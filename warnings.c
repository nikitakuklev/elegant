#include "mdb.h"
#include "track.h"

static long warnings = 0;
static char **warningText = NULL;
static long *warningCount = NULL;
static FILE *fpWarn = NULL;

void setWarningFilePointer(FILE *fp)
{
  fpWarn = fp;
}

void printWarning(char *text,  char *detail)
{
  printWarningWithContext(NULL, NULL, text, detail);
}

void printWarningForTracking(char *text, char *detail)
/* extracts element info from tracking context */
{
  TRACKING_CONTEXT trackingContext;
  getTrackingContext(&trackingContext);
  if (!strlen(trackingContext.elementName) || trackingContext.elementOccurrence || !trackingContext.element)
    printWarning(text, detail);
  printWarningWithContext(entity_name[trackingContext.element->type], trackingContext.elementName,
                          text, detail);
}

void printWarningWithContext(char *context1, char  *context2, char *text,  char *detail)
{
  long i;
  char buffer[32768];
  for (i=0; i<warnings; i++) 
    if (strcmp(warningText[i], buffer)==0)
      break;
  if (i==warnings)  {
    if (!(warningText = SDDS_Realloc(warningText, sizeof(*warningText)*(warnings+1))))
      bombElegant("Memory allocation error in printWarning\n", NULL);
    if (!(warningCount = SDDS_Realloc(warningCount, sizeof(*warningCount)*(warnings+1))))
      bombElegant("Memory allocation error in printWarning\n", NULL);
    cp_str(&warningText[warnings], buffer);
    warningCount[warnings] = 1;
    warnings++;
  } else 
    warningCount[i] += 1;

  if (!fpWarn)
    fpWarn = stdout;
  if (context1 && strlen(context1)) {
    if (context2 && strlen(context2))
      /* context1 is typically the element type name, context2 is typically the element name */
      snprintf(buffer, 32768, "%s %s: %s", context1, context2, text);
    else
      /* context1 is typically the command name or subroutine name */
      snprintf(buffer, 32768, "%s: %s", context1, text);
  }
  else
    strncpy(buffer, text, 32768);
  if (warningCount[i]<=warningCountLimit) {
    if (detail && strlen(detail))
      fprintf(fpWarn, "*** Warning: %s---%s\n", buffer, detail);
    else
      fprintf(fpWarn, "*** Warning: %s\n", buffer);
    if (warningCount[i]==warningCountLimit)
      fprintf(fpWarn, "Further warnings of this type will be suppressed.\n");
    fflush(fpWarn);
  }
}

void summarizeWarnings()
{
  long i;
  if (!fpWarn)
    fpWarn = stdout;
  if (warnings) {
    fprintf(fpWarn, "\n************************** Summary of warnings ****************************\n");
    fprintf(fpWarn, "*** NB: the warning summary is still in development and covers only some warnings.\n\n");
    fprintf(fpWarn, "%ld types of warnings were recorded:\n", warnings);
    for (i=0; i<warnings; i++) {
      fprintf(fpWarn, "%ld* %s\n",
             warningCount[i], warningText[i]);
      free(warningText[i]);
      warningCount[i] = 0;
      warningText[i] = NULL;
    }
    fprintf(fpWarn, "\n*** NB: the warning summary is still in development and covers only some warnings.\n");
    fprintf(fpWarn, "*****************************************************************************\n\n");

    warnings = 0;
  }
}
