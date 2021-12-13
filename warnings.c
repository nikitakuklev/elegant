#include "mdb.h"
#include "track.h"

static long warnings = 0;
static char **warningText = NULL;
static long *warningCount = NULL;
static FILE *fpWarn = NULL;
static htab *hash_table = NULL;

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
  char warningBuffer[1024];
  getTrackingContext(&trackingContext);
  if (!strlen(trackingContext.elementName) || !trackingContext.elementOccurrence || !trackingContext.element)
    printWarning(text, detail);
  else {
    snprintf(warningBuffer, 1024, "%s#%ld", trackingContext.elementName, trackingContext.elementOccurrence);
    printWarningWithContext(entity_name[trackingContext.element->type], warningBuffer, 
                            text, detail);
  }
}

void printWarningWithContext(char *context1, char  *context2, char *text,  char *detail)
{
  long *counterPtr;

  if (!hash_table)
    hash_table = hcreate(12);
  if (hcount(hash_table)==0 ||
      hfind(hash_table, text, strlen(text))==FALSE) {
    if (!(warningText = SDDS_Realloc(warningText, sizeof(*warningText)*(warnings+1))))
      bombElegant("Memory allocation error in printWarning\n", NULL);
    if (!(warningCount = SDDS_Realloc(warningCount, sizeof(*warningCount)*(warnings+1))))
      bombElegant("Memory allocation error in printWarning\n", NULL);
    cp_str(&warningText[warnings], text);
    warningCount[warnings] = 1;
    counterPtr = &warningCount[warnings];
    hadd(hash_table, warningText[warnings], strlen(warningText[warnings]), 
         (void*)&warningCount[warnings]);
    warnings++;
  } else {
    counterPtr = hstuff(hash_table);
    *counterPtr += 1;
  }

  if (!fpWarn)
    fpWarn = stdout;
  if ((*counterPtr)<=warningCountLimit) {
    fprintf(fpWarn, "*** Warning: %s", text);
    if (context1 && strlen(context1)) {
      if (context2 && strlen(context2)) {
        /* context1 is typically the element type name, context2 is typically the element name */
        if (detail && strlen(detail))
          fprintf(fpWarn, " --- %s %s: %s\n", context1, context2, detail);
        else 
          fprintf(fpWarn, " --- %s %s\n", context1, context2);
      } else {
        /* context1 is typically the command name or subroutine name */
        if (detail && strlen(detail))
          fprintf(fpWarn, " --- %s: %s\n", context1, detail);
        else
          fprintf(fpWarn, " --- %s\n", context1);
      }
    } else {
      if (detail && strlen(detail))
        fprintf(fpWarn, " --- %s\n", detail);
      else
        fprintf(fpWarn, "\n");
    }
    if ((*counterPtr)==warningCountLimit)
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
    fprintf(fpWarn, "%ld types of warnings were recorded:\n", warnings);
    for (i=0; i<warnings; i++) {
      fprintf(fpWarn, "%ld* %s\n",
             warningCount[i], warningText[i]);
      free(warningText[i]);
      warningCount[i] = 0;
      warningText[i] = NULL;
    }
    fprintf(fpWarn, "*****************************************************************************\n\n");

    warnings = 0;
    hdestroy(hash_table);
    hash_table = NULL;
    free(warningText);
    free(warningCount);
  }
}
