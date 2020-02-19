#include "mdb.h"
#include "track.h"

static long warnings = 0;
static char **warningText = NULL;
static long *warningCount = NULL;

void printWarning(char *text,  char *detail)
{
  long i;
  if (detail && strlen(detail))
    printf("*** Warning: %s%s\n", text, detail);
  else
    printf("*** Warning: %s.\n", text);
  fflush(stdout);
  for (i=0; i<warnings; i++) 
    if (strcmp(warningText[i], text)==0)
      break;
  if (i==warnings)  {
    if (!(warningText = SDDS_Realloc(warningText, sizeof(*warningText)*(warnings+1))))
      bombElegant("Memory allocation error in printWarning\n", NULL);
    if (!(warningCount = SDDS_Realloc(warningCount, sizeof(*warningCount)*(warnings+1))))
      bombElegant("Memory allocation error in printWarning\n", NULL);
    cp_str(&warningText[warnings], text);
    warningCount[warnings] = 1;
    warnings++;
  } else 
    warningCount[i] += 1;
}

void summarizeWarnings()
{
  long i;
  if (warnings) {
    printf("\n************************** Summary of warnings ****************************\n");
    printf("*** NB: the warning summary is still in development and covers only some warnings.\n\n");
    printf("%ld types of warnings were recorded:\n", warnings);
    for (i=0; i<warnings; i++) {
      printf("%ld* %s\n",
             warningCount[i], warningText[i]);
      free(warningText[i]);
      warningCount[i] = 0;
      warningText[i] = NULL;
    }
    printf("\n*** NB: the warning summary is still in development and covers only some warnings.\n");
    printf("*****************************************************************************\n\n");

    warnings = 0;
  }
}
