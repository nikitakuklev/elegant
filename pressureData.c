/*************************************************************************\
* Copyright (c) 2017 The University of Chicago, as Operator of Argonne
* National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: ionEffects.c
 * purpose: simulation of ion interaction with the beam
 *
 * Joe Calvey, Michael Borland 2017
 */
#if defined(SOLARIS) && !defined(__GNUC__)
#include <sunmath.h>
#endif

#include "mdb.h"
#include "track.h"
#include "ionEffects.h"
#include "constants.h"
#include "pressureData.h"

void readGasPressureData(char *filename, PRESSURE_DATA *pressureData)
{
  /* Assumed file structure:
   * Parameters: Gasses --- SDDS_STRING giving comma- or space-separated list of gas species, e.g., "H2O H2 N2 O2 CO2 CO CH4"
   * Columns:
   * s         --- SDDS_FLOAT or SDDS_DOUBLE giving location in the lattice
   * <gasName> --- SDDS_FLOAT or SDDS_DOUBLE giving pressure of <gasName> in Torr or nT
   */

  SDDS_DATASET SDDSin;
  char *gasColumnList, *ptr;
  long i;
  double dsMin, dsMax, ds, pressureMultiplier;

  if (!SDDS_InitializeInputFromSearchPath(&SDDSin, filename)) {
    printf("Problem opening pressure data file %s\n", filename);
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
  }

  if (!check_sdds_column(&SDDSin, "s", "m"))
    bombElegantVA("Column 's' is missing or does not have units of 'm' in %s", filename);
  if (SDDS_CheckParameter(&SDDSin, "Gasses", NULL, SDDS_STRING, stdout)!=SDDS_CHECK_OK)
    bombElegantVA("Parameters \"Gasses\" is missing or not string type in %s", filename);
  
  if (SDDS_ReadPage(&SDDSin)<=0) 
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);

  if (!SDDS_GetParameter(&SDDSin, "Gasses", &gasColumnList))
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
  
  pressureData->nGasses = 0;
  pressureData->gasName = NULL;
  while ((ptr=get_token(gasColumnList))!=NULL) {
    pressureData->gasName = (char**)SDDS_Realloc(pressureData->gasName, sizeof(*(pressureData->gasName))*(pressureData->nGasses+1));
    cp_str(&pressureData->gasName[pressureData->nGasses], ptr);
    pressureData->nGasses += 1;
  }
  printf("%ld gasses listed in %s: ", pressureData->nGasses, filename);
  for (i=0; i<pressureData->nGasses; i++)
    printf("%s%c", pressureData->gasName[i], i==(pressureData->nGasses-1)?'\n':' ');
  free(gasColumnList);

  pressureData->nLocations = SDDS_RowCount(&SDDSin);
  printf("Gas data provided at %ld s locations\n", pressureData->nLocations);
  if (!(pressureData->s=SDDS_GetColumnInDoubles(&SDDSin, "s")))
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);

  pressureData->pressure = (double**)czarray_2d(sizeof(double), pressureData->nLocations, pressureData->nGasses);
  for (i=0; i<pressureData->nGasses; i++) {
    pressureMultiplier = 1;
    if (!check_sdds_column(&SDDSin, pressureData->gasName[i], "Torr") && !check_sdds_column(&SDDSin, pressureData->gasName[i], "T")) {
      pressureMultiplier = 1e-9;
      if (!check_sdds_column(&SDDSin, pressureData->gasName[i], "nT") && !check_sdds_column(&SDDSin, pressureData->gasName[i], "nTorr"))
        bombElegantVA("Column \"%s\" is missing, not floating-point type, or does not have units of \"Torr\" or \"nT\" in %s", 
                      pressureData->gasName[i], filename);
    }
    if (!(pressureData->pressure[i] = SDDS_GetColumnInDoubles(&SDDSin, pressureData->gasName[i]))) {
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
    }
    if (pressureMultiplier!=1) {
      /* Convert to Torr */
      long j;
      for (j=0; j<pressureData->nLocations; j++)
        pressureData->pressure[i][j] *= pressureMultiplier;
    }
  }

  dsMax = -(dsMin = DBL_MAX);
  for (i=1; i<pressureData->nLocations; i++) {
    ds = pressureData->s[i] - pressureData->s[i-1];
    if (ds<=0)
      bombElegantVA("s data is not monotonically increasing in pressure data file %s (%le, %le)", filename, pressureData->s[i-1], pressureData->s[i]);
    if (dsMin>ds)
      dsMin = ds;
    if (dsMax<ds)
      dsMax = ds;
  }
  if (fabs(1-dsMin/dsMax)>1e-3)
      bombElegantVA("s data is not uniformly spaced to within desired 0.1% in pressure data file %s", filename);

  printf("Finished reading pressure data file %s\n", filename);
  fflush(stdout);
}

void computeAverageGasPressures(double sStart, double sEnd, double *pressure, PRESSURE_DATA *pressureData)
{
  double sum;
  long iGas, iLocation, iStart, iEnd;

  /* Could improve this by interpolating the pressure at sStart and sEnd */

  /* Find the indices spanning the desired region */
  iStart = iEnd = -1;
  for (iLocation=0; iLocation<pressureData->nLocations; iLocation++) {
    if (pressureData->s[iLocation]>=sStart && iStart==-1)
      iStart = iLocation;
    if (pressureData->s[iLocation]<=sEnd)
      iEnd = iLocation;
    else
      break;
  }
  if (iStart==-1 || iEnd==-1 || iEnd<=iStart)
    bombElegantVA("Failed to find indices corresponding to pressure region s:[%le, %le] m\n",
                  sStart, sEnd);

  for (iGas=0; iGas<pressureData->nGasses; iGas++) {
    sum = 0;
    for (iLocation=iStart; iLocation<=iEnd; iLocation++) 
      sum += pressureData->pressure[iGas][iLocation];
    pressure[iGas] = sum/(iEnd-iStart+1);
  }
}
