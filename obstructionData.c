#include "track.h"
#include "math.h"
#include "stdio.h"

OBSTRUCTION_DATASETS obstructionDataSets = {0, 0, 0, NULL, 0};

void readObstructionInput(NAMELIST_TEXT *nltext, RUN *run)
{
  SDDS_DATASET SDDSin;
  char s[16384];
  long code;
  
#include "obstructionData.h"

  set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
  set_print_namelist_flags(0);
  if (processNamelist(&obstruction_data, nltext)==NAMELIST_ERROR)
    bombElegant(NULL, NULL);
  if (echoNamelists) print_namelist(stdout, &obstruction_data);

  resetObstructionData(&obstructionDataSets);

  if (disable)
    return;
  
  if (!SDDS_InitializeInputFromSearchPath(&SDDSin, input)) {
    sprintf(s, "Problem opening aperture input file %s", input);
    SDDS_SetError(s);
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
  }
  
  if (!check_sdds_column(&SDDSin, "Z", "m") ||
      !check_sdds_column(&SDDSin, "X", "m")) {
    printf("Necessary data quantities (Z, X) have wrong units or are not present in %s\n",
            input);
    printf("Note that units must be \"m\" on all quantities\n");
    fflush(stdout);
    exitElegant(1);
  }

  obstructionDataSets.periodic = periodic;
  obstructionDataSets.superperiodicity = superperiodicity;
  
  while ((code=SDDS_ReadPage(&SDDSin))>0) {
    obstructionDataSets.data = SDDS_Realloc(obstructionDataSets.data,
                                            sizeof(*(obstructionDataSets.data))*(obstructionDataSets.nDataSets+1));;
    if ((obstructionDataSets.data[obstructionDataSets.nDataSets].points=SDDS_RowCount(&SDDSin))<3) {
      printf("Obstruction input file %s has fewer than 3 rows on page %ld", input, code);
      fflush(stdout);
      exitElegant(1);
    }
    if (!(obstructionDataSets.data[obstructionDataSets.nDataSets].Z = SDDS_GetColumnInDoubles(&SDDSin, "Z")) ||
        !(obstructionDataSets.data[obstructionDataSets.nDataSets].X = SDDS_GetColumnInDoubles(&SDDSin, "X"))) {
      sprintf(s, "Problem getting data from page %ld of obstruction input file %s", code, input);
      SDDS_SetError(s);
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
    }
    obstructionDataSets.nDataSets += 1;
  }
  if (code==0 || obstructionDataSets.nDataSets<1) {
    sprintf(s, "Problem reading obstruction  input file %s---seems to be empty", input);
    SDDS_SetError(s);
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
  }
  printf("%ld datasets read from obstruction file %s\n", obstructionDataSets.nDataSets, input);

  obstructionDataSets.initialized = 1;

  return;
}

/* Returns 0 if not obstructed */
long insideObstruction(double *part, long segment, long nSegments)
{
  TRACKING_CONTEXT context;
  ELEMENT_LIST *eptr;
  /*
  static FILE *fpObs = NULL;
  static ELEMENT_LIST *lastEptr = NULL;
  */
  long ic;
  double Z, X, Y;

  /*
  if (!fpObs) {
    fpObs = fopen("globalPart.sdds", "w");
    fprintf(fpObs, "SDDS1\n&column name=Z type=double units=m &end\n");
    fprintf(fpObs, "&column name=X type=double units=m &end\n");
    fprintf(fpObs, "&column name=particleID type=long &end\n");
    fprintf(fpObs, "&data mode=ascii no_row_counts=1 &end\n");
  }
  */

  if (!obstructionDataSets.initialized) return 0;

  getTrackingContext(&context);
  if (!(eptr=context.element)) return 0;

  /*
  if (eptr!=lastEptr) {
    printf("%s#%04ld: Z=%le, X=%le, theta=%le\n", 
           eptr->name, eptr->occurence, eptr->floorCoord[2], eptr->floorCoord[0], eptr->floorAngle[0]);
    lastEptr = eptr;
  }
  */

  convertLocalCoordinatesToGlobal(&Z, &X, &Y, part, eptr, segment, nSegments);
  /* 
  fprintf(fpObs, "%le %le %ld\n", Z, X, (long)part[6]); 
  */

  for (ic=0; ic<obstructionDataSets.nDataSets; ic++) {
    if (pointIsInsideContour(Z, X, 
                             obstructionDataSets.data[ic].Z, 
                             obstructionDataSets.data[ic].X, 
                             obstructionDataSets.data[ic].points)) {
      return 1;
    }
  }
  return 0;
}

long insideObstruction_xy(double x, double y, double xyTilt, long particleID, long segment, long nSegments)
{
  double part[7] ={0,0,0,0,0,0,0};
  double sin_tilt, cos_tilt;
  sin_tilt = sin(xyTilt);
  cos_tilt = cos(xyTilt);
  part[0] =  x*cos_tilt - y*sin_tilt;
  part[2] =  x*sin_tilt + y*cos_tilt;
  part[6] = particleID;
  return insideObstruction(part, segment, nSegments);
}

long insideObstruction_xy_dz(double x, double y, double xyTilt, long particleID, double dz)
{
  TRACKING_CONTEXT context;
  ELEMENT_LIST *eptr;
  long ic;
  double Z, X, Y;
  double part[7] ={0,0,0,0,0,0,0};
  double sin_tilt, cos_tilt;

  if (!obstructionDataSets.initialized) return 0;

  getTrackingContext(&context);
  if (!(eptr=context.element)) return 0;

  sin_tilt = sin(xyTilt);
  cos_tilt = cos(xyTilt);
  part[0] =  x*cos_tilt - y*sin_tilt;
  part[2] =  x*sin_tilt + y*cos_tilt;
  part[6] = particleID;

  convertLocalCoordinatesToGlobal(&Z, &X, &Y, part, eptr, 0, 1);
  Z += dz;

  for (ic=0; ic<obstructionDataSets.nDataSets; ic++) {
    if (pointIsInsideContour(Z, X, 
                             obstructionDataSets.data[ic].Z, 
                             obstructionDataSets.data[ic].X, 
                             obstructionDataSets.data[ic].points))
      return 1;
  }
  return 0;
}


long filterParticlesWithObstructions(double **coord, long np, double **accepted, double z, double Po,
                                     long segment, long nSegments)
{
  long ip, itop;
  itop = np - 1;
  for (ip=0; ip<=itop; ip++) {
    if (insideObstruction(coord[ip], segment, nSegments)) {
      TRACKING_CONTEXT context;
      getTrackingContext(&context);
      if (ip!=itop)
        swapParticles(coord[ip], coord[itop]);
      coord[itop][4] = z;
      coord[itop][5] = Po*(1+coord[itop][5]);
      if (accepted)
        swapParticles(accepted[ip], accepted[itop]);
      --itop;
      --ip;
    } 
  }
  return itop+1;
}

void resetObstructionData(OBSTRUCTION_DATASETS *obsData)
{
  long i;
  if (obsData->initialized) {
    for (i=0; i<obsData->nDataSets; i++)  {
      free(obsData->data[i].X);
      free(obsData->data[i].Z);
    }
    free(obsData->data);
    obsData->data = NULL;
    obsData->initialized = 0;
  }
}


