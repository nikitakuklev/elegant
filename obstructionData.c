#include "track.h"
#include "math.h"
#include "stdio.h"

OBSTRUCTION_DATASETS obstructionDataSets = {0, 0, 0, NULL, 0};

static long obstructionsInForce = 1;
void setObstructionsMode(long state) 
{
  obstructionsInForce = state;
}

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
long insideObstruction(double *part, short mode, double dz, long segment, long nSegments)
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

  if (!obstructionDataSets.initialized || !obstructionsInForce) return 0;

  getTrackingContext(&context);
  if (!(eptr=context.element)) return 0;

  /*
  if (eptr!=lastEptr) {
    printf("%s#%04ld: Z=%le, X=%le, theta=%le\n", 
           eptr->name, eptr->occurence, eptr->floorCoord[2], eptr->floorCoord[0], eptr->floorAngle[0]);
    lastEptr = eptr;
  }
  */

  convertLocalCoordinatesToGlobal(&Z, &X, &Y, mode, part, eptr, dz, segment, nSegments);
  /*
    printf("Checking obstruction for x=%le, y=%le, s=%le, segment=%ld/%ld: Z=%le, X=%le, Y=%le\n",
         part[0], part[2], part[4], segment, nSegments, Z, X, Y);
  */

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

long insideObstruction_xyz
(
 double x, /* local x coordinate */
 double y, /* local y coordinate */
 double xyTilt, /* tilt of element */
 short mode,
 double dz,
 long segment,  /* for segmented elements, the segment index */
 long nSegments /* For segmented elements, the number of segments.
                * If non-positive, the s value is used to determine the
                * longitudinal position more accurately.
                */
 )
{
  double part[7] ={0,0,0,0,0,0,0};
  double sin_tilt, cos_tilt;

  if (!obstructionDataSets.initialized || !obstructionsInForce) return 0;

  if (xyTilt) {
    sin_tilt = sin(xyTilt);
    cos_tilt = cos(xyTilt);
  } else {
    sin_tilt = 0;
    cos_tilt = 1;
  }
  part[0] =  x*cos_tilt - y*sin_tilt;
  part[2] =  x*sin_tilt + y*cos_tilt;
  part[4] = 0;
  part[6] = 0;
  return insideObstruction(part, mode, dz, segment, nSegments);
}

long insideObstruction_XYZ
/* Used for elements with internal Cartesian coordinate system that
 * may be offset and rotated relative to global coordinates.
 * E.g., BRAT element.
 */
(
 /* magnet-frame coordinates of particle */
 double X  , double Y  , double Z  , 
 /* coordinate offsets of the nominal entrance */
 double dXi, double dYi, double dZi, 
 /* angle of the internal Z axis w.r.t. nominal incoming trajectory */
 double thetai,
 /* return of global loss coordinates (X, Y, Z) */
 double *lossCoordinates
 )
{
  TRACKING_CONTEXT context;
  double C, S;
  double X1, Y1, Z1;
  long ic;
#ifdef DEBUG
  static FILE *fp = NULL;

  if (!fp) {
    fp = fopen("obstructPath.sdds", "w");
    fprintf(fp, "SDDS1\n&column name=Z type=double units=m &end\n");
    fprintf(fp, "&column name=X type=double units=m &end\n");
    fprintf(fp, "&column name=thetai type=double &end\n");
    fprintf(fp, "&data mode=ascii no_row_counts=1 &end\n");
  }
#endif

  if (!obstructionDataSets.initialized || !obstructionsInForce) return 0;

  X -= dXi;
  Y -= dYi;
  Z -= dZi;

  getTrackingContext(&context);
  if (context.element->pred)
    thetai -= context.element->pred->floorAngle[0];

  C = cos(thetai);
  S = sin(thetai);
  X1 =  C*X - S*Z;
  Y1 = Y;
  Z1 =  S*X + C*Z;
  if (context.element->pred) {
    X1 += context.element->pred->floorCoord[0];
    Y1 += context.element->pred->floorCoord[1];
    Z1 += context.element->pred->floorCoord[2];
  }

#ifdef DEBUG
  fprintf(fp, "%le %le %le\n", Z1, X1, thetai);
  return 0;
#endif

  for (ic=0; ic<obstructionDataSets.nDataSets; ic++) {
    if (pointIsInsideContour(Z1, X1, 
                             obstructionDataSets.data[ic].Z, 
                             obstructionDataSets.data[ic].X, 
                             obstructionDataSets.data[ic].points)) {
      if (lossCoordinates) {
        lossCoordinates[0] = X1;
        lossCoordinates[1] = Y1;
        lossCoordinates[2] = Z1;
      }
      /*
      printf("Lost on obstruction: %le, %le, %le\n", X, Y, Z);
      fflush(stdout);
      */
      return 1;
    }
  }
  return 0;
}

long filterParticlesWithObstructions(double **coord, long np, double **accepted, double z, double Po)
{
  long ip, itop;
  itop = np - 1;
  for (ip=0; ip<=itop; ip++) {
    if (insideObstruction(coord[ip], GLOBAL_LOCAL_MODE_END, 0.0, 1, 1)) {
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


