#include "track.h"
#include "mdb.h"
#include "math.h"
#include "stdio.h"

OBSTRUCTION_DATASETS obstructionDataSets = {0, 0, 0, {0.0, 0.0}, {-10.0, 10.0}, NULL, 0};

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
  obstructionDataSets.yLimit[0] = yLimit[0];
  obstructionDataSets.yLimit[1] = yLimit[1];

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
  if (!check_sdds_parameter(&SDDSin, "ZCenter", "m") ||
      !check_sdds_parameter(&SDDSin, "XCenter", "m") ||
      !check_sdds_parameter(&SDDSin, "Superperiodicity", NULL) ) {
    printf("Necessary data quantities (ZCenter, XCenter, Superperiodicity) have wrong type or units, or are not present in %s\n",
            input);
    fflush(stdout);
    exitElegant(1);
  }

  obstructionDataSets.periods = periods;

  while ((code=SDDS_ReadPage(&SDDSin))>0) {
    if (code==1) {
      int32_t superperiodicity;
      if (!SDDS_GetParameterAsDouble(&SDDSin, "ZCenter", &obstructionDataSets.center[0]) ||
	  !SDDS_GetParameterAsDouble(&SDDSin, "XCenter", &obstructionDataSets.center[1]) ||
	  !SDDS_GetParameterAsLong(&SDDSin, "Superperiodicity", &superperiodicity)) {
	sprintf(s, "Problem getting data from page %ld of obstruction input file %s", code, input);
	SDDS_SetError(s);
	SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
      }
      obstructionDataSets.superperiodicity = superperiodicity;
      printf("ZCenter = %le, XCenter = %le, Superperiodicity = %ld\n", 
	     obstructionDataSets.center[0], obstructionDataSets.center[1], obstructionDataSets.superperiodicity);
    }
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

#define DEBUG 1
void logInside(double X, double Z, long particleID, short where)
{
#ifdef DEBUG
  static FILE *fpInside = NULL;
  TRACKING_CONTEXT context;
  if (!fpInside) {
    char buffer[1024];
#if USE_MPI
    snprintf(buffer, 1024, "insideObstruction-%04d.sdds", myid);
#else
    snprintf(buffer, 1024, "insideObstruction.sdds");
#endif
    fpInside = fopen(buffer, "w");
    fprintf(fpInside, "SDDS1\n&column name=Z type=double units=m &end\n");
    fprintf(fpInside, "&column name=X type=double units=m &end\n");
    fprintf(fpInside, "&column name=particleID type=long &end\n");
    fprintf(fpInside, "&column name=call type=short &end\n");
    fprintf(fpInside, "&column name=ElementName type=string &end\n");
    fprintf(fpInside, "&column name=ElementType type=string &end\n");
    fprintf(fpInside, "&data mode=ascii no_row_counts=1 &end\n");
  }
  getTrackingContext(&context);
  fprintf(fpInside, "%le %le %ld %hd %s %s\n", Z, X, particleID, where, 
	  context.element->name, entity_name[context.element->type]);
  fflush(fpInside);
#endif
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
  long ic, iperiod, lost;
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

  if (!obstructionDataSets.initialized) {
    /*
    printf("insideObstruction: obstructions not initialized\n");
    */
    return 0;
  }
  
  if (!obstructionsInForce) {
    /*
    printf("insideObstruction: obstructions not in force\n");
    */
    return 0;
  }

  getTrackingContext(&context);
  if (!(eptr=context.element)) {
    printf("No element pointer in insideObstruction()\n");
    return 0;
  }

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

  if (obstructionDataSets.yLimit[0]<obstructionDataSets.yLimit[1] && 
      (Y<obstructionDataSets.yLimit[0] || Y>obstructionDataSets.yLimit[1]))
    lost = 1;
  else {
    lost = 0;
    for (iperiod=0; iperiod<obstructionDataSets.periods && !lost; iperiod++) {
      for (ic=0; ic<obstructionDataSets.nDataSets; ic++) {
	if (pointIsInsideContour(Z, X, 
				 obstructionDataSets.data[ic].Z, 
				 obstructionDataSets.data[ic].X, 
				 obstructionDataSets.data[ic].points,
				 obstructionDataSets.center,
				 (iperiod*PIx2)/obstructionDataSets.superperiodicity)) {
	  lost = 1;
	  break;
	}
      }
    }
  }
  if (lost) {
    logInside(X, Z, (long)part[particleIDIndex], 1);
    if (globalLossCoordOffset!=-1) {
      part[globalLossCoordOffset+0] = X;
      part[globalLossCoordOffset+1] = Y;
      part[globalLossCoordOffset+2] = Z;
    }
    return 1;
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
  double part[MAX_PROPERTIES_PER_PARTICLE];
  double sin_tilt, cos_tilt;

  if (!obstructionDataSets.initialized || !obstructionsInForce) return 0;

  if (xyTilt) {
    sin_tilt = sin(xyTilt);
    cos_tilt = cos(xyTilt);
  } else {
    sin_tilt = 0;
    cos_tilt = 1;
  }
  memset(&part[0], 0, sizeof(double)*MAX_PROPERTIES_PER_PARTICLE);
  part[0] =  x*cos_tilt - y*sin_tilt;
  part[2] =  x*sin_tilt + y*cos_tilt;
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
  long ic, iperiod, lost;

  /*
  static FILE *fp = NULL;
  if (!fp) {
    fp = fopen("obstructPath.sdds", "w");
    fprintf(fp, "SDDS1\n&column name=Z type=double units=m &end\n");
    fprintf(fp, "&column name=X type=double units=m &end\n");
    fprintf(fp, "&column name=dZi type=double units=m &end\n");
    fprintf(fp, "&column name=dXi type=double units=m &end\n");
    fprintf(fp, "&column name=Z1 type=double units=m &end\n");
    fprintf(fp, "&column name=X1 type=double units=m &end\n");
    fprintf(fp, "&column name=thetai type=double &end\n");
    fprintf(fp, "&data mode=ascii no_row_counts=1 &end\n");
  }
  fprintf(fp, "%21.15e %21.15e %21.15e %21.15e ", Z, X, dZi, dXi);
  */

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
    
    /*
      fprintf(fp, "%21.15e %21.15e %21.15e\n", Z1, X1, thetai);
  */

  if (obstructionDataSets.yLimit[0]<obstructionDataSets.yLimit[1] && 
      (Y<obstructionDataSets.yLimit[0] || Y>obstructionDataSets.yLimit[1]))
    lost = 1;
  else {
    lost = 0;
    for (iperiod=0; iperiod<obstructionDataSets.periods && !lost; iperiod++) {
      for (ic=0; ic<obstructionDataSets.nDataSets; ic++) {
	if (pointIsInsideContour(Z1, X1, 
				 obstructionDataSets.data[ic].Z, 
				 obstructionDataSets.data[ic].X, 
				 obstructionDataSets.data[ic].points,
				 obstructionDataSets.center, 
				 (iperiod*PIx2)/obstructionDataSets.superperiodicity)) {
	  lost = 1;
	  break;
	}
      }
    }
  }
  if (lost) {
    if (lossCoordinates) {
      lossCoordinates[0] = X1;
      lossCoordinates[1] = Y1;
      lossCoordinates[2] = Z1;
    }
    return 1;
  }
  return 0;
}

long filterParticlesWithObstructions(double **coord, long np, double **accepted, double z, double Po)
{
  long ip, itop;
  itop = np - 1;
  for (ip=0; ip<=itop; ip++) {
    /* printf("filterParticlesWithObstructions: ip=%ld, itop=%ld\n", ip, itop); */
    if (insideObstruction(coord[ip], GLOBAL_LOCAL_MODE_END, 0.0, 0, 1)) {
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
    obsData->yLimit[0] = -10;
    obsData->yLimit[1] = 10;
  }
}


