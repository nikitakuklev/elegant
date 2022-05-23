/************************************************************************* \
* Copyright (c) 2022 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2022 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/
#define DEBUG 

/* routine: lgbend()
 * purpose: tracking through canonical bend in cartesian coordinates
 * 
 * Michael Borland, 2022
 */
#include "mdb.h"
#include "track.h"
#include "multipole.h"

static LGBEND lgbendCopy;
static ELEMENT_LIST *eptrCopy = NULL;
static double PoCopy, lastRho, lastX, lastXp;

#ifdef DEBUG
static FILE *fpDeb = NULL;
#endif

void readLGBendConfiguration(LGBEND *lgbend, ELEMENT_LIST *eptr);
void flipLGBEND(LGBEND *lgbend);

void switchLgbendPlane(double **particle, long n_part, double dx, double alpha, double po);
void lgbendFringe(double **particle, long n_part, double alpha, double rho0, double K1, double K2,
                  LGBEND_SEGMENT *segment, short angleSign, short isExit);
int integrate_kick_KnL(double *coord, double dx, double dy, 
                      double Po, double rad_coef, double isr_coef,
                      double *KnL, long nTerms,
                      long integration_order, long n_parts, long iPart, long iFinalSlice,
                      double drift,
                      MULTIPOLE_DATA *multData, MULTIPOLE_DATA *edge1MultData, MULTIPOLE_DATA *edge2MultData, 
                      MULT_APERTURE_DATA *apData, double *dzLoss, double *sigmaDelta2);
double lgbend_trajectory_error(double *value, long *invalid);

long track_through_lgbend
(
 double **particle,   /* initial/final phase-space coordinates */
 long n_part,         /* number of particles */
 ELEMENT_LIST *eptr,
 LGBEND *lgbend,
 double Po,
 double **accepted,
 double z_start,
 double *sigmaDelta2, /* use for accumulation of energy spread for radiation matrix computation */
 char *rootname,
 MAXAMP *maxamp,
 APCONTOUR *apContour,
 APERTURE_DATA *apFileData,
 /* If iPart non-negative, we do one step. The caller is responsible 
  * for handling the coordinates appropriately outside this routine. 
  * The element must have been previously optimized to determine FSE and X offsets.
  */
 long iPart,
 /* If iFinalSlice is positive, we terminate integration inside the magnet. The caller is responsible 
  * for handling the coordinates appropriately outside this routine. 
  * The element must have been previously optimized to determine FSE and X offsets.
  */
 long iFinalSlice
 )
{
  double KnL[9];
  long iTerm, nTerms, iSegment;
  double dx, dy, dz; /* offsets of the multipole center */
  long nSlices, integ_order;
  long i_part, i_top;
  double *coef;
  double fse, etilt, tilt, rad_coef, isr_coef, dzLoss=0;
  double rho0, length, angle, angleSign, extraTilt;
  double entryAngle, exitAngle, entryPosition, exitPosition;
  MULT_APERTURE_DATA apertureData;

#ifdef DEBUG
  if (!fpDeb) {
    char buffer[1024];
    TRACKING_CONTEXT context;
    getTrackingContext(&context);
    snprintf(buffer, 1024, "%s.lgb", context.rootname);
    fpDeb = fopen(buffer, "w");
    fprintf(fpDeb, "SDDS1\n&column name=segment type=short &end\n");
    fprintf(fpDeb, "&column name=z type=double units=m &end\n");
    fprintf(fpDeb, "&column name=s type=double units=m &end\n");
    fprintf(fpDeb, "&column name=x type=double units=m &end\n");
    fprintf(fpDeb, "&column name=xp type=double &end\n");
    fprintf(fpDeb, "&column name=y type=double units=m &end\n");
    fprintf(fpDeb, "&column name=yp type=double &end\n");
    fprintf(fpDeb, "&column name=part type=string &end\n");
    fprintf(fpDeb, "&data mode=ascii no_row_counts=1 &end\n");
  }
#endif

  if (!particle)
    bombTracking("particle array is null (track_through_lgbend)");

  if (iPart>=0 && lgbend->optimized!=1)
    bombTracking("Programming error: one-step mode invoked for unoptimized LGBEND.");
  if (iFinalSlice>0 && iPart>=0)
    bombTracking("Programming error: partial integration mode and one-step mode invoked together for LGBEND.");
  if (iFinalSlice>0 && lgbend->optimized!=1)
    bombTracking("Programming error: partial integration mode invoked for unoptimized LGBEND.");

  if (N_LGBEND_FRINGE_INT!=8)
    bombTracking("Coding error detected: number of fringe integrals for LGBEND is different than expected.");

  if (!lgbend->initialized)
    readLGBendConfiguration(lgbend, eptr);

  if (lgbend->optimizeFse && !lgbend->optimized && lgbend->angle!=0) {
    double acc;
    double startValue[2], stepSize[2], lowerLimit[2], upperLimit[2];
    short disable[2];

    PoCopy = lastRho = lastX = lastXp = 0;
    if (iPart>=0)
      bombTracking("Programming error: oneStep mode is incompatible with optmization for LGBEND.");
    startValue[0] = lgbend->segment[0].fse;
    startValue[1] = lgbend->segment[lgbend->nSegments-1].fse;
    lgbend->optimized = -1; /* flag to indicate calls to track_through_lgbend will be for FSE optimization */
    memcpy(&lgbendCopy, lgbend, sizeof(lgbendCopy));
    eptrCopy = eptr;
    lgbendCopy.fse = lgbendCopy.dx = lgbendCopy.dy = lgbendCopy.dz = 
      lgbendCopy.etilt = lgbendCopy.tilt = lgbendCopy.isr = lgbendCopy.synch_rad = lgbendCopy.isr1Particle = 0;

    PoCopy = Po;
    stepSize[0] = stepSize[1] = 1e-3;
    lowerLimit[0] = lowerLimit[1] = -1;
    upperLimit[0] = upperLimit[1] = 1;
    disable[0] = disable[1] = 0;
    if (simplexMin(&acc, startValue, stepSize, lowerLimit, upperLimit, disable, 2, 
                   fabs(1e-15*lgbend->length), fabs(1e-16*lgbend->length),
                   lgbend_trajectory_error, NULL, 1500, 3, 12, 3.0, 1.0, 0)<0) {
      bombElegantVA("failed to find FSE and x offset to center trajectory for lgbend. accuracy acheived was %le.", acc);
    }
    lgbend->segment[0].fse = startValue[0];
    lgbend->segment[lgbend->nSegments-1].fse = startValue[1];

    lgbend->optimized = 1;
    if (lgbend->verbose) {
      printf("LGBEND %s#%ld optimized: FSE[0]=%le, FSE[%ld]=%le, accuracy=%le\n",
             eptr?eptr->name:"?", eptr?eptr->occurence:-1, 
             lgbend->segment[0].fse, lgbend->nSegments-1, lgbend->segment[lgbend->nSegments-1].fse, acc);
      fflush(stdout);
    }
  }

  rad_coef = isr_coef = 0;

  if (lgbend->synch_rad)
    rad_coef = sqr(particleCharge)*pow3(Po)/(6*PI*epsilon_o*sqr(c_mks)*particleMass); 
  isr_coef = particleRadius*sqrt(55.0/(24*sqrt(3))*pow5(Po)*137.0359895);
  if (!lgbend->isr || (lgbend->isr1Particle==0 && n_part==1))
    /* minus sign indicates we accumulate into sigmadelta^2 only, don't perturb particles */
    isr_coef *= -1;
  if (lgbend->length<1e-6 && (lgbend->isr || lgbend->synch_rad)) {
    rad_coef = isr_coef = 0;  /* avoid unphysical results */
  }

  if (!(coef = expansion_coefficients(0)))
    bombTracking("expansion_coefficients(0) returned NULL pointer (track_through_lgbend)");
  if (!(coef = expansion_coefficients(1)))
    bombTracking("expansion_coefficients(1) returned NULL pointer (track_through_lgbend)");
  if (!(coef = expansion_coefficients(2)))
    bombTracking("expansion_coefficients(2) returned NULL pointer (track_through_lgbend)");

  nSlices = lgbend->nSlices;
  integ_order = lgbend->integration_order;
  if (nSlices<=0)
    bombTracking("nSlices<=0 in track_lgbend()");
  if (integ_order!=2 && integ_order!=4 && integ_order!=6) 
    bombTracking("multipole integration_order must be 2, 4, or 6");

  fse = lgbend->fse;

  exactDrift(particle, n_part, lgbend->predrift);

  for (iSegment=0; iSegment<lgbend->nSegments; iSegment++) {
    for (iTerm=0; iTerm<9; iTerm++)
      KnL[iTerm] = 0;
    angle = lgbend->segment[iSegment].angle;
    length = lgbend->segment[iSegment].length;
    entryAngle = lgbend->segment[iSegment].entryAngle;
    entryPosition = lgbend->segment[iSegment].entryX;
    exitAngle = lgbend->segment[iSegment].exitAngle;
    exitPosition = lgbend->segment[iSegment].exitX;
    rho0 = length/(sin(entryAngle) + sin(angle-entryAngle));
    KnL[0] = (1+fse+lgbend->segment[iSegment].fse)/rho0*length;
    KnL[1] = (1+fse+lgbend->segment[iSegment].fse)*lgbend->segment[iSegment].K1*length/(1-lgbend->segment[iSegment].KnDelta);
    KnL[2] = (1+fse+lgbend->segment[iSegment].fse)*lgbend->segment[iSegment].K2*length/(1-lgbend->segment[iSegment].KnDelta);
#ifdef DEBUG
    /*
    printf("segment %ld: angle=%le, length=%le, rho0=%le, entryX=%le, entryAngle=%le, exitX=%le, exitAngle=%le, K1L=%le, K2L=%le\n",
           iSegment, angle, length, rho0, entryPosition, entryAngle, exitPosition, exitAngle, KnL[1], KnL[2]);
    printf("             K1 = %le, K2 = %le, FSE = %le, FSE1 = %le, KnDelta = %le\n",
           lgbend->segment[iSegment].K1, lgbend->segment[iSegment].K2, fse, lgbend->segment[iSegment].fse,
           lgbend->segment[iSegment].KnDelta);
    */
#endif
    /*
    KnL[3] = (1+fse)*lgbend->K3*length/(1-lgbend->KnDelta);
    KnL[4] = (1+fse)*lgbend->K4*length/(1-lgbend->KnDelta);
    KnL[5] = (1+fse)*lgbend->K5*length/(1-lgbend->KnDelta);
    KnL[6] = (1+fse)*lgbend->K6*length/(1-lgbend->KnDelta);
    KnL[7] = (1+fse)*lgbend->K7*length/(1-lgbend->KnDelta);
    KnL[8] = (1+fse)*lgbend->K8*length/(1-lgbend->KnDelta);
    */
    if (angle<0) {
      angleSign = -1;
      for (iTerm=0; iTerm<9; iTerm+=2)
        KnL[iTerm] *= -1;
      angle = -angle;
      rho0 = -rho0;
      entryAngle = -entryAngle;
      exitAngle = -exitAngle;
      extraTilt = PI;
      entryPosition = -entryPosition;
    } else {
      angleSign = 1;
      extraTilt = 0;
    }
    
    tilt = lgbend->tilt + extraTilt;
    etilt = lgbend->etilt;
    dx = lgbend->dx;
    dy = lgbend->dy;
    dz = lgbend->dz*(lgbend->edgeFlip?-1:1);
    
    if (tilt) {
      /* this is needed because the DX and DY offsets will be applied after the particles are
       * tilted into the magnet reference frame
       */
      dx =  lgbend->dx*cos(tilt) + lgbend->dy*sin(tilt);
      dy = -lgbend->dx*sin(tilt) + lgbend->dy*cos(tilt);
    }

    setupMultApertureData(&apertureData, -tilt, apContour, maxamp, apFileData, z_start+length/2);

#ifdef DEBUG
    if (lgbend->optimized!=-1)
      fprintf(fpDeb, "%ld %le %le %le %le %le %le segstart\n",
              iSegment, (iSegment==0?0:lgbend->segment[iSegment-1].zAccumulated), particle[0][4], 
              particle[0][0], particle[0][1], 
              particle[0][2], particle[0][3]);
#endif

    if (iPart<=0) {
      if (tilt)
        rotateBeamCoordinatesForMisalignment(particle, n_part, tilt);
      if (iSegment==0) {
        switchLgbendPlane(particle, n_part, entryPosition, entryAngle, Po);
#ifdef DEBUG
        if (lgbend->optimized!=-1)
          fprintf(fpDeb, "%ld %le %le %le %le %le %le switch-lgbend-plane\n",
                  iSegment, (iSegment==0?0:lgbend->segment[iSegment-1].zAccumulated), particle[0][4], 
                  particle[0][0], particle[0][1], 
                  particle[0][2], particle[0][3]);
#endif
        if (dx || dy || dz)
          offsetBeamCoordinatesForMisalignment(particle, n_part, dx, dy, dz);
        if (etilt)
          rotateBeamCoordinatesForMisalignment(particle, n_part, etilt);
      }
      if (lgbend->segment[iSegment].has1) {
        lgbendFringe(particle, n_part, entryAngle, rho0, KnL[1]/length, KnL[2]/length,
                     lgbend->segment+iSegment, angleSign, 0);
#ifdef DEBUG
        if (lgbend->optimized!=-1)
          fprintf(fpDeb, "%ld %le %le %le %le %le %le fringe1\n",
                  iSegment, (iSegment==0?0:lgbend->segment[iSegment-1].zAccumulated), particle[0][4], 
                  particle[0][0], particle[0][1], 
                  particle[0][2], particle[0][3]);
#endif
      }
    }

    nTerms = 0;
    for (iTerm=8; iTerm>=0; iTerm--)
      if (KnL[iTerm]) {
        nTerms = iTerm+1;
        break;
      }
    
    if (sigmaDelta2)
      *sigmaDelta2 = 0;
    i_top = n_part-1;
    for (i_part=0; i_part<=i_top; i_part++) {
      if (!integrate_kick_KnL(particle[i_part], dx, dy, Po, rad_coef, isr_coef, KnL, nTerms,
                              integ_order, nSlices, iPart, iFinalSlice, length, NULL, NULL, NULL,
                              &apertureData, &dzLoss, sigmaDelta2)) {
        swapParticles(particle[i_part], particle[i_top]);
        if (accepted)
          swapParticles(accepted[i_part], accepted[i_top]);
        particle[i_top][4] = z_start+dzLoss;
        particle[i_top][5] = Po*(1+particle[i_top][5]);
        i_top--;
        i_part--;
        continue;
      }
    }
#ifdef DEBUG
    if (lgbend->optimized!=-1)
      fprintf(fpDeb, "%ld %le %le %le %le %le %le integrated\n",
              iSegment, lgbend->segment[iSegment].zAccumulated, particle[0][4], 
              particle[0][0], particle[0][1], 
              particle[0][2], particle[0][3]);
#endif
    multipoleKicksDone += (i_top+1)*nSlices;
    if (sigmaDelta2)
      *sigmaDelta2 /= i_top+1;
    
    if ((iPart<0 || iPart==(lgbend->nSlices-1)) && iFinalSlice<=0) {
        /*
          printf("output before adjustments: %16.10le %16.10le %16.10le %16.10le %16.10le %16.10le\n",
          particle[0][0], particle[0][1], particle[0][2],
          particle[0][3], particle[0][4], particle[0][5]);
        */
      if (lgbend->segment[iSegment].has2)  {
        lgbendFringe(particle, i_top+1, -exitAngle, rho0, KnL[1]/length, KnL[2]/length,
                     lgbend->segment+iSegment, angleSign, 1);
#ifdef DEBUG
        if (lgbend->optimized!=-1)
          fprintf(fpDeb, "%ld %le %le %le %le %le %le fringe2\n",
                  iSegment, lgbend->segment[iSegment].zAccumulated, particle[0][4], 
                  particle[0][0], particle[0][1], 
                  particle[0][2], particle[0][3]);
#endif
      }
      if (iSegment==(lgbend->nSegments-1)) {
        if (etilt)
          rotateBeamCoordinatesForMisalignment(particle, n_part, -etilt);
        if (dx || dy || dz)
          offsetBeamCoordinatesForMisalignment(particle, i_top+1, -dx, -dy, -dz);
        switchLgbendPlane(particle, i_top+1, -exitPosition, -exitAngle, Po);
        exactDrift(particle, i_top+1, lgbend->postdrift);
#ifdef DEBUG
        if (lgbend->optimized!=-1)
          fprintf(fpDeb, "%ld %le %le %le %le %le %le switch-lgbend-plane\n",
                  iSegment, lgbend->segment[iSegment].zAccumulated, particle[0][4], 
                  particle[0][0], particle[0][1], 
                  particle[0][2], particle[0][3]);
#endif
        if (tilt)
          /* use n_part here so lost particles get rotated back */
          rotateBeamCoordinatesForMisalignment(particle, n_part, -tilt);
        /*
          if (lgbend->optimized) {
          for (i_part=0; i_part<=i_top; i_part++)
          particle[i_part][4] += lgbend->lengthCorrection;
          }
        */
        /*
          printf("output after adjustments: %16.10le %16.10le %16.10le %16.10le %16.10le %16.10le\n",
          particle[0][0], particle[0][1], particle[0][2],
          particle[0][3], particle[0][4], particle[0][5]);
          fflush(stdout);
        */
      }
    } else if (iFinalSlice>0) {
      if (tilt)
        /* use n_part here so lost particles get rotated back */
        rotateBeamCoordinatesForMisalignment(particle, n_part, -tilt);
    }

    if (angleSign<0) {
      lastRho *= -1;
      lastX *= -1;
      lastXp *= -1;
    }
  }

  log_exit("track_through_lgbend");

  return(i_top+1);
}


#ifdef COMPILE_THIS
VMATRIX *determinePartialLgbendLinearMatrix(LGBEND *lgbend, double *startingCoord, double pCentral, long iFinalSlice)
/* This routine is used for getting the linear transport matrix from the start of the LGBEND to some interior slice.
 * We need this to compute radiation integrals.
 */
{
  double **coord;
  long n_track, i, j;
  VMATRIX *M;
  double **R, *C;
  double stepSize[6] = {1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5};
  long ltmp1, ltmp2;
  double angle0[2];

#if USE_MPI
  long notSinglePart_saved = notSinglePart;

  /* All the particles should do the same thing for this routine. */	
  notSinglePart = 0;
#endif
   		 
  coord = (double**)czarray_2d(sizeof(**coord), 1+6*4, totalPropertiesPerParticle);

  n_track = 4*6+1;
  for (j=0; j<6; j++)
    for (i=0; i<n_track; i++)
      coord[i][j] = startingCoord ? startingCoord[j] : 0;

  /* particles 0 and 1 are for d/dx */
  coord[0][0] += stepSize[0] ;
  coord[1][0] -= stepSize[0] ;
  /* particles 2 and 3 are for d/dxp */
  coord[2][1] += stepSize[1];
  coord[3][1] -= stepSize[1];
  /* particles 4 and 5 are for d/dy */
  coord[4][2] += stepSize[2] ;
  coord[5][2] -= stepSize[2] ;
  /* particles 6 and 7 are for d/dyp */
  coord[6][3] += stepSize[3] ;
  coord[7][3] -= stepSize[3] ;
  /* particles 8 and 9 are for d/ds */
  coord[8][4] += stepSize[4] ;
  coord[9][4] -= stepSize[4] ;
  /* particles 10 and 11 are for d/delta */
  coord[10][5] += stepSize[5];
  coord[11][5] -= stepSize[5];

  /* particles 12 and 13 are for d/dx */
  coord[12][0] += 3*stepSize[0] ;
  coord[13][0] -= 3*stepSize[0] ;
  /* particles 14 and 15 are for d/dxp */
  coord[14][1] += 3*stepSize[1];
  coord[15][1] -= 3*stepSize[1];
  /* particles 16 and 17 are for d/dy */
  coord[16][2] += 3*stepSize[2] ;
  coord[17][2] -= 3*stepSize[2] ;
  /* particles 18 and 19 are for d/dyp */
  coord[18][3] += 3*stepSize[3] ;
  coord[19][3] -= 3*stepSize[3] ;
  /* particles 20 and 21 are for d/ds */
  coord[20][4] += 3*stepSize[4] ;
  coord[21][4] -= 3*stepSize[4] ;
  /* particles 22 and 23 are for d/delta */
  coord[22][5] += 3*stepSize[5];
  coord[23][5] -= 3*stepSize[5];
  /* particle n_track-1 is the reference particle (coordinates set above) */

  /* save radiation-related parameters */
  ltmp1 = lgbend->isr;
  ltmp2 = lgbend->synch_rad;

  lgbend->isr = lgbend->synch_rad = 0;
  track_through_lgbend(coord, n_track, NULL, lgbend, pCentral, NULL, 0.0, NULL, NULL, NULL, NULL, NULL, -1, iFinalSlice);

  lgbend->isr = ltmp1;
  lgbend->synch_rad = ltmp2;
  
  M = tmalloc(sizeof(*M));
  M->order = 1;
  initialize_matrices(M, M->order);
  R = M->R;
  C = M->C;

  /* Rotate the coordinates into the local frame by assuming that the reference particle
   * is on the reference trajectory. Not valid if errors are present, but not a bad approximation presumably.
   */
  for (j=0; j<2; j++)
    angle0[j] = atan(coord[n_track-1][2*j+1]);
  for (i=0; i<n_track; i++) {
    for (j=0; j<2; j++)
      coord[i][2*j] = coord[i][2*j] - coord[n_track-1][2*j];
    for (j=0; j<2; j++)
      coord[i][2*j+1] = tan(atan(coord[i][2*j+1])-angle0[j]);
    coord[i][4] -= coord[n_track-1][4];
  }
  
  for (i=0; i<6; i++) {
    /* i indexes the dependent quantity */

    /* Determine C[i] */
    C[i] = coord[n_track-1][i];

    /* Compute R[i][j] */
    for (j=0; j<6; j++) {
      /* j indexes the initial coordinate value */
      R[i][j] = 
        (27*(coord[2*j][i]-coord[2*j+1][i])-(coord[2*j+12][i]-coord[2*j+13][i]))/(48*stepSize[j]);
    }
  }

  free_czarray_2d((void**)coord, 1+4*6, totalPropertiesPerParticle);

#if USE_MPI
  notSinglePart = notSinglePart_saved;
#endif
  return M;
}
#endif

#undef DEBUG

#ifdef COMPILE_THIS

void addLgbendRadiationIntegrals(LGBEND *lgbend, double *startingCoord, double pCentral,
                                 double eta0, double etap0, double beta0, double alpha0,
                                 double *I1, double *I2, double *I3, double *I4, double *I5, ELEMENT_LIST *elem)
{
  long iSlice;
  VMATRIX *M;
  double gamma0, K1;
  double eta1, beta1, alpha1, etap1;
  double eta2, beta2, alpha2, etap2;
  double C, S, Cp, Sp, ds, H1, H2;
#ifdef DEBUG
  double s0;
  static FILE *fpcr = NULL;
  if (fpcr==NULL) {
    fpcr = fopen("lgbend-RI.sdds", "w");
    fprintf(fpcr, "SDDS1\n&column name=Slice type=short &end\n&column name=s type=double units=m &end\n");
    fprintf(fpcr, "&column name=betax type=double units=m &end\n");
    fprintf(fpcr, "&column name=etax type=double units=m &end\n");
    fprintf(fpcr, "&column name=etaxp type=double &end\n");
    fprintf(fpcr, "&column name=alphax type=double &end\n");
    fprintf(fpcr, "&column name=ElementName type=string &end\n");
    fprintf(fpcr, "&data mode=ascii no_row_counts=1 &end\n");
  }
  s0 = elem->end_pos - lgbend->length;
  fprintf(fpcr, "0 %le %le %le %le %le %s\n", s0, beta0, eta0, etap0, alpha0, elem->name);
#endif
  
  if (lgbend->tilt)
    bombElegant("Can't add radiation integrals for tilted LGBEND\n", NULL);

  gamma0 = (1+alpha0*alpha0)/beta0;

  beta1 = beta0;
  eta1 = eta0;
  etap1 = etap0;
  alpha1 = alpha0;
  H1 = (eta1*eta1 + sqr(beta1*etap1 + alpha1*eta1))/beta1;
  ds = lgbend->length/lgbend->nSlices; /* not really right... */
  for (iSlice=1; iSlice<=lgbend->nSlices; iSlice++) {
    /* Determine matrix from start of element to exit of slice iSlice */
    M = determinePartialLgbendLinearMatrix(lgbend, startingCoord, pCentral, iSlice);
    C = M->R[0][0];
    Cp = M->R[1][0];
    S = M->R[0][1];
    Sp = M->R[1][1];

    /* propagate lattice functions from start of element to exit of this slice */
    beta2 = (sqr(C)*beta0 - 2*C*S*alpha0 + sqr(S)*gamma0);
    alpha2 = -C*Cp*beta0 + (Sp*C+S*Cp)*alpha0 - S*Sp*gamma0;
    eta2 = C*eta0 + S*etap0 + M->R[0][5];
    etap2 = Cp*eta0 + Sp*etap0 + M->R[1][5];
#ifdef DEBUG
    s0 += ds;
    fprintf(fpcr, "%ld %le %le %le %le %le %s\n", iSlice, s0, beta2, eta2, etap2, alpha2, elem->name);
#endif

    /* Compute contributions to radiation integrals in this slice.
     * lastRho is saved by the routine track_through_lgbend(). Since determinePartialLgbendLinearMatrix()
     * puts the reference particle first, it is appropriate to the central trajectory. 
     */
    *I1 += ds*(eta1+eta2)/2/lastRho;
    *I2 += ds/sqr(lastRho);
    *I3 += ds/ipow(fabs(lastRho), 3);
    /* Compute effective K1 including the sextupole effect plus rotation.
     * lastX and lastXp are saved by track_through_lgbend().
     */
    K1 = (lgbend->K1+(lastX-lgbend->dxOffset)*lgbend->K2)*cos(atan(lastXp));
    *I4 += ds*(eta1+eta2)/2*(1/ipow(lastRho, 3) + 2*K1/lastRho); 
    H2 = (eta2*eta2 + sqr(beta2*etap2 + alpha2*eta2))/beta2;
    *I5 += ds/ipow(fabs(lastRho), 3)*(H1+H2)/2;

    /* Save lattice functions as values at start of next slice */
    beta1 = beta2;
    alpha1 = alpha2;
    eta1 = eta2;
    etap1 = etap2;
    H1 = H2;
    free_matrices(M);
    free(M);
  }
}
#endif

void readLGBendConfiguration(LGBEND *lgbend, ELEMENT_LIST *eptr)
{
  SDDS_DATASET SDDSin;
  long readCode, i, index;
  short postdriftSeen, predriftSeen;
  uint64_t iRow, rows;
  double *parameterValue;
  char **parameterName;
  /* These names need to have the same order as the items in LGBEND_SEGMENT */
#define NLGBEND (8+2*N_LGBEND_FRINGE_INT+2)
  static char *knownParameterName[NLGBEND] = {
    "LONGIT_L", "K1", "K2", "ANGLE", "ENTRY_X", "ENTRY_ANGLE", "EXIT_X", "EXIT_ANGLE", 
    "FRINGE1K0","FRINGE1I0", "FRINGE1K2","FRINGE1I1", "FRINGE1K4", "FRINGE1K5", "FRINGE1K6", "FRINGE1K7",
    "FRINGE2K0","FRINGE2I0", "FRINGE2K2","FRINGE2I1", "FRINGE2K4", "FRINGE2K5", "FRINGE2K6", "FRINGE2K7",
    "PREDRIFT", "POSTDRIFT",
  };
  short required[NLGBEND] = {
    2, 0, 0, 1, 1, 1, 2, 2, 
    0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 
    1, 0
  };
  short disallowed[NLGBEND] = {
    0, 0, 0, 0, 0, 0, 0, 0, 
    1, 1, 1, 1, 1, 1, 1, 1, 
    0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0
  };
  short provided[NLGBEND];
  short offset[NLGBEND];
  LGBEND_SEGMENT segmentExample;

  offset[0] = 0;
  offset[1] = (short)(&(segmentExample.K1) - &(segmentExample.length));
  offset[2] = (short)(&(segmentExample.K2) -  &(segmentExample.length));
  offset[3] = (short)(&(segmentExample.angle) -  &(segmentExample.length));
  offset[4] = (short)(&(segmentExample.entryX) -  &(segmentExample.length));
  offset[5] = (short)(&(segmentExample.entryAngle) -  &(segmentExample.length));
  offset[6] = (short)(&(segmentExample.exitX) -  &(segmentExample.length));
  offset[7] = (short)(&(segmentExample.exitAngle) -  &(segmentExample.length));
  offset[8] = (short)(&(segmentExample.fringeInt1K0) -  &(segmentExample.length));
  offset[9] = (short)(&(segmentExample.fringeInt1I0) -  &(segmentExample.length));
  offset[10] = (short)(&(segmentExample.fringeInt1K2) -  &(segmentExample.length));
  offset[11] = (short)(&(segmentExample.fringeInt1I1) -  &(segmentExample.length));
  offset[12] = (short)(&(segmentExample.fringeInt1K4) -  &(segmentExample.length));
  offset[13] = (short)(&(segmentExample.fringeInt1K5) -  &(segmentExample.length));
  offset[14] = (short)(&(segmentExample.fringeInt1K6) -  &(segmentExample.length));
  offset[15] = (short)(&(segmentExample.fringeInt1K7) -  &(segmentExample.length));
  offset[16] = (short)(&(segmentExample.fringeInt2K0) -  &(segmentExample.length));
  offset[17] = (short)(&(segmentExample.fringeInt2I0) -  &(segmentExample.length));
  offset[18] = (short)(&(segmentExample.fringeInt2K2) -  &(segmentExample.length));
  offset[19] = (short)(&(segmentExample.fringeInt2I1) -  &(segmentExample.length));
  offset[20] = (short)(&(segmentExample.fringeInt2K4) -  &(segmentExample.length));
  offset[21] = (short)(&(segmentExample.fringeInt2K5) -  &(segmentExample.length));
  offset[22] = (short)(&(segmentExample.fringeInt2K6) -  &(segmentExample.length));
  offset[23] = (short)(&(segmentExample.fringeInt2K7) -  &(segmentExample.length));
  offset[24] = -1;
  offset[25] = -1;
  for (i=0; i<NLGBEND; i++)
    offset[i] *= sizeof(double);

  if (!SDDS_InitializeInputFromSearchPath(&SDDSin, lgbend->configuration)) {
    fprintf(stderr, "Error: unable to find or read LGBEND configuration file %s\n", lgbend->configuration);
    exitElegant(1);
  }
  if (SDDS_CheckColumn(&SDDSin, "ParameterName", NULL, SDDS_STRING, stdout)!=SDDS_CHECK_OK ||
      SDDS_CheckColumn(&SDDSin, "ParameterValue", NULL, SDDS_ANY_FLOATING_TYPE, stdout)!=SDDS_CHECK_OK) {
    fprintf(stderr, "Error: required columns with required types not present in LGBEND configuration file %s\n", 
            lgbend->configuration);
    exitElegant(1);
  }

  lgbend->angle = lgbend->length = 0;
  postdriftSeen = 0;
  while ((readCode=SDDS_ReadPage(&SDDSin))>=1) {
    if (postdriftSeen) {
      fprintf(stderr, "Error: POSTDRIFT parameter seen for interior segment reading LGBEND configuration file %s\n", 
              lgbend->configuration);
      exitElegant(1);
    }
    if (!(lgbend->segment = SDDS_Realloc(lgbend->segment, (lgbend->nSegments+1)*sizeof(*(lgbend->segment))))) {
      fprintf(stderr, "Error: memory allocation failure reading LGBEND configuration file %s\n", lgbend->configuration);
      exitElegant(1);
    }
    memset(&(lgbend->segment[lgbend->nSegments]), 0, sizeof(*(lgbend->segment)));
    parameterName = NULL;
    parameterValue = NULL;
    if ((rows = SDDS_RowCount(&SDDSin))<=0 ||
        !(parameterName=SDDS_GetColumn(&SDDSin, "ParameterName")) ||
        !(parameterValue=SDDS_GetColumn(&SDDSin, "ParameterValue"))) {
      fprintf(stderr, "Error: problem getting data from LGBEND configuration file %s\n", lgbend->configuration);
      exitElegant(1);
    }
    memset(provided, 0, sizeof(provided[0])*(NLGBEND));
    postdriftSeen = predriftSeen = 0;
    for (iRow=0; iRow<rows; iRow++) {
      if ((index=match_string(parameterName[iRow], knownParameterName, NLGBEND, EXACT_MATCH))<0) {
        fprintf(stdout, "Error: unrecognized parameter name \"%s\" in LGBEND configuration file %s (page %ld, row %ld)\n", 
                parameterName[iRow], lgbend->configuration,
                readCode, iRow);
        exitElegant(1);
      }
      provided[index] = 1;
      if (index==24) {
        predriftSeen += 1;
        lgbend->predrift = parameterValue[iRow];
      } else if (index==25) {
        postdriftSeen += 1;
        lgbend->postdrift = parameterValue[iRow];
      } else
        *((double*)(((char*)&(lgbend->segment[lgbend->nSegments]))+offset[index])) = parameterValue[iRow];
    }
    free(parameterValue);
    SDDS_FreeStringArray(parameterName, rows);
    if (predriftSeen && readCode!=1) {
      fprintf(stderr, "Error: PREDRIFT parameter seen for interior segment %ld reading LGBEND configuration file %s\n", 
              readCode, lgbend->configuration);
      exitElegant(1);
    }
    if (!predriftSeen && readCode==1) {
      fprintf(stderr, "Error: PREDRIFT parameter not seen for initial segment reading LGBEND configuration file %s\n", 
              lgbend->configuration);
      exitElegant(1);
    }
    for (i=0; i<NLGBEND; i++)  {
      if ((readCode==1 && required[i] && !provided[i]) ||
          (readCode>1 && required[i]==2 && !provided[i])) {
        fprintf(stdout, "Error: parameter %s is required but was not provided in page %ld of configuration file %s for LGBEND %s#%ld\n",
                knownParameterName[i], readCode, lgbend->configuration, eptr->name, eptr->occurence);
        exitElegant(1);
      }
      if (readCode>1 && provided[i] && disallowed[i]) {
        fprintf(stdout, "Error: parameter %s is provided in page %ld of configuration file %s for LGBEND %s#%ld. Not meaningful.\n",
                knownParameterName[i], readCode, lgbend->configuration, eptr->name, eptr->occurence);
        exitElegant(1);
      }
    }
    lgbend->segment[lgbend->nSegments].has2 = 1;
    lgbend->segment[lgbend->nSegments].has1 = readCode==1?1:0;
    if (lgbend->nSegments>0)  {
      lgbend->segment[lgbend->nSegments].entryX = lgbend->segment[lgbend->nSegments-1].exitX;
      lgbend->segment[lgbend->nSegments].entryAngle = lgbend->segment[lgbend->nSegments-1].exitAngle;
    }
    lgbend->angle += lgbend->segment[lgbend->nSegments].angle;
    lgbend->length += lgbend->segment[lgbend->nSegments].length;
    lgbend->segment[lgbend->nSegments].zAccumulated = lgbend->length;
    lgbend->nSegments ++;
  }
  SDDS_Terminate(&SDDSin);
  lgbend->initialized = 1;
  if (!postdriftSeen) {
    fprintf(stderr, "Error: POSTDRIFT parameter not seen for final segment reading LGBEND configuration file %s\n", 
            lgbend->configuration);
    exitElegant(1);
  }

  if (lgbend->edgeFlip) {
    flipLGBEND(lgbend);
    lgbend->edgeFlip = 0;
  }

  if (lgbend->verbose>10) {
    printf("%ld segments found for LGBEND %s#%ld in file %s\n",
           lgbend->nSegments, eptr->name, eptr->occurence, lgbend->configuration);
    printf("total angle is %le, total length is %le\n", lgbend->angle, lgbend->length);
    fflush(stdout);
    for (i=0; i<lgbend->nSegments; i++) {
      printf("Segment %ld\nLONGIT_L = %le\nK1 = %le\nK2 = %le\n", 
             i, 
             lgbend->segment[i].length, lgbend->segment[i].K1, lgbend->segment[i].K2);
      printf("ANGLE = %le\nENTRY_X = %le\nENTRY_ANGLE = %le\nEXIT_X = %le\nEXIT_ANGLE = %le\n",
             lgbend->segment[i].angle,
             lgbend->segment[i].entryX, lgbend->segment[i].entryAngle,
             lgbend->segment[i].exitX, lgbend->segment[i].exitAngle);
      if (lgbend->segment[i].has1)
        printf("FRINGE1K0 = %le\nFRINGE1I0 = %le\nFRINGE1K2 = %le\nFRINGE1I1 = %le\nFRINGE1K4 = %le\nFRINGE1K5 = %le\nFRINGE1K6 = %le\nFRINGE1K7 = %le\n",
               lgbend->segment[i].fringeInt1K0,
               lgbend->segment[i].fringeInt1I0,
               lgbend->segment[i].fringeInt1K2,
               lgbend->segment[i].fringeInt1I1,
               lgbend->segment[i].fringeInt1K4,
               lgbend->segment[i].fringeInt1K5,
               lgbend->segment[i].fringeInt1K6,
               lgbend->segment[i].fringeInt1K7);
      if (lgbend->segment[i].has2)
        printf("FRINGE2K0 = %le\nFRINGE2I0 = %le\nFRINGE2K2 = %le\nFRINGE2I1 = %le\nFRINGE2K4 = %le\nFRINGE2K5 = %le\nFRINGE2K6 = %le\nFRINGE2K7 = %le\n",
               lgbend->segment[i].fringeInt2K0,
               lgbend->segment[i].fringeInt2I0,
               lgbend->segment[i].fringeInt2K2,
               lgbend->segment[i].fringeInt2I1,
               lgbend->segment[i].fringeInt2K4,
               lgbend->segment[i].fringeInt2K5,
               lgbend->segment[i].fringeInt2K6,
               lgbend->segment[i].fringeInt2K7);
    }
  }

}

void flipLGBEND(LGBEND *lgbend)
{
  long i, j;
  LGBEND_SEGMENT *segment;

  SWAP_DOUBLE(lgbend->xEntry, lgbend->xExit);
  lgbend->zEntry *= -1;
  lgbend->zExit *= -1;
  SWAP_DOUBLE(lgbend->zEntry, lgbend->zExit);
  lgbend->zVertex *= -1;

  /* need to reverse the order of the segments and exchange edge integrals */
  /* - copy the data to temporary memory */
  segment = tmalloc(sizeof(*segment)*lgbend->nSegments);
  for (i=0; i<lgbend->nSegments; i++)
    memcpy(segment+i, lgbend->segment+i, sizeof(*segment));
  /* - copy back to caller's memory, but in reversed order */
  for (i=0, j=lgbend->nSegments-1; i<lgbend->nSegments; i++, j--) {
    if (i!=j)
      memcpy(lgbend->segment+i, segment+j, sizeof(*segment));
  }
  free(segment);
  /* - swap edge data and adjust signs as needed */
  for (i=0; i<lgbend->nSegments; i++) {
    SWAP_DOUBLE(lgbend->segment[i].fringeInt1K0, lgbend->segment[i].fringeInt2K0);
    SWAP_DOUBLE(lgbend->segment[i].fringeInt1I0, lgbend->segment[i].fringeInt2I0);
    SWAP_DOUBLE(lgbend->segment[i].fringeInt1K2, lgbend->segment[i].fringeInt2K2);
    SWAP_DOUBLE(lgbend->segment[i].fringeInt1I1, lgbend->segment[i].fringeInt2I1);
    SWAP_DOUBLE(lgbend->segment[i].fringeInt1K4, lgbend->segment[i].fringeInt2K4);
    SWAP_DOUBLE(lgbend->segment[i].fringeInt1K5, lgbend->segment[i].fringeInt2K5);
    SWAP_DOUBLE(lgbend->segment[i].fringeInt1K6, lgbend->segment[i].fringeInt2K6);
    SWAP_DOUBLE(lgbend->segment[i].fringeInt1K7, lgbend->segment[i].fringeInt2K7);
    lgbend->segment[i].fringeInt1K0 *= -1;
    lgbend->segment[i].fringeInt1I1 *= -1;
    lgbend->segment[i].fringeInt1K5 *= -1;
    lgbend->segment[i].fringeInt2K0 *= -1;
    lgbend->segment[i].fringeInt2I1 *= -1;
    lgbend->segment[i].fringeInt2K5 *= -1;
    SWAP_SHORT(lgbend->segment[i].has1, lgbend->segment[i].has2);
    SWAP_DOUBLE(lgbend->segment[i].entryX, lgbend->segment[i].exitX);
    if (i==0) {
      lgbend->segment[i].entryAngle = -lgbend->segment[i].exitAngle;
      lgbend->segment[i].exitAngle = lgbend->segment[i].entryAngle - lgbend->segment[i].angle;
    }
  }
  for (i=1; i<lgbend->nSegments; i++) {
    lgbend->segment[i].entryAngle = lgbend->segment[i-1].exitAngle;
    lgbend->segment[i].exitAngle = lgbend->segment[i].entryAngle - lgbend->segment[i].angle;
  }
  SWAP_DOUBLE(lgbend->predrift, lgbend->postdrift);
}

void lgbendFringe
(
 double **particle, long n_part, 
 double alpha, // edge angle relative to beam path
 double rho0,  // bending radius (always positive)
 double K1,    // interior gradient
 double K2,    // interior sextupole
 LGBEND_SEGMENT *segment,
 short angleSign,      // -1 or 1
 short isExit
 )
{
  double x1, px1, y1, py1, tau1, delta;
  double x2, px2, y2, py2, tau2, tau;
  
  double intK0, intK2, intK4, intK5, intK6, intK7, intI0, intI1;
  double invRhoPlus, invRhoMinus, K1plus, K1minus;
  double tant, sect, sect3, sint, temp;
  double focX0, focXd, focY0, focYd, invP;
  
  long i;
  if (isExit) {
    // exit fringe
    alpha = -alpha;
    intK0 = segment->fringeInt2K0;
    intI0 = segment->fringeInt2I0;
    intK2 = segment->fringeInt2K2;
    intI1 = segment->fringeInt2I1;
    intK4 = segment->fringeInt2K4;
    intK5 = segment->fringeInt2K5;
    intK6 = segment->fringeInt2K6;
    intK7 = segment->fringeInt2K7;
    K1minus = K1;
    K1plus = 0;
    invRhoMinus = 1/rho0;
    invRhoPlus = 0;
  } else {
    // entry fringe
    intK0 = segment->fringeInt1K0;
    intI0 = segment->fringeInt1I0;
    intK2 = segment->fringeInt1K2;
    intI1 = segment->fringeInt1I1;
    intK4 = segment->fringeInt1K4;
    intK5 = segment->fringeInt1K5;
    intK6 = segment->fringeInt1K6;
    intK7 = segment->fringeInt1K7;
    K1plus = K1;
    K1minus = 0;
    invRhoPlus = 1/rho0;
    invRhoMinus = 0;
  }

  // Per R. Lindberg, some of the fringe integrals change sign with the sign of the bending angle.
  intK0 *= angleSign;
  intK4 *= angleSign;
  intK5 *= angleSign;
  intK6 *= angleSign;
  
  tant = tan(alpha);
  sint = sin(alpha);
  sect = 1.0/cos(alpha);
  sect3 = sect*sect*sect;
  focX0 = -tant*intK5 - 0.5*intI0*(2.0 - tant*tant);
  focXd = sect3*intK7;
  focY0 = -tant*(invRhoPlus - invRhoMinus) + tant*sect*sect*intK5 + 0.5*intI0*(2.0 + tant*tant);
  focYd = sect3*((1.0+sint*sint)*intK2 - intK7);
  for (i=0; i<n_part; i++) {
    x1  = particle[i][0];
    temp = sqrt( 1.0 + sqr(particle[i][1]) + sqr(particle[i][3]) );
    px1 = (1.0 + particle[i][5])*particle[i][1]/temp;
    y1  = particle[i][2];
    py1 = (1.0 + particle[i][5])*particle[i][3]/temp;
    tau1 = -particle[i][4] + sint*x1;
    delta = particle[i][5];
    px1 = px1 - (1.0 + delta)*sint;
    invP = 1.0/(1.0 + delta);
      
    tau = tau1;
    temp = 1.0 - 0.5*sect3*intK5*x1*invP;
    x2  = x1/temp;
    y2  = y1;
    px2 = px1*temp*temp;
    py2 = py1;
    tau -= 0.5*sect3*intK5*x1*x1*px1*invP*invP;
      
    temp = sect*intK5*invP;
    x1  = x2;
    y1  = y2*exp(-temp*x2);
    px1 = px2 + temp*y2*py2;
    py1 = py2*exp(temp*x2);
    tau += temp*x2*y2*py2*invP;

    temp = 0.5*sect3*(intK5 - (invRhoPlus-invRhoMinus))*invP;
    x2  = x1 - temp*y1*y1;
    y2  = y1;
    px2 = px1;
    py2 = py1 + 2.0*temp*px1*y1;
    tau += temp*y1*y1*px1*invP;

    x1 = x2;
    y1 = y2;
    px1 = px2 + (focX0 + focXd*invP)*x2
      - 0.25*tant*(K1plus - K1minus)*(x2*x2 + y2*y2) - 0.5*intK6*(x2*x2 - sect*sect*y2*y2)
      + ( 0.5*(sect*intK5*focY0 + 0.5*sect3*(intK5 - (invRhoPlus-invRhoMinus))*focX0)*y2*y2
          - 0.75*sect3*intK5*focX0*x2*x2 )*invP;
    py1 = py2 + (focY0 + focYd*invP)*y2
      + (sect*sect*intK6 - 0.5*tant*(K1plus - K1minus))*x2*y2
      + ( sect*intK5*focY0 + 0.5*sect3*focX0*(intK5 - (invRhoPlus-invRhoMinus)) )*x2*y2*invP;

    tau += 0.5*( (focXd*x2*x2 + focYd*y2*y2) + sect*intK5*focY0*x2*y2*y2
                 + 0.5*sect3*( (intK5 - (invRhoPlus-invRhoMinus))*y2*y2 - intK5*x2*x2 )*focX0*x2 )*invP*invP;

    x2  = x1 - sect3*intK0*invP;
    y2  = y1;
    px2 = px1;
    py2 = py1;
    //      tau2 = tau + sect3*intK0*px2*invP*invP;
    tau += sect3*intK0*px2*invP*invP;

    temp = sect*intI1*invP;
    tau2 = tau + temp*invP*(py2*y2 - px2*x2);
    temp = exp(temp);
    // temp = exp(sect*intI1);
    x2 = x2*temp;
    y2 = y2/temp;
    px2 = px2/temp;
    py2 = py2*temp;	

    /* x2 = x1;
       y2 = y1;
       px2 = px1;
       py2 = py1;
       tau2 = tau1;

       x2 += 0.5*sect3*( intK5*(x1*x1 - y1*y1) + y1*y1*(invRhoPlus-invRhoMinus) )*invP;
       y2 += -sect*intK5*x1*y1*invP;
       px2 += -tant*intI1*0.0 + (focX0 + focXd*invP)*x1
       - 0.25*tant*(K1plus - K1minus)*(x1*x1 + y1*y1) - 0.5*intK6*(x1*x1 - sect*sect*y1*y1)
       + intK5*(sect*py1*y1 - sect3*px1*x1)*invP; // hereit
       py2 += (focY0 + focYd*invP)*y1
       + (sect*sect*intK6 - 0.5*tant*(K1plus - K1minus))*x1*y1
       + ( intK5*(sect*py1*x1 + sect3*px1*y1) - sect3*(invRhoPlus - invRhoMinus)*px1*y1 )*invP;
       tau2 += ( 0.5*(focXd*x1*x1 + focYd*y1*y1) + sect*intK5*py1*x1*y1
       + sect3*(intK5*px1*(y1*y1-x1*x1) - px1*y1*y1*(invRhoPlus-invRhoMinus)) )*invP*invP;
       x2 += -sect3*intK0*invP;
       tau2 += sect3*intK0*px2*invP*invP; */

    particle[i][0] = x2;
    particle[i][2] = y2;
    px2 = px2 + (1.0 + delta)*sint;
    temp = sqrt( sqr(1+delta) - sqr(px2) - sqr(py2) );
    particle[i][1] = px2/temp;
    particle[i][3] = py2/temp;
    particle[i][4] = -tau2 + sint*x2;
  }
}

void switchLgbendPlane(double **particle, long n_part, double dx, double alpha, double po)
/* transforms the reference plane to one that is at an angle alpha relative to the
 * initial plane. 
 */
{
  long i;
  double s, *coord, sin_alpha, cos_alpha, tan_alpha;
  double d, qx0, qy0, qz0, qx, qy, qz;

  tan_alpha = tan(alpha);
  cos_alpha = cos(alpha);
  sin_alpha = sin(alpha);

  for (i=0; i<n_part; i++) {
    coord = particle[i];
    s = coord[0]*tan_alpha/(1-coord[1]*tan_alpha);
    coord[0] = (coord[0]+coord[1]*s)/cos_alpha;
    coord[2] = (coord[2]+coord[3]*s);
    coord[4] += s;
    d = sqrt(1+sqr(coord[1])+sqr(coord[3]));
    qx0 = coord[1]*(1+coord[5])/d;
    qy0 = coord[3]*(1+coord[5])/d;
    qz0 = (1+coord[5])/d;
    qx =  qx0*cos_alpha + qz0*sin_alpha;
    qy = qy0;
    qz = -qx0*sin_alpha + qz0*cos_alpha;
    coord[1] = qx/qz;
    coord[3] = qy/qz;
    coord[0] += dx;
  }
}

double lgbend_trajectory_error(double *value, long *invalid)
{
  static double **particle = NULL;
  double result;

  *invalid = 0;
  if ((lgbendCopy.segment[0].fse=value[0])>=1 || value[0]<=-1) {
    *invalid = 1;
    return DBL_MAX;
  }
  if ((lgbendCopy.segment[lgbendCopy.nSegments-1].fse=value[1])>=1 || value[1]<=-1) {
    *invalid = 1;
    return DBL_MAX;
  }
  if (!particle) 
    particle = (double**)czarray_2d(sizeof(**particle), 1, totalPropertiesPerParticle);
  memset(particle[0], 0, totalPropertiesPerParticle*sizeof(**particle));
  if (lgbendCopy.compensateKn) {
    lgbendCopy.segment[0].KnDelta = -lgbendCopy.segment[0].fse;
    lgbendCopy.segment[lgbendCopy.nSegments-1].KnDelta = -lgbendCopy.segment[lgbendCopy.nSegments-1].fse;
  }
  if (!track_through_lgbend(particle, 1, eptrCopy, &lgbendCopy, PoCopy, NULL, 0.0, NULL, NULL, NULL, NULL, NULL, -1, -1)) {
    *invalid = 1;
    return DBL_MAX;
  }
  result = fabs(particle[0][0]) + fabs(particle[0][1]);
  return result;
}
