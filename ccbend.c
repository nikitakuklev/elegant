/************************************************************************* \
* Copyright (c) 2018 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2018 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* routine: crbend()
 * purpose: tracking through canonical bend in cartesian coordinates
 * 
 * Michael Borland, 2018
 */
#include "mdb.h"
#include "track.h"
#include "multipole.h"
#define DEBUG

static CCBEND ccbendCopy;
static double PoCopy, xMin, xFinal, xAve, xMax, xError, xpError, lastRho, lastX, lastXp;
#define OPTIMIZE_X  0x01UL
#define OPTIMIZE_XP 0x02UL
static unsigned long optimizationFlags = 0;
#ifdef DEBUG
static short logHamiltonian = 0;
static FILE *fpHam = NULL;
#endif

void switchRbendPlane(double **particle, long n_part, double alpha, double Po);
void verticalRbendFringe(double **particle, long n_part, double alpha, double rho0, double K1, double K2, double gK, long order);
int integrate_kick_K012(double *coord, double dx, double dy, 
                        double Po, double rad_coef, double isr_coef,
                        double K0L, double K1L, double K2L,
                        long integration_order, long n_parts, long iPart, long iFinalSlice,
                        double drift,
                        MULTIPOLE_DATA *multData, MULTIPOLE_DATA *edgeMultData, 
                        MULT_APERTURE_DATA *apData, double *dzLoss, double *sigmaDelta2);
int integrate_kick_KnL(double *coord, double dx, double dy, 
                      double Po, double rad_coef, double isr_coef,
                      double *KnL, long nTerms,
                      long integration_order, long n_parts, long iPart, long iFinalSlice,
                      double drift,
                      MULTIPOLE_DATA *multData, MULTIPOLE_DATA *edge1MultData, MULTIPOLE_DATA *edge2MultData, 
                      MULT_APERTURE_DATA *apData, double *dzLoss, double *sigmaDelta2);
double ccbend_trajectory_error(double *value, long *invalid);

long track_through_ccbend(
                          double **particle,   /* initial/final phase-space coordinates */
                          long n_part,         /* number of particles */
			  ELEMENT_LIST *eptr,
                          CCBEND *ccbend,
                          double Po,
                          double **accepted,
                          double z_start,
                          double *sigmaDelta2, /* use for accumulation of energy spread for radiation matrix computation */
                          char *rootname,
                          MAXAMP *maxamp,
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
  long iTerm, nTerms;
  double dx, dy, dz; /* offsets of the multipole center */
  long n_kicks, integ_order;
  long i_part, i_top;
  double *coef;
  double fse, tilt, rad_coef, isr_coef, dzLoss=0;
  double rho0, arcLength, length, angle, yaw, angleSign;
  MULTIPOLE_DATA *multData = NULL, *edge1MultData = NULL, *edge2MultData = NULL;
  long freeMultData=0;
  MULT_APERTURE_DATA apertureData;
  double referenceKnL = 0;
  double gK[2];

  if (!particle)
    bombTracking("particle array is null (track_through_ccbend)");

  if (iPart>=0 && ccbend->optimized!=1)
    bombTracking("Programming error: one-step mode invoked for unoptimized CCBEND.");
  if (iFinalSlice>0 && iPart>=0)
    bombTracking("Programming error: partial integration mode and one-step mode invoked together for CCBEND.");
  if (iFinalSlice>0 && ccbend->optimized!=1)
    bombTracking("Programming error: partial integration mode invoked for unoptimized CCBEND.");

  if ((ccbend->optimizeFse || ccbend->optimizeDx) && ccbend->optimized!=-1 && ccbend->angle!=0) {
    if (ccbend->optimized==0 || 
        ccbend->length!=ccbend->referenceData[0] ||
        ccbend->angle!=ccbend->referenceData[1] ||
        ccbend->K1!=ccbend->referenceData[2] ||
        ccbend->K2!=ccbend->referenceData[3] ||
        ccbend->yaw!=ccbend->referenceData[4]) {
      double acc;
      double startValue[2], stepSize[2], lowerLimit[2], upperLimit[2];
      short disable[2] = {0, 0};
      /*
      printf("Reoptimizing CCBEND due to changes: (optimized=%hd)\n", ccbend->optimized);
      printf("delta L: %le\n", ccbend->length-ccbend->referenceData[0]);
      printf("delta ANGLE: %le\n", ccbend->angle-ccbend->referenceData[1]);
      printf("delta K1: %le\n", ccbend->K1-ccbend->referenceData[2]);
      printf("delta K2: %le\n", ccbend->K2-ccbend->referenceData[3]);
      fflush(stdout);
      */
      if (iPart>=0)
        bombTracking("Programming error: oneStep mode is incompatible with optmization for CCBEND.");
      optimizationFlags = OPTIMIZE_X|OPTIMIZE_XP;
      if (ccbend->optimized) {
        startValue[0] = ccbend->fseOffset;
        startValue[1] = ccbend->dxOffset;
        if (ccbend->optimizeFseOnce) {
          disable[0] = 1;
          optimizationFlags &= ~OPTIMIZE_XP;
        }
        if (ccbend->optimizeDxOnce) {
          disable[1] = 1;
          optimizationFlags &= ~OPTIMIZE_X;
        }
      } else {
        startValue[0] = startValue[1] = 0;
        if (!ccbend->optimizeFse) {
          disable[0] = 1;
          optimizationFlags &= ~OPTIMIZE_XP;
        }
        if (!ccbend->optimizeDx) {
          disable[1] = 1;
          optimizationFlags &= ~OPTIMIZE_X;
        }
      }
      if (!disable[0] || !disable[1]) {
        ccbend->optimized = -1; /* flag to indicate calls to track_through_ccbend will be for FSE optimization */
        memcpy(&ccbendCopy, ccbend, sizeof(ccbendCopy));
        ccbendCopy.fse = ccbendCopy.dx = ccbendCopy.dy = ccbendCopy.dz = ccbendCopy.etilt = ccbendCopy.tilt = 
          ccbendCopy.isr = ccbendCopy.synch_rad = ccbendCopy.isr1Particle = ccbendCopy.KnDelta = 0;
        ccbendCopy.systematic_multipoles = ccbendCopy.edge_multipoles = ccbendCopy.edge1_multipoles = 
          ccbendCopy.edge2_multipoles = ccbendCopy.random_multipoles = NULL;
        PoCopy = Po;
        stepSize[0] = 1e-3; /* FSE */
        stepSize[1] = 1e-4; /* X */
        lowerLimit[0] = lowerLimit[1] = -1;
        upperLimit[0] = upperLimit[1] = 1;
        if (simplexMin(&acc, startValue, stepSize, lowerLimit, upperLimit, disable, 2, 
                       fabs(1e-14*ccbend->length), fabs(1e-16*ccbend->length),
                       ccbend_trajectory_error, NULL, 1500, 3, 12, 3.0, 1.0, 0)<0) {
          bombElegantVA("failed to find FSE and x offset to center trajectory for ccbend. accuracy acheived was %le.", acc);
        }
        ccbend->fseOffset = startValue[0];
        ccbend->dxOffset = startValue[1];
        ccbend->xAdjust = xFinal;
        ccbend->KnDelta = ccbendCopy.KnDelta;
        ccbend->referenceData[0] = ccbend->length;
        ccbend->referenceData[1] = ccbend->angle;
        ccbend->referenceData[2] = ccbend->K1;
        ccbend->referenceData[3] = ccbend->K2;
        ccbend->referenceData[4] = ccbend->yaw;
	if (ccbend->verbose) {
	  printf("CCBEND %s#%ld optimized: FSE=%le, dx=%le, accuracy=%le\n", 
		 eptr?eptr->name:"?", eptr?eptr->occurence:-1, ccbend->fseOffset, ccbend->dxOffset, acc);
          printf("xMin = %le, xMax = %le, xAve = %le, xFinal = %le, xpError = %le\n", 
                 xMin, xMax, xAve, xFinal, xpError);
	  fflush(stdout);
	}
        ccbend->optimized = 1;
      }
    }
  }

#ifdef DEBUG
  if (ccbend->optimized==1) {
    char filename[1000];
    logHamiltonian = 1;
    snprintf(filename, 1000, "ccbend-%s.sdds", eptr->name);
    fpHam = fopen(filename, "w");
    fprintf(fpHam, "SDDS1\n&column name=z type=double units=m &end\n");
    fprintf(fpHam, "&column name=x type=double units=m &end\n");
    fprintf(fpHam, "&column name=y type=double units=m &end\n");
    fprintf(fpHam, "&column name=qx type=double &end\n");
    fprintf(fpHam, "&column name=qy type=double &end\n");
    fprintf(fpHam, "&column name=dH type=double &end\n");
    fprintf(fpHam, "&data mode=ascii no_row_counts=1 &end\n");
  }
#endif

  rad_coef = isr_coef = 0;

  n_kicks = ccbend->n_kicks;
  arcLength = ccbend->length;
  if (ccbend->optimized!=0)
    fse = ccbend->fse + ccbend->fseOffset;
  else 
    fse = ccbend->fse;
  for (iTerm=0; iTerm<9; iTerm++)
    KnL[iTerm] = 0;
  length = rho0 = 0; /* prevent compiler warnings */
  if ((angle = ccbend->angle)!=0) {
    rho0 = arcLength/angle;
    length = 2*rho0*sin(angle/2);
    KnL[0] = (1+fse)/rho0*length;
  }
  else
    bombTracking("Can't have zero ANGLE for CCBEND.");
  KnL[1] = (1+fse)*ccbend->K1*length/(1-ccbend->KnDelta);
  KnL[2] = (1+fse)*ccbend->K2*length/(1-ccbend->KnDelta);
  KnL[3] = (1+fse)*ccbend->K3*length/(1-ccbend->KnDelta);
  KnL[4] = (1+fse)*ccbend->K4*length/(1-ccbend->KnDelta);
  KnL[5] = (1+fse)*ccbend->K5*length/(1-ccbend->KnDelta);
  KnL[6] = (1+fse)*ccbend->K6*length/(1-ccbend->KnDelta);
  KnL[7] = (1+fse)*ccbend->K7*length/(1-ccbend->KnDelta);
  KnL[8] = (1+fse)*ccbend->K8*length/(1-ccbend->KnDelta);
  if (angle<0) {
    angleSign = -1;
    for (iTerm=0; iTerm<9; iTerm+=2)
      KnL[iTerm] *= -1;
    angle = -angle;
    rho0 = -rho0;
    yaw = -ccbend->yaw*ccbend->edgeFlip;
    rotateBeamCoordinates(particle, n_part, PI);
  } else {
    angleSign = 1;
    yaw = ccbend->yaw*ccbend->edgeFlip;
  }
  if (ccbend->systematic_multipoles || ccbend->edge_multipoles || ccbend->random_multipoles ||
      ccbend->edge1_multipoles || ccbend->edge2_multipoles) {
    /* Note that with referenceOrder=0, referenceKnL will always be positive if angle is nonzero */
    if (ccbend->referenceOrder==0 && (referenceKnL=KnL[0])==0)
      bombElegant("REFERENCE_ORDER=0 but CCBEND ANGLE is zero", NULL);
    if (ccbend->referenceOrder==1 && (referenceKnL=KnL[1])==0)
      bombElegant("REFERENCE_ORDER=1 but CCBEND K1 is zero", NULL);
    if (ccbend->referenceOrder==2 && (referenceKnL=KnL[2])==0)
      bombElegant("REFERENCE_ORDER=2 but CCBEND K2 is zero", NULL);
    if (ccbend->referenceOrder<0 || ccbend->referenceOrder>2)
      bombElegant("REFERENCE_ORDER must be 0, 1, or 2 for CCBEND", NULL);
  }
  if (ccbend->edgeFlip==1) {
    gK[0] = 2*ccbend->fint1*ccbend->hgap;
    gK[1] = 2*ccbend->fint2*ccbend->hgap;
  } else {
    gK[1] = 2*ccbend->fint1*ccbend->hgap;
    gK[0] = 2*ccbend->fint2*ccbend->hgap;
  }

  integ_order = ccbend->integration_order;
  if (ccbend->synch_rad)
    rad_coef = sqr(particleCharge)*pow3(Po)/(6*PI*epsilon_o*sqr(c_mks)*particleMass); 
  isr_coef = particleRadius*sqrt(55.0/(24*sqrt(3))*pow5(Po)*137.0359895);
  if (!ccbend->isr || (ccbend->isr1Particle==0 && n_part==1))
    /* minus sign indicates we accumulate into sigmadelta^2 only, don't perturb particles */
    isr_coef *= -1;
  if (ccbend->length<1e-6 && (ccbend->isr || ccbend->synch_rad)) {
    rad_coef = isr_coef = 0;  /* avoid unphysical results */
  }
  if (!ccbend->multipolesInitialized) {
    /* read the data files for the error multipoles */
    readErrorMultipoleData(&(ccbend->systematicMultipoleData), ccbend->systematic_multipoles, 0);
    if ((ccbend->edge1_multipoles && !ccbend->edgeFlip) || (ccbend->edge2_multipoles && ccbend->edgeFlip))
        readErrorMultipoleData(&(ccbend->edge1MultipoleData), 
                               ccbend->edgeFlip?ccbend->edge2_multipoles:ccbend->edge1_multipoles, 0);
    else 
      readErrorMultipoleData(&(ccbend->edge1MultipoleData), ccbend->edge_multipoles, 0);
    if ((ccbend->edge2_multipoles && !ccbend->edgeFlip) || (ccbend->edge1_multipoles && ccbend->edgeFlip))
        readErrorMultipoleData(&(ccbend->edge2MultipoleData), 
                               ccbend->edgeFlip?ccbend->edge1_multipoles:ccbend->edge2_multipoles, 0);
    else 
      readErrorMultipoleData(&(ccbend->edge2MultipoleData), ccbend->edge_multipoles, 0);
    readErrorMultipoleData(&(ccbend->randomMultipoleData), ccbend->random_multipoles, 0);
    ccbend->multipolesInitialized = 1;
  }
  if (!ccbend->totalMultipolesComputed) {
    computeTotalErrorMultipoleFields(&(ccbend->totalMultipoleData),
                                     &(ccbend->systematicMultipoleData), ccbend->systematicMultipoleFactor, 
                                     &(ccbend->edge1MultipoleData), &(ccbend->edge2MultipoleData),
                                     &(ccbend->randomMultipoleData), ccbend->randomMultipoleFactor,
                                     NULL, 0.0,
                                     referenceKnL, 
                                     ccbend->referenceOrder, 0
                                     );
    if (angleSign<0) {
      long i;
      for (i=0; i<ccbend->systematicMultipoleData.orders; i++) {
        if (ccbend->totalMultipoleData.order[i]%2) {
          ccbend->totalMultipoleData.KnL[i] *= -1;
          ccbend->totalMultipoleData.JnL[i] *= -1;
        }
      }
      for (i=0; i<ccbend->edge1MultipoleData.orders; i++) {
        if (ccbend->totalMultipoleData.order[i]%2) {
          ccbend->edge1MultipoleData.KnL[i] *= -1;
          ccbend->edge1MultipoleData.JnL[i] *= -1;
        }
      }
      for (i=0; i<ccbend->edge2MultipoleData.orders; i++) {
        if (ccbend->totalMultipoleData.order[i]%2) {
          ccbend->edge2MultipoleData.KnL[i] *= -1;
          ccbend->edge2MultipoleData.JnL[i] *= -1;
        }
      }
    }
    ccbend->totalMultipolesComputed = 1;
  }
  multData = &(ccbend->totalMultipoleData);
  edge1MultData = &(ccbend->edge1MultipoleData);
  edge2MultData = &(ccbend->edge2MultipoleData);

  if (multData && !multData->initialized)
    multData = NULL;

  if (n_kicks<=0)
    bombTracking("n_kicks<=0 in track_ccbend()");
  if (integ_order!=2 && integ_order!=4) 
    bombTracking("multipole integration_order must be 2 or 4");

  if (!(coef = expansion_coefficients(0)))
    bombTracking("expansion_coefficients(0) returned NULL pointer (track_through_ccbend)");
  if (!(coef = expansion_coefficients(1)))
    bombTracking("expansion_coefficients(1) returned NULL pointer (track_through_ccbend)");
  if (!(coef = expansion_coefficients(2)))
    bombTracking("expansion_coefficients(2) returned NULL pointer (track_through_ccbend)");

  tilt = ccbend->tilt;
  dx = ccbend->dx;
  dy = ccbend->dy;
  dz = ccbend->dz;

  setupMultApertureData(&apertureData, maxamp, tilt, apFileData, z_start+length/2);

  if (iPart<=0) {
    xpError = particle[0][1];
    if (tilt)
      rotateBeamCoordinates(particle, n_part, tilt);
    switchRbendPlane(particle, n_part, angle/2-yaw, Po);
    if (dx || dy || dz)
      offsetBeamCoordinates(particle, n_part, dx, dy, dz);
    if (ccbend->optimized)
      offsetBeamCoordinates(particle, n_part, ccbend->dxOffset, 0, 0);
    verticalRbendFringe(particle, n_part, angle/2-yaw, rho0, KnL[1]/length, KnL[2]/length, gK[0], ccbend->edgeOrder);
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
                             integ_order, n_kicks, iPart, iFinalSlice, length, multData, edge1MultData, edge2MultData, 
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
  multipoleKicksDone += (i_top+1)*n_kicks;
  if (multData)
    multipoleKicksDone += (i_top+1)*n_kicks*multData->orders;
  if (sigmaDelta2)
    *sigmaDelta2 /= i_top+1;

  if ((iPart<0 || iPart==(ccbend->n_kicks-1)) && iFinalSlice<=0) {
    verticalRbendFringe(particle, i_top+1, angle/2+yaw, rho0, KnL[1]/length, KnL[2]/length, gK[1], ccbend->edgeOrder);
    if (ccbend->optimized)
      offsetBeamCoordinates(particle, i_top+1, ccbend->xAdjust, 0, 0);
    if (dx || dy || dz)
      offsetBeamCoordinates(particle, i_top+1, -dx, -dy, -dz);
    switchRbendPlane(particle, i_top+1, angle/2+yaw, Po);
    if (tilt)
      rotateBeamCoordinates(particle, i_top+1, -tilt);
    xpError = fabs(xpError+particle[0][1]);
  }

  if (angleSign<0)
    /* note that we use n_part here so lost particles get rotated back as well */
    rotateBeamCoordinates(particle, n_part, -PI);

  if (freeMultData && !multData->copy) {
    if (multData->order)
      free(multData->order);
    if (multData->KnL)
      free(multData->KnL);
    free(multData);
  }

  if (angleSign<0) {
    lastRho *= -1;
    lastX *= -1;
    lastXp *= -1;
  }

  if (fpHam) {
    fclose(fpHam);
    fpHam = NULL;
  }

  log_exit("track_through_ccbend");
  return(i_top+1);
}


/* beta is 2^(1/3) */
#define BETA 1.25992104989487316477

int integrate_kick_K012(double *coord, /* coordinates of the particle */
                        double dx, double dy,  /* misalignments, needed for aperture checks */
                        double Po, double rad_coef, double isr_coef, /* radiation effects */
                        double K0L, double K1L, double K2L, /* strength parameters for full magnet */
                        long integration_order, /* 2 or 4 */
                        long n_parts, /* N_KICKS */
                        long iPart,   /* If <0, integrate the full magnet. If >=0, integrate just a single part and return.
                                       * This is needed to allow propagation of the radiation matrix. */
                        long iFinalSlice, /* If >0, integrate to the indicated slice. Needed to allow extracting the
                                           * interior matrix from tracking data. */
                        double drift, /* length of the full element */
                        MULTIPOLE_DATA *multData, MULTIPOLE_DATA *edgeMultData, /* error multipoles */
                        MULT_APERTURE_DATA *apData,  /* aperture */
                        double *dzLoss,      /* if particle is loss, offset from start of element where this occurs */
                        double *sigmaDelta2  /* accumulate the energy spread increase for propagation of radiation matrix */
                        )
{
  double p, qx, qy, denom, beta0, beta1, dp, s;
  double x, y, xp, yp, sum_Fx, sum_Fy;
  double xSum;
  long i_kick, step, steps, iMult, nSum;
  double dsh;
  long maxOrder;
  double *xpow, *ypow;
  static double driftFrac[4] = {
    0.5/(2-BETA),  (1-BETA)/(2-BETA)/2,  (1-BETA)/(2-BETA)/2,  0.5/(2-BETA)
    } ;
  static double kickFrac[4] = {
    1./(2-BETA),  -BETA/(2-BETA),  1/(2-BETA),  0
    } ;
  if (integration_order==2) {
    driftFrac[0] = driftFrac[1] = 0.5;
    kickFrac[0] = 1;
    kickFrac[1] = 0;
    steps = 2;
  } else
    steps = 4;

  drift = drift/n_parts;
  K0L /= n_parts;
  K1L /= n_parts;
  K2L /= n_parts;

  x = coord[0];
  xp = coord[1];
  y = coord[2];
  yp = coord[3];
  s  = 0;
  dp = coord[5];
  p = Po*(1+dp);
  beta0 = p/sqrt(sqr(p)+1);

#if defined(ieee_math)
  if (isnan(x) || isnan(xp) || isnan(y) || isnan(yp)) {
    return 0;
  }
#endif
  if (fabs(x)>COORD_LIMIT || fabs(y)>COORD_LIMIT ||
      fabs(xp)>SLOPE_LIMIT || fabs(yp)>SLOPE_LIMIT) {
    return 0;
  }

  /* calculate initial canonical momenta */
  qx = (1+dp)*xp/(denom=sqrt(1+sqr(xp)+sqr(yp)));
  qy = (1+dp)*yp/denom;

  maxOrder = findMaximumOrder(1, 2, edgeMultData, multData, NULL);
  xpow = tmalloc(sizeof(*xpow)*(maxOrder+1));
  ypow = tmalloc(sizeof(*ypow)*(maxOrder+1));

  if (iPart<=0 && edgeMultData && edgeMultData->orders) {
    fillPowerArray(x, xpow, maxOrder);
    fillPowerArray(y, ypow, maxOrder);
    for (iMult=0; iMult<edgeMultData->orders; iMult++) {
      apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, xpow, ypow, 
                                      edgeMultData->order[iMult], 
                                      edgeMultData->KnL[iMult], 0);
      apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, xpow, ypow, 
                                      edgeMultData->order[iMult], 
                                      edgeMultData->JnL[iMult], 1);
    }
  }
  /* we must do this to avoid numerical precision issues that may subtly change the results
   * when edge multipoles are enabled to disabled
   */
  if ((denom=sqr(1+dp)-sqr(qx)-sqr(qy))<=0) {
    coord[0] = x;
    coord[2] = y;
    free(xpow); free(ypow);
    return 0;
  }
  denom = sqrt(denom);
  xp = qx/denom;
  yp = qy/denom;

  *dzLoss = 0;
  if (iFinalSlice<=0)
    iFinalSlice = n_parts;
  xMin = DBL_MAX;
  xMax = -DBL_MAX;
  xSum = x;
  nSum = 1;
  for (i_kick=0; i_kick<iFinalSlice; i_kick++) {
#ifdef DEBUG
    double H0;
    if (logHamiltonian && fpHam) {
      double H;
      H = -sqrt(ipow(1+dp, 2)-qx*qx-qy*qy) +
        K0L/drift*x +
        K1L/(2*drift)*(x*x - y*y) +
        K2L/(6*drift)*(ipow(x, 3) - 3*x*y*y);
      if (i_kick==0) 
        H0 = H;
      fprintf(fpHam, "%e %e %e %e %le %le\n", i_kick*drift, x, y, qx, qy, H-H0);
      if (i_kick==n_parts)
        fputs("\n", fpHam);
    }
#endif
    if (apData && !checkMultAperture(x+dx, y+dy, apData))  {
      coord[0] = x;
      coord[2] = y;
      free(xpow); free(ypow);
      return 0;
    }
    for (step=0; step<steps; step++) {
      if (drift) {
        dsh = drift*driftFrac[step];
        x += xp*dsh;
        y += yp*dsh;
        s += dsh*sqrt(1 + sqr(xp) + sqr(yp));
        *dzLoss += dsh;
      }

      if (!kickFrac[step])
        break;

      fillPowerArray(x, xpow, maxOrder);
      fillPowerArray(y, ypow, maxOrder);

      sum_Fx = sum_Fy = 0;

      if (K0L)
        apply_canonical_multipole_kicks(&qx, &qy, &sum_Fx, &sum_Fy, xpow, ypow, 
                                        0, K0L*kickFrac[step], 0);
      if (K1L)
        apply_canonical_multipole_kicks(&qx, &qy, &sum_Fx, &sum_Fy, xpow, ypow, 
                                        1, K1L*kickFrac[step], 0);
      if (K2L)
        apply_canonical_multipole_kicks(&qx, &qy, &sum_Fx, &sum_Fy, xpow, ypow, 
                                        2, K2L*kickFrac[step], 0);

      if (multData) {
        /* do kicks for spurious multipoles */
        for (iMult=0; iMult<multData->orders; iMult++) {
          if (multData->KnL && multData->KnL[iMult]) {
            apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, xpow, ypow, 
                                            multData->order[iMult], 
                                            multData->KnL[iMult]*kickFrac[step]/n_parts,
                                            0);
          }
          if (multData->JnL && multData->JnL[iMult]) {
            apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, xpow, ypow, 
                                            multData->order[iMult], 
                                            multData->JnL[iMult]*kickFrac[step]/n_parts,
                                            1);
          }
        }
      }
      if ((denom=sqr(1+dp)-sqr(qx)-sqr(qy))<=0) {
        coord[0] = x;
        coord[2] = y;
        free(xpow); free(ypow);
        return 0;
      }
      xp = qx/(denom=sqrt(denom));
      yp = qy/denom;
      /* these three quantities are needed for radiation integrals */
      lastRho = 1/(K0L/drift + x*(K1L/drift) + x*x*(K2L/drift)/2);
      lastX = x;
      lastXp = xp;
      if ((rad_coef || isr_coef) && drift) {
	double deltaFactor, F2, dsFactor, dsIsrFactor;
        qx /= (1+dp);
        qy /= (1+dp);
	deltaFactor = sqr(1+dp);
	F2 = (sqr(sum_Fy)+sqr(sum_Fx))/sqr(lastRho);
	dsFactor = sqrt(1+sqr(xp)+sqr(yp));
	dsIsrFactor = dsFactor*drift/3;   /* recall that kickFrac may be negative */
	dsFactor *= drift*kickFrac[step]; /* that's ok here, since we don't take sqrt */
	if (rad_coef)
	  dp -= rad_coef*deltaFactor*F2*dsFactor;
	if (isr_coef>0)
	  dp -= isr_coef*deltaFactor*pow(F2, 0.75)*sqrt(dsIsrFactor)*gauss_rn_lim(0.0, 1.0, srGaussianLimit, random_2);
        if (sigmaDelta2)
          *sigmaDelta2 += sqr(isr_coef*deltaFactor)*pow(F2, 1.5)*dsFactor;
        qx *= (1+dp);
        qy *= (1+dp);
      }
    }
    xSum += x;
    nSum ++;
    if (x>xMax) xMax = x;
    if (x<xMin) xMin = x;
    if (iPart>=0)
      break;
  }
  xFinal = x;
  xAve = xSum/nSum;
  xError = xMax + xMin;
  /*
  printf("x init, min, max, fin = %le, %le, %le, %le\n", x0, xMin, xMax, x);
  printf("xp init, fin = %le, %le\n", xp0, xp);
  */

  if (apData && !checkMultAperture(x+dx, y+dy, apData))  {
    coord[0] = x;
    coord[2] = y;
    free(xpow); free(ypow);
    return 0;
  }
  
  if ((iPart<0 || iPart==n_parts) && (iFinalSlice==n_parts-1) && edgeMultData && edgeMultData->orders) {
    fillPowerArray(x, xpow, maxOrder);
    fillPowerArray(y, ypow, maxOrder);
    for (iMult=0; iMult<edgeMultData->orders; iMult++) {
      apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, xpow, ypow, 
                                      edgeMultData->order[iMult], 
                                      edgeMultData->KnL[iMult], 0);
      apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, xpow, ypow, 
                                      edgeMultData->order[iMult], 
                                      edgeMultData->JnL[iMult], 1);
    }
  }
  if ((denom=sqr(1+dp)-sqr(qx)-sqr(qy))<=0) {
    coord[0] = x;
    coord[2] = y;
    free(xpow); free(ypow);
    return 0;
  }
  denom = sqrt(denom);
  xp = qx/denom;
  yp = qy/denom;

  free(xpow);
  free(ypow);

  coord[0] = x;
  coord[1] = xp;
  coord[2] = y;
  coord[3] = yp;
  if (rad_coef) {
    p = Po*(1+dp);
    beta1 = p/sqrt(sqr(p)+1);
    coord[4] = beta1*(coord[4]/beta0 + 2*s/(beta0+beta1));
  }
  else 
    coord[4] += s;
  coord[5] = dp;

#if defined(ieee_math)
  if (isnan(x) || isnan(xp) || isnan(y) || isnan(yp)) {
    return 0;
  }
#endif
  if (fabs(x)>COORD_LIMIT || fabs(y)>COORD_LIMIT ||
      fabs(xp)>SLOPE_LIMIT || fabs(yp)>SLOPE_LIMIT) {
    return 0;
  }
  return 1;
}

int integrate_kick_KnL(double *coord, /* coordinates of the particle */
                       double dx, double dy,  /* misalignments, needed for aperture checks */
                       double Po, double rad_coef, double isr_coef, /* radiation effects */
                       double *KnLFull, 
                       long nTerms,
                       long integration_order, /* 2 or 4 */
                       long n_parts, /* N_KICKS */
                       long iPart,   /* If <0, integrate the full magnet. If >=0, integrate just a single part and return.
                                      * This is needed to allow propagation of the radiation matrix. */
                       long iFinalSlice, /* If >0, integrate to the indicated slice. Needed to allow extracting the
                                          * interior matrix from tracking data. */
                       double drift, /* length of the full element */
                       /* error multipoles */
                       MULTIPOLE_DATA *multData,       /* body terms */
                       MULTIPOLE_DATA *edge1MultData,  /* entrance */
                       MULTIPOLE_DATA *edge2MultData,  /* exit */
                       MULT_APERTURE_DATA *apData,  /* aperture */
                       double *dzLoss,      /* if particle is loss, offset from start of element where this occurs */
                       double *sigmaDelta2  /* accumulate the energy spread increase for propagation of radiation matrix */
                       )
{
  double p, qx, qy, denom, beta0, beta1, dp, s;
  double x, y, xp, yp, sum_Fx, sum_Fy;
  double xSum;
  long i_kick, step, steps, iMult, nSum;
  double dsh;
  long maxOrder, iTerm;
  double *xpow, *ypow, *KnL;
  static double driftFrac[4] = {
    0.5/(2-BETA),  (1-BETA)/(2-BETA)/2,  (1-BETA)/(2-BETA)/2,  0.5/(2-BETA)
    } ;
  static double kickFrac[4] = {
    1./(2-BETA),  -BETA/(2-BETA),  1/(2-BETA),  0
    } ;
  if (integration_order==2) {
    driftFrac[0] = driftFrac[1] = 0.5;
    kickFrac[0] = 1;
    kickFrac[1] = 0;
    steps = 2;
  } else
    steps = 4;

  drift = drift/n_parts;

  KnL = tmalloc(sizeof(*KnL)*nTerms);
  for (iTerm=0; iTerm<nTerms; iTerm++)
    KnL[iTerm] = KnLFull[iTerm]/n_parts;

  x = coord[0];
  xp = coord[1];
  y = coord[2];
  yp = coord[3];
  s  = 0;
  dp = coord[5];
  p = Po*(1+dp);
  beta0 = p/sqrt(sqr(p)+1);

#if defined(ieee_math)
  if (isnan(x) || isnan(xp) || isnan(y) || isnan(yp)) {
    return 0;
  }
#endif
  if (fabs(x)>COORD_LIMIT || fabs(y)>COORD_LIMIT ||
      fabs(xp)>SLOPE_LIMIT || fabs(yp)>SLOPE_LIMIT) {
    return 0;
  }

  /* calculate initial canonical momenta */
  qx = (1+dp)*xp/(denom=sqrt(1+sqr(xp)+sqr(yp)));
  qy = (1+dp)*yp/denom;

  maxOrder = findMaximumOrder(1, nTerms, edge1MultData, edge2MultData, multData);
  xpow = tmalloc(sizeof(*xpow)*(maxOrder+1));
  ypow = tmalloc(sizeof(*ypow)*(maxOrder+1));

  if (iPart<=0 && edge1MultData && edge1MultData->orders) {
    fillPowerArray(x, xpow, maxOrder);
    fillPowerArray(y, ypow, maxOrder);
    for (iMult=0; iMult<edge1MultData->orders; iMult++) {
      apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, xpow, ypow, 
                                      edge1MultData->order[iMult], 
                                      edge1MultData->KnL[iMult], 0);
      apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, xpow, ypow, 
                                      edge1MultData->order[iMult], 
                                      edge1MultData->JnL[iMult], 1);
    }
  }
  /* we must do this to avoid numerical precision issues that may subtly change the results
   * when edge multipoles are enabled to disabled
   */
  if ((denom=sqr(1+dp)-sqr(qx)-sqr(qy))<=0) {
    coord[0] = x;
    coord[2] = y;
    free(xpow); free(ypow);
    free(KnL);
    return 0;
  }
  denom = sqrt(denom);
  xp = qx/denom;
  yp = qy/denom;

  *dzLoss = 0;
  if (iFinalSlice<=0)
    iFinalSlice = n_parts;
  xMin = DBL_MAX;
  xMax = -DBL_MAX;
  xSum = x;
  nSum = 1;
  for (i_kick=0; i_kick<iFinalSlice; i_kick++) {
#ifdef DEBUG
    double H0;
    if (logHamiltonian && fpHam) {
      double H;
      H = -sqrt(ipow(1+dp, 2)-qx*qx-qy*qy) +
        KnL[0]/drift*x +
        KnL[1]/(2*drift)*(x*x - y*y) +
        KnL[2]/(6*drift)*(ipow(x, 3) - 3*x*y*y);
      if (i_kick==0) 
        H0 = H;
      fprintf(fpHam, "%e %e %e %e %le %le\n", i_kick*drift, x, y, qx, qy, H-H0);
      if (i_kick==n_parts)
        fputs("\n", fpHam);
    }
#endif
    if (apData && !checkMultAperture(x+dx, y+dy, apData))  {
      coord[0] = x;
      coord[2] = y;
      free(xpow); free(ypow);
      free(KnL);
      return 0;
    }
    for (step=0; step<steps; step++) {
      if (drift) {
        dsh = drift*driftFrac[step];
        x += xp*dsh;
        y += yp*dsh;
        s += dsh*sqrt(1 + sqr(xp) + sqr(yp));
        *dzLoss += dsh;
      }

      if (!kickFrac[step])
        break;

      fillPowerArray(x, xpow, maxOrder);
      fillPowerArray(y, ypow, maxOrder);

      sum_Fx = sum_Fy = 0;

      for (iTerm=0; iTerm<nTerms; iTerm++)
        if (KnL[iTerm])
          apply_canonical_multipole_kicks(&qx, &qy, &sum_Fx, &sum_Fy, xpow, ypow, 
                                          iTerm, KnL[iTerm]*kickFrac[step], 0);

      if (multData) {
        /* do kicks for spurious multipoles */
        for (iMult=0; iMult<multData->orders; iMult++) {
          if (multData->KnL && multData->KnL[iMult]) {
            apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, xpow, ypow, 
                                            multData->order[iMult], 
                                            multData->KnL[iMult]*kickFrac[step]/n_parts,
                                            0);
          }
          if (multData->JnL && multData->JnL[iMult]) {
            apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, xpow, ypow, 
                                            multData->order[iMult], 
                                            multData->JnL[iMult]*kickFrac[step]/n_parts,
                                            1);
          }
        }
      }
      if ((denom=sqr(1+dp)-sqr(qx)-sqr(qy))<=0) {
        coord[0] = x;
        coord[2] = y;
        free(xpow); free(ypow);
        free(KnL);
        return 0;
      }
      xp = qx/(denom=sqrt(denom));
      yp = qy/denom;
      /* these three quantities are needed for radiation integrals */
      lastRho = 1/(KnL[0]/drift + x*(KnL[1]/drift) + x*x*(KnL[2]/drift)/2);
      lastX = x;
      lastXp = xp;
      if ((rad_coef || isr_coef) && drift) {
	double deltaFactor, F2, dsFactor, dsIsrFactor;
        qx /= (1+dp);
        qy /= (1+dp);
	deltaFactor = sqr(1+dp);
	F2 = (sqr(sum_Fy)+sqr(sum_Fx))/sqr(lastRho);
	dsFactor = sqrt(1+sqr(xp)+sqr(yp));
	dsIsrFactor = dsFactor*drift/3;   /* recall that kickFrac may be negative */
	dsFactor *= drift*kickFrac[step]; /* that's ok here, since we don't take sqrt */
	if (rad_coef)
	  dp -= rad_coef*deltaFactor*F2*dsFactor;
	if (isr_coef>0)
	  dp -= isr_coef*deltaFactor*pow(F2, 0.75)*sqrt(dsIsrFactor)*gauss_rn_lim(0.0, 1.0, srGaussianLimit, random_2);
        if (sigmaDelta2)
          *sigmaDelta2 += sqr(isr_coef*deltaFactor)*pow(F2, 1.5)*dsFactor;
        qx *= (1+dp);
        qy *= (1+dp);
      }
    }
    xSum += x;
    nSum ++;
    if (x>xMax) xMax = x;
    if (x<xMin) xMin = x;
    if (iPart>=0)
      break;
  }
  xFinal = x;
  xAve = xSum/nSum;
  xError = xMax + xMin;
  /*
  printf("x init, min, max, fin = %le, %le, %le, %le\n", x0, xMin, xMax, x);
  printf("xp init, fin = %le, %le\n", xp0, xp);
  */

  if (apData && !checkMultAperture(x+dx, y+dy, apData))  {
    coord[0] = x;
    coord[2] = y;
    free(xpow); free(ypow);
    free(KnL);
    return 0;
  }
  
  if ((iPart<0 || iPart==n_parts) && (iFinalSlice==n_parts-1) && edge2MultData && edge2MultData->orders) {
    fillPowerArray(x, xpow, maxOrder);
    fillPowerArray(y, ypow, maxOrder);
    for (iMult=0; iMult<edge2MultData->orders; iMult++) {
      apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, xpow, ypow, 
                                      edge2MultData->order[iMult], 
                                      edge2MultData->KnL[iMult], 0);
      apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, xpow, ypow, 
                                      edge2MultData->order[iMult], 
                                      edge2MultData->JnL[iMult], 1);
    }
  }
  if ((denom=sqr(1+dp)-sqr(qx)-sqr(qy))<=0) {
    coord[0] = x;
    coord[2] = y;
    free(xpow); free(ypow);
    free(KnL);
    return 0;
  }
  denom = sqrt(denom);
  xp = qx/denom;
  yp = qy/denom;

  free(xpow);
  free(ypow);
  free(KnL);

  coord[0] = x;
  coord[1] = xp;
  coord[2] = y;
  coord[3] = yp;
  if (rad_coef) {
    p = Po*(1+dp);
    beta1 = p/sqrt(sqr(p)+1);
    coord[4] = beta1*(coord[4]/beta0 + 2*s/(beta0+beta1));
  }
  else 
    coord[4] += s;
  coord[5] = dp;

#if defined(ieee_math)
  if (isnan(x) || isnan(xp) || isnan(y) || isnan(yp)) {
    return 0;
  }
#endif
  if (fabs(x)>COORD_LIMIT || fabs(y)>COORD_LIMIT ||
      fabs(xp)>SLOPE_LIMIT || fabs(yp)>SLOPE_LIMIT) {
    return 0;
  }
  return 1;
}


void switchRbendPlane(double **particle, long n_part, double alpha, double po)
/* transforms the reference plane to one that is at an angle alpha relative to the
 * initial plane. use alpha>0 at the entrance to an rbend, alpha<0 at the exit.
 * alpha = theta/2, where theta is the total bend angle.
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
  }
}

void verticalRbendFringe(double **particle, long n_part, double alpha, double rho0, double K1, double K2, double gK, long order)
{
  long i;
  double c, d, e;
  if (order<1)
    return;
  c = d = e = 0;
  if (gK!=0) 
    alpha -= gK/fabs(rho0)/cos(alpha)*(1+sqr(sin(alpha)));
  c = sin(alpha)/rho0;
  if (order>1)
    d = sin(alpha)*K1;
  if (order>2)
    e = sin(alpha)*K2/2;
  for (i=0; i<n_part; i++) {
    double x, y;
    x = particle[i][0];
    y = particle[i][2];
    particle[i][3] -= y*(c + d*x + e*(x*x-y*y))/(1+particle[i][5]);
  }
}

double ccbend_trajectory_error(double *value, long *invalid)
{
  static double **particle = NULL;
  double result;

  *invalid = 0;
  if ((ccbendCopy.fseOffset = value[0])>=1 || value[0]<=-1) {
    *invalid = 1;
    return DBL_MAX;
  }
  if (value[1]>=1 || value[1]<=-1) {
    *invalid = 1;
    return DBL_MAX;
  }
  if (!particle) 
    particle = (double**)czarray_2d(sizeof(**particle), 1, COORDINATES_PER_PARTICLE);
  memset(particle[0], 0, COORDINATES_PER_PARTICLE*sizeof(**particle));
  ccbendCopy.dxOffset = value[1];
  if (ccbendCopy.compensateKn)
    ccbendCopy.KnDelta = -ccbendCopy.fseOffset;
  /* printf("** fse = %le, dx = %le, x[0] = %le\n", value[0], value[1], particle[0][0]); fflush(stdout);  */
  if (!track_through_ccbend(particle, 1, NULL, &ccbendCopy, PoCopy, NULL, 0.0, NULL, NULL, NULL, NULL, -1, -1)) {
    *invalid = 1;
    return DBL_MAX;
  }
  result = 0;
  if (optimizationFlags&OPTIMIZE_X)
    result += fabs(xError);
  if (optimizationFlags&OPTIMIZE_XP)
    result += fabs(xpError/ccbendCopy.angle);
  /*
  printf("particle[0][0] = %le, particle[0][1] = %le, xError = %le, xpError = %le, result = %le\n", 
         particle[0][0], particle[0][1], xError, xpError, result); fflush(stdout);
  */
  return result;
}

VMATRIX *determinePartialCcbendLinearMatrix(CCBEND *ccbend, double *startingCoord, double pCentral, long iFinalSlice)
/* This routine is used for getting the linear transport matrix from the start of the CCBEND to some interior slice.
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
   		 
  coord = (double**)czarray_2d(sizeof(**coord), 1+6*4, COORDINATES_PER_PARTICLE);

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
  ltmp1 = ccbend->isr;
  ltmp2 = ccbend->synch_rad;

  ccbend->isr = ccbend->synch_rad = 0;
  track_through_ccbend(coord, n_track, NULL, ccbend, pCentral, NULL, 0.0, NULL, NULL, NULL, NULL, -1, iFinalSlice);

  ccbend->isr = ltmp1;
  ccbend->synch_rad = ltmp2;
  
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

  free_czarray_2d((void**)coord, 1+4*6, COORDINATES_PER_PARTICLE);

#if USE_MPI
  notSinglePart = notSinglePart_saved;
#endif
  return M;
}

#undef DEBUG

void addCcbendRadiationIntegrals(CCBEND *ccbend, double *startingCoord, double pCentral,
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
    fpcr = fopen("ccbend-RI.sdds", "w");
    fprintf(fpcr, "SDDS1\n&column name=s type=double units=m &end\n");
    fprintf(fpcr, "&column name=betax type=double units=m &end\n");
    fprintf(fpcr, "&column name=etax type=double units=m &end\n");
    fprintf(fpcr, "&column name=etaxp type=double &end\n");
    fprintf(fpcr, "&column name=alphax type=double &end\n");
    fprintf(fpcr, "&column name=ElementName type=string &end\n");
    fprintf(fpcr, "&data mode=ascii no_row_counts=1 &end\n");
  }
  s0 = elem->end_pos - ccbend->length;
  fprintf(fpcr, "%le %le %le %le %le %s\n", s0, beta0, eta0, etap0, alpha0, elem->name);
#endif
  
  if (ccbend->tilt)
    bombElegant("Can't add radiation integrals for tilted CCBEND\n", NULL);

  gamma0 = (1+alpha0*alpha0)/beta0;

  beta1 = beta0;
  eta1 = eta0;
  etap1 = etap0;
  alpha1 = alpha0;
  H1 = (eta1*eta1 + sqr(beta1*etap1 + alpha1*eta1))/beta1;
  ds = ccbend->length/ccbend->n_kicks; /* not really right... */
  for (iSlice=1; iSlice<=ccbend->n_kicks; iSlice++) {
    /* Determine matrix from start of element to exit of slice iSlice */
    M = determinePartialCcbendLinearMatrix(ccbend, startingCoord, pCentral, iSlice);
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
    fprintf(fpcr, "%le %le %le %le %le %s\n", s0, beta2, eta2, etap2, alpha2, elem->name);
#endif

    /* Compute contributions to radiation integrals in this slice.
     * lastRho is saved by the routine track_through_ccbend(). Since determinePartialCcbendLinearMatrix()
     * puts the reference particle first, it is appropriate to the central trajectory. 
     */
    *I1 += ds*(eta1+eta2)/2/lastRho;
    *I2 += ds/sqr(lastRho);
    *I3 += ds/ipow(fabs(lastRho), 3);
    /* Compute effective K1 including the sextupole effect plus rotation.
     * lastX and lastXp are saved by track_through_ccbend().
     */
    K1 = (ccbend->K1+(lastX-ccbend->dxOffset)*ccbend->K2)*cos(atan(lastXp));
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
