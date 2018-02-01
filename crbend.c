/*************************************************************************\
* Copyright (c) 2018 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2018 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* routine: crbend()
 * purpose: tracking through canonical rectangular bend
 * 
 * Michael Borland, 2018
 */
#include "mdb.h"
#include "track.h"
#include "multipole.h"

static CRBEND *crbendCopy;
static double PoCopy, xMax, xError, xpError;

void switchRbendPlane(double **particle, long n_part, double alpha, double Po);
void verticalRbendFringe(double **particle, long n_part, double alpha, double rho0);
int integrate_kick_K012(double *coord, double dx, double dy, 
                        double Po, double rad_coef, double isr_coef,
                        double K0L, double K1L, double K2L,
                        long integration_order, long n_parts, double drift,
                        MULTIPOLE_DATA *multData, MULTIPOLE_DATA *edgeMultData, 
                        MULT_APERTURE_DATA *apData, double *dzLoss, double *sigmaDelta2);
double crbend_trajectory_error(double *value, long *invalid);

long track_through_crbend(
                          double **particle,   /* initial/final phase-space coordinates */
                          long n_part,         /* number of particles */
                          CRBEND *crbend,
                          double Po,
                          double **accepted,
                          double z_start,
                          double *sigmaDelta2,
                          char *rootname,
                          MAXAMP *maxamp,
                          APERTURE_DATA *apFileData
                          )
{
  double K0L, K1L, K2L;
  double dx, dy, dz; /* offsets of the multipole center */
  long n_kicks, integ_order;
  long i_part, i_top;
  double *coef;
  double fse, tilt, rad_coef, isr_coef, dzLoss=0, multSign = 1;
  double rho0, arcLength, length, angle;
  MULTIPOLE_DATA *multData = NULL, *edgeMultData = NULL;
  long freeMultData=0;
  MULT_APERTURE_DATA apertureData;
  double referenceKnL;

  if (!particle)
    bombTracking("particle array is null (track_through_crbend)");
  
  if (crbend->optimized!=-1 && crbend->angle!=0) {
    if (crbend->optimized==0 ||
        crbend->length!=crbend->referenceData[0] ||
        crbend->angle!=crbend->referenceData[1] ||
        crbend->K1!=crbend->referenceData[2] ||
        crbend->K2!=crbend->referenceData[3]) {
      double acc;
      double startValue[2], stepSize[2], lowerLimit[2], upperLimit[2];
      crbend->optimized = -1; /* flag to indicate calls to track_through_crbend will be for FSE optimization */
      crbendCopy = crbend;
      PoCopy = Po;
      startValue[0] = startValue[1] = 0;
      stepSize[0] = 1e-3; /* FSE */
      stepSize[1] = 1e-4; /* X */
      lowerLimit[0] = lowerLimit[1] = -1;
      upperLimit[0] = upperLimit[1] = 1;
      if (simplexMin(&acc, startValue, stepSize, lowerLimit, upperLimit, NULL, 2, 
                     fabs(1e-14*crbend->length/crbend->angle), fabs(1e-16*crbend->length/crbend->angle),
                     crbend_trajectory_error, NULL, 1500, 3, 12, 3.0, 1.0, 0)<0) {
        bombElegantVA("failed to find FSE and x offset to center trajectory for crbend. accuracy acheived was %le.", acc);
      }
      crbend->fseOffset = startValue[0];
      crbend->dxOffset = startValue[1];
      crbend->referenceData[0] = crbend->length;
      crbend->referenceData[1] = crbend->angle;
      crbend->referenceData[2] = crbend->K1;
      crbend->referenceData[3] = crbend->K2;
      printf("CRBEND optimized FSE=%le, x=%le, accuracy=%le\n", crbend->fseOffset, crbend->dxOffset, acc);
      crbend->optimized = 1;
    }
  }

  rad_coef = isr_coef = 0;

  n_kicks = crbend->n_kicks;
  arcLength = crbend->length;
  if (crbend->optimized!=0)
    fse = crbend->fse + crbend->fseOffset;
  else 
    fse = crbend->fse;
  if ((angle = crbend->angle)!=0) {
    rho0 = arcLength/angle;
    length = 2*rho0*sin(angle/2);
    K0L = (1+fse)/rho0*length;
  }
  else {
    rho0 = DBL_MAX;
    length = arcLength;
    K0L = 0;
  }
  K1L = (1+fse)*crbend->K1*length;
  K2L = (1+fse)*crbend->K2*length;
  if (crbend->systematic_multipoles || crbend->edge_multipoles || crbend->random_multipoles) {
    if (crbend->referenceOrder==0 && (referenceKnL=K0L)==0)
        bombElegant("REFERENCE_ORDER=0 but CRBEND ANGLE is zero", NULL);
    if (crbend->referenceOrder==1 && (referenceKnL=K1L)==0)
      bombElegant("REFERENCE_ORDER=1 but CRBEND K1 is zero", NULL);
    if (crbend->referenceOrder==2 && (referenceKnL=K2L)==0)
      bombElegant("REFERENCE_ORDER=2 but CRBEND K2 is zero", NULL);
    if (crbend->referenceOrder<0 || crbend->referenceOrder>2)
      bombElegant("REFERENCE_ORDER must be 0, 1, or 2 for CRBEND", NULL);
  }
  if (angle<0) {
    K0L *= -1;
    K1L *= -1;
    K2L *= -1;
    multSign = -1;
    rotateBeamCoordinates(particle, n_part, PI);
  }

  integ_order = crbend->integration_order;
  if (crbend->synch_rad)
    rad_coef = sqr(particleCharge)*pow3(Po)/(6*PI*epsilon_o*sqr(c_mks)*particleMass); 
  isr_coef = particleRadius*sqrt(55.0/(24*sqrt(3))*pow5(Po)*137.0359895);
  if (!crbend->isr || (crbend->isr1Particle==0 && n_part==1))
    /* minus sign indicates we accumulate into sigmadelta^2 only, don't perturb particles */
    isr_coef *= -1;
  if (crbend->length<1e-6 && (crbend->isr || crbend->synch_rad)) {
    rad_coef = isr_coef = 0;  /* avoid unphysical results */
  }
  if (!crbend->multipolesInitialized) {
    /* read the data files for the error multipoles */
    readErrorMultipoleData(&(crbend->systematicMultipoleData), crbend->systematic_multipoles, 0);
    readErrorMultipoleData(&(crbend->edgeMultipoleData), crbend->edge_multipoles, 0);
    readErrorMultipoleData(&(crbend->randomMultipoleData), crbend->random_multipoles, 0);
    crbend->multipolesInitialized = 1;
  }
  if (!crbend->totalMultipolesComputed) {
    computeTotalErrorMultipoleFields(&(crbend->totalMultipoleData),
                                     &(crbend->systematicMultipoleData), multSign*crbend->systematicMultipoleFactor, 
                                     &(crbend->edgeMultipoleData),
                                     &(crbend->randomMultipoleData), multSign*crbend->randomMultipoleFactor,
                                     NULL, 0.0,
                                     referenceKnL, crbend->referenceOrder);
    crbend->totalMultipolesComputed = 1;
  }
  multData = &(crbend->totalMultipoleData);
  edgeMultData = &(crbend->edgeMultipoleData);

  if (multData && !multData->initialized)
    multData = NULL;

  if (crbend->optimized==-1) {
    /* ignore multipole errors during fse optimization */
    multData = NULL;
    edgeMultData = NULL;
    isr_coef = 0;
    rad_coef = 0;
  }

  if (n_kicks<=0)
    bombTracking("n_kicks<=0 in track_crbend()");
  if (integ_order!=2 && integ_order!=4) 
    bombTracking("multipole integration_order must be 2 or 4");

  if (!(coef = expansion_coefficients(0)))
    bombTracking("expansion_coefficients(0) returned NULL pointer (track_through_crbend)");
  if (!(coef = expansion_coefficients(1)))
    bombTracking("expansion_coefficients(1) returned NULL pointer (track_through_crbend)");
  if (!(coef = expansion_coefficients(2)))
    bombTracking("expansion_coefficients(2) returned NULL pointer (track_through_crbend)");

  tilt = crbend->tilt;
  dx = crbend->dx + (crbend->optimized ? crbend->dxOffset : 0);
  dy = crbend->dy;
  dz = crbend->dz;

  setupMultApertureData(&apertureData, maxamp, tilt, apFileData, z_start+length/2);

  if (angle!=0) {
    switchRbendPlane(particle, n_part, fabs(angle/2), Po);
    verticalRbendFringe(particle, n_part, fabs(angle/2), rho0);
  }

  if (dx || dy || dz)
    offsetBeamCoordinates(particle, n_part, dx, dy, dz);
  if (tilt)
    rotateBeamCoordinates(particle, n_part, tilt);

  if (sigmaDelta2)
    *sigmaDelta2 = 0;
  i_top = n_part-1;
  for (i_part=0; i_part<=i_top; i_part++) {
    if (!integrate_kick_K012(particle[i_part], dx, dy, Po, rad_coef, isr_coef, K0L, K1L, K2L, 
                             integ_order, n_kicks, length, multData, edgeMultData, 
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

  if (tilt)
    rotateBeamCoordinates(particle, n_part, -tilt);
  if (dx || dy || dz)
    offsetBeamCoordinates(particle, n_part, -dx, -dy, -dz);
  if (crbend->optimized!=-1) {
    if (angle!=0) {
      verticalRbendFringe(particle, n_part, fabs(angle/2), rho0);
      switchRbendPlane(particle, n_part, fabs(angle/2), Po);
    }
  }

  if (angle<0)
    /* note that we use n_part here so lost particles get rotated back as well */
    rotateBeamCoordinates(particle, n_part, -PI);

  if (freeMultData && !multData->copy) {
    if (multData->order)
      free(multData->order);
    if (multData->KnL)
      free(multData->KnL);
    free(multData);
  }

  log_exit("track_through_crbend");
  return(i_top+1);
}


/* beta is 2^(1/3) */
#define BETA 1.25992104989487316477

int integrate_kick_K012(double *coord, double dx, double dy, 
                        double Po, double rad_coef, double isr_coef,
                        double K0L, double K1L, double K2L,
                        long integration_order, long n_parts, double drift,
                        MULTIPOLE_DATA *multData, MULTIPOLE_DATA *edgeMultData, 
                        MULT_APERTURE_DATA *apData, double *dzLoss, double *sigmaDelta2)
{
  double p, qx, qy, denom, beta0, beta1, dp, s;
  double x, y, xp, yp, sum_Fx, sum_Fy;
  double xMin, x0, xp0;
  long i_kick, step, steps, iMult;
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

  if (edgeMultData && edgeMultData->orders) {
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
    return 0;
  }
  denom = sqrt(denom);
  xp = qx/denom;
  yp = qy/denom;

  *dzLoss = 0;
  xMax = xMin = x0 = x;
  xp0 = xp;
  for (i_kick=0; i_kick<n_parts; i_kick++) {
    if (apData && !checkMultAperture(x+dx, y+dy, apData))  {
      coord[0] = x;
      coord[2] = y;
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
      if (x>xMax)
        xMax = x;
      if (x<xMin)
        xMin = x;

      if (!kickFrac[step])
        break;

      fillPowerArray(x, xpow, maxOrder);
      fillPowerArray(y, ypow, maxOrder);

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
        return 0;
      }
      xp = qx/(denom=sqrt(denom));
      yp = qy/denom;
      if ((rad_coef || isr_coef) && drift) {
	double deltaFactor, F2, dsFactor, dsIsrFactor;
        qx /= (1+dp);
        qy /= (1+dp);
	deltaFactor = sqr(1+dp);
	F2 = (sqr(sum_Fy)+sqr(sum_Fx))*sqr(K0L/drift);
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
  }
  xError = fabs(x - x0) + fabs(x + xMax);
  xpError =  fabs(xp0+xp);
  /*
  printf("x init, min, max, fin = %le, %le, %le, %le\n", x0, xMin, xMax, x);
  printf("xp init, fin = %le, %le\n", xp0, xp);
  */

  if (apData && !checkMultAperture(x+dx, y+dy, apData))  {
    coord[0] = x;
    coord[2] = y;
    return 0;
  }
  
  if (edgeMultData && edgeMultData->orders) {
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

void verticalRbendFringe(double **particle, long n_part, double alpha, double rho0)
{
  long i;
  double c;
  c = sin(alpha)/rho0;
  for (i=0; i<n_part; i++) 
    particle[i][3] -= particle[i][2]*c/(1+particle[i][5]);
}

double crbend_trajectory_error(double *value, long *invalid)
{
  static double **particle = NULL;
  double result;

  *invalid = 0;
  if ((crbendCopy->fseOffset = value[0])>=1 || value[0]<=-1) {
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
  crbendCopy->dxOffset = value[1];
  /* printf("** fse = %le, dx = %le, x[0] = %le\n", value[0], value[1], particle[0][0]); */
  if (!track_through_crbend(particle, 1, crbendCopy, PoCopy, NULL, 0.0, NULL, NULL, NULL, NULL)) {
    *invalid = 1;
    return DBL_MAX;
  }
  result = xError + xpError/fabs(crbendCopy->angle);
  /* printf("particle[0][0] = %le, particle[0][1] = %le, result = %le\n", particle[0][0], particle[0][1], result); */
  return result;
}
