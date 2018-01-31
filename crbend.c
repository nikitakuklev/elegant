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

void switchRbendPlane(double **particle, long n_part, double alpha, double Po);
void verticalRbendFringe(double **particle, long n_part, double alpha, double rho0);
int integrate_kick_K012(double *coord, double dx, double dy, 
                        double Po, double rad_coef, double isr_coef,
                        double K0L, double K1L, double K2L,
                        long integration_order, long n_parts, double drift,
                        MULTIPOLE_DATA *multData, MULTIPOLE_DATA *edgeMultData, 
                        MULT_APERTURE_DATA *apData, double *dzLoss, double *sigmaDelta2);

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
  double tilt, rad_coef, isr_coef, dzLoss=0, multSign = 1;
  double rho0, arcLength, length, angle;
  MULTIPOLE_DATA *multData = NULL, *edgeMultData = NULL;
  long freeMultData=0;
  MULT_APERTURE_DATA apertureData;
  
  if (!particle)
    bombTracking("particle array is null (crbend_tracking)");

  rad_coef = isr_coef = 0;

  n_kicks = crbend->n_kicks;
  arcLength = crbend->length;
  if ((angle = crbend->angle)!=0) {
    rho0 = arcLength/angle;
    length = 2*rho0*sin(angle/2);
    K0L = (1+crbend->fse)/rho0*length;
  }
  else {
    rho0 = DBL_MAX;
    length = arcLength;
    K0L = 0;
  }
  K1L = (1+crbend->fse)*crbend->K1*length;
  K2L = (1+crbend->fse)*crbend->K2*length;
  if (angle<0) {
    K0L *= -1;
    K1L *= -1;
    K2L *= -1;
    multSign = -1;
    rotateBeamCoordinates(particle, n_part, PI);
  }
  tilt = crbend->tilt;
  dx = crbend->dx;
  dy = crbend->dy;
  dz = crbend->dz;

  integ_order = crbend->integration_order;
  if (crbend->synch_rad)
    rad_coef = sqr(particleCharge)*pow3(Po)/(6*PI*epsilon_o*sqr(c_mks)*particleMass); 
  isr_coef = particleRadius*sqrt(55.0/(24*sqrt(3))*pow5(Po)*137.0359895);
  if (!crbend->isr || (crbend->isr1Particle==0 && n_part==1))
    /* Minus sign indicates we accumulate into sigmaDelta^2 only, don't perturb particles */
    isr_coef *= -1;
  if (crbend->length<1e-6 && (crbend->isr || crbend->synch_rad)) {
    rad_coef = isr_coef = 0;  /* avoid unphysical results */
  }
  if (!crbend->multipolesInitialized) {
    /* read the data files for the error multipoles */
    readErrorMultipoleData(&(crbend->systematicMultipoleData),
                           crbend->systematic_multipoles, 0);
    readErrorMultipoleData(&(crbend->edgeMultipoleData),
                           crbend->edge_multipoles, 0);
    readErrorMultipoleData(&(crbend->randomMultipoleData),
                           crbend->random_multipoles, 0);
    crbend->multipolesInitialized = 1;
  }
  if (!crbend->totalMultipolesComputed) {
    computeTotalErrorMultipoleFields(&(crbend->totalMultipoleData),
                                     &(crbend->systematicMultipoleData), multSign*crbend->systematicMultipoleFactor, 
                                     &(crbend->edgeMultipoleData),
                                     &(crbend->randomMultipoleData), multSign*crbend->randomMultipoleFactor,
                                     NULL, 0.0,
                                     K0L, 1);
    crbend->totalMultipolesComputed = 1;
  }
  multData = &(crbend->totalMultipoleData);
  edgeMultData = &(crbend->edgeMultipoleData);

  if (multData && !multData->initialized)
    multData = NULL;
  
  if (n_kicks<=0)
    bombTracking("n_kicks<=0 in track_crbend()");
  if (integ_order!=2 && integ_order!=4) 
    bombTracking("multipole integration_order must be 2 or 4");

  if (!(coef = expansion_coefficients(0)))
    bombTracking("expansion_coefficients(0) returned null pointer (crbend_tracking)");
  if (!(coef = expansion_coefficients(1)))
    bombTracking("expansion_coefficients(1) returned null pointer (crbend_tracking)");
  if (!(coef = expansion_coefficients(2)))
    bombTracking("expansion_coefficients(2) returned null pointer (crbend_tracking)");

  setupMultApertureData(&apertureData, maxamp, tilt, apFileData, z_start+length/2);

  printf("Need to implement misalignments in crbend_tracking\n");
  /*
  if (dx || dy || dz)
    offsetBeamCoordinates(particle, n_part, dx, dy, dz);
  if (tilt)
    rotateBeamCoordinates(particle, n_part, tilt);
  */

  if (angle!=0) {
    switchRbendPlane(particle, n_part, fabs(angle/2), Po);
    verticalRbendFringe(particle, n_part, fabs(angle/2), rho0);
  }

  if (sigmaDelta2)
    *sigmaDelta2 = 0;
  i_top = n_part-1;
  for (i_part=0; i_part<=i_top; i_part++) {
    if (!integrate_kick_K012(particle[i_part], dx, dy, Po, rad_coef, isr_coef, K0L, K1L, K2L, integ_order, n_kicks, length, 
                            multData, edgeMultData, &apertureData, &dzLoss, sigmaDelta2)) {
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

  if (angle!=0) {
    verticalRbendFringe(particle, n_part, fabs(angle/2), rho0);
    switchRbendPlane(particle, n_part, fabs(angle/2), Po);
  }

  printf("Need to implement misalignments in crbend_tracking\n");
  /*
  if (tilt)
    rotateBeamCoordinates(particle, n_part, -tilt);
  if (dx || dy || dz)
    offsetBeamCoordinates(particle, n_part, -dx, -dy, -dz);
  */

  if (angle<0)
    /* Note that we use n_part here so lost particles get rotated back as well */
    rotateBeamCoordinates(particle, n_part, -PI);

  if (freeMultData && !multData->copy) {
    if (multData->order)
      free(multData->order);
    if (multData->KnL)
      free(multData->KnL);
    free(multData);
  }

  log_exit("crbend_tracking");
  return(i_top+1);
}


/* BETA is 2^(1/3) */
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
  long i_kick, step, steps, imult;
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

#if defined(IEEE_MATH)
  if (isnan(x) || isnan(xp) || isnan(y) || isnan(yp)) {
    return 0;
  }
#endif
  if (FABS(x)>COORD_LIMIT || FABS(y)>COORD_LIMIT ||
      FABS(xp)>SLOPE_LIMIT || FABS(yp)>SLOPE_LIMIT) {
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
    for (imult=0; imult<edgeMultData->orders; imult++) {
      apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, xpow, ypow, 
                                      edgeMultData->order[imult], 
                                      edgeMultData->KnL[imult], 0);
      apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, xpow, ypow, 
                                      edgeMultData->order[imult], 
                                      edgeMultData->JnL[imult], 1);
    }
  }
  /* We must do this to avoid numerical precision issues that may subtly change the results
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
        for (imult=0; imult<multData->orders; imult++) {
          if (multData->KnL && multData->KnL[imult]) {
            apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, xpow, ypow, 
                                            multData->order[imult], 
                                            multData->KnL[imult]*kickFrac[step]/n_parts,
                                            0);
          }
          if (multData->JnL && multData->JnL[imult]) {
            apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, xpow, ypow, 
                                            multData->order[imult], 
                                            multData->JnL[imult]*kickFrac[step]/n_parts,
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
	double deltaFactor, F2, dsFactor, dsISRFactor;
        qx /= (1+dp);
        qy /= (1+dp);
	deltaFactor = sqr(1+dp);
	F2 = (sqr(sum_Fy)+sqr(sum_Fx))*sqr(K0L/drift);
	dsFactor = sqrt(1+sqr(xp)+sqr(yp));
	dsISRFactor = dsFactor*drift/3;   /* recall that kickFrac may be negative */
	dsFactor *= drift*kickFrac[step]; /* that's ok here, since we don't take sqrt */
	if (rad_coef)
	  dp -= rad_coef*deltaFactor*F2*dsFactor;
	if (isr_coef>0)
	  dp -= isr_coef*deltaFactor*pow(F2, 0.75)*sqrt(dsISRFactor)*gauss_rn_lim(0.0, 1.0, srGaussianLimit, random_2);
        if (sigmaDelta2)
          *sigmaDelta2 += sqr(isr_coef*deltaFactor)*pow(F2, 1.5)*dsFactor;
        qx *= (1+dp);
        qy *= (1+dp);
      }
    }
  }
  
  if (apData && !checkMultAperture(x+dx, y+dy, apData))  {
    coord[0] = x;
    coord[2] = y;
    return 0;
  }
  
  if (edgeMultData && edgeMultData->orders) {
    fillPowerArray(x, xpow, maxOrder);
    fillPowerArray(y, ypow, maxOrder);
    for (imult=0; imult<edgeMultData->orders; imult++) {
      apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, xpow, ypow, 
                                      edgeMultData->order[imult], 
                                      edgeMultData->KnL[imult], 0);
      apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, xpow, ypow, 
                                      edgeMultData->order[imult], 
                                      edgeMultData->JnL[imult], 1);
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

#if defined(IEEE_MATH)
  if (isnan(x) || isnan(xp) || isnan(y) || isnan(yp)) {
    return 0;
  }
#endif
  if (FABS(x)>COORD_LIMIT || FABS(y)>COORD_LIMIT ||
      FABS(xp)>SLOPE_LIMIT || FABS(yp)>SLOPE_LIMIT) {
    return 0;
  }
  return 1;
}

void switchRbendPlane(double **particle, long n_part, double alpha, double Po)
/* Transforms the reference plane to one that is at an angle alpha relative to the
 * initial plane. Use alpha>0 at the entrance to an RBEND, alpha<0 at the exit.
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
  double C;
  C = sin(alpha)/rho0;
  for (i=0; i<n_part; i++) 
    particle[i][3] -= particle[i][2]*C/(1+particle[i][5]);
}
