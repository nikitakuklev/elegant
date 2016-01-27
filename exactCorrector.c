/*************************************************************************\
* Copyright (c) 2015 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2015 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: exactCorrector.c
 * contents:  track_through_exact_corrector()
 *
 *
 * Michael Borland, 2015
 */
#include "mdb.h"
#include "track.h"

long trackThroughExactCorrector(double **part, long n_part, ELEMENT_LIST  *eptr, double Po, double **accepted, double z_start, double *sigmaDelta2)
{
  long i_part, i_top;
  double xkick, ykick, kick, tilt, length, *coord, rho0, theta0, rho, theta, alpha, arg;
  long isr, sr;
  EHCOR *ehcor;
  EVCOR *evcor;
  EHVCOR *ehvcor;
  void *pePtr;

  pePtr = eptr->p_elem;
  
  xkick = ykick = kick = tilt = length = 0;
  isr = sr = 0;
  switch (eptr->type) {
  case T_EHCOR:
    ehcor = (EHCOR*) pePtr;
    length = ehcor->length;
    xkick = ehcor->kick*ehcor->calibration;
    tilt = ehcor->tilt;
    isr = ehcor->isr;
    sr = ehcor->synchRad;
    break;
  case T_EVCOR:
    evcor = (EVCOR*) pePtr;
    length = evcor->length;
    ykick = evcor->kick*evcor->calibration;
    tilt = evcor->tilt;
    isr = evcor->isr;
    sr = evcor->synchRad;
    break;
  case T_EHVCOR:
    ehvcor = (EHVCOR*) pePtr;
    length = ehvcor->length;
    xkick = ehvcor->xkick*ehvcor->xcalibration;
    ykick = ehvcor->ykick*ehvcor->ycalibration;
    tilt = ehvcor->tilt;
    isr = ehvcor->isr;
    sr = ehvcor->synchRad;
    break;
  default:
    break;
  }

  /* xkick and ykick are actually angles, not kicks */
  theta0 = atan(sqrt(sqr(tan(xkick)) + sqr(tan(ykick))));
  if (theta0==0) {
    exactDrift(part, n_part, length);
    return n_part;
  }
  /* the tilt will account for the direction of bending */
  /* theta0 and rho0 are positive values */
  tilt += atan2(tan(ykick), tan(xkick));
  rho0 = length/sin(theta0);
  
  if (sigmaDelta2)
    *sigmaDelta2 = 0;

  i_top = -1;
  if (isSlave || !notSinglePart) {
    i_top = n_part-1;

    /* Rotate particle coordinates into frame such that bending is all in the "x" coordinates */
    rotateBeamCoordinates(part, n_part, tilt);

    for (i_part=0; i_part<=i_top; i_part++) {
      if (!(coord = part[i_part])) {
        fprintf(stdout, "error: null coordinate pointer for particle %ld (trackThroughExactCorrector)\n", i_part);
        fflush(stdout);
        abort();
      }

      if (rho0==0 || length==0)
        coord[1] = tan(theta0/(1+coord[5]) + atan(coord[1]));
      else {
        alpha = atan(coord[1]);
        rho = rho0*(1+coord[5]);
        if (rho>0 && (arg = length/rho + sin(alpha))<=1 && arg >= -1) {
          double dyf;
          theta = asin(arg) - alpha;
          dyf = coord[3]/sqrt(1 + sqr(coord[1]) - sqr(coord[3]));
          coord[2] += dyf*theta*rho;
          coord[4] += theta*rho*sqrt(1 + sqr(dyf));
          coord[0] += rho*(cos(alpha) - cos(theta+alpha)); 
          coord[1]  = tan(theta + alpha);
        } else {
          /* particle is lost */
          if (accepted && !accepted[i_part]) {
            fprintf(stdout, "error: null accepted particle pointer for particle %ld (track_through_csbend)\n", i_part);
            fflush(stdout);
            abort();
          }
          if (!part[i_top]) {
            fprintf(stdout, "error: couldn't swap particles %ld and %ld--latter is null pointer (track_through_csbend)\n",
                    i_part, i_top);
            fflush(stdout);
            abort();
          }
          swapParticles(part[i_part], part[i_top]);
          if (accepted) {
            if (!accepted[i_top]) {
              fprintf(stdout, 
                      "error: couldn't swap acceptance data for particles %ld and %ld--latter is null pointer (track_through_csbend)\n",
                      i_part, i_top);
              fflush(stdout);
              abort();
            }
            swapParticles(accepted[i_part], accepted[i_top]);
          }
          part[i_top][4] = z_start;
          part[i_top][5] = Po*(1+part[i_top][5]);
          i_top--;
          i_part--;
        }
      }
    }
    
    /* note that we use n_part here instead of i_top, so the lost particles get rotated also */
    rotateBeamCoordinates(part, n_part, -tilt);
  }
  
  if (sr)
    addCorrectorRadiationKick(part, i_top, eptr, eptr->type, Po, sigmaDelta2, 0);
  
  return i_top+1;
}

