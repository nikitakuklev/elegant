/************************************************************************* \
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

#include <ctype.h>
#include "mdb.h"
#include "track.h"
#include "sort.h"
#ifdef HAVE_GPU
#include <gpu_limit_amplitudes.h>
#endif

long evaluateLostWithOpenSides(long code, double dx, double dy, double xsize, double ysize);

/* routine: rectangular_collimator()
 * purpose: eliminate particles that hit the walls of a non-zero length
 *          rectangular collimator
 *
 * Michael Borland, 1989
 */

#ifdef VAX_VMS
#define isnan(x) 0
#define isinf(x) 0
#endif

long rectangular_collimator(
                            double **initial, RCOL *rcol, long np, double **accepted, double z,
                            double Po
                            )
{
  double length, *ini;
  long ip, itop, is_out, lost, openCode;
  double xsize, ysize;
  double x_center, y_center;
  double x1, y1, zx, zy, dx, dy;

#ifdef HAVE_GPU
   if(getElementOnGpu()){
      startGpuTimer();
      ip = gpu_rectangular_collimator(rcol, np, accepted, z, Po);
#ifdef GPU_VERIFY     
      startCpuTimer();
      rectangular_collimator(initial, rcol, np, accepted, z, Po);
      compareGpuCpu(np, "rectangular_collimator");
#endif /* GPU_VERIFY */
      return ip;
    }
#endif /* HAVE_GPU */

  xsize  = rcol->x_max;
  ysize  = rcol->y_max;
  x_center = rcol->dx;
  y_center = rcol->dy;

  /*
  if (rcol->invert && rcol->length) {
    TRACKING_CONTEXT tc;
    getTrackingContext(&tc);
    fprintf(stderr, "Problem for %s#%ld:\n", tc.elementName, tc.elementOccurrence);
    bombElegant("Cannot have invert=1 and non-zero length for RCOL", NULL);
  }
  */

  if (xsize<=0 && ysize<=0) {
    exactDrift(initial, np, rcol->length);
    return(np);
  }
  openCode = determineOpenSideCode(rcol->openSide);
  itop = np-1;
  for (ip=0; ip<np; ip++) {
    ini = initial[ip];
    dx = ini[0] - x_center;
    dy = ini[2] - y_center;
    lost = 0;
    if ((xsize>0 && fabs(dx) > xsize) ||
        (ysize>0 && fabs(dy) > ysize)) {
      lost = openCode ? evaluateLostWithOpenSides(openCode, dx, dy, xsize, ysize) : 1;
    } else if (isinf(ini[0]) || isinf(ini[2]) ||
               isnan(ini[0]) || isnan(ini[2]) )
      lost = 1;
    if (rcol->invert)
      lost = !lost;
    if (lost) {
      swapParticles(initial[ip], initial[itop]);
     if (accepted)
        swapParticles(accepted[ip], accepted[itop]);
      initial[itop][4] = z; /* record position of particle loss */
      initial[itop][5] = Po*(1+initial[itop][5]);
      --itop;
      --ip;
      --np;
    }
  }
  if (np==0 || (length=rcol->length)<=0) {
    return(np);
  }

  itop = np-1;
  for (ip=0; ip<np; ip++) {
    ini = initial[ip];
    x1 = ini[0] + length*ini[1];
    y1 = ini[2] + length*ini[3];
    dx = x1 - x_center;
    dy = y1 - y_center;
    is_out = 0;
    if (xsize>0 && fabs(dx)>xsize)
      is_out += 1*(openCode?evaluateLostWithOpenSides(openCode, dx, 0, xsize, ysize):1);
    if (ysize>0 && fabs(dy)>ysize)
      is_out += 2*(openCode?evaluateLostWithOpenSides(openCode, 0, dy, xsize, ysize):1);
    if (isinf(x1) || isinf(y1) || isnan(x1) || isnan(y1) )
      is_out += 4;
    if (is_out&4) {
      ini[4] = z+length;
      ini[0] = x1;
      ini[2] = y1;
      ini[5] = Po*(1+ini[5]);
      swapParticles(initial[ip], initial[itop]);
      if (accepted)
        swapParticles(accepted[ip], accepted[itop]);
      --itop;
      --ip;
      --np;
    } else if ((is_out && !rcol->invert) || (!is_out && rcol->invert)) {
      if (!openCode) {
        zx = zy = DBL_MAX;
        if (is_out&1 && ini[1]!=0)
          zx = (SIGN(ini[1])*xsize-(ini[0]-x_center))/ini[1];
        if (is_out&2 && ini[3]!=0)
          zy = (SIGN(ini[3])*ysize-(ini[2]-y_center))/ini[3];
        if (zx<zy) {
          ini[0] += ini[1]*zx;
          ini[2] += ini[3]*zx;
          ini[4] = z+zx; 
        }
        else {
          ini[0] += ini[1]*zy;
          ini[2] += ini[3]*zy;
          ini[4] = z+zy;
        }
      }
      ini[5] = Po*(1+ini[5]);
      swapParticles(initial[ip], initial[itop]);
      if (accepted)
        swapParticles(accepted[ip], accepted[itop]);
      --itop;
      --ip;
      --np;
    }
    else {
      ini[4] += length*sqrt(1+sqr(ini[1])+sqr(ini[3]));
      ini[0] = x1;
      ini[2] = y1;
    }
  }
  return(np);
}


/* routine: limit_amplitudes()
 * purpose: eliminate particles with (x,y) larger than given values
 *
 * Michael Borland, 1989
 */

long limit_amplitudes(
    double **coord, double xmax, double ymax, long np, double **accepted,
    double z, double Po, long extrapolate_z, long openCode)
{
    long ip, itop, is_out;
    double *part;
    double dz, dzx, dzy;

#ifdef HAVE_GPU
   if(getElementOnGpu()){
      startGpuTimer();
      ip = gpu_limit_amplitudes(xmax, ymax, np, accepted, z, Po, 
                                extrapolate_z, openCode);
#ifdef GPU_VERIFY     
      startCpuTimer();
      limit_amplitudes(coord, xmax, ymax, np, accepted, z, Po,
                       extrapolate_z, openCode);
      compareGpuCpu(np, "limit_amplitudes");
#endif /* GPU_VERIFY */
      return ip;
    }
#endif /* HAVE_GPU */

    if (xmax<0 && ymax<0) {
      return(np);
    }

    itop = np-1;

    for (ip=0; ip<np; ip++) {
        part = coord[ip];
        is_out = 0;
        if (xmax>0 && fabs(part[0])>xmax)
            is_out += 1;
        if (ymax>0 && fabs(part[2])>ymax)
            is_out += 2;
        if (openCode)
          is_out *= evaluateLostWithOpenSides(openCode, part[0], part[2], xmax, ymax);
        if (isinf(part[0]) || isinf(part[2]) || isnan(part[0]) || isnan(part[2]) )
            is_out += 4;
        dz = 0;
        if (is_out && !(is_out&4) && !openCode && extrapolate_z) {
            /* find the actual position of loss, assuming a drift preceded with 
             * the same aperture 
             */
            dzx = dzy = -DBL_MAX;
            if (is_out&1 && part[1]!=0)
                dzx = (part[0]-SIGN(part[1])*xmax)/part[1];
            if (is_out&2 && part[3]!=0)
                dzy = (part[2]-SIGN(part[3])*ymax)/part[3];
            if (dzx>dzy) 
                dz = -dzx;
            else
                dz = -dzy;
            if (dz==-DBL_MAX)
                dz = 0;
            part[0] += dz*part[1];
            part[2] += dz*part[3];
            }
        if (is_out) {
          if (ip!=itop)
            swapParticles(coord[ip], coord[itop]);
          coord[itop][4] = z+dz;  /* record position of loss */
          coord[itop][5] = Po*(1+coord[itop][5]);
          if (accepted)
            swapParticles(accepted[ip], accepted[itop]);
          --itop;
          --ip;
          np--;
        }
      }
    return(np);
  }


long removeInvalidParticles(
    double **coord, long np, double **accepted,
    double z, double Po)
{
    long ip, itop, is_out, ic;
    double *part;

#ifdef HAVE_GPU
   if(getElementOnGpu()){
      startGpuTimer();
      ip = gpu_removeInvalidParticles(np, accepted, z, Po);
#ifdef GPU_VERIFY     
      startCpuTimer();
      removeInvalidParticles(coord, np, accepted, z, Po);
      compareGpuCpu(np, "removeInvalidParticles");
#endif /* GPU_VERIFY */
      return ip;
    }
#endif /* HAVE_GPU */

    itop = np-1;
    for (ip=0; ip<np; ip++) {
        part = coord[ip];
        is_out = 0;
        for (ic=0; ic<6; ic++)
          if (isnan(part[ic])) {
            is_out = 1;
            break;
          }
        if (part[5]<=-1)
          is_out = 1;
        if (is_out) {
          swapParticles(coord[ip], coord[itop]);
          coord[itop][4] = z; 
          coord[itop][5] = Po*(1+coord[itop][5]);
          if (accepted)
            swapParticles(accepted[ip], accepted[itop]);
          --itop;
          --ip;
          np--;
        }
      }
    return(np);
  }

            
/* routine: elliptical_collimator()
 * purpose: eliminate particles that hit the walls of a non-zero length
 *          elliptical collimator
 *
 * Michael Borland, 1989
 */

long elliptical_collimator(
                           double **initial, ECOL *ecol, long np, double **accepted, double z,
                           double Po)
{
  double length, *ini;
  long ip, itop, lost, openCode;
  double a2, b2;
  double dx, dy, xo, yo, xsize, ysize;
  TRACKING_CONTEXT context;
  long xe, ye;
  
#ifdef HAVE_GPU
   if(getElementOnGpu()){
      startGpuTimer();
      ip = gpu_elliptical_collimator(ecol, np, accepted, z, Po);
#ifdef GPU_VERIFY     
      startCpuTimer();
      elliptical_collimator(initial, ecol, np, accepted, z, Po);
      compareGpuCpu(np, "elliptical_collimator");
#endif /* GPU_VERIFY */
      return ip;
    }
#endif /* HAVE_GPU */

  xsize = ecol->x_max;
  ysize = ecol->y_max;
  if ((xe=ecol->exponent)<2 || xe%2) {
    getTrackingContext(&context);
    fprintf(stderr, "Error for %s: exponent=%ld is not valid.  Give even integer >=2\n",
            context.elementName, xe);
    exitElegant(1);
  }
  ye = xe;
  if (ecol->yExponent) {
    if ((ye=ecol->yExponent)<2 || ye%2) {
      getTrackingContext(&context);
      fprintf(stderr, "Error for %s: exponent=%ld is not valid.  Give even integer >=2\n",
              context.elementName, ye);
      exitElegant(1);
    }
  }
  
  a2 = ipow(ecol->x_max, xe);
  b2 = ipow(ecol->y_max, ye);
  dx = ecol->dx;
  dy = ecol->dy;

  if (ecol->x_max<=0 || ecol->y_max<=0) {
    /* At least one of x_max or y_max is non-positive */
    if (ecol->x_max>0 || ecol->y_max>0) {
      /* One of x_max or y_max is positive. Use rectangular collimator routine to implement this. */
      RCOL rcol;
      rcol.length = ecol->length;
      rcol.x_max = ecol->x_max;
      rcol.y_max = ecol->y_max;
      rcol.dx = ecol->dx;
      rcol.dy = ecol->dy;
      rcol.invert = ecol->invert;
      rcol.openSide = ecol->openSide;
      return rectangular_collimator(initial, &rcol, np, accepted, z, Po);
    } 
    exactDrift(initial, np, ecol->length);
    return(np);
  }
  openCode = determineOpenSideCode(ecol->openSide);

  itop = np-1;
  for (ip=0; ip<np; ip++) {
    ini = initial[ip];
    lost = 0;
    xo = ini[0] - dx;
    yo = ini[2] - dy;
    if ((ipow(xo, xe)/a2 + ipow(yo, ye)/b2)>1)
      lost = openCode ? evaluateLostWithOpenSides(openCode, xo, yo, xsize, ysize) : 1;
    else if (isinf(ini[0]) || isinf(ini[2]) ||
             isnan(ini[0]) || isnan(ini[2]) )
      lost = 1;
    if (ecol->invert)
      lost = !lost;
    if (lost) {
      swapParticles(initial[ip], initial[itop]);
      initial[itop][4] = z;
      initial[itop][5] = sqrt(sqr(Po*(1+initial[itop][5]))+1);
      if (accepted)
        swapParticles(accepted[ip], accepted[itop]);
      --itop;
      --ip;
      np--;
    }
  }

  if (np==0 || (length=ecol->length)<=0)
    return(np);
  
  itop = np-1;
  for (ip=0; ip<np; ip++) {
    ini = initial[ip];
    lost = 0;
    ini[0] += length*ini[1];
    ini[2] += length*ini[3];
    xo = (ini[0]-dx)/xsize;
    yo = (ini[2]-dy)/ysize;
    if ((ipow(xo, ecol->exponent) + ipow(yo, ecol->exponent))>1)
      lost = openCode ? evaluateLostWithOpenSides(openCode, xo, yo, 1, 1) : 1;
    else if (isinf(ini[0]) || isinf(ini[2]) ||
             isnan(ini[0]) || isnan(ini[2]) )
      lost = 1;
    if (ecol->invert)
      lost = !lost;
    if (lost) {
      swapParticles(initial[ip], initial[itop]);
      initial[itop][4] = z + length;
      initial[itop][5] = sqrt(sqr(Po*(1+initial[itop][5]))+1);
      if (accepted)
        swapParticles(accepted[ip], accepted[itop]);
      --itop;
      --ip;
      np--;
    }
    else 
      ini[4] += length*sqrt(1+sqr(ini[1])+sqr(ini[3]));
  }

  return(np);
}



/* routine: elimit_amplitudes()
 * purpose: eliminate particles outside an ellipse with given semi-major
 *          and semi-minor axes.
 *
 * Michael Borland, 1989
 */

long elimit_amplitudes(
    double **coord,
    double xmax,       /* half-axis in x direction */ 
    double ymax,       /* half-axis in y direction */
    long np,
    double **accepted,
    double z,
    double Po,
    long extrapolate_z,
    long openCode,
    long exponent,
    long yexponent
    )
{
    long ip, itop, lost;
    double *part;
    double a2, b2, c1, c2, c0, dz, det;
    TRACKING_CONTEXT context;
    long xe, ye;
    
#ifdef HAVE_GPU
   if(getElementOnGpu()){
      startGpuTimer();
      ip = gpu_elimit_amplitudes(xmax, ymax, np, accepted, z, Po, 
                                extrapolate_z, openCode, exponent, yexponent);
#ifdef GPU_VERIFY     
      startCpuTimer();
      elimit_amplitudes(coord, xmax, ymax, np, accepted, z, Po,
                       extrapolate_z, openCode, exponent, yexponent);
      compareGpuCpu(np, "elimit_amplitudes");
#endif /* GPU_VERIFY */
      return ip;
    }
#endif /* HAVE_GPU */

    if ((xe=exponent)<2 || xe%2) {
      getTrackingContext(&context);
      fprintf(stderr, "Error for %s: exponent=%ld is not valid.  Give even integer >=2\n",
              context.elementName, xe);
      exitElegant(1);
    }
    ye = xe;
    if (yexponent) {
      if ((ye=yexponent)<2 || ye%2) {
        getTrackingContext(&context);
        fprintf(stderr, "Error for %s: exponent=%ld is not valid.  Give even integer >=2\n",
                context.elementName, ye);
        exitElegant(1);
      }
    }

    if (xmax<=0 || ymax<=0) {
      /* At least one of the dimensions is non-positive and therefore ignored */
      if (xmax>0 || ymax>0) 
        return limit_amplitudes(coord, xmax, ymax, np, accepted, z, Po, extrapolate_z, openCode);
      return(np);
    }

    a2 = ipow(xmax, xe);
    b2 = ipow(ymax, ye);
    itop = np-1;

    for (ip=0; ip<np; ip++) {
        part = coord[ip];
        if (isinf(part[0]) || isinf(part[2]) || isnan(part[0]) || isnan(part[2]) ) {
            swapParticles(coord[ip], coord[itop]);
            coord[itop][4] = z;
            coord[itop][5] = sqrt(sqr(Po*(1+coord[itop][5]))+1);
            if (accepted)
                swapParticles(accepted[ip], accepted[itop]);
            --itop;
            --ip;
            np--;
            continue;
            }
        lost = 0;
        if ((ipow(part[0], xe)/a2 + ipow(part[2], ye)/b2) > 1)
          lost = 1;
        else
          continue;
        if (openCode)
          lost *= evaluateLostWithOpenSides(openCode, part[0], part[2], xmax, ymax);
        if (lost) {
          dz = 0;
          if (extrapolate_z && !openCode && xe==2 && ye==2) {
            c0 = sqr(part[0])/a2 + sqr(part[2])/b2 - 1;
            c1 = 2*(part[0]*part[1]/a2 + part[2]*part[3]/b2);
            c2 = sqr(part[1])/a2 + sqr(part[3])/b2;
            det = sqr(c1)-4*c0*c2;
            if (z>0 && c2 && det>=0) {
              if ((dz = (-c1+sqrt(det))/(2*c2))>0)
                dz = (-c1-sqrt(det))/(2*c2);
              if ((z+dz)<0) 
                dz = -z;
              part[0] += dz*part[1];
              part[2] += dz*part[3];
            }
          }
          swapParticles(coord[ip], coord[itop]);
          coord[itop][4] = z+dz;
          coord[itop][5] = sqrt(sqr(Po*(1+coord[itop][5]))+1);
          if (accepted)
            swapParticles(accepted[ip], accepted[itop]);
          --itop;
          --ip;
          np--;
        }
      }
    log_exit("elimit_amplitudes");
    return(np);
  }

#define SQRT_3 1.7320508075688772
            
long beam_scraper(
                  double **initial, SCRAPER *scraper, long np, double **accepted, double z,
                  double Po
                  )
{
  double length, *ini;
  long do_x, do_y, ip, idir, hit;
  long dsign[2] = {1, -1};
  long dflag[2] = {0, 0};
  MATTER matter;

  if (!((scraper->direction = interpretScraperDirection(scraper->insert_from, scraper->oldDirection))&(DIRECTION_X|DIRECTION_Y))) {
    exactDrift(initial, np, scraper->length);
    return np;
  }

#ifdef HAVE_GPU
   if(getElementOnGpu()){
      startGpuTimer();
      ip = gpu_beam_scraper(scraper, np, accepted, z, Po); 
#ifdef GPU_VERIFY     
      startCpuTimer();
      beam_scraper(initial, scraper, np, accepted, z, Po); 
      compareGpuCpu(np, "beam_scraper");
#endif /* GPU_VERIFY */
      return ip;
    }
#endif /* HAVE_GPU */

  log_entry("beam_scraper");

  do_x = do_y = 0;
  if (scraper->direction&DIRECTION_X) {
    /* For now, x and +x are the same */
    do_x = 1;
    dflag[0] = scraper->direction&DIRECTION_PLUS_X ? 1 : 0;
    dflag[1] = scraper->direction&DIRECTION_MINUS_X ? 1 : 0;
  }
  else if (scraper->direction&DIRECTION_Y) {
    do_y = 1;
    dflag[0] = scraper->direction&DIRECTION_PLUS_Y ? 1 : 0;
    dflag[1] = scraper->direction&DIRECTION_MINUS_Y ? 1 : 0;
  } else
    return np;

  if (scraper->length && (scraper->Xo || scraper->Z)) {
    /* scraper has material properties that scatter beam and
     * absorb energy
     */
    matter.length = scraper->length;
    matter.lEffective = 0;
    matter.Xo = scraper->Xo;
    matter.energyDecay = scraper->energyDecay;
    matter.energyStraggle = scraper->energyStraggle;
    matter.nuclearBremsstrahlung = scraper->nuclearBremsstrahlung;
    matter.electronRecoil = scraper->electronRecoil;
    matter.Z = scraper->Z;
    matter.A = scraper->A;
    matter.rho = scraper->rho;
    matter.pLimit = scraper->pLimit;
    matter.width = matter.spacing = matter.tilt = matter.center = 0;
    matter.nSlots = 0;
    matter.startPass = matter.endPass = -1;

    for (ip=0; ip<np; ip++) {
      ini = initial[ip];
      hit = 0;
      if (dflag[0] && dflag[1]) {
        if ((do_x && fabs(ini[0]-scraper->dx)>scraper->position) ||
            (do_y && fabs(ini[2]-scraper->dy)>scraper->position)) {
          /* scatter and/or absorb energy */
          hit = 1;
          if (!track_through_matter(&ini, 1, 0, &matter, Po, NULL, z))
            ini[5] = -1;
        }
      } else {
        for (idir=0; idir<2 && !hit; idir++) {
          if (!dflag[idir]) continue;
          if ((do_x && dsign[idir]*(ini[0]-scraper->dx)>dsign[idir]*scraper->position) ||
              (do_y && dsign[idir]*(ini[2]-scraper->dy)>dsign[idir]*scraper->position)) {
            /* scatter and/or absorb energy */
            hit = 1;
            if (!track_through_matter(&ini, 1, 0, &matter, Po, NULL, z))
              ini[5] = -1;
          }
        }
      }
      if (!hit) {
        ini[0] = ini[0] + ini[1]*scraper->length;
        ini[2] = ini[2] + ini[3]*scraper->length;
        ini[4] += scraper->length*sqrt(1+sqr(ini[1])+sqr(ini[3]));
      }
    }
    log_exit("beam_scraper");
    return(np);
  }

  if (dflag[0] && dflag[1]) {
    /* Scraper from both sides is like an RCOL */
    RCOL rcol;
    if (scraper->position<=0) {
      for (ip=0; ip<np; ip++) {
        initial[ip][4] = z;
        initial[ip][5] = Po*(initial[ip][5]);
      }
      return 0;
    }
    rcol.length = scraper->length;
    rcol.x_max = rcol.y_max = 0;
    if (do_x)
      rcol.x_max = scraper->position;
    else
      rcol.y_max = scraper->position;
    rcol.dx = scraper->dx;
    rcol.dy = scraper->dy;
    rcol.invert = 0;
    rcol.openSide = NULL;
    return rectangular_collimator(initial, &rcol, np, accepted, z, Po);
  }

  /* come here for idealized on-sided scraper that just absorbs particles */
  for (ip=0; ip<np; ip++) {
    ini = initial[ip];
    hit = 0;
    for (idir=0; idir<2 && !hit; idir++) {
      if (!dflag[idir]) continue;
      if ((do_x && dsign[idir]*(ini[0]-scraper->dx)>dsign[idir]*scraper->position) ||
          (do_y && dsign[idir]*(ini[2]-scraper->dy)>dsign[idir]*scraper->position)) {
        swapParticles(initial[ip], initial[np-1]);
        if (accepted)
          swapParticles(accepted[ip], accepted[np-1]);
        initial[np-1][4] = z; /* record position of particle loss */
        initial[np-1][5] = Po*(1+initial[np-1][5]);
        --ip;
        --np;
        hit = 1;
      }
    }
  }
  if (np==0 || (length=scraper->length)<=0) {
    log_exit("beam_scraper");
    return(np);
  }

  z += length;
  for (ip=0; ip<np; ip++) {
    ini = initial[ip];
    ini[0] += length*ini[1];
    ini[2] += length*ini[3];
    hit = 0;
    for (idir=0; idir<2 && !hit; idir++) {
      if (!dflag[idir]) continue;
      if ((do_x && dsign[idir]*(ini[0]-scraper->dx)>dsign[idir]*scraper->position) ||
          (do_y && dsign[idir]*(ini[2]-scraper->dy)>dsign[idir]*scraper->position)) {
        double dz;
        if (do_x) {
          double dx;
          dx = ini[0]-scraper->dx-scraper->position;
          dz = dx/ini[1];
        } else {
          double dy;
          dy = ini[2]-scraper->dy-scraper->position;
          dz = dy/ini[3];
        }
        ini[4] = z-dz;
        ini[0] -= ini[1]*dz;
        ini[2] -= ini[3]*dz;
        ini[5] = Po*(1+ini[5]);
        swapParticles(initial[ip], initial[np-1]);
        if (accepted)
          swapParticles(accepted[ip], accepted[np-1]);
        --ip;
        --np;
        hit = 1;
      }
    }
    if (!hit)
      ini[4] += length*sqrt(1+sqr(ini[1])+sqr(ini[3]));
  }
  log_exit("beam_scraper");
  return(np);
}

long track_through_pfilter(
    double **initial, PFILTER *pfilter, long np, double **accepted, double z,
    double Po
    )
{
  long ip, itop;
  static double *deltaBuffer=NULL;
  static long maxBuffer = 0;
  double reference;
#ifdef USE_KAHAN
  double error = 0.0; 
#endif
  
  itop = np-1;
  if ((pfilter->lowerFraction || pfilter->upperFraction) && !pfilter->limitsFixed) {
    double level[2]={-1,-1}, limit[2], upper[2];
    long count = 0, i;
    if (maxBuffer<np &&
        !(deltaBuffer=SDDS_Realloc(deltaBuffer, sizeof(*deltaBuffer)*(maxBuffer=np))))
      SDDS_Bomb("memory allocation failure");  
    if (isSlave || !notSinglePart)
    for (ip=0; ip<np; ip++)
      deltaBuffer[ip] = initial[ip][5];
    /* eliminate lowest lowerfraction of particles and highest
       upperfraction */
    if (pfilter->lowerFraction>0 && pfilter->lowerFraction<1) {
      upper[count] = 0;
      level[count++] = pfilter->lowerFraction*100;
    }
    if (pfilter->upperFraction>0 && pfilter->upperFraction<1) {
      upper[count] = 1;
      level[count++] = 100-pfilter->upperFraction*100;
    }
#if SDDS_MPI_IO
    if (notSinglePart)
      approximate_percentiles_p(limit, level, count, deltaBuffer, np, pfilter->bins); 
    else
      compute_percentiles(limit, level, count, deltaBuffer, np);
#else     
    compute_percentiles(limit, level, count, deltaBuffer, np);
#endif

    itop = np-1;
    for (i=0; i<2; i++) {
      if (level[i]<0)
        break;
      if (upper[i]) {
        pfilter->pUpper = (1+limit[i])*Po;
        pfilter->hasUpper = 1;
      }
      else {
        pfilter->pLower = (1+limit[i])*Po;
        pfilter->hasLower = 1;
      }
      if (!pfilter->fixPLimits) {
	/* filter in next block so there are no discrepancies due to
	 * small numerical differences
	 */
	if (isSlave || !notSinglePart)
	for (ip=0; ip<=itop; ip++) {
	  if ((upper[i] && initial[ip][5]>limit[i]) ||
	      (!upper[i] && initial[ip][5]<limit[i])) {
	    swapParticles(initial[ip], initial[itop]);
	    initial[itop][4] = z;  /* record position of particle loss */
	    initial[itop][5] = Po*(1+initial[itop][5]);  /* momentum at loss */
	    if (accepted)
	      swapParticles(accepted[ip], accepted[itop]);
	    --itop;
	    --ip;
	  }
	}
      }
    }
    if (pfilter->fixPLimits)
      pfilter->limitsFixed = 1;
  }
  if (pfilter->limitsFixed) {
    double p;
    if(isSlave || !notSinglePart)
    for (ip=0; ip<=itop; ip++) {
      p = (1+initial[ip][5])*Po;
      if ((pfilter->hasUpper && p>pfilter->pUpper) ||
	  (pfilter->hasLower && p<pfilter->pLower)) {
	swapParticles(initial[ip], initial[itop]);
	initial[itop][4] = z;  /* record position of particle loss */
	initial[itop][5] = p;  /* momentum at loss */
	if (accepted)
	  swapParticles(accepted[ip], accepted[itop]);
	--itop;
	--ip;
      }
    }
  }
  
  if (pfilter->deltaLimit<0) {
#if defined(MINIMIZE_MEMORY)
    free(deltaBuffer);
    deltaBuffer = NULL;
    maxBuffer= 0;
#endif
    return itop+1;
  }
  reference = 0.0;
  if (pfilter->beamCentered) {
    if (isSlave)
    for (ip=0; ip<=itop; ip++) {
#ifndef USE_KAHAN
      reference += initial[ip][5];
#else
      reference = KahanPlus(reference, initial[ip][5], &error);
#endif
    }
#if USE_MPI
    if (notSinglePart) {
      if (USE_MPI) {
	long itop_total;
	double reference_total;
	MPI_Allreduce(&itop, &itop_total, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&reference, &reference_total, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	reference = reference_total/(itop_total+n_processors);
      }
    } else
      reference /= (itop+1);
#else
    reference /= (itop+1);
#endif
  }
  if (isSlave || !notSinglePart)
  for (ip=0; ip<=itop; ip++) {
    if (fabs(initial[ip][5]-reference)<pfilter->deltaLimit)
      continue;
    swapParticles(initial[ip], initial[itop]);
    initial[itop][4] = z;  /* record position of particle loss */
    initial[itop][5] = Po*(1+initial[itop][5]);  /* momentum at loss */
    if (accepted)
      swapParticles(accepted[ip], accepted[itop]);
    --itop;
    --ip;
  }
#if defined(MINIMIZE_MEMORY)
  free(deltaBuffer);
  deltaBuffer = NULL;
  maxBuffer= 0;
#endif
  return itop+1;
}

long remove_outlier_particles(
                              double **initial, CLEAN *clean, long np, double **accepted, double z,
                              double Po
                              )
{
  double *ini, beta, p, *sSave;
  long ip, itop, is_out, j, mode;
  double limit[6], centroid[6], stDev[6];
  long count;
#define CLEAN_STDEV 0
#define CLEAN_ABSDEV 1
#define CLEAN_ABSVAL 2
#define CLEAN_MODES 3
  static char *modeName[CLEAN_MODES] = {"stdeviation", "absdeviation", "absvalue"};

  switch ((mode=match_string(clean->mode, modeName, 3, 0))) {
  case CLEAN_STDEV:
  case CLEAN_ABSDEV:
  case CLEAN_ABSVAL:
    break;
  default:
    fprintf(stderr, "Error: mode for CLEAN element must be one of the following:\n");
    for (j=0; j<CLEAN_MODES; j++)
      fprintf(stderr, "%s  ", modeName[j]);
    fprintf(stderr, "\n");
    exitElegant(1);
    break;
  }
  
  if (!(clean->xLimit>0 || clean->xpLimit>0 ||
        clean->yLimit>0 || clean->ypLimit>0 ||
        clean->tLimit>0 || clean->deltaLimit>0)) {
    return(np);
  }

  /* copy limits into array */
  limit[0] = clean->xLimit;
  limit[1] = clean->xpLimit;
  limit[2] = clean->yLimit;
  limit[3] = clean->ypLimit;
  limit[4] = clean->tLimit;
  limit[5] = clean->deltaLimit;

  if (!(sSave = malloc(sizeof(*sSave)*np)))
    bombElegant("memory allocation failure", NULL);

  /* compute centroids for each coordinate */
  for (j=0; j<6; j++) {
    centroid[j] = stDev[j] = 0;
  }
  for (ip=count=0; ip<np; ip++) {
    ini = initial[ip];
    sSave[ip] = ini[4];
    for (j=0; j<6; j++)
      if (isnan(ini[j]) || isinf(ini[j]))
        break;
    if (j!=6)
      continue;
    /* compute time of flight and store in place of s */
    p = Po*(1+ini[5]);
    beta = p/sqrt(p*p+1);
    ini[4] /= beta*c_mks;
    count++;
    for (j=0; j<6; j++) {
      centroid[j] += ini[j];
    }
  }

#if USE_MPI
  if (isSlave) {
    double sum[6];
    long count_total;

    MPI_Allreduce(centroid,sum,6,MPI_DOUBLE,MPI_SUM,workers);
    memcpy(centroid,sum,sizeof(double)*6);
    MPI_Allreduce(&count,&count_total,1,MPI_LONG,MPI_SUM,workers);
    count = count_total;
  }
#endif

  if (!count) {
    for (ip=0; ip<np; ip++) {
      initial[ip][4] = z; /* record position of particle loss */
      initial[ip][5] = Po*(1+initial[ip][5]);
    }
    free(sSave);
    return 0;
  }
  for (j=0; j<6; j++)
    centroid[j] /= count;

  /* compute standard deviation of coordinates, if needed */
  if (mode==CLEAN_STDEV) {
    for (ip=0; ip<np; ip++) {
      ini = initial[ip];
      for (j=0; j<6; j++) {
        if (!isnan(ini[j]) && !isinf(ini[j]))
          stDev[j] += sqr(centroid[j] - ini[j]);
      }
    }
#if USE_MPI
  if (isSlave) {
    double sum[6];

    MPI_Allreduce(stDev,sum,6,MPI_DOUBLE,MPI_SUM,workers);
    memcpy(stDev,sum,sizeof(double)*6);
  }
#endif
    for (j=0; j<6; j++)
      stDev[j] = sqrt(stDev[j]/count);
  }
  
  itop = np-1;
  for (ip=0; ip<np; ip++) {
    ini = initial[ip];
    for (j=is_out=0; j<6; j++) {
      if (limit[j]<=0)
        continue;
      switch (mode) {
      case CLEAN_STDEV:
        if (fabs(ini[j]-centroid[j])/stDev[j]>limit[j])
          is_out = 1;
        break;
      case CLEAN_ABSDEV:
        if (fabs(ini[j]-centroid[j])>limit[j])
          is_out = 2;
        break;
      case CLEAN_ABSVAL:
        if (fabs(ini[j])>limit[j])
          is_out = 3;
        break;
      default:
        fprintf(stderr, "invalid mode in remove_outlier_particles---programming error!\n");
        exitElegant(1);
        break;
      }
      if (isnan(ini[j]) || isinf(ini[j]))
	is_out = 1;
      if (is_out)
        break;
    }
    if (is_out) {
      swapParticles(initial[ip], initial[itop]);
      SWAP_DOUBLE(sSave[ip], sSave[itop]);
      if (accepted)
        swapParticles(accepted[ip], accepted[itop]);
      initial[itop][4] = z; /* record position of particle loss */
      initial[itop][5] = Po*(1+initial[itop][5]);
      --itop;
      --ip;
      --np;
    }
  }

  /* restore s data to particle array */
  for (ip=0; ip<np; ip++)
    initial[ip][4] = sSave[ip];
  free(sSave);

  return(np);
}


long evaluateLostWithOpenSides(long code, double dx, double dy, double xsize, double ysize)
{
  long lost = 1;
  switch (code) {
  case OPEN_PLUS_X:
    if (dx>0 && fabs(dy)<ysize)
      lost = 0;
    break;
  case OPEN_MINUS_X:
    if (dx<0 && fabs(dy)<ysize)
      lost = 0;
    break;
  case OPEN_PLUS_Y:
    if (dy>0 && fabs(dx)<xsize)
      lost = 0;
    break;
          case OPEN_MINUS_Y:
    if (dy<0 && fabs(dx)<xsize)
      lost = 0;
    break;
  default:
    break;
  }
  return lost;
}

long determineOpenSideCode(char *openSide)
{
  TRACKING_CONTEXT context;
  long value = -1;
  
  if (!openSide)
    return 0;
  if (strlen(openSide)==2 &&
      (openSide[0]=='+' || openSide[0]=='-')) {
    switch (openSide[1]) {
    case 'x':
    case 'X':
    case 'h':
    case 'H':
      value = (openSide[0]=='+' ? OPEN_PLUS_X : OPEN_MINUS_X);
      break;
    case 'y':
    case 'Y':
    case 'v':
    case 'V':
      value = (openSide[0]=='+' ? OPEN_PLUS_Y : OPEN_MINUS_Y);
      break;
    default:
      value = -1;
      break;
    }
  }
  if (value==-1) {
    getTrackingContext(&context);
    fprintf(stderr, "Error for %s: open_side=%s is not valid\n",
            context.elementName, openSide);
    exitElegant(1);
  }
  return value;
}


/* routine: imposeApertureData()
 * purpose: eliminate particles with (x,y) larger than values given in aperture input file.
 *
 * Michael Borland, 2005
 */

#define DEBUG_APERTURE 0

long interpolateApertureData(double z, APERTURE_DATA *apData,
                             double *xCenter, double *yCenter, double *xSize, double *ySize)
{
  double z0, period;
  long iz;
 
  z0 = z;
  if (apData->s[apData->points-1]<z) {
    if (apData->periodic) {
      period = apData->s[apData->points-1]-apData->s[0];
      z0 -= period*((long) ((z0 - apData->s[0])/period));
    } else
      return 0;
  }
  if ((iz = binaryArraySearch(apData->s, sizeof(apData->s[0]), apData->points, &z0, 
                              double_cmpasc, 1))<0)
    return 0;
  if (iz==apData->points-1)
    iz -= 1;

  if (apData->s[iz]==apData->s[iz+1]) {
    *xCenter = apData->dx[iz];
    *yCenter = apData->dy[iz];
    *xSize = apData->xMax[iz];
    *ySize = apData->yMax[iz];
  }
  else {
    *xCenter = INTERPOLATE(apData->dx[iz], apData->dx[iz+1], apData->s[iz], apData->s[iz+1], z0);
    *yCenter = INTERPOLATE(apData->dy[iz], apData->dy[iz+1], apData->s[iz], apData->s[iz+1], z0);
    *xSize = INTERPOLATE(apData->xMax[iz], apData->xMax[iz+1], apData->s[iz], apData->s[iz+1], z0);
    *ySize = INTERPOLATE(apData->yMax[iz], apData->yMax[iz+1], apData->s[iz], apData->s[iz+1], z0);
  }
  return 1;
}


long imposeApertureData(
                        double **initial, long np, double **accepted,
                        double z, double Po, APERTURE_DATA *apData)
{
  long ip, itop, lost;
  double *ini;
  double xSize, ySize;
  double xCenter, yCenter;
  double dx, dy;

#ifdef HAVE_GPU
   if(getElementOnGpu()){
      startGpuTimer();
      ip = gpu_imposeApertureData(np, accepted, z, Po, apData);
#ifdef GPU_VERIFY     
      startCpuTimer();
      imposeApertureData(initial, np, accepted, z, Po, apData);
      compareGpuCpu(np, "imposeApertureData");
#endif /* GPU_VERIFY */
      return ip;
    }
#endif /* HAVE_GPU */

#if DEBUG_APERTURE
  static FILE *fp = NULL;
  if (!fp) {
    TRACKING_CONTEXT tcontext;
    char s[1000];
    getTrackingContext(&tcontext);
    sprintf(s, "%s.aplos", tcontext.rootname);
    fflush(stdout);
    if (!(fp = fopen(s, "w")))
      bombElegant("unable to open debug file for aperture losses", NULL);
    fprintf(fp, "SDDS1\n");
    fprintf(fp, "&column name=z type=double units=m &end\n");
    fprintf(fp, "&column name=x type=double units=m &end\n");
    fprintf(fp, "&column name=y type=double units=m &end\n");
    fprintf(fp, "&data mode=ascii no_row_counts=1 &end\n");
  }  
#endif

  if (!interpolateApertureData(z, apData, 
                               &xCenter, &yCenter, &xSize, &ySize))
    return np;
  
  itop = np-1;
  for (ip=0; ip<np; ip++) {
    ini = initial[ip];
    dx = ini[0] - xCenter;
    dy = ini[2] - yCenter;
    lost = 0;
    if ((xSize && fabs(dx) > xSize) ||
        (ySize && fabs(dy) > ySize))
      lost = 1;
    if (lost) {
#if DEBUG_APERTURE
      fprintf(fp, "%e %e %e\n", z, initial[ip][0], initial[ip][2]);
#endif
      swapParticles(initial[ip], initial[itop]);
      if (accepted)
        swapParticles(accepted[ip], accepted[itop]);
      initial[itop][4] = z; /* record position of particle loss */
      initial[itop][5] = Po*(1+initial[itop][5]);
      --itop;
      --ip;
      --np;
    }
  }
  
#if DEBUG_APERTURE
  fflush(fp);
#endif

  return(np);
}

long track_through_speedbump(double **initial, SPEEDBUMP *speedbump, long np, double **accepted, double z,
                             double Po
                             )
{
  double *ini, radius, xiHit, yiHit;
  long iplane, ip, idir, hit, phit;
  long dsign[2] = {1, -1};
  long dflag[2] = {0, 0};

  if (!((speedbump->direction = interpretScraperDirection(speedbump->insertFrom, -1))&(DIRECTION_X|DIRECTION_Y))) {
    exactDrift(initial, np, speedbump->length);
    return np;
  }

  iplane = 0;
  if (speedbump->direction&DIRECTION_X) {
    /* For now, x and +x are the same */
    iplane = 0;
    dflag[0] = speedbump->direction&DIRECTION_PLUS_X ? 1 : 0;
    dflag[1] = speedbump->direction&DIRECTION_MINUS_X ? 1 : 0;
  }
  else if (speedbump->direction&DIRECTION_Y) {
    iplane = 2;
    dflag[0] = speedbump->direction&DIRECTION_PLUS_Y ? 1 : 0;
    dflag[1] = speedbump->direction&DIRECTION_MINUS_Y ? 1 : 0;
  } else
    return np;

  if (speedbump->chord<0 || speedbump->height<0)
    bombElegant("SPEEDBUMP element needs non-negative CHORD and HEIGHT", NULL);

  if ((speedbump->length/2-speedbump->chord/2+speedbump->dzCenter)<0)
    bombElegant("SPEEDBUMP element is displaced too far upstream---the element, reduce the chord, or reduce |dzCenter|", NULL);
  if ((speedbump->length/2+speedbump->chord/2-speedbump->dzCenter)<0)
    bombElegant("SPEEDBUMP element is displaced too far downstream---lengthen the element, reduce the chord, or reduce |dzCenter|", NULL);

  if (speedbump->chord==0 || speedbump->height==0)
    radius = 0; /* doesn't matter */
  else
    radius = (sqr(speedbump->chord) + 4*sqr(speedbump->height))/(8*speedbump->height);

  /* idealized speedbump that just absorbs particles */
  for (ip=0; ip<np; ip++) {
    double xi=0, yi=0;
    ini = initial[ip];
    xiHit = yiHit = DBL_MAX; /* no hit */
    phit = 0;
    for (idir=0; idir<2; idir++) {
      /* Variables for circle-line intersection.
       * See Circle-Line Intersection article on Wolfram MathWorld 
       * x is the longitudinal coordinate, y is the transverse coordinate 
       */
      double x1, x2, y1, y2, dx, dy, dr, D, disc; 
      double dxPlane, offset;
      if (!dflag[idir]) continue;
      offset = iplane==0 ? speedbump->dx : speedbump->dy;
      hit = 0;
      dxPlane = 2*speedbump->length; /* by default, we miss */
      /* 1. Check for intersection of the ray with the substrate plane upstream of the bump */
      if (dsign[idir]*ini[iplane]>=(speedbump->position+speedbump->height+offset*dsign[idir])) {
        /* hit the leading edge */
        dxPlane = 0;
        xi = 0;
        yi = ini[iplane];
        hit = 1;
      } else if (ini[iplane+1]) {
        dxPlane = (dsign[idir]*(speedbump->position+speedbump->height+offset*dsign[idir])-ini[iplane])/ini[iplane+1];
        if (dxPlane>=0 && dxPlane<=(speedbump->length/2 + speedbump->dzCenter - speedbump->chord/2)) {
          hit = 1;
          xi = dxPlane;
          yi = (speedbump->position+speedbump->height)*dsign[idir] + offset;
        }
      }
      if (!hit && radius) {
        /* 2. Check for a hit on the bump */
        x1 = -(speedbump->length/2 + speedbump->dzCenter);
        dx = speedbump->length;
        x2 = x1 + dx;
        y1 = speedbump->position + radius + offset*dsign[idir] - dsign[idir]*ini[iplane];
        dy = -speedbump->length*dsign[idir]*ini[iplane+1];
        y2 = y1 + dy;
        D = x1*y2 - x2*y1;
        dr = sqrt(dx*dx + dy*dy);
        if ((disc = sqr(radius*dr) - sqr(D))>=0) {
          if (disc==0) {
            /* tangent hit */
            xi = D*dy/sqr(dr);
            yi = -D*dx/sqr(dr);
            hit = 1;
            /* put into accelerator coordinates */
            xi += speedbump->length/2 + speedbump->dzCenter;
            yi = (yi - (speedbump->position + radius + offset*dsign[idir]))/(-dsign[idir]);
          } else {
            /* secant line */
            double xi1, xi2, yi1, yi2;
            xi1 = (D*dy + SIGN(dy)*dx*sqrt(disc))/sqr(dr);
            yi1 = (-D*dx + fabs(dy)*sqrt(disc))/sqr(dr);
            xi2 = (D*dy - SIGN(dy)*dx*sqrt(disc))/sqr(dr);
            yi2 = (-D*dx - fabs(dy)*sqrt(disc))/sqr(dr);
            if (xi1<xi2) {
              xi = xi1;
              yi = yi1;
            } else {
              xi = xi2;
              yi = yi2;
            }
            hit = 1;
            xi += speedbump->length/2 + speedbump->dzCenter;
            yi = (yi - (speedbump->position + radius + offset*dsign[idir]))/(-dsign[idir]);
          }
        }
        if (!hit && dxPlane>=(speedbump->length/2+speedbump->dzCenter) && dxPlane<=speedbump->length) {
          /* we missed the bump but hit the plane downstream */
          xi = dxPlane;
          yi = (speedbump->position+speedbump->height)*dsign[idir] + offset;
          hit = 1;
        }
      }
      if (hit && xi<xiHit) {
        xiHit = xi;
        yiHit = yi;
        phit = 1;
      }
    }
    if (!phit) {
      /* Particle survives */
      ini[0] += speedbump->length*ini[1];
      ini[2] += speedbump->length*ini[3];
      ini[4] += speedbump->length*sqrt(1+sqr(ini[1])+sqr(ini[3]));
    } else {
      swapParticles(initial[ip], initial[np-1]);
      if (accepted)
        swapParticles(accepted[ip], accepted[np-1]);
      initial[np-1][iplane==0 ? 0 : 2] = yiHit;
      initial[np-1][iplane==0 ? 2 : 0] += xiHit*initial[np-1][iplane==0 ? 3 : 1];
      initial[np-1][4] = z+xiHit;
      initial[np-1][5] = Po*(1+initial[np-1][5]);
      --ip;
      --np;
    }
  }

  return(np);
}

int pointIsInsideContour(double x0, double y0, double *x, double *y, long n);

long trackThroughApContour(double **coord, APCONTOUR *apcontour, long np, double **accepted, double z,
                           double Po
                           )
{
  long ip, i_top;
  double z0, z1, zLost;
  short lost0, lost1, lost2;
  int lossCode;
  
  if (!apcontour->initialized) {
    SDDS_DATASET SDDSin;
    SDDSin.parallel_io = 0;
    if (apcontour->x) free(apcontour->x);
    if (apcontour->y) free(apcontour->y);
    apcontour->x = apcontour->y = NULL;
    apcontour->nPoints = 0;
    if (!apcontour->filename || !strlen(apcontour->filename))
      bombElegantVA("Error: No filename given for APCONTOUR\n", apcontour->filename);
    if (!apcontour->xColumn || !strlen(apcontour->xColumn))
      bombElegantVA("Error: No XCOLUMN given for APCONTOUR\n", apcontour->xColumn);
    if (!apcontour->yColumn || !strlen(apcontour->yColumn))
      bombElegantVA("Error: No YCOLUMN given for APCONTOUR\n", apcontour->yColumn);
    if (!SDDS_InitializeInputFromSearchPath(&SDDSin, apcontour->filename) ||
	SDDS_ReadPage(&SDDSin)!=1 ||
	(apcontour->nPoints = SDDS_RowCount(&SDDSin))<0 ||
	apcontour->nPoints<3)
      bombElegantVA("Error: APCONTOUR file %s is unreadable, or has insufficient data (<3 points)\n", apcontour->filename);
    if (SDDS_CheckColumn(&SDDSin, apcontour->xColumn, "m", SDDS_ANY_FLOATING_TYPE, stdout)!=SDDS_CHECK_OK) 
      bombElegantVA("Error: problem with x column (%s) for APCONTOUR file %s---check existence, units, and type\n", 
		    apcontour->xColumn, apcontour->filename);
    if (SDDS_CheckColumn(&SDDSin, apcontour->yColumn, "m", SDDS_ANY_FLOATING_TYPE, stdout)!=SDDS_CHECK_OK) 
      bombElegantVA("Error: problem with y column (%s) for APCONTOUR file %s---check existence, units, and type\n", 
		    apcontour->yColumn, apcontour->filename);
    if (!(apcontour->x = SDDS_GetColumnInDoubles(&SDDSin, apcontour->xColumn)) ||
	!(apcontour->y = SDDS_GetColumnInDoubles(&SDDSin, apcontour->yColumn)))
      bombElegantVA("Error: failed to get x or y data from APCONTOUR file %s\n", apcontour->filename);
    SDDS_Terminate(&SDDSin);
    if (apcontour->resolution<=0)
      bombElegant("Error: APCONTOUR has RESOLUTION<=0", NULL);
    if (apcontour->x[0]!=apcontour->x[apcontour->nPoints-1] ||
	apcontour->y[0]!=apcontour->y[apcontour->nPoints-1])
      bombElegantVA("Error: contour provided in file %s for APCONTOUR is not a closed shape\n", apcontour->filename);
    printf("Read aperture contour data from file %s\n", apcontour->filename);
    fflush(stdout);
    apcontour->initialized = 1;
  }

  /* misalignments */
  if (apcontour->dx || apcontour->dy || apcontour->dz)
    offsetBeamCoordinates(coord, np, apcontour->dx, apcontour->dy, apcontour->dz);
  if (apcontour->tilt)
    rotateBeamCoordinates(coord, np, apcontour->tilt);

  lossCode = 0;
  if (apcontour->invert)
    lossCode = 1;
  i_top = np - 1;
  for (ip=0; ip<=i_top; ip++) {
    z0 = zLost = 0;
    z1 = apcontour->length;
    lost0 = lost1 = 0;
    if (pointIsInsideContour(coord[ip][0]+coord[ip][1]*z0,
			     coord[ip][2]+coord[ip][3]*z0,
			     apcontour->x, apcontour->y, apcontour->nPoints)==lossCode) {
      lost0 = 1;
    } else if (apcontour->length>0) {
      if (pointIsInsideContour(coord[ip][0]+coord[ip][1]*z1,
			       coord[ip][2]+coord[ip][3]*z1,
			       apcontour->x, apcontour->y, apcontour->nPoints)==lossCode) {
	lost1 = 1;
	while ((z1-z0)>apcontour->resolution) {
	  zLost = (z0+z1)/2;
	  lost2 = pointIsInsideContour(coord[ip][0]+coord[ip][1]*zLost,
				       coord[ip][2]+coord[ip][3]*zLost,
				       apcontour->x, apcontour->y, apcontour->nPoints)==lossCode;
	  if (lost2==lost1)
	    z1 = zLost;
	  else if (lost2==lost0)
	    z0 = zLost;
	}
	zLost = (z0+z1)/2;
      }
    }
    if (lost0+lost1) {
      exactDrift(coord+ip, 1, zLost);
      coord[ip][4] = z + zLost;
      coord[ip][5] = Po*(1+coord[ip][5]);
      swapParticles(coord[ip], coord[i_top]);
      if (accepted)
        swapParticles(accepted[ip], accepted[i_top]);
      --i_top;
      --ip;
    } else
      exactDrift(coord+ip, 1, apcontour->length);
  }

  /* misalignments */
  if (apcontour->tilt)
    rotateBeamCoordinates(coord, np, -apcontour->tilt);
  if (apcontour->dx || apcontour->dy || apcontour->dz)
    offsetBeamCoordinates(coord, np, -apcontour->dx, -apcontour->dy, -apcontour->dz);

  return i_top+1;
}

long trackThroughTaperApCirc(double **initial, TAPERAPC *taperApC, long np, double **accepted, double z,
                             double Po)
{
  long ip, itop, isLost;
  double *coord, determinant, rho, rho2, r, rStart, rStart2, r2Limit, dz, dx, dy;
  double x, y, xp, yp;

  rho = rho2 = rStart = rStart2 = r2Limit = 0;
  if (taperApC->length>0) {
    rho = (taperApC->r[taperApC->e2Index]-taperApC->r[taperApC->e1Index])/taperApC->length;
    rho2 = sqr(rho);
    rStart = taperApC->r[taperApC->e1Index];
    rStart2 = sqr(taperApC->r[taperApC->e1Index]);
  } else {
    rho = 0;
    r = MIN(taperApC->r[taperApC->e1Index], taperApC->r[taperApC->e2Index]);
    r2Limit = sqr(r);
  }
  dx = taperApC->dx;
  dy = taperApC->dy;

  itop = np-1;
  for (ip=0; ip<=itop; ip++) {
    coord = initial[ip];
    isLost = 0;
    x = coord[0] - dx;
    y = coord[2] - dy;
    xp = coord[1];
    yp = coord[3];
    dz = 0;
    if ((sqr(x)+sqr(y))>=rStart2) 
      isLost = 1;
    else if (rho) {
      double a, b, c;
      a = sqr(xp) + sqr(yp) - rho2;
      b = 2*x*xp + 2*y*yp - 2*rStart*rho;
      c = sqr(x) + sqr(y) - rStart2;
      if ((determinant = sqr(b)-4*a*c)>=0 && (dz = (-b+sqrt(determinant))/(2*a))>=0 && dz<=taperApC->length)
        isLost = 1;
    } else if ((sqr(x)+sqr(y))>=r2Limit)
      isLost = 1;
    if (isLost) {
      coord[0] += dz*xp;
      coord[2] += dz*yp;
      swapParticles(initial[ip], initial[itop]);
      if (accepted)
        swapParticles(accepted[ip], accepted[itop]);
      initial[itop][4] = z+dz;
      initial[itop][5] = Po*(1+initial[itop][5]);
      --itop;
      --ip;
    } else
      exactDrift(initial+ip, 1, taperApC->length);
  }
  return itop+1;
}

int insideTaperedEllipse(double x0, double xp, double y0, double yp, double z,
			 double a0, double dadz, double b0, double dbdz,
			 long xExponent, long yExponent)
{
  double x, y, r2, theta, ct, st, a, b;
  double xe, ye, re2;
  
  x = x0 + xp*z;
  y = y0 + yp*z;
  r2 = sqr(x) + sqr(y);
  theta = atan2(y, x);
  a = dadz*z + a0;
  b = dbdz*z + b0;
  ct = cos(theta);
  st = sin(theta);
  xe = a*pow(fabs(ct), 2/xExponent)*SIGN(ct);
  ye = b*pow(fabs(st), 2/yExponent)*SIGN(st);
  re2 = sqr(xe)+sqr(ye);
  if (re2<=r2)
    return 0;
  return 1;
}

long trackThroughTaperApElliptical(double **initial, TAPERAPE *taperApE, long np, double **accepted, double zStartElem,
                                   double Po)
{
  long ip, itop, isLost0, isLost1, isLost2;
  double z0, z1, x0, y0, xp, yp, zLost, dadz, dbdz;
  double *coord;

  if (taperApE->length<0)
    bombElegant("TAPERAPE has negative length, which is not allowed", NULL);
  if (taperApE->length>0) {
    dadz = (taperApE->a[taperApE->e2Index]-taperApE->a[taperApE->e1Index])/taperApE->length;
    dbdz = (taperApE->b[taperApE->e2Index]-taperApE->b[taperApE->e1Index])/taperApE->length;
  } else {
    dadz = dbdz = 0;
  }

  if (taperApE->dx || taperApE->dy)
    offsetBeamCoordinates(initial, np, taperApE->dx, taperApE->dy, 0);
  if (taperApE->tilt)
    rotateBeamCoordinates(initial, np, taperApE->tilt);

  itop = np-1;
  for (ip=0; ip<=itop; ip++) {
    coord = initial[ip];
    x0 = coord[0];
    y0 = coord[2];
    xp = coord[1];
    yp = coord[3];
    z0 = zLost = 0;
    z1 = taperApE->length;
    isLost0 = isLost1 = 0;
    if (!insideTaperedEllipse(x0, xp, y0, yp, z0,
			      taperApE->a[taperApE->e1Index], dadz, 
			      taperApE->b[taperApE->e1Index], dbdz, 
			      taperApE->xExponent, taperApE->yExponent)) {
      isLost0 = 1;
      zLost = z0;
    } else if (taperApE->length>0) {
      if (!insideTaperedEllipse(x0, xp, y0, yp, z1,
				taperApE->a[taperApE->e1Index], dadz, 
				taperApE->b[taperApE->e1Index],	dbdz,
				taperApE->xExponent, taperApE->yExponent)) {
	isLost1 = 1;
	while ((z1-z0)>taperApE->resolution) {
	  zLost = (z0+z1)/2;
	  isLost2 = !insideTaperedEllipse(x0, xp, y0, yp, zLost,
					 taperApE->a[taperApE->e1Index], dadz, 
					 taperApE->b[taperApE->e1Index], dbdz,
					 taperApE->xExponent, taperApE->yExponent);
	  if (isLost2==isLost1)
	    z1 = zLost;
	  else if (isLost2==isLost0)
	    z0 = zLost;
	}
	zLost = (z0+z1)/2;
      }
    }

    if (isLost0+isLost1) {
      coord[0] += zLost*xp;
      coord[2] += zLost*yp;
      swapParticles(initial[ip], initial[itop]);
      if (accepted)
        swapParticles(accepted[ip], accepted[itop]);
      initial[itop][4] = zStartElem+zLost;
      initial[itop][5] = Po*(1+initial[itop][5]);
      --itop;
      --ip;
    } else
      exactDrift(initial+ip, 1, taperApE->length);
  }

  if (taperApE->tilt)
    rotateBeamCoordinates(initial, np, -taperApE->tilt);
  if (taperApE->dx || taperApE->dy)
    offsetBeamCoordinates(initial, np, -taperApE->dx, -taperApE->dy, 0);

  return itop+1;
}

long trackThroughTaperApRectangular(double **initial, TAPERAPR *taperApR, long np, double **accepted, double zStartElem,
                                    double Po)
{
  long ip, itop, isLost;
  double *coord;
  double x0, y0, xp, yp, zLost, zLost0;

  if (taperApR->length<0)
    bombElegant("TAPERAPR has negative length, which is not allowed", NULL);

  if (taperApR->dx || taperApR->dy)
    offsetBeamCoordinates(initial, np, taperApR->dx, taperApR->dy, 0);
  if (taperApR->tilt)
    rotateBeamCoordinates(initial, np, taperApR->tilt);

  itop = np-1;
  for (ip=0; ip<=itop; ip++) {
    coord = initial[ip];
    x0 = coord[0];
    y0 = coord[2];
    xp = coord[1];
    yp = coord[3];

    isLost = 0;
    zLost = zLost0 = DBL_MAX;

    if (taperApR->xmax[taperApR->e1Index]>0 && taperApR->xmax[taperApR->e2Index]>0) {
      /* x>0 plane */
      if (x0>taperApR->xmax[taperApR->e1Index]) {
        isLost = 1;
        zLost = 0;
      } else {
        /* find intersection */
        zLost0 = (x0 - taperApR->xmax[taperApR->e1Index])*taperApR->length/
          (taperApR->xmax[taperApR->e2Index] - taperApR->xmax[taperApR->e1Index] - xp*taperApR->length);
        if (zLost0>=0 && zLost0<=taperApR->length) {
          isLost = 1;
          zLost = zLost0;
        }
      }
      if (!isLost || zLost>0) {
        /* x<0 plane */
        if (x0 < -taperApR->xmax[taperApR->e1Index]) {
          isLost = 1;
          zLost = 0;
        } else {
          /* find intersection */
          zLost0 = (x0 + taperApR->xmax[taperApR->e1Index])*taperApR->length/
            (-taperApR->xmax[taperApR->e2Index] + taperApR->xmax[taperApR->e1Index] - xp*taperApR->length);
          if (zLost0>=0 && zLost0<=taperApR->length) {
            isLost = 1;
            if (zLost0<zLost)
              zLost = zLost0;
          }
        }
      }
    }
    if (taperApR->ymax[taperApR->e1Index]>0 && taperApR->ymax[taperApR->e2Index]>0) {
      if (!isLost || zLost>0) {
        /* y>0 plane */
        if (y0>taperApR->ymax[taperApR->e1Index]) {
          isLost = 1;
          zLost = 0;
        } else {
          /* find intersection */
          zLost0 = (y0 - taperApR->ymax[taperApR->e1Index])*taperApR->length/
            (taperApR->ymax[taperApR->e2Index] - taperApR->ymax[taperApR->e1Index] - yp*taperApR->length);
          if (zLost0>=0 && zLost0<=taperApR->length) {
            isLost = 1;
            if (zLost0<zLost)
              zLost = zLost0;
          }
        }
        if (!isLost || zLost>0) {
          /* y<0 plane */
          if (y0 < -taperApR->ymax[taperApR->e1Index]) {
            isLost = 1;
            zLost = 0;
          } else {
            /* find intersection */
            zLost0 = (y0 + taperApR->ymax[taperApR->e1Index])*taperApR->length/
              (-taperApR->ymax[taperApR->e2Index] + taperApR->ymax[taperApR->e1Index] - yp*taperApR->length);
            if (zLost0>=0 && zLost0<=taperApR->length) {
              isLost = 1;
              if (zLost0<zLost)
                zLost = zLost0;
            }
          }
        }
      }
    }

    if (isLost) {
      coord[0] += zLost*xp;
      coord[2] += zLost*yp;
      swapParticles(initial[ip], initial[itop]);
      if (accepted)
        swapParticles(accepted[ip], accepted[itop]);
      initial[itop][4] = zStartElem+zLost;
      initial[itop][5] = Po*(1+initial[itop][5]);
      --itop;
      --ip;
    } else
      exactDrift(initial+ip, 1, taperApR->length);
  }

  if (taperApR->tilt)
    rotateBeamCoordinates(initial, np, -taperApR->tilt);
  if (taperApR->dx || taperApR->dy)
    offsetBeamCoordinates(initial, np, -taperApR->dx, -taperApR->dy, 0);

  return itop+1;
}

