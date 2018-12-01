/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* routine: polynomial_kicks()
 * purpose: apply kicks due to a polynomial dependence on x or y
 * 
 * Michael Borland, 1992 
 */
#include "mdb.h"
#include "track.h"

long polynomial_kicks(
    double **particle,  /* initial/final phase-space coordinates */
    long n_part,        /* number of particles */
    KPOLY *kpoly,       /* kick-polynomial structure */
    double p_error,     /* p_nominal/p_central */
    double Po,
    double **accepted,
    double z_start
    )
{
    double KL0;         /* on-momentum integrated strength at x=1 or y=1 */
    double KL;
    double dx, dy, dz;  /* offsets of the multipole center */
    long order;         /* order (n) */
    long i_part, i_top, yplane;
    double *coord;
    double cos_tilt, sin_tilt;
    double x, xp, y, yp;

    log_entry("polynomial_kicks");

    if (!particle)
        bombElegant("particle array is null (polynomial_kicks)", NULL);

    if (!kpoly)
        bombElegant("null KPOLY pointer (polynomial_kicks)", NULL);

    if ((order=kpoly->order)<0)
        bombElegant("order < 0 for KPOLY element (polynomial_kicks)", NULL);

    if (kpoly->plane && (kpoly->plane[0]=='y' || kpoly->plane[0]=='Y'))
        yplane = 1;
    else {
        if (kpoly->plane && !(kpoly->plane[0]=='x' || kpoly->plane[0]=='X')) {
            fputs("warning: KPOLY plane not recognized--x plane assumed.", stdout);
            cp_str(&kpoly->plane, "x");
            }
        yplane = 0;
        }

    KL0 = kpoly->coefficient*kpoly->factor;

    cos_tilt = cos(kpoly->tilt);
    sin_tilt = sin(kpoly->tilt);
    dx = kpoly->dx;
    dy = kpoly->dy;
    dz = kpoly->dz;

    i_top = n_part-1;
    for (i_part=0; i_part<=i_top; i_part++) {
        if (!(coord = particle[i_part])) {
            printf("null coordinate pointer for particle %ld (polynomial_kicks)", i_part);
            fflush(stdout);
            abort();
            }
        if (accepted && !accepted[i_part]) {
            printf("null accepted coordinates pointer for particle %ld (polynomial_kicks)", i_part);
            fflush(stdout);
            abort();
            }

        /* adjust strength for momentum offset--good to all orders */
        KL = KL0/(1+coord[5]);

        /* calculate coordinates in rotated and offset frame */
        coord[4] += dz*sqrt(1 + sqr(coord[1]) + sqr(coord[3]));
        coord[0] += -dx + dz*coord[1];
        coord[2] += -dy + dz*coord[3];

        x  =   cos_tilt*coord[0] + sin_tilt*coord[2];
        y  = - sin_tilt*coord[0] + cos_tilt*coord[2];
        xp =   cos_tilt*coord[1] + sin_tilt*coord[3];
        yp = - sin_tilt*coord[1] + cos_tilt*coord[3];

#if defined(IEEE_MATH)
        if (isnan(x) || isnan(xp) || isnan(y) || isnan(yp)) {
            swapParticles(particle[i_part], particle[i_top]);
            if (accepted)
                swapParticles(accepted[i_part], accepted[i_top]);
            particle[i_top][4] = z_start;
            particle[i_top][5] = Po*(1+particle[i_top][5]);
            i_top--;
            i_part--;
            continue;
            }
#endif
        if (FABS(x)>COORD_LIMIT || FABS(y)>COORD_LIMIT ||
            FABS(xp)>SLOPE_LIMIT || FABS(yp)>SLOPE_LIMIT) {
            swapParticles(particle[i_part], particle[i_top]);
            if (accepted)
                swapParticles(accepted[i_part], accepted[i_top]);
            particle[i_top][4] = z_start;
            particle[i_top][5] = Po*(1+particle[i_top][5]);
            i_top--;
            i_part--;
            continue;
            }

        if (yplane)
            yp += KL*ipow(y, order);
        else
            xp += KL*ipow(x, order);
                
        /* undo the rotation and store in place of initial coordinates */
        /* don't need to change coord[0] or coord[2] since x and y are unchanged */
        coord[1] = cos_tilt*xp - sin_tilt*yp;
        coord[3] = sin_tilt*xp + cos_tilt*yp;

        /* remove the coordinate offsets */
        coord[0] += dx - coord[1]*dz;
        coord[2] += dy - coord[3]*dz;
        coord[4] -= dz*sqrt(1+ sqr(coord[1]) + sqr(coord[3]));
        }
    log_exit("polynomial_kicks");
    return(i_top+1);
    }

long polynomial_hamiltonian(
    double **particle,  /* initial/final phase-space coordinates */
    long n_part,        /* number of particles */
    HKPOLY *hkpoly,     /* kick-polynomial structure */
    double p_error,     /* p_nominal/p_central */
    double Po,
    double **accepted,
    double z_start
    )
{
  double dx, dy, dz;  /* offsets of the multipole center */
  long i_part, i_top, ic, ix, iy;
  double *coord, kick;
  double cos_tilt, sin_tilt;
  double x, xp, y, yp, qx, qy;
  double dl;
  long ik, nk;

  if (!particle)
    bombElegant("particle array is null (polynomial_hamiltonian)", NULL);
  
  if (!hkpoly)
    bombElegant("null HKPOLY pointer (polynomial_hamiltonian)", NULL);
  if (hkpoly->nKicks<=0)
    bombElegant("HKPOLY N_KICKS must be positive (polynomial_hamiltonian)", NULL);
  if (hkpoly->length<0)
    bombElegant("HKPOLY length (L) must be non-negative (polynomial_hamiltonian)", NULL);

  cos_tilt = cos(hkpoly->tilt);
  sin_tilt = sin(hkpoly->tilt);
  dx = hkpoly->dx;
  dy = hkpoly->dy;
  dz = hkpoly->dz;

  nk = hkpoly->nKicks;
  if ((dl = hkpoly->length/hkpoly->nKicks)==0)
    nk = 1;
  
  i_top = n_part-1;
  for (i_part=0; i_part<=i_top; i_part++) {
    if (!(coord = particle[i_part])) {
      printf("null coordinate pointer for particle %ld (polynomial_hamiltonian)", i_part);
      fflush(stdout);
      abort();
    }
    if (accepted && !accepted[i_part]) {
      printf("null accepted coordinates pointer for particle %ld (polynomial_hamiltonian)", i_part);
      fflush(stdout);
      abort();
    }
    
    /* calculate coordinates in rotated and offset frame */
    coord[4] += dz*sqrt(1 + sqr(coord[1]) + sqr(coord[3]));
    coord[0] += -dx + dz*coord[1];
    coord[2] += -dy + dz*coord[3];
    
    x  =   cos_tilt*coord[0] + sin_tilt*coord[2];
    y  = - sin_tilt*coord[0] + cos_tilt*coord[2];
    xp =   cos_tilt*coord[1] + sin_tilt*coord[3];
    yp = - sin_tilt*coord[1] + cos_tilt*coord[3];
    
    if (FABS(x)>COORD_LIMIT || FABS(y)>COORD_LIMIT ||
        FABS(xp)>SLOPE_LIMIT || FABS(yp)>SLOPE_LIMIT) {
      swapParticles(particle[i_part], particle[i_top]);
      if (accepted)
        swapParticles(accepted[i_part], accepted[i_top]);
      particle[i_top][4] = z_start;
      particle[i_top][5] = Po*(1+particle[i_top][5]);
      i_top--;
      i_part--;
      continue;
    }

#if defined(IEEE_MATH)
    if (isnan(x) || isnan(xp) || isnan(y) || isnan(yp)) {
      swapParticles(particle[i_part], particle[i_top]);
      if (accepted)
        swapParticles(accepted[i_part], accepted[i_top]);
      particle[i_top][4] = z_start;
      particle[i_top][5] = Po*(1+particle[i_top][5]);
      i_top--;
      i_part--;
      continue;
    }
#endif

    convertSlopesToMomenta(&qx, &qy, xp, yp, coord[5]);
    for (ik=0; ik<nk; ik++) {
      if (dl) {
        x += xp*dl/2;
        y += yp*dl/2;
      }
      for (ic=0; ic<24; ic++) {
        if (hkpoly->coefficient[ic]) {
          ix = (ic+1)%5;
          iy = (ic+1)/5;
          if (ix || iy) {
            if (ix)
              qx -= hkpoly->factor*hkpoly->coefficient[ic]*ix/(1+coord[5])*ipow(x, ix-1)*ipow(y, iy  )/hkpoly->nKicks;
            if (iy)
              qy -= hkpoly->factor*hkpoly->coefficient[ic]*iy/(1+coord[5])*ipow(x, ix  )*ipow(y, iy-1)/hkpoly->nKicks;
            convertMomentaToSlopes(&xp, &yp, qx, qy, coord[5]);
          }
        }
      }
      if (dl) {
        x += xp*dl/2;
        y += yp*dl/2;
      }
    }
         
    /* undo the rotation and store in place of initial coordinates */
    /* don't need to change coord[0] or coord[2] since x and y are unchanged */
    coord[1] = cos_tilt*xp - sin_tilt*yp;
    coord[3] = sin_tilt*xp + cos_tilt*yp;

    /* remove the coordinate offsets */
    coord[0] += dx - coord[1]*dz;
    coord[2] += dy - coord[3]*dz;
    coord[4] -= dz*sqrt(1+ sqr(coord[1]) + sqr(coord[3]));
  }
  return(i_top+1);
}
