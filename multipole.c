/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* routine: multipole()
 * purpose: apply kicks due to an arbitrary multipole
 * 
 * Michael Borland, 1991. 
 */
#include "mdb.h"
#include "track.h"

double *expansion_coefficients(long n);
void apply_canonical_multipole_kicks(double *qx, double *qy,
                                   double *sum_Fx, double *sum_Fy,
                                   double x, double y,
                                   long order, double KnL, long skew);
int integrate_kick_multipole_ord2(double *coord, double cos_tilt, double sin_tilt,
                                  double dx, double dy, double xkick, double ykick,
                                  double Po, double rad_coef,
                                  long order, double KnL, long n_kicks, double drift,
                                  MULTIPOLE_DATA *multData, MULTIPOLE_DATA *steeringMultData);
int integrate_kick_multipole_ord4(double *coord, double cos_tilt, double sin_tilt,
                                  double dx, double dy, double xkick, double ykick,
                                  double Po, double rad_coef,
                                  long order, double KnL, long n_kicks, double drift,
                                  MULTIPOLE_DATA *multData, MULTIPOLE_DATA *steeringMultData);
void computeTotalErrorMultipoleFields(MULTIPOLE_DATA *totalMult,
                                      MULTIPOLE_DATA *systematicMult,
                                      MULTIPOLE_DATA *randomMult,
                                      MULTIPOLE_DATA *steeringMult,
                                      double KmL, long rootOrder);
void randomizeErrorMultipoleFields(MULTIPOLE_DATA *randomMult);

unsigned long multipoleKicksDone = 0;

#define ODD(j) ((j)%2)

void readErrorMultipoleData(MULTIPOLE_DATA *multData,
                               char *multFile, long steering)
{
  SDDS_DATASET SDDSin;
  char buffer[1024];
  if (!multFile || !strlen(multFile)) {
    multData->orders = 0;
    multData->initialized = 0;
    return;
  }
  if (multData->initialized)
    return;
  if (!SDDS_InitializeInputFromSearchPath(&SDDSin, multFile)) {
    fprintf(stdout, "Problem opening file %s\n", multFile);
    fflush(stdout);
    exit(1);
  }
  if (SDDS_CheckColumn(&SDDSin, "order", NULL, SDDS_ANY_INTEGER_TYPE, stdout)!=SDDS_CHECK_OK ||
      SDDS_CheckColumn(&SDDSin, "an", NULL, SDDS_ANY_FLOATING_TYPE, stdout)!=SDDS_CHECK_OK ||
      (!steering && SDDS_CheckColumn(&SDDSin, "bn", NULL, SDDS_ANY_FLOATING_TYPE, stdout)!=SDDS_CHECK_OK) ||
      SDDS_CheckParameter(&SDDSin, "referenceRadius", "m", SDDS_ANY_FLOATING_TYPE, stdout)!=SDDS_CHECK_OK) {
    fprintf(stdout, "Problems with data in multipole file %s\n", multFile);
    fflush(stdout);
    exit(1);
  }
  if (steering && SDDS_CheckColumn(&SDDSin, "bn", NULL, SDDS_ANY_FLOATING_TYPE, NULL)==SDDS_CHECK_OK) {
    fprintf(stdout, "Warning: Steering multipole file %s should not have bn.\n",
            multFile);
    fprintf(stdout, "Use an to specify multipole content for a horizontal steerer.\n");
    fprintf(stdout, "Multipole content for vertical steerer is deduced from this.\n");
    fflush(stdout);
  }
  if (SDDS_ReadPage(&SDDSin)!=1)  {
    sprintf(buffer, "Problem reading multipole file %s\n", multFile);
    SDDS_SetError(buffer);
    SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors);
    exit(1);
  }
  if ((multData->orders = SDDS_RowCount(&SDDSin))<=0) {
    fprintf(stdout, "Warning: no data in multipole file %s\n", multFile);
    fflush(stdout);
    SDDS_Terminate(&SDDSin);
  }
  if (!SDDS_GetParameterAsDouble(&SDDSin, "referenceRadius", &multData->referenceRadius) ||
      !(multData->order=SDDS_GetColumnInLong(&SDDSin, "order")) ||
      !(multData->an=SDDS_GetColumnInDoubles(&SDDSin, "an")) || 
      (!steering && !(multData->bn=SDDS_GetColumnInDoubles(&SDDSin, "bn")))) {
    sprintf(buffer, "Unable to read multipole data for file %s\n", multFile);
    SDDS_SetError(buffer);
    SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors);
    exit(1);
  }    
  if (steering &&
      !(multData->bn=SDDS_Malloc(sizeof(*(multData->bn))*multData->orders))) {
    fprintf(stdout, "Memory allocation failure (readErrorMultipoleData)\n");
    exit(1);
  }
  if (SDDS_ReadPage(&SDDSin)==2) {
    fprintf(stdout, "Warning: multipole file %s has multiple pages, which are ignored\n",
            multFile);
    fflush(stdout);
  }
  SDDS_Terminate(&SDDSin);
  if (steering) {
    long i, j;
    /* check for disallowed multipoles */
    for (i=0; i<multData->orders; i++) {
      if (ODD(multData->order[i])) {
        fprintf(stdout, "Error: steering multipole file %s has disallowed odd orders.\n",
                multFile);
        exit(1);
      }
    }
    /* normalize to n=0 if present */
    /* find i such that order[i] is 0 (dipole) */
    for (i=0; i<multData->orders; i++) {
      if (multData->order[i]==0)
        break;
    }
    if (multData->orders>1 && i!=multData->orders) {
      /* dipole present */
      if (!multData->an[i] || !multData->bn[i]) {
        fprintf(stdout, "Steering multipole data in %s is invalid: an or bn is zero for order=0\n",
                multFile);
        exit(1);
      }
      /* normalize to dipole for normal and skew separately */
      for (j=0; j<multData->orders; j++)
        if (j!=i)
          multData->an[j] /= multData->an[i];
      /* remove the dipole data */
      for (j=i+1; j<multData->orders; j++) {
        multData->an[j-1] = multData->an[j];
        multData->order[j-1] = multData->order[j];
      }
      multData->orders -= 1;
      for (i=0; i<multData->orders; i++)
        multData->bn[i] = multData->an[i]*ipow(-1.0, multData->order[i]/2);
#ifdef DEBUG
      fprintf(stdout, "Steering multipole data: \n");
      for (i=0; i<multData->orders; i++)
        fprintf(stdout, "%ld: %e %e\n", multData->order[i],
                multData->an[i], multData->bn[i]);
#endif
    }
  }
  multData->initialized = 1;
}

void initialize_fmultipole(FMULT *multipole)
{
  SDDS_DATASET SDDSin;
  char buffer[1024];
  MULTIPOLE_DATA *multData;
  
  multData = &(multipole->multData);
  if (multData->initialized)
    return;
  if (!multipole->filename)
    bomb("FMULT element doesn't have filename", NULL);
  if (!SDDS_InitializeInputFromSearchPath(&SDDSin, multipole->filename)) {
    sprintf(buffer, "Problem opening file %s (FMULT)\n", multipole->filename);
    SDDS_SetError(buffer);
    SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors);
    exit(1);
  }
  if (SDDS_CheckColumn(&SDDSin, "order", NULL, SDDS_ANY_INTEGER_TYPE, stdout)!=SDDS_CHECK_OK ||
      SDDS_CheckColumn(&SDDSin, "KnL", NULL, SDDS_ANY_FLOATING_TYPE, stdout)!=SDDS_CHECK_OK ||
      SDDS_CheckColumn(&SDDSin, "JnL", NULL, SDDS_ANY_FLOATING_TYPE, stdout)!=SDDS_CHECK_OK)
    bomb("problems with data in FMULT input file", NULL);
  if (SDDS_ReadPage(&SDDSin)!=1)  {
    sprintf(buffer, "Problem reading FMULT file %s\n", multipole->filename);
    SDDS_SetError(buffer);
    SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors);
    exit(1);
  }
  if ((multData->orders = SDDS_RowCount(&SDDSin))<=0) {
    fprintf(stdout, "Warning: no data in FMULT file %s\n", multipole->filename);
    fflush(stdout);
    SDDS_Terminate(&SDDSin);
    return;
  }
  multData->JnL = NULL;
  if (!(multData->order=SDDS_GetColumnInLong(&SDDSin, "order")) ||
      !(multData->KnL=SDDS_GetColumnInDoubles(&SDDSin, "KnL")) ||
      !(multData->JnL=SDDS_GetColumnInDoubles(&SDDSin, "JnL"))) {
    sprintf(buffer, "Unable to read data for FMULT file %s\n", multipole->filename);
    SDDS_SetError(buffer);
    SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors);
    exit(1);
  }    
  if (SDDS_ReadPage(&SDDSin)==2) {
    fprintf(stdout, "Warning: FMULT file %s has multiple pages, which are ignored\n",
            multipole->filename);
    fflush(stdout);
  }
  SDDS_Terminate(&SDDSin);
  multData->initialized = 1;
}


long fmultipole_tracking(
                         double **particle,  /* initial/final phase-space coordinates */
                         long n_part,        /* number of particles */
                         FMULT *multipole,   /* multipole structure */
                         double p_error,     /* p_nominal/p_central */
                         double Po,
                         double **accepted,
                         double z_start
                         )
{
  double dx, dy, dz;  /* offsets of the multipole center */
  long n_kicks;       /* number of kicks to split multipole into */
  long i_part, i_top, is_lost=0, i_order;
  double *coord;
  double drift, cos_tilt, sin_tilt;
  double x=0.0, xp=0.0, y=0.0, yp=0.0;
  double rad_coef;
  MULTIPOLE_DATA multData;
  
  if (!particle)
    bomb("particle array is null (fmultipole_tracking)", NULL);

  if (!multipole)
    bomb("null MULT pointer (fmultipole_tracking)", NULL);
  
  if (!multipole->multData.initialized)
    initialize_fmultipole(multipole);
  
  if ((n_kicks=multipole->n_kicks)<=0)
    bomb("n_kicks<=0 in multipole()", NULL);

  drift = multipole->length;

  cos_tilt = cos(multipole->tilt);
  sin_tilt = sin(multipole->tilt);
  dx = multipole->dx;
  dy = multipole->dy;
  dz = multipole->dz;
  if (multipole->synch_rad)
    rad_coef = sqr(e_mks)*pow3(Po)/(6*PI*epsilon_o*sqr(c_mks)*me_mks);
  else
    rad_coef = 0;

  multData = multipole->multData;
  for (i_order=0; i_order<multData.orders; i_order++) {
    multData.KnL[i_order] *= (1+multipole->fse);
    if (multData.JnL)
      multData.JnL[i_order] *= (1+multipole->fse);
  }
  
  i_top = n_part-1;
  multipoleKicksDone += (i_top+1)*multData.orders*n_kicks*4;
  for (i_part=0; i_part<=i_top; i_part++) {
    if (!(coord = particle[i_part])) {
      fprintf(stdout, "null coordinate pointer for particle %ld (fmultipole_tracking)", i_part);
      fflush(stdout);
      abort();
    }
    if (accepted && !accepted[i_part]) {
      fprintf(stdout, "null accepted coordinates pointer for particle %ld (fmultipole_tracking)", i_part);
      fflush(stdout);
      abort();
    }

    coord[4] += dz*sqrt(1 + sqr(coord[1]) + sqr(coord[3]));
    coord[0]  = coord[0] + dz*coord[1];
    coord[2]  = coord[2] + dz*coord[3];

    if (!integrate_kick_multipole_ord4(coord, cos_tilt, sin_tilt, dx, dy, 0.0, 0.0, Po, rad_coef,
                                       1, 0.0, n_kicks, drift, &multData, NULL)) {
      is_lost = 1;
      break;
    }
    
    if (!is_lost) {
      coord[4] -= dz*sqrt(1 + sqr(coord[1]) + sqr(coord[3]));
      x = coord[0]  = coord[0] - dz*coord[1];
      y = coord[2]  = coord[2] - dz*coord[3];
      xp = coord[1];
      yp = coord[3];

#if defined(IEEE_MATH)
      if (isnan(x) || isnan(xp) || isnan(y) || isnan(yp)) {
        SWAP_PTR(particle[i_part], particle[i_top]);
        if (accepted)
          SWAP_PTR(accepted[i_part], accepted[i_top]);
        particle[i_top][4] = z_start;
        particle[i_top][5] = Po*(1+particle[i_top][5]);
        i_top--;
        i_part--;
        continue;
      }
#endif
    }
    if (is_lost || FABS(x)>COORD_LIMIT || FABS(y)>COORD_LIMIT ||
        FABS(xp)>SLOPE_LIMIT || FABS(yp)>SLOPE_LIMIT) {
      SWAP_PTR(particle[i_part], particle[i_top]);
      if (accepted)
        SWAP_PTR(accepted[i_part], accepted[i_top]);
      particle[i_top][4] = z_start;
      particle[i_top][5] = Po*(1+particle[i_top][5]);
      i_top--;
      i_part--;
      continue;
    }
  }
  log_exit("fmultipole_tracking");
  return(i_top+1);
}


long multipole_tracking(
    double **particle,  /* initial/final phase-space coordinates */
    long n_part,        /* number of particles */
    MULT *multipole,    /* multipole structure */
    double p_error,     /* p_nominal/p_central */
    double Po,
    double **accepted,
    double z_start
    )
{
    double KnL;         /* integrated strength = L/(B.rho)*(Dx^n(By))_o for central momentum */
    double dx, dy, dz;  /* offsets of the multipole center */
    long order;         /* order (n) */
    long n_kicks;       /* number of kicks to split multipole into */
    long i_part, i_kick, i, i_top, is_lost;
    double sum_Fx, sum_Fy, xypow, denom, qx, qy;
    double *coord;
    double drift, cos_tilt, sin_tilt;
    double *coef;
    double x, xp, y, yp, s, dp;
    double ratio, rad_coef;
    double beta0, beta1, p;

    log_entry("multipole_tracking");

    if (!particle)
        bomb("particle array is null (multipole_tracking)", NULL);

    if (!multipole)
        bomb("null MULT pointer (multipole_tracking)", NULL);

    if ((n_kicks=multipole->n_kicks)<=0)
        bomb("n_kicks<=0 in multipole()", NULL);

    if ((order=multipole->order)<0)
        bomb("order < 0 in multipole()", NULL);

    if (!(coef = expansion_coefficients(order)))
        bomb("expansion_coefficients() returned null pointer (multipole_tracking)", NULL);

    drift = multipole->length/n_kicks/2;
    if (multipole->bore)
        /* KnL = d^nB/dx^n * L/(B.rho) = n! B(a)/a^n * L/(B.rho) */
        KnL = dfactorial(multipole->order)*multipole->BTipL/ipow(multipole->bore, multipole->order)*
              (e_mks/(me_mks*c_mks*Po))*multipole->factor;
    else
      KnL = multipole->KnL*multipole->factor/n_kicks;

    cos_tilt = cos(multipole->tilt);
    sin_tilt = sin(multipole->tilt);
    dx = multipole->dx;
    dy = multipole->dy;
    dz = multipole->dz;
    if (multipole->synch_rad)
        rad_coef = sqr(e_mks)*pow3(Po)/(6*PI*epsilon_o*sqr(c_mks)*me_mks);
    else
        rad_coef = 0;

    i_top = n_part-1;
    multipoleKicksDone += (i_top+1)*n_kicks*4;
    for (i_part=0; i_part<=i_top; i_part++) {
        if (!(coord = particle[i_part])) {
            fprintf(stdout, "null coordinate pointer for particle %ld (multipole_tracking)", i_part);
            fflush(stdout);
            abort();
            }
        if (accepted && !accepted[i_part]) {
            fprintf(stdout, "null accepted coordinates pointer for particle %ld (multipole_tracking)", i_part);
            fflush(stdout);
            abort();
            }
        if (KnL==0) {
            coord[4] += multipole->length*sqrt(1+sqr(coord[1])+sqr(coord[3]));
            coord[0] += multipole->length*coord[1];
            coord[2] += multipole->length*coord[3];
            continue;
            }

        /* calculate coordinates in rotated and offset frame */
        coord[4] += dz*sqrt(1 + sqr(coord[1]) + sqr(coord[3]));
        coord[0]  = coord[0] - dx + dz*coord[1];
        coord[2]  = coord[2] - dy + dz*coord[3];

        x  =   cos_tilt*coord[0] + sin_tilt*coord[2];
        y  = - sin_tilt*coord[0] + cos_tilt*coord[2];
        xp =   cos_tilt*coord[1] + sin_tilt*coord[3];
        yp = - sin_tilt*coord[1] + cos_tilt*coord[3];
        s  = 0;
        dp = coord[5];
        p = Po*(1+dp);
        beta0 = p/sqrt(sqr(p)+1);

#if defined(IEEE_MATH)
        if (isnan(x) || isnan(xp) || isnan(y) || isnan(yp)) {
            SWAP_PTR(particle[i_part], particle[i_top]);
            if (accepted)
                SWAP_PTR(accepted[i_part], accepted[i_top]);
            particle[i_top][4] = z_start;
            particle[i_top][5] = Po*(1+particle[i_top][5]);
            i_top--;
            i_part--;
            continue;
            }
#endif
        if (FABS(x)>COORD_LIMIT || FABS(y)>COORD_LIMIT ||
            FABS(xp)>SLOPE_LIMIT || FABS(yp)>SLOPE_LIMIT) {
            SWAP_PTR(particle[i_part], particle[i_top]);
            if (accepted)
                SWAP_PTR(accepted[i_part], accepted[i_top]);
            particle[i_top][4] = z_start;
            particle[i_top][5] = Po*(1+particle[i_top][5]);
            i_top--;
            i_part--;
            continue;
            }

        /* calculate initial canonical momenta */
        qx = (1+dp)*xp/(denom=sqrt(1+sqr(xp)+sqr(yp)));
        qy = (1+dp)*yp/denom;
        is_lost = 0;
        for (i_kick=0; i_kick<n_kicks; i_kick++) {
            if (drift) {
                x += xp*drift*(i_kick?2:1);
                y += yp*drift*(i_kick?2:1);
                s += (i_kick?2:1)*drift*sqrt(1+sqr(xp)+sqr(yp));
                }
            if (x==0) {
                xypow = ipow(y, order);
                ratio = 0;
                i = order;
                }
            else {
                xypow = ipow(x, order);
                ratio = y/x;
                i = 0;
                }
            /* now sum up the terms for the multipole expansion */
            for (sum_Fx=sum_Fy=0; i<=order; i++) {
                if (ODD(i))
                    sum_Fx += coef[i]*xypow;
                else
                    sum_Fy += coef[i]*xypow;
                xypow *= ratio;
                }
            /* apply kicks canonically */
            qx -= KnL*sum_Fy;
            qy += KnL*sum_Fx;
            if ((denom=sqr(1+dp)-sqr(qx)-sqr(qy))<=0) {
                is_lost = 1;
                break;
                }
            xp = qx/(denom=sqrt(denom));
            yp = qy/denom;
            if (rad_coef && drift) {
                qx /= (1+dp);
                qy /= (1+dp);
                dp -= rad_coef*sqr(KnL*(1+dp))*(sqr(sum_Fy)+sqr(sum_Fx))*sqrt(1+sqr(xp)+sqr(yp))/(2*drift);
                qx *= (1+dp);
                qy *= (1+dp);
                }
            }
        if (drift && !is_lost) {
            /* go through final drift */
            x += xp*drift;
            y += yp*drift;
            s += drift*sqrt(1 + sqr(xp) + sqr(yp));
            }
        /* undo the rotation and store in place of initial coordinates */
        coord[0] = cos_tilt*x  - sin_tilt*y ;
        coord[2] = sin_tilt*x  + cos_tilt*y ;
        coord[1] = cos_tilt*xp - sin_tilt*yp;
        coord[3] = sin_tilt*xp + cos_tilt*yp;
        if (rad_coef) {
            p = Po*(1+dp);
            beta1 = p/sqrt(sqr(p)+1);
            coord[4] = beta1*(coord[4]/beta0 + 2*s/(beta0+beta1));
            }
        else 
            coord[4] += s;
        coord[5] = dp;

        /* remove the coordinate offsets */
        coord[0] += dx - coord[1]*dz;
        coord[2] += dy - coord[3]*dz;
        coord[4] -= dz*sqrt(1+ sqr(coord[1]) + sqr(coord[3]));

#if defined(IEEE_MATH)
        if (isnan(x) || isnan(xp) || isnan(y) || isnan(yp)) {
            SWAP_PTR(particle[i_part], particle[i_top]);
            if (accepted)
                SWAP_PTR(accepted[i_part], accepted[i_top]);
            particle[i_top][4] = z_start;
            particle[i_top][5] = Po*(1+particle[i_top][5]);
            i_top--;
            i_part--;
            continue;
            }
#endif
        if (FABS(x)>COORD_LIMIT || FABS(y)>COORD_LIMIT ||
            FABS(xp)>SLOPE_LIMIT || FABS(yp)>SLOPE_LIMIT || is_lost) {
            SWAP_PTR(particle[i_part], particle[i_top]);
            if (accepted)
                SWAP_PTR(accepted[i_part], accepted[i_top]);
            particle[i_top][4] = z_start;
            particle[i_top][5] = Po*(1+particle[i_top][5]);
            i_top--;
            i_part--;
            continue;
            }

        }
    log_exit("multipole_tracking");
    return(i_top+1);
    }


double *expansion_coefficients(long n)
{
    static double **expansion_coef=NULL;
    static long *order=NULL;
    static long n_expansions = 0;
    long i;

    log_entry("expansion_coefficients");

    /* look and see if this is already stored */
    for (i=0; i<n_expansions; i++)
        if (n==order[i]) {
            log_exit("expansion_coefficients");
            return(expansion_coef[i]);
            }

    expansion_coef = trealloc(expansion_coef, sizeof(*expansion_coef)*(n_expansions+1));
    order         = trealloc(order, sizeof(*order)*(n_expansions+1));
    expansion_coef[n_expansions] = tmalloc(sizeof(**expansion_coef)*(n+1));
    order[n_expansions] = n;

    /* calculate expansion coefficients with signs for (x+iy)^n/n! */
#ifdef DEBUG
    fprintf(stdout, "coefficients of expansion for multipole of order %ld\n", n);
    fflush(stdout);
#endif
    for (i=0; i<=n; i++) {
        expansion_coef[n_expansions][i] = (ODD(i/2)?-1.0:1.0)/(dfactorial(i)*dfactorial(n-i));
#ifdef DEBUG
        fprintf(stdout, "%.16lf*%sx^%ld*y^%ld \n", expansion_coef[n_expansions][i]*dfactorial(n),
                (ODD(i)?"i*":""), n-i, i);
        fflush(stdout);
#endif
        }             
    log_exit("expansion_coefficients");
    return(expansion_coef[n_expansions++]);
    }

long multipole_tracking2(
                         double **particle,   /* initial/final phase-space coordinates */
                         long n_part,         /* number of particles */
                         ELEMENT_LIST *elem,  /* element pointer */
                         double p_error,      /* p_nominal/p_central */
                         double Po,
                         double **accepted,
                         double z_start
                         )
{
  double KnL;         /* integrated strength = L/(B.rho)*(Dx^n(By))_o for central momentum */
  double dx, dy;      /* offsets of the multipole center */
  long order;         /* order (n) */
  long n_kicks, integ_order;
  long i_part, i_top, n_parts;
  double *coef, *coord;
  double drift, cos_tilt, sin_tilt;
  double tilt, rad_coef, xkick, ykick;
  KQUAD *kquad;
  KSEXT *ksext;
  MULTIPOLE_DATA *multData, *steeringMultData;
  
  log_entry("multipole_tracking2");

  if (!particle)
    bomb("particle array is null (multipole_tracking)", NULL);

  if (!elem)
    bomb("null element pointer (multipole_tracking2)", NULL);
  if (!elem->p_elem)
    bomb("null p_elem pointer (multipole_tracking2)", NULL);

  rad_coef = xkick = ykick = 0;
  steeringMultData = NULL;

  switch (elem->type) {
  case T_KQUAD:
    kquad = ((KQUAD*)elem->p_elem);
    n_kicks = kquad->n_kicks;
    order = 1;
    if (kquad->bore)
      /* KnL = d^nB/dx^n * L/(B.rho) = n! B(a)/a^n * L/(B.rho) * (1+FSE) */
      KnL = kquad->B/kquad->bore*(e_mks/(me_mks*c_mks*Po))*kquad->length*(1+kquad->fse);
    else
      KnL = kquad->k1*kquad->length*(1+kquad->fse);
    drift = kquad->length;
    tilt = kquad->tilt;
    dx = kquad->dx;
    dy = kquad->dy;
    xkick = kquad->xkick;
    ykick = kquad->ykick;
    integ_order = kquad->integration_order;
    if (kquad->synch_rad)
      rad_coef = sqr(e_mks)*pow3(Po)/(6*PI*epsilon_o*sqr(c_mks)*me_mks); 
    if (!kquad->multipolesInitialized) {
      /* read the data files for the error multipoles */
      readErrorMultipoleData(&(kquad->systematicMultipoleData),
                             kquad->systematic_multipoles, 0);
      readErrorMultipoleData(&(kquad->randomMultipoleData),
                             kquad->random_multipoles, 0);
      readErrorMultipoleData(&(kquad->steeringMultipoleData), 
                             kquad->steering_multipoles, 1);
      kquad->multipolesInitialized = 1;
    }
    computeTotalErrorMultipoleFields(&(kquad->totalMultipoleData),
                                     &(kquad->systematicMultipoleData),
                                     &(kquad->randomMultipoleData),
                                     &(kquad->steeringMultipoleData),
                                     KnL, 1);
    multData = &(kquad->totalMultipoleData);
    steeringMultData = &(kquad->steeringMultipoleData);
    break;
  case T_KSEXT:
    ksext = ((KSEXT*)elem->p_elem);
    n_kicks = ksext->n_kicks;
    order = 2;
    if (ksext->bore)
      /* KnL = d^nB/dx^n * L/(B.rho) = n! B(a)/a^n * L/(B.rho) * (1+FSE) */
      KnL = 2*ksext->B/sqr(ksext->bore)*(e_mks/(me_mks*c_mks*Po))*ksext->length*(1+ksext->fse);
    else
      KnL = ksext->k2*ksext->length*(1+ksext->fse);
    drift = ksext->length;
    tilt = ksext->tilt;
    dx = ksext->dx;
    dy = ksext->dy;
    integ_order = ksext->integration_order;
    if (ksext->synch_rad)
      rad_coef = sqr(e_mks)*pow3(Po)/(6*PI*epsilon_o*sqr(c_mks)*me_mks);
    if (!ksext->multipolesInitialized) {
      /* read the data files for the error multipoles */
      readErrorMultipoleData(&(ksext->systematicMultipoleData),
                             ksext->systematic_multipoles, 0);
      readErrorMultipoleData(&(ksext->randomMultipoleData),
                             ksext->random_multipoles, 0);
      ksext->multipolesInitialized = 1;
    }
    computeTotalErrorMultipoleFields(&(ksext->totalMultipoleData),
                                     &(ksext->systematicMultipoleData),
                                     &(ksext->randomMultipoleData),
                                     NULL,
                                     KnL, 2);
    multData = &(ksext->totalMultipoleData);
    break;
  default:
    fprintf(stdout, "error: multipole_tracking2() called for element %s--not supported!\n", elem->name);
    fflush(stdout);
    exit(1);
    break;
  }
  if (!multData->initialized)
    multData = NULL;
  
  if (n_kicks<=0)
    bomb("n_kicks<=0 in multipole()", NULL);
  if (order<=0)
    bomb("order <= 0 in multipole()", NULL);
  if (integ_order!=2 && integ_order!=4) 
    bomb("multipole integration_order must be 2 or 4", NULL);
  
  if (!(coef = expansion_coefficients(order)))
    bomb("expansion_coefficients() returned null pointer (multipole_tracking)", NULL);

  cos_tilt = cos(tilt);
  sin_tilt = sin(tilt);

  i_top = n_part-1;
  if (integ_order==4) {
    if ((n_parts = ceil(n_kicks/4.0))<1)
      n_parts = 1;
    n_kicks = n_parts*4;
  }
  else
    n_parts = n_kicks;
  multipoleKicksDone += (i_top+1)*n_kicks;
  if (multData)
    multipoleKicksDone += (i_top+1)*n_kicks*multData->orders;
  for (i_part=0; i_part<=i_top; i_part++) {
    if (!(coord = particle[i_part])) {
      fprintf(stdout, "null coordinate pointer for particle %ld (multipole_tracking)", i_part);
      fflush(stdout);
      abort();
    }
    if (accepted && !accepted[i_part]) {
      fprintf(stdout, "null accepted coordinates pointer for particle %ld (multipole_tracking)", i_part);
      fflush(stdout);
      abort();
    }

    if ((integ_order==4 &&
         !integrate_kick_multipole_ord4(coord, cos_tilt, sin_tilt, dx, dy, xkick, ykick,
                                        Po, rad_coef, order, KnL, n_parts, drift, multData, steeringMultData)) ||
        (integ_order==2 &&
         !integrate_kick_multipole_ord2(coord, cos_tilt, sin_tilt, dx, dy, xkick, ykick,
                                        Po, rad_coef, order, KnL, n_parts, drift, multData, steeringMultData))) {
      SWAP_PTR(particle[i_part], particle[i_top]);
      if (accepted)
        SWAP_PTR(accepted[i_part], accepted[i_top]);
      particle[i_top][4] = z_start;
      particle[i_top][5] = Po*(1+particle[i_top][5]);
      i_top--;
      i_part--;
      continue;
    }
  }
  log_exit("multipole_tracking2");
  return(i_top+1);
}

int integrate_kick_multipole_ord2(double *coord, double cos_tilt, double sin_tilt,
                                  double dx, double dy, double xkick, double ykick,
                                  double Po, double rad_coef,
                                  long order, double KnL, long n_kicks, double drift,
                                  MULTIPOLE_DATA *multData, 
                                  MULTIPOLE_DATA *steeringMultData) 
{
  double p, qx, qy, denom, beta0, beta1, dp, s;
  double x, y, xp, yp, sum_Fx, sum_Fy;
  long i_kick, imult;
  
  drift = drift/n_kicks/2.0;
  KnL = KnL/n_kicks;
  
  /* calculate coordinates in rotated and offset frame */
  coord[0] -= dx;
  coord[2] -= dy;
  x  =   cos_tilt*coord[0] + sin_tilt*coord[2];
  y  = - sin_tilt*coord[0] + cos_tilt*coord[2];
  xp =   cos_tilt*coord[1] + sin_tilt*coord[3];
  yp = - sin_tilt*coord[1] + cos_tilt*coord[3];
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

  /* apply steering corrector kick */
  xp += xkick/(1+dp)/2;
  yp += ykick/(1+dp)/2;

  /* calculate initial canonical momenta */
  qx = (1+dp)*xp/(denom=sqrt(1+sqr(xp)+sqr(yp)));
  qy = (1+dp)*yp/denom;

  if (steeringMultData) {
    /* apply steering corrector multipoles */
    for (imult=0; imult<steeringMultData->orders; imult++) {
      apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, x, y, 
                                      steeringMultData->order[imult], 
                                      steeringMultData->KnL[imult]*xkick/2, 0);
      apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, x, y, 
                                      steeringMultData->order[imult], 
                                      steeringMultData->JnL[imult]*ykick/2, 1);
    }
    if ((denom=sqr(1+dp)-sqr(qx)-sqr(qy))<=0) {
      return 0;
    }
    xp = qx/(denom=sqrt(denom));
    yp = qy/denom;
  }

  for (i_kick=0; i_kick<n_kicks; i_kick++) {
    if (drift) {
      x += xp*drift*(i_kick?2:1);
      y += yp*drift*(i_kick?2:1);
      s += drift*(i_kick?2:1)*sqrt(1 + sqr(xp) + sqr(yp));
    }
    apply_canonical_multipole_kicks(&qx, &qy, &sum_Fx, &sum_Fy, x, y, order, KnL, 0);

    /* do kicks for spurious multipoles */
    if (multData) {
      for (imult=0; imult<multData->orders; multData++) {
        if (multData->KnL && multData->KnL[imult]) 
          apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, x, y, 
                                          multData->order[imult], 
                                          multData->KnL[imult]/n_kicks, 0);
        if (multData->JnL && multData->JnL[imult]) 
          apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, x, y, 
                                          multData->order[imult], 
                                          multData->JnL[imult]/n_kicks, 1);
      }
    }

    if ((denom=sqr(1+dp)-sqr(qx)-sqr(qy))<=0) {
      return 0;
    }
    xp = qx/(denom=sqrt(denom));
    yp = qy/denom;
    if (rad_coef && drift) {
      qx /= (1+dp);
      qy /= (1+dp);
      dp -= rad_coef*sqr(KnL*(1+dp))*(sqr(sum_Fy)+sqr(sum_Fx))*sqrt(1+sqr(xp)+sqr(yp))/(2*drift);
      qx *= (1+dp);
      qy *= (1+dp);
    }
  }
  if (drift) {
    /* go through final drift */
    x += xp*drift;
    y += yp*drift;
    s += drift*sqrt(1 + sqr(xp) + sqr(yp));
  }

  if (steeringMultData) {
    /* apply steering corrector multipoles */
    for (imult=0; imult<steeringMultData->orders; imult++) {
      apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, x, y, 
                                      steeringMultData->order[imult], 
                                      steeringMultData->KnL[imult]*xkick/2, 0);
      apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, x, y, 
                                      steeringMultData->order[imult], 
                                      steeringMultData->JnL[imult]*ykick/2, 1);
    }
    if ((denom=sqr(1+dp)-sqr(qx)-sqr(qy))<=0) {
      return 0;
    }
    xp = qx/(denom=sqrt(denom));
    yp = qy/denom;
  }

  /* apply steering corrector kick */
  xp += xkick/(1+dp)/2;
  yp += ykick/(1+dp)/2;

  /* undo the rotation and store in place of initial coordinates */
  coord[0] = cos_tilt*x  - sin_tilt*y ;
  coord[2] = sin_tilt*x  + cos_tilt*y ;
  coord[1] = cos_tilt*xp - sin_tilt*yp;
  coord[3] = sin_tilt*xp + cos_tilt*yp;
  if (rad_coef) {
    p = Po*(1+dp);
    beta1 = p/sqrt(sqr(p)+1);
    coord[4] = beta1*(coord[4]/beta0 + 2*s/(beta0+beta1));
  }
  else 
    coord[4] += s;
  coord[5] = dp;

  /* remove the coordinate offsets */
  coord[0] += dx;
  coord[2] += dy;
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


/* BETA is 2^(1/3) */
#define BETA 1.25992104989487316477

int integrate_kick_multipole_ord4(double *coord, double cos_tilt, double sin_tilt,
                                  double dx, double dy, double xkick, double ykick,
                                  double Po, double rad_coef,
                                  long order, double KnL, long n_parts, double drift,
                                  MULTIPOLE_DATA *multData, MULTIPOLE_DATA *steeringMultData) 
{
  double p, qx, qy, denom, beta0, beta1, dp, s;
  double x, y, xp, yp, sum_Fx, sum_Fy;
  long i_kick, step, imult;
  double dsh;
  static double driftFrac[4] = {
    0.5/(2-BETA),  (1-BETA)/(2-BETA)/2,  (1-BETA)/(2-BETA)/2,  0.5/(2-BETA)
    } ;
  static double kickFrac[4] = {
    1./(2-BETA),  -BETA/(2-BETA),  1/(2-BETA),  0
    } ;
  
  drift = drift/n_parts;
  KnL = KnL/n_parts;
  
  /* calculate coordinates in rotated and offset frame */
  coord[0] -= dx;
  coord[2] -= dy;
  x  =   cos_tilt*coord[0] + sin_tilt*coord[2];
  y  = - sin_tilt*coord[0] + cos_tilt*coord[2];
  xp =   cos_tilt*coord[1] + sin_tilt*coord[3];
  yp = - sin_tilt*coord[1] + cos_tilt*coord[3];
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

  /* apply steering corrector kick */
  xp += xkick/(1+dp)/2;
  yp += ykick/(1+dp)/2;

  /* calculate initial canonical momenta */
  qx = (1+dp)*xp/(denom=sqrt(1+sqr(xp)+sqr(yp)));
  qy = (1+dp)*yp/denom;

  if (steeringMultData) {
    /* apply steering corrector multipoles */
    for (imult=0; imult<steeringMultData->orders; imult++) {
      if (steeringMultData->KnL[imult])
        apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, x, y, 
                                        steeringMultData->order[imult], 
                                        steeringMultData->KnL[imult]*xkick/2, 0);
      if (steeringMultData->JnL[imult]) 
        apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, x, y, 
                                        steeringMultData->order[imult], 
                                        steeringMultData->JnL[imult]*ykick/2, 1);
    }
    if ((denom=sqr(1+dp)-sqr(qx)-sqr(qy))<=0) {
      return 0;
    }
    xp = qx/(denom=sqrt(denom));
    yp = qy/denom;
  }

  for (i_kick=0; i_kick<n_parts; i_kick++) {
    for (step=0; step<4; step++) {
      if (drift) {
        dsh = drift*driftFrac[step];
        x += xp*dsh;
        y += yp*dsh;
        s += dsh*sqrt(1 + sqr(xp) + sqr(yp));
      }
      if (!kickFrac[step])
        break;
      apply_canonical_multipole_kicks(&qx, &qy, &sum_Fx, &sum_Fy, x, y, 
                                      order, KnL*kickFrac[step], 0);
      if (multData) {
        /* do kicks for spurious multipoles */
        for (imult=0; imult<multData->orders; imult++) {
          if (multData->KnL && multData->KnL[imult]) {
            apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, x, y, 
                                            multData->order[imult], 
                                            multData->KnL[imult]*kickFrac[step]/n_parts,
                                            0);
          }
          if (multData->JnL && multData->JnL[imult]) {
            apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, x, y, 
                                            multData->order[imult], 
                                            multData->JnL[imult]*kickFrac[step]/n_parts,
                                            1);
          }
        }
      }
      if ((denom=sqr(1+dp)-sqr(qx)-sqr(qy))<=0) {
        return 0;
      }
      xp = qx/(denom=sqrt(denom));
      yp = qy/denom;
      if (rad_coef && drift) {
        qx /= (1+dp);
        qy /= (1+dp);
        dp -= rad_coef*sqr(KnL*kickFrac[step]*(1+dp))*
          (sqr(sum_Fy)+sqr(sum_Fx))*sqrt(1+sqr(xp)+sqr(yp))/(2*drift*kickFrac[step]);
        qx *= (1+dp);
        qy *= (1+dp);
      }
    }
  }
  
  if (steeringMultData) {
    /* apply steering corrector multipoles */
    for (imult=0; imult<steeringMultData->orders; imult++) {
      if (steeringMultData->KnL[imult]) 
        apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, x, y, 
                                        steeringMultData->order[imult], 
                                        steeringMultData->KnL[imult]*xkick/2, 0);
      if (steeringMultData->JnL[imult]) 
        apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, x, y, 
                                        steeringMultData->order[imult], 
                                        steeringMultData->JnL[imult]*ykick/2, 1);
    }
    if ((denom=sqr(1+dp)-sqr(qx)-sqr(qy))<=0) {
      return 0;
    }
    xp = qx/(denom=sqrt(denom));
    yp = qy/denom;
  }

  /* apply steering corrector kick */
  xp += xkick/(1+dp)/2;
  yp += ykick/(1+dp)/2;

  /* undo the rotation and store in place of initial coordinates */
  coord[0] = cos_tilt*x  - sin_tilt*y ;
  coord[2] = sin_tilt*x  + cos_tilt*y ;
  coord[1] = cos_tilt*xp - sin_tilt*yp;
  coord[3] = sin_tilt*xp + cos_tilt*yp;
  if (rad_coef) {
    p = Po*(1+dp);
    beta1 = p/sqrt(sqr(p)+1);
    coord[4] = beta1*(coord[4]/beta0 + 2*s/(beta0+beta1));
  }
  else 
    coord[4] += s;
  coord[5] = dp;

  /* remove the coordinate offsets */
  coord[0] += dx;
  coord[2] += dy;
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

void apply_canonical_multipole_kicks(double *qx, double *qy, 
                                     double *sum_Fx_return, double *sum_Fy_return,
                                     double x, double y,
                                     long order, double KnL, long skew)
{
  long i;
  double sum_Fx, sum_Fy, xypow, ratio;
  double *coef;
  coef = expansion_coefficients(order);
  if (x==0) {
    if (y==0)
      return;
    xypow = ipow(y, order);
    i = order;
    ratio = 0;
  }
  else {
    xypow = ipow(x, order);
    ratio = y/x;
    i = 0;
  }
  /* now sum up the terms for the multipole expansion */
  for (sum_Fx=sum_Fy=0; i<=order; i++) {
    if (ODD(i))
      sum_Fx += coef[i]*xypow;
    else
      sum_Fy += coef[i]*xypow;
    xypow *= ratio;
  }
  if (skew) {
    SWAP_DOUBLE(sum_Fx, sum_Fy);
    sum_Fx = -sum_Fx;
  }
  /* add the kicks */
  *qx -= KnL*sum_Fy;
  *qy += KnL*sum_Fx;
  if (sum_Fx_return)
    *sum_Fx_return = sum_Fx;
  if (sum_Fy_return)
    *sum_Fy_return = sum_Fy;
}

void randomizeErrorMultipoleFields(MULTIPOLE_DATA *randomMult)
{
  long i;
  double nFactorial, rpow, rn1, rn2;

  if (!randomMult || randomMult->randomized)
    return;
  for (i=0; i<randomMult->orders; i++) {
    nFactorial = dfactorial(randomMult->order[i]);
    rpow = ipow(randomMult->referenceRadius, randomMult->order[i]);
    rn1 = gauss_rn_lim(0.0, 1.0, 2.0, random_1);
    rn2 = gauss_rn_lim(0.0, 1.0, 2.0, random_1);
    randomMult->anMod[i] = randomMult->an[i]*nFactorial*rn1/rpow;
    randomMult->bnMod[i] = randomMult->bn[i]*nFactorial*rn2/rpow;
  }
  randomMult->randomized = 1;
}

void computeTotalErrorMultipoleFields(MULTIPOLE_DATA *totalMult,
                                      MULTIPOLE_DATA *systematicMult,
                                      MULTIPOLE_DATA *randomMult,
                                      MULTIPOLE_DATA *steeringMult,
                                      double KmL, long rootOrder)
{
  long i;
  double sFactor=0.0, rFactor=0.0;
  
  if (!totalMult->initialized) {
    totalMult->initialized = 1;
    /* make a list of unique orders for random and systematic multipoles */
    if (systematicMult->orders && randomMult->orders &&
        systematicMult->orders!=randomMult->orders)
      bomb("The number of systematic and random multipole error orders must be the same for any give element", NULL);
    if (systematicMult->orders)
      totalMult->orders = systematicMult->orders;
    else
      totalMult->orders = randomMult->orders;
    if (!(totalMult->order=SDDS_Malloc(sizeof(*totalMult->order)*(totalMult->orders))))
      bomb("memory allocation failure (computeTotalMultipoleFields)", NULL);
    if (systematicMult->orders &&
        (!(systematicMult->anMod=SDDS_Malloc(sizeof(*systematicMult->anMod)*systematicMult->orders)) ||
         !(systematicMult->bnMod=SDDS_Malloc(sizeof(*systematicMult->bnMod)*systematicMult->orders)) ||
         !(systematicMult->KnL=SDDS_Malloc(sizeof(*systematicMult->KnL)*systematicMult->orders)) ||
         !(systematicMult->JnL=SDDS_Malloc(sizeof(*systematicMult->JnL)*systematicMult->orders))))
      bomb("memory allocation failure (computeTotalMultipoleFields)", NULL);
    if (randomMult->orders &&
        (!(randomMult->anMod=SDDS_Malloc(sizeof(*randomMult->anMod)*randomMult->orders)) ||
         !(randomMult->bnMod=SDDS_Malloc(sizeof(*randomMult->bnMod)*randomMult->orders)) ||
         !(randomMult->KnL=SDDS_Malloc(sizeof(*randomMult->KnL)*randomMult->orders)) ||
         !(randomMult->JnL=SDDS_Malloc(sizeof(*randomMult->JnL)*randomMult->orders))))
      bomb("memory allocation failure (computeTotalMultipoleFields", NULL);
    if (!(totalMult->KnL = SDDS_Malloc(sizeof(*totalMult->KnL)*totalMult->orders)) ||
        !(totalMult->JnL = SDDS_Malloc(sizeof(*totalMult->JnL)*totalMult->orders)) )
      bomb("memory allocation failure (computeTotalMultipoleFields)", NULL);
    if (steeringMult && steeringMult->orders) {
      if (!(steeringMult->KnL = SDDS_Malloc(sizeof(*steeringMult->KnL)*steeringMult->orders)) ||
          !(steeringMult->JnL = SDDS_Malloc(sizeof(*steeringMult->JnL)*steeringMult->orders)) )
        bomb("memory allocation failure (computeTotalMultipoleFields)", NULL);
    }
    for (i=0; i<totalMult->orders; i++) {
      if (systematicMult->orders && randomMult->orders &&
          systematicMult->order[i]!=randomMult->order[i])
        bomb("multipole orders in systematic and random lists must match up for any given element.",
             NULL);
      if (systematicMult->orders) {
        totalMult->order[i] = systematicMult->order[i] ;
        systematicMult->anMod[i] = systematicMult->an[i]*dfactorial(systematicMult->order[i])/
          ipow(systematicMult->referenceRadius, systematicMult->order[i]);
        systematicMult->bnMod[i] = systematicMult->bn[i]*dfactorial(systematicMult->order[i])/
          ipow(systematicMult->referenceRadius, systematicMult->order[i]);
      } else {
        totalMult->order[i] = randomMult->order[i];
        /* anMod and bnMod will be computed later for randomized multipoles */
      }
    }
  }

  if (randomMult->orders)
    randomizeErrorMultipoleFields(randomMult);
  
  /* compute normal (KnL) and skew (JnL) from an and bn
   * KnL = an*n!/r^n*(KmL*r^m/m!), 
   * JnL = bn*n!/r^n*(KmL*r^m/m!), where m is the root order 
   * of the magnet with strength KmL
   */
  if (systematicMult->orders)
    sFactor = KmL/dfactorial(rootOrder)*ipow(systematicMult->referenceRadius, rootOrder);
  if (randomMult->orders)
    rFactor = KmL/dfactorial(rootOrder)*ipow(randomMult->referenceRadius, rootOrder);
  for (i=0; i<totalMult->orders; i++) {
    totalMult->KnL[i] = totalMult->JnL[i] = 0;
    if (systematicMult->orders) {
      totalMult->KnL[i] += sFactor*systematicMult->anMod[i];
      totalMult->JnL[i] += sFactor*systematicMult->bnMod[i];
    }
    if (randomMult->orders) {
      totalMult->KnL[i] += rFactor*randomMult->anMod[i];
      totalMult->JnL[i] += rFactor*randomMult->bnMod[i];
    }
  }
  if (steeringMult) {
    /* same for steering multipoles, but compute KnL/theta and JnL/theta (in this case m=0) */
    for (i=0; i<steeringMult->orders; i++) {
      steeringMult->KnL[i] = 
        -1*steeringMult->an[i]*dfactorial(steeringMult->order[i])/ipow(steeringMult->referenceRadius, steeringMult->order[i]);
      steeringMult->JnL[i] = 
        -1*steeringMult->bn[i]*dfactorial(steeringMult->order[i])/ipow(steeringMult->referenceRadius, steeringMult->order[i]);
    }
  }
}

