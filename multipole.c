/* Copyright 1994 by Michael Borland and Argonne National Laboratory,
 * all rights reserved.
 */
/* routine: multipole()
 * purpose: apply kicks due to an arbitrary multipole
 * 
 * Michael Borland, 1991. 
 */
#include "mdb.h"
#include "track.h"

double *expansion_coefficients(long n);

#define ODD(j) (j%2)

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
        KnL = dfactorial(multipole->order)*multipole->BnL/ipow(multipole->bore, multipole->order)*
              (e_mks/(me_mks*c_mks*Po))*multipole->factor;
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
    for (i_part=0; i_part<=i_top; i_part++) {
        if (!(coord = particle[i_part])) {
            fprintf(stderr, "null coordinate pointer for particle %ld (multipole_tracking)", i_part);
            abort();
            }
        if (accepted && !accepted[i_part]) {
            fprintf(stderr, "null accepted coordinates pointer for particle %ld (multipole_tracking)", i_part);
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
#if DEBUG
    fprintf(stderr, "coefficients of expansion for multipole of order %ld\n", n);
#endif
    for (i=0; i<=n; i++) {
        expansion_coef[n_expansions][i] = (ODD(i/2)?-1.0:1.0)/(dfactorial(i)*dfactorial(n-i));
#if DEBUG
        fprintf(stderr, "%.16lf*%sx^%ld*y^%ld \n", expansion_coef[n_expansions][i]*dfactorial(n),
                (ODD(i)?"i*":""), n-i, i);
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
    long n_kicks;       /* number of kicks to split multipole into */
    long i_part, i_kick, i, i_top, is_lost;
    double sum_Fx, sum_Fy, xypow, denom, qx, qy;
    double *coord;
    double drift, cos_tilt, sin_tilt;
    double *coef;
    double x, xp, y, yp, dp, s;
    double ratio, tilt;
    double rad_coef;
    double beta0, beta1, p;
    KQUAD *kquad;
    KSEXT *ksext;

    log_entry("multipole_tracking2");

    if (!particle)
        bomb("particle array is null (multipole_tracking)", NULL);

    if (!elem)
        bomb("null element pointer (multipole_tracking2)", NULL);
    if (!elem->p_elem)
        bomb("null p_elem pointer (multipole_tracking2)", NULL);

    rad_coef = 0;
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
        if (kquad->synch_rad)
            rad_coef = sqr(e_mks)*pow3(Po)/(6*PI*epsilon_o*sqr(c_mks)*me_mks);
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
        if (ksext->synch_rad)
            rad_coef = sqr(e_mks)*pow3(Po)/(6*PI*epsilon_o*sqr(c_mks)*me_mks);
        break;
      default:
        fprintf(stderr, "error: multipole_tracking2() called for element %s--not supported!\n", elem->name);
        exit(1);
        break;
        }

    if (n_kicks<=0)
        bomb("n_kicks<=0 in multipole()", NULL);

    if (order<=0)
        bomb("order <= 0 in multipole()", NULL);

    if (!(coef = expansion_coefficients(order)))
        bomb("expansion_coefficients() returned null pointer (multipole_tracking)", NULL);

    drift = drift/n_kicks/2;
    KnL = KnL/n_kicks;

    cos_tilt = cos(tilt);
    sin_tilt = sin(tilt);

    i_top = n_part-1;
    for (i_part=0; i_part<=i_top; i_part++) {
        if (!(coord = particle[i_part])) {
            fprintf(stderr, "null coordinate pointer for particle %ld (multipole_tracking)", i_part);
            abort();
            }
        if (accepted && !accepted[i_part]) {
            fprintf(stderr, "null accepted coordinates pointer for particle %ld (multipole_tracking)", i_part);
            abort();
            }

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
                s += drift*(i_kick?2:1)*sqrt(1 + sqr(xp) + sqr(yp));
                }
            if (x==0) {
                if (y==0)
                    continue;
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
        coord[0] += dx;
        coord[2] += dy;
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
    log_exit("multipole_tracking2");
    return(i_top+1);
    }
  
