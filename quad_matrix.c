/* Copyright 1994 by Michael Borland and Argonne National Laboratory,
 * all rights reserved.
 */
/* contents: quadrupole_matrix(), qfringe_matrix(), quad_fringe(),
 *           qfringe_R_matrix(), qfringe_T_matrix()
 *
 * Michael Borland, 1989.
 */
#include "track.h"
#include "mdb.h"

static double swap_tmp;
#define swap_double(x, y) (swap_tmp=(x),(x)=(y),(y)=swap_tmp)

VMATRIX *quadrupole_matrix(double K1, double l, long maximum_order,
                           double tilt, double ffringe, double fse)
{
    VMATRIX *M;
    double *C, **R, ***T;
    double kl, k, sin_kl, cos_kl, cosh_kl, sinh_kl;
    double ltot;

    log_entry("quadrupole_matrix");

    K1 *= (1+fse);

    if (K1==0 || l==0) {
        log_exit("quadrupole_matrix");
        return(drift_matrix(l, maximum_order));
        }

    M = tmalloc(sizeof(*M));
    initialize_matrices(M, M->order = MIN(2,maximum_order));
    R = M->R;
    C = M->C;

    ltot = l;
    if (ffringe) 
        l = ltot*(1-ffringe);

    kl = (k=sqrt(fabs(K1)))*l;
    sin_kl  = sin(kl);
    cos_kl  = cos(kl);
    cosh_kl = cosh(kl);
    sinh_kl = sinh(kl);

    R[4][4] = R[5][5] = 1;
    C[4] = l;
    if (K1>0) {
        /* focussing in horizontal plane */
        R[0][0] = R[1][1] = cos_kl;
        R[0][1] = sin_kl/k;
        R[1][0] = -k*sin_kl;
        R[2][2] = R[3][3] = cosh_kl;
        R[2][3] = sinh_kl/k;
        R[3][2] = k*sinh_kl;
        
        if (M->order>=2) {
            T = M->T;
            T[0][5][0] = T[1][5][1] = kl*sin_kl/2;
            T[0][5][1] = sin_kl/(2*k) - l*cos_kl/2;
            T[1][5][0] = k/2*(kl*cos_kl + sin_kl);
            T[2][5][2] = -kl/2*sinh_kl;
            T[2][5][3] = (sinh_kl/k - l*cosh_kl)/2;
            T[3][5][2] = -k/2*(kl*cosh_kl + sinh_kl);
            T[3][5][3] = -kl/2*sinh_kl;
            T[4][0][0] = sqr(k)*(l - sin_kl/k*cos_kl)/4;
            T[4][1][0] = -sqr(sin_kl)/2;
            T[4][1][1] = (l + sin_kl/k*cos_kl)/4;
            T[4][2][2] = -sqr(k)*(l - sinh_kl/k*cosh_kl)/4;
            T[4][3][2] = sqr(sinh_kl)/2;
            T[4][3][3] = (l + sinh_kl/k*cosh_kl)/4;
            }
        }
    else {
        /* defocussing in horizontal plane */
        R[2][2] = R[3][3] = cos_kl;
        R[2][3] = sin_kl/k;
        R[3][2] = -k*sin_kl;
        R[0][0] = R[1][1] = cosh_kl;
        R[0][1] = sinh_kl/k;
        R[1][0] = k*sinh_kl;
        
        if (M->order>=2) {
            T = M->T;
            T[2][5][2] = T[3][5][3] = kl*sin_kl/2;
            T[2][5][3] = sin_kl/(2*k) - l*cos_kl/2;
            T[3][5][2] = k/2*(kl*cos_kl + sin_kl);
            T[0][5][0] = T[1][5][1] = -kl/2*sinh_kl;
            T[0][5][1] = (sinh_kl/k - l*cosh_kl)/2;
            T[1][5][0] = -k/2*(kl*cosh_kl + sinh_kl);
            T[4][0][0] = -sqr(k)*(l - sinh_kl/k*cosh_kl)/4;
            T[4][1][0] = sqr(sinh_kl)/2;
            T[4][1][1] = (l + sinh_kl/k*cosh_kl)/4;
            T[4][2][2] = sqr(k)*(l - sin_kl/k*cos_kl)/4;
            T[4][3][2] = -sqr(sin_kl)/2;
            T[4][3][3] = (l + sin_kl/k*cos_kl)/4;
            }
        }
    
    if (l = ltot*ffringe/2) {
        VMATRIX *Mfringe, *Mtot, *tmp;

        Mtot = tmalloc(sizeof(*Mtot));
        initialize_matrices(Mtot, M->order);

        /* entrance fringe fields */
        Mfringe = quad_fringe(l, K1, M->order, 0, 0.0);
        concat_matrices(Mtot, M, Mfringe);
        tmp  = Mtot;
        Mtot = M;
        M    = tmp;
        free_matrices(Mfringe); tfree(Mfringe); Mfringe = NULL;

        /* exit fringe fields */
        Mfringe = quad_fringe(l, K1, M->order, 1, 0.0);
        concat_matrices(Mtot, Mfringe, M);
        tmp  = Mtot;
        Mtot = M;
        M    = tmp;
        free_matrices(Mfringe); tfree(Mfringe); Mfringe = NULL;
        free_matrices(Mtot); tfree(Mtot); Mtot = NULL;
        }

    C[4] = ltot;

    tilt_matrices(M, tilt);
    log_exit("quadrupole_matrix");
    return(M);
    }


VMATRIX *quad_fringe(double l, double ko, long order, long reverse, double fse)
{
    VMATRIX *M;

    log_entry("quad_fringe");

    ko *= (1+fse);
    M  = tmalloc(sizeof(*M));
    initialize_matrices(M, order);

    M->C[4] = l;
    M->R[5][5] = M->R[4][4] = 1;

    qfringe_R_matrix(
        &M->R[0][0], &M->R[0][1], &M->R[1][0], &M->R[1][1],  ko/l, l);
    qfringe_R_matrix(
        &M->R[2][2], &M->R[2][3], &M->R[3][2], &M->R[3][3], -ko/l, l);
    if (reverse) {
        swap_double(M->R[0][0], M->R[1][1]);
        swap_double(M->R[2][2], M->R[3][3]);
        }

    if (order>=2) {
        qfringe_T_matrix(
            &M->T[0][5][0], &M->T[0][5][1], &M->T[1][5][0], &M->T[1][5][1],
            &M->T[4][0][0], &M->T[4][1][0], &M->T[4][1][1], 
            ko/l, l, reverse);
        qfringe_T_matrix(
            &M->T[2][5][2], &M->T[2][5][3], &M->T[3][5][2], &M->T[3][5][3], 
            &M->T[4][2][2], &M->T[4][3][2], &M->T[4][3][3],
            -ko/l, l, reverse);
        }

    log_exit("quad_fringe");
    return(M);
    }

void qfringe_R_matrix(
    double *R11, double *R12,
    double *R21, double *R22,
    double dk_dz, double l
    )
{
    double term, l3;
    long n;

    log_entry("qfringe_R_matrix");

    if (!l) {
        *R11 = *R22 = 1;
        *R21 = *R12 = 0;
        log_exit("qfringe_R_matrix");
        return;
        }
    if (!dk_dz) {
        *R11 = *R22 = 1;
        *R21 = 0;
        *R12 = l;
        log_exit("qfringe_R_matrix");
        return;
        }

    l3 = pow3(l);

    /* compute R11, R21 */
    *R11 = *R21 = 0;
    term = 1;
    n = 0;
    do {
        *R11 += term;
        *R21 += n*term;
        term *= -l3*dk_dz/((n+3)*(n+2));
        n += 3;
        } while (FABS(term/(*R11))>1e-16);
    *R21 /= l;

    /* compute R12, R22 */
    *R12 = 0;
    *R22 = 0;
    term = l;
    n = 1;
    do {
        *R12 += term;
        *R22 += n*term;
        term *= -l3*dk_dz/((n+3)*(n+2));
        n += 3;
        } while (FABS(term/(*R12))>1e-16);
    *R22 /= l;
    log_exit("qfringe_R_matrix");
    }

void qfringe_T_matrix(
    double *T116, double *T126, double *T216, double *T226, 
    double *T511, double *T512, double *T522,
    double dk_dz, double l, long reverse
    )
{
    double term, l3, ko, mult;
    long n;

    log_entry("qfringe_T_matrix");

    if (!l || !dk_dz) {
        *T116 = *T226 = *T216 = *T126 = *T511 = *T512 = *T522 = 0;
        log_exit("qfringe_T_matrix");
        return;
        }

    l3 = pow3(l);

    /* compute T116, T216 */
    *T116 = *T216 = 0;
    term = 1;
    n = 0;
    do {
        *T116 += -n*term/3;
        *T216 += -sqr(n)*term/3;
        term *= -l3*dk_dz/((n+3)*(n+2));
        n += 3;
        } while (FABS(term/(*T116?*T116:1))>1e-16);
    *T216 /= l;

    /* compute T126, T226 */
    *T126 = 0;
    *T226 = 0;
    term = l;
    n = 1;
    do {
        *T126 += -(n-1)*term/3;
        *T226 += -n*(n-1)*term/3;
        term *= -l3*dk_dz/((n+3)*(n+2));
        n += 3;
        } while (FABS(term/(*T126?*T126:1))>1e-16);
    *T226 /= l;

    if (reverse) 
        swap_double(*T116, *T226);

    /* compute path-length terms using truncated series (truncated at z^12) */
    if (!reverse) {
        /* entrance fringe field */
        ko = dk_dz*l;
        /* T511 terms */
        term = sqr(ko)*pow3(l);
        if ((mult = ko*sqr(l))>1)
            printf("warning: path-length terms for qfringe may be inaccurate: ko*lf^2>1\n");
        *T511  = 1./20.*term;         term *= mult;
        *T511 += -1./240.*term;       term *= mult; 
        *T511 +=  13./79200.*term;    term *= mult;
        *T511 += -19./4989600.*term;  term *= mult;
        *T511 += 19./325721088.*term;
        /* T522 terms */
        *T522  = term = l;            term *= mult;
        *T522 += -1./6.*term;         term *= mult;
        *T522 += 5./252.*term;        term *= mult;
        *T522 += -11./11340.*term;    term *= mult;
        *T522 += 187./7076160.*term;   term *= mult;
        *T522 += -391./849139200.*term; 
        *T522 /= 2;
        /* T512 terms */
        term = mult;
        *T512  = -1./6.*term;             term *= mult;
        *T512 += 1./30.*term;             term *= mult;
        *T512 += -1./480.*term;           term *= mult;
        *T512 += 1./14784.*term;          term *= mult;
        *T512 += -1./739200.*term; term *= mult;
        *T512 += 1./54454400.*term;
        }
    else {
        /* exit fringe field */
        ko = dk_dz*l;
        /* T511 terms */
        term = sqr(ko)*pow3(l);
        if ((mult = ko*sqr(l))>1)
            printf("warning: path-length terms for qfringe may be inaccurate: ko*lf^2>1\n");
        *T511  = 2./15.*term;        term *= mult;
        *T511 += -1./80.*term;       term *= mult; 
        *T511 +=  1./1848.*term;     term *= mult;
        *T511 += -1./73920*term;     term *= mult;
        *T511 += 3./13613600.*term;
        /* T522 terms */
        *T522  = term = l;            term *= mult;
        *T522 += -1./6.*term;          term *= mult;
        *T522 += 1./70.*term;          term *= mult;
        *T522 += -1./1680.*term;       term *= mult;
        *T522 += 1./68640.*term;   term *= mult;
        *T522 += -3./12812800.*term; 
        *T522 /= 2;
        /* T512 terms */
        term = mult;
        *T512  = -1./3.*term;             term *= mult;
        *T512 += 1./20.*term;             term *= mult;
        *T512 += -1./336.*term;           term *= mult;
        *T512 += 1./10560.*term;          term *= mult;
        *T512 += -3/1601600.*term; term *= mult;
        *T512 += 1./39603200.*term;
        }
    log_exit("qfringe_T_matrix");
    }

/* This subroutine returns a stand-alone quadrpole fringe
 * matrix, in the event the user asks for a qfringe element
 * separate from a quadrupole
 */

VMATRIX *qfringe_matrix(
    double K1, double l, double tilt, long direction, long order, double fse
    )
{
    VMATRIX *M;

    log_entry("qfringe_matrix");

    M = quad_fringe(l, K1, order, (direction==-1?1:0), fse);

    M->C[4] = l;
    M->R[5][5] = M->R[4][4] = 1;

    tilt_matrices(M, tilt);

    log_exit("qfringe_matrix");
    return(M);
    }

