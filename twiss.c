/* Copyright 1994 by Michael Borland and Argonne National Laboratory,
 * all rights reserved.
 */
/* file: twiss.c
 * purpose: computation of Twiss parameters
 *
 * Michael Borland, 1989
 */
#include "mdb.h"
#include "track.h"


void copy_doubles(double *target, double *source, long n);
double find_acceptance(ELEMENT_LIST *elem, long plane, RUN *run, char **name, double *end_pos);
void modify_rfca_matrices(ELEMENT_LIST *eptr, long order);
void compute_twiss_statistics(LINE_LIST *beamline, TWISS *twiss_ave, TWISS *twiss_min, TWISS *twiss_max);

#define TWISS_Y_OFFSET 6
#define TWISS_CENT_OFFSET 12

VMATRIX *compute_periodic_twiss(
    double *betax, double *alphax, double *etax, double *etapx, double *NUx,
    double *betay, double *alphay, double *etay, double *etapy, double *NUy,
    ELEMENT_LIST *elem, double *clorb, RUN *run)
{
    VMATRIX *M, *M1;
    double cos_phi, sin_phi, **R, beta[2], alpha[2], eta[2], etap[2], phi[2];
    double det;
    long i, j;

    log_entry("compute_periodic_twiss");

    if ((i = fill_in_matrices(elem, run)))
        fprintf(stderr, "%ld matrices recomputed for periodic Twiss parameter computation\n", i);
    
    modify_rfca_matrices(elem, run->default_order);  /* replace rf cavities with drifts */
    if (clorb) {
        /* use the closed orbit to compute the on-orbit R matrix */
        M1 = tmalloc(sizeof(*M1));
        initialize_matrices(M1, 1);
        for (i=0; i<6; i++) {
            M1->C[i] = clorb[i];
            M1->R[i][i] = 1;
            }
        /* use 3rd order here to get as close to tracked closed orbit as possible */
        M = append_full_matrix(elem, run, M1, 3);
#if defined(DEBUG)
        fprintf(stderr, "matrix concatenation for periodic Twiss computation:\n");
        fprintf(stderr, "closed orbit at input:\n  ");
        for (i=0; i<6; i++)
            fprintf(stderr, "%14.6e ", clorb[i]);
        fprintf(stderr, "\nclosed orbit at output:\n  ");
        for (i=0; i<6; i++)
            fprintf(stderr, "%14.6e ", M->C[i]);
        fprintf(stderr, "\nR matrix:\n");
        for (i=0; i<6; i++) {
            fprintf(stderr, "  ");
            for (j=0; j<6; j++)
                fprintf(stderr, "%14.6e ", M->R[i][j]);
            fprintf(stderr, "\n");
            }
#endif
        }
    else
        M = full_matrix(elem, run, run->default_order);
    R = M->R;

    for (i=0; i<4; i+=2 ) {
        if (fabs(cos_phi = (R[i][i] + R[i+1][i+1])/2)>1) {
            fprintf(stderr, "warning: beamline unstable for %c plane--can't match beta functions.\n", 'x'+i/2);
            sin_phi = 1e-6;
            }
        beta[i/2] = fabs(R[i][i+1]/sin(acos(cos_phi)));
        sin_phi   = R[i][i+1]/beta[i/2];
        phi[i/2]  = atan2(sin_phi, cos_phi);
        if (phi[i/2]<0)
            phi[i/2] += PIx2;
        alpha[i/2] = (R[i][i]-R[i+1][i+1])/(2*sin_phi);
        if ((det = (2-R[i][i]-R[i+1][i+1]))<=0) {
            fprintf(stderr, "error: beamline unstable for %c plane--can't match dispersion functions.\n",
                'x'+i/2);
            det = 1e-6;
            }
        eta[i/2] = ((1-R[i+1][i+1])*R[i][5]+R[i][i+1]*R[i+1][5])/det;
        etap[i/2] = (R[i+1][i]*R[i][5] + (1-R[i][i])*R[i+1][5])/det;
        }
    *betax = beta[0];
    *betay = beta[1];
    *alphax = alpha[0];
    *alphay = alpha[1];
    *etax = eta[0];
    *etay = eta[1];
    *etapx = etap[0];
    *etapy = etap[1];
    *NUx = phi[0]/PIx2;
    *NUy = phi[1]/PIx2;

    log_exit("compute_periodic_twiss");
    return(M);
    }

void propagate_twiss_parameters(TWISS *twiss0, double *tune, RADIATION_INTEGRALS *radIntegrals, ELEMENT_LIST *elem, 
                                long plane, RUN *run, double *traj
        )
{
    double beta0, alpha0, phi0, eta0, etap0, gamma0;
    double *func, *path, detR;
    double **R, C, S, Cp, Sp, D, Dp, sin_dphi, cos_dphi, dphi;
    double centroid;
    long n_mat_computed, i;
    VMATRIX *M1, *M2;
    double func0, func1;
    ELEMENT_LIST *elem0;
    
    log_entry("propagate_twiss_parameters");

    if (!twiss0)
        bomb("initial Twiss parameters not given (propagate_twiss_parameters())", NULL);

    switch (plane) {
        case 0:
            beta0  = twiss0->betax;
            alpha0 = twiss0->alphax;
            phi0   = twiss0->phix;
            eta0   = twiss0->etax;
            etap0  = twiss0->etapx;
            break;
        default:
            beta0  = twiss0->betay;
            alpha0 = twiss0->alphay;
            phi0   = twiss0->phiy;
            eta0   = twiss0->etay;
            etap0  = twiss0->etapy;
            break;
        }

#if defined(DEBUG)
    fprintf(stderr, "now in propagate_twiss_parameters\n");
    fprintf(stderr, "initial twiss parameters for %c plane:\nbeta = %e, alpha = %e, eta = %e, etap = %e\n",
           plane==0?'x':'y', beta0, alpha0, eta0, etap0);
#endif

    path = tmalloc(sizeof(*path)*6);
    M1 = tmalloc(sizeof(*M1));
    M2 = tmalloc(sizeof(*M2));
    initialize_matrices(M1, 1);
    initialize_matrices(M2, run->default_order);
    if (traj) {
        for (i=0; i<6; i++)
            path[i] = traj[i];
        }
    else {
        for (i=0; i<6; i++)
            path[i] = 0;
        }

    n_mat_computed = 0;

    if (radIntegrals) {
      for (i=0; i<6; i++) 
        radIntegrals->I[i] = 0;
    }
    
    while (elem) {
        gamma0 = (1+sqr(alpha0))/beta0;
        if (entity_description[elem->type].flags&HAS_MATRIX) {
            if (!elem->matrix || !(elem->matrix->R) || (elem->pred && elem->pred->Pref_output!=elem->Pref_input)) {
                elem->matrix = compute_matrix(elem, run, NULL);
                n_mat_computed++;
                }
            /* use matrix concatenation to include effect of beam path */
            for (i=0; i<6; i++) {
                M1->C[i] = path[i];
                M1->R[i][i] = 1;
                }
            concat_matrices(M2, elem->matrix, M1);
            R = M2->R;
            /* record new centroids for beam path */
            for (i=0; i<6; i++)
                path[i] = M2->C[i];
            switch (plane) {
              case 0:
                C  = R[0][0]; S  = R[0][1];
                Cp = R[1][0]; Sp = R[1][1];
                D  = R[0][5]; Dp = R[1][5];
                break;
              default:
                C  = R[2][2]; S  = R[2][3];
                Cp = R[3][2]; Sp = R[3][3];
                D  = R[2][5]; Dp = R[3][5];
                break;
                }
            }
        else {
            if (elem->pred)
                elem->Pref_input = elem->Pref_output = elem->pred->Pref_output;
            if (entity_description[elem->type].flags&HAS_LENGTH) {
                C = Sp = 1;
                Cp = 0;
                D = Dp = 0;
                S = *((double*)elem->p_elem);
                }
            else {
                C = Sp = 1;
                Cp = S = 0;
                D = Dp = 0;
                }
            }
        detR = C*Sp - Cp*S;
        if (!elem->twiss)
            elem->twiss = tmalloc(sizeof(*elem->twiss));
        func = ((double*)elem->twiss) + (plane?TWISS_Y_OFFSET:0);
        /* store centroid position */
        *(((double*)elem->twiss) + (plane?1:0) + TWISS_CENT_OFFSET) = path[plane?2:0];

        if (radIntegrals && plane==0) {
          /* compute contribution to radiation integrals */
          long isBend;
          BEND *bptr;
          KSBEND *kbptr;
          CSBEND *cbptr;
          double length, angle, E1, E2, K1;
          double k2, rho, k, kl, gamma1;
          double I1, I2, I3, I4, I5;
          double alpha1, gamam1, etap1, eta2, sin_kl, cos_kl;
          double etaAve, etaK1_rhoAve, HAve;
          
          isBend = 1;
          switch (elem->type) {
          case T_SBEN:
          case T_RBEN:
            bptr = (BEND*)(elem->p_elem);
            length = bptr->length;
            angle = bptr->angle;
            E1 = bptr->e1*(bptr->edge1_effects?1:0);
            E2 = bptr->e2*(bptr->edge2_effects?1:0);
            K1 = bptr->k1;
            break;
          case T_KSBEND:
            kbptr = (KSBEND*)(elem->p_elem);
            length = kbptr->length;
            angle = kbptr->angle;
            E1 = kbptr->e1*(kbptr->edge1_effects?1:0);
            E2 = kbptr->e2*(kbptr->edge2_effects?1:0);
            K1 = kbptr->k1;
            break;
          case T_CSBEND:
            cbptr = (CSBEND*)(elem->p_elem);
            length = cbptr->length;
            angle = cbptr->angle;
            E1 = cbptr->e1*(cbptr->edge1_effects?1:0);
            E2 = cbptr->e2*(cbptr->edge2_effects?1:0);
            K1 = cbptr->k1;
            break;
          default:
            isBend = 0;
            break;
          }
          if (isBend && angle!=0) {
            rho = length/angle;
            k2 = K1+1./(rho*rho);
            /* equations are from SLAC 1193 */
            if (k2<0) {
              k = sqrt(-k2);
              kl = k*length;
              cos_kl = cosh(kl);
              sin_kl = sinh(kl);
            } else {
              k = sqrt(k2);
              kl = k*length;
              sin_kl = sin(kl);
              cos_kl = cos(kl);
            }
            etap1 = etap0 + eta0/rho*tan(E1);
            eta2  = eta0*cos_kl + etap1*sin_kl/k + (1-cos_kl)/(rho*k2);
            alpha1 = alpha0 - beta0/rho*tan(E1);
            gamma1 = (1+sqr(alpha1))/beta0;
            etaAve = eta0*sin_kl/kl + etap1*(1-cos_kl)/(k2*length) +
              (kl-sin_kl)/(k2*kl*rho);
            etaK1_rhoAve =  -etaAve*K1/rho + (eta0*tan(E1)+eta2*tan(E2))/(2*length*sqr(rho));
            HAve = gamma1*sqr(eta0) + 2*alpha1*eta0*etap1 + beta0*sqr(etap1) 
              + 2*angle*( -(gamma1*eta0+alpha1*etap1)*(kl-sin_kl)/(k*k2*sqr(length)) +
                         (alpha1*eta0+beta0*etap1)*(1-cos_kl)/(k2*sqr(length))
                         )
                + sqr(angle)*(gamma1*(3*kl-4*sin_kl+sin_kl*cos_kl)/(2*k*k2*k2*ipow(length,3)) 
                              - alpha1*sqr(1-cos_kl)/(k2*k2*ipow(length,3))
                              + beta0*(kl-cos_kl*sin_kl)/(2*kl*k2*sqr(length)));
            I1 = etaAve*length/rho;
            I2 = length/sqr(rho);
            I3 = I2/fabs(rho);
            I4 = I2/rho*etaAve - 2*length*etaK1_rhoAve;
            I5 = HAve*I3;
            radIntegrals->I[0] += I1;
            radIntegrals->I[1] += I2;
            radIntegrals->I[2] += I3;
            radIntegrals->I[3] += I4;
            radIntegrals->I[4] += I5;
          }
        }
        
        /* calculate new beta and alpha */
        func[0] = (sqr(C)*beta0 - 2*C*S*alpha0 + sqr(S)*gamma0)/detR;
        func[1] = (-C*Cp*beta0 + (Sp*C+S*Cp)*alpha0 - S*Sp*gamma0)/detR;
#if defined(DEBUG)
        fprintf(stderr, "** working on element %s at z=%e\n", elem->name, elem->end_pos);
        fprintf(stderr, "initial twiss parameters: beta = %e, alpha = %e, gamma = %e, eta = %e, etap = %e\n",
               beta0, alpha0, gamma0, eta0, etap0);
        fprintf(stderr, "matrix: C = %e, S = %e, C' = %e, S' = %e\n",
            C, S, Cp, Sp);
        fprintf(stderr, "func[0] = %e, func[1] = %e \n", func[0], func[1]);
        fprintf(stderr, "func0   = %e, func1   = %e \n",
          sqr(C)*beta0 - 2*C*S*alpha0 + sqr(S)*gamma0,
          -C*Cp*beta0 + (Sp*C+S*Cp)*alpha0 - S*Sp*gamma0);
#endif
        /* use R12=S to find sin(dphi) */
        if ((sin_dphi=S/sqrt(beta0*func[0]))>1) {
            fprintf(stderr, "warning: argument of asin > 1 by %f (propagate_twiss)\n", sin_dphi-1);
            fprintf(stderr, "element is %s at z=%em\n", elem->name, elem->end_pos);
            fprintf(stderr, "%c-plane matrix:  C = %e,  S = %e,  C' = %e,  S' = %e\n",
                   (plane==0?'x':'y'), C, S, Cp, Sp);
            fprintf(stderr, "beta0 = %e, func[0] = %e\n", beta0, func[0]);
            sin_dphi = 1;
            cos_dphi = 0;
            }
        else if (sin_dphi<-1) {
            fprintf(stderr, "warning: argument of asin < -1 by %f (propagate_twiss)\n", sin_dphi+1);
            fprintf(stderr, "element is %s at z=%em\n", elem->name, elem->end_pos);
            sin_dphi = -1;
            cos_dphi = 0;
            }
        else {
            /* use R11=C to find cos(dphi) */
            cos_dphi = sqrt(beta0/func[0])*C - alpha0*sin_dphi;
          }
        if ((dphi = atan2(sin_dphi, cos_dphi))<0) 
          dphi += PIx2;
        phi0   = func[2] = phi0 + dphi;
        beta0  = func[0];
        alpha0 = func[1];
        /* compute new dispersion function and slope */
        func[3] = eta0*C  + etap0*S  + D;
        func[4] = eta0*Cp + etap0*Sp + Dp;
        eta0  = func[3];
        etap0 = func[4];
        elem = elem->succ;
#if defined(DEBUG)
        fprintf(stderr, "new twiss parameters: beta = %e, alpha = %e, eta = %e, etap = %e\n",
               beta0, alpha0, eta0, etap0);
#endif
        }
    *tune = phi0/PIx2;

    if (traj) {
        tfree(path);
        free_matrices(M1); tfree(M1);
        free_matrices(M2); tfree(M2);
        }
    log_exit("propagate_twiss_parameters");
    }


static long twiss_initialized = 0;
static long SDDS_twiss_initialized = 0;
static SDDS_TABLE SDDS_twiss;
static long twiss_count = 0;

#define IC_S 0
#define IC_BETAX 1
#define IC_ALPHAX 2
#define IC_PHIX 3
#define IC_ETAX 4
#define IC_ETAPX 5
#define IC_APX 6
#define IC_BETAY 7
#define IC_ALPHAY 8
#define IC_PHIY 9
#define IC_ETAY 10
#define IC_ETAPY 11
#define IC_APY 12
#define IC_PCENTRAL 13
#define N_DOUBLE_COLUMNS 14
#define IC_ELEMENT 14
#define IC_OCCURENCE 15
#define IC_TYPE 16
#define N_COLUMNS 17
static SDDS_DEFINITION column_definition[N_COLUMNS] = {
  {"s", "&column name=s, type=double, units=m, description=Distance &end"},
  {"betax", "&column name=betax, type=double, units=m, symbol=\"$gb$r$bx$n\", description=\"Horizontal beta-function\" &end"},
  {"alphax", "&column name=alphax, type=double, symbol=\"$ga$r$bx$n\", description=\"Horizontal alpha-function\" &end"},
  {"psix", "&column name=psix, type=double, units=rad, symbol=\"$gy$r$bx$n\", description=\"Horizontal phase advance\" &end"},
  {"etax", "&column name=etax, type=double, units=m, symbol=\"$gc$r$bx$n\", description=\"Horizontal dispersion\" &end"},
  {"etaxp", "&column name=etaxp, type=double, symbol=\"$gc$r$bx$n$a'$n\", description=\"Slope of horizontal dispersion\" &end"},
  {"xAperture", "&column name=xAperture, type=double, units=m, symbol=\"a$bx$n\", description=\"Horizontal aperture\" &end"},
  {"betay", "&column name=betay, type=double, units=m, symbol=\"$gb$r$by$n\", description=\"Vertical beta-function\" &end"},
  {"alphay", "&column name=alphay, type=double, symbol=\"$ga$r$by$n\", description=\"Vertical alpha-function\" &end"},
  {"psiy", "&column name=psiy, type=double, units=rad, symbol=\"$gy$r$by$n\", description=\"Vertical phase advance\" &end"},
  {"etay", "&column name=etay, type=double, units=m, symbol=\"$gc$r$by$n\", description=\"Vertical dispersion\" &end"},
  {"etayp", "&column name=etayp, type=double, symbol=\"$gc$r$by$n$a'$n\", description=\"Slope of vertical dispersion\" &end"},
  {"yAperture", "&column name=yAperture, type=double, units=m, symbol=\"a$by$n\", description=\"Vertical aperture\" &end"},
  {"pCentral", "&column name=pCentral, type=double, units=\"m$be$nc\", symbol=\"p$bcent$n\", description=\"Central momentum\" &end"},
  {"ElementName", "&column name=ElementName, type=string, description=\"Element name\", format_string=%10s &end"},
  {"ElementOccurence", "&column name=ElementOccurence, type=long, description=\"Occurence of element\", format_string=%6ld &end"},
  {"ElementType", "&column name=ElementType, type=string, description=\"Element-type name\", format_string=%10s &end"},
  };

#define IP_STEP 0
#define IP_NUX 1
#define IP_DNUXDP 2
#define IP_AX 3
#define IP_NUY 4
#define IP_DNUYDP 5
#define IP_AY 6
#define IP_STAGE 7
#define IP_ALPHAC 8
#define IP_I1 9
#define IP_I2 10
#define IP_I3 11
#define IP_I4 12
#define IP_I5 13
#define IP_EX0 14
#define IP_TAUX 15
#define IP_JX 16
#define IP_TAUY 17
#define IP_JY 18
#define IP_SIGMADELTA 19
#define IP_TAUDELTA 20
#define IP_JDELTA 21
#define IP_U0 22
#define N_PARAMETERS 23
static SDDS_DEFINITION parameter_definition[N_PARAMETERS] = {
    {"Step", "&parameter name=Step, type=long, description=\"Simulation step\" &end"},
    {"nux", "&parameter name=nux, symbol=\"$gn$r$bx$n\", type=double, units=\"1/(2$gp$r)\", description=\"Horizontal tune\" &end"},
    {"dnux/dp", "&parameter name=dnux/dp, symbol=\"$gx$r$bx$n\", type=double, units=\"1/(2$gp$r)\", description=\"Horizontal chromaticity\" &end"},
    {"Ax", "&parameter name=Ax, symbol=\"A$bx$n\", type=double, units=\"$gp$rm\", description=\"Horizontal acceptance\" &end"},
    {"nuy", "&parameter name=nuy, symbol=\"$gn$r$by$n\", type=double, units=\"1/(2$gp$r)\", description=\"Vertical tune\" &end"},
    {"dnuy/dp", "&parameter name=dnuy/dp, symbol=\"$gx$r$by$n\", type=double, units=\"1/(2$gp$r)\", description=\"Vertical chromaticity\" &end"},
    {"Ay", "&parameter name=Ay, symbol=\"A$by$n\", type=double, units=\"$gp$rm\", description=\"Vertical acceptance\" &end"},
    {"Stage", "&parameter name=Stage, type=string, description=\"Stage of computation\" &end"},
    {"alphac", "&parameter name=alphac, symbol=\"$ga$r$bc$n\", type=double, description=\"Momentum compaction factor\" &end"},
  {"I1", "&parameter name=I1, type=double, description=\"Radiation integral 1\", units=m &end"} ,
  {"I2", "&parameter name=I2, type=double, description=\"Radiation integral 2\", units=1/m &end"} ,
  {"I3", "&parameter name=I3, type=double, description=\"Radiation integral 3\", units=1/m$a2$n &end"} ,
  {"I4", "&parameter name=I4, type=double, description=\"Radiation integral 4\", units=1/m &end"} ,
  {"I5", "&parameter name=I5, type=double, description=\"Radiation integral 5\", units=1/m &end"} ,
  {"ex0", "&parameter name=ex0, type=double, description=\"Damped horizontal emittance\", units=$gp$rm &end"},
  {"taux", "&parameter name=taux, type=double, description=\"Horizontal damping time\", units=s &end"},
  {"Jx", "&parameter name=Jx, type=double, description=\"Horizontal damping partition number\" &end"},
  {"tauy", "&parameter name=tauy, type=double, description=\"Vertical damping time\", units=s &end"},
  {"Jy", "&parameter name=Jy, type=double, description=\"Vertical damping partition number\" &end"},
  {"Sdelta0", "&parameter name=Sdelta0, type=double, description=\"RMS fractional energy spread\" &end"},
  {"taudelta", "&parameter name=taudelta, type=double, description=\"Longitudinal damping time\", units=s &end"},
  {"Jdelta", "&parameter name=Jdelta, type=double, description=\"Longitudinal damping partition number\" &end"},
  {"U0", "&parameter name=U0, type=double, units=MeV, description=\"Energy loss per turn\" &end"},
    } ;

void dump_twiss_parameters(
    TWISS *twiss0,
    ELEMENT_LIST *elem,
    long n_elem,
    double *tune,
    RADIATION_INTEGRALS *radIntegrals,                           
    double *chromaticity,
    double *acceptance,
    double alphac,
    long final_values_only,
    long tune_corrected,
    RUN *run
    )
{
    double *data;
    long i, j, row_count;
    char *stage;

    log_entry("dump_twiss_parameters");

    if (!twiss0)
        bomb("Twiss data not computed prior to dump_twiss_parameters() call (1)", NULL);

    data = tmalloc(sizeof(*data)*N_DOUBLE_COLUMNS);

    if (tune_corrected==1)
        stage = "tunes corrected";
    else
        stage = "tunes uncorrected";

    if (!SDDS_StartTable(&SDDS_twiss, final_values_only?1:n_elem+1)) {
        SDDS_SetError("Problem starting SDDS table (dump_twiss_parameters)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
    if (!SDDS_SetParameters(&SDDS_twiss, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE,
                            IP_STEP, twiss_count, IP_STAGE, stage, 
                            IP_NUX, tune[0], IP_DNUXDP, chromaticity[0], IP_AX, acceptance[0],
                            IP_NUY, tune[1], IP_DNUYDP, chromaticity[1], IP_AY, acceptance[1], 
                            IP_ALPHAC, alphac, -1)) {
        SDDS_SetError("Problem setting SDDS parameters (dump_twiss_parameters 1)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
    if (radIntegrals) {
      if (!SDDS_SetParameters(&SDDS_twiss, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE,
                              IP_I1, radIntegrals->I[0],
                              IP_I2, radIntegrals->I[1],
                              IP_I3, radIntegrals->I[2],
                              IP_I4, radIntegrals->I[3],
                              IP_I5, radIntegrals->I[4],
                              IP_EX0, radIntegrals->ex0,
                              IP_TAUX, radIntegrals->taux,
                              IP_JX, radIntegrals->Jx,
                              IP_TAUY, radIntegrals->tauy,
                              IP_JY, radIntegrals->Jy,
                              IP_SIGMADELTA, radIntegrals->sigmadelta,
                              IP_TAUDELTA, radIntegrals->taudelta,
                              IP_JDELTA, radIntegrals->Jdelta,
                              IP_U0, radIntegrals->Uo, -1)) {
        SDDS_SetError("Problem setting SDDS parameters (dump_twiss_parameters 2)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      }
    }
    
    if (!final_values_only) {
        row_count = 0;
        data[0] = 0;     /* position */
        copy_doubles(data+1, (double*)twiss0, N_DOUBLE_COLUMNS-2);
        data[N_DOUBLE_COLUMNS-1] = elem->Pref_input;
        for (j=0; j<N_DOUBLE_COLUMNS; j++)
            if (!SDDS_SetRowValues(&SDDS_twiss, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, row_count, j, data[j], -1)) {
                SDDS_SetError("Problem setting SDDS rows (dump_twiss_parameters)");
                SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
                }
        if (!SDDS_SetRowValues(&SDDS_twiss, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, row_count++, 
                               IC_ELEMENT, "_BEG_", IC_OCCURENCE, (long)1, IC_TYPE, "MARK", -1)) {
            SDDS_SetError("Problem setting SDDS rows (dump_twiss_parameters)");
            SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
            }

        i = 0;
        while (elem) {
            data[0] = elem->end_pos;     /* position */
            data[N_DOUBLE_COLUMNS-1] = elem->Pref_output;
            if (!elem->twiss)
                bomb("Twiss data not computed prior to dump_twiss_parameters() call (2)", NULL);
            copy_doubles(data+1, (double*)elem->twiss, N_DOUBLE_COLUMNS-2);
            for (j=0; j<N_DOUBLE_COLUMNS; j++)
                if (!SDDS_SetRowValues(&SDDS_twiss, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, row_count, j, data[j], -1)) {
                    SDDS_SetError("Problem setting SDDS rows (dump_twiss_parameters)");
                    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
                    }
            if (!SDDS_SetRowValues(&SDDS_twiss, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, row_count++, 
                                   IC_ELEMENT, elem->name, IC_OCCURENCE, elem->occurence, 
                                   IC_TYPE, entity_name[elem->type], -1)) {
                SDDS_SetError("Problem setting SDDS rows (dump_twiss_parameters)");
                SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
                }
            i++;
            elem = elem->succ;
            }
        if (i!=n_elem)
            bomb("element count error in dump_twiss_parameters()", NULL);
        }
    else {
        /* find final element */
        i = 0;
        while (1) {
            if (!elem->twiss)
                bomb("Twiss data not computed prior to dump_twiss_parameters() call (2)", NULL);
            i++;
            if (!elem->succ)
                break;
            elem = elem->succ;
            }
        if (i!=n_elem)
            bomb("element count error in dump_twiss_parameters()", NULL);
        data[0] = elem->end_pos;     /* position */
        data[N_DOUBLE_COLUMNS-1] = elem->Pref_output;
        copy_doubles(data+1, (double*)elem->twiss, N_DOUBLE_COLUMNS-2);
        for (j=0; j<N_DOUBLE_COLUMNS; j++)
            if (!SDDS_SetRowValues(&SDDS_twiss, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, 0, j, data[j], -1)) {
                SDDS_SetError("Problem setting SDDS rows (dump_twiss_parameters)");
                SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
                }
        if (!SDDS_SetRowValues(&SDDS_twiss, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, 0, 
                               IC_ELEMENT, elem->name, IC_OCCURENCE, elem->occurence, 
                               IC_TYPE, entity_name[elem->type], -1)) {
            SDDS_SetError("Problem setting SDDS rows (dump_twiss_parameters)");
            SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
            }
        }

    if (!SDDS_WriteTable(&SDDS_twiss)) {
        SDDS_SetError("Unable to write Twiss parameter data (dump_twiss_parameters)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
    if (!SDDS_EraseData(&SDDS_twiss)) {
        SDDS_SetError("Unable to erase Twiss parameter data (dump_twiss_parameters)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }

    free(data);
    if (tune_corrected)
        twiss_count++;

    log_exit("dump_twiss_parameters");
    }

#include "twiss.h"

long get_twiss_mode(long *mode, double *x_twiss, double *y_twiss)
{
    if (!twiss_initialized)
        return(0);
    *mode = matched;
    x_twiss[0] = beta_x;
    x_twiss[1] = alpha_x;
    x_twiss[2] = eta_x;
    x_twiss[3] = etap_x;
    x_twiss[4] = 0;
    y_twiss[0] = beta_y;
    y_twiss[1] = alpha_y;
    y_twiss[2] = eta_y;
    y_twiss[3] = etap_y;
    y_twiss[4] = 0;
    return(1);
    }

void setup_twiss_output(NAMELIST_TEXT *nltext, RUN *run, LINE_LIST *beamline, long *do_twiss_output)
{
    char label[200];

    log_entry("setup_twiss_output");

    /* process namelist input */
    set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
    set_print_namelist_flags(0);
    process_namelist(&twiss_output, nltext);
    print_namelist(stderr, &twiss_output);

    if (filename)
        filename = compose_filename(filename, run->rootname);

    *do_twiss_output = output_at_each_step;

    if (!matched) {
        if (beta_x<=0 || beta_y<=0)
            bomb("invalid initial beta-functions given in twiss_output namelist", NULL);
        }

    if (filename) {
        SDDS_ElegantOutputSetup(&SDDS_twiss, filename, SDDS_BINARY, 1, "Twiss parameters",
                                run->runfile, run->lattice, parameter_definition, 
                                (radiation_integrals?N_PARAMETERS:IP_ALPHAC+1),
                                column_definition, N_COLUMNS, "setup_twiss_output",
                                SDDS_EOS_NEWFILE|SDDS_EOS_COMPLETE);
        SDDS_twiss_initialized = 1;
        twiss_count = 0;
        }
    else
        SDDS_twiss_initialized = 0;
    twiss_initialized = 1;

    beamline->flags |= BEAMLINE_TWISS_WANTED;
    if (radiation_integrals)
      beamline->flags |= BEAMLINE_RADINT_WANTED;
    
    log_exit("setup_twiss_output");
    }

void finish_twiss_output(void)
{
    log_entry("finish_twiss_output");
    if (SDDS_twiss_initialized && !SDDS_Terminate(&SDDS_twiss)) {
        SDDS_SetError("Problem terminating SDDS output (finish_twiss_output)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
    SDDS_twiss_initialized = twiss_count = 0;
    log_exit("finish_twiss_output");
    }

void run_twiss_output(RUN *run, LINE_LIST *beamline, double *starting_coord, long tune_corrected)
{
    ELEMENT_LIST *eptr, *elast;
    long n_elem, last_n_elem;
    TWISS twiss_ave, twiss_min, twiss_max;

/*
    if (beamline->flags&BEAMLINE_TWISS_CURRENT)
        return;
 */

    log_entry("run_twiss_output");

#ifdef DEBUG
    fprintf(stderr, "now in run_twiss_output\n");
#endif

    if (tune_corrected==0 && !output_before_tune_correction) {
        log_exit("run_twiss_output");
        return;
        }

    eptr = beamline->elem_twiss = &(beamline->elem);
    n_elem = last_n_elem = beamline->n_elems;
    while (eptr) {
        if (eptr->type==T_RECIRC) {
            last_n_elem = n_elem;
            beamline->elem_twiss = beamline->elem_recirc = eptr;
            }
        eptr = eptr->succ;
        n_elem --;
        }
    n_elem = last_n_elem;
    
    compute_twiss_parameters(run, beamline, starting_coord, matched, radiation_integrals,
                             beta_x, alpha_x, eta_x, etap_x,
                             beta_y, alpha_y, eta_y, etap_y);
    elast = beamline->elast;

    if (run->default_order>=2) {
#ifdef DEBUG
        fprintf(stderr, "computing chromaticities\n");
#endif
        fprintf(stderr, "%s Twiss parameters (chromaticity valid for fully second-order calculation only!):\n",
            matched?"periodic":"final");
        fprintf(stderr, "         beta          alpha           nu           eta          eta'       dnu/d(dp/p)     accept.\n");
        fprintf(stderr, "          m                          1/2pi           m                         1/2pi        mm-mrad\n");
        fprintf(stderr, "------------------------------------------------------------------------------------------------------\n");
        fprintf(stderr, "  x: %13.6le %13.6le %13.6le %13.6le %13.6le %13.6le %13.6le\n",  
               elast->twiss->betax, elast->twiss->alphax, beamline->tune[0], elast->twiss->etax, elast->twiss->etapx,
               beamline->chromaticity[0], 1e6*beamline->acceptance[0]);
        if (statistics) {
            compute_twiss_statistics(beamline, &twiss_ave, &twiss_min, &twiss_max);
            fprintf(stderr, "ave: %13.6le %13.6le %-13s %13.6le %13.6le\n",
                   twiss_ave.betax, twiss_ave.alphax, "", twiss_ave.etax, twiss_ave.etapx);
            fprintf(stderr, "min: %13.6le %13.6le %-13s %13.6le %13.6le\n",
                   twiss_min.betax, twiss_min.alphax, "", twiss_min.etax, twiss_min.etapx);
            fprintf(stderr, "max: %13.6le %13.6le %-13s %13.6le %13.6le\n",
                   twiss_max.betax, twiss_max.alphax, "", twiss_max.etax, twiss_max.etapx);
            }
        fprintf(stderr, "  y: %13.6le %13.6le %13.6le %13.6le %13.6le %13.6le %13.6le\n",  
               elast->twiss->betay, elast->twiss->alphay, beamline->tune[1], elast->twiss->etay, elast->twiss->etapy,
               beamline->chromaticity[1], 1e6*beamline->acceptance[1]);
        if (statistics) {
            fprintf(stderr, "ave: %13.6le %13.6le %-13s %13.6le %13.6le\n",
                   twiss_ave.betay, twiss_ave.alphay, "", twiss_ave.etay, twiss_ave.etapy);
            fprintf(stderr, "min: %13.6le %13.6le %-13s %13.6le %13.6le\n",
                   twiss_min.betay, twiss_min.alphay, "", twiss_min.etay, twiss_min.etapy);
            fprintf(stderr, "max: %13.6le %13.6le %-13s %13.6le %13.6le\n",
                   twiss_max.betay, twiss_max.alphay, "", twiss_max.etay, twiss_max.etapy);
            }
        }
    else {
        fprintf(stderr, "%s Twiss parameters:\n", matched?"periodic":"final");
        fprintf(stderr, "         beta          alpha           nu           eta          eta'        accept.\n");
        fprintf(stderr, "          m                          1/2pi           m                       mm-mrad\n");
        fprintf(stderr, "---------------------------------------------------------------------------------------\n");
        fprintf(stderr, "  x: %13.6le %13.6le %13.6le %13.6le %13.6le %13.6le\n",  
               elast->twiss->betax, elast->twiss->alphax, beamline->tune[0], elast->twiss->etax, elast->twiss->etapx,
               1e6*beamline->acceptance[0]);
        if (statistics) {
            compute_twiss_statistics(beamline, &twiss_ave, &twiss_min, &twiss_max);
            fprintf(stderr, "ave: %13.6le %13.6le %-13s %13.6le %13.6le\n",
                   twiss_ave.betax, twiss_ave.alphax, "", twiss_ave.etax, twiss_ave.etapx);
            fprintf(stderr, "min: %13.6le %13.6le %-13s %13.6le %13.6le\n",
                   twiss_min.betax, twiss_min.alphax, "", twiss_min.etax, twiss_min.etapx);
            fprintf(stderr, "max: %13.6le %13.6le %-13s %13.6le %13.6le\n",
                   twiss_max.betax, twiss_max.alphax, "", twiss_max.etax, twiss_max.etapx);
            }
        fprintf(stderr, "  y: %13.6le %13.6le %13.6le %13.6le %13.6le %13.6le\n",  
               elast->twiss->betay, elast->twiss->alphay, beamline->tune[1], elast->twiss->etay, elast->twiss->etapy,
               1e6*beamline->acceptance[1]);
        if (statistics) {
            fprintf(stderr, "ave: %13.6le %13.6le %-13s %13.6le %13.6le\n",
                   twiss_ave.betay, twiss_ave.alphay, "", twiss_ave.etay, twiss_ave.etapy);
            fprintf(stderr, "min: %13.6le %13.6le %-13s %13.6le %13.6le\n",
                   twiss_min.betay, twiss_min.alphay, "", twiss_min.etay, twiss_min.etapy);
            fprintf(stderr, "max: %13.6le %13.6le %-13s %13.6le %13.6le\n",
                   twiss_max.betay, twiss_max.alphay, "", twiss_max.etay, twiss_max.etapy);
            }
        }

    if (beamline->acc_limit_name[0])
        fprintf(stderr, "x acceptance limited by %s ending at %e m\n", beamline->acc_limit_name[0], beamline->acceptance[2]);
    if (beamline->acc_limit_name[1])
        fprintf(stderr, "y acceptance limited by %s ending at %e m\n", beamline->acc_limit_name[1], beamline->acceptance[3]);

    if (SDDS_twiss_initialized) {
      double alphac = 0;
      if (beamline->matrix->C[4]!=0) {
        alphac = (beamline->matrix->R[4][5] + 
                  beamline->matrix->R[4][0]*elast->twiss->etax +
                  beamline->matrix->R[4][1]*elast->twiss->etapx +
                  beamline->matrix->R[4][2]*elast->twiss->etay +
                  beamline->matrix->R[4][3]*elast->twiss->etapy)/beamline->matrix->C[4];
      }
      dump_twiss_parameters(beamline->twiss0, beamline->elem_twiss, n_elem,
                            beamline->tune, 
                            radiation_integrals?&(beamline->radIntegrals):NULL, 
                            beamline->chromaticity, beamline->acceptance, 
                            alphac,
                            final_values_only, tune_corrected, run);
    }

    log_exit("run_twiss_output");
    }

void compute_twiss_parameters(RUN *run, LINE_LIST *beamline, double *starting_coord, long periodic,
                              long radiation_integrals,
                              double betax, double alphax, double etax, double etapx, 
                              double betay, double alphay, double etay, double etapy)
{
    VMATRIX *M;
    double chromx, chromy;
    double x_acc_z, y_acc_z;
    ELEMENT_LIST *eptr, *elast;
    char *x_acc_name, *y_acc_name;
    long i;
    
    log_entry("compute_twiss_parameters");

    if (!beamline->twiss0)
        beamline->twiss0 = tmalloc(sizeof(*beamline->twiss0));

    eptr = beamline->elem_twiss = &(beamline->elem);

    elast = eptr;
    while (eptr) {
        if (eptr->type==T_RECIRC)
            beamline->elem_twiss = beamline->elem_recirc = eptr;
        elast = eptr;
        eptr = eptr->succ;
        }
    beamline->elast = elast;

    if (periodic) {
        beamline->matrix = compute_periodic_twiss(&betax, &alphax, &etax, &etapx, beamline->tune,
            &betay, &alphay, &etay, &etapy, beamline->tune+1, beamline->elem_twiss, starting_coord, run);
#ifdef DEBUG
        fprintf(stderr, "matched parameters computed--returned to compute_twiss_parameters\n");
        fprintf(stderr, "beamline matrix has order %ld\n", beamline->matrix->order);
        fprintf(stderr, "beamline matrix pointers:  %x, %x, %x, %x\n",
               beamline->matrix, beamline->matrix->C, beamline->matrix->R, beamline->matrix->T);
#endif
        if (run->default_order>=2 && !(beamline->matrix->T))
            bomb("logic error: T matrix is NULL on return from compute_periodic_twiss", NULL);
        }
    else {
        VMATRIX *M1;
        if (run->default_order>1 && starting_coord) {
            M1 = tmalloc(sizeof(*M1));
            initialize_matrices(M1, 1);
            M1->C[0] = starting_coord[0];
            M1->C[1] = starting_coord[0];
            M1->C[2] = starting_coord[0];
            M1->C[3] = starting_coord[0];
            M1->C[4] = starting_coord[0];
            M1->C[5] = starting_coord[0];
            for (i=0; i<6; i++)
                M1->R[i][i] = 1;
            fill_in_matrices(beamline->elem_twiss, run);
            beamline->matrix = append_full_matrix(beamline->elem_twiss, run, M1, run->default_order);
            }
        else
            beamline->matrix = full_matrix(beamline->elem_twiss, run, run->default_order);
        }

    beamline->twiss0->betax  = betax;
    beamline->twiss0->alphax = alphax;
    beamline->twiss0->phix   = 0;
    beamline->twiss0->etax   = etax;
    beamline->twiss0->etapx  = etapx;
    beamline->twiss0->betay  = betay;
    beamline->twiss0->alphay = alphay;
    beamline->twiss0->phiy   = 0;
    beamline->twiss0->etay   = etay;
    beamline->twiss0->etapy  = etapy;

    if (run->default_order>=2 && !(beamline->matrix->T))
        bomb("logic error: beamline T matrix is NULL in compute_twiss_parameters", NULL);

#ifdef DEBUG
    fprintf(stderr, "propagating parameters\n");
    fprintf(stderr, "beamline matrix pointers:  %x, %x, %x, %x\n",
           beamline->matrix, beamline->matrix->C, beamline->matrix->R, beamline->matrix->T);
#endif

    propagate_twiss_parameters(beamline->twiss0, beamline->tune, 
                               (radiation_integrals?&(beamline->radIntegrals):NULL),
                               beamline->elem_twiss, 0, run, starting_coord);
    propagate_twiss_parameters(beamline->twiss0, beamline->tune+1, NULL,
                               beamline->elem_twiss, 1, run, starting_coord);
    
    if (radiation_integrals) {
      double Rce, gamma;
      RADIATION_INTEGRALS *RI;
      RI = &(beamline->radIntegrals);
      gamma = sqrt(sqr(run->p_central)+1);
      Rce = sqr(e_mks)/(1e7*me_mks);
      RI->Uo = me_mev*Rce*RI->I[1]*2./3.*ipow(gamma,4);
      RI->Jx = 1 - RI->I[3]/RI->I[1];
      RI->Jdelta = 3 - RI->Jx;
      RI->Jy = 1;
      RI->tauy = 1./(Rce/3*ipow(gamma,3)*c_mks/beamline->revolution_length*RI->I[1]);
      RI->taux = RI->tauy*RI->Jy/RI->Jx;
      RI->taudelta = RI->tauy*RI->Jy/RI->Jdelta;
      RI->sigmadelta = gamma*sqrt(55./32./sqrt(3.)*hbar_mks/(me_mks*c_mks)*RI->I[2]/(2*RI->I[1]+RI->I[3]));
      RI->ex0 = sqr(gamma)*55./32./sqrt(3.)*hbar_mks/(me_mks*c_mks)*RI->I[4]/(RI->I[1]-RI->I[3]);
    }
    
#ifdef DEBUG
    fprintf(stderr, "finding acceptance\n");
#endif
    beamline->acceptance[0] = find_acceptance(beamline->elem_twiss, 0, run, &x_acc_name, &x_acc_z);
    beamline->acceptance[1] = find_acceptance(beamline->elem_twiss, 1, run, &y_acc_name, &y_acc_z);
    beamline->acceptance[2] = x_acc_z;
    beamline->acceptance[3] = y_acc_z;
    if (x_acc_name)
        cp_str(&beamline->acc_limit_name[0], x_acc_name);
    else
        beamline->acc_limit_name[0] = NULL;
    if (y_acc_name)
        cp_str(&beamline->acc_limit_name[1], y_acc_name);
    else
        beamline->acc_limit_name[1] = NULL;

    beamline->twiss0->apx = beamline->elem_twiss->twiss->apx;
    beamline->twiss0->apy = beamline->elem_twiss->twiss->apy;

#ifdef DEBUG
    fprintf(stderr, "computing chromaticities\n");
#endif

    chromx = chromy = 0;
    if (periodic) {
        if (!(M = beamline->matrix))
            bomb("logic error: revolution matrix is NULL in compute_twiss_parameters", NULL);

        if (run->default_order>=2) {
#ifdef DEBUG
            fprintf(stderr, "computing chromaticities\n");
#endif
            if (!(M->T))
                bomb("logic error: T matrix is NULL in compute_twiss_parameters", NULL);
            computeChromaticities(&chromx, &chromy, beamline->twiss0, M);
            }
        }
    else {
        M = beamline->matrix;
        elast = beamline->elast;

        if (!elast)
            bomb("logic error in compute_twiss_parameters--elast pointer is NULL", NULL);
        if (!elast->twiss)
            bomb("logic error in compute_twiss_parameters--elast->twiss pointer is NULL", NULL);

        betax = elast->twiss->betax;
        alphax = elast->twiss->alphax;
        etax = elast->twiss->etax;
        etapx = elast->twiss->etapx;
        betay = elast->twiss->betay;
        alphay = elast->twiss->alphay;
        etay = elast->twiss->etay;
        etapy = elast->twiss->etapy;

        if (run->default_order>=2) {
            if (!(M->T))
                bomb("logic error: T matrix is NULL in compute_twiss_parameters", NULL);
            computeChromaticities(&chromx, &chromy, beamline->twiss0, M);
#ifdef DEBUG
            fprintf(stderr, "chomaticities: %e, %e\n", chromx, chromy);
#endif
            }
        }
    beamline->chromaticity[0] = chromx;
    beamline->chromaticity[1] = chromy;

    beamline->flags |= BEAMLINE_TWISS_DONE+BEAMLINE_TWISS_CURRENT;
    if (radiation_integrals) 
      beamline->flags |= BEAMLINE_RADINT_DONE+BEAMLINE_RADINT_CURRENT;
    
    log_exit("compute_twiss_parameters");
    }

void update_twiss_parameters(RUN *run, LINE_LIST *beamline)
{
    log_entry("update_twiss_parameters");
    compute_twiss_parameters(run, beamline, 
                             beamline->closed_orbit?beamline->closed_orbit->centroid:NULL, matched, 
                             radiation_integrals,
                             beta_x, alpha_x, eta_x, etap_x, beta_y, alpha_y, eta_y, etap_y);
    log_exit("update_twiss_parameters");
    }

void copy_doubles(double *target, double *source, long n)
{
    log_entry("copy_doubles");
    while (n--)
        *target++ = *source++;
    log_exit("copy_doubles");
    }


long has_aperture(ELEMENT_LIST *elem);

double find_acceptance(
        ELEMENT_LIST *elem, long plane, RUN *run, char **name, double *z
        )
{
    double beta, acceptance, tmp, last_z, *last_aperture;
    double tube_aperture, aperture, centroid, room;
    double other_centroid, a, b, a_tube, b_tube;
    SCRAPER *scraper;
    ELEMENT_LIST *ap_elem;

    log_entry("find_acceptance");

    acceptance = tube_aperture = a_tube = b_tube = 0;
    last_z = 0;
    last_aperture = NULL;
    *name = NULL;
    *z = 0;
    while (elem) {
        beta = *(((double*)elem->twiss) + (plane?TWISS_Y_OFFSET:0));
        centroid = *(((double*)elem->twiss) + (plane?1:0) + TWISS_CENT_OFFSET);
        other_centroid = *(((double*)elem->twiss) + (plane?0:1) + TWISS_CENT_OFFSET);
        aperture = 0;
        if (!has_aperture(elem) && has_aperture(elem->succ))
            ap_elem = elem->succ;
        else
            ap_elem = elem;
        switch (ap_elem->type) {
          case T_MAXAMP:
            if (((MAXAMP*)ap_elem->p_elem)->elliptical==0) {
                if (plane)
                    tube_aperture = aperture = ((MAXAMP*)ap_elem->p_elem)->y_max;
                else
                    tube_aperture = aperture = ((MAXAMP*)ap_elem->p_elem)->x_max;
                a_tube = b_tube = 0;
                }
            else {
                if (plane) {
                    a_tube = ((MAXAMP*)ap_elem->p_elem)->y_max;
                    b_tube = ((MAXAMP*)ap_elem->p_elem)->x_max;
                    }
                else {
                    a_tube = ((MAXAMP*)ap_elem->p_elem)->x_max;
                    b_tube = ((MAXAMP*)ap_elem->p_elem)->y_max;
                    }
                if ((aperture = sqr(a_tube)*(1 - sqr(other_centroid/b_tube)))<0)
                    aperture = 0;
                else
                    aperture = sqrt(aperture);
                }
            break;
          case T_RCOL:
            if (plane) {
                centroid -= ((RCOL*)ap_elem->p_elem)->dy;
                other_centroid -= ((RCOL*)ap_elem->p_elem)->dx;
                }
            else {
                centroid -= ((RCOL*)ap_elem->p_elem)->dx;
                other_centroid -= ((RCOL*)ap_elem->p_elem)->dy;
                }
            if (plane)
                aperture = ((RCOL*)ap_elem->p_elem)->y_max;
            else
                aperture = ((RCOL*)ap_elem->p_elem)->x_max;
            break;
          case T_ECOL:
            if (plane) {
                centroid -= ((ECOL*)ap_elem->p_elem)->dy;
                other_centroid -= ((ECOL*)ap_elem->p_elem)->dx;
                }
            else {
                centroid -= ((ECOL*)ap_elem->p_elem)->dx;
                other_centroid -= ((ECOL*)ap_elem->p_elem)->dy;
                }
            if (plane) {
                a = ((ECOL*)ap_elem->p_elem)->y_max;
                b = ((ECOL*)ap_elem->p_elem)->x_max;
                }
            else {
                a = ((ECOL*)ap_elem->p_elem)->x_max;
                b = ((ECOL*)ap_elem->p_elem)->y_max;
                }
            if ((aperture = sqr(a)*(1 - sqr(other_centroid/b)))<0)
                aperture = 1/HUGE;
            else
                aperture = sqrt(aperture);
            break;
          case T_SCRAPER:
            scraper = (SCRAPER*)ap_elem->p_elem;
            if (plane) {
                centroid -= scraper->dy;
                other_centroid -= scraper->dx;
                }
            else {
                centroid -= scraper->dx;
                other_centroid -= scraper->dy;
                }
            if (plane) {
                if (scraper->direction==1 || scraper->direction==3)
                    aperture = fabs(scraper->position);
                }
            else {
                if (scraper->direction==0 || scraper->direction==2)
                    aperture = fabs(scraper->position);
                }
            break;
          default:
            if (a_tube && b_tube) {
                if ((aperture = sqr(a_tube)*(1-sqr(other_centroid/b_tube)))<0)
                    aperture = 1/HUGE;
                else
                    aperture = sqrt(aperture);
                }
            else {
                aperture = tube_aperture;
                }
            break;
            }
        if (aperture) {
            if ((room=aperture-fabs(centroid))>0) {
                if (((tmp=sqr(room)/beta)<acceptance || !acceptance)) {
                    *name = elem->name;
                    *z = elem->end_pos;
                    acceptance = tmp;
                    }
                }
            else {
                *name = elem->name;
                *z = elem->end_pos;
                acceptance = 0;
                break;
                }
            }
        
        if (last_aperture) {
            if (last_z==elem->end_pos) {
                if (aperture==0)
                    aperture = *last_aperture;
                else if (*last_aperture==0)
                    *last_aperture = aperture;
                }
            }
        
        *(last_aperture = (((double*)elem->twiss) + (plane?TWISS_Y_OFFSET:0) + 5)) = aperture;
        last_z = elem->end_pos;
        elem = elem->succ;
        }
    log_exit("find_acceptance");
    return(acceptance);
    }

long has_aperture(ELEMENT_LIST *elem)
{
    if (!elem)
        return(0);
    switch (elem->type) {
      case T_MAXAMP: case T_RCOL: case T_ECOL: case T_SCRAPER:
        return(1);
        break;
      default:
        return(0);
        break;
        }
    }

void modify_rfca_matrices(ELEMENT_LIST *eptr, long order)
/* Replace the matrices for rf cavities with drift matrices. 
   Used prior to twiss parameter computations. */
{
    while (eptr) {
        if (eptr->type==T_RFCA || eptr->type==T_MODRF) {
            if (eptr->matrix) {
                free_matrices(eptr->matrix);
                tfree(eptr->matrix);
                }
            if (eptr->type==T_RFCA)
                eptr->matrix = drift_matrix(((RFCA*)eptr->p_elem)->length, order);
            else
                eptr->matrix = drift_matrix(((MODRF*)eptr->p_elem)->length, order);
            }
        eptr = eptr->succ;
        }
    }

#define ASSIGN_MINMAX(min, max, value) ((min>value?min=value:1),(max<value?max=value:1))

void compute_twiss_statistics(LINE_LIST *beamline, TWISS *twiss_ave, TWISS *twiss_min, TWISS *twiss_max)
{
    ELEMENT_LIST *eptr;
    double dz, end_pos;

    if (!twiss_ave) {
        fprintf(stderr, "error: NULL twiss_ave pointer in compute_twiss_statistics\n");
        abort();
        }
    if (!twiss_min) {
        fprintf(stderr, "error: NULL twiss_min pointer in compute_twiss_statistics\n");
        abort();
        }
    if (!twiss_max) {
        fprintf(stderr, "error: NULL twiss_max pointer in compute_twiss_statistics\n");
        abort();
        }
        
    zero_memory(twiss_ave, sizeof(*twiss_ave));
    twiss_min->betax = twiss_min->alphax = twiss_min->etax = twiss_min->etapx = HUGE;
    twiss_min->betay = twiss_min->alphay = twiss_min->etay = twiss_min->etapy = HUGE;
    twiss_max->betax = twiss_max->alphax = twiss_max->etax = twiss_max->etapx = -HUGE;
    twiss_max->betay = twiss_max->alphay = twiss_max->etay = twiss_max->etapy = -HUGE;

    eptr = beamline->elem_twiss;
    while (eptr) {
        ASSIGN_MINMAX(twiss_min->betax, twiss_max->betax, eptr->twiss->betax);
        ASSIGN_MINMAX(twiss_min->alphax, twiss_max->alphax, eptr->twiss->alphax);
        ASSIGN_MINMAX(twiss_min->etax, twiss_max->etax, eptr->twiss->etax);
        ASSIGN_MINMAX(twiss_min->etapx, twiss_max->etapx, eptr->twiss->etapx);
        ASSIGN_MINMAX(twiss_min->betay, twiss_max->betay, eptr->twiss->betay);
        ASSIGN_MINMAX(twiss_min->alphay, twiss_max->alphay, eptr->twiss->alphay);
        ASSIGN_MINMAX(twiss_min->etay, twiss_max->etay, eptr->twiss->etay);
        ASSIGN_MINMAX(twiss_min->etapy, twiss_max->etapy, eptr->twiss->etapy);
        if (eptr->pred) {
           dz = (eptr->end_pos - eptr->pred->end_pos)/2;
           twiss_ave->betax += (eptr->pred->twiss->betax + eptr->twiss->betax)*dz;
           twiss_ave->alphax += (eptr->pred->twiss->alphax + eptr->twiss->alphax)*dz;
           twiss_ave->etax += (eptr->pred->twiss->etax + eptr->twiss->etax)*dz;
           twiss_ave->etapx += (eptr->pred->twiss->etapx + eptr->twiss->etapx)*dz;
           twiss_ave->betay += (eptr->pred->twiss->betay + eptr->twiss->betay)*dz;
           twiss_ave->alphay += (eptr->pred->twiss->alphay + eptr->twiss->alphay)*dz;
           twiss_ave->etay += (eptr->pred->twiss->etay + eptr->twiss->etay)*dz;
           twiss_ave->etapy += (eptr->pred->twiss->etapy + eptr->twiss->etapy)*dz;
           }
        end_pos = eptr->end_pos;
        eptr = eptr->succ;
        }
    if (end_pos) {
        twiss_ave->betax  /= end_pos;
        twiss_ave->alphax /= end_pos;
        twiss_ave->etax   /= end_pos;
        twiss_ave->etapx  /= end_pos;
        twiss_ave->betay  /= end_pos;
        twiss_ave->alphay /= end_pos;
        twiss_ave->etay   /= end_pos;
        twiss_ave->etapy  /= end_pos;
        }
    }
