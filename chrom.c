/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: chrom.c
 * purpose: chromaticity correction by adjustment of sextupoles.
 *
 * Michael Borland, 1992
 */
#if defined(SOLARIS) && !defined(__GNUC__)
#include <sunmath.h>
#endif
#include "mdb.h"
#include "track.h"
#include "chromDefs.h"

static FILE *fp_sl = NULL;
static long alter_defined_values;

void setup_chromaticity_correction(NAMELIST_TEXT *nltext, RUN *run, LINE_LIST *beamline, CHROM_CORRECTION *chrom)
{
    VMATRIX *M;
    ELEMENT_LIST *eptr, *elast;
#include "chrom.h"
    unsigned long unstable;
    
    log_entry("setup_chromaticity_correction");

    cp_str(&sextupoles, "sf sd");

    /* process namelist input */
    set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
    set_print_namelist_flags(0);
    process_namelist(&chromaticity, nltext);
    str_toupper(sextupoles);
    print_namelist(stdout, &chromaticity);

    if (run->default_order<2)
        bomb("default order must be >= 2 for chromaticity correction", NULL);

    if (chrom->name)
        tfree(chrom->name);
    chrom->name = tmalloc(sizeof(*chrom->name)*(chrom->n_families=1));
    while ((chrom->name[chrom->n_families-1]=get_token(sextupoles)))
        chrom->name = trealloc(chrom->name, sizeof(*chrom->name)*(chrom->n_families+=1));
    if ((--chrom->n_families)<1)
        bomb("too few sextupoles given for chromaticity correction", NULL);
    chrom->chromx = dnux_dp;
    chrom->chromy = dnuy_dp;
    chrom->n_iterations = n_iterations;
    chrom->correction_fraction = correction_fraction;
    alter_defined_values = change_defined_values;
    chrom->strengthLimit = strength_limit;
    chrom->use_perturbed_matrix = use_perturbed_matrix;
    chrom->sextupole_tweek = sextupole_tweek;
    chrom->tolerance = tolerance;
    
    if (!use_perturbed_matrix) {
      if (!beamline->twiss0 || !beamline->matrix) {
        double beta_x, alpha_x, eta_x, etap_x;
        double beta_y, alpha_y, eta_y, etap_y;

        fprintf(stdout, "Computing periodic Twiss parameters.\n");
        fflush(stdout);

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
        
        M = beamline->matrix = compute_periodic_twiss(&beta_x, &alpha_x, &eta_x, &etap_x, beamline->tune,
                                                      &beta_y, &alpha_y, &eta_y, &etap_y, beamline->tune+1, 
                                                      beamline->elem_twiss, NULL, run, 
                                                      &unstable, NULL, NULL);
        beamline->twiss0->betax  = beta_x;
        beamline->twiss0->alphax = alpha_x;
        beamline->twiss0->phix   = 0;
        beamline->twiss0->etax   = eta_x;
        beamline->twiss0->etapx  = etap_x;
        beamline->twiss0->betay  = beta_y;
        beamline->twiss0->alphay = alpha_y;
        beamline->twiss0->phiy   = 0;
        beamline->twiss0->etay   = eta_y;
        beamline->twiss0->etapy  = etap_y;
        
        propagate_twiss_parameters(beamline->twiss0, beamline->tune, beamline->waists,
                                   NULL, beamline->elem_twiss, run, NULL);
      }

      if (!(M=beamline->matrix) || !M->C || !M->R || !M->T)
        bomb("something wrong with transfer map for beamline (setup_chromaticity_correction)", NULL);

      computeChromCorrectionMatrix(run, beamline, chrom);
    }
    

    if (strength_log) {
        strength_log = compose_filename(strength_log, run->rootname);
        fp_sl = fopen_e(strength_log, "w", 0);
        fprintf(fp_sl, "SDDS1\n&column name=Step, type=long, description=\"Simulation step\" &end\n");
        fprintf(fp_sl, "&column name=K2, type=double, units=\"1/m$a2$n\" &end\n");
        fprintf(fp_sl, "&column name=SextupoleName, type=string  &end\n");
        fprintf(fp_sl, "&data mode=ascii, no_row_counts=1 &end\n");
        fflush(fp_sl);
        }

    log_exit("setup_chromaticity_correction");
    
}

void computeChromCorrectionMatrix(RUN *run, LINE_LIST *beamline, CHROM_CORRECTION *chrom)
{    
    VMATRIX *M;
    double chromx, chromy;
    double chromx0, chromy0;
    double K2=0.0, *K2ptr;
    ELEMENT_LIST *context;
    long i, count, K2_param=0;
    MATRIX *C, *Ct, *CtC, *inv_CtC;

    m_alloc(&C, 2, chrom->n_families);
    m_alloc(&Ct, chrom->n_families, 2);
    m_alloc(&CtC, chrom->n_families, chrom->n_families);
    m_alloc(&inv_CtC, chrom->n_families, chrom->n_families);
    m_alloc(&(chrom->T), chrom->n_families, 2);
    m_alloc(&(chrom->dK2), chrom->n_families, 1);
    m_alloc(&(chrom->dchrom), 2, 1);

    fprintf(stdout, "Computing chromaticity influence matrix for all named sextupoles.\n");
    fflush(stdout);
    computeChromaticities(&chromx0, &chromy0, 
                          NULL, NULL, NULL, NULL, 
			  beamline->twiss0, M=beamline->matrix);

    for (i=0; i<chrom->n_families; i++) {
        count = 0;
        context = NULL;
        while ((context=find_element(chrom->name[i], &context, beamline->elem_twiss))) {
            if (!count && !(K2_param=confirm_parameter("K2", context->type))) {
                fprintf(stdout, "error: element %s does not have K2 parameter\n", 
                        context->name);
                fflush(stdout);
                exit(1);
                }
            if (!(K2ptr = (double*)(context->p_elem + entity_description[context->type].parameter[K2_param].offset)))
                bomb("K2ptr NULL in setup_chromaticity_correction", NULL);
            if (count==0)
                K2 = *K2ptr;
            *K2ptr = K2 + chrom->sextupole_tweek;
            if (context->matrix)
                free_matrices(context->matrix);
            compute_matrix(context, run, NULL);
            count++;
            }
        if (count==0) {
            fprintf(stdout, "error: element %s is not in the beamline.\n",  chrom->name[i]);
            fflush(stdout);
            exit(1);
            }
        if (beamline->links)
            assert_element_links(beamline->links, run, beamline, STATIC_LINK+DYNAMIC_LINK);
        M = full_matrix(beamline->elem_twiss, run, 2);
        computeChromaticities(&chromx, &chromy, 
                              NULL, NULL, NULL, NULL, beamline->twiss0, M);

        C->a[0][i] = (chromx-chromx0)/chrom->sextupole_tweek;
        C->a[1][i] = (chromy-chromy0)/chrom->sextupole_tweek;
        if (C->a[0][i]==0 || C->a[1][i]==0) {
            fprintf(stdout, "error: element %s does not change the chromaticity!\n", chrom->name[i]);
            fflush(stdout);
            exit(1);
            }
        count = 0;
        context = NULL;
        while ((context=find_element(chrom->name[i], &context, beamline->elem_twiss))) {
            if (!count && !(K2_param=confirm_parameter("K2", context->type))) {
                fprintf(stdout, "error: element %s does not have K2 parameter\n", 
                        context->name);
                fflush(stdout);
                exit(1);
                }
            if (!(K2ptr = (double*)(context->p_elem + entity_description[context->type].parameter[K2_param].offset)))
                bomb("K2ptr NULL in setup_chromaticity_correction", NULL);
            if (!K2ptr)
                bomb("K2ptr NULL in setup_chromaticity_correction", NULL);
            *K2ptr = K2;
            if (context->matrix)
                free_matrices(context->matrix);
            compute_matrix(context, run, NULL);
            count++;
            }
        }
    beamline->matrix = full_matrix(beamline->elem_twiss, run, run->default_order);

    fprintf(stdout, "\nfamily           dCHROMx/dK2        dCHROMy/dK2\n");
    fflush(stdout);
    for (i=0; i<chrom->n_families; i++)
       fprintf(stdout, "%10s:    %14.7e     %14.7e\n", chrom->name[i], C->a[0][i], C->a[1][i]);
       fflush(stdout);

    m_trans(Ct, C);
    m_mult(CtC, Ct, C);
    m_invert(inv_CtC, CtC);
    m_mult(chrom->T, inv_CtC, Ct);

    fprintf(stdout, "\nfamily           dK2/dCHROMx        dK2/dCHROMy\n");
    fflush(stdout);
    for (i=0; i<chrom->n_families; i++)
       fprintf(stdout, "%10s:    %14.7e     %14.7e\n", chrom->name[i], chrom->T->a[i][0], chrom->T->a[i][1]);
       fflush(stdout);
    fprintf(stdout, "\n");
    fflush(stdout);

    m_free(&C);
    m_free(&Ct);
    m_free(&CtC);
    m_free(&inv_CtC);
  }



long do_chromaticity_correction(CHROM_CORRECTION *chrom, RUN *run, LINE_LIST *beamline,
            double *clorb, long step, long last_iteration)
{
    VMATRIX *M;
    double chromx0, chromy0, dchromx, dchromy;
    double K2=0.0, *K2ptr;
    ELEMENT_LIST *context;
    long i, K2_param=0, type=0, iter, count;
    double beta_x, alpha_x, eta_x, etap_x;
    double beta_y, alpha_y, eta_y, etap_y;
    unsigned long unstable;
    double lastError, presentError;
    
    log_entry("do_chromaticity_correction");

    if (!beamline->matrix || !beamline->twiss0) {
        if (!beamline->twiss0)
            beamline->twiss0 = tmalloc(sizeof(*beamline->twiss0));
        if (!beamline->elem_twiss) {
            ELEMENT_LIST *eptr;
            eptr = beamline->elem_twiss = &(beamline->elem);
            while (eptr) {
                if (eptr->type==T_RECIRC)
                    beamline->elem_twiss = beamline->elem_recirc = eptr;
                eptr = eptr->succ;
                }
            }
    
        M = beamline->matrix = compute_periodic_twiss(&beta_x, &alpha_x, &eta_x, &etap_x, beamline->tune,
                                                      &beta_y, &alpha_y, &eta_y, &etap_y, beamline->tune+1, 
                                                      beamline->elem_twiss, clorb, run, &unstable, NULL, NULL);
        beamline->twiss0->betax  = beta_x;
        beamline->twiss0->alphax = alpha_x;
        beamline->twiss0->phix   = 0;
        beamline->twiss0->etax   = eta_x;
        beamline->twiss0->etapx  = etap_x;
        beamline->twiss0->betay  = beta_y;
        beamline->twiss0->alphay = alpha_y;
        beamline->twiss0->phiy   = 0;
        beamline->twiss0->etay   = eta_y;
        beamline->twiss0->etapy  = etap_y;
        
        propagate_twiss_parameters(beamline->twiss0, beamline->tune, beamline->waists,
                                   NULL, beamline->elem_twiss, run, clorb);
        }
    else if (beamline->matrix->order<2)
        beamline->matrix = full_matrix(beamline->elem_twiss, run, 2);

    if (!(M = beamline->matrix) || !M->C || !M->R || !M->T)
        bomb("something wrong with transfer map for beamline (do_chromaticity_correction.1)", NULL);

    computeChromaticities(&chromx0, &chromy0, 
			  NULL, NULL, NULL, NULL, beamline->twiss0, M);

    fprintf(stdout, "\nAdjusting chromaticities:\n");
    fflush(stdout);
    fprintf(stdout, "initial chromaticities:  %e  %e\n", chromx0, chromy0);
    fflush(stdout);

    presentError = DBL_MAX;
    for (iter=0; iter<chrom->n_iterations; iter++) {
        dchromx = chrom->chromx - chromx0;
        dchromy = chrom->chromy - chromy0;
        if (chrom->tolerance>0 &&
            chrom->tolerance>fabs(dchromx) &&
            chrom->tolerance>fabs(dchromy) )
          break;

        lastError = presentError;
        presentError = sqr(dchromx)+sqr(dchromy);
        if (iter && presentError>lastError) {
          fprintf(stdout, "Error increasing---reducing gain\n");
          fflush(stdout);
          chrom->correction_fraction /= 10;
        }

        if (chrom->use_perturbed_matrix)
          computeChromCorrectionMatrix(run, beamline, chrom);
        chrom->dchrom->a[0][0] = dchromx;
        chrom->dchrom->a[1][0] = dchromy;

        m_mult(chrom->dK2, chrom->T, chrom->dchrom);
        for (i=0; i<chrom->n_families; i++) {
          if (isnan(chrom->correction_fraction*chrom->dK2->a[i][0]) ||
              isinf(chrom->correction_fraction*chrom->dK2->a[i][0]))
            break;
        }
        if (i!=chrom->n_families) {
          fprintf(stdout, "Unable to correct chromaticity---diverged\n");
          fflush(stdout);
          return 0;
        }
        for (i=0; i<chrom->n_families; i++) {
            context = NULL;
            count = 0;
            while ((context=find_element(chrom->name[i], &context, beamline->elem_twiss))) {
                if (count==0 && (K2_param = confirm_parameter("K2", context->type))<0) {
                    fprintf(stdout, "error: element %s doesn't have K2 parameter\n",
                            context->name);
                    fflush(stdout);
                    exit(1);
                    }
                if (!(K2ptr = (double*)(context->p_elem + entity_description[context->type].parameter[K2_param].offset)))
                    bomb("K2ptr NULL in setup_chromaticity_correction", NULL);
                K2 = (*K2ptr += chrom->correction_fraction*chrom->dK2->a[i][0]);
                if (chrom->strengthLimit>0 && chrom->strengthLimit<fabs(K2)) {
                  K2 = *K2ptr = SIGN(K2)*chrom->strengthLimit;
                }
                if (context->matrix)
                    free_matrices(context->matrix);
                compute_matrix(context, run, NULL);
                type = context->type;
                count++;
                /* fprintf(stdout, "new value of %s#%ld[K2] is  %.15lg 1/m^3\n", 
                       chrom->name[i], context->occurence, K2);
                   fflush(stdout);
                       */
              }
            fprintf(stdout, "Change for family %ld (%ld sextupoles): %e\n",
                    i, count, chrom->correction_fraction*chrom->dK2->a[i][0]);
            fflush(stdout);
            if (alter_defined_values)
              change_defined_parameter(chrom->name[i], K2_param, type, K2, NULL, LOAD_FLAG_ABSOLUTE);
            }    

        if (beamline->links)
            assert_element_links(beamline->links, run, beamline, 
                                 STATIC_LINK+DYNAMIC_LINK+(alter_defined_values?LINK_ELEMENT_DEFINITION:0));

        M = beamline->matrix = compute_periodic_twiss(&beta_x, &alpha_x, &eta_x, &etap_x, beamline->tune,
                                                      &beta_y, &alpha_y, &eta_y, &etap_y, beamline->tune+1, 
                                                      beamline->elem_twiss, clorb, run, &unstable, NULL, NULL);
        beamline->twiss0->betax  = beta_x;
        beamline->twiss0->alphax = alpha_x;
        beamline->twiss0->phix   = 0;
        beamline->twiss0->etax   = eta_x;
        beamline->twiss0->etapx  = etap_x;
        beamline->twiss0->betay  = beta_y;
        beamline->twiss0->alphay = alpha_y;
        beamline->twiss0->phiy   = 0;
        beamline->twiss0->etay   = eta_y;
        beamline->twiss0->etapy  = etap_y;
        
        if (!M || !M->C || !M->R || !M->T)
            bomb("something wrong with transfer map for beamline (do_chromaticity_correction.2)", NULL);
        computeChromaticities(&chromx0, &chromy0, 
			      NULL, NULL, NULL, NULL, beamline->twiss0, M);
        beamline->chromaticity[0] = chromx0;
        beamline->chromaticity[1] = chromy0;
        fprintf(stdout, "resulting chromaticities:  %e  %e\n", chromx0, chromy0);
        fflush(stdout);
        }

    if (fp_sl && last_iteration) {
        for (i=0; i<chrom->n_families; i++) {
            context = NULL;
            while ((context=find_element(chrom->name[i], &context, beamline->elem_twiss))) {
              if (( K2_param = confirm_parameter("K2", context->type))<0)
                bomb("confirm_parameter doesn't return offset for K2 parameter.\n", NULL);
              fprintf(fp_sl, "%ld %e %s\n", 
                      step,
                      *((double*)(context->p_elem + entity_description[context->type].parameter[K2_param].offset)),
                      chrom->name[i]);
            }
          }
        fflush(fp_sl);
      }
    

    propagate_twiss_parameters(beamline->twiss0, beamline->tune, 
                               beamline->waists, NULL, beamline->elem_twiss, run, clorb);
    log_exit("do_chromaticity_correction");
    return 1;
    }

void computeChromaticities(double *chromx, double *chromy, 
                           double *dbetax, double *dbetay,
                           double *dalphax, double *dalphay,
                           TWISS *twiss, VMATRIX *M)
{
  double computeChromaticDerivRElem(long i, long j, TWISS *twiss, VMATRIX *M);
  double dR11, dR22, dR12, dR33, dR34, dR44;

  dR11 = computeChromaticDerivRElem(1, 1, twiss, M);
  dR12 = computeChromaticDerivRElem(1, 2, twiss, M);
  dR22 = computeChromaticDerivRElem(2, 2, twiss, M);
  dR33 = computeChromaticDerivRElem(3, 3, twiss, M);
  dR34 = computeChromaticDerivRElem(3, 4, twiss, M);
  dR44 = computeChromaticDerivRElem(4, 4, twiss, M);

  *chromx = -(dR11 + dR22)/M->R[0][1]*twiss->betax/(2*PIx2);
  *chromy = -(dR33 + dR44)/M->R[2][3]*twiss->betay/(2*PIx2);

  if (dbetax) {
    if (M->R[0][1])
      *dbetax = twiss->betax*(dR12 - PI*twiss->betax*(*chromx)*(M->R[0][0]+M->R[1][1]))/M->R[0][1];
    else
      *dbetax = DBL_MAX;
  }
  if (dalphax) {
    if (M->R[0][1])
      *dalphax = twiss->betax/(2*M->R[0][1])*
        (-twiss->alphax*(M->R[0][0]+M->R[1][1])*PIx2*(*chromx) + dR11-dR22);
    else
      *dalphax = DBL_MAX;
  }
  if (dbetay) {
    if (M->R[2][3])
      *dbetay = twiss->betay*(dR34 - PI*twiss->betay*(*chromy)*(M->R[2][2]+M->R[3][3]))/M->R[2][3];
    else
      *dbetay = DBL_MAX;
  }
  if (dalphay) {
    if (M->R[2][3])
      *dalphay = twiss->betay/(2*M->R[2][3])*
        (-twiss->alphay*(M->R[2][2]+M->R[3][3])*PIx2*(*chromy) + dR33-dR44);
    else
      *dalphay = DBL_MAX;
  }
}

void computeHigherOrderChromaticities(LINE_LIST *beamline, double *clorb, RUN *run,
				      long concatOrder)
{
#define NDELTA_VALUES 5   /* must be at least 5 */
  double trace[2][NDELTA_VALUES], delta[NDELTA_VALUES];
  double eta[6];
  double coef[NDELTA_VALUES], sCoef[NDELTA_VALUES], chi;
  long i, p;
  double c1;
  VMATRIX *M1, M0;

  beamline->chrom2[0] = beamline->chrom2[1] = 0;
  beamline->chrom3[0] = beamline->chrom3[1] = 0;
  initialize_matrices(&M0, 1);
  eta[0] = beamline->twiss0->etax;
  eta[1] = beamline->twiss0->etapx;
  eta[2] = beamline->twiss0->etay;
  eta[3] = beamline->twiss0->etapy;
  eta[4] = 0;
  eta[5] = 1;
  for (p=0; p<NDELTA_VALUES; p++) {
    delta[p] = (p-3)*1e-4;
    for (i=0; i<6; i++) {
      M0.C[i] = (clorb?clorb[i]:0)+delta[p]*eta[i] +
	(i<4? sqr(delta[p])*beamline->eta2[i] + pow3(delta[p])*beamline->eta3[i] : 0);
      M0.R[i][i] = 1;
    }
    M1 = append_full_matrix(beamline->elem_twiss, run, &M0, concatOrder);
    for (i=0; i<2; i++)
      /* Tr[0,1][p] is the trace for x,y plane for point p */
      trace[i][p] = M1->R[2*i][2*i] + M1->R[2*i+1][2*i+1];
    free_matrices(M1);
  }
  for (i=0; i<2; i++) {
    lsfn(delta, trace[i], NULL, NDELTA_VALUES, NDELTA_VALUES-2, coef, sCoef, &chi, NULL);
    c1 = -coef[1]/(4*PI*sin(PIx2*beamline->tune[i]));
    beamline->chrom2[i] = (
		 2*coef[2]
		 + 8*sqr(PI)*sqr(c1)*cos(PIx2*beamline->tune[i])
		 )/(-4*PI*sin(PIx2*beamline->tune[i]));
    beamline->chrom3[i] = (
		 6*coef[3] 
		 - 16*pow3(PI*c1)*sin(PIx2*beamline->tune[i]) 
		 + 24*sqr(PI)*c1*beamline->chrom2[i]*cos(PIx2*beamline->tune[i])
		 )/(-4*PI*sin(PIx2*beamline->tune[i]));
  }
  free_matrices(&M0);
}

double computeChromaticDerivRElem(long i, long j, TWISS *twiss, VMATRIX *M)
{
  long k;
  double sum, eta[6] = {0,0,0,0,0,0};

  if (!(M->T))
    return 0.0;
  
  eta[0] = twiss->etax;
  eta[1] = twiss->etapx;
  eta[2] = twiss->etay;
  eta[3] = twiss->etapy;
  eta[5] = 1;
  i--;
  j--;
  for (k=sum=0; k<6; k++) {
    if (k>j) {
      sum += eta[k]*M->T[i][k][j];
    } else if (k==j) {
      sum += eta[k]*M->T[i][k][k]*2;
    } else {
      sum += eta[k]*M->T[i][j][k];
    }
  }
  return sum;
}  

double computeChromaticDeriv2RElem(long i, long m, TWISS *twiss, VMATRIX *M,
                                   double *eta2q)
{
  long j, k, l, count;
  double sum, eta[6] = {0,0,0,0,0,0};
  double eta2[6] = {0,0,0,0,0,0};
   
  fprintf(stderr, "Computing d2R%ld%ld:\n", i, m);
  eta[0] = twiss->etax;
  eta[1] = twiss->etapx;
  eta[2] = twiss->etay;
  eta[3] = twiss->etapy;
  eta[5] = 1;
  for (j=0; j<4; j++)
    eta2[j] = eta2q[j];

  i--;
  m--;
  sum = 0;
  count = 0;
  if (M->T)
    for (k=0; k<6; k++) {
      if (k>m) {
        sum += eta2[k]*M->T[i][k][m];
      } else if (k==m) {
        sum += eta2[k]*M->T[i][k][k]*2;
      } else {
        sum += eta2[k]*M->T[i][m][k];
      }
      count ++;
    }
  if (M->Q) {
    for (j=0; j<6; j++)
      for (k=0; k<=j; k++)
	for (l=0; l<=k; l++) {
	  if (j==m) {
	    sum += M->Q[i][m][k][l]*eta[k]*eta[l];
	    count ++;
	  }
	  if (k==m) {
	    sum += M->Q[i][j][m][l]*eta[j]*eta[l];
	    count ++;
	  }
	  if (l==m) {
	    sum += M->Q[i][j][k][l]*eta[j]*eta[k];
	    count ++;
	  }
	}
  }
  fprintf(stderr, "sum is %e from %ld terms\n",
	  sum, count);
  return sum;
}  

double computeChromaticDeriv3RElem(long i, long m, TWISS *twiss, VMATRIX *M,
                                   double *eta2q, double *eta3q)
{
  long j, k, l, count;
  double sum, eta[6] = {0,0,0,0,0,0};
  double eta2[6] = {0,0,0,0,0,0};
  double eta3[6] = {0,0,0,0,0,0};
   
  fprintf(stderr, "Computing d3R%ld%ld:\n", i, m);
  eta[0] = twiss->etax;
  eta[1] = twiss->etapx;
  eta[2] = twiss->etay;
  eta[3] = twiss->etapy;
  eta[5] = 1;
  for (j=0; j<4; j++) {
    eta2[j] = eta2q[j];
    eta3[j] = eta3q[j];
  }
  for (j=0; j<6; j++) {
    fprintf(stderr, "eta[%ld] = %e, %e, ", j, eta[j], eta2[j]);
    fprintf(stderr, "%e\n", eta3[j]);
  }
  i--;
  m--;
  sum = count = 0;
  if (M->T)
    for (k=0; k<6; k++) {
      if (k>m) {
	count ++;
        sum += eta3[k]*M->T[i][k][m];
	fprintf(stderr, "term eta3[%ld]*T[%ld][%ld][%ld] = %e*%e\n",
		k, i, k, m, eta3[k], M->T[i][k][m]);
      } else if (k==m) {
	count ++;
        sum += eta3[k]*M->T[i][k][k]*2;
	fprintf(stderr, "term eta3[%ld]*T[%ld][%ld][%ld] = %e*%e\n",
		k, i, k, k, eta3[k], M->T[i][k][k]);
      } else {
	count ++;
        sum += eta3[k]*M->T[i][m][k];
	fprintf(stderr, "term eta3[%ld]*T[%ld][%ld][%ld] = %e*%e\n",
		k, i, m, k, eta3[k], M->T[i][m][k]);
      }
    }
  if (M->Q) {
    for (j=0; j<6; j++)
      for (k=0; k<=j; k++)
	for (l=0; l<=k; l++) {
	  if (j==m) {
	    count ++;
	    sum += M->Q[i][j][k][l]*(eta2[k]*eta[l]+eta[k]*eta2[l]);
	    fprintf(stderr, "term eta*U[%ld][%ld][%ld][%ld] = %e*%e\n",
		    i, j, k, l, (eta2[k]*eta[l]+eta[k]*eta2[l]), 
		    M->Q[i][j][k][l]);
	  }
	  if (k==m) {
	    count ++;
	    sum += M->Q[i][j][k][l]*(eta2[j]*eta[l]+eta[j]*eta2[l]);
	    fprintf(stderr, "term eta*U[%ld][%ld][%ld][%ld] = %e*%e\n",
		    i, j, k, l, (eta2[j]*eta[l]+eta[j]*eta2[l]),
		    M->Q[i][j][k][l]);
	  }
	  if (l==m) {
	    count ++;
	    sum += M->Q[i][j][k][l]*(eta2[j]*eta[k]+eta[j]*eta2[k]);
	    fprintf(stderr, "term eta*U[%ld][%ld][%ld][%ld] = %e*%e\n",
		    i, j, k, l, (eta2[j]*eta[k]+eta[j]*eta2[k]),
		    M->Q[i][j][k][l]);
	  }
	}
  }
  fprintf(stderr, "sum of %ld terms is %e\n", count, sum);
  return sum;
}  
