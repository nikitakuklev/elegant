/* Copyright 1994 by Michael Borland and Argonne National Laboratory,
 * all rights reserved.
 */
/* file: chrom.c
 * purpose: chromaticity correction by adjustment of sextupoles.
 *
 * Michael Borland, 1992
 */
#include "mdb.h"
#include "track.h"
#include "chromDefs.h"

static FILE *fp_sl = NULL;
static long alter_defined_values;
void computeChromaticities(double *chromx, double *chromy, TWISS *twiss, VMATRIX *M);

void setup_chromaticity_correction(NAMELIST_TEXT *nltext, RUN *run, LINE_LIST *beamline, CHROM_CORRECTION *chrom)
{
    VMATRIX *M;
    static long initial_call = 1;
    ELEMENT_LIST *eptr, *elast, *context;
    long i, n_elem, last_n_elem, count, K2_param;
    MATRIX *C, *Ct, *CtC, *inv_CtC;
#include "chrom.h"
    unsigned long unstable;
    
    log_entry("setup_chromaticity_correction");

    cp_str(&sextupoles, "sf sd");

    /* process namelist input */
    set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
    set_print_namelist_flags(0);
    process_namelist(&chromaticity, nltext);
    str_toupper(sextupoles);
    print_namelist(stderr, &chromaticity);

    if (run->default_order<2)
        bomb("default order must be >= 2 for chromaticity correction", NULL);

    if (chrom->name)
        tfree(chrom->name);
    chrom->name = tmalloc(sizeof(*chrom->name)*(chrom->n_families=1));
    while (chrom->name[chrom->n_families-1]=get_token(sextupoles))
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
    
    if (!beamline->twiss0 || !beamline->matrix) {
        double beta_x, alpha_x, eta_x, etap_x;
        double beta_y, alpha_y, eta_y, etap_y;

        fprintf(stderr, "Computing periodic Twiss parameters.\n");

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
                                                      &unstable, NULL);
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
        
        propagate_twiss_parameters(beamline->twiss0, beamline->tune, NULL, beamline->elem_twiss, run, NULL);
        }

    if (!(M=beamline->matrix) || !M->C || !M->R || !M->T)
        bomb("something wrong with transfer map for beamline (setup_chromaticity_correction)", NULL);

    if (!use_perturbed_matrix)
      computeChromCorrectionMatrix(run, beamline, chrom);

    if (strength_log) {
        strength_log = compose_filename(strength_log, run->rootname);
        fp_sl = fopen_e(strength_log, "w", FOPEN_SAVE_IF_EXISTS);
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
    double K2, *K2ptr;
    ELEMENT_LIST *eptr, *elast, *context;
    long i, n_elem, last_n_elem, count, K2_param;
    MATRIX *C, *Ct, *CtC, *inv_CtC;

    m_alloc(&C, 2, chrom->n_families);
    m_alloc(&Ct, chrom->n_families, 2);
    m_alloc(&CtC, chrom->n_families, chrom->n_families);
    m_alloc(&inv_CtC, chrom->n_families, chrom->n_families);
    m_alloc(&(chrom->T), chrom->n_families, 2);
    m_alloc(&(chrom->dK2), chrom->n_families, 1);
    m_alloc(&(chrom->dchrom), 2, 1);

    fprintf(stderr, "Computing chromaticity influence matrix for all named sextupoles.\n");
    computeChromaticities(&chromx0, &chromy0, beamline->twiss0, M=beamline->matrix);

    for (i=0; i<chrom->n_families; i++) {
        count = 0;
        context = NULL;
        while (context=find_element(chrom->name[i], &context, beamline->elem_twiss)) {
            if (!count && !(K2_param=confirm_parameter("K2", context->type))) {
                fprintf(stderr, "error: element %s does not have K2 parameter\n", 
                        context->name);
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
            fprintf(stderr, "error: element %s is not in the beamline.\n",  chrom->name[i]);
            exit(1);
            }
        if (beamline->links)
            assert_element_links(beamline->links, run, beamline, STATIC_LINK+DYNAMIC_LINK);
        M = full_matrix(beamline->elem_twiss, run, 2);
        computeChromaticities(&chromx, &chromy, beamline->twiss0, M);

        C->a[0][i] = (chromx-chromx0)/chrom->sextupole_tweek;
        C->a[1][i] = (chromy-chromy0)/chrom->sextupole_tweek;
        if (C->a[0][i]==0 || C->a[1][i]==0) {
            fprintf(stderr, "error: element %s does not change the chromaticity!\n", chrom->name[i]);
            exit(1);
            }
        count = 0;
        context = NULL;
        while (context=find_element(chrom->name[i], &context, beamline->elem_twiss)) {
            if (!count && !(K2_param=confirm_parameter("K2", context->type))) {
                fprintf(stderr, "error: element %s does not have K2 parameter\n", 
                        context->name);
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

    fprintf(stderr, "\nfamily           dCHROMx/dK2        dCHROMy/dK2\n");
    for (i=0; i<chrom->n_families; i++)
       fprintf(stderr, "%10s:    %14.7e     %14.7e\n", chrom->name[i], C->a[0][i], C->a[1][i]);

    m_trans(Ct, C);
    m_mult(CtC, Ct, C);
    m_invert(inv_CtC, CtC);
    m_mult(chrom->T, inv_CtC, Ct);

    fprintf(stderr, "\nfamily           dK2/dCHROMx        dK2/dCHROMy\n");
    for (i=0; i<chrom->n_families; i++)
       fprintf(stderr, "%10s:    %14.7e     %14.7e\n", chrom->name[i], chrom->T->a[i][0], chrom->T->a[i][1]);
    fprintf(stderr, "\n");

    m_free(&C);
    m_free(&Ct);
    m_free(&CtC);
    m_free(&inv_CtC);
  }



long do_chromaticity_correction(CHROM_CORRECTION *chrom, RUN *run, LINE_LIST *beamline,
            double *clorb, long step, long last_iteration)
{
    VMATRIX *M;
    double chromx0, chromy0;
    double K2, *K2ptr;
    ELEMENT_LIST *context;
    long i, K2_param, type, iter, count;
    double beta_x, alpha_x, eta_x, etap_x;
    double beta_y, alpha_y, eta_y, etap_y;
    unsigned long unstable;

    
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
                                                      beamline->elem_twiss, clorb, run, &unstable, NULL);
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
        
        propagate_twiss_parameters(beamline->twiss0, beamline->tune, NULL, beamline->elem_twiss, run, clorb);
        }
    else if (beamline->matrix->order<2)
        beamline->matrix = full_matrix(beamline->elem_twiss, run, 2);

    if (!(M = beamline->matrix) || !M->C || !M->R || !M->T)
        bomb("something wrong with transfer map for beamline (do_chromaticity_correction.1)", NULL);

    if (chrom->use_perturbed_matrix)
      computeChromCorrectionMatrix(run, beamline, chrom);
    
    computeChromaticities(&chromx0, &chromy0, beamline->twiss0, M);

    fprintf(stderr, "\nAdjusting chromaticities:\n");
    fprintf(stderr, "initial chromaticities:  %e  %e\n", chromx0, chromy0);

    for (iter=0; iter<chrom->n_iterations; iter++) {
        chrom->dchrom->a[0][0] = chrom->chromx - chromx0;
        chrom->dchrom->a[1][0] = chrom->chromy - chromy0;
        if (chrom->tolerance>0 &&
            chrom->tolerance>fabs(chrom->dchrom->a[0][0]) &&
            chrom->tolerance>fabs(chrom->dchrom->a[1][0]) )
          break;

        m_mult(chrom->dK2, chrom->T, chrom->dchrom);
        for (i=0; i<chrom->n_families; i++) {
          if (isnan(chrom->correction_fraction*chrom->dK2->a[i][0]) ||
              isinf(chrom->correction_fraction*chrom->dK2->a[i][0]))
            break;
        }
        if (i!=chrom->n_families) {
          fprintf(stderr, "Unable to correct chromaticity---diverged\n");
          return 0;
        }
        for (i=0; i<chrom->n_families; i++) {
            context = NULL;
            count = 0;
            while (context=find_element(chrom->name[i], &context, beamline->elem_twiss)) {
                if (count==0 && (K2_param = confirm_parameter("K2", context->type))<0) {
                    fprintf(stderr, "error: element %s doesn't have K2 parameter\n",
                            context->name);
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
                fprintf(stderr, "new value of %s#%ld[K2] is  %.15lg 1/m^3\n", 
                       chrom->name[i], context->occurence, K2);
                }
            if (alter_defined_values)
              change_defined_parameter(chrom->name[i], K2_param, type, K2, NULL, LOAD_FLAG_ABSOLUTE);
            }    

        if (beamline->links)
            assert_element_links(beamline->links, run, beamline, 
                                 STATIC_LINK+DYNAMIC_LINK+(alter_defined_values?LINK_ELEMENT_DEFINITION:0));

        M = beamline->matrix = compute_periodic_twiss(&beta_x, &alpha_x, &eta_x, &etap_x, beamline->tune,
                                                      &beta_y, &alpha_y, &eta_y, &etap_y, beamline->tune+1, 
                                                      beamline->elem_twiss, clorb, run, &unstable, NULL);
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
        computeChromaticities(&chromx0, &chromy0, beamline->twiss0, M);
        beamline->chromaticity[0] = chromx0;
        beamline->chromaticity[1] = chromy0;
        fprintf(stderr, "resulting chromaticities:  %e  %e\n", chromx0, chromy0);
        }

    if (fp_sl && last_iteration) {
        for (i=0; i<chrom->n_families; i++) {
            context = NULL;
            while (context=find_element(chrom->name[i], &context, beamline->elem_twiss)) {
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
    

    propagate_twiss_parameters(beamline->twiss0, beamline->tune, NULL, beamline->elem_twiss, run, clorb);
    log_exit("do_chromaticity_correction");
    return 1;
    }

void computeChromaticities(double *chromx, double *chromy, TWISS *twiss, VMATRIX *M)
{
    double dR11, dR22, dR33, dR44;

    dR11 = M->T[0][5][0] + 2*M->T[0][0][0]*twiss->etax  + M->T[0][1][0]*twiss->etapx +
        M->T[0][2][0]*twiss->etay + M->T[0][3][0]*twiss->etapy;
    dR22 = M->T[1][5][1] + 2*M->T[1][1][1]*twiss->etapx + M->T[1][1][0]*twiss->etax  +
        M->T[1][2][1]*twiss->etay + M->T[1][3][1]*twiss->etapy;
    *chromx = -(dR11 + dR22)/M->R[0][1]*twiss->betax/(2*PIx2);
    dR33 = M->T[2][5][2] + 2*M->T[2][2][2]*twiss->etay  + M->T[2][3][2]*twiss->etapy +
        M->T[2][2][0]*twiss->etax + M->T[2][2][1]*twiss->etapx;
    dR44 = M->T[3][5][3] + 2*M->T[3][3][3]*twiss->etapy + M->T[3][3][2]*twiss->etay  +
        M->T[3][3][0]*twiss->etax + M->T[3][3][1]*twiss->etapx;
    *chromy = -(dR33 + dR44)/M->R[2][3]*twiss->betay/(2*PIx2);
    }
