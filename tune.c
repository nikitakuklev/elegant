/* Copyright 1994 by Michael Borland and Argonne National Laboratory,
 * all rights reserved.
 */
/* file: tune.c
 * purpose: tune correction by adjustment of quadrupoles
 *
 * Michael Borland, 1992
 */
#include "mdb.h"
#include "track.h"
#include "tuneDefs.h"

static FILE *fp_sl = NULL;
static long alter_defined_values;

void setup_tune_correction(NAMELIST_TEXT *nltext, RUN *run, LINE_LIST *beamline, TUNE_CORRECTION *tune)
{
    VMATRIX *M;
    ELEMENT_LIST *eptr, *elast, *context;
    double beta_x, alpha_x, eta_x, etap_x;
    double beta_y, alpha_y, eta_y, etap_y;
    long i, n_elem, last_n_elem, count;

#include "tune.h"

    log_entry("setup_tune_correction");

    if (fp_sl) {
        fclose(fp_sl);
        fp_sl = NULL;
        }

    /* process namelist input */
    set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
    set_print_namelist_flags(0);
    process_namelist(&correct_tunes, nltext);
    str_toupper(quadrupoles);
    print_namelist(stderr, &correct_tunes);

    if (tune->name)
        tfree(tune->name);
    tune->name = tmalloc(sizeof(*tune->name)*(tune->n_families=1));
    while (tune->name[tune->n_families-1]=get_token(quadrupoles))
        tune->name = trealloc(tune->name, sizeof(*tune->name)*(tune->n_families+=1));
    if ((--tune->n_families)<2)
        bomb("too few quadrupoles given for tune correction", NULL);
    tune->tunex = tune_x;
    tune->tuney = tune_y;
    tune->gain = correction_fraction;
    tune->n_iterations = n_iterations;
    tune->use_perturbed_matrix = use_perturbed_matrix;
    tune->tolerance = tolerance;
    tune->maximum_gain = max_correction_fraction;
    tune->step_up_interval = step_up_interval;
    tune->delta_gain = delta_correction_fraction;
    alter_defined_values = change_defined_values;

    if (!beamline->twiss0 || !beamline->matrix) {
        if (!beamline->twiss0)
            beamline->twiss0 = tmalloc(sizeof(*beamline->twiss0));

        eptr = beamline->elem_twiss = &(beamline->elem);
        n_elem = last_n_elem = beamline->n_elems;
        elast = eptr;
        while (eptr) {
            if (eptr->type==T_RECIRC) {
                last_n_elem = n_elem;
                beamline->elem_twiss = beamline->elem_recirc = eptr;
                }
            elast = eptr;
            eptr = eptr->succ;
            n_elem --;
            }
        n_elem = last_n_elem;
    
        }

    fprintf(stderr, "Computing periodic Twiss parameters.\n");
    M = beamline->matrix = compute_periodic_twiss(&beta_x, &alpha_x, &eta_x, &etap_x, beamline->tune,
            &beta_y, &alpha_y, &eta_y, &etap_y, beamline->tune+1, beamline->elem_twiss, NULL, run);
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

    if (tune->tunex<0)
        fprintf(stderr, "horizontal tune will be held at %f\n", tune->tunex = beamline->tune[0]);
    if (tune->tuney<0)
        fprintf(stderr, "vertical tune will be held at %f\n", tune->tuney = beamline->tune[1]);
    
    if (strength_log) {
        strength_log = compose_filename(strength_log, run->rootname);
        fp_sl = fopen_e(strength_log, "w", FOPEN_SAVE_IF_EXISTS);
        fprintf(fp_sl, "SDDS1\n&column name=Step, type=long, description=\"Simulation step\" &end\n");
        fprintf(fp_sl, "&column name=K1, type=double, units=\"1/m$a2$n\" &end\n");
        fprintf(fp_sl, "&column name=QuadrupoleName, type=string  &end\n");
        fprintf(fp_sl, "&data mode=ascii, no_row_counts=1 &end\n");
        fflush(fp_sl);
        }

    if (!tune->use_perturbed_matrix)
      computeTuneCorrectionMatrix(run, beamline, tune, 1);
}

void computeTuneCorrectionMatrix(RUN *run, LINE_LIST *beamline, TUNE_CORRECTION *tune, long printout)
{
    MATRIX *C, *Ct, *CtC, *inv_CtC;
    VMATRIX *M;
    long i, n_elem, last_n_elem, count;
    double betax_L_sum, betay_L_sum;
    ELEMENT_LIST *eptr, *elast, *context;
    
    if (!(M=beamline->matrix) || !M->C || !M->R)
        bomb("something wrong with transfer map for beamline (setup_tune_correction)", NULL);

    m_alloc(&C, 2, tune->n_families);
    m_alloc(&Ct, tune->n_families, 2);
    m_alloc(&CtC, tune->n_families, tune->n_families);
    m_alloc(&inv_CtC, tune->n_families, tune->n_families);
    m_alloc(&(tune->T), tune->n_families, 2);
    m_alloc(&(tune->dK1), tune->n_families, 1);
    m_alloc(&(tune->dtune), 2, 1);

    fprintf(stderr, "Computing tune influence matrix for all named quadrupoles.\n");

    for (i=0; i<tune->n_families; i++) {
        count = 0;
        context = NULL;
        betax_L_sum = betay_L_sum = 0;
        while (context=find_element(tune->name[i], &context, &(beamline->elem))) {
            if (count==0) {
                if (context->type!=T_QUAD && context->type!=T_KQUAD) {
                    fprintf(stderr, "%s is not a QUAD or KQUAD element!\n", context->name);
                    exit(1);
                    }
                }
            betax_L_sum += context->twiss->betax*((QUAD*)context->p_elem)->length;
            betay_L_sum += context->twiss->betay*((QUAD*)context->p_elem)->length;
            count++;
            }
        if (count==0) {
            fprintf(stderr, "error: element %s is not part of the beamline\n", tune->name[i]);
            exit(1);
            }
        if (printout)
          fprintf(stderr, "%ld instances of %s found\n", count, tune->name[i]); 
        
        C->a[0][i] = betax_L_sum/(4*PI);
        C->a[1][i] = -betay_L_sum/(4*PI);
        if (C->a[0][i]==0 || C->a[1][i]==0) {
            fprintf(stderr, "error: element %s does not change the tune!\n", tune->name[i]);
            exit(1);
            }
        }

    if (printout) {
      fprintf(stderr, "\nfamily               dNUx/dK1                  dNUy/dK1\n");

      for (i=0; i<tune->n_families; i++)
        fprintf(stderr, "%10s:    %22.15le     %22.15le\n", tune->name[i], C->a[0][i], C->a[1][i]);
    }
    
    m_trans(Ct, C);
    m_mult(CtC, Ct, C);
    m_invert(inv_CtC, CtC);
    m_mult(tune->T, inv_CtC, Ct);

    if (printout) {
      fprintf(stderr, "\nfamily               dK1/dNUx                  dK1/dNUy\n");
      for (i=0; i<tune->n_families; i++)
        fprintf(stderr, "%10s:    %22.15le     %22.15le\n", tune->name[i], tune->T->a[i][0], tune->T->a[i][1]);
      fprintf(stderr, "\n");
    }
    
    m_free(&C);
    m_free(&Ct);
    m_free(&CtC);
    m_free(&inv_CtC);

    log_exit("setup_tune_correction");
    }


long do_tune_correction(TUNE_CORRECTION *tune, RUN *run, LINE_LIST *beamline, 
                        double *clorb, long step, long last_iteration)
{
  VMATRIX *M;
  double K1, gain, LastMsError, MsError;
  long steps_since_gain_change;
  ELEMENT_LIST *context;
  long i, K1_param, type, iter;
  double beta_x, alpha_x, eta_x, etap_x;
  double beta_y, alpha_y, eta_y, etap_y;
  static long tunes_saved=0;
  static double nux_orig, nuy_orig;
  
  if (tune->use_perturbed_matrix)
    computeTuneCorrectionMatrix(run, beamline, tune, 0);
  
  M = beamline->matrix = compute_periodic_twiss(&beta_x, &alpha_x, &eta_x, &etap_x, beamline->tune,
                                                &beta_y, &alpha_y, &eta_y, &etap_y, beamline->tune+1, 
                                                beamline->elem_twiss, clorb, run);
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

  if (!M || !M->C || !M->R)
    bomb("something wrong with transfer map for beamline (do_tune_correction.1)", NULL);

  fprintf(stderr, "\nAdjusting tunes:\n");
  fprintf(stderr, "initial tunes:  %e  %e\n", beamline->tune[0], beamline->tune[1]);
  MsError = sqr(beamline->tune[0]-tune->tunex)+sqr(beamline->tune[1]-tune->tuney);
  
  if (!tunes_saved) {
    nux_orig = beamline->tune[0];
    nuy_orig = beamline->tune[1];
    tunes_saved = 1;
  }

  steps_since_gain_change = 0;
  gain = tune->gain;
  
  for (iter=0; iter<tune->n_iterations; iter++) {
    LastMsError = MsError;
    MsError = sqr(beamline->tune[0]-tune->tunex)+sqr(beamline->tune[1]-tune->tuney);
    if (MsError>LastMsError) {
      /* reset to minimum gain */
      fprintf(stderr, "Warning: tune correction diverging---gain reset to minimum\n");
      gain = tune->gain;
      steps_since_gain_change = 0;
    }
    
    tune->dtune->a[0][0] = tune->tunex - beamline->tune[0];
    tune->dtune->a[1][0] = tune->tuney - beamline->tune[1];
    if (tune->tolerance>0 &&
        tune->tolerance>fabs(tune->dtune->a[0][0]) &&
        tune->tolerance>fabs(tune->dtune->a[1][0]))
      break;
    
    if (( K1_param = confirm_parameter("K1", T_QUAD))<0)
      bomb("confirm_parameter doesn't return offset for K1 parameter of quadrupole!\n", NULL);
    
    m_mult(tune->dK1, tune->T, tune->dtune);
    m_scmul(tune->dK1, tune->dK1, gain);
    for (i=0; i<tune->n_families; i++) {
      if (isnan(tune->dK1->a[i][0]) || isinf(tune->dK1->a[i][0]))
        break;
    }
    if (i!=tune->n_families) {
      fprintf(stderr, "Unable to correct tune---diverged.\n");
      return 0;
    }
    
    for (i=0; i<tune->n_families; i++) {
      context = NULL;
      while (context=find_element(tune->name[i], &context, &(beamline->elem))) {
        K1 = (((QUAD*)context->p_elem)->k1 += tune->dK1->a[i][0]);
        if (context->matrix)
          free_matrices(context->matrix);
        compute_matrix(context, run, NULL);
        type = context->type;
      }
      fprintf(stderr, "change of %s[K1] is  %.15g 1/m^3\n", tune->name[i], tune->dK1->a[i][0]);
      if (alter_defined_values) {
        fprintf(stderr, "new value of %s[K1] is  %.15g 1/m^3\n", tune->name[i], K1);
        change_defined_parameter(tune->name[i], K1_param, type, K1, NULL, LOAD_FLAG_ABSOLUTE);
      }
    }    
    if (++steps_since_gain_change==tune->step_up_interval) {
      if ((gain += tune->delta_gain)>tune->maximum_gain)
        gain = tune->maximum_gain;
      steps_since_gain_change = 0;
    }

    if (beamline->links)
      assert_element_links(beamline->links, run, beamline, 
                           STATIC_LINK+DYNAMIC_LINK+(alter_defined_values?LINK_ELEMENT_DEFINITION:0));

    M = beamline->matrix = compute_periodic_twiss(&beta_x, &alpha_x, &eta_x, &etap_x, beamline->tune,
                                                  &beta_y, &alpha_y, &eta_y, &etap_y, beamline->tune+1, beamline->elem_twiss, clorb, run);

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
    
    if (tune->use_perturbed_matrix)
      computeTuneCorrectionMatrix(run, beamline, tune, 0);
  
    if (!M || !M->C || !M->R)
      bomb("something wrong with transfer map for beamline (do_tune_correction.2)", NULL);

    propagate_twiss_parameters(beamline->twiss0, beamline->tune, NULL, beamline->elem_twiss, run, clorb);
    fprintf(stderr, "new tunes: %e %e\n", beamline->tune[0], beamline->tune[1]);
  }

  if (fp_sl && last_iteration) {
    tunes_saved = 0;
    for (i=0; i<tune->n_families; i++) {
      context = NULL;
      while (context=find_element(tune->name[i], &context, &(beamline->elem)))
        fprintf(fp_sl, "%ld %21.15e %s\n", step, ((QUAD*)context->p_elem)->k1, tune->name[i]);
    }    
    fflush(fp_sl);
  }

  return 1;
}
