/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: correct.c
 * purpose: trajectory/orbit correction for elegant
 * 
 * Michael Borland, 1991. 
 */
#include "mdb.h"
#include "mdbsun.h"
#include "track.h"
#include "match_string.h"
#include "correctDefs.h"

#define N_CORRECTION_MODES 2
char *correction_mode[N_CORRECTION_MODES] = {
    "trajectory", "orbit"
    } ;

#define GLOBAL_CORRECTION 0
#define ONE_TO_ONE_CORRECTION 1
#define N_CORRECTION_METHODS 2
char *correction_method[N_CORRECTION_METHODS] = {
    "global", "one-to-one"
    } ;

/* For trajectory correction:
 *   C      :     NMON x NCOR matrix of dQ_mon/dkick_cor 
 *   Qo     :     NMON x 1    matrix of initial monitor positions, Q stands for either x or y 
 *   Q      :     NMON x 1    matrix of final monitor positions
 *   dK     :     NCOR x 1    matrix of changes to correctors 
 *
 *   Q = Qo + C*dK
 *   Want to minimize S = SUM((Qi*wi)^2), where wi is the weight for the ith monitor.
 *   Solution is dK = -Inverse(Trans(C) W C) Trans(C) W Qo, or dK = T Qo.
 *   W is NMON x NMON, and W(i,j) = delta(i,j)*wi.  A monitor with a zero or negative weight is
 *       ignored.
 *   T is NCOR x NMON, and is computed before any errors are added and before any parameters are varied
 */

long global_trajcor_plane(CORMON_DATA *CM, STEERING_LIST *SL, long coord, TRAJECTORY **traject, long n_iterations, 
                          RUN *run, LINE_LIST *beamline, double *starting_coord, BEAM *beam);
void one_to_one_trajcor_plane(CORMON_DATA *CM, STEERING_LIST *SL, long coord, TRAJECTORY **traject, long n_iterations, RUN *run, 
            LINE_LIST *beamline, double *starting_coord, BEAM *beam);
long orbcor_plane(CORMON_DATA *CM, STEERING_LIST *SL, long coord, TRAJECTORY **orbit, 
                  long n_iterations, double accuracy,
                  long clorb_iter, double clorb_iter_frac,
                  RUN *run, LINE_LIST *beamline, double *closed_orbit, double *Cdp);
ELEMENT_LIST *find_useable_moni_corr(long *nmon, long *ncor, long **mon_index,
            ELEMENT_LIST ***umoni, ELEMENT_LIST ***ucorr, double **kick_coef, long **sl_index,
            long plane, STEERING_LIST *SL, RUN *run, LINE_LIST *beamline, long recircs);
ELEMENT_LIST *next_element_of_type(ELEMENT_LIST *elem, long type);
ELEMENT_LIST *next_element_of_type2(ELEMENT_LIST *elem, long type1, long type2);
ELEMENT_LIST *next_element_of_types(ELEMENT_LIST *elem, long *type, long n_types, long *index);
long find_parameter_offset(char *param_name, long elem_type);
long zero_correctors_one_plane(ELEMENT_LIST *elem, RUN *run, STEERING_LIST *SL);
long zero_correctors(ELEMENT_LIST *elem, RUN *run, CORRECTION *correct);
long zero_hcorrectors(ELEMENT_LIST *elem, RUN *run, CORRECTION *correct);
long zero_vcorrectors(ELEMENT_LIST *elem, RUN *run, CORRECTION *correct);
double rms_value(double *data, long n_data);
long steering_corrector(ELEMENT_LIST *eptr, STEERING_LIST *SL);
void zero_closed_orbit(TRAJECTORY *clorb, long n);
long find_index(long key, long *list, long n_listed);
void add_steer_elem_to_lists(STEERING_LIST *SL, long plane, char *name, char *item, 
                             char *element_type, double tweek, double limit,
                             LINE_LIST *beamline, RUN *run);
void add_steer_type_to_lists(STEERING_LIST *SL, long plane, long type, char *item, double tweek, double limit,
    LINE_LIST *beamline, RUN *run);
double compute_kick_coefficient(ELEMENT_LIST *elem, long plane, long type, double corr_tweek, char *name, char *item, RUN *run);
double noise_value(double xamplitude, double xcutoff, long xerror_type);
void do_response_matrix_output(char *filename, char *type, RUN *run, char *beamline_name, CORMON_DATA *CM, 
                               STEERING_LIST *SL, long plane);
long findFixedLengthClosedOrbit(TRAJECTORY *clorb, double clorb_acc, long clorb_iter, LINE_LIST *beamline, 
                                VMATRIX *M, RUN *run, double dp, long start_from_recirc, double *starting_point, 
                                double change_fraction, double *deviation);

static long rpn_x_mem= -1, rpn_y_mem= -1;
static long usePerturbedMatrix = 0, fixedLengthMatrix = 0;

double getMonitorWeight(ELEMENT_LIST *elem);
double getMonitorCalibration(ELEMENT_LIST *elem, long coord);
double getCorrectorCalibration(ELEMENT_LIST *elem, long coord);

#define UNIFORM_ERRORS 0
#define GAUSSIAN_ERRORS 1
#define PLUS_OR_MINUS_ERRORS 2
#define N_ERROR_TYPES 3
static char *known_error_type[N_ERROR_TYPES] = {
    "uniform", "gaussian", "plus_or_minus"
    };


void correction_setup(
    CORRECTION *_correct,     /* the underscore is to avoid conflicts with the namelist "correct" */
    NAMELIST_TEXT *nltext, RUN *run, LINE_LIST *beamline
    )
{
#include "correct.h"
    char *item;

    log_entry("correction_setup");

    if (_correct->traj) {
        free_zarray_2d((void**)_correct->traj, 3, beamline->n_elems+1);
        _correct->traj = NULL;
        }

    if (_correct->CMx) {
        m_free(&_correct->CMx->T);
        m_free(&_correct->CMx->dK);
        m_free(&_correct->CMx->Qo);
        _correct->CMx->T = _correct->CMx->dK = _correct->CMx->Qo = _correct->CMx->C = NULL;
        free_zarray_2d((void**)_correct->CMx->kick, _correct->n_iterations+1, _correct->CMx->ncor);
        free_zarray_2d((void**)_correct->CMx->posi, _correct->n_iterations+1, _correct->CMx->nmon);
        tfree(_correct->CMx->mon_index); _correct->CMx->mon_index = NULL;
        tfree(_correct->CMx); _correct->CMx = NULL;
        _correct->CMx->ncor = _correct->CMx->nmon = _correct->CMx->inverse_computed = 0;
        }

    if (_correct->CMy) {
        m_free(&_correct->CMy->T);
        m_free(&_correct->CMy->dK);
        m_free(&_correct->CMy->Qo);
        _correct->CMy->T = _correct->CMy->dK = _correct->CMy->Qo = _correct->CMy->C = NULL;
        free_zarray_2d((void**)_correct->CMy->kick, _correct->n_iterations+1, _correct->CMy->ncor);
        free_zarray_2d((void**)_correct->CMy->posi, _correct->n_iterations+1, _correct->CMy->nmon);
        tfree(_correct->CMy->mon_index); _correct->CMy->mon_index = NULL;
        tfree(_correct->CMy); _correct->CMy = NULL;
        _correct->CMy->ncor = _correct->CMy->nmon = _correct->CMy->inverse_computed = 0;
        }

    /* process the namelist text */
    set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
    set_print_namelist_flags(0);
    process_namelist(&correct, nltext);
    print_namelist(stdout, &correct);
    usePerturbedMatrix = use_perturbed_matrix;
    fixedLengthMatrix = fixed_length_matrix;
    
    /* check for valid input data */
    if ((_correct->mode=match_string(mode, correction_mode, N_CORRECTION_MODES, 0))<0)
        bomb("invalid correction mode", NULL);
    if ((_correct->method=match_string(method, correction_method, N_CORRECTION_METHODS, 0))<0)
        bomb("invalid correction method", NULL);
    if (corrector_tweek[0]==0 || corrector_tweek[0]==0)
        bomb("invalid corrector tweek(s)", NULL);
    if (corrector_limit[0]<0 || corrector_limit[1]<0 ||
        (corrector_limit[0] && corrector_limit[0]<corrector_tweek[0]) ||
        (corrector_limit[1] && corrector_limit[1]<corrector_tweek[1]) )
        bomb("invalid corrector_limit(s)", NULL);
    if (correction_fraction[0]<=0 || correction_fraction[1]<=0)
        bomb("invalid correction_fraction(s)", NULL);
    if (correction_accuracy[0]<0 || correction_accuracy[1]<0)
        bomb("invalid correction_accuracy(s)", NULL);
    if (trajectory_output)
        setup_orb_traj_output(trajectory_output=compose_filename(trajectory_output, run->rootname),
                    correction_mode[_correct->mode], run);
    if (statistics)
        setup_cormon_stats(statistics=compose_filename(statistics, run->rootname), run);
    if (corrector_output)
        setup_corrector_output(corrector_output=compose_filename(corrector_output, run->rootname), run);
    if (closed_orbit_accuracy<=0)
        bomb("closed_orbit_accuracy must be > 0", NULL);
    if (closed_orbit_iteration_fraction<=0 ||
        closed_orbit_iteration_fraction>1)
      bomb("closed_orbit_iteration_fraction must be on (0, 1]", NULL);
    _correct->clorb_accuracy = closed_orbit_accuracy;
    _correct->clorb_iter_fraction = closed_orbit_iteration_fraction;
    _correct->verbose = verbose;
    _correct->track_before_and_after = track_before_and_after;
    _correct->prezero_correctors = prezero_correctors;
    _correct->start_from_centroid = start_from_centroid;
    _correct->use_actual_beam = use_actual_beam;
    if ((_correct->clorb_iterations=closed_orbit_iterations)<=0)
        bomb("closed_orbit_iterations <= 0", NULL);
    if ((_correct->n_iterations = n_iterations)<0)
        bomb("n_iterations < 0", NULL);
    if ((_correct->n_xy_cycles = n_xy_cycles)<0)
        bomb("n_xy_cycles < 0", NULL);

    _correct->CMx = tmalloc(sizeof(*_correct->CMx));
    _correct->CMy = tmalloc(sizeof(*_correct->CMy));

    _correct->CMx->default_tweek = corrector_tweek[0];
    _correct->CMy->default_tweek = corrector_tweek[1];
    _correct->CMx->corr_limit = corrector_limit[0];
    _correct->CMy->corr_limit = corrector_limit[1];
    _correct->CMx->corr_fraction = correction_fraction[0];
    _correct->CMy->corr_fraction = correction_fraction[1];
    _correct->CMx->corr_accuracy = correction_accuracy[0];
    _correct->CMy->corr_accuracy = correction_accuracy[1];
    _correct->CMx->bpm_noise = bpm_noise[0];
    _correct->CMy->bpm_noise = bpm_noise[1];
    if ((_correct->CMx->bpm_noise_distribution 
         = match_string(bpm_noise_distribution[0], known_error_type, N_ERROR_TYPES, 0))<0)
      bomb("unknown noise distribution type", NULL);
    if ((_correct->CMy->bpm_noise_distribution 
         = match_string(bpm_noise_distribution[1], known_error_type, N_ERROR_TYPES, 0))<0)
      bomb("unknown noise distribution type", NULL);
    _correct->CMx->bpm_noise_cutoff = bpm_noise_cutoff[0];
    _correct->CMy->bpm_noise_cutoff = bpm_noise_cutoff[1];
    _correct->CMx->fixed_length = _correct->CMy->fixed_length = fixed_length;
    _correct->response_only = n_iterations==0;

    /* find correction matrices Qo, T, C, and W for all monitors using all correctors */
    if (verbose)
      fputs("finding correctors/monitors and/or computing correction matrices\n", stdout);
    
    if (_correct->SLx.n_corr_types==0) {
      cp_str(&item, "KICK");
      add_steer_type_to_lists(&_correct->SLx, 0, T_HCOR, item, _correct->CMx->default_tweek, 
                              _correct->CMx->corr_limit, beamline, run);
      cp_str(&item, "HKICK");
      add_steer_type_to_lists(&_correct->SLx, 0, T_HVCOR, item, _correct->CMx->default_tweek, 
                              _correct->CMy->corr_limit, beamline, run);
    }
    if (_correct->SLy.n_corr_types==0) {
      cp_str(&item, "KICK");
      add_steer_type_to_lists(&_correct->SLy, 2, T_VCOR, item, _correct->CMy->default_tweek, 
                              _correct->CMx->corr_limit, beamline, run);
      cp_str(&item, "VKICK");
      add_steer_type_to_lists(&_correct->SLy, 2, T_HVCOR, item, _correct->CMy->default_tweek, 
                              _correct->CMy->corr_limit, beamline, run);
    }

    if (_correct->mode==TRAJECTORY_CORRECTION) {
      compute_trajcor_matrices(_correct->CMx, &_correct->SLx, 0, run, beamline, 0,
                               !_correct->response_only);
      compute_trajcor_matrices(_correct->CMy, &_correct->SLy, 2, run, beamline, 0, 
                               !_correct->response_only);
    }
    else if (_correct->mode==ORBIT_CORRECTION) {
      compute_orbcor_matrices(_correct->CMx, &_correct->SLx, 0, run, beamline, 0, 
                              !_correct->response_only, fixed_length_matrix);
      compute_orbcor_matrices(_correct->CMy, &_correct->SLy, 2, run, beamline, 0, 
                              !_correct->response_only, fixed_length_matrix);
    }
    else
      bomb("something impossible happened (correction_setup)", NULL);

    if (n_iterations!=0) {
      /* allocate space to store before/after data for correctors and monitors */
      _correct->CMx->kick = (double**)zarray_2d(sizeof(**(_correct->CMx->kick)), _correct->n_iterations+1, _correct->CMx->ncor);
      _correct->CMx->posi = (double**)zarray_2d(sizeof(**(_correct->CMx->posi)), _correct->n_iterations+1, _correct->CMx->nmon);
      _correct->CMy->kick = (double**)zarray_2d(sizeof(**(_correct->CMy->kick)), _correct->n_iterations+1, _correct->CMy->ncor);
      _correct->CMy->posi = (double**)zarray_2d(sizeof(**(_correct->CMy->posi)), _correct->n_iterations+1, _correct->CMy->nmon);
      
      /* Allocate space to store before/after trajectories/closed orbits.
       * After each correction pass through x and y, the first trajectory is the initial one,
       * the second is after x correction, and the third is after x and y correction.
       */
      _correct->traj = (TRAJECTORY**)zarray_2d(sizeof(**_correct->traj), 3, beamline->n_elems+1);
      
    }
    if (verbose) {
      fprintf(stdout, "there are %ld useable horizontal monitors and %ld useable horizontal correctors\n",
              _correct->CMx->nmon, _correct->CMx->ncor);
      fflush(stdout);
      fprintf(stdout, "there are %ld useable   vertical monitors and %ld useable   vertical correctors\n",
              _correct->CMy->nmon, _correct->CMy->ncor);
      fflush(stdout);
    }
    rpn_x_mem = rpn_create_mem("x");
    rpn_y_mem = rpn_create_mem("y");

    log_exit("correction_setup");
  }

void add_steering_element(CORRECTION *correct, LINE_LIST *beamline, RUN *run, NAMELIST_TEXT *nltext)
{
#include "steer_elem.h"

  /* process the namelist text */
  set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
  set_print_namelist_flags(0);
  process_namelist(&steering_element, nltext);
  print_namelist(stdout, &steering_element);

  if (limit && (limit<tweek || limit<0))
    bomb("invalid limit specified for steering element", NULL);

  if (plane[0]=='h' || plane[0]=='H') 
    add_steer_elem_to_lists(&correct->SLx, 0, name, item, element_type, tweek, limit, beamline, run);
  else if (plane[0]=='v' || plane[0]=='V')
    add_steer_elem_to_lists(&correct->SLy, 2, name, item, element_type, tweek, limit, beamline, run);
  else
    bomb("invalid plane specified for steering element", NULL);
}

void add_steer_type_to_lists(STEERING_LIST *SL, long plane, long type, char *item, double tweek, double limit,
                             LINE_LIST *beamline, RUN *run)
{
  ELEMENT_LIST *context;
  context = &(beamline->elem);
  while (context && (context=next_element_of_type(context, type))) {
    add_steer_elem_to_lists(SL, plane, context->name, item, NULL, tweek, limit, beamline, run);
    context = context->succ;
  }
}

void add_steer_elem_to_lists(STEERING_LIST *SL, long plane, char *name, char *item, 
                             char *element_type, double tweek, double limit, 
                             LINE_LIST *beamline, RUN *run)
{
  ELEMENT_LIST *context;
  long param_number, i, found;

  if (SL->n_corr_types==0) {
    if (SL->corr_name)    tfree(SL->corr_name);
    SL->corr_name = NULL;
    if (SL->corr_type)    tfree(SL->corr_type);
    SL->corr_type = NULL;
    if (SL->corr_param)   tfree(SL->corr_param);
    SL->corr_param = NULL;
    if (SL->corr_tweek)   tfree(SL->corr_tweek);
    SL->corr_tweek = NULL;
    if (SL->corr_limit)   tfree(SL->corr_limit);
    SL->corr_limit = NULL;
    if (SL->param_offset) tfree(SL->param_offset);
    SL->param_offset = NULL;
    if (SL->param_index) tfree(SL->param_index);
    SL->param_index = NULL;
  }

  if (!name && !element_type)
    bomb("NULL name and element_type passed to add_steer_elem_to_list", NULL);
  if (name) {
    str_toupper(name);
    if (has_wildcards(name) && strchr(name, '-'))
      name = expand_ranges(name);
  }
  else 
    cp_str(&name, "*");
  if (element_type) {
    str_toupper(element_type);
    if (has_wildcards(element_type) && strchr(element_type, '-'))
      element_type = expand_ranges(element_type);
  }

  if (!item)
    bomb("NULL item passed to add_steer_elem_to_list", NULL);
  str_toupper(item);

  context = NULL;
  found = 0;
  
  while ((context=wfind_element(name, &context, &(beamline->elem)))) {
    if (element_type &&
        !wild_match(entity_name[context->type], element_type))
      continue;

    for (i=0; i<SL->n_corr_types; i++)
      if (strcmp(context->name, SL->corr_name[i])==0)
        return;

    SL->corr_name    = trealloc(SL->corr_name, (SL->n_corr_types+1)*sizeof(*SL->corr_name));
    SL->corr_type    = trealloc(SL->corr_type, (SL->n_corr_types+1)*sizeof(*SL->corr_type));
    SL->corr_param   = trealloc(SL->corr_param, (SL->n_corr_types+1)*sizeof(*SL->corr_param));
    SL->corr_tweek   = trealloc(SL->corr_tweek, (SL->n_corr_types+1)*sizeof(*SL->corr_tweek));
    SL->corr_limit   = trealloc(SL->corr_limit, (SL->n_corr_types+1)*sizeof(*SL->corr_limit));
    SL->param_offset = trealloc(SL->param_offset, (SL->n_corr_types+1)*sizeof(*SL->param_offset));
    SL->param_index  = trealloc(SL->param_index , (SL->n_corr_types+1)*sizeof(*SL->param_index ));
    
    SL->corr_type[SL->n_corr_types] = context->type;
    cp_str(SL->corr_name+SL->n_corr_types, context->name);
    
    cp_str(SL->corr_param+SL->n_corr_types, item);
    SL->corr_tweek[SL->n_corr_types] = tweek;
    SL->corr_limit[SL->n_corr_types] = limit;
    
    if ((SL->param_index[SL->n_corr_types]=param_number=confirm_parameter(item, context->type))<0 ||
        entity_description[context->type].parameter[param_number].type!=IS_DOUBLE ||
        (SL->param_offset[SL->n_corr_types]=find_parameter_offset(item, SL->corr_type[SL->n_corr_types]))<0) {
      fprintf(stderr, "No such floating-point parameter (%s) for %s (add_steer_elem_to_lists)\n", 
              item, context->name);
      exit(1);
    }
    
    SL->n_corr_types += 1;    
    found = 1;
  }
  if (!found)
    bomb("no match to give name or element type", NULL);
}


double compute_kick_coefficient(ELEMENT_LIST *elem, long plane, long type, double corr_tweek, char *name, char *item, RUN *run)
{
  double value, coef;
  long param_offset=0, param_number;

  if ((param_number=confirm_parameter(item, type))<0 || (param_offset=find_parameter_offset(item, type))<0)
    bomb("invalid parameter or element type (compute_kick_coefficient)", NULL);

  switch (type) {
  case T_HCOR:
    if (plane==0 && param_offset==find_parameter_offset("KICK", type))
      return(1.0);
    break;
  case T_VCOR:
    if (plane!=0 && param_offset==find_parameter_offset("KICK", type))
      return(1.0);
    break;
  case T_HVCOR:
    if (plane==0 && param_offset==find_parameter_offset("HKICK", type))
      return(1.0);
    if (plane!=0 && param_offset==find_parameter_offset("VKICK", type))
      return(1.0);
    break;
  default:
    if (entity_description[elem->type].flags&HAS_MATRIX) {
      VMATRIX *M1, *M2, *M0;
      M0 = elem->matrix;
      value = *((double*)(elem->p_elem+param_offset));
      *((double*)(elem->p_elem+param_offset)) += corr_tweek;
      M1 = compute_matrix(elem, run, NULL);
      *((double*)(elem->p_elem+param_offset)) -= 2*corr_tweek;
      M2 = compute_matrix(elem, run, NULL);
      if (plane==0) 
        coef = (M1->C[1]-M2->C[1])/(2*corr_tweek);
      else
        coef = (M1->C[3]-M2->C[3])/(2*corr_tweek);
      /*
        fprintf(stdout, "computed kick coefficient for %s.%s: %g rad/%s\n",
        name, item, coef, entity_description[type].parameter[param_number].unit);
        fflush(stdout);
        */
      free_matrices(M1); tfree(M1); M1 = NULL;
      free_matrices(M2); tfree(M2); M2 = NULL;
      *((double*)(elem->p_elem+param_offset)) = value;
      elem->matrix = M0;
      return(coef);
    }
    else {
      fprintf(stdout, "error: element %s does not have a matrix--can't be used for steering.\n",
              elem->name);
      fflush(stdout);
      exit(1);
    }
    break;
  }
  return(0.0);
}

long do_correction(CORRECTION *correct, RUN *run, LINE_LIST *beamline, double *starting_coord, 
                   BEAM *beam, long sim_step, long initial_correction)
{
  long i, i_cycle, x_failed, y_failed, n_iter_taken, bombed=0, final_traj;
  double *closed_orbit, rms_before, rms_after, *Cdp;

  log_entry("do_correction");

  if (correct->response_only)
    return 1;

  closed_orbit = starting_coord;   /* for return of closed orbit at starting point */

  if (correct->verbose && correct->n_iterations>=1 && correct->n_xy_cycles>0) {
    if (correct->CMx->fixed_length && correct->mode==ORBIT_CORRECTION) {
      fputs(" plane   iter     rms kick     rms pos.    max kick     max pos.   mom.offset\n", stdout);
      fputs("                    mrad          mm         mrad          mm           %\n", stdout);
      fputs("-------  ----    ----------   ---------   ----------   ----------  ----------\n", stdout);
    }
    else {
      fputs(" plane   iter     rms kick     rms pos.    max kick     max pos.\n", stdout);
      fputs("                    mrad          mm         mrad          mm\n", stdout);
      fputs("-------  ----    ----------   ---------   ----------   ----------\n", stdout);
    }
  }

  if (correct->prezero_correctors && initial_correction) {
    long changed;
    if (beamline->elem_recirc)
      changed = zero_correctors(beamline->elem_recirc, run, correct);
    else
      changed = zero_correctors(&(beamline->elem), run, correct);
    if (changed && beamline->matrix) {
      free_matrices(beamline->matrix);
      tfree(beamline->matrix);
      beamline->matrix = NULL;
    }
  }

  correct->CMx->n_cycles_done = correct->CMy->n_cycles_done = 0;
  final_traj = 1;
  switch (correct->mode) {
  case TRAJECTORY_CORRECTION:
    x_failed = y_failed = 0;
    if (usePerturbedMatrix) {
      compute_trajcor_matrices(correct->CMx, &correct->SLx, 0, run, beamline, 0,
                               !correct->response_only);
      compute_trajcor_matrices(correct->CMy, &correct->SLy, 0, run, beamline, 0,
                               !correct->response_only);
    }
    for (i_cycle=0; i_cycle<correct->n_xy_cycles; i_cycle++) {
      final_traj = 1;
      if (!x_failed && correct->CMx->ncor) {
        if (correct->method==GLOBAL_CORRECTION) {
          if (!global_trajcor_plane(correct->CMx, &correct->SLx, 0, correct->traj, correct->n_iterations, 
                                    run, beamline, starting_coord, (correct->use_actual_beam?beam:NULL)))
            return 0;
        } else
          one_to_one_trajcor_plane(correct->CMx, &correct->SLx, 0, correct->traj, correct->n_iterations, 
                                   run, beamline, starting_coord, (correct->use_actual_beam?beam:NULL));
        correct->CMx->n_cycles_done = i_cycle+1;
        if (correct->n_iterations<1)
          break;
        rms_before = rms_value(correct->CMx->posi[0], correct->CMx->nmon);
        rms_after  = rms_value(correct->CMx->posi[correct->n_iterations], correct->CMx->nmon);
#if defined(IEEE_MATH)
        if (isnan(rms_before) || isnan(rms_after) || isinf(rms_before) || isinf(rms_after)) {
          x_failed = 1;
          if (correct->verbose)
            fputs("horizontal trajectory diverged--setting correctors to zero\n", stdout);
          zero_hcorrectors(&(beamline->elem), run, correct);
        }
        else
#endif
          if (rms_before<=rms_after+correct->CMx->corr_accuracy) {
            x_failed = 1;
            if (correct->verbose)
              fputs("trajectory not improved--discontinuing horizontal correction\n", stdout);
          }
        dump_cormon_stats(correct->verbose, 0, correct->CMx->kick, 
                          correct->CMx->ncor, correct->CMx->posi, correct->CMx->nmon, NULL, 
                          correct->n_iterations, i_cycle, i_cycle==correct->n_xy_cycles-1 || x_failed,
                          sim_step, initial_correction);
        if (!initial_correction && (i_cycle==correct->n_xy_cycles-1 || x_failed))
          dump_corrector_data(correct->CMx, &correct->SLx, correct->n_iterations, "horizontal", sim_step);
      }
      if (!y_failed && correct->CMy->ncor) {                    
        final_traj = 2;
        if (correct->method==GLOBAL_CORRECTION) {
          if (!global_trajcor_plane(correct->CMy, &correct->SLy, 2, correct->traj+1, correct->n_iterations, 
                                    run, beamline, starting_coord, (correct->use_actual_beam?beam:NULL)))
            return 0;
        }
        else
          one_to_one_trajcor_plane(correct->CMy, &correct->SLy, 2, correct->traj+1, correct->n_iterations, 
                                   run, beamline, starting_coord, (correct->use_actual_beam?beam:NULL));
        correct->CMy->n_cycles_done = i_cycle+1;
        rms_before = rms_value(correct->CMy->posi[0], correct->CMy->nmon);
        rms_after  = rms_value(correct->CMy->posi[correct->n_iterations], correct->CMy->nmon);
#if defined(IEEE_MATH)
        if (isnan(rms_before) || isnan(rms_after) || isinf(rms_before) || isinf(rms_after)) {
          y_failed = 1;
          if (correct->verbose)
            fputs("vertical trajectory diverged--setting correctors to zero\n", stdout);
          zero_vcorrectors(&(beamline->elem), run, correct);
        }
        else
#endif
          if (rms_before<=rms_after+correct->CMy->corr_accuracy) {
            y_failed = 1;
            if (correct->verbose)
              fputs("trajectory not improved--discontinuing vertical correction\n", stdout);
          }
        dump_cormon_stats(correct->verbose, 2, correct->CMy->kick, 
                          correct->CMy->ncor, correct->CMy->posi, correct->CMy->nmon, NULL, 
                          correct->n_iterations, i_cycle, i_cycle==correct->n_xy_cycles-1 || y_failed,
                          sim_step, initial_correction);
        if (!initial_correction && (i_cycle==correct->n_xy_cycles-1 || y_failed))
          dump_corrector_data(correct->CMy, &correct->SLy, correct->n_iterations, "vertical", sim_step);
      }
      if (initial_correction && i_cycle==0)
        dump_orb_traj(correct->traj[0], beamline->n_elems, "uncorrected", sim_step);
      if (x_failed && y_failed) {
        if (correct->verbose)
          fputs("trajectory correction discontinued\n", stdout);
        break;
      }
    }
    if (!initial_correction && correct->n_iterations>=1)
      dump_orb_traj(correct->traj[final_traj], beamline->n_elems, "corrected", sim_step);
    if (starting_coord)
      for (i=0; i<6; i++)
        starting_coord[i] = 0;  /* don't want to seem to be returning a closed orbit here */
    bombed = 0;
    break;
  case ORBIT_CORRECTION:
    if (correct->CMx->fixed_length)
      Cdp = tmalloc(sizeof(*Cdp)*(correct->n_iterations+1));
    else
      Cdp = NULL;
    x_failed = y_failed = bombed = 0;
    if (usePerturbedMatrix) {
      compute_orbcor_matrices(correct->CMx, &correct->SLx, 0, run, beamline, 0, 
                              !correct->response_only, fixedLengthMatrix);
      compute_orbcor_matrices(correct->CMy, &correct->SLy, 2, run, beamline, 0, 
                              !correct->response_only, fixedLengthMatrix);
    }
    for (i_cycle=0; i_cycle<correct->n_xy_cycles; i_cycle++) {
      final_traj = 1;
      if (!x_failed && correct->CMx->ncor) {
        if ((n_iter_taken = orbcor_plane(correct->CMx, &correct->SLx, 0, correct->traj, 
                                         correct->n_iterations, correct->clorb_accuracy, 
                                         correct->clorb_iterations, 
                                         correct->clorb_iter_fraction,
                                         run, beamline, 
                                         closed_orbit, Cdp))<0) {
          fprintf(stdout, "Horizontal correction has failed.\n");
          fflush(stdout);
          bombed = 1; 
          break;
        }
        correct->CMx->n_cycles_done = i_cycle+1;
        if (correct->n_iterations<1)
          break;
        rms_before = rms_value(correct->CMx->posi[0], correct->CMx->nmon);
        rms_after  = rms_value(correct->CMx->posi[n_iter_taken], correct->CMx->nmon);
#if defined(IEEE_MATH)
        if (isnan(rms_before) || isnan(rms_after) || isinf(rms_before) || isinf(rms_after)) {
          x_failed = 1;
          if (correct->verbose)
            fputs("horizontal orbit diverged--setting correctors to zero\n", stdout);
          zero_hcorrectors(&(beamline->elem), run, correct);
        }
        else
#endif
          if (rms_before<=rms_after+correct->CMx->corr_accuracy) {
            x_failed = 1;
            if (correct->verbose)
              fputs("orbit not improved--discontinuing horizontal correction\n", stdout);
          }
        dump_cormon_stats(correct->verbose, 0, correct->CMx->kick, 
                          correct->CMx->ncor, correct->CMx->posi, correct->CMx->nmon, Cdp, 
                          correct->n_iterations, i_cycle, i_cycle==correct->n_xy_cycles-1 || x_failed,
                          sim_step, initial_correction);
        if (!initial_correction && (i_cycle==correct->n_xy_cycles-1 || x_failed)) 
          dump_corrector_data(correct->CMx, &correct->SLx, correct->n_iterations, "horizontal", sim_step);
      }
      if (!y_failed && correct->CMy->ncor) {
        final_traj = 2;
        if ((n_iter_taken = orbcor_plane(correct->CMy, &correct->SLy, 2, correct->traj+1, 
                                         correct->n_iterations, correct->clorb_accuracy, 
                                         correct->clorb_iterations, 
                                         correct->clorb_iter_fraction,
                                         run, beamline, 
                                         closed_orbit, Cdp))<0) {
          bombed = 1;
          fprintf(stdout, "Vertical correction has failed.\n");
          fflush(stdout);
          break;
        }
        correct->CMy->n_cycles_done = i_cycle+1;
        rms_before = rms_value(correct->CMy->posi[0], correct->CMy->nmon);
        rms_after  = rms_value(correct->CMy->posi[n_iter_taken], correct->CMy->nmon);
#if defined(IEEE_MATH)
        if (isnan(rms_before) || isnan(rms_after) || isinf(rms_before) || isinf(rms_after)) {
          y_failed = 1;
          if (correct->verbose)
            fputs("vertical orbit diverged--setting correctors to zero\n", stdout);
          zero_vcorrectors(&(beamline->elem), run, correct);
        }
        else
#endif
          if (rms_before<=rms_after+correct->CMx->corr_accuracy) {
            y_failed = 1;
            if (correct->verbose)
              fputs("orbit not improved--discontinuing vertical correction\n", stdout);
          }
        dump_cormon_stats(correct->verbose, 2, correct->CMy->kick, 
                          correct->CMy->ncor, correct->CMy->posi, correct->CMy->nmon, Cdp, 
                          correct->n_iterations, i_cycle, i_cycle==correct->n_xy_cycles-1 || y_failed,
                          sim_step, initial_correction);
        if (!initial_correction && (i_cycle==correct->n_xy_cycles-1 || y_failed))
          dump_corrector_data(correct->CMy, &correct->SLy, correct->n_iterations, "vertical", sim_step);
      }
      if (initial_correction && i_cycle==0) 
        dump_orb_traj(correct->traj[0], beamline->n_elems, "uncorrected", sim_step);
      if (x_failed && y_failed) {
        if (correct->verbose)
          fputs("orbit correction discontinued\n", stdout);
        break;
      }
    }
    if (!initial_correction && !bombed && correct->n_iterations>=1)
      dump_orb_traj(correct->traj[final_traj], beamline->n_elems, "corrected", sim_step);
    break;
  }

  beamline->closed_orbit = correct->traj[final_traj];

  log_exit("do_correction");
  return(!bombed);
}

void compute_trajcor_matrices(CORMON_DATA *CM, STEERING_LIST *SL, long coord, RUN *run, LINE_LIST *beamline, long find_only, long invert)
{
#ifdef DEBUG
  long i_debug;
  FILE *fpdeb;
  char s[100];
#endif
  ELEMENT_LIST *corr, *start;
  TRAJECTORY *traj0, *traj1;
  long kick_offset, i_corr, i_moni, i, equalW;
  long n_part;
  static MATRIX *I1=NULL, *I2=NULL, *I3=NULL, *I4=NULL, *W=NULL;
  double **one_part, p, p0, kick0, kick1, corr_tweek, corrCalibration, *moniCalibration, W0=0.0;
  VMATRIX *save;
  long i_type;

  log_entry("compute_trajcor_matrices");

  start = find_useable_moni_corr(&CM->nmon, &CM->ncor, &CM->mon_index,
                                 &CM->umoni, &CM->ucorr, &CM->kick_coef, &CM->sl_index, coord, SL, run, beamline, 0);
  if (CM->nmon<CM->ncor) {
    fprintf(stdout, "*** Warning: more correctors than monitors for %c plane.\n",  (coord==0?'x':'y'));
    fprintf(stdout, "*** Correction will probably be unstable!\n");
    fflush(stdout);
  }
  if (CM->ncor==0) {
    fprintf(stdout, "Warning: no correctors for %c plane.  No correction done.\n",  (coord==0?'x':'y'));
    fflush(stdout);
    return;
  }
  if (CM->nmon==0) {
    fprintf(stdout, "Warning: no monitors for %c plane.  No correction done.\n",  (coord==0?'x':'y'));
    fflush(stdout);
    CM->ncor = 0;
    return;
  }

  if (find_only) {
    log_exit("compute_trajcor_matrices");
    return;
  }

  fprintf(stdout, "computing response matrix...\n");
  fflush(stdout);
  report_stats(stdout, "start");

  /* allocate correction matrix for this plane, plus others: dK = T*Qo */
  m_alloc1(&CM->T , CM->ncor, CM->nmon);
  m_alloc1(&CM->Qo, CM->nmon, 1);
  m_alloc1(&CM->dK, CM->ncor, 1);
  m_alloc1(&CM->C , CM->nmon, CM->ncor);

  /* intermediate matrices for computations: T = -I4*I2 */
  m_alloc1(&W , CM->nmon, CM->nmon);
  m_alloc1(&I1, CM->ncor, CM->nmon);        /* I1 = TRANS(C) */
  m_alloc1(&I2, CM->ncor, CM->nmon);        /* I2 = TRANS(C).W */
  m_alloc1(&I3, CM->ncor, CM->ncor);        /* I3 = TRANS(C).W.C */
  m_alloc1(&I4, CM->ncor, CM->ncor);        /* I4 = INVERSE(I3) */

  /* arrays for trajectory data */
  traj0 = tmalloc(sizeof(*traj0)*beamline->n_elems);
  traj1 = tmalloc(sizeof(*traj1)*beamline->n_elems);
  one_part = (double**)zarray_2d(sizeof(**one_part), 1, 7);

  /* find initial trajectory */
  p = p0 = sqrt(sqr(run->ideal_gamma)-1);
  n_part = 1;
  fill_double_array(*one_part, 7, 0.0);
  if (!do_tracking(one_part, &n_part, NULL, beamline, &p, (double**)NULL, (BEAM_SUMS**)NULL, (long*)NULL,
                   traj0, run, 0, TEST_PARTICLES+TIME_DEPENDENCE_OFF, 1, 0, NULL, NULL, NULL))
    bomb("tracking failed for test particle (compute_trajcor_matrices())", NULL);

#if  DEBUG
  i_debug = 0;
  sprintf(s, "traj%c-%ld.deb", (coord==0?'x':'y'), i_debug++);
  fpdeb = fopen_e(s, "w", 0);
  fprintf(fpdeb, "z (m)\n%c (m)\ninitial trajectory computed in compute_trajcor_matrices\n\n%ld\n", 
          (coord==0?'x':'y'), beamline->n_elems);
  for (i=0; i<beamline->n_elems; i++)
    fprintf(fpdeb, "%e\t%e\n", traj0[i].elem->end_pos, traj0[i].centroid[coord]);
  fclose(fpdeb);
#endif

  /* set up weight matrix and monitor calibration array */
  moniCalibration = tmalloc(sizeof(*moniCalibration)*CM->nmon);
  equalW = 1;
  for (i_moni=0; i_moni<CM->nmon; i_moni++) {
    long j;
    for (j=0; j<CM->nmon; j++) 
      W->a[i_moni][j] = 0.0;
    W->a[i_moni][i_moni] = getMonitorWeight(CM->umoni[i_moni]);
    if (!i_moni)
      W0 = W->a[i_moni][i_moni];
    else if (W0!=W->a[i_moni][i_moni])
      equalW = 0;
    moniCalibration[i_moni] = getMonitorCalibration(CM->umoni[i_moni], coord);
  }

  for (i_corr = 0; i_corr<CM->ncor; i_corr++) {
    corr = CM->ucorr[i_corr];

    if ((i_type = find_index(corr->type, SL->corr_type, SL->n_corr_types))<0)
      bomb("failed to find corrector type in type list", NULL);

    kick_offset = SL->param_offset[i_type];
    corr_tweek  = SL->corr_tweek[i_type];

    /* record value of corrector */
    kick0 = *((double*)(corr->p_elem+kick_offset));

    /* change the corrector by corr_tweek and compute the new matrix for the corrector */
    kick1 = *((double*)(corr->p_elem+kick_offset)) = kick0 + corr_tweek;
#ifdef DEBUG
    fprintf(stdout, "corrector %s tweeked to %e\n", corr->name, *((double*)(corr->p_elem+kick_offset)));
    fflush(stdout);
#endif

    if (corr->matrix)
      save = corr->matrix;
    else
      save = NULL;
    compute_matrix(corr, run, NULL);
#ifdef DEBUG
    print_matrices(stdout, "*** corrector matrix:", corr->matrix);
#endif

    /* track with positively-tweeked corrector */
    p = p0;
    n_part = 1;
    fill_double_array(*one_part, 7, 0.0);
    if (!do_tracking(one_part, &n_part, NULL, beamline, &p, (double**)NULL, (BEAM_SUMS**)NULL, (long*)NULL,
                     traj1, run, 0, TEST_PARTICLES+TIME_DEPENDENCE_OFF, 1, 0, NULL, NULL, NULL))
      bomb("tracking failed for test particle (compute_trajcor_matrices())", NULL);

#ifdef DEBUG
    sprintf(s, "traj%c-%ld.deb", (coord==0?'x':'y'), i_debug++);
    fpdeb = fopen_e(s, "w", 0);
    fprintf(fpdeb, "z (m)\n%c (m)\ntrajectory computed in compute_trajcor_matrices\n%s = %e\n%ld\n", 
            (coord==0?'x':'y'), corr->name, kick1, beamline->n_elems);
    for (i=0; i<beamline->n_elems; i++)
      fprintf(fpdeb, "%e\t%e\n", traj1[i].elem->end_pos, traj1[i].centroid[coord]);
    fclose(fpdeb);
#endif


#if TWO_POINT_TRAJRESPONSE
    /* change the corrector by -corr_tweek and compute the new matrix for the corrector */
    kick1 = *((double*)(corr->p_elem+kick_offset)) = kick0 - corr_tweek;
    if (beamline->links)
      assert_element_links(beamline->links, run, beamline, DYNAMIC_LINK);
#ifdef DEBUG
    fprintf(stdout, "corrector %s tweeked to %e\n", corr->name, *((double*)(corr->p_elem+kick_offset)));
    fflush(stdout);
#endif
    free_matrices(corr->matrix); corr->matrix = NULL;
    compute_matrix(corr, run, NULL);
#ifdef DEBUG
    print_matrices(stdout, "*** corrector matrix:", corr->matrix);
#endif

    /* track with tweeked corrector */
    p = p0;
    n_part = 1;
    fill_double_array(*one_part, 7, 0.0);
    if (!do_tracking(one_part, &n_part, NULL, beamline, &p, (double**)NULL, (BEAM_SUMS**)NULL, (long*)NULL,
                     traj0, run, 0, TEST_PARTICLES+TIME_DEPENDENCE_OFF, 1, 0, NULL, NULL, NULL))
      bomb("tracking failed for test particle (compute_trajcor_matrices())", NULL);

    /* compute coefficients of array C that are driven by this corrector */
    corrCalibration = getCorrectorCalibration(CM->ucorr[i_corr], coord)/(2*corr_tweek);
    for (i_moni=0; i_moni<CM->nmon; i_moni++) {
      i = CM->mon_index[i_moni];
      CM->C->a[i_moni][i_corr] = moniCalibration[i_moni]*corrCalibration*
        (traj1[i].centroid[coord] - traj0[i].centroid[coord]);
    }
#else

    /* compute coefficients of array C that are driven by this corrector */
    corrCalibration = getCorrectorCalibration(CM->ucorr[i_corr], coord)/corr_tweek;
    for (i_moni=0; i_moni<CM->nmon; i_moni++) {
      i = CM->mon_index[i_moni];
      CM->C->a[i_moni][i_corr] = moniCalibration[i_moni]*corrCalibration*
        (traj1[i].centroid[coord] - traj0[i].centroid[coord]);
    }
#endif

    /* change the corrector back */
    *((double*)(corr->p_elem+kick_offset)) = kick0;
    if (beamline->links)
      assert_element_links(beamline->links, run, beamline, DYNAMIC_LINK);
    if (corr->matrix) {
      free_matrices(corr->matrix);
      corr->matrix = NULL;
    }
    if (save)
      corr->matrix = save;
    else 
      compute_matrix(corr, run, NULL);
  } 
  free(moniCalibration);
  tfree(traj0); traj0 = NULL;
  tfree(traj1); traj1 = NULL;
  free_zarray_2d((void**)one_part, 1, 7); one_part = NULL;
#ifdef DEBUG
  m_show(CM->C    , "%13.6le ", "influence matrix\n", stdout);
#endif
  report_stats(stdout, "done");

  if (invert && !CM->inverse_computed) {
    fprintf(stdout, "computing correction matrix...\n");
    fflush(stdout);
    report_stats(stdout, "start");
    /* compute correction matrix T */
    CM->inverse_computed = 1;
    if (equalW) {
      m_trans(I1, CM->C);
      m_mult(I3, I1, CM->C);
      m_invert(I4, I3);
      m_mult(CM->T, I4, I1);
      m_scmul(CM->T, CM->T, -W0); 
    }
    else {
      m_trans(I1, CM->C);
      m_mult(I2, I1, W);
      m_mult(I3, I2, CM->C);
      m_invert(I4, I3);
      m_mult(CM->T, I4, I2);
      m_scmul(CM->T, CM->T, -1.0); 
    }
    report_stats(stdout, "done");
#ifdef DEBUG
    m_show(CM->T, "%13.6le ", "correction matrix\n", stdout);
#endif
  }

  log_exit("compute_trajcor_matrices");
}

long global_trajcor_plane(CORMON_DATA *CM, STEERING_LIST *SL, long coord, TRAJECTORY **traject, long n_iterations, 
                          RUN *run, LINE_LIST *beamline, double *starting_coord, BEAM *beam)
{
  ELEMENT_LIST *corr, *eptr;
  TRAJECTORY *traj;
  long iteration, kick_offset;
  long i_moni, i_corr;
  long n_part, i, tracking_flags, sl_index;
  double **particle;
  double p, x, y, reading, fraction, minFraction, param, change;

  log_entry("global_trajcor_plane");

  if (!m_check(CM->T) || !m_check(CM->Qo) || !m_check(CM->dK))
    bomb("corrupted correction matrix detected (global_trajcor_plane)", NULL);
  if (!CM->mon_index)
    bomb("monitor index array is NULL (global_trajcor_plane)", NULL);
  if (!CM->posi)
    bomb("monitor readout array is NULL (global_trajcor_plane)", NULL);
  if (!CM->kick)
    bomb("corrector value array is NULL (global_trajcor_plane)", NULL);
  if (!traject)
    bomb("no trajectory arrays supplied (global_trajcor_plane)", NULL);
  if (!CM->inverse_computed)
    bomb("no inverse matrix computed (global_trajcor_plane)", NULL);

  if (!beam) {
    particle = (double**)zarray_2d(sizeof(**particle), 1, 7);
    tracking_flags = TEST_PARTICLES;
  }
  else {
    if (beam->n_to_track==0)
      bomb("no particles to track in global_trajcor_plane()", NULL);
    particle = (double**)zarray_2d(sizeof(**particle), beam->n_to_track, 7);
    tracking_flags = TEST_PARTICLES+TEST_PARTICLE_LOSSES;
  }

  if (CM->nmon<CM->ncor) {
    fprintf(stdout, "*** Warning: more correctors than monitors for %c plane.\n",  (coord==0?'x':'y'));
    fprintf(stdout, "*** Correction will probably be unstable!\n");
    fflush(stdout);
  }
  for (iteration=0; iteration<=n_iterations; iteration++) {
    if (!CM->posi[iteration])
      bomb("monitor readout array for this iteration is NULL (global_trajcor_plane)", NULL);
    if (!CM->kick[iteration])
      bomb("corrector value array for this iteration is NULL (global_trajcor_plane)", NULL);
    if (!(traj = traject[iteration?1:0]))
      bomb("trajectory array for this iteration is NULL (global_trajcor_plane)", NULL);

    /* find trajectory */
    p = sqrt(sqr(run->ideal_gamma)-1);

    if (!beam) {
      if (!starting_coord)
        fill_double_array(*particle, 7, 0.0);
      else {
        for (i=0; i<6; i++)
          particle[0][i] = starting_coord[i];
      }
      n_part = 1;
    }
    else
      copy_particles(particle, beam->particle, n_part=beam->n_to_track);

    n_part = do_tracking(particle, &n_part, NULL, beamline, &p, (double**)NULL, 
                         (BEAM_SUMS**)NULL, (long*)NULL,
                         traj, run, 0, tracking_flags, 1, 0, NULL, NULL, NULL);
    if (beam) {
      fprintf(stdout, "%ld particles survived tracking", n_part);
      fflush(stdout);
      if (n_part==0) {
        for (i=0; i<beamline->n_elems+1; i++)
          if (traj[i].n_part==0)
            break;
        if (i!=0 && i<beamline->n_elems+1)
          fprintf(stdout, "---all beam lost before z=%em (element %s)",
                  traj[i].elem->end_pos, traj[i].elem->name);
        fflush(stdout);
        fputc('\n', stdout);
        return 0;
      }
      fputc('\n', stdout);
    }

    /* find readings at monitors and add in reading errors */
    for (i_moni=0; i_moni<CM->nmon; i_moni++) {
      if (!(eptr=traj[CM->mon_index[i_moni]].elem))
        bomb("invalid element pointer in trajectory array (global_trajcor_plane)", NULL);
      x = traj[CM->mon_index[i_moni]].centroid[0];
      y = traj[CM->mon_index[i_moni]].centroid[2];
      reading = computeMonitorReading(eptr, coord, x, y, 0);
      if (isnan(reading) || isinf(reading)) 
        return 0;
      CM->posi[iteration][i_moni] = CM->Qo->a[i_moni][0] =  
        reading + (CM->bpm_noise?noise_value(CM->bpm_noise, CM->bpm_noise_cutoff, CM->bpm_noise_distribution):0);
    }
    
    if (iteration==n_iterations)
      break;

    /* solve for the corrector changes */
    m_mult(CM->dK, CM->T, CM->Qo);

#ifdef DEBUG
    m_show(CM->Qo, "%13.6le ", "traj matrix\n", stdout);
    m_show(CM->dK, "%13.6le ", "kick matrix\n", stdout);
#endif

    /* step through beamline find any kicks that are over their limits */
    minFraction = 1;
    for (i_corr=0; i_corr<CM->ncor; i_corr++) {
      corr = CM->ucorr[i_corr];
      sl_index = CM->sl_index[i_corr];
      kick_offset = SL->param_offset[sl_index];
      param = fabs(*((double*)(corr->p_elem+kick_offset)) +
                   (change=CM->dK->a[i_corr][0]/CM->kick_coef[i_corr]*CM->corr_fraction));
      if (SL->corr_limit[sl_index] && param>SL->corr_limit[sl_index]) {
        fraction = fabs((SL->corr_limit[sl_index]-fabs(*((double*)(corr->p_elem+kick_offset))))/change);
        if (fraction<minFraction)
          minFraction = fraction;
      }
    }
    fraction = minFraction*CM->corr_fraction;

#if defined(DEBUG)
    fprintf(stdout, "Changing correctors:");
    fflush(stdout);
#endif
    /* step through beamline and change correctors */
    for (i_corr=0; i_corr<CM->ncor; i_corr++) {
      corr = CM->ucorr[i_corr];
      sl_index = CM->sl_index[i_corr];
      kick_offset = SL->param_offset[sl_index];
#if defined(DEBUG)
      fprintf(stdout, "before = %e, ", *((double*)(corr->p_elem+kick_offset)));
      fflush(stdout);
#endif
      if (iteration==0)
        CM->kick[iteration][i_corr] = *((double*)(corr->p_elem+kick_offset))*CM->kick_coef[i_corr];
      *((double*)(corr->p_elem+kick_offset)) += CM->dK->a[i_corr][0]*fraction;
      CM->kick[iteration+1][i_corr] = *((double*)(corr->p_elem+kick_offset))*CM->kick_coef[i_corr];
#if defined(DEBUG)
      fprintf(stdout, "after = %e\n", *((double*)(corr->p_elem+kick_offset)));
      fflush(stdout);
#endif
      if (corr->matrix) {
        free_matrices(corr->matrix);
        tfree(corr->matrix);
        corr->matrix = NULL;
      }
      compute_matrix(corr, run, NULL);
    }
    if (beamline->links)
      assert_element_links(beamline->links, run, beamline, DYNAMIC_LINK);
  }

  /* indicate that beamline concatenation and Twiss parameter computation (if wanted) are not current */
  beamline->flags &= ~BEAMLINE_CONCAT_CURRENT;
  beamline->flags &= ~BEAMLINE_TWISS_CURRENT;

  if (!beam)
    free_zarray_2d((void**)particle, 1, 7);
  else
    free_zarray_2d((void**)particle, beam->n_to_track, 7);
  particle = NULL;
  log_exit("global_trajcor_plane");
  return(1);
}

void one_to_one_trajcor_plane(CORMON_DATA *CM, STEERING_LIST *SL, long coord, TRAJECTORY **traject, long n_iterations, RUN *run,
                              LINE_LIST *beamline, double *starting_coord, BEAM *beam)
{
  ELEMENT_LIST *corr, *eptr;
  TRAJECTORY *traj;
  long iteration, kick_offset;
  long i_moni, i_corr, sl_index;
  long n_part, i, tracking_flags;
  double **particle, param, fraction;
  double p, x, y, reading;
  
  log_entry("one_to_one_trajcor_plane");
  
  if (!m_check(CM->T) || !m_check(CM->Qo) || !m_check(CM->dK))
    bomb("corrupted correction matrix detected (one_to_one_trajcor_plane)", NULL);
  if (!CM->mon_index)
    bomb("monitor index array is NULL (one_to_one_trajcor_plane)", NULL);
  if (!CM->posi)
    bomb("monitor readout array is NULL (one_to_one_trajcor_plane)", NULL);
  if (!CM->kick)
    bomb("corrector value array is NULL (one_to_one_trajcor_plane)", NULL);
  if (!traject)
    bomb("no trajectory arrays supplied (one_to_one_trajcor_plane)", NULL);
  
  if (!beam) {
    particle = (double**)zarray_2d(sizeof(**particle), 1, 7);
    tracking_flags = TEST_PARTICLES;
  }
  else {
    if (beam->n_to_track==0)
      bomb("no particles to track in one_to_one_trajcor_plane()", NULL);
    particle = (double**)zarray_2d(sizeof(**particle), beam->n_to_track, 7);
    tracking_flags = TEST_PARTICLES+TEST_PARTICLE_LOSSES;
  }
  
  for (iteration=0; iteration<=n_iterations; iteration++) {
    if (!CM->posi[iteration])
      bomb("monitor readout array for this iteration is NULL (one_to_one_trajcor_plane)", NULL);
    if (!CM->kick[iteration])
      bomb("corrector value array for this iteration is NULL (one_to_one_trajcor_plane)", NULL);
    if (!(traj = traject[iteration?1:0]))
      bomb("trajectory array for this iteration is NULL (one_to_one_trajcor_plane)", NULL);
    
    for (i_moni=i_corr=0; i_moni<CM->nmon; i_moni++) {
      /* find trajectory */
      p = sqrt(sqr(run->ideal_gamma)-1);
      
      if (!beam) {
        if (!starting_coord)
          fill_double_array(*particle, 7, 0.0);
        else {
          for (i=0; i<7; i++)
            particle[0][i] = starting_coord[i];
        }
        n_part = 1;
      }
      else
        copy_particles(particle, beam->particle, n_part=beam->n_to_track);
      
      n_part = do_tracking(particle, &n_part, NULL, beamline, &p, (double**)NULL, 
                           (BEAM_SUMS**)NULL, (long*)NULL,
                           traj, run, 0, tracking_flags, 1, 0, NULL, NULL, NULL);
      if (beam) {
        fprintf(stdout, "%ld particles survived tracking", n_part);
        fflush(stdout);
        if (n_part==0) {
          for (i=0; i<beamline->n_elems+1; i++)
            if (traj[i].n_part==0)
              break;
          if (i!=0 && i<beamline->n_elems+1)
            fprintf(stdout, "---all beam lost before z=%em (element %s)",
                    traj[i].elem->end_pos, traj[i].elem->name);
          fflush(stdout);
        }
        fputc('\n', stdout);
      }
      
      if (!(eptr=traj[CM->mon_index[i_moni]].elem))
        bomb("invalid element pointer in trajectory array (one_to_one_trajcor_plane)", NULL);
      x = traj[CM->mon_index[i_moni]].centroid[0];
      y = traj[CM->mon_index[i_moni]].centroid[2];
      reading = computeMonitorReading(eptr, coord, x, y, 0);
      reading += (CM->bpm_noise?noise_value(CM->bpm_noise, CM->bpm_noise_cutoff, CM->bpm_noise_distribution):0);
      CM->posi[iteration][i_moni] = reading;

#if defined(DEBUG)
      fprintf(stdout, "i_moni = %ld, i_corr = %ld, reading = %e", i_moni, i_corr, reading);
      fflush(stdout);
#endif
      if (iteration==n_iterations || (i_corr>=CM->ncor || CM->C->a[i_moni][i_corr]==0)) {
#if defined(DEBUG)
        fprintf(stdout, "\n");
        fflush(stdout);
#endif
        continue;
      }

      /* change corrector */
      corr = CM->ucorr[i_corr];
      sl_index = CM->sl_index[i_corr];
      kick_offset = SL->param_offset[sl_index];
      if (iteration==0)
        CM->kick[iteration][i_corr] = *((double*)(corr->p_elem+kick_offset))*CM->kick_coef[i_corr];
      param = *((double*)(corr->p_elem+kick_offset)) - reading/CM->C->a[i_moni][i_corr]*CM->corr_fraction;
      fraction = 1;
      if (param && SL->corr_limit[sl_index])
        if ((fraction = SL->corr_limit[sl_index]/fabs(param))>1)
          fraction = 1;
      *((double*)(corr->p_elem+kick_offset)) = param*fraction;
      CM->kick[iteration+1][i_corr] = *((double*)(corr->p_elem+kick_offset))*CM->kick_coef[i_corr];
#if defined(DEBUG)
      fprintf(stdout, ", param = %e, fraction=%e\n", param, fraction);
      fflush(stdout);
#endif
      if (corr->matrix) {
        free_matrices(corr->matrix);
        tfree(corr->matrix);
        corr->matrix = NULL;
      }
      compute_matrix(corr, run, NULL);
      if (beamline->links)
        assert_element_links(beamline->links, run, beamline, DYNAMIC_LINK);
      
      /* indicate that beamline concatenation and Twiss parameter computation (if wanted) are not current */
      beamline->flags &= ~BEAMLINE_CONCAT_CURRENT;
      beamline->flags &= ~BEAMLINE_TWISS_CURRENT;
      i_corr++;
    }
  }

  if (!beam)
    free_zarray_2d((void**)particle, 1, 7);
  else
    free_zarray_2d((void**)particle, beam->n_to_track, 7);
  particle = NULL;
  log_exit("one_to_one_trajcor_plane");
}


ELEMENT_LIST *find_useable_moni_corr(long *nmon, long *ncor, long **mon_index,
                                     ELEMENT_LIST ***umoni, ELEMENT_LIST ***ucorr,
                                     double **kick_coef, long **sl_index, long plane, STEERING_LIST *SL, 
                                     RUN *run, LINE_LIST *beamline, long recircs)
{
  ELEMENT_LIST *corr, *moni, *start;
  long i_elem, moni_follows, corr_seen, index;
  long moni_type_1, moni_type_2;

  log_entry("find_useable_moni_corr");

  switch (plane) {
  case 0:
    moni_type_1 = T_HMON;
    moni_type_2 = T_MONI;
    break;
  case 2:
    moni_type_1 = T_VMON;
    moni_type_2 = T_MONI;
    break;
  default:
    fprintf(stdout, "error: invalid coordinate for correction: %ld\n", plane);
    fflush(stdout);
    exit(1);
    break;
  }

  start = &(beamline->elem);
  if (recircs) {
    /* advance to recirculation point, if there is one */
    while (start) {
      if (start->type==T_RECIRC)
        break;
      start = start->succ;
    }
    if (!start)
      start = &(beamline->elem);
  }
  if (*ncor && *nmon) {
    log_exit("find_useable_moni_corr");
    return start;
  }

  
  *ncor = *nmon = 0;
  *mon_index = NULL;
  *umoni = NULL;
  *ucorr = NULL;
  *kick_coef = NULL;
  *sl_index = NULL;

  /* first count correctors */
  corr = start;
  do {
    /* advance to position of next corrector */
    if (!(corr = next_element_of_types(corr, SL->corr_type, SL->n_corr_types, &index)))
      break;
    if (steering_corrector(corr, SL) && match_string(corr->name, SL->corr_name, SL->n_corr_types, EXACT_MATCH)>=0) {
      *ucorr = trealloc(*ucorr, sizeof(**ucorr)*(*ncor+1));
      (*ucorr)[*ncor] = corr;
      *sl_index = trealloc(*sl_index, sizeof(**sl_index)*(*ncor+1));
      (*sl_index)[*ncor] = index;
      *kick_coef = trealloc(*kick_coef, sizeof(**kick_coef)*(*ncor+1));
      if (!((*kick_coef)[*ncor] = compute_kick_coefficient(corr, plane, corr->type, 
                                                           SL->corr_tweek[index], corr->name, SL->corr_param[index], run))) {
        fprintf(stdout, "error: changing %s.%s does not kick the beam--can't use for steering.\n",
                corr->name, SL->corr_param[index]);
        fflush(stdout);
        exit(1);
      }
      
      if (!recircs) {
        /* Since the beamline doesn't recirculate, advance through remainder of beamline
         * to make sure there are subsequent monitors */
        moni = corr->succ;
        moni_follows = 0;
        do {
          while (moni && !((moni->type==moni_type_1 || moni->type==moni_type_2) &&
                           ((MONI*)moni->p_elem)->weight>0) )
            moni = moni->succ;
          if (moni)
            moni_follows = 1;
        } while (moni && (moni=moni->succ) && !moni_follows);
        if (!moni_follows)
          break;        /* ignore correctors with no monitors following */
      }
      *ncor += 1;
    }
  } while ((corr=corr->succ));

  if (!recircs) {
    /* count all monitors with one or more correctors upstream and non-negative weight */
    moni = start;
    corr_seen = 0;
    i_elem = 0;
    while (moni) {
      if (find_index(moni->type, SL->corr_type, SL->n_corr_types)!=-1 && 
          (steering_corrector(moni, SL) || match_string(moni->name, SL->corr_name, SL->n_corr_types, EXACT_MATCH)>=0))
        corr_seen = 1;
      if ((moni->type==moni_type_1 || moni->type==moni_type_2) 
          && ((MONI*)moni->p_elem)->weight>0 && corr_seen) {
        *nmon += 1;
        *umoni = trealloc(*umoni, *nmon*sizeof(**umoni));
        (*umoni)[*nmon-1] = moni;
        *mon_index = trealloc(*mon_index, *nmon*sizeof(**mon_index));
        (*mon_index)[*nmon-1] = i_elem;
      }
      moni = moni->succ;
      i_elem++;
    }
  }
  else {
    /* find all monitors */
    moni = start;
    i_elem = 0;
    while (moni) {
      if ((moni->type==moni_type_1 || moni->type==moni_type_2) 
          && ((MONI*)moni->p_elem)->weight>0) {
        *nmon += 1;
        *umoni = trealloc(*umoni, *nmon*sizeof(**umoni));
        (*umoni)[*nmon-1] = moni;
        *mon_index = trealloc(*mon_index, *nmon*sizeof(**mon_index));
        (*mon_index)[*nmon-1] = i_elem;
      }
      moni = moni->succ;
      i_elem++;
    }
  } 
  log_exit("find_useable_moni_corr");
  return(start);
}

void compute_orbcor_matrices(CORMON_DATA *CM, STEERING_LIST *SL, long coord, RUN *run, LINE_LIST *beamline, 
                             long find_only, long invert, long fixed_length)
{
  ELEMENT_LIST *start;
  long i_corr, i_moni, equalW;
  double coef, htune, moniFactor, *corrFactor, *corrFactorFL, coefFL, W0=0.0;
  static MATRIX *I1=NULL, *I2=NULL, *I3=NULL, *I4=NULL, *W=NULL;

  log_entry("compute_orbcor_matrices");

  start = find_useable_moni_corr(&CM->nmon, &CM->ncor, &CM->mon_index, &CM->umoni, &CM->ucorr, 
                                 &CM->kick_coef, &CM->sl_index, coord, SL, run, beamline, 1);
#ifdef DEBUG
  fprintf(stdout, "finding twiss parameters beginning at %s.\n", start->name);
  fflush(stdout);
#endif
  if (!(beamline->flags&BEAMLINE_TWISS_CURRENT)) {
    fprintf(stdout, "updating twiss parameters...");
    fflush(stdout);
    update_twiss_parameters(run, beamline, NULL);
#ifdef DEBUG
    fprintf(stdout, "Tunes: %e, %e\n", beamline->tune[0], beamline->tune[1]);
    fflush(stdout);
    fprintf(stdout, "Initial eta: %e, %e\n", 
            beamline->elem.twiss->etax, 
            beamline->elem.twiss->etay);
    fflush(stdout);
    fprintf(stdout, "Final eta: %e, %e\n", 
            beamline->elast->twiss->etax, 
            beamline->elast->twiss->etay);
    fflush(stdout);
#endif
    report_stats(stdout, "\ndone: ");
  }

#ifdef DEBUG
  fprintf(stdout, "monitors: ");
  fflush(stdout);
  for (i_moni=0; i_moni<CM->nmon; i_moni++)
    fprintf(stdout, "%s ", CM->umoni[i_moni]->name);
  fflush(stdout);
  fprintf(stdout, "\ncorrectors: ");
  fflush(stdout);
  for (i_corr=0; i_corr<CM->ncor; i_corr++)
    fprintf(stdout, "%s ", CM->ucorr[i_corr]->name);
  fflush(stdout);
  fputc('\n', stdout);
#endif

  if (CM->nmon<CM->ncor) {
    fprintf(stdout, "*** Warning: more correctors than monitors for %c plane.\n",  (coord==0?'x':'y'));
    fprintf(stdout, "*** Correction will probably be unstable!\n");
    fflush(stdout);
  }
  if (CM->ncor==0) {
    fprintf(stdout, "Warning: no correctors for %c plane.  No correction done.\n",  (coord==0?'x':'y'));
    fflush(stdout);
    return;
  }
  if (CM->nmon==0) {
    fprintf(stdout, "Warning: no monitors for %c plane.  No correction done.\n",  (coord==0?'x':'y'));
    fflush(stdout);
    CM->ncor = 0;
    return;
  }

  log_entry("compute_orbcor_matrices.1");
  /* allocate correction matrix for this plane, plus others: dK = T*Qo */
  m_alloc1(&CM->T , CM->ncor, CM->nmon);
  m_alloc1(&CM->Qo, CM->nmon, 1);
  m_alloc1(&CM->dK, CM->ncor, 1);
  m_alloc1(&CM->C , CM->nmon, CM->ncor);

  if (find_only) {
    log_exit("compute_orbcor_matrices");
    return;
  }

  /* intermediate matrices for computations: T = -I4*I2 */
  m_alloc1(&W , CM->nmon, CM->nmon);
  m_alloc1(&I1, CM->ncor, CM->nmon);        /* I1 = TRANS(C) */
  m_alloc1(&I2, CM->ncor, CM->nmon);        /* I2 = TRANS(C).W */
  m_alloc1(&I3, CM->ncor, CM->ncor);        /* I3 = TRANS(C).W.C */
  m_alloc1(&I4, CM->ncor, CM->ncor);        /* I4 = INVERSE(I3) */
  log_exit("compute_orbcor_matrices.1");

  log_entry("compute_orbcor_matrices.2");
  /* set up weight matrix */
  equalW = 1;
  for (i_moni=0; i_moni<CM->nmon; i_moni++) {
    long j;
    for (j=0; j<CM->nmon; j++)
      W->a[i_moni][j] = 0;
    W->a[i_moni][i_moni] = getMonitorWeight(CM->umoni[i_moni]);
    if (!i_moni)
      W0 = W->a[i_moni][i_moni];
    else if (W->a[i_moni][i_moni]!=W0)
      equalW = 0;
  }

  /* find transfer matrix from correctors to monitors */
  if ((coef = 2*sin(htune=PIx2*beamline->tune[coord?1:0]/2))==0)
    bomb("can't compute response matrix--beamline unstable", NULL);
  coefFL = beamline->matrix->R[4][5];
  
  corrFactor   = tmalloc(sizeof(*corrFactor)*CM->ncor);
  corrFactorFL = tmalloc(sizeof(*corrFactorFL)*CM->ncor);
  fprintf(stdout, "computing orbit response matrix...");
  fflush(stdout);
  switch (coord) {
  case 0:
    for (i_corr=0; i_corr<CM->ncor; i_corr++) {
      double betax, etax;
      betax = CM->ucorr[i_corr]->twiss->betax;
      etax = CM->ucorr[i_corr]->twiss->etax;
      if (CM->ucorr[i_corr]->pred) {
        betax = (betax + CM->ucorr[i_corr]->pred->twiss->betax)/2;
        etax = (etax + CM->ucorr[i_corr]->pred->twiss->etax)/2;
      }
      corrFactor[i_corr] = 
        getCorrectorCalibration(CM->ucorr[i_corr], coord)*sqrt(betax)/coef;
      corrFactorFL[i_corr] = 
        getCorrectorCalibration(CM->ucorr[i_corr], coord)*etax/coefFL;
    }
    for (i_moni=0; i_moni<CM->nmon; i_moni++) {
      moniFactor = 
        getMonitorCalibration(CM->umoni[i_moni], coord)*sqrt(CM->umoni[i_moni]->twiss->betax);
      for (i_corr=0; i_corr<CM->ncor; i_corr++) {
        double phi;
        phi = CM->ucorr[i_corr]->twiss->phix;
        if (CM->ucorr[i_corr]->pred)
          phi = (phi + CM->ucorr[i_corr]->pred->twiss->phix)/2;
        CM->C->a[i_moni][i_corr] 
          = moniFactor*corrFactor[i_corr]*
            cos(htune-fabs(CM->umoni[i_moni]->twiss->phix - phi));
        if (fixed_length)
          CM->C->a[i_moni][i_corr] += CM->umoni[i_moni]->twiss->etax*corrFactorFL[i_corr];
      }
    }
    break;
  default:
    for (i_corr=0; i_corr<CM->ncor; i_corr++) {
      double betay;
      betay = CM->ucorr[i_corr]->twiss->betay;
      if (CM->ucorr[i_corr]->pred)
        betay = (betay + CM->ucorr[i_corr]->pred->twiss->betay)/2;
      corrFactor[i_corr] = 
        getCorrectorCalibration(CM->ucorr[i_corr], coord)*sqrt(betay)/coef;
    }
    for (i_moni=0; i_moni<CM->nmon; i_moni++) {
      moniFactor = 
        getMonitorCalibration(CM->umoni[i_moni], coord)*sqrt(CM->umoni[i_moni]->twiss->betay);
      for (i_corr=0; i_corr<CM->ncor; i_corr++) {
        double phi;
        phi = CM->ucorr[i_corr]->twiss->phiy;
        if (CM->ucorr[i_corr]->pred)
          phi = (phi + CM->ucorr[i_corr]->pred->twiss->phiy)/2;
        CM->C->a[i_moni][i_corr] 
          = moniFactor*corrFactor[i_corr]*
            cos(htune-fabs(CM->umoni[i_moni]->twiss->phiy - phi));
      }
    }
  }
  free(corrFactor);
  report_stats(stdout, "\ndone");
#ifdef DEBUG
  m_show(CM->C    , "%13.6le ", "influence matrix\n", stdout);
#endif
  log_exit("compute_orbcor_matrices.2");

  if (invert && !CM->inverse_computed) {
    /* compute correction matrix T */
    CM->inverse_computed = 1;
    fprintf(stdout, "computing correction matrix...");
    fflush(stdout);
    if (equalW) {
      m_trans(I1, CM->C);
      m_mult(I3, I1, CM->C);
      m_invert(I4, I3);
      m_mult(CM->T, I4, I1);
      m_scmul(CM->T, CM->T, -W0); 
    } else {
      m_trans(I1, CM->C);
      m_mult(I2, I1, W);
      m_mult(I3, I2, CM->C);
      m_invert(I4, I3);
      m_mult(CM->T, I4, I2);
      m_scmul(CM->T, CM->T, -1.0); 
    }
    report_stats(stdout, "\ndone");
#ifdef DEBUG
    m_show(CM->T, "%13.6le ", "correction matrix\n", stdout);
#endif
  }

  log_exit("compute_orbcor_matrices");
}

long orbcor_plane(CORMON_DATA *CM, STEERING_LIST *SL, long coord, TRAJECTORY **orbit, long n_iterations, 
                  double clorb_acc, long clorb_iter, double clorb_iter_frac, RUN *run, LINE_LIST *beamline, double *closed_orbit, double *Cdp)
{
  ELEMENT_LIST *corr, *eptr;
  TRAJECTORY *clorb;
  VMATRIX *M;
  long iteration, kick_offset;
  long i_moni, i_corr, i, sl_index;
  double dp, x, y, reading;
  double last_rms_pos, best_rms_pos, rms_pos, corr_fraction;
  double fraction, minFraction, param, change;
  
  log_entry("orbcor_plane");

  if (!CM)
    bomb("NULL CORMON_DATA pointer passed to orbcor_plane", NULL);
  if (!orbit)
    bomb("NULL TRAJECTORY pointer passed to orbcor_plane", NULL);
  if (!run)
    bomb("NULL RUN pointer passed to orbcor_plane", NULL);
  if (!beamline)
    bomb("NULL LINE_LIST pointer passed to orbcor_plane", NULL);
  if (!CM->inverse_computed)
    bomb("No inverse matrix computed prior to orbcor_plane", NULL);

  if (closed_orbit)
    dp = closed_orbit[5];
  else
    dp = 0;

  clorb = orbit[0];
  if (closed_orbit) {
    for (i=0; i<4; i++)
      clorb[0].centroid[i] = closed_orbit[i];
  }
  else {
    for (i=0; i<4; i++)
      clorb[0].centroid[i] = 0;
  }

  if (CM->nmon<CM->ncor) {
    fprintf(stdout, "*** Warning: more correctors than monitors for %c plane.\n",  (coord==0?'x':'y'));
    fprintf(stdout, "*** Correction will probably be unstable!\n");
    fflush(stdout);
  }

  best_rms_pos = rms_pos = DBL_MAX/4;
  corr_fraction = CM->corr_fraction;
  for (iteration=0; iteration<=n_iterations; iteration++) {
    clorb = orbit[iteration?1:0];
    delete_phase_references();
    if (!(M = beamline->matrix)) {
      if (beamline->elem_twiss)
        M = beamline->matrix = full_matrix(beamline->elem_twiss, run, 1);
      else
        M = beamline->matrix = full_matrix(&(beamline->elem), run, 1);
    }
    if (!M || !M->R)
      bomb("problem calculating full matrix (orbcor_plane)", NULL);
    if (iteration==1)
      for (i=0; i<6; i++)
        orbit[1][0].centroid[i] = orbit[0][0].centroid[i];
    if (!find_closed_orbit(clorb, clorb_acc, clorb_iter, beamline, M, run, dp, 1, CM->fixed_length, NULL, 
                           clorb_iter_frac, NULL)) {
      fprintf(stdout, "Failed to find closed orbit.\n");
      fflush(stdout);
      return(-1);
    }
    if (Cdp)
      Cdp[iteration] = clorb[0].centroid[5];
    if (closed_orbit) {
      for (i=0; i<6; i++) 
        closed_orbit[i] = clorb[0].centroid[i];
    }

    /* find readings at monitors and add in reading errors */
    i = 1;
    last_rms_pos = rms_pos;
    if (best_rms_pos>rms_pos)
      best_rms_pos = rms_pos;
    rms_pos = 0;
    for (i_moni=0; i_moni<CM->nmon; i_moni++, i++) {
      while (CM->umoni[i_moni]!=clorb[i].elem) {
        if (clorb[i].elem->succ)
          i++;
        else
          bomb("missing monitor in closed orbit", NULL);
      }

      x = clorb[i].centroid[0];
      y = clorb[i].centroid[2];
      eptr = clorb[i].elem;
      reading = computeMonitorReading(eptr, coord, x, y, 0);

      CM->Qo->a[i_moni][0] = CM->posi[iteration][i_moni] = reading + 
        (CM->bpm_noise?noise_value(CM->bpm_noise, CM->bpm_noise_cutoff, CM->bpm_noise_distribution):0.0);
      rms_pos += sqr(CM->Qo->a[i_moni][0]);
      if (!clorb[i].elem->succ)
        break;
    }
    rms_pos = sqrt(rms_pos/CM->nmon);
    if (iteration==0 && rms_pos>1e9) {
      /* if the closed orbit has RMS > 1e9m, I assume correction won't work and routine bombs */
      fprintf(stdout, "Orbit beyond 10^9 m.  Aborting correction.\n");
      fflush(stdout);
      return(-1);
    }
    if (rms_pos>best_rms_pos*1.01) {
      if (corr_fraction==1 || corr_fraction==0)
        break;            
      fprintf(stdout, "orbit diverging on iteration %ld: RMS=%e m (was %e m)--redoing with correction fraction of %e\n", 
              iteration, rms_pos, best_rms_pos, corr_fraction = sqr(corr_fraction));
      best_rms_pos = rms_pos;
      fflush(stdout);
#if 0
      if (iteration>=1) {
        /* step through beamline and change correctors back to last values */
        for (i_corr=0; i_corr<CM->ncor; i_corr++) {
          corr = CM->ucorr[i_corr];
          sl_index = CM->sl_index[i_corr];
          kick_offset = SL->param_offset[sl_index];
          *((double*)(corr->p_elem+kick_offset)) = CM->kick[iteration-1][i_corr]/CM->kick_coef[i_corr];
          if (corr->matrix) {
            free_matrices(corr->matrix);
            tfree(corr->matrix);
            corr->matrix = NULL;
          }
          compute_matrix(corr, run, NULL);
        }
      }
      if (beamline->links)
        assert_element_links(beamline->links, run, beamline, DYNAMIC_LINK);

      /* free concatenated matrix to force it to be recomputed with new corrector matrices */
      free_matrices(beamline->matrix);
      tfree(beamline->matrix);
      beamline->matrix = NULL;
      iteration -= 1;
      rms_pos = last_rms_pos;
      continue;
#endif
    }

    if (CM->fixed_length)
      dp = clorb[0].centroid[5];

    if (iteration==n_iterations)
      break;

    /* solve for the corrector kicks */
    m_mult(CM->dK, CM->T, CM->Qo);

#ifdef DEBUG
    m_show(CM->Qo, "%13.6le ", "traj matrix\n", stdout);
    m_show(CM->dK, "%13.6le ", "kick matrix\n", stdout);
#endif

    /* see if any correctors are over their limit */
    minFraction = 1;
    for (i_corr=0; i_corr<CM->ncor; i_corr++) {
      corr = CM->ucorr[i_corr];
      sl_index = CM->sl_index[i_corr];
      kick_offset = SL->param_offset[sl_index];
      param = fabs(*((double*)(corr->p_elem+kick_offset)) +
                   (change=CM->dK->a[i_corr][0]/CM->kick_coef[i_corr]*corr_fraction));
      if (SL->corr_limit[sl_index] && param>SL->corr_limit[sl_index]) {
        fraction = fabs((SL->corr_limit[sl_index]-fabs(*((double*)(corr->p_elem+kick_offset))))/change);
        if (fraction<minFraction)
          minFraction = fraction;
      }
    }
    fraction = minFraction*corr_fraction;
    
    /* step through beamline and change correctors */
    for (i_corr=0; i_corr<CM->ncor; i_corr++) {
      corr = CM->ucorr[i_corr];
      sl_index = CM->sl_index[i_corr];
      kick_offset = SL->param_offset[sl_index];
      if (iteration==0) 
        CM->kick[iteration][i_corr] = *((double*)(corr->p_elem+kick_offset))*CM->kick_coef[i_corr];
      *((double*)(corr->p_elem+kick_offset)) += CM->dK->a[i_corr][0]/CM->kick_coef[i_corr]*fraction;
      CM->kick[iteration+1][i_corr] = *((double*)(corr->p_elem+kick_offset))*CM->kick_coef[i_corr];
      if (corr->matrix) {
        free_matrices(corr->matrix);
        tfree(corr->matrix);
        corr->matrix = NULL;
      }
      compute_matrix(corr, run, NULL);
    }
    if (beamline->links)
      assert_element_links(beamline->links, run, beamline, DYNAMIC_LINK);

    /* free concatenated matrix to force it to be recomputed with new corrector matrices */
    free_matrices(beamline->matrix);
    tfree(beamline->matrix);
    beamline->matrix = NULL;
  }

  if (rms_pos>1e9) {
    /* if the final closed orbit has RMS > 1e9m, I assume correction didn't work and routine bombs */
    fprintf(stdout, "Orbit beyond 1e9 m.  Aborting correction.\n");
    fflush(stdout);
    return(-1);
  }

  if (rms_pos>best_rms_pos+CM->corr_accuracy && (corr_fraction==1 || corr_fraction==0)) {
    fprintf(stdout, "orbit not improving--iteration terminated at iteration %ld with last result of %e m RMS\n",
            iteration, rms_pos);
    fflush(stdout);
    /* step through beamline and change correctors back to last values */
    for (i_corr=0; i_corr<CM->ncor; i_corr++) {
      corr = CM->ucorr[i_corr];
      sl_index = CM->sl_index[i_corr];
      kick_offset = SL->param_offset[sl_index];
      *((double*)(corr->p_elem+kick_offset)) = CM->kick[iteration][i_corr]/CM->kick_coef[i_corr];
      if (corr->matrix) {
        free_matrices(corr->matrix);
        tfree(corr->matrix);
        corr->matrix = NULL;
      }
      compute_matrix(corr, run, NULL);
    }
    if (beamline->links)
      assert_element_links(beamline->links, run, beamline, DYNAMIC_LINK);

    /* free concatenated matrix to force it to be recomputed with new corrector matrices */
    free_matrices(beamline->matrix);
    tfree(beamline->matrix);
    beamline->matrix = NULL;
  }

  /* indicate that beamline concatenation and Twiss computation (if wanted) are not current */
  beamline->flags &= ~BEAMLINE_CONCAT_CURRENT;
  beamline->flags &= ~BEAMLINE_TWISS_CURRENT;

  log_exit("orbcor_plane");
  return(iteration);
}

long find_closed_orbit(TRAJECTORY *clorb, double clorb_acc, long clorb_iter, LINE_LIST *beamline, VMATRIX *M, RUN *run, 
                       double dp, long start_from_recirc, long fixed_length, double *starting_point, double change_fraction,
                       double *deviation)
{
  static MATRIX *R, *ImR, *INV_ImR, *INV_R, *C, *co, *diff, *change;
  static double **one_part;
  static long initialized = 0, been_warned = 0;
  long i, j, n_iter = 0, bad_orbit;
  long n_part;
  double p, error, last_error;

  log_entry("find_closed_orbit");

  if (fixed_length)
    return findFixedLengthClosedOrbit(clorb, clorb_acc, clorb_iter, beamline, M, run, dp,
                                      start_from_recirc, starting_point, change_fraction, deviation);
  
  /* method for finding closed orbit: 
   * 1. solve co[i] = C[i] + R[i][j]*co[j] for co[i]:
   *        co = INV(I-R)*C
   * 2. use Newton's method iteration starting with this solution:
   *        dco = INV(R)*(-co + F(co))
   *    where F(co) returns the coordinates at the end for starting
   *    coordinates co.
   */

  if (!M)
    bomb("no transport matrix passed to find_closed_orbit()", NULL);
  if (!M->R)
    bomb("faulty transport matrix passed to find_closed_orbit()", NULL);

  if (!initialized) {
    m_alloc(&ImR, 4, 4);
    m_alloc(&R, 4, 4);
    m_alloc(&INV_ImR, 4, 4);
    m_alloc(&INV_R, 4, 4);
    m_alloc(&C, 4, 1);
    m_alloc(&co, 4, 1);
    m_alloc(&diff, 4, 1);
    m_alloc(&change, 4, 1);
    one_part = (double**)zarray_2d(sizeof(**one_part), 1, 7);
    initialized = 1;
  }

  for (i=0; i<4; i++) {
    C->a[i][0] = M->C[i];
    for (j=0; j<4; j++) {
      R->a[i][j]   = M->R[i][j];
      ImR->a[i][j] = (i==j?1:0)-R->a[i][j];
    }
  }

  if (!m_invert(INV_ImR, ImR)) {
    fprintf(stdout, "error: unable to invert matrix to find closed orbit!\nThe R matrix is:");
    fflush(stdout);
    for (i=0; i<4; i++) {
      fprintf(stdout, "R[%ld]: ", i+1);
      fflush(stdout);
      for (j=0; j<4; j++)
        fprintf(stdout, "%14.6e ", R->a[i][j]);
      fflush(stdout);
      fputc('\n', stdout);
    }
    for (i=0; i<4; i++) {
      for (j=0; j<4; j++) {
        INV_ImR->a[i][j] = 0;
      }
    }
  }

  if (!starting_point) {
    if (!m_mult(co, INV_ImR, C))
      bomb("unable to solve for closed orbit--matrix multiplication error", NULL);
    
    one_part[0][0] = co->a[0][0];
    one_part[0][1] = co->a[1][0];
    one_part[0][2] = co->a[2][0];
    one_part[0][3] = co->a[3][0];
    one_part[0][4] = 0;
    one_part[0][5] = dp;
  }
  else {
    one_part[0][0] = co->a[0][0] = starting_point[0];
    one_part[0][1] = co->a[1][0] = starting_point[1];
    one_part[0][2] = co->a[2][0] = starting_point[2];
    one_part[0][3] = co->a[3][0] = starting_point[3];
    one_part[0][4] = 0;
    one_part[0][5] = dp;
  }

  p = run->p_central;
  error = DBL_MAX/4;
  bad_orbit = 0;
  if (deviation)
    deviation[4] = deviation[5] = 0;
  do {
    n_part = 1;
    do_tracking(one_part, &n_part, NULL, beamline, &p, (double**)NULL, (BEAM_SUMS**)NULL, (long*)NULL,
                clorb+1, run, 0, 
                TEST_PARTICLES+TIME_DEPENDENCE_OFF+(start_from_recirc?BEGIN_AT_RECIRC:0), 1, 0,
                NULL, NULL, NULL);
    if (n_part==0)
      bomb("tracking failed for closed-orbit test particle", NULL);
#ifdef DEBUG
    fprintf(stdout, "particle coordinates at end of closed-orbit tracking:\n%e %e %e %e %e %e\n",
            one_part[0][0], one_part[0][1], one_part[0][2], one_part[0][3],
            one_part[0][4], one_part[0][5]);
    fflush(stdout);
#endif
    for (i=0; i<4; i++) {
      diff->a[i][0] = one_part[0][i] - co->a[i][0];
      if (deviation)
        deviation[i] = diff->a[i][0];
    }
    last_error = error;
    if ((error = sqrt(sqr(diff->a[0][0]) + sqr(diff->a[1][0]) + sqr(diff->a[2][0]) + sqr(diff->a[3][0])))<clorb_acc)
      break;
    if (error>2*last_error) {
      fprintf(stdout, "warning: closed orbit diverging--iteration stopped\n");
      fflush(stdout);
      fprintf(stdout, "last error was %e, current is %e\n", last_error, error);
      fflush(stdout);
      n_iter = clorb_iter;
      break;
    }
    if (change_fraction) {
      m_mult(change, INV_ImR, diff);
      if (change_fraction!=1)
        m_scmul(change, change, change_fraction);
      m_add(co, co, change);
      for (i=0; i<4; i++)
        one_part[0][i] = co->a[i][0];
    }
    else {
      for (i=0; i<4; i++) {
        co->a[i][0] = (co->a[i][0]+one_part[0][i])/2;
        one_part[0][i] = co->a[i][0];
      }
    }
    one_part[0][4] = 0;
    one_part[0][5] = dp;
  } while (++n_iter<clorb_iter);
  if (n_iter==clorb_iter)  {
    fprintf(stdout, "warning: closed orbit did not converge to better than %e after %ld iterations\n",
            error, n_iter);
    fflush(stdout);
    if (isnan(error) || isinf(error)) {
      return 0;
    }
    fflush(stdout);
  }
#ifdef DEBUG
  fprintf(stdout, "final closed-orbit after %ld iterations:\n%e %e %e %e %e %e\n",
          n_iter, one_part[0][0], one_part[0][1], one_part[0][2], one_part[0][3],
          one_part[0][4], one_part[0][5]);
  fflush(stdout);
#endif
  clorb[0].centroid[0] = one_part[0][0];
  clorb[0].centroid[1] = one_part[0][1];
  clorb[0].centroid[2] = one_part[0][2];
  clorb[0].centroid[3] = one_part[0][3];
  clorb[0].centroid[4] = 0;
  clorb[0].centroid[5] = dp;

  log_exit("find_closed_orbit");
  if (bad_orbit)
    return(0);
  return(1);
}

long findFixedLengthClosedOrbit(TRAJECTORY *clorb, double clorb_acc, long clorb_iter, LINE_LIST *beamline, VMATRIX *M, RUN *run, 
                                double dp, long start_from_recirc, double *starting_point, double change_fraction,
                                double *deviation)
{
  long nElems, iterationsLeft, i, iterationsDone;
  double error, ds, last_dp, last_ds;
  
  nElems = beamline->n_elems;
  iterationsLeft = clorb_iter/10+10;
  last_ds = last_dp = sqrt(DBL_MAX/10);
  iterationsDone = 0;
  while (iterationsDone<iterationsLeft) {
    if (!find_closed_orbit(clorb, clorb_acc, clorb_iter, beamline, M, run, dp, start_from_recirc,
                           0, starting_point, change_fraction, deviation))
      return 0;
    ds = clorb[nElems].centroid[4] - beamline->revolution_length;
    if (deviation)
      deviation[4] = ds;
    for (i=error=0; i<4; i++)
      error += sqr(clorb[nElems].centroid[i]-clorb[0].centroid[i]);
    /* The error for delta is scaled by 1/change_fraction so that it
     * gets a chance to converge even with a small change_fraction.
     * Otherwise, the delta iteration may not converge even though the
     * individual orbits are converging.
     */
    error = sqrt(error + sqr((last_ds-ds)/(beamline->revolution_length*change_fraction)));
    /* 
      fprintf(stdout, "orbit error for dp=%le is %le, ds=%le:\n", dp, error, ds);
      for (i=0; i<6; i++)
      fprintf(stdout, "%10.3e ", clorb[0].centroid[i]);
      fprintf(stdout, "\n");
      for (i=0; i<6; i++)
      fprintf(stdout, "%10.3e ", clorb[nElems].centroid[i]-(i==4?beamline->revolution_length:0));
      fprintf(stdout, "\n");
      */
    if (error<clorb_acc)
      break;
    last_ds = ds;
    last_dp = dp;
    dp -= change_fraction*ds/M->R[4][5];
    iterationsDone++;
  }
#if DEBUG
  fprintf(stdout, "%ld iterations done for delta in fixed-length orbit computation\ndelta convergence error was %le\ndelta=%le, length error was %le\n", 
          iterationsDone, last_dp-dp, dp, ds);
#endif
  if (iterationsDone<iterationsLeft)
    return 1;
  fprintf(stdout, "Warning: fixed length orbit iteration didn't converge (error is %le)\n", error);
  fprintf(stdout, "dp = %le, %le\n", dp, last_dp);
  for (i=0; i<6; i++)
    fprintf(stdout, "%10.3e ", clorb[0].centroid[i]);
  fprintf(stdout, "\n");
  for (i=0; i<6; i++)
    fprintf(stdout, "%10.3e ", clorb[nElems].centroid[i]-(i==4?beamline->revolution_length:0));
  fprintf(stdout, "\n");
  
  return 0;
}


void zero_closed_orbit(TRAJECTORY *clorb, long n)
{
  long i, j;

  for (i=0; i<n; i++)
    for (j=0; j<6; j++)
      clorb[i].centroid[j] = 0;
}


ELEMENT_LIST *next_element_of_type(ELEMENT_LIST *elem, long type)
{
  while (elem) {
    if (elem->type==type)
      return(elem);
    elem = elem->succ;
  }
  return(elem);
}

ELEMENT_LIST *next_element_of_type2(ELEMENT_LIST *elem, long type1, long type2)
{
  while (elem) {
    if (elem->type==type1 || elem->type==type2)
      return(elem);
    elem = elem->succ;
  }
  return(elem);
}

ELEMENT_LIST *next_element_of_types(ELEMENT_LIST *elem, long *type, long n_types, long *index)
{
  register long i;
  while (elem) {
    for (i=0; i<n_types; i++) {
      if (elem->type==type[i]) {
        *index = i;
        return(elem);
      }
    }
    elem = elem->succ;
  }
  return(elem);
}

long find_parameter_offset(char *param_name, long elem_type)
{
  long param;
  if ((param=confirm_parameter(param_name, elem_type))<0)
    return(-1);
  return(entity_description[elem_type].parameter[param].offset);
}


long zero_correctors(ELEMENT_LIST *elem, RUN *run, CORRECTION *correct)
{
  return 
    zero_hcorrectors(elem, run, correct) +
      zero_vcorrectors(elem, run, correct);
}

long zero_correctors_one_plane(ELEMENT_LIST *elem, RUN *run, STEERING_LIST *SL)
{  
  long n_zeroed = 0, i;
  long paramOffset;

  while (elem) {
    paramOffset = -1;
    for (i=0; i<SL->n_corr_types; i++)
      if (strcmp(elem->name, SL->corr_name[i])==0)
        break;
    if (i!=SL->n_corr_types) {
      paramOffset = SL->param_offset[i];
      *((double*)(elem->p_elem+paramOffset)) = 0;
      if (elem->matrix) {
        free_matrices(elem->matrix);
        tfree(elem->matrix);
        elem->matrix = NULL;
      }
      compute_matrix(elem, run, NULL);
      n_zeroed++;
    }
    elem = elem->succ;
  }
  return(n_zeroed);
}

long zero_hcorrectors(ELEMENT_LIST *elem, RUN *run, CORRECTION *correct)
{
  long nz;
  nz = zero_correctors_one_plane(elem, run, &(correct->SLx));
  fprintf(stderr, "%ld H correctors set to zero\n", nz);
  return nz;
}

long zero_vcorrectors(ELEMENT_LIST *elem, RUN *run, CORRECTION *correct)
{
  long nz;
  nz = zero_correctors_one_plane(elem, run, &(correct->SLy));
  fprintf(stderr, "%ld V correctors set to zero\n", nz);
  return nz;
}

void rotate_xy(double *x, double *y, double angle)
{
    static double x0, y0, ca, sa;

    x0 = *x;
    y0 = *y;
    *x  =  x0*(ca=cos(angle)) + y0*(sa=sin(angle));
    *y  = -x0*sa              + y0*ca;
    }

double rms_value(double *data, long n_data)
{
    double sum2;
    long i;

    for (i=sum2=0; i<n_data; i++) 
        sum2 += sqr(data[i]);
    if (n_data)
        return(sqrt(sum2/n_data));
    else
        return(0.0);
    }

long steering_corrector(ELEMENT_LIST *eptr, STEERING_LIST *SL)
{
  long i;
  for (i=0; i<SL->n_corr_types; i++)
    if (strcmp(eptr->name, SL->corr_name[i])==0) {
      switch (eptr->type) {
      case T_HVCOR:
        return ((HVCOR*)(eptr->p_elem))->steering;
        break;
      case T_HCOR:
        return ((HCOR*)(eptr->p_elem))->steering;
        break;
      case T_VCOR:
        return ((VCOR*)(eptr->p_elem))->steering;
        break;
      default:
        return 1;
        break;
      }
    }
  return 0;
}

long find_index(long key, long *list, long n_listed)
{
    register long i;
    for (i=0; i<n_listed; i++)
        if (key==list[i])
            return(i);
    return(-1);
    }

double noise_value(double xamplitude, double xcutoff, long xerror_type)
{
    switch (xerror_type) {
        case UNIFORM_ERRORS:
            return(2*xamplitude*(random_3(0)-0.5));
        case GAUSSIAN_ERRORS:
            return(gauss_rn_lim(0.0, xamplitude, xcutoff, random_3));
        case PLUS_OR_MINUS_ERRORS:
            /* return a number on [-x-1, x-1]), which is added to 1 in the calling routine
             * (since these are implemented as fractional errors)
             */
            return(xamplitude*(random_3(0)>0.5?1.0:-1.0) - 1);
        default:
            bomb("unknown error type in perturbation()", NULL);
            exit(1);
            break;
        }
    return(0.0);
    }

double computeMonitorReading(ELEMENT_LIST *elem, long coord, double x, double y,
                             unsigned long flags)
/* coord = 0 is x, otherwise y */
{
  double calibration, tilt, reading;
  char *equation;

  switch (elem->type) {
  case T_MONI:  
    x -= ((MONI*)(elem->p_elem))->dx;
    y -= ((MONI*)(elem->p_elem))->dy;
    if (coord==0)
      calibration = ((MONI*)(elem->p_elem))->xcalibration;
    else 
      calibration = ((MONI*)(elem->p_elem))->ycalibration;
    tilt = ((MONI*)(elem->p_elem))->tilt;
    if (flags&COMPUTEMONITORREADING_TILT_0)
      tilt = 0;
    if (tilt)
      rotate_xy(&x, &y, tilt);   
    if (coord==0)
      equation = ((MONI*)(elem->p_elem))->x_readout; 
    else
      equation = ((MONI*)(elem->p_elem))->y_readout;
    break;
  case T_HMON: 
    x -= ((HMON*)(elem->p_elem))->dx;
    y -= ((HMON*)(elem->p_elem))->dy;
    calibration = ((HMON*)(elem->p_elem))->calibration;
    tilt = ((HMON*)(elem->p_elem))->tilt;
    if (flags&COMPUTEMONITORREADING_TILT_0)
      tilt = 0;
    if (tilt)
      rotate_xy(&x, &y, tilt);   
    equation = ((HMON*)(elem->p_elem))->readout; 
    if (coord!=0)
      bomb("element in horizontal monitor list is not a vertical monitor--internal logic error", NULL);
    break;
  case T_VMON:
    x -= ((VMON*)(elem->p_elem))->dx;
    y -= ((VMON*)(elem->p_elem))->dy;
    calibration = ((VMON*)(elem->p_elem))->calibration;
    tilt = ((VMON*)(elem->p_elem))->tilt;
    if (flags&COMPUTEMONITORREADING_TILT_0)
      tilt = 0;
    if (tilt)
      rotate_xy(&x, &y, tilt);   
    equation = ((VMON*)(elem->p_elem))->readout; 
    if (!coord)
      bomb("element in vertical monitor list is not a vertical monitor--internal logic error", NULL);
    break;
  default:
    fprintf(stdout, "error: element %s found in monitor list--internal logic error\n", 
            elem->name);
    fflush(stdout);
    abort();
    break;
  }
  
  if (flags&COMPUTEMONITORREADING_CAL_1)
    calibration = 1;

  if (equation) {
    rpn_clear();
    rpn_store(x, rpn_x_mem);
    rpn_store(y, rpn_y_mem);
    reading = rpn(equation)*calibration;
    if (rpn_check_error()) exit(1);
  }
  else {
    switch (coord) {
    case 0: 
      reading = x*calibration;
      break;
    default:
      reading = y*calibration;
      break;
    }
  }
  return reading;
}

double getMonitorWeight(ELEMENT_LIST *elem)
{
    switch (elem->type) {
      case T_MONI:
        return ((MONI*)(elem->p_elem))->weight;
      case T_HMON:
        return ((HMON*)(elem->p_elem))->weight;
      case T_VMON:
        return ((VMON*)(elem->p_elem))->weight;
        }
    bomb("invalid element type in getMonitorWeight()", NULL);
	return(1);
    }

double getMonitorCalibration(ELEMENT_LIST *elem, long coord)
{
    switch (elem->type) {
      case T_HMON:
        return ((HMON*)(elem->p_elem))->calibration;
      case T_VMON:
        return ((VMON*)(elem->p_elem))->calibration;
      case T_MONI:
        if (coord)
            return ((MONI*)(elem->p_elem))->ycalibration;
        return ((MONI*)(elem->p_elem))->xcalibration;
        }
    bomb("invalid element type in getMonitorCalibration()", NULL);
	return(1);
    }

void setMonitorCalibration(ELEMENT_LIST *elem, double calib, long coord)
{
  switch (elem->type) {
  case T_HMON:
    ((HMON*)(elem->p_elem))->calibration = calib;
    break;
  case T_VMON:
    ((VMON*)(elem->p_elem))->calibration = calib;
    break;
  case T_MONI:
    if (coord)
      ((MONI*)(elem->p_elem))->ycalibration = calib;
    else 
      ((MONI*)(elem->p_elem))->xcalibration = calib;
    break;
  default:
    bomb("invalid element type in setMonitorCalibration()", NULL);
    break;
  }
}


double getCorrectorCalibration(ELEMENT_LIST *elem, long coord)
{
    switch (elem->type) {
      case T_HCOR:
        return ((HCOR*)(elem->p_elem))->calibration;
      case T_VCOR:
        return ((VCOR*)(elem->p_elem))->calibration;
      case T_HVCOR:
        if (coord)
            return ((HVCOR*)(elem->p_elem))->ycalibration;
        return ((HVCOR*)(elem->p_elem))->xcalibration;
      default:
        return 1;
        }
    }
