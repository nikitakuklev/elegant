/* Copyright 1994 by Michael Borland and Argonne National Laboratory,
 * all rights reserved.
 */
/* routine: do_tracking()
 * purpose: track a collection of particles through a beamline
 *
 * Michael Borland, 1989
 */
#include "mdb.h"
#include "mdbsun.h"
#include "track.h"
/* #include "smath.h" */

ELEMENT_LIST *findBeamlineMatrixElement(ELEMENT_LIST *eptr);
void trackLongitudinalOnlyRing(double **part, long np, VMATRIX *M, double *alpha);
void store_fitpoint_matrix_values(MARK *fpt, char *name, long occurence, VMATRIX *M);
void trackWithChromaticLinearMatrix(double **particle, long particles,
                                    ELEMENT_LIST *eptr,
                                    TWISS *twiss0,
                                    double *tune0,
                                    double *chrom,
                                    double *dbeta_dPoP, 
                                    double *dalpha_dPoP);

static TRACKING_CONTEXT trackingContext;

double beta_from_delta(double p, double delta)
{
  p *= 1+delta;
  return( p/sqrt(p*p+1));
}

long do_tracking(
                 double **coord,
                 long *n_original,
                 long *effort,
                 LINE_LIST *beamline,
                 double *P_central,    /* beta*gamma for central particle */
                 double **accepted,
                 BEAM_SUMS **sums_vs_z,
                 long *n_z_points,
                 TRAJECTORY *traj_vs_z,
                 RUN *run,
                 long step,
                 unsigned long flags,
                 long n_passes,
                 SASEFEL_OUTPUT *sasefel,
		 SLICE_OUTPUT *sliceAnalysis,
                 double *finalCharge
                 )
{
  RFMODE *rfmode; TRFMODE *trfmode;
  WATCH *watch;
  HISTOGRAM *histogram;
  ENERGY *energy;
  MAXAMP *maxamp;
  MALIGN *malign;
  ELEMENT_LIST *eptr, *eptrCLMatrix=NULL;
  long n_left, show_dE;
  double dgamma, dP[3], z, z_recirc, last_z;
  long i, j, i_traj=0, i_sums, n_to_track, i_pass, isConcat;
  long i_sums_recirc, saveISR=0;
  long watch_pt_seen;
  double sum, x_max, y_max;
  long elliptical;
  double et1, et2;
  long is_batch = 0, last_type;
  static long is_ansi_term = -1;
  char s[100], *name;
  long check_nan, sums_allocated = 0;
  long elementsTracked, sliceAnDone = 0;
  CHARGE *charge;
  static long warnedAboutChargePosition = 0;
  

  if (!coord)
    bomb("Null particle coordinate array! (do_tracking)", NULL);

#ifdef WATCH_MEMORY
  fprintf(stdout, "start do_tracking():  CPU: %6.2lf  PF: %6ld  MEM: %6ld\n",
          cpu_time()/100.0, page_faults(), memory_count());
  fflush(stdout);
#endif
  
#if defined(UNIX) || defined(_WIN32)
  if (is_ansi_term==-1) {
    char *ptr;
    is_ansi_term = 1;
    if (!(ptr=getenv("TERM")))
      is_ansi_term = 0;
    else if (strcmp(ptr, "emacs")==0)
      is_ansi_term = 0;
  }
#endif
  
  if (accepted)
    copy_particles(accepted, coord, *n_original);

#ifdef VAX_VMS
  is_batch = job_mode(getpid())==2?1:0;
#endif
  
  z = z_recirc = last_z =  0;
  i_sums = i_sums_recirc = 0;
  x_max = y_max = 0;
  n_to_track = *n_original;
  et1 = -2.0;
  elliptical = isConcat = 0;
  n_left = n_to_track;
  watch_pt_seen = 0;
  
  check_nan = 1;
  eptr = &(beamline->elem);

  if (flags&FIRST_BEAM_IS_FIDUCIAL && !(flags&FIDUCIAL_BEAM_SEEN)) {
    /* this is required just in case rf elements etc. have previously
     * been activated by computation of correction matrices or trajectory
     * correction.
     */
    delete_phase_references();
    reset_special_elements(beamline, 1);
  }
  reset_driftCSR();
  
  while (eptr) {
    if (check_nan) 
      switch (eptr->type) {
      case T_RCOL: case T_ECOL: case T_MAXAMP:
        check_nan = 0;
        break;
      default:
        break;
      }
    if (flags&FIRST_BEAM_IS_FIDUCIAL && !(flags&FIDUCIAL_BEAM_SEEN))
      eptr->Pref_output_fiducial = 0;
    eptr = eptr->succ;
  }
  if (!(flags&FIDUCIAL_BEAM_SEEN) && flags&PRECORRECTION_BEAM)
    flags &= ~FIRST_BEAM_IS_FIDUCIAL; 
  if (flags&FIRST_BEAM_IS_FIDUCIAL && !(flags&FIDUCIAL_BEAM_SEEN) && !(flags&SILENT_RUNNING)) {
    fprintf(stdout, "This step establishes energy profile vs s (fiducial beam).\n");
    fflush(stdout);
  }
  
  log_exit("do_tracking.1");
  log_entry("do_tracking.2");
  name = "_BEG_";
  last_type = sums_allocated = 0;
  charge = NULL;
  if (finalCharge)
    *finalCharge = 0;
  for (i_pass=0; i_pass<n_passes; i_pass++) {
    log_entry("do_tracking.2.1");
    
    if (beamline->links) {
      sprintf(s, "%.15e sto p_central  %ld sto turn", *P_central, i_pass);
      rpn(s);
      rpn_clear();
      if (rpn_check_error()) exit(1);
      if (assert_element_links(beamline->links, run, beamline, TURN_BY_TURN_LINK)) {
        beamline->flags &= ~BEAMLINE_CONCAT_CURRENT;
        beamline->flags &= ~BEAMLINE_TWISS_CURRENT;
      }
    }
    if (beamline->flags&BEAMLINE_TWISS_WANTED && !(beamline->flags&BEAMLINE_TWISS_CURRENT)
        && !(flags&TEST_PARTICLES)) {
      update_twiss_parameters(run, beamline, NULL);
    }
    if (run->concat_order && !(flags&TEST_PARTICLES) && 
        !(beamline->flags&BEAMLINE_CONCAT_CURRENT) ) {
      /* form concatenated beamline for tracking */
      concatenate_beamline(beamline, run);
    }

    if (run->concat_order && beamline->flags&BEAMLINE_CONCAT_DONE &&
        !(flags&TEST_PARTICLES)) {
      if (beamline->ecat_recirc && (i_pass || flags&BEGIN_AT_RECIRC))
        eptr = beamline->ecat_recirc;
      else
        eptr = &(beamline->ecat);
      isConcat = 1;
    }
    else if (beamline->elem_recirc && (i_pass || flags&BEGIN_AT_RECIRC))
      eptr = beamline->elem_recirc;
    else
      eptr = &(beamline->elem);
    
    if (flags&LINEAR_CHROMATIC_MATRIX) {
      if (!isConcat) {
        fprintf(stdout, "Error: in order to use the \"linear chromatic matrix\" for\n");
        fflush(stdout);
        fprintf(stdout, "tracking, you must ask for matrix concatenation in the run_setup.\n");
        fflush(stdout);
        exit(1);
      }
      eptrCLMatrix = findBeamlineMatrixElement(eptr);
    }
    if (flags&LONGITUDINAL_RING_ONLY) {
      if (!isConcat) {
        fprintf(stdout, "Error: in order to use the \"longitudinal ring\" mode of\n");
        fflush(stdout);
        fprintf(stdout, "tracking, you must ask for matrix concatenation in the run_setup.\n");
        fflush(stdout);
        exit(1);
      }
      eptrCLMatrix = findBeamlineMatrixElement(eptr);
    }
    
    if (sums_vs_z && n_z_points) {
      if (!sums_allocated || !*sums_vs_z) {
        /* allocate storage for beam sums */
        if (!isConcat)
          *n_z_points = beamline->n_elems + 1 + 
            (run->wrap_around?0:(n_passes-1)*(beamline->n_elems-beamline->i_recirc));
        else
          *n_z_points = beamline->ncat_elems + 1 +
            (run->wrap_around?0:(n_passes-1)*(beamline->ncat_elems-beamline->i_recirc));
        if (flags&FINAL_SUMS_ONLY)
          *n_z_points = 0;
        *sums_vs_z = tmalloc(sizeof(**sums_vs_z)*(*n_z_points+1));
        zero_beam_sums(*sums_vs_z, *n_z_points+1);
        sums_allocated = 1;
      }
      else if (!run->combine_bunch_statistics && i_pass==0)
        zero_beam_sums(*sums_vs_z, *n_z_points+1);
    }

    if (run->wrap_around) {
      i_sums = i_sums_recirc;  /* ==0 for i_pass==0 */
      z = z_recirc;            /* ditto */
      last_z = z;
    }
    log_exit("do_tracking.2.1");
    log_entry("do_tracking.2.2");
    if (check_nan) {
      n_left = n_to_track = limit_amplitudes(coord, DBL_MAX, DBL_MAX, n_to_track, accepted, z, *P_central, 0);
    }
    if (!(flags&SILENT_RUNNING) && !is_batch && n_passes!=1 && !(flags&TEST_PARTICLES)
        && !(run->tracking_updates==0)) {
#if defined(VAX_VMS)
      sprintf(s, "%ld particles left after pass %ld        ",
              n_to_track, i_pass);
      fputs(s, stdout);
      if (is_ansi_term)
        backspace(strlen(s));
      else
        fputc('\n', stdout);
      fflush(stdout);
      et1 = et2;
#endif
#if defined(UNIX) || defined(_WIN32)
      if ((et2=delapsed_time())-et1>2.0) {
        sprintf(s, "%ld particles left after pass %ld        ", 
                n_to_track, i_pass);
        fputs(s, stdout);
        if (is_ansi_term)
          backspace(strlen(s));
        else
          fputc('\n', stdout);
        fflush(stdout);
        et1 = et2;
      }
#else
      sprintf(s, "%ld particles left after pass %ld        ", 
              n_to_track, i_pass);
      fputs(s, stdout);
      if (is_ansi_term)
        backspace(strlen(s));
      else
        fputc('\n', stdout);
      fflush(stdout);
#endif 
    }
    elementsTracked = -1;
    while (eptr && n_to_track) {
      elementsTracked++;
      log_entry("do_tracking.2.2.0");
      if (!eptr->name) {
        fprintf(stdout, "error: element ending at %em has NULL name pointer\n", eptr->end_pos);
        fflush(stdout);
        if (eptr->pred && eptr->pred->name) {
          fprintf(stdout, "previous element is %s\n", eptr->pred->name);
          fflush(stdout);
        }
        else if (eptr->succ && eptr->succ->name) {
          fprintf(stdout, "next element is %s\n", eptr->succ->name);
          fflush(stdout);
        }
        abort();
      }
      if (!eptr->p_elem && !run->concat_order) {
        fprintf(stdout, "element %s has NULL p_elem pointer", eptr->name);
        fflush(stdout);
        exit(1);
      }
      if (eptr->type<=0 || eptr->type>=N_TYPES) {
        fprintf(stdout, "element %s has type %ld--not recognized/not allowed\n", eptr->name, eptr->type);
        fflush(stdout);
        exit(1);
      }
      log_exit("do_tracking.2.2.0");

      log_entry("do_tracking.2.2.1");
      if (sums_vs_z && *sums_vs_z && !(flags&FINAL_SUMS_ONLY) && !(flags&TEST_PARTICLES)) {
        if (i_sums<0)
          bomb("attempt to accumulate beam sums with negative index!", NULL);
        accumulate_beam_sums(*sums_vs_z+i_sums, coord, n_to_track, *P_central);
        (*sums_vs_z)[i_sums].z = z;
#if defined(BEAM_SUMS_DEBUG)
        fprintf(stdout, "beam sums accumulated in slot %ld for %s at z=%em, sx=%e\n", 
                i_sums, name, z, sqrt((*sums_vs_z)[i_sums].sum2[0]/n_left));
        fflush(stdout);
#endif
        i_sums++;
      }
      name = eptr->name;
      last_z = z;
      if (entity_description[eptr->type].flags&HAS_LENGTH && eptr->p_elem)
        z += ((DRIFT*)eptr->p_elem)->length;
      else {
        if (eptr->pred)
          z += eptr->end_pos - eptr->pred->end_pos;
        else
          z += eptr->end_pos;
      }
      /* fill a structure that can be used to pass to other routines 
       * information on the tracking context 
       */
      strcpy(trackingContext.elementName, eptr->name);
      trackingContext.elementOccurrence = eptr->occurence;
      trackingContext.sliceAnalysis = sliceAnalysis;
      trackingContext.zStart = last_z;
      trackingContext.zEnd = z;
      trackingContext.step = step;

      log_exit("do_tracking.2.2.1");
      if (eptr->p_elem || eptr->matrix) {
#ifdef VAX_VMS
        if (run->print_statistics && !(flags&TEST_PARTICLES))
          fprintf(stdout, "tracking through %s%c", eptr->name, ' ');
          fflush(stdout);
#else
        if (run->print_statistics && !(flags&TEST_PARTICLES))
          fprintf(stdout, "tracking through %s%c", eptr->name, '\n');
          fflush(stdout);
#endif
        n_left = n_to_track;
        show_dE = 0;
        if (eptr==eptrCLMatrix)  {
          /* This element is the place-holder for the chromatic linear matrix or
           * the longitudinal-only matrix 
           */
          if (flags&LINEAR_CHROMATIC_MATRIX) 
            trackWithChromaticLinearMatrix(coord, n_to_track,
                                           eptrCLMatrix,
                                           beamline->twiss0,
                                           beamline->tune,
                                           beamline->chromaticity,
                                           beamline->dbeta_dPoP,
                                           beamline->dalpha_dPoP);
          else 
            trackLongitudinalOnlyRing(coord, n_to_track, 
                                      eptrCLMatrix->matrix,
                                      beamline->alpha);
        }
        else if (entity_description[eptr->type].flags&MATRIX_TRACKING) {
          if (!(entity_description[eptr->type].flags&HAS_MATRIX))
            bomb("attempt to matrix-multiply for element with no matrix!",  NULL);
          if (!eptr->matrix) {
            if (!(eptr->matrix=compute_matrix(eptr, run, NULL)))
              bomb("no matrix for element that must have matrix", NULL);
          }
          if (eptr->matrix->C[5]!=0) {
            fprintf(stdout, "Warning: matrix with C5!=0 detected in matrix multiplier--this shouldn't happen!\nAll particles considered lost!\n");
            fprintf(stdout, "Element in question is %s, C5=%le\n", eptr->name, eptr->matrix->C[5]);
            fflush(stdout);
            n_left = 0;
          } else {
            if (run->print_statistics>1 && !(flags&TEST_PARTICLES)) {
              fprintf(stdout, "Tracking matrix for %s\n", eptr->name);
              fflush(stdout);
              print_elem(stdout, eptr);
              print_matrices(stdout, "", eptr->matrix);
            }
            track_particles(coord, eptr->matrix, coord, n_to_track);
          }
        }
        else {
          if (run->print_statistics>1 && !(flags&TEST_PARTICLES)) {
            fprintf(stdout, "Tracking element: ");
            fflush(stdout);
            print_elem(stdout, eptr);
          }
          switch (eptr->type) {
          case T_CHARGE:
            if (i_pass==0) {
              if (elementsTracked!=0 && !warnedAboutChargePosition) {
                warnedAboutChargePosition = 1;
                fprintf(stdout, "Warning: the CHARGE element is not at the start of the beamline.\n");
                fflush(stdout);
              }
              if (charge!=NULL) {
                fprintf(stdout, "Fatal error: multipole CHARGE elements in one beamline.\n");
                fflush(stdout);
                exit(1);
              }
              charge = (CHARGE*)eptr->p_elem;
              charge->macroParticleCharge = 0;
              if (*n_original)
                charge->macroParticleCharge = charge->charge/(*n_original);
              if (charge->chargePerParticle)
                charge->macroParticleCharge = charge->chargePerParticle;
            }
            break;
          case T_MARK:
            if (((MARK*)eptr->p_elem)->fitpoint && i_pass==n_passes-1) {
              if (beamline->flags&BEAMLINE_TWISS_WANTED) {
                if (!(beamline->flags&BEAMLINE_TWISS_DONE))
                  update_twiss_parameters(run, beamline, NULL);
                store_fitpoint_twiss_parameters((MARK*)eptr->p_elem, eptr->name, 
                                                eptr->occurence, eptr->twiss);
              }
              store_fitpoint_matrix_values((MARK*)eptr->p_elem, eptr->name, 
                                           eptr->occurence, eptr->accumMatrix);
              store_fitpoint_beam_parameters((MARK*)eptr->p_elem, eptr->name,eptr->occurence, 
                                             coord, n_to_track, *P_central); 
            }
            break;
          case T_RECIRC:
            /* Recognize and record recirculation point.  */
            if (i_pass==0) {
              i_sums_recirc = i_sums-1;
              z_recirc = last_z;
            }
            break;
          case T_RFDF:
            if (!(flags&TIME_DEPENDENCE_OFF))
              track_through_rf_deflector(coord, (RFDF*)eptr->p_elem,
                                         coord, n_to_track, *P_central);
            else
              drift_beam(coord, n_to_track, ((RFDF*)eptr->p_elem)->length, run->default_order);
            break;
          case T_RMDF:
            if (!(flags&TIME_DEPENDENCE_OFF))
              track_through_ramped_deflector(coord, (RMDF*)eptr->p_elem,
                                             coord, n_to_track, *P_central);
            else
              drift_beam(coord, n_to_track, ((RMDF*)eptr->p_elem)->length, run->default_order);
            break;
          case T_RFTMEZ0:
          case T_TMCF:
          case T_CEPL:
          case T_TWPL:
            if (!(flags&TIME_DEPENDENCE_OFF)) {
              n_left = motion(coord, n_to_track, eptr->p_elem, eptr->type, P_central, 
                              &dgamma, dP, accepted, last_z);
              show_dE = 1;
            }
            else
              drift_beam(coord, n_to_track, ((TW_LINAC*)eptr->p_elem)->length, run->default_order);
            break;
          case T_MAPSOLENOID:
            n_left = motion(coord, n_to_track, eptr->p_elem, eptr->type, P_central, 
                            &dgamma, dP, accepted, last_z);
            break;
          case T_TWLA:
          case T_TWMTA:
            n_left = motion(coord, n_to_track, eptr->p_elem, eptr->type, P_central, 
                            &dgamma, dP, accepted, last_z);
            show_dE = 1;
            break;
          case T_RCOL:
            if (flags&TEST_PARTICLES && !(flags&TEST_PARTICLE_LOSSES))
              drift_beam(coord, n_to_track, ((RCOL*)eptr->p_elem)->length, run->default_order);
            else {
              n_left = rectangular_collimator(coord, (RCOL*)eptr->p_elem, n_to_track, accepted, last_z, *P_central);
            }
            break;
          case T_ECOL:
            if (flags&TEST_PARTICLES && !(flags&TEST_PARTICLE_LOSSES))
              drift_beam(coord, n_to_track, ((RCOL*)eptr->p_elem)->length, run->default_order);
            else
              n_left = elliptical_collimator(coord, (ECOL*)eptr->p_elem, n_to_track, accepted, z, *P_central);
            break;
          case T_CLEAN:
            if (!(flags&TEST_PARTICLES && !(flags&TEST_PARTICLE_LOSSES)))
              n_left = remove_outlier_particles(coord, (CLEAN*)eptr->p_elem, 
                                                n_to_track, accepted, z, *P_central);
            break;
          case T_SCRAPER:
            if (!(flags&TEST_PARTICLES && !(flags&TEST_PARTICLE_LOSSES))) {
              n_left = beam_scraper(coord, (SCRAPER*)eptr->p_elem, n_to_track, accepted, z, *P_central);
            }
            break;
          case T_PFILTER:
            if (!(flags&TEST_PARTICLES && !(flags&TEST_PARTICLE_LOSSES)))
              n_left = track_through_pfilter(coord, (PFILTER*)eptr->p_elem, n_to_track, 
                                             accepted, z, *P_central);
            break;
          case T_CENTER:
            center_beam(coord, (CENTER*)eptr->p_elem, n_to_track);
            break;
          case T_REMCOR:
            remove_correlations(coord, (REMCOR*)eptr->p_elem, n_to_track);
            break;
          case T_RFCA:
            n_left = simple_rf_cavity(coord, n_to_track, (RFCA*)eptr->p_elem, accepted, P_central, z);
            break;
          case T_RFCW:
            n_left = track_through_rfcw(coord, n_to_track, (RFCW*)eptr->p_elem, accepted, P_central, z,
                                        run, i_pass, charge);
            break;
          case T_MODRF:
            modulated_rf_cavity(coord, n_to_track, (MODRF*)eptr->p_elem, *P_central, z);
            break;
          case T_WATCH:
            watch_pt_seen = 1;
            if (!(flags&TEST_PARTICLES) && !(flags&INHIBIT_FILE_OUTPUT)) {
              watch = (WATCH*)eptr->p_elem;
              if (!watch->initialized) 
                set_up_watch_point(watch, run);
              if (i_pass==0 && (n_passes/watch->interval)==0)
                fprintf(stdout, "warning: n_passes = %ld and WATCH interval = %ld--no output will be generated!\n",
                        n_passes, watch->interval);
                fflush(stdout);
              if (i_pass>=watch->start_pass && (i_pass-watch->start_pass)%watch->interval==0) {
                switch (watch->mode_code) {
                case WATCH_COORDINATES:
                  dump_watch_particles(watch, step, i_pass, coord, n_to_track, *P_central,
                                       beamline->revolution_length);
                  break;
                case WATCH_PARAMETERS:
                case WATCH_CENTROIDS:
                  dump_watch_parameters(watch, step, i_pass, n_passes, coord, n_to_track, *n_original, *P_central,
                                        beamline->revolution_length);
                  break;
                case WATCH_FFT:
                  dump_watch_FFT(watch, step, i_pass, n_passes, coord, n_to_track, *n_original, *P_central);
                  break;
                }
              }
            }
            break;
          case T_HISTOGRAM:
            watch_pt_seen = 1;   /* yes, this should be here */
            if (!(flags&TEST_PARTICLES) && !(flags&INHIBIT_FILE_OUTPUT)) {
              histogram = (HISTOGRAM*)eptr->p_elem;
              if (!histogram->initialized) 
                set_up_histogram(histogram, run);
              if (i_pass==0 && (n_passes/histogram->interval)==0)
                fprintf(stdout, "warning: n_passes = %ld and HISTOGRAM interval = %ld--no output will be generated!\n",
                        n_passes, histogram->interval);
                fflush(stdout);
              if (i_pass>=histogram->startPass && (i_pass-histogram->startPass)%histogram->interval==0) {
                dump_particle_histogram(histogram, step, i_pass, coord, n_to_track, *P_central,
                                       beamline->revolution_length);
              }
            }
            break;
          case T_MALIGN:
            malign = (MALIGN*)eptr->p_elem;
            if (malign->on_pass==-1 || malign->on_pass==i_pass)
              offset_beam(coord, n_to_track, (MALIGN*)eptr->p_elem, *P_central);
            break;
          case T_PEPPOT:
            n_left = pepper_pot_plate(coord, (PEPPOT*)eptr->p_elem, n_to_track, accepted);
            break;
          case T_ENERGY:
            energy = (ENERGY*)eptr->p_elem;
            if (energy->match_beamline) {
              if ((flags&FIDUCIAL_BEAM_SEEN) && eptr->Pref_output_fiducial>0)
                /* Beamline momentum is defined.  Change particle reference momentum to match. */
                set_central_momentum(coord, n_to_track, eptr->Pref_output_fiducial, P_central);
              else
                /* Compute new central momentum to match the average momentum of the particles. */
                do_match_energy(coord, n_to_track, P_central, 0);
              if (energy->match_particles)
                bomb("can't match_beamline AND match_particles for ENERGY element", NULL);
            }
            else if (energy->match_particles) {
              /* change the particle momenta so that the centroid is the central momentum */
              do_match_energy(coord, n_to_track, P_central, 1);
            }
            else if (energy->central_energy)
              /* Change particle reference momentum to match the given energy */
              set_central_momentum(coord, n_to_track, sqrt(sqr(energy->central_energy+1)-1), 
                                   P_central);
            else if (energy->central_momentum)
              /* Change particle reference momentum to match the given value */
              set_central_momentum(coord, n_to_track, energy->central_momentum, P_central);
            break;
          case T_MAXAMP:
            maxamp = (MAXAMP*) eptr->p_elem;
            x_max = maxamp->x_max;
            y_max = maxamp->y_max;
            elliptical = maxamp->elliptical;
            break;
          case T_TRCOUNT:
            *n_original = n_left;
            if (accepted && i_pass==0)
              copy_particles(accepted, coord, *n_original);
            break;
          case T_ALPH:
            if (!eptr->matrix && !(eptr->matrix=compute_matrix(eptr, run, NULL)))
              bomb("no matrix for alpha magnet", NULL);
            n_left = alpha_magnet_tracking(coord, eptr->matrix, (ALPH*)eptr->p_elem, n_to_track,
                                           accepted, *P_central, z);
            break;
          case T_MULT:
            n_left = multipole_tracking(coord, n_to_track, (MULT*)eptr->p_elem, 0.0,
                                        *P_central, accepted, z);
            break;
          case T_FMULT:
            n_left = fmultipole_tracking(coord, n_to_track, (FMULT*)eptr->p_elem, 0.0,
                                         *P_central, accepted, z);
            break;
          case T_KICKER:
            if (flags&TIME_DEPENDENCE_OFF)
              drift_beam(coord, n_to_track, ((KICKER*)eptr->p_elem)->length, run->default_order);
            else
              track_through_kicker(coord, n_to_track, (KICKER*)eptr->p_elem, *P_central, i_pass, run->default_order);
            break;
          case T_KSBEND:
            n_left = track_through_kick_sbend(coord, n_to_track, (KSBEND*)eptr->p_elem, 0.0,
                                              *P_central, accepted, z);
            break;
          case T_CSBEND:
            if (flags&TEST_PARTICLES) {
              saveISR = ((CSBEND*)eptr->p_elem)->isr;
              ((CSBEND*)eptr->p_elem)->isr = 0;
            }
            n_left = track_through_csbend(coord, n_to_track, (CSBEND*)eptr->p_elem, 0.0,
                                          *P_central, accepted, last_z);
            if (flags&TEST_PARTICLES)
              ((CSBEND*)eptr->p_elem)->isr = saveISR;
            break;
          case T_CSRCSBEND:
            if (flags&TEST_PARTICLES) {
              saveISR = ((CSRCSBEND*)eptr->p_elem)->isr;
              ((CSRCSBEND*)eptr->p_elem)->isr = 0;
            }
            n_left = track_through_csbendCSR(coord, n_to_track, (CSRCSBEND*)eptr->p_elem, 0.0,
                                             *P_central, accepted, last_z, z, charge, run->rootname);
            if (flags&TEST_PARTICLES)
              ((CSRCSBEND*)eptr->p_elem)->isr = saveISR;
            break;
          case T_CSRDRIFT:
            n_left = track_through_driftCSR(coord, n_to_track, (CSRDRIFT*)eptr->p_elem,
                                            *P_central, accepted, last_z, run->rootname);
            break;
          case T_TUBEND:
            n_left = track_through_tubend(coord, n_to_track, 
                                          (TUBEND*)eptr->p_elem, 0.0,
                                          *P_central, accepted, z);
            break;
          case T_KQUAD:
          case T_KSEXT:
            n_left = multipole_tracking2(coord, n_to_track, eptr, 0.0,
                                         *P_central, accepted, z);
            break;
          case T_SAMPLE:
            if (!(flags&TEST_PARTICLES))
              n_left = sample_particles(coord, (SAMPLE*)eptr->p_elem, n_to_track, accepted, z, *P_central);
            break;
          case T_SCATTER:
            if (!(flags&TEST_PARTICLES))
              scatter(coord, n_to_track, *P_central, (SCATTER*)eptr->p_elem);
            break;
          case T_NIBEND:
            n_left = lorentz(coord, n_to_track, (NIBEND*)eptr->p_elem, T_NIBEND, *P_central, accepted);
            break;
          case T_NISEPT:
            n_left = lorentz(coord, n_to_track, (NISEPT*)eptr->p_elem, T_NISEPT, *P_central, accepted);
            break;
          case T_BMAPXY:
            n_left = lorentz(coord, n_to_track, (BMAPXY*)eptr->p_elem, T_BMAPXY, *P_central, accepted);
            break;
          case T_KPOLY:
            n_left = polynomial_kicks(coord, n_to_track, (KPOLY*)eptr->p_elem, 0.0,
                                      *P_central, accepted, z);
            break;
          case T_RAMPRF:
            ramped_rf_cavity(coord, n_to_track, (RAMPRF*)eptr->p_elem, *P_central, beamline->revolution_length,
                             eptr->end_pos, i_pass);
            break;
          case T_RAMPP:
            ramp_momentum(coord, n_to_track, (RAMPP*)eptr->p_elem, P_central, i_pass);
            break;
          case T_SOLE:
            if (((SOLE*)eptr->p_elem)->B) {
              SOLE *sptr;
              double ks;
              sptr = (SOLE*)eptr->p_elem;
              if ((ks = -sptr->B/(*P_central*me_mks*c_mks/e_mks))!=sptr->ks) {
                sptr->ks = ks;
                if (eptr->matrix)
                  free_matrices(eptr->matrix);
                if (!(eptr->matrix = compute_matrix(eptr, run, NULL)))
                  bomb("no matrix for element that must have matrix", NULL);
              }
              sptr->ks = 0;  /* reset so it is clear that B is fundamental quantity */
            }
            if (!eptr->matrix) {
              if (!(eptr->matrix=compute_matrix(eptr, run, NULL)))
                bomb("no matrix for element that must have matrix", NULL);
            }
            track_particles(coord, eptr->matrix, coord, n_to_track);
            break;
          case T_MATTER:
            track_through_matter(coord, n_to_track, (MATTER*)eptr->p_elem, *P_central);
            break;
          case T_RFMODE:
            rfmode = (RFMODE*)eptr->p_elem;
            if (!rfmode->initialized)
              set_up_rfmode(rfmode, eptr->name, eptr->end_pos, n_passes, run, 
                            *n_original, *P_central,
                            beamline->revolution_length);
            track_through_rfmode(coord, n_to_track, (RFMODE*)eptr->p_elem, *P_central,
                                 eptr->name, eptr->end_pos, i_pass, n_passes,
                                 charge);
            break;
          case T_TRFMODE:
            trfmode = (TRFMODE*)eptr->p_elem;
            if (!trfmode->initialized)
              set_up_trfmode(trfmode, eptr->name, eptr->end_pos, n_passes, run, *n_original);
            track_through_trfmode(coord, n_to_track, (TRFMODE*)eptr->p_elem, *P_central,
                                  eptr->name, eptr->end_pos, i_pass, n_passes,
                                  charge);
            break;
          case T_ZLONGIT:
            track_through_zlongit(coord, n_to_track, (ZLONGIT*)eptr->p_elem, *P_central, run, i_pass,
                                  charge);
            break;
          case T_ZTRANSVERSE:
            track_through_ztransverse(coord, n_to_track, (ZTRANSVERSE*)eptr->p_elem, *P_central, run, i_pass,
                                  charge);
            break;
          case T_WAKE:
            track_through_wake(coord, n_to_track, (WAKE*)eptr->p_elem, P_central, run, i_pass,
                               charge);
            break;
          case T_TRWAKE:
            track_through_trwake(coord, n_to_track, (TRWAKE*)eptr->p_elem, *P_central, run, i_pass, 
                                 charge);
            break;
          case T_SREFFECTS:
            if (!(flags&TEST_PARTICLES))
              track_SReffects(coord, n_to_track, (SREFFECTS*)eptr->p_elem, *P_central, eptr->twiss, &(beamline->radIntegrals));
            break;
          case T_IBSCATTER:
            if (!(flags&TEST_PARTICLES))
              track_IBS(coord, n_to_track, (IBSCATTER*)eptr->p_elem,
                        *P_central, &(beamline->elem), &(beamline->radIntegrals),
                        charge);
            break;
          default:
            fprintf(stdout, "programming error: no tracking statements for element %s (type %s)\n",
                    eptr->name, entity_name[eptr->type]);
            fflush(stdout);
            exit(1);
            break;
          }
        }
        if (!(flags&TEST_PARTICLES && !(flags&TEST_PARTICLE_LOSSES)) && (x_max || y_max)) {
          if (!elliptical) 
            n_left = limit_amplitudes(coord, x_max, y_max, n_left, accepted, z, *P_central, 
                                      eptr->type==T_DRIF || eptr->type==T_STRAY);
          else
            n_left = elimit_amplitudes(coord, x_max, y_max, n_left, accepted, z, *P_central, 
                                       eptr->type==T_DRIF || eptr->type==T_STRAY);
        }
        if (run->print_statistics && !(flags&TEST_PARTICLES)) {
          report_stats(stdout, ": ");
          /*
            if (show_dE && n_left) {
            fprintf(stdout, "average energy imparted: %e MeV\n",
            dgamma*me_mev/n_left);
            fflush(stdout);
            fprintf(stdout, "average x,y,z momentum imparted: %e, %e, %e MeV/c\n",
            dP[0]*me_mev/n_left, dP[1]*me_mev/n_left,
            dP[2]*me_mev/n_left);
            fflush(stdout);
            }
            */
          fprintf(stdout, "central momentum is %e    zstart = %em  zend = %em\n", *P_central, last_z, z);
          fflush(stdout);
          if (n_left!=n_to_track)
            fprintf(stdout, "%ld particles left\n", n_left);
            fflush(stdout);
        }
      }
      else if (!(flags&TEST_PARTICLES)) {
        fprintf(stdout, "element %s was ignored in tracking.\n",
                eptr->name);
        fflush(stdout);
      }
      if (flags&FIRST_BEAM_IS_FIDUCIAL && !(flags&FIDUCIAL_BEAM_SEEN)) {
        if (!(flags&RESTRICT_FIDUCIALIZATION) ||
            (entity_description[eptr->type].flags&MAY_CHANGE_ENERGY)) {
          do_match_energy(coord, n_left, P_central, 0);
        }
        eptr->Pref_output_fiducial = *P_central;
      } else if (flags&FIDUCIAL_BEAM_SEEN) {
        if (*P_central!=eptr->Pref_output_fiducial)
          set_central_momentum(coord, n_left, eptr->Pref_output_fiducial, P_central);
      }
      else if (run->always_change_p0)
        do_match_energy(coord, n_left, P_central, 0);
      if (i_pass==0 && traj_vs_z) {
        /* collect trajectory data--used mostly by trajectory correction routines */
        if (!traj_vs_z[i_traj].centroid) {
          fprintf(stdout, "error: the trajectory centroid array for %s is NULL (do_tracking)",
                  eptr->name);
          fflush(stdout);
          exit(1);
        }
        traj_vs_z[i_traj].elem = eptr;
        if (!(traj_vs_z[i_traj].n_part=n_left)) {
          for (i=0; i<6; i++)
            traj_vs_z[i_traj].centroid[i] = 0;
        }
        else {
          for (i=0; i<6; i++) {
            for (j=sum=0; j<n_to_track; j++)
              sum += coord[j][i];
            traj_vs_z[i_traj].centroid[i] = sum/n_left;
          }
        }
        i_traj++;
      }
      if (!(flags&TEST_PARTICLES)) {
	performSliceAnalysisOutput(sliceAnalysis, coord, n_to_track, 
				   !sliceAnDone, step, 
				   *P_central, 
				   charge?charge->macroParticleCharge*n_to_track:0.0, 
				   eptr->name, eptr->end_pos, 0);
	sliceAnDone = 1;
      }

      last_type = eptr->type;
      eptr = eptr->succ;
      n_to_track = n_left;
    }
    
    log_entry("do_tracking.2.2.3");
    if (effort)
      *effort += n_left;

    if (sums_vs_z && (*sums_vs_z) && !(flags&FINAL_SUMS_ONLY) && !(flags&TEST_PARTICLES) &&
        (run->wrap_around || i_pass==n_passes-1)) {
      if (i_sums<0)
        bomb("attempt to accumulate beam sums with negative index!", NULL);
      accumulate_beam_sums(*sums_vs_z+i_sums, coord, n_to_track, *P_central);
      (*sums_vs_z)[i_sums].z = z;
#if defined(BEAM_SUMS_DEBUG)
      fprintf(stdout, "beam sums accumulated in slot %ld for %s at z=%em, sx=%e\n", 
              i_sums, name, z, sqrt((*sums_vs_z)[i_sums].sum2[0]/n_left));
      fflush(stdout);
#endif
      i_sums++;
    }
    log_exit("do_tracking.2.2.3");
    
    log_exit("do_tracking.2.2");
#ifdef WATCH_MEMORY
    fprintf(stdout, "main tracking loop done: CPU: %6.2lf  PF: %6ld  MEM: %6ld\n",
            cpu_time()/100.0, page_faults(), memory_count());
    fflush(stdout);
#endif
    
    if (i_pass==0 || watch_pt_seen) {
      /* if eptr is not NULL, then all particles have been lost */
      /* some work still has to be done, however. */
      while (eptr) {
        if (entity_description[eptr->type].flags&HAS_LENGTH && eptr->p_elem)
          z += ((DRIFT*)eptr->p_elem)->length;
        else {
          if (eptr->pred)
            z += eptr->end_pos - eptr->pred->end_pos;
          else
            z += eptr->end_pos;
        }
        if (sums_vs_z && *sums_vs_z && !(flags&FINAL_SUMS_ONLY) && !(flags&TEST_PARTICLES)) {
          if (i_sums<0)
            bomb("attempt to accumulate beam sums with negative index!", NULL);
          accumulate_beam_sums(*sums_vs_z+i_sums, coord, n_to_track, *P_central);
          (*sums_vs_z)[i_sums].z = z;
          i_sums++;
        }
        switch (eptr->type) {
        case T_WATCH:
          if (!(flags&TEST_PARTICLES) && !(flags&INHIBIT_FILE_OUTPUT)) {
            watch = (WATCH*)eptr->p_elem;
            if (i_pass%watch->interval==0) {
              switch (watch->mode_code) {
              case WATCH_COORDINATES:
                break;
              case WATCH_PARAMETERS:
              case WATCH_CENTROIDS:
                dump_watch_parameters(watch, step, i_pass, n_passes, coord, n_to_track, *n_original, *P_central,
                                      beamline->revolution_length);
                break;
              case WATCH_FFT:
                dump_watch_FFT(watch, step, i_pass, n_passes, coord, n_to_track, *n_original, *P_central);
                break;
              }
            }
          }
          break;
        default:
          break;
        }
        if (i_pass==0 && traj_vs_z) {
          /* collect trajectory data--used mostly by trajectory correction routines */
          if (!traj_vs_z[i_traj].centroid) {
            fprintf(stdout, "error: the trajectory centroid array for %s is NULL (do_tracking)",
                    eptr->name);
            fflush(stdout);
            exit(1);
          }
          traj_vs_z[i_traj].elem = eptr;
          traj_vs_z[i_traj].n_part = 0;
          for (i=0; i<6; i++)
            traj_vs_z[i_traj].centroid[i] = 0;
          i_traj++;
        }
        eptr = eptr->succ;
      }
    }
  }

  /* do this here to get report of CSR drift normalization */
  reset_driftCSR();

  log_exit("do_tracking.2");
  log_entry("do_tracking.3");
  
  if (n_left && sums_vs_z && *sums_vs_z && !(flags&TEST_PARTICLES)) {
    if (flags&FINAL_SUMS_ONLY) {
      log_entry("do_tracking.3.1");
      i_sums = 0;
      accumulate_beam_sums(*sums_vs_z+i_sums, coord, n_to_track, *P_central);
      (*sums_vs_z)[i_sums].z = z;
#if defined(BEAM_SUMS_DEBUG)
      fprintf(stdout, "beam sums accumulated in slot %ld for final sums at z=%em, sx=%e\n", 
              i_sums, z, sqrt((*sums_vs_z)[i_sums].sum2[0]/n_left));
      fflush(stdout);
#endif
      log_exit("do_tracking.3.1");
    }
    else if (run->wrap_around) {
      log_entry("do_tracking.3.2");
      if (i_sums<0)
        bomb("attempt to accumulate beam sums with negative index!", NULL);
      /* accumulate sums for final output */
      accumulate_beam_sums(*sums_vs_z+i_sums, coord, n_to_track, *P_central);
#if defined(BEAM_SUMS_DEBUG)
      fprintf(stdout, "beam sums accumulated in slot %ld for final sums at z=%em, sx=%e\n", 
              i_sums, z, sqrt((*sums_vs_z)[i_sums].sum2[0]/n_left));
      fflush(stdout);
#endif
      log_exit("do_tracking.3.2");
    }
    else {
      log_entry("do_tracking.3.3");
      if (i_sums<0)
        bomb("attempt to accumulate beam sums with negative index!", NULL);
      copy_beam_sums(*sums_vs_z+i_sums, *sums_vs_z+i_sums-1);
#if defined(BEAM_SUMS_DEBUG)
      fprintf(stdout, "beam sums copied to slot %ld from slot %ld for final sums at z=%em, sx=%e\n", 
              i_sums, i_sums-1, z, (*sums_vs_z)[i_sums].sum2[0]);
      fflush(stdout);
#endif
      log_exit("do_tracking.3.3");
    }
  }

  if (sasefel && sasefel->active) {
    if (!charge) {
      fprintf(stdout, "Can't compute SASE FEL---no CHARGE element seen");
      fflush(stdout);
      exit(1);
    }
    computeSASEFELAtEnd(sasefel, coord, n_to_track, *P_central, charge->macroParticleCharge*n_to_track);
  }
  
  log_exit("do_tracking.3");
  log_entry("do_tracking.4");
  if (!(flags&SILENT_RUNNING) && !is_batch && n_passes!=1 && !(flags&TEST_PARTICLES)) {
    fprintf(stdout, "%ld particles left after pass %ld        \n", 
            n_to_track, i_pass);
    fflush(stdout);
  }

  log_exit("do_tracking.4");

  log_exit("do_tracking");
  if (charge && finalCharge)
    *finalCharge = n_to_track*charge->macroParticleCharge;
  return(n_to_track);
}

void offset_beam(
                 double **coord,
                 long n_to_track, 
                 MALIGN *offset,
                 double P_central
                 )
{
  long i_part;
  double *part, pc, beta, gamma, t;
  double ds;
  
  log_entry("offset_beam");
  
  for (i_part=n_to_track-1; i_part>=0; i_part--) {
    part = coord[i_part];
    if (offset->dz)
      ds = offset->dz*sqrt(1+sqr(part[1])+sqr(part[3]));
    else
      ds = 0;
    part[0] += offset->dx + offset->dz*part[1];
    part[1] += offset->dxp;
    part[2] += offset->dy + offset->dz*part[3];
    part[3] += offset->dyp;
    part[4] += ds;
    if (offset->dt || offset->dp || offset->de) {
      pc = P_central*(1+part[5]);
      beta = pc/(gamma=sqrt(1+pc*pc));
      t = part[4]/(beta*c_mks) + offset->dt;
      if (offset->dp) {
        part[5] += offset->dp;
        pc = P_central*(1+part[5]);
        beta = pc/sqrt(1+pc*pc);
      }
      if (offset->de) {
        gamma += offset->de*gamma;
        pc = sqrt(gamma*gamma-1);
        beta = pc/gamma;
        part[5] = (pc-P_central)/P_central;
      }
      part[4] = t*beta*c_mks;
    }
  }
  log_exit("offset_beam");
}

void do_match_energy(
                     double **coord, 
                     long np,
                     double *P_central,
                     long change_beam
                     )
{
  long ip;
  double P_average, dP_centroid, P, t;
  
  log_entry("do_match_energy");
  
  if (!np) {
    log_exit("do_match_energy");
    return;
  }
  
  if (!change_beam) {
    /* change the central momentum so that it matches the beam's centroid */
    P_average = 0;
    for (ip=0; ip<np; ip++) {
      P_average += (*P_central)*(1+coord[ip][5]);
    }
    P_average /= np;
    if (P_average!= *P_central) {
      for (ip=0; ip<np; ip++)
        coord[ip][5] = ((1+coord[ip][5])*(*P_central) - P_average)/ P_average;
      *P_central =  P_average;
    }
  }
  else {
    /* change the particle momenta so that the centroid is the central momentum */
    /* the path length is adjusted so that the time-of-flight at the current
       velocity is fixed */
    P_average = 0;
    for (ip=0; ip<np; ip++) 
      P_average += (*P_central*(1+coord[ip][5]));
    P_average /= np;
    dP_centroid =  *P_central - P_average;
    for (ip=0; ip<np; ip++) {
      P = (1+coord[ip][5])*(*P_central);
      t = coord[ip][4]/(P/sqrt(P*P+1));
      P += dP_centroid;
      coord[ip][5] = (P - *P_central)/ (*P_central);
      coord[ip][4] = t*(P/sqrt(P*P+1));
#if defined(IEEE_MATH)
      if (isnan(coord[ip][4]) || isinf(coord[ip][4])) {
        long i;
        fprintf(stdout, "error: bad time coordinate for particle %ld\n", ip);
        fflush(stdout);
        for (i=0; i<6; i++)
          fprintf(stdout, "%15.8e ", coord[ip][i]);
          fflush(stdout);
        fputc('\n', stdout);
        fprintf(stdout, "P_average = %e  P_central = %e  t = %e  dP_centroid = %e\n",
                P_average, *P_central, t, dP_centroid);
        fflush(stdout);
        abort();
      }
#endif
    }
  }

  log_exit("do_match_energy");
}

void set_central_energy(
                        double **coord, 
                        long np,
                        double new_energy,  /* new central gamma - 1*/
                        double *P_central
                        )
{
  
  log_entry("set_central_energy");
  set_central_momentum(coord, np, sqrt(sqr(new_energy+1)-1), P_central);
  log_exit("set_central_energy");
}

void set_central_momentum(
                          double **coord, 
                          long np,
                          double  P_new,  /* new central beta*gamma */
                          double *P_central
                          )
{
  long ip;
  
  if (!np) {
    *P_central =  P_new;
    return;
  }
  if (*P_central != P_new) {
    for (ip=0; ip<np; ip++)
      coord[ip][5] = ((1+coord[ip][5])*(*P_central) - P_new)/P_new;
    *P_central =  P_new;
  }
}

void remove_correlations(double **part, REMCOR *remcor, long np)
{
  double sumxy, sumy2, ratio;
  long ip, ic, wc;
  long removeFrom[4];
  
  if (!np) 
    return;
  
  removeFrom[0] = remcor->x;
  removeFrom[1] = remcor->xp;
  removeFrom[2] = remcor->y;
  removeFrom[3] = remcor->yp;
  wc = remcor->with-1;

  for (ip=sumy2=0; ip<np; ip++)
    sumy2 += part[ip][wc]*part[ip][wc];
  if (!sumy2)
    return;
  
  for (ic=0; ic<4; ic++) {
    if (!removeFrom[ic] || ic==wc)
      continue;
    if (!remcor->ratioSet[ic]) {
      for (ip=sumxy=0; ip<np; ip++)
        sumxy += part[ip][ic]*part[ip][wc];
      ratio = sumxy/sumy2;
      if (remcor->onceOnly) {
        remcor->ratio[ic] = ratio;
        remcor->ratioSet[ic] = 1;
      }
    }
    else 
      ratio = remcor->ratio[ic];
    for (ip=0; ip<np; ip++)
      part[ip][ic] -= part[ip][wc]*ratio;
  }
}


void center_beam(double **part, CENTER *center, long np)
{
  double sum, offset;
  long i, ic;
  long centerCoord[4];

  if (!np) {
    return;
  }

  centerCoord[0] = center->x;
  centerCoord[1] = center->xp;
  centerCoord[2] = center->y;
  centerCoord[3] = center->yp;

  for (ic=0; ic<4; ic++) {
    if (centerCoord[ic]) {
      if (!center->deltaSet[ic]) {
        for (i=sum=0; i<np; i++)
          sum += part[i][ic];
        center->delta[ic] = offset = sum/np;
        if (center->onceOnly)
          center->deltaSet[ic] = 1;
      } else 
        offset = center->delta[ic];
      for (i=0; i<np; i++)
        part[i][ic] -= offset;
    }
  }
}


void drift_beam(double **part, long np, double length, long order)
{
  VMATRIX *M;
  
  log_entry("drift_beam");
  
  if (length) {
    M = drift_matrix(length, order);
    track_particles(part, M, part, np);
    free_matrices(M);
    tfree(M);
    M = NULL;
  }
  log_exit("drift_beam");
}

void scatter(double **part, long np, double Po, SCATTER *scat)
{
  long i, ip;
  double t, P, beta;
  double sigma[4];

  if (!np)
    return;
  
  log_entry("scatter");
  sigma[0] = scat->x;
  sigma[1] = scat->xp;
  sigma[2] = scat->y;
  sigma[3] = scat->yp;
  for (i=0; i<4; i++) {
    if (!sigma[i])
      continue;
    for (ip=0; ip<np; ip++)
      part[ip][i] += gauss_rn(0, random_2)*sigma[i];
  }

  if (scat->dp) {
    for (ip=0; ip<np; ip++) {
      P = (1+part[ip][5])*Po;
      beta = P/sqrt(sqr(P)+1);
      t = part[ip][4]/beta;
      part[ip][5] += scat->dp*gauss_rn(0, random_2);
      P = (1+part[ip][5])*Po;
      beta = P/sqrt(sqr(P)+1);
      part[ip][4] = t*beta;
    }
  }

  log_exit("scatter");
}

void store_fitpoint_matrix_values(MARK *fpt, char *name, long occurence, VMATRIX *M)
{
  char buffer[1000];
  long i, j, k, count;

  if (!M) 
    return;
  if (!M->R)
    bomb("NULL R matrix passed to store_fitpoint_matrix_values", NULL);

  if (!(fpt->init_flags&8)) {
    if (M->order==1) {
      if (!(fpt->matrix_mem = malloc(sizeof(*(fpt->matrix_mem))*36)))
        bomb("memory allocation failure (store_fitpoint_matrix_values)", NULL);
    } else {
      if (!(fpt->matrix_mem = malloc(sizeof(*(fpt->matrix_mem))*(36+126))))
        bomb("memory allocation failure (store_fitpoint_matrix_values)", NULL);
    }
    for (i=count=0; i<6; i++) {
      for (j=0; j<6; j++) {
        sprintf(buffer, "%s#%ld.R%ld%ld", name, occurence, i+1, j+1);
        fpt->matrix_mem[count++] = rpn_create_mem(buffer);
      }
    }
    if (M->order>1) {
      for (i=0; i<6; i++) {
        for (j=0; j<6; j++) {
          for (k=0; k<=j; k++) {
            sprintf(buffer, "%s#%ld.T%ld%ld%ld", name, occurence, i+1, j+1, k+1);
            fpt->matrix_mem[count++] = rpn_create_mem(buffer);
          }
        }
      }
    }
    fpt->init_flags |= 8;
  }
  
  for (i=count=0; i<6; i++)
    for (j=0; j<6; j++)
      rpn_store(M->R[i][j], fpt->matrix_mem[count++]);
  if (M->order>1)
    for (i=0; i<6; i++)
      for (j=0; j<6; j++)
        for (k=0; k<=j; k++) 
          rpn_store(M->T[i][j][k], fpt->matrix_mem[count++]);
}


void store_fitpoint_beam_parameters(MARK *fpt, char *name, long occurence, double **coord, long np, double Po)
{
  long i;
  static double centroid[6], sigma[6], emit[3];
  static char *centroid_name_suffix[8] = {
    "Cx", "Cxp", "Cy", "Cyp", "Cs", "Cdelta", "pCentral", "Particles" };
  static char *sigma_name_suffix[6] = {
    "Sx", "Sxp", "Sy", "Syp", "Ss", "Sdelta" };
  static char *emit_name_suffix[3] = {
    "ex", "ey", "es"};
  static char s[1000];

  compute_centroids(centroid, coord, np);
  compute_sigmas(emit, sigma, centroid, coord, np);
  if (!(fpt->init_flags&2)) {
    fpt->centroid_mem = tmalloc(sizeof(*fpt->centroid_mem)*8);
    fpt->sigma_mem = tmalloc(sizeof(*fpt->sigma_mem)*6);
    fpt->emit_mem = tmalloc(sizeof(*fpt->emit_mem)*3);
    for (i=0; i<8; i++) {
      sprintf(s, "%s#%ld.%s", name, occurence, centroid_name_suffix[i]);
      fpt->centroid_mem[i] = rpn_create_mem(s);
    }
    for (i=0; i<6; i++) {
      sprintf(s, "%s#%ld.%s", name, occurence, sigma_name_suffix[i]);
      fpt->sigma_mem[i] = rpn_create_mem(s);
    }
    for (i=0; i<3; i++) {
      sprintf(s, "%s#%ld.%s", name, occurence, emit_name_suffix[i]);
      fpt->emit_mem[i] = rpn_create_mem(s);
    }
    fpt->init_flags |= 2;
  }
  for (i=0; i<6; i++) {
    rpn_store(centroid[i], fpt->centroid_mem[i]);
    rpn_store(sigma[i], fpt->sigma_mem[i]);
  }
  for (i=0; i<3; i++)
    rpn_store(emit[i], fpt->emit_mem[i]);
  rpn_store(Po, fpt->centroid_mem[6]);
  rpn_store((double)np, fpt->centroid_mem[7]);
}

void store_fitpoint_twiss_parameters(MARK *fpt, char *name, long occurence,TWISS *twiss)
{
  long i;
  static char *twiss_name_suffix[12] = {
    "betax", "alphax", "nux", "etax", "etapx", "etaxp",
    "betay", "alphay", "nuy", "etay", "etapy", "etaxp",
    } ;
  static char s[100];
  if (!(fpt->init_flags&1)) {
    fpt->twiss_mem = tmalloc(sizeof(*(fpt->twiss_mem))*12);
    fpt->init_flags |= 1;
    for (i=0; i<12; i++) {
      sprintf(s, "%s#%ld.%s", name, occurence, twiss_name_suffix[i]);
      fpt->twiss_mem[i] = rpn_create_mem(s);
    }
  }
  if (!twiss) {
    fprintf(stdout, "twiss parameter pointer unexpectedly NULL\n");
    fflush(stdout);
    abort();
  }
  for (i=0; i<5; i++) {
    rpn_store(*((&twiss->betax)+i)/(i==2?PIx2:1), fpt->twiss_mem[i]);
    rpn_store(*((&twiss->betay)+i)/(i==2?PIx2:1), fpt->twiss_mem[i+6]);
  }
  /* store etaxp and etayp in under two names each: etapx and etaxp */
  i = 4;
  rpn_store(*((&twiss->betax)+i), fpt->twiss_mem[i+1]);
  rpn_store(*((&twiss->betay)+i), fpt->twiss_mem[i+7]);
}


ELEMENT_LIST *findBeamlineMatrixElement(ELEMENT_LIST *eptr)
{
  ELEMENT_LIST *eptr0=NULL;
  long matrixSeen = 0;
  while (eptr) {
    if ((eptr->p_elem || eptr->matrix) && 
        entity_description[eptr->type].flags&MATRIX_TRACKING) {
      if (!(entity_description[eptr->type].flags&HAS_MATRIX))
        bomb("attempt to matrix track for element with no matrix!",  NULL);
      eptr0 = eptr;
      matrixSeen = 1;
      eptr = eptr->succ;
      break;
    }
    eptr = eptr->succ;
  }
  if (!matrixSeen)
    bomb("Can't do \"linear chromatic\" or \"longitudinal-only\" matrix tracking---no matrices!", NULL);
  while (eptr) {
    if ((eptr->p_elem || eptr->matrix) &&
        entity_description[eptr->type].flags&MATRIX_TRACKING) {
      bomb("Can't do \"linear chromatic\" or \"longitudinal-only\" matrix tracking---concatenation resulted in more than one matrix",
           NULL);
    }
    eptr = eptr->succ;
  }
  return eptr0;
}

void trackWithChromaticLinearMatrix(double **particle, long particles,
                                    ELEMENT_LIST *eptr,
                                    TWISS *twiss,
                                    double *tune0,
                                    double *chrom,
                                    double *dbeta_dPoP, 
                                    double *dalpha_dPoP)
{
  long ip, plane, offset, i, j;
  double *coord, deltaPoP, tune2pi, sin_phi, cos_phi;
  double alpha[2], beta[2], beta1, alpha1;
  double R11, R22, R12;
  static VMATRIX *M1 = NULL;
  double lastDPoP = DBL_MAX, det;
  if (!M1) {
    M1 = tmalloc(sizeof(*M1));
    initialize_matrices(M1, 1);
  }
  beta[0] = twiss->betax;
  beta[1] = twiss->betay;
  alpha[0] = twiss->alphax;
  alpha[1] = twiss->alphay;
  for (i=0; i<6; i++) {
    M1->C[i] = eptr->matrix->C[i];
    for (j=0; j<6; j++)
      M1->R[i][j] = eptr->matrix->R[i][j];
  }
  for (ip=0; ip<particles; ip++) {
    coord = particle[ip];
    deltaPoP = coord[5];
    if (deltaPoP!=lastDPoP) {
      for (plane=0; plane<2; plane++) {
        tune2pi = PIx2*(tune0[plane] + chrom[plane]*deltaPoP);
        offset = 2*plane;
        /* R11=R22 or R33=R44 */
        sin_phi = sin(tune2pi);
        cos_phi = cos(tune2pi);
        beta1 = beta[plane]+dbeta_dPoP[plane]*deltaPoP;
        alpha1 = alpha[plane]+dalpha_dPoP[plane]*deltaPoP;
        /* R11 or R33 */
        R11 = M1->R[0+offset][0+offset] = cos_phi + alpha1*sin_phi;
        /* R22 or R44 */
        R22 = M1->R[1+offset][1+offset] = cos_phi - alpha1*sin_phi;
        /* R12 or R34 */
        if ((R12 = M1->R[0+offset][1+offset] = beta1*sin_phi)) {
          /* R21 or R43 */
          M1->R[1+offset][0+offset] = 
            (R11*R22-1)/R12;
        }
        else {
          bomb("divided by zero in trackWithChromaticLinearMatrix", NULL);
        }
        det = M1->R[0+offset][0+offset]*M1->R[1+offset][1+offset] -
          M1->R[0+offset][1+offset]*M1->R[1+offset][0+offset];
        if (fabs(det-1)>1e-6) {
          fprintf(stdout, "Determinant is suspect: %e\n", det);
          fflush(stdout);
          exit(1);
        }
      }
    } 
    lastDPoP = deltaPoP;
    track_particles(&coord, M1, &coord, 1);
  }
}

void trackLongitudinalOnlyRing(double **part, long np, VMATRIX *M, double *alpha)
{
  long ip;
  double *coord, length, alpha1, alpha2;
  
  alpha1 = alpha[0];
  alpha2 = alpha[1];
  length = M->C[4];
  for (ip=0; ip<np; ip++) {
    coord = part[ip];
    coord[0] = coord[1] = coord[2] = coord[3] = 0;
    coord[4] += length*(1+(alpha1+alpha2*coord[5])*coord[5]);
  }
}

void getTrackingContext(TRACKING_CONTEXT *trackingContext0) 
{
  memcpy(trackingContext0, &trackingContext, sizeof(trackingContext));
}

