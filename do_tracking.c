/* Copyright 1994 by Michael Borland and Argonne National Laboratory,
 * all rights reserved.
 */
/* routine: do_tracking()
 * purpose: track a collection of particles through a beamline
 *
 * Michael Borland, 1989
 */
#include "mdb.h"
#include "track.h"
/* #include "smath.h" */

#pragma inline beta_from_delta

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
                 long flags,
                 long n_passes
                 )
{
    RFMODE *rfmode; TRFMODE *trfmode; ZLONGIT *zlongit;
    WATCH *watch;
    ENERGY *energy;
    MAXAMP *maxamp;
    MALIGN *malign;
    ELEMENT_LIST *eptr;
    long n_left, show_dE;
    double dgamma, dP[3], z, z_recirc, last_z;
    long i, j, i_traj=0, i_sums, n_to_track, i_pass, isConcat;
    long i_sums_recirc;
    long watch_pt_seen;
    double sum, x_max, y_max;
    long elliptical;
    double et1, et2;
    long is_batch = 0, last_type;
    static long is_ansi_term = -1;
    char s[100], *name;
    long check_nan, links_asserted=0, sums_allocated = 0;

    log_entry("do_tracking");
    
    log_entry("do_tracking.1");

#ifdef WATCH_MEMORY
    printf("start do_tracking():  CPU: %6.2lf  PF: %6ld  MEM: %6ld\n",
           cpu_time()/100.0, page_faults(), memory_count());
#endif
    
#if defined(UNIX)
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
    while (eptr && check_nan) {
        switch (eptr->type) {
          case T_RCOL: case T_ECOL: case T_MAXAMP:
            check_nan = 0;
            break;
          default:
            break;
            }
        eptr = eptr->succ;
        }

    log_exit("do_tracking.1");
    log_entry("do_tracking.2");
    name = "_BEG_";
    last_type = sums_allocated = 0;
    for (i_pass=0; i_pass<n_passes; i_pass++) {
        log_entry("do_tracking.2.1");
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
#if defined(UNIX)
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
            update_twiss_parameters(run, beamline);
            }
        if (run->concat_order && !(flags&TEST_PARTICLES) && !(beamline->flags&BEAMLINE_CONCAT_CURRENT) ) {
            /* form concatenated beamline for tracking */
            printf("Concatenating lattice into matrices of order %ld where possible.\n", run->concat_order);
            concatenate_beamline(beamline, run);
            }

        if (run->concat_order && beamline->flags&BEAMLINE_CONCAT_DONE && !(flags&TEST_PARTICLES)) {
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
        while (eptr && n_to_track) {
            log_entry("do_tracking.2.2.0");
            if (!eptr->name) {
                fprintf(stderr, "error: element ending at %em has NULL name pointer\n", eptr->end_pos);
                if (eptr->pred && eptr->pred->name)
                    fprintf(stderr, "previous element is %s\n", eptr->pred->name);
                else if (eptr->succ && eptr->succ->name)
                    fprintf(stderr, "next element is %s\n", eptr->succ->name);
                abort();
                }
            if (!eptr->p_elem && !run->concat_order) {
                printf("element %s has NULL p_elem pointer", eptr->name);
                exit(1);
                }
            if (eptr->type<=0 || eptr->type>=N_TYPES) {
                printf("element %s has type %ld--not recognized/not allowed\n", eptr->name, eptr->type);
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
                printf("beam sums accumulated in slot %ld for %s at z=%em, sx=%e\n", 
                       i_sums, name, z, sqrt((*sums_vs_z)[i_sums].sum2[0]/n_left));
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
            log_exit("do_tracking.2.2.1");
            if (eptr->p_elem || eptr->matrix) {
#ifdef VAX_VMS
                if (run->print_statistics && !(flags&TEST_PARTICLES))
                    printf("tracking through %s%c", eptr->name, ' ');
#else
                if (run->print_statistics && !(flags&TEST_PARTICLES))
                    printf("tracking through %s%c", eptr->name, '\n');
#endif
                n_left = n_to_track;
                show_dE = 0;
                if (entity_description[eptr->type].flags&MATRIX_TRACKING) {
                    if (!(entity_description[eptr->type].flags&HAS_MATRIX))
                        bomb("attempt to matrix-multiply for element with no matrix!",  NULL);
                    if (!eptr->matrix) {
                        if (!(eptr->matrix=compute_matrix(eptr, run, NULL)))
                            bomb("no matrix for element that must have matrix", NULL);
                        }
                    if (eptr->matrix->C[5]!=0)
                        bomb("matrix with C5!=0 detected in matrix multiplier--this shouldn't happen", NULL);
                    track_particles(coord, eptr->matrix, coord, n_to_track);
                    }
                else {
                    switch (eptr->type) {
                      case T_MARK:
                        if (((MARK*)eptr->p_elem)->fitpoint && i_pass==n_passes-1) {
                            if (beamline->flags&BEAMLINE_TWISS_WANTED)
                                store_fitpoint_twiss_parameters((MARK*)eptr->p_elem, eptr->name, eptr->occurence, eptr->twiss);
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
                      case T_RFTM:
                        bomb("RFTM not implemented yet", NULL);
                        break;
                      case T_TMCF:
                      case T_CEPL:
                      case T_TWPL:
                        if (!(flags&TIME_DEPENDENCE_OFF)) {
                            n_left = motion(coord, n_to_track, eptr->p_elem, eptr->type, *P_central, 
                                            &dgamma, dP, accepted, last_z);
                            show_dE = 1;
                            }
                        else
                            drift_beam(coord, n_to_track, ((TW_LINAC*)eptr->p_elem)->length, run->default_order);
                        break;
                      case T_TWLA:
                      case T_TWMTA:
                        n_left = motion(coord, n_to_track, eptr->p_elem, eptr->type, *P_central, 
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
                      case T_SCRAPER:
                        if (!(flags&TEST_PARTICLES && !(flags&TEST_PARTICLE_LOSSES))) {
                            n_left = beam_scraper(coord, (SCRAPER*)eptr->p_elem, n_to_track, accepted, z, *P_central);
                            }
                        break;
                      case T_CENTER:
                        center_beam(coord, (CENTER*)eptr->p_elem, n_to_track);
                        break;
                      case T_RFCA:
                        n_left = simple_rf_cavity(coord, n_to_track, (RFCA*)eptr->p_elem, accepted, P_central, z);
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
                                printf("warning: n_passes = %ld and WATCH interval = %ld--no output will be generated!\n",
                                       n_passes, watch->interval);
                            if (i_pass%watch->interval==0) {
                                switch (watch->mode_code) {
                                  case WATCH_COORDINATES:
                                    dump_watch_particles(watch, step, i_pass, coord, n_to_track, *P_central,
                                                         beamline->revolution_length);
                                    break;
                                  case WATCH_PARAMETERS:
                                  case WATCH_CENTROIDS:
                                    dump_watch_parameters(watch, step, i_pass, n_passes, coord, n_to_track, *n_original, *P_central);
                                    break;
                                  case WATCH_FFT:
                                    dump_watch_FFT(watch, step, i_pass, n_passes, coord, n_to_track, *n_original, *P_central);
                                    break;
                                    }
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
                            do_match_energy(coord, n_to_track, P_central, 0);
                            if (energy->match_particles)
                                bomb("can't match_beamline AND match_particles for ENERGY element", NULL);
                            }
                        else if (energy->match_particles) {
                            do_match_energy(coord, n_to_track, P_central, 1);
                            }
                        else if (energy->central_energy)
                            set_central_momentum(coord, n_to_track, sqrt(sqr(energy->central_energy+1)-1), 
                                                 P_central);
                        else if (energy->central_momentum)
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
                        n_left = track_through_csbend(coord, n_to_track, (CSBEND*)eptr->p_elem, 0.0,
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
                            set_up_rfmode(rfmode, eptr->name, eptr->end_pos, n_passes, run, *n_original, *P_central,
                                          beamline->revolution_length);
                        track_through_rfmode(coord, n_to_track, (RFMODE*)eptr->p_elem, *P_central,
                                             eptr->name, eptr->end_pos, i_pass, n_passes);
                        break;
                      case T_TRFMODE:
                        trfmode = (TRFMODE*)eptr->p_elem;
                        if (!trfmode->initialized)
                            set_up_trfmode(trfmode, eptr->name, eptr->end_pos, n_passes, run, *n_original);
                        track_through_trfmode(coord, n_to_track, (TRFMODE*)eptr->p_elem, *P_central,
                                              eptr->name, eptr->end_pos, i_pass, n_passes);
                        break;
                      case T_ZLONGIT:
                        track_through_zlongit(coord, n_to_track, (ZLONGIT*)eptr->p_elem, *P_central, run, i_pass);
                        break;
                      case T_SREFFECTS:
                        track_SReffects(coord, n_to_track, (SREFFECTS*)eptr->p_elem, *P_central, eptr->twiss);
                        break;
                      default:
                        fprintf(stderr, "programming error: no tracking statements for element %s (type %s)\n",
                                eptr->name, entity_name[eptr->type]);
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
                      printf("average energy imparted: %e MeV\n",
                      dgamma*me_mev/n_left);
                      printf("average x,y,z momentum imparted: %e, %e, %e MeV/c\n",
                      dP[0]*me_mev/n_left, dP[1]*me_mev/n_left,
                      dP[2]*me_mev/n_left);
                      }
                      */
                    printf("central momentum is %e    za = %em   ze = %em\n", *P_central, z, eptr->end_pos);
                    if (n_left!=n_to_track)
                        printf("%ld particles left\n", n_left);
                    fflush(stdout);
                    }
                }
            else if (!(flags&TEST_PARTICLES)) {
                printf("element %s was ignored in tracking.\n",
                       eptr->name);
                fflush(stdout);
                }
            if (i_pass==0 && traj_vs_z) {
                /* collect trajectory data--used mostly by trajectory correction routines */
                if (!traj_vs_z[i_traj].centroid) {
                    fprintf(stderr, "error: the trajectory centroid array for %s is NULL (do_tracking)",
                            eptr->name);
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
            printf("beam sums accumulated in slot %ld for %s at z=%em, sx=%e\n", 
                   i_sums, name, z, sqrt((*sums_vs_z)[i_sums].sum2[0]/n_left));
#endif
            i_sums++;
            }
        log_exit("do_tracking.2.2.3");
        
        log_exit("do_tracking.2.2");
#ifdef WATCH_MEMORY
        printf("main tracking loop done: CPU: %6.2lf  PF: %6ld  MEM: %6ld\n",
               cpu_time()/100.0, page_faults(), memory_count());
#endif
        
        if (i_pass==0 || watch_pt_seen) {
            while (eptr) {
                /* some elements need to be executed even if there are no particles */
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
                                dump_watch_parameters(watch, step, i_pass, n_passes, coord, n_to_track, *n_original, *P_central);
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
                        fprintf(stderr, "error: the trajectory centroid array for %s is NULL (do_tracking)",
                                eptr->name);
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
    log_exit("do_tracking.2");
    log_entry("do_tracking.3");
    
    if (n_left && sums_vs_z && *sums_vs_z && !(flags&TEST_PARTICLES)) {
        if (flags&FINAL_SUMS_ONLY) {
            log_entry("do_tracking.3.1");
            i_sums = 0;
            accumulate_beam_sums(*sums_vs_z+i_sums, coord, n_to_track, *P_central);
            (*sums_vs_z)[i_sums].z = z;
#if defined(BEAM_SUMS_DEBUG)
            printf("beam sums accumulated in slot %ld for final sums at z=%em, sx=%e\n", 
                   i_sums, z, sqrt((*sums_vs_z)[i_sums].sum2[0]/n_left));
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
            printf("beam sums accumulated in slot %ld for final sums at z=%em, sx=%e\n", 
                   i_sums, z, sqrt((*sums_vs_z)[i_sums].sum2[0]/n_left));
#endif
            log_exit("do_tracking.3.2");
            }
        else {
            log_entry("do_tracking.3.3");
            if (i_sums<0)
                bomb("attempt to accumulate beam sums with negative index!", NULL);
            copy_beam_sums(*sums_vs_z+i_sums, *sums_vs_z+i_sums-1);
#if defined(BEAM_SUMS_DEBUG)
            printf("beam sums copied to slot %ld from slot %ld for final sums at z=%em, sx=%e\n", 
                   i_sums, i_sums-1, z, (*sums_vs_z)[i_sums].sum2[0]);
#endif
            log_exit("do_tracking.3.3");
            }
        }
    
    log_exit("do_tracking.3");
    log_entry("do_tracking.4");
    if (!(flags&SILENT_RUNNING) && !is_batch && n_passes!=1 && !(flags&TEST_PARTICLES)) {
        printf("%ld particles left after pass %ld        \n", 
               n_to_track, i_pass);
        fflush(stdout);
        }

    log_exit("do_tracking.4");

    log_exit("do_tracking");
    return(n_to_track);
    }

void do_element_misalignment(ELEMENT_LIST *elem, double **coord, long n, long mode)
{
    QUAD *quad; BEND *bend; SEXT *sext; 
    ECOL *ecol; RCOL *rcol; ALPH *alph;
    SOLE *sole; QFRING *qfring;
    double dx, dy;
    long i;
    
    log_entry("do_element_misalignment");
    
    dx = dy = 0;
    
    switch (elem->type) {
      case T_RBEN: case T_SBEN:
        bend = (BEND*)elem->p_elem;
        dx = bend->dx;
        dy = bend->dy;
        break;
      case T_QUAD:
        quad = (QUAD*)elem->p_elem;
        dx = quad->dx;
        dy = quad->dy;
        break;
      case T_SEXT:
        sext = (SEXT*)elem->p_elem;
        dx   = sext->dx;
        dy   = sext->dy;
        break;
      case T_ALPH:
        alph = (ALPH*)elem->p_elem;
        dx   = alph->dx;
        dy   = alph->dy;
        break;  
      case T_SOLE:
        sole = (SOLE*) elem->p_elem;
        dx = sole->dx;
        dy = sole->dy;
        break;
      case T_ECOL:
        ecol = (ECOL*) elem->p_elem;
        dx = ecol->dx;
        dy = ecol->dy;
        break;
      case T_RCOL:
        rcol = (RCOL*) elem->p_elem;
        dx = rcol->dx;
        dy = rcol->dy;
        break;
      case T_QFRING:
        qfring = (QFRING*)elem->p_elem;
        dx = qfring->dx;
        dy = qfring->dy;
        break;
      default:
        break;
        }
    
    if (dx) {
        if (mode<0)
            dx = -dx;
        for (i=0; i<n; i++)
            coord[i][0] -= dx;
        }
    
    if (dy) {
        if (mode<0)
            dy = -dy;
        for (i=0; i<n; i++)
            coord[i][2] -= dy;
        }
    
    log_exit("do_element_misalignment");
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
        for (ip=0; ip<np; ip++) 
            P_average += (coord[ip][5] = *P_central*(1+coord[ip][5]));
        P_average /= np;
        for (ip=0; ip<np; ip++)
            coord[ip][5] = (coord[ip][5]- P_average)/ P_average;
        *P_central =  P_average;
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
                printf("error: bad time coordinate for particle %ld\n", ip);
                for (i=0; i<6; i++)
                    printf("%15.8e ", coord[ip][i]);
                putchar('\n');
                printf("P_average = %e  P_central = %e  t = %e  dP_centroid = %e\n",
                       P_average, P_central, t, dP_centroid);
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
    
    log_entry("set_central_momentum");
    
    if (!np) {
        *P_central =  P_new;
        log_exit("set_central_momentum");
        return;
        }
    
    for (ip=0; ip<np; ip++)
        coord[ip][5] = (1+coord[ip][5])*(*P_central)/ P_new-1;
    
    *P_central =  P_new;
    log_exit("set_central_momentum");
    }

void center_beam(double **part, CENTER *center, long np)
{
    double sum, offset;
    long i;
    
    log_entry("center_beam");
    
    if (!np) {
        log_exit("center_beam");
        return;
        }
    
    if (center->x) {
        for (i=sum=0; i<np; i++)
            sum += part[i][0];
        offset = sum/np;
        for (i=0; i<np; i++)
            part[i][0] -= offset;
        }
    
    if (center->xp) {
        for (i=sum=0; i<np; i++)
            sum += part[i][1];
        offset = sum/np;
        for (i=0; i<np; i++)
            part[i][1] -= offset;
        }
    
    if (center->y) {
        for (i=sum=0; i<np; i++)
            sum += part[i][2];
        offset = sum/np;
        for (i=0; i<np; i++)
            part[i][2] -= offset;
        }
    
    if (center->yp) {
        for (i=sum=0; i<np; i++)
            sum += part[i][3];
        offset = sum/np;
        for (i=0; i<np; i++)
            part[i][3] -= offset;
        }
    log_exit("center_beam");
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

void store_fitpoint_beam_parameters(MARK *fpt, char *name, long occurence, double **coord, long np, double Po)
{
    long i;
    static double centroid[6], sigma[6];
    static char *centroid_name_prefix[8] = {
        "Cx", "Cxp", "Cy", "Cyp", "Cs", "Cdelta", "pCentral", "Particles" };
    static char *sigma_name_prefix[6] = {
        "Sx", "Sxp", "Sy", "Syp", "Ss", "Sdelta" };
    static char s[100];

    compute_centroids(centroid, coord, np);
    compute_sigmas(sigma, centroid, coord, np);
    if (!(fpt->init_flags&2)) {
        fpt->centroid_mem = tmalloc(sizeof(*fpt->centroid_mem)*8);
        fpt->sigma_mem = tmalloc(sizeof(*fpt->sigma_mem)*6);
        for (i=0; i<8; i++) {
            sprintf(s, "%s#%ld.%s", name, occurence, centroid_name_prefix[i]);
            fpt->centroid_mem[i] = rpn_create_mem(s);
            }
        for (i=0; i<6; i++) {
            sprintf(s, "%s#%ld.%s", name, occurence, sigma_name_prefix[i]);
            fpt->sigma_mem[i] = rpn_create_mem(s);
            }
        fpt->init_flags |= 2;
        }
    for (i=0; i<6; i++) {
        rpn_store(centroid[i], fpt->centroid_mem[i]);
        rpn_store(sigma[i], fpt->sigma_mem[i]);
        }
    rpn_store(Po, fpt->centroid_mem[6]);
    rpn_store((double)np, fpt->centroid_mem[7]);
    }

void store_fitpoint_twiss_parameters(MARK *fpt, char *name, long occurence,TWISS *twiss)
{
    long i;
    static char *twiss_name_prefix[10] = {
        "betax", "alphax", "nux", "etax", "etapx",
        "betay", "alphay", "nuy", "etay", "etapy"
        } ;
    static char s[100];
    if (!(fpt->init_flags&1)) {
        fpt->twiss_mem = tmalloc(sizeof(*(fpt->twiss_mem))*10);
        fpt->init_flags |= 1;
        for (i=0; i<10; i++) {
            sprintf(s, "%s#%ld.%s", name, occurence, twiss_name_prefix[i]);
            fpt->twiss_mem[i] = rpn_create_mem(s);
            }
        }
    if (!twiss) {
        fprintf(stderr, "twiss parameter pointer unexpectedly NULL\n");
        abort();
        }
    for (i=0; i<5; i++) {
        rpn_store(*((&twiss->betax)+i)/(i==2?PIx2:1), fpt->twiss_mem[i]);
        rpn_store(*((&twiss->betay)+i)/(i==2?PIx2:1), fpt->twiss_mem[i+5]);
        }
    }


