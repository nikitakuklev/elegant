/*************************************************************************\
 * Copyright (c) 2002 The University of Chicago, as Operator of Argonne
 * National Laboratory.
 * Copyright (c) 2002 The Regents of the University of California, as
 * Operator of Los Alamos National Laboratory.
 * This file is distributed subject to a Software License Agreement found
 * in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* routine: do_tracking()
 * purpose: track a collection of particles through a beamline
 *
 * Michael Borland, 1989
 */
#include "mdb.h"
#include "mdbsun.h"
#include "track.h"
/* #include "smath.h" */
void set_up_frfmode(FRFMODE *rfmode, char *element_name, double element_z, long n_passes,  RUN *run, long n_particles, double Po, double total_length);
void track_through_frfmode(double **part, long np, FRFMODE *rfmode, double Po,char *element_name, double element_z, long pass, long n_passes,CHARGE *charge);
void set_up_ftrfmode(FTRFMODE *rfmode, char *element_name, double element_z, long n_passes,RUN *run, long n_particles,double Po, double total_length);
void track_through_ftrfmode(double **part, long np, FTRFMODE *trfmode, double Po,char *element_name, double element_z, long pass, long n_passes,CHARGE *charge);

ELEMENT_LIST *findBeamlineMatrixElement(ELEMENT_LIST *eptr);
void trackLongitudinalOnlyRing(double **part, long np, VMATRIX *M, double *alpha);
void store_fitpoint_matrix_values(MARK *fpt, char *name, long occurence, VMATRIX *M);
long trackWithChromaticLinearMatrix(double **particle, long particles,
                                    double **accepted, double Po, double z,
                                    ELEMENT_LIST *eptr,
                                    TWISS *twiss0, double *tune0,
                                    double *chrom, double *chrom2, double *chrom3,
                                    double *dbeta_dPoP, double *dalpha_dPoP,
                                    double *alphac, double *eta2, double *eta3);
void matr_element_tracking(double **coord, VMATRIX *M, MATR *matr,
                           long np, double z);
void ematrix_element_tracking(double **coord, VMATRIX *M, EMATRIX *matr,
                              long np, double z);
long transformBeamWithScript(SCRIPT *script, double pCentral, CHARGE *charge, BEAM *beam, double **part, 
                             long np, long nLost, char *mainRootname, long iPass, long driftOrder);
void distributionScatter(double **part, long np, double Po, DSCATTER *scat, long iPass);
void recordLossPass(long *lostOnPass, long *nLost, long nLeft, long nMaximum, long pass, int myid,
                    int lostSinceSeqMode);
void storeMonitorOrbitValues(ELEMENT_LIST *eptr, double **part, long np);

#if USE_MPI
typedef enum pMode {notParallel, initialMode, trueParallel} parallelMode;
typedef enum balanceMode {badBalance, startMode, goodBalance} balance;
void scatterParticles(double **coord, long *nToTrack, double **accepted,
                      long n_processors, int myid, balance balanceStatus, 
                      double my_rate, double nParPerElements, double round,
                      int lostSinceSeqMod,int *distributed);
void gatherParticles(double **coord, long *lostOnPass, long *nToTrack, 
                     long *nLost, double **accepted, long n_processors, 
                     int myid, double *round);
balance checkBalance(double my_wtime, int myid, long n_processors);
#endif

#ifdef SORT   
int comp_IDs(const void *coord1, const void *coord2);
#endif
static TRACKING_CONTEXT trackingContext;

double beta_from_delta(double p, double delta)
{
  p *= 1+delta;
  return( p/sqrt(p*p+1));
}

long do_tracking(
                 /* Either the beam pointer or the coord pointer must be supplied, but not both */
                 BEAM *beam,  
                 double **coord,
                 long nOriginal,   /* Used only if coord is supplied */
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
                 long passOffset,
                 SASEFEL_OUTPUT *sasefel,
		 SLICE_OUTPUT *sliceAnalysis,
                 double *finalCharge,
		 long *lostOnPass
                 )
{
  RFMODE *rfmode; TRFMODE *trfmode;
  FRFMODE *frfmode; FTRFMODE *ftrfmode;
  WATCH *watch;
  STRAY *stray;
  HISTOGRAM *histogram;
  ENERGY *energy;
  MAXAMP *maxamp;
  MALIGN *malign;
  ELEMENT_LIST *eptr, *eptrPred, *eptrCLMatrix=NULL;
  long nToTrack;  /* number of particles being tracked */
  long nLeft;     /* number of those that are left after a tracking routine returns */
  long nLost=0;     /* accumulated number lost */
  long nMaximum=0;  /* maximum number of particles seen */
  long show_dE, maxampOpenCode=0, maxampExponent=0;
  double dgamma, dP[3], z, z_recirc, last_z;
  long i, j, i_traj=0, i_sums, i_pass, isConcat;
  long i_sums_recirc, saveISR=0;
  long watch_pt_seen, feedbackDriverSeen;
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
  unsigned long classFlags;
  long nParticlesStartPass = 0;
  int myid, active = 1, lostSinceSeqMode = 0;
#if USE_MPI 
  long old_nToTrack = 0, nParElements, nElements = 0;
  long n_processors = run->n_processors; 
  int checkFlags;
  double my_wtime, start_wtime, end_wtime, nParPerElements, my_rate;
  double round = 0.5;
  parallelMode parallelStatus; 
  parallelStatus = initialMode;
  balance balanceStatus;
  balanceStatus = startMode;
  int distributed = 0; /* identify if the particles have been scattered */ 
  MPI_Comm_rank(MPI_COMM_WORLD, &myid); /* get ID number for each processor */
  trackingContext.myid = myid;
#endif 

  if (!coord && !beam)
    bomb("Null particle coordinate array and null beam pointer! (do_tracking)", NULL);
  if (coord && beam)
    bomb("Particle coordinate array and beam pointer both supplied!  (do_tracking)", NULL);
  if (beam) {
    coord = beam->particle;
    nOriginal = beam->n_to_track;  /* used only for computing macroparticle charge */
  }
  
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
    copy_particles(accepted, coord, nOriginal);

#ifdef VAX_VMS
  is_batch = job_mode(getpid())==2?1:0;
#endif
  
  z = z_recirc = last_z =  0;
  i_sums = i_sums_recirc = 0;
  x_max = y_max = 0;
  nToTrack = nLeft = nMaximum = nOriginal;
  et1 = -2.0;
  elliptical = isConcat = 0;
  watch_pt_seen = feedbackDriverSeen = 0;

#ifdef SORT
  int nToTrackAtLastSort = nToTrack;
#endif
  
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

  if (check_nan) {
      nLeft = nToTrack = limit_amplitudes(coord, DBL_MAX, DBL_MAX, nToTrack, accepted, z, *P_central, 0,
					  0);
      recordLossPass(lostOnPass, &nLost, nLeft, nMaximum, 0, myid, lostSinceSeqMode);
  }

#if USE_MPI
  if (nToTrack<(n_processors-1)) {
    printf("************************************************************************\n");
    printf("* The number of particles can't be less than the number of processors! *\n");
    printf("************************************************************************\n");
    MPI_Barrier (MPI_COMM_WORLD);
    MPI_Abort (MPI_COMM_WORLD, 0);
  }
#endif

  for (i_pass=passOffset; i_pass<n_passes+passOffset; i_pass++) {
    log_entry("do_tracking.2.1");

    ResetNoiseGroupValues();

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

    if (i_pass==0) {
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
    }
      
    if (sums_vs_z && n_z_points) {
      if (!sums_allocated && !*sums_vs_z) {
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
    if (run->final_pass && sums_vs_z && n_z_points)
      zero_beam_sums(*sums_vs_z, *n_z_points+1);

    log_exit("do_tracking.2.1");
    log_entry("do_tracking.2.2");
    if (!(flags&SILENT_RUNNING) && !is_batch && n_passes!=1 && !(flags&TEST_PARTICLES)
        && !(run->tracking_updates==0)) {
#if defined(VAX_VMS)
      sprintf(s, "%ld particles present after pass %ld        ",
              nToTrack, i_pass);
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
        sprintf(s, "%ld particles present after pass %ld        ", 
                nToTrack, i_pass);
        fputs(s, stdout);
        if (is_ansi_term)
          backspace(strlen(s));
        else
          fputc('\n', stdout);
        fflush(stdout);
        et1 = et2;
      }
#else
      sprintf(s, "%ld particles present after pass %ld        ", 
              nToTrack, i_pass);
      fputs(s, stdout);
      if (is_ansi_term)
        backspace(strlen(s));
      else
        fputc('\n', stdout);
      fflush(stdout);
#endif 
    }
    elementsTracked = -1;
    eptrPred = eptr;
#if USE_MPI
    my_wtime = 0.0;
    nParElements = 0;
#endif
    nParticlesStartPass = nToTrack;

    if (getSCMULTSpecCount()) 
      /* prepare space charge effects calculation  */
      initializeSCMULT(eptr, coord, nToTrack, *P_central, i_pass);

    while (eptr && (nToTrack || USE_MPI)) {
      classFlags = entity_description[eptr->type].flags;
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
        accumulate_beam_sums(*sums_vs_z+i_sums, coord, nToTrack, *P_central);
        (*sums_vs_z)[i_sums].z = z;
#if defined(BEAM_SUMS_DEBUG)
        fprintf(stdout, "beam sums accumulated in slot %ld for %s at z=%em, sx=%e\n", 
                i_sums, name, z, sqrt((*sums_vs_z)[i_sums].sum2[0]/nLeft));
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
      strncpy(trackingContext.elementName, eptr->name, CONTEXT_BUFSIZE);
      trackingContext.elementOccurrence = eptr->occurence;
      trackingContext.sliceAnalysis = sliceAnalysis?
	(sliceAnalysis->finalValuesOnly?NULL:sliceAnalysis):NULL;
      trackingContext.zStart = last_z;
      trackingContext.zEnd = z;
      trackingContext.step = step;
      strncpy(trackingContext.rootname, run->rootname, CONTEXT_BUFSIZE);

#if USE_MPI
      active = 0;
      if (classFlags&UNIPROCESSOR) {
        /* This element cannot be done in parallel. Only the master CPU will work. */
        if (myid == 0)
          active = 1;
        else 
          active = 0;
        if (parallelStatus==trueParallel) {
          gatherParticles(coord, lostOnPass, &nToTrack, &nLost, accepted, 
                          n_processors, myid, &round);

          /* update the nMaximum for recording the nLost on all the slave processors */
          if (myid!=0)
          nMaximum = nToTrack;

          /* The element will change the state of particles. Scatter is required
             for parallel computation */
          if (!(classFlags&(UNIDIAGNOSTIC&(~UNIPROCESSOR)))) {
            parallelStatus = notParallel;           
          }
        }    
      } 
      else {
        /* This element can be done in parallel. Only the slave CPUS will work. */
        if (myid != 0)
          active = 1; 
        else
          active = 0;
        nElements++;
        if ((balanceStatus==badBalance) && (parallelStatus==trueParallel)) {
          gatherParticles(coord, lostOnPass, &nToTrack, &nLost, accepted, n_processors, myid, &round);
          if (myid != 0) {
            /* update the nMaximum for recording the nLost on all the slave processors */
            nMaximum = nToTrack;
          }
        } 
        /* Particles will be scattered in startMode, bad balancing status or notParallel state */  
        if ((balanceStatus==badBalance) || (parallelStatus==notParallel)) { 
          scatterParticles(coord, &nToTrack, accepted, n_processors, myid,
                           balanceStatus, my_rate, nParPerElements, round, lostSinceSeqMode, &distributed);
          if (myid != 0) {
           /* update the nMaximum for recording the nLost on all the slave processors */
           nMaximum = nToTrack;  
          }
          if (balanceStatus!=startMode)
            balanceStatus = goodBalance;
        }
        else if (balanceStatus==startMode) { 
	  /* For the first pass, scatter when it is not in parallel mode */
	  if (parallelStatus!=trueParallel) {
	    scatterParticles(coord, &nToTrack, accepted, n_processors, myid,
                             balanceStatus, my_rate, nParPerElements, round, lostSinceSeqMode, &distributed); 
            if (myid != 0) {
              /* update the nMaximum for recording the nLost on all the slave processors */
              nMaximum = nToTrack; 
            }
          }
	}
	parallelStatus = trueParallel;         
        if (myid != 0) {
          /* count the total number of particles tracked by all of the elements for each pass */
          nParElements = nParElements+nToTrack;
        }
        lostSinceSeqMode = 0;
        start_wtime = MPI_Wtime();        
      }
#endif

#ifdef SORT
      if (!USE_MPI)
	if (nToTrackAtLastSort != nToTrack) {/* indicates more particles are lost, need sort */
          qsort(coord[0], nToTrack, COORDINATES_PER_PARTICLE*sizeof(double), comp_IDs);
          if (accepted!=NULL)
            qsort(accepted[0], nToTrack, COORDINATES_PER_PARTICLE*sizeof(double), comp_IDs);
	  nToTrackAtLastSort = nToTrack;
	}   
#endif   
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
        show_dE = 0;
        nLeft = nToTrack;  /* in case it isn't set by the element tracking */
        if (eptr==eptrCLMatrix) {
          /* This element is the place-holder for the chromatic linear matrix or
           * the longitudinal-only matrix 
           */
          if (!USE_MPI || (USE_MPI && (myid!=0))) {
            /* Only the slave CPUs will work on this part */ 
            if (flags&LINEAR_CHROMATIC_MATRIX) 
              nLeft
	        = trackWithChromaticLinearMatrix(coord, nToTrack, accepted,
		       			         *P_central, z, eptrCLMatrix,
					         beamline->twiss0, beamline->tune,
					         beamline->chromaticity,
					         beamline->chrom2, beamline->chrom3,
					         beamline->dbeta_dPoP, beamline->dalpha_dPoP,
					         beamline->alpha, beamline->eta2, beamline->eta3);
            else 
              trackLongitudinalOnlyRing(coord, nToTrack, 
                                        eptrCLMatrix->matrix,
                                        beamline->alpha);
          }
	}
        else if (entity_description[eptr->type].flags&MATRIX_TRACKING &&
		 !(flags&IBS_ONLY_TRACKING)) {
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
            nLeft = 0;
          } else {
            if (run->print_statistics>1 && !(flags&TEST_PARTICLES)) {
              fprintf(stdout, "Tracking matrix for %s\n", eptr->name);
              fflush(stdout);
              print_elem(stdout, eptr);
              print_matrices(stdout, "", eptr->matrix);
            }
            if (flags&CLOSED_ORBIT_TRACKING) {
              switch (eptr->type) {
              case T_MONI:
              case T_HMON:
              case T_VMON:
                storeMonitorOrbitValues(eptr, coord, nToTrack);
                break;
              default:
                break;
              }
            }
            if (!USE_MPI || (USE_MPI && (myid!=0)))  /* Only the slave CPUs will track */ 
              track_particles(coord, eptr->matrix, coord, nToTrack);
          }
        }
        else {
          long type;
          if (run->print_statistics>1 && !(flags&TEST_PARTICLES)) {
            fprintf(stdout, "Tracking element: ");
            fflush(stdout);
            print_elem(stdout, eptr);
          }
	  type = eptr->type;
	  if (flags&IBS_ONLY_TRACKING) {
	    switch (type) {
	    case T_IBSCATTER:
	    case T_WATCH:
	    case T_CLEAN:
	    case T_RCOL:
	    case T_CHARGE:
	      break;
	    default:
	      type = -1;
	      break;
	    }
	  }
          
	  if (active && ((!USE_MPI && nParticlesStartPass) || nToTrack || 
              entity_description[eptr->type].flags&RUN_ZERO_PARTICLES)) {
	    switch (type) {
	    case -1:
	      break;
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
		if (nOriginal)
		  charge->macroParticleCharge = charge->charge/(nOriginal);
		if (charge->chargePerParticle)
		  charge->macroParticleCharge = charge->chargePerParticle;
	      }
	      break;
	    case T_MARK:
	      if (((MARK*)eptr->p_elem)->fitpoint && i_pass==n_passes-1) {
		/*
		  if (beamline->flags&BEAMLINE_TWISS_WANTED) {
		  if (!(beamline->flags&BEAMLINE_TWISS_DONE))
                  update_twiss_parameters(run, beamline, NULL);
		  store_fitpoint_twiss_parameters((MARK*)eptr->p_elem, eptr->name, 
		  eptr->occurence, eptr->twiss);
		  }
		*/
		store_fitpoint_matrix_values((MARK*)eptr->p_elem, eptr->name, 
					     eptr->occurence, eptr->accumMatrix);
		store_fitpoint_beam_parameters((MARK*)eptr->p_elem, eptr->name,eptr->occurence, 
					       coord, nToTrack, *P_central); 
		if (flags&CLOSED_ORBIT_TRACKING)
		  storeMonitorOrbitValues(eptr, coord, nToTrack);
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
					   coord, nToTrack, *P_central);
	      else
		drift_beam(coord, nToTrack, ((RFDF*)eptr->p_elem)->length, run->default_order);
	      break;
	    case T_RFTM110:
	      if (!(flags&TIME_DEPENDENCE_OFF))
		track_through_rftm110_deflector(coord, (RFTM110*)eptr->p_elem,
						coord, nToTrack, *P_central,
						beamline->revolution_length, eptr->end_pos,
						i_pass);
	      break;
	    case T_RMDF:
	      if (!(flags&TIME_DEPENDENCE_OFF))
		track_through_ramped_deflector(coord, (RMDF*)eptr->p_elem,
					       coord, nToTrack, *P_central);
	      else
		drift_beam(coord, nToTrack, ((RMDF*)eptr->p_elem)->length, run->default_order);
	      break;
	    case T_RFTMEZ0:
              nLeft = motion(coord, nToTrack, eptr->p_elem, eptr->type, P_central, 
			     &dgamma, dP, accepted, last_z);
              show_dE = 1;
	      break;
	    case T_TMCF:
	    case T_CEPL:
	    case T_TWPL:
	      if (!(flags&TIME_DEPENDENCE_OFF)) {
		nLeft = motion(coord, nToTrack, eptr->p_elem, eptr->type, P_central, 
			       &dgamma, dP, accepted, last_z);
		show_dE = 1;
	      }
	      else
		drift_beam(coord, nToTrack, ((TW_LINAC*)eptr->p_elem)->length, run->default_order);
	      break;
	    case T_MAPSOLENOID:
	      nLeft = motion(coord, nToTrack, eptr->p_elem, eptr->type, P_central, 
			     &dgamma, dP, accepted, last_z);
	      break;
	    case T_TWLA:
	    case T_TWMTA:
	      nLeft = motion(coord, nToTrack, eptr->p_elem, eptr->type, P_central, 
			     &dgamma, dP, accepted, last_z);
	      show_dE = 1;
	      break;
	    case T_RCOL:
	      if (flags&TEST_PARTICLES && !(flags&TEST_PARTICLE_LOSSES))
		drift_beam(coord, nToTrack, ((RCOL*)eptr->p_elem)->length, run->default_order);
	      else {
		nLeft = rectangular_collimator(coord, (RCOL*)eptr->p_elem, nToTrack, accepted, last_z, *P_central);
	      }
	      break;
	    case T_ECOL:
	      if (flags&TEST_PARTICLES && !(flags&TEST_PARTICLE_LOSSES))
		drift_beam(coord, nToTrack, ((RCOL*)eptr->p_elem)->length, run->default_order);
	      else
		nLeft = elliptical_collimator(coord, (ECOL*)eptr->p_elem, nToTrack, accepted, z, *P_central);
	      break;
	    case T_CLEAN:
	      if (!(flags&TEST_PARTICLES && !(flags&TEST_PARTICLE_LOSSES)))
		nLeft = remove_outlier_particles(coord, (CLEAN*)eptr->p_elem, 
						 nToTrack, accepted, z, *P_central);
	      break;
	    case T_SCRAPER:
	      if (!(flags&TEST_PARTICLES && !(flags&TEST_PARTICLE_LOSSES))) {
		nLeft = beam_scraper(coord, (SCRAPER*)eptr->p_elem, nToTrack, accepted, z, *P_central);
	      }
	      break;
	    case T_PFILTER:
	      if (!(flags&TEST_PARTICLES && !(flags&TEST_PARTICLE_LOSSES)))
		nLeft = track_through_pfilter(coord, (PFILTER*)eptr->p_elem, nToTrack, 
					      accepted, z, *P_central);
	      break;
	    case T_CENTER:
	      center_beam(coord, (CENTER*)eptr->p_elem, nToTrack, i_pass);
	      break;
	    case T_REMCOR:
	      remove_correlations(coord, (REMCOR*)eptr->p_elem, nToTrack);
	      break;
	    case T_RFCA:
	      nLeft = simple_rf_cavity(coord, nToTrack, (RFCA*)eptr->p_elem, accepted, P_central, z);
	      break;
	    case T_RFCW:
	      nLeft = track_through_rfcw(coord, nToTrack, (RFCW*)eptr->p_elem, accepted, P_central, z,
					 run, i_pass, charge);
	      break;
	    case T_MODRF:
	      modulated_rf_cavity(coord, nToTrack, (MODRF*)eptr->p_elem, *P_central, z);
	      break;
	    case T_WATCH:
	      if (!(flags&TEST_PARTICLES) && !(flags&INHIBIT_FILE_OUTPUT)) {
	        watch = (WATCH*)eptr->p_elem;
	        if (!watch->disable) {
	          watch_pt_seen = 1;
	          if (!watch->initialized) 
	            set_up_watch_point(watch, run);
	          if (i_pass==0 && (n_passes/watch->interval)==0)
	            fprintf(stdout, "warning: n_passes = %ld and WATCH interval = %ld--no output will be generated!\n",
	     	     n_passes, watch->interval);
		  fflush(stdout);
		  if (i_pass>=watch->start_pass && (i_pass-watch->start_pass)%watch->interval==0) {
	            switch (watch->mode_code) {
	            case WATCH_COORDINATES:
		      dump_watch_particles(watch, step, i_pass, coord, nToTrack, *P_central,
			        	   beamline->revolution_length, 
					   charge?charge->macroParticleCharge*nToTrack:0.0, z);
		      break;
		    case WATCH_PARAMETERS:
		    case WATCH_CENTROIDS:
		      dump_watch_parameters(watch, step, i_pass, n_passes, coord, nToTrack, nOriginal, *P_central,
					    beamline->revolution_length);
		      break;
		    case WATCH_FFT:
		      dump_watch_FFT(watch, step, i_pass, n_passes, coord, nToTrack, nOriginal, *P_central);
		      break;
		    }
		  }
		}
	      }
	      break;
	    case T_HISTOGRAM:
	      if (!(flags&TEST_PARTICLES) && !(flags&INHIBIT_FILE_OUTPUT)) {
		histogram = (HISTOGRAM*)eptr->p_elem;
		if (!histogram->disable) {
		  watch_pt_seen = 1;   /* yes, this should be here */
		  if (!histogram->initialized) 
		    set_up_histogram(histogram, run);
		  if (i_pass==0 && (n_passes/histogram->interval)==0)
		    fprintf(stdout, "warning: n_passes = %ld and HISTOGRAM interval = %ld--no output will be generated!\n",
			    n_passes, histogram->interval);
		  fflush(stdout);
		  if (i_pass>=histogram->startPass && (i_pass-histogram->startPass)%histogram->interval==0) {
		    dump_particle_histogram(histogram, step, i_pass, coord, nToTrack, *P_central,
					    beamline->revolution_length, 
					    charge?charge->macroParticleCharge*nToTrack:0.0, z);
		  }
		}
	      }
	      break;
	    case T_MALIGN:
	      malign = (MALIGN*)eptr->p_elem;
	      if (malign->on_pass==-1 || malign->on_pass==i_pass)
		offset_beam(coord, nToTrack, (MALIGN*)eptr->p_elem, *P_central);
	      break;
	    case T_PEPPOT:
	      nLeft = pepper_pot_plate(coord, (PEPPOT*)eptr->p_elem, nToTrack, accepted);
	      break;
	    case T_ENERGY:
	      energy = (ENERGY*)eptr->p_elem;
	      if (energy->match_beamline) {
		if ((flags&FIDUCIAL_BEAM_SEEN) && eptr->Pref_output_fiducial>0)
		  /* Beamline momentum is defined.  Change particle reference momentum to match. */
		  set_central_momentum(coord, nToTrack, eptr->Pref_output_fiducial, P_central);
		else
		  /* Compute new central momentum to match the average momentum of the particles. */
		  do_match_energy(coord, nToTrack, P_central, 0);
		if (energy->match_particles)
		  bomb("can't match_beamline AND match_particles for ENERGY element", NULL);
	      }
	      else if (energy->match_particles) {
		/* change the particle momenta so that the centroid is the central momentum */
		do_match_energy(coord, nToTrack, P_central, 1);
	      }
	      else if (energy->central_energy)
		/* Change particle reference momentum to match the given energy */
		set_central_momentum(coord, nToTrack, sqrt(sqr(energy->central_energy+1)-1), 
				     P_central);
	      else if (energy->central_momentum)
		/* Change particle reference momentum to match the given value */
		set_central_momentum(coord, nToTrack, energy->central_momentum, P_central);
	      break;
	    case T_MAXAMP:
	      maxamp = (MAXAMP*) eptr->p_elem;
	      x_max = maxamp->x_max;
	      y_max = maxamp->y_max;
	      elliptical = maxamp->elliptical;
	      maxampOpenCode = determineOpenSideCode(maxamp->openSide);
	      maxampExponent = maxamp->exponent;
	      break;
	    case T_TRCOUNT:
	      /* >>>>> Needs to be updated */
	      /* 
	       *n_original = nLeft; 
	       if (accepted && i_pass==0)
               copy_particles(accepted, coord, *n_original);
              */
	      break;
	    case T_ALPH:
	      if (!eptr->matrix && !(eptr->matrix=compute_matrix(eptr, run, NULL)))
		bomb("no matrix for alpha magnet", NULL);
	      nLeft = alpha_magnet_tracking(coord, eptr->matrix, (ALPH*)eptr->p_elem, nToTrack,
					    accepted, *P_central, z);
	      break;
	    case T_MATR:
	      if (!eptr->matrix)
		eptr->matrix = &(((MATR*)eptr->p_elem)->M);
	      matr_element_tracking(coord, eptr->matrix, (MATR*)eptr->p_elem, nToTrack,
				    z);
	      break;
	    case T_EMATRIX:
	      if (!eptr->matrix)
		eptr->matrix = compute_matrix(eptr, run, NULL);
	      ematrix_element_tracking(coord, eptr->matrix, (EMATRIX*)eptr->p_elem, nToTrack,
				       z);
	      break;
	    case T_MULT:
	      nLeft = multipole_tracking(coord, nToTrack, (MULT*)eptr->p_elem, 0.0,
					 *P_central, accepted, z);
	      break;
	    case T_FMULT:
	      nLeft = fmultipole_tracking(coord, nToTrack, (FMULT*)eptr->p_elem, 0.0,
					  *P_central, accepted, z);
	      break;
	    case T_TAYLORSERIES:
	      nLeft = taylorSeries_tracking(coord, nToTrack, (TAYLORSERIES*)eptr->p_elem, 0.0,
					    *P_central, accepted, z);
	      break;
	    case T_KICKER:
	      if (flags&TIME_DEPENDENCE_OFF)
		drift_beam(coord, nToTrack, ((KICKER*)eptr->p_elem)->length, run->default_order);
	      else
		track_through_kicker(coord, nToTrack, (KICKER*)eptr->p_elem, *P_central, i_pass, run->default_order);
	      break;
	    case T_KSBEND:
	      nLeft = track_through_kick_sbend(coord, nToTrack, (KSBEND*)eptr->p_elem, 0.0,
					       *P_central, accepted, z);
	      break;
	    case T_CSBEND:
	      if (flags&TEST_PARTICLES) {
		saveISR = ((CSBEND*)eptr->p_elem)->isr;
		((CSBEND*)eptr->p_elem)->isr = 0;
	      }
	      nLeft = track_through_csbend(coord, nToTrack, (CSBEND*)eptr->p_elem, 0.0,
					   *P_central, accepted, last_z);
	      if (flags&TEST_PARTICLES)
		((CSBEND*)eptr->p_elem)->isr = saveISR;	  
	      break;
	    case T_CSRCSBEND:
	      if (flags&TEST_PARTICLES) {
		saveISR = ((CSRCSBEND*)eptr->p_elem)->isr;
		((CSRCSBEND*)eptr->p_elem)->isr = 0;
	      }
	      nLeft = track_through_csbendCSR(coord, nToTrack, (CSRCSBEND*)eptr->p_elem, 0.0,
					      *P_central, accepted, last_z, z, charge, run->rootname);
	      if (flags&TEST_PARTICLES)
		((CSRCSBEND*)eptr->p_elem)->isr = saveISR;
	      break;
	    case T_CSRDRIFT:
	      nLeft = track_through_driftCSR(coord, nToTrack, (CSRDRIFT*)eptr->p_elem,
					     *P_central, accepted, last_z, 
					     beamline->revolution_length,
					     run->rootname);
	      break;
	    case T_LSCDRIFT:
	      track_through_lscdrift(coord, nToTrack, (LSCDRIFT*)eptr->p_elem, *P_central, charge);
	      break;
	    case T_SCMULT:
	      trackThroughSCMULT(coord, nToTrack, eptr);
	      break;
	    case T_EDRIFT:
	      exactDrift(coord, nToTrack, ((EDRIFT*)eptr->p_elem)->length);
	      break;
	    case T_TUBEND:
	      nLeft = track_through_tubend(coord, nToTrack, 
					   (TUBEND*)eptr->p_elem, 0.0,
					   *P_central, accepted, z);
	      break;
	    case T_KQUAD:
	    case T_KSEXT:
	      nLeft = multipole_tracking2(coord, nToTrack, eptr, 0.0,
					    *P_central, accepted, z);
	      break;
	    case T_SAMPLE:
	      if (!(flags&TEST_PARTICLES))
		nLeft = sample_particles(coord, (SAMPLE*)eptr->p_elem, nToTrack, accepted, z, *P_central);
	      break;
	    case T_SCATTER:
	      if (!(flags&TEST_PARTICLES))
		scatter(coord, nToTrack, *P_central, (SCATTER*)eptr->p_elem);
	      break;
	    case T_DSCATTER:
	      if (!(flags&TEST_PARTICLES))
		distributionScatter(coord, nToTrack, *P_central, (DSCATTER*)eptr->p_elem, i_pass);
	      break;
	    case T_NIBEND:
	      nLeft = lorentz(coord, nToTrack, (NIBEND*)eptr->p_elem, T_NIBEND, *P_central, accepted);
	      break;
	    case T_NISEPT:
	      nLeft = lorentz(coord, nToTrack, (NISEPT*)eptr->p_elem, T_NISEPT, *P_central, accepted);
	      break;
	    case T_BMAPXY:
	      nLeft = lorentz(coord, nToTrack, (BMAPXY*)eptr->p_elem, T_BMAPXY, *P_central, accepted);
	      break;
	    case T_KPOLY:
	      nLeft = polynomial_kicks(coord, nToTrack, (KPOLY*)eptr->p_elem, 0.0,
				       *P_central, accepted, z);
	      break;
	    case T_RAMPRF:
	      ramped_rf_cavity(coord, nToTrack, (RAMPRF*)eptr->p_elem, *P_central, beamline->revolution_length,
			       eptr->end_pos, i_pass);
	      break;
	    case T_RAMPP:
	      ramp_momentum(coord, nToTrack, (RAMPP*)eptr->p_elem, P_central, i_pass);
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
	      track_particles(coord, eptr->matrix, coord, nToTrack);
	      break;
	    case T_MATTER:
	      track_through_matter(coord, nToTrack, (MATTER*)eptr->p_elem, *P_central);
	      break;
	    case T_RFMODE:
	      rfmode = (RFMODE*)eptr->p_elem;
	      if (!rfmode->initialized)
		set_up_rfmode(rfmode, eptr->name, eptr->end_pos, n_passes, run, 
			      nOriginal, *P_central,
			      beamline->revolution_length);
	      track_through_rfmode(coord, nToTrack, (RFMODE*)eptr->p_elem, *P_central,
				   eptr->name, eptr->end_pos, i_pass, n_passes,
				   charge);
	      break;
	    case T_FRFMODE:
	      frfmode = (FRFMODE*)eptr->p_elem;
	      if (!frfmode->initialized)
		set_up_frfmode(frfmode, eptr->name, eptr->end_pos, n_passes, run, 
			       nOriginal, *P_central,
			       beamline->revolution_length);
	      track_through_frfmode(coord, nToTrack, frfmode, *P_central,
				    eptr->name, eptr->end_pos, i_pass, n_passes,
				    charge);
	      break;
	    case T_TRFMODE:
	      trfmode = (TRFMODE*)eptr->p_elem;
	      if (!trfmode->initialized)
		set_up_trfmode(trfmode, eptr->name, eptr->end_pos, n_passes, run, nOriginal);
	      track_through_trfmode(coord, nToTrack, (TRFMODE*)eptr->p_elem, *P_central,
				    eptr->name, eptr->end_pos, i_pass, n_passes,
				    charge);
	      break;
	    case T_FTRFMODE:
	      ftrfmode = (FTRFMODE*)eptr->p_elem;
	      if (!ftrfmode->initialized)
		set_up_ftrfmode(ftrfmode, eptr->name, eptr->end_pos, n_passes, run, 
				nOriginal, *P_central,
				beamline->revolution_length);
	      track_through_ftrfmode(coord, nToTrack, ftrfmode, *P_central,
				     eptr->name, eptr->end_pos, i_pass, n_passes,
				     charge);
	      break;
	    case T_ZLONGIT:
	      track_through_zlongit(coord, nToTrack, (ZLONGIT*)eptr->p_elem, *P_central, run, i_pass,
				    charge);
	      break;
	    case T_ZTRANSVERSE:
	      track_through_ztransverse(coord, nToTrack, (ZTRANSVERSE*)eptr->p_elem, *P_central, run, i_pass,
					charge);
	      break;
	    case T_WAKE:
	      track_through_wake(coord, nToTrack, (WAKE*)eptr->p_elem, P_central, run, i_pass,
				 charge);
	      break;
	    case T_TRWAKE:
	      track_through_trwake(coord, nToTrack, (TRWAKE*)eptr->p_elem, *P_central, run, i_pass, 
				   charge);
	      break;
	    case T_SREFFECTS:
	      if (!(flags&TEST_PARTICLES))
		track_SReffects(coord, nToTrack, (SREFFECTS*)eptr->p_elem, *P_central, eptr->twiss, &(beamline->radIntegrals));
	      break;
	    case T_IBSCATTER:
	      if (!(flags&TEST_PARTICLES))
		track_IBS(coord, nToTrack, (IBSCATTER*)eptr->p_elem,
			  *P_central, &(beamline->elem), &(beamline->radIntegrals),
			  charge);
	      break;
	    case T_SCRIPT:
	      nLeft = transformBeamWithScript((SCRIPT*)eptr->p_elem, *P_central, charge, 
					      beam, coord, nToTrack, nLost, 
					      run->rootname, i_pass, run->default_order);
	      /* 
		 fprintf(stderr, "nLost=%ld, beam->n_particle=%ld, nLeft=%ld\n",
		 nLost, beam->n_particle, nLeft);
	      */
	      if (beam && coord!=beam->particle) {
		/* particles were created and so the particle array was changed */
		coord = beam->particle;
		if (nLost != (beam->n_particle - nLeft)) {
		  fprintf(stderr, "Particle accounting problem after return from script.\n");
		  fprintf(stderr, "nLost=%ld, beam->n_particle=%ld, nLeft=%ld\n",
			  nLost, beam->n_particle, nLeft);
		  exit(1);
		}
	      }

	      nToTrack = nLeft;
	      lostOnPass = beam->lostOnPass;
	      nMaximum = beam->n_particle;
	      break;
	    case T_FLOORELEMENT:
	      break;
	    case T_TFBPICKUP:
	      if (!(flags&TEST_PARTICLES))
		transverseFeedbackPickup((TFBPICKUP*)eptr->p_elem, coord, nToTrack, i_pass);
	      break;
	    case T_STRAY:
	      if (eptr->matrix)
		free_matrices(eptr->matrix);
	      stray = (STRAY*)eptr->p_elem;
	      eptr->matrix = stray_field_matrix(stray->length, &stray->lBx, &stray->gBx, 
						eptr->end_theta, stray->order?stray->order:run->default_order,
						*P_central, 
						stray->Wi);
	      track_particles(coord, eptr->matrix, coord, nToTrack);
	      break;
	    case T_TFBDRIVER:
	      if (!(flags&TEST_PARTICLES))
		transverseFeedbackDriver((TFBDRIVER*)eptr->p_elem, coord, nToTrack, beamline, i_pass, n_passes, run->rootname);
	      feedbackDriverSeen = 1;
	      break;
	    case T_LSRMDLTR:
	      nLeft = motion(coord, nToTrack, eptr->p_elem, eptr->type, P_central, 
			     &dgamma, dP, accepted, last_z);
	      show_dE = 1;
	      break;
	    case T_CWIGGLER:
	      GWigSymplecticPass(coord, nToTrack, *P_central, (CWIGGLER*)eptr->p_elem);
	      break;
	    default:
	      fprintf(stdout, "programming error: no tracking statements for element %s (type %s)\n",
		      eptr->name, entity_name[eptr->type]);
	      fflush(stdout);
	      exit(1);
	      break;
	    }
	  }
	}
        if (!(flags&TEST_PARTICLES && !(flags&TEST_PARTICLE_LOSSES)) && (x_max || y_max)) {
          if (!elliptical) 
            nLeft = limit_amplitudes(coord, x_max, y_max, nLeft, accepted, z, *P_central, 
				     eptr->type==T_DRIF || eptr->type==T_STRAY,
				     maxampOpenCode);
          else
            nLeft = elimit_amplitudes(coord, x_max, y_max, nLeft, accepted, z, *P_central, 
				      eptr->type==T_DRIF || eptr->type==T_STRAY,
				      maxampOpenCode, maxampExponent);
        }
        if (run->print_statistics && !(flags&TEST_PARTICLES)) {
          report_stats(stdout, ": ");
          /*
            if (show_dE && nLeft) {
            fprintf(stdout, "average energy imparted: %e MeV\n",
            dgamma*me_mev/nLeft);
            fflush(stdout);
            fprintf(stdout, "average x,y,z momentum imparted: %e, %e, %e MeV/c\n",
            dP[0]*me_mev/nLeft, dP[1]*me_mev/nLeft,
            dP[2]*me_mev/nLeft);
            fflush(stdout);
            }
	  */
          fprintf(stdout, "central momentum is %e    zstart = %em  zend = %em\n", *P_central, last_z, z);
          fflush(stdout);
          if (nLeft!=nToTrack)
            fprintf(stdout, "%ld particles left\n", nLeft);
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
          do_match_energy(coord, nLeft, P_central, 0);
        }
        eptr->Pref_output_fiducial = *P_central;
      } else if (flags&FIDUCIAL_BEAM_SEEN) {
        if (*P_central!=eptr->Pref_output_fiducial)
          set_central_momentum(coord, nLeft, eptr->Pref_output_fiducial, P_central);
      }
      else if (run->always_change_p0)
        do_match_energy(coord, nLeft, P_central, 0);
      if (i_pass==0 && traj_vs_z) {
        /* collect trajectory data--used mostly by trajectory correction routines */
        if (!traj_vs_z[i_traj].centroid) {
          fprintf(stdout, "error: the trajectory centroid array for %s is NULL (do_tracking)",
                  eptr->name);
          fflush(stdout);
          exit(1);
        }
        traj_vs_z[i_traj].elem = eptr;
        if (!(traj_vs_z[i_traj].n_part=nLeft)) {
          for (i=0; i<6; i++)
            traj_vs_z[i_traj].centroid[i] = 0;
        }
        else {
          for (i=0; i<6; i++) {
            for (j=sum=0; j<nToTrack; j++)
              sum += coord[j][i];
            traj_vs_z[i_traj].centroid[i] = sum/nLeft;
          }
        }
        i_traj++;
      }
      if (!(flags&TEST_PARTICLES) && sliceAnalysis && !sliceAnalysis->finalValuesOnly) {
	performSliceAnalysisOutput(sliceAnalysis, coord, nToTrack, 
				   !sliceAnDone, step, 
				   *P_central, 
				   charge?charge->macroParticleCharge*nToTrack:0.0, 
				   eptr->name, eptr->end_pos, 0); 
	sliceAnDone = 1;
      }
#if USE_MPI
      if (!(classFlags&UNIPROCESSOR)) {
        end_wtime = MPI_Wtime();
        my_wtime = my_wtime+end_wtime-start_wtime; 
      }
      else if (!(classFlags&(UNIDIAGNOSTIC&(~UNIPROCESSOR)))) { 
        /* a non-diagnostic uniprocessor element */
        if ((myid == 0) && (nMaximum!=(nLeft+nLost)))         
          /* there are additional losses occurred */
          lostSinceSeqMode = 1;
        MPI_Bcast (&lostSinceSeqMode, 1, MPI_INT, 0, MPI_COMM_WORLD);
      }
#endif
      if ((!USE_MPI) || (USE_MPI && active)) {
        nLeft = limit_amplitudes(coord, DBL_MAX, DBL_MAX, nLeft, accepted, z, *P_central, 0,
			       0);
        recordLossPass(lostOnPass, &nLost, nLeft, nMaximum, i_pass, myid, lostSinceSeqMode);
      }  

      if (getSCMULTSpecCount() && entity_description[eptr->type].flags&HAS_LENGTH) {
	/* calcaulate beam size at exit of element for use in space  charge calculation with SCMULT */
	/* need special care for element with 0 length but phase space rotation */
      	if (((DRIFT*)eptr->p_elem)->length > 0.0) 
        	accumulateSCMULT(coord, nToTrack, eptr);
      }

      last_type = eptr->type;
      eptrPred = eptr;
      eptr = eptr->succ;
      nToTrack = nLeft;
    } /* end of the while loop */
    
    if (!(flags&TEST_PARTICLES) && sliceAnalysis && sliceAnalysis->finalValuesOnly) {
      performSliceAnalysisOutput(sliceAnalysis, coord, nToTrack, 
				 !sliceAnDone, step, 
				 *P_central, 
				 charge?charge->macroParticleCharge*nToTrack:0.0, 
				 eptrPred->name, eptrPred->end_pos, 0);
      
      sliceAnDone = 1;
    }

    
    log_entry("do_tracking.2.2.3");
    if (effort)
      *effort += nLeft;

    if (sums_vs_z && (*sums_vs_z) && !(flags&FINAL_SUMS_ONLY) && !(flags&TEST_PARTICLES) &&
        (run->wrap_around || i_pass==n_passes-1)) {
      if (i_sums<0)
        bomb("attempt to accumulate beam sums with negative index!", NULL);
      accumulate_beam_sums(*sums_vs_z+i_sums, coord, nToTrack, *P_central);
      (*sums_vs_z)[i_sums].z = z;
#if defined(BEAM_SUMS_DEBUG)
      fprintf(stdout, "beam sums accumulated in slot %ld for %s at z=%em, sx=%e\n", 
              i_sums, name, z, sqrt((*sums_vs_z)[i_sums].sum2[0]/nLeft));
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
    
    if (!USE_MPI && (i_pass==0 || watch_pt_seen || feedbackDriverSeen)) {
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
          accumulate_beam_sums(*sums_vs_z+i_sums, coord, nToTrack, *P_central);
          (*sums_vs_z)[i_sums].z = z;
          i_sums++;
        }
        switch (eptr->type) {
        case T_TFBDRIVER:
          flushTransverseFeedbackDriverFiles((TFBDRIVER *)(eptr->p_elem));
          break;
        case T_WATCH:
          if (!(flags&TEST_PARTICLES) && !(flags&INHIBIT_FILE_OUTPUT)) {
            watch = (WATCH*)eptr->p_elem;
	    if (!watch->disable) {
	      if (i_pass%watch->interval==0) {
		switch (watch->mode_code) {
		case WATCH_COORDINATES:
		  break;
		case WATCH_PARAMETERS:
		case WATCH_CENTROIDS:
		  dump_watch_parameters(watch, step, i_pass, n_passes, coord, nToTrack, nOriginal, *P_central,
					beamline->revolution_length);
		  break;
		case WATCH_FFT:
		  dump_watch_FFT(watch, step, i_pass, n_passes, coord, nToTrack, nOriginal, *P_central);
		  break;
		}
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
 
#if USE_MPI
    if (balanceStatus==startMode) { 
      balanceStatus = checkBalance (my_wtime, myid, n_processors);  
      /* calculate the rate for all of the slave processors */
      if (myid==0) {
        my_rate = 0.0;
        nParPerElements = 0.0;
      }
      else {
        nParPerElements = (double)nParElements/(double)nElements; 
        my_rate = nParPerElements/my_wtime;
      } 
      lostSinceSeqMode = 1; /* set flag to distribute jobs according to  the speed */
    }
    else { /* The workload balancing will be checked for every pass by default.
              If user defined the CHECKFLAGS, the balance will be checked only 
              when the nToTrack is changed. */
  #ifdef CHECKFLAGS 
      if (myid==0) {
        if (old_nToTrack!=nToTrack) {
          checkFlags = 1;
	}
        else
          checkFlags = 0;
      }
      MPI_Bcast(&checkFlags, 1, MPI_INT, 0, MPI_COMM_WORLD);
  #else
      checkFlags = 1; /* the default option */
  #endif
      if (checkFlags)
        balanceStatus = checkBalance (my_wtime, myid, n_processors);   
      if (balanceStatus == badBalance) {
	if (myid==0) {
	  my_rate = 0.0;
	  nParPerElements = 0.0;
	}
	else {
	  nParPerElements = (double)nParElements/(double)nElements; 
	  my_rate = nParPerElements/my_wtime;
	} 
      }
    }
    if (myid==0)
      old_nToTrack = nToTrack;
#endif
  } /* end of the for loop for n_passes*/

#if USE_MPI
  /* change back to sequential mode before leaving the do_tracking function */
  if (parallelStatus!=notParallel) {
    gatherParticles(coord, lostOnPass, &nToTrack, &nLost, accepted, n_processors, myid, &round);
    parallelStatus = notParallel ;
  }
#endif

  /* do this here to get report of CSR drift normalization */
  reset_driftCSR();

  log_exit("do_tracking.2");
  log_entry("do_tracking.3");
  
  if (nLeft && sums_vs_z && *sums_vs_z && !(flags&TEST_PARTICLES)) {
    if (flags&FINAL_SUMS_ONLY) {
      log_entry("do_tracking.3.1");
      i_sums = 0;
      accumulate_beam_sums(*sums_vs_z+i_sums, coord, nToTrack, *P_central);
      (*sums_vs_z)[i_sums].z = z;
#if defined(BEAM_SUMS_DEBUG)
      fprintf(stdout, "beam sums accumulated in slot %ld for final sums at z=%em, sx=%e\n", 
              i_sums, z, sqrt((*sums_vs_z)[i_sums].sum2[0]/nLeft));
      fflush(stdout);
#endif
      log_exit("do_tracking.3.1");
    }
    else if (run->wrap_around) {
      log_entry("do_tracking.3.2");
      if (i_sums<0)
        bomb("attempt to accumulate beam sums with negative index!", NULL);
      /* accumulate sums for final output */
      accumulate_beam_sums(*sums_vs_z+i_sums, coord, nToTrack, *P_central);
#if defined(BEAM_SUMS_DEBUG)
      fprintf(stdout, "beam sums accumulated in slot %ld for final sums at z=%em, sx=%e\n", 
              i_sums, z, sqrt((*sums_vs_z)[i_sums].sum2[0]/nLeft));
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
    computeSASEFELAtEnd(sasefel, coord, nToTrack, *P_central, charge->macroParticleCharge*nToTrack);
  }
  
  log_exit("do_tracking.3");
  log_entry("do_tracking.4");
  if (!(flags&SILENT_RUNNING) && !is_batch && n_passes!=1 && !(flags&TEST_PARTICLES)) {
    fprintf(stdout, "%ld particles present after pass %ld        \n", 
            nToTrack, i_pass);
    fflush(stdout);
  }

#ifdef MPI_DEBUG  
  #ifdef CHECKFLAGS 
    printf("Balance is checked for the first pass and when particles are lost only.\n"); 
    fflush(stdout);
  #else
    printf("Balance is checked for every pass.\n"); 
    fflush(stdout);
  #endif
#endif

  log_exit("do_tracking.4");

  log_exit("do_tracking");
  if (charge && finalCharge)
    *finalCharge = nToTrack*charge->macroParticleCharge;
  return(nToTrack);
}

void offset_beam(
                 double **coord,
                 long nToTrack, 
                 MALIGN *offset,
                 double P_central
                 )
{
  long i_part;
  double *part, pc, beta, gamma, t;
  double ds;
  
  log_entry("offset_beam");
  
  for (i_part=nToTrack-1; i_part>=0; i_part--) {
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


void center_beam(double **part, CENTER *center, long np, long iPass)
{
  double sum, offset;
  long i, ic;
  long centerCoord[4];

  if (!np) {
    return;
  }

  if (center->onPass>=0 && iPass!=center->onPass)
    return;

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
  for (ip=0; ip<np; ip++) {
    if (scat->probability<1 && random_2(1)>scat->probability)
      continue;
    for (i=0; i<4; i++) {
      if (!sigma[i])
        continue;
      part[ip][i] += gauss_rn(0, random_2)*sigma[i];
    }
    if (scat->dp) {
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
        fpt->matrix_mem[count++] = rpn_create_mem(buffer, 0);
      }
    }
    if (M->order>1) {
      for (i=0; i<6; i++) {
        for (j=0; j<6; j++) {
          for (k=0; k<=j; k++) {
            sprintf(buffer, "%s#%ld.T%ld%ld%ld", name, occurence, i+1, j+1, k+1);
            fpt->matrix_mem[count++] = rpn_create_mem(buffer, 0);
          }
        }
      }
    }
    fpt->init_flags |= 8;
  }
  
  for (i=count=0; i<6; i++)
    for (j=0; j<6; j++)
      rpn_store(M->R[i][j], NULL, fpt->matrix_mem[count++]);
  if (M->order>1)
    for (i=0; i<6; i++)
      for (j=0; j<6; j++)
        for (k=0; k<=j; k++) 
          rpn_store(M->T[i][j][k], NULL, fpt->matrix_mem[count++]);
}


void store_fitpoint_beam_parameters(MARK *fpt, char *name, long occurence, double **coord, long np, double Po)
{
  long i, j, k;
  static double emit[3], sigma[6], centroid[6];
  static BEAM_SUMS sums;
  static char *centroid_name_suffix[8] = {
    "Cx", "Cxp", "Cy", "Cyp", "Cs", "Cdelta", "pCentral", "Particles" };
  static char *sigma_name_suffix[6] = {
    "Sx", "Sxp", "Sy", "Syp", "Ss", "Sdelta" };
  static char *emit_name_suffix[3] = {
    "ex", "ey", "es"};
  static char s[1000];

  zero_beam_sums(&sums, 1);
  accumulate_beam_sums(&sums, coord, np, Po);
  for (i=0; i<6; i++) {
    centroid[i] = sums.centroid[i];
    sigma[i] = sqrt(sums.sigma[i][i]);
    if ((emit[i/2] = sums.sigma[i][i]*sums.sigma[i+1][i+1] - sqr(sums.sigma[i][i+1]))>0)
      emit[i/2] = sqrt(emit[i/2]);
    else
      emit[i/2] = 0;
  }
  
  if (!(fpt->init_flags&2)) {
    fpt->centroid_mem = tmalloc(sizeof(*fpt->centroid_mem)*8);
    fpt->sigma_mem = tmalloc(sizeof(*fpt->sigma_mem)*6);
    fpt->emit_mem = tmalloc(sizeof(*fpt->emit_mem)*3);
    fpt->sij_mem = tmalloc(sizeof(*fpt->sigma_mem)*15);
    for (i=0; i<8; i++) {
      sprintf(s, "%s#%ld.%s", name, occurence, centroid_name_suffix[i]);
      fpt->centroid_mem[i] = rpn_create_mem(s, 0);
    }
    for (i=0; i<6; i++) {
      sprintf(s, "%s#%ld.%s", name, occurence, sigma_name_suffix[i]);
      fpt->sigma_mem[i] = rpn_create_mem(s, 0);
    }
    for (i=0; i<3; i++) {
      sprintf(s, "%s#%ld.%s", name, occurence, emit_name_suffix[i]);
      fpt->emit_mem[i] = rpn_create_mem(s, 0);
    }
    for (i=k=0; i<6; i++) {
      for (j=i+1; j<6; j++, k++) {
        sprintf(s, "%s#%ld.s%ld%ld", name, occurence, i+1, j+1);
        fpt->sij_mem[k] = rpn_create_mem(s, 0);
      }
    }
    fpt->init_flags |= 2;
  }
  for (i=0; i<6; i++) {
    rpn_store(centroid[i], NULL, fpt->centroid_mem[i]);
    rpn_store(sigma[i], NULL, fpt->sigma_mem[i]);
  }
  for (i=0; i<3; i++)
    rpn_store(emit[i], NULL, fpt->emit_mem[i]);
  for (i=k=0; i<6; i++)
    for (j=i+1; j<6; j++, k++)
      rpn_store(sums.sigma[i][j], NULL, fpt->sij_mem[k]);
  rpn_store(Po, NULL, fpt->centroid_mem[6]);
  rpn_store((double)np, NULL, fpt->centroid_mem[7]);
}

ELEMENT_LIST *findBeamlineMatrixElement(ELEMENT_LIST *eptr)
{
  ELEMENT_LIST *eptr0=NULL, *eptrPassed;
  long matrixSeen = 0;
  eptrPassed = eptr;
  while (eptr) {
    if (eptr->type==T_MATR) {
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
    if ((eptr->p_elem || eptr->matrix) && eptr->type==T_MATR) {
      fprintf(stderr, "***** WARNING ****\n");
      fprintf(stderr, "Possible problem with \"linear chromatic\" or \"longitudinal-only\" matrix tracking\n");
      fprintf(stderr, "Concatenation resulted in more than one matrix.  Make sure the additional matrices do\n");
      fprintf(stderr, "not affect the revolution matrix!\n");
      print_elem_list(stderr, eptrPassed);
      fprintf(stderr, "***** WARNING ****\n");
      break;
    }
    eptr = eptr->succ;
  }
  return eptr0;
}

long trackWithChromaticLinearMatrix(double **particle, long particles, double **accepted,
                                    double Po, double z, ELEMENT_LIST *eptr,
                                    TWISS *twiss, double *tune0,
                                    double *chrom,    /* d   nu /ddelta   */
                                    double *chrom2,   /* d^2 nu /ddelta^2 */
                                    double *chrom3,   /* d^3 nu /ddelta^3 */
                                    double *dbeta_dPoP, 
                                    double *dalpha_dPoP,
                                    double *alphac,   /* Cs = Cs(0) + delta*alphac[0] + delta^2*alphac[1] */
                                    double *eta2,
                                    double *eta3      /* x = x(0) + eta*delta + eta2*delta^2 + eta3*delta^3 */
                                    )
{
  long ip, plane, offset, i, j, itop, is_lost;
  double *coord, deltaPoP, tune2pi, sin_phi, cos_phi;
  double alpha[2], beta[2], eta[4], beta1, alpha1;
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
  eta[0] = twiss->etax;
  eta[1] = twiss->etapx;
  eta[2] = twiss->etay;
  eta[3] = twiss->etapy;
  for (i=0; i<6; i++) {
    M1->C[i] = eptr->matrix->C[i];
    for (j=0; j<6; j++)
      M1->R[i][j] = i==j?1:0;
  }
  itop = particles-1;
  for (ip=0; ip<particles; ip++) {
    coord = particle[ip];
    deltaPoP = coord[5];
    is_lost = 0;
    if (deltaPoP!=lastDPoP) {
      for (plane=0; !is_lost && plane<2; plane++) {
        tune2pi = PIx2*(tune0[plane] + 
                        deltaPoP*(chrom[plane] +
                                  deltaPoP/2*(chrom2[plane] + 
                                              deltaPoP/3*chrom3[plane])));
        offset = 2*plane;
        if ((beta1 = beta[plane]+dbeta_dPoP[plane]*deltaPoP)<=0) {
          fprintf(stdout, "nonpositive beta function for particle with delta=%le\n",
                  deltaPoP);
          fprintf(stdout, "particle is lost\n");
          lastDPoP = DBL_MAX;
          is_lost = 1;
          continue;
        }
        if (fabs( ((long)(2*tune2pi/PIx2)) - ((long)(2*tune0[plane]))) != 0) {
          fprintf(stdout, "particle with delta=%le crossed integer or half-integer resonance\n",
                  deltaPoP);
          fprintf(stdout, "particle is lost\n");
          lastDPoP = DBL_MAX;
          is_lost = 1;
          continue;
        }
        /* R11=R22 or R33=R44 */
        sin_phi = sin(tune2pi);
        cos_phi = cos(tune2pi);
        alpha1 = alpha[plane]+dalpha_dPoP[plane]*deltaPoP;
        /* R11 or R33 */
        R11 = M1->R[0+offset][0+offset] = cos_phi + alpha1*sin_phi;
        /* R22 or R44 */
        R22 = M1->R[1+offset][1+offset] = cos_phi - alpha1*sin_phi;
        /* R12 or R34 */
        if ((R12 = M1->R[0+offset][1+offset] = beta1*sin_phi)) {
          /* R21 or R43 */
          M1->R[1+offset][0+offset] = (R11*R22-1)/R12;
        }
        else {
          bomb("divided by zero in trackWithChromaticLinearMatrix", NULL);
        }
        det = M1->R[0+offset][0+offset]*M1->R[1+offset][1+offset] -
          M1->R[0+offset][1+offset]*M1->R[1+offset][0+offset];
        if (fabs(det-1)>1e-6) {
          fprintf(stdout, "Determinant is suspect for particle with delta=%e\n", deltaPoP);
          fprintf(stdout, "particle is lost\n");
          lastDPoP = DBL_MAX;
          is_lost = 1;
          continue;
        }
      }
    }
    if (is_lost) {
      swapParticles(particle[ip], particle[itop]);
      if (accepted)
        swapParticles(accepted[ip], accepted[itop]);
      particle[itop][4] = z;
      particle[itop][5] = Po*(1+deltaPoP);
      --itop;
      --ip;
      --particles;
    } else {
      lastDPoP = deltaPoP;
      for (plane=0; !is_lost && plane<2; plane++) {
        /* remove the dispersive orbit from the particle coordinates */
        coord[2*plane]   -= deltaPoP*(eta[2*plane]   + deltaPoP*eta2[2*plane]);
        coord[2*plane+1] -= deltaPoP*(eta[2*plane+1] + deltaPoP*eta2[2*plane+1]);
      }
      /* momentum-dependent pathlength --- note that other path-length terms are ignored ! */
      M1->C[4] = eptr->matrix->C[4]*(1 + deltaPoP*(alphac[0] + deltaPoP*alphac[1]));
      coord[5] -= lastDPoP;
      track_particles(&coord, M1, &coord, 1);
      coord[5] += lastDPoP;
      for (plane=0; !is_lost && plane<2; plane++) {
        /* add back the dispersive orbit to the particle coordinates */
        coord[2*plane]   += deltaPoP*eta[2*plane]   + sqr(deltaPoP)*eta2[2*plane];
        coord[2*plane+1] += deltaPoP*eta[2*plane+1] + sqr(deltaPoP)*eta2[2*plane+1];
      }
    }
  }
  return particles;
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

void matr_element_tracking(double **coord, VMATRIX *M, MATR *matr,
                           long np, double z)
/* subtract off <s> prior to using a user-supplied matrix to avoid possible
 * problems with R5? and T?5? elements
 */
{
  long i;
  if (!np)
    return;
  if (!matr) {
    track_particles(coord, M, coord, np);
  } else {
    if (!matr->fiducialSeen) {
      double sum = 0;
      for (i=0; i<np; i++)
        sum += coord[i][4];
      matr->sReference = sum/np;
      matr->fiducialSeen = 1;
    }
    for (i=0; i<np; i++)
      coord[i][4] -= matr->sReference;
    track_particles(coord, M, coord, np);
    for (i=0; i<np; i++)
      coord[i][4] += matr->sReference;
  }
}

void ematrix_element_tracking(double **coord, VMATRIX *M, EMATRIX *matr,
			      long np, double z)
/* subtract off <s> prior to using a user-supplied matrix to avoid possible
 * problems with R5? and T?5? elements
 */
{
  long i;
  if (!np)
    return;
  if (!matr) {
    fprintf(stderr, "ematrix_element_tracking: matr=NULL, tracking with M (%ld order)\n",
            M->order);
    track_particles(coord, M, coord, np);
  } else {
    if (!matr->fiducialSeen) {
      double sum = 0;
      for (i=0; i<np; i++)
        sum += coord[i][4];
      matr->sReference = sum/np;
      matr->fiducialSeen = 1;
    }
    for (i=0; i<np; i++)
      coord[i][4] -= matr->sReference;
    track_particles(coord, M, coord, np);
    for (i=0; i<np; i++)
      coord[i][4] += matr->sReference;
  }
}

long transformBeamWithScript(SCRIPT *script, double pCentral, CHARGE *charge, 
                             BEAM *beam, double **part, long np, long nLost,
                             char *mainRootname, long iPass, long driftOrder)
{
  char *rootname, *input, *output;
  char *cmdBuffer0, *cmdBuffer1;
  SDDS_DATASET SDDSout, SDDSin;
  double *data;
  char *dataname[6] = {"x","xp","y","yp","t","p"};
  long i, j, npNew, nameLength, doDrift;

  if (!script->rootname || !strlen(script->rootname)) {
    /* generate random rootname */
    if (!(rootname = tmpname(NULL)))
      bomb("problem generating temporary filename for script", NULL);
  } else 
    rootname = compose_filename(script->rootname, mainRootname);
  nameLength = (script->directory?strlen(script->directory):0) + \
    strlen(rootname) + strlen(script->inputExtension) +
    strlen(script->outputExtension) + 4;
  if (!(input = malloc(sizeof(*input)*nameLength)) ||
      !(output = malloc(sizeof(*output)*nameLength)))
    bomb("problem generating temporary filename for script", NULL);

  doDrift = 0;
  if (script->onPass>=0) {
    if (script->onPass!=iPass)
      doDrift = 1;
  } else if (script->startPass>=0) {
    if (script->startPass>iPass)
      doDrift = 1;
  }

  if (doDrift) {
    drift_beam(part, np, script->length, driftOrder);
    return np;
  }

  /* prepare command */
  if (script->directory && strlen(script->directory)) {
#if defined(_WIN32)
    sprintf(input, "%s\\%s.%s", script->directory, rootname, script->inputExtension);
    sprintf(output, "%s\\%s.%s", script->directory, rootname, script->outputExtension);
#else
    sprintf(input, "%s/%s.%s", script->directory, rootname, script->inputExtension);
    sprintf(output, "%s/%s.%s", script->directory, rootname, script->outputExtension);
#endif
  } else {
    sprintf(input, "%s.%s", rootname, script->inputExtension);
    sprintf(output, "%s.%s", rootname, script->outputExtension);
  }
  if (rootname!=script->rootname)
    free(rootname);

  if (!(cmdBuffer0=malloc(sizeof(char)*(strlen(script->command)+10*strlen(input)+10*strlen(output)))) ||
      !(cmdBuffer1=malloc(sizeof(char)*(strlen(script->command)+10*strlen(input)+10*strlen(output)))))
    bomb("memory allocation failure making command buffer for script", NULL);
  replaceString(cmdBuffer0, script->command, "%i", input, 9, 0);
  replaceString(cmdBuffer1, cmdBuffer0, "%o", output, 9, 0);

  /* substitute numerical parameters */
  for (i=0; i<10; i++) {
    long count = 0;
    char tag[10], value[25], *ptr;
    sprintf(tag, "%%np%ld", i);
    ptr = cmdBuffer1;
    while ((ptr=strstr(ptr, tag))) {
      count ++;
      ptr += 3;
    }
    if (!count) continue;
    sprintf(value, "%21.15e", script->NP[i]);
    if (!(cmdBuffer0 = SDDS_Realloc(cmdBuffer0, sizeof(*cmdBuffer1)*(strlen(cmdBuffer1)+count*25+1))) ||
        !(cmdBuffer1 = SDDS_Realloc(cmdBuffer1, sizeof(*cmdBuffer1)*(strlen(cmdBuffer1)+count*25+1))))
      SDDS_Bomb("memory allocation failure");
    replaceString(cmdBuffer0, cmdBuffer1, tag, value, count, 0);
    strcpy(cmdBuffer1, cmdBuffer0);
  }
  /* substitute string parameters */
  for (i=0; i<10; i++) {
    long count = 0;
    char tag[10], *ptr;
    if (!script->SP[i] || strlen(script->SP[i])==0)
      continue;
    sprintf(tag, "%%sp%ld", i);
    ptr = cmdBuffer1;
    while ((ptr=strstr(ptr, tag))) {
      count ++;
      ptr += 3;
    }
    if (!count) continue;
    if (!(cmdBuffer0 = 
          SDDS_Realloc(cmdBuffer0, sizeof(*cmdBuffer1)*(strlen(cmdBuffer1)+count*strlen(script->SP[i])+1))) ||
        !(cmdBuffer1 = 
          SDDS_Realloc(cmdBuffer1, sizeof(*cmdBuffer1)*(strlen(cmdBuffer1)+count*strlen(script->SP[i])+1))))
      SDDS_Bomb("memory allocation failure");
    replaceString(cmdBuffer0, cmdBuffer1, tag, script->SP[i], count, 0);
    strcpy(cmdBuffer1, cmdBuffer0);
  }

  if (script->verbosity>0) {
    fprintf(stdout, "%s\n", cmdBuffer1);
    fflush(stdout);
  }
  
  /* dump the data to script input file */
  SDDS_ForceInactive(&SDDSout);
  SDDS_PhaseSpaceSetup(&SDDSout, input, SDDS_BINARY, 1, "script input", 
		       "unknown", "unknown",
                       "transformBeamWithScript");
  dump_phase_space(&SDDSout, part, np, 0, pCentral, charge?charge->macroParticleCharge*np:0.0);

  if (!SDDS_Terminate(&SDDSout))
    SDDS_Bomb("problem terminating script input file");
  
#if defined(CONDOR_COMPILE)
  _condor_ckpt_disable();
#endif

  /* run the script */
  if (script->useCsh)
    executeCshCommand(cmdBuffer1);
  else 
    system(cmdBuffer1);

#if defined(CONDOR_COMPILE)
  _condor_ckpt_enable();
#endif

  if (script->verbosity>0) {
    fprintf(stdout, "Command completed\n");
    fflush(stdout);
  }
  
  /* read the data from script output file */
  if (!fexists(output)) 
    SDDS_Bomb("unable to find script output file");
  if (!SDDS_InitializeInput(&SDDSin, output)) {
    SDDS_SetError("Unable to read script output file");
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
  }
  if (!check_sdds_column(&SDDSin, "x", "m") ||
      !check_sdds_column(&SDDSin, "y", "m") ||
      !check_sdds_column(&SDDSin, "xp", NULL) ||
      !check_sdds_column(&SDDSin, "yp", NULL) ||
      !check_sdds_column(&SDDSin, "p", "m$be$nc") ||
      !check_sdds_column(&SDDSin, "t", "s")) {
    if (!check_sdds_column(&SDDSin, "p", "m$be$nc") &&
        check_sdds_column(&SDDSin, "p", NULL)) {
      fprintf(stdout, "Warning: p has no units in script output file.  Expected m$be$nc\n");
      fflush(stdout);
    } else {
      fprintf(stdout, 
              "necessary data quantities (x, x', y, y', t, p) have the wrong units or are not present in script output");
      fflush(stdout);
      exit(1);
    }
  }
  
  if (SDDS_ReadPage(&SDDSin)!=1) {
    SDDS_SetError("Unable to read script output file");
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
  }
  npNew = SDDS_RowCount(&SDDSin);
  if (script->verbosity>0) {
    fprintf(stdout, "%ld particles in script output file\n", npNew);
    fflush(stdout);
  }
  if (!npNew) {
    return 0;
  }
  if (npNew>np) {
    /* We have to resize the arrays in the BEAM structure */
    /*
      fprintf(stdout, "Increasing number of particles from %ld (%ld active) to %ld (%ld active)\n",
      np+nLost, np, npNew+nLost, npNew);
    */
    if (!beam) {
      fprintf(stderr, "Error: script element increased the number of particles from %ld to %ld\n.",
              np, npNew);
      fprintf(stderr, "This happened (apparently) during a pre-tracking stage, which isn't allowed\n");
      fprintf(stderr, "in this version of elegant.\n");
      exit(1);
    }
    if ((np+nLost)!=beam->n_particle) {
      fprintf(stderr, "Particle accounting problem in SCRIPT element:\n");
      fprintf(stderr, "np = %ld, nLost = %ld, n_particle = %ld\n",
              np, nLost, beam->n_particle);
      exit(1);
    }
    
    if (beam->original==beam->particle) {
      /* This means, oddly enough, that the particle array and original array are the same because the
       * separate original array wasn't needed.  n_original gives the size of both arrays (including
       * live and lost particles).  To avoid confusion, we'll copy the data to a new array before
       * doing anything else, even though it means the original array is not used for anything and
       * contains a useless frozen copy of the present beam.
       * Use n_original since that's the size of the array, including lost particles. 
       */
      beam->particle = (double**)czarray_2d(sizeof(double), beam->n_original, 7);
      copy_particles(beam->particle, beam->original, beam->n_original);
    }
    /* resize the particle array, leaving space for the lost particle data at the top */
    if (!(beam->particle = realloc(beam->particle, sizeof(*(beam->particle))*(npNew+nLost))) ||
        !(beam->lostOnPass = realloc(beam->lostOnPass, sizeof(beam->lostOnPass)*(npNew+nLost)))) {
      fprintf(stderr, "Memory allocation failure increasing particle array size to %ld\n",
              npNew+nLost);
      exit(1);
    }
    /* allocate space for new particle pointers */
    for (i=np+nLost; i<npNew+nLost; i++) 
      beam->particle[i] = tmalloc(sizeof(**(beam->particle))*7);
    /* move lost particles into the upper part of the arrays */
    /* fprintf(stdout, "Moving %ld lost particles higher into buffer\n",
       nLost);
    */
    for (i=nLost-1; i>=0; i--) {
      swapParticles(beam->particle[np+i], beam->particle[npNew+i]);
      SWAP_LONG(beam->lostOnPass[np+i], beam->lostOnPass[npNew+i]);
    }
    if (beam->accepted)  {
      /* this data is invalid when particles are added */
      free_czarray_2d((void**)beam->accepted, np+nLost, 7);
      beam->accepted = NULL;
    }
    beam->n_particle = npNew+nLost;
    beam->n_to_track = npNew;
    /* fprintf(stdout, "beam->n_particle = %ld, beam->n_to_track = %ld\n",
       beam->n_particle, beam->n_to_track);
    */
    part = beam->particle;
  }
  for (i=0; i<6; i++) {
    if (!(data = SDDS_GetColumnInDoubles(&SDDSin, dataname[i]))) {
      SDDS_SetError("Unable to read script output file");
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
    }
    for (j=0; j<npNew; j++)
      part[j][i] = data[j];
    free(data);
  }
  /* assign new particle IDs ??? */
  for (j=0; j<npNew; j++)
    part[j][6] = j;

  if (charge) {
    double totalCharge;
    if (!SDDS_GetParameterAsDouble(&SDDSin, "Charge", &totalCharge)) {
      SDDS_SetError("Unable to read Charge parameter from script output file");
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
    }
    charge->charge = totalCharge;
    charge->macroParticleCharge = 0;
    if (npNew)
      charge->macroParticleCharge = totalCharge/npNew;
  }
  if (SDDS_ReadPage(&SDDSin)!=-1)
    SDDS_Bomb("Script output file has multiple pages");
  if (!SDDS_Terminate(&SDDSin))
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
  if (script->verbosity) {
    fprintf(stdout, "done with file\n");
    fflush(stdout);
  }
    
  /* convert (t, p) data to (s, delta) */
  for (j=0; j<npNew; j++) {
    double p, beta;
    p = part[j][5];
    part[j][5] = (p-pCentral)/pCentral;
    beta = p/sqrt(sqr(p)+1);
    part[j][4] *= beta*c_mks;
  }

  if (!script->keepFiles) {
    /* delete the input and output files */
    remove(input);
    remove(output);
  }

  /* clean up */
  free(cmdBuffer0);
  free(cmdBuffer1);
  
  return npNew;
}

void distributionScatter(double **part, long np, double Po, DSCATTER *scat, long iPass)
{
  static DSCATTER_GROUP *dscatterGroup = NULL;
  static long dscatterGroups = 0;
  long i, ip, interpCode, nScattered, nLeftThisPass;
  double t, P, beta, amplitude, cdf;
  TRACKING_CONTEXT context;
    
  if (!np)
    return;
  getTrackingContext(&context);

  if (!scat->initialized) {
    SDDS_DATASET SDDSin;
    static char *planeName[3] = {"xp", "yp", "dp"};
    static short planeIndex[3] = {1, 3, 5};
    scat->initialized = 1;
    scat->indepData = scat->cdfData = NULL;
    scat->groupIndex = -1;
    if ((i=match_string(scat->plane, planeName, 3, MATCH_WHOLE_STRING))<0) {
      fprintf(stderr, "Error for %s: plane is not valid.  Give xp, yp, or dp\n",
              context.elementName);
      exit(1);
    }
    scat->iPlane = planeIndex[i];
    if (!SDDS_InitializeInputFromSearchPath(&SDDSin, scat->fileName) ||
        SDDS_ReadPage(&SDDSin)!=1) {
      fprintf(stderr, "Error for %s: file is not valid.\n", context.elementName);
      exit(1);
    }
    if ((scat->nData=SDDS_RowCount(&SDDSin))<2) {
      fprintf(stderr, "Error for %s: file contains insufficient data.\n", context.elementName);
      exit(1);
    }
    /* Get independent data */
    if (!(scat->indepData=SDDS_GetColumnInDoubles(&SDDSin, scat->valueName))) {
      fprintf(stderr, "Error for %s: independent variable data is invalid.\n",
              context.elementName);
      exit(1);
    }
    /* Check that independent data is monotonically increasing */
    for (i=1; i<scat->nData; i++)
      if (scat->indepData[i]<=scat->indepData[i-1]) {
        fprintf(stderr, "Error for %s: independent variable data is not monotonically increasing.\n",
                context.elementName);
        exit(1);
      }
    /* Get CDF or PDF data */
    if (!(scat->cdfData=SDDS_GetColumnInDoubles(&SDDSin, 
                                                scat->cdfName?scat->cdfName:scat->pdfName))) {
      fprintf(stderr, "Error for %s: CDF/PDF data is invalid.\n",
              context.elementName);
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
      exit(1);
    }
    SDDS_Terminate(&SDDSin);
    if (!(scat->cdfName)) {
      /* must integrate to get the CDF */
      double *integral, *ptr;
      if (!(integral=malloc(sizeof(*integral)*scat->nData))) {
        fprintf(stderr, "Error for %s: memory allocation failure.\n", context.elementName);
        exit(1);
      }
      trapazoidIntegration1(scat->indepData, scat->cdfData, scat->nData, integral);
      ptr = scat->cdfData;
      scat->cdfData = integral;
      free(ptr);
    }
    /* Check that CDF is monotonically increasing */
    for (i=0; i<scat->nData; i++) {
      if (scat->cdfData[i]<0 || (i && scat->cdfData[i]<=scat->cdfData[i-1])) {
        long j;
        fprintf(stderr, "Error for %s: CDF not monotonically increasing at index %ld.\n", context.elementName,
                i);
        for (j=0; j<=i+1 && j<scat->nData; j++)
          fprintf(stderr, "%ld %21.15e\n", j, scat->cdfData[i]);
        exit(1);
      }
    }
    /* Normalize CDF to 1 */
    for (i=0; i<scat->nData; i++) 
      scat->cdfData[i] /= scat->cdfData[scat->nData-1];

    if (scat->oncePerParticle) {
      /* Figure out scattering group assignment */
      if (scat->group>=0) {
        for (i=0; i<dscatterGroups; i++) {
          if (dscatterGroup[i].group==scat->group)
            break;
        }
      } else
        i = dscatterGroups;
      if (i==dscatterGroups) {
        /* make a new group */
        fprintf(stderr, "Making new scattering group structure (#%ld, group %ld) for %s #%ld\n",
                i, scat->group, context.elementName, context.elementOccurrence);
        if (!(dscatterGroup = SDDS_Realloc(dscatterGroup, sizeof(*dscatterGroup)*(dscatterGroups+1)))) {
          fprintf(stderr, "Error for %s: memory allocation failure.\n", context.elementName);
          exit(1);
        }
        dscatterGroup[i].particleIDScattered = NULL;
        dscatterGroup[i].nParticles = 0;
        dscatterGroup[i].group = scat->group;
        scat->firstInGroup = 1;
        dscatterGroups++;
      } else {
        fprintf(stderr, "Linking to scattering group structure (#%ld, group %ld) for %s #%ld\n",
                i, scat->group, context.elementName, context.elementOccurrence);
        scat->firstInGroup = 0;
      }
      scat->groupIndex = i;
    }
    if (scat->startOnPass<0)
      scat->startOnPass = 0;
  }


  if (iPass==0 && scat->oncePerParticle && scat->firstInGroup) {
    /* Initialize data for remembering which particles have been scattered */
    dscatterGroup[scat->groupIndex].nParticles = np;
    if (dscatterGroup[scat->groupIndex].particleIDScattered)
      free(dscatterGroup[scat->groupIndex].particleIDScattered);
    if (!(dscatterGroup[scat->groupIndex].particleIDScattered 
          = malloc(sizeof(*(dscatterGroup[scat->groupIndex].particleIDScattered))*np))) {
      fprintf(stderr, "Error for %s: memory allocation failure.\n", context.elementName);
      exit(1);
    }
    dscatterGroup[scat->groupIndex].nScattered = 0;
    dscatterGroup[scat->groupIndex].allScattered = 0;
  }

  if (iPass<scat->startOnPass) 
    return;
  if (scat->oncePerParticle && dscatterGroup[scat->groupIndex].allScattered)
    return;

  if (iPass==scat->startOnPass) {
    scat->nLeft = scat->limitTotal>=0 ? scat->limitTotal : np;
    if (scat->nLeft>np)
      scat->nLeft = np;
  }
  nLeftThisPass = np;
  if (scat->limitTotal>=0)
    nLeftThisPass = scat->nLeft;
  if (scat->limitPerPass>=0) {
    if (scat->limitPerPass<nLeftThisPass)
      nLeftThisPass = scat->limitPerPass;
  }
  if (nLeftThisPass==0)
    return;

  if (scat->oncePerParticle && dscatterGroup[scat->groupIndex].nParticles<np) {
    fprintf(stderr, "Error for %s: number of particles is greater than the size of the particle ID array.\n",
            context.elementName);
    exit(1);
  }

  nScattered = 0;
  for (ip=0; ip<np; ip++) {
    if (scat->probability<1 && random_2(1)>scat->probability)
      continue;
    if (nLeftThisPass==0)
      break;
    if (scat->oncePerParticle) {
      short found = 0;
      if (dscatterGroup[scat->groupIndex].nScattered>=dscatterGroup[scat->groupIndex].nParticles) {
        fprintf(stderr, "All particles scattered for group %ld (nscattered=%ld, np0=%ld)\n",
                scat->group, dscatterGroup[scat->groupIndex].nScattered, 
		dscatterGroup[scat->groupIndex].nParticles);
	dscatterGroup[scat->groupIndex].allScattered = 1;
	return ;
      }
      for (i=0; i<dscatterGroup[scat->groupIndex].nScattered; i++) {
        if (dscatterGroup[scat->groupIndex].particleIDScattered[i]==part[ip][6]) {
          found = 1;
          break;
        }
      }
      if (found)
        continue;
      dscatterGroup[scat->groupIndex].particleIDScattered[i] = part[ip][6];
      dscatterGroup[scat->groupIndex].nScattered += 1;
    }
    nScattered++;
    nLeftThisPass--;
    cdf = random_2(1);
    amplitude = scat->factor*interp(scat->indepData, scat->cdfData, scat->nData, cdf, 0, 1, &interpCode);
    if (scat->randomSign)
      amplitude *= random_2(1)>0.5 ? 1 : -1;
    if (!interpCode)
      fprintf(stderr, "Warning: interpolation error for %s.  cdf=%e\n",
              context.elementName, cdf);
    if (scat->iPlane==5) {
      /* momentum scattering */
      P = (1+part[ip][5])*Po;
      beta = P/sqrt(sqr(P)+1);
      t = part[ip][4]/beta;
      part[ip][5] += amplitude;
      P = (1+part[ip][5])*Po;
      beta = P/sqrt(sqr(P)+1);
      part[ip][4] = t*beta;
    } else
      part[ip][scat->iPlane] += amplitude;
  }
  scat->nLeft -= nScattered;
  /*
    fprintf(stderr, "%ld particles scattered by %s#%ld, group %ld, total=%ld\n", 
    nScattered, context.elementName, context.elementOccurrence, 
    scat->group, dscatterGroup[scat->groupIndex].nScattered);
  */
}

void recordLossPass(long *lostOnPass, long *nLost, long nLeft, long nMaximum, 
                    long pass, int myid, int lostSinceSeqMode)
{
  long ip;
  if (!lostOnPass || !nLost)
    return;
  if (nMaximum==(nLeft+*nLost))
    /* no additional losses occurred */
    return;
  if ((!USE_MPI) || (USE_MPI && (myid != 0)) || (USE_MPI && (myid == 0) && lostSinceSeqMode)) {
    /* The information has been recorded in the gathering procedure for the master processor
       if it is a diagnostic element */
    for (ip=nLeft; ip<nMaximum-*nLost; ip++) {
      lostOnPass[ip] = pass;
    }
  }
  *nLost = nMaximum - nLeft;
}

void storeMonitorOrbitValues(ELEMENT_LIST *eptr, double **part, long np)
{
  MONI *moni;
  HMON *hmon;
  VMON *vmon;
  MARK *mark;
  char s[1000];
  double centroid[6];
  long i;
  
  if (!np)
    return;
  
  switch (eptr->type) {
  case T_MONI:
    moni = (MONI*)eptr->p_elem;
    if (!moni->coFitpoint) 
      return;
    compute_centroids(centroid, part, np);
    if (!moni->initialized) {
      sprintf(s, "%s#%ld.xco", eptr->name, eptr->occurence);
      moni->coMemoryNumber[0] = rpn_create_mem(s, 0);
      sprintf(s, "%s#%ld.xpco", eptr->name, eptr->occurence);
      moni->coMemoryNumber[1] = rpn_create_mem(s, 0);
      sprintf(s, "%s#%ld.yco", eptr->name, eptr->occurence);
      moni->coMemoryNumber[2] = rpn_create_mem(s, 0);
      sprintf(s, "%s#%ld.ypco", eptr->name, eptr->occurence);
      moni->coMemoryNumber[3] = rpn_create_mem(s, 0);
      moni->initialized = 1;
    }
    for (i=0; i<4; i++) 
      rpn_store(centroid[i], NULL, moni->coMemoryNumber[i]);
    break;
  case T_MARK:
    mark = (MARK*)eptr->p_elem;
    if (!mark->fitpoint) 
      return;
    compute_centroids(centroid, part, np);
    if (mark->co_mem==NULL) {
      mark->co_mem = tmalloc(4*sizeof(*(mark->co_mem)));
      sprintf(s, "%s#%ld.xco", eptr->name, eptr->occurence);
      mark->co_mem[0] = rpn_create_mem(s, 0);
      sprintf(s, "%s#%ld.xpco", eptr->name, eptr->occurence);
      mark->co_mem[1] = rpn_create_mem(s, 0);
      sprintf(s, "%s#%ld.yco", eptr->name, eptr->occurence);
      mark->co_mem[2] = rpn_create_mem(s, 0);
      sprintf(s, "%s#%ld.ypco", eptr->name, eptr->occurence);
      mark->co_mem[3] = rpn_create_mem(s, 0);
    }
    for (i=0; i<4; i++) 
      rpn_store(centroid[i], NULL, mark->co_mem[i]);
    break;
  case T_HMON:
    hmon = (HMON*)eptr->p_elem;
    if (!hmon->coFitpoint) 
      return;
    compute_centroids(centroid, part, np);
    if (!hmon->initialized) {
      sprintf(s, "%s#%ld.xco", eptr->name, eptr->occurence);
      hmon->coMemoryNumber[0] = rpn_create_mem(s, 0);
      sprintf(s, "%s#%ld.xpco", eptr->name, eptr->occurence);
      hmon->coMemoryNumber[1] = rpn_create_mem(s, 0);
      hmon->initialized = 1;
    }
    for (i=0; i<2; i++) 
      rpn_store(centroid[i], NULL, hmon->coMemoryNumber[i]);
    break;
  case T_VMON:
    vmon = (VMON*)eptr->p_elem;
    if (!vmon->coFitpoint) 
      return;
    compute_centroids(centroid, part, np);
    if (!vmon->initialized) {
      sprintf(s, "%s#%ld.yco", eptr->name, eptr->occurence);
      vmon->coMemoryNumber[0] = rpn_create_mem(s, 0);
      sprintf(s, "%s#%ld.ypco", eptr->name, eptr->occurence);
      vmon->coMemoryNumber[1] = rpn_create_mem(s, 0);
      vmon->initialized = 1;
    }
    for (i=0; i<2; i++) 
      rpn_store(centroid[i+2], NULL, vmon->coMemoryNumber[i]);
    break;
  default:
    return ;
    break;
  }
}
 
#if USE_MPI
void scatterParticles(double **coord, long *nToTrack, double **accepted,
                      long n_processors, int myid, balance balanceStatus, 
                      double my_rate, double nParPerElements, double round, 
                      int lostSinceSeqMode, int *distributed)
{
  long work_processors = n_processors-1; 
  int root = 0, i;
  int my_nToTrack, nItems, nToTrackCounts[n_processors], nRemainParticles;
  double total_rate, constTime;
  MPI_Status status;

  /* The particles will be distributed to slave processors evenly for the first pass */
  if ((balanceStatus==startMode) && (!*distributed))  {
    if (myid==0) 
      my_nToTrack = 0;
    else {
      my_nToTrack = *nToTrack/work_processors; 
      if (myid<=(*nToTrack%work_processors)) 
	my_nToTrack++;
    } 
    *distributed = 1; 
  }
  else if ((balanceStatus == badBalance) || lostSinceSeqMode) { 
   /* calculating the number of jobs to be sent according to the speed of each processors */
    MPI_Reduce(&my_rate, &total_rate, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (myid==0) {
      if (fabs(total_rate-0.0)>1e-12)
        constTime = *nToTrack/total_rate;
      else
        constTime = 0.0;
      my_nToTrack = 0;
#ifdef MPI_DEBUG
      printf("total_rate=%lf, nToTrack=%ld\n", total_rate, *nToTrack );
#endif
    }
    MPI_Bcast(&constTime, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (myid!=0) {
      *nToTrack = my_rate*constTime+round; /* round to the nearest integer */
      my_nToTrack = *nToTrack;
    }
  }
  else { /* keep the nToTrack unchanged */
    if (myid==0)
      my_nToTrack = 0;
    else       
      my_nToTrack =*nToTrack;
  }
 
  /* gather the number of particles to be sent to each processor */ 
  MPI_Gather(&my_nToTrack, 1, MPI_INT, nToTrackCounts, 1, MPI_INT, root, MPI_COMM_WORLD);

  /* scatter particles to all of the slave processors */
  if (myid==0) {
    for (i=1; i<work_processors; i++) {
#ifdef MPI_DEBUG
      if ((balanceStatus == badBalance) || lostSinceSeqMode)
	printf("%d will be computed on %d\n",nToTrackCounts[i],i);
#endif
      /* calculate the number of elements that will be sent to each processor */
      nItems = nToTrackCounts[i]*COORDINATES_PER_PARTICLE;
      MPI_Send (&coord[my_nToTrack][0], nItems, MPI_DOUBLE, i, 104, MPI_COMM_WORLD); 
      if (accepted!=NULL)
        MPI_Send (&accepted[my_nToTrack][0], nItems, MPI_DOUBLE, i, 105, MPI_COMM_WORLD); 
      /* count the total number of particles that have been scattered */
      my_nToTrack = my_nToTrack+nToTrackCounts[i];      
    }
  } 
  else if (myid!=work_processors) {
    MPI_Recv (&coord[0][0], my_nToTrack*COORDINATES_PER_PARTICLE, MPI_DOUBLE, 0,
              104, MPI_COMM_WORLD, &status); 
    if (accepted!=NULL)
      MPI_Recv (&accepted[0][0], my_nToTrack*COORDINATES_PER_PARTICLE, 
                MPI_DOUBLE, 0, 105, MPI_COMM_WORLD, &status);
  }
  
  /* The last processor will be responsible for all of the remaining particles */
  if (myid==0) {
    nRemainParticles = *nToTrack-my_nToTrack;
#ifdef MPI_DEBUG
    if (((balanceStatus==startMode) && (!*distributed)) || 
	(balanceStatus == badBalance) || lostSinceSeqMode)
      printf("%d will be computed on %ld\n",nRemainParticles,work_processors);
#endif
    nItems = nRemainParticles*COORDINATES_PER_PARTICLE;
    MPI_Send (&nRemainParticles, 1, MPI_INT, work_processors, 106, MPI_COMM_WORLD);
    MPI_Send (&coord[my_nToTrack][0], nItems, MPI_DOUBLE, work_processors, 104, MPI_COMM_WORLD); 
    if (accepted!=NULL)
      MPI_Send (&accepted[my_nToTrack][0], nItems, MPI_DOUBLE, work_processors, 105, MPI_COMM_WORLD);
  }
  else {
    if (myid==work_processors) {
      MPI_Recv (&my_nToTrack, 1, MPI_INT, 0, 106, MPI_COMM_WORLD, &status);
      MPI_Recv (&coord[0][0], my_nToTrack*COORDINATES_PER_PARTICLE, MPI_DOUBLE, 
                0, 104, MPI_COMM_WORLD, &status); 
      if (accepted!=NULL)
        MPI_Recv (&accepted[0][0], my_nToTrack*COORDINATES_PER_PARTICLE, 
                  MPI_DOUBLE, 0, 105, MPI_COMM_WORLD, &status);  
    }
    *nToTrack = my_nToTrack;
  }
}

void gatherParticles(double **coord, long *lostOnPass, long *nToTrack, long *nLost, double **accepted, long n_processors, int myid, double *round)
{
  long work_processors = n_processors-1;
  int root = 0, i, nItems, displs ;
  int my_nToTrack, my_nLost, nToTrackCounts[n_processors], nLostCounts[n_processors], current_nLost = 0; 
  MPI_Status status;

  if (myid==0) {
    my_nToTrack = 0;  
    my_nLost = 0;
  }
  else {
    my_nToTrack = *nToTrack;
    my_nLost = *nLost;
  }

  /* gather nToTrack and nLost from all of the slave processors to the master processors */ 
  MPI_Gather(&my_nToTrack, 1, MPI_INT, nToTrackCounts, 1, MPI_INT, root, MPI_COMM_WORLD);
  MPI_Gather(&my_nLost, 1, MPI_INT, nLostCounts, 1, MPI_INT, root, MPI_COMM_WORLD);

#ifdef SORT
  if ((myid!=0) && (my_nLost != 0)) {/* indicates more particles are lost, need sort */
    qsort(coord[0], my_nToTrack, COORDINATES_PER_PARTICLE*sizeof(double), comp_IDs);  
    if (accepted!=NULL)
      qsort(accepted[0], my_nToTrack, COORDINATES_PER_PARTICLE*sizeof(double), comp_IDs);  
  }
#endif

  if (myid==0) {
    for (i=1; i<=work_processors; i++) {
      /* the number of elements that are received from each processor (for root only) */
      nItems = nToTrackCounts[i]*COORDINATES_PER_PARTICLE;
      /* collect information for the left particles */
      MPI_Recv (&coord[my_nToTrack][0], nItems, MPI_DOUBLE, i, 100, MPI_COMM_WORLD, &status); 

      /* count the total number of particles to track and the total number of lost after the most recent scattering */
      my_nToTrack = my_nToTrack+nToTrackCounts[i];
      current_nLost = current_nLost+nLostCounts[i];      
    }
    *nLost = *nLost+current_nLost;
    *nToTrack = my_nToTrack;
  } 
  else {
    MPI_Send (&coord[0][0], my_nToTrack*COORDINATES_PER_PARTICLE, MPI_DOUBLE, 0, 100, MPI_COMM_WORLD); 
  }

  /* collect information for the lost particles and gather the accepted and lostOnPass arrays
     only when some particles are lost */
  
  MPI_Bcast(&current_nLost, 1, MPI_INT, root, MPI_COMM_WORLD);

  if (current_nLost>0) {
    if (myid==0) {
      /* set up the displacement array and the number of elements that are received from each processor */ 
      nLostCounts[0] = 0;
      displs = my_nToTrack;
      my_nToTrack = 0;
      for (i=1; i<=work_processors; i++) {
        /* gather information for lost particles */  
  	displs = displs+nLostCounts[i-1];
        nItems = nLostCounts[i]*COORDINATES_PER_PARTICLE;
        MPI_Recv (&coord[displs][0], nItems, MPI_DOUBLE, i, 102, MPI_COMM_WORLD, &status);
        if (accepted!=NULL){
          MPI_Recv (&accepted[my_nToTrack][0], nToTrackCounts[i]*COORDINATES_PER_PARTICLE, MPI_DOUBLE, i, 101, MPI_COMM_WORLD, &status); 
          MPI_Recv (&accepted[displs][0], nItems, MPI_DOUBLE, i, 103, MPI_COMM_WORLD, &status);
        my_nToTrack = my_nToTrack+nToTrackCounts[i];
	}
        if (lostOnPass!=NULL) {
          MPI_Recv (&lostOnPass[displs], nLostCounts[i], MPI_LONG, i, 107, 
                    MPI_COMM_WORLD, &status);
        }  
      }
      /* update the round parameter to avoid more particles are distributed 
         than the available particles */
      if ((my_nToTrack/work_processors)<work_processors)
        *round = 0.0; 
    }
    else {
      /* send information for lost particles */
      MPI_Send (&coord[my_nToTrack][0], my_nLost*COORDINATES_PER_PARTICLE, MPI_DOUBLE, root, 102, MPI_COMM_WORLD);  
      if (accepted!=NULL) {
        MPI_Send (&accepted[0][0], my_nToTrack*COORDINATES_PER_PARTICLE, MPI_DOUBLE, root, 101, MPI_COMM_WORLD);  
        MPI_Send (&accepted[my_nToTrack][0], my_nLost*COORDINATES_PER_PARTICLE, MPI_DOUBLE, root, 103, MPI_COMM_WORLD); 
      }
      if (lostOnPass!=NULL)
        MPI_Send (&lostOnPass[my_nToTrack], my_nLost, MPI_LONG, root, 107, MPI_COMM_WORLD);  
      /* The number of lost particles on the slave processor will be set to 0 after gathering */ 
      *nLost = 0;     
    }
    MPI_Bcast (round, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
  } 
}

balance checkBalance (double my_wtime, int myid, long n_processors)
{
  double maxTime, minTime, time[n_processors];
  int i, balanceFlag = -1; 

  MPI_Gather (&my_wtime, 1, MPI_DOUBLE, time, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
  if (myid==0) {
    maxTime = minTime = time[1];
    for (i=2; i<n_processors; i++) {
      if (maxTime<time[i])
        maxTime = time[i];
      if (minTime>time[i])
        minTime = time[i]; 
    } 
    if ((maxTime-minTime)/minTime>0.10)
      balanceFlag = 0;        /* the workload balancing is bad */  
    else             
      balanceFlag = 1;  /* the workload balancing is good */
  }
  MPI_Bcast (&balanceFlag, 1, MPI_INT, 0, MPI_COMM_WORLD);     

#ifdef MPI_DEBUG
  if (balanceFlag==1) 
    printf("The balance is in good status.\n");
  else
    printf("The balance is in bad status.\n");
#endif

  if (balanceFlag==1) 
    return goodBalance;
  else
    return badBalance;
}   

#endif

#ifdef SORT
  int comp_IDs(const void *coord1, const void *coord2)
  {
    if (((double*) coord1)[6] <((double*) coord2)[6]) 
      return -1;
    else if (((double*) coord1)[6] > ((double*) coord2)[6]) 
      return  1;
    else 
      return 0;
  }
#endif
