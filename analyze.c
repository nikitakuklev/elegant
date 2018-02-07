/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: analyze.c
 * purpose: Do tracking to find transfer matrix and tunes.
 *          See file analyze.nl for input parameters.
 *
 * Michael Borland, 1989
 */
#include "mdb.h"
#include "track.h"
#include "analyze.h"
#include "matlib.h"

#define CMATRIX_OFFSET 0
#define RMATRIX_OFFSET CMATRIX_OFFSET+6
#define X_BETA_OFFSET RMATRIX_OFFSET+36
#define Y_BETA_OFFSET X_BETA_OFFSET+6
#define DETR_OFFSET Y_BETA_OFFSET+6
#define X_ETA_OFFSET DETR_OFFSET+1
#define Y_ETA_OFFSET X_ETA_OFFSET+2
#define CLORB_ETA_OFFSET Y_ETA_OFFSET+2
#define N_ANALYSIS_COLUMNS CLORB_ETA_OFFSET+4


SDDS_DEFINITION analysis_column[N_ANALYSIS_COLUMNS] = {
    {"C1", "&column name=C1, units=m, symbol=\"C$b1$n\", type=double &end"},
    {"C2", "&column name=C2, symbol=\"C$b2$n\", type=double &end"},
    {"C3", "&column name=C3, units=m, symbol=\"C$b3$n\", type=double &end"},
    {"C4", "&column name=C4, symbol=\"C$b4$n\", type=double &end"},
    {"C5", "&column name=C5, units=m, symbol=\"C$b5$n\", type=double &end"},
    {"C6", "&column name=C6, symbol=\"C$b6$n\", type=double &end"},
    {"R11", "&column name=R11, symbol=\"R$b11$n\", type=double &end"},
    {"R12", "&column name=R12, units=m, symbol=\"R$b12$n\", type=double &end"},
    {"R13", "&column name=R13, symbol=\"R$b13$n\", type=double &end"},
    {"R14", "&column name=R14, units=m, symbol=\"R$b14$n\", type=double &end"},
    {"R15", "&column name=R15, symbol=\"R$b15$n\", type=double &end"},
    {"R16", "&column name=R16, units=m, symbol=\"R$b16$n\", type=double &end"},
    {"R21", "&column name=R21, units=1/m, symbol=\"R$b21$n\", type=double &end"},
    {"R22", "&column name=R22, symbol=\"R$b22$n\", type=double &end"},
    {"R23", "&column name=R23, units=1/m, symbol=\"R$b23$n\", type=double &end"},
    {"R24", "&column name=R24, symbol=\"R$b24$n\", type=double &end"},
    {"R25", "&column name=R25, units=1/m, symbol=\"R$b25$n\", type=double &end"},
    {"R26", "&column name=R26, symbol=\"R$b26$n\", type=double &end"},
    {"R31", "&column name=R31, symbol=\"R$b31$n\", type=double &end"},
    {"R32", "&column name=R32, units=m, symbol=\"R$b32$n\", type=double &end"},
    {"R33", "&column name=R33, symbol=\"R$b33$n\", type=double &end"},
    {"R34", "&column name=R34, units=m, symbol=\"R$b34$n\", type=double &end"},
    {"R35", "&column name=R35, symbol=\"R$b35$n\", type=double &end"},
    {"R36", "&column name=R36, units=m, symbol=\"R$b36$n\", type=double &end"},
    {"R41", "&column name=R41, units=1/m, symbol=\"R$b41$n\", type=double &end"},
    {"R42", "&column name=R42, symbol=\"R$b42$n\", type=double &end"},
    {"R43", "&column name=R43, units=1/m, symbol=\"R$b43$n\", type=double &end"},
    {"R44", "&column name=R44, symbol=\"R$b44$n\", type=double &end"},
    {"R45", "&column name=R45, units=1/m, symbol=\"R$b45$n\", type=double &end"},
    {"R46", "&column name=R46, symbol=\"R$b46$n\", type=double &end"},
    {"R51", "&column name=R51, symbol=\"R$b51$n\", type=double &end"},
    {"R52", "&column name=R52, units=m, symbol=\"R$b52$n\", type=double &end"},
    {"R53", "&column name=R53, symbol=\"R$b53$n\", type=double &end"},
    {"R54", "&column name=R54, units=m, symbol=\"R$b54$n\", type=double &end"},
    {"R55", "&column name=R55, symbol=\"R$b55$n\", type=double &end"},
    {"R56", "&column name=R56, units=m, symbol=\"R$b56$n\", type=double &end"},
    {"R61", "&column name=R61, units=1/m, symbol=\"R$b61$n\", type=double &end"},
    {"R62", "&column name=R62, symbol=\"R$b62$n\", type=double &end"},
    {"R63", "&column name=R63, units=1/m, symbol=\"R$b63$n\", type=double &end"},
    {"R64", "&column name=R64, symbol=\"R$b64$n\", type=double &end"},
    {"R65", "&column name=R65, units=1/m, symbol=\"R$b65$n\", type=double &end"},
    {"R66", "&column name=R66, symbol=\"R$b66$n\", type=double &end"},
    {"betax", "&column name=betax, units=m, symbol=\"$gb$r$bx$n\", type=double &end"},
    {"alphax", "&column name=alphax, units=m, symbol=\"$gb$r$bx$n\", type=double &end"},
    {"nux", "&column name=nux, symbol=\"$gn$r$bx$n\", type=double &end"},
    {"dbetax/dp", "&column name=dbetax/dp, units=m, type=double, description=\"Derivative of betax with momentum offset\" &end"},
    {"dalphax/dp", "&column name=dalphax/dp, type=double, description=\"Derivative of alphax with momentum offset\" &end"},
    {"dnux/dp", "&column name=dnux/dp, symbol=\"$gx$r$bx$n\", type=double, description=\"Horizontal chromaticity\" &end"},
    {"betay", "&column name=betay, units=m, symbol=\"$gb$r$by$n\", type=double &end"},
    {"alphay", "&column name=alphay, units=m, symbol=\"$gb$r$bx$n\", type=double &end"},
    {"nuy", "&column name=nuy, symbol=\"$gn$r$by$n\", type=double &end"},
    {"dbetay/dp", "&column name=dbetay/dp, units=m, type=double, description=\"Derivative of betay with momentum offset\" &end"},
    {"dalphay/dp", "&column name=dalphay/dp, type=double, description=\"Derivative of alphay with momentum offset\" &end"},
    {"dnuy/dp", "&column name=dnuy/dp, symbol=\"$gx$r$by$n\", type=double, description=\"Vertical chromaticity\" &end"},
    {"detR", "&column name=detR, symbol=\"detR\", type=double &end"},
    {"etax", "&column name=etax, units=m, symbol=\"$gc$r$bx$n\", type=double &end"},
    {"etapx", "&column name=etapx, symbol=\"$gc$r$bx$n$a'$n\", type=double &end"},
    {"etay", "&column name=etay, units=m, symbol=\"$gc$ry\", type=double &end"},
    {"etapy", "&column name=etapy, symbol=\"$gc$r$by$n$a'$n\", type=double &end"},
    {"cetax", "&column name=cetax, units=m, symbol=\"$gc$r$bxc$n\", type=double &end"},
    {"cetapx", "&column name=cetapx, symbol=\"$gc$r$bxc$n$a'$n\", type=double &end"},
    {"cetay", "&column name=cetay, units=m, symbol=\"$gc$r$byc$n\", type=double &end"},
    {"cetapy", "&column name=cetapy, symbol=\"$gc$r$byc$n$a'$n\", type=double &end"},
    } ;

#define IP_STEP 0
#define IP_SUBSTEP 1
#define N_PARAMETERS 2
static SDDS_DEFINITION parameter_definition[N_PARAMETERS] = {
    {"Step", "&parameter name=Step, type=long, description=\"Simulation step\" &end"},
    {"Substep", "&parameter name=Substep, type=long, description=\"Simulation substep\" &end"},
    } ;

static long SDDS_analyze_initialized = 0;
static SDDS_TABLE SDDS_analyze;
static FILE *fpPrintout = NULL;

typedef struct {
  double tune1[2], beta1[2], alpha1[2];
} CHROM_DERIVS;

void copyParticles(double ***coordCopy, double **coord, long np);
void performChromaticAnalysisFromMap(VMATRIX *M, TWISS *twiss, CHROM_DERIVS *chromDeriv);
void printMapAnalysisResults(FILE *fp, long printoutOrder, VMATRIX *M, TWISS *twiss, CHROM_DERIVS *chromDeriv, double *data);
void propagateTwissParameters(TWISS *twiss1, TWISS *twiss0, VMATRIX *M);
  
void setup_transport_analysis(
    NAMELIST_TEXT *nltext,
    RUN *run,
    VARY *control,
    ERRORVAL *errcon
    )
{
    log_entry("setup_transport_analysis");

    /* process namelist input */
    set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
    set_print_namelist_flags(0);
    if (processNamelist(&analyze_map, nltext)==NAMELIST_ERROR)
      bombElegant(NULL, NULL);
    if (echoNamelists) print_namelist(stdout, &analyze_map);

    /* check for data errors */
    if (!output && !printout)
        bombElegant("no output filename or printout filename specified", NULL);
    if (printout_order>3) {
      printf("warning: maximum printout_order is 3\n");
      printout_order = 3;
    }
    
#if USE_MPI
      if (myid==0) {
	/* In MPI mode, all output is handled by the master processor */
#endif
    if (SDDS_analyze_initialized) {
      if (SDDS_IsActive(&SDDS_analyze) && !SDDS_Terminate(&SDDS_analyze)) {
        SDDS_SetError("Problem terminating SDDS output (finish_transport_analysis)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      }
      SDDS_analyze_initialized = 0;
    }

    if (output) {
      output = compose_filename(output, run->rootname);
      SDDS_ElegantOutputSetup(&SDDS_analyze, output, SDDS_BINARY, 1, "transport analysis", 
			      run->runfile, run->lattice, parameter_definition, N_PARAMETERS,
			      analysis_column, N_ANALYSIS_COLUMNS, "setup_transport_analysis", 
			      SDDS_EOS_NEWFILE);
      SDDS_analyze_initialized = 1;
      
      if (!SDDS_DefineSimpleColumns(&SDDS_analyze, control->n_elements_to_vary,
				    control->varied_quan_name, control->varied_quan_unit, SDDS_DOUBLE) ||
	  !SDDS_DefineSimpleColumns(&SDDS_analyze, errcon->n_items, errcon->quan_name, errcon->quan_unit, 
				    SDDS_DOUBLE)) {
        SDDS_SetError("Unable to define additional SDDS columns (setup_transport_analysis)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      }
      
      if (!SDDS_SaveLayout(&SDDS_analyze) || !SDDS_WriteLayout(&SDDS_analyze)) {
        SDDS_SetError("Unable to write SDDS layout for transport analysis");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      }
    }

    if (fpPrintout)
      fclose(fpPrintout);
    fpPrintout = NULL;
    if (printout) {
      printout = compose_filename(printout, run->rootname);
      fpPrintout = fopen(printout, "w");
    }
#if USE_MPI
      }
#endif
    
    log_exit("setup_transport_analysis");
    }

void do_transport_analysis(
    RUN *run,
    VARY *control,
    ERRORVAL *errcon,
    LINE_LIST *beamline,
    double *orbit
    )
{
    double **initialCoord, **finalCoord, **coordError, **finalCoordCopy;
    double *data, *offset, *orbit_p, *orbit_m;
    MATRIX *R;
    double p_central;
    long n_left, n_track, i, j, k, index, code;
    TRAJECTORY *clorb=NULL;
    double stepSize[6], maximumValue[6], sum2Difference, maxAbsDifference, difference;
    VMATRIX *M;
    TWISS twiss;
    CHROM_DERIVS chromDeriv;
#if USE_MPI
    long *nToTrackCounts, my_nTrack, my_offset, n_leftTotal, nItems;
    MPI_Status status;
    notSinglePart = 0;
#endif

    if (center_on_orbit && !orbit)
        bombElegant("you've asked to center the analysis on the closed orbit, but you didn't issue a closed_orbit command", NULL);

    data  = (double*)tmalloc(sizeof(*data)*SDDS_analyze.layout.n_columns);
    offset = (double*)tmalloc(sizeof(*offset)*6);
    m_alloc(&R , 6, 6);
    orbit_p = tmalloc(sizeof(*orbit_p)*6);
    orbit_m = tmalloc(sizeof(*orbit_m)*6);

    if (center_on_orbit)
        clorb = tmalloc(sizeof(*clorb)*(beamline->n_elems+1));

    stepSize[0] = delta_x;
    stepSize[1] = delta_xp;
    stepSize[2] = delta_y;
    stepSize[3] = delta_yp;
    stepSize[4] = delta_s;
    stepSize[5] = delta_dp;

    n_track = 
      makeInitialParticleEnsemble(&initialCoord,
				  (orbit && center_on_orbit? orbit: NULL),
				  &finalCoord, &coordError,
				  7, stepSize);
#if DEBUG
    if (1) {
      long i, j;
      FILE *fp;
      fp = fopen("analyzeMap.in", "w");
      fprintf(fp, "SDDS1\n");
      fprintf(fp, "&column name=x0 type=double units=m &end\n");
      fprintf(fp, "&column name=xp0 type=double &end\n");
      fprintf(fp, "&column name=y0 type=double units=m &end\n");
      fprintf(fp, "&column name=yp0 type=double &end\n");
      fprintf(fp, "&column name=s0 type=double units=m &end\n");
      fprintf(fp, "&column name=delta0 type=double &end\n");
      fprintf(fp, "&data mode=ascii no_row_counts=1 &end\n");
      for (i=0; i<n_track; i++) {
        for (j=0; j<6; j++)
          fprintf(fp, "%21.15e ", initialCoord[i][j]);
        fprintf(fp, "\n");
      }
    }    
#endif

    if (canonical_variables)
      /* Assume particles are generated in canonical variables (x, px, y, py, -s, delta) and
       * convert to (x, x', y, y', s, delta) coordinates for tracking.
       */
      convertFromCanonicalCoordinates(finalCoord, n_track, run->p_central, 1);

    /* Track the reference particle for fiducialization. In MPI mode, all cores do this */
    beamline->fiducial_flag = 0;
    p_central = run->p_central;
    if (verbosity>0) 
      printf("Tracking fiducial particle\n");
    code = do_tracking(NULL, finalCoord, 1, NULL, beamline, &p_central, 
                       NULL, NULL, NULL, NULL, run, control->i_step, 
                       FIRST_BEAM_IS_FIDUCIAL+(verbosity>1?0:SILENT_RUNNING)+INHIBIT_FILE_OUTPUT,
		       1, 0, NULL, NULL, NULL, NULL, NULL);

    if (!code) {
      printf("Fiducial particle lost. Don't know what to do.\n");
      exitElegant(1);
    }
    if (verbosity>0)
      printf("Fiducialization completed\n");
    
    /* Track the other particles. */
    p_central = run->p_central;
#if USE_MPI
    /* We partition by processor in MPI mode.
       The arrays are oversized, but not so large that it will hurt. 
    */
    n_left = n_track - 1;
    my_nTrack = n_left/n_processors+0.5;
    my_offset = myid*my_nTrack + 1;
    if (myid==(n_processors-1))
      my_nTrack = n_left - my_offset + 1;
#if MPI_DEBUG
    printf("Tracking %ld particles, offset %ld, processor %d\n",
	   my_nTrack, my_offset, myid);
    fflush(stdout);
#endif    
    n_left = do_tracking(NULL, finalCoord+my_offset, my_nTrack, NULL, beamline, &p_central, 
			 NULL, NULL, NULL, NULL, run, control->i_step, 
			 FIDUCIAL_BEAM_SEEN+FIRST_BEAM_IS_FIDUCIAL+SILENT_RUNNING+INHIBIT_FILE_OUTPUT,
			 1, 0, NULL, NULL, NULL, NULL, NULL);
#if MPI_DEBUG
    printf("Tracking done\n"); fflush(stdout);
#endif
    MPI_Allreduce(&n_left, &n_leftTotal, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
    if (n_leftTotal!=(n_track-1)) {
      printf("warning: %ld particle(s) lost during transport analysis--continuing with next step", n_track-1-n_leftTotal);
      return;
    }
#else
    n_left = do_tracking(NULL, finalCoord+1, n_track-1, NULL, beamline, &p_central, 
			 NULL, NULL, NULL, NULL, run, control->i_step, 
			 FIDUCIAL_BEAM_SEEN+FIRST_BEAM_IS_FIDUCIAL+SILENT_RUNNING+INHIBIT_FILE_OUTPUT,
			 1, 0, NULL, NULL, NULL, NULL, NULL);
    if (n_left!=(n_track-1)) {
      fputs("warning: particle(s) lost during transport analysis--continuing with next step", stdout);
      return;
    }
#endif

#if USE_MPI
    /* Gather final particles back to master */
    MPI_Barrier(MPI_COMM_WORLD);
    nToTrackCounts = tmalloc(sizeof(long)*n_processors);
    MPI_Gather(&my_nTrack, 1, MPI_LONG, nToTrackCounts, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    if (myid==0) {
      /* Copy data from each slave */
      my_nTrack += my_offset; /* to account for the fiducial particle */
      for (i=1; i<n_processors; i++) {
	if (verbosity>2) {
	  printf("Pulling %ld particles from processor %ld\n", nToTrackCounts[i], i);
	  fflush(stdout);
	}
	nItems = nToTrackCounts[i]*COORDINATES_PER_PARTICLE;
	MPI_Recv(&finalCoord[my_nTrack][0], nItems, MPI_DOUBLE, i, 100, MPI_COMM_WORLD, &status); 
	my_nTrack += nToTrackCounts[i];
      }
    } else {
      /* Send data to master */
      MPI_Send (&finalCoord[my_offset][0], my_nTrack*COORDINATES_PER_PARTICLE, MPI_DOUBLE, 0, 100, MPI_COMM_WORLD);
    }
    free(nToTrackCounts);
#endif
    
#if USE_MPI
    if (myid==0) {
      /* In MPI mode, only master does analysis and output */
#endif
    
#if DEBUG
    if (1) {
      long i, j;
      FILE *fp;
      fp = fopen("analyzeMap.out", "w");
      fprintf(fp, "SDDS1\n");
      fprintf(fp, "&column name=x type=double units=m &end\n");
      fprintf(fp, "&column name=xp type=double &end\n");
      fprintf(fp, "&column name=y type=double units=m &end\n");
      fprintf(fp, "&column name=yp type=double &end\n");
      fprintf(fp, "&column name=s type=double units=m &end\n");
      fprintf(fp, "&column name=delta type=double &end\n");
      fprintf(fp, "&data mode=ascii no_row_counts=1 &end\n");
      for (i=0; i<n_track; i++) {
        for (j=0; j<6; j++)
          fprintf(fp, "%21.15e ", finalCoord[i][j]);
        fprintf(fp, "\n");
      }
    }
#endif
    if (canonical_variables)
      /* convert back to (x, px, y, py, -s, delta) for analysis */
      convertToCanonicalCoordinates(finalCoord, n_track, run->p_central, 1);
    
    copyParticles(&finalCoordCopy, finalCoord, n_track);
    if (p_central!=run->p_central && verbosity>0)
      printf("Central momentum changed from %e to %e\n", run->p_central, p_central);

    if (verbosity>0){
        if (orbit && center_on_orbit) {
            printf("closed orbit: \n");
            fflush(stdout);
            for (i=0; i<6; i++)
                printf("%16.8e ", orbit[i]);
            fputc('\n', stdout);
            }
        printf("final coordinates of refence particle: \n");
        fflush(stdout);
        for (i=0; i<6; i++)
            printf("%16.8e ", finalCoord[0][i]);
        fputc('\n', stdout);
        fflush(stdout);
        if (verbosity>1) {
          for (i=0; i<n_track-1; i++) {
            printf("Particle %ld start/end/delta coordinates:\n", i);
            for (j=0; j<6; j++)
              printf("%21.15e%c", initialCoord[i][j], j==5?'\n':' ');
            for (j=0; j<6; j++)
              printf("%21.15e%c", finalCoord[i][j], j==5?'\n':' ');
            for (j=0; j<6; j++)
              printf("%21.15e%c", finalCoord[i][j]-finalCoord[0][j], j==5?'\n':' ');
          }
          fflush(stdout);
        }
      }
    

    for (i=0; i<6; i++) {
      maximumValue[i] = -DBL_MAX;
      for (j=0; j<n_track; j++) {
	if (maximumValue[i]<initialCoord[j][i])
	  maximumValue[i] = initialCoord[j][i];
      }
    }

    /* Set errors as fraction of the absolute maximum coordinate */
    for (i=0; i<6; i++) {
      double min, max;
      max = -(min = DBL_MAX);
      for (j=0; j<n_track; j++) {
        if (max<finalCoord[j][i])
          max = finalCoord[j][i];
        if (min>finalCoord[j][i])
          min = finalCoord[j][i];
      }
      max = fabs(max);
      min = fabs(min);
      max = max>min ? max : min;
      for (j=0; j<n_track; j++)
        coordError[j][i] = max*accuracy_factor;
    }
    
    M = computeMatricesFromTracking(stdout, initialCoord, finalCoord, coordError, stepSize,
				    maximumValue, 7, n_track, 8, verbosity>1?1:0);
    
    for (i=0; i<6; i++) 
      data[i+CMATRIX_OFFSET] = M->C[i];
    for (i=k=0; i<6; i++)
      for (j=0; j<6; j++, k++)
	data[k+RMATRIX_OFFSET] = M->R[i][j];

    performChromaticAnalysisFromMap(M, &twiss, &chromDeriv);

    data[X_BETA_OFFSET  ] = twiss.betax;
    data[X_BETA_OFFSET+1] = twiss.alphax;
    data[X_BETA_OFFSET+2] = twiss.phix/PIx2;
    data[X_BETA_OFFSET+3] = chromDeriv.beta1[0];
    data[X_BETA_OFFSET+4] = chromDeriv.alpha1[0];
    data[X_BETA_OFFSET+5] = chromDeriv.tune1[0];
    data[X_ETA_OFFSET   ] = twiss.etax;
    data[X_ETA_OFFSET+1 ] = twiss.etapx;

    data[Y_BETA_OFFSET  ] = twiss.betay;
    data[Y_BETA_OFFSET+1] = twiss.alphay;
    data[Y_BETA_OFFSET+2] = twiss.phiy/PIx2;
    data[Y_BETA_OFFSET+3] = chromDeriv.beta1[1];
    data[Y_BETA_OFFSET+4] = chromDeriv.alpha1[1];
    data[Y_BETA_OFFSET+5] = chromDeriv.tune1[1];
    data[Y_ETA_OFFSET   ] = twiss.etay;
    data[Y_ETA_OFFSET+1 ] = twiss.etapy;

    /* compute determinant of R */
    for (i=0; i<6; i++)
      for (j=0; j<6; j++)
	R->a[i][j] = M->R[i][j];
    data[DETR_OFFSET] = m_det(R);

    if (delta_dp && center_on_orbit) {
        if (!beamline->matrix)
            beamline->matrix = full_matrix(&(beamline->elem), run, 1);
        find_closed_orbit(clorb, 1e-12, 20, beamline, beamline->matrix, run, delta_dp, 1, 0, NULL, 0.5,
                          1.05, 5, NULL, 0);
        for (i=0; i<6; i++)
            orbit_p[i] = clorb[0].centroid[i];
        find_closed_orbit(clorb, 1e-12, 20, beamline, beamline->matrix, run, -delta_dp, 1, 0, NULL, 0.5,
                          1.05, 5, NULL, 0);
        for (i=0; i<6; i++)
            orbit_m[i] = clorb[0].centroid[i];
        for (i=0; i<4; i++)
            data[CLORB_ETA_OFFSET+i] = (orbit_p[i]-orbit_m[i])/(2*delta_dp);
        }
    else
        for (i=0; i<4; i++)
            data[CLORB_ETA_OFFSET+i] = 0;

    index = N_ANALYSIS_COLUMNS;
    for (i=0; i<control->n_elements_to_vary; i++, index++)
        data[index] = control->varied_quan_value[i];
    for (i=0 ; i<errcon->n_items; i++, index++)
        data[index] = errcon->error_value[i];

    if (!SDDS_StartTable(&SDDS_analyze, 1)) {
        printf("Unable to start SDDS table (do_transport_analysis)");
        fflush(stdout);
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
        exitElegant(1);
        }
    if (!SDDS_SetParameters(&SDDS_analyze, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, 
                            IP_STEP, control->i_step, IP_SUBSTEP, control->i_vary, -1)) {
        SDDS_SetError("Unable to set SDDS parameter values (do_transport_analysis)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
    for (i=0; i<SDDS_ColumnCount(&SDDS_analyze); i++)
        if (!SDDS_SetRowValues(&SDDS_analyze, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, 0,
                               i, data[i], -1)) {
            printf("Unable to set SDDS column %s (do_transport_analysis)\n", 
                    analysis_column[i].name);
            fflush(stdout);
            SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
            exitElegant(1);
            }
    if (!SDDS_WriteTable(&SDDS_analyze)) {
        printf("Unable to write SDDS table (do_transport_analysis)");
        fflush(stdout);
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
        exitElegant(1);
        }
    if (!inhibitFileSync)
      SDDS_DoFSync(&SDDS_analyze);

    if (fpPrintout)
      printMapAnalysisResults(fpPrintout, printout_order, M, &twiss, &chromDeriv, data);
    if (verbosity>0)
      printMapAnalysisResults(stdout, printout_order, M, &twiss, &chromDeriv, data);

    /* check accuracy of matrix in reproducing coordinates */
    /* all processors do this for all particles, which shouldn't take much time */
    track_particles(initialCoord, M, initialCoord, n_track);
    for (k=0; k<6; k++) {
      sum2Difference = maxAbsDifference = 0;
      for (i=0; i<n_track; i++) {
	difference = fabs(initialCoord[i][k] - finalCoordCopy[i][k]);
	if (difference > maxAbsDifference)
	  maxAbsDifference = difference;
	sum2Difference += sqr(difference);
      }
      printf("Error for coordinate %ld: maximum = %le, rms = %le\n",
	     k, maxAbsDifference, sqrt(sum2Difference/n_track));
    }

#if USE_MPI
    }
#endif
    
#if USE_MPI
    if (myid==0) {
#endif
      free_matrices(M);
      free(M);
      free_czarray_2d((void**)finalCoordCopy, n_track, COORDINATES_PER_PARTICLE);
#if USE_MPI
    }
#endif
    
    free_czarray_2d((void**)initialCoord, n_track, COORDINATES_PER_PARTICLE);
    free_czarray_2d((void**)finalCoord, n_track, COORDINATES_PER_PARTICLE);
    free_czarray_2d((void**)coordError, n_track, COORDINATES_PER_PARTICLE);
    free(data);
    free(offset);
    free(orbit_p);
    free(orbit_m);
    m_free(&R);
    }


void printMapAnalysisResults(FILE *fp, long printoutOrder, VMATRIX *M, TWISS *twiss, CHROM_DERIVS *chromDeriv, double *data)
{
  long saveOrder;
  saveOrder = M->order;
  M->order = printout_order>3 ? 3 : printout_order;
  print_matrices(fp, "Matrix from fitting:", M);
  M->order = saveOrder;
  fprintf(fp, "determinant of R = 1 + %14.6e\n", data[DETR_OFFSET]-1);
  if (delta_dp && center_on_orbit) 
    fprintf(fp, "dispersion functions from closed-orbit calculations:\nx: %e m    %e\ny: %e m    %e\n",
	    data[CLORB_ETA_OFFSET  ], data[CLORB_ETA_OFFSET+1],
	    data[CLORB_ETA_OFFSET+2], data[CLORB_ETA_OFFSET+3]);
  if (periodic)
    fprintf(fp, "Lattice functions computed on assumption of periodic system:\n");
  else
    fprintf(fp, "Lattice functions computed on assumption of non-periodic system:\n");
  fprintf(fp, " horizontal:   tune = %14.6e  beta = %14.6e alpha = %14.6e  eta = %14.6e  eta' = %14.6e\n",
          twiss->phix/PIx2, twiss->betax, twiss->alphax, twiss->etax, twiss->etapx);
  fprintf(fp, "               dnu/dp = %14.6e  dbeta/dp = %14.6e  dalpha/dp = %14.6e\n",
          chromDeriv->tune1[0], chromDeriv->beta1[0], chromDeriv->alpha1[0]);
  fflush(fp);
  fprintf(fp, " vertical  :   tune = %14.6e  beta = %14.6e alpha = %14.6e  eta = %14.6e  eta' = %14.6e\n",
          twiss->phiy/PIx2, twiss->betay, twiss->alphay, twiss->etay, twiss->etapy);
  fprintf(fp, "               dnu/dp = %14.6e  dbeta/dp = %14.6e  dalpha/dp = %14.6e\n",
          chromDeriv->tune1[1], chromDeriv->beta1[1], chromDeriv->alpha1[1]);
  fflush(fp);
}


void finish_transport_analysis(
    RUN *run,
    VARY *control,
    ERRORVAL *errcon,
    LINE_LIST *beamline
    )
{
    if (SDDS_IsActive(&SDDS_analyze) && !SDDS_Terminate(&SDDS_analyze)) {
        SDDS_SetError("Problem terminating SDDS output (finish_transport_analysis)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }        
    SDDS_analyze_initialized = 0;
    if (fpPrintout)
      fclose(fpPrintout);
    }


VMATRIX *determineMatrix(RUN *run, ELEMENT_LIST *eptr, double *startingCoord, double *stepSize)
{
  double **coord;
  long n_track, i, j;
  VMATRIX *M;
  double **R, *C;
  double defaultStep[6] = {1e-4, 1e-4, 1e-4, 1e-4, 1e-4, 1e-4};
  long ltmp1, ltmp2;
  double dgamma, dtmp1, dP[3];
 
#if USE_MPI
  long notSinglePart_saved = notSinglePart;

  /* All the particles should do the same thing for this routine. */	
  notSinglePart = 0;
#endif
   		 
  coord = (double**)czarray_2d(sizeof(**coord), 1+6*4, COORDINATES_PER_PARTICLE);

  if (stepSize==NULL)
    stepSize = defaultStep;
  
  n_track = 4*6+1;
  for (j=0; j<COORDINATES_PER_PARTICLE; j++)
    for (i=0; i<n_track; i++)
      coord[i][j] = startingCoord ? startingCoord[j] : 0;

  /* particles 0 and 1 are for d/dx */
  coord[0][0] += stepSize[0] ;
  coord[1][0] -= stepSize[0] ;
  /* particles 2 and 3 are for d/dxp */
  coord[2][1] += stepSize[1];
  coord[3][1] -= stepSize[1];
  /* particles 4 and 5 are for d/dy */
  coord[4][2] += stepSize[2] ;
  coord[5][2] -= stepSize[2] ;
  /* particles 6 and 7 are for d/dyp */
  coord[6][3] += stepSize[3] ;
  coord[7][3] -= stepSize[3] ;
  /* particles 8 and 9 are for d/ds */
  coord[8][4] += stepSize[4] ;
  coord[9][4] -= stepSize[4] ;
  /* particles 10 and 11 are for d/delta */
  coord[10][5] += stepSize[5];
  coord[11][5] -= stepSize[5];

  /* particles 12 and 13 are for d/dx */
  coord[12][0] += 3*stepSize[0] ;
  coord[13][0] -= 3*stepSize[0] ;
  /* particles 14 and 15 are for d/dxp */
  coord[14][1] += 3*stepSize[1];
  coord[15][1] -= 3*stepSize[1];
  /* particles 16 and 17 are for d/dy */
  coord[16][2] += 3*stepSize[2] ;
  coord[17][2] -= 3*stepSize[2] ;
  /* particles 18 and 19 are for d/dyp */
  coord[18][3] += 3*stepSize[3] ;
  coord[19][3] -= 3*stepSize[3] ;
  /* particles 20 and 21 are for d/ds */
  coord[20][4] += 3*stepSize[4] ;
  coord[21][4] -= 3*stepSize[4] ;
  /* particles 22 and 23 are for d/delta */
  coord[22][5] += 3*stepSize[5];
  coord[23][5] -= 3*stepSize[5];
  /* particle n_track-1 is the reference particle (coordinates set above) */

  switch (eptr->type) {
  case T_CSBEND:
    ltmp1 = ((CSBEND*)eptr->p_elem)->isr;
    ltmp2 = ((CSBEND*)eptr->p_elem)->synch_rad;
    ((CSBEND*)eptr->p_elem)->isr = ((CSBEND*)eptr->p_elem)->synch_rad = 0;
    track_through_csbend(coord, n_track, (CSBEND*)eptr->p_elem, 0.0, run->p_central, NULL, 0.0,
                         NULL, NULL, NULL, NULL);
    ((CSBEND*)eptr->p_elem)->isr = ltmp1;
    ((CSBEND*)eptr->p_elem)->synch_rad = ltmp2;
   break;
  case T_CRBEND:
    ltmp1 = ((CRBEND*)eptr->p_elem)->isr;
    ltmp2 = ((CRBEND*)eptr->p_elem)->synch_rad;
    ((CRBEND*)eptr->p_elem)->isr = ((CRBEND*)eptr->p_elem)->synch_rad = 0;
    track_through_crbend(coord, n_track, (CRBEND*)eptr->p_elem, run->p_central, NULL, 0.0,
                         NULL, NULL, NULL, NULL, -1, -1);
    ((CRBEND*)eptr->p_elem)->isr = ltmp1;
    ((CRBEND*)eptr->p_elem)->synch_rad = ltmp2;
   break;
  case T_CWIGGLER:
    ltmp1 = ((CWIGGLER*)eptr->p_elem)->isr;
    ltmp2 = ((CWIGGLER*)eptr->p_elem)->sr;
    ((CWIGGLER*)eptr->p_elem)->isr = ((CWIGGLER*)eptr->p_elem)->sr = 0;
    GWigSymplecticPass(coord, n_track, run->p_central, (CWIGGLER*)eptr->p_elem, NULL, 0, NULL);
    ((CWIGGLER*)eptr->p_elem)->isr = ltmp1;
    ((CWIGGLER*)eptr->p_elem)->sr = ltmp2;
    break;
  case T_APPLE:
    ltmp1 = ((APPLE*)eptr->p_elem)->isr;
    ltmp2 = ((APPLE*)eptr->p_elem)->sr;
    ((APPLE*)eptr->p_elem)->isr = ((APPLE*)eptr->p_elem)->sr = 0;
    APPLE_Track(coord, n_track, run->p_central, (APPLE*)eptr->p_elem);
    ((APPLE*)eptr->p_elem)->isr = ltmp1;
    ((APPLE*)eptr->p_elem)->sr = ltmp2;
    break;
  case T_UKICKMAP:
    ltmp1 = ((UKICKMAP*)eptr->p_elem)->isr;
    ltmp2 = ((UKICKMAP*)eptr->p_elem)->synchRad;
    ((UKICKMAP*)eptr->p_elem)->isr = ((UKICKMAP*)eptr->p_elem)->synchRad = 0;
    if (trackUndulatorKickMap(coord, NULL, n_track, run->p_central, (UKICKMAP*)eptr->p_elem, 0)!=n_track) {
      printf("*** Error: particles lost in determineMatrix call for UKICKMAP\n");
      exitElegant(1);
    }
    ((UKICKMAP*)eptr->p_elem)->isr = ltmp1;
    ((UKICKMAP*)eptr->p_elem)->synchRad = ltmp2;
    break;
  case T_BGGEXP:
    ltmp1 = ((BGGEXP*)eptr->p_elem)->isr;
    ltmp2 = ((BGGEXP*)eptr->p_elem)->synchRad;
    ((BGGEXP*)eptr->p_elem)->isr = ((BGGEXP*)eptr->p_elem)->synchRad = 0;
    trackBGGExpansion(coord, n_track, (BGGEXP*)eptr->p_elem, run->p_central, NULL, NULL);
    ((BGGEXP*)eptr->p_elem)->isr = ltmp1;
    ((BGGEXP*)eptr->p_elem)->synchRad = ltmp2;
    break;
  case T_SCRIPT:
#if USE_MPI
    if (myid==0) {
#if MPI_DEBUG
      printf("Calling transformBeamWithScript with force serial=1 \n");
      fflush(stdout);
#endif
      transformBeamWithScript((SCRIPT*)eptr->p_elem, run->p_central, NULL, NULL, coord, n_track, 
                              NULL, 0, 2, -1.0, 1, eptr->occurence);
    } else  {
#if MPI_DEBUG
      printf("myid=%d in determineMatrix for SCRIPT\n", myid);
      fflush(stdout);
#endif
    }
#if MPI_DEBUG
    printf("Preparing MPI_Bcast coordinates from serial tracking\n");
    fflush(stdout);
#endif
    MPI_Bcast(&(coord[0][0]), n_track*COORDINATES_PER_PARTICLE, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#if MPI_DEBUG
    printf("Finished MPI_Bcast coordinates from serial tracking\n");
    fflush(stdout);
    if (1) {
      for (i=0; i<n_track; i++) 
        printf("%le %le %le %le %le %le\n",
               coord[i][0], coord[i][1], coord[i][2], coord[i][3], coord[i][4], coord[i][5]);
      fflush(stdout);
    }
#endif
#else
    /* Serial version */
    transformBeamWithScript((SCRIPT*)eptr->p_elem, run->p_central, NULL, NULL, coord, n_track, 
                            NULL, 0, 2, -1, 0, eptr->occurence);
#endif
    break;
  case T_TWMTA:
  case T_MAPSOLENOID:
  case T_TWLA:
    motion(coord, n_track, eptr->p_elem, eptr->type, &run->p_central, &dgamma, dP, NULL, 0.0);
    break;
  case T_RFDF:
    /* Don't actually use this */
    track_through_rf_deflector(coord, (RFDF*)eptr->p_elem,
			       coord, n_track, run->p_central, 0, eptr->end_pos, 0);
    break;
  case T_LSRMDLTR:
    ltmp1 = ((LSRMDLTR*)eptr->p_elem)->isr;
    ltmp2 = ((LSRMDLTR*)eptr->p_elem)->synchRad;
    dtmp1 = ((LSRMDLTR*)eptr->p_elem)->laserPeakPower;
    ((LSRMDLTR*)eptr->p_elem)->isr = ((LSRMDLTR*)eptr->p_elem)->synchRad = 0;
    ((LSRMDLTR*)eptr->p_elem)->laserPeakPower = 0;
    motion(coord, n_track, eptr->p_elem, eptr->type, &run->p_central, &dgamma, dP, NULL, 0.0);
    ((LSRMDLTR*)eptr->p_elem)->isr = ltmp1;
    ((LSRMDLTR*)eptr->p_elem)->synchRad = ltmp2;
    ((LSRMDLTR*)eptr->p_elem)->laserPeakPower = dtmp1;
    break;
  case T_FTABLE:
    field_table_tracking(coord, n_track, (FTABLE*)eptr->p_elem, run->p_central, run);
    break;
  default:
    printf("*** Error: determineMatrix called for element that is not supported!\n");
    printf("***        Seek professional help!\n");
    exitElegant(1);
    break;
  }
  
  M = tmalloc(sizeof(*M));
  M->order = 1;
  initialize_matrices(M, M->order);
  R = M->R;
  C = M->C;
  
  for (i=0; i<6; i++) {
    /* i indexes the dependent quantity */

    /* Determine C[i] */
    C[i] = coord[n_track-1][i];

    /* Compute R[i][j] */
    for (j=0; j<6; j++) {
      /* j indexes the initial coordinate value */
      /*
      if (n_points==2) 
        R[i][j] = (coord[2*j][i]-coord[2*j+1][i])/(2*stepSize[j]);
      else
      */
        R[i][j] = 
          (27*(coord[2*j][i]-coord[2*j+1][i])-(coord[2*j+12][i]-coord[2*j+13][i]))/(48*stepSize[j]);
    }
  }

  free_czarray_2d((void**)coord, 1+4*6, COORDINATES_PER_PARTICLE);

  /*
  if (strlen(eptr->name)<900)
    sprintf(s, "\nElement %s#%ld matrix determined from tracking:\n", eptr->name, eptr->occurence);
  else
    sprintf(s, "\nElement matrix determined from tracking:\n");
  print_matrices(stdout, s, M);
  */
#if USE_MPI
  notSinglePart = notSinglePart_saved;
#endif
  return M;
}

VMATRIX *determineMatrixHigherOrder(RUN *run, ELEMENT_LIST *eptr, double *startingCoord, double *stepSize, long order)
{
  double **initialCoord, **finalCoord, **coordError;
  long n_track, i, j, n_left;
  VMATRIX *M;
  double defaultStep[6] = {1e-4, 1e-4, 1e-4, 1e-4, 1e-4, 1e-4};
  double maximumValue[6];
  long ltmp1, ltmp2;
  double dgamma, dtmp1, dP[3];
  long nPoints1 = 5;
  long maxFitOrder = 4;
  /* We'll store some of the matrices to avoid recomputing them */
#define MAX_N_STORED_MATRICES 100
  static long nStoredMatrices = 0, iStoredMatrices = 0;
  static ELEMENT_LIST **storedElement=NULL;
  static VMATRIX **storedMatrix=NULL;

  /*
  if (storedElement==NULL) {
    storedElement = tmalloc(sizeof(*storedElement)*MAX_N_STORED_MATRICES);
    storedMatrix = tmalloc(sizeof(*storedMatrix)*MAX_N_STORED_MATRICES);
  }

  for (i=0; i<nStoredMatrices; i++) {
    CRBEND *crbptr0, *crbptr1;
    if (eptr->type==storedElement[i]->type && 
	memcmp(storedElement[i]->p_elem, eptr->p_elem, entity_description[eptr->type].user_structure_size)==0) {
      M = tmalloc(sizeof(*M));
      copy_matrices(M, storedMatrix[i]);
      switch (eptr->type) {
      case T_CRBEND:
	crbptr0 = (CRBEND*)storedElement[i]->p_elem;
	crbptr1 = (CRBEND*)eptr->p_elem;
	crbptr1->optimized = crbptr0->optimized;
	crbptr1->fseOffset = crbptr0->fseOffset;
	crbptr1->dxOffset = crbptr0->dxOffset;
	crbptr1->KnDelta = crbptr0->KnDelta;
	crbptr1->referenceData[0] = crbptr0->referenceData[0];
	crbptr1->referenceData[1] = crbptr0->referenceData[1];
	crbptr1->referenceData[2] = crbptr0->referenceData[2];
	crbptr1->referenceData[3] = crbptr0->referenceData[3];
	break;
      default:
	break;
      }
      return M;
    }
  }
  */

  if (stepSize==NULL)
    stepSize = defaultStep;

  n_track = makeInitialParticleEnsemble(&initialCoord, startingCoord, &finalCoord, &coordError, nPoints1, stepSize);
  n_left = n_track;

  switch (eptr->type) {
  case T_CSBEND:
    ltmp1 = ((CSBEND*)eptr->p_elem)->isr;
    ltmp2 = ((CSBEND*)eptr->p_elem)->synch_rad;
    ((CSBEND*)eptr->p_elem)->isr = ((CSBEND*)eptr->p_elem)->synch_rad = 0;
    n_left = track_through_csbend(finalCoord, n_track, (CSBEND*)eptr->p_elem, 0.0, run->p_central, NULL, 0.0,
                         NULL, NULL, NULL, NULL);
    ((CSBEND*)eptr->p_elem)->isr = ltmp1;
    ((CSBEND*)eptr->p_elem)->synch_rad = ltmp2;
   break;
  case T_CRBEND:
    ltmp1 = ((CRBEND*)eptr->p_elem)->isr;
    ltmp2 = ((CRBEND*)eptr->p_elem)->synch_rad;
    ((CRBEND*)eptr->p_elem)->isr = ((CRBEND*)eptr->p_elem)->synch_rad = 0;
    n_left = track_through_crbend(finalCoord, n_track, (CRBEND*)eptr->p_elem, run->p_central, NULL, 0.0,
                                  NULL, NULL, NULL, NULL, -1, -1);
    ((CRBEND*)eptr->p_elem)->isr = ltmp1;
    ((CRBEND*)eptr->p_elem)->synch_rad = ltmp2;
   break;
  case T_CWIGGLER:
    ltmp1 = ((CWIGGLER*)eptr->p_elem)->isr;
    ltmp2 = ((CWIGGLER*)eptr->p_elem)->sr;
    ((CWIGGLER*)eptr->p_elem)->isr = ((CWIGGLER*)eptr->p_elem)->sr = 0;
    GWigSymplecticPass(finalCoord, n_track, run->p_central, (CWIGGLER*)eptr->p_elem, NULL, 0, NULL);
    ((CWIGGLER*)eptr->p_elem)->isr = ltmp1;
    ((CWIGGLER*)eptr->p_elem)->sr = ltmp2;
    break;
  case T_APPLE:
    ltmp1 = ((APPLE*)eptr->p_elem)->isr;
    ltmp2 = ((APPLE*)eptr->p_elem)->sr;
    ((APPLE*)eptr->p_elem)->isr = ((APPLE*)eptr->p_elem)->sr = 0;
    APPLE_Track(finalCoord, n_track, run->p_central, (APPLE*)eptr->p_elem);
    ((APPLE*)eptr->p_elem)->isr = ltmp1;
    ((APPLE*)eptr->p_elem)->sr = ltmp2;
    break;
  case T_UKICKMAP:
    ltmp1 = ((UKICKMAP*)eptr->p_elem)->isr;
    ltmp2 = ((UKICKMAP*)eptr->p_elem)->synchRad;
    ((UKICKMAP*)eptr->p_elem)->isr = ((UKICKMAP*)eptr->p_elem)->synchRad = 0;
    if (trackUndulatorKickMap(finalCoord, NULL, n_track, run->p_central, (UKICKMAP*)eptr->p_elem, 0)!=n_track) {
      printf("*** Error: particles lost in determineMatrix call for UKICKMAP\n");
      exitElegant(1);
    }
    ((UKICKMAP*)eptr->p_elem)->isr = ltmp1;
    ((UKICKMAP*)eptr->p_elem)->synchRad = ltmp2;
    break;
  case T_BGGEXP:
    ltmp1 = ((BGGEXP*)eptr->p_elem)->isr;
    ltmp2 = ((BGGEXP*)eptr->p_elem)->synchRad;
    ((BGGEXP*)eptr->p_elem)->isr = ((BGGEXP*)eptr->p_elem)->synchRad = 0;
    trackBGGExpansion(finalCoord, n_track, (BGGEXP*)eptr->p_elem, run->p_central, NULL, NULL);
    ((BGGEXP*)eptr->p_elem)->isr = ltmp1;
    ((BGGEXP*)eptr->p_elem)->synchRad = ltmp2;
    break;
  case T_TWMTA:
  case T_MAPSOLENOID:
  case T_TWLA:
    motion(finalCoord, n_track, eptr->p_elem, eptr->type, &run->p_central, &dgamma, dP, NULL, 0.0);
    break;
  case T_RFDF:
    /* Don't actually use this */
    track_through_rf_deflector(finalCoord, (RFDF*)eptr->p_elem,
			       finalCoord, n_track, run->p_central, 0, eptr->end_pos, 0);
    break;
  case T_LSRMDLTR:
    ltmp1 = ((LSRMDLTR*)eptr->p_elem)->isr;
    ltmp2 = ((LSRMDLTR*)eptr->p_elem)->synchRad;
    dtmp1 = ((LSRMDLTR*)eptr->p_elem)->laserPeakPower;
    ((LSRMDLTR*)eptr->p_elem)->isr = ((LSRMDLTR*)eptr->p_elem)->synchRad = 0;
    ((LSRMDLTR*)eptr->p_elem)->laserPeakPower = 0;
    motion(finalCoord, n_track, eptr->p_elem, eptr->type, &run->p_central, &dgamma, dP, NULL, 0.0);
    ((LSRMDLTR*)eptr->p_elem)->isr = ltmp1;
    ((LSRMDLTR*)eptr->p_elem)->synchRad = ltmp2;
    ((LSRMDLTR*)eptr->p_elem)->laserPeakPower = dtmp1;
    break;
  case T_FTABLE:
    field_table_tracking(finalCoord, n_track, (FTABLE*)eptr->p_elem, run->p_central, run);
    break;
  default:
    printf("*** Error: determineMatrixHigherOrder called for element that is not supported!\n");
    printf("***        Seek professional help!\n");
    exitElegant(1);
    break;
  }
  if (n_left!=n_track) {
    bombElegantVA("lost particles when tracking to determine matrix for %s\n", eptr->name);
  }

  /*
  if (1) {
    FILE *fplog;
    fplog = fopen("matrix.log", "w");
    fprintf(fplog, "SDDS1\n");
    fprintf(fplog, "&column name=x0 type=double &end\n");
    fprintf(fplog, "&column name=x1 type=double &end\n");
    fprintf(fplog, "&column name=xp0 type=double &end\n");
    fprintf(fplog, "&column name=xp1 type=double &end\n");
    fprintf(fplog, "&column name=y0 type=double &end\n");
    fprintf(fplog, "&column name=y1 type=double &end\n");
    fprintf(fplog, "&column name=yp0 type=double &end\n");
    fprintf(fplog, "&column name=yp1 type=double &end\n");
    fprintf(fplog, "&column name=s0 type=double &end\n");
    fprintf(fplog, "&column name=s1 type=double &end\n");
    fprintf(fplog, "&column name=delta0 type=double &end\n");
    fprintf(fplog, "&column name=delta1 type=double &end\n");
    fprintf(fplog, "&data mode=ascii no_row_counts=1 &end\n");
    for (i=0; i<n_track; i++) {
      for (j=0; j<6; j++)
        fprintf(fplog, "%13.6le %13.6le ", initialCoord[i][j], finalCoord[i][j]);
      fprintf(fplog, "\n");
    }
    fclose(fplog);
  }
  */

  /* Set errors as fraction of the absolute maximum coordinate */
  for (i=0; i<6; i++) {
    double min, max;
    max = -(min = DBL_MAX);
    maximumValue[i] = -DBL_MAX;
    for (j=0; j<n_track; j++) {
      if (max<finalCoord[j][i])
        max = finalCoord[j][i];
      if (min>finalCoord[j][i])
        min = finalCoord[j][i];
      if (maximumValue[i]<initialCoord[j][i])
        maximumValue[i] = initialCoord[j][i];
    }
    max = fabs(max);
    min = fabs(min);
    max = max>min ? max : min;
    for (j=0; j<n_track; j++)
      coordError[j][i] = max*accuracy_factor;
  }
  
  M = computeMatricesFromTracking(stdout, initialCoord, finalCoord, coordError, stepSize,
                                  maximumValue, nPoints1, n_track, maxFitOrder, 0);

  /*
  storedElement[iStoredMatrices] = eptr;
  storedMatrix[iStoredMatrices] = M;
  if (nStoredMatrices<MAX_N_STORED_MATRICES) {
    iStoredMatrices++;
    nStoredMatrices++;
  } else {
    iStoredMatrices++;
    if (iStoredMatrices==nStoredMatrices)
      iStoredMatrices = 0;
  }
  */

  /*
  for (i=0; i<6; i++) {
    printf("R[%ld]: ", i);
    for (j=0; j<6; j++) {
      printf("%le ", M->R[i][j]);
    }
    printf("\n");
  }
  */

  free_czarray_2d((void**)initialCoord, n_track, COORDINATES_PER_PARTICLE);
  free_czarray_2d((void**)finalCoord, n_track, COORDINATES_PER_PARTICLE);
  free_czarray_2d((void**)coordError, n_track, COORDINATES_PER_PARTICLE);

  free_matrices_above_order(M, order);

  /*
  for (i=0; i<6; i++) {
    if (fabs(M->C[i])<1e-10)
      M->C[i] = 0;
    for (j=0; j<6; j++) {
      long k;
      if (fabs(M->R[i][j])<1e-8)
        M->R[i][j] = 0;
      for (k=0; k<=j; k++) 
        if (fabs(M->T[i][j][k])<1e-6)
          M->T[i][j][k] = 0;
    }
  }
  */

  /* print_matrices(stdout, "Matrix from fitting:", M); */

  return M;
}

/* FILE *fpdeb = NULL; */

void determineRadiationMatrix(VMATRIX *Mr, RUN *run, ELEMENT_LIST *eptr, double *startingCoord, double *Dr, long nSlices, long sliceEtilted, long order)
{
  CSBEND csbend; CSRCSBEND *csrcsbend; BEND *sbend; WIGGLER *wig; CRBEND crbend;
  KQUAD kquad;  QUAD *quad; CWIGGLER cwig; BGGEXP bggexp;
  KSEXT ksext; SEXT *sext;
  HCOR hcor; VCOR vcor; HVCOR hvcor;
  EHCOR ehcor; EVCOR evcor; EHVCOR ehvcor;
  double length, z;
  long i, j, k, slice, nSlices0;
  double *accumD1, *accumD2, *dtmp;
  double post_xkick, post_ykick;
  VMATRIX *M1, *M2, *Ml1, *Mtmp;
  ELEMENT_LIST elem;
  MATRIX *Ms;
  long ignoreRadiation = 0;
  
/*
  fpdeb = fopen("analyze.sdds", "w");
  fprintf(fpdeb, "SDDS1\n&column name=z type=double &end\n");
  fprintf(fpdeb, "&column name=x type=double &end\n&column name=xp type=double &end\n");
  fprintf(fpdeb, "&column name=y type=double &end\n&column name=yp type=double &end\n");
  fprintf(fpdeb, "&column name=s type=double &end\n&column name=p type=double &end\n");
  fprintf(fpdeb, "&data mode=ascii no_row_counts=1 &end\n");
*/

  /* Accumulated diffusion matrix */
  accumD1 = tmalloc(21*sizeof(*(accumD1)));
  memset(accumD1, 0, 21*sizeof(*(accumD1)));
  accumD2 = tmalloc(21*sizeof(*(accumD2)));
  memset(accumD2, 0, 21*sizeof(*(accumD2)));

  /* Matrices for C, R matrix propagation: */
  initialize_matrices(M1=tmalloc(sizeof(*M1)), order);
  initialize_matrices(M2=tmalloc(sizeof(*M2)), order);
  /* Temporary variable for linear matrix with radiation: */
  initialize_matrices(Ml1=tmalloc(sizeof(*Ml1)), 1);
  for (i=0; i<6; i++) {
    M1->R[i][i] = 1;
    M1->C[i] = startingCoord[i];
  }
  /* Matrix for sigma matrix (D matrix) propagation */
  m_alloc(&Ms, 21, 21);

  nSlices0 = nSlices;

  /* NB: CWIGGLER will be used in single-step mode, which is why this slicing ends up
   * sub-dividing the periods 
   */
  if (eptr->type==T_CWIGGLER)
    nSlices = ((CWIGGLER*)eptr->p_elem)->periods*((CWIGGLER*)eptr->p_elem)->stepsPerPeriod;
  else if (eptr->type==T_WIGGLER)
    nSlices *= (((WIGGLER*)eptr->p_elem)->poles/2);
  else if (eptr->type==T_CRBEND) {
    if (nSlices > ((CRBEND*)eptr->p_elem)->n_kicks)
      nSlices = ((CRBEND*)eptr->p_elem)->n_kicks;
  }
  z = 0;
  
  elem.end_pos = eptr->end_pos;
  elem.name = NULL;
  elem.occurence = 0;
  elem.type = eptr->type;
  length = 0;
  for (slice=0; slice<nSlices; slice++) {
    post_xkick = post_ykick = 0; /* use this to handle pre- and post-KQUAD kicks */
    switch (eptr->type) {
    case T_CSBEND:
      memcpy(&csbend, (CSBEND*)eptr->p_elem, sizeof(CSBEND));
      if (csbend.etilt && !sliceEtilted) {
        nSlices = 1;
      }
      csbend.isr = 0;
      csbend.angle /= nSlices;
      length = (csbend.length /= nSlices);
      csbend.n_kicks = fabs(csbend.angle/0.005) + 1;
      csbend.refTrajectoryChangeSet = 0;
      csbend.refLength = 0;
      csbend.refAngle = 0;
      csbend.refTrajectoryChange = NULL;
      csbend.referenceCorrection = 0;
      if (slice!=0)
        csbend.edgeFlags &= ~BEND_EDGE1_EFFECTS;
      if (slice!=nSlices-1) 
        csbend.edgeFlags &= ~BEND_EDGE2_EFFECTS;
      elem.type = T_CSBEND;
      elem.p_elem = (void*)&csbend;
      break;
    case T_CRBEND:
      if (slice==0) {
        memcpy(&crbend, (CRBEND*)eptr->p_elem, sizeof(CRBEND));
        crbend.isr = 0;
        crbend.n_kicks = nSlices;
        elem.type = T_CRBEND;
        elem.p_elem = (void*)&crbend;
      }
      break;
    case T_SBEN:
      if (slice!=0)
        csbend.edge_effects[0] = 0;
      if (slice!=nSlices-1)
        csbend.edge_effects[0] = 0;
      elem.type = T_CSBEND;
      elem.p_elem = (void*)&csbend;
      sbend = (BEND*)eptr->p_elem;
      memset(&csbend, 0, sizeof(csbend));
      csbend.isr = 0;
      csbend.synch_rad = 1;
      if (sbend->etilt && !sliceEtilted) {
        nSlices = 1;
      }
      length = csbend.length = sbend->length/nSlices;
      csbend.angle = sbend->angle/nSlices;
      csbend.k1 = sbend->k1;
      csbend.e[0] = sbend->e[0];
      csbend.e[1] = sbend->e[1];
      csbend.k2 = sbend->k2;
      csbend.h[0] = sbend->h[0];
      csbend.h[1] = sbend->h[1];
      csbend.hgap = sbend->hgap;
      csbend.fint = sbend->fint;
      csbend.dx = sbend->dx;
      csbend.dy = sbend->dy;
      csbend.dz = sbend->dz;
      csbend.fse = sbend->fse;
      csbend.tilt = sbend->tilt;
      csbend.etilt = sbend->etilt;
      csbend.edge_effects[0] = sbend->edge_effects[0];
      csbend.edge_effects[1] = sbend->edge_effects[1];
      csbend.edge_order = sbend->edge_order;
      csbend.edgeFlags = sbend->edgeFlags;
      csbend.refTrajectoryChangeSet = 0;
      csbend.refLength = 0;
      csbend.refAngle = 0;
      csbend.refTrajectoryChange = NULL;
      csbend.e1Index = sbend->e1Index;
      csbend.e2Index = sbend->e2Index;
      if (slice!=0)
        csbend.edgeFlags &= ~BEND_EDGE1_EFFECTS;
      if (slice!=nSlices-1)
        csbend.edgeFlags &= ~BEND_EDGE2_EFFECTS;
      csbend.k1 = sbend->k1;
      csbend.k2 = sbend->k2;
      csbend.n_kicks = fabs(csbend.angle/0.005) + 1;
      csbend.integration_order = 4;
      break;
    case T_CSRCSBEND:
      if (slice!=0)
        csbend.edge_effects[0] = 0;
      if (slice!=nSlices-1)
        csbend.edge_effects[1] = 0;
      elem.type = T_CSBEND;
      elem.p_elem = (void*)&csbend;
      csrcsbend = (CSRCSBEND*)eptr->p_elem;
      if (csrcsbend->etilt && !sliceEtilted) {
        nSlices = 1;
      }
      memset(&csbend, 0, sizeof(csbend));
      csbend.isr = 0;
      csbend.synch_rad = 1;
      length = csbend.length = csrcsbend->length/nSlices;
      csbend.angle = csrcsbend->angle/nSlices;
      csbend.k1 = csrcsbend->k1;
      csbend.e[0] = csrcsbend->e[0];
      csbend.e[1] = csrcsbend->e[1];
      csbend.k2 = csrcsbend->k2;
      csbend.h[0] = csrcsbend->h[0];
      csbend.h[1] = csrcsbend->h[1];
      csbend.hgap = csrcsbend->hgap;
      csbend.fint = csrcsbend->fint;
      csbend.dx = csrcsbend->dx;
      csbend.dy = csrcsbend->dy;
      csbend.dz = csrcsbend->dz;
      csbend.fse = csrcsbend->fse;
      csbend.tilt = csrcsbend->tilt;
      csbend.etilt = csrcsbend->etilt;
      csbend.edge_effects[0] = csrcsbend->edge_effects[0];
      csbend.edge_effects[1] = csrcsbend->edge_effects[1];
      csbend.edge_order = csrcsbend->edge_order;
      csbend.edgeFlags = csrcsbend->edgeFlags;
      csbend.refTrajectoryChangeSet = 0;
      csbend.refLength = 0;
      csbend.refAngle = 0;
      csbend.refTrajectoryChange = NULL;
      csbend.e1Index = csrcsbend->e1Index;
      csbend.e2Index = csrcsbend->e2Index;
      if (slice!=0)
        csbend.edgeFlags &= ~BEND_EDGE1_EFFECTS;
      if (slice!=nSlices-1)
        csbend.edgeFlags &= ~BEND_EDGE2_EFFECTS;
      csbend.k1 = csrcsbend->k1;
      csbend.k2 = csrcsbend->k2;
      csbend.n_kicks = fabs(csbend.angle/0.005) + 1;
      csbend.integration_order = 4;
      break;
    case T_KQUAD:
      memcpy(&kquad, (KQUAD*)eptr->p_elem, sizeof(KQUAD));
      kquad.isr = 0;
      if (slice==0) {
        M1->C[1] += kquad.xkick/2;
        M1->C[3] += kquad.ykick/2;
      }
      if (slice==nSlices-1) {
        post_xkick = kquad.xkick/2;
        post_ykick = kquad.ykick/2;
      }
      kquad.xkick = kquad.ykick = 0;
      length = (kquad.length /= nSlices);
      kquad.n_kicks = 4 + (long)(fabs(kquad.k1)*sqr(kquad.length));
      if (slice!=0)
	kquad.edge1_effects = 0;
      if (slice!=nSlices-1)
	kquad.edge2_effects = 0;
      elem.type = T_KQUAD;
      elem.p_elem = (void*)&kquad;
      break;
    case T_QUAD:
      quad = (QUAD*)eptr->p_elem;
      memset(&kquad, 0, sizeof(KQUAD));
      kquad.isr = 0;
      kquad.synch_rad = 1;
      length = (kquad.length = quad->length/nSlices);
      kquad.k1 = quad->k1;
      kquad.tilt = quad->tilt;
      if (quad->ffringe)
        bombElegant("Can't perform radiation matrix calculations when QUAD has nonzero FFRINGE parameter", NULL);
      kquad.dx = quad->dx;
      kquad.dy = quad->dy;
      kquad.dz = quad->dz;
      kquad.fse = quad->fse;
      kquad.xkick = quad->xkick;
      kquad.ykick = quad->ykick;
      kquad.xKickCalibration = quad->xKickCalibration;
      kquad.yKickCalibration = quad->yKickCalibration;
      kquad.n_kicks = 4 + (long)(fabs(kquad.k1)*sqr(kquad.length));
      kquad.integration_order = 4;
      if (slice!=0)
	kquad.edge1_effects = 0;
      if (slice!=nSlices-1)
	kquad.edge2_effects = 0;
      elem.type = T_KQUAD;
      elem.p_elem = (void*)&kquad;
      break;
    case T_KSEXT:
      memcpy(&ksext, (KSEXT*)eptr->p_elem, sizeof(KSEXT));
      ksext.isr = 0;
      ksext.n_kicks = 4;
      length = (ksext.length /= nSlices);
      if (length<1e-6) {
	ignoreRadiation = 1;
	ksext.synch_rad = 0;
      }
      elem.type = T_KSEXT;
      elem.p_elem = (void*)&ksext;
      break;
    case T_SEXT:
      sext = (SEXT*)eptr->p_elem;
      memset(&ksext, 0, sizeof(KSEXT));
      ksext.isr = 0;
      ksext.synch_rad = 1;
      length = (ksext.length = sext->length/nSlices);
      if (length<1e-6) {
	ignoreRadiation = 1;
	ksext.synch_rad = 0;
      }
      ksext.k2 = sext->k2;
      ksext.tilt = sext->tilt;
      ksext.dx = sext->dx;
      ksext.dy = sext->dy;
      ksext.dz = sext->dz;
      ksext.fse = sext->fse;
      ksext.n_kicks = 4;
      ksext.integration_order = 4;
      elem.type = T_KSEXT;
      elem.p_elem = (void*)&ksext;
      break;
    case T_RFCA:
      nSlices = 1;
      elem.type = T_RFCA;
      elem.p_elem = (void*)eptr->p_elem;
      length = ((RFCA*)eptr->p_elem)->length;
      break;
    case T_TWLA:
      nSlices = 1;
      elem.type = T_TWLA;
      elem.p_elem = (void*)eptr->p_elem;
      length = ((TW_LINAC*)eptr->p_elem)->length/nSlices;
      break;
    case T_HCOR:
      memcpy(&elem, eptr, sizeof(elem));
      elem.matrix = NULL;
      elem.p_elem = (void*)&hcor;
      memcpy(&hcor, eptr->p_elem, sizeof(hcor));
      length = (hcor.length /= nSlices);
      hcor.lEffRad /= nSlices;
      hcor.kick /= nSlices;
      hcor.isr = 0;
      break;
    case T_VCOR:
      memcpy(&elem, eptr, sizeof(elem));
      elem.matrix = NULL;
      elem.p_elem = (void*)&vcor;
      memcpy(&vcor, eptr->p_elem, sizeof(vcor));
      length = (vcor.length /= nSlices);
      vcor.lEffRad /= nSlices;
      vcor.kick /= nSlices;
      vcor.isr = 0;
      break;
    case T_HVCOR:
      memcpy(&elem, eptr, sizeof(elem));
      elem.matrix = NULL;
      elem.p_elem = (void*)&hvcor;
      memcpy(&hvcor, eptr->p_elem, sizeof(hvcor));
      length = (hvcor.length /= nSlices);
      hvcor.lEffRad /= nSlices;
      hvcor.xkick /= nSlices;
      hvcor.ykick /= nSlices;
      hvcor.isr = 0;
      break;
    case T_EHCOR:
      memcpy(&elem, eptr, sizeof(elem));
      elem.matrix = NULL;
      elem.p_elem = (void*)&ehcor;
      memcpy(&ehcor, eptr->p_elem, sizeof(ehcor));
      length = (ehcor.length /= nSlices);
      ehcor.lEffRad /= nSlices;
      ehcor.kick /= nSlices;
      ehcor.isr = 0;
      break;
    case T_EVCOR:
      memcpy(&elem, eptr, sizeof(elem));
      elem.matrix = NULL;
      elem.p_elem = (void*)&evcor;
      memcpy(&evcor, eptr->p_elem, sizeof(evcor));
      length = (evcor.length /= nSlices);
      evcor.lEffRad /= nSlices;
      evcor.kick /= nSlices;
      evcor.isr = 0;
      break;
    case T_EHVCOR:
      memcpy(&elem, eptr, sizeof(elem));
      elem.matrix = NULL;
      elem.p_elem = (void*)&ehvcor;
      memcpy(&ehvcor, eptr->p_elem, sizeof(ehvcor));
      length = (ehvcor.length /= nSlices);
      ehvcor.lEffRad /= nSlices;
      ehvcor.xkick /= nSlices;
      ehvcor.ykick /= nSlices;
      ehvcor.isr = 0;
      break;
    case T_CWIGGLER:
      if (slice==0) {
        memcpy(&cwig, (CWIGGLER*)eptr->p_elem, sizeof(CWIGGLER));
        cwig.isr = 0;
        elem.type = T_CWIGGLER;
        elem.p_elem = (void*)&cwig;
        length = cwig.length/nSlices;
      }
      break;
    case T_WIGGLER:
      if (slice==0) {
        wig = (WIGGLER*)eptr->p_elem;
        memset(&cwig, 0, sizeof(cwig));
        cwig.sinusoidal = 1;
        cwig.length = wig->length;
        cwig.tilt = wig->tilt;
        cwig.dx = wig->dx;
        cwig.dy = wig->dy;
        cwig.dz = wig->dz;
        cwig.periods = wig->poles/2;
        cwig.stepsPerPeriod = nSlices0;
        cwig.sr = 1;
        cwig.integrationOrder = 4;
        if (wig->K) {
          double lambda;
          lambda = cwig.length/cwig.periods;
          cwig.ByMax = wig->K/(UNDULATOR_K_FACTOR*lambda);
        } else if (wig->radius) {
          double H;
          H = run->p_central*RIGIDITY_FACTOR;
          cwig.ByMax = H/wig->radius;
        } else
        cwig.ByMax = wig->B;
        elem.type = T_CWIGGLER;
        elem.p_elem = (void*)&cwig;
        length = cwig.length/nSlices;
      }
      break;
    case T_BGGEXP:
      if (slice==0) {
        nSlices = 1;
        memcpy(&bggexp, (BGGEXP*)eptr->p_elem, sizeof(BGGEXP));
        bggexp.isr = 0;
        elem.type = T_BGGEXP;
        elem.p_elem = (void*)&bggexp;
      }
      break;
    default:
      printf("*** Error: determineRadiationMatrix called for element (%s) that is not supported!\n", eptr->name);
      printf("***        Seek professional help!\n");
      exit(1);
      break;
    }

    /* Step 1: determine effective R matrix for this element, as well as the diffusion matrix */
    determineRadiationMatrix1(Ml1, run, &elem, M1->C, accumD2, ignoreRadiation, &z, slice); 
    Ml1->C[1] += post_xkick;
    Ml1->C[3] += post_ykick;
    /* printf("z = %le ", z);
     print_matrices(stdout, "matrix1:", Ml1);*/

    /* Step 2: Propagate the diffusion matrix */
    fillSigmaPropagationMatrix(Ms->a, Ml1->R);
    for (i=0; i<21; i++) 
      for (j=0; j<21; j++)
        accumD2[i] += Ms->a[i][j]*accumD1[j];
    dtmp    = accumD1;
    accumD1 = accumD2;
    accumD2 = dtmp;
    
    /* Step 3: Copy the propagated C vector */
    memcpy(M2->C, Ml1->C, 6*sizeof(*(M2->C)));

    /* Step 4: Multiply the R matrices */
    for (i=0; i<6; i++)
      for (j=0; j<6; j++) {
        M2->R[i][j] = 0;
        for (k=0; k<6; k++) 
          M2->R[i][j] += Ml1->R[i][k]*M1->R[k][j];
      }

    Mtmp = M2;
    M2   = M1;
    M1   = Mtmp;

    elem.end_pos += length;
  }
  
  /* Copy the matrix to the caller's buffer */
  for (i=0; i<6; i++) {
    Mr->C[i] = M1->C[i];
    for (j=0; j<6; j++) {
      Mr->R[i][j] = M1->R[i][j];
    }
  }
  for (i=0; i<21; i++)
    Dr[i] = accumD1[i];

  free(accumD1);
  free(accumD2);
  free_matrices(M2); tfree(M2);
  free_matrices(M1); tfree(M1);
  free_matrices(Ml1); tfree(Ml1);
  M1 = M2 = Ml1 = NULL;
  m_free(&Ms);
}


void determineRadiationMatrix1(VMATRIX *Mr, RUN *run, ELEMENT_LIST *elem, double *startingCoord, double *D, long ignoreRadiation, double *z,
                               long iSlice)
{
  CSBEND *csbend;
  CRBEND *crbend;
  KQUAD *kquad;
  KSEXT *ksext;
  double **coord, pCentral;
  long n_track, i, j;
  double **R, *C, Cs0;
  double stepSize[6] = {1e-5, 1e-5, 1e-5, 1e-5, 1e-3, 1e-5};
  double sigmaDelta2;
  double dP[3], dgamma;
  VMATRIX *matrix;

  coord = (double**)czarray_2d(sizeof(**coord), 1+6*2, COORDINATES_PER_PARTICLE);

  Cs0 = startingCoord ? startingCoord[4]: 0;
  n_track = 2*6+1;
  for (j=0; j<6; j++)
    for (i=0; i<n_track; i++)
      coord[i][j] = (startingCoord && j!=4) ? startingCoord[j] : 0;

  /* particles 0 and 1 are for d/dx, 2 and 3 are for d/dxp, etc */
  /* particle n_track-1 is the reference particle (coordinates set above) */
  for (i=0; i<6; i++) {
    coord[2*i+0][i] += stepSize[i] ;
    coord[2*i+1][i] -= stepSize[i] ;
  }
  
  if (D) {
    for (i=0; i<21; i++)
        D[i] = 0;
  }
  sigmaDelta2 = 0;
  setTrackingContext(elem->name, elem->occurence, elem->type, run->rootname);
  switch (elem->type) {
  case T_CSBEND:
    csbend = (CSBEND*)elem->p_elem;
    track_through_csbend(coord, n_track, csbend, 0, run->p_central, NULL, elem->end_pos-csbend->length, &sigmaDelta2, run->rootname, NULL, NULL);
    break;
  case T_CRBEND:
    crbend = (CRBEND*)elem->p_elem;
    track_through_crbend(coord, n_track, crbend, run->p_central, NULL, elem->end_pos-crbend->length, &sigmaDelta2, 
                         run->rootname, NULL, NULL, iSlice, -1);
    break;
  case T_SBEN:
    track_particles(coord, elem->matrix, coord, n_track);
    break;
  case T_KQUAD:
    kquad = (KQUAD*)elem->p_elem;
    multipole_tracking2(coord, n_track, elem, 0.0, run->p_central, NULL, elem->end_pos-kquad->length,
                        NULL, NULL, &sigmaDelta2);
    break;
  case T_KSEXT:
    ksext = (KSEXT*)elem->p_elem;
    multipole_tracking2(coord, n_track, elem, 0.0, run->p_central, NULL, elem->end_pos-ksext->length,
                        NULL, NULL, &sigmaDelta2);
    break;
  case T_RFCA:
    pCentral = run->p_central;
    for (i=0; i<n_track; i++)
      coord[i][4] += Cs0;
    Cs0 = 0;
    simple_rf_cavity(coord, n_track, (RFCA*)elem->p_elem, NULL, &pCentral, elem->end_pos);
    break;
  case T_HCOR:
  case T_VCOR:
  case T_HVCOR:
    matrix = compute_matrix(elem, run, NULL);
    track_particles(coord, matrix, coord, n_track);
    addCorrectorRadiationKick(coord, n_track, elem, elem->type, run->p_central, &sigmaDelta2, 1);
    free_matrices(matrix); free(matrix); matrix = NULL;
    break;
  case T_EHCOR:
  case T_EVCOR:
  case T_EHVCOR:
    trackThroughExactCorrector(coord, n_track, elem, run->p_central, NULL, 0, &sigmaDelta2);
    break;
  case T_TWLA:
    pCentral = run->p_central;
    motion(coord, n_track, elem->p_elem, elem->type, &pCentral, &dgamma, dP, NULL, 0.0);
    break;
  case T_CWIGGLER:
    GWigSymplecticPass(coord, n_track, run->p_central, (CWIGGLER*)elem->p_elem, &sigmaDelta2, 1, z);
    break;
  case T_BGGEXP:
    trackBGGExpansion(coord, n_track, (BGGEXP*)elem->p_elem, run->p_central, NULL, &sigmaDelta2);
    break;
  default:
    printf("*** Error: determineRadiationMatrix1 called for element (%s) that is not supported!\n", elem->name);
    printf("***        Seek professional help!\n");
    exitElegant(1);
    break;
  }
  if (ignoreRadiation)
    sigmaDelta2 = 0;

  /*  fprintf(fpdeb, "%le %le %le %le %le %le %le\n",
          *z,
          (coord[0][0]+coord[1][0])/2,
          (coord[0][1]+coord[1][1])/2,
          (coord[0][2]+coord[1][2])/2,
          (coord[0][3]+coord[1][3])/2,
          (coord[0][4]+coord[1][4])/2,
          (coord[0][5]+coord[1][5])/2);
          */

  R = Mr->R;
  C = Mr->C;
  for (i=0; i<6; i++) {
    /* i indexes the dependent quantity */
    
    /* Determine C[i] */
    C[i] = coord[n_track-1][i];
    
    /* Compute R[i][j] */
    for (j=0; j<6; j++) {
      /* j indexes the initial coordinate value */
      R[i][j] = (coord[2*j][i]-coord[2*j+1][i])/(2*stepSize[j]);
    }
  }

  D[sigmaIndex3[5][5]] = sigmaDelta2;
  
  C[4] += Cs0;
  free_czarray_2d((void**)coord, 1+2*6, COORDINATES_PER_PARTICLE);

}

void copyParticles(double ***coordCopy, double **coord, long np)
{
  long i;
  *coordCopy = (double**)zarray_2d(sizeof(***coordCopy), np, COORDINATES_PER_PARTICLE);
  for (i=0; i<np; i++)
    memcpy((*coordCopy)[i], coord[i], sizeof(double)*COORDINATES_PER_PARTICLE);
}

void performChromaticAnalysisFromMap(VMATRIX *M, TWISS *twiss, CHROM_DERIVS *chromDeriv)
{
  long i;
  double beta, alpha, cos_phi, sin_phi, phi, det, eta, etap;

  if (periodic) {
    for (i=0; i<4; i+=2) {
      /* find nu, beta, alpha */
      if (fabs(cos_phi = (M->R[i][i] + M->R[i+1][i+1])/2)>1) {
        printf("warning: beamline unstable for %c plane\n", i==0?'x':'y');
        fflush(stdout);
        beta = alpha = phi = 0;
      } else {
        beta = fabs(M->R[i][i+1]/sin(acos(cos_phi)));
        sin_phi = M->R[i][i+1]/beta;
        if ((phi = atan2(sin_phi, cos_phi))<0)
          phi += PIx2;
        alpha = (M->R[i][i]-M->R[i+1][i+1])/(2*sin_phi);
      }
      
      /* compute eta and eta' */
      if ((det = (2 - M->R[i][i] - M->R[i+1][i+1]))<=0) {
        printf("error: beamline unstable for %c plane--can't match dispersion functions\n",
		i==0?'x':'y');
        fflush(stdout);
        det = 1e-6;
      }
      eta  = ((1-M->R[i+1][i+1])*M->R[i][5]+M->R[i][i+1]*M->R[i+1][5])/det;
      etap = (M->R[i+1][i]*M->R[i][5] + (1-M->R[i][i])*M->R[i+1][5])/det;
      if (i==0) {
        twiss->betax = beta;
        twiss->alphax = alpha;
        twiss->phix = phi;
        twiss->etax = eta;
        twiss->etapx = etap;
      } else {
        twiss->betay = beta;
        twiss->alphay = alpha;
        twiss->phiy = phi;
        twiss->etay = eta;
        twiss->etapy = etap;
      }
    }
    
    twiss->periodic = 1;
    computeChromaticities(&(chromDeriv->tune1[0]), &(chromDeriv->tune1[1]),
                          &(chromDeriv->beta1[0]), &(chromDeriv->beta1[1]),
                          &(chromDeriv->alpha1[0]), &(chromDeriv->alpha1[1]),
                          twiss, twiss, M);
  } else {
    TWISS twiss0;
    twiss0.betax = beta_x;
    twiss0.alphax = alpha_x;
    twiss0.phix = 0;
    twiss0.etax = eta_x;
    twiss0.etapx = etap_x;
    twiss0.betay = beta_y;
    twiss0.alphay = alpha_y;
    twiss0.phiy = 0;
    twiss0.etay = eta_y;
    twiss0.etapy = etap_y;

    
    propagateTwissParameters(twiss, &twiss0, M);
    
    twiss->periodic = 0;
    
    computeChromaticities(&(chromDeriv->tune1[0]), &(chromDeriv->tune1[1]),
                          &(chromDeriv->beta1[0]), &(chromDeriv->beta1[1]),
                          &(chromDeriv->alpha1[0]), &(chromDeriv->alpha1[1]),
                          &twiss0, twiss, M);
  }
}

void propagateTwissParameters(TWISS *twiss1, TWISS *twiss0, VMATRIX *M)
{
  size_t offset;
  double C, S, Cp, Sp;
  double beta1;
  //double alpha1;
  double beta0, alpha0, eta0, etap0, gamma0;
  double cos_dphi, sin_dphi, dphi;
  long plane;
  
  
  for (plane=offset=0; plane<2; plane++) {
    beta0  = *(&(twiss0->betax )+offset);
    alpha0 = *(&(twiss0->alphax)+offset);
    eta0   = *(&(twiss0->etax  )+offset);
    etap0  = *(&(twiss0->etapx )+offset);
    gamma0 = (1+alpha0*alpha0)/beta0;
    C  = M->R[plane*2  ][plane*2  ];
    S  = M->R[plane*2  ][plane*2+1];
    Cp = M->R[plane*2+1][plane*2  ];
    Sp = M->R[plane*2+1][plane*2+1];
    beta1  = *(&(twiss1->betax )+offset) = C*C*beta0 - 2*S*C*alpha0 + S*S*gamma0;
    //alpha1 = *(&(twiss1->alphax)+offset) = -C*Cp*beta0 + (Sp*C+S*Cp)*alpha0 - S*Sp*gamma0;
    *(&(twiss1->alphax)+offset) = -C*Cp*beta0 + (Sp*C+S*Cp)*alpha0 - S*Sp*gamma0;
    *(&(twiss1->etax )+offset) = eta0*C  + etap0*S  + M->R[plane*2  ][5];
    *(&(twiss1->etapx)+offset) = eta0*Cp + etap0*Sp + M->R[plane*2+1][5];
    if ((sin_dphi = S/sqrt(beta0*beta1))>1) {
      sin_dphi = 1;
      cos_dphi = 0;
    } else if (sin_dphi<-1) {
      sin_dphi = -1;
      cos_dphi = 0;
    } else 
      cos_dphi = sqrt(beta0/beta1)*C - alpha0*sin_dphi;
    if ((dphi = atan2(sin_dphi, cos_dphi))<0)
      dphi += PIx2;
    *(&(twiss1->phix)+offset) = dphi + *(&(twiss0->phix)+offset);
    offset = TWISS_Y_OFFSET;
  }
}


