/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: sdds_beam.c
 * purpose: Do tracking for beams from sdds tables.
 *          See file sdds_beam.nl for input parameters.
 *
 * Michael Borland, 1994
 */
#include "mdb.h"
#include "track.h"
#include "sdds_beam.h"

static long input_type_code;
#define ELEGANT_BEAM 0
#define SPIFFE_BEAM 1
#define N_SDDS_INPUT_TYPES 2
static char *input_type_name[N_SDDS_INPUT_TYPES] = {
  "elegant", "spiffe"
  } ;

#ifdef VAX_VMS
#define isnan(x) 0
#define isinf(x) 0
#endif

  /* SDDS routines will be asked to deliver data in the order given in
   * this string, which matches the order of the #define's 
   */
static char *spiffe_columns = "r pr pz t pphi";
#define ISC_R 0
#define ISC_PR 1
#define ISC_PZ 2
#define ISC_T 3
#define ISC_PPHI 4

/* SDDS routines will be asked to deliver data in the order given in
 * this string, which matches the order of the #define's 
 */
static char *elegant_columns = "x xp y yp t p";
#define IEC_X 0
#define IEC_XP 1
#define IEC_Y 2
#define IEC_YP 3
#define IEC_T 4
#define IEC_P 5

long get_sdds_particles(double ***particle, long one_dump, long n_skip);

static SDDS_TABLE SDDS_input;
static long input_initialized = 0, has_been_read = 0;
char **inputFile = NULL;
long inputFiles = 0;
long inputFileIndex = 0;

void setup_sdds_beam(
                     BEAM *beam,
                     NAMELIST_TEXT *nltext,
                     RUN *run, 
                     VARY *control,
                     ERRORVAL *errcon,
                     OPTIM_VARIABLES *optim,
                     OUTPUT_FILES *output,
                     LINE_LIST *beamline,
                     long n_elements,
                     long save_original
                     )
{
  static long initial_call = 1;

  log_entry("setup_sdds_beam");

  if (!beam)
    bomb("BEAM pointer is null in setup_sdds_beam", NULL);
  if (!nltext)
    bomb("NAMELIST_TEXT pointer is null in setup_sdds_beam", NULL);
  if (!run)
    bomb("RUN pointer is null in setup_sdds_beam", NULL);
  if (!control)
    bomb("VARY pointer is null in setup_sdds_beam", NULL);
  if (!errcon)
    bomb("ERROR pointer is null in setup_sdds_beam", NULL);
  if (!output)
    bomb("OUTPUT_FILES pointer is null in setup_sdds_beam", NULL);
  if (!beamline)
    bomb("beamline pointer is null in setup_sdds_beam", NULL);
  if (!run->runfile || !run->lattice)
    bomb("null runfile or lattice pointer in RUN structure in setup_sdds_beam", NULL);

  if (!initial_call)
    get_sdds_particles(NULL, 0, 0);
  else
    initial_call = 0;

  /* process namelist input */
  set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
  set_print_namelist_flags(0);
  process_namelist(&sdds_beam, nltext);
  print_namelist(stdout, &sdds_beam);
  fflush(stdout);

  /* check for validity of namelist inputs */
  if (input==NULL && input_list==NULL)
    bomb("no input file given in namelist sdds_beam", NULL);
  if ((selection_parameter && !selection_string) || (!selection_parameter && selection_string))
    bomb("must specify selection_parameter and selection_string together or not at all", NULL);

  if (inputFile)
    free(inputFile);
  inputFile = NULL;
  if (input) {
    inputFiles = 1;
    inputFile = tmalloc(sizeof(*inputFile)*1);
    cp_str(inputFile, input);
    inputFile[0] = compose_filename(inputFile[0], run->rootname);
  }
  else {
    char *ptr;
    inputFiles = 0;
    while ((ptr=get_token(input_list))) {
      inputFiles += 1;
      if (!(inputFile = SDDS_Realloc(inputFile, sizeof(*inputFile)*inputFiles)))
        bomb("memory allocation failure", NULL);
      cp_str(inputFile+inputFiles-1, ptr);
      inputFile[inputFiles-1] = compose_filename(inputFile[inputFiles-1], run->rootname);
    }
  }
  
  inputFileIndex = input_initialized = has_been_read = 0;


  if ((input_type_code=match_string(input_type, input_type_name, N_SDDS_INPUT_TYPES, 0))<0)
    bomb("unknown sdds input type", NULL);
  if (input_type_code==SPIFFE_BEAM && n_particles_per_ring<=0)
    bomb("n_particles_per_ring is invalid", NULL);
  if (input_type_code!=SPIFFE_BEAM)
    n_particles_per_ring = 1;
  if (n_particles_per_ring!=1 && one_random_bunch) 
    bomb("must have n_particles_per_ring==1 for one_random_bunch!=0", NULL);
  if (p_lower>p_upper)
    bomb("p_lower and p_upper are invalid", NULL);
  if (sample_interval<1)
    bomb("sample_interval < 1", NULL);
  if (sample_fraction>1)
    bomb("sample_fraction > 1", NULL);
  if (sample_fraction<1 && sample_interval>1)
    bomb("either sample_fraction or sample_interval must be 1", NULL);
  if (save_initial_coordinates && !reuse_bunch)
    save_initial_coordinates = 0;

  beam->original = beam->particle = beam->accepted = NULL;
  beam->n_original = beam->n_to_track = beam->n_accepted = beam->n_saved = beam->n_particle = 0;
  save_initial_coordinates = save_original || save_initial_coordinates;
  
  log_exit("setup_sdds_beam");
}

long new_sdds_beam(
                   BEAM *beam,
                   RUN *run,
                   VARY *control,
                   OUTPUT_FILES *output,
                   long flags
                   )
{    
  double p, t_offset, gamma, beta, p_central;
#ifdef VAX_VMS
  char s[100];
#endif
  long i_store, i, j;
  static long new_particle_data, generate_new_bunch;
  double r, theta, pr, pz, pphi, delta;
  double sin_theta, cos_theta, path, *pti;

  log_entry("new_sdds_beam");

  if (flags&TRACK_PREVIOUS_BUNCH) {
    /* retracking bunch that has already been set up */
    if (!save_initial_coordinates)
      bomb("logic error---initial beam coordinates not saved", NULL);
    if (beam->original==NULL)
      bomb("can't retrack with previous bunch--there isn't one!", NULL);
    if (n_particles_per_ring!=1) {
      fputs("Warning: can't do retracking with previous bunch when n_particles_per_ring!=1\n", stdout);
      fputs("Will use a new bunch generated from previously read data.\n", stdout);
      generate_new_bunch = 1;
    }
    else
      generate_new_bunch = 0;
  }
  else {
    if (!prebunched) {
      /* The beam in the input file is to be treated as a single bunch,
       * even though it may be spread over several pages.
       * Read in all the particles from the input file and allocate arrays
       * for storing initial coordinates and coordinates of accepted 
       * particles. */
      if (beam->original==NULL) {
        /* no beam has been read before, or else it was purged from memory to save RAM */
        /* free any arrays we may have from previous pass */
        if (beam->particle)
          free_zarray_2d((void**)beam->particle, beam->n_particle, 7);
        if (beam->accepted)
          free_zarray_2d((void**)beam->accepted, beam->n_particle, 7);
        beam->particle = beam->accepted = beam->original = NULL;
        /* read the particle data */
        if ((beam->n_original=get_sdds_particles(&beam->original, prebunched, 0))<0) {
          if (has_been_read) 
            return -1;
          bomb("no particles in input file", NULL);
        }
        has_been_read = 1;
        if (save_initial_coordinates || n_particles_per_ring!=1)
          beam->particle = (double**)zarray_2d
            (sizeof(double), beam->n_particle=n_particles_per_ring*beam->n_original, 7);
        else {
          beam->particle = beam->original;
          beam->n_particle = beam->n_original;
        }
        if (run->acceptance)
          beam->accepted = (double**)zarray_2d
            (sizeof(double), beam->n_particle, 7);
        new_particle_data = 1;
      }
      else
        /* since we've read the whole file already and saved it, there is no new data */
        new_particle_data = 0;
    }
    else {
      /* Each page in the input file is to be treated as a separate
       * bunch.  Read in the next page for tracking.
       */
      /* Free arrays from previous pass */
      if (beam->particle)
        free_zarray_2d((void**)beam->particle, beam->n_particle, 7);
      if (beam->accepted)
        free_zarray_2d((void**)beam->accepted, beam->n_particle, 7);
      if (beam->original && beam->original!=beam->particle)
        free_zarray_2d((void**)beam->original, beam->n_original, 7);
      beam->particle = beam->accepted = beam->original = NULL;
      /* read the new page */
      if ((beam->n_original=get_sdds_particles(&beam->original, prebunched, n_tables_to_skip))>=0) { 
        n_tables_to_skip = 0;    /* use the user's parameter only the first time */
        if (save_initial_coordinates || n_particles_per_ring!=1)
          beam->particle = (double**)zarray_2d
            (sizeof(double), beam->n_particle=n_particles_per_ring*beam->n_original, 7);
        else {
          beam->particle = beam->original;
          beam->n_particle = beam->n_original;
        }
        if (run->acceptance) 
          beam->accepted = (double**)zarray_2d
            (sizeof(double), beam->n_particle, 7);
      }
      else { 
        log_exit("new_sdds_beam");
        return(-1);
      }
      has_been_read = 1;
      new_particle_data = 1;
    }
  }

  t_offset = (control->bunch_frequency?(control->i_step-1)/control->bunch_frequency:0);

  p_central = beam->p0_original = run->p_central;
  if (new_particle_data || generate_new_bunch || 
      (input_type_code==SPIFFE_BEAM && !one_random_bunch && !(flags&TRACK_PREVIOUS_BUNCH))) {
    /* Create the initial distribution from the beam->original particle data 
     * or generate a new distribution from those data 
     */
    if (input_type_code==SPIFFE_BEAM) {
      if (!beam->original)
        bomb("beam->original array is NULL (new_sdds_beam-2)", NULL);
      if (!beam->particle)
        bomb("beam->particle array is NULL (new_sdds_beam-2)", NULL);
      for (i=i_store=0; i<beam->n_original; i+=sample_interval) {
        if (!beam->original[i]) {
          fprintf(stdout, "error: beam->original[%ld] is NULL (new_sdds_beam-2)\n", i);;
          fflush(stdout);
          exit(1);
        }
        if (sample_fraction!=1 && random_4(1)>sample_fraction)
          continue;
        pti = beam->original[i];
        pz = pti[ISC_PZ];
        pr = pti[ISC_PR];
        pphi = pti[ISC_PPHI];
        p = sqrt(sqr(pz) + sqr(pr) + sqr(pphi));
        gamma = sqrt(sqr(p)+1);
        if (p_lower && (p_lower>p || p_upper<p))
          continue;
        r      = pti[ISC_R];
        path   = (t_offset + pti[ISC_T])*c_mks*(beta = p/gamma) ;
        delta  = (p-p_central)/p_central;
        theta  = PIx2*random_4(1);
        for (j=0; j<n_particles_per_ring; j++, i_store++) {
          sin_theta = sin(theta);
          cos_theta = cos(theta);
          theta += PIx2/n_particles_per_ring;
          if (!beam->particle[i_store]) {
            fprintf(stdout, "error: beam->particle[%ld] is NULL (new_sdds_beam-2)\n", i_store);
            fflush(stdout);
            exit(1);
          }
          beam->particle[i_store][0] = r*cos_theta;
          beam->particle[i_store][1] = (pr*cos_theta - pphi*sin_theta)/pz;
          beam->particle[i_store][2] = r*sin_theta;
          beam->particle[i_store][3] = (pr*sin_theta + pphi*cos_theta)/pz;
          beam->particle[i_store][4] = path;
          beam->particle[i_store][5] = delta;
          beam->particle[i_store][6] = particleID++;
        }
      }
      beam->n_to_track = i_store;
#ifndef VAX_VMS
      for (i_store=0; i_store<beam->n_to_track; i_store++) {
        for (i=0; i<6; i++) {
          if (!beam->particle[i_store]) {
            fprintf(stdout, "error: beam->particle[%ld] is NULL\n", i_store);
            fflush(stdout);
            exit(1);
          }
          if (isnan(beam->particle[i_store][i]) || isinf(beam->particle[i_store][i])) {
            fprintf(stdout, "error: NaN or Infinity detected in initial particle data, coordinate %ld\n", i);
            fflush(stdout);
            exit(1);
          }
        }
      }
#endif
      if (center_transversely) {
        zero_centroid(beam->particle, beam->n_to_track, 0);
        zero_centroid(beam->particle, beam->n_to_track, 1);
        zero_centroid(beam->particle, beam->n_to_track, 2);
        zero_centroid(beam->particle, beam->n_to_track, 3);
      }
      if (center_arrival_time || reverse_t_sign)
        adjust_arrival_time_data(beam->particle, beam->n_to_track, p_central, 
                                 center_arrival_time, reverse_t_sign);
    }
    else {
      /* In this case, the data is for point particles already,
       * so I just copy the data for the most part, except for sampling.
       */
      if (!beam->original && beam->n_original)
        bomb("beam->original array is NULL (new_sdds_beam)", NULL);
      if (!beam->particle)
        bomb("beam->particle array is NULL (new_sdds_beam)", NULL);
      for (i=i_store=0; i<beam->n_original; i_store++,i+=sample_interval) {
        if (sample_fraction!=1 && random_4(1)>sample_fraction) {
          i_store--;
          continue;
        }
        if (!beam->original[i]) {
          fprintf(stdout, "error: beam->original[%ld] is NULL (new_sdds_beam.2)\n", i);
          fflush(stdout);
          exit(1);
        }
        gamma = sqrt(sqr(p=beam->original[i][IEC_P])+1);
        if (p_lower && (p_lower>p || p_upper<p)) {
          i_store--;
          continue;
        }
        if (!beam->particle[i_store]) {
          fprintf(stdout, "error: beam->particle[%ld] is NULL (new_sdds_beam.2)\n", i_store);
          fflush(stdout);
          exit(1);
        }
        for (j=0; j<4; j++)
          beam->particle[i_store][j] = beam->original[i][j];
        /* convert time to path-length */
        beam->particle[i_store][4] = (t_offset+beam->original[i][IEC_T])*c_mks*p/gamma;
        /* convert energy to dp/p */
        beam->particle[i_store][5] = (p-p_central)/p_central;
        beam->particle[i_store][6] = beam->original[i][6];
      }
      beam->n_to_track = i_store;
#ifndef VAX_VMS
      for (i_store=0; i_store<beam->n_to_track; i_store++) {
        for (i=0; i<6; i++) {
          if (isnan(beam->particle[i_store][i]) || isinf(beam->particle[i_store][i])) {
            fprintf(stdout, "error: NaN or Infinity detected in initial particle data, coordinate %ld\n", i);
            fflush(stdout);
            exit(1);
          }
        }
      }
#endif
      if (center_transversely) {
        zero_centroid(beam->particle, beam->n_to_track, 0);
        zero_centroid(beam->particle, beam->n_to_track, 1);
        zero_centroid(beam->particle, beam->n_to_track, 2);
        zero_centroid(beam->particle, beam->n_to_track, 3);
      }
      if (center_arrival_time || reverse_t_sign)
        adjust_arrival_time_data(beam->particle, beam->n_to_track, p_central, 
                                 center_arrival_time, reverse_t_sign);
    }
  }
  else {
    /* use (x, x', y, x', s, dp/p) saved in beam->original[] */
    if (!(beam->n_to_track = beam->n_saved)) {
      log_exit("new_sdds_beam");
      return -1;
    }
    if (!beam->original)
      bomb("beam->original is NULL (new_sdds_beam.3)", NULL);
    if (!beam->particle)
      bomb("beam->particle is NULL (new_sdds_beam.3)", NULL);
    for (i=0; i<beam->n_saved; i++) {
      if (!beam->original[i]) {
        fprintf(stdout, "error: beam->original[%ld] is NULL (new_sdds_beam.3)\n", i);
        fflush(stdout);
        exit(1);
      }
      if (!beam->particle[i]) {
        fprintf(stdout, "error: beam->particle[%ld] is NULL (new_sdds_beam.3)\n", i);
        fflush(stdout);
        exit(1);
      }
      for (j=0; j<7; j++) 
        beam->particle[i][j] = beam->original[i][j];
      p = p_central*(1+beam->particle[i][5]);
      beta = p/sqrt(p*p+1);
      beam->particle[i][4] += t_offset*beta*c_mks;
    }
    new_particle_data = 0;
  }

  if (new_particle_data && save_initial_coordinates && 
      (one_random_bunch || (reuse_bunch && input_type_code!=SPIFFE_BEAM))) {
    /* Copy the new "initial" data into original[] in case it needs to be reused,  
       but only if the stuff already in original[] is not going to be needed again
       to generate new beams.
       */
    if (beam->original==beam->particle) 
      bomb("logic error in new_sdds_beam: array for original coordinates is missing", NULL);
    if (!beam->original)
      bomb("beam->original is NULL (new_sdds_beam.4)", NULL);
    if (!beam->particle)
      bomb("beam->particle is NULL (new_sdds_beam.4)", NULL);
    for (i=0; i<beam->n_to_track; i++) {
      if (!beam->original[i]) {
        fprintf(stdout, "error: beam->original[%ld] is NULL (new_sdds_beam.4)\n", i);
        fflush(stdout);
        exit(1);
      }
      if (!beam->particle[i]) {
        fprintf(stdout, "error: beam->particle[%ld] is NULL (new_sdds_beam.4)\n", i);
        fflush(stdout);
        exit(1);
      }
      for (j=0; j<7; j++)
        beam->original[i][j] = beam->particle[i][j];
    }
    beam->n_saved = beam->n_to_track;
    new_particle_data = 0;
  }

  if (!save_initial_coordinates) {
    if (reuse_bunch) {
      /* close the SDDS file to free memory and read again from scratch next time */
      if (!SDDS_Terminate(&SDDS_input)) {
        SDDS_SetError("Problem terminate sdds beam input");
        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
      }
      input_initialized = 0;
    }
    /* free the 'original' particle data array */
    if (beam->original && beam->original!=beam->particle) {
      free_zarray_2d((void**)beam->original, beam->n_original, 7);
      beam->original = NULL;
      beam->n_original = 0;
    }
  if (beam->original==beam->particle)
      beam->original = NULL;
  }
  
  log_exit("new_sdds_beam");
  return(beam->n_to_track);
}

/* get_sdds_particles reads data from all of the input files and puts it into the particle array as
 * follows:
 *     for spiffe input:   (*particle)[i] = (r, pr, pz, pphi, t) for ith particle
 *     for elegant input:  (*particle)[i] = (x, xp, y, yp, t, p) for ith particle
 */
long get_sdds_particles(double ***particle, 
                        long one_dump,   /* read only one page */
                        long n_skip      /* number of pages to skip */
                        )
{
  long i, np_max, np, np_new, rows, dump_rejected;
  long retval, data_seen;
  double **data=NULL, **new_data;
  static char s[200];
  long indexID = -1;

  log_entry("get_sdds_particles");

  if (particle==NULL) {
    /* reset for reading again */
    if (input_initialized && !SDDS_Terminate(&SDDS_input)) {
      sprintf(s, "Problem terminating SDDS beam input from file %s", inputFile[inputFileIndex]);
      SDDS_SetError(s);
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
    }
    input_initialized = 0;
    return 0;
  }
  if (inputFileIndex>=inputFiles) {
    fprintf(stdout, "No particles left in input file(s)\n");
    fflush(stdout);
    return -1;
  }
  
  retval = data_seen = np = 0;
  while (inputFileIndex<inputFiles) {
    if (!input_initialized) {
      if (!SDDS_InitializeInputFromSearchPath(&SDDS_input, inputFile[inputFileIndex])) {
        sprintf(s, "Problem opening beam input file %s", inputFile[inputFileIndex]);
        SDDS_SetError(s);
        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
      }
      input_initialized = 1;
      if (selection_parameter) {
        if ((i=SDDS_GetParameterIndex(&SDDS_input, selection_parameter))<0)
          fprintf(stdout, "warning: SDDS beam file %s does not contain the selection parameter %s\n",
                  inputFile[inputFileIndex], selection_parameter);
        fflush(stdout);
        if (SDDS_GetParameterType(&SDDS_input, i)!=SDDS_STRING) {
          sprintf(s, "SDDS beam file %s contains parameter %s, but parameter is not a string", 
                  inputFile[inputFileIndex], selection_parameter);
          SDDS_SetError(s);
          SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
        }
      }
      if (input_type_code==SPIFFE_BEAM) {
        if (!check_sdds_column(&SDDS_input, "r", "m") || 
            !check_sdds_column(&SDDS_input, "pr", "m$be$nc") ||
            !check_sdds_column(&SDDS_input, "pz", "m$be$nc") ||
            !check_sdds_column(&SDDS_input, "pphi", "m$be$nc") ||
            !check_sdds_column(&SDDS_input, "t", "s")) {
          fprintf(stdout, 
                  "necessary data quantities (r, pr, pz, t) have the wrong units or are not present in %s", 
                  inputFile[inputFileIndex]);
          fflush(stdout);
          exit(1);
        }
      }
      else {
        if (!check_sdds_column(&SDDS_input, "x", "m") ||
            !check_sdds_column(&SDDS_input, "y", "m") ||
            !check_sdds_column(&SDDS_input, "xp", NULL) ||
            !check_sdds_column(&SDDS_input, "yp", NULL) ||
            !check_sdds_column(&SDDS_input, "t", "s")) {
          fprintf(stdout, 
                  "necessary data quantities (x, x', y, y', t, p) have the wrong units or are not present in %s\n", 
                  inputFile[inputFileIndex]);
          fflush(stdout);
          exit(1);
        }
        if (!check_sdds_column(&SDDS_input, "p", "m$be$nc")) {
          if (check_sdds_column(&SDDS_input, "p", NULL)) {
            fprintf(stdout, "Warning: p has no units.  Expected m$be$nc\n");
            fflush(stdout);
          }
        }
      }
      fprintf(stdout, "File %s opened and checked.\n", inputFile[inputFileIndex]);
      fflush(stdout);
    }
    
     
    np_max = np = 0;
    data = NULL;
    data_seen = 1;
    while (data_seen) {
      data_seen = 0;
      if ((retval=SDDS_ReadTable(&SDDS_input))==0)
        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
      else if (retval!=-1) 
        data_seen = 1;
      else 
        /* end of file */
        break;
      if (one_dump && n_skip>0) {
        n_skip--;
        continue;
      }
      dump_rejected = 0;
      if (selection_parameter) {
        char *value;
        if (!SDDS_GetParameter(&SDDS_input, selection_parameter, &value)) {
          sprintf(s, "Problem getting value of parameter %s from file %s", selection_parameter, inputFile[inputFileIndex]);
          SDDS_SetError(s);
          SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
        }
        if (!wild_match(value, selection_string)) {
          dump_rejected = 1;
          break;
        }
      }
      if ((rows = SDDS_CountRowsOfInterest(&SDDS_input))<=0) {
        if (rows==-1) {
          sprintf(s, "Problem counting rows of interest for file %s", inputFile[inputFileIndex]);
          SDDS_SetError(s);
          SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
        }
        if (!one_dump)
          continue;
      }
      SDDS_SetColumnFlags(&SDDS_input, 0);
      if (!SDDS_SetColumnsOfInterest(&SDDS_input, SDDS_NAMES_STRING, 
                                     input_type_code==SPIFFE_BEAM?spiffe_columns:elegant_columns)) {
        sprintf(s, "Problem setting columns of interest for file %s", inputFile[inputFileIndex]);
        SDDS_SetError(s);
        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
      }
      if ((np_new=np+rows)>np_max) {
        /* must reallocate to get more space */
        np_max = np + 2*rows;
        data = trealloc(data, np_max*sizeof(*data));
      }
      if (!(new_data=SDDS_GetCastMatrixOfRows(&SDDS_input, &i, SDDS_DOUBLE))) {
        sprintf(s, "Problem getting matrix of rows for file %s", inputFile[inputFileIndex]);
        SDDS_SetError(s);
        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
      }
      if (i!=rows) {
        sprintf(s, "Row count mismatch for file %s", inputFile[inputFileIndex]);
        SDDS_SetError(s);
        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
      }
      for (i=np; i<np_new; i++)
        /* will want to use this storage for 6D+1 phase space latter */
        data[i] = trealloc(new_data[i-np], sizeof(**new_data)*7);
      if ((indexID=SDDS_GetColumnIndex(&SDDS_input, "particleID"))>=0) {
        double *index;
        if (!(index=SDDS_GetColumnInDoubles(&SDDS_input, "particleID"))) {
          sprintf(s, "Problem reading particleID column for file %s", inputFile[inputFileIndex]);
          SDDS_SetError(s);
          SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
        }
        for (i=np; i<np_new; i++)
          data[i][6] = index[i-np];
        free(index);
      }
      else if (input_type_code!=SPIFFE_BEAM)
        for (i=np; i<np_new; i++)
          data[i][6] = particleID++;
      free(new_data);
      np = np_new;
      
      if (one_dump && !dump_rejected)
        break;
    }
    if (retval==-1) {
      /* go to next file */
      if (!SDDS_Terminate(&SDDS_input)) {
        sprintf(s, "Problem terminating SDDS beam input from file %s", inputFile[inputFileIndex]);
        SDDS_SetError(s);
        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
      }
      fprintf(stdout, "File %s was used up and closed.\n", inputFile[inputFileIndex]);
      fflush(stdout);
      inputFileIndex ++;
      input_initialized = 0;
    }
    if (np)
      break;
    fprintf(stdout, "Checking next file\n");
    fflush(stdout);
  }    
  
  if (input_initialized && !SDDS_ShortenTable(&SDDS_input, 1)) {
    SDDS_SetError("Problem releasing table memory when reading SDDS beam file.");
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
  }
  
  if (!data_seen && one_dump)
    return -1;
  
  fprintf(stdout, "a total of %ld data points were read\n\n", np);
  fflush(stdout);
  *particle = data;

  log_exit("get_sdds_particles");
  return(np);
}            

void adjust_arrival_time_data(double **coord, long np, double Po, long center_t, long flip_t)
{
  long ip;
  double P, beta;
  double tave;

  if (np<1)
    return;

  if (flip_t) 
    for (ip=0; ip<np; ip++)
      coord[ip][4] *= -1;
  
  if (center_t) {
    for (ip=tave=0; ip<np; ip++) {
      P = Po*(1+coord[ip][5]);
      beta = P/sqrt(sqr(P)+1);
      tave += coord[ip][4]/beta;
    }
    tave /= np;
    for (ip=0; ip<np; ip++) {
      P = Po*(1+coord[ip][5]);
      beta = P/sqrt(sqr(P)+1);
      coord[ip][4] -= beta*tave;
    }
  }
}
