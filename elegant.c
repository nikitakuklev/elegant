/* Copyright 1994 by Michael Borland and Argonne National Laboratory,
 * all rights reserved.
 */
/* program: elegant
 * purpose: accelerator simulation
 * Michael Borland, 1989-1994
 */
#include "mdb.h"
#include "mdbsun.h"
#include "track.h"
#include "elegant.h"
#include "scan.h"
#include <ctype.h>
#include "match_string.h"
#include <signal.h>
#include <time.h>
#if defined(UNIX) || defined(_WIN32)
#include <malloc.h>
#endif

#include "chromDefs.h"
#include "correctDefs.h"
#include "tuneDefs.h"

void traceback_handler(int code);

#define DESCRIBE_INPUT 0
#define DEFINE_MACRO 1
#define N_OPTIONS  2
char *option[N_OPTIONS] = {
    "describeinput",
    "macro",
        };
char *USAGE="elegant <inputfile> [-macro=<tag>=<value>,[...]]\n\nProgram by Michael Borland. (This is version 14.5, "__DATE__".)";

char *GREETING="This is elegant, by Michael Borland. (This is version 14.5, "__DATE__".)";

#define RUN_SETUP        0
#define RUN_CONTROL      1
#define VARY_ELEMENT     2
#define ERROR_CONTROL    3
#define ERROR_ELEMENT    4
#define SET_AWE_BEAM     5
#define SET_BUNCHED_BEAM 6
#define CORRECTION_SETUP 7
#define MATRIX_OUTPUT    8
#define TWISS_OUTPUT     9
#define TRACK           10
#define STOP            11
#define OPTIMIZATION_SETUP  12
#define OPTIMIZE        13
#define OPTIMIZATION_VARIABLE   14
#define OPTIMIZATION_CONSTRAINT 15
#define OPTIMIZATION_COVARIABLE 16
#define SAVE_LATTICE    17
#define RPN_EXPRESSION  18
#define PROGRAM_TRACE   19
#define CHROMATICITY 20
#define CLOSED_ORBIT    21
#define FIND_APERTURE   22
#define ANALYZE_MAP     23
#define CORRECT_TUNES   24
#define LINK_CONTROL    25
#define LINK_ELEMENTS   26
#define STEERING_ELEMENT 27
#define AMPLIF_FACTORS 28
#define PRINT_DICTIONARY 29
#define FLOOR_COORDINATES 30
#define CORRECTION_MATRIX_OUTPUT 31
#define LOAD_PARAMETERS 32
#define SET_SDDS_BEAM   33
#define SUBPROCESS      34
#define FIT_TRACES      35
#define SASEFEL_AT_END  36
#define ALTER_ELEMENTS  37
#define OPTIMIZATION_TERM 38
#define N_COMMANDS      39

char *command[N_COMMANDS] = {
    "run_setup", "run_control", "vary_element", "error_control", "error_element", "awe_beam", "bunched_beam",
    "correct", "matrix_output", "twiss_output", "track", "stop", 
    "optimization_setup", "optimize", "optimization_variable", "optimization_constraint",
    "optimization_covariable", "save_lattice", "rpn_expression", "trace", "chromaticity", "closed_orbit",
    "find_aperture", "analyze_map", "correct_tunes", "link_control", "link_elements",
    "steering_element", "amplification_factors", "print_dictionary", "floor_coordinates", "correction_matrix_output",
    "load_parameters", "sdds_beam", "subprocess", "fit_traces", "sasefel", "alter_elements",
    "optimization_term",
        } ;

char *description[N_COMMANDS] = {
    "run_setup                   defines lattice input file, primary output files, tracking order, etc.",
    "run_control                 defines number of steps, number of passes, number of indices, etc.",
    "vary_element                defines element family and item to vary with an index",
    "error_control               sets up and control random errors",
    "error_element               defines an element family and item to add errors to",
    "awe_beam                    defines name of input beam data file, type of data, and some preprocessing",
    "bunched_beam                defines beam distribution",
    "correct                     requests orbit or trajectory correction and define parameters",
    "matrix_output               requests awe-format output or printed output of the matrix",
    "twiss_output                requests output of Twiss parameters, chromaticity, and acceptance",
    "track                       command to begin tracking",
    "stop                        command to stop reading input file and end the run",
    "optimization_setup          requests running of optimization mode and sets it up",
    "optimize                    command to begin optimization",
    "optimization_variable       defines an element family and item to vary for optimization",
    "optimization_constraint     defines a constraint on optimization",
    "optimization_covariable     defines an element family and item to compute from optimization variables",
    "optimization_term           specifies an individual term in the optimization equation",
    "save_lattice                command to save the current lattice",
    "rpn_expression              command to execute an rpn expression (useful for optimization)",
    "trace                       requests tracing of program calls and defines parameters of trace",
    "chromaticity                requests correction of the chromaticity",
    "closed_orbit                requests output of the closed orbit",
    "find_aperture               command to do aperture searches",
    "analyze_map                 command to do map analysis",
    "correct_tunes               requests correction of the tunes",
    "link_control                sets up and control element links",
    "link_elements               defines a link between two sets of elements",
    "steering_element            defines an element (group) and item for steering the beam",
    "amplification_factors       computes orbit/trajectory amplification factors",
    "print_dictionary            prints a list of accelerator element types and allowed parameters",
    "floor_coordinates           computes the floor coordinates for the ends of elements",
    "correction_matrix_output    prints response matrices and their inverses",
    "load_parameters             sets up loading of parameter values for elements",
    "sdds_beam                   defines name of input beam data file",
    "subprocess                  executes a string in a sub-shell",
    "fit_traces                  obtains a lattice model by fitting to multiple tracks through a beamline",
    "sasefel                     computes parameters of SASE FEL at end of system",
    "alter_elements              alters a common parameter for one or more elements",
        } ;

void initialize_structures(RUN *run_conditions, VARY *run_control, ERRORVAL *error_control, CORRECTION *correct, 
                           BEAM *beam, OUTPUT_FILES *output_data, OPTIMIZATION_DATA *optimize,
                           CHROM_CORRECTION *chrom_corr_data, TUNE_CORRECTION *tune_corr_data,
                           ELEMENT_LINKS *links);
void free_beamdata(BEAM *beam);

#define NAMELIST_BUFLEN 65536

#define DEBUG 0

int main(argc, argv)
int argc;
char **argv;
{
  char **macroTag, **macroValue;
  long macros;
  LINE_LIST *beamline;        /* pointer to root of linked list */
  FILE *fp_in;
  char *inputfile;
  SCANNED_ARG *scanned;
  char s[NAMELIST_BUFLEN], *ptr;
  long i;
  RUN run_conditions;
  VARY run_control;
  CORRECTION correct;
  ERRORVAL error_control;
  BEAM beam;
  OUTPUT_FILES output_data;
  OPTIMIZATION_DATA optimize;
  CHROM_CORRECTION chrom_corr_data;
  TUNE_CORRECTION tune_corr_data;
  ELEMENT_LINKS links;
  char *saved_lattice = NULL;
  long correction_setuped, run_setuped, run_controled, error_controled, beam_type;
  long do_chromatic_correction = 0, do_twiss_output = 0, fl_do_tune_correction = 0;
  long do_closed_orbit = 0, do_matrix_output = 0, do_response_output = 0;
  long last_default_order = 0, new_beam_flags, links_present, twiss_computed = 0;
  long correctionDone;
  double *starting_coord, finalCharge;
  long namelists_read = 0, failed, firstPass;
                    
#if defined(UNIX) || defined(_WIN32)
    signal(SIGINT, traceback_handler);
    signal(SIGILL, traceback_handler);
    signal(SIGABRT, traceback_handler);
    signal(SIGFPE, traceback_handler);
    signal(SIGSEGV, traceback_handler);
#endif
#if defined(UNIX)
    signal(SIGHUP, traceback_handler);
    signal(SIGQUIT, traceback_handler);
    signal(SIGTRAP, traceback_handler);
    signal(SIGBUS, traceback_handler);
#endif
    
    log_entry("main");
    if (!SDDS_CheckTableStructureSize(sizeof(SDDS_TABLE))) {
        fprintf(stderr, "table structure size is inconsistent\n");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
        exit(1);
        }

    compute_offsets();
    set_max_name_length(12);
    macros = 0;
    macroTag = macroValue = NULL;
  
#if defined(VAX_VMS) || defined(UNIX) || defined(_WIN32)
    init_stats();
#endif

/*
#ifdef SUNOS4
    malloc_debug(1);
      mallopt(M_MXFAST, 1024);
      mallopt(M_NLBLKS, 8);
#endif
 */
    
    argc = scanargs(&scanned, argc, argv);
    if (argc<2) {
        fprintf(stdout, "usage: %s\n", USAGE);
        fflush(stdout);
        link_date();
        exit(1);
        }
    
    fprintf(stdout, "%s\n", GREETING);
    fflush(stdout);
    link_date();
    if (getenv("RPN_DEFNS")) {
        rpn(getenv("RPN_DEFNS"));
        if (rpn_check_error()) exit(1);
        }

    inputfile = NULL;
    for (i=1; i<argc; i++) {
        if (scanned[i].arg_type==OPTION) {
            switch (match_string(scanned[i].list[0], option, N_OPTIONS, 0)) {
              case DESCRIBE_INPUT:
                show_namelists_fields(stdout, namelist_pointer, namelist_name, n_namelists);
                if (argc==2)
                    exit(0);
                break;
              case DEFINE_MACRO:
                if ((scanned[i].n_items-=1)<1)
                  bomb("invalid -macro syntax", USAGE);
                if (!(macroTag=SDDS_Realloc(macroTag, sizeof(*macroTag)*(macros+scanned[i].n_items))) ||
                    !(macroValue=SDDS_Realloc(macroValue, sizeof(*macroValue)*(macros+scanned[i].n_items))))
                  bomb("memory allocation failure (-macro)", NULL);
                else {
                  long j;
                  for (j=0; j<scanned[i].n_items; j++) {
                    macroTag[macros] = scanned[i].list[j+1];
                    if (!(macroValue[macros] = strchr(macroTag[macros], '=')))
                      bomb("invalid -macro syntax", USAGE);
                    macroValue[macros][0] = NULL;
                    macroValue[macros] += 1;
                    macros++;
                  }
                }
                break;
              default:
                bomb("unknown option given.", USAGE);
                break;
                }
            }
        else { 
            /* filenames */
            if (!inputfile)
                fp_in = fopen_e(inputfile = scanned[i].list[0], "r", 0);
            else 
                bomb("too many file names listed.", USAGE);
            }
        }
    
    if (!inputfile)
        bomb("no input file was given", USAGE);
    
    initialize_structures(&run_conditions, &run_control, &error_control, &correct, &beam, &output_data,
                          &optimize, &chrom_corr_data, &tune_corr_data, &links);
    
    run_setuped = run_controled = error_controled = correction_setuped = 0;
    
    starting_coord = tmalloc(sizeof(*starting_coord)*7);
    
    beam_type = -1;

    while (get_namelist(s, NAMELIST_BUFLEN, fp_in)) {
        substituteTagValue(s, NAMELIST_BUFLEN, macroTag, macroValue, macros);
#if DEBUG
      fprintf(stderr, "%s\n", s);
#endif
#if defined(VAX_VMS) || defined(UNIX) || defined(_WIN32)
        report_stats(stdout, "statistics: ");
        fflush(stdout);
#endif
        if (namelists_read)
            free_namelist_text(&namelist_text);
        scan_namelist(&namelist_text, s);
        namelists_read = 1;
        switch (match_string(namelist_text.group_name, command, N_COMMANDS, EXACT_MATCH)) {
          case RUN_SETUP:
            beam_type = -1;
            
            initialize_structures(NULL, &run_control, &error_control, &correct, &beam, &output_data,
                          &optimize, &chrom_corr_data, &tune_corr_data, &links);
            run_setuped = run_controled = error_controled = correction_setuped = 0;
            
            run_setuped = run_controled = error_controled = correction_setuped = do_closed_orbit = do_chromatic_correction = 
                fl_do_tune_correction = 0;
            do_twiss_output = do_matrix_output = do_response_output = 0;

            set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
            set_print_namelist_flags(0);
            process_namelist(&run_setup, &namelist_text);
            print_namelist(stdout, &run_setup);
            setSearchPath(search_path);
            /* check for validity of namelist inputs */
            if (lattice==NULL) {
                if (!saved_lattice)
                    bomb("no lattice supplied", NULL);
                if (default_order!=last_default_order)
                    delete_matrix_data(NULL);
                }
            else {
                /* free previous lattice info */
                free_elements(NULL);
                free_beamlines(NULL);
                saved_lattice = lattice;
                }
            free_beamdata(&beam);
            if (default_order<1 || default_order>3)
                bomb("default_order is out of range", NULL);
            if (concat_order>3)
                bomb("concat_order is out of range", NULL);
            if (p_central<=0 && !expand_for)
                bomb("p_central <= 0", NULL);
            if (expand_for)
                p_central = find_beam_p_central(expand_for);
            if (random_number_seed==0) {
                random_number_seed = (long)time((time_t)0);
                random_number_seed = 2*(random_number_seed/2) + 1;
                fprintf(stdout, "clock-generated random_number_seed = %ld\n", random_number_seed);
                fflush(stdout);
                }
            
            /* seed random number generators.  Note that random_1 seeds random_2, random_3,
             * and random_4.
             * random_1 is used for beamline errors.  
             * random_4 is used for beam generation 
             * random_2 is used for random scraping/sampling/scattering.  
             * random_3 is used for BPM noise.
             */
            random_1(-FABS(random_number_seed));

            /* copy run data into run_conditions structure */
            run_conditions.ideal_gamma = sqrt(sqr(p_central)+1);
            run_conditions.p_central = p_central;
            run_conditions.default_order = default_order;
            run_conditions.concat_order = concat_order;
            run_conditions.print_statistics = print_statistics;
            run_conditions.combine_bunch_statistics = combine_bunch_statistics;
            run_conditions.wrap_around = wrap_around;
            run_conditions.tracking_updates = tracking_updates;
            run_conditions.always_change_p0 = always_change_p0;
            
            /* extract the root filename from the input filename */
            strcpy(s, inputfile);
            if (rootname==NULL) {
                clean_filename(s);
                if ((ptr=strrchr(s, '.')))
                    *ptr = 0;
                cp_str(&rootname, s);
                run_conditions.rootname = rootname;
                run_conditions.runfile  = inputfile;
                }
            else {
                run_conditions.rootname = rootname;
                run_conditions.runfile  = compose_filename(inputfile, rootname);
                }
            run_conditions.acceptance = compose_filename(acceptance, rootname);
            run_conditions.centroid   = compose_filename(centroid, rootname);
            run_conditions.sigma      = compose_filename(sigma, rootname);
            run_conditions.final      = compose_filename(final, rootname);
            run_conditions.output     = compose_filename(output, rootname);
            run_conditions.losses     = compose_filename(losses, rootname);
            magnets                   = compose_filename(magnets, rootname);
            semaphore_file            = compose_filename(semaphore_file, rootname);
            parameters                = compose_filename(parameters, rootname);
            
            if (semaphore_file && fexists(semaphore_file))
              remove(semaphore_file);

            /* parse the lattice file and create the beamline */
            run_conditions.lattice = compose_filename(saved_lattice, rootname);
            beamline = get_beamline(lattice, use_beamline, p_central, echo_lattice, element_divisions);
            fprintf(stdout, "length of beamline %s per pass: %21.15e m\n", beamline->name, beamline->revolution_length);
            fflush(stdout);
            lattice = saved_lattice;
            
            /* output the magnet layout */
            if (magnets)
                output_magnets(magnets, lattice, beamline);

            delete_phase_references();    /* necessary for multi-step runs */
            reset_special_elements(beamline, 1);
            reset_driftCSR();
            last_default_order = default_order;
            run_setuped = 1;
            break;
          case RUN_CONTROL:
            if (!run_setuped)
                bomb("run_setup must precede run_control namelist", NULL);
            vary_setup(&run_control, &namelist_text, &run_conditions, beamline);
            run_controled = 1;
            break;
          case VARY_ELEMENT:
            if (!run_controled)
                bomb("run_control must precede vary_element namelists", NULL);
            if (beam_type!=-1)
                bomb("vary_element statements must come before beam definition", NULL);
            add_varied_element(&run_control, &namelist_text, &run_conditions, beamline);
            break;
          case ERROR_CONTROL:
            if (!run_setuped || !run_controled)
                bomb("run_setup and run_control must precede error_control namelist", NULL);
            if (beam_type!=-1)
                bomb("error specifications must be completed before beam type is specified", NULL);
            error_setup(&error_control, &namelist_text, &run_conditions, beamline);
            error_controled = 1;
            break;                    
          case ERROR_ELEMENT:
            if (beam_type!=-1)
                bomb("error_element statements must come before beam definition", NULL);
            if (!error_controled)
                bomb("error_control namelist must precede error_element namelists", NULL);
            add_error_element(&error_control, &namelist_text, beamline);
            break;
          case CORRECTION_SETUP:
            if (!run_setuped)
                bomb("run_setup must precede correction", NULL);
            if (beam_type!=-1)
                bomb("beam setup (bunched_beam or sdds_beam) must follow correction setup", NULL);
            correction_setuped = 1;
            correction_setup(&correct, &namelist_text, &run_conditions, beamline); 
            delete_phase_references();
            reset_special_elements(beamline, 1);
            reset_driftCSR();
            break;
          case SET_AWE_BEAM: 
            fprintf(stdout, "This program no longer supports awe-format files.\n");
            fflush(stdout);
            fprintf(stdout, "Use awe2sdds to convert your data files, and use\n");
            fflush(stdout);
            fprintf(stdout, "the sdds_beam command instead of awe_beam.\n");
            fflush(stdout);
            break;
          case SET_BUNCHED_BEAM:
            if (!run_setuped || !run_controled)
                bomb("run_setup and run_control must precede bunched_beam namelist", NULL);
            setup_bunched_beam(&beam, &namelist_text, &run_conditions, &run_control, &error_control, &optimize.variables,
                               &output_data, beamline, beamline->n_elems,
                               correct.mode!=-1 && 
                               (correct.track_before_and_after || correct.start_from_centroid));
            setup_output(&output_data, &run_conditions, &run_control, &error_control, &optimize.variables, beamline);
            beam_type = SET_BUNCHED_BEAM;
            break;
          case SET_SDDS_BEAM: 
            if (!run_setuped || !run_controled)
                bomb("run_setup and run_control must precede sdds_beam namelist", NULL);
            setup_sdds_beam(&beam, &namelist_text, &run_conditions, &run_control, &error_control, 
                            &optimize.variables, &output_data, beamline, beamline->n_elems,
                            correct.mode!=-1 && 
                            (correct.track_before_and_after || correct.start_from_centroid));
            setup_output(&output_data, &run_conditions, &run_control, &error_control, &optimize.variables, beamline);
            beam_type = SET_SDDS_BEAM;
            break;
          case TRACK:
            set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
            set_print_namelist_flags(0);
            process_namelist(&track, &namelist_text);
            print_namelist(stdout, &track);
            if ((use_linear_chromatic_matrix || longitudinal_ring_only) 
                && (!twiss_computed && !do_twiss_output))
              bomb("you must compute twiss parameters to do linear chromatic matrix tracking or longitudinal ring tracking", NULL);
            if (use_linear_chromatic_matrix && longitudinal_ring_only)
              bomb("can't do linear chromatic tracking and longitudinal-only tracking together", NULL);
            if (beam_type==-1)
                bomb("beam must be defined prior to tracking", NULL);
            new_beam_flags = 0;
            firstPass = 1;
            while (vary_beamline(&run_control, &error_control, &run_conditions, beamline)) {
              fill_double_array(starting_coord, 6, 0.0);
              correctionDone = 0;
              if (correct.mode!= -1) {
                if (correct.track_before_and_after || correct.start_from_centroid) {
                  if (beam_type==SET_AWE_BEAM) {
                    bomb("beam type of SET_AWE_BEAM in main routine--this shouldn't happen", NULL);
                  }
                  else if (beam_type==SET_SDDS_BEAM) {
                    if (new_sdds_beam(&beam, &run_conditions, &run_control, &output_data, 0)<0)
                      break;
                  }
                  else
                    new_bunched_beam(&beam, &run_conditions, &run_control, &output_data, 0);
                  new_beam_flags = TRACK_PREVIOUS_BUNCH;
                  if (correct.start_from_centroid)
                    compute_centroids(starting_coord, beam.particle, beam.n_to_track);
                  if (correct.track_before_and_after) {
                    track_beam(&run_conditions, &run_control, &error_control, &optimize.variables, 
                               beamline, &beam, &output_data, 
                               PRECORRECTION_BEAM, 0, &finalCharge);
                  }
                }
                else if (correct.use_actual_beam) {
                  if (beam_type==SET_AWE_BEAM) {
                    bomb("beam type of SET_AWE_BEAM in main routine--this shouldn't happen", NULL);
                  }
                  else if (beam_type==SET_SDDS_BEAM) {
                    if (new_sdds_beam(&beam, &run_conditions, &run_control, &output_data, 0)<0)
                      break;
                  }
                  else
                    new_bunched_beam(&beam, &run_conditions, &run_control, &output_data, 0);
                  new_beam_flags = TRACK_PREVIOUS_BUNCH;
                }
                if (!do_correction(&correct, &run_conditions, beamline, starting_coord, &beam, 
                                   run_control.i_step, 
                                   (fl_do_tune_correction || do_chromatic_correction)) ) {
                  fputs("warning: orbit correction failed--continuing with next step\n", stdout);
                  continue;
                }
                correctionDone = 1;
              }
              if (beam_type==SET_AWE_BEAM) {
                bomb("beam type of SET_AWE_BEAM in main routine--this shouldn't happen", NULL);
              }
              else if (beam_type==SET_SDDS_BEAM) {
                if (new_sdds_beam(&beam, &run_conditions, &run_control, &output_data, new_beam_flags)<0)
                  break;
              }
              else
                new_bunched_beam(&beam, &run_conditions, &run_control, &output_data, new_beam_flags);
              if (do_closed_orbit && 
                  !run_closed_orbit(&run_conditions, beamline, starting_coord, &beam, 0) &&
                  !soft_failure) {
                fprintf(stdout, "Closed orbit not found---continuing to next step\n");
                fflush(stdout);
                continue;
              }
              if (do_twiss_output && 
                  !run_twiss_output(&run_conditions, beamline, starting_coord, 0) &&
                  !soft_failure) {
                fprintf(stdout, "Twiss parameters not defined---continuing to next step\n");
                fflush(stdout);
                continue;
              }
              if (do_response_output)
                run_response_output(&run_conditions, beamline, &correct, 0);
              run_matrix_output(&run_conditions, beamline);
              for (i=failed=0; (fl_do_tune_correction || do_chromatic_correction) && i<correction_iterations; i++) {
                correctionDone = 0;
                if (correction_iterations>1) {
                  fprintf(stdout, "\nTune/chromaticity correction iteration %ld\n", i+1);
                  fflush(stdout);
                }
                if (fl_do_tune_correction) {
                  if (do_closed_orbit && 
                      !run_closed_orbit(&run_conditions, beamline, starting_coord, &beam, 0) &&
                      !soft_failure) {
                    fprintf(stdout, "Closed orbit not found---continuing to next step\n");
                    fflush(stdout);
                    failed = 1;
                    break;
                  }
                  if (!do_tune_correction(&tune_corr_data, &run_conditions, beamline, starting_coord,
                                          run_control.i_step, i==correction_iterations-1) &&
                      !soft_failure) {
                    fprintf(stdout, "Tune correction failed---continuing to next step\n");
                    fflush(stdout);
                    failed = 1;
                    break;
                  }
                }
                if (do_chromatic_correction) {
                  if (do_closed_orbit && 
                      !run_closed_orbit(&run_conditions, beamline, starting_coord, &beam, 0) &&
                      !soft_failure) {
                    fprintf(stdout, "Closed orbit not found---continuing to next step\n");
                    fflush(stdout);
                    failed = 1;
                    break;
                  }
                  if (!do_chromaticity_correction(&chrom_corr_data, &run_conditions, beamline, starting_coord,
                                                  run_control.i_step, i==correction_iterations-1) &&
                      !soft_failure) {
                    fprintf(stdout, "Chromaticity correction failed---continuing to next step\n");
                    fflush(stdout);
                    failed = 1;
                    break;
                  }
                }
              }
              if (failed)
                continue;
              perturb_beamline(&run_control, &error_control, &run_conditions, beamline); 
              if (correct.mode!=-1 && !correctionDone &&
                  !do_correction(&correct, &run_conditions, beamline, starting_coord, &beam, run_control.i_step, 0) &&
                  !soft_failure) {
                fputs("warning: orbit correction failed--continuing with next step\n", stdout);
                continue;
              }
              if (do_closed_orbit && 
                  !run_closed_orbit(&run_conditions, beamline, starting_coord, &beam, 1) &&
                  !soft_failure) {
                fprintf(stdout, "Closed orbit not found---continuing to next step\n");
                fflush(stdout);
                continue;
              }
              if (do_twiss_output && 
                  !run_twiss_output(&run_conditions, beamline, starting_coord, 1) &&
                  !soft_failure) {
                fprintf(stdout, "Twiss parameters not defined---continuing to next step\n");
                fflush(stdout);
                continue;
              }
              if (do_response_output)
                run_response_output(&run_conditions, beamline, &correct, 1);
              if (center_on_orbit)
                center_beam_on_coords(beam.particle, beam.n_to_track, starting_coord, center_momentum_also);
              if (firstPass) {
                /* prevent fiducialization of RF etc. by correction etc. */
                delete_phase_references();
                reset_special_elements(beamline, 1);
                reset_driftCSR();
              }
              firstPass = 0;
              track_beam(&run_conditions, &run_control, &error_control, &optimize.variables, 
                         beamline, &beam, &output_data, 
                         (use_linear_chromatic_matrix?LINEAR_CHROMATIC_MATRIX:0)+
                         (longitudinal_ring_only?LONGITUDINAL_RING_ONLY:0), 0, &finalCharge);
              if (parameters)
                dumpLatticeParameters(parameters, &run_conditions, beamline);
            }
            finish_output(&output_data, &run_conditions, &run_control, &error_control, &optimize.variables, 
                          beamline, beamline->n_elems, &beam, finalCharge);
            free_beamdata(&beam);
            if (do_closed_orbit)
              finish_clorb_output();
            if (do_twiss_output)
              finish_twiss_output();
            if (do_response_output)
              finish_response_output();
            if (parameters)
              finishLatticeParametersFile();
#ifdef SUNOS4
            check_heap();
#endif
            fprintf(stdout, "Finished tracking.\n");
            fflush(stdout);
            /* reassert defaults for namelist run_setup */
            lattice = use_beamline = acceptance = centroid = sigma = final = output = rootname = losses = 
              parameters = NULL;
            combine_bunch_statistics = 0;
            random_number_seed = 987654321;
            wrap_around = 1;
            default_order = 2;
            concat_order = 0;
            tracking_updates = 1;
            concat_order = print_statistics = p_central = 0;
            run_setuped = run_controled = error_controled = correction_setuped = do_chromatic_correction =
              fl_do_tune_correction = do_closed_orbit = do_twiss_output = do_response_output = 0;
            break;
          case MATRIX_OUTPUT:
            if (!run_setuped)
                bomb("run_setup must precede matrix_output namelist", NULL);
            setup_matrix_output(&namelist_text, &run_conditions, beamline);
            do_matrix_output = 1;
            break;
          case TWISS_OUTPUT:
            if (!run_setuped)
                bomb("run_setup must precede twiss_output namelist", NULL);
            setup_twiss_output(&namelist_text, &run_conditions, beamline, &do_twiss_output,
                               run_conditions.default_order);
            if (!do_twiss_output) {
                twiss_computed = 1;
                run_twiss_output(&run_conditions, beamline, NULL, -1);
                delete_phase_references();
                reset_special_elements(beamline, 1);
                reset_driftCSR();
                finish_twiss_output();
            }
            break;
          case STOP:
            lorentz_report();
            finish_load_parameters();
            if (semaphore_file && !fexists(semaphore_file)) {
              fclose(fopen(semaphore_file, "w"));
            }
            free_beamdata(&beam);
            exit(0);
            break;
          case OPTIMIZATION_SETUP:
            if (beam_type!=-1)
                bomb("optimization statements must come before beam definition", NULL);
            do_optimization_setup(&optimize, &namelist_text, &run_conditions, beamline);
            break;
          case OPTIMIZE:
            if (beam_type==-1)
                bomb("beam definition must come before optimize command", NULL);
            while (vary_beamline(&run_control, &error_control, &run_conditions, beamline)) {
                do_optimize(&namelist_text, &run_conditions, &run_control, &error_control, beamline, &beam,
                            &output_data, &optimize, beam_type);
                if (parameters)
                  dumpLatticeParameters(parameters, &run_conditions, beamline);
                }
            if (parameters)
              finishLatticeParametersFile();
            break;
          case OPTIMIZATION_VARIABLE:
            if (beam_type!=-1)
                bomb("optimization statements must come before beam definition", NULL);
            add_optimization_variable(&optimize, &namelist_text, &run_conditions, beamline);
            break;
          case OPTIMIZATION_CONSTRAINT:
            if (beam_type!=-1)
                bomb("optimization statements must come before beam definition", NULL);
            add_optimization_constraint(&optimize, &namelist_text, &run_conditions, beamline);
            break;
          case OPTIMIZATION_COVARIABLE:
            if (beam_type!=-1)
                bomb("optimization statements must come before beam definition", NULL);
            add_optimization_covariable(&optimize, &namelist_text, &run_conditions, beamline);
            break;
          case OPTIMIZATION_TERM:
            if (beam_type!=-1)
                bomb("optimization statements must come before beam definition", NULL);
            add_optimization_term(&optimize, &namelist_text, &run_conditions, beamline);
            break;
          case SAVE_LATTICE:
            do_save_lattice(&namelist_text, &run_conditions, beamline);
            break;
          case RPN_EXPRESSION:
            run_rpn_expression(&namelist_text);
            break;
          case PROGRAM_TRACE:
            process_trace_request(&namelist_text);
            break;
          case CHROMATICITY:
            setup_chromaticity_correction(&namelist_text, &run_conditions, beamline, &chrom_corr_data);
            do_chromatic_correction = 1;
            break;
          case CORRECT_TUNES:
            setup_tune_correction(&namelist_text, &run_conditions, beamline, &tune_corr_data);
            fl_do_tune_correction = 1;
            break;
          case CLOSED_ORBIT:
            setup_closed_orbit(&namelist_text, &run_conditions, beamline);
            do_closed_orbit = 1;
            if (correction_setuped)
                fprintf(stdout, "warning: you've asked to do both closed-orbit calculation and orbit correction.\nThis may duplicate effort.\n");
                fflush(stdout);
            break;
          case FIND_APERTURE:
            setup_aperture_search(&namelist_text, &run_conditions, &run_control);
            while (vary_beamline(&run_control, &error_control, &run_conditions, beamline)) {
              fill_double_array(starting_coord, 6, 0.0);
              if (correct.mode!= -1) {
                if (!do_correction(&correct, &run_conditions, beamline, starting_coord, &beam, 
                                   run_control.i_step, 1) ) {
                  fputs("warning: orbit correction failed--continuing with next step\n", stdout);
                  continue;
                }
              }
              if (do_closed_orbit && 
                  !run_closed_orbit(&run_conditions, beamline, starting_coord, &beam, 0) &&
                  !soft_failure) {
                fprintf(stdout, "Closed orbit not found---continuing to next step\n");
                fflush(stdout);
                continue;
              }
              if (do_twiss_output && 
                  !run_twiss_output(&run_conditions, beamline, starting_coord, 0) &&
                  !soft_failure) {
                fprintf(stdout, "Twiss parameters not defined---continuing to next step\n");
                fflush(stdout);
                continue;
              }
              if (do_response_output)
                run_response_output(&run_conditions, beamline, &correct, 0);
              run_matrix_output(&run_conditions, beamline);
              for (i=failed=0; (fl_do_tune_correction || do_chromaticity_correction) && i<correction_iterations; i++) {
                if (correction_iterations>1)
                  fprintf(stdout, "\nTune/chromaticity correction iteration %ld\n", i+1);
                  fflush(stdout);
                if (fl_do_tune_correction) {
                  if (do_closed_orbit && 
                      !run_closed_orbit(&run_conditions, beamline, starting_coord, &beam, 0) &&
                      !soft_failure) {
                    fprintf(stdout, "Closed orbit not found---continuing to next step\n");
                    fflush(stdout);
                    failed = 1;
                    break;
                  }
                  if (!do_tune_correction(&tune_corr_data, &run_conditions, beamline, starting_coord,
                                          run_control.i_step, i==correction_iterations-1) &&
                      !soft_failure) {
                    fprintf(stdout, "Tune correction failed---continuing to next step\n");
                    fflush(stdout);
                    failed = 1;
                    break;
                  }
                }
                if (do_chromatic_correction) {
                  if (do_closed_orbit && 
                      !run_closed_orbit(&run_conditions, beamline, starting_coord, &beam, 0) &&
                      !soft_failure) {
                    fprintf(stdout, "Closed orbit not found---continuing to next step\n");
                    fflush(stdout);
                    failed = 1;
                    break;
                  }
                  if (!do_chromaticity_correction(&chrom_corr_data, &run_conditions, beamline, starting_coord,
                                                  run_control.i_step, i==correction_iterations-1) &&
                      !soft_failure) {
                    fprintf(stdout, "Chromaticity correction failed---continuing to next step\n");
                    fflush(stdout);
                    failed = 1;
                    break;
                  }
                }
              }
              if (failed)
                continue;
              perturb_beamline(&run_control, &error_control, &run_conditions, beamline); 
              if (correct.mode!=-1 &&
                  !do_correction(&correct, &run_conditions, beamline, starting_coord, &beam, run_control.i_step, 0) &&
                  !soft_failure) {
                fputs("warning: orbit correction failed--continuing with next step\n", stdout);
                continue;
              }
              if (do_closed_orbit && 
                  !run_closed_orbit(&run_conditions, beamline, starting_coord, &beam, 1) &&
                  !soft_failure) {
                fprintf(stdout, "Closed orbit not found---continuing to next step\n");
                fflush(stdout);
                continue;
              }
              if (do_twiss_output && 
                  !run_twiss_output(&run_conditions, beamline, starting_coord, 1) &&
                  !soft_failure) {
                fprintf(stdout, "Twiss parameters not defined---continuing to next step\n");
                fflush(stdout);
                continue;
              }
              if (do_response_output)
                run_response_output(&run_conditions, beamline, &correct, 1);
              do_aperture_search(&run_conditions, &run_control, &error_control, beamline);
            }
            fprintf(stdout, "Finished all tracking steps.\n"); fflush(stdout);
            fflush(stdout);
            finish_aperture_search(&run_conditions, &run_control, &error_control, beamline);
            if (do_closed_orbit)
                finish_clorb_output();
#ifdef SUNOS4
            check_heap();
#endif
            fprintf(stdout, "Finished dynamic aperture search.\n");
            fflush(stdout);
            /* reassert defaults for namelist run_setup */
            lattice = use_beamline = acceptance = centroid = sigma = final = output = rootname = losses =
              parameters = NULL;
            combine_bunch_statistics = 0;
            random_number_seed = 987654321;
            wrap_around = 1;
            default_order = 2;
            concat_order = 0;
            tracking_updates = 1;
            concat_order = print_statistics = p_central = 0;
            run_setuped = run_controled = error_controled = correction_setuped = do_chromatic_correction =
                fl_do_tune_correction = do_closed_orbit = do_twiss_output = do_response_output = 0;
            break;
          case ANALYZE_MAP:
            setup_transport_analysis(&namelist_text, &run_conditions, &run_control, &error_control);
            while (vary_beamline(&run_control, &error_control, &run_conditions, beamline)) {
                if (run_control.i_vary<=1) {
                    fill_double_array(starting_coord, 6, 0.0);
                    if (correct.mode!= -1) {
                        if (!do_correction(&correct, &run_conditions, beamline, starting_coord, NULL, run_control.i_step, 1)) {
                            fputs("warning: correction failed--continuing with next step", stdout);
                            continue;
                            }
                        }
                    }
                if (do_twiss_output)
                    run_twiss_output(&run_conditions, beamline, starting_coord, 0);
                if (do_response_output)
                    run_response_output(&run_conditions, beamline, &correct, 0);
                for (i=0; (fl_do_tune_correction || do_chromaticity_correction) && i<correction_iterations; i++) {
                    if (correction_iterations>1)
                        fprintf(stdout, "\nTune/chromaticity correction iteration %ld\n", i+1);
                        fflush(stdout);
                    if (fl_do_tune_correction) {
                        if (do_closed_orbit)
                            run_closed_orbit(&run_conditions, beamline, starting_coord, NULL, 0);
                        do_tune_correction(&tune_corr_data, &run_conditions, beamline, starting_coord,
                                           run_control.i_step, i==correction_iterations-1);
                        }
                    if (do_chromatic_correction) {
                        if (do_closed_orbit)
                            run_closed_orbit(&run_conditions, beamline, starting_coord, NULL, 0);
                        do_chromaticity_correction(&chrom_corr_data, &run_conditions, beamline, starting_coord,
                                                   run_control.i_step, i==correction_iterations-1);
                        }
                    }
                perturb_beamline(&run_control, &error_control, &run_conditions, beamline);
                if (correct.mode!=-1 &&
                    !do_correction(&correct, &run_conditions, beamline, starting_coord, &beam, run_control.i_step, 0)) {
                  fputs("warning: orbit correction failed--continuing with next step\n", stdout);
                  continue;
                }
                if (do_closed_orbit)
                    run_closed_orbit(&run_conditions, beamline, starting_coord, NULL, 1);
                if (do_twiss_output)
                    run_twiss_output(&run_conditions, beamline, starting_coord, 1);
                if (do_response_output)
                    run_response_output(&run_conditions, beamline, &correct, 1);
                do_transport_analysis(&run_conditions, &run_control, &error_control, beamline, 
                                      (do_closed_orbit || correct.mode!=-1?starting_coord:NULL));
                }
            finish_transport_analysis(&run_conditions, &run_control, &error_control, beamline);
            if (do_closed_orbit)
                finish_clorb_output();
            if (do_twiss_output)
                finish_twiss_output();
            if (do_response_output)
                finish_response_output();
#ifdef SUNOS4
            check_heap();
#endif
            fprintf(stdout, "Finished transport analysis.\n");
            fflush(stdout);
            /* reassert defaults for namelist run_setup */
            lattice = use_beamline = acceptance = centroid = sigma = final = output = rootname = losses = 
              parameters = NULL;
            combine_bunch_statistics = 0;
            random_number_seed = 987654321;
            wrap_around = 1;
            default_order = 2;
            concat_order = 0;
            tracking_updates = 1;
            concat_order = print_statistics = p_central = 0;
            run_setuped = run_controled = error_controled = correction_setuped = do_chromatic_correction =
                fl_do_tune_correction = do_closed_orbit = do_twiss_output = do_matrix_output = do_twiss_output = 0;
            break;
          case LINK_CONTROL:
            if (!run_setuped || !run_controled)
                bomb("run_setup and run_control must precede link_control namelist", NULL);
            element_link_control(&links, &namelist_text, &run_conditions, beamline);
            break;                    
          case LINK_ELEMENTS:
            if (!run_setuped || !run_controled)
                bomb("run_setup and run_control must precede link_elements namelist", NULL);
            if (!beamline)
                bomb("beamline not defined--can't add element links", NULL);
            add_element_links(&links, &namelist_text, beamline);
            links_present = 1;
            beamline->links = &links;
            break;
          case STEERING_ELEMENT:
            if (correction_setuped)
                bomb("you must define steering elements prior to giving the 'correct' namelist", NULL);
            add_steering_element(&correct, beamline, &run_conditions, &namelist_text);
            break;
          case AMPLIF_FACTORS:
            compute_amplification_factors(&namelist_text, &run_conditions, &correct, do_closed_orbit, beamline);
            break;
          case PRINT_DICTIONARY:
            set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
            set_print_namelist_flags(0);
            process_namelist(&print_dictionary, &namelist_text);
            print_namelist(stdout, &print_dictionary);
            do_print_dictionary(filename, latex_form);
            break;
          case FLOOR_COORDINATES:
            output_floor_coordinates(&namelist_text, &run_conditions, beamline);
            break;
          case CORRECTION_MATRIX_OUTPUT:
            if (!run_setuped)
                bomb("run setup must precede correction_matrix_output namelist", NULL);
            setup_correction_matrix_output(&namelist_text, &run_conditions, beamline, &correct,
                                           &do_response_output, 
                                           do_twiss_output+do_matrix_output+twiss_computed);
            if (!do_response_output) {
                run_response_output(&run_conditions, beamline, &correct, -1);
                delete_phase_references();
                reset_special_elements(beamline, 1);
                reset_driftCSR();
                finish_response_output();
                }
            break;
          case LOAD_PARAMETERS:
            if (!run_setuped)
                bomb("run_setup must precede load_parameters namelists", NULL);
            if (run_controled)
                bomb("load_parameters namelists must precede run_control namelist", NULL);
            if (error_controled)
                bomb("load_parameters namelists must precede error_control and error namelists", NULL);
            if (setup_load_parameters(&namelist_text, &run_conditions, beamline) && magnets)
              /* make sure the magnet output is right in case loading parameters changed something */
              output_magnets(magnets, lattice, beamline);
            break;
          case SUBPROCESS:
            run_subprocess(&namelist_text, &run_conditions);
            break;
          case FIT_TRACES:
            do_fit_trace_data(&namelist_text, &run_conditions, beamline);
            if (parameters) {
              dumpLatticeParameters(parameters, &run_conditions, beamline);
              finishLatticeParametersFile();
            }
            /* reassert defaults for namelist run_setup */
            lattice = use_beamline = acceptance = centroid = sigma = final = output = rootname = losses =
              parameters = NULL;
            combine_bunch_statistics = 0;
            random_number_seed = 987654321;
            wrap_around = 1;
            default_order = 2;
            concat_order = 0;
            tracking_updates = 1;
            concat_order = print_statistics = p_central = 0;
            run_setuped = run_controled = error_controled = correction_setuped = do_chromatic_correction =
              fl_do_tune_correction = do_closed_orbit = do_twiss_output = do_response_output = 0;
            break;
          case SASEFEL_AT_END:
            if (!run_setuped)
              bomb("run_setup must precede sasefel namelist", NULL);
            if (beam_type!=-1)
              bomb("sasefel namelist must precede beam definition", NULL);
            setupSASEFELAtEnd(&namelist_text, &run_conditions, &output_data);
            break;
          case ALTER_ELEMENTS:
            if (!run_setuped)
              bomb("run_setup must precede alter_element namelist", NULL);
            do_alter_element(&namelist_text, &run_conditions, beamline);
            break;
          default:
            fprintf(stdout, "unknown namelist %s given.  Known namelists are:\n", namelist_text.group_name);
            fflush(stdout);
            for (i=0; i<N_COMMANDS; i++)
                fprintf(stdout, "%s\n", description[i]);
                fflush(stdout);
            exit(1);
            break;
            }
#ifdef SUNOS4
        check_heap();
#endif
        }
    fprintf(stdout, "End of input data encountered.\n"); fflush(stdout);
    fflush(stdout);
    lorentz_report();
    finish_load_parameters();
    if (semaphore_file && !fexists(semaphore_file)) {
      fclose(fopen(semaphore_file, "w"));
    }
    free_beamdata(&beam);
#if defined(VAX_VMS) || defined(UNIX) || defined(_WIN32)
    report_stats(stdout, "statistics: ");
    fflush(stdout);
#endif
    log_exit("main");
    return(0);
    }

double find_beam_p_central(char *input)
{
  SDDS_DATASET SDDSin;
  char s[SDDS_MAXLINE];
  double *p, psum;
  long i, rows;
  
  if (!SDDS_InitializeInputFromSearchPath(&SDDSin, input) || !SDDS_ReadPage(&SDDSin)) {
    sprintf(s, "Problem opening beam input file %s", input);
    SDDS_SetError(s);
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
  }
  if ((rows=SDDS_RowCount(&SDDSin))<=0 || !(p=SDDS_GetColumnInDoubles(&SDDSin, "p"))) {
    sprintf(s, "No data in input file %s", input);
    SDDS_SetError(s);
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
  }
  for (i=psum=0; i<rows; i++) 
    psum += p[i];
  if (SDDS_ReadPage(&SDDSin)>0)
    fprintf(stdout, "Warning: file %s has multiple pages.  Only the first is used for expand_for.\n",
            input);
    fflush(stdout);
  SDDS_Terminate(&SDDSin);
  fprintf(stdout, "Expanding about p = %21.15e\n", psum/rows);
  fflush(stdout);
  return psum/rows;
}

#ifdef SUNOS4
#include <malloc.h>

void check_heap() 
{
    struct mallinfo info;
    
    fprintf(stdout, "Performing memory heap verification..."); fflush(stdout);
    fflush(stdout);
    if (malloc_verify())
        fputs("okay.", stdout);
    else
        fputs("errors, detected.", stdout);
    fflush(stdout);
#if defined(VAX_VMS) || defined(UNIX) || defined(_WIN32)
    report_stats(stdout, "statistics: ");
    fflush(stdout);
#endif
    return;
    info.arena = info.ordblks = info.smblks = info.hblks = info.hblkhd = info.usmblks = 0;   
    info.fsmblks = info.uordblks = info.fordblks = info.keepcost = info.allocated = info.treeoverhead = 0;
    info = mallinfo();
    fprintf(stdout, "memory allocation information:\n");
    fflush(stdout);
    fprintf(stdout, "  total space in arena: %ld\n", info.arena);     
    fflush(stdout);
    fprintf(stdout, "  number of ordinary blocks: %ld\n", info.ordblks);   
    fflush(stdout);
    fprintf(stdout, "  number of small blocks: %ld\n", info.smblks);    
    fflush(stdout);
    fprintf(stdout, "  number of holding blocks: %ld\n", info.hblks);     
    fflush(stdout);
    fprintf(stdout, "  space in holding block headers: %ld\n", info.hblkhd);    
    fflush(stdout);
    fprintf(stdout, "  space in small blocks in use: %ld\n", info.usmblks);   
    fflush(stdout);
    fprintf(stdout, "  space in free small blocks: %ld\n", info.fsmblks);   
    fflush(stdout);
    fprintf(stdout, "  space in ordinary blocks in use: %ld\n", info.uordblks);  
    fflush(stdout);
    fprintf(stdout, "  space in free ordinary blocks: %ld\n", info.fordblks);  
    fflush(stdout);
    fprintf(stdout, "  cost of enabling keep option: %ld\n", info.keepcost);  
    fflush(stdout);
    fprintf(stdout, "  number of ordinary blocks allocated: %ld\n", info.allocated);
    fflush(stdout);
    fprintf(stdout, "  bytes used in maintaining the free tree: %ld\n", info.treeoverhead);
    fflush(stdout);
    }
#endif

void center_beam_on_coords(double **part, long np, double *coord, long center_dp)
{
    double sum, offset;
    long i, j, lim;
    
    if (!np)
        return;
    
    if (center_dp)
        lim = 5;
    else
        lim = 4;
    
    for (j=0; j<=lim; j++) {
        for (i=sum=0; i<np; i++)
            sum += part[i][j];
        offset = sum/np - coord[j];
        for (i=0; i<np; i++)
            part[i][j] -= offset;
        }
    }

typedef struct {
  char *elementName;
  long index;
} DICTIONARY_ENTRY;

int dictionaryEntryCmp(const void *p1, const void *p2)
{
  DICTIONARY_ENTRY *de1, *de2;
  de1 = (DICTIONARY_ENTRY *)p1;
  de2 = (DICTIONARY_ENTRY *)p2;
  return strcmp(de1->elementName, de2->elementName);
}

void do_print_dictionary(char *filename, long latex_form)
{
  FILE *fp;
  long i;
  DICTIONARY_ENTRY *dictList;
  
  if (!filename)
    bomb("filename invalid (do_print_dictionary)", NULL);
  if (!(dictList = malloc(sizeof(*dictList)*N_TYPES)))
    bomb("memory allocation failure", NULL);
  for (i=1; i<N_TYPES; i++) {
    dictList[i-1].elementName = entity_name[i];
    dictList[i-1].index = i;
  }
  qsort((void*)dictList, N_TYPES-1, sizeof(*dictList), dictionaryEntryCmp);
  fp = fopen_e(filename, "w", 0);
  if (latex_form) {
    fprintf(fp, "\\newlength{\\descwidth}\n");
    fprintf(fp, "\\setlength{\\descwidth}{2in}\n");
  }
  for (i=0; i<N_TYPES-1; i++)
    print_dictionary_entry(fp, dictList[i].index, latex_form);
  free(dictList);
  fclose(fp);
}

#define PRINTABLE_NULL(s) (s?s:"NULL")
char *translateUnitsToTex(char *source);
char *makeTexSafeString(char *source);

void print_dictionary_entry(FILE *fp, long type, long latex_form)
{
  char *type_name[3] = {"double", "long", "STRING", };
  long j, texLines;
  if (latex_form) {
    fprintf(fp, "\\begin{latexonly}\n\\newpage\n\\begin{center}{\\Large\\verb|%s|}\\end{center}\n\\end{latexonly}\\subsection{%s}\n", 
            entity_name[type], entity_name[type]);
    fprintf(fp, "%s\n\\\\\n", makeTexSafeString(entity_text[type]));
    fprintf(fp, "\\begin{tabular}{|l|l|l|l|p{\\descwidth}|} \\hline\n");
    fprintf(fp, "Parameter Name & Units & Type & Default & Description \\\\ \\hline \n");
  }
  else  {
    fprintf(fp, "***** element type %s:\n", entity_name[type]);
    fprintf(fp, "%s\n", entity_text[type]);
  }
  for (j=texLines=0; j<entity_description[type].n_params; j++) {
    /* 35 lines fits on a latex page */
    if (latex_form && texLines>35) {
      texLines = 0;
      fprintf(fp, "\\end{tabular}\n\n");
      fprintf(fp, "\\begin{latexonly}\n\\newpage\n\\begin{center}{\\Large\\verb|%s| continued}\\end{center}\n\\end{latexonly}\n", 
              entity_name[type], entity_name[type]);
      fprintf(fp, "%s\n\\\\\n", makeTexSafeString(entity_text[type]));
      fprintf(fp, "\\begin{tabular}{|l|l|l|l|p{\\descwidth}|} \\hline\n");
      fprintf(fp, "Parameter Name & Units & Type & Default & Description \\\\ \\hline \n");
    }
    if (!latex_form)
      fprintf(fp, "%20s %20s %10s", 
              PRINTABLE_NULL(entity_description[type].parameter[j].name), 
              PRINTABLE_NULL(entity_description[type].parameter[j].unit),
              PRINTABLE_NULL(type_name[entity_description[type].parameter[j].type-1]));
    else {
      fprintf(fp, "%s ",
              makeTexSafeString(PRINTABLE_NULL(entity_description[type].parameter[j].name)));
      fprintf(fp, "& %s ",
              translateUnitsToTex(PRINTABLE_NULL(entity_description[type].parameter[j].unit)));
      fprintf(fp, "& %s & ",
              makeTexSafeString(PRINTABLE_NULL(type_name[entity_description[type].parameter[j].type-1])));
    }
    switch (entity_description[type].parameter[j].type) {
    case IS_DOUBLE:
      if (latex_form && entity_description[type].parameter[j].number==0)
        fprintf(fp, " 0.0");
      else 
        fprintf(fp, "  %.15g", entity_description[type].parameter[j].number);
      break;
    case IS_LONG:
      if (latex_form && entity_description[type].parameter[j].integer==0)
        fprintf(fp, " \\verb|0|");
      else      
        fprintf(fp, "  %-15ld", entity_description[type].parameter[j].integer);
      break;
    case IS_STRING:
      fprintf(fp, "  %-15s", 
              PRINTABLE_NULL(entity_description[type].parameter[j].string));
      break;
    default:
      fprintf(stdout, "Invalid parameter type for %s item of %s\n",
              PRINTABLE_NULL(entity_description[type].parameter[j].name),
              entity_name[type]);
      fflush(stdout);
      exit(1);
    }
    if (latex_form) {
      char *ptr0, buffer[1024];
      strcpy(buffer, entity_description[type].parameter[j].description);
      if (strlen(ptr0 = buffer)) {
        /* don't need splitting of strings since the p tabular code 
           is used.
           */
        /*
        while (strlen(ptr0)>20) {
          ptr1 = ptr0+20;
          while (*ptr1 && *ptr1!=' ')
            ptr1++;
          if (!*ptr1)
            break;
          *ptr1 = 0;
          fprintf(fp, " & %s ", makeTexSafeString(ptr0));
          fprintf(fp, "\\\\ \n & & & ");
          ptr0 = ptr1+1;
          texLines++;
        }
        */
        /* add to lines counter based on estimate of wrap-around lines
           in the latex parbox. 28 is approximate number of char. in 2 in */
        texLines += strlen(ptr0) / 28;
        if (*ptr0) {
          fprintf(fp, " & %s ", 
                  makeTexSafeString(ptr0));
          fprintf(fp, " \\\\ \\hline \n");
          texLines++;
        }
      } else {
        fprintf(fp, " & \\\\ \\hline \n");
        texLines++;
      }
    }
    else 
      fprintf(fp, "  %s\n", entity_description[type].parameter[j].description);
  }
  if (latex_form) {
    char physicsFile[1024];
    
    fprintf(fp, "\\end{tabular}\n\n");

    sprintf(physicsFile, "%s.tex", entity_name[type]);
    str_tolower(physicsFile);
    if (fexists(physicsFile)) {
      fprintf(stderr, "Including file %s\n", physicsFile);
      fprintf(fp, "\\vspace*{0.5in}\n\\input{%s}\n", physicsFile);
    }
  }
}

#include <memory.h>

void initialize_structures(RUN *run_conditions, VARY *run_control, ERRORVAL *error_control, CORRECTION *correct, 
                           BEAM *beam, OUTPUT_FILES *output_data, OPTIMIZATION_DATA *optimize,
                           CHROM_CORRECTION *chrom_corr_data, TUNE_CORRECTION *tune_corr_data,
                           ELEMENT_LINKS *links)
{
    if (run_conditions)
        memset((void*)run_conditions, 0, sizeof(*run_conditions));
    if (run_control)
        memset((void*)run_control, 0, sizeof(*run_control));
    run_control->bunch_frequency = 2856e6;
    if (error_control)
        memset((void*)error_control, 0, sizeof(*error_control));
    if (correct)
        memset((void*)correct, 0, sizeof(*correct));
    correct->mode = correct->method = -1;
    if (beam)
        memset((void*)beam, 0, sizeof(*beam));
    if (output_data)
        memset((void*)output_data, 0, sizeof(*output_data));
    if (optimize)
        memset((void*)optimize, 0, sizeof(*optimize));
    if (chrom_corr_data)
        memset((void*)chrom_corr_data, 0, sizeof(*chrom_corr_data));
    if (tune_corr_data)
        memset((void*)tune_corr_data, 0, sizeof(*tune_corr_data));
    if (links)
        memset((void*)links, 0, sizeof(*links));
    }

char *makeTexSafeString(char *source)
{
  static char buffer[1024];
  long index = 0;
  buffer[0] = 0;
  while (*source) {
    if (*source=='_' || *source=='^' || *source=='{' || *source=='}') {
      buffer[index++] = '\\';
      buffer[index++] = *source++;
    }
    else if  (*source=='<' || *source=='>' || *source=='|') {
      buffer[index++] = '$';
      buffer[index++] = *source++;
      buffer[index++] = '$';
    }
    else
      buffer[index++] = *source++;
  }
  buffer[index] = 0;
  return buffer;
}

char *translateUnitsToTex(char *source)
{
  static char buffer[1024];
  char *ptr, *ptr1;

  if (strlen(source)==0)
    return source;
  buffer[0] = '$';
  buffer[1] = 0;
  ptr = source;
  while (*ptr) {
    if (*ptr=='$') {
      switch (*(ptr1=ptr+1)) {
      case 'a':
      case 'b':
      case 'A':
      case 'B':
        while (*ptr1 && *ptr1!='$')
          ptr1++;
        if (*ptr1 && (*(ptr1+1)=='n' || *(ptr1+1)=='N')) {
          if (*(ptr+1)=='a' || *(ptr+1)=='A')
            strcat(buffer, "^{");
          else
            strcat(buffer, "_{");
          strncat(buffer, ptr+2, (ptr1-1)-(ptr+2)+1);
          strcat(buffer, "}");
          ptr = ptr1+1;
        } else
          strncat(buffer, ptr, 1);
        break;
      default:
        fprintf(stdout, "Unrecognized $ sequence: %s\n", ptr);
        fflush(stdout);
        strncat(buffer, ptr, 1);
        break;
      }
    } else 
      strncat(buffer, ptr, 1);
    ptr++;
  }
  strcat(buffer, "$");
  return buffer;
}

void free_beamdata(BEAM *beam)
{
  if (beam->particle)
    free_zarray_2d((void**)beam->particle, beam->n_particle, 7);
  if (beam->accepted)
    free_zarray_2d((void**)beam->accepted, beam->n_particle, 7);
  if (beam->original && beam->original!=beam->particle)
    free_zarray_2d((void**)beam->original, beam->n_original, 7);
  beam->particle = beam->accepted = beam->original = NULL;
  beam->n_original = beam->n_to_track = beam->n_accepted = beam->p0 = beam->n_saved = beam->n_particle = 0;
}  


