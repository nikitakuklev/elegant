/* Copyright 1994 by Michael Borland and Argonne National Laboratory,
 * all rights reserved.
 */
/* program: elegant
 * purpose: accelerator simulation
 * Michael Borland, 1989-1994
 */
#include "mdb.h"
#include "track.h"
#include "elegant.h"
#include "scan.h"
#include <ctype.h>
#include "match_string.h"
#include <signal.h>
#include <time.h>
#ifdef UNIX
#include <malloc.h>
#endif

#include "chromDefs.h"
#include "correctDefs.h"
#include "tuneDefs.h"

void traceback_handler(int code);

#define DESCRIBE_INPUT 0
#define N_OPTIONS  1
char *option[N_OPTIONS] = {
    "describeinput",
        };
char *USAGE="elegant <inputfile> [-describeInput]\n\nProgram by Michael Borland. (This is version 13.10, February 1999.)";

char *GREETING="This is elegant version 13.10, by Michael Borland. (This is version 13.10, February 1999.)";

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
#define N_COMMANDS      36

char *command[N_COMMANDS] = {
    "run_setup", "run_control", "vary_element", "error_control", "error_element", "awe_beam", "bunched_beam",
    "correct", "matrix_output", "twiss_output", "track", "stop", 
    "optimization_setup", "optimize", "optimization_variable", "optimization_constraint",
    "optimization_covariable", "save_lattice", "rpn_expression", "trace", "chromaticity", "closed_orbit",
    "find_aperture", "analyze_map", "correct_tunes", "link_control", "link_elements",
    "steering_element", "amplification_factors", "print_dictionary", "floor_coordinates", "correction_matrix_output",
    "load_parameters", "sdds_beam", "subprocess", "fit_traces",
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
        } ;

void initialize_structures(RUN *run_conditions, VARY *run_control, ERROR *error_control, CORRECTION *correct, 
                           BEAM *beam, OUTPUT_FILES *output_data, OPTIMIZATION_DATA *optimize,
                           CHROM_CORRECTION *chrom_corr_data, TUNE_CORRECTION *tune_corr_data,
                           ELEMENT_LINKS *links);

#define NAMELIST_BUFLEN 4096

main(argc, argv)
int argc;
char **argv;
{
    LINE_LIST *beamline;        /* pointer to root of linked list */
    FILE *fp_in;
    char *inputfile;
    SCANNED_ARG *scanned;
    char s[NAMELIST_BUFLEN], *ptr;
    long i;
    RUN run_conditions;
    VARY run_control;
    CORRECTION correct;
    ERROR error_control;
    BEAM beam;
    OUTPUT_FILES output_data;
    OPTIMIZATION_DATA optimize;
    CHROM_CORRECTION chrom_corr_data;
    TUNE_CORRECTION tune_corr_data;
    ELEMENT_LINKS links;
    char *saved_lattice = NULL;
    long correction_setuped, run_setuped, run_controled, error_controled, beam_type;
    long do_chromatic_correction = 0, do_twiss_output = 0, fl_do_tune_correction = 0;
    long do_closed_orbit = 0, do_matrix_output = 0, do_response_output = 0, step;
    long last_default_order = 0, new_beam_flags, links_present, twiss_computed = 0;
    double *starting_coord;
    long namelists_read = 0, failed;

#if defined(UNIX)
    signal(SIGHUP, traceback_handler);
    signal(SIGINT, traceback_handler);
    signal(SIGQUIT, traceback_handler);
    signal(SIGILL, traceback_handler);
    signal(SIGABRT, traceback_handler);
    signal(SIGTRAP, traceback_handler);
    signal(SIGFPE, traceback_handler);
    signal(SIGBUS, traceback_handler);
    signal(SIGSEGV, traceback_handler);
#endif
    
    log_entry("main");
    compute_offsets();
    set_max_name_length(12);

#if defined(VAX_VMS) || defined(UNIX)
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
    if (argc<2 || argc>(2+N_OPTIONS)) {
        fprintf(stderr, "usage: %s\n", USAGE);
        link_date();
        exit(1);
        }
    
    fprintf(stderr, "%s\n", GREETING);
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
                show_namelists_fields(stderr, namelist_pointer, namelist_name, n_namelists);
                if (argc==2)
                    exit(0);
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
            print_namelist(stderr, &run_setup);
            
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
                fprintf(stderr, "clock-generated random_number_seed = %ld\n", random_number_seed);
                }
            
            /* seed random number generators.  Note that random_1 seeds random_2, and random_3.
             * random_1 is used for beamline errors.  random_2 is used for beam generation and 
             * random scraping/sampling/scattering.  random_3 is used for BPM noise.
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
            
            /* extract the root filename from the input filename */
            strcpy(s, inputfile);
            if (rootname==NULL) {
                clean_filename(s);
                if (ptr=strrchr(s, '.'))
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
            
            /* parse the lattice file and create the beamline */
            run_conditions.lattice = compose_filename(saved_lattice, rootname);
            beamline = get_beamline(lattice, use_beamline, p_central);
            fprintf(stderr, "length of beamline %s per pass: %21.15le m\n", beamline->name, beamline->revolution_length);
            lattice = saved_lattice;
            
            /* output the magnet layout in mpl format */
            if (magnets)
                output_magnets(magnets, lattice, beamline);

            delete_phase_references();    /* necessary for multi-step runs */
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
            correction_setuped = 1;
            correction_setup(&correct, &namelist_text, &run_conditions, beamline); 
            break;
          case SET_AWE_BEAM: 
            fprintf(stderr, "This program no longer supports awe-format files.\n");
            fprintf(stderr, "Use awe2sdds to convert your data files, and use\n");
            fprintf(stderr, "the sdds_beam command instead of awe_beam.\n");
            break;
          case SET_BUNCHED_BEAM:
            if (!run_setuped || !run_controled)
                bomb("run_setup and run_control must precede bunched_beam namelist", NULL);
            setup_bunched_beam(&beam, &namelist_text, &run_conditions, &run_control, &error_control, &optimize.variables,
                               &output_data, beamline, beamline->n_elems);
            setup_output(&output_data, &run_conditions, &run_control, &error_control, &optimize.variables, beamline);
            beam_type = SET_BUNCHED_BEAM;
            break;
          case SET_SDDS_BEAM: 
            if (!run_setuped || !run_controled)
                bomb("run_setup and run_control must precede sdds_beam namelist", NULL);
            setup_sdds_beam(&beam, &namelist_text, &run_conditions, &run_control, &error_control, 
                           &optimize.variables, &output_data, beamline, beamline->n_elems);
            setup_output(&output_data, &run_conditions, &run_control, &error_control, &optimize.variables, beamline);
            beam_type = SET_SDDS_BEAM;
            break;
          case TRACK:
            set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
            set_print_namelist_flags(0);
            process_namelist(&track, &namelist_text);
            print_namelist(stderr, &track);
            if (use_linear_chromatic_matrix && !twiss_computed)
              bomb("you must compute twiss parameters to do linear chromatic matrix tracking", NULL);
            
            if (beam_type==-1)
                bomb("beam must be defined prior to tracking", NULL);
            new_beam_flags = 0;
            while (vary_beamline(&run_control, &error_control, &run_conditions, beamline)) {
              fill_double_array(starting_coord, 6, 0.0);
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
                               beamline, &beam, &output_data, 0);
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
                                   run_control.i_step, 1) ) {
                  fputs("warning: orbit correction failed--continuing with next step\n", stderr);
                  continue;
                }
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
                fprintf(stderr, "Closed orbit not found---continuing to next step\n");
                continue;
              }
              if (do_twiss_output && 
                  !run_twiss_output(&run_conditions, beamline, starting_coord, 0) &&
                  !soft_failure) {
                fprintf(stderr, "Twiss parameters not defined---continuing to next step\n");
                continue;
              }
              if (do_response_output)
                run_response_output(&run_conditions, beamline, &correct, 0);
              run_matrix_output(&run_conditions, beamline);
              for (i=failed=0; (fl_do_tune_correction || do_chromaticity_correction) && i<correction_iterations; i++) {
                if (correction_iterations>1)
                  fprintf(stderr, "\nTune/chromaticity correction iteration %ld\n", i+1);
                if (fl_do_tune_correction) {
                  if (do_closed_orbit && 
                      !run_closed_orbit(&run_conditions, beamline, starting_coord, &beam, 0) &&
                      !soft_failure) {
                    fprintf(stderr, "Closed orbit not found---continuing to next step\n");
                    failed = 1;
                    break;
                  }
                  if (!do_tune_correction(&tune_corr_data, &run_conditions, beamline, starting_coord,
                                          run_control.i_step, i==correction_iterations-1) &&
                      !soft_failure) {
                    fprintf(stderr, "Tune correction failed---continuing to next step\n");
                    failed = 1;
                    break;
                  }
                }
                if (do_chromatic_correction) {
                  if (do_closed_orbit && 
                      !run_closed_orbit(&run_conditions, beamline, starting_coord, &beam, 0) &&
                      !soft_failure) {
                    fprintf(stderr, "Closed orbit not found---continuing to next step\n");
                    failed = 1;
                    break;
                  }
                  if (!do_chromaticity_correction(&chrom_corr_data, &run_conditions, beamline, starting_coord,
                                                  run_control.i_step, i==correction_iterations-1) &&
                      !soft_failure) {
                    fprintf(stderr, "Chromaticity correction failed---continuing to next step\n");
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
                fputs("warning: orbit correction failed--continuing with next step\n", stderr);
                continue;
              }
              if (do_closed_orbit && 
                  !run_closed_orbit(&run_conditions, beamline, starting_coord, &beam, 1) &&
                  !soft_failure) {
                fprintf(stderr, "Closed orbit not found---continuing to next step\n");
                continue;
              }
              if (do_twiss_output && 
                  !run_twiss_output(&run_conditions, beamline, starting_coord, 1) &&
                  !soft_failure) {
                fprintf(stderr, "Twiss parameters not defined---continuing to next step\n");
                continue;
              }
              if (do_response_output)
                run_response_output(&run_conditions, beamline, &correct, 1);
              if (center_on_orbit)
                center_beam_on_coords(beam.particle, beam.n_to_track, starting_coord, center_momentum_also);
              track_beam(&run_conditions, &run_control, &error_control, &optimize.variables, 
                         beamline, &beam, &output_data, 
                         use_linear_chromatic_matrix?LINEAR_CHROMATIC_MATRIX:0);
            }
            finish_output(&output_data, &run_conditions, &run_control, &error_control, &optimize.variables, 
                          beamline, beamline->n_elems, &beam);
            if (do_closed_orbit)
              finish_clorb_output();
            if (do_twiss_output)
              finish_twiss_output();
            if (do_response_output)
              finish_response_output();
#ifdef SUNOS4
            check_heap();
#endif
            fprintf(stderr, "Finished tracking.\n");
            /* reassert defaults for namelist run_setup */
            lattice = use_beamline = acceptance = centroid = sigma = final = output = rootname = losses = NULL;
            combine_bunch_statistics = 0;
            random_number_seed = 987654321;
            wrap_around = 1;
            default_order = 2;
            concat_order = 0;
            tracking_updates = 1;
            concat_order = print_statistics = p_central = 0;
#if defined(VAX_VMS) || defined(UNIX)
            report_stats(stderr, "statistics: ");
            fflush(stderr);
#endif
            run_setuped = run_controled = error_controled = correction_setuped = do_chromatic_correction =
              fl_do_tune_correction = do_closed_orbit = do_twiss_output = do_response_output = 0;
            break;
          case MATRIX_OUTPUT:
            if (!run_setuped)
                bomb("run_setup must precede matrix_output namelist", NULL);
            setup_matrix_output(&namelist_text, &run_conditions, beamline);
            break;
          case TWISS_OUTPUT:
            if (!run_setuped)
                bomb("run_setup must precede twiss_output namelist", NULL);
            setup_twiss_output(&namelist_text, &run_conditions, beamline, &do_twiss_output);
            if (!do_twiss_output) {
                twiss_computed = 1;
                run_twiss_output(&run_conditions, beamline, NULL, -1);
                finish_twiss_output();
                }
            break;
          case STOP:
            lorentz_report();
            finish_load_parameters();
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
                }
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
                fprintf(stderr, "warning: you've asked to do both closed-orbit calculation and orbit correction.\nThis may duplicate effort.\n");
            break;
          case FIND_APERTURE:
            setup_aperture_search(&namelist_text, &run_conditions, &run_control);
            while (vary_beamline(&run_control, &error_control, &run_conditions, beamline)) {
              fill_double_array(starting_coord, 6, 0.0);
              if (correct.mode!= -1) {
                if (!do_correction(&correct, &run_conditions, beamline, starting_coord, &beam, 
                                   run_control.i_step, 1) ) {
                  fputs("warning: orbit correction failed--continuing with next step\n", stderr);
                  continue;
                }
              }
              if (do_closed_orbit && 
                  !run_closed_orbit(&run_conditions, beamline, starting_coord, &beam, 0) &&
                  !soft_failure) {
                fprintf(stderr, "Closed orbit not found---continuing to next step\n");
                continue;
              }
              if (do_twiss_output && 
                  !run_twiss_output(&run_conditions, beamline, starting_coord, 0) &&
                  !soft_failure) {
                fprintf(stderr, "Twiss parameters not defined---continuing to next step\n");
                continue;
              }
              if (do_response_output)
                run_response_output(&run_conditions, beamline, &correct, 0);
              run_matrix_output(&run_conditions, beamline);
              for (i=failed=0; (fl_do_tune_correction || do_chromaticity_correction) && i<correction_iterations; i++) {
                if (correction_iterations>1)
                  fprintf(stderr, "\nTune/chromaticity correction iteration %ld\n", i+1);
                if (fl_do_tune_correction) {
                  if (do_closed_orbit && 
                      !run_closed_orbit(&run_conditions, beamline, starting_coord, &beam, 0) &&
                      !soft_failure) {
                    fprintf(stderr, "Closed orbit not found---continuing to next step\n");
                    failed = 1;
                    break;
                  }
                  if (!do_tune_correction(&tune_corr_data, &run_conditions, beamline, starting_coord,
                                          run_control.i_step, i==correction_iterations-1) &&
                      !soft_failure) {
                    fprintf(stderr, "Tune correction failed---continuing to next step\n");
                    failed = 1;
                    break;
                  }
                }
                if (do_chromatic_correction) {
                  if (do_closed_orbit && 
                      !run_closed_orbit(&run_conditions, beamline, starting_coord, &beam, 0) &&
                      !soft_failure) {
                    fprintf(stderr, "Closed orbit not found---continuing to next step\n");
                    failed = 1;
                    break;
                  }
                  if (!do_chromaticity_correction(&chrom_corr_data, &run_conditions, beamline, starting_coord,
                                                  run_control.i_step, i==correction_iterations-1) &&
                      !soft_failure) {
                    fprintf(stderr, "Chromaticity correction failed---continuing to next step\n");
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
                fputs("warning: orbit correction failed--continuing with next step\n", stderr);
                continue;
              }
              if (do_closed_orbit && 
                  !run_closed_orbit(&run_conditions, beamline, starting_coord, &beam, 1) &&
                  !soft_failure) {
                fprintf(stderr, "Closed orbit not found---continuing to next step\n");
                continue;
              }
              if (do_twiss_output && 
                  !run_twiss_output(&run_conditions, beamline, starting_coord, 1) &&
                  !soft_failure) {
                fprintf(stderr, "Twiss parameters not defined---continuing to next step\n");
                continue;
              }
              if (do_response_output)
                run_response_output(&run_conditions, beamline, &correct, 1);
              do_aperture_search(&run_conditions, &run_control, &error_control, beamline);
            }
            fprintf(stderr, "Finished all tracking steps.\n"); fflush(stderr);
            finish_aperture_search(&run_conditions, &run_control, &error_control, beamline);
            if (do_closed_orbit)
                finish_clorb_output();
#ifdef SUNOS4
            check_heap();
#endif
            fprintf(stderr, "Finished dynamic aperture search.\n");
            /* reassert defaults for namelist run_setup */
            lattice = use_beamline = acceptance = centroid = sigma = final = output = rootname = losses = NULL;
            combine_bunch_statistics = 0;
            random_number_seed = 987654321;
            wrap_around = 1;
            default_order = 2;
            concat_order = 0;
            tracking_updates = 1;
            concat_order = print_statistics = p_central = 0;
#if defined(VAX_VMS) || defined(UNIX)
            report_stats(stderr, "statistics: ");
            fflush(stderr);
#endif
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
                            fputs("warning: correction failed--continuing with next step", stderr);
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
                        fprintf(stderr, "\nTune/chromaticity correction iteration %ld\n", i+1);
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
                  fputs("warning: orbit correction failed--continuing with next step\n", stderr);
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
            fprintf(stderr, "Finished transport analysis.\n");
            /* reassert defaults for namelist run_setup */
            lattice = use_beamline = acceptance = centroid = sigma = final = output = rootname = losses = NULL;
            combine_bunch_statistics = 0;
            random_number_seed = 987654321;
            wrap_around = 1;
            default_order = 2;
            concat_order = 0;
            tracking_updates = 1;
            concat_order = print_statistics = p_central = 0;
#if defined(VAX_VMS) || defined(UNIX)
            report_stats(stderr, "statistics: ");
            fflush(stderr);
#endif
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
            print_namelist(stderr, &print_dictionary);
            do_print_dictionary(filename);
            break;
          case FLOOR_COORDINATES:
            output_floor_coordinates(&namelist_text, &run_conditions, beamline);
            break;
          case CORRECTION_MATRIX_OUTPUT:
            if (!run_setuped)
                bomb("run setup must precede correction_matrix_output namelist", NULL);
            setup_correction_matrix_output(&namelist_text, &run_conditions, beamline, &correct, &do_response_output);
            if (!do_response_output) {
                run_response_output(&run_conditions, beamline, &correct, -1);
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
            setup_load_parameters(&namelist_text, &run_conditions, beamline);
            break;
          case SUBPROCESS:
            run_subprocess(&namelist_text, &run_conditions);
            break;
          case FIT_TRACES:
            do_fit_trace_data(&namelist_text, &run_conditions, beamline);
            break;
          default:
            fprintf(stderr, "unknown namelist %s given.  Known namelists are:\n", namelist_text.group_name);
            for (i=0; i<N_COMMANDS; i++)
                fprintf(stderr, "%s\n", description[i]);
            exit(1);
            break;
            }
#ifdef SUNOS4
        check_heap();
#endif
        }
    fprintf(stderr, "End of input data encountered.\n"); fflush(stderr);
    lorentz_report();
    finish_load_parameters();
    log_exit("main");
    exit(0);
    }

double find_beam_p_central(char *input)
{
    bomb("the expand_for feature is disabled at this time", NULL);
    }

#ifdef SUNOS4
#include <malloc.h>

void check_heap() 
{
    struct mallinfo info;
    
    fprintf(stderr, "Performing memory heap verification..."); fflush(stderr);
    if (malloc_verify())
        fputs("okay.", stderr);
    else
        fputs("errors, detected.", stderr);
    fflush(stderr);
#if defined(VAX_VMS) || defined(UNIX)
    report_stats(stderr, "statistics: ");
    fflush(stderr);
#endif
    return;
    info.arena = info.ordblks = info.smblks = info.hblks = info.hblkhd = info.usmblks = 0;   
    info.fsmblks = info.uordblks = info.fordblks = info.keepcost = info.allocated = info.treeoverhead = 0;
    info = mallinfo();
    fprintf(stderr, "memory allocation information:\n");
    fprintf(stderr, "  total space in arena: %ld\n", info.arena);     
    fprintf(stderr, "  number of ordinary blocks: %ld\n", info.ordblks);   
    fprintf(stderr, "  number of small blocks: %ld\n", info.smblks);    
    fprintf(stderr, "  number of holding blocks: %ld\n", info.hblks);     
    fprintf(stderr, "  space in holding block headers: %ld\n", info.hblkhd);    
    fprintf(stderr, "  space in small blocks in use: %ld\n", info.usmblks);   
    fprintf(stderr, "  space in free small blocks: %ld\n", info.fsmblks);   
    fprintf(stderr, "  space in ordinary blocks in use: %ld\n", info.uordblks);  
    fprintf(stderr, "  space in free ordinary blocks: %ld\n", info.fordblks);  
    fprintf(stderr, "  cost of enabling keep option: %ld\n", info.keepcost);  
    fprintf(stderr, "  number of ordinary blocks allocated: %ld\n", info.allocated);
    fprintf(stderr, "  bytes used in maintaining the free tree: %ld\n", info.treeoverhead);
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

void do_print_dictionary(char *filename)
{
    FILE *fp;
    long i, j;

    if (!filename)
        bomb("filename invalid (do_print_dictionary)", NULL);
    fp = fopen_e(filename, "w", 0);
    for (i=1; i<N_TYPES; i++)
        print_dictionary_entry(fp, i);
    }

#define PRINTABLE_NULL(s) (s?s:"NULL")
void print_dictionary_entry(FILE *fp, long type)
{
    char *type_name[3] = {"double", "long", "STRING", };
    long j;
    fprintf(fp, "***** element type %s:\n", entity_name[type]);
    fprintf(fp, "%s\n", entity_text[type]);
    for (j=0; j<entity_description[type].n_params; j++) {
        fprintf(fp, "%20s  %20s  %10s", 
                PRINTABLE_NULL(entity_description[type].parameter[j].name), 
                PRINTABLE_NULL(entity_description[type].parameter[j].unit),
                PRINTABLE_NULL(type_name[entity_description[type].parameter[j].type-1]));
        switch (entity_description[type].parameter[j].type) {
          case IS_DOUBLE:
            fprintf(fp, "  %.15g\n", entity_description[type].parameter[j].number);
            break;
          case IS_LONG:
            fprintf(fp, "  %-15ld\n", entity_description[type].parameter[j].integer);
            break;
          case IS_STRING:
            fprintf(fp, "  %-15s\n", 
                    PRINTABLE_NULL(entity_description[type].parameter[j].string));
            break;
          default:
            fprintf(stderr, "Invalid parameter type for %s item of %s\n",
                    PRINTABLE_NULL(entity_description[type].parameter[j].name),
                    entity_name[type]);
            exit(1);
            }
        }
    }

#include <memory.h>

void initialize_structures(RUN *run_conditions, VARY *run_control, ERROR *error_control, CORRECTION *correct, 
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


