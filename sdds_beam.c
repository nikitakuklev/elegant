/* Copyright 1994 by Michael Borland and Argonne National Laboratory,
 * all rights reserved.
 */
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

static long n_input_files = 0;
static char **input_file = NULL;
static SDDS_TABLE *SDDS_input = NULL;
static long *input_initialized = NULL;

#ifdef VAX_VMS
#define isnan(x) 0
#define isinf(x) 0
#endif

/* SDDS routines will be asked to deliver data in the order given in
 * this string, which matches the order of the #define's 
 */
static char *spiffe_columns = "r pr pz t";
#define ISC_R 0
#define ISC_PR 1
#define ISC_PZ 2
#define ISC_T 3

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
long check_sdds_beam_column(SDDS_TABLE *SDDS_table, char *name, char *units);

void setup_sdds_beam(
    BEAM *beam,
    NAMELIST_TEXT *nltext,
    RUN *run, 
    VARY *control,
    ERROR *errcon,
    OPTIM_VARIABLES *optim,
    OUTPUT_FILES *output,
    LINE_LIST *beamline,
    long n_elements
    )
{
    char t[1024], *ptr;
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
    print_namelist(stderr, &sdds_beam);
    fflush(stderr);

    /* check for validity of namelist inputs */
    if (input==NULL)
        bomb("no input file given in namelist sdds_beam", NULL);
    if ((selection_parameter && !selection_string) || (!selection_parameter && selection_string))
        bomb("must specify selection_parameter and selection_string together or not at all", NULL);

    n_input_files = 0;
    while (ptr=get_token(input)) {
        input_file = trealloc(input_file, sizeof(*input_file)*(n_input_files+1));
        SDDS_input = trealloc(SDDS_input, sizeof(*SDDS_input)*(n_input_files+1));
        input_initialized = trealloc(input_initialized, sizeof(*input_initialized)*(n_input_files+1));
        cp_str(input_file+n_input_files, ptr);
        SDDS_ZeroMemory(SDDS_input+n_input_files, sizeof(*SDDS_input));
        input_initialized[n_input_files] = 0;
        if (!fexists(input_file[n_input_files])) {
            fprintf(stderr, "error: input file %s does not exist\n", input_file[n_input_files]);
            exit(1);
            }
        n_input_files++;
        }
    if ((input_type_code=match_string(input_type, input_type_name, N_SDDS_INPUT_TYPES, 0))<0)
        bomb("unknown sdds input type", NULL);
    if (input_type_code==SPIFFE_BEAM && n_particles_per_ring<=0)
        bomb("n_particles_per_ring is invalid", NULL);
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

    beam->original = beam->particle = beam->accepted = NULL;
    beam->n_original = beam->n_to_track = beam->n_accepted = 0;

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
    double r, theta, rp, delta;
    double sin_theta, cos_theta, path, *pti;

    log_entry("new_sdds_beam");

    if (flags&TRACK_PREVIOUS_BUNCH) {
        if (beam->n_original==0)
            bomb("can't retrack with previous bunch--there isn't one!", NULL);
        if (n_particles_per_ring!=1) {
            fputs("Warning: can't do retracking with previous bunch when n_particles_per_ring!=1\n", stderr);
            fputs("Will use a new bunch generated from previously read data.\n", stderr);
            generate_new_bunch = 1;
            }
        else
            generate_new_bunch = 0;
        }
    else {
        if (!prebunched) {
            /* The beam in the input file is to be treated as a single bunch,
             * even though it may be spread over several tables.
             * Read in all the particles from the input file and allocate arrays
             * for storing initial coordinates and coordinates of accepted 
             * particles. */
            if (beam->n_original==0) {
                /* there is no bunch in memory */
                if (!(beam->n_original=get_sdds_particles(&beam->original, prebunched, 0))) 
                    bomb("no particles in input file", NULL);
                beam->particle = (double**)zarray_2d(sizeof(double), beam->n_original, 7);
                if (run->acceptance)
                    beam->accepted = (double**)zarray_2d(sizeof(double), beam->n_original, 7);
                new_particle_data = 1;
                }
            else
                new_particle_data = 0;
            }
        else {
            /* Each dump in the input file is to be treated as a separate
             * bunch.  Read in the next dump for tracking.
             */
            if (beam->n_original) {
                if (beam->particle)
                    free_zarray_2d((void**)beam->particle, beam->n_original, 7);
                if (beam->accepted)
                    free_zarray_2d((void**)beam->accepted, beam->n_original, 7);
                if (beam->original)
                    free_zarray_2d((void**)beam->original, beam->n_original, 7);
                }
            beam->particle = beam->accepted = beam->original = NULL;
            if ((beam->n_original=get_sdds_particles(&beam->original, prebunched, n_tables_to_skip))>=0) { 
                n_tables_to_skip = 0;    /* use the user's parameter only the first time */
                beam->particle = (double**)zarray_2d(sizeof(double), beam->n_original, 7);
                if (run->acceptance) 
                    beam->accepted = (double**)zarray_2d(sizeof(double), beam->n_original, 7);
                }
            else { 
                log_exit("new_sdds_beam");
                return(-1);
                }
            new_particle_data = 1;
            }
        }

    t_offset = (control->bunch_frequency?(control->i_step-1)/control->bunch_frequency:0);

    p_central = beam->p0_original = run->p_central;
    /* Create the initial distribution from the beam->original particle data */
    if (new_particle_data || generate_new_bunch || 
        (input_type_code==SPIFFE_BEAM && !one_random_bunch && !(flags&TRACK_PREVIOUS_BUNCH))) {
        if (input_type_code==SPIFFE_BEAM) {
            /* Since the spiffe particles represent rings of charge,
             * it is necessary to approximate each ring by a number of 
             * point particles. */
            if (!beam->original)
                bomb("beam->original array is NULL (new_sdds_beam-2)", NULL);
            if (!beam->particle)
                bomb("beam->particle array is NULL (new_sdds_beam-2)", NULL);
            for (i=i_store=0; i<beam->n_original; i+=sample_interval, i_store++) {
                if (!beam->original[i]) {
                    fprintf(stderr, "error: beam->original[%ld] is NULL (new_sdds_beam-2)\n", i);;
                    exit(1);
                    }
                if (sample_fraction!=1 && random_2(1)>sample_fraction) {
                    i_store--;
                    continue;
                    }
                pti = beam->original[i];
                p = sqrt(sqr(pti[ISC_PZ]) + sqr(pti[ISC_PR]));
                gamma = sqrt(sqr(p)+1);
                if (p_lower && (p_lower>p || p_upper<p)) {
                    i_store--;
                    continue;
                    }
                r      = pti[ISC_R];
                path   = (t_offset + pti[ISC_T])*c_mks*(beta = p/gamma) ;
                delta  = (p-p_central)/p_central;
                rp     = pti[ISC_PR]/pti[ISC_PZ];
                theta  = PIx2*random_2(1);
                sin_theta = sin(theta);
                cos_theta = cos(theta);
                if (!beam->particle[i_store]) {
                    fprintf(stderr, "error: beam->particle[%ld] is NULL (new_sdds_beam-2)\n", i_store);
                    exit(1);
                    }
                beam->particle[i_store][0] = r*cos_theta;
                beam->particle[i_store][1] = rp*cos_theta;
                beam->particle[i_store][2] = r*sin_theta;
                beam->particle[i_store][3] = rp*sin_theta;
                beam->particle[i_store][4] = path;
                beam->particle[i_store][5] = delta;
                beam->particle[i_store][6] = particleID++;
                }
            beam->n_to_track = i_store;
#ifndef VAX_VMS
            for (i_store=0; i_store<beam->n_to_track; i_store++) {
                for (i=0; i<6; i++) {
                    if (!beam->particle[i_store]) {
                        fprintf(stderr, "error: beam->particle[%ld] is NULL\n", i_store);
                        exit(1);
                        }
                    if (isnan(beam->particle[i_store][i]) || isinf(beam->particle[i_store][i])) {
                        fprintf(stderr, "error: NaN or Infinity detected in initial particle data, coordinate %ld\n", i);
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
            if (center_arrival_time)
                adjust_arrival_time_data(beam->particle, beam->n_to_track, p_central);
            }
        else {
            /* In this case, the data is for point particles already,
             * so I just copy the data for the most part. */
            if (!beam->original && beam->n_original)
                bomb("beam->original array is NULL (new_sdds_beam)", NULL);
            if (!beam->particle)
                bomb("beam->particle array is NULL (new_sdds_beam)", NULL);
            for (i=i_store=0; i<beam->n_original; i_store++,i+=sample_interval) {
                if (sample_fraction!=1 && random_2(1)>sample_fraction) {
                    i_store--;
                    continue;
                    }
                if (!beam->original[i]) {
                    fprintf(stderr, "error: beam->original[%ld] is NULL (new_sdds_beam.2)\n", i);
                    exit(1);
                    }
                gamma = sqrt(sqr(p=beam->original[i][IEC_P])+1);
                if (p_lower && (p_lower>p || p_upper<p)) {
                    i_store--;
                    continue;
                    }
                if (!beam->particle[i_store]) {
                    fprintf(stderr, "error: beam->particle[%ld] is NULL (new_sdds_beam.2)\n", i_store);
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
                        fprintf(stderr, "error: NaN or Infinity detected in initial particle data, coordinate %ld\n", i);
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
            if (center_arrival_time)
                adjust_arrival_time_data(beam->particle, beam->n_to_track, p_central);
            }
        }
    else {
        /* use (x, x', y, x', s, dp/p) saved in beam->original[] */
        beam->n_to_track = beam->n_saved;
        if (!beam->original)
            bomb("beam->original is NULL (new_sdds_beam.3)", NULL);
        if (!beam->particle)
            bomb("beam->particle is NULL (new_sdds_beam.3)", NULL);
        for (i=0; i<beam->n_saved; i++) {
            if (!beam->original[i]) {
                fprintf(stderr, "error: beam->original[%ld] is NULL (new_sdds_beam.3)\n", i);
                exit(1);
                }
            if (!beam->particle[i]) {
                fprintf(stderr, "error: beam->particle[%ld] is NULL (new_sdds_beam.3)\n", i);
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

     if (new_particle_data && (one_random_bunch || reuse_bunch)) {
        /* Copy the new "initial" data into original[] in case it needs to be reused,  
           but only if the stuff already in original[] is not going to be needed again.
         */

        if (!beam->original)
            bomb("beam->original is NULL (new_sdds_beam.4)", NULL);
        if (!beam->particle)
            bomb("beam->particle is NULL (new_sdds_beam.4)", NULL);
        for (i=0; i<beam->n_to_track; i++) {
            if (!beam->original[i]) {
                fprintf(stderr, "error: beam->original[%ld] is NULL (new_sdds_beam.4)\n", i);
                exit(1);
                }
            if (!beam->particle[i]) {
                fprintf(stderr, "error: beam->particle[%ld] is NULL (new_sdds_beam.4)\n", i);
                exit(1);
                }
            for (j=0; j<7; j++)
                beam->original[i][j] = beam->particle[i][j];
            }
        beam->n_saved = beam->n_to_track;
        new_particle_data = 0;
        }

    log_exit("new_sdds_beam");
    return(beam->n_to_track);
    }

/* get_sdds_particles reads data from all of the input files and puts it into the particle array as
 * follows:
 *     for spiffe input:   (*particle)[i] = (r, pr, pz, t) for ith particle
 *     for elegant input:  (*particle)[i] = (x, xp, y, yp, t, p) for ith particle
 */
long get_sdds_particles(double ***particle, long one_dump, long n_skip)
{
    long i, np_max, np, np_new, rows, ifile, dump_rejected;
    long files_initialized, retval, data_seen;
    double **data, **new_data;
    static char s[200];
    long indexID = -1;

    log_entry("get_sdds_particles");

    if (particle==NULL) {
        /* reset for reading again */
        for (ifile=0; ifile<n_input_files; ifile++)  {
            if (input_initialized[ifile] && !SDDS_Terminate(SDDS_input+ifile)) {
                sprintf(s, "Problem terminating SDDS beam input from file %s", input_file[ifile]);
                SDDS_SetError(s);
                SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
                }
            input_initialized[ifile] = 0;
            }
        log_exit("get_sdds_particles");
        return(-1);
        }

    files_initialized = 0;
    for (ifile=0; ifile<n_input_files; ifile++) {
        if (input_initialized[ifile])
            continue;
        files_initialized++;
        input_initialized[ifile] = 1;
        if (!SDDS_InitializeInput(SDDS_input+ifile, input_file[ifile])) {
            sprintf(s, "Problem opening beam input file %s", input_file[ifile]);
            SDDS_SetError(s);
            SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
            }
        if (selection_parameter) {
            if ((i=SDDS_GetParameterIndex(SDDS_input+ifile, selection_parameter))<0)
                fprintf(stderr, "warning: SDDS beam file %s does not contain the selection parameter %s\n",
                        selection_parameter);
            if (SDDS_GetParameterType(SDDS_input+ifile, i)!=SDDS_STRING) {
              sprintf(s, "SDDS beam file %s contains parameter %s, but parameter is not a string", 
                      selection_parameter);
              SDDS_SetError(s);
              SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
            }
          }
        if (input_type_code==SPIFFE_BEAM) {
            if (!check_sdds_beam_column(SDDS_input+ifile, "r", "m") || 
                !check_sdds_beam_column(SDDS_input+ifile, "pr", "m$be$nc") ||
                !check_sdds_beam_column(SDDS_input+ifile, "pz", "m$be$nc") ||
                !check_sdds_beam_column(SDDS_input+ifile, "t", "s")) {
                fprintf(stderr, 
                        "necessary data quantities (r, pr, pz, t) have the wrong units or are not present in %s", 
                        input_file[ifile]);
                exit(1);
                }
            }
        else {
            if (!check_sdds_beam_column(SDDS_input+ifile, "x", "m") ||
                !check_sdds_beam_column(SDDS_input+ifile, "y", "m") ||
                !check_sdds_beam_column(SDDS_input+ifile, "xp", NULL) ||
                !check_sdds_beam_column(SDDS_input+ifile, "yp", NULL) ||
                !check_sdds_beam_column(SDDS_input+ifile, "p", "m$be$nc") ||
                !check_sdds_beam_column(SDDS_input+ifile, "t", "s")) {
                fprintf(stderr, 
                        "necessary data quantities (x, x', y, y', t, p) have the wrong units or are not present in %s", 
                        input_file[ifile]);
                exit(1);
                }
            }
        }
    if (files_initialized && files_initialized!=n_input_files)
        bomb("the number of SDDS beam files initialized isn't equal to the number of files in use", NULL);
        
    np_max = np = 0;
    data = NULL;
    data_seen = 1;
    while (data_seen) {
        data_seen = 0;
        for (ifile=0; ifile<n_input_files; ifile++) {
            if ((retval=SDDS_ReadTable(SDDS_input+ifile))==0)
                SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
            else if (retval==-1) 
                continue;
            data_seen = 1;
            }
        if (!data_seen)
            break;
        if (one_dump && n_skip>0) {
            n_skip--;
            continue;
            }
        dump_rejected = 0;
        for (ifile=0; ifile<n_input_files; ifile++) {
            if (selection_parameter) {
                char *value;
                if (!SDDS_GetParameter(SDDS_input+ifile, selection_parameter, &value)) {
                    sprintf(s, "Problem getting value of parameter %s from file %s", selection_parameter, input_file[ifile]);
                    SDDS_SetError(s);
                    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
                    }
                if (!wild_match(value, selection_string)) {
                    dump_rejected = 1;
                    break;
                    }
                }
            SDDS_SetColumnFlags(SDDS_input+ifile, 0);
            if (input_type_code==SPIFFE_BEAM && 
                !SDDS_SetColumnsOfInterest(SDDS_input+ifile, SDDS_NAMES_STRING, spiffe_columns)) {
                sprintf(s, "Problem setting columns of interest for file %s", input_file[ifile]);
                SDDS_SetError(s);
                SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
                }
            else if (input_type_code!=SPIFFE_BEAM &&
                     !SDDS_SetColumnsOfInterest(SDDS_input+ifile, SDDS_NAMES_STRING, elegant_columns)) {
                sprintf(s, "Problem setting columns of interest for file %s", input_file[ifile]);
                SDDS_SetError(s);
                SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
                }
            if ((rows = SDDS_CountRowsOfInterest(SDDS_input+ifile))<=0) {
                if (rows==-1) {
                    sprintf(s, "Problem counting rows of interest for file %s", input_file[ifile]);
                    SDDS_SetError(s);
                    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
                    }
                break;
                }
            if ((np_new=np+rows)>np_max) {
                /* must reallocate to get more space */
                np_max = np + 2*rows;
                data = trealloc(data, np_max*sizeof(*data));
                }
            if (!(new_data=SDDS_GetCastMatrixOfRows(SDDS_input+ifile, &i, SDDS_DOUBLE))) {
                sprintf(s, "Problem getting matrix of rows for file %s", input_file[ifile]);
                SDDS_SetError(s);
                SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
                }
            if (i!=rows) {
                sprintf(s, "Row count mismatch for file %s", input_file[ifile]);
                SDDS_SetError(s);
                SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
                }
            for (i=np; i<np_new; i++)
                /* will want to use this storage for 6D+1 phase space latter */
                data[i] = trealloc(new_data[i-np], sizeof(**new_data)*7);
            if ((indexID=SDDS_GetColumnIndex(SDDS_input+ifile, "particleID"))>=0) {
                double *index;
                if (!(index=SDDS_GetColumnInDoubles(SDDS_input+ifile, "particleID"))) {
                    sprintf(s, "Problem reading particleID column for file %s", input_file[ifile]);
                    SDDS_SetError(s);
                    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
                    }
                for (i=np; i<np_new; i++)
                    data[i][6] = index[i-np];
                }
            else if (input_type_code!=SPIFFE_BEAM)
                for (i=np; i<np_new; i++)
                    data[i][6] = particleID++;
            free(new_data);
            np = np_new;
            }
        if (one_dump && !dump_rejected)
            break;
        }
    
    fprintf(stderr, "a total of %ld data points were read\n\n", np);
    *particle = data;

    log_exit("get_sdds_particles");
    return(np);
    }            

long check_sdds_beam_column(SDDS_TABLE *SDDS_table, char *name, char *units)
{
    char *units1;
    if (SDDS_GetColumnIndex(SDDS_table, name)<0)
        return(0);
    if (SDDS_GetColumnInformation(SDDS_table, "units", &units1, SDDS_GET_BY_NAME, name)!=SDDS_STRING) {
        SDDS_SetError("units field of column has wrong data type!");
        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
        }
    if (!units) {
        if (!units1)
            return(1);
        if (SDDS_StringIsBlank(units1)) {
            free(units1);
            return(1);
            }
        return(0);
        }
    if (!units1)
        return(0);
    if (strcmp(units, units1)==0) {
        free(units1);
        return(1);
        }
    free(units1);
    return(0);
    }

void adjust_arrival_time_data(double **coord, long np, double Po)
{
    long ip;
    double P, beta, time;
    double tave;

    if (np<1)
        return;

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
