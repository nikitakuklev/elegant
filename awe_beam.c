/* Copyright 1994 by Michael Borland and Argonne National Laboratory,
 * all rights reserved.
 */
/* file: awe_beam.c
 * purpose: Do tracking for beams from awe dumps.
 *          See file awe_beam.nl for input parameters.
 *
 * Michael Borland, 1989
 */
#include "mdb.h"
#include "track.h"
#include "awe.h"
#include "awe_beam.h"
#include "match_string.h"

long read_dump_for_list(DUMP *dump, DUMP_SPECS *dump_specs, FILE **fp, long n_files);

static long input_type_code;
#define ELEGANT_BEAM 0
#define MASK_BEAM 1
#define SPIFFE_BEAM 2
#define N_AWE_INPUT_TYPES 3
static char *input_type_name[N_AWE_INPUT_TYPES] = {
    "elegant", "mask", "spiffe"
    } ;


/* indices of data from MASK output in arrays returned by get_particles()
 */
static long i_x1, i_x2, i_beta1, i_beta2, i_time, i_gamma, i_p;

/* additional indices for spiffe beam data */
static long i_r, i_pz, i_pr;

/* additional indices of data from elegant output in arrays returned by get_particles()
 */
static long i_transverse[4];

static long n_input_files = 0;
static char **input_file = NULL;

#ifdef VAX_VMS
#define isnan(x) 0
#define isinf(x) 0
#endif

void setup_awe_beam(
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

    log_entry("setup_awe_beam");

    if (!beam)
        bomb("BEAM pointer is null in setup_awe_beam", NULL);
    if (!nltext)
        bomb("NAMELIST_TEXT pointer is null in setup_awe_beam", NULL);
    if (!run)
        bomb("RUN pointer is null in setup_awe_beam", NULL);
    if (!control)
        bomb("VARY pointer is null in setup_awe_beam", NULL);
    if (!errcon)
        bomb("ERROR pointer is null in setup_awe_beam", NULL);
    if (!output)
        bomb("OUTPUT_FILES pointer is null in setup_awe_beam", NULL);
    if (!beamline)
        bomb("beamline pointer is null in setup_awe_beam", NULL);
    if (!run->runfile || !run->lattice)
        bomb("null runfile or lattice pointer in RUN structure in setup_awe_beam", NULL);

    if (!initial_call) {
        get_particles(NULL, NULL, 0, 0, 0);
        }
    else
        initial_call = 0;

    /* process namelist input */
    set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
    set_print_namelist_flags(0);
    process_namelist(&awe_beam, nltext);
    print_namelist(stderr, &awe_beam);
    fflush(stderr);

    /* check for validity of namelist inputs */
    if (input==NULL)
        bomb("no input file given in namelist awe_beam", NULL);
    n_input_files = 0;
    while (ptr=get_token(input)) {
        input_file = trealloc(input_file, sizeof(*input_file)*(n_input_files+1));
        cp_str(input_file+n_input_files, ptr);
        if (!fexists(input_file[n_input_files])) {
            fprintf(stderr, "error: input file %s does not exist\n", input_file[n_input_files]);
            exit(1);
            }
        n_input_files++;
        }
    if ((input_type_code=match_string(input_type, input_type_name, N_AWE_INPUT_TYPES, 0))<0)
        bomb("unknown awe input type", NULL);
    if (input_type_code==MASK_BEAM && n_particles_per_ring<=0)
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

    log_exit("setup_awe_beam");
    }

long new_awe_beam(
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

    log_entry("new_awe_beam");

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
             * even though it may be spread over several dumps.
             * Read in all the particles from the input file and allocate arrays
             * for storing initial coordinates and coordinates of accepted 
             * particles. */
            if (beam->n_original==0) {
                /* there is no bunch in memory */
                if (!(beam->n_original=get_particles(&beam->original, input_file, n_input_files, prebunched, 0))) 
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
            if (beam->n_original=get_particles(&beam->original, input_file, n_input_files, prebunched, n_dumps_to_skip)) { 
                n_dumps_to_skip = 0;    /* use the user's parameter only the first time */
                beam->particle = (double**)zarray_2d(sizeof(double), beam->n_original, 7);
                if (run->acceptance) 
                    beam->accepted = (double**)zarray_2d(sizeof(double), beam->n_original, 7);
                }
            else { 
                log_exit("new_awe_beam");
                return(0);
                }
            new_particle_data = 1;
            }
        }

    t_offset = (control->bunch_frequency?(control->i_step-1)/control->bunch_frequency:0);

    p_central = beam->p0_original = run->p_central;
    /* Create the initial distribution from the beam->original particle
     * data, then track */
    if (new_particle_data || generate_new_bunch || 
        ((input_type_code==MASK_BEAM || input_type_code==SPIFFE_BEAM) && !one_random_bunch && !(flags&TRACK_PREVIOUS_BUNCH))) {
        if (input_type_code==MASK_BEAM) {
            /* Since the MASK particles represent rings of charge,
             * it is necessary to approximate each ring by a number of 
             * point particles. */
            if (!beam->original)
                bomb("beam->original array is NULL (new_awe_beam-1)", NULL);
            if (!beam->particle)
                bomb("beam->particle array is NULL (new_awe_beam-1)", NULL);
            for (i=i_store=0; i<beam->n_original; i+=sample_interval, i_store++) {
                if (!beam->original[i]) {
                    fprintf(stderr, "error: beam->original[%ld] is NULL (new_awe_beam-1)\n", i);;
                    exit(1);
                    }
                if (sample_fraction!=1 && random_2(1)>sample_fraction) {
                    i_store--;
                    continue;
                    }
                if (i_gamma>=0)
                    p      = sqrt(sqr(gamma=(pti=beam->original[i])[i_gamma]) - 1);
                else
                    gamma  = sqrt(sqr(p=(pti=beam->original[i])[i_p]) + 1);
                if (p_lower && (p_lower>p || p_upper<p)) {
                    i_store--;
                    continue;
                    }
                r      = pti[i_x2];
                path   = (t_offset + pti[i_time])*c_mks*(beta = p/gamma) ;
                delta  = (p-p_central)/p_central;
                rp     = pti[i_beta2]/pti[i_beta1];
                theta  = PIx2*random_2(1);
                sin_theta = sin(theta);
                cos_theta = cos(theta);
                if (!beam->particle[i_store]) {
                    fprintf(stderr, "error: beam->particle[%ld] is NULL (new_awe_beam-1)\n", i_store);
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
        else if (input_type_code==SPIFFE_BEAM) {
            /* Since the spiffe particles represent rings of charge,
             * it is necessary to approximate each ring by a number of 
             * point particles. */
            
            if (!beam->original)
                bomb("beam->original array is NULL (new_awe_beam-2)", NULL);
            if (!beam->particle)
                bomb("beam->particle array is NULL (new_awe_beam-2)", NULL);
            for (i=i_store=0; i<beam->n_original; i+=sample_interval, i_store++) {
                if (!beam->original[i]) {
                    fprintf(stderr, "error: beam->original[%ld] is NULL (new_awe_beam-2)\n", i);;
                    exit(1);
                    }
                if (sample_fraction!=1 && random_2(1)>sample_fraction) {
                    i_store--;
                    continue;
                    }
                pti = beam->original[i];
                p = sqrt(sqr(pti[i_pz]) + sqr(pti[i_pr]));
                gamma = sqrt(sqr(p)+1);
                if (p_lower && (p_lower>p || p_upper<p)) {
                    i_store--;
                    continue;
                    }
                r      = pti[i_r];
                path   = (t_offset + pti[i_time])*c_mks*(beta = p/gamma) ;
                delta  = (p-p_central)/p_central;
                rp     = pti[i_pr]/pti[i_pz];
                theta  = PIx2*random_2(1);
                sin_theta = sin(theta);
                cos_theta = cos(theta);
                if (!beam->particle[i_store]) {
                    fprintf(stderr, "error: beam->particle[%ld] is NULL (new_awe_beam-2)\n", i_store);
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
            if (!beam->original)
                bomb("beam->original array is NULL (new_awe_beam)", NULL);
            if (!beam->particle)
                bomb("beam->particle array is NULL (new_awe_beam)", NULL);
            for (i=i_store=0; i<beam->n_original; i_store++,i+=sample_interval) {
                if (sample_fraction!=1 && random_2(1)>sample_fraction) {
                    i_store--;
                    continue;
                    }
                if (!beam->original[i]) {
                    fprintf(stderr, "error: beam->original[%ld] is NULL (new_awe_beam.2)\n", i);
                    exit(1);
                    }
                if (i_gamma>=0)
                    p     = sqrt(sqr(gamma = beam->original[i][i_gamma])-1);
                else
                    gamma = sqrt(sqr(p=beam->original[i][i_p])+1);
                if (p_lower && (p_lower>p || p_upper<p)) {
                    i_store--;
                    continue;
                    }
                if (!beam->particle[i_store]) {
                    fprintf(stderr, "error: beam->particle[%ld] is NULL (new_awe_beam.2)\n", i_store);
                    exit(1);
                    }
                for (j=0; j<4; j++)
                    beam->particle[i_store][j] = beam->original[i][i_transverse[j]];
                /* convert time to path-length */
                beam->particle[i_store][4] = (t_offset+beam->original[i][i_time])*c_mks*p/gamma;
                /* convert energy to dp/p */
                beam->particle[i_store][5] = (p-p_central)/p_central;
                beam->particle[i_store][6] = particleID++;
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
            bomb("beam->original is NULL (new_awe_beam.3)", NULL);
        if (!beam->particle)
            bomb("beam->particle is NULL (new_awe_beam.3)", NULL);
        for (i=0; i<beam->n_saved; i++) {
            if (!beam->original[i]) {
                fprintf(stderr, "error: beam->original[%ld] is NULL (new_awe_beam.3)\n", i);
                exit(1);
                }
            if (!beam->particle[i]) {
                fprintf(stderr, "error: beam->particle[%ld] is NULL (new_awe_beam.3)\n", i);
                exit(1);
                }
            for (j=0; j<6; j++) 
                beam->particle[i][j] = beam->original[i][j];
            p = p_central*(1+beam->particle[i][5]);
            beta = p/sqrt(p*p+1);
            beam->particle[i][4] += (control->bunch_frequency?1./control->bunch_frequency:0)*beta*c_mks;
            beam->particle[i][6] = particleID++;
            }
        new_particle_data = 0;
        }

     if (new_particle_data && (one_random_bunch || reuse_bunch)) {
        /* Copy the new "initial" data into original[] in case it needs to be reused,  
           but only if the stuff already in original[] is not going to be needed again.
         */

        if (!beam->original)
            bomb("beam->original is NULL (new_awe_beam.4)", NULL);
        if (!beam->particle)
            bomb("beam->particle is NULL (new_awe_beam.4)", NULL);
        for (i=0; i<beam->n_to_track; i++) {
            if (!beam->original[i]) {
                fprintf(stderr, "error: beam->original[%ld] is NULL (new_awe_beam.4)\n", i);
                exit(1);
                }
            if (!beam->particle[i]) {
                fprintf(stderr, "error: beam->particle[%ld] is NULL (new_awe_beam.4)\n", i);
                exit(1);
                }
            for (j=0; j<6; j++)
                beam->original[i][j] = beam->particle[i][j];
            }
        beam->n_saved = beam->n_to_track;
        new_particle_data = 0;
        }

    log_exit("new_awe_beam");
    return(beam->n_to_track);
    }

long get_particles(double ***particle, char **input_file, long n_input_files, long one_dump, long n_skip)
{
    static FILE **fp_in=NULL;
    static DUMP_SPECS *input_specs = NULL;
    static DUMP *input_dump = NULL;
    long i, np_max, np, np_new, ifile, dump_rejected;
    static long last_nif = 0;
    double **data;

    log_entry("get_particles");

#ifdef DEBUG
    fprintf(stderr, "heap verification: get_particles.1\n");
    malloc_verify();
#endif

    if (particle==NULL) {
        /* reset for reading again */
        for (ifile=0; ifile<last_nif; ifile++) 
            if (fp_in[ifile]) {
                fclose(fp_in[ifile]);
                fp_in[ifile] = NULL;
                }
        last_nif = 0;
        log_exit("get_particles");
        return(0);
        }
        
    if (!last_nif) {
#ifdef DEBUG
        fprintf(stderr, "heap verification: get_particles.2\n");
        malloc_verify();
#endif
        last_nif = n_input_files;
        fp_in = trealloc(fp_in, sizeof(*fp_in)*n_input_files);
        input_specs = trealloc(input_specs, sizeof(*input_specs)*n_input_files);
        input_dump  = trealloc(input_dump, sizeof(*input_dump)*n_input_files);
        for (ifile=0; ifile<n_input_files; ifile++) {
            fp_in[ifile] = fopen_e(input_file[ifile], "r", 0);
            read_dump_specs(input_specs+ifile, fp_in[ifile]);
            if ((input_type_code==MASK_BEAM && input_specs[ifile].n_quantities<8) ||
                (input_type_code==ELEGANT_BEAM && input_specs[ifile].n_quantities<7) ||
                (input_type_code==SPIFFE_BEAM && input_specs[ifile].n_quantities<3) )
                bomb("too few quantities in dump file.", NULL);
            fprintf(stderr, "dump specifications:\nrun identifier = %s\n%ld quantities per point:\n",
                   input_specs[ifile].run_identifier, input_specs[ifile].n_quantities);
            for (i=0; i<input_specs[ifile].n_quantities; i++) 
                fprintf(stderr, "    %s in %s in %s precision\n", 
                       input_specs[ifile].quantity_name[i], input_specs[ifile].quantity_unit[i],
                       (input_specs[ifile].quantity_precision[i]==DOUBLE_PRECISION?
                        "double":"single"));
            if (input_type_code==MASK_BEAM) {
                /* find indices of data */
                i_x1    = match_string(   "x1", input_specs[ifile].quantity_name, input_specs[ifile].n_quantities, EXACT_MATCH);
                i_x2    = match_string(   "x2", input_specs[ifile].quantity_name, input_specs[ifile].n_quantities, EXACT_MATCH);
                i_beta1 = match_string("beta1", input_specs[ifile].quantity_name, input_specs[ifile].n_quantities, EXACT_MATCH);
                i_beta2 = match_string("beta2", input_specs[ifile].quantity_name, input_specs[ifile].n_quantities, EXACT_MATCH);
                i_time  = match_string(    "t", input_specs[ifile].quantity_name, input_specs[ifile].n_quantities, EXACT_MATCH);
                i_gamma = match_string("gamma", input_specs[ifile].quantity_name, input_specs[ifile].n_quantities, EXACT_MATCH);
                i_p     = match_string(    "p", input_specs[ifile].quantity_name, input_specs[ifile].n_quantities, EXACT_MATCH);
                if (i_x1<0 || i_x2<0 || i_beta1<0 || i_beta2<0 || i_time<0 ||
                    (i_gamma<0 && i_p<0))
                    bomb("necessary data quantities not present in input file", NULL); 
                }
            else if (input_type_code==SPIFFE_BEAM) {
                /* find indices of data */
                i_r    = match_string( "r", input_specs[ifile].quantity_name, input_specs[ifile].n_quantities, EXACT_MATCH);
                i_pr   = match_string("pr", input_specs[ifile].quantity_name, input_specs[ifile].n_quantities, EXACT_MATCH);
                i_pz   = match_string("pz", input_specs[ifile].quantity_name, input_specs[ifile].n_quantities, EXACT_MATCH);
                i_time = match_string( "t", input_specs[ifile].quantity_name, input_specs[ifile].n_quantities, EXACT_MATCH);
                if (i_r<0 || i_pr<0 || i_pz<0 || i_time<0) 
                    bomb("necessary data quantities (r, pr, pz, t) not present in input file", NULL); 
                }
            else {
                i_transverse[0] = match_string( "x", 
                                               input_specs[ifile].quantity_name, input_specs[ifile].n_quantities, EXACT_MATCH);
                i_transverse[1] = match_string("xp", 
                                               input_specs[ifile].quantity_name, input_specs[ifile].n_quantities, EXACT_MATCH);
                i_transverse[2] = match_string( "y", 
                                               input_specs[ifile].quantity_name, input_specs[ifile].n_quantities, EXACT_MATCH);
                i_transverse[3] = match_string("yp", 
                                               input_specs[ifile].quantity_name, input_specs[ifile].n_quantities, EXACT_MATCH);
                i_time  = match_string(    "t", input_specs[ifile].quantity_name, input_specs[ifile].n_quantities, EXACT_MATCH);
                i_gamma = match_string("gamma", input_specs[ifile].quantity_name, input_specs[ifile].n_quantities, EXACT_MATCH);
                i_p     = match_string(    "p", input_specs[ifile].quantity_name, input_specs[ifile].n_quantities, EXACT_MATCH);
                if (i_transverse[0]<0 || i_transverse[1]<0 || i_transverse[2]<0 ||
                    i_transverse[3]<0 || i_time<0 || (i_gamma<0 && i_p<0))
                    bomb("necessary data quantities not present in input file", NULL); 
                }
            }
        }

    np_max = np = 0;
    data = NULL;
    while (read_dump_for_list(input_dump, input_specs, fp_in, n_input_files)) {
        if (one_dump && n_skip>0) {
            n_skip--;
            continue;
            }
#ifdef DEBUG
        fprintf(stderr, "heap verification: get_particles.3\n");
        malloc_verify();
#endif
        dump_rejected = 0;
        for (ifile=0; ifile<n_input_files; ifile++) {
            if (dump_selection_string) {
                if (!input_dump[ifile].dump_label || !wild_match(input_dump[ifile].dump_label, dump_selection_string)) {
                    dump_rejected = 1;
                    break;
                    }
                }
            if ((np_new=np+input_dump[ifile].n_points)>np_max) {
#ifdef DEBUG
                fprintf(stderr, "heap verification: get_particles.4\n");
                malloc_verify();
#endif
                /* must reallocate to get more space */
                np_max = np + 2*input_dump[ifile].n_points;
                data = trealloc(data, np_max*sizeof(*data));
                }
            if (input_specs[ifile].n_quantities>=6)
                for (i=np; i<np_new; i++)
                    data[i] = input_dump[ifile].point[i-np];
            else 
                /* will want to use this array for 6-d phase space coordinates latter on */
                for (i=np; i<np_new; i++) 
                    data[i] = trealloc(input_dump[ifile].point[i-np], sizeof(**data)*7);
            np = np_new;
#ifdef DEBUG
            fprintf(stderr, "heap verification: get_particles.5\n");
            malloc_verify();
#endif
            tfree(input_dump[ifile].point);
            input_dump[ifile].point = NULL;
#ifdef DEBUG
            fprintf(stderr, "heap verification: get_particles.6\n");
            malloc_verify();
#endif
#ifdef DEBUG
            fprintf(stderr, "heap verification: get_particles.7\n");
            malloc_verify();
#endif
            }
        if (one_dump && !dump_rejected)
            break;
        }
    
    fprintf(stderr, "a total of %ld data points were read\n\n", np);
    *particle = data;
#ifdef DEBUG
    fprintf(stderr, "heap verification: get_particles.8\n");
    malloc_verify();
#endif
    log_exit("get_particles");
    return(np);
    }            

long read_dump_for_list(DUMP *dump, DUMP_SPECS *dump_specs, FILE **fp, long n_files)
{
    long ifile, data_seen;

    log_entry("read_dump_for_list");

    if (!dump)
        bomb("dump pointer is NULL (read_dump_for_list)", NULL);
    if (!dump_specs)
        bomb("dump_specs pointer is NULL (read_dump_for_list)", NULL);
    if (!fp)
        bomb("file pointer is NULL (read_dump_for_list)", NULL);
    
    data_seen = 0;
    for (ifile=0; ifile<n_files; ifile++) {
        if (fp[ifile]) {
            if (read_dump(dump+ifile, dump_specs+ifile, fp[ifile])) {
                data_seen = 1;
                }
            else {
                fclose(fp[ifile]);
                fp[ifile] = NULL;
                }
            }
        }
    log_exit("read_dump_for_list");
    return(data_seen);
    }
