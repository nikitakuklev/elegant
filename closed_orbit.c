/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: closed_orbit.c
 * purpose: computation closed orbits
 *
 * Michael Borland, 1992
 */
#include "mdb.h"
#include "track.h"

static long SDDS_clorb_initialized = 0;
static SDDS_TABLE SDDS_clorb;
static long clorb_count = 0;

#define IC_S 0
#define IC_X 1
#define IC_Y 2
#define IC_ELEMENT 3
#define IC_OCCURENCE 4
#define IC_TYPE 5
#define N_COLUMNS 6
static SDDS_DEFINITION column_definition[N_COLUMNS] = {
    {"s", "&column name=s, units=m, type=double, description=\"Distance\" &end"},
    {"x", "&column name=x, units=m, type=double, description=\"Horizontal position\" &end"},
    {"y", "&column name=y, units=m, type=double, description=\"Vertical position\" &end"},
    {"ElementName", "&column name=ElementName, type=string, description=\"Element name\", format_string=%10s &end"},
    {"ElementOccurence", 
         "&column name=ElementOccurence, type=long, description=\"Occurence of element\", format_string=%6ld &end"},
    {"ElementType", "&column name=ElementType, type=string, description=\"Element-type name\", format_string=%10s &end"},
    } ;

#define IP_STEP 0
#define IP_XERROR 1
#define IP_YERROR 2
#define IP_DELTA 3
#define N_PARAMETERS 4
static SDDS_DEFINITION parameter_definition[N_PARAMETERS] = {
    {"Step", "&parameter name=Step, type=long, description=\"Simulation step\" &end"},
    {"xError", "&parameter name=xError, type=double, units=m, description=\"Horizontal closed orbit convergence error\" &end"},
    {"yError", "&parameter name=yError, type=double, units=m, description=\"Vertical closed orbit convergence error\" &end"},
    {"delta", "&parameter name=delta, symbol=\"$gd$r\", type=double, description=\"Fractional energy offset of closed orbit\" &end"},
    } ;

#include "closed_orbit.h"

static TRAJECTORY *clorb = NULL;


void setup_closed_orbit(NAMELIST_TEXT *nltext, RUN *run, LINE_LIST *beamline)
{

    log_entry("setup_closed_orbit");

    if (clorb)
        free(clorb);
    clorb = tmalloc(sizeof(*clorb)*(beamline->n_elems+1));

    /* process namelist input */
    set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
    set_print_namelist_flags(0);
    process_namelist(&closed_orbit, nltext);
    print_namelist(stdout, &closed_orbit);

    if (output)
        output = compose_filename(output, run->rootname);
    if (closed_orbit_accuracy<=0)
        bomb("closed_orbit_accuracy <= 0", NULL);
    if (closed_orbit_iterations<1)
        bomb("closed_orbit_iterations < 1", NULL);
    if (iteration_fraction<0 || iteration_fraction>1)
        bomb("iteration_fraction must be on [0, 1]", NULL);
    if (output) {
        SDDS_ElegantOutputSetup(&SDDS_clorb, output, SDDS_BINARY, 1, "closed orbit", 
                                run->runfile, run->lattice, parameter_definition, N_PARAMETERS,
                                column_definition, N_COLUMNS, "setup_closed_orbit", 
                                SDDS_EOS_NEWFILE|SDDS_EOS_COMPLETE);
        SDDS_clorb_initialized = 1;
        }

    log_exit("setup_closed_orbit");
    }

long run_closed_orbit(RUN *run, LINE_LIST *beamline, double *starting_coord, BEAM *beam, long do_output)
{
    double dp, deviation[4];
    long i, bad_orbit;
    VMATRIX *M;
    
    if (!starting_coord)
        bomb("starting_coord array is NULL (run_closed_orbit)", NULL);

    if (start_from_centroid) {
        if (!beam)
            bomb("no beam present for closed-orbit calculation starting from centroid", NULL);
        compute_centroids(starting_coord, beam->particle, beam->n_to_track);
        dp = starting_coord[5];
        }
    else
        dp = 0;

    if (!clorb)
        bomb("TRAJECTORY array for closed orbit not allocated (run_closed_orbit)", NULL);
    beamline->closed_orbit = clorb;

    if (beamline->elem_recirc)
        M = full_matrix(beamline->elem_recirc, run, 1);
    else 
        M = full_matrix(&(beamline->elem), run, 1);

    bad_orbit = !find_closed_orbit(clorb, closed_orbit_accuracy, closed_orbit_iterations, beamline, M, 
                                   run, dp, start_from_recirc, fixed_length, 
                                   (start_from_centroid?starting_coord:NULL), iteration_fraction,
                                   deviation);
    free_matrices(M); tfree(M); M = NULL;
    
    /* return closed orbit at the beginning of the ring */
    for (i=0; i<6; i++)
        starting_coord[i] = clorb[0].centroid[i];

    /* do output, if required */
    if (verbosity && !bad_orbit) {
        fprintf(stdout, "closed orbit: \n");
        fflush(stdout);
        for (i=0; i<7; i++)
            fprintf(stdout, "%.8e ", starting_coord[i]);
            fflush(stdout);
        fputc('\n', stdout);
        }

    if (do_output && SDDS_clorb_initialized)
        dump_closed_orbit(clorb, beamline->n_elems, clorb_count++, deviation);

    return !bad_orbit;
    }
        
void finish_clorb_output(void)
{
    log_entry("finish_clorb_output");
    if (SDDS_IsActive(&SDDS_clorb) && !SDDS_Terminate(&SDDS_clorb)) {
        SDDS_SetError("Problem terminating SDDS output (finish_clorb_output)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
    SDDS_clorb_initialized = clorb_count = 0;
    log_exit("finish_clorb_output");
    }


void dump_closed_orbit(TRAJECTORY *traj, long n_elems, long step, double *deviation)
{
    long i, n, occurence, row;
    double position;
    char *name;

    log_entry("dump_closed_orbit");

    if (!SDDS_clorb_initialized)
        return;

    /* count number of trajectory elements actually used */
    for (i=1; i<n_elems+1; i++) {
        if (!traj[i].elem)
            break;
        }
    n = i;

    if (!SDDS_StartTable(&SDDS_clorb, n)) {
        SDDS_SetError("Unable to start SDDS table (dump_closed_orbit)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }

    if (!SDDS_SetParameters(&SDDS_clorb, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE,
                            IP_STEP, step,
                            IP_XERROR, deviation[0],
                            IP_YERROR, deviation[2],
                            IP_DELTA, traj[0].centroid[5],
                            -1)) {
        SDDS_SetError("Unable to set SDDS parameters (dump_closed_orbit)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }

    position = traj[1].elem->end_pos - 
              (entity_description[traj[1].elem->type].flags&HAS_LENGTH ?
                  *((double*)traj[1].elem->p_elem):0.0);
    name = "_BEG_";
    occurence = 1;

    for (i=row=0; i<n; i++) {
        if (i) {
            position = traj[i].elem->end_pos;
            name = traj[i].elem->name;
            occurence = traj[i].elem->occurence;
            }
        if (output_monitors_only &&
            (i==0 ||
             !(traj[i].elem->type==T_MONI || traj[i].elem->type==T_HMON || traj[i].elem->type==T_VMON)))
            continue;
        if (!SDDS_SetRowValues(&SDDS_clorb, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, row++,
                               IC_S, position, IC_X, traj[i].centroid[0], IC_Y, traj[i].centroid[2],
                               IC_ELEMENT, name, IC_OCCURENCE, occurence, 
                               IC_TYPE, i==0?"MARK":entity_name[traj[i].elem->type], -1)) {
            fprintf(stdout, "Unable to set row %ld values (dump_closed_orbit)\n", i);
            fflush(stdout);
            SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
            exit(1);
            }
        }

    if (!SDDS_WriteTable(&SDDS_clorb)) {
        SDDS_SetError("Unable to write closed orbit data (dump_closed_orbit)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
    if (!SDDS_EraseData(&SDDS_clorb)) {
        SDDS_SetError("Unable to erase closed orbit data (dump_closed_orbit)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }

    log_exit("dump_closed_orbit");
    }

