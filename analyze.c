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
#include "correctDefs.h"

#define CMATRIX_OFFSET 0
#define RMATRIX_OFFSET 6
#define X_BETA_OFFSET RMATRIX_OFFSET+36
#define Y_BETA_OFFSET X_BETA_OFFSET+2
#define DETR_OFFSET Y_BETA_OFFSET+2
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
    {"nux", "&column name=nux, symbol=\"$gn$r$bx$n\", type=double &end"},
    {"betay", "&column name=betay, units=m, symbol=\"$gb$r$by$n\", type=double &end"},
    {"nuy", "&column name=nuy, symbol=\"$gn$r$by$n\", type=double &end"},
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
    process_namelist(&analyze_map, nltext);
    print_namelist(stdout, &analyze_map);

    /* check for data errors */
    if (!output)
        bomb("no output filename specified", NULL);
    if (n_points!=2 && n_points!=4)
        bomb("n_points must be either 2 or 4", NULL);

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
    static double **coord;
    static double *data;
    static double *offset;
    static long initialized = 0;
    static MATRIX *R, *Rc;
    double p_central, beta, tune;
    long n_track, n_trpoint, i, j, effort, index;
    double sin_phi, cos_phi, det;
    static double *orbit_p, *orbit_m;
    TRAJECTORY *clorb=NULL;

    log_entry("do_transport_analysis");
        
    if (center_on_orbit && !orbit)
        bomb("you've asked to center the analysis on the closed orbit, but you didn't issue a closed_orbit command", NULL);

    if (!initialized) {
        coord = (double**)zarray_2d(sizeof(**coord), 1+6*n_points, 7);
        data  = (double*)tmalloc(sizeof(*data)*SDDS_analyze.layout.n_columns);
        offset = (double*)tmalloc(sizeof(*offset)*6);
        m_alloc(&R , 6, 6);
        m_alloc(&Rc, 6, 6);
        orbit_p = tmalloc(sizeof(*orbit_p)*6);
        orbit_m = tmalloc(sizeof(*orbit_m)*6);
        initialized = 1;
        }
    if (center_on_orbit)
        clorb = tmalloc(sizeof(*clorb)*(beamline->n_elems+1));

    n_track = n_trpoint = n_points*6+1;
    for (j=0; j<7; j++)
        for (i=0; i<n_track; i++)
            coord[i][j] = (orbit && center_on_orbit?orbit[j]:0);


    /* particles 0 and 1 are for d/dx */
    offset[0] = delta_x ;
    coord[0][0] += delta_x ;
    coord[1][0] -= delta_x ;
    /* particles 2 and 3 are for d/dxp */
    offset[1] = delta_xp;
    coord[2][1] += delta_xp;
    coord[3][1] -= delta_xp;
    /* particles 4 and 5 are for d/dy */
    offset[2] = delta_y ;
    coord[4][2] += delta_y ;
    coord[5][2] -= delta_y ;
    /* particles 6 and 7 are for d/dy */
    offset[3] = delta_y ;
    coord[6][3] += delta_y ;
    coord[7][3] -= delta_y ;
    /* particles 8 and 9 are for d/ds */
    offset[4] = delta_s ;
    coord[8][4] += delta_s ;
    coord[9][4] -= delta_s ;
    /* particles 10 and 11 are for d/dp */
    offset[5] = delta_dp;
    coord[10][5] += delta_dp;
    coord[11][5] -= delta_dp;
    if (n_points==4) {
        /* particles 12 and 13 are for d/dx */
        coord[12][0] += 3*delta_x ;
        coord[13][0] -= 3*delta_x ;
        /* particles 14 and 15 are for d/dxp */
        coord[14][1] += 3*delta_xp;
        coord[15][1] -= 3*delta_xp;
        /* particles 16 and 17 are for d/dy */
        coord[16][2] += 3*delta_y ;
        coord[17][2] -= 3*delta_y ;
        /* particles 18 and 19 are for d/dy */
        coord[18][3] += 3*delta_y ;
        coord[19][3] -= 3*delta_y ;
        /* particles 20 and 21 are for d/ds */
        coord[20][4] += 3*delta_s ;
        coord[21][4] -= 3*delta_s ;
        /* particles 22 and 23 are for d/dp */
        coord[22][5] += 3*delta_dp;
        coord[23][5] -= 3*delta_dp;
        }
    /* particle n_track-1 is the reference particle */

    effort = 0;
    p_central = run->p_central;
    if (do_tracking(coord, &n_trpoint, &effort, beamline, &p_central, 
                    NULL, NULL, NULL, NULL, run, control->i_step, SILENT_RUNNING, control->n_passes, 0,
                    NULL, NULL, NULL)!=n_track) {
        fputs("warning: particle(s) lost during transport analysis--continuing with next step", stdout);
        log_exit("do_transport_analysis");
        return;
        }

    if (verbosity>0){
        if (orbit && center_on_orbit) {
            fprintf(stdout, "closed orbit: \n");
            fflush(stdout);
            for (i=0; i<6; i++)
                fprintf(stdout, "%15.8e ", orbit[i]);
            fputc('\n', stdout);
            }
        fprintf(stdout, "final coordinates of refence particle: \n");
        fflush(stdout);
        for (i=0; i<6; i++)
            fprintf(stdout, "%15.8e ", coord[n_track-1][i]);
        fputc('\n', stdout);
        fflush(stdout);
        if (verbosity>1) {
          for (i=0; i<n_track-1; i++) {
            fprintf(stdout, "Particle %ld end coordinates:\n", i);
            for (j=0; j<4; j++)
              fprintf(stdout, "%e%c", coord[i][j], j==3?'\n':' ');
          }
          fflush(stdout);
        }
      }
    

    for (i=0; i<6; i++) 
        data[i+CMATRIX_OFFSET] = coord[n_track-1][i];

    for (i=0; i<6; i++) {
        /* Compute R[i][j] --> data[6*i+j] */
        /* i indexes the dependent quantity */
        for (j=0; j<6; j++) {
            /* j indexes the initial coordinate value */
            if (!offset[j]) 
                R->a[i][j] = data[RMATRIX_OFFSET+6*i+j] = (i==j?1:0);
            else {
                if (n_points==2) 
                    R->a[i][j] = data[RMATRIX_OFFSET+6*i+j] = (coord[2*j][i]-coord[2*j+1][i])/(2*offset[j]);
                else
                    R->a[i][j] = data[RMATRIX_OFFSET+6*i+j] = 
                        (27*(coord[2*j][i]-coord[2*j+1][i])-(coord[2*j+12][i]-coord[2*j+13][i]))/(48*offset[j]);
                }
            }
        }
    /* find NUx */
    if (fabs(cos_phi = (R->a[0][0]+R->a[1][1])/2)>=1) {
        fprintf(stdout, "warning: beamline unstable for x plane\n");
        fflush(stdout);
        beta = 0;
        tune = 0;
        }
    else {
        beta    = fabs(R->a[0][1]/sin(acos(cos_phi)));
        sin_phi = R->a[0][1]/beta;
        if ((tune = atan2(sin_phi, cos_phi)/PIx2)<0)
            tune += 1;
        }
    data[X_BETA_OFFSET  ] = beta;
    data[X_BETA_OFFSET+1] = tune;
    /* find NUy */
    if (fabs(cos_phi = (R->a[2][2]+R->a[3][3])/2)>=1) {
        fprintf(stdout, "warning: beamline unstable for y plane\n");
        fflush(stdout);
        beta = 0;
        tune = 0;
        }
    else {
        beta    = fabs(R->a[2][3]/sin(acos(cos_phi)));
        sin_phi = R->a[2][3]/beta;
        if ((tune = atan2(sin_phi, cos_phi)/PIx2)<0)
            tune += 1;
        }
    data[Y_BETA_OFFSET  ] = beta;
    data[Y_BETA_OFFSET+1] = tune;
    m_copy(Rc, R);
    data[DETR_OFFSET] = m_det(R);

    /* compute etax and etax' */
    if ((det = (2 - R->a[0][0] - R->a[1][1]))<=0) {
        fprintf(stdout, "error: beamline unstable for x plane--can't match dispersion functions\n");
        fflush(stdout);
        det = 1e-6;
        }
    data[X_ETA_OFFSET  ] = ((1-R->a[1][1])*R->a[0][5]+R->a[0][1]*R->a[1][5])/det;
    data[X_ETA_OFFSET+1] = (R->a[1][0]*R->a[0][5] + (1-R->a[0][0])*R->a[1][5])/det;

    /* compute etay and etay' */
    if ((det = (2 - R->a[2][2] - R->a[3][3]))<=0) {
        fprintf(stdout, "error: beamline unstable for y plane--can't match dispersion functions\n");
        fflush(stdout);
        det = 1e-6;
        }
    data[Y_ETA_OFFSET  ] = ((1-R->a[3][3])*R->a[2][5]+R->a[2][3]*R->a[3][5])/det;
    data[Y_ETA_OFFSET+1] = (R->a[3][2]*R->a[2][5] + (1-R->a[2][2])*R->a[3][5])/det;

    if (delta_dp && center_on_orbit) {
        if (!beamline->matrix)
            beamline->matrix = full_matrix(&(beamline->elem), run, 1);
        find_closed_orbit(clorb, 1e-12, 20, beamline, beamline->matrix, run, delta_dp, 1, 0, NULL, 1.0,
                          NULL);
        for (i=0; i<6; i++)
            orbit_p[i] = clorb[0].centroid[i];
        find_closed_orbit(clorb, 1e-12, 20, beamline, beamline->matrix, run, -delta_dp, 1, 0, NULL, 1.0,
                          NULL);
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
        fprintf(stdout, "Unable to start SDDS table (do_transport_analysis)");
        fflush(stdout);
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
        exit(1);
        }
    if (!SDDS_SetParameters(&SDDS_analyze, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, 
                            IP_STEP, control->i_step, IP_SUBSTEP, control->i_vary, -1)) {
        SDDS_SetError("Unable to set SDDS parameter values (do_transport_analysis)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
    for (i=0; i<SDDS_ColumnCount(&SDDS_analyze); i++)
        if (!SDDS_SetRowValues(&SDDS_analyze, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, 0,
                               i, data[i], -1)) {
            fprintf(stdout, "Unable to set SDDS column %s (do_transport_analysis)\n", 
                    analysis_column[i].name);
            fflush(stdout);
            SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
            exit(1);
            }
    if (!SDDS_WriteTable(&SDDS_analyze)) {
        fprintf(stdout, "Unable to write SDDS table (do_transport_analysis)");
        fflush(stdout);
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
        exit(1);
        }
    if (!SDDS_DoFSync(&SDDS_analyze))
      fprintf(stdout, "Warning: problem fsync'ing map analysis output file\n");
    
    if (verbosity>0) {
        for (i=0; i<6; i++) {
            fprintf(stdout, "R%ld: ", i+1);
            fflush(stdout);
            for (j=0; j<6; j++) 
                fprintf(stdout, "%20.13e ", R->a[i][j]);
                fflush(stdout);
            fputc('\n', stdout);
            }
        fprintf(stdout, "horizontal:   tune = %20.13e  beta = %20.13e  eta = %20.13e  eta' = %20.13e\n",
            data[X_BETA_OFFSET+1], data[X_BETA_OFFSET], data[X_ETA_OFFSET], data[X_ETA_OFFSET+1]);
        fflush(stdout);
        fprintf(stdout, "vertical  :   tune = %20.13e  beta = %20.13e  eta = %20.13e  eta' = %20.13e\n",
            data[Y_BETA_OFFSET+1], data[Y_BETA_OFFSET], data[Y_ETA_OFFSET], data[Y_ETA_OFFSET+1]);
        fflush(stdout);
        fprintf(stdout, "determinant of R = 1 + %20.13e\n", data[DETR_OFFSET]-1);
        fflush(stdout);
        fprintf(stdout, "dispersion functions from closed-orbit calculations:\nx: %e m    %e\ny: %e m    %e\n",
            data[CLORB_ETA_OFFSET  ], data[CLORB_ETA_OFFSET+1],
            data[CLORB_ETA_OFFSET+2], data[CLORB_ETA_OFFSET+3]);
        fflush(stdout);
        }

    log_exit("do_transport_analysis");
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
    }

