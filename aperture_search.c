/* Copyright 1994 by Michael Borland and Argonne National Laboratory,
 * all rights reserved.
 */
/* file: aperture_search.c
 * purpose: Do tracking to find machine aperture.
 *          See file aperture_search.nl for input parameters.
 *
 * Michael Borland, 1989
 */
#include "mdb.h"
#include "track.h"
#include "aperture_search.h"

#define IC_X 0
#define IC_Y 1
#define N_COLUMNS 2
static SDDS_DEFINITION column_definition[N_COLUMNS] = {
    {"x", "&column name=x, symbol=x, units=m, type=double &end"},
    {"y", "&column name=y, symbol=y, units=m, type=double &end"},
    } ;

#define IP_STEP 0
#define N_PARAMETERS 1
static SDDS_DEFINITION parameter_definition[N_PARAMETERS] = {
    {"Step", "&parameter name=Step, type=long, description=\"Simulation step\" &end"},
    } ;

static SDDS_TABLE SDDS_aperture;

#define N_SEARCH_MODES 2
static char *search_mode[N_SEARCH_MODES] = {
    "many-particle", "single-particle"
    } ;
static long mode_code = 0;

void setup_aperture_search(
    NAMELIST_TEXT *nltext,
    RUN *run,
    VARY *control
    )
{

    log_entry("setup_aperture_search");

    /* process namelist input */
    set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
    set_print_namelist_flags(0);
    process_namelist(&find_aperture, nltext);
    print_namelist(stderr, &find_aperture);

    /* check for data errors */
    if (!output)
        bomb("no output filename specified", NULL);
    if (xmin>=xmax)
        bomb("xmin >= xmax", NULL);
    if (ymin>=ymax)
        bomb("ymin >= ymax", NULL);
    if (nx<3)
        bomb("nx < 3", NULL);
    if (ny<2)
        bomb("ny < 2", NULL);
    if (n_splits && n_splits<1)
        bomb("n_splits is non-zero, but less than 1", NULL);
    if (n_splits) {
        if (split_fraction<=0 || split_fraction>=1)
            bomb("split_fraction must be greater than 0 and less than 1", NULL);
        if (desired_resolution<=0 || desired_resolution>=1)
            bomb("desired_resolution must be greater than 0 and less than 1", NULL);
        if ((desired_resolution *= (xmax-xmin))>(xmax-xmin)/(nx-1))
            bomb("desired_resolution is larger than coarse mesh", NULL);
        }
    if ((mode_code=match_string(mode, search_mode, N_SEARCH_MODES, 0))<0)
        bomb("unknown search mode", NULL);

    output = compose_filename(output, run->rootname);
    SDDS_ElegantOutputSetup(&SDDS_aperture, output, SDDS_BINARY, 1, 
                            (mode_code==0?"multi-particle aperture search":"single-particle aperture search"),
                            run->runfile, run->lattice, parameter_definition, N_PARAMETERS,
                            column_definition, N_COLUMNS, "setup_aperture_search", SDDS_EOS_NEWFILE);

    if (control->n_elements_to_vary) 
        if (!SDDS_DefineSimpleParameters(&SDDS_aperture, control->n_elements_to_vary,
                                         control->varied_quan_name, control->varied_quan_unit, SDDS_DOUBLE)) {
            SDDS_SetError("Unable to define additional SDDS parameters (setup_aperture_search)");
            SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
            }

    if (!SDDS_WriteLayout(&SDDS_aperture)) {
        SDDS_SetError("Unable to write SDDS layout for aperture search");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }

    if (boundary) {
        FILE *fp;
        boundary = compose_filename(boundary, run->rootname);
        fp = fopen_e(boundary, "w", 0);
        fputs("SDDS1\n&column name=x, units=m, type=double &end\n", fp);
        fputs("&column name=y, units=m, type=double &end\n", fp);
        fprintf(fp, "&parameter name=MplTitle, type=string, fixed_value=\"Aperture search boundary for run %s\", &end\n",
                run->runfile);
        fputs("&data mode=ascii, no_row_counts=1 &end\n", fp);
        fprintf(fp, "%e\t%e\n", xmin, ymin);
        fprintf(fp, "%e\t%e\n", xmin, ymax);
        fprintf(fp, "%e\t%e\n", xmax, ymax);
        fprintf(fp, "%e\t%e\n", xmax, ymin);
        fprintf(fp, "%e\t%e\n", xmin, ymin);
        fclose(fp);
        }

    log_exit("setup_aperture_search");
    }


long do_aperture_search(
    RUN *run,
    VARY *control,
    ERROR *errcon,
    LINE_LIST *beamline
    )
{    
    long retcode;

    log_entry("do_aperture_search");
    if (mode_code==0)
        retcode = do_aperture_search_mp(run, control, errcon, beamline);
    else
        retcode = do_aperture_search_sp(run, control, errcon, beamline);
    log_exit("do_aperture_search");
    return(retcode);
    }

/* many-particle search routine */

long do_aperture_search_mp(
    RUN *run,
    VARY *control,
    ERROR *errcon,
    LINE_LIST *beamline
    )
{
    double **coord, **accepted;
    double y, dx, dy;
    double **xy_left, **xy_right;
    long *found;
    long n_left, n_right, n_survived;
    double p_central;
    long n_trpoint, ix, iy, is, ny1;
    long effort, n_stable;

    log_entry("do_aperture_search_mp");

    log_entry("do_aperture_search_mp.1");
    coord     = (double**)zarray_2d(sizeof(**coord), ny+1, 7);
    coord[ny] = NULL;
    accepted  = (double**)zarray_2d(sizeof(**accepted), ny+1, 7);
    accepted[ny] = NULL;
    xy_left   = (double**)zarray_2d(sizeof(**xy_left), ny+1, 2); 
    xy_left[ny] = NULL;
    xy_right  = (double**)zarray_2d(sizeof(**xy_right), ny+1, 2); 
    xy_right[ny] = NULL;
    found     = (long*)tmalloc(sizeof(*found)*ny);
    n_left = n_right = 0;
    log_exit("do_aperture_search_mp.1");

    log_entry("do_aperture_search_mp.2");
    dx  = (xmax-xmin)/(nx-1);
    dy = (ymax-ymin)/(ny-1);
    effort = 0;
    n_stable = 0;

    for (iy=0, y=ymin; iy<ny; iy++, y+=dy) {
        xy_left[iy][1] = xy_right[iy][1] = y;
        xy_left[iy][0]  = xmin - dx;
        xy_right[iy][0] = xmax + dx;
        }

    ny1 = ny;
    fill_long_array(found, ny, 0L);
    log_exit("do_aperture_search_mp.2");

    log_entry("do_aperture_search_mp.3");
    while (ny1) {
        /* search from left */
        ny1 = 0;
        for (iy=0; iy<ny; iy++) {
            if (!found[iy]) {
                if ((xy_left[iy][0] += dx)>xmax) 
                    found[iy] = -1;
                else {
                    coord[ny1][0] = xy_left[iy][0];
                    coord[ny1][2] = xy_left[iy][1];
                    coord[ny1][1] = coord[ny1][3] = coord[ny1][4] = coord[ny1][5] = 0;
                    ny1++;
                    }
                }
            }
        if (!ny1)
            break;
        if (verbosity>1) {
            fprintf(stderr, "tracking %ld particles with x = %e:  \n", ny1, coord[0][0]);
            if (verbosity>2) {
                for (iy=0; iy<ny1; iy++) 
                    fprintf(stderr, "    y = %e ", coord[iy][2]);
                fputc('\n', stderr);
                }
            }
        p_central = run->p_central;
        n_trpoint = ny1;
        n_survived = do_tracking(coord, &n_trpoint, &effort, beamline, &p_central, 
                                 accepted, NULL, NULL, NULL, run, control->i_step, 
                                 SILENT_RUNNING, control->n_passes, NULL, NULL);
        if (verbosity>1) {
            fprintf(stderr, "    %ld particles survived\n", n_survived);
            }
        
        for (is=0; is<n_survived; is++) {
            if (verbosity>2)
                fprintf(stderr, "survivor: x = %e, y = %e\n", accepted[is][0], accepted[is][2]);
            for (iy=0; iy<ny; iy++) {
                if (!found[iy] && accepted[is][2]==xy_left[iy][1])
                    break;
                }
            if (iy!=ny) {
                iy = (accepted[is][2]-ymin)/dy + 0.5;
                if (iy<0 || iy>=ny)
                    bomb("invalid index (do_aperture_search.1)", NULL);
                found[iy] = 1;
                }
            else {
                fprintf(stderr, "error: data handling error (do_aperture_search)\n");
                fprintf(stderr, "no match for particle with y = %e\n", accepted[is][2]);
                abort();
                }
            }
        n_stable += n_survived;
        ny1 -= n_survived;
        }
    n_left = ny;
    for (iy=0; iy<n_left; iy++) {
        if (found[iy]!=1) {
            for (ix=iy+1; ix<ny; ix++) {
                found[ix-1] = found[ix];
                xy_left[ix-1][0] = xy_left[ix][0];
                xy_left[ix-1][1] = xy_left[ix][1];
                }
            iy--;
            n_left--;
            }
        }
    if (verbosity>1) {
        fprintf(stderr, "results for scan from left\n");
        for (iy=0; iy<n_left; iy++)
            fprintf(stderr, "    stable particle at x=%e, y=%e\n", xy_left[iy][0], xy_left[iy][1]);
        }
    log_exit("do_aperture_search_mp.3");

    log_entry("do_aperture_search_mp.4");
    ny1 = ny;
    fill_long_array(found, ny, 0L);
    while (ny1) {
        /* search from right */
        ny1 = 0;
        for (iy=0; iy<ny; iy++) {
            if (!found[iy]) {
                if ((xy_right[iy][0] -= dx)<xmin) 
                    found[iy] = -1;
                else {
                    coord[ny1][0] = xy_right[iy][0];
                    coord[ny1][2] = xy_right[iy][1];
                    coord[ny1][1] = coord[ny1][3] = coord[ny1][4] = coord[ny1][5] = 0;
                    ny1++;
                    }
                }
            }
        if (!ny1)
            break;
        if (verbosity>1) {
            fprintf(stderr, "tracking %ld particles with x = %e:  \n", ny1, coord[0][0]);
            if (verbosity>2) {
                for (iy=0; iy<ny1; iy++) 
                    fprintf(stderr, "    y = %e ", coord[iy][2]);
                fputc('\n', stderr);
                }
            }
        p_central = run->p_central;
        n_trpoint = ny1;
        n_survived = do_tracking(coord, &n_trpoint, &effort, beamline, &p_central, 
                                 accepted, NULL, NULL, NULL, run, control->i_step, 
                                 SILENT_RUNNING, control->n_passes, NULL, NULL);
        if (verbosity>1) {
            fprintf(stderr, "    %ld particles survived\n", n_survived);
            }
        
        for (is=0; is<n_survived; is++) {
            if (verbosity>2)
                fprintf(stderr, "survivor: x = %e, y = %e\n", accepted[is][0], accepted[is][2]);
            for (iy=0; iy<ny; iy++) {
                if (!found[iy] && accepted[is][2]==xy_right[iy][1])
                    break;
                }
            if (iy!=ny) {
                iy = (accepted[is][2]-ymin)/dy + 0.5;
                if (iy<0 || iy>=ny)
                    bomb("invalid index (do_aperture_search.1)", NULL);
                found[iy] = 1;
                }
            else {
                fprintf(stderr, "error: data handling error (do_aperture_search)\n");
                fprintf(stderr, "no match for particle with y = %e\n", accepted[is][2]);
                abort();
                }
            }
        n_stable += n_survived;
        ny1 -= n_survived;
        }
    n_right = ny;
    for (iy=0; iy<n_right; iy++) {
        if (found[iy]!=1) {
            for (ix=iy+1; ix<ny; ix++) {
                found[ix-1] = found[ix];
                if (!xy_right[ix-1]) {
                    fprintf(stderr, "error: xy_right[%ld] is NULL\n", ix-1);
                    abort();
                    }
                if (!xy_right[ix]) {
                    fprintf(stderr, "error: xy_right[%ld] is NULL\n", ix);
                    abort();
                    }
                xy_right[ix-1][0] = xy_right[ix][0];
                xy_right[ix-1][1] = xy_right[ix][1];
                }
            iy--;
            n_right--;
            }
        }
    if (verbosity>1) {
        fprintf(stderr, "results for scan from right\n");
        for (iy=0; iy<n_right; iy++)
            fprintf(stderr, "    stable particle at x=%e, y=%e\n", xy_right[iy][0], xy_right[iy][1]);
        }
    log_exit("do_aperture_search_mp.4");

    if (verbosity>0) {
        fprintf(stderr, "total effort:  %ld particle-turns   %ld stable particles were tracked\n", effort, n_stable);
        }

    log_entry("do_aperture_search_mp.5");

    if (!SDDS_StartTable(&SDDS_aperture, n_left+n_right)) {
        SDDS_SetError("Unable to start SDDS table (do_aperture_search)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
    SDDS_SetParameters(&SDDS_aperture, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, 0, control->i_step, -1);
    if (control->n_elements_to_vary) {
        for (ix=0; ix<control->n_elements_to_vary; ix++)
            if (!SDDS_SetParameters(&SDDS_aperture, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, ix+1,
                                    control->varied_quan_value[ix], -1))
                break;
        }
    if (SDDS_NumberOfErrors()) {
        SDDS_SetError("Problem setting SDDS parameter values (do_aperture_search)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }

    log_entry("do_aperture_search_mp.6");
    for (iy=0; iy<n_left; iy++) 
        if (!SDDS_SetRowValues(&SDDS_aperture, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, iy,
                               IC_X, xy_left[iy][0], IC_Y, xy_left[iy][1], -1)) {
            SDDS_SetError("Problem setting SDDS row values (do_aperture_search)");
            SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
            }
    for (iy=0; iy<n_right; iy++) {
        if (!SDDS_SetRowValues(&SDDS_aperture, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, iy+n_left,
                               IC_X, xy_right[n_right-iy-1][0], IC_Y, xy_right[n_right-iy-1][1], -1)) {
            SDDS_SetError("Problem setting SDDS row values (do_aperture_search)");
            SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
            }
        }
    log_exit("do_aperture_search_mp.6");
    if (!SDDS_WriteTable(&SDDS_aperture)) {
        SDDS_SetError("Problem writing SDDS table (do_aperture_search)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
        
    log_exit("do_aperture_search_mp.5");

    log_entry("do_aperture_search_mp.8");
    free_zarray_2d((void**)coord, ny, 7);
    free_zarray_2d((void**)accepted, ny, 7);
    free_zarray_2d((void**)xy_left, ny, 2);
    free_zarray_2d((void**)xy_right, ny, 2);
    free(found);
    log_exit("do_aperture_search_mp.8");

    log_exit("do_aperture_search_mp");
    return(1);
    }


long do_aperture_search_sp(
    RUN *run,
    VARY *control,
    ERROR *errcon,
    LINE_LIST *beamline
    )
{    
    double **coord;
    double x, y, dx, dy;
    double **xy_left, **xy_right;
    long n_left, n_right;
    double last_x_left, last_x_right, x1, x2;
    double p_central;
    long n_trpoint, ix, iy, is;
    long effort, n_stable;

    log_entry("do_aperture_search_sp");

    coord = (double**)zarray_2d(sizeof(**coord), 1, 7);
    xy_left   = (double**)zarray_2d(sizeof(**xy_left), ny, 2); 
    xy_right  = (double**)zarray_2d(sizeof(**xy_right), ny, 2); 
    n_left = n_right = 0;

    dx  = (xmax-xmin)/(nx-1);
    dy = (ymax-ymin)/(ny-1);
    last_x_left  = xmin;
    last_x_right = xmax;
    effort = 0;
    n_stable = 0;
    for (iy=0, y=ymin; iy<ny; iy++, y+=dy) {
        if (verbosity>0)
            fprintf(stderr, "searching for aperture for y = %e m\n", y);
        /* search from left */
        if (verbosity>1) {
            fprintf(stderr, "    searching from left to right\n");
            }
        if (assume_nonincreasing && iy!=0) {
            x = last_x_left;
            ix = (x - xmin)/dx + 0.5;
            if (ix>nx-1)
                ix = nx-1;
            }
        else {
            x = xmin;
            ix = 0;
            }
        for ( ; ix<nx; ix++, x+=dx) {
            if (verbosity>1) {
                fprintf(stderr, "    tracking for x = %e m\n", x);
                }
            coord[0][0] = x;
            coord[0][2] = y;
            p_central = run->p_central;
            coord[0][1] = coord[0][3] = coord[0][4] = coord[0][5] = 0;
            n_trpoint = 1;
            if (do_tracking(coord, &n_trpoint, &effort, beamline, &p_central, 
                            NULL, NULL, NULL, NULL, run, control->i_step, 
                            SILENT_RUNNING, control->n_passes, NULL, NULL))
                break;
            }
        if (ix!=nx) {
            n_stable++;
            last_x_left = x;
            if (ix!=0 && n_splits) {
                /* do secondary search */
                x1 = x;            /* stable   */
                x2 = x - dx;       /* unstable */
                for (is=0; is<n_splits; is++) {
                    if (fabs(x1-x2)<desired_resolution)
                        break;
                    x = (1-split_fraction)*x1 + split_fraction*x2;
                    if (verbosity>1) {
                        fprintf(stderr, "    splitting:  %e, %e --> %e \n", x1, x2, x);
                        }
                    coord[0][0] = x;
                    coord[0][2] = y;
                    p_central = run->p_central;
                    coord[0][1] = coord[0][3] = coord[0][4] = coord[0][5] = 0;
                    n_trpoint = 1;
                    if (do_tracking(coord, &n_trpoint, &effort, beamline, &p_central, 
                                    NULL, NULL, NULL, NULL, run, control->i_step, 
                                    SILENT_RUNNING, control->n_passes, NULL, NULL)) {
                        n_stable++;
                        x1 = x;    /* stable */
                        }
                    else
                        x2 = x;    /* unstable */
                    }
                x = x1;
                }
            xy_left[n_left][0] = x;
            xy_left[n_left][1] = y;
            if (verbosity>0) {
                fprintf(stderr, "    x = %e m is stable\n", x);
                }
            n_left++;
            }
        else {
            if (verbosity>0) {
                fprintf(stderr, "    no stable particles seen\n");
                }
            continue;
            }
        /* search from right */
        if (verbosity>1) {
            fprintf(stderr, "    searching from right to left\n");
            }
        if (assume_nonincreasing && iy!=0) {
            x = last_x_right;
            ix = (xmax - x)/dx + 0.5;
            if (ix>nx-1)
                ix = nx-1;
            }
        else {
            x = xmax;
            ix = 0;
            }
        for ( ; ix<nx; ix++, x-=dx) {
            if (verbosity>1) {
                fprintf(stderr, "    tracking for x = %e m\n", x);
                }
            coord[0][0] = x;
            coord[0][2] = y;
            p_central = run->p_central;
            coord[0][1] = coord[0][3] = coord[0][4] = coord[0][5] = 0;
            n_trpoint = 1;
            if (do_tracking(coord, &n_trpoint, &effort, beamline, &p_central, 
                            NULL, NULL, NULL, NULL, run, control->i_step, 
                            SILENT_RUNNING, control->n_passes, NULL, NULL))
                break;
            }
        if (ix!=nx) {
            n_stable++;
            last_x_right = x;
            if (ix!=0 && n_splits) {
                /* do secondary search */
                x1 = x;            /* stable   */
                x2 = x + dx;       /* unstable */
                for (is=0; is<n_splits; is++) {
                    if (fabs(x1-x2)<desired_resolution)
                        break;
                    x = (1-split_fraction)*x1 + split_fraction*x2;
                    if (verbosity>1) {
                        fprintf(stderr, "    splitting:  %e, %e --> %e \n", x1, x2, x);
                        }
                    coord[0][0] = x;
                    coord[0][2] = y;
                    p_central = run->p_central;
                    coord[0][1] = coord[0][3] = coord[0][4] = coord[0][5] = 0;
                    n_trpoint = 1;
                    if (do_tracking(coord, &n_trpoint, &effort, beamline, &p_central, 
                                    NULL, NULL, NULL, NULL, run, control->i_step, 
                                    SILENT_RUNNING, control->n_passes, NULL, NULL)) {
                        n_stable++;
                        x1 = x;    /* stable */
                        }
                    else
                        x2 = x;    /* unstable */
                    }
                x = x1;
                }
            xy_right[n_right][0] = x;
            xy_right[n_right][1] = y;
            if (verbosity>0) {
                fprintf(stderr, "    x = %e m is stable\n", x);
                }
            n_right++;
            }
        }
    if (verbosity>0) {
        fprintf(stderr, "total effort:  %ld particle-turns   %ld stable particles were tracked\n", effort, n_stable);
        }

    if (!SDDS_StartTable(&SDDS_aperture, n_left+n_right)) {
        SDDS_SetError("Unable to start SDDS table (do_aperture_search)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
    SDDS_SetParameters(&SDDS_aperture, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, 0, control->i_step, -1);
    if (control->n_elements_to_vary) {
        for (ix=0; ix<control->n_elements_to_vary; ix++)
            if (!SDDS_SetParameters(&SDDS_aperture, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, ix+1,
                                    control->varied_quan_value[ix], -1))
                break;
        }
    if (SDDS_NumberOfErrors()) {
        SDDS_SetError("Problem setting SDDS parameter values (do_aperture_search)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }

    for (iy=0; iy<n_left; iy++) 
        if (!SDDS_SetRowValues(&SDDS_aperture, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, iy,
                               IC_X, xy_left[iy][0], IC_Y, xy_left[iy][1], -1)) {
            SDDS_SetError("Problem setting SDDS row values (do_aperture_search)");
            SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
            }
    for (iy=0; iy<n_right; iy++) {
        if (!SDDS_SetRowValues(&SDDS_aperture, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, iy+n_left,
                               IC_X, xy_right[n_right-iy-1][0], IC_Y, xy_right[n_right-iy-1][1], -1)) {
            SDDS_SetError("Problem setting SDDS row values (do_aperture_search)");
            SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
            }
        }
    if (!SDDS_WriteTable(&SDDS_aperture)) {
        SDDS_SetError("Problem writing SDDS table (do_aperture_search)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
        

    free_zarray_2d((void**)coord, 1, 7);
    free_zarray_2d((void**)xy_left, ny, 2);
    free_zarray_2d((void**)xy_right, ny, 2);

    log_exit("do_aperture_search_sp");
    return(1);
    }

void finish_aperture_search(
    RUN *run,
    VARY *control,
    ERROR *errcon,
    LINE_LIST *beamline
    )
{
    if (SDDS_IsActive(&SDDS_aperture) && !SDDS_Terminate(&SDDS_aperture)) {
        SDDS_SetError("Problem terminating SDDS output (finish_aperture_search)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
    }

