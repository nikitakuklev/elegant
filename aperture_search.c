/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

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

static SDDS_DATASET SDDS_aperture;
static FILE *fpSearchOutput = NULL;

#define MP_MODE 0
#define SP_MODE 1
#define LINE_MODE 2
#define N_SEARCH_MODES 3
static char *search_mode[N_SEARCH_MODES] = {
  "many-particle", "single-particle", "particle-line"
    } ;
static long mode_code = 0;

long do_aperture_search_line(RUN *run, VARY *control, ERRORVAL *errcon, LINE_LIST *beamline);

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
    print_namelist(stdout, &find_aperture);

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

    fpSearchOutput = NULL;
    if (search_output) {
      if (mode_code!=SP_MODE) {
        fprintf(stdout, "Error: search_output field can only be used with single-particle mode\n");
        exit(1);
      }
      search_output = compose_filename(search_output, run->rootname);
      fpSearchOutput = fopen_e(search_output, "w", 0);
      fputs("SDDS1\n&parameter name=Step, type=long &end\n", fpSearchOutput);
      fputs("&parameter name=x0, type=double, units=m &end\n", fpSearchOutput);
      fputs("&parameter name=y0, type=double, units=m &end\n", fpSearchOutput);
      fputs("&parameter name=SearchFromRight, type=short &end\n", fpSearchOutput);
      fputs("&parameter name=IsStable, type=short &end\n", fpSearchOutput);
      fputs("&data mode=ascii no_row_counts=1 &end\n", fpSearchOutput);
    }
    
    log_exit("setup_aperture_search");
    }


long do_aperture_search(
    RUN *run,
    VARY *control,
    ERRORVAL *errcon,
    LINE_LIST *beamline
    )
{    
    long retcode;

    log_entry("do_aperture_search");
    switch (mode_code) {
    case MP_MODE:
      retcode = do_aperture_search_mp(run, control, errcon, beamline);
      break;
    case LINE_MODE:
      retcode = do_aperture_search_line(run, control, errcon, beamline);
      break;
    case SP_MODE:
    default:
      retcode = do_aperture_search_sp(run, control, errcon, beamline);
      break;
    }
    return(retcode);
}

/* many-particle search routine */

long do_aperture_search_mp(
    RUN *run,
    VARY *control,
    ERRORVAL *errcon,
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
            fprintf(stdout, "tracking %ld particles with x = %e:  \n", ny1, coord[0][0]);
            fflush(stdout);
            if (verbosity>2) {
                for (iy=0; iy<ny1; iy++) 
                    fprintf(stdout, "    y = %e ", coord[iy][2]);
                    fflush(stdout);
                fputc('\n', stdout);
                }
            }
        p_central = run->p_central;
        n_trpoint = ny1;
        n_survived = do_tracking(NULL, coord, n_trpoint, &effort, beamline, &p_central, 
                                 accepted, NULL, NULL, NULL, run, control->i_step, 
                                 SILENT_RUNNING, control->n_passes, 0, NULL, NULL, NULL, NULL);
        if (verbosity>1) {
            fprintf(stdout, "    %ld particles survived\n", n_survived);
            fflush(stdout);
            }
        
        for (is=0; is<n_survived; is++) {
            if (verbosity>2)
                fprintf(stdout, "survivor: x = %e, y = %e\n", accepted[is][0], accepted[is][2]);
                fflush(stdout);
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
                fprintf(stdout, "error: data handling error (do_aperture_search)\n");
                fflush(stdout);
                fprintf(stdout, "no match for particle with y = %e\n", accepted[is][2]);
                fflush(stdout);
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
        fprintf(stdout, "results for scan from left\n");
        fflush(stdout);
        for (iy=0; iy<n_left; iy++)
            fprintf(stdout, "    stable particle at x=%e, y=%e\n", xy_left[iy][0], xy_left[iy][1]);
            fflush(stdout);
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
            fprintf(stdout, "tracking %ld particles with x = %e:  \n", ny1, coord[0][0]);
            fflush(stdout);
            if (verbosity>2) {
                for (iy=0; iy<ny1; iy++) 
                    fprintf(stdout, "    y = %e ", coord[iy][2]);
                    fflush(stdout);
                fputc('\n', stdout);
                }
            }
        p_central = run->p_central;
        n_trpoint = ny1;
        n_survived = do_tracking(NULL, coord, n_trpoint, &effort, beamline, &p_central, 
                                 accepted, NULL, NULL, NULL, run, control->i_step, 
                                 SILENT_RUNNING, control->n_passes, 0, NULL, NULL, NULL, NULL);
        if (verbosity>1) {
            fprintf(stdout, "    %ld particles survived\n", n_survived);
            fflush(stdout);
            }
        
        for (is=0; is<n_survived; is++) {
            if (verbosity>2)
                fprintf(stdout, "survivor: x = %e, y = %e\n", accepted[is][0], accepted[is][2]);
                fflush(stdout);
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
                fprintf(stdout, "error: data handling error (do_aperture_search)\n");
                fflush(stdout);
                fprintf(stdout, "no match for particle with y = %e\n", accepted[is][2]);
                fflush(stdout);
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
                    fprintf(stdout, "error: xy_right[%ld] is NULL\n", ix-1);
                    fflush(stdout);
                    abort();
                    }
                if (!xy_right[ix]) {
                    fprintf(stdout, "error: xy_right[%ld] is NULL\n", ix);
                    fflush(stdout);
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
        fprintf(stdout, "results for scan from right\n");
        fflush(stdout);
        for (iy=0; iy<n_right; iy++)
            fprintf(stdout, "    stable particle at x=%e, y=%e\n", xy_right[iy][0], xy_right[iy][1]);
            fflush(stdout);
        }
    log_exit("do_aperture_search_mp.4");

    if (verbosity>0) {
        fprintf(stdout, "total effort:  %ld particle-turns   %ld stable particles were tracked\n", effort, n_stable);
        fflush(stdout);
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
    SDDS_DoFSync(&SDDS_aperture);
        
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
                           ERRORVAL *errcon,
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
      fprintf(stdout, "searching for aperture for y = %e m\n", y);
    fflush(stdout);
    /* search from left */
    if (verbosity>1) {
      fprintf(stdout, "    searching from left to right\n");
      fflush(stdout);
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
        fprintf(stdout, "    tracking for x = %e m\n", x);
        fflush(stdout);
      }
      coord[0][0] = x;
      coord[0][2] = y;
      p_central = run->p_central;
      coord[0][1] = coord[0][3] = coord[0][4] = coord[0][5] = 0;
      n_trpoint = 1;
      if (do_tracking(NULL, coord, n_trpoint, &effort, beamline, &p_central, 
                      NULL, NULL, NULL, NULL, run, control->i_step, 
                      SILENT_RUNNING, control->n_passes, 0, NULL, NULL, NULL, NULL)) {
        /* stable */
        if (fpSearchOutput)
          fprintf(fpSearchOutput, "%ld\n%le\n%le\n0\n1\n", control->i_step, x, y);
        break;
      } else {
        /* unstable */
        if (fpSearchOutput)
          fprintf(fpSearchOutput, "%ld\n%le\n%le\n0\n0\n", control->i_step, x, y);
      }
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
            fprintf(stdout, "    splitting:  %e, %e --> %e \n", x1, x2, x);
            fflush(stdout);
          }
          coord[0][0] = x;
          coord[0][2] = y;
          p_central = run->p_central;
          coord[0][1] = coord[0][3] = coord[0][4] = coord[0][5] = 0;
          n_trpoint = 1;
          if (do_tracking(NULL, coord, n_trpoint, &effort, beamline, &p_central, 
                          NULL, NULL, NULL, NULL, run, control->i_step, 
                          SILENT_RUNNING, control->n_passes, 0, NULL, NULL, NULL, NULL)) {
            n_stable++;
            x1 = x;    /* stable */
            if (fpSearchOutput)
              fprintf(fpSearchOutput, "%ld\n%le\n%le\n0\n1\n", control->i_step, x, y);
          }
          else {
            x2 = x;    /* unstable */
            if (fpSearchOutput)
              fprintf(fpSearchOutput, "%ld\n%le\n%le\n0\n0\n", control->i_step, x, y);
          }
        }
        x = x1;
      }
      xy_left[n_left][0] = x;
      xy_left[n_left][1] = y;
      if (verbosity>0) {
        fprintf(stdout, "    x = %e m is stable\n", x);
        fflush(stdout);
      }
      n_left++;
    }
    else {
      if (verbosity>0) {
        fprintf(stdout, "    no stable particles seen\n");
        fflush(stdout);
      }
      continue;
    }
    if (fpSearchOutput)
      fflush(fpSearchOutput);
    /* search from right */
    if (verbosity>1) {
      fprintf(stdout, "    searching from right to left\n");
      fflush(stdout);
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
        fprintf(stdout, "    tracking for x = %e m\n", x);
        fflush(stdout);
      }
      coord[0][0] = x;
      coord[0][2] = y;
      p_central = run->p_central;
      coord[0][1] = coord[0][3] = coord[0][4] = coord[0][5] = 0;
      n_trpoint = 1;
      if (do_tracking(NULL, coord, n_trpoint, &effort, beamline, &p_central, 
                      NULL, NULL, NULL, NULL, run, control->i_step, 
                      SILENT_RUNNING, control->n_passes, 0, NULL, NULL, NULL, NULL)) {
        /* stable */
        if (fpSearchOutput)
          fprintf(fpSearchOutput, "%ld\n%le\n%le\n1\n1\n", control->i_step, x, y);
        break;
      } else {
        /* unstable */
        if (fpSearchOutput)
          fprintf(fpSearchOutput, "%ld\n%le\n%le\n1\n0\n", control->i_step, x, y);
      }
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
            fprintf(stdout, "    splitting:  %e, %e --> %e \n", x1, x2, x);
            fflush(stdout);
          }
          coord[0][0] = x;
          coord[0][2] = y;
          p_central = run->p_central;
          coord[0][1] = coord[0][3] = coord[0][4] = coord[0][5] = 0;
          n_trpoint = 1;
          if (do_tracking(NULL, coord, n_trpoint, &effort, beamline, &p_central, 
                          NULL, NULL, NULL, NULL, run, control->i_step, 
                          SILENT_RUNNING, control->n_passes, 0, NULL, NULL, NULL, NULL)) {
            n_stable++;
            x1 = x;    /* stable */
            if (fpSearchOutput)
              fprintf(fpSearchOutput, "%ld\n%le\n%le\n1\n1\n", control->i_step, x, y);
          }
          else {
            x2 = x;    /* unstable */
            if (fpSearchOutput)
              fprintf(fpSearchOutput, "%ld\n%le\n%le\n1\n0\n", control->i_step, x, y);
          }
        }
        x = x1;
      }
      xy_right[n_right][0] = x;
      xy_right[n_right][1] = y;
      if (verbosity>0) {
        fprintf(stdout, "    x = %e m is stable\n", x);
        fflush(stdout);
      }
      n_right++;
    }
  }
  if (fpSearchOutput)
    fflush(fpSearchOutput);
  if (verbosity>0) {
    fprintf(stdout, "total effort:  %ld particle-turns   %ld stable particles were tracked\n", effort, n_stable);
    fflush(stdout);
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
  SDDS_DoFSync(&SDDS_aperture);
  
  free_zarray_2d((void**)coord, 1, 7);
  free_zarray_2d((void**)xy_left, ny, 2);
  free_zarray_2d((void**)xy_right, ny, 2);

  log_exit("do_aperture_search_sp");
  return(1);
}

void finish_aperture_search(
                            RUN *run,
                            VARY *control,
                            ERRORVAL *errcon,
                            LINE_LIST *beamline
                            )
{
  if (SDDS_IsActive(&SDDS_aperture) && !SDDS_Terminate(&SDDS_aperture)) {
    SDDS_SetError("Problem terminating SDDS output (finish_aperture_search)");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  if (fpSearchOutput) {
    fclose(fpSearchOutput);
    fpSearchOutput = NULL;
  }
}

/* line search routine */

long do_aperture_search_line(
    RUN *run,
    VARY *control,
    ERRORVAL *errcon,
    LINE_LIST *beamline
    )
{
    double **coord;
    double x0, y0, dx, dy;
    double p_central;
    long index, iSplit, nSteps;
    long effort, n_trpoint;
    double xSurvived, ySurvived ;

    coord     = (double**)zarray_2d(sizeof(**coord), 1, 7);

    effort = 0;
    xSurvived = ySurvived = -1;
    dx = dy = 0;

    for (iSplit=0; iSplit<=n_splits; iSplit++) {
      if (iSplit==0) {
	x0 = y0 = 0;
	dx  = xmax/(nx-1);
	dy = ymax/(nx-1);
	nSteps = nx;
      } else {
	if ((x0 = xSurvived)<=0 || (y0 = ySurvived)<=0) 
	  x0 = y0 = 0;
	dx *= split_fraction;
	dy *= split_fraction;
	x0 += dx;
	y0 += dy;
	nSteps = nx/split_fraction;
	if (verbosity>=1) {
	  printf("divided search interval to %e, %e\n", dx, dy);
	  fflush(stdout);
	}
      }

      for (index=0; index<nx; index++) {
	coord[0][1] = coord[0][3] = coord[0][4] = coord[0][5] = 0;
	coord[0][0] = index*dx + x0;
	coord[0][2] = index*dy + y0;
	
	p_central = run->p_central;
	n_trpoint = 1;
	if (do_tracking(NULL, coord, n_trpoint, &effort, beamline, &p_central, 
			  NULL, NULL, NULL, NULL, run, control->i_step, 
			  SILENT_RUNNING, control->n_passes, 0, NULL, NULL, NULL, NULL)!=1) {
	  if (verbosity>=2) {
	    fprintf(stdout, "particle lost for x=%e, y=%e\n", index*dx + x0, index*dy + y0);
	    fflush(stdout);
	  }
	  break;
	}
	if (verbosity>=2) {
	  fprintf(stdout, "particle survived for x=%e, y=%e\n", x0+index*dx, y0+index*dy);
	  fflush(stdout);
	}
	if (xSurvived<(x0+index*dx)) {
	  xSurvived = x0+index*dx;
	  ySurvived = y0+index*dy;      
	}
      }
      
      if (verbosity>=1) {
        fprintf(stdout, "Sweep done, particle survived up to x=%e, y=%e\n", xSurvived, ySurvived);
        fflush(stdout);
      }
    }
    
    if (!SDDS_StartTable(&SDDS_aperture, 1)) {
      SDDS_SetError("Unable to start SDDS table (do_aperture_search)");
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
    SDDS_SetParameters(&SDDS_aperture, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, 0, control->i_step, -1);
    if (control->n_elements_to_vary) {
      for (index=0; index<control->n_elements_to_vary; index++)
	if (!SDDS_SetParameters(&SDDS_aperture, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, index+1,
				control->varied_quan_value[index], -1))
	  break;
    }
    if (SDDS_NumberOfErrors()) {
      SDDS_SetError("Problem setting SDDS parameter values (do_aperture_search)");
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
    
    if (!SDDS_SetRowValues(&SDDS_aperture, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, 0,
			   IC_X, xSurvived, IC_Y, ySurvived, -1)) {
      SDDS_SetError("Problem setting SDDS row values (do_aperture_search)");
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
    
    if (!SDDS_WriteTable(&SDDS_aperture)) {
      SDDS_SetError("Problem writing SDDS table (do_aperture_search)");
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
    SDDS_DoFSync(&SDDS_aperture);
    
    free_zarray_2d((void**)coord, 1, 7);
    return(1);
}

