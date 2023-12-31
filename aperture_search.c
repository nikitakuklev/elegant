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
#define IC_XC 2
#define IC_YC 3
#define IC_SLOST 4
#define IC_XLOST 5
#define IC_YLOST 6
#define N_COLUMNS 7
static SDDS_DEFINITION column_definition[N_COLUMNS] = {
  {"x", "&column name=x, symbol=x, units=m, type=double &end"},
  {"y", "&column name=y, symbol=y, units=m, type=double &end"},
  {"xClipped", "&column name=xClipped, symbol=xClipped, units=m, type=double &end"},
  {"yClipped", "&column name=yClipped, symbol=yClipped, units=m, type=double &end"},
  {"sLost", "&column name=sLost, units=m, type=double &end"},
  {"xLost", "&column name=xLost, units=m, type=double &end"},
  {"yLost", "&column name=yLost, units=m type=double &end"},
};

static SDDS_DEFINITION columnp_definition[N_COLUMNS] = {
  {"xp", "&column name=xp, type=double &end"},
  {"yp", "&column name=yp, type=double &end"},
  {"xpClipped", "&column name=xpClipped, type=double &end"},
  {"ypClipped", "&column name=ypClipped, type=double &end"},
  {"sLost", "&column name=sLost, units=m, type=double &end"},
  {"xLost", "&column name=xLost, units=m, type=double &end"},
  {"yLost", "&column name=yLost, units=m type=double &end"},
};

#define IP_STEP 0
#define IP_SVN 1
#define IP_AREA 2
#define N_PARAMETERS 3
static SDDS_DEFINITION parameter_definition[N_PARAMETERS] = {
  {"Step", "&parameter name=Step, type=long, description=\"Simulation step\" &end"},
  {"SVNVersion", "&parameter name=SVNVersion, type=string, description=\"SVN version number\", fixed_value=" SVN_VERSION " &end"},
  {"Area", "&parameter name=Area, type=double, units=\"m$a2$n\" &end"},
};
static SDDS_DEFINITION parameterp_definition[N_PARAMETERS] = {
  {"Step", "&parameter name=Step, type=long, description=\"Simulation step\" &end"},
  {"SVNVersion", "&parameter name=SVNVersion, type=string, description=\"SVN version number\", fixed_value=" SVN_VERSION " &end"},
  {"Area", "&parameter name=Area, type=double, &end"},
};

static SDDS_DATASET SDDS_aperture;
static FILE *fpSearchOutput = NULL;

#define MP_MODE 0
#define SP_MODE 1
#define ONE_LINE_MODE 2
/* This one needs to be the first line mode */
#define TWO_LINE_MODE 3
#define THREE_LINE_MODE 4
#define FIVE_LINE_MODE 5
#define SEVEN_LINE_MODE 6
#define NINE_LINE_MODE 7
#define ELEVEN_LINE_MODE 8
#define N_LINE_MODE 9
/* This one needs to be the last line mode */
#define LINE_MODE 10
#define N_SEARCH_MODES 11
static char *search_mode[N_SEARCH_MODES] = {
  "many-particle",
  "single-particle",
  "one-line",
  "two-line",
  "three-line",
  "five-line",
  "seven-line",
  "nine-line",
  "eleven-line",
  "n-line",
  "particle-line",
};
static long mode_code = 0;
static long slope_mode = 0;

#if USE_MPI
long do_aperture_search_line_p(RUN *run, VARY *control, double *referenceCoord, ERRORVAL *errcon,
                               LINE_LIST *beamline, long lines, double *returnValue);
#endif

void setup_aperture_search(
  NAMELIST_TEXT *nltext,
  RUN *run,
  VARY *control,
  long *optimizationMode) {
  char description[200];

  log_entry("setup_aperture_search");

  /* process namelist input */
  set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
  set_print_namelist_flags(0);
  if (processNamelist(&find_aperture, nltext) == NAMELIST_ERROR)
    bombElegant(NULL, NULL);
  if (echoNamelists)
    print_namelist(stdout, &find_aperture);

  /* check for data errors */
  if (!output && !optimization_mode)
    bombElegant("no output filename specified (required if optimization_mode=0)", NULL);
  if (xmin == xmax) {
    slope_mode = 1;
    if (xpmin == xpmax || ypmin == ypmax)
      bombElegant("if xmin==xmax, must not have xpmin==xpmax or ypmin==ypmax", NULL);
    if (ymin != ymax)
      bombElegant("if xmin==xmax, must also have ymin==ymax", NULL);
  } else {
    if (xmin >= xmax)
      bombElegant("xmin >= xmax", NULL);
    if (ymin >= ymax)
      bombElegant("ymin >= ymax", NULL);
    if (xpmin != xpmax || ypmin != ypmax)
      bombElegant("if xmin!=xmax, must have xpmin==xpmax and ypmin==ypmax", NULL);
  }

  if (nx < 3)
    bombElegant("nx < 3", NULL);
  if (ny < 2)
    bombElegant("ny < 2", NULL);
  if (n_splits && n_splits < 1)
    bombElegant("n_splits is non-zero, but less than 1", NULL);
  if (n_splits) {
    if (split_fraction <= 0 || split_fraction >= 1)
      bombElegant("split_fraction must be greater than 0 and less than 1", NULL);
    if (desired_resolution <= 0 || desired_resolution >= 1)
      bombElegant("desired_resolution must be greater than 0 and less than 1", NULL);
    if ((desired_resolution *= (xmax - xmin)) > (xmax - xmin) / (nx - 1))
      bombElegant("desired_resolution is larger than coarse mesh", NULL);
  }
  if ((mode_code = match_string(mode, search_mode, N_SEARCH_MODES, 0)) < 0)
    bombElegant("unknown search mode", NULL);
  if (slope_mode && mode_code != N_LINE_MODE)
    bombElegant("slope-based aperture search only supported for mode=\"n-line\"", NULL);
#if USE_MPI
  if (mode_code != N_LINE_MODE) {
    bombElegant("only mode=\"n-line\" is supported by parallel elegant", NULL);
  }
#endif
  if (mode_code < ONE_LINE_MODE) {
    printWarning("find_aperture: using one of the line-based search methods is strongly recommended",
                 "other methods, which search from large amplitude to small amplitude, may overstate the aperture because of stable islands");
  }
  if (full_plane != 0 && mode_code < ONE_LINE_MODE)
    bombElegant("full_plane=1 is only supported for line-based modes at present", NULL);

  if (optimization_mode && (mode_code <= TWO_LINE_MODE || (mode_code == LINE_MODE && n_lines < 3)))
    bombElegant("dynamic aperture optimization requires use of n-line mode with at least 3 lines", NULL);
  if (offset_by_orbit && mode_code == SP_MODE)
    bombElegant("can't presently offset_by_orbit for that mode", NULL);

#if USE_MPI
  watch_not_allowed = 1;
  if (isMaster) {
#endif
    if (!optimization_mode) {
      output = compose_filename(output, run->rootname);
      sprintf(description, "%s aperture search", search_mode[mode_code]);
      SDDS_ElegantOutputSetup(&SDDS_aperture, output, SDDS_BINARY, 1,
                              description, run->runfile, run->lattice, slope_mode ? parameterp_definition : parameter_definition,
                              N_PARAMETERS - (mode_code >= ONE_LINE_MODE && mode_code <= LINE_MODE ? 0 : 1),
                              slope_mode ? columnp_definition : column_definition,
                              N_COLUMNS - (mode_code >= ONE_LINE_MODE && mode_code <= LINE_MODE ? 0 : 3),
                              "setup_aperture_search", SDDS_EOS_NEWFILE);
      if (control->n_elements_to_vary)
        if (!SDDS_DefineSimpleParameters(&SDDS_aperture, control->n_elements_to_vary,
                                         control->varied_quan_name, control->varied_quan_unit, SDDS_DOUBLE)) {
          SDDS_SetError("Unable to define additional SDDS parameters (setup_aperture_search)");
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors | SDDS_EXIT_PrintErrors);
        }

      if (!SDDS_WriteLayout(&SDDS_aperture)) {
        SDDS_SetError("Unable to write SDDS layout for aperture search");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors | SDDS_EXIT_PrintErrors);
      }

      if (boundary && (mode_code == SP_MODE || mode_code == MP_MODE)) {
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
        if (mode_code != SP_MODE) {
          printf("Error: search_output field can only be used with single-particle mode\n");
          exitElegant(1);
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
    }
#if USE_MPI
  } else /* set output to NULL for the slave processors */
    output = NULL;

  if (optimization_mode) {
    notSinglePart = 0; /* All the processors will track independently */
    lessPartAllowed = 1;
  }
#endif

  *optimizationMode = optimization_mode;

  log_exit("setup_aperture_search");
}

long do_aperture_search(
  RUN *run,
  VARY *control,
  double *referenceCoord,
  ERRORVAL *errcon,
  LINE_LIST *beamline,
  double *returnValue) {
  long retcode;
  *returnValue = 0;

  log_entry("do_aperture_search");
  switch (mode_code) {
  case N_LINE_MODE:
    if (n_lines <= 2)
      bombElegant("n_lines must be greater than 2 for aperture search", NULL);
    retcode = do_aperture_search_line(run, control, referenceCoord, errcon, beamline, n_lines, returnValue);
    break;
  case MP_MODE:
    retcode = do_aperture_search_mp(run, control, referenceCoord, errcon, beamline);
    break;
  case ONE_LINE_MODE:
  case LINE_MODE:
    retcode = do_aperture_search_line(run, control, referenceCoord, errcon, beamline, 1, returnValue);
    break;
  case TWO_LINE_MODE:
    retcode = do_aperture_search_line(run, control, referenceCoord, errcon, beamline, 2, returnValue);
    break;
  case THREE_LINE_MODE:
    retcode = do_aperture_search_line(run, control, referenceCoord, errcon, beamline, 3, returnValue);
    break;
  case FIVE_LINE_MODE:
    retcode = do_aperture_search_line(run, control, referenceCoord, errcon, beamline, 5, returnValue);
    break;
  case SEVEN_LINE_MODE:
    retcode = do_aperture_search_line(run, control, referenceCoord, errcon, beamline, 7, returnValue);
    break;
  case NINE_LINE_MODE:
    retcode = do_aperture_search_line(run, control, referenceCoord, errcon, beamline, 9, returnValue);
    break;
  case ELEVEN_LINE_MODE:
    retcode = do_aperture_search_line(run, control, referenceCoord, errcon, beamline, 11, returnValue);
    break;
  case SP_MODE:
  default:
    retcode = do_aperture_search_sp(run, control, referenceCoord, errcon, beamline);
    break;
  }
  return (retcode);
}

/* many-particle search routine */

long do_aperture_search_mp(
  RUN *run,
  VARY *control,
  double *referenceCoord,
  ERRORVAL *errcon,
  LINE_LIST *beamline) {
  double **coord, **accepted;
  double y, dx, dy;
  double **xy_left, **xy_right;
  long *found;
  long n_left, n_right, n_survived, ic;
  double p_central;
  long n_trpoint, ix, iy, is, ny1;
  long effort, n_stable;
  double orbit[6] = {0, 0, 0, 0, 0, 0};

  log_entry("do_aperture_search_mp");

  log_entry("do_aperture_search_mp.1");
  coord = (double **)czarray_2d(sizeof(**coord), ny + 1, totalPropertiesPerParticle);
  coord[ny] = NULL;
  accepted = (double **)czarray_2d(sizeof(**accepted), ny + 1, totalPropertiesPerParticle);
  accepted[ny] = NULL;
  xy_left = (double **)czarray_2d(sizeof(**xy_left), ny + 1, 2);
  xy_left[ny] = NULL;
  xy_right = (double **)czarray_2d(sizeof(**xy_right), ny + 1, 2);
  xy_right[ny] = NULL;
  found = (long *)tmalloc(sizeof(*found) * ny);
  n_left = n_right = 0;

  if (offset_by_orbit)
    memcpy(orbit, referenceCoord, sizeof(*referenceCoord) * 6);

  log_exit("do_aperture_search_mp.1");

  log_entry("do_aperture_search_mp.2");
  dx = (xmax - xmin) / (nx - 1);
  dy = (ymax - ymin) / (ny - 1);
  effort = 0;
  n_stable = 0;

  for (iy = 0, y = ymin; iy < ny; iy++, y += dy) {
    xy_left[iy][1] = xy_right[iy][1] = y;
    xy_left[iy][0] = xmin - dx;
    xy_right[iy][0] = xmax + dx;
  }

  ny1 = ny;
  fill_long_array(found, ny, 0L);
  log_exit("do_aperture_search_mp.2");

  log_entry("do_aperture_search_mp.3");
  while (ny1) {
    /* search from left */
    ny1 = 0;
    for (iy = 0; iy < ny; iy++) {
      if (!found[iy]) {
        if ((xy_left[iy][0] += dx) > xmax)
          found[iy] = -1;
        else {
          coord[ny1][0] = xy_left[iy][0];
          coord[ny1][2] = xy_left[iy][1];
          coord[ny1][1] = coord[ny1][3] = coord[ny1][4] = coord[ny1][5] = 0;
          coord[ny1][6] = iy;
          ny1++;
        }
      }
    }
    if (!ny1)
      break;
    if (verbosity > 1) {
      printf("tracking %ld particles with x = %e:  \n", ny1, coord[0][0]);
      fflush(stdout);
      if (verbosity > 2) {
        for (iy = 0; iy < ny1; iy++)
          printf("    y = %e ", coord[iy][2]);
        fflush(stdout);
        fputc('\n', stdout);
      }
    }
    p_central = run->p_central;
    n_trpoint = ny1;
    for (iy = 0; iy < ny1; iy++)
      for (ic = 0; ic < 6; ic++)
        coord[iy][ic] += orbit[ic];
    n_survived = do_tracking(NULL, coord, n_trpoint, &effort, beamline, &p_central,
                             accepted, NULL, NULL, NULL, run, control->i_step,
                             SILENT_RUNNING, control->n_passes, 0, NULL, NULL, NULL, NULL, NULL);
    if (verbosity > 1) {
      printf("    %ld particles survived\n", n_survived);
      fflush(stdout);
    }

    for (is = 0; is < n_survived; is++) {
      if (verbosity > 2)
        printf("survivor: x = %e, y = %e\n",
               accepted[is][0] - orbit[0], accepted[is][2] - orbit[2]);
      fflush(stdout);

      iy = accepted[is][6];
      if (iy < 0 || iy >= ny)
        bombElegant("invalid index (do_aperture_search.1)", NULL);
      found[iy] = 1;
    }
    n_stable += n_survived;
    ny1 -= n_survived;
  }
  n_left = ny;
  for (iy = 0; iy < n_left; iy++) {
    if (found[iy] != 1) {
      for (ix = iy + 1; ix < ny; ix++) {
        found[ix - 1] = found[ix];
        xy_left[ix - 1][0] = xy_left[ix][0];
        xy_left[ix - 1][1] = xy_left[ix][1];
      }
      iy--;
      n_left--;
    }
  }
  if (verbosity > 1) {
    printf("results for scan from left\n");
    fflush(stdout);
    for (iy = 0; iy < n_left; iy++)
      printf("    stable particle at x=%e, y=%e\n", xy_left[iy][0], xy_left[iy][1]);
    fflush(stdout);
  }
  log_exit("do_aperture_search_mp.3");

  log_entry("do_aperture_search_mp.4");
  ny1 = ny;
  fill_long_array(found, ny, 0L);
  while (ny1) {
    /* search from right */
    ny1 = 0;
    for (iy = 0; iy < ny; iy++) {
      if (!found[iy]) {
        if ((xy_right[iy][0] -= dx) < xmin)
          found[iy] = -1;
        else {
          coord[ny1][0] = xy_right[iy][0];
          coord[ny1][2] = xy_right[iy][1];
          coord[ny1][1] = coord[ny1][3] = coord[ny1][4] = coord[ny1][5] = 0;
          coord[ny1][6] = iy;
          ny1++;
        }
      }
    }
    if (!ny1)
      break;
    if (verbosity > 1) {
      printf("tracking %ld particles with x = %e:  \n", ny1, coord[0][0]);
      fflush(stdout);
      if (verbosity > 2) {
        for (iy = 0; iy < ny1; iy++)
          printf("    y = %e ", coord[iy][2]);
        fflush(stdout);
        fputc('\n', stdout);
      }
    }
    p_central = run->p_central;
    n_trpoint = ny1;
    for (iy = 0; iy < ny1; iy++)
      for (ic = 0; ic < 6; ic++)
        coord[iy][ic] += orbit[ic];
    n_survived = do_tracking(NULL, coord, n_trpoint, &effort, beamline, &p_central,
                             accepted, NULL, NULL, NULL, run, control->i_step,
                             SILENT_RUNNING, control->n_passes, 0, NULL, NULL, NULL, NULL, NULL);
    if (verbosity > 1) {
      printf("    %ld particles survived\n", n_survived);
      fflush(stdout);
    }

    for (is = 0; is < n_survived; is++) {
      if (verbosity > 2)
        printf("survivor: x = %e, y = %e\n",
               accepted[is][0] - orbit[0], accepted[is][2] - orbit[2]);
      fflush(stdout);

      iy = accepted[is][6];
      if (iy < 0 || iy >= ny)
        bombElegant("invalid index (do_aperture_search.1)", NULL);
      found[iy] = 1;
    }
    n_stable += n_survived;
    ny1 -= n_survived;
  }
  n_right = ny;
  for (iy = 0; iy < n_right; iy++) {
    if (found[iy] != 1) {
      for (ix = iy + 1; ix < ny; ix++) {
        found[ix - 1] = found[ix];
        if (!xy_right[ix - 1]) {
          printf("error: xy_right[%ld] is NULL\n", ix - 1);
          fflush(stdout);
          abort();
        }
        if (!xy_right[ix]) {
          printf("error: xy_right[%ld] is NULL\n", ix);
          fflush(stdout);
          abort();
        }
        xy_right[ix - 1][0] = xy_right[ix][0];
        xy_right[ix - 1][1] = xy_right[ix][1];
      }
      iy--;
      n_right--;
    }
  }
  if (verbosity > 1) {
    printf("results for scan from right\n");
    fflush(stdout);
    for (iy = 0; iy < n_right; iy++)
      printf("    stable particle at x=%e, y=%e\n", xy_right[iy][0], xy_right[iy][1]);
    fflush(stdout);
  }
  log_exit("do_aperture_search_mp.4");

  if (verbosity > 0) {
    printf("total effort:  %ld particle-turns   %ld stable particles were tracked\n", effort, n_stable);
    fflush(stdout);
  }

  log_entry("do_aperture_search_mp.5");

  if (!SDDS_StartTable(&SDDS_aperture, n_left + n_right)) {
    SDDS_SetError("Unable to start SDDS table (do_aperture_search)");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors | SDDS_EXIT_PrintErrors);
  }
  SDDS_SetParameters(&SDDS_aperture, SDDS_SET_BY_INDEX | SDDS_PASS_BY_VALUE, IP_STEP, control->i_step, -1);
  if (control->n_elements_to_vary) {
    for (ix = 0; ix < control->n_elements_to_vary; ix++)
      if (!SDDS_SetParameters(&SDDS_aperture, SDDS_SET_BY_INDEX | SDDS_PASS_BY_VALUE, ix + N_PARAMETERS,
                              control->varied_quan_value[ix], -1))
        break;
  }
  if (SDDS_NumberOfErrors()) {
    SDDS_SetError("Problem setting SDDS parameter values (do_aperture_search)");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors | SDDS_EXIT_PrintErrors);
  }

  log_entry("do_aperture_search_mp.6");
  for (iy = 0; iy < n_left; iy++)
    if (!SDDS_SetRowValues(&SDDS_aperture, SDDS_SET_BY_INDEX | SDDS_PASS_BY_VALUE, iy,
                           IC_X, xy_left[iy][0], IC_Y, xy_left[iy][1], -1)) {
      SDDS_SetError("Problem setting SDDS row values (do_aperture_search)");
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors | SDDS_EXIT_PrintErrors);
    }
  for (iy = 0; iy < n_right; iy++) {
    if (!SDDS_SetRowValues(&SDDS_aperture, SDDS_SET_BY_INDEX | SDDS_PASS_BY_VALUE, iy + n_left,
                           IC_X, xy_right[n_right - iy - 1][0], IC_Y, xy_right[n_right - iy - 1][1], -1)) {
      SDDS_SetError("Problem setting SDDS row values (do_aperture_search)");
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors | SDDS_EXIT_PrintErrors);
    }
  }
  log_exit("do_aperture_search_mp.6");
  if (!SDDS_WriteTable(&SDDS_aperture)) {
    SDDS_SetError("Problem writing SDDS table (do_aperture_search)");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors | SDDS_EXIT_PrintErrors);
  }
  if (!inhibitFileSync)
    SDDS_DoFSync(&SDDS_aperture);
  log_exit("do_aperture_search_mp.5");

  log_entry("do_aperture_search_mp.8");
  free_czarray_2d((void **)coord, ny, totalPropertiesPerParticle);
  free_czarray_2d((void **)accepted, ny, totalPropertiesPerParticle);
  free_czarray_2d((void **)xy_left, ny, 2);
  free_czarray_2d((void **)xy_right, ny, 2);
  free(found);
  log_exit("do_aperture_search_mp.8");

  log_exit("do_aperture_search_mp");
  return (1);
}

long do_aperture_search_sp(
  RUN *run,
  VARY *control,
  double *referenceCoord,
  ERRORVAL *errcon,
  LINE_LIST *beamline) {
  double **coord;
  double x, y, dx, dy;
  double **xy_left, **xy_right;
  long n_left, n_right;
  double last_x_left, last_x_right, x1, x2;
  double p_central;
  long n_trpoint, ix, iy, is;
  long effort, n_stable;

  log_entry("do_aperture_search_sp");

  coord = (double **)czarray_2d(sizeof(**coord), 1, totalPropertiesPerParticle);
  xy_left = (double **)czarray_2d(sizeof(**xy_left), ny, 2);
  xy_right = (double **)czarray_2d(sizeof(**xy_right), ny, 2);
  n_left = n_right = 0;

  dx = (xmax - xmin) / (nx - 1);
  dy = (ymax - ymin) / (ny - 1);
  last_x_left = xmin;
  last_x_right = xmax;
  effort = 0;
  n_stable = 0;
  for (iy = 0, y = ymin; iy < ny; iy++, y += dy) {
#if USE_MPI
    /* A y value will be searched with one CPU */
    if (myid != iy % n_processors)
      continue;
#endif
    if (verbosity > 0)
      printf("searching for aperture for y = %e m\n", y);
    fflush(stdout);
    if (verbosity > 1) {
      printf("    searching from left to right\n");
      fflush(stdout);
    }
    if (assume_nonincreasing && iy != 0) {
      x = last_x_left;
      ix = (x - xmin) / dx + 0.5;
      if (ix > nx - 1)
        ix = nx - 1;
    } else {
      x = xmin;
      ix = 0;
    }
    for (; ix < nx; ix++, x += dx) {
      if (verbosity > 1) {
        printf("    tracking for x = %e m\n", x);
        fflush(stdout);
      }
      coord[0][0] = x;
      coord[0][2] = y;
      p_central = run->p_central;
      coord[0][1] = coord[0][3] = coord[0][4] = coord[0][5] = 0;
      n_trpoint = 1;
      if (do_tracking(NULL, coord, n_trpoint, &effort, beamline, &p_central,
                      NULL, NULL, NULL, NULL, run, control->i_step,
                      SILENT_RUNNING, control->n_passes, 0, NULL, NULL, NULL, NULL, NULL)) {
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
    if (ix != nx) {
      n_stable++;
      last_x_left = x;
      if (ix != 0 && n_splits) {
        /* do secondary search */
        x1 = x;      /* stable   */
        x2 = x - dx; /* unstable */
        for (is = 0; is < n_splits; is++) {
          if (fabs(x1 - x2) < desired_resolution)
            break;
          x = (1 - split_fraction) * x1 + split_fraction * x2;
          if (verbosity > 1) {
            printf("    splitting:  %e, %e --> %e \n", x1, x2, x);
            fflush(stdout);
          }
          coord[0][0] = x;
          coord[0][2] = y;
          p_central = run->p_central;
          coord[0][1] = coord[0][3] = coord[0][4] = coord[0][5] = 0;
          n_trpoint = 1;
          if (do_tracking(NULL, coord, n_trpoint, &effort, beamline, &p_central,
                          NULL, NULL, NULL, NULL, run, control->i_step,
                          SILENT_RUNNING, control->n_passes, 0, NULL, NULL, NULL, NULL, NULL)) {
            n_stable++;
            x1 = x; /* stable */
            if (fpSearchOutput)
              fprintf(fpSearchOutput, "%ld\n%le\n%le\n0\n1\n", control->i_step, x, y);
          } else {
            x2 = x; /* unstable */
            if (fpSearchOutput)
              fprintf(fpSearchOutput, "%ld\n%le\n%le\n0\n0\n", control->i_step, x, y);
          }
        }
        x = x1;
      }
      xy_left[n_left][0] = x;
      xy_left[n_left][1] = y;
      if (verbosity > 0) {
        printf("    x = %e m is stable\n", x);
        fflush(stdout);
      }
      n_left++;
    } else {
      if (verbosity > 0) {
        printf("    no stable particles seen\n");
        fflush(stdout);
      }
      continue;
    }
    if (fpSearchOutput)
      fflush(fpSearchOutput);
    /* search from right */
    if (verbosity > 1) {
      printf("    searching from right to left\n");
      fflush(stdout);
    }
    if (assume_nonincreasing && iy != 0) {
      x = last_x_right;
      ix = (xmax - x) / dx + 0.5;
      if (ix > nx - 1)
        ix = nx - 1;
    } else {
      x = xmax;
      ix = 0;
    }
    for (; ix < nx; ix++, x -= dx) {
      if (verbosity > 1) {
        printf("    tracking for x = %e m\n", x);
        fflush(stdout);
      }
      coord[0][0] = x;
      coord[0][2] = y;
      p_central = run->p_central;
      coord[0][1] = coord[0][3] = coord[0][4] = coord[0][5] = 0;
      n_trpoint = 1;
      if (do_tracking(NULL, coord, n_trpoint, &effort, beamline, &p_central,
                      NULL, NULL, NULL, NULL, run, control->i_step,
                      SILENT_RUNNING, control->n_passes, 0, NULL, NULL, NULL, NULL, NULL)) {
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
    if (ix != nx) {
      n_stable++;
      last_x_right = x;
      if (ix != 0 && n_splits) {
        /* do secondary search */
        x1 = x;      /* stable   */
        x2 = x + dx; /* unstable */
        for (is = 0; is < n_splits; is++) {
          if (fabs(x1 - x2) < desired_resolution)
            break;
          x = (1 - split_fraction) * x1 + split_fraction * x2;
          if (verbosity > 1) {
            printf("    splitting:  %e, %e --> %e \n", x1, x2, x);
            fflush(stdout);
          }
          coord[0][0] = x;
          coord[0][2] = y;
          p_central = run->p_central;
          coord[0][1] = coord[0][3] = coord[0][4] = coord[0][5] = 0;
          n_trpoint = 1;
          if (do_tracking(NULL, coord, n_trpoint, &effort, beamline, &p_central,
                          NULL, NULL, NULL, NULL, run, control->i_step,
                          SILENT_RUNNING, control->n_passes, 0, NULL, NULL, NULL, NULL, NULL)) {
            n_stable++;
            x1 = x; /* stable */
            if (fpSearchOutput)
              fprintf(fpSearchOutput, "%ld\n%le\n%le\n1\n1\n", control->i_step, x, y);
          } else {
            x2 = x; /* unstable */
            if (fpSearchOutput)
              fprintf(fpSearchOutput, "%ld\n%le\n%le\n1\n0\n", control->i_step, x, y);
          }
        }
        x = x1;
      }
      xy_right[n_right][0] = x;
      xy_right[n_right][1] = y;
      if (verbosity > 0) {
        printf("    x = %e m is stable\n", x);
        fflush(stdout);
      }
      n_right++;
    }
  }
  if (fpSearchOutput)
    fflush(fpSearchOutput);
#if !USE_MPI
  if (verbosity > 0) {
    printf("total effort:  %ld particle-turns   %ld stable particles were tracked\n", effort, n_stable);
    fflush(stdout);
  }
#else
  /* Gather all the simulation result to master to write into a file */
  if (USE_MPI) {
    int *n_vector = (int *)tmalloc(n_processors * sizeof(*n_vector));
    int *offset = (int *)tmalloc(n_processors * sizeof(*offset));
    long i;

    MPI_Allgather(&n_left, 1, MPI_LONG, n_vector, 1, MPI_INT, MPI_COMM_WORLD);
    offset[0] = 0;
    for (i = 0; i < n_processors - 1; i++) {
      n_vector[i] *= 2;
      offset[i + 1] = offset[i] + n_vector[i];
    }
    n_vector[n_processors - 1] *= 2;
    MPI_Allgatherv(xy_left[0], 2 * n_left, MPI_DOUBLE, xy_left[0], n_vector, offset, MPI_DOUBLE, MPI_COMM_WORLD);

    MPI_Allgather(&n_right, 1, MPI_LONG, n_vector, 1, MPI_INT, MPI_COMM_WORLD);
    offset[0] = 0;
    for (i = 0; i < n_processors - 1; i++) {
      n_vector[i] *= 2;
      offset[i + 1] = offset[i] + n_vector[i];
    }
    n_vector[n_processors - 1] *= 2;
    MPI_Allgatherv(xy_right[0], 2 * n_right, MPI_DOUBLE, xy_right[0], n_vector, offset, MPI_DOUBLE, MPI_COMM_WORLD);
  }

  if (USE_MPI) {
    long effort_total, n_stable_total;
    long tmp[4], total[4];

    tmp[0] = effort;
    tmp[1] = n_stable;
    tmp[2] = n_left;
    tmp[3] = n_right;
    MPI_Reduce(tmp, &total, 4, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    if (isMaster) {
      effort_total = total[0];
      n_stable_total = total[1];
      n_left = total[2];
      n_right = total[3];
      if (verbosity > 0) {
        printf("total effort:  %ld particle-turns   %ld stable particles were tracked\n", effort_total, n_stable_total);
        fflush(stdout);
      }
    }
  }
  if (isMaster) {
#endif
  if (!SDDS_StartTable(&SDDS_aperture, n_left + n_right)) {
    SDDS_SetError("Unable to start SDDS table (do_aperture_search)");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors | SDDS_EXIT_PrintErrors);
  }
  SDDS_SetParameters(&SDDS_aperture, SDDS_SET_BY_INDEX | SDDS_PASS_BY_VALUE, IP_STEP, control->i_step, -1);
  if (control->n_elements_to_vary) {
    for (ix = 0; ix < control->n_elements_to_vary; ix++)
      if (!SDDS_SetParameters(&SDDS_aperture, SDDS_SET_BY_INDEX | SDDS_PASS_BY_VALUE, ix + N_PARAMETERS,
                              control->varied_quan_value[ix], -1))
        break;
  }
  if (SDDS_NumberOfErrors()) {
    SDDS_SetError("Problem setting SDDS parameter values (do_aperture_search)");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors | SDDS_EXIT_PrintErrors);
  }

  for (iy = 0; iy < n_left; iy++)
    if (!SDDS_SetRowValues(&SDDS_aperture, SDDS_SET_BY_INDEX | SDDS_PASS_BY_VALUE, iy,
                           IC_X, xy_left[iy][0], IC_Y, xy_left[iy][1], -1)) {
      SDDS_SetError("Problem setting SDDS row values (do_aperture_search)");
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors | SDDS_EXIT_PrintErrors);
    }
  for (iy = 0; iy < n_right; iy++) {
    if (!SDDS_SetRowValues(&SDDS_aperture, SDDS_SET_BY_INDEX | SDDS_PASS_BY_VALUE, iy + n_left,
                           IC_X, xy_right[n_right - iy - 1][0], IC_Y, xy_right[n_right - iy - 1][1], -1)) {
      SDDS_SetError("Problem setting SDDS row values (do_aperture_search)");
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors | SDDS_EXIT_PrintErrors);
    }
  }
  if (!SDDS_WriteTable(&SDDS_aperture)) {
    SDDS_SetError("Problem writing SDDS table (do_aperture_search)");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors | SDDS_EXIT_PrintErrors);
  }
  if (!inhibitFileSync)
    SDDS_DoFSync(&SDDS_aperture);
#if USE_MPI
}
#endif
free_czarray_2d((void **)coord, 1, totalPropertiesPerParticle);
free_czarray_2d((void **)xy_left, ny, 2);
free_czarray_2d((void **)xy_right, ny, 2);

log_exit("do_aperture_search_sp");
return (1);
}

void finish_aperture_search(
  RUN *run,
  VARY *control,
  ERRORVAL *errcon,
  LINE_LIST *beamline) {
  if (output) {
    if (SDDS_IsActive(&SDDS_aperture) && !SDDS_Terminate(&SDDS_aperture)) {
      SDDS_SetError("Problem terminating SDDS output (finish_aperture_search)");
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors | SDDS_EXIT_PrintErrors);
    }
    if (fpSearchOutput) {
      fclose(fpSearchOutput);
      fpSearchOutput = NULL;
    }
  }
}

#ifdef USE_MPI
int comp_index(const void *idx1, const void *idx2) {
  double a1 = *((long *)idx1), a2 = *((long *)idx2);

  if (a1 < a2)
    return -1;
  else if (a1 > a2)
    return 1;
  else
    return 0;
}
#endif

/* line search routine */

long do_aperture_search_line(
  RUN *run,
  VARY *control,
  double *referenceCoord,
  ERRORVAL *errcon,
  LINE_LIST *beamline,
  long lines,
  double *returnValue) {
  double **coord;
  double x0, y0, dx, dy;
  double p_central;
  long index, split, nSteps;
  long effort, n_trpoint, line;
  double xSurvived, ySurvived, area, dtheta;
  double orbit[6] = {0, 0, 0, 0, 0, 0};
  double *dxFactor, *dyFactor;
  double *xLimit, *yLimit;
  double xLost, yLost, sLost;
  long originStable = 0;
#if USE_MPI
  long break_index; /* The index for which the particle is lost on each CPU */
  long last_index;  /* The index for which the last particle is survived */
#endif

#if USE_MPI
  return do_aperture_search_line_p(run, control, referenceCoord, errcon, beamline, lines, returnValue);
#endif

  coord = (double **)czarray_2d(sizeof(**coord), 1, totalPropertiesPerParticle);

  dxFactor = tmalloc(sizeof(*dxFactor) * lines);
  dyFactor = tmalloc(sizeof(*dyFactor) * lines);
  xLimit = tmalloc(sizeof(*xLimit) * lines);
  yLimit = tmalloc(sizeof(*yLimit) * lines);

  switch (lines) {
  case 1:
    dxFactor[0] = 1;
    dyFactor[0] = 1;
    break;
  case 2:
    dxFactor[0] = -1;
    dyFactor[0] = 1;
    dxFactor[1] = 1;
    dyFactor[1] = 1;
    break;
  case 3:
    dxFactor[0] = 0;
    dyFactor[0] = 1;
    dxFactor[1] = dyFactor[1] = 1 / sqrt(2);
    dxFactor[2] = 1;
    dyFactor[2] = 0;
    break;
  default:
    if (full_plane == 0)
      dtheta = PI / (lines - 1);
    else
      dtheta = PIx2 / (lines - 1);
    for (line = 0; line < lines; line++) {
      dxFactor[line] = sin(-PI / 2 + dtheta * line);
      dyFactor[line] = cos(-PI / 2 + dtheta * line);
      if (fabs(dxFactor[line]) < 1e-6)
        dxFactor[line] = 0;
      if (fabs(dyFactor[line]) < 1e-6)
        dyFactor[line] = 0;
    }
    break;
  }

  effort = 0;
  xSurvived = ySurvived = -1;
  dx = dy = 0;
  if (offset_by_orbit) {
    /* N.B.: for an off-momentum orbit that is created with an initial
     * MALIGN element, the momentum offset will not appear in the
     * referenceCoord array.  So this works if the user sets ON_PASS=0
     * for the MALIGN.
     */
    memcpy(orbit, referenceCoord, sizeof(*referenceCoord) * 6);
    /*
      fprintf(stderr, "offseting by orbit: %e, %e, %e, %e, %e, %e \n",
      orbit[0], orbit[1], orbit[2],
	    orbit[3], orbit[4], orbit[5]);
    */
  }

  if (verbosity >= 1) {
    printf("** Starting %ld-line aperture search\n", lines);
    fflush(stdout);
  }

  if (output) {
    if (!SDDS_StartTable(&SDDS_aperture, lines)) {
      SDDS_SetError("Unable to start SDDS table (do_aperture_search)");
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors | SDDS_EXIT_PrintErrors);
    }
    SDDS_SetParameters(&SDDS_aperture, SDDS_SET_BY_INDEX | SDDS_PASS_BY_VALUE, IP_STEP, control->i_step, -1);
    if (control->n_elements_to_vary) {
      for (index = 0; index < control->n_elements_to_vary; index++)
        if (!SDDS_SetParameters(&SDDS_aperture, SDDS_SET_BY_INDEX | SDDS_PASS_BY_VALUE, index + N_PARAMETERS,
                                control->varied_quan_value[index], -1))
          break;
    }
    if (SDDS_NumberOfErrors()) {
      SDDS_SetError("Problem setting SDDS parameter values (do_aperture_search)");
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors | SDDS_EXIT_PrintErrors);
    }
  }

  for (line = 0; line < lines; line++) {
    if (dxFactor[line] > 0)
      xSurvived = ySurvived = -1;
    else
      xSurvived = -(ySurvived = -1);
    printf("* Searching line %ld\n", line);
    fflush(stdout);
    xLost = yLost = sLost = DBL_MAX;
    for (split = 0; split <= n_splits; split++) {
      if (split == 0) {
        if (!slope_mode) {
          dx = xmax / (nx - 1) * dxFactor[line];
          dy = ymax / (nx - 1) * dyFactor[line];
        } else {
          dx = xpmax / (nx - 1) * dxFactor[line];
          dy = ypmax / (nx - 1) * dyFactor[line];
        }
        x0 = y0 = 0;
        nSteps = nx;
      } else {
        x0 = xSurvived;
        y0 = ySurvived;
        dx *= split_fraction;
        dy *= split_fraction;
        x0 += dx;
        y0 += dy;
        nSteps = 1 / split_fraction - 0.5;
        if (nSteps < 1)
          nSteps = 1;
        if (verbosity >= 1) {
          printf("divided search interval to %e, %e\n", dx, dy);
          fflush(stdout);
        }
      }
#if USE_MPI
      break_index = nSteps; /* Initialization, no particle has been lost */
      last_index = -1;
      if ((verbosity >= 1) && isMaster && split == 0) {
        if (nSteps <= n_processors) {
          char buffer[16384];
          snprintf(buffer, 16384, ". Using %ld or fewer CPUs is suggested.", nSteps - 1);
          printWarning("find_aperture: please reduce the number of CPUs or use a finer grid to avoid wasting resources", buffer);
        }
      }
#endif
      for (index = 0; index < nSteps; index++) {
#if USE_MPI
        if (myid == index % n_processors || !originStable) { /* All CPUs will track the first point */
#endif
          if (index != 0 || split != 0 || !originStable) {
            memcpy(coord[0], orbit, sizeof(*orbit) * 6);
            if (!slope_mode) {
              coord[0][0] = index * dx + x0 + orbit[0];
              coord[0][2] = index * dy + y0 + orbit[2];
            } else {
              coord[0][1] = index * dx + x0 + orbit[1];
              coord[0][3] = index * dy + y0 + orbit[3];
            }

            p_central = run->p_central;
            n_trpoint = 1;
            if (do_tracking(NULL, coord, n_trpoint, &effort, beamline, &p_central,
                            NULL, NULL, NULL, NULL, run, control->i_step,
                            SILENT_RUNNING, control->n_passes, 0, NULL, NULL, NULL, NULL, NULL) != 1) {
              if (verbosity >= 2) {
                if (!slope_mode)
                  printf("particle lost for x=%e, y=%e\n", index * dx + x0, index * dy + y0);
                else
                  printf("particle lost for xp=%e, yp=%e\n", index * dx + x0, index * dy + y0);
                fflush(stdout);
              }
              xLost = coord[0][0];
              yLost = coord[0][2];
              sLost = coord[0][4];
#if USE_MPI
              break_index = index;
#endif
              break;
            }
          }
          if (index == 0 && split == 0)
            originStable = 1;
          if (verbosity >= 2) {
            if (!slope_mode)
              printf("particle survived for x=%e, y=%e\n", x0 + index * dx, y0 + index * dy);
            else
              printf("particle survived for xp=%e, yp=%e\n", x0 + index * dx, y0 + index * dy);
            fflush(stdout);
          }
          if (dxFactor[line]) {
            if ((dxFactor[line] > 0 && (xSurvived < (x0 + index * dx))) ||
                (dxFactor[line] < 0 && (xSurvived > (x0 + index * dx)))) {
              xSurvived = x0 + index * dx;
              ySurvived = y0 + index * dy;
            }
          } else {
            if (ySurvived < (y0 + index * dy)) {
              xSurvived = x0 + index * dx;
              ySurvived = y0 + index * dy;
            }
          }

#if USE_MPI
        }
#endif
      }

#if USE_MPI
      /* find the global extreme value and save it on master */
      if (USE_MPI) {
        long *index_array;
        double *xLost_array, *yLost_array, *sLost_array;

        if (!(index_array = (long *)malloc(n_processors * sizeof(*index_array))) ||
            !(xLost_array = (double *)malloc(n_processors * sizeof(*xLost_array))) ||
            !(yLost_array = (double *)malloc(n_processors * sizeof(*yLost_array))) ||
            !(sLost_array = (double *)malloc(n_processors * sizeof(*sLost_array))))
          bombElegant("memory allocation failure gathering data from processors", NULL);

        MPI_Gather(&break_index, 1, MPI_LONG, index_array, 1, MPI_LONG, 0, MPI_COMM_WORLD);
        MPI_Gather(&xLost, 1, MPI_DOUBLE, xLost_array, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Gather(&yLost, 1, MPI_DOUBLE, yLost_array, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Gather(&sLost, 1, MPI_DOUBLE, sLost_array, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        if (isMaster) {
          long imin, imax;
          imin = 0;
          index_min_max_long(&imin, &imax, index_array, n_processors);
          if (index_array[imin] > 0) {
            if (index_array[imin] != nSteps)      /* A particle is lost */
              last_index = index_array[imin] - 1; /* Get the index for the survived particle */
            else                                  /* The last one will be the index we need */
              last_index = nSteps - 1;
          } else /* A special case where only the first particle survived */
            last_index = 0;

          xLost = xLost_array[imin];
          yLost = yLost_array[imin];
          sLost = sLost_array[imin];
          if ((last_index != -1) && (last_index != 0)) {
            xSurvived = x0 + last_index * dx;
            ySurvived = y0 + last_index * dy;
          }
        }
        free(index_array);
        free(xLost_array);
        free(yLost_array);
        free(sLost_array);

        /* Broadcasts the last survived particle coordincates to all the processors */
        MPI_Bcast(&xSurvived, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&ySurvived, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&xLost, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&yLost, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&sLost, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      }
#endif
      if (verbosity >= 1) {
        printf("Sweep done, particle survived up to x=%e, y=%e\n", xSurvived, ySurvived);
        fflush(stdout);
      }
    }
    xLimit[line] = xSurvived;
    yLimit[line] = ySurvived;
#if USE_MPI
    if (isMaster)
#endif
      if (output) {
        if (!SDDS_SetRowValues(&SDDS_aperture, SDDS_SET_BY_INDEX | SDDS_PASS_BY_VALUE, line,
                               IC_X, xSurvived, IC_Y, ySurvived,
                               IC_XLOST, xLost, IC_YLOST, yLost, IC_SLOST, sLost, -1)) {
          SDDS_SetError("Problem setting SDDS row values (do_aperture_search)");
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors | SDDS_EXIT_PrintErrors);
        }
      }
  }

  area = 0;
#if USE_MPI
  if (isMaster) {
#endif
    if (originStable && lines > 1) {
      /* compute the area */

      area = trimApertureSearchResult(lines, xLimit, yLimit, dxFactor, dyFactor, full_plane);
    }
    *returnValue = area;

    if (output) {
      if (!SDDS_SetColumn(&SDDS_aperture, SDDS_SET_BY_INDEX, xLimit, lines, IC_XC) ||
          !SDDS_SetColumn(&SDDS_aperture, SDDS_SET_BY_INDEX, yLimit, lines, IC_YC) ||
          !SDDS_SetParameters(&SDDS_aperture, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE,
                              "Area", area, NULL)) {
        SDDS_SetError("Problem setting parameters values in SDDS table (do_aperture_search)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors | SDDS_EXIT_PrintErrors);
      }
      if (control->n_elements_to_vary) {
        long i;
        for (i = 0; i < control->n_elements_to_vary; i++)
          if (!SDDS_SetParameters(&SDDS_aperture, SDDS_SET_BY_INDEX | SDDS_PASS_BY_VALUE, i + N_PARAMETERS,
                                  control->varied_quan_value[i], -1))
            break;
      }

      if (!SDDS_WriteTable(&SDDS_aperture)) {
        SDDS_SetError("Problem writing SDDS table (do_aperture_search)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors | SDDS_EXIT_PrintErrors);
      }
      if (!inhibitFileSync)
        SDDS_DoFSync(&SDDS_aperture);
    }
#if USE_MPI
  }
#endif

  free_czarray_2d((void **)coord, 1, totalPropertiesPerParticle);
  free(dxFactor);
  free(dyFactor);
  free(xLimit);
  free(yLimit);

  return (1);
}

#if USE_MPI

/* Line aperture search that makes use of larger parallel resources by tracking all particles for all
   lines at the same time (if enough processors are available
*/

long do_aperture_search_line_p(
  RUN *run,
  VARY *control,
  double *referenceCoord,
  ERRORVAL *errcon,
  LINE_LIST *beamline,
  long lines,
  double *returnValue) {
  double **coord;
  double p_central;
  long index, split, step, nSteps;
  long effort, n_trpoint, line;
  double area, dtheta;
  double orbit[6] = {0, 0, 0, 0, 0, 0};
  double *x0, *y0, *dx, *dy, *dxFactor, *dyFactor;
  double *xLimit, *yLimit, *xLost, *yLost, *sLost;
  double **xLost2, **yLost2, **sLost2;
  double **survived, **buffer;
  long maxSteps;
#  if DEBUG
  FILE *fpd, *fpd2, *fpd3;
  char s[100];
  sprintf(s, "pda-%03d.sdds", myid);
  fpd = fopen(s, "w");
  fprintf(fpd, "SDDS1\n&column name=x type=double &end\n&column name=y type=double &end\n&column name=survived type=short &end\n");
  fprintf(fpd, "&column name=split type=short &end\n&column name=line type=short &end\n&column name=step type=short &end\n");
  fprintf(fpd, "&column name=index type=short &end\n&column name=ID type=short &end\n&data mode=ascii no_row_counts=1 &end\n");
  sprintf(s, "pda-%03d.txt", myid);
  fpd2 = fopen(s, "w");
  fpd3 = fopen("pda.sdds", "w");
  fprintf(fpd3, "SDDS1\n&column name=x type=double &end\n&column name=y type=double &end\n&column name=survived type=short &end\n");
  fprintf(fpd3, "&column name=split type=short &end\n&column name=line type=short &end\n&column name=step type=short &end\n");
  fprintf(fpd3, "&column name=index type=short &end\n&column name=ID type=short &end\n&data mode=ascii no_row_counts=1 &end\n");
#  endif

  maxSteps = 1 / split_fraction + 0.5;
  if (nx > maxSteps)
    maxSteps = nx;

  coord = (double **)czarray_2d(sizeof(**coord), 1, totalPropertiesPerParticle);

  /* These arrays are used to store results for each line and each step on the line */
  survived = (double **)czarray_2d(sizeof(**survived), lines, maxSteps);
  xLost2 = (double **)czarray_2d(sizeof(**xLost2), lines, maxSteps);
  yLost2 = (double **)czarray_2d(sizeof(**yLost2), lines, maxSteps);
  sLost2 = (double **)czarray_2d(sizeof(**sLost2), lines, maxSteps);
  /* buffer for sharing data among processors */
  buffer = (double **)czarray_2d(sizeof(**buffer), lines, maxSteps);

  /* Some of these were just scalars in the previous version.  We use arrays to avoid recomputing values */
  dx = tmalloc(sizeof(*dx) * lines);
  dy = tmalloc(sizeof(*dy) * lines);
  x0 = tmalloc(sizeof(*x0) * lines);
  y0 = tmalloc(sizeof(*y0) * lines);
  xLost = tmalloc(sizeof(*xLost) * lines);
  yLost = tmalloc(sizeof(*yLost) * lines);
  sLost = tmalloc(sizeof(*sLost) * lines);
  dxFactor = tmalloc(sizeof(*dxFactor) * lines);
  dyFactor = tmalloc(sizeof(*dyFactor) * lines);
  xLimit = calloc(sizeof(*xLimit), lines);
  yLimit = calloc(sizeof(*yLimit), lines);

  for (line = 0; line < lines; line++)
    xLost[line] = yLost[line] = sLost[line] = DBL_MAX;

  switch (lines) {
  case 1:
    dxFactor[0] = 1;
    dyFactor[0] = 1;
    break;
  case 2:
    dxFactor[0] = -1;
    dyFactor[0] = 1;
    dxFactor[1] = 1;
    dyFactor[1] = 1;
    break;
  case 3:
    dxFactor[0] = 0;
    dyFactor[0] = 1;
    dxFactor[1] = dyFactor[1] = 1 / sqrt(2);
    dxFactor[2] = 1;
    dyFactor[2] = 0;
    break;
  default:
    if (full_plane == 0)
      dtheta = PI / (lines - 1);
    else
      dtheta = PIx2 / (lines - 1);
    for (line = 0; line < lines; line++) {
      dxFactor[line] = sin(-PI / 2 + dtheta * line);
      dyFactor[line] = cos(-PI / 2 + dtheta * line);
      if (fabs(dxFactor[line]) < 1e-6)
        dxFactor[line] = 0;
      if (fabs(dyFactor[line]) < 1e-6)
        dyFactor[line] = 0;
    }
    break;
  }

  effort = 0;
  if (offset_by_orbit) {
    /* N.B.: for an off-momentum orbit that is created with an initial
     * MALIGN element, the momentum offset will not appear in the
     * referenceCoord array.  So this works if the user sets ON_PASS=0
     * for the MALIGN.
     */
    memcpy(orbit, referenceCoord, sizeof(*referenceCoord) * 6);
  }

  if (verbosity >= 1) {
    printf("** Starting %ld-line aperture search\n", lines);
    fflush(stdout);
  }

  if (isMaster && output) {
    if (!SDDS_StartTable(&SDDS_aperture, lines)) {
      SDDS_SetError("Unable to start SDDS table (do_aperture_search)");
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors | SDDS_EXIT_PrintErrors);
    }
    SDDS_SetParameters(&SDDS_aperture, SDDS_SET_BY_INDEX | SDDS_PASS_BY_VALUE, IP_STEP, control->i_step, -1);
    if (control->n_elements_to_vary) {
      for (index = 0; index < control->n_elements_to_vary; index++)
        if (!SDDS_SetParameters(&SDDS_aperture, SDDS_SET_BY_INDEX | SDDS_PASS_BY_VALUE, index + N_PARAMETERS,
                                control->varied_quan_value[index], -1))
          break;
    }
    if (SDDS_NumberOfErrors()) {
      SDDS_SetError("Problem setting SDDS parameter values (do_aperture_search)");
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors | SDDS_EXIT_PrintErrors);
    }
  }

  /* For each split, loop over lines and steps on lines to evaluate which particles survive */
  for (split = 0; split <= n_splits; split++) {
    if (split == 0)
      nSteps = nx;
    else
      nSteps = 1 / split_fraction - 0.5;
    if (nSteps < 1)
      nSteps = 1;
    /* Mark all particles as not surviving */
    for (line = 0; line < lines; line++) {
      for (step = 0; step < nSteps; step++)
        survived[line][step] = 0;
    }
    /* loop over lines.  'index' counts the particles to be tracked */
    for (index = line = 0; line < lines; line++) {
      /* Set start and delta values for the line */
      if (split == 0) {
        if (!slope_mode) {
          dx[line] = xmax / (nx - 1) * dxFactor[line];
          dy[line] = ymax / (nx - 1) * dyFactor[line];
        } else {
          dx[line] = xpmax / (nx - 1) * dxFactor[line];
          dy[line] = ypmax / (nx - 1) * dyFactor[line];
        }
        x0[line] = y0[line] = 0;
      } else {
        x0[line] = xLimit[line];
        y0[line] = yLimit[line];
        dx[line] *= split_fraction;
        dy[line] *= split_fraction;
        x0[line] += dx[line];
        y0[line] += dy[line];
      }
      /* step along the line */
      for (step = 0; step < nSteps; step++, index++) {
        /* initialize 2d arrays to zero even if we won't track this particle */
        xLost2[line][step] = 0;
        yLost2[line][step] = 0;
        sLost2[line][step] = 0;
        survived[line][step] = 0;
#  if DEBUG
        fprintf(fpd2, "myid=%d, split=%ld, line=%ld, step=%ld, index=%ld, n_processors=%d---",
                myid, split, line, step, index, n_processors);
        fflush(fpd2);
#  endif

        /* decide if we are going to track this particle */
        if (myid != index % n_processors) {
#  if DEBUG
          fprintf(fpd2, "not tracking\n");
          fflush(fpd2);
#  endif
          continue;
        }
#  if DEBUG
        fprintf(fpd2, "tracking---");
        fflush(fpd2);
#  endif

        /* Track a particle */
        memcpy(coord[0], orbit, sizeof(*orbit) * 6);
        if (!slope_mode) {
          coord[0][0] = step * dx[line] + x0[line] + orbit[0];
          coord[0][2] = step * dy[line] + y0[line] + orbit[2];
        } else {
          coord[0][1] = step * dx[line] + x0[line] + orbit[1];
          coord[0][3] = step * dy[line] + y0[line] + orbit[3];
        }
        p_central = run->p_central;
        n_trpoint = 1;
        if (do_tracking(NULL, coord, n_trpoint, &effort, beamline, &p_central,
                        NULL, NULL, NULL, NULL, run, control->i_step,
                        SILENT_RUNNING, control->n_passes, 0, NULL, NULL, NULL, NULL, NULL) != 1) {
          /* Particle lost, so record information */
          xLost2[line][step] = coord[0][0];
          yLost2[line][step] = coord[0][2];
          sLost2[line][step] = coord[0][4];
          survived[line][step] = 0;
#  if DEBUG
          fprintf(fpd2, "lost\n");
          fflush(fpd2);
#  endif
        } else {
          /* Particle survived */
          survived[line][step] = 1;
#  if DEBUG
          fprintf(fpd2, "survived\n");
          fflush(fpd2);
#  endif
        }
#  if DEBUG
        fprintf(fpd, "%le %le %.0f %ld %ld %ld %ld %d\n", step * dx[line] + x0[line], step * dy[line] + y0[line], survived[line][step], split, line, step,
                index, myid);
        fflush(fpd);
#  endif
      }
    }
    /* Wait for all processors to exit the loop */
#  if MPI_DEBUG
    printf("Waiting on barrier after loop over all lines\n");
    fflush(stdout);
#  endif
    MPI_Barrier(MPI_COMM_WORLD);

    /* Sum values of arrays over all processors */
    MPI_Allreduce(survived[0], buffer[0], lines * maxSteps, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    memcpy(survived[0], buffer[0], sizeof(double) * lines * maxSteps);

    MPI_Allreduce(xLost2[0], buffer[0], lines * maxSteps, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    memcpy(xLost2[0], buffer[0], sizeof(double) * lines * maxSteps);

    MPI_Allreduce(yLost2[0], buffer[0], lines * maxSteps, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    memcpy(yLost2[0], buffer[0], sizeof(double) * lines * maxSteps);

    MPI_Allreduce(sLost2[0], buffer[0], lines * maxSteps, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    memcpy(sLost2[0], buffer[0], sizeof(double) * lines * maxSteps);

#  if DEBUG
    for (line = 0; line < lines; line++)
      for (step = 0; step < nSteps; step++) {
        fprintf(fpd3, "%le %le %.0f %ld %ld %ld %ld %d\n", step * dx[line] + x0[line], step * dy[line] + y0[line], survived[line][step], split, line, step,
                index, myid);
        fflush(fpd3);
      }
#  endif

    /* Scan array to determine x and y stability limits */
    for (line = 0; line < lines; line++) {
#  if DEBUG
      printf("Checking line=%ld\n", line);
      fflush(stdout);
#  endif
      if (split == 0) {
        if (!survived[line][0])
          bombElegant("Particle unstable at origin", NULL);
        for (step = 1; step < nSteps; step++) {
          if (!survived[line][step])
            break;
        }
#  if DEBUG
        printf("split=0, lost at step=%ld of %ld\n", step, nSteps);
#  endif
      } else {
        for (step = 0; step < nSteps; step++) {
          if (!survived[line][step])
            break;
        }
#  if DEBUG
        printf("split=%ld, lost at step=%ld of %ld\n", split, step, nSteps);
        fflush(stdout);
#  endif
        if (step == 0)
          continue;
      }
      step--;

#  if DEBUG
      printf("line=%ld, split=%ld, x0=%le, y0=%le, dx=%le, dy=%le\n", line, split,
             x0[line], y0[line], dx[line], dy[line]);
      printf("Particle survived at step=%ld, x=%le, y=%le\n", step, x0[line] + step * dx[line], y0[line] + step * dy[line]);
      fflush(stdout);
#  endif

      if ((dx[line] > 0 && xLimit[line] < (step * dx[line] + x0[line])) ||
          (dx[line] == 0 && yLimit[line] < (step * dy[line] + y0[line])) ||
          (dx[line] < 0 && xLimit[line] > (step * dx[line] + x0[line]))) {
#  if DEBUG
        printf("Change limit from (%le, %le) to (%le, %le)\n",
               xLimit[line], yLimit[line], step * dx[line] + x0[line], step * dy[line] + y0[line]);
        fflush(stdout);
#  endif
        xLimit[line] = step * dx[line] + x0[line];
        yLimit[line] = step * dy[line] + y0[line];
        if (step < (nSteps - 1)) {
          xLost[line] = xLost2[line][step + 1];
          yLost[line] = yLost2[line][step + 1];
          sLost[line] = sLost2[line][step + 1];
        }
      }
    }
  }

  if (isMaster) {
    /* Analyze and store results */

    if (output) {
      for (line = 0; line < lines; line++)
        if (!SDDS_SetRowValues(&SDDS_aperture, SDDS_SET_BY_INDEX | SDDS_PASS_BY_VALUE, line,
                               IC_X, xLimit[line], IC_Y, yLimit[line],
                               IC_XLOST, xLost[line], IC_YLOST, yLost[line], IC_SLOST, sLost[line], -1)) {
          SDDS_SetError("Problem setting SDDS row values (do_aperture_search)");
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors | SDDS_EXIT_PrintErrors);
        }
    }

    area = 0;
    if (lines > 1) {
      /* compute the area */
      area = trimApertureSearchResult(lines, xLimit, yLimit, dxFactor, dyFactor, full_plane);
    }
    *returnValue = area;

    if (output) {
      if (!SDDS_SetColumn(&SDDS_aperture, SDDS_SET_BY_INDEX, xLimit, lines, IC_XC) ||
          !SDDS_SetColumn(&SDDS_aperture, SDDS_SET_BY_INDEX, yLimit, lines, IC_YC) ||
          !SDDS_SetParameters(&SDDS_aperture, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE,
                              "Area", area, NULL)) {
        SDDS_SetError("Problem setting parameters values in SDDS table (do_aperture_search)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors | SDDS_EXIT_PrintErrors);
      }
      if (control->n_elements_to_vary) {
        long i;
        for (i = 0; i < control->n_elements_to_vary; i++)
          if (!SDDS_SetParameters(&SDDS_aperture, SDDS_SET_BY_INDEX | SDDS_PASS_BY_VALUE, i + N_PARAMETERS,
                                  control->varied_quan_value[i], -1))
            break;
      }

      if (!SDDS_WriteTable(&SDDS_aperture)) {
        SDDS_SetError("Problem writing SDDS table (do_aperture_search)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors | SDDS_EXIT_PrintErrors);
      }
      if (!inhibitFileSync)
        SDDS_DoFSync(&SDDS_aperture);
    }
  }

  free_czarray_2d((void **)coord, 1, totalPropertiesPerParticle);
  free_czarray_2d((void **)survived, lines, maxSteps);
  free_czarray_2d((void **)xLost2, lines, maxSteps);
  free_czarray_2d((void **)yLost2, lines, maxSteps);
  free_czarray_2d((void **)sLost2, lines, maxSteps);
  free_czarray_2d((void **)buffer, lines, maxSteps);

  free(dx);
  free(dy);
  free(x0);
  free(y0);
  free(xLost);
  free(yLost);
  free(dxFactor);
  free(dyFactor);
  free(xLimit);
  free(yLimit);

  return (1);
}

#endif
