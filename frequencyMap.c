/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: frequencyMap.c
 * purpose: Do frequency map tracking and analysis.
 *          See file frequencyMap.nl for input parameters.
 *
 * Michael Borland, 2004
 */
#include "mdb.h"
#include "track.h"
#include "frequencyMap.h"

#define IC_X 0
#define IC_Y 1
#define IC_NUX 2
#define IC_NUY 3
#define IC_DNUX 4
#define IC_DNUY 5
#define N_COLUMNS 6
static SDDS_DEFINITION column_definition[N_COLUMNS] = {
    {"x", "&column name=x, symbol=x, units=m, type=double &end"},
    {"y", "&column name=y, symbol=y, units=m, type=double &end"},
    {"nux", "&column name=nux, symbol=$gn$r$bx$n, type=double &end"},
    {"nuy", "&column name=nuy, symbol=$gn$r$by$n, type=double &end"},
    {"dnux", "&column name=dnux, symbol=$gDn$r$bx$n, type=double &end"},
    {"dnuy", "&column name=dnuy, symbol=$gDn$r$by$n, type=double &end"},
    } ;

#define IP_STEP 0
#define N_PARAMETERS 1
static SDDS_DEFINITION parameter_definition[N_PARAMETERS] = {
    {"Step", "&parameter name=Step, type=long, description=\"Simulation step\" &end"},
    } ;

static SDDS_DATASET SDDS_fmap;

void setupFrequencyMap(
    NAMELIST_TEXT *nltext,
    RUN *run,
    VARY *control
    )
{
  /* process namelist input */
  set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
  set_print_namelist_flags(0);
  process_namelist(&frequency_map, nltext);
  print_namelist(stdout, &frequency_map);
  
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
  
  output = compose_filename(output, run->rootname);
  SDDS_ElegantOutputSetup(&SDDS_fmap, output, SDDS_BINARY, 1, "frequency map analysis",
                          run->runfile, run->lattice, parameter_definition, N_PARAMETERS,
                          column_definition, N_COLUMNS, "setup_frequencyMap", SDDS_EOS_NEWFILE);
  
  if (control->n_elements_to_vary) 
    if (!SDDS_DefineSimpleParameters(&SDDS_fmap, control->n_elements_to_vary,
                                     control->varied_quan_name, control->varied_quan_unit, SDDS_DOUBLE)) {
      SDDS_SetError("Unable to define additional SDDS parameters (setup_aperture_search)");
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
  
  if (!SDDS_WriteLayout(&SDDS_fmap)) {
    SDDS_SetError("Unable to write SDDS layout for aperture search");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  
}


long doFrequencyMap(
                    RUN *run,
                    VARY *control,
                    double *referenceCoord,
                    ERRORVAL *errcon,
                    LINE_LIST *beamline
                    )
{
  double firstTune[2], secondTune[2], startingCoord[6], endingCoord[6];
  double dx, dy, x, y;
  long ix, iy, ip;
  
  if (!SDDS_StartPage(&SDDS_fmap, nx*ny) || 
      !SDDS_SetParameters(&SDDS_fmap, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, 0, control->i_step, -1)) {
    SDDS_SetError("Unable to start SDDS page (do_frequencyMap)");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  if (control->n_elements_to_vary) {
    for (ip=0; ip<control->n_elements_to_vary; ip++)
      if (!SDDS_SetParameters(&SDDS_fmap, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, ip+1,
                              control->varied_quan_value[ip], -1)) {
        SDDS_SetError("Unable to start SDDS page (do_frequencyMap)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      }
  }

  dx  = (xmax-xmin)/(nx-1);
  dy = (ymax-ymin)/(ny-1);
  ip = 0;
  for (ix=0; ix<nx; ix++) {
    x = xmin + ix*dx;
    for (iy=0; iy<ny; iy++) {
      y = ymin + iy*dy;
      memcpy(startingCoord, referenceCoord, sizeof(*startingCoord)*6);
      if (!computeTunesFromTracking(firstTune, beamline->matrix, beamline, run,
                                    startingCoord, x, y, control->n_passes/2, 
                                    0, endingCoord)) 
        continue;
      memcpy(startingCoord, endingCoord, sizeof(*startingCoord)*6);
      if (!computeTunesFromTracking(secondTune, beamline->matrix, beamline, run,
                                    startingCoord, 0.0, 0.0, control->n_passes/2, 
                                    0, endingCoord)) 
        continue;
      if (!SDDS_SetRowValues(&SDDS_fmap, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, ip,
                             IC_X, x, IC_Y, y, 
                             IC_NUX, firstTune[0], 
                             IC_NUY, firstTune[1], 
                             IC_DNUX, fabs(secondTune[0]-firstTune[0]), 
                             IC_DNUY, fabs(secondTune[1]-firstTune[1]), -1)) {
        SDDS_SetError("Problem setting SDDS row values (doFrequencyMap)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      }
      ip ++;
      if (verbosity) {
        fprintf(stdout, "Done with particle %ld of %ld\n",
                ip, nx*ny);
        fflush(stdout);
      }
    }
  }
  SDDS_DoFSync(&SDDS_fmap);
  if (!SDDS_WriteTable(&SDDS_fmap)) {
    SDDS_SetError("Problem writing SDDS table (doFrequencyMap)");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  
  return(1);
}

void finishFrequencyMap()
{
  if (SDDS_IsActive(&SDDS_fmap) && !SDDS_Terminate(&SDDS_fmap)) {
    SDDS_SetError("Problem terminating SDDS output (finish_aperture_search)");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
}

