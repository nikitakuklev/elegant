/* Copyright 1994 by Michael Borland and Argonne National Laboratory,
 * all rights reserved.
 */
/* file: sdds_support.c
 * purpose: routines for working with SDDS files
 * 
 * Michael Borland, 1993.
 */
#include "mdb.h"
#include "fftpackC.h"
#include "track.h"
#include "SDDS.h"

void SDDS_ElegantOutputSetup(SDDS_TABLE *SDDS_table, char *filename, long mode, long lines_per_row,
                             char *contents, char *command_file, char *lattice_file, SDDS_DEFINITION *parameter_definition,
                             long n_parameters, SDDS_DEFINITION *column_definition, long n_columns,
                             char *caller, long flags)
{
    static char description[SDDS_MAXLINE];
    long i, last_index, index;

    log_entry("SDDS_ElegantOutputSetup");
    last_index = -1;

    if (flags&SDDS_EOS_NEWFILE) {
        if (SDDS_IsActive(SDDS_table)==1) {
            if (!SDDS_Terminate(SDDS_table)) {
                fprintf(stdout, "Unable to terminate SDDS output (%s)\n", caller);
                fflush(stdout);
                SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
                exit(1);
                }
            }
        sprintf(description, "%s--input: %s  lattice: %s", contents, command_file, lattice_file);
        if (!SDDS_InitializeOutput(SDDS_table, mode, 1, description, contents, filename)) {
            SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
            exit(1);
            }
        }

    /* define SDDS parameters */
    for (i=0; i<n_parameters; i++) {
        if (!SDDS_ProcessParameterString(SDDS_table, parameter_definition[i].text, 0) ||
            (index=SDDS_GetParameterIndex(SDDS_table, parameter_definition[i].name))<0) {
            fprintf(stdout, "Unable to define SDDS parameter for %s--string was:\n%s\n", caller,
                    parameter_definition[i].text);
            fflush(stdout);
            SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
            exit(1);
            }
        if (last_index==-1 && flags&SDDS_EOS_NEWFILE && index!=0)
            fprintf(stdout, "\7\7\7WARNING: first defined parameter index for SDDS file %s is not 0--this will probably cause unexpected results\n", filename);
            fflush(stdout);
        if (last_index!=-1 && index!=(last_index+1))
            fprintf(stdout, "\7\7\7WARNING: parameter indices for SDDS file %s are not sequential--this will probably cause unexpected results\n", filename);
            fflush(stdout);
        last_index = index;
        }

    last_index = -1;
    /* define SDDS columns */
    for (i=0; i<n_columns; i++) {
        if (!SDDS_ProcessColumnString(SDDS_table, column_definition[i].text, 0) ||
            (index=SDDS_GetColumnIndex(SDDS_table, column_definition[i].name))<0) {
            fprintf(stdout, "Unable to define SDDS column for %s--string was:\n%s\n", caller,
                    column_definition[i].text);
            fflush(stdout);
            SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
            exit(1);
            }
        if (last_index==-1 && flags&SDDS_EOS_NEWFILE && index!=0)
            fprintf(stdout, "\7\7\7WARNING: first defined column index for SDDS file %s is not 0--this will probably cause unexpected results\n", filename);
            fflush(stdout);
        if (last_index!=-1 && index!=(last_index+1))
            fprintf(stdout, "\7\7\7WARNING: column indices for SDDS file %s are not sequential--this will probably cause unexpected results\n", filename);
            fflush(stdout);
        last_index = index;
        }

    if (flags&SDDS_EOS_COMPLETE) {
        if (!SDDS_WriteLayout(SDDS_table)) {
            fprintf(stdout, "Unable to write SDDS layout for file %s (%s)\n", filename, caller);
            fflush(stdout);
            SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
            exit(1);
            }
        }

    log_exit("SDDS_ElegantOutputSetup");
    }


#define STANDARD_PARAMETERS 1
static SDDS_DEFINITION standard_parameter[STANDARD_PARAMETERS] = {
    {"Step", "&parameter name=Step, type=long, description=\"Simulation step\" &end"},
    } ;

#define ELEMENT_COLUMNS 4
static SDDS_DEFINITION element_column[ELEMENT_COLUMNS] = {
    {"s", "&column name=s, units=m, type=double, description=\"Distance\" &end"},
    {"ElementName", "&column name=ElementName, type=string, description=\"Element name\", format_string=%10s &end"},
    {"ElementOccurence", 
         "&column name=ElementOccurence, type=long, description=\"Occurence of element\", format_string=%6ld &end"},
    {"ElementType", "&column name=ElementType, type=string, description=\"Element-type name\", format_string=%10s &end"},
    } ;

#define PHASE_SPACE_COLUMNS 7
static SDDS_DEFINITION phase_space_column[PHASE_SPACE_COLUMNS] = {
    {"x", "&column name=x, units=m, type=double &end"},
    {"xp", "&column name=xp, symbol=\"x'\", type=double &end"},
    {"y", "&column name=y, units=m, type=double &end"},
    {"yp", "&column name=yp, symbol=\"y'\", type=double &end"},
    {"t", "&column name=t, units=s, type=double &end"},
    {"p", "&column name=p, units=\"m$be$nc\", type=double &end"},
    {"particleID", "&column name=particleID, type=long &end"},
    } ;

#define PHASE_SPACE_PARAMETERS 2
static SDDS_DEFINITION phase_space_parameter[PHASE_SPACE_PARAMETERS] = {
  {"Step", "&parameter name=Step, type=long, description=\"Simulation step\" &end"},
  {"pCentral", "&parameter name=pCentral, symbol=\"p$bcen$n\", units=\"m$be$nc\", type=double &end"},
};

 

  
void SDDS_PhaseSpaceSetup(SDDS_TABLE *SDDS_table, char *filename, long mode, long lines_per_row, char *contents, 
                          char *command_file, char *lattice_file, char *caller)
{
    log_entry("SDDS_PhaseSpaceSetup");
    SDDS_ElegantOutputSetup(SDDS_table, filename, mode, lines_per_row, contents, command_file, lattice_file,
                            phase_space_parameter, PHASE_SPACE_PARAMETERS, phase_space_column, PHASE_SPACE_COLUMNS,
                            caller, SDDS_EOS_NEWFILE|SDDS_EOS_COMPLETE);
    log_exit("SDDS_PhaseSpaceSetup");
    }


#define BEAM_LOSS_COLUMNS 7
static SDDS_DEFINITION beam_loss_column[PHASE_SPACE_COLUMNS] = {
    {"x", "&column name=x, units=m, type=double &end"},
    {"xp", "&column name=xp, symbol=\"x'\", type=double &end"},
    {"y", "&column name=y, units=m, type=double &end"},
    {"yp", "&column name=yp, symbol=\"y'\", type=double &end"},
    {"s", "&column name=s, units=m, type=double &end"},
    {"p", "&column name=p, units=\"m$be$nc\", type=double &end"},
    {"particleID", "&column name=particleID, type=long &end"},
    } ;

void SDDS_BeamLossSetup(SDDS_TABLE *SDDS_table, char *filename, long mode, long lines_per_row, char *contents, 
                          char *command_file, char *lattice_file, char *caller)
{
    log_entry("SDDS_BeamLossSetup");
    SDDS_ElegantOutputSetup(SDDS_table, filename, mode, lines_per_row, contents, command_file, lattice_file,
                            standard_parameter, STANDARD_PARAMETERS, beam_loss_column, BEAM_LOSS_COLUMNS,
                            caller, SDDS_EOS_NEWFILE|SDDS_EOS_COMPLETE);
    log_exit("SDDS_BeamLossSetup");
    }


#define CENTROID_COLUMNS 8
static SDDS_DEFINITION centroid_column[CENTROID_COLUMNS] = {
    {"Cx", "&column name=Cx, symbol=\"<x>\", units=m, type=double &end"},
    {"Cxp", "&column name=Cxp, symbol=\"<x'>\", type=double &end"},
    {"Cy", "&column name=Cy, symbol=\"<y>\", units=m, type=double &end"},
    {"Cyp", "&column name=Cyp, symbol=\"<y'>\", type=double &end"},
    {"Cs", "&column name=Cs, symbol=\"<s>\", units=m, type=double &end"},
    {"Cdelta", "&column name=Cdelta, symbol=\"<$gd$r>\", type=double &end"},
    {"Particles", "&column name=Particles, description=\"Number of particles\", type=long &end"},
    {"pCentral", "&column name=pCentral, symbol=\"p$bcen$n\", units=\"m$be$nc\", type=double &end"},
    } ;

void SDDS_CentroidOutputSetup(SDDS_TABLE *SDDS_table, char *filename, long mode, long lines_per_row, char *contents, 
                          char *command_file, char *lattice_file, char *caller)
{
    log_entry("SDDS_CentroidOutputSetup");
    SDDS_ElegantOutputSetup(SDDS_table, filename, mode, lines_per_row, contents, command_file, lattice_file,
                            standard_parameter, STANDARD_PARAMETERS, element_column, ELEMENT_COLUMNS,
                            caller, SDDS_EOS_NEWFILE);
    SDDS_ElegantOutputSetup(SDDS_table, NULL, 0, 0, NULL, NULL, NULL, NULL, 0, 
                            centroid_column, CENTROID_COLUMNS, caller, SDDS_EOS_COMPLETE);
    log_exit("SDDS_CentroidOutputSetup");
    }

#define SIGMA_MATRIX_COLUMNS 35
static SDDS_DEFINITION sigma_matrix_column[SIGMA_MATRIX_COLUMNS] = {
    {"s1",    "&column name=s1, symbol=\"$gs$r$b1$n\", units=m, type=double &end"},
    {"s12",    "&column name=s12, symbol=\"$gs$r$b12$n\", units=m, type=double &end"},
    {"s13",    "&column name=s13, symbol=\"$gs$r$b13$n\", units=\"m$a2$n\", type=double &end"},
    {"s14",    "&column name=s14, symbol=\"$gs$r$b14$n\", units=m, type=double &end"},
    {"s15",    "&column name=s15, symbol=\"$gs$r$b15$n\", units=\"m$a2$n\", type=double &end"},
    {"s16",    "&column name=s16, symbol=\"$gs$r$b16$n\", units=m, type=double &end"},
    {"s2",    "&column name=s2, symbol=\"$gs$r$b2$n\", type=double &end"},
    {"s23",    "&column name=s23, symbol=\"$gs$r$b23$n\", units=m, type=double &end"},
    {"s24",    "&column name=s24, symbol=\"$gs$r$b24$n\", type=double &end"},
    {"s25",    "&column name=s25, symbol=\"$gs$r$b25$n\", units=m, type=double &end"},
    {"s26",    "&column name=s26, symbol=\"$gs$r$b26$n\", type=double &end"},
    {"s3",    "&column name=s3, symbol=\"$gs$r$b3$n\", units=m, type=double &end"},
    {"s34",    "&column name=s34, symbol=\"$gs$r$b34$n\", units=m, type=double &end"},
    {"s35",    "&column name=s35, symbol=\"$gs$r$b35$n\", units=\"m$a2$n\", type=double &end"},
    {"s36",    "&column name=s36, symbol=\"$gs$r$b36$n\", units=m, type=double &end"},
    {"s4",    "&column name=s4, symbol=\"$gs$r$b4$n\", type=double &end"},
    {"s45",    "&column name=s45, symbol=\"$gs$r$b45$n\", units=m, type=double &end"},
    {"s46",    "&column name=s46, symbol=\"$gs$r$b46$n\", type=double &end"},
    {"s5",    "&column name=s5, symbol=\"$gs$r$b5$n\", units=m, type=double &end"},
    {"s56",    "&column name=s56, symbol=\"$gs$r$b56$n\", units=m, type=double &end"},
    {"s6",    "&column name=s6, symbol=\"$gs$r$b6$n\", type=double &end"},
    {"ma1",    "&column name=ma1, symbol=\"max$sb$ex$sb$e\", units=m, type=double &end"},
    {"ma2",    "&column name=ma2, symbol=\"max$sb$ex'$sb$e\", type=double &end"},
    {"ma3",    "&column name=ma3, symbol=\"max$sb$ey$sb$e\", units=m, type=double &end"},
    {"ma4",    "&column name=ma4, symbol=\"max$sb$ey'$sb$e\", type=double &end"},
    {"Sx",    "&column name=Sx, symbol=\"$gs$r$bx$n\", units=m, type=double &end"},
    {"Sxp",    "&column name=Sxp, symbol=\"$gs$r$bx'$n\", type=double &end"},
    {"Sy",    "&column name=Sy, symbol=\"$gs$r$by$n\", units=m, type=double &end"},
    {"Syp",    "&column name=Syp, symbol=\"$gs$r$by'$n\", type=double &end"},
    {"Ss",    "&column name=Ss, units=m, symbol=\"$gs$r$bs$n\", type=double &end"},
    {"Sdelta",    "&column name=Sdelta, symbol=\"$gs$bd$n$r\", type=double &end"},
    {"ex", "&column name=ex, symbol=\"$ge$r$bx$n\", units=\"$gp$rm\", type=double &end"},
    {"enx", "&column name=enx, symbol=\"$ge$r$bx,n$n\", type=double, units=\"$gp$rm$be$ncm\"  &end"},
    {"ey", "&column name=ey, symbol=\"$ge$r$by$n\", units=\"$gp$rm\", type=double &end"},
    {"eny", "&column name=eny, symbol=\"$ge$r$by,n$n\", type=double, units=\"$gp$rm$be$ncm\"  &end"},
    } ;

void SDDS_SigmaOutputSetup(SDDS_TABLE *SDDS_table, char *filename, long mode, long lines_per_row,
                           char *command_file, char *lattice_file, char *caller)
{
    log_entry("SDDS_SigmaOutputSetup");
    SDDS_ElegantOutputSetup(SDDS_table, filename, mode, lines_per_row, "sigma matrix", command_file, lattice_file, 
                            standard_parameter, STANDARD_PARAMETERS, element_column, ELEMENT_COLUMNS,
                            caller, SDDS_EOS_NEWFILE);
    SDDS_ElegantOutputSetup(SDDS_table, filename, mode, lines_per_row, NULL, NULL, NULL,
                            NULL, 0, sigma_matrix_column, SIGMA_MATRIX_COLUMNS,
                            caller, SDDS_EOS_COMPLETE);
    log_exit("SDDS_SigmaOutputSetup");
    }


#define WATCH_PARAMETER_MODE_COLUMNS 25
#define WATCH_CENTROID_MODE_COLUMNS 15
static SDDS_DEFINITION watch_parameter_mode_column[WATCH_PARAMETER_MODE_COLUMNS] = {
    {"Step", "&column name=Step, type=long &end"},
    {"Pass", "&column name=Pass, type=long &end"},
    {"Cx", "&column name=Cx, symbol=\"<x>\", units=m, type=double &end"},
    {"Cxp", "&column name=Cxp, symbol=\"<x'>\", type=double &end"},
    {"Cy", "&column name=Cy, symbol=\"<y>\", units=m, type=double &end"},
    {"Cyp", "&column name=Cyp, symbol=\"<y'>\", type=double &end"},
    {"Cs", "&column name=Cs, symbol=\"<s>\", units=m, type=double &end"},
    {"Cdelta", "&column name=Cdelta, symbol=\"<$gd$r>\", type=double &end"},
    {"Ct", "&column name=Ct, symbol=\"<t>\", units=s, type=double &end"},
    {"dCt", "&column name=dCt, symbol=\"$gD$r<t>\", units=s, type=double &end"},
    {"Particles", "&column name=Particles, description=\"Number of particles\", type=long &end"},
    {"Transmission", "&column name=Transmission, description=Transmission, type=double &end"},
    {"pCentral", "&column name=pCentral, symbol=\"p$bcen$n\", units=\"m$be$nc\", type=double &end"},
    {"pAverage", "&column name=pAverage, symbol=\"p$bave$n\", units=\"m$be$nc\", type=double &end"},
    {"KAverage",  "&column name=KAverage, symbol=\"K$bave$n\", units=MeV, type=double &end"},
    {"Sx", "&column name=Sx, symbol=\"$gs$r$bx$n\", units=m, type=double &end"},
    {"Sxp", "&column name=Sxp, symbol=\"$gs$r$bx'$n\", type=double &end"},
    {"Sy", "&column name=Sy, symbol=\"$gs$r$by$n\", units=m, type=double &end"},
    {"Syp", "&column name=Syp, symbol=\"$gs$r$by'$n\", type=double &end"},
    {"Ss", "&column name=Ss, symbol=\"$gs$r$bs$n\", units=m, type=double &end"},
    {"Sdelta", "&column name=Sdelta, symbol=\"$gs$bd$n$r\", type=double &end"},
    {"St", "&column name=St, symbol=\"$gs$r$bt$n\", units=s, type=double &end"},
    {"ex", "&column name=ex, symbol=\"$ge$r$bx$n\", units=\"$gp$rm\", type=double &end"},
    {"ey", "&column name=ey, symbol=\"$ge$r$by$n\", units=\"$gp$rm\", type=double &end"},
    {"el", "&column name=el, symbol=\"$ge$r$bl$n\", units=s, type=double &end"},
    } ;

#define WATCH_POINT_FFT_COLUMNS 4
static SDDS_DEFINITION watch_point_fft_column[WATCH_POINT_FFT_COLUMNS] = {
    {"f", "&column name=f, symbol=\"f/f$bsample$n\", description=\"Normalized frequency\", type=double &end"},
    {"Cx", "&column name=Cx, symbol=\"FFT <x>\", units=m, type=double &end"},
    {"Cy", "&column name=Cy, symbol=\"FFT <y>\", units=m, type=double &end"},
    {"Cdelta", "&column name=Cdelta, symbol=\"FFT <$gd$r>\", type=double &end"},
    } ;

void SDDS_WatchPointSetup(WATCH *watch, long mode, long lines_per_row,
                          char *command_file, char *lattice_file, char *caller, char *qualifier)
{
  char s[100];
  SDDS_TABLE *SDDS_table;
  char *filename;
  long watch_mode, columns;
  
  SDDS_table = &watch->SDDS_table;
  filename = watch->filename;
  watch_mode = watch->mode_code;
  switch (watch_mode) {
    case WATCH_COORDINATES:
    SDDS_ElegantOutputSetup(SDDS_table, filename, mode, lines_per_row, "watch-point phase space",
                            command_file, lattice_file,
                            standard_parameter, STANDARD_PARAMETERS,
                            NULL, 0, caller, SDDS_EOS_NEWFILE);
    columns = 0;
    if (watch->xData) {
      if ((watch->xIndex[0]=SDDS_DefineColumn(SDDS_table, "x", NULL, "m", NULL, NULL, SDDS_DOUBLE, 0))<0 ||
          (watch->xIndex[1]=SDDS_DefineColumn(SDDS_table, "xp", NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0))<0) {
        fprintf(stdout, "Unable to define SDDS columns x and xp for file %s (%s)\n",
                filename, caller);
        fflush(stdout);
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
        exit(1);
      }
    }
    if (watch->yData) {
      if ((watch->yIndex[0]=SDDS_DefineColumn(SDDS_table, "y", NULL, "m", NULL, NULL, SDDS_DOUBLE, 0))<0 ||
          (watch->yIndex[1]=SDDS_DefineColumn(SDDS_table, "yp",  NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0))<0) {
              fprintf(stdout, "Unable to define SDDS columns y and yp for file %s (%s)\n",
                filename, caller);
              fflush(stdout);
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
        exit(1);
      }
    }
    if (watch->longitData) {
      if ((watch->longitIndex[0]
           =SDDS_DefineColumn(SDDS_table, "t", NULL, "s", NULL, NULL, SDDS_DOUBLE, 0))<0 ||
          (watch->longitIndex[1]
           =SDDS_DefineColumn(SDDS_table, "p", NULL, "m$be$nc", NULL, NULL, SDDS_DOUBLE, 0))<0 ||
          (watch->longitIndex[2]
           =SDDS_DefineColumn(SDDS_table, "dt", NULL, "s", NULL, NULL, SDDS_DOUBLE, 0))<0) {
        fprintf(stdout, "Unable to define SDDS columns t, dt, and p for file %s (%s)\n",
                filename, caller);
        fflush(stdout);
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
        exit(1);
      }
    }
    if ((watch->IDIndex = SDDS_DefineColumn(SDDS_table, "particleID", NULL, NULL, NULL, NULL, SDDS_LONG, 0))<0) {
      fprintf(stdout, "Unable to define SDDS columns t, dt, and p for file %s (%s)\n",
              filename, caller);
      fflush(stdout);
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      exit(1);
    }
    
    if (!SDDS_DefineSimpleParameter(SDDS_table, "Pass", NULL, SDDS_LONG) ||
        !SDDS_DefineSimpleParameter(SDDS_table, "Particles", NULL, SDDS_LONG) ||
        !SDDS_DefineSimpleParameter(SDDS_table, "pCentral", "m$be$nc", SDDS_DOUBLE) ||
        !SDDS_DefineSimpleParameter(SDDS_table, "PassLength", "m", SDDS_DOUBLE) ||
        !SDDS_DefineSimpleParameter(SDDS_table, "PassCentralTime", "s", SDDS_DOUBLE)) {
      fprintf(stdout, "Unable define SDDS parameter for file %s (%s)\n", filename, caller);
      fflush(stdout);
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      exit(1);
    }
    if (!SDDS_WriteLayout(SDDS_table)) {
      fprintf(stdout, "Unable to write SDDS layout for file %s (%s)\n", filename, caller);
      fflush(stdout);
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      exit(1);
    }
    break;
  case WATCH_CENTROIDS:
    SDDS_ElegantOutputSetup(SDDS_table, filename, mode, lines_per_row, "watch-point centroids",
                            command_file, lattice_file,
                            standard_parameter, STANDARD_PARAMETERS,
                            watch_parameter_mode_column, WATCH_CENTROID_MODE_COLUMNS,
                            caller, SDDS_EOS_NEWFILE|SDDS_EOS_COMPLETE);
    break;
  case WATCH_PARAMETERS:
    SDDS_ElegantOutputSetup(SDDS_table, filename, mode, lines_per_row, "watch-point parameters", 
                            command_file, lattice_file,
                            standard_parameter, STANDARD_PARAMETERS, 
                            watch_parameter_mode_column, WATCH_PARAMETER_MODE_COLUMNS,
                            caller, SDDS_EOS_NEWFILE|SDDS_EOS_COMPLETE);
    break;
  case WATCH_FFT:
    if (!qualifier)
      watch->window_code = FFT_HANNING;
    else if ((watch->window_code=match_string(qualifier, fft_window_name, N_FFT_WINDOWS, 0))<0) 
      bomb("invalid window mode for WATCH fft", NULL);
    sprintf(s,  "watch-point centroid FFT (%s window)", fft_window_name[watch->window_code]);
    SDDS_ElegantOutputSetup(SDDS_table, filename, mode, lines_per_row, s,
                            command_file, lattice_file, standard_parameter, STANDARD_PARAMETERS, 
                            watch_point_fft_column, WATCH_POINT_FFT_COLUMNS,
                            caller, SDDS_EOS_NEWFILE|SDDS_EOS_COMPLETE);
    break;
  default:
    break;
  }
}

void SDDS_HistogramSetup(HISTOGRAM *histogram, long mode, long lines_per_row,
                          char *command_file, char *lattice_file, char *caller)
{
  SDDS_TABLE *SDDS_table;
  char *filename;
  long column, columns;
  
  SDDS_table = &histogram->SDDS_table;
  filename = histogram->filename;
  SDDS_ElegantOutputSetup(SDDS_table, filename, mode, lines_per_row, "histograms of phase-space coordinates",
                          command_file, lattice_file,
                          standard_parameter, STANDARD_PARAMETERS,
                          NULL, 0, caller, SDDS_EOS_NEWFILE);
  for (column=0; column<7; column++)
    histogram->columnIndex[column][0] = 
      histogram->columnIndex[column][1] = -1;
  columns = 0;
  if (histogram->xData) {
    if ((histogram->columnIndex[0][0]
         =SDDS_DefineColumn(SDDS_table, "x", NULL, "m", NULL, NULL, SDDS_DOUBLE, 0))<0 ||
        (histogram->columnIndex[0][1]
         =SDDS_DefineColumn(SDDS_table, "xFrequency", NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0))<0 ||
        (histogram->columnIndex[1][0]
         =SDDS_DefineColumn(SDDS_table, "xp", NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0))<0 ||
        (histogram->columnIndex[1][1]
         =SDDS_DefineColumn(SDDS_table, "xpFrequency", NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0))<0) {
      fprintf(stdout, "Unable to define SDDS columns x and xp for file %s (%s)\n",
              filename, caller);
      fflush(stdout);
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      exit(1);
    }
    columns += 4;
  }
  if (histogram->yData) {
    if ((histogram->columnIndex[2][0]
         =SDDS_DefineColumn(SDDS_table, "y", NULL, "m", NULL, NULL, SDDS_DOUBLE, 0))<0 ||
        (histogram->columnIndex[2][1]
         =SDDS_DefineColumn(SDDS_table, "yFrequency", NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0))<0 ||
        (histogram->columnIndex[3][0]
         =SDDS_DefineColumn(SDDS_table, "yp", NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0))<0 ||
        (histogram->columnIndex[3][1]
         =SDDS_DefineColumn(SDDS_table, "ypFrequency", NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0))<0) {
      fprintf(stdout, "Unable to define SDDS columns y and yp for file %s (%s)\n",
              filename, caller);
      fflush(stdout);
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      exit(1);
    }
    columns += 4;
  }
  if (histogram->longitData) {
    if ((histogram->columnIndex[4][0]
         =SDDS_DefineColumn(SDDS_table, "t", NULL, "s", NULL, NULL, SDDS_DOUBLE, 0))<0 ||
        (histogram->columnIndex[4][1]
         =SDDS_DefineColumn(SDDS_table, "tFrequency", NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0))<0 ||
        (histogram->columnIndex[5][0]
         =SDDS_DefineColumn(SDDS_table, "delta", NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0))<0 ||
        (histogram->columnIndex[5][1]
         =SDDS_DefineColumn(SDDS_table, "deltaFrequency", NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0))<0 ||
        (histogram->columnIndex[6][0]
         =SDDS_DefineColumn(SDDS_table, "dt", NULL, "s", NULL, NULL, SDDS_DOUBLE, 0))<0 ||
        (histogram->columnIndex[6][1]
         =SDDS_DefineColumn(SDDS_table, "dtFrequency", NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0))<0 ) {
      fprintf(stdout, "Unable to define SDDS columns t, dt, and p for file %s (%s)\n",
              filename, caller);
      fflush(stdout);
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      exit(1);
    }
    columns += 6;
  }
  if (!columns)
    bomb("no output selected for histogram", NULL);
  if (!SDDS_DefineSimpleParameter(SDDS_table, "Pass", NULL, SDDS_LONG) ||
      !SDDS_DefineSimpleParameter(SDDS_table, "Particles", NULL, SDDS_LONG) ||
      !SDDS_DefineSimpleParameter(SDDS_table, "pCentral", "m$be$nc", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleParameter(SDDS_table, "PassLength", "m", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleParameter(SDDS_table, "PassCentralTime", "s", SDDS_DOUBLE)) {
    fprintf(stdout, "Unable define SDDS parameter for file %s (%s)\n", filename, caller);
    fflush(stdout);
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    exit(1);
  }
  if (!SDDS_WriteLayout(SDDS_table)) {
    fprintf(stdout, "Unable to write SDDS layout for file %s (%s)\n", filename, caller);
    fflush(stdout);
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    exit(1);
  }
}

void dump_watch_particles(WATCH *watch, long step, long pass, double **particle, long particles, 
                          double Po, double length)
{
  long i, row;
  double p, t0, t;

  log_entry("dump_watch_particles");
  if (!watch->initialized)
    bomb("uninitialized watch-point (coordinate mode) encountered (dump_watch_particles)", NULL);
  if (!particle)
    bomb("NULL coordinate pointer passed to dump_watch_particles", NULL);
  if (watch->fraction>1)
    bomb("logic error--fraction>1 in dump_watch_particles", NULL);
  for (i=0; i<particles; i++)
    if (!particle[i]) {
      fprintf(stdout, "error: coordinate slot %ld is NULL (dump_watch_particles)\n", i);
      fflush(stdout);
      abort();
    }

  if (!SDDS_StartTable(&watch->SDDS_table, particles)) {
    SDDS_SetError("Problem starting SDDS table (dump_watch_particles)");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  row = 0;
  t0 = pass*length*sqrt(Po*Po+1)/(c_mks*(Po+1e-32));
  for (i=0; i<particles; i++) {
    if (watch->fraction==1 || random_2(0)<watch->fraction) {
      if (watch->xData && 
          !SDDS_SetRowValues(&watch->SDDS_table, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, row,
                             watch->xIndex[0], particle[i][0], watch->xIndex[1], particle[i][1],
                             -1)) {
        SDDS_SetError("Problem setting SDDS row values (dump_watch_particles)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      }
      if (watch->yData &&
          !SDDS_SetRowValues(&watch->SDDS_table, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, row,
                             watch->yIndex[0], particle[i][2], watch->yIndex[1], particle[i][3],
                             -1)) {
        SDDS_SetError("Problem setting SDDS row values (dump_watch_particles)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      }
      if (watch->longitData) {
        p = Po*(1+particle[i][5]);
        t = particle[i][4]/(c_mks*p/sqrt(sqr(p)+1));
        if (!SDDS_SetRowValues(&watch->SDDS_table, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, row,
                               watch->longitIndex[0], t,
                               watch->longitIndex[1], p,
                               watch->longitIndex[2], t-t0,
                               -1)) {
          SDDS_SetError("Problem setting SDDS row values (dump_watch_particles)");
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
      }
      if (!SDDS_SetRowValues(&watch->SDDS_table, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, row++,
                             watch->IDIndex, (long)particle[i][6], -1)) {
        SDDS_SetError("Problem setting SDDS row values (dump_watch_particles)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      }
    }
  }
  if (!SDDS_SetParameters(&watch->SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, 
                          "Step", step, "Pass", pass, "Particles", row, "pCentral", Po,
                          "PassLength", length, 
                          "PassCentralTime", t0, 
                          NULL)) {
    SDDS_SetError("Problem setting SDDS parameters (dump_watch_particles)");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  if (!SDDS_WriteTable(&watch->SDDS_table)) {
    SDDS_SetError("Problem writing SDDS table (dump_watch_particles)");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  } 
  log_exit("dump_watch_particles");
}

double tmp_safe_sqrt;
#define SAFE_SQRT(x) ((tmp_safe_sqrt=(x))<0?(double)0.0:sqrt(tmp_safe_sqrt))

void dump_watch_parameters(WATCH *watch, long step, long pass, long n_passes, double **particle, long particles, 
                           long original_particles,  double Po)
{
    long sample, i, j;
    double tc, last_tc, p_sum, gamma_sum, sum, p;

    log_entry("dump_watch_parameters");

    if (!watch->initialized)
        bomb("uninitialized watch-point (coordinate mode) encountered (dump_watch_parameters)", NULL);
    if (!particle)
        bomb("NULL coordinate pointer passed to dump_watch_parameters", NULL);
    if (watch->fraction>1)
        bomb("logic error--fraction>1 in dump_watch_parameters", NULL);

    for (i=0; i<particles; i++)
        if (!particle[i]) {
            fprintf(stdout, "error: coordinate slot %ld is NULL (dump_watch_parameters)\n", i);
            fflush(stdout);
            abort();
            }
    if (pass==0 && !SDDS_StartTable(&watch->SDDS_table, n_passes/watch->interval+1)) {
        SDDS_SetError("Problem starting SDDS table (dump_watch_parameters)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }

    sample = pass/watch->interval;
    if (particles) {
        /* compute centroids, sigmas, and emittances for x, y, and s */
        for (i=0; i<6; i+=2 ) {
            double x, xp, cx, cxp, sum_x, sum_xp, sum_x2, sum_xxp, sum_xp2;
            for (j=sum_x=sum_xp=0; j<particles; j++) {
                sum_x   += particle[j][i];
                sum_xp  += particle[j][i+1];
                }
            cx  = sum_x/particles;
            cxp = sum_xp/particles;
            if (!SDDS_SetRowValues(&watch->SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, sample,
                                   i==0?"Cx" :(i==2?"Cy" :"Cs" ), cx,
                                   i==0?"Cxp":(i==2?"Cyp":"Cdelta"), cxp, NULL)) {
                SDDS_SetError("Problem setting row values for SDDS table (dump_watch_parameters)");
                SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
                }
            if (watch->mode_code==WATCH_PARAMETERS) {
                for (j=sum_x2=sum_xp2=sum_xxp=0; j<particles; j++) {
                    x  = particle[j][i]   - cx;
                    xp = particle[j][i+1] - cxp;
                    sum_x2  += sqr(x);
                    sum_xxp += x*xp;
                    sum_xp2 += sqr(xp);
                    }
                if (!SDDS_SetRowValues(&watch->SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, sample,
                                       i==0?"Sx" :(i==2?"Sy" :"Ss" ), sqrt(sum_x2/particles),
                                       i==0?"Sxp":(i==2?"Syp":"Sdelta"), sqrt(sum_xp2/particles), NULL)) {
                    SDDS_SetError("Problem setting row values for SDDS table (dump_watch_parameters)");
                    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
                    }
                if (i!=4)
                    if (!SDDS_SetRowValues(&watch->SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, sample,
                                           i==0?"ex":"ey", SAFE_SQRT(sum_x2*sum_xp2-sqr(sum_xxp))/particles, NULL)) {
                        SDDS_SetError("Problem setting row values for SDDS table (dump_watch_parameters)");
                        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
                        }
                }
            }
        if (watch->mode_code==WATCH_PARAMETERS)
            if (!SDDS_SetRowValues(&watch->SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, sample,
                                   "el", rms_longitudinal_emittance(particle, particles, Po), NULL)) {
                SDDS_SetError("Problem setting row values for SDDS table (dump_watch_parameters)");
                SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
                }

        /* time centroid and sigma */
        for (i=p_sum=gamma_sum=sum=0; i<particles; i++) {
            p = Po*(1+particle[i][5]);
            p_sum     += p;
            gamma_sum += sqrt(sqr(p)+1);
            sum += particle[i][4]/(p/sqrt(sqr(p)+1)*c_mks) ;
            }
        tc = sum/particles;
        last_tc = 0;
        if (sample && !SDDS_GetValue(&watch->SDDS_table, "Ct", sample-1, &last_tc)) {
            SDDS_SetError("Problem retrieving row value from SDDS table (dump_watch_parameters)");
            SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
            }
        if (!SDDS_SetRowValues(&watch->SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, sample,
                               "Ct", tc, "dCt", tc-last_tc,
                               "pAverage", p_sum/particles, "pCentral", Po, 
                               "KAverage", (gamma_sum/particles-1)*me_mev, NULL)) {
            SDDS_SetError("Problem setting row values for SDDS table (dump_watch_parameters)");
            SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
            }

        if (watch->mode_code==WATCH_PARAMETERS) {
            for (i=sum=0; i<particles; i++)
                sum += sqr(particle[i][4]/(p/sqrt(sqr(p)+1)*c_mks) - tc);
            if (!SDDS_SetRowValues(&watch->SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, sample,
                               "St", sqrt(sum/particles), NULL)) {
                SDDS_SetError("Problem setting row values for SDDS table (dump_watch_parameters)");
                SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
                }
            }
        }

    /* number of particles */
    if (!SDDS_SetRowValues(&watch->SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, sample,
                           "Particles", particles, 
                           "Transmission", (original_particles?((double)particles)/original_particles:(double)0.0),
                           "Pass", sample*watch->interval, 
                           "Step", step, NULL)) {
        SDDS_SetError("Problem setting row values for SDDS table (dump_watch_parameters)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
    if (!SDDS_SetParameters(&watch->SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, 
                            "Step", step, NULL)) {
        SDDS_SetError("Problem setting parameter values for SDDS table (dump_watch_parameters)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }

    if (sample==(n_passes-1)/watch->interval) {
        if (!SDDS_WriteTable(&watch->SDDS_table)) {
            SDDS_SetError("Problem writing data for SDDS table (dump_watch_parameters)");
            SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
            }
        }

    log_exit("dump_watch_parameters");
    }


void dump_watch_FFT(WATCH *watch, long step, long pass, long n_passes, double **particle, long particles,
                           long original_particles,  double Po)
{
    long sample_index, i, samples;
    double sum_x, sum_y, sum_dp, *sample[4];

    log_entry("dump_watch_FFT");

    samples = n_passes/watch->interval;
    if (!power_of_2(samples)) 
        bomb("number of samples must be a power of 2 for FFT (dump_watch_FFT)", NULL);
    if (!watch->initialized)
        bomb("uninitialized watch-point (coordinate mode) encountered (dump_watch_FFT)", NULL);
    if (!particle)
        bomb("NULL coordinate pointer passed to dump_watch_FFT", NULL);
    for (i=0; i<particles; i++)
        if (!particle[i]) {
            fprintf(stdout, "error: coordinate slot %ld is NULL (dump_watch_FFT)\n", i);
            fflush(stdout);
            abort();
            }

    if (pass==0 && !SDDS_StartTable(&watch->SDDS_table, 2*samples)) {
        SDDS_SetError("Problem starting SDDS table (dump_watch_FFT)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }

    sample_index = pass/watch->interval;
    for (i=sum_x=sum_y=sum_dp=0; i<particles; i++) {
        sum_x  += particle[i][0];
        sum_y  += particle[i][2];
        sum_dp += particle[i][5];
        }
    if (particles) {
        sum_x  /= particles;
        sum_y  /= particles;
        sum_dp /= particles;
        }
    /* SDDS is used to store the centroid vs pass data prior to FFT */
    if (!SDDS_SetRowValues(&watch->SDDS_table, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, sample_index,
                           0, ((double)sample_index)/samples, 1, sum_x, 2, sum_y, 3, sum_dp, -1)) {
        SDDS_SetError("Problem setting row values for SDDS table (dump_watch_FFT)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }

    if (sample_index==samples-1) {
        /* sample/samples stored temporarily in frequency column: */
        if (!(sample[0]=(double*)SDDS_GetColumn(&watch->SDDS_table, "f"))  || 
            !(sample[1]=(double*)SDDS_GetColumn(&watch->SDDS_table, "Cx")) ||
            !(sample[2]=(double*)SDDS_GetColumn(&watch->SDDS_table, "Cy")) ||
            !(sample[3]=(double*)SDDS_GetColumn(&watch->SDDS_table, "Cdelta"))) {
            SDDS_SetError("Problem retrieving sample values from SDDS table (dump_watch_FFT)");
            SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
            }
        do_watch_FFT(sample, samples, 1, watch->window_code);
        do_watch_FFT(sample, samples, 2, watch->window_code);
        do_watch_FFT(sample, samples, 3, watch->window_code);
        if (!SDDS_EraseData(&watch->SDDS_table) ||
            !SDDS_SetColumn(&watch->SDDS_table, SDDS_SET_BY_INDEX, sample[0], samples/2, 0) ||
            !SDDS_SetColumn(&watch->SDDS_table, SDDS_SET_BY_INDEX, sample[1], samples/2, 1) ||
            !SDDS_SetColumn(&watch->SDDS_table, SDDS_SET_BY_INDEX, sample[2], samples/2, 2) ||
            !SDDS_SetColumn(&watch->SDDS_table, SDDS_SET_BY_INDEX, sample[3], samples/2, 3) ) {
            SDDS_SetError("Problem setting columns of SDDS table (dump_watch_FFT)");
            SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
            }
        for (i=0; i<4; i++)
            free(sample[i]);
        if (!SDDS_SetParameters(&watch->SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, 
                                "Step", step, NULL)) {
            SDDS_SetError("Problem setting parameter values for SDDS table (dump_watch_FFT)");
            SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
            }
        if (!SDDS_WriteTable(&watch->SDDS_table)) {
            SDDS_SetError("Problem writing data to SDDS file (dump_watch_FFT)");
            SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
            }
        }
    log_exit("dump_watch_FFT");
    }
 
void do_watch_FFT(double **data, long n_data, long slot, long window_code) 
{
    double *real_imag, r1, r2;
    long i;

    log_entry("do_watch_FFT");
    if (!data[slot])
        bomb("data pointer for specified slot is NULL (do_watch_FFT)", NULL);

    /* load into temporary array */
    real_imag = tmalloc(sizeof(*real_imag)*(n_data+2));
    for (i=0; i<n_data; i++)
        real_imag[i] = data[slot][i];

    /* apply the window */
    switch (window_code) {
      case FFT_HANNING:
        r1 = PIx2/(n_data-1);
        for (i=0; i<n_data; i++) 
            real_imag[i] *= (1-cos(i*r1))/2;
        break;
      case FFT_PARZEN:
        r1 = (n_data-1)/2.0;
        for (i=0; i<n_data; i++) 
            real_imag[i] *= 1-FABS((i-r1)/r1);
        break;
      case FFT_WELCH:
        r1 = (n_data-1)/2.0;
        r2 = sqr((n_data+1)/2.0);
        for (i=0; i<n_data; i++) 
            real_imag[i] *= 1 - sqr(i-r1)/r2;
        break;
      case FFT_UNIFORM:
        break;
      default:
        break;
        }

    /* do the FFT */
    if (!realFFT2(real_imag, real_imag, n_data, 0)) {
        fprintf(stdout, "error: FFT fails for slot %ld of watch point\n", slot);
        fflush(stdout);
        abort();
        }

    /* put the amplitudes back into the data array.  Multiply the non-DC
     * terms by 2 to account for amplitude in negative frequencies. 
     */
    data[slot][0] = fabs(real_imag[0]);
    for (i=1; i<n_data/2; i++) 
        data[slot][i] = sqrt(sqr(real_imag[2*i]) + (i==0?0:sqr(real_imag[2*i+1])))*2;
    
    free(real_imag); 
    log_exit("do_watch_FFT");
    }


void dump_particle_histogram(HISTOGRAM *histogram, long step, long pass, double **particle, long particles, 
                          double Po, double length)
{
  long icoord, ipart, ibin;
  double p, t0;
  static long maxBins = 0, maxParticles = 0;
  static double *frequency=NULL, *coordinate=NULL, *histData=NULL;
  double center, range, lower, upper;
  
  log_entry("dump_particle_histogram");
  if (!histogram->initialized)
    bomb("uninitialized histogram encountered (dump_particle_histogram)", NULL);
  if (!particle)
    bomb("NULL coordinate pointer passed to dump_particle_histogram", NULL);
  for (ipart=0; ipart<particles; ipart++)
    if (!particle[ipart]) {
      fprintf(stdout, "error: coordinate slot %ld is NULL (dump_particle_histogram)\n", ipart);
      fflush(stdout);
      abort();
    }

  if (!SDDS_StartTable(&histogram->SDDS_table, histogram->bins)) {
    SDDS_SetError("Problem starting SDDS table (dump_particle_histogram)");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  if (histogram->bins>maxBins) {
    maxBins = histogram->bins;
    if (!(frequency=SDDS_Realloc(frequency, sizeof(*frequency)*maxBins)) ||
        !(coordinate=SDDS_Realloc(coordinate, sizeof(*coordinate)*maxBins))) {
      SDDS_Bomb("Memory allocation failure (dump_particle_histogram)");
    }
  }
  if (particles>maxParticles) {
    maxParticles = particles;
    if (!(histData = SDDS_Realloc(histData, sizeof(*histData)*maxParticles)))
      SDDS_Bomb("Memory allocation failure (dump_particle_histogram)");
  }
  
  t0 = pass*length*sqrt(Po*Po+1)/(c_mks*(Po+1e-32));
  for (icoord=0; icoord<7; icoord++) {
    upper = -(lower=DBL_MAX);
    if (histogram->columnIndex[icoord][0]==-1)
      continue;
    if (icoord==4) {
      /* t */
      for (ipart=0; ipart<particles; ipart++) {
        p = Po*(1+particle[ipart][5]);
        histData[ipart] = particle[ipart][4]/(c_mks*p/sqrt(sqr(p)+1));
      }
    } else if (icoord==6) {
      /* dt (could reuse computations for t, but would need to store another array) */
      for (ipart=0; ipart<particles; ipart++) {
        p = Po*(1+particle[ipart][5]);
        histData[ipart] = particle[ipart][4]/(c_mks*p/sqrt(sqr(p)+1)) - t0;
      }
    } else {
      for (ipart=0; ipart<particles; ipart++)
        histData[ipart] = particle[ipart][icoord];
    }
    find_min_max(&lower, &upper, histData, particles);
    if (lower==upper) {
      lower -= 1e-16;
      upper += 1e-16;
    }
    if (histogram->count==0 || !histogram->fixedBinSize) {
      range = histogram->binSizeFactor*(upper-lower);
      histogram->binSize[icoord] = range/histogram->bins;
    } else
      range = histogram->binSize[icoord]*histogram->bins;
    center = (lower+upper)/2;
    lower = center-range/2;
    upper = center+range/2;
    for (ibin=0; ibin<histogram->bins; ibin++) 
      coordinate[ibin] = (ibin+0.5)*histogram->binSize[icoord] + lower;
    make_histogram(frequency, histogram->bins, lower, upper, histData, particles, 1);
    if (!SDDS_SetColumn(&histogram->SDDS_table, SDDS_SET_BY_INDEX, 
                        coordinate, histogram->bins, histogram->columnIndex[icoord][0]) ||
        !SDDS_SetColumn(&histogram->SDDS_table, SDDS_SET_BY_INDEX,
                        frequency, histogram->bins, histogram->columnIndex[icoord][1])) {
      SDDS_Bomb("Problem setting column values (dump_particle_histogram)");
    }
  }
  if (!SDDS_SetParameters(&histogram->SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, 
                          "Step", step, "Pass", pass, "Particles", particles, "pCentral", Po,
                          "PassLength", length, 
                          "PassCentralTime", t0, 
                          NULL)) {
    SDDS_SetError("Problem setting SDDS parameters (dump_particle_histogram)");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  if (!SDDS_WriteTable(&histogram->SDDS_table)) {
    SDDS_SetError("Problem writing SDDS table (dump_particle_histogram)");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  } 
  histogram->count++;
}


void dump_phase_space(SDDS_TABLE *SDDS_table, double **particle, long particles, long step, double Po)
{
    long i;
    double p;

    log_entry("dump_phase_space");
    if (!particle)
        bomb("NULL coordinate pointer passed to dump_phase_space", NULL);
    for (i=0; i<particles; i++)
        if (!particle[i]) {
            fprintf(stdout, "error: coordinate slot %ld is NULL (dump_phase_space)\n", i);
            fflush(stdout);
            abort();
            }

    if (!SDDS_StartTable(SDDS_table, particles)) {
        SDDS_SetError("Problem starting SDDS table (dump_phase_space)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
    if (!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, 
                            "Step", step, "pCentral", Po, NULL)) {
        SDDS_SetError("Problem setting parameter values for SDDS table (dump_phase_space)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
    for (i=0; i<particles; i++) {
        p = Po*(1+particle[i][5]);
        if (!SDDS_SetRowValues(SDDS_table, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, i,
                               0, particle[i][0], 1, particle[i][1], 2, particle[i][2], 3, particle[i][3],
                               4, particle[i][4]/(c_mks*p/sqrt(sqr(p)+1)), 5, p, 
                               6, (long)particle[i][6], -1)) {
            SDDS_SetError("Problem setting SDDS row values (dump_phase_space)");
            SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
            }
        }
    if (!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "Step", step, NULL)) {
        SDDS_SetError("Problem setting SDDS parameters (dump_phase_space)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
    if (!SDDS_WriteTable(SDDS_table)) {
        SDDS_SetError("Problem writing SDDS table (dump_phase_space)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        } 

    log_exit("dump_phase_space");
    }

void dump_lost_particles(SDDS_TABLE *SDDS_table, double **particle, long particles, long step)
{
    long i;

    log_entry("dump_lost_particles");
    if (!particle)
        bomb("NULL coordinate pointer passed to dump_lost_particles", NULL);
    for (i=0; i<particles; i++)
        if (!particle[i]) {
            fprintf(stdout, "error: coordinate slot %ld is NULL (dump_lost_particles)\n", i);
            fflush(stdout);
            abort();
            }

    if (!SDDS_StartTable(SDDS_table, particles)) {
        SDDS_SetError("Problem starting SDDS table (dump_lost_particles)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
    for (i=0; i<particles; i++) {
        if (!SDDS_SetRowValues(SDDS_table, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, i,
                               0, particle[i][0], 1, particle[i][1], 2, particle[i][2], 3, particle[i][3],
                               4, particle[i][4], 5, particle[i][5],
                               6, (long)particle[i][6], -1)) {
            SDDS_SetError("Problem setting SDDS row values (dump_lost_particles)");
            SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
            }
        }
    if (!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "Step", step, NULL)) {
        SDDS_SetError("Problem setting SDDS parameters (dump_lost_particles)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
    if (!SDDS_WriteTable(SDDS_table)) {
        SDDS_SetError("Problem writing SDDS table (dump_lost_particles)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        } 

    log_exit("dump_lost_particles");
    }

void dump_centroid(SDDS_TABLE *SDDS_table, BEAM_SUMS *sums, LINE_LIST *beamline, long n_elements, long step,
                          double p_central)
{
    long i, j;
    BEAM_SUMS *beam;
    ELEMENT_LIST *eptr;
    double *cent;
    char *name, *type_name;
    long s_index, Cx_index, occurence;

    log_entry("dump_centroid");

    if (!SDDS_table)
        bomb("SDDS_TABLE pointer is NULL (dump_centroid)", NULL);
    if (!sums)
        bomb("BEAM_SUMS pointer is NULL (dump_centroid)", NULL);
    if (!beamline)
        bomb("LINE_LIST pointer is NULL (dump_centroid)", NULL);
    if ((s_index=SDDS_GetColumnIndex(SDDS_table, "s"))<0 ||
        (Cx_index=SDDS_GetColumnIndex(SDDS_table, "Cx"))<0) {
        SDDS_SetError("Problem getting index of SDDS columns (dump_centroid)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }

    if (!SDDS_StartTable(SDDS_table, n_elements+1)) {
        SDDS_SetError("Problem starting SDDS table (dump_centroid)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
    if (!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, 
                            "Step", step, NULL)) {
        SDDS_SetError("Problem setting parameter values for SDDS table (dump_centroid)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
 
    if (!(beamline->flags&BEAMLINE_CONCAT_DONE))
        eptr = &(beamline->elem);
    else 
        eptr = &(beamline->ecat);
    name = "_BEG_";
    type_name = "MARK";
    occurence = 1;
    cent = tmalloc(sizeof(*cent)*6);
    for (i=0; i<n_elements; i++) {
        if (!eptr) {
            fprintf(stdout, "element pointer is NULL, i=%ld (dump_centroid)", i);
            fflush(stdout);
            abort();
            }
        beam = sums+i;
        if (!beam->sum) {
            fprintf(stdout, "beam->sum is NULL, i=%ld (dump_centroid)\n", i);
            fflush(stdout);
            abort();
            }
        if (beam->n_part)
            for (j=0; j<6; j++) {
                cent[j] = beam->sum[j]/beam->n_part;
                if (isnan(cent[j]) || isinf(cent[j]))
                    cent[j] = DBL_MAX;
                }
        else
            for (j=0; j<6; j++)
                cent[j] = 0;
        if (!SDDS_SetRowValues(SDDS_table, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, i,
                               Cx_index  , cent[0], Cx_index+1, cent[1], Cx_index+2, cent[2],
                               Cx_index+3, cent[3], Cx_index+4, cent[4], Cx_index+5, cent[5],
                               Cx_index+6, beam->n_part, Cx_index+7, beam->p0, 
                               s_index, beam->z, s_index+1, name, s_index+2, occurence,
                               s_index+3, type_name, -1)) {
            SDDS_SetError("Problem setting row values for SDDS table (dump_centroid)");
            SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
            }
        
        name = eptr->name;
        type_name = entity_name[eptr->type];
        occurence = eptr->occurence;
        if (eptr->succ)
            eptr = eptr->succ;
        else if (beamline->elem_recirc)
            eptr = beamline->elem_recirc;
        else
            eptr = &(beamline->elem);
        }

    if (!SDDS_WriteTable(SDDS_table)) {
        SDDS_SetError("Unable to write centroid data (dump_centroid)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
    if (!SDDS_EraseData(SDDS_table)) {
        SDDS_SetError("Unable to erase centroid data (dump_centroid)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
    free(cent);
    log_exit("dump_centroid");
    }

void dump_sigma(SDDS_TABLE *SDDS_table, BEAM_SUMS *sums, LINE_LIST *beamline, long n_elements, long step,
                double p_central)
{
  long i, j, ie, offset, plane, index;
  BEAM_SUMS *beam;
  ELEMENT_LIST *eptr;
  double Sigma[6], sigma[6][6], centroid[6], emit, emitNorm;
  char *name, *type_name;
  long s_index, ma1_index, Sx_index, occurence, ex_index;
  long sNIndex[6];

  if (!SDDS_table)
    bomb("SDDS_TABLE pointer is NULL (dump_sigma)", NULL);
  if (!sums)
    bomb("BEAM_SUMS pointer is NULL (dump_centroid)", NULL);
  if (!beamline)
    bomb("LINE_LIST pointer is NULL (dump_centroid)", NULL);
  if ((s_index=SDDS_GetColumnIndex(SDDS_table, "s"))<0 ||
      (sNIndex[0]=SDDS_GetColumnIndex(SDDS_table, "s1"))<0 ||
      (sNIndex[1]=SDDS_GetColumnIndex(SDDS_table, "s2"))<0 ||
      (sNIndex[2]=SDDS_GetColumnIndex(SDDS_table, "s3"))<0 ||
      (sNIndex[3]=SDDS_GetColumnIndex(SDDS_table, "s4"))<0 ||
      (sNIndex[4]=SDDS_GetColumnIndex(SDDS_table, "s5"))<0 ||
      (sNIndex[5]=SDDS_GetColumnIndex(SDDS_table, "s6"))<0 ||
      (ma1_index=SDDS_GetColumnIndex(SDDS_table, "ma1"))<0 ||
      (Sx_index=SDDS_GetColumnIndex(SDDS_table, "Sx"))<0 ||
      (ex_index=SDDS_GetColumnIndex(SDDS_table, "ex"))<0) {
    SDDS_SetError("Problem getting index of SDDS columns (dump_sigma)");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }

  if (!SDDS_StartTable(SDDS_table, n_elements+1)) {
    SDDS_SetError("Problem starting SDDS table (dump_sigma)");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  if (!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, 
                          "Step", step, NULL)) {
    SDDS_SetError("Problem setting parameter values for SDDS table (dump_sigma)");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }

  if (!(beamline->flags&BEAMLINE_CONCAT_DONE))
    eptr = &(beamline->elem);
  else 
    eptr = &(beamline->ecat);
  name = "_BEG_";
  type_name = "MARK";
  occurence = 1;
  for (ie=0; ie<n_elements; ie++) {
    beam  = sums+ie;
    if (!beam->sum || !beam->sum2) {
      fprintf(stdout, "error: sum or sum2 element of BEAM_SUMS is undefined for element %ld (%s)\n",
              ie, name);
      fflush(stdout);
      exit(1);
    }
    if (beam->n_part) {
      /* compute centroids, sigma matrix, and beam sizes (with centroid term removed) */
      for (i=0; i<6; i++)
        centroid[i] = beam->sum[i]/beam->n_part;
      for (i=0; i<6; i++) {
        sigma[i][i] = beam->sum2[i][i]/beam->n_part;
        Sigma[i] = SAFE_SQRT(beam->sum2[i][i]/beam->n_part-sqr(centroid[i]));
#if DO_NORMEMIT_SUMS
        if (i<4) {
          pcentroid[i] = beam->psum[i]/beam->n_part;
          psigma[i][i] = beam->psum2[i][i]/beam->n_part;
          pSigma[i] = SAFE_SQRT(beam->psum2[i][i]/beam->n_part-sqr(pcentroid[i])); 
        }
#endif
        if (!SDDS_SetRowValues(SDDS_table, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, ie, 
                               sNIndex[i], sqrt(sigma[i][i]), 
                               Sx_index+i, Sigma[i],
                               -1)) {
          SDDS_SetError("Problem setting SDDS row values (dump_sigma)");
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
        offset = 1;
        for (j=i+1; j<6; j++, offset++) {
          sigma[i][j] = beam->sum2[i][j]/beam->n_part - centroid[i]*centroid[j]; 
#if DO_NORMEMIT_SUMS
          if (i<4 && j<4)
            psigma[i][j] = beam->psum2[i][j]/beam->n_part; 
#endif
          if (!SDDS_SetRowValues(SDDS_table, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, ie, 
                                 sNIndex[i]+offset, sigma[i][j], -1)) {
            SDDS_SetError("Problem setting SDDS row values (dump_sigma)");
            SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
          }
        }
      }
      
      /* Set values for maximum amplitudes of transverse coordinates */
      for (i=0; i<4; i++)
        if (!SDDS_SetRowValues(SDDS_table, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, ie, ma1_index+i, beam->maxabs[i], -1)) {
          SDDS_SetError("Problem setting SDDS row values (dump_sigma)");
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
      for (plane=0; plane<=2; plane+=2) {
        /* emittance */
        emit = SAFE_SQRT(sqr(Sigma[0+plane]*Sigma[1+plane]) 
                         - sqr(sigma[0+plane][1+plane]));
#if DO_NORMEMIT_SUMS
        emitNorm = SAFE_SQRT(sqr(pSigma[0+plane]*pSigma[1+plane]) 
                         - sqr(psigma[0+plane][1+plane]-pcentroid[0+plane]*pcentroid[1+plane]));
#else
        emitNorm = emit*beam->p0;
#endif
        if (!SDDS_SetRowValues(SDDS_table, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, ie, 
                               ex_index+plane, emit, 
                               ex_index+1+plane, emitNorm,
                               -1)) {
          SDDS_SetError("Problem setting SDDS row values (dump_sigma)");
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
      }
    }
    else {
      for (index=sNIndex[0]; index<sNIndex[0]+21; index++) 
        if (!SDDS_SetRowValues(SDDS_table, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, ie, index, (double)0.0, -1)) {
          SDDS_SetError("Problem setting SDDS row values (dump_sigma)");
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
      for (index=ma1_index; index<ma1_index+4; index++) 
        if (!SDDS_SetRowValues(SDDS_table, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, ie, index, (double)0.0, -1)) {
          SDDS_SetError("Problem setting SDDS row values (dump_sigma)");
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
      for (index=Sx_index; index<Sx_index+6; index++) 
        if (!SDDS_SetRowValues(SDDS_table, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, ie, index, (double)0.0, -1)) {
          SDDS_SetError("Problem setting SDDS row values (dump_sigma)");
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
      for (index=ex_index; index<ex_index+4; index++)
        if (!SDDS_SetRowValues(SDDS_table, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, ie, index, (double)0.0, -1)) {
          SDDS_SetError("Problem setting SDDS row values (dump_sigma)");
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
    }
    
    if (!SDDS_SetRowValues(SDDS_table,  SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, ie, 
                           s_index, beam->z, s_index+1, name, s_index+2, occurence,
                           s_index+3, type_name, -1)) {
      SDDS_SetError("Problem setting row values for SDDS table (dump_sigma)");
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
    
    name = eptr->name;
    type_name = entity_name[eptr->type];
    occurence = eptr->occurence;
    if (eptr->succ)
      eptr = eptr->succ;
    else if (beamline->elem_recirc)
      eptr = beamline->elem_recirc;
    else
      eptr = &(beamline->elem);
  }

  if (!SDDS_WriteTable(SDDS_table)) {
    SDDS_SetError("Unable to write sigma data (dump_sigma)");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  if (!SDDS_EraseData(SDDS_table)) {
    SDDS_SetError("Unable to erase sigma data (dump_sigma)");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
}

