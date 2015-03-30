/*************************************************************************\
* Copyright (c) 2015 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2015 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* program: sddsbunchingfactor
 * purpose: compute bunching factor from time data in elegant output file.
 *
 * Michael Borland, 2015
 *
 */
#include "mdb.h"
#include "scan.h"
#include "SDDS.h"

#define SET_PIPE 0
#define SET_OMEGARANGE 1
#define SET_POINTS 2
#define SET_MODE 3
#define N_OPTIONS 4

char *option[N_OPTIONS] = {
  "pipe", "omegarange", "points", "mode",
} ;

char *modeOption[2] = {
  "linear", "logarithmic"
};

char *USAGE="sddsbunchingfactor [-pipe=[input][,output]] [<SDDSinputfile>] [<SDDSoutputfile>]\n\
  [-omegaRange=<lower>,<upper>] [-points=<number>] [-mode={linear|logarithmic}]\n\
Computes the bunching factor vs angular frequency using time coordinates of particles in the input.\n\
-pipe         The standard SDDS pipe option.\n\
-omegaRange   Lower and upper limits of angular frequency in Hz.\n\
-points       Number of values of omega.\n\
-mode         Linear or logarithmic spacing of omega points?\n\n\
Program by Michael Borland.  (This is version 1, March 30, 2015)\n";

long SetUpOutputFile(SDDS_DATASET *SDDSout, char *outputfile, SDDS_DATASET *SDDSin);
long check_sdds_beam_column(SDDS_TABLE *SDDS_table, char *name, char *units);

int main(int argc, char **argv)
{
  SDDS_DATASET SDDSin, SDDSout;
  char *inputfile, *outputfile;
  long i_arg, readCode, omegaMode;
  long iw, ip, particles;
  SCANNED_ARG *s_arg;
  unsigned long pipeFlags;
  double *t, omega, cosSum, sinSum;
  double omegaLower, omegaUpper;
  long omegaPoints, tmpFileUsed;
  
  SDDS_RegisterProgramName(argv[0]);
  argc = scanargs(&s_arg, argc, argv);
  if (argc<2) 
    bomb(NULL, USAGE);

  inputfile = outputfile = NULL;
  pipeFlags = 0;
  omegaPoints = -1;
  omegaLower = omegaUpper = 0;
  omegaMode = 0;
  
  for (i_arg=1; i_arg<argc; i_arg++) {
    if (s_arg[i_arg].arg_type==OPTION) {
      switch (match_string(s_arg[i_arg].list[0], option, N_OPTIONS, 0)) {
      case SET_PIPE:
        if (!processPipeOption(s_arg[i_arg].list+1, s_arg[i_arg].n_items-1, &pipeFlags))
          SDDS_Bomb("invalid -pipe syntax");
        break;
      case SET_OMEGARANGE:
	if (s_arg[i_arg].n_items!=3 ||
	    sscanf(s_arg[i_arg].list[1], "%lf", &omegaLower)!=1 ||
	    sscanf(s_arg[i_arg].list[2], "%lf", &omegaUpper)!=1 ||
            omegaLower>=omegaUpper)
	  SDDS_Bomb("invalid -omegaRange syntax");
        break;
      case SET_POINTS:
	if (s_arg[i_arg].n_items!=2 ||
	    sscanf(s_arg[i_arg].list[1], "%ld", &omegaPoints)!=1 ||
            omegaPoints<=0)
	  SDDS_Bomb("invalid -omegaPoints syntax");
        break;
      case SET_MODE:
        if (s_arg[i_arg].n_items!=2 ||
            (omegaMode=match_string(s_arg[i_arg].list[1], modeOption, 2, 0))<0)
	  SDDS_Bomb("invalid -mode syntax");
        break;
      default:
        fprintf(stdout, "error: unknown switch: %s\n", s_arg[i_arg].list[0]);
        fflush(stdout);
        exit(1);
        break;
      }
    }
    else {
      if (inputfile==NULL)
        inputfile = s_arg[i_arg].list[0];
      else if (outputfile==NULL)
        outputfile = s_arg[i_arg].list[0];
      else
        SDDS_Bomb("too many filenames");
    }
  }

  processFilenames("sddsbunchingfactor", &inputfile, &outputfile, pipeFlags, 0, &tmpFileUsed);
  if (tmpFileUsed)
    SDDS_Bomb("can't overwrite input file");
  
  if (!SDDS_InitializeInput(&SDDSin, inputfile))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);

  /* check the input file for valid data */
  if (!check_sdds_beam_column(&SDDSin, "t", "s")) {
    fprintf(stderr, 
            "sddsbunchingfactor: column t is not present or has the wrong units (should be \"s\") in file %s\n", 
            inputfile);
    exit(1);
  }

  if (!SetUpOutputFile(&SDDSout, outputfile, &SDDSin)) 
    SDDS_Bomb("problem setting up output file");
  
  while ((readCode=SDDS_ReadPage(&SDDSin))>0) {
    if (!SDDS_StartPage(&SDDSout, omegaPoints) || !SDDS_CopyParameters(&SDDSout, &SDDSin))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    if ((particles=SDDS_RowCount(&SDDSin))>1) {
      if (!(t =  SDDS_GetColumnInDoubles(&SDDSin, "t")))
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      for (iw=0; iw<omegaPoints; iw++) {
        if (omegaMode==0) {
          /* linear */
          omega = omegaLower + iw*(omegaUpper-omegaLower)/(omegaPoints-1);
        } else {
          /* logarithmic */
          omega = omegaLower*pow(10, iw*log10(omegaUpper/omegaLower)/(omegaPoints-1.));
        }
        cosSum = sinSum = 0;
        for (ip=0; ip<particles; ip++) {
          cosSum += cos(omega*t[ip]);
          sinSum += sin(omega*t[ip]);
        }
        if (!SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, iw, 
                               "omega", omega,
                               "BunchingFactor", sqrt(cosSum*cosSum+sinSum*sinSum)/particles, 
                               NULL)) {
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
      }
      free(t);
      t = NULL;
    }
    if (!SDDS_SetParameters(&SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "Particles", particles, NULL) ||
        !SDDS_WritePage(&SDDSout)) {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
  }
  if (readCode==0)
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  if (!SDDS_Terminate(&SDDSin) || !SDDS_Terminate(&SDDSout))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  return 0;
}

long SetUpOutputFile(SDDS_DATASET *SDDSout, char *outputfile, SDDS_DATASET *SDDSin)
{
  if (!SDDS_InitializeOutput(SDDSout, SDDS_BINARY, 1, NULL, NULL, outputfile) ||
      SDDS_DefineColumn(SDDSout, "omega", "$gw$r", "1/s", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineColumn(SDDSout, "BunchingFactor", NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineParameter(SDDSout, "Particles", NULL, NULL, NULL, NULL, SDDS_LONG, 0)<0 ||
      !SDDS_TransferAllParameterDefinitions(SDDSout, SDDSin, SDDS_TRANSFER_KEEPOLD)) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  if (!SDDS_SaveLayout(SDDSout) || !SDDS_WriteLayout(SDDSout)) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  return(1);
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

