/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* program: sddsanalyzebeam.c
 * purpose: take particle coordinate input and analyze to give
 *          moments and twiss parameters
 *
 * Michael Borland, 2000
 *
 $Log: not supported by cvs2svn $
 Revision 1.6  2001/10/15 20:37:06  soliday
 Cleaned up for Linux.

 Revision 1.5  2001/10/15 15:42:33  soliday
 Cleaned up for WIN32.

 Revision 1.4  2001/08/06 13:56:40  borland
 Fixed memory leak.

 Revision 1.3  2001/07/02 18:09:46  borland
 Usage and error messages now mention the requirement for 't' column.

 Revision 1.2  2000/12/06 02:14:47  borland
 Now provides "corrected" and "uncorrected" values for beta and alpha.
 The corrected values are the best ones if one really knows the dispersion.

 Revision 1.1  2000/08/14 21:44:57  borland
 First version.

 */
#include "mdb.h"
#include "scan.h"
#include "SDDS.h"

#define DEBUG 0

#define SET_PIPE 0
#define SET_NOWARNINGS 1
#define N_OPTIONS 2

char *option[N_OPTIONS] = {
  "pipe", "nowarnings"
} ;

char *USAGE="sddsanalyzebeam [-pipe=[input][,output]] [<SDDSinputfile>] [<SDDSoutputfile>]\n\
  [-nowarnings]\n\
Computes Twiss parameters and other properties of a particle beam.\n\
The input file must have columns x, xp, y, yp, t, and p; for example, an elegant\n\
beam output file is acceptable.\n\
Program by Michael Borland.  (This is version 1, August 2000.)\n";


long SetUpOutputFile(SDDS_DATASET *SDDSout, char *outputfile);
long check_sdds_beam_column(SDDS_TABLE *SDDS_table, char *name, char *units);

char *CenName[6];
char *CorName[6][6];
char *psName[6] = {"x", "xp", "y", "yp", "t", "delta" };
  
int main(int argc, char **argv)
{
  SDDS_DATASET SDDSin, SDDSout;
  char *inputfile, *outputfile;
  long iPart, particles, i_arg, readCode, noWarnings, row;
  long i, j;
  SCANNED_ARG *s_arg;
  unsigned long pipeFlags;
  double *x, *xp, *y, *yp, *p=NULL, *t;
  double pAve=0.0, sum;
  double S[6][6], C[6], beta[2], alpha[2], eta[4], emit[2], emitcor[2], beamsize[6];
  double betacor[2], alphacor[2];
  double *data[6], Sbeta[4][4];
  long tmpFileUsed;
  
  SDDS_RegisterProgramName(argv[0]);
  argc = scanargs(&s_arg, argc, argv);
  if (argc<2) 
    bomb(NULL, USAGE);

  inputfile = outputfile = NULL;
  pipeFlags = noWarnings = 0;
  
  for (i_arg=1; i_arg<argc; i_arg++) {
    if (s_arg[i_arg].arg_type==OPTION) {
      switch (match_string(s_arg[i_arg].list[0], option, N_OPTIONS, 0)) {
      case SET_PIPE:
        if (!processPipeOption(s_arg[i_arg].list+1, s_arg[i_arg].n_items-1, &pipeFlags))
          SDDS_Bomb("invalid -pipe syntax");
        break;
      case SET_NOWARNINGS:
        noWarnings = 1;
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

  processFilenames("sddsanalyzebeam", &inputfile, &outputfile, pipeFlags, noWarnings, &tmpFileUsed);
  if (tmpFileUsed)
    SDDS_Bomb("can't overwrite input file");
  
  if (!SDDS_InitializeInput(&SDDSin, inputfile))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);

  /* check the input file for valid data */
  if (!check_sdds_beam_column(&SDDSin, "x", "m") || !check_sdds_beam_column(&SDDSin, "y", "m") ||
      !check_sdds_beam_column(&SDDSin, "xp", NULL) || !check_sdds_beam_column(&SDDSin, "yp", NULL) ||
      !check_sdds_beam_column(&SDDSin, "t", "s") ||
      (!check_sdds_beam_column(&SDDSin, "p", "m$be$nc") && !check_sdds_beam_column(&SDDSin, "p", NULL))) {
    fprintf(stderr, 
            "sddsanalyzebeam: one or more of (x, xp, y, yp, t, p) have the wrong units or are not present in %s", 
            inputfile);
    exit(1);
  }

  if (!SetUpOutputFile(&SDDSout, outputfile))
    SDDS_Bomb("problem setting up output file");

  for (i=0; i<6; i++)
    data[i] = NULL;
  
  row = 0;
  while ((readCode=SDDS_ReadPage(&SDDSin))>0) {
    if (readCode!=1 && !SDDS_LengthenTable(&SDDSout, 1))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    for (i=0; i<2; i++)
      beta[i] = alpha[i] = emit[i] = emitcor[i] = eta[i] = eta[i+2] = 0;
    for (i=0; i<6; i++) {
      for (j=C[i]=0; j<=i; j++)
        S[i][j] = 0;
      if (data[i])
        free(data[i]);
    }
    if ((particles=SDDS_RowCount(&SDDSin))>2) {
      if (!(data[0] = x = SDDS_GetColumnInDoubles(&SDDSin, "x")) ||
          !(data[1] = xp = SDDS_GetColumnInDoubles(&SDDSin, "xp")) ||
          !(data[2] = y = SDDS_GetColumnInDoubles(&SDDSin, "y")) ||
          !(data[3] = yp = SDDS_GetColumnInDoubles(&SDDSin, "yp")) ||
          !(data[4] = t =  SDDS_GetColumnInDoubles(&SDDSin, "t")) ||
          !(data[5] = p = SDDS_GetColumnInDoubles(&SDDSin, "p")))
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      /* compute and subtract off average values */
      for (i=0; i<5; i++) {
        C[i] = arithmeticAverage(data[i], particles);
        for (j=0; j<particles; j++)
          data[i][j] -= C[i];
      }
      /* convert momentum to (p-po)/po */
      C[6] = 0;
      pAve = arithmeticAverage(p, particles);
      for (i=0; i<particles; i++)
        p[i] = p[i]/pAve-1;

      /* compute correlations */
      for (i=0; i<6; i++) {
        for (j=0; j<=i; j++) {
          for (iPart=sum=0; iPart<particles; iPart++)
            sum += data[i][iPart]*data[j][iPart];
          S[j][i] = S[i][j] = sum/particles;
        }
      }
      for (i=0; i<6; i++)
        beamsize[i] = sqrt(S[i][i]);

      /* compute correlations with energy correlations removed */
      if (S[5][5])
        for (i=0; i<4; i++) 
          eta[i] = S[i][5]/S[5][5];
      else 
        for (i=0; i<4; i++) 
          eta[i] = 0;
      for (i=0; i<4; i++) {
        for (iPart=0; iPart<particles; iPart++)
          data[i][iPart] -= eta[i]*data[5][iPart];
      }
      for (i=0; i<4; i++) {
        for (j=0; j<4; j++) {
          for (iPart=sum=0; iPart<particles; iPart++)
            sum += data[i][iPart]*data[j][iPart];
          Sbeta[j][i] = Sbeta[i][j] = sum/particles;
        }
      }
      /* compute beta functions etc */
      for (i=0; i<2; i++) {
        emitcor[i] = emit[i] = beta[i] = alpha[i] = 0;
        if ((emit[i] = S[2*i+0][2*i+0]*S[2*i+1][2*i+1]-sqr(S[2*i+0][2*i+1]))>0) {
          emit[i] = sqrt(emit[i]);
          beta[i] = S[2*i+0][2*i+0]/emit[i];
          alpha[i] = -S[2*i+0][2*i+1]/emit[i];
        }
        if ((emitcor[i] = Sbeta[2*i+0][2*i+0]*Sbeta[2*i+1][2*i+1]-sqr(Sbeta[2*i+0][2*i+1]))>0) {
          emitcor[i] = sqrt(emitcor[i]);
          betacor[i] = Sbeta[2*i+0][2*i+0]/emitcor[i];
          alphacor[i] = -Sbeta[2*i+0][2*i+1]/emitcor[i];
        }
      }
    }
    /* set centroids and sigmas */
    for (i=0; i<6; i++) {
      if (!SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,
                             row, CenName[i], C[i], NULL)) {
        SDDS_SetError("Problem setting centroid value");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      }
      for (j=0; j<=i; j++) {
        if (!SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,
                               row, CorName[i][j], S[i][j], NULL)) {
          SDDS_SetError("Problem setting sigma value");
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
      }
    }
    if (!SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, row,
                           "Sx", beamsize[0], "Sxp", beamsize[1],
                           "Sy", beamsize[2], "Syp", beamsize[3],
                           "St", beamsize[4], "Sdelta", beamsize[5], NULL)) {
      SDDS_SetError("Problem setting sigma value");
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
    
    /* set twiss parameter, etc */
    if (!SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, row, 
                           "betax", beta[0], "alphax", alpha[0], "etax", eta[0],
                           "etaxp", eta[1], "ex", emit[0], "ecx", emitcor[0],
                           "enx", emit[0]*pAve, "ecnx", emitcor[0]*pAve,
                           "betacx", betacor[0], "alphacx", alphacor[0],
                           "betay", beta[1], "alphay", alpha[1], "etay", eta[2],
                           "etayp", eta[3], "ey", emit[1], 
                           "ecy", emitcor[1], "eny", emit[1]*pAve, "ecny", emitcor[1]*pAve,
                           "betacy", betacor[1], "alphacy", alphacor[1],
                           "pAverage", pAve, NULL)) {
      SDDS_SetError("Problem setting Twiss values");
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
    row++;
  }
  if (readCode==0)
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  if (!SDDS_WritePage(&SDDSout) || !SDDS_Terminate(&SDDSin) ||
      !SDDS_Terminate(&SDDSout))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  return 0;
}


long SetUpOutputFile(SDDS_DATASET *SDDSout, char *outputfile)
{
  char units[128], name[128];
  char *ppUnits[6] = {"m", "", "m", "", "s", ""} ;
  long i, j;
  
  if (!SDDS_InitializeOutput(SDDSout, SDDS_BINARY, 1, NULL, NULL, outputfile) ||
      !SDDS_DefineSimpleColumn(SDDSout, "ex", "m", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(SDDSout, "enx", "m", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(SDDSout, "betax", "m", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(SDDSout, "alphax", NULL, SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(SDDSout, "ecx", "m", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(SDDSout, "ecnx", "m", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(SDDSout, "betacx", "m", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(SDDSout, "alphacx", NULL, SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(SDDSout, "etax", "m", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(SDDSout, "etaxp", NULL, SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(SDDSout, "ey", "m", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(SDDSout, "eny", "m", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(SDDSout, "betay", "m", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(SDDSout, "alphay", NULL, SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(SDDSout, "ecy", "m", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(SDDSout, "ecny", "m", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(SDDSout, "betacy", "m", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(SDDSout, "alphacy", NULL, SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(SDDSout, "etay", "m", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(SDDSout, "etayp", NULL, SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(SDDSout, "Sx", "m", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(SDDSout, "Sxp", NULL, SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(SDDSout, "Sy", "m", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(SDDSout, "Syp", NULL, SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(SDDSout, "St", "s", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(SDDSout, "Sdelta", NULL, SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(SDDSout, "pAverage", NULL, SDDS_DOUBLE)) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  for (i=0; i<6; i++) {
    sprintf(name, "C%s", psName[i]);
    if (!SDDS_CopyString(&CenName[i], name))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    if (!SDDS_DefineSimpleColumn(SDDSout, name, ppUnits[i], SDDS_DOUBLE))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    for (j=0; j<=i; j++) {
      sprintf(name, "s%ld%ld", j+1, i+1);
      units[0] = 0;
      if (ppUnits[i]) {
        if (ppUnits[j]) {
          if (strcmp(ppUnits[i], ppUnits[j])==0 && strlen(ppUnits[i]))
            sprintf(units, "%s$a2$n", ppUnits[i]);
          else if (strlen(ppUnits[i]))
            sprintf(units, "%s %s", ppUnits[i], ppUnits[j]);
        }
        else {
          strcpy(units, ppUnits[i]);
        }
      } else if (ppUnits[j])
        strcpy(units, ppUnits[j]);
      if (!SDDS_DefineSimpleColumn(SDDSout, name, units, SDDS_DOUBLE))
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      if (!SDDS_CopyString(&CorName[i][j], name))
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
  }
  if (!SDDS_SaveLayout(SDDSout) || !SDDS_WriteLayout(SDDSout) ||
      !SDDS_StartPage(SDDSout, 1)) {
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

