/*************************************************************************\
* Copyright (c) 2009 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2009 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* program: astra2elegant.c
 * purpose: Convert ASTRA phase space output to a form acceptable to elegant
 *
 * Michael Borland, 2009
 *
 $Log: not supported by cvs2svn $
 Revision 1.1  2009/07/28 14:15:51  borland
 First version.

 
 */
#include "mdb.h"
#include "scan.h"
#include "SDDS.h"

#define MAX_ROW_INCREMENT 100000

#define SET_PIPE 0
#define N_OPTIONS 1

char *option[N_OPTIONS] = {
  "pipe", 
} ;

char *USAGE="astra2elegant [-pipe=[input][,output]] [<ASTRAoutputfile>] [<SDDSoutputfile>]\n\
Converts ASTRA phase space output to a form acceptable to elegant.\n\
Program by Michael Borland.  (This is version 1, July 2009.)\n";

long SetUpOutputFile(SDDS_DATASET *SDDSout, char *outputfile);
long addParticleToFile(SDDS_DATASET *SDDSout, long *row, long *maxRows,
                       double x, double y, double z, double px, double py, double pz);

int main(int argc, char **argv)
{
  SDDS_DATASET SDDSout;
  char *inputfile, *outputfile;
  long row, maxRows;
  SCANNED_ARG *s_arg;
  unsigned long pipeFlags;
  FILE *fpin;
  double x, y, z, px, py, pzRef, dpz, clock, charge, totalCharge;
  long index, flag, i_arg, tmpFileUsed;
  
  SDDS_RegisterProgramName(argv[0]);
  argc = scanargs(&s_arg, argc, argv);
  if (argc<2) 
    bomb(NULL, USAGE);

  inputfile = outputfile = NULL;
  pipeFlags = 0;
  
  for (i_arg=1; i_arg<argc; i_arg++) {
    if (s_arg[i_arg].arg_type==OPTION) {
      switch (match_string(s_arg[i_arg].list[0], option, N_OPTIONS, 0)) {
      case SET_PIPE:
        if (!processPipeOption(s_arg[i_arg].list+1, s_arg[i_arg].n_items-1, &pipeFlags))
          SDDS_Bomb("invalid -pipe syntax");
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

  if (pipeFlags&USE_STDIN)
    SDDS_Bomb("give input filename");
  processFilenames("astra2elegant", &inputfile, &outputfile, pipeFlags, 0, &tmpFileUsed);
  if (tmpFileUsed)
    SDDS_Bomb("can't overwrite input file");
  
  if (!(fpin=fopen(inputfile, "r")))
    SDDS_Bomb("problem opening input file");
  if (fscanf(fpin, "%le %le %le %le %le %le %le %le %ld %ld",
             &x, &y, &z, &px, &py, &pzRef, &clock, &charge, &index, &flag)!=10 ||
      flag!=5)
    SDDS_Bomb("problem reading reference particle data");
  totalCharge = charge;
  
  if (!SetUpOutputFile(&SDDSout, outputfile) ||
      !SDDS_StartTable(&SDDSout, maxRows=MAX_ROW_INCREMENT))
    SDDS_Bomb("problem setting up output file");

  row = 0;
  if (!addParticleToFile(&SDDSout, &row, &maxRows, x, y, 0.0, px, py, pzRef))
    return 1;

  while (!feof(fpin) &&
         fscanf(fpin, "%le %le %le %le %le %le %le %le %ld %ld",
                &x, &y, &z, &px, &py, &dpz, &clock, &charge, &index, &flag)==10) {
    if (flag==3 || flag==5) {
      if (!addParticleToFile(&SDDSout, &row, &maxRows, x, y, z, px, py, pzRef+dpz))
        return 1;
      totalCharge += charge;
    }
  }
  
  if (!SDDS_SetParameters(&SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,
                          "Charge", fabs(totalCharge/1e9), NULL))
    SDDS_Bomb("problem setting Charge value");

  if (row%maxRows && !SDDS_WritePage(&SDDSout))
    SDDS_Bomb("problem writing rows with particle data");

  if (!SDDS_Terminate(&SDDSout))
    SDDS_Bomb("problem terminating output file");
    
  return 0;
}


long SetUpOutputFile(SDDS_DATASET *SDDSout, char *outputfile)
{
  if (!SDDS_InitializeOutput(SDDSout, SDDS_BINARY, 1, NULL, NULL, outputfile) ||
      SDDS_DefineColumn(SDDSout, "x", NULL, "m", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineColumn(SDDSout, "xp", NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineColumn(SDDSout, "y", NULL, "m", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineColumn(SDDSout, "yp", NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineColumn(SDDSout, "t", NULL, "s", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineColumn(SDDSout, "p", NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineParameter(SDDSout, "Charge", NULL, "C", NULL, NULL, SDDS_DOUBLE, NULL)<0 ||
      !SDDS_WriteLayout(SDDSout)) 
    return 0;
  return 1;
}

#define K1 (e_mks/me_mks/sqr(c_mks))

long addParticleToFile(SDDS_DATASET *SDDSout, long *row, long *maxRows,
                       double x, double y, double z, double px, double py, double pz)
{
  double xp, yp;
  double p2, gamma, betaz;
  
  if (!pz)
    return 0;
  xp = px/pz;
  yp = py/pz;

  x -= z*xp;
  y -= z*yp;
  p2 = (sqr(px)+sqr(py)+sqr(pz))*sqr(K1);
  gamma = sqrt(p2+1);
  betaz = (pz*K1)/gamma;
  if (!SDDS_SetRowValues(SDDSout, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, 
                         *row,
                         0, x, 1, xp, 2, y, 3, yp, 4, -z/(betaz*c_mks), 5, sqrt(p2), -1)) {
    SDDS_Bomb("problem setting rows with particle data");
    return 0;
  }
  *row += 1;
  if (*row==*maxRows) {
    *maxRows += MAX_ROW_INCREMENT;
    if (!SDDS_LengthenTable(SDDSout, *maxRows)) {
      SDDS_Bomb("problem writing rows with particle data");
      return 0;
    }
  }
  return 1;
}


