/*CopyrightNotice001

*****************************************************************
                          COPYRIGHT NOTIFICATION
*****************************************************************

THE FOLLOWING IS A NOTICE OF COPYRIGHT, AVAILABILITY OF THE CODE,
AND DISCLAIMER WHICH MUST BE INCLUDED IN THE PROLOGUE OF THE CODE
AND IN ALL SOURCE LISTINGS OF THE CODE.
 
(C)  COPYRIGHT 1998 UNIVERSITY OF CHICAGO
 
Argonne National Laboratory (ANL), with facilities in the States of 
Illinois and Idaho, is owned by the United States Government, and
operated by the University of Chicago under provision of a contract
with the Department of Energy.

Portions of this material resulted from work developed under a U.S.
Government contract and are subject to the following license:  For
a period of five years from December 30, 1998, the Government is
granted for itself and others acting on its behalf a paid-up,
nonexclusive, irrevocable worldwide license in this computer
software to reproduce, prepare derivative works, and perform
publicly and display publicly.  With the approval of DOE, this
period may be renewed for two additional five year periods. 
Following the expiration of this period or periods, the Government
is granted for itself and others acting on its behalf, a paid-up,
nonexclusive, irrevocable worldwide license in this computer
software to reproduce, prepare derivative works, distribute copies
to the public, perform publicly and display publicly, and to permit
others to do so.

*****************************************************************
                                DISCLAIMER
*****************************************************************

NEITHER THE UNITED STATES GOVERNMENT NOR ANY AGENCY THEREOF, NOR
THE UNIVERSITY OF CHICAGO, NOR ANY OF THEIR EMPLOYEES OR OFFICERS,
MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL
LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR
USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT, OR PROCESS
DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE PRIVATELY
OWNED RIGHTS.  

*****************************************************************
LICENSING INQUIRIES MAY BE DIRECTED TO THE INDUSTRIAL TECHNOLOGY
DEVELOPMENT CENTER AT ARGONNE NATIONAL LABORATORY.

CopyrightNotice001*/

/* program: sddsmatchtwiss.c
 * purpose: take particle coordinate input and convert the Twiss parameters
 *          to user-supplied values.
 *
 * Michael Borland, 2000
 *
 $Log: not supported by cvs2svn $
 Revision 1.1  2000/08/14 21:44:56  borland
 First version.

 */
#include "mdb.h"
#include "scan.h"
#include "SDDS.h"

#define DEBUG 0

#define SET_PIPE 0
#define SET_XPLANE 1
#define SET_YPLANE 2
#define SET_NOWARNINGS 3
#define N_OPTIONS 4

char *option[N_OPTIONS] = {
  "pipe", "xplane", "yplane", "nowarnings"
} ;

char *USAGE="sddsmatchtwiss [-pipe=[input][,output]] [<SDDSinputfile>] [<SDDSoutputfile>]\n\
  [-xPlane=[beta=<meters>,alpha=<value>][,etaValue=<meters>][,etaSlope=<value>]]\n\
  [-yPlane=[beta=<meters>,alpha=<value>][,etaValue=<meters>][,etaSlope=<value>]]\n\
  [-nowarnings]\n\
The input file must have columns x, xp, y, yp, and p; for example, an elegant\n\
beam output file is acceptable.\n\
beta and alpha must be given together, or omitted together.\n\
etaValue and etaSlope may be given individually or together.  If etaValue is not given,\n\
the coordinates have no dispersion adjustment.  If etaSlope is not given, then\n\
slopes have no dispersion adjustment.\n\n\
Program by Michael Borland.  (This is version 1, August 2000.)\n";

typedef struct {
  double beta, alpha, eta, etap;
  unsigned long flags;
#define BETA_GIVEN  0x0001UL
#define ALPHA_GIVEN 0x0002UL
#define ETA_GIVEN   0x0004UL
#define ETAP_GIVEN  0x0008UL
} PLANE_SPEC;

long PerformTransformation(double *x, double *xp, double *p, long rows, PLANE_SPEC match);
long check_sdds_beam_column(SDDS_TABLE *SDDS_table, char *name, char *units);


main(int argc, char **argv)
{
  SDDS_DATASET SDDSin, SDDSout;
  char *inputfile, *outputfile;
  long i_arg, row, rows, readCode, noWarnings;
  SCANNED_ARG *s_arg;
  unsigned long pipeFlags;
  PLANE_SPEC xSpec, ySpec;
  double *x, *xp, *y, *yp, *p;
  
  SDDS_RegisterProgramName(argv[0]);
  argc = scanargs(&s_arg, argc, argv);
  if (argc<2) 
    bomb(NULL, USAGE);

  inputfile = outputfile = NULL;
  pipeFlags = noWarnings = 0;
  xSpec.flags = ySpec.flags = 0;
  
  for (i_arg=1; i_arg<argc; i_arg++) {
    if (s_arg[i_arg].arg_type==OPTION) {
      switch (match_string(s_arg[i_arg].list[0], option, N_OPTIONS, 0)) {
      case SET_XPLANE:
        xSpec.flags = 0;
        s_arg[i_arg].n_items--;
        if (!scanItemList(&xSpec.flags, s_arg[i_arg].list+1, &s_arg[i_arg].n_items, 0,
                          "beta", SDDS_DOUBLE, &xSpec.beta, 1, BETA_GIVEN,
                          "alpha", SDDS_DOUBLE, &xSpec.alpha, 1, ALPHA_GIVEN,
                          "etavalue", SDDS_DOUBLE, &xSpec.eta, 1, ETA_GIVEN,
                          "etaslope", SDDS_DOUBLE, &xSpec.etap, 1, ETAP_GIVEN,
                          NULL) ||
            (xSpec.flags&BETA_GIVEN && !(xSpec.flags&ALPHA_GIVEN)) ||
            (!(xSpec.flags&BETA_GIVEN) && xSpec.flags&ALPHA_GIVEN) ||
            (xSpec.flags&BETA_GIVEN && xSpec.beta<=0))
          SDDS_Bomb("invalid -xPlane syntax/values---watch out for abbreviations of etaValue and etaSlope");
        break;
      case SET_YPLANE:
        ySpec.flags = 0;
        s_arg[i_arg].n_items--;
        if (!scanItemList(&ySpec.flags, s_arg[i_arg].list+1, &s_arg[i_arg].n_items, 0,
                          "beta", SDDS_DOUBLE, &ySpec.beta, 1, BETA_GIVEN,
                          "alpha", SDDS_DOUBLE, &ySpec.alpha, 1, ALPHA_GIVEN,
                          "etavalue", SDDS_DOUBLE, &ySpec.eta, 1, ETA_GIVEN,
                          "etaslope", SDDS_DOUBLE, &ySpec.etap, 1, ETAP_GIVEN,
                          NULL) ||
            (ySpec.flags&BETA_GIVEN && !(ySpec.flags&ALPHA_GIVEN)) ||
            (!(ySpec.flags&BETA_GIVEN) && ySpec.flags&ALPHA_GIVEN) ||
            (ySpec.flags&BETA_GIVEN && ySpec.beta<=0))
          SDDS_Bomb("invalid -yPlane syntax/values---watch out for abbreviations of etaValue and etaSlope");
        break;
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

  processFilenames("sddsmatchbeam", &inputfile, &outputfile, pipeFlags, noWarnings, NULL);

  if (!SDDS_InitializeInput(&SDDSin, inputfile))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);

  /* check the input file for valid data */
  if (!check_sdds_beam_column(&SDDSin, "x", "m") || !check_sdds_beam_column(&SDDSin, "y", "m") ||
      !check_sdds_beam_column(&SDDSin, "xp", NULL) || !check_sdds_beam_column(&SDDSin, "yp", NULL) ||
      (!check_sdds_beam_column(&SDDSin, "p", "m$be$nc") && !check_sdds_beam_column(&SDDSin, "p", NULL))) {
    fprintf(stderr, 
            "sddsmatchtwiss: one or more data quantities (x, xp, y, yp, p) have the wrong units or are not present in %s", 
            inputfile);
    exit(1);
  }

  if (!SDDS_InitializeCopy(&SDDSout, &SDDSin, outputfile, "w") ||
      !SDDS_SaveLayout(&SDDSout) || !SDDS_WriteLayout(&SDDSout))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  
  while ((readCode=SDDS_ReadPage(&SDDSin))>0) {
    if ((rows=SDDS_RowCount(&SDDSin))==0) {
      if (!SDDS_CopyPage(&SDDSin, &SDDSout) || !SDDS_WritePage(&SDDSout))
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      continue;
    }
    if (!(x = SDDS_GetColumnInDoubles(&SDDSin, "x")) ||
        !(xp = SDDS_GetColumnInDoubles(&SDDSin, "xp")) ||
        !(y = SDDS_GetColumnInDoubles(&SDDSin, "y")) ||
        !(yp = SDDS_GetColumnInDoubles(&SDDSin, "yp")) ||
        !(p = SDDS_GetColumnInDoubles(&SDDSin, "p")))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    PerformTransformation(x, xp, p, rows, xSpec);
    PerformTransformation(y, yp, p, rows, ySpec);
    if (!SDDS_CopyPage(&SDDSout, &SDDSin) ||
        !(SDDS_SetColumnFromDoubles(&SDDSout, SDDS_SET_BY_NAME, x, rows, "x")) ||
        !(SDDS_SetColumnFromDoubles(&SDDSout, SDDS_SET_BY_NAME, y, rows, "y")) ||
        !(SDDS_SetColumnFromDoubles(&SDDSout, SDDS_SET_BY_NAME, xp, rows, "xp")) ||
        !(SDDS_SetColumnFromDoubles(&SDDSout, SDDS_SET_BY_NAME, yp, rows, "yp")) ||
        !(SDDS_SetColumnFromDoubles(&SDDSout, SDDS_SET_BY_NAME, p, rows, "p")) ||
        !SDDS_WritePage(&SDDSout))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  if (readCode==0)
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
}

long PerformTransformation(double *x, double *xp, double *p, long rows, PLANE_SPEC match)
{
  long i;
  double pAve;
  double S11, S22, S12, S16, S26, S66;
  double R11, R12, R21, R22;
  double x0, xp0;
  double emit, beta1, beta2, alpha1, alpha2;
  double eta1, etap1, eta2, etap2;
  
  pAve = arithmeticAverage(p, rows);
  for (i=0; i<rows; i++)
    p[i] = p[i]/pAve - 1;
  computeCorrelations(&S11, &S16, &S66, x, p, rows);
  eta1 = S16/S66;
  for (i=0; i<rows; i++)
    x[i] -= p[i]*eta1;
  computeCorrelations(&S22, &S26, &S66, xp, p, rows);
  etap1 = S26/S66;
  for (i=0; i<rows; i++)
    xp[i] -= p[i]*etap1;
  
  if (match.flags&BETA_GIVEN) {
    computeCorrelations(&S11, &S12, &S22, x, xp,rows);
    if ((emit = S11*S22-sqr(S12))<0) {
      SDDS_Bomb("emittance is zero");
    }
    else
      emit = sqrt(emit);
    beta1 = S11/emit;
    alpha1 = -S12/emit;
    beta2 = match.beta;
    alpha2 = match.alpha;
    R11 = beta2/sqrt(beta1*beta2);
    R12 = 0;
    R21 = (alpha1-alpha2)/sqrt(beta1*beta2);
    R22 = beta1/sqrt(beta1*beta2);
    for (i=0; i<rows; i++) {
      x0 = x[i];
      xp0 = xp[i];
      x[i] = R11*x0 + R12*xp0;
      xp[i] = R21*x0 + R22*xp0;
    }
  }

  if (match.flags&ETA_GIVEN)
    eta2 = match.eta;
  else
    eta2 = eta1;
  if (match.flags&ETAP_GIVEN)
    etap2 = match.etap;
  else
    etap2 = etap1;
  for (i=0; i<rows; i++) {
    x[i] += eta2*p[i];
    xp[i] += etap2*p[i];
    p[i] = (p[i]+1)*pAve;
  }
  return 1;
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

