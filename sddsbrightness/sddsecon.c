/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

#include "mdb.h"
#include "SDDS.h"
#include "scan.h"

#define SET_PIPE 0
#define SET_DE1 1
#define SET_OUTPUTCOLUMN 2
#define N_OPTIONS 3

#define NPPSIGMA 6

char *option[N_OPTIONS] = {"pipe", "de1", "outputColumn"};

char *USAGE = "sddsecon [<input>] [<output>] [-pipe=[in][,out]] -de1=<value> -outputColumn=<name>\n\
<input>        The input file is expected to contain two columns called 'e' and 'spec'.\n\
               'e' is the photon energy (ev,keV)\n\
               'spec' is the spectrum \n\
-de1           Detector energy resolution at lower energy (ev,keV)\n\
-outputColumn  The name of the column containing the results.\n\
RESTRICTIONS:\n\
 The energy points need to be equally spaced in ascending order.\n\
 Due to the convolution performed, some data points will be set to\n\
 zero at each end of the spectrum. The energy range with zeroed data\n\
 will be bigger at the end of the spectrum due to the transformation used,\n\
 and it depends on the value of the beam energy spread. The larger the\n\
 beam energy spread, the larger the range with zeroed data.\n\
PROCEDURE:\n\
 The spectrum is convolved with a Gaussian function with a constant width.\n\n\
Program by Robert Soliday. ("__DATE__")\n";

double* econ(double *e, double *spec, long rows, double de1);
long checkMonotonicity(double *indepValue, long rows);

int main(int argc, char **argv)
{
  SCANNED_ARG *s_arg;
  long i_arg;

  double de1=0;
  char *input=NULL, *output=NULL, *outputColumn=NULL;
  unsigned long pipeFlags=0;
  long noWarnings=0, tmpfile_used=0, rows;

  SDDS_DATASET SDDS_input, SDDS_output;
  double *e, *spec, *specResults=NULL;

  SDDS_RegisterProgramName(argv[0]);
  argc = scanargs(&s_arg, argc, argv);
  if (argc<5) {
    fprintf(stderr, "%s", USAGE);
    return(1);
  }
  /* parse the command line */
  for (i_arg=1; i_arg<argc; i_arg++) {
    if (s_arg[i_arg].arg_type==OPTION) {
      switch (match_string(s_arg[i_arg].list[0], option, N_OPTIONS, 0)) {
      case SET_DE1:
	if (s_arg[i_arg].n_items!=2) {
	  fprintf(stderr, "error: invalid -de1 syntax\n");
	  return(1);
	}
	if (sscanf(s_arg[i_arg].list[1], "%lf", &de1) != 1) {
	  fprintf(stderr, "error: invalid -de1 syntax or value\n");
	  return(1);
	}
	break;
      case SET_OUTPUTCOLUMN:
	if (s_arg[i_arg].n_items!=2) {
	  fprintf(stderr, "error: invalid -de1 syntax\n");
	  return(1);
	}
        outputColumn = s_arg[i_arg].list[1];
	break;
      default:
	fprintf(stderr, "error: unknown switch: %s\n", s_arg[i_arg].list[0]);
	return(1);
      }
    } else {
      if (input == NULL) {
	input = s_arg[i_arg].list[0];
      } else if (output == NULL) {
	output = s_arg[i_arg].list[0];
      } else {
	fprintf(stderr, "too many filenames\n");
	return(1);
      }
    }
  }
  processFilenames("sddsecon", &input, &output, pipeFlags, noWarnings, &tmpfile_used);

  /* open the input file */
  if (!SDDS_InitializeInput(&SDDS_input, input)) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    return(1);
  }

  if (SDDS_ReadPage(&SDDS_input) != 1) {
    fprintf(stderr, "error: no data found in input file\n");
    return(1);
  }
  rows = SDDS_RowCount(&SDDS_input);

  if (!(e = SDDS_GetColumnInDoubles(&SDDS_input, "e")) ||
      !(spec = SDDS_GetColumnInDoubles(&SDDS_input, "spec"))) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    return(1);
  }
  /* close the input file */
  if (!SDDS_Terminate(&SDDS_input)) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    return(1);
  }

  specResults = econ(e, spec, rows, de1);
  if (specResults == NULL) {
    return(1);
  }
  /* open the output file */
  if (!SDDS_InitializeOutput(&SDDS_output, SDDS_BINARY, 0, NULL, NULL, output)) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    return(1);
  }
  if (!SDDS_DefineSimpleColumn(&SDDS_output, "e", NULL, SDDS_DOUBLE)) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    return(1);
  }
  if (!SDDS_DefineSimpleColumn(&SDDS_output, "spec", NULL, SDDS_DOUBLE)) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    return(1);
  }
  if (!SDDS_DefineSimpleColumn(&SDDS_output, outputColumn, NULL, SDDS_DOUBLE)) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    return(1);
  }
  if (!SDDS_WriteLayout(&SDDS_output)) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    return(1);
  }
  if (!SDDS_StartTable(&SDDS_output, rows)) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    return(1);
  }
  if (!SDDS_SetColumn(&SDDS_output, SDDS_SET_BY_NAME, e, rows, "e")) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    return(1);
  }
  if (!SDDS_SetColumn(&SDDS_output, SDDS_SET_BY_NAME, spec, rows, "spec")) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    return(1);
  }
  if (!SDDS_SetColumn(&SDDS_output, SDDS_SET_BY_NAME, specResults, rows, outputColumn)) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    return(1);
  }
  if (!SDDS_WriteTable(&SDDS_output)) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    return(1);
  }
  /* close the input file */
  if (!SDDS_Terminate(&SDDS_output)) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    return(1);
  }

  free(e);
  free(spec);
  free(specResults);
  free_scanargs(&s_arg,argc);
  return(0);
}

double* econ(double *e, double *spec, long rows, double de1) {
  double *cf, *cg, *specg;
  long ng, i, monotonicity;

  double q, sigf;
  unsigned long interpCode=0;
  OUTRANGE_CONTROL aboveRange, belowRange;

  long nsigma, nsigma2, n_pts, t, k2;
  double conv, xvar, *gaus, gausSum, *specc, *specResults;

  aboveRange.flags = belowRange.flags = OUTRANGE_ABORT;

  cf = malloc(sizeof(double) * rows);
  q = (rows - 1) / log(e[rows-1] / e[0]);
  for (i = 0; i < rows; i++) {
    cf[i] = q * log(e[i]/e[0]);
  }
  sigf = 2.0 * de1 * q;
  if (sigf < .01) {
    /*Return original spec array here*/
    specResults = malloc(sizeof(double) * rows);
    for (i = 0; i < rows; i++) {
      specResults[i] = spec[i];
    }
    return(specResults);
  }
  ng = ceil(NPPSIGMA * rows / sigf);
  if (rows > ng)
    ng = rows;
  cg = malloc(sizeof(double) * ng);
  specg = malloc(sizeof(double) * ng);
  for (i = 0; i < ng; i++) {
    cg[i] = i * (rows - 1) / (ng - 1);
  }
  if ((monotonicity=checkMonotonicity(cf, rows))==0) {
    fprintf(stderr, "independent data values do not change monotonically or repeated independent values exist\n");
    return(NULL);
  }
  for (i = 0; i < ng; i++) {
    specg[i] = interpolate(spec, cf, rows, cg[i], &belowRange, &aboveRange, 1, &interpCode, monotonicity);
    if (interpCode) {
      if (interpCode&OUTRANGE_ABORT) {
        fprintf(stderr, "error: can't do interpolation, value %lf is out of range.\n", cg[i]);
        return(NULL);
      }
    }
  }
  nsigma = 3.0;
  nsigma2 = nsigma * 2;
  conv = (cg[ng-1] - cg[0]) / (ng - 1);
  n_pts = ceil(nsigma2 * sigf / conv);
  if (2 > n_pts)
    n_pts = 2;
  if ((ng - 2) < n_pts)
    n_pts = ng - 2;
  if ((long)(n_pts/(long)2 * 2) == n_pts)
    n_pts += 1;
  gaus = malloc(sizeof(double) * n_pts);
  gausSum = 0;
  for (i = 0; i < n_pts; i++) {
    xvar = (i / (n_pts - 1.0) - .5) * n_pts;
    gaus[i] = exp(-.5 * pow((xvar / (sigf / conv)),2));
    gausSum += gaus[i];
  }

  specc = malloc(sizeof(double) * ng);
  k2 = n_pts / (long)2;
  for (t=0; t < ng; t++) {
    specc[t]=0;
    for (i=0; i <= n_pts - 1; i++) {
      if ((k2 <= t) && (t <= ng - k2 - 1)) {
        specc[t] += specg[t + i - k2] * gaus[i];
      }
    }
    specc[t] = specc[t] / gausSum;
  }

  specResults = malloc(sizeof(double) * rows);
  if ((monotonicity=checkMonotonicity(cg, ng))==0) {
    fprintf(stderr, "independent data values do not change monotonically or repeated independent values exist\n");
    return(NULL);
  }
  for (i = 0; i < rows; i++) {
    specResults[i] = interpolate(specc, cg, ng, cf[i], &belowRange, &aboveRange, 1, &interpCode, monotonicity);
    if (interpCode) {
      if (interpCode&OUTRANGE_ABORT) {
        fprintf(stderr, "error: can't do interpolation, value %lf is out of range.\n", cf[i]);
        return(NULL);
      }
    }
  }
  free(cf);
  free(cg);
  free(specg);
  free(gaus);
  free(specc);
  return(specResults);
}

long checkMonotonicity(double *indepValue, long rows)
{
  if (rows==1)
    return 1;
  if (indepValue[rows-1]>indepValue[0]) {
    while (--rows>0)
      if (indepValue[rows]<=indepValue[rows-1])
        return 0;
    return 1;
  }
  else {
    while (--rows>0)
      if (indepValue[rows]>=indepValue[rows-1])
        return 0;
    return -1;
  }
}
