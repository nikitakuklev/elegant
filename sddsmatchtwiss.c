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
 Revision 1.8  2001/10/15 15:42:33  soliday
 Cleaned up for WIN32.

 Revision 1.7  2001/05/30 15:53:27  borland
 -zPlane option now allows specifying alphaz or correlation.

 Revision 1.6  2001/05/15 16:52:54  borland
 Added to -zPlane option the ability to change the central momentum.

 Revision 1.5  2001/05/08 14:07:22  borland
 Added -ztransform option, which permits changing the longitudinal parameters.

 Revision 1.4  2001/04/06 14:44:55  borland
 Added ability to change the emittance of the beam in both planes.

 Revision 1.3  2001/01/08 20:22:16  borland
 Added -oneTransform option, which allows computing the transform only for
 the first page, and using it for all pages.
 Also, added hidden -verbose option.

 Revision 1.2  2000/08/17 13:45:56  borland
 Fixed errors in usage message.

 Revision 1.1  2000/08/14 21:44:56  borland
 First version.

 */
#include "mdb.h"
#include "scan.h"
#include "SDDS.h"

#define SET_PIPE 0
#define SET_XPLANE 1
#define SET_YPLANE 2
#define SET_NOWARNINGS 3
#define SET_ONETRANSFORM 4
#define SET_VERBOSE 5
#define SET_ZPLANE 6
#define N_OPTIONS 7

char *option[N_OPTIONS] = {
  "pipe", "xplane", "yplane", "nowarnings", "oneTransform", "verbose", 
  "zplane",
} ;

char *USAGE="sddsmatchtwiss [-pipe=[input][,output]] [<SDDSinputfile>] [<SDDSoutputfile>]\n\
  [-xPlane=[beta=<meters>,alpha=<value>][nemittance=<meters>,][,etaValue=<meters>][,etaSlope=<value>]]\n\
  [-yPlane=[beta=<meters>,alpha=<value>][nemittance=<meters>,][,etaValue=<meters>][,etaSlope=<value>]]\n\
  [-zPlane=[deltaStDev=<value>][,tStDev=<seconds>][,{correlation=<seconds>|alpha=<value>}][,betaGamma=<central-value>]]\n\
  [-nowarnings] [-oneTransform]\n\
The input file must have columns x, xp, y, yp, and p; for example, an elegant\n\
beam output file is acceptable.\n\
beta and alpha must be given together, or omitted together.\n\
etaValue and etaSlope may be given individually or together.  If etaValue is not given,\n\
the coordinates have no dispersion adjustment.  If etaSlope is not given, then\n\
slopes have no dispersion adjustment.\n\
If -oneTransform is given, then the transformation is computed for the first page only,\n\
then reused for all subsequent pages.\n\n\
Program by Michael Borland.  (This is version 4, May 2001.)\n";

typedef struct {
  double beta, alpha, eta, etap, normEmittance;
  unsigned long flags;
  double R11, R12, R21, R22;
  double etaBeam, etapBeam;
#define BETA_GIVEN  0x0001UL
#define ALPHA_GIVEN 0x0002UL
#define ETA_GIVEN   0x0004UL
#define ETAP_GIVEN  0x0008UL
#define NEMIT_GIVEN 0x0010UL
} PLANE_SPEC;

typedef struct {
  double deltaStDev, tStDev, correlation, alpha, betaGamma;
  unsigned long flags;
  double R11, R12, R21, R22;
#define DELTASTDEV_GIVEN  0x0001UL
#define TSTDEV_GIVEN      0x0002UL
#define CORRELATION_GIVEN 0x0004UL
#define BETAGAMMA_GIVEN   0x0008UL
#define ALPHAZ_GIVEN      0x0010UL
} ZPLANE_SPEC;

long PerformTransformation(double *x, double *xp, double *p, long rows, PLANE_SPEC *match,
                           long compute);
long PerformZTransformation(double *t, double *p, double *x, double *xp,
                            double *y, double *yp, long rows, ZPLANE_SPEC *match, long compute);

long check_sdds_beam_column(SDDS_TABLE *SDDS_table, char *name, char *units);


int main(int argc, char **argv)
{
  SDDS_DATASET SDDSin, SDDSout;
  char *inputfile, *outputfile;
  long i_arg, rows, readCode, noWarnings;
  SCANNED_ARG *s_arg;
  unsigned long pipeFlags;
  PLANE_SPEC xSpec, ySpec;
  ZPLANE_SPEC zSpec;
  double *x, *xp=NULL, *y=NULL, *yp=NULL, *p=NULL, *t=NULL;
  long oneTransform, verbose;
  
  SDDS_RegisterProgramName(argv[0]);
  argc = scanargs(&s_arg, argc, argv);
  if (argc<2) 
    bomb(NULL, USAGE);

  inputfile = outputfile = NULL;
  pipeFlags = noWarnings = 0;
  xSpec.flags = ySpec.flags = zSpec.flags = 0;
  verbose = oneTransform = 0;
  
  for (i_arg=1; i_arg<argc; i_arg++) {
    if (s_arg[i_arg].arg_type==OPTION) {
      switch (match_string(s_arg[i_arg].list[0], option, N_OPTIONS, 0)) {
      case SET_XPLANE:
        xSpec.flags = 0;
        s_arg[i_arg].n_items--;
        if (!scanItemList(&xSpec.flags, s_arg[i_arg].list+1, &s_arg[i_arg].n_items, 0,
                          "beta", SDDS_DOUBLE, &xSpec.beta, 1, BETA_GIVEN,
                          "alpha", SDDS_DOUBLE, &xSpec.alpha, 1, ALPHA_GIVEN,
                          "nemittance", SDDS_DOUBLE, &xSpec.normEmittance, 1, NEMIT_GIVEN,
                          "etavalue", SDDS_DOUBLE, &xSpec.eta, 1, ETA_GIVEN,
                          "etaslope", SDDS_DOUBLE, &xSpec.etap, 1, ETAP_GIVEN,
                          NULL) ||
            (xSpec.flags&BETA_GIVEN && !(xSpec.flags&ALPHA_GIVEN)) ||
            (!(xSpec.flags&BETA_GIVEN) && xSpec.flags&ALPHA_GIVEN) ||
            (xSpec.flags&BETA_GIVEN && xSpec.beta<=0) || 
            (xSpec.flags&NEMIT_GIVEN && xSpec.normEmittance<=0))
          SDDS_Bomb("invalid -xPlane syntax/values---watch out for abbreviations of etaValue and etaSlope");
        break;
      case SET_YPLANE:
        ySpec.flags = 0;
        s_arg[i_arg].n_items--;
        if (!scanItemList(&ySpec.flags, s_arg[i_arg].list+1, &s_arg[i_arg].n_items, 0,
                          "beta", SDDS_DOUBLE, &ySpec.beta, 1, BETA_GIVEN,
                          "alpha", SDDS_DOUBLE, &ySpec.alpha, 1, ALPHA_GIVEN,
                          "nemittance", SDDS_DOUBLE, &ySpec.normEmittance, 1, NEMIT_GIVEN,
                          "etavalue", SDDS_DOUBLE, &ySpec.eta, 1, ETA_GIVEN,
                          "etaslope", SDDS_DOUBLE, &ySpec.etap, 1, ETAP_GIVEN,
                          NULL) ||
            (ySpec.flags&BETA_GIVEN && !(ySpec.flags&ALPHA_GIVEN)) ||
            (!(ySpec.flags&BETA_GIVEN) && ySpec.flags&ALPHA_GIVEN) ||
            (ySpec.flags&BETA_GIVEN && ySpec.beta<=0) ||
            (ySpec.flags&NEMIT_GIVEN && ySpec.normEmittance<=0))
          SDDS_Bomb("invalid -yPlane syntax/values---watch out for abbreviations of etaValue and etaSlope");
        break;
      case SET_ZPLANE:
        zSpec.flags = 0;
        s_arg[i_arg].n_items--;
        if (!scanItemList(&zSpec.flags, s_arg[i_arg].list+1, &s_arg[i_arg].n_items, 0,
                          "tStDev", SDDS_DOUBLE, &zSpec.tStDev, 1, TSTDEV_GIVEN,
                          "deltaStDev", SDDS_DOUBLE, &zSpec.deltaStDev, 1, DELTASTDEV_GIVEN,
                          "correlation", SDDS_DOUBLE, &zSpec.correlation, 1, CORRELATION_GIVEN,
                          "alpha", SDDS_DOUBLE, &zSpec.alpha, 1, ALPHAZ_GIVEN,
                          "betagamma", SDDS_DOUBLE, &zSpec.betaGamma, 1, BETAGAMMA_GIVEN,
                          NULL) ||
            bitsSet(zSpec.flags)<1 ||
            (zSpec.flags&ALPHAZ_GIVEN && zSpec.flags&CORRELATION_GIVEN) ||
            (zSpec.flags&TSTDEV_GIVEN &&  zSpec.tStDev<0) ||
            (zSpec.flags&DELTASTDEV_GIVEN && zSpec.deltaStDev<0) ||
            (zSpec.flags&BETAGAMMA_GIVEN && zSpec.betaGamma<=0))
          SDDS_Bomb("invalid -zPlane syntax/values");
        break;
      case SET_PIPE:
        if (!processPipeOption(s_arg[i_arg].list+1, s_arg[i_arg].n_items-1, &pipeFlags))
          SDDS_Bomb("invalid -pipe syntax");
        break;
      case SET_NOWARNINGS:
        noWarnings = 1;
        break;
      case SET_ONETRANSFORM:
        oneTransform = 1;
        break;
      case SET_VERBOSE:
        verbose = 1;
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
        !(t = SDDS_GetColumnInDoubles(&SDDSin, "t")) ||
        !(p = SDDS_GetColumnInDoubles(&SDDSin, "p")))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    PerformTransformation(x, xp, p, rows, &xSpec, readCode==1?1:!oneTransform);
    PerformTransformation(y, yp, p, rows, &ySpec, readCode==1?1:!oneTransform);
    PerformZTransformation(t, p, x, xp, y, yp, rows, &zSpec, readCode==1?1:!oneTransform);
    if (verbose) {
      if (xSpec.flags)
        fprintf(stderr, "x transformation: %le, %le, %le, %le, %le, %le\n",
                xSpec.R11, xSpec.R12, xSpec.R21, xSpec.R22, xSpec.etaBeam,
                xSpec.etapBeam);
      if (ySpec.flags)
        fprintf(stderr, "y transformation: %le, %le, %le, %le, %le, %le\n",
                ySpec.R11, ySpec.R12, ySpec.R21, ySpec.R22, ySpec.etaBeam,
                ySpec.etapBeam);
      if (zSpec.flags)
        fprintf(stderr, "z transformation: %le, %le, %le, %le\n",
                zSpec.R11, zSpec.R12, zSpec.R21, zSpec.R22);
    }
    if (!SDDS_CopyPage(&SDDSout, &SDDSin) ||
        !(SDDS_SetColumnFromDoubles(&SDDSout, SDDS_SET_BY_NAME, x, rows, "x")) ||
        !(SDDS_SetColumnFromDoubles(&SDDSout, SDDS_SET_BY_NAME, y, rows, "y")) ||
        !(SDDS_SetColumnFromDoubles(&SDDSout, SDDS_SET_BY_NAME, xp, rows, "xp")) ||
        !(SDDS_SetColumnFromDoubles(&SDDSout, SDDS_SET_BY_NAME, yp, rows, "yp")) ||
        !(SDDS_SetColumnFromDoubles(&SDDSout, SDDS_SET_BY_NAME, t, rows, "t")) ||
        !(SDDS_SetColumnFromDoubles(&SDDSout, SDDS_SET_BY_NAME, p, rows, "p")) ||
        !SDDS_WritePage(&SDDSout))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  if (readCode==0)
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  return 0;
}

long PerformTransformation(double *x, double *xp, double *p, long rows, PLANE_SPEC *match,
                           long computeTransform)
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

  if (computeTransform) {
    computeCorrelations(&S11, &S16, &S66, x, p, rows);
    match->etaBeam = eta1 = S16/S66;
    computeCorrelations(&S22, &S26, &S66, xp, p, rows);
    match->etapBeam = etap1 = S26/S66;
  } else {
    eta1 = match->etaBeam;
    etap1 = match->etapBeam;
  }
  
  for (i=0; i<rows; i++)
    x[i] -= p[i]*eta1;
  for (i=0; i<rows; i++)
    xp[i] -= p[i]*etap1;
  
  if (match->flags&BETA_GIVEN || match->flags&NEMIT_GIVEN) {
    if (computeTransform || match->flags&NEMIT_GIVEN) {
      computeCorrelations(&S11, &S12, &S22, x, xp,rows);
      if ((emit = S11*S22-sqr(S12))<=0) {
        SDDS_Bomb("emittance is zero");
      }
      else
        emit = sqrt(emit);
      beta1 = S11/emit;
      alpha1 = -S12/emit;
      if (match->flags&NEMIT_GIVEN) {
        match->R11 = R11 = match->R22 = R22 = 
          sqrt(match->normEmittance/(emit*pAve));
        match->R12 = R12 = match->R21 = R21 = 0;
        beta2 = beta1;
        alpha2 = alpha1;
      } else {
        beta2 = match->beta;
        alpha2 = match->alpha;
        match->R11 = R11 = beta2/sqrt(beta1*beta2);
        match->R12 = R12 = 0;
        match->R21 = R21 = (alpha1-alpha2)/sqrt(beta1*beta2);
        match->R22 = R22 = beta1/sqrt(beta1*beta2);
      }
    }
    else {
      R11 = match->R11;
      R12 = match->R12;
      R21 = match->R21;
      R22 = match->R22;
    }
    for (i=0; i<rows; i++) {
      x0 = x[i];
      xp0 = xp[i];
      x[i] = R11*x0 + R12*xp0;
      xp[i] = R21*x0 + R22*xp0;
    }
  }

  if (match->flags&ETA_GIVEN)
    eta2 = match->eta;
  else
    eta2 = eta1;
  if (match->flags&ETAP_GIVEN)
    etap2 = match->etap;
  else
    etap2 = etap1;
  for (i=0; i<rows; i++) {
    x[i] += eta2*p[i];
    xp[i] += etap2*p[i];
    p[i] = (p[i]+1)*pAve;
  }
  return 1;
}

long PerformZTransformation(double *t, double *p, 
                            double *x, double *xp,
                            double *y, double *yp,
                            long rows, ZPLANE_SPEC *match,
                            long computeTransform)
{
  long i;
  double pAve, tAve;
  double S11, S22, S12;
  double R11, R12, R21, R22;
  double delta0, t0;
  double emit1, beta1, beta2, alpha1, alpha2;
  double emit2, ratio;

  pAve = arithmeticAverage(p, rows);
  for (i=0; i<rows; i++)
    p[i] = p[i]/pAve - 1;
  tAve = arithmeticAverage(t, rows);
  for (i=0; i<rows; i++)
    t[i] = t[i] - tAve;

  if (computeTransform) {
    computeCorrelations(&S11, &S12, &S22, t, p, rows);
    if ((emit1 = S11*S22-sqr(S12))<=0)
      SDDS_Bomb("longitudinal emittance is zero");
    else
      emit1 = sqrt(emit1);
    beta1 = S11/emit1;
    alpha1 = -S12/emit1;

    if (match->flags&TSTDEV_GIVEN)
      S11 = sqr(match->tStDev);
    if (match->flags&DELTASTDEV_GIVEN)
      S22 = sqr(match->deltaStDev);
    if (match->flags&CORRELATION_GIVEN) 
      S12 = match->correlation*sqrt(S11*S22);
    if (match->flags&ALPHAZ_GIVEN)
      S12 = -match->alpha*sqrt(S11*S22)/sqrt(1+sqr(match->alpha));
    if ((emit2 = S11*S22-sqr(S12))<=0)
      SDDS_Bomb("longitudinal emittance is zero");
    else
      emit2 = sqrt(emit2);
    beta2 = S11/emit2;
    if (!(match->flags&ALPHAZ_GIVEN))
      alpha2 = -S12/emit2;
    else
      alpha2 = match->alpha;
    ratio = sqrt(emit2/emit1);
    
    match->R11 = R11 = ratio*beta2/sqrt(beta1*beta2);
    match->R12 = R12 = 0;
    match->R21 = R21 = ratio*(alpha1-alpha2)/sqrt(beta1*beta2);
    match->R22 = R22 = ratio*beta1/sqrt(beta1*beta2);
  }
  else {
    R11 = match->R11;
    R12 = match->R12;
    R21 = match->R21;
    R22 = match->R22;
  }

  for (i=0; i<rows; i++) {
    t0 = t[i];
    delta0 = p[i];
    t[i] = R11*t0 + R12*delta0;
    p[i] = R21*t0 + R22*delta0;
  }

  if (match->flags&BETAGAMMA_GIVEN) {
    double ratio;
    ratio = sqrt(pAve/match->betaGamma);
    for (i=0; i<rows; i++) {
      x[i] *= ratio;
      xp[i] *= ratio;
      y[i] *= ratio;
      yp[i] *= ratio;
      t[i] += tAve;
      p[i] = p[i]*pAve + match->betaGamma;
    }
  } else {
    for (i=0; i<rows; i++) {
      p[i] = (p[i]+1)*pAve;
      t[i] += tAve;
    }
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

