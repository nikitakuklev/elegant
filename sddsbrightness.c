/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* program: sddsbrightness.c
 * purpose: take twiss output from elegant and compute undulator
 *          brightness curves
 *
 * Michael Borland, 2002
 *
 $Log: not supported by cvs2svn $
 Revision 1.4  2002/05/07 20:02:24  shang
 modified the computation of brightness by using convolution when the width of
 sinc() function is wider than that of energ spread (gaussian func)

 Revision 1.3  2002/04/24 14:26:38  borland
 Refined calculation of radiation size and divergence.  Added
 -noSpectralBroadening switch to allow turning off this part of the
 calculation.

 Revision 1.2  2002/04/22 21:01:01  borland
 Added factors of 2 for square radiation opening angle and square of
 radiation size, per Dejus.
 Refined energy-broadening effect by fitting a gaussian to the sinc function
 to get a width of xn=0.36 for sigma.

 Revision 1.1  2002/04/18 23:41:57  borland
 First version, with assistence of L. Emery.

 */
#include "mdb.h"
#include "scan.h"
#include "SDDS.h"

#define DEBUG 0

#define SET_PIPE 0
#define SET_HARMONICS 1
#define SET_KRANGE 2
#define SET_CURRENT 3
#define SET_TOTALLENGTH 4
#define SET_PERIODLENGTH 5
#define SET_EMITTANCERATIO 6
#define SET_COUPLING 7
#define SET_NOSPECTRALBROADENING 8
#define N_OPTIONS 9

char *option[N_OPTIONS] = {
  "pipe", "harmonics", "krange", "current", "totallength", "periodlength", 
  "emittanceratio", "coupling", "nospectralbroadening",
} ;

char *USAGE="sddsbrightness [-pipe=[input][,output]] [<twissFile>] [<SDDSoutputfile>]\n\
 -harmonics=<integer> -Krange=start=<value>,end=<value>,points=<integer>\n\
 -current=<Amps> -totalLength=<meters> -periodLength=<meters>\n\
 [-emittanceRatio=<value> | -coupling=<value>] [-noSpectralBroadening]\n\n\
harmonics        number of harmonics to compute\n\
Krange           range and number of undulator K parameter to evaluate\n\
current          beam current in amperes\n\
totalLength      total length of the undulator in meters\n\
periodLength     length of the undulator period in meters\n\
emittanceRatio   ratio of y emittance to x emittance.  x emittance is\n\
                 ex0 from input file. y emittance is ratio*ex0\n\
coupling         x emittance is ex0/(1+coupling), while y emittance is\n\
                 coupling*ex0/(1+coupling).\n\
noSpectralBroadening\n\
                 Turns off the default inclusion of spectral broadening in\n\
                 the calculation.  Gives an over-estimate of the brightness.\n\n\
Computes on-axis brightness for an undulator centered on the end of the beamline\n\
the Twiss parameters for which are in the input file.  You should generate the\n\
input file using elegant's twiss_output command with radiation_integrals=1 .\n\n\
Program by Michael Borland.  (This is version 2, April 2002.)\n";

long SetUpOutputFile(SDDS_DATASET *SDDSout, SDDS_DATASET *SDDSin, char *outputfile, long harmonics);
double ComputeBrightness(double period, long Nu, double K, long n,
                         double gamma, double ex0, double Sdelta0, 
                         double coupling, double emitRatio, 
                         double current, 
                         double betax, double alphax, double etax, double etaxp,
                         double betay, double alphay, double etay, double etayp,
                         double *lambda, double *Fn, 
                         short spectralBroadening);
long GetTwissValues(SDDS_DATASET *SDDSin, 
                    double *betax, double *alphax, double *etax, double *etaxp, 
                    double *betay, double *alphay, double *etay, double *etayp, 
                    double *ex0, double *Sdelta0, double *pCentral);
double computeFactorOfConvolution(long periods, long harmonic, double Sdelta0);
double convolutionFunc(double x);
double delta0,sincNu; /*two constants used in convolutionFunc() */

int main(int argc, char **argv)
{
  SDDS_DATASET SDDSin, SDDSout;
  char *inputfile, *outputfile;
  SCANNED_ARG *s_arg;
  unsigned long pipeFlags;
  double current, totalLength, periodLength, KStart, KEnd, coupling, emittanceRatio, dK;
  long KPoints, harmonics, tmpFileUsed, iK, i_arg, poles, readCode, h, ih, periods;
  unsigned long dummyFlags;
  double betax, alphax, betay, alphay, etax, etaxp, etay, etayp, lambda, energy;
  double pCentral, ex0, Bn, Fn, K, Sdelta0, conFactor=1.0; /*the factor from convolution of
                                                             sinc() and gaussian() */
  short spectralBroadening;
  
  SDDS_RegisterProgramName(argv[0]);
  argc = scanargs(&s_arg, argc, argv);
  if (argc<2) 
    bomb(NULL, USAGE);

  inputfile = outputfile = NULL;
  pipeFlags = 0;
  spectralBroadening = 1;
  current = totalLength = periodLength = 0;
  KStart = KEnd = coupling = emittanceRatio = 0;
  harmonics = KPoints = 0;
  
  for (i_arg=1; i_arg<argc; i_arg++) {
    if (s_arg[i_arg].arg_type==OPTION) {
      switch (match_string(s_arg[i_arg].list[0], option, N_OPTIONS, 0)) {
      case SET_PIPE:
        if (!processPipeOption(s_arg[i_arg].list+1, s_arg[i_arg].n_items-1, &pipeFlags))
          SDDS_Bomb("invalid -pipe syntax");
        break;
      case SET_HARMONICS:
        if (s_arg[i_arg].n_items!=2 ||
            !sscanf(s_arg[i_arg].list[1], "%ld", &harmonics) ||
            harmonics<1)
          SDDS_Bomb("invalid -harmonics value");
        break;
      case SET_KRANGE:
        s_arg[i_arg].n_items--;
        if (!scanItemList(&dummyFlags, s_arg[i_arg].list+1, &s_arg[i_arg].n_items, 0,
                          "start", SDDS_DOUBLE, &KStart, 1, 0,
                          "end", SDDS_DOUBLE, &KEnd, 1, 0,
                          "points", SDDS_LONG, &KPoints, 1, 0,
                          NULL) ||
            KStart>=KEnd || KStart<0 || KPoints<2)
          SDDS_Bomb("invalid -Krange syntax");
        break;
      case SET_CURRENT:
        if (s_arg[i_arg].n_items!=2 ||
            !sscanf(s_arg[i_arg].list[1], "%le", &current) ||
            current<=0)
          SDDS_Bomb("invalid -current value");
        break;
      case SET_TOTALLENGTH:
        if (s_arg[i_arg].n_items!=2 ||
            !sscanf(s_arg[i_arg].list[1], "%le", &totalLength) ||
            totalLength<=0)
          SDDS_Bomb("invalid -totalLength value");
        break;
      case SET_PERIODLENGTH:
        if (s_arg[i_arg].n_items!=2 ||
            !sscanf(s_arg[i_arg].list[1], "%le", &periodLength) ||
            periodLength<=0)
          SDDS_Bomb("invalid -periodLength value");
        break;
      case SET_EMITTANCERATIO:
        if (s_arg[i_arg].n_items!=2 ||
            !sscanf(s_arg[i_arg].list[1], "%le", &emittanceRatio) ||
            emittanceRatio<=0)
          SDDS_Bomb("invalid -emittanceRatio value");
        break;
      case SET_COUPLING:
        if (s_arg[i_arg].n_items!=2 ||
            !sscanf(s_arg[i_arg].list[1], "%le", &coupling) ||
            coupling<=0)
          SDDS_Bomb("invalid -coupling value");
        break;
      case SET_NOSPECTRALBROADENING:
        spectralBroadening = 0;
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

  if (coupling && emittanceRatio)
    SDDS_Bomb("give only one of -coupling or -emittanceRatio");
  if (!harmonics)
    SDDS_Bomb("you must specify the number of harmonics to compute");
  if (!KPoints)
    SDDS_Bomb("you must specify the range and number of K values");
  if (!current)
    SDDS_Bomb("you must specify the current");
  if (!totalLength)
    SDDS_Bomb("you must specify the total undulator length");
  if (!periodLength)
    SDDS_Bomb("you must specify the undulator period length");
  /*poles = totalLength/periodLength; rename poles as periods, which makes sense.*/
  periods=totalLength/periodLength;
  if (periods<1)
    SDDS_Bomb("period lengths is shorter than undulator length!");
  
  processFilenames("sddsbrightness", &inputfile, &outputfile, pipeFlags, 0, &tmpFileUsed);
  if (tmpFileUsed)
    SDDS_Bomb("can't overwrite input file");
  
  if (!SDDS_InitializeInput(&SDDSin, inputfile))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);

  /* check the input file for valid data */
  if (SDDS_CheckColumn(&SDDSin, "betax", "m", SDDS_ANY_FLOATING_TYPE, stderr)!=SDDS_CHECK_OK || 
      SDDS_CheckColumn(&SDDSin, "alphax", "", SDDS_ANY_FLOATING_TYPE, stderr)!=SDDS_CHECK_OK || 
      SDDS_CheckColumn(&SDDSin, "etax", "m", SDDS_ANY_FLOATING_TYPE, stderr)!=SDDS_CHECK_OK ||
      SDDS_CheckColumn(&SDDSin, "etaxp", "", SDDS_ANY_FLOATING_TYPE, stderr)!=SDDS_CHECK_OK ||
      SDDS_CheckColumn(&SDDSin, "betay", "m", SDDS_ANY_FLOATING_TYPE, stderr)!=SDDS_CHECK_OK || 
      SDDS_CheckColumn(&SDDSin, "alphay", "", SDDS_ANY_FLOATING_TYPE, stderr)!=SDDS_CHECK_OK ||
      SDDS_CheckColumn(&SDDSin, "etay", "m", SDDS_ANY_FLOATING_TYPE, stderr)!=SDDS_CHECK_OK ||
      SDDS_CheckColumn(&SDDSin, "etayp", "", SDDS_ANY_FLOATING_TYPE, stderr)!=SDDS_CHECK_OK ) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    SDDS_Bomb("something wrong with twiss parameter columns");
  }
  if (SDDS_CheckParameter(&SDDSin, "ex0", "$gp$rm", SDDS_ANY_FLOATING_TYPE, stderr)!=SDDS_CHECK_OK || 
      SDDS_CheckParameter(&SDDSin, "pCentral", "m$be$nc", SDDS_ANY_FLOATING_TYPE, stderr)!=SDDS_CHECK_OK ||
      SDDS_CheckParameter(&SDDSin, "Sdelta0", "", SDDS_ANY_FLOATING_TYPE, stderr)!=SDDS_CHECK_OK) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    SDDS_Bomb("something wrong with ex0, pCentral, or Sdelta0 parameters");
  }

  if (!SetUpOutputFile(&SDDSout, &SDDSin, outputfile, harmonics))
    SDDS_Bomb("problem setting up output file");

  dK = (KEnd-KStart)/(KPoints-1);
  while ((readCode=SDDS_ReadPage(&SDDSin))>0) {
    if (!GetTwissValues(&SDDSin, 
                        &betax, &alphax, &etax, &etaxp, 
                        &betay, &alphay, &etay, &etayp, 
                        &ex0, &Sdelta0, &pCentral))
      SDDS_Bomb("problem getting twiss parameters and other values from input file");
    if (!SDDS_StartPage(&SDDSout, KPoints))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    
    for (ih=0; ih<harmonics; ih++) {
      h = ih*2+1;
      /*convolution factor is same for all K s, it varies on periods,h and Sdelta0 only*/
      conFactor=1.0; 
      if (spectralBroadening && Sdelta0!=0 && 0.36/periods/h>Sdelta0) {
        conFactor=computeFactorOfConvolution(periods,h,Sdelta0);
      } 
      for (K=KStart, iK=0; iK<KPoints; iK++, K+=dK) {
        Bn = conFactor*ComputeBrightness(periodLength, periods, K, h,
                               pCentral, ex0, Sdelta0,
                               coupling, emittanceRatio, current,
                               betax, alphax, etax, etayp, 
                               betay, alphay, etay, etaxp, 
                               &lambda, &Fn, spectralBroadening);
        energy = 12.39/(lambda*1e10);
        if (!SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, iK,
                               ih*4+1, Bn,
                               ih*4+2, Fn,
                               ih*4+3, lambda, 
                               ih*4+4, energy, -1) ||
            (h==1 &&
             !SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, iK, 0, K, -1)))
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      }
    }
    if (!SDDS_WritePage(&SDDSout))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  
  if (!SDDS_Terminate(&SDDSin) || !SDDS_Terminate(&SDDSout))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  free_scanargs(&s_arg,argc);
  
  return 0;
}

long SetUpOutputFile(SDDS_DATASET *SDDSout, SDDS_DATASET *SDDSin, char *outputfile, long harmonics)
{
  long h;
  char buffer[1024];
  
  if (!SDDS_InitializeOutput(SDDSout, SDDS_BINARY, 1, NULL, NULL, outputfile) ||
      !SDDS_TransferAllParameterDefinitions(SDDSout, SDDSin, 0) ||
      !SDDS_DefineSimpleColumn(SDDSout, "K", NULL, SDDS_DOUBLE))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  for (h=1; h<2*harmonics; h+=2) {
    sprintf(buffer, "Brightness%ld", h);
    if (!SDDS_DefineSimpleColumn(SDDSout, buffer, "photons/s/mrad$a2$n/mm$a2$n0.1%BW", SDDS_DOUBLE))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    sprintf(buffer, "F%ld", h);
    if (!SDDS_DefineSimpleColumn(SDDSout, buffer, NULL, SDDS_DOUBLE))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    sprintf(buffer, "wavelength%ld", h);
    if (!SDDS_DefineSimpleColumn(SDDSout, buffer, "m", SDDS_DOUBLE))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    sprintf(buffer, "photonEnergy%ld", h);
    if (!SDDS_DefineSimpleColumn(SDDSout, buffer, "keV", SDDS_DOUBLE))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  if (!SDDS_WriteLayout(SDDSout))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  return 1;
}

double ComputeBrightness(double period, long Nu, double K, long n, 
                         double gamma, double ex0, double Sdelta0, 
                         double coupling, double emitRatio, 
                         double current, 
                         double betax, double alphax, double etax, double etaxp,
                         double betay, double alphay, double etay, double etayp,
                         double *lambdaOut, double *FnOut,
                         short spectralBroadening)
{
  double JArg, Nn;
  double Sx, Sy, Sxp, Syp, Srp2, Sr2;
  double ex, ey, lambda, Fn;
  double gammax, gammay, length, broadening;
  double sincWidth;
  
  JArg = n*K*K/(4+2*K*K);
  Fn = sqr(n*K/(1+K*K/2)*(jn((n+1)/2, JArg) - jn((n-1)/2, JArg)));
  
  if (coupling) {
    ex = ex0/(1+coupling);
    ey = coupling*ex;
  } else {
    ex = ex0;
    ey = ex0*emitRatio;
  }
  lambda = period/(2*n*sqr(gamma))*(1+sqr(K)/2);

  /* 0.001 is for 0.1% bandwidth */
  Nn = 0.5*PI/137.04*Nu*(current/e_mks)*0.001*Fn*(1+sqr(K)/2)/n;
  length = Nu*period;
  sincWidth=0.36/Nu/n;
  
  broadening = 1;
  if (spectralBroadening) {
    /* this factor includes the loss of peak brightness due to spectral
     * broadening from the energy spread. The factor of 2 in front of Sdelta0
     * is because wavelength ~ 1/gamma^2 .  The 0.36 is from a gaussian fit
     * to the sinc function, which gives sigma(omega)/omega = 0.36/(n*Nu)
     */
    broadening = sqrt( 1 + sqr(2*Sdelta0*n*Nu/0.36) );
  }

  /* radiation divergence and size, squared.
   * numerical factors are from "Radiation Properties of an Undulator",
   * SLS-Note 4/95, by W. Joho
   */
  Srp2 = sqr(0.58)*lambda/length;   /* Joho has *broadening here, which I think is wrong */
  Sr2  = sqr(0.12/0.58)*lambda*length;

  gammax = (1+sqr(alphax))/betax;
  gammay = (1+sqr(alphay))/betay;
  Sxp = sqrt(ex*gammax + sqr(Sdelta0*etaxp) + Srp2);
  Syp = sqrt(ey*gammay + sqr(Sdelta0*etayp) + Srp2);
  Sx = sqrt(ex*betax + sqr(Sdelta0*etax) + Sr2);
  Sy = sqrt(ey*betay + sqr(Sdelta0*etay) + Sr2);
  
  *lambdaOut = lambda;
  *FnOut = Fn;
  /* 1e-12 is to convert to 1/(mm^2*mrad^2) units */
  /* Joho does not put the broadening here, which I think is wrong */
  if (Sdelta0!=0 && sincWidth>Sdelta0) 
    broadening=1.0; /*calculate the factor by convolution instead */
  return Nn/(sqr(PIx2)*Sx*Sxp*Sy*Syp)/broadening*1e-12;
}

/*this function is to get the brightness by the convolution of sinc() and
  gaussian function (the energy spreed).
  since sinc() is symmetric, the peak position of the convolution is the same
  as sinc(), that is, at x=0, thus, the convolution is calculated only at
  one point: x=0, and returned as the brightness */
/* S(x)=(sin(N*pi*x)/(N*pi*x))^2
   g(x)=Cexp(-1/8 * (x/Sdelta0)^2) */
/* g(x)=1.0*10^-7 while x=+- 4*Sdelta0*log(1e-7/C) (the range of convoluation */
/* convoluation=integral of (g(x)*S(-x)dx, using gaussianQuadrature() to calculate it */
double computeFactorOfConvolution(long periods, long harmonic, double Sdelta0)
{
  double startX,C, conFactor=1.0;
  C=1/sqrt(PI*2)*0.5/Sdelta0; /*theorectical value for C=1/sqrt(2*pi)*1/(2*Sdelat0) */
  startX=-(4*Sdelta0*sqrt(fabs(log(1e-7/C)))); 
  sincNu=harmonic*periods*PI; /*sincNu and delta0 are two globals*/
  delta0=Sdelta0;
  gaussianQuadrature(convolutionFunc,startX,fabs(startX),1,0.01,&conFactor); 
  /*fprintf(stderr,"integrated result=%.9f\n",conFactor); */
  return conFactor;
}
/*integral function of gaussianQuadrature() */
double convolutionFunc(double x)
{
  double multResult,C,sx,gx;
  C=1/sqrt(PI*2)*0.5/delta0;
  if (x==0)
    sx=1;
  else
    sx=sqr(sin(sincNu*x)/(sincNu*x));
  gx=exp(-(x*x)/(8*delta0*delta0));
  multResult=C*sx*gx;
  return multResult;
}


long GetTwissValues(SDDS_DATASET *SDDSin, 
                    double *betax, double *alphax, double *etax, double *etaxp, 
                    double *betay, double *alphay, double *etay, double *etayp, 
                    double *ex0, double *Sdelta0, double *pCentral)
{
  double *data;
  long rows;
  
  if (!(rows=SDDS_RowCount(SDDSin)))
    return 0;

  if (!(data=SDDS_GetColumnInDoubles(SDDSin, "betax")))
    SDDS_Bomb("unable to get betax");
  *betax = data[rows-1];
  free(data);

  if (!(data=SDDS_GetColumnInDoubles(SDDSin, "alphax")))
    SDDS_Bomb("unable to get alphax");
  *alphax = data[rows-1];
  free(data);
  
  if (!(data=SDDS_GetColumnInDoubles(SDDSin, "etax")))
    SDDS_Bomb("unable to get etax");
  *etax = data[rows-1];
  free(data);

  if (!(data=SDDS_GetColumnInDoubles(SDDSin, "etaxp")))
    SDDS_Bomb("unable to get etax");
  *etaxp = data[rows-1];
  free(data);

  if (!(data=SDDS_GetColumnInDoubles(SDDSin, "betay")))
    SDDS_Bomb("unable to get betay");
  *betay = data[rows-1];
  free(data);

  if (!(data=SDDS_GetColumnInDoubles(SDDSin, "alphay")))
    SDDS_Bomb("unable to get alphay");
  *alphay = data[rows-1];
  free(data);

  if (!(data=SDDS_GetColumnInDoubles(SDDSin, "etay")))
    SDDS_Bomb("unable to get etay");
  *etay = data[rows-1];
  free(data);

  if (!(data=SDDS_GetColumnInDoubles(SDDSin, "etayp")))
    SDDS_Bomb("unable to get etay");
  *etayp = data[rows-1];
  free(data);

  if (!SDDS_GetParameterAsDouble(SDDSin, "ex0", ex0) ||
      !SDDS_GetParameterAsDouble(SDDSin, "Sdelta0", Sdelta0) ||
      !SDDS_GetParameterAsDouble(SDDSin, "pCentral", pCentral) )
    SDDS_Bomb("unable to get parameters from twiss file");
  return 1;
}

