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
 Revision 1.5  2002/08/14 20:23:48  soliday
 Added Open License

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
#include "sddsbrightness.h"

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
#define SET_METHOD 9
#define N_OPTIONS 10

char *option[N_OPTIONS] = {
  "pipe", "harmonics", "krange", "current", "totallength", "periodlength", 
  "emittanceratio", "coupling", "nospectralbroadening","method",
} ;

char *USAGE="sddsbrightness [-pipe=[input][,output]] [<twissFile>] [<SDDSoutputfile>]\n\
 -harmonics=<integer> -Krange=start=<value>,end=<value>,points=<integer>\n\
 -current=<Amps> -totalLength=<meters> -periodLength=<meters>\n\
 [-emittanceRatio=<value> | -coupling=<value>] [-noSpectralBroadening]\n\
 [-method=<string value>,device=<string value>,neks=<value>]] \n\n\
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
                 the calculation.  Gives an over-estimate of the brightness.\n\
method           choose method for calculating brightness \n\
                 method=borland  Michael Borland's approximation method. \n\
                 method=dejus    Non-zero emittance, \n\
                                 infinite-N +convolution (Dejus' approach) (default) \n\
                 method=walkerinfinite        Non-zero emittance, \n\
                                              infinite-N +convolution (Walker's approach) \n\
                 method=walkerfinite          Non-zero emittance; finite-N (Walker's) \n\
                 device=planar or helical      undulator type. \n\
                 neks=<value> number of points for peaking search. \n\n\
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

/*following functions are needed for calculating brightness using Dejus's method */
void FindPeak(double *E,double *spec,double *ep,double *sp,long n);
int Gauss_Convolve(double *E,double *spec,long *ns,double sigmaE);
/* note that sigmaE=Sdelta0 */
void Dejus_CalculateBrightness(double ENERGY, double current,long nE,
                               double period_mks, long nP, long device,
                               long ihMin,long ihMax,long ihStep,double sigmaE,
                               double gamma, double ex0,double coupling, double emitRatio, 
                               double betax, double alphax, double etax, double etaxp,
                               double betay, double alphay, double etay, double etayp,
                               long minNEKS, long maxNEKS, long neks,
                               double kMin,double kMax, long method,
                               double *sigmax,double *sigmay,double *sigmaxp,double *sigmayp,
                               double **K, double ***FnOut,
                               double ***Energy, double ***Brightness, double ***LamdarOut);
void ComputeBeamSize(double period, long Nu, double ex0, double Sdelta0, 
                     double coupling, double emitRatio,
                     double betax, double alphax, double etax, double etaxp,
                     double betay, double alphay, double etay, double etayp,
                     double *Sx, double *Sy, double *Sxp, double *Syp);


/*fortran subroutine*/
void usb_();

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
  long method, device,nE,ihMin,ihMax;
  double ikMin,ikMax;
  double ENERGY=7.0; /*storage ring energy Gev */
  double *KK,**FnOut,**Energy,**Brightness,**LamdarOut;
  long minNEKS,maxNEKS,neks;
  double sigmax,sigmay,sigmaxp,sigmayp;
  char *deviceOption;

  SDDS_RegisterProgramName(argv[0]);
  argc = scanargs(&s_arg, argc, argv);
  if (argc<2) 
    bomb(NULL, USAGE);

  KK=NULL;
  FnOut=Energy=Brightness=LamdarOut=NULL;
  deviceOption=NULL;

  inputfile = outputfile = NULL;
  pipeFlags = dummyFlags=0;
  spectralBroadening = 1;
  current = totalLength = periodLength = 0;
  KStart = KEnd = coupling = emittanceRatio = 0;
  harmonics = KPoints = 0;
  method=0;
  device=0;
  minNEKS=100;
  maxNEKS=500;
  neks=100;
  
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
      case SET_METHOD:
        if (s_arg[i_arg].n_items<2) 
          SDDS_Bomb("invalid -dejus syntax/values");
        switch (match_string(s_arg[i_arg].list[1], method_option, METHOD_OPTIONS, 0)) {
        case BORLAND:
          method=0;
          break;
        case DEJUS:
          method=1;
          break;
        case WALKER_INF:
          method=2;
          break;
        case WALKER_FIN:
          method=3;
          break;
        default:
          SDDS_Bomb("Invalid method given!");
          break;
        }
        s_arg[i_arg].n_items -=2;
        if (s_arg[i_arg].n_items>0 &&
            !scanItemList(&dummyFlags, s_arg[i_arg].list+2, &s_arg[i_arg].n_items, 0,
                          "device", SDDS_STRING, &deviceOption, 1, 0,
                          "neks", SDDS_LONG, &neks, 1, 0, NULL))
          SDDS_Bomb("invalid -method syntax/values");
        if (deviceOption) {
          switch (match_string(deviceOption, device_option, DEVICE_OPTIONS, 0)) {
          case PLANAR:
            device=0;
            break;
          case HELICAL:
            /*note that helical device has only one harmonics */
            harmonics=1;
            device=1;
            break;
          default:
            SDDS_Bomb("unknow device given!");
            break;
          }
          break;
        default:
          fprintf(stdout, "error: unknown switch: %s\n", s_arg[i_arg].list[0]);
          fflush(stdout);
          exit(1);
          break;
        }
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
    if (method) {
      ihMin=1;
      ihMax=2*(harmonics-1)+1;
      nE=KPoints;
      current=current*1.0e3; /*change unit from A to mA for dejus's method */
      Dejus_CalculateBrightness(ENERGY,current,nE,periodLength, periods,device,ihMin,ihMax,2,Sdelta0,
                                pCentral,ex0,coupling,emittanceRatio, 
                                betax,alphax,etax,etaxp,betay,alphay,etay,etayp,
                                minNEKS,maxNEKS,neks,KStart,KEnd,method,
                                &sigmax,&sigmay,&sigmaxp,&sigmayp,
                                &KK, &FnOut,&Energy,&Brightness,&LamdarOut);
      for (ih=0; ih<harmonics; ih++) {
        h = ih*2+1;
        if (h==1 && !SDDS_SetColumn(&SDDSout,SDDS_SET_BY_INDEX,KK,nE,0))
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        if (!SDDS_SetColumn(&SDDSout,SDDS_SET_BY_INDEX,Brightness[ih],nE,ih*4+1) ||
            !SDDS_SetColumn(&SDDSout,SDDS_SET_BY_INDEX,FnOut[ih],nE,ih*4+2) ||
            !SDDS_SetColumn(&SDDSout,SDDS_SET_BY_INDEX,LamdarOut[ih],nE,ih*4+3) ||
            !SDDS_SetColumn(&SDDSout,SDDS_SET_BY_INDEX,Energy[ih],nE,ih*4+4))
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      }
      if (!SDDS_SetParameters(&SDDSout,SDDS_BY_NAME|SDDS_PASS_BY_VALUE,"current",current,
                             "EnergySpread",Sdelta0,"sigmax",sigmax,"sigmay",sigmay,
                             "sigmaxprime",sigmaxp,"sigmayprime",sigmayp,
                             "period",periodLength,"numberOfPeriods",periods,NULL))
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    } else {    
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
    }
    
    if (!SDDS_WritePage(&SDDSout))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  
  if (!SDDS_Terminate(&SDDSin) || !SDDS_Terminate(&SDDSout))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  free_scanargs(&s_arg,argc);
  if (method) {
    for (ih=0;ih<harmonics;ih++) {
      free(Energy[ih]);
      free(FnOut[ih]);
      free(Brightness[ih]);
      free(LamdarOut[ih]);
    }
    free(Energy);
    free(FnOut);
    free(Brightness);
    free(LamdarOut);
    free(KK);
  }
  return 0;
}

long SetUpOutputFile(SDDS_DATASET *SDDSout, SDDS_DATASET *SDDSin, char *outputfile, long harmonics)
{
  long h;
  char buffer[1024];
  
  if (!SDDS_InitializeOutput(SDDSout, SDDS_BINARY, 1, NULL, NULL, outputfile) ||
      !SDDS_DefineSimpleColumn(SDDSout, "K", NULL, SDDS_DOUBLE))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  if (SDDS_DefineParameter(SDDSout,"current",NULL, "mA", NULL,NULL,SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineParameter(SDDSout,"EnergySpread",NULL, NULL, NULL,NULL,SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineParameter(SDDSout,"sigmax",NULL, "mm", NULL,NULL,SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineParameter(SDDSout,"sigmay",NULL, "mm", NULL,NULL,SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineParameter(SDDSout,"sigmaxprime",NULL,"mrad", NULL,NULL,SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineParameter(SDDSout,"sigmayprime",NULL,"mrad", NULL,NULL,SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineParameter(SDDSout,"period",NULL, "m", NULL,NULL,SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineParameter(SDDSout,"numberOfPeriods",NULL, NULL, NULL,NULL,SDDS_LONG, 0)<0)
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

void FindPeak(double *E,double *spec,double *ep,double *sp,long n)
{
  long i;
  
  *sp=spec[0];
  *ep=E[0];
  for (i=1; i<n; i++) {
    if (*sp<spec[i]) {
      *sp=spec[i];
      *ep=E[i];
    }
  }
}

int Gauss_Convolve(double *E,double *spec,long *ns,double sigmaE) 
{
  long nSigma=3,nppSigma=6,e_size=10000, ne1,ne2,ns1,np;
  
  int i,j;
  double ep,sp,sigp,de,sum, *gs,x, *spec2;
	
  ns1=*ns;
  gs=spec2=NULL;
  
  if (!E || !spec) {
    fprintf(stderr,"No energy or spectra points!\n");
    return 1;
  }
  FindPeak(E,spec,&ep,&sp,ns1);
  /*generate Gaussian with correct sigma in units of x-axis */
  de=E[1]-E[0];
  sigp=2.0*sigmaE*ep/de; /*sigma in x-axis units */
  
  if (sigp < (nppSigma-1)) {
    fprintf(stderr,"too few data points for Gaussian convolution\n");
    return 1; 
  }
  np=(2*nSigma)*sigp+1;
  if (np%2==0) np=np+1; /* make odd */
  gs=(double*)calloc(np,sizeof(*gs));
  spec2=(double*)calloc(ns1,sizeof(*spec2));
  sum=0.0;
  for (i=0;i<np;i++) {
    x=i*1.0-0.5*(np-1);
    gs[i]=exp(-x*x/2.0/(sigp*sigp));
  /*  fprintf(stderr,"x=%e, gs=%e\n",x,gs[i]); */
    sum=sum+gs[i];
  }
  /*fprintf(stderr,"sigp=%e,nSigma=%d\n",sigp,nSigma); */
  
  
  /*make convolution */
  ne1=np/2;
  ne2=ns1-ne1-1;
  if (ne2<0) {
    fprintf(stderr,"Error: Check the number of peak search points\n");
    return 1;
  }
  for (i=ne1;i<=ne2;i++) {
    spec2[i]=0.0;
    for (j=0;j<np;j++)
      spec2[i]=spec2[i]+spec[i+ne1-j]*gs[j];
  }
 /* fprintf(stderr,"np=%d, sum=%e, ne1=%d, ne2=%d\n",np,sum,ne1,ne2); */
  /*retun in original array and make adjustment of array sizes */
  *ns=ne2-ne1+1;
  for (i=ne1;i<=ne2;i++) {
    E[i-ne1]=E[i];
    spec[i-ne1]=spec2[i]/sum;
  }
  free(spec2);
  free(gs);
  return 0;
}

/*compute the beam size from twiss parameters
  output: Sx --- sigmaX
          Sy --- sigmaY
          Sxp --- sigmaX prime
          Syp --- sigmaY prime
*/
void ComputeBeamSize(double period, long Nu, double ex0, double Sdelta0, 
                         double coupling, double emitRatio,
                         double betax, double alphax, double etax, double etaxp,
                         double betay, double alphay, double etay, double etayp,
                         double *Sx, double *Sy, double *Sxp, double *Syp)
{
  double ex, ey;
  double gammax, gammay, length;
  
  if (coupling) {
    ex = ex0/(1+coupling);
    ey = coupling*ex;
  } else {
    ex = ex0;
    ey = ex0*emitRatio;
  }
  length = Nu*period;
  
  gammax = (1+sqr(alphax))/betax;
  gammay = (1+sqr(alphay))/betay;
  *Sxp = sqrt(ex*gammax + sqr(Sdelta0*etaxp))*1.0e3;
  *Syp = sqrt(ey*gammay + sqr(Sdelta0*etayp))*1.0e3;
  *Sx = sqrt(ex*betax + sqr(Sdelta0*etax))*1.0e3;
  *Sy = sqrt(ey*betay + sqr(Sdelta0*etay))*1.0e3;
}

/*caculate brightness by calling fortran subroutine usb(), written by Dejus */
/* input: nE number of points for calculation 
          ENERGY         Storage ring energy              (GeV)
          current        Storage ring current             (mA)
          sigmaE	 Energy spread (sigma(E)/E)
          

  output:
          sigmax:        RMS beam size (horizontal)       (mm
          sigmay:        RMS beam size (vertical)         (mm)
          sigmaxp:       RMS beam divergence (horizontal) (mrad)
          sigmayp:       RMS beam divergence (vertical)   (mrad)
         double *K,      array of K values
         double **FnOut  array of F factor fb[H][nE], H is the number of harmonics 
         double **Emergy  array of energy, eb[H][nE] in [kev]
         double **Brightness  array of brightness, sb[H][nE]
         double **LamdarOut array of lamda, LamdarOut[H][nE]
*/  

void Dejus_CalculateBrightness(double ENERGY, double current,long nE,
                               double period_mks, long nP, long device,
                               long ihMin,long ihMax,long ihStep,double sigmaE,
                               double gamma, double ex0,double coupling, double emitRatio, 
                               double betax, double alphax, double etax, double etaxp,
                               double betay, double alphay, double etay, double etayp,
                               long minNEKS, long maxNEKS, long neks,
                               double kMin,double kMax, long method,
                               double *sigmax,double *sigmay,double *sigmaxp,double *sigmayp,
                               double **K, double ***FnOut,
                               double ***Energy, double ***Brightness, double ***LamdarOut)
{
  double gamma1,lamdar,reducedE,kx,ky,eMin,eMax,ekMin,ekMax,ep,sp,dep1,dep2,fc,fc2,de,smax;
  long ih,ik,i,j,ie,ir,je,errorFlag=0;
  long nSigma=3,nppSigma=6,nek,ns,exitLoop=0;
  double JArg,sigmaEE,gk,dek,ek;
  double *tmpE,*tmpSpec,**ei,*ptot,*pd,*kyb,**eb,**sb;
  double e,k;
  double sigmaX,sigmaX1,sigmaY,sigmaY1,period;
  
  tmpE=tmpSpec=pd=ptot=NULL;
  ei=eb=sb=NULL;
  period=period_mks*1.0e2; /*use cm as units */
  if (*K) free(*K);
  if (*FnOut) free(*FnOut);
  if (*Energy) free(*Energy);
  if (*Brightness) free(*Brightness);
  if (*LamdarOut) free(*LamdarOut);
  *K=NULL;
  *FnOut=*Energy=*Brightness=*LamdarOut=NULL;
  
  if (neks<=nE) neks=nE+50;
  
  tmpE=(double*)calloc(MAXIMUM_E,sizeof(*tmpE));
  tmpSpec=(double*)calloc(MAXIMUM_E,sizeof(*tmpSpec));
  kyb=(double*)calloc(neks+100,sizeof(*kyb));
  ptot=(double*)calloc(neks+100,sizeof(*ptot));
  pd=(double*)calloc(neks+100,sizeof(*pd));
 
  ei=(double**)malloc(sizeof(*ei)*(neks+100));
  eb=(double**)malloc(sizeof(*eb)*(neks+100));
  sb=(double**)malloc(sizeof(*sb)*(neks+100));
  
  for (i=0;i<(neks+100);i++) {
    ei[i]=(double*)calloc(MAXIMUM_H,sizeof(**ei));
    eb[i]=(double*)calloc(MAXIMUM_H,sizeof(**eb));
    sb[i]=(double*)calloc(MAXIMUM_H,sizeof(**sb));
  }
  
  gamma1=ENERGY/me_mev*1.0E3;
  lamdar=period*1.0E8/(2.0*gamma1*gamma1); /*reduced wavelength A */
  reducedE=c_evang/lamdar; /*reduced energy ev */
  kx=0.0;
  eMin=reducedE/(1+kMax*kMax/2.0);
  eMax=reducedE/(1+kMin*kMin/2.0);
  
  if (device==HELICAL) {
    kMin=kMin/sqrt(2.0);
    kMax=kMax/sqrt(2.0);
  } 
  /*determin peak shifts for first and second harmonics at kMin */
  ky=kMin;
  nek=maxNEKS;
  /*compute electron beam size */
  ComputeBeamSize(period_mks,nP, ex0,sigmaE,coupling,emitRatio,
                  betax,alphax,etax,etaxp,
                  betay, alphay,etay,etayp,&sigmaX, &sigmaY, &sigmaX1, &sigmaY1);
  *sigmax=sigmaX;
  *sigmay=sigmaY;
  *sigmaxp=sigmaX1;
  *sigmayp=sigmaY1;
  
  for (i=1;i<3;i++) {
    if (i==1) {
      ekMin=0.95*i*eMax;
      ekMax=1.01*i*eMax;
    } else {
      ekMin=0.820*i*eMax;
      ekMax=1.002*i*eMax;
    }
    /*call fortran subroutin usb */
    usb_(&ENERGY,&current,&sigmaX,&sigmaY,&sigmaX1,&sigmaY1,&period,&nP,&kx,&ky,
        &ekMin,&ekMax,&nek,&method,tmpE,tmpSpec,&ns,&errorFlag);
    if (errorFlag) {
      fprintf(stderr,"error occurred in calling fortran subroutine usb\n");
      exit(1);
    }
    /*find the peak */
    FindPeak(tmpE,tmpSpec,&ep,&sp,ns);
    if (i==1)
      dep1=eMax*i-ep;
    else
      dep2=eMax*i-ep;
  }
  /*Main loop over harmonics and K-values */
  ih=0;
  de=(eMax-eMin)/(nE-1);
  for (i=ihMin;i<=ihMax;i=i+ihStep) {
    /*Try Kmin initially to find fc for each harmonic  
      Omit Brilliances < SPECLIM (will be set to 0.0) and peak positions will be
      calculated from the zero emittance formula. */
    smax=0;
    je=nE;
    nek=neks;
    do {
      je--;
      ek=eMin+je*de;
      ky=sqrt(2.0*(reducedE/ek-1.0));
      if (device==HELICAL) {
        ky=ky/sqrt(2.0);
        kx=ky;
      }
      if (ih==0) {
        kyb[je]=ky;
        ptot[je]=ptot_fac*nP*(kx*kx+ky*ky)*(ENERGY*ENERGY)*current*1.0e-3/(period*1.0e-2);
        if (device==HELICAL)
          gk=HELICAK(ky);
        else
          gk=PLANARK(ky);
        pd[je]=pd_fac*nP*ky*gk*pow(ENERGY,4)*current*1.0e-3/(period*1.0e-2);
      }
      if (i%2) {
        /*odd harmonics */
        ekMin=i*ek-i*dep1;
        ekMax=i*ek+i*dep1/2.0;
        if (i==1) ekMin=ekMin-dep1;
        if (i> (ek/dep1)) {
          fprintf(stderr,"Warning: overlapping range for initial peak search for harmonic %d\n",i);
          exitLoop=1;
          break;
        }
      } else {
        /*even harmonics */
        ekMin=i*ek-4.0*dep2;
        ekMax=i*ek;
      }
      /*call usb */
      usb_(&ENERGY,&current,&sigmaX,&sigmaY,&sigmaX1,&sigmaY1,&period,&nP,&kx,&ky,&ekMin,&ekMax,&nek,
          &method,tmpE,tmpSpec,&ns,&errorFlag);
      if (errorFlag) {        
        fprintf(stderr,"error occurred in calling fortran subroutine usb\n");
        exit(1);
      }      
      smax=0;
      for (j=0;j<ns;j++)
        if (tmpSpec[j]>smax) smax=tmpSpec[j];
      ei[je][ih]=i*ek;
      eb[je][ih]=i*ek;
      sb[je][ih]=0.0;
      
    } while (smax<minB && je>=0);
    if (exitLoop)
      break;
    if (smax < minB ) {
      fprintf(stderr,"Warning, Harmonic intensity too small, for harmonic number %d\n",i);
      break;
    }
    FindPeak(tmpE,tmpSpec,&ep,&sp,ns);
    /*define fc */
    fc=0.985*ep/ek/i;
    if (i > 1/(1-fc) ) {
      fprintf(stderr,"Warning: overlapping range for peak search for harmonics %d\n",i);
      break;
    }
    je=nE-1;
    for (j=0;j<=je;j++) {
      ek=eMin+j*de;
      ky=sqrt(2.0*(reducedE/ek-1.0));
      if (device==HELICAL) {
        ky=ky/sqrt(2.0);
        kx=ky;
      }
      if (ih==0) {
        kyb[j]=ky;
        ptot[j]=ptot_fac*nP*(kx*kx+ky*ky)*(ENERGY*ENERGY)*current*1.0e-3/(period*1.0e-2);
        if (device==HELICAL)
          gk=HELICAK(ky);
        else
          gk=PLANARK(ky);
        pd[j]=pd_fac*nP*ky*gk*pow(ENERGY,4)*current*1.0e-3/(period*1.0e-2);
      }
      
      if (i%2) 
        fc2=1.002;
      else
        fc2=1.000;
      ekMin=fc*i*ek;
      ekMax=fc2*i*ek;
      /* adjust ekmin, ekmax, and number of points if beam energy spread applied */
      if (sigmaE > 0) {
        nek=neks;
        dek=(ekMax-ekMin)/nek;
        sigmaEE=2.0*sigmaE*i*ek; /*  estimated width (eV) */
        ekMin=ekMin-sigmaEE*nSigma; /* adjust for convolution */
        ekMax=ekMax+sigmaEE*nSigma;  /* adjust for convolution */
        if (sigmaEE/nppSigma < dek) dek=sigmaEE/nppSigma;
        nek=(ekMax-ekMin)/dek+1; /*number of points */
        if (nek>MAXIMUM_E) {
          fprintf(stderr,"Energy points out of boundary (constrainted by FORTRAN usb subroutine).\n");
          exit(1);
        }
      }
      /* get undulator on-axis brilliance for given K value in energy range ekmin to ekmax
         returns energy array e (eV) and spec (ph/s/mrad^2/mm^2/0.1%bw) */
      usb_(&ENERGY,&current,&sigmaX,&sigmaY,&sigmaX1,&sigmaY1,&period,
          &nP,&kx,&ky,&ekMin,&ekMax,&nek,&method,tmpE,tmpSpec,&ns,&errorFlag);
      if (errorFlag) {
        fprintf(stderr,"error occurred in calling fortran subroutine usb\n");
        exit(1);
      } 
      FindPeak(tmpE,tmpSpec,&ep,&sp,ns);
     /* fprintf(stderr,"j=%d,points %d,maxE %e,minE %e,peakE %e,peakB %e\n",j,nek,ekMax,ekMin,ep,sp); */
      if (sigmaE>0) {
        /* gauss convolve */
        if (Gauss_Convolve(tmpE,tmpSpec,&ns,sigmaE))
          exit(1);
      }
      FindPeak(tmpE,tmpSpec,&ep,&sp,ns);
     /* fprintf(stderr,"after gauss convolve, peakE %e, peakB %e\n",ep,sp); */
      
      ei[j][ih]=i*ek;
      eb[j][ih]=ep;
      sb[j][ih]=sp;
    }
   /* fprintf(stderr,"Harmonics %d completed.\n",i); */
    ih++;
  } /*end for harmonics loop */
  
  /*output the result */
  *K=(double*)calloc(nE,sizeof(**K));
  *FnOut=(double**)malloc(sizeof(**FnOut)*ih);
  *Energy=(double**)malloc(sizeof(**Energy)*ih);
  *Brightness=(double**)malloc(sizeof(**Brightness)*ih);
  *LamdarOut=(double**)malloc(sizeof(**LamdarOut)*ih);
  for (i=0;i<ih;i++) {
    (*FnOut)[i]=calloc(nE,sizeof(***FnOut));
    (*Energy)[i]=calloc(nE,sizeof(***Energy));
    (*Brightness)[i]=calloc(nE,sizeof(***Brightness));
    (*LamdarOut)[i]=calloc(nE,sizeof(***LamdarOut));
  }
  for (j=0;j<nE;j++) {
    if (device==HELICAL)
      (*K)[j]=kyb[nE-1-j]*sqrt(2);
    else
      (*K)[j]=kyb[nE-1-j];
  }
  
  ih=0;
  for (i=ihMin;i<=ihMax;i=i+ihStep) {
    for (j=0;j<nE;j++) {
      k=(*K)[j];
      JArg = i*k*k/(4+2*k*k);
      (*FnOut)[ih][j]=pow(i*k/(1+k*k/2)*(jn((i+1)/2, JArg) - jn((i-1)/2, JArg)),2);
      (*Energy)[ih][j]=eb[nE-1-j][ih]*1.0e-3;
      (*Brightness)[ih][j]=sb[nE-1-j][ih];
      e=(*Energy)[ih][j];
      (*LamdarOut)[ih][j]=12.39/(e*1.0e10);
    }
    ih++;
  }
  free(tmpE);
  free(tmpSpec);
  free(kyb);
  free(pd);
  free(ptot);
  
  for (i=0;i<(neks+100);i++) {
    free(ei[i]);
    free(eb[i]);
    free(sb[i]);
  }
  free(ei);
  free(eb);
  free(sb);
  return;
}
