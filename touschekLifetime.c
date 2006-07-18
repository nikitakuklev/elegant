/*************************************************************************\
 * Copyright (c) 2006 The University of Chicago, as Operator of Argonne
 * National Laboratory.
 * This file is distributed subject to a Software License Agreement found
 * in the file LICENSE that is included with this distribution. 
 \*************************************************************************/

/* 
 * $Log: not supported by cvs2svn $

 * sdds program to return Touschek lifetime.
 * Using A. Piwinski's formula, DESY 98-179/ISSN 0418-9833.
 * Input files are elegant twiss file with radiation integrals
 * parameters, and a momentum aperture file.
 * Calculated parameters will go in an output sdds file.
 * 
 * Original author: Aimin Xiao, ANL.
 */

#include <stdio.h>
#include "mdb.h"
#include "scan.h"
#include "match_string.h"
#include "SDDS.h"
#include "constants.h"

static char *USAGE = "touschekLifetime <resultsFile>\n\
 -twiss=<twissFile> -aperture=<momentumApertureFile>\n\
 {-charge=<nC>|-particles=<value>} -coupling=<value>\n\
 {-RF=Voltage=<MV>,harmonic=<value>|-length=<mm>}\n\
 [-emitxInput=<value>] [-deltaInput=<value>] [-verbose=<value>]\n\n\
 Computes Touschek lifetime using Piwinski's method, DESY 98-179/ISSN 0418-09833\n\n\
 This is version 1, A. Xiao (ANL/APS).";

#define VERBOSE 0
#define CHARGE 1
#define PARTICLES 2
#define COUPLING 3
#define RF 4
#define LENGTH 5
#define EMITXINPUT 6
#define DELTAINPUT 7
#define TWISSFILE 8
#define APERFILE 9
#define N_OPTIONS 10

char *option[N_OPTIONS] = {
  "verbose",
  "charge",
  "particles",
  "coupling",
  "rf",
  "length",
  "emitxinput",
  "deltainput",
  "twiss",
  "aperture",
};

void TouschekLifeCalc();  
void FIntegral(double *tm, double *B1, double *B2, double *F, long index); 
double Fvalue (double t, double tm, double b1, double b2);
double linear_interpolation(double *y, double *t, long n, double t0, long *i);

/* global varibles */
long plane_orbit = 0;
double *s, *s2, *dpp, *dpm;
double *betax, *alphax, *etax, *etaxp;
double *betay, *alphay, *etay, *etayp;
double *tm, *B1, *B2, *F, *coeff;
double tLife;
long elements, elem2;
long *eOccur1;
double pCentral, sz, sigmap; 
double NP, emitx, emity; 
char **eName1, **eName2, **eType1;
#define NDIV 1000;

int main( int argc, char **argv)
{
  SCANNED_ARG *scanned;
  char *inputfile1, *inputfile2, *outputfile;
  SDDS_DATASET twissPage, aperPage, resultsPage;
  long verbosity;
  double etaymin, etaymax;
  long i;
  unsigned long dummyFlags;
  double emitxInput, sigmaDeltaInput, rfVoltage, rfHarmonic;
  double alphac, U0, circumference, EMeV;
  double coupling, emitx0, charge;

  /****************************************************\
   * read from command line                           *
   \****************************************************/

  SDDS_RegisterProgramName(argv[0]);
  argc  =  scanargs(&scanned, argc, argv);
  if (argc == 1)
    bomb(NULL, USAGE);

  inputfile1  =  NULL;
  inputfile2  =  NULL;  
  outputfile  =  NULL;
  verbosity = 0;
  NP = 0;
  charge = 0;
  coupling = 0;
  sz = 0;
  emitxInput = 0;
  sigmaDeltaInput = 0;
  rfVoltage = rfHarmonic = 0;
  
  for (i = 1; i<argc; i++) {
    if (scanned[i].arg_type == OPTION) {
      delete_chars(scanned[i].list[0], "_");
      switch(match_string(scanned[i].list[0], option, N_OPTIONS, UNIQUE_MATCH)) {
      case VERBOSE:
        if (scanned[i].n_items > 1 ) {
          get_long(&verbosity, scanned[i].list[1]);
        } else {
          verbosity=1;
        }
        break;
      case CHARGE:
        if (scanned[i].n_items != 2 )
          bomb("invalid -charge syntax/values", "-charge=<nC>");
        get_double(&charge, scanned[i].list[1]);
        break;
      case EMITXINPUT:
        if (scanned[i].n_items != 2 ) 
          bomb("invalid -exitxInput syntax/values", "-emitxInput=<value>");        
        get_double(&emitxInput, scanned[i].list[1]);
        break;
      case DELTAINPUT:
        if (scanned[i].n_items != 2 ) 
          bomb("invalid -deltaInput syntax/values", "-deltaInput=<value>");        
        get_double(&sigmaDeltaInput, scanned[i].list[1]);
        break;
      case LENGTH:
        if (scanned[i].n_items != 2 )
          bomb("invalid -length syntax/values", "-length=<mm>");        
        get_double(&sz, scanned[i].list[1]);
        sz /= 1000;   /* convert input length from mm to m */
        break;
      case COUPLING:
        if (scanned[i].n_items != 2 )
          bomb("invalid -coupling syntax/values", "-coupling=<value>");        
        get_double(&coupling, scanned[i].list[1]);
        break;
      case PARTICLES:
        if (scanned[i].n_items != 2 )
          bomb("invalid -particles syntax/values", "-particles=<value>");        
        get_double(&NP, scanned[i].list[1]);
        break;
      case RF:
        if (scanned[i].n_items<2)
          bomb("invalid -rf syntax", NULL);
        scanned[i].n_items--;
        if (!scanItemList(&dummyFlags, scanned[i].list+1, &scanned[i].n_items, 0,
                          "voltage", SDDS_DOUBLE, &rfVoltage, 1, 0,
                          "harmonic", SDDS_DOUBLE, &rfHarmonic, 1, 0,
                          NULL) ||
            rfVoltage<=0 || rfHarmonic<=0)
          bomb("invalid -rf syntax/values", "-rf=voltage=MV,harmonic=<value>");
        break;
      case TWISSFILE:
        if (scanned[i].n_items<2)
          bomb("invalid -twiss syntax", NULL);
        inputfile1  =  scanned[i].list[1];
        break;
      case APERFILE:
        if (scanned[i].n_items<2)
          bomb("invalid -twiss syntax", NULL);
        inputfile2  =  scanned[i].list[1];
        break;
      default:
        bomb("unknown option given.", NULL);  
      }
    }
    else {
      if (!outputfile) 
        outputfile =  scanned[i].list[0];
      else
        bomb("too many filenames given", NULL);
    }
  }
  if (charge && NP) {
    bomb("Options charge and particles cannot be both specified.",NULL);
  }
  if (!charge) 
    charge = NP * e_mks;
  if (!NP) {
    /* command line input value is in units of nC */
    charge /= 1e9; 
    NP = charge/ e_mks;
  }
  if (!coupling) 
    bomb("Coupling value not specified.",NULL);
  if (!sz && !rfVoltage) 
    bomb("Specify either the bunch length or the rf voltage.", NULL);
  
  /****************************************************\
   * Check input twissfile                            *
   \****************************************************/
  if (verbosity)
    fprintf( stdout, "Opening \"%s\" for checking presence of parameters.\n", inputfile1);
  if (!SDDS_InitializeInput(&twissPage, inputfile1))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  /* Check presence of first radiation integral */
  SDDS_ReadPage(&twissPage);
  switch(SDDS_CheckParameter(&twissPage, "I1", NULL, SDDS_DOUBLE, verbosity?stdout:NULL)) {
  case SDDS_CHECK_NONEXISTENT:
    if (verbosity)
      fprintf( stdout, "\tParameter I1 not found in input file.\n");
    exit(1);
    break;
  case SDDS_CHECK_WRONGTYPE:
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    exit(1);
    break;
  case SDDS_CHECK_OKAY:
    break;
  default:
    fprintf( stdout, "Unexpected result from SDDS_CheckParameter routine.\n");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    exit(1);
    break;
  }
  /****************************************************\
   * Check input aperturefile                         *
   \****************************************************/
  if (verbosity)
    fprintf( stdout, "Opening \"%s\" for checking presence of parameters.\n", inputfile2);
  if (!SDDS_InitializeInput(&aperPage, inputfile2))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  /* Check presence of momentum aperture */
  SDDS_ReadPage(&aperPage);
  switch(SDDS_CheckColumn(&aperPage, "deltaPositive", NULL, SDDS_DOUBLE, verbosity?stdout:NULL)) {
  case SDDS_CHECK_NONEXISTENT:
    if (verbosity)
      fprintf( stdout, "\tColumn deltaPositive is not found in input file.\n");
    exit(1);
    break;
  case SDDS_CHECK_WRONGTYPE:
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    exit(1);
    break;
  case SDDS_CHECK_OKAY:
    break;
  default:
    fprintf( stdout, "Unexpected result from SDDS_CheckColumn routine.\n");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    exit(1);
    break;
  }
  /****************************************************\
   * Check output file and write file layout          *
   \****************************************************/
  if (verbosity)
    fprintf( stdout, "Opening \"%s\" for writing...\n", outputfile);
  if (!SDDS_InitializeOutput(&resultsPage, SDDS_BINARY, 1, "Touschek lifetime calculation",
                             "Touschek lifetime calculation", outputfile))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);

  if (!SDDS_TransferParameterDefinition(&resultsPage, &twissPage, "pCentral", NULL) ||
      !SDDS_TransferParameterDefinition(&resultsPage, &twissPage, "Sdelta0", NULL) ||
      !SDDS_TransferParameterDefinition(&resultsPage, &twissPage, "ex0", NULL))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);

  if (0>SDDS_DefineParameter(&resultsPage, "coupling", NULL, NULL, 
                             "Coupling", NULL, SDDS_DOUBLE, NULL) ||
      0>SDDS_DefineParameter(&resultsPage, "emitx", "$ge$r$bx$n", "m", 
                             "Horizontal emittance with coupling", NULL, SDDS_DOUBLE, NULL) ||
      0>SDDS_DefineParameter(&resultsPage, "emity", "$ge$r$by$n", "m", 
                             "Vertical emittance with coupling", NULL, SDDS_DOUBLE, NULL) ||
      0>SDDS_DefineParameter(&resultsPage, "Particles", NULL, NULL, 
                             "Particles", NULL, SDDS_DOUBLE, NULL) ||
      0>SDDS_DefineParameter(&resultsPage, "Charge", NULL, "nC", 
                             "Charge", NULL, SDDS_DOUBLE, NULL) ||
      0>SDDS_DefineParameter(&resultsPage, "sigmaz", "$gs$r$bz$n", "m", 
                             "Bunch length", NULL, SDDS_DOUBLE, NULL) ||
      0>SDDS_DefineParameter(&resultsPage, "tLifetime", NULL, "s", 
                             "Touschek half lifetime", NULL, SDDS_DOUBLE, NULL))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  
  if (!SDDS_TransferColumnDefinition(&resultsPage, &twissPage, "s", NULL) ||
      !SDDS_TransferColumnDefinition(&resultsPage, &twissPage, "betax", NULL) ||
      !SDDS_TransferColumnDefinition(&resultsPage, &twissPage, "alphax", NULL) ||
      !SDDS_TransferColumnDefinition(&resultsPage, &twissPage, "etax", NULL) ||
      !SDDS_TransferColumnDefinition(&resultsPage, &twissPage, "etaxp", NULL) ||
      !SDDS_TransferColumnDefinition(&resultsPage, &twissPage, "betay", NULL) ||
      !SDDS_TransferColumnDefinition(&resultsPage, &twissPage, "alphay", NULL) ||
      !SDDS_TransferColumnDefinition(&resultsPage, &twissPage, "etay", NULL) ||
      !SDDS_TransferColumnDefinition(&resultsPage, &twissPage, "etayp", NULL) ||
      !SDDS_TransferColumnDefinition(&resultsPage, &twissPage, "ElementName", NULL) ||
      !SDDS_TransferColumnDefinition(&resultsPage, &twissPage, "ElementType", NULL) ||
      !SDDS_TransferColumnDefinition(&resultsPage, &twissPage, "ElementOccurence", NULL))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  if (0>SDDS_DefineColumn(&resultsPage, "tm", NULL, NULL, 
                          "Local momentum aperture", NULL, SDDS_DOUBLE, 0) ||
      0>SDDS_DefineColumn(&resultsPage, "B1", NULL, NULL, 
                          "Piwinski's parameter B1", NULL, SDDS_DOUBLE, 0) || 
      0>SDDS_DefineColumn(&resultsPage, "B2", NULL, NULL, 
                          "Piwinski's parameter B2", NULL, SDDS_DOUBLE, 0) || 
      0>SDDS_DefineColumn(&resultsPage, "c0", NULL, NULL, 
                          "1/T=c0*F", NULL, SDDS_DOUBLE, 0) || 
      0>SDDS_DefineColumn(&resultsPage, "F", NULL, NULL, 
                          "Piwinski's parameter F", NULL, SDDS_DOUBLE, 0))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);

  if (!SDDS_WriteLayout(&resultsPage) )
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);

  /****************************************************\
   * read from twiss file                             *
   \****************************************************/

  if (!SDDS_GetParameters(&twissPage,
                          "pCentral", &pCentral,
                          "ex0", &emitx0,
                          "Sdelta0", &sigmap,                            
                          "alphac", &alphac,
                          "U0", &U0,
                          NULL) )
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);

  elements = SDDS_CountRowsOfInterest(&twissPage);
  s = SDDS_GetColumnInDoubles(&twissPage, "s");
  elem2 = SDDS_CountRowsOfInterest(&aperPage);
  s2 = SDDS_GetColumnInDoubles(&aperPage, "s");
  if (s[0]<s2[0] || s[elements-1]>s2[elem2-1])
    bomb("aperture file s range does not cover twiss file s range", NULL);
  emitx = emitx0/ ( 1 + coupling);
  if (emitxInput) 
    emitx = emitxInput;
  emity = emitx * coupling;
  if (sigmaDeltaInput) 
    sigmap = sigmaDeltaInput;
  circumference = s[elements-1];
  EMeV = sqrt(sqr(pCentral) + 1) * me_mev;
  if (!sz) {
    /* compute length in m from rf voltage, energy spread, etc */
    sz = 
      circumference*sigmap*
        sqrt(alphac*EMeV/(PIx2*rfHarmonic*sqrt(sqr(rfVoltage)-sqr(U0))));
  }

  betax = SDDS_GetColumnInDoubles(&twissPage, "betax");
  betay = SDDS_GetColumnInDoubles(&twissPage, "betay");
  alphax = SDDS_GetColumnInDoubles(&twissPage, "alphax");
  alphay = SDDS_GetColumnInDoubles(&twissPage, "alphay");
  etax = SDDS_GetColumnInDoubles(&twissPage, "etax");
  etaxp = SDDS_GetColumnInDoubles(&twissPage, "etaxp");
  etay = SDDS_GetColumnInDoubles(&twissPage, "etay");
  etayp = SDDS_GetColumnInDoubles(&twissPage, "etayp");
  eName1 = SDDS_GetColumn(&twissPage, "ElementName");
  eType1 = SDDS_GetColumn(&twissPage, "ElementType");
  eOccur1 = SDDS_GetColumn(&twissPage, "ElementOccurence");

  eName2 = SDDS_GetColumn(&aperPage, "ElementName");
  dpp = SDDS_GetColumnInDoubles(&aperPage, "deltaPositive");
  dpm = SDDS_GetColumnInDoubles(&aperPage, "deltaNegative");
  
  /****************************************************\
   * calculate Touschek Lifetime                      *
   \****************************************************/
  find_min_max(&etaymin, &etaymax, etay, elements);
  if ((etaymax-etaymin) < 1e-6)
    plane_orbit=1; 
  
  if (!(tm = SDDS_Malloc(sizeof(*tm)*elements)) ||
      !(B1 = SDDS_Malloc(sizeof(*B1)*elements)) ||
      !(B2 = SDDS_Malloc(sizeof(*B2)*elements)) ||
      !(F = SDDS_Malloc(sizeof(*F)*elements)) ||
      !(coeff = SDDS_Malloc(sizeof(*coeff)*elements)))
    bomb("memory allocation failure (integration arrays)", NULL);
  
  TouschekLifeCalc(verbosity);
  
  /****************************************************\
   * Write output file                                *
   \****************************************************/
  if (0>SDDS_StartPage(&resultsPage, elements) ||
      !SDDS_SetParameters(&resultsPage, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,
                          "pCentral", pCentral, 
                          "Sdelta0", sigmap,
                          "ex0", emitx0,
                          "coupling", coupling, 
                          "emitx", emitx,
                          "emity", emity,
                          "Particles", NP,
                          "Charge", (1e9 * charge),
                          "tLifetime", tLife,
                          "sigmaz", sz, NULL) ||
      !SDDS_SetColumn(&resultsPage, SDDS_SET_BY_NAME, s, elements, "s") ||
      !SDDS_SetColumn(&resultsPage, SDDS_SET_BY_NAME, betax, elements, "betax") ||
      !SDDS_SetColumn(&resultsPage, SDDS_SET_BY_NAME, alphax, elements, "alphax") ||
      !SDDS_SetColumn(&resultsPage, SDDS_SET_BY_NAME, etax, elements, "etax") ||
      !SDDS_SetColumn(&resultsPage, SDDS_SET_BY_NAME, etaxp, elements, "etaxp") ||
      !SDDS_SetColumn(&resultsPage, SDDS_SET_BY_NAME, betay, elements, "betay") ||
      !SDDS_SetColumn(&resultsPage, SDDS_SET_BY_NAME, alphay, elements, "alphay") ||
      !SDDS_SetColumn(&resultsPage, SDDS_SET_BY_NAME, etay, elements, "etay") ||
      !SDDS_SetColumn(&resultsPage, SDDS_SET_BY_NAME, etayp, elements, "etayp") ||
      !SDDS_SetColumn(&resultsPage, SDDS_SET_BY_NAME, eName1, elements, "ElementName") ||
      !SDDS_SetColumn(&resultsPage, SDDS_SET_BY_NAME, eType1, elements, "ElementType") ||
      !SDDS_SetColumn(&resultsPage, SDDS_SET_BY_NAME, eOccur1, elements, "ElementOccurence") ||
      !SDDS_SetColumn(&resultsPage, SDDS_SET_BY_NAME, tm, elements, "tm") ||
      !SDDS_SetColumn(&resultsPage, SDDS_SET_BY_NAME, B1, elements, "B1") ||
      !SDDS_SetColumn(&resultsPage, SDDS_SET_BY_NAME, B2, elements, "B2") ||
      !SDDS_SetColumn(&resultsPage, SDDS_SET_BY_NAME, coeff, elements, "c0") ||
      !SDDS_SetColumn(&resultsPage, SDDS_SET_BY_NAME, F, elements, "F") ||

      !SDDS_WritePage(&resultsPage))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  if (SDDS_ReadPage(&twissPage)>0)
    fprintf( stdout, "The code doesn't support multi twiss pages.\n");
  if (SDDS_ReadPage(&aperPage)>0)
    fprintf( stdout, "The code doesn't support multi aperture pages.\n");
  
  if (!SDDS_Terminate(&twissPage) || 
      !SDDS_Terminate(&aperPage) ||
      !SDDS_Terminate(&resultsPage))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  
  return(0);
  
}

void TouschekLifeCalc(long verbosity) 
{
  long i, j;
  double sp2, sp4, beta2, betagamma2, sh2;
  double sx2, sxb2, dx2, dx_2;
  double sy2, syb2, dy2, dy_2;
  double a0, c0, c1, c2, c3;
  double pm, pp;
  
  sp2 = sqr(sigmap);
  sp4 = sqr(sp2);
  beta2 = sqr(pCentral)/(sqr(pCentral)+1);
  betagamma2 = sqr(pCentral);
  a0 = sqr(re_mks)*c_mks*NP/(4*sqrt(PI)*(sqr(pCentral)+1)*sz);
  
  i=j=0;
  for (i = 0; i<elements; i++) {
    pp = linear_interpolation(dpp, s2, elem2, s[i], &j); 
    pm = linear_interpolation(dpm, s2, elem2, s[i], &j); 

    tm[i] = fabs(pp)>fabs(pm)?fabs(pm):fabs(pp);
    tm[i] = beta2*tm[i]*tm[i];

    sxb2 = betax[i]*emitx;
    syb2 = betay[i]*emity;
    dx2 = ipow(etax[i], 2);
    sx2  = sxb2 + dx2*sp2;
    
    dx_2 = ipow(alphax[i]*etax[i]+betax[i]*etaxp[i], 2);
    c1 = sqr(betax[i])/(2*betagamma2*sxb2);
    c2 = sqr(betay[i])/(2*betagamma2*syb2);
    
    if (plane_orbit) {
      sh2 = 1/(1/sp2+(dx2+dx_2)/sxb2);
      B1[i] = c1*(1-sh2*dx_2/sxb2)+c2;
      c0 = sqrt(sh2)/(sigmap*betagamma2*emitx*emity);
      c3 = sx2*syb2;
      B2[i] = sqr(B1[i])-sqr(c0)*c3;   
      if (B2[i]<0) {
        fprintf(stdout, "B2^2<0 at \"%s\" occurence %ld", eName1[i], eOccur1[i]);
        exit(1);
      }
      B2[i]=sqrt(B2[i]);
    }
    
    if (!plane_orbit) {
      dy2 = ipow(etay[i], 2);
      sy2  = syb2 + dy2*sp2;
      dy_2 = ipow(alphay[i]*etay[i]+betay[i]*etayp[i], 2);      
      sh2 = 1/(1/sp2+(dx2+dx_2)/sxb2+(dy2+dy_2)/syb2);
      c0 = sqrt(sh2)/(sigmap*betagamma2*emitx*emity);
      c3 = sx2*sy2-sp4*dx2*dy2;
      B1[i] = c1*(1-sh2*dx_2/sxb2)+c2*(1-sh2*dy_2/syb2);
      B2[i] = sqr(B1[i])-sqr(c0)*c3;   	  
      if (B2[i]<0) {
        fprintf(stdout, "B2^2<0 at \"%s\" occurence %ld", eName1[i], eOccur1[i]);
        exit(1);
      }
      B2[i]=sqrt(B2[i]);   	  
    }
    
    coeff[i] = a0*c0;
    
    if (i==0) FIntegral(tm, B1, B2, F, i);
    if (i>0 && s[i]>s[i-1]) 
      FIntegral(tm, B1, B2, F, i);
    else
      F[i] = F[i-1];
  }
  
  tLife = 0;  
  for (i = 1; i<elements; i++) {
    if (s[i]>s[i-1]) {
      tLife += (s[i]-s[i-1])*(coeff[i]*F[i]+coeff[i-1]*F[i-1])/2;
    }
  }
  tLife /= s[elements-1];
  tLife = 1/tLife;
  return;
} 

void FIntegral(double *tm, double *B1, double *B2, double *F, long index) 
{
  double f0, f1, sum;
  double HPI, step;
  long converge=0;
  double t1, k0, k1;
  double tstart, km, b1, b2;
  
  HPI = PI/2;
  step = HPI/NDIV;
  tstart = tm[index];
  b1 = B1[index];
  b2 = B2[index];
  k0 = km = atan(sqrt(tstart));
  f0 = Fvalue(tstart, tstart, b1, b2);
  sum = 0;
  f1 = 0;
  
  while (!converge) {
    k1 = k0 + step;
    t1 = sqr(tan(k1));
    if (isnan(t1)) {
      fprintf(stdout, "Integration failed at Element# %ld", index);
      break;
    }
    
    f1 = Fvalue(t1, tstart, b1, b2);

    if (abs(f1*(HPI-k1))<1e-3*sum)
      converge = 1;
    
    sum +=(f1+f0)/2*step;
    k0 = k1;
    f0 = f1;    
  }
  F[index] = sum;
  return;
}

double Fvalue (double t, double tm, double b1, double b2)
{
  double c0, c1, c2, result;
  
  c0 = (sqr(2*t+1)*(t/tm/(1+t)-1)/t+t-sqrt(t*tm*(1+t))-(2+1/2/t)*log(t/tm/(1+t)))*sqrt(1+t);
  c1 = exp(-b1*t);
  c2 = dbesi0(b2*t);
  result = c0 * c1 * c2;
  /* If overflow/underflow use approximate equation for modified bessel function. */
  if (isnan(result) || result>MAXFLOAT) {
    result=c0*exp(b2*t-b1*t)/sqrt(PIx2*b2*t);
  } 
  return result;
}

/* Only linear_interpolate from i=0 to i=n-2. For i<0 using i=0. For i>n-2, using i=n-1 */
double linear_interpolation(double *y, double *t, long n, double t0, long *iStart)
{
  long i;
  i = *iStart;
  if (i<0)
    i = 0;
  if (i>=n-1)
    i = n-2;
  while (i<=n-2 && t0>t[i+1])
    i++;
  if (i==n-1) {
    *iStart = n-1;
    return(y[n-1]);
  }
  while (i>=0 && t0<t[i])
    i--;
  if (i==-1) {
    *iStart =  0;
    return(y[0]);
  }
  if (!(t0>=t[i] && t0<=t[i+1])) {
    fprintf(stdout, "failure to bracket point in time array: t0=%e, t[0] = %e, t[n-1] = %e\n",
            t0, t[0], t[n-1]);
    fflush(stdout);
    abort();
  }

  *iStart = i;

  /* this handles zero length elements */ 
  if (t[i+1]==t[i])
    return y[i];
  
  return( y[i] + (y[i+1]-y[i])/(t[i+1]-t[i])*(t0-t[i]) );
}
