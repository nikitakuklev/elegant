/* $Log: not supported by cvs2svn $
 * Revision 1.1  1999/01/25 10:51:05  emery
 * First installation of ibsEmittance.c zibs.c
 * for intra-beam scattering rate calculation
 * and equilibrium emittance calculation.
 *
 * sdds program to return value of emittance
 * with intra-beam scattering included.
 * Input file is an elegant twiss file with
 * radiation integrals parameters.
 * Commmand line argument is energy in MeV.
 * Calculated parameters will go in an output
 * sdds file.
*/
#include <stdio.h>
#include "mdb.h"
#include "scan.h"
#include "match_string.h"
#include "SDDS.h"
#include "constants.h"
static char *USAGE = "ibsEmittance <twissFile> <resultsFile>\n\
 {-charge=<nC>|-particles=<value>} -coupling=<value>\n\
 [-emitxInput=<value>] [-deltaInput=<value>] \n\
 [-growthRatesOnly] \n\
 [-superperiods=<value>]\n\
 {-RF=Voltage=<MV>,harmonic=<value>|-length=<mm>}\n\
 [-energyRange=initial=<value>,final=<value>,delta=<value>,file=<string>]";

#define SET_ENERGY 0
#define VERBOSE 1
#define CHARGE 2
#define PARTICLES 3
#define COUPLING 4
#define RF 5
#define LENGTH 6
#define SUPERPERIOD 7
#define METHOD 8
#define EMITXINPUT 9
#define DELTAINPUT 10
#define GROWTHRATESONLY 11
#define N_OPTIONS 12
char *option[N_OPTIONS] = {
  "energyRange",
  "verbose",
  "charge",
  "particles",
  "coupling",
  "rf",
  "length",
  "superperiod",
  "method",
  "emitxInput",
  "deltaInput",
  "growthRatesOnly"
  };

#include "zibs.h"
double IBSequations(double *x, long *invalid);
void IBSsimplexReport(double ymin, double *xmin, long pass, long evals, long dims);

/* global variables */
double *s, *betax, *alphax, *betay, *alphay, *etax, *etaxp;

int main( int argc, char **argv)
{
  SCANNED_ARG *scanned;
  char *inputfile, *outputfile;
  SDDS_DATASET twissPage, resultsPage;
  double energy, particles, charge, length;
  long verbosity, i, j, k, elements, superperiods, growthRatesOnly;
  double pCentral, I1, I2, I3, I4, I5, taux, taudelta;
  double EMeV;
  double emitx0, emitx, emitxInput, emityInput, emity, coupling, sigmaz0, sigmaz;
  double sigmaDelta0, sigmaDelta, sigmaDeltaInput, xGrowthRate, yGrowthRate, zGrowthRate;
  double xGrowthRateOld, yGrowthRateOld, zGrowthRateOld;
  double Gx, Gy, Gz;
  double a, b, c, emitxOld, emityOld, sigmaDeltaOld, emittanceTolerance = 1e-5;
  double d, e, f, dfde, dfdd, dgde, dgdd, determinant, func1, func2;
  double x, dfdx, dgdx;
  long method, converged, loopLimit = 50, totalSteps;
/* used in simplex minimization */
  double yReturn, *xGuess, *dxGuess, *xLowerLimit, *xUpperLimit;
  short *disable;
  long dimensions = 14, maxEvaluations = 500, maxPasses = 2;
  double target = 0.0, tolerance = 0.01;
  
  SDDS_RegisterProgramName(argv[0]);
  argc  =  scanargs(&scanned, argc, argv);
  if (argc == 1)
    bomb(NULL, USAGE);

  inputfile  =  NULL;
  outputfile  =  NULL;
/*  energy = 0; */
  verbosity = 0;
  particles = 0;
  charge = 0;
  coupling = 0;
  length = 0;
  superperiods=1;
  method = 0;
  emitxInput = 0;
  sigmaDeltaInput = 0;
  growthRatesOnly = 0;
  for (i = 1; i<argc; i++) {
    if (scanned[i].arg_type == OPTION) {
      delete_chars(scanned[i].list[0], "_");
      switch(match_string(scanned[i].list[0], option, N_OPTIONS, UNIQUE_MATCH)) {
      case VERBOSE:
        if(scanned[i].n_items > 1 ) {
          get_long(&verbosity, scanned[i].list[1]);
        } else {
          verbosity=1;
        }
        break;
      case CHARGE:
        get_double(&charge, scanned[i].list[1]);
        break;
      case EMITXINPUT:
        get_double(&emitxInput, scanned[i].list[1]);
        break;
      case DELTAINPUT:
        get_double(&sigmaDeltaInput, scanned[i].list[1]);
        break;
      case LENGTH:
        get_double(&length, scanned[i].list[1]);
        length /= 1000; /* convert input length from mm to m */
        break;
      case COUPLING:
        get_double(&coupling, scanned[i].list[1]);
        break;
      case PARTICLES:
        get_double(&particles, scanned[i].list[1]);
        break;
      case SUPERPERIOD:
        get_long(&superperiods, scanned[i].list[1]);
        break;
      case METHOD:
        get_long(&method, scanned[i].list[1]);
        break;
      case GROWTHRATESONLY:
        growthRatesOnly = 1;
        break;
      case RF:
      case SET_ENERGY:
        bomb("Option %s not implemented.\n", scanned[i].list[0]);
        break;
      default:
        bomb("unknown option given.", NULL);  
      }
    }
    else {
      if (!inputfile)
        inputfile  =  scanned[i].list[0];
      else if (!outputfile) 
        outputfile =  scanned[i].list[0];
      else
        bomb("too many filenames given", NULL);
    }
  }
  if (charge && particles) {
    bomb("Options charge and particles cannot be both specified.",NULL);
  }
  if (!charge) 
    charge = particles * e_mks;
  if (!particles) {
    /* command line input value is in units of nC */
    charge /= 1e9; 
    particles = charge/ e_mks;
  }
  if (!coupling) 
    bomb("Coupling value not specified.",NULL);
  if (!length) 
    bomb("Bunch length value not specified.",NULL);

  /***************************************************\
   * get parameter information from first input file  *
   \***************************************************/
  if (verbosity)
    fprintf( stderr, "Opening \"%s\" for checking presence of parameters.\n", inputfile);
  if (!SDDS_InitializeInput(&twissPage, inputfile))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  /* read first page of input file to get parameters 
     I1 I2 I3 I4 I5.
     Check presence of first radiation integral.
     */
  SDDS_ReadPage(&twissPage);
  /* parameter Type */
  switch(SDDS_CheckParameter(&twissPage, "I1", NULL, SDDS_DOUBLE, verbosity?stderr:NULL)) {
  case SDDS_CHECK_NONEXISTENT:
    if (verbosity)
      fprintf( stderr, "\tParameter I1 not found in input file.\n");
    break;
  case SDDS_CHECK_WRONGTYPE:
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    exit(1);
    break;
  case SDDS_CHECK_OKAY:
    break;
  default:
    fprintf( stderr, "Unexpected result from SDDS_CheckParameter routine while checking parameter Type.\n");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    exit(1);
    break;
  }
  if (verbosity)
    fprintf( stderr, "Opening \"%s\" for writing...\n", outputfile);
  if (!SDDS_InitializeOutput(&resultsPage, SDDS_BINARY, 1, "Intra-beam scattering rates",
                             "Intra-beam scattering rates", outputfile))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  if (!SDDS_TransferParameterDefinition(&resultsPage, &twissPage, "I1", NULL) ||
      !SDDS_TransferParameterDefinition(&resultsPage, &twissPage, "I2", NULL) ||
      !SDDS_TransferParameterDefinition(&resultsPage, &twissPage, "I3", NULL) ||
      !SDDS_TransferParameterDefinition(&resultsPage, &twissPage, "I4", NULL) ||             
      !SDDS_TransferParameterDefinition(&resultsPage, &twissPage, "I5", NULL) ||             
      !SDDS_TransferParameterDefinition(&resultsPage, &twissPage, "pCentral", NULL) ||
      !SDDS_TransferParameterDefinition(&resultsPage, &twissPage, "taux", NULL) ||
      !SDDS_TransferParameterDefinition(&resultsPage, &twissPage, "taudelta", NULL) )
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);

  if (0>SDDS_DefineParameter(&resultsPage, "Superperiods", NULL, NULL, "Superperiods", NULL, SDDS_LONG, NULL) ||
      0>SDDS_DefineParameter(&resultsPage, "Energy", "E", "MeV", "Total Energy", NULL, SDDS_DOUBLE, NULL) ||
      0>SDDS_DefineParameter(&resultsPage, "Particles", NULL, NULL, "Particles", NULL, SDDS_DOUBLE, NULL) ||
      0>SDDS_DefineParameter(&resultsPage, "Charge", NULL, "nC", "Charge", NULL, SDDS_DOUBLE, NULL) ||
      0>SDDS_DefineParameter(&resultsPage, "PeakCurrent", "I$bp$n", "A", "Peak Current", NULL, SDDS_DOUBLE, NULL) ||
      0>SDDS_DefineParameter(&resultsPage, "xGrowthRate", "g$bIBS,x$n", "1/s", "IBS emittance growth rate in the horizontal plane", NULL, SDDS_DOUBLE, NULL) ||
      0>SDDS_DefineParameter(&resultsPage, "yGrowthRate", "g$bIBS,y$n", "1/s", "IBS emittance growth rate in the vertical plane", NULL, SDDS_DOUBLE, NULL) ||
      0>SDDS_DefineParameter(&resultsPage, "zGrowthRate", "g$bIBS,z$n", "1/s", "IBS emittance growth rate in the longitudinal plane", NULL, SDDS_DOUBLE, NULL))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);

  if (0>SDDS_DefineParameter(&resultsPage, "Convergence", NULL, NULL, "Convergence state of emittance calculations.", NULL, SDDS_STRING, NULL) ||
      0>SDDS_DefineParameter(&resultsPage, "emitx0", "$ge$r$bx,0$n", "m-rad", "Horizontal emittance with no coupling and no IBS.", NULL, SDDS_DOUBLE, NULL) ||
      0>SDDS_DefineParameter(&resultsPage, "emitxInput", "$ge$r$bx,Input$n", "m-rad", "Horizontal emittance with coupling and no IBS.", NULL, SDDS_DOUBLE, NULL) ||
      0>SDDS_DefineParameter(&resultsPage, "emityInput", "$ge$r$by,Input$n", "m-rad", "Vertical emittance with coupling and no IBS.", NULL, SDDS_DOUBLE, NULL) ||
      0>SDDS_DefineParameter(&resultsPage, "emitx", "$ge$r$bx$n", "m-rad", "Horizontal emittance with coupling and with IBS.", NULL, SDDS_DOUBLE, NULL) ||
      0>SDDS_DefineParameter(&resultsPage, "emity", "$ge$r$by$n", "m-rad", "Vertical emittance with coupling and with IBS.", NULL, SDDS_DOUBLE, NULL))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);

  if (0>SDDS_DefineParameter(&resultsPage, "sigmaDelta0", "$gs$bd,0$n", NULL, "Relative momentum spread without IBS.", NULL, SDDS_DOUBLE, NULL) ||
      0>SDDS_DefineParameter(&resultsPage, "sigmaDelta", "$gs$bd$n", NULL, "Relative momentum spread with IBS.", NULL, SDDS_DOUBLE, NULL) ||
      0>SDDS_DefineParameter(&resultsPage, "sigmaz0", "$gs$r$bz,0$n", "m", "Bunch length without IBS.", NULL, SDDS_DOUBLE, NULL) ||
      0>SDDS_DefineParameter(&resultsPage, "sigmaz", "$gs$r$bz$n", "m", "Bunch length with IBS.", NULL, SDDS_DOUBLE, NULL))
    SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  
  if (verbosity)
    fprintf( stderr, "Opening for reading \"%s\"\n", inputfile);
  if (!SDDS_InitializeInput(&twissPage, inputfile))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  if (!SDDS_WriteLayout(&resultsPage) )
    SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  
  while(SDDS_ReadPage(&twissPage)>0) {
    if (!SDDS_GetParameters(&twissPage,
                            "pCentral", &pCentral,
                            "I1", &I1,
                            "I2", &I2,
                            "I3", &I3,
                            "I4", &I4,
                            "I5", &I5,
                            "taux", &taux,
                            "taudelta", &taudelta,
                            NULL) )
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    EMeV = sqrt(sqr(pCentral) + 1) * me_mev;
    elements = SDDS_CountRowsOfInterest(&twissPage);
    s = SDDS_GetColumnInDoubles(&twissPage, "s");
    betax = SDDS_GetColumnInDoubles(&twissPage, "betax");
    betay = SDDS_GetColumnInDoubles(&twissPage, "betay");
    alphax = SDDS_GetColumnInDoubles(&twissPage, "alphax");
    alphay = SDDS_GetColumnInDoubles(&twissPage, "alphay");
    etax = SDDS_GetColumnInDoubles(&twissPage, "etax");
    etaxp = SDDS_GetColumnInDoubles(&twissPage, "etaxp");

    /* emitx0 and sigmaDelta0 are unperturbed quantities 
       (i.e. no coupling and no IBS) that
       zibs requires to internally calculate the quantum excitation.
       (zibs doesn't use the radiation intergrals but should!) 
       */
    emitx0 = 55.0/ (32.*sqrt(3.)) * hbar_mks * sqr(pCentral)/ (me_mks * c_mks)
      * I5 / (I2 - I4);
    sigmaDelta0 = sqrt(55.0/ (32.*sqrt(3.)) * hbar_mks * sqr(pCentral)/ (me_mks * c_mks)
      * I3 / (2 * I2 + I4));
    /* use unperturbed quantities in no input supplied. */
    if (!sigmaDeltaInput)
      sigmaDeltaInput = sigmaDelta0;
    if (!emitxInput)
      emitxInput = emitx0/ ( 1 + coupling);
    emityInput = emitxInput * coupling;
    sigmaz0 = length;
    sigmaz = sigmaz0;
    sigmaDelta = sigmaDeltaInput;
    emity = emityInput;
    emitx = emitxInput;
    
    if (!growthRatesOnly) {
      if (verbosity > 1) {
        fprintf (stderr, "Starting values:\nemitx: %10.5lg sigmaDelta %10.5lg.\n", emitx, sigmaDelta);
      }
      emitxOld = emitx;
      sigmaDeltaOld = sigmaDelta;
      xGuess = (double*) malloc(sizeof(double)*dimensions);
      dxGuess = (double*) malloc(sizeof(double)*dimensions);
      xLowerLimit = (double*) malloc(sizeof(double)*dimensions);
      xUpperLimit = (double*) malloc(sizeof(double)*dimensions);
      disable = (short*) malloc(sizeof(short)*dimensions);
      xGuess[0] = emitx;
      xGuess[1] = sigmaDelta;
      dxGuess[0] = emitx * 0.1;
      dxGuess[1] = sigmaDelta * 0.1;
      xLowerLimit[0] = emitx0/ (1 + coupling);
      xLowerLimit[1] = sigmaDelta0;
      xUpperLimit[0] = emitx0/ (1 + coupling) * 200;
      xUpperLimit[1] = MIN(sigmaDelta0 * 100, 1.0);
      /* assign other variables to array which are not supoosed
         to be varied by simplex minimization
         */
      xGuess[2] = pCentral;
      xGuess[3] = emity;
      xGuess[4] = sigmaz0;
      xGuess[5] = particles;
      xGuess[6] = emitx0;
      xGuess[7] = sigmaDelta0;
      xGuess[8] = taux;
      xGuess[9] = taudelta;
      xGuess[10] = coupling;
      xGuess[11] = elements;
      xGuess[12] = superperiods;
      xGuess[13] = verbosity;
      xLowerLimit[2] = pCentral;
      xLowerLimit[3] = emity;
      xLowerLimit[4] = sigmaz0;
      xLowerLimit[5] = particles;
      xLowerLimit[6] = emitx0;
      xLowerLimit[7] = sigmaDelta0;
      xLowerLimit[8] = taux;
      xLowerLimit[9] = taudelta;
      xLowerLimit[10] = coupling;
      xLowerLimit[11] = elements;
      xLowerLimit[12] = superperiods;
      xLowerLimit[13] = verbosity;
      xUpperLimit[2] = pCentral;
      xUpperLimit[3] = emity;
      xUpperLimit[4] = sigmaz0;
      xUpperLimit[5] = particles;
      xUpperLimit[6] = emitx0;
      xUpperLimit[7] = sigmaDelta0;
      xUpperLimit[8] = taux;
      xUpperLimit[9] = taudelta;
      xUpperLimit[10] = coupling;
      xUpperLimit[11] = elements;
      xUpperLimit[12] = superperiods;
      xUpperLimit[13] = verbosity;
      disable[0] = 0;
      disable[1] = 0;
      for (i=2 ; i<dimensions ; i++) {
        dxGuess[i] = 0.0;
        disable[i] = 1;
      }
      if (verbosity) {
        fprintf( stderr, "Doing simplex minimization...\n");
      }
      simplexMin( &yReturn, xGuess, dxGuess, xLowerLimit, xUpperLimit, disable, dimensions,
                 target, tolerance, IBSequations, verbosity?IBSsimplexReport:NULL, 
                 maxEvaluations, maxPasses);
      /* final answers */
      emitx = xGuess[0];
      sigmaDelta = xGuess[1];
      emity = emitx * coupling;
      sigmaz = sigmaz0 * (sigmaDelta/ sigmaDelta0);
    }
    IBSGrowthRates( pCentral, emitx, emity, sigmaDelta, sigmaz, particles,
                   emitx0, sigmaDelta0, 2./taux, 2./taudelta, coupling,
                   s, betax, alphax, betay, alphay, etax, etaxp, elements,
                   superperiods, 1,
                   &xGrowthRate, &yGrowthRate, &zGrowthRate);
    converged = 1;
    if (0>SDDS_StartPage(&resultsPage, 0) ||
        !SDDS_SetParameters(&resultsPage, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,
                            "Convergence", converged?"Emittance converged":"Emittance did not converge",
                            "pCentral", pCentral,
                            "I1", I1,
                            "I2", I2,
                            "I3", I3,
                            "I4", I4,
                            "I5", I5,
                            "taux", taux,
                            "taudelta", taudelta,
                            "Energy", EMeV,
                            "Particles", particles,
                            "Charge", (1e9 * charge),
                            "PeakCurrent", (charge*c_mks/(sqrt(2*PI)*sigmaz)),
                            "Superperiods", superperiods,
                            "emitx0", emitx0,
                            "emitxInput", emitxInput,
                            "emityInput", emityInput,
                            "xGrowthRate", xGrowthRate,
                            "yGrowthRate", yGrowthRate,
                            "zGrowthRate", zGrowthRate,
                            "emitx", emitx,
                            "emity", emity,
                            "sigmaDelta0", sigmaDelta0,
                            "sigmaDelta", sigmaDelta,
                            "sigmaz0", sigmaz0,
                            "sigmaz", sigmaz, NULL))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    if (!SDDS_WritePage(&resultsPage))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  
  if (!SDDS_Terminate(&twissPage) || !SDDS_Terminate(&resultsPage))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  
  exit(0);
}

double IBSequations(double *x, long *invalid) {
  double emitx, sigmaDelta;
  double pCentral, emity, sigmaz, sigmaz0, particles, emitx0, 
  sigmaDelta0, taux, taudelta, coupling;
  long elements, superperiods, verbosity;
  double xGrowthRate, yGrowthRate, zGrowthRate;
  double Gx, Gy, Gz;
  double a, b, c, d, e, f, func1, func2;
  long findEquilibrium;
  
  emitx = x[0];
  sigmaDelta = x[1];
  pCentral = x[2];
  sigmaz0 = x[4];
  particles = x[5];
  emitx0 = x[6];
  sigmaDelta0 = x[7];
  taux = x[8];
  taudelta = x[9];
  coupling = x[10];
  elements = x[11];
  superperiods = x[12];
  verbosity = x[13];

    /* zap code requires damping rate for horizontal and longitudinal emittances
       which is twice the damping rate for one coordinate.
       The quantities emitx0 and 2/taux are used to determined
       the quantum excitation inside the zap algortihm.
     */

    /* During iterations maintain the ratios
       sigmaz/sigmaDelta and emity/emitx constant.
       Their original values are
       sigmaz0/sigmaDelta0 and emityInput/emitxInput respectively.       
       */

    /* equations to solve (using ZAP manual notation) and throwing in coupling
       terms.
       damping term     quantum excitation      IBS terms
                       (a constant)                                              
        SR             SR   0    1              IBS       1         IBS       k  
       g   e        = g    e   -----       +   g    e   -----  +   g    e   -----
        x   x          x    x  1 + k            x    x  1 + k       y    y  1 + k
     
       The quantum excitation term is a constant and is reduced in the x plane 
       because of coupling. Also the IBS growth rate in x is reduced because of
       coupling. 
     
       In the y-plane:
       damping term     quantum excitation      IBS terms
                       (a constant)                                              
        SR             SR   0    k              IBS       1         IBS       k  
       g   e        = g    e   -----       +   g    e   -----  +   g    e   -----
        y   y          y    y  1 + k            y    y  1 + k       x    x  1 + k
     
       In the longitudinal plane,
       damping term     quantum excitation      IBS term
                       (a constant)                                              
        SR      2      SR        2              IBS      2
       g   delta    = g    delta0          +   g    delta
        z              z                        z    
     
       Assume that g^IBS will have the approximate dependence which will help finding
       solutions:
                    Input Input          2  
                   e    e     deltaInput          
        IBS         x    y                          1
       g    = G' ----------------------- = G ----------------
                               2                           2               
                   e  e   delta                e  e   delta 
                    x  y                        x  y        
       so that G = g^IBS from the very first calculation.
       Correspondence with our variables:
       g^SR_x    -> 2/taux
       g^IBS_x   -> xGrowthRate
       e^0_x     -> emitx0
       k         -> coupling

       One can ignore the y equation and simply put emity = coupling * emitx 
       during and after the calculation of emitx.
       */
    
  emity = emitx * coupling;
  sigmaz = sigmaz0 * (sigmaDelta/ sigmaDelta0);
  IBSGrowthRates( pCentral, emitx, emity, sigmaDelta, sigmaz, particles,
                 emitx0, sigmaDelta0, 2./taux, 2./taudelta, coupling,
                 s, betax, alphax, betay, alphay, etax, etaxp, elements,
                 superperiods, verbosity,
                 &xGrowthRate, &yGrowthRate, &zGrowthRate);
  Gx = xGrowthRate * emitx * emity * sqr(sigmaDelta);
  Gy = yGrowthRate * emitx * emity * sqr(sigmaDelta);
  Gz = zGrowthRate * emitx * emity * sqr(sigmaDelta);
  a = -2./taux;
  b = 2./taux * emitx0/(1+coupling);
  c = Gx/ (coupling * (1 + coupling));
  d = -2./taudelta;
  e = 2./taudelta * sqr(sigmaDelta0);
  f = Gz/ (coupling );
  func1 = a * emitx + b + c/ (emitx * sqr(sigmaDelta));
  func2 = d * sqr(sigmaDelta) + e + f/ sqr(emitx);
  *invalid = 0;
  return (sqr(func1/b) + sqr(func2/e));
}

void IBSsimplexReport(double ymin, double *xmin, long pass, long evals, long dims) {
  fprintf( stderr, "IBS Simplex Report:\nMinimum value obtained: %le.\n", ymin);
  fprintf( stderr, "    %ld passes, %ld evaluations.\n", pass, evals);
  fprintf( stderr, "    emitx = %le   sigmaDelta = %le.\n", xmin[0], xmin[1]);
  return;
}
