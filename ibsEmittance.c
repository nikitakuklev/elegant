/* 
 * $Log: not supported by cvs2svn $
 * Revision 1.9  2001/05/16 19:02:58  borland
 * Modified calls to simplexMin() to accomodate new argument.
 * Added simplified Rosenzweig/Serafini model for rf focusing in the
 * body of RFCA and RFCW elements.
 *
 * Revision 1.8  2000/10/23 18:57:08  borland
 * Implemented -energy option.
 *
 * Revision 1.7  1999/10/12 21:49:55  borland
 * All printouts now go to the stdout rather than stderr.  fflush statements,
 * some unnecessary, were added in a mostly automated fashion.
 *
 * Revision 1.6  1999/08/05 15:35:52  soliday
 * Added WIN32 and Linux support
 *
 * Revision 1.5  1999/03/18 21:02:44  borland
 * Implemented the -rf option.
 *
 * Revision 1.4  1999/02/15 14:41:26  borland
 * Added integration vs time.
 *
 * Revision 1.3  1999/02/11 20:46:35  borland
 * Now uses target parameter of simplexMin to accept any solution with
 * 1e-4 or better convergence.
 *
 * Revision 1.2  1999/01/27 17:55:15  borland
 * Removed prototypes from ibsEmittance.c and zibs.c and put them in common
 * file zibs.h
 *
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
 [-energy=<MeV>]\n\
 [-integrate=turns=<number>[,stepSize=<number>]]";

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
#define SET_TARGET 12
#define SET_INTEGRATE 13
#define N_OPTIONS 14
char *option[N_OPTIONS] = {
  "energy",
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
  "growthRatesOnly",
  "target",
  "integrate",
  };

#include "zibs.h"
double IBSequations(double *x, long *invalid);
void IBSsimplexReport(double ymin, double *xmin, long pass, long evals, long dims);

void IBSIntegrate(double *exInteg, double *eyInteg, double *elInteg, long *passInteg,
                  double *xRateInteg, double *yRateInteg, double *zRateInteg,
                  long integTurns, long integStepSize, 
                  double P, double emitx, double emity,
                  double sigmaDelta, double sigmaz,
                  double particles,
                  double emitx0, double sigmaDelta0, 
                  double transSRdampRate, double longSRdampRate,
                  double coupling,
                  double *s, double *betax, double *alphax, double *betay, 
                  double *alphay, double *etax, double *etaxp, long elements, 
                  long superperiods);

/* global variables */
double *s, *betax, *alphax, *betay, *alphay, *etax, *etaxp;

int main( int argc, char **argv)
{
  SCANNED_ARG *scanned;
  char *inputfile, *outputfile;
  SDDS_DATASET twissPage, resultsPage;
  double particles, charge, length;
  long verbosity, i, elements, superperiods, growthRatesOnly;
  double pCentral, I1, I2, I3, I4, I5, taux, taudelta;
  double EMeV;
  double emitx0, emitx, emitxInput, emityInput, emity, coupling, sigmaz0, sigmaz;
  double sigmaDelta0, sigmaDelta, sigmaDeltaInput, xGrowthRate, yGrowthRate, zGrowthRate;
  double emitxOld, sigmaDeltaOld;
  long method, converged;
/* used in simplex minimization */
  double yReturn, *xGuess, *dxGuess, *xLowerLimit, *xUpperLimit;
  short *disable;
  long dimensions = 14, maxEvaluations = 500, maxPasses = 2;
  double target = 1e-4, tolerance = 1e-6;
  long integrationTurns, integrationStepSize, integrationPoints = 0;
  double *exInteg, *eyInteg, *elInteg, *xRateInteg, *yRateInteg, *zRateInteg;
  long *passInteg;
  unsigned long dummyFlags;
  double rfVoltage, rfHarmonic;
  double alphac, U0, circumference, energy;
  SDDS_RegisterProgramName(argv[0]);
  argc  =  scanargs(&scanned, argc, argv);
  if (argc == 1)
    bomb(NULL, USAGE);

  inputfile  =  NULL;
  outputfile  =  NULL;
  energy = 0;
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
  integrationTurns = 0;
  rfVoltage = rfHarmonic = 0;
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
      case SET_TARGET:
        if (scanned[i].n_items!=2 ||
            !get_double(&target, scanned[i].list[1]) ||
            target<0)
          bomb("invalid -target syntax", NULL);
        break;
      case RF:
        if (scanned[i].n_items<2)
          bomb("invalid -rf syntax", NULL);
        scanned[i].n_items--;
        rfVoltage = rfHarmonic = 0;
        if (!scanItemList(&dummyFlags, scanned[i].list+1, &scanned[i].n_items, 0,
                          "voltage", SDDS_DOUBLE, &rfVoltage, 1, 0,
                          "harmonic", SDDS_DOUBLE, &rfHarmonic, 1, 0,
                          NULL) ||
            rfVoltage<=0 || rfHarmonic<=0)
          bomb("invalid -rf syntax/values", "-rf=voltage=MV,harmonic=<value>");
        break;
      case SET_ENERGY:
        if (scanned[i].n_items!=2)
          bomb("invalid -energy syntax", NULL);
        if (!sscanf(scanned[i].list[1], "%lf", &energy) || energy<=0)
          bomb("invalid -energy syntax/values", "-energy=<MeV>");
        break;
      case SET_INTEGRATE:
        if (scanned[i].n_items<2)
          bomb("invalid -integrate syntax", NULL);
        integrationTurns = 0;
        integrationStepSize = 1;
        scanned[i].n_items--;
        if (!scanItemList(&dummyFlags, scanned[i].list+1, &scanned[i].n_items, 0,
                          "turns", SDDS_LONG, &integrationTurns, 1, 0,
                          "stepsize", SDDS_LONG, &integrationStepSize, 1, 0,
                          NULL) ||
            integrationTurns<=0 || integrationStepSize<1) 
          bomb("invalid -integrate syntax", NULL);
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
  if (!length && !rfVoltage) 
    bomb("Specify either the bunch length or the rf voltage.", NULL);

  /***************************************************\
   * get parameter information from first input file  *
   \***************************************************/
  if (verbosity)
    fprintf( stdout, "Opening \"%s\" for checking presence of parameters.\n", inputfile);
  if (!SDDS_InitializeInput(&twissPage, inputfile))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  /* read first page of input file to get parameters 
     I1 I2 I3 I4 I5.
     Check presence of first radiation integral.
     */
  SDDS_ReadPage(&twissPage);
  /* parameter Type */
  switch(SDDS_CheckParameter(&twissPage, "I1", NULL, SDDS_DOUBLE, verbosity?stdout:NULL)) {
  case SDDS_CHECK_NONEXISTENT:
    if (verbosity)
      fprintf( stdout, "\tParameter I1 not found in input file.\n");
    break;
  case SDDS_CHECK_WRONGTYPE:
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    exit(1);
    break;
  case SDDS_CHECK_OKAY:
    break;
  default:
    fprintf( stdout, "Unexpected result from SDDS_CheckParameter routine while checking parameter Type.\n");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    exit(1);
    break;
  }
  if (verbosity)
    fprintf( stdout, "Opening \"%s\" for writing...\n", outputfile);
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
      0>SDDS_DefineParameter(&resultsPage, "emitx0", "$ge$r$bx,0$n", "$gp$rm", "Horizontal emittance with no coupling and no IBS.", NULL, SDDS_DOUBLE, NULL) ||
      0>SDDS_DefineParameter(&resultsPage, "emitxInput", "$ge$r$bx,Input$n", "$gp$rm", "Horizontal emittance with coupling and no IBS.", NULL, SDDS_DOUBLE, NULL) ||
      0>SDDS_DefineParameter(&resultsPage, "emityInput", "$ge$r$by,Input$n", "$gp$rm", "Vertical emittance with coupling and no IBS.", NULL, SDDS_DOUBLE, NULL) ||
      0>SDDS_DefineParameter(&resultsPage, "emitx", "$ge$r$bx$n", "$gp$rm", "Horizontal emittance with coupling and with IBS.", NULL, SDDS_DOUBLE, NULL) ||
      0>SDDS_DefineParameter(&resultsPage, "emity", "$ge$r$by$n", "$gp$rm", "Vertical emittance with coupling and with IBS.", NULL, SDDS_DOUBLE, NULL))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);

  if (0>SDDS_DefineParameter(&resultsPage, "sigmaDelta0", "$gs$r$bd,0$n", NULL, "Relative momentum spread without IBS.", NULL, SDDS_DOUBLE, NULL) ||
      0>SDDS_DefineParameter(&resultsPage, "sigmaDelta", "$gs$r$bd$n", NULL, "Relative momentum spread with IBS.", NULL, SDDS_DOUBLE, NULL) ||
      0>SDDS_DefineParameter(&resultsPage, "sigmaz0", "$gs$r$bz,0$n", "m", "Bunch length without IBS.", NULL, SDDS_DOUBLE, NULL) ||
      0>SDDS_DefineParameter(&resultsPage, "sigmaz", "$gs$r$bz$n", "m", "Bunch length with IBS.", NULL, SDDS_DOUBLE, NULL))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);

  if (integrationTurns) {
    if (SDDS_DefineColumn(&resultsPage, "ex", "$ge$r$bx$n", "$gp$rm", "Horizontal Emittance", NULL, SDDS_DOUBLE, 0)<0 ||
        SDDS_DefineColumn(&resultsPage, "ey", "$ge$r$by$n", "$gp$rm", "Vertical Emittance", NULL, SDDS_DOUBLE, 0)<0 ||
        SDDS_DefineColumn(&resultsPage, "el", "$ge$r$bl$n", "s", "Longitudinal Emittance", NULL, SDDS_DOUBLE, 0)<0 ||
        SDDS_DefineColumn(&resultsPage, "IBSRatex", NULL, "s", "Horizontal IBS Emittance Growth Rate", NULL, SDDS_DOUBLE, 0)<0 ||
        SDDS_DefineColumn(&resultsPage, "IBSRatey", NULL, "s", "Vertical IBS Emittance Growth Rate", NULL, SDDS_DOUBLE, 0)<0 ||
        SDDS_DefineColumn(&resultsPage, "IBSRatel", NULL, "s", "Longitudinal IBS Emittance Growth Rate", NULL, SDDS_DOUBLE, 0)<0 ||
        SDDS_DefineColumn(&resultsPage, "Pass", NULL, NULL, NULL, NULL, SDDS_LONG, 0)<0)
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    integrationPoints = integrationTurns/integrationStepSize+1;
    if (!(exInteg = SDDS_Malloc(sizeof(*exInteg)*integrationPoints)) ||
        !(eyInteg = SDDS_Malloc(sizeof(*eyInteg)*integrationPoints)) ||
        !(elInteg = SDDS_Malloc(sizeof(*elInteg)*integrationPoints)) ||
        !(xRateInteg = SDDS_Malloc(sizeof(*xRateInteg)*integrationPoints)) ||
        !(yRateInteg = SDDS_Malloc(sizeof(*yRateInteg)*integrationPoints)) ||
        !(zRateInteg = SDDS_Malloc(sizeof(*zRateInteg)*integrationPoints)) ||
        !(passInteg = SDDS_Malloc(sizeof(*passInteg)*integrationPoints)))
      bomb("memory allocation failure (integration arrays)", NULL);
  }
  
  if (verbosity)
    fprintf( stdout, "Opening for reading \"%s\"\n", inputfile);
  if (!SDDS_InitializeInput(&twissPage, inputfile))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  if (!SDDS_WriteLayout(&resultsPage) )
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  
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
                            "alphac", &alphac,
                            "U0", &U0,
                            NULL) )
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    EMeV = sqrt(sqr(pCentral) + 1) * me_mev;
    elements = SDDS_CountRowsOfInterest(&twissPage);
    s = SDDS_GetColumnInDoubles(&twissPage, "s");
    circumference = s[elements-1]*superperiods;
    U0 *= superperiods;
    if (energy!=0) {
      /* scale to new energy */
      pCentral = sqrt(sqr(energy/me_mev)-1);
      taux /= ipow(energy/EMeV, 3);
      taudelta /= ipow(energy/EMeV, 3);
      U0 *= ipow(energy/EMeV, 4);
      EMeV = energy;
    }
    
    if (!length && U0>rfVoltage)
      bomb("energy loss per turn is greater than rf voltage", NULL);
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
    sigmaDelta = sigmaDeltaInput;
    if (length)
      sigmaz0 = length;
    else {
      /* compute length in m from rf voltage, energy spread, etc */
      sigmaz0 = 
        circumference*sigmaDelta*
          sqrt(alphac*EMeV/(PIx2*rfHarmonic*sqrt(sqr(rfVoltage)-sqr(U0))));
    }
    sigmaz = sigmaz0;
    emity = emityInput;
    emitx = emitxInput;

    if (integrationPoints) {
      IBSIntegrate(exInteg, eyInteg, elInteg, passInteg,
                   xRateInteg, yRateInteg, zRateInteg,
                   integrationTurns, integrationStepSize, 
                   pCentral, emitx, emity, sigmaDelta, sigmaz, particles,
                   emitx0, sigmaDelta0, 2./taux, 2./taudelta, coupling,
                   s, betax, alphax, betay, alphay, etax, etaxp, elements,
                   superperiods);
    }
    
    if (!growthRatesOnly) {
      if (verbosity > 1) {
        fprintf (stdout, "Starting values:\nemitx: %10.5g sigmaDelta %10.5g.\n", emitx, sigmaDelta);
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
        fprintf( stdout, "Doing simplex minimization...\n");
      }
      simplexMin( &yReturn, xGuess, dxGuess, xLowerLimit, xUpperLimit, disable, dimensions,
                 target, tolerance, IBSequations, verbosity?IBSsimplexReport:NULL, 
                 maxEvaluations, maxPasses, 12, 0);
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
    if (0>SDDS_StartPage(&resultsPage, integrationPoints) ||
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
                            "sigmaz", sigmaz, NULL) ||
        (integrationPoints && 
         (!SDDS_SetColumn(&resultsPage, SDDS_SET_BY_NAME, exInteg, integrationPoints, "ex") ||
          !SDDS_SetColumn(&resultsPage, SDDS_SET_BY_NAME, eyInteg, integrationPoints, "ey") ||
          !SDDS_SetColumn(&resultsPage, SDDS_SET_BY_NAME, elInteg, integrationPoints, "el") ||
          !SDDS_SetColumn(&resultsPage, SDDS_SET_BY_NAME, xRateInteg, integrationPoints, "IBSRatex") ||
          !SDDS_SetColumn(&resultsPage, SDDS_SET_BY_NAME, yRateInteg, integrationPoints, "IBSRatey") ||
          !SDDS_SetColumn(&resultsPage, SDDS_SET_BY_NAME, zRateInteg, integrationPoints, "IBSRatel") ||
          !SDDS_SetColumn(&resultsPage, SDDS_SET_BY_NAME, passInteg, integrationPoints, "Pass"))) ||
        !SDDS_WritePage(&resultsPage))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  
  if (!SDDS_Terminate(&twissPage) || !SDDS_Terminate(&resultsPage))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  
  return(0);
}

double IBSequations(double *x, long *invalid) {
  double emitx, sigmaDelta;
  double pCentral, emity, sigmaz, sigmaz0, particles, emitx0, 
  sigmaDelta0, taux, taudelta, coupling;
  long elements, superperiods, verbosity;
  double xGrowthRate, yGrowthRate, zGrowthRate;
  double Gx, Gy, Gz;
  double a, b, c, d, e, f, func1, func2;
  
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
  fprintf( stdout, "IBS Simplex Report:\nMinimum value obtained: %e.\n", ymin);
  fprintf( stdout, "    %ld passes, %ld evaluations.\n", pass, evals);
  fprintf( stdout, "    emitx = %e   sigmaDelta = %e.\n", xmin[0], xmin[1]);
  return;
}

void IBSIntegrate(double *exInteg, double *eyInteg, double *elInteg, long *passInteg,
                  double *xRateInteg, double *yRateInteg, double *zRateInteg,
                  long integTurns, long integStepSize, 
                  double P, double emitx, double emity,
                  double sigmaDelta, double sigmaz,
                  double particles,
                  double emitx0, double sigmaDelta0, 
                  double transSRdampRate, double longitSRdampRate,
                  double coupling,
                  double *s, double *betax, double *alphax, double *betay, 
                  double *alphay, double *etax, double *etaxp, long elements, 
                  long superperiods)
{
  long turn, slot;
  double dT, gamma, vz, emitz, zRatio;
  double xGrowthRate, yGrowthRate, zGrowthRate, emitz0;
  
  gamma = sqrt(sqr(P)+1);
  vz = c_mks*P/gamma;
  dT = s[elements-1]*superperiods*integStepSize/vz;

  emitz = sigmaDelta*sigmaz;
  emitz0 = sigmaDelta0*sigmaz;
  zRatio = sigmaDelta/sigmaz;
  for (turn=slot=0; turn<integTurns; turn+=integStepSize, slot++) {
    exInteg[slot] = emitx;
    eyInteg[slot] = emity;
    elInteg[slot] = emitz/vz;
    passInteg[slot] = turn;
    IBSGrowthRates(P, emitx, emity, sigmaDelta, sigmaz, particles,
                   emitx0, sigmaDelta0, transSRdampRate, longitSRdampRate, coupling,
                   s, betax, alphax, betay, alphay, etax, etaxp, elements,
                   superperiods, 0,
                   &xGrowthRate, &yGrowthRate, &zGrowthRate);
    xRateInteg[slot] = xGrowthRate;
    yRateInteg[slot] = yGrowthRate;
    zRateInteg[slot] = zGrowthRate;
    emitx += (xGrowthRate-transSRdampRate)*emitx*dT+transSRdampRate*emitx0*dT/(1+coupling);
    emity += (yGrowthRate-transSRdampRate)*emity*dT+transSRdampRate*emitx0*coupling*dT/(1+coupling);
    emitz += (zGrowthRate-longitSRdampRate)*emitz*dT+longitSRdampRate*emitz0*dT;
    sigmaDelta = sqrt(emitz*zRatio);
    sigmaz = emitz/sigmaDelta;
  }
  exInteg[slot] = emitx;
  eyInteg[slot] = emity;
  elInteg[slot] = emitz/vz;
  passInteg[slot] = turn;
  xRateInteg[slot] = xGrowthRate;
  yRateInteg[slot] = yGrowthRate;
  zRateInteg[slot] = zGrowthRate;
}


