/************************************************************************* \
* Copyright (c) 2017 The University of Chicago, as Operator of Argonne
* National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: ionEffects.c
 * purpose: simulation of ion interaction with the beam
 *
 * Joe Calvey, Michael Borland 2017
 */
#if defined(SOLARIS) && !defined(__GNUC__)
#include <sunmath.h>
#endif

#include <complex>
#include "mdb.h"
#include "track.h"
#include "ionEffects.h"
#include "constants.h"
#include "pressureData.h"
//#include <algorithm>

#define ION_FIELD_GAUSSIAN 0
#define ION_FIELD_BIGAUSSIAN 1
#define ION_FIELD_BILORENTZIAN 2
#define ION_FIELD_TRIGAUSSIAN 3
#define ION_FIELD_TRILORENTZIAN 4
#define ION_FIELD_AUTO 5
#define N_ION_FIELD_METHODS 6
static char *ionFieldMethodOption[N_ION_FIELD_METHODS] = {
  (char*)"gaussian",
  (char*)"bigaussian",
  (char*)"bilorentzian",
  (char*)"trigaussian",
  (char*)"trilorentzian",
  (char*)"auto"
};
static long ionFieldMethod = -1;

#define ION_FIT_RESIDUAL_SUM_ABS_DEV 0
#define ION_FIT_RESIDUAL_RMS_DEV 1
#define ION_FIT_RESIDUAL_MAX_ABS_DEV 2
#define ION_FIT_RESIDUAL_MAX_PLUS_RMS_DEV 3
#define ION_FIT_RESIDUAL_SUM_ABS_PLUS_RMS_DEV 4
#define ION_FIT_RESIDUAL_RMS_DEV_PLUS_ABS_DEV_SUM 5
#define ION_FIT_RESIDUAL_SUM_ABS_PLUS_ABS_DEV_SUM 6
#define ION_FIT_RESIDUAL_RMS_DEV_PLUS_CENTROID 7
#define N_ION_FIT_RESIDUAL_OPTIONS 8
static char *ionFitResidualOption[N_ION_FIT_RESIDUAL_OPTIONS] = {
  (char*)"sum-ad",
  (char*)"rms-dev",
  (char*)"max-ad",
  (char*)"max-ad-plus-rms-dev",
  (char*)"sum-ad-plus-rms-dev",
  (char*)"rms-dev-plus-ad-sum",
  (char*)"sum-ad-plus-ad-sum",
  (char*)"rms-dev-plus-centroid",
};

static long residualType = -1;

static long ionsInitialized = 0;

static PRESSURE_DATA pressureData;

typedef struct {
  long nSpecies;
  char **ionName;
  double *mass, *chargeState;     /* Atomic mass, charge state */
  long *sourceGasIndex;           /* If charge state=1, index in pressureData of the source gas */
  long *sourceIonIndex;           /* If charge state>1, index in this list of the ion that sources this ion */
  double *crossSection;
} ION_PROPERTIES;

static ION_PROPERTIES ionProperties;
void readIonProperties(char *filename);

void addIons(IONEFFECTS *ionEffects, long iSpecies, long nToAdd, double qToAdd, double centroid[2], double sigma[2], long symmetrize);

void addIon_point(IONEFFECTS *ionEffects, long iSpecies, double qToAdd,  double x, double y);

void gaussianBeamKick(double *coord, double center[2], double sigma[2], double kick[2], double charge, 
		      double ionMass, double ionCharge);

void roundGaussianBeamKick(double *coord, double center[2], double sigma[2], double kick[2], double charge, 
		      double ionMass, double ionCharge);

void makeIonHistograms(IONEFFECTS *ionEffects, long nSpecies, double *beamCentroid, double *beamSigma, 
		       double *ionCentroid, double *ionSigma);
double findIonBinningRange(IONEFFECTS *ionEffects, long iPlane, long nSpecies);

#if USE_MPI
void shareIonHistograms(IONEFFECTS *ionEffects);
#endif
void determineOffsetAndActiveBins(double *histogram, long nBins, long *binOffset, long *activeBins);

static SDDS_DATASET *SDDS_beamOutput = NULL;
static SDDS_DATASET *SDDS_ionDensityOutput = NULL;
static SDDS_DATASET *SDDS_ionHistogramOutput = NULL;
static long ionHistogramOutputInterval, ionHistogramMinOutputBins, ionHistogramMaxBins;
static double ionHistogramOutput_sStart, ionHistogramOutput_sEnd;

static double sStartFirst = -1;
static long iIonEffectsElement = -1, nIonEffectsElements = 0, iBeamOutput, iIonDensityOutput;
static IONEFFECTS *firstIonEffects = NULL; /* first in the lattice */
#if USE_MPI
static long leftIonCounter = 0;
extern void find_global_min_index (double *min, int *processor_ID, MPI_Comm comm);
#endif

//for fit (e.g., bi-gaussian)
static double *xData=NULL, *yData=NULL, *yFit=NULL, yDataSum;
static long nData = 0;
static long nFunctions = 2; /* should be 2 or 3 */
static long mFunctions;

short multipleWhateverFit(double beamSigma[2], double beamCentroid[2], double paramValueX[9], 
                    double paramValueY[9], IONEFFECTS *ionEffects, double ionSigma[2], double ionCentroid[2]);
double multiGaussianFunction(double *param, long *invalid);
double multiLorentzianFunction(double *param, long *invalid);

//void report();
void report(double res, double *a, long pass, long n_eval, long n_dimen);

#if USE_MPI
void findGlobalMinIndex (double *min, int *processor_ID, MPI_Comm comm) {
    struct {
      double val;
      int rank;
    } in, out;
    in.val = *min;
    MPI_Comm_rank(comm, &(in.rank));
    MPI_Allreduce(&in, &out, 1, MPI_DOUBLE_INT, MPI_MINLOC, comm);
    *min = out.val;
    *processor_ID = out.rank;
}
void findGlobalMaxIndex (double *max, int *processor_ID, MPI_Comm comm) {
    struct {
      double val;
      int rank;
    } in, out;
    in.val = *max;
    MPI_Comm_rank(comm, &(in.rank));
    MPI_Allreduce(&in, &out, 1, MPI_DOUBLE_INT, MPI_MAXLOC, comm);
    *max = out.val;
    *processor_ID = out.rank;
}
#endif

char speciesNameBuffer[100];
char *makeSpeciesName(const char *prefix, char *suffix)
{
  strcpy(speciesNameBuffer, prefix);
  strcat(speciesNameBuffer, suffix);
  return &speciesNameBuffer[0];
}

  
void closeIonEffectsOutputFiles() {
  if (SDDS_beamOutput) {
    SDDS_Terminate(SDDS_beamOutput);
    SDDS_beamOutput = NULL;
  }
  if (SDDS_ionDensityOutput) {
    SDDS_Terminate(SDDS_ionDensityOutput);
    SDDS_ionDensityOutput = NULL;
  }
}

void setUpIonEffectsOutputFiles(long nPasses) 
{
  long iSpecies;
  closeIonEffectsOutputFiles(); /* Shouldn't be needed, but won't hurt */

  if (beam_output) {
    /* Setup the beam output file */
#if USE_MPI 
    if (myid==0) {
#endif
      if (!SDDS_beamOutput) {
        SDDS_beamOutput = (SDDS_DATASET*)tmalloc(sizeof(*SDDS_beamOutput));
        if (!SDDS_InitializeOutput(SDDS_beamOutput, SDDS_BINARY, 1, "electron beam output", NULL, beam_output)) {
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
          exitElegant(1);
        }
        if (!SDDS_DefineSimpleColumn(SDDS_beamOutput, "t", "s", SDDS_DOUBLE) ||
            !SDDS_DefineSimpleColumn(SDDS_beamOutput, "Pass", NULL, SDDS_LONG) ||
            !SDDS_DefineSimpleColumn(SDDS_beamOutput, "Bunch", NULL, SDDS_LONG) ||
            !SDDS_DefineSimpleColumn(SDDS_beamOutput, "qBunch", "C", SDDS_DOUBLE) ||
            !SDDS_DefineSimpleColumn(SDDS_beamOutput, "npBunch", NULL, SDDS_LONG) ||
            !SDDS_DefineSimpleColumn(SDDS_beamOutput, "s", "m", SDDS_DOUBLE) ||
            !SDDS_DefineSimpleColumn(SDDS_beamOutput, "Sx", "m", SDDS_DOUBLE) ||
            !SDDS_DefineSimpleColumn(SDDS_beamOutput, "Sy", "m", SDDS_DOUBLE) ||
            !SDDS_DefineSimpleColumn(SDDS_beamOutput, "Cx", "m", SDDS_DOUBLE) ||
            !SDDS_DefineSimpleColumn(SDDS_beamOutput, "Cy", "m", SDDS_DOUBLE) ||
            !SDDS_SaveLayout(SDDS_beamOutput) || !SDDS_WriteLayout(SDDS_beamOutput)) {
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
          exitElegant(1);
        }
      }
#if USE_MPI
    }
#endif
  }
  
  if (ion_density_output) {
    /* Setup the ion density file */
#if USE_MPI 
    if (myid==0) {
#endif
      if (!SDDS_ionDensityOutput) {
        SDDS_ionDensityOutput = (SDDS_DATASET*)tmalloc(sizeof(*SDDS_ionDensityOutput));
        if (!SDDS_InitializeOutput(SDDS_ionDensityOutput, SDDS_BINARY, 1, "ion density output", NULL, ion_density_output)) {
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
          exitElegant(1);
        }
        if (!SDDS_DefineSimpleColumn(SDDS_ionDensityOutput, "Pass", NULL, SDDS_LONG) ||
            !SDDS_DefineSimpleColumn(SDDS_ionDensityOutput, "Bunch", NULL, SDDS_LONG) ||
            !SDDS_DefineSimpleColumn(SDDS_ionDensityOutput, "t", "s", SDDS_DOUBLE) ||
            !SDDS_DefineSimpleColumn(SDDS_ionDensityOutput, "s", "m", SDDS_DOUBLE) ||
            !SDDS_DefineSimpleColumn(SDDS_ionDensityOutput, "qIons", "C", SDDS_DOUBLE) ||
            !SDDS_DefineSimpleColumn(SDDS_ionDensityOutput, "nMacroIons", NULL, SDDS_LONG) ||
            !SDDS_DefineSimpleColumn(SDDS_ionDensityOutput, "Sx", "m", SDDS_DOUBLE) ||
            !SDDS_DefineSimpleColumn(SDDS_ionDensityOutput, "Sy", "m", SDDS_DOUBLE) ||
            !SDDS_DefineSimpleColumn(SDDS_ionDensityOutput, "Cx", "m", SDDS_DOUBLE) ||
            !SDDS_DefineSimpleColumn(SDDS_ionDensityOutput, "Cy", "m", SDDS_DOUBLE)) {
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
          exitElegant(1);
        }
        if (ion_species_output) {
          for (iSpecies=0; iSpecies<ionProperties.nSpecies; iSpecies++) {
            if (!SDDS_DefineSimpleColumn(SDDS_ionDensityOutput,
                                         makeSpeciesName("qIons_", ionProperties.ionName[iSpecies]),
                                         "C", SDDS_DOUBLE) ||
                !SDDS_DefineSimpleColumn(SDDS_ionDensityOutput,
                                         makeSpeciesName("nMacroIons_", ionProperties.ionName[iSpecies]),
                                         NULL, SDDS_LONG) ||
                !SDDS_DefineSimpleColumn(SDDS_ionDensityOutput,
                                         makeSpeciesName("Cx_", ionProperties.ionName[iSpecies]),
                                         "m", SDDS_DOUBLE) ||
                !SDDS_DefineSimpleColumn(SDDS_ionDensityOutput,
                                         makeSpeciesName("Cy_", ionProperties.ionName[iSpecies]),
                                         "m", SDDS_DOUBLE) ||
                !SDDS_DefineSimpleColumn(SDDS_ionDensityOutput,
                                         makeSpeciesName("Sx_", ionProperties.ionName[iSpecies]),
                                         "m", SDDS_DOUBLE) ||
                !SDDS_DefineSimpleColumn(SDDS_ionDensityOutput,
                                         makeSpeciesName("Sy_", ionProperties.ionName[iSpecies]),
                                         "m", SDDS_DOUBLE) 
                ) {
              SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
              exitElegant(1);
            }
          }
        }
#if USE_MPI
        if (!SDDS_DefineSimpleColumn(SDDS_ionDensityOutput, "nMacroIonsMin", NULL, SDDS_LONG) ||
          !SDDS_DefineSimpleColumn(SDDS_ionDensityOutput, "nMacroIonsMax", NULL, SDDS_LONG)) {
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
          exitElegant(1);
        }
#endif
        if (!SDDS_SaveLayout(SDDS_ionDensityOutput) || !SDDS_WriteLayout(SDDS_ionDensityOutput)) {
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
          exitElegant(1);
        }
      }
#if USE_MPI
    }
#endif
  }

  if (ion_histogram_output) {
    /* Setup the ion histogram output file */
#if USE_MPI 
    if (myid==0) {
#endif
      if (!SDDS_ionHistogramOutput) {
        SDDS_ionHistogramOutput = (SDDS_DATASET*)tmalloc(sizeof(*SDDS_ionHistogramOutput));
        if (!SDDS_InitializeOutput(SDDS_ionHistogramOutput, SDDS_BINARY, 1,
				   "ion histogram output", NULL, ion_histogram_output)) {
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
          exitElegant(1);
        }
        if (!SDDS_DefineSimpleParameter(SDDS_ionHistogramOutput, "Pass", NULL, SDDS_LONG) ||
            !SDDS_DefineSimpleParameter(SDDS_ionHistogramOutput, "Bunch", NULL, SDDS_LONG) ||
            !SDDS_DefineSimpleParameter(SDDS_ionHistogramOutput, "t", "s", SDDS_DOUBLE) ||
            !SDDS_DefineSimpleParameter(SDDS_ionHistogramOutput, "s", "m", SDDS_DOUBLE) ||
	    !SDDS_DefineSimpleParameter(SDDS_ionHistogramOutput, "Plane", NULL, SDDS_STRING) ||
            !SDDS_DefineSimpleParameter(SDDS_ionHistogramOutput, "qIonsOutside", "C", SDDS_DOUBLE) ||
            !SDDS_DefineSimpleParameter(SDDS_ionHistogramOutput, "fractionIonChargeOutside", NULL, SDDS_DOUBLE) ||
            !SDDS_DefineSimpleParameter(SDDS_ionHistogramOutput, "binSize", "m", SDDS_DOUBLE) ||
            !SDDS_DefineSimpleParameter(SDDS_ionHistogramOutput, "binRange", "m", SDDS_DOUBLE) ||
            !SDDS_DefineSimpleParameter(SDDS_ionHistogramOutput, "nBins", NULL, SDDS_LONG) ||
	    !SDDS_DefineSimpleColumn(SDDS_ionHistogramOutput, "Position", "m", SDDS_DOUBLE) ||
	    !SDDS_DefineSimpleColumn(SDDS_ionHistogramOutput, "Charge", "C", SDDS_DOUBLE)) {
	  SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
	  exitElegant(1);
        }
	if (ionFieldMethod!=ION_FIELD_GAUSSIAN &&
	    (!SDDS_DefineSimpleColumn(SDDS_ionHistogramOutput, "ChargeFit", "C", SDDS_DOUBLE) ||
             !SDDS_DefineSimpleParameter(SDDS_ionHistogramOutput, "fitType", NULL, SDDS_STRING) ||
	     !SDDS_DefineSimpleParameter(SDDS_ionHistogramOutput, "fitResidual", NULL, SDDS_DOUBLE) ||
#if USE_MPI
	     !SDDS_DefineSimpleParameter(SDDS_ionHistogramOutput, "nEvaluationsBest", NULL, SDDS_LONG) ||
	     !SDDS_DefineSimpleParameter(SDDS_ionHistogramOutput, "nEvaluationsMin", NULL, SDDS_LONG) ||
	     !SDDS_DefineSimpleParameter(SDDS_ionHistogramOutput, "nEvaluationsMax", NULL, SDDS_LONG) ||
#else
	     !SDDS_DefineSimpleParameter(SDDS_ionHistogramOutput, "nEvaluations", NULL, SDDS_LONG) ||
#endif
	     !SDDS_DefineSimpleParameter(SDDS_ionHistogramOutput, "sigma1", "m", SDDS_DOUBLE) ||
	     !SDDS_DefineSimpleParameter(SDDS_ionHistogramOutput, "centroid1", "m", SDDS_DOUBLE) ||
	     !SDDS_DefineSimpleParameter(SDDS_ionHistogramOutput, "q1", "C", SDDS_DOUBLE) ||
	     !SDDS_DefineSimpleParameter(SDDS_ionHistogramOutput, "sigma2", "m", SDDS_DOUBLE) ||
	     !SDDS_DefineSimpleParameter(SDDS_ionHistogramOutput, "centroid2", "m", SDDS_DOUBLE) ||
	     !SDDS_DefineSimpleParameter(SDDS_ionHistogramOutput, "q2", "C", SDDS_DOUBLE))) {
	  SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
	  exitElegant(1);
	}
        if ((ionFieldMethod==ION_FIELD_TRIGAUSSIAN || ionFieldMethod==ION_FIELD_TRILORENTZIAN) &&
            (!SDDS_DefineSimpleParameter(SDDS_ionHistogramOutput, "sigma3", "m", SDDS_DOUBLE) ||
	     !SDDS_DefineSimpleParameter(SDDS_ionHistogramOutput, "centroid3", "m", SDDS_DOUBLE) ||
	     !SDDS_DefineSimpleParameter(SDDS_ionHistogramOutput, "q3", "C", SDDS_DOUBLE))) {
	  SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
	  exitElegant(1);
        }
        if (!SDDS_SaveLayout(SDDS_ionHistogramOutput) || !SDDS_WriteLayout(SDDS_ionHistogramOutput)) {
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
          exitElegant(1);
        }
      }
#if USE_MPI
    }
#endif
  }  
}

void setupIonEffects(NAMELIST_TEXT *nltext, VARY *control, RUN *run)
{
  /* process namelist input */
  set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
  set_print_namelist_flags(0);
  if (processNamelist(&ion_effects, nltext)==NAMELIST_ERROR)
    bombElegant(NULL, NULL);
  if (echoNamelists) print_namelist(stdout, &ion_effects);

  /* Basic check of input values */
  if (macro_ions<=0)
    bombElegant("macro_ions must be positive", NULL);
  if (generation_interval<=0)
    bombElegant("generation_interval must be positive", NULL);
  if (!pressure_profile || !strlen(pressure_profile))
    bombElegant("pressure_profile undefined", NULL);
  if (!ion_properties || !strlen(ion_properties))
    bombElegant("ion_properties undefined", NULL);
  if (beam_output)
    beam_output = compose_filename(beam_output, run->rootname);
  if (ion_density_output)
    ion_density_output = compose_filename(ion_density_output, run->rootname);
  if (ion_histogram_output) {
    ion_histogram_output = compose_filename(ion_histogram_output, run->rootname);
    if (ion_histogram_output_s_start>ion_histogram_output_s_end)
      bombElegantVA((char*)"ion_histogram_s_start (%le) is > ion_histogram_s_end (%le)", ion_histogram_output_s_start,
		    ion_histogram_output_s_end);
    if (ion_histogram_output_interval<=0)
      bombElegantVA((char*)"ion_histogram_s_interval (%ld) is <=0", ion_histogram_output_interval);
    ionHistogramOutputInterval = ion_histogram_output_interval;
    ionHistogramOutput_sStart = ion_histogram_output_s_start;
    ionHistogramOutput_sEnd= ion_histogram_output_s_end;
    ionHistogramMinOutputBins = ion_histogram_min_output_bins;
    ionHistogramMaxBins = ion_histogram_max_bins;
  }
  if (!field_calculation_method || !strlen(field_calculation_method))
    bombElegant("field_calculation_method undefined", NULL);
  if ((ionFieldMethod = match_string(field_calculation_method, ionFieldMethodOption, N_ION_FIELD_METHODS, EXACT_MATCH))<0)
    bombElegantVA((char*)"field_calculation_method=\"%s\" not recognized", field_calculation_method);
  if (ionFieldMethod==ION_FIELD_AUTO)
    bombElegantVA((char*)"field_calculation_method=\"%s\" is not yet implemented", field_calculation_method);
  if (ionFieldMethod==ION_FIELD_BIGAUSSIAN || ionFieldMethod==ION_FIELD_BILORENTZIAN)
    nFunctions = 2;
  else if (ionFieldMethod==ION_FIELD_TRIGAUSSIAN || ionFieldMethod==ION_FIELD_TRILORENTZIAN)
    nFunctions = 3;

  if (!fit_residual_type || !strlen(fit_residual_type))
    residualType = ION_FIT_RESIDUAL_RMS_DEV_PLUS_ABS_DEV_SUM;
  else if ((residualType = match_string(fit_residual_type, ionFitResidualOption, N_ION_FIT_RESIDUAL_OPTIONS, EXACT_MATCH))<0)
    bombElegantVA((char*)"fit_residual_type=\"%s\" not recognized", fit_residual_type);

  for (int iPlane=0; iPlane<2; iPlane++) {
    if (ion_span[iPlane]<=0)
      bombElegantVA("ion_span must be positive for both planes---%le given for %s plane\n", 
                    ion_span[iPlane], iPlane?"y":"x");

    if (ion_bin_divisor[iPlane]<=0)
      bombElegantVA("ion_bin_divisor must be positive for both planes---%le given for %s plane\n", 
                    ion_bin_divisor[iPlane], iPlane?"y":"x");

  }

  readGasPressureData(pressure_profile, &pressureData);

  readIonProperties(ion_properties);

  setUpIonEffectsOutputFiles(control->n_passes);
}

void readIonProperties(char *filename)
{
  /* Assumed file structure:
   * Columns:
   * IonName       --- SDDS_STRING, Name of the ion, e.g., "H2O+", "CO++"
   * Mass          --- SDDS_FLOAT or SDDS_DOUBLE, in AMU
   * ChargeState   --- SDDS_LONG or SDDS_SHORT, ion charge state (positive integer)
   * SourceName    --- SDDS_STRING, Name of the source molecule for this ion, e.g., "H2O", "CO+"
   * CrossSection  --- SDDS_FLOAT or SDDS_DOUBLE, Cross section for producing ion from source, in Mb (megabarns)
   */

  SDDS_DATASET SDDSin;
  long i;
  char **sourceName;

  if (!SDDS_InitializeInputFromSearchPath(&SDDSin, filename)) {
    printf("Problem opening ion properties data file %s\n", filename);
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
  }
  if (!check_sdds_column(&SDDSin, (char*)"Mass", (char*)"AMU"))
    bombElegantVA((char*)"Column \"Mass\" is missing, not floating-point type, or does not have units of \"AMU\" in %s\n",
                  filename);
  if (!check_sdds_column(&SDDSin, (char*)"CrossSection", (char*)"Mb") &&
      !check_sdds_column(&SDDSin, (char*)"CrossSection", (char*)"MBarns") && 
      !check_sdds_column(&SDDSin, (char*)"CrossSection", (char*)"megabarns"))
    bombElegantVA((char*)"Column \"CrossSection\" is missing, not floating-point type, or does not have units of megabarns (or Mbarns or Mb)",
                  filename);

  if (SDDS_ReadPage(&SDDSin)<=0) 
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
  
  if ((ionProperties.nSpecies = SDDS_RowCount(&SDDSin))<=0)
    bombElegantVA((char*)"Ion properties file %s appears to have no rows.\n", filename);

  if (!(ionProperties.ionName = (char**)SDDS_GetColumn(&SDDSin, (char*)"IonName")) ||
      !(ionProperties.mass = SDDS_GetColumnInDoubles(&SDDSin, (char*)"Mass")) ||
      !(ionProperties.chargeState = SDDS_GetColumnInDoubles(&SDDSin, (char*)"ChargeState")) ||
      !(ionProperties.crossSection = SDDS_GetColumnInDoubles(&SDDSin, (char*)"CrossSection")) ||
      !(sourceName = (char**)SDDS_GetColumn(&SDDSin, (char*)"SourceName")))
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);

  if (!(ionProperties.sourceGasIndex = (long*)tmalloc(sizeof(*(ionProperties.sourceGasIndex))*ionProperties.nSpecies)) ||
      !(ionProperties.sourceIonIndex = (long*)tmalloc(sizeof(*(ionProperties.sourceIonIndex))*ionProperties.nSpecies)))
    bombElegantVA((char*)"Memory allocation failure allocating arrays for %ld ion species.\n", ionProperties.nSpecies);

  /* Figure out the source gas or source ion indices */
  for (i=0; i<ionProperties.nSpecies; i++) {
    ionProperties.sourceGasIndex[i] = ionProperties.sourceIonIndex[i] = -1;
    if (ionProperties.chargeState[i]<=0) 
      bombElegantVA((char*)"Ion %s has non-positive charge state", ionProperties.ionName[i]);
    if (ionProperties.chargeState[i]==1) {
      if ((ionProperties.sourceGasIndex[i] = match_string(sourceName[i], pressureData.gasName, pressureData.nGasses, EXACT_MATCH))<0) {
	if ((ionProperties.sourceIonIndex[i] = match_string(sourceName[i], ionProperties.ionName, ionProperties.nSpecies, EXACT_MATCH))<0) {
	  bombElegantVA((char*)"Unable to find match to gas source \"%s\" for species \"%s\"", 
			sourceName[i], ionProperties.ionName[i]);
	}
      }
    } else {
      if ((ionProperties.sourceIonIndex[i] = match_string(sourceName[i], ionProperties.ionName, ionProperties.nSpecies, EXACT_MATCH))<0) 
        bombElegantVA((char*)"Unable to find match to ion source \"%s\" for species \"%s\"", 
                      sourceName[i], ionProperties.ionName[i]);
    }
  }

  printf("Finished reading ion properties file %s\n", filename);
  fflush(stdout);
}


void completeIonEffectsSetup(RUN *run, LINE_LIST *beamline)
{
  /* scan through the beamline, find IONEFFECTS elements, set parameters */
  ELEMENT_LIST *eptr, *eptrLast;
  IONEFFECTS *ionEffects;
  short chargeSeen = 0;
  long iPlane;
  
  eptr = &(beamline->elem);
  eptrLast = eptr;
  if (eptr->type == T_IONEFFECTS) 
    bombElegant("ION_EFFECTS element cannot be the first element in the beamline", NULL);
  nIonEffectsElements = 0;
  sStartFirst = 0;
  while (eptr) {
    if (eptr->type == T_CHARGE)
      chargeSeen = 1;
    if (eptr->type == T_IONEFFECTS) {
      if (!chargeSeen)
        bombElegant("ION_EFFECTS element preceeds the CHARGE element", NULL);
      /* Set the start of the s range for this element */
      ionEffects = (IONEFFECTS*)eptr->p_elem;
      ionEffects->sLocation = eptr->end_pos;
      ionEffects->sStart = (eptrLast->end_pos + eptr->end_pos)/2;
      /* in case this is the last element in the beamline, set s so the range ends here */
      ionEffects->sEnd = eptr->end_pos;
      if (eptrLast && eptrLast->type == T_IONEFFECTS) {
        /* set the s range for the previous ion effects element */
        ionEffects = (IONEFFECTS*)eptrLast->p_elem;
        ionEffects->sEnd = (eptrLast->end_pos + eptr->end_pos)/2;
      }
      eptrLast = eptr;
      if (nIonEffectsElements==0)
        sStartFirst = ionEffects->sStart;
      for (iPlane=0; iPlane<2; iPlane++)
        ionEffects->xyIonHistogram[iPlane] = ionEffects->ionHistogram[iPlane] = 
          ionEffects->ionHistogramFit[iPlane] = NULL;
      nIonEffectsElements++;
    } 
    eptr = eptr->succ;
  }

  eptr = &(beamline->elem);
  firstIonEffects = NULL;
  while (eptr) {
    if (eptr->type == T_IONEFFECTS) {
      ionEffects = (IONEFFECTS*)eptr->p_elem;
      if (!firstIonEffects)
        firstIonEffects = ionEffects;
      if (verbosity>10) 
        printf("IONEFFECTS element %s#%ld at s=%le m spans s:[%le, %le] m\n",
               eptr->name, eptr->occurence, eptr->end_pos, ionEffects->sStart, ionEffects->sEnd);

      /* Determine the average pressure for each gas */
      ionEffects->pressure = (double*)tmalloc(sizeof(*(ionEffects->pressure))*pressureData.nGasses);
      computeAverageGasPressures(ionEffects->sStart, ionEffects->sEnd, ionEffects->pressure, &pressureData);
      if (verbosity>20) {
        long i;
        printf("Average pressures over s:[%le, %le] m\n", ionEffects->sStart, ionEffects->sEnd);
        for (i=0; i<pressureData.nGasses; i++)
          printf("%s:%.2f nT  ", pressureData.gasName[i], ionEffects->pressure[i]*1e9);
        printf("\n");
        fflush(stdout);
      }

      /* Allocate arrays (should really clean up the arrays first in case this is one in a series of runs) */
      ionEffects->coordinate = (double***)calloc(ionProperties.nSpecies, sizeof(*(ionEffects->coordinate)));
      ionEffects->nIons = (long*)calloc(ionProperties.nSpecies, sizeof(*(ionEffects->nIons)));
      ionEffects->t = 0;
      if (ionEffects->macroIons<=0)
        ionEffects->macroIons = macro_ions;
      if (ionEffects->generationInterval<=0)
        ionEffects->generationInterval = generation_interval;
      for (iPlane=0; iPlane<2; iPlane++) {
        if (ionEffects->span[iPlane]<=0)
          ionEffects->span[iPlane] = ion_span[iPlane];
        if (ionEffects->binDivisor[iPlane]<=0)
          ionEffects->binDivisor[iPlane] = ion_bin_divisor[iPlane];
	ionEffects->rangeMultiplier[iPlane] = ion_range_multiplier[iPlane];
	ionEffects->sigmaLimitMultiplier[iPlane] = ion_sigma_limit_multiplier[iPlane];
      }
    }
    eptr = eptr->succ;
  }

  ionsInitialized = 1;

}

void trackWithIonEffects
(
 double **part0,          /* part0[i][j] is the jth coordinate (x,x',y,y',t,delta) for the ith particle */
 long np0,                /* number of particles (on this processor) */
 IONEFFECTS *ionEffects,  /* ion effects element data */
 double Po,               /* central momentum (beta*gamma) */
 long iPass,              /* pass number */ 
 long nPasses,            /* number of passes */
 CHARGE *charge           /* beam charge structure */
 )
{
  long ip, iSpecies, iIon;
  long iBunch, nBunches=0;
  double *time0 = NULL;          /* array to record arrival time of each particle */
  double **part = NULL;          /* particle buffer for working bunch */
  double *time = NULL;           /* array to record arrival time of each particle in working bunch */
  long *ibParticle = NULL;       /* array to record which bunch each particle is in */
  long **ipBunch = NULL;         /* array to record particle indices in part0 array for all particles in each bunch */
  long *npBunch = NULL;          /* array to record how many particles are in each bunch */
  long np, npTotal, max_np = 0;
  /* properties of the electron beam */
  double centroid[2], sigma[2], tNow, qBunch;
  //double sigmatemp[2];
  /* properties of the ion cloud */
  double ionCentroid[2], ionSigma[2], qIon;
  double **speciesCentroid=NULL, *speciesCharge=NULL, **speciesSigma=NULL;
  long *speciesCount=NULL;
  double unitsFactor; /* converts Torr to 1/m^3 and mBarns to m^2 */
#if USE_MPI
  MPI_Status mpiStatus;
#endif

  if (verbosity>30) {
    printf("Running ION_EFFECTS\n");
    fflush(stdout);
  }

  if (!ionsInitialized) {
    bombElegant("IONEFFECTS element seen, but ion_effects command was not given to initialize ion modeling.", NULL);
  }
  if (ionEffects->disable) {
    if (verbosity>30) {
      printf("ION_EFFECTS disabled, returning\n");
      fflush(stdout);
    }
    return ;
  }

  if (ionEffects->startPass<0)
    ionEffects->startPass = 0;
  if (iPass==0)
    ionEffects->xyFitSet[0] = ionEffects->xyFitSet[1] = 0;

  if ((ionEffects->startPass>=0 && iPass<ionEffects->startPass) ||
      (ionEffects->endPass>=0 && iPass>ionEffects->endPass) ||
      (ionEffects->passInterval>=1 && (iPass-ionEffects->startPass)%ionEffects->passInterval!=0)) {
    if (verbosity>30) {
      printf("ION_EFFECTS out of pass range (pass:%ld, start:%ld, end:%ld, interval:%ld), returning\n",
	     iPass, ionEffects->startPass, ionEffects->endPass, ionEffects->passInterval);
      fflush(stdout);
    }
    return;
  }
  if (verbosity>40) {
    printf("ION_EFFECTS within range (pass:%ld, start:%ld, end:%ld, interval:%ld), returning\n",
	   iPass, ionEffects->startPass, ionEffects->endPass, ionEffects->passInterval);
    fflush(stdout);
  }

  /* converts Torr to 1/m^3 and mBarns to m^2 */
  unitsFactor = 1e-22/(7.5006e-3*k_boltzmann_mks*pressureData.temperature);
    
  if (isSlave || !notSinglePart) {
    /* Determine which bunch each particle is in */
    determine_bucket_assignments(part0, np0, charge?charge->idSlotsPerBunch:0, Po, &time0, &ibParticle, &ipBunch, &npBunch, &nBunches, -1);
#if USE_MPI
    if (mpiAbort)
      return;
#endif
  }

#if USE_MPI
  /* Share the number of bunches with the master node */
  MPI_Barrier(MPI_COMM_WORLD);
  if (myid==1) 
    MPI_Send(&nBunches, 1, MPI_LONG, 0, 1, MPI_COMM_WORLD);
  if (myid==0) 
    MPI_Recv(&nBunches, 1, MPI_LONG, 1, 1, MPI_COMM_WORLD, &mpiStatus);
#endif
  if (verbosity>30) {
    printf("Running ION_EFFECTS with %ld bunches\n", nBunches);
    fflush(stdout);
  }
    
#if USE_MPI
  if (myid==0) {
#endif
    if (iPass==0 && ionEffects->sStart==sStartFirst) {
      iIonEffectsElement = 0;
      if (SDDS_beamOutput) {
	if (verbosity>10) {
	  printf("Starting page (%ld rows) for ion-related electron beam output\n", 
		 nPasses*(beam_output_all_locations?nIonEffectsElements:1)*nBunches);
	  fflush(stdout);
	}
        if (!SDDS_StartPage(SDDS_beamOutput, nPasses*(beam_output_all_locations?nIonEffectsElements:1)*nBunches)) {
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
          exitElegant(1);
        }
        iBeamOutput = 0;
      }
      if (SDDS_ionDensityOutput) {
	if (verbosity>10) {
	  printf("Starting page (%ld rows) for ion density output\n", 
		 nPasses*(beam_output_all_locations?nIonEffectsElements:1)*nBunches);
	  fflush(stdout);
	}
        if (!SDDS_StartPage(SDDS_ionDensityOutput, nPasses*(ion_output_all_locations?nIonEffectsElements:1)*nBunches)) {
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
          exitElegant(1);
        }
        iIonDensityOutput = 0;
      }
    }
#if USE_MPI
  }
#endif

  for (iBunch=0; iBunch<nBunches; iBunch++) {
    np = 0;
    /* Loop over all bunches */
    if (verbosity>20) {
      printf("Working on bunch %ld\n", iBunch);
      fflush(stdout);
    }
    if (isSlave || !notSinglePart) {
      if (nBunches==1) {
        time = time0;
        part = part0;
        np = np0;
      } else {
        if (npBunch)
          np = npBunch[iBunch];
        else 
          np = 0;
        if (np && (!ibParticle || !ipBunch || !time0)) {
#if USE_MPI
          mpiAbort = MPI_ABORT_BUCKET_ASSIGNMENT_ERROR;
          return;
#else
          printf("Problem in determine_bucket_assignments. Seek professional help.\n");
          exitElegant(1);
#endif
        }
      
	if (np>max_np) {
	  if (part)
	    free_czarray_2d((void**)part, max_np, COORDINATES_PER_PARTICLE);
	  part = (double**)czarray_2d(sizeof(double), np, COORDINATES_PER_PARTICLE);
	  time = (double*)trealloc(time, sizeof(*time)*np);
	  max_np = np;
	}
	if (np>0) {
	  for (ip=0; ip<np; ip++) {
	    time[ip] = time0[ipBunch[iBunch][ip]];
	    memcpy(part[ip], part0[ipBunch[iBunch][ip]], sizeof(double)*COORDINATES_PER_PARTICLE);
	  }
	}
      }
    }

#if USE_MPI
    if (myid==0)
      np = 0;
#endif

    /* Compute Sx, Cx, Sy, Cy, <t> (in parallel over all cores if needed) */


#if USE_MPI
    sigma[0] = computeRmsCoordinate_p(part, 0, np, centroid+0, &npTotal, MPALGORITHM);
    sigma[1] = computeRmsCoordinate_p(part, 2, np, centroid+1, &npTotal, MPALGORITHM);
    tNow = computeAverage_p(time, np, MPI_COMM_WORLD);
#else
    sigma[0] = computeRmsCoordinate(part, 0, np, centroid+0, &npTotal);
    sigma[1] = computeRmsCoordinate(part, 2, np, centroid+1, NULL);
    compute_average(&tNow, time, np);
#endif
    qBunch = npTotal*charge->macroParticleCharge;

#if DEBUG
    if ((sigma[0] > 0.01) || (sigma[1] > 0.01)) {
      printf("beam sigma too large: bunch %ld, Pass %ld, s=%f \n", iBunch, iPass, ionEffects->sLocation);
    }
#endif


    if (verbosity>30) {
      printf("np: %ld, <t>: %le, sigma x,y: %le, %le,  centroid x,y: %le, %le,  q: %le\n",
             npTotal, tNow, sigma[0], sigma[1], centroid[0], centroid[1], qBunch);
      fflush(stdout);
    }

#ifdef DEBUG
    if ((verbosity > 4) &&  (ionEffects->sStart < 10)) {
      FILE * fbeam;
      fbeam = fopen("beam_info.dat", "a");
      fprintf(fbeam, "%le  %le  %le  %le  %le  %le \n",
              tNow, sigma[0], sigma[1], centroid[0], centroid[1], qBunch);
      fclose(fbeam);
    }
#endif
#if USE_MPI
    if (myid==0) {
#endif
      if (SDDS_beamOutput) {
        if ((beam_output_all_locations || ionEffects==firstIonEffects) &&
            !SDDS_SetRowValues(SDDS_beamOutput, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, iBeamOutput++,
                               "t", tNow, "Pass", iPass,
                               "Bunch", iBunch, "qBunch", qBunch, "npBunch", npTotal,
                               "s", ionEffects->sLocation,
                               "Sx", sigma[0], "Sy", sigma[1], "Cx", centroid[0], "Cy", centroid[1],
                               NULL)) {
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
          exitElegant(1);
        }
      }
#if USE_MPI
    }
#endif

    if (isSlave || !notSinglePart) {
      /*** Advance the ion positions */
      for (iSpecies=0; iSpecies<ionProperties.nSpecies; iSpecies++) {
        for (iIon=0; iIon<ionEffects->nIons[iSpecies]; iIon++) {
          ionEffects->coordinate[iSpecies][iIon][0] += (tNow-ionEffects->t)*ionEffects->coordinate[iSpecies][iIon][1];
          ionEffects->coordinate[iSpecies][iIon][2] += (tNow-ionEffects->t)*ionEffects->coordinate[iSpecies][iIon][3];
        }
      }
      ionEffects->t = tNow;
      
      if (ionEffects->span[0] || ionEffects->span[1]) {
        /*** Eliminate ions that are outside the simulation region */
        for (iSpecies=0; iSpecies<ionProperties.nSpecies; iSpecies++) {
          for (iIon=0; iIon<ionEffects->nIons[iSpecies]; iIon++) {
            if ((ionEffects->span[0] && fabs(ionEffects->coordinate[iSpecies][iIon][0])>ionEffects->span[0]) ||
                (ionEffects->span[1] && fabs(ionEffects->coordinate[iSpecies][iIon][2])>ionEffects->span[1])) {
              if (ionEffects->nIons[iSpecies]>1) {
                /* Move ion at the top of the array into this slot in the array */
                long k;
                for (k=0; k<5; k++) {
                  ionEffects->coordinate[iSpecies][iIon][k] = ionEffects->coordinate[iSpecies][ionEffects->nIons[iSpecies]-1][k];
#if DEBUG
                  if (isnan(ionEffects->coordinate[iSpecies][iIon][k])) {
                    printf("nan");
                  }
#endif
		}
              }
              ionEffects->nIons[iSpecies] -= 1;
              iIon -= 1;
            }
          }
        }
      }
      

      
      if (((iPass-ionEffects->startPass)*nBunches+iBunch)%ionEffects->generationInterval==0) {
        /*** Generate ions */
        for (iSpecies=0; iSpecies<ionProperties.nSpecies; iSpecies++) {
          long nToAdd = 0, index;
          double qToAdd = 0;
          if ((index=ionProperties.sourceGasIndex[iSpecies])>=0) {
            /* this is a singly-ionized molecule, so use source gas 
               nToAdd =  someFunctionOfPressure(ionEffects->pressure[index], ...);
            */
#if USE_MPI
            /* The macroIons parameter is the number for all processors, so we need to 
             * apportion the ions among the working processors 
             */
            nToAdd = ionEffects->macroIons/(n_processors-1.0);
            long nLeft = ionEffects->macroIons - nToAdd*(n_processors-1);
            for (long iLeft=0; iLeft<nLeft; iLeft++) {
              if (leftIonCounter%(n_processors-1)==(myid-1))
                nToAdd ++;
              leftIonCounter++; /* This counter will be the same on all processors */
            }
#else
            nToAdd = ionEffects->macroIons;
#endif

            if (nToAdd) {
              qToAdd = unitsFactor * qBunch * ionEffects->pressure[index] * ionEffects->generationInterval * \
                ionProperties.crossSection[iSpecies] * (ionEffects->sEnd - ionEffects->sStart) / ionEffects->macroIons;
              if (symmetrize) {
                nToAdd *= 2;
                qToAdd /= 2;
              }
              addIons(ionEffects, iSpecies, nToAdd, qToAdd, centroid, sigma, symmetrize);
            }
          } else if (((index=ionProperties.sourceIonIndex[iSpecies])>=0) && \
		     (((iPass-ionEffects->startPass)*nBunches+iBunch)%multiple_ionization_interval == 0)) {
            /* This is a multiply-ionized molecule, so use source ion density.
             * Relevant quantities:
             * ionEffects->nIons[index] --- Number of ions of the source species
             * ionEffects->coordinate[index][j][k] --- kth coordinate of jth source ion 
             * ionProperties.crossSection[iSpecies] --- Cross section for producing new ion from the source ions
             */
            /* 
               nToAdd = someFunctionOfExistingNumberOfIons(...); 
            */


	    double beamFact = 0, jx = 0, jy = 0, Pmi = 0, rnd = 0;
	    beamFact = multiple_ionization_interval * 1e-22 * qBunch / e_mks / (2*PI * sigma[0] * sigma[1]);
	    for (int jMacro = 0; jMacro < ionEffects->nIons[index]; jMacro++) {
	      jx = ionEffects->coordinate[index][jMacro][0] - centroid[0];
	      jy = ionEffects->coordinate[index][jMacro][2] - centroid[1];
	      Pmi = beamFact * ionProperties.crossSection[iSpecies] * \
		exp(-sqr(jx) / (2*sqr(sigma[0])) - sqr(jy) / (2*sqr(sigma[1])));

	      rnd = random_2(0);
	      if (rnd < Pmi) { //multiple ionization occurs
		double qToAdd, mx, my;
		qToAdd = ionEffects->coordinate[index][jMacro][4];
		mx = ionEffects->coordinate[index][jMacro][0];
		my = ionEffects->coordinate[index][jMacro][2];
		addIon_point(ionEffects, iSpecies, qToAdd, mx, my); //add multiply ionized ion

		// Initial kinetic energy
		double vmag, ionMass, vx, vy, rangle, Emi;	      
		ionMass = 1.672621898e-27 * ionProperties.mass[iSpecies]; 
		Emi = fabs(gauss_rn_lim(20, 10, 3, random_4));
		//Emi = 0;
		vmag = sqrt(2 * Emi * e_mks / ionMass);
		rangle = random_2(0) * 2 * PI;
		vx = vmag * cos(rangle);
		vy = vmag * sin(rangle);
		ionEffects->coordinate[iSpecies][ionEffects->nIons[iSpecies]-1][1] = vx;
		ionEffects->coordinate[iSpecies][ionEffects->nIons[iSpecies]-1][3] = vy;

		//delete source ion
		long k;
                for (k=0; k<5; k++) {
                  ionEffects->coordinate[index][jMacro][k] = ionEffects->coordinate[index][ionEffects->nIons[index]-1][k];
		}
		ionEffects->nIons[index] -= 1;
		jMacro -= 1;
              }
              
	    } 
	    
            //bombElegant("Multiple ionization not implemented at this time", NULL);
          }
        }
      }


     

          
      /*** Determine and apply kicks from beam to ions */
      for (iSpecies=0; iSpecies<ionProperties.nSpecies; iSpecies++) {
        double ionMass, ionCharge, *coord, kick[2]={0,0};
      
        /* Relevant quantities:
         * ionProperties.chargeState[iSpecies] --- Charge state of the ion (integer)
         * ionProperties.mass[iSpecies] --- Mass of the ion (AMUs)
         * qBunch --- bunch charge (C)
         * sigma[0], sigma[1] --- x, y sigma (m)
         * centroid[0], centroid[1] --- x , y centroid (m)
         * ionEffects->nIons[index] --- Number of ions of the source species
         * ionEffects->coordinate[index][j][k] --- kth coordinate of jth source ion 
         */
        
	ionMass = 1.672621898e-27 * ionProperties.mass[iSpecies]; 
	ionCharge = (double)ionProperties.chargeState[iSpecies];


	double tempkick[2], maxkick[2], tempart[4];
	tempart[0] = sigma[0] + centroid[0];
	tempart[2] = 0;
	gaussianBeamKick(tempart, centroid, sigma, tempkick, qBunch, ionMass, ionCharge);
	maxkick[0] = 2*abs(tempkick[0]);
	
	tempart[2] = sigma[1] + centroid[1];
	tempart[0] = 0;
	gaussianBeamKick(tempart, centroid, sigma, tempkick, qBunch, ionMass, ionCharge);
	maxkick[1] = 2*abs(tempkick[1]);
        
        for (iIon=0; iIon<ionEffects->nIons[iSpecies]; iIon++) {
          /*
            Get effective fields in V/m
            BeamFieldFunction(qBunch, sigma, centroid, &(ionEffects->coordinate[iSpecies][iIon][0]), &Ex, &Ey);
          */
	  coord = ionEffects->coordinate[iSpecies][iIon];
          
	  gaussianBeamKick(coord, centroid, sigma, kick, qBunch, ionMass, ionCharge);


          
#if DEBUG
	  if (isnan(kick[0])) {
	    printf("kick nan");
	  }
#endif

	  if (abs(kick[0]) < maxkick[0] && abs(kick[1]) < maxkick[1]) {
	    ionEffects->coordinate[iSpecies][iIon][1] += kick[0];
	    ionEffects->coordinate[iSpecies][iIon][3] += kick[1];
	  } 

	  //          ionEffects->coordinate[iSpecies][iIon][1] += 
          //  Ex*ionProperties.chargeState[iSpecies]/ionProperties.mass[iSpecies]; /* plus other constants */
          //ionEffects->coordinate[iSpecies][iIon][3] += Ey*ionProperties.chargeState[iSpecies]/ionProperties.mass[iSpecies]; /* plus other constants */
	  
        }
      }
    }  /* branch for workers */
    
    /* Master should also have valid values for these since we may want to log them */
    ionCentroid[0] = 0;
    ionCentroid[1] = 0;
    ionSigma[0] = 0;
    ionSigma[1] = 0;
    qIon = 0;
    long mTot = 0;
    double bx1, bx2, by1, by2;
    bx1 = centroid[0] - 3*sigma[0];
    bx2 = centroid[0] + 3*sigma[0];
    by1 = centroid[1] - 3*sigma[1];
    by2 = centroid[1] + 3*sigma[1];


    if (ion_species_output && iBunch==0) {
      speciesCentroid = (double**)zarray_2d(sizeof(double), ionProperties.nSpecies, 2);
      speciesSigma = (double**)zarray_2d(sizeof(double), ionProperties.nSpecies, 2);
      speciesCharge = (double*)tmalloc(sizeof(*speciesCharge)*ionProperties.nSpecies);
      speciesCount = (long*)tmalloc(sizeof(*speciesCount)*ionProperties.nSpecies);
    }
    
    /*** Compute field due to ions */
    if (isSlave || !notSinglePart) {
      
      /* Compute charge-weighted centroids */
      for (iSpecies=0; iSpecies<ionProperties.nSpecies; iSpecies++) {
	/* Relevant quantities:
	 * ionProperties.chargeState[iSpecies] --- Charge state of the ion (integer)
	 * ionEffects->nIons[index] --- Number of ions of the source species
	 * ionEffects->coordinate[index][j][k] --- kth coordinate of jth source ion 
	 */
	
	int jMacro = 0;
	if (ion_species_output) {
	  speciesCentroid[iSpecies][0] = speciesCentroid[iSpecies][1] = 0;
	  speciesCharge[iSpecies] = 0;
	  speciesCount[iSpecies] = 0;
	}
	for (jMacro=0; jMacro < ionEffects->nIons[iSpecies]; jMacro++) {
	  if ((ionEffects->coordinate[iSpecies][jMacro][0] > bx1) && (ionEffects->coordinate[iSpecies][jMacro][0] < bx2) &&
	      (ionEffects->coordinate[iSpecies][jMacro][2] > by1) && (ionEffects->coordinate[iSpecies][jMacro][2] < by2)) {
	    if (ion_species_output) {
	      speciesCentroid[iSpecies][0] += ionEffects->coordinate[iSpecies][jMacro][0]*ionEffects->coordinate[iSpecies][jMacro][4];
	      speciesCentroid[iSpecies][1] += ionEffects->coordinate[iSpecies][jMacro][2]*ionEffects->coordinate[iSpecies][jMacro][4];
	      speciesCharge[iSpecies] += ionEffects->coordinate[iSpecies][jMacro][4];
	      speciesCount[iSpecies] += 1;
	    }
	    ionCentroid[0] += ionEffects->coordinate[iSpecies][jMacro][0]*ionEffects->coordinate[iSpecies][jMacro][4];
	    ionCentroid[1] += ionEffects->coordinate[iSpecies][jMacro][2]*ionEffects->coordinate[iSpecies][jMacro][4];
	    qIon += ionEffects->coordinate[iSpecies][jMacro][4];
	    mTot++;
	  }
	}
      }
    } else {
      for (iSpecies=0; iSpecies<ionProperties.nSpecies; iSpecies++) {
	if (ion_species_output) {
	  speciesCentroid[iSpecies][0] = speciesCentroid[iSpecies][1] = 0;
	  speciesCharge[iSpecies] = 0;
	  speciesCount[iSpecies] = 0;
	}
      }
    }
    
    long mTotTotal = 0;
#if USE_MPI
    /* Sum ion centroid and charge data over all nodes */
    long mTotMin, mTotMax;
    double qIonTotal, ionCentroidTotal[2];
    /* Use buffers to perform reduce operations to improve efficiency */
    double inBuffer[4], sumBuffer[4];
    inBuffer[0] = mTot;
    inBuffer[1] = qIon;
    inBuffer[2] = ionCentroid[0];
    inBuffer[3] = ionCentroid[1];
    MPI_Allreduce(inBuffer, sumBuffer, 4, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    mTotTotal = sumBuffer[0];
    qIonTotal = sumBuffer[1];
    ionCentroidTotal[0] = sumBuffer[2];
    ionCentroidTotal[1] = sumBuffer[3];
    if (myid==0)
      mTot = LONG_MAX;
    MPI_Allreduce(&mTot, &mTotMin, 1, MPI_LONG, MPI_MIN, MPI_COMM_WORLD);
    if (myid==0)
      mTot = LONG_MIN;
    MPI_Allreduce(&mTot, &mTotMax, 1, MPI_LONG, MPI_MAX, MPI_COMM_WORLD);
    qIon = qIonTotal;
    if (qIon!=0) {
      ionCentroid[0] = ionCentroidTotal[0]/qIon;
      ionCentroid[1] = ionCentroidTotal[1]/qIon;
    }
    if (ion_species_output) {
      for (iSpecies=0; iSpecies<ionProperties.nSpecies; iSpecies++) {
        if (myid!=0) {
          inBuffer[0] = speciesCentroid[iSpecies][0];
          inBuffer[1] = speciesCentroid[iSpecies][1];
          inBuffer[2] = speciesCharge[iSpecies];
          inBuffer[3] = speciesCount[iSpecies];
        } else {
          inBuffer[0] = inBuffer[1] = inBuffer[2] = inBuffer[3] = 0;
        }
        MPI_Allreduce(inBuffer, sumBuffer, 4, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        speciesCentroid[iSpecies][0] = sumBuffer[0];
        speciesCentroid[iSpecies][1] = sumBuffer[1];
        speciesCharge[iSpecies] = sumBuffer[2];
        speciesCount[iSpecies] = sumBuffer[3];
        if (speciesCharge[iSpecies]) {
          speciesCentroid[iSpecies][0] /= speciesCharge[iSpecies];
          speciesCentroid[iSpecies][1] /= speciesCharge[iSpecies];
        }
      }
    }
#else
    if (qIon) {
      ionCentroid[0] = ionCentroid[0]/qIon;
      ionCentroid[1] = ionCentroid[1]/qIon;
      mTotTotal = mTot;
      if (ion_species_output) {
        for (iSpecies=0; iSpecies<ionProperties.nSpecies; iSpecies++) {
          if (speciesCharge[iSpecies]) {
            speciesCentroid[iSpecies][0] /= speciesCharge[iSpecies];
            speciesCentroid[iSpecies][1] /= speciesCharge[iSpecies];
          }
        }
      }
    }
#endif

    if (isSlave || !notSinglePart) {
      /* Compute charge-weighted rms size */
      for (iSpecies=0; iSpecies<ionProperties.nSpecies; iSpecies++) {
        /* Relevant quantities:
         * ionProperties.chargeState[iSpecies] --- Charge state of the ion (integer)
         * ionEffects->nIons[index] --- Number of ions of the source species
         * ionEffects->coordinate[index][j][k] --- kth coordinate of jth source ion 
         */

	int jMacro = 0;
        if (ion_species_output) 
          speciesSigma[iSpecies][0] = speciesSigma[iSpecies][1] = 0;
	for (jMacro=0; jMacro < ionEffects->nIons[iSpecies]; jMacro++) {
	  if ((ionEffects->coordinate[iSpecies][jMacro][0] > bx1) && (ionEffects->coordinate[iSpecies][jMacro][0] < bx2) &&
	      (ionEffects->coordinate[iSpecies][jMacro][2] > by1) && (ionEffects->coordinate[iSpecies][jMacro][2] < by2)) {
	    ionSigma[0] += sqr(ionEffects->coordinate[iSpecies][jMacro][0]-ionCentroid[0])*ionEffects->coordinate[iSpecies][jMacro][4];
	    ionSigma[1] += sqr(ionEffects->coordinate[iSpecies][jMacro][2]-ionCentroid[1])*ionEffects->coordinate[iSpecies][jMacro][4];
	    if (ion_species_output) {
	      speciesSigma[iSpecies][0] += sqr(ionEffects->coordinate[iSpecies][jMacro][0]-speciesCentroid[iSpecies][0])*ionEffects->coordinate[iSpecies][jMacro][4];
	      speciesSigma[iSpecies][1] += sqr(ionEffects->coordinate[iSpecies][jMacro][2]-speciesCentroid[iSpecies][1])*ionEffects->coordinate[iSpecies][jMacro][4];
	    }
	  }
	}
      }
    } else {
      for (iSpecies=0; iSpecies<ionProperties.nSpecies; iSpecies++) {
        if (ion_species_output) 
          speciesSigma[iSpecies][0] = speciesSigma[iSpecies][1] = 0;
      }
    }
    
#if USE_MPI
    double ionSigmaTotal[2];
    MPI_Allreduce(ionSigma, ionSigmaTotal, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    if (qIon) {
      ionSigma[0] = sqrt(ionSigmaTotal[0]/qIon);
      ionSigma[1] = sqrt(ionSigmaTotal[1]/qIon);
    }
    if (ion_species_output) {
      for (iSpecies=0; iSpecies<ionProperties.nSpecies; iSpecies++) {
        MPI_Allreduce(speciesSigma[iSpecies], ionSigmaTotal, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        if (speciesCharge[iSpecies]) {
          speciesSigma[iSpecies][0] = sqrt(ionSigmaTotal[0]/speciesCharge[iSpecies]);
          speciesSigma[iSpecies][1] = sqrt(ionSigmaTotal[1]/speciesCharge[iSpecies]);
        }
      }
    }
#else
    if (qIon) {
      ionSigma[0] = sqrt(ionSigma[0]/qIon);
      ionSigma[1] = sqrt(ionSigma[1]/qIon);
    }
    if (ion_species_output) {
      for (iSpecies=0; iSpecies<ionProperties.nSpecies; iSpecies++) {
        if (speciesCharge[iSpecies]) {
          speciesSigma[iSpecies][0] = sqrt(speciesSigma[iSpecies][0]/speciesCharge[iSpecies]);
          speciesSigma[iSpecies][1] = sqrt(speciesSigma[iSpecies][1]/speciesCharge[iSpecies]);
        }
      }
    }
#endif

#if DEBUG
    if ((ionSigma[0] > 0.01) || (ionSigma[1] > 0.01)) {
      //      printf("ion sigma too large");  
      printf("ion sigma too large (sx=%e, sy=%e): bunch %d, Pass %d, s=%f \n", ionSigma[0], ionSigma[1], iBunch, iPass, ionEffects->sLocation);
    }
#endif


#if USE_MPI
    if (myid==0) {
#endif
      if ((SDDS_ionDensityOutput) && (iPass-ionEffects->startPass+iBunch)%ionEffects->generationInterval==0) {
      if (ion_output_all_locations || ionEffects==firstIonEffects) {
        if (!SDDS_SetRowValues(SDDS_ionDensityOutput, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, iIonDensityOutput,
                               "t", tNow, "Pass", iPass, "Bunch", iBunch,  "s", ionEffects->sLocation,
                               "qIons", qIon, "Sx", ionSigma[0], "Sy", ionSigma[1],
                               "Cx", ionCentroid[0], "Cy", ionCentroid[1], "nMacroIons", mTotTotal,
#if USE_MPI
                               "nMacroIonsMin", mTotMin, 
                               "nMacroIonsMax", mTotMax, 
#endif
                               NULL)) {
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
          exitElegant(1);
        }
        if (ion_species_output) {
          for (iSpecies=0; iSpecies<ionProperties.nSpecies; iSpecies++) {
            if (!SDDS_SetRowValues(SDDS_ionDensityOutput, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, iIonDensityOutput,
                                   makeSpeciesName("qIons_", ionProperties.ionName[iSpecies]),
                                   speciesCharge[iSpecies], NULL) ||
                !SDDS_SetRowValues(SDDS_ionDensityOutput, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, iIonDensityOutput,
                                   makeSpeciesName("nMacroIons_", ionProperties.ionName[iSpecies]),
                                   speciesCount[iSpecies], NULL) ||
                !SDDS_SetRowValues(SDDS_ionDensityOutput, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, iIonDensityOutput,
                                   makeSpeciesName("Cx_", ionProperties.ionName[iSpecies]),
                                   speciesCentroid[iSpecies][0], NULL) ||
                !SDDS_SetRowValues(SDDS_ionDensityOutput, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, iIonDensityOutput,
                                   makeSpeciesName("Cy_", ionProperties.ionName[iSpecies]),
                                   speciesCentroid[iSpecies][1], NULL) ||
                !SDDS_SetRowValues(SDDS_ionDensityOutput, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, iIonDensityOutput,
                                   makeSpeciesName("Sx_", ionProperties.ionName[iSpecies]),
                                   speciesSigma[iSpecies][0], NULL) ||
                !SDDS_SetRowValues(SDDS_ionDensityOutput, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, iIonDensityOutput,
                                   makeSpeciesName("Sy_", ionProperties.ionName[iSpecies]),
                                   speciesSigma[iSpecies][1], NULL)) {
              SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
              exitElegant(1);
            }
          }
        }
        iIonDensityOutput++;
      }
    }
#if USE_MPI
  }
#endif
    
    //for multi-gaussian and multi-lorentzian kick
    double paramValueX[9], paramValueY[9];
    double tempCentroid[9][2], tempSigma[9][2], tempkick[2], tempQ[9];
    double normX, normY;
    
    ionEffects->ionFieldMethod = ionFieldMethod;

    if (ionFieldMethod != ION_FIELD_GAUSSIAN) {
      makeIonHistograms(ionEffects, ionProperties.nSpecies, centroid, sigma, ionCentroid, ionSigma);
#if USE_MPI && MPI_DEBUG
      printf("Running bigaussian/bilorentzian fit for s=%le, bunch %ld\n", ionEffects->sLocation, iBunch);
      printf("Ion centroid = (%le, %le), sigma = (%le, %le)\n", ionCentroid[0], ionCentroid[1], 
	     ionSigma[0], ionSigma[1]);
      fflush(stdout);
#endif
      /* We take a return value here for future improvement in which the fitting function is automatically selected. */
      ionEffects->ionFieldMethod = 
        multipleWhateverFit(sigma, centroid, paramValueX, paramValueY, ionEffects, ionSigma, ionCentroid);
      
      /* these factors needed because we fit charge histograms instead of charge densities */
      normX = ionEffects->ionRange[0] / ionEffects->ionBins[0];
      normY = ionEffects->ionRange[1] / ionEffects->ionBins[1];
      
      if (ionEffects->ionFieldMethod==ION_FIELD_BIGAUSSIAN || ionEffects->ionFieldMethod==ION_FIELD_TRIGAUSSIAN) {
        /* paramValueX[0..8] = sigma1, centroid1, height1, sigma2, centroid2, height2, [sigma3, centroid3, height3] */
        /* paramValueY[0..8] = sigma1, centroid1, height1, sigma2, centroid2, height2, [sigma3, centroid3, height3] */
        for (int ix=0; ix<nFunctions; ix++) {
          for (int iy=0; iy<nFunctions; iy++) {
            tempCentroid[ix+iy*nFunctions][0] = paramValueX[1+ix*3];
            tempCentroid[ix+iy*nFunctions][1] = paramValueY[1+iy*3];
            tempSigma[ix+iy*nFunctions][0] = paramValueX[0+ix*3];
            tempSigma[ix+iy*nFunctions][1] = paramValueY[0*iy*3];

            /* We need 2*Pi*SigmaX*Sigmay here because in biGaussianFunction() the factor 1/(sqrt(2*pi)*sigma) is
             * hidden in the height parameter. We need to remove it for use in the B-E formula.
             * Dividing by qIon is needed because we are making a product of two functions, rho(x) and rho(y),
             * each of which will integrate to qIon.
             */
            tempQ[ix+iy*nFunctions] = 
              paramValueX[2+ix*3] / normX * paramValueY[2+iy*3] / normY * 
              2 * PI * paramValueX[0+ix*3] * paramValueY[0+iy*3] / qIon;
          }
        }
      } else if (ionEffects->ionFieldMethod==ION_FIELD_BILORENTZIAN || ionEffects->ionFieldMethod==ION_FIELD_TRILORENTZIAN) {
        /* paramValueX[0..8] = ax1, centroidx1, heightx1, ax2, centroidx2, heightx2, [ax3, centroidx3, heightx3] */
        /* paramValueY[0..8] = ay1, centroidy1, heighty1, ay2, centroidy2, heighty2, [ay3, centroidy3, heighty3] */
        for (int ix=0; ix<nFunctions; ix++) {
          for (int iy=0; iy<nFunctions; iy++) {
            tempCentroid[ix+iy*nFunctions][0] = paramValueX[1+ix*3];
            tempCentroid[ix+iy*nFunctions][1] = paramValueY[1+iy*3];
            tempSigma[ix+iy*nFunctions][0] = paramValueX[0+ix*3];
            tempSigma[ix+iy*nFunctions][1] = paramValueY[0*iy*3];

            /* Here we account for the bin sizes (normX and normY) and convert the height parameters to those
             * used in a standard Lorentzian, PI*L(0)*a. We also divide out the total
             * charge since otherwise the 2d integral will be qIon^2, instead of qIon. 
             */
            tempQ[ix+iy*nFunctions] = 
              paramValueX[2+ix*3]*PI*paramValueX[0+ix*3]/normX *
              paramValueY[2+iy*3]*PI*paramValueY[0+iy*3]/normY / qIon;
          }
        }
      } else
        bombElegant("invalid field method used for ION_EFFECTS, seek professional help", NULL);
    }

    if (isSlave || !notSinglePart) {
      /*** Determine and apply kicks to beam from the total ion field */
#if MPI_DEBUG
      printf("Applying kicks to electron beam\n");
#endif
      if (qIon && ionSigma[0]>0 && ionSigma[1]>0 && mTotTotal>10) {
        for (ip=0; ip<np; ip++) {
          double kick[2] = {0,0};

          switch (ionEffects->ionFieldMethod) {
          case ION_FIELD_GAUSSIAN:
	    gaussianBeamKick(part[ip], ionCentroid, ionSigma, kick, qIon, me_mks, 1);
            part[ip][1] += kick[0] / c_mks / Po;
            part[ip][3] += kick[1] / c_mks / Po; 
            break;
          case ION_FIELD_BIGAUSSIAN:
          case ION_FIELD_TRIGAUSSIAN:
            kick[0] = kick[1] = 0;
            for (int i=0; i<nFunctions*nFunctions; i++)  {
              gaussianBeamKick(part[ip], tempCentroid[i], tempSigma[i], tempkick, tempQ[i], me_mks, 1);
              kick[0] += tempkick[0];
              kick[1] += tempkick[1];
            }
            part[ip][1] += kick[0] / c_mks / Po;
            part[ip][3] += kick[1] / c_mks / Po; 
            break;
          case ION_FIELD_BILORENTZIAN:
          case ION_FIELD_TRILORENTZIAN:
            kick[0] = kick[1] = 0;
            for (int i=0; i<nFunctions*nFunctions; i++) {
              evaluateVoltageFromLorentzian(tempkick, 
                                           tempSigma[i][0], tempSigma[i][1], 
                                           part[ip][0] - tempCentroid[i][0], part[ip][2] - tempCentroid[i][1]);
              kick[0] += tempQ[i]*tempkick[0];
              kick[1] += tempQ[i]*tempkick[1];
            }
            part[ip][1] -= kick[0] / (Po*particleMassMV*1e6*particleRelSign);
            part[ip][3] -= kick[1] / (Po*particleMassMV*1e6*particleRelSign);
            break;
          default:
            bombElegant("invalid field method used for ION_EFFECTS, seek professional help", NULL);
            break;
	  }
        }
      }

      if (nBunches!=1) {
        /*** Copy bunch coordinates back to original array */
        for (ip=0; ip<np; ip++)
          memcpy(part0[ipBunch[iBunch][ip]], part[ip], sizeof(double)*7);
      }
    }

    if ((ionFieldMethod==ION_FIELD_BIGAUSSIAN || ionFieldMethod==ION_FIELD_BILORENTZIAN ||
         ionFieldMethod==ION_FIELD_TRIGAUSSIAN || ionFieldMethod==ION_FIELD_TRILORENTZIAN) &&
	(SDDS_ionHistogramOutput) && (iPass%ionHistogramOutputInterval == 0)
	&& (ionEffects->sLocation >= ionHistogramOutput_sStart) 
	&& (ionEffects->sLocation <= ionHistogramOutput_sEnd)) {
      long iPlane;
      /* output ion density histogram */
#if USE_MPI
      if (myid==0) {
#endif
	for (iPlane=0; iPlane<2; iPlane++) {
	  long binOffset, activeBins;
	  determineOffsetAndActiveBins(ionEffects->ionHistogram[iPlane], ionEffects->ionBins[iPlane],
				       &binOffset, &activeBins);
	  if (!SDDS_StartPage(SDDS_ionHistogramOutput, activeBins) ||
	      !SDDS_SetParameters(SDDS_ionHistogramOutput, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,
				  "Pass", iPass, "Bunch", iBunch, "t", tNow, "s", ionEffects->sLocation,
				  "qIonsOutside", ionEffects->ionHistogramMissed[iPlane],
				  "fractionIonChargeOutside", 
				  ionEffects->qTotal ? ionEffects->ionHistogramMissed[iPlane]/ionEffects->qTotal : -1,
				  "binSize", ionEffects->ionDelta[iPlane],
				  "binRange", ionEffects->ionRange[iPlane],
				  "nBins", ionEffects->ionBins[iPlane],
				  "Plane", iPlane==0?"x":"y", 
				  NULL) ||
	      !SDDS_SetColumn(SDDS_ionHistogramOutput, SDDS_SET_BY_NAME,
			      ionEffects->xyIonHistogram[iPlane]+binOffset, activeBins, "Position") ||
	      !SDDS_SetColumn(SDDS_ionHistogramOutput, SDDS_SET_BY_NAME,
			      ionEffects->ionHistogram[iPlane]+binOffset, activeBins, "Charge")) {
	    SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors);
	    SDDS_Bomb((char*)"Problem writing ion histogram data");
	  }
	  if (!SDDS_SetColumn(SDDS_ionHistogramOutput, SDDS_SET_BY_NAME,
			      ionEffects->ionHistogramFit[iPlane]+binOffset, activeBins, "ChargeFit") ||
	      !SDDS_SetParameters(SDDS_ionHistogramOutput, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,
				  "fitResidual", ionEffects->xyFitResidual[iPlane],
#if USE_MPI
				  "nEvaluationsBest", ionEffects->nEvaluationsBest[iPlane],
				  "nEvaluationsMin", ionEffects->nEvaluationsMin[iPlane],
				  "nEvaluationsMax", ionEffects->nEvaluationsMax[iPlane],
#else
				  "nEvaluations", ionEffects->nEvaluations[iPlane],
#endif
				  "sigma1", ionEffects->xyFitParameter2[iPlane][0],
				  "sigma2", ionEffects->xyFitParameter2[iPlane][3],
				  "centroid1", ionEffects->xyFitParameter2[iPlane][1],
				  "centroid2", ionEffects->xyFitParameter2[iPlane][4],
				  "q1", ionEffects->xyFitParameter2[iPlane][2],
				  "q2", ionEffects->xyFitParameter2[iPlane][5],
				  NULL)) {
	    SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors);
	    SDDS_Bomb((char*)"Problem writing ion histogram data");
	  }
          if ((ionFieldMethod==ION_FIELD_TRIGAUSSIAN || ionFieldMethod==ION_FIELD_TRILORENTZIAN) &&
	      !SDDS_SetParameters(SDDS_ionHistogramOutput, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,
				  "sigma3", ionEffects->xyFitParameter3[iPlane][6],
				  "centroid3", ionEffects->xyFitParameter3[iPlane][7],
				  "q3", ionEffects->xyFitParameter3[iPlane][8],
				  NULL)) {
	    SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors);
	    SDDS_Bomb((char*)"Problem writing ion histogram data");
	  }
	  if (!SDDS_WritePage(SDDS_ionHistogramOutput)) {
	    SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors);
	    SDDS_Bomb((char*)"Problem writing ion histogram data");
	  }
	}
#if USE_MPI
      }
#endif	
    }

#if USE_MPI
#ifdef DEBUG
      printf("Preparing to wait on barrier at end of loop for bucket %ld\n", iBunch);
      fflush(stdout);
#endif
      MPI_Barrier(MPI_COMM_WORLD);
#endif
#ifdef DEBUG
      printf("Done with bunch %ld\n", iBunch);
      fflush(stdout);
#endif


#if DEBUG
      if (verbosity > 5) {
	// add up total charge at this element
	double qTotal = 0;
	int jMacro = 0;
	double qPart[ionProperties.nSpecies];

	for (iSpecies=0; iSpecies<ionProperties.nSpecies; iSpecies++) {
	  qPart[iSpecies] = 0;
	  for (jMacro=0; jMacro < ionEffects->nIons[iSpecies]; jMacro++) {
	    qPart[iSpecies] += ionEffects->coordinate[iSpecies][jMacro][4];
	    qTotal += ionEffects->coordinate[iSpecies][jMacro][4];
	  }
	}
	
	FILE * fion;
	fion = fopen("ion_densities.dat", "a");
	// time, element s, charge
	fprintf(fion, "%e  %f  %e %e %e %e %e \n", tNow, ionEffects->sStart, qTotal, qPart[0], qPart[1], qPart[2], qPart[3]);
	fclose(fion);
  }
#endif


#if DEBUG
      // write out coordinates of each ion
      if ((verbosity > 20) && (iBunch == 323)) {
	double xtemp, ytemp, qtemp; 
	int jMacro = 0;
	FILE * fion;
	fion = fopen("ion_coord_all.dat", "a");
	for (iSpecies=0; iSpecies<ionProperties.nSpecies; iSpecies++) {
	  for (jMacro=0; jMacro < ionEffects->nIons[iSpecies]; jMacro++) {
	    xtemp = ionEffects->coordinate[iSpecies][jMacro][0];
	    ytemp = ionEffects->coordinate[iSpecies][jMacro][2];
	    qtemp = ionEffects->coordinate[iSpecies][jMacro][4];
	    fprintf(fion, "%f  %f  %f  %e  %d \n",  ionEffects->sStart, xtemp, ytemp, qtemp, iSpecies);
	  }
	}
	fclose(fion);
      }
#endif


    } /* End of loop over bunches */


#if USE_MPI
    if (myid==0) {
#endif
      iIonEffectsElement++;
      if (iIonEffectsElement==nIonEffectsElements) {
	if (verbosity>10 && SDDS_beamOutput) {
	  printf("Pass %ld, updating output for ion-related electron beam output file following last of %ld IONEFFECTS elements\n", 
		 iPass, nIonEffectsElements);
	  fflush(stdout);
	}
        if (SDDS_beamOutput && !SDDS_UpdatePage(SDDS_beamOutput, FLUSH_TABLE)) {
          SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors);
          SDDS_Bomb((char*)"problem flushing data for ion_effects beam parameters output file");
        }
	if (verbosity>10 && SDDS_ionDensityOutput) {
	  printf("Pass %ld, updating output for ion density output file following last of %ld IONEFFECTS elements\n", 
		 iPass, nIonEffectsElements);
	  fflush(stdout);
	}
        if (SDDS_ionDensityOutput && !SDDS_UpdatePage(SDDS_ionDensityOutput, FLUSH_TABLE)) {
          SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors);
          SDDS_Bomb((char*)"problem flushing data for ion_effects ion density output file");
        }
        iIonEffectsElement = 0;
      }
#if USE_MPI
    }
#endif

#if DEBUG
    // add up total charge at this element
    double qTotal = 0;
    long mTotal = 0;
    int jMacro = 0;
    for (iSpecies=0; iSpecies<ionProperties.nSpecies; iSpecies++) {
      mTotal += ionEffects->nIons[iSpecies];
      for (jMacro=0; jMacro < ionEffects->nIons[iSpecies]; jMacro++) {
	qTotal += ionEffects->coordinate[iSpecies][jMacro][4];
      }
    }

    FILE * fion;
    fion = fopen("ion_info.dat", "a");
    // turn, element s, total macros, total charge
    fprintf(fion, "%d  %f  %d  %e \n", iPass, ionEffects->sStart, mTotal, qTotal);
    fclose(fion);
#endif

#if DEBUG
    // write out coordinates of each ion
    //if (((iPass == 0) || (iPass == 99)) || (verbosity > 9)) {
    if ((iPass % 100 == 0) || (verbosity > 9)) {
      double xtemp, ytemp, qtemp; 
      fion = fopen("ion_coord.dat", "a");
      for (iSpecies=0; iSpecies<ionProperties.nSpecies; iSpecies++) {
	for (jMacro=0; jMacro < ionEffects->nIons[iSpecies]; jMacro++) {
	  xtemp = ionEffects->coordinate[iSpecies][jMacro][0];
	  ytemp = ionEffects->coordinate[iSpecies][jMacro][2];
	  qtemp = ionEffects->coordinate[iSpecies][jMacro][4];
	  fprintf(fion, "%f  %e  %e  %e  %d \n",  ionEffects->sStart, xtemp, ytemp, qtemp, iSpecies);
	}
      }
      fclose(fion);
    }
#endif

  if (time && time!=time0)
    free(time);
  if (isSlave || !notSinglePart)
    free_bucket_assignment_memory(time0, ibParticle, ipBunch, npBunch, nBunches);
  if (speciesCentroid)
    free_zarray_2d((void**)speciesCentroid, ionProperties.nSpecies, 2);
  if (speciesSigma)
    free_zarray_2d((void**)speciesSigma, ionProperties.nSpecies, 2);
  if (speciesCharge)
    free(speciesCharge);
  if (speciesCount)
    free(speciesCount);
}

void addIons(IONEFFECTS *ionEffects, long iSpecies, long nToAdd, double qToAdd,  double centroid[2], double sigma[2], long symmetrize)
{
  long iNew;
  
  /* Allocate space for ion coordinates */
  if (ionEffects->coordinate[iSpecies]==NULL)
    ionEffects->coordinate[iSpecies] 
      = (double**)czarray_2d(sizeof(**(ionEffects->coordinate[iSpecies])), nToAdd, COORDINATES_PER_ION);
  else
    ionEffects->coordinate[iSpecies]
      = (double**)resize_czarray_2d((void**)ionEffects->coordinate[iSpecies], 
                                    sizeof(**(ionEffects->coordinate[iSpecies])), 
                                    ionEffects->nIons[iSpecies]+nToAdd, COORDINATES_PER_ION);

  iNew = ionEffects->nIons[iSpecies];
  ionEffects->nIons[iSpecies] += nToAdd;
  for ( ; iNew<ionEffects->nIons[iSpecies]; iNew++) {
    ionEffects->coordinate[iSpecies][iNew][0] = gauss_rn_lim(centroid[0], sigma[0], 3, random_4) ; /* initial x position */
    ionEffects->coordinate[iSpecies][iNew][1] = 0 ; /* initial x velocity */
    ionEffects->coordinate[iSpecies][iNew][2] =  gauss_rn_lim(centroid[1], sigma[1], 3, random_4) ; /* initial y position */
    ionEffects->coordinate[iSpecies][iNew][3] = 0 ; /* initial y velocity */
    ionEffects->coordinate[iSpecies][iNew][4] = qToAdd ; /* macroparticle charge */
    if (symmetrize) {
      iNew ++;
      ionEffects->coordinate[iSpecies][iNew][0] = centroid[0] - (ionEffects->coordinate[iSpecies][iNew-1][0]-centroid[0]);
      ionEffects->coordinate[iSpecies][iNew][1] = 0 ;
      ionEffects->coordinate[iSpecies][iNew][2] = centroid[1] - (ionEffects->coordinate[iSpecies][iNew-1][2]-centroid[1]);
      ionEffects->coordinate[iSpecies][iNew][3] = 0 ;
      ionEffects->coordinate[iSpecies][iNew][4] = qToAdd ;
    }

    //ionEffects->qIon[iSpecies][iNew] = qToAdd;
  }
  
}

void makeIonHistograms(IONEFFECTS *ionEffects, long nSpecies, double *beamCentroid, double *beamSigma, 
		       double *ionCentroid, double *ionSigma)
{
  long iSpecies, iBin, iPlane;
  long iIon;
  double qTotal;

  for (iPlane=0; iPlane<2; iPlane++) {
    if (beamSigma[iPlane]>0 && ionEffects->binDivisor[iPlane]>1)
      ionEffects->ionDelta[iPlane] = beamSigma[iPlane]/ionEffects->binDivisor[iPlane];
    else
      ionEffects->ionDelta[iPlane] = 1e-3;
    if (ionEffects->rangeMultiplier[iPlane]<0)
      ionEffects->ionRange[iPlane] = -ionSigma[iPlane]*ionEffects->rangeMultiplier[iPlane];
    else
      ionEffects->ionRange[iPlane] = findIonBinningRange(ionEffects, iPlane, nSpecies);
    ionEffects->ionBins[iPlane] = ionEffects->ionRange[iPlane]/ionEffects->ionDelta[iPlane]+0.5;
    if (ionEffects->ionBins[iPlane]>ionHistogramMaxBins) {
      ionEffects->ionBins[iPlane] = ionHistogramMaxBins;
      ionEffects->ionDelta[iPlane] = 2*ionEffects->ionRange[iPlane]/(ionHistogramMaxBins-1.0);
    }

    /* allocate and set x or y coordinates of histogram (for output to file) */
    if (ionEffects->xyIonHistogram[iPlane])
      free(ionEffects->xyIonHistogram[iPlane]);
    ionEffects->xyIonHistogram[iPlane] = (double*)tmalloc(sizeof(*ionEffects->xyIonHistogram[iPlane])*ionEffects->ionBins[iPlane]);
    for (int ixy=0; ixy<ionEffects->ionBins[iPlane]; ixy++)
      ionEffects->xyIonHistogram[iPlane][ixy] = ixy*ionEffects->ionDelta[iPlane] - ionEffects->ionRange[iPlane]/2;

    if (ionEffects->ionHistogram[iPlane])
      free(ionEffects->ionHistogram[iPlane]);
    ionEffects->ionHistogram[iPlane] = (double*)tmalloc(sizeof(*ionEffects->ionHistogram[iPlane])*ionEffects->ionBins[iPlane]);
    memset(ionEffects->ionHistogram[iPlane], 0, sizeof(*ionEffects->ionHistogram[iPlane])*ionEffects->ionBins[iPlane]);

    if (ionEffects->ionHistogramFit[iPlane])
      free(ionEffects->ionHistogramFit[iPlane]);
    ionEffects->ionHistogramFit[iPlane] = (double*)tmalloc(sizeof(*ionEffects->ionHistogramFit[iPlane])*ionEffects->ionBins[iPlane]);

    ionEffects->ionHistogramMissed[iPlane] = 0;

    double delta, xyStart;
    delta = ionEffects->ionDelta[iPlane];
    xyStart = -ionEffects->ionRange[iPlane]/2;

    /* histogram ion charge */
    for (iSpecies=qTotal=0; iSpecies<nSpecies; iSpecies++) {
      for (iIon=0; iIon<ionEffects->nIons[iSpecies]; iIon++) {
	iBin = floor((ionEffects->coordinate[iSpecies][iIon][2*iPlane] - xyStart)/delta);
	if (iBin<0 || iBin>(ionEffects->ionBins[iPlane]-1))
	  ionEffects->ionHistogramMissed[iPlane] += ionEffects->coordinate[iSpecies][iIon][4];
	else
	  ionEffects->ionHistogram[iPlane][iBin] += ionEffects->coordinate[iSpecies][iIon][4];
	qTotal += ionEffects->coordinate[iSpecies][iIon][4];
      }
    }
  }
  ionEffects->qTotal = qTotal;

#if USE_MPI
  shareIonHistograms(ionEffects);
#endif
}

double findIonBinningRange(IONEFFECTS *ionEffects, long iPlane, long nSpecies)
{
  double *histogram;
  double min, max, hrange, delta;
  long i, quickBins = 0, nIons, nIonsMissed;

  max = -(min=DBL_MAX);
  nIons = 0;
  for (long iSpecies=0; iSpecies<nSpecies; iSpecies++) {
    nIons += ionEffects->nIons[iSpecies];
    for (long iIon=0; iIon<ionEffects->nIons[iSpecies]; iIon++) {
      if (ionEffects->coordinate[iSpecies][iIon][2*iPlane]<min)
	min = ionEffects->coordinate[iSpecies][iIon][2*iPlane];
      if (ionEffects->coordinate[iSpecies][iIon][2*iPlane]>max)
	max = ionEffects->coordinate[iSpecies][iIon][2*iPlane];
    }
  }
#if USE_MPI
  double gmin, gmax;
  long nIonsGlobal;
  MPI_Allreduce(&min, &gmin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&max, &gmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&nIons, &nIonsGlobal, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  min = gmin;
  max = gmax;
  nIons = nIonsGlobal;
#endif
  if (abs(max)>abs(min))
    hrange = abs(max) + ionEffects->ionDelta[iPlane];
  else
    hrange = abs(min) + ionEffects->ionDelta[iPlane];
  /*  printf("nIons = %ld, min = %le, max = %le, hrange = %le\n", nIons, min, max, hrange); */
  if (ionEffects->rangeMultiplier[iPlane]==0) {
    /* printf("Using full range [%le, %le] for ion binning\n", -hrange, hrange); fflush(stdout); */
    return 2*hrange;
  }
  delta = 10*ionEffects->ionDelta[iPlane];
  quickBins = (2*hrange)/delta+0.5;
  if (quickBins<50 || nIons<8) {
    /* printf("Using full range [%le, %le] for ion binning\n", -hrange, hrange); fflush(stdout); */
    return 2*hrange;
  }

  /* printf("Using %ld bins, delta = %le for quick ion binning over range [%le, %le]\n",
     quickBins, delta, -hrange, hrange); fflush(stdout); */

  histogram = (double*)calloc(quickBins, sizeof(*histogram));

  /* make charge-weighted histogram */
  nIonsMissed = 0;
  for (long iSpecies=0; iSpecies<nSpecies; iSpecies++) {
    for (long iIon=0; iIon<ionEffects->nIons[iSpecies]; iIon++) {
      long iBin;
      iBin = floor((ionEffects->coordinate[iSpecies][iIon][2*iPlane] + hrange)/delta);
      if (iBin>=0 && iBin<quickBins)
	histogram[iBin] += 1; /* ionEffects->coordinate[iSpecies][iIon][4]; */
      else
	nIonsMissed += 1;
    }
  }

#if USE_MPI
  double *histogramGlobal;
  long nIonsMissedGlobal;
  histogramGlobal = (double*)calloc(quickBins, sizeof(*histogramGlobal));
  MPI_Allreduce(histogram, histogramGlobal, quickBins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&nIonsMissed, &nIonsMissedGlobal, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  /* printf("nIonsMissedGlobal = %ld\n", nIonsMissedGlobal); fflush(stdout); */
  for (i=0; i<quickBins; i++)
    histogram[i] = histogramGlobal[i];
  free(histogramGlobal);
#endif
  
  /* find cumulative distribution */
  /*
  for (i=0; i<quickBins; i++)
    printf("histogram[%ld] = %le\n", i, histogram[i]);
  */
  for (i=1; i<quickBins; i++)
    histogram[i] += histogram[i-1];
  for (i=0; i<quickBins; i++)
    histogram[i] /= histogram[quickBins-1];
  /* 
  for (i=0; i<quickBins; i++)
    printf("cdf[%ld] = %le\n", i, histogram[i]);
  */

  /* find 10% and 90% points */
  min = -hrange;
  max = hrange;
  for (i=0; i<quickBins; i++) {
    if (histogram[i]>0.1) {
      min = -hrange + (i-1)*delta;
      break;
    }
  }
  for (   ; i<quickBins; i++) {
    if (histogram[i]>0.9) {
      max = -hrange + (i-1)*delta;
      break;
    }
  }

  if (abs(max)>abs(min))
    hrange = abs(max);
  else
    hrange = abs(min);
  
  hrange *= abs(ionEffects->rangeMultiplier[iPlane]);
  if (hrange>ionEffects->span[iPlane])
    hrange = ionEffects->span[iPlane];
  /* printf("Using hrange=%le for ion binning\n", hrange); fflush(stdout); */

  free(histogram);
  return 2*hrange;
}

#if USE_MPI
void shareIonHistograms(IONEFFECTS *ionEffects)
{
  double *buffer, partBuffer3[3], sumBuffer3[3];
  long iPlane;

  MPI_Barrier(MPI_COMM_WORLD);

  partBuffer3[0] = ionEffects->qTotal;
  partBuffer3[1] = ionEffects->ionHistogramMissed[0];
  partBuffer3[2] = ionEffects->ionHistogramMissed[1];
  MPI_Allreduce(partBuffer3, sumBuffer3, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  ionEffects->qTotal = sumBuffer3[0];
  ionEffects->ionHistogramMissed[0] = sumBuffer3[1];
  ionEffects->ionHistogramMissed[1] = sumBuffer3[2];

  if (ionEffects->ionBins[0]<=0)
    return;

  for (iPlane=0; iPlane<2; iPlane++) {
    buffer = (double*)calloc(sizeof(*buffer), ionEffects->ionBins[iPlane]);
    MPI_Allreduce(ionEffects->ionHistogram[iPlane], buffer, ionEffects->ionBins[iPlane],
		  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    memcpy(ionEffects->ionHistogram[iPlane], buffer, sizeof(*buffer)*ionEffects->ionBins[iPlane]);
    free(buffer);
  }
}
#endif
  
void determineOffsetAndActiveBins(double *histogram, long nBins, long *binOffset, long *activeBins)
{
  long i, j;
  double min, max, threshold;

  if (nBins<ionHistogramMinOutputBins || nBins<20) {
    *binOffset = 0;
    *activeBins = nBins;
    return;
  }

  find_min_max(&min, &max, histogram, nBins);
  threshold = max/1000;

  for (i=0; i<nBins; i++)
    if (histogram[i]>threshold)
      break;
  if ((*binOffset = i-1)<0)
    *binOffset = 0;

  for (j=nBins-1; j>i; j--)
    if (histogram[j]>threshold)
      break;
  *activeBins = j - *binOffset;
  if (*activeBins<=0) {
    *binOffset = 0;
    *activeBins = nBins;
    return;
  }

  if (*activeBins<ionHistogramMinOutputBins) {
    j = ionHistogramMinOutputBins - *activeBins;
    *binOffset -= j/2;
    *activeBins += j;
    if (*binOffset<0 || (*binOffset+*activeBins)>=nBins) {
      *binOffset = 0;
      *activeBins = ionHistogramMinOutputBins;
      return;
    }
  }
  return;
}


void addIon_point(IONEFFECTS *ionEffects, long iSpecies, double qToAdd,  double x, double y)
{
  long iNew;
  
  /* Allocate space for ion coordinates */
  if (ionEffects->coordinate[iSpecies]==NULL)
    ionEffects->coordinate[iSpecies] 
      = (double**)czarray_2d(sizeof(**(ionEffects->coordinate[iSpecies])), 1, COORDINATES_PER_ION);
  else
    ionEffects->coordinate[iSpecies]
      = (double**)resize_czarray_2d((void**)ionEffects->coordinate[iSpecies], 
                                    sizeof(**(ionEffects->coordinate[iSpecies])), 
                                    ionEffects->nIons[iSpecies]+1, COORDINATES_PER_ION);

  iNew = ionEffects->nIons[iSpecies];
  ionEffects->nIons[iSpecies] += 1;
  ionEffects->coordinate[iSpecies][iNew][0] = x;
  ionEffects->coordinate[iSpecies][iNew][1] = 0 ; /* initial x velocity */
  ionEffects->coordinate[iSpecies][iNew][2] =  y;
  ionEffects->coordinate[iSpecies][iNew][3] = 0 ; /* initial y velocity */
  ionEffects->coordinate[iSpecies][iNew][4] = qToAdd ; /* macroparticle charge */
}


void gaussianBeamKick(double *coord, double center[2], double sigma[2], double kick[2], double charge, 
		      double ionMass, double ionCharge) 
{
  // calculate beam kick on ion, assuming Gaussian beam
  double sx, sy, x, y, sd, Fx, Fy, C1, C2, C3, ay;
  std::complex <double> Fc, w1, w2, Fc0, erf1, erf2;
  long flag, flag2;


  kick[0] = 0;
  kick[1] = 0;
  //  return;
  
  sx = sigma[0];
  sy = sigma[1];
  x = coord[0] - center[0];
  y = coord[2] - center[1];


  C1 = c_mks * charge * re_mks * me_mks * ionCharge / e_mks;

  if (sx > sy) {
    ay = abs(y);
    sd = sqrt(2.0*(sqr(sx)-sqr(sy)));
    w1 = std::complex <double> (x/sd, ay/sd);
    w2 = std::complex <double> (x/sd*sy/sx, ay/sd*sx/sy);

    C2 = sqrt(2*PI / (sqr(sx)-sqr(sy)));
    C3 = exp(-sqr(x)/(2*sqr(sx))-sqr(y)/(2*sqr(sy)));

    erf1 = complexErf(w1, &flag);
    erf2 = complexErf(w2, &flag2);

    Fc = C1 * C2 * (erf1 - C3*erf2);

    Fx = Fc.imag();
    if (y > 0) Fy = Fc.real();
    else Fy = -Fc.real();


  } else {
    sd = sqrt(2.0*(sqr(sy)-sqr(sx)));
    w1 = std::complex <double> (y/sd, abs(x)/sd);
    w2 = std::complex <double> (y/sd*sx/sy, abs(x)/sd*sy/sx);

    Fc0 = std::complex <double> (0, C1 * sqrt(2*PI / (sqr(sy)-sqr(sx))) );

    erf1 = complexErf(w1, &flag);
    erf2 = complexErf(w2, &flag2);

    C3 = exp(-sqr(x)/(2*sqr(sx))-sqr(y)/(2*sqr(sy)));

    Fc = -Fc0 * (erf1 - C3*erf2);

    if (x > 0) Fx = -Fc.imag();
    else Fx = Fc.imag();
    
    Fy = Fc.real();

  }

  kick[0] = -Fx / ionMass;
  kick[1] = -Fy / ionMass;

}



void roundGaussianBeamKick(double *coord, double center[2], double sigma[2], double kick[2], double charge, 
		      double ionMass, double ionCharge) 
{
  // calculate beam kick on ion, assuming round Gaussian beam
  double sx, sy, x, y, sig, r, C1, dp, theta;


  kick[0] = 0;
  kick[1] = 0;
  //  return;
  
  sx = sigma[0];
  sy = sigma[1];
  sig = (sx + sy) / 2;

  x = coord[0] - center[0];
  y = coord[2] - center[1];
  r = sqrt(sqr(x) + sqr(y));

  C1 = 2 * c_mks * charge * re_mks * me_mks * ionCharge / e_mks;

  dp = C1 / r * (1 - exp(-sqr(r)/(2*sqr(sig))));

  theta = atan2(y,x);
  kick[0] = -dp * cos(theta) / ionMass;
  kick[1] = -dp * sin(theta) / ionMass;

}


short multipleWhateverFit(double beamSigma[2], double beamCentroid[2], double *paramValueX, double *paramValueY, 
                   IONEFFECTS *ionEffects, double ionSigma[2], double ionCentroid[2]) {
  double result = 0;
  double paramValue[9], paramDelta[9], lowerLimit[9], upperLimit[9];
  int32_t nEvalMax=distribution_fit_evaluations, nPassMax=distribution_fit_passes;
  unsigned long simplexFlags = SIMPLEX_NO_1D_SCANS;
  double peakVal, minVal, xMin, xMax;
  long fitReturn, dummy, nEvaluations;
#if USE_MPI
  double bestResult, lastBestResult;
  int min_location;
#else
  double lastResult = DBL_MAX;
#endif
  long verbosity = 0, plane;

  for (int i=0; i<3*nFunctions; i++)
    paramValueX[i] = paramValueY[i] = 0;

  for (plane=0; plane<2; plane++) {
    fitReturn = 0;
#if MPI_DEBUG
    printf("Performing fit for %c plane\n", plane?'y':'x');
    fflush(stdout);
#endif

    for (int i=0; i<9; i++)
      paramValue[i] = paramDelta[i] = lowerLimit[i] = upperLimit[i] = 0;

    nData = ionEffects->ionBins[plane];
    for (mFunctions=2; mFunctions<=nFunctions; mFunctions++) {
      xData = ionEffects->xyIonHistogram[plane];
      yData = ionEffects->ionHistogram[plane];
      yFit = ionEffects->ionHistogramFit[plane];
      yDataSum = 0;
      for (int i=0; i<nData; i++)
        yDataSum += yData[i];
      result = find_min_max(&minVal, &peakVal, yData, nData);
      find_min_max(&xMin, &xMax, xData, nData);
      
      /* smaller sigma is close to the beam size, larger is close to ion sigma */
      if (mFunctions==2) {
	if (ionEffects->xyFitSet[plane]&0x01) {
	  paramValue[0] = ionEffects->xyFitParameter2[plane][0];
	  paramValue[1] = ionEffects->xyFitParameter2[plane][1];
	  paramValue[2] = ionEffects->xyFitParameter2[plane][2];
	} else {
	  paramValue[0] = beamSigma[plane];
	  paramValue[1] = beamCentroid[plane];
	  paramValue[2] = peakVal/2;
	}
        paramDelta[0] = paramValue[0]/2;
        paramDelta[1] = abs(beamCentroid[plane])/2;
        paramDelta[2] = peakVal/4;
        lowerLimit[0] = paramValue[0]/100;
        if (ionEffects->sigmaLimitMultiplier[plane]>0 && lowerLimit[0]<(ionEffects->sigmaLimitMultiplier[plane]*ionEffects->ionDelta[plane]))
          lowerLimit[0] = ionEffects->sigmaLimitMultiplier[plane]*ionEffects->ionDelta[plane];
        lowerLimit[1] = xMin/10;
        lowerLimit[2] = peakVal/20;
        upperLimit[0] = paramValue[0]*10;
        if (upperLimit[0]<lowerLimit[0])
          upperLimit[0] = 2*lowerLimit[0];
        upperLimit[1] = xMax/10;
        upperLimit[2] = peakVal;
	
	if (ionEffects->xyFitSet[plane]&0x01) {
	  paramValue[3] = ionEffects->xyFitParameter2[plane][3];
	  paramValue[4] = ionEffects->xyFitParameter2[plane][4];
	  paramValue[5] = ionEffects->xyFitParameter2[plane][5];
	} else {
	  paramValue[3] = ionSigma[plane];
	  paramValue[4] = ionCentroid[plane];
	  paramValue[5] = paramValue[2]/3;
	}
        paramDelta[3] = paramValue[3]/2;
        paramDelta[4] = abs(ionCentroid[plane])/2;
        paramDelta[5] = peakVal/4;
        lowerLimit[3] = paramValue[3]/100;
        if (ionEffects->sigmaLimitMultiplier[plane]>0 && lowerLimit[3]<(ionEffects->sigmaLimitMultiplier[plane]*ionEffects->ionDelta[plane]))
          lowerLimit[3] = ionEffects->sigmaLimitMultiplier[plane]*ionEffects->ionDelta[plane];
        lowerLimit[4] = xMin/10;
        lowerLimit[5] = paramDelta[5]/5;
        upperLimit[3] = paramValue[3]*10;
        if (upperLimit[3]<lowerLimit[3])
          upperLimit[3] = 2*lowerLimit[3];
        upperLimit[4] = xMax/10;
        upperLimit[5] = peakVal;
      } else if (mFunctions==3) {
	if (ionEffects->xyFitSet[plane]&0x02) {
#if USE_MPI
	  /* In this case, the odd processors will use the held-over values from 2-function fit */
	  if (myid%2==0) {
	    memcpy(paramValue, ionEffects->xyFitParameter3[plane], 9*sizeof(double));
	  } else {
	    paramValue[6] = 10*ionSigma[plane];
	    paramValue[7] = 0;
	    paramValue[8] = 0;
	  }
#else
	  memcpy(paramValue, ionEffects->xyFitParameter3[plane], 9*sizeof(double));
#endif
	} else {
	  /* paramValue[0-5] are held over from 2-function fit */
	  paramValue[6] = 10*ionSigma[plane];
	  paramValue[7] = 0;
	  paramValue[8] = 0;
	}
        paramDelta[6] = paramValue[6]/2;
        paramDelta[7] = abs(ionCentroid[plane])/2;
        paramDelta[8] = peakVal/20;
        lowerLimit[6] = paramValue[6]/100;
        if (ionEffects->sigmaLimitMultiplier[plane]>0 && 
            lowerLimit[6]<(ionEffects->sigmaLimitMultiplier[plane]*ionEffects->ionDelta[plane]))
          lowerLimit[6] = ionEffects->sigmaLimitMultiplier[plane]*ionEffects->ionDelta[plane];
        lowerLimit[7] = xMin/10;
        lowerLimit[8] = 0;
        upperLimit[6] = paramValue[6]*10;
        if (upperLimit[6]<lowerLimit[6])
          upperLimit[6] = 2*lowerLimit[6];
        upperLimit[7] = xMax/10;
        upperLimit[8] = peakVal;
      }

    for (int i=0; i<3*mFunctions; i++) {
      if (lowerLimit[i]>paramValue[i]) {
	if (paramValue[i]<0)
	  lowerLimit[i] = 10*paramValue[i];
	else
	  lowerLimit[i] = paramValue[i]/10;
      }
      if (upperLimit[i]<paramValue[i]) {
	if (paramValue[i]<0)
	  upperLimit[i] = paramValue[i]/10;
	else
	  upperLimit[i] = paramValue[i]*10;
      }
    }

#if USE_MPI
    /* Randomize step sizes and starting points */
    for (int i=0; i<3*mFunctions; i++) {
      if (myid!=0 && (mFunctions!=3 || myid!=1)) {
	paramValue[i] *= (1+(random_2(0)-0.5)/5);
	if (paramValue[i]<lowerLimit[i])
	  paramValue[i] = lowerLimit[i];
	if (paramValue[i]>upperLimit[i])
	  paramValue[i] = upperLimit[i];
	paramDelta[i] *= random_2(0)*9.9+0.1;
      }
    }
    lastBestResult = DBL_MAX;
#endif
    int nTries = distribution_fit_restarts;
    nEvaluations = 0;
#if !USE_MPI
    lastResult = DBL_MAX;
#endif
    while (nTries--) {
      fitReturn +=  simplexMin(&result, paramValue, paramDelta, lowerLimit, upperLimit,
			       NULL, mFunctions*3, distribution_fit_target, distribution_fit_tolerance/10.0, 
			       ionEffects->ionFieldMethod==ION_FIELD_BIGAUSSIAN || ionEffects->ionFieldMethod==ION_FIELD_TRIGAUSSIAN
			       ?multiGaussianFunction:multiLorentzianFunction,
			       (verbosity>0?report:NULL) , nEvalMax, nPassMax, 12, 3, 1.0, simplexFlags);
      if (fitReturn>=0)
	nEvaluations += fitReturn;
#if USE_MPI
#if MPI_DEBUG
      printf("Waiting on barrier after simplexMin, return=%ld, result=%le\n", fitReturn, result); fflush(stdout);
#endif
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Allreduce(&result, &bestResult, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#if MPI_DEBUG
      double worstResult;
      MPI_Allreduce(&result, &worstResult, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      printf("nTries = %d, bestResult = %le, worstResult = %le\n", nTries, bestResult, worstResult); fflush(stdout);
#endif
      if (bestResult<distribution_fit_target || (lastBestResult - bestResult)<distribution_fit_tolerance) {
	break;
      }
      lastBestResult = bestResult;
#else
      if (result<distribution_fit_target || (lastResult-result)<distribution_fit_tolerance)
	break;
      lastResult = result;
#endif
    
#if USE_MPI
      if (nTries!=0) {
        MPI_Barrier(MPI_COMM_WORLD);
	min_location = 1;
	findGlobalMinIndex(&result, &min_location, MPI_COMM_WORLD);
#if MPI_DEBUG
        printf("distributing and randomizing bestResult from processor %d\n", min_location);
        fflush(stdout);
#endif
	MPI_Bcast(paramValue, 3*mFunctions, MPI_DOUBLE, min_location, MPI_COMM_WORLD);
	for (int i=0; i<3*mFunctions; i++) {
          paramValue[i] *= (1+(random_2(0)-0.5)/20);
	  if (paramValue[i]<lowerLimit[i])
	    paramValue[i] = lowerLimit[i];
	  if (paramValue[i]>upperLimit[i])
	    paramValue[i] = upperLimit[i];
        }
      }
#endif
    }

#if USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
    min_location = 1;
    findGlobalMinIndex(&result, &min_location, MPI_COMM_WORLD);
#if MPI_DEBUG
    printf("distributing best result from processor %d\n", min_location);
    fflush(stdout);
#endif
    MPI_Bcast(paramValue, 3*mFunctions, MPI_DOUBLE, min_location, MPI_COMM_WORLD);
    MPI_Bcast(&fitReturn, 1, MPI_LONG, min_location, MPI_COMM_WORLD);
#endif
    
    if (fitReturn>0) {
      ionEffects->xyFitSet[plane] |= mFunctions==2 ? 0x01 : 0x02;
      for (int i=0; i<3*mFunctions; i++) {
	if (mFunctions<3)
	  ionEffects->xyFitParameter2[plane][i] = paramValue[i];
	else
	  ionEffects->xyFitParameter3[plane][i] = paramValue[i];
      }
    }
    
    if (result<distribution_fit_target)
      break;
    }

#if USE_MPI
#if MPI_DEBUG
    printf("Waiting on barrier after optimization loop, result=%le, fitReturn=%ld, nEvaluations=%ld\n", 
	   result, fitReturn, nEvaluations); 
    fflush(stdout);
#endif
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allreduce(&nEvaluations, &ionEffects->nEvaluationsMin[plane], 1, MPI_LONG, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&nEvaluations, &ionEffects->nEvaluationsMax[plane], 1, MPI_LONG, MPI_MAX, MPI_COMM_WORLD);
    min_location = -1;
    findGlobalMinIndex(&result, &min_location, MPI_COMM_WORLD);
#if MPI_DEBUG
    int max_location;
    findGlobalMaxIndex(&result, &max_location, MPI_COMM_WORLD);
    printf("min,max_location = %d,%d after findGlobalMinIndex\n", min_location, max_location);
#endif
    MPI_Bcast(&nEvaluations, 1, MPI_LONG, min_location, MPI_COMM_WORLD);
    ionEffects->nEvaluationsBest[plane] = nEvaluations;
#else
    ionEffects->nEvaluations[plane] = nEvaluations;
#endif
    if (ionEffects->ionFieldMethod==ION_FIELD_BIGAUSSIAN || ionEffects->ionFieldMethod==ION_FIELD_TRIGAUSSIAN)
      ionEffects->xyFitResidual[plane] = multiGaussianFunction(paramValue, &dummy);
    else
      ionEffects->xyFitResidual[plane] = multiLorentzianFunction(paramValue, &dummy);
    
#if USE_MPI && MPI_DEBUG
    printf("residual is %le\n", ionEffects->xyFitResidual[plane]);
#endif
    
    if (fitReturn>0) {
      for (int i=0; i<nData; i++)
	ionEffects->ionHistogramFit[plane][i] = yFit[i];
    } else {
      ionEffects->xyFitSet[plane] = 0;
      for (int i=0; i<nData; i++)
	ionEffects->ionHistogramFit[plane][i] = -1;
    }
    
    for (int i=0; i<3*mFunctions; i++) {
      if (plane==0)
	paramValueX[i] = paramValue[i];
      else
	paramValueY[i] = paramValue[i];
    }
  }
  
  return ionEffects->ionFieldMethod;
}


double multiGaussianFunction(double *param, long *invalid) 
{
  double sum = 0, sum2 = 0, tmp = 0, result, max = 0, yFitSum = 0;
  double wSumData = 0, wSumFit = 0;

  *invalid = 0;

  //param[3*j+0] = sig[j]
  //param[3*j+1] = cen[j]
  //param[3*j+2] = h[j]

  for (int i=0; i<nData; i++) {
    double z;
    wSumData += xData[i]*yData[i];
    yFit[i] = 0;
    for (int j=0; j<mFunctions; j++) {
      z = (xData[i]-param[3*j+1])/param[3*j+0];
      if (z<6 && z>-6)
        yFit[i] += param[3*j+2] * exp(-z*z/2);
    }
    tmp = abs(yFit[i]-yData[i]);
    sum += tmp;
    sum2 += sqr(tmp);
    yFitSum += yFit[i];
    wSumFit += xData[i]*yFit[i];
    if (tmp>max)
      max = tmp;
  }
  
  switch (residualType) {
  case ION_FIT_RESIDUAL_RMS_DEV:
    result = sqrt(sum2)/yDataSum;
    break;
  case ION_FIT_RESIDUAL_MAX_ABS_DEV:
    result = max/yDataSum;
    break;
  case ION_FIT_RESIDUAL_MAX_PLUS_RMS_DEV:
    result = sqrt(sum2)/yDataSum + max/yDataSum;
    break;
  case ION_FIT_RESIDUAL_SUM_ABS_PLUS_RMS_DEV:
    result = sqrt(sum2)/yDataSum + sum/yDataSum;
    break;
  case ION_FIT_RESIDUAL_RMS_DEV_PLUS_ABS_DEV_SUM:
    result = sqrt(sum2)/yDataSum + abs(yDataSum-yFitSum)/yDataSum;
    break;
  case ION_FIT_RESIDUAL_SUM_ABS_PLUS_ABS_DEV_SUM:
    result = sum/yDataSum + abs(yDataSum-yFitSum)/yDataSum;
    break;
  case ION_FIT_RESIDUAL_RMS_DEV_PLUS_CENTROID:
    if (yFitSum)
      result = sqrt(sum2)/yDataSum + abs(wSumFit/yFitSum-wSumData/yDataSum);
    else
      result = sqrt(sum2)/yDataSum + abs(wSumData/yDataSum);
    break;
  case ION_FIT_RESIDUAL_SUM_ABS_DEV:
    result = sum/yDataSum;
    break;
  default:
    result = 0;
    bombElegant("Invalid residual code in multiLorentzianFunction---seek professional help", NULL);
    break;
  }
    
  return result;
}

double multiLorentzianFunction(double *param, long *invalid) 
{
  double sum = 0, sum2 = 0, tmp = 0, result, max = 0, yFitSum = 0;
  double wSumData = 0, wSumFit = 0;

  *invalid = 0;

  /* The parameters are different from the standard Lorentzian.
   * Instead of A/(pi*a*(1 + (x/a)^2) we use peak/(1 + (x/a)^2)
   */
  //param[3*j+0] = a[j]
  //param[3*j+1] = cen[j]
  //param[3*j+2] = peak[j]
  
  for (int i=0; i<nData; i++) {
    double z;
    wSumData += xData[i]*yData[i];
    yFit[i] = 0;
    for (int j=0; j<mFunctions; j++) {
      z = (xData[i]-param[3*j+1])/param[3*j+0];
      yFit[i] += param[3*j+2]/(1+sqr(z));
    }
    tmp = abs(yFit[i]-yData[i]);
    sum += tmp;
    sum2 += sqr(tmp);
    yFitSum += yFit[i];
    wSumFit += xData[i]*yFit[i];
    if (tmp>max)
      max = tmp;
  }
  
  switch (residualType) {
  case ION_FIT_RESIDUAL_RMS_DEV:
    result = sqrt(sum2)/yDataSum;
    break;
  case ION_FIT_RESIDUAL_MAX_ABS_DEV:
    result = max/yDataSum;
    break;
  case ION_FIT_RESIDUAL_MAX_PLUS_RMS_DEV:
    result = sqrt(sum2)/yDataSum + max/yDataSum;
    break;
  case ION_FIT_RESIDUAL_SUM_ABS_PLUS_RMS_DEV:
    result = sqrt(sum2)/yDataSum + sum/yDataSum;
    break;
  case ION_FIT_RESIDUAL_RMS_DEV_PLUS_ABS_DEV_SUM:
    result = sqrt(sum2)/yDataSum + abs(yDataSum-yFitSum)/yDataSum;
    break;
  case ION_FIT_RESIDUAL_SUM_ABS_PLUS_ABS_DEV_SUM:
    result = sum/yDataSum + abs(yDataSum-yFitSum)/yDataSum;
    break;
  case ION_FIT_RESIDUAL_RMS_DEV_PLUS_CENTROID:
    if (yFitSum)
      result = sqrt(sum2)/yDataSum + abs(wSumFit/yFitSum-wSumData/yDataSum);
    else
      result = sqrt(sum2)/yDataSum + abs(wSumData/yDataSum);
    break;
  case ION_FIT_RESIDUAL_SUM_ABS_DEV:
    result = sum/yDataSum;
    break;
  default:
    result = 0;
    bombElegant("Invalid residual code in multiLorentzianFunction---seek professional help", NULL);
    break;
  }
    
  return result;
}

void report(double y, double *x, long pass, long nEval, long n_dimen)
{
  long i;

  fprintf(stderr, "pass %ld, after %ld evaluations: result = %.16e\na = ", pass, nEval, y);
  for (i=0; i<n_dimen; i++)
    fprintf(stderr, "%.8e ", x[i]);
  fputc('\n', stderr);
}
