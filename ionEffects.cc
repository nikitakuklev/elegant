/*************************************************************************\
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

#define ION_FIELD_GAUSSIAN 0
#define N_ION_FIELD_METHODS 1
static char *ionFieldMethodOption[N_ION_FIELD_METHODS] = {
  (char*)"gaussian"
};

static long ionFieldMethod = -1;

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

void addIons(IONEFFECTS *ionEffects, long iSpecies, long nToAdd, double qToAdd, double centroid[2], double sigma[2]);

void gaussianBeamKick(double *coord, double center[2], double sigma[2], double kick[2], double charge, 
		      double ionMass, double ionCharge);

static SDDS_DATASET *SDDS_beamOutput = NULL;
static SDDS_DATASET *SDDS_ionDensityOutput = NULL;

static double sStartFirst = -1;
static long iIonEffectsElement = -1, nIonEffectsElements = 0, iBeamOutput, iIonDensityOutput;
static IONEFFECTS *firstIonEffects = NULL; /* first in the lattice */
#if USE_MPI
static long leftIonCounter = 0;
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
  
}

void setupIonEffects(NAMELIST_TEXT *nltext, VARY *control, RUN *run)
{
  cp_str(&field_calculation_method, (char*)"gaussian");

  /* process namelist input */
  set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
  set_print_namelist_flags(0);
  if (processNamelist(&ion_effects, nltext)==NAMELIST_ERROR)
    bombElegant(NULL, NULL);
  if (echoNamelists) print_namelist(stdout, &ion_effects);

  /* Basic check of input values */
  if (!pressure_profile || !strlen(pressure_profile))
    bombElegant("pressure_profile undefined", NULL);
  if (!ion_properties || !strlen(ion_properties))
    bombElegant("ion_properties undefined", NULL);
  if (beam_output)
    beam_output = compose_filename(beam_output, run->rootname);
  if (ion_density_output)
    ion_density_output = compose_filename(ion_density_output, run->rootname);
  if (!field_calculation_method || !strlen(field_calculation_method))
    bombElegant("field_calculation_method undefined", NULL);
  if ((ionFieldMethod = match_string(field_calculation_method, ionFieldMethodOption, N_ION_FIELD_METHODS, EXACT_MATCH))<0)
    bombElegantVA((char*)"field_calculation_method=\"%s\" not recognized", field_calculation_method);

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
        bombElegantVA((char*)"Unable to find match to gas source \"%s\" for species \"%s\"", 
                      sourceName[i], ionProperties.ionName[i]);
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
      ionEffects->macroIons = macro_ions;
      ionEffects->xSpan = x_span;
      ionEffects->ySpan = y_span;
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
  /* properties of the ion cloud */
  double ionCentroid[2], ionSigma[2], qIon;
  double **speciesCentroid=NULL, *speciesCharge=NULL, **speciesSigma=NULL;
  long *speciesCount=NULL;
  double unitsFactor; /* converts Torr to 1/m^3 and mBarns to m^2 */
#if USE_MPI
  MPI_Status mpiStatus;
#endif

  if (!ionsInitialized) {
    bombElegant("IONEFFECTS element seen, but ion_effects command was not given to initialize ion modeling.", NULL);
  }
  if (ionEffects->disable)
    return ;

  if (ionEffects->startPass<0)
    ionEffects->startPass = 0;

  if ((ionEffects->startPass>=0 && iPass<ionEffects->startPass) ||
      (ionEffects->endPass>=0 && iPass>ionEffects->endPass) ||
      (ionEffects->passInterval>=1 && (iPass-ionEffects->startPass)%ionEffects->passInterval!=0))
    return;

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
  
#if USE_MPI
  if (myid==0) {
#endif
    if (iPass==0 && ionEffects->sStart==sStartFirst) {
      iIonEffectsElement = 0;
      if (SDDS_beamOutput) {
        if (!SDDS_StartPage(SDDS_beamOutput, nPasses*(beam_output_all_locations?nIonEffectsElements:1)*nBunches)) {
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
          exitElegant(1);
        }
        iBeamOutput = 0;
      }
      if (SDDS_ionDensityOutput) {
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
      
      if (ionEffects->xSpan || ionEffects->ySpan) {
        /*** Eliminate ions that are outside the simulation region */
        for (iSpecies=0; iSpecies<ionProperties.nSpecies; iSpecies++) {
          for (iIon=0; iIon<ionEffects->nIons[iSpecies]; iIon++) {
            if ((ionEffects->xSpan && fabs(ionEffects->coordinate[iSpecies][iIon][0])>ionEffects->xSpan) ||
                (ionEffects->ySpan && fabs(ionEffects->coordinate[iSpecies][iIon][2])>ionEffects->ySpan)) {
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
      
      /*
      for (iSpecies=0; iSpecies<ionProperties.nSpecies; iSpecies++) {
        for (iIon=0; iIon<ionEffects->nIons[iSpecies]; iIon++) {
          if ((ionEffects->xSpan && fabs(ionEffects->coordinate[iSpecies][iIon][0])>ionEffects->xSpan) ||
              (ionEffects->ySpan && fabs(ionEffects->coordinate[iSpecies][iIon][2])>ionEffects->ySpan)) {
            long k;
            k = 0;
          }
        }
      }
      */
      
      /*** Generate ions */
      for (iSpecies=0; iSpecies<ionProperties.nSpecies; iSpecies++) {
        long nToAdd = 0, index;
	double qToAdd = 0;
        if ((index=ionProperties.sourceGasIndex[iSpecies])>=0) {
          /* this is a singly-ionized molecule, so use source gas 
             nToAdd =  someFunctionOfPressure(ionEffects->pressure[index], ...);
          */
          /* Shouldn't there be some statistics here ? -- MB */
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
            qToAdd = unitsFactor * qBunch * ionEffects->pressure[index] *      \
              ionProperties.crossSection[iSpecies] * (ionEffects->sEnd - ionEffects->sStart) / ionEffects->macroIons;
            
            addIons(ionEffects, iSpecies, nToAdd, qToAdd, centroid, sigma);
          }
        } else if ((index=ionProperties.sourceIonIndex[iSpecies])>=0) {
          /* This is a multiply-ionized molecule, so use source ion density.
           * Relevant quantities:
           * ionEffects->nIons[index] --- Number of ions of the source species
           * ionEffects->coordinate[index][j][k] --- kth coordinate of jth source ion 
           * ionProperties.crossSection[iSpecies] --- Cross section for producing new ion from the source ions
           */
          /* 
             nToAdd = someFunctionOfExistingNumberOfIons(...); 
          */
          bombElegant("Multiple ionization not implemented at this time", NULL);
        }
      }
          
      /*** Determine and apply kicks from beam to ions */
      for (iSpecies=0; iSpecies<ionProperties.nSpecies; iSpecies++) {
        double ionMass, ionCharge, *coord, kick[2];
      
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

	  ionEffects->coordinate[iSpecies][iIon][1] += kick[0];
	  ionEffects->coordinate[iSpecies][iIon][3] += kick[1];

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
	  ionSigma[0] += sqr(ionEffects->coordinate[iSpecies][jMacro][0]-ionCentroid[0])*ionEffects->coordinate[iSpecies][jMacro][4];
	  ionSigma[1] += sqr(ionEffects->coordinate[iSpecies][jMacro][2]-ionCentroid[1])*ionEffects->coordinate[iSpecies][jMacro][4];
          if (ion_species_output) {
            speciesSigma[iSpecies][0] += sqr(ionEffects->coordinate[iSpecies][jMacro][0]-speciesCentroid[iSpecies][0])*ionEffects->coordinate[iSpecies][jMacro][4];
            speciesSigma[iSpecies][1] += sqr(ionEffects->coordinate[iSpecies][jMacro][2]-speciesCentroid[iSpecies][1])*ionEffects->coordinate[iSpecies][jMacro][4];
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
    if (SDDS_ionDensityOutput) {
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

    if (isSlave || !notSinglePart) {

      if (ionSigma[0]/ionSigma[1]>0.9 && ionSigma[0]/ionSigma[1] < 1.0) {
#if DEBUG
	printf("ion sigx = sigy (sx=%e, sy=%e): bunch %d, Pass %d, s=%f \n",  ionSigma[0], ionSigma[1], iBunch, iPass, ionEffects->sLocation);
#endif
	ionSigma[0] = 0.9 * ionSigma[1];
      } else if (ionSigma[0]/ionSigma[1]>1.0 && ionSigma[0]/ionSigma[1] < 1.1) {
	ionSigma[0] = 1.1 * ionSigma[1];
	//	printf("ion sigx = sigy");
      }

      /*** Determine and apply kicks to beam from the total ion field */
      if (qIon && ionSigma[0]>0 && ionSigma[1]>0 && mTotTotal>10) {
        for (ip=0; ip<np; ip++) {
          double kick[2];
          gaussianBeamKick(part[ip], ionCentroid, ionSigma, kick, qIon, me_mks, 1);
          part[ip][1] += kick[0] / c_mks / Po;
          part[ip][3] += kick[1] / c_mks / Po; 

	  //if (kick[0] > 1000) {
	  //  printf("kick wrong");
	  //}

        }
      }

      if (nBunches!=1) {
        /*** Copy bunch coordinates back to original array */
        for (ip=0; ip<np; ip++)
          memcpy(part0[ipBunch[iBunch][ip]], part[ip], sizeof(double)*7);
      }
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
      if (verbosity > 20) {
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
        if (SDDS_beamOutput && !SDDS_UpdatePage(SDDS_beamOutput, FLUSH_TABLE)) {
          SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors);
          SDDS_Bomb((char*)"problem flushing data for ion_effects beam parameters output file");
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

  if (time0) 
    free(time0);
  if (time && time!=time0)
    free(time);
  if (part && part!=part0)
    free_czarray_2d((void**)part, max_np, COORDINATES_PER_PARTICLE);
  if (ibParticle) 
    free(ibParticle);
  if (ipBunch)
    free_czarray_2d((void**)ipBunch, nBunches, np0);
  if (npBunch)
    free(npBunch);
  if (speciesCentroid)
    free_zarray_2d((void**)speciesCentroid, ionProperties.nSpecies, 2);
  if (speciesSigma)
    free_zarray_2d((void**)speciesSigma, ionProperties.nSpecies, 2);
  if (speciesCharge)
    free(speciesCharge);
  if (speciesCount)
    free(speciesCount);
}

void addIons(IONEFFECTS *ionEffects, long iSpecies, long nToAdd, double qToAdd,  double centroid[2], double sigma[2])
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
    
    //ionEffects->qIon[iSpecies][iNew] = qToAdd;
  }


  
}


void gaussianBeamKick(double *coord, double center[2], double sigma[2], double kick[2], double charge, 
		      double ionMass, double ionCharge) 
{
  // calculate beam kick on ion, assuming Gaussian beam
  double sx, sy, x, y, sd, Fx, Fy, C1, C2, C3;
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
    sd = sqrt(2.0*(sqr(sx)-sqr(sy)));
    w1 = std::complex <double> (x/sd, y/sd);
    w2 = std::complex <double> (x/sd*sy/sx, y/sd*sx/sy);

    C2 = sqrt(2*PI / (sqr(sx)-sqr(sy)));
    C3 = exp(-sqr(x)/(2*sqr(sx))-sqr(y)/(2*sqr(sy)));

    erf1 = complexErf(w1, &flag);
    erf2 = complexErf(w2, &flag2);

    Fc = C1 * C2 * (erf1 - C3*erf2);

    //    Fc = C1 * sqrt(2*PI / (sqr(sx)-sqr(sy))) *  (complexErf(w1, &flag) - 
    //					  exp(-sqr(x)/(2*sqr(sx))-sqr(y)/(2*sqr(sy))) * complexErf(w2, &flag2));

  } else {
    sd = sqrt(2.0*(sqr(sy)-sqr(sx)));
    w1 = std::complex <double> (y/sd, -x/sd);
    w2 = std::complex <double> (y/sd*sx/sy, -x/sd*sy/sx);

    Fc0 = std::complex <double> (0, C1 * sqrt(2*PI / (sqr(sy)-sqr(sx))) );

    erf1 = complexErf(w1, &flag);
    erf2 = complexErf(w2, &flag2);

    C3 = exp(-sqr(x)/(2*sqr(sx))-sqr(y)/(2*sqr(sy)));

    Fc = -Fc0 * (erf1 - C3*erf2);

    //Fc = Fc0 * (complexErf(w1, &flag) -  exp(-sqr(x)/(2*sqr(sx))-sqr(y)/(2*sqr(sy))) * complexErf(w2, &flag2));

  }

  Fx = Fc.imag();
  Fy = Fc.real();

  kick[0] = -Fx / ionMass;
  kick[1] = -Fy / ionMass;
  


}
