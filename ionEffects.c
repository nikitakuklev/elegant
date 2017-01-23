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
#include "mdb.h"
#include "track.h"
#include "ionEffects.h"

#define ION_FIELD_GAUSSIAN 0
#define N_ION_FIELD_METHODS 1
static char *ionFieldMethodOption[N_ION_FIELD_METHODS] = {
  "gaussian"
};

static long ionFieldMethod = -1;

static long ionsInitialized = 0;

typedef struct {
  long nGasses;
  char **gasName;
  long nLocations;
  double *s;         /* s[j] is the location of the jth set of pressure samples */
  double **pressure; /* pressure[i][j] is the pressure of the ith species at the jth location */
} PRESSURE_DATA;

static PRESSURE_DATA pressureData;
void readGasPressureData(char *filename);
void computeAverageGasPressures(double sStart, double sEnd, double *pressure);

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

void addIons(IONEFFECTS *ionEffects, long iSpecies, long nToAdd);

void setupIonEffects(NAMELIST_TEXT *nltext, RUN *run)
{
  cp_str(&field_calculation_method, "gaussian");

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
  if (!field_calculation_method || !strlen(field_calculation_method))
    bombElegant("field_calculation_method undefined", NULL);
  if ((ionFieldMethod = match_string(field_calculation_method, ionFieldMethodOption, N_ION_FIELD_METHODS, EXACT_MATCH))<0)
    bombElegantVA("field_calculation_method=\"%s\" not recognized", field_calculation_method);

  readGasPressureData(pressure_profile);

  readIonProperties(ion_properties);
}

void readGasPressureData(char *filename)
{
  /* Assumed file structure:
   * Parameters: Gasses --- SDDS_STRING giving comma- or space-separated list of gas species, e.g., "H2O H2 N2 O2 CO2 CO CH4"
   * Columns:
   * s         --- SDDS_FLOAT or SDDS_DOUBLE giving location in the lattice
   * <gasName> --- SDDS_FLOAT or SDDS_DOUBLE giving pressure of <gasName> in Torr or nT
   */

  SDDS_DATASET SDDSin;
  char *gasColumnList, *ptr;
  long i;
  double dsMin, dsMax, ds, pressureMultiplier;

  if (!SDDS_InitializeInputFromSearchPath(&SDDSin, filename)) {
    printf("Problem opening pressure data file %s\n", filename);
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
  }

  if (!check_sdds_column(&SDDSin, "s", "m"))
    bombElegantVA("Column 's' is missing or does not have units of 'm' in %s", filename);
  if (SDDS_CheckParameter(&SDDSin, "Gasses", NULL, SDDS_STRING, stdout)!=SDDS_CHECK_OK)
    bombElegantVA("Parameters \"Gasses\" is missing or not string type in %s", filename);
  
  if (SDDS_ReadPage(&SDDSin)<=0) 
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);

  if (!SDDS_GetParameter(&SDDSin, "Gasses", &gasColumnList))
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
  
  pressureData.nGasses = 0;
  pressureData.gasName = NULL;
  while ((ptr=get_token(gasColumnList))!=NULL) {
    pressureData.gasName = SDDS_Realloc(pressureData.gasName, sizeof(*(pressureData.gasName))*(pressureData.nGasses+1));
    cp_str(&pressureData.gasName[pressureData.nGasses], ptr);
    pressureData.nGasses += 1;
  }
  printf("%ld gasses listed in %s: ", pressureData.nGasses, filename);
  for (i=0; i<pressureData.nGasses; i++)
    printf("%s%c", pressureData.gasName[i], i==(pressureData.nGasses-1)?'\n':' ');
  free(gasColumnList);

  pressureData.nLocations = SDDS_RowCount(&SDDSin);
  printf("Gas data provided at %ld s locations\n", pressureData.nLocations);
  if (!(pressureData.s=SDDS_GetColumnInDoubles(&SDDSin, "s")))
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);

  pressureData.pressure = (double**)czarray_2d(sizeof(double), pressureData.nLocations, pressureData.nGasses);
  for (i=0; i<pressureData.nGasses; i++) {
    pressureMultiplier = 1;
    if (!check_sdds_column(&SDDSin, pressureData.gasName[i], "Torr") && !check_sdds_column(&SDDSin, pressureData.gasName[i], "T")) {
      pressureMultiplier = 1e-9;
      if (!check_sdds_column(&SDDSin, pressureData.gasName[i], "nT") && !check_sdds_column(&SDDSin, pressureData.gasName[i], "nTorr"))
        bombElegantVA("Column \"%s\" is missing, not floating-point type, or does not have units of \"Torr\" or \"nT\" in %s", 
                      pressureData.gasName[i], filename);
    }
    if (!(pressureData.pressure[i] = SDDS_GetColumnInDoubles(&SDDSin, pressureData.gasName[i]))) {
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
    }
    if (pressureMultiplier!=1) {
      /* Convert to Torr */
      long j;
      for (j=0; j<pressureData.nLocations; j++)
        pressureData.pressure[i][j] *= pressureMultiplier;
    }
  }

  dsMax = -(dsMin = DBL_MAX);
  for (i=1; i<pressureData.nLocations; i++) {
    ds = pressureData.s[i] - pressureData.s[i-1];
    if (ds<=0)
      bombElegantVA("s data is not monotonically increasing in pressure data file %s (%le, %le)", filename, pressureData.s[i-1], pressureData.s[i]);
    if (dsMin>ds)
      dsMin = ds;
    if (dsMax<ds)
      dsMax = ds;
  }
  if (fabs(1-dsMin/dsMax)>1e-3)
      bombElegantVA("s data is not uniformly spaced to within desired 0.1% in pressure data file %s", filename);

  printf("Finished reading pressure data file %s\n", filename);
  fflush(stdout);
}

void readIonProperties(char *filename)
{
  /* Assumed file structure:
   * Columns:
   * IonName       --- SDDS_STRING, Name of the ion, e.g., "H2O+", "CO++"
   * Mass          --- SDDS_FLOAT or SDDS_DOUBLE, in AMU
   * ChargeState   --- SDDS_LONG or SDDS_SHORT, ion charge state (positive integer)
   * SourceName    --- SDDS_STRING, Name of the source molecule for this ion, e.g., "H2O", "CO+"
   * CrossSection  --- SDDS_FLOAT or SDDS_DOUBLE, Cross section for producing ion from source, in m^2 (m$a2$n)
   */

  SDDS_DATASET SDDSin;
  long i;
  char **sourceName;

  if (!SDDS_InitializeInputFromSearchPath(&SDDSin, filename)) {
    printf("Problem opening ion properties data file %s\n", filename);
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
  }
  if (!check_sdds_column(&SDDSin, "Mass", "AMU"))
    bombElegantVA("Column \"Mass\" is missing, not floating-point type, or does not have units of \"AMU\" in %s\n",
                  filename);
  if (!check_sdds_column(&SDDSin, "CrossSection", "m$a2$n") &&
      !check_sdds_column(&SDDSin, "CrossSection", "m^2"))
    bombElegantVA("Column \"CrossSection\" is missing, not floating-point type, or does not have units of \"m$a2$n\" in %s\n",
                  filename);

  if (SDDS_ReadPage(&SDDSin)<=0) 
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
  
  if ((ionProperties.nSpecies = SDDS_RowCount(&SDDSin))<=0)
    bombElegantVA("Ion properties file %s appears to have no rows.\n", filename);

  if (!(ionProperties.ionName = (char**)SDDS_GetColumn(&SDDSin, "IonName")) ||
      !(ionProperties.mass = SDDS_GetColumnInDoubles(&SDDSin, "Mass")) ||
      !(ionProperties.chargeState = SDDS_GetColumnInDoubles(&SDDSin, "ChargeState")) ||
      !(ionProperties.crossSection = SDDS_GetColumnInDoubles(&SDDSin, "CrossSection")) ||
      !(sourceName = SDDS_GetColumn(&SDDSin, "SourceName")))
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);

  if (!(ionProperties.sourceGasIndex = tmalloc(sizeof(*(ionProperties.sourceGasIndex))*ionProperties.nSpecies)) ||
      !(ionProperties.sourceIonIndex = tmalloc(sizeof(*(ionProperties.sourceIonIndex))*ionProperties.nSpecies)))
    bombElegantVA("Memory allocation failure allocating arrays for %ld ion species.\n", ionProperties.nSpecies);

  /* Figure out the source gas or source ion indices */
  for (i=0; i<ionProperties.nSpecies; i++) {
    ionProperties.sourceGasIndex[i] = ionProperties.sourceIonIndex[i] = -1;
    if (ionProperties.chargeState[i]<=0) 
      bombElegantVA("Ion %s has non-positive charge state", ionProperties.ionName[i]);
    if (ionProperties.chargeState[i]==1) {
      if ((ionProperties.sourceGasIndex[i] = match_string(sourceName[i], pressureData.gasName, pressureData.nGasses, EXACT_MATCH))<0) {
        bombElegantVA("Unable to find match to gas source \"%s\" for species \"%s\"", 
                      sourceName[i], ionProperties.ionName[i]);
      }
    } else {
      if ((ionProperties.sourceIonIndex[i] = match_string(sourceName[i], ionProperties.ionName, ionProperties.nSpecies, EXACT_MATCH))<0) 
        bombElegantVA("Unable to find match to ion source \"%s\" for species \"%s\"", 
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
  while (eptr) {
    if (eptr->type == T_CHARGE)
      chargeSeen = 1;
    if (eptr->type == T_IONEFFECTS) {
      if (!chargeSeen)
        bombElegant("ION_EFFECTS element preceeds the CHARGE element", NULL);
      /* Set the start of the s range for this element */
      ionEffects = (IONEFFECTS*)eptr->p_elem;
      ionEffects->sStart = (eptrLast->end_pos + eptr->end_pos)/2;
      /* in case this is the last element in the beamline, set s so the range ends here */
      ionEffects->sEnd = eptr->end_pos;
      if (eptrLast && eptrLast->type == T_IONEFFECTS) {
        /* set the s range for the previous ion effects element */
        ionEffects = (IONEFFECTS*)eptrLast->p_elem;
        ionEffects->sEnd = (eptrLast->end_pos + eptr->end_pos)/2;
      }
      eptrLast = eptr;
    } 
    eptr = eptr->succ;
  }

  eptr = &(beamline->elem);
  while (eptr) {
    if (eptr->type == T_IONEFFECTS) {
      ionEffects = (IONEFFECTS*)eptr->p_elem;
      if (verbosity>1) 
        printf("IONEFFECTS element %s#%ld at s=%le m spans s:[%le, %le] m\n",
               eptr->name, eptr->occurence, eptr->end_pos, ionEffects->sStart, ionEffects->sEnd);

      /* Determine the average pressure for each gas */
      ionEffects->pressure = tmalloc(sizeof(*(ionEffects->pressure))*pressureData.nGasses);
      computeAverageGasPressures(ionEffects->sStart, ionEffects->sEnd, ionEffects->pressure);
      if (verbosity>2) {
        long i;
        printf("Average pressures over s:[%le, %le] m\n", ionEffects->sStart, ionEffects->sEnd);
        for (i=0; i<pressureData.nGasses; i++)
          printf("%s:%.2f nT  ", pressureData.gasName[i], ionEffects->pressure[i]*1e9);
        printf("\n");
        fflush(stdout);
      }

      /* Allocate arrays (should really clean up the arrays first in case this is one in a series of runs) */
      ionEffects->coordinate = calloc(ionProperties.nSpecies, sizeof(*(ionEffects->coordinate)));
      ionEffects->nIons = calloc(ionProperties.nSpecies, sizeof(*(ionEffects->nIons)));
      ionEffects->t = 0;
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
 CHARGE *charge           /* beam charge structure */
 )
{
  long ip, iSpecies, iIon;
  long iBunch, nBunches;
  double *time0 = NULL;          /* array to record arrival time of each particle */
  double **part = NULL;          /* particle buffer for working bunch */
  double *time = NULL;           /* array to record arrival time of each particle in working bunch */
  long *ibParticle = NULL;       /* array to record which bunch each particle is in */
  long **ipBunch = NULL;         /* array to record particle indices in part0 array for all particles in each bunch */
  long *npBunch = NULL;          /* array to record how many particles are in each bunch */
  long np, npTotal, max_np = 0;
  /* properties of the electron beam */
  double centroid[2], sigma[2], tNow, qBunch;

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
    
  if (isSlave || !notSinglePart) {
    /* Determine which bunch each particle is in */
    determine_bucket_assignments(part0, np0, charge?charge->idSlotsPerBunch:0, Po, &time0, &ibParticle, &ipBunch, &npBunch, &nBunches, -1);
#if USE_MPI
    if (mpiAbort)
      return;
#endif

    for (iBunch=0; iBunch<nBunches; iBunch++) {
      if (verbosity>2) {
        printf("Working on bunch %ld\n", iBunch);
        fflush(stdout);
      }
      /* Loop over all bunches */
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

#if !USE_MPI
        if (np==0)
          continue;
#endif

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

      /* Compute Sx, Cx, Sy, Cy, <t> */
      sigma[0] = computeRmsCoordinate(part, 0, np, centroid+0, &npTotal);
      sigma[1] = computeRmsCoordinate(part, 2, np, centroid+1, NULL);
#if USE_MPI
      tNow = computeAverage_p(time, np, workers);
#else
      compute_average(&tNow, time, np);
#endif
      /* total bunch charge */
      qBunch = npTotal*charge->macroParticleCharge;
      if (verbosity>3) {
        printf("np: %ld, <t>: %le, sigma x,y: %le, %le,  centroid x,y: %le, %le,  q: %le\n",
               npTotal, tNow, sigma[0], sigma[1], centroid[0], centroid[1], qBunch);
        fflush(stdout);
      }

      /*** Advance the ion positions */
      for (iSpecies=0; iSpecies<ionProperties.nSpecies; iSpecies++) {
        for (iIon=0; iIon<ionEffects->nIons[iIon]; iIon++) {
          ionEffects->coordinate[iSpecies][iIon][0] += (tNow-ionEffects->t)*ionEffects->coordinate[iSpecies][iIon][1];
          ionEffects->coordinate[iSpecies][iIon][2] += (tNow-ionEffects->t)*ionEffects->coordinate[iSpecies][iIon][3];
        }
      }
      ionEffects->t = tNow;

      if (ionEffects->xSpan || ionEffects->ySpan) {
        /*** Eliminate ions that are outside the simulation region */
        for (iSpecies=0; iSpecies<ionProperties.nSpecies; iSpecies++) {
          for (iIon=0; iIon<ionEffects->nIons[iIon]; iIon++) {
            if ((ionEffects->xSpan && fabs(ionEffects->coordinate[iSpecies][iIon][0])>ionEffects->xSpan) ||
                (ionEffects->ySpan && fabs(ionEffects->coordinate[iSpecies][iIon][2])>ionEffects->ySpan)) {
              if (ionEffects->nIons[iIon]>1) {
                /* Move ion at the top of the array into this slot in the array */
                long k;
                for (k=0; k<4; k++) 
                  ionEffects->coordinate[iSpecies][iIon][k] = ionEffects->coordinate[iSpecies][ionEffects->nIons[iIon]-1][k];
              }
              ionEffects->nIons[iSpecies] -= 1;
              iIon -= 1;
            }
          }
        }
      }
      
      /*** Generate ions */
      for (iSpecies=0; iSpecies<ionProperties.nSpecies; iSpecies++) {
        long nToAdd = 0, index;
        if ((index=ionProperties.sourceGasIndex[iSpecies])>=0) {
          /* this is a singly-ionized molecule, so use source gas 
             nToAdd =  someFunctionOfPressure(ionEffects->pressure[index], ...);
          */
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
        }
        addIons(ionEffects, iSpecies, nToAdd);
      }

      /*** Determine and apply kicks from beam to ions */
      for (iSpecies=0; iSpecies<ionProperties.nSpecies; iSpecies++) {
        double Ex = 0, Ey = 0;
        /* Relevant quantities:
         * ionProperties.chargeState[iSpecies] --- Charge state of the ion (integer)
         * ionProperties.mass[iSpecies] --- Mass of the ion (AMUs)
         * qBunch --- bunch charge (C)
         * sigma[0], sigma[1] --- x, y sigma (m)
         * centroid[0], centroid[1] --- x , y centroid (m)
         * ionEffects->nIons[index] --- Number of ions of the source species
         * ionEffects->coordinate[index][j][k] --- kth coordinate of jth source ion 
         */
        for (iIon=0; iIon<ionEffects->nIons[iIon]; iIon++) {
          /*
            Get effective fields in V/m
            BeamFieldFunction(qBunch, sigma, centroid, &(ionEffects->coordinate[iSpecies][iIon][0]), &Ex, &Ey);
          */
          ionEffects->coordinate[iSpecies][iIon][1] += 
            Ex*ionProperties.chargeState[iSpecies]/ionProperties.mass[iSpecies]; /* plus other constants */
          ionEffects->coordinate[iSpecies][iIon][3] += Ey*ionProperties.chargeState[iSpecies]/ionProperties.mass[iSpecies]; /* plus other constants */
        }
      }

      /*** Compute field due to ions */
      for (iSpecies=0; iSpecies<ionProperties.nSpecies; iSpecies++) {
        /* Relevant quantities:
         * ionProperties.chargeState[iSpecies] --- Charge state of the ion (integer)
         * ionEffects->nIons[index] --- Number of ions of the source species
         * ionEffects->coordinate[index][j][k] --- kth coordinate of jth source ion 
         */
      }

      /*** Determine and apply kicks to beam from the total ion field */
      for (ip=0; ip<np; ip++) {
        /* Here, for illustration, just give a 1 urad x kick and a -0.5 urad y kick */
        part[ip][1] += 1e-6;
        part[ip][3] -= 0.5e-6;
      }
      
      if (nBunches!=1) {
        /*** Copy bunch coordinates back to master array */
        for (ip=0; ip<np; ip++)
          memcpy(part0[ipBunch[iBunch][ip]], part[ip], sizeof(double)*7);
      }

#if USE_MPI
#ifdef DEBUG
      printf("Preparing to wait on barrier at end of loop for bucket %ld\n", iBunch);
      fflush(stdout);
#endif
      MPI_Barrier(workers);
#endif
#ifdef DEBUG
      printf("Done with bunch %ld\n", iBunch);
      fflush(stdout);
#endif
    } /* End of loop over bunches */
  } /* End of branch restricting execution to worker nodes */
}

void computeAverageGasPressures(double sStart, double sEnd, double *pressure)
{
  double sum;
  long iGas, iLocation, iStart, iEnd;

  /* Find the indices spanning the desired region */
  iStart = iEnd = -1;
  for (iLocation=0; iLocation<pressureData.nLocations; iLocation++) {
    if (pressureData.s[iLocation]>=sStart && iStart==-1)
      iStart = iLocation;
    if (pressureData.s[iLocation]<=sEnd)
      iEnd = iLocation;
  }
  if (iStart==-1 || iEnd==-1 || iEnd<=iStart)
    bombElegantVA("Failed to find indices corresponding to pressure region s:[%le, %le] m\n",
                  sStart, sEnd);

  for (iGas=0; iGas<pressureData.nGasses; iGas++) {
    sum = 0;
    for (iLocation=iStart; iLocation<=iEnd; iLocation++) 
      sum += pressureData.pressure[iGas][iLocation];
    pressure[iGas] = sum/(iEnd-iStart+1);
  }
}

void addIons(IONEFFECTS *ionEffects, long iSpecies, long nToAdd)
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
    ionEffects->coordinate[iSpecies][iNew][0] = 0 ; /* initial x position */
    ionEffects->coordinate[iSpecies][iNew][1] = 0 ; /* initial x velocity */
    ionEffects->coordinate[iSpecies][iNew][2] = 0 ; /* initial y position */
    ionEffects->coordinate[iSpecies][iNew][3] = 0 ; /* initial y velocity */
  }
  
}
