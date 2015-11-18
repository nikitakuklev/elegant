/*************************************************************************\
* Copyright (c) 2013 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2013 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

#include "mdb.h"
#include "track.h"
#include "table.h"
#include "fftpackC.h"

void set_up_lrwake(LRWAKE *wakeData, RUN *run, long pass, long particles, CHARGE *charge, long nBuckets);
void convolveArrays(double *output, long outputs, 
                    double *a1, long n1,
                    double *a2, long n2);

void determine_bucket_assignments(double **part, long np, long idSlotsPerBunch, double P0, double **time, long **ibParticle, long ***ipBucket, long **npBucket, 
                                  long *nBuckets, 
                                  long lastNBuckets /* Supply this only if the calling algorithm insists that the number of buckets not change */
                                  )
{
  long ip, ib;
  long ibMin, ibMax;
#if USE_MPI
  long ibMinGlobal, ibMaxGlobal;
  long nBucketsRange[2];
  MPI_Status mpiStatus;
#endif

#ifdef DEBUG
  printf("determine_bucket_assignments called, np=%ld, idSlotsPerBunch=%ld, lastNBuckets=%ld\n", 
         np, idSlotsPerBunch, lastNBuckets);
  printf("pointer check: part %s, time %s, ibParticle %s, ipBucket %s, npBucket %s, nBuckets %s\n",
         part ? "ok" : "NULL",
         time ? "ok" : "NULL",
         ibParticle ? "ok" : "NULL",
         ipBucket ? "ok" : "NULL",
         npBucket ? "ok" : "NULL",
         nBuckets ? "ok" : "NULL");
  fflush(stdout);
#endif

  if (idSlotsPerBunch<=0) {
    ibMin = 0;
    *nBuckets = 1;
#ifdef DEBUG
    printf("determine_bucket_assignments: only one bunch\n");
    fflush(stdout);
#endif
  } else {
    ibMin = LONG_MAX;
    ibMax = LONG_MIN;
    for (ip=0; ip<np; ip++) {
      ib = (part[ip][6]-1)/idSlotsPerBunch;
      if (ib<ibMin)
        ibMin = ib;
      if (ib>ibMax)
        ibMax = ib;
    }
#if USE_MPI
#ifdef DEBUG
    printf("Sharing bucket min/max data: %ld, %ld\n", ibMin, ibMax);
    fflush(stdout);
#endif
    MPI_Allreduce(&ibMin, &ibMinGlobal, 1, MPI_LONG, MPI_MIN, workers);
    MPI_Allreduce(&ibMax, &ibMaxGlobal, 1, MPI_LONG, MPI_MAX, workers);
    ibMin = ibMinGlobal;
    ibMax = ibMaxGlobal;
#endif
    if (ibMin==LONG_MAX || ibMax==LONG_MIN)
      *nBuckets = 0;
    else
      *nBuckets = (ibMax-ibMin)+1;
#ifdef DEBUG
    printf("nPPB=%ld, ibMin = %ld, ibMax = %ld, nBuckets = %ld\n", idSlotsPerBunch, ibMin, ibMax, *nBuckets);
    fflush(stdout);
#endif
  }
  if (*nBuckets) {
    /* To prevent problems in LRWAKE, need to ensure that number of buckets does not increase */
    if (lastNBuckets>0 && *nBuckets>lastNBuckets) {
#ifdef DEBUG
      printf("Error: lastNBuckets = %ld, *nBuckets = %ld\n", lastNBuckets, *nBuckets);
      fflush(stdout);
#endif
      bombElegant("Error: number of bunches has increased.", NULL);
    }

#ifdef DEBUG
    fprintf(stdout, "Performing bucket assignment, nBuckets=%ld\n", *nBuckets);
    fflush(stdout);
#endif
#if USE_MPI
    if (isSlave || !notSinglePart) {
#ifdef DEBUG
      fprintf(stdout, "...performing bucket assignment\n");
      fflush(stdout);
#endif
#endif
      if (np) {
#ifdef DEBUG
        printf("Doing branch for np!=0 (np=%ld)\n", np);
        fflush(stdout);
#endif
        /* Compute time coordinate of each particle */
        if (!(*time = tmalloc(sizeof(**time)*np)))
          bombElegant("Memory allocation problem in determine_bucket_assignments\n", NULL);
        computeTimeCoordinatesOnly(*time, P0, part, np);
#ifdef DEBUG
        printf("Computed time coordinates\n");
        fflush(stdout);
#endif
        *ibParticle = tmalloc(sizeof(**ibParticle)*np);
        if (ipBucket && npBucket) {
#ifdef DEBUG
          printf("Allocating cross-reference arrays\n");
          fflush(stdout);
#endif
          *ipBucket = (long**) czarray_2d(sizeof(***ipBucket), *nBuckets, np);
          *npBucket = (long*) tmalloc(sizeof(*npBucket)*(*nBuckets));
          for (ib=0; ib<(*nBuckets); ib++)
            (*npBucket)[ib] = 0;
        }
#ifdef DEBUG
        printf("Looping over all particles\n");
        fflush(stdout);
#endif
        for (ip=0; ip<np; ip++) {
          if (*nBuckets==1)
            ib = 0;
          else
            ib = (part[ip][6]-1)/idSlotsPerBunch - ibMin;
          if (ib<0 || ib>=(*nBuckets)) {
            fprintf(stdout, "Error: particle outside bunch: ib=%ld, nBuckets=%ld, particleID=%ld\n", ib, *nBuckets, (long)(part[ip][6]));
            exitElegant(1);
          }
          (*ibParticle)[ip] = ib;
          if (ipBucket && npBucket) {
            (*ipBucket)[ib][(*npBucket)[ib]] = ip;
            (*npBucket)[ib] += 1;
          }
        }
      }
      
#if USE_MPI
    }
#endif

  }

#if USE_MPI
  printf("Waiting on barrier at end of determine_bucket_assignment\n");
  fflush(stdout);
  MPI_Barrier(workers);
#endif

#ifdef DEBUG
  printf("Leaving determine_bucket_assignment\n");
  fflush(stdout);
#endif
}

void track_through_lrwake(double **part, long np, LRWAKE *wakeData, double *P0Input,
			  RUN *run, long i_pass, CHARGE *charge
                        )
{
  double *tBucket = NULL;     /* array for <t> for buckets */
  double *xBucket = NULL;     /* array for Q*<x> for bnuches */
  double *yBucket = NULL;     /* array for Q*<y> for buckets */
  double *QBucket = NULL;     /* array for Q for buckets */
  static long lastNBuckets = -1;
  
  double *VzBucket=NULL, *VxBucket=NULL, *VyBucket=NULL; /* arrays of voltage at each bucket */
  double *time = NULL;        /* array to record arrival time of each particle */
  long ib, ip, nBuckets;
  double factor, P0, rampFactor;
  long *ibParticle = NULL;
#if USE_MPI
  double *buffer = NULL;
#endif

#ifdef DEBUG
#if DEBUG>1
  static FILE *fph = NULL, *fpw = NULL;
  if (fph==NULL) {
    fph = fopen("lrwake.history", "w");
    fprintf(fph, "SDDS1\n");
    fprintf(fph, "&parameter name=Pass type=long &end\n");
    fprintf(fph, "&column name=ib type=long &end\n");
    fprintf(fph, "&column name=tMean units=s type=double &end\n");
    fprintf(fph, "&column name=QSum units=C type=double &end\n");
    fprintf(fph, "&data mode=ascii no_row_counts=1 &end\n");
  }
  if (fpw==NULL) {
    fpw = fopen("lrwake.wake", "w");
    fprintf(fpw, "SDDS1\n");
    fprintf(fpw, "&parameter name=Pass type=long &end\n");
    fprintf(fpw, "&column name=t units=s type=double &end\n");
    fprintf(fpw, "&column name=Vz units=V type=double &end\n");
    fprintf(fpw, "&data mode=ascii no_row_counts=1 &end\n");
  }
#endif
#endif

  if (isSlave || !notSinglePart) {
#ifdef DEBUG
    fprintf(stdout, "Running track_through_lrwake, isSlave=%ld, notSinglePart=%ld\n", isSlave, notSinglePart);
#endif

  P0 = *P0Input;
  determine_bucket_assignments(part, np, charge?charge->idSlotsPerBunch:0, P0, &time, &ibParticle, NULL, NULL, &nBuckets, lastNBuckets);
  lastNBuckets = nBuckets;
    
#ifdef DEBUG
  fputs("determine_bucket_assignment returned\n", stdout);
  printf("np = %ld, nBuckets = %ld\n", np, nBuckets);
#endif

  set_up_lrwake(wakeData, run, i_pass, np, charge, nBuckets);

#ifdef DEBUG
  fputs("set_up_lrwake returned\n", stdout);
#endif

  rampFactor = 0;
  if (i_pass>=(wakeData->rampPasses-1))
    rampFactor = 1;
  else
    rampFactor = (i_pass+1.0)/wakeData->rampPasses;

  tBucket = tmalloc(sizeof(*tBucket)*nBuckets);
  xBucket = tmalloc(sizeof(*xBucket)*nBuckets);
  yBucket = tmalloc(sizeof(*yBucket)*nBuckets);
  QBucket = tmalloc(sizeof(*QBucket)*nBuckets);
  for (ib=0; ib<nBuckets; ib++)
    xBucket[ib] = yBucket[ib] = tBucket[ib] = QBucket[ib] = 0;

  for (ip=0; ip<np; ip++) {
    ib = ibParticle[ip];
    tBucket[ib] += time[ip];
    QBucket[ib] += 1;
    xBucket[ib] += part[ip][0];
    yBucket[ib] += part[ip][2];
  }

#if USE_MPI
    /* sum data across processors */
#ifdef MPI_DEBUG
  printf("Summing data across processors\n");
  fflush(stdout);
#endif

    buffer = tmalloc(sizeof(*buffer)*nBuckets);
    MPI_Allreduce(tBucket, buffer, nBuckets, MPI_DOUBLE, MPI_SUM, workers);
    memcpy(tBucket, buffer, sizeof(*tBucket)*nBuckets);

    MPI_Allreduce(QBucket, buffer, nBuckets, MPI_DOUBLE, MPI_SUM, workers);
    memcpy(QBucket, buffer, sizeof(*QBucket)*nBuckets);

    if (wakeData->xFactor && wakeData->WColumn[1]) {
      MPI_Allreduce(xBucket, buffer, nBuckets, MPI_DOUBLE, MPI_SUM, workers);
      memcpy(xBucket, buffer, sizeof(*xBucket)*nBuckets);
    }

    if (wakeData->yFactor && wakeData->WColumn[2]) {
      MPI_Allreduce(yBucket, buffer, nBuckets, MPI_DOUBLE, MPI_SUM, workers);
      memcpy(yBucket, buffer, sizeof(*yBucket)*nBuckets);
    }

    free(buffer);
    buffer = NULL;
#ifdef MPI_DEBUG
    printf("Computing values within buckets\n");
    fflush(stdout);
#endif
#endif

  /* Compute mean time and space position of particles within buckets, multiplied by charge */
  for (ib=0; ib<nBuckets; ib++) {
    if (QBucket[ib]) {
      /* QBucket[ib] holds particle count at this point */
      tBucket[ib] /= QBucket[ib]; 
      xBucket[ib] /= QBucket[ib];
      yBucket[ib] /= QBucket[ib];
      /* multiply by macro-particle charge to QBucket means what it says */
      QBucket[ib] *= charge->macroParticleCharge; 
    }
  }

#ifdef MPI_DEBUG
    printf("Moving history data back\n");
    fflush(stdout);
#endif

  /* Move the existing history data back */
  memmove(wakeData->tHistory, wakeData->tHistory+nBuckets, sizeof(*(wakeData->tHistory))*(wakeData->nHistory-nBuckets));
  memmove(wakeData->QHistory, wakeData->QHistory+nBuckets, sizeof(*(wakeData->QHistory))*(wakeData->nHistory-nBuckets));
  memmove(wakeData->xHistory, wakeData->xHistory+nBuckets, sizeof(*(wakeData->xHistory))*(wakeData->nHistory-nBuckets));
  memmove(wakeData->yHistory, wakeData->yHistory+nBuckets, sizeof(*(wakeData->yHistory))*(wakeData->nHistory-nBuckets));

#ifdef MPI_DEBUG
    printf("Moving new data into buffer\n");
    fflush(stdout);
#endif
  /* Move the new history data into the buffer */
  memmove(wakeData->tHistory+wakeData->nHistory-nBuckets, tBucket, sizeof(*tBucket)*nBuckets);
  memmove(wakeData->QHistory+wakeData->nHistory-nBuckets, QBucket, sizeof(*QBucket)*nBuckets);
  memmove(wakeData->xHistory+wakeData->nHistory-nBuckets, xBucket, sizeof(*xBucket)*nBuckets);
  memmove(wakeData->yHistory+wakeData->nHistory-nBuckets, yBucket, sizeof(*yBucket)*nBuckets);

#ifdef DEBUG
#if DEBUG>1
  fprintf(fph, "%ld\n", i_pass); 
  for (ib=0; ib<wakeData->nHistory; ib++) 
    fprintf(fph, "%ld %e %e\n", ib, wakeData->tHistory[ib], wakeData->QHistory[ib]);
  fprintf(fph, "\n");
  fflush(fph);
#endif
#endif

#ifdef MPI_DEBUG
  printf("Computing wake function at each new bucket\n");
  fflush(stdout);
#endif
  /* Compute the wake function at each new bucket */
  VxBucket = tmalloc(sizeof(*VxBucket)*nBuckets);
  VyBucket = tmalloc(sizeof(*VyBucket)*nBuckets);
  VzBucket = tmalloc(sizeof(*VzBucket)*nBuckets);
  for (ib=0; ib<nBuckets; ib++) {
    long ibh;
    VxBucket[ib] = VyBucket[ib] = VzBucket[ib ] = 0;
#ifdef DEBUG
    fprintf(stdout, "bin %ld: summing from %ld to %ld\n", ib, 0L, (wakeData->nHistory-nBuckets+ib)-1);
#endif
    for (ibh=0; ibh<wakeData->nHistory-nBuckets+ib; ibh++) {
      if (wakeData->QHistory[ibh]) {
	long it;
	double dt;
	dt = tBucket[ib] - wakeData->tHistory[ibh];
	/* interpolate wake functions */
	it = find_nearby_array_entry(wakeData->W[0], wakeData->wakePoints, dt);
	/*
#ifdef DEBUG
	fprintf(stdout, "ib=%ld, ibh=%ld, dt=%le, it=%ld\n", ib, ibh, dt, it);
#endif
	*/
	if (wakeData->W[1])
	  VxBucket[ib] += wakeData->xFactor*wakeData->QHistory[ibh]*wakeData->xHistory[ibh]*
	    linear_interpolation(wakeData->W[1], wakeData->W[0], wakeData->wakePoints, dt, it);
	if (wakeData->W[2])
	  VyBucket[ib] += wakeData->yFactor*wakeData->QHistory[ibh]*wakeData->yHistory[ibh]*
	    linear_interpolation(wakeData->W[2], wakeData->W[0], wakeData->wakePoints, dt, it);
	if (wakeData->W[3])
	  VzBucket[ib] += wakeData->zFactor*wakeData->QHistory[ibh]*
	    linear_interpolation(wakeData->W[3], wakeData->W[0], wakeData->wakePoints, dt, it);
      }
    }
  }

#ifdef DEBUG
#if DEBUG>1
  fprintf(fpw, "%ld\n", i_pass);
  for (ib=0; ib<nBuckets; ib++)
    fprintf(fpw, "%e %e\n", tBucket[ib], VzBucket[ib]);
  fprintf(fpw, "\n");
  fflush(fpw);
#endif
#endif

#ifdef MPI_DEBUG
  printf("Imparting kicks to particles\n");
  fflush(stdout);
#endif
  for (ip=0; ip<np; ip++) {
    /* Impart kick for each particle */
    double dgam, pz;
    factor = wakeData->factor*rampFactor;
    ib = ibParticle[ip];
    dgam = VzBucket[ib]/(1e6*particleMassMV)*factor;
    add_to_particle_energy(part[ip], time[ip], P0, -dgam);
    pz = P0*(1+part[ip][5])/sqrt(1+sqr(part[ip][1])+sqr(part[ip][3]));
    part[ip][1] += VxBucket[ib]*factor/(1e6*particleMassMV)/pz;
    part[ip][3] += VyBucket[ib]*factor/(1e6*particleMassMV)/pz;
  }

  /* Free memory */
  if (time)
    free(time);
  if (tBucket)
    free(tBucket);
  if (xBucket)
    free(xBucket);
  if (yBucket)
    free(yBucket);
  if (QBucket)
    free(QBucket);
  if (ibParticle)
    free(ibParticle);
  if (VxBucket)
    free(VxBucket);
  if (VyBucket)
    free(VyBucket);
  if (VzBucket)
    free(VzBucket);
  }

#if USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif

}

typedef struct {
  char *filename;
  long points;
  int32_t nColumns;
  char **columnName, **units;
  double **data;
} LRWAKE_DATA;

static LRWAKE_DATA *storedWake = NULL;
static long storedWakes = 0;
static char *expectedUnits[4] = {"s", "V/C/m", "V/C/m", "V/C"};

void set_up_lrwake(LRWAKE *wakeData, RUN *run, long pass, long particles, CHARGE *charge, long nBuckets)
{
  SDDS_DATASET SDDSin;
  double tmin, tmax;
  long iw, icol, ic;
#if SDDS_MPI_IO 
/* All the processes will read the wake file, but not in parallel.
   Zero the Memory when call  SDDS_InitializeInput */
  SDDSin.parallel_io = 0; 
#endif
  if (!charge)
    bombElegant("LRWAKE given but no charge element in beamline.", NULL);
  
  if (wakeData->initialized)
    return;

  wakeData->initialized = 1;

  if (!wakeData->inputFile || !strlen(wakeData->inputFile))
    bombElegant("supply inputFile for LRWAKE element", NULL);

  for (icol=0; icol<4; icol++) 
    if (wakeData->WColumn[icol] && !strlen(wakeData->WColumn[icol]))
      wakeData->WColumn[icol] = NULL;
  if (!wakeData->WColumn[0])
    bombElegant("supply tColumn for LRWAKE element", NULL);
  if (!wakeData->WColumn[1] && !wakeData->WColumn[2] && !wakeData->WColumn[3])
    bombElegant("supply at least one of WxColumn, WyColumn, or WzColumn for LRWAKE element", NULL);
    
  for (icol=0; icol<4; icol++)
    wakeData->W[icol] = NULL;

  for (iw=0; iw<storedWakes; iw++) {
    if (strcmp(storedWake[iw].filename, wakeData->inputFile)==0)
      break;
  }

  if (iw==storedWakes) {
    /* read new wake file */
#ifdef DEBUG
    fputs("Reading new wake file for LRWAKE\n", stdout);
#endif
    if (!(storedWake = SDDS_Realloc(storedWake, sizeof(*storedWake)*(storedWakes+1)))) {
      fprintf(stdout, "Error: memory allocation failur storing wake data from LRWAKE file %s\n", wakeData->inputFile);
      exitElegant(1);
    }
    cp_str(&(storedWake[storedWakes].filename), wakeData->inputFile);
    if (!SDDS_InitializeInputFromSearchPath(&SDDSin, wakeData->inputFile) || SDDS_ReadPage(&SDDSin)!=1) {
      fprintf(stdout, "Error: unable to open or read LRWAKE file %s\n", wakeData->inputFile);
      exitElegant(1);
    }
    if ((storedWake[storedWakes].points=SDDS_RowCount(&SDDSin))<0) {
      fprintf(stdout, "Error: no data in LRWAKE file %s\n",  wakeData->inputFile);
      exitElegant(1);
    }
    if (storedWake[storedWakes].points<2) {
      fprintf(stdout, "Error: too little data in LRWAKE file %s\n",  wakeData->inputFile);
      exitElegant(1);
    }
    if (!(storedWake[storedWakes].columnName = SDDS_GetColumnNames(&SDDSin, &storedWake[storedWakes].nColumns)) ||
	storedWake[storedWakes].nColumns<1) {
      fprintf(stdout, "Error: too few columns in LRWAKE file %s\n",  wakeData->inputFile);
      exitElegant(1);
    }
    storedWake[storedWakes].data = tmalloc(sizeof(*(storedWake[storedWakes].data))*storedWake[storedWakes].nColumns);
    storedWake[storedWakes].units = tmalloc(sizeof(*(storedWake[storedWakes].units))*storedWake[storedWakes].nColumns);
#ifdef DEBUG
    fputs("Reading new wake file for LRWAKE (1)\n", stdout);
#endif
    for (ic=0; ic<storedWake[storedWakes].nColumns; ic++) {
      storedWake[storedWakes].data[ic] = NULL;
      storedWake[storedWakes].units[ic] = NULL;
      if (SDDS_CheckColumn(&SDDSin, storedWake[storedWakes].columnName[ic], NULL, SDDS_ANY_FLOATING_TYPE, stdout)!=SDDS_CHECK_WRONGTYPE &&
	  (!(storedWake[storedWakes].data[ic] = 
	     SDDS_GetColumnInDoubles(&SDDSin, storedWake[storedWakes].columnName[ic])) ||
	   !(SDDS_GetColumnInformation(&SDDSin, "units", &(storedWake[storedWakes].units[ic]), SDDS_BY_NAME, storedWake[storedWakes].columnName[ic])))) {
	fprintf(stdout, "Error: problem reading column %s from LRWAKE file %s\n",  storedWake[storedWakes].columnName[ic], wakeData->inputFile);
	exitElegant(1);
      }
    }
#ifdef DEBUG
    fputs("Reading new wake file for LRWAKE (2)\n", stdout);
#endif

    SDDS_Terminate(&SDDSin);
    storedWakes++;
  } 

  /* search stored data for the required columns */
  for (iw=0; iw<storedWakes; iw++) {
    /* Scan over all stored wake input files */
    if (strcmp(storedWake[iw].filename, wakeData->inputFile)==0) {
#ifdef DEBUG
      fprintf(stdout, "Using wakefile %ld (%s)\n", iw, storedWake[iw].filename);
#endif
      /* Found the right file */
      wakeData->wakePoints = storedWake[iw].points;
      for (icol=0; icol<4; icol++) {
#ifdef DEBUG
	fprintf(stdout, "Looking for column %ld \n", icol);
#endif
	/* check for data for each type of data (t, Wx, Wy, Wz) */
	if (wakeData->WColumn[icol]  && strlen(wakeData->WColumn[icol])) {
#ifdef DEBUG
	  fprintf(stdout, "Looking for column %s \n", wakeData->WColumn[icol]);
#endif
	  /* user has requested this type of data */
	  long ic;
	  for (ic=0; ic<storedWake[iw].nColumns; ic++) {
#ifdef DEBUG
	    fprintf(stdout, "Comparing data column %ld (%s) \n", ic, storedWake[iw].columnName[ic]);
#endif
	    /* scan over all the columns for this file */
	    if (strcmp(storedWake[iw].columnName[ic], wakeData->WColumn[icol])==0) {
	      /* found the matching column, so check units */
#ifdef DEBUG
	      fprintf(stdout, "Found match\n");
#endif
	      if (!storedWake[iw].units[ic] || !strlen(storedWake[iw].units[ic]) || 
		  strcmp(storedWake[iw].units[ic], expectedUnits[icol])!=0) {
#ifdef DEBUG
		fprintf(stdout, "Units don't match\n");
#endif
		fprintf(stdout, "Expected units %s for column %s, but found %s, for LRWAKE file %s\n", 
			expectedUnits[icol], storedWake[iw].columnName[ic], storedWake[iw].units[ic], wakeData->inputFile);
		exitElegant(1);
	      }
	      /* copy the data pointer */
	      wakeData->W[icol] = storedWake[iw].data[ic];
#ifdef DEBUG
	      fprintf(stdout, "Stored %x for slot %ld\n", wakeData->W[icol], icol);
#endif
	      break;
	    }
	  }
	}
      }
      break;
    }
  }
#ifdef DEBUG
  fprintf(stdout, "Completed seerch loop\n");
#endif

  for (icol=0; icol<4; icol++) {
    if (wakeData->WColumn[icol] && strlen(wakeData->WColumn[icol]) && !wakeData->W[icol]) {
      fprintf(stdout, "Data for %s not found in LRWAKE file %s\n", wakeData->WColumn[icol], wakeData->inputFile);
      exitElegant(1);
    }
  }

  find_min_max(&tmin, &tmax, wakeData->W[0], wakeData->wakePoints);
  if (tmin>=tmax) {
    fprintf(stdout, "Error: zero or negative time span in LRWAKE file %s\n",  wakeData->inputFile);
    exitElegant(1);
  }
  if (tmin!=0) {
    fprintf(stdout, "Error: LRWAKE function does not start at t=0 for file %s\n",  wakeData->inputFile);
    exitElegant(1);
  }
  wakeData->dt = (tmax-tmin)/(wakeData->wakePoints-1);

  wakeData->nHistory = wakeData->turnsToKeep*nBuckets;

  wakeData->tHistory = tmalloc(sizeof(*(wakeData->tHistory))*wakeData->nHistory);
  wakeData->xHistory = tmalloc(sizeof(*(wakeData->xHistory))*wakeData->nHistory);
  wakeData->yHistory = tmalloc(sizeof(*(wakeData->yHistory))*wakeData->nHistory);
  wakeData->QHistory = tmalloc(sizeof(*(wakeData->QHistory))*wakeData->nHistory);
  for (iw=0; iw<wakeData->nHistory; iw++) 
    wakeData->tHistory[iw] = wakeData->xHistory[iw] = wakeData->yHistory[iw] = wakeData->QHistory[iw] = 0;
#ifdef DEBUG
  fprintf(stdout, "Returing from lrwake setup\n");
#endif
}

