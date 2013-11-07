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

void set_up_lrwake(LRWAKE *wakeData, RUN *run, long pass, long particles, CHARGE *charge);
void convolveArrays(double *output, long outputs, 
                    double *a1, long n1,
                    double *a2, long n2);

/* #define DEBUG 1  */

void determine_bucket_assignments(double **part, double P0, double **time, long **ibParticle, long np, long nBuckets, double revolutionLength) 
{
  double tmin, tmax, beta0;
  double Trev, dtBucket, tStart;
  long ip, ib;

  /* Compute time coordinate of each particle */
  *time = tmalloc(sizeof(**time)*np);
  computeTimeCoordinates(*time, P0, part, np);

  find_min_max(&tmin, &tmax, *time, np);

  /* define boundaries of buckets */
  beta0 = P0/sqrt(P0*P0+1);
  Trev = revolutionLength/(beta0*c_mks);
  dtBucket = Trev/nBuckets;
  tStart = tmin-dtBucket/2;

  *ibParticle = tmalloc(sizeof(**ibParticle)*np);
  for (ip=0; ip<np; ip++) {
    ib = ((*time)[ip]-tStart)/dtBucket;
    if (ib<0 || ib>=nBuckets) {
      fprintf(stderr, "Error: particle outside LRWAKE region, ib=%ld, nBuckets=%ld\n", ib, nBuckets);
      fprintf(stderr, "t=%21.15e, tStart=%21.15e, dtBucket=%21.15e\n", (*time)[ip], tStart, dtBucket);
      exitElegant(1);
    }
    (*ibParticle)[ip] = ib;
#ifdef DEBUG
    fprintf(stderr, "Particle %ld, t=%le assigned to bucket %ld\n", ip, (*time)[ip], ib);
#endif
  }
}

void track_through_lrwake(double **part, long np, LRWAKE *wakeData, double *P0Input,
			  RUN *run, long i_pass, double revolutionLength, CHARGE *charge
                        )
{
  double *tBucket = NULL;     /* array for <t> for buckets */
  double *xBucket = NULL;     /* array for Q*<x> for bnuches */
  double *yBucket = NULL;     /* array for Q*<y> for buckets */
  double *QBucket = NULL;     /* array for Q for buckets */

  double *VzBucket, *VxBucket, *VyBucket; /* arrays of voltage at each bucket */
  double *time = NULL;        /* array to record arrival time of each particle */
  long ib, ip;
  double factor, P0, rampFactor;
  long *ibParticle;
#ifdef DEBUG
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

  set_up_lrwake(wakeData, run, i_pass, np, charge);
  P0 = *P0Input;
  determine_bucket_assignments(part, P0, &time, &ibParticle, np, wakeData->nBuckets, revolutionLength);
#ifdef DEBUG
  fputs("set_up_lrwake returned\n", stderr);
#endif
  rampFactor = 0;
  if (i_pass>=(wakeData->rampPasses-1))
    rampFactor = 1;
  else
    rampFactor = (i_pass+1.0)/wakeData->rampPasses;

  tBucket = tmalloc(sizeof(*tBucket)*wakeData->nBuckets);
  xBucket = tmalloc(sizeof(*xBucket)*wakeData->nBuckets);
  yBucket = tmalloc(sizeof(*yBucket)*wakeData->nBuckets);
  QBucket = tmalloc(sizeof(*QBucket)*wakeData->nBuckets);
  for (ib=0; ib<wakeData->nBuckets; ib++)
    xBucket[ib] = yBucket[ib] = tBucket[ib] = QBucket[ib] = 0;

  for (ip=0; ip<np; ip++) {
    ib = ibParticle[ip];
    tBucket[ib] += time[ip];
    QBucket[ib] += 1;
    xBucket[ib] += part[ip][0];
    yBucket[ib] += part[ip][2];
  }
  /* Compute mean time and space position of particles within buckets, multiplied by charge */
  for (ib=0; ib<wakeData->nBuckets; ib++) {
    if (QBucket[ib]) {
      /* QBucket[ib] holds particle count at this point */
      tBucket[ib] /= QBucket[ib]; 
      xBucket[ib] /= QBucket[ib];
      yBucket[ib] /= QBucket[ib];
      /* multiply by macro-particle charge to QBucket means what it says */
      QBucket[ib] *= charge->macroParticleCharge; 
    }
  }

  /* Move the existing history data back */
  memmove(wakeData->tHistory, wakeData->tHistory+wakeData->nBuckets, sizeof(*(wakeData->tHistory))*(wakeData->nHistory-wakeData->nBuckets));
  memmove(wakeData->QHistory, wakeData->QHistory+wakeData->nBuckets, sizeof(*(wakeData->QHistory))*(wakeData->nHistory-wakeData->nBuckets));
  memmove(wakeData->xHistory, wakeData->xHistory+wakeData->nBuckets, sizeof(*(wakeData->xHistory))*(wakeData->nHistory-wakeData->nBuckets));
  memmove(wakeData->yHistory, wakeData->yHistory+wakeData->nBuckets, sizeof(*(wakeData->yHistory))*(wakeData->nHistory-wakeData->nBuckets));

  /* Move the new history data into the buffer */
  memmove(wakeData->tHistory+wakeData->nHistory-wakeData->nBuckets, tBucket, sizeof(*tBucket)*wakeData->nBuckets);
  memmove(wakeData->QHistory+wakeData->nHistory-wakeData->nBuckets, QBucket, sizeof(*QBucket)*wakeData->nBuckets);
  memmove(wakeData->xHistory+wakeData->nHistory-wakeData->nBuckets, xBucket, sizeof(*xBucket)*wakeData->nBuckets);
  memmove(wakeData->yHistory+wakeData->nHistory-wakeData->nBuckets, yBucket, sizeof(*yBucket)*wakeData->nBuckets);

#ifdef DEBUG
  fprintf(fph, "%ld\n", i_pass); 
  for (ib=0; ib<wakeData->nHistory; ib++) 
    fprintf(fph, "%ld %e %e\n", ib, wakeData->tHistory[ib], wakeData->QHistory[ib]);
  fprintf(fph, "\n");
  fflush(fph);
#endif

  /* Compute the wake function at each new bucket */
  VxBucket = tmalloc(sizeof(*VxBucket)*wakeData->nBuckets);
  VyBucket = tmalloc(sizeof(*VyBucket)*wakeData->nBuckets);
  VzBucket = tmalloc(sizeof(*VzBucket)*wakeData->nBuckets);
  for (ib=0; ib<wakeData->nBuckets; ib++) {
    long ibh;
    VxBucket[ib] = VyBucket[ib] = VzBucket[ib ] = 0;
#ifdef DEBUG
    fprintf(stderr, "bin %ld: summing from %ld to %ld\n", ib, 0, (wakeData->nHistory-wakeData->nBuckets+ib)-1);
#endif
    for (ibh=0; ibh<wakeData->nHistory-wakeData->nBuckets+ib; ibh++) {
      if (wakeData->QHistory[ibh]) {
	long it;
	double dt;
	dt = tBucket[ib] - wakeData->tHistory[ibh];
	/* interpolate wake functions */
	it = find_nearby_array_entry(wakeData->W[0], wakeData->wakePoints, dt);
	/*
#ifdef DEBUG
	fprintf(stderr, "ib=%ld, ibh=%ld, dt=%le, it=%ld\n", ib, ibh, dt, it);
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
  fprintf(fpw, "%ld\n", i_pass);
  for (ib=0; ib<wakeData->nBuckets; ib++)
    fprintf(fpw, "%e %e\n", tBucket[ib], VzBucket[ib]);
  fprintf(fpw, "\n");
  fflush(fpw);
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
  free(time);
  free(tBucket);
  free(xBucket);
  free(yBucket);
  free(QBucket);
  free(ibParticle);
  free(VxBucket);
  free(VyBucket);
  free(VzBucket);
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

void set_up_lrwake(LRWAKE *wakeData, RUN *run, long pass, long particles, CHARGE *charge)
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

  if (wakeData->nBuckets<1)
    bombElegant("n_buckets must be >1 for LRWAKE element", NULL);

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
    fputs("Reading new wake file for LRWAKE\n", stderr);
#endif
    if (!(storedWake = SDDS_Realloc(storedWake, sizeof(*storedWake)*(storedWakes+1)))) {
      fprintf(stderr, "Error: memory allocation failur storing wake data from LRWAKE file %s\n", wakeData->inputFile);
      exitElegant(1);
    }
    cp_str(&(storedWake[storedWakes].filename), wakeData->inputFile);
    if (!SDDS_InitializeInputFromSearchPath(&SDDSin, wakeData->inputFile) || SDDS_ReadPage(&SDDSin)!=1) {
      fprintf(stderr, "Error: unable to open or read LRWAKE file %s\n", wakeData->inputFile);
      exitElegant(1);
    }
    if ((storedWake[storedWakes].points=SDDS_RowCount(&SDDSin))<0) {
      fprintf(stderr, "Error: no data in LRWAKE file %s\n",  wakeData->inputFile);
      exitElegant(1);
    }
    if (storedWake[storedWakes].points<2) {
      fprintf(stderr, "Error: too little data in LRWAKE file %s\n",  wakeData->inputFile);
      exitElegant(1);
    }
    if (!(storedWake[storedWakes].columnName = SDDS_GetColumnNames(&SDDSin, &storedWake[storedWakes].nColumns)) ||
	storedWake[storedWakes].nColumns<1) {
      fprintf(stderr, "Error: too few columns in LRWAKE file %s\n",  wakeData->inputFile);
      exitElegant(1);
    }
    storedWake[storedWakes].data = tmalloc(sizeof(*(storedWake[storedWakes].data))*storedWake[storedWakes].nColumns);
    storedWake[storedWakes].units = tmalloc(sizeof(*(storedWake[storedWakes].units))*storedWake[storedWakes].nColumns);
#ifdef DEBUG
    fputs("Reading new wake file for LRWAKE (1)\n", stderr);
#endif
    for (ic=0; ic<storedWake[storedWakes].nColumns; ic++) {
      storedWake[storedWakes].data[ic] = NULL;
      storedWake[storedWakes].units[ic] = NULL;
      if (SDDS_CheckColumn(&SDDSin, storedWake[storedWakes].columnName[ic], NULL, SDDS_ANY_FLOATING_TYPE, stderr)!=SDDS_CHECK_WRONGTYPE &&
	  (!(storedWake[storedWakes].data[ic] = 
	     SDDS_GetColumnInDoubles(&SDDSin, storedWake[storedWakes].columnName[ic])) ||
	   !(SDDS_GetColumnInformation(&SDDSin, "units", &(storedWake[storedWakes].units[ic]), SDDS_BY_NAME, storedWake[storedWakes].columnName[ic])))) {
	fprintf(stderr, "Error: problem reading column %s from LRWAKE file %s\n",  storedWake[storedWakes].columnName[ic], wakeData->inputFile);
	exitElegant(1);
      }
    }
#ifdef DEBUG
    fputs("Reading new wake file for LRWAKE (2)\n", stderr);
#endif

    SDDS_Terminate(&SDDSin);
    storedWakes++;
  } 

  /* search stored data for the required columns */
  for (iw=0; iw<storedWakes; iw++) {
    /* Scan over all stored wake input files */
    if (strcmp(storedWake[iw].filename, wakeData->inputFile)==0) {
#ifdef DEBUG
      fprintf(stderr, "Using wakefile %ld (%s)\n", iw, storedWake[iw].filename);
#endif
      /* Found the right file */
      wakeData->wakePoints = storedWake[iw].points;
      for (icol=0; icol<4; icol++) {
#ifdef DEBUG
	fprintf(stderr, "Looking for column %ld \n", icol);
#endif
	/* check for data for each type of data (t, Wx, Wy, Wz) */
	if (wakeData->WColumn[icol]  && strlen(wakeData->WColumn[icol])) {
#ifdef DEBUG
	  fprintf(stderr, "Looking for column %s \n", wakeData->WColumn[icol]);
#endif
	  /* user has requested this type of data */
	  long ic;
	  for (ic=0; ic<storedWake[iw].nColumns; ic++) {
#ifdef DEBUG
	    fprintf(stderr, "Comparing data column %ld (%s) \n", ic, storedWake[iw].columnName[ic]);
#endif
	    /* scan over all the columns for this file */
	    if (strcmp(storedWake[iw].columnName[ic], wakeData->WColumn[icol])==0) {
	      /* found the matching column, so check units */
#ifdef DEBUG
	      fprintf(stderr, "Found match\n");
#endif
	      if (!storedWake[iw].units[ic] || !strlen(storedWake[iw].units[ic]) || 
		  strcmp(storedWake[iw].units[ic], expectedUnits[icol])!=0) {
#ifdef DEBUG
		fprintf(stderr, "Units don't match\n");
#endif
		fprintf(stderr, "Expected units %s for column %s, but found %s, for LRWAKE file %s\n", 
			expectedUnits[icol], storedWake[iw].columnName[ic], storedWake[iw].units[ic], wakeData->inputFile);
		exitElegant(1);
	      }
	      /* copy the data pointer */
	      wakeData->W[icol] = storedWake[iw].data[ic];
#ifdef DEBUG
	      fprintf(stderr, "Stored %x for slot %ld\n", wakeData->W[icol], icol);
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
  fprintf(stderr, "Completed seerch loop\n");
#endif

  for (icol=0; icol<4; icol++) {
    if (wakeData->WColumn[icol] && strlen(wakeData->WColumn[icol]) && !wakeData->W[icol]) {
      fprintf(stderr, "Data for %s not found in LRWAKE file %s\n", wakeData->WColumn[icol], wakeData->inputFile);
      exitElegant(1);
    }
  }

  find_min_max(&tmin, &tmax, wakeData->W[0], wakeData->wakePoints);
#if USE_MPI
  if (isSlave && notSinglePart)
    find_global_min_max(&tmin, &tmax, wakeData->wakePoints, workers);      
#endif
  if (tmin>=tmax) {
    fprintf(stderr, "Error: zero or negative time span in WAKE file %s\n",  wakeData->inputFile);
    exitElegant(1);
  }
  if (tmin!=0) {
    fprintf(stderr, "Error: WAKE function does not start at t=0 for file %s\n",  wakeData->inputFile);
    exitElegant(1);
  }
  wakeData->dt = (tmax-tmin)/(wakeData->wakePoints-1);

  wakeData->nHistory = wakeData->turnsToKeep*wakeData->nBuckets;

  wakeData->tHistory = tmalloc(sizeof(*(wakeData->tHistory))*wakeData->nHistory);
  wakeData->xHistory = tmalloc(sizeof(*(wakeData->xHistory))*wakeData->nHistory);
  wakeData->yHistory = tmalloc(sizeof(*(wakeData->yHistory))*wakeData->nHistory);
  wakeData->QHistory = tmalloc(sizeof(*(wakeData->QHistory))*wakeData->nHistory);
  for (iw=0; iw<wakeData->nHistory; iw++) 
    wakeData->tHistory[iw] = wakeData->xHistory[iw] = wakeData->yHistory[iw] = wakeData->QHistory[iw] = 0;
}

