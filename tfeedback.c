/*************************************************************************\
* Copyright (c) 2003 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2003 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: tfeedback.c
 *
 * Michael Borland, 2003
 */
#include "mdb.h"
#include "track.h"

void transverseFeedbackPickup(TFBPICKUP *tfbp, double **part0, long np0, long pass, double Po, long idSlotsPerBunch)
{
  double sum, position, output;
  long i, j;
  long np;
  double *time0 = NULL;           /* array to record arrival time of each particle */
  long *ibParticle = NULL;        /* array to record which bucket each particle is in */
  long **ipBucket = NULL;                /* array to record particle indices in part0 array for all particles in each bucket */
  long *npBucket = NULL;                 /* array to record how many particles are in each bucket */
  long iBucket, nBuckets;
#if USE_MPI
  long npTotal;
  double sumTotal;
  MPI_Status mpiStatus;
#endif
  
  if (tfbp->initialized==0)
    initializeTransverseFeedbackPickup(tfbp);

  if (isSlave || !notSinglePart) 
    determine_bucket_assignments(part0, np0, tfbp->bunchedBeamMode?idSlotsPerBunch:0, Po, &time0, &ibParticle, &ipBucket, &npBucket, &nBuckets, -1);

#if USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
  if (myid==0)
    MPI_Recv(&nBuckets, 1, MPI_LONG, 1, 1, MPI_COMM_WORLD, &mpiStatus);
  else if (myid==1)
    MPI_Send(&nBuckets, 1, MPI_LONG, 0, 1, MPI_COMM_WORLD);
#endif
#ifdef DEBUG
  printf("TFBPICKUP: %ld bunches\n", nBuckets);
#endif


  if (tfbp->nBunches==0) {
    tfbp->nBunches = nBuckets;
    tfbp->data = SDDS_Realloc(tfbp->data, sizeof(*tfbp->data)*nBuckets);
    tfbp->filterOutput = SDDS_Realloc(tfbp->filterOutput, sizeof(*tfbp->filterOutput)*nBuckets);
    for (i=0; i<nBuckets; i++) 
      if (!(tfbp->data[i] = calloc(TFB_FILTER_LENGTH, sizeof(**tfbp->data))))
        bombElegant("Memory allocation problem for TFBPICKUP", NULL);
  } else if (tfbp->nBunches!=nBuckets)
    bombElegant("Number of bunches has changed while using TFBPICKUP", NULL);
  
  for (iBucket=0; iBucket<nBuckets; iBucket++) {
#if USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    sum = position = np = 0;
    if (isSlave || !notSinglePart) {
      j = tfbp->yPlane?2:0;
      if (nBuckets==1) {
        np = np0;
        for (i=0; i<np0; i++)
          sum += part0[i][j];
      } else {
        np = npBucket[iBucket];
        for (i=0; i<np; i++)
          sum += part0[ipBucket[iBucket][i]][j];
      }
    }
    
#if USE_MPI
    if (myid==0)
      np = 0;
    MPI_Allreduce(&sum, &sumTotal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&np,  &npTotal, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
    if (npTotal>0)
      position = sumTotal/npTotal;
#else
    if (np>0)
      position = sum/np;
#endif

    if (tfbp->rmsNoise) {
      double dposition;
      dposition = gauss_rn_lim(0.0, tfbp->rmsNoise, 2, random_3);
#if USE_MPI
      MPI_Bcast(&dposition, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
      position += dposition;
  }
#ifdef DEBUG
    fprintf(stdout, "TFBPICKUP: Putting value %e in slot %ld for bunch %ld\n",
            position, pass%tfbp->filterLength, iBucket);
#endif
    tfbp->data[iBucket][pass%tfbp->filterLength] = position;
    if (pass<tfbp->filterLength) {
      tfbp->filterOutput[iBucket] = 0;
    } else {
      for (i=output=0, j=pass; i<tfbp->filterLength; i++, j--) {
        output += tfbp->a[i]*tfbp->data[iBucket][j%tfbp->filterLength];
      }
#ifdef DEBUG
      fprintf(stdout, "TFBPICKUP: filter output is %e\n", output);
#endif
      tfbp->filterOutput[iBucket] = output;
    }
  }

  if (time0) 
    free(time0);
  if (ibParticle) 
    free(ibParticle);
  if (ipBucket)
    free_czarray_2d((void**)ipBucket, nBuckets, np0);
  if (npBucket)
    free(npBucket);
}


void initializeTransverseFeedbackPickup(TFBPICKUP *tfbp)
{
  long i;
  double sum;

  if (tfbp->ID==NULL || !strlen(tfbp->ID))
    bombElegant("you must give an ID string for TFBPICKUP", NULL);

  for (i=sum=0; i<TFB_FILTER_LENGTH; i++)
    sum += tfbp->a[i];
  if (fabs(sum)>1e-6)
    fprintf(stdout, "Warning: sum of a[i] is nonzero for TFBPICKUP\n");

  for (i=TFB_FILTER_LENGTH-1; i>=0; i--) {
    if (tfbp->a[i]!=0)
      break;
  }
  if (i<0)
    bombElegant("All filter coefficients are zero for TFBPICKUP", NULL);
  tfbp->filterLength = i+1;

  if (strcmp(tfbp->plane, "x")==0 || strcmp(tfbp->plane, "X")==0) 
    tfbp->yPlane = 0;
  else {
    if (!(strcmp(tfbp->plane, "y")==0 || strcmp(tfbp->plane, "Y")==0))
      bombElegant("PLANE must be x or y for TFBPICKUP", NULL);
    tfbp->yPlane = 1;
  }

  if (tfbp->data && tfbp->nBunches) {
    for (i=0; i<tfbp->nBunches; i++)
      if (tfbp->data)
        free(tfbp->data);
    tfbp->data = NULL;
  }
  tfbp->nBunches = 0;
  
  tfbp->initialized = 1;
}

void transverseFeedbackDriver(TFBDRIVER *tfbd, double **part0, long np0, LINE_LIST *beamline, long pass, long nPasses, char *rootname, double Po, long idSlotsPerBunch)
{
  double kick;
  long i, j;
  double *time0 = NULL;           /* array to record arrival time of each particle */
  long *ibParticle = NULL;        /* array to record which bucket each particle is in */
  long **ipBucket = NULL;                /* array to record particle indices in part0 array for all particles in each bucket */
  long *npBucket = NULL;                 /* array to record how many particles are in each bucket */
  long iBucket, nBuckets;
#if USE_MPI
  MPI_Status mpiStatus;
#endif

  if (isSlave || !notSinglePart) 
    determine_bucket_assignments(part0, np0, tfbd->bunchedBeamMode?idSlotsPerBunch:0, Po, &time0, &ibParticle, &ipBucket, &npBucket, &nBuckets, -1);

#if USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
  if (myid==0)
    MPI_Recv(&nBuckets, 1, MPI_LONG, 1, 1, MPI_COMM_WORLD, &mpiStatus);
  else if (myid==1)
    MPI_Send(&nBuckets, 1, MPI_LONG, 0, 1, MPI_COMM_WORLD);
#endif
#ifdef DEBUG
  printf("TFBDRIVER: %ld bunches\n", nBuckets);
#endif
  
  if (tfbd->initialized==0)
    initializeTransverseFeedbackDriver(tfbd, beamline, nPasses*nBuckets, rootname);

  if (pass==0)
    tfbd->dataWritten = 0;
  
  if (tfbd->nBunches==0) {
    tfbd->nBunches = nBuckets;
    if (!(tfbd->driverSignal = tmalloc(sizeof(*tfbd->driverSignal)*nBuckets)))
      bombElegant("memory allocation failure (transverseFeedbackDriver)", NULL);
    for (iBucket=0; iBucket<nBuckets; iBucket++) 
      tfbd->driverSignal[iBucket] = calloc((tfbd->delay+1+TFB_FILTER_LENGTH), sizeof(**tfbd->driverSignal));
  } else if (tfbd->nBunches!=nBuckets)
    bombElegant("number of bunches has changed while using TFBDRIVER", NULL);

  if (tfbd->nBunches!=tfbd->pickup->nBunches)
    bombElegant("mismatch in number of buckets between TFBDRIVER and TFBPICKUP", NULL);
  
  for (iBucket=0; iBucket<nBuckets; iBucket++) {
#if USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    kick = tfbd->pickup->filterOutput[iBucket]*tfbd->strength;
    if (tfbd->kickLimit>0 && fabs(kick)>tfbd->kickLimit)
      kick = SIGN(kick)*tfbd->kickLimit;

#ifdef DEBUG
    fprintf(stdout, "TFBDRIVER: pass %ld\nstoring kick %e in slot %ld based on filter output of %e\n",
            pass, kick, pass%(tfbd->delay+tfbd->filterLength), tfbd->pickup->filterOutput[iBucket]);
#endif
  
    tfbd->driverSignal[iBucket][pass%(tfbd->delay+tfbd->filterLength)] = kick;
  
    if (pass<tfbd->delay+tfbd->filterLength) {
#ifdef DEBUG
      fprintf(stdout, "TFBDRIVER: no kick applied for pass %ld due to delay of %ld\n",
              pass, tfbd->delay);
#endif
      kick = 0;
    }
    else {
      kick = 0;
      for (i=0; i<tfbd->filterLength; i++) {
#ifdef DEBUG
        fprintf(stdout, "TFBDRIVER: adding term a[%ld]=%e  *   %e\n",
                i, tfbd->a[i], tfbd->driverSignal[iBucket][(pass - tfbd->delay - i)%(tfbd->delay+tfbd->filterLength)]);
        fflush(stdout);
#endif
        kick += tfbd->a[i]*tfbd->driverSignal[iBucket][(pass - tfbd->delay - i)%(tfbd->delay+tfbd->filterLength)];
      }
#ifdef DEBUG
      fprintf(stdout, "TFBDRIVER: kick = %le\n", kick);
      fflush(stdout);
#endif

      if (isSlave || !notSinglePart) {
        j = tfbd->pickup->yPlane?3:1;
        if (nBuckets==1) {
          for (i=0; i<np0; i++)
            part0[i][j] += kick;
        } else {
          for (i=0; i<npBucket[iBucket]; i++) {
            part0[ipBucket[iBucket][i]][j] += kick;
          }
        }
      }
      
#ifdef DEBUG
      fprintf(stdout, "TFBDRIVER: kick applied for bunch %ld\n", iBucket);
      fflush(stdout);
#endif
    }
    
#if USE_MPI
    if (myid==0) 
#endif
      if (tfbd->outputFile) {
        if (!SDDS_SetRowValues(&tfbd->SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, 
                               tfbd->outputIndex++, "Bunch", iBucket, "Pass", pass,
                               "PickupOutput", tfbd->pickup->filterOutput[iBucket], 
                               "DriverOutput", kick, NULL)) {
          SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors);
          SDDS_Bomb("problem writing data for TFBDRIVER output file");
        }
        if ((pass+1)==nPasses) {
          if (!SDDS_WritePage(&tfbd->SDDSout)) {
            SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors);
            SDDS_Bomb("problem writing data for TFBDRIVER output file");
          }
          tfbd->outputIndex = 0;
          tfbd->dataWritten = 1;
        }
      }
  }
  
  if (time0) 
    free(time0);
  if (ibParticle) 
    free(ibParticle);
  if (ipBucket)
    free_czarray_2d((void**)ipBucket, nBuckets, np0);
  if (npBucket)
    free(npBucket);

}

void flushTransverseFeedbackDriverFiles(TFBDRIVER *tfbd)
{
#if USE_MPI
  if (myid!=0)
    return;
#endif
  if (tfbd->initialized && !(tfbd->dataWritten)) {
    if (!SDDS_WritePage(&tfbd->SDDSout)) {
      SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors);
      SDDS_Bomb("problem writing data for TFBDRIVER output file (flushTransverseFeedbackDriverFiles)");
    }
    tfbd->dataWritten = 1;
  }
  tfbd->outputIndex = 0;
}

void initializeTransverseFeedbackDriver(TFBDRIVER *tfbd, LINE_LIST *beamline, long nPasses, char *rootname)
{
  ELEMENT_LIST *eptr;
  long pickupFound = 0, i;

  if (tfbd->ID==NULL || !strlen(tfbd->ID))
    bombElegant("you must give an ID string for TFBDRIVER", NULL);
  
  for (i=TFB_FILTER_LENGTH-1; i>=0; i--) {
    if (tfbd->a[i]!=0)
      break;
  }
  if (i<0)
    bombElegant("All filter coefficients are zero for TFBDRIVER", NULL);
  tfbd->filterLength = i+1;

  eptr = &(beamline->elem);
  while (eptr) {
    if (eptr->type==T_TFBPICKUP && strcmp(tfbd->ID, ((TFBPICKUP*)eptr->p_elem)->ID)==0) {
      pickupFound = 1;
      tfbd->pickup = ((TFBPICKUP*)eptr->p_elem);
      break;
    }
    eptr = eptr->succ;
  }
  if (!pickupFound) 
    bombElegant("pickup not found for TFBDRIVER", NULL);
  
  if (tfbd->delay<0)
    bombElegant("TFBDRIVER delay is negative", NULL);

#if USE_MPI
  if (myid==0)
#endif
  if (tfbd->outputFile) {
    tfbd->outputFile = compose_filename(tfbd->outputFile, rootname);
    if (!SDDS_InitializeOutput(&tfbd->SDDSout, SDDS_BINARY, 1, NULL, NULL, tfbd->outputFile) ||
        !SDDS_DefineSimpleColumn(&tfbd->SDDSout, "Pass", NULL, SDDS_LONG) ||
        !SDDS_DefineSimpleColumn(&tfbd->SDDSout, "Bunch", NULL, SDDS_LONG) ||
        !SDDS_DefineSimpleColumn(&tfbd->SDDSout, "PickupOutput", NULL, SDDS_DOUBLE) ||
        !SDDS_DefineSimpleColumn(&tfbd->SDDSout, "DriverOutput", "rad", SDDS_DOUBLE) ||
        !SDDS_WriteLayout(&tfbd->SDDSout) || !SDDS_StartPage(&tfbd->SDDSout, nPasses)) {
      SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors);
      SDDS_Bomb("Problem setting up TFBDRIVER output file");
    }
  }
  tfbd->dataWritten = tfbd->outputIndex = 0;

  if (tfbd->driverSignal) {
    long i;
    for (i=0;i<tfbd->nBunches; i++)
      free(tfbd->driverSignal[i]);
    free(tfbd->driverSignal);
    tfbd->driverSignal = NULL;
  }
  tfbd->nBunches = 0;

  tfbd->initialized = 1;
}


