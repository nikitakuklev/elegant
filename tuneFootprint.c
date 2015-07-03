/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: tuneFootprint.c
 * purpose: Do frequency map tracking and analysis to determine the tune footprint.
 *          See file tuneFootprint.nl for input parameters.
 *
 * Michael Borland, 2004
 */
#include "mdb.h"
#include "track.h"
#include "tuneFootprint.h"

typedef struct {
  short ix, iy, used;
  double diffusionRate, position[2], nu[2];
} XY_TF_DATA;

typedef struct {
  short idelta, used;
  double diffusionRate, delta, nu[2];
} DELTA_TF_DATA;

int xyTfDataCompare(const void *xy1c, const void *xy2c)
{
  short dix, diy;
  XY_TF_DATA *xy1, *xy2;
  xy1 = (XY_TF_DATA*)xy1c;
  xy2 = (XY_TF_DATA*)xy2c;
  if (xy2->used==0)
    return -1; /* unused slot is always greater than anything else */
  if (xy1->used==0) 
    return 1; /* unused slot is always greater than anything else */
  dix = xy1->ix - xy2->ix;
  diy = xy1->iy - xy2->iy;
  if (diy!=0)
    return diy;
  return dix;
}

int deltaTfDataCompare(const void *delta1c, const void *delta2c)
{
  DELTA_TF_DATA *delta1, *delta2;
  delta1 = (DELTA_TF_DATA*)delta1c;
  delta2 = (DELTA_TF_DATA*)delta2c;
  if (delta2->used==0)
    return -1; /* unused slot is always greater than anything else */
  if (delta1->used==0)
    return 1; /* unused slot is always greater than anything else */
  return delta1->idelta - delta2->idelta;
}

void determineDeltaTuneFootprint(DELTA_TF_DATA *deltaTfData, long nDelta, double *tuneRange, double *deltaRange, double *diffusionRateMax, double *nuxLimit, double *nuyLimit)
{
  long id, id0, id1, id2, coord;
  double delta0, delta1, delta2, nu0[2];
  double nuMin, nuMax;

  /* Find dividing point between positive and negative delta */
  id0 = -1;
  for (id=0; id<nDelta; id++) {
    if (deltaTfData[id].delta>=0) {
      id0 = id;
      nu0[0] = deltaTfData[id].nu[0];
      nu0[1] = deltaTfData[id].nu[1];
      delta0 = deltaTfData[id].delta;
      break;
    }
  }
  if (id0==-1 || 
      (deltaTfData[id].nu[0]<=0 || deltaTfData[id].nu[0]>=1 ||
       deltaTfData[id].nu[1]<=0 || deltaTfData[id].nu[1]>=1 || deltaTfData[id].diffusionRate>diffusion_rate_limit)) {
    tuneRange[0] = tuneRange[1] = deltaRange[0] = deltaRange[1]= 0;
    return;
  }
#ifdef DEBUG
  printf("id0 = %ld, nux = %le, nuy = %le\n", id0, deltaTfData[id0].nu[0], deltaTfData[id0].nu[1]);
#endif

  for (coord=0; coord<2; coord++) {
    id1 = id2 = -1;
    delta1 = delta2 = 0; 
    for (id=id0; id>=0; id--) {
      /* scan to negative delta side from center until we find a place where the 
       * tune is invalid or crosses a half integer */
      if (deltaTfData[id].nu[coord]<=0 || deltaTfData[id].nu[coord]>=1 || 
          (!ignore_half_integer && ((long)(2*deltaTfData[id].nu[coord])!=(long)(2*nu0[coord]))) || 
          (compute_diffusion && deltaTfData[id].diffusionRate>diffusion_rate_limit)) {
#ifdef DEBUG
        printf("coord = %ld, bad: id = %ld, delta = %le, nu = %le, nu0 = %le, dR = %le\n",
             coord, id, deltaTfData[id].delta, deltaTfData[id].nu[coord], nu0[coord], deltaTfData[id].diffusionRate);
#endif
        break;
      }
#ifdef DEBUG
      printf("coord = %ld, ok: id = %ld, delta = %le, nu = %le, nu0 = %le, dR = %le\n",
             coord, id, deltaTfData[id].delta, deltaTfData[id].nu[coord], nu0[coord], deltaTfData[id].diffusionRate);
#endif
      delta1 = deltaTfData[id].delta;
      id1 = id;
    }
    for (id=id0; id<nDelta; id++) {
      /* scan to positive delta side from center until we find a place where the 
       * tune is invalid or crosses a half integer */
      if (deltaTfData[id].nu[coord]<=0 || deltaTfData[id].nu[coord]>=1 || 
          (!ignore_half_integer && ((long)(2*deltaTfData[id].nu[coord])!=(long)(2*nu0[coord]))) ||
          (compute_diffusion && deltaTfData[id].diffusionRate>diffusion_rate_limit)) {
#ifdef DEBUG
        printf("coord = %ld, bad: id = %ld, delta = %le, nu = %le, nu0 = %le, dR = %le\n",
             coord, id, deltaTfData[id].delta, deltaTfData[id].nu[coord], nu0[coord], deltaTfData[id].diffusionRate);
#endif
        break;
      }
#ifdef DEBUG
      printf("coord = %ld, ok: id = %ld, delta = %le, nu = %le, nu0 = %le, dR = %le\n",
             coord, id, deltaTfData[id].delta, deltaTfData[id].nu[coord], nu0[coord], deltaTfData[id].diffusionRate);
#endif
      delta2 = deltaTfData[id].delta;
      id2 = id;
    }
    deltaRange[coord] = MIN(fabs(delta1), fabs(delta2));
#ifdef DEBUG
    printf("coord = %ld, delta1 = %le, id1 = %ld, delta2 = %le, id2 = %ld, deltaRange = %le\n",
           coord, delta1, id1, delta2, id2, deltaRange[coord]);
#endif

    for (id=0; id<id1; id++)
      deltaTfData[id].used = 0;
    for (id=id2+1; id<nDelta; id++)
      deltaTfData[id].used = 0;

    nuMin = 1;
    nuMax = 0;
    for (id=id1; id<=id2; id++) {
      if (deltaTfData[id].nu[coord]<nuMin)
        nuMin = deltaTfData[id].nu[coord];
      if (deltaTfData[id].nu[coord]>nuMax)
        nuMax = deltaTfData[id].nu[coord];
    }
    if (coord==0) {
      if (nuxLimit) {
	nuxLimit[0] = nuMin;
	nuxLimit[1] = nuMax;
      }
    } else {
      if (nuyLimit) {
	nuyLimit[0] = nuMin;
	nuyLimit[1] = nuMax;
      }
    }
    tuneRange[coord] = nuMax - nuMin;
  }

  *diffusionRateMax = -DBL_MAX;
  for (id1=0; id<=id2; id++) {
    if (*diffusionRateMax < deltaTfData[id].diffusionRate)
      *diffusionRateMax = deltaTfData[id].diffusionRate;
  }

}

void determineXyTuneFootprint(XY_TF_DATA *xyTfData, long nx, long ny, double *tuneRange, double *positionRange, double *diffusionRateMax, double *area,
			      double *nuxLimit, double *nuyLimit)
{
  long id;
  double nuMin, nuMax, pMin, pMax;
  double distance, bestDistance;
  long ix0, ix1, ix2, ix, iy, iy1, iy2, ic;
  
  /* check for indexing issues */
  for (ix=0; ix<nx; ix++) {
    for (iy=0; iy<ny; iy++) {
      id = ix + nx*iy;
      if (ix!=xyTfData[id].ix || xyTfData[id].iy!=iy) {
        fprintf(stderr, "indexing problem in determineXyTuneFootprint: ix=%ld, iy=%ld, id=%ld, ix[id]=%d, iy[id]=%d\n",
                ix, iy, id, xyTfData[id].ix, xyTfData[id].iy);
        bombElegant("Indexing bug", NULL);
      }
      if (xyTfData[id].nu[0]<=0 || xyTfData[id].nu[0]>=1 ||
          xyTfData[id].nu[1]<=0 || xyTfData[id].nu[1]>=1 || 
	  (compute_diffusion && xyTfData[id].diffusionRate>diffusion_rate_limit))
        xyTfData[id].used = 0;
    }
  }

  /* Find point closest to origin (iy=0) */
  bestDistance = DBL_MAX;
  ix0 = -1;
  for (ix=0; ix<nx; ix++) {
    if (xyTfData[ix].used) {
      distance = fabs(xyTfData[ix].position[0]);
      if (distance<bestDistance) {
        ix0 = ix;
        bestDistance = distance;
      }
    }
  }
  if (ix0==-1) {
    tuneRange[0] = tuneRange[1] = positionRange[0] = positionRange[1]= 0;
    return;
  }

  for (iy=0; iy<ny; iy++) {
    /* find limiting indices ix1 and ix2 of valid data for iy */
    ix1 = ix2 = ix0;
    for (ix=ix0; ix>=0; ix--) {
      id = ix + nx*iy;
      if (!xyTfData[id].used)
        break;
      ix1 = ix;
    }
    for (ix=ix0; ix<nx; ix++) {
      id = ix + nx*iy;
      if (!xyTfData[id].used)
        break;
      ix2 = ix;
    }
    /* Mark all points at or above this iy value with ix<ix1 or ix>ix2 as bad */
    for (iy1=iy; iy1<ny; iy1++) {
      for (ix=0; ix<ix1; ix++)
        xyTfData[ix + nx*iy1].used = 0;
      for (ix=ix2+1; ix<nx; ix++)
        xyTfData[ix + nx*iy1].used = 0;
    }
  }

  /* find the area for points with acceptable characteristics */
  for (ix1=0; ix1<nx; ix1++) {
    if (xyTfData[ix1].used)
      break;
  }
  for (ix2=nx-1; ix2>ix1; ix2--) {
    if (xyTfData[ix2].used)
      break;
  }
  iy1 = 0;
  for (iy=0; iy<ny; iy++) {
    id = ix + nx*iy;
    if (!xyTfData[id].used)
      break;
    iy1 = iy;
  }
  *area = 0;
  for (ix=ix1+1; ix<=ix2; ix++) {    
    iy2 = 0;
    for (iy=0; iy<ny; iy++) {
      id = ix + nx*iy;
      if (!xyTfData[id].used)
        break;
      iy2 = iy;
    }
    *area += (xyTfData[ix + nx*iy2].position[1] + xyTfData[ix-1 + nx*iy1].position[1])*
      (xyTfData[ix].position[0] - xyTfData[ix-1].position[0])/2;
    iy1 = iy2;
  }
  
  /* find the spread in tune and position for points with acceptable characteristics */
  for (ic=0; ic<2; ic++) {
    nuMin = 1;
    nuMax = 0;
    pMax = -(pMin = DBL_MAX);
    for (ix=0; ix<nx; ix++) {
      for (iy=0; iy<ny; iy++) {
        id = ix + nx*iy;
        if (xyTfData[id].used) {
          if (xyTfData[id].nu[ic]<nuMin)
            nuMin = xyTfData[id].nu[ic];
          if (xyTfData[id].nu[ic]>nuMax)
            nuMax = xyTfData[id].nu[ic];
          if (xyTfData[id].position[ic]<pMin)
            pMin = xyTfData[id].position[ic];
          if (xyTfData[id].position[ic]>pMax)
            pMax = xyTfData[id].position[ic];
        }
      }
    }
    tuneRange[ic] = nuMax - nuMin;
    if (ic==0) {
      if (nuxLimit) {
	nuxLimit[0] = nuMin;
	nuxLimit[1] = nuMax;
      }
    } else {
      if (nuyLimit) {
	nuyLimit[0] = nuMin;
	nuyLimit[1] = nuMax;
      }
    }

    positionRange[ic] = pMax - pMin;
  }

  /* Find the maximum diffusion rate for points with acceptable characteristics */
  *diffusionRateMax = -DBL_MAX;
  for (ix=0; ix<nx; ix++) {
    for (iy=0; iy<ny; iy++) {
      id = ix + nx*iy;
      if (xyTfData[id].used && *diffusionRateMax < xyTfData[id].diffusionRate)
        *diffusionRateMax = xyTfData[id].diffusionRate;
    }
  }
}

#if USE_MPI
MPI_Datatype xyTfDataType, deltaTfDataType, tfReturnDataType;
void setupTuneFootprintDataTypes ()
{
  MPI_Datatype oldType[5];
  int blockLength[5];
  MPI_Aint offset[5];
  XY_TF_DATA xyTfExample;
  DELTA_TF_DATA deltaTfExample;
  TUNE_FOOTPRINTS tfExample;
  long i;
  
  oldType[0] = MPI_SHORT;
  oldType[1] = MPI_DOUBLE;

  blockLength[0] = 3;
  blockLength[1] = 5;
  MPI_Get_address(&xyTfExample.ix, &offset[0]);
  MPI_Get_address(&xyTfExample.diffusionRate, &offset[1]);
  for (i=1; i>=0; i--)
    offset[i] -= offset[0];
  MPI_Type_create_struct(2, blockLength, offset, oldType, &xyTfDataType);
  MPI_Type_commit(&xyTfDataType);
  
  blockLength[0] = 2;
  blockLength[1] = 4;
  MPI_Get_address(&deltaTfExample.idelta, &offset[0]);
  MPI_Get_address(&deltaTfExample.diffusionRate, &offset[1]);
  for (i=1; i>=0; i--)
    offset[i] -= offset[0];
  MPI_Type_struct(2, blockLength, offset, oldType, &deltaTfDataType);
  MPI_Type_commit(&deltaTfDataType);

  oldType[0] = oldType[1] = oldType[2] = oldType[3] = oldType[4] = MPI_DOUBLE;
  blockLength[0] = 2;
  blockLength[1] = 3;
  blockLength[2] = 2;
  blockLength[3] = 2;
  blockLength[4] = 3;
  MPI_Get_address(&tfExample.chromaticTuneRange, &offset[0]);
  MPI_Get_address(&tfExample.deltaRange, &offset[1]);
  MPI_Get_address(&tfExample.amplitudeTuneRange, &offset[2]);
  MPI_Get_address(&tfExample.positionRange, &offset[3]);
  MPI_Get_address(&tfExample.chromaticDiffusionMaximum, &offset[4]);
  for (i=4; i>=0; i--)
    offset[i] -= offset[0];
  MPI_Type_struct(5, blockLength, offset, oldType, &tfReturnDataType);
  MPI_Type_commit(&tfReturnDataType);
}
#endif

static short tuneFootprintOn = 0;

long setupTuneFootprint(
    NAMELIST_TEXT *nltext,
    RUN *run,
    VARY *control
    )
{
  /* process namelist input */
  set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
  set_print_namelist_flags(0);
  if (processNamelist(&tune_footprint, nltext)==NAMELIST_ERROR)
    bombElegant(NULL, NULL);
  if (echoNamelists) print_namelist(stdout, &tune_footprint);
  
  /* check for data errors */
  if (xmin>xmax)
    bombElegant("xmin > xmax", NULL);
  if (ymin>ymax)
    bombElegant("ymin > ymax", NULL);
  if (ymin<0)
    ymin = 0;
  if (delta_min>delta_max)
    bombElegant("delta_min > delta_max", NULL);
  if (quadratic_spacing) {
    if (xmin<0)
      xmin = 0;
    if (ymin<0)
      ymin = 0;
    if (delta_min<0)
      delta_min = 0;
  }

  if (delta_output)
    delta_output = compose_filename(delta_output, run->rootname);
  if (xy_output)
    xy_output = compose_filename(xy_output, run->rootname);

  tuneFootprintOn = 1;
  return immediate;
}

/* Used to save results of last run for output after optimization */
static XY_TF_DATA *allXyTfData = NULL;
static DELTA_TF_DATA *allDeltaTfData = NULL;

long doTuneFootprint(
                    RUN *run,
                    VARY *control,
                    double *referenceCoord,
                    LINE_LIST *beamline,
                    TUNE_FOOTPRINTS *tfReturn
                    )
{
  double firstTune[2], secondTune[2], startingCoord[6], endingCoord[6];
  double firstAmplitude[2], secondAmplitude[2];
  double dx, dy, ddelta, x, y, delta;
  long ix, iy, ixy, idelta, turns;
  long my_ixy, my_idelta;
  static double **one_part;
  double p, oldPercentage=0;
  long n_part, lost;
  XY_TF_DATA *xyTfData;
  DELTA_TF_DATA *deltaTfData;
  long my_nxy, my_ndelta;
  double chromTuneRange[2], chromDeltaRange[2], xyTuneRange[2], xyPositionRange[2], diffusionRateMax, xyArea;
  double nuxLimit[2], nuyLimit[2];

#ifdef DEBUG
  FILE *fpdebug = NULL;
#endif
#if USE_MPI 
  MPI_Status mpiStatus;
#endif

  memset(tfReturn, 0, sizeof(*tfReturn));

  if (!tuneFootprintOn)
    return 0;
  if (allXyTfData) {
    free(allXyTfData);
    allXyTfData = NULL;
  }
  if (allDeltaTfData) {
    free(allDeltaTfData);
    allDeltaTfData = NULL;
  }
  
#if USE_MPI
  setupTuneFootprintDataTypes();
  /* Note that the master is a working processor for this algorithm */
  if (ndelta) {
    if (ndelta<n_processors) {
      ndelta = n_processors;
      if (myid==0) {
	printf("NB: ndelta increased to equal the number of working processors (%d)\n", n_processors);
      }
    }
    if (ndelta%n_processors!=0) {
      /* ensure that all processors have an equal number of particles to work on */
      ndelta = (ndelta/n_processors+1)*n_processors;
      if (myid==0) {
	printf("NB: ndelta increased to equal a multiple (%ld) of the number of working processors (%d)\n", ndelta, n_processors);
      }
    }
  }
  if (nx && ny) {
    if (nx*ny<n_processors && myid==0) {
      printf("NB: number of x, y grid points less than the number of cores, which is a waste of resources\n");
    }
  }
  my_nxy = (nx*ny)/n_processors+1;
  my_ndelta = ndelta/n_processors+1;
#else
  my_nxy = nx*ny;
  my_ndelta = ndelta;
#endif
  
  if (my_nxy) 
    xyTfData = calloc(my_nxy, sizeof(*xyTfData));
  else
    xyTfData = NULL;
  if (my_ndelta)
    deltaTfData = calloc(my_ndelta, sizeof(*deltaTfData));
  else
    deltaTfData = NULL;

  /* Perform fiducialization by tracking one turn */
  if (!one_part)
    one_part = (double**)czarray_2d(sizeof(**one_part), 1, 7);
  n_part = 1;
  if (referenceCoord) {
    long i;
    for (i=0; i<6; i++)
      one_part[0][i] = referenceCoord[i];
  }
  p = run->p_central;
  if (!do_tracking(NULL, one_part, n_part, NULL, beamline, &p, (double**)NULL, (BEAM_SUMS**)NULL, (long*)NULL,
                   NULL, run, 0, TEST_PARTICLES, 1, 0,
                   NULL, NULL, NULL, NULL, NULL)) {
    printf("Error: lost particle when fiducializing\n");
    exitElegant(1);
  }

  turns = control->n_passes/2;

  if (ndelta) {
    /* First perform delta scan */
    x = x_for_delta;
    y = y_for_delta;
    if (ndelta>1)
      ddelta = (delta_max-delta_min)/(ndelta-1);
    else
      ddelta = 0;
    for (idelta=my_idelta=0; idelta<ndelta; idelta++) {
      deltaTfData[my_idelta].used = 0;
      if (!quadratic_spacing)
	delta = delta_min + idelta*ddelta;
      else {
	if (idelta<(ndelta-1.)/2)
	  delta = -((delta_max-delta_min)*sqrt(fabs((idelta-(ndelta-1)/2.)/((ndelta-1)/2.)))+delta_min);
	else
	  delta = ((delta_max-delta_min)*sqrt(fabs((idelta-(ndelta-1)/2.)/((ndelta-1)/2.)))+delta_min);
      }
      memcpy(startingCoord, referenceCoord, sizeof(*startingCoord)*6);
#if USE_MPI
      if (myid == idelta%n_processors) /* Partition the job according to particle ID */
#endif
	{
	  if (verbosity>=2) {
#if USE_MPI
	    if (myid==0)
#endif
	      printf("computing tune for delta = %le\n", delta);
	  }
	  lost = 0;
	  if (!computeTunesFromTracking(firstTune, firstAmplitude,
					beamline->matrix, beamline, run,
					startingCoord, x, y, delta, turns, 0,
					0, endingCoord, NULL, NULL, 1, 1) ||
	      firstTune[0]>1.0 || firstTune[0]<0 || firstTune[1]>1.0 || firstTune[1]<0) {
	    lost = 1;
	  } else if (compute_diffusion) {
	    memcpy(startingCoord, endingCoord, sizeof(*startingCoord)*6);
	    if (!computeTunesFromTracking(secondTune, secondAmplitude,
					  beamline->matrix, beamline, run,
					  startingCoord, 0.0, 0.0, 0.0, turns, turns,
					  0, endingCoord, NULL, NULL, 1, 1) || 
		secondTune[0]>1.0 || secondTune[0]<0 || secondTune[1]>1.0 || secondTune[1]<0) {
	      lost = 1;
	    }
	  }             
#if USE_MPI
	  if (my_idelta>=my_ndelta) {
	    fprintf(stderr, "delta index too large on processor %d\n", myid);
	    exit(1);
	  }
#endif
	  deltaTfData[my_idelta].idelta = idelta;
	  deltaTfData[my_idelta].delta = delta;
	  deltaTfData[my_idelta].used = 1;
	  if (lost) {
	    deltaTfData[my_idelta].nu[0] = deltaTfData[my_idelta].nu[1] = -2;
	    deltaTfData[my_idelta].diffusionRate = DBL_MAX;
	  } else {
	    deltaTfData[my_idelta].nu[0] = firstTune[0];
	    deltaTfData[my_idelta].nu[1] = firstTune[1];
	    if (compute_diffusion)
	      deltaTfData[my_idelta].diffusionRate = log10((sqr(secondTune[0] - firstTune[0]) + sqr(secondTune[1] - firstTune[1]))/turns);
	    else
	      deltaTfData[my_idelta].diffusionRate = -DBL_MAX;
	  }
	  my_idelta ++;
	  if (verbosity>=2) {
#if USE_MPI
	    if (myid==0) {
	      double newPercentage = (100.0*idelta)/ndelta;
	      if ((newPercentage-oldPercentage)>=1) {
		fprintf(stdout, "About %.1f%% done with energy scan\n", newPercentage);
		oldPercentage = newPercentage;
		fflush(stdout);
	      }
	    }
#else
	    fprintf(stdout, "Done with particle %ld of %ld for energy scan\n",
		    idelta+1, ndelta);
	    fflush(stdout);
#endif
	  }
	}
    }

#if USE_MPI
    /* Collect data onto the master */
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid==0) {
      long id, iTotal;
      MPI_Status mpiStatus;
      allDeltaTfData  = calloc(my_ndelta*n_processors, sizeof(*allDeltaTfData));
      memcpy(allDeltaTfData, deltaTfData, sizeof(*deltaTfData)*my_ndelta);
      /* receive data */
      iTotal = my_ndelta;
      for (id=1; id<n_processors; id++) {
	if (MPI_Recv(allDeltaTfData+iTotal, my_ndelta, deltaTfDataType, id, 1, MPI_COMM_WORLD, &mpiStatus)!=MPI_SUCCESS) {
	  printf("Error: MPI_Recv returns error retrieving data from processor %ld\n", id);
	  bombElegant("Communication error", NULL);
	}
	iTotal += my_ndelta;
      }
      my_ndelta = iTotal;
    } else {
      /* send data */
      MPI_Send(deltaTfData, my_ndelta, deltaTfDataType, 0, 1, MPI_COMM_WORLD);
      allDeltaTfData = NULL;
    }
    MPI_Barrier(MPI_COMM_WORLD);
#else
    allDeltaTfData = deltaTfData;
#endif
    
#if USE_MPI
    if (myid==0) {
#endif
      qsort(allDeltaTfData, my_ndelta, sizeof(*allDeltaTfData), deltaTfDataCompare);
      for (idelta=0; idelta<my_ndelta; idelta++) {
	if (allDeltaTfData[idelta].used==0) {
	  my_ndelta = idelta;
	  break;
	}
      }
#ifdef DEBUG
      fpdebug = fopen("tfDelta.sdds", "w");
      fprintf(fpdebug, "SDDS1\n&column name=idelta type=short &end\n");
      fprintf(fpdebug, "&column name=delta type=double &end\n");
      fprintf(fpdebug, "&column name=nux type=double &end\n");
      fprintf(fpdebug, "&column name=nuy type=double &end\n");
      fprintf(fpdebug, "&column name=diffusionRate type=double &end\n");
      fprintf(fpdebug, "&data mode=ascii no_row_counts=1 &end\n");
      for (idelta=0; idelta<my_ndelta; idelta++) {
	fprintf(fpdebug, "%d %e %e %e %e\n",
		allDeltaTfData[idelta].idelta, 
		allDeltaTfData[idelta].delta,
		allDeltaTfData[idelta].nu[0], allDeltaTfData[idelta].nu[1],
		allDeltaTfData[idelta].diffusionRate);
      }
      fclose(fpdebug);
#endif
      
      determineDeltaTuneFootprint(allDeltaTfData, my_ndelta, chromTuneRange, chromDeltaRange, &diffusionRateMax, nuxLimit, nuyLimit);
      if (verbosity)
	printf("nux/chromatic: tune range=%le (%le, %le) delta range = %le\nnuy/chromatic: tune range=%le (%le, %le) delta range = %le\n",
	       chromTuneRange[0], nuxLimit[0], nuxLimit[1], chromDeltaRange[0],
	       chromTuneRange[1], nuyLimit[0], nuyLimit[1], chromDeltaRange[1]);
      if (tfReturn) {
	memcpy(tfReturn->chromaticTuneRange, chromTuneRange, sizeof(*chromTuneRange)*2);
	memcpy(tfReturn->deltaRange, chromDeltaRange, sizeof(*chromDeltaRange)*2);
	memcpy(tfReturn->nuxChromLimit, nuxLimit, sizeof(*nuxLimit)*2);
	memcpy(tfReturn->nuyChromLimit, nuyLimit, sizeof(*nuyLimit)*2);
	tfReturn->chromaticDiffusionMaximum = diffusionRateMax;
	tfReturn->deltaRange[2] = MIN(tfReturn->deltaRange[0], tfReturn->deltaRange[1]);
      }
#if USE_MPI
    }
#endif
  }

  if (nx!=0 && ny!=0) {
    dx = dy = delta = 0;
    if (!quadratic_spacing) {
      if (nx>1)
	dx  = (xmax-xmin)/(nx-1);
      if (ny>1)
	dy = (ymax-ymin)/(ny-1);
    }
    my_ixy = 0;
    oldPercentage = 0;
    for (ix=0; ix<nx; ix++) {
      if (quadratic_spacing) {
	if (ix<nx/2) 
	  x = -((xmax-xmin)*sqrt(fabs((ix-(nx-1)/2.)/((nx-1)/2.))) + xmin);
	else
	  x = ((xmax-xmin)*sqrt(fabs((ix-(nx-1)/2.)/((nx-1)/2.))) + xmin);
      } else {
	x = xmin + ix*dx;
      }
      for (iy=0; iy<ny; iy++) {
	if (quadratic_spacing) {
	  y = (ymax-ymin)*sqrt(iy/(ny-1.)) + ymin;
	} else {
	  y = ymin + iy*dy;
	}
	xyTfData[my_ixy].used = 0;
	memcpy(startingCoord, referenceCoord, sizeof(*startingCoord)*6);
#if USE_MPI
	if (myid == (ix*ny+iy)%n_processors) /* Partition the job according to particle ID */
#endif
	  {
            lost = 0;
	    if (!computeTunesFromTracking(firstTune, firstAmplitude,
					  beamline->matrix, beamline, run,
					  startingCoord, x, y, delta, turns, 0,
					  0, endingCoord, NULL, NULL, 1, 1) ||
		firstTune[0]>1.0 || firstTune[0]<0 || firstTune[1]>1.0 || firstTune[1]<0) {
              lost = 1;
            } else if (compute_diffusion) {
              memcpy(startingCoord, endingCoord, sizeof(*startingCoord)*6);
              if (!computeTunesFromTracking(secondTune, secondAmplitude,
                                            beamline->matrix, beamline, run,
                                            startingCoord, 0.0, 0.0, 0.0, turns, turns,
                                            0, endingCoord, NULL, NULL, 1, 1) || 
                  secondTune[0]>1.0 || secondTune[0]<0 || secondTune[1]>1.0 || secondTune[1]<0) {
                lost = 1;
              }
            }
            xyTfData[my_ixy].ix = ix;
            xyTfData[my_ixy].iy = iy;
            xyTfData[my_ixy].position[0] = x;
            xyTfData[my_ixy].position[1] = y;
            xyTfData[my_ixy].used = 1;
            if (lost) {
              xyTfData[my_ixy].nu[0] = xyTfData[my_ixy].nu[1] = -2;
              xyTfData[my_ixy].diffusionRate = DBL_MAX;
            } else {
              xyTfData[my_ixy].nu[0] = firstTune[0];
              xyTfData[my_ixy].nu[1] = firstTune[1];
	      if (compute_diffusion)
		xyTfData[my_ixy].diffusionRate = log10((sqr(secondTune[0] - firstTune[0]) + sqr(secondTune[1] - firstTune[1]))/turns);
	      else
		xyTfData[my_ixy].diffusionRate = -DBL_MAX;
            }
            my_ixy ++;

	    if (verbosity>=2) {
#if USE_MPI
	      if (myid==0) {
		double newPercentage = 100*(ix*ny+iy+1.0)/(nx*ny);
		if ((newPercentage-oldPercentage)>=5) {
		  fprintf(stdout, "About %.1f%% done with x, y scan\n", newPercentage);
		  oldPercentage = newPercentage;
		  fflush(stdout);
		}
	      }
#else
	      fprintf(stdout, "Done with particle %ld of %ld\n",
		      ix*ny+iy+1, nx*ny);
	      fflush(stdout);
#endif
	    }
	  }
      }
    }
    
#if USE_MPI
    /* Collect data onto the master */
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid==0) {
      long id, iTotal;
      MPI_Status mpiStatus;
      allXyTfData  = calloc(my_nxy*n_processors, sizeof(*allXyTfData));
      memcpy(allXyTfData, xyTfData, sizeof(*xyTfData)*my_nxy);
      /* receive data */
      iTotal = my_nxy;
      for (id=1; id<n_processors; id++) {
	if (MPI_Recv(allXyTfData+iTotal, my_nxy, xyTfDataType, id, 1, MPI_COMM_WORLD, &mpiStatus)!=MPI_SUCCESS) {
	  printf("Error: MPI_Recv returns error retrieving data from processor %ld\n", id);
	  bombElegant("Communication error", NULL);
	}
	iTotal += my_nxy;
      }
      my_nxy = iTotal;
    } else {
      /* send data */
      MPI_Send(xyTfData, my_nxy, xyTfDataType, 0, 1, MPI_COMM_WORLD);
      allXyTfData = NULL;
    }
    MPI_Barrier(MPI_COMM_WORLD);
#else
    allXyTfData = xyTfData;
#endif

#if USE_MPI
    if (myid==0) {
#endif    
      qsort(allXyTfData, my_nxy, sizeof(*allXyTfData), xyTfDataCompare);
      my_ixy = 0;
      for (ixy=0; ixy<my_nxy; ixy++) {
	if (!allXyTfData[ixy].used) {
	  my_nxy = ixy;
	  break;
	}
      }
      if (my_nxy!=nx*ny) {
	fprintf(stderr, "my_nxy = %ld, nx*ny = %ld\n", my_nxy, nx*ny);
	bombElegant("Error: counting for nx*ny didn't work", NULL);
      }
      
      determineXyTuneFootprint(allXyTfData, nx, ny, xyTuneRange, xyPositionRange, &diffusionRateMax, &xyArea, nuxLimit, nuyLimit);
      if (verbosity)
	printf("nux/amplitude: tune range=%le  position range = %le\nnuy/amplitude: tune range=%le position range = %le\nArea = %le\n",
	       xyTuneRange[0], xyPositionRange[0],
	       xyTuneRange[1], xyPositionRange[1], xyArea);
      if (tfReturn) {
	memcpy(tfReturn->amplitudeTuneRange, xyTuneRange, sizeof(*xyTuneRange)*2);
	memcpy(tfReturn->positionRange, xyPositionRange, sizeof(*xyPositionRange)*2);
	memcpy(tfReturn->nuxAmpLimit, nuxLimit, sizeof(*nuxLimit)*2);
	memcpy(tfReturn->nuyAmpLimit, nuyLimit, sizeof(*nuyLimit)*2);
	tfReturn->xyArea = xyArea;
	tfReturn->amplitudeDiffusionMaximum = diffusionRateMax;
      }
      
#ifdef DEBUG
      fpdebug = fopen("tfXy.sdds", "w");
      fprintf(fpdebug, "SDDS1\n");
      fprintf(fpdebug, "&column name=ix type=short &end\n");
      fprintf(fpdebug, "&column name=iy type=short &end\n");
      fprintf(fpdebug, "&column name=x type=double &end\n");
      fprintf(fpdebug, "&column name=y type=double &end\n");
      fprintf(fpdebug, "&column name=nux type=double &end\n");
      fprintf(fpdebug, "&column name=nuy type=double &end\n");
      fprintf(fpdebug, "&column name=diffusionRate type=double &end\n");
      fprintf(fpdebug, "&data mode=ascii no_row_counts=1 &end\n");
      for (ix=0; ix<nx; ix++) {
	for (iy=0; iy<ny; iy++) {
	  ixy = ix + iy*nx;
	  fprintf(fpdebug, "%d %d %e %e %e %e %e\n",
		  allXyTfData[ixy].ix,
		  allXyTfData[ixy].iy,
		  allXyTfData[ixy].position[0], allXyTfData[ixy].position[1], 
		  allXyTfData[ixy].nu[0], allXyTfData[ixy].nu[1], 
		  allXyTfData[ixy].diffusionRate);
	}
      }
      fclose(fpdebug);
#endif
#if USE_MPI
    }
#endif
  }

#if USE_MPI
  /* Share tfReturn with other ranks */
  if (tfReturn) {
    if (myid==0) {
      long id;
      for (id=1; id<n_processors; id++)
	MPI_Send(tfReturn, 1, tfReturnDataType, id, 1, MPI_COMM_WORLD);
    } else {
      MPI_Recv(tfReturn, 1, tfReturnDataType, 0, 1, MPI_COMM_WORLD, &mpiStatus);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
#endif

  /* Free arrays of data, but only if they are not being saved */
  if (xyTfData && xyTfData!=allXyTfData)
    free(xyTfData);
  xyTfData = NULL;
  if (deltaTfData && deltaTfData!=allDeltaTfData)
    free(deltaTfData);
  deltaTfData = NULL;
  
  return 1;
}

void outputTuneFootprint() 
{
  /* Output previously-saved tune footprint data */
  /* For MPI, only the master writes */
  SDDS_DATASET sddsOut;
  long id, ix, iy, iRow;
  
  if (!tuneFootprintOn)
    return;
  
  if (delta_output && allDeltaTfData) {
#if USE_MPI
    if (myid==0) {
#endif
      if (!SDDS_InitializeOutput(&sddsOut, SDDS_BINARY, 0, NULL, NULL, delta_output) ||
          SDDS_DefineColumn(&sddsOut, "delta", "$gd$r", NULL, NULL, NULL, SDDS_DOUBLE, 0)<0 ||
          SDDS_DefineColumn(&sddsOut, "nux", "$gn$r$bx$n", NULL, NULL, NULL, SDDS_DOUBLE, 0)<0 ||
          SDDS_DefineColumn(&sddsOut, "nuy", "$gn$r$by$n", NULL, NULL, NULL, SDDS_DOUBLE, 0)<0 ||
          SDDS_DefineColumn(&sddsOut, "diffusionRate", "log$b10$n(($gDn$r$bx$n$a2$n+$gDn$r$bx$n$a2$n)/Turns)", NULL, NULL, NULL, SDDS_DOUBLE, 0)<0 ||
          !SDDS_WriteLayout(&sddsOut) ||
          !SDDS_StartPage(&sddsOut, ndelta)) {
        SDDS_SetError("Problem setting up chromatic tune footprint output file");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      }
      for (id=iRow=0; id<ndelta; id++) {
        if (filtered_output && !allDeltaTfData[id].used)
          continue;
        if (!SDDS_SetRowValues(&sddsOut, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, iRow,
                               0, allDeltaTfData[id].delta,
                               1, allDeltaTfData[id].nu[0],
                               2, allDeltaTfData[id].nu[1],
                               3, allDeltaTfData[id].diffusionRate,
                               -1)) {
          SDDS_SetError("Problem setting up chromatic tune footprint output values");
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
        iRow++;
      }
      if (!SDDS_WriteTable(&sddsOut) || !SDDS_Terminate(&sddsOut)) {
          SDDS_SetError("Problem writing chromatic tune footprint output file");
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      }
#if USE_MPI
    }
#endif
  }
  if (xy_output && allXyTfData) {
#if USE_MPI
    if (myid==0) {
#endif
      if (!SDDS_InitializeOutput(&sddsOut, SDDS_BINARY, 0, NULL, NULL, xy_output) ||
          SDDS_DefineColumn(&sddsOut, "x", NULL, "m", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
          SDDS_DefineColumn(&sddsOut, "y", NULL, "m", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
          SDDS_DefineColumn(&sddsOut, "nux", "$gn$r$bx$n", NULL, NULL, NULL, SDDS_DOUBLE, 0)<0 ||
          SDDS_DefineColumn(&sddsOut, "nuy", "$gn$r$by$n", NULL, NULL, NULL, SDDS_DOUBLE, 0)<0 ||
          SDDS_DefineColumn(&sddsOut, "diffusionRate", "log$b10$n(($gDn$r$bx$n$a2$n+$gDn$r$bx$n$a2$n)/Turns)", NULL, NULL, NULL, SDDS_DOUBLE, 0)<0 ||
          !SDDS_WriteLayout(&sddsOut) ||
          !SDDS_StartPage(&sddsOut, nx*ny)) {
        SDDS_SetError("Problem setting up amplitude tune footprint output file");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      }
      for (ix=iRow=0; ix<nx; ix++) {
        for (iy=0; iy<ny; iy++) {
          id = ix + iy*nx;
          if (filtered_output && !allXyTfData[id].used)
            continue;
          if (!SDDS_SetRowValues(&sddsOut, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, iRow,
                                 0, allXyTfData[id].position[0],
                                 1, allXyTfData[id].position[1],
                                 2, allXyTfData[id].nu[0],
                                 3, allXyTfData[id].nu[1],
                                 4, allXyTfData[id].diffusionRate,
                               -1)) {
            SDDS_SetError("Problem setting up amplitude tune footprint output values");
            SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
          }
          iRow++;
        }
      }
      if (!SDDS_WriteTable(&sddsOut) || !SDDS_Terminate(&sddsOut)) {
          SDDS_SetError("Problem writing amplitude tune footprint output file");
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      }
#if USE_MPI
    }
#endif
  }

#if USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
}


void finishTuneFootprint()
{
  tuneFootprintOn = 0;
}

