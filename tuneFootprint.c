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

void determineDeltaTuneFootprint(DELTA_TF_DATA *deltaTfData, long nDelta, double *tuneRange, double *deltaRange)
{
  long id, id0, id1, id2, coord;
  double delta0, delta1, delta2, nu0[2];
  double nuMin, nuMax;

  /* Find dividing point between positive and negative delta */
  id0 = -1;
  for (id=0; id<nDelta; id++) {
    if (deltaTfData[id].nu[0]<=0 || deltaTfData[id].nu[0]>=1 ||
        deltaTfData[id].nu[1]<=0 || deltaTfData[id].nu[1]>=1 || deltaTfData[id].diffusionRate>diffusion_rate_limit)
      continue;
    if (deltaTfData[id].delta>=0) {
      id0 = id;
      nu0[0] = deltaTfData[id].nu[0];
      nu0[1] = deltaTfData[id].nu[1];
      delta0 = deltaTfData[id].delta;
      break;
    }
  }
  if (id0==-1) {
    tuneRange[0] = tuneRange[1] = deltaRange[0] = deltaRange[1]= 0;
    return;
  }

  id1 = id2 = -1;
  delta1 = delta2 = 0;
  for (coord=0; coord<2; coord++) {
    for (id=id0; id>=0; id--) {
      /* scan to negative delta side from center until we find a place where the 
       * tune is invalid or crosses a half integer */
      if (deltaTfData[id].nu[coord]<=0 || deltaTfData[id].nu[coord]>=1 || 
          ((long)(2*deltaTfData[id].nu[coord])!=(long)(2*nu0[coord])) || deltaTfData[id].diffusionRate>diffusion_rate_limit)
        break;
      delta1 = deltaTfData[id].delta;
      id1 = id;
    }  
    for (id=id0; id<nDelta; id++) {
      /* scan to positive delta side from center until we find a place where the 
       * tune is invalid or crosses a half integer */
      if (deltaTfData[id].nu[coord]<=0 || deltaTfData[id].nu[coord]>=1 || 
          ((long)(2*deltaTfData[id].nu[coord])!=(long)(2*nu0[coord])))
        break;
      delta2 = deltaTfData[id].delta;
      id2 = id;
    }
    deltaRange[coord] = MIN(fabs(delta1), fabs(delta2));

    nuMin = 1;
    nuMax = 0;
    for (id=id1; id<=id2; id++) {
      if (deltaTfData[id].nu[coord]<nuMin)
        nuMin = deltaTfData[id].nu[coord];
      if (deltaTfData[id].nu[coord]>nuMax)
        nuMax = deltaTfData[id].nu[coord];
    }
    tuneRange[coord] = nuMax - nuMin;
  }
  
}

void determineXyTuneFootprint(XY_TF_DATA *xyTfData, long nx, long ny, double *tuneRange, double *positionRange)
{
  long id;
  double nuMin, nuMax, pMin, pMax;
  double distance, bestDistance;
  long ix0, ix1, ix2, ix, iy, iy1, ic;
  
  /* check for indexing issues */
  for (ix=0; ix<nx; ix++) {
    for (iy=0; iy<ny; iy++) {
      id = ix + nx*iy;
      if (ix!=xyTfData[id].ix || xyTfData[id].iy!=iy) {
        fprintf(stderr, "indexing problem in determineXyTuneFootprint: ix=%ld, iy=%ld, id=%ld, ix[id]=%d, iy[id]=%d\n",
                ix, iy, id, xyTfData[id].ix, xyTfData[id].iy);
        bombElegant("Indexing bug", NULL);
      }
    }
  }

  /* Find point closest to origin (iy=0) */
  bestDistance = DBL_MAX;
  ix0 = -1;
  for (ix=0; ix<nx; ix++) {
    if (xyTfData[ix].nu[0]<=0 || xyTfData[ix].nu[0]>=1 ||
        xyTfData[ix].nu[1]<=0 || xyTfData[ix].nu[1]>=1 || xyTfData[ix].diffusionRate>diffusion_rate_limit)
      continue;
    distance = fabs(xyTfData[ix].position[0]);
    if (distance<bestDistance) {
      ix0 = ix;
      bestDistance = distance;
    }
  }
  if (ix0==-1) {
    tuneRange[0] = tuneRange[1] = positionRange[0] = positionRange[1]= 0;
    return;
  }

  ix1 = ix2 = -1;
  for (iy=0; iy<ny; iy++) {
    /* find limiting indices ix1 and ix2 of valid data for iy */
    for (ix=ix0; ix>=0; ix--) {
      id = ix + nx*iy;
      if (xyTfData[id].nu[0]<=0 || xyTfData[id].nu[0]>=1 ||
          xyTfData[id].nu[1]<=0 || xyTfData[id].nu[1]>=1 || xyTfData[id].diffusionRate>diffusion_rate_limit)
        break;
      ix1 = ix;
    }
    for (ix=ix0; ix<nx; ix++) {
      id = ix + nx*iy;
      if (xyTfData[id].nu[0]<=0 || xyTfData[id].nu[0]>=1 ||
          xyTfData[id].nu[1]<=0 || xyTfData[id].nu[1]>=1 || xyTfData[id].diffusionRate>diffusion_rate_limit)
        break;
      ix2 = ix;
    }
    /* Mark all points at or above this iy value with ix<ix1 or ix>ix2 as bad */
    for (iy1=iy; iy1<ny; iy1++) {
      for (ix=0; ix<nx; ix++) {
        id = ix + nx*iy1;
        if (ix<ix1 || ix>ix2)
          xyTfData[id].diffusionRate = diffusion_rate_limit+1;
      }
    }
  }

  /* find the spread in tune and position for points with acceptable characteristics */
  for (ic=0; ic<2; ic++) {
    nuMin = 1;
    nuMax = 0;
    pMax = -(pMin = DBL_MAX);
    for (ix=0; ix<nx; ix++) {
      for (iy=0; iy<ny; iy++) {
        id = ix + nx*iy;
        if (xyTfData[id].nu[0]<=0 || xyTfData[id].nu[0]>=1 ||
            xyTfData[id].nu[1]<=0 || xyTfData[id].nu[1]>=1 || xyTfData[id].diffusionRate>diffusion_rate_limit)
          continue;
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
    tuneRange[ic] = nuMax - nuMin;
    positionRange[ic] = pMax - pMin;
  }
}

#if USE_MPI
MPI_Datatype xyTfDataType, deltaTfDataType;
void setupTuneFootprintDataTypes ()
{
  MPI_Datatype oldType[2];
  int blockLength[2];
  MPI_Aint offset[2];
  XY_TF_DATA xyTfExample;
  DELTA_TF_DATA deltaTfExample;
  long i;
  
  oldType[0] = MPI_SHORT;
  oldType[1] = MPI_DOUBLE;

  blockLength[0] = 3;
  MPI_Get_address(&xyTfExample.ix, &offset[0]);
  MPI_Get_address(&xyTfExample.diffusionRate, &offset[1]);
  for (i=1; i>=0; i--)
    offset[i] -= offset[0];
  blockLength[1] = 5;
  MPI_Type_create_struct(2, blockLength, offset, oldType, &xyTfDataType);
  MPI_Type_commit(&xyTfDataType);
  
  blockLength[0] = 2;
  MPI_Get_address(&deltaTfExample.idelta, &offset[0]);
  MPI_Get_address(&deltaTfExample.diffusionRate, &offset[1]);
  for (i=1; i>=0; i--)
    offset[i] -= offset[0];
  blockLength[1] = 4;
  MPI_Type_struct(2, blockLength, offset, oldType, &deltaTfDataType);
  MPI_Type_commit(&deltaTfDataType);
}
#endif

void setupTuneFootprint(
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
  if (delta_min>delta_max)
    bombElegant("delta_min > delta_max", NULL);
  if (quadratic_spacing) {
    if (xmin<0)
      xmin = 0;
    if (ymin<0)
      ymin = 0;
  }
  if (nx<1)
    nx = 1;
  if (ny<1)
    ny = 1;
  if (ndelta<1)
    ndelta = 1;

}


long doTuneFootprint(
                    RUN *run,
                    VARY *control,
                    double *referenceCoord,
                    ERRORVAL *errcon,
                    LINE_LIST *beamline,
                    TUNE_FOOTPRINTS *tfOutput
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
  XY_TF_DATA *xyTfData, *allXyTfData;
  DELTA_TF_DATA *deltaTfData, *allDeltaTfData;
  long my_nxy, my_ndelta;
  double chromTuneRange[2], chromDeltaRange[2], xyTuneRange[2], xyPositionRange[2];
#define DEBUG 1
#ifdef DEBUG
  FILE *fpdebug = NULL;
#endif
#if USE_MPI  
  MPI_Status mpiStatus;
  setupTuneFootprintDataTypes();
  /* Note that the master is a working processor for this algorithm */
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
  if (nx*ny<n_processors && myid==0) {
    printf("NB: number of x, y grid points less than the number of cores, which is a waste of resources\n");
  }
  my_nxy = (nx*ny)/n_processors+1;
  my_ndelta = ndelta/n_processors+1;
#else
  my_nxy = nx*ny;
  my_ndelta = ndelta;
#endif
  
  xyTfData = calloc(my_nxy, sizeof(*xyTfData));
  deltaTfData = calloc(my_ndelta, sizeof(*deltaTfData));

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

  /* First perform delta scan */
  x = x_for_delta;
  y = y_for_delta;
  if (ndelta>1)
    ddelta = (delta_max-delta_min)/(ndelta-1);
  else
    ddelta = 0;
  for (idelta=my_idelta=0; idelta<ndelta; idelta++) {
    delta = delta_min + idelta*ddelta;
    memcpy(startingCoord, referenceCoord, sizeof(*startingCoord)*6);
#if USE_MPI
    if (myid == idelta%n_processors) /* Partition the job according to particle ID */
#endif
	  {
            if (verbosity>=1) {
#if USE_MPI
              if (myid==0)
#endif
                printf("computing tune for delta = %le\n", delta);
            }
            lost = 0;
	    if (!computeTunesFromTracking(firstTune, firstAmplitude,
					  beamline->matrix, beamline, run,
					  startingCoord, x, y, delta, turns,
					  0, endingCoord, NULL, NULL, 1, 1) ||
		firstTune[0]>1.0 || firstTune[0]<0 || firstTune[1]>1.0 || firstTune[1]<0) {
              lost = 1;
	    } else {
              memcpy(startingCoord, endingCoord, sizeof(*startingCoord)*6);
              if (!computeTunesFromTracking(secondTune, secondAmplitude,
                                            beamline->matrix, beamline, run,
                                            startingCoord, 0.0, 0.0, 0.0, turns,
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
            if (lost)
              deltaTfData[my_idelta].nu[0] = deltaTfData[my_idelta].nu[1] = 
                deltaTfData[my_idelta].diffusionRate = -1;
            else {
              deltaTfData[my_idelta].nu[0] = firstTune[0];
              deltaTfData[my_idelta].nu[1] = firstTune[1];
              deltaTfData[my_idelta].diffusionRate = log10((sqr(secondTune[0] - firstTune[0]) + sqr(secondTune[1] - firstTune[1]))/turns);
            }
            my_idelta ++;
	    if (verbosity) {
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
    allDeltaTfData  = calloc(my_ndelta*n_processors, sizeof(*allDeltaTfData));
    memcpy(allDeltaTfData, deltaTfData, sizeof(*deltaTfData)*my_ndelta);
    /* receive data */
    iTotal = my_ndelta;
    for (id=1; id<n_processors; id++) {
      if (MPI_Recv(allDeltaTfData+iTotal, my_ndelta, deltaTfDataType, id, 1, MPI_COMM_WORLD, &mpiStatus)!=MPI_SUCCESS) {
        printf("error receiving delta data from processor %ld\n", id);
        exit(1);
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
        my_ndelta = idelta+1;
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

    determineDeltaTuneFootprint(allDeltaTfData, my_ndelta, chromTuneRange, chromDeltaRange);
    printf("nux/chromatic: tune range=%le  delta range = %le\nnuy/chromatic: tune range=%le  delta range = %le\n",
           chromTuneRange[0], chromDeltaRange[0],
           chromTuneRange[1], chromDeltaRange[1]);
    if (tfOutput) {
      memcpy(tfOutput->chromaticTuneRange, chromTuneRange, sizeof(*chromTuneRange)*2);
      memcpy(tfOutput->deltaRange, chromDeltaRange, sizeof(*chromDeltaRange)*2);
      tfOutput->deltaRange[2] = MIN(tfOutput->deltaRange[0], tfOutput->deltaRange[1]);
    }
#if USE_MPI
  }
#endif

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
      x = xmin + (xmax-xmin)*sqrt((ix+1.)/nx);
    } else {
      x = xmin + ix*dx;
    }
    for (iy=0; iy<ny; iy++) {
      if (quadratic_spacing) {
        y = ymin + (ymax-ymin)*sqrt((iy+1.)/ny);
      } else {
        y = ymin + iy*dy;
      }
      memcpy(startingCoord, referenceCoord, sizeof(*startingCoord)*6);
#if USE_MPI
      if (myid == (ix*ny+iy)%n_processors) /* Partition the job according to particle ID */
#endif
	  {
            lost = 0;
	    if (!computeTunesFromTracking(firstTune, firstAmplitude,
					  beamline->matrix, beamline, run,
					  startingCoord, x, y, delta, turns,
					  0, endingCoord, NULL, NULL, 1, 1) ||
		firstTune[0]>1.0 || firstTune[0]<0 || firstTune[1]>1.0 || firstTune[1]<0) {
              lost = 1;
            } else {
              memcpy(startingCoord, endingCoord, sizeof(*startingCoord)*6);
              if (!computeTunesFromTracking(secondTune, secondAmplitude,
                                            beamline->matrix, beamline, run,
                                            startingCoord, 0.0, 0.0, 0.0, turns,
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
              xyTfData[my_ixy].nu[0] = xyTfData[my_ixy].nu[1] = 
                xyTfData[my_ixy].diffusionRate = -1;
            } else {
              xyTfData[my_ixy].nu[0] = firstTune[0];
              xyTfData[my_ixy].nu[1] = firstTune[1];
              xyTfData[my_ixy].diffusionRate = log10((sqr(secondTune[0] - firstTune[0]) + sqr(secondTune[1] - firstTune[1]))/turns);
            }
            my_ixy ++;

	    if (verbosity) {
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
    long iTotal;
    allXyTfData  = calloc(my_nxy*n_processors, sizeof(*allXyTfData));
    memcpy(allXyTfData, xyTfData, sizeof(*xyTfData)*my_nxy);
    /* receive data */
    iTotal = my_nxy;
    for (idelta=1; idelta<n_processors; idelta++) {
      if (MPI_Recv(allXyTfData+iTotal, my_nxy, xyTfDataType, idelta, 1, MPI_COMM_WORLD, &mpiStatus)!=MPI_SUCCESS) {
        printf("error receiving xy data from processor %ld\n", idelta);
        exit(1);
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

    determineXyTuneFootprint(allXyTfData, nx, ny, xyTuneRange, xyPositionRange);
    printf("nux/amplitude: tune range=%le  position range = %le\nnuy/amplitude: tune range=%le position range = %le\n",
           xyTuneRange[0], xyPositionRange[0],
           xyTuneRange[1], xyPositionRange[1]);
    if (tfOutput) {
      memcpy(tfOutput->amplitudeTuneRange, xyTuneRange, sizeof(*xyTuneRange)*2);
      memcpy(tfOutput->positionRange, xyPositionRange, sizeof(*xyPositionRange)*2);
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

  free(xyTfData);
  free(deltaTfData);

  return(1);
}

void finishTuneFootprint()
{
  
}

