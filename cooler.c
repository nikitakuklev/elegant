/*************************************************************************\
 * Copyright (c) 2003 The University of Chicago, as Operator of Argonne
 * National Laboratory.
 * Copyright (c) 2003 The Regents of the University of California, as
 * Operator of Los Alamos National Laboratory.
 * This file is distributed subject to a Software License Agreement found
 * in the file LICENSE that is included with this distribution. 
 \*************************************************************************/

/*  file: cooler.cc 
 *  based on: tfeedback.c
 *
 */

#if defined(SOLARIS) && !defined(__GNUC__)
#  include <sunmath.h>
#endif

#include "mdb.h"
#include "track.h"
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>
//#include <algorithm>

#define cms c_mks
#define twopi PIx2

void coolerPickup(CPICKUP *cpickup, double **part0, long np0, long pass, double Po, long idSlotsPerBunch) {

  /*Initialize some variables */
  double sum;
  long i;
  double *time0 = NULL;    /* array to record arrival time of each particle */
  long *ibParticle = NULL; /* array to record which bucket each particle is in */
  long **ipBucket = NULL;  /* array to record particle indices in part0 array for all particles in each bucket */
  long *npBucket = NULL;   /* array to record how many particles are in each bucket */
  long nBuckets = 0;
  long ib;

  // MPI Required vars
#if USE_MPI
  long npTotal;
  double sumTotal;

  if (!partOnMaster) {
    // this element does nothing in single particle mode (e.g., trajectory, orbit, ..)
    if (notSinglePart == 0)
      return;
  }
#endif

#ifdef DEBUG
  printf("Continuing cooler pickup\n");
  fflush(stdout);
#endif

  // makes sure element is initialized
  if (cpickup->initialized == 0)
    initializeCoolerPickup(cpickup);
  else {
    for (ib = 0; ib < cpickup->nBunches; ib++) {
      free(cpickup->long_coords[ib]);
      free(cpickup->horz_coords[ib]);
      free(cpickup->vert_coords[ib]);
      free(cpickup->pid[ib]);
      hdestroy(cpickup->pidHashTable[ib]);
      free(cpickup->index[ib]);
    }
    free(cpickup->long_coords);
    free(cpickup->horz_coords);
    free(cpickup->vert_coords);
    free(cpickup->pid);
    free(cpickup->npBunch);
    free(cpickup->pidHashTable);
    free(cpickup->index);

    cpickup->long_coords = cpickup->horz_coords = cpickup->vert_coords = NULL;
    cpickup->pid = NULL;
    cpickup->npBunch = NULL;
    cpickup->pidHashTable = NULL;
    cpickup->index = NULL;
  }

  // runs only between user specified passes
  if ((cpickup->startPass > 0 && pass < cpickup->startPass) || (cpickup->endPass > 0 && pass > cpickup->endPass))
    return;
  if (cpickup->updateInterval > 1 && pass % cpickup->updateInterval != 0)
    return;
  cpickup->lastPass = pass;

#ifdef DEBUG
  printf("Initialized cooler pickup\n");
  fflush(stdout);
#endif

  // assign buckets normal
  if (isSlave || !notSinglePart
#if USE_MPI
      || (partOnMaster && myid == 0)
#endif
  )
    index_bunch_assignments(part0, np0, cpickup->bunchedBeamMode ? idSlotsPerBunch : 0, Po,
                            &time0, &ibParticle, &ipBucket, &npBucket, &nBuckets, -1);

#ifdef DEBUG
  printf("indexed bunch assignments\n");
  fflush(stdout);
#endif

#if USE_MPI
  if (partOnMaster || myid != 0) {
#endif
#ifdef DEBUG
    printf("%ld buckets found\n", nBuckets);
    fflush(stdout);
#endif

    if (cpickup->nBunches == 0 || cpickup->nBunches != nBuckets) {
      if (cpickup->nBunches != nBuckets) {
        printf("Number of bunches has changed, re-initializing cooler pickup\n");
        fflush(stdout);
      }
      cpickup->nBunches = nBuckets;
      cpickup->tAverage = (double *)SDDS_Realloc(cpickup->tAverage, sizeof(*cpickup->tAverage) * nBuckets);
      for (i = 0; i < nBuckets; i++) {
        cpickup->tAverage[i] = 0;
      }
    }

    // Allocate space for particle coordinates, which may be needed by the CKICKER
    cpickup->long_coords = (double **)malloc(sizeof(double *) * nBuckets);
    cpickup->horz_coords = (double **)malloc(sizeof(double *) * nBuckets);
    cpickup->vert_coords = (double **)malloc(sizeof(double *) * nBuckets);
    cpickup->pid = (long **)malloc(sizeof(long *) * nBuckets);
    cpickup->npBunch = (long *)malloc(sizeof(long) * nBuckets);
    cpickup->pidHashTable = (htab **)malloc(sizeof(htab *) * nBuckets);
    cpickup->index = (long **)malloc(sizeof(long *) * nBuckets);

    // Important stuff happens here
    //    tAverage takes average time of arrival
    //    part_coords takes each particles individual time of arival

    for (ib = 0; ib < nBuckets; ib++) {
      cpickup->long_coords[ib] = (double *)malloc(sizeof(double) * npBucket[ib]);
      cpickup->horz_coords[ib] = (double *)malloc(sizeof(double) * npBucket[ib]);
      cpickup->vert_coords[ib] = (double *)malloc(sizeof(double) * npBucket[ib]);
      cpickup->pid[ib] = (long *)malloc(sizeof(long) * npBucket[ib]);
      cpickup->npBunch[ib] = npBucket[ib];
      cpickup->pidHashTable[ib] = hcreate(12);
      cpickup->index[ib] = (long *)malloc(sizeof(long) * npBucket[ib]);

      sum = 0; // Used to average all particle arival times in the pickup
      for (i = 0; i < npBucket[ib]; i++) {
        char pidString[32];
        sum += time0[ipBucket[ib][i]];
        // store the particle time in the coord array to be used in kicker
        cpickup->long_coords[ib][i] = time0[ipBucket[ib][i]];
        // store x and y for use in computing transverse effects
        cpickup->horz_coords[ib][i] = part0[ipBucket[ib][i]][0] + cpickup->dx;
        cpickup->vert_coords[ib][i] = part0[ipBucket[ib][i]][2] + cpickup->dy;
        // will use this to double-check ordering of particles between pickup and kicker
        cpickup->pid[ib][i] = part0[ipBucket[ib][i]][particleIDIndex];
        snprintf(pidString, 32, "%ld", cpickup->pid[ib][i]);
        cpickup->index[ib][i] = i;
        if (!hadd(cpickup->pidHashTable[ib], pidString, strlen(pidString), &(cpickup->index[ib][i]))) {
          bombElegantVA("Problem creating PID hash table: duplicate PID %ld\n", cpickup->pid[ib][i]);
        }
      }

#ifdef DEBUG
      printf("Finished looping over %ld particles, sum = %le\n", npBucket[ib], sum);
      fflush(stdout);
#endif

#if USE_MPI
      if (!partOnMaster) {
#  ifdef DEBUG
        printf("Reducing sum data\n");
        fflush(stdout);
#  endif
        MPI_Allreduce(&npBucket[ib], &npTotal, 1, MPI_LONG, MPI_SUM, workers);
#  ifdef DEBUG
        printf("npTotal = %ld\n", npTotal);
        fflush(stdout);
#  endif
        MPI_Allreduce(&sum, &sumTotal, 1, MPI_DOUBLE, MPI_SUM, workers);
#  ifdef DEBUG
        printf("sumTotal = %le\n", sumTotal);
        fflush(stdout);
#  endif
        cpickup->tAverage[ib] = sumTotal / npTotal;
      } else
        cpickup->tAverage[ib] = sum / npBucket[ib];
#else
    cpickup->tAverage[ib] = sum / npBucket[ib];
#endif
    }

#if USE_MPI
  } /* if (!partOnMaster || myid==0) */
#endif

  //memory managment
  if (isSlave || !notSinglePart
#if USE_MPI
      || (partOnMaster && myid == 0)
#endif
  )
    free_bunch_index_memory(time0, ibParticle, ipBucket, npBucket, nBuckets);

#if USE_MPI
#  ifdef DEBUG
  printf("Waiting on barrier before returning from coolerPickup\n");
  fflush(stdout);
#  endif
  if (!partOnMaster)
    MPI_Barrier(MPI_COMM_WORLD);
#endif

#ifdef DEBUG
  printf("Returning from coolerPickup\n");
  fflush(stdout);
#endif
}

// Pickup initialization function{
void initializeCoolerPickup(CPICKUP *cpickup) {
  if (cpickup->ID == NULL || !strlen(cpickup->ID))
    bombElegant("you must give an ID string for CPICKUP", NULL);

  cpickup->nBunches = 0;

  if (cpickup->updateInterval < 1)
    cpickup->updateInterval = 1;

  cpickup->initialized = 1;
  cpickup->lastPass = -1;
  cpickup->long_coords = cpickup->horz_coords = cpickup->vert_coords = NULL;
  cpickup->pid = NULL;
  cpickup->npBunch = NULL;
}

struct E_params {
  double x;
  double y;
  double gamma;
  double lambda;
};

double Ex(double theta, void *params) {
  //unpacking parameters
  struct E_params *p = (struct E_params *)params;
  double x = (p->x);
  double y = (p->y);
  double gamma = (p->gamma);
  double lambda = (p->lambda);

  // define some parameters for calc
  double k0 = twopi / lambda;
  double rho = sqrt(x * x + y * y);
  double phi = atan2(y, x);

  double J0_term = gsl_sf_bessel_Jn(0, rho * k0 * theta / (1 + pow(gamma * theta, 2)));
  double J2_term = gsl_sf_bessel_Jn(2, rho * k0 * theta / (1 + pow(gamma * theta, 2))) * pow(gamma * theta, 2) * cos(2 * phi);
  double I0_term = theta / pow(1 + pow(gamma * theta, 2), 4);

  return (J0_term + J2_term) * I0_term;
}

typedef struct {
  int index;
  double t;
  double x;
  double y;
} indexed_coord;

int compare_indexed_coord(const void *arg1, const void *arg2) {
  indexed_coord const *lhs = (indexed_coord const *)(arg1);
  indexed_coord const *rhs = (indexed_coord const *)(arg2);
  return (lhs->t < rhs->t) ? -1 : ((rhs->t < lhs->t) ? 1 : 0);
}

void coolerKicker(CKICKER *ckicker, double **part0, long np0, LINE_LIST *beamline, long pass, long nPasses, char *rootname, double Po, long idSlotsPerBunch) {
  long i, j;
  double *time0 = NULL;    /* array to record arrival time of each particle */
  long *ibParticle = NULL; /* array to record which bucket each particle is in */
  long **ipBucket = NULL;  /* array to record particle indices in part0 array for all particles in each bucket */
  long *npBucket = NULL;   /* array to record how many particles are in each bucket */
  long nBuckets = 0;
  long updateInterval;
  double Ex0;
  long ib;

#if USE_MPI
  double sumTotal, s_sumTotal;
  long npTotal;

  if (!partOnMaster) {
    if (notSinglePart == 0)
      /* this element does nothing in single particle mode (e.g., trajectory, orbit, ..) */
      return;
  }
#endif

  if ((ckicker->startPass > 0 && pass < ckicker->startPass) || (ckicker->endPass > 0 && pass > ckicker->endPass))
    return;

#ifdef DEBUG
  printf("Continuing cooler kicker\n");
  fflush(stdout);
#endif

  if (ckicker->transverseMode != 0) {
    // integration variables
    double error;
    size_t neval;
    gsl_function F;

    struct E_params params = {0, 0, Po, ckicker->lambda_rad}; // {x, y, gamma, lambda}

    F.function = &Ex;
    F.params = &params;

    // (gsl_function F, a, b, epsabs, epsrel, result, error, neval)
    gsl_integration_qng(&F, 0, ckicker->angle_rad / Po, 0, 1e-7, &Ex0, &error, &neval);
  }

  if (ckicker->initialized == 0)
    initializeCoolerKicker(ckicker, beamline, nPasses * nBuckets, rootname);

#ifdef DEBUG
  printf("Initialized cooler kicker\n");
  fflush(stdout);
#endif

  if (ckicker->startPass > 0 && ckicker->startPass != ckicker->pickup->startPass)
    bombElegantVA((char *)"CKICKER linked to CPICKUP with different START_PASS value (%ld vs %ld).",
                  ckicker->startPass, ckicker->pickup->startPass);
  if (ckicker->endPass > 0 && ckicker->endPass != ckicker->pickup->endPass)
    bombElegantVA((char *)"CKICKER linked to CPICKUP with different END_PASS value (%ld vs %ld).",
                  ckicker->endPass, ckicker->pickup->endPass);

  if ((updateInterval = ckicker->pickup->updateInterval * ckicker->updateInterval) <= 0)
    bombElegantVA((char *)"CKICKER and CPICKUP with ID=%s have UPDATE_INTERVAL product of %d", ckicker->ID, updateInterval);
  if (pass % updateInterval != 0)
    return;

  if (ckicker->pickup->lastPass != pass)
    bombElegantVA("CKICKER and CPICKUP with ID=%s are not synchronized to the same pass (%ld vs %ld)\n",
                  pass, ckicker->pickup->lastPass);

  if (isSlave || !notSinglePart
#if USE_MPI
      || (partOnMaster && myid == 0)
#endif
  )
    index_bunch_assignments(part0, np0, ckicker->bunchedBeamMode ? idSlotsPerBunch : 0, Po,
                            &time0, &ibParticle, &ipBucket, &npBucket, &nBuckets, -1);

#ifdef DEBUG
  printf("indexed bunch assignments\n");
  fflush(stdout);
#endif

#if USE_MPI
  if (partOnMaster || myid != 0) {
#endif

    if (ckicker->nBunches == 0 || ckicker->nBunches != nBuckets) {
      if (ckicker->nBunches != nBuckets) {
        printf("Number of bunches has changed, re-initializing cooler driver.\n");
        fflush(stdout);
      }
      ckicker->nBunches = nBuckets;
    }

    if (ckicker->nBunches != ckicker->pickup->nBunches)
      bombElegantVA("mismatch in number of buckets between CKICKER (%ld) and CPICKUP (%ld)",
                    ckicker->nBunches, ckicker->pickup->nBunches);

    for (ib = 0; ib < nBuckets; ib++) {
      // ================================= //
      // This is where the kick is applied //
      // ================================= //
      double sum, nom_kick, env_term;
      double t_avg, t_i, x_i, y_i;
      char pidString[32];
      long *sitmp, *sourceIndex = NULL;
      indexed_coord *incoherent_ar = NULL;
      double s_sum;
      // double s_avg;

      sourceIndex = (long *)malloc(sizeof(long) * npBucket[ib]);
      for (i = 0; i < npBucket[ib]; i++) {
        snprintf(pidString, 32, "%ld", (long)part0[ipBucket[ib][i]][particleIDIndex]);
        if (!hfind(ckicker->pickup->pidHashTable[ib], pidString, strlen(pidString)))
          bombElegantVA("PID of particle in CKICKER not present at CPICKUP! Cooling ID=%s\n", ckicker->ID);
        if (!(sitmp = (long *)hstuff(ckicker->pickup->pidHashTable[ib])))
          bombElegantVA("Problem retrieving index of PID %ld. Cooling ID=%s\n", (long)part0[ipBucket[ib][i]][particleIDIndex],
                        ckicker->ID);
        sourceIndex[i] = *sitmp;
      }

      if (ckicker->incoherentMode != 0)
        incoherent_ar = (indexed_coord *)malloc(sizeof(indexed_coord) * npBucket[ib]);
      sum = 0;
      s_sum = 0;
      for (i = 0; i < npBucket[ib]; i++) {
        sum += time0[ipBucket[ib][i]];
        s_sum += part0[ipBucket[ib][i]][4];

        if (ckicker->incoherentMode != 0) {
          incoherent_ar[i].index = i;
          incoherent_ar[i].x = ckicker->pickup->horz_coords[ib][sourceIndex[i]];
          incoherent_ar[i].y = ckicker->pickup->vert_coords[ib][sourceIndex[i]];
          incoherent_ar[i].t = ckicker->pickup->long_coords[ib][sourceIndex[i]];
        }
      }

      if (ckicker->incoherentMode != 0) {
        // sort into arrival-time order, earliest particles first
        qsort(incoherent_ar, npBucket[ib], sizeof(indexed_coord), compare_indexed_coord);
      }

      // t_avg equals average arrival time at kicker - average time at pickup
#if USE_MPI
      if (!partOnMaster) {
        MPI_Allreduce(&sum, &sumTotal, 1, MPI_DOUBLE, MPI_SUM, workers);
        MPI_Allreduce(&npBucket[ib], &npTotal, 1, MPI_LONG, MPI_SUM, workers);
        t_avg = sumTotal / npTotal - ckicker->pickup->tAverage[ib];

        MPI_Allreduce(&s_sum, &s_sumTotal, 1, MPI_DOUBLE, MPI_SUM, workers);
        // s_avg = s_sumTotal / npTotal;

      } else {
        t_avg = sum / npBucket[ib] - ckicker->pickup->tAverage[ib];
        // s_avg = s_sum / npBucket[ib];
      }
#else
    t_avg = sum / npBucket[ib] - ckicker->pickup->tAverage[ib];
    // s_avg = s_sum / npBucket[ib];
#endif

#ifdef DEBUG
      printf("bunch %ld has %ld particles, t_avg = %le\n", ib, npBucket[ib], t_avg);
      fflush(stdout);
#endif

      if (ckicker->incoherentMode != 0) {
        for (i = 0; i < npBucket[ib]; i++) {
          long index;
          index = incoherent_ar[i].index; /* index of target particle in the unsorted arrays */
          for (j = 0; j < npBucket[ib]; j++) {
            if (i != j) {
              if (fabs(incoherent_ar[j].t - incoherent_ar[i].t) < (ckicker->Nu * ckicker->lambda_rad) / c_mks) {

                double incoherent_phase, coherent_phase, incoherent_strength, kick;

                t_i = time0[ipBucket[ib][incoherent_ar[i].index]] - ckicker->pickup->long_coords[ib][sourceIndex[incoherent_ar[i].index]];

                coherent_phase = (t_avg - t_i) * c_mks * (twopi / ckicker->lambda_rad);
                incoherent_phase = (incoherent_ar[i].t - incoherent_ar[j].t) * c_mks * (twopi / ckicker->lambda_rad);

                double phi = (coherent_phase + incoherent_phase);
                double nu_2 = ckicker->Nu / 2;

                incoherent_strength = (1 - abs(phi / twopi) / (ckicker->Nu)) * (2.0 / twopi) * (100 * atan((phi + nu_2)) + atan(-100 * (phi - nu_2)));

                kick = -(ckicker->strength * incoherent_strength) * sin(phi + ckicker->phase * twopi);

                part0[ipBucket[ib][index]][5] += kick;
              }
            }
          }
        }
        free(incoherent_ar);
      }

      //printf("%li\n",pass);
      for (i = 0; i < npBucket[ib]; i++) {
        // particle coords
        t_i = time0[ipBucket[ib][i]] - ckicker->pickup->long_coords[ib][sourceIndex[i]];
        x_i = part0[ipBucket[ib][i]][0] - ckicker->magnification * ckicker->pickup->horz_coords[ib][sourceIndex[i]];
        y_i = part0[ipBucket[ib][i]][2] - ckicker->magnification * ckicker->pickup->vert_coords[ib][sourceIndex[i]];

        // on-axis kick = sin(c*(difference in arrival time)*2*Pi/lambda + phase)

        double nu_2 = ckicker->Nu / 2;
        double phi = (t_avg - t_i) * (c_mks / ckicker->lambda_rad);

        env_term = (2.0 / twopi) * (atan(100 * (phi + (2 * nu_2))) + atan(-100 * (phi - (2 * nu_2)))) * (1 - abs(phi) / (2 * nu_2));

        nom_kick = -(ckicker->strength) * env_term * sin((t_avg - t_i) * c_mks * (twopi / (ckicker->lambda_rad)) + ckicker->phase * twopi);

        if (ckicker->transverseMode != 0) {
          double Exi, error;
          size_t neval;
          gsl_function F;
          struct E_params params = {x_i, y_i, Po, ckicker->lambda_rad}; // {x, y, gamma, lambda}
          //struct E_params params = {0.0 , 0.0, Po, ckicker->lambda_rad};
          F.function = &Ex;
          F.params = &params;

          gsl_integration_qng(&F, 0, ckicker->angle_rad / Po, 0, 1e-7, &Exi, &error, &neval);

          part0[ipBucket[ib][i]][5] += nom_kick * Exi / Ex0;
        } else {
          part0[ipBucket[ib][i]][5] += nom_kick;
        }
      }
      free(sourceIndex);
      sourceIndex = NULL;
    }

    if (isSlave || !notSinglePart)
      free_bunch_index_memory(time0, ibParticle, ipBucket, npBucket, nBuckets);

#if USE_MPI
  } /* if (partOnMaster || myid!=0) */
  if (!partOnMaster)
    MPI_Barrier(MPI_COMM_WORLD);
#endif
}

// Kicker initialization function

void initializeCoolerKicker(CKICKER *ckicker, LINE_LIST *beamline, long nPasses, char *rootname) {
  ELEMENT_LIST *eptr;
  long pickupFound = 0;

  if (ckicker->ID == NULL || !strlen(ckicker->ID))
    bombElegant("you must give an ID string for CKICKER", NULL);

  eptr = beamline->elem;
  while (eptr) {
    if (eptr->type == T_CPICKUP && strcmp(ckicker->ID, ((CPICKUP *)eptr->p_elem)->ID) == 0) {
      pickupFound = 1;
      ckicker->pickup = ((CPICKUP *)eptr->p_elem);
      break;
    }
    eptr = eptr->succ;
  }
  if (!pickupFound)
    bombElegant("pickup not found for CKICKER", NULL);

  ckicker->nBunches = 0;

  if (ckicker->updateInterval < 1)
    ckicker->updateInterval = 1;
  ckicker->initialized = 1;
}

void setCoolingTrackingMode(ELEMENT_LIST *eptr) {
  short incoherentCooling = 0;
  if ((entity_description[T_CKICKER].flags & UNIPROCESSOR) &&
      (entity_description[T_CPICKUP].flags & UNIPROCESSOR))
    return;
  while (eptr && !incoherentCooling) {
    if (eptr->type == T_CKICKER)
      incoherentCooling = ((CKICKER *)(eptr->p_elem))->incoherentMode;
    eptr = eptr->succ;
  }
  if (incoherentCooling) {
    /* Force cooling into serial model */
    printWarning("CPICKUP and CKICKER elements running in serial mode", "Parallel cooling incompatible with inclusion of incoherent effects.");
    entity_description[T_CKICKER].flags &= ~MPALGORITHM;
    entity_description[T_CKICKER].flags |= UNIPROCESSOR;
    entity_description[T_CPICKUP].flags &= ~MPALGORITHM;
    entity_description[T_CPICKUP].flags |= UNIPROCESSOR;
  }
}
