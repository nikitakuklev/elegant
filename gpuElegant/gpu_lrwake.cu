#include <gpu_track.h>
#include <gpu_base.h>
#include <gpu_particle_template_function.hcu>
#include <gpu_reductions.h>
#include <gpu_trwake.h>                // gpu_computeTimeCoordinates
#include <gpu_bin_time_distribution.h> // gpu_binTimeDistribution_and_countBinned

class gpu_index_bunch_assignments_kernel1
{
 public:
  unsigned int idSlotsPerBunch;
  double *d_ib;

 gpu_index_bunch_assignments_kernel1(unsigned int idSlotsPerBunch,
                                          double *d_ib) : idSlotsPerBunch(idSlotsPerBunch), d_ib(d_ib){};

  __device__ void operator()(gpuParticleAccessor &coord)
  {
    unsigned int tid = coord.getParticleIndex();
    //d_ib[tid] = (coord[6] - 1) / idSlotsPerBunch;
    d_ib[tid] = coord[bunchIndex];
  }
};

class gpu_index_bunch_assignments_kernel2
{
 public:
  unsigned int idSlotsPerBunch, nBunches, ibMin;
  double *d_ibParticle;

 gpu_index_bunch_assignments_kernel2(unsigned int idSlotsPerBunch,
                                          unsigned int nBunches, unsigned int ibMin, double *d_ibParticle) : idSlotsPerBunch(idSlotsPerBunch), nBunches(nBunches), ibMin(ibMin),
    d_ibParticle(d_ibParticle){};

  __device__ void operator()(gpuParticleAccessor &coord)
  {
    unsigned int tid = coord.getParticleIndex();
    unsigned int ib;
    if (nBunches == 1)
      ib = 0;
    else {
      //ib = (coord[6] - 1) / idSlotsPerBunch - ibMin;
      ib = coord[bunchIndex] - ibMin;
    }
    d_ibParticle[tid] = ib;
  }
};

void gpu_index_bunch_assignments(long np, long idSlotsPerBunch,
                                      double P0, double **d_time, long **npBunch, long *nBunches, long lastNBunches)
{
  long ibMin, ibMax, ib;
  double fibMin, fibMax;
  double *fnpBunch = NULL;
  struct GPUBASE *gpuBase = getGpuBase();
  unsigned int particlePitch = gpuBase->gpu_array_pitch;
  double *d_ibParticle, *d_npBunch;
#if USE_MPI
  long ibMinGlobal, ibMaxGlobal;
#endif

#ifdef DEBUG
  printf("gpu_index_bunch_assignments called, np=%ld, idSlotsPerBunch=%ld, lastNBunches=%ld\n",
         np, idSlotsPerBunch, lastNBunches);
  printf("pointer check: time %s, npBunch %s, nBunches %s\n",
         d_time ? "ok" : "NULL",
         npBunch ? "ok" : "NULL",
         nBunches ? "ok" : "NULL");
  fflush(stdout);
#endif

  if (trajectoryTracking)
    {
      /* when tracking for trajectory, ignore bunch assignments */
      idSlotsPerBunch = -1;
    }

  if (idSlotsPerBunch <= 0)
    {
      ibMin = 0;
      *nBunches = 1;
#ifdef DEBUG
      printf("gpu_index_bunch_assignments: only one bunch\n");
      fflush(stdout);
#endif
    }
  else
    {
      sortByPID(np); // sort particles
      d_ibParticle = gpuBase->d_temp_particles + 4 * particlePitch;
      ibMin = LONG_MAX;
      ibMax = LONG_MIN;
      gpuDriver(np,
                gpu_index_bunch_assignments_kernel1(idSlotsPerBunch, d_ibParticle));
      gpuReduceMinMax(d_ibParticle, np, &fibMin, &fibMax);
      ibMin = (unsigned int)fibMin;
      ibMax = (unsigned int)fibMax;
#if USE_MPI
#  ifdef DEBUG
      printf("Sharing bunch min/max data: %ld, %ld\n", ibMin, ibMax);
      fflush(stdout);
#  endif
      MPI_Allreduce(&ibMin, &ibMinGlobal, 1, MPI_LONG, MPI_MIN, workers);
      MPI_Allreduce(&ibMax, &ibMaxGlobal, 1, MPI_LONG, MPI_MAX, workers);
      ibMin = ibMinGlobal;
      ibMax = ibMaxGlobal;
#endif
      if (ibMin == LONG_MAX || ibMax == LONG_MIN)
        *nBunches = 0;
      else
        *nBunches = (ibMax - ibMin) + 1;
#ifdef DEBUG
      printf("nPPB=%ld, ibMin = %ld, ibMax = %ld, nBunches = %ld\n", idSlotsPerBunch, ibMin, ibMax, *nBunches);
      fflush(stdout);
#endif
    }
  /* set pointers after the sort as the temp array is swapped. */
  d_npBunch = gpuBase->d_temp_particles + 5 * particlePitch;
  *d_time = gpuBase->d_temp_particles + 6 * particlePitch;

  if (*nBunches)
    {
      /* To prevent problems in LRWAKE, need to ensure that number of bunches does not increase */
      if (lastNBunches > 0 && *nBunches > lastNBunches)
        {
#ifdef DEBUG
          printf("Error: lastNBunches = %ld, *nBunches = %ld\n", lastNBunches, *nBunches);
          fflush(stdout);
#endif
#if USE_MPI
          mpiAbort = MPI_ABORT_BUNCH_ASSIGNMENT_ERROR;
          return;
#else
          bombElegant("Error: number of bunches has increased.", NULL);
#endif
        }

#ifdef DEBUG
      printf("Performing bunch assignment, nBunches=%ld\n", *nBunches);
      fflush(stdout);
#endif
#if USE_MPI
      if (isSlave || !notSinglePart)
        {
#  ifdef DEBUG
          printf("...performing bunch assignment\n");
          fflush(stdout);
#  endif
#endif
          if (np)
            {
#ifdef DEBUG
              printf("Doing branch for np!=0 (np=%ld)\n", np);
              fflush(stdout);
#endif
              /* Compute time coordinate of each particle */
              gpu_computeTimeCoordinates(np, *d_time, P0);
#ifdef DEBUG
              printf("Computed time coordinates\n");
              fflush(stdout);
#endif
              if (npBunch)
                {
#ifdef DEBUG
                  printf("Allocating npBunch array\n");
                  fflush(stdout);
#endif
                  *npBunch = (long *)malloc(sizeof(*npBunch) * (*nBunches));
                  /* No kernel launches for 1 bunch */
                  if (*nBunches == 1)
                    {
                      *npBunch[0] = np;
                    }
                  else
                    {
                      gpuDriver(np, gpu_index_bunch_assignments_kernel2(idSlotsPerBunch,
                                                                             *nBunches, ibMin, d_ibParticle));
                      gpu_binTimeDistribution_and_countBinned(d_npBunch, d_ibParticle, np, 0.0,
                                                              1.0, *nBunches);
                      fnpBunch = (double *)malloc(sizeof(double) * (*nBunches));
                      cudaMemcpy(fnpBunch, d_npBunch, sizeof(double) * (*nBunches), cudaMemcpyDeviceToHost);
                      for (ib = 0; ib < (*nBunches); ib++)
                        *npBunch[ib] = (long)fnpBunch[ib];
                      free(fnpBunch);
                    }
                }
            }

#if USE_MPI
        }
#endif
    }

#ifdef DEBUG
  printf("%ld bunches found\n", *nBunches);
  for (ib = 0; ib < *nBunches; ib++)
    printf("npBunch[%ld] = %ld\n", ib, (*npBunch)[ib]);
#endif

#if USE_MPI
#  ifdef DEBUG
  printf("Waiting on barrier at end of index_bunch_assignment\n");
  fflush(stdout);
#  endif
  if (notSinglePart)
    MPI_Barrier(workers);
  else
    MPI_Barrier(MPI_COMM_WORLD);
#endif

#ifdef DEBUG
  printf("Leaving gpu_index_bunch_assignment\n");
  fflush(stdout);
#endif
}
