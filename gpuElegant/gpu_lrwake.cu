#include <gpu_track.h>
#include <gpu_base.h>
#include <gpu_particle_template_function.hcu>
#include <gpu_reductions.h>
#include <gpu_trwake.h> // gpu_computeTimeCoordinates
#include <gpu_bin_time_distribution.h> // gpu_binTimeDistribution_and_countBinned

class gpu_determine_bucket_assignments_kernel1{
public:
  unsigned int idSlotsPerBunch;
  double* d_ib;

  gpu_determine_bucket_assignments_kernel1(unsigned int idSlotsPerBunch,
    double* d_ib) : idSlotsPerBunch(idSlotsPerBunch), d_ib(d_ib) {};

  __device__ void operator()(gpuParticleAccessor& coord){
    unsigned int tid = coord.getParticleIndex();
    d_ib[tid]=(coord[6]-1)/idSlotsPerBunch;
  }
};

class gpu_determine_bucket_assignments_kernel2{
public:
  unsigned int idSlotsPerBunch, nBuckets, ibMin;
  double* d_ibParticle;

  gpu_determine_bucket_assignments_kernel2(unsigned int idSlotsPerBunch,
    unsigned int nBuckets, unsigned int ibMin, double* d_ibParticle) : 
    idSlotsPerBunch(idSlotsPerBunch), nBuckets(nBuckets), ibMin(ibMin),
    d_ibParticle(d_ibParticle) {};

  __device__ void operator()(gpuParticleAccessor& coord){
    unsigned int tid = coord.getParticleIndex();
    unsigned int ib;
    if (nBuckets==1)
      ib = 0;
    else
      ib = (coord[6]-1)/idSlotsPerBunch - ibMin;
    d_ibParticle[tid] = ib;
  }
};

void gpu_determine_bucket_assignments(long np, long idSlotsPerBunch, 
       double P0, double **d_time, long **npBucket, long *nBuckets, long lastNBuckets)
{
  long ibMin, ibMax, ib;
  double fibMin, fibMax;
  double* fnpBucket = NULL;
  struct GPUBASE* gpuBase = getGpuBase();
  unsigned int particlePitch = gpuBase->gpu_array_pitch;
  double* d_ibParticle, *d_npBucket;
#if USE_MPI
  long ibMinGlobal, ibMaxGlobal;
#endif

#ifdef DEBUG
  printf("gpu_determine_bucket_assignments called, np=%ld, idSlotsPerBunch=%ld, lastNBuckets=%ld\n",
         np, idSlotsPerBunch, lastNBuckets);
  printf("pointer check: time %s, npBucket %s, nBuckets %s\n",
         d_time ? "ok" : "NULL",
         npBucket ? "ok" : "NULL",
         nBuckets ? "ok" : "NULL");
  fflush(stdout);
#endif

  if (trajectoryTracking) {
    /* when tracking for trajectory, ignore bunch assignments */
    idSlotsPerBunch = -1;
  }

  if (idSlotsPerBunch<=0) {
    ibMin = 0;
    *nBuckets = 1;
#ifdef DEBUG
    printf("gpu_determine_bucket_assignments: only one bunch\n");
    fflush(stdout);
#endif
  } else {
    sortByPID(np); // sort particles 
    d_ibParticle = gpuBase->d_temp_particles + 4*particlePitch;
    ibMin = UINT_MAX;
    ibMax = 0;
    gpuDriver(np,
      gpu_determine_bucket_assignments_kernel1(idSlotsPerBunch, d_ibParticle));
    gpuReduceMinMax(d_ibParticle, np, &fibMin, &fibMax);
    ibMin = (unsigned int)fibMin; 
    ibMax = (unsigned int)fibMax; 
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
  /* set pointers after the sort as the temp array is swapped. */
  d_npBucket = gpuBase->d_temp_particles + 5*particlePitch;
  *d_time =    gpuBase->d_temp_particles + 6*particlePitch;

  if (*nBuckets) {
    /* To prevent problems in LRWAKE, need to ensure that number of buckets does not increase */
    if (lastNBuckets>0 && *nBuckets>lastNBuckets) {
#ifdef DEBUG
      printf("Error: lastNBuckets = %ld, *nBuckets = %ld\n", lastNBuckets, *nBuckets);
      fflush(stdout);
#endif
#if USE_MPI
      mpiAbort = MPI_ABORT_BUCKET_ASSIGNMENT_ERROR;
      return;
#else
      bombElegant("Error: number of bunches has increased.", NULL);
#endif
    }

#ifdef DEBUG
    printf("Performing bucket assignment, nBuckets=%ld\n", *nBuckets);
    fflush(stdout);
#endif
#if USE_MPI
    if (isSlave || !notSinglePart) {
#ifdef DEBUG
      printf("...performing bucket assignment\n");
      fflush(stdout);
#endif
#endif
      if (np) {
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
        if (npBucket) {
#ifdef DEBUG
          printf("Allocating npBucket array\n");
          fflush(stdout);
#endif
          *npBucket = (long*) malloc(sizeof(*npBucket)*(*nBuckets));
          /* No kernel launches for 1 bucket */
          if (*nBuckets==1) {
            *npBucket[0] = np;
          } else {
            gpuDriver(np, gpu_determine_bucket_assignments_kernel2(idSlotsPerBunch,
                                              *nBuckets, ibMin, d_ibParticle));
            gpu_binTimeDistribution_and_countBinned(d_npBucket, d_ibParticle, np, 0.0,
                                                1.0, *nBuckets);
            fnpBucket = (double*) malloc(sizeof(double)*(*nBuckets));
            cudaMemcpy(fnpBucket, d_npBucket, sizeof(double)*(*nBuckets), cudaMemcpyDeviceToHost);
            for (ib=0; ib<(*nBuckets); ib++)
              *npBucket[ib] = (long)fnpBucket[ib];
            free(fnpBucket);
          }
        }
      }

#if USE_MPI
    }
#endif

  }

#ifdef DEBUG
  printf("%ld buckets found\n", *nBuckets);
  for (ib=0; ib<*nBuckets; ib++) 
    printf("npBucket[%ld] = %ld\n", ib, (*npBucket)[ib]);
#endif

#if USE_MPI
#ifdef DEBUG
  printf("Waiting on barrier at end of determine_bucket_assignment\n");
  fflush(stdout);
#endif
  if (notSinglePart)
    MPI_Barrier(workers);
  else
    MPI_Barrier(MPI_COMM_WORLD);
#endif

#ifdef DEBUG
    printf("Leaving gpu_determine_bucket_assignment\n");
    fflush(stdout);
#endif
}
