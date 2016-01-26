#include <gpu_track.h>

#include <gpu_base.h>
#include <gpu_particle_template_function.hcu>
#include <gpu_particle_reduction.hcu>
#include <gpu_bin_transverse_distribution.h>
#include <gpu_reductions.h>
#include <gpu_trwake.h>
#include <gpu_convolve_arrays.h> // gpuConvolveArrays
#include <gpu_funcs.h> // gpu_rotateBeamCoordinates
#include <gpu_lrwake.h> // gpu_determine_bucket_assignments

#define c_mks   (2.99792458e8)

class gpuApplyTransverseWakeKicks{
public:
  double* d_time;
  double* d_pz;
  double Po;
  unsigned int plane;
  double* d_Vtime;
  unsigned int nb;
  double tmin;
  double dt;
  int interpolate;
  int exponent;
  double particleMassMV;
  double particleRelSign;
  double factor;

  gpuApplyTransverseWakeKicks(double* d_time, double* d_pz, double Po,
      unsigned int plane, double* d_Vtime, unsigned int nb, double tmin,
      double dt, int interpolate, int exponent, double particleMassMV,
      double particleRelSign, double factor) : d_time(d_time), d_pz(d_pz),
      Po(Po), plane(plane), d_Vtime(d_Vtime), nb(nb), tmin(tmin), dt(dt),
      interpolate(interpolate), exponent(exponent),
      particleMassMV(particleMassMV), particleRelSign(particleRelSign),
      factor(factor) {}

  __device__ void inline operator()(gpuParticleAccessor& particle){
    int offset;
   
    unsigned int ip = particle.getParticleIndex();

    offset = 2*plane+1;  /* xp or yp index */

    /* Bin CENTERS are at tmin+ib*dt */
    int ib = (d_time[ip]-tmin)/dt+0.5;
    double Vinterp = 0;
    if( ib <= nb -1 && ib >= 0){
      if (interpolate) {
        double dt1 = d_time[ip]-(tmin+dt*ib); /* distance to bin center */
        if ((dt1<0 && ib) || ib==nb-1) {
          ib--;
          dt1 += dt;
        }
        Vinterp = factor*d_Vtime[ib]+(factor*d_Vtime[ib+1]-factor*d_Vtime[ib])/dt*dt1;
      } else
        Vinterp = factor*d_Vtime[ib];
      if (exponent)
        Vinterp *= pow(particle[offset-1], exponent);
      if (Vinterp)
        particle[offset] += Vinterp/(1e6*particleMassMV*particleRelSign)/d_pz[ip];
    }//ib<=nb
  }
};

extern "C" {

void set_up_trwake(TRWAKE *wakeData, RUN *run, long pass, long particles, CHARGE *charge);

void gpu_track_through_trwake(long np0,
			      TRWAKE *wakeData, double Po,
			      RUN *run, long i_pass, CHARGE *charge)
{
  double *posItime[2] = {NULL, NULL};     /* array for histogram of particle density times x, y*/
  double *Vtime = NULL;           /* array for voltage acting on each bin */
  long max_n_bins = 0;
  long *npBucket = NULL;          /* array to record how many particles are in each bucket */
  long nb=0, n_binned=0, plane, nb_max=0;
  long iBucket, nBuckets, np;
  double factor, tmin, tmean=0, tmax, dt=0, rampFactor=1;
  unsigned int offset=0;
#if USE_MPI
  double *buffer;
#endif
  double* d_particles, *d_wakeData, *d_pz, *d_time, *d_Vtime;  
  unsigned int particlePitch = getGpuBase()->gpu_array_pitch;
  double* d_posItime[2] = {NULL, NULL};
  bool allocGpuMem=false;

  set_up_trwake(wakeData, run, i_pass, np0, charge);
  cudaMalloc( (void**)&d_wakeData, sizeof(double)*wakeData->wakePoints);

  if (i_pass>=(wakeData->rampPasses-1))
    rampFactor = 1;
  else
    rampFactor = (i_pass+1.0)/wakeData->rampPasses;

  if (isSlave || !notSinglePart) {
    gpu_determine_bucket_assignments(np0, 
      (charge && wakeData->bunchedBeamMode)?charge->idSlotsPerBunch:0,
      Po, &d_time, &npBucket, &nBuckets, -1);
    gpuErrorHandler("gpu_track_through_trwake::gpu_determine_bucket_assignments");
    /* set pointers after potential sort in gpu_determine_bucket_assignments */
    d_particles =   getGpuBase()->d_particles;
    d_pz =          getGpuBase()->d_temp_particles + 2*particlePitch;
    d_Vtime =       getGpuBase()->d_temp_particles + 3*particlePitch;
    d_posItime[0] = getGpuBase()->d_temp_particles + 4*particlePitch;
    d_posItime[1] = getGpuBase()->d_temp_particles + 5*particlePitch;

    offset=0;
    for (iBucket=0; iBucket<nBuckets; iBucket++) {
      if (nBuckets==1) {
        np = np0;
        tmean = gpuReduceAdd(d_time, np);
        if (np>0)
          tmean /= np;
      } else {
        if ((np = npBucket[iBucket])==0)
          continue;
#ifdef DEBUG
        printf("WAKE: copying data to work array, iBucket=%ld, np=%ld\n", iBucket, np);
        fflush(stdout);
#endif
        selectParticleSubset(offset);
        tmean = gpuReduceAdd(d_time+offset, np);
        if (np>0)
          tmean /= np;
      }

      gpuReduceMinMax(d_time+offset, np, &tmin, &tmax);
      gpuErrorHandler("gpu_track_through_trwake::gpuReduceMinMax");
#ifdef DEBUG
      printf("WAKE: tmin=%21.15le, tmax=%21.15le, tmax-tmin=%21.15le, np=%ld\n", tmin, tmax, tmax-tmin, np);
      fflush(stdout);
#endif
#if USE_MPI
      if (isSlave && notSinglePart)
        find_global_min_max(&tmin, &tmax, np, workers);
#endif
      if (isSlave ||  !notSinglePart) {
        if ((tmax-tmin) > (wakeData->t[wakeData->wakePoints-1]-wakeData->t[0])) {
          fprintf(stderr, "The beam is longer than the transverse wake function.\nThis would produce unphysical results.\n");
          fprintf(stderr, "The beam length is %le s, while the wake length is %le s\n",
                  tmax-tmin, wakeData->t[wakeData->wakePoints-1]-wakeData->t[0]);
          exitElegant(1);
        }

        dt = wakeData->dt;

        if (wakeData->n_bins) {
          tmin = tmean-dt*wakeData->n_bins/2.0;
          nb = wakeData->n_bins;
        }
        else {
          nb = (tmax-tmin)/dt+3;
          tmin -= dt;
          tmax += dt;
        }

        if (tmin>tmax || nb<=0) {
          fprintf(stdout, "Problem with time coordinates in TRWAKE.  Po=%le\n", Po);
          exitElegant(1);
        }

        if (nb>max_n_bins) {
          posItime[0] = (double*) trealloc(posItime[0], sizeof(**posItime)*(max_n_bins=nb));
          posItime[1] = (double*) trealloc(posItime[1], sizeof(**posItime)*(max_n_bins=nb));
          Vtime = (double*) trealloc(Vtime, sizeof(*Vtime)*(max_n_bins+1));
        }

        if (nb > particlePitch && nb > nb_max) {
          if (allocGpuMem) {
            cudaFree(d_Vtime);
            cudaFree(d_posItime[0]);
            cudaFree(d_posItime[1]);
          }
          cudaMalloc( (void**)&d_Vtime, sizeof(double)*nb);
          cudaMalloc( (void**)&(d_posItime[0]), sizeof(double)*nb);
          cudaMalloc( (void**)&(d_posItime[1]), sizeof(double)*nb);
          allocGpuMem=true;
          nb_max=nb;
        }

        if (wakeData->tilt)
          gpu_rotateBeamCoordinates(np, wakeData->tilt);

        n_binned = gpu_binTransverseDistribution(d_posItime[0], d_posItime[1],
                     nb, d_particles+offset, particlePitch, np, d_time, tmin, dt, Po,
                     wakeData->dx, wakeData->dy, wakeData->xDriveExponent, 
                     wakeData->yDriveExponent, d_pz);
        gpuErrorHandler("gpu_track_through_trwake::gpu_binTransverseDistribution");
#if (!USE_MPI)
#ifdef DEBUG
        //dumpTransverseTimeDistributions(posItime, nb); /* no gpu implementation */
#endif
#endif
      }
#if (!USE_MPI)
      if (n_binned!=np) {
        fprintf(stdout, "warning: only %ld of %ld particles where binned (TRWAKE)\n", n_binned, np);
        fprintf(stdout, "consider setting n_bins=0 in TRWAKE definition to invoke autoscaling\n");
        fflush(stdout);
      }
#else
      if (notSinglePart) {
        if (isSlave) {
          int all_binned, result = 1;

          result = ((n_binned==np) ? 1 : 0);
          MPI_Allreduce(&result, &all_binned, 1, MPI_INT, MPI_LAND, workers);
          if (!all_binned) {
            if (myid==1) {
              /* This warning will be given only if the flag MPI_DEBUG is defined for the Pelegant */
              fprintf(stdout, "warning: Not all of %ld particles were binned (WAKE)\n", np);
              fprintf(stdout, "consider setting n_bins=0 in WAKE definition to invoke autoscaling\n");
              fflush(stdout);
            }
          }
        }
      } else {
        if (n_binned!=np) {
          fprintf(stdout, "warning: only %ld of %ld particles where binned (TRWAKE)\n", n_binned, np);
          fprintf(stdout, "consider setting n_bins=0 in TRWAKE definition to invoke autoscaling\n");
          fflush(stdout);
        }
      }
#endif
      if (isSlave || !notSinglePart) {
        for (plane=0; plane<2; plane++) {
          if (!wakeData->W[plane])
            continue;

          cudaMemcpy(posItime[plane], d_posItime[plane], sizeof(double)*nb,
                      cudaMemcpyDeviceToHost);
#if USE_MPI 
          if (isSlave && notSinglePart) {
            buffer = (double*) malloc(sizeof(double) * nb);
            MPI_Allreduce(posItime[plane], buffer, nb, MPI_DOUBLE, MPI_SUM, workers);
            memcpy(posItime[plane], buffer, sizeof(double)*nb);
            free(buffer); buffer = NULL;
          }
#endif
          factor = wakeData->macroParticleCharge*particleRelSign*wakeData->factor;
          if (plane==0)
            factor *= wakeData->xfactor*rampFactor;
          else
            factor *= wakeData->yfactor*rampFactor;
          if (!factor)
            continue;

          if (wakeData->smoothing && nb>=(2*wakeData->SGHalfWidth+1)) {
            if (!SavitzyGolaySmooth(posItime[plane], nb, wakeData->SGOrder,
                                    wakeData->SGHalfWidth, wakeData->SGHalfWidth, 0)) {
              fprintf(stderr, "Problem with smoothing for TRWAKE element (file %s)\n",
                      wakeData->inputFile);
              fprintf(stderr, "Parameters: nbins=%ld, order=%ld, half-width=%ld\n",
                      nb, wakeData->SGOrder, wakeData->SGHalfWidth);
              exitElegant(1);
            }

          }

          /* Do the convolution of the particle density and the wake function,
             V(T) = Integral[W(T-t)*I(t)dt, t={-infinity, T}]
             Note that T<0 is the head of the bunch.
             For the wake, the argument is the normal convention wherein larger
             arguments are later times.
             */
          cudaMemcpy(d_posItime[plane], posItime[plane], sizeof(double)*nb,
                     cudaMemcpyHostToDevice);
          cudaMemcpy(d_wakeData, wakeData->W[plane],
                     sizeof(double)*wakeData->wakePoints, cudaMemcpyHostToDevice);
          gpuErrorHandler("gpu_track_through_trwake::cudaMemcpy d_wakeData");
          gpuConvolveArrays(d_Vtime, nb, d_posItime[plane], nb, d_wakeData,
                            wakeData->wakePoints);
          gpuErrorHandler("gpu_track_through_trwake::gpuConvolveArrays");

          /* change particle transverse momenta to reflect voltage in relevant bin */
           gpuDriver(np,
               gpuApplyTransverseWakeKicks(d_time, d_pz, Po, plane, 
               d_Vtime, nb, tmin, dt, wakeData->interpolate,
               plane==0?wakeData->xProbeExponent:wakeData->yProbeExponent,
               particleMassMV, particleRelSign, factor) );
          gpuErrorHandler("gpu_track_through_trwake::gpuApplyTransverseWakeKicks");
        }

        if (wakeData->tilt)
          gpu_rotateBeamCoordinates(np, -wakeData->tilt);

      }
      /* No copy back with gpu code as we set offset */
      offset += np; // for next iteration
    }
  }

#if USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  if (npBucket)
    free(npBucket);
  if (Vtime)
    free(Vtime);
  if (posItime[0])
    free(posItime[0]);
  if (posItime[1])
    free(posItime[1]);

  resetToFullParticleArray();
  cudaFree(d_wakeData);
  if (allocGpuMem) {
    cudaFree(d_Vtime);
    cudaFree(d_posItime[0]);
    cudaFree(d_posItime[1]);
  }
}

} // extern "C"

class gpuComputeTimeCoordinatesFunctor {
public:
  double* d_time;
  double Po;
  gpuComputeTimeCoordinatesFunctor(double* time, double _Po) : 
    d_time(time), Po(_Po) {}
  
  __device__ inline double operator()(gpuParticleAccessor& particle){
    unsigned int tid = particle.getParticleIndex();

    double P = Po * (particle[5] + 1.0);
    double r_time = particle[4] * __dsqrt_rn((P*P)+1)/(c_mks*P);
    d_time[tid] = r_time;
    return r_time;
  }
};
  
void gpu_computeTimeCoordinates(unsigned int np, double* d_time, double Po) {
  gpu_computeTimeCoordinatesAndMinMax(np, d_time, Po, NULL, NULL);
}

void gpu_computeTimeCoordinatesAndMinMax(unsigned int np, double* d_time, 
    double Po, double* tmin, double* tmax){

  if (tmin==NULL || tmax==NULL) 
    gpuDriver(np, 
      gpuComputeTimeCoordinatesFunctor(d_time, Po));
  else  
    gpuParticleDoubleReduction(np,
      gpuComputeTimeCoordinatesFunctor(d_time, Po), Min<double>(), 
      Max<double>(), tmin, tmax);
}

double gpu_computeTimeCoordinatesAndMean(long np, double* d_time, 
                                         double Po) {
  double tmean;
#ifdef USE_KAHAN
  double error=0.0;
#endif

#if (!USE_MPI)
#ifndef USE_KAHAN
  tmean = gpuParticleReduction(np,
            gpuComputeTimeCoordinatesFunctor(d_time, Po), Add<double>());
#else
  tmean = gpuParticleKahanReduction(np, &error,
            gpuComputeTimeCoordinatesFunctor(d_time, Po));
#endif
  return tmean/np;
#else /* USE_MPI */
  if (!partOnMaster){
    long np_total;
    double tmean_total;
    if (isSlave || !notSinglePart) {
#ifndef USE_KAHAN
      tmean = gpuParticleReduction(np,
                gpuComputeTimeCoordinatesFunctor(d_time, Po), Add<double>());
#else
      tmean = gpuParticleKahanReduction(np, &error,
                gpuComputeTimeCoordinatesFunctor(d_time, Po));
#endif
    }
    tmean_total = DBL_MAX;
    if (notSinglePart) {
      if (isMaster) {
        tmean = 0;
        np = 0;
      }
      MPI_Allreduce(&np, &np_total, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
#ifndef USE_KAHAN
      MPI_Allreduce(&tmean, &tmean_total, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
      tmean_total = KahanParallel (tmean, error, MPI_COMM_WORLD);
#endif
    }
    if (tmean_total==DBL_MAX)
      bombElegant("invalid value of tmean_total in computeTimeCoordinates. Seek professional help!", NULL);
    if (np_total>0)
    return tmean_total/np_total;
    return (double)0;
  } 
  else { /* This serial part can be removed after all the upper level functions (e.g. wake) are parallelized */
#ifndef USE_KAHAN
    tmean = gpuParticleReduction(np,
              gpuComputeTimeCoordinatesFunctor(d_time, Po), Add<double>());
#else
    tmean = gpuParticleKahanReduction(np, &error,
              gpuComputeTimeCoordinatesFunctor(d_time, Po));
#endif
    return tmean/np;   
  }
#endif
}

class gpuComputeDistanceCoordinatesFunctor {
public:
  double* d_time;
  double Po;
  gpuComputeDistanceCoordinatesFunctor(double* d_time, double Po) : 
    d_time(d_time), Po(Po) {}
  
  __device__ inline void operator()(gpuParticleAccessor& particle) {
    unsigned int tid = particle.getParticleIndex();

    double P = Po*(particle[5]+1);
    double beta = P/sqrt(P*P+1);
    particle[4] = c_mks*beta*d_time[tid];
  }
};

void gpu_computeDistanceCoordinates(double *d_time, double Po, long np) {
  gpuDriver(np, gpuComputeDistanceCoordinatesFunctor(d_time, Po));
}
