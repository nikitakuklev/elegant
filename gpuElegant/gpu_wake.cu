#include <gpu_track.h>
#include <table.h>
#include <fftpackC.h>

#include <gpu_base.h>
#include <gpu_particle_template_function.hcu>
#include <gpu_bin_time_distribution.h>
#include <gpu_reductions.h>
#include <gpu_wake.h>
#include <gpu_convolve_arrays.h> // gpuConvolveArrays
#include <gpu_funcs.h>           // gpu_do_match_energy
#include <gpu_lsc.h>             // gpu_applyLongitudinalWakeKicksWithFactor
#include <gpu_trwake.h>          // gpu_computeTimeCoordinatesAndMinMax
#include <gpu_lrwake.h>          // gpu_determine_bucket_assignments
#include <gpu_simple_rfca.h>     // gpu_add_to_particle_energy

extern "C"
{

  void set_up_wake(WAKE *wakeData, RUN *run, long pass, long particles, CHARGE *charge);

  void gpu_track_through_wake(unsigned int np0, WAKE *wakeData, double *PoInput,
                              RUN *run, long i_pass, CHARGE *charge)
  {
    double *d_Itime = NULL; /* array for histogram of particle density */
    double *Itime = NULL;   /* array for histogram of particle density */
    double *d_Vtime = NULL; /* array for voltage acting on each bin */
    double *d_time = NULL;  /* array to record arrival time of each particle */
    long *npBucket = NULL;  /* array to record how many particles are in each bucket */
    short shortBunchWarning = 0;
    long nb = 0, n_binned = 0, nb_max = 0;
    long iBucket, nBuckets, np;
    double factor, tmin, tmax, tmean = 0, dt = 0, Po, rampFactor;
    unsigned int offset = 0;
#if USE_MPI
    double *buffer;
#  if MPI_DEBUG
    printf("myid=%d, np0=%ld\n", myid, np0);
    fflush(stdout);
#  endif
#endif

    struct GPUBASE *gpuBase = getGpuBase();
    double *d_particles, *d_temp_particles, *d_WakeData_W;
    unsigned int particlePitch = gpuBase->gpu_array_pitch;
    bool allocGpuMem = false;

    set_up_wake(wakeData, run, i_pass, np0, charge);
    cudaMalloc((void **)&d_WakeData_W, sizeof(double) * wakeData->wakePoints);
    cudaMemcpy(d_WakeData_W, wakeData->W, sizeof(double) * wakeData->wakePoints, cudaMemcpyHostToDevice);

    if (isSlave || !notSinglePart)
      {
        rampFactor = 0;
        if (i_pass >= (wakeData->rampPasses - 1))
          rampFactor = 1;
        else
          rampFactor = (i_pass + 1.0) / wakeData->rampPasses;
        Po = *PoInput;

        gpu_determine_bucket_assignments(np0,
                                         (charge && wakeData->bunchedBeamMode) ? charge->idSlotsPerBunch : 0,
                                         Po, &d_time, &npBucket, &nBuckets, -1);
        gpuErrorHandler("gpu_track_through_wake::gpu_determine_bucket_assignments");
        /* set pointers after potential sort in gpu_determine_bucket_assignments */
        d_particles = gpuBase->d_particles;
        d_temp_particles = gpuBase->d_temp_particles;
        d_Itime = d_temp_particles + 2 * particlePitch;
        d_Vtime = d_temp_particles + 3 * particlePitch;

#ifdef DEBUG
        if (nBuckets > 1)
          {
            printf("%ld buckets\n", nBuckets);
            fflush(stdout);
            for (iBucket = 0; iBucket < nBuckets; iBucket++)
              {
                printf("bucket %ld: %ld particles\n", iBucket, npBucket[iBucket]);
                fflush(stdout);
              }
          }
#endif

        offset = 0;
        for (iBucket = 0; iBucket < nBuckets; iBucket++)
          {
            if (nBuckets == 1)
              {
                np = np0;
                tmean = gpuReduceAdd(d_time, np);
                if (np > 0)
                  tmean /= np;
              }
            else
              {
                if ((np = npBucket[iBucket]) == 0)
                  continue;
#ifdef DEBUG
                printf("WAKE: pointing to data work array, iBucket=%ld, np=%ld\n", iBucket, np);
                fflush(stdout);
#endif
                selectParticleSubset(offset);
                tmean = gpuReduceAdd(d_time + offset, np);
                if (np > 0)
                  tmean /= np;
              }

            gpuReduceMinMax(d_time + offset, np, &tmin, &tmax);
#ifdef DEBUG
            printf("WAKE: tmin=%21.15le, tmax=%21.15le, np=%ld\n", tmin, tmax, np);
            fflush(stdout);
#endif
#if USE_MPI
            if (isSlave && notSinglePart)
              find_global_min_max(&tmin, &tmax, np, workers);
#  ifdef DEBUG
            printf("WAKE: global tmin=%21.15le, tmax=%21.15le, np=%ld\n", tmin, tmax, np);
            fflush(stdout);
#  endif
#endif
            if (isSlave || !notSinglePart)
              {
                if ((tmax - tmin) > (wakeData->t[wakeData->wakePoints - 1] - wakeData->t[0]))
                  {
                    if (!wakeData->allowLongBeam)
                      {
                        fprintf(stderr, "Error: The beam is longer than the longitudinal wake function.\nThis may produce unphysical results.\n");
                        fprintf(stderr, "The beam length is %le s, while the wake length is %le s\n",
                                tmax - tmin, wakeData->t[wakeData->wakePoints - 1] - wakeData->t[0]);
                        exit(1);
                      }
                    printf("Warning: The beam is longer than the longitudinal wake function.\nThis may produce unphysical results.\n");
                    printf("The beam length is %le s, while the wake length is %le s\n",
                           tmax - tmin, wakeData->t[wakeData->wakePoints - 1] - wakeData->t[0]);
                    /*
                      if (abs(tmax-tmean)<abs(tmin-tmean)) 
                      tmin = tmax - (wakeData->t[wakeData->wakePoints-1]-wakeData->t[0]);
                      else
                      tmax = tmin + (wakeData->t[wakeData->wakePoints-1]-wakeData->t[0]);
                    */
                  }

                dt = wakeData->dt;
                if (np > 1 && (tmax - tmin) < 20 * dt && !shortBunchWarning)
                  {
                    printf("Warning: The beam is shorter than 20*DT, where DT is the spacing of the wake points.\n");
                    printf("         Depending on the longitudinal distribution and shape of the wake, this may produce poor results.\n");
                    printf("         Consider using a wake with finer time spacing in WAKE elements.\n");
                    fflush(stdout);
                    shortBunchWarning = 1;
                  }

                if (wakeData->n_bins)
                  {
                    nb = wakeData->n_bins;
                    tmin = tmean - dt * nb / 2.0;
                  }
                else
                  {
                    nb = (tmax - tmin) / dt + 3;
                    tmin -= dt;
                    tmax += dt;
                  }
                if (nb > particlePitch && nb > nb_max)
                  {
                    if (allocGpuMem)
                      {
                        cudaFree(d_Itime);
                        cudaFree(d_Vtime);
                      }
                    cudaMalloc((void **)&d_Itime, sizeof(double) * nb);
                    cudaMalloc((void **)&d_Vtime, sizeof(double) * nb);
                    allocGpuMem = true;
                    nb_max = nb;
                  }
                if (nb <= 0)
                  {
                    printf("Warning: Number of wake bins is 0 or negative\n");
                    printf("probably indicating an extremely long bunch\n");
                    printf("Wake ignored!\n");
                    return;
                  }
#ifdef DEBUG
                printf("WAKE: dt=%le, nb=%ld\n", dt, nb);
                fflush(stdout);
#endif

                n_binned = gpu_binTimeDistribution_and_countBinned(d_Itime, d_time + offset, np, tmin, dt, nb);
                gpuErrorHandler("gpu_track_through_wake::gpu_binTimeDistribution_and_countBinned");
              }

            if (USE_MPI || isSlave)
              {
                Itime = (double *)malloc(sizeof(double) * nb);
                cudaMemcpy(Itime, d_Itime, sizeof(double) * nb, cudaMemcpyDeviceToHost);
                gpuErrorHandler("gpu_track_through_wake::cudaMemcpy d_Itime");
              }

            if (!USE_MPI || !notSinglePart)
              {
                if (n_binned != np)
                  {
                    printf("warning: only %ld of %ld particles were binned (WAKE)\n", n_binned, np);
                    printf("consider setting n_bins=0 in WAKE definition to invoke autoscaling\n");
                    fflush(stdout);
                    return;
                  }
              }
#if USE_MPI
            else if (isSlave)
              {
                int all_binned, result = 1;

                result = ((n_binned == np) ? 1 : 0);
                MPI_Allreduce(&result, &all_binned, 1, MPI_INT, MPI_LAND, workers);
                if (!all_binned)
                  {
                    if (myid == 1)
                      {
                        /* This warning will be given only if the flag MPI_DEBUG is defined for the Pelegant */
                        printf("warning: Not all of %ld particles were binned (WAKE)\n", np);
                        printf("consider setting n_bins=0 in WAKE definition to invoke autoscaling\n");
                        fflush(stdout);
                      }
                  }
              }

            if (isSlave && notSinglePart)
              {
                buffer = (double *)malloc(sizeof(double) * nb);
                MPI_Allreduce(Itime, buffer, nb, MPI_DOUBLE, MPI_SUM, workers);
                memcpy(Itime, buffer, sizeof(double) * nb);
                free(buffer);
              }
#endif
            if (isSlave || !notSinglePart)
              {
                if (wakeData->smoothing && nb >= (2 * wakeData->SGHalfWidth + 1))
                  {
                    if (!SavitzyGolaySmooth(Itime, nb, wakeData->SGOrder, wakeData->SGHalfWidth, wakeData->SGHalfWidth, 0))
                      {
                        fprintf(stderr, "Problem with smoothing for WAKE element (file %s)\n",
                                wakeData->inputFile);
                        fprintf(stderr, "Parameters: nbins=%ld, order=%ld, half-width=%ld\n",
                                nb, wakeData->SGOrder, wakeData->SGHalfWidth);
                        exit(1);
                      }
                  }
              }

            if (USE_MPI || isSlave)
              {
                cudaMemcpy(d_Itime, Itime, sizeof(double) * nb, cudaMemcpyHostToDevice);
                gpuErrorHandler("gpu_track_through_wake::cudaMemcpy d_Itime (2)");
                free(Itime);
                Itime = NULL;
              }

            /* Do the convolution of the particle density and the wake function,
               V(T) = Integral[W(T-t)*I(t)dt, t={-infinity, T}]
               Note that T<0 is the head of the bunch.
               For the wake, the argument is the normal convention wherein larger
               arguments are later times.
            */
            if (isSlave || !notSinglePart)
              {
                if (wakeData->i0 == 0) 
                  {
                    gpuConvolveArrays(d_Vtime, nb, d_Itime, nb, d_WakeData_W, wakeData->wakePoints);
                    gpuErrorHandler("gpu_track_through_wake::gpuConvolveArrays");
                  } 
                else 
                  {
                    /*FIX THIS*/
                    exit(0);
                  }
                
                factor = wakeData->macroParticleCharge * particleRelSign * wakeData->factor * rampFactor;

                gpu_applyLongitudinalWakeKicksWithFactor(d_particles + offset, particlePitch, np,
                                                         d_time + offset, Po, d_Vtime, nb, tmin, dt, wakeData->interpolate,
                                                         particleMassMV * particleRelSign, factor);
                gpuErrorHandler("gpu_track_through_wake::gpu_applyLongitudinalWakeKicksWithFactor");
              }

            /* No copy back with gpu code as we set offset */
            offset += np; // for next iteration
          }

        if (nBuckets > 1)
          {
            free(npBucket);
            npBucket = NULL;
          }
      }

#if USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    resetToFullParticleArray();
    gpuErrorHandler("gpu_track_through_wake::resetToFullParticleArray");
    cudaFree(d_WakeData_W);

    if (wakeData->change_p0)
      gpu_do_match_energy(np, PoInput, 0);
    gpuErrorHandler("gpu_track_through_wake::gpu_do_match_energy");

    if (npBucket)
      free(npBucket);
    if (Itime)
      free(Itime);

    if (allocGpuMem)
      {
        cudaFree(d_Itime);
        cudaFree(d_Vtime);
      }
  }

} // extern "C"

class gpuApplyLongitudinalWakeKicks : public gpuFunctionIfc
{
  double *d_time;
  double *d_Vtime;
  unsigned int nb;
  double tmin;
  double dt;
  double Po;
  bool interpolate;
  double particleMassMvTimesParticleRelSign;
  double factor;

 public:
  gpuApplyLongitudinalWakeKicks(double *d_time, double *d_Vtime,
                                unsigned int nb, double tmin, double dt, double Po, bool interpolate,
                                double particleMassMvTimesParticleRelSign, double factor = 1.0)
    {
      this->d_time = d_time;
      this->d_Vtime = d_Vtime;
      this->nb = nb;
      this->tmin = tmin;
      this->dt = dt;
      this->Po = Po;
      this->interpolate = interpolate;
      this->particleMassMvTimesParticleRelSign = particleMassMvTimesParticleRelSign;
      this->factor = factor;
    }
  __device__ inline void operator()(gpuParticleAccessor &particle)
  {
    double r_time = d_time[particle.getParticleIndex()];
    unsigned int ib = (r_time - tmin) / dt + 0.5;
    double dgam = 0;

    if (ib < nb)
      {
        if (interpolate)
          {
            double dt1 = r_time - (tmin + dt * ib);
            if ((dt1 < 0 && ib) || ib == nb - 1)
              {
                ib--;
                dt1 += dt;
              }
            dgam = (factor * d_Vtime[ib] + (factor * d_Vtime[ib + 1] - factor * d_Vtime[ib]) / dt * dt1) / (1e6 * particleMassMvTimesParticleRelSign);
          }
        else
          {
            dgam = factor * d_Vtime[ib] / (1e6 * particleMassMvTimesParticleRelSign);
          }
        if (dgam)
          {
            gpu_add_to_particle_energy(particle, r_time, Po, -dgam);
          }
      }
  }
};

void gpu_applyLongitudinalWakeKicks(double *d_particles,
                                    unsigned int particlePitch, unsigned int np, double *d_time, double Po,
                                    double *d_Vtime, unsigned int nb, double tmin, double dt, int interpolate,
                                    double particleMassMvTimesParticleRelSign)
{
  gpuDriver(np, gpuApplyLongitudinalWakeKicks(
                                              d_time, d_Vtime, nb, tmin, dt, Po, interpolate,
                                              particleMassMvTimesParticleRelSign));
}

void gpu_applyLongitudinalWakeKicksWithFactor(double *d_particles,
                                              unsigned int particlePitch, unsigned int np, double *d_time, double Po,
                                              double *d_Vtime, unsigned int nb, double tmin, double dt, int interpolate,
                                              double particleMassMvTimesParticleRelSign, double factor)
{
  gpuDriver(np, gpuApplyLongitudinalWakeKicks(
                                              d_time, d_Vtime, nb, tmin, dt, Po, interpolate,
                                              particleMassMvTimesParticleRelSign, factor));
}

void gpu_applyLongitudinalWakeKicksAndDrift(double *d_particles,
                                            unsigned int particlePitch, unsigned int np, double *d_time, double Po,
                                            double *d_Vtime, unsigned int nb, double tmin, double dt, int interpolate,
                                            double particleMassMvTimesParticleRelSign, double length)
{
  gpuDriver(np, gpuApplyLongitudinalWakeKicks(d_time, d_Vtime, nb, tmin, dt, Po, interpolate, particleMassMvTimesParticleRelSign), gpuExactDrift(length));
}
