#ifndef GPU_WAKE_H
#define GPU_WAKE_H

#include <gpu_track.h>
#include <gpu_base.h>

#ifdef __cplusplus
extern "C" {
#endif

void gpu_track_through_wake(unsigned int np0, WAKE *wakeData, double *PoInput,
		                        RUN *run, long i_pass, CHARGE *charge);

#ifdef __cplusplus
}
#endif

void gpu_applyLongitudinalWakeKicks(double* d_particles,
       unsigned int particlePitch, unsigned int np, double* d_time, double Po,
       double* d_Vtime, unsigned int nb, double tmin, double dt, int interpolate,
       double particleMassMvTimesParticleRelSign);

void gpu_applyLongitudinalWakeKicksWithFactor(double* d_particles, 
       unsigned int particlePitch, unsigned int np, double* d_time, double Po, 
       double* d_Vtime, unsigned int nb, double tmin, double dt, int interpolate, 
       double particleMassMvTimesParticleRelSign, double factor);

void gpu_applyLongitudinalWakeKicksAndDrift(double* d_particles,
       unsigned int particlePitch, unsigned int np, double* d_time, double Po,
       double* d_Vtime, unsigned int nb, double tmin, double dt, int interpolate,
       double particleMassMvTimesParticleRelSign, double length);

#endif /* GPU_WAKE_H */
