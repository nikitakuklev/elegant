#ifndef GPU_LSC_H
#define GPU_LSC_H

#include <gpu_track.h>
#include <gpu_base.h>

#ifdef __cplusplus
extern "C" {
#endif

void gpu_track_through_lscdrift(unsigned int np, LSCDRIFT *LSC, double Po,
       CHARGE *charge);

void gpu_addLSCKick(double* d_particles, unsigned int particlePitch,
       unsigned int np, LSCKICK *LSC, double Po, CHARGE *charge, 
       double lengthScale, double dgammaOverGamma, double* d_temp_particles);

#ifdef __cplusplus
}
#endif

#endif /* GPU_LSC_H */
