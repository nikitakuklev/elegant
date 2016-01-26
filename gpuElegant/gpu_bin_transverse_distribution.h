#ifndef GPU_BIN_TRANSVERSE_DISTRIBUTION_H
#define GPU_BIN_TRANSVERSE_DISTRIBUTION_H

#include <gpu_track.h>

#ifdef __cplusplus /* if it's C++ use C linkage */
extern "C" {
#endif
unsigned int gpu_binTransverseDistribution(double* d_Itime0, double* d_Itime1, const unsigned int& nb, double* d_particles, unsigned int particlePitch, const unsigned int& np,
			      double* d_time, const double& tmin, const double& dt, const double& Po, const double& dx, const double& dy, 
					   const int& xPower, const int& yPower, double* d_pz);
#ifdef __cplusplus
}
#endif
#endif
