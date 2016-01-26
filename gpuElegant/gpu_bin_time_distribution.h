
#ifndef GPU_BIN_TIME_DISTRIBUTION_H
#define GPU_BIN_TIME_DISTRIBUTION_H

#include <gpu_track.h>

#ifdef __cplusplus /* if it's C++ use C linkage */
extern "C" {
#endif
  unsigned int cpu_binParticles_and_countBinned(double* h_histogram, double* h_particles, unsigned int particlePitch, unsigned int np, unsigned int indexToBin, double min, double delta, unsigned int nbins);

  void gpu_binTimeDistribution_and_findMax(double* d_Itime, double* d_time, const unsigned int np,  const double tmin, const double dt, const unsigned int nbins, double* Imax);

  unsigned int gpu_binParticles_and_countBinned(double* d_histogram, double* d_particles, unsigned int particlePitch, unsigned int np, unsigned int indexToBin, double min, double delta, unsigned int nbins);

  unsigned int gpu_binTimeDistribution_and_countBinned(double* d_Itime, double* d_time, const unsigned int np,  const double tmin, const double dt, const unsigned int nbins);

#ifdef __cplusplus
}
#endif
#endif
