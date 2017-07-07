#ifndef GPU_LIMIT_AMPLITUDES_H
#define GPU_LIMIT_AMPLITUDES_H

#include <gpu_base.h>
#include <gpu_track.h>

#ifdef __cplusplus /* if it's C++ use C linkage */
extern "C" {
#endif

long gpu_rectangular_collimator(RCOL *rcol, long np, double **accepted, double z,
                                double Po);

long gpu_limit_amplitudes(
    double xmax, double ymax, long np, double **accepted,
    double z, double Po, long extrapolate_z, long openCode);

long gpu_removeInvalidParticles(long np, double **accepted, double z, double Po);

long gpu_elliptical_collimator(ECOL *ecol, long np, double **accepted,
                               double z, double Po);

long gpu_beam_scraper(SCRAPER *scraper, long np, double **accepted,
                      double z, double Po);

long gpu_imposeApertureData(long np, double **accepted,
                            double z, double Po, APERTURE_DATA *apData);

long gpu_elimit_amplitudes(double xmax, double ymax, long np, 
                           double **accepted, double z, double Po, long extrapolate_z, 
                           long openCode, long exponent, long yexponent);

#ifdef __cplusplus
}
#endif
#endif

