#ifndef GPU_TRWAKE_H
#define GPU_TRWAKE_H

#include <gpu_track.h>
#include <gpu_base.h>

#ifdef __cplusplus
extern "C" {
#endif

void gpu_track_through_trwake(long np, TRWAKE *wakeData, double Po,
			      RUN *run, long i_pass, CHARGE *charge);

#ifdef __cplusplus
}
#endif

void gpu_computeTimeCoordinates(unsigned int np, double* d_time, double Po);

void gpu_computeTimeCoordinatesAndMinMax(unsigned int np, double* d_time,
       double Po, double* tmin, double* tmax);

double gpu_computeTimeCoordinatesAndMean(long np, double* d_time, 
                                         double Po);

void gpu_computeDistanceCoordinates(double *d_time, double Po, long np);

#endif /* GPU_TRWAKE_H */
