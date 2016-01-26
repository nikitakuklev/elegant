#ifndef GPU_KAHAN_H
#define GPU_KAHAN_H

#include <gpu_track.h>

#ifdef __cplusplus /* if it's C++ use C linkage */
extern "C" {
#endif

void gpuKahanAsync(double* d_a, double* sum, double* error, long n);

double gpuKahan(double* d_a, double* error, long n);

double gpu_P_average_kahanSum(double* d_coord5, double* P_central, 
                double* error, long n);

#ifdef __cplusplus
}
#endif
#endif
