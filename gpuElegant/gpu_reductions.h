#ifndef GPU_REDUCTIONS_H
#define GPU_REDUCTIONS_H

#include <gpu_track.h>

#ifdef __CUDACC__
template<class T, class OPERATOR>
void gpuReduction(T* d_a, unsigned int n, OPERATOR op);
#endif

#ifdef __cplusplus /* if it's C++ use C linkage */
extern "C" {
#endif

void gpuFindMinIndex(double* d_a, unsigned int n, double* minValue, 
    unsigned int* indexOfMinValue);

void gpuFindMaxIndex(double* d_a, unsigned int n, double* maxValue, 
    unsigned int* indexOfMaxValue);

double gpuReduceAdd(double* d_a, unsigned int n);

void gpuReduceAddAsync(double* d_a, unsigned int n, double* sum);

long gpuReduceAddLong(long* d_a, unsigned int n, long* d_templ=0);

double gpuReduceMin(double* d_a, unsigned int n);

double gpuReduceMax(double* d_a, unsigned int n);

void gpuReduceMinMax(double* d_a, unsigned int n, double* min, double* max);

void gpuReduceMinMaxAsync(double* d_a, unsigned int n, double* min, double* max);

void gpuReduceUMinMax(unsigned int* d_a, unsigned int n, unsigned int* min, 
                     unsigned int* max);

double thrustReduceAdd(double* d_a, unsigned int n); 

#ifdef __cplusplus
}
#endif 

#endif // GPU_REDUCTIONS_H
