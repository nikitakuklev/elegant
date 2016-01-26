#ifndef GPU_CONVOLVE_ARRAYS_H
#define GPU_CONVOLVE_ARRAYS_H

#include <gpu_track.h>

#ifdef __cplusplus /* if it's C++ use C linkage */
extern "C" {
#endif

  void gpuConvolveArrays(double* d_output, unsigned int nOutputs, double* d_a1, unsigned int n1, double* d_a2, unsigned int n2);
  

#ifdef __cplusplus
}
#endif 

#endif 
