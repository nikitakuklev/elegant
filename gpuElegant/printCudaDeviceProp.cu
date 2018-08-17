#include <printCudaDeviceProp.h>
#include <cuda.h>
#include <stdio.h>

void printCudaDeviceProp()
{
  int n = 10;

  double *d_v;
  cudaMalloc((void **)&d_v, n * sizeof(double));
  cudaMemset(d_v, 0, n * sizeof(double));
  cudaFree(d_v);
  cudaDeviceProp prop;

  cudaGetDeviceProperties(&prop, 0);

  printf("CUDA DEVICE:     name:   %s\n", prop.name);
  printf("CUDA DEVICE:   mem mb:   %g\n", (prop.totalGlobalMem) / 1024. / 1024.);
  printf("CUDA DEVICE:    major:   %d\n", prop.major);
  printf("CUDA DEVICE:    minor:   %d\n", prop.minor);
  printf("CUDA DEVICE:  MPcount:   %d\n", prop.multiProcessorCount);
  printf("CUDA DEVICE:   integrated:   %d\n", prop.integrated);
  printf("CUDA DEVICE:   warpSize:   %d\n", prop.warpSize);
}
