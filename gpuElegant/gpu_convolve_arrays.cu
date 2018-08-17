#include <iostream>
#include <stdio.h>

#include <gpu_convolve_arrays.h>

#ifdef max
#  undef max
#endif
#ifdef min
#  undef min
#endif

inline __device__ double atomicAddOAG(double *address, double val)
{

  double old = *address, assumed;
  do
    {
      assumed = old;
      old = __longlong_as_double(atomicCAS((unsigned long long int *)address,
                                           __double_as_longlong(assumed),
                                           __double_as_longlong(val + assumed)));
    }
  while (assumed != old);
  return old;
}

template <int BUFFER_X, int BUFFER_Y, int NTHREADS>
  __global__ void gpuConvolveArraysKernel(double *d_output,
                                          unsigned int nOutputs, double *d_a1, unsigned int n1, double *d_a2,
                                          unsigned int n2)
{

  __shared__ double s_a2_buffer[BUFFER_X];

  __shared__ double s_sums[BUFFER_Y][NTHREADS];

  int blockRow = blockIdx.y;
  int blockCol = blockIdx.x;

  if (blockCol * BUFFER_X > blockRow * BUFFER_Y)
    return;

  // row index is the result index
  int ibStart = blockRow * BUFFER_Y;

  // col index is what we're buffering for a2
  int ib2Start = blockCol * BUFFER_X;

  for (int tx = threadIdx.x; tx < BUFFER_X; tx += NTHREADS)
    {
      if (ib2Start + tx < n2)
        s_a2_buffer[tx] = d_a2[ib2Start + tx];
      else
        s_a2_buffer[tx] = 0;
    }
  __syncthreads();

  for (int sumIndex = 0; sumIndex < BUFFER_Y; sumIndex++)
    {
      int ib = ibStart + sumIndex;

      s_sums[sumIndex][threadIdx.x] = 0;

      // There's a bug on the cpu side when n1 < nOutputs
      // the CPU code reads past n1 during the convolution.
      if (ib < nOutputs)
        {
          for (int tx = threadIdx.x; tx < BUFFER_X; tx += NTHREADS)
            {
              int ib2 = ib2Start + tx;
              int ib1 = ib - ib2;
              if (ib2 <= ib && ib1 < n1)
                {
                  double r_a1 = d_a1[ib1];
                  double r_a2 = s_a2_buffer[ib2 - ib2Start];
                  s_sums[sumIndex][threadIdx.x] += r_a1 * r_a2;
                }
            }
        }
    }
  __syncthreads();

  unsigned int blkStep = NTHREADS;
  blkStep >>= 1;
  while (blkStep >= 1)
    {
      if (threadIdx.x < blkStep && threadIdx.x + blkStep < NTHREADS)
        for (int ib = 0; ib < BUFFER_Y; ib++)
          s_sums[ib][threadIdx.x] =
            s_sums[ib][threadIdx.x] + s_sums[ib][threadIdx.x + blkStep];
      __syncthreads();
      blkStep >>= 1;
    }

  for (unsigned int tx = threadIdx.x; tx < BUFFER_Y; tx += NTHREADS)
    {
      unsigned int ib = ibStart + tx;
      if (ib < nOutputs)
        {
          atomicAddOAG((double *)(d_output + ib), s_sums[tx][0]);
        }
    }
}

void gpuConvolveArrays(double *d_output, unsigned int nOutputs, double *d_a1,
                       unsigned int n1, double *d_a2, unsigned int n2)
{

  cudaMemset(d_output, 0, nOutputs * sizeof(double));

  const unsigned int nTx = 64;
  const unsigned int BUFFER_X = 512;
  const unsigned int BUFFER_Y = 8;

  int buffered_x_blocks = (nOutputs + BUFFER_X - 1) / BUFFER_X;
  int buffered_y_blocks = (nOutputs + BUFFER_Y - 1) / BUFFER_Y;
  dim3 blocks(buffered_x_blocks, buffered_y_blocks, 1);

  gpuConvolveArraysKernel<BUFFER_X, BUFFER_Y, nTx><<<blocks, nTx>>>(d_output, nOutputs, d_a1, n1, d_a2, n2);
}
