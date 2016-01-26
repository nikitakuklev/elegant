#include <gpu_particle_template_function.hcu>
#include <gpu_particle_reduction.hcu>

#define USE_SWAP
#ifdef USE_SWAP
/**
 * Setup a mapping to killed particle indicies from d_linearIndex
 * which is the result of an inclusive scan of a killed particle boolean
 * @param d_map (output) the map array
 * @param d_linearIndex (input) the scanned killed particle boolean array
 * @param np number of particles
 */
__global__ void setMap(unsigned int* d_map, unsigned int* d_linearIndex, 
                       unsigned int np) {
  unsigned int ind;
  extern volatile __shared__ unsigned int s_linInd[];

  for(ind = blockDim.x * blockIdx.x + threadIdx.x; ind < np; 
      ind += blockDim.x * gridDim.x) {
    s_linInd[threadIdx.x+1] = d_linearIndex[ind];
    if (threadIdx.x == 0) {
      if (ind > 0)
        s_linInd[0] = d_linearIndex[ind-1];
      else 
        s_linInd[0] = 0;
    }

    __syncthreads();

    if (s_linInd[threadIdx.x+1]-s_linInd[threadIdx.x])
      d_map[s_linInd[threadIdx.x]] = ind;
  }
}

/**
 * swap killed particles to the end of the array
 * @param d_part device particles array
 * @param d_map map to killed particle indicies
 * @param particlePitch particle SoA array pitch
 * @param nLeft number of particles remaining
 * @param nKilled number of particles lost
 */
__global__ void swapKilledInPlace(double* d_part, unsigned int* d_map,
    unsigned int particlePitch, unsigned int nLeft, unsigned int nKilled) {
  unsigned int idx, ikilled;
  double new_part[7];

  for(idx = blockDim.x*blockIdx.x+threadIdx.x; idx < nKilled; 
      idx += blockDim.x*gridDim.x) {
    if (d_map[idx]>=nLeft) return;
    ikilled = nLeft + idx;
    #pragma unroll 7
    for(unsigned int ii=0; ii<7; ii++) {
      new_part[ii]=d_part[ikilled+ii*particlePitch];
      d_part[ikilled+ii*particlePitch]=d_part[d_map[idx]+ii*particlePitch];
      d_part[d_map[idx]+ii*particlePitch]=new_part[ii];
    }
  }
}
#endif /* USE_SWAP */
