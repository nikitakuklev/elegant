#include <iostream>

#include <gpu_bin_time_distribution.h>
#include <gpu_reductions.hcu>
#include <gpu_base.h>
#include <gpu_particle_accessor.hcu>

__inline__ __device__ double ipow(double x, long p)
{
  register double hp;
  register long n;

  if (!x)  {
    return(p==0?1.:0.);
  }

  if (p<0)
    return(1./ipow(x, -p));

  switch (p) {
  case 0: 
    return(1.);
  case 1: 
    return(x);
  case 2: 
    return(x*x);
  case 3: 
    hp = x*x;
    return(hp*x);
  case 4: 
    hp = x*x;
    return(hp*hp);
  case 5: 
    hp = x*x;
    return(hp*hp*x);
  case 6: 
    hp = x*x;
    return(hp*hp*hp);
  case 7: 
    hp = x*x*x;
    return(hp*hp*x);
  case 8: 
    hp = x*x;
    hp = hp*hp;
    return(hp*hp);
  default:
    n = p/2;
    hp = ipow(x, n);
    switch (p-2*n) {
    case 0:
      return(hp*hp);
    case 1:
      return(hp*hp*x);
    }
    break;
  }
  return(0.);  /* never actually executed--keeps compiler happy */
}


inline __device__ double atomicAdd(volatile double* address, double val) {
  double old = *address, assumed;
  do {
    assumed = old;
    old = __longlong_as_double( atomicCAS((unsigned long long int*)address,
					  __double_as_longlong(assumed),
					  __double_as_longlong(val + assumed)));
  } while (assumed != old);
  return old;
} 

#define sharedWarpBlockPitchOffset 0
#define sharedWarpMemoryTarget 12288
#define maxSharedMemoryInBytes 24576
#define nb_decomp_threshold (maxSharedMemoryInBytes / 16 )

#include <gpu_reductions.hcu>

#ifndef sqr
#define sqr(x) (x*x)
#endif

__global__ void gpu_binTransverseKernel(double* d_Itime0, double* d_Itime1,
    const unsigned int nb, double* d_particles, unsigned int particlePitch,
    const unsigned int np, double* d_time, const double tmin, const double dt,
    const double Po, const double dx, const double dy, const int xPower,
    const int yPower, double* d_pz, unsigned int* d_unbinned, 
    const unsigned int nb_pitch, const unsigned int nb_height,
    unsigned int* retirementCount) {
  // shared memory to calculate per-block histograms.  
  // Of size nb_height * nb_pitch, sized to maximize use of shared memory
  // Computing nb_height sub-histograms per thread block to reduce thread 
  // contention in atomics
  extern volatile __shared__ double s_data[]; 
  //[ (sizeof(double))(2* nb_height * nb_pitch) + sizeof(unsigned int)]
  volatile __shared__ double * s_blockHistogram_0; //[nb_height * nb_pitch];
  volatile __shared__ double * s_blockHistogram_1; //[nb_height * nb_pitch];

  // column each thread indexes into for sub-histogram calculation
  unsigned int nb_lane = threadIdx.x % nb_height;
  if(threadIdx.x==0){
    s_blockHistogram_0 = s_data;
    s_blockHistogram_1 = s_data + nb_height * nb_pitch;
  }
  __syncthreads();

  // zero array
  for(unsigned int tx=threadIdx.x; tx<nb_height*nb_pitch; tx+=blockDim.x){
    s_blockHistogram_0[tx] = 0;
    s_blockHistogram_1[tx] = 0;
  }
  __syncthreads();

  double* posItime0 = (double*)(s_blockHistogram_0 + nb_lane*nb_pitch);
  double* posItime1 = (double*)(s_blockHistogram_1 + nb_lane*nb_pitch);

  unsigned int unbinned = 0;// this thread

  // loop over particles, compute histogram, and perform
  // increment(histogram[lane][bin])
  for(unsigned int tid= threadIdx.x + blockIdx.x*blockDim.x; tid < np;
      tid += blockDim.x * gridDim.x) {
    gpuParticleAccessor part(d_particles, tid, particlePitch);
    
    /* Bin CENTERS are at tmin+ib*dt */
    int ib = (d_time[tid]-tmin)/dt+0.5;
    d_pz[tid] = Po*(1+part[5])/sqrt(1+sqr(part[1])+sqr(part[3]));

    if (ib<0){
      unbinned++;
      continue;
    }
    if (ib>nb - 1){
      unbinned++;
      continue;
    }
    if (xPower==1)
      atomicAdd(posItime0+ib,part[0]-dx);
    else if (xPower<1)
      atomicAdd(posItime0+ib, 1);
    else
      atomicAdd(posItime0+ib, ipow(part[0]-dx, xPower));
    if (yPower==1)
      atomicAdd(posItime1+ib, part[2]-dy);
    else if (yPower<1)
      atomicAdd(posItime1+ib, 1);
    else
      atomicAdd(posItime1+ib, ipow(part[2]-dy, yPower));
  }
  
  // synchronize shared memory
  __syncthreads();

  // combine sub-histograms in shared memory
  for(unsigned int w = 1; w < nb_height; w++){    
    for(unsigned int tx=threadIdx.x; tx < nb; tx+= blockDim.x){
      s_blockHistogram_0[tx] += s_blockHistogram_0[w * nb_pitch + tx];
      s_blockHistogram_1[tx] += s_blockHistogram_1[w * nb_pitch + tx];
    }
  }

  __syncthreads();

  // atomic add to global memory to calculate final histogram
  for ( uint tx = threadIdx.x; tx < nb; tx+=blockDim.x){
    double r_temp_0 = s_blockHistogram_0[tx];
    double r_temp_1 = s_blockHistogram_1[tx];
    
    unsigned int ib = tx + blockIdx.y * nb_pitch;
    if(r_temp_0)
      atomicAdd( d_Itime0+ib, r_temp_0);    
    if(r_temp_1)
      atomicAdd( d_Itime1+ib, r_temp_1);    
    
  }
  
  __syncthreads();

  s_blockHistogram_0[threadIdx.x] = unbinned;
    
  __syncthreads();

  double unbinnedf = unbinned;
  reduceBlock<double, 2048, Add<double> >
      (s_blockHistogram_0, unbinnedf, threadIdx.x, blockDim.x, Add<double>() );
    
  if(threadIdx.x==0 && s_blockHistogram_0[0]){
    atomicAdd(d_unbinned, (unsigned int)s_blockHistogram_0[0]);
  }
}

template<unsigned int NB_PER_BLOCK>
__global__ void gpu_binTransverseKernel_binDecomp(double* d_Itime0,
    double* d_Itime1, const unsigned int nb, double* d_particles,
    unsigned int particlePitch, const unsigned int np, double* d_time,
    const double tmin, const double dt, const double Po, const double dx,
    const double dy, const int xPower, const int yPower, double* d_pz,
    unsigned int* d_unbinned, unsigned int* retirementCount) {
  // shared memory to calculate per-block sub-histograms.  
  volatile __shared__ double s_blockHistogram_0[NB_PER_BLOCK];
  volatile __shared__ double s_blockHistogram_1[NB_PER_BLOCK];

  // zero array
  for(unsigned int tx=threadIdx.x; tx<NB_PER_BLOCK; tx+=blockDim.x){
    s_blockHistogram_0[tx] = 0;
    s_blockHistogram_1[tx] = 0;
  }
  __syncthreads();

  int ib_min = blockIdx.y * NB_PER_BLOCK;
  int ib_max = (blockIdx.y + 1) * NB_PER_BLOCK;
 
  unsigned int unbinned = 0;// this thread

  // loop over particles, compute histogram, and perform
  // increment(histogram[lane][bin])
  for(unsigned int tid= threadIdx.x + blockIdx.x*blockDim.x; tid < np;
      tid += blockDim.x * gridDim.x){
    gpuParticleAccessor part(d_particles, tid, particlePitch);
    
    /* Bin CENTERS are at tmin+ib*dt */
    int ib = (d_time[tid]-tmin)/dt+0.5;
    if(blockIdx.y == 0)
      d_pz[tid] = Po*(1+part[5])/sqrt(1+sqr(part[1])+sqr(part[3]));

    if(ib >= ib_min && ib < ib_max){
      int local_ib = ib - ib_min;
      
      if (ib<0){
	unbinned++;
	continue;
      }
      if (ib>nb - 1){
	unbinned++;
	continue;
      }
      if (xPower==1)
	atomicAdd(s_blockHistogram_0+local_ib,part[0]-dx);
      else if (xPower<1)
	atomicAdd(s_blockHistogram_0+local_ib, 1);
      else
	atomicAdd(s_blockHistogram_0+local_ib, ipow(part[0]-dx, xPower));
      if (yPower==1)
	atomicAdd(s_blockHistogram_1+local_ib, part[2]-dy);
      else if (yPower<1)
	atomicAdd(s_blockHistogram_1+local_ib, 1);
      else
	atomicAdd(s_blockHistogram_1+local_ib, ipow(part[2]-dy, yPower));
    }
  }
  // synchronize shared memory
  __syncthreads();
 
  // atomic add to global memory to calculate final histogram
  for ( uint tx = threadIdx.x; tx < NB_PER_BLOCK; tx+=blockDim.x){
    unsigned int ib = tx + ib_min;
    if(ib < nb){
      double r_temp_0 = s_blockHistogram_0[tx];
      double r_temp_1 = s_blockHistogram_1[tx];

      if(r_temp_0)
	atomicAdd( d_Itime0+ib, r_temp_0);    
      if(r_temp_1)
	atomicAdd( d_Itime1+ib, r_temp_1);    
    }
  }
  
  __syncthreads();

  s_blockHistogram_0[threadIdx.x] = unbinned;
    
  __syncthreads();

  double unbinnedf = unbinned;
  reduceBlock<double, NB_PER_BLOCK, Add<double> >
      (s_blockHistogram_0, unbinnedf, threadIdx.x, blockDim.x, Add<double>() );
    
  if(threadIdx.x==0 && s_blockHistogram_0[0]) {
    atomicAdd(d_unbinned, (unsigned int)s_blockHistogram_0[0]);
  }
}



static bool isSet_binTransverse = false;
void setCacheConfig_binTransverse(){
  if(isSet_binTransverse == false){
    isSet_binTransverse = true;    
    cudaFuncSetCacheConfig(gpu_binTransverseKernel,cudaFuncCachePreferShared);
    cudaFuncSetCacheConfig(gpu_binTransverseKernel_binDecomp<nb_decomp_threshold>,
                           cudaFuncCachePreferShared);
  }
}

extern "C" {

unsigned int gpu_binTransverseDistribution(double* d_Itime0,
    double* d_Itime1, const unsigned int& nb, double* d_particles,
    unsigned int particlePitch, const unsigned int& np, double* d_time,
    const double& tmin, const double& dt, const double& Po, const double& dx,
    const double& dy, const int& xPower, const int& yPower, double* d_pz) {
  if(!np) return 0;
  
  setCacheConfig_binTransverse();

  GPUBASE* gpuBase = getGpuBase();

  unsigned int* d_unbinned = gpuBase->d_tempu_alpha;
  unsigned int* d_retirementCount = gpuBase->d_retirementCount;  

  cudaMemset(d_Itime0,0,sizeof(double)*nb);
  cudaMemset(d_Itime1,0,sizeof(double)*nb);
  cudaMemset(d_unbinned, 0, sizeof(unsigned int) );
  cudaMemset(d_retirementCount, 0, sizeof(unsigned int));
  
  if(nb > nb_decomp_threshold){
    unsigned int nThreads = gpuBase->nReductionThreads;
  
    unsigned int nBx = (np + nThreads - 1)/nThreads;
    if(nBx > gpuBase->nReductionBlocks) nBx = gpuBase->nReductionBlocks;
    unsigned int nBy = (nb + nb_decomp_threshold - 1)/nb_decomp_threshold;

    dim3 gridDim(nBx,nBy,1);
    dim3 blockDim(nThreads,1,1);
    gpu_binTransverseKernel_binDecomp<nb_decomp_threshold>
      <<<gridDim,blockDim>>>(d_Itime0, d_Itime1, nb, d_particles,
          particlePitch, np, d_time, tmin, dt, Po, dx, dy, xPower, yPower,
          d_pz, d_unbinned, d_retirementCount);
    gpuErrorHandler("gpu_binTransverseDistribution::gpu_binTransverseKernel_binDecomp");
  } else {
    // find the nearest value of number of bins to 32 (rounded up), and
    // thus the pitch of the shared memory arrays
    unsigned int nb_rounded = (nb + 32 - 1) / 32;
    nb_rounded = 32 * nb_rounded;
  
    unsigned int nb_pitch = nb_rounded + sharedWarpBlockPitchOffset;
    unsigned int nThreads = gpuBase->nReductionThreads;
    unsigned int sharedMemoryInBytes = sharedWarpMemoryTarget; // bytes
    unsigned int nb_height = sharedMemoryInBytes / nb_pitch / (2*sizeof(double));

    // there is an oddly large dependence on making nb_height an odd number.  
    if(nb_height % 2 == 0){
      nb_height = nb_height + 1;
    }
    sharedMemoryInBytes = 2* nb_pitch * nb_height * sizeof(double);

    unsigned int nBlocks = (np + nThreads -1 )/nThreads;
    if(nBlocks > gpuBase->nReductionBlocks) nBlocks = gpuBase->nReductionBlocks;
    dim3 dimGrid(nBlocks, 1, 1); 
    dim3 dimBlock(nThreads, 1, 1);

    gpu_binTransverseKernel<<<dimGrid, dimBlock, sharedMemoryInBytes>>>
      ( d_Itime0,  d_Itime1,  nb,  d_particles,  particlePitch,  np,
        d_time,  tmin,  dt,  Po,  dx,  dy, xPower, yPower,  d_pz,
        d_unbinned, nb_pitch, nb_height, d_retirementCount);
    gpuErrorHandler("gpu_binTransverseDistribution::gpu_binTransverseKernel");
  }

  unsigned int unbinned= 0;
  cudaMemcpy(&unbinned, d_unbinned, sizeof(unsigned int), cudaMemcpyDeviceToHost);
  unsigned int nbinned = np - unbinned;
  return nbinned;
}

} // extern "C"
