#include <iostream>

#include <gpu_bin_time_distribution.h>
#include <gpu_reductions.hcu>
#include <gpu_base.h>

inline __device__ double atomicAddOAG(double* address, double val) {
  double old = *address, assumed;
  do {
    assumed = old;
    old = __longlong_as_double( atomicCAS((unsigned long long int*)address,
					  __double_as_longlong(assumed),
					  __double_as_longlong(val + assumed)));
  } while (assumed != old);
  return old;
} 

class gpuBinParticleDistribution {
public:
  gpuBinParticleDistribution(double* d_particles, unsigned int particlePitch, 
      unsigned int index, double min, double delta){
    this->d_particles = d_particles;
    this->particlePitch = particlePitch;
    this->index = index;
    this->min = min;
    this->delta = delta;
  }
  inline __device__ bool calculateBin(unsigned int& ib, unsigned int tid, 
      const unsigned int& nb){
    ib = (d_particles[particlePitch*index + tid] - min)/delta;
    return (ib < nb);
  }
  double* d_particles;
  unsigned int particlePitch;
  unsigned int index;
  double min;
  double delta;
};

class gpuBinTimeDistribution {
public:
  gpuBinTimeDistribution(double* d_particleTime, double tmin, double dt){
    this->d_particleTime = d_particleTime;
    this->tmin = tmin;
    this->dt = dt;
  }
  inline __device__ bool calculateBin(unsigned int& ib, unsigned int tid, 
      const unsigned int& nb){
    ib = (d_particleTime[tid] - tmin)/dt + 0.5;
    return (ib < nb);
  }
  double* d_particleTime;
  double tmin;
  double dt;
};

class reduceNothing {
public:
  reduceNothing() { }
  inline __device__ __host__ bool doReduction(){ return false; }
  inline __device__ __host__ unsigned int operation(unsigned int r1, 
      unsigned int r2) {  return r1; }

  inline __device__ __host__ void reduce(volatile unsigned int* s_data, 
      unsigned int r1, const unsigned int& tx, const unsigned int& blockSize){
  }
};

template<unsigned int NTHREADS>
class reduceAdd {
public:
  reduceAdd(unsigned int* d_nBinned){
    this->d_nBinned = d_nBinned;
  }
  inline __device__ bool doReduction(){ return true; }

  inline __device__ unsigned int operation(unsigned int r1, unsigned int r2){
    return r1 + r2;
  }
  inline __device__ void reduce(volatile unsigned int* s_data, unsigned int r1,
      const unsigned int& tx, const unsigned int& blockSize){
    reduceBlock<unsigned int, NTHREADS, Add<unsigned int> >
        (s_data, r1, tx, blockSize, Add<unsigned int>());
    if(tx==0){
      *d_nBinned = s_data[0];
    }
  }
  unsigned int* d_nBinned;
};

template<unsigned int NTHREADS>
class reduceMax {
public:
  reduceMax(unsigned* d_globalMax){ 
    this->d_globalMax = d_globalMax;
  }
  inline __device__ bool doReduction(){ return true; }

  inline __device__ unsigned int operation(unsigned int r1, unsigned int r2){
    return max(r1, r2);
  }
  inline __device__ void reduce(volatile unsigned int* s_data, unsigned int r1, 
      const unsigned int& tx, const unsigned int& blockSize){
    reduceBlock<unsigned int, NTHREADS, Max<unsigned int> >
        (s_data, r1, tx, blockSize, Max<unsigned int>());
    if(tx == 0)
      *d_globalMax = s_data[0];
  }
  unsigned* d_globalMax;
};

#define sharedWarpBlockPitchOffset 0
#define sharedWarpMemoryTarget 12288
#define maxSharedMemoryInBytes 24576
#define nb_decomp_threshold (maxSharedMemoryInBytes / 8 )

template<class BIN_FUNCTION, class REDUCTION_OPERATION>
__global__ void gpu_binDistribution_kernel(double* d_doubleHistogram, 
    unsigned int* d_unsignedHistogram, const unsigned int nb, 
    const unsigned int nb_pitch, const unsigned int nb_height,
    const unsigned int np, unsigned int* retirementCount, 
    BIN_FUNCTION binFunction, REDUCTION_OPERATION reductionOperation){
  // shared memory to calculate per-block histograms.  
  // Of size nb_height * nb_pitch, sized to maximize use of shared memory
  // Computing nb_height sub-histograms per thread block to reduce thread 
  // contention in atomics
  extern volatile __shared__ unsigned int s_blockHistogram[]; //[nb_height * nb_pitch];
  // column each thread indexes into for sub-histogram calculation
  unsigned int nb_lane = threadIdx.x % nb_height;
  
  // zero array
  for(unsigned int tx=threadIdx.x; tx<nb_height*nb_pitch; tx+=blockDim.x){
    s_blockHistogram[tx] = 0;
  }
  __syncthreads();
  
  // loop over particles, compute histogram, and perform 
  // increment(histogram[lane][bin])
  for(unsigned int tid= threadIdx.x + blockIdx.x*blockDim.x; tid < np; 
      tid += blockDim.x * gridDim.x){
    unsigned int ib;
    bool validBin = binFunction.calculateBin(ib, tid, nb);
    if(validBin){
      atomicInc((unsigned int*)(s_blockHistogram + nb_lane*nb_pitch + ib),np);
    }
  }
  
  // synchronize shared memory
  __syncthreads();

  // combine sub-histograms in shared memory
  for(unsigned int w = 1; w < nb_height; w++){    
    for(unsigned int tx=threadIdx.x; tx < nb; tx+= blockDim.x){
      s_blockHistogram[tx] += s_blockHistogram[w * nb_pitch + tx];
    }
  }

  __syncthreads();

  // atomic add to global memory to calculate final histogram
  for ( uint tx = threadIdx.x; tx < nb; tx+=blockDim.x){
    unsigned int r_temp = s_blockHistogram[tx];

       if(r_temp){
      unsigned int oldt = atomicAddOAG((double*)(d_unsignedHistogram + tx), r_temp);    
    }
  }
  
  // code taken and modified from the threadFenceReduction CUDA SDK example
  if(reductionOperation.doReduction() ){
    __shared__ volatile bool amLast;
    __threadfence();
    // each block takes a ticket, the last block to finish the calculation 
    // does the reduction
    if(threadIdx.x==0){
      uint ticket = atomicInc(retirementCount, gridDim.x);
      amLast = (ticket == gridDim.x - 1);
    }
    __syncthreads();    

    if( amLast ){

      uint r_op = 0;

      for( uint tx = threadIdx.x; tx < nb; tx += blockDim.x){
	unsigned int temp = d_unsignedHistogram[tx] ;

	d_doubleHistogram[tx] = temp;
      	r_op = reductionOperation.operation(r_op, temp);
      }
      
      s_blockHistogram[threadIdx.x] = r_op;
      // find the maximum or summed value of histogram
      reductionOperation.reduce(s_blockHistogram, r_op, threadIdx.x,
                                blockDim.x);
    } //amLast
  } // gridDim.x > 1
}

unsigned int cpu_binParticles_and_countBinned(double* h_histogram, 
    double* h_particles, unsigned int particlePitch, unsigned int np,
    unsigned int indexToBin, double min, double delta, unsigned int nbins){

  unsigned int nBinned = 0;
  memset(h_histogram, 0, nbins*sizeof(double) );
  for(unsigned int i = 0; i< np; i++){
    unsigned int ib = (h_particles[i + particlePitch*indexToBin] - min)/delta;
    if(ib < nbins){
      h_histogram[ib]++;
      nBinned++;
    }
  }
  return nBinned;
}
      
#include <gpu_reductions.h>

template<unsigned int NB_PER_BLOCK, class BIN_FUNCTION, class REDUCTION_OPERATION>
__global__ void gpu_binDistribution_decompKernel(double* d_doubleHistogram, 
    unsigned int* d_unsignedHistogram, const unsigned int nb, 
    const unsigned int np, unsigned int* retirementCount, 
    BIN_FUNCTION binFunction, REDUCTION_OPERATION reductionOperation){
  // shared memory to calculate per-block sub-histograms.  
  volatile __shared__ unsigned int s_blockHistogram[NB_PER_BLOCK];

  // zero array
  for(unsigned int tx=threadIdx.x; tx<NB_PER_BLOCK; tx+=blockDim.x){
    s_blockHistogram[tx] = 0;
  }
  __syncthreads();

  int ib_min = blockIdx.y * NB_PER_BLOCK;
  int ib_max = (blockIdx.y + 1) * NB_PER_BLOCK;
 
  // loop over particles, compute histogram, and perform
  // increment(histogram[lane][bin])
  for(unsigned int tid= threadIdx.x + blockIdx.x*blockDim.x; tid < np;
      tid += blockDim.x * gridDim.x){
    unsigned int ib;
    bool validBin = binFunction.calculateBin(ib, tid, nb);
    if(validBin){
      if(ib >= ib_min && ib < ib_max){
	int local_ib = ib - ib_min;
	atomicInc((unsigned int*)(s_blockHistogram + local_ib),np);
      }
    }
  }
  // synchronize shared memory
  __syncthreads();
 
  // atomic add to global memory to calculate final histogram
  for ( uint tx = threadIdx.x; tx < NB_PER_BLOCK; tx+=blockDim.x){
    unsigned int ib = tx + ib_min;
    if(ib < nb){
      unsigned int r_temp = s_blockHistogram[tx];
      if(r_temp)
	atomicAddOAG((double*)(d_unsignedHistogram + ib), r_temp);    
    }
  }
  
  __syncthreads();
  
  // code taken and modified from the threadFenceReduction CUDA SDK example
  if(reductionOperation.doReduction() ){
    __shared__ volatile bool amLast;
    __threadfence();
    // each block takes a ticket, the last block to finish the calculation
    // does the reduction
    if(threadIdx.x==0){
      uint ticket = atomicInc(retirementCount, gridDim.x*gridDim.y);
      amLast = (ticket == gridDim.x*gridDim.y - 1);
    }
    __syncthreads();    

    if( amLast ){

      uint r_op = 0;

      for( uint tx = threadIdx.x; tx < nb; tx += blockDim.x){
	unsigned int temp = d_unsignedHistogram[tx] ;

	d_doubleHistogram[tx] = temp;
      	r_op = reductionOperation.operation(r_op, temp);
      }
      
      s_blockHistogram[threadIdx.x] = r_op;
      // find the maximum or summed value of histogram
      reductionOperation.reduce(s_blockHistogram, r_op, threadIdx.x,
                                blockDim.x);
    } //amLast
  }
}

/*
static bool isSet_binParticles_reduceAdd = false;
void setCacheConfig_binParticles_reduceAdd(){
  if(isSet_binParticles_reduceAdd == false){
    isSet_binParticles_reduceAdd = true;  
    cudaFuncSetCacheConfig(gpu_binDistribution_kernel<gpuBinParticleDistribution,
        reduceAdd<256>>, cudaFuncCachePreferShared);  
    cudaFuncSetCacheConfig(gpu_binDistribution_kernel<gpuBinParticleDistribution,
        reduceAdd<512>>, cudaFuncCachePreferShared);  
    cudaFuncSetCacheConfig(gpu_binDistribution_kernel<gpuBinParticleDistribution,
        reduceAdd<1024>>, cudaFuncCachePreferShared);  
    cudaFuncSetCacheConfig(gpu_binDistribution_kernel<gpuBinParticleDistribution,
        reduceAdd<2048>>, cudaFuncCachePreferShared);  
    cudaFuncSetCacheConfig(
        gpu_binDistribution_decompKernel<nb_decomp_threshold, 
        gpuBinParticleDistribution, reduceAdd>,cudaFuncCachePreferShared);
  }
}

static bool isSet_binTime_reduceAdd = false;
void setCacheConfig_binTime_reduceAdd(){
  if(isSet_binTime_reduceAdd == false){
    isSet_binTime_reduceAdd = true;    
    cudaFuncSetCacheConfig(gpu_binDistribution_kernel<gpuBinParticleDistribution,
        reduceAdd<256>>, cudaFuncCachePreferShared);  
    cudaFuncSetCacheConfig(gpu_binDistribution_kernel<gpuBinParticleDistribution,
        reduceAdd<512>>, cudaFuncCachePreferShared);  
    cudaFuncSetCacheConfig(gpu_binDistribution_kernel<gpuBinParticleDistribution,
        reduceAdd<1024>>, cudaFuncCachePreferShared);  
    cudaFuncSetCacheConfig(gpu_binDistribution_kernel<gpuBinParticleDistribution,
        reduceAdd<2048>>, cudaFuncCachePreferShared);  
    cudaFuncSetCacheConfig(
        gpu_binDistribution_decompKernel<nb_decomp_threshold,
        gpuBinTimeDistribution, reduceAdd>,cudaFuncCachePreferShared);
  }
}

static bool isSet_binTime_reduceMax= false;
void setCacheConfig_binTime_reduceMax(){
  if(isSet_binTime_reduceMax== false){
    isSet_binTime_reduceMax = true;    
    cudaFuncSetCacheConfig(gpu_binDistribution_kernel<gpuBinParticleDistribution,
        reduceMax<256>>, cudaFuncCachePreferShared);  
    cudaFuncSetCacheConfig(gpu_binDistribution_kernel<gpuBinParticleDistribution,
        reduceMax<512>>, cudaFuncCachePreferShared);  
    cudaFuncSetCacheConfig(gpu_binDistribution_kernel<gpuBinParticleDistribution,
        reduceMax<1024>>, cudaFuncCachePreferShared);  
    cudaFuncSetCacheConfig(gpu_binDistribution_kernel<gpuBinParticleDistribution,
        reduceMax<2048>>, cudaFuncCachePreferShared);  
    cudaFuncSetCacheConfig(
        gpu_binDistribution_decompKernel<nb_decomp_threshold,
        gpuBinTimeDistribution, reduceMax>,cudaFuncCachePreferShared);
  }
}
*/

template<class TRANSFORM, class REDUCTION>
void gpuBinDriver(double* d_histogram, unsigned int nbins, unsigned int np,
    TRANSFORM& transform, REDUCTION& reduction){
  if (!np) return;
 
  GPUBASE* gpuBase = getGpuBase();
 
  unsigned int* d_tempu = gpuBase->d_tempu_alpha;
  unsigned int* d_retirementCount = gpuBase->d_retirementCount;  
 
  cudaMemset(d_retirementCount, 0, sizeof(unsigned int));
   
  // find the nearest value of number of bins to 32 (rounded up), and thus 
  // the pitch of the shared memory arrays
  unsigned int nb_rounded = (nbins + 32 - 1) / 32;
  nb_rounded = 32 * nb_rounded;
  cudaMemset(d_tempu, 0, sizeof(unsigned int)*nb_rounded);
 
  if(nbins > nb_decomp_threshold){    
    unsigned int nThreads = gpuBase->nReductionThreads;
     
    unsigned int nBx = (np + nThreads - 1)/nThreads;
    if(nBx > gpuBase->nReductionBlocks) nBx = gpuBase->nReductionBlocks;
    unsigned int nBy = (nbins + nb_decomp_threshold - 1)/nb_decomp_threshold;
 
    dim3 gridDim(nBx,nBy,1);
    dim3 blockDim(nThreads,1,1);
    gpu_binDistribution_decompKernel
      <nb_decomp_threshold, TRANSFORM, REDUCTION><<<gridDim,blockDim>>>
      (d_histogram, d_tempu, nbins, np, d_retirementCount, transform, reduction);

  } else {
   
    unsigned int nb_pitch = nb_rounded + sharedWarpBlockPitchOffset;
    unsigned int nThreads = gpuBase->nReductionThreads;
    unsigned int sharedMemoryInBytes = sharedWarpMemoryTarget; // bytes
    unsigned int nb_height = sharedMemoryInBytes / nb_pitch / sizeof(unsigned int);
 
    // there is an oddly large dependence on making nb_height an odd number.  
    if(nb_height % 2 == 0)
      nb_height = nb_height + 1;
 
    sharedMemoryInBytes = nb_pitch * nb_height * sizeof(unsigned int);
 
    unsigned int nBlocks = (np + nThreads -1 )/nThreads;
    if(nBlocks > gpuBase->nReductionBlocks) nBlocks = gpuBase->nReductionBlocks;
    dim3 dimGrid(nBlocks, 1, 1); 
    dim3 dimBlock(nThreads, 1, 1);

    gpu_binDistribution_kernel<TRANSFORM, REDUCTION>
      <<<dimGrid, dimBlock, sharedMemoryInBytes>>>(d_histogram, d_tempu, nbins,
          nb_pitch, nb_height, np, d_retirementCount, transform, reduction);
  }
}

unsigned int gpu_binParticles_and_countBinned(double* d_histogram,
    double* d_particles, unsigned int particlePitch, unsigned int np,
    unsigned int indexToBin, double min, double delta, unsigned int nbins) {
  GPUBASE* gpuBase = getGpuBase();
  unsigned int* d_tempu = gpuBase->d_blockTempu;

  gpuBinParticleDistribution transform(d_particles, particlePitch, indexToBin, min, delta); 
  unsigned int nTx = gpuBase->nReductionThreads; 
  if(nTx==256) {
    reduceAdd<256> reduction(d_tempu);
    gpuBinDriver(d_histogram, nbins, np, transform, reduction);
  }
  else if(nTx==512) {
    reduceAdd<512> reduction(d_tempu);
    gpuBinDriver(d_histogram, nbins, np, transform, reduction);
  }
  else if(nTx==1024) {
    reduceAdd<1024> reduction(d_tempu);
    gpuBinDriver(d_histogram, nbins, np, transform, reduction);
  }
  else if(nTx==2048) {
    reduceAdd<2048> reduction(d_tempu);
    gpuBinDriver(d_histogram, nbins, np, transform, reduction);
  }
  //setCacheConfig_binParticles_reduceAdd();

  //gpuBinDriver(d_histogram, nbins, np, transform, reduction);

  unsigned int nBinned;
  cudaMemcpy(&nBinned, d_tempu, sizeof(unsigned int), cudaMemcpyDeviceToHost);
  return nBinned;
}
  
unsigned int gpu_binTimeDistribution_and_countBinned(double* d_Itime, 
    double* d_time, const uint np,  const double tmin, const double dt,
    const uint nbins){
  GPUBASE* gpuBase = getGpuBase();
  unsigned int* d_tempu = gpuBase->d_blockTempu;

  gpuBinTimeDistribution transform(d_time, tmin, dt);
  unsigned int nTx = gpuBase->nReductionThreads; 
  if(nTx==256) {
    reduceAdd<256> reduction(d_tempu);
    gpuBinDriver(d_Itime, nbins, np, transform, reduction);
  }
  else if(nTx==512) {
    reduceAdd<512> reduction(d_tempu);
    gpuBinDriver(d_Itime, nbins, np, transform, reduction);
  }
  else if(nTx==1024) {
    reduceAdd<1024> reduction(d_tempu);
    gpuBinDriver(d_Itime, nbins, np, transform, reduction);
  }
  else if(nTx==2048) {
    reduceAdd<2048> reduction(d_tempu);
    gpuBinDriver(d_Itime, nbins, np, transform, reduction);
  }

  //setCacheConfig_binTime_reduceAdd();
  //gpuBinDriver(d_Itime, nbins, np, transform, reduction);

  unsigned int nBinned = 0;
  cudaMemcpy(&nBinned, d_tempu, sizeof(unsigned int), cudaMemcpyDeviceToHost);
  return nBinned;
}

void gpu_binTimeDistribution_and_findMax(double* d_Itime, double* d_time, 
    const uint np,  const double tmin, const double dt, const uint nbins, 
    double* Imax){
  GPUBASE* gpuBase = getGpuBase();
  unsigned int* d_tempu = gpuBase->d_blockTempu;

  gpuBinTimeDistribution transform(d_time, tmin, dt);
  unsigned int nTx = gpuBase->nReductionThreads; 
  if(nTx==256) {
    reduceMax<256> reduction(d_tempu);
    gpuBinDriver(d_Itime, nbins, np, transform, reduction);
  }
  else if(nTx==512) {
    reduceMax<512> reduction(d_tempu);
    gpuBinDriver(d_Itime, nbins, np, transform, reduction);
  }
  else if(nTx==1024) {
    reduceMax<1024> reduction(d_tempu);
    gpuBinDriver(d_Itime, nbins, np, transform, reduction);
  }
  else if(nTx==2048) {
    reduceMax<2048> reduction(d_tempu);
    gpuBinDriver(d_Itime, nbins, np, transform, reduction);
  }

  //setCacheConfig_binTime_reduceMax();
  //gpuBinDriver(d_Itime, nbins, np, transform, reduction);

  unsigned int uImax = 0;
  cudaMemcpy(&uImax, d_tempu, sizeof(unsigned int), cudaMemcpyDeviceToHost);

  *Imax = (double)uImax;
}


