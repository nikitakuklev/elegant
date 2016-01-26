#include <gpu_kahan.h>
#include <gpu_kahan.hcu>

template<class FLOATTYPE, class TRANSFORM, unsigned int NTHREADS>
__global__ void sumReductionTrueKahanThreadfenceKernel(TRANSFORM transform,
    unsigned int n, FLOATTYPE* d_blockSum, FLOATTYPE* d_blockError, 
    unsigned int* retirementCount){

  // perform reduction in shared memory first
  __shared__ FLOATTYPE s_sum[NTHREADS];
  __shared__ FLOATTYPE s_error[NTHREADS];
  FLOATTYPE sum = 0;
  FLOATTYPE c = 0;
  
  // loop over lots of particles per thread to minimize tree reductions.
  for(unsigned int tid = threadIdx.x + blockIdx.x*blockDim.x; tid < n;
      tid += gridDim.x * blockDim.x){
    kahanSumDevice(sum, c, transform(tid));
  }
  FLOATTYPE error = -c;
  s_sum[threadIdx.x] = sum;
  s_error[threadIdx.x] = error;
  __syncthreads();

  // perform shared memory reduction
  reduceBlockKahanAdd<FLOATTYPE, NTHREADS>
      (s_sum, s_error, sum, error, threadIdx.x, blockDim.x);

  // write out per-block results
  if(threadIdx.x==0){
    d_blockSum[blockIdx.x] = sum = s_sum[0];
    d_blockError[blockIdx.x] = error = s_error[0];
  }

  // the last block reads the results of all previous blocks, performs a reduction,
  // and writes the result back to global code taken and modified from the
  // threadFenceReduction CUDA SDK example
  
  if(gridDim.x > 1){
    __shared__ volatile bool amLast;
    __threadfence();
  
    if(threadIdx.x==0){
      uint ticket = atomicInc(retirementCount, gridDim.x);
      amLast = (ticket == gridDim.x - 1);
    }
    __syncthreads();
  
    if( amLast ){
      FLOATTYPE sum = 0;
      FLOATTYPE error = 0;
      if(threadIdx.x < gridDim.x){
        kahanSumTwo(sum, d_blockSum[threadIdx.x], error, d_blockError[threadIdx.x]);
      }
  
      reduceBlockKahanAdd<FLOATTYPE, NTHREADS>
          (s_sum, s_error, sum, error, threadIdx.x, gridDim.x);
      
      if(threadIdx.x==0){
	d_blockSum[0] = sum = s_sum[0];
	d_blockError[0] = error = s_error[0];
      }  
      
    }//amLast
  }// gridDim.x > 1    
  
}

#include <gpu_base.h>

template<class TRANSFORM>
double gpuKahan(TRANSFORM transform, double* error, long n){
  GPUBASE* gpuBase = getGpuBase();
  double* d_tempf = gpuBase->d_blockTempf;
  unsigned int* d_retirementCount = gpuBase->d_retirementCount;  
  cudaMemset(d_retirementCount, 0, 128);
  
  unsigned int nTx = gpuBase->nReductionThreads;
  unsigned int nBx = (n + nTx -1 )/nTx;
  if(nBx > gpuBase->nReductionBlocks) nBx = gpuBase->nReductionBlocks;

  double* d_tempf_alpha = d_tempf;
  double* d_tempf_beta = d_tempf + gpuBase->nReductionBlocks;

  dim3 dimBlock(nTx,1,1);
  dim3 dimGrid(nBx,1,1);

  if (nTx==256)
    sumReductionTrueKahanThreadfenceKernel<double, TRANSFORM, 256>
      <<<dimGrid,dimBlock>>>(transform, n, d_tempf_alpha,
                             d_tempf_beta, d_retirementCount);
  else if (nTx==512)
    sumReductionTrueKahanThreadfenceKernel<double, TRANSFORM, 512>
      <<<dimGrid,dimBlock>>>(transform, n, d_tempf_alpha,
                             d_tempf_beta, d_retirementCount);
  else if (nTx==1024)
    sumReductionTrueKahanThreadfenceKernel<double, TRANSFORM, 1024>
      <<<dimGrid,dimBlock>>>(transform, n, d_tempf_alpha,
                             d_tempf_beta, d_retirementCount);
  else if (nTx==2048)
    sumReductionTrueKahanThreadfenceKernel<double, TRANSFORM, 2048>
      <<<dimGrid,dimBlock>>>(transform, n, d_tempf_alpha,
                             d_tempf_beta, d_retirementCount);

  double sum = 0;
  cudaMemcpy(&sum, d_tempf_alpha, sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(error, d_tempf_beta, sizeof(double), cudaMemcpyDeviceToHost);

  return sum;
}

template<class TRANSFORM>
void gpuKahanAsync(TRANSFORM transform, double* sum, double* error, long n){
  GPUBASE* gpuBase = getGpuBase();
  double *d_tempf=NULL, *pinsum=NULL, *pinerr=NULL;
  unsigned int* d_retirementCount=NULL;
  cudaStream_t* stream = (cudaStream_t*) getNextReductionStream(&d_tempf, 
                         &d_retirementCount, sum, &pinsum, error, &pinerr);
  
  unsigned int nTx = gpuBase->nReductionThreads;
  unsigned int nBx = (n + nTx -1 )/nTx;
  if(nBx > gpuBase->nReductionBlocks) nBx = gpuBase->nReductionBlocks;

  double* d_tempf_alpha = d_tempf;
  double* d_tempf_beta = d_tempf + gpuBase->nReductionBlocks;

  dim3 dimBlock(nTx,1,1);
  dim3 dimGrid(nBx,1,1);

  if (nTx==256)
    sumReductionTrueKahanThreadfenceKernel<double, TRANSFORM, 256>
      <<<dimGrid,dimBlock,0,*stream>>>(transform, n, d_tempf_alpha, 
      d_tempf_beta, d_retirementCount);
  else if (nTx==512)
    sumReductionTrueKahanThreadfenceKernel<double, TRANSFORM, 512>
      <<<dimGrid,dimBlock,0,*stream>>>(transform, n, d_tempf_alpha, 
      d_tempf_beta, d_retirementCount);
  else if (nTx==1024)
    sumReductionTrueKahanThreadfenceKernel<double, TRANSFORM, 1024>
      <<<dimGrid,dimBlock,0,*stream>>>(transform, n, d_tempf_alpha, 
      d_tempf_beta, d_retirementCount);
  else if (nTx==2048)
    sumReductionTrueKahanThreadfenceKernel<double, TRANSFORM, 2048>
      <<<dimGrid,dimBlock,0,*stream>>>(transform, n, d_tempf_alpha, 
      d_tempf_beta, d_retirementCount);

  cudaMemcpyAsync(pinsum, d_tempf_alpha, sizeof(double), 
                  cudaMemcpyDeviceToHost, *stream);
  cudaMemcpyAsync(pinerr, d_tempf_beta, sizeof(double),
                  cudaMemcpyDeviceToHost, *stream);
}

class basicKahanTransform{
public:
  double* d_a;
  basicKahanTransform(double* d_a):d_a(d_a){};

  __inline__ __device__ double operator()(const unsigned int& tid){
    return d_a[tid];
  }
};

double gpuKahan(double* d_a, double* error, long n){
  basicKahanTransform transform(d_a);
  return gpuKahan(transform, error, n);
}

void gpuKahanAsync(double* d_a, double* sum, double* error, long n){
  basicKahanTransform transform(d_a);
  gpuKahanAsync(transform, sum, error, n);
}

class P_average_transform{
public:
  double* d_coord5;
  double P_central;
  P_average_transform(double* d_coord5, double P_central):
    d_coord5(d_coord5),P_central(P_central){};
  
  __inline__ __device__ double operator()(const unsigned int& tid){
    double P_average = P_central*(1 + d_coord5[tid]);
    return P_average;
  }
};

double gpu_P_average_kahanSum(double* d_coord5, double* P_central,
                              double* error, long n){
  P_average_transform transform(d_coord5, *P_central);
  return gpuKahan(transform, error, n);
}
