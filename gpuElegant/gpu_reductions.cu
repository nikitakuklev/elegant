#include <typeinfo>

#include <gpu_reductions.h>
#include <gpu_reductions.hcu>
#include <gpu_base.h>

template<class T, class OPERATOR>
T gpuReduction(T* d_a, unsigned int n, OPERATOR op){
  GPUBASE* gpuBase = getGpuBase();

  unsigned int* d_retirementCount = gpuBase->d_retirementCount;  
  cudaMemset(d_retirementCount, 0, 128);

  unsigned int nTx = gpuBase->nReductionThreads;
  unsigned int nBx = (n + nTx -1 )/nTx;
  if(nBx > gpuBase->nReductionBlocks) nBx = gpuBase->nReductionBlocks;

  T* d_temp;
  if( typeid(T) == typeid(double)){
    d_temp = (T*)gpuBase->d_blockTempf;  
  }else if(typeid(T) == typeid(unsigned int)){      
    d_temp = (T*)gpuBase->d_blockTempu;
  }else{
    exit(1);
  }

  if (nTx==256)
    gpuReduceKernel<T, 256, OPERATOR><<<nBx,nTx>>>(d_a, d_temp, n,
      d_retirementCount, op );
  else if (nTx==512)
    gpuReduceKernel<T, 512, OPERATOR><<<nBx,nTx>>>(d_a, d_temp, n,
      d_retirementCount, op );
  else if (nTx==1024)
    gpuReduceKernel<T, 1024, OPERATOR><<<nBx,nTx>>>(d_a, d_temp, n,
      d_retirementCount, op );
  else if (nTx==2048)
    gpuReduceKernel<T, 2048, OPERATOR><<<nBx,nTx>>>(d_a, d_temp, n,
      d_retirementCount, op );

  T result = 0;
  cudaMemcpy(&result, d_temp, sizeof(T), cudaMemcpyDeviceToHost); 

  return result;
}

// instantiations
template double 
gpuReduction<double, Max<double> >(double* d_a, unsigned int a, Max<double> op);
template unsigned int 
gpuReduction<unsigned int, Max<unsigned int> >(unsigned int*, unsigned int, 
    Max<unsigned int> );

void gpuFindMinIndex(double* d_a, unsigned int n, double* minValue, 
       unsigned int* indexOfMinValue){
  GPUBASE* gpuBase = getGpuBase();
  double* d_tempf = gpuBase->d_blockTempf;
  unsigned int* d_tempu = gpuBase->d_blockTempu;
  unsigned int* d_retirementCount = gpuBase->d_retirementCount;  
  cudaMemset(d_retirementCount, 0, 128);
  
  unsigned int nTx = gpuBase->nReductionThreads;
  unsigned int nBx = (n + nTx -1 )/nTx;
  if(nBx > gpuBase->nReductionBlocks) nBx = gpuBase->nReductionBlocks;

  if (nTx==256)
    gpuComparatorReduceKernel<double, unsigned int, 
      LessThan<double, unsigned int>, 256 ><<<nBx, nTx>>>(d_a, 0, d_tempf,
      d_tempu, n, d_retirementCount, LessThan<double, unsigned int>() );
  else if (nTx==512)
    gpuComparatorReduceKernel<double, unsigned int, 
      LessThan<double, unsigned int>, 512 ><<<nBx, nTx>>>(d_a, 0, d_tempf,
      d_tempu, n, d_retirementCount, LessThan<double, unsigned int>() );
  else if (nTx==1024)
    gpuComparatorReduceKernel<double, unsigned int, 
      LessThan<double, unsigned int>, 1024 ><<<nBx, nTx>>>(d_a, 0, d_tempf,
      d_tempu, n, d_retirementCount, LessThan<double, unsigned int>() );
  else if (nTx==2048)
    gpuComparatorReduceKernel<double, unsigned int, 
      LessThan<double, unsigned int>, 2048 ><<<nBx, nTx>>>(d_a, 0, d_tempf,
      d_tempu, n, d_retirementCount, LessThan<double, unsigned int>() );

  cudaMemcpy(minValue, d_tempf, sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(indexOfMinValue, d_tempu, sizeof(unsigned int), cudaMemcpyDeviceToHost);
}

void gpuFindMaxIndex(double* d_a, unsigned int n, double* maxValue, 
    unsigned int* indexOfMaxValue){
  GPUBASE* gpuBase = getGpuBase();
  double* d_tempf = gpuBase->d_blockTempf;
  unsigned int* d_tempu = gpuBase->d_blockTempu;
  unsigned int* d_retirementCount = gpuBase->d_retirementCount;  
  cudaMemset(d_retirementCount, 0, 128);
  
  unsigned int nTx = gpuBase->nReductionThreads;
  unsigned int nBx = (n + nTx -1 )/nTx;
  if(nBx > gpuBase->nReductionBlocks) nBx = gpuBase->nReductionBlocks;

  if (nTx==256)
    gpuComparatorReduceKernel<double, unsigned int,
      GreaterThan<double, unsigned int>, 256 ><<<nBx, nTx>>>
      (d_a, 0, d_tempf, d_tempu, n, d_retirementCount,
      GreaterThan<double, unsigned int>() );
  else if (nTx==512)
    gpuComparatorReduceKernel<double, unsigned int,
      GreaterThan<double, unsigned int>, 512 ><<<nBx, nTx>>>
      (d_a, 0, d_tempf, d_tempu, n, d_retirementCount,
      GreaterThan<double, unsigned int>() );
  else if (nTx==1024)
    gpuComparatorReduceKernel<double, unsigned int,
      GreaterThan<double, unsigned int>, 1024 ><<<nBx, nTx>>>
      (d_a, 0, d_tempf, d_tempu, n, d_retirementCount,
      GreaterThan<double, unsigned int>() );
  else if (nTx==2048)
    gpuComparatorReduceKernel<double, unsigned int,
      GreaterThan<double, unsigned int>, 2048 ><<<nBx, nTx>>>
      (d_a, 0, d_tempf, d_tempu, n, d_retirementCount,
      GreaterThan<double, unsigned int>() );

  cudaMemcpy(maxValue, d_tempf, sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(indexOfMaxValue, d_tempu, sizeof(unsigned int), cudaMemcpyDeviceToHost);
}

#include <iostream>
double gpuReduceAdd(double* d_a, unsigned int n){
  
  GPUBASE* gpuBase = getGpuBase();
  double* d_tempf = gpuBase->d_blockTempf;

  unsigned int* d_retirementCount = gpuBase->d_retirementCount;  
  cudaMemset(d_retirementCount, 0, 128);

  unsigned int nTx = gpuBase->nReductionThreads;
  unsigned int nBx = (n + nTx -1 )/nTx;
  if(nBx > gpuBase->nReductionBlocks) nBx = gpuBase->nReductionBlocks;

  if (nTx==256)
    gpuReduceKernel<double, 256, Add<double> ><<<nBx,nTx>>>(d_a,
      d_tempf, n, d_retirementCount, Add<double>() );
  else if (nTx==512)
    gpuReduceKernel<double, 512, Add<double> ><<<nBx,nTx>>>(d_a,
      d_tempf, n, d_retirementCount, Add<double>() );
  else if (nTx==1024)
    gpuReduceKernel<double, 1024, Add<double> ><<<nBx,nTx>>>(d_a,
      d_tempf, n, d_retirementCount, Add<double>() );
  else if (nTx==2048)
    gpuReduceKernel<double, 2048, Add<double> ><<<nBx,nTx>>>(d_a,
      d_tempf, n, d_retirementCount, Add<double>() );

  double result = 0;
  cudaMemcpy(&result, d_tempf, sizeof(double), cudaMemcpyDeviceToHost); 
  return result;
}

void gpuReduceAddAsync(double* d_a, unsigned int n, double* sum){
  GPUBASE* gpuBase = getGpuBase();
  double* d_tempf=NULL, *pinsum=NULL;
  unsigned int* d_retirementCount=NULL;
  cudaStream_t* stream = (cudaStream_t*)
    getNextReductionStream(&d_tempf, &d_retirementCount, sum, &pinsum, 
      NULL, NULL);

  unsigned int nTx = gpuBase->nReductionThreads;
  unsigned int nBx = (n + nTx -1 )/nTx;
  if(nBx > gpuBase->nReductionBlocks) nBx = gpuBase->nReductionBlocks;

  if (nTx==256)
    gpuReduceKernel<double, 256, Add<double> ><<<nBx,nTx,0,*stream>>>
      (d_a, d_tempf, n, d_retirementCount, Add<double>() );
  else if (nTx==512)
    gpuReduceKernel<double, 512, Add<double> ><<<nBx,nTx,0,*stream>>>
      (d_a, d_tempf, n, d_retirementCount, Add<double>() );
  else if (nTx==1024)
    gpuReduceKernel<double, 1024, Add<double> ><<<nBx,nTx,0,*stream>>>
      (d_a, d_tempf, n, d_retirementCount, Add<double>() );
  else if (nTx==2048)
    gpuReduceKernel<double, 2048, Add<double> ><<<nBx,nTx,0,*stream>>>
      (d_a, d_tempf, n, d_retirementCount, Add<double>() );

  cudaMemcpyAsync(pinsum, d_tempf, sizeof(double), 
      cudaMemcpyDeviceToHost, *stream); 
}

long gpuReduceAddLong(long* d_a, unsigned int n, long* d_templ){
  GPUBASE* gpuBase = getGpuBase();

  unsigned int* d_retirementCount = gpuBase->d_retirementCount;  
  cudaMemset(d_retirementCount, 0, 128);
  
  unsigned int nTx = gpuBase->nReductionThreads;
  unsigned int nBx = (n + nTx -1 )/nTx;
  if(nBx > gpuBase->nReductionBlocks) nBx = gpuBase->nReductionBlocks;

  bool allocate_templ = d_templ ? false: true;
  if(allocate_templ) cudaMalloc( (void**)&d_templ, sizeof(long)*nBx);

  if (nTx==256)
    gpuReduceKernel<long, 256, Add<long> ><<<nBx,nTx,nTx*sizeof(long)>>>(d_a, 
      d_templ, n, d_retirementCount, Add<long>() );
  else if (nTx==512)
    gpuReduceKernel<long, 512, Add<long> ><<<nBx,nTx,nTx*sizeof(long)>>>(d_a, 
      d_templ, n, d_retirementCount, Add<long>() );
  else if (nTx==1024)
    gpuReduceKernel<long, 1024, Add<long> ><<<nBx,nTx,nTx*sizeof(long)>>>(d_a, 
      d_templ, n, d_retirementCount, Add<long>() );
  else if (nTx==2048)
    gpuReduceKernel<long, 2048, Add<long> ><<<nBx,nTx,nTx*sizeof(long)>>>(d_a, 
      d_templ, n, d_retirementCount, Add<long>() );

  long result = 0;
  cudaMemcpy(&result, d_templ, sizeof(long), cudaMemcpyDeviceToHost);

  if(allocate_templ) cudaFree(d_templ);
  return result;
}

double gpuReduceMin(double* d_a, unsigned int n){
  GPUBASE* gpuBase = getGpuBase();
  double* d_tempf = gpuBase->d_blockTempf;
  unsigned int* d_retirementCount = gpuBase->d_retirementCount;  
  cudaMemset(d_retirementCount, 0, 128);
  
  unsigned int nTx = gpuBase->nReductionThreads;
  unsigned int nBx = (n + nTx -1 )/nTx;
  if(nBx > gpuBase->nReductionBlocks) nBx = gpuBase->nReductionBlocks;

  if (nTx==256)
    gpuReduceKernel<double, 256, Min<double> ><<<nBx,nTx>>>(d_a, d_tempf, n,
      d_retirementCount, Min<double>() );
  else if (nTx==512)
    gpuReduceKernel<double, 512, Min<double> ><<<nBx,nTx>>>(d_a, d_tempf, n,
      d_retirementCount, Min<double>() );
  else if (nTx==1024)
    gpuReduceKernel<double, 1024, Min<double> ><<<nBx,nTx>>>(d_a, d_tempf, n,
      d_retirementCount, Min<double>() );
  else if (nTx==2048)
    gpuReduceKernel<double, 2048, Min<double> ><<<nBx,nTx>>>(d_a, d_tempf, n,
      d_retirementCount, Min<double>() );

  double result = 0;
  cudaMemcpy(&result, d_tempf, sizeof(double), cudaMemcpyDeviceToHost); 
  return result;
}

double gpuReduceMax(double* d_a,unsigned int n){
  GPUBASE* gpuBase = getGpuBase();
  double* d_tempf = gpuBase->d_blockTempf;

  unsigned int* d_retirementCount = gpuBase->d_retirementCount;  
  cudaMemset(d_retirementCount, 0, 128);
  
  unsigned int nTx = gpuBase->nReductionThreads;
  unsigned int nBx = (n + nTx -1 )/nTx;
  if(nBx > gpuBase->nReductionBlocks) nBx = gpuBase->nReductionBlocks;

  if (nTx==256)
    gpuReduceKernel<double, 256, Max<double> ><<<nBx,nTx>>>(d_a, d_tempf, n,
      d_retirementCount, Max<double>() );
  else if (nTx==512)
    gpuReduceKernel<double, 512, Max<double> ><<<nBx,nTx>>>(d_a, d_tempf, n,
      d_retirementCount, Max<double>() );
  else if (nTx==1024)
    gpuReduceKernel<double, 1024, Max<double> ><<<nBx,nTx>>>(d_a, d_tempf, n,
      d_retirementCount, Max<double>() );
  else if (nTx==2048)
    gpuReduceKernel<double, 2048, Max<double> ><<<nBx,nTx>>>(d_a, d_tempf, n,
      d_retirementCount, Max<double>() );

  double result = 0;
  cudaMemcpy(&result, d_tempf, sizeof(double), cudaMemcpyDeviceToHost); 
  return result;
}

#ifdef min
#undef min
#endif
#ifdef max
#undef max
#endif

void gpuReduceMinMax(double* d_a, unsigned int n, double* min, double* max){
  GPUBASE* gpuBase = getGpuBase();
  double* d_tempf = gpuBase->d_blockTempf;
  unsigned int* d_retirementCount = gpuBase->d_retirementCount;  
  cudaMemset(d_retirementCount, 0, 128);
  
  unsigned int nTx = gpuBase->nReductionThreads;
  unsigned int nBx = (n + nTx -1 )/nTx;
  if(nBx > gpuBase->nReductionBlocks) nBx = gpuBase->nReductionBlocks;

  double* d_min = d_tempf;
  double* d_max = d_tempf + gpuBase->nReductionBlocks;
 
  if (nTx==256)
    gpuDoubleReduceKernel<double, 256, nullTransform<double>, Min<double>, 
      Max<double> ><<<nBx,nTx>>>(d_min, d_max, n, d_retirementCount,
      nullTransform<double>(d_a), Min<double>(), Max<double>() );
  else if (nTx==512)
    gpuDoubleReduceKernel<double, 512, nullTransform<double>, Min<double>, 
      Max<double> ><<<nBx,nTx>>>(d_min, d_max, n, d_retirementCount,
      nullTransform<double>(d_a), Min<double>(), Max<double>() );
  else if (nTx==1024)
    gpuDoubleReduceKernel<double, 1024, nullTransform<double>, Min<double>, 
      Max<double> ><<<nBx,nTx>>>(d_min, d_max, n, d_retirementCount,
      nullTransform<double>(d_a), Min<double>(), Max<double>() );
  else if (nTx==2048)
    gpuDoubleReduceKernel<double, 2048, nullTransform<double>, Min<double>, 
      Max<double> ><<<nBx,nTx>>>(d_min, d_max, n, d_retirementCount,
      nullTransform<double>(d_a), Min<double>(), Max<double>() );

  cudaMemcpy(min,d_min,sizeof(double),cudaMemcpyDeviceToHost);
  cudaMemcpy(max,d_max,sizeof(double),cudaMemcpyDeviceToHost);  
}

void gpuReduceMinMaxAsync(double* d_a, unsigned int n, double* min, double* max){
  GPUBASE* gpuBase = getGpuBase();
  double* d_tempf=NULL, *pinmin=NULL, *pinmax=NULL;
  unsigned int* d_retirementCount=NULL;
  cudaStream_t* stream = (cudaStream_t*)
    getNextReductionStream(&d_tempf, &d_retirementCount, min, &pinmin, 
                           max, &pinmax);
  
  unsigned int nTx = gpuBase->nReductionThreads;
  unsigned int nBx = (n + nTx -1 )/nTx;
  if(nBx > gpuBase->nReductionBlocks) nBx = gpuBase->nReductionBlocks;

  double* d_min = d_tempf;
  double* d_max = d_tempf + gpuBase->nReductionBlocks;

  if (nTx==256)
    gpuDoubleReduceKernel<double, 256, nullTransform<double>, Min<double>, 
      Max<double> ><<<nBx,nTx,0,*stream>>>(d_min, d_max, n, d_retirementCount, 
      nullTransform<double>(d_a), Min<double>(), Max<double>() );
  else if (nTx==512)
    gpuDoubleReduceKernel<double, 512, nullTransform<double>, Min<double>, 
      Max<double> ><<<nBx,nTx,0,*stream>>>(d_min, d_max, n, d_retirementCount, 
      nullTransform<double>(d_a), Min<double>(), Max<double>() );
  else if (nTx==1024)
    gpuDoubleReduceKernel<double, 1024, nullTransform<double>, Min<double>, 
      Max<double> ><<<nBx,nTx,0,*stream>>>(d_min, d_max, n, d_retirementCount, 
      nullTransform<double>(d_a), Min<double>(), Max<double>() );
  else if (nTx==2048)
    gpuDoubleReduceKernel<double, 2048, nullTransform<double>, Min<double>, 
      Max<double> ><<<nBx,nTx,0,*stream>>>(d_min, d_max, n, d_retirementCount, 
      nullTransform<double>(d_a), Min<double>(), Max<double>() );

  cudaMemcpyAsync(pinmin,d_min,sizeof(double),cudaMemcpyDeviceToHost,*stream);
  cudaMemcpyAsync(pinmax,d_max,sizeof(double),cudaMemcpyDeviceToHost,*stream);  
}

void gpuReduceUMinMax(unsigned int* d_a, unsigned int n, unsigned int* min, 
                     unsigned int* max){
  GPUBASE* gpuBase = getGpuBase();
  unsigned int* d_tempu = gpuBase->d_blockTempu;
  unsigned int* d_retirementCount = gpuBase->d_retirementCount;  
  cudaMemset(d_retirementCount, 0, 128);
  
  unsigned int nTx = gpuBase->nReductionThreads;
  unsigned int nBx = (n + nTx -1 )/nTx;
  if(nBx > gpuBase->nReductionBlocks) nBx = gpuBase->nReductionBlocks;

  unsigned int* d_min = d_tempu;
  unsigned int* d_max = d_tempu + gpuBase->nReductionBlocks;
 
  if (nTx==256)
    gpuDoubleReduceKernel<unsigned int, 256, nullTransform<unsigned int>, Min<unsigned int>, 
      Max<unsigned int> ><<<nBx,nTx>>>(d_min, d_max, n, d_retirementCount,
      nullTransform<unsigned int>(d_a), Min<unsigned int>(), Max<unsigned int>() );
  else if (nTx==512)
    gpuDoubleReduceKernel<unsigned int, 512, nullTransform<unsigned int>, Min<unsigned int>, 
      Max<unsigned int> ><<<nBx,nTx>>>(d_min, d_max, n, d_retirementCount,
      nullTransform<unsigned int>(d_a), Min<unsigned int>(), Max<unsigned int>() );
  else if (nTx==1024)
    gpuDoubleReduceKernel<unsigned int, 1024, nullTransform<unsigned int>, Min<unsigned int>, 
      Max<unsigned int> ><<<nBx,nTx>>>(d_min, d_max, n, d_retirementCount,
      nullTransform<unsigned int>(d_a), Min<unsigned int>(), Max<unsigned int>() );
  else if (nTx==2048)
    gpuDoubleReduceKernel<unsigned int, 2048, nullTransform<unsigned int>, Min<unsigned int>, 
      Max<unsigned int> ><<<nBx,nTx>>>(d_min, d_max, n, d_retirementCount,
      nullTransform<unsigned int>(d_a), Min<unsigned int>(), Max<unsigned int>() );

  cudaMemcpy(min,d_min,sizeof(double),cudaMemcpyDeviceToHost);
  cudaMemcpy(max,d_max,sizeof(double),cudaMemcpyDeviceToHost);  
}

#include <thrust/device_ptr.h>
#include <thrust/reduce.h>
#include <thrust/functional.h>

double thrustReduceAdd(double* d_a, unsigned int n){
  thrust::device_ptr<double> tdv_a(d_a);

  return thrust::reduce(tdv_a, tdv_a + n);
}
