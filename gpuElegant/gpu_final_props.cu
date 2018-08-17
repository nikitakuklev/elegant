#include <gpu_track.h>

#include <gpu_base.h>
#include <gpu_reductions.hcu>
#include <gpu_csbend.h> // gpu_binParticleCoordinate

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

template <bool S11, bool S12, bool S22, unsigned int NTHREADS>
  __global__ void gpu_rms_kernel(double *d_particles, unsigned int particlePitch,
                                 unsigned int i1, unsigned int i2, double *d_result, unsigned int n)
{

  unsigned int tx = threadIdx.x;
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

  __shared__ volatile double s_v1[NTHREADS];
  __shared__ volatile double s_v2[NTHREADS];
  __shared__ volatile double s_s11[S11 ? NTHREADS : 1];
  __shared__ volatile double s_s12[S12 ? NTHREADS : 1];
  __shared__ volatile double s_s22[S22 ? NTHREADS : 1];

  s_v1[tx] = 0;
  s_v2[tx] = 0;
  if (S11)
    s_s11[tx] = 0;
  if (S12)
    s_s12[tx] = 0;
  if (S22)
    s_s22[tx] = 0;

  for (; tid < n; tid += gridDim.x * blockDim.x)
    {
      double v1 = 0.;
      double v2 = 0.;
      v1 = d_particles[tid + particlePitch * i1];
      s_v1[tx] += v1;
      v2 = d_particles[tid + particlePitch * i2];
      s_v2[tx] += v2;
      if (S11)
        s_s11[tx] += v1 * v1;
      if (S12)
        s_s12[tx] += v1 * v2;
      if (S22)
        s_s22[tx] += v2 * v2;
    }
  __syncthreads();

  double r_v1 = s_v1[tx];
  reduceBlock<double, NTHREADS, Add<double> >(s_v1, r_v1, tx, NTHREADS, Add<double>());
  double r_v2 = s_v2[tx];
  reduceBlock<double, NTHREADS, Add<double> >(s_v2, r_v2, tx, NTHREADS, Add<double>());
  if (S11)
    {
      double r_s11 = s_s11[tx];
      reduceBlock<double, NTHREADS, Add<double> >(s_s11, r_s11, tx, NTHREADS, Add<double>());
    }
  if (S12)
    {
      double r_s12 = s_s12[tx];
      reduceBlock<double, NTHREADS, Add<double> >(s_s12, r_s12, tx, NTHREADS, Add<double>());
    }
  if (S22)
    {
      double r_s22 = s_s22[tx];
      reduceBlock<double, NTHREADS, Add<double> >(s_s22, r_s22, tx, NTHREADS, Add<double>());
    }

  if (tx == 0)
    {
      atomicAddOAG((double *)(d_result), s_v1[0]);
      atomicAddOAG((double *)(d_result + 1), s_v2[0]);
      if (S11)
        atomicAddOAG((double *)(d_result + 2), s_s11[0]);
      if (S12)
        atomicAddOAG((double *)(d_result + 3), s_s12[0]);
      if (S22)
        atomicAddOAG((double *)(d_result + 4), s_s22[0]);
    }
}

#define NTHREADS 256

void gpu_rms_emittance(unsigned int i1, unsigned int i2, unsigned int n,
                       double *s11, double *s12, double *s22)
{
  GPUBASE *gpuBase = getGpuBase();
  double *d_particles = gpuBase->d_particles;
  unsigned int particlePitch = gpuBase->gpu_array_pitch;
  double *d_tempf = gpuBase->d_blockTempf;
  unsigned int *d_retirementCount = gpuBase->d_retirementCount;
  cudaMemset(d_retirementCount, 0, 128);

  const unsigned int nTx = NTHREADS;
  unsigned int nBx = (n + nTx - 1) / nTx;
  if (nBx > gpuBase->nReductionBlocks)
    nBx = gpuBase->nReductionBlocks;

  dim3 grid(nBx, 1, 1);
  dim3 block(nTx, 1, 1);

  cudaMemset(d_tempf, 0, sizeof(double) * 16);

  if (s11 && s12 && s22)
    gpu_rms_kernel<true, true, true, nTx><<<grid, block>>>(d_particles, particlePitch, i1, i2, d_tempf, n);
  else if (s11 && s12 && !s22)
    gpu_rms_kernel<true, true, false, nTx><<<grid, block>>>(d_particles, particlePitch, i1, i2, d_tempf, n);
  else if (s11 && !s12 && s22)
    gpu_rms_kernel<true, false, true, nTx><<<grid, block>>>(d_particles, particlePitch, i1, i2, d_tempf, n);
  else if (s11 && !s12 && !s22)
    gpu_rms_kernel<true, false, false, nTx><<<grid, block>>>(d_particles, particlePitch, i1, i2, d_tempf, n);
  else if (!s11 && s12 && s22)
    gpu_rms_kernel<false, true, true, nTx><<<grid, block>>>(d_particles, particlePitch, i1, i2, d_tempf, n);
  else if (!s11 && s12 && !s22)
    gpu_rms_kernel<false, true, false, nTx><<<grid, block>>>(d_particles, particlePitch, i1, i2, d_tempf, n);
  else if (!s11 && !s12 && s22)
    gpu_rms_kernel<false, false, true, nTx><<<grid, block>>>(d_particles, particlePitch, i1, i2, d_tempf, n);
  else
    gpu_rms_kernel<false, false, false, nTx><<<grid, block>>>(d_particles, particlePitch, i1, i2, d_tempf, n);

  double tempResult[5];
  cudaMemcpy(tempResult, d_tempf, sizeof(double) * 5, cudaMemcpyDeviceToHost);

  double sum_v1 = tempResult[0];
  double sum_v2 = tempResult[1];
  double sum_v1v1 = tempResult[2];
  double sum_v1v2 = tempResult[3];
  double sum_v2v2 = tempResult[4];

  if (s11)
    {
      *s11 = (sum_v1v1 - sum_v1 * sum_v1 / n) / n;
    }
  if (s12)
    {
      *s12 = (sum_v1v2 - sum_v1 * sum_v2 / n) / n;
    }
  if (s22)
    {
      *s22 = (sum_v2v2 - sum_v2 * sum_v2 / n) / n;
    }
}

#if USE_MPI
void gpu_rms_emittance_p(long i1, long i2, long n,
                         double *s11, double *s12, double *s22)
{
  long n_total = 0;

  if (notSinglePart)
    {
      if (isMaster)
        n = 0; /* The master will not contribute anything in this routine */
      MPI_Allreduce(&n, &n_total, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
    }
  else
    {
      if (!n)
        return;
    }

  if (!n_total)
    return;

  GPUBASE *gpuBase = getGpuBase();
  double *d_particles = gpuBase->d_particles;
  unsigned int particlePitch = gpuBase->gpu_array_pitch;
  double *d_tempf = gpuBase->d_blockTempf;
  unsigned int *d_retirementCount = gpuBase->d_retirementCount;
  cudaMemset(d_retirementCount, 0, 128);

  const unsigned int nTx = NTHREADS;
  unsigned int nBx = (n + nTx - 1) / nTx;
  if (nBx > gpuBase->nReductionBlocks)
    nBx = gpuBase->nReductionBlocks;

  dim3 grid(nBx, 1, 1);
  dim3 block(nTx, 1, 1);

  cudaMemset(d_tempf, 0, sizeof(double) * 16);

  if (s11 && s12 && s22)
    gpu_rms_kernel<true, true, true, nTx><<<grid, block>>>(d_particles, particlePitch, i1, i2, d_tempf, n);
  else if (s11 && s12 && !s22)
    gpu_rms_kernel<true, true, false, nTx><<<grid, block>>>(d_particles, particlePitch, i1, i2, d_tempf, n);
  else if (s11 && !s12 && s22)
    gpu_rms_kernel<true, false, true, nTx><<<grid, block>>>(d_particles, particlePitch, i1, i2, d_tempf, n);
  else if (s11 && !s12 && !s22)
    gpu_rms_kernel<true, false, false, nTx><<<grid, block>>>(d_particles, particlePitch, i1, i2, d_tempf, n);
  else if (!s11 && s12 && s22)
    gpu_rms_kernel<false, true, true, nTx><<<grid, block>>>(d_particles, particlePitch, i1, i2, d_tempf, n);
  else if (!s11 && s12 && !s22)
    gpu_rms_kernel<false, true, false, nTx><<<grid, block>>>(d_particles, particlePitch, i1, i2, d_tempf, n);
  else if (!s11 && !s12 && s22)
    gpu_rms_kernel<false, false, true, nTx><<<grid, block>>>(d_particles, particlePitch, i1, i2, d_tempf, n);
  else
    gpu_rms_kernel<false, false, false, nTx><<<grid, block>>>(d_particles, particlePitch, i1, i2, d_tempf, n);

  double tempResult[5];
  cudaMemcpy(tempResult, d_tempf, sizeof(double) * 5, cudaMemcpyDeviceToHost);

  double sum_v1 = tempResult[0];
  double sum_v2 = tempResult[1];
  double sum_v1v1 = tempResult[2];
  double sum_v1v2 = tempResult[3];
  double sum_v2v2 = tempResult[4];

  double total_sum_v1 = 0;
  double total_sum_v2 = 0;
  double total_sum_v1v1 = 0;
  double total_sum_v1v2 = 0;
  double total_sum_v2v2 = 0;

  if (s11 || s12)
    MPI_Allreduce(&sum_v1, &total_sum_v1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  if (s22 || s12)
    MPI_Allreduce(&sum_v2, &total_sum_v2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  if (s11)
    MPI_Allreduce(&sum_v1v1, &total_sum_v1v1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  if (s12)
    MPI_Allreduce(&sum_v1v2, &total_sum_v1v2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  if (s22)
    MPI_Allreduce(&sum_v2v2, &total_sum_v2v2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  if (s11)
    *s11 = (total_sum_v1v1 - total_sum_v1 * total_sum_v1 / n_total) / n_total;
  if (s22)
    *s22 = (total_sum_v2v2 - total_sum_v2 * total_sum_v2 / n_total) / n_total;
  if (s12)
    *s12 = (total_sum_v1v2 - total_sum_v1 * total_sum_v2 / n_total) / n_total;
}
#endif

#define ANALYSIS_BINS 10000

double gpu_approximateBeamWidth(double fraction, double *d_hist, long nPart,
                                long iCoord)
{
  double *hist, *cdf;
  long maxBins = ANALYSIS_BINS, bins = ANALYSIS_BINS, i50, iLo, iHi, i;
  double xMin, xMax, dx;
#if USE_MPI
  double *buffer;
#endif
  xMin = xMax = dx = 0;

  /* make histogram of the coordinate */
  hist = (double *)tmalloc(sizeof(*hist) * bins);
  gpu_binParticleCoordinate(hist, d_hist, &maxBins, &xMin, &xMax, &dx, &bins,
                            1.01, nPart, iCoord);

#if USE_MPI
  if (notSinglePart)
    { /* Master needs to know the information to write the result */
      buffer = (double *)malloc(sizeof(double) * bins);
      MPI_Allreduce(hist, buffer, bins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      memcpy(hist, buffer, sizeof(double) * bins);
      free(buffer);
    }
#endif
  /* sum histogram to get CDF */
  cdf = hist;
  for (i = 1; i < bins; i++)
    cdf[i] += cdf[i - 1];
  /* normalize CDF and find 50% point */
  i50 = bins / 2;
  for (i = 0; i < bins; i++)
    {
      cdf[i] /= cdf[bins - 1];
      if (cdf[i] < 0.50)
        i50 = i;
    }
  /* find locations containing half the indicated area around the 50% point */
  iLo = iHi = i50;
  for (i = i50; i < bins; i++)
    {
      if ((cdf[i] - 0.5) < fraction / 2)
        iHi = i;
      else
        break;
    }
  for (i = i50; i >= 0; i--)
    {
      if ((0.5 - cdf[i]) < fraction / 2)
        iLo = i;
      else
        break;
    }
  free(hist);
  return (iHi - iLo) * dx;
}
