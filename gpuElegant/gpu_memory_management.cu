#include <stdio.h>
#include <iostream>

#include <gpu_track.h>

#include <gpu_base.h>

#define MINIMUM_ARRAY_SIZE 1024
/**
 * Allocate temporary arrays of unsigned ints for use in elegant
 * @param d_array double pointer to array.  
 */
void allocateGpuUnsignedIntArray(unsigned int **ptr_to_d_tempu)
{
  unsigned int *d_tempu;
  struct GPUBASE *gpubase = getGpuBase();
  unsigned int array_pitch = gpubase->gpu_array_pitch;
  size_t length = MINIMUM_ARRAY_SIZE;
  if (2 * array_pitch > length)
    length = 2 * array_pitch;
  cudaError_t err = cudaMalloc((void **)&d_tempu, sizeof(unsigned int) * length);
  gpuErrorHandler(err, "after allocating unsigned int arrays");

  *ptr_to_d_tempu = d_tempu;
}

/**
 * Allocate arrays for particles, in structure of arrays format
 * @param d_array pointer to an unallocated array
 */
void allocateGpuParticles(double **ptr_to_d_array)
{
  cudaError_t err;
  struct GPUBASE *gpubase = getGpuBase();
  unsigned int *array_pitch = &gpubase->gpu_array_pitch;
  unsigned int np = gpubase->nOriginal;
  unsigned int n_comp = gpubase->n_comp;
  if (np < MINIMUM_ARRAY_SIZE)
    np = MINIMUM_ARRAY_SIZE;
  size_t pitch = 0;
  double *d_array;
  err = cudaMallocPitch(&d_array, &pitch, np * sizeof(double), n_comp);
  *array_pitch = pitch / sizeof(double);
  gpuErrorHandler(err, "after cudaMallocPitch of SoA ");
  if (array_pitch == 0)
    {
      bombElegant("cudaMallocPitch for particles returned pitch = 0", NULL);
    }
  *ptr_to_d_array = d_array;
}

/**
 * Free GPU particles
 * @param d_array array to be free'd
 */
void freeGpuParticles(double **d_array)
{
  cudaError_t err;
  err = cudaFree(*d_array);
  gpuErrorHandler(err, "after freeing cuda array");
}

#ifdef NTRANSPOSE_PARTICLES
#  undef NTRANSPOSE_PARTICLES
#endif
#define NTRANSPOSE_PARTICLES 32

/**
 * use shared memory to both read and write in a coalesced way, from AoS to SoA
 * @param d_SoA array in structures-of-arrays format (pitched by n_part)
 * @param d_SoA_pitch pitch of d_SoA (pitched by n_part)
 * @param d_AoS array in array-of-structures format (pitched by n_comp)
 * @param d_AoS_pitch pitch of d_AoS (pitched by n_comp)
 * @param n_part number of particles
 * @param n_comp number of componenets per particles
 */
__global__ void sharedTransposeKernel(double *d_SoA, unsigned int d_SoA_pitch,
                                      double *d_AoS, unsigned int d_AoS_pitch, unsigned int n_part,
                                      unsigned int n_comp)
{
  // read in AoS order
  unsigned int rp = threadIdx.x / d_AoS_pitch;
  unsigned int rc = threadIdx.x - rp * d_AoS_pitch;
  // write in SoA order
  unsigned int wc = threadIdx.x / NTRANSPOSE_PARTICLES;
  unsigned int wp = threadIdx.x - wc * NTRANSPOSE_PARTICLES;
  // size: NTRANSPOSE_PARTICLES * d_AoS_pitch * sizeof(double)
  extern __shared__ double s_AoS[];

  for (unsigned int p_base = blockIdx.x * NTRANSPOSE_PARTICLES;
       p_base < n_part; p_base += gridDim.x * NTRANSPOSE_PARTICLES)
    {
      if (p_base + rp < n_part && rc < n_comp)
        {
          s_AoS[rp * d_AoS_pitch + rc] = d_AoS[(p_base + rp) * d_AoS_pitch + rc];
        }
      __syncthreads();
      if (p_base + wp < n_part && wc < n_comp)
        {
          d_SoA[p_base + wp + wc * d_SoA_pitch] = s_AoS[wp * d_AoS_pitch + wc];
        }
      __syncthreads();
    }
}

void sharedTranspose(double *d_SoA, unsigned int d_SoA_pitch, double *d_AoS,
                     unsigned int d_AoS_pitch, unsigned int n_part, unsigned int n_comp)
{
  struct GPUBASE *gpuBase = getGpuBase();
  unsigned int nTx = NTRANSPOSE_PARTICLES * d_AoS_pitch;
  unsigned int nBx = (n_part + NTRANSPOSE_PARTICLES - 1) / NTRANSPOSE_PARTICLES;
  if (nBx > gpuBase->maxBlocks)
    nBx = gpuBase->maxBlocks;
  dim3 dimBlock(nTx, 1, 1);
  dim3 dimGrid(nBx, 1, 1);
  size_t sharedMemorySize = nTx * sizeof(double);
  sharedTransposeKernel<<<dimGrid, dimBlock, sharedMemorySize>>>(d_SoA, d_SoA_pitch, d_AoS, d_AoS_pitch, n_part, n_comp);
}

/**
 * use shared memory to both read and write in a coalesced way, from SoA to AoS
 * @param d_AoS array in array-of-structures format (pitched by n_comp)
 * @param d_AoS_pitch pitch of d_AoS (pitched by n_comp)
 * @param d_SoA array in structures-of-arrays format (pitched by n_part)
 * @param d_SoA_pitch pitch of d_SoA (pitched by n_part)
 * @param n_part number of particles
 * @param n_comp number of componenets per particles
 */
__global__ void sharedInverseTransposeKernel(double *d_AoS,
                                             unsigned int d_AoS_pitch, double *d_SoA, unsigned int d_SoA_pitch,
                                             unsigned int n_part, unsigned int n_comp)
{
  // read in SoA format:  component c of particle p is indexed by d_SoA[p + c*SoA_pitch]
  unsigned int rc = threadIdx.x / NTRANSPOSE_PARTICLES;
  unsigned int rp = threadIdx.x - rc * NTRANSPOSE_PARTICLES;
  // write in AoS format:  component c of particle p is indexed by d_AoS[p*AoS_pitch + c]
  unsigned int wp = threadIdx.x / d_AoS_pitch;
  unsigned int wc = threadIdx.x - wp * d_AoS_pitch;
  // size: NTRANSPOSE_PARTICLES * d_AoS_pitch * sizeof(double)
  extern __shared__ double s_SoA[];

  for (unsigned int p_base = blockIdx.x * NTRANSPOSE_PARTICLES;
       p_base < n_part; p_base += gridDim.x * NTRANSPOSE_PARTICLES)
    {
      if (p_base + rp < n_part)
        {
          s_SoA[rp + rc * NTRANSPOSE_PARTICLES] = d_SoA[rp + p_base + rc * d_SoA_pitch];
        }
      __syncthreads();
      if (p_base + wp < n_part)
        {
          d_AoS[(p_base + wp) * d_AoS_pitch + wc] = s_SoA[wp + wc * NTRANSPOSE_PARTICLES];
        }
      __syncthreads();
    }
}

void sharedInverseTranspose(double *d_AoS, unsigned int d_AoS_pitch,
                            double *d_SoA, unsigned int d_SoA_pitch, unsigned int n_part,
                            unsigned int n_comp)
{
  struct GPUBASE *gpuBase = getGpuBase();
  unsigned int nTx = NTRANSPOSE_PARTICLES * d_AoS_pitch;
  unsigned int nBx = (n_part + NTRANSPOSE_PARTICLES - 1) / NTRANSPOSE_PARTICLES;
  if (nBx > gpuBase->maxBlocks)
    nBx = gpuBase->maxBlocks;
  dim3 dimBlock(nTx, 1, 1);
  dim3 dimGrid(nBx, 1, 1);
  size_t sharedMemorySize = NTRANSPOSE_PARTICLES * d_AoS_pitch * sizeof(double);
  //size_t sharedMemorySize = nTx * d_AoS_pitch * sizeof(double);
  sharedInverseTransposeKernel<<<dimGrid, dimBlock, sharedMemorySize>>>(d_AoS, d_AoS_pitch, d_SoA, d_SoA_pitch, n_part, n_comp);
}

#undef NTRANSPOSE_PARTICLES

/**
 * copy cpu particles to gpu and perform relevant GPU-based transpose
 * @param n_part number of particles
 */
void copyParticlesFromCpuToGpu(unsigned int n_part)
{
  struct GPUBASE *gpubase = getGpuBase();
  double *d_particles = gpubase->d_particles;
  double *d_temp_particles = gpubase->d_temp_particles;
  unsigned int array_pitch = gpubase->gpu_array_pitch;
  double **h_particles = gpubase->coord;
  unsigned int n_comp = gpubase->n_comp;
  static bool copyAccepted = false;
  cudaError_t err;

  if (gpubase->particlesOnGpu == 1)
    return;
  if (n_part <= 0)
    return;

  if (copyAccepted)
    {
      h_particles = gpubase->accepted;
      d_particles = gpubase->d_accepted;
    }

  err = cudaMemset(d_particles, 0, sizeof(double) * array_pitch * n_comp);
  gpuErrorHandler(err, "copyParticlesFromCpuToGpu::cudaMemset");
  err = cudaMemset(d_temp_particles, 0, sizeof(double) * array_pitch * n_comp);
  gpuErrorHandler(err, "copyParticlesFromCpuToGpu::cudaMemset (2)");

  err = cudaMemcpy(d_temp_particles, h_particles[0],
                   sizeof(double) * n_comp * n_part, cudaMemcpyHostToDevice);
  gpuErrorHandler(err, "copyParticlesFromCpuToGpu::cudaMemcpy");

  sharedTranspose(d_particles, array_pitch, d_temp_particles, n_comp, n_part, n_comp);
  gpuErrorHandler("copyParticlesFromCpuToGpu::sharedTranspose");

  if (gpubase->accepted && !copyAccepted)
    {
      copyAccepted = true;
      copyParticlesFromCpuToGpu(n_part);
    }

  if (copyAccepted)
    {
      copyAccepted = false;
    }
  else
    {
      gpubase->particlesOnGpu = 1;
      //std::cout << "Copied particles TO GPU." << std::endl;
    }
}

/**
 * Perform GPU-based transpose and copy particles to CPU
 * @param n_part number of particles
 */
void copyParticlesFromGpuToCpu(unsigned int n_part)
{
  struct GPUBASE *gpubase = getGpuBase();
  double *d_particles = gpubase->d_particles;
  double *d_temp_particles = gpubase->d_temp_particles;
  unsigned int array_pitch = gpubase->gpu_array_pitch;
  double **h_particles = gpubase->coord;
  unsigned int n_comp = gpubase->n_comp;
  static bool copyAccepted = false;
  cudaError_t err;

  if (n_part <= 0)
    return;

  /* Allow for device on host particles for comparison and accepted */
  if (gpubase->doh_particles != NULL)
    {
      h_particles = gpubase->doh_particles;
    }
  else if (copyAccepted)
    {
      h_particles = gpubase->accepted;
      d_particles = gpubase->d_accepted;
    }
  else if (!gpubase->particlesOnGpu && !gpubase->copyForValidation)
    return;

  err = cudaMemset(d_temp_particles, 0, sizeof(double) * array_pitch * n_comp);
  gpuErrorHandler(err, "after trying to memset d_temp_particles");

  sharedInverseTranspose(d_temp_particles, n_comp, d_particles,
                         array_pitch, n_part, n_comp);
  gpuErrorHandler("after shared inverse transpose from d_particles to d_temp_particles");

  err = cudaMemcpy(h_particles[0], d_temp_particles,
                   sizeof(double) * n_comp * n_part, cudaMemcpyDeviceToHost);
  gpuErrorHandler(err, "after linear copy from d_temp_particles to h_particles[0]");

  if (gpubase->accepted && !copyAccepted)
    {
      copyAccepted = true;
      copyParticlesFromGpuToCpu(n_part);
    }

  if (copyAccepted)
    {
      copyAccepted = false;
    }
  else if (gpubase->doh_particles != NULL)
    {
      //std::cout << "Copied particles FROM GPU (validation comparison)." << std::endl;
    }
  else if (gpubase->copyForValidation)
    {
      //std::cout << "Copied particles FROM GPU (validation initialize)." << std::endl;
      gpubase->copyForValidation = 0;
    }
  else if (gpubase->isMaster == -1)
    {
      //std::cout << "Copied particles FROM GPU for finalize." << std::endl;
      gpubase->particlesOnGpu = 0;
    }
  else if (!gpubase->elementOnGpu)
    {
      //std::cout << "Copied particles FROM GPU as element, "
      //          << gpubase->elemName << ", is of type "
      //          << gpubase->elemType << "." << std::endl;
      gpubase->particlesOnGpu = 0;
    }
}

/**
 * copy the component of a particle (particleIndex) from device to host
 * @param particleIndex the particle of interest
 * @param component the component to be copied.
 */
double getParticleData(unsigned int particleIndex, unsigned int component)
{
  struct GPUBASE *gpubase = getGpuBase();
  double *d_particles = gpubase->d_particles;
  unsigned int particlePitch = gpubase->gpu_array_pitch;
  double r;
  cudaMemcpy(&r, d_particles + particleIndex + particlePitch * component,
             sizeof(double), cudaMemcpyDeviceToHost);
  return r;
}

/**
 * fill an already allocated particle array on the host with a particle from the GPU
 * @param h_coord a seven-long array on the host.
 * @param d_particles pointer to particle array on device (SoA)
 * @param particlePitch pitch of SoA particle array
 * @param particleIndex the particle of interest
 */
void fillParticleData(double *h_coord, double *d_particles,
                      unsigned int particlePitch, unsigned int particleIndex)
{
  for (unsigned int component = 0; component < 7; component++)
    {
      cudaMemcpy(h_coord + component,
                 d_particles + particleIndex + particlePitch * component,
                 sizeof(double), cudaMemcpyDeviceToHost);
    }
}
