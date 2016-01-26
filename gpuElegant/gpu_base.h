#ifndef GPU_BASE_H
#define GPU_BASE_H
/**
 * Base singleton struct with pointers to data used by the GPU routines.
 *
 * GPU routines are included in the CPU routines of the same type 
 * with three basic statements:
 *
 *   #ifdef HAVE_GPU
 *     copyParticlesFromCpuToGpu(n_part);
 *     gpu_function_call(...);
 *     copyParticlesFromGpuToCpu(n_part);
 *     return;
 *   #endif
 *
 * However, the actual implementation is slightly more complicated as
 * GPU validation is also included. See gpu_verify.h.
 *
 * GPU singleton data is accessed through a statement such as
 *   struct GPUBASE* gpubase = getGpuBase();
 *   double** h_particles = gpubase->coord;
 *   double* d_particles = gpubase->d_particles;
 * etc.
 */

#include <gpu_memory_management.h>
#ifdef GPU_VERIFY
#include <gpu_verify.h>
#endif /* GPU_VERIFY */

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Struct to hold GPU stream data.
 */

/**
 * Struct to hold common data for GPU routines.
 */
struct GPUBASE{
  double* d_particles;             /* Device particles */
  double* d_accepted;              /* Accepted device particles */
  double* d_lostParticles;         /* Lost device particles */
  double* d_temp_particles;
  double* d_csbend;                /* csbend coeffecient storage */
  unsigned int* d_tempu_alpha;
  unsigned int* d_retirementCount; /* Used for reduction-type operations */
  double* d_blockTempf;            /* used for reduction-type operations */
  unsigned int* d_blockTempu;      /* used for reduction-type operations*/
  unsigned int nThreads;          /* number of cuda threads/block */
  unsigned int maxBlocks;         /* maximum number of cuda blocks */
  unsigned int nReductionThreads; /* number of cuda threads/block for reductions */
  unsigned int nReductionBlocks;  /* maximum number of cuda blocks for reductions */
  double** coord;                 /* Host particles */
  double** accepted;              /* Host accepted particles */
  double** lostParticles;         /* Host lost particles */
  double** doh_particles;         /* Device on host particles (for validation) */
  unsigned int gpu_array_pitch;   /* Device array nparticles + padding */
  unsigned int particlesOnGpu;    /* Particle location */
  unsigned int nOriginal;         /* original # of particles */
  unsigned int n_comp;            /* # of components per particle */
  unsigned int elementOnGpu;      /* boolean to run gpu routine */
  unsigned int copyForValidation; /* boolean to initialize cpu for validation */
  unsigned int nBins;             /* number of currently allocated bins */
  unsigned int gpuFtnLevel;       /* boolean to avoid recursive timings */
  long elemType;                  /* next element type */
  char* elemName;                 /* next element name */
  long isMaster;                  /* isMaster from elegant */
  double* full_d_particles;       /* pointer to full device particles if subset is selected */
  unsigned int particlesUnsorted; /* boolean for particles sort status */
};

/**
 * Returns a Singleton instance of gpubase. 
 */
struct GPUBASE* getGpuBase();

/**
 * Provides quick access to elementOnGpu within the gpubase singleton.
 */
unsigned int getElementOnGpu();

/**
 * Initialize the gpubase singleton struct
 * @param coord host particle pointer
 * @param nOriginal original number of particles
 * @param accepted accepted host pointer
 * @param lostParticles lostParticles host pointer
 * @param isMaster isMaster from elegant.c
 */
void gpuBaseInit(double** coord, unsigned int nOriginal, double** accepted, 
                 double** lostParticles, long isMaster);

/**
 * Deallocate the gpubase singleton struct
 */
void gpuBaseDealloc();

/**
 * Set GPU data for this element
 * @param veptr pointer to the eptr linked list of elements
 */
void setElementGpuData(void* veptr);

/**
 * Send particles to the CPU if on the GPU
 * @param routine rountine name for transfer
 * @return returns the host coord array
 */
double** forceParticlesToCpu(char* routine);

/**
 * Synchornize the lost particles array with the cpu.
 * @param nToTrack nparticles before the element
 * @param nLeft nparticles after the element
 */
void syncWithGpuLossBuffer(long nToTrack, long nLeft);

/**
 * Get the next reduction stream
 * @param d_tempf (output) global block reduction array for this stream
 * @param d_retirementCount (output) counter for this stream
 * @param val (input) pointer to host output location
 * @param pinval (output) output pointer to pinned val memory
 * @param val2 optional (input) pointer to second host output location
 * @param pinval2 optional (output) pointer to second host output location
 */
void* getNextReductionStream(double** d_tempf,
    unsigned int** d_retirementCount, double* value, double** pinval,
    double* value2, double** pinval2);

/**
 * Finish all reductions and write to host memory
 */
void finishReductionStreams();

/**
 * Generate a random normal distribution with a specified mean and sigma.
 * The seed is only used for the initial generation, after that
 * the curand generator is responsible for changing the seed.
 * @param d_ranarr pointer to the random number array stored in the 
 *                 global device memory
 * @param n_num number of random numbers to generate
 * @param mean normal distribution mean
 * @param sigma normal distribution standard deviation
 * @param seed random number generator seed
 * @note if n_num is odd, d_ranarr should be of size n_num+1
 */
void gpu_d_gauss_rn(double* d_ranarr, unsigned int n_num, double mean,
                       double sigma, double seed);

/**
 * Get the cuRAND random state array of size nOriginal.
 * This is useful if there is a need for on-demand thread random
 * number generation.
 * @param d_ranarr pointer to the random number array stored in the 
 *                 global device memory
 * @param n_num number of random numbers to generate
 * @param seed random number generator seed
 */
void* gpu_get_rand_state(double* d_ranarr, unsigned int n_num, double seed);

/**
 * Setup the curandState_t array for thread-based random number generation.
 * If a sufficiently sized state vector has already been created, this 
 * function does nothing. The seed is only used for the initial generation, 
 * after that the curand generator is responsible for changing the seed.
 * @param d_ranarr pointer to the random number array stored in the 
 *                 global device memory
 * @param n_num number of random numbers to generate
 * @param seed random number generator seed
 */
void gpu_init_rand_state(double* d_ranarr, unsigned int n_num,
                         double seed);

/**
 * Generate a random normal distribution within a limit.
 * The seed is only used for the initial generation, after that
 * the curand generator is responsible for changing the seed.
 * @param d_ranarr pointer to the random number array stored in the 
 *                 global device memory
 * @param n_num number of random numbers to generate
 * @param mean normal distribution mean
 * @param sigma normal distribution standard deviation
 * @param limit_in_sigmas cutoff in stdevs of the normal distribtion
 * @param seed random number generator seed
 * @note if n_num is odd, d_ranarr should be of size n_num+1
 */
void gpu_d_gauss_rn_lim(double* d_ranarr, unsigned int n_num, double mean,
                        double sigma, double limit_in_sigmas, double seed);

/**
 * Set the d_particles array pointer with an offset
 * @param offset offset index
 */
void selectParticleSubset(unsigned int offset);

/**
 * Reset the d_particles array pointer
 */
void resetToFullParticleArray();

/**
 * Sort particles by ID
 * @param np number of particles
 */
void sortByPID(unsigned int np);

#ifndef GPU_VERIFY
/**
 * Dummy function
 */
void startGpuTimer();
#endif

#ifdef __cplusplus
} /* extern C */
#endif /* __cplusplus */

#ifdef __CUDACC__
#include <string>
#if DEBUG || CUDA_DEBUG
#include <iostream>
#include <gpu_track.h>
#endif
/**
 * Check for CUDA errors
 */
inline void gpuErrorHandler(std::string diagnosis){
#if DEBUG || CUDA_DEBUG
  cudaError_t err = cudaDeviceSynchronize();
  if (err) {
    std::cout << "CUDA Error string at " << diagnosis << " is " 
              << cudaGetErrorString(err) << std::endl;
    bombElegant("CUDA error", NULL);
  }
#endif
}

/**
 * Check for CUDA errors
 */
inline void gpuErrorHandler(cudaError_t err, std::string diagnosis){
#if DEBUG || CUDA_DEBUG
  if (err) {
    std::cout << "CUDA Error string at " << diagnosis << " is " 
              << cudaGetErrorString(err) << std::endl;
    bombElegant("CUDA error", NULL);
  }
#endif
}
#endif /* __CUDACC__ */

#endif /* GPU_BASE_H */
