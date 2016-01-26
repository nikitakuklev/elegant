#ifdef GPU_VERIFY
#ifndef GPU_VERIFY_H
#define GPU_VERIFY_H
/**
 * Routines to verify the GPU implementations. When compiled
 * with -DENABLE_GPU_VERIFY these routine run the both the CPU and
 * GPU versions of each routine with the same initial state, and
 * then compare the results.
 *   CPU routines (e.g. track_particles), include relevant GPU
 * header, and, if HAVE_GPU is enabled and getElementOnGpu() 
 * returns true (the element is on the GPU), call the GPU function
 * then immediately return. However, if verify is true, 
 * additional validation steps and timing steps are taken.
 *   General structure for validation:
 * 1) The call to startGpuTimer() starts the GPU timer if validation
 *    is enabled.
 * 2) The GPU function is executed.
 * 3) After the GPU function call, is the validation code block.
 *    This code stops the GPU timer, starts the CPU timer, executes the 
 *    CPU version of the routine, stops the CPU timer, and compares the 
 *    resulting particles from both routines. A recursive call of the 
 *    routine is done to minimize the number of locations the cpu code 
 *    must be modified. getElementOnGpu returns false for this recursive
 *    call.
 *
 * For example, the following is the hook implementation of exactDrift
 * in csbend.c: 
 *#ifdef HAVE_GPU
 *  if(getElementOnGpu()){
 *    startGpuTimer();
 *    gpu_exactDrift(np, length);
 *#ifdef GPU_VERIFY     
 *    startCpuTimer();
 *    exactDrift(part, np, length);
 *    compareGpuCpu(np, "exactDrift");
 *#endif
 *    return;
 *  }
 *#endif
 */

#ifdef __cplusplus /* if it's C++ use C linkage */
extern "C" {
#endif

/**
 * Initialize the timer
 */
void initTimer();

/**
 * Deallocate the timer list
 */
void finalizeTimer();

/**
 * Start the GPU timer
 * @param synchronize if true the CPU and GPU use the same initial paricle state
 */
void startGpuTimer();

/**
 * Start the CPU timer (record the GPU time)
 */
void startCpuTimer();

/**
 * Stop the timer and record timings
 * @param name timed routine name
 */
void endTimer(char* name);

/**
 * Record timings
 * @param name timed routine name
 * @param gputime gputime of routine
 * @param cputime cputime of routine
 */
void recordTimings(char* name, float gputime, float cputime);

/**
 * Display the timings
 */
void displayTimings();

/**
 * Compare GPU and CPU results and aggregate timings
 * @param n_part number of particles
 * @param name compared routine name
 */
void compareGpuCpu(unsigned int n_part, char* name);

/**
 * Get the GPU beam sums array and sync (copy) the CPU beam sums array
 * @param cpu_beam_sums input beam sums array to copy 
 */
void* getGpuBeamSums(void* cpu_beam_sums);

/**
 * Copy and store CUDA computed reduction arrays.
 * @centroid input centroid array
 * @sigma input sigma array
 */
void copyReductionArrays(double* centroid, double* sigma);

/**
 * Compare passed CPU arrays with the stored GPU results.
 * @centroid CPU centroid array (can be NULL)
 * @sigma CPU sigma array (can be NULL)
 * @vsums CPU beam sums array (can be NULL)
 * @name name of the routine being verified
 */
void compareReductionArrays(double* centroid, double* sigma, void* vsums,
                            char* name);

/**
 * Compare two CSR_LAST_WAKE structs 
 * @vgpuCsrWake the gpu instance of csrWake
 * @vcpuCsrWake the cpu instance of csrWake
 */
void compareCSR_LAST_WAKE(void* vgpuCsrWake, void* vcpuCsrWake);

#ifdef __cplusplus
}
#endif

#endif /* GPU_VERIFY_H */
#endif /* GPU_VERIFY */
