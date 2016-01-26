#ifndef GPU_MEMORY_MANAGEMENT_H
#define GPU_MEMORY_MANAGEMENT_H

#ifdef __cplusplus /* if it's C++ use C linkage */
extern "C" {
#endif

/**
 * Allocate temporary arrays of unsigned ints for use in elegant
 * @param d_array double pointer to array.  
 * @param array_pitch pitch of pre-existing particle arrays.
 */
void allocateGpuUnsignedIntArray(unsigned int** ptr_to_d_tempu);

/**
 * Allocate GPU particles
 * @param d_array array to be allocated
 */
void allocateGpuParticles(double** ptr_to_d_array);

/**
 * Free GPU particles
 * @param d_array array to be free'd
 */
void freeGpuParticles(double** ptr_to_d_array);

/**
 * Perform GPU-based transpose and copy particles to GPU
 * @param n_part number of particles
 */
void copyParticlesFromCpuToGpu(unsigned int n_part);

/**
 * Perform GPU-based transpose and copy particles to CPU
 * @param n_part number of particles
 */
void copyParticlesFromGpuToCpu(unsigned int n_part);

/**
 * copy the component of a particle (particleIndex) from device to host
 * @param particleIndex the particle of interest
 * @param component the component to be copied.
 */
double getParticleData(unsigned int particleIndex, unsigned int component);

/**
 * fill an already allocated particle array on the host with a particle from the GPU
 * @param h_coord a seven-long array on the host.
 * @param d_particles pointer to particle array on device (SoA)
 * @param particlePitch pitch of SoA particle array
 * @param particleIndex the particle of interest
 */
void fillParticleData(double* h_coord, double* d_particles, unsigned int particlePitch, unsigned int particleIndex);

#ifdef __cplusplus
}
#endif
#endif

