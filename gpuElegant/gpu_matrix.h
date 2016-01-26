#ifndef GPU_MATRIX_H
#define GPU_MATRIX_H

#include <gpu_track.h>
#include <gpu_base.h>

#ifdef __cplusplus /* if it's C++ use C linkage */
extern "C" {
#endif

void gpu_track_particles(VMATRIX *M, unsigned int n_part);

#ifdef __cplusplus
}
#endif
#endif
