#ifndef GPU_COMPUTE_CENTROIDS_H
#define GPU_COMPUTE_CENTROIDS_H

#include <gpu_track.h>
#include <gpu_base.h>

#ifdef __cplusplus /* if it's C++ use C linkage */
extern "C" {
#endif

void gpu_compute_centroids(double *centroid, long np);

void gpu_accumulate_beam_sums(BEAM_SUMS *sums, long n_part, double p_central, 
                              double mp_charge, long startPID, long endPID, 
                              unsigned long flags);

#ifdef __cplusplus
}
#endif 

#endif/* GPU_COMPUTE_CENTROIDS_H*/
