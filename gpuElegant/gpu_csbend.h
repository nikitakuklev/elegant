#ifndef _GPU_CSBEND_H_
#define _GPU_CSBEND_H_

#include <csbend.h>

#ifdef __cplusplus
extern "C" {
#endif

long gpu_track_through_csbend(long n_part, CSBEND *csbend,
      double p_error, double Po, double **accepted, double z_start,
      double *sigmaDelta2, char *rootname, MAXAMP *maxamp, APERTURE_DATA *apFileData);

long gpu_track_through_csbendCSR(long n_part, CSRCSBEND *csbend,
       double p_error, double Po, double **accepted, double z_start,
       double z_end, CHARGE *charge, char *rootname, MAXAMP *maxamp, APERTURE_DATA *apFileData);

long gpu_track_through_driftCSR(long np, CSRDRIFT *csrDrift, double Po,
       double **accepted, double zStart, double revolutionLength,
       CHARGE *charge, char *rootname);

long gpu_track_through_driftCSR_Stupakov(long np, CSRDRIFT *csrDrift,
           double Po, double **accepted, double zStart, CHARGE *charge, char *rootname);

void gpu_exactDrift(long np, double length);

#ifdef __cplusplus
}
#endif

long gpu_binParticleCoordinate(double *hist, double *d_hist, long *maxBins, 
       double *lower, double *upper, double *binSize, long *bins, 
       double expansionFactor, long nParticles, long coordinateIndex);

#endif /* _GPU_CSBEND_H_ */
