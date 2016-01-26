#ifndef _GPU_MATTER_H_
#define _GPU_MATTER_H_

#ifdef __cplusplus
extern "C" {
#endif

/* setup */
void gpu_set_track_through_matter(long np, MATTER *matter, double Po,
       int impulseMode, int *multipleScattering, double *Nrad, double *L,
       double *theta_rms, int *sections0, double *prob, double *L1, double *K2,
       double *probBS, double *probER, double *dGammaFactor);

/* full function */
long gpu_track_through_matter(long np, MATTER *matter, double Po,
                              double **accepted, double z0);

#ifdef __cplusplus
}
#endif

#endif /* _GPU_MATTER_H_ */
