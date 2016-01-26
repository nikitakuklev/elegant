#ifndef _GPU_FINAL_PROPS_H_
#define _GPU_FINAL_PROPS_H_

void gpu_rms_emittance(unsigned int i1, unsigned int i2, unsigned int n,
                       double* s11, double* s12, double* s22);

#if USE_MPI
void gpu_rms_emittance_p(long i1, long i2, long n,
                         double* s11, double* s12, double* s22);
#endif

double gpu_approximateBeamWidth(double fraction, double *d_hist, long nPart,
                                long iCoord);

#endif /* _GPU_FINAL_PROPS_H_ */
