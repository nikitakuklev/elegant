#ifndef GPU_SIMPLE_RFCA_H
#define GPU_SIMPLE_RFCA_H

#include <gpu_base.h>
#include <gpu_track.h>

#ifdef __cplusplus /* if it's C++ use C linkage */
extern "C" {
#endif

long gpu_trackRfCavityWithWakes (long np, RFCA *rfca, double **accepted,
		      double *P_central, double zEnd, long iPass, RUN *run, 
		      CHARGE *charge, WAKE *wake,
		      TRWAKE *trwake, LSCKICK *LSCKick, long wakesAtEnd);

double gpu_findFiducialTime(long np, double s0, double sOffset,
                            double p0, unsigned long mode);

long gpu_track_through_rfcw(long np, RFCW *rfcw, double **accepted,
    double *P_central, double zEnd, RUN *run, long i_pass, CHARGE *charge);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#ifdef __CUDACC__
#ifdef sqr
#undef sqr
#define sqr(x) (x*x)
#endif
__inline__ __device__ void 
  gpu_add_to_particle_energy(gpuParticleAccessor& coord, double timeOfFlight, 
                             double Po, double dgamma)
{
  double gamma, gamma1, PRatio, P, P1, Pz1, Pz;

  P = Po*(1+coord[5]);                    /* old momentum */
  gamma1 = (gamma=sqrt(P*P+1)) + dgamma;  /* new gamma */
  if (gamma1<=1)
    gamma1 = 1+1e-7;
  P1 = sqrt(gamma1*gamma1-1);             /* new momentum */
  coord[5] = (P1-Po)/Po;                  

  /* adjust s for the new particle velocity */
  coord[4] = timeOfFlight*c_mks*P1/gamma1;

  /* adjust slopes so that Px and Py are conserved */
  Pz = P/sqrt(1+sqr(coord[1])+sqr(coord[3]));
  Pz1 = sqrt(Pz*Pz + gamma1*gamma1 - gamma*gamma);
  PRatio = Pz/Pz1;
  coord[1] *= PRatio;
  coord[3] *= PRatio;
}
#endif /* __CUDACC__ */

#endif /* GPU_SIMPLE_RFCA_H */
