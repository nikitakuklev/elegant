#ifndef _MULTIPOLE_H
#define _MULTIPOLE_H

#ifdef __cplusplus
extern "C" {
#endif

double *expansion_coefficients(long n);
void applyRadialCanonicalMultipoleKicks(double *qx, double *qy, 
					double *sum_Fx_return, double *sum_Fy_return,
					double x, double y,
					long order, double KnL, long skew);

void apply_canonical_multipole_kicks(double *qx, double *qy,
                                   double *sum_Fx, double *sum_Fy,
                                   double x, double y,
                                   long order, double KnL, long skew);

void computeTotalErrorMultipoleFields(MULTIPOLE_DATA *totalMult,
                                      MULTIPOLE_DATA *systematicMult,
                                      MULTIPOLE_DATA *edgeMult,
                                      MULTIPOLE_DATA *randomMult,
                                      MULTIPOLE_DATA *steeringMult,
                                      double KmL, long rootOrder);

void randomizeErrorMultipoleFields(MULTIPOLE_DATA *randomMult);

#ifdef __cplusplus
}
#endif

#define EXSQRT(value, order) (order==0?sqrt(value):(1+0.5*((value)-1)))

#ifndef _GPU_MULTIPOLE_H_
extern unsigned long multipoleKicksDone;
#endif

#define ODD(j) ((j)%2)

#endif
