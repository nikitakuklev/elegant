#ifndef _MULTIPOLE_H
#define _MULTIPOLE_H

#ifdef __cplusplus
extern "C" {
#endif

double *expansion_coefficients(long n);
void computeTotalErrorMultipoleFields(MULTIPOLE_DATA *totalMult,
                                      MULTIPOLE_DATA *systematicMult,
				      double systematicMultFactor,
                                      MULTIPOLE_DATA *edge1Mult,
                                      MULTIPOLE_DATA *edge2Mult,
                                      MULTIPOLE_DATA *randomMult,
				      double randomMultFactor,
                                      MULTIPOLE_DATA *steeringMult,
				      double steeringMultFactor,
                                      double KmL, long defaultOrder,
                                      long orderOverride, short minOrder, short maxOrder);

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
