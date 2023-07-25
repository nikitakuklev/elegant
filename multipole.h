#ifndef _MULTIPOLE_H
#define _MULTIPOLE_H

#ifdef __cplusplus
extern "C" {
#endif

#include "manual.h"


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
                                      long orderOverride, short *minOrder, short *maxOrder);

void randomizeErrorMultipoleFields(MULTIPOLE_DATA *randomMult);

#define ODD(j) ((j)%2)

#if TURBO_MODE >= 22
//static double expansion_coef_7[8] = {+1./(1*5040), +1/(1*720), -1/(2*120), -1./(6*24),
//                                     +1./(24*6), +1./(120*2), -1./(720*1), -1./(5040*1)};
static double expansion_coef_6[7] = {+1./(1*720), +1./(1*120), -1./(2*24), -1./(6*6),
                                     +1./(24*2), +1./(120*1), -1./(720*1)};
static double expansion_coef_5[6] = {+1./(1*120), +1./(1*24), -1./(2*6), -1./(6*2), +1./(24*1), +1./(120*1)};
static double expansion_coef_4[5] = {+1./(1*24), +1./(1*6), -1./(2*2), -1./(6*1), +1./(24*1)};
static double expansion_coef_3[4] = {+1./(1*6), +1./(1*2), -1./(2*1), -1./(6*1)};
static double expansion_coef_2[3] = {+1./(1*1), +1./(1*2), -1./(2*1)};
static double expansion_coef_1[2] = {+1./(1*1), +1./(1*1)};
static double expansion_coef_0[] = {+1./(1*1)};

static double *expansion_coefficients(const long n) {
  switch (n) {
    case 0:
      return expansion_coef_0;
    case 1:
      return expansion_coef_1;
    case 2:
      return expansion_coef_2;
    case 3:
      return expansion_coef_3;
    case 4:
      return expansion_coef_4;
    case 5:
      return expansion_coef_5;
    case 6:
      return expansion_coef_6;
//    case 7:
//      return expansion_coef_7;
    default:
      printf("Request for field expansion of order %ld seems too high???\n", n);
      return expansion_coefficients_dynamic(n);
  }
}

static double *expansions_coefficients_dynamic(const long n) {
  static double **expansion_coef = NULL;
  static long *orderDone = NULL;
  static long maxOrder = -1;
  long i;

  if (n <= maxOrder && orderDone[n])
    return (expansion_coef[n]);

  if (n > maxOrder) {
    expansion_coef = SDDS_Realloc(expansion_coef, sizeof(*expansion_coef) * (n + 1));
    orderDone = SDDS_Realloc(orderDone, sizeof(*orderDone) * (n + 1));
    for (i = maxOrder + 1; i <= n; i++)
      orderDone[i] = 0;
    maxOrder = n;
  }

  expansion_coef[n] = tmalloc(sizeof(**expansion_coef) * (n + 1));

  /* calculate expansion coefficients with signs for (x+iy)^n/n! */
  for (i = 0; i <= n; i++) {
    expansion_coef[n][i] = (ODD(i / 2) ? -1.0 : 1.0) / (dfactorial(i) * dfactorial(n - i));
  }
  orderDone[n] = 1;

  return (expansion_coef[n]);
}

static inline void fillPowerArray(const double x, double *xpow, const long order) {
  long i;
  xpow[0] = 1;
  for (i = 1; i <= order; i++) {
    xpow[i] = xpow[i - 1] * x;
  }
}

static inline void fillPowerArrayReverse(const double x, double *xpow, const long order) {
  long i;
  xpow[order] = 1;
  for (i = order-1; i >=0; i--) {
    xpow[i] = xpow[i + 1] * x;
  }
}

static void mkicks_fast(double *qx, double *qy,
                        double *delta_qx_return, double *delta_qy_return,
                        double *xpow, double *ypow,
                        const long order, const double KnL, const long skew) {
  double sum_Fx = 0, sum_Fy = 0;
  double *coef;
  int i;
  // Do we ever go above 8?
  // Assume that xpow is reversed such that highest order is at index 0
  printf("mkicks_fast order (%ld) skew (%ld)\n", order, skew);
  switch (order) {
    case 1:
      coef = expansion_coef_1;
      sum_Fy += coef[0] * xpow[0] * ypow[0];
      sum_Fx += coef[1] * xpow[1] * ypow[1];
    case 2:
      coef = expansion_coef_2;
      sum_Fy += coef[0] * xpow[0] * ypow[0];
      sum_Fx += coef[1] * xpow[1] * ypow[1];
      sum_Fy += coef[2] * xpow[2] * ypow[2];
    case 3:
      coef = expansion_coef_3;
      sum_Fy += coef[0] * xpow[0] * ypow[0];
      sum_Fx += coef[1] * xpow[1] * ypow[1];
      sum_Fy += coef[2] * xpow[2] * ypow[2];
      sum_Fx += coef[3] * xpow[3] * ypow[3];
    case 4:
      coef = expansion_coef_4;
      sum_Fy += coef[0] * xpow[0] * ypow[0];
      sum_Fx += coef[1] * xpow[1] * ypow[1];
      sum_Fy += coef[2] * xpow[2] * ypow[2];
      sum_Fx += coef[3] * xpow[3] * ypow[3];
      sum_Fy += coef[4] * xpow[4] * ypow[4];
    default:
      coef = expansion_coefficients(order);

      //odd
      for (i = 1; i <= order; i += 2)
        sum_Fx += coef[i] * xpow[i] * ypow[i];
      //even
      for (i = 0; i <= order; i += 2)
        sum_Fy += coef[i] * xpow[i] * ypow[i];
  }
  if (skew) {
    SWAP_DOUBLE(sum_Fx, sum_Fy);
    sum_Fx = -sum_Fx;
  }
  /* add the kicks */
  *qx -= KnL * sum_Fy;
  *qy += KnL * sum_Fx;
  if (delta_qx_return)
    *delta_qx_return -= KnL * sum_Fy;
  if (delta_qy_return)
    *delta_qy_return += KnL * sum_Fx;
}

#else
#if TURBO_MODE < 2
double *expansion_coefficients(const long n);
void fillPowerArray(double x, double *xpow, long order);
void apply_canonical_multipole_kicks(double *qx, double *qy,
                                     double *delta_qx_return,
                                     double *delta_qy_return,
                                     double *xpow, double *ypow,
                                     long order, double KnL, long skew);
#else
static double *expansion_coefficients(const long n) {
  static double **expansion_coef = NULL;
  static long *orderDone = NULL;
  static long maxOrder = -1;
  long i;

  if (n <= maxOrder && orderDone[n])
    return (expansion_coef[n]);

  if (n > maxOrder) {
    expansion_coef = SDDS_Realloc(expansion_coef, sizeof(*expansion_coef) * (n + 1));
    orderDone = SDDS_Realloc(orderDone, sizeof(*orderDone) * (n + 1));
    for (i = maxOrder + 1; i <= n; i++)
      orderDone[i] = 0;
    maxOrder = n;
  }

  expansion_coef[n] = tmalloc(sizeof(**expansion_coef) * (n + 1));

  /* calculate expansion coefficients with signs for (x+iy)^n/n! */
  for (i = 0; i <= n; i++) {
    expansion_coef[n][i] = (ODD(i / 2) ? -1.0 : 1.0) / (dfactorial(i) * dfactorial(n - i));
  }
  orderDone[n] = 1;

  return (expansion_coef[n]);
}

static void fillPowerArray(const double x, double *xpow, const long order) {
  long i;

  if (!xpow)
    bombElegant("Error: NULL pointer passed to fillPowerArray---Seek expert help!", NULL);

  xpow[0] = 1;
  for (i = 1; i <= order; i++) {
    xpow[i] = xpow[i - 1] * x;
  }
}

static void apply_canonical_multipole_kicks(double *qx, double *qy,
                                     double *delta_qx_return, double *delta_qy_return,
                                     double *xpow, double *ypow,
                                     const long order, const double KnL, const long skew) {
  long i;
  double sum_Fx, sum_Fy;
  double *coef;

  coef = expansion_coefficients(order);
  sum_Fx = sum_Fy = 0;

  //odd
  for (i = 1; i <= order; i += 2)
    sum_Fx += coef[i] * xpow[order - i] * ypow[i];
  //even
  for (i = 0; i <= order; i += 2)
    sum_Fy += coef[i] * xpow[order - i] * ypow[i];
  if (skew) {
    SWAP_DOUBLE(sum_Fx, sum_Fy);
    sum_Fx = -sum_Fx;
  }
  /* add the kicks */
  *qx -= KnL * sum_Fy;
  *qy += KnL * sum_Fx;
  if (delta_qx_return)
    *delta_qx_return -= KnL * sum_Fy;
  if (delta_qy_return)
    *delta_qy_return += KnL * sum_Fx;
}
#endif

#endif


#if TURBO_MODE >= 7
extern unsigned short expandHamiltonian;

static inline int convertSlopesToMomenta(double *qx, double *qy, double xp, double yp, double delta) {
  if (expandHamiltonian) {
    *qx = (1 + delta) * xp;
    *qy = (1 + delta) * yp;
  } else {
    double denom;
    denom = sqrt(1 + sqr(xp) + sqr(yp));
    *qx = (1 + delta) * xp / denom;
    *qy = (1 + delta) * yp / denom;
  }
  return 1;
}

static inline int convertMomentaToSlopes(double *xp, double *yp, double qx, double qy, double delta) {
  if (expandHamiltonian) {
    *xp = qx / (1 + delta);
    *yp = qy / (1 + delta);
  } else {
    double denom;
    if ((denom = sqr(1 + delta) - sqr(qx) - sqr(qy)) <= 0) {
      printWarningForTracking("Particle acquired undefined slopes when integrating through kick multipole.", NULL);
      return 0;
    }
    denom = sqrt(denom);
    *xp = qx / denom;
    *yp = qy / denom;
  }
  return 1;
}

#endif

#ifdef __cplusplus
}
#endif

#define EXSQRT(value, order) (order==0?sqrt(value):(1+0.5*((value)-1)))

#ifndef _GPU_MULTIPOLE_H_
extern unsigned long multipoleKicksDone;
#endif

#endif
