#include <curand.h>
#include <curand_kernel.h>

#include <gpu_killParticles.hcu>
#include <gpu_base.h>
#include <gpu_particle_template_function.hcu>
#include <gpu_reductions.h>
#include <gpu_funcs.h> // gpu_offsetBeamCoordinates, gpu_rotateBeamCoordinates
#include <gpu_multipole.h>

#include <multipole.h>
#include <gpu_track.h>

#ifdef __cplusplus
extern "C" {
#endif
void gpu_exactDrift(long np, double length);
#ifdef __cplusplus
}
#endif

#ifdef sqr
#undef sqr
#endif
#define sqr(x) x*x
// Redefine these to avoid compiler errors (from gpu_track.h)
#undef SLOPE_LIMIT
#define SLOPE_LIMIT 1.0
#undef COORD_LIMIT
#define COORD_LIMIT 10.0

#define expansionOrderMax 10
// From multData and steeringMultData 
__device__ __constant__ int32_t multDataOrder[expansionOrderMax];
__device__ __constant__ double multDataKnL[expansionOrderMax];
__device__ __constant__ double multDataJnL[expansionOrderMax];
__device__ __constant__ int32_t steeringMultDataOrder[expansionOrderMax];
__device__ __constant__ double steeringMultDataKnL[expansionOrderMax];
__device__ __constant__ double steeringMultDataJnL[expansionOrderMax];
__device__ __constant__ int32_t edgeMultDataOrder[expansionOrderMax];
__device__ __constant__ double edgeMultDataKnL[expansionOrderMax];
__device__ __constant__ double edgeMultDataJnL[expansionOrderMax];
__device__ __constant__ double d_coef[expansionOrderMax*expansionOrderMax];

// Defined below
__device__ void gpu_apply_canonical_multipole_kicks(double *qx, double *qy,
    double *sum_Fx_return, double *sum_Fy_return, double x, double y,
    int order, double KnL, int skew);
__device__ void gpu_applyRadialCanonicalMultipoleKicks(double *qx, double *qy,
    double *sum_Fx_return, double *sum_Fy_return, double x, double y, 
    int order, double KnL, int skew);

class gpu_multipole_tracking2_kernel {
public:

  unsigned int *d_sortIndex;
  double *d_sigmaDelta2;
  curandState_t *state;
  double dx, dy, xkick, ykick, Po, rad_coef, isr_coef, KnL, drift, z_start;
  int order, sqrtOrder, n_parts, integ_order;
  // Associated MULTIPOLE_DATA is also in constant memory
  // From MULTIPOLE_DATA (*multData) 
  int multDataOrders; // set to -1 if no multData
  // From MULTIPOLE_DATA (*edgeMultData) 
  int edgeMultDataOrders; // set to -1 if no edgeMultData
  // From MULTIPOLE_DATA (*steeringMultData) 
  int steeringMultDataOrders; // set to -1 if no steeringMultData
  // From MULT_APERTURE_DATA (*apData)
  int present; 
  double xMax, xCen, yMax, yCen, srGaussianLimit;
  int radial;

  gpu_multipole_tracking2_kernel(unsigned int *d_sortIndex,
      double *d_sigmaDelta2, curandState_t *state, double dx, double dy,
      double xkick, double ykick, double Po, double rad_coef, double isr_coef,
      double KnL, double drift, double z_start, int order, int sqrtOrder,
      int n_parts, int integ_order, int multDataOrders, int edgeMultDataOrders,
      int steeringMultDataOrders, int present, double xMax, double xCen,
      double yMax, double yCen, int radial, double srGaussianLimit) :
    d_sortIndex(d_sortIndex), d_sigmaDelta2(d_sigmaDelta2),
    state(state), dx(dx), dy(dy), xkick(xkick), ykick(ykick), Po(Po),
    rad_coef(rad_coef), isr_coef(isr_coef), KnL(KnL), drift(drift),
    z_start(z_start), order(order), sqrtOrder(sqrtOrder), n_parts(n_parts), 
    integ_order(integ_order), multDataOrders(multDataOrders),
    edgeMultDataOrders(edgeMultDataOrders),
    steeringMultDataOrders(steeringMultDataOrders), present(present),
    xMax(xMax), xCen(xCen), yMax(yMax), yCen(yCen), radial(radial),
    srGaussianLimit(srGaussianLimit) {};

  __device__ unsigned int operator()(gpuParticleAccessor& coord) {
    unsigned int tid = coord.getParticleIndex();

    int particle_lost;
    double dzLoss;
    double *tSigmaDelta2 = NULL;
    if (d_sigmaDelta2) tSigmaDelta2 = &d_sigmaDelta2[tid];
    if (integ_order==4) {
      particle_lost = !gpu_integrate_kick_multipole_ord4(coord, KnL,
                         drift, &dzLoss, tSigmaDelta2, radial, srGaussianLimit);
    } else if (integ_order==2) {
      particle_lost = !gpu_integrate_kick_multipole_ord2(coord, KnL,
                         drift, &dzLoss, tSigmaDelta2, radial, srGaussianLimit);
    }
    if (particle_lost) {
      coord[4] = z_start+dzLoss;
      coord[5] = Po*(1+coord[5]);
      d_sortIndex[tid] = tid+coord.getParticlePitch();
      return 0;
    }
    d_sortIndex[tid] = tid;
    return 1;
  }

#define x  coord[0]    
#define xp coord[1]    
#define y  coord[2]    
#define yp coord[3]    
#define dp coord[5]    

  __device__ int
  gpu_integrate_kick_multipole_ord2(gpuParticleAccessor& coord, double KnL,
      double drift, double *dzLoss, double *sigmaDelta2, int radial,
      double srGaussianLimit) {
    double p, qx, qy, denom, beta0, beta1, s;
    double sum_Fx, sum_Fy;
    int i_kick, imult;
    
    drift = drift/n_parts/2.0;
    KnL = KnL/n_parts;
   
    // prepocessor defines 
    //x = coord[0];
    //xp = coord[1];
    //y = coord[2];
    //yp = coord[3];
    s  = 0;
    //dp = coord[5];
    p = Po*(1+dp);
    beta0 = p/sqrt(sqr(p)+1);
  
#if defined(IEEE_MATH)
    if (isnan(x) || isnan(xp) || isnan(y) || isnan(yp)) {
      return 0;
    }
#endif
    if (fabs(x)>COORD_LIMIT || fabs(y)>COORD_LIMIT ||
        fabs(xp)>SLOPE_LIMIT || fabs(yp)>SLOPE_LIMIT) {
      return 0;
    }
  
    /* apply steering corrector kick */
    xp += xkick/(1+dp)/2;
    yp += ykick/(1+dp)/2;
  
    /* calculate initial canonical momenta */
    denom = 1+sqr(xp)+sqr(yp);
    denom = EXSQRT(denom, sqrtOrder);
    qx = (1+dp)*xp/denom;
    qy = (1+dp)*yp/denom;
  
    /* apply steering corrector multipoles */
    if (steeringMultDataOrders>=0) {
      for (imult=0; imult<steeringMultDataOrders; imult++) {
        gpu_apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, x, y, 
                                        steeringMultDataOrder[imult], 
                                        steeringMultDataKnL[imult]*xkick/2, 0);
        gpu_apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, x, y, 
                                        steeringMultDataOrder[imult], 
                                        steeringMultDataJnL[imult]*ykick/2, 1);
      }
    }
    if (edgeMultDataOrders>=0) {
      for (imult=0; imult<edgeMultDataOrders; imult++) {
        gpu_apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, x, y, 
                                        edgeMultDataOrder[imult], 
                                        edgeMultDataKnL[imult], 0);
        gpu_apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, x, y, 
                                        edgeMultDataOrder[imult], 
                                        edgeMultDataJnL[imult], 1);
      }
    }
    // for numerical accuracy
    denom=(1+dp)*(1+dp)-(qx*qx+qy*qy);
    if (denom <= 0) return 0;
    //if ((denom=sqr(1+dp)-sqr(qx)-sqr(qy))<=0) {
    //  return 0;
    //}
    denom = EXSQRT(denom, sqrtOrder);
    xp = qx/denom;
    yp = qy/denom;
 
    *dzLoss = 0;
    for (i_kick=0; i_kick<n_parts; i_kick++) {
      if (drift) {
        x += xp*drift*(i_kick?2:1);
        y += yp*drift*(i_kick?2:1);
        s += drift*(i_kick?2:1)*sqrt(1 + sqr(xp) + sqr(yp));
        *dzLoss += drift*(i_kick?2:1);
      }
      if (present &&
          ((xMax && fabs(x + dx - xCen)>xMax) ||
           (yMax && fabs(y + dy - yCen)>yMax) )) {
        return 0;
      }
      
      if (!radial)
        gpu_apply_canonical_multipole_kicks(&qx, &qy, &sum_Fx, &sum_Fy,
            x, y, order, KnL, 0);
      else 
        gpu_applyRadialCanonicalMultipoleKicks(&qx, &qy, &sum_Fx, &sum_Fy,
            x, y, order, KnL, 0);
  
      /* do kicks for spurious multipoles */
      for (imult=0; imult<multDataOrders; imult++) {
        if (multDataKnL[imult]) 
          gpu_apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, x, y, 
                                          multDataOrder[imult], 
                                          multDataKnL[imult]/n_parts, 0);
        if (multDataJnL[imult]) 
          gpu_apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, x, y, 
                                          multDataOrder[imult], 
                                          multDataJnL[imult]/n_parts, 1);
      }
  
      // for numerical accuracy
      denom=(1+dp)*(1+dp)-(qx*qx+qy*qy);
      if (denom <= 0) return 0;
      //if ((denom=sqr(1+dp)-sqr(qx)-sqr(qy))<=0) {
      //  return 0;
      //}
      denom = EXSQRT(denom, sqrtOrder);
      xp = qx/denom;
      yp = qy/denom;
      if ((rad_coef || isr_coef) && drift) {
        double deltaFactor, F2, dsFactor;
        deltaFactor = sqr(1+dp);
        F2 = (sqr(sum_Fy)+sqr(sum_Fx))*sqr(KnL/(2*drift));
        dsFactor = xp*xp+yp*yp;
        dsFactor = EXSQRT(1+dsFactor, sqrtOrder)*2*drift;
        // comment out useless q->q transform
        //qx /= (1+dp);
        //qy /= (1+dp);
        if (rad_coef)
          dp -= rad_coef*deltaFactor*F2*dsFactor;
        if (isr_coef>0) {
          // Get a srGaussianLimit sigma (1) limited random number 
          double rand = curand_uniform_double(state);
          while(fabs(rand)>srGaussianLimit) 
            rand=curand_uniform_double(state);
          dp -= isr_coef*deltaFactor*pow(F2, 0.75)*sqrt(dsFactor)*rand;
        }
        if (sigmaDelta2)
          *sigmaDelta2 += sqr(isr_coef*deltaFactor)*pow(F2, 1.5)*dsFactor;
        //qx *= (1+dp);
        //qy *= (1+dp);
      }
    }
    if (drift) {
      /* go through final drift */
      x += xp*drift;
      y += yp*drift;
      s += drift*EXSQRT(1 + sqr(xp) + sqr(yp), sqrtOrder);
      *dzLoss += drift;
    }
    if (present &&
        ((xMax && fabs(x + dx - xCen)>xMax) ||
         (yMax && fabs(y + dy - yCen)>yMax) )) {
      return 0;
    }
  
    if (edgeMultDataOrders>=0) {
      for (imult=0; imult<edgeMultDataOrders; imult++) {
        gpu_apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, x, y, 
                                        edgeMultDataOrder[imult], 
                                        edgeMultDataKnL[imult], 0);
        gpu_apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, x, y, 
                                        edgeMultDataOrder[imult], 
                                        edgeMultDataJnL[imult], 1);
      }
    }

    /* apply steering corrector multipoles */
    if (steeringMultDataOrders>=0) {
      for (imult=0; imult<steeringMultDataOrders; imult++) {
        gpu_apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, x, y, 
                                        steeringMultDataOrder[imult], 
                                        steeringMultDataKnL[imult]*xkick/2, 0);
        gpu_apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, x, y, 
                                        steeringMultDataOrder[imult], 
                                        steeringMultDataJnL[imult]*ykick/2, 1);
      }
    }
    // for numerical accuracy
    denom=(1+dp)*(1+dp)-(qx*qx+qy*qy);
    if (denom <= 0) return 0;
    //if ((denom=sqr(1+dp)-sqr(qx)-sqr(qy))<=0) {
    //  return 0;
    //}
    denom = EXSQRT(denom, sqrtOrder);
    xp = qx/denom;
    yp = qy/denom;
  
    /* apply steering corrector kick */
    xp += xkick/(1+dp)/2;
    yp += ykick/(1+dp)/2;
  
    //coord[0] = x;
    //coord[1] = xp;
    //coord[2] = y;
    //coord[3] = yp;
    if (rad_coef || isr_coef) {
      p = Po*(1+dp);
      beta1 = p/sqrt(sqr(p)+1);
      coord[4] = beta1*(coord[4]/beta0 + 2*s/(beta0+beta1));
    }
    else 
      coord[4] += s;
    //coord[5] = dp;
  
#if defined(IEEE_MATH)
    if (isnan(x) || isnan(xp) || isnan(y) || isnan(yp)) {
      return 0;
    }
#endif
    if (fabs(x)>COORD_LIMIT || fabs(y)>COORD_LIMIT ||
        fabs(xp)>SLOPE_LIMIT || fabs(yp)>SLOPE_LIMIT) {
      return 0;
    }
    return 1;
  }

/* BETA is 2^(1/3) */
#define BETA 1.25992104989487316477

  __device__ int
  gpu_integrate_kick_multipole_ord4(gpuParticleAccessor& coord, double KnL,
      double drift, double *dzLoss, double *sigmaDelta2, int radial,
      double srGaussianLimit) {
    double p, qx, qy, denom, beta0, beta1, s;
    double sum_Fx, sum_Fy;
    int i_kick, step, imult;
    double dsh;
    const double driftFrac[4] = { 0.5/(2-BETA),  (1-BETA)/(2-BETA)/2,
                                 (1-BETA)/(2-BETA)/2,  0.5/(2-BETA)  } ;
    const double kickFrac[4] = { 1./(2-BETA),  -BETA/(2-BETA),
                                 1/(2-BETA),  0 } ;
    
    drift = drift/n_parts;
    KnL = KnL/n_parts;
  
    //x = coord[0];
    //xp = coord[1];
    //y = coord[2];
    //yp = coord[3];
    s  = 0;
    //dp = coord[5];
    p = Po*(1+dp);
    beta0 = p/sqrt(sqr(p)+1);
  
#if defined(IEEE_MATH)
    if (isnan(x) || isnan(xp) || isnan(y) || isnan(yp)) {
      return 0;
    }
#endif
    if (fabs(x)>COORD_LIMIT || fabs(y)>COORD_LIMIT ||
        fabs(xp)>SLOPE_LIMIT || fabs(yp)>SLOPE_LIMIT) {
      return 0;
    }
  
    /* apply steering corrector kick */
    xp += xkick/(1+dp)/2;
    yp += ykick/(1+dp)/2;
  
    /* calculate initial canonical momenta */
    qx = (1+dp)*xp/(denom=EXSQRT(1+sqr(xp)+sqr(yp), sqrtOrder));
    qy = (1+dp)*yp/denom;
  
    /* apply steering corrector multipoles */
    if (steeringMultDataOrders>=0) {
      for (imult=0; imult<steeringMultDataOrders; imult++) {
        if (steeringMultDataKnL[imult])
          gpu_apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, x, y, 
                                          steeringMultDataOrder[imult], 
                                          steeringMultDataKnL[imult]*xkick/2, 0);
        if (steeringMultDataJnL[imult]) 
          gpu_apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, x, y, 
                                          steeringMultDataOrder[imult], 
                                          steeringMultDataJnL[imult]*ykick/2, 1);
      }
    }

    if (edgeMultDataOrders>=0) {
      for (imult=0; imult<edgeMultDataOrders; imult++) {
        if (edgeMultDataKnL[imult])
          gpu_apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, x, y, 
                                          edgeMultDataOrder[imult], 
                                          edgeMultDataKnL[imult], 0);
        if (edgeMultDataJnL[imult]) 
          gpu_apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, x, y, 
                                          edgeMultDataOrder[imult], 
                                          edgeMultDataJnL[imult], 1);
      }
    }
    // for numerical accuracy
    denom=(1+dp)*(1+dp)-(qx*qx+qy*qy);
    if (denom <= 0) return 0;
    //if ((denom=sqr(1+dp)-sqr(qx)-sqr(qy))<=0) {
    //  return 0;
    //}
    xp = qx/(denom=EXSQRT(denom, sqrtOrder));
    yp = qy/denom;

    *dzLoss = 0;
    for (i_kick=0; i_kick<n_parts; i_kick++) {
      for (step=0; step<4; step++) {
        if (drift) {
          dsh = drift*driftFrac[step];
          x += xp*dsh;
          y += yp*dsh;
          s += dsh*EXSQRT(1 + sqr(xp) + sqr(yp), sqrtOrder);
          *dzLoss += dsh;
        }
        if (present &&
            ((xMax && fabs(x + dx - xCen)>xMax) ||
             (yMax && fabs(y + dy - yCen)>yMax) )) {
          return 0;
        }
  
        if (!kickFrac[step])
          break;

        if (!radial)
          gpu_apply_canonical_multipole_kicks(&qx, &qy, &sum_Fx, &sum_Fy, x, y, 
                                              order, KnL*kickFrac[step], 0);
        else 
          gpu_applyRadialCanonicalMultipoleKicks(&qx, &qy, &sum_Fx, &sum_Fy, x, y,
                                                 order, KnL*kickFrac[step], 0);

        /* do kicks for spurious multipoles */
        for (imult=0; imult<multDataOrders; imult++) {
          if (multDataKnL[imult]) {
            gpu_apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, x, y, 
                                            multDataOrder[imult], 
                                            multDataKnL[imult]*kickFrac[step]/n_parts,
                                            0);
          }
          if (multDataJnL[imult]) {
            gpu_apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, x, y, 
                                            multDataOrder[imult], 
                                            multDataJnL[imult]*kickFrac[step]/n_parts,
                                            1);
          }
        }
        // for numerical accuracy
        denom=(1+dp)*(1+dp)-(qx*qx+qy*qy);
        if (denom <= 0) return 0;
        //if ((denom=sqr(1+dp)-sqr(qx)-sqr(qy))<=0) {
        //  return 0;
        //}
        xp = qx/(denom=EXSQRT(denom, sqrtOrder));
        yp = qy/denom;
        if ((rad_coef || isr_coef) && drift) {
          double deltaFactor, F2, dsFactor, dsISRFactor;
          // comment out useless q->q transform
          //qx /= (1+dp);
          //qy /= (1+dp);
          deltaFactor = sqr(1+dp);
          F2 = (sqr(sum_Fy)+sqr(sum_Fx))*sqr(KnL/drift);
          dsFactor = xp*xp+yp*yp;
          dsFactor = EXSQRT(1+dsFactor, sqrtOrder);
          dsISRFactor = dsFactor*drift/3;   /* recall that kickFrac may be negative */
          dsFactor *= drift*kickFrac[step]; /* that's ok here, since we don't take sqrt */
          if (rad_coef)
            dp -= rad_coef*deltaFactor*F2*dsFactor;
          if (isr_coef>0) {
            // Get a srGaussianLimit sigma (sigma=1) limited random number 
            double rand = curand_uniform_double(state);
            while(fabs(rand)>srGaussianLimit) 
              rand=curand_uniform_double(state);
            dp -= isr_coef*deltaFactor*pow(F2, 0.75)*sqrt(dsISRFactor)*rand;
          }
          if (sigmaDelta2)
            *sigmaDelta2 += sqr(isr_coef*deltaFactor)*pow(F2, 1.5)*dsFactor;
          //qx *= (1+dp);
          //qy *= (1+dp);
        }
      }
    }
    
    if (present &&
        ((xMax && fabs(x + dx - xCen)>xMax) ||
         (yMax && fabs(y + dy - yCen)>yMax) )) {
      return 0;
    }
    
    if (edgeMultDataOrders>=0) {
      for (imult=0; imult<edgeMultDataOrders; imult++) {
        if (edgeMultDataKnL[imult]) 
          gpu_apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, x, y, 
                                          edgeMultDataOrder[imult], 
                                          edgeMultDataKnL[imult], 0);
        if (edgeMultDataJnL[imult]) 
          gpu_apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, x, y, 
                                          edgeMultDataOrder[imult], 
                                          edgeMultDataJnL[imult], 1);
      }
    }

    /* apply steering corrector multipoles */
    if (steeringMultDataOrders>=0) {
      for (imult=0; imult<steeringMultDataOrders; imult++) {
        if (steeringMultDataKnL[imult]) 
          gpu_apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, x, y, 
                                          steeringMultDataOrder[imult], 
                                          steeringMultDataKnL[imult]*xkick/2, 0);
        if (steeringMultDataJnL[imult]) 
          gpu_apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, x, y, 
                                          steeringMultDataOrder[imult], 
                                          steeringMultDataJnL[imult]*ykick/2, 1);
      }
    }
    // for numerical accuracy
    denom=(1+dp)*(1+dp)-(qx*qx+qy*qy);
    if (denom <= 0) return 0;
    //if ((denom=sqr(1+dp)-sqr(qx)-sqr(qy))<=0) {
    //  return 0;
    //}
    xp = qx/(denom=EXSQRT(denom, sqrtOrder));
    yp = qy/denom;
  
    /* apply steering corrector kick */
    xp += xkick/(1+dp)/2;
    yp += ykick/(1+dp)/2;
  
    //coord[0] = x;
    //coord[1] = xp;
    //coord[2] = y;
    //coord[3] = yp;
    if (rad_coef) {
      p = Po*(1+dp);
      beta1 = p/sqrt(sqr(p)+1);
      coord[4] = beta1*(coord[4]/beta0 + 2*s/(beta0+beta1));
    }
    else 
      coord[4] += s;
    //coord[5] = dp;
  
#if defined(IEEE_MATH)
    if (isnan(x) || isnan(xp) || isnan(y) || isnan(yp)) {
      return 0;
    }
#endif
    if (fabs(x)>COORD_LIMIT || fabs(y)>COORD_LIMIT ||
        fabs(xp)>SLOPE_LIMIT || fabs(yp)>SLOPE_LIMIT) {
      return 0;
    }
    return 1;
  }

#undef x
#undef xp
#undef y
#undef yp
#undef dp
};

extern unsigned long multipoleKicksDone;

extern "C" {

long gpu_multipole_tracking2(long n_part, ELEMENT_LIST *elem,
       double p_error, double Po, double **accepted, double z_start, 
       double x_max, double y_max, long elliptical, APERTURE_DATA *apFileData,
       double *sigmaDelta2) {
  double KnL(0);      /* integrated strength = L/(B.rho)*(Dx^n(By))_o for central momentum */
  double dx(0), dy(0), dz(0);  /* offsets of the multipole center */
  long order;         /* order (n) */
  long n_kicks(0), integ_order(0);
  long i_top, n_parts;
  double *coef;
  double drift(0);
  double tilt, rad_coef, isr_coef, xkick, ykick;
  KQUAD *kquad = NULL;
  KSEXT *ksext = NULL;
  KQUSE *kquse = NULL;
  KOCT *koct = NULL;
  static long sextWarning = 0, quadWarning = 0, octWarning = 0, quseWarning = 0;
  double lEffective = -1, lEnd = 0;
  short doEndDrift = 0;
  
  MULTIPOLE_DATA *multData = NULL, *steeringMultData = NULL, *edgeMultData = NULL;
  long sqrtOrder, freeMultData=0;
  MULT_APERTURE_DATA apertureData;
  double K2L;
  
  log_entry("multipole_tracking2");

  if (!elem)
    bombTracking("null element pointer (multipole_tracking2)");
  if (!elem->p_elem)
    bombTracking("null p_elem pointer (multipole_tracking2)");

  rad_coef = xkick = ykick = isr_coef = 0;
  sqrtOrder = 0;

  switch (elem->type) {
  case T_KQUAD:
    kquad = ((KQUAD*)elem->p_elem);
    n_kicks = kquad->n_kicks;
    order = 1;
    if ((lEffective = kquad->lEffective)<=0)
      lEffective = kquad->length;
    else {
      lEnd = (kquad->length-lEffective)/2;
      doEndDrift = 1;
    }
    if (kquad->bore)
      /* KnL = d^nB/dx^n * L/(B.rho) = n! B(a)/a^n * L/(B.rho) * (1+FSE) */
      KnL = kquad->B/kquad->bore*(particleCharge/(particleMass*c_mks*Po))*lEffective*(1+kquad->fse);
    else
      KnL = kquad->k1*lEffective*(1+kquad->fse);
    drift = lEffective;
    tilt = kquad->tilt;
    dx = kquad->dx;
    dy = kquad->dy;
    dz = kquad->dz;
    xkick = kquad->xkick*kquad->xKickCalibration;
    ykick = kquad->ykick*kquad->yKickCalibration;
    integ_order = kquad->integration_order;
    sqrtOrder = kquad->sqrtOrder?1:0;
    if (kquad->synch_rad)
      rad_coef = sqr(particleCharge)*pow(Po,3)/(6*PI*epsilon_o*sqr(c_mks)*particleMass); 
    isr_coef = particleRadius*sqrt(55.0/(24*sqrt(3))*pow(Po,5)*137.0359895);
    if (!kquad->isr || (kquad->isr1Particle==0 && n_part==1))
      /* Minus sign indicates we accumulate into sigmaDelta^2 only, don't perturb particles */
      isr_coef *= -1;
    if (kquad->length<1e-6 && (kquad->isr || kquad->synch_rad)) {
      rad_coef = isr_coef = 0;  /* avoid unphysical results */
      if (!quadWarning) {
        printf("**** Warning: one or more quadrupoles with length < 1e-6 have had SYNCH_RAD=0 and ISR=0 forced to avoid unphysical results.\n");
	quadWarning = 1;
      }
    }
    if (!kquad->multipolesInitialized) {
      /* read the data files for the error multipoles */
      readErrorMultipoleData(&(kquad->systematicMultipoleData),
                             kquad->systematic_multipoles, 0);
      readErrorMultipoleData(&(kquad->edgeMultipoleData),
                             kquad->edge_multipoles, 0);
      readErrorMultipoleData(&(kquad->randomMultipoleData),
                             kquad->random_multipoles, 0);
      readErrorMultipoleData(&(kquad->steeringMultipoleData), 
                             kquad->steering_multipoles, 1);
      kquad->multipolesInitialized = 1;
    }
    if (!kquad->totalMultipolesComputed) {
      computeTotalErrorMultipoleFields(&(kquad->totalMultipoleData),
                                       &(kquad->systematicMultipoleData),
                                       &(kquad->edgeMultipoleData),
                                       &(kquad->randomMultipoleData),
                                       &(kquad->steeringMultipoleData),
                                       KnL, 1);
      kquad->totalMultipolesComputed = 1;
    }
    multData = &(kquad->totalMultipoleData);
    edgeMultData = &(kquad->edgeMultipoleData);
    steeringMultData = &(kquad->steeringMultipoleData);
    break;
  case T_KSEXT:
    ksext = ((KSEXT*)elem->p_elem);
    n_kicks = ksext->n_kicks;
    order = 2;
    if (ksext->bore)
      /* KnL = d^nB/dx^n * L/(B.rho) = n! B(a)/a^n * L/(B.rho) * (1+FSE) */
      KnL = 2*ksext->B/sqr(ksext->bore)*(particleCharge/(particleMass*c_mks*Po))*ksext->length*(1+ksext->fse);
    else
      KnL = ksext->k2*ksext->length*(1+ksext->fse);
    drift = ksext->length;
    tilt = ksext->tilt;
    dx = ksext->dx;
    dy = ksext->dy;
    dz = ksext->dz;
    xkick = ksext->xkick*ksext->xKickCalibration;
    ykick = ksext->ykick*ksext->yKickCalibration;
    integ_order = ksext->integration_order;
    sqrtOrder = ksext->sqrtOrder?1:0;
    if (ksext->synch_rad)
      rad_coef = sqr(particleCharge)*pow(Po,3)/(6*PI*epsilon_o*sqr(c_mks)*particleMass);
    isr_coef = particleRadius*sqrt(55.0/(24*sqrt(3))*pow(Po,5)*137.0359895);
    if (!ksext->isr || (ksext->isr1Particle==0 && n_part==1))
      /* Minus sign indicates we accumulate into sigmaDelta^2 only, don't perturb particles */
      isr_coef *= -1;
    if (ksext->length<1e-6 && (ksext->isr || ksext->synch_rad)) {
      rad_coef = isr_coef = 0;  /* avoid unphysical results */
      if (!sextWarning) {
        printf("**** Warning: one or more sextupoles with length < 1e-6 have had SYNCH_RAD=0 and ISR=0 forced to avoid unphysical results.\n");
	sextWarning = 1;
      }
    }
    if (!ksext->multipolesInitialized) {
      /* read the data files for the error multipoles */
      readErrorMultipoleData(&(ksext->systematicMultipoleData),
                             ksext->systematic_multipoles, 0);
      readErrorMultipoleData(&(ksext->edgeMultipoleData),
                             ksext->edge_multipoles, 0);
      readErrorMultipoleData(&(ksext->randomMultipoleData),
                             ksext->random_multipoles, 0);
      readErrorMultipoleData(&(ksext->steeringMultipoleData), 
                             ksext->steering_multipoles, 1);
      ksext->multipolesInitialized = 1;
    }
    if (!ksext->totalMultipolesComputed) {
      computeTotalErrorMultipoleFields(&(ksext->totalMultipoleData),
                                       &(ksext->systematicMultipoleData),
                                       &(ksext->edgeMultipoleData),
                                       &(ksext->randomMultipoleData),
                                       &(ksext->steeringMultipoleData),
                                       KnL, 2);
      ksext->totalMultipolesComputed = 1;
    }
    multData = &(ksext->totalMultipoleData);
    edgeMultData = &(ksext->edgeMultipoleData);
    steeringMultData = &(ksext->steeringMultipoleData);
    break;
  case T_KOCT:
    koct = ((KOCT*)elem->p_elem);
    n_kicks = koct->n_kicks;
    order = 3;
    if (koct->bore)
      /* KnL = d^nB/dx^n * L/(B.rho) = n! B(a)/a^n * L/(B.rho) * (1+FSE) */
      KnL = 6*koct->B/pow(koct->bore, 3)*(particleCharge/(particleMass*c_mks*Po))*koct->length*(1+koct->fse);
    else
      KnL = koct->k3*koct->length*(1+koct->fse);
    drift = koct->length;
    tilt = koct->tilt;
    dx = koct->dx;
    dy = koct->dy;
    dz = koct->dz;
    integ_order = koct->integration_order;
    sqrtOrder = koct->sqrtOrder?1:0;
    if (koct->synch_rad)
      rad_coef = sqr(particleCharge)*pow(Po,3)/(6*PI*epsilon_o*sqr(c_mks)*particleMass);
    isr_coef = particleRadius*sqrt(55.0/(24*sqrt(3))*pow(Po,5)*137.0359895);
    if (!koct->isr || (koct->isr1Particle==0 && n_part==1))
      /* Minus sign indicates we accumulate into sigmaDelta^2 only, don't perturb particles */
      isr_coef *= -1;
    if (koct->length<1e-6 && (koct->isr || koct->synch_rad)) {
      rad_coef = isr_coef = 0;  /* avoid unphysical results */
      if (!octWarning) {
        printf("**** Warning: one or more octupoles with length < 1e-6 have had SYNCH_RAD=0 and ISR=0 forced to avoid unphysical results.\n");
	octWarning = 1;
      }
    }
    if (!koct->multipolesInitialized) {
      /* read the data files for the error multipoles */
      readErrorMultipoleData(&(koct->systematicMultipoleData),
                             koct->systematic_multipoles, 0);
      readErrorMultipoleData(&(koct->randomMultipoleData),
                             koct->random_multipoles, 0);
      koct->multipolesInitialized = 1;
    }
    if (!koct->totalMultipolesComputed) {
      computeTotalErrorMultipoleFields(&(koct->totalMultipoleData),
                                       &(koct->systematicMultipoleData),
                                       NULL,
                                       &(koct->randomMultipoleData),
                                       NULL,
                                       KnL, 3);
      koct->totalMultipolesComputed = 1;
    }
    multData = &(koct->totalMultipoleData);
    break;
  case T_KQUSE:
    /* Implemented as a quadrupole with sextupole as a secondary multipole */
    kquse = ((KQUSE*)elem->p_elem);
    n_kicks = kquse->n_kicks;
    order = 1;
    KnL = kquse->k1*kquse->length*(1+kquse->fse1);
    drift = kquse->length;
    tilt = kquse->tilt;
    dx = kquse->dx;
    dy = kquse->dy;
    dz = kquse->dz;
    integ_order = kquse->integration_order;
    sqrtOrder = 0;
    if (kquse->synch_rad)
      rad_coef = sqr(particleCharge)*pow(Po,3)/(6*PI*epsilon_o*sqr(c_mks)*particleMass); 
    isr_coef = particleRadius*sqrt(55.0/(24*sqrt(3))*pow(Po,5)*137.0359895);
    if (!kquse->isr || (kquse->isr1Particle==0 && n_part==1))
      /* Minus sign indicates we accumulate into sigmaDelta^2 only, don't perturb particles */
      isr_coef *= -1;
    if (kquse->length<1e-6 && (kquse->isr || kquse->synch_rad)) {
      rad_coef = isr_coef = 0;  /* avoid unphysical results */
      if (!quseWarning) {
        printf("**** Warning: one or more KQUSE's with length < 1e-6 have had SYNCH_RAD=0 and ISR=0 forced to avoid unphysical results.\n");
	quseWarning = 1;
      }
    }
    K2L = kquse->k2*kquse->length*(1+kquse->fse2);
    if (K2L) {
      multData = (MULTIPOLE_DATA*)tmalloc(sizeof(*multData));
      multData->orders = multData->initialized = 1;
      multData->randomized = 0;
      multData->order = (int32_t*)tmalloc(sizeof(*(multData->order))*1);
      multData->order[0] = 2;
      multData->KnL = (double*)tmalloc(sizeof(*(multData->KnL))*1);
      multData->KnL[0] = K2L;
      multData->JnL = NULL;
      freeMultData = 1;
    }
    break;
  default:
    fprintf(stdout, "error: multipole_tracking2() called for element %s--not supported!\n", elem->name);
    fflush(stdout);
    KnL = dx = dy = dz = tilt = drift = 0;
    integ_order = order = n_kicks = 0;
    exitElegant(1);
    break;
  }
  if (multData && !multData->initialized)
    multData = NULL;
  
  if (n_kicks<=0)
    bombTracking("n_kicks<=0 in multipole()");
  if (order<=0)
    bombTracking("order <= 0 in multipole()");
  if (integ_order!=2 && integ_order!=4) 
    bombTracking("multipole integration_order must be 2 or 4");
  
  if (!(coef = expansion_coefficients(order)))
    bombTracking("expansion_coefficients() returned null pointer (multipole_tracking)");

  i_top = n_part-1;
  if (integ_order==4) {
    if ((n_parts = ceil(n_kicks/4.0))<1)
      n_parts = 1;
    n_kicks = n_parts*4;
  }
  else
    n_parts = n_kicks;
  
  multipoleKicksDone += (i_top+1)*n_kicks;
  if (multData)
    multipoleKicksDone += (i_top+1)*n_kicks*multData->orders;

  setupMultApertureData(&apertureData, x_max, y_max, elliptical, tilt, apFileData, z_start+drift/2);
  
  struct GPUBASE* gpuBase = getGpuBase();
  unsigned int particlePitch = gpuBase->gpu_array_pitch;
  unsigned int* d_sortIndex = gpuBase->d_tempu_alpha;

  if (dx || dy || dz)
    gpu_offsetBeamCoordinates(n_part, dx, dy, dz);
  if (tilt)
    gpu_rotateBeamCoordinates(n_part, tilt);

  if (doEndDrift) {
    gpu_exactDrift(n_part, lEnd);
  }

  /* Fringe treatment, if any */
  switch (elem->type) {
  case T_KQUAD:
    if (kquad->edge1_effects>0)
      //TODO quadFringe
      std::cout << "Implement gpu_quadFringe" << std::endl;
      //quadFringe(particle, n_part, kquad->k1, kquad->fringeIntM, kquad->fringeIntP, -1, kquad->edge1_effects-1, kquad->edge1Linear, kquad->edge1NonlinearFactor);
    break;
  default:
    break;
  }

  // Copy multData and steeringMultData to device
  if (multData) {
    cudaMemcpyToSymbol(multDataOrder, multData->order, 
        sizeof(multDataOrder)*multData->orders, 0, cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol(multDataKnL, multData->KnL, 
        sizeof(multDataKnL)*multData->orders, 0, cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol(multDataJnL, multData->JnL, 
        sizeof(multDataJnL)*multData->orders, 0, cudaMemcpyHostToDevice);
  } 
  if (steeringMultData) {
    cudaMemcpyToSymbol(steeringMultDataOrder, steeringMultData->order, 
        sizeof(steeringMultDataOrder)*steeringMultData->orders, 0,
        cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol(steeringMultDataKnL, steeringMultData->KnL, 
        sizeof(steeringMultDataKnL)*steeringMultData->orders, 0,
        cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol(steeringMultDataJnL, steeringMultData->JnL, 
        sizeof(steeringMultDataJnL)*steeringMultData->orders, 0,
        cudaMemcpyHostToDevice);
  }
  if (edgeMultData) {
    cudaMemcpyToSymbol(edgeMultDataOrder, edgeMultData->order, 
        sizeof(edgeMultDataOrder)*edgeMultData->orders, 0, cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol(edgeMultDataKnL, edgeMultData->KnL, 
        sizeof(edgeMultDataKnL)*edgeMultData->orders, 0, cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol(edgeMultDataJnL, edgeMultData->JnL, 
        sizeof(edgeMultDataJnL)*edgeMultData->orders, 0, cudaMemcpyHostToDevice);
  } 

  // Ensure expansionOrderMax is sufficiently large
  int maxorder = order;
  for (int iorder=0; iorder < expansionOrderMax-1; ++iorder) {
    if (multData && multData->orders > iorder &&  
        multData->order[iorder] > maxorder)
      maxorder = multData->order[iorder];
    if (steeringMultData && steeringMultData->orders > iorder &&
        steeringMultData->order[iorder] > maxorder) 
      maxorder = steeringMultData->order[iorder];
    if (edgeMultData && edgeMultData->orders > iorder &&  
        edgeMultData->order[iorder] > maxorder)
      maxorder = edgeMultData->order[iorder];
  }
  if (maxorder > expansionOrderMax-1) 
    bombTracking("gpu_multipole_tracking2: Error: increase expansionOrderMax");

  // Copy expansion coefficients to constant memory
  for (int iorder=0; iorder <= maxorder; ++iorder) {
    coef = expansion_coefficients(iorder);
    cudaMemcpyToSymbol(d_coef, coef, sizeof(double)*(iorder+1),
                       sizeof(double)*(iorder*expansionOrderMax),
                       cudaMemcpyHostToDevice);
    gpuErrorHandler("gpu_multipole_tracking2::cudaMemcpyToSymbol");
  }
  
  double* d_gauss_rn = gpuBase->d_temp_particles + particlePitch;
  curandState_t* state = NULL;
  if (isr_coef>0)
    state =
      (curandState_t*)gpu_get_rand_state(d_gauss_rn, n_part, random_2(0));

  double* d_sigmaDelta2 = NULL;
  if (sigmaDelta2) {
    d_sigmaDelta2 = gpuBase->d_temp_particles;
    cudaMemset(d_sigmaDelta2, 0, particlePitch);
  }

  /*
  printf(" dx=%f dy=%f xkick=%f ykick=%f\n Po=%f rad_coef=%f isr_coef=%f KnL=%f\n",
         dx, dy, xkick, ykick, Po, rad_coef, isr_coef, KnL);
  printf(" drift=%f z_start=%f order=%ld\n sqrtOrder=%ld n_parts=%ld integ_order=%ld\n",
         drift, z_start, order, sqrtOrder, n_parts, integ_order);
  printf("multData->orders=%ld steeringMultData->orders=%ld\n",
         multData?multData->orders:-1, steeringMultData?steeringMultData->orders:-1);
  printf("apertureData.present=%d apertureData.xMax=%f apertureData.xCen=%f\n",
         apertureData.present, apertureData.xMax, apertureData.xCen);
  printf("apertureData.yMax=%f apertureData.yCen=%f\n",
         apertureData.yMax, apertureData.yCen);
  */

  n_part = killParticles(n_part, d_sortIndex, accepted,
    gpu_multipole_tracking2_kernel(d_sortIndex, d_sigmaDelta2, state,
    dx, dy, xkick, ykick, Po, rad_coef, isr_coef, KnL, drift, z_start, order,
    sqrtOrder, n_parts, integ_order, multData?multData->orders:-1, 
    edgeMultData?edgeMultData->orders:-1,
    steeringMultData?steeringMultData->orders:-1, apertureData.present,
    apertureData.xMax, apertureData.xCen, apertureData.yMax,
    apertureData.yCen, elem->type==T_KQUAD?kquad->radial:0, srGaussianLimit));
  gpuErrorHandler("gpu_multipole_tracking2::gpu_multipole_tracking2_kernel");

  if (sigmaDelta2) {
    *sigmaDelta2 = gpuReduceAdd(d_sigmaDelta2, n_part);
    *sigmaDelta2 /= n_part;
  }

  /* Fringe treatment, if any */
  switch (elem->type) {
  case T_KQUAD:
    if (kquad->edge2_effects>0)
      bombTracking("gpu_multipole_tracking2: quadFringe not implemented");
      //TODO quadFringe
      //quadFringe(particle, n_part, kquad->k1, kquad->fringeIntM, kquad->fringeIntP, 1, kquad->edge2_effects-1, kquad->edge2Linear, kquad->edge2NonlinearFactor);
    break;
  default:
    break;
  }
  if (doEndDrift) {
    gpu_exactDrift(n_part, lEnd);
  }
  
  if (tilt)
    gpu_rotateBeamCoordinates(n_part, -tilt);
  if (dx || dy || dz)
    gpu_offsetBeamCoordinates(n_part, -dx, -dy, -dz);

  if (freeMultData && !multData->copy) {
    if (multData->order)
      free(multData->order);
    if (multData->KnL)
      free(multData->KnL);
    free(multData);
  }

  log_exit("multipole_tracking2");
  return n_part;
}

} // extern "C"


__device__ void gpu_apply_canonical_multipole_kicks(double *qx, double *qy,
    double *sum_Fx_return, double *sum_Fy_return, double x, double y,
    int order, double KnL, int skew) {
  int i;
  double sum_Fx, sum_Fy;
  if (sum_Fx_return)
    *sum_Fx_return = 0;
  if (sum_Fy_return)
    *sum_Fy_return = 0;
  /* sum up the terms for the multipole expansion */
  for (i=sum_Fx=sum_Fy=0; i<=order; i++) {
    if (ODD(i))
      sum_Fx += d_coef[expansionOrderMax*order+i]*pow(x, order-i)*pow(y, i);
    else
      sum_Fy += d_coef[expansionOrderMax*order+i]*pow(x, order-i)*pow(y, i);
  }
  if (skew) {
    SWAP_DOUBLE(sum_Fx, sum_Fy);
    sum_Fx = -sum_Fx;
  }
  /* add the kicks */
  *qx -= KnL*sum_Fy;
  *qy += KnL*sum_Fx;
  if (sum_Fx_return)
    *sum_Fx_return = sum_Fx;
  if (sum_Fy_return)
    *sum_Fy_return = sum_Fy;
}

__device__ void gpu_applyRadialCanonicalMultipoleKicks(double *qx, double *qy,
     double *sum_Fx_return, double *sum_Fy_return, double x, double y, 
     int order, double KnL, int skew)
{
  int i;
  double sum_Fx, sum_Fy, xypow, ratio;
  if (sum_Fx_return)
    *sum_Fx_return = 0;
  if (sum_Fy_return)
    *sum_Fy_return = 0;
  if (x==0) {
    if (y==0)
      return;
    xypow = pow(y, order);
    i = order;
    ratio = 0;
  }
  else {
    xypow = pow(x, order);
    ratio = y/x;
    i = 0;
  }
  /* now sum up the terms for the multipole expansion */
  for (sum_Fx=sum_Fy=0; i<=order; i++) {
    if (ODD(i))
      sum_Fx -= d_coef[expansionOrderMax*order+i]*xypow;
    else
      sum_Fy += d_coef[expansionOrderMax*order+i]*xypow;
    xypow *= ratio;
  }
  if (skew) {
    SWAP_DOUBLE(sum_Fx, sum_Fy);
    sum_Fx = -sum_Fx;
  }
  /* add the kicks */
  *qx -= KnL*sum_Fy;
  *qy += KnL*sum_Fx;
  if (sum_Fx_return)
    *sum_Fx_return = sum_Fx;
  if (sum_Fy_return)
    *sum_Fy_return = sum_Fy;
}
