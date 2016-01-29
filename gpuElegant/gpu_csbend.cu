#include <gpu_track.h>

#include <curand.h>
#include <curand_kernel.h>

#include <gpu_killParticles.hcu> 
#include <gpu_base.h>
#include <gpu_particle_template_function.hcu>
#include <gpu_particle_reduction.hcu>
#include <gpu_reductions.h>
#include <gpu_bin_time_distribution.h>
#include <gpu_matrix.h> // gpu_track_particles
#include <gpu_final_props.h> // gpu_rms_emittance
#include <gpu_fringe.hcu> // gpu_dipoleFringe
#include <gpu_multipole.h> // gpu_multipole_tracking2
#include <gpu_lsc.h> // gpu_addLSCKick

#include <gpu_csbend.h>
#include "csbend.h"

#ifdef sqr
#undef sqr
#endif
#define sqr(x) pow(x, 2)

#define expansionOrderMax 10
#define xOrderMax2 expansionOrderMax*expansionOrderMax

__device__ __constant__ double c_Fx_xy[xOrderMax2];
__device__ __constant__ double c_Fy_xy[xOrderMax2];
__device__ __constant__ double c_ksiTable[200];
__device__ __constant__ double c_FTable[200];

// device kernel declarations (defined below)
__device__ void 
gpu_addRadiationKick(double *Qx, double *Qy, double *dPoP, 
    double *sigmaDelta2, int sqrtOrder, double x, double h0, double Fx, 
    double Fy, double ds, double radCoef, double dsISR, double isrCoef, 
    int distributionBased, int includeOpeningAngle, 
    double meanPhotonsPerMeter, double normalizedCriticalEnergy, double Po,
    double *d_gauss_rn, curandState_t *state, double srGaussianLimit);
__device__ double gpu_pickNormalizedPhotonEnergy(double RN);
__device__ void 
gpu_integrate_csbend_ord2(double *Qf, double *Qi, double *sigmaDelta2, 
    double s, int n, int sqrtOrder, double rho0, 
    double p0, int &particle_lost, double &s_lost,
    double d_rho_actual, double d_rad_coef, double d_isrConstant, 
    int d_distributionBasedRadiation, int d_includeOpeningAngle,
    double d_meanPhotonsPerMeter0, double d_normalizedCriticalEnergy0,
    int d_expansionOrder, double *d_gauss_rn, curandState_t *state,
    double srGaussianLimit, double *d_refTrajectoryData, double d_refTrajectoryMode);
__device__ void 
gpu_integrate_csbend_ord4(double *Qf, double *Qi, double *sigmaDelta2, 
    double s, int n, int sqrtOrder, double rho0, double p0,
    int &particle_lost, double &s_lost,
    double d_rho_actual, double d_rad_coef, double d_isrConstant, 
    int d_distributionBasedRadiation, int d_includeOpeningAngle,
    double d_meanPhotonsPerMeter0, double d_normalizedCriticalEnergy0,
    int d_expansionOrder, double *d_gauss_rn, curandState_t *state,
    double srGaussianLimit, double *d_refTrajectoryData, double d_refTrajectoryMode);
__device__ void 
gpu_apply_edge_effects(double *x, double *xp, double *y, double *yp, 
                       double rho, double n, double beta, double he,
                       double psi, int which_edge);
__device__ void 
gpu_dipoleFringeSym(double *x, double *xp, double *y, double *yp,
                    double *dp, double rho, double inFringe, long edgeOrder, 
                    double K1, double edge, double gap, double fint, double Rhe);



// global kernel declarations (defined below)
void gpu_convertFromCSBendCoords(long np, double rho0, double cos_ttilt,
                                 double sin_ttilt, long ctMode);
void gpu_convertToCSBendCoords(long np, double rho0, double cos_ttilt, 
                               double sin_ttilt, long ctMode);
void gpu_initializeInterpolationTables();

// cpu function declarations
extern "C" {

#define RECORD_TRAJECTORY 1
#define SUBTRACT_TRAJECTORY 2
void convolveArrays1(double *output, long n, double *a1, double *a2);
void computeCSBENDFieldCoefficients(double *b, double h, long nonlinear,
                                    long expansionOrder);
void readWakeFilterFile(long *values, double **freq, double **real,
                        double **imag,  char *freqName, char *realName,
                        char *imagName, char *filename);

long correctDistribution(double *array, long npoints, double desiredSum);

void applyFilterTable(double *function, long bins, double dt, long fValues,
                      double *fFreq, double *fReal, double *fImag);

}

__device__ void 
gpu_convertToDipoleCanonicalCoordinates(double *Qi, double rho, int sqrtOrder)
{
  double f;
  f = (1 + Qi[5])/EXSQRT(1 + sqr(Qi[1]) + sqr(Qi[3]), sqrtOrder);
  Qi[1] *= f;
  Qi[3] *= f;
}

__device__ void 
gpu_convertFromDipoleCanonicalCoordinates(double *Qi, double rho, int sqrtOrder)
{
  double f;
  f = 1/EXSQRT(sqr(1+Qi[5])-sqr(Qi[1])-sqr(Qi[3]), sqrtOrder);
  Qi[1] *= f;
  Qi[3] *= f;
}

__device__ int gpu_inversePoissonCDF(double mu, double C)
{
  double sum, expMinusMu, term;
  int r, rMax;

  r = 0;
  if ((rMax = 50*mu)<10)
    rMax = 10;
  expMinusMu = exp(-mu);
  term = sum = expMinusMu;
  while (r<=rMax && C>=sum) {
    term *= mu/(++r);
    sum += term;
  }
  return r;
}

__device__ void gpu_computeCSBENDFields(double *Fx, double *Fy, 
    double x, double y, double d_expansionOrder)
{
  double xp[11], yp[11];
  double sumFx=0, sumFy=0;
  int i, j;
  
  xp[0] = yp[0] = 1;
  for (i=1; i<d_expansionOrder; i++) {
    xp[i] = xp[i-1]*x;
    yp[i] = yp[i-1]*y;
  }
  
  for (i=0; i<d_expansionOrder; i++)
    for (j=1; j<d_expansionOrder-i; j+=2) {
      if (c_Fx_xy[i*expansionOrderMax+j]) {
        sumFx += c_Fx_xy[i*expansionOrderMax+j]*xp[i]*yp[j];
      }
    }
  *Fx = sumFx;

  for (i=0; i<d_expansionOrder; i++)
    for (j=0; j<d_expansionOrder-i; j+=2) {
      if (c_Fy_xy[i*expansionOrderMax+j]) {
        sumFy += c_Fy_xy[i*expansionOrderMax+j]*xp[i]*yp[j];
      }
    }
  *Fy = sumFy;
}

class gpu_track_through_csbend_kernel{
public:
  unsigned int *d_sortIndex;
  double *d_sigmaDelta2;
  double *d_gauss_rn;
  curandState_t *state;
  double Po, z_start, srGaussianLimit;
  double dxi, dyi, dzi, dxf, dyf, dzf;
  double cos_ttilt, sin_ttilt;
  double e1, e2, he1, he2, psi1, psi2, n;
  double dcet0, dcet1, dcet2, dcet3;
  double e1_kick_limit, e2_kick_limit;
  // From CSBEND struct
  double length, b0, hgap, fint, h1, h2;
  int edge1_effects, edge2_effects, edge_order, n_kicks;
  int integration_order, d_expansionOrder, sqrtOrder;
  unsigned long edgeFlags;
  // From csbend.h globals
  double d_rho0, d_rho_actual, d_rad_coef, d_isrConstant;
  double d_meanPhotonsPerMeter0, d_normalizedCriticalEnergy0;
  int d_distributionBasedRadiation, d_includeOpeningAngle;
  double *d_refTrajectoryData;
  double d_refTrajectoryMode;

  gpu_track_through_csbend_kernel(unsigned int *d_sortIndex,
      double *d_sigmaDelta2, double *d_gauss_rn, curandState_t *state,
      double Po, double z_start, double dxi,
      double dyi, double dzi, double dxf, double dyf, double dzf,
      double cos_ttilt, double sin_ttilt, double e1, double e2,
      double he1, double he2, double psi1, double psi2, double n,
      double dcet0, double dcet1, double dcet2, double dcet3,
      double e1_kick_limit, double e2_kick_limit, double sqrtOrder,
      double length, double b0, double hgap, double fint, double h1, 
      double h2, int edge1_effects, int edge2_effects, 
      int edge_order, int n_kicks, int integration_order, 
      unsigned long edgeFlags, double d_rho0, double d_rho_actual,
      double d_rad_coef, double d_isrConstant, int d_expansionOrder,
      double d_meanPhotonsPerMeter0, double d_normalizedCriticalEnergy0,
      int d_distributionBasedRadiation, int d_includeOpeningAngle,
      double srGaussianLimit, double *d_refTrajectoryData, double d_refTrajectoryMode) :
    d_sortIndex(d_sortIndex), d_sigmaDelta2(d_sigmaDelta2),
    d_gauss_rn(d_gauss_rn), state(state), Po(Po), z_start(z_start), 
    dxi(dxi), dyi(dyi), dzi(dzi), dxf(dxf), dyf(dyf),
    dzf(dzf), cos_ttilt(cos_ttilt), sin_ttilt(sin_ttilt), e1(e1), e2(e2),
    he1(he1), he2(he2), psi1(psi1), psi2(psi2), n(n), dcet0(dcet0), 
    dcet1(dcet1), dcet2(dcet2), dcet3(dcet3), e1_kick_limit(e1_kick_limit),
    e2_kick_limit(e2_kick_limit), sqrtOrder(sqrtOrder), length(length),
    b0(b0), hgap(hgap), fint(fint), h1(h1), h2(h2), 
    edge1_effects(edge1_effects), edge2_effects(edge2_effects),
    edge_order(edge_order), n_kicks(n_kicks),
    integration_order(integration_order), edgeFlags(edgeFlags), d_rho0(d_rho0),
    d_rho_actual(d_rho_actual), d_rad_coef(d_rad_coef),
    d_isrConstant(d_isrConstant), d_expansionOrder(d_expansionOrder),
    d_meanPhotonsPerMeter0(d_meanPhotonsPerMeter0),
    d_normalizedCriticalEnergy0(d_normalizedCriticalEnergy0),
    d_distributionBasedRadiation(d_distributionBasedRadiation),
    d_includeOpeningAngle(d_includeOpeningAngle),
    srGaussianLimit(srGaussianLimit), d_refTrajectoryData(d_refTrajectoryData),
    d_refTrajectoryMode(d_refTrajectoryMode) {};

  __device__ unsigned int operator()(gpuParticleAccessor& coord){
    unsigned int tid = coord.getParticleIndex();

    double rho, s, x, xp, y, yp, dp, dp0, delta_xp;
    double Fx, Fy, dp_prime;
    double Qi[6], Qf[6];

    coord[4] += dzi*EXSQRT(1 + sqr(coord[1]) + sqr(coord[3]), sqrtOrder);
    coord[0]  = coord[0] + dxi + dzi*coord[1];
    coord[2]  = coord[2] + dyi + dzi*coord[3];

    x  =  coord[0]*cos_ttilt + coord[2]*sin_ttilt;
    y  = -coord[0]*sin_ttilt + coord[2]*cos_ttilt;
    xp =  coord[1]*cos_ttilt + coord[3]*sin_ttilt;
    yp = -coord[1]*sin_ttilt + coord[3]*cos_ttilt;
    s  = coord[4];
    dp = dp0 = coord[5];

    if (edgeFlags&BEND_EDGE1_EFFECTS) {
      rho = (1+dp)*d_rho_actual;
      if (edge_order<=1 || edge1_effects==1) {
        /* apply edge focusing */
        delta_xp = tan(e1)/rho*x;
        if (e1_kick_limit>0 && fabs(delta_xp)>e1_kick_limit)
          delta_xp = SIGN(delta_xp)*e1_kick_limit;
        xp += delta_xp;
        yp -= tan(e1-psi1/(1+dp))/rho*y;
      } else if (edge_order>=2 && edge1_effects==1) {
        gpu_apply_edge_effects(&x, &xp, &y, &yp, rho, n, e1, he1, psi1*(1+dp), -1);
      } else if (edge1_effects==2) {
        rho = (1+dp)*d_rho_actual;
        gpu_dipoleFringeSym(&x, &xp, &y, &yp, &dp, d_rho_actual, -1., edge_order, b0/d_rho0, e1, 2*hgap, fint, h1);
      }
    }

    /* load input coordinates into arrays */
    Qi[0] = x;  Qi[1] = xp;  Qi[2] = y;  Qi[3] = yp;  Qi[4] = 0;  Qi[5] = dp;

    if (edgeFlags&BEND_EDGE1_EFFECTS && e1!=0 && d_rad_coef) {
      /* pre-adjust dp/p to anticipate error made by integrating over entire sector */
      gpu_computeCSBENDFields(&Fx, &Fy, x, y, d_expansionOrder);

      dp_prime = -d_rad_coef*(sqr(Fx)+sqr(Fy))*sqr(1+dp)*EXSQRT(sqr(1+x/d_rho0)+sqr(xp)+sqr(yp), sqrtOrder);
      Qi[5] -= dp_prime*x*tan(e1);
    }

    gpu_convertToDipoleCanonicalCoordinates(Qi, d_rho0, sqrtOrder);
 
    int particle_lost = 0;
    double s_lost;
    double *tSigmaDelta2 = NULL;
    if (d_sigmaDelta2) tSigmaDelta2 = &d_sigmaDelta2[tid];
    double *t_gauss_rn = NULL;
    if (d_gauss_rn) t_gauss_rn = &d_gauss_rn[tid];
    curandState_t* t_state = NULL;
    if (state) t_state = &state[tid];
    if (integration_order==4)
      gpu_integrate_csbend_ord4(Qf, Qi, tSigmaDelta2, length, 
          n_kicks, sqrtOrder, d_rho0, Po, particle_lost, s_lost,
          d_rho_actual, d_rad_coef, d_isrConstant, 
          d_distributionBasedRadiation, d_includeOpeningAngle,
          d_meanPhotonsPerMeter0, d_normalizedCriticalEnergy0,
          d_expansionOrder, t_gauss_rn, t_state, srGaussianLimit,
          d_refTrajectoryData, d_refTrajectoryMode);
    else
      gpu_integrate_csbend_ord2(Qf, Qi, tSigmaDelta2, length,
          n_kicks, sqrtOrder, d_rho0, Po, particle_lost, s_lost,
          d_rho_actual, d_rad_coef, d_isrConstant, 
          d_distributionBasedRadiation, d_includeOpeningAngle,
          d_meanPhotonsPerMeter0, d_normalizedCriticalEnergy0,
          d_expansionOrder, t_gauss_rn, t_state, srGaussianLimit,
          d_refTrajectoryData, d_refTrajectoryMode);

    gpu_convertFromDipoleCanonicalCoordinates(Qf, d_rho0, sqrtOrder);

    if (particle_lost) {
      // TODO: accepted particles
      coord[4] = z_start + s_lost;
      coord[5] = Po*(1+coord[5]);
      d_sortIndex[tid] = tid+coord.getParticlePitch();
      return 0;
    }

    if (edgeFlags&BEND_EDGE2_EFFECTS && e2!=0 && d_rad_coef) {
      /* post-adjust dp/p to correct error made by integrating over entire sector */
      x = Qf[0];
      xp = Qf[1];
      y = Qf[2];
      yp = Qf[3];
      dp = Qf[5];

      gpu_computeCSBENDFields(&Fx, &Fy, x, y, d_expansionOrder);

      dp_prime = -d_rad_coef*(sqr(Fx)+sqr(Fy))*sqr(1+dp)*EXSQRT(sqr(1+x/d_rho0)+sqr(xp)+sqr(yp), sqrtOrder);
      Qf[5] -= dp_prime*x*tan(e2);
    }

    /* get final coordinates */
    if (d_rad_coef || d_isrConstant) {
      double p0, p1;
      double beta0, beta1;
      /* fix previous distance information to reflect new velocity--since distance
       * is really time-of-flight at the current velocity 
       */
      p0 = Po*(1+dp0);
      beta0 = p0/sqrt(sqr(p0)+1);
      p1 = Po*(1+Qf[5]);
      beta1 = p1/sqrt(sqr(p1)+1);
      s = beta1*s/beta0 + Qf[4];
    }
    else
      s += Qf[4];
    x = Qf[0];  xp = Qf[1];  y = Qf[2];  yp = Qf[3];  dp = Qf[5];

    if (edgeFlags&BEND_EDGE2_EFFECTS) {
      /* apply edge focusing */
      rho = (1+dp)*d_rho_actual;
      if (edge_order<=1 || edge2_effects==1) {
        delta_xp = tan(e2)/rho*x;
        if (e2_kick_limit>0 && fabs(delta_xp)>e2_kick_limit)
          delta_xp = SIGN(delta_xp)*e2_kick_limit;
        xp += delta_xp;
        yp -= tan(e2-psi2/(1+dp))/rho*y;
      } else if (edge_order>=2 && edge2_effects==1) {
        gpu_apply_edge_effects(&x, &xp, &y, &yp, rho, n, e2, he2, psi2*(1+dp), 1);
      } else if (edge2_effects==2) {
        //rho = (1+dp)*d_rho_actual;
        gpu_dipoleFringeSym(&x, &xp, &y, &yp, &dp, d_rho_actual, 1., edge_order, b0/d_rho0, e2, 2*hgap, fint, h2);
      }
    }
    
    coord[0] =  x*cos_ttilt -  y*sin_ttilt + dcet0; // dcoord_etilt[0];
    coord[2] =  x*sin_ttilt +  y*cos_ttilt + dcet2; // dcoord_etilt[2];
    coord[1] = xp*cos_ttilt - yp*sin_ttilt + dcet1; // dcoord_etilt[1];
    coord[3] = xp*sin_ttilt + yp*cos_ttilt + dcet3; // dcoord_etilt[3];
    coord[4] = s;
    coord[5] = dp;

    coord[0] += dxf + dzf*coord[1];
    coord[2] += dyf + dzf*coord[3];
    coord[4] += dzf*EXSQRT(1+ sqr(coord[1]) + sqr(coord[3]), sqrtOrder);

    d_sortIndex[tid] = tid;
    return 1;
  }
};

extern "C" {

long gpu_track_through_csbend(long n_part, CSBEND *csbend, 
       double p_error, double Po, double **accepted, double z_start, 
       double *sigmaDelta2, char *rootname)
{
  double h;
  double n, fse;
  double tilt, etilt, cos_ttilt, sin_ttilt, ttilt;
  double angle, e1, e2, Kg;
  double psi1, psi2, he1, he2;
  double dcoord_etilt[6];
  double dxi, dyi, dzi;
  double dxf, dyf, dzf;
  double e1_kick_limit, e2_kick_limit;
  static long largeRhoWarning = 0;

  struct GPUBASE* gpuBase = getGpuBase();

  if (!csbend)
    bombElegant("null CSBEND pointer (track_through_csbend)", NULL);

  if (csbend->edge_order>1 && (csbend->edge1_effects==2 || csbend->edge2_effects==2) && csbend->hgap==0)
    bombElegant("CSBEND has EDGE_ORDER>1 and EDGE[12]_EFFECTS==2, but HGAP=0. This gives undefined results.", NULL);
  
  if (csbend->referenceCorrection) {
    if (csbend->refTrajectoryChangeSet==0 || csbend->refLength!=csbend->length || csbend->refAngle!=csbend->angle || csbend->refKicks!=csbend->n_kicks) {
      /* Figure out the reference trajectory offsets to suppress inaccuracy in the integrator */
      CSBEND csbend0;
      double **part0;
      TRACKING_CONTEXT tcontext;

      getTrackingContext(&tcontext);
      if (tcontext.elementOccurrence>0) {
	printf("Determining reference trajectory for CSBEND %s#%ld at s=%e\n", tcontext.elementName, tcontext.elementOccurrence, tcontext.zStart);
      }
      
      if (csbend->refTrajectoryChange && csbend->refKicks) {
        free_czarray_2d((void**)csbend->refTrajectoryChange, csbend->refKicks, 5);
        csbend->refTrajectoryChange = NULL;
        csbend->refKicks = 0;
      }
      
      part0 = (double**)czarray_2d(sizeof(double), 1, 7);
      memset(part0[0], 0, sizeof(**part0)*7);
      memcpy(&csbend0, csbend, sizeof(*csbend));
      csbend0.dx = csbend0.dy = csbend0.dz = csbend0.fse = csbend0.etilt = csbend0.isr = csbend0.synch_rad = 0;
      
      csbend0.refTrajectoryChange = csbend->refTrajectoryChange = (double**)czarray_2d(sizeof(double), csbend->n_kicks, 5);
      refTrajectoryPoints = csbend->n_kicks;
      csbend0.refLength = csbend0.length;
      csbend0.refAngle = csbend0.angle;
      csbend0.refKicks = csbend0.n_kicks;
      /* This forces us into the next branch on the next call to this routine */
      csbend0.refTrajectoryChangeSet = 1;
      setTrackingContext((char*)"csbend0", 0, T_CSBEND, (char*)"none");
      // keep single particle csbend on CPU
      gpuBase->elementOnGpu=0;
      track_through_csbend(part0, 1, &csbend0, p_error, Po, NULL, 0, NULL, rootname);
      gpuBase->elementOnGpu=1;
      csbend->refTrajectoryChangeSet = 2;  /* indicates that reference trajectory has been determined */

      csbend->refKicks = csbend->n_kicks;
      csbend->refLength = csbend->length;
      csbend->refAngle = csbend->angle;
      free_czarray_2d((void**)part0, 1, 7);

      refTrajectoryData = csbend->refTrajectoryChange;
      refTrajectoryPoints = csbend->refKicks;
      refTrajectoryMode = SUBTRACT_TRAJECTORY;
    } else if (csbend->refTrajectoryChangeSet==1) {
      /* indicates reference trajectory is about to be determined */
      refTrajectoryData = csbend->refTrajectoryChange;
      refTrajectoryPoints = csbend->n_kicks;
      refTrajectoryMode = RECORD_TRAJECTORY;
      csbend->refTrajectoryChangeSet = 2;
    } else {
      /* assume that reference trajectory already determined */
      refTrajectoryData = csbend->refTrajectoryChange;
      refTrajectoryPoints = csbend->refKicks;
      refTrajectoryMode = SUBTRACT_TRAJECTORY;
    }
  } else
    refTrajectoryMode = 0;

  if (csbend->angle==0) {
    gpu_exactDrift(n_part, csbend->length);
    return n_part;
  }
  
  if (!(csbend->edgeFlags&BEND_EDGE_DETERMINED)) 
    bombElegant("CSBEND element doesn't have edge flags set.", NULL);
  
  if (csbend->integration_order!=2 && csbend->integration_order!=4)
    bombElegant("CSBEND integration_order is invalid--must be either 2 or 4", NULL);

  rho0 =  csbend->length/csbend->angle;
  if (csbend->use_bn) {
    csbend->b[0] = csbend->b1;
    csbend->b[1] = csbend->b2;
    csbend->b[2] = csbend->b3;
    csbend->b[3] = csbend->b4;
    csbend->b[4] = csbend->b5;
    csbend->b[5] = csbend->b6;
    csbend->b[6] = csbend->b7;
    csbend->b[7] = csbend->b8;
  } else {
    csbend->b[0] = csbend->k1*rho0;
    csbend->b[1] = csbend->k2*rho0;
    csbend->b[2] = csbend->k3*rho0;
    csbend->b[3] = csbend->k4*rho0;
    csbend->b[4] = csbend->k5*rho0;
    csbend->b[5] = csbend->k6*rho0;
    csbend->b[6] = csbend->k7*rho0;
    csbend->b[7] = csbend->k8*rho0;
  }
  if (csbend->xReference>0) {
    double term = 1/csbend->xReference, f[8];
    long i;
    f[0] = csbend->f1;
    f[1] = csbend->f2;
    f[2] = csbend->f3;
    f[3] = csbend->f4;
    f[4] = csbend->f5;
    f[5] = csbend->f6;
    f[6] = csbend->f7;
    f[7] = csbend->f8;
    for (i=0; i<=8; i++) {
      csbend->b[i] += f[i]*term;
      term *= (i+2)/csbend->xReference;
    }
  }
  
  he1 = csbend->h1;
  he2 = csbend->h2;
  if (csbend->angle<0) {
    long i;
    angle = -csbend->angle;
    e1    = -csbend->e1;
    e2    = -csbend->e2;
    etilt = csbend->etilt;
    tilt  = csbend->tilt + PI;      /* work in rotated system */
    rho0  = csbend->length/angle;
    for (i=0; i<8; i+=2)
      csbend->b[i] *= -1;
  }
  else {
    angle = csbend->angle;
    e1    = csbend->e1;
    e2    = csbend->e2;
    etilt = csbend->etilt;
    tilt  = csbend->tilt;
    rho0  = csbend->length/angle;
  }


  if (rho0>1e6) {
    if (csbend->k2!=0)
      bombElegant("Error: One or more CSBENDs have radius > 1e6 but non-zero K2. Best to convert this to KQUSE or KSEXT.\n", NULL);
    if (csbend->k1!=0) {
      ELEMENT_LIST elem;
      KQUAD kquad;
      static short largeRhoWarningK1 = 0;
      if (!largeRhoWarningK1) {
#if USE_MPI
	if (myid==1)
	  dup2(fd,fileno(stdout)); /* Let the first slave processor write the output */
#endif
        printf("Warning: One or more CSBENDs have radius > 1e6 but non-zero K1.  Treated as KQUAD.\n");
        printf("*** All higher multipoles are ignored for these elements!\n");
        largeRhoWarningK1 = 1;
#if USE_MPI
	if (myid==1) {
#if defined(_WIN32)
	  freopen("NUL","w",stdout); 
#else
	  freopen("/dev/null","w",stdout); 
#endif
	}
#endif  
      }
      memset(&elem, 0, sizeof(elem));
      memset(&kquad, 0, sizeof(kquad));
      elem.p_elem = (char*)&kquad;
      elem.type = T_KQUAD;
      kquad.length = csbend->length;
      kquad.k1 = csbend->k1;
      kquad.tilt = csbend->tilt;
      kquad.dx = csbend->dx;
      kquad.dy = csbend->dy;
      kquad.dz = csbend->dz;
      kquad.synch_rad = csbend->synch_rad;
      kquad.isr = csbend->isr;
      kquad.isr1Particle = csbend->isr1Particle;
      kquad.n_kicks = csbend->n_kicks;
      kquad.integration_order = csbend->integration_order;
      return gpu_multipole_tracking2(n_part, &elem, p_error, Po, accepted, z_start, 1, 1, 0, NULL, sigmaDelta2);
    } else {
      if (!largeRhoWarning) {
#if USE_MPI
	if (myid==1)
	  dup2(fd,fileno(stdout)); /* Let the first slave processor write the output */
#endif
        printf("Warning: One or more CSBENDs have radius > 1e6.  Treated as EDRIFT.\n");
        printf("*** All higher multipoles are ignored for these elements!\n");
#if USE_MPI
	if (myid==1) {
#if defined(_WIN32)
	  freopen("NUL","w",stdout); 
#else
	  freopen("/dev/null","w",stdout); 
#endif
	}
#endif  
        largeRhoWarning = 1;
      }
      gpu_exactDrift(n_part, csbend->length);
      return n_part;
    }
  }


  fse = csbend->fse;
  h = 1/rho0;
  n = -csbend->b[0]/h;
  if (fse>-1)
    rho_actual = 1/((1+fse)*h);
  else
    rho_actual = 1e16/h;

  e1_kick_limit = csbend->edge1_kick_limit;
  e2_kick_limit = csbend->edge2_kick_limit;
  if (csbend->kick_limit_scaling) {
    e1_kick_limit *= rho0/rho_actual;
    e2_kick_limit *= rho0/rho_actual;
  }
  if (e1_kick_limit>0 || e2_kick_limit>0)
    fprintf(stdout, "rho0=%e  rho_a=%e fse=%e e1_kick_limit=%e e2_kick_limit=%e\n",
            rho0, rho_actual, csbend->fse, e1_kick_limit, e2_kick_limit);
    fflush(stdout);
  
  /* angles for fringe-field effects */
  Kg   = 2*csbend->hgap*csbend->fint;
  psi1 = Kg/rho_actual/cos(e1)*(1+sqr(sin(e1)));
  psi2 = Kg/rho_actual/cos(e2)*(1+sqr(sin(e2)));

  /* rad_coef is d((P-Po)/Po)/ds for the on-axis, on-momentum particle, where po is the momentum of
   * the central particle.
   */
  if (csbend->synch_rad)
    rad_coef = sqr(particleCharge)*pow3(Po)*sqr(1+fse)/(6*PI*epsilon_o*sqr(c_mks)*particleMass*sqr(rho0));
  else
    rad_coef = 0;
  /* isrConstant is the RMS increase in dP/P per meter due to incoherent SR.  */
  isrConstant = particleRadius*sqrt(55.0/(24*sqrt(3))*pow5(Po)*
                            137.0359895/pow3(fabs(rho_actual)));
  if (!csbend->isr || (csbend->isr1Particle==0 && n_part==1))
    /* Minus sign here indicates that we accumulate ISR into sigmaDelta^2 but don't apply it to particles. */
    isrConstant *= -1; 

  if ((distributionBasedRadiation = csbend->distributionBasedRadiation)) {
    /* Sands 5.15 */
    meanPhotonsPerRadian0 = 5.0/(2.0*sqrt(3))*Po/137.0359895;  
    meanPhotonsPerMeter0 = (5*c_mks*Po*particleMass*particleRadius)/(2*sqrt(3)*hbar_mks*rho_actual);
    /* Critical energy normalized to beam energy, Sands 5.9 */
    normalizedCriticalEnergy0 = 3.0/2*hbar_mks*c_mks*pow3(Po)/fabs(rho_actual)/(Po*particleMass*sqr(c_mks));
    /* fprintf(stderr, "Mean photons per radian expected: %le   ECritical/E: %le\n", 
            meanPhotonsPerRadian0, normalizedCriticalEnergy0);
    */
    includeOpeningAngle = csbend->includeOpeningAngle;
  }
    
  computeCSBENDFieldCoefficients(csbend->b, h, csbend->nonlinear, csbend->expansionOrder);

  ttilt = tilt + etilt;
  if (ttilt==0) {
    cos_ttilt = 1;
    sin_ttilt = 0;
  }
  else if (fabs(fabs(ttilt)-PI)<1e-12) {
    cos_ttilt = -1;
    sin_ttilt = 0;
  }
  else if (fabs(ttilt-PIo2)<1e-12) {
    cos_ttilt = 0;
    sin_ttilt = 1;
  }
  else if (fabs(ttilt+PIo2)<1e-12) {
    cos_ttilt = 0;
    sin_ttilt = -1;
  }
  else {
    cos_ttilt = cos(ttilt);
    sin_ttilt = sin(ttilt);
  }

  computeEtiltCentroidOffset(dcoord_etilt, rho0, angle, etilt, tilt);

  dxi = -csbend->dx;
  dzi =  csbend->dz;
  dyi = -csbend->dy;
  
  /* must use the original angle here because the translation is done after
   * the final rotation back
   */
  dxf =  csbend->dx*cos(csbend->angle) + csbend->dz*sin(csbend->angle);
  dzf =  csbend->dx*sin(csbend->angle) - csbend->dz*cos(csbend->angle);
  dyf = csbend->dy;

#if !defined(PARALLEL)
  multipoleKicksDone += n_part*csbend->n_kicks*(csbend->integration_order==4?4:1);
#endif

  if (sigmaDelta2)
    *sigmaDelta2 = 0;
  double* d_sigmaDelta2 = NULL;

  if (isSlave || !notSinglePart) {
    unsigned int particlePitch = gpuBase->gpu_array_pitch;
    unsigned int* d_sortIndex = gpuBase->d_tempu_alpha;

    /* Copy Fx_xy, Fy_xy and table coefficients if needed */
    double Fall_xy[2*xOrderMax2];
    for (int ii=0; ii<expansionOrderMax; ii++) {
      for (int jj=0; jj<expansionOrderMax; jj++) {
        Fall_xy[ii*expansionOrderMax+jj]            = Fx_xy[ii][jj];
        Fall_xy[ii*expansionOrderMax+jj+xOrderMax2] = Fy_xy[ii][jj];
      }
    }
    cudaMemcpyToSymbol(c_Fx_xy, Fall_xy, sizeof(double)*xOrderMax2,
                       0, cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol(c_Fy_xy, &(Fall_xy[xOrderMax2]), sizeof(double)*xOrderMax2, 
                       0, cudaMemcpyHostToDevice);
    gpu_initializeInterpolationTables();

    /* d_temp_particles conflicts with killParticles */
    if (sigmaDelta2)
      cudaMalloc( (void**)&d_sigmaDelta2, sizeof(double)*particlePitch);

    /* initalize random numbers as needed. */
    curandState_t* state = NULL;
    double* d_gauss_rn =  gpuBase->d_temp_particles + particlePitch;
    if (rad_coef || isrConstant)
      if (!distributionBasedRadiation && isrConstant)
        gpu_d_gauss_rn_lim(d_gauss_rn, n_part, 0.0, 1.0, 3.0, random_2(0));
      else if (distributionBasedRadiation)
        state =
          (curandState_t*)gpu_get_rand_state(d_gauss_rn, n_part, random_2(0));
    
    /* copy trajectory data if needed */
    double* d_refTrajectoryData =  gpuBase->d_temp_particles + 2*particlePitch;
    if (refTrajectoryMode) {
      cudaMemcpy(refTrajectoryData[0], d_refTrajectoryData, 
                 5*refTrajectoryPoints, cudaMemcpyHostToDevice);
    }

    n_part = killParticles(n_part, d_sortIndex, accepted,
      gpu_track_through_csbend_kernel(d_sortIndex, d_sigmaDelta2, 
        d_gauss_rn, state,
        Po, z_start, dxi, dyi, dzi, dxf, dyf, dzf, cos_ttilt, 
        sin_ttilt, e1, e2, he1, he2, psi1, psi2, n, dcoord_etilt[0],
        dcoord_etilt[1], dcoord_etilt[2], dcoord_etilt[3], 
        e1_kick_limit, e2_kick_limit, csbend->sqrtOrder, 
        csbend->length, csbend->b[0], csbend->hgap, csbend->fint, 
        csbend->h1, csbend->h2, csbend->edge1_effects, csbend->edge2_effects, 
        csbend->edge_order, csbend->n_kicks, csbend->integration_order,
        csbend->edgeFlags, rho0, rho_actual, rad_coef, isrConstant,
        expansionOrder1, meanPhotonsPerMeter0, normalizedCriticalEnergy0,
        distributionBasedRadiation, includeOpeningAngle, srGaussianLimit,
        d_refTrajectoryData, refTrajectoryMode)); 
    gpuErrorHandler("gpu_track_through_csbend: gpu_track_through_csbend_kernel");
  }
  if (distributionBasedRadiation) {
    radiansTotal += fabs(csbend->angle);
    distributionBasedRadiation = 0;
  }

  if (sigmaDelta2) {
    *sigmaDelta2 = gpuReduceAdd(d_sigmaDelta2, n_part);
    /* Return average value for all particles */
    *sigmaDelta2 /= n_part;
    cudaFree(d_sigmaDelta2);
  }

  return(n_part);
}

} // extern "C"

__device__ void 
gpu_integrate_csbend_ord2(double *Qf, double *Qi, double *sigmaDelta2, 
    double s, int n, int sqrtOrder, double rho0, 
    double p0, int &particle_lost, double &s_lost,
    double d_rho_actual, double d_rad_coef, double d_isrConstant, 
    int d_distributionBasedRadiation, int d_includeOpeningAngle,
    double d_meanPhotonsPerMeter0, double d_normalizedCriticalEnergy0,
    int d_expansionOrder, double *d_gauss_rn, curandState_t *state,
    double srGaussianLimit,  double *d_refTrajectoryData, double d_refTrajectoryMode)
{
  int i;
  double factor, f, phi, ds, dsh, dist;
  double Fx, Fy, x, y;
  double sine, cosi, tang;
  double sin_phi, cos_phi;
  
#define X0 Qi[0]
#define XP0 Qi[1]
#define Y0 Qi[2]
#define YP0 Qi[3]
#define S0 Qi[4]
#define DPoP0 Qi[5]

#define X Qf[0]
#define QX Qf[1]
#define Y Qf[2]
#define QY Qf[3]
#define S Qf[4]
#define DPoP Qf[5]

//  if (!Qf)
//    bombElegant("NULL final coordinates pointer ()", NULL);
//  if (!Qi)
//    bombElegant("NULL initial coordinates pointer (integrate_csbend_ord2)", NULL);
//  if (n<1)
//    bombElegant("invalid number of steps (integrate_csbend_ord2)", NULL);

  /* calculate canonical momenta (scaled to central momentum) */
/*  dp = DPoP0;
  f = (1+dp)/EXSQRT(sqr(1+X0/rho0) + sqr(XP0) + sqr(YP0), sqrtOrder);
  QX = XP0*f;
  QY = YP0*f;
*/

  memcpy(Qf, Qi, sizeof(*Qi)*6);
  
  X = X0;
  Y = Y0;
  S = S0;
  DPoP = DPoP0;

  ds = s/n;
  dsh = ds/2;
  dist = 0;
  
  for (i=0; i<n; i++) {
    if (i==0) {
      /* do half-length drift */
      if ((f=sqr(1+DPoP)-sqr(QY))<=0) {
        particle_lost = 1;
        s_lost = dist;
        return;
      }
      f = EXSQRT(f, sqrtOrder);
      if (fabs(QX/f)>1) {
        particle_lost = 1;
        s_lost = dist;
        return;
      }
      phi = asin(sin_phi=QX/f);
      sine = sin(dsh/rho0+phi);
      if ((cosi = cos(dsh/rho0+phi))==0) {
        particle_lost = 1;
        s_lost = dist;
        return;
      }
      tang = sine/cosi;
      cos_phi = cos(phi);
      QX = f*sine;
      Y += QY*(factor=(rho0+X)*cos_phi/f*(tang-sin_phi/cos_phi));
      dist += factor*(1+DPoP);
      f = cos_phi/cosi;
      X  = rho0*(f-1) + f*X;
    }

    /* calculate the scaled fields */
    x = X;
    y = Y;

    gpu_computeCSBENDFields(&Fx, &Fy, x, y, d_expansionOrder);

    /* do kicks */
    QX += -ds*(1+X/rho0)*Fy/d_rho_actual;
    QY += ds*(1+X/rho0)*Fx/d_rho_actual;
    if (d_rad_coef || d_isrConstant)
      gpu_addRadiationKick(&QX, &QY, &DPoP, sigmaDelta2, sqrtOrder, 
		       X, 1./rho0, Fx, Fy, 
		       ds, d_rad_coef, ds, d_isrConstant, 
                       d_distributionBasedRadiation, d_includeOpeningAngle,
                       d_meanPhotonsPerMeter0, d_normalizedCriticalEnergy0, 
                       p0, d_gauss_rn, state, srGaussianLimit);
    
    if (i==n-1) {
      /* do half-length drift */
      if ((f=sqr(1+DPoP)-sqr(QY))<=0) {
        particle_lost = 1;
        s_lost = dist;
        return;
      }
      f = EXSQRT(f, sqrtOrder);
      if (fabs(QX/f)>1) {
        particle_lost = 1;
        s_lost = dist;
        return;
      }
      phi = asin(sin_phi=QX/f);
      sine = sin(dsh/rho0+phi);
      if ((cosi = cos(dsh/rho0+phi))==0) {
        particle_lost = 1;
        s_lost = dist;
        return;
      }
      tang = sine/cosi;
      cos_phi = cos(phi);
      QX = f*sine;
      Y += QY*(factor=(rho0+X)*cos_phi/f*(tang-sin_phi/cos_phi));
      dist += factor*(1+DPoP);
      f = cos_phi/cosi;
      X  = rho0*(f-1) + f*X;
    }
    else {
      /* do full-length drift */
      if ((f=sqr(1+DPoP)-sqr(QY))<=0) {
        particle_lost = 1;
        s_lost = dist;
        return;
      }
      f = EXSQRT(f, sqrtOrder);
      if (fabs(QX/f)>1) {
        particle_lost = 1;
        s_lost = dist;
        return;
      }
      phi = asin(sin_phi=QX/f);
      sine = sin(ds/rho0+phi);
      if ((cosi = cos(ds/rho0+phi))==0) {
        particle_lost = 1;
        s_lost = dist;
        return;
      }
      tang = sine/cosi;
      cos_phi = cos(phi);
      QX = f*sine;
      Y += QY*(factor=(rho0+X)*cos_phi/f*(tang-sin_phi/cos_phi));
      dist += factor*(1+DPoP);
      f = cos_phi/cosi;
      X  = rho0*(f-1) + f*X;
    }

    if (d_refTrajectoryMode==RECORD_TRAJECTORY) {
      d_refTrajectoryData[i*5+0] = X;
      d_refTrajectoryData[i*5+1] = QX;
      d_refTrajectoryData[i*5+2] = Y;
      d_refTrajectoryData[i*5+3] = QY;
      d_refTrajectoryData[i*5+4] = dist - ds;
      X = QX = Y = QY = dist = 0;
    }
    if (d_refTrajectoryMode==SUBTRACT_TRAJECTORY) {
      X -= d_refTrajectoryData[i*5+0];
      QX -= d_refTrajectoryData[i*5+1];
      Y -= d_refTrajectoryData[i*5+2];
      QY -= d_refTrajectoryData[i*5+3];
      dist -= d_refTrajectoryData[i*5+4];
    }
  }

  /* convert back to slopes */
/*
  f = (1+X/rho0)/EXSQRT(sqr(1+DPoP)-sqr(QX)-sqr(QY), sqrtOrder);
  Qf[1] *= f;
  Qf[3] *= f;
*/

  Qf[4] += dist;
}

__device__ void 
gpu_integrate_csbend_ord4(double *Qf, double *Qi, double *sigmaDelta2, 
    double s, int n, int sqrtOrder, double rho0, double p0,
    int &particle_lost, double &s_lost,
    double d_rho_actual, double d_rad_coef, double d_isrConstant, 
    int d_distributionBasedRadiation, int d_includeOpeningAngle,
    double d_meanPhotonsPerMeter0, double d_normalizedCriticalEnergy0,
    int d_expansionOrder, double *d_gauss_rn, curandState_t *state,
    double srGaussianLimit, double *d_refTrajectoryData, double d_refTrajectoryMode)
{
  int i;
  double factor, f, phi, ds, dsh, dist;
  double Fx, Fy, x, y;
  double sine, cosi, tang;
  double sin_phi, cos_phi;
  
#define X0 Qi[0]
#define XP0 Qi[1]
#define Y0 Qi[2]
#define YP0 Qi[3]
#define S0 Qi[4]
#define DPoP0 Qi[5]

#define X Qf[0]
#define QX Qf[1]
#define Y Qf[2]
#define QY Qf[3]
#define S Qf[4]
#define DPoP Qf[5]

  /* BETA is 2^(1/3) */
#define BETA 1.25992104989487316477

//  if (!Qf)
//    bombElegant("NULL final coordinates pointer ()", NULL);
//  if (!Qi)
//    bombElegant("NULL initial coordinates pointer (integrate_csbend_ord4)", NULL);
//  if (n<1)
//    bombElegant("invalid number of steps (integrate_csbend_ord4)", NULL);

  /* calculate canonical momenta (scaled to central momentum) */
/*
  dp = DPoP0;
  f = (1+dp)/EXSQRT(sqr(1+X0/rho0) + sqr(XP0) + sqr(YP0), sqrtOrder);
  QX = XP0*f;
  QY = YP0*f;

  X = X0;
  Y = Y0;
  S = S0;
  DPoP = DPoP0;
*/
  
  memcpy(Qf, Qi, sizeof(*Qi)*6);

  dist = 0;

  s /= n;
  for (i=0; i<n; i++) {
    
    /* do first drift */
    dsh = s/2/(2-BETA);
    if ((f=sqr(1+DPoP)-sqr(QY))<=0) {
      particle_lost = 1;
      s_lost = dist;
      return;
    }
    f = EXSQRT(f, sqrtOrder);
    if (fabs(QX/f)>1) {
      particle_lost = 1;
      s_lost = dist;
      return;
    }
    phi = asin(sin_phi=QX/f);
    sine = sin(dsh/rho0+phi);
    if ((cosi = cos(dsh/rho0+phi))==0) {
      particle_lost = 1;
      s_lost = dist;
      return;
    }
    tang = sine/cosi;
    cos_phi = cos(phi);
    QX = f*sine;
    Y += QY*(factor=(rho0+X)*cos_phi/f*(tang-sin_phi/cos_phi));
    dist += factor*(1+DPoP);
    f = cos_phi/cosi;
    X  = rho0*(f-1) + f*X;
 
    /* do first kick */
    ds = s/(2-BETA);
    /* -- calculate the scaled fields */
    x = X;
    y = Y;

    gpu_computeCSBENDFields(&Fx, &Fy, x, y, d_expansionOrder);
    
    /* --do kicks */
    QX += -ds*(1+X/rho0)*Fy/d_rho_actual;
    QY += ds*(1+X/rho0)*Fx/d_rho_actual;
    if (d_rad_coef || d_isrConstant) {
      gpu_addRadiationKick(&QX, &QY, &DPoP, sigmaDelta2, sqrtOrder,
		       X, 1./rho0, Fx, Fy, 
		       ds, d_rad_coef, s/3, d_isrConstant,
                       d_distributionBasedRadiation, d_includeOpeningAngle,
                       d_meanPhotonsPerMeter0, d_normalizedCriticalEnergy0,
                       p0, d_gauss_rn, state, srGaussianLimit);
    }

    /* do second drift */
    dsh = s*(1-BETA)/(2-BETA)/2;
    if ((f=sqr(1+DPoP)-sqr(QY))<=0) {
      particle_lost = 1;
      s_lost = dist;
      return;
    }
    f = EXSQRT(f, sqrtOrder);
    if (fabs(QX/f)>1) {
      particle_lost = 1;
      s_lost = dist;
      return;
    }
    phi = asin(sin_phi=QX/f);
    sine = sin(dsh/rho0+phi);
    if ((cosi = cos(dsh/rho0+phi))==0) {
      particle_lost = 1;
      s_lost = dist;
      return;
    }
    tang = sine/cosi;
    cos_phi = cos(phi);
    QX = f*sine;
    Y += QY*(factor=(rho0+X)*cos_phi/f*(tang-sin_phi/cos_phi));
    dist += factor*(1+DPoP);
    f = cos_phi/cosi;
    X  = rho0*(f-1) + f*X;
    
    /* do second kick */
    ds = -s*BETA/(2-BETA);
    /* -- calculate the scaled fields */
    x = X;
    y = Y;
    gpu_computeCSBENDFields(&Fx, &Fy, x, y, d_expansionOrder);

    /* --do kicks */
    QX += -ds*(1+X/rho0)*Fy/d_rho_actual;
    QY += ds*(1+X/rho0)*Fx/d_rho_actual;
    if (d_rad_coef || d_isrConstant)
      gpu_addRadiationKick(&QX, &QY, &DPoP, sigmaDelta2, sqrtOrder,
		       X, 1./rho0, Fx, Fy, 
		       ds, d_rad_coef, s/3, d_isrConstant,
                       d_distributionBasedRadiation, d_includeOpeningAngle,
                       d_meanPhotonsPerMeter0, d_normalizedCriticalEnergy0,
                       p0, d_gauss_rn, state, srGaussianLimit);

    /* do third drift */
    dsh = s*(1-BETA)/(2-BETA)/2;
    if ((f=sqr(1+DPoP)-sqr(QY))<=0) {
      particle_lost = 1;
      s_lost = dist;
      return;
    }
    f = EXSQRT(f, sqrtOrder);
    if (fabs(QX/f)>1) {
      particle_lost = 1;
      return;
    }
    phi = asin(sin_phi=QX/f);
    sine = sin(dsh/rho0+phi);
    if ((cosi = cos(dsh/rho0+phi))==0) {
      particle_lost = 1;
      s_lost = dist;
      return;
    }
    tang = sine/cosi;
    cos_phi = cos(phi);
    QX = f*sine;
    Y += QY*(factor=(rho0+X)*cos_phi/f*(tang-sin_phi/cos_phi));
    dist += factor*(1+DPoP);
    f = cos_phi/cosi;
    X  = rho0*(f-1) + f*X;
    
    /* do third kick */
    ds = s/(2-BETA);
    /* -- calculate the scaled fields */
    x = X;
    y = Y;
    gpu_computeCSBENDFields(&Fx, &Fy, x, y, d_expansionOrder);

    /* --do kicks */
    QX += -ds*(1+X/rho0)*Fy/d_rho_actual;
    QY += ds*(1+X/rho0)*Fx/d_rho_actual;
    if (d_rad_coef || d_isrConstant) 
      gpu_addRadiationKick(&QX, &QY, &DPoP, sigmaDelta2, sqrtOrder,
		       X, 1./rho0, Fx, Fy, 
		       ds, d_rad_coef, s/3, d_isrConstant,
                       d_distributionBasedRadiation, d_includeOpeningAngle,
                       d_meanPhotonsPerMeter0, d_normalizedCriticalEnergy0,
                       p0, d_gauss_rn, state, srGaussianLimit);
    
    /* do fourth drift */
    dsh = s/2/(2-BETA);
    if ((f=sqr(1+DPoP)-sqr(QY))<=0) {
      particle_lost = 1;
      s_lost = dist;
      return;
    }
    f = EXSQRT(f, sqrtOrder);
    if (fabs(QX/f)>1) {
      particle_lost = 1;
      s_lost = dist;
      return;
    }
    phi = asin(sin_phi=QX/f);
    sine = sin(dsh/rho0+phi);
    if ((cosi = cos(dsh/rho0+phi))==0) {
      particle_lost = 1;
      s_lost = dist;
      return;
    }
    tang = sine/cosi;
    cos_phi = cos(phi);
    QX = f*sine;
    Y += QY*(factor=(rho0+X)*cos_phi/f*(tang-sin_phi/cos_phi));
    dist += factor*(1+DPoP);
    f = cos_phi/cosi;
    X  = rho0*(f-1) + f*X;

    if (d_refTrajectoryMode==RECORD_TRAJECTORY) {
      d_refTrajectoryData[i*5+0] = X;
      d_refTrajectoryData[i*5+1] = QX;
      d_refTrajectoryData[i*5+2] = Y;
      d_refTrajectoryData[i*5+3] = QY;
      d_refTrajectoryData[i*5+4] = dist - s;
      X = QX = Y = QY = dist = 0;
    }
    if (d_refTrajectoryMode==SUBTRACT_TRAJECTORY) {
      X -= d_refTrajectoryData[i*5+0];
      QX -= d_refTrajectoryData[i*5+1];
      Y -= d_refTrajectoryData[i*5+2];
      QY -= d_refTrajectoryData[i*5+3];
      dist -= d_refTrajectoryData[i*5+4];
    }
  }

  /* convert back to slopes */
/*
  f = (1+X/rho0)/EXSQRT(sqr(1+DPoP)-sqr(QX)-sqr(QY), sqrtOrder);
  Qf[1] *= f;
  Qf[3] *= f;
*/

  Qf[4] += dist;
}

class gpu_track_through_csbendCSR_kernel1{
public:
  double length;
  gpu_track_through_csbendCSR_kernel1(double length) : length(length) {};

  __device__ inline void operator()(gpuParticleAccessor& part){

    part[0] += length*part[1];
    part[2] += length*part[3];
    part[4] += length;
  }
};

class gpu_track_through_csbendCSR_kernel2{
public:
  double *d_beta0, *d_x, *d_y;
  double Po;
  double dxi, dyi, dzi;
  double cos_ttilt, sin_ttilt;

  gpu_track_through_csbendCSR_kernel2(
    double *d_beta0, double *d_x, double *d_y, double Po, double dxi, 
    double dyi, double dzi, double cos_ttilt, double sin_ttilt) : 
    d_beta0(d_beta0), d_x(d_x), d_y(d_y), Po(Po), dxi(dxi), dyi(dyi), 
    dzi(dzi), cos_ttilt(cos_ttilt), sin_ttilt(sin_ttilt) {};

  /* check particle data, transform coordinates, and handle edge effects */
  __device__ inline void operator()(gpuParticleAccessor& coord){
    unsigned int tid = coord.getParticleIndex();

    double xp, yp, p0;

    /* adjust for element offsets */
    coord[4] += dzi*sqrt(1 + sqr(coord[1]) + sqr(coord[3]));
    coord[0]  = coord[0] + dxi + dzi*coord[1];
    coord[2]  = coord[2] + dyi + dzi*coord[3];
  
    /* perform tilt transformations and save some data */
    d_x[tid]  =  coord[0]*cos_ttilt + coord[2]*sin_ttilt;
    d_y[tid]  = -coord[0]*sin_ttilt + coord[2]*cos_ttilt;
    coord[0] = d_x[tid];
    coord[2] = d_y[tid];
    xp =  coord[1]*cos_ttilt + coord[3]*sin_ttilt;
    yp = -coord[1]*sin_ttilt + coord[3]*cos_ttilt;
    coord[1] = xp;
    coord[3] = yp;
    p0 = Po*(1+coord[5]);
    d_beta0[tid] = p0/sqrt(p0*p0+1);
    coord[4] /= d_beta0[tid];
  }
};

#undef X
#undef Y
#define X coord[0]
#define Y coord[2]
#define XP coord[1]
#define YP coord[3]
#define CT coord[4]
#define DP coord[5]

class gpu_track_through_csbendCSR_kernel3{
public:
  double n, e1, psi1, he1;
  // From CSRCSBEND struct
  int edge_order, edge1_effects;
  double b0, hgap, fint, h1;
  // From csbend.h globals
  double d_rho_actual, d_rho0;

  double length;
  gpu_track_through_csbendCSR_kernel3(double d_rho_actual, double d_rho0, 
     double n, double e1, double psi1, double he1, int edge_order, 
     int edge1_effects, double b0, double hgap, double fint, double h1) :
     d_rho_actual(d_rho_actual), d_rho0(d_rho0), n(n), e1(e1), psi1(psi1), he1(he1), 
     edge_order(edge_order), edge1_effects(edge1_effects), b0(b0),
     hgap(hgap), fint(fint), h1(h1) {};

  __device__ inline void operator()(gpuParticleAccessor& coord){

    double rho = (1+DP)*d_rho_actual;
    if (edge_order<=1 && edge1_effects==1) {
      double delta_xp = tan(e1)/rho*X;
      XP += delta_xp;
      YP -= tan(e1-psi1/(1+DP))/rho*Y;
    } else if (edge_order>=2 && edge1_effects==1) 
      gpu_apply_edge_effects(&X, &XP, &Y, &YP, rho, n, e1, he1, psi1*(1+DP), -1);
    else if (edge1_effects>=2) {
      //rho = (1+DP)*d_rho_actual;
      gpu_dipoleFringeSym(&X, &XP, &Y, &YP, &DP, d_rho_actual, -1., edge_order, b0/d_rho0, e1, 2*hgap, fint, h1);
    }
  }
};

class gpu_track_through_csbendCSR_kernel4{
public:
  double *d_x, *d_y;
  double e1;
  // From CSRCSBEND struct
  int d_expansionOrder;
  // From csbend.h globals
  double d_rho0, d_rad_coef;

  gpu_track_through_csbendCSR_kernel4(double *d_x, double *d_y, double e1,
      int d_expansionOrder, double d_rho0, double d_rad_coef) : d_x(d_x), d_y(d_y), e1(e1), 
      d_expansionOrder(d_expansionOrder), 
      d_rho0(d_rho0), d_rad_coef(d_rad_coef) {};

  __device__ inline void operator()(gpuParticleAccessor& coord){
    unsigned int tid = coord.getParticleIndex();

    double Fx, Fy;

    /* pre-adjust dp/p to anticipate error made by integrating over entire sector */
    gpu_computeCSBENDFields(&Fx, &Fy, d_x[tid], d_y[tid], d_expansionOrder);
  
    double dp_prime = -d_rad_coef*(sqr(Fx)+sqr(Fy))*sqr(1+DP)*
      sqrt(sqr(1+X/d_rho0)+sqr(XP)+sqr(YP));
    DP -= dp_prime*X*tan(e1);
  }
};

class gpu_track_through_csbendCSR_kernel5{
public:
  unsigned int *d_sortIndex;
  double *d_beta0, *d_gauss_rn;
  curandState_t *state;
  double Po, srGaussianLimit;
  // From CSRCSBEND struct
  int integration_order;
  double length, n_kicks;
  int d_expansionOrder;
  // From csbend.h globals
  double d_rho0, d_rho_actual, d_rad_coef, d_isrConstant;
  double d_meanPhotonsPerMeter0, d_normalizedCriticalEnergy0;
  int d_distributionBasedRadiation, d_includeOpeningAngle;

  gpu_track_through_csbendCSR_kernel5(unsigned int *d_sortIndex, double* d_beta0,
    double *d_gauss_rn, curandState_t *state, double Po, int integration_order,
    double length, double n_kicks, int d_expansionOrder, double d_rho0,
    double d_rho_actual, double d_rad_coef, double d_isrConstant,
    double d_meanPhotonsPerMeter0, double d_normalizedCriticalEnergy0,
    int d_distributionBasedRadiation, int d_includeOpeningAngle, 
    double srGaussianLimit) :
    d_sortIndex(d_sortIndex), d_beta0(d_beta0), d_gauss_rn(d_gauss_rn), 
    state(state), Po(Po), integration_order(integration_order), length(length), 
    n_kicks(n_kicks), d_expansionOrder(d_expansionOrder), d_rho0(d_rho0),
    d_rho_actual(d_rho_actual), d_rad_coef(d_rad_coef), 
    d_isrConstant(d_isrConstant),
    d_meanPhotonsPerMeter0(d_meanPhotonsPerMeter0),
    d_normalizedCriticalEnergy0(d_normalizedCriticalEnergy0),
    d_distributionBasedRadiation(d_distributionBasedRadiation),
    d_includeOpeningAngle(d_includeOpeningAngle),
    srGaussianLimit(srGaussianLimit) {};

  __device__ inline void operator()(gpuParticleAccessor& coord){
    unsigned int tid = coord.getParticleIndex();

    double p1, beta1;
    double Qi[6], Qf[6];

    /* load input coordinates into arrays */
    Qi[0] = X;
    Qi[1] = XP;
    Qi[2] = Y;
    Qi[3] = YP;
    Qi[4] = 0;  
    Qi[5] = DP;
    gpu_convertToDipoleCanonicalCoordinates(Qi, d_rho0, 0);
  
    int particle_lost = 0;
    double s_lost; /* unused here */
    double *t_gauss_rn = NULL;
    if (d_gauss_rn) t_gauss_rn = &d_gauss_rn[tid];
    curandState_t* t_state = NULL;
    if (state) t_state = &state[tid];
    if (integration_order==4)
      gpu_integrate_csbend_ord4(Qf, Qi, NULL, length/n_kicks,
         1, 0, d_rho0, Po, particle_lost, s_lost,
         d_rho_actual, d_rad_coef, d_isrConstant,
         d_distributionBasedRadiation, d_includeOpeningAngle,
         d_meanPhotonsPerMeter0, d_normalizedCriticalEnergy0,
         d_expansionOrder, t_gauss_rn, t_state, srGaussianLimit,
         NULL, 0);
    else
      gpu_integrate_csbend_ord2(Qf, Qi, NULL, length/n_kicks,
         1, 0, d_rho0, Po, particle_lost, s_lost,
         d_rho_actual, d_rad_coef, d_isrConstant,
         d_distributionBasedRadiation, d_includeOpeningAngle,
         d_meanPhotonsPerMeter0, d_normalizedCriticalEnergy0,
         d_expansionOrder, t_gauss_rn, t_state, srGaussianLimit,
         NULL, 0);

    /* Lost particles are killed in kernel7 and just marked here */
    if (particle_lost)
      d_sortIndex[tid] = tid + coord.getParticlePitch();
    else 
      d_sortIndex[tid] = tid;  
  
    /* retrieve coordinates from arrays */
    gpu_convertFromDipoleCanonicalCoordinates(Qf, d_rho0, 0);
    X  = Qf[0];  
    XP = Qf[1];  
    Y  = Qf[2];  
    YP = Qf[3];  
    DP = Qf[5];
  
    if (d_rad_coef || d_isrConstant) {
      /* convert additional distance traveled to ct using mean velocity */
      p1 = Po*(1+DP);
      beta1 = p1/sqrt(p1*p1+1);
      CT += Qf[4]*2/(d_beta0[tid]+beta1);
      d_beta0[tid] = beta1;
    } else {
      CT += Qf[4]/d_beta0[tid];  
    }
  }
};

class gpu_track_through_csbendCSR_kernel6{
public:
  unsigned int *d_sortIndex;
  double *d_dGamma;
  double Po, dct, ctLower, d_rho0;
  long nBins1; // nBins - 1

  gpu_track_through_csbendCSR_kernel6(unsigned int *d_sortIndex,
    double* d_dGamma, double Po, double dct, double ctLower, long nBins1,
    double d_rho0) :
    d_sortIndex(d_sortIndex), d_dGamma(d_dGamma), Po(Po), dct(dct), 
    ctLower(ctLower), nBins1(nBins1), d_rho0(d_rho0) {};

  __device__ inline void operator()(gpuParticleAccessor& coord){
    unsigned int tid = coord.getParticleIndex();

    long iBin;
    /* if particle is not lost */
    if (d_sortIndex[tid] < coord.getParticlePitch()) {
      double f;
      /* apply CSR kick */
      iBin = (f=(CT-ctLower)/dct);
      f -= iBin;
      if (iBin>=0 && iBin<nBins1)
	DP += ((1-f)*d_dGamma[iBin]+f*d_dGamma[iBin+1])/Po*(1+X/d_rho0);
    }
  }
};

class gpu_track_through_csbendCSR_kernel7{
public:
  unsigned int *d_sortIndex;
  double Po, e2, z_start;
  // From CSRCSBEND struct
  int d_expansionOrder;
  unsigned long edgeFlags;
  // From csbend.h globals
  double d_rho0, d_rad_coef;

  gpu_track_through_csbendCSR_kernel7(unsigned int *d_sortIndex, double Po,
    double e2, double z_start, int d_expansionOrder,
    unsigned long edgeFlags, double d_rho0, double d_rad_coef) : 
    d_sortIndex(d_sortIndex), Po(Po), e2(e2), z_start(z_start), 
    d_expansionOrder(d_expansionOrder),
    edgeFlags(edgeFlags), d_rho0(d_rho0), d_rad_coef(d_rad_coef) {};

  __device__ unsigned int operator()(gpuParticleAccessor& coord){
    unsigned int tid = coord.getParticleIndex();

    double dp_prime, p1, beta1;
    double Fx, Fy;

    if (edgeFlags&BEND_EDGE2_EFFECTS && e2!=0 && d_rad_coef) {
      /* post-adjust dp/p to correct error made by integrating over entire sector */
      gpu_computeCSBENDFields(&Fx, &Fy, X, Y, d_expansionOrder);
      
      dp_prime = -d_rad_coef*(sqr(Fx)+sqr(Fy))*sqr(1+DP)*
        sqrt(sqr(1+X/d_rho0)+sqr(XP)+sqr(YP));
      DP -= dp_prime*X*tan(e2);
    }

    /* convert CT to distance traveled at final velocity */
    p1 = Po*(1+DP);
    beta1 = p1/sqrt(sqr(p1)+1);
    coord[4] = CT*beta1;

    /* if particle is lost (can be lost in kernel5) */
    if (d_sortIndex[tid] > coord.getParticlePitch() || p1<=0) {
      // TODO s_lost should be stored from kernel5, this is
      // probably a bug on the CPU side as well.
      coord[4] = z_start;// + s_lost;
      coord[5] = Po*(1+coord[5]);
      d_sortIndex[tid] = tid + coord.getParticlePitch();
      return 0;
    }

    d_sortIndex[tid] = tid;
    return 1;
  }
};

class gpu_track_through_csbendCSR_kernel8{
public:
  double e2, psi2, n, he2;
  // From CSRCSBEND struct
  int edge_order, edge2_effects;
  double b0, hgap, fint, h2;
  // From csbend.h globals
  double d_rho_actual, d_rho0;

  gpu_track_through_csbendCSR_kernel8(double d_rho_actual, double d_rho0, 
     double n, double e2, double psi2, double he2, int edge_order, 
     int edge2_effects, double b0, double hgap, double fint, double h2) :
     d_rho_actual(d_rho_actual), d_rho0(d_rho0), n(n), e2(e2), psi2(psi2), he2(he2), 
     edge_order(edge_order), edge2_effects(edge2_effects), b0(b0),
     hgap(hgap), fint(fint), h2(h2) {};

  __device__ inline void operator()(gpuParticleAccessor& coord){

    double rho, delta_xp;

    /* apply edge focusing */
    rho = (1+DP)*d_rho_actual;
    if (edge_order<=1 && edge2_effects==1) {
      delta_xp = tan(e2)/rho*X;
      XP += delta_xp;
      YP -= tan(e2-psi2/(1+DP))/rho*Y;
    } else if (edge_order>=2 && edge2_effects==1) {
      gpu_apply_edge_effects(&X, &XP, &Y, &YP, rho, n, e2, he2, psi2*(1+DP), 1);
    } else if (edge2_effects>=2) {
      //rho = (1+DP)*rho_actual;
      gpu_dipoleFringeSym(&X, &XP, &Y, &YP, &DP, d_rho_actual, 1., edge_order, b0/d_rho0, e2, 2*hgap, fint, h2);
    }

  }
};

class gpu_track_through_csbendCSR_kernel9{
public:
  double dxf, dyf, dzf;
  double dcet0, dcet1, dcet2, dcet3;
  double cos_ttilt, sin_ttilt;

  gpu_track_through_csbendCSR_kernel9(double dxf, double dyf, double dzf,
      double dcet0, double dcet1, double dcet2, double dcet3,
      double cos_ttilt, double sin_ttilt) : dxf(dxf), dyf(dyf), dzf(dzf), 
      dcet0(dcet0), dcet1(dcet1), dcet2(dcet2), dcet3(dcet3),
      cos_ttilt(cos_ttilt), sin_ttilt(sin_ttilt) {};

  __device__ inline void operator()(gpuParticleAccessor& coord){

    double x, y, xp, yp;

    x  =  X*cos_ttilt -  Y*sin_ttilt + dcet0; // dcoord_etilt[0];
    y  =  X*sin_ttilt +  Y*cos_ttilt + dcet2; // dcoord_etilt[2];
    xp = XP*cos_ttilt - YP*sin_ttilt + dcet1; // dcoord_etilt[1];
    yp = XP*sin_ttilt + YP*cos_ttilt + dcet3; // dcoord_etilt[3];
    X  = x;
    Y  = y;
    XP = xp;
    YP = yp;
    coord[0] += dxf + dzf*coord[1];
    coord[2] += dyf + dzf*coord[3];
    coord[4] += dzf*sqrt(1+ sqr(coord[1]) + sqr(coord[3]));
  }
};


extern "C" {

static char *derbenevCriterionOption[N_DERBENEV_CRITERION_OPTIONS]
  = {(char*)"disable", (char*)"evaluate", (char*)"enforce"};

long gpu_track_through_csbendCSR(long n_part, CSRCSBEND *csbend, 
       double p_error, double Po, double **accepted, double z_start,
       double z_end, CHARGE *charge, char *rootname)
{
  double h, n, he1, he2;
  static long csrWarning = 0;
  static double *ctHist=NULL, *ctHistDeriv=NULL;
  static double *dGamma=NULL, *T1=NULL, *T2=NULL, *denom=NULL, *chik=NULL, *grnk=NULL;
  static long maxBins = 0 ;
  double ctLower(0), ctUpper(0), dct(0), slippageLength, phiBend, slippageLength13;
  long diSlippage, diSlippage4;
  long nBins, nBinned = 0;
  long kick;
  double fse;
  double tilt, etilt, cos_ttilt, sin_ttilt, ttilt;
  double angle, e1, e2, Kg;
  double psi1, psi2;
  double dcoord_etilt[6];
  double dxi, dyi, dzi;
  double dxf, dyf, dzf;
  double macroParticleCharge, CSRConstant, gamma2, gamma3;
  long iBin, iBinBehind;
  long csrInhibit = 0, largeRhoWarning = 0;
  double derbenevRatio = 0;
  long n_partMoreThanOne = 0;
  TRACKING_CONTEXT tContext;
  VMATRIX *Msection=NULL, *Me1=NULL, *Me2=NULL;
  static double accumulatedAngle = 0;
  short accumulatingAngle = 1;
#if USE_MPI
  double *buffer;
#endif
#ifdef DEBUG_IGF
  FILE *fpdeb;
  fpdeb = fopen("csr.sdds","w");
  fprintf(fpdeb, "SDDS1\n&parameter name = Slice, type=long &end\n");
  fprintf(fpdeb, "&column name=s, type=double, units=m &end\n");
  fprintf(fpdeb, "&column name=iBin, type=long &end\n");
  fprintf(fpdeb, "&column name=Chi, type=double &end\n");
  fprintf(fpdeb, "&column name=G, units=V/m, type=double &end\n");
  fprintf(fpdeb, "&column name=dGamma, type=double &end\n");
  fprintf(fpdeb, "&data mode=ascii &end\n");
#endif

  struct GPUBASE* gpuBase = getGpuBase();
  unsigned int particlePitch = gpuBase->gpu_array_pitch;
  unsigned int* d_sortIndex = gpuBase->d_tempu_alpha;
  curandState_t* state = NULL;
  double* d_gauss_rn = gpuBase->d_temp_particles;
  double* d_beta0    = gpuBase->d_temp_particles + 1*particlePitch;
  double* d_dGamma   = gpuBase->d_temp_particles + 2*particlePitch;
  double* d_x        = gpuBase->d_temp_particles + 3*particlePitch;
  double* d_y        = gpuBase->d_temp_particles + 4*particlePitch;
  double* d_ctHist;
  bool allocGpuMem=false;
  /* always alloc d_ctHist to avoid conflict with killParticles */
  cudaMalloc( (void**)&d_ctHist, sizeof(double)*csbend->bins);
  if ((isSlave || !notSinglePart) && particlePitch < csbend->bins) {
    cudaMalloc( (void**)&d_dGamma, sizeof(double)*csbend->bins);
    allocGpuMem=true;
  }

  gamma2 = Po*Po+1;
  gamma3 = pow(gamma2, 3./2);

#if USE_MPI
  if (notSinglePart)
    n_partMoreThanOne = 1; /* This is necessary to solve synchronization issue in parallel version*/
  else
    if (n_part > 1) n_partMoreThanOne = 1;	
#else
  if (n_part > 1) n_partMoreThanOne = 1;
#endif

  if (!(csbend->edgeFlags&SAME_BEND_PRECEDES))
    accumulatedAngle = accumulatingAngle = 0;
  

  csrWake.valid = 0;
  refTrajectoryMode = 0;
  if (isSlave || !notSinglePart) 
    reset_driftCSR();

  getTrackingContext(&tContext);
  
  if (!csbend)
    bombElegant("null CSRCSBEND pointer (track_through_csbend)", NULL);
  if (csbend->integratedGreensFunction && !csbend->steadyState) 
    bombElegant("CSRCSBEND requires STEADYSTATE=1 if IGF=1.", NULL);
  if (csbend->edge_order>1 && (csbend->edge1_effects==2 || csbend->edge2_effects==2) && csbend->hgap==0)
    bombElegant("CSRCSBEND has EDGE_ORDER>1 and EDGE[12]_EFFECTS==2, but HGAP=0. This gives undefined results.", NULL);

  if (csbend->angle==0) {
    if (!csbend->useMatrix)
      gpu_exactDrift(n_part, csbend->length); 
    else {
      if (isSlave || !notSinglePart) {
        gpuDriver(n_part, gpu_track_through_csbendCSR_kernel1(csbend->length));
        gpuErrorHandler("gpu_track_through_csbendCSR: gpu_track_through_csbendCSR_kernel1");
      }
    }
    return n_part;
  }

  if (csbend->integration_order!=2 && csbend->integration_order!=4)
    bombElegant("CSBEND integration_order is invalid--must be either 2 or 4", NULL);

  macroParticleCharge = 0;
  if (charge) {
    macroParticleCharge = charge->macroParticleCharge;
  } else if (csbend->bins && !csrWarning && csbend->csr) {
    fprintf(stdout, "Warning: you asked for CSR on CSBEND but didn't give a CHARGE element\n");
    fflush(stdout);
    csrWarning = 1;
  }
  
  if ((nBins=csbend->bins)<2)
    bombElegant("Less than 2 bins for CSR!", NULL);

  if (csbend->SGDerivHalfWidth<=0)
    csbend->SGDerivHalfWidth = csbend->SGHalfWidth;
  if (csbend->SGDerivHalfWidth<=0)
    csbend->SGDerivHalfWidth = 1;

  if (csbend->SGDerivOrder<=0)
    csbend->SGDerivOrder = csbend->SGOrder;
  if (csbend->SGDerivOrder<=0)
    csbend->SGDerivOrder = 1;
  
  rho0 = csbend->length/csbend->angle;
  if (csbend->use_bn) {
    csbend->b[0] = csbend->b1;
    csbend->b[1] = csbend->b2;
    csbend->b[2] = csbend->b3;
    csbend->b[3] = csbend->b4;
    csbend->b[4] = csbend->b5;
    csbend->b[5] = csbend->b6;
    csbend->b[6] = csbend->b7;
    csbend->b[7] = csbend->b8;
  } else {
    csbend->b[0] = csbend->k1*rho0;
    csbend->b[1] = csbend->k2*rho0;
    csbend->b[2] = csbend->k3*rho0;
    csbend->b[3] = csbend->k4*rho0;
    csbend->b[4] = csbend->k5*rho0;
    csbend->b[5] = csbend->k6*rho0;
    csbend->b[6] = csbend->k7*rho0;
    csbend->b[7] = csbend->k8*rho0;
  }

  he1 = csbend->h1;
  he2 = csbend->h2;
  if (csbend->angle<0) {
    long i;
    angle = -csbend->angle;
    e1    = -csbend->e1;
    e2    = -csbend->e2;
    etilt = csbend->etilt;
    tilt  = csbend->tilt + PI;
    rho0  = csbend->length/angle;
    for (i=0; i<8; i+=2)
      csbend->b[i] *= -1;
  }
  else {
    angle = csbend->angle;
    e1    = csbend->e1;
    e2    = csbend->e2;
    etilt = csbend->etilt;
    tilt  = csbend->tilt;
    rho0  = csbend->length/angle;
  }
  
  if (rho0>1e6) {
    if (!largeRhoWarning) {
      printf("Warning: One or more CSRCSBENDs have radius > 1e6.  Treated as drift.\n");
      largeRhoWarning = 1;
    }
    gpu_exactDrift(n_part, csbend->length);
    return n_part;
  }
  
  h = 1/rho0;
  n = -csbend->b[0]/h;
  fse = csbend->fse;
  if (fse>-1)
    rho_actual = 1/((1+fse)*h);
  else
    rho_actual = 1e16/h;

  /* angles for fringe-field effects */
  Kg   = 2*csbend->hgap*csbend->fint;
  psi1 = Kg/rho_actual/cos(e1)*(1+sqr(sin(e1)));
  psi2 = Kg/rho_actual/cos(e2)*(1+sqr(sin(e2)));

  /* rad_coef is d((P-Po)/Po)/ds for the on-axis, on-momentum particle, where po is the momentum of
   * the central particle.
   */
  if (csbend->synch_rad)
    rad_coef = sqr(particleCharge)*pow3(Po)*sqr(1+fse)/(6*PI*epsilon_o*sqr(c_mks)*particleMass*sqr(rho0));
  else
    rad_coef = 0;
  /* isrConstant is the RMS increase in dP/P per meter due to incoherent SR.  */
  if (csbend->isr && (n_part>1 || !csbend->isr1Particle)) 
    isrConstant = particleRadius*sqrt(55.0/(24*sqrt(3))*pow5(Po)*
                              137.0359895/pow3(fabs(rho_actual)));
  else
    isrConstant = 0;

  distributionBasedRadiation = 0;
  
  if (csbend->useMatrix) {
    csbend->nonlinear = 0;
    Me1 = edge_matrix(e1, 1./(rho0/(1+csbend->fse)), 0.0, n, -1, Kg, 1, 0, 0);
    Msection = bend_matrix(csbend->length/csbend->n_kicks, 
                                   angle/csbend->n_kicks, 0.0, 0.0, 
                                   0.0, 0.0, csbend->b[0]*h,  0.0,
                                   0.0, 0.0, 0.0, csbend->fse, csbend->etilt, 1, 1, 0, 0);
    Me2 = edge_matrix(e2, 1./(rho0/(1+csbend->fse)), 0.0, n, 1, Kg, 1, 0, 0);
  }
  computeCSBENDFieldCoefficients(csbend->b, h, csbend->nonlinear, csbend->expansionOrder);

  ttilt = tilt + etilt;
  if (ttilt==0) {
    cos_ttilt = 1;
    sin_ttilt = 0;
  }
  else if (fabs(fabs(ttilt)-PI)<1e-12) {
    cos_ttilt = -1;
    sin_ttilt = 0;
  }
  else if (fabs(ttilt-PIo2)<1e-12) {
    cos_ttilt = 0;
    sin_ttilt = 1;
  }
  else if (fabs(ttilt+PIo2)<1e-12) {
    cos_ttilt = 0;
    sin_ttilt = -1;
  }
  else {
    cos_ttilt = cos(ttilt);
    sin_ttilt = sin(ttilt);
  }


  if (etilt) {
    /* compute final offsets due to error-tilt of the magnet */
    /* see pages 90-93 of notebook 1 about this */
    double q1a, q2a, q3a;
    double q1b, q2b, q3b;
    double qp1, qp2, qp3; 
    double dz, tan_alpha, k;

    q1a = (1-cos(angle))*rho0*(cos(etilt)-1);
    q2a = 0;
    q3a = (1-cos(angle))*rho0*sin(etilt);
    qp1 = sin(angle)*cos(etilt);
    qp2 = cos(angle);
    k = sqrt(sqr(qp1)+sqr(qp2));
    qp1 /= k;
    qp2 /= k;
    qp3 = sin(angle)*sin(etilt)/k;
    tan_alpha = 1./tan(angle)/cos(etilt);
    q1b = q1a*tan_alpha/(tan(angle)+tan_alpha);
    q2b = -q1b*tan(angle);
    dz  = sqrt(sqr(q1b-q1a)+sqr(q2b-q2a));
    q3b = q3a + qp3*dz;

    dcoord_etilt[0] = sqrt(sqr(q1b) + sqr(q2b));
    dcoord_etilt[1] = tan(atan(tan_alpha)-(PIo2-angle));
    dcoord_etilt[2] = q3b;
    dcoord_etilt[3] = qp3;
    dcoord_etilt[4] = dz*sqrt(1+sqr(qp3));
    dcoord_etilt[5] = 0;

    /* rotate by tilt to get into same frame as bend equations. */
    rotate_coordinates(dcoord_etilt, tilt);
  }
  else
    fill_double_array(dcoord_etilt, 6L, 0.0);

  dxi = -csbend->dx;
  dzi =  csbend->dz;
  dyi = -csbend->dy;

  /* must use the original angle here because the translation is done after
   * the final rotation back
   */
  dxf =  csbend->dx*cos(csbend->angle) + csbend->dz*sin(csbend->angle);
  dzf =  csbend->dx*sin(csbend->angle) - csbend->dz*cos(csbend->angle);
  dyf = csbend->dy;

  if (isMaster) {
  if (csbend->particleOutputFile && strlen(csbend->particleOutputFile) && !csbend->particleFileActive) {
    /* set up SDDS output file for particle coordinates inside bend */
    csbend->particleFileActive = 1;
    csbend->particleOutputFile = compose_filename(csbend->particleOutputFile, rootname);
    if (!SDDS_InitializeOutput(&csbend->SDDSpart, SDDS_BINARY, 1, 
                               NULL, NULL, csbend->particleOutputFile) ||
        !SDDS_DefineSimpleParameter(&csbend->SDDSpart, "Pass", NULL, SDDS_LONG) ||
        !SDDS_DefineSimpleParameter(&csbend->SDDSpart, "Kick", NULL, SDDS_LONG) ||
        !SDDS_DefineSimpleParameter(&csbend->SDDSpart, "pCentral", "m$be$nc", SDDS_DOUBLE) ||
        !SDDS_DefineSimpleParameter(&csbend->SDDSpart, "Angle", NULL, SDDS_DOUBLE) ||
        (csbend->xIndex=SDDS_DefineColumn(&csbend->SDDSpart, "x", NULL, "m", 
                                          NULL, NULL, SDDS_DOUBLE, 0 ))<0 ||
        (csbend->xpIndex=SDDS_DefineColumn(&csbend->SDDSpart, "xp", NULL, NULL, 
                                           NULL, NULL, SDDS_DOUBLE, 0))<0 ||
        (csbend->tIndex=SDDS_DefineColumn(&csbend->SDDSpart, "t", NULL, "s", 
                                          NULL, NULL, SDDS_DOUBLE, 0))<0 ||
        (csbend->pIndex=SDDS_DefineColumn(&csbend->SDDSpart, "p", NULL, "m$be$nc", 
                                          NULL, NULL, SDDS_DOUBLE, 0))<0 ||
        !SDDS_WriteLayout(&csbend->SDDSpart)) {
      SDDS_SetError((char*)"Problem setting up particle output file for CSR");
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
    }
  }
  }
  
  if (isMaster) { 
  if (csbend->histogramFile && strlen(csbend->histogramFile) && !csbend->wakeFileActive) {
    /* set up SDDS output file for CSR monitoring */
    csbend->wakeFileActive = 1;
    csbend->histogramFile = compose_filename(csbend->histogramFile, rootname);
    if (!SDDS_InitializeOutput(&csbend->SDDSout, SDDS_BINARY, 1, NULL, NULL, csbend->histogramFile) ||
        !SDDS_DefineSimpleParameter(&csbend->SDDSout, "Pass", NULL, SDDS_LONG) ||
        !SDDS_DefineSimpleParameter(&csbend->SDDSout, "Kick", NULL, SDDS_LONG) ||
        !SDDS_DefineSimpleParameter(&csbend->SDDSout, "pCentral", "m$be$nc", SDDS_DOUBLE) ||
        !SDDS_DefineSimpleParameter(&csbend->SDDSout, "Angle", NULL, SDDS_DOUBLE) ||
        !SDDS_DefineSimpleParameter(&csbend->SDDSout, "SlippageLength", "m", SDDS_DOUBLE) ||
        !SDDS_DefineSimpleParameter(&csbend->SDDSout, "TotalBunchLength", "m", SDDS_DOUBLE) ||
        !SDDS_DefineSimpleParameter(&csbend->SDDSout, "BinSize", "m", SDDS_DOUBLE) ||
        !SDDS_DefineSimpleParameter(&csbend->SDDSout, "dsKick", "m", SDDS_DOUBLE) ||
        !SDDS_DefineSimpleParameter(&csbend->SDDSout, "DerbenevRatio", NULL, SDDS_DOUBLE) ||
        !SDDS_DefineSimpleColumn(&csbend->SDDSout, "s", "m", SDDS_DOUBLE) ||
        !SDDS_DefineSimpleColumn(&csbend->SDDSout, "LinearDensity", "C/s", SDDS_DOUBLE) ||
        !SDDS_DefineSimpleColumn(&csbend->SDDSout, "LinearDensityDeriv", "C/s$a2$n", SDDS_DOUBLE) ||
        !SDDS_DefineSimpleColumn(&csbend->SDDSout, "DeltaGamma", NULL, SDDS_DOUBLE) ||
        !SDDS_DefineSimpleColumn(&csbend->SDDSout, "GammaDeriv", "1/m", SDDS_DOUBLE) ||
        !SDDS_DefineSimpleColumn(&csbend->SDDSout, "DeltaGammaT1", NULL, SDDS_DOUBLE) ||
        !SDDS_DefineSimpleColumn(&csbend->SDDSout, "DeltaGammaT2", NULL, SDDS_DOUBLE) ||
        !SDDS_WriteLayout(&csbend->SDDSout)) {
      SDDS_SetError((char*)"Problem setting up wake output file for CSR");
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
    }
  }
  }
  if (csbend->wakeFilterFile && strlen(csbend->wakeFilterFile) && !csbend->wffValues) 
    readWakeFilterFile(&csbend->wffValues,
                       &csbend->wffFreqValue, &csbend->wffRealFactor, &csbend->wffImagFactor, 
                       csbend->wffFreqColumn, csbend->wffRealColumn, csbend->wffImagColumn,
                       csbend->wakeFilterFile);
  
  /*  prepare arrays for CSR integrals */
  nBins = csbend->bins;
  if (!(ctHist=(double*)SDDS_Realloc(ctHist, sizeof(*ctHist)*nBins)) ||
      !(ctHistDeriv=(double*)SDDS_Realloc(ctHistDeriv, sizeof(*ctHistDeriv)*nBins)) ||
      !(denom=(double*)SDDS_Realloc(denom, sizeof(*denom)*nBins)) ||
      !(T1=(double*)SDDS_Realloc(T1, sizeof(*T1)*nBins)) ||
      !(T2=(double*)SDDS_Realloc(T2, sizeof(*T2)*nBins)) ||
      !(dGamma=(double*)SDDS_Realloc(dGamma, sizeof(*dGamma)*nBins)))
    bombElegant("memory allocation failure (track_through_csbendCSR)", NULL);

  /* prepare some data for CSRDRIFT */
  csrWake.dGamma = dGamma;
  csrWake.bins = nBins;
  csrWake.ds0 = csbend->length/csbend->n_kicks;
  csrWake.zLast = csrWake.z0 = z_end;
  csrWake.highFrequencyCutoff0 = csbend->highFrequencyCutoff0;
  csrWake.highFrequencyCutoff1 = csbend->highFrequencyCutoff1;
  csrWake.lowFrequencyCutoff0 = csbend->lowFrequencyCutoff0;
  csrWake.lowFrequencyCutoff1 = csbend->lowFrequencyCutoff1;
  csrWake.clipNegativeBins = csbend->clipNegativeBins;
  csrWake.wffValues = csbend->wffValues;
  csrWake.wffFreqValue = csbend->wffFreqValue;
  csrWake.wffRealFactor = csbend->wffRealFactor;
  csrWake.wffImagFactor = csbend->wffImagFactor;
  
#if !defined(PARALLEL)  
  multipoleKicksDone += n_part*csbend->n_kicks*(csbend->integration_order==4?4:1);
#endif

  if (isSlave || !notSinglePart) {
    gpuDriver(n_part, gpu_track_through_csbendCSR_kernel2(d_beta0, d_x, d_y, Po, 
              dxi, dyi, dzi, cos_ttilt, sin_ttilt));
    gpuErrorHandler("gpu_track_through_csbendCSR: gpu_track_through_csbendCSR_kernel2");
    if (csbend->edgeFlags&BEND_EDGE1_EFFECTS) {
      /* apply edge focusing */
      if (csbend->useMatrix)
        gpu_track_particles(Me1, n_part);
      else {
        gpuDriver(n_part, gpu_track_through_csbendCSR_kernel3(rho_actual, rho0, n, e1,
                  psi1, he1, csbend->edge_order, csbend->edge1_effects, 
                  csbend->b[0], csbend->hgap, csbend->fint, csbend->h1));
        gpuErrorHandler("gpu_track_through_csbendCSR: gpu_track_through_csbendCSR_kernel3");
      }
    }
    if (csbend->edgeFlags&BEND_EDGE1_EFFECTS && e1!=0 && rad_coef) {
      gpuDriver(n_part, gpu_track_through_csbendCSR_kernel4(d_x, d_y, e1,
                expansionOrder1, rho0, rad_coef));
      gpuErrorHandler("gpu_track_through_csbendCSR: gpu_track_through_csbendCSR_kernel4");
    }
  }
  if (csbend->csr && n_partMoreThanOne)
    CSRConstant = 2*macroParticleCharge*particleCharge/pow(3*rho0*rho0, 1./3.)/(4*PI*epsilon_o*particleMass*sqr(c_mks));
  else
    CSRConstant = 0;
  /* Now do the body of the sector dipole */
  phiBend = accumulatedAngle;
  if (isSlave || !notSinglePart) {
    /* Copy Fx_xy, Fy_xy and table coefficients */
    double Fall_xy[2*xOrderMax2];
    for (int ii=0; ii<expansionOrderMax; ii++) {
      for (int jj=0; jj<expansionOrderMax; jj++) {
        Fall_xy[ii*expansionOrderMax+jj]            = Fx_xy[ii][jj];
        Fall_xy[ii*expansionOrderMax+jj+xOrderMax2] = Fy_xy[ii][jj];
      }
    }
    cudaMemcpyToSymbol(c_Fx_xy, Fall_xy, sizeof(double)*xOrderMax2,
                       0, cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol(c_Fy_xy, &(Fall_xy[xOrderMax2]), sizeof(double)*xOrderMax2,
                       0, cudaMemcpyHostToDevice);
    gpu_initializeInterpolationTables();

    /* initalize random numbers as needed. */
    if (rad_coef || isrConstant)
      if (!distributionBasedRadiation && isrConstant)
        gpu_d_gauss_rn_lim(d_gauss_rn, n_part, 0.0, 1.0, srGaussianLimit, random_2(0));
      else if (distributionBasedRadiation)
        state =
          (curandState_t*)gpu_get_rand_state(d_gauss_rn, n_part, random_2(0));
  }
  for (kick=0; kick<csbend->n_kicks; kick++) {
    if (isSlave || !notSinglePart) {
      if (csbend->useMatrix) {
        gpu_track_particles(Msection, n_part);
      } else {
        gpuDriver(n_part,
          gpu_track_through_csbendCSR_kernel5(d_sortIndex, d_beta0, 
          d_gauss_rn, state, Po, csbend->integration_order, csbend->length, 
          csbend->n_kicks, expansionOrder1, rho0, rho_actual, rad_coef,
          isrConstant, meanPhotonsPerMeter0, normalizedCriticalEnergy0,
          distributionBasedRadiation, includeOpeningAngle, srGaussianLimit));
        gpuErrorHandler("gpu_track_through_csbendCSR: gpu_track_through_csbendCSR_kernel5");
      }
    }

    if (n_partMoreThanOne && csbend->derbenevCriterionMode) {
      /* evaluate Derbenev criterion from TESLA-FEL 1995-05: sigma_x/sigma_z << (R/sigma_z)^(1/3) */
      long code;
      double Sz, Sx;
      switch (code=match_string(csbend->derbenevCriterionMode, derbenevCriterionOption, N_DERBENEV_CRITERION_OPTIONS, 0)) {
      case DERBENEV_CRITERION_DISABLE:
	break;
      case DERBENEV_CRITERION_EVAL:
      case DERBENEV_CRITERION_ENFORCE:
#if !USE_MPI
	gpu_rms_emittance(4, 5, n_part, &Sz, NULL, NULL);
	gpu_rms_emittance(0, 1, n_part, &Sx, NULL, NULL);
#else
     if (notSinglePart) {
        /* The master will get the result from the rms_emittance routine */
	gpu_rms_emittance_p(4, 5, n_part, &Sz, NULL, NULL);
	gpu_rms_emittance_p(0, 1, n_part, &Sx, NULL, NULL);
     } else {
        gpu_rms_emittance(4, 5, n_part, &Sz, NULL, NULL);
        gpu_rms_emittance(0, 1, n_part, &Sx, NULL, NULL);
     }
#endif
	Sz = sqrt(Sz);
	Sx = sqrt(Sx);
	derbenevRatio = (Sx/Sz)/pow(rho0/Sz, 1./3.);
	if (derbenevRatio>0.1) {
	  if (code==DERBENEV_CRITERION_EVAL)
	    fprintf(stderr, "Warning: Using 1-D CSR formalism but Derbenev criterion not satisfied (%le > 0.1).\n",
		    derbenevRatio);
	  else {
	    csrInhibit = 1;
	    fprintf(stderr, "Warning: Derbenev criterion not satisfied (%le > 0.1)---not applying CSR\n",
		    derbenevRatio);
	  }
	}
	break;
      default:
	fprintf(stderr, "Error: invalid value for DERBENEV_CRITERION_MODE. Give 'disable', 'evaluate', or 'enforce'\n");
	exit(1);
	break;
      }
    }
    

#if (!USE_MPI)
    if (n_partMoreThanOne && !csrInhibit) {
#else
      if (!csrInhibit && (notSinglePart || (!notSinglePart && n_partMoreThanOne))) { /* n_part could be 0 for some processors, which could cause synchronization problem */
#endif
      /* compute CSR potential function */
      if (kick==0 || !csbend->binOnce) {
        /* - first make a density histogram */
        ctLower = ctUpper = dct = 0;
        nBinned = gpu_binParticleCoordinate(ctHist, d_ctHist, &maxBins, &ctLower,
                    &ctUpper, &dct, &nBins,
                    csbend->binRangeFactor<1.1?1.1:csbend->binRangeFactor,
                    n_part, 4);
 
#if (!USE_MPI) 
	if (nBinned != n_part) {
          fprintf(stdout, "Only %ld of %ld particles binned for CSRCSBEND (z0=%le, kick=%ld, BRF=%le)\n", 
		  nBinned, n_part, z_start, kick, csbend->binRangeFactor<1.1?1.1:csbend->binRangeFactor);
	  fprintf(stdout, "ct min, max = %21.15e, %21.15e, dct = %21.15e, nBins=%ld, maxBins=%ld\n",
		  ctLower, ctUpper, dct, nBins, maxBins);
          fflush(stdout);
        }
#else
     if (notSinglePart) {
	if (USE_MPI) {
	  long all_binned, result = 1, nBinned_total;

          if (isSlave || !notSinglePart) {
	    result = ((nBinned==n_part) ? 1 : 0);
	  }
	  MPI_Allreduce(&result, &all_binned, 1, MPI_LONG, MPI_LAND, MPI_COMM_WORLD);
	  MPI_Allreduce(&nBinned, &nBinned_total, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
          nBinned = nBinned_total; 
	  if (!all_binned && isMaster) {
	    fprintf(stdout, "Not all particles binned for CSRCSBEND (z0=%le, kick=%ld, BRF=%le)\n", 
		    z_start, kick,
		  csbend->binRangeFactor<1.1?1.1:csbend->binRangeFactor);
	  fprintf(stdout, "ct min, max = %21.15e, %21.15e, dct = %21.15e, nBins=%ld, maxBins=%ld\n",
		  ctLower, ctUpper, dct, nBins, maxBins);
          fflush(stdout);
        }
        }

	if (USE_MPI) {  /* Master needs to know the information to write the result */
	  buffer = (double*)malloc(sizeof(double) * nBins);
	  MPI_Allreduce(ctHist, buffer, nBins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	  memcpy(ctHist, buffer, sizeof(double)*nBins);
	  free(buffer);
	}
     }
#endif
        
        /* - smooth the histogram, normalize to get linear density, and 
           copy in preparation for taking derivative
           */
        if (csbend->highFrequencyCutoff0>0 || csbend->lowFrequencyCutoff0>=0) {
          long nz;
          nz = applyLHPassFilters(ctHist, nBins, 
                                  csbend->lowFrequencyCutoff0, csbend->lowFrequencyCutoff1,
                                  csbend->highFrequencyCutoff0, csbend->highFrequencyCutoff1,
                                  csbend->clipNegativeBins);
          if (nz && negativeWarningsLeft) {
	    fprintf(stdout, "Warning: low pass filter resulted in negative values in %ld bins\n", nz);
            if (--negativeWarningsLeft==0)
              fprintf(stdout, "         Further warnings will be suppressed for this run.\n");
            fflush(stdout);
          }
        }
        if (csbend->SGHalfWidth>0) {
          SavitzyGolaySmooth(ctHist, nBins, csbend->SGOrder, csbend->SGHalfWidth, csbend->SGHalfWidth,  0);
          correctDistribution(ctHist, nBins, 1.0*nBinned);
        }
        for (iBin=0; iBin<nBins; iBin++) {
          denom[iBin] = pow(dct*iBin, 1./3.);
          ctHistDeriv[iBin] = (ctHist[iBin] /= dct);
        }
        /* - compute derivative with smoothing.  The deriv is w.r.t. index number and
         * I won't scale it now as it will just fall out in the integral 
         */
        SavitzyGolaySmooth(ctHistDeriv, nBins, csbend->SGDerivOrder, 
                           csbend->SGDerivHalfWidth, csbend->SGDerivHalfWidth, 1);
      } else {
        ctLower += rho0*angle/csbend->n_kicks;
        ctUpper += rho0*angle/csbend->n_kicks;
      }
      
      
      phiBend += angle/csbend->n_kicks;
      slippageLength = rho0*pow3(phiBend)/24.0;
      slippageLength13 = pow(slippageLength, 1./3.);
      diSlippage = slippageLength/dct;
      diSlippage4 = 4*slippageLength/dct;
      if (kick==0 || !csbend->binOnce) {
        if (csbend->integratedGreensFunction) {
          /* Integrated Greens function method */
          double const2;
          double z, xmu, a, b, frac, const1;
          if (kick==0) {
            if (!csbend->steadyState)
              bombElegant("Must have STEADY_STATE=1 when IGF=1\n", NULL);
            if (!(grnk=(double*)SDDS_Realloc(grnk, sizeof(*grnk)*nBins)) ||
                !(chik=(double*)SDDS_Realloc(chik, sizeof(*chik)*nBins)))
              bombElegant("memory allocation failure (track_through_csbendCSR)", NULL);
          }
          frac = 9.0/16.0;
          const1 = 6.0-log(27.0/4.0);
          for (iBin=0; iBin<nBins; iBin++) {
            z   = iBin*dct;
            xmu = 3.0*gamma3*z/(2.0*rho0);
            a   = sqrt(xmu*xmu+1.0);
            b   = a+xmu;
            if (xmu < 1e-3)
              chik[iBin] = frac*const1 + 0.50*sqr(xmu)-(7.0/54.0)*pow4(4)+(140.0/2187.0)*pow6(xmu);
            else
              chik[iBin] = frac*( 3.0*( -2.0*xmu*pow(b,1.0/3.0) + pow(b,2.0/3.0) + pow(b,4.0/3.0) ) +
                                 log( pow((1-pow(b,2.0/3.0))/xmu,2)  / (1+pow(b,2.0/3.0)+pow(b,4.0/3.0)) ) );
          }
          const2 = (16.0/27.0)*(particleCharge/(4*PI*epsilon_o))/(gamma2*dct);
          grnk[0] = const2*(chik[1]-chik[0]);
          for (iBin=1; iBin<nBins-1; iBin++)
            grnk[iBin] = const2*(chik[iBin+1] - 2.0*chik[iBin] + chik[iBin-1] );
          grnk[nBins-1] = 0;
        } else {
	  for (iBin=0; iBin<nBins; iBin++) {
	    double term1, term2;
	    long count;
	    T1[iBin] = T2[iBin] = 0;
	    term1 = term2 = 0;
	    if (CSRConstant) {
	      if (csbend->steadyState) {
                if (!csbend->integratedGreensFunction) {
	          if (!csbend->trapazoidIntegration) {
	  	  for (iBinBehind=iBin+1; iBinBehind<nBins; iBinBehind++)
	  	    T1[iBin] += ctHistDeriv[iBinBehind]/denom[iBinBehind-iBin];
	          }
	          else {
	  	  if ((iBinBehind=iBin+1)<nBins)
	  	    term1 = ctHistDeriv[iBinBehind]/denom[iBinBehind-iBin];
	  	  for (count=0, iBinBehind=iBin+1; iBinBehind<nBins; iBinBehind++, count++)
	  	    T1[iBin] += (term2=ctHistDeriv[iBinBehind]/denom[iBinBehind-iBin]);
	  	  if ((iBin+1)<nBins)
	  	    T1[iBin] += 0.3*sqr(denom[1])*(2*ctHistDeriv[iBin+1]+3*ctHistDeriv[iBin])/dct;
	  	  if (count>1)
	  	    T1[iBin] -= (term1+term2)/2;
	          }
                }
	      } else {
                /* Transient CSR */
	        if (!csbend->trapazoidIntegration) {
	  	for (iBinBehind=iBin+1; iBinBehind<=(iBin+diSlippage) && iBinBehind<nBins; iBinBehind++)
	  	  T1[iBin] += ctHistDeriv[iBinBehind]/denom[iBinBehind-iBin];
	        }
	        else {
	  	if ((iBinBehind = iBin+1)<nBins && iBinBehind<=(iBin+diSlippage))
	  	  term1 = ctHistDeriv[iBinBehind]/denom[iBinBehind-iBin]/2;
	  	for (count=0, iBinBehind = iBin+1; iBinBehind<=(iBin+diSlippage) && iBinBehind<nBins; 
	  	     count++, iBinBehind++)
	  	  T1[iBin] += (term2=ctHistDeriv[iBinBehind]/denom[iBinBehind-iBin]);
	  	if (diSlippage>0 && (iBin+1)<nBins)
	  	  T1[iBin] += 0.3*sqr(denom[1])*(2*ctHistDeriv[iBin+1]+3*ctHistDeriv[iBin])/dct;
	  	if (count>1)
	  	  T1[iBin] -= (term1+term2)/2;
	        }
	        if ((iBin+diSlippage)<nBins)
	  	T2[iBin] += ctHist[iBin+diSlippage];
	        if ((iBin+diSlippage4)<nBins)
	  	T2[iBin] -= ctHist[iBin+diSlippage4];
	      }
	      /* there is no negative sign here because my derivative is w.r.t. -s
	         in notation of Saldin, et. al. */
	      T1[iBin] *= CSRConstant*csbend->length/csbend->n_kicks; 
	      /* keep the negative sign on this term, which has no derivative */
	      T2[iBin] *= -CSRConstant*csbend->length/csbend->n_kicks/slippageLength13;
	    }
	    dGamma[iBin] = T1[iBin]+T2[iBin];
	  }
        }

        if (csbend->integratedGreensFunction) {
          convolveArrays1(dGamma, nBins, ctHist, grnk);
          for (iBin=0; iBin<nBins; iBin++)
            dGamma[iBin] *= -macroParticleCharge/(particleMass*sqr(c_mks))*csbend->length/csbend->n_kicks;
#ifdef DEBUG_IGF
          fprintf(fpdeb, "%ld\n%ld\n", kick, nBins);
          for (iBin=0; iBin<nBins; iBin++)
            fprintf(fpdeb, "%le %ld %le %le %le\n", iBin*dct, iBin, chik[iBin], grnk[iBin], dGamma[iBin]);
#endif
        }


	if (csbend->wffValues) 
	  applyFilterTable(dGamma, nBins, dct/c_mks, csbend->wffValues, csbend->wffFreqValue,
			   csbend->wffRealFactor, csbend->wffImagFactor);
      }
      if (isSlave || !notSinglePart) {
	if (CSRConstant) {
          cudaMemcpy((void*)d_dGamma, (void *)dGamma, sizeof(double)*nBins,
                     cudaMemcpyHostToDevice);
          gpuDriver(n_part,
            gpu_track_through_csbendCSR_kernel6(d_sortIndex, d_dGamma,
            Po, dct, ctLower, nBins-1, rho0));
          gpuErrorHandler("gpu_track_through_csbendCSR: gpu_track_through_csbendCSR_kernel6");
	}
      }
  
      if (csbend->particleFileActive && kick%csbend->particleOutputInterval==0) {
	if (isMaster) {
          // TODO
          printf ("GPUelegant does not support particle output with csbend_csr.");
      //  long ip;
      //  /* dump particle data at this location */
      //  if (!SDDS_StartPage(&csbend->SDDSpart, n_part) ||
      //      !SDDS_SetParameters(&csbend->SDDSpart, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, 
      //                          "Pass", -1, "Kick", kick, "pCentral", Po, "Angle", phiBend, 
      //                          NULL))
      //    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
      //  convertFromCSBendCoords(part, n_part, rho0, cos_ttilt, sin_ttilt, 1);
      //  for (ip=0; ip<n_part; ip++) {
      //    if (!SDDS_SetRowValues(&csbend->SDDSpart, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE,
      //                           ip, 
      //                           csbend->xIndex, part[ip][0],
      //                           csbend->xpIndex, part[ip][1],
      //                           csbend->tIndex, part[ip][4],
      //                           csbend->pIndex, Po*(1+part[ip][5]),
      //                           -1)) 
      //      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
      //  }
      //  convertToCSBendCoords(part, n_part, rho0, cos_ttilt, sin_ttilt, 1);
      //  if (!SDDS_WritePage(&csbend->SDDSpart))
      //    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
      //  if (!inhibitFileSync)
      //    SDDS_DoFSync(&csbend->SDDSpart);
	}
      }

      if (tContext.sliceAnalysis && tContext.sliceAnalysis->active &&
	  kick!=(csbend->n_kicks-1) &&
	  (csbend->sliceAnalysisInterval==0 ||
	   kick%csbend->sliceAnalysisInterval==0)) {
#if (!USE_MPI)
	//convertFromCSBendCoords(part, n_part, rho0, cos_ttilt, sin_ttilt, 1);
	//performSliceAnalysisOutput(tContext.sliceAnalysis, part, n_part, 
	//			   0, tContext.step, Po, 
	//			   macroParticleCharge*n_part,
	//			   tContext.elementName, 
	//			   z_start + (kick*(z_end-z_start))/(csbend->n_kicks-1),
	//			   1);
	//convertToCSBendCoords(part, n_part, rho0, cos_ttilt, sin_ttilt, 1);
      if (isMaster) 
	printf ("GPUelegant does not support slice analysis output inside an element now.");
#else 
      if (isMaster) 
	printf ("Pelegant does not support slice analysis output inside an element now.");
    
#endif
      }

      if (csbend->wakeFileActive && 
          ((!csbend->outputLastWakeOnly && kick%csbend->outputInterval==0) ||
           (csbend->outputLastWakeOnly && kick==(csbend->n_kicks-1)))) {
        /* scale the linear density and its derivative to get C/s and C/s^2 
         * ctHist is already normalized to dct, but ctHistDeriv requires an additional factor
         */
        for (iBin=0; iBin<nBins; iBin++) {
          ctHist[iBin] *= macroParticleCharge*c_mks;
          ctHistDeriv[iBin] *= macroParticleCharge*sqr(c_mks)/dct;
        }
 
	if (isMaster) {
        if (!SDDS_StartPage(&csbend->SDDSout, nBins) ||
            !SDDS_SetColumn(&csbend->SDDSout, SDDS_SET_BY_NAME, dGamma, nBins, "DeltaGamma") ||
            !SDDS_SetColumn(&csbend->SDDSout, SDDS_SET_BY_NAME, T1, nBins, "DeltaGammaT1") ||
            !SDDS_SetColumn(&csbend->SDDSout, SDDS_SET_BY_NAME, T2, nBins, "DeltaGammaT2") ||
            !SDDS_SetColumn(&csbend->SDDSout, SDDS_SET_BY_NAME, ctHist, nBins, "LinearDensity") ||
            !SDDS_SetColumn(&csbend->SDDSout, SDDS_SET_BY_NAME, ctHistDeriv, nBins, "LinearDensityDeriv") ||
            !SDDS_SetParameters(&csbend->SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, 
                                "Pass", -1, "Kick", kick, "dsKick", csbend->length/csbend->n_kicks,
                                "pCentral", Po, "Angle", phiBend, "SlippageLength", slippageLength,
                                "TotalBunchLength", ctUpper-ctLower,
                                "BinSize", dct, 
                                "DerbenevRatio", derbenevRatio, NULL))
          SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
	}
        if (csbend->binOnce) {
          /* fix these arrays so they can be used again */
          ctHist[iBin] /= macroParticleCharge*c_mks;
          ctHistDeriv[iBin] /= macroParticleCharge*sqr(c_mks)/dct;
        }
        /* use T1 array to output s and T2 to output dGamma/ds */
        for (iBin=0; iBin<nBins; iBin++) {
          T1[iBin] = ctLower-(ctLower+ctUpper)/2.0+dct*(iBin+0.5);
          T2[iBin] = dGamma[iBin]/(csbend->length/csbend->n_kicks);
        }
	if (isMaster){
        if (!SDDS_SetColumn(&csbend->SDDSout, SDDS_SET_BY_NAME, T1, nBins, "s") ||
            !SDDS_SetColumn(&csbend->SDDSout, SDDS_SET_BY_NAME, T2, nBins, "GammaDeriv") ||
            !SDDS_WritePage(&csbend->SDDSout))
          SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
        if (!inhibitFileSync)
          SDDS_DoFSync(&csbend->SDDSout);
      }
    }
  }
  }

  if (!csbend->binOnce && n_partMoreThanOne && !csrInhibit && !csbend->csrBlock) {
    /* prepare some data for use by CSRDRIFT element */
    csrWake.dctBin = dct;
    ctLower = ctUpper = dct = 0;

    nBinned = gpu_binParticleCoordinate(ctHist, d_ctHist, &maxBins, &ctLower,
                &ctUpper, &dct, &nBins,
                csbend->binRangeFactor<1.1?1.1:csbend->binRangeFactor,
                n_part, 4);

#if (!USE_MPI)
    if (nBinned!=n_part) {
      fprintf(stdout, "Only %ld of %ld particles binned for CSRCSBEND (z0=%le, end, BRF=%le)\n", 
	      nBinned, n_part, z_start, csbend->binRangeFactor<1.1?1.1:csbend->binRangeFactor);
      fprintf(stdout, "ct min, max = %21.15e, %21.15e, dct = %21.15e, nBins=%ld, maxBins=%ld\n",
	      ctLower, ctUpper, dct, nBins, maxBins);
      fflush(stdout);
    }
#else
    if (USE_MPI && notSinglePart) {
      long all_binned, result = 1, nBinned_total;

      if (isSlave || !notSinglePart) {
	result = ((nBinned==n_part) ? 1 : 0);
      }
      else
	nBinned = 0;
      MPI_Allreduce(&result, &all_binned, 1, MPI_LONG, MPI_LAND, MPI_COMM_WORLD);
      MPI_Allreduce(&nBinned, &nBinned_total, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
      nBinned = nBinned_total; 
      if (!all_binned && isMaster) {
	fprintf(stdout, "Not all particles binned for CSRCSBEND (z0=%le, kick=%ld, BRF=%le)\n", 
		z_start, kick,
	      csbend->binRangeFactor<1.1?1.1:csbend->binRangeFactor);
      fprintf(stdout, "ct min, max = %21.15e, %21.15e, dct = %21.15e, nBins=%ld, maxBins=%ld\n",
	      ctLower, ctUpper, dct, nBins, maxBins);
      fflush(stdout);
      }
      if (notSinglePart) {  /* Master needs to know the information to write the result */
	buffer = (double*)malloc(sizeof(double) * nBins);
	MPI_Allreduce(ctHist, buffer, nBins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	memcpy(ctHist, buffer, sizeof(double)*nBins);
	free(buffer);
      }
    }
#endif
    csrWake.s0 = ctLower + dzf;
  } else {
    ctLower = ctUpper = dct = 0;
    csrWake.dctBin = dct;
    csrWake.s0 = ctLower + dzf;
  }

  if (isSlave || !notSinglePart) {
    /* remove lost particles, handle edge effects, and transform coordinates */    
    n_part = killParticles(n_part, d_sortIndex, accepted,
               gpu_track_through_csbendCSR_kernel7(d_sortIndex, Po,
               e2, z_start, expansionOrder1,
               csbend->edgeFlags, rho0, rad_coef));
    gpuErrorHandler("gpu_track_through_csbendCSR: gpu_track_through_csbendCSR_kernel7");
    if (csbend->edgeFlags&BEND_EDGE2_EFFECTS) {
      if (csbend->useMatrix)
        gpu_track_particles(Me2, n_part);
      else
        gpuDriver(n_part, gpu_track_through_csbendCSR_kernel8(rho_actual, rho0, n, e2,
                  psi2, he2, csbend->edge_order, csbend->edge2_effects, 
                  csbend->b[0], csbend->hgap, csbend->fint, csbend->h2));
        gpuErrorHandler("gpu_track_through_csbendCSR: gpu_track_through_csbendCSR_kernel8");
    }
    gpuDriver(n_part,
      gpu_track_through_csbendCSR_kernel9(dxf, dyf, dzf, 
      dcoord_etilt[0], dcoord_etilt[1], dcoord_etilt[2], dcoord_etilt[3],
      cos_ttilt, sin_ttilt));
    gpuErrorHandler("gpu_track_through_csbendCSR: gpu_track_through_csbendCSR_kernel9");
  }

  if (n_partMoreThanOne && !csbend->csrBlock) {
    /* prepare more data for CSRDRIFT */
    long imin, imax;
    double S55;

#if !USE_MPI   
    gpu_rms_emittance(0, 1, n_part, &csrWake.S11, &csrWake.S12, &csrWake.S22);
    gpu_rms_emittance(4, 5, n_part, &S55, NULL, NULL);
#else
    if (notSinglePart) {	
    	gpu_rms_emittance_p(0, 1, n_part, &csrWake.S11, &csrWake.S12, &csrWake.S22);
    	gpu_rms_emittance_p(4, 5, n_part, &S55, NULL, NULL);
    } else {
     	gpu_rms_emittance(0, 1, n_part, &csrWake.S11, &csrWake.S12, &csrWake.S22);
    	gpu_rms_emittance(4, 5, n_part, &S55, NULL, NULL);
    }
#endif

    csrWake.perc68BunchLength = 
      gpu_approximateBeamWidth(0.6826, d_ctHist, n_part, 4)/2;
    csrWake.perc90BunchLength = 
      gpu_approximateBeamWidth(0.9, d_ctHist, n_part, 4)/2;
	
    csrWake.rmsBunchLength = sqrt(S55);


#ifdef DEBUG
    fprintf(stderr, "rms bunch length = %le, percentile bunch length (68, 90) = %le, %le\n",
            csrWake.rmsBunchLength, csrWake.perc68BunchLength,
            csrWake.perc90BunchLength);
#endif
    if (macroParticleCharge) {
      index_min_max(&imin, &imax, csrWake.dGamma, csrWake.bins);
      csrWake.peakToPeakWavelength = 2*fabs(1.0*imax-imin)*dct;
    } else {
      csrWake.peakToPeakWavelength = csrWake.perc68BunchLength;
    }

    csrWake.valid = 1;
    csrWake.rho = rho_actual;
    csrWake.bendingAngle = accumulatingAngle ? fabs(phiBend) : fabs(angle);
    csrWake.Po = Po;
    csrWake.SGOrder = csbend->SGOrder;
    csrWake.SGDerivOrder = csbend->SGDerivOrder;
    csrWake.SGHalfWidth = csbend->SGHalfWidth;
    csrWake.SGDerivHalfWidth = csbend->SGDerivHalfWidth;
    csrWake.GSConstant = CSRConstant*pow(3*rho0*rho0, 1./3.)/2;  /* used for G. Stupakov's drift formulae */
    csrWake.MPCharge = macroParticleCharge;
    csrWake.binRangeFactor = csbend->binRangeFactor;
    csrWake.trapazoidIntegration = csbend->trapazoidIntegration;
    if (csbend->useMatrix) {
      free_matrices(Msection);
      free_matrices(Me1);
      free_matrices(Me2);
      free(Msection);
      free(Me1);
      free(Me2);
      Msection = Me1 = Me2 = NULL;
    }
  }

  if (csbend->csrBlock)
    accumulatedAngle = 0;
  else
    /* accumulate the bending angle just in case the same type of dipole follows */
    accumulatedAngle += fabs(angle);
    
#if defined(MINIMIZE_MEMORY)
  /* leave dGamma out of this because that memory is used by CSRDRIFT */
  free(ctHist);
  free(ctHistDeriv);
  free(T1);
  free(T2);
  free(denom);
  if (grnk)
    free(grnk);
  if (chik)
    free(chik);
  ctHist = ctHistDeriv = T1 = T2 = denom = NULL;
  maxBins = 0;
#endif

  cudaFree(d_ctHist);
  if (allocGpuMem) {
    cudaFree(d_dGamma);
  }

  return(n_part);
}
#undef DEBUG_IGF

} // extern "C"

long gpu_binParticleCoordinate(double *hist, double *d_hist, long *maxBins,
       double *lower, double *upper, double *binSize, long *bins,
       double expansionFactor, long nParticles, long coordinateIndex)
{
  long iBin, nBinned;

  if (*binSize<=0 && *bins<1)
    return -1;
  if (*binSize>0 && *bins>1)
    return -2;

  struct GPUBASE* gpuBase = getGpuBase();
  double *d_particles = gpuBase->d_particles;
  unsigned int particlePitch = gpuBase->gpu_array_pitch;

  /* if (*lower==*upper)  This condition will be removed */
  if (isSlave || !notSinglePart) {
    /* find range of points */
    gpuReduceMinMax(d_particles + particlePitch*coordinateIndex,
                    nParticles, lower, upper);
  }

#if USE_MPI
  /* find the global maximum and minimum */
  if (notSinglePart) {
    if (isMaster)
      nParticles = 0;
    find_global_min_max(lower, upper, nParticles, MPI_COMM_WORLD);
  }
#endif

  if (expansionFactor>1) {
    double center, range;
    center = (*lower+*upper)/2;
    range = (*upper-*lower)*expansionFactor;
    *lower = center-range/2;
    *upper = center+range/2;
  }

  if (*binSize>0)
    /* bin size given, so determine the number of bins */
    *bins = (*upper-*lower)/(*binSize);
  *binSize = (*upper-*lower)/(*bins);

  if (isMaster && USE_MPI) {
    for (iBin=0; iBin<*bins; iBin++)
      hist[iBin] = 0;
    nBinned = 0;
  } else {
    nBinned = gpu_binParticles_and_countBinned(d_hist, d_particles, particlePitch,
                nParticles, coordinateIndex, *lower, *binSize, *bins);
    cudaMemcpy(hist, d_hist, sizeof(double) * (*bins),
               cudaMemcpyDeviceToHost);
  }
  return nBinned;
}

class gpu_track_through_driftCSR_kernel1{
public:
  double dz, Po;
  // from CSRDRIFT
  int linearOptics;

  gpu_track_through_driftCSR_kernel1(double dz, double Po,
    int linearOptics) : dz(dz), Po(Po), linearOptics(linearOptics) {};

  __device__ inline void operator()(gpuParticleAccessor& coord){
    double p, beta;

    coord[0] += coord[1]*dz;
    coord[2] += coord[3]*dz;
    p = Po*(1+coord[5]);
    beta = p/sqrt(p*p+1);
    if (linearOptics)
      coord[4] = (coord[4]+dz)/beta;
    else
      coord[4] = (coord[4]+dz*sqrt(1+sqr(coord[1])+sqr(coord[3])))/beta;
  }
};

class gpu_track_through_driftCSR_kernel2{
public:
  double *d_dGamma;
  double Po, ctLower, dct, factor;
  unsigned int nBins1;

  gpu_track_through_driftCSR_kernel2(double *d_dGamma, double Po, 
    double ctLower, double dct, double factor, unsigned int nBins1) :
    d_dGamma(d_dGamma), Po(Po), ctLower(ctLower), dct(dct),
    factor(factor), nBins1(nBins1) {};

  __device__ unsigned int operator()(gpuParticleAccessor& coord){

    double f, p, beta;
    unsigned int iBin, binned(0);

    if(d_dGamma) {
      iBin = (f=(coord[4]-ctLower)/dct);
      f -= iBin;
      if (iBin<nBins1) {
        coord[5] += ((1-f)*d_dGamma[iBin] + f*d_dGamma[iBin+1])/Po*factor;
        binned = 1;
      }
      //} else {
      //  No message on GPU
      //  fprintf(stdout, "Particle out of bin range---not kicked: ct-ctLower=%21.15e, dct=%21.15e, iBin=%ld\n", coord[4]-ctLower, dct, iBin);
      //}
    }
    p = (1+coord[5])*Po;
    beta = p/sqrt(p*p+1);
    coord[4] = beta*coord[4];

    return binned;
  }
};

class gpu_track_through_driftCSR_kernel3{
public:
  double dz;
  // from CSRDRIFT
  int linearOptics;

  gpu_track_through_driftCSR_kernel3(double dz, int linearOptics) : 
    dz(dz), linearOptics(linearOptics) {};

  __device__ inline void operator()(gpuParticleAccessor& coord){

    coord[0] += dz*coord[1];
    coord[2] += dz*coord[3];
    if (linearOptics)
      coord[4] += dz;
    else
      coord[4] += dz*sqrt(1+sqr(coord[1])+sqr(coord[3]));
  }
};

extern "C" {

// From csbend.c
void computeSaldinFdNorm(double **FdNorm, double **x, long *n,
       double sMax, long ns, double Po, double radius, double angle,
       double dx, char *normMode);

long gpu_track_through_driftCSR(long np, CSRDRIFT *csrDrift, double Po,
       double **accepted, double zStart, double revolutionLength, 
       CHARGE *charge, char *rootname)
{
  long iKick, iBin, binned=0, nKicks, iSpreadMode=0;
  double dz, ct0=0.0, factor, dz0, dzFirst;
  double ctmin, ctmax, dct;
  double zTravel, attenuationLength, thetaRad=0.0, sigmaZ, overtakingLength, criticalWavelength, wavelength=0.0;
  static char *spreadMode[3] = {(char*)"full", (char*)"simple", (char*)"radiation-only"};
  static char *wavelengthMode[3] = {(char*)"sigmaz", (char*)"bunchlength", (char*)"peak-to-peak"};
  static char *bunchlengthMode[3] = {(char*)"rms", (char*)"68-percentile", (char*)"90-percentile"};
  unsigned long mode;
  static long warned = 0, incrementWarningsLeft=100;
  long nBins1;
  TRACKING_CONTEXT tContext;
#if USE_MPI 
  long np_total=1, np_tmp=np, binned_total;
#endif
  struct GPUBASE* gpuBase = getGpuBase();
  unsigned int particlePitch = gpuBase->gpu_array_pitch;
  double* d_dGamma   = gpuBase->d_temp_particles + 3*particlePitch;
  if ((isSlave || !notSinglePart) && particlePitch < csrWake.bins) {
    std::stringstream msg;
    msg << "driftCSR: GPUElegant requires more particles than bins!\n"
        << "particlePitch = " << particlePitch 
        << " number of bins = " << csrWake.bins;
    char* cmsg = new char[msg.str().size()+1];
    strcpy(cmsg, msg.str().c_str());
    bombElegant(cmsg, NULL);
  }
  
  getTrackingContext(&tContext);

#if (!USE_MPI)
  if (np<=1 || !csrWake.valid || !csrDrift->csr) {
#else
  if (notSinglePart){
    if (isMaster) 
      np_tmp = 0;  
    MPI_Allreduce(&np_tmp, &np_total, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);   
  } else
    np_total = np;

  if (np_total<=1 || !csrWake.valid || !csrDrift->csr) 	{
    if (isSlave||!notSinglePart) {
#endif
      if (csrDrift->linearOptics) {
        gpuDriver(np,
            gpu_track_through_csbendCSR_kernel1(csrDrift->length));
        gpuErrorHandler("gpu_track_through_driftCSR: gpu_track_through_csbendCSR_kernel1");
      }
      else
	gpu_exactDrift(np, csrDrift->length);
#if (USE_MPI)
    }
#endif
    return np;
  }	
  nBins1 = csrWake.bins - 1;

  mode = 
    (csrDrift->spread?CSRDRIFT_SPREAD:0) +
      (csrDrift->useOvertakingLength?CSRDRIFT_OVERTAKINGLENGTH:0) +
        (csrDrift->useSaldin54?CSRDRIFT_SALDIN54:0) +
          (csrDrift->attenuationLength>0?CSRDRIFT_ATTENUATIONLENGTH:0) +
            (csrDrift->useStupakov?CSRDRIFT_STUPAKOV:0) ;
  while (zStart+1.e-12<csrWake.zLast) {
    if (incrementWarningsLeft) {
      fprintf(stdout, "*** Warning: incrementing zStart by revolution length for CSRDRIFT (%s #%ld).\n",
              tContext.elementName, tContext.elementOccurrence);
      fprintf(stdout, "    If you are not simulating a ring, this could be a problem!\n");
      incrementWarningsLeft --;
    }
    zStart += revolutionLength;
  }
  if (bitsSet(mode)>1) {
    fprintf(stdout, "Error: Too many modes set for CSRDRIFT.\n");
    exitElegant(1);
  }
  if (csrWake.lastMode && csrWake.lastMode!=mode) {
    fprintf(stdout, "Error: CSRDRIFT mode changed between dipoles. Pick one mode following each dipole.\n");
    exitElegant(1);
  }
  csrWake.lastMode = mode;
  
  if (mode&CSRDRIFT_STUPAKOV)
    return gpu_track_through_driftCSR_Stupakov(np, csrDrift, Po, accepted, zStart, charge, rootname);

  if (!warned) {
    fprintf(stdout, "Warning: USE_STUPAKOV=1 is recommended for CSRDRIFT elements.\n");
    fprintf(stdout, "This is the most physical model available at this time in elegant.\n");
    warned = 1;
  }
  
  dct = csrWake.dctBin;
  if (csrDrift->dz>0) {
    if ((nKicks = csrDrift->length/csrDrift->dz)<1)
      nKicks = 1;
  } else 
    nKicks = csrDrift->nKicks;
  if (nKicks<=0)
    bombElegant("nKicks=0 in CSR drift.", NULL);
  dz = (dz0=csrDrift->length/nKicks)/2;
  
  sigmaZ = 0;
  switch (match_string(csrDrift->bunchlengthMode, bunchlengthMode, 3, 0)) {
  case 0:
    sigmaZ = csrWake.rmsBunchLength;
    break;
  case 1:
    sigmaZ = csrWake.perc68BunchLength;
    break;
  case 2:
    sigmaZ = csrWake.perc90BunchLength;
    break;
  default:
    bombElegant("invalid bunchlength_mode for CSRDRIFT.  Use rms or percentile.", NULL);
  }
  
  overtakingLength = pow(24*sigmaZ*csrWake.rho*csrWake.rho, 1./3.);

  if (mode&CSRDRIFT_OVERTAKINGLENGTH)
    attenuationLength = overtakingLength*csrDrift->overtakingLengthMultiplier;
  else
    attenuationLength = csrDrift->attenuationLength;
  
  if (mode&CSRDRIFT_SPREAD) {
    iSpreadMode = 0;
    if (csrDrift->spreadMode && 
        (iSpreadMode=match_string(csrDrift->spreadMode, spreadMode, 3, 0))<0)
      bombElegant("invalid spread_mode for CSR DRIFT.  Use full, simple, or radiation-only", NULL);
    switch (match_string(csrDrift->wavelengthMode, wavelengthMode, 3, 0)) {
    case 0:
    case 1:
      /* bunch length */
      wavelength = sigmaZ;
      break;
    case 2:
      /* peak-to-peak */
      wavelength = csrWake.peakToPeakWavelength;
      break;
    default:
      bombElegant("invalid wavelength_mode for CSR DRIFT.  Use sigmaz or peak-to-peak", NULL);
      break;
    }
    criticalWavelength = 4.19/pow3(csrWake.Po)*csrWake.rho;
    if (!particleIsElectron)
      bombElegant("CSRDRIFT spread mode is not supported for particles other than electrons", NULL);
    thetaRad = 0.5463e-3/(csrWake.Po*0.511e-3)/pow(criticalWavelength/wavelength, 1./3.);
  }

  if (mode&CSRDRIFT_SALDIN54) {
    if (csrWake.FdNorm==NULL) {
      if (csrDrift->nSaldin54Points<20) 
        csrDrift->nSaldin54Points = 20;
      computeSaldinFdNorm(&csrWake.FdNorm, &csrWake.xSaldin, &csrWake.nSaldin,
                          2*sigmaZ, csrDrift->nSaldin54Points, csrWake.Po, csrWake.rho, csrWake.bendingAngle, dz,
                          csrDrift->normMode);
      if (csrDrift->Saldin54Output)  {
        long ix;
        if (!csrDrift->fpSaldin) {
          csrDrift->Saldin54Output = compose_filename(csrDrift->Saldin54Output, rootname);
          csrDrift->fpSaldin = fopen(csrDrift->Saldin54Output, "w");
          fprintf(csrDrift->fpSaldin, "SDDS1\n&column name=z, type=double &end\n&column name=Factor, type=double &end\n");
          fprintf(csrDrift->fpSaldin, "&data mode=ascii no_row_counts=1 &end\n");
        } else
          fprintf(csrDrift->fpSaldin, "\n");
        for (ix=0; ix<csrWake.nSaldin; ix++) 
          fprintf(csrDrift->fpSaldin, "%le %le\n", csrWake.xSaldin[ix], csrWake.FdNorm[ix]);
        fflush(csrDrift->fpSaldin);
      }
    }
  }

  dzFirst = zStart - csrWake.zLast;
  zTravel = zStart-csrWake.z0;  /* total distance traveled by radiation to reach this point */
#ifdef DEBUG
  fprintf(stdout, "CSR in drift:\n");
  fprintf(stdout, "zStart = %21.15le, zLast = %21.15le, zTravel = %21.15le\n", zStart, csrWake.zLast,
          zTravel);
  fprintf(stdout, "dzFirst = %21.15e, s0 = %21.15e\n", dzFirst, csrWake.s0);
#endif

  for (iKick=0; iKick<nKicks; iKick++) {
    /* first drift is dz=dz0/2, others are dz0 */
    if (iKick==1)
      dz = dz0;
    zTravel += dz;

    ctmin = DBL_MAX;
    ctmax = -DBL_MAX;

    /* propagate particles forward, converting s to c*t=s/beta */
    if (isSlave || !notSinglePart) {
      gpuDriver(np,
          gpu_track_through_driftCSR_kernel1(dz, Po, csrDrift->linearOptics));
        gpuErrorHandler("gpu_track_through_driftCSR: gpu_track_through_driftCSR_kernel1");
#ifdef DEBUG
      gpuReduceMinMax(getGpuBase()->d_particles + 4*particlePitch, np, &ctmin, &ctmax);
#endif
    }

    factor = 1;
    if (csrWake.dGamma) {
      /* propagate wake forward */
      csrWake.s0 += dz+dzFirst;   /* accumulates position of back end of the radiation pulse */
      ct0 = csrWake.s0;
      
      if (attenuationLength>0) {
        /* attenuate wake */
        if ((factor = exp(-(dz+dzFirst)/attenuationLength))<1) {
          for (iBin=0; iBin<csrWake.bins; iBin++)
            csrWake.dGamma[iBin] *= factor;
        }
      }
      /* factor to account for difference in drift lengths here and in
       * csrcsbend integration.  Use dz0 here because that is the
       * length integrated by each kick.  Add dzFirst to account for any
       * length we may have missed due to intervening non-drift elements.
       */
      factor = (dz0+dzFirst)/csrWake.ds0;
    }
    if (mode&CSRDRIFT_SPREAD) {
      /* compute loss of on-axis field due to spread of beam using a simple-minded
       * computation of beam sizes */
      switch (iSpreadMode) {
      case 0:  /* full */
        factor *= sqrt(csrWake.S11/(csrWake.S11 + 
                                     2*zTravel*csrWake.S12 + 
                                     zTravel*zTravel*(sqr(thetaRad)+csrWake.S22)));
        break;
      case 1: /* simple */
        factor *= sqrt(csrWake.S11/(csrWake.S11 + zTravel*zTravel*(sqr(thetaRad)+csrWake.S22)));
        break;
      case 2: /* radiation only */
        factor *= sqrt(csrWake.S11/(csrWake.S11 + sqr(zTravel*thetaRad)));
        break;
      default:
        bombElegant("invalid spread code---programming error!", NULL);
        break;
      }
    }
    
    if (mode&CSRDRIFT_SALDIN54) {
      long code;
      double f0 = 0;
      if (zTravel<=csrWake.xSaldin[csrWake.nSaldin-1]) 
        factor *= (f0=interp(csrWake.FdNorm, csrWake.xSaldin, csrWake.nSaldin, zTravel, 0, 1, &code));
      else 
        factor = 0;
      csrWake.lastFdNorm = f0;
#ifdef DEBUG
      fprintf(csrWake.fpSaldin, "%le %le\n", zTravel, f0);
      fflush(csrWake.fpSaldin);
#endif
      if (!code) {
        fprintf(stderr, "Warning: interpolation failure for Saldin eq. 54\n");
        fprintf(stderr, "zTravel = %le,  csrWake available up to %le\n",
                zTravel, csrWake.xSaldin[csrWake.nSaldin-1]);
        factor = 0;
      }
    }
    
    dzFirst = 0;

    /* apply kick to each particle and convert back to normal coordinates */
    if (isSlave || !notSinglePart) {
      if (csrWake.dGamma)
        cudaMemcpy(d_dGamma, (void*)csrWake.dGamma,
                   sizeof(double)*csrWake.bins, cudaMemcpyHostToDevice);
      else
        d_dGamma = NULL;
      binned = gpuUnsignedIntParticleReduction(np,
                 gpu_track_through_driftCSR_kernel2(d_dGamma, Po,
                 ct0, dct, factor, nBins1), Add<unsigned int>() );
    }
#if USE_MPI
    if (isSlave && notSinglePart) {
      MPI_Allreduce(&binned, &binned_total, 1, MPI_LONG, MPI_SUM, workers);
    }
    if ((myid==1) && (csrWake.dGamma && np_total!=binned_total)) {
      dup2(fd,fileno(stdout)); /* Let the first slave processor write the output */
      fprintf(stdout, "only %ld of %ld particles binned for CSR drift %s (track_through_driftCSR)\n",
              binned_total, np_total, tContext.elementName);
#else
    if (csrWake.dGamma && np!=binned) {
      fprintf(stdout, "only %ld of %ld particles binned for CSR drift %s (track_through_driftCSR)\n",
              binned, np, tContext.elementName);
#endif
      fprintf(stdout, "beam ct min, max = %21.15e, %21.15e\n",
              ctmin, ctmax);
      fprintf(stdout, "wake ct0 = %21.15e, ct1 = %21.15e\n",
              ct0, ct0+csrWake.dctBin*csrWake.bins);
      fflush(stdout);
#if USE_MPI
#if defined(_WIN32)
    freopen("NUL","w",stdout); 
#else
    freopen("/dev/null","w",stdout); 
#endif
#endif  
    }
  }
  /* do final drift of dz0/2 */
  dz = dz0/2;
  if (isSlave || !notSinglePart) {
    gpuDriver(np,
              gpu_track_through_driftCSR_kernel3(dz, csrDrift->linearOptics));
    gpuErrorHandler("gpu_track_through_driftCSR: gpu_track_through_driftCSR_kernel3");
  }
  csrWake.zLast = zStart+csrDrift->length;
  
  if (csrWake.dGamma) {
    /* propagate wake forward */
    csrWake.s0 += dz;
    ct0 = csrWake.s0;
    
    if (attenuationLength>0) {
      /* attenuate wake */
      if ((factor = exp(-dz/attenuationLength))<1) {
        for (iBin=0; iBin<csrWake.bins; iBin++)
            csrWake.dGamma[iBin] *= factor;
        }
    }
  }

  return np;
}

void gpu_exactDrift(long np, double length){
  gpuDriver(np, gpuExactDrift(length) );
  gpuErrorHandler("gpu_exactDrift: gpuExactDrift");
}

// From csbend.c
double SolveForPhiStupakov(double x, double ds, double phim);
void DumpStupakovOutput(char *filename, SDDS_DATASET *SDDSout, long *active,
                        double zTravel, double *ctHist, double *ctHistDeriv,
                        double *dGamma, long nBins, double dct, 
                        double MPCharge, double dz,
                        long nCaseC, long nCaseD1,long nCaseD2,
                        double x, double dsMax, double phi0, double phi1);

// From csbend.c
extern double SolveForPhiStupakovDiffSum;
extern long SolveForPhiStupakovDiffCount;

long gpu_track_through_driftCSR_Stupakov(long np, CSRDRIFT *csrDrift, 
       double Po, double **accepted, double zStart, CHARGE *charge, char *rootname)
{
  long iKick, iBin, binned=0, nKicks;
  long nCaseC, nCaseD1, nCaseD2;
  double ctLower, ctUpper, ds;
  long nBins, maxBins, nBinned, diBin;
  double dz, factor, dz0, dzFirst;
  double zTravel, dct, zOutput;
  double *ctHist=NULL, *ctHistDeriv=NULL, *phiSoln=NULL;
  double length;
  long nBins1;
  double dsMax, x;
  TRACKING_CONTEXT tContext;
  LSCKICK lscKick;
#if USE_MPI
  long binned_total=1, np_total=1;
  double *buffer;
#endif

  struct GPUBASE* gpuBase = getGpuBase();
  unsigned int particlePitch = gpuBase->gpu_array_pitch;
  double* d_ctHist = gpuBase->d_temp_particles + 2*particlePitch;
  double* d_dGamma = gpuBase->d_temp_particles + 3*particlePitch;
  if ((isSlave || !notSinglePart) && particlePitch < csrWake.bins) {
    std::stringstream msg;
    msg << "driftCSR_Stupakov: GPUElegant requires more particles than bins!\n"
        << "particlePitch = " << particlePitch 
        << " number of bins = " << csrWake.bins;
    char* cmsg = new char[msg.str().size()+1];
    strcpy(cmsg, msg.str().c_str());
    bombElegant(cmsg, NULL);
  }

  getTrackingContext(&tContext);

  SolveForPhiStupakovDiffCount = 0;
  SolveForPhiStupakovDiffSum = 0;
  
  length = csrDrift->length;
  if (zStart!=csrWake.zLast) {
    length += (dzFirst = zStart-csrWake.zLast);
    /* propagate beam back so we can tranverse the missing length including CSR
     */
    if (isSlave || !notSinglePart) {
      gpuDriver(np,
                gpu_track_through_driftCSR_kernel3(-dzFirst,
                csrDrift->linearOptics));
      gpuErrorHandler("gpu_track_through_driftCSR_Stupakov: gpu_track_through_driftCSR_kernel3");
    }
    zStart = csrWake.zLast;
  }
  zOutput = zStart;  /* absolute coordinate used for output of data vs z or s */
  
  if (csrDrift->dz>0) {
    if ((nKicks = length/csrDrift->dz+0.5)<1)
      nKicks = 1;
  } else 
    nKicks = csrDrift->nKicks;
  if (nKicks<=0)
    bombElegant("nKicks=0 in CSR drift.", NULL);
  dz = (dz0=length/nKicks)/2;
  
  zTravel = zStart-csrWake.z0;  /* total distance traveled by radiation to reach this point */

  maxBins = nBins = csrWake.bins;
  nBins1 = nBins-1;
  if (!(ctHist=(double*)SDDS_Malloc(sizeof(*ctHist)*nBins)) ||
      !(ctHistDeriv=(double*)SDDS_Malloc(sizeof(*ctHistDeriv)*nBins)) ||
      !(phiSoln=(double*)SDDS_Malloc(sizeof(*phiSoln)*nBins)))
    bombElegant("memory allocation failure (track_through_driftCSR)", NULL);

  if ((lscKick.bins = csrDrift->LSCBins)>0) {
    lscKick.interpolate = csrDrift->LSCInterpolate;
    lscKick.radiusFactor = csrDrift->LSCRadiusFactor;
    lscKick.lowFrequencyCutoff0 = csrDrift->LSCLowFrequencyCutoff0;
    lscKick.lowFrequencyCutoff1 = csrDrift->LSCLowFrequencyCutoff1;
    lscKick.highFrequencyCutoff0 = csrDrift->LSCHighFrequencyCutoff0;
    lscKick.highFrequencyCutoff1 = csrDrift->LSCHighFrequencyCutoff1;
  } 
  for (iKick=0; iKick<nKicks; iKick++) {
    /* first drift is dz=dz0/2, others are dz0 */
    if (iKick==1)
      dz = dz0;
    zTravel += dz;
    zOutput += dz;
    
    x = zTravel/csrWake.rho;
    dsMax = csrWake.rho/24*pow(csrWake.bendingAngle, 3)
      *(csrWake.bendingAngle+4*x)/(csrWake.bendingAngle+x);
    /* propagate particles forward, converting s to c*t=s/beta */
    if (isSlave || !notSinglePart) {
      gpuDriver(np,
        gpu_track_through_driftCSR_kernel1(dz, Po, csrDrift->linearOptics));
      gpuErrorHandler("gpu_track_through_driftCSR_Stupakov: gpu_track_through_driftCSR_kernel1");
    }
    /* bin the particle distribution */
    ctLower = ctUpper = dct = 0;
    nBinned = gpu_binParticleCoordinate(ctHist, d_ctHist, &maxBins, &ctLower,
                &ctUpper, &dct, &nBins,
                csrWake.binRangeFactor<1.1?1.1:csrWake.binRangeFactor,
                np, 4);

#if USE_MPI
  if (notSinglePart) {
    if (isSlave)
      MPI_Allreduce(&np, &np_total, 1, MPI_LONG, MPI_SUM, workers);    
    MPI_Allreduce(&nBinned, &binned_total, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  }

  if (notSinglePart) {  /* Master needs to know the information to write the result */
    buffer = (double*)malloc(sizeof(double) * nBins);
    MPI_Allreduce(ctHist, buffer, nBins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    memcpy(ctHist, buffer, sizeof(double)*nBins);
    free(buffer);
  }
  if ((myid==1) && (np_total!=binned_total)) {
    dup2(fd,fileno(stdout)); /* Let the first slave processor write the output */
    fprintf(stdout, "Only %ld of %ld particles binned for CSRDRIFT (%s, BRF=%le, Stupakov)\n", 
	    binned_total, np_total,
	    tContext.elementName, csrWake.binRangeFactor);
    fflush(stdout);
#else
  if (nBinned!=np) {
    fprintf(stdout, "Only %ld of %ld particles binned for CSRDRIFT (%s, BRF=%le, Stupakov)\n", 
	    nBinned, np,
	    tContext.elementName, csrWake.binRangeFactor);
#endif
    fprintf(stdout, "ct min, max = %21.15e, %21.15e, dct = %21.15e, nBins=%ld, maxBins=%ld\n",
	    ctLower, ctUpper, dct, nBins, maxBins);
    fflush(stdout);
#if USE_MPI
#if defined(_WIN32)
    freopen("NUL","w",stdout); 
#else
    freopen("/dev/null","w",stdout); 
#endif
#endif 

    if (csrDrift->LSCBins>0)
      gpu_addLSCKick(gpuBase->d_particles, particlePitch, np, &lscKick, Po, 
                     charge, dz, 0.0, gpuBase->d_temp_particles);
  }
      
    /* - smooth the histogram, normalize to get linear density, and 
       copy in preparation for taking derivative
       */
    if (csrWake.highFrequencyCutoff0>0 || csrWake.lowFrequencyCutoff0>=0) {
      long nz;
      nz = applyLHPassFilters(ctHist, nBins, 
                              csrWake.lowFrequencyCutoff0, csrWake.lowFrequencyCutoff1,
                              csrWake.highFrequencyCutoff0, csrWake.highFrequencyCutoff1,
                              csrWake.clipNegativeBins);
      if (nz && negativeWarningsLeft) {
        fprintf(stdout, "Warning: low pass filter resulted in negative values in %ld bins\n",
                nz);
        if (--negativeWarningsLeft==0)
          fprintf(stdout, "         Further warnings will be suppressed for this run.\n");
        fflush(stdout);
      }
    }

    if (csrWake.SGHalfWidth>0) {
      SavitzkyGolaySmooth(ctHist, nBins, csrWake.SGOrder, csrWake.SGHalfWidth, csrWake.SGHalfWidth,  0);
#if (!USE_MPI)
      correctDistribution(ctHist, nBins, 1.0*nBinned);
#else
      if (notSinglePart)
	correctDistribution(ctHist, nBins, 1.0*binned_total);
      else
	correctDistribution(ctHist, nBins, 1.0*nBinned);
#endif
    }
    for (iBin=0; iBin<nBins; iBin++)
      ctHistDeriv[iBin] = (ctHist[iBin] /= dct);
    /* - compute derivative with smoothing.  The deriv is w.r.t. index number and
     * I won't scale it now as it will just fall out in the integral 
     */
    SavitzkyGolaySmooth(ctHistDeriv, nBins, csrWake.SGDerivOrder, 
                       csrWake.SGDerivHalfWidth, csrWake.SGDerivHalfWidth, 1);

    /* Case C */ 
    nCaseC = 0;
    nCaseD1 = 0;
    nCaseD2 = 0;
    for (iBin=0; iBin<nBins; iBin++) {
      double f;
      ds = csrWake.rho/6*sqr(csrWake.bendingAngle)*(csrWake.bendingAngle + 3*x);
      diBin = ds/dct;
      if (iBin+diBin<nBins) {
        f = -1/(csrWake.bendingAngle+2*x); 
        csrWake.dGamma[iBin] = f*ctHist[iBin+diBin];
        nCaseC++;
      } else
        csrWake.dGamma[iBin] = 0;
    }
    /* Case D */
    for (iBin=0; iBin<nBins; iBin++) {
      phiSoln[iBin] = -1;
      if ((ds = iBin*dct)>dsMax)
        break;
      phiSoln[iBin] = SolveForPhiStupakov(x, iBin*dct/csrWake.rho, csrWake.bendingAngle);
    }
    for (iBin=0; iBin<nBins; iBin++) {
      long jBin, first, count;
      double term1=0, term2=0;
      diBin = dsMax/dct;
      if (iBin+diBin<nBins) {
        nCaseD1 ++;
        csrWake.dGamma[iBin] += ctHist[iBin+diBin]/(csrWake.bendingAngle+2*x);
      }
      first = 1;
      count = 0;
      for (jBin=iBin; jBin<nBins; jBin++) {
        double phi;
        if ((phi = phiSoln[jBin-iBin])>=0) {
          /* I put in a negative sign here because my s is opposite in direction to 
           * Saldin et al. and Stupakov, so my derivative has the opposite sign.
           * Note lack of ds factor here as I use the same one in my unnormalized derivative.
           */
          if (phi>0) {
            /* ^^^ If I test phi+2*x here, I get noisy, unphysical results very close
             * to the dipole exit 
             */
            term2 = ctHistDeriv[jBin]/(phi+2*x);
            csrWake.dGamma[iBin] -= term2;
            if (first) {
              term1 = term2;
              first = 0;
            }
            count++;
            nCaseD2++;
          }
        } else
          break;
      }
      if (count>1 && csrWake.trapazoidIntegration)
        /* trapazoid rule correction for ends */
        csrWake.dGamma[iBin] += (term1+term2)/2;
    }
    /* the minus sign adjusts for Stupakov using wake<0 to indicate energy gain
     */
    factor = -4/csrWake.rho*csrWake.GSConstant*dz0;
    for (iBin=0; iBin<nBins; iBin++)
      csrWake.dGamma[iBin] *= factor;

    if (csrWake.wffValues) 
      applyFilterTable(csrWake.dGamma, nBins, dct/c_mks, csrWake.wffValues, csrWake.wffFreqValue,
                       csrWake.wffRealFactor, csrWake.wffImagFactor);

    if ((csrDrift->StupakovOutput || csrWake.StupakovFileActive) && 
        (csrDrift->StupakovOutputInterval<2 || iKick%csrDrift->StupakovOutputInterval==0)) {
      double x, dsMax, phi0, phi1;
      if (!csrWake.StupakovFileActive) {
        if (!SDDS_CopyString(&csrWake.StupakovOutput, csrDrift->StupakovOutput))
          bombElegant("string copying problem preparing Stupakov output for CSRDRIFT", NULL);
        csrWake.StupakovOutput = compose_filename(csrWake.StupakovOutput, rootname);
      }
      x = zTravel/csrWake.rho;
      dsMax = csrWake.rho/24*pow(csrWake.bendingAngle, 3)
        *(csrWake.bendingAngle+4*x)/(csrWake.bendingAngle+x);
      phi0 = SolveForPhiStupakov(x, 0.0, csrWake.bendingAngle);
      phi1 = SolveForPhiStupakov(x, dsMax/csrWake.rho*0.999, csrWake.bendingAngle);
      
      /* note that the contents of ctHist and ctHistDeriv are corrupted by this operation */
      DumpStupakovOutput(csrWake.StupakovOutput, &csrWake.SDDS_Stupakov, 
                         &csrWake.StupakovFileActive, zTravel,
                         ctHist, ctHistDeriv, csrWake.dGamma, nBins, dct, csrWake.MPCharge,
                         dz0, nCaseC, nCaseD1, nCaseD2,
                         x, dsMax/csrWake.rho, phi0, phi1);
    }
    
    /* apply kick to each particle and convert back to normal coordinates */
    if (isSlave || !notSinglePart) {
      cudaMemcpy(d_dGamma, (void*)csrWake.dGamma, sizeof(double)*nBins,
                 cudaMemcpyHostToDevice);
      binned = gpuUnsignedIntParticleReduction(np,
                 gpu_track_through_driftCSR_kernel2(d_dGamma, Po,
                 ctLower, dct, 1.0, nBins1), Add<unsigned int>() );
    }

    if (tContext.sliceAnalysis && tContext.sliceAnalysis->active &&
	(csrDrift->sliceAnalysisInterval==0 ||
	 iKick%csrDrift->sliceAnalysisInterval==0)) {
#if USE_MPI
      /* This function will be parallelized in the future */
      fprintf(stdout, "performSliceAnalysisOutput is not supported in parallel mode currently.\n");
      MPI_Barrier(MPI_COMM_WORLD); /* Make sure the information can be printed before aborting */
      MPI_Abort(MPI_COMM_WORLD, 1); 
#endif
      // TODO
      bombElegant("performSliceAnalysisOutput is not supported on the GPU", NULL);
	//performSliceAnalysisOutput(tContext.sliceAnalysis, part, np, 
	//			   0, tContext.step, Po, 
	//			   csrWake.MPCharge*np,
	//			   tContext.elementName, 
	//			   zOutput, 0);
    }
#if USE_MPI
    if (isSlave && notSinglePart) {
      MPI_Allreduce(&binned, &binned_total, 1, MPI_LONG, MPI_SUM, workers);
    }
    if ((myid==1) && (np_total!=binned_total)) {
      dup2(fd,fileno(stdout)); /* Let the first slave processor write the output */
      fprintf(stdout, "Only %ld of %ld particles kicked for CSRDRIFT (%s, BRF=%le, Stupakov)\n", 
	      binned_total, np_total,
	      tContext.elementName, csrWake.binRangeFactor);
#else
    if (np!=binned) {
      fprintf(stdout, "Only %ld of %ld particles kicked for CSRDRIFT (%s, BRF=%le, Stupakov)\n", 
	      binned, np,
	      tContext.elementName, csrWake.binRangeFactor);
#endif
      fprintf(stdout, "ct min, max = %21.15e, %21.15e, dct = %21.15e, nBins=%ld, maxBins=%ld\n",
	      ctLower, ctUpper, dct, nBins, maxBins);
      fflush(stdout);
#if USE_MPI
#if defined(_WIN32)
    freopen("NUL","w",stdout); 
#else
    freopen("/dev/null","w",stdout); 
#endif
#endif 
    }
  }
  
  /* do final drift of dz0/2 */
  dz = dz0/2;
  if (isSlave || !notSinglePart) {
    gpuDriver(np,
              gpu_track_through_driftCSR_kernel3(dz,
              csrDrift->linearOptics));
    gpuErrorHandler("gpu_track_through_driftCSR_Stupakov: gpu_track_through_driftCSR_kernel3");
  }

  if (csrDrift->LSCBins>0)
    gpu_addLSCKick(gpuBase->d_particles, particlePitch, np, &lscKick, Po,
                   charge, dz, 0.0, gpuBase->d_temp_particles);

  csrWake.zLast = zStart + length;
  free(ctHist);
  free(ctHistDeriv);
  free(phiSoln);
#if DEBUG
  if (SolveForPhiStupakovDiffCount)
    fprintf(stdout, "Phi solution accuracy for %ld solutions: %le\n",
            SolveForPhiStupakovDiffCount, SolveForPhiStupakovDiffSum/SolveForPhiStupakovDiffCount);
#endif
  return np;
}

} // extern "C"

__device__ void 
gpu_apply_edge_effects(double *x, double *xp, double *y, double *yp, 
                       double rho, double n, double beta, double he,
                       double psi, int which_edge)
{
  double h, tan_beta, tan2_beta, sec_beta, sec2_beta, h2;
  double R21, R43;
  double T111, T133, T211, T441, T331, T221, T233, T243, T431, T432;
  double x0, xp0, y0, yp0;

  h = 1/rho;
  R21 = h*(tan_beta=tan(beta));
  R43 = -h*tan(beta-psi);

  h2 = sqr(h);
  T111 = which_edge*h/2*(tan2_beta=sqr(tan_beta));
  T133 = -which_edge*h/2*(sec2_beta=sqr(sec_beta=1./cos(beta)));
  T211 = which_edge==-1?
    -n*h2*tan_beta:
    -h2*(n+tan2_beta/2)*tan_beta;
  T441 = -(T331 = T221 = -which_edge*h*tan2_beta);
  T233 =  which_edge==-1?
    h2*(n+.5+tan2_beta)*tan_beta:
    h2*(n-tan2_beta/2)*tan_beta;
  T243 = which_edge*h*tan2_beta;
  T431 = h2*(2*n+(which_edge==1?sec2_beta:0))*tan_beta;
  T432 = which_edge*h*sec2_beta;
  if (he!=0) {
    double term;
    term = h/2*he*sec2_beta*sec_beta;
    T211 += term;
    T233 -= term;
    T431 -= 2*term;
  }

  x0 = *x;  xp0 = *xp;  y0 = *y;  yp0 = *yp;
  *x  = x0  + T111*sqr(x0) + T133*sqr(y0);
  *xp = xp0 + R21*x0 + T211*sqr(x0) + T221*x0*xp0 + T233*sqr(y0) + T243*y0*yp0;
  *y  = y0  + T331*x0*y0;
  *yp = yp0 + R43*y0 + T441*yp0*x0 + T431*x0*y0 + T432*xp0*y0;
}

__inline__ __device__ double ipow(double x, long p)
{
  register double hp;
  register long n;

  if (!x)  {
    return(p==0?1.:0.);
  }

  if (p<0)
    return(1./ipow(x, -p));

  switch (p) {
  case 0:
    return(1.);
  case 1:
    return(x);
  case 2:
    return(x*x);
  case 3:
    hp = x*x;
    return(hp*x);
  case 4:
    hp = x*x;
    return(hp*hp);
  case 5:
    hp = x*x;
    return(hp*hp*x);
  case 6:
    hp = x*x;
    return(hp*hp*hp);
  case 7:
    hp = x*x*x;
    return(hp*hp*x);
  case 8:
    hp = x*x;
    hp = hp*hp;
    return(hp*hp);
  default:
    n = p/2;
    hp = ipow(x, n);
    switch (p-2*n) {
    case 0:
      return(hp*hp);
    case 1:
      return(hp*hp*x);
    }
    break;
  }
  return(0.);  /* never actually executed--keeps compiler happy */
}

/* dipole fringe effects symplectic tracking, based on work of Kilean Hwang */
__device__ void 
gpu_dipoleFringeSym(double *x, double *xp, double *y, double *yp,
                    double *dp, double rho, double inFringe, long edgeOrder, 
                    double K1, double edge, double gap, double fint, double Rhe)
{
  double dx, dpx, dy, dpy;
  double tan_edge, sin_edge, sec_edge, cos_edge;
  double x0, xp0, y0, yp0, dp0;
  /*double psi;*/
  double k0, k3, k2;
  /*double Kg;*/
  double k4, k5, k6;

  k0 = sqr(PI)/6.;
  k2 = fint;
  k3 = 1.0*1./6.;
  /*Kg = gap*fint;*/
  k4 = -1.0*sqr(PI)/3.;
  k5 = 0.0;
  k6 = -1.0;

  x0 = *x;  xp0 = *xp;  y0 = *y;  yp0 = *yp; dp0 = *dp;
  dx = dpx = dy = dpy = 0;
  /*psi = Kg/rho/cos(edge)*(1+sqr(sin(edge)));*/

  sec_edge=1./cos(edge);
  tan_edge=tan(edge);
  sin_edge=sin(edge);
  cos_edge=cos(edge);
  

  if (edgeOrder>1) {

    /* entrance */
    if (inFringe==-1.) {
      dx  =   inFringe*ipow(sec_edge,2)*ipow(gap,2)*k0/rho/(1+dp0)
        + inFringe*ipow(x0,2)*ipow(tan_edge,2)/2/rho/(1+dp0) 
          - inFringe*ipow(y0,2)*ipow(sec_edge,2)/2/rho/(1+dp0);
      dy  =  -inFringe*x0*y0*ipow(tan_edge,2)/rho/(1+dp0);
      dpx  =  -1.*ipow(sec_edge,3)*sin_edge*ipow(gap,2)*k0/rho/rho/(1+dp0)
        +tan_edge*x0/rho
          +ipow(y0,2)/2*(2*ipow(tan_edge,3))/ipow(rho,2)/(1+dp0)
            +ipow(y0,2)/2*(ipow(tan_edge,1))/ipow(rho,2)/(1+dp0)
              -inFringe*(x0*xp0-y0*yp0)*ipow(tan_edge,2)/rho
		+k4*ipow(sin_edge,2)*ipow(gap,2)/2/ipow(cos_edge,3)/rho*Rhe
                  -k5*x0*ipow(sin_edge,1)*ipow(gap,1)/ipow(cos_edge,3)/rho*Rhe
                    +k6*(y0*y0-x0*x0)/2/ipow(cos_edge,3)/rho*Rhe;
      dpy  =  -1.*tan_edge*y0/rho 
        +k2*y0*(1+ipow(sin_edge,2))*gap/(1+dp0)/ipow(rho,2)/ipow(cos_edge,3)
          +inFringe*(x0*yp0+y0*xp0)*ipow(tan_edge,2)/rho
            +inFringe*y0*xp0/rho
              +k3*ipow(y0,3)*(2./3./cos_edge-4./3./ipow(cos_edge,3))/(1+dp0)/rho/rho/gap
		+k6*x0*y0/ipow(cos_edge,3)/rho*Rhe;
    }
    /* exit */
    if (inFringe==1.) {
      dx  =   inFringe*ipow(sec_edge,2)*ipow(gap,2)*k0/rho/(1+dp0)
        + inFringe*ipow(x0,2)*ipow(tan_edge,2)/2/rho/(1+dp0) 
          - inFringe*ipow(y0,2)*ipow(sec_edge,2)/2/rho/(1+dp0);
      dy  =  -inFringe*x0*y0*ipow(tan_edge,2)/rho/(1+dp0);
      dpx  =  tan_edge*x0/rho
          -ipow(y0,2)/2*(1*ipow(tan_edge,3))/ipow(rho,2)/(1+dp0)
            -ipow(x0,2)/2*(1*ipow(tan_edge,3))/ipow(rho,2)/(1+dp0)
              -inFringe*(x0*xp0-y0*yp0)*ipow(tan_edge,2)/rho
		+k4*ipow(sin_edge,2)*ipow(gap,2)/2/ipow(cos_edge,3)/rho*Rhe
                  -k5*x0*ipow(sin_edge,1)*ipow(gap,1)/ipow(cos_edge,3)/rho*Rhe
                    +k6*(y0*y0-x0*x0)/2/ipow(cos_edge,3)/rho*Rhe;
      dpy  =  -1.*tan_edge*y0/rho 
        +k2*y0*(1+ipow(sin_edge,2))*gap/(1+dp0)/ipow(rho,2)/ipow(cos_edge,3)
          +inFringe*(x0*yp0+y0*xp0)*ipow(tan_edge,2)/rho
            +inFringe*y0*xp0/rho
              +x0*y0*ipow(sec_edge,2)*tan_edge/ipow(rho,2)/(1+dp0)
		+k3*ipow(y0,3)*(2./3./cos_edge-4./3./ipow(cos_edge,3))/(1+dp0)/rho/rho/gap
                  -k5*y0*ipow(sin_edge,1)*ipow(gap,1)/ipow(cos_edge,3)/rho*Rhe
                    +k6*x0*y0/ipow(cos_edge,3)/rho*Rhe;
    }
    
  } else {
    /* linear terms in transverse coordinates only */

    /* entrance */
    if (inFringe==-1.) {
      dx  =   inFringe*ipow(sec_edge,2)*ipow(gap,2)*k0/rho/(1+dp0);
      dy  =  0;
      dpx  =  -1.*ipow(sec_edge,3)*sin_edge*ipow(gap,2)*k0/rho/rho/(1+dp0)
        +tan_edge*x0/rho
          +k4*ipow(sin_edge,2)*ipow(gap,2)/2/ipow(cos_edge,3)/rho*Rhe
            -k5*x0*ipow(sin_edge,1)*ipow(gap,1)/ipow(cos_edge,3)/rho*Rhe;
      dpy  =  -1.*tan_edge*y0/rho 
        +k2*y0*(1+ipow(sin_edge,2))*gap/(1+dp0)/ipow(rho,2)/ipow(cos_edge,3);
    }

    /* exit */
    if (inFringe==1.) {
      dx  =   inFringe*ipow(sec_edge,2)*ipow(gap,2)*k0/rho/(1+dp0);
      dy  =  0;
      dpx  =  tan_edge*x0/rho
        +k4*ipow(sin_edge,2)*ipow(gap,2)/2/ipow(cos_edge,3)/rho*Rhe
          -k5*x0*ipow(sin_edge,1)*ipow(gap,1)/ipow(cos_edge,3)/rho*Rhe;
      dpy  =  -1.*tan_edge*y0/rho 
        +k2*y0*(1+ipow(sin_edge,2))*gap/(1+dp0)/ipow(rho,2)/ipow(cos_edge,3)
          -k5*y0*ipow(sin_edge,1)*ipow(gap,1)/ipow(cos_edge,3)/rho*Rhe;
    }
  }
  
  *x  = x0  + dx;
  *xp = xp0 + dpx/(1+dp0);
  *y  = y0  + dy;
  *yp = yp0 + dpy/(1+dp0);
  /*  printf("x %f y %f xp %f yp %f dp0 %f\n", *x, *y, *xp, *yp, dp0); */
}


//TODO: add for particle output in CSBEND_CSR
/* this is used solely to convert coordinates inside the element for
 * the purpose of generating output.  It ignores misalignments.
 */
//void gpu_convertFromCSBendCoords(long np, double rho0, double cos_ttilt,
//                                 double sin_ttilt, long ctMode)
//{
//  long ip;
//  double x, y, xp, yp, *coord;
//
//  for (ip=0; ip<np; ip++) {
//    coord = part[ip];
//
//    x  =  X*cos_ttilt -  Y*sin_ttilt;
//    y  =  X*sin_ttilt +  Y*cos_ttilt;
//    xp = XP*cos_ttilt - YP*sin_ttilt;
//    yp = XP*sin_ttilt + YP*cos_ttilt;
//
//    X  = x;
//    Y  = y;
//    XP = xp;
//    YP = yp;
//
//    if (ctMode)
//      coord[4] /= c_mks;
//  }
//}

//TODO: add for particle output in CSBEND_CSR
/* this is used solely to undo the transformation done by 
 * convertFromCSBendCoords
 */
//void gpu_convertToCSBendCoords(long np, double rho0, double cos_ttilt, 
//                               double sin_ttilt, long ctMode)
//{
//  long ip;
//  double x, y, xp, yp, *coord;
//
//  for (ip=0; ip<np; ip++) {
//    coord = part[ip];
//
//    x  =   X*cos_ttilt +  Y*sin_ttilt;
//    y  =  -X*sin_ttilt +  Y*cos_ttilt;
//    xp =  XP*cos_ttilt + YP*sin_ttilt;
//    yp = -XP*sin_ttilt + YP*cos_ttilt;
//
//    X  = x;
//    Y  = y;
//    XP = xp;
//    YP = yp;
//
//    if (ctMode)
//      coord[4] *= c_mks;
//  }
//}

__device__ void 
gpu_addRadiationKick(double *Qx, double *Qy, double *dPoP, 
    double *sigmaDelta2, int sqrtOrder, double x, double h0, double Fx, 
    double Fy, double ds, double radCoef, double dsISR, double isrCoef, 
    int distributionBased, int includeOpeningAngle, 
    double meanPhotonsPerMeter, double normalizedCriticalEnergy0, 
    double Po, double *d_gauss_rn, curandState_t *state, 
    double srGaussianLimit)
{
  double f, xp, yp, F2, F, deltaFactor, dsFactor;
  double nMean, dDelta, thetaRms;
  int i, nEmitted;
  double y, logy;
  double normalizedCriticalEnergy;
  
  f = (1+x*h0)/EXSQRT(sqr(1+*dPoP)-sqr(*Qx)-sqr(*Qy), sqrtOrder);
  xp = *Qx*f;
  yp = *Qy*f;
  dsFactor = EXSQRT(sqr(1+x*h0)+sqr(xp)+sqr(yp), sqrtOrder);
  F2 = sqr(Fx)+sqr(Fy);

  if (!distributionBased) {
    deltaFactor = sqr(1 + *dPoP);
    *Qx /= (1 + *dPoP);
    *Qy /= (1 + *dPoP);
    if (radCoef)
      *dPoP -= radCoef*deltaFactor*F2*ds*dsFactor;
    if (isrCoef>0) {
      /* The minus sign is for consistency with the previous version. */
      *dPoP -=
        isrCoef*deltaFactor*pow(F2,0.75)*sqrt(dsISR*dsFactor)*(*d_gauss_rn);
    }
    if (sigmaDelta2)
      *sigmaDelta2 += sqr(isrCoef*deltaFactor)*pow(F2,1.5)*dsISR*dsFactor;
    *Qx *= (1 + *dPoP);
    *Qy *= (1 + *dPoP);
  } else {
    F = sqrt(F2);
    /* Compute the mean number of photons emitted = meanPhotonsPerMeter*meters */
    nMean = meanPhotonsPerMeter*dsISR*dsFactor*F;
    /* Pick the actual number of photons emitted from Poisson distribution */
    nEmitted = gpu_inversePoissonCDF(nMean, curand_uniform_double(state));
    /* Adjust normalized critical energy to local field strength 
       (FSE is already included via rho_actual) */
    /* Variation with energy offset included below. */
    normalizedCriticalEnergy = normalizedCriticalEnergy0*(1+*dPoP)*F;
    /* For each photon, pick its energy and emission angles */
    for (i=0; i<nEmitted; i++) {
      /* Pick photon energy normalized to critical energy */
      y=gpu_pickNormalizedPhotonEnergy(curand_uniform_double(state));
      /* Multiply by critical energy normalized to beam energy, adjusting for variation with
       * individual electron energy offset. */
      dDelta = normalizedCriticalEnergy*(1 + *dPoP)*y;
      // GPU: these reductions are for commented-out print diagnostics only
      //photonCount ++;
      //energyCount += y;
      /* Change the total electron momentum */
      *dPoP -= dDelta;
      if (includeOpeningAngle) {
        /* Compute rms spread in electron angle = (rms photon angle)*dDelta */
        logy = log10(y);
        thetaRms = dDelta*pow(10., -2.418673276661232e-01
                              + logy*(-4.472680955382907e-01+logy*(-4.535350424882360e-02
                              -logy*6.181818621278201e-03)))/Po;
        /* Compute change in electron angle due to photon angle */
        double dran = curand_normal_double(state);
        while (fabs(dran) > srGaussianLimit) dran = curand_normal_double(state);
        xp += thetaRms*dran;
        dran = curand_normal_double(state);
        while (fabs(dran) > srGaussianLimit) dran = curand_normal_double(state);
        yp += thetaRms*dran;
      }
    }
    f = (1 + *dPoP)/EXSQRT(sqr(1+x*h0)+sqr(xp)+sqr(yp), sqrtOrder);
    *Qx = xp*f;
    *Qy = yp*f;
  }
  
}

void gpu_initializeInterpolationTables() {
  /* load interpolation tables into constant memory */
  static int interpTablesInitialized=0;
  if (interpTablesInitialized) return;
  interpTablesInitialized=1;

  double ksiTable[200] = {
1.000000000000000e-07, 1.103351074554523e-07, 1.217383310646075e-07, 1.343200559586096e-07, 1.482020747927429e-07,
1.635189113578160e-07, 1.804187271717404e-07, 1.990651067932036e-07, 2.196385495378513e-07, 2.423384085535426e-07,
2.673842889192374e-07, 2.950186137528527e-07, 3.255088868939687e-07, 3.591505332405233e-07, 3.962690515521040e-07,
4.372236992560780e-07, 4.824109227537678e-07, 5.322685173365927e-07, 5.872789374272963e-07, 6.479745823257918e-07,
7.149429940622619e-07, 7.888329490825732e-07, 8.703595415651299e-07, 9.603117573208774e-07, 1.059560346180786e-06,
1.169066741218290e-06, 1.289890855361064e-06, 1.423201921186706e-06, 1.570290408531488e-06, 1.732581081736759e-06,
1.911644944229382e-06, 2.109214732888093e-06, 2.327202951554652e-06, 2.567720983189001e-06, 2.833097369770374e-06,
3.125899929822278e-06, 3.448963034857157e-06, 3.805415579829938e-06, 4.198708930670648e-06, 4.632648449870741e-06,
5.111434736502657e-06, 5.639704554113462e-06, 6.222573518736223e-06, 6.865680963315421e-06, 7.575252265320169e-06,
8.358158724774669e-06, 9.221982709850737e-06, 1.017508139438384e-05, 1.122668092062936e-05, 1.238696398672945e-05,
1.366716918178818e-05, 1.507968136669955e-05, 1.663817388115531e-05, 1.835773664119261e-05, 2.025502748788774e-05,
2.234840008736815e-05, 2.465811862080412e-05, 2.720654509486478e-05, 3.001836976776838e-05, 3.312079173855694e-05,
3.654384295607505e-05, 4.032066206311618e-05, 4.448784496925588e-05, 4.908569920644159e-05, 5.415873276357582e-05,
5.975605436563109e-05, 6.593190648243325e-05, 7.274602260436053e-05, 8.026436452046312e-05, 8.855971070475770e-05,
9.771244913354602e-05, 1.078111117026042e-04, 1.189534517802125e-04, 1.312473286188584e-04, 1.448118705606887e-04,
1.597782996280689e-04, 1.762914808166189e-04, 1.945112638414439e-04, 2.146141870658145e-04, 2.367947478384840e-04,
2.612676279365202e-04, 2.882697280504828e-04, 3.180626634504074e-04, 3.509347182037930e-04, 3.872040386295810e-04,
4.272217166086854e-04, 4.713754443547046e-04, 5.200925166529275e-04, 5.738444084509028e-04, 6.331514454110419e-04,
6.985881553542330e-04, 7.707878748140519e-04, 8.504493031356108e-04, 9.383435740575402e-04, 1.035322090321881e-03,
1.142323583631410e-03, 1.260383487096962e-03, 1.390644637378946e-03, 1.534368742954572e-03, 1.692947192375836e-03,
1.867914439234830e-03, 2.060964191854564e-03, 2.273966197872862e-03, 2.508982774009835e-03, 2.768287894108917e-03,
3.054391670052556e-03, 3.370064913211098e-03, 3.718364390303022e-03, 4.102660002531562e-03, 4.526671789275104e-03,
4.994505863796759e-03, 5.510692966986859e-03, 6.080227102658105e-03, 6.708621448554406e-03, 7.401960924834442e-03,
8.166960997457158e-03, 9.011022498794720e-03, 9.942316075907171e-03, 1.096985907697749e-02, 1.210360516825099e-02,
1.335452196639184e-02, 1.473471854468308e-02, 1.625755786686624e-02, 1.793779326078131e-02, 1.979167811637735e-02,
2.183715833862031e-02, 2.409403673367596e-02, 2.658418068941497e-02, 2.933167681381903e-02, 3.236312132092089e-02,
3.570786033002121e-02, 3.939830569107346e-02, 4.347015239952341e-02, 4.796281661196050e-02, 5.291978786613488e-02,
5.838910396032759e-02, 6.442366592223334e-02, 7.108188831787103e-02, 7.842822357081021e-02, 8.653385973120035e-02,
9.547720673026440e-02, 1.053448316706813e-01, 1.162322544203413e-01, 1.282449697899047e-01, 1.414991970199129e-01,
1.561232236887442e-01, 1.722586122799903e-01, 1.900616973883389e-01, 2.097047392306722e-01, 2.313778527830923e-01,
2.552908366685073e-01, 2.816753658394816e-01, 3.107867644741113e-01, 3.429067722728065e-01, 3.783463152734567e-01,
4.174487166108054e-01, 4.605924179038830e-01, 5.081949415377639e-01, 5.607170865377987e-01, 6.186676294134220e-01,
6.826074965088431e-01, 7.531554325044472e-01, 8.309943513916955e-01, 9.168782178988778e-01, 1.011638437307961e+00,
1.116191954411967e+00, 1.231550862322391e+00, 1.358832474428039e+00, 1.499269100004806e+00, 1.654219599620752e+00,
1.825183916868001e+00, 2.013817817791925e+00, 2.221947831729780e+00, 2.451587713539253e+00, 2.704960411634972e+00,
2.984519634347505e+00, 3.292972664221137e+00, 3.633303772560254e+00, 4.008807416353689e+00, 4.423119788364888e+00,
4.880253623874576e+00, 5.384631444881934e+00, 5.941135706944927e+00, 6.555154946882635e+00, 7.232636842499024e+00,
7.980135322277263e+00, 8.804886289535018e+00, 9.714875109915180e+00, 1.071891743371295e+01, 1.182672581369469e+01,
1.304902401296616e+01, 1.439764568247248e+01, 1.588565738231238e+01, 1.752745256838863e+01, 1.933892408641842e+01,
2.133760842747432e+01, 2.354287285156119e+01, 2.597604764912193e+01, 2.866068635656761e+01, 3.162277660168377e+01,
  };
  double FTable[200] = {
0.000000000000000e+00, 1.916076787477782e-04, 3.896006996482199e-04, 5.941918318862451e-04, 8.056009324383097e-04,
1.024055848381587e-03, 1.249790750550654e-03, 1.483048166648730e-03, 1.724078746036354e-03, 1.973142196708657e-03,
2.230505581886648e-03, 2.496445345396121e-03, 2.771247236692068e-03, 3.055207274452791e-03, 3.348630028390361e-03,
3.651830594751359e-03, 3.965134731031564e-03, 4.288879835176022e-03, 4.623413241840414e-03, 4.969094094184835e-03,
5.326293748409966e-03, 5.695396753910589e-03, 6.076799201662367e-03, 6.470910424341261e-03, 6.878153743490802e-03,
7.298967431515415e-03, 7.733803160081558e-03, 8.183127443537616e-03, 8.647422816544402e-03, 9.127188753348749e-03,
9.622940276393589e-03, 1.013520903749105e-02, 1.066454503150837e-02, 1.121151744173274e-02, 1.177671348365064e-02,
1.236073899369794e-02, 1.296422080554008e-02, 1.358780748228750e-02, 1.423216848805338e-02, 1.489799412460975e-02,
1.558599872144761e-02, 1.629692120583610e-02, 1.703152471578071e-02, 1.779059568966145e-02, 1.857494805017956e-02,
1.938542355142036e-02, 2.022289197174982e-02, 2.108824911944652e-02, 2.198242222447622e-02, 2.290636998738409e-02,
2.386108352008151e-02, 2.484758298036615e-02, 2.586692442491193e-02, 2.692019946744418e-02, 2.800853716992024e-02,
2.913309895920114e-02, 3.029508724043934e-02, 3.149574455091229e-02, 3.273635664445640e-02, 3.401824527268329e-02,
3.534277891084329e-02, 3.671137126806664e-02, 3.812548585471972e-02, 3.958662612641901e-02, 4.109634873515033e-02,
4.265626122542306e-02, 4.426802842804741e-02, 4.593335936539433e-02, 4.765402350963371e-02, 4.943184769110832e-02,
5.126872354805775e-02, 5.316659282366442e-02, 5.512746484065590e-02, 5.715341374017695e-02, 5.924658623150298e-02,
6.140918640941097e-02, 6.364349310252398e-02, 6.595185823803101e-02, 6.833671466591608e-02, 7.080056062143326e-02,
7.334597653614160e-02, 7.597562485484224e-02, 7.869225776043576e-02, 8.149870143512994e-02, 8.439787179667740e-02,
8.739277614110098e-02, 9.048652052440445e-02, 9.368229400251835e-02, 9.698338259579789e-02, 1.003931731738744e-01,
1.039151601553213e-01, 1.075529299247917e-01, 1.113101721273651e-01, 1.151906862423360e-01, 1.191983871317649e-01,
1.233372898347266e-01, 1.276115170379167e-01, 1.320253087906072e-01, 1.365830262538971e-01, 1.412891370418868e-01,
1.461482173909104e-01, 1.511649654280279e-01, 1.563442022251273e-01, 1.616908577721921e-01, 1.672099659551390e-01,
1.729066816736924e-01, 1.787862779608226e-01, 1.848541324986933e-01, 1.911157129897795e-01, 1.975765981791920e-01,
2.042424693223850e-01, 2.111190968366426e-01, 2.182123129823302e-01, 2.255280364093208e-01, 2.330722555939911e-01,
2.408510146807283e-01, 2.488703695479055e-01, 2.571364147677684e-01, 2.656552557248533e-01, 2.744329918358574e-01,
2.834756510338721e-01, 2.927892168723024e-01, 3.023795848158047e-01, 3.122525396990814e-01, 3.224136624350130e-01,
3.328683532505147e-01, 3.436217660758924e-01, 3.546787751835536e-01, 3.660438466342482e-01, 3.777210511874600e-01,
3.897139689097260e-01, 4.020256374046105e-01, 4.146583795759362e-01, 4.276137967068266e-01, 4.408926345108801e-01,
4.544946994027770e-01, 4.684186411623745e-01, 4.826619073829413e-01, 4.972205662043195e-01, 5.120891723276200e-01,
5.272605134477985e-01, 5.427255018050126e-01, 5.584729557729631e-01, 5.744894096219003e-01, 5.907588309168045e-01,
6.072624513895929e-01, 6.239785205673239e-01, 6.408820784452591e-01, 6.579446823895981e-01, 6.751342192611512e-01,
6.924146905220633e-01, 7.097460264751286e-01, 7.270839321656986e-01, 7.443798099419868e-01, 7.615807315005799e-01,
7.786294876113157e-01, 7.954647879517510e-01, 8.120215670459850e-01, 8.282314411609915e-01, 8.440233375507971e-01,
8.593244148619452e-01, 8.740611430503441e-01, 8.881606055288961e-01, 9.015520524341847e-01, 9.141687939214495e-01,
9.259502059677074e-01, 9.368438193470870e-01, 9.468075155760183e-01, 9.558117659236037e-01, 9.638415906288208e-01,
9.708980284014210e-01, 9.769991556010101e-01, 9.821804314564566e-01, 9.864941142754962e-01, 9.900075473722704e-01,
9.928006021230259e-01, 9.949622418923878e-01, 9.965864066620354e-01, 9.977674409043544e-01, 9.985957390441925e-01,
9.991538996011060e-01, 9.995138052403904e-01, 9.997348433078885e-01, 9.998634876235176e-01, 9.999340435105117e-01,
9.999702897374164e-01, 9.999876125809349e-01, 9.999952573365697e-01, 9.999983471293699e-01, 9.999994807411241e-01,
9.999998545219742e-01, 9.999999640793696e-01, 9.999999922833868e-01, 9.999999985785893e-01, 9.999999997790499e-01,
9.999999999715266e-01, 9.999999999970198e-01, 9.999999999997560e-01, 9.999999999999872e-01, 1.000000000000000e+00,
  };

  cudaMemcpyToSymbol(c_ksiTable, ksiTable, sizeof(double)*200,
                     0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(c_FTable, FTable, sizeof(double)*200, 
                     0, cudaMemcpyHostToDevice);
}

#include <gpu_interp.hcu> // gpu_interp

/* Return randomly-chosen photon energy normalized to the critical energy */
__device__ double gpu_pickNormalizedPhotonEnergy(double RN)
{
  int interpCode;
  double value;

  value = gpu_interp(c_ksiTable, c_FTable, 200, RN, 0, 2, &interpCode);
  if (!interpCode)
    return c_ksiTable[0];
  return value;
}
