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
extern "C"
{
#endif
  void gpu_exactDrift(long np, double length);
#ifdef __cplusplus
}
#endif

#ifdef sqr
#  undef sqr
#endif
#define sqr(x) x *x
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
__device__ __constant__ double d_coef[expansionOrderMax * expansionOrderMax];

// Defined below
__device__ int
gpu_integrate_kick_multipole_ordn(gpuParticleAccessor &coord, long *order, double *KnL, short *skew, long n_parts, long i_part, long iSlice,
                                      double drift, long integration_order, double *dzLoss, double *sigmaDelta2, int radial, double refTilt,
                                      double srGaussianLimit, short expandHamiltonian);
__device__ void gpu_apply_canonical_multipole_kicks(double *qx, double *qy,
                                                    double *sum_Fx_return, double *sum_Fy_return, double x, double y,
                                                    int order, double KnL, int skew);
__device__ void gpu_applyRadialCanonicalMultipoleKicks(double *qx, double *qy,
                                                       double *sum_Fx_return, double *sum_Fy_return, double x, double y,
                                                       int order, double KnL, int skew);
__device__ int gpu_convertSlopesToMomenta(double *qx, double *qy, double xp, double yp, double delta, short expandHamiltonian);
__device__ int gpu_convertMomentaToSlopes(double *xp, double *yp, double qx, double qy, double delta, short expandHamiltonian);

class gpu_multipole_tracking2_kernel
{
 public:
  unsigned int *d_sortIndex;
  double *d_sigmaDelta2;
  curandState_t *state;
  double dx, dy, xkick, ykick, Po, rad_coef, isr_coef;
  double KnL0, KnL1, KnL2;
  double drift, z_start;
  long order0, order1, order2;
  short skew0, skew1, skew2;
  int n_parts, i_part, integ_order;
  short expandHamiltonian;
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
  double refTilt;

 gpu_multipole_tracking2_kernel(unsigned int *d_sortIndex,
                                double *d_sigmaDelta2, curandState_t *state, double dx, double dy,
                                double xkick, double ykick, double Po, double rad_coef, double isr_coef,
                                double KnL0, double KnL1, double KnL2, double drift, double z_start, 
                                long order0, long order1, long order2, short skew0, short skew1, short skew2,
                                int n_parts, int i_part, int integ_order, int multDataOrders, int edgeMultDataOrders,
                                int steeringMultDataOrders, int present, double xMax, double xCen,
                                double yMax, double yCen, int radial, double refTilt, double srGaussianLimit, short expandHamiltonian) : d_sortIndex(d_sortIndex), d_sigmaDelta2(d_sigmaDelta2),
    state(state), dx(dx), dy(dy), xkick(xkick), ykick(ykick), Po(Po),
    rad_coef(rad_coef), isr_coef(isr_coef), KnL0(KnL0), KnL1(KnL1), KnL2(KnL2), drift(drift),
    z_start(z_start), order0(order0), order1(order1), order2(order2), 
    skew0(skew0), skew1(skew1), skew2(skew2), n_parts(n_parts), i_part(i_part),
    integ_order(integ_order), multDataOrders(multDataOrders),
    edgeMultDataOrders(edgeMultDataOrders),
    steeringMultDataOrders(steeringMultDataOrders), present(present),
    xMax(xMax), xCen(xCen), yMax(yMax), yCen(yCen), radial(radial), refTilt(refTilt),
    srGaussianLimit(srGaussianLimit), expandHamiltonian(expandHamiltonian){};

  __device__ unsigned int operator()(gpuParticleAccessor &coord)
  {
    unsigned int tid = coord.getParticleIndex();

    int particle_lost;
    double dzLoss;
    double *tSigmaDelta2 = NULL;
    double KnL[3];
    long order[3];
    short skew[3];
    KnL[0] = KnL0;
    KnL[1] = KnL1;
    KnL[2] = KnL2;
    order[0] = order0;
    order[1] = order1;
    order[2] = order2;
    skew[0] = skew0;
    skew[1] = skew1;
    skew[2] = skew2;
    if (d_sigmaDelta2)
      tSigmaDelta2 = &d_sigmaDelta2[tid];

    particle_lost = !gpu_integrate_kick_multipole_ordn(coord, order, KnL, skew, n_parts, i_part,
                                                       drift, integ_order, &dzLoss, tSigmaDelta2, radial, refTilt, srGaussianLimit, expandHamiltonian);
    if (particle_lost)
      {
        coord[4] = z_start + dzLoss;
        coord[5] = Po * (1 + coord[5]);
        d_sortIndex[tid] = tid + coord.getParticlePitch();
        return 0;
      }
    d_sortIndex[tid] = tid;
    return 1;
  }

#define x coord[0]
#define xp coord[1]
#define y coord[2]
#define yp coord[3]
#define dp coord[5]


  /* BETA is 2^(1/3) */
#define BETA 1.25992104989487316477

  __device__ int
    gpu_integrate_kick_multipole_ordn(gpuParticleAccessor &coord, long *order, double *KnL, short *skew, long n_parts, long i_part,
                                      double drift, long integration_order, double *dzLoss, double *sigmaDelta2, int radial, double refTilt,
                                      double srGaussianLimit, short expandHamiltonian)
  {
    double p, qx, qy, beta0, beta1, s;
    double delta_qx, delta_qy;
    int i_kick, step, imult;
    double dsh;
    
    static double driftFrac2[2] = {
      0.5, 0.5
    };
    static double kickFrac2[2] = {
      1.0, 0.0
    };
    static double driftFrac4[4] = {
      0.5/(2-BETA),  (1-BETA)/(2-BETA)/2,  (1-BETA)/(2-BETA)/2,  0.5/(2-BETA)
    };
    static double kickFrac4[4] = {
      1./(2-BETA),  -BETA/(2-BETA),  1/(2-BETA),  0
    };

    /* From AOP-TN-2020-064 */
    static double driftFrac6[8] = {
      0.39225680523878, 0.5100434119184585, -0.47105338540975655, 0.0687531682525181,
      0.0687531682525181, -0.47105338540975655, 0.5100434119184585, 0.39225680523878,
    } ;
    static double kickFrac6[8] = {
      0.784513610477560, 0.235573213359357, -1.17767998417887, 1.3151863206839063,
      -1.17767998417887,  0.235573213359357, 0.784513610477560, 0
    } ;

    double *driftFrac = NULL, *kickFrac = NULL;
    long nSubsteps = 0;

    switch (integration_order) {
    case 2:
      nSubsteps = 2;
      driftFrac = driftFrac2;
      kickFrac = kickFrac2;
      break;
    case 4:
      nSubsteps = 4;
      driftFrac = driftFrac4;
      kickFrac = kickFrac4;
      break;
    case 6:
      nSubsteps = 8;
      driftFrac = driftFrac6;
      kickFrac = kickFrac6;
      break;
    default:
      printf("invalid order %ld given for symplectic integrator", integration_order);
      return(0);
    }


    drift = drift / n_parts;
    //KnL = KnL / n_parts;
    //KnL2 = KnL2 / n_parts;
    xkick = xkick / n_parts;
    ykick = ykick / n_parts;

    //x = coord[0];
    //xp = coord[1];
    //y = coord[2];
    //yp = coord[3];
    s = 0;
    //dp = coord[5];
    p = Po * (1 + dp);
    beta0 = p / sqrt(sqr(p) + 1);

#if defined(IEEE_MATH)
    if (isnan(x) || isnan(xp) || isnan(y) || isnan(yp))
      {
        return 0;
      }
#endif
    if (fabs(x) > COORD_LIMIT || fabs(y) > COORD_LIMIT ||
        fabs(xp) > SLOPE_LIMIT || fabs(yp) > SLOPE_LIMIT)
      {
        return 0;
      }

    /* calculate initial canonical momenta */
    gpu_convertSlopesToMomenta(&qx, &qy, xp, yp, dp, expandHamiltonian);

    if (i_part<=0) {
    if (edgeMultDataOrders >= 0)
      {
        for (imult = 0; imult < edgeMultDataOrders; imult++)
          {
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
    }
    /* We must do this in case steering or edge multipoles were run. We do it even if not in order
     * to avoid numerical precision issues that may subtly change the results
     */
    if (!gpu_convertMomentaToSlopes(&xp, &yp, qx, qy, dp, expandHamiltonian))
      return 0;

    *dzLoss = 0;
    for (i_kick = 0; i_kick < n_parts; i_kick++)
      {
        if (present &&
            ((xMax && fabs(x + dx - xCen) > xMax) ||
             (yMax && fabs(y + dy - yCen) > yMax)))
          {
            return 0;
          }
        /* obstructions are not implemented in GPU code
        if (insideObstruction_xyz(x, xp, y, yp, coord[particleIDIndex], 
			                            globalLossCoordOffset>0?coord+globalLossCoordOffset:NULL, 
			                            refTilt,  GLOBAL_LOCAL_MODE_SEG, 0.0, i_kick, n_parts))
         {
            return 0;
          }
        */
        delta_qx = delta_qy = 0;
        for (step = 0; step < nSubsteps; step++)
          {
            if (drift)
              {
                dsh = drift * driftFrac[step];
                x += xp * dsh;
                y += yp * dsh;
                if (expandHamiltonian)
                  {
                    s += dsh * (1 + (sqr(xp) + sqr(yp)) / 2);
                  }
                else
                  {
                    s += dsh * sqrt(1 + sqr(xp) + sqr(yp));
                  }
                *dzLoss += dsh;
              }

            if (!kickFrac[step])
              break;

            delta_qx = delta_qy = 0;

            if (!radial) {
              for (int iOrder=0; iOrder<3; iOrder++) {
                if (KnL[iOrder]) {
                  gpu_apply_canonical_multipole_kicks(&qx, &qy, &delta_qx, &delta_qy, x, y, 
                                                      order[iOrder], KnL[iOrder]/n_parts*kickFrac[step], skew[iOrder]);
                }
              }
            } else 
              gpu_applyRadialCanonicalMultipoleKicks(&qx, &qy, &delta_qx, &delta_qy, x, y, 
                                                     order[0], KnL[0]/n_parts*kickFrac[step], 0);

            if (xkick)
              gpu_apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, x, y, 0, -xkick * kickFrac[step], 0);
            if (ykick)
              gpu_apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, x, y, 0, -ykick * kickFrac[step], 1);

            /* apply steering corrector multipoles */
            if (steeringMultDataOrders >= 0)
              {
                for (imult = 0; imult < steeringMultDataOrders; imult++)
                  {
                    if (steeringMultDataKnL[imult])
                      gpu_apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, x, y,
                                                          steeringMultDataOrder[imult],
                                                          steeringMultDataKnL[imult] * xkick * kickFrac[step], 0);
                    if (steeringMultDataJnL[imult])
                      gpu_apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, x, y,
                                                          steeringMultDataOrder[imult],
                                                          steeringMultDataJnL[imult] * ykick * kickFrac[step], 1);
                  }
              }

            /* do kicks for spurious multipoles */
            for (imult = 0; imult < multDataOrders; imult++)
              {
                if (multDataKnL[imult])
                  {
                    gpu_apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, x, y,
                                                        multDataOrder[imult],
                                                        multDataKnL[imult] * kickFrac[step] / n_parts,
                                                        0);
                  }
                if (multDataJnL[imult])
                  {
                    gpu_apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, x, y,
                                                        multDataOrder[imult],
                                                        multDataJnL[imult] * kickFrac[step] / n_parts,
                                                        1);
                  }
              }
            if (!gpu_convertMomentaToSlopes(&xp, &yp, qx, qy, dp, expandHamiltonian))
              return 0;
            if ((rad_coef || isr_coef) && drift)
              {
                double deltaFactor, F2, dsFactor, dsISRFactor;
                qx /= (1 + dp);
                qy /= (1 + dp);
                deltaFactor = sqr(1 + dp);
                /* delta_qx and delta_qy are for the last step and have kickFrac[step-1] included, so remove it */
                delta_qx /= kickFrac[step];
                delta_qy /= kickFrac[step];
                F2 = sqr(delta_qx/drift-xkick/drift)+sqr(delta_qy/drift+ykick/drift);
                delta_qx = 0;
                delta_qy = 0;
                dsFactor = xp * xp + yp * yp;
                dsFactor = sqrt(1 + dsFactor);
                dsISRFactor = dsFactor * drift / (nSubsteps - 1); /* recall that kickFrac may be negative */
                dsFactor *= drift * kickFrac[step]; /* that's ok here, since we don't take sqrt */
                if (rad_coef)
                  dp -= rad_coef * deltaFactor * F2 * dsFactor;
                if (isr_coef > 0)
                  {
                    // Get a srGaussianLimit sigma (sigma=1) limited random number
                    double rand = curand_uniform_double(state);
                    while (fabs(rand) > srGaussianLimit)
                      rand = curand_uniform_double(state);
                    dp -= isr_coef * deltaFactor * pow(F2, 0.75) * sqrt(dsISRFactor) * rand;
                  }
                if (sigmaDelta2)
                  *sigmaDelta2 += sqr(isr_coef * deltaFactor) * pow(F2, 1.5) * dsISRFactor;
                qx *= (1 + dp);
                qy *= (1 + dp);
                if (!gpu_convertMomentaToSlopes(&xp, &yp, qx, qy, dp, expandHamiltonian))
                  return 0;
              }
          }
        if (i_part>=0)
          break;
      }
    if (present &&
        ((xMax && fabs(x + dx - xCen) > xMax) ||
         (yMax && fabs(y + dy - yCen) > yMax)))
      {
        return 0;
      }
    /* obstructions are not implemented in GPU code
    if (insideObstruction_xy(x, xp, y, yp, coord[particleIDIndex], 
			                       globalLossCoordOffset>0?coord+globalLossCoordOffset:NULL, 
			                       refTilt,  GLOBAL_LOCAL_MODE_SEG, 0.0, i_kick, n_parts))
      {
        return 0;
      }
    */
  if (i_part<0 || i_part==(n_parts-1)) {
    if (edgeMultDataOrders >= 0)
      {
        for (imult = 0; imult < edgeMultDataOrders; imult++)
          {
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
  }
    if (!gpu_convertMomentaToSlopes(&xp, &yp, qx, qy, dp, expandHamiltonian))
      return 0;

    //coord[0] = x;
    //coord[1] = xp;
    //coord[2] = y;
    //coord[3] = yp;
    if (rad_coef)
      {
        p = Po * (1 + dp);
        beta1 = p / sqrt(sqr(p) + 1);
        coord[4] = beta1 * (coord[4] / beta0 + 2 * s / (beta0 + beta1));
      }
    else
      coord[4] += s;
    //coord[5] = dp;

#if defined(IEEE_MATH)
    if (isnan(x) || isnan(xp) || isnan(y) || isnan(yp))
      {
        return 0;
      }
#endif
    if (fabs(x) > COORD_LIMIT || fabs(y) > COORD_LIMIT ||
        fabs(xp) > SLOPE_LIMIT || fabs(yp) > SLOPE_LIMIT)
      {
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

__device__ extern unsigned long multipoleKicksDone;
unsigned short host_expandHamiltonian;

extern "C"
{

  long gpu_multipole_tracking2(long n_part, ELEMENT_LIST *elem,
                               double p_error, double Po, double **accepted, double z_start,
                               MAXAMP *maxamp, APCONTOUR *apcontour, APERTURE_DATA *apFileData,
                               double *sigmaDelta2, long iSlice)
  {
    double *KnL;  /* integrated strength = L/(B.rho)*(Dx^n(By))_o for central momentum */
    long *order;  /* order (n) */
    short *skew;


    /*double KnL(0), KnL2(0);     
      long order, order2;    */     
    double dx(0), dy(0), dz(0); /* offsets of the multipole center */
    long n_kicks(0), integ_order(0), iOrder;
    long i_top, nSlices;
    double *coef;
    double drift(0);
    double tilt, pitch, yaw, rad_coef, isr_coef, xkick, ykick;
    KQUAD *kquad = NULL;
    KSEXT *ksext = NULL;
    KQUSE *kquse = NULL;
    KOCT *koct = NULL;
    double lEffective = -1, lEnd = 0;
    short doEndDrift = 0, malignMethod;
    unsigned long d_multipoleKicksDone;

    MULTIPOLE_DATA *multData = NULL, *steeringMultData = NULL, *edgeMultData = NULL;
    long freeMultData = 0;
    MULT_APERTURE_DATA apertureData;
    
    KnL = (double*)(calloc(sizeof(double), sizeof(double) * 3));
    order = (long*)(calloc(sizeof(long), sizeof(long) * 3));
    skew = (short*)(calloc(sizeof(short), sizeof(short) * 3));

    log_entry("multipole_tracking2");

    if (!elem)
      bombTracking("null element pointer (multipole_tracking2)");
    if (!elem->p_elem)
      bombTracking("null p_elem pointer (multipole_tracking2)");

    rad_coef = xkick = ykick = isr_coef = 0;
    pitch = yaw = tilt = 0;
    dx = dy = dz = 0;
    malignMethod = 0;
    /*order2 = 0;
      KnL2 = 0;*/

    switch (elem->type)
      {
      case T_KQUAD:
        kquad = ((KQUAD *)elem->p_elem);
        nSlices = kquad->nSlices;
        n_kicks = kquad->n_kicks;
        host_expandHamiltonian = kquad->expandHamiltonian;
        order[0] = 1;
        if ((lEffective = kquad->lEffective) <= 0)
          lEffective = kquad->length;
        else
          {
            lEnd = (kquad->length - lEffective) / 2;
            doEndDrift = 1;
          }
        if (kquad->bore)
          /* KnL = d^nB/dx^n * L/(B.rho) = n! B(a)/a^n * L/(B.rho) * (1+FSE) */
          KnL[0] = kquad->B / kquad->bore * (particleCharge / (particleMass * c_mks * Po)) * lEffective * (1 + kquad->fse);
        else
          KnL[0] = kquad->k1 * lEffective * (1 + kquad->fse);
        drift = lEffective;
        tilt = kquad->tilt;
        pitch = kquad->pitch;
        yaw = kquad->yaw;
        dx = kquad->dx;
        dy = kquad->dy;
        dz = kquad->dz;
        malignMethod = kquad->malignMethod;
        xkick = kquad->xkick * kquad->xKickCalibration;
        ykick = kquad->ykick * kquad->yKickCalibration;
        integ_order = kquad->integration_order;
        if (kquad->synch_rad)
          rad_coef = sqr(particleCharge) * pow(Po, 3.) / (6 * PI * epsilon_o * sqr(c_mks) * particleMass);
        isr_coef = particleRadius * sqrt(55.0 / (24 * sqrt(3)) * pow(Po, 5.) * 137.0359895);
        if (!kquad->isr || (kquad->isr1Particle == 0 && n_part == 1))
          /* Minus sign indicates we accumulate into sigmaDelta^2 only, don't perturb particles */
          isr_coef *= -1;
        if (kquad->length < 1e-6 && (kquad->isr || kquad->synch_rad))
          {
            rad_coef = isr_coef = 0; /* avoid unphysical results */
            printWarningForTracking((char*)"Quadrupole with length < 1e-6 has SYNCH_RAD=0 and ISR=0 forced to avoid unphysical results.",
                                    NULL);
          }
        if (!kquad->multipolesInitialized)
          {
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
        if (!kquad->totalMultipolesComputed)
          {
            computeTotalErrorMultipoleFields(&(kquad->totalMultipoleData),
                                             &(kquad->systematicMultipoleData),
                                             kquad->systematicMultipoleFactor,
                                             &(kquad->edgeMultipoleData),
                                             NULL,
                                             &(kquad->randomMultipoleData),
                                             kquad->randomMultipoleFactor,
                                             &(kquad->steeringMultipoleData),
                                             kquad->steeringMultipoleFactor,
                                             KnL[0], 1, 1,
                                             kquad->minMultipoleOrder, kquad->maxMultipoleOrder);
            kquad->totalMultipolesComputed = 1;
          }
        multData = &(kquad->totalMultipoleData);
        edgeMultData = &(kquad->edgeMultipoleData);
        steeringMultData = &(kquad->steeringMultipoleData);
        break;
      case T_KSEXT:
        ksext = ((KSEXT *)elem->p_elem);
        n_kicks = ksext->n_kicks;
        nSlices = ksext->nSlices;
        host_expandHamiltonian = ksext->expandHamiltonian;
        order[0] = 2;
        if (ksext->bore)
          /* KnL = d^nB/dx^n * L/(B.rho) = n! B(a)/a^n * L/(B.rho) * (1+FSE) */
          KnL[0] = 2 * ksext->B / sqr(ksext->bore) * (particleCharge / (particleMass * c_mks * Po)) * ksext->length * (1 + ksext->fse);
        else
          KnL[0] = ksext->k2 * ksext->length * (1 + ksext->fse);
        drift = ksext->length;
        tilt = ksext->tilt;
        pitch = ksext->pitch;
        yaw = ksext->yaw;
        dx = ksext->dx;
        dy = ksext->dy;
        dz = ksext->dz;
        malignMethod = ksext->malignMethod;
        xkick = ksext->xkick * ksext->xKickCalibration;
        ykick = ksext->ykick * ksext->yKickCalibration;
        integ_order = ksext->integration_order;
        if (ksext->k1)
          {
            KnL[1] = ksext->k1 * ksext->length;
            order[1] = 1;
          }
        if (ksext->j1)
          {
            KnL[2] = ksext->j1 * ksext->length;
            order[2] = 1;
            skew[2] = 1;
          }
        if (ksext->synch_rad)
          rad_coef = sqr(particleCharge) * pow(Po, 3.) / (6 * PI * epsilon_o * sqr(c_mks) * particleMass);
        isr_coef = particleRadius * sqrt(55.0 / (24 * sqrt(3)) * pow(Po, 5.) * 137.0359895);
        if (!ksext->isr || (ksext->isr1Particle == 0 && n_part == 1))
          /* Minus sign indicates we accumulate into sigmaDelta^2 only, don't perturb particles */
          isr_coef *= -1;
        if (ksext->length < 1e-6 && (ksext->isr || ksext->synch_rad))
          {
            rad_coef = isr_coef = 0; /* avoid unphysical results */
            printWarningForTracking((char*)"Sextupole with length < 1e-6 has SYNCH_RAD=0 and ISR=0 forced to avoid unphysical results.",
                                    NULL);
          }
        if (!ksext->multipolesInitialized)
          {
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
        if (!ksext->totalMultipolesComputed)
          {
            computeTotalErrorMultipoleFields(&(ksext->totalMultipoleData),
                                             &(ksext->systematicMultipoleData),
                                             ksext->systematicMultipoleFactor,
                                             &(ksext->edgeMultipoleData),
                                             NULL,
                                             &(ksext->randomMultipoleData),
                                             ksext->randomMultipoleFactor,
                                             &(ksext->steeringMultipoleData),
                                             ksext->steeringMultipoleFactor,
                                             KnL[0], 2, 1,
                                       ksext->minMultipoleOrder, ksext->maxMultipoleOrder);
            ksext->totalMultipolesComputed = 1;
          }
        multData = &(ksext->totalMultipoleData);
        edgeMultData = &(ksext->edgeMultipoleData);
        steeringMultData = &(ksext->steeringMultipoleData);
        break;
      case T_KOCT:
        koct = ((KOCT *)elem->p_elem);
        n_kicks = koct->n_kicks;
        nSlices = koct->nSlices;
        host_expandHamiltonian = koct->expandHamiltonian;
        order[0] = 3;
        if (koct->bore)
          /* KnL = d^nB/dx^n * L/(B.rho) = n! B(a)/a^n * L/(B.rho) * (1+FSE) */
          KnL[0] = 6 * koct->B / pow(koct->bore, 3.) * (particleCharge / (particleMass * c_mks * Po)) * koct->length * (1 + koct->fse);
        else
          KnL[0] = koct->k3 * koct->length * (1 + koct->fse);
        drift = koct->length;
        tilt = koct->tilt;
        pitch = koct->pitch;
        yaw = koct->yaw;
        dx = koct->dx;
        dy = koct->dy;
        dz = koct->dz;
        malignMethod = koct->malignMethod;
        integ_order = koct->integration_order;
        if (koct->synch_rad)
          rad_coef = sqr(particleCharge) * pow(Po, 3.) / (6 * PI * epsilon_o * sqr(c_mks) * particleMass);
        isr_coef = particleRadius * sqrt(55.0 / (24 * sqrt(3)) * pow(Po, 5.) * 137.0359895);
        if (!koct->isr || (koct->isr1Particle == 0 && n_part == 1))
          /* Minus sign indicates we accumulate into sigmaDelta^2 only, don't perturb particles */
          isr_coef *= -1;
        if (koct->length < 1e-6 && (koct->isr || koct->synch_rad))
          {
            rad_coef = isr_coef = 0; /* avoid unphysical results */
            printWarningForTracking((char*)"Octupole with length < 1e-6 has SYNCH_RAD=0 and ISR=0 forced to avoid unphysical results.", 
                                    NULL);
          }
        if (!koct->multipolesInitialized)
          {
            /* read the data files for the error multipoles */
            readErrorMultipoleData(&(koct->systematicMultipoleData),
                                   koct->systematic_multipoles, 0);
            readErrorMultipoleData(&(koct->randomMultipoleData),
                                   koct->random_multipoles, 0);
            koct->multipolesInitialized = 1;
          }
        if (!koct->totalMultipolesComputed)
          {
            computeTotalErrorMultipoleFields(&(koct->totalMultipoleData),
                                             &(koct->systematicMultipoleData),
                                             1.0,
                                             NULL,
                                             NULL,
                                             &(koct->randomMultipoleData),
                                             1.0,
                                             NULL,
                                             1.0,
                                             KnL[0], 3, 1,
                                             NULL, NULL);
            koct->totalMultipolesComputed = 1;
          }
        multData = &(koct->totalMultipoleData);
        break;
      case T_KQUSE:
        /* Implemented as a quadrupole with sextupole as a secondary multipole */
        kquse = ((KQUSE *)elem->p_elem);
        n_kicks = kquse->n_kicks;
        nSlices = kquse->nSlices;
        host_expandHamiltonian = kquse->expandHamiltonian;
        order[0] = 1;
        KnL[0] = kquse->k1 * kquse->length * (1 + kquse->fse1);
        drift = kquse->length;
        tilt = kquse->tilt;
        dx = kquse->dx;
        dy = kquse->dy;
        dz = kquse->dz;
        malignMethod = 0;
        integ_order = kquse->integration_order;
        if (kquse->synch_rad)
          rad_coef = sqr(particleCharge) * pow(Po, 3.) / (6 * PI * epsilon_o * sqr(c_mks) * particleMass);
        isr_coef = particleRadius * sqrt(55.0 / (24 * sqrt(3)) * pow(Po, 5.) * 137.0359895);
        if (!kquse->isr || (kquse->isr1Particle == 0 && n_part == 1))
          /* Minus sign indicates we accumulate into sigmaDelta^2 only, don't perturb particles */
          isr_coef *= -1;
        if (kquse->length < 1e-6 && (kquse->isr || kquse->synch_rad))
          {
            rad_coef = isr_coef = 0; /* avoid unphysical results */
            printWarningForTracking((char*)"KQUSE with length < 1e-6 has SYNCH_RAD=0 and ISR=0 forced to avoid unphysical results.", NULL);
          }
        KnL[1] = kquse->k2 * kquse->length * (1 + kquse->fse2);
        order[1] = 2;
        break;
      default:
        printf("error: multipole_tracking2() called for element %s--not supported!\n", elem->name);
        fflush(stdout);
        KnL[0] = dx = dy = dz = tilt = drift = 0;
        integ_order = order[0] = n_kicks = nSlices = 0;
        exitElegant(1);
        break;
      }
    if (multData && !multData->initialized)
      multData = NULL;

    if (order[0] <= 0)
      bombTracking("order <= 0 in multipole()");
    if (integ_order != 2 && integ_order != 4 && integ_order != 6)
      bombTracking("multipole integration_order must be 2, 4 or 6");

    for (iOrder=0; iOrder<3; iOrder++)
      {
        if (KnL[iOrder] && !expansion_coefficients(order[iOrder]))
          bombTracking("expansion_coefficients() returned null pointer (multipole_tracking2)");
      }

    i_top = n_part - 1;

    if (n_kicks<=0) {
      if (nSlices<=0)
        bombTracking("N_KICKS<=0 and N_SLICES<=0 in multipole tracking");
    } else {
      if (integ_order>2) {
        if ((nSlices = ceil(n_kicks/(1.0*integ_order)))<1)
          nSlices = 1;
        n_kicks = nSlices*integ_order;
      }
      else
        nSlices = n_kicks;
    }

    //Copy multipoleKicksDone from host cpu to gpu device
    cudaMemcpyToSymbol(&d_multipoleKicksDone, &multipoleKicksDone, sizeof(long), 0, cudaMemcpyHostToDevice);

    d_multipoleKicksDone += (i_top + 1) * n_kicks;
    if (multData)
      d_multipoleKicksDone += (i_top + 1) * n_kicks * multData->orders;
    //Copy d_multipoleKicksDone from gpu device to host cpu
    cudaMemcpyToSymbol(&multipoleKicksDone, &d_multipoleKicksDone, sizeof(long), 0, cudaMemcpyDeviceToHost);

    setupMultApertureData(&apertureData, -tilt, apcontour, maxamp, apFileData, NULL, z_start + drift / 2, elem);

    struct GPUBASE *gpuBase = getGpuBase();
    unsigned int particlePitch = gpuBase->gpu_array_pitch;
    unsigned int *d_sortIndex = gpuBase->d_tempu_alpha;

    if (iSlice<=0) {
      if (malignMethod!=0) {
        if (dx || dy || dz || tilt || pitch || yaw) {
          bombTracking("gpu_multipole_tracking2: Non zero MALIGN_METHOD not supported in the GPU version.");
          /*
            if (malignMethod==1) {
            gpu_offsetParticlesForEntranceCenteredMisalignmentExact
            (n_part, dx, dy, dz, pitch, yaw, tilt, 0.0, 0.0, drift, 1);
            }
            else {
            gpu_offsetParticlesForBodyCenteredMisalignmentExact
            (n_part, dx, dy, dz, pitch, yaw, tilt, 0.0, 0.0, drift, 1);
            }
          */
        }
      }
    else {
      if (dx || dy || dz)
        gpu_offsetBeamCoordinatesForMisalignment(n_part, dx, dy, dz);
      if (tilt)
        gpu_rotateBeamCoordinatesForMisalignment(n_part, tilt);
    }

    if (doEndDrift)
      {
        gpu_exactDrift(n_part, lEnd);
      }

    /* Fringe treatment, if any */
    switch (elem->type)
      {
      case T_KQUAD:
        if (kquad->edge1_effects > 0)
          //TODO quadFringe
          std::cout << "Implement gpu_quadFringe" << std::endl;
        //quadFringe(particle, n_part, kquad->k1, kquad->fringeIntM, kquad->fringeIntP, -1, kquad->edge1_effects, kquad->edge1Linear, kquad->edge1NonlinearFactor);
        break;
      default:
        break;
      }
    }

    // Copy multData and steeringMultData to device
    if (multData)
      {
        cudaMemcpyToSymbol(multDataOrder, multData->order,
                           sizeof(multDataOrder) * multData->orders, 0, cudaMemcpyHostToDevice);
        cudaMemcpyToSymbol(multDataKnL, multData->KnL,
                           sizeof(multDataKnL) * multData->orders, 0, cudaMemcpyHostToDevice);
        cudaMemcpyToSymbol(multDataJnL, multData->JnL,
                           sizeof(multDataJnL) * multData->orders, 0, cudaMemcpyHostToDevice);
      }
    if (steeringMultData)
      {
        cudaMemcpyToSymbol(steeringMultDataOrder, steeringMultData->order,
                           sizeof(steeringMultDataOrder) * steeringMultData->orders, 0,
                           cudaMemcpyHostToDevice);
        cudaMemcpyToSymbol(steeringMultDataKnL, steeringMultData->KnL,
                           sizeof(steeringMultDataKnL) * steeringMultData->orders, 0,
                           cudaMemcpyHostToDevice);
        cudaMemcpyToSymbol(steeringMultDataJnL, steeringMultData->JnL,
                           sizeof(steeringMultDataJnL) * steeringMultData->orders, 0,
                           cudaMemcpyHostToDevice);
      }
    if (edgeMultData)
      {
        cudaMemcpyToSymbol(edgeMultDataOrder, edgeMultData->order,
                           sizeof(edgeMultDataOrder) * edgeMultData->orders, 0, cudaMemcpyHostToDevice);
        cudaMemcpyToSymbol(edgeMultDataKnL, edgeMultData->KnL,
                           sizeof(edgeMultDataKnL) * edgeMultData->orders, 0, cudaMemcpyHostToDevice);
        cudaMemcpyToSymbol(edgeMultDataJnL, edgeMultData->JnL,
                           sizeof(edgeMultDataJnL) * edgeMultData->orders, 0, cudaMemcpyHostToDevice);
      }

    // Ensure expansionOrderMax is sufficiently large
    int maxorder = order[0];
    for (int iorder = 0; iorder < expansionOrderMax - 1; ++iorder)
      {
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
    if (maxorder > expansionOrderMax - 1)
      bombTracking("gpu_multipole_tracking2: Error: increase expansionOrderMax");

    // Copy expansion coefficients to constant memory
    for (int iorder = 0; iorder <= maxorder; ++iorder)
      {
        coef = expansion_coefficients(iorder);
        cudaMemcpyToSymbol(d_coef, coef, sizeof(double) * (iorder + 1),
                           sizeof(double) * (iorder * expansionOrderMax),
                           cudaMemcpyHostToDevice);
        gpuErrorHandler("gpu_multipole_tracking2::cudaMemcpyToSymbol");
      }

    double *d_gauss_rn = gpuBase->d_temp_particles + particlePitch;
    curandState_t *state = NULL;
    if (isr_coef > 0)
      state =
        (curandState_t *)gpu_get_rand_state(d_gauss_rn, n_part, random_2(0));

    double *d_sigmaDelta2 = NULL;
    if (sigmaDelta2)
      {
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
                                                          dx, dy, xkick, ykick, Po, rad_coef, isr_coef, KnL[0], KnL[1], KnL[2], drift, z_start, 
                                                          order[0], order[1], order[2], skew[0], skew[1], skew[2],
                                                          nSlices, iSlice, integ_order, multData ? multData->orders : -1,
                                                          edgeMultData ? edgeMultData->orders : -1,
                                                          steeringMultData ? steeringMultData->orders : -1, apertureData.present,
                                                          apertureData.xMax, apertureData.xCen, apertureData.yMax,
                                                          apertureData.yCen, elem->type == T_KQUAD ? kquad->radial : 0, tilt, srGaussianLimit, host_expandHamiltonian));
    gpuErrorHandler("gpu_multipole_tracking2::gpu_multipole_tracking2_kernel");

    if (sigmaDelta2)
      {
        *sigmaDelta2 = gpuReduceAdd(d_sigmaDelta2, n_part);
        *sigmaDelta2 /= n_part;
      }

  if (iSlice<0 || iSlice==(nSlices-1)) {
    /* Fringe treatment, if any */
    switch (elem->type)
      {
      case T_KQUAD:
        if (kquad->edge2_effects > 0)
          bombTracking("gpu_multipole_tracking2: quadFringe not implemented");
        //TODO quadFringe
        //quadFringe(particle, n_part, kquad->k1, kquad->fringeIntM, kquad->fringeIntP, 1, kquad->edge2_effects, kquad->edge2Linear, kquad->edge2NonlinearFactor);
        break;
      default:
        break;
      }
    if (doEndDrift)
      {
        gpu_exactDrift(n_part, lEnd);
      }

    if (malignMethod!=0) {
      if (dx || dy || dz || tilt || pitch || yaw)  {
        bombTracking("gpu_multipole_tracking2: Non zero MALIGN_METHOD not supported in the GPU version.");
        /*
          if (malignMethod==1) {
          gpu_offsetParticlesForEntranceCenteredMisalignmentExact
          (n_part, dx, dy, dz, pitch, yaw, tilt, 0.0, 0.0, drift, 2);
          }
          else {
          gpu_offsetParticlesForBodyCenteredMisalignmentExact
          (n_part, dx, dy, dz, pitch, yaw, tilt, 0.0, 0.0, drift, 2);
          }
        */
      }
    } else {
      if (tilt)
        gpu_rotateBeamCoordinatesForMisalignment(n_part, -tilt);
      if (dx || dy || dz)
        gpu_offsetBeamCoordinatesForMisalignment(n_part, -dx, -dy, -dz);
    }
  }

    if (freeMultData && !multData->copy)
      {
        if (multData->order)
          free(multData->order);
        if (multData->KnL)
          free(multData->KnL);
        free(multData);
      }
    free(KnL);
    free(order);
    free(skew);

    log_exit("multipole_tracking2");
    host_expandHamiltonian = 0;
    return n_part;
  }

} // extern "C"

__device__ void gpu_apply_canonical_multipole_kicks(double *qx, double *qy,
                                                    double *delta_qx_return, double *delta_qy_return, double x, double y,
                                                    int order, double KnL, int skew)
{
  int i;
  double sum_Fx, sum_Fy;

  /* sum up the terms for the multipole expansion */
  for (i = sum_Fx = sum_Fy = 0; i <= order; i++)
    {
      if (ODD(i))
        sum_Fx += d_coef[expansionOrderMax * order + i] * pow(x, 1. * order - i) * pow(y, 1. * i);
      else
        sum_Fy += d_coef[expansionOrderMax * order + i] * pow(x, 1. * order - i) * pow(y, 1. * i);
    }
  if (skew)
    {
      SWAP_DOUBLE(sum_Fx, sum_Fy);
      sum_Fx = -sum_Fx;
    }
  /* add the kicks */
  *qx -= KnL * sum_Fy;
  *qy += KnL * sum_Fx;
  if (delta_qx_return)
    *delta_qx_return -= KnL*sum_Fy;
  if (delta_qy_return)
    *delta_qy_return += KnL*sum_Fx;
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
  if (x == 0)
    {
      if (y == 0)
        return;
      xypow = pow(y, 1. * order);
      i = order;
      ratio = 0;
    }
  else
    {
      xypow = pow(x, 1. * order);
      ratio = y / x;
      i = 0;
    }
  /* now sum up the terms for the multipole expansion */
  for (sum_Fx = sum_Fy = 0; i <= order; i++)
    {
      if (ODD(i))
        sum_Fx -= d_coef[expansionOrderMax * order + i] * xypow;
      else
        sum_Fy += d_coef[expansionOrderMax * order + i] * xypow;
      xypow *= ratio;
    }
  if (skew)
    {
      SWAP_DOUBLE(sum_Fx, sum_Fy);
      sum_Fx = -sum_Fx;
    }
  /* add the kicks */
  *qx -= KnL * sum_Fy;
  *qy += KnL * sum_Fx;
  if (sum_Fx_return)
    *sum_Fx_return = sum_Fx;
  if (sum_Fy_return)
    *sum_Fy_return = sum_Fy;
}

__device__ int gpu_convertSlopesToMomenta(double *qx, double *qy, double xp, double yp, double delta, short expandHamiltonian)
{
  if (expandHamiltonian)
    {
      *qx = (1 + delta) * xp;
      *qy = (1 + delta) * yp;
    }
  else
    {
      double denom;
      denom = sqrt(1 + sqr(xp) + sqr(yp));
      *qx = (1 + delta) * xp / denom;
      *qy = (1 + delta) * yp / denom;
    }
  return 1;
}

__device__ int gpu_convertMomentaToSlopes(double *xp, double *yp, double qx, double qy, double delta, short expandHamiltonian)
{
  static short warningCounter = 100;
  if (expandHamiltonian)
    {
      *xp = qx / (1 + delta);
      *yp = qy / (1 + delta);
    }
  else
    {
      double denom;
      if ((denom = (1 + delta) * (1 + delta) - (qx * qx + qy * qy)) <= 0)
        {
          if (warningCounter) 
            {
              printf("Warning: particle acquired undefined slopes when integrating through kick multipole\n");
              if (--warningCounter==0)
                printf("         No further warnings of this type will be issued.\n");
              /*fflush(stdout);*/
            }
          return 0;
        }
      denom = sqrt(denom);
      *xp = qx / denom;
      *yp = qy / denom;
    }
  return 1;
}

