#include <float.h> // DBL_MAX

#include <gpu_track.h>

#include <gpu_particle_template_function.hcu>
#include <gpu_csbend.h> // gpu_exactDrift
#include <gpu_matter.hcu> // gpu_track_through_matter_dfunc, gpu_setMatterTables
#include <gpu_matter.h> // gpu_set_track_through_matter
#undef align
#include <gpu_killParticles.hcu>

#ifdef SIGN
#undef SIGN
#endif
#define SIGN(x) ((x)<0?-1:((x)>0?1:0))
#ifdef sqr
#undef sqr
#endif
#define sqr(x) (x*x)

__device__ inline unsigned int evaluateLostWithOpenSidesDevice(unsigned int code, 
    double dx, double dy, double xsize, double ysize)
{
  long lost = 1;
  switch (code) {
  case OPEN_PLUS_X:
    if (dx>0 && fabs(dy)<ysize)
      lost = 0;
    break;
  case OPEN_MINUS_X:
    if (dx<0 && fabs(dy)<ysize)
      lost = 0;
    break;
  case OPEN_PLUS_Y:
    if (dy>0 && fabs(dx)<xsize)
      lost = 0;
    break;
  case OPEN_MINUS_Y:
    if (dy<0 && fabs(dx)<xsize)
      lost = 0;
    break;
  default:
    break;
  }
  return lost;
}

class gpu_rectangular_collimator_kernel1{
public:
  unsigned int * d_sortIndex;
  double xsize, ysize, x_center, y_center, z, Po;
  unsigned int invert;
  long openCode;
  gpu_rectangular_collimator_kernel1(unsigned int* d_sortIndex, double xsize, 
      double ysize, double x_center, double y_center, unsigned int invert, 
      double z, double Po, long openCode): 
    d_sortIndex(d_sortIndex), xsize(xsize), ysize(ysize), x_center(x_center), 
    y_center(y_center), invert(invert), z(z), Po(Po), openCode(openCode) {};

  __device__ unsigned int operator()(gpuParticleAccessor& part){
    unsigned int tid = part.getParticleIndex();
    double dx, dy;
    unsigned int lost;
    dx = part[0] - x_center;
    dy = part[2] - y_center;
    lost = 0;
    if ((xsize>0 && fabs(dx) > xsize) ||
        (ysize>0 && fabs(dy) > ysize)) {
      lost = openCode ? evaluateLostWithOpenSidesDevice(openCode, 
                                             dx, dy, xsize, ysize) : 1;
    } else if (isinf(part[0]) || isinf(part[2]) ||
               isnan(part[0]) || isnan(part[2]) )
      lost = 1;
    if (invert)
      lost = !lost;
    if (lost) {
      d_sortIndex[tid] = tid+part.getParticlePitch();
      part[4] = z; /* record position of particle loss */
      part[5] = Po*(1+part[5]);
      return 0;
    } else {
      d_sortIndex[tid] = tid;
      return 1;
    }
  }
};

class gpu_rectangular_collimator_kernel2{
public:
  unsigned int * d_sortIndex;
  double xsize, ysize, x_center, y_center, length, z, Po;
  int invert;
  long openCode;
  gpu_rectangular_collimator_kernel2(unsigned int* d_sortIndex, double xsize, 
      double ysize, double x_center, double y_center, double length, 
      double z, double Po, long openCode, int invert): 
    d_sortIndex(d_sortIndex), xsize(xsize), ysize(ysize), x_center(x_center), 
    y_center(y_center), length(length), z(z), Po(Po), openCode(openCode),
    invert(invert) {};

  __device__ unsigned int operator()(gpuParticleAccessor& part){
    unsigned int tid = part.getParticleIndex();
    double x1, y1, dx, dy, zx, zy;
    unsigned int is_out;
    
    x1 = part[0] + length*part[1];
    y1 = part[2] + length*part[3];
    dx = x1 - x_center;
    dy = y1 - y_center;
    is_out = 0;
    if (xsize>0 && fabs(dx)>xsize)
      is_out += 1*(openCode?evaluateLostWithOpenSidesDevice(openCode, 
            dx, 0, xsize, ysize):1);
    if (ysize>0 && fabs(dy)>ysize)
      is_out += 2*(openCode?evaluateLostWithOpenSidesDevice(openCode, 
            0, dy, xsize, ysize):1);
    if (isinf(x1) || isinf(y1) || isnan(x1) || isnan(y1) )
      is_out += 4;
    if (is_out&4) {
      part[4] = z+length;
      part[0] = x1;
      part[2] = y1;
      part[5] = Po*(1+part[5]);
      d_sortIndex[tid] = tid+part.getParticlePitch();
      return 0;
    } else if ((is_out && !invert) || (!is_out && invert)) {
      if (!openCode) {
        zx = zy = DBL_MAX;
        if (is_out&1 && part[1]!=0)
          zx = (SIGN(part[1])*xsize-(part[0]-x_center))/part[1];
        if (is_out&2 && part[3]!=0)
          zy = (SIGN(part[3])*ysize-(part[2]-y_center))/part[3];
        if (zx<zy) {
          part[0] += part[1]*zx;
          part[2] += part[3]*zx;
          part[4] = z+zx;
        }
        else {
          part[0] += part[1]*zy;
          part[2] += part[3]*zy;
          part[4] = z+zy;
        }
      }
      part[5] = Po*(1+part[5]);
      d_sortIndex[tid] = tid+part.getParticlePitch();
      return 0;
    } else {
      part[4] += length*sqrt(1+sqr(part[1])+sqr(part[3]));
      part[0] = x1;
      part[2] = y1;
      d_sortIndex[tid] = tid;
      return 1;
    }
  }
};

extern "C" {

/* routine: gpu_rectangular_collimator()
 * purpose: eliminate particles that hit the walls of a non-zero length
 *          rectangular collimator
 */
long gpu_rectangular_collimator(RCOL *rcol, long np, double **accepted, double z,
                                double Po) {
  long openCode;
  double xsize, ysize, x_center, y_center, length;

  xsize  = rcol->x_max;
  ysize  = rcol->y_max;
  x_center = rcol->dx;
  y_center = rcol->dy;
 
  /*
  if (rcol->invert && rcol->length) {
    TRACKING_CONTEXT tc;
    getTrackingContext(&tc);
    fprintf(stderr, "Problem for %s#%ld:\n", tc.elementName, tc.elementOccurrence);
    bombElegant("Cannot have invert=1 and non-zero length for RCOL", NULL);
  }
  */

  if (xsize<=0 && ysize<=0) {
    gpu_exactDrift(np, rcol->length);
    return(np);
  }
  openCode = determineOpenSideCode(rcol->openSide);
    
  unsigned int* d_sortIndex = getGpuBase()->d_tempu_alpha;
  np = killParticles(np, d_sortIndex, accepted,
    gpu_rectangular_collimator_kernel1(d_sortIndex, xsize, ysize, 
    x_center, y_center, rcol->invert, z, Po, openCode));
  gpuErrorHandler("gpu_rectangular_collimator::gpu_rectangular_collimator_kernel1");
  if (np==0 || (length=rcol->length)<=0) {
    return(np);
  }

  return killParticles(np, d_sortIndex, accepted,
    gpu_rectangular_collimator_kernel2(d_sortIndex, xsize, ysize, 
    x_center, y_center, length, z, Po, openCode, rcol->invert));
}

} // extern "C"

class gpu_limit_amplitudes_kernel{
public:
  unsigned int * d_sortIndex;
  double xmax, ymax, z, Po;
  long extrapolate_z, openCode;
  gpu_limit_amplitudes_kernel(unsigned int* d_sortIndex, double xmax, 
      double ymax, double z, double Po, long extrapolate_z, long openCode): 
    d_sortIndex(d_sortIndex), xmax(xmax), ymax(ymax), z(z), Po(Po),
    extrapolate_z(extrapolate_z), openCode(openCode) {};

  __device__ unsigned int operator()(gpuParticleAccessor& part){
    unsigned int tid = part.getParticleIndex();
    double dz, dzx, dzy;
    unsigned int is_out=0;

    if (xmax>0 && fabs(part[0])>xmax)
      is_out += 1;
    if (ymax>0 && fabs(part[2])>ymax)
      is_out += 2;
    if (openCode)
      is_out *= evaluateLostWithOpenSidesDevice(openCode, part[0], 
                                                part[2], xmax, ymax);
    dz = 0;
    if (is_out && !(is_out&4) && !openCode && extrapolate_z) {
      /* find the actual position of loss, assuming a drift preceded with 
       * the same aperture 
       */
      dzx = dzy = -DBL_MAX;
      if (is_out&1 && part[1]!=0)
        dzx = (part[0]-SIGN(part[1])*xmax)/part[1];
      if (is_out&2 && part[3]!=0)
        dzy = (part[2]-SIGN(part[3])*ymax)/part[3];
      if (dzx>dzy)
        dz = -dzx;
      else
        dz = -dzy;
      if (dz==-DBL_MAX)
        dz = 0;
      part[0] += dz*part[1];
      part[2] += dz*part[3];
    }
    if (is_out) {
      d_sortIndex[tid] = tid + part.getParticlePitch();
      part[4] = z+dz;  /* record position of loss */
      part[5] = Po*(1+part[5]);
      return 0;
    } else {
      d_sortIndex[tid] = tid;
      return 1;
    }
  }
};

extern "C" {

/* routine: gpu_limit_amplitudes()
 * purpose: eliminate particles with (x,y) larger than given values
 */
long gpu_limit_amplitudes(
    double xmax, double ymax, long np, double **accepted,
    double z, double Po, long extrapolate_z, long openCode) {

  if (xmax<0 && ymax<0) {
    return(np);
  }

  unsigned int* d_sortIndex = getGpuBase()->d_tempu_alpha;
  return killParticles(np, d_sortIndex, accepted,
    gpu_limit_amplitudes_kernel(d_sortIndex, xmax, ymax, z, Po,
    extrapolate_z, openCode));
}

} // extern "C" 

class gpu_removeInvalidParticles_kernel{
public:
  unsigned int * d_sortIndex;
  double xmax, ymax, z, Po;
  long extrapolate_z, openCode;
  gpu_removeInvalidParticles_kernel(unsigned int* d_sortIndex, 
                                    double z, double Po): 
    d_sortIndex(d_sortIndex), z(z), Po(Po) {};

  __device__ unsigned int operator()(gpuParticleAccessor& part){
    unsigned int tid = part.getParticleIndex();
    unsigned int is_out, ic;

    is_out = 0;
    for (ic=0; ic<6; ic++)
      if (isnan(part[ic])) {
        is_out = 1;
        break;
      }
    if (part[5]<=-1)
      is_out = 1;
    if (is_out) {
      d_sortIndex[tid] = tid + part.getParticlePitch();
      part[4] = z;
      part[5] = Po*(1+part[5]);
      return 0;
    } else {
      d_sortIndex[tid] = tid;
      return 1;
    }
  }
};

extern "C" {

long gpu_removeInvalidParticles(long np, double **accepted, double z, double Po)
{
  unsigned int* d_sortIndex = getGpuBase()->d_tempu_alpha;
  return killParticles(np, d_sortIndex, accepted,
    gpu_removeInvalidParticles_kernel(d_sortIndex, z, Po));
}

} // extern "C" 

class gpu_elliptical_collimator_kernel1{
public:
  unsigned int * d_sortIndex;
  double xsize, ysize, dx, dy, xe, ye, a2, b2, z, Po;
  long openCode;
  gpu_elliptical_collimator_kernel1(unsigned int* d_sortIndex, double z, 
      double Po, double xsize, double ysize, double dx, double dy, long xe,
      long ye, double a2, double b2, long openCode): 
    d_sortIndex(d_sortIndex), z(z), Po(Po), xsize(xsize), ysize(ysize),
    dx(dx), dy(dy), xe((double)xe), ye((double)ye), a2(a2), b2(b2), 
    openCode(openCode) {};

  __device__ unsigned int operator()(gpuParticleAccessor& part){
    unsigned int tid = part.getParticleIndex();
    unsigned int lost;
    double xo, yo;

    lost = 0;
    xo = part[0] - dx;
    yo = part[2] - dy;
    if ((pow(xo, xe)/a2 + pow(yo, ye)/b2)>1)
      lost = openCode ? 
        evaluateLostWithOpenSidesDevice(openCode, xo, yo, xsize, ysize) : 1;
    else if (isinf(part[0]) || isinf(part[2]) ||
             isnan(part[0]) || isnan(part[2]) )
      lost = 1;
    if (lost) {
      d_sortIndex[tid] = tid + part.getParticlePitch();
      part[4] = z;
      part[5] = sqrt(sqr(Po*(1+part[5]))+1);
      return 0;
    } else {
      d_sortIndex[tid] = tid;
      return 1;
    }
  }
};

class gpu_elliptical_collimator_kernel2{
public:
  unsigned int * d_sortIndex;
  double z, Po, xsize, ysize, length, exponent, dx, dy;
  long openCode;
  int invert;
  gpu_elliptical_collimator_kernel2(unsigned int* d_sortIndex, double z, 
      double Po, double xsize, double ysize, double length, double dx,
      double dy, long exponent, long openCode, int invert): 
    d_sortIndex(d_sortIndex), z(z), Po(Po), xsize(xsize), ysize(ysize),
    length(length), dx(dx), dy(dy), exponent((double)exponent),
    openCode(openCode), invert(invert) {};

  __device__ unsigned int operator()(gpuParticleAccessor& part){
    unsigned int tid = part.getParticleIndex();
    unsigned int lost;
    double xo, yo;

    lost = 0;
    part[0] += length*part[1];
    part[2] += length*part[3];
    xo = (part[0]-dx)/xsize;
    yo = (part[2]-dy)/ysize;
    if ((pow(xo, exponent) + pow(yo, exponent))>1)
      lost = openCode ? evaluateLostWithOpenSidesDevice(openCode, xo, yo, 1, 1) : 1;
    else if (isinf(part[0]) || isinf(part[2]) ||
             isnan(part[0]) || isnan(part[2]) )
      lost = 1;
    if (invert) 
      lost = !lost;
    if (lost) {
      d_sortIndex[tid] = tid + part.getParticlePitch();
      part[4] = z + length;
      part[5] = sqrt(sqr(Po*(1+part[5]))+1);
      return 0;
    } else {
      part[4] += length*sqrt(1+sqr(part[1])+sqr(part[3]));
      d_sortIndex[tid] = tid;
      return 1;
    }
  }
};

extern "C" {

/* routine: gpu_elliptical_collimator()
 * purpose: eliminate particles that hit the walls of a non-zero length
 *          elliptical collimator
 */
long gpu_elliptical_collimator(ECOL *ecol, long np, double **accepted, 
                               double z, double Po) {
  double length;
  long openCode;
  double a2, b2;
  double dx, dy, xsize, ysize;
  TRACKING_CONTEXT context;
  long xe, ye;
  unsigned int* d_sortIndex = getGpuBase()->d_tempu_alpha;

  xsize = ecol->x_max;
  ysize = ecol->y_max;
  if ((xe=ecol->exponent)<2 || xe%2) {
    getTrackingContext(&context);
    fprintf(stderr, "Error for %s: exponent=%ld is not valid.  Give even integer >=2\n",
            context.elementName, xe);
    exitElegant(1);
  }
  ye = xe;
  if (ecol->yExponent) {
    if ((ye=ecol->yExponent)<2 || ye%2) {
      getTrackingContext(&context);
      fprintf(stderr, "Error for %s: exponent=%ld is not valid.  Give even integer >=2\n",
              context.elementName, ye);
      exitElegant(1);
    }
  }

  a2 = std::pow(ecol->x_max, xe);
  b2 = std::pow(ecol->y_max, ye);
  dx = ecol->dx;
  dy = ecol->dy;

  if (ecol->x_max<=0 || ecol->y_max<=0) {
    /* At least one of x_max or y_max is non-positive */
    if (ecol->x_max>0 || ecol->y_max>0) {
      /* One of x_max or y_max is positive. Use rectangular collimator routine to implement this. */
      RCOL rcol;
      rcol.length = ecol->length;
      rcol.x_max = ecol->x_max;
      rcol.y_max = ecol->y_max;
      rcol.dx = ecol->dx;
      rcol.dy = ecol->dy;
      rcol.invert = ecol->invert;
      rcol.openSide = ecol->openSide;
      return gpu_rectangular_collimator(&rcol, np, accepted, z, Po);
    }
    gpu_exactDrift(np, ecol->length);
    return(np);
  }
  openCode = determineOpenSideCode(ecol->openSide);

  np = killParticles(np, d_sortIndex, accepted,
    gpu_elliptical_collimator_kernel1(d_sortIndex, z, Po,
    xsize, ysize, dx, dy, xe, ye, a2, b2, openCode));

  if (np==0 || (length=ecol->length)<=0)
    return(np);

  return killParticles(np, d_sortIndex, accepted,
    gpu_elliptical_collimator_kernel2(d_sortIndex, z, Po,
    xsize, ysize, length, dx,dy, ecol->exponent, openCode, ecol->invert));
}

} // extern "C"

class gpu_elimit_amplitudes_kernel{
public:
  unsigned int * d_sortIndex;
  double xmax, ymax, dx, dy, xe, ye, a2, b2, z, Po;
  long extrapolate_z, openCode;
  gpu_elimit_amplitudes_kernel(unsigned int* d_sortIndex, double z, 
      double Po, double xmax, double ymax, long xe, long ye, 
      double a2, double b2, long extrapolate_z, long openCode): 
    d_sortIndex(d_sortIndex), z(z), Po(Po), xmax(xmax), ymax(ymax),
    xe((double)xe), ye((double)ye), a2(a2), b2(b2),
    extrapolate_z(extrapolate_z), openCode(openCode) {};

  __device__ unsigned int operator()(gpuParticleAccessor& part){
    unsigned int tid = part.getParticleIndex();
    unsigned int lost;
    double c1, c2, c0, dz, det;
 
    if (isinf(part[0]) || isinf(part[2]) || isnan(part[0]) || isnan(part[2]) ) {
      d_sortIndex[tid] = tid + part.getParticlePitch();
      part[4] = z;
      part[5] = sqrt(sqr(Po*(1+part[5]))+1);
      return 0;
    }
    lost = 0;
    if ((pow(part[0], xe)/a2 + pow(part[2], ye)/b2) > 1) {
      lost = 1;
    } else {
      d_sortIndex[tid] = tid;
      return 1;
    }
    if (openCode)
      lost *= evaluateLostWithOpenSidesDevice(openCode, part[0], part[2], xmax, ymax);
    if (lost) {
      dz = 0;
      if (extrapolate_z && !openCode && fabs(xe-2)<1e-12 && fabs(ye-2)<1e-12) {
        c0 = sqr(part[0])/a2 + sqr(part[2])/b2 - 1;
        c1 = 2*(part[0]*part[1]/a2 + part[2]*part[3]/b2);
        c2 = sqr(part[1])/a2 + sqr(part[3])/b2;
        det = sqr(c1)-4*c0*c2;
        if (z>0 && c2 && det>=0) {
          if ((dz = (-c1+sqrt(det))/(2*c2))>0)
            dz = (-c1-sqrt(det))/(2*c2);
          if ((z+dz)<0)
            dz = -z;
          part[0] += dz*part[1];
          part[2] += dz*part[3];
        }
      }
      d_sortIndex[tid] = tid + part.getParticlePitch();
      part[4] = z+dz;
      part[5] = sqrt(sqr(Po*(1+part[5]))+1);
      return 0;
    } else {
      d_sortIndex[tid] = tid;
      return 1;
    }
  }
};

extern "C" {

/* routine: gpu_elimit_amplitudes()
 * purpose: eliminate particles outside an ellipse with given semi-major
 *          and semi-minor axes.
 */
long gpu_elimit_amplitudes(double xmax, double ymax, long np, 
    double **accepted, double z, double Po, long extrapolate_z, 
    long openCode, long exponent, long yexponent) { 
  double a2, b2;
  TRACKING_CONTEXT context;
  long xe, ye;
  unsigned int* d_sortIndex = getGpuBase()->d_tempu_alpha;

  if ((xe=exponent)<2 || xe%2) {
    getTrackingContext(&context);
    fprintf(stderr, "Error for %s: exponent=%ld is not valid.  Give even integer >=2\n",
            context.elementName, xe);
    exitElegant(1);
  }
  ye = xe;
  if (yexponent) {
    if ((ye=yexponent)<2 || ye%2) {
      getTrackingContext(&context);
      fprintf(stderr, "Error for %s: exponent=%ld is not valid.  Give even integer >=2\n",
              context.elementName, ye);
      exitElegant(1);
    }
  }

  if (xmax<=0 || ymax<=0) {
    /* At least one of the dimensions is non-positive and therefore ignored */
    if (xmax>0 || ymax>0)
      return gpu_limit_amplitudes(xmax, ymax, np, accepted, 
                                  z, Po, extrapolate_z, openCode);
    return(np);
  }

  a2 = std::pow(xmax, xe);
  b2 = std::pow(ymax, ye);

  np = killParticles(np, d_sortIndex, accepted,
    gpu_elimit_amplitudes_kernel(d_sortIndex, z, Po,
    xmax, ymax, xe, ye, a2, b2, extrapolate_z, openCode));

  log_exit("gpu_elimit_amplitudes");
  return(np);
}

} // extern "C"

class gpu_beam_scraper_kernel1{
public:
  double dx, dy, limit, length;
  double Nrad, L, theta_rms, prob, L1, K2, probBS, probER, Po, dGammaFactor, z0;
  int do_x, do_y, sections0, multipleScattering;
  unsigned int *d_sortIndex;
  curandState_t *d_state;
  /* from MATTER struct */
  int energyDecay, nuclearBremsstrahlung, energyStraggle;
  int nSlots;
  double spacing,  width,  tilt,  center;

  gpu_beam_scraper_kernel1(int do_x, int do_y, double dx, double dy,
    double limit, double length, int energyDecay, int nuclearBremsstrahlung,
    int energyStraggle, int nSlots, double spacing, double width, double tilt,
    double center, int multipleScattering, double Nrad, double L,
    double theta_rms, int sections0, double prob, double L1, double K2,
    double probBS, double probER, double Po, double dGammaFactor, double z0,
    unsigned int *d_sortIndex, curandState_t *d_state) :
    do_x(do_x), do_y(do_y), dx(dx), dy(dy), limit(limit), length(length),
    energyDecay(energyDecay), nuclearBremsstrahlung(nuclearBremsstrahlung),
    energyStraggle(energyStraggle), nSlots(nSlots), spacing(spacing), width(width),
    tilt(tilt), center(center), multipleScattering(multipleScattering),
    Nrad(Nrad), L(L), theta_rms(theta_rms), sections0(sections0),
    prob(prob), L1(L1), K2(K2), probBS(probBS), probER(probER), Po(Po),
    dGammaFactor(dGammaFactor), z0(z0), d_sortIndex(d_sortIndex),
    d_state(d_state) {};

  __device__ unsigned int operator()(gpuParticleAccessor& part){
    unsigned int tid = part.getParticleIndex();

    if ((do_x && do_x*(part[0]-dx)>limit) ||
        (do_y && do_y*(part[2]-dy)>limit)) {
      /* scatter and/or absorb energy */
      //if (!track_through_matter(&part, 1, &matter, Po, NULL, z))
      if (!gpu_track_through_matter_dfunc(part, energyDecay,
             nuclearBremsstrahlung, energyStraggle, nSlots, spacing,
             width, tilt, center, multipleScattering,
             Nrad, L, theta_rms, sections0, prob, L1, K2, probBS, probER,
             Po, dGammaFactor, z0, d_sortIndex, d_state)) {
        part[5] = -1;
        return 0;
      }
    } else {
      part[0] = part[0] + part[1]*length;
      part[2] = part[2] + part[3]*length;
      part[4] += length*sqrt(1+sqr(part[1])+sqr(part[3]));
    }
    return 1;
  }
};

class gpu_beam_scraper_kernel2{
public:
  int do_x, do_y;
  double dx, dy, limit, Po, z;
  unsigned int *d_sortIndex;
  gpu_beam_scraper_kernel2(int do_x, int do_y, double dx, double dy,
    double limit, double Po, double z, unsigned int *d_sortIndex) :
    do_x(do_x), do_y(do_y), dx(dx), dy(dy), limit(limit), Po(Po), z(z),
    d_sortIndex(d_sortIndex) {};

  __device__ unsigned int operator()(gpuParticleAccessor& part){
    unsigned int tid = part.getParticleIndex();

    if ((do_x && do_x*(part[0]-dx) > limit) ||
        (do_y && do_y*(part[2]-dy) > limit) || part[5]<=-1) {
      part[4] = z; /* record position of particle loss */
      part[5] = Po*(1+part[5]);
      d_sortIndex[tid] = tid + part.getParticlePitch();
      return 0;
    }
    d_sortIndex[tid] = tid;
    return 1;
  }
};

class gpu_beam_scraper_kernel3{
public:
  int do_x, do_y;
  double dx, dy, limit, length, Po, z;
  unsigned int *d_sortIndex;
  gpu_beam_scraper_kernel3(int do_x, int do_y, double dx, double dy,
    double limit, double Po, double z, unsigned int *d_sortIndex) : 
    do_x(do_x), do_y(do_y), dx(dx), dy(dy), limit(limit), length(length), Po(Po), 
    z(z), d_sortIndex(d_sortIndex) {};

  __device__ unsigned int operator()(gpuParticleAccessor& part){
    unsigned int tid = part.getParticleIndex();

    part[0] += length*part[1];
    part[2] += length*part[3];
    if ((do_x && do_x*(part[0]-dx) > limit) ||
        (do_y && do_y*(part[2]-dy) > limit) ) {
      part[4] = z; /* record position of particle loss */
      part[5] = Po*(1+part[5]);
      d_sortIndex[tid] = tid + part.getParticlePitch();
      return 0;
    }
    else
      part[4] += length*sqrt(1+sqr(part[1])+sqr(part[3]));
    return 1;
  }
};

extern "C" {

long gpu_beam_scraper(SCRAPER *scraper, long np, double **accepted,
                      double z, double Po) {
  double length;
  long do_x, do_y;
  double limit;
  struct GPUBASE* gpuBase = getGpuBase();
  unsigned int particlePitch = gpuBase->gpu_array_pitch;
  unsigned int* d_sortIndex = gpuBase->d_tempu_alpha;
 
  log_entry("beam_scraper");

  if (scraper->direction<0 || scraper->direction>3)
    return np;
  
  if (scraper->direction==0 || scraper->direction==2) {
    do_x = scraper->direction==0 ? 1 : -1;
    do_y = 0;
    limit = scraper->position*do_x;
  }
  else {
    do_x = 0;
    do_y = scraper->direction==1 ? 1 : -1;
    limit = scraper->position*do_y;
  }    

  if (scraper->length && (scraper->Xo || scraper->Z)) {
    /* scraper has material properties that scatter beam and
     * absorb energy
     */
    MATTER matter;
    matter.length = scraper->length;
    matter.lEffective = 0;
    matter.Xo = scraper->Xo;
    matter.energyDecay = scraper->energyDecay;
    matter.energyStraggle = scraper->energyStraggle;
    matter.nuclearBremsstrahlung = scraper->nuclearBremsstrahlung;
    matter.electronRecoil = scraper->electronRecoil;
    matter.Z = scraper->Z;
    matter.A = scraper->A;
    matter.rho = scraper->rho;
    matter.pLimit = scraper->pLimit;
    matter.width = matter.spacing = matter.tilt = matter.center = 0;
    matter.nSlots = 0;
    
    /* setup matter variables */
    double L, Nrad, theta_rms=0;
    double dGammaFactor;
    double K2=0.0;
    int sections0=1;
    double L1, prob, probBS, probER;
    int multipleScattering = 0, impulseMode;

    if (matter.length!=0) {
      L = matter.length;
      impulseMode = 0;
    } else if (matter.lEffective!=0) {
      L = matter.lEffective;
      impulseMode = 1;
    }
    else
      return np;

    gpu_set_track_through_matter(np, &matter, Po, impulseMode, &multipleScattering, &Nrad,
           &L, &theta_rms, &sections0, &prob, &L1, &K2, &probBS, &probER, &dGammaFactor);
    gpu_setMatterTables();

    double* d_gauss_rn =  gpuBase->d_temp_particles + particlePitch;
    curandState_t* d_state = (curandState_t*)gpu_get_rand_state(d_gauss_rn, np, random_2(0));
    np = killParticles(np, d_sortIndex, accepted,
           gpu_beam_scraper_kernel1(do_x, do_y, scraper->dx, scraper->dy, limit,
             scraper->length, matter.energyDecay, matter.nuclearBremsstrahlung,
             matter.energyStraggle, matter.nSlots, matter.spacing, matter.width,
             matter.tilt, matter.center, multipleScattering, Nrad, L, theta_rms, sections0,
             prob, L1, K2, probBS, probER, Po, dGammaFactor, z, d_sortIndex, d_state));
    gpuErrorHandler("gpu_beam_scraper::gpu_beam_scraper_kernel1");

    log_exit("beam_scraper");
    return(np);
  }

  /* come here for idealized scraper that just absorbs particles */
  np = killParticles(np, d_sortIndex, accepted,
         gpu_beam_scraper_kernel2(do_x, do_y, scraper->dx, scraper->dy, limit,
           Po, z, d_sortIndex));
  gpuErrorHandler("gpu_beam_scraper::gpu_beam_scraper_kernel2");

  if (np==0 || (length=scraper->length)<=0) {
    log_exit("beam_scraper");
    return(np);
  }

  z += length;

  np = killParticles(np, d_sortIndex, accepted,
         gpu_beam_scraper_kernel3(do_x, do_y, scraper->dx, scraper->dy, limit,
           Po, z, d_sortIndex));
  gpuErrorHandler("gpu_beam_scraper::gpu_beam_scraper_kernel3");

  log_exit("beam_scraper");
  return(np);
}

} // extern "C"

class gpu_imposeApertureData_kernel{
public:
  double xCenter, yCenter, xSize, ySize, z, Po;
  unsigned int *d_sortIndex;
  gpu_imposeApertureData_kernel(double xCenter, double yCenter, double xSize,
    double ySize, double z, double Po, unsigned int *d_sortIndex) :
    xCenter(xCenter), yCenter(yCenter), xSize(xSize), ySize(ySize),
    z(z), Po(Po), d_sortIndex(d_sortIndex) {};

  __inline__ __device__ unsigned int operator()(gpuParticleAccessor& ini){
    unsigned int tid = ini.getParticleIndex();
    double dx, dy;

    dx = ini[0] - xCenter;
    dy = ini[2] - yCenter;
    if ((xSize && fabs(dx) > xSize) ||
        (ySize && fabs(dy) > ySize)) {
      ini[4] = z; /* record position of particle loss */
      ini[5] = Po*(1+ini[5]);
      d_sortIndex[tid] = tid + ini.getParticlePitch();
      return 0;
    }
    d_sortIndex[tid] = tid;
    return 1;
  }
};

extern "C" {

long interpolateApertureData(double z, APERTURE_DATA *apData,
                             double *xCenter, double *yCenter, double *xSize, double *ySize);

long gpu_imposeApertureData(long np, double **accepted,
                            double z, double Po, APERTURE_DATA *apData) {
  double xSize, ySize;
  double xCenter, yCenter;

#if DEBUG_APERTURE
  static FILE *fp = NULL;
  if (!fp) {
    TRACKING_CONTEXT tcontext;
    char s[1000];
    getTrackingContext(&tcontext);
    sprintf(s, "%s.aplos", tcontext.rootname);
    fflush(stdout);
    if (!(fp = fopen(s, "w")))
      bombElegant("unable to open debug file for aperture losses", NULL);
    fprintf(fp, "SDDS1\n");
    fprintf(fp, "&column name=z type=double units=m &end\n");
    fprintf(fp, "&column name=x type=double units=m &end\n");
    fprintf(fp, "&column name=y type=double units=m &end\n");
    fprintf(fp, "&data mode=ascii no_row_counts=1 &end\n");
  }
#endif

  if (!interpolateApertureData(z, apData,
                               &xCenter, &yCenter, &xSize, &ySize))
    return np;

  struct GPUBASE* gpuBase = getGpuBase();
  unsigned int* d_sortIndex = gpuBase->d_tempu_alpha;
  np = killParticles(np, d_sortIndex, accepted,
         gpu_imposeApertureData_kernel(xCenter, yCenter, xSize, ySize, 
           z, Po, d_sortIndex));

#if DEBUG_APERTURE
  fflush(fp);
#endif

  return(np);
}

} //extern "C"
