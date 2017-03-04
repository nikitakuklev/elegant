#include <gpu_track.h>

#include <curand.h>
#include <curand_kernel.h>

#undef align
#include <gpu_killParticles.hcu>
#include <gpu_particle_template_function.hcu>
#include <gpu_matter.hcu> // device code shared with scraper in limit_amplitudes
#include <gpu_base.h>

#define AMU (1.6605e-27)
#define ALPHA (1./137.036)

class gpu_track_through_matter_kernel{
public:
  double Nrad, L, theta_rms, prob, L1, K2, probBS, probER, Po, dGammaFactor, z0;
  int sections0, multipleScattering;
  unsigned int *d_sortIndex;
  curandState_t *d_state;
  /* from MATTER struct */
  int energyDecay, nuclearBremsstrahlung, energyStraggle;
  int nSlots;
  double spacing,  width,  tilt,  center;

 gpu_track_through_matter_kernel(
    int energyDecay, int nuclearBremsstrahlung, int energyStraggle,
    int multipleScattering, int nSlots, double spacing, double width,
    double tilt, double center, double Nrad, double L, double theta_rms,
    int sections0, double prob, double L1, double K2,
    double probBS, double probER, double Po, double dGammaFactor,
    double z0, unsigned int *d_sortIndex, curandState_t *d_state) :
    energyDecay(energyDecay), nuclearBremsstrahlung(nuclearBremsstrahlung),
    energyStraggle(energyStraggle), nSlots(nSlots), spacing(spacing), width(width),
    tilt(tilt), center(center), multipleScattering(multipleScattering),
    Nrad(Nrad), L(L), theta_rms(theta_rms), sections0(sections0),
    prob(prob), L1(L1), K2(K2), probBS(probBS), probER(probER), Po(Po),
    dGammaFactor(dGammaFactor), z0(z0), d_sortIndex(d_sortIndex),
    d_state(d_state) {};

  __device__ unsigned int operator()(gpuParticleAccessor& part){
    unsigned int tid = part.getParticleIndex();

    return
      gpu_track_through_matter_dfunc(part, energyDecay,
        nuclearBremsstrahlung, energyStraggle, multipleScattering,
        nSlots, spacing,  width,  tilt,  center,
        Nrad, L, theta_rms, sections0, prob, L1, K2, probBS, probER,
        Po, dGammaFactor, z0, d_sortIndex, d_state);
  }
};

extern "C" {

double radiationLength(long Z, double A, double rho);

void gpu_set_track_through_matter(long np, MATTER *matter, double Po,
       int impulseMode, int *multipleScattering, double *Nrad, double *L,
       double *theta_rms, int *sections0, double *prob, double *L1, double *K2,
       double *probBS, double *probER, double *dGammaFactor){

  double sigmaTotal, probScatter=0.0;
  double K1, Xo, probBSScatter = 0, probERScatter=0, beta;

  if (particleIsElectron==0)
    bombElegant("MATTER element doesn't work for particles other than electrons", NULL);

  if (matter->energyDecay && (matter->nuclearBremsstrahlung || matter->electronRecoil))
    bombElegant("ENERGY_DECAY=1 and NUCLEAR_BREHMSSTRAHLUNG=1 or ELECTRON_RECOIL=1 options to MATTER/SCATTER element are mutually exclusive", NULL);

  beta = Po/sqrt(sqr(Po)+1);
  if (matter->Xo==0) {
    if (matter->Z<1 || matter->A<1 || matter->rho==0)
      bombElegant("XO=0 but Z, A, or rho invalid for MATTER element", NULL);
    Xo = radiationLength(matter->Z, matter->A, matter->rho);
    /* printf("Computed radiation length for Z=%ld, A=%le, rho=%le is %le m\n",
           matter->Z, matter->A, matter->rho, Xo);
    */
  } else 
    Xo = matter->Xo;
  
  *Nrad = *L/Xo;
  *dGammaFactor = 1-exp(-(*Nrad));
  *prob = *probBS = *probER = 0;
  if (*Nrad<1e-3 || matter->nuclearBremsstrahlung || matter->electronRecoil) {
    if (matter->Z<1 || matter->A<1 || matter->rho<=0)
      bombElegant("MATTER element is too thin or requests special features---provide Z, A, and rho for single-scattering calculation.", NULL);
    K1 = 4*matter->Z*(matter->Z+1)*sqr(particleRadius/(beta*Po));
    *K2 = sqr(pow(matter->Z, 1./3.)*ALPHA/Po);
    sigmaTotal = K1*pow(PI, 3)/(sqr(*K2)+(*K2)*SQR_PI);
    probScatter = matter->rho/(AMU*matter->A)*(*L)*sigmaTotal;
    /* printf("K1=%le, K2=%le, mean expected number of scatters is %le\n", K1, K2, probScatter); */
    probBSScatter = 0;
    if (matter->nuclearBremsstrahlung) {
      probBSScatter = 4*(*L)/(3*Xo)*(-log(BS_Y0)-(1-BS_Y0)+3./8.*(1-BS_Y0*BS_Y0));
    }
    if (matter->electronRecoil) {
      probERScatter = (*L)*matter->rho/(AMU*matter->A)*PIx2*matter->Z*sqr(re_mks)/Po*(1/BS_Y0-1);
    }
    *sections0 = probScatter/matter->pLimit+1;
    *Nrad /= *sections0;
    *multipleScattering = 0;    
    *L1 = (*L)/(*sections0);
    *prob = probScatter/(*sections0);
    *probBS = probBSScatter/(*sections0);
    *probER = probERScatter/(*sections0);
    printf("Sections=%d, L1 = %le, probIS = %le, probBS = %le, probER = %le\n", *sections0, *L1, *prob, *probBS, *probER);
  } else {
    *multipleScattering = 1;
    *theta_rms = 13.6/particleMassMV/Po/sqr(beta)*sqrt(*Nrad)*(1+0.038*log(*Nrad));
  }
  
  if (impulseMode)
    *L = *L1 = 0;

  if (matter->spacing>0 && matter->width>0) {
    if (matter->spacing<=matter->width)
        bombElegant("MATTER SPACING parameter is less than WIDTH parameter", NULL);
    if (matter->spacing<=0)
      bombElegant("MATTER SPACING parameter is <= 0", NULL); 
  }
}

long gpu_track_through_matter(long np, MATTER *matter, double Po, 
                              double **accepted, double z0)
{
  double L, Nrad, theta_rms=0;
  double dGammaFactor;
  double K2=0.0;
  int sections0=1;
  double L1, prob, probBS, probER;  
  int multipleScattering = 0, impulseMode;
  
  log_entry("track_through_matter");
  
  if (matter->length!=0) {
    L = matter->length;
    impulseMode = 0;
  } else if (matter->lEffective!=0) {
    L = matter->lEffective;
    impulseMode = 1;
  }
  else 
    return np;

  gpu_set_track_through_matter(np, matter, Po, impulseMode, &multipleScattering, &Nrad, 
       &L, &theta_rms, &sections0, &prob, &L1, &K2, &probBS, &probER, &dGammaFactor);
  gpu_setMatterTables();

  struct GPUBASE* gpuBase = getGpuBase();
  unsigned int particlePitch = gpuBase->gpu_array_pitch;
  double* d_gauss_rn =  gpuBase->d_temp_particles + particlePitch;
  curandState_t* d_state = (curandState_t*)gpu_get_rand_state(d_gauss_rn, np, random_2(0));
  unsigned int* d_sortIndex = gpuBase->d_tempu_alpha;
  np = killParticles(np, d_sortIndex, accepted,
    gpu_track_through_matter_kernel(matter->energyDecay, matter->nuclearBremsstrahlung,
    matter->energyStraggle, matter->nSlots, matter->spacing, matter->width, 
    matter->tilt, matter->center,
    multipleScattering, Nrad, L, theta_rms, sections0,
    prob, L1, K2, probBS, probER, Po, dGammaFactor, z0, d_sortIndex, d_state));
  
  log_exit("track_through_matter");
  return np;
}

} // extern C
