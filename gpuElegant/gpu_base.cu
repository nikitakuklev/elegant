/**
 * Base singleton struct with pointers to data used by the GPU routines
 * and some smaller GPU routines used in do_tracking.c
 */
#include <iostream>
#include <string>

#include <curand.h>
#include <curand_kernel.h>

#include <thrust/device_ptr.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/gather.h>
#include <thrust/remove.h>
#include <thrust/sort.h>
#include <thrust/scan.h>

#include <gpu_track.h>

#include <gpu_base.h>

extern "C" {

#define NSTREAMS 7 // # of CUDA streams

static bool deallocGpuBase=false;
static bool deallocGpuLocal=false;
// for elemName 
static char initName[5] = "init";
static char noneName[5] = "none";
static char unknownName[8] = "unknown";

struct GPUBASE* getGpuBase(){
  static struct GPUBASE* gpuBase = NULL;
  if (!gpuBase && !deallocGpuBase) {
    gpuBase = (GPUBASE*)malloc(sizeof(GPUBASE));
    gpuBase->d_particles = NULL;
    gpuBase->d_accepted = NULL;
    gpuBase->d_lostParticles = NULL;
    gpuBase->d_temp_particles = NULL;
    gpuBase->d_tempu_alpha = NULL;
    gpuBase->d_blockTempf = NULL;
    gpuBase->d_blockTempu = NULL;
    gpuBase->d_retirementCount = NULL;
    gpuBase->d_csbend=NULL;
    gpuBase->coord = NULL;
    gpuBase->accepted = NULL;
    gpuBase->lostParticles = NULL;
    gpuBase->doh_particles = NULL;
    gpuBase->gpu_array_pitch = 0;
    gpuBase->particlesOnGpu = 0;
    gpuBase->nOriginal = 0;
    gpuBase->n_comp = 0;
    gpuBase->elementOnGpu = 0;
    gpuBase->copyForValidation = 0;
    gpuBase->gpuFtnLevel = 0;
    gpuBase->elemType = -100;
    gpuBase->elemName = initName;
    gpuBase->nBins = 0;
    gpuBase->isMaster = 0;
    gpuBase->nThreads = 1;
    gpuBase->maxBlocks = 1;
    gpuBase->nReductionThreads = 1;
    gpuBase->nReductionBlocks = 1;
    gpuBase->full_d_particles = NULL;
    gpuBase->particlesUnsorted = 0;
  } else if (gpuBase && deallocGpuBase) {
    free(gpuBase);
    gpuBase=NULL;
    deallocGpuBase=false;
  }
  return gpuBase;
}

unsigned int getElementOnGpu(){
  return getGpuBase()->elementOnGpu;
}

/**
 * Boolean to indicate if it is a top level GPU function call for timing.
 * This is used to avoid recursive timing comparisons.
 */
unsigned int topLevelGpuFtn(){
  return (getGpuBase()->gpuFtnLevel==1);
}

void gpuBaseInit(double** coord, unsigned int nOriginal, double** accepted,
                 double** lostParticles, long isMaster){
  struct GPUBASE* gpuBase = getGpuBase();
  gpuBase->isMaster=isMaster;

  /* Adjust global thread/block decomposition based on device prop
     only 256, 512, 1024 and 2048 are supported for nReductionBlocks */
  cudaDeviceProp props;
  cudaGetDeviceProperties(&props, 0);
  if (props.major>=3) {
    gpuBase->nThreads = 512;
    gpuBase->maxBlocks = 4*props.multiProcessorCount;
    gpuBase->nReductionThreads = 256;
    gpuBase->nReductionBlocks = 3*props.multiProcessorCount;
  } else {
    gpuBase->nThreads = 256;
    gpuBase->maxBlocks = 8*props.multiProcessorCount;
    gpuBase->nReductionThreads = 256;
    gpuBase->nReductionBlocks = 3*props.multiProcessorCount;
  }
  static bool report=true;
  if(gpuBase->isMaster && report) {
    std::cout << "gpuBaseInit: Using cuda device of compute capability " 
      << props.major << "." << props.minor 
      << ".\ngpuBaseInit: Particle kernels use " << gpuBase->nThreads 
      << " threads and " << gpuBase->maxBlocks << " blocks"
      << ".\ngpuBaseInit: Reductions use " << gpuBase->nReductionThreads 
      << " threads and " << gpuBase->nReductionBlocks << " blocks." 
      << std::endl;
  }

  /* Use separate Devices for each MPI process.
     The master node does not own any particles and thus does 
     not use a device */
#if USE_MPI
  cudaError_t err;
  int deviceCount, rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
#ifdef CUDA_DEVICE_COUNT
  deviceCount = CUDA_DEVICE_COUNT;
#else
  err = cudaGetDeviceCount(&deviceCount); 
  gpuErrorHandler(err,"gpuBaseInit: cudaGetDeviceCount");
#endif
  if(gpuBase->isMaster && report) {
    std::cout << "gpuBaseInit: Number of devices per node is "
              << deviceCount << "." << std::endl;
    report=false;
  }
  if(nOriginal==0) {
    return;
  }
  int deviceNum = rank % deviceCount; 
  if (deviceNum>=0) err = cudaSetDevice(deviceNum);
  gpuErrorHandler(err,"gpuBaseInit: cudaSetDevice");
#else  /* USE_MPI */
  report=false;
#endif /* USE_MPI */

  if (props.major<3)
    cudaSetDeviceFlags(cudaDeviceScheduleBlockingSync);

  gpuBase->coord = coord;
  gpuBase->accepted = accepted;
  gpuBase->lostParticles = lostParticles;
  gpuBase->nOriginal = nOriginal;
  if(coord && nOriginal > 1){
    gpuBase->n_comp = (unsigned int) ( (coord[1] - coord[0]) );
    if(gpuBase->n_comp != 7) {
      std::cout << "gpuElegant: unknown n_comp=" << gpuBase->n_comp 
                << ", np=" << nOriginal << std::endl;
      exit(-1);
    }
  }else{
    gpuBase->n_comp = 7;
  }
  allocateGpuParticles(&(gpuBase->d_particles));
  if (accepted)
    allocateGpuParticles(&(gpuBase->d_accepted));
  if (lostParticles)
    allocateGpuParticles(&(gpuBase->d_lostParticles));
  allocateGpuParticles(&(gpuBase->d_temp_particles));
  allocateGpuUnsignedIntArray(&(gpuBase->d_tempu_alpha));
  
  // 2x because double reductions need double memory
  cudaMalloc( (void**)&(gpuBase->d_blockTempf), 
      2*NSTREAMS*sizeof(double)*gpuBase->nReductionBlocks);
  cudaMalloc( (void**)&(gpuBase->d_blockTempu), 
      2*NSTREAMS*sizeof(unsigned int)*gpuBase->nReductionBlocks);
  // 128 bytes for alignment 
  cudaMalloc( (void**)&(gpuBase->d_retirementCount), 128*NSTREAMS);

#if GPU_VERIFY
  // initialize random particles to avoid initialization costs in timings
  gpu_get_rand_state(gpuBase->d_temp_particles, gpuBase->nOriginal, 10101010.);
  // initialize streams to avoid initialization costs in timings
  finishReductionStreams();
#endif
}

void deallocateGpuLocal();

void gpuBaseDealloc(){
  struct GPUBASE* gpuBase = getGpuBase();
  cudaError_t err;

  if (gpuBase->d_particles==NULL) {
#if USE_MPI
    if (gpuBase->nOriginal==0) return;
#endif
    std::cout << "ERROR: gpuBaseDealloc() called without gpuBaseInit()"  
      << std::endl;
    return;
  }

  /* Copy particles back to the CPU before deallocation. Use isMaster 
     as a flag to indicate this is a finalize step */
  gpuBase->isMaster=-1;
  if (gpuBase->particlesOnGpu && gpuBase->coord)
    copyParticlesFromGpuToCpu(gpuBase->nOriginal);
 

  err = cudaFree( (void*)gpuBase->d_retirementCount);
  gpuErrorHandler(err,"gpuBaseDealloc: cudaFree d_retirementCount");
  err = cudaFree( (void*)gpuBase->d_blockTempu);
  gpuErrorHandler(err,"gpuBaseDealloc: cudaFree d_blockTempu");
  err = cudaFree( (void*)gpuBase->d_blockTempf);
  gpuErrorHandler(err,"gpuBaseDealloc: cudaFree d_blockTempf");
  err = cudaFree( (void*)gpuBase->d_tempu_alpha);
  gpuErrorHandler(err,"gpuBaseDealloc: cudaFree d_tempu_alpha");
  freeGpuParticles(&(gpuBase->d_temp_particles));
  freeGpuParticles(&(gpuBase->d_particles));
  if (gpuBase->d_accepted)
    freeGpuParticles(&(gpuBase->d_accepted));
  if (gpuBase->d_lostParticles)
    freeGpuParticles(&(gpuBase->d_lostParticles));

  if(gpuBase->d_csbend) {
    err = cudaFree( (void*)gpuBase->d_csbend);
    gpuErrorHandler(err,"gpuBaseDealloc: cudaFree d_csbend");
  }

#ifdef GPU_VERIFY
  finalizeTimer();
#endif
  deallocateGpuLocal();
  
  deallocGpuBase=true;
  getGpuBase(); // deallocate
}

/**
 * Determine if element can be computed on the GPU
 * Return 1/0 for true/false and 2 if the element does not
 * involve the particle arrays and thus particles should 
 * remain on either the host or device.
 * @param veptr pointer to the eptr linked list of elements
 */
unsigned int isElementOnGpu(void* veptr){

  ELEMENT_LIST* eptr = (ELEMENT_LIST*)veptr;
  if( eptr == NULL || eptr == 0 )
    return 0;

  switch(eptr->type){ // TYPEDEF #
  case T_QUAD: return 1;          // 1
  case T_SBEN: return 1;          // 2
  case T_DRIF: return 1;          // 4
  case T_SEXT: return 1;          // 5
  case T_SOLE: return 1;          // 8
  case T_VCOR: return 1;          // 9
  case T_HCOR: return 1;          // 10
  case T_RFCA: return 1;          // 11
  case T_HMON: return 1;          // 13
  case T_VMON: return 1;          // 14
  case T_MONI: return 1;          // 15
  case T_RCOL: return 1;          // 16
  case T_ECOL: return 1;          // 17
  case T_MARK: return 1;          // 18
  case T_MATR: return 1;          // 19
  case T_MALIGN: return 1;        // 28
  case T_ENERGY: return 1;        // 31
  case T_MAXAMP: return 2;        // 32
  case T_RECIRC: return 2;        // 35
  case T_SCRAPER: return 1;       // 37
  case T_CENTER: return 1;        // 38
  case T_KSEXT: return 1;         // 40
  case T_KQUAD: return 1;         // 42
  case T_HVCOR: return 1;         // 45
  case T_CSBEND: return 1;        // 53
  case T_MATTER: return 1;        // 55
  case T_WAKE: return 1;          // 65
  case T_TRWAKE: return 1;        // 66
  case T_CHARGE: return 2;        // 68
  case T_CSRCSBEND: return 1;     // 71
  case T_CSRDRIFT: return 1;      // 72
  case T_RFCW: return 1;          // 73
  case T_FLOORELEMENT: return 2;  // 81
  case T_EMATRIX: return 1;       // 84
  case T_LSCDRIFT: return 1;      // 89
  case T_EDRIFT: return 1;        // 95
  case T_KQUSE: return 1;         // 99
  case T_KOCT: return 1;          // 105
  case T_MRADINTEGRALS: return 2; // 106
  default: return 0;
  }
}

void setElementGpuData(void* veptr) {
  ELEMENT_LIST* eptr = (ELEMENT_LIST*)veptr;
  struct GPUBASE* gpuBase = getGpuBase();

  if(eptr==NULL) {
    gpuBase->elementOnGpu = 0;
    gpuBase->elemType = T_RETURN;
    gpuBase->elemName = noneName;
    return;
  }

  unsigned int nextElementFlag = isElementOnGpu(veptr);
  if (nextElementFlag<2) gpuBase->elementOnGpu = nextElementFlag;
  //if (nOriginal<500) gpuBase->elementOnGpu = 0; //particle limit
#if DEBUG || CUDA_DEBUG
  if (gpuBase->elementOnGpu)
    std::cout << "Next GPU element is " << eptr->type << std::endl;
  else
    std::cout << "Next CPU element is " << eptr->type << std::endl;
#endif
  gpuBase->elemType = eptr->type;
  gpuBase->elemName = eptr->name?eptr->name:unknownName;

  if(gpuBase->elementOnGpu && !gpuBase->particlesOnGpu)
    copyParticlesFromCpuToGpu(gpuBase->nOriginal);
  else if(!gpuBase->elementOnGpu && gpuBase->particlesOnGpu)
    copyParticlesFromGpuToCpu(gpuBase->nOriginal);
}

double** forceParticlesToCpu(char* routine) {
  struct GPUBASE* gpuBase = getGpuBase();
  if(gpuBase->particlesOnGpu) {
    copyParticlesFromGpuToCpu(gpuBase->nOriginal);
    //std::cout << "Copied particles to the CPU for " << routine << std::endl;
    gpuBase->elementOnGpu = 0;
    gpuBase->particlesOnGpu = 0;
  }
  return gpuBase->coord;
}

void syncWithGpuLossBuffer(long nToTrack, long nLeft) {
  struct GPUBASE* gpuBase = getGpuBase();
  unsigned int particlePitch = gpuBase->gpu_array_pitch;
  double** lossBuffer = gpuBase->lostParticles;
  double* d_lossBuffer = gpuBase->d_lostParticles;
  unsigned int n_comp = gpuBase->n_comp;
  cudaError_t err;
  
  long nLostSync = nToTrack - nLeft;
  double *tLossBuffer = (double*)malloc(sizeof(double)*nLostSync*n_comp);

  err = cudaMemcpy(tLossBuffer, d_lossBuffer + nLeft,
                   sizeof(double)*n_comp*nLostSync, cudaMemcpyDeviceToHost);
  gpuErrorHandler(err,"syncWithGpuLossBuffer::cudaMemcpy");
  
  for (long ip=nLeft; ip<nToTrack; ip++)  
    for (int j=0; j<7; j++)
      lossBuffer[ip][j] = tLossBuffer[ip-nLeft + j*particlePitch];

  free(tLossBuffer);
}

/**
 * Struct to hold GPU data local to this file.
 */
struct GPULOCAL{
  curandState_t* cuRandState;
  curandGenerator_t* cuRandGen;
  unsigned int nstreams;
  unsigned int nextStream;
  cudaStream_t* streams;
  double* pinnedVal;
  double* pinnedVal2;
  unsigned int* streamActive;
  double** outputValue;
  double** outputValue2;
};

/**
 * Returns a Singleton instance of gpulocal. Initializes if needed.
 */
struct GPULOCAL* getGpuLocal(){
  static struct GPULOCAL* gpuLocal = NULL;
  unsigned int ii;
  cudaError_t err;

  if (!gpuLocal && !deallocGpuLocal) {
    gpuLocal = (GPULOCAL*)malloc(sizeof(GPULOCAL));
    gpuLocal->cuRandState = NULL;
    gpuLocal->cuRandGen = NULL;
    gpuLocal->nstreams = NSTREAMS; 
    unsigned int nstreams = gpuLocal->nstreams;
    gpuLocal->nextStream = 0;
    gpuLocal->streams = (cudaStream_t*)malloc(nstreams*sizeof(cudaStream_t));
    for(ii = 0; ii < nstreams; ii++) {
      err = cudaStreamCreate(&(gpuLocal->streams[ii]));
      gpuErrorHandler(err, "getGpuLocal: cudaStreamCreate failed");
    }
    err = cudaHostAlloc(&gpuLocal->pinnedVal, nstreams*sizeof(double),
                        cudaHostAllocPortable);
    gpuErrorHandler(err, "getGpuLocal: cudaHostAlloc pinnedVal failed");
    err = cudaHostAlloc(&gpuLocal->pinnedVal2, nstreams*sizeof(double),
                        cudaHostAllocPortable);
    gpuErrorHandler(err, "getGpuLocal: cudaHostAlloc pinnedVal2 failed");
    gpuLocal->streamActive = 
      (unsigned int*)malloc(nstreams*sizeof(unsigned int));
    for(ii = 0; ii < nstreams; ii++) gpuLocal->streamActive[ii]=0;
    gpuLocal->outputValue = (double**)malloc(nstreams*sizeof(double*));
    gpuLocal->outputValue2 = (double**)malloc(nstreams*sizeof(double*));
    for(ii = 0; ii < nstreams; ii++) gpuLocal->outputValue2[ii] = NULL;
  } else if (gpuLocal && deallocGpuLocal) {
    free(gpuLocal);
    gpuLocal=NULL;
    deallocGpuLocal=false;
  }
  return gpuLocal;
}

void deallocateGpuLocal() {
  struct GPULOCAL* gpuLocal = getGpuLocal();
  cudaError_t err;
  for(unsigned int i = 0; i < gpuLocal->nstreams; i++) {
    err = cudaStreamDestroy(gpuLocal->streams[i]);
    gpuErrorHandler(err, "deallocateGpuLocal: cudaStreamDestroy failed");
  }
  free(gpuLocal->streams);
  err = cudaFreeHost(gpuLocal->pinnedVal);
  gpuErrorHandler(err, "deallocateGpuLocal: cudaFreeHost pinnedVal failed");
  err = cudaFreeHost(gpuLocal->pinnedVal2);
  gpuErrorHandler(err, "deallocateGpuLocal: cudaFreeHost pinnedVal2 failed");
  free(gpuLocal->streamActive);
  free(gpuLocal->outputValue);
  free(gpuLocal->outputValue2);
  if (gpuLocal->cuRandGen) {
    free(gpuLocal->cuRandGen);
  }
  if (gpuLocal->cuRandState) {
    err = cudaFree( gpuLocal->cuRandState);
    gpuErrorHandler(err, "gpuBaseDealloc: cudaFree cuRandState");
  }

  deallocGpuLocal=true;
  getGpuLocal(); // deallocate
}

void* getNextReductionStream(double** d_tempf, 
    unsigned int** d_retirementCount, double* value, double** pinval,
    double* value2, double** pinval2){
  GPUBASE* gpuBase = getGpuBase();
  GPULOCAL* gpuLocal = getGpuLocal();
  unsigned int istream = gpuLocal->nextStream; 
  gpuLocal->nextStream++;
  gpuLocal->nextStream%=gpuLocal->nstreams;
  cudaStream_t* stream = &gpuLocal->streams[istream];
  cudaError_t err;

  // Finish and copy from pinned memory if stream is active
  if(gpuLocal->streamActive[istream]) {
    err = cudaStreamSynchronize(*stream);
    gpuErrorHandler(err, 
        "getNextReductionStream: cudaStreamSynchronize failed");
    *(gpuLocal->outputValue[istream]) = gpuLocal->pinnedVal[istream];
    if(gpuLocal->outputValue2[istream]) {
      *(gpuLocal->outputValue2[istream]) = gpuLocal->pinnedVal2[istream];
      gpuLocal->outputValue2[istream] = NULL;
    }
  }

  // Setup next stream arrays
  *d_tempf = gpuBase->d_blockTempf+2*istream*gpuBase->nReductionBlocks;
  *d_retirementCount = gpuBase->d_retirementCount+istream*128/sizeof(unsigned int);
  err = cudaMemsetAsync(*d_retirementCount, 0, 128, *stream);
  gpuErrorHandler(err, "getNextReductionStream: cudaMemsetAsync failed");

  *pinval = &(gpuLocal->pinnedVal[istream]);
  gpuLocal->outputValue[istream] = value;
  if(value2) {
    *pinval2 = &(gpuLocal->pinnedVal2[istream]);
    gpuLocal->outputValue2[istream] = value2;
  }

  gpuLocal->streamActive[istream] = 1;
  return (void*)stream;
}

void finishReductionStreams() {
  GPULOCAL* gpuLocal = getGpuLocal();
  for (unsigned int istream = 0; istream < gpuLocal->nstreams; istream++) {
    if(gpuLocal->streamActive[istream]) {
      gpuLocal->streamActive[istream] = 0;
      cudaStream_t* stream = &gpuLocal->streams[istream];
      cudaError_t err = cudaStreamSynchronize(*stream);
      gpuErrorHandler(err, 
          "getNextReductionStream: cudaStreamSynchronize failed");
      // Copy results from pinned memory
      *(gpuLocal->outputValue[istream]) = gpuLocal->pinnedVal[istream];
      if(gpuLocal->outputValue2[istream]) {
        *(gpuLocal->outputValue2[istream]) = gpuLocal->pinnedVal2[istream];
        gpuLocal->outputValue2[istream] = NULL;
      }
    }
  }
}

/**
 * Handle cuRAND errors
 */
inline void gpuCuRandErrorHandler(int err, std::string diagnosis){
  if(err != 0)
    std::cout << "cuRand Kernel " << diagnosis 
      << " returned with error " << err << std::endl;
}

/**
 * Initialize the cuRAND random number generator
 * @param seed random number generator seed
 */
void gpu_init_rand_generator(double seed) {
  struct GPULOCAL* gpuLocal = getGpuLocal();
  curandGenerator_t* generator = gpuLocal->cuRandGen;
  if (generator) return;
  
  /* Do not modify the bits on conversion to a long long for the seed */
  unsigned long long iseed=*reinterpret_cast<unsigned long long*>(&seed);

  curandStatus_t errRand;
  generator = (curandGenerator_t*)malloc(sizeof(curandGenerator_t));
  errRand = curandCreateGenerator(generator, CURAND_RNG_PSEUDO_DEFAULT);
  gpuCuRandErrorHandler((int)errRand,
      "gpu_init_rand_generator: curandCreateGenerator");
  errRand = curandSetPseudoRandomGeneratorSeed(*generator, iseed);
  gpuCuRandErrorHandler((int)errRand,
      "gpu_init_rand_generator: curandSetPseudoRandomGeneratorSeed");
  errRand = curandSetGeneratorOrdering(*generator, 
                                       CURAND_ORDERING_PSEUDO_SEEDED);
  gpuCuRandErrorHandler((int)errRand,
      "gpu_init_rand_generator: curandSetGeneratorOrdering");
  gpuLocal->cuRandGen = generator;
}

void gpu_d_gauss_rn(double* d_ranarr, unsigned int n_num, double mean,
                    double sigma, double seed) {

  gpu_init_rand_generator(seed);

  struct GPULOCAL* gpuLocal = getGpuLocal();
  curandGenerator_t* generator = gpuLocal->cuRandGen;
  curandStatus_t errRand;

  /* requires an even n_num */
  if (n_num%2>0) n_num+=1;

  errRand = curandGenerateNormalDouble(*generator, d_ranarr,
                                       n_num, mean, sigma);
  gpuCuRandErrorHandler((int)errRand,
      "gpu_d_gauss_rn: curandGenerateNormalDouble");
}

/**
 * Setup curand state vector.
 * @param state array of curandState_t to initialize
 */
__global__ void setupState(double* d_ranarr, unsigned int np, 
                           curandState_t* state) {
  for(unsigned int tid=threadIdx.x + blockIdx.x*blockDim.x; tid < np; 
      tid += blockDim.x*gridDim.x) {
    /* Do not modify the bits on conversion to a long long for the seed */
    curand_init(*(unsigned long long*)&d_ranarr[tid],tid,0,state+tid);
  }
}

/**
 * Initialize the random cuRAND random state array
 * The seed is only used for the initial generation, after that
 * the curand generator is responsible for changing the seed.
 * @param d_ranarr pointer to the random number array stored in the 
 *                 global device memory
 * @param n_num number of random numbers to generate
 * @param seed random number generator seed
 */
void gpu_init_rand_state(double* d_ranarr, unsigned int n_num, 
                         double seed) {
  struct GPULOCAL* gpuLocal = getGpuLocal();
  curandState_t* state = gpuLocal->cuRandState;
  if (state) return;

  struct GPUBASE* gpuBase = getGpuBase();
  unsigned int nOriginal = gpuBase->nOriginal;

  // Compute the state seed
  gpu_init_rand_generator(seed);
  curandStatus_t errRand;
  curandGenerator_t* generator = gpuLocal->cuRandGen;
  errRand = curandGenerateUniformDouble(*generator, d_ranarr,
                                        nOriginal);
  gpuCuRandErrorHandler((int)errRand,
      "gpu_init_rand_state: curandGenerateUniformDouble");

  // setup the state
  cudaMalloc((void **)&state, nOriginal*sizeof(curandState));
  unsigned int nTx = 256;
  unsigned int nBx = (n_num + nTx - 1) / nTx;
  if(nBx > gpuBase->nReductionBlocks) nBx = gpuBase->nReductionBlocks;
  setupState<<<nBx,nTx>>>(d_ranarr, nOriginal, state);
  gpuLocal->cuRandState = state;

  gpuErrorHandler("gpu_init_rand_state: setupState");
}

void* gpu_get_rand_state(double* d_ranarr, unsigned int n_num, double seed) {
  struct GPUBASE* gpuBase = getGpuBase();
  unsigned int nOriginal = gpuBase->nOriginal;
  if (n_num > nOriginal) 
    std::cout << "Warning: n_num too large in gpu_init_rand_state "
              << "- state array sized at nOriginal."
              << std::endl;

  struct GPULOCAL* gpuLocal = getGpuLocal();
  curandState_t* state = gpuLocal->cuRandState;
  if (state) return (void*)state;

  gpu_init_rand_state(d_ranarr, n_num, seed);
  state = gpuLocal->cuRandState;

  return (void*)state;
}

/**
 * Ensure random numbers are bounded within a limit.
 * @param n_num number of random numbers to generate
 */
__global__ void enforceLimit(double* d_ranarr, unsigned int np, double mean, 
                             double sigma, double limit, curandState_t* state){
  for(unsigned int tid=threadIdx.x + blockIdx.x*blockDim.x; tid < np; 
      tid += blockDim.x*gridDim.x) {
    while(fabs(d_ranarr[tid]-mean)>limit) {
      d_ranarr[tid]=sigma*curand_normal_double(&state[tid])+mean;
    }
  }
}

/**
 * Generate a random normal distribution within a limit.
 * The seed is only used for the initial generation, after that
 * the curand generator is responsible for changing the seed.
 * @param d_ranarr pointer to the random number array stored in the 
 *                 global device memory
 * @param n_num number of random numbers to generate
 * @param mean normal distribution mean
 * @param sigma normal distribution standard deviation
 * @param limit_in_sigmas cutoff in stdevs of the normal distribtion
 * @param seed random number generator seed
 */
void gpu_d_gauss_rn_lim(double* d_ranarr, unsigned int n_num, double mean, 
                        double sigma, double limit_in_sigmas, double seed) {

  curandState_t* state =
    (curandState_t*)gpu_get_rand_state(d_ranarr, n_num, seed);

  gpu_d_gauss_rn(d_ranarr, n_num, mean, sigma, seed);

  unsigned int nTx = 256;
  unsigned int nBx = (n_num + nTx - 1) / nTx;
  struct GPUBASE* gpuBase = getGpuBase();
  if(nBx > gpuBase->nReductionBlocks) nBx = gpuBase->nReductionBlocks;
  enforceLimit<<<nBx,nTx>>>(d_ranarr, n_num, mean, 
                            sigma, limit_in_sigmas*sigma, state);
  gpuErrorHandler("gpu_d_gauss_rn_lim: enforceLimit");
}

/**
 * Set the d_particles array pointer with an offset
 * @param offset offset index
 */
void selectParticleSubset(unsigned int offset) {
  struct GPUBASE* gpuBase = getGpuBase();
  if (gpuBase->full_d_particles == NULL) {
    gpuBase->full_d_particles = gpuBase->d_particles;
    gpuBase->d_particles = gpuBase->d_particles + offset;
  } else {
    gpuBase->d_particles = gpuBase->full_d_particles + offset;
  }
}

/**
 * Reset the d_particles array pointer
 */
void resetToFullParticleArray() {
  struct GPUBASE* gpuBase = getGpuBase();
  if (gpuBase->full_d_particles != NULL) {
    gpuBase->d_particles = gpuBase->full_d_particles;
    gpuBase->full_d_particles = NULL;
  }
}

/**
 * Sort particles by ID
 * @param np number of particles
 */
void sortByPID(unsigned int np) {
  struct GPUBASE* gpuBase = getGpuBase();
  double* d_particles = gpuBase->d_particles;
  double* d_accepted = gpuBase->d_accepted;
  unsigned int particlePitch = gpuBase->gpu_array_pitch;
  double* d_temp_particles = gpuBase->d_temp_particles;
  unsigned int* d_linearIndex = gpuBase->d_tempu_alpha;

  /* Don't sort unless we have to */
  if (!gpuBase->particlesUnsorted) return;

  /* Create iterators over the device arrays */
  thrust::device_ptr<double> dp_0(d_particles + 0*particlePitch);
  thrust::device_ptr<double> dp_1(d_particles + 1*particlePitch);
  thrust::device_ptr<double> dp_2(d_particles + 2*particlePitch);
  thrust::device_ptr<double> dp_3(d_particles + 3*particlePitch);
  thrust::device_ptr<double> dp_4(d_particles + 4*particlePitch);
  thrust::device_ptr<double> dp_5(d_particles + 5*particlePitch);
  thrust::device_ptr<double> dp_6(d_particles + 6*particlePitch);

  thrust::device_ptr<double> dtp_0(d_temp_particles + 0*particlePitch);
  thrust::device_ptr<double> dtp_1(d_temp_particles + 1*particlePitch);
  thrust::device_ptr<double> dtp_2(d_temp_particles + 2*particlePitch);
  thrust::device_ptr<double> dtp_3(d_temp_particles + 3*particlePitch);
  thrust::device_ptr<double> dtp_4(d_temp_particles + 4*particlePitch);
  thrust::device_ptr<double> dtp_5(d_temp_particles + 5*particlePitch);
  thrust::device_ptr<double> dtp_6(d_temp_particles + 6*particlePitch);

  typedef thrust::zip_iterator< thrust::tuple< thrust::device_ptr<double>,
          thrust::device_ptr<double>,thrust::device_ptr<double>,
          thrust::device_ptr<double>,thrust::device_ptr<double>,
          thrust::device_ptr<double>,thrust::device_ptr<double> > >
          dp7_zip_iterator;
 
  dp7_zip_iterator zipInit =  make_zip_iterator( 
           thrust::make_tuple(dp_0, dp_1, dp_2, dp_3, dp_4, dp_5, dp_6) );
  dp7_zip_iterator zipFinal = make_zip_iterator( 
           thrust::make_tuple(dtp_0, dtp_1, dtp_2, dtp_3, dtp_4, dtp_5, dtp_6) );

  /* thrust sort particles */
  thrust::device_ptr<double> tdp_sortIndex(d_particles + 6*particlePitch);
  thrust::device_ptr<unsigned int> tdp_linearIndex(d_linearIndex);
  thrust::counting_iterator<unsigned int> counter(0);
  thrust::copy(counter, counter + np, tdp_linearIndex);
  thrust::sort_by_key(tdp_sortIndex, tdp_sortIndex + np, tdp_linearIndex);

  /* apply the sort (copies particles to d_temp_particles) */
  thrust::gather(tdp_linearIndex, tdp_linearIndex + np, zipInit, zipFinal);
  gpuErrorHandler("sortByPID: thrust::gather");

  /* swap the device arrays */
  std::swap(gpuBase->d_particles, gpuBase->d_temp_particles);

  /* sort accepted particles to maintain corresponding orderings */
  if (gpuBase->d_accepted) {
    thrust::device_ptr<double> dap_0(d_accepted + 0*particlePitch);
    thrust::device_ptr<double> dap_1(d_accepted + 1*particlePitch);
    thrust::device_ptr<double> dap_2(d_accepted + 2*particlePitch);
    thrust::device_ptr<double> dap_3(d_accepted + 3*particlePitch);
    thrust::device_ptr<double> dap_4(d_accepted + 4*particlePitch);
    thrust::device_ptr<double> dap_5(d_accepted + 5*particlePitch);
    thrust::device_ptr<double> dap_6(d_accepted + 6*particlePitch);
    
    dp7_zip_iterator zipAccept =  make_zip_iterator(
           thrust::make_tuple(dap_0, dap_1, dap_2, dap_3, dap_4, dap_5, dap_6) );

    /* apply the sort (copies accepted to zipInit which is now the temp array) */
    thrust::gather(tdp_linearIndex, tdp_linearIndex + np, zipAccept, zipInit);
    gpuErrorHandler("sortByPID: thrust::gather (2)");

    /* swap the device arrays */
    std::swap(gpuBase->d_accepted, gpuBase->d_temp_particles);
  }

  /* mark particles as sorted */
  gpuBase->particlesUnsorted = 0;
}

#ifndef GPU_VERIFY
/**
 * Dummy function
 */
void startGpuTimer() {
}
#endif

} // extern "C"
