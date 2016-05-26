/**
 * Routines to verify the GPU implementations. 
 */

#ifdef GPU_VERIFY
#include <string>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <fstream>

#include <gpu_track.h>
#include <gpu_base.h>

#define WARN_TOL 1e-10

/**
 * Handle cuda errors
 * @param err CUDA error code
 * @param diagnosis location string for error
 */
inline void gpuDriverErrorHandler(cudaError_t err, std::string diagnosis){
  if(err != 0)
    std::cout << "Kernel " << diagnosis << " returned with error type: " 
              << cudaGetErrorString(err) << std::endl;
}

/**
 * Singleton class to store timers
 */
class timerType {
public:
  cudaEvent_t start, stop;
  float gputime, cputime;
  timerType() {
    cudaError_t err;
    err = cudaEventCreate(&this->start);
    gpuDriverErrorHandler(err,"at cudaEventCreate in timerType");
    err = cudaEventCreate(&this->stop);
    gpuDriverErrorHandler(err,"at cudaEventCreate in timerType");
  }
  virtual ~timerType() {
    cudaError_t err;
    err = cudaEventDestroy(this->start);
    gpuDriverErrorHandler(err,"at cudaEventDestroy in ~timerType");
    err = cudaEventDestroy(this->stop);
    gpuDriverErrorHandler(err,"at cudaEventDestroy in ~timerType");
  }
};

/**
 * Get singleton timer storage class instance
 */
timerType* getTimer() {
  static timerType *timer = NULL;
  if (timer==NULL) timer = new timerType();
  return timer;
}

/**
 * Start the GPU timer
 */
void startGpuTimer(){
  cudaError_t err;
  timerType* timer = getTimer();
  struct GPUBASE* gpuBase = getGpuBase();
  gpuBase->gpuFtnLevel++;
  /* If you see this error message it means a GPU function called a CPU 
     function instead of the GPU analogue. */
  if (gpuBase->gpuFtnLevel>1) {
    //double *segfault=NULL;
    //*segfault=1;
    std::cout << "CATASTROPHIC ERROR: recursive gpu timings detected!" 
      << std::endl;
    exit(-1);
  }
  const bool synchronize = true;
  if (synchronize) {
    gpuBase->copyForValidation=1;
    copyParticlesFromGpuToCpu(gpuBase->nOriginal); 
  }
  err = cudaEventRecord(timer->start, 0);
  gpuDriverErrorHandler(err,"at cudaEventRecord in startGpuTimer");
}

/**
 * Start the CPU timer (record the GPU time)
 */
void startCpuTimer(){
  cudaError_t err;
  timerType* timer = getTimer();
  struct GPUBASE* gpuBase = getGpuBase();
  err = cudaEventRecord(timer->stop, 0);
  gpuDriverErrorHandler(err,"at cudaEventRecord in startCpuTimer");
  err = cudaEventSynchronize(timer->stop);
  gpuDriverErrorHandler(err,"at cudaEventSynchronize in startCpuTimer");
  err = cudaEventElapsedTime(&timer->gputime, timer->start, timer->stop);
  gpuDriverErrorHandler(err,"at cudaEventElapsedTime in startCpuTimer");
  gpuBase->elementOnGpu=0; /* run routines on the CPU */
  err = cudaEventRecord(timer->start, 0);
  gpuDriverErrorHandler(err,"at cudaEventRecord in startCpuTimer");
}

/**
 * Stop the timer and record timings
 * @param name timed routine name
 */
void endTimer(char* name){
  cudaError_t err;
  timerType* timer = getTimer();
  struct GPUBASE* gpuBase = getGpuBase();
  err = cudaEventRecord(timer->stop, 0);
  gpuDriverErrorHandler(err,"at cudaEventRecord in endTimer");
  err = cudaEventSynchronize(timer->stop);
  gpuDriverErrorHandler(err,"at cudaEventSynchronize in endTimer");
  err = cudaEventElapsedTime(&timer->cputime, timer->start, timer->stop);
  gpuDriverErrorHandler(err,"at cudaEventElapsedTime in endTimer");
  recordTimings(name,timer->gputime,timer->cputime);
  gpuBase->elementOnGpu=1;
  gpuBase->gpuFtnLevel--;
}

/**
 * Singleton class to store function timings and linked list
 */
class timerListType {
public:
  std::string* name;
  int ncalls;
  float gpuTime;
  float cpuTime;
  timerListType* next;
  timerListType() {
    this->name = new std::string;
    *this->name = "none"; 
    this->ncalls = 0;
    this->gpuTime = 0; 
    this->cpuTime = 0;
    this->next = NULL;
  }
  virtual ~timerListType() {
    delete this->name;
  }
};

/**
 * Get the first instance of the singleton linked list for timings
 */
timerListType* getTimerList(bool reset=false) {
  static timerListType* timer_list = NULL;
  if (reset) {
    finalizeTimer();
    timer_list = NULL;
  }
  if (timer_list == NULL) timer_list = new timerListType;
  return timer_list;
}

/**
 * Deallocate the timer list
 */
void finalizeTimer(){
  timerListType* tlist;
  timerListType* timer_list = getTimerList();
  while (true) {
    if (timer_list->next != NULL) tlist=timer_list->next;
    else tlist=NULL;
    delete timer_list;
    if (tlist == NULL) return;
    else timer_list=tlist;
  }
}

/**
 * Record timings
 * @param name timed routine name
 * @param gputime gputime of routine
 * @param cputime cputime of routine
 */
void recordTimings(char* name, float gputime, float cputime){

  timerListType* timer_list = getTimerList();
  if (timer_list->name->compare("none")==0) *timer_list->name=name;

  timerListType* tlist=timer_list;
  while (tlist->name->compare(name)) {
    if (tlist->next==NULL) {
      tlist->next = new timerListType;
      tlist = tlist->next;
      *tlist->name = name; 
    } else {
      tlist=tlist->next;
    }
  }
 
  tlist->ncalls++;
  tlist->gpuTime+=gputime;
  tlist->cpuTime+=cputime;
}

/**
 * This class handles MPI file based output for validation
 */
class mpiStream {
private:
  std::fstream coss;
public:
#if USE_MPI
  mpiStream() {
    char filename[25];
    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    sprintf(filename, "gpuval.%i", rank);
    coss.open(filename, std::fstream::out | std::fstream::app);
  };
  virtual ~mpiStream() {
    if(coss.is_open()) coss.close();
  };
#else
  mpiStream() {}
  virtual ~mpiStream() {} 
#endif
  template <class T> mpiStream& operator<< (T val);
  mpiStream& operator<< (std::ostream& (*pfun)(std::ostream&));
};

/**
 * Wrapper for MPI file based writes of a value
 */
template <class T>
mpiStream& mpiStream::operator<< (T val){
#if USE_MPI
  this->coss << val;
#endif
  std::cout << val;
  return *this;
}

/**
 * Wrapper for MPI file based writes of a function pointer
 */
mpiStream& mpiStream::operator<< (std::ostream& (*pfun)(std::ostream&)){
#if USE_MPI
  pfun(this->coss);
#endif
  pfun(std::cout);
  return *this;
}

/**
 * Display the timings
 */
void displayTimings() {
  timerListType* timer_list = getTimerList();
  timerListType* tlist = timer_list;

#if USE_MPI
  int err=0;
  int rank, nprocs;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  if (rank==0 && nprocs>1) {
    /* Master node: collect timers from other processes and average */
    // 1) Determine number of function types called
    int minNFtns(0), maxNFtns;
    err += MPI_Allreduce(&minNFtns, &maxNFtns, 1, MPI_INT, MPI_MAX,
                         MPI_COMM_WORLD);
    err += MPI_Allreduce(&maxNFtns, &minNFtns, 1, MPI_INT, MPI_MIN,
                         MPI_COMM_WORLD);
    if (maxNFtns != minNFtns) {
      std::cout << "displayTimings: error: unbalanced gpu kernel calls."
                << std::endl;
      return;
    }

    // 2) Collect and average each function type
    int izero=0;
    float fzero=0;
    int csize;
    timer_list = getTimerList(true); // Reset my timer
    tlist = timer_list;
    for(int iFtn=0; iFtn < maxNFtns; iFtn++) {
      err += MPI_Recv(&csize, 1, MPI_INT, 1, 0, MPI_COMM_WORLD,
                      MPI_STATUS_IGNORE);
      char* tname = new char[csize];
      err += MPI_Recv(tname, csize, MPI_CHAR, 1,
                      0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      *tlist->name=tname;
      delete tname;
      err += MPI_Reduce(&izero, &tlist->ncalls, 1, MPI_INT, 
                        MPI_SUM, 0, MPI_COMM_WORLD);
      err += MPI_Reduce(&fzero, &tlist->gpuTime, 1, MPI_FLOAT,
                        MPI_SUM, 0, MPI_COMM_WORLD);
      err += MPI_Reduce(&fzero, &tlist->cpuTime, 1, MPI_FLOAT,
                        MPI_SUM, 0, MPI_COMM_WORLD);
      tlist->ncalls/=float(nprocs-1);
      tlist->gpuTime/=float(nprocs-1);
      tlist->cpuTime/=float(nprocs-1);
      if (iFtn < maxNFtns - 1) {
        tlist->next = new timerListType();
        tlist = tlist->next;
      }
    }

    if (err != 0) {
      std::cout << "displayTimings: MPI communication error." << std::endl;
      return;
    }

    // 3) Print statistics (below)
    tlist=timer_list;

  } else {
    /* Slave nodes: send timings to the master node*/
    // 1) Determine number of function types called
    int myNFtns(0), maxNFtns, minNFtns;
    while(tlist != NULL) {
      myNFtns++;
      tlist = tlist->next;
    }
    tlist = timer_list;
    err += MPI_Allreduce(&myNFtns, &maxNFtns, 1, MPI_INT, 
                         MPI_MAX, MPI_COMM_WORLD);
    err += MPI_Allreduce(&myNFtns, &minNFtns, 1, MPI_INT,
                         MPI_MIN, MPI_COMM_WORLD);
    if (maxNFtns != minNFtns) return;

    // 2) Collect and average each function type
    int csize;
    for(int iFtn=0; iFtn < maxNFtns; iFtn++) {
      if (rank==1) {
        csize = 2*( (int)tlist->name->length() + 2);
        err += MPI_Send(&csize, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        err += MPI_Send((void*)tlist->name->c_str(), csize, MPI_CHAR, 0, 0,
                        MPI_COMM_WORLD);
      }
      err += MPI_Reduce(&tlist->ncalls, NULL, 1, MPI_INT, 
                        MPI_SUM, 0, MPI_COMM_WORLD);
      err += MPI_Reduce(&tlist->gpuTime, NULL, 1, MPI_FLOAT, 
                        MPI_SUM, 0, MPI_COMM_WORLD);
      err += MPI_Reduce(&tlist->cpuTime, NULL, 1, MPI_FLOAT, 
                        MPI_SUM, 0, MPI_COMM_WORLD);
      tlist = tlist->next;
    }   

    // 3) Return
    return;
  }
#endif

  std::cout << "\nGPU vs CPU average timings (ms) \n"
            << "Routine                   #Calls   GPUtime      CPUtime      SpeedUp \n";
  while(tlist != NULL) {
    std::cout << std::setw(25) << std::left << tlist->name->c_str() << " "
              << std::setw(8)  << std::left << tlist->ncalls << " "
              << std::setw(12) << std::left << tlist->gpuTime/tlist->ncalls << " "
              << std::setw(12) << std::left << tlist->cpuTime/tlist->ncalls << " "
              << std::setw(8)  << std::left << tlist->cpuTime/tlist->gpuTime 
              << "\n";
    tlist=(timerListType*)tlist->next;
  }
  std::cout << std::endl;
}

/**
 * Compare two particles based on the last index.
 * @param p1 pointer to particle 1
 * @param p2 pointer to particle 2
 */
int indexCompare(const void* p1, const void* p2) {
#define COMPIND 6 /* use the last component */
  if (((double*)p1)[COMPIND] < ((double*)p2)[COMPIND]) return -1;
  if (((double*)p1)[COMPIND] > ((double*)p2)[COMPIND]) return  1;
  return 0;
}

/**
 * Compare GPU and CPU results and aggregate timings
 * @param n_part number of particles
 * @param name compared routine name
 */
void compareGpuCpu(unsigned int n_part, char* name){
  struct GPUBASE* gpuBase = getGpuBase();
  double** h_particles = gpuBase->coord;
  unsigned int n_comp = gpuBase->n_comp;

  endTimer(name);
  if(n_part==0) return;
  //std::cout << "Comparing reduction routine " << name << std::endl;

  /* Copy device particles to the host */
  gpuBase->doh_particles = (double**)czarray_2d(sizeof(**gpuBase->doh_particles),
                                               n_part, n_comp);
  double** doh_particles = gpuBase->doh_particles;
  copyParticlesFromGpuToCpu(n_part);

  /* sort based on the last index before comparison */
  std::qsort(  *h_particles,n_part,sizeof(double)*n_comp,indexCompare);
  std::qsort(*doh_particles,n_part,sizeof(double)*n_comp,indexCompare);

  /* compare */
  double residual, avgdiff(0), maxdiff(0);
  double const tol=WARN_TOL;
  int maxip(0), maxic(0), ndiff(0);
  double mindbl = DBL_MIN/tol*1.e3;
  bool gpuNaN = false;
  bool cpuNaN = false;
  for(unsigned int ip=0; ip<n_part; ip++) {
    for(unsigned int ic=0; ic<n_comp; ic++) {
      if (isnan(h_particles[ip][ic])) cpuNaN = true;
      if (isnan(doh_particles[ip][ic])) gpuNaN = true;
      if (std::fabs(h_particles[ip][ic])>mindbl) {
        residual=std::fabs((doh_particles[ip][ic]-h_particles[ip][ic])
                           /h_particles[ip][ic]);
        if (residual > tol) {
          ndiff++;
          avgdiff+=residual;
          if (residual > maxdiff) {
            maxip=ip; maxic=ic;
            maxdiff=residual;
          }
        }
      }
    }
  }

  mpiStream fcout;

  if (ndiff > 0) {
    avgdiff=avgdiff/(double)ndiff;
    fcout << "WARNING: GPU/CPU comparison failed for routine " << name
              << ", element " << gpuBase->elemName 
              << "\n # failed/total particle components: " 
              << ndiff << " / " << n_part*n_comp
              << "\n average relative difference: " << avgdiff
              << "\n maximum relative difference: " << maxdiff
              << " for particle " << maxip << " component " << maxic
              << "\n  host max:";
    for(unsigned int ic=0; ic<n_comp; ic++)
      fcout << std::fixed << std::setprecision(20) << " " << h_particles[maxip][ic];
    fcout << "\n  dev. max:";
    for(unsigned int ic=0; ic<n_comp; ic++)
      fcout << std::fixed << std::setprecision(20) << " " << doh_particles[maxip][ic];
    fcout << std::endl;
  }
  if (cpuNaN && gpuNaN) {
    fcout << "WARNING: Host and device NaNs found after routine " << name
          << ", element " << gpuBase->elemName << std::endl;
  } else if (cpuNaN) {
    fcout << "WARNING: Host NaNs found after routine " << name
          << ", element " << gpuBase->elemName << std::endl;
  } else if (gpuNaN) {
    fcout << "*WARNING*: Device NaNs found after routine " << name
          << ", element " << gpuBase->elemName << std::endl;
  }

  free_czarray_2d((void**)gpuBase->doh_particles, n_part, n_comp);
  gpuBase->doh_particles=NULL;
}

#define NCOMP 7
/**
 * Set pointers to stored reduction arrays computed by the CUDA code
 */
void setGpuReductionArrays(double** pgpu_centroid, double** pgpu_sigma, 
                           BEAM_SUMS** pgpu_sums) {
  static double gpu_centroid[NCOMP]={0.,0.,0.,0.,0.,0.,0.};
  static double gpu_sigma[NCOMP]={0.,0.,0.,0.,0.,0.,0.};
  static BEAM_SUMS gpu_sums;

  if (pgpu_centroid) *pgpu_centroid = gpu_centroid;
  if (pgpu_sigma) *pgpu_sigma = gpu_sigma;
  if (pgpu_sums) *pgpu_sums = &gpu_sums;
}



/**
 * Get the GPU beam sums array
 */
void* getGpuBeamSums(void* cpu_beam_sums) {
  BEAM_SUMS* gpu_sums;
  BEAM_SUMS* sums = (BEAM_SUMS*)cpu_beam_sums;
  setGpuReductionArrays(NULL, NULL, &gpu_sums);

  /* copy the beam sums */
  long i, j;
  for (i=0; i<7; i++) {
    gpu_sums->maxabs[i] = sums->maxabs[i];
    gpu_sums->max[i] = sums->max[i];
    gpu_sums->min[i] = sums->min[i];
  }
  for (i=0; i<7; i++) {
    gpu_sums->centroid[i] = sums->centroid[i];
    for (j=i; j<7; j++)
      gpu_sums->sigma[i][j] = sums->sigma[i][j];
  }
  gpu_sums->n_part = sums->n_part;
  gpu_sums->z      = sums->z;
  gpu_sums->p0     = sums->p0;

  return (void*)gpu_sums;
}

/**
 * Copy and store CUDA computed reduction arrays.
 */
void copyReductionArrays(double* centroid, double* sigma) {
  double* gpu_centroid, *gpu_sigma;
  setGpuReductionArrays(&gpu_centroid, &gpu_sigma, NULL);

  for (unsigned int ii=0; ii<NCOMP; ii++){
    if (centroid) gpu_centroid[ii] = centroid[ii];
    if (sigma)    gpu_sigma[ii]    = sigma[ii];
  }
}

/**
 * Compare two values and write output if not within tolerances
 */
inline void compareValues(double cpuval, double gpuval,
                          std::string name, int index=-1, int index2=-1) {
  const double tol = WARN_TOL;
  const double mindbl = DBL_MIN/tol*1.e3;
  double reldiff;

  if(std::fabs(cpuval) > mindbl) {
    reldiff = std::fabs((cpuval-gpuval)/cpuval);
    if (reldiff > tol) { 
      if (index2 > -1)
        std::cout << "WARNING: " << name << "[" << index  << "][" << index2
                  << "] differs: CPU= " << cpuval << " GPU= " 
                  << gpuval << " reldiff= " << reldiff << std::endl;
      else if (index > -1)
        std::cout << "WARNING: " << name << "[" << index 
                  << "] differs: CPU= " << cpuval << " GPU= " 
                  << gpuval << " reldiff= " << reldiff << std::endl;
      else
        std::cout << "WARNING: " << name << " differs: CPU= " 
                  << cpuval << " GPU= " << gpuval << " reldiff= " 
                  << reldiff << std::endl;
    }
  }
}

/**
 * Compare passed CPU arrays with the stored GPU results.
 */
void compareReductionArrays(double* centroid, double* sigma, void* vsums,
                            char* name) {
  BEAM_SUMS* sums = (BEAM_SUMS*)vsums;
  double* gpu_centroid, *gpu_sigma;
  BEAM_SUMS* gpu_sums;
  long ii, jj;

  endTimer(name);
//  std::cout << "Comparing reduction routine " << name << std::endl;

  setGpuReductionArrays(&gpu_centroid, &gpu_sigma, &gpu_sums);
  if (sums) {
    compareValues(sums->p0, gpu_sums->p0, "p0");
    for (ii=0; ii<7; ii++) {
      compareValues(sums->maxabs[ii], gpu_sums->maxabs[ii], "maxabs", ii);
      compareValues(sums->max[ii], gpu_sums->max[ii], "max", ii);
      compareValues(sums->min[ii], gpu_sums->min[ii], "min", ii);
      compareValues(sums->centroid[ii], gpu_sums->centroid[ii], "centroid", ii);
      for (jj=ii; jj<7; jj++) 
        compareValues(sums->sigma[ii][jj], gpu_sums->sigma[ii][jj],
                      "sigma", ii, jj);
    }
  }

  if (centroid) 
    for (ii=0; ii<NCOMP; ii++)
      compareValues(centroid[ii], gpu_centroid[ii], "-centroid", ii);

  if (sigma) 
    for (ii=0; ii<NCOMP; ii++)
      compareValues(sigma[ii], gpu_sigma[ii], "-sigma", ii);

}

#include "csbend.h"

void compareCSR_LAST_WAKE(void* vgpuCsrWake, void* vcpuCsrWake) {
  CSR_LAST_WAKE *gpuCsrWake = (CSR_LAST_WAKE*) vgpuCsrWake;
  CSR_LAST_WAKE *cpuCsrWake = (CSR_LAST_WAKE*) vcpuCsrWake;

  compareValues(gpuCsrWake->valid, cpuCsrWake->valid, "valid");
  compareValues(gpuCsrWake->bins, cpuCsrWake->bins, "bins");
  compareValues(gpuCsrWake->dctBin, cpuCsrWake->dctBin, "dctBin");
  compareValues(gpuCsrWake->s0, cpuCsrWake->s0, "s0");
  compareValues(gpuCsrWake->ds0, cpuCsrWake->ds0, "ds0");
  compareValues(gpuCsrWake->zLast, cpuCsrWake->zLast, "zLast");
  compareValues(gpuCsrWake->z0, cpuCsrWake->z0, "z0");
  compareValues(gpuCsrWake->S11, cpuCsrWake->S11, "S11");
  compareValues(gpuCsrWake->S12, cpuCsrWake->S12, "S12");
  compareValues(gpuCsrWake->S22, cpuCsrWake->S22, "S22");
  compareValues(gpuCsrWake->nSaldin, cpuCsrWake->nSaldin, "nSaldin");
  if (gpuCsrWake->nSaldin == cpuCsrWake->nSaldin) {
    for (int ii = 0; ii < gpuCsrWake->nSaldin; ii++) {
      compareValues(gpuCsrWake->FdNorm[ii], cpuCsrWake->FdNorm[ii],
          "FdNorm", ii);
      compareValues(gpuCsrWake->xSaldin[ii], cpuCsrWake->xSaldin[ii],
          "xSaldin", ii);
    }
  }
  compareValues(gpuCsrWake->rho, cpuCsrWake->rho, "rho");
  compareValues(gpuCsrWake->bendingAngle, cpuCsrWake->bendingAngle,
      "bendingAngle");
  compareValues(gpuCsrWake->Po, cpuCsrWake->Po, "Po");
  compareValues(gpuCsrWake->SGOrder, cpuCsrWake->SGOrder, "SGOrder");
  compareValues(gpuCsrWake->SGDerivOrder, cpuCsrWake->SGDerivOrder,
      "SGDerivOrder");
  compareValues(gpuCsrWake->SGHalfWidth, cpuCsrWake->SGHalfWidth,
      "SGHalfWidth");
  compareValues(gpuCsrWake->SGDerivHalfWidth, cpuCsrWake->SGDerivHalfWidth,
      "SGDerivHalfWidth");
  compareValues(gpuCsrWake->GSConstant, cpuCsrWake->GSConstant, "GSConstant");
  compareValues(gpuCsrWake->MPCharge, cpuCsrWake->MPCharge, "MPCharge");
  compareValues(gpuCsrWake->binRangeFactor, cpuCsrWake->binRangeFactor,
      "binRangeFactor");
  compareValues(gpuCsrWake->trapazoidIntegration,
      cpuCsrWake->trapazoidIntegration, "trapazoidIntegration");
}

#endif /* GPU_VERIFY */
