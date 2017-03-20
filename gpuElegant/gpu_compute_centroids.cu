#include <gpu_track.h>

#include <gpu_base.h>
#include <gpu_particle_template_function.hcu>
#include <gpu_particle_reduction.hcu>
#include <gpu_reductions.hcu>
#include <gpu_reductions.h>
#include <gpu_kahan.hcu>
#include <gpu_kahan.h>
#include <gpu_trwake.h> // gpu_computeTimeCoordinates

extern "C" {

void gpu_compute_centroids(double *centroid, long n_part)
{
  struct GPUBASE* gpubase = getGpuBase();
  double* d_particles = gpubase->d_particles;
  unsigned int particlePitch = gpubase->gpu_array_pitch;
  long i_coord;
  double sum[6];
  long active = 1; 
#ifdef USE_KAHAN
  double error[6];
#endif

#if USE_MPI  /* In the non-parallel mode, it will be same with the serial version */ 
  long n_total;
#ifdef USE_KAHAN
  long j;
  double error_sum=0.0, error_total=0.0,
    **sumMatrix, **errorMatrix,
    *sumArray, *errorArray;
  sumMatrix = (double**)czarray_2d(sizeof(**sumMatrix), n_processors, 6); 
  errorMatrix = (double**)czarray_2d(sizeof(**errorMatrix), n_processors, 6);
  sumArray = (double*)malloc(sizeof(double) * n_processors);
  errorArray = (double*)malloc(sizeof(double) * n_processors);
#endif
  if (notSinglePart) {
    if (((parallelStatus==trueParallel) && isSlave) || ((parallelStatus!=trueParallel) && isMaster))
      active = 1;
    else 
      active = 0;
  }
#endif
 
  for (i_coord=0; i_coord<6; i_coord++) {
    sum[i_coord] = centroid[i_coord] = 0;
#ifdef USE_KAHAN
    error[i_coord] = 0.0;
#endif
  }

  if (active) {
    for (i_coord=0; i_coord<6; i_coord++)
#ifndef USE_KAHAN
      gpuReduceAddAsync(d_particles+particlePitch*i_coord, n_part, &sum[i_coord]);
#else
      gpuKahanAsync(d_particles+particlePitch*i_coord, &sum[i_coord],
                    &error[i_coord], n_part);
#endif
  }
  finishReductionStreams();

  if (!USE_MPI || !notSinglePart) {
    if (n_part)
      for (i_coord=0; i_coord<6; i_coord++)
	centroid[i_coord] = sum[i_coord]/n_part;

  }
#if USE_MPI
  if (notSinglePart) {
    if (parallelStatus!=trueParallel) {
      if (isMaster && n_part)
	for (i_coord=0; i_coord<6; i_coord++)
	  centroid[i_coord] = sum[i_coord]/n_part;
    }
    else {
      if (isMaster) {
	n_part = 0;
      }
      /* compute centroid sum over processors */
#ifndef USE_KAHAN  
      MPI_Allreduce(sum,centroid,6,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#else
      MPI_Allgather(sum,6,MPI_DOUBLE,&sumMatrix[0][0],6,MPI_DOUBLE,MPI_COMM_WORLD);
      /* compute error sum over processors */
      MPI_Allgather(error,6,MPI_DOUBLE,&errorMatrix[0][0],6,MPI_DOUBLE,MPI_COMM_WORLD);

      for (i_coord=0; i_coord<6; i_coord++) {
        error_sum = 0.0;
	/* extract the columnwise array from the matrix */
	for (j=0; j<n_processors; j++) {         
	  sumArray[j] = sumMatrix[j][i_coord];             
	  errorArray[j] = errorMatrix[j][i_coord];
	}
	centroid[i_coord] = Kahan(n_processors-1,&sumArray[1],&error_sum);
	error_total = Kahan(n_processors-1,&errorArray[1],&error_sum);
	centroid[i_coord] += error_total; 
      }
#endif 
      /* compute total number of particles over processors */
      MPI_Allreduce(&n_part, &n_total, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);  
      if (n_total)
	for (i_coord=0; i_coord<6; i_coord++)
	  centroid[i_coord] /= n_total;

    }
  }
#endif  
   
#if USE_MPI  /* In the non-parallel mode, it will be same with the serial version */ 
#ifdef USE_KAHAN
  free_czarray_2d((void**)sumMatrix, n_processors, 6); 
  free_czarray_2d((void**)errorMatrix, n_processors, 6);
  free(sumArray);
  free(errorArray);
#endif
#endif

}

} // extern "C"

class gpuPZ {
public:
  gpuPZ(double* d_pz, double p_central){
    this->d_pz=d_pz;
    this->p_central=p_central;
  }
  __device__ void inline operator()(gpuParticleAccessor& particle){
    unsigned int ind = particle.getParticleIndex();
    d_pz[ind] = p_central*(1+particle[5])
                /sqrt(1 + particle[1]*particle[1] + particle[3]*particle[3]);
  }
  double* d_pz;
  double p_central;
};

class gpuCenter {
public:
  gpuCenter(double* d_center, double* d_input, double center){
    this->d_center=d_center;
    this->d_input=d_input;
    this->center=center;
  }
  __device__ void inline operator()(gpuParticleAccessor& particle){
    unsigned int ind = particle.getParticleIndex();
    d_center[ind]=d_input[ind]-center;
  }
  double* d_center;
  double* d_input;
  double center;
};

class gpuPZmult {
public:
  gpuPZmult(double* d_pz, int ii){
    this->d_pz=d_pz;
    this->ii=ii;
  }
  __device__ double inline operator()(gpuParticleAccessor& particle){
    unsigned int ind = particle.getParticleIndex();
    return particle[ii]*d_pz[ind];
  }
  double* d_pz;
  int ii;
};

class gpuSij {
public:
  gpuSij(double* d_arr1, double center1, double* d_arr2, double center2) {
    this->d_arr1=d_arr1;
    this->center1=center1;
    this->d_arr2=d_arr2;
    this->center2=center2;
  }
  __device__ double inline operator()(gpuParticleAccessor& particle){
    unsigned int ind = particle.getParticleIndex();
    return (d_arr1[ind]-center1)*(d_arr2[ind]-center2);
  }
  double *d_arr1, *d_arr2;
  double center1, center2;
};

class gpuSijn {
public:
  gpuSijn(int ii, double center1, int jj, double center2, double* d_pz) {
    this->ii=ii;
    this->center1=center1;
    this->jj=jj;
    this->center2=center2;
    this->d_pz=d_pz;
  }
  __device__ double inline operator()(gpuParticleAccessor& particle){
    unsigned int ind = particle.getParticleIndex();
    double Sijn;
    if (ii==1 || ii==3)
      if (jj==1 || jj==3)
        Sijn = (particle[ii]*d_pz[ind]-center2)*(particle[jj]*d_pz[ind]-center2);
      else
        Sijn = (particle[ii]*d_pz[ind]-center2)*(particle[jj]-center2);
    else
      if (jj==1 || jj==3)
        Sijn = (particle[ii]-center2)*(particle[jj]*d_pz[ind]-center2);
      else
        Sijn = (particle[ii]-center2)*(particle[jj]-center2);
    return Sijn;
  }
  int ii, jj;
  double *d_pz;
  double center1, center2;
};

class gpuChosen {
public:
  gpuChosen(double startPID, double endPID, double* d_stPID) : 
    startPID(startPID), endPID(endPID), d_stPID(d_stPID) {}
  __device__ int inline operator()(gpuParticleAccessor& particle){
    unsigned int ind = particle.getParticleIndex();
    /* use to mark the startPID */
    if (particle[6]==startPID)
      d_stPID[ind] = -1;
    else
      d_stPID[ind] = 1;
    /* return 1 if chosen */
    if (particle[6]>=startPID && particle[6]<=endPID) 
      return 1;
    return 0;
  }
  double startPID, endPID;
  double* d_stPID;
};

extern "C" {

#if !USE_MPI
void gpu_accumulate_beam_sums(
                          BEAM_SUMS *sums,
                          long n_part,
                          double p_central,
			  double mp_charge,
                          double *timeValue, double tMin, double tMax,
                          long startPID, long endPID,
                          unsigned long flags
                          )
{
  long i, j;
  double centroid[7], centroidn[7], pmin[7], pmax[7];
  double Sijarr[7][7], Sijnarr[7][7];
  double Sij, Sijn;
  long npCount=0;
  short sparse[7][7] = {
    { 1, 1, 0, 0, 0, 1, 0 },
    { 0, 1, 0, 0, 0, 1, 0 },
    { 0, 0, 1, 1, 0, 1, 0 },
    { 0, 0, 0, 1, 0, 1, 0 },
    { 0, 0, 0, 0, 1, 1, 1 },
    { 0, 0, 0, 0, 0, 1, 0 },
    { 0, 0, 0, 0, 0, 0, 1 }
  };
  
#ifdef USE_KAHAN
  double errorCen[7], errorCenn[7], errorSig, errorSign;
#endif

  if (startPID<endPID && n_part) {
    sortByPID(n_part);
  }
  /* set gpu pointers after the sort */
  struct GPUBASE* gpubase = getGpuBase();
  double* d_particles = gpubase->d_particles;
  unsigned int particlePitch = gpubase->gpu_array_pitch;
  double* d_temp_particles = gpubase->d_temp_particles;
  double* d_timeCoord = d_temp_particles;
  double* d_center  = d_temp_particles + 1*particlePitch;
  double* d_center4 = d_temp_particles + 2*particlePitch;
  double* d_pz      = d_temp_particles + 3*particlePitch;

  if (timeValue != NULL) {
    fprintf(stderr, "timeValue in gpu_accumulate_beam_sums is not implemented yet\n");
    exit;
  }
  gpu_computeTimeCoordinates(n_part, d_timeCoord, p_central);
  
  
  if (exactNormalizedEmittance) {
    gpuDriver(n_part, gpuPZ(d_pz, p_central));
  }

  if (!sums->n_part)
    sums->p0 = p_central;

  if (n_part) {
    if (startPID<endPID && n_part) {
      npCount =
        gpuParticleReduction(n_part, gpuChosen(startPID, endPID, d_center),
                             Add<double>());
      /* adjust to subset of particles */
      double minVal;
      unsigned int partOffset;
      gpuFindMinIndex(d_center, n_part, &minVal, &partOffset);
      selectParticleSubset(partOffset);
      d_particles = d_particles + partOffset;
      d_timeCoord = d_timeCoord + partOffset;
      d_pz = d_pz + partOffset;
    } else {
      npCount = n_part;
    }
    if (!(flags&BEAM_SUMS_NOMINMAX)) {
      /* maximum amplitudes */
      for (i=0; i<6; i++) {
        if (i==4)
          continue;  /* done below */
        gpuReduceMinMaxAsync(d_particles+particlePitch*i, npCount, &pmin[i],
                             &pmax[i]);
      }
    }
    
    /* compute centroids for present beam and add in to existing centroid data */
    for (i=0; i<7; i++) {
#ifdef USE_KAHAN
      errorCen[i] = errorCenn[i] = 0.0;
#endif
      centroid[i]=0;
      centroidn[i]=0;
#ifndef USE_KAHAN
      if (exactNormalizedEmittance && (i==1 || i==3))
        gpuParticleReductionAsync(npCount, &centroidn[i],
                                  gpuPZmult(d_pz, i), Add<double>());
      if (i<6) gpuReduceAddAsync(d_particles+particlePitch*i, npCount,
                                 &centroid[i]);
      else     gpuReduceAddAsync(d_timeCoord, npCount, &centroid[i]);
#else
      if (exactNormalizedEmittance && (i==1 || i==3))
        gpuParticleKahanReductionAsync(npCount, &centroidn[i], &errorCenn[i],
                                       gpuPZmult(d_pz, i));
      if (i<6) gpuKahanAsync(d_particles+particlePitch*i, &centroid[i],
                             &errorCenn[i], npCount);
      else     gpuKahanAsync(d_timeCoord, &centroid[i], &errorCen[i], npCount);
#endif
    }
    finishReductionStreams();
    for (i=0; i<7; i++) {
      if (npCount) {	
        sums->centroid[i] = (sums->centroid[i]*sums->n_part+centroid[i])/(sums->n_part+npCount);
        centroid[i] /= npCount;
        centroidn[i] /= npCount;
      }
    }
    for (i=0; i<6; i++) {
      if (i==4)
        continue;  /* done below */
      if(pmax[i] > -pmin[i]) sums->maxabs[i] = pmax[i];
      else sums->maxabs[i] = -pmin[i];
      if(pmax[i] > sums->max[i]) sums->max[i] = pmax[i];
      if(pmin[i] < sums->min[i]) sums->min[i] = pmin[i];
      if (i!=1 && i!=3)
        centroidn[i] = centroid[i];
    }

    
    if (!(flags&BEAM_SUMS_NOMINMAX)) {
      i = 4;
      gpuDriver(npCount, gpuCenter(d_center4, d_particles+particlePitch*i, centroid[i]));
      gpuReduceMinMaxAsync(d_center4, npCount, &pmin[i], &pmax[i]);
      i = 6;  /* time coordinate */
      gpuDriver(npCount, gpuCenter(d_center, d_timeCoord, centroid[i]));
      gpuReduceMinMaxAsync(d_center, npCount, &pmin[i], &pmax[i]);
      finishReductionStreams();
      for (i = 4; i<=6; i+=2) {
        if(pmax[i] > -pmin[i]) sums->maxabs[i] = pmax[i];
        else sums->maxabs[i] = -pmin[i];
        if(pmax[i] > sums->max[i]) sums->max[i] = pmax[i];
        if(pmin[i] < sums->min[i]) sums->min[i] = pmin[i];
      }
    }
    
    /* compute Sigma[i][j] for present beam and add to existing data */
    for (i=0; i<7; i++) {
      for (j=i; j<7; j++) {
#ifdef USE_KAHAN
        errorSig=errorSign=0.0;
#endif
        Sij = Sijn = 0;
        if (flags&BEAM_SUMS_SPARSE && !sparse[i][j]) {
          /* Only compute the diagonal blocks and dispersive correlations */
          sums->sigma[j][i] = sums->sigma[i][j] = 0;
          sums->sigman[j][i] = sums->sigman[i][j] = 0;
          continue;
        }
#ifndef USE_KAHAN
        if (i<6 && j<6)
          gpuParticleReductionAsync(npCount,
              &Sijarr[i][j], gpuSij(d_particles + particlePitch*i,
              centroid[i], d_particles + particlePitch*j, centroid[j]),
              Add<double>());
        else if (i<6)
          gpuParticleReductionAsync(npCount,
              &Sijarr[i][j], gpuSij(d_particles + particlePitch*i,
                centroid[i], d_center, 0.), Add<double>());
        else if (j<6)
          gpuParticleReductionAsync(npCount,
              &Sijarr[i][j], gpuSij(d_center, 0.,
              d_particles + particlePitch*j, centroid[j]), Add<double>());
        else
          gpuParticleReductionAsync(npCount,
              &Sijarr[i][j], gpuSij(d_center, 0., d_center, 0.),
              Add<double>());
        if (exactNormalizedEmittance && i<4 && j<4)
          gpuParticleReductionAsync(npCount, &Sijnarr[i][j],
              gpuSijn(i, centroidn[i], j, centroidn[j], d_pz));
#else
        if (i<6 && j<6)
          gpuParticleKahanReductionAsync(npCount,
              &Sijarr[i][j], &errorSig,
              gpuSij(d_particles + particlePitch*i, centroid[i],
              d_particles + particlePitch*j, centroid[j]));
        else if (i<6)
          gpuParticleKahanReductionAsync(npCount,
              &Sijarr[i][j], &errorSig,
              gpuSij(d_particles + particlePitch*i, centroid[i],
              d_center, 0.));
        else if (j<6)
          gpuParticleKahanReductionAsync(npCount,
              &Sijarr[i][j], &errorSig, gpuSij(d_center, 0.,
              d_particles + particlePitch*j, centroid[i]));
        else
          gpuParticleKahanReductionAsync(npCount,
              &Sijarr[i][j], &errorSig,
              gpuSij(d_center, 0., d_center, 0.));
        if (exactNormalizedEmittance && i<4 && j<4)
          gpuParticleKahanReductionAsync(npCount,
              &Sijnarr[i][j], &errorSign,
              gpuSijn(i, centroidn[i], j, centroidn[j], d_pz));
#endif
      }
    }
    finishReductionStreams();
    for (i=0; i<7; i++) {
      for (j=i; j<7; j++) {
        sums->sigma[j][i] = (sums->sigma[i][j] = (sums->sigma[i][j]*sums->n_part+Sijarr[i][j])/(sums->n_part+npCount));
        if (exactNormalizedEmittance) {
          sums->sigman[j][i] = (sums->sigman[i][j] = (sums->sigman[i][j]*sums->n_part+Sijnarr[i][j])/(sums->n_part+npCount));
        }
      }
    }
  }
  
  sums->charge = mp_charge*npCount;
  sums->n_part += npCount;

  resetToFullParticleArray();
}

#else
/* USE_MPI=1 */
/*TODO: This procedure probably doesn't work because it was not tested */
void gpu_accumulate_beam_sums(
                          BEAM_SUMS *sums,
                          long n_part,
                          double p_central, 
			  double mp_charge,
                          double *timeValue, double tMin, double tMax,
                          long startPID, long endPID,
                          unsigned long flags
                          )
{
  long i, j;
  double centroid[7], pmin[7], pmax[7];
  double Sijarr[7][7], Sijnarr[7][7];
  long active = 1;
  long npCount=0, npCount_total = 0;
  short sparse[7][7] = {
    { 1, 1, 0, 0, 0, 1, 0 },
    { 0, 1, 0, 0, 0, 1, 0 },
    { 0, 0, 1, 1, 0, 1, 0 },
    { 0, 0, 0, 1, 0, 1, 0 },
    { 0, 0, 0, 0, 1, 1, 1 },
    { 0, 0, 0, 0, 0, 1, 0 },
    { 0, 0, 0, 0, 0, 0, 1 }
  };
#ifdef USE_KAHAN
  double errorCen[7], errorSig[28];
#endif

  double buffer[7], Sij_p[28], Sijn_p[28], Sij_total[28];
  long offset=0, index;  
#ifdef USE_KAHAN
  double error_sum=0.0, error_total=0.0,
    **sumMatrixCen, **errorMatrixCen,
    **sumMatrixSig, **errorMatrixSig,
    *sumArray, *errorArray;
  long k;
  sumMatrixCen = (double**)czarray_2d(sizeof(**sumMatrixCen), n_processors, 7); 
  errorMatrixCen = (double**)czarray_2d(sizeof(**errorMatrixCen), n_processors, 7);
  sumMatrixSig = (double**)czarray_2d(sizeof(**sumMatrixSig), n_processors, 28);
  errorMatrixSig = (double**)czarray_2d(sizeof(**errorMatrixSig), n_processors, 28);
  sumArray = (double*)malloc(sizeof(double) * n_processors);
  errorArray = (double*)malloc(sizeof(double) * n_processors);
#endif

  if (notSinglePart) {
    if (((parallelStatus==trueParallel) && isSlave) || ((parallelStatus!=trueParallel) && isMaster))
      active = 1;
    else 
      active = 0;
  }

  if (startPID<endPID && n_part) {
    sortByPID(n_part);
  }
  /* set gpu pointers after the sort */
  struct GPUBASE* gpubase = getGpuBase();
  double* d_particles = gpubase->d_particles;
  unsigned int particlePitch = gpubase->gpu_array_pitch;
  double* d_temp_particles = gpubase->d_temp_particles;
  double* d_timeCoord = d_temp_particles;
  double* d_center  = d_temp_particles + 1*particlePitch;
  double* d_center4 = d_temp_particles + 2*particlePitch;
  double* d_pz      = d_temp_particles + 3*particlePitch;

  if (timeValue != NULL) {
    fprintf(stderr, "timeValue in gpu_accumulate_beam_sums is not implemented yet\n");
    exit;
  }
  gpu_computeTimeCoordinates(n_part, d_timeCoord, p_central);

  if (!sums->n_part) 
    sums->p0 = p_central;
  
  if (exactNormalizedEmittance) {
    gpuDriver(n_part, gpuPZ(d_pz, p_central));
  }

  if (startPID<endPID && n_part) {
    npCount =
      gpuParticleReduction(n_part, gpuChosen(startPID, endPID, d_center),
                           Add<double>());
    /* adjust to subset of particles */
    double minVal;
    unsigned int partOffset;
    gpuFindMinIndex(d_center, n_part, &minVal, &partOffset);
    selectParticleSubset(partOffset);
    d_particles = d_particles + partOffset;
    d_timeCoord = d_timeCoord + partOffset;
    d_pz = d_pz + partOffset;
  } else {
    npCount = n_part;
  }
  if (active) {
    if (!(flags&BEAM_SUMS_NOMINMAX) && !(flags&BEAM_SUMS_EXACTEMIT)) {
      /* maximum amplitudes */
      for (i=0; i<6; i++) {
        if (i==4)
          continue;  /* done below */
        gpuReduceMinMaxAsync(d_particles+particlePitch*i, npCount, &pmin[i],
                             &pmax[i]);
      }
    }
    /* compute centroids for present beam and add in to existing centroid data */
    for (i=0; i<7; i++) {
#ifdef USE_KAHAN
      errorCen[i] = 0.0;
#endif
      centroid[i]=0;
#ifndef USE_KAHAN
      if (flags&BEAM_SUMS_EXACTEMIT && (i==1 || i==3))
        gpuParticleReductionAsync(npCount, &centroid[i],
                                  gpuPZmult(d_pz, i), Add<double>());
      if (i<6) gpuReduceAddAsync(d_particles+particlePitch*i, npCount,
                                 &centroid[i]);
      else     gpuReduceAddAsync(d_timeCoord, npCount, &centroid[i]);
#else
      if (flags&BEAM_SUMS_EXACTEMIT && (i==1 || i==3))
        gpuParticleKahanReductionAsync(npCount, &centroid[i], &errorCenn[i],
                                       gpuPZmult(d_pz, i));
      if (i<6) gpuKahanAsync(d_particles+particlePitch*i, &centroid[i],
                             &errorCenn[i], npCount);
      else     gpuKahanAsync(d_timeCoord, &centroid[i], &errorCen[i], npCount);
#endif
    }
    finishReductionStreams();
    for (i=0; i<7; i++) {
      if (!notSinglePart && npCount) {	
        /* single-particle mode and this process has a particle */
        if (!(flags&BEAM_SUMS_EXACTEMIT))
          sums->centroid[i] = (sums->centroid[i]*sums->n_part+centroid[i])/(sums->n_part+npCount);
        centroid[i] /= npCount;
      } 
      if (notSinglePart && (parallelStatus!=trueParallel) && isMaster) {
        /* multi-particle mode and, but not true parallel mode, so master has the particles */
        if (npCount) {
          if (!(flags&BEAM_SUMS_EXACTEMIT))
            sums->centroid[i] = (sums->centroid[i]*sums->n_part+centroid[i])/(sums->n_part+npCount);
          centroid[i] /= npCount;
        }
      }
    }
    for (i=0; i<6; i++) {
      if (i==4)
        continue;  /* done below */
      if(pmax[i] > -pmin[i]) sums->maxabs[i] = pmax[i];
      else sums->maxabs[i] = -pmin[i];
      if(pmax[i] > sums->max[i]) sums->max[i] = pmax[i];
      if(pmin[i] < sums->min[i]) sums->min[i] = pmin[i];
    }
  }
  
  if (notSinglePart) {
    if (parallelStatus==trueParallel) {
      if (isMaster) { 
        n_part = 0;  /* All the particles have been distributed to the slave processors */
        npCount = 0;
        memset(centroid, 0.0,  sizeof(double)*6);
      }
      /* compute centroid sum over processors */
#ifndef USE_KAHAN 
      MPI_Allreduce(centroid, buffer, 7, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      memcpy(centroid, buffer, sizeof(double)*7);
#else
      MPI_Allgather(centroid,7,MPI_DOUBLE,&sumMatrixCen[0][0],7,MPI_DOUBLE,MPI_COMM_WORLD);
      /* compute error sum over processors */
      MPI_Allgather(errorCen,7,MPI_DOUBLE,&errorMatrixCen[0][0],7,MPI_DOUBLE,MPI_COMM_WORLD); 
      for (i=0; i<7; i++) {
        error_sum = 0.0;
        centroid[i] = 0.0;
        /* extract the columnwise array from the matrix */
        for (j=1; j<n_processors; j++) {
          centroid[i] = KahanPlus(centroid[i], sumMatrixCen[j][i], &error_sum);
          centroid[i] = KahanPlus(centroid[i], errorMatrixCen[j][i], &error_sum);
        }
      }
#endif 
      /* compute total number of particles over processors */
      MPI_Allreduce(&npCount, &npCount_total, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
      if (npCount_total)
        for (i=0; i<7; i++) {
          if (!(flags&BEAM_SUMS_EXACTEMIT))
            sums->centroid[i] = (sums->centroid[i]*sums->n_part+centroid[i])/(sums->n_part+npCount_total);
          centroid[i] /= npCount_total;
        }     
    }
    else
      if (isMaster) {
        npCount_total = npCount;
      }
  }
  

  if (active) {	        
    if (!(flags&BEAM_SUMS_NOMINMAX) && !(flags&BEAM_SUMS_EXACTEMIT)) {
      i = 4;
      gpuDriver(npCount, gpuCenter(d_center4, d_particles+particlePitch*i, centroid[i]));
      gpuReduceMinMaxAsync(d_center4, npCount, &pmin[i], &pmax[i]);
      i = 6;  /* time coordinate */
      gpuDriver(npCount, gpuCenter(d_center, d_timeCoord, centroid[i]));
      gpuReduceMinMaxAsync(d_center, npCount, &pmin[i], &pmax[i]);
      finishReductionStreams();
      for (i = 4; i<=6; i+=2) {
        if(pmax[i] > -pmin[i]) sums->maxabs[i] = pmax[i];
        else sums->maxabs[i] = -pmin[i];
        if(pmax[i] > sums->max[i]) sums->max[i] = pmax[i];
        if(pmin[i] < sums->min[i]) sums->min[i] = pmin[i];
      }
    }
  }
  
  if (!(flags&BEAM_SUMS_NOMINMAX) && !(flags&BEAM_SUMS_EXACTEMIT)) {
    if (notSinglePart) {
      if (parallelStatus==trueParallel) {
        /* compute sums->maxabs over processors*/
        MPI_Allreduce(sums->maxabs, buffer, 7, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD); 
        memcpy(sums->maxabs, buffer, sizeof(double)*7);     
        /* compute sums->max over processors */
        MPI_Allreduce(sums->max, buffer, 7, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD); 
        memcpy(sums->max, buffer, sizeof(double)*7);      
        /* compute sums->min over processors */
        MPI_Allreduce(sums->min, buffer, 7, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD); 
        memcpy(sums->min, buffer, sizeof(double)*7);  
      }    
    }
  }

  /* compute Sigma[i][j] for present beam and add to existing data */
  if (active) {
    for (i=0; i<7; i++) {

      if ((parallelStatus==trueParallel) && notSinglePart)
        if (i>=1)        
          offset += i-1;

      for (j=i; j<7; j++) {
#ifdef USE_KAHAN
        errorSig[j]=0.0;
#endif
        if (flags&BEAM_SUMS_SPARSE && !sparse[i][j]) {
          /* Only compute the diagonal blocks and dispersive correlations */
          if (flags&BEAM_SUMS_EXACTEMIT)
            sums->sigman[j][i] = sums->sigman[i][j] = 0;
          else
            sums->sigma[j][i] = sums->sigma[i][j] = 0;
          continue;
        }

        if (notSinglePart) {
          if (parallelStatus==trueParallel) {
            index = 6*i+j-offset;
            Sij_p[index] = Sijn_p[index] = 0;
#ifdef USE_KAHAN
            errorSig[index] = 0.0;
#endif
#ifndef USE_KAHAN
            if (i<6 && j<6)
              gpuParticleReductionAsync(npCount,
                  &Sij_p[index], gpuSij(d_particles + particlePitch*i,
                  centroid[i], d_particles + particlePitch*j, centroid[j]),
                  Add<double>());
            else if (i<6)
              gpuParticleReductionAsync(npCount,
                  &Sij_p[index], gpuSij(d_particles + particlePitch*i,
                    centroid[i], d_center, 0.), Add<double>());
            else if (j<6)
              gpuParticleReductionAsync(npCount,
                  &Sij_p[index], gpuSij(d_center, 0.,
                  d_particles + particlePitch*j, centroid[j]), Add<double>());
            else
              gpuParticleReductionAsync(npCount,
                  &Sij_p[index], gpuSij(d_center, 0., d_center, 0.),
                  Add<double>());
	    if (flags&BEAM_SUMS_EXACTEMIT && i<4 && j<4)
	      gpuParticleReductionAsync(npCount, &Sijn_p[index],
		  gpuSijn(i, centroidn[i], j, centroidn[j], d_pz));
#else
            if (i<6 && j<6)
              gpuParticleKahanReductionAsync(npCount,
                  &Sij_p[index], &errorSig[j],
                  gpuSij(d_particles + particlePitch*i, centroid[i],
                  d_particles + particlePitch*j, centroid[j]));
            else if (i<6)
              gpuParticleKahanReductionAsync(npCount,
                  &Sij_p[index], &errorSig[j],
                  gpuSij(d_particles + particlePitch*i, centroid[i],
                  d_center, 0.));
            else if (j<6)
              gpuParticleKahanReductionAsync(npCount,
                  &Sij_p[index], &errorSig[j], gpuSij(d_center, 0.,
                  d_particles + particlePitch*j, centroid[i]));
            else
              gpuParticleKahanReductionAsync(npCount,
                  &Sij_p[index], &errorSig[j],
                  gpuSij(d_center, 0., d_center, 0.));
	    if (flags&BEAM_SUMS_EXACTEMIT && i<4 && j<4)
	      gpuParticleKahanReductionAsync(npCount,
		  &Sijn_p[index], &errorSign,
                  gpuSijn(i, centroidn[i], j, centroidn[j], d_pz));
#endif
          }
          else if (isMaster) {
#ifndef USE_KAHAN
            if (i<6 && j<6)
              gpuParticleReductionAsync(npCount,
                  &Sijarr[i][j], gpuSij(d_particles + particlePitch*i,
                  centroid[i], d_particles + particlePitch*j, centroid[j]),
                  Add<double>());
            else if (i<6)
              gpuParticleReductionAsync(npCount,
                  &Sijarr[i][j], gpuSij(d_particles + particlePitch*i,
                  centroid[i], d_center, 0.), Add<double>());
            else if (j<6)
              gpuParticleReductionAsync(npCount,
                  &Sijarr[i][j], gpuSij(d_center, 0.,
                  d_particles + particlePitch*j, centroid[j]), Add<double>());
            else
              gpuParticleReductionAsync(npCount,
                  &Sijarr[i][j], gpuSij(d_center, 0., d_center, 0.),
                  Add<double>());
	    if (flags&BEAM_SUMS_EXACTEMIT && i<4 && j<4)
	      gpuParticleReductionAsync(npCount, &Sijnarr[i][j],
		  gpuSijn(i, centroidn[i], j, centroidn[j], d_pz));
#else
            if (i<6 && j<6)
              gpuParticleKahanReductionAsync(npCount,
                  &Sijarr[i][j], &errorSig[j],
                  gpuSij(d_particles + particlePitch*i, centroid[i],
                  d_particles + particlePitch*j, centroid[j]));
            else if (i<6)
              gpuParticleKahanReductionAsync(npCount,
                  &Sijarr[i][j], &errorSig[j],
                  gpuSij(d_particles + particlePitch*i, centroid[i],
                  d_center, 0.));
            else if (j<6)
              gpuParticleKahanReductionAsync(npCount,
                  &Sijarr[i][j], &errorSig[j], gpuSij(d_center, 0.,
                  d_particles + particlePitch*j, centroid[i]));
            else
              gpuParticleKahanReductionAsync(npCount,
                  &Sijarr[i][j], &errorSig[j],
                  gpuSij(d_center, 0., d_center, 0.));
	    if (flags&BEAM_SUMS_EXACTEMIT && i<4 && j<4)
	      gpuParticleKahanReductionAsync(npCount,
		  &Sijnarr[i][j], &errorSign,
                  gpuSijn(i, centroidn[i], j, centroidn[j], d_pz));
#endif
            }
          }
          else { /* Single particle case */
#ifndef USE_KAHAN
            if (i<6 && j<6)
              gpuParticleReductionAsync(npCount,
                  &Sijarr[i][j], gpuSij(d_particles + particlePitch*i,
                  centroid[i], d_particles + particlePitch*j, centroid[j]),
                  Add<double>());
            else if (i<6)
              gpuParticleReductionAsync(npCount,
                  &Sijarr[i][j], gpuSij(d_particles + particlePitch*i,
                  centroid[i], d_center, 0.), Add<double>());
            else if (j<6)
              gpuParticleReductionAsync(npCount,
                  &Sijarr[i][j], gpuSij(d_center, 0.,
                  d_particles + particlePitch*j, centroid[j]), Add<double>());
            else
              gpuParticleReductionAsync(npCount,
                  &Sijarr[i][j], gpuSij(d_center, 0., d_center, 0.),
                  Add<double>());
	    if (flags&BEAM_SUMS_EXACTEMIT && i<4 && j<4)
	      gpuParticleReductionAsync(npCount, &Sijnarr[i][j],
		  gpuSijn(i, centroidn[i], j, centroidn[j], d_pz));
#else
            if (i<6 && j<6)
              gpuParticleKahanReductionAsync(npCount,
                  &Sijarr[i][j], &errorSig[j],
                  gpuSij(d_particles + particlePitch*i, centroid[i],
                  d_particles + particlePitch*j, centroid[j]));
            else if (i<6)
              gpuParticleKahanReductionAsync(npCount,
                  &Sijarr[i][j], &errorSig[j],
                  gpuSij(d_particles + particlePitch*i, centroid[i],
                  d_center, 0.));
            else if (j<6)
              gpuParticleKahanReductionAsync(npCount,
                  &Sijarr[i][j], &errorSig[j], gpuSij(d_center, 0.,
                  d_particles + particlePitch*j, centroid[i]));
            else
              gpuParticleKahanReductionAsync(npCount,
                  &Sijarr[i][j], &errorSig[j],
                  gpuSij(d_center, 0., d_center, 0.));
	    if (flags&BEAM_SUMS_EXACTEMIT && i<4 && j<4)
	      gpuParticleKahanReductionAsync(npCount,
		  &Sijnarr[i][j], &errorSign,
                  gpuSijn(i, centroidn[i], j, centroidn[j], d_pz));
#endif
        }
      }
    }
    finishReductionStreams();
    for (i=0; i<7; i++) {
      for (j=i; j<7; j++) {
        if (notSinglePart) {
          if (parallelStatus==trueParallel) {
            /* nothing to do */
          }
          else if (isMaster) {
            if (npCount) {
              sums->sigma[j][i] = (sums->sigma[i][j] = (sums->sigma[i][j]*sums->n_part+Sijarr[i][j])/(sums->n_part+npCount));
              if (flags&BEAM_SUMS_EXACTEMIT) {
		sums->sigman[j][i] = (sums->sigman[i][j] = (sums->sigman[i][j]*sums->n_part+Sijnarr[i][j])/(sums->n_part+npCount));
	      }
	    }
          }
        }
        else { /* Single particle case */
          if (n_part) {
	    if (flags&BEAM_SUMS_EXACTEMIT) {
	      sums->sigman[j][i] = (sums->sigman[i][j] = (sums->sigman[i][j]*sums->n_part-npCount+Sijnarr[i][j])/(sums->n_part));
	    } else {
	      sums->sigma[j][i] = (sums->sigma[i][j] = (sums->sigma[i][j]*sums->n_part+Sijarr[i][j])/(sums->n_part+npCount));
	    }
	  }
        }
      }
    }
  }

  if (notSinglePart) {
    if (parallelStatus==trueParallel) {
      if (isMaster) {
        memset(Sij_p, 0,  sizeof(double)*28);
#ifdef USE_KAHAN
        memset(errorSig, 0,  sizeof(double)*28);
#endif
      }
      /* compute Sij sum over processors */
#ifndef USE_KAHAN
      MPI_Allreduce(Sij_p, Sij_total, 28, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
      MPI_Allgather(Sij_p, 28, MPI_DOUBLE, &sumMatrixSig[0][0], 28, MPI_DOUBLE, MPI_COMM_WORLD);
      /* compute error sum over processors */
      MPI_Allgather(errorSig, 28, MPI_DOUBLE, &errorMatrixSig[0][0], 28, MPI_DOUBLE,MPI_COMM_WORLD);
      offset = 0;
      for (i=0; i<7; i++) {
        if (i>=1)        
          offset += i-1;
        for (j=i; j<7; j++) {
          index = 6*i+j-offset;
          if (flags&BEAM_SUMS_SPARSE && !sparse[i][j]) {
            /* Only compute the diagonal blocks and dispersive correlations */
            Sij_total[index] = 0;
            continue;
          }
          error_sum = 0.0;
          /* extract the columnwise array from the matrix */
          for (k=0; k<n_processors; k++) {         
            sumArray[k] = sumMatrixSig[k][index];             
            errorArray[k] = errorMatrixSig[k][index];
          }
          Sij_total[index] = Kahan(n_processors-1,&sumArray[1],&error_sum);
          error_total = Kahan(n_processors-1,&errorArray[1],&error_sum);
          Sij_total[index] += error_total;
        } 
      }
      
#endif
      if (npCount_total) {
        offset = 0; 
        for (i=0; i<7; i++) {
          if (i>=1)
            offset += i-1;
          for (j=i; j<7; j++) {
            if (flags&BEAM_SUMS_SPARSE && !sparse[i][j]) {
              /* Only compute the diagonal blocks and dispersive correlations */
              if (flags&BEAM_SUMS_EXACTEMIT)
                sums->sigman[j][i] = sums->sigman[i][j] = 0;
              else
		sums->sigma[j][i] = sums->sigma[i][j] = 0;
              continue;
            }
            index = 6*i+j-offset;
            if (flags&BEAM_SUMS_EXACTEMIT)
                /* sums->n_part was updated already, so we need to *subtract* npCount, not add it */
              sums->sigman[j][i] = (sums->sigman[i][j] = (sums->sigman[i][j]*(sums->n_part-npCount_total)+Sij_total[index])/sums->n_part);
            else
	      sums->sigma[j][i] = (sums->sigma[i][j] = (sums->sigma[i][j]*sums->n_part+Sij_total[index])/(sums->n_part+npCount_total));
          }
        }
      }
    }
  }

  //sums->charge = mp_charge*npCount_total;
  
  if (!(flags&BEAM_SUMS_EXACTEMIT)) {
    if (!notSinglePart) {
      sums->n_part += npCount; 
      sums->charge = mp_charge*npCount;
    } else if (!SDDS_MPI_IO) {
      if (parallelStatus==trueParallel) {
	sums->n_part += npCount_total;
	sums->charge = mp_charge*npCount_total;
      } else if (isMaster) {
	sums->n_part += npCount;
	sums->charge = mp_charge*npCount;
      }
    } else {
      if (isMaster) {
	sums->n_part += npCount_total;
	sums->charge = mp_charge*npCount_total;
      } else {
	sums->n_part += npCount;
	sums->charge = mp_charge*npCount;
      }
    }
  }

#ifdef USE_KAHAN
  free_czarray_2d((void**)sumMatrixCen, n_processors, 7); 
  free_czarray_2d((void**)errorMatrixCen, n_processors, 7);
  free_czarray_2d((void**)sumMatrixSig, n_processors, 28);
  free_czarray_2d((void**)errorMatrixSig, n_processors, 28);
  free(sumArray);
  free(errorArray);
#endif

  resetToFullParticleArray();
}
#endif

} // extern "C"
