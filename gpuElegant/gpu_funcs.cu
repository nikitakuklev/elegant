#include <gpu_track.h>

#include <gpu_base.h>
#include <gpu_particle_template_function.hcu>
#include <gpu_particle_reduction.hcu>
#include <gpu_reductions.h>
#include <gpu_trwake.h> // gpu_computeTimeCoordinatesAndMean, gpu_computeDistanceCoordinates
#include <gpu_matrix.h> // gpu_track_particles

class gpu_rotateBeamCoordinates_kernel {
public:

  double cos_a, sin_a;

  gpu_rotateBeamCoordinates_kernel(double cos_a, double sin_a) :
    cos_a(cos_a), sin_a(sin_a) {};

  __device__ void operator()(gpuParticleAccessor& coord) {

    double x, xp, y, yp;
    x = coord[0];
    xp = coord[1];
    y = coord[2];
    yp = coord[3];
    coord[0] =   x*cos_a + y*sin_a;
    coord[2] =  -x*sin_a + y*cos_a;
    coord[1] =  xp*cos_a + yp*sin_a;
    coord[3] = -xp*sin_a + yp*cos_a;
  }
};

/* The name is a misnomer: we rotate the coordinates of a particle
 * to emulate rotation of the upcoming element by the given angle
 */

void gpu_rotateBeamCoordinates(long np, double angle)
{
  double sin_a, cos_a;

  if (!angle || fabs(fabs(angle)-PIx2)<1e-12)
    return;
  if (fabs(fabs(angle)-PI)<1e-12) {
    cos_a = -1;
    sin_a = 0;
  }
  else if (fabs(angle-PIo2)<1e-12) {
    cos_a = 0;
    sin_a = 1;
  }
  else if (fabs(angle+PIo2)<1e-12) {
    cos_a = 0;
    sin_a = -1;
  }
  else {
    cos_a = cos(angle);
    sin_a = sin(angle);
  }
  
  gpuDriver(np,
      gpu_rotateBeamCoordinates_kernel(cos_a, sin_a));

}

#ifdef sqr 
#undef sqr
#endif
#define sqr(x) x*x

class gpu_offsetBeamCoordinates_kernel {
public:

  double dx, dy, dz;

  gpu_offsetBeamCoordinates_kernel(double dx, double dy, double dz) :
    dx(dx), dy(dy), dz(dz) {};

  __device__ void operator()(gpuParticleAccessor& part) {

    part[4] += dz*sqrt(1 + sqr(part[1]) + sqr(part[3]));
    part[0]  = part[0] - dx + dz*part[1];
    part[2]  = part[2] - dy + dz*part[3];
  }
};

/* The name is a misnomer: we offset the coordinates of a beam
 * to emulate offsetting of the upcoming element by the values
 */

void gpu_offsetBeamCoordinates(long np, double dx, double dy, double dz)
{ 
  gpuDriver(np, gpu_offsetBeamCoordinates_kernel(dx, dy, dz));
}

class gpu_center_beam_kernel {
public:
  double offset[7];
  double *d_timeCoord;

  gpu_center_beam_kernel(double o0, double o1, double o2, double o3,
      double o4, double o5, double o6, double* d_timeCoord) : 
      d_timeCoord(d_timeCoord) {
    offset[0]=o0; offset[1]=o1; offset[2]=o2; offset[3]=o3;
    offset[4]=o4; offset[5]=o5; offset[6]=o6;
  };

  __device__ void operator()(gpuParticleAccessor& coord) {
    unsigned int tid = coord.getParticleIndex();

    for (int ic=0; ic<6; ic++)
      if (offset[ic])
        coord[ic] -= offset[ic];
    if (offset[6])
      d_timeCoord[tid] -= offset[6];
  }
};

extern "C" {

void gpu_center_beam(CENTER *center, long np, long iPass, double p0)
{
  double sum[6], delta6(0);
  long ic;
  struct GPUBASE* gpubase = getGpuBase();
  double* d_particles = gpubase->d_particles;
  unsigned int particlePitch = gpubase->gpu_array_pitch;
  double* d_timeCoord = gpubase->d_temp_particles;

  if (!np) {
    return;
  }

  if (center->onPass>=0 && iPass!=center->onPass)
    return;

  for (ic=0; ic<6; ic++) {
    sum[ic]=0;
    if (center->doCoord[ic]) {
      if (!center->deltaSet[ic]) {
        gpuReduceAddAsync(d_particles+particlePitch*ic, np, &sum[ic]);
      }
    }
  }
  finishReductionStreams();
  for (ic=0; ic<6; ic++) {
    if (center->doCoord[ic]) {
      if (!center->deltaSet[ic]) {
#if USE_MPI
	if (notSinglePart) {
	  double sum_total;
          long np_total;

          MPI_Allreduce (&sum[ic], &sum_total, 1, MPI_DOUBLE, MPI_SUM, workers);
          sum[ic] = sum_total;
	  MPI_Allreduce (&np, &np_total, 1, MPI_LONG, MPI_SUM, workers);
          center->delta[ic] = sum[ic]/np_total;
	} else
	  center->delta[ic] = sum[ic]/np;
#else
        center->delta[ic] = sum[ic]/np;
#endif
        if (center->onceOnly)
          center->deltaSet[ic] = 1;
      } 
    }
  }

  if (center->doCoord[ic=6]) {
    /* Special treatment for time coordinate */
    delta6 = gpu_computeTimeCoordinatesAndMean(np, d_timeCoord, p0);
    if (center->deltaSet[ic])
      delta6 = center->delta[ic];
  }

  gpuDriver(np, 
    gpu_center_beam_kernel(center->delta[0], center->delta[1], center->delta[2],
    center->delta[3], center->delta[4], center->delta[5], delta6, 
    d_timeCoord));

  if (center->doCoord[ic=6]) {
    gpu_computeDistanceCoordinates(d_timeCoord, p0, np);
    if (center->onceOnly && !center->deltaSet[ic]) {
      center->delta[ic] = delta6;
      center->deltaSet[ic] = 1;
    }
  }
}

} // extern "C"

class gpuAddCorrectorRadiationKick {
public:
  gpuAddCorrectorRadiationKick(short isr, double Po, double F2, double length,
      double isrCoef, double radCoef, short do_reduce, double* d_gauss_rn){
    this->isr = isr;
    this->Po = Po;
    this->F2 = F2;
    this->length = length;
    this->isrCoef = isrCoef;
    this->radCoef = radCoef;
    this->do_reduce = do_reduce;
    this->d_gauss_rn = d_gauss_rn;
  }
  __device__ double inline operator()(gpuParticleAccessor& particle){
    double dp, p, beta0, beta1, deltaFactor;
    unsigned int ind = particle.getParticleIndex();

    dp = particle[5];
    p = Po*(1+dp);
    beta0 = p/sqrt(pow(p,2)+1);
    deltaFactor = pow(1+dp,2);
    dp -= radCoef*deltaFactor*F2*length;
    if (isr)
      dp += isrCoef*deltaFactor*pow(F2, 0.75)*sqrt(length)*d_gauss_rn[ind];
    p = Po*(1+dp);
    beta1 = p/sqrt(pow(p,2)+1);
    particle[5] = dp;
    particle[4] = beta1*particle[4]/beta0;
    if (do_reduce) return pow(isrCoef*deltaFactor,2)*pow(F2, 1.5)*length;
    else return 0;
  }
  short isr, do_reduce;
  double Po, F2, length, isrCoef, radCoef;
  double* d_gauss_rn;
};

extern "C" {

void gpu_addCorrectorRadiationKick(long np, ELEMENT_LIST *elem, long type, 
       double Po, double *sigmaDelta2, long disableISR){
  double F2;
  double kick, length;
  double isrCoef, radCoef;
  short isr, sr;

  if (!np)
    return;

  isr = sr = 0;

  switch (type) {
  case T_HCOR:
    kick = ((HCOR*)elem->p_elem)->kick;
    if ((length = ((HCOR*)elem->p_elem)->length)==0)
      length = ((HCOR*)elem->p_elem)->lEffRad;
    if (((HCOR*)elem->p_elem)->synchRad) {
      sr = 1;
      if (((HCOR*)elem->p_elem)->isr)
        isr = 1;
    }
    break;
  case T_VCOR:
    kick = ((VCOR*)elem->p_elem)->kick;
    if ((length = ((VCOR*)elem->p_elem)->length)==0)
      length = ((VCOR*)elem->p_elem)->lEffRad;
    if (((VCOR*)elem->p_elem)->synchRad) {
      sr = 1;
      if (((VCOR*)elem->p_elem)->isr)
        isr = 1;
    }
    break;
  case T_HVCOR:
    kick = sqrt(sqr(((HVCOR*)elem->p_elem)->xkick)+sqr(((HVCOR*)elem->p_elem)->ykick));
    if ((length = ((HVCOR*)elem->p_elem)->length)==0)
      length = ((HVCOR*)elem->p_elem)->lEffRad;
    if (((HVCOR*)elem->p_elem)->synchRad) {
      sr = 1;
      if (((HVCOR*)elem->p_elem)->isr)
        isr = 1;
    }
    break;
  case T_EHCOR:
    kick = ((EHCOR*)elem->p_elem)->kick;
    if ((length = ((EHCOR*)elem->p_elem)->length)==0)
      length = ((EHCOR*)elem->p_elem)->lEffRad;
    if (((EHCOR*)elem->p_elem)->synchRad) {
      sr = 1;
      if (((EHCOR*)elem->p_elem)->isr)
        isr = 1;
    }
    break;
  case T_EVCOR:
    kick = ((EVCOR*)elem->p_elem)->kick;
    if ((length = ((EVCOR*)elem->p_elem)->length)==0)
      length = ((EVCOR*)elem->p_elem)->lEffRad;
    if (((EVCOR*)elem->p_elem)->synchRad) {
      sr = 1;
      if (((EVCOR*)elem->p_elem)->isr)
        isr = 1;
    }
    break;
  case T_EHVCOR:
    kick = sqrt(sqr(((EHVCOR*)elem->p_elem)->xkick)+sqr(((EHVCOR*)elem->p_elem)->ykick));
    if ((length = ((EHVCOR*)elem->p_elem)->length)==0)
      length = ((EHVCOR*)elem->p_elem)->lEffRad;
    if (((EHVCOR*)elem->p_elem)->synchRad) {
      sr = 1;
      if (((EHVCOR*)elem->p_elem)->isr)
        isr = 1;
    }
    break;
  }
  if (sr==0 || length==0)
    return ;
  if (disableISR)
    isr = 0;
  radCoef = sqr(particleCharge)*pow3(Po)/(6*PI*epsilon_o*sqr(c_mks)*particleMass);
  isrCoef = particleRadius*sqrt(55.0/(24*sqrt(3))*pow5(Po)*137.0359895);

  F2 = sqr(kick/length);

  struct GPUBASE* gpuBase = getGpuBase();
  double* d_gauss_rn=gpuBase->d_temp_particles;

  if (isr>0) 
    gpu_d_gauss_rn_lim(d_gauss_rn, np, 0.0, 1.0, 3.0, random_2(0));

  if (sigmaDelta2) 
    *sigmaDelta2 = gpuParticleReduction(np,
      gpuAddCorrectorRadiationKick(isr, Po, F2, length, isrCoef, radCoef,
      1, d_gauss_rn), Add<double>());
  else
    gpuDriver(np, 
      gpuAddCorrectorRadiationKick(isr, Po, F2, length, isrCoef, radCoef, 
      0, d_gauss_rn));
}

} // extern "C"

class gpu_offset_beam_kernel{
public:
  double dx, dy, dz, dxp, dyp, dt, dp, de, P_central;
  bool allParticles;
  double startPID, endPID;

  gpu_offset_beam_kernel(double dx, double dy, double dz,
    double dxp, double dyp, double dt, double dp, double de, 
    double P_central, bool allParticles, long startPID, long endPID) :
    dx(dx), dy(dy), dz(dz), dxp(dxp), dyp(dyp), dt(dt), de(de),
    P_central(P_central), allParticles(allParticles), 
    startPID(startPID), endPID(endPID) {};

  __inline__ __device__ void operator()(gpuParticleAccessor& part){
    double pc, beta, gamma, t;
    double ds;

    if (!allParticles && (part[6]<startPID || part[6]>endPID))
      return;
    if (dz)
      ds = dz*sqrt(1+sqr(part[1])+sqr(part[3]));
    else
      ds = 0;
    part[0] += dx + dz*part[1];
    part[1] += dxp;
    part[2] += dy + dz*part[3];
    part[3] += dyp;
    part[4] += ds;
    if (dt || dp || de) {
      pc = P_central*(1+part[5]);
      beta = pc/(gamma=sqrt(1+pc*pc));
      t = part[4]/(beta*c_mks) + dt;
      if (dp) {
        part[5] += dp;
        pc = P_central*(1+part[5]);
        beta = pc/sqrt(1+pc*pc);
      }
      if (de) {
        gamma += de*gamma;
        pc = sqrt(gamma*gamma-1);
        beta = pc/gamma;
        part[5] = (pc-P_central)/P_central;
      }
      part[4] = t*beta*c_mks;
    }
  }
};

extern "C" {

void gpu_offset_beam(long nToTrack, MALIGN *offset, double P_central) {
  bool allParticles;
  log_entry("offset_beam");

  if (offset->startPID>=0 && offset->startPID>offset->endPID)
    bombElegantVA("Error: startPID (%ld) greater than endPID (%ld) for MALIGN element (offset_beam)\n", offset->startPID, offset->endPID);
  if ((offset->endPID>=0 && offset->startPID<0) || (offset->startPID>=0 && offset->endPID<0))
    bombElegantVA("Error: Invalid startPID (%ld) and endPID (%ld) in MALIGN element (offset_beam)\n", offset->startPID, offset->endPID);

  allParticles = (offset->startPID==-1) && (offset->endPID==-1);

  gpuDriver(nToTrack,
    gpu_offset_beam_kernel(offset->dx, offset->dy, offset->dz, offset->dxp,
      offset->dyp, offset->dt, offset->dp, offset->de, P_central,
      allParticles, offset->startPID, offset->endPID));

  log_exit("offset_beam");
}

} // extern "C"

class do_match_energy_P_average_transform{
public:
  double P_central;
  do_match_energy_P_average_transform(double P_central) : 
    P_central(P_central) {};

  __inline__ __device__ double operator()(gpuParticleAccessor& coord){
    double P_average = P_central * (1 + coord[5]);
    return P_average;
  }
};

class do_match_energy_adjust_to_P_average{
 public:
  double P_central, P_average, dp, dr;

  do_match_energy_adjust_to_P_average(double P_central, double P_average,
    double dp, double dr) : P_central(P_central), P_average(P_average), 
    dp(dp), dr(dr) {};

  __inline__ __device__ void operator()(gpuParticleAccessor& coord){
    double dp = (P_central - P_average)/P_average;
    double dr = P_central/P_average;
    //coord[5] = ((1+coord[5])*P_central - P_average)/P_average;
    coord[5] = dp + coord[5]*dr;
  }
};

class do_match_energy_match_energy{
public:
  double dP_centroid;
  double P_central;
  do_match_energy_match_energy(double dP_centroid, double P_central) : 
    dP_centroid(dP_centroid), P_central(P_central) {};

  __inline__ __device__ void operator()(gpuParticleAccessor& coord){    
    double P = (1+coord[5])*(P_central);
    double t = coord[4]/(P/sqrt(P*P+1));
    P += dP_centroid;
    coord[5] = (P - P_central)/ (P_central);
    coord[4] = t*(P/sqrt(P*P+1));

#if defined(IEEE_MATH)
    if(!isfinite(coord[4])){
      printf("GPU Error: bad time coordinate for particle %d:  %g %g %g  %g %g %g, due to P_central %g, t %g, dP_centroid %g\n",
	     coord.getParticleIndex(),coord[0],coord[1],coord[2],coord[3],coord[4],coord[5],P_central, t, dP_centroid);
      //assert(isfinite(coord[4]));
    }
#endif
    
  }
};

extern "C" {

void gpu_do_match_energy(long np, double* P_central, long change_beam) {

  double P_average, dP_centroid;
  long active = 1;
#ifdef USE_KAHAN
  double error = 0.0;
#endif
#if USE_MPI
  long np_total;
  double P_total = 0.0;
  if (notSinglePart) {
    if (((parallelStatus==trueParallel) && isSlave) || ((parallelStatus!=trueParallel) && isMaster))
      active = 1;
    else 
      active = 0;
  }  
#endif

  log_entry("do_match_energy");

#if (!USE_MPI)  
  if (!np) {
    log_exit("do_match_energy");
    return;
  }
#else
  if (notSinglePart) {
    if (parallelStatus!=trueParallel) {
      if (!np) {
	log_exit("do_match_energy");   
	return;   
      }
    }
    else {
      if (isMaster) 
	np = 0; /* All the particles have been distributed to the slave processors */    
      MPI_Allreduce(&np, &np_total, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);        
    }    
  }
  else if (!np) {
    log_exit("do_match_energy");
    return;
  }
#endif
  
  if (!change_beam) {
    /* change the central momentum so that it matches the beam's centroid */
    P_average = 0;
    if (active) {

      do_match_energy_P_average_transform P_average_transform(*P_central);
#ifndef USE_KAHAN	     
      P_average = gpuParticleReduction(np, 
                    P_average_transform, Add<double>() );
#else

      P_average = gpuParticleKahanReduction(np, 
                    &error, P_average_transform);
#endif	

    }
#if (!USE_MPI)
    P_average /= np;
#else
    if (notSinglePart) {
      if (parallelStatus!=trueParallel) {
	if (isMaster)            
	  P_average /= np; 
      }
      else {
#ifndef USE_KAHAN    
	MPI_Allreduce(&P_average, &P_total, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );  
#else
        P_total = KahanParallel (P_average, error, MPI_COMM_WORLD);

#endif 
	P_average = P_total/np_total;
      }
    }
    else /* Single particle case, all the processors will do the same as in serial version */
      P_average /= np;
#endif 

#ifdef DEBUG_FIDUCIALIZATION
      fprintf(stdout, "Changing reference momentum from %e to %e to match beam, ratio %e\n",
              *P_central, P_average, (*P_central/P_average-1)*active);
#endif
    if (fabs(P_average-(*P_central))/(*P_central)>1e-14){ 
     /* if (P_average!= *P_central) { */
      double dp = (*P_central - P_average)/P_average;
      double dr = (*P_central)/P_average;
      if (active) {
	do_match_energy_adjust_to_P_average adjust_to_P_average(*P_central,
           P_average, dp, dr);
	gpuDriver(np, adjust_to_P_average);
      }
      *P_central =  P_average;
    }
  }
  else {
    /* change the particle momenta so that the centroid is the central momentum */
    /* the path length is adjusted so that the time-of-flight at the current
       velocity is fixed */
    P_average = 0;
    if (active) {
      do_match_energy_P_average_transform P_average_transform(*P_central);
#ifndef USE_KAHAN	     	
      P_average = gpuParticleReductionDriver(np, P_average_transform, Add<double>());
#else
      P_average = gpuParticleKahanReduction(np, &error, P_average_transform);      
#endif	
      
    }
#if (!USE_MPI)
    P_average /= np;
#else
    if (notSinglePart) { 
      if (parallelStatus!=trueParallel) {
	if (isMaster)
	  P_average /= np; 
      }
      else {
#ifndef USE_KAHAN    
	MPI_Allreduce(&P_average, &P_total, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );  
#else
        P_total = KahanParallel (P_average, error, MPI_COMM_WORLD);
#endif
	P_average = P_total/np_total;
      }
    }
    else
      P_average /= np;
#endif       

    if (active) {

      dP_centroid =  *P_central - P_average;

      do_match_energy_match_energy match_energy(dP_centroid, *P_central);
      
      gpuDriver(np, match_energy);

      cudaError_t err = cudaThreadSynchronize();
      if(err!=0){
	fprintf(stdout,"GPU Error or caught assert(isfinite(coord[4])) after trying to match energy\n");
	fprintf(stdout, "P_average = %e  P_central = %e  dP_centroid = %e\n",
		P_average, *P_central, dP_centroid);
	fflush(stdout);
#if (USE_MPI)
	if (active)
	  MPI_Abort(MPI_COMM_WORLD, 1);
#endif    
	abort();
      
      }
    }
  }
  log_exit("do_match_energy");
}

} // extern "C"

class apply_offset4{
public:
  double offset;
  apply_offset4(double offset) : offset(offset) {};

  __inline__ __device__ void operator()(gpuParticleAccessor& coord){
    coord[4] += offset;
  }
};

extern "C" {

void gpu_matr_element_tracking(VMATRIX *M, MATR *matr,
                               long np, double z)
/* subtract off <s> prior to using a user-supplied matrix to avoid possible
 * problems with R5? and T?5? elements
 */
{
  struct GPUBASE* gpuBase = getGpuBase();
#if !USE_MPI
  if (!np)
    return;
#endif
  if (!matr) {
    gpu_track_particles(M, np);
  } else {
    if (!matr->fiducialSeen) {
      double sum = 0;
      sum = gpuReduceAdd(gpuBase->d_particles+4*gpuBase->gpu_array_pitch, np);
#if !USE_MPI
      matr->sReference = sum/np;
#else
      if (notSinglePart) {
	if (isSlave) {
	  double sum_total;
	  long np_total;

	  MPI_Allreduce(&np, &np_total, 1, MPI_LONG, MPI_SUM, workers);
	  MPI_Allreduce(&sum, &sum_total, 1, MPI_DOUBLE, MPI_SUM, workers);	      
	  matr->sReference = sum_total/np_total;
	}
      } else
	matr->sReference = sum/np;
#endif
      matr->fiducialSeen = 1;
    }
    gpuDriver(np, apply_offset4(-matr->sReference));
    gpu_track_particles(M, np);
    gpuDriver(np, apply_offset4(matr->sReference));
  }
}

void gpu_ematrix_element_tracking(VMATRIX *M, EMATRIX *matr,
			      long np, double z, double *P_central)
/* subtract off <s> prior to using a user-supplied matrix to avoid possible
 * problems with R5? and T?5? elements
 */
{
  struct GPUBASE* gpuBase = getGpuBase();
#if !USE_MPI
  if (!np)
    return;
#endif
  if (!matr) {
    fprintf(stderr, "ematrix_element_tracking: matr=NULL, tracking with M (%ld order)\n",
            M->order);
    gpu_track_particles(M, np);
  } else {
    if (!matr->fiducialSeen) {
      double sum = 0;
      sum = gpuReduceAdd(gpuBase->d_particles+4*gpuBase->gpu_array_pitch, np);
#if !USE_MPI
      matr->sReference = sum/np;
#else
      if (notSinglePart) {
	if (isSlave) {
	  double sum_total;
	  long np_total;
	
	  MPI_Allreduce(&np, &np_total, 1, MPI_LONG, MPI_SUM, workers);
	  MPI_Allreduce(&sum, &sum_total, 1, MPI_DOUBLE, MPI_SUM, workers);	      
	  matr->sReference = sum_total/np_total;
	}
      } else
	matr->sReference = sum/np;
#endif
      matr->fiducialSeen = 1;
    }

    gpuDriver(np, apply_offset4(-matr->sReference));
    gpu_track_particles(M, np);
    gpuDriver(np, apply_offset4(matr->sReference));
  }
  if (matr->deltaP)
    *P_central += matr->deltaP;
}

void gpu_collect_trajectory_data(double* centroid, long np) {
  struct GPUBASE* gpubase = getGpuBase();
  double* d_particles = gpubase->d_particles;
  unsigned int particlePitch = gpubase->gpu_array_pitch;
  double sums[6];
  for (int i_coord=0; i_coord<6; i_coord++) 
    gpuReduceAddAsync(d_particles+particlePitch*i_coord, np, &sums[i_coord]);
  finishReductionStreams();
  for (int i_coord=0; i_coord<6; i_coord++) 
    centroid[i_coord] = sums[i_coord]/np;
}

} //extern "C"

class gpu_set_central_momentum_kernel{
public:
  double P_central, P_new;
  gpu_set_central_momentum_kernel(double P_central, double P_new) : 
    P_central(P_central), P_new(P_new) {};

  __inline__ __device__ void operator()(gpuParticleAccessor& coord){
    coord[5] = ((1+coord[5])*P_central - P_new)/P_new;
  }
};

extern "C" {

void gpu_set_central_momentum(long np, double  P_new, double *P_central) {

#if (!USE_MPI)  
  if (!np) {
    *P_central =  P_new;
    return;
  }
  if (*P_central != P_new) {
    gpuDriver(np, gpu_set_central_momentum_kernel(*P_central, P_new));
    *P_central =  P_new;
  }
#else
  if (notSinglePart) {
    if (!np)
      *P_central =  P_new;

    if (*P_central != P_new) {
#ifdef DEBUG_FIDUCIALIZATION
      fprintf(stdout, "Changing reference momentum from %e to %e in %s at %e to match beam\n",
              *P_central, P_new, trackingContext.elementName, trackingContext.zEnd);
#endif
      if (((parallelStatus==trueParallel) && isSlave) || ((parallelStatus!=trueParallel) && isMaster)) {
        gpuDriver(np, gpu_set_central_momentum_kernel(*P_central, P_new));
      }
      *P_central =  P_new;
    }
  }
  else {
    if (!np) {
      *P_central =  P_new;
      return;
    }
    if (*P_central != P_new) {
#ifdef DEBUG_FIDUCIALIZATION
      fprintf(stdout, "Changing reference momentum from %e to %e in %s\n",
              *P_central, P_new, trackingContext.elementName);
#endif
      gpuDriver(np, gpu_set_central_momentum_kernel(*P_central, P_new));
      *P_central =  P_new;
    }
  }
#endif
}

} //extern "C"
