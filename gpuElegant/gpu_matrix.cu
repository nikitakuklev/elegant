#include <iostream>
#include <string>
#include <stdio.h>

#include <gpu_matrix.h>

__device__ __constant__ double constData_C[6];
__device__ __constant__ double constData_R[36];
__device__ __constant__ double constData_T[216];
__device__ __constant__ double constData_Q[6*6*6*6];

void gpu_matrix_copyCRTtoGPU(double* C, double* R){
  cudaError_t err;
  err = cudaMemcpyToSymbol(constData_C, C, sizeof(double)*6, 0, cudaMemcpyHostToDevice);
  gpuErrorHandler(err,"after memcpy to symbol constData_C");
  err = cudaMemcpyToSymbol(constData_R, R, sizeof(double)*6*6, 0, cudaMemcpyHostToDevice);
  gpuErrorHandler(err,"after memcpy to symbol constData_R");
}

void gpu_matrix_copyCRTtoGPU(double* C, double* R, double* T){
  cudaError_t err;
  err = cudaMemcpyToSymbol(constData_C, C, sizeof(double)*6, 0, cudaMemcpyHostToDevice);
  gpuErrorHandler(err,"after memcpy to symbol constData_C");
  err = cudaMemcpyToSymbol(constData_R, R, sizeof(double)*6*6, 0, cudaMemcpyHostToDevice);
  gpuErrorHandler(err,"after memcpy to symbol constData_R");
  err = cudaMemcpyToSymbol(constData_T, T, sizeof(double)*6*6*6, 0, cudaMemcpyHostToDevice);
  gpuErrorHandler(err,"after memcpy to symbol constData_T");
}

void gpu_matrix_copyCRTtoGPU(double* C, double* R, double* T, double* Q){
  cudaError_t err;
  err = cudaMemcpyToSymbol(constData_C, C, sizeof(double)*6, 0, cudaMemcpyHostToDevice);
  gpuErrorHandler(err,"after memcpy to symbol constData_C");
  err = cudaMemcpyToSymbol(constData_R, R, sizeof(double)*6*6, 0, cudaMemcpyHostToDevice);
  gpuErrorHandler(err,"after memcpy to symbol constData_R");
  err = cudaMemcpyToSymbol(constData_T, T, sizeof(double)*6*6*6, 0, cudaMemcpyHostToDevice);
  gpuErrorHandler(err,"after memcpy to symbol constData_T");
  err = cudaMemcpyToSymbol(constData_Q, Q, sizeof(double)*6*6*6*6, 0, cudaMemcpyHostToDevice);
  gpuErrorHandler(err,"after memcpy to symbol constData_Q");
}

__global__ void gpu_dipole_matrix_track_particles_CRT_Constant_Pragma_kernel(
    double *d_final_0, double *d_final_1, double *d_final_2, double *d_final_3,
    double *d_final_4, double *d_final_5, double *d_initial_0, 
    double *d_initial_1, double *d_initial_2, double *d_initial_3,
    double *d_initial_4, double *d_initial_5, int n_part) {
 
  int tx = threadIdx.x;
  int bx = blockIdx.x;
  int tid = bx*blockDim.x + tx;

  double t_initial[6];
  double t_final[6];
  for(; tid < n_part; tid += gridDim.x*blockDim.x){

    t_initial[0] = d_initial_0[tid];
    t_initial[1] = d_initial_1[tid];
    t_initial[2] = d_initial_2[tid];
    t_initial[3] = d_initial_3[tid];
    t_initial[4] = d_initial_4[tid];
    t_initial[5] = d_initial_5[tid];
  
#pragma unroll 6
    for(int i=0; i<6; i++){
      t_final[i]= constData_C[i];
    }
#pragma unroll 6
    for(int i=0; i<6; i++){
      double alpha = 0;
#pragma unroll 6
      for(int j=0; j<6; j++){
	alpha += constData_R[i*6+j]*t_initial[j];
      }
      t_final[i] += alpha;
    }
    d_final_0[tid] = t_final[0];
    d_final_1[tid] = t_final[1];
    d_final_2[tid] = t_final[2];
    d_final_3[tid] = t_final[3];
    d_final_4[tid] = t_final[4];
    d_final_5[tid] = t_final[5];
  }//for
  
}

void gpu_dipole_matrix_track_particles_launch_kernel(double *d_final_0, 
    double *d_final_1, double *d_final_2, double *d_final_3, double *d_final_4,
    double *d_final_5, double *d_initial_0, double *d_initial_1,
    double *d_initial_2, double *d_initial_3, double *d_initial_4,
    double *d_initial_5, int n_part) {

  struct GPUBASE* gpuBase = getGpuBase();
  unsigned int nTx = gpuBase->nThreads;
  unsigned int nBx = (n_part + nTx - 1) / nTx;
  if(nBx > gpuBase->maxBlocks) nBx = gpuBase->maxBlocks;
  dim3 dimBlock(nTx,1,1);
  dim3 dimGrid(nBx,1,1);
  
  gpu_dipole_matrix_track_particles_CRT_Constant_Pragma_kernel
    <<<dimGrid,dimBlock>>>(d_final_0, d_final_1, d_final_2, d_final_3,
        d_final_4, d_final_5, d_initial_0, d_initial_1, d_initial_2, 
        d_initial_3, d_initial_4, d_initial_5, n_part);

  gpuErrorHandler("gpu_dipole_matrix_track_particles_launch_kernel");
}

void gpu_dipole_matrix_track_particles(VMATRIX *M, unsigned int n_part){
  struct GPUBASE* gpubase = getGpuBase();
  double* d_particles_initial = gpubase->d_particles;
  double* d_particles_final = gpubase->d_particles;
  unsigned int gpu_pitch = gpubase->gpu_array_pitch;
  double h_C[6];
  double h_R[6*6];
  log_entry("gpu dipole matrix track particles");
  
  for(unsigned int i=0; i < 6; i++){
    h_C[i] = M->C[i];
    //    printf("h_C[%d] = %g\n",i,h_C[i]);
    for(unsigned int j=0; j< 6; j++){
      h_R[i*6 + j] = M->R[i][j];
      //      printf("h_R[%d] = %g\n",i*6+j,h_C[i*6+j]);
    }
  }
  gpu_matrix_copyCRTtoGPU(h_C, h_R);

  double* d_initial_0 = d_particles_initial + 0 * gpu_pitch;
  double* d_initial_1 = d_particles_initial + 1 * gpu_pitch;
  double* d_initial_2 = d_particles_initial + 2 * gpu_pitch;
  double* d_initial_3 = d_particles_initial + 3 * gpu_pitch;
  double* d_initial_4 = d_particles_initial + 4 * gpu_pitch;
  double* d_initial_5 = d_particles_initial + 5 * gpu_pitch;

  double* d_final_0 = d_particles_final + 0 * gpu_pitch;
  double* d_final_1 = d_particles_final + 1 * gpu_pitch;
  double* d_final_2 = d_particles_final + 2 * gpu_pitch;
  double* d_final_3 = d_particles_final + 3 * gpu_pitch;
  double* d_final_4 = d_particles_final + 4 * gpu_pitch;
  double* d_final_5 = d_particles_final + 5 * gpu_pitch;
  
  gpu_dipole_matrix_track_particles_launch_kernel(
      d_final_0, d_final_1, d_final_2, d_final_3, d_final_4, d_final_5,
      d_initial_0, d_initial_1, d_initial_2, d_initial_3, d_initial_4,
      d_initial_5, n_part);

  log_exit("gpu dipole matrix track particles");
}
  
__global__ void gpu_track_particles_CRT_Constant_Pragma_kernel(
    double *d_final_0, double *d_final_1, double *d_final_2,
    double *d_final_3, double *d_final_4, double *d_final_5,
    double *d_initial_0, double *d_initial_1, double *d_initial_2,
    double *d_initial_3, double *d_initial_4, double *d_initial_5, int n_part){
 
  int tx = threadIdx.x;
  int bx = blockIdx.x;
  int tid = bx*blockDim.x + tx;

  double t_initial[6];
  double t_final[6];
  for(; tid < n_part; tid += gridDim.x*blockDim.x){

    t_initial[0] = d_initial_0[tid];
    t_initial[1] = d_initial_1[tid];
    t_initial[2] = d_initial_2[tid];
    t_initial[3] = d_initial_3[tid];
    t_initial[4] = d_initial_4[tid];
    t_initial[5] = d_initial_5[tid];
  
         
#pragma unroll 6
    for(int i=0; i<6; i++){
      t_final[i]= constData_C[i];
    }
#pragma unroll 6
    for(int i=0; i<6; i++){
      double alpha = 0;
#pragma unroll 6
      for(int j=0; j<6; j++){
	alpha += constData_R[i*6+j]*t_initial[j];
      }
      t_final[i] += alpha;
    }

#pragma unroll 6
    for(int i=0; i<6; i++){
      double alpha = 0;
      // j = 0, k = 0
      alpha+= constData_T[i*36 + 0*6 + 0]*t_initial[0]*t_initial[0];
      // j = 1, k = 0 to 1
      alpha+= constData_T[i*36 + 1*6 + 0]*t_initial[0]*t_initial[1];
      alpha+= constData_T[i*36 + 1*6 + 1]*t_initial[1]*t_initial[1];
      // j = 2, k = 0 to 2
      alpha+= constData_T[i*36 + 2*6 + 0]*t_initial[0]*t_initial[2];
      alpha+= constData_T[i*36 + 2*6 + 1]*t_initial[1]*t_initial[2];
      alpha+= constData_T[i*36 + 2*6 + 2]*t_initial[2]*t_initial[2];
      // j = 3, k = 0 to 3
      alpha+= constData_T[i*36 + 3*6 + 0]*t_initial[0]*t_initial[3];
      alpha+= constData_T[i*36 + 3*6 + 1]*t_initial[1]*t_initial[3];
      alpha+= constData_T[i*36 + 3*6 + 2]*t_initial[2]*t_initial[3];
      alpha+= constData_T[i*36 + 3*6 + 3]*t_initial[3]*t_initial[3];
      // j = 4, k = 0 to 4
      alpha+= constData_T[i*36 + 4*6 + 0]*t_initial[0]*t_initial[4];
      alpha+= constData_T[i*36 + 4*6 + 1]*t_initial[1]*t_initial[4];
      alpha+= constData_T[i*36 + 4*6 + 2]*t_initial[2]*t_initial[4];
      alpha+= constData_T[i*36 + 4*6 + 3]*t_initial[3]*t_initial[4];
      alpha+= constData_T[i*36 + 4*6 + 4]*t_initial[4]*t_initial[4];
      // j = 5, k = 0 to 5
      alpha+= constData_T[i*36 + 5*6 + 0]*t_initial[0]*t_initial[5];
      alpha+= constData_T[i*36 + 5*6 + 1]*t_initial[1]*t_initial[5];
      alpha+= constData_T[i*36 + 5*6 + 2]*t_initial[2]*t_initial[5];
      alpha+= constData_T[i*36 + 5*6 + 3]*t_initial[3]*t_initial[5];
      alpha+= constData_T[i*36 + 5*6 + 4]*t_initial[4]*t_initial[5];
      alpha+= constData_T[i*36 + 5*6 + 5]*t_initial[5]*t_initial[5];

      t_final[i] += alpha;
    }// for i
    d_final_0[tid] = t_final[0];
    d_final_1[tid] = t_final[1];
    d_final_2[tid] = t_final[2];
    d_final_3[tid] = t_final[3];
    d_final_4[tid] = t_final[4];
    d_final_5[tid] = t_final[5];
  }//for
  
}

void gpu_quad_matrix_track_particles_launch_kernel(double *d_final_0,
    double *d_final_1, double *d_final_2, double *d_final_3, double *d_final_4,
    double *d_final_5, double *d_initial_0, double *d_initial_1,
    double *d_initial_2, double *d_initial_3, double *d_initial_4,
    double *d_initial_5, int n_part){

  struct GPUBASE* gpuBase = getGpuBase();
  unsigned int nTx = gpuBase->nThreads;
  unsigned int nBx = (n_part + nTx - 1) / nTx;
  if(nBx > gpuBase->maxBlocks) nBx = gpuBase->maxBlocks;
  dim3 dimBlock(nTx,1,1);
  dim3 dimGrid(nBx,1,1);

  gpu_track_particles_CRT_Constant_Pragma_kernel<<<dimGrid,dimBlock>>>(
      d_final_0, d_final_1, d_final_2, d_final_3, d_final_4, d_final_5,
      d_initial_0, d_initial_1, d_initial_2, d_initial_3, d_initial_4,
      d_initial_5, n_part);

  gpuErrorHandler("gpu_quad_matrix_track_particles_launch_kernel");
}

void gpu_quad_matrix_track_particles(VMATRIX *M, unsigned int n_part){
  struct GPUBASE* gpubase = getGpuBase();
  double* d_particles_initial = gpubase->d_particles;
  double* d_particles_final = gpubase->d_particles;
  unsigned int gpu_pitch = gpubase->gpu_array_pitch;
  double h_C[6];
  double h_R[36];
  double h_T[6*6*6];
  log_entry("gpu quad matrix track particles");  
  
  for(unsigned int i=0; i < 6; i++){
    h_C[i] = M->C[i];
    //    printf("h_C[%d] = %g\n",i,h_C[i]);
    for(unsigned int j=0; j< 6; j++){
      h_R[i*6 + j] = M->R[i][j];
      //      printf("h_R[%d] = %g\n",i*6+j,h_C[i*6+j]);
      for(unsigned int k =0; k<= j; k++){	
	h_T[i*6*6 + j*6 + k] = M->T[i][j][k];
	//	printf("h_T[%d] = %g\n",i*6*6 + j*6 +k, h_T[i*6*6 + j*6 + k]);
      }
    }
  }
  gpu_matrix_copyCRTtoGPU(h_C, h_R, h_T);

  double* d_initial_0 = d_particles_initial + 0 * gpu_pitch;
  double* d_initial_1 = d_particles_initial + 1 * gpu_pitch;
  double* d_initial_2 = d_particles_initial + 2 * gpu_pitch;
  double* d_initial_3 = d_particles_initial + 3 * gpu_pitch;
  double* d_initial_4 = d_particles_initial + 4 * gpu_pitch;
  double* d_initial_5 = d_particles_initial + 5 * gpu_pitch;

  double* d_final_0 = d_particles_final + 0 * gpu_pitch;
  double* d_final_1 = d_particles_final + 1 * gpu_pitch;
  double* d_final_2 = d_particles_final + 2 * gpu_pitch;
  double* d_final_3 = d_particles_final + 3 * gpu_pitch;
  double* d_final_4 = d_particles_final + 4 * gpu_pitch;
  double* d_final_5 = d_particles_final + 5 * gpu_pitch;
  
  gpu_quad_matrix_track_particles_launch_kernel(
      d_final_0, d_final_1, d_final_2, d_final_3, d_final_4, d_final_5,
      d_initial_0, d_initial_1, d_initial_2, d_initial_3, d_initial_4,
      d_initial_5, n_part);

  log_exit("gpu quad matrix track particles");
}
  
__global__ void gpu_hexpole_matrix_track_particles_CRT_Constant_Pragma_kernel(
    double *d_final_0, double *d_final_1, double *d_final_2, double *d_final_3,
    double *d_final_4, double *d_final_5, double *d_initial_0,
    double *d_initial_1, double *d_initial_2, double *d_initial_3,
    double *d_initial_4, double *d_initial_5, int n_part){
 
  int tx = threadIdx.x;
  int bx = blockIdx.x;
  int tid = bx*blockDim.x + tx;

  double t_initial[6];
  double t_final[6];
  for(; tid < n_part; tid += gridDim.x*blockDim.x){

    t_initial[0] = d_initial_0[tid];
    t_initial[1] = d_initial_1[tid];
    t_initial[2] = d_initial_2[tid];
    t_initial[3] = d_initial_3[tid];
    t_initial[4] = d_initial_4[tid];
    t_initial[5] = d_initial_5[tid];
  
#pragma unroll 6
    for(int i=0; i<6; i++){
      t_final[i]= constData_C[i];
    }

    double alpha, alpha1, alpha2;
    int ii, iij, iijk;
#pragma unroll 6
    for(int i=0; i<6; i++){
      ii = i*6;
      alpha = 0.0;

#pragma unroll 6
      for(int j=0; j<6; j++){
        iij = ii*6 + j*6;
        alpha1  = 0.0;
        alpha  += constData_R[ii+j]*t_initial[j];

        // The the parts of the matrices which are unused in the CPU implementation
        // are set to zero in this implementation enabling pragma unroll
#pragma unroll 6
        for(int k=0; k<6; ++k) {
          iijk = iij*6 + k*6;
          alpha2  = 0.0;
          alpha1 += constData_T[iij + k]*t_initial[k];

#pragma unroll 6
          for(int l=0; l<6; ++l) {
            alpha2 += constData_Q[iijk + l]*t_initial[l];
          } // for k

          alpha1 += alpha2*t_initial[k];
        } // for k

        alpha += alpha1*t_initial[j];
      } //for j

      t_final[i] += alpha;
    } // for i
    d_final_0[tid] = t_final[0];
    d_final_1[tid] = t_final[1];
    d_final_2[tid] = t_final[2];
    d_final_3[tid] = t_final[3];
    d_final_4[tid] = t_final[4];
    d_final_5[tid] = t_final[5];
  } // for tid
  
}

void gpu_hexpole_matrix_track_particles_launch_kernel(double *d_final_0,
    double *d_final_1, double *d_final_2, double *d_final_3,
    double *d_final_4, double *d_final_5, double *d_initial_0,
    double *d_initial_1, double *d_initial_2, double *d_initial_3,
    double *d_initial_4, double *d_initial_5, int n_part){

  struct GPUBASE* gpuBase = getGpuBase();
  unsigned int nTx = gpuBase->nThreads;
  unsigned int nBx = (n_part + nTx - 1) / nTx;
  if(nBx > gpuBase->maxBlocks) nBx = gpuBase->maxBlocks;
  dim3 dimBlock(nTx,1,1);
  dim3 dimGrid(nBx,1,1);

  gpu_hexpole_matrix_track_particles_CRT_Constant_Pragma_kernel
    <<<dimGrid,dimBlock>>>(d_final_0, d_final_1, d_final_2, d_final_3,
        d_final_4, d_final_5, d_initial_0,d_initial_1, d_initial_2,
        d_initial_3, d_initial_4, d_initial_5, n_part);

  gpuErrorHandler("gpu_hexpole_matrix_track_particles_launch_kernel");
}

void gpu_hexpole_matrix_track_particles(VMATRIX *M, unsigned int n_part){
  struct GPUBASE* gpubase = getGpuBase();
  double* d_particles_initial = gpubase->d_particles;
  double* d_particles_final = gpubase->d_particles;
  unsigned int gpu_pitch = gpubase->gpu_array_pitch;
  double h_C[6];
  double h_R[6*6];
  double h_T[6*6*6];
  double h_Q[6*6*6*6];
  log_entry("gpu hexpole matrix track particles");

  // Fill lower diagonal portion of the matrix
  for(unsigned int i=0; i < 6; i++){
    h_C[i] = M->C[i];
    //    printf("h_C[%d] = %g\n",i,h_C[i]);
    for(unsigned int j=0; j< 6; j++){
      h_R[i*6 + j] = M->R[i][j];
      //      printf("h_R[%d] = %g\n",i*6+j,h_C[i*6+j]);
      for(unsigned int k =0; k<= j; k++){
        h_T[i*6*6 + j*6 + k] = M->T[i][j][k];
        //  printf("h_T[%d] = %g\n",i*6*6 + j*6 +k, h_T[i*6*6 + j*6 + k]);
        for(unsigned int l =0; l<= k; l++){
          h_Q[i*6*6*6 + j*6*6 + k*6 + l] = M->Q[i][j][k][l];
          //  printf("h_T[%d] = %g\n",i*6*6 + j*6 +k, h_T[i*6*6 + j*6 + k]);
        }
        // Fill remaining parts corresponding to this k of Q matrix with zeros
        for(unsigned int l=(k+1); l<6; l++){
          h_Q[i*6*6*6 + j*6*6 + k*6 + l] = 0.0;
        }
      }
      // Fill remaining parts of T & Q matrices with zeros
      for(unsigned int k=(j+1); k<6; k++){
        h_T[i*6*6 + j*6 + k] = 0.0;
        for(unsigned int l=0; l<6; l++){
          h_Q[i*6*6*6 + j*6*6 + k*6 + l] = 0.0;
        }
      }
    }
  }

  gpu_matrix_copyCRTtoGPU(h_C, h_R, h_T, h_Q);

  double* d_initial_0 = d_particles_initial + 0 * gpu_pitch;
  double* d_initial_1 = d_particles_initial + 1 * gpu_pitch;
  double* d_initial_2 = d_particles_initial + 2 * gpu_pitch;
  double* d_initial_3 = d_particles_initial + 3 * gpu_pitch;
  double* d_initial_4 = d_particles_initial + 4 * gpu_pitch;
  double* d_initial_5 = d_particles_initial + 5 * gpu_pitch;

  double* d_final_0 = d_particles_final + 0 * gpu_pitch;
  double* d_final_1 = d_particles_final + 1 * gpu_pitch;
  double* d_final_2 = d_particles_final + 2 * gpu_pitch;
  double* d_final_3 = d_particles_final + 3 * gpu_pitch;
  double* d_final_4 = d_particles_final + 4 * gpu_pitch;
  double* d_final_5 = d_particles_final + 5 * gpu_pitch;
  
  gpu_hexpole_matrix_track_particles_launch_kernel(
      d_final_0, d_final_1, d_final_2, d_final_3, d_final_4, d_final_5,
      d_initial_0, d_initial_1, d_initial_2, d_initial_3, d_initial_4,
      d_initial_5, n_part);

  log_exit("gpu hexpole matrix track particles");
}
 
extern "C" {

void gpu_track_particles(VMATRIX *M, unsigned int n_part) {
  switch (M->order) {
    case 1:
      gpu_dipole_matrix_track_particles(M, n_part);
      return;
    case 2:
      gpu_quad_matrix_track_particles(M, n_part);
      return;
    case 3:
      gpu_hexpole_matrix_track_particles(M, n_part);
      return;
    default:
      fprintf(stdout, "invalid order: %ld  (gpu_track_particle)\n", M->order);
      fflush(stdout);
      exitElegant(1);
      break;
  }
}

} // extern "C"
