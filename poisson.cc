#include <iostream>
#include <fstream>
#include <complex>
#include <algorithm>
#include <string>
#include <time.h>
#include <iomanip>

#include "poisson.h"

#include <fftw3.h>
#include <mdb.h>

using namespace std;


const double pi = 3.14159265358979323846;


void meshgrid_wavenumber(int N_x, int N_y, double x_domain, double y_domain, double **k_square) {
  //vector<double> v_1d_x(N_x, 0.0);
  double *v_1d_x;
  v_1d_x = (double *)tmalloc(sizeof(*v_1d_x) * N_x);
  
  for ( int i=0; i<N_x; i++ ) {
    if (i<N_x/2) {
      v_1d_x[i] = i;
    } 
    else {
      v_1d_x[i] = -(N_x-i);
    }
  } 
  
  //vector<double> v_1d_y(N_y, 0.0);
  double* v_1d_y;
  v_1d_y = (double *)tmalloc(sizeof(*v_1d_y) * N_y);
  for ( int i=0; i<N_y; i++ ) {
    if (i<N_y/2) {
      v_1d_y[i] = i;
    } 
    else {
      v_1d_y[i] = -(N_y-i);
    }
  } 

  //vector<vector<double> > kx(N_x, vector<double> (N_y, 0));
  //vector<vector<double> > ky(N_x, vector<double> (N_y, 0));

  double **kx, **ky;
  kx = (double **)czarray_2d(sizeof(double), N_x, N_y);
  ky = (double **)czarray_2d(sizeof(double), N_x, N_y);
  
  for ( int i=0; i<N_x; i++ ) {
    for ( int j=0; j<N_y; j++ ) {
      kx[i][j] = v_1d_x[i] / x_domain;
      ky[i][j] = v_1d_y[j] / y_domain;
    }
  }
  
  
  //vector<vector<double> > k_square(N_x, vector<double> (N_y, 0));
  //k_square = (double **)czarray_2d(sizeof(double), N_x, N_y);
  for ( int i=0; i<N_x; i++ ) {
    for ( int j=0; j<N_y; j++ ) {
      k_square[i][j] = -(kx[i][j]*kx[i][j] + ky[i][j]*ky[i][j]) * (4.0*pi*pi);
    }
  }
  

  k_square[0][0] = -1.0;
  free(v_1d_x);
  free(v_1d_y);
  free_czarray_2d((void **)kx, N_x, N_y);
  free_czarray_2d((void **)ky, N_x, N_y);
  //free(kx);
  //free(ky);
  //return k_square;
}



void subtract_const(int N_x, int N_y, double **vec_c) {
  double c00;

  //vector<vector<double> > vec_x(N_x, vector<double> (N_y, 0));
  c00 = vec_c[0][0];
  for ( int i=0; i<N_x; i++ ) {
    for ( int j=0; j<N_y; j++ ) {
      vec_c[i][j] = (vec_c[i][j] - c00) / N_x / N_y;
    }
  }
  //return vec_x;
}



void poisson_solver(double **vec_rhs, double x_domain, double y_domain, int N_x, int N_y, double **u_fft) {

  // Save the mesh grid for re-use, since it is expensive to generate.
  static double **k_square = NULL;
  static double last_x_domain, last_y_domain;
  static int last_N_x, last_N_y;
  int parametersChanged = 1;

  if (N_x != last_N_x || N_y != last_N_y || x_domain != last_x_domain || y_domain != last_y_domain)
    parametersChanged = 1;
  else
    parametersChanged = 0;
  if (parametersChanged || !k_square) {
    if (k_square)
      free_czarray_2d((void **)k_square, N_x, N_y);
    k_square = (double **)czarray_2d(sizeof(double), N_x, N_y);
    meshgrid_wavenumber(N_x, N_y, x_domain, y_domain, k_square);
  }
  
  static fftw_complex *fftwIn = NULL;
  if (!fftwIn || parametersChanged) {
    int n0 = N_x * (N_y / 2 + 1);
    if (fftwIn)
      fftw_free(fftwIn);
    fftwIn = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n0);
  }
  //fftwOut = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N_x * N_y);
    
  /*
  for (int i=0; i<N_x; i++) {
     for (int j=0; j<N_y; j++) {
       fftwIn[i*N_y+j][0] = vec_rhs[i][j];
       fftwIn[i*N_y+j][1] = 0;
     }
  }
  */

  static short p1Created = 0;
  static fftw_plan p1;
  static double **p1Buffer = NULL;
  if (!p1Created || parametersChanged) {
    p1Buffer = (double**)czarray_2d(sizeof(double), N_x, N_y);
    memcpy(&(p1Buffer[0][0]), &(vec_rhs[0][0]), sizeof(double)*N_x*N_y);
    p1 = fftw_plan_dft_r2c_2d(N_x, N_y, &(p1Buffer[0][0]), fftwIn, FFTW_MEASURE);
    p1Created = 1;
    printf("FFTW_MEASURE plan 1 created\n"); fflush(stdout);
  }
  memcpy(&(p1Buffer[0][0]), &(vec_rhs[0][0]), sizeof(double)*N_x*N_y);

  fftw_execute(p1);
  
  divide_fftw_vec(N_x, N_y, fftwIn, k_square);
  
  //p = fftw_plan_dft_2d(N_x, N_y, fftwIn, fftwIn,
  //		       FFTW_BACKWARD, FFTW_MEASURE);


  static short p2Created = 0;
  static fftw_plan p2;
  static double **p2Buffer = NULL;
  if (!p2Created || parametersChanged) {
    p2Buffer = (double**)czarray_2d(sizeof(double), N_x, N_y);
    fftw_complex *fftwTmp = NULL;
    int n0 = N_x * (N_y / 2 + 1);
    fftwTmp = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n0);
    memcpy(&(fftwTmp[0][0]), &(fftwIn[0][0]), sizeof(fftw_complex)*n0);
    p2 = fftw_plan_dft_c2r_2d(N_x, N_y, fftwIn, &(p2Buffer[0][0]), FFTW_MEASURE);
    p2Created = 1;
    memcpy(&(fftwIn[0][0]), &(fftwTmp[0][0]), sizeof(fftw_complex)*n0);
    fftw_free(fftwTmp);
    printf("FFTW_MEASURE plan 2 created\n"); fflush(stdout);
  }

  fftw_execute(p2);
  memcpy(&(u_fft[0][0]), &(p2Buffer[0][0]), sizeof(double)*N_x*N_y);

  /*
  for (int i=0; i<N_x; i++) {
    for (int j=0; j<N_y; j++) {
       u_fft[i][j] = fftwIn[i*N_y+j][0];
     }
  }
  */

  //fftw_free(fftwOut);
  //fftw_free(fftwOut2);
  
  subtract_const(N_x, N_y, u_fft);
  
  last_N_x = N_x;
  last_N_y = N_y;
  last_x_domain = x_domain;
  last_y_domain = y_domain;
  
  //return u_fft;
}


void divide_fftw_vec(int N_x, int N_y, fftw_complex *fftwIn, double **vec_x) {

  int jmax = N_y / 2 + 1;
   for ( int i=0; i<N_x; i++ ) {
     for ( int j=0; j<jmax; j++ ) {
       fftwIn[i*jmax+j][0] = fftwIn[i*jmax+j][0] / vec_x[i][j];
       fftwIn[i*jmax+j][1] = fftwIn[i*jmax+j][1] / vec_x[i][j];
     }
   }
   return;
}




