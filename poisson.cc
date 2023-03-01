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

  //vector<vector<double> > vec_x(N_x, vector<double> (N_y, 0));
  for ( int i=0; i<N_x; i++ ) {
    for ( int j=0; j<N_y; j++ ) {
      vec_c[i][j] = (vec_c[i][j] - vec_c[0][0]) / N_x / N_y;
    }
  }
  //return vec_x;
}



void poisson_solver(double **vec_rhs, double x_domain, double y_domain, int N_x, int N_y, double **u_fft) {

  
  double **k_square;
  k_square = (double **)czarray_2d(sizeof(double), N_x, N_y);
  meshgrid_wavenumber(N_x, N_y, x_domain, y_domain, k_square);
  
  fftw_complex *fftwIn;
  fftw_plan p;
  int n0 = N_x * (N_y / 2 + 1);
  fftwIn = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n0);
  //fftwOut = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N_x * N_y);

  /*
  for (int i=0; i<N_x; i++) {
     for (int j=0; j<N_y; j++) {
       fftwIn[i*N_y+j][0] = vec_rhs[i][j];
       fftwIn[i*N_y+j][1] = 0;
     }
  }
  */
  
  p = fftw_plan_dft_r2c_2d(N_x, N_y, &(vec_rhs[0][0]), fftwIn, FFTW_ESTIMATE);
  
 
  // p = fftw_plan_dft_2d(N_x, N_y, fftwIn, fftwIn,
  //		       FFTW_FORWARD, FFTW_ESTIMATE);
  
  fftw_execute(p);
  
  divide_fftw_vec(N_x, N_y, fftwIn, k_square);
  
  //p = fftw_plan_dft_2d(N_x, N_y, fftwIn, fftwIn,
  //		       FFTW_BACKWARD, FFTW_ESTIMATE);

  p = fftw_plan_dft_c2r_2d(N_x, N_y, fftwIn, &(u_fft[0][0]), FFTW_ESTIMATE);
  
  fftw_execute(p);

  /*
  for (int i=0; i<N_x; i++) {
    for (int j=0; j<N_y; j++) {
       u_fft[i][j] = fftwIn[i*N_y+j][0];
     }
  }
  */

  fftw_destroy_plan(p);
  fftw_free(fftwIn);
  free_czarray_2d((void **)k_square, N_x, N_y);
  //free(k_square);
  //fftw_free(fftwOut);
  //fftw_free(fftwOut2);
  
  subtract_const(N_x, N_y, u_fft);
  
  
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




