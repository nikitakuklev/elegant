/* Copyright 1994 by Michael Borland and Argonne National Laboratory,
 * all rights reserved.
 */
/* file: zlongit.c
 * contents: track_through_zlongit()
 *
 * Michael Borland, 1993
 */
#include "mdb.h"
#include "track.h"
#include "table.h"
#include "fftpackC.h"

void set_up_ztransverse(ZTRANSVERSE *ztransverse, RUN *run, long pass, long particles, CHARGE *charge);
double *getTransverseImpedance(SDDS_DATASET *SDDSin, char *ZName);

void track_through_ztransverse(double **part, long np, ZTRANSVERSE *ztransverse, double Po,
                               RUN *run, long i_pass, CHARGE *charge
                               )
{
  static double *posItime[2] = {NULL, NULL};     /* array for histogram of particle density times x, y*/
  static double *Vtime = NULL;           /* array for voltage acting on each bin */
  static long max_n_bins = 0;
  static long *pbin = NULL;              /* array to record which bin each particle is in */
  static double *time = NULL;            /* array to record arrival time of each particle */
  static double *pz = NULL;
  static long max_np = 0;
  double *posIfreq, *Vfreq, *iZ;
  long ib, nb, n_binned, nfreq, iReal, iImag, plane;
  double factor, tmin, tmax, tmean, dt;
#if defined(DEBUG)
  static long first_call = 1;
  FILE *fp;
#endif

  set_up_ztransverse(ztransverse, run, i_pass, np, charge);
  nb = ztransverse->n_bins;
  dt = ztransverse->bin_size;

  if (ztransverse->n_bins>max_n_bins) {
    posItime[0] = trealloc(posItime[0], 2*sizeof(**posItime)*(max_n_bins=ztransverse->n_bins));
    posItime[1] = trealloc(posItime[1], 2*sizeof(**posItime)*(max_n_bins=ztransverse->n_bins));
    Vtime = trealloc(Vtime, 2*sizeof(*Vtime)*(max_n_bins+1));
  }

  if (np>max_np) {
    pbin = trealloc(pbin, sizeof(*pbin)*(max_np=np));
    time = trealloc(time, sizeof(*time)*max_np);
    pz = trealloc(pz, sizeof(*pz)*max_np);
  }

  for (ib=0; ib<ztransverse->n_bins; ib++)
    posItime[0][2*ib] = posItime[0][2*ib+1] = 
      posItime[1][2*ib] = posItime[1][2*ib+1] = 0;

  /* Compute time coordinate of each particle */
  tmean = computeTimeCoordinates(time, Po, part, np);
  find_min_max(&tmin, &tmax, time, np);
  tmin -= dt;
  tmax -= dt;
  
  /* make arrays of I(t)*x and I(t)*y */
  n_binned = binTransverseTimeDistribution(posItime, pz, pbin, tmin, dt, nb, time, part, Po, np,
                                           ztransverse->dx, ztransverse->dy);

  if (n_binned!=np)
    fprintf(stdout, "warning: only %ld of %ld particles where binned (ZTRANSVERSE)\n", n_binned, np);
    fflush(stdout);

  for (plane=0; plane<2; plane++) {
    if (ztransverse->smoothing)
      SavitzyGolaySmooth(posItime[plane], nb, ztransverse->SGOrder,
                         ztransverse->SGHalfWidth, ztransverse->SGHalfWidth, 0);

    /* Take the FFT of (x*I)(t) to get (x*I)(f) */
    realFFT(posItime[plane], nb, 0);
    posIfreq = posItime[plane]; 
    
    /* Compute V(f) = i*Z(f)*(x*I)(f), putting in a factor 
     * to normalize the current waveform
     */
    Vfreq = Vtime;
    factor = ztransverse->macroParticleCharge/dt;
    iZ = ztransverse->iZ[plane];
    Vfreq[0] = posIfreq[0]*iZ[0]*factor;
    nfreq = nb/2 + 1;
    if (nb%2==0)
      /* Nyquist term */
      Vfreq[nb-1] = posIfreq[nb-1]*iZ[nb-1]*factor;
    for (ib=1; ib<nfreq-1; ib++) {
      iImag = (iReal = 2*ib-1)+1;
      /* The signs are chosen here to get agreement with TRFMODE.
         In particular, test particles following closely behind the 
         drive particle get defocused.
         */
      Vfreq[iReal] =  (posIfreq[iReal]*iZ[iImag] + posIfreq[iImag]*iZ[iReal])*factor; 
      Vfreq[iImag] = -(posIfreq[iReal]*iZ[iReal] - posIfreq[iImag]*iZ[iImag])*factor;
    }

    /* Compute inverse FFT of V(f) to get V(t) */
    realFFT(Vfreq, nb, INVERSE_FFT);
    Vtime = Vfreq;
            
    /* change particle transverse momenta to reflect voltage in relevant bin */
    applyTransverseWakeKicks(part, time, pz, pbin, np, 
                             Po, plane, 
                             Vtime, nb, tmin, dt, ztransverse->interpolate);
    
  }

#if defined(MIMIMIZE_MEMORY)
  free(posItime[0]);
  free(posItime[1]);
  free(Vtime);
  free(pbin);
  free(time);
  free(pz);
  posItime[0] = posItime[1] = Vtime = time = pz = NULL;
  pbin = NULL;
  max_n_bins = max_np = 0 ;
#endif

}


void set_up_ztransverse(ZTRANSVERSE *ztransverse, RUN *run, long pass, long particles, CHARGE *charge)
{
  long i, nfreq;
  double df, t_range;

  if (charge) {
    ztransverse->macroParticleCharge = charge->macroParticleCharge;
  } else if (pass==0) {
    ztransverse->macroParticleCharge = 0;
    if (particles)
      ztransverse->macroParticleCharge = ztransverse->charge/particles;
  }
  
  if (ztransverse->initialized)
    return;

  if (ztransverse->bin_size<=0)
    bomb("bin_size must be positive for ZTRANSVERSE element", NULL);
  if (ztransverse->broad_band) {
    /* Use impedance Z = j*wr/w*Rs/(1 + j*Q(w/wr-wr/w))
       */
    double term;
    if (ztransverse->n_bins%2!=0)
      bomb("ZTRANSVERSE element must have n_bins divisible by 2", NULL);
    if (ztransverse->ZxReal || ztransverse->ZxImag ||
        ztransverse->ZyReal || ztransverse->ZyImag )
      bomb("can't specify both broad_band impedance and Z(f) files for ZTRANSVERSE element", NULL);
    nfreq = ztransverse->n_bins/2 + 1;
    ztransverse->iZ[0] = tmalloc(sizeof(**(ztransverse->iZ))*ztransverse->n_bins);
    ztransverse->iZ[1] = tmalloc(sizeof(**(ztransverse->iZ))*ztransverse->n_bins);
    /* df is the frequency spacing normalized to the resonant frequency */
    df = 1/(ztransverse->n_bins*ztransverse->bin_size)/(ztransverse->freq);
    /* DC term of iZ is zero */
    ztransverse->iZ[0][0] = ztransverse->iZ[1][0] = 0;
    /* Nyquist term--real part of iZ only */
    for (i=1; i<nfreq-1; i++) {
      term = ztransverse->Q*(i*df-1.0/(i*df));
      /* real part of i*Z */
      ztransverse->iZ[0][2*i-1] =  
        ztransverse->iZ[1][2*i-1] =  
          ztransverse->Rs/(i*df)/(1+term*term);
      /* imaginary part of i*Z is -Real[i*Z]*term */
      ztransverse->iZ[0][2*i] = 
        ztransverse->iZ[1][2*i] = 
          -term*ztransverse->iZ[0][2*i-1];
    }
    term = ztransverse->Q*(1.0/(nfreq*df)-nfreq*df);
    ztransverse->iZ[0][ztransverse->n_bins-1] = 
      ztransverse->iZ[1][ztransverse->n_bins-1] = 
        ztransverse->Rs/(nfreq*df)/(1+term*term);
    df *= ztransverse->freq;
  } else {
    double *ZReal[2], *ZImag[2], *freqData;
    double df_spect;
    long n_spect;
    SDDS_DATASET SDDSin;
    if (ztransverse->n_bins<1)
      bomb("ZTRANSVERSE element must have n_bins>=1", NULL);
    if (!ztransverse->freqColumn || !ztransverse->inputFile)
      bomb("you must give an inputFile and freqColumn, or use a broad band model (ZTRANSVERSE)", NULL);
    if (!ztransverse->ZxReal && !ztransverse->ZxImag &&
        !ztransverse->ZyReal && !ztransverse->ZxImag)
      bomb("you must either give broad_band=1, or Z[xy]Real and/or Z[xy]Imag (ZTRANSVERSE)", NULL);
    if (!SDDS_InitializeInput(&SDDSin, ztransverse->inputFile) || !SDDS_ReadPage(&SDDSin)) {
      fprintf(stdout, "unable to read file %s\n", ztransverse->inputFile);
      fflush(stdout);
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors); 
      exit(1);
    }
    if ((n_spect=SDDS_RowCount(&SDDSin))<4) {
      fprintf(stdout, "too little data in %s\n", ztransverse->inputFile);
      fflush(stdout);
      exit(1);
    }
    if (!power_of_2(n_spect-1))
      bomb("number of spectrum points must be 2^n+1, n>1 (ZTRANSVERSE)", NULL);
    ZReal[0] = getTransverseImpedance(&SDDSin, ztransverse->ZxReal);
    ZImag[0] = getTransverseImpedance(&SDDSin, ztransverse->ZxImag);
    ZReal[1] = getTransverseImpedance(&SDDSin, ztransverse->ZyReal);
    ZImag[1] = getTransverseImpedance(&SDDSin, ztransverse->ZyImag);
    if (!(freqData=SDDS_GetColumnInDoubles(&SDDSin, ztransverse->freqColumn)))
      fprintf(stdout, "Unable to read column %s (ZTRANSVERSE)\n", ztransverse->freqColumn);
      fflush(stdout);
    if ((df_spect = (freqData[n_spect-1]-freqData[0])/(n_spect-1))<=0) {
      fprintf(stdout, "Zero or negative frequency spacing in %s (ZTRANSVERSE)\n",
              ztransverse->inputFile);
      fflush(stdout);
      exit(1);
    }
    df = df_spect;
    nfreq = n_spect;
    t_range = ztransverse->n_bins*ztransverse->bin_size;
    ztransverse->n_bins = 2*(n_spect-1);
    ztransverse->bin_size = 1.0/(ztransverse->n_bins*df_spect);
    if (t_range>ztransverse->n_bins*ztransverse->bin_size) {
      fprintf(stdout, "error for ZTRANSVERSE element:\nimpedance-spectrum-equivalent binning range not sufficient.\n");
      fflush(stdout);
      fprintf(stdout, "consider padding the impedance spectrum\n");
      fflush(stdout);
      exit(1);
    }
    if (!SDDS_Terminate(&SDDSin)) {
      fprintf(stdout, "Error closing data set %s\n",
              ztransverse->inputFile);
      fflush(stdout);
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      exit(1);
    }
    if (!(ztransverse->iZ[0] =
          calloc(sizeof(*ztransverse->iZ[0]), n_spect*2)) ||
        !(ztransverse->iZ[1] =
          calloc(sizeof(*ztransverse->iZ[1]), n_spect*2)))
      bomb("memory allocation failure (ZTRANSVERSE)", NULL);
    for (i=0; i<n_spect; i++) {
      if (i==0) {
        /* DC term */
        ztransverse->iZ[0][i] = -ZImag[0][i];
        ztransverse->iZ[1][i] = -ZImag[1][i];
      } else if (i==n_spect-1 && ztransverse->n_bins%2==0) {
        /* Nyquist term */
        ztransverse->iZ[0][n_spect-1] = -ZImag[0][i];
        ztransverse->iZ[1][n_spect-1] = -ZImag[1][i];
      } else {
        /* real part of iZ */
        ztransverse->iZ[0][2*i-1] = -ZImag[0][i];
        ztransverse->iZ[1][2*i-1] = -ZImag[1][i];
        /* imaginary part of iZ */
        ztransverse->iZ[0][2*i  ] = ZReal[0][i];
        ztransverse->iZ[1][2*i  ] = ZReal[1][i];
      }
    }
    free(ZReal[0]);
    free(ZReal[1]);
    free(ZImag[0]);
    free(ZImag[1]);
  }
#if 0
  if (!ztransverse->initialized) {
      FILE *fp;
      fp = fopen_e("ztransverse.sdds", "w", 0);
      fprintf(fp, "SDDS1\n&column name=f units=Hz type=double &end\n");
      fprintf(fp, "&column name=ZReal type=double &end\n");
      fprintf(fp, "&column name=ZImag type=double &end\n");
      fprintf(fp, "&data mode=ascii no_row_counts=1 &end\n");
      for (i=0; i<nfreq; i++) 
        fprintf(fp, "%21.15e %21.15e %21.15e\n",
                i*df, ztransverse->iZ[0][2*i], -ztransverse->iZ[0][2*i-1]);
      fclose(fp);
    }
#endif

  ztransverse->initialized = 1;
}

double *getTransverseImpedance(SDDS_DATASET *SDDSin, 
                               char *ZName)
{
  long rows;
  double *Z;
  rows = SDDS_RowCount(SDDSin);
  if (!ZName || !strlen(ZName)) {
    Z = calloc(sizeof(*Z), rows);
    return Z;
  }
  if (!(Z=SDDS_GetColumnInDoubles(SDDSin, ZName))) {
    fprintf(stdout, "Unable to read column %s (ZTRANSVERSE)\n",
            ZName);
    fflush(stdout);
    exit(1);
  }
  return Z;
}

