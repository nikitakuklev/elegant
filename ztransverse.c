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
  double *posIfreq, *Vfreq, *iZ, Vinterp;
  long ip, ib, nb, n_binned, nfreq, iReal, iImag, plane, offset;
  double factor, tmin, tmax, tmean, dt, dt1, P, dgam, gam, frac;
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

  tmin = DBL_MAX;
  tmean = 0;
  for (ip=0; ip<np; ip++) {
    P = Po*(part[ip][5]+1);
    time[ip] = part[ip][4]*sqrt(sqr(P)+1)/(c_mks*P);
    if (time[ip]<tmin)
      tmin = time[ip];
    tmean += time[ip];
  }
  tmin -= dt;
  tmax = tmin + ztransverse->bin_size*ztransverse->n_bins;
  tmean /= np;

  for (ib=0; ib<ztransverse->n_bins; ib++)
    posItime[2*ib][0] = posItime[2*ib+1][0] = 
      posItime[2*ib][1] = posItime[2*ib+1][1] = 0;

  /* make arrays of I(t)*x and I(t)*y */
  for (ip=n_binned=0; ip<np; ip++) {
    pbin[ip] = -1;
    ib = (time[ip]-tmin)/dt;
    if (ib<0 || ib>nb - 1)
      continue;
    posItime[0][ib] += part[ip][0]-ztransverse->dx;
    posItime[1][ib] += part[ip][2]-ztransverse->dy;
    pbin[ip] = ib;
    pz[ip] = Po*(1+part[ip][5])/sqrt(1+sqr(part[ip][1])+sqr(part[ip][3]));
    n_binned++;
  }
  if (n_binned!=np)
    fprintf(stderr, "warning: only %ld of %ld particles where binned (ZTRANSVERSE)\n", n_binned, np);

  for (plane=0; plane<1; plane++) {
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
    factor = ztransverse->macroParticleCharge/((tmax-tmin)/nb);
    iZ = ztransverse->iZ[plane];
    Vfreq[0] = posIfreq[0]*iZ[0]*factor;
    nfreq = nb/2 + 1;
    if (nb%2==0)
      /* Nyquist term */
      Vfreq[nb-1] = posIfreq[nb-1]*iZ[nb-1]*factor;
    for (ib=1; ib<nfreq-1; ib++) {
      iImag = (iReal = 2*ib-1)+1;
      Vfreq[iReal] = -(posIfreq[iReal]*iZ[iImag] + posIfreq[iImag]*iZ[iReal])*factor; 
      Vfreq[iImag] =  (posIfreq[iReal]*iZ[iReal] - posIfreq[iImag]*iZ[iImag])*factor;
    }

    /* Compute inverse FFT of V(f) to get V(t) */
    realFFT(Vfreq, nb, INVERSE_FFT);
    Vtime = Vfreq;

    /* put zero voltage in Vtime[nb] for use in interpolation */
    Vtime[nb] = 0;
    /* change particle transverse momenta to reflect voltage in relevant bin */
    offset = 2*plane+1;
    for (ip=0; ip<np; ip++) {
      if ((ib=pbin[ip])>=1 && ib<=nb-1) {
        if (ztransverse->interpolate) {
          dt1 = time[ip]-(tmin+dt*ib);
          if (dt1<0)
            dt1 = dt-dt1;
          else 
            ib += 1;
          Vinterp = Vtime[ib-1]+(Vtime[ib]-Vtime[ib-1])/dt*dt1;
        }
        else
          Vinterp = Vtime[ib];
        if (Vinterp)
          part[ip][offset] += Vinterp/(1e6*me_mev)/pz[ip] ;
      }
    }
  }
}


void set_up_ztransverse(ZTRANSVERSE *ztransverse, RUN *run, long pass, long particles, CHARGE *charge)
{
  long i, nfreq;
  double df, t_range;

  if (charge) {
    ztransverse->macroParticleCharge = charge->macroParticleCharge;
    ztransverse->charge = charge->macroParticleCharge*particles;
  } else if (pass==0) {
    ztransverse->macroParticleCharge = 0;
    if (particles)
      ztransverse->macroParticleCharge = ztransverse->charge/particles;
  }
  
  if (ztransverse->initialized)
    return;
  ztransverse->initialized = 1;

  if (ztransverse->bin_size<=0)
    bomb("bin_size must be positive for ZTRANSVERSE element", NULL);
  if (ztransverse->broad_band) {
    double term;
    if (ztransverse->n_bins%2!=0)
      bomb("ZTRANSVERSE element must have n_bins divisible by 2", NULL);
    if (ztransverse->ZxReal || ztransverse->ZxImag ||
        ztransverse->ZyReal || ztransverse->ZyImag )
      bomb("can't specify both broad_band impedance and Z(f) files for ZTRANSVERSE element", NULL);
    nfreq = ztransverse->n_bins/2 + 1;
    ztransverse->iZ[0] = tmalloc(sizeof(**(ztransverse->iZ))*ztransverse->n_bins);
    ztransverse->iZ[1] = tmalloc(sizeof(**(ztransverse->iZ))*ztransverse->n_bins);
    df = 1/(ztransverse->n_bins*ztransverse->bin_size)/(ztransverse->freq);
    /* DC term of iZ is pure real */
    ztransverse->iZ[0][0] = ztransverse->iZ[1][0] = 
      -(ztransverse->Rs)/ztransverse->Q;
    /* Nyquist term--real part of iZ only */
    term = ztransverse->Q*(1.0/(nfreq*df)-nfreq*df);
    ztransverse->iZ[0][ztransverse->n_bins-1] = 
      ztransverse->iZ[1][ztransverse->n_bins-1] = 
        (ztransverse->Rs)*ztransverse->freq/(nfreq*df)/(1-sqr(term))*term;
    for (i=1; i<nfreq-1; i++) {
      term = ztransverse->Q*(1.0/(i*df)-i*df);
      /* imaginary part of i*Z */
      ztransverse->iZ[0][2*i] = 
        ztransverse->iZ[1][2*i] = (ztransverse->Rs)*ztransverse->freq/(i*df)/(1-sqr(term));  
      /* real part is Imag(i*Z)*Q*(fr/f-f/fr) */
      ztransverse->iZ[0][2*i-1] =  
        ztransverse->iZ[1][2*i-1] =  ztransverse->iZ[0][2*i]*term;
    }
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
    if (!SDDS_InitializeInput(&SDDSin, ztransverse->inputFile)) {
      fprintf(stderr, "unable to read file %s\n", ztransverse->inputFile);
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors); 
      exit(1);
    }
    if ((n_spect=SDDS_RowCount(&SDDSin))<4) {
      fprintf(stderr, "too little data in %s\n", ztransverse->inputFile);
      exit(1);
    }
    if (!power_of_2(n_spect-1))
      bomb("number of spectrum points must be 2^n+1, n>1 (ZTRANSVERSE)", NULL);
    ZReal[0] = getTransverseImpedance(&SDDSin, ztransverse->ZxReal);
    ZImag[0] = getTransverseImpedance(&SDDSin, ztransverse->ZxImag);
    ZReal[1] = getTransverseImpedance(&SDDSin, ztransverse->ZyReal);
    ZImag[1] = getTransverseImpedance(&SDDSin, ztransverse->ZyImag);
    if (!(freqData=SDDS_GetColumnInDoubles(&SDDSin, ztransverse->freqColumn)))
      fprintf(stderr, "Unable to read column %s (ZTRANSVERSE)\n", ztransverse->freqColumn);
    if ((df_spect = (freqData[n_spect-1]-freqData[0])/(n_spect-1))<=0) {
      fprintf(stderr, "Zero or negative frequency spacing in %s (ZTRANSVERSE)\n",
              ztransverse->inputFile);
      exit(1);
    }
    t_range = ztransverse->n_bins*ztransverse->bin_size;
    ztransverse->n_bins = 2*(n_spect-1);
    ztransverse->bin_size = 1.0/(ztransverse->n_bins*df_spect);
    if (t_range>ztransverse->n_bins*ztransverse->bin_size) {
      fprintf(stderr, "error for ZTRANSVERSE element:\nimpedance-spectrum-equivalent binning range not sufficient.\n");
      fprintf(stderr, "consider padding the impedance spectrum\n");
      exit(1);
    }
    if (!SDDS_Terminate(&SDDSin)) {
      fprintf(stderr, "Error closing data set %s\n",
              ztransverse->inputFile);
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
        ztransverse->iZ[0][2*i-1] = -ZImag[0][i];
        ztransverse->iZ[1][2*i-1] = -ZImag[1][i];
        ztransverse->iZ[0][2*i  ] = ZReal[0][i];
        ztransverse->iZ[1][2*i  ] = ZReal[1][i];
      }
    }
    free(ZReal[0]);
    free(ZReal[1]);
    free(ZImag[0]);
    free(ZImag[1]);
  }
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
    fprintf(stderr, "Unable to read column %s (ZTRANSVERSE)\n",
            ZName);
    exit(1);
  }
  return Z;
}

