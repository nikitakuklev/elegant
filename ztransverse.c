/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: zlongit.c
 * contents: track_through_zlongit()
 *
 * Michael Borland, 1993
 */
#include "mdb.h"
#include "track.h"
#include "table.h"
#include "fftpackC.h"

#define WAKE_COLUMNS 5
static SDDS_DEFINITION wake_column[WAKE_COLUMNS] = {
    {"Deltat", "&column name=Deltat, symbol=\"$gD$rt\", units=s, type=double, description=\"Time after head of bunch\" &end"},
    {"Ix", "&column name=Ix, symbol=\"<I*x>\", units=C*m/s, type=double, description=\"Transverse horizontal moment\" &end"},
    {"Wx", "&column name=Wx, symbol=\"W$bx$n\", units=V/m, type=double, description=\"Transverse horizontal wake\" &end"},
    {"Iy", "&column name=Iy, symbol=\"<I*y>\", units=C*m/s, type=double, description=\"Transverse vertical moment\" &end"},
    {"Wy", "&column name=Wy, symbol=\"W$by$n\", units=V/m, type=double, description=\"Transverse vertical wake\" &end"},
    };

#define WAKE_PARAMETERS 5
#define BB_WAKE_PARAMETERS 5
#define NBB_WAKE_PARAMETERS 2
static SDDS_DEFINITION wake_parameter[WAKE_PARAMETERS] = {
    {"Pass", "&parameter name=Pass, type=long &end"},
    {"q", "&parameter name=q, units=C, type=double, description=\"Total charge\" &end"},
    {"Rs", "&parameter name=Rs, symbol=\"R$bs$n\", units=\"$gW$r\", type=double, description=\"Broad-band impedance\" &end"},
    {"fo", "&parameter name=fo, symbol=\"f$bo$n\", units=Hz, type=double, description=\"Frequency of BB resonator\" &end"},
    {"Deltaf", "&parameter name=Deltaf, symbol=\"$gD$rf\", units=Hz, type=double, description=\"Frequency sampling interval\" &end"},
    } ;


void set_up_ztransverse(ZTRANSVERSE *ztransverse, RUN *run, long pass, long particles, CHARGE *charge);
double *getTransverseImpedance(SDDS_DATASET *SDDSin, char *ZName);

void track_through_ztransverse(double **part, long np, ZTRANSVERSE *ztransverse, double Po,
                               RUN *run, long i_pass, CHARGE *charge
                               )
{
  static double *posItime[2] = {NULL, NULL};     /* array for particle density times x, y*/
  static double *posIfreq;                       /* array for FFT of particle density times x or y*/
  static double *Vtime = NULL;           /* array for voltage acting on each bin */
  static long max_n_bins = 0;
  static long *pbin = NULL;              /* array to record which bin each particle is in */
  static double *time = NULL;            /* array to record arrival time of each particle */
  static double *pz = NULL;
  static long max_np = 0;
  double *Vfreq, *iZ;
  long ib, nb, n_binned, nfreq, iReal, iImag, plane, first;
  double factor, tmin, tmax, tmean, dt, userFactor[2];
  static long not_first_call = -1;
#if defined(DEBUG)
  FILE *fp;
#endif
  not_first_call += 1;
  
  set_up_ztransverse(ztransverse, run, i_pass, np, charge);
  nb = ztransverse->n_bins;
  dt = ztransverse->bin_size;

  if (nb>max_n_bins) {
    posItime[0] = trealloc(posItime[0], 2*sizeof(**posItime)*(max_n_bins=nb));
    posItime[1] = trealloc(posItime[1], 2*sizeof(**posItime)*(max_n_bins=nb));
    posIfreq = trealloc(posIfreq, 2*sizeof(*posIfreq)*(max_n_bins=nb));
    Vtime = trealloc(Vtime, 2*sizeof(*Vtime)*(max_n_bins+1));
  }

  if (np>max_np) {
    pbin = trealloc(pbin, sizeof(*pbin)*(max_np=np));
    time = trealloc(time, sizeof(*time)*max_np);
    pz = trealloc(pz, sizeof(*pz)*max_np);
  }

  for (ib=0; ib<nb; ib++)
    posItime[0][2*ib] = posItime[0][2*ib+1] = 
      posItime[1][2*ib] = posItime[1][2*ib+1] = 0;

  /* Compute time coordinate of each particle */
  tmean = computeTimeCoordinates(time, Po, part, np);
  find_min_max(&tmin, &tmax, time, np);
  tmin -= dt;
  tmax -= dt;
  
  /* make arrays of I(t)*x and I(t)*y */
  n_binned = binTransverseTimeDistribution(posItime, pz, pbin, tmin, dt, nb, time, part, Po, np,
                                           ztransverse->dx, ztransverse->dy, 1, 1);

  if (n_binned!=np) {
    fprintf(stdout, "Warning: only %ld of %ld particles were binned (ZTRANSVERSE)!\n", n_binned, np);
    if (!not_first_call) {
      fprintf(stdout, "*** This may produce unphysical results.  Your wake needs smaller frequency\n");
      fprintf(stdout, "    spacing to cover a longer time span.\n");
    }
    fflush(stdout);
  }  

  userFactor[0] = ztransverse->factor*ztransverse->xfactor;
  userFactor[1] = ztransverse->factor*ztransverse->yfactor;
  first = 1;
  for (plane=0; plane<2; plane++) {
    if (userFactor[plane]==0) {
      for (ib=0; ib<nb; ib++)
	Vtime[ib] = 0;
    } else {
      if (ztransverse->smoothing)
	SavitzyGolaySmooth(posItime[plane], nb, ztransverse->SGOrder,
			   ztransverse->SGHalfWidth, ztransverse->SGHalfWidth, 0);
      
      /* Take the FFT of (x*I)(t) to get (x*I)(f) */
      memcpy(posIfreq, posItime[plane], 2*nb*sizeof(*posIfreq));
      realFFT(posIfreq, nb, 0);
      
      /* Compute V(f) = i*Z(f)*(x*I)(f), putting in a factor 
       * to normalize the current waveform
       */
      Vfreq = Vtime;
      factor = ztransverse->macroParticleCharge/dt*userFactor[plane];
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

    if (ztransverse->SDDS_wake_initialized && ztransverse->wakes) {
      /* wake potential output */
      factor = ztransverse->macroParticleCharge/dt;
      if (ztransverse->wake_interval<=0 || (i_pass%ztransverse->wake_interval)==0) {
	if (first && !SDDS_StartTable(&ztransverse->SDDS_wake, nb)) {
	  SDDS_SetError("Problem starting SDDS table for wake output (track_through_ztransverse)");
	  SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
	}
	for (ib=0; ib<nb; ib++) {
	  if (!SDDS_SetRowValues(&ztransverse->SDDS_wake, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, ib,
				 0, ib*dt, 
				 1+plane*2, posItime[plane][ib]*factor,  
				 2+plane*2, Vtime[ib], -1)) {
	    SDDS_SetError("Problem setting rows of SDDS table for wake output (track_through_ztransverse)");
	    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
	  }
	}
	if (!SDDS_SetParameters(&ztransverse->SDDS_wake, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,
				"Pass", i_pass, "q", ztransverse->macroParticleCharge*np, NULL)) {
	  SDDS_SetError("Problem setting parameters of SDDS table for wake output (track_through_ztransverse)");
	  SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
	}
	if (ztransverse->broad_band) {
	  if (!SDDS_SetParameters(&ztransverse->SDDS_wake, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,
				  "Rs", ztransverse->Rs, "fo", ztransverse->freq, 
				  "Deltaf", ztransverse->bin_size, NULL)) {
	    SDDS_SetError("Problem setting parameters of SDDS table for wake output (track_through_ztransverse)");
	    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
	  }
	}
	if (!first) {
	  if (!SDDS_WriteTable(&ztransverse->SDDS_wake)) {
	    SDDS_SetError("Problem writing SDDS table for wake output (track_through_ztransverse)");
	    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
	  }
	  SDDS_DoFSync(&ztransverse->SDDS_wake);
	}
      }
    }
    first = 0;
  }

#if defined(MIMIMIZE_MEMORY)
  free(posItime[0]);
  free(posItime[1]);
  free(posIfreq);
  free(Vtime);
  free(pbin);
  free(time);
  free(pz);
  posItime[0] = posItime[1] = Vtime = time = pz = posIfreq = NULL;
  pbin = NULL;
  max_n_bins = max_np = 0 ;
#endif

}


void set_up_ztransverse(ZTRANSVERSE *ztransverse, RUN *run, long pass, long particles, CHARGE *charge)
{
  long i, nfreq;
  double df;

  if (charge) {
    ztransverse->macroParticleCharge = charge->macroParticleCharge;
  } else if (pass==0) {
    ztransverse->macroParticleCharge = 0;
    if (particles)
      ztransverse->macroParticleCharge = ztransverse->charge/particles;
  }

  if (ztransverse->initialized)
    return;

  ztransverse->SDDS_wake_initialized = 0;

  if (ztransverse->broad_band) {
    /* Use impedance Z = -i*wr/w*Rs/(1 + i*Q(w/wr-wr/w))
       */
    double term;
    if (ztransverse->bin_size<=0)
      bomb("bin_size must be positive for ZTRANSVERSE element", NULL);
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
    /* DC term of iZ is 0  */
    ztransverse->iZ[0][0] = ztransverse->iZ[1][0] = 0;
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
    /* Nyquist term--real part of iZ only */
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
    if (!ztransverse->freqColumn || !ztransverse->inputFile)
      bomb("you must give an inputFile and freqColumn, or use a broad band model (ZTRANSVERSE)", NULL);
    if (!ztransverse->ZxReal && !ztransverse->ZxImag &&
        !ztransverse->ZyReal && !ztransverse->ZxImag)
      bomb("you must either give broad_band=1, or Z[xy]Real and/or Z[xy]Imag (ZTRANSVERSE)", NULL);
    if (!SDDS_InitializeInputFromSearchPath(&SDDSin, ztransverse->inputFile) || !SDDS_ReadPage(&SDDSin)) {
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
    if (!(freqData=SDDS_GetColumnInDoubles(&SDDSin, ztransverse->freqColumn))) {
      fprintf(stdout, "Unable to read column %s (ZTRANSVERSE)\n", ztransverse->freqColumn);
      fflush(stdout);
      exit(1);
    }
    if (!checkPointSpacing(freqData, n_spect, 1e-6)) {
      fprintf(stdout, "Frequency values are not equispaced (ZTRANSVERSE)\n");
      fflush(stdout);
      exit(1);
    }
    if ((df_spect = (freqData[n_spect-1]-freqData[0])/(n_spect-1))<=0) {
      fprintf(stdout, "Zero or negative frequency spacing in %s (ZTRANSVERSE)\n",
              ztransverse->inputFile);
      fflush(stdout);
      exit(1);
    }
    df = df_spect;
    nfreq = n_spect;
    ztransverse->n_bins = 2*(n_spect-1);
    ztransverse->bin_size = 1.0/(ztransverse->n_bins*df_spect);
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

  if (ztransverse->wakes) {
    ztransverse->wakes = compose_filename(ztransverse->wakes, run->rootname);
    if (ztransverse->broad_band) 
      SDDS_ElegantOutputSetup(&ztransverse->SDDS_wake, ztransverse->wakes, SDDS_BINARY, 
			      1, "transverse wake",
			      run->runfile, run->lattice, wake_parameter, BB_WAKE_PARAMETERS,
			      wake_column, WAKE_COLUMNS, "set_up_ztransverse", 
			      SDDS_EOS_NEWFILE|SDDS_EOS_COMPLETE);
    else {
      SDDS_ElegantOutputSetup(&ztransverse->SDDS_wake, ztransverse->wakes, SDDS_BINARY, 
			      1, "transverse wake",
			      run->runfile, run->lattice, wake_parameter, NBB_WAKE_PARAMETERS,
			      wake_column, WAKE_COLUMNS, "set_up_ztransverse", 
			      SDDS_EOS_NEWFILE|SDDS_EOS_COMPLETE);
    }
    ztransverse->SDDS_wake_initialized = 1;
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
    fprintf(stdout, "Unable to read column %s (ZTRANSVERSE)\n",
            ZName);
    fflush(stdout);
    exit(1);
  }
  return Z;
}

