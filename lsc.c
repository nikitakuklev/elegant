/*************************************************************************\
* Copyright (c) 2003 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2003 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: lsc.c 
 *
 * Michael Borland, 2003
 */
#include "mdb.h"
#include "track.h"

#include "fftpackC.h"

void track_through_lscdrift(double **part, long np, LSCDRIFT *LSC, double Po, CHARGE *charge)
{
  static double *Itime = NULL;           /* array for histogram of particle density */
  static double *Ifreq = NULL;           /* array for FFT of histogram of particle density */
  static double *Vtime = NULL;           /* array for voltage acting on each bin */
  static long max_n_bins = 0;
  static long *pbin = NULL;              /* array to record which bin each particle is in */
  static double *time = NULL;            /* array to record arrival time of each particle */
  static long max_np = 0;
  double *Vfreq, ZImag;
  long ib, nb, n_binned, nfreq, iReal, iImag;
  double factor, tmin, tmax, tmean, dt, df, dk, a1, a2;
  double lengthLeft, Imin, Imax, kSC, Zmax;
  double Ia = 17045, Z0, length, k;
  double S11, S33, beamRadius;
#if DEBUG
  FILE *fpd = NULL;
#endif

  Z0 = sqrt(mu_o/epsilon_o);
  nb = LSC->bins;

#if DEBUG
  fprintf(stdout, "%ld bins for LSC\n", nb);
  fflush(stdout);
#endif

  if (nb>max_n_bins) {
    max_n_bins = nb;
    Itime = trealloc(Itime, 2*sizeof(*Itime)*(max_n_bins+1));
    Ifreq = trealloc(Ifreq, 2*sizeof(*Ifreq)*(max_n_bins+1));
    Vtime = trealloc(Vtime, 2*sizeof(*Vtime)*(max_n_bins+1));
  }

  if (np>max_np) {
    pbin = trealloc(pbin, sizeof(*pbin)*(max_np=np));
    time = trealloc(time, sizeof(*time)*max_np);
  }

  lengthLeft = LSC->length;
  while (lengthLeft>0) {
    /* compute time coordinates and make histogram */
    tmean = computeTimeCoordinates(time, Po, part, np);
    find_min_max(&tmin, &tmax, time, np);
    dt = (tmax-tmin)/(nb-3);
#if DEBUG
    fprintf(stdout, "tmean=%e, tmin=%e, tmax=%e, dt=%e\n",
            tmean, tmin, tmax, dt);
    fflush(stdout);
#endif
    n_binned = binTimeDistribution(Itime, pbin, tmin, dt, nb, time, part, Po, np);
#if DEBUG
    fprintf(stdout, "%ld of %ld particles binned\n", n_binned, np);
    fflush(stdout);
#endif
    if (n_binned!=np) {
      fprintf(stdout, "Warning: only %ld of %ld particles were binned (LSCDRIFT)!\n", n_binned, np);
      fprintf(stdout, "This shouldn't happen.\n");
      fflush(stdout);
    }
    if (LSC->smoothing) {
      SavitzyGolaySmooth(Itime, nb, LSC->SGOrder, LSC->SGHalfWidth, LSC->SGHalfWidth, 0);
#if DEBUG
      fprintf(stdout, "Smoothing completed\n");
      fflush(stdout);
#endif
    }

    /* Compute kSC and length to drift */
    /* - find maximum current */
    find_min_max(&Imin, &Imax, Itime, nb);
#if DEBUG
    fprintf(stdout, "Maximum particles/bin: %e    Q/MP: %e C    Imax: %e A\n", 
            Imax, charge->macroParticleCharge, Imax*charge->macroParticleCharge/dt);
  fflush(stdout);
#endif
    Imax *= charge->macroParticleCharge/dt;
    /* - compute beam radius as the average rms beam size in x and y */
    rms_emittance(part, 0, 2, np, &S11, NULL, &S33);
    if ((beamRadius = (sqrt(S11)+sqrt(S33))/2)==0) {
      fprintf(stdout, "Error: beam radius is zero in LSCDRIFT\n");
      exit(1);
    }
    /* - compute kSC */
    kSC = 2/beamRadius*sqrt(Imax/ipow(Po,3)/Ia);
    /* - compute length to drift */
    length = 0.1/kSC;
    if (length>lengthLeft)
      length = lengthLeft;
    /* - compute values for computing impedance */
    df = 1./(dt*nb);
    dk = df*PIx2/c_mks;

#if DEBUG
    fprintf(stdout, "rb: %e m   LSC: I0=%e A    kSC=%e 1/m\nlength = %e m   dt = %e s    df = %e Hz   dk = %e 1/m\n",
            beamRadius, Imax, kSC, length, dt, df, dk);
  fflush(stdout);
#endif

    /* Take the FFT of I(t) to get I(f) */
    memcpy(Ifreq, Itime, 2*nb*sizeof(*Ifreq));
    realFFT(Ifreq, nb, 0);
    nfreq = nb/2 + 1;

    if (LSC->highFrequencyCutoff0>0) {
      /* apply low-pass filter */
      long i, i1, i2;
      double dfraction, fraction;
      i1 = LSC->highFrequencyCutoff0*nfreq;
      if (i1<1)
	i1 = 1;
      i2 = LSC->highFrequencyCutoff1*nfreq;
      if (i2>=nfreq)
	i2 = nfreq-1;
      dfraction = i1==i2 ? 0 : 1./(i2-i1);
      fraction = 1;
      for (i=i1; i<i2; i++) {
	Ifreq[2*i-1] *= fraction;
	Ifreq[2*i  ] *= fraction;
	if ((fraction -= dfraction)<0)
	  fraction = 0;
      }
      for ( ; i<nfreq-1; i++) {
	Ifreq[2*i-1] = 0;
	Ifreq[2*i  ] = 0;
      }
      /* kill the Nyquist term */
      Ifreq[nb-1] = 0;
    }
    
    /* Compute V(f) = Z(f)*I(f), putting in a factor 
     * to normalize the current waveform.
     */
    Vfreq = Vtime;

    /* impedance is zero at DC */
    Vfreq[0] = 0;
    
    /* Nyquist term for current histogram is pure imaginary.
     * Since impedance is also imaginary, the product is pure real.
     * Hence, the Nyquist term we store is zero.
     */
    if (nb%2==0)
      Vfreq[nb-1] = 0;
    
    factor = charge->macroParticleCharge/dt;
    a2 = Z0/(PI*sqr(beamRadius))*length;
#if DEBUG
    fprintf(stdout, "nfreq = %ld   a2 = %e Ohms/m\n", nfreq, a2);
    fflush(stdout);

    if (!fpd) {
      fpd = fopen_e("lscZ.sdds", "w", 0);
      fprintf(fpd, "SDDS1\n&column name=k type=double units=m &end\n");
      fprintf(fpd, "&column name=ZImag, type=double, units=Ohms &end\n");
      fprintf(fpd, "&data no_row_counts=1 mode=ascii &end\n");
    } else
      fprintf(fpd, "\n");
#endif
    Zmax = 0;
    for (ib=1; ib<nfreq-1; ib++) {
      k = ib*dk;
      a1 = k*beamRadius/Po;        
      ZImag = a2/k*(1-a1*dbesk1(a1));
#if DEBUG
      fprintf(fpd, "%e %e\n", k, ZImag);
#endif
      if (ZImag>Zmax)
        Zmax = ZImag;
      iImag = (iReal = 2*ib-1)+1;
      /* There is a minus sign here because I use t<0 for the head */
      Vfreq[iReal] = Ifreq[iImag]*ZImag*factor;
      Vfreq[iImag] = -Ifreq[iReal]*ZImag*factor;
    }
#if DEBUG
    fprintf(stdout, "Maximum |Z| = %e Ohm\n", Zmax);
    fflush(stdout);
#endif

    /* Compute inverse FFT of V(f) to get V(t) */
    realFFT(Vfreq, nb, INVERSE_FFT);
    Vtime = Vfreq;
    
    /* put zero voltage in Vtime[nb] for use in interpolation */
    Vtime[nb] = 0;
    applyLongitudinalWakeKicks(part, time, pbin, np, Po, Vtime, 
                             nb, tmin, dt, LSC->interpolate);

    /* advance particles to the next step */
    for (ib=0; ib<np; ib++) {
      part[ib][4] += length;
    }
    
    lengthLeft -= length;
  }

#if DEBUG  
  if (fpd)
    fclose(fpd);
#endif

#if defined(MINIMIZE_MEMORY)
  free(Itime);
  free(Vtime);
  free(pbin);
  free(time);
  Itime = Vtime = time = NULL;
  pbin = NULL;
  max_np = max_n_bins = 0;
#endif

}

void addLSCKick(double **part, long np, LSCKICK *LSC, double Po, CHARGE *charge, 
                double lengthScale, double dgammaOverGamma)
{
  static double *Itime = NULL;           /* array for histogram of particle density */
  static double *Ifreq = NULL;           /* array for FFT of histogram of particle density */
  static double *Vtime = NULL;           /* array for voltage acting on each bin */
  static long max_n_bins = 0;
  static long *pbin = NULL;              /* array to record which bin each particle is in */
  static double *time = NULL;            /* array to record arrival time of each particle */
  static long max_np = 0;
  double *Vfreq, ZImag;
  long ib, nb, n_binned, nfreq, iReal, iImag;
  double factor, tmin, tmax, tmean, dt, df, dk, a1, a2;
  double Imin, Imax, kSC, Zmax;
  double Ia = 17045, Z0, length, k;
  double S11, S33, beamRadius;
  
  Z0 = sqrt(mu_o/epsilon_o);
  nb = LSC->bins;

#if DEBUG
  fprintf(stdout, "%ld bins for LSC\n", nb);
  fflush(stdout);
#endif

  if (nb>max_n_bins) {
    max_n_bins = nb;
    Itime = trealloc(Itime, 2*sizeof(*Itime)*(max_n_bins+1));
    Ifreq = trealloc(Ifreq, 2*sizeof(*Ifreq)*(max_n_bins+1));
    Vtime = trealloc(Vtime, 2*sizeof(*Vtime)*(max_n_bins+1));
  }

  if (np>max_np) {
    pbin = trealloc(pbin, sizeof(*pbin)*(max_np=np));
    time = trealloc(time, sizeof(*time)*max_np);
  }

  /* compute time coordinates and make histogram */
  tmean = computeTimeCoordinates(time, Po, part, np);
  find_min_max(&tmin, &tmax, time, np);
  dt = (tmax-tmin)/(nb-3);
#if DEBUG
  fprintf(stdout, "tmean=%e, tmin=%e, tmax=%e, dt=%e\n",
          tmean, tmin, tmax, dt);
  fflush(stdout);
#endif
  n_binned = binTimeDistribution(Itime, pbin, tmin, dt, nb, time, part, Po, np);
#if DEBUG
  fprintf(stdout, "%ld of %ld particles binned\n", n_binned, np);
  fflush(stdout);
#endif
  if (n_binned!=np) {
    fprintf(stdout, "Warning: only %ld of %ld particles were binned (LSCDRIFT)!\n", n_binned, np);
    fprintf(stdout, "This shouldn't happen.\n");
    fflush(stdout);
  }

  /* Compute kSC and length to drift */
  /* - compute values involved in binning and FFTs */
  df = 1./(dt*nb);
  dk = df*PIx2/c_mks;
  /* - find maximum current */
  find_min_max(&Imin, &Imax, Itime, nb);
#if DEBUG
  fprintf(stdout, "Maximum particles/bin: %e    Q/MP: %e C    Imax: %e A\n", 
          Imax, charge->macroParticleCharge, Imax*charge->macroParticleCharge/dt);
  fflush(stdout);
#endif
  Imax *= charge->macroParticleCharge/dt;
  /* - compute beam radius as the average rms beam size in x and y */
  rms_emittance(part, 0, 2, np, &S11, NULL, &S33);
  if ((beamRadius = (sqrt(S11)+sqrt(S33))/2)==0) {
    fprintf(stdout, "Error: beam radius is zero in LSCDRIFT\n");
    exit(1);
  }
  /* - compute kSC */
  kSC = 2/beamRadius*sqrt(Imax/ipow(Po,3)/Ia);

  /* - compute maximum length that we should be traveling between kicks */
#if DEBUG
  fprintf(stdout, "rb=%e m   I0=%e A    kSC=%e 1/m    dt=%e s    df=%e Hz   dk=%e 1/m\n",
          beamRadius, Imax, kSC, dt, df, dk);
  fprintf(stdout, "lengthScale=%e m   dgamma/gamma=%e\n", lengthScale, dgammaOverGamma);
  fflush(stdout);
#endif
  length = 1/kSC;
  if (dgammaOverGamma) {
    double length2;
    length2 = fabs(lengthScale/dgammaOverGamma);
    if (length2<length)
      length = length2;
  }
  length /= 10;
  
  if (length<lengthScale) {
    /* length scale used by calling routine is too large, so refuse to do it */
    TRACKING_CONTEXT context;
    getTrackingContext(&context);
    fprintf(stdout, "Error: distance between LSC kicks for %s at z=%e is too large.\n",
            context.elementName, context.zStart);
    fprintf(stdout, "Suggest reducing distance between kicks by factor %e\n",
            lengthScale/length);
    exit(1);
  }    

  /* Take the FFT of I(t) to get I(f) */
  memcpy(Ifreq, Itime, 2*nb*sizeof(*Ifreq));
  realFFT(Ifreq, nb, 0);
  nfreq = nb/2 + 1;

  if (LSC->highFrequencyCutoff0>0) {
    /* apply low-pass filter */
    long i, i1, i2;
    double dfraction, fraction;
    i1 = LSC->highFrequencyCutoff0*nfreq;
    if (i1<1)
      i1 = 1;
    i2 = LSC->highFrequencyCutoff1*nfreq;
    if (i2>=nfreq)
      i2 = nfreq-1;
    dfraction = i1==i2 ? 0 : 1./(i2-i1);
    fraction = 1;
    for (i=i1; i<i2; i++) {
      Ifreq[2*i-1] *= fraction;
      Ifreq[2*i  ] *= fraction;
      if ((fraction -= dfraction)<0)
        fraction = 0;
    }
    for ( ; i<nfreq-1; i++) {
      Ifreq[2*i-1] = 0;
      Ifreq[2*i  ] = 0;
    }
    /* kill the Nyquist term */
    Ifreq[nb-1] = 0;
  }
  
  /* Compute V(f) = Z(f)*I(f), putting in a factor 
   * to normalize the current waveform.
   */
  Vfreq = Vtime;

  /* impedance is zero at DC */
  Vfreq[0] = 0;
  
  /* Nyquist term for current histogram is pure imaginary.
   * Since impedance is also imaginary, the product is pure real.
   * Hence, the Nyquist term we store is zero.
   */
  if (nb%2==0)
    Vfreq[nb-1] = 0;
  
  factor = charge->macroParticleCharge/dt;
  a2 = Z0/(PI*sqr(beamRadius))*lengthScale;
#if DEBUG
  fprintf(stdout, "nfreq = %ld   a2 = %e Ohms/m\n", nfreq, a2);
  fflush(stdout);
#endif
  Zmax = 0;
  for (ib=1; ib<nfreq-1; ib++) {
    k = ib*dk;
    a1 = k*beamRadius/Po;        
    ZImag = a2/k*(1-a1*dbesk1(a1));
    if (ZImag>Zmax)
      Zmax = ZImag;
    iImag = (iReal = 2*ib-1)+1;
    /* There is a minus sign here because I use t<0 for the head */
    Vfreq[iReal] = Ifreq[iImag]*ZImag*factor;
    Vfreq[iImag] = -Ifreq[iReal]*ZImag*factor;
  }
#if DEBUG
  fprintf(stdout, "Maximum |Z| = %e Ohm\n", Zmax);
  fflush(stdout);
#endif

  /* Compute inverse FFT of V(f) to get V(t) */
  realFFT(Vfreq, nb, INVERSE_FFT);
  Vtime = Vfreq;
    
  /* put zero voltage in Vtime[nb] for use in interpolation */
  Vtime[nb] = 0;
  applyLongitudinalWakeKicks(part, time, pbin, np, Po, Vtime, 
                             nb, tmin, dt, LSC->interpolate);
  
#if defined(MINIMIZE_MEMORY)
  free(Itime);
  free(Vtime);
  free(pbin);
  free(time);
  Itime = Vtime = time = NULL;
  pbin = NULL;
  max_np = max_n_bins = 0;
#endif

}

