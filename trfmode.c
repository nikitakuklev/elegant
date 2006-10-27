/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: trfmode.c
 * contents: track_through_trfmode()
 *
 * Michael Borland, 1999
 */
#include "mdb.h"
#include "track.h"

void runBinlessTrfMode(double **part, long np, TRFMODE *trfmode, double Po,
                       char *element_name, double element_z, long pass, long n_passes,
                       CHARGE *charge);

#define DEBUG 0

void track_through_trfmode(
                           double **part, long np, TRFMODE *trfmode, double Po,
                           char *element_name, double element_z, long pass, long n_passes,
                           CHARGE *charge
                           )
{
  static double *xsum = NULL;              /* sum of x coordinate in each bin = N*<x> */
  static double *ysum = NULL;              /* sum of y coordinate in each bin = N*<y> */
  static double *Vxbin = NULL;             /* array for voltage acting on each bin MV */
  static double *Vybin = NULL;             /* array for voltage acting on each bin MV */
  static long max_n_bins = 0;
  static long *pbin = NULL;                /* array to record which bin each particle is in */
  static double *time = NULL;              /* array to record arrival time of each particle */
  static long max_np = 0;
  long ip, ib;
  double tmin, tmax, tmean, dt, P;
  double Vxb, Vyb, V, omega, phase, t, k, damping_factor, tau;
  double Px, Py, Pz;
  double Q, Qrp;
  long n_binned, lastBin;
  static long been_warned = 0;
#if DEBUG
  static FILE *fpdeb = NULL;
  static long debugPass = 0;
#endif

  if (trfmode->binless) {
    runBinlessTrfMode(part, np, trfmode, Po, element_name, element_z, pass, n_passes, charge);
    return;
  }
  
  if (charge) {
    trfmode->mp_charge = charge->macroParticleCharge;
  } else if (pass==0) {
    trfmode->mp_charge = 0;
    if (np)
      trfmode->mp_charge = trfmode->charge/np;
  }

#if DEBUG
  if (!fpdeb) {
    fpdeb = fopen("trfmode.debug", "w");
    fprintf(fpdeb, "SDDS1\n&parameter name=Pass type=long &end\n");
    fprintf(fpdeb, "&parameter name=nBinned type=long &end\n");
    fprintf(fpdeb, "&column name=Bin , type=double &end\n");
    fprintf(fpdeb, "&column name=xSum , type=double &end\n");
    fprintf(fpdeb, "&column name=ySum , type=double &end\n");
    fprintf(fpdeb, "&column name=xVoltage , type=double &end\n");
    fprintf(fpdeb, "&column name=yVoltage , type=double &end\n");
    fprintf(fpdeb, "&data mode=ascii &end\n");
  }
#endif

  omega = PIx2*trfmode->freq;
  if ((Q = trfmode->Q/(1+trfmode->beta))<=0.5) {
    fprintf(stdout, "The effective Q<=0.5 for TRFMODE.  Use the ZTRANSVERSE element.\n");
    fflush(stdout);
    exit(1);
  }
  tau = 2*Q/omega;
  Qrp = sqrt(Q*Q - 0.25);
  k = omega/4*trfmode->RaInternal/trfmode->Q;

  if (!trfmode->doX && !trfmode->doY)
    bomb("x and y turned off for TRFMODE---this shouldn't happen", NULL);
  
  if (!been_warned) {        
    if (trfmode->freq<1e3 && trfmode->freq)  {
      fprintf(stdout, "\7\7\7warning: your TRFMODE frequency is less than 1kHz--this may be an error\n");
      fflush(stdout);
      been_warned = 1;
    }
    if (been_warned) {
      fprintf(stdout, "units of parameters for TRFMODE are as follows:\n");
      fflush(stdout);
      print_dictionary_entry(stdout, T_TRFMODE, 0, 0);
    }
  }

  if (!trfmode->initialized)
    bomb("track_through_trfmode called with uninitialized element", NULL);

  if (trfmode->n_bins>max_n_bins) {
    max_n_bins = trfmode->n_bins;
    xsum = trealloc(xsum, sizeof(*xsum)*max_n_bins);
    ysum = trealloc(ysum, sizeof(*ysum)*max_n_bins);
    Vxbin = trealloc(Vxbin, sizeof(*Vxbin)*max_n_bins);
    Vybin = trealloc(Vybin, sizeof(*Vybin)*max_n_bins);
  }

  if (np>max_np) {
    pbin = trealloc(pbin, sizeof(*pbin)*(max_np=np));
    time = trealloc(time, sizeof(*time)*max_np);
  }

  tmean = 0;
  for (ip=0; ip<np; ip++) {
    P = Po*(part[ip][5]+1);
    time[ip] = part[ip][4]*sqrt(sqr(P)+1)/(c_mks*P);
    tmean += time[ip];
  }
  tmean /= np;
  tmin = tmean - trfmode->bin_size*trfmode->n_bins/2.;
  tmax = tmean + trfmode->bin_size*trfmode->n_bins/2.;

  for (ib=0; ib<trfmode->n_bins; ib++)
    xsum[ib] = ysum[ib] = 0;
  dt = (tmax - tmin)/trfmode->n_bins;
  n_binned = 0;
  lastBin = -1;
  for (ip=0; ip<np; ip++) {
    pbin[ip] = -1;
    ib = (time[ip]-tmin)/dt;
    if (ib<0)
      continue;
    if (ib>trfmode->n_bins - 1)
      continue;
    xsum[ib] += part[ip][0]-trfmode->dx;
    ysum[ib] += part[ip][2]-trfmode->dy;
    pbin[ip] = ib;
    if (ib>lastBin)
      lastBin = ib;
    n_binned++;
  }
  if (n_binned!=np) {
    fprintf(stdout, "Warning: only %ld of %ld particles binned (TRFMODE)\n",
            n_binned, np);
    fflush(stdout);
  }
  
  /* These adjustments per Zotter and Kheifets, 3.2.4, 3.3.2 */
  k *= Q/Qrp;
  omega *= Qrp/Q;

  if (pass <= (trfmode->rampPasses-1)) 
    k *= (pass+1.0)/trfmode->rampPasses;
    
  if (trfmode->single_pass) {
    trfmode->Vx = trfmode->Vy = 0;
    trfmode->last_t = tmin + ib*dt;
    trfmode->last_xphase = trfmode->last_yphase = 0;
  }

  for (ib=0; ib<=lastBin; ib++) {
    if (!xsum[ib] && !ysum[ib])
      continue;

    t = tmin+(ib+0.5)*dt;           /* middle arrival time for this bin */
    
    /* advance cavity to this time */
    damping_factor = exp(-(t-trfmode->last_t)/tau);
    if (trfmode->doX) {
      /* -- x plane */
      phase = trfmode->last_xphase + omega*(t - trfmode->last_t);
      V = trfmode->Vx*damping_factor;
      trfmode->Vxr = V*cos(phase);
      trfmode->Vxi = V*sin(phase);
      trfmode->last_xphase = phase;
    }
    if (trfmode->doY) {
      /* -- y plane */
      phase = trfmode->last_yphase + omega*(t - trfmode->last_t);
      V = trfmode->Vy*damping_factor;
      trfmode->Vyr = V*cos(phase);
      trfmode->Vyi = V*sin(phase);
      trfmode->last_yphase = phase;
    }
    
    trfmode->last_t = t;

    /* compute beam-induced voltage for this bin */
    if (trfmode->doX) {
      /* -- x plane */
      Vxb = 2*k*trfmode->mp_charge*xsum[ib]*trfmode->xfactor;
      Vxbin[ib] = trfmode->Vxr;
      /* add beam-induced voltage to cavity voltage---it is imaginary as
       * the voltage is 90deg out of phase 
       */
      trfmode->Vxi -= Vxb;
      if (trfmode->Vxi==0 && trfmode->Vxr==0)
        trfmode->last_xphase = 0;
      else
        trfmode->last_xphase = atan2(trfmode->Vxi, trfmode->Vxr);
      trfmode->Vx = sqrt(sqr(trfmode->Vxr)+sqr(trfmode->Vxi));
    }
    if (trfmode->doY) {
      /* -- y plane */
      Vyb = 2*k*trfmode->mp_charge*ysum[ib]*trfmode->yfactor;
      Vybin[ib] = trfmode->Vyr;
      /* add beam-induced voltage to cavity voltage---it is imaginary as
       * the voltage is 90deg out of phase 
       */
      trfmode->Vyi -= Vyb;
      if (trfmode->Vyi==0 && trfmode->Vyr==0)
        trfmode->last_yphase = 0;
      else
        trfmode->last_yphase = atan2(trfmode->Vyi, trfmode->Vyr);
      trfmode->Vy = sqrt(sqr(trfmode->Vyr)+sqr(trfmode->Vyi));
    }    
  }
  
#if DEBUG
  fprintf(fpdeb, "%ld\n%ld\n%ld\n", 
          debugPass++, n_binned, lastBin+1);
  for (ib=0; ib<=lastBin; ib++) {
    fprintf(fpdeb, "%ld %e %e %e %e\n",
            ib, xsum[ib], ysum[ib], 
            xsum[ib]?Vxbin[ib]:0.0,
            ysum[ib]?Vybin[ib]:0.0);
  }
  fflush(fpdeb);
#endif

  /* change particle slopes to reflect voltage in relevant bin */
  for (ip=0; ip<np; ip++) {
    if (pbin[ip]>=0) {
      P = Po*(1+part[ip][5]);
      Pz = P/sqrt(1+sqr(part[ip][1])+sqr(part[ip][3]));
      Px = part[ip][1]*Pz + Vxbin[pbin[ip]]/(1e6*me_mev);
      Py = part[ip][3]*Pz + Vybin[pbin[ip]]/(1e6*me_mev);
      P  = sqrt(Pz*Pz+Px*Px+Py*Py);
      part[ip][1] = Px/Pz;
      part[ip][3] = Py/Pz;
      part[ip][5] = (P-Po)/Po;
      part[ip][4] = time[ip]*c_mks*P/sqrt(sqr(P)+1);
    }
  }

#if defined(MINIMIZE_MEMORY)
  free(xsum);
  free(ysum);
  free(Vxbin);
  free(Vybin);
  free(pbin);
  free(time);
  xsum = ysum = Vxbin = Vybin = time = NULL;
  pbin = NULL;
  max_n_bins =  max_np = 0;
#endif

}


#include "complex.h"

void set_up_trfmode(TRFMODE *trfmode, char *element_name, double element_z, 
                    long n_passes, RUN *run, long n_particles)
{
  double T;

  if (trfmode->initialized)
    return;
  
  trfmode->initialized = 1;
  
  if (n_particles<1)
    bomb("too few particles in set_up_trfmode()", NULL);
  if (trfmode->n_bins<2)
    bomb("too few bins for TRFMODE", NULL);
  if (trfmode->bin_size<=0)
    bomb("bin_size must be positive for TRFMODE", NULL);
  if (trfmode->Ra && trfmode->Rs) 
    bomb("TRFMODE element may have only one of Ra or Rs nonzero.  Ra is just 2*Rs", NULL);
  if (trfmode->Ra)
    trfmode->RaInternal = trfmode->Ra;
  else
    trfmode->RaInternal = 2*trfmode->Rs;
  if (trfmode->bin_size*trfmode->freq>0.1) {
    T = trfmode->bin_size*trfmode->n_bins;
    trfmode->bin_size = 0.1/trfmode->freq;
    trfmode->n_bins = T/trfmode->bin_size;
    fprintf(stdout, "The TRFMODE %s bin size is too large--setting to %e and increasing to %ld bins\n",
            element_name, trfmode->bin_size, trfmode->n_bins);
    fflush(stdout);
  }
  trfmode->last_t = element_z/c_mks;
  trfmode->Vxr = trfmode->Vxi = trfmode->Vx = 0;
  trfmode->Vyr = trfmode->Vyi = trfmode->Vy = 0;
  trfmode->doX = trfmode->doY = 0;
  if (strcmp_ci(trfmode->plane, "BOTH")==0)
    trfmode->doX = trfmode->doY = 1;
  else if (strcmp_ci(trfmode->plane, "X")==0)
    trfmode->doX = 1;
  else if (strcmp_ci(trfmode->plane, "Y")==0)
    trfmode->doY = 1;
  if (!trfmode->doX && !trfmode->doY) 
    bomb("No planes selected for TRFMODE", NULL);
}

void runBinlessTrfMode(
                       double **part, long np, TRFMODE *trfmode, double Po,
                       char *element_name, double element_z, long pass, long n_passes,
                       CHARGE *charge
                       )
{
  static TIMEDATA *tData;
  static long max_np = 0;
  long ip, ip0;
  double P;
  double Vxb, Vyb, dVxb, dVyb, V, omega, phase, t, k, damping_factor, tau;
  double Px, Py, Pz;
  double Q, Qrp;
  double x, y;
  static long been_warned = 0;

  if (np==0)
    return;
  
  if (charge) {
    trfmode->mp_charge = charge->macroParticleCharge;
  } else if (pass==0) {
    trfmode->mp_charge = 0;
    if (np)
      trfmode->mp_charge = trfmode->charge/np;
  }

  omega = PIx2*trfmode->freq;
  if ((Q = trfmode->Q/(1+trfmode->beta))<=0.5) {
    fprintf(stdout, "The effective Q<=0.5 for TRFMODE.  Use the ZTRANSVERSE element.\n");
    fflush(stdout);
    exit(1);
  }
  tau = 2*Q/omega;
  Qrp = sqrt(Q*Q - 0.25);
  k = omega/4*trfmode->RaInternal/trfmode->Q;

  if (!trfmode->doX && !trfmode->doY)
    bomb("x and y turned off for TRFMODE---this shouldn't happen", NULL);
  
  if (!been_warned) {        
    if (trfmode->freq<1e3 && trfmode->freq)  {
      fprintf(stdout, "\7\7\7warning: your TRFMODE frequency is less than 1kHz--this may be an error\n");
      fflush(stdout);
      been_warned = 1;
    }
    if (been_warned) {
      fprintf(stdout, "units of parameters for TRFMODE are as follows:\n");
      fflush(stdout);
      print_dictionary_entry(stdout, T_TRFMODE, 0, 0);
    }
  }

  if (!trfmode->initialized)
    bomb("track_through_trfmode called with uninitialized element", NULL);

  if (np>max_np) 
    tData = trealloc(tData, sizeof(*tData)*(max_np=np));

  for (ip=0; ip<np; ip++) {
    P = Po*(part[ip][5]+1);
    tData[ip].t = part[ip][4]*sqrt(sqr(P)+1)/(c_mks*P);
    tData[ip].ip = ip;
  }
  qsort(tData, np, sizeof(*tData), compTimeData);
  
  /* These adjustments per Zotter and Kheifets, 3.2.4, 3.3.2 */
  k *= Q/Qrp;
  omega *= Qrp/Q;

  if (pass <= (trfmode->rampPasses-1)) 
    k *= (pass+1.0)/trfmode->rampPasses;
    
  if (trfmode->single_pass) {
    trfmode->Vx = trfmode->Vy = 0;
    trfmode->last_t = tData[0].t;
    trfmode->last_xphase = trfmode->last_yphase = 0;
  }

  for (ip0=0; ip0<np; ip0++) {
    ip = tData[ip0].ip;
    x = part[ip][0];
    y = part[ip][2];
    if (x==0 && y==0)
      continue;
    
    t = tData[ip0].t;
    
    /* advance cavity to this time */
    damping_factor = exp(-(t-trfmode->last_t)/tau);
    if (trfmode->doX) {
      /* -- x plane */
      phase = trfmode->last_xphase + omega*(t - trfmode->last_t);
      V = trfmode->Vx*damping_factor;
      trfmode->Vxr = V*cos(phase);
      trfmode->Vxi = V*sin(phase);
      trfmode->last_xphase = phase;
    }
    if (trfmode->doY) {
      /* -- y plane */
      phase = trfmode->last_yphase + omega*(t - trfmode->last_t);
      V = trfmode->Vy*damping_factor;
      trfmode->Vyr = V*cos(phase);
      trfmode->Vyi = V*sin(phase);
      trfmode->last_yphase = phase;
    }
    
    trfmode->last_t = t;

    /* compute beam-induced voltage for this bin */
    if (trfmode->doX) {
      /* -- x plane */
      dVxb = 2*k*trfmode->mp_charge*x*trfmode->xfactor;
      Vxb = trfmode->Vxr;
      /* add beam-induced voltage to cavity voltage---it is imaginary as
       * the voltage is 90deg out of phase 
       */
      trfmode->Vxi -= dVxb;
      if (trfmode->Vxi==0 && trfmode->Vxr==0)
        trfmode->last_xphase = 0;
      else
        trfmode->last_xphase = atan2(trfmode->Vxi, trfmode->Vxr);
      trfmode->Vx = sqrt(sqr(trfmode->Vxr)+sqr(trfmode->Vxi));
    }
    if (trfmode->doY) {
      /* -- y plane */
      dVyb = 2*k*trfmode->mp_charge*y*trfmode->yfactor;
      Vyb = trfmode->Vyr;
      /* add beam-induced voltage to cavity voltage---it is imaginary as
       * the voltage is 90deg out of phase 
       */
      trfmode->Vyi -= dVyb;
      if (trfmode->Vyi==0 && trfmode->Vyr==0)
        trfmode->last_yphase = 0;
      else
        trfmode->last_yphase = atan2(trfmode->Vyi, trfmode->Vyr);
      trfmode->Vy = sqrt(sqr(trfmode->Vyr)+sqr(trfmode->Vyi));
    }    
    
    /* change particle slopes to reflect voltage in relevant bin */
    P = Po*(1+part[ip][5]);
    Pz = P/sqrt(1+sqr(part[ip][1])+sqr(part[ip][3]));
    Px = part[ip][1]*Pz + Vxb/(1e6*me_mev);
    Py = part[ip][3]*Pz + Vyb/(1e6*me_mev);
    P  = sqrt(Pz*Pz+Px*Px+Py*Py);
    part[ip][1] = Px/Pz;
    part[ip][3] = Py/Pz;
    part[ip][5] = (P-Po)/Po;
    part[ip][4] = tData[ip0].t*c_mks*P/sqrt(sqr(P)+1);
  }

#if defined(MINIMIZE_MEMORY)
  free(tData);
  tData = NULL;
  max_np = 0;
#endif

}

int compTimeData(void *tv1, void *tv2)
{
  double diff;
  diff = ((TIMEDATA*)tv1)->t - ((TIMEDATA*)tv2)->t;
  if (diff<0)
    return -1;
  if (diff>0)
    return 1;
  return 0;
}

