/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: ftrfmode.c
 * contents: track_through_ftrfmode(), set_up_ftrfmode(),
 *
 * Michael Borland, 1993, 2003
 */
#include "mdb.h"
#include "track.h"

void track_through_ftrfmode(
                            double **part, long np, FTRFMODE *trfmode, double Po,
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
  long lastBin, imode;

  if (charge)
    trfmode->mp_charge = charge->macroParticleCharge;
  else
    bomb("CHARGE element required to use FTRFMODE", NULL);
  if (trfmode->mp_charge==0 || (trfmode->xfactor==0 && trfmode->yfactor==0))
    return;
  
  if (!trfmode->initialized)
    bomb("track_through_ftrfmode called with uninitialized element", NULL);

  
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
  }

  for (ib=0; ib<=lastBin; ib++) {
    if (!xsum[ib] && !ysum[ib])
      continue;
    t = tmin+(ib+0.5)*dt;           /* middle arrival time for this bin */
    for (imode=0; imode<trfmode->modes; imode++) {
      if (trfmode->cutoffFrequency>0 && (trfmode->omega[imode] > PIx2*trfmode->cutoffFrequency))
        continue;

      omega = trfmode->omega[imode];
      Q = trfmode->Q[imode]/(1+trfmode->beta[imode]);
      tau = 2*Q/omega;
      Qrp = sqrt(Q*Q - 0.25);
      k = omega/2*trfmode->Rs[imode]/trfmode->Q[imode];
      /* These adjustments per Zotter and Kheifets, 3.2.4, 3.3.2 */
      k *= Q/Qrp;
      omega *= Qrp/Q;
      
      if (!trfmode->doX[imode] && !trfmode->doY[imode])
        bomb("x and y turned off for FTRFMODE---this shouldn't happen", NULL);
      
      /* advance cavity to this time */
      damping_factor = exp(-(t-trfmode->last_t)/tau);
      if (trfmode->doX[imode]) {
        /* -- x plane */
        /* advance the phasor */
        phase = trfmode->lastPhasex[imode] + omega*(t - trfmode->last_t);
        V = trfmode->Vx[imode]*damping_factor;
        trfmode->Vxr[imode] = V*cos(phase);
        trfmode->Vxi[imode] = V*sin(phase);
        trfmode->lastPhasex[imode] = phase;
        /* add this cavity's contribution to this bin */
        Vxbin[ib] += trfmode->Vxr[imode];
        /* compute beam-induced voltage for this bin */
        Vxb = 2*k*trfmode->mp_charge*xsum[ib]*trfmode->xfactor; 
        /* add beam-induced voltage to cavity voltage---it is imaginary as
         * the voltage is 90deg out of phase 
         */
        trfmode->Vxi[imode] -= Vxb;
        if (trfmode->Vxi[imode]==0 && trfmode->Vxr[imode]==0)
          trfmode->lastPhasex[imode] = 0;
        else
          trfmode->lastPhasex[imode] = atan2(trfmode->Vxi[imode], trfmode->Vxr[imode]);
        trfmode->Vx[imode] = sqrt(sqr(trfmode->Vxr[imode])+sqr(trfmode->Vxi[imode]));
      }
      if (trfmode->doY[imode]) {
        /* -- y plane */
        /* advance the phasor */
        phase = trfmode->lastPhasey[imode] + omega*(t - trfmode->last_t);
        V = trfmode->Vy[imode]*damping_factor;
        trfmode->Vyr[imode] = V*cos(phase);
        trfmode->Vyi[imode] = V*sin(phase);
        trfmode->lastPhasey[imode] = phase;
        /* add this cavity's contribution to this bin */
        Vybin[ib] += trfmode->Vyr[imode];
        /* compute beam-induced voltage for this bin */
        Vyb = 2*k*trfmode->mp_charge*ysum[ib]*trfmode->yfactor;
        /* add beam-induced voltage to cavity voltage---it is imaginary as
         * the voltage is 90deg out of phase 
         */
        trfmode->Vyi[imode] -= Vyb;
        if (trfmode->Vyi[imode]==0 && trfmode->Vyr[imode]==0)
          trfmode->lastPhasey[imode] = 0;
        else
          trfmode->lastPhasey[imode] = atan2(trfmode->Vyi[imode], trfmode->Vyr[imode]);
        trfmode->Vy[imode] = sqrt(sqr(trfmode->Vyr[imode])+sqr(trfmode->Vyi[imode]));
      }
    }
    trfmode->last_t = t;
  }
  
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

void set_up_ftrfmode(FTRFMODE *rfmode, char *element_name, double element_z, long n_passes, 
                     RUN *run, long n_particles,
                     double Po, double total_length)
{
  long i, n, imode;
  double T;
  SDDS_DATASET SDDSin;
  
  if (rfmode->initialized)
    return;
  if (n_particles<1)
    bomb("too few particles in set_up_ftrfmode()", NULL);
  if (rfmode->n_bins<2)
    bomb("too few bins for FTRFMODE", NULL);
  if (!rfmode->filename ||
      !SDDS_InitializeInput(&SDDSin, rfmode->filename))
    bomb("unable to open file for FTRFMODE element", NULL);
  /* check existence and properties of required columns  */
  if (SDDS_CheckColumn(&SDDSin, "Frequency", "Hz", SDDS_ANY_FLOATING_TYPE,
                       stdout)!=SDDS_CHECK_OK) {
    fprintf(stdout, "Error: problem with Frequency column for FTRFMODE file %s.  Check existence, type, and units.\n", rfmode->filename);
    exit(1);
  }
  if (SDDS_CheckColumn(&SDDSin, "Frequency", "Hz", SDDS_ANY_FLOATING_TYPE,
                       stdout)!=SDDS_CHECK_OK) {
    fprintf(stdout, "Error: problem with Frequency column for FTRFMODE file %s.  Check existence, type, and units.\n", rfmode->filename);
    exit(1);
  }
  if (SDDS_CheckColumn(&SDDSin, "Q", NULL, SDDS_ANY_FLOATING_TYPE,
                       stdout)!=SDDS_CHECK_OK) {
    fprintf(stdout, "Error: problem with Q column for FTRFMODE file %s.  Check existence, type, and units.\n", rfmode->filename);
    exit(1);
  }
  if (rfmode->useSymmData) {
    if (SDDS_CheckColumn(&SDDSin, "ShuntImpedanceSymm", "$gW$r/m", SDDS_ANY_FLOATING_TYPE,
                         NULL)!=SDDS_CHECK_OK &&
        SDDS_CheckColumn(&SDDSin, "ShuntImpedanceSymm", "Ohms/m", SDDS_ANY_FLOATING_TYPE,
                         NULL)!=SDDS_CHECK_OK) {
      fprintf(stdout, "Error: problem with ShuntImpedanceSymm column for FTRFMODE file %s.  Check existence, type, and units.\n", rfmode->filename);
      exit(1);
    }
  }
  else {
    if (SDDS_CheckColumn(&SDDSin, "ShuntImpedance", "$gW$r/m", SDDS_ANY_FLOATING_TYPE,
                         NULL)!=SDDS_CHECK_OK &&
        SDDS_CheckColumn(&SDDSin, "ShuntImpedance", "Ohms/m", SDDS_ANY_FLOATING_TYPE,
                         NULL)!=SDDS_CHECK_OK) {
      fprintf(stdout, "Error: problem with ShuntImpedance column for FTRFMODE file %s.  Check existence, type, and units.\n", rfmode->filename);
      exit(1);
    }
  }

  if (!SDDS_ReadPage(&SDDSin))
    SDDS_Bomb("unable to read page from file for FTRFMODE element");
  if ((rfmode->modes = SDDS_RowCount(&SDDSin))<1) {
    fprintf(stdout, "Error: no data in FTRFMODE file %s\n", rfmode->filename);
    exit(1);
  }
  if (SDDS_CheckColumn(&SDDSin, "beta", NULL, SDDS_ANY_FLOATING_TYPE,
                       NULL)!=SDDS_CHECK_NONEXISTENT) {
    if (SDDS_CheckColumn(&SDDSin, "beta", NULL, SDDS_ANY_FLOATING_TYPE,
                         NULL)!=SDDS_CHECK_OK) {
      fprintf(stdout, "Error: problem with \"beta\" column for FRFMODE file %s.  Check type and units.\n", rfmode->filename);
      exit(1);
    }
  }
  if (SDDS_CheckColumn(&SDDSin, "doX", NULL, SDDS_ANY_INTEGER_TYPE,
                       NULL)!=SDDS_CHECK_NONEXISTENT) {
    if (SDDS_CheckColumn(&SDDSin, "doX", NULL, SDDS_ANY_INTEGER_TYPE,
                         NULL)!=SDDS_CHECK_OK) {
      fprintf(stdout, "Error: problem with \"doX\" column for FTRFMODE file %s.  Check type and units.\n", rfmode->filename);
      exit(1);
    }
  }
  if (SDDS_CheckColumn(&SDDSin, "doY", NULL, SDDS_ANY_INTEGER_TYPE,
                       NULL)!=SDDS_CHECK_NONEXISTENT) {
    if (SDDS_CheckColumn(&SDDSin, "doY", NULL, SDDS_ANY_INTEGER_TYPE,
                         NULL)!=SDDS_CHECK_OK) {
      fprintf(stdout, "Error: problem with \"doY\" column for FTRFMODE file %s.  Check type and units.\n", rfmode->filename);
      exit(1);
    }
  }
  if (!(rfmode->omega = SDDS_GetColumnInDoubles(&SDDSin, "Frequency")) ||
      !(rfmode->Q = SDDS_GetColumnInDoubles(&SDDSin, "Q")) ||
      (rfmode->useSymmData &&
       !(rfmode->Rs = SDDS_GetColumnInDoubles(&SDDSin, "ShuntImpedanceSymm"))) ||
      (!rfmode->useSymmData &&
       !(rfmode->Rs = SDDS_GetColumnInDoubles(&SDDSin, "ShuntImpedance")))) 
    SDDS_Bomb("Problem getting data from FTRFMODE file");
  if (!(rfmode->beta = SDDS_GetColumnInDoubles(&SDDSin, "beta"))) {
    if (!(rfmode->beta = malloc(sizeof(*(rfmode->beta))*rfmode->modes)))
      bomb("memory allocation failure (FTRFMODE)", NULL);
    for (imode=0; imode<rfmode->modes; imode++)
      rfmode->beta[imode] = 0;
  }
  if (!(rfmode->doX = SDDS_GetColumnInLong(&SDDSin, "xMode"))) {
    if (!(rfmode->doX = malloc(sizeof(*(rfmode->doX))*rfmode->modes)))
      bomb("memory allocation failure (FTRFMODE)", NULL);
    for (imode=0; imode<rfmode->modes; imode++)
      rfmode->doX[imode] = 1;
  }
  if (!(rfmode->doY = SDDS_GetColumnInLong(&SDDSin, "yMode"))) {
    if (!(rfmode->doY = malloc(sizeof(*(rfmode->doY))*rfmode->modes)))
      bomb("memory allocation failure (FTRFMODE)", NULL);
    for (imode=0; imode<rfmode->modes; imode++)
      rfmode->doY[imode] = 1;
  }

  if (!(rfmode->Vx  = malloc(sizeof(*(rfmode->Vx ))*rfmode->modes)) ||
      !(rfmode->Vxr = malloc(sizeof(*(rfmode->Vxr))*rfmode->modes)) ||
      !(rfmode->Vxi = malloc(sizeof(*(rfmode->Vxi))*rfmode->modes)) ||
      !(rfmode->Vy  = malloc(sizeof(*(rfmode->Vy ))*rfmode->modes)) ||
      !(rfmode->Vyr = malloc(sizeof(*(rfmode->Vyr))*rfmode->modes)) ||
      !(rfmode->Vyi = malloc(sizeof(*(rfmode->Vyi))*rfmode->modes)) ||
      !(rfmode->lastPhasex = malloc(sizeof(*(rfmode->lastPhasex))*rfmode->modes)) ||
      !(rfmode->lastPhasey = malloc(sizeof(*(rfmode->lastPhasey))*rfmode->modes)))
    bomb("memory allocation failure (FTRFMODE)", NULL);
  
  for (imode=0; imode<rfmode->modes; imode++) {
    rfmode->omega[imode] *= PIx2;
    rfmode->Vx[imode] = rfmode->Vxr[imode] = rfmode->Vxi[imode] = 0;
    rfmode->lastPhasex[imode] = 0;
    rfmode->Vy[imode] = rfmode->Vyr[imode] = rfmode->Vyi[imode] = 0;
    rfmode->lastPhasey[imode] = 0;
  }

  for (imode=0; imode<rfmode->modes; imode++) {
    if (rfmode->bin_size*rfmode->omega[imode]/PIx2>0.1) {
      T = rfmode->bin_size*rfmode->n_bins;
      rfmode->bin_size = 0.1/(rfmode->omega[imode]/PIx2);
      rfmode->n_bins = T/rfmode->bin_size+1;
      rfmode->bin_size = T/rfmode->n_bins;
      fprintf(stdout, "The FTRFMODE %s bin size is too large for mode %ld--setting to %e and increasing to %ld bins\n",
              element_name, imode, rfmode->bin_size, rfmode->n_bins);
      fflush(stdout);
    }
  }

  for (imode=0; imode<rfmode->modes; imode++) {
    if (rfmode->bin_size*rfmode->omega[imode]/PIx2>0.1) {
      fprintf(stdout, "Error: FTRFMODE bin size adjustment failed\n");
      exit(1);
    }
  }

  rfmode->last_t = element_z/c_mks;

  rfmode->initialized = 1;
}

