/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: frfmode.c
 * contents: track_through_frfmode(), set_up_frfmode()
 *
 * Michael Borland, 1993, 2003
 */
#include "mdb.h"
#include "track.h"

void track_through_frfmode(
                           double **part, long np, FRFMODE *rfmode, double Po,
                           char *element_name, double element_z, long pass, long n_passes,
                           CHARGE *charge
                           )
{
  static long *Ihist = NULL;               /* array for histogram of particle density */
  static double *Vbin = NULL;              /* array for voltage acting on each bin */
  static long max_n_bins = 0;
  static long *pbin = NULL;                /* array to record which bin each particle is in */
  static double *time = NULL;              /* array to record arrival time of each particle */
  static long max_np = 0;
  long ip, ib, nb2, lastBin, n_binned;
  double tmin, tmax, tmean, dt, P;
  double Vb, V, omega, phase, t, k, damping_factor, tau;
  double V_sum, Vr_sum, phase_sum;
  double Vc, Vcr, dgamma;
  long max_hist, n_occupied, imode;
  double Qrp, VbImagFactor, Q;
  
  if (charge)
    rfmode->mp_charge = charge->macroParticleCharge;
  else
    bomb("CHARGE element required to use FRFMODE", NULL);
  
  if (rfmode->mp_charge==0 || rfmode->factor==0)
    return ;
  
  if (!rfmode->initialized)
    bomb("track_through_rfmode called with uninitialized element", NULL);

  if (rfmode->n_bins>max_n_bins) {
    Ihist = trealloc(Ihist, sizeof(*Ihist)*(max_n_bins=rfmode->n_bins));
    Vbin = trealloc(Vbin, sizeof(*Vbin)*max_n_bins);
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
  tmin = tmean - rfmode->bin_size*rfmode->n_bins/2.;
  tmax = tmean + rfmode->bin_size*rfmode->n_bins/2.;

  for (ib=0; ib<rfmode->n_bins; ib++) {
    Ihist[ib] = 0;
    Vbin[ib] = 0;
  }
  
  dt = (tmax - tmin)/rfmode->n_bins;
  n_binned = lastBin = 0;
  for (ip=0; ip<np; ip++) {
    pbin[ip] = -1;
    ib = (time[ip]-tmin)/dt;
    if (ib<0)
      continue;
    if (ib>rfmode->n_bins - 1)
      continue;
    Ihist[ib] += 1;
    pbin[ip] = ib;
    if (ib>lastBin)
      lastBin = ib;
    n_binned++;
  }

  for (ib=0; ib<=lastBin; ib++) {
    if (!Ihist[ib])
      continue;
    
    for (imode=0; imode<rfmode->modes; imode++) {
      if (rfmode->cutoffFrequency>0 && (rfmode->omega[imode] > PIx2*rfmode->cutoffFrequency))
        continue;
      
      V_sum = Vr_sum = phase_sum = Vc = Vcr = 0;
      max_hist = n_occupied = 0;
      nb2 = rfmode->n_bins/2;
      
      omega = rfmode->omega[imode];
      Q = rfmode->Q[imode]/(1+rfmode->beta[imode]);
      tau = 2*Q/omega;
      k = omega/2*rfmode->Rs[imode]*rfmode->factor/rfmode->Q[imode];

      /* These adjustments per Zotter and Kheifets, 3.2.4 */
      Qrp = sqrt(Q*Q - 0.25);
      VbImagFactor = 1/(2*Qrp);
      omega *= Qrp/Q;

      t = tmin+(ib+0.5)*dt;           /* middle arrival time for this bin */
      
      /* advance cavity to this time */
      phase = rfmode->last_phase[imode] + omega*(t - rfmode->last_t);
      damping_factor = exp(-(t-rfmode->last_t)/tau);
      rfmode->last_phase[imode] = phase;
      V = rfmode->V[imode]*damping_factor;
      rfmode->Vr[imode] = V*cos(phase);
      rfmode->Vi[imode] = V*sin(phase);

      /* compute beam-induced voltage for this bin */
      Vb = 2*k*rfmode->mp_charge*Ihist[ib]; 
      Vbin[ib] += rfmode->Vr[imode] - Vb/2;
      
      /* add beam-induced voltage to cavity voltage */
      rfmode->Vr[imode] -= Vb;
      rfmode->Vi[imode] -= Vb*VbImagFactor;
      rfmode->last_phase[imode] = atan2(rfmode->Vi[imode], rfmode->Vr[imode]);
      rfmode->V[imode] = sqrt(sqr(rfmode->Vr[imode])+sqr(rfmode->Vi[imode]));
      
      V_sum  += Ihist[ib]*rfmode->V[imode];
      Vr_sum += Ihist[ib]*rfmode->Vr[imode];
      phase_sum += Ihist[ib]*rfmode->last_phase[imode];
      
    }
    rfmode->last_t = t;
  }
  
  if (rfmode->rigid_until_pass<=pass) {
    /* change particle momentum offsets to reflect voltage in relevant bin */
    /* also recompute slopes for new momentum to conserve transverse momentum */
    for (ip=0; ip<np; ip++) {
      if (pbin[ip]>=0) {
        /* compute new momentum and momentum offset for this particle */
        dgamma = Vbin[pbin[ip]]/(1e6*me_mev);
        add_to_particle_energy(part[ip], time[ip], Po, dgamma);
      }
    }
  }
  
#if defined(MINIMIZE_MEMORY)
  free(Ihist);
  free(Vbin);
  free(pbin);
  free(time);
  Ihist = pbin = NULL;
  Vbin = time = NULL;
  max_n_bins = max_np = 0;
#endif

}


#include "complex.h"

void set_up_frfmode(FRFMODE *rfmode, char *element_name, double element_z, long n_passes, 
                    RUN *run, long n_particles,
                    double Po, double total_length)
{
  long i, n, imode;
  double T;
  SDDS_DATASET SDDSin;
  
  if (rfmode->initialized)
    return;
  if (n_particles<1)
    bomb("too few particles in set_up_frfmode()", NULL);
  if (rfmode->n_bins<2)
    bomb("too few bins for FRFMODE", NULL);
  if (!rfmode->filename ||
      !SDDS_InitializeInput(&SDDSin, rfmode->filename))
    bomb("unable to open file for FRFMODE element", NULL);
  /* check existence and properties of required columns  */
  if (SDDS_CheckColumn(&SDDSin, "Frequency", "Hz", SDDS_ANY_FLOATING_TYPE,
                       stdout)!=SDDS_CHECK_OK) {
    fprintf(stdout, "Error: problem with Frequency column for FRFMODE file %s.  Check existence, type, and units.\n", rfmode->filename);
    exit(1);
  }
  if (SDDS_CheckColumn(&SDDSin, "Frequency", "Hz", SDDS_ANY_FLOATING_TYPE,
                       stdout)!=SDDS_CHECK_OK) {
    fprintf(stdout, "Error: problem with Frequency column for FRFMODE file %s.  Check existence, type, and units.\n", rfmode->filename);
    exit(1);
  }
  if (SDDS_CheckColumn(&SDDSin, "Q", NULL, SDDS_ANY_FLOATING_TYPE,
                       stdout)!=SDDS_CHECK_OK) {
    fprintf(stdout, "Error: problem with Q column for FRFMODE file %s.  Check existence, type, and units.\n", rfmode->filename);
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
  if (rfmode->useSymmData) {
    if (SDDS_CheckColumn(&SDDSin, "ShuntImpedanceSymm", "$gW$r", SDDS_ANY_FLOATING_TYPE,
                         NULL)!=SDDS_CHECK_OK &&
        SDDS_CheckColumn(&SDDSin, "ShuntImpedanceSymm", "Ohms", SDDS_ANY_FLOATING_TYPE,
                         NULL)!=SDDS_CHECK_OK) {
      fprintf(stdout, "Error: problem with ShuntImpedanceSymm column for FRFMODE file %s.  Check existence, type, and units.\n", rfmode->filename);
      exit(1);
    }
  }
  else {
    if (SDDS_CheckColumn(&SDDSin, "ShuntImpedance", "$gW$r", SDDS_ANY_FLOATING_TYPE,
                         NULL)!=SDDS_CHECK_OK &&
        SDDS_CheckColumn(&SDDSin, "ShuntImpedance", "Ohms", SDDS_ANY_FLOATING_TYPE,
                         NULL)!=SDDS_CHECK_OK) {
      fprintf(stdout, "Error: problem with ShuntImpedance column for FRFMODE file %s.  Check existence, type, and units.\n", rfmode->filename);
      exit(1);
    }
  }

  if (!SDDS_ReadPage(&SDDSin))
    SDDS_Bomb("unable to read page from file for FRFMODE element");
  if ((rfmode->modes = SDDS_RowCount(&SDDSin))<1) {
    fprintf(stdout, "Error: no data in FRFMODE file %s\n", rfmode->filename);
    exit(1);
  }
  if (!(rfmode->omega = SDDS_GetColumnInDoubles(&SDDSin, "Frequency")) ||
      !(rfmode->Q = SDDS_GetColumnInDoubles(&SDDSin, "Q")) ||
      (rfmode->useSymmData &&
       !(rfmode->Rs = SDDS_GetColumnInDoubles(&SDDSin, "ShuntImpedanceSymm"))) ||
      (!rfmode->useSymmData &&
       !(rfmode->Rs = SDDS_GetColumnInDoubles(&SDDSin, "ShuntImpedance")))) 
    SDDS_Bomb("Problem getting data from FRFMODE file");
  if (!(rfmode->beta = SDDS_GetColumnInDoubles(&SDDSin, "beta"))) {
    if (!(rfmode->beta = malloc(sizeof(*(rfmode->beta))*rfmode->modes)))
      bomb("memory allocation failure (FRFMODE)", NULL);
    for (imode=0; imode<rfmode->modes; imode++)
      rfmode->beta[imode] = 0;
  }

  if (!(rfmode->V  = malloc(sizeof(*(rfmode->V ))*rfmode->modes)) ||
      !(rfmode->Vr = malloc(sizeof(*(rfmode->Vr))*rfmode->modes)) ||
      !(rfmode->Vi = malloc(sizeof(*(rfmode->Vi))*rfmode->modes)) ||
      !(rfmode->last_phase = malloc(sizeof(*(rfmode->last_phase))*rfmode->modes)))
    bomb("memory allocation failure (FRFMODE)", NULL);
  
  for (imode=0; imode<rfmode->modes; imode++) {
    rfmode->omega[imode] *= PIx2;
    rfmode->V[imode] = rfmode->Vr[imode] = rfmode->Vi[imode] = 0;
    rfmode->last_phase[imode] = 0;
  }


  for (imode=0; imode<rfmode->modes; imode++) {
    if (rfmode->bin_size*rfmode->omega[imode]/PIx2>0.1) {
      T = rfmode->bin_size*rfmode->n_bins;
      rfmode->bin_size = 0.1/(rfmode->omega[imode]/PIx2);
      rfmode->n_bins = T/rfmode->bin_size+1;
      rfmode->bin_size = T/rfmode->n_bins;
      fprintf(stdout, "The FRFMODE %s bin size is too large for mode %ld--setting to %e and increasing to %ld bins\n",
              element_name, imode, rfmode->bin_size, rfmode->n_bins);
      fflush(stdout);
    }
  }

  for (imode=0; imode<rfmode->modes; imode++) {
    if (rfmode->bin_size*rfmode->omega[imode]/PIx2>0.1) {
      fprintf(stdout, "Error: FRFMODE bin size adjustment failed\n");
      exit(1);
    }
  }
  

  rfmode->last_t = element_z/c_mks;

  rfmode->initialized = 1;
}

