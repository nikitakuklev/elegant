/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: rfmode.c
 * contents: track_through_rfmode()
 *
 * Michael Borland, 1993
 */
#include "mdb.h"
#include "track.h"

void track_through_rfmode(
                          double **part, long np, RFMODE *rfmode, double Po,
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
    double Vc, Vcr, Q_sum, dgamma;
    long n_summed, max_hist, n_occupied;
    static long been_warned = 0;
    double Qrp, VbImagFactor, Q;
    
    if (charge) {
      rfmode->mp_charge = charge->macroParticleCharge;
    } else if (pass==0) {
      rfmode->mp_charge = 0;
      if (np)
        rfmode->mp_charge = rfmode->charge/np;
    }

    if (pass%rfmode->pass_interval)
      return;
      
    omega = PIx2*rfmode->freq;
    if ((Q = rfmode->Q/(1+rfmode->beta))<=0.5) {
      fprintf(stdout, "The effective Q<=0.5 for RFMODE.  Use the ZLONGIT element.\n");
      fflush(stdout);
      exit(1);
    }
    tau = 2*Q/omega;
    k = omega/4*(rfmode->Ra)/rfmode->Q;

    /* These adjustments per Zotter and Kheifets, 3.2.4 */
    Qrp = sqrt(Q*Q - 0.25);
    VbImagFactor = 1/(2*Qrp);
    omega *= Qrp/Q;

    if (!been_warned) {        
        if (rfmode->freq<1e3 && rfmode->freq)  {
            fprintf(stdout, "\7\7\7warning: your RFMODE frequency is less than 1kHz--this may be an error\n");
            fflush(stdout);
            been_warned = 1;
            }
        if (been_warned) {
            fprintf(stdout, "units of parameters for RFMODE are as follows:\n");
            fflush(stdout);
            print_dictionary_entry(stdout, T_RFMODE, 0);
            }
        }

    if (rfmode->mp_charge==0) {
      return ;
    }
    if (rfmode->detuned_until_pass>pass) {
      return ;
    }
    
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

    for (ib=0; ib<rfmode->n_bins; ib++)
        Ihist[ib] = 0;

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

    V_sum = Vr_sum = phase_sum = Vc = Vcr = Q_sum = 0;
    n_summed = max_hist = n_occupied = 0;
    nb2 = rfmode->n_bins/2;
    if (rfmode->single_pass)
        rfmode->V = 0;
    for (ib=0; ib<=lastBin; ib++) {
        if (!Ihist[ib])
            continue;
        if (Ihist[ib]>max_hist)
            max_hist = Ihist[ib];
        n_occupied++;

        t = tmin+(ib+0.5)*dt;           /* middle arrival time for this bin */
        
        /* advance cavity to this time */
        phase = rfmode->last_phase + omega*(t - rfmode->last_t);
        damping_factor = exp(-(t-rfmode->last_t)/tau);
        rfmode->last_t = t;
        rfmode->last_phase = phase;
        V = rfmode->V*damping_factor;
        rfmode->Vr = V*cos(phase);
        rfmode->Vi = V*sin(phase);

        /* compute beam-induced voltage for this bin */
        Vb = 2*k*rfmode->mp_charge*rfmode->pass_interval*Ihist[ib]; 
        Vbin[ib] = rfmode->Vr - Vb/2;
        
        /* add beam-induced voltage to cavity voltage */
        rfmode->Vr -= Vb;
        rfmode->Vi -= Vb*VbImagFactor;
        rfmode->last_phase = atan2(rfmode->Vi, rfmode->Vr);
        rfmode->V = sqrt(sqr(rfmode->Vr)+sqr(rfmode->Vi));
 
        V_sum  += Ihist[ib]*rfmode->V;
        Vr_sum += Ihist[ib]*rfmode->Vr;
        phase_sum += Ihist[ib]*rfmode->last_phase;
        Q_sum += Ihist[ib]*rfmode->mp_charge;
        n_summed  += Ihist[ib];
        
        if (ib==nb2) {
            Vc = rfmode->V;
            Vcr = rfmode->Vr;
            }
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
    
    if (rfmode->record) {
      if ((pass%rfmode->sample_interval)==0 && 
          (!SDDS_SetRowValues(&rfmode->SDDSrec, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,
			      (pass/rfmode->sample_interval),
			      "Pass", pass, "NumberOccupied", n_occupied,
			      "FractionBinned", np?(1.0*n_binned)/np:0.0,
			      "VPostBeam", rfmode->V, "PhasePostBeam", rfmode->last_phase,
			      "tPostBeam", rfmode->last_t,
			      "V", n_summed?V_sum/n_summed:0.0,
			      "VReal", n_summed?Vr_sum/n_summed:0.0,
			      "Phase", n_summed?phase_sum/n_summed:0.0, NULL) ||
	   !SDDS_UpdatePage(&rfmode->SDDSrec, 0))) {
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
        SDDS_Bomb("problem setting up data for RFMODE record file");
      }
      if (pass==n_passes-1 && !SDDS_Terminate(&rfmode->SDDSrec)) {
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
        SDDS_Bomb("problem writing data for RFMODE record file");
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

void set_up_rfmode(RFMODE *rfmode, char *element_name, double element_z, long n_passes, 
                   RUN *run, long n_particles,
                   double Po, double total_length)
{
  long n;
  double T;

  if (rfmode->initialized)
    return;
  
  rfmode->initialized = 1;
  if (rfmode->pass_interval<=0)
    bomb("pass_interval <= 0 for RFMODE", NULL);
  if (n_particles<1)
    bomb("too few particles in set_up_rfmode()", NULL);
  if (rfmode->n_bins<2)
    bomb("too few bins for RFMODE", NULL);
  if (rfmode->bin_size<=0)
    bomb("bin_size must be positive for RFMODE", NULL);
  if (rfmode->Ra && rfmode->Rs) 
    bomb("RFMODE element may have only one of Ra or Rs nonzero.  Ra is just 2*Rs", NULL);
  if (!rfmode->Ra)
    rfmode->Ra = 2*rfmode->Rs;
  if (rfmode->bin_size*rfmode->freq>0.1) {
    T = rfmode->bin_size*rfmode->n_bins;
    rfmode->bin_size = 0.1/rfmode->freq;
    rfmode->n_bins = T/rfmode->bin_size+1;
    rfmode->bin_size = T/rfmode->n_bins;
    fprintf(stdout, "The RFMODE %s bin size is too large--setting to %e and increasing to %ld bins\n",
            element_name, rfmode->bin_size, rfmode->n_bins);
    fprintf(stdout, "Total span changed from %le to %le\n",
            T, rfmode->n_bins*rfmode->bin_size);
    fflush(stdout);
  }
  if (rfmode->sample_interval<=0)
    rfmode->sample_interval = 1;
  if (rfmode->record) {
    rfmode->record = compose_filename(rfmode->record, run->rootname);
    n = n_passes/rfmode->sample_interval;
    if (!SDDS_InitializeOutput(&rfmode->SDDSrec, SDDS_BINARY, 1, NULL, NULL, rfmode->record) ||
        !SDDS_DefineSimpleColumn(&rfmode->SDDSrec, "Pass", NULL, SDDS_LONG) ||
        !SDDS_DefineSimpleColumn(&rfmode->SDDSrec, "NumberOccupied", NULL, SDDS_LONG) ||
        !SDDS_DefineSimpleColumn(&rfmode->SDDSrec, "FractionBinned", NULL, SDDS_DOUBLE) ||
        !SDDS_DefineSimpleColumn(&rfmode->SDDSrec, "V", NULL, SDDS_DOUBLE) ||
        !SDDS_DefineSimpleColumn(&rfmode->SDDSrec, "VReal", NULL, SDDS_DOUBLE) ||
        !SDDS_DefineSimpleColumn(&rfmode->SDDSrec, "Phase", NULL, SDDS_DOUBLE) ||
        !SDDS_DefineSimpleColumn(&rfmode->SDDSrec, "VPostBeam", NULL, SDDS_DOUBLE) ||
        !SDDS_DefineSimpleColumn(&rfmode->SDDSrec, "PhasePostBeam", NULL, SDDS_DOUBLE) ||
        !SDDS_DefineSimpleColumn(&rfmode->SDDSrec, "tPostBeam", NULL, SDDS_DOUBLE) ||
        !SDDS_WriteLayout(&rfmode->SDDSrec) ||
        !SDDS_StartPage(&rfmode->SDDSrec, n+1)) {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      SDDS_Bomb("problem setting up RFMODE record file");
    }
  }
  if (rfmode->preload && rfmode->charge) {
    double Vb, omega, To, tau;
    COMPLEX Vc;
    To = total_length/(Po*c_mks/sqrt(sqr(Po)+1));
    omega = rfmode->freq*PIx2;
    Vb = 2 * omega/4*(rfmode->Ra)/rfmode->Q * rfmode->charge * rfmode->preload_factor;
    tau = 2*rfmode->Q/(omega*(1+rfmode->beta));

    Vc = cmulr(
               cdiv(cassign(1, 0), 
                    cadd(cassign(1, 0), 
                         cmulr(cexpi(omega*To), -exp(-To/tau)))),
               -Vb);
    rfmode->V = sqrt(sqr(Vc.r)+sqr(Vc.i));
    rfmode->last_phase = atan2(Vc.i, Vc.r);
    fprintf(stdout, "RFMODE %s at z=%fm preloaded:  V = (%e, %e) V  =  %eV at %fdeg \n",
            element_name, element_z, Vc.r, Vc.i,
            rfmode->V, rfmode->last_phase*180/PI);
    fflush(stdout);
    fprintf(stdout, "omega=%21.15e To=%21.15es, Vb = %21.15eV, tau = %21.15es\n", omega, To, Vb, tau);
    fflush(stdout);
    rfmode->last_t = element_z/c_mks - To;
  }
  else if (rfmode->initial_V) {
    rfmode->last_phase = rfmode->initial_phase;
    rfmode->V  = rfmode->initial_V;
    rfmode->last_t = rfmode->initial_t;
  }
  else 
    rfmode->last_t = element_z/c_mks;
}

