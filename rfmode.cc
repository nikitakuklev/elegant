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
#if USE_MPI 
#include "mpi.h" /* Defines the interface to MPI allowing the use of all MPI functions. */
#if USE_MPE
#include "mpe.h" /* Defines the MPE library */ 
#endif
#endif
#include <complex>
#include "mdb.h"
#include "track.h"

#ifdef __cplusplus
extern "C" {
#endif
long find_nearby_array_entry(double *entry, long n, double key);
double linear_interpolation(double *y, double *t, long n, double t0, long i);
#ifdef __cplusplus
}
#endif

void runBinlessRfMode(double **part, long np, RFMODE *rfmode, double Po,
		      char *element_name, double element_z, long pass, long n_passes,
		      CHARGE *charge);
void fillBerencABMatrices(MATRIX *A, MATRIX *B, RFMODE *rfmode, double dt);

void track_through_rfmode(
                          double **part0, long np0, RFMODE *rfmode, double Po,
                          char *element_name, double element_z, long pass, long n_passes,
                          CHARGE *charge
    )
{
    long *Ihist = NULL;               /* array for histogram of particle density */
    double *Vbin = NULL;              /* array for voltage acting on each bin */
    long max_n_bins = 0;
    long max_np = 0;
    long *pbin = NULL;                /* array to record which bin each particle is in */
    double *time0 = NULL;             /* array to record arrival time of each particle */
    double *time = NULL;              /* array to record arrival time of each particle */
    double **part = NULL;             /* particle buffer for working bucket */
    long *ibParticle = NULL;          /* array to record which bucket each particle is in */
    long **ipBucket = NULL;           /* array to record particle indices in part0 array for all particles in each bucket */
    long *npBucket = NULL;            /* array to record how many particles are in each bucket */
    long iBucket, nBuckets, np;
    static FILE *fpdeb = NULL;
    static FILE *fpdeb2 = NULL;
    double phig;
    
    long ip, ib, lastBin=0, firstBin=0, n_binned=0;
    double tmin=0, tmax, last_tmax, tmean, dt=0, P;
    double Vb, V, omega=0, phase, t, k, damping_factor, tau;
    double VPrevious, tPrevious, phasePrevious;
    double V_sum, Vr_sum, phase_sum, Vg_sum, phase_g_sum, Vc_sum;
    double Q_sum, dgamma;
    long n_summed, max_hist, n_occupied;
    static long been_warned = 0;
    double Qrp, VbImagFactor, Q=0;
    long deltaPass;
#if USE_MPI
    long nonEmptyBins = 0;
    long np_total;
#endif

    /*
    if (!fpdeb) {
      fpdeb = fopen("rfmode.sdds", "w");
      fprintf(fpdeb, "SDDS1\n");
      fprintf(fpdeb, "&column name=t type=double units=s &end\n");
      fprintf(fpdeb, "&column name=VI type=double units=V &end\n");
      fprintf(fpdeb, "&column name=VQ type=double units=V &end\n");
      fprintf(fpdeb, "&column name=V type=double units=V &end\n");
      fprintf(fpdeb, "&column name=phase type=double &end\n");
      fprintf(fpdeb, "&column name=dV type=double units=V &end\n");
      fprintf(fpdeb, "&column name=dPhase type=double &end\n");
      fprintf(fpdeb, "&data mode=ascii no_row_counts=1 &end\n");
    }
    if (!fpdeb2) {
      fpdeb2 = fopen("rfmode2.sdds", "w");
      fprintf(fpdeb2, "SDDS1\n");
      fprintf(fpdeb2, "&column name=Pass type=long &end\n");
      fprintf(fpdeb2, "&column name=t type=double units=s &end\n");
      fprintf(fpdeb2, "&column name=dt type=double units=s &end\n");
      fprintf(fpdeb2, "&column name=V type=double units=V &end\n");
      fprintf(fpdeb2, "&column name=phase type=double units=V &end\n");
      fprintf(fpdeb2, "&column name=VReal type=double units=V &end\n");
      fprintf(fpdeb2, "&data mode=ascii no_row_counts=1 &end\n");
    }
    */

    if (rfmode->binless) { /* This can't be done in parallel mode */
#if USE_MPI
    fprintf(stdout, (char*)"binless in rfmode is not supported in the current parallel version.\n");
    fprintf(stdout, (char*)"Please use serial version.\n");
    fflush(stdout);
    MPI_Barrier (MPI_COMM_WORLD);
    MPI_Abort(MPI_COMM_WORLD, T_RFMODE);
#endif
    runBinlessRfMode(part0, np0, rfmode, Po, element_name, element_z, pass, n_passes, charge);
    return;
    }
    
    if (charge) {
      rfmode->mp_charge = charge->macroParticleCharge;
    } else if (pass==0) {
      rfmode->mp_charge = 0;
      if (rfmode->charge<0)
        bombElegant((char*)"RFMODE charge parameter should be non-negative. Use change_particle to set particle charge state.", NULL);
#if (!USE_MPI) 
      if (np0)
        rfmode->mp_charge = rfmode->charge/np0;
#else
      if (notSinglePart) {
	MPI_Allreduce(&np0, &np_total, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
	if (np_total)
	  rfmode->mp_charge = rfmode->charge/np_total; 
      } else {
        if (np0)
          rfmode->mp_charge = rfmode->charge/np0;
      }
#endif
    }
#ifdef DEBUG
    printf("RFMODE: np0=%ld, charge=%le, mp_charge=%le\n", np0, rfmode->charge, rfmode->mp_charge);
#endif

    if (pass%rfmode->pass_interval)
      return;
    
    if (isMaster && (!been_warned)) {        
      if (rfmode->freq<1e3 && rfmode->freq)  {
        fprintf(stdout, (char*)"\7\7\7warning: your RFMODE frequency is less than 1kHz--this may be an error\n");
        fflush(stdout);
        been_warned = 1;
      }
      if (been_warned) {
        fprintf(stdout, (char*)"units of parameters for RFMODE are as follows:\n");
        fflush(stdout);
        print_dictionary_entry(stdout, T_RFMODE, 0, 0);
      }
    }
    
    if (rfmode->mp_charge==0) {
#ifdef DEBUG
      printf("RFMODE: mp_charge=0, returning\n");
#endif
      return ;
    }
    if (rfmode->detuned_until_pass>pass) {
      return ;
    }

    if (!rfmode->initialized)
        bombElegant((char*)"track_through_rfmode called with uninitialized element", NULL);
    if (rfmode->Ra)
      rfmode->RaInternal = rfmode->Ra;
    else
      rfmode->RaInternal = 2*rfmode->Rs;

    if (rfmode->n_bins>max_n_bins) {
      Ihist = (long*)trealloc(Ihist, sizeof(*Ihist)*(max_n_bins=rfmode->n_bins));
      Vbin = (double*)trealloc(Vbin, sizeof(*Vbin)*max_n_bins);
    }

    if (isSlave || !notSinglePart) {
#ifdef DEBUG
      printf("RFMODE: Determining bucket assignments\n");
#endif
      determine_bucket_assignments(part0, np0, (charge && rfmode->bunchedBeamMode)?charge->idSlotsPerBunch:0, Po, &time0, &ibParticle, &ipBucket, &npBucket, &nBuckets, -1);
#ifdef DEBUG
      printf("RFMODE: Done determining bucket assignments\n");
      fflush(stdout);
#endif 
    } else 
      nBuckets = 1;

    if (rfmode->fileInitialized) {
      long rowsNeeded = nBuckets*(n_passes/rfmode->sample_interval+1);
      if (pass==0 && (!SDDS_StartPage(&rfmode->SDDSrec, rowsNeeded))) {
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
        SDDS_Bomb((char*)"problem setting up RFMODE record file");
      }
      if ((pass!=0 || rowsNeeded>rfmode->SDDSrec.n_rows_allocated) && !SDDS_LengthenTable(&rfmode->SDDSrec, rowsNeeded-rfmode->SDDSrec.n_rows_allocated)) {
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
        SDDS_Bomb((char*)"problem setting up RFMODE record file");
      }
      if (pass==0)
        rfmode->sample_counter = 0;
    }

#if USE_MPI
    /* Master needs to know the number of buckets */
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allreduce(&nBuckets, &iBucket, 1, MPI_LONG, MPI_MAX, MPI_COMM_WORLD);
    if (myid==0)
      nBuckets = iBucket;
#endif
#ifdef DEBUG
    printf("RFMODE: nBuckets = %ld\n", nBuckets);
    fflush(stdout);
#endif
    
    for (iBucket=0; iBucket<nBuckets; iBucket++) {
#ifdef DEBUG
      printf("iBucket = %ld\n", iBucket);
      fflush(stdout);
#endif
#if USE_MPI
      /* Master needs to know if this bucket has particles */
      if (isSlave || !notSinglePart) {
        if (nBuckets==1)
          np = np0;
        else
          np = npBucket[iBucket];
      } else {
        np = 0;
      }
      MPI_Allreduce(&np, &np_total, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
      if (np_total==0)
        continue;
#endif
      
      if (isSlave || !notSinglePart) {
        if (nBuckets==1) {
          time = time0;
          part = part0;
          np = np0;
          pbin = (long*)trealloc(pbin, sizeof(*pbin)*(max_np=np));
        } else {
          if ((np = npBucket[iBucket])==0)
            continue;
          if (part)
            free_czarray_2d((void**)part, max_np, 7);
          part = (double**)czarray_2d(sizeof(double), np, 7);
          time = (double*)trealloc(time, sizeof(*time)*np);
          pbin = (long*)trealloc(pbin, sizeof(*pbin)*np);
          max_np = np;
          for (ip=0; ip<np; ip++) {
            time[ip] = time0[ipBucket[iBucket][ip]];
            memcpy(part[ip], part0[ipBucket[iBucket][ip]], sizeof(double)*7);
          }
        }
#ifdef DEBUG
        printf("Working on bucket %ld of %ld, %ld particles\n", iBucket, nBuckets, np);
        fflush(stdout);
#endif
        tmean = 0;
        if (isSlave) {
          for (ip=0; ip<np; ip++) {
            tmean += time[ip];
          }
        }
#if USE_MPI
        if (notSinglePart) {
          if (isSlave) {
            double t_total;
            MPI_Allreduce(&np, &np_total, 1, MPI_LONG, MPI_SUM, workers);
            MPI_Allreduce(&tmean, &t_total, 1, MPI_DOUBLE, MPI_SUM, workers);
            tmean = t_total;
          }
          tmean /= np_total;
        } else
          tmean /= np;
#else
        tmean /= np;
#endif
#ifdef DEBUG
        printf("computed tmean = %21.15le\n", tmean);
        fflush(stdout);
#endif
        tmin = tmean - rfmode->bin_size*rfmode->n_bins/2.;
        tmax = tmean + rfmode->bin_size*rfmode->n_bins/2.;
        if (iBucket>0 && tmin<last_tmax) {
#if USE_MPI
          if (myid==0)
#endif
            bombElegant("Error: time range overlap between buckets\n", NULL);
        }
        last_tmax = tmax;

        if (rfmode->driveFrequency>0) {
          /* handle generator voltage, cavity feedback */
          if (!rfmode->fbRunning) {
            /* This next statement phases the generator to the first bunch at the desired phase */
            rfmode->tGenerator = rfmode->fbLastTickTime = tmean - 0.5/rfmode->driveFrequency;
            rfmode->fbNextTickTime = rfmode->fbLastTickTime + rfmode->updateInterval/rfmode->driveFrequency;
            rfmode->fbRunning = 1;
          }
          if (tmean > rfmode->fbNextTickTime) {
            /* Need to advance the generator phasors to the next sample time before handling this bunch */
            long nTicks;
            double IgAmp, IgPhase, VI, VQ, omegaDrive, omegaRes, dt;
            
            nTicks = (tmean-rfmode->fbLastTickTime)/(rfmode->updateInterval/rfmode->driveFrequency)-0.5;
            while (nTicks--) {
              /* Update the voltage using the cavity state-space model */
              m_mult(rfmode->Mt1, rfmode->A, rfmode->Viq);
              m_mult(rfmode->Mt2, rfmode->B, rfmode->Iiq);
              m_add(rfmode->Viq, rfmode->Mt1, rfmode->Mt2);

              rfmode->fbLastTickTime = rfmode->fbNextTickTime;
              rfmode->fbNextTickTime += rfmode->updateInterval/rfmode->driveFrequency;

              /** Do feedback **/

              /* Calculate the net voltage and phase at this time */
              omegaRes = PIx2*rfmode->freq;
              omegaDrive = PIx2*rfmode->driveFrequency;
              /* - Calculate beam-induced voltage components. */
              dt = rfmode->fbLastTickTime - rfmode->last_t;
              damping_factor = exp(-dt/tau);
              /*
              VI = damping_factor*(rfmode->Vr*cos((omegaRes-omegaDrive)*dt) - rfmode->Vi*sin((omegaRes-omegaDrive)*dt));
              VQ = damping_factor*(rfmode->Vr*sin((omegaRes-omegaDrive)*dt) + rfmode->Vi*cos((omegaRes-omegaDrive)*dt));
              */
              VI = damping_factor*(rfmode->Vr*cos(omegaDrive*dt) - rfmode->Vi*sin(omegaDrive*dt));
              VQ = damping_factor*(rfmode->Vr*sin(omegaDrive*dt) + rfmode->Vi*cos(omegaDrive*dt));
              /* - Add generator voltage components */
              dt = rfmode->fbLastTickTime - rfmode->tGenerator;
              VI += rfmode->Viq->a[0][0]*cos(omegaDrive*dt) - rfmode->Viq->a[1][0]*sin(omegaDrive*dt);
              VQ += rfmode->Viq->a[0][0]*sin(omegaDrive*dt) + rfmode->Viq->a[1][0]*cos(omegaDrive*dt);

              /* - Compute total voltage amplitude and phase */
              V = sqrt(VI*VI+VQ*VQ);
              phase = atan2(VQ, VI);
              
              /* Calculate updated generator amplitude and phase
                 - Compute errors for voltage amplitude and phase
                 - Run these through the IIR filters
                 - Add to nominal generator amplitude and phase
              */

              /*
                fprintf(fpdeb, "%21.15le %le %le %le %le %le %le\n",
                      rfmode->fbLastTickTime, VI, VQ, V, phase, rfmode->voltageSetpoint-V, rfmode->phaseg - phase);
              */

              IgAmp = sqrt(sqr(rfmode->Ig0->a[0][0])+sqr(rfmode->Ig0->a[1][0]))
                + applyIIRFilter(rfmode->amplitudeFilter, rfmode->nAmplitudeFilters, rfmode->lambdaA*(rfmode->voltageSetpoint - V));
              IgPhase = atan2(rfmode->Ig0->a[1][0], rfmode->Ig0->a[0][0]) 
                + applyIIRFilter(rfmode->phaseFilter, rfmode->nPhaseFilters, rfmode->phaseg - phase);

              /* Calculate updated I/Q components for generator current */
              rfmode->Iiq->a[0][0] = IgAmp*cos(IgPhase);
              rfmode->Iiq->a[1][0] = IgAmp*sin(IgPhase);
            }
          }
        }

        if (isSlave) {
#ifdef DEBUG
          printf("tmin = %21.15le, tmax = %21.15le, tmean = %21.15le\n", tmin, tmax, tmean);
          fflush(stdout);
#endif
          
          for (ib=0; ib<rfmode->n_bins; ib++)
            Ihist[ib] = 0;
        
          dt = (tmax - tmin)/rfmode->n_bins;
          n_binned = lastBin = 0;
          firstBin = rfmode->n_bins;
#if USE_MPI
          nonEmptyBins = 0;
#endif
          for (ip=0; ip<np; ip++) {
            pbin[ip] = -1;
            ib = (long)((time[ip]-tmin)/dt);
            if (ib<0)
              continue;
            if (ib>rfmode->n_bins - 1)
              continue;
            Ihist[ib] += 1;
#if USE_MPI
            if (Ihist[ib]==1)
              nonEmptyBins++;
#endif
            pbin[ip] = ib;
            if (ib>lastBin)
              lastBin = ib;
            if (ib<firstBin)    
              firstBin = ib;
            n_binned++;
          }
          if (n_binned!=np) {
#if USE_MPI
            dup2(fd,fileno(stdout)); 
            printf("%ld of %ld particles outside of binning region in RFMODE. Consider increasing number of bins.\n", 
                   np-n_binned, np);
            printf("Also, check particleID assignments for bunch identification. Bunches should be on separate pages of the input file.\n");
            fflush(stdout);
            close(fd);
            mpiAbort = 1;
            MPI_Abort(MPI_COMM_WORLD, T_RFMODE);
#else 
            bombElegant("some particles  outside of binning region in RFMODE. Consider increasing number of bins. Also, particleID assignments should be checked.", NULL);
#endif
          }
          V_sum = Vr_sum = phase_sum = Q_sum = Vg_sum = Vc_sum = phase_g_sum = 0;
          n_summed = max_hist = n_occupied = 0;
    
          /* find frequency and Q at this time */
          omega = PIx2*rfmode->freq;
          if (rfmode->nFreq) {
            double omegaFactor;
            ib = find_nearby_array_entry(rfmode->tFreq, rfmode->nFreq, tmean);
            omega *= (omegaFactor=linear_interpolation(rfmode->fFreq, rfmode->tFreq, rfmode->nFreq, tmean, ib));
            /* keeps stored energy constant for constant R/Q */
            rfmode->V *= sqrt(omegaFactor);
          }
          Q = rfmode->Q/(1+rfmode->beta);
          if (rfmode->nQ) {
            ib = find_nearby_array_entry(rfmode->tQ, rfmode->nQ, tmean);
            Q *= linear_interpolation(rfmode->fQ, rfmode->tQ, rfmode->nQ, tmean, ib);
          }
        }
#if (!USE_MPI)
        if (Q<0.5) {
          fprintf(stdout, (char*)"The effective Q<=0.5 for RFMODE.  Use the ZLONGIT element.\n");
          fflush(stdout);
          exitElegant(1);
        }
#else
        if (myid == 1) { /* Let the first slave processor write the output */
          if (Q<0.5) {
            dup2(fd,fileno(stdout)); 
            fprintf(stdout, (char*)"The effective Q<=0.5 for RFMODE.  Use the ZLONGIT element.\n");
            fflush(stdout);
            close(fd);
            MPI_Abort(MPI_COMM_WORLD, T_RFMODE);
          }
        }
#endif
        tau = 2*Q/omega;
        k = omega/4*(rfmode->RaInternal)/rfmode->Q;
        if ((deltaPass = (pass-rfmode->detuned_until_pass)) <= (rfmode->rampPasses-1)) 
          k *= (deltaPass+1.0)/rfmode->rampPasses;
      
        /* These adjustments per Zotter and Kheifets, 3.2.4 */
        Qrp = sqrt(Q*Q - 0.25);
        VbImagFactor = 1/(2*Qrp);
        omega *= Qrp/Q;
      
        if (rfmode->single_pass) {
          rfmode->V = rfmode->last_phase = 0;
          rfmode->last_t = tmin + 0.5*dt;
        }
      }
    
      
#if USE_MPI
        if (nBuckets==0) {
          /* Since nBuckets can never be 0, this code is not called. */
          /* There are unresolved issues with histogram_sums() introducing noise into the results. */
          histogram_sums(nonEmptyBins, firstBin, &lastBin, Ihist);
        } else {
          long firstBin_global, lastBin_global;
          if (myid==0) {
            firstBin = rfmode->n_bins;
            lastBin = 0;
          }
          MPI_Allreduce(&lastBin, &lastBin_global, 1, MPI_LONG, MPI_MAX, MPI_COMM_WORLD);
          MPI_Allreduce(&firstBin, &firstBin_global, 1, MPI_LONG, MPI_MIN, MPI_COMM_WORLD);
          firstBin = firstBin_global;
          lastBin = lastBin_global;
          if (isSlave || !notSinglePart) { 
            double *buffer;
            buffer = (double*)calloc(lastBin-firstBin+1, sizeof(double));
            MPI_Allreduce(&Ihist[firstBin], buffer, lastBin-firstBin+1, MPI_DOUBLE, MPI_SUM, workers);
            memcpy(Ihist+firstBin, buffer, sizeof(double)*(lastBin-firstBin+1));
            free(buffer);
          }
#ifdef DEBUG
          printf("bucket = %ld, firstBin = %ld, lastBin = %ld\n",  iBucket, firstBin_global, lastBin_global);
          if (myid!=0) {
            printf("%ld particles binned\n", n_binned);
            for (ib=firstBin; ib<=lastBin; ib++) 
              printf("%ld %ld\n", ib, Ihist[ib]);
            fflush(stdout);
          }
#endif
        }
#else 
#ifdef DEBUG
        printf("%ld particles binned\n", n_binned);
        for (ib=firstBin; ib<=lastBin; ib++) 
          printf("%ld %ld\n", ib, Ihist[ib]);
        fflush(stdout);
#endif
#endif

      if (isSlave || !notSinglePart) {
        /* These values are fixed and can be used to compute the effect on the beam of
         * the "long-range" fields (previous turns) only
         */
        VPrevious = rfmode->V;
        tPrevious = rfmode->last_t;
        phasePrevious = rfmode->last_phase;

        if (rfmode->driveFrequency>0) {
          /* compute generator voltage I and Q envelopes at bunch center */
          fillBerencABMatrices(rfmode->At, rfmode->Bt, rfmode, tmean-rfmode->fbLastTickTime);
          m_mult(rfmode->Mt1, rfmode->At, rfmode->Viq);
          m_mult(rfmode->Mt2, rfmode->Bt, rfmode->Iiq);
          m_add(rfmode->Mt3, rfmode->Mt1, rfmode->Mt2);
        }
          
        for (ib=firstBin; ib<=lastBin; ib++) {
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
          Vb = 2*k*rfmode->mp_charge*particleRelSign*rfmode->pass_interval*Ihist[ib];
          if (rfmode->long_range_only)
            Vbin[ib] = VPrevious*exp(-(t-tPrevious)/tau)*cos(phasePrevious + omega*(t - tPrevious));
          else 
            Vbin[ib] = rfmode->Vr - Vb/2;

          if (rfmode->driveFrequency>0) {
            /* add generator-induced voltage */
            double Vr, Vi, dt;
            Vg_sum += Ihist[ib]*sqrt(sqr(rfmode->Mt3->a[0][0])+sqr(rfmode->Mt3->a[1][0]));
            dt = t - rfmode->tGenerator;
            Vr = rfmode->Mt3->a[0][0]*cos(PIx2*rfmode->driveFrequency*dt) - rfmode->Mt3->a[1][0]*sin(PIx2*rfmode->driveFrequency*dt);
            Vi = rfmode->Mt3->a[0][0]*sin(PIx2*rfmode->driveFrequency*dt) + rfmode->Mt3->a[1][0]*cos(PIx2*rfmode->driveFrequency*dt);
            phase_g_sum += Ihist[ib]*atan2(Vi, Vr);
            Vbin[ib] += Vr;
            Vc_sum += Ihist[ib]*sqrt(sqr(Vr+rfmode->Vr-Vb/2)+sqr(Vi+rfmode->Vi));
            /* fprintf(fpdeb2, "%ld %21.15le %21.15le %le %le %le\n", pass, t, dt, sqrt(Vi*Vi+Vr*Vr), atan2(Vi, Vr), Vr); */
          }
          
          /* add beam-induced voltage to cavity voltage */
          rfmode->Vr -= Vb;
          rfmode->Vi -= Vb*VbImagFactor;
          rfmode->last_phase = atan2(rfmode->Vi, rfmode->Vr);
          rfmode->V = sqrt(sqr(rfmode->Vr)+sqr(rfmode->Vi));
        
          V_sum  += Ihist[ib]*rfmode->V;
          Vr_sum += Ihist[ib]*rfmode->Vr;
          phase_sum += Ihist[ib]*rfmode->last_phase;
          Q_sum += Ihist[ib]*rfmode->mp_charge*particleRelSign;
          n_summed  += Ihist[ib];
        }
        /*
        fprintf(fpdeb, "\n");
        fflush(fpdeb);
        */
#ifdef DEBUG
        printf("Computed voltage values in bins\n");
        fflush(stdout);
#endif

        if (rfmode->rigid_until_pass<=pass) {
          /* change particle momentum offsets to reflect voltage in relevant bin */
          /* also recompute slopes for new momentum to conserve transverse momentum */
          for (ip=0; ip<np; ip++) {
            if (pbin[ip]>=0) {
              /* compute new momentum and momentum offset for this particle */
              dgamma = rfmode->n_cavities*Vbin[pbin[ip]]/(1e6*particleMassMV*particleRelSign);
              add_to_particle_energy(part[ip], time[ip], Po, dgamma);
            }
          }
        }

#ifdef DEBUG
        printf("Applied voltages to particles\n");
        fflush(stdout);
#endif

        if (rfmode->record) {
#if (USE_MPI)
          if (myid == 1) {
	    /* We let the first slave to dump the parameter */
#endif
#ifdef DEBUG
            printf("Writing record file\n");
            fflush(stdout);
#endif
            if ((pass%rfmode->sample_interval)==0) {
              if (!SDDS_SetRowValues(&rfmode->SDDSrec, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,
                                     rfmode->sample_counter++,
                                     (char*)"Bunch", iBucket, 
                                     (char*)"Pass", pass, (char*)"NumberOccupied", n_occupied,
                                     (char*)"FractionBinned", np?(1.0*n_binned)/np:0.0,
                                     (char*)"VPostBeam", rfmode->V, (char*)"PhasePostBeam", rfmode->last_phase,
                                     (char*)"tPostBeam", rfmode->last_t,
                                     (char*)"V", n_summed?V_sum/n_summed:0.0,
                                     (char*)"VReal", n_summed?Vr_sum/n_summed:0.0,
                                     (char*)"Phase", n_summed?phase_sum/n_summed:0.0, 
                                     (char*)"Charge", rfmode->mp_charge*np, 
                                     (char*)"VGenerator", n_summed?Vg_sum/n_summed:0.0,
                                     (char*)"PhaseGenerator", n_summed?phase_g_sum/n_summed:0.0,
                                     (char*)"VCavity", n_summed?Vc_sum/n_summed:0.0,
                                     NULL)) {
                SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
                printf("Warning: problem setting up data for RFMODE record file, row %ld\n", rfmode->sample_counter);
              }
              if (iBucket==nBuckets-1 && !SDDS_UpdatePage(&rfmode->SDDSrec, 0)) {
                SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
                printf("Warning: problem writing data for RFMODE record file, row %ld\n", rfmode->sample_counter);
              }
            }
            if (pass==n_passes-1) {
	      if (!SDDS_UpdatePage(&rfmode->SDDSrec, 0)) {
		SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
		SDDS_Bomb((char*)"problem writing data for RFMODE record file");
	      }
	    }
#ifdef DEBUG
          printf("Done writing record file\n");
          fflush(stdout);
#endif
#if USE_MPI
	  }
#endif
	}

        if (nBuckets!=1) {
          for (ip=0; ip<np; ip++)
            memcpy(part0[ipBucket[iBucket][ip]], part[ip], sizeof(double)*7);
        }
      }
#if USE_MPI
#ifdef DEBUG
      printf("Preparing to wait on barrier (1)\n");
      fflush(stdout);
#endif
      MPI_Barrier(MPI_COMM_WORLD);
#ifdef DEBUG
      printf("Finished waiting on barrier (1)\n");
      fflush(stdout);
#endif
#endif
    }

#if USE_MPI
#ifdef DEBUG
      printf("Preparing to wait on barrier (2)\n");
      fflush(stdout);
#endif
      MPI_Barrier(MPI_COMM_WORLD);
#endif
#ifdef DEBUG
      printf("Finished waiting on barrier (2)\n");
      fflush(stdout);
#endif

#ifdef DEBUG
    printf("RFMODE: Exited bunch loop\n");
    fflush(stdout);
#endif

    if (Ihist) free(Ihist);
    if (Vbin) free(Vbin);
    if (part && part!=part0)
      free_czarray_2d((void**)part, max_np, 7);
    if (time && time!=time0) 
      free(time);
    if (time0) 
      free(time0);
    if (pbin)
      free(pbin);
    if (ibParticle) 
      free(ibParticle);
    if (ipBucket)
      free_czarray_2d((void**)ipBucket, nBuckets, np0);
    if (npBucket)
      free(npBucket);
  }

void set_up_rfmode(RFMODE *rfmode, char *element_name, double element_z, long n_passes, 
                   RUN *run, long n_particles,
                   double Po, double total_length)
{
  long i, n;
  double T;
  TABLE data;
  
  if (rfmode->initialized)
    return;

  rfmode->initialized = 1;
  if (rfmode->pass_interval<=0)
    bombElegant((char*)"pass_interval <= 0 for RFMODE", NULL);
  if (rfmode->long_range_only) {
    if (rfmode->binless)
      bombElegantVA((char*)"Error: binless and long-range modes are incompatible in RFMODE (element %s)\n", element_name);
    if (rfmode->single_pass)
      bombElegantVA((char*)"single-pass and long-range modes are incompatible in RFMODE (element %s)\n", element_name);
  }      
  if (rfmode->driveFrequency>0) {
    if (rfmode->Qwaveform || rfmode->fwaveform)
      bombElegantVA("Error: Unfortunately, can't do Q_WAVEFORM or FREQ_WAVEFORM with RF feedback for RFMODE (element %s)\n", element_name);
    if (rfmode->binless)
      bombElegantVA("Error: Unfortunately, can't use BINLESS mode with RF feedback for RFMODE (element %s)\n", element_name);
  }

#if !SDDS_MPI_IO
  if (n_particles<1)
    bombElegant((char*)"too few particles in set_up_rfmode()", NULL);
#endif
  if (rfmode->n_bins<2)
    bombElegant((char*)"too few bins for RFMODE", NULL);
  if (rfmode->bin_size<=0 && !rfmode->binless)
    bombElegant((char*)"bin_size must be positive for RFMODE", NULL);
  if (rfmode->Ra && rfmode->Rs) 
    bombElegant((char*)"RFMODE element may have only one of Ra or Rs nonzero.  Ra is just 2*Rs", NULL);
  if (rfmode->Ra)
    rfmode->RaInternal = rfmode->Ra;
  else
    rfmode->RaInternal = 2*rfmode->Rs;
  if (rfmode->bin_size*rfmode->freq>0.1) {
    T = rfmode->bin_size*rfmode->n_bins;
    rfmode->bin_size = 0.1/rfmode->freq;
    rfmode->n_bins = ((long)(T/rfmode->bin_size+1));
    rfmode->bin_size = T/rfmode->n_bins;
    fprintf(stdout, (char*)"The RFMODE %s bin size is too large--setting to %e and increasing to %ld bins\n",
            element_name, rfmode->bin_size, rfmode->n_bins);
    fprintf(stdout, (char*)"Total span changed from %le to %le\n",
            T, rfmode->n_bins*rfmode->bin_size);
    fflush(stdout);
  }
  if (rfmode->sample_interval<=0)
    rfmode->sample_interval = 1;
  if (rfmode->record && !(rfmode->fileInitialized)) {
    rfmode->record = compose_filename(rfmode->record, run->rootname);
#if (USE_MPI)
    if (myid == 1) /* We let the first slave to dump the parameter */
#endif
    if (!SDDS_InitializeOutput(&rfmode->SDDSrec, SDDS_BINARY, 1, NULL, NULL, rfmode->record) ||
        !SDDS_DefineSimpleColumn(&rfmode->SDDSrec, (char*)"Bunch", NULL, SDDS_LONG) ||
        !SDDS_DefineSimpleColumn(&rfmode->SDDSrec, (char*)"Pass", NULL, SDDS_LONG) ||
        !SDDS_DefineSimpleColumn(&rfmode->SDDSrec, (char*)"NumberOccupied", NULL, SDDS_LONG) ||
        !SDDS_DefineSimpleColumn(&rfmode->SDDSrec, (char*)"FractionBinned", NULL, SDDS_DOUBLE) ||
        !SDDS_DefineSimpleColumn(&rfmode->SDDSrec, (char*)"V", NULL, SDDS_DOUBLE) ||
        !SDDS_DefineSimpleColumn(&rfmode->SDDSrec, (char*)"VReal", NULL, SDDS_DOUBLE) ||
        !SDDS_DefineSimpleColumn(&rfmode->SDDSrec, (char*)"Phase", NULL, SDDS_DOUBLE) ||
        !SDDS_DefineSimpleColumn(&rfmode->SDDSrec, (char*)"VPostBeam", NULL, SDDS_DOUBLE) ||
        !SDDS_DefineSimpleColumn(&rfmode->SDDSrec, (char*)"PhasePostBeam", NULL, SDDS_DOUBLE) ||
        !SDDS_DefineSimpleColumn(&rfmode->SDDSrec, (char*)"tPostBeam", NULL, SDDS_DOUBLE) ||
        !SDDS_DefineSimpleColumn(&rfmode->SDDSrec, (char*)"Charge", NULL, SDDS_DOUBLE) ||
        !SDDS_DefineSimpleColumn(&rfmode->SDDSrec, (char*)"VGenerator", NULL, SDDS_DOUBLE) ||
        !SDDS_DefineSimpleColumn(&rfmode->SDDSrec, (char*)"PhaseGenerator", NULL, SDDS_DOUBLE) ||
        !SDDS_DefineSimpleColumn(&rfmode->SDDSrec, (char*)"VCavity", NULL, SDDS_DOUBLE) ||
        !SDDS_WriteLayout(&rfmode->SDDSrec) ||
        !SDDS_StartPage(&rfmode->SDDSrec, n+1)) {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      SDDS_Bomb((char*)"problem setting up RFMODE record file");
    }
    rfmode->fileInitialized = 1;
  }
  if (rfmode->preload && rfmode->charge) {
    double Vb, omega, To, tau;
    std::complex <double> Vc;
    if (rfmode->fwaveform || rfmode->Qwaveform) {
      printf((char*)"Warning: preloading of RFMODE doesn't work properly with frequency or Q waveforms\n");
      printf((char*)"unless the initial values of the frequency and Q factors are 1.\n");
    }
    To = total_length/(Po*c_mks/sqrt(sqr(Po)+1));
    omega = rfmode->freq*PIx2;
    Vb = 2 * omega/4*(rfmode->Ra)/rfmode->Q * rfmode->charge * rfmode->preload_factor * particleRelSign;
    tau = 2*rfmode->Q/(omega*(1+rfmode->beta));

    Vc = -Vb/(1.0-cexpi(omega*To)*exp(-To/tau));
    rfmode->V = std::abs<double>(Vc);
    rfmode->last_phase = atan2(Vc.imag(), Vc.real());
    fprintf(stdout, (char*)"RFMODE %s at z=%fm preloaded:  V = (%e, %e) V  =  %eV at %fdeg \n",
            element_name, element_z, Vc.real(), Vc.imag(),
            rfmode->V, rfmode->last_phase*180/PI);
    fflush(stdout);
    fprintf(stdout, (char*)"omega=%21.15e To=%21.15es, Vb = %21.15eV, tau = %21.15es\n", omega, To, Vb, tau);
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
  rfmode->tFreq = rfmode->fFreq = NULL;
  rfmode->nFreq = 0;
  if (rfmode->fwaveform) {
    if (!getTableFromSearchPath(&data, rfmode->fwaveform, 1, 0)) 
      bombElegant((char*)"unable to read frequency waveform for RFMODE", NULL);
    if (data.n_data<=1)
      bombElegant((char*)"RFMODE frequency table contains less than 2 points", NULL);
    rfmode->tFreq = data.c1;
    rfmode->fFreq = data.c2;
    rfmode->nFreq = data.n_data;
    for (i=1; i<rfmode->nFreq; i++)
      if (rfmode->tFreq[i-1]>rfmode->tFreq[i])
        bombElegant((char*)"time values are not monotonically increasing in RFMODE frequency waveform", NULL);
    tfree(data.xlab); tfree(data.ylab); tfree(data.title); tfree(data.topline);
    data.xlab = data.ylab = data.title = data.topline = NULL;
    data.c1 = data.c2 = NULL;
  }
  rfmode->tQ = rfmode->fQ = NULL;
  rfmode->nQ = 0;
  if (rfmode->Qwaveform) {
    if (!getTableFromSearchPath(&data, rfmode->Qwaveform, 1, 0)) 
      bombElegant((char*)"unable to read Q waveform for RFMODE", NULL);
    if (data.n_data<=1)
      bombElegant((char*)"RFMODE Q table contains less than 2 points", NULL);
    rfmode->tQ = data.c1;
    rfmode->fQ = data.c2;
    rfmode->nQ = data.n_data;
    for (i=1; i<rfmode->nQ; i++)
      if (rfmode->tQ[i-1]>rfmode->tQ[i])
        bombElegant((char*)"time values are not monotonically increasing in RFMODE Q waveform", NULL);
    tfree(data.xlab); tfree(data.ylab); tfree(data.title); tfree(data.topline);
    data.xlab = data.ylab = data.title = data.topline = NULL;
    data.c1 = data.c2 = NULL;
  }

  if (rfmode->driveFrequency>0) {
    /* Set up the generator and related data. */
    /* See T. Berenc, RF-TN-2015-001 */
    MATRIX *C;
    double deltaOmega, decrement, sigma, QL, Tsample, k, b1, alpha, beta, sin1, cos1;

    rfmode->lambdaA = 2*(rfmode->beta+1)/rfmode->RaInternal;
    m_alloc(&rfmode->A, 2, 2);
    m_alloc(&rfmode->B, 2, 2);
    m_alloc(&rfmode->At, 2, 2);
    m_alloc(&rfmode->Bt, 2, 2);
    m_alloc(&rfmode->Mt1, 2, 1);
    m_alloc(&rfmode->Mt2, 2, 1);
    m_alloc(&rfmode->Mt3, 2, 1);
    m_alloc(&rfmode->Ig0, 2, 1);
    m_alloc(&rfmode->Viq, 2, 1);
    m_alloc(&rfmode->Iiq, 2, 1);

    fillBerencABMatrices(rfmode->A, rfmode->B, rfmode, rfmode->updateInterval/rfmode->driveFrequency);

    /* convert from V*sin(phi) to V*cos(phi) convention and account for fact that the
     * feedback is performed in the nominally empty part of each bucket (i.e., 180 degrees before 
     * the nominal beam phase)
     */
    rfmode->phaseg = PI/180*rfmode->phaseSetpoint - PI/2 - PI;
    rfmode->Viq->a[0][0] = rfmode->voltageSetpoint*cos(rfmode->phaseg);
    rfmode->Viq->a[1][0] = rfmode->voltageSetpoint*sin(rfmode->phaseg);

    /* Compute nominal generator current */
    QL = rfmode->Q/(1+rfmode->beta);
    deltaOmega = PIx2*(rfmode->freq - rfmode->driveFrequency);
    sigma = PIx2*rfmode->freq/(2*QL);
    Tsample = rfmode->updateInterval/rfmode->driveFrequency;
    decrement = exp(-sigma*Tsample);
    k = PIx2*rfmode->freq/4*(rfmode->RaInternal/rfmode->Q);
    m_alloc(&C, 2, 2);
    C->a[0][0] = C->a[1][1] = sigma/k;
    C->a[1][0] = -(C->a[0][1] = deltaOmega/k);
    m_mult(rfmode->Ig0, C, rfmode->Viq);
    m_copy(rfmode->Iiq, rfmode->Ig0);
    m_free(&C);

    rfmode->fbRunning = 0;
    rfmode->fbNextTickTime = 0;
    rfmode->fbLastTickTime = 0;
    
    /* Read FB filters */
    rfmode->nAmplitudeFilters = rfmode->nPhaseFilters = 0;
    if (rfmode->amplitudeFilterFile && !(rfmode->nAmplitudeFilters=readIIRFilter(rfmode->amplitudeFilter, 4, rfmode->amplitudeFilterFile)))
      bombElegantVA("Error: problem reading amplitude filter file for RFMODE %s\n", element_name);
    if (rfmode->phaseFilterFile && !(rfmode->nPhaseFilters=readIIRFilter(rfmode->phaseFilter, 4, rfmode->phaseFilterFile))) 
      bombElegantVA("Error: problem reading phase filter file for RFMODE %s\n", element_name);
  }
   
}

void runBinlessRfMode(
		      double **part, long np, RFMODE *rfmode, double Po,
		      char *element_name, double element_z, long pass, long n_passes,
		      CHARGE *charge
		      )
{
  static TIMEDATA *tData;
  static long max_np = 0;
  long ip, ip0, ib;
  double P;
  double Vb, V, omega, phase, t, k, damping_factor, tau;
  double V_sum, Vr_sum, phase_sum;
  double Vc, Vcr, Q_sum, dgamma, Vb_sum;
  long n_summed, max_hist, n_occupied;
  static long been_warned = 0;
  double Qrp, VbImagFactor, Q;
  double tmean;
  long deltaPass;
#if USE_MPI
  long np_total;

  if (isSlave)
    MPI_Allreduce(&np, &np_total, 1, MPI_LONG, MPI_SUM, workers);
#endif

  if (charge) {
    rfmode->mp_charge = charge->macroParticleCharge;
  } else if (pass==0) {
    rfmode->mp_charge = 0;
    if (rfmode->charge<0)
      bombElegant((char*)"RFMODE charge parameter should be non-negative. Use change_particle to set particle charge state.", NULL);
#if !USE_MPI
    if (np)
      rfmode->mp_charge = rfmode->charge/np;
#else
    if (notSinglePart) {
      if (np_total)
	rfmode->mp_charge = rfmode->charge/np_total;
    } else {
      if (np)
	rfmode->mp_charge = rfmode->charge/np;  
    }
#endif
  }

  if (pass%rfmode->pass_interval)
    return;
    
  if (!been_warned) {        
    if (rfmode->freq<1e3 && rfmode->freq)  {
      fprintf(stdout, (char*)"\7\7\7warning: your RFMODE frequency is less than 1kHz--this may be an error\n");
      fflush(stdout);
      been_warned = 1;
    }
    if (been_warned) {
      fprintf(stdout, (char*)"units of parameters for RFMODE are as follows:\n");
      fflush(stdout);
      print_dictionary_entry(stdout, T_RFMODE, 0, 0);
    }
  }

  if (rfmode->mp_charge==0) {
    return ;
  }
  if (rfmode->detuned_until_pass>pass) {
    return ;
  }

  if (!rfmode->initialized)
    bombElegant((char*)"track_through_rfmode called with uninitialized element", NULL);

  if (np>max_np) 
    tData = (TIMEDATA*)trealloc(tData, sizeof(*tData)*(max_np=np));

  tmean = 0;
  for (ip=0; ip<np; ip++) {
    P = Po*(part[ip][5]+1);
    tData[ip].t = part[ip][4]*sqrt(sqr(P)+1)/(c_mks*P);
    tData[ip].ip = ip;
    if (isSlave) {
      tmean += tData[ip].t;
    }
  }
  qsort(tData, np, sizeof(*tData), compTimeData);

#if USE_MPI
  if (notSinglePart) {
    if (isSlave) {
      double t_total;

      MPI_Allreduce(&tmean, &t_total, 1, MPI_DOUBLE, MPI_SUM, workers);
      tmean = t_total;
    }
    tmean /= np_total;
  } else
    tmean /= np;
#else
  tmean /= np;
#endif

  V_sum = Vr_sum = phase_sum = Vc = Vcr = Q_sum = Vb_sum = 0;
  n_summed = max_hist = n_occupied = 0;
    
  /* find frequency and Q at this time */
  omega = PIx2*rfmode->freq;
  if (rfmode->nFreq) {
    double omegaFactor;
    ib = find_nearby_array_entry(rfmode->tFreq, rfmode->nFreq, tmean);
    omega *= (omegaFactor=linear_interpolation(rfmode->fFreq, rfmode->tFreq, rfmode->nFreq, tmean, ib));
    /* keeps stored energy constant for constant R/Q */
    rfmode->V *= sqrt(omegaFactor);
  }
  Q = rfmode->Q/(1+rfmode->beta);
  if (rfmode->nQ) {
    ib = find_nearby_array_entry(rfmode->tQ, rfmode->nQ, tmean);
    Q *= linear_interpolation(rfmode->fQ, rfmode->tQ, rfmode->nQ, tmean, ib);
  }
  if (Q<0.5) {
    fprintf(stdout, (char*)"The effective Q<=0.5 for RFMODE.  Use the ZLONGIT element.\n");
    fflush(stdout);
    exitElegant(1);
  }
  tau = 2*Q/omega;
  k = omega/4*(rfmode->RaInternal)/rfmode->Q;
  if ((deltaPass = (pass-rfmode->detuned_until_pass)) <= (rfmode->rampPasses-1)) 
    k *= (deltaPass+1.0)/rfmode->rampPasses;

  /* These adjustments per Zotter and Kheifets, 3.2.4 */
  Qrp = sqrt(Q*Q - 0.25);
  VbImagFactor = 1/(2*Qrp);
  omega *= Qrp/Q;

  if (rfmode->single_pass) {
    rfmode->V = rfmode->last_phase = 0;
    rfmode->last_t = tData[0].t;
  }
    
  for (ip0=0; ip0<np; ip0++) {
    ip = tData[ip0].ip;
    t = tData[ip0].t;
        
    /* advance cavity to this time */
    phase = rfmode->last_phase + omega*(t - rfmode->last_t);
    damping_factor = exp(-(t-rfmode->last_t)/tau);
    rfmode->last_t = t;
    rfmode->last_phase = phase;
    V = rfmode->V*damping_factor;
    rfmode->Vr = V*cos(phase);
    rfmode->Vi = V*sin(phase);
        
    /* compute beam-induced voltage from this particle */
    Vb = 2*k*rfmode->mp_charge*particleRelSign*rfmode->pass_interval;
    /* compute voltage seen by this particle */
    V = rfmode->Vr - Vb/2;
        
    /* add beam-induced voltage to cavity voltage */
    rfmode->Vr -= Vb;
    rfmode->Vi -= Vb*VbImagFactor;
    rfmode->last_phase = atan2(rfmode->Vi, rfmode->Vr);
    rfmode->V = sqrt(sqr(rfmode->Vr)+sqr(rfmode->Vi));
        
    V_sum  += rfmode->V;
    Vb_sum += Vb;
    Vr_sum += rfmode->Vr;
    phase_sum += rfmode->last_phase;
    Q_sum += rfmode->mp_charge*particleRelSign;
    n_summed  += 1;
        
    if (rfmode->rigid_until_pass<=pass) {
      /* change the particle energy */
      dgamma = rfmode->n_cavities*V/(1e6*particleMassMV*particleRelSign);
      add_to_particle_energy(part[ip], tData[ip0].t, Po, dgamma);
    }
  }

  if (rfmode->record) {
    if ((pass%rfmode->sample_interval)==0 && 
	(!SDDS_SetRowValues(&rfmode->SDDSrec, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,
			    (pass/rfmode->sample_interval),
			    (char*)"Pass", pass, (char*)"NumberOccupied", n_occupied,
			    (char*)"FractionBinned", 1.0,
			    (char*)"VPostBeam", rfmode->V, (char*)"PhasePostBeam", rfmode->last_phase,
			    (char*)"tPostBeam", rfmode->last_t,
			    (char*)"V", n_summed?V_sum/n_summed:0.0,
			    (char*)"VReal", n_summed?Vr_sum/n_summed:0.0,
			    (char*)"Phase", n_summed?phase_sum/n_summed:0.0, 
			    NULL) ||
	 !SDDS_UpdatePage(&rfmode->SDDSrec, 0))) {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      SDDS_Bomb((char*)"problem setting up data for RFMODE record file");
    }
    if (pass==n_passes-1 && !SDDS_Terminate(&rfmode->SDDSrec)) {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      SDDS_Bomb((char*)"problem writing data for RFMODE record file");
    }
  }
 
#if defined(MINIMIZE_MEMORY)
  free(tData);
  tData = NULL;
  max_np = 0;
#endif

}


#if USE_MPI
typedef struct {
  long index;  /* Records the real index in the whole histogram */
  long sum;
} SUB_HISTOGRAM;

void histogram_sums(long nonEmptyBins, long firstBin, long *lastBin, long *his)
{
  static long *nonEmptyArray = NULL;
  long lastBin_global, firstBin_global;
  long nonEmptyBins_total = 0, offset = 0;
  long i, j, map_j, ib;
  static SUB_HISTOGRAM *subHis; /* a compressed histogram with non-zero bins only */
  static long max_nonEmptyBins_total = 0;
  MPI_Status status;

#ifdef DEBUG
  printf("histogram_sums 1\n"); fflush(stdout);
#endif

  MPI_Reduce(lastBin, &lastBin_global, 1, MPI_LONG, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&firstBin, &firstBin_global, 1, MPI_LONG, MPI_MIN, 0, MPI_COMM_WORLD);
  if (isMaster) {
    *lastBin = lastBin_global;
    firstBin = firstBin_global; 
  }
#ifdef DEBUG
  printf("histogram_sums 2\n"); fflush(stdout);
#endif

  if (!nonEmptyArray)
    nonEmptyArray =(long*) tmalloc(sizeof(*nonEmptyArray)*n_processors);

#ifdef DEBUG
  printf("histogram_sums 3\n"); fflush(stdout);
#endif

  MPI_Gather(&nonEmptyBins,1,MPI_LONG,nonEmptyArray,1,MPI_LONG,0,MPI_COMM_WORLD);
  if (isMaster){
    for (i=1; i<n_processors; i++) 
      nonEmptyBins_total += nonEmptyArray[i];
  }
  MPI_Bcast(&nonEmptyBins_total, 1, MPI_LONG, 0, MPI_COMM_WORLD);
#ifdef DEBUG
  printf("histogram_sums 4: nonEmptyBins_total = %ld\n", nonEmptyBins_total); fflush(stdout);
#endif
 
  if (nonEmptyBins_total>max_nonEmptyBins_total) {
    if (!(subHis = (SUB_HISTOGRAM*) trealloc(subHis, sizeof(*subHis)*(max_nonEmptyBins_total=nonEmptyBins_total))))
      bomb ("Memory allocation failure in track_through_ftrfmod", NULL);
  }

#ifdef DEBUG
  printf("histogram_sums 5\n"); fflush(stdout);
#endif
  if (isSlave) {
#ifdef DEBUG
    printf("histogram_sums 6, firstBin=%ld, lastBin=%ld\n", firstBin, *lastBin); fflush(stdout);
#endif
    for (i=0,ib=firstBin; ib<=*lastBin; ib++) {
      if (his[ib]){        
	subHis[i].index = ib;
	subHis[i].sum = his[ib]; 
	i++;
      }
    }
#ifdef DEBUG
    printf("histogram_sums 6.1\n"); fflush(stdout);
#endif
    MPI_Send(subHis, nonEmptyBins*sizeof(*subHis), MPI_BYTE, 0, 108, MPI_COMM_WORLD);
#ifdef DEBUG
    printf("histogram_sums 6.2\n"); fflush(stdout);
#endif
  }
  else {
#ifdef DEBUG
    printf("histogram_sums 7\n"); fflush(stdout);
#endif
    for (i=1; i<n_processors; i++) {
      if (i>1)
	offset += nonEmptyArray[i-1];
      MPI_Recv (&subHis[offset], nonEmptyArray[i]*sizeof(*subHis), MPI_BYTE, i, 108, MPI_COMM_WORLD, &status); 
      for (j=offset; j<nonEmptyArray[i]+offset; j++) {
        #define current subHis[j]
	map_j = current.index;
	his[map_j] += current.sum; 
      }
    }
#ifdef DEBUG
    printf("histogram_sums 7.1\n"); fflush(stdout);
#endif
      
    for (i=0, ib=firstBin; ib<=*lastBin; ib++) { 
      if (his[ib]) {
	subHis[i].index = ib;
	subHis[i].sum = his[ib];
	i++;  
      }
    } 
#ifdef DEBUG
    printf("histogram_sums 7.2\n"); fflush(stdout);
#endif
    /* If there are overlapped bins between different processors, the number should be less than the original
       nonEmptyBins_total */
    nonEmptyBins_total = i;  
  }
#ifdef DEBUG
  printf("histogram_sums 8\n"); fflush(stdout);
#endif
  MPI_Bcast (&nonEmptyBins_total, 1, MPI_LONG, 0, MPI_COMM_WORLD);
  MPI_Bcast (subHis, nonEmptyBins_total*sizeof(*subHis), MPI_BYTE, 0, MPI_COMM_WORLD);
  MPI_Bcast (lastBin, 1, MPI_LONG, 0, MPI_COMM_WORLD);

  if (isSlave) {
    for (i=0; i<nonEmptyBins_total; i++) {      
      his[subHis[i].index] = subHis[i].sum; 
    }
  } 
#ifdef DEBUG
  printf("histogram_sums 9\n"); fflush(stdout);
#endif
}
#endif

void fillBerencABMatrices(MATRIX *A, MATRIX *B, RFMODE *rfmode, double dt)
{
  double deltaOmega, decrement, sigma, QL, Tsample, k, b1, alpha, beta, sin1, cos1;

  QL = rfmode->Q/(1+rfmode->beta);
  deltaOmega = PIx2*(rfmode->freq - rfmode->driveFrequency);
  sigma = PIx2*rfmode->freq/(2*QL);
  Tsample = dt;
  decrement = exp(-sigma*Tsample);
  
  sin1 = sin(deltaOmega*Tsample);
  cos1 = cos(deltaOmega*Tsample);
  A->a[0][0] = A->a[1][1] = decrement*cos1;
  A->a[0][1] = -(A->a[1][0] = decrement*sin1);
  
  k = PIx2*rfmode->freq/4*(rfmode->RaInternal/rfmode->Q);
  b1 = k/(sqr(sigma) + sqr(deltaOmega));
  alpha = deltaOmega*decrement*sin1 - sigma*decrement*cos1 + sigma;
  beta = sigma*decrement*sin1 + deltaOmega*decrement*cos1 - deltaOmega;
  B->a[0][0] = B->a[1][1] = b1*alpha;
  B->a[1][0] = -(B->a[0][1] = b1*beta);
}

