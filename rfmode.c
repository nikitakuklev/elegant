/* Copyright 1994 by Michael Borland and Argonne National Laboratory,
 * all rights reserved.
 */
/* file: rfmode.c
 * contents: track_through_rfmode()
 *
 * Michael Borland, 1993
 */
#include "mdb.h"
#include "track.h"

void track_through_rfmode(
    double **part, long np, RFMODE *rfmode, double Po,
    char *element_name, double element_z, long pass, long n_passes
    )
{
    static long *Ihist = NULL;               /* array for histogram of particle density */
    static double *Vbin = NULL;              /* array for voltage acting on each bin */
    static long max_n_bins = 0;
    static long *pbin = NULL;                /* array to record which bin each particle is in */
    static double *time = NULL;              /* array to record arrival time of each particle */
    static long max_np = 0;
    long ip, ib, nb2;
    double tmin, tmax, tmean, dt, P;
    double Vb, V, omega, phase, t, k, damping_factor, tau;
    double V_sum, Vr_sum, phase_sum;
    double Vc, Vcr, Q_sum, dgamma, gamma;
    long n_summed, max_hist, n_occupied;
    static long been_warned = 0;

    log_entry("track_through_rfmode");

    if (!been_warned) {        
        if (rfmode->freq<1e3 && rfmode->freq)  {
            printf("\7\7\7warning: your RFMODE frequency is less than 1kHz--this may be an error\n");
            been_warned = 1;
            }
        if (been_warned) {
            printf("units of parameters for RFMODE are as follows:\n");
            print_dictionary_entry(stdout, T_RFMODE);
            }
        }

    if (rfmode->charge==0) {
        log_exit("track_through_rfmode");
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
    for (ip=0; ip<np; ip++) {
        pbin[ip] = -1;
        ib = (time[ip]-tmin)/dt;
        if (ib<0)
            continue;
        if (ib>rfmode->n_bins - 1)
            continue;
        Ihist[ib] += 1;
        pbin[ip] = ib;
        }

    omega = PIx2*rfmode->freq;
    tau = 2*rfmode->Q/(omega*(1+rfmode->beta));
    k = omega/4*(rfmode->Ra)/rfmode->Q;
    
    V_sum = Vr_sum = phase_sum = Vc = Vcr = Q_sum = 0;
    n_summed = max_hist = n_occupied = 0;
    nb2 = rfmode->n_bins/2;
    if (rfmode->single_pass)
        rfmode->V = 0;
    for (ib=0; ib<rfmode->n_bins; ib++) {
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
        Vb = 2*k*rfmode->mp_charge*Ihist[ib]; 
        Vbin[ib] = rfmode->Vr - Vb/2;
        
        /* add beam-induced voltage to cavity voltage */
        rfmode->Vr -= Vb;
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
        double PRatio, P1;
        for (ip=0; ip<np; ip++) {
            if (pbin[ip]>=0) {
                dgamma = Vbin[pbin[ip]]/(1e6*me_mev);
                P = Po*(1+part[ip][5]);
                gamma = sqrt(P*P+1) + dgamma;
                P1 = sqrt(gamma*gamma-1);
                part[ip][5] = (P1-Po)/Po;
                part[ip][4] = time[ip]*c_mks*P1/gamma;
                /* strictly speaking should use Pz/Pz1 in the next line */
                part[ip][1] *= (PRatio=P/P1);
                part[ip][3] *= PRatio;
                }
            }
        }

    if (rfmode->record && (pass%rfmode->sample_interval)==0) {
        long longdata[3];
        double realdata[6];
        longdata[0] = pass;
        longdata[1] = max_hist;
        longdata[2] = n_occupied;
        realdata[3] = Vcr;
        realdata[4] = Vc;
        realdata[5] = Q_sum;
        if (n_summed) {
            realdata[0] = V_sum/n_summed;
            realdata[1] = Vr_sum/n_summed;
            realdata[2] = 180/PI*phase_sum/n_summed;
            } 
        else
            realdata[0] = realdata[1] = realdata[2] = 0;
        fwrite(longdata, sizeof(*longdata), 3, rfmode->fprec);
        fwrite(realdata, sizeof(*realdata), 6, rfmode->fprec);
        fflush(rfmode->fprec);
        }
    log_exit("track_through_rfmode");
    }


#include "complex.h"

void set_up_rfmode(RFMODE *rfmode, char *element_name, double element_z, long n_passes, RUN *run, long n_particles,
                   double Po, double total_length)
{
    long n;

    rfmode->initialized = 1;
    if (n_particles<1)
        bomb("too few particles in set_up_rfmode()", NULL);
    if (rfmode->n_bins<2)
        bomb("too few bins for RFMODE", NULL);
    if (rfmode->bin_size<=0)
        bomb("bin_size must be positive for RFMODE", NULL);
    if (rfmode->sample_interval<=0)
        rfmode->sample_interval = 1;
    if (rfmode->record) {
        rfmode->record = compose_filename(rfmode->record, run->rootname);
        rfmode->fprec = fopen_e(rfmode->record, "w", 0);
        fprintf(rfmode->fprec, "SDDS1\n&column name=Pass, type=long &end\n");
        fprintf(rfmode->fprec, "&column name=MaximumBin, description=\"Maximum Number of Particles in any Bin\", symbol=\"H$bmax$n\", type=long &end\n");
        fprintf(rfmode->fprec, "&column name=NumberBinned, description=\"Number of Occupied Bins\", symbol=\"N$bocc$n\", type=long &end\n");
        fprintf(rfmode->fprec, "&column name=V, description=\"Average Voltage\", type=double &end\n");
        fprintf(rfmode->fprec, "&column name=VReal, description=\"Average Real[V]\", symbol=\"V$br$n\", type=double &end\n");
        fprintf(rfmode->fprec, "&column name=Phase, description=\"Average RF Phase\", symbol=\"$gF$r$brf$n\", type=double &end\n");
        fprintf(rfmode->fprec, "&column name=VRealCentral, description=\"Central Real[V]\", symbol=\"V$brc$n\", type=double &end\n");
        fprintf(rfmode->fprec, "&column name=VCentral, description=\"Central Voltage\", symbol=\"V$bc$n\", type=double &end\n");
        fprintf(rfmode->fprec, "&column name=QBinned, description=\"Total Charge Binned\", symbol=\"Q$bbin$n\", type=double &end\n");
        fprintf(rfmode->fprec, "&parameter name=mplTopline, type=string, fixed_value=\"RFMODE %s at z=%fm\" &end\n", 
                element_name, element_z);
        fprintf(rfmode->fprec, "&data mode=binary &end\n");
        n = n_passes/rfmode->sample_interval;
        fwrite(&n, sizeof(n), 1, rfmode->fprec);
        }
    if (rfmode->n_bins%2==0) {
        printf("warning: number of bins for RFMODE %s increased from %d to %d (odd number preferred)\n",
               element_name, rfmode->n_bins, rfmode->n_bins+1);
        rfmode->n_bins += 1;
        }
    rfmode->mp_charge = rfmode->charge/n_particles;
    if (rfmode->preload && rfmode->charge) {
        double Vb, omega, To, tau;
        COMPLEX Vc;
        To = total_length/(Po*c_mks/sqrt(sqr(Po)+1));
        omega = rfmode->freq*PIx2;
        Vb = 2 * omega/4*(rfmode->Ra)/rfmode->Q * rfmode->charge * rfmode->preload_factor;
        tau = 2*rfmode->Q/(omega*(1+rfmode->beta));
        Vc = cmulr(cdiv(cassign(1, 0), cadd(cassign(1, 0), cmulr(cexpi(omega*To), -exp(-To/tau)))), -Vb*exp(-To/tau));
        rfmode->V = sqrt(sqr(Vc.r)+sqr(Vc.i));
        rfmode->last_phase = atan2(Vc.i, Vc.r);
        printf("RFMODE %s at z=%fm preloaded:  V = (%e, %e) V  =  %eV at %fdeg \n",
               element_name, element_z, Vc.r, Vc.i,
               rfmode->V, rfmode->last_phase*180/PI);
        printf("To = %es, Vb = %eV, tau = %es\n", To, Vb, tau);
        }
    else {
        /* calculate phasor for specified initial voltage--convert to V=Vo*sin(phase) convention */
        rfmode->last_phase = rfmode->initial_phase*PI/180 - PIo2;
        rfmode->V  = rfmode->initial_V;
        }
    rfmode->last_t = element_z/c_mks;
    }

