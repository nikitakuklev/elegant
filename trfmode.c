/* Copyright 1994 by Michael Borland and Argonne National Laboratory,
 * all rights reserved.
 */
/* file: trfmode.c
 * contents: track_through_trfmode()
 *
 * Michael Borland, 1993
 */
#include "mdb.h"
#include "track.h"

void track_through_trfmode(
    double **part, long np, TRFMODE *trfmode, double Po,
    char *element_name, double element_z, long pass, long n_passes
    )
{
    static long *Ihist = NULL;               /* array for histogram of particle density */
    static double *xave = NULL;              /* average x coordinate in each bin */
    static double *yave = NULL;              /* average y coordinate in each bin */
    static double *Vxbin = NULL;             /* array for voltage acting on each bin (MV/m) */
    static double *Vybin = NULL;             /* array for voltage acting on each bin (MV/m) */
    static long max_n_bins = 0;
    static long *pbin = NULL;                /* array to record which bin each particle is in */
    static double *time = NULL;              /* array to record arrival time of each particle */
    static long max_np = 0;
    long ip, ib, nb2;
    double tmin, tmax, tmean, dt, P;
    double Vxb, Vyb, V, omega, phase, t, k, damping_factor, tau;
    double V_sum, Vr_sum, phase_sum;
    double Vc, Vcr, Q_sum;
    double Px, Py, Pz;
    long n_summed, max_hist, n_occupied;
    static long been_warned = 0;

    log_entry("track_through_trfmode");

    if (!been_warned) {        
        if (trfmode->freq<1e3 && trfmode->freq)  {
            printf("\7\7\7warning: your TRFMODE frequency is less than 1kHz--this may be an error\n");
            been_warned = 1;
            }
        if (been_warned) {
            printf("units of parameters for TRFMODE are as follows:\n");
            print_dictionary_entry(stdout, T_TRFMODE);
            }
        }

    if (!trfmode->initialized)
        bomb("track_through_trfmode called with uninitialized element", NULL);

    if (trfmode->n_bins>max_n_bins) {
       Ihist = trealloc(Ihist, sizeof(*Ihist)*(max_n_bins=trfmode->n_bins));
       xave = trealloc(xave, sizeof(*xave)*max_n_bins);
       yave = trealloc(yave, sizeof(*yave)*max_n_bins);
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
        Ihist[ib] = xave[ib] = yave[ib] = 0;

    dt = (tmax - tmin)/trfmode->n_bins;
    for (ip=0; ip<np; ip++) {
        pbin[ip] = -1;
        ib = (time[ip]-tmin)/dt;
        if (ib<0)
            continue;
        if (ib>trfmode->n_bins - 1)
            continue;
        Ihist[ib] += 1;
        xave[ib] += part[ip][0];
        yave[ib] += part[ip][2];
        pbin[ip] = ib;
        }
    for (ib=0; ib<trfmode->n_bins; ib++) 
        if (Ihist[ib]) {
            xave[ib] /= Ihist[ib];
            yave[ib] /= Ihist[ib];
            }

    omega = PIx2*trfmode->freq;
    tau = 2*trfmode->Q/(omega*(1+trfmode->beta));
    k = omega/4*(trfmode->Ra)/trfmode->Q;

    V_sum = Vr_sum = phase_sum = Vc = Vcr = Q_sum = 0;
    n_summed = max_hist = n_occupied = 0;
    nb2 = trfmode->n_bins/2;
    if (trfmode->single_pass)
        trfmode->Vx = trfmode->Vy = 0;
    for (ib=0; ib<trfmode->n_bins; ib++) {
        if (!Ihist[ib])
            continue;
        if (Ihist[ib]>max_hist)
            max_hist = Ihist[ib];
        n_occupied++;

        t = tmin+(ib+0.5)*dt;           /* middle arrival time for this bin */
        
        /* advance cavity to this time */
        /* -- x plane */
        phase = trfmode->last_xphase + omega*(t - trfmode->last_t);
        damping_factor = exp(-(t-trfmode->last_t)/tau);
        trfmode->last_xphase = phase;
        V = trfmode->Vx*damping_factor;
        trfmode->Vxr = V*cos(phase);
        trfmode->Vxi = V*sin(phase);
        /* -- y plane */
        phase = trfmode->last_yphase + omega*(t - trfmode->last_t);
        damping_factor = exp(-(t-trfmode->last_t)/tau);
        trfmode->last_yphase = phase;
        V = trfmode->Vy*damping_factor;
        trfmode->Vyr = V*cos(phase);
        trfmode->Vyi = V*sin(phase);

        trfmode->last_t = t;

        /* compute beam-induced voltage for this bin */
        /* -- x plane */
        Vxb = 2*k*trfmode->mp_charge*Ihist[ib]*xave[ib]; 
        Vxbin[ib] = trfmode->Vxr;
        /* -- y plane */
        Vyb = 2*k*trfmode->mp_charge*Ihist[ib]*yave[ib]; 
        Vybin[ib] = trfmode->Vyr;
        
        /* add beam-induced voltage to cavity voltage */
        /* -- x plane */
        trfmode->Vxi -= Vxb;
        if (trfmode->Vxi==0 && trfmode->Vxr==0)
            trfmode->last_xphase = 0;
        else
            trfmode->last_xphase = atan2(trfmode->Vxi, trfmode->Vxr);
        trfmode->Vx = sqrt(sqr(trfmode->Vxr)+sqr(trfmode->Vxi));
        /* -- y plane */
        trfmode->Vyi -= Vyb;
        if (trfmode->Vyi==0 && trfmode->Vyr==0)
            trfmode->last_yphase = 0;
        else
            trfmode->last_yphase = atan2(trfmode->Vyi, trfmode->Vxr);
        trfmode->Vy = sqrt(sqr(trfmode->Vyr)+sqr(trfmode->Vyi));

        V_sum  += Ihist[ib]*trfmode->Vx;
        Vr_sum += Ihist[ib]*trfmode->Vxr;
        n_summed  += Ihist[ib];
        
        if (ib==nb2) {
            Vc = trfmode->Vx;
            Vcr = trfmode->Vxr;
            }
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

    if (trfmode->record && (pass%trfmode->sample_interval)==0) {
        if (n_summed)
            fprintf(trfmode->fprec, "%d\t%e\t%e\t%e\t%e\t%ld\t%ld\n", 
                    pass, V_sum/n_summed, Vr_sum/n_summed, 
                    Vcr, Vc, max_hist, n_occupied);
        else
            fprintf(trfmode->fprec, "%d\t%e\t%e\t%e\t%e\t%ld\t%ld\n", 
                    pass, 0.0, 0.0,
                    Vcr, Vc, max_hist, n_occupied);
        fflush(trfmode->fprec);
        }
    log_exit("track_through_trfmode");
    }


#include "complex.h"

void set_up_trfmode(TRFMODE *trfmode, char *element_name, double element_z, long n_passes, RUN *run, long n_particles)
{
    trfmode->initialized = 1;
    if (n_particles<1)
        bomb("too few particles in set_up_trfmode()", NULL);
    if (trfmode->n_bins<2)
        bomb("too few bins for TRFMODE", NULL);
    if (trfmode->bin_size<=0)
        bomb("bin_size must be positive for TRFMODE", NULL);
    if (trfmode->sample_interval<=0)
        trfmode->sample_interval = 1;
    if (trfmode->record) {
        trfmode->record = compose_filename(trfmode->record, run->rootname);
        trfmode->fprec = fopen_e(trfmode->record, "w", 0);
        fprintf(trfmode->fprec, "7 1 0\npass\nV\\V\\Average Voltage\nV$br$n\\V\\Average Real[V]\n");
        fprintf(trfmode->fprec, "V$br,c$n\\V\\Central Real[V]\nV$bc$n\\V\\Central Voltage\n");
        fprintf(trfmode->fprec, "H$bmax$n\\ \\Maximum Number of Particles in Any Bin\n");
        fprintf(trfmode->fprec, "N$bocc$n\\ \\Number of Occupied Bins\n");
        fprintf(trfmode->fprec, "TRFMODE %s at z=%fm\n\n%10d\n", element_name, element_z, n_passes/trfmode->sample_interval+1);
        }
    if (trfmode->n_bins%2==0) {
        printf("warning: number of bins for TRFMODE %s increased from %d to %d (odd number preferred)\n",
               element_name, trfmode->n_bins, trfmode->n_bins+1);
        trfmode->n_bins += 1;
        }
    trfmode->mp_charge = trfmode->charge/n_particles;
    trfmode->last_t = element_z/c_mks;
    trfmode->Vxr = trfmode->Vxi = trfmode->Vx = 0;
    trfmode->Vyr = trfmode->Vyi = trfmode->Vy = 0;
    }

