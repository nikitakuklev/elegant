/* Copyright 1994 by Michael Borland and Argonne National Laboratory,
 * all rights reserved.
 */
/* file: ramped_rfca.c
 * contents: ramped_rf_cavity()
 *
 * Michael Borland, 1992, 1993.
 */
#include "mdb.h"
#include "track.h"
#include "table.h"

void set_up_ramped_rfca(RAMPRF *ramprf);
long find_nearby_array_entry(double *entry, long n, double key);
double linear_interpolation(double *y, double *t, long n, double t0, long i);

long ramped_rf_cavity(
    double **part, long np, RAMPRF *ramprf, double P_central, double L_central, double zEnd, long pass
    )
{
    long ip, i_volt, i_phase, i_freq, i;
    double PRatio, gamma1;
    double P, dP, gamma, beta, dgamma, phase, length, volt;
    double *coord, t, t0, Ts, omega, beta_i, beta_f;
    long fixed_freq;
    static long been_warned = 0;
#if DEBUG
    static SDDS_TABLE debugTable;
    static long debugInitialized, debugCount = 0, debugLength;
#endif

    log_entry("ramped_rf_cavity");

    if (!been_warned) {        
        if (ramprf->freq<1e3 && ramprf->freq)  {
            printf("\7\7\7warning: your RAMPRF frequency is less than 1kHz--this may be an error\n");
            been_warned = 1;
            }
        if (fabs(ramprf->volt)<100 && ramprf->volt) {
            printf("\7\7\7warning: your RAMPRF voltage is less than 100V--this may be an error\n");
            been_warned = 1;
            }
        if (been_warned) {
            printf("units of parameters for RAMPRF are as follows:\n");
            print_dictionary_entry(stdout, T_RAMPRF);
            }
        }

    if (np<=0) {
        log_exit("ramped_rf_cavity");
        return(np);
        }

    length = ramprf->length;

    if (ramprf->volt==0) {
        if (ramprf->length) {
            for (ip=0; ip<np; ip++) {
                coord = part[ip];
                coord[0] += coord[1]*length;
                coord[2] += coord[3]*length;
                coord[4] += length;
                }
            }
        log_exit("ramped_rf_cavity");
        return(np);
        }

    if (!ramprf->t_Vf)
        set_up_ramped_rfca(ramprf);

    if (!ramprf->t_Vf || !ramprf->Vfactor || !ramprf->n_Vpts)
        bomb("no (valid) waveform data for RAMPRF", NULL);
    
    gamma = sqrt(sqr(P_central)+1);
    beta  = P_central/gamma;
    if (pass==0)
        ramprf->Ts = (zEnd - length)/(beta*c_mks);
    else
        ramprf->Ts += L_central/(beta*c_mks);
    t0 = ramprf->Ts;

    /* find position within voltage and phase ramp arrays */
    i_volt = find_nearby_array_entry(ramprf->t_Vf, ramprf->n_Vpts, t0);
    fixed_freq = 0;
    if (!ramprf->pwaveform)
        fixed_freq = 1;
    if (!fixed_freq) {
        i_phase = find_nearby_array_entry(ramprf->t_dP, ramprf->n_Ppts, t0);
        i_freq = find_nearby_array_entry(ramprf->t_ff, ramprf->n_fpts, t0);
        }

    volt = ramprf->volt/(1e6*me_mev)*linear_interpolation(ramprf->Vfactor, ramprf->t_Vf, ramprf->n_Vpts, t0, i_volt);
    if (!fixed_freq) {
        phase = (ramprf->phase + linear_interpolation(ramprf->dPhase, ramprf->t_dP, ramprf->n_Ppts, t0, i_phase))*PI/180.;
        omega = PIx2*ramprf->freq*linear_interpolation(ramprf->ffactor, ramprf->t_ff, ramprf->n_fpts, t0, i_freq);
        }
    else {
        omega = PIx2*ramprf->freq;
        if (ramprf->phase_reference==0) 
            ramprf->phase_reference = unused_phase_reference();
        
        switch (get_phase_reference(&phase, ramprf->phase_reference)) {
          case REF_PHASE_RETURNED:
            break;
          case REF_PHASE_NOT_SET:
          case REF_PHASE_NONEXISTENT:
            if (!ramprf->fiducial_seen) {
                unsigned long mode;
                if (!(mode = parseFiducialMode(ramprf->fiducial)))
                    bomb("invalid fiducial mode for RAMPRF element", NULL);
                t0 = findFiducialTime(part, np, zEnd-length, length/2, P_central, mode);
                ramprf->phase_fiducial = -omega*t0;
                ramprf->fiducial_seen = 1;
                }
            set_phase_reference(ramprf->phase_reference, phase=ramprf->phase_fiducial);
            break;
          default:
            bomb("unknown return value from get_phase_reference()", NULL);
            break;
            }
        }

    if (omega) {
        t0 = -ramprf->phase_fiducial/omega;
        }
    else
        t0 = 0;
    phase += ramprf->phase*PI/180.0;
    
#ifdef DEBUG
    if (!debugInitialized) {
        debugInitialized = 1;
        debugCount = 0;
        if (!SDDS_InitializeOutput(&debugTable, SDDS_BINARY, 0, NULL, NULL, "rampedrf.debug") ||
            SDDS_DefineColumn(&debugTable, "t", NULL, "s", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
            SDDS_DefineColumn(&debugTable, "phase", NULL, "rad", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
            SDDS_DefineColumn(&debugTable, "V", NULL, "V", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
            SDDS_DefineColumn(&debugTable, "freq", NULL, "Hz", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
            SDDS_DefineColumn(&debugTable, "RF", NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0)<0 ||
            !SDDS_WriteLayout(&debugTable) || !SDDS_StartTable(&debugTable, debugLength=1024)) 
            SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
        }
#endif

    for (ip=0; ip<np; ip++) {
        coord  = part[ip];
        /* apply initial drift */
        coord[0] += coord[1]*length/2;
        coord[2] += coord[3]*length/2;
        coord[4] += length/2*sqrt(1+sqr(coord[1])+sqr(coord[3]));

        /* compute energy kick */
        P     = P_central*(1+coord[5]);
        beta_i = P/(gamma=sqrt(sqr(P)+1));
        t     = coord[4]/(c_mks*beta_i);
        dgamma = volt*sin(omega*t+phase);
        
        /* apply energy kick */
        add_to_particle_energy(coord, t, P_central, dgamma);

        /* apply final drift */
        coord[0] += coord[1]*length/2;
        coord[2] += coord[3]*length/2;
        coord[4] += length/2.0*sqrt(1+sqr(coord[1])+sqr(coord[3]));

#if defined(IEEE_MATH)
        for (i=0; i<6; i++)
            if (isnan(coord[i]) || isinf(coord[i]))
                break;
        if (i!=6) {
            fprintf(stderr, "error: bad coordinate for particle %ld in RAMPRF\n", i);
            for (i=0; i<6; i++)
                fprintf(stderr, "%15.8e ", coord[i]);
            fprintf(stderr, "\nP = %15.8e  t = %15.8e\n", P, t);
            abort();
            }

#endif
        
#ifdef DEBUG
        if (ip==0) {
            double redPhase1;
            if ((debugCount+1)>debugLength && !SDDS_LengthenTable(&debugTable, (debugLength+=1024)))
                SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
            if (!SDDS_SetRowValues(&debugTable, SDDS_BY_INDEX|SDDS_PASS_BY_VALUE,
                                   debugCount, 
                                   0, t, 
                                   1, fmod(omega*(t-t0)+phase, PIx2), 2, volt, 3, omega/PIx2,
                                   4, volt*sin(omega*(t-t0)+phase), -1) ||
                !SDDS_UpdateTable(&debugTable))
                SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
            debugCount++;
            }
#endif
        }

    log_exit("ramped_rf_cavity");
    return(np);
    }

void set_up_ramped_rfca(RAMPRF *ramprf)
{
    TABLE data;
    long i;

    log_entry("set_up_ramprf");

    ramprf->Ts = 0;
    ramprf->fiducial_seen = 0;

    if (!ramprf->vwaveform)
        bomb("no voltage waveform filename given for ramprf", NULL);

    if (!get_table(&data, ramprf->vwaveform, 1, 0))
        bomb("unable to read voltage waveform for ramprf", NULL);

    if (data.n_data<=1)
        bomb("ramprf voltage waveform contains less than 2 points", NULL);

    ramprf->t_Vf    = data.c1;
    ramprf->Vfactor = data.c2;
    ramprf->n_Vpts  = data.n_data;
    for (i=0; i<ramprf->n_Vpts-1; i++)
        if (ramprf->t_Vf[i]>ramprf->t_Vf[i+1])
            bomb("time values are not monotonically increasing in ramprf voltage waveform", NULL);
    tfree(data.xlab); tfree(data.ylab); tfree(data.title); tfree(data.topline);
    data.xlab = data.ylab = data.title = data.topline = NULL;

    if ((ramprf->fwaveform && !ramprf->pwaveform) || (!ramprf->fwaveform && ramprf->pwaveform))
        bomb("you must give both a freq_waveform and a phase_waveform, or else give neither (RAMPRF)", NULL);

    if (ramprf->pwaveform) {
        if (!get_table(&data, ramprf->pwaveform, 1, 0))
            bomb("unable to read phase waveform for ramprf", NULL);
        
        if (data.n_data<=1)
            bomb("ramprf phase waveform contains less than 2 points", NULL);
        
        ramprf->t_dP    = data.c1;
        ramprf->dPhase  = data.c2;
        ramprf->n_Ppts  = data.n_data;
        for (i=0; i<ramprf->n_Ppts-1; i++)
            if (ramprf->t_dP[i]>ramprf->t_dP[i+1])
                bomb("time values are not monotonically increasing in phase ramprf waveform", NULL);
        tfree(data.xlab); tfree(data.ylab); tfree(data.title); tfree(data.topline);
        data.xlab = data.ylab = data.title = data.topline = NULL;
        }

    if (ramprf->fwaveform) {
        if (!get_table(&data, ramprf->fwaveform, 1, 0))
            bomb("unable to read frequency waveform for ramprf", NULL);
        
        if (data.n_data<=1)
            bomb("ramprf frequency waveform contains less than 2 points", NULL);
        
        ramprf->t_ff    = data.c1;
        ramprf->ffactor = data.c2;
        ramprf->n_fpts  = data.n_data;
        for (i=0; i<ramprf->n_fpts-1; i++)
            if (ramprf->t_ff[i]>ramprf->t_ff[i+1])
                bomb("time values are not monotonically increasing in frequency ramprf waveform", NULL);
        tfree(data.xlab); tfree(data.ylab); tfree(data.title); tfree(data.topline);
        data.xlab = data.ylab = data.title = data.topline = NULL;
        }

    log_exit("set_up_ramprf");
    }

long find_nearby_array_entry(double *entry, long n, double key)
{
    long  hi, mid, lo;

    lo = 0;
    hi = n-1;
    if (entry[lo]>key)
        return(-1);
    if (entry[hi]<key)
        return(n-1);
    while ((hi-lo)>1) {
        mid = (hi+lo)/2;
        if (key<entry[mid])
            hi = mid;
        else
            lo = mid;
        }
    return(lo);
    }

double linear_interpolation(double *y, double *t, long n, double t0, long i)
{
    if (i<0)
        i = 0;
    if (i>=n-1)
        i = n-2;
    while (i<=n-2 && t0>t[i+1])
        i++;
    if (i==n-1)
        return(y[n-1]);
    while (i>=0 && t0<t[i])
        i--;
    if (i==-1)
        return(y[0]);
    if (!(t0>=t[i] && t0<=t[i+1])) {
        printf("failure to bracket point in time array: t0=%e, t[0] = %e, t[n-1] = %e\n",
               t0, t[0], t[n-1]);
        abort();
        }
    return( y[i] + (y[i+1]-y[i])/(t[i+1]-t[i])*(t0-t[i]) );
    }
