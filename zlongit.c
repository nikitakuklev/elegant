/* Copyright 1994 by Michael Borland and Argonne National Laboratory,
 * all rights reserved.
 */
/* file: zlongit.c
 * contents: track_through_zlongit()
 *
 * Michael Borland, 1993
 */
#include "mdb.h"
#include "track.h"
#include "table.h"
#include "fftpackC.h"

#define WAKE_COLUMNS 2
static SDDS_DEFINITION wake_column[WAKE_COLUMNS] = {
    {"Deltat", "&column name=Deltat, symbol=\"$gD$rt\", units=s, type=double, description=\"Time after head of bunch\" &end"},
    {"Wz", "&column name=Wz, symbol=\"W$bz$n\", units=V, type=double, description=\"Longitudinal wake\" &end"},
    };

#define WAKE_PARAMETERS 5
#define BB_WAKE_PARAMETERS 5
#define NBB_WAKE_PARAMETERS 2
static SDDS_DEFINITION wake_parameter[WAKE_PARAMETERS] = {
    {"Pass", "&parameter name=Pass, type=long &end"},
    {"q", "&parameter name=q, units=C, type=double, description=\"Total charge\" &end"},
    {"Ra", "&parameter name=Ra, symbol=\"R$ba$n\", units=\"$gW$r\", type=double, description=\"Broad-band impedance\" &end"},
    {"fo", "&parameter name=fo, symbol=\"f$bo$n\", units=Hz, type=double, description=\"Frequency of BB resonator\" &end"},
    {"Deltaf", "&parameter name=Deltaf, symbol=\"$gD$rf\", units=Hz, type=double, description=\"Frequency sampling interval\" &end"},
    } ;

void track_through_zlongit(double **part, long np, ZLONGIT *zlongit, double Po,
    RUN *run, long i_pass
    )
{
    static double *Itime = NULL;           /* array for histogram of particle density */
    static double *Vtime = NULL;           /* array for voltage acting on each bin */
    static long max_n_bins = 0;
    static long *pbin = NULL;              /* array to record which bin each particle is in */
    static double *time = NULL;            /* array to record arrival time of each particle */
    static long max_np = 0;
    double *Ifreq, *Vfreq, *Z;
    long ip, ib, nb, n_binned, nfreq, iReal, iImag;
    double factor, tmin, tmax, tmean, dt, dt1, P, dgam, gam, frac;
#if defined(DEBUG)
    static long first_call = 1;
    FILE *fp;
#endif

    log_entry("track_through_zlongit");

    if (!zlongit->initialized)
        set_up_zlongit(zlongit, run);
    nb = zlongit->n_bins;
    dt = zlongit->bin_size;

    if (zlongit->n_bins>max_n_bins) {
       Itime = trealloc(Itime, 2*sizeof(*Itime)*(max_n_bins=zlongit->n_bins));
       Vtime = trealloc(Vtime, 2*sizeof(*Vtime)*(max_n_bins+1));
       }

    if (np>max_np) {
        pbin = trealloc(pbin, sizeof(*pbin)*(max_np=np));
        time = trealloc(time, sizeof(*time)*max_np);
        }

    tmin = HUGE;
    tmean = 0;
    for (ip=0; ip<np; ip++) {
        P = Po*(part[ip][5]+1);
        time[ip] = part[ip][4]*sqrt(sqr(P)+1)/(c_mks*P);
        if (time[ip]<tmin)
            tmin = time[ip];
        tmean += time[ip];
        }
    tmin -= dt;
    tmax = tmin + zlongit->bin_size*zlongit->n_bins;
    tmean /= np;

    for (ib=0; ib<zlongit->n_bins; ib++)
        Itime[2*ib] = Itime[2*ib+1] = 0;

    for (ip=n_binned=0; ip<np; ip++) {
        pbin[ip] = -1;
        ib = (time[ip]-tmin)/dt;
        if (zlongit->area_weight) {
            frac = (time[ip]-tmin)/dt-ib;
            if (ib==0) {
                ib = 1;
                frac = 0;
                }
            }
        else
            frac = 0;
        if (ib<0)
            continue;
        if (ib>nb - 1)
            continue;
        if (ib<nb-1) {
            Itime[ib]   += (1-frac);
            Itime[ib+1] += frac;
            }
        else
            Itime[ib] += 1;
        pbin[ip] = ib;
        n_binned++;
        }
    if (n_binned!=np)
        printf("warning: only %ld of %ld particles where binned (ZLONGIT)\n", n_binned, np);

    /* Take the FFT of I(t) to get I(f) */
    realFFT(Itime, nb, 0);
    Ifreq = Itime;

    /* Compute V(f) = Z(f)*I(f), putting in a factor 
     * to normalize the current waveform
     */
    Vfreq = Vtime;
    factor = -(zlongit->charge/np)/((tmax-tmin)/nb);
    Z = zlongit->Z;
    Vfreq[0] = Ifreq[0]*Z[0]*factor;
    nfreq = nb/2 + 1;
    if (nb%2==0)
        /* Nyquist term */
        Vfreq[nb-1] = Ifreq[nb-1]*Z[nb-1]*factor;
    for (ib=1; ib<nfreq-1; ib++) {
        iImag = (iReal = 2*ib-1)+1;
        Vfreq[iReal] = (Ifreq[iReal]*Z[iReal] - Ifreq[iImag]*Z[iImag])*factor;
        Vfreq[iImag] = (Ifreq[iReal]*Z[iImag] + Ifreq[iImag]*Z[iReal])*factor; 
        }

    /* Compute inverse FFT of V(f) to get V(t) */
    realFFT(Vfreq, nb, INVERSE_FFT);
    Vtime = Vfreq;

    if (zlongit->SDDS_wake_initialized && zlongit->wakes) {
        char s[100];
        /* wake potential output */
        if (zlongit->wake_interval<=0 || (i_pass%zlongit->wake_interval)==0) {
            if (!SDDS_StartTable(&zlongit->SDDS_wake, nb)) {
                SDDS_SetError("Problem starting SDDS table for wake output (track_through_zlongit)");
                SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
                }
            for (ib=0; ib<nb; ib++) {
                if (!SDDS_SetRowValues(&zlongit->SDDS_wake, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, ib,
                                       0, ib*dt, 1, Vtime[ib], -1)) {
                    SDDS_SetError("Problem setting rows of SDDS table for wake output (track_through_zlongit)");
                    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
                    }
                }
            if (!SDDS_SetParameters(&zlongit->SDDS_wake, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,
                                    "Pass", i_pass, "q", zlongit->charge, NULL)) {
                SDDS_SetError("Problem setting parameters of SDDS table for wake output (track_through_zlongit)");
                SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
                }
            if (zlongit->broad_band) {
                if (!SDDS_SetParameters(&zlongit->SDDS_wake, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,
                                        "Ra", zlongit->Ra, "fo", zlongit->freq, "Deltaf", zlongit->bin_size, NULL)) {
                    SDDS_SetError("Problem setting parameters of SDDS table for wake output (track_through_zlongit)");
                    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
                    }
                }
            if (!SDDS_WriteTable(&zlongit->SDDS_wake)) {
                SDDS_SetError("Problem writing SDDS table for wake output (track_through_zlongit)");
                SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
                }
            }
        }

    /* put zero voltage in Vtime[nb] for use in interpolation */
    Vtime[nb] = 0;
    /* change particle momentum offsets to reflect voltage in relevant bin */
    for (ip=0; ip<np; ip++) {
        if ((ib=pbin[ip])>=1 && ib<=nb-1) {
            if (zlongit->interpolate) {
                dt1 = time[ip]-(tmin+dt*ib);
                if (dt1<0)
                    dt1 = dt-dt1;
                else 
                    ib += 1;
                dgam = (Vtime[ib-1]+(Vtime[ib]-Vtime[ib-1])/dt*dt1)/(1e6*me_mev);
                }
            else
                dgam = Vtime[ib]/(1e6*me_mev);
            if (dgam) {
                P = Po*(1+part[ip][5]);
                gam = sqrt(sqr(P)+1);
                P = sqrt(sqr(gam+dgam)-1);
                part[ip][5] = (P-Po)/Po;
                part[ip][4] = time[ip]*c_mks*P/gam;
                }
            }
        }

    log_exit("track_through_zlongit");
    }



void set_up_zlongit(ZLONGIT *zlongit, RUN *run)
{
    long i, nfreq;
    double df, t_range;
    static char associate[SDDS_MAXLINE];

    log_entry("set_up_zlongit");

    zlongit->initialized = 1;
    if (zlongit->bin_size<=0)
        bomb("bin_size must be positive for ZLONGIT element", NULL);
    if (zlongit->broad_band) {
        double term;
        if (zlongit->n_bins%2!=0)
            bomb("ZLONGIT element must have n_bins divisible by 2", NULL);
        if (zlongit->Zreal  || zlongit->Zimag) 
            bomb("can't specify both broad_band impedance and Z(f) files for ZLONGIT element", NULL);
        nfreq = zlongit->n_bins/2 + 1;
        zlongit->Z = tmalloc(sizeof(*(zlongit->Z))*zlongit->n_bins);
        df = 1/(zlongit->n_bins*zlongit->bin_size)/(zlongit->freq);
        zlongit->Z[0] = 0;
        zlongit->Z[zlongit->n_bins-1] = 0;    /* Nyquist term */
        for (i=1; i<nfreq-1; i++) {
            term = zlongit->Q*(1.0/(i*df)-i*df);
            zlongit->Z[2*i-1] =  zlongit->Ra/2/(1+sqr(term));
            zlongit->Z[2*i  ] =  zlongit->Z[2*i-1]*term;
            }
        if (1) {
            FILE *fp;
            fp = fopen("zbb.debug", "w");
            fputs("SDDS1\n&column name=Index, type=long &end\n", fp);
            fputs("&column name=zReal, type=double &end\n", fp);
            fputs("&column name=zImag, type=double &end\n", fp);
            fputs("&data mode=ascii, no_row_counts=1 &end\n", fp);
            fprintf(fp, "0 %e 0\n", zlongit->Z[0]);
            for (i=1; i<nfreq-1; i++)
                fprintf(fp, "%ld %e %e\n",
                        i, zlongit->Z[2*i-1], zlongit->Z[2*i]);
            fprintf(fp, "%ld %e 0\n", i, zlongit->Z[zlongit->n_bins-1]);
            fclose(fp);
            }
        }
    else {
        TABLE Zr_data, Zi_data;
        double *Zr, *Zi;
        double df_spect;
        long n_spect;
        if (zlongit->n_bins<1)
            bomb("ZLONGIT element must have n_bins>=1", NULL);
        if (!zlongit->Zreal && !zlongit->Zimag)
            bomb("you must either give broad_band=1, or Zreal and/or Zimag (ZLONGIT)", NULL);
        if (zlongit->Zreal && !get_table(&Zr_data, zlongit->Zreal, 1, 0))
            bomb("unable to read real impedance function (ZLONGIT)", NULL);
        if (zlongit->Zimag && !get_table(&Zi_data, zlongit->Zimag, 1, 0))
            bomb("unable to read imaginary impedance function (ZLONGIT)", NULL);
        if (zlongit->Zreal && !zlongit->Zimag) {
            Zr = Zr_data.c2;
            if ((n_spect = Zr_data.n_data)<2)
                bomb("too little data in real impedance input file (ZLONGIT)", NULL);
            df_spect = Zr_data.c1[1]-Zr_data.c1[0];
            Zi = tmalloc(sizeof(*Zi)*n_spect);
            for (i=0; i<n_spect; i++)
                Zi[i] = 0;
            }
        else if (zlongit->Zimag && !zlongit->Zreal) {
            Zi = Zi_data.c2;
            if ((n_spect = Zi_data.n_data)<2)
                bomb("too little data in imaginary impedance input file (ZLONGIT)", NULL);
            df_spect = Zi_data.c1[1]-Zi_data.c1[0];
            Zr = tmalloc(sizeof(*Zr)*n_spect);
            for (i=0; i<n_spect; i++)
                Zr[i] = 0;
            }
        else if (zlongit->Zimag && zlongit->Zreal) {
            if (Zi_data.n_data!=Zr_data.n_data)
                bomb("real and imaginary impedance files have different amounts of data (ZLONGIT)", NULL);
            n_spect = Zi_data.n_data;
            df_spect = Zi_data.c1[1]-Zi_data.c1[0];
            if (df_spect!=(Zi_data.c1[1]-Zi_data.c1[0]))
                bomb("real and imaginary impedance files have different frequency spacing (ZLONGIT)", NULL);
            Zi = Zi_data.c2;
            Zr = Zr_data.c2;
            }
        if (Zi[0])
            bomb("impedance spectrum has non-zero imaginary DC term (ZLONGIT)", NULL);
        if (!power_of_2(n_spect-1))
            bomb("number of spectrum points must be 2^n+1, n>1 (ZLONGIT)", NULL);
        t_range = zlongit->n_bins*zlongit->bin_size;
        zlongit->n_bins = 2*(n_spect-1);
        zlongit->bin_size = 1.0/(zlongit->n_bins*df_spect);
        if (t_range>zlongit->n_bins*zlongit->bin_size) {
            fprintf(stderr, "error for ZLONGIT element:\nimpedance-spectrum-equivalent binning range not sufficient.\n");
            fprintf(stderr, "consider padding the impedance spectrum\n");
            exit(1);
            }
        zlongit->Z = tmalloc(sizeof(*zlongit->Z)*2*zlongit->n_bins);
        for (i=0; i<n_spect; i++) {
            if (i==0)
                /* DC term */
                zlongit->Z[0] = Zr[0];
            else if (i==n_spect-1 && zlongit->n_bins%2==0)
                /* Nyquist term */
                zlongit->Z[n_spect-1] = Zr[i];
            else {
                zlongit->Z[2*i-1] = Zr[i];
                zlongit->Z[2*i  ] = Zi[i];
                }
            }
        }

    if (zlongit->SDDS_wake_initialized && !SDDS_Terminate(&zlongit->SDDS_wake)) {
        SDDS_SetError("Problem terminating SDDS output (set_up_zlongit)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
    zlongit->SDDS_wake_initialized = 0;

    if (zlongit->wakes) {
        zlongit->wakes = compose_filename(zlongit->wakes, run->rootname);
        if (zlongit->broad_band) 
            SDDS_ElegantOutputSetup(&zlongit->SDDS_wake, zlongit->wakes, SDDS_BINARY, 1, "longitudinal wake",
                                    run->runfile, run->lattice, wake_parameter, BB_WAKE_PARAMETERS,
                                    wake_column, WAKE_COLUMNS, "set_up_zlongit", SDDS_EOS_NEWFILE|SDDS_EOS_COMPLETE);
        else {
            SDDS_ElegantOutputSetup(&zlongit->SDDS_wake, zlongit->wakes, SDDS_BINARY, 1, "longitudinal wake",
                                    run->runfile, run->lattice, wake_parameter, NBB_WAKE_PARAMETERS,
                                    wake_column, WAKE_COLUMNS, "set_up_zlongit", SDDS_EOS_NEWFILE);
            if (zlongit->Zreal) {
                sprintf(associate, 
                        "&associate filename=\"%s\", path=\"%s\", description=\"Real part of impedance spectrum\",\
 contents=\"real impedance spectrum, parent\", sdds=0 &end", 
                        zlongit->Zreal, getenv("PWD"));
                if (!SDDS_ProcessAssociateString(&zlongit->SDDS_wake, associate)) {
                    fprintf(stderr, "Unable to define SDDS associate (set_up_zlongit)--string was:%s\n", associate);
                    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
                    exit(1);
                    }
                sprintf(associate, 
                        "&associate filename=\"%s\", path=\"%s\", description=\"Imaginary part of impedance spectrum\",\
 contents=\"real impedance spectrum, parent\", sdds=0 &end", 
                        zlongit->Zimag, getenv("PWD"));
                if (!SDDS_ProcessAssociateString(&zlongit->SDDS_wake, associate)) {
                    fprintf(stderr, "Unable to define SDDS associate (set_up_zlongit)--string was:%s\n", associate);
                    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
                    exit(1);
                    }
                }
            if (!SDDS_WriteLayout(&zlongit->SDDS_wake)) {
                SDDS_SetError("Unable to write SDDS layout (set_up_zlongit)");
                SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
                exit(1);
                }
            }
        zlongit->SDDS_wake_initialized = 1;
        }
    zlongit->initialized = 1;
    log_exit("set_up_zlongit");
    }
