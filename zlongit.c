/* Copyright 1994 by Michael Borland and Argonne National Laboratory,
 * all rights reserved.
 */
/* file: zlongit.c
 * contents: track_through_zlongit()
 *
 * Michael Borland, 1993
 */

/* Sign conventions in this module:
 * V = I*Z 
 * V(t) = V(w)*exp(i*w*t)  I(t) = I(w)*exp(i*w*t)
 * For inductor:   Z = i*w*L
 * For capacitor:  Z = -i/(w*C)
 * For resistor:   Z = R
 * For resonator:  Z = R/(1+i*Q*(w/wr-wr/w))
 *
 * The minus sign to get energy loss instead of energy gain is 
 * used internally and must not be included in the impedance.
 *
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

#define DEBUG 1

void set_up_zlongit(ZLONGIT *zlongit, RUN *run, long pass, long particles, CHARGE *charge);

void track_through_zlongit(double **part, long np, ZLONGIT *zlongit, double Po,
    RUN *run, long i_pass, CHARGE *charge
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

    set_up_zlongit(zlongit, run, i_pass, np, charge);
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

    tmean = computeTimeCoordinates(time, Po, part, np);
    tmin = tmean - dt*zlongit->n_bins/2.0;
    
    for (ib=0; ib<zlongit->n_bins; ib++)
        Itime[2*ib] = Itime[2*ib+1] = 0;

    for (ip=n_binned=0; ip<np; ip++) {
      pbin[ip] = -1;
      ib = (time[ip]-tmin)/dt;
      if (ib<0)
        continue;
      if (ib>nb - 1)
        continue;
      if (zlongit->area_weight && ib>1 && ib<(nb-1)) {
        double dist;
        dist = (time[ip]-((ib+0.5)*dt+tmin))/dt;
        Itime[ib] += 0.5;
        Itime[ib-1] += 0.25-0.5*dist;
        Itime[ib+1] += 0.25+0.5*dist;
      }
      else 
        Itime[ib] += 1;
      pbin[ip] = ib;
      n_binned++;
    }
    if (n_binned!=np)
        fprintf(stderr, "warning: only %ld of %ld particles where binned (ZLONGIT)\n", n_binned, np);
    if (zlongit->smoothing)
      SavitzyGolaySmooth(Itime, nb, zlongit->SGOrder, 
                         zlongit->SGHalfWidth, zlongit->SGHalfWidth, 0);
    
#if DEBUG 
    /* Output the time-binned data */
    if (1) {
      FILE *fp;
      fp = fopen("zlongit.tbin", "w");
      fprintf(fp, "SDDS1\n&column name=t type=double units=s &end\n&column name=I type=double &end\n&data mode=ascii &end\n");
      fprintf(fp, "%ld\n", nb);
      for (ib=0; ib<nb; ib++) 
        fprintf(fp, "%e %e\n",
                ib*dt+tmin, Itime[ib]*zlongit->macroParticleCharge/dt);
      fclose(fp);
    }
#endif

    /* Take the FFT of I(t) to get I(f) */
    realFFT(Itime, nb, 0);
    Ifreq = Itime;

    /* Compute V(f) = Z(f)*I(f), putting in a factor 
     * to normalize the current waveform.
     */
    Vfreq = Vtime;
    factor = zlongit->macroParticleCharge/dt;
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
      if ((ib=pbin[ip])>=0 && ib<=nb-1) {
        if (zlongit->interpolate && ib>0) {
          /* dt/2 offset is so that center of bin is location where
           * particle sees voltage for that bin only
           */
          dt1 = time[ip]-(tmin+dt/2.0+dt*ib);
          if (dt1<0)
            dt1 += dt;
          else
            ib += 1;
          if (ib<nb)
            dgam = (Vtime[ib-1]+(Vtime[ib]-Vtime[ib-1])/dt*dt1)/(1e6*me_mev);
          else
            continue;
        }
        else
          dgam = Vtime[ib]/(1e6*me_mev);
        if (dgam) {
          /* Put in minus sign here as the voltage decelerates the beam */
          add_to_particle_energy(part[ip], time[ip], Po, -dgam);
        }
      }
    }
  }

void set_up_zlongit(ZLONGIT *zlongit, RUN *run, long pass, long particles, CHARGE *charge)
{
    long i, nfreq;
    double df, t_range;
    static char associate[SDDS_MAXLINE];

    if (charge) {
      zlongit->macroParticleCharge = charge->macroParticleCharge;
      zlongit->charge = charge->macroParticleCharge*particles;
    } else if (pass==0) {
      zlongit->macroParticleCharge = 0;
      if (particles)
        zlongit->macroParticleCharge = zlongit->charge/particles;
    }

    if (zlongit->initialized)
      return ;
    zlongit->initialized = 1;

    if (zlongit->bin_size<=0)
        bomb("bin_size must be positive for ZLONGIT element", NULL);
    if (zlongit->broad_band) {
      /* compute impedance for a resonator.  Recall that I use V(t) = Vo*exp(i*w*t) convention,
       * so the impedance is Z(w) = (Ra/2)*(1 + i*T)/(1+T^2), where T=Q*(wo/w-w/wo).
       * The imaginary and real parts are positive for small w.
       */
        double term, factor1, factor2, factor;
        if (zlongit->Ra && zlongit->Rs) 
          bomb("ZLONGIT element broad-band resonator may have only one of Ra or Rs nonzero.  Ra is just 2*Rs", NULL);
        if (!zlongit->Ra)
          zlongit->Ra = 2*zlongit->Rs;
        if (zlongit->n_bins%2!=0)
            bomb("ZLONGIT element must have n_bins divisible by 2", NULL);
        if (zlongit->Zreal  || zlongit->Zimag) 
            bomb("can't specify both broad_band impedance and Z(f) files for ZLONGIT element", NULL);
        if (2/(zlongit->freq*zlongit->bin_size)<10) {
          /* want maximum frequency in Z > 10*fResonance */
          fprintf(stderr, "ZLONGIT has excessive bin size for given frequency\n", NULL);
          zlongit->bin_size = 0.2/zlongit->freq;
          fprintf(stderr, "  Bin size adjusted to %e\n", zlongit->bin_size);
        }
        /* Want frequency resolution < fResonance/10 */
        factor1 = 10/(zlongit->n_bins*zlongit->bin_size*zlongit->freq);
        /* Want frequency resolution < fResonanceWidth/10 */
        factor2 = 10/(zlongit->n_bins*zlongit->bin_size*zlongit->freq/zlongit->Q);
        factor = MAX(factor1, factor2);
        if (factor>1) {
          fprintf(stderr, "ZLONGIT has too few bins or excessive bin size for given frequency\n");
          zlongit->n_bins = pow(2, (long)(log(zlongit->n_bins*factor)/log(2)+1));
          fprintf(stderr, "  Number of bins adjusted to %ld\n",
                  zlongit->n_bins);
        }
        df = 1/(zlongit->n_bins*zlongit->bin_size)/(zlongit->freq);
        nfreq = zlongit->n_bins/2 + 1;
        fprintf(stderr, "ZLONGIT has %ld frequency points with df=%le\n",
                nfreq, df);
        zlongit->Z = tmalloc(sizeof(*(zlongit->Z))*zlongit->n_bins);
        zlongit->Z[0] = 0;
        zlongit->Z[zlongit->n_bins-1] = 0;    /* Nyquist term */
        for (i=1; i<nfreq-1; i++) {
            term = zlongit->Q*(1.0/(i*df)-i*df);
            zlongit->Z[2*i-1] =  zlongit->Ra/2/(1+sqr(term));
            zlongit->Z[2*i  ] =  zlongit->Z[2*i-1]*term;
            }
        if (0) {
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
        /* Recalculate the bin size and number of bins for consistency with the given
         * frequency range and spacing.  Try to maintain sufficient range of time binning. 
         */
        t_range = zlongit->n_bins*zlongit->bin_size;
        zlongit->n_bins = 2*(n_spect-1);
        zlongit->bin_size = 1.0/(zlongit->n_bins*df_spect);
        if (t_range>zlongit->n_bins*zlongit->bin_size) {
            fprintf(stderr, "error for ZLONGIT element:\nimpedance-spectrum-equivalent binning range not sufficient.\n");
            fprintf(stderr, "consider padding the impedance spectrum\n");
            exit(1);
            }
        fprintf(stderr, "Using Nb=%ld and dt=%le in ZLONGIT\n",
                zlongit->n_bins, zlongit->bin_size);
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
  }

