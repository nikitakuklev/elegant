/* Copyright 1999 by Michael Borland and Argonne National Laboratory,
 * all rights reserved.
 */

#include "mdb.h"
#include "track.h"
#include "table.h"
#include "fftpackC.h"

void set_up_wake(WAKE *wakeData, RUN *run, long pass, long particles, CHARGE *charge);
void convolveArrays(double *output, long outputs, 
                    double *a1, long n1,
                    double *a2, long n2);


void track_through_wake(double **part, long np, WAKE *wakeData, double *PoInput,
                        RUN *run, long i_pass, CHARGE *charge
                        )
{
  static double *Itime = NULL;           /* array for histogram of particle density */
  static double *Vtime = NULL;           /* array for voltage acting on each bin */
  static long max_n_bins = 0;
  static long *pbin = NULL;              /* array to record which bin each particle is in */
  static double *time = NULL;            /* array to record arrival time of each particle */
  static long max_np = 0;
  long ib, nb, n_binned;
  double factor, tmin, tmax, tmean, dt, Po;

  set_up_wake(wakeData, run, i_pass, np, charge);
  Po = *PoInput;
  
  if (np>max_np) {
    pbin = trealloc(pbin, sizeof(*pbin)*(max_np=np));
    time = trealloc(time, sizeof(*time)*max_np);
  }

  /* Compute time coordinate of each particle */
  tmean = computeTimeCoordinates(time, Po, part, np);
  find_min_max(&tmin, &tmax, time, np);
  if ((tmax-tmin) > (wakeData->t[wakeData->wakePoints-1]-wakeData->t[0])) {
    fprintf(stderr, "The beam is longer than the longitudinal wake function.\nThis may produce unphysical results.\n");
    fprintf(stderr, "The beam length is %le s, while the wake length is %le s\n",
            tmax-tmin, wakeData->t[wakeData->wakePoints-1]-wakeData->t[0]);
    if (!wakeData->allowLongBeam) 
      exit(1);
    if (abs(tmax-tmean)<abs(tmin-tmean)) 
      tmin = tmax - (wakeData->t[wakeData->wakePoints-1]-wakeData->t[0]);
    else
      tmax = tmin + (wakeData->t[wakeData->wakePoints-1]-wakeData->t[0]);
  }

  dt = wakeData->dt;
  if (wakeData->n_bins) {
    nb = wakeData->n_bins;
    tmin = tmean-dt*nb/2.0;
  }
  else {
    nb = (tmax-tmin)/dt+3;
    tmin -= dt;
    tmax += dt;
  }
  
  if (nb>max_n_bins) {
    Itime = trealloc(Itime, 2*sizeof(*Itime)*(max_n_bins=nb));
    Vtime = trealloc(Vtime, 2*sizeof(*Vtime)*(max_n_bins+1));
  }

  n_binned = binTimeDistribution(Itime, pbin, tmin, dt, nb, time, part, Po, np);
  if (n_binned!=np) {
    fprintf(stdout, "warning: only %ld of %ld particles where binned (WAKE)\n", n_binned, np);
    fprintf(stdout, "consider setting n_bins=0 in WAKE definition to invoke autoscaling\n");
    fflush(stdout);
  }
  
  if (wakeData->smoothing && nb>=(2*wakeData->SGHalfWidth+1)) {
    if (!SavitzyGolaySmooth(Itime, nb, wakeData->SGOrder, wakeData->SGHalfWidth, wakeData->SGHalfWidth, 0)) {
      fprintf(stderr, "Problem with smoothing for WAKE element (file %s)\n",
              wakeData->inputFile);
      fprintf(stderr, "Parameters: nbins=%ld, order=%ld, half-width=%ld\n",
              nb, wakeData->SGOrder, wakeData->SGHalfWidth);
      exit(1);
    }
  }
  
  
  /* Do the convolution of the particle density and the wake function,
     V(T) = Integral[W(T-t)*I(t)dt, t={-infinity, T}]
     Note that T<0 is the head of the bunch.
     For the wake, the argument is the normal convention wherein larger
     arguments are later times.
     */
  Vtime[nb] = 0;
  convolveArrays(Vtime, nb,
                 Itime, nb, 
                 wakeData->W, wakeData->wakePoints);

  factor = wakeData->macroParticleCharge*wakeData->factor;
  for (ib=0; ib<nb; ib++)
    Vtime[ib] *= factor;
  
  applyLongitudinalWakeKicks(part, time, pbin, np, Po, 
                             Vtime, nb, tmin, dt, wakeData->interpolate);

  if (wakeData->change_p0)
    do_match_energy(part, np, PoInput, 0);
  
#if defined(MINIMIZE_MEMORY)
  free(Itime);
  free(Vtime);
  free(pbin);
  free(time);
  Itime = Vtime = time = NULL;
  pbin = NULL;
  max_n_bins =  max_np= 0;
#endif
}

void applyLongitudinalWakeKicks(double **part, double *time, long *pbin, long np, double Po,
                                double *Vtime, long nb, double tmin, double dt,
                                long interpolate)
{
  long ip, ib;
  double dt1, dgam;
  
  /* change particle momentum offsets to reflect voltage in relevant bin */
  for (ip=0; ip<np; ip++) {
    if ((ib=pbin[ip])>=0 && ib<=nb-1) {
      if (interpolate) {
        dt1 = time[ip]-(tmin+dt*ib);  /* distance to bin center */
        if ((dt1<0 && ib) || ib==nb-1) {
          ib--;
          dt1 += dt;
        }
        dgam = (Vtime[ib]+(Vtime[ib+1]-Vtime[ib])/dt*dt1)/(1e6*me_mev);
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

typedef struct {
  char *filename;
  long points;
  double *t, *W;
} WAKE_DATA;

static WAKE_DATA *storedWake = NULL;
static long storedWakes = 0;

void set_up_wake(WAKE *wakeData, RUN *run, long pass, long particles, CHARGE *charge)
{
  SDDS_DATASET SDDSin;
  double tmin, tmax;
  long iw;
  
  if (charge) {
    wakeData->macroParticleCharge = charge->macroParticleCharge;
  } else if (pass==0) {
    wakeData->macroParticleCharge = 0;
    if (particles)
      wakeData->macroParticleCharge = wakeData->charge/particles;
  }
  
  if (wakeData->initialized)
    return;
  wakeData->initialized = 1;
  wakeData->W = wakeData->t = NULL;

  if (wakeData->n_bins<2 && wakeData->n_bins!=0)
    bomb("n_bins must be >=2 or else 0 (autoscale) for WAKE element", NULL);

  if (!wakeData->inputFile || !strlen(wakeData->inputFile))
    bomb("supply inputFile for WAKE element", NULL);
  if (!wakeData->tColumn || !strlen(wakeData->tColumn))
    bomb("supply tColumn for WAKE element", NULL);
  if (!wakeData->WColumn || !strlen(wakeData->WColumn))
    bomb("supply WColumn for WAKE element", NULL);
  
  for (iw=0; iw<storedWakes; iw++) {
    if (strcmp(storedWake[iw].filename, wakeData->inputFile)==0)
      break;
  }
  
  if (iw==storedWakes) {
    /* read in a new wake */
    if (!SDDS_InitializeInput(&SDDSin, wakeData->inputFile) || SDDS_ReadPage(&SDDSin)!=1) {
      fprintf(stderr, "Error: unable to open or read WAKE file %s\n", wakeData->inputFile);
      exit(1);
    }
    if ((wakeData->wakePoints=SDDS_RowCount(&SDDSin))<0) {
      fprintf(stderr, "Error: no data in WAKE file %s\n",  wakeData->inputFile);
      exit(1);
    }
    if (wakeData->wakePoints<2) {
      fprintf(stderr, "Error: too little data in WAKE file %s\n",  wakeData->inputFile);
      exit(1);
    }
    if (SDDS_CheckColumn(&SDDSin, wakeData->tColumn, "s", SDDS_ANY_FLOATING_TYPE, 
                         stdout)!=SDDS_CHECK_OK) {
      fprintf(stderr, "Error: problem with time column in WAKE file %s.  Check existence, type, and units.\n",  wakeData->inputFile);
      exit(1);
    }
    if (!(wakeData->t=SDDS_GetColumnInDoubles(&SDDSin, wakeData->tColumn))) {
      fprintf(stderr, "Error: problem retrieving time data from WAKE file %s\n",  wakeData->inputFile);
      exit(1);
    }
    if (SDDS_CheckColumn(&SDDSin, wakeData->WColumn, "V/C", SDDS_ANY_FLOATING_TYPE, 
                         stdout)!=SDDS_CHECK_OK) {
      fprintf(stderr, "Error: problem with wake column in WAKE file %s.  Check existence, type, and units.\n",  wakeData->inputFile);
      exit(1);
    }
    if (!(wakeData->W=SDDS_GetColumnInDoubles(&SDDSin, wakeData->WColumn))) {
      fprintf(stderr, "Error: problem retrieving wake data from WAKE file %s\n",  wakeData->inputFile);
      exit(1);
    }
    SDDS_Terminate(&SDDSin);

    /* record in wake storage */
    if (!(storedWake=SDDS_Realloc(storedWake, sizeof(*storedWake)*(storedWakes+1))) || 
        !SDDS_CopyString(&storedWake[storedWakes].filename, wakeData->inputFile))
      SDDS_Bomb("Memory allocation failure (WAKE)");
    storedWake[storedWakes].t = wakeData->t;
    storedWake[storedWakes].W = wakeData->W;
    storedWake[storedWakes].points = wakeData->wakePoints;
    wakeData->isCopy = 0;
    storedWakes++;
  }
  else {
    /* point to an existing wake */
    wakeData->t = storedWake[iw].t;
    wakeData->W = storedWake[iw].W;
    wakeData->wakePoints = storedWake[iw].points;
    wakeData->isCopy = 1;
  }
  find_min_max(&tmin, &tmax, wakeData->t, wakeData->wakePoints);
  if (tmin>=tmax) {
    fprintf(stderr, "Error: zero or negative time span in WAKE file %s\n",  wakeData->inputFile);
    exit(1);
  }
  if (tmin!=0) {
    fprintf(stderr, "Error: WAKE function does not start at t=0 for file %s\n",  wakeData->inputFile);
    exit(1);
  }
  wakeData->dt = (tmax-tmin)/(wakeData->wakePoints-1);
}

void convolveArrays(double *output, long outputs, 
                    double *a1, long n1,
                    double *a2, long n2)
{
  long ib, ib1, ib2;
  for (ib=0; ib<outputs; ib++) {
    output[ib] = 0;
    ib2 = ib;
    ib1 = 0;
    if (ib2>=n2) {
      ib2 = n2-1;
      ib1 = ib-ib2;
      if (ib1>=n1)
        continue;
    }
    for (; ib1<=ib; ib1++, ib2--)
      output[ib] += a1[ib1]*a2[ib2];
  }
}

long binTimeDistribution(double *Itime, long *pbin, double tmin,
                         double dt, long nb, double *time, double **part, double Po, long np)
{
  long ib, ip, n_binned;
  
  for (ib=0; ib<nb; ib++)
    Itime[ib] = 0;

  for (ip=n_binned=0; ip<np; ip++) {
    pbin[ip] = -1;
    /* Bin CENTERS are at tmin+ib*dt */
    ib = (time[ip]-tmin)/dt+0.5;
    if (ib<0)
      continue;
    if (ib>nb - 1)
      continue;
    Itime[ib] += 1;
    pbin[ip] = ib;
    n_binned++;
  }
  return n_binned;
}
