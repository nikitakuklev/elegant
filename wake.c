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


void track_through_wake(double **part, long np, WAKE *wakeData, double Po,
                        RUN *run, long i_pass, CHARGE *charge
                        )
{
  static double *Itime = NULL;           /* array for histogram of particle density */
  static double *Vtime = NULL;           /* array for voltage acting on each bin */
  static long max_n_bins = 0;
  static long *pbin = NULL;              /* array to record which bin each particle is in */
  static double *time = NULL;            /* array to record arrival time of each particle */
  static long max_np = 0;
  long ip, ib, ib1, ib2, nb, n_binned;
  double factor, tmin, tmean, dt, dt1, P, dgam, gam, frac;

  set_up_wake(wakeData, run, i_pass, np, charge);
  nb = wakeData->n_bins;
  dt = wakeData->dt;

  if (nb>max_n_bins) {
    Itime = trealloc(Itime, 2*sizeof(*Itime)*(max_n_bins=nb));
    Vtime = trealloc(Vtime, 2*sizeof(*Vtime)*(max_n_bins+1));
  }

  if (np>max_np) {
    pbin = trealloc(pbin, sizeof(*pbin)*(max_np=np));
    time = trealloc(time, sizeof(*time)*max_np);
  }

  /* Compute time coordinate of each particle */
  tmean = computeTimeCoordinates(time, Po, part, np);
  tmin = tmean-dt*nb/2.0;
  
  n_binned = binTimeDistribution(Itime, pbin, tmin, dt, nb, time, part, Po, np);
  if (n_binned!=np)
    fprintf(stderr, "warning: only %ld of %ld particles where binned (WAKE)\n", n_binned, np);

  if (wakeData->smoothing)
    SavitzyGolaySmooth(Itime, nb, wakeData->SGOrder, wakeData->SGHalfWidth, wakeData->SGHalfWidth, 0);
  
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

void set_up_wake(WAKE *wakeData, RUN *run, long pass, long particles, CHARGE *charge)
{
  SDDS_DATASET SDDSin;
  double tmin, tmax;
  
  if (charge) {
    wakeData->macroParticleCharge = charge->macroParticleCharge;
    wakeData->charge = charge->macroParticleCharge*particles;
  } else if (pass==0) {
    wakeData->macroParticleCharge = 0;
    if (particles)
      wakeData->macroParticleCharge = wakeData->charge/particles;
  }
  
  if (wakeData->initialized)
    return;
  wakeData->initialized = 1;
  
  if (wakeData->n_bins<2)
    bomb("n_bins must be >=2 for WAKE element", NULL);

  if (!wakeData->inputFile || !wakeData->tColumn || !wakeData->WColumn ||
      !strlen(wakeData->inputFile) || !strlen(wakeData->tColumn) || !strlen(wakeData->WColumn))
    bomb("supply inputFile, tColumn, and WColumn for WAKE element", NULL);
  
  if (!SDDS_InitializeInput(&SDDSin, wakeData->inputFile) || SDDS_ReadPage(&SDDSin)!=1)
    SDDS_Bomb("unable to read WAKE file");
  if ((wakeData->wakePoints=SDDS_RowCount(&SDDSin))<0)
    bomb("no data in WAKE file", NULL);
  if (wakeData->wakePoints<2)
    bomb("too little data in WAKE file", NULL);
  
  if (SDDS_CheckColumn(&SDDSin, wakeData->tColumn, "s", SDDS_ANY_FLOATING_TYPE, 
                       stderr)!=SDDS_CHECK_OK)
    bomb("problem with time column for WAKE file---check existence, units, and type", NULL);
  if (!(wakeData->t=SDDS_GetColumnInDoubles(&SDDSin, wakeData->tColumn)))
    SDDS_Bomb("unable to read WAKE file");
  
  if (SDDS_CheckColumn(&SDDSin, wakeData->WColumn, "V/C", SDDS_ANY_FLOATING_TYPE, 
                       stderr)!=SDDS_CHECK_OK)
    bomb("problem with wake column for WAKE file---check existence, units, and type", NULL);
  if (!(wakeData->W=SDDS_GetColumnInDoubles(&SDDSin, wakeData->WColumn)))
    SDDS_Bomb("unable to read WAKE file");
  SDDS_Terminate(&SDDSin);
  
  find_min_max(&tmin, &tmax, wakeData->t, wakeData->wakePoints);
  if (tmin==tmax)
    bomb("no time span in WAKE data", NULL);
  if (tmin!=0)
    fprintf(stderr, "warning: WAKE function does not start at t=0.  Offseting the function!\n");
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
