/* Copyright 1999 by Michael Borland and Argonne National Laboratory,
 * all rights reserved.
 */

#include "mdb.h"
#include "track.h"
#include "table.h"
#include "fftpackC.h"

void set_up_trwake(TRWAKE *wakeData, RUN *run, long pass, long particles, CHARGE *charge);

void track_through_trwake(double **part, long np, TRWAKE *wakeData, double Po,
                        RUN *run, long i_pass, CHARGE *charge
                        )
{
  static double *posItime[2] = {NULL, NULL};     /* array for histogram of particle density times x, y*/
  static double *Vtime = NULL;           /* array for voltage acting on each bin */
  static long max_n_bins = 0;
  static long *pbin = NULL;              /* array to record which bin each particle is in */
  static double *time = NULL;            /* array to record arrival time of each particle */
  static double *pz = NULL;
  static long max_np = 0;
  long ib, nb, n_binned, plane;
  double factor, tmin, tmean, dt;

  log_entry("track_through_wake");

  set_up_trwake(wakeData, run, i_pass, np, charge);

  nb = wakeData->n_bins;
  dt = wakeData->dt;

  if (wakeData->n_bins>max_n_bins) {
    posItime[0] = trealloc(posItime[0], sizeof(**posItime)*(max_n_bins=wakeData->n_bins));
    posItime[1] = trealloc(posItime[1], sizeof(**posItime)*(max_n_bins=wakeData->n_bins));
    Vtime = trealloc(Vtime, sizeof(*Vtime)*(max_n_bins+1));
  }

  if (np>max_np) {
    pbin = trealloc(pbin, sizeof(*pbin)*(max_np=np));
    time = trealloc(time, sizeof(*time)*max_np);
    pz = trealloc(pz, sizeof(*pz)*max_np);
  }

  /* Compute time coordinate of each particle */
  tmean = computeTimeCoordinates(time, Po, part, np);
  tmin = tmean-dt*wakeData->n_bins/2.0;
  
  n_binned = binTransverseTimeDistribution(posItime, pz, pbin, tmin, dt, nb, time, part, Po, np,
                                           wakeData->dx, wakeData->dy);
  if (n_binned!=np)
    fprintf(stdout, "warning: only %ld of %ld particles where binned (WAKE)\n", n_binned, np);
    fflush(stdout);

  for (plane=0; plane<2; plane++) {
    if (!wakeData->W[plane])
      continue;
    
    if (wakeData->smoothing)
      SavitzyGolaySmooth(posItime[plane], nb, wakeData->SGOrder, 
                         wakeData->SGHalfWidth, wakeData->SGHalfWidth, 0);

    /* Do the convolution of the particle density and the wake function,
       V(T) = Integral[W(T-t)*I(t)dt, t={-infinity, T}]
       Note that T<0 is the head of the bunch.
       For the wake, the argument is the normal convention wherein larger
       arguments are later times.
       */
    Vtime[nb] = 0;
    convolveArrays(Vtime, wakeData->n_bins, 
                   posItime[plane], wakeData->n_bins,
                   wakeData->W[plane], wakeData->wakePoints);

    factor = wakeData->macroParticleCharge*wakeData->factor;
    for (ib=0; ib<wakeData->n_bins; ib++)
      Vtime[ib] *= factor;

    /* change particle transverse momenta to reflect voltage in relevant bin */
    applyTransverseWakeKicks(part, time, pz, pbin, np, 
                             Po, plane, 
                             Vtime, nb, tmin, dt, wakeData->interpolate);
    
  }
  log_exit("track_through_wake");
}

void applyTransverseWakeKicks(double **part, double *time, double *pz, long *pbin, long np,
                              double Po, long plane,
                              double *Vtime, long nb, double tmin, double dt, 
                              long interpolate)
{
  long ip, ib, offset;
  double dt1, Vinterp;
  
  offset = 2*plane+1;  /* xp or yp index */
  for (ip=0; ip<np; ip++) {
    if ((ib=pbin[ip])>=0 && ib<=nb-1) {
      if (interpolate) {
        dt1 = time[ip]-(tmin+dt*ib); /* distance to bin center */
        if ((dt1<0 && ib) || ib==nb-1) {
          ib--;
          dt1 += dt;
        }
        Vinterp = Vtime[ib]+(Vtime[ib+1]-Vtime[ib])/dt*dt1;
      } else
        Vinterp = Vtime[ib];
      if (Vinterp)
        part[ip][offset] += Vinterp/(1e6*me_mev)/pz[ip];
    }
  }
}

void set_up_trwake(TRWAKE *wakeData, RUN *run, long pass, long particles, CHARGE *charge)
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
    bomb("n_bins must be >=2 for TRWAKE element", NULL);

  if (!wakeData->inputFile || !wakeData->tColumn || 
      (!wakeData->WxColumn && !wakeData->WyColumn) ||
      !strlen(wakeData->inputFile) || !strlen(wakeData->tColumn))
    bomb("supply inputFile, tColumn, and WxColumn and/or WyColumn for TRWAKE element", NULL);
  if (wakeData->WxColumn && !strlen(wakeData->WxColumn))
    bomb("supply valid column name for WxColumn for TRWAKE element", NULL);
  if (wakeData->WyColumn && !strlen(wakeData->WyColumn))
    bomb("supply valid column name for WyColumn for TRWAKE element", NULL);
  
  if (!SDDS_InitializeInput(&SDDSin, wakeData->inputFile) || SDDS_ReadPage(&SDDSin)!=1)
    SDDS_Bomb("unable to read TRWAKE file");
  if ((wakeData->wakePoints=SDDS_RowCount(&SDDSin))<0)
    bomb("no data in TRWAKE file", NULL);
  if (wakeData->wakePoints<2)
    bomb("too little data in TRWAKE file", NULL);
  
  if (SDDS_CheckColumn(&SDDSin, wakeData->tColumn, "s", SDDS_ANY_FLOATING_TYPE, 
                       stdout)!=SDDS_CHECK_OK)
    bomb("problem with time column for TRWAKE file---check existence, units, and type", NULL);
  if (!(wakeData->t=SDDS_GetColumnInDoubles(&SDDSin, wakeData->tColumn)))
    SDDS_Bomb("unable to read TRWAKE file");

  wakeData->W[0] = wakeData->W[1] = NULL;
  if (wakeData->WxColumn) {
    if (SDDS_CheckColumn(&SDDSin, wakeData->WxColumn, "V/C/m", SDDS_ANY_FLOATING_TYPE, 
                         stdout)!=SDDS_CHECK_OK)
      bomb("problem with Wx wake column for TRWAKE file---check existence, units, and type", NULL);
    if (!(wakeData->W[0]=SDDS_GetColumnInDoubles(&SDDSin, wakeData->WxColumn)))
      SDDS_Bomb("unable to read WAKE file");
  }
  if (wakeData->WyColumn) {
    if (SDDS_CheckColumn(&SDDSin, wakeData->WyColumn, "V/C/m", SDDS_ANY_FLOATING_TYPE, 
                         stdout)!=SDDS_CHECK_OK)
      bomb("problem with Wy wake column for TRWAKE file---check existence, units, and type", NULL);
    if (!(wakeData->W[1]=SDDS_GetColumnInDoubles(&SDDSin, wakeData->WyColumn)))
      SDDS_Bomb("unable to read TRWAKE file");
  }
  SDDS_Terminate(&SDDSin);

  if (!wakeData->W[0] && !wakeData->W[1])
    bomb("no valid wake data for TRWAKE element", NULL);
  
  find_min_max(&tmin, &tmax, wakeData->t, wakeData->wakePoints);
  if (tmin==tmax)
    bomb("no time span in TRWAKE data", NULL);
  if (tmin!=0)
    fprintf(stdout, "warning: TRWAKE function does not start at t=0.  Offseting the function!\n");
    fflush(stdout);
  wakeData->dt = (tmax-tmin)/(wakeData->wakePoints-1);
}

double computeTimeCoordinates(double *time, double Po, double **part, long np)
{
  double tmean, P;
  long ip;
  for (ip=tmean=0; ip<np; ip++) {
    P = Po*(part[ip][5]+1);
    tmean += (time[ip] = part[ip][4]*sqrt(sqr(P)+1)/(c_mks*P));
  }
  return tmean/np;
}

long binTransverseTimeDistribution(double **posItime, double *pz, long *pbin, double tmin, double dt, long nb,
                                   double *time, double **part, double Po, long np,
                                   double dx, double dy)
{
  long ip, ib, n_binned;
  for (ib=0; ib<nb; ib++)
    posItime[0][ib] = posItime[1][ib] = 0;
  for (ip=n_binned=0; ip<np; ip++) {
    pbin[ip] = -1;
    /* Bin CENTERS are at tmin+ib*dt */
    ib = (time[ip]-tmin)/dt+0.5;
    if (ib<0)
      continue;
    if (ib>nb - 1)
      continue;
    posItime[0][ib] += part[ip][0]-dx;
    posItime[1][ib] += part[ip][2]-dy;
    pbin[ip] = ib;
    pz[ip] = Po*(1+part[ip][5])/sqrt(1+sqr(part[ip][1])+sqr(part[ip][3]));
    n_binned++;
  }
  return n_binned;
}
  
