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
  double factor, tmin, tmean, tmax, dt;

  set_up_trwake(wakeData, run, i_pass, np, charge);

  if (np>max_np) {
    pbin = trealloc(pbin, sizeof(*pbin)*(max_np=np));
    time = trealloc(time, sizeof(*time)*max_np);
    pz = trealloc(pz, sizeof(*pz)*max_np);
  }

  /* Compute time coordinate of each particle */
  tmean = computeTimeCoordinates(time, Po, part, np);
  find_min_max(&tmin, &tmax, time, np);
  
  if ((tmax-tmin) > (wakeData->t[wakeData->wakePoints-1]-wakeData->t[0])) {
    fprintf(stderr, "The beam is longer than the transverse wake function.\nThis would produce unphysical results.\n");
    fprintf(stderr, "The beam length is %le s, while the wake length is %le s\n",
            tmax-tmin, wakeData->t[wakeData->wakePoints-1]-wakeData->t[0]);
    exit(1);
  }

  dt = wakeData->dt;
  if (wakeData->n_bins) {
    tmin = tmean-dt*wakeData->n_bins/2.0;
    nb = wakeData->n_bins;
  }
  else {
    nb = (tmax-tmin)/dt+3;
    tmin -= dt;
    tmax += dt;
  }

  if (tmin>tmax || nb<=0) {
    fprintf(stdout, "Problem with time coordinates in TRWAKE.  Po=%le\n", Po);
    exit(1);
  }
  
  if (nb>max_n_bins) {
    posItime[0] = trealloc(posItime[0], sizeof(**posItime)*(max_n_bins=nb));
    posItime[1] = trealloc(posItime[1], sizeof(**posItime)*(max_n_bins=nb));
    Vtime = trealloc(Vtime, sizeof(*Vtime)*(max_n_bins+1));
  }

  n_binned = binTransverseTimeDistribution(posItime, pz, pbin, tmin, dt, nb, time, part, Po, np,
                                           wakeData->dx, wakeData->dy);
  if (n_binned!=np) {
    fprintf(stdout, "warning: only %ld of %ld particles where binned (TRWAKE)\n", n_binned, np);
    fprintf(stdout, "consider setting n_bins=0 in TRWAKE definition to invoke autoscaling\n");
    fflush(stdout);
  }
  
  for (plane=0; plane<2; plane++) {
    if (!wakeData->W[plane])
      continue;
    
    if (wakeData->smoothing && 
        !SavitzyGolaySmooth(posItime[plane], nb, wakeData->SGOrder, 
                            wakeData->SGHalfWidth, wakeData->SGHalfWidth, 0)) {
      fprintf(stderr, "Problem with smoothing for TRWAKE element (file %s)\n",
              wakeData->inputFile);
      fprintf(stderr, "Parameters: nbins=%ld, order=%ld, half-width=%ld\n",
              nb, wakeData->SGOrder, wakeData->SGHalfWidth);
      exit(1);
    }

    /* Do the convolution of the particle density and the wake function,
       V(T) = Integral[W(T-t)*I(t)dt, t={-infinity, T}]
       Note that T<0 is the head of the bunch.
       For the wake, the argument is the normal convention wherein larger
       arguments are later times.
       */
    Vtime[nb] = 0;
    convolveArrays(Vtime, nb, 
                   posItime[plane], nb,
                   wakeData->W[plane], wakeData->wakePoints);

    factor = wakeData->macroParticleCharge*wakeData->factor;
    for (ib=0; ib<nb; ib++)
      Vtime[ib] *= factor;

    /* change particle transverse momenta to reflect voltage in relevant bin */
    applyTransverseWakeKicks(part, time, pz, pbin, np, 
                             Po, plane, 
                             Vtime, nb, tmin, dt, wakeData->interpolate);
    
  }
#if defined(MINIMIZE_MEMORY)
  free(posItime[0]);
  free(posItime[1]);
  free(Vtime);
  free(pbin);
  free(time);
  free(pz);
  Vtime = time = pz = posItime[0] = posItime[1] = NULL;
  pbin = NULL;
  max_n_bins = max_np = 0;
#endif
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

typedef struct {
  char *filename;
  long points;
  double *t, *Wx, *Wy;
} WAKE_DATA;

static WAKE_DATA *storedWake = NULL;
static long storedWakes = 0;

void set_up_trwake(TRWAKE *wakeData, RUN *run, long pass, long particles, CHARGE *charge)
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
  wakeData->W[0] = wakeData->W[1] = wakeData->t = NULL;
  
  if (wakeData->n_bins<2 && wakeData->n_bins!=0)
    bomb("n_bins must be >=2 or 0 (autoscale) for TRWAKE element", NULL);

  if (!wakeData->inputFile || !wakeData->tColumn || 
      (!wakeData->WxColumn && !wakeData->WyColumn) ||
      !strlen(wakeData->inputFile) || !strlen(wakeData->tColumn))
    bomb("supply inputFile, tColumn, and WxColumn and/or WyColumn for TRWAKE element", NULL);
  if (wakeData->WxColumn && !strlen(wakeData->WxColumn))
    bomb("supply valid column name for WxColumn for TRWAKE element", NULL);
  if (wakeData->WyColumn && !strlen(wakeData->WyColumn))
    bomb("supply valid column name for WyColumn for TRWAKE element", NULL);
  
  for (iw=0; iw<storedWakes; iw++) {
    if (strcmp(storedWake[iw].filename, wakeData->inputFile)==0)
      break;
  }
  if (iw==storedWakes) {
    if (!SDDS_InitializeInput(&SDDSin, wakeData->inputFile) || SDDS_ReadPage(&SDDSin)!=1 ||
        (wakeData->wakePoints=SDDS_RowCount(&SDDSin))<0 ||
        wakeData->wakePoints<2) {
      fprintf(stderr, "Error: TRWAKE file is unreadable, or has insufficient data.\n", wakeData->inputFile);
      exit(1);
    }
    if (SDDS_CheckColumn(&SDDSin, wakeData->tColumn, "s", SDDS_ANY_FLOATING_TYPE, 
                         stdout)!=SDDS_CHECK_OK) {
      fprintf(stderr, "Error: problem with time column for TRWAKE file %s---check existence, units, and type", 
              wakeData->inputFile);
      exit(1);
    }
    if (!(wakeData->t=SDDS_GetColumnInDoubles(&SDDSin, wakeData->tColumn))) {
      fprintf(stderr, "Error: unable to retrieve time data from TRWAKE file %s\n", wakeData->inputFile);
      exit(1);
    }

    wakeData->W[0] = wakeData->W[1] = NULL;
    if (wakeData->WxColumn) {
      if (SDDS_CheckColumn(&SDDSin, wakeData->WxColumn, "V/C/m", SDDS_ANY_FLOATING_TYPE, 
                           stdout)!=SDDS_CHECK_OK) {
        fprintf(stderr, "Error: problem with Wx wake column for TRWAKE file %s---check existence, units, and type", 
                wakeData->inputFile);
        exit(1);
      }
      if (!(wakeData->W[0]=SDDS_GetColumnInDoubles(&SDDSin, wakeData->WxColumn))) {
        fprintf(stderr, "Error: unable to retrieve Wx data from TRWAKE file %s\n", wakeData->inputFile);
        exit(1);
      }
    }
    if (wakeData->WyColumn) {
      if (SDDS_CheckColumn(&SDDSin, wakeData->WyColumn, "V/C/m", SDDS_ANY_FLOATING_TYPE, 
                           stdout)!=SDDS_CHECK_OK) {
        fprintf(stderr, "Error: problem with Wy wake column for TRWAKE file %s---check existence, units, and type", 
                wakeData->inputFile);
        exit(1);
      }
      if (!(wakeData->W[1]=SDDS_GetColumnInDoubles(&SDDSin, wakeData->WyColumn))) {
        fprintf(stderr, "Error: unable to retrieve Wy data from TRWAKE file %s\n", wakeData->inputFile);
        exit(1);
      }
    }

    SDDS_Terminate(&SDDSin);

    /* record in wake storage */
    if (!(storedWake=SDDS_Realloc(storedWake, sizeof(*storedWake)*(storedWakes+1))) || 
        !SDDS_CopyString(&storedWake[storedWakes].filename, wakeData->inputFile))
      SDDS_Bomb("Memory allocation failure (WAKE)");
    storedWake[storedWakes].t = wakeData->t;
    storedWake[storedWakes].Wx = wakeData->W[0];
    storedWake[storedWakes].Wy = wakeData->W[1];
    storedWake[storedWakes].points = wakeData->wakePoints;
    wakeData->isCopy = 0;
    storedWakes++;
  }
  else {
    /* point to an existing wake */
    wakeData->t = storedWake[iw].t;
    wakeData->W[0] = storedWake[iw].Wx;
    wakeData->W[1] = storedWake[iw].Wy;
    wakeData->wakePoints = storedWake[iw].points;
    wakeData->isCopy = 1;
  }
  
  if (!wakeData->W[0] && !wakeData->W[1])
    bomb("no valid wake data for TRWAKE element", NULL);
  
  find_min_max(&tmin, &tmax, wakeData->t, wakeData->wakePoints);
  if (tmin==tmax)
    bomb("no time span in TRWAKE data", NULL);
  if (tmin!=0)
    bomb("TRWAKE function does not start at t=0.\n", NULL);
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
  
