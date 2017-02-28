/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* routine: multipole()
 * purpose: apply kicks due to an arbitrary multipole
 * 
 * Michael Borland, 1991. 
 */
#include "mdb.h"
#include "track.h"
#include "multipole.h"
#ifdef HAVE_GPU
#include "gpu_multipole.h"
#endif

unsigned long multipoleKicksDone = 0;

#define ODD(j) ((j)%2)

typedef struct {
  char *filename;
  MULTIPOLE_DATA data;
  short steering;
} STORED_MULTIPOLE_DATA;
static STORED_MULTIPOLE_DATA *storedMultipoleData = NULL;
static long nMultipoleDataSets = 0;

long searchForStoredMultipoleData(char *multFile)
{
  long i;
  /* should use a hash table ! */
  for (i=0; i<nMultipoleDataSets; i++)
    if (strcmp(multFile, storedMultipoleData[i].filename)==0) {
      return i;
    }
  return -1;
}

void copyMultipoleDataset(MULTIPOLE_DATA *multData, long index)
{
  memcpy(multData, &storedMultipoleData[index].data, sizeof(*multData));
  multData->copy = 1;
}
 
void addMultipoleDatasetToStore(MULTIPOLE_DATA *multData, char *filename)
{
  printf("Adding file %s to multipole data store\n", filename, multData->filename);
  fflush(stdout);
  storedMultipoleData = SDDS_Realloc(storedMultipoleData, sizeof(*storedMultipoleData)*(nMultipoleDataSets+1));
  memcpy(&storedMultipoleData[nMultipoleDataSets].data, multData, sizeof(*multData));
  storedMultipoleData[nMultipoleDataSets].filename = tmalloc(sizeof(*filename)*(strlen(filename)+1));
  storedMultipoleData[nMultipoleDataSets].data.copy = 1;
  strcpy(storedMultipoleData[nMultipoleDataSets].filename, filename);
  nMultipoleDataSets++;
}

void readErrorMultipoleData(MULTIPOLE_DATA *multData,
                               char *multFile, long steering)
{
  SDDS_DATASET SDDSin;
  char buffer[1024];
  short anCheck, bnCheck, normalCheck, skewCheck;
  long index;

  if (!multFile || !strlen(multFile)) {
    multData->orders = 0;
    multData->initialized = 0;
    return;
  }
  if (multData->initialized)
    return;
  if ((index=searchForStoredMultipoleData(multFile))>=0) {
    copyMultipoleDataset(multData, index);
    return;
  }
  cp_str(&(multData->filename), multFile);
  if (!SDDS_InitializeInputFromSearchPath(&SDDSin, multFile)) {
    printf("Problem opening file %s\n", multFile);
    fflush(stdout);
    exitElegant(1);
  }
  if (SDDS_CheckColumn(&SDDSin, "order", NULL, SDDS_ANY_INTEGER_TYPE, NULL)!=SDDS_CHECK_OK ||
      SDDS_CheckParameter(&SDDSin, "referenceRadius", "m", SDDS_ANY_FLOATING_TYPE, NULL)!=SDDS_CHECK_OK) {
    printf("Problems with data in multipole file %s\n", multFile);
    fflush(stdout);
    exitElegant(1);
  }
  anCheck = SDDS_CheckColumn(&SDDSin, "an", NULL, SDDS_ANY_FLOATING_TYPE, NULL) == SDDS_CHECK_OK;
  normalCheck = SDDS_CheckColumn(&SDDSin, "normal", NULL, SDDS_ANY_FLOATING_TYPE, NULL) == SDDS_CHECK_OK;
  if (!anCheck && !normalCheck) {
    printf("Problems with data in multipole file %s: neither \"an\" nor \"normal\" column found\n", multFile);
    exitElegant(1);
  }
  if (anCheck && normalCheck) {
    printf("*** Warning: multipole file %s has both \"an\" and \"normal\" columns. \"normal\" used.\n", multFile);
    anCheck = 0;
  }

  bnCheck = SDDS_CheckColumn(&SDDSin, "bn", NULL, SDDS_ANY_FLOATING_TYPE, NULL) == SDDS_CHECK_OK;
  skewCheck = SDDS_CheckColumn(&SDDSin, "skew", NULL, SDDS_ANY_FLOATING_TYPE, NULL) == SDDS_CHECK_OK;
  if (!steering) {
    if (!bnCheck && !skewCheck) {
      printf("Problems with data in multipole file %s: neither \"bn\" nor \"skew\" column found\n", multFile);
      exitElegant(1);
    }
    if (bnCheck && skewCheck) {
      printf("*** Warning: multipole file %s has both \"bn\" and \"skew\" columns. \"skew\" used.\n", multFile);
      bnCheck = 0;
    }
  } else {
    if (bnCheck || skewCheck) {
      printf("*** Warning: Steering multipole file %s should not have bn or skew columns.\n",
            multFile);
      printf("Use \"normal\" column to specify multipole content for a horizontal steerer.\n");
      printf("Multipole content for vertical steerer is deduced from this.\n");
      fflush(stdout);
    }
  }
  
  if (SDDS_ReadPage(&SDDSin)!=1)  {
    sprintf(buffer, "Problem reading multipole file %s\n", multFile);
    SDDS_SetError(buffer);
    SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors);
    exitElegant(1);
  }
  if ((multData->orders = SDDS_RowCount(&SDDSin))<=0) {
    printf("Warning: no data in multipole file %s\n", multFile);
    fflush(stdout);
    SDDS_Terminate(&SDDSin);
  }
  if (!SDDS_GetParameterAsDouble(&SDDSin, "referenceRadius", &multData->referenceRadius) ||
      !(multData->order=SDDS_GetColumnInLong(&SDDSin, "order")) ||
      !(multData->an=SDDS_GetColumnInDoubles(&SDDSin, anCheck?"an":"normal")) || 
      (!steering && !(multData->bn=SDDS_GetColumnInDoubles(&SDDSin, bnCheck?"bn":"skew")))) {
    sprintf(buffer, "Unable to read multipole data for file %s\n", multFile);
    SDDS_SetError(buffer);
    SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors);
    exitElegant(1);
  }    
  if (steering &&
      !(multData->bn=SDDS_Malloc(sizeof(*(multData->bn))*multData->orders))) {
    printf("Memory allocation failure (readErrorMultipoleData)\n");
    exitElegant(1);
  }
  if (SDDS_ReadPage(&SDDSin)==2) {
    printf("Warning: multipole file %s has multiple pages, which are ignored\n",
            multFile);
    fflush(stdout);
  }
  SDDS_Terminate(&SDDSin);
  if (steering) {
    long i, j;
    /* check for disallowed multipoles */
    for (i=0; i<multData->orders; i++) {
      if (ODD(multData->order[i])) {
        printf("Error: steering multipole file %s has disallowed odd orders.\n",
                multFile);
        exitElegant(1);
      }
    }
    /* normalize to n=0 if present */
    /* find i such that order[i] is 0 (dipole) */
    for (i=0; i<multData->orders; i++) {
      if (multData->order[i]==0)
        break;
    }
    if (multData->orders>1 && i!=multData->orders) {
      /* dipole present */
      if (!multData->an[i]) {
        printf("Steering multipole data in %s is invalid: an is zero for order=0\n",
                multFile);
        exitElegant(1);
      }
      /* normalize to dipole for normal and skew separately */
      for (j=0; j<multData->orders; j++)
        if (j!=i)
          multData->an[j] /= multData->an[i];
      /* remove the dipole data */
      for (j=i+1; j<multData->orders; j++) {
        multData->an[j-1] = multData->an[j];
        multData->order[j-1] = multData->order[j];
      }
      multData->orders -= 1;
    }
    for (i=0; i<multData->orders; i++)
      multData->bn[i] = multData->an[i]*ipow(-1.0, multData->order[i]/2);
#ifdef DEBUG
    printf("Steering multipole data: \n");
    for (i=0; i<multData->orders; i++)
      printf("%ld: %e %e\n", multData->order[i],
              multData->an[i], multData->bn[i]);
#endif
  }

  if (steering) {
    long i;
    if (!(multData->KnL = SDDS_Malloc(sizeof(*multData->KnL)*multData->orders)) ||
        !(multData->JnL = SDDS_Malloc(sizeof(*multData->JnL)*multData->orders)) )
      bombTracking("memory allocation failure (readErrorMultipoleData)");
    for (i=0; i<multData->orders; i++) {
      multData->KnL[i] = 
        -1*multData->an[i]*dfactorial(multData->order[i])/ipow(multData->referenceRadius, multData->order[i]);
      multData->JnL[i] = 
        -1*multData->bn[i]*dfactorial(multData->order[i])/ipow(multData->referenceRadius, multData->order[i]);
    }
  }


  multData->initialized = 1;
  multData->copy = 0;

  addMultipoleDatasetToStore(multData, multFile);
}

void fillPowerArray(double x, double *xpow, long order)
{
  long i;
  xpow[0] = 1;
  for (i=1; i<=order; i++) {
    xpow[i] = xpow[i-1]*x;
  }
}

void initialize_fmultipole(FMULT *multipole)
{
  SDDS_DATASET SDDSin;
  char buffer[1024];
  MULTIPOLE_DATA *multData;
  
  multData = &(multipole->multData);
  if (multData->initialized)
    return;
  if (!multipole->filename)
    bombTracking("FMULT element doesn't have filename");
  if (!SDDS_InitializeInputFromSearchPath(&SDDSin, multipole->filename)) {
    sprintf(buffer, "Problem opening file %s (FMULT)\n", multipole->filename);
    SDDS_SetError(buffer);
    SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors);
    exitElegant(1);
  }
  if (SDDS_CheckColumn(&SDDSin, "order", NULL, SDDS_ANY_INTEGER_TYPE, stdout)!=SDDS_CHECK_OK ||
      SDDS_CheckColumn(&SDDSin, "KnL", NULL, SDDS_ANY_FLOATING_TYPE, stdout)!=SDDS_CHECK_OK ||
      SDDS_CheckColumn(&SDDSin, "JnL", NULL, SDDS_ANY_FLOATING_TYPE, stdout)!=SDDS_CHECK_OK)
    bombTracking("problems with data in FMULT input file");
  if (SDDS_ReadPage(&SDDSin)!=1)  {
    sprintf(buffer, "Problem reading FMULT file %s\n", multipole->filename);
    SDDS_SetError(buffer);
    SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors);
    exitElegant(1);
  }
  if ((multData->orders = SDDS_RowCount(&SDDSin))<=0) {
    printf("Warning: no data in FMULT file %s\n", multipole->filename);
    fflush(stdout);
    SDDS_Terminate(&SDDSin);
    return;
  }
  multData->JnL = NULL;
  if (!(multData->order=SDDS_GetColumnInLong(&SDDSin, "order")) ||
      !(multData->KnL=SDDS_GetColumnInDoubles(&SDDSin, "KnL")) ||
      !(multData->JnL=SDDS_GetColumnInDoubles(&SDDSin, "JnL"))) {
    sprintf(buffer, "Unable to read data for FMULT file %s\n", multipole->filename);
    SDDS_SetError(buffer);
    SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors);
    exitElegant(1);
  }    
  if (SDDS_ReadPage(&SDDSin)==2) {
    printf("Warning: FMULT file %s has multiple pages, which are ignored\n",
            multipole->filename);
    fflush(stdout);
  }
  SDDS_Terminate(&SDDSin);
  multData->initialized = 1;
}


long fmultipole_tracking(
                         double **particle,  /* initial/final phase-space coordinates */
                         long n_part,        /* number of particles */
                         FMULT *multipole,   /* multipole structure */
                         double p_error,     /* p_nominal/p_central */
                         double Po,
                         double **accepted,
                         double z_start
                         )
{
  double dummy;
  double dzLoss;
  long n_kicks;       /* number of kicks to split multipole into */
  long i_part, i_top, is_lost=0, i_order;
  double *coord;
  double drift;
  double x=0.0, xp=0.0, y=0.0, yp=0.0;
  double rad_coef;
  MULTIPOLE_DATA multData;
  
  if (!particle)
    bombTracking("particle array is null (fmultipole_tracking)");

  if (!multipole)
    bombTracking("null MULT pointer (fmultipole_tracking)");
  
  if (!multipole->multData.initialized)
    initialize_fmultipole(multipole);
  
  if ((n_kicks=multipole->n_kicks)<=0)
    bombTracking("n_kicks<=0 in fmultipole_tracking()");

  drift = multipole->length;

  if (multipole->synch_rad)
    rad_coef = sqr(particleCharge)*pow3(Po)/(6*PI*epsilon_o*sqr(c_mks)*particleMass);
  else
    rad_coef = 0;

  multData = multipole->multData;
  for (i_order=0; i_order<multData.orders; i_order++) {
    multData.KnL[i_order] *= (1+multipole->fse);
    if (multData.JnL)
      multData.JnL[i_order] *= (1+multipole->fse);
  }
  
  if (multipole->dx || multipole->dy || multipole->dz)
    offsetBeamCoordinates(particle, n_part, multipole->dx, multipole->dy, multipole->dz);
  if (multipole->tilt)
    rotateBeamCoordinates(particle, n_part, multipole->tilt);

  i_top = n_part-1;
  multipoleKicksDone += (i_top+1)*multData.orders*n_kicks*4;
  for (i_part=0; i_part<=i_top; i_part++) {
    if (!(coord = particle[i_part])) {
      printf("null coordinate pointer for particle %ld (fmultipole_tracking)", i_part);
      fflush(stdout);
      abort();
    }
    if (accepted && !accepted[i_part]) {
      printf("null accepted coordinates pointer for particle %ld (fmultipole_tracking)", i_part);
      fflush(stdout);
      abort();
    }

    is_lost = 0;
    if (!integrate_kick_multipole_ord4(coord, multipole->dx, multipole->dy, 0.0, 0.0, Po, rad_coef, 0.0,
                                       1, multipole->sqrtOrder, 0.0, n_kicks, drift, &multData, NULL, NULL, NULL,
                                       &dzLoss, NULL, 0))
      is_lost = 1;
    
    x = coord[0];
    y = coord[2];
    xp = coord[1];
    yp = coord[3];
    
    if (is_lost || isnan(x) || isnan(y) || isnan(xp) || isnan(yp) ||
	FABS(x)>COORD_LIMIT || FABS(y)>COORD_LIMIT ||
        FABS(xp)>SLOPE_LIMIT || FABS(yp)>SLOPE_LIMIT) {
      swapParticles(particle[i_part], particle[i_top]);
      if (accepted)
        swapParticles(accepted[i_part], accepted[i_top]);
      particle[i_top][4] = z_start + dzLoss;
      particle[i_top][5] = Po*(1+particle[i_top][5]);
      i_top--;
      i_part--;
    }
  }

  if (multipole->tilt)
    rotateBeamCoordinates(particle, n_part, -multipole->tilt);
  if (multipole->dx || multipole->dy || multipole->dz)
    offsetBeamCoordinates(particle, n_part, -multipole->dx, -multipole->dy, -multipole->dz);

  log_exit("fmultipole_tracking");
  return(i_top+1);
}


long multipole_tracking(
    double **particle,  /* initial/final phase-space coordinates */
    long n_part,        /* number of particles */
    MULT *multipole,    /* multipole structure */
    double p_error,     /* p_nominal/p_central */
    double Po,
    double **accepted,
    double z_start
    )
{
    double KnL;         /* integrated strength = L/(B.rho)*(Dx^n(By))_o for central momentum */
    long order;         /* order (n) */
    long n_kicks;       /* number of kicks to split multipole into */
    long i_part, i_kick, i, i_top, is_lost;
    double sum_Fx, sum_Fy, denom, qx, qy;
    double *coord;
    double drift;
    double *coef;
    double x, xp, y, yp, s, dp;
    double ratio, rad_coef;
    double beta0, beta1, p;
    static long maxOrder = -1;
    static double *xpow = NULL, *ypow = NULL;
    
    log_entry("multipole_tracking");

    if (!particle)
        bombTracking("particle array is null (multipole_tracking)");

    if (!multipole)
        bombTracking("null MULT pointer (multipole_tracking)");

    if ((n_kicks=multipole->n_kicks)<=0) 
      bombTracking("n_kicks<=0 in multipole_tracking()");

    if ((order=multipole->order)<0)
      bombTracking("order < 0 in multipole_tracking()");
    if (order>maxOrder || maxOrder==-1) {
      xpow = SDDS_Realloc(xpow, sizeof(*xpow)*(order+1));
      ypow = SDDS_Realloc(ypow, sizeof(*ypow)*(order+1));
      maxOrder = order;
    }

    if (!(coef = expansion_coefficients(order)))
      bombTracking("expansion_coefficients() returned null pointer (multipole_tracking)");

    drift = multipole->length/n_kicks/2;
    if (multipole->bore)
        /* KnL = d^nB/dx^n * L/(B.rho) = n! B(a)/a^n * L/(B.rho) */
        KnL = dfactorial(multipole->order)*multipole->BTipL/ipow(multipole->bore, multipole->order)*
              (particleCharge/(particleMass*c_mks*Po))*multipole->factor;
    else
      KnL = multipole->KnL*multipole->factor/n_kicks;

    if (KnL==0) {
      if (multipole->length==0)
        return n_part;
      exactDrift(particle, n_part, multipole->length);
      return n_part;
    }
    
    if (multipole->synch_rad)
        rad_coef = sqr(particleCharge)*pow3(Po)/(6*PI*epsilon_o*sqr(c_mks)*particleMass);
    else
        rad_coef = 0;

    if (multipole->dx || multipole->dy || multipole->dz)
      offsetBeamCoordinates(particle, n_part, multipole->dx, multipole->dy, multipole->dz);
    if (multipole->tilt)
      rotateBeamCoordinates(particle, n_part, multipole->tilt);

    i_top = n_part-1;
    multipoleKicksDone += (i_top+1)*n_kicks*4;
    for (i_part=0; i_part<=i_top; i_part++) {
        if (!(coord = particle[i_part])) {
            printf("null coordinate pointer for particle %ld (multipole_tracking)", i_part);
            fflush(stdout);
            abort();
            }
        if (accepted && !accepted[i_part]) {
            printf("null accepted coordinates pointer for particle %ld (multipole_tracking)", i_part);
            fflush(stdout);
            abort();
            }
        if (KnL==0) {
            coord[4] += multipole->length*sqrt(1+sqr(coord[1])+sqr(coord[3]));
            coord[0] += multipole->length*coord[1];
            coord[2] += multipole->length*coord[3];
            continue;
            }

        x = coord[0];
        xp = coord[1];
        y = coord[2];
        yp = coord[3];
        s  = 0;
        dp = coord[5];
        p = Po*(1+dp);
        beta0 = p/sqrt(sqr(p)+1);

#if defined(IEEE_MATH)
        if (isnan(x) || isnan(xp) || isnan(y) || isnan(yp)) {
            swapParticles(particle[i_part], particle[i_top]);
            if (accepted)
                swapParticles(accepted[i_part], accepted[i_top]);
            particle[i_top][4] = z_start;
            particle[i_top][5] = Po*(1+particle[i_top][5]);
            i_top--;
            i_part--;
            continue;
            }
#endif
        if (FABS(x)>COORD_LIMIT || FABS(y)>COORD_LIMIT ||
            FABS(xp)>SLOPE_LIMIT || FABS(yp)>SLOPE_LIMIT) {
            swapParticles(particle[i_part], particle[i_top]);
            if (accepted)
                swapParticles(accepted[i_part], accepted[i_top]);
            particle[i_top][4] = z_start;
            particle[i_top][5] = Po*(1+particle[i_top][5]);
            i_top--;
            i_part--;
            continue;
            }

        /* calculate initial canonical momenta */
        qx = (1+dp)*xp/(denom=sqrt(1+sqr(xp)+sqr(yp)));
        qy = (1+dp)*yp/denom;
        is_lost = 0;
        for (i_kick=0; i_kick<n_kicks; i_kick++) {
            if (drift) {
                x += xp*drift*(i_kick?2:1);
                y += yp*drift*(i_kick?2:1);
                s += (i_kick?2:1)*drift*sqrt(1+sqr(xp)+sqr(yp));
                }
            fillPowerArray(x, xpow, order);
            fillPowerArray(y, ypow, order);
            /* now sum up the terms for the multipole expansion */
            for (i=sum_Fx=sum_Fy=0; i<=order; i++) {
                if (ODD(i))
                    sum_Fx += coef[i]*xpow[order-i]*ypow[i];
                else
                    sum_Fy += coef[i]*xpow[order-i]*ypow[i];
                }
            /* apply kicks canonically */
            qx -= KnL*sum_Fy;
            qy += KnL*sum_Fx;
            if ((denom=sqr(1+dp)-sqr(qx)-sqr(qy))<=0) {
                is_lost = 1;
                break;
                }
            xp = qx/(denom=sqrt(denom));
            yp = qy/denom;
            if (rad_coef && drift) {
                qx /= (1+dp);
                qy /= (1+dp);
                dp -= rad_coef*sqr(KnL*(1+dp))*(sqr(sum_Fy)+sqr(sum_Fx))*sqrt(1+sqr(xp)+sqr(yp))/(2*drift);
                qx *= (1+dp);
                qy *= (1+dp);
                }
            }
        if (drift && !is_lost) {
            /* go through final drift */
            x += xp*drift;
            y += yp*drift;
            s += drift*sqrt(1 + sqr(xp) + sqr(yp));
            }

        coord[0] = x;
        coord[1] = xp;
        coord[2] = y;
        coord[3] = yp;
        coord[5] = dp;

        if (rad_coef) {
            p = Po*(1+dp);
            beta1 = p/sqrt(sqr(p)+1);
            coord[4] = beta1*(coord[4]/beta0 + 2*s/(beta0+beta1));
            }
        else 
            coord[4] += s;

#if defined(IEEE_MATH)
        if (isnan(x) || isnan(xp) || isnan(y) || isnan(yp)) {
            swapParticles(particle[i_part], particle[i_top]);
            if (accepted)
                swapParticles(accepted[i_part], accepted[i_top]);
            particle[i_top][4] = z_start;
            particle[i_top][5] = Po*(1+particle[i_top][5]);
            i_top--;
            i_part--;
            continue;
            }
#endif
        if (FABS(x)>COORD_LIMIT || FABS(y)>COORD_LIMIT ||
            FABS(xp)>SLOPE_LIMIT || FABS(yp)>SLOPE_LIMIT || is_lost) {
          swapParticles(particle[i_part], particle[i_top]);
          if (accepted)
            swapParticles(accepted[i_part], accepted[i_top]);
          particle[i_top][4] = z_start;
          particle[i_top][5] = Po*(1+particle[i_top][5]);
          i_top--;
          i_part--;
          continue;
        }

      }
    
    if (multipole->tilt)
      rotateBeamCoordinates(particle, n_part, -multipole->tilt);
    if (multipole->dx || multipole->dy || multipole->dz)
      offsetBeamCoordinates(particle, n_part, -multipole->dx, -multipole->dy, -multipole->dz);
    
    log_exit("multipole_tracking");
    return(i_top+1);
    }



double *expansion_coefficients(long n)
{
  static double **expansion_coef=NULL;
  static long *orderDone=NULL;
  static long maxOrder = -1;
  long i;
  
  if (n<=maxOrder && orderDone[n]) 
    return(expansion_coef[n]);

  if (n>maxOrder) {
    expansion_coef = trealloc(expansion_coef, sizeof(*expansion_coef)*(n+1));
    orderDone      = trealloc(orderDone, sizeof(*orderDone)*(n+1));
    for (i=maxOrder+1; i<=n; i++)
      orderDone[i] = 0;
    maxOrder = n;
  }

  expansion_coef[n] = tmalloc(sizeof(**expansion_coef)*(n+1));
  
  /* calculate expansion coefficients with signs for (x+iy)^n/n! */
  for (i=0; i<=n; i++) {
    expansion_coef[n][i] = (ODD(i/2)?-1.0:1.0)/(dfactorial(i)*dfactorial(n-i));
  }             
  orderDone[n] = 1;
  
  return(expansion_coef[n]);
}

long multipole_tracking2(
                         double **particle,   /* initial/final phase-space coordinates */
                         long n_part,         /* number of particles */
                         ELEMENT_LIST *elem,  /* element pointer */
                         double p_error,      /* p_nominal/p_central */
                         double Po,
                         double **accepted,
                         double z_start,
                         MAXAMP *maxamp,
                         /* from aperture_data command */
                         APERTURE_DATA *apFileData,
                         /* For return of accumulated change in sigmaDelta^2 */
                         double *sigmaDelta2
                         )
{
  double KnL;         /* integrated strength = L/(B.rho)*(Dx^n(By))_o for central momentum */
  double dx, dy, dz;  /* offsets of the multipole center */
  long order;         /* order (n) */
  long n_kicks, integ_order;
  long i_part, i_top, n_parts;
  double *coef, *coord;
  double drift;
  double tilt, rad_coef, isr_coef, xkick, ykick, dzLoss;
  KQUAD *kquad = NULL;
  KSEXT *ksext;
  KQUSE *kquse;
  KOCT *koct;
  static long sextWarning = 0, quadWarning = 0, octWarning = 0, quseWarning = 0;
  double lEffective = -1, lEnd = 0;
  short doEndDrift = 0;
  
  MULTIPOLE_DATA *multData = NULL, *steeringMultData = NULL, *edgeMultData = NULL;
  long sqrtOrder, freeMultData=0;
  MULT_APERTURE_DATA apertureData;
  double K2L;
  
#ifdef HAVE_GPU
  if(getElementOnGpu()){
    startGpuTimer();
    i_part = gpu_multipole_tracking2(n_part, elem, p_error, Po, accepted, z_start, maxamp, apFileData, sigmaDelta2);
#ifdef GPU_VERIFY     
    startCpuTimer();
    multipole_tracking2(particle, n_part, elem, p_error, Po, accepted, z_start, maxamp, apFileData, sigmaDelta2);
    compareGpuCpu(n_part, "multipole_tracking2");
#endif /* GPU_VERIFY */
    return i_part;
  }
#endif /* HAVE_GPU */


  log_entry("multipole_tracking2");

  if (!particle)
    bombTracking("particle array is null (multipole_tracking)");

  if (!elem)
    bombTracking("null element pointer (multipole_tracking2)");
  if (!elem->p_elem)
    bombTracking("null p_elem pointer (multipole_tracking2)");

  rad_coef = xkick = ykick = isr_coef = 0;
  sqrtOrder = 0;

  switch (elem->type) {
  case T_KQUAD:
    kquad = ((KQUAD*)elem->p_elem);
    n_kicks = kquad->n_kicks;
    order = 1;
    if ((lEffective = kquad->lEffective)<=0)
      lEffective = kquad->length;
    else {
      lEnd = (kquad->length-lEffective)/2;
      doEndDrift = 1;
    }
    if (kquad->bore)
      /* KnL = d^nB/dx^n * L/(B.rho) = n! B(a)/a^n * L/(B.rho) * (1+FSE) */
      KnL = kquad->B/kquad->bore*(particleCharge/(particleMass*c_mks*Po))*lEffective*(1+kquad->fse);
    else
      KnL = kquad->k1*lEffective*(1+kquad->fse);
    drift = lEffective;
    tilt = kquad->tilt;
    dx = kquad->dx;
    dy = kquad->dy;
    dz = kquad->dz;
    xkick = kquad->xkick*kquad->xKickCalibration;
    ykick = kquad->ykick*kquad->yKickCalibration;
    integ_order = kquad->integration_order;
    sqrtOrder = kquad->sqrtOrder?1:0;
    if (kquad->synch_rad)
      rad_coef = sqr(particleCharge)*pow3(Po)/(6*PI*epsilon_o*sqr(c_mks)*particleMass); 
    isr_coef = particleRadius*sqrt(55.0/(24*sqrt(3))*pow5(Po)*137.0359895);
    if (!kquad->isr || (kquad->isr1Particle==0 && n_part==1))
      /* Minus sign indicates we accumulate into sigmaDelta^2 only, don't perturb particles */
      isr_coef *= -1;
    if (kquad->length<1e-6 && (kquad->isr || kquad->synch_rad)) {
      rad_coef = isr_coef = 0;  /* avoid unphysical results */
      if (!quadWarning) {
        printf("**** Warning: one or more quadrupoles with length < 1e-6 have had SYNCH_RAD=0 and ISR=0 forced to avoid unphysical results.\n");
	quadWarning = 1;
      }
    }
    if (!kquad->multipolesInitialized) {
      /* read the data files for the error multipoles */
      readErrorMultipoleData(&(kquad->systematicMultipoleData),
                             kquad->systematic_multipoles, 0);
      readErrorMultipoleData(&(kquad->edgeMultipoleData),
                             kquad->edge_multipoles, 0);
      readErrorMultipoleData(&(kquad->randomMultipoleData),
                             kquad->random_multipoles, 0);
      readErrorMultipoleData(&(kquad->steeringMultipoleData), 
                             kquad->steering_multipoles, 1);
      kquad->multipolesInitialized = 1;
    }
    if (!kquad->totalMultipolesComputed) {
      computeTotalErrorMultipoleFields(&(kquad->totalMultipoleData),
                                       &(kquad->systematicMultipoleData),
                                       &(kquad->edgeMultipoleData),
                                       &(kquad->randomMultipoleData),
                                       &(kquad->steeringMultipoleData),
                                       KnL, 1);
      kquad->totalMultipolesComputed = 1;
    }
    multData = &(kquad->totalMultipoleData);
    edgeMultData = &(kquad->edgeMultipoleData);
    steeringMultData = &(kquad->steeringMultipoleData);
    break;
  case T_KSEXT:
    ksext = ((KSEXT*)elem->p_elem);
    n_kicks = ksext->n_kicks;
    order = 2;
    if (ksext->bore)
      /* KnL = d^nB/dx^n * L/(B.rho) = n! B(a)/a^n * L/(B.rho) * (1+FSE) */
      KnL = 2*ksext->B/sqr(ksext->bore)*(particleCharge/(particleMass*c_mks*Po))*ksext->length*(1+ksext->fse);
    else
      KnL = ksext->k2*ksext->length*(1+ksext->fse);
    drift = ksext->length;
    tilt = ksext->tilt;
    dx = ksext->dx;
    dy = ksext->dy;
    dz = ksext->dz;
    xkick = ksext->xkick*ksext->xKickCalibration;
    ykick = ksext->ykick*ksext->yKickCalibration;
    integ_order = ksext->integration_order;
    sqrtOrder = ksext->sqrtOrder?1:0;
    if (ksext->synch_rad)
      rad_coef = sqr(particleCharge)*pow3(Po)/(6*PI*epsilon_o*sqr(c_mks)*particleMass);
    isr_coef = particleRadius*sqrt(55.0/(24*sqrt(3))*pow5(Po)*137.0359895);
    if (!ksext->isr || (ksext->isr1Particle==0 && n_part==1))
      /* Minus sign indicates we accumulate into sigmaDelta^2 only, don't perturb particles */
      isr_coef *= -1;
    if (ksext->length<1e-6 && (ksext->isr || ksext->synch_rad)) {
      rad_coef = isr_coef = 0;  /* avoid unphysical results */
      if (!sextWarning) {
        printf("**** Warning: one or more sextupoles with length < 1e-6 have had SYNCH_RAD=0 and ISR=0 forced to avoid unphysical results.\n");
	sextWarning = 1;
      }
    }
    if (!ksext->multipolesInitialized) {
      /* read the data files for the error multipoles */
      readErrorMultipoleData(&(ksext->systematicMultipoleData),
                             ksext->systematic_multipoles, 0);
      readErrorMultipoleData(&(ksext->edgeMultipoleData),
                             ksext->edge_multipoles, 0);
      readErrorMultipoleData(&(ksext->randomMultipoleData),
                             ksext->random_multipoles, 0);
      readErrorMultipoleData(&(ksext->steeringMultipoleData), 
                             ksext->steering_multipoles, 1);
      ksext->multipolesInitialized = 1;
    }
    if (!ksext->totalMultipolesComputed) {
      computeTotalErrorMultipoleFields(&(ksext->totalMultipoleData),
                                       &(ksext->systematicMultipoleData),
                                       &(ksext->edgeMultipoleData),
                                       &(ksext->randomMultipoleData),
                                       &(ksext->steeringMultipoleData),
                                       KnL, 2);
      ksext->totalMultipolesComputed = 1;
    }
    multData = &(ksext->totalMultipoleData);
    edgeMultData = &(ksext->edgeMultipoleData);
    steeringMultData = &(ksext->steeringMultipoleData);
    break;
  case T_KOCT:
    koct = ((KOCT*)elem->p_elem);
    n_kicks = koct->n_kicks;
    order = 3;
    if (koct->bore)
      /* KnL = d^nB/dx^n * L/(B.rho) = n! B(a)/a^n * L/(B.rho) * (1+FSE) */
      KnL = 6*koct->B/ipow(koct->bore, 3)*(particleCharge/(particleMass*c_mks*Po))*koct->length*(1+koct->fse);
    else
      KnL = koct->k3*koct->length*(1+koct->fse);
    drift = koct->length;
    tilt = koct->tilt;
    dx = koct->dx;
    dy = koct->dy;
    dz = koct->dz;
    integ_order = koct->integration_order;
    sqrtOrder = koct->sqrtOrder?1:0;
    if (koct->synch_rad)
      rad_coef = sqr(particleCharge)*pow3(Po)/(6*PI*epsilon_o*sqr(c_mks)*particleMass);
    isr_coef = particleRadius*sqrt(55.0/(24*sqrt(3))*pow5(Po)*137.0359895);
    if (!koct->isr || (koct->isr1Particle==0 && n_part==1))
      /* Minus sign indicates we accumulate into sigmaDelta^2 only, don't perturb particles */
      isr_coef *= -1;
    if (koct->length<1e-6 && (koct->isr || koct->synch_rad)) {
      rad_coef = isr_coef = 0;  /* avoid unphysical results */
      if (!octWarning) {
        printf("**** Warning: one or more octupoles with length < 1e-6 have had SYNCH_RAD=0 and ISR=0 forced to avoid unphysical results.\n");
	octWarning = 1;
      }
    }
    if (!koct->multipolesInitialized) {
      /* read the data files for the error multipoles */
      readErrorMultipoleData(&(koct->systematicMultipoleData),
                             koct->systematic_multipoles, 0);
      readErrorMultipoleData(&(koct->randomMultipoleData),
                             koct->random_multipoles, 0);
      koct->multipolesInitialized = 1;
    }
    if (!koct->totalMultipolesComputed) {
      computeTotalErrorMultipoleFields(&(koct->totalMultipoleData),
                                       &(koct->systematicMultipoleData),
                                       NULL,
                                       &(koct->randomMultipoleData),
                                       NULL,
                                       KnL, 3);
      koct->totalMultipolesComputed = 1;
    }
    multData = &(koct->totalMultipoleData);
    break;
  case T_KQUSE:
    /* Implemented as a quadrupole with sextupole as a secondary multipole */
    kquse = ((KQUSE*)elem->p_elem);
    n_kicks = kquse->n_kicks;
    order = 1;
    KnL = kquse->k1*kquse->length*(1+kquse->fse1);
    drift = kquse->length;
    tilt = kquse->tilt;
    dx = kquse->dx;
    dy = kquse->dy;
    dz = kquse->dz;
    integ_order = kquse->integration_order;
    sqrtOrder = 0;
    if (kquse->synch_rad)
      rad_coef = sqr(particleCharge)*pow3(Po)/(6*PI*epsilon_o*sqr(c_mks)*particleMass); 
    isr_coef = particleRadius*sqrt(55.0/(24*sqrt(3))*pow5(Po)*137.0359895);
    if (!kquse->isr || (kquse->isr1Particle==0 && n_part==1))
      /* Minus sign indicates we accumulate into sigmaDelta^2 only, don't perturb particles */
      isr_coef *= -1;
    if (kquse->length<1e-6 && (kquse->isr || kquse->synch_rad)) {
      rad_coef = isr_coef = 0;  /* avoid unphysical results */
      if (!quseWarning) {
        printf("**** Warning: one or more KQUSE's with length < 1e-6 have had SYNCH_RAD=0 and ISR=0 forced to avoid unphysical results.\n");
	quseWarning = 1;
      }
    }
    K2L = kquse->k2*kquse->length*(1+kquse->fse2);
    if (K2L) {
      multData = tmalloc(sizeof(*multData));
      multData->orders = multData->initialized = 1;
      multData->randomized = 0;
      multData->order = tmalloc(sizeof(*(multData->order))*1);
      multData->order[0] = 2;
      multData->KnL = tmalloc(sizeof(*(multData->KnL))*1);
      multData->KnL[0] = K2L;
      multData->JnL = NULL;
      freeMultData = 1;
    }
    break;
  default:
    printf("error: multipole_tracking2() called for element %s--not supported!\n", elem->name);
    fflush(stdout);
    KnL = dx = dy = dz = tilt = drift = 0;
    integ_order = order = n_kicks = 0;
    exitElegant(1);
    break;
  }
  if (multData && !multData->initialized)
    multData = NULL;
  
  if (n_kicks<=0)
    bombTracking("n_kicks<=0 in multipole()");
  if (order<=0)
    bombTracking("order <= 0 in multipole()");
  if (integ_order!=2 && integ_order!=4) 
    bombTracking("multipole integration_order must be 2 or 4");
  
  if (!(coef = expansion_coefficients(order)))
    bombTracking("expansion_coefficients() returned null pointer (multipole_tracking)");

  i_top = n_part-1;
  if (integ_order==4) {
    if ((n_parts = ceil(n_kicks/4.0))<1)
      n_parts = 1;
    n_kicks = n_parts*4;
  }
  else
    n_parts = n_kicks;
  
  multipoleKicksDone += (i_top+1)*n_kicks;
  if (multData)
    multipoleKicksDone += (i_top+1)*n_kicks*multData->orders;

  setupMultApertureData(&apertureData, maxamp, tilt, apFileData, z_start+drift/2);
  
  if (dx || dy || dz)
    offsetBeamCoordinates(particle, n_part, dx, dy, dz);
  if (tilt)
    rotateBeamCoordinates(particle, n_part, tilt);

  if (doEndDrift) {
    exactDrift(particle, n_part, lEnd);
  }

  /* Fringe treatment, if any */
  switch (elem->type) {
  case T_KQUAD:
    if (kquad->edge1_effects>0)
      quadFringe(particle, n_part, kquad->k1, kquad->fringeIntM, kquad->fringeIntP, -1, kquad->edge1_effects-1,
                 kquad->edge1Linear, kquad->edge1NonlinearFactor);
    break;
  default:
    break;
  }
  
  if (sigmaDelta2)
    *sigmaDelta2 = 0;
  for (i_part=0; i_part<=i_top; i_part++) {
    if (!(coord = particle[i_part])) {
      printf("null coordinate pointer for particle %ld (multipole_tracking)", i_part);
      fflush(stdout);
      abort();
    }
    if (accepted && !accepted[i_part]) {
      printf("null accepted coordinates pointer for particle %ld (multipole_tracking)", i_part);
      fflush(stdout);
      abort();
    }

    if ((integ_order==4 &&
         !integrate_kick_multipole_ord4(coord, dx, dy, xkick, ykick,
                                        Po, rad_coef, isr_coef, order, sqrtOrder, KnL,
                                        n_parts, drift, 
                                        multData, edgeMultData, steeringMultData,
                                        &apertureData, &dzLoss, sigmaDelta2,
					elem->type==T_KQUAD?kquad->radial:0)) ||
        (integ_order==2 &&
         !integrate_kick_multipole_ord2(coord, dx, dy, xkick, ykick,
                                        Po, rad_coef, isr_coef, order, sqrtOrder, KnL, 
                                        n_parts, drift,
                                        multData, edgeMultData, steeringMultData,
                                        &apertureData, &dzLoss, sigmaDelta2,
					elem->type==T_KQUAD?kquad->radial:0))) {
      swapParticles(particle[i_part], particle[i_top]);
      if (accepted)
        swapParticles(accepted[i_part], accepted[i_top]);
      particle[i_top][4] = z_start+dzLoss;
      particle[i_top][5] = Po*(1+particle[i_top][5]);
      i_top--;
      i_part--;
      continue;
    }
  }
  if (sigmaDelta2)
    *sigmaDelta2 /= i_top+1;

  /* Fringe treatment, if any */
  switch (elem->type) {
  case T_KQUAD:
    if (kquad->edge2_effects>0)
      quadFringe(particle, n_part, kquad->k1, kquad->fringeIntM, kquad->fringeIntP, 1, kquad->edge2_effects-1,
                 kquad->edge2Linear, kquad->edge2NonlinearFactor);
    break;
  default:
    break;
  }

  if (doEndDrift) {
    exactDrift(particle, n_part, lEnd);
  }
  
  if (tilt)
    rotateBeamCoordinates(particle, n_part, -tilt);
  if (dx || dy || dz)
    offsetBeamCoordinates(particle, n_part, -dx, -dy, -dz);

  if (freeMultData && !multData->copy) {
    if (multData->order)
      free(multData->order);
    if (multData->KnL)
      free(multData->KnL);
    free(multData);
  }

  log_exit("multipole_tracking2");
  return(i_top+1);
}

int integrate_kick_multipole_ord2(double *coord, double dx, double dy, double xkick, double ykick,
                                  double Po, double rad_coef, double isr_coef,
                                  long order, long sqrtOrder, double KnL, long n_kicks, double drift,
                                  MULTIPOLE_DATA *multData, 
                                  MULTIPOLE_DATA *edgeMultData, 
                                  MULTIPOLE_DATA *steeringMultData,
                                  MULT_APERTURE_DATA *apData, double *dzLoss, double *sigmaDelta2,
				  long radial) 
{
  double p, qx, qy, denom, beta0, beta1, dp, s;
  double x, y, xp, yp, sum_Fx, sum_Fy;
  long i_kick, imult;
  
  drift = drift/n_kicks/2.0;
  KnL = KnL/n_kicks;
  
  x = coord[0];
  xp = coord[1];
  y = coord[2];
  yp = coord[3];
  s  = 0;
  dp = coord[5];
  p = Po*(1+dp);
  beta0 = p/sqrt(sqr(p)+1);

#if defined(IEEE_MATH)
  if (isnan(x) || isnan(xp) || isnan(y) || isnan(yp)) {
    return 0;
  }
#endif
  if (FABS(x)>COORD_LIMIT || FABS(y)>COORD_LIMIT ||
      FABS(xp)>SLOPE_LIMIT || FABS(yp)>SLOPE_LIMIT) {
    return 0;
  }

  /* apply steering corrector kick */
  xp += xkick/(1+dp)/2;
  yp += ykick/(1+dp)/2;

  /* calculate initial canonical momenta */
  denom = 1+sqr(xp)+sqr(yp);
  denom = EXSQRT(denom, sqrtOrder);
  qx = (1+dp)*xp/denom;
  qy = (1+dp)*yp/denom;

  if (steeringMultData && steeringMultData->orders) {
    /* apply steering corrector multipoles */
    for (imult=0; imult<steeringMultData->orders; imult++) {
      apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, x, y, 
                                      steeringMultData->order[imult], 
                                      steeringMultData->KnL[imult]*xkick/2, 0);
      apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, x, y, 
                                      steeringMultData->order[imult], 
                                      steeringMultData->JnL[imult]*ykick/2, 1);
    }
  }
  if (edgeMultData && edgeMultData->orders) {
    for (imult=0; imult<edgeMultData->orders; imult++) {
      apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, x, y, 
                                      edgeMultData->order[imult], 
                                      edgeMultData->KnL[imult], 0);
      apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, x, y, 
                                      edgeMultData->order[imult], 
                                      edgeMultData->JnL[imult], 1);
    }
  }
  /* We must do this in case steering or edge multipoles were run. We do it even if not in order
   * to avoid numerical precision issues that may subtly change the results
   */
  if ((denom=sqr(1+dp)-sqr(qx)-sqr(qy))<=0) {
    coord[0] = x;
    coord[2] = y;
    return 0;
  }
  denom = EXSQRT(denom, sqrtOrder);
  xp = qx/denom;
  yp = qy/denom;

  
  *dzLoss = 0;
  for (i_kick=0; i_kick<n_kicks; i_kick++) {
    if (drift) {
      x += xp*drift*(i_kick?2:1);
      y += yp*drift*(i_kick?2:1);
      s += drift*(i_kick?2:1)*sqrt(1 + sqr(xp) + sqr(yp));
      *dzLoss += drift*(i_kick?2:1);
    }
    if (apData && !checkMultAperture(x+dx, y+dy, apData))  {
      coord[0] = x;
      coord[2] = y;
      return 0;
    }

    if (!radial)
      apply_canonical_multipole_kicks(&qx, &qy, &sum_Fx, &sum_Fy, x, y, order, KnL, 0);
    else
      applyRadialCanonicalMultipoleKicks(&qx, &qy, &sum_Fx, &sum_Fy, x, y, order, KnL, 0);

    /* do kicks for spurious multipoles */
    if (multData) {
      for (imult=0; imult<multData->orders; imult++) {
        if (multData->KnL && multData->KnL[imult])
          apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, x, y, 
                                          multData->order[imult], 
                                          multData->KnL[imult]/n_kicks, 0);
        if (multData->JnL && multData->JnL[imult]) 
          apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, x, y, 
                                          multData->order[imult], 
                                          multData->JnL[imult]/n_kicks, 1);
      }
    }

    if ((denom=sqr(1+dp)-sqr(qx)-sqr(qy))<=0) {
      coord[0] = x;
      coord[2] = y;
      return 0;
    }
    denom = EXSQRT(denom, sqrtOrder);
    xp = qx/denom;
    yp = qy/denom;
    if ((rad_coef || isr_coef) && drift) {
      double deltaFactor, F2, dsFactor;
      deltaFactor = sqr(1+dp);
      F2 = (sqr(sum_Fy)+sqr(sum_Fx))*sqr(KnL/(2*drift));
      dsFactor = EXSQRT(1+sqr(xp)+sqr(yp), sqrtOrder)*2*drift;
      qx /= (1+dp);
      qy /= (1+dp);
      if (rad_coef)
	dp -= rad_coef*deltaFactor*F2*dsFactor;
      if (isr_coef>0)
	dp -= isr_coef*deltaFactor*pow(F2, 0.75)*sqrt(dsFactor)*gauss_rn_lim(0.0, 1.0, srGaussianLimit, random_2);
      if (sigmaDelta2)
        *sigmaDelta2 += sqr(isr_coef*deltaFactor)*pow(F2, 1.5)*dsFactor;
      qx *= (1+dp);
      qy *= (1+dp);
    }
  }
  if (drift) {
    /* go through final drift */
    x += xp*drift;
    y += yp*drift;
    s += drift*EXSQRT(1 + sqr(xp) + sqr(yp), sqrtOrder);
    *dzLoss += drift;
  }
  if (apData && !checkMultAperture(x+dx, y+dy, apData))  {
    coord[0] = x;
    coord[2] = y;
    return 0;
  }

  if (edgeMultData && edgeMultData->orders) {
    for (imult=0; imult<edgeMultData->orders; imult++) {
      apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, x, y, 
                                      edgeMultData->order[imult], 
                                      edgeMultData->KnL[imult], 0);
      apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, x, y, 
                                      edgeMultData->order[imult], 
                                      edgeMultData->JnL[imult], 1);
    }
  }
  if (steeringMultData && steeringMultData->orders) {
    /* apply steering corrector multipoles */
    for (imult=0; imult<steeringMultData->orders; imult++) {
      apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, x, y, 
                                      steeringMultData->order[imult], 
                                      steeringMultData->KnL[imult]*xkick/2, 0);
      apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, x, y, 
                                      steeringMultData->order[imult], 
                                      steeringMultData->JnL[imult]*ykick/2, 1);
    }
  }
  if ((denom=sqr(1+dp)-sqr(qx)-sqr(qy))<=0) {
    coord[0] = x;
    coord[2] = y;
    return 0;
  }
  denom = EXSQRT(denom, sqrtOrder);
  xp = qx/denom;
  yp = qy/denom;

  /* apply steering corrector kick */
  xp += xkick/(1+dp)/2;
  yp += ykick/(1+dp)/2;

  coord[0] = x;
  coord[1] = xp;
  coord[2] = y;
  coord[3] = yp;
  if (rad_coef || isr_coef) {
    p = Po*(1+dp);
    beta1 = p/sqrt(sqr(p)+1);
    coord[4] = beta1*(coord[4]/beta0 + 2*s/(beta0+beta1));
  }
  else 
    coord[4] += s;
  coord[5] = dp;

#if defined(IEEE_MATH)
  if (isnan(x) || isnan(xp) || isnan(y) || isnan(yp)) {
    return 0;
  }
#endif
  if (FABS(x)>COORD_LIMIT || FABS(y)>COORD_LIMIT ||
      FABS(xp)>SLOPE_LIMIT || FABS(yp)>SLOPE_LIMIT) {
    return 0;
  }
  return 1;
}


/* BETA is 2^(1/3) */
#define BETA 1.25992104989487316477

int integrate_kick_multipole_ord4(double *coord, double dx, double dy, double xkick, double ykick,
                                  double Po, double rad_coef, double isr_coef,
                                  long order, long sqrtOrder, double KnL, long n_parts, double drift,
                                  MULTIPOLE_DATA *multData, MULTIPOLE_DATA *edgeMultData, MULTIPOLE_DATA *steeringMultData,
                                  MULT_APERTURE_DATA *apData, double *dzLoss, double *sigmaDelta2,
				  long radial) 
{
  double p, qx, qy, denom, beta0, beta1, dp, s;
  double x, y, xp, yp, sum_Fx, sum_Fy;
  long i_kick, step, imult;
  double dsh;
  static double driftFrac[4] = {
    0.5/(2-BETA),  (1-BETA)/(2-BETA)/2,  (1-BETA)/(2-BETA)/2,  0.5/(2-BETA)
    } ;
  static double kickFrac[4] = {
    1./(2-BETA),  -BETA/(2-BETA),  1/(2-BETA),  0
    } ;
  
  drift = drift/n_parts;
  KnL = KnL/n_parts;

  x = coord[0];
  xp = coord[1];
  y = coord[2];
  yp = coord[3];
  s  = 0;
  dp = coord[5];
  p = Po*(1+dp);
  beta0 = p/sqrt(sqr(p)+1);

#if defined(IEEE_MATH)
  if (isnan(x) || isnan(xp) || isnan(y) || isnan(yp)) {
    return 0;
  }
#endif
  if (FABS(x)>COORD_LIMIT || FABS(y)>COORD_LIMIT ||
      FABS(xp)>SLOPE_LIMIT || FABS(yp)>SLOPE_LIMIT) {
    return 0;
  }

  /* apply steering corrector kick */
  xp += xkick/(1+dp)/2;
  yp += ykick/(1+dp)/2;

  /* calculate initial canonical momenta */
  qx = (1+dp)*xp/(denom=EXSQRT(1+sqr(xp)+sqr(yp), sqrtOrder));
  qy = (1+dp)*yp/denom;

  if (steeringMultData && steeringMultData->orders) {
    /* apply steering corrector multipoles */
    for (imult=0; imult<steeringMultData->orders; imult++) {
      if (steeringMultData->KnL[imult])
        apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, x, y, 
                                        steeringMultData->order[imult], 
                                        steeringMultData->KnL[imult]*xkick/2, 0);
      if (steeringMultData->JnL[imult]) 
        apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, x, y, 
                                        steeringMultData->order[imult], 
                                        steeringMultData->JnL[imult]*ykick/2, 1);
    }
  }
  if (edgeMultData && edgeMultData->orders) {
    for (imult=0; imult<edgeMultData->orders; imult++) {
      apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, x, y, 
                                      edgeMultData->order[imult], 
                                      edgeMultData->KnL[imult], 0);
      apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, x, y, 
                                      edgeMultData->order[imult], 
                                      edgeMultData->JnL[imult], 1);
    }
  }
  /* We must do this in case steering or edge multipoles were run. We do it even if not in order
   * to avoid numerical precision issues that may subtly change the results
   */
  if ((denom=sqr(1+dp)-sqr(qx)-sqr(qy))<=0) {
    coord[0] = x;
    coord[2] = y;
    return 0;
  }
  denom = EXSQRT(denom, sqrtOrder);
  xp = qx/denom;
  yp = qy/denom;

  *dzLoss = 0;
  for (i_kick=0; i_kick<n_parts; i_kick++) {
    if (apData && !checkMultAperture(x+dx, y+dy, apData))  {
      coord[0] = x;
      coord[2] = y;
      return 0;
    }
    for (step=0; step<4; step++) {
      if (drift) {
        dsh = drift*driftFrac[step];
        x += xp*dsh;
        y += yp*dsh;
        s += dsh*EXSQRT(1 + sqr(xp) + sqr(yp), sqrtOrder);
        *dzLoss += dsh;
      }

      if (!kickFrac[step])
        break;

      if (!radial)
	apply_canonical_multipole_kicks(&qx, &qy, &sum_Fx, &sum_Fy, x, y, 
					order, KnL*kickFrac[step], 0);
      else 
	applyRadialCanonicalMultipoleKicks(&qx, &qy, &sum_Fx, &sum_Fy, x, y, 
					   order, KnL*kickFrac[step], 0);
	
      if (multData) {
        /* do kicks for spurious multipoles */
        for (imult=0; imult<multData->orders; imult++) {
          if (multData->KnL && multData->KnL[imult]) {
            apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, x, y, 
                                            multData->order[imult], 
                                            multData->KnL[imult]*kickFrac[step]/n_parts,
                                            0);
          }
          if (multData->JnL && multData->JnL[imult]) {
            apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, x, y, 
                                            multData->order[imult], 
                                            multData->JnL[imult]*kickFrac[step]/n_parts,
                                            1);
          }
        }
      }
      if ((denom=sqr(1+dp)-sqr(qx)-sqr(qy))<=0) {
        coord[0] = x;
        coord[2] = y;
        return 0;
      }
      xp = qx/(denom=EXSQRT(denom, sqrtOrder));
      yp = qy/denom;
      if ((rad_coef || isr_coef) && drift) {
	double deltaFactor, F2, dsFactor, dsISRFactor;
        qx /= (1+dp);
        qy /= (1+dp);
	deltaFactor = sqr(1+dp);
	F2 = (sqr(sum_Fy)+sqr(sum_Fx))*sqr(KnL/drift);
	dsFactor = EXSQRT(1+sqr(xp)+sqr(yp), sqrtOrder);
	dsISRFactor = dsFactor*drift/3;   /* recall that kickFrac may be negative */
	dsFactor *= drift*kickFrac[step]; /* that's ok here, since we don't take sqrt */
	if (rad_coef)
	  dp -= rad_coef*deltaFactor*F2*dsFactor;
	if (isr_coef>0)
	  dp -= isr_coef*deltaFactor*pow(F2, 0.75)*sqrt(dsISRFactor)*gauss_rn_lim(0.0, 1.0, srGaussianLimit, random_2);
        if (sigmaDelta2)
          *sigmaDelta2 += sqr(isr_coef*deltaFactor)*pow(F2, 1.5)*dsFactor;
        qx *= (1+dp);
        qy *= (1+dp);
      }
    }
  }
  
  if (apData && !checkMultAperture(x+dx, y+dy, apData))  {
    coord[0] = x;
    coord[2] = y;
    return 0;
  }
  
  if (edgeMultData && edgeMultData->orders) {
    for (imult=0; imult<edgeMultData->orders; imult++) {
      apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, x, y, 
                                      edgeMultData->order[imult], 
                                      edgeMultData->KnL[imult], 0);
      apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, x, y, 
                                      edgeMultData->order[imult], 
                                      edgeMultData->JnL[imult], 1);
    }
  }
  if (steeringMultData && steeringMultData->orders) {
    /* apply steering corrector multipoles */
    for (imult=0; imult<steeringMultData->orders; imult++) {
      if (steeringMultData->KnL[imult]) 
        apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, x, y, 
                                        steeringMultData->order[imult], 
                                        steeringMultData->KnL[imult]*xkick/2, 0);
      if (steeringMultData->JnL[imult]) 
        apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, x, y, 
                                        steeringMultData->order[imult], 
                                        steeringMultData->JnL[imult]*ykick/2, 1);
    }
  }
  if ((denom=sqr(1+dp)-sqr(qx)-sqr(qy))<=0) {
    coord[0] = x;
    coord[2] = y;
    return 0;
  }
  denom = EXSQRT(denom, sqrtOrder);
  xp = qx/denom;
  yp = qy/denom;

  /* apply steering corrector kick */
  xp += xkick/(1+dp)/2;
  yp += ykick/(1+dp)/2;

  coord[0] = x;
  coord[1] = xp;
  coord[2] = y;
  coord[3] = yp;
  if (rad_coef) {
    p = Po*(1+dp);
    beta1 = p/sqrt(sqr(p)+1);
    coord[4] = beta1*(coord[4]/beta0 + 2*s/(beta0+beta1));
  }
  else 
    coord[4] += s;
  coord[5] = dp;

#if defined(IEEE_MATH)
  if (isnan(x) || isnan(xp) || isnan(y) || isnan(yp)) {
    return 0;
  }
#endif
  if (FABS(x)>COORD_LIMIT || FABS(y)>COORD_LIMIT ||
      FABS(xp)>SLOPE_LIMIT || FABS(yp)>SLOPE_LIMIT) {
    return 0;
  }
  return 1;
}

void apply_canonical_multipole_kicks(double *qx, double *qy, 
                                     double *sum_Fx_return, double *sum_Fy_return,
                                     double x, double y,
                                     long order, double KnL, long skew)
{
  long i;
  double sum_Fx, sum_Fy;
  double *coef;
  static long maxOrder = -1;
  static double *xpow = NULL, *ypow = NULL;
  if (order>maxOrder || maxOrder==-1) {
    xpow = SDDS_Realloc(xpow, sizeof(*xpow)*(order+1));
    ypow = SDDS_Realloc(ypow, sizeof(*ypow)*(order+1));
    maxOrder = order;
  }
  if (sum_Fx_return)
    *sum_Fx_return = 0;
  if (sum_Fy_return)
    *sum_Fy_return = 0;
  coef = expansion_coefficients(order);

  fillPowerArray(x, xpow, order);
  fillPowerArray(y, ypow, order);
  
  /* sum up the terms for the multipole expansion */
  for (i=sum_Fx=sum_Fy=0; i<=order; i++) {
    /*
    if (ODD(i))
      sum_Fx += coef[i]*ipow(x, order-i)*ipow(y, i);
    else
      sum_Fy += coef[i]*ipow(x, order-i)*ipow(y, i);
      */
    if (ODD(i))
      sum_Fx += coef[i]*xpow[order-i]*ypow[i];
    else
      sum_Fy += coef[i]*xpow[order-i]*ypow[i];
  }
  if (skew) {
    SWAP_DOUBLE(sum_Fx, sum_Fy);
    sum_Fx = -sum_Fx;
  }
  /* add the kicks */
  *qx -= KnL*sum_Fy;
  *qy += KnL*sum_Fx;
  if (sum_Fx_return)
    *sum_Fx_return = sum_Fx;
  if (sum_Fy_return)
    *sum_Fy_return = sum_Fy;
}

void applyRadialCanonicalMultipoleKicks(double *qx, double *qy, 
					double *sum_Fx_return, double *sum_Fy_return,
					double x, double y,
					long order, double KnL, long skew)
{
  long i;
  double sum_Fx, sum_Fy;
  double *coef;
  static long maxOrder = -1;
  static double *xpow = NULL, *ypow = NULL;

  if (order>maxOrder || maxOrder==-1) {
    xpow = SDDS_Realloc(xpow, sizeof(*xpow)*(order+1));
    ypow = SDDS_Realloc(ypow, sizeof(*ypow)*(order+1));
    maxOrder = order;
  }
  fillPowerArray(x, xpow, order);
  fillPowerArray(y, ypow, order);
  
  if (sum_Fx_return)
    *sum_Fx_return = 0;
  if (sum_Fy_return)
    *sum_Fy_return = 0;
  coef = expansion_coefficients(order);
  i = 0;
  /* now sum up the terms for the multipole expansion */
  for (sum_Fx=sum_Fy=0; i<=order; i++) {
    if (ODD(i))
      sum_Fx -= coef[i-1]*xpow[order-i]*ypow[i];
    else
      sum_Fy += coef[i]*xpow[order-i]*ypow[i];
  }
  if (skew) {
    SWAP_DOUBLE(sum_Fx, sum_Fy);
    sum_Fx = -sum_Fx;
  }
  /* add the kicks */
  *qx -= KnL*sum_Fy;
  *qy += KnL*sum_Fx;
  if (sum_Fx_return)
    *sum_Fx_return = sum_Fx;
  if (sum_Fy_return)
    *sum_Fy_return = sum_Fy;
}

void randomizeErrorMultipoleFields(MULTIPOLE_DATA *randomMult)
{
  long i;
  double nFactorial, rpow, rn1, rn2;
  if (!randomMult || randomMult->randomized)
    return;
  for (i=0; i<randomMult->orders; i++) {
    nFactorial = dfactorial(randomMult->order[i]);
    rpow = ipow(randomMult->referenceRadius, randomMult->order[i]);
    rn1 = gauss_rn_lim(0.0, 1.0, 2.0, random_1_elegant);
    rn2 = gauss_rn_lim(0.0, 1.0, 2.0, random_1_elegant);
    randomMult->anMod[i] = randomMult->an[i]*nFactorial*rn1/rpow;
    randomMult->bnMod[i] = randomMult->bn[i]*nFactorial*rn2/rpow;
  }
#ifdef DEBUG_RANDOMIZE
  printf("randomized multipoles\n");
#endif
  randomMult->randomized = 1;
}

void computeTotalErrorMultipoleFields(MULTIPOLE_DATA *totalMult,
                                      MULTIPOLE_DATA *systematicMult,
                                      MULTIPOLE_DATA *edgeMult,
                                      MULTIPOLE_DATA *randomMult,
                                      MULTIPOLE_DATA *steeringMult,
                                      double KmL, long rootOrder)
{
  long i;
  double sFactor=0.0, rFactor=0.0;
  
  if (!totalMult->initialized) {
    totalMult->initialized = 1;
    /* make a list of unique orders for random and systematic multipoles */
    if (systematicMult->orders && randomMult->orders &&
        systematicMult->orders!=randomMult->orders) {
      fprintf(stderr, "Issue with files %s and %s\n", 
	      systematicMult->filename, randomMult->filename);
      bombTracking("The number of systematic and random multipole error orders must be the same for any given element");
    }
    if (systematicMult->orders)
      totalMult->orders = systematicMult->orders;
    else
      totalMult->orders = randomMult->orders;
    if (totalMult->orders) {
      if (!(totalMult->order=SDDS_Malloc(sizeof(*totalMult->order)*(totalMult->orders))))
        bombTracking("memory allocation failure (computeTotalMultipoleFields)");
      if (systematicMult->orders &&
          (!(systematicMult->anMod=SDDS_Malloc(sizeof(*systematicMult->anMod)*systematicMult->orders)) ||
           !(systematicMult->bnMod=SDDS_Malloc(sizeof(*systematicMult->bnMod)*systematicMult->orders)) ||
           !(systematicMult->KnL=SDDS_Malloc(sizeof(*systematicMult->KnL)*systematicMult->orders)) ||
           !(systematicMult->JnL=SDDS_Malloc(sizeof(*systematicMult->JnL)*systematicMult->orders))))
        bombTracking("memory allocation failure (computeTotalMultipoleFields)");
      if (randomMult->orders &&
          (!(randomMult->anMod=SDDS_Malloc(sizeof(*randomMult->anMod)*randomMult->orders)) ||
           !(randomMult->bnMod=SDDS_Malloc(sizeof(*randomMult->bnMod)*randomMult->orders)) ||
           !(randomMult->KnL=SDDS_Malloc(sizeof(*randomMult->KnL)*randomMult->orders)) ||
           !(randomMult->JnL=SDDS_Malloc(sizeof(*randomMult->JnL)*randomMult->orders))))
        bombTracking("memory allocation failure (computeTotalMultipoleFields");
      if (!(totalMult->KnL = SDDS_Malloc(sizeof(*totalMult->KnL)*totalMult->orders)) ||
          !(totalMult->JnL = SDDS_Malloc(sizeof(*totalMult->JnL)*totalMult->orders)) )
        bombTracking("memory allocation failure (computeTotalMultipoleFields)");
      for (i=0; i<totalMult->orders; i++) {
        if (systematicMult->orders && randomMult->orders &&
            systematicMult->order[i]!=randomMult->order[i])
          bombTracking("multipole orders in systematic and random lists must match up for any given element.");
        if (systematicMult->orders) {
          totalMult->order[i] = systematicMult->order[i] ;
          systematicMult->anMod[i] = systematicMult->an[i]*dfactorial(systematicMult->order[i])/
            ipow(systematicMult->referenceRadius, systematicMult->order[i]);
          systematicMult->bnMod[i] = systematicMult->bn[i]*dfactorial(systematicMult->order[i])/
            ipow(systematicMult->referenceRadius, systematicMult->order[i]);
        } else {
          totalMult->order[i] = randomMult->order[i];
          /* anMod and bnMod will be computed later for randomized multipoles */
        }
      }
    }
    if (edgeMult && edgeMult->orders) {
      if (!(edgeMult->anMod=SDDS_Malloc(sizeof(*edgeMult->anMod)*edgeMult->orders)) ||
          !(edgeMult->bnMod=SDDS_Malloc(sizeof(*edgeMult->bnMod)*edgeMult->orders)) ||
          !(edgeMult->KnL=SDDS_Malloc(sizeof(*edgeMult->KnL)*edgeMult->orders)) ||
          !(edgeMult->JnL=SDDS_Malloc(sizeof(*edgeMult->JnL)*edgeMult->orders)))
        bombTracking("memory allocation failure (computeTotalMultipoleFields");
      for (i=0; i<edgeMult->orders; i++) {
        edgeMult->anMod[i] = edgeMult->an[i]*dfactorial(edgeMult->order[i])/
          ipow(edgeMult->referenceRadius, edgeMult->order[i]);
        edgeMult->bnMod[i] = edgeMult->bn[i]*dfactorial(edgeMult->order[i])/
          ipow(edgeMult->referenceRadius, edgeMult->order[i]);
      }
    }
  }
  
  if (randomMult->orders)
    randomizeErrorMultipoleFields(randomMult);
  
  /* body multipoles:
   * compute normal (KnL) and skew (JnL) from an and bn
   * KnL = an*n!/r^n*(KmL*r^m/m!), 
   * JnL = bn*n!/r^n*(KmL*r^m/m!), where m is the root order 
   * of the magnet with strength KmL
   */
  if (systematicMult->orders)
    sFactor = KmL/dfactorial(rootOrder)*ipow(systematicMult->referenceRadius, rootOrder);
  if (randomMult->orders)
    rFactor = KmL/dfactorial(rootOrder)*ipow(randomMult->referenceRadius, rootOrder);
  for (i=0; i<totalMult->orders; i++) {
    totalMult->KnL[i] = totalMult->JnL[i] = 0;
    if (systematicMult->orders) {
      totalMult->KnL[i] += sFactor*systematicMult->anMod[i];
      totalMult->JnL[i] += sFactor*systematicMult->bnMod[i];
    }
    if (randomMult->orders) {
      totalMult->KnL[i] += rFactor*randomMult->anMod[i];
      totalMult->JnL[i] += rFactor*randomMult->bnMod[i];
    }
  }

  if (edgeMult && edgeMult->orders) {
    /* edge multipoles: 
     * compute normal (KnL) and skew (JnL) from an and bn
     * KnL = an*n!/r^n*(KmL*r^m/m!), 
     * JnL = bn*n!/r^n*(KmL*r^m/m!), where m is the root order 
     * of the magnet with strength KmL
     */
    sFactor = KmL/dfactorial(rootOrder)*ipow(edgeMult->referenceRadius, rootOrder);
    for (i=0; i<edgeMult->orders; i++) {
      edgeMult->KnL[i] = sFactor*edgeMult->anMod[i];
      edgeMult->JnL[i] = sFactor*edgeMult->bnMod[i];
    }
  }
  
  if (steeringMult) {
    /* same for steering multipoles, but compute KnL/theta and JnL/theta (in this case m=0) */
    for (i=0; i<steeringMult->orders; i++) {
      steeringMult->KnL[i] = 
        -1*steeringMult->an[i]*dfactorial(steeringMult->order[i])/ipow(steeringMult->referenceRadius, steeringMult->order[i]);
      steeringMult->JnL[i] = 
        -1*steeringMult->bn[i]*dfactorial(steeringMult->order[i])/ipow(steeringMult->referenceRadius, steeringMult->order[i]);
    }
  }

}

void setupMultApertureData(MULT_APERTURE_DATA *apertureData, MAXAMP *maxamp, double tilt, 
     APERTURE_DATA *apFileData, double zPosition)
{
  double x_max, y_max;
  /* zPosition=_start+drift/2 */
  apertureData->xCen = apertureData->yCen = 0;
  x_max = y_max = apertureData->xMax = apertureData->yMax = 0;
  apertureData->elliptical = 0;
  apertureData->present = apertureData->openSide = 0;
  if (maxamp) {
    x_max = apertureData->xMax = maxamp->x_max;
    y_max = apertureData->yMax = maxamp->y_max;
    apertureData->elliptical = maxamp->elliptical;
    apertureData->present = x_max>0 || y_max>0;
    apertureData->xExponent = 
      apertureData->yExponent = maxamp->exponent;
    if (maxamp->yExponent)
      apertureData->yExponent = maxamp->yExponent;
    apertureData->openSide = determineOpenSideCode(maxamp->openSide);
  }
  if (apFileData && apFileData->initialized) {
    /* If there is file-based aperture data, it may override MAXAMP data. */
    double xCenF, yCenF, xMaxF, yMaxF;
    apertureData->present = 1;
    if (interpolateApertureData(zPosition,  apFileData, 
                                &xCenF, &yCenF, &xMaxF, &yMaxF)) {
      if (x_max<=0 || (x_max>fabs(xCenF+xMaxF) && x_max>fabs(xCenF-xMaxF))) {
        apertureData->xMax = xMaxF;
        apertureData->xCen = xCenF;
        apertureData->elliptical = 0;
        apertureData->openSide = 0;
      }
      if (y_max<=0 || (y_max>fabs(yCenF+yMaxF) && y_max>fabs(yCenF-yMaxF))) {
        apertureData->yMax = yMaxF;
        apertureData->yCen = yCenF;
        apertureData->elliptical = 0;
        apertureData->openSide = 0;
      }
    }
  }
  if (fabs(tilt)>0.1) {
    /* If rotation is greater than 100 mrad, disable aperture inside the element */
    /* Prevents unexpected results with skew elements */
    apertureData->present = 0;
  }
}

long checkMultAperture(double x, double y, MULT_APERTURE_DATA *apData) 
{
  double xa, yb;
  if (!apData || !apData->present)
    return 1;

  x -= apData->xCen;
  y -= apData->yCen;

  if (apData->elliptical==0 || apData->xMax<=0 || apData->yMax<=0) {
    /* rectangular or one-dimensional */
    if ((apData->xMax>0 && fabs(x)>apData->xMax) ||
        (apData->yMax>0 && fabs(y)>apData->yMax)) {
      if (apData->openSide==0 ||
          evaluateLostWithOpenSides(apData->openSide, x, y, apData->xMax, apData->yMax))
        return 0;
    }
    return 1;
  }

  /* Elliptical or super-elliptical */
  xa = x/apData->xMax;
  yb = y/apData->yMax;
  if ((ipow(xa, apData->xExponent) + ipow(yb, apData->yExponent))>=1) {
    if (apData->openSide==0 ||
        evaluateLostWithOpenSides(apData->openSide, x, y, apData->xMax, apData->yMax))
      return 0;
  }
  return 1;
}
