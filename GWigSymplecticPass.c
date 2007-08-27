/* GWigSymplecticPass.c for
   elegant
*/

/*
 *---------------------------------------------------------------------------
 * Modification Log:
 * -----------------
 * .03  2004-11-01     M. Borland, APS, borland@aps.anl.gov
 *                              Converted from AT to elegant.
 *
 * .02  2003-06-18     J. Li, jing@fel.duke.edu
 *				Cleanup the code
 *
 * .01  2003-04-20     YK Wu, wu@fel.duke.edu
 *				GWiggler interface
 *
 *---------------------------------------------------------------------------
 *  Accelerator Physics Group, Duke FEL Lab, www.fel.duke.edu  
 */

#include <stdlib.h>
#include <math.h>
#include "gwig.h"

/******************************************************************************/
/* PHYSICS SECTION ************************************************************/

void GWigInit(struct gwig *Wig, 
	      double Ltot, /* total length of wiggler */
	      double Lw,   /* wiggler period (m) */
	      double Bmax, /* peak magnetic field (Tesla) */
	      int Nstep,   /* number of integration steps (per period?) */
	      int Nmeth,   /* integration method (2 or 4 for integration order) */
	      int NHharm,  /* number of horizontal harmonics (By) */
	      int NVharm,  /* number of vertical harmonics (Bx) */
	      double *pBy, /* data for horizontal harmonics (By) */ 
	                   /* (harmonic, strength, kx/kw, ky/kw, kz/kw, phase, ...) */
	      double *pBx, /* data for vertical harmonics (Bx) */
	      double *T1,  /* unused */
	      double *T2,  /* unused */
	      double *R1,  /* unused */
	      double *R2,  /* unused */
	      double pCentral, /* central momentum (beta*gamma) */
              long synchRad,   /* classical radiation ? */
              long isr         /* quantum (incoherent) radiation ? */
	      )
{
  double *tmppr;
  int    i;
  double kw;

  Wig->Po = pCentral;
  Wig->E0 = pCentral*XMC2;
  Wig->Pmethod = Nmeth;
  Wig->PN = Nstep;
  Wig->Nw = (int)(Ltot / Lw + 0.5);  /* In case Lw is slightly inaccurate */
  Wig->NHharm = NHharm;
  Wig->NVharm = NVharm;
  Wig->PB0 = Bmax;
  Wig->Lw  = Lw;
  Wig->srCoef = 0;

  if (Wig->sr = synchRad)
    Wig->srCoef = sqr(e_mks)*ipow(pCentral, 3)/(6*PI*epsilon_o*me_mks*sqr(c_mks));

  Wig->isrCoef = 0;
  if (Wig->isr = isr)
    Wig->isrCoef = re_mks*sqrt(55/(24.0*sqrt(3))*ipow(pCentral, 5)*137.0359895);

  kw = 2.0e0*PI/(Wig->Lw);
  Wig->Zw = 0.0;
  Wig->Aw = 0.0;
  tmppr = pBy;
  for (i = 0; i < NHharm; i++){
    tmppr++;
    Wig->HCw[i] = 0.0;
    Wig->HCw_raw[i] = *tmppr;

    tmppr++;
    Wig->Hkx[i]     = (*tmppr) * kw;

    tmppr++;
    Wig->Hky[i]     = (*tmppr) * kw;

    tmppr++;
    Wig->Hkz[i]     = (*tmppr) * kw;

    tmppr++;
    Wig->Htz[i]     =  *tmppr;

    tmppr++;
  }

  tmppr = pBx;
  for (i = 0; i < NVharm; i++){
    tmppr++;
    Wig->VCw[i] = 0.0;
    Wig->VCw_raw[i] = *tmppr;

    tmppr++;
    Wig->Vkx[i]     = (*tmppr) * kw;

    tmppr++;
    Wig->Vky[i]     = (*tmppr) * kw;

    tmppr++;
    Wig->Vkz[i]     = (*tmppr) * kw;

    tmppr++;
    Wig->Vtz[i]     =  *tmppr;

    tmppr++;
  }
  
  for (i = NHharm ; i< WHmax; i++) {
    Wig->HCw[i] = 0.0;
    Wig->HCw_raw[i] = 0.0;
    Wig->Hkx[i] = 0.0;
    Wig->Hky[i] = 0.0;
    Wig->Hkz[i] = 0.0;
    Wig->Htz[i] = 0.0;
  }
  for (i = NVharm ; i< WHmax; i++) {
    Wig->VCw[i] = 0.0;
    Wig->VCw_raw[i] = 0.0;
    Wig->Vkx[i] = 0.0;
    Wig->Vky[i] = 0.0;
    Wig->Vkz[i] = 0.0;
    Wig->Vtz[i] = 0.0;
  }
}

#define second 2
#define fourth 4

void GWigSymplecticPass(double **coord, long num_particles, double pCentral,
			CWIGGLER *cwiggler)
{	

  int c;
  double r6[6], denom;
  struct gwig Wig;
  MALIGN malign;

  if (!cwiggler->initialized) 
    InitializeCWiggler(cwiggler);

/*
  if ((cwiggler->sr || cwiggler->isr) && cwiggler->integrationOrder==fourth) {
    printf("Error: Can't presently include synchrotron radiation effects for fourth-order integration of CWIGGLER\n");
    exit(1);
  }
*/

  GWigInit(&Wig, cwiggler->length, cwiggler->length/cwiggler->periods, 
	   cwiggler->BMax, cwiggler->stepsPerPeriod, 
	   cwiggler->integrationOrder,
	   cwiggler->ByHarmonics, 
	   cwiggler->BxHarmonics,
	   cwiggler->ByData,
	   cwiggler->BxData,
	   NULL, NULL, NULL, NULL, pCentral,
           cwiggler->sr, cwiggler->isr);

  if (cwiggler->tilt)
    rotateBeamCoordinates(coord, num_particles, cwiggler->tilt);
  if (cwiggler->dx || cwiggler->dy || cwiggler->dz) {
    memset(&malign, 0, sizeof(malign));
    malign.dx = -cwiggler->dx;
    malign.dy = -cwiggler->dy;
    malign.dz = cwiggler->dz;
    offset_beam(coord, num_particles, &malign, pCentral);
  }

  for (c=0; c<num_particles; c++) {	
    /* convert from (x, x', y, y', s, delta) to Canonical coordinates 
     * (x, qx, y, qy, delta, s) 
     * d =  sqrt(1+sqr(xp)+sqr(yp))
     * qx = (1+delta)*xp/d
     * qy = (1+delta)*yp/d
     */
    r6[0] = coord[c][0];
    r6[2] = coord[c][2];
    r6[1] = (1+coord[c][5])*coord[c][1]/(denom=sqrt(1+sqr(coord[c][1])+sqr(coord[c][3])));
    r6[3] = (1+coord[c][5])*coord[c][3]/denom;
    /* For some reason, they swap the order here */
    r6[4] = coord[c][5]; 
    r6[5] = coord[c][4];

    /* Track through the wiggler */
    switch (cwiggler->integrationOrder) {
      case second :
	GWigPass_2nd(&Wig, r6);
	break;
      case fourth:
	GWigPass_4th(&Wig, r6);
	break;
      default:
	printf("Error: Invalid method integration order for CWIGGLER (use 2 or 4)\n");
	exit(1);
	break;
    }

    /* convert back to elegant coordinates */
    coord[c][0] = r6[0];
    coord[c][2] = r6[2];
    coord[c][5] = r6[4]; 
    coord[c][4] = r6[5];
    /* d = sqrt(sqr(1+delta)-sqr(qx)-sqr(qy))
     * xp = qx/d, yp=qy/d
     */
    denom = sqrt(sqr(1+coord[c][5])-sqr(r6[1])-sqr(r6[3]));
    coord[c][1] = r6[1]/denom;
    coord[c][3] = r6[3]/denom;
  }

  if (cwiggler->dx || cwiggler->dy || cwiggler->dz) {
    memset(&malign, 0, sizeof(malign));
    malign.dx = cwiggler->dx;
    malign.dy = cwiggler->dy;
    malign.dz = -cwiggler->dz;
    offset_beam(coord, num_particles, &malign, pCentral);
  }
  if (cwiggler->tilt)
    rotateBeamCoordinates(coord, num_particles, -cwiggler->tilt);
}

void InitializeCWiggler(CWIGGLER *cwiggler)
{
  double sumCmn2[2] = {0,0};
  long i;
  if (cwiggler->initialized)
    return;
  if (cwiggler->sinusoidal) {
    if (cwiggler->BxFile || cwiggler->ByFile)
      printf("*** Warning: CWIGGLER element has SINUSOIDAL=1, but also has filenames\n");
    cwiggler->BxHarmonics = 0;
    cwiggler->BxData = NULL;
    cwiggler->ByHarmonics = 1;
    cwiggler->ByData = tmalloc(sizeof(*(cwiggler->ByData))*6);
    cwiggler->ByData[0] = 0;  /* row */
    cwiggler->ByData[1] = 1;  /* Cmn */
    cwiggler->ByData[2] = 0;  /* kx */
    cwiggler->ByData[3] = 1;  /* ky */
    cwiggler->ByData[4] = 1;  /* kz */
    cwiggler->ByData[5] = 0;  /* phase */
  } else {
    ReadCWigglerHarmonics(&cwiggler->ByData, &cwiggler->ByHarmonics, 
                          cwiggler->ByFile, "By");
    ReadCWigglerHarmonics(&cwiggler->BxData, &cwiggler->BxHarmonics, 
                          cwiggler->BxFile, "Bx");
  }
  for (i=0; i<cwiggler->ByHarmonics; i++)
    sumCmn2[0] += sqr(cwiggler->ByData[6*i+1]);
  for (i=0; i<cwiggler->BxHarmonics; i++)
    sumCmn2[1] += sqr(cwiggler->BxData[6*i+1]);
  cwiggler->sumCmn2 = MAX(sumCmn2[0], sumCmn2[1]);
  cwiggler->initialized = 1;
}


long ReadCWigglerHarmonics(double **BData, long *harmonics, char *file, char *name)
{
  SDDS_DATASET SDDSin;
  double *Cmn, *kx, *ky, *kz, *phase;
  long row, rows;

  if (!file) {
    *harmonics = 0;
    return 0;
  }
  if (!SDDS_InitializeInput(&SDDSin, file) || !SDDS_ReadPage(&SDDSin)) {
    printf("Error: problem initializing file %s\n", file);
    exit(1);
  }
  if (!(rows=SDDS_RowCount(&SDDSin))) {
    printf("Error: no rows in file %s\n", file);
    exit(1);
  }
  if (!(Cmn=SDDS_GetColumnInDoubles(&SDDSin, "Cmn")) ||
      !(kx=SDDS_GetColumnInDoubles(&SDDSin, "KxOverKw")) ||
      !(ky=SDDS_GetColumnInDoubles(&SDDSin, "KyOverKw")) ||
      !(kz=SDDS_GetColumnInDoubles(&SDDSin, "KzOverKw")) ||
      !(phase=SDDS_GetColumnInDoubles(&SDDSin, "Phase"))) {
    printf("Error: problem reading file %s\n", file);
    printf("Check for existence of Cmn, KxOverKw, KyOverKw, KzOverKw, and Phase\n");
    exit(1);
  }
  
  *harmonics = rows;
  if (!(*BData = calloc(rows*6, sizeof(**BData)))) {
    printf("ReadCWigglerHarmonics: memory allocation failure (%ld harmonics)\n",
	   rows);
    exit(1);
  }

  for (row=0; row<rows; row++) {
    (*BData)[row*6]   = row;
    (*BData)[row*6+1] = Cmn[row];
    (*BData)[row*6+2] = kx[row];
    (*BData)[row*6+3] = ky[row];
    (*BData)[row*6+4] = kz[row];
    (*BData)[row*6+5] = phase[row];
  }
  if (!SDDS_Terminate(&SDDSin)) {
    printf("Warning: problem terminating CWIGGLER input file\n");
  }
  return rows;
}
