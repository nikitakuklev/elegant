/* Copyright 2016 by Michael Borland, Ryan Lindberg, and Argonne National Laboratory,
 * all rights reserved.
 */
#include "mdb.h"
#include "SDDS.h"
#include "track.h"


typedef struct {
  double radius;       /* reference radius */
  long nz;             /* number of z points */
  double dz;           /* z spacing */
  long nm;             /* number of values of m (angular harmonic) */
  long *m;             /* value of m */
  double nGradients;   /* number of gradient functions per m */
  /* See M. Venturini and A. Dragt, NIM A 427 (1999) 387-392, Eq. 9. */
  double ***Cmns;      /* generalized gradient: Cnms[im][in][iz] */
  double ***dCmns_dz;  /* z derivative of generalized gradient */
} STORED_BGGEXP_DATA;

static STORED_BGGEXP_DATA *storedBGGExpData = NULL;
static long nBGGExpDataSets = 0;
static htab *fileHashTable = NULL;

#define BUFSIZE 16834

long addBGGExpData(char *filename)
{
  SDDS_DATASET SDDSin;
  TRACKING_CONTEXT tcontext;
  char buffer[BUFSIZE];
  long im, ic, nc, readCode, nz;
  short m;
  
  if (!fileHashTable)
    fileHashTable = hcreate(12);

  if (hfind(fileHashTable, filename, strlen(filename))) {
    return *((long*)hstuff(fileHashTable));
  }

  getTrackingContext(&tcontext);
  printf("Adding BGGEXP data from file %s for element %s #%ld\n", filename, tcontext.elementName, tcontext.elementOccurrence);
  fflush(stdout);
  
  if (!SDDS_InitializeInputFromSearchPath(&SDDSin, filename))
    bombElegantVA("Unable to read file %s for BGGEXP %s #%ld\n", filename, tcontext.elementName, tcontext.elementOccurrence);

  if (SDDS_CheckParameter(&SDDSin, "m", NULL, SDDS_SHORT, stderr)!=SDDS_CHECK_OK) 
    bombElegantVA("Unable to find short-integer parameter \"m\" in file %s for BGGEXP %s #%ld\n", filename, tcontext.elementName, tcontext.elementOccurrence);
  /* 
  if (SDDS_CheckParameter(&SDDSin, "referenceRadius", "m", SDDS_ANY_FLOATING_TYPE, stderr)!=SDDS_CHECK_OK) 
    bombElegantVA("Unable to find floating point parameter \"referenceRadius\" with units \"m\" in file %s for BGGEXP %s #%ld\n", filename, tcontext.elementName, tcontext.elementOccurrence);
  */

  /* Check presence of z column */
  if (SDDS_CheckColumn(&SDDSin, "z", "m", SDDS_ANY_FLOATING_TYPE, stderr)!=SDDS_CHECK_OK) 
    bombElegantVA("Unable to find floating-point column \"z\" with units \"m\" in file %s for BGGEXP %s #%ld\n", filename, tcontext.elementName, tcontext.elementOccurrence);

  /* Check presence of Cnm* columns */
  ic = 0;
  while (1) {
    snprintf(buffer, BUFSIZE, "Cnm%ld", 2*ic);
    if (SDDS_CheckColumn(&SDDSin, buffer, NULL, SDDS_ANY_FLOATING_TYPE, NULL)!=SDDS_CHECK_OK)
      break;
    ic ++;
  }
  if (ic==0) 
    bombElegantVA("Unable to find any floating-point columns Cnm* in file %s for BGGEXP %s #%ld\n", filename, tcontext.elementName, tcontext.elementOccurrence);
  nc = ic;
  printf("Found %ld Cnm* columns\n", nc);

  /* Check for presence of matching dCnmXXX/dz columns */
  for (ic=0; ic<nc; ic++) {
    snprintf(buffer, BUFSIZE, "dCnm%ld/dz", 2*ic);
    if (SDDS_CheckColumn(&SDDSin, buffer, NULL, SDDS_ANY_FLOATING_TYPE, stderr)!=SDDS_CHECK_OK)
      break;
  }
  if (ic!=nc) 
    bombElegantVA("Unable to find matching floating-point columns dCnm*/dz in file %s for BGGEXP %s #%ld\n", filename, tcontext.elementName, tcontext.elementOccurrence);

  if (!(storedBGGExpData = SDDS_Realloc(storedBGGExpData, sizeof(*storedBGGExpData)*(nBGGExpDataSets+1))))
    bombElegantVA("Memory allocation error reading data from file %s for BGGEXP %s #%ld\n", filename, tcontext.elementName, tcontext.elementOccurrence);

  storedBGGExpData[nBGGExpDataSets].nm = storedBGGExpData[nBGGExpDataSets].nz = 0;
  storedBGGExpData[nBGGExpDataSets].nGradients = nc;
  storedBGGExpData[nBGGExpDataSets].m = NULL;
  storedBGGExpData[nBGGExpDataSets].Cmns = NULL;
  storedBGGExpData[nBGGExpDataSets].dCmns_dz = NULL;

  im = nz = 0;
  while ((readCode=SDDS_ReadPage(&SDDSin))>0) {
    if (!SDDS_GetParameter(&SDDSin, "m", &m) || m<=1)
      bombElegantVA("Problem with value of m for page %ld of file %s for BGGEXP %s #%ld\n", readCode, filename, tcontext.elementName, tcontext.elementOccurrence);
    if (readCode==1) {
      long iz;
      double dz0, dz, *z;
      if ((nz = SDDS_RowCount(&SDDSin))<=1)
        bombElegantVA("Too few z values in file %s for BGGEXP %s #%ld\n", filename, tcontext.elementName, tcontext.elementOccurrence);
      if (!(z=SDDS_GetColumnInDoubles(&SDDSin, "z")))
        bombElegantVA("Problem reading column z from %s for BGGEXP %s #%ld\n", buffer, filename, tcontext.elementName, tcontext.elementOccurrence);
      dz0 = z[1] - z[0];
      for (iz=1; iz<nz; iz++) {
        dz = z[iz] - z[iz-1];
        if (dz<=0 || fabs(dz0/dz-1)>1e-6)
          bombElegantVA("Data not uniformly and monotonically increasing in z column from %s for BGGEXP %s #%ld\n", buffer, filename, tcontext.elementName, tcontext.elementOccurrence);
      }
      free(z);
      storedBGGExpData[nBGGExpDataSets].dz = dz0;
    } else {
      if (nz != SDDS_RowCount(&SDDSin))
        bombElegantVA("Inconsistent number of z values in file %s for BGGEXP %s #%ld\n", filename, tcontext.elementName, tcontext.elementOccurrence);
    }
    if (!(storedBGGExpData[nBGGExpDataSets].m = SDDS_Realloc(storedBGGExpData[nBGGExpDataSets].m, 
                                                             sizeof(*storedBGGExpData[nBGGExpDataSets].m)*(im+1))) ||
        !(storedBGGExpData[nBGGExpDataSets].Cmns = SDDS_Realloc(storedBGGExpData[nBGGExpDataSets].Cmns, 
                                                          sizeof(*storedBGGExpData[nBGGExpDataSets].Cmns)*(im+1))) ||
        !(storedBGGExpData[nBGGExpDataSets].dCmns_dz = SDDS_Realloc(storedBGGExpData[nBGGExpDataSets].dCmns_dz, 
                                                                    sizeof(*storedBGGExpData[nBGGExpDataSets].dCmns_dz)*(im+1))))
      bombElegantVA("Memory allocation failure (1) loading data from file %s for BGGEXP %s #%ld\n", filename, tcontext.elementName, tcontext.elementOccurrence);      
    storedBGGExpData[nBGGExpDataSets].Cmns[im] = NULL;
    storedBGGExpData[nBGGExpDataSets].dCmns_dz[im] = NULL;
    if (!(storedBGGExpData[nBGGExpDataSets].Cmns[im] = malloc(sizeof(*storedBGGExpData[nBGGExpDataSets].Cmns[im])*nc)) ||
        !(storedBGGExpData[nBGGExpDataSets].dCmns_dz[im] = malloc(sizeof(*storedBGGExpData[nBGGExpDataSets].dCmns_dz[im])*nc)))
      bombElegantVA("Memory allocation failure (2) loading data from file %s for BGGEXP %s #%ld\n", filename, tcontext.elementName, tcontext.elementOccurrence);      

    for (ic=0; ic<nc; ic++) {
      snprintf(buffer, BUFSIZE, "Cnm%ld", 2*ic);
      if (!(storedBGGExpData[nBGGExpDataSets].Cmns[im][ic] = SDDS_GetColumnInDoubles(&SDDSin, buffer))) 
        bombElegantVA("Problem reading column %s from %s for BGGEXP %s #%ld\n", buffer, filename, tcontext.elementName, tcontext.elementOccurrence);
      snprintf(buffer, BUFSIZE, "dCnm%ld/dz", 2*ic);
      if (!(storedBGGExpData[nBGGExpDataSets].dCmns_dz[im][ic] = SDDS_GetColumnInDoubles(&SDDSin, buffer))) 
        bombElegantVA("Problem reading column %s from %s for BGGEXP %s #%ld\n", buffer, filename, tcontext.elementName, tcontext.elementOccurrence);
    }

    storedBGGExpData[nBGGExpDataSets].m[im] = m;
    im ++;
  }
  SDDS_Terminate(&SDDSin);
  
  storedBGGExpData[nBGGExpDataSets].nm = im;
  storedBGGExpData[nBGGExpDataSets].nz = nz;
  
  hadd(fileHashTable, filename, strlen(filename), nBGGExpDataSets);

  printf("Done adding BGGEXP data from file %s\n", filename);
  fflush(stdout);

  return nBGGExpDataSets++;
}

long trackBGGExpansion(double **part, long np, BGGEXP *bgg, double pCentral, double **accepted)
{
  long ip, ig, im, iz, m, igLimit, izLast;
  STORED_BGGEXP_DATA *bggData;
  double ds, x, y, xp, yp, delta, s, r, phi, denom;
  double gamma, step,  length, fieldLength;

#ifdef DEBUG
  static FILE *fpdebug = NULL;
  if (!fpdebug) {
    fpdebug = fopen("bggexp.deb", "w");
    fprintf(fpdebug, "SDDS1\n");
    fprintf(fpdebug, "&column name=particleID type=long units=m &end\n");
    fprintf(fpdebug, "&column name=ds type=double units=m &end\n");
    fprintf(fpdebug, "&column name=x type=double units=m &end\n");
    fprintf(fpdebug, "&column name=y type=double units=m &end\n");
    fprintf(fpdebug, "&column name=z type=double units=m &end\n");
    fprintf(fpdebug, "&column name=Bx type=double units=T &end\n");
    fprintf(fpdebug, "&column name=By type=double units=T &end\n");
    fprintf(fpdebug, "&column name=Bz type=double units=T &end\n");
    fprintf(fpdebug, "&column name=px type=double &end\n");
    fprintf(fpdebug, "&column name=py type=double &end\n");
    fprintf(fpdebug, "&column name=pz type=double &end\n");
    fprintf(fpdebug, "&column name=dpx type=double &end\n");
    fprintf(fpdebug, "&column name=dpy type=double &end\n");
    fprintf(fpdebug, "&column name=dpz type=double &end\n");
    fprintf(fpdebug, "&data mode=ascii no_row_counts=1 &end\n");
  }
#endif

  TRACKING_CONTEXT tcontext;
  getTrackingContext(&tcontext);
  
  if (!bgg->initialized) {
    bgg->initialized = 1;
    if (!(bgg->filename) || !strlen(bgg->filename)) {
      bombElegantVA("No filename given for BGGEXP %s #%ld\n", tcontext.elementName, tcontext.elementOccurrence);
    }
    bgg->dataIndex = addBGGExpData(bgg->filename);
  }
  bggData = storedBGGExpData+bgg->dataIndex;

  length = bgg->length; 
  fieldLength = bgg->fieldLength;
  if (fieldLength<0)
    fieldLength = length;
  
  step = bgg->fieldLength/(bggData->nz-1);
  if (fabs(step/bggData->dz-1)>1e-6) 
    bombElegantVA("Length mismatch for BGGEXP %s #%ld: %le vs %le\n", tcontext.elementName, tcontext.elementOccurrence, step, bggData->dz);

  /* adjust for insertion length differing from field length */
  if (length!=fieldLength)
    exactDrift(part, np, (length-fieldLength)/2);
  
  /* Do misalignments */
  if (bgg->dx || bgg->dy || bgg->dz)
    offsetBeamCoordinates(part, np, bgg->dx, bgg->dy, bgg->dz);
  if (bgg->tilt)
    rotateBeamCoordinates(part, np, bgg->tilt);

  igLimit = bggData->nGradients;
  if (bgg->maximum2n>=0) {
    igLimit = bgg->maximum2n/2+1;
    if (igLimit>bggData->nGradients)
      igLimit = bggData->nGradients;
  }
  
  if (bgg->zInterval<=0) 
    bombElegantVA("zInterval %ld is invalid for BGGEXP %s #%ld\n", bgg->zInterval, tcontext.elementName, tcontext.elementOccurrence);
  izLast = bggData->nz-bgg->zInterval;

  if (bgg->symplectic) {
    long iImpLoop;
    double xMid, yMid,xNext, yNext, xLoop, yLoop, delta_s;
    double px, py, pxNext, pyNext, pxLoop, pyLoop, ux, uy;
    double delta, r, phi, denom, scaleA;
    double Ax, dAx_dx, dAx_dy, Ay, dAy_dx, dAy_dy, dAz_dx, dAz_dy;
    double GenGrad_s, dGenGrad_s;
    double epsImplConverge = 1.e-14*1.0e-3;

    scaleA = -bgg->strength*particleCharge*particleRelSign/(pCentral*particleMass*c_mks);  /** [factor in parentheses of a = (q/p_0)*A] **/
    /* Element body */
    for (ip=0; ip<np; ip++) {
      x = part[ip][0];
      y = part[ip][2];
      delta = part[ip][5];
      denom = 1.0/sqrt(1.0 + part[ip][1]*part[ip][1] + part[ip][3]*part[ip][3]);
      px = part[ip][1]*(1.0 + delta)*denom;
      py = part[ip][3]*(1.0 + delta)*denom;
      s = part[ip][4];
      delta_s = 0.0;
      
      /* Integrate through the magnet */
      for (iz=0; iz<bggData->nz; iz+=bgg->zInterval) {
        r = sqrt(sqr(x)+sqr(y));
        phi = atan2(y, x);
        
        /** Calculate vector potential A and its relevant derivatives from the generalized gradients **/
        Ax = dAx_dx = dAx_dy = Ay = dAy_dx = dAy_dy = dAz_dx = dAz_dy = 0.0;
        /* Compute fields */
        for (im=0; im<bggData->nm; im++) {
          double m_1fact, term, sin_mphi, cos_mphi;
          m = bggData->m[im];
          if (bgg->mMaximum>0 && m>bgg->mMaximum)
            continue;
          m_1fact = dfactorial(m-1);
          sin_mphi = sin(m*phi);
          cos_mphi = cos(m*phi);
          for (ig=0; ig<igLimit; ig++) {
            term  = ipow(-1, ig)*m_1fact*ipow(r, 2*ig+m-2)/(ipow(4, ig)*factorial(ig)*factorial(ig+m));
            dGenGrad_s = bggData->dCmns_dz[im][ig][iz];
            GenGrad_s = bggData->Cmns[im][ig][iz];
            /** Assume skew components Cmnc, dCmnc_dz = 0 **/
            /** Ax += term*( cos_mphi*GenGrad_s -  sin_mphi*GenGrad_c ); **/
            Ax += term*cos_mphi*dGenGrad_s;
            /** dAx_dx += term*( (((2*ig+m+1)*x*x + y*y)*cos_mphi + m*x*y*sin_mphi)*dGenGrad_s - (((2*ig+m+1)*x*x + y*y)*sin_mphi - m*x*y*cos_mphi )*dGenGrad_c ); **/
            dAx_dx += term*( ((2*ig+m+1)*x*x + y*y)*cos_mphi + m*x*y*sin_mphi )*dGenGrad_s;
            /** dAx_dy += term*( ((2*ig+m)*x*y*cos_mphi - m*x*x*sin_mphi)*dGenGrad_s - ((2*ig+m)*x*y*sin_mphi + m*x*x*cos_mphi )*dGenGrad_c ); **/
            dAx_dy += term*( (2*ig+m)*x*y*cos_mphi - m*x*x*sin_mphi )*dGenGrad_s;
            /** Ay += term*( cos_mphi*dGenGrad_s -  sin_mphi*dGenGrad_c ); **/
            Ay += term*cos_mphi*dGenGrad_s;
            /** dAy_dx += term*( ((2*ig+m)*x*y*cos_mphi + m*y*y*sin_mphi)*dGenGrad_s - ((2*ig+m)*x*y*sin_mphi - m*y*y*cos_mphi )*dGenGrad_c ); **/
            dAy_dx += term*( (2*ig+m)*x*y*cos_mphi + m*y*y*sin_mphi )*dGenGrad_s;
            /** dAy_dy += term*( (((2*ig+m+1)*y*y + x*x)*cos_mphi - m*x*y*sin_mphi)*dGenGrad_s - (((2*ig+m+1)*y*y + x*x)*sin_mphi + m*x*y*cos_mphi )*dGenGrad_c ); **/
            dAy_dy += term*( ((2*ig+m+1)*y*y + x*x)*cos_mphi - m*x*y*sin_mphi )*dGenGrad_s;
            /** dAz_dx += (2*ig+m)*term*( ((2*ig+m)*x*sin_mphi - m*y*cos_mphi)*dGenGrad_c - ((2*ig+m)*x*cos_mphi + m*y*sin_mphi)*dGenGrad_s ); **/
            dAz_dx +=-(2*ig+m)*term*( (2*ig+m)*x*cos_mphi + m*y*sin_mphi )*GenGrad_s;
            /** dAz_dy += (2*ig+m)*term*( ((2*ig+m)*y*sin_mphi + m*x*cos_mphi)*dGenGrad_c - ((2*ig+m)*y*cos_mphi - m*x*sin_mphi)*dGenGrad_s ); **/
            dAz_dy +=-(2*ig+m)*term*( (2*ig+m)*y*cos_mphi - m*x*sin_mphi )*GenGrad_s;
          }
        }
        Ax = x*r*r*Ax;
        Ay = y*r*r*Ay;
        
        /** Start with first order guess for the 'Next' coordinates **/
        ux = px - scaleA*Ax;
        uy = py - scaleA*Ay;
        denom = (1.0 + delta)*(1.0 + delta) - ux*ux - uy*uy;
        denom = 1.0/sqrt(denom);
        xNext = x + step*bgg->zInterval*ux*denom;
        yNext = y + step*bgg->zInterval*uy*denom;
        pxNext = px + step*bgg->zInterval*scaleA*( (ux*dAx_dx + uy*dAy_dx)*denom + dAz_dx );
        pyNext = py + step*bgg->zInterval*scaleA*( (ux*dAx_dy + uy*dAy_dy)*denom + dAz_dy );
        iImpLoop=0;
        do {
          /** Convergence determined when updated 'Next' coordinates match previous 'Loop' coords.  **/
          xLoop = xNext;
          yLoop = yNext;
          pxLoop = pxNext;
          pyLoop = pyNext;
          xMid = 0.5*(x + xLoop);
          yMid = 0.5*(y + yLoop);
          r = sqrt(sqr(xMid)+sqr(yMid));
          phi = atan2(yMid, xMid);
          /** Calculate vector potential A and its relevant derivatives from the generalized gradients **/
          Ax = dAx_dx = dAx_dy = Ay = dAy_dx = dAy_dy = dAz_dx = dAz_dy = 0.0;
          for (im=0; im<bggData->nm; im++) {
            double m_1fact, term, sin_mphi, cos_mphi;
            m = bggData->m[im];
            if (bgg->mMaximum>0 && m>bgg->mMaximum)
              continue;
            m_1fact = dfactorial(m-1);
            sin_mphi = sin(m*phi);
            cos_mphi = cos(m*phi);
            for (ig=0; ig<igLimit; ig++) {
              term  = ipow(-1, ig)*m_1fact*ipow(r, 2*ig+m-2)/(ipow(4, ig)*factorial(ig)*factorial(ig+m));
              dGenGrad_s = bggData->dCmns_dz[im][ig][iz];
              GenGrad_s = bggData->Cmns[im][ig][iz];
              /** Assume skew components Cmnc, dCmnc_dz = 0 **/
              /** Ax += term*( cos_mphi*dGenGrad_s -  sin_mphi*dGenGrad_c ); **/
              Ax += term*cos_mphi*dGenGrad_s;
              /** dAx_dx += term*( (((2*ig+m+1)*xMid*xMid + yMid*yMid)*cos_mphi + m*xMid*yMid*sin_mphi)*dGenGrad_s - (((2*ig+m+1)*xMid*xMid + yMid*yMid)*sin_mphi - m*xMid*yMid*cos_mphi )*dGenGrad_c ); **/
              dAx_dx += term*( ((2*ig+m+1)*xMid*xMid + yMid*yMid)*cos_mphi + m*xMid*yMid*sin_mphi )*dGenGrad_s;
              /** dAx_dy += term*( ((2*ig+m)*xMid*yMid*cos_mphi - m*xMid*xMid*sin_mphi)*dGenGrad_s - ((2*ig+m)*xMid*yMid*sin_mphi + m*xMid*xMid*cos_mphi )*dGenGrad_c ); **/
              dAx_dy += term*( (2*ig+m)*xMid*yMid*cos_mphi - m*xMid*xMid*sin_mphi )*dGenGrad_s;
              /** Ay += term*( cos_mphi*dGenGrad_s -  sin_mphi*dGenGrad_c ); **/
              Ay += term*cos_mphi*dGenGrad_s;
              /** dAy_dx += term*( ((2*ig+m)*xMid*yMid*cos_mphi + m*yMid*yMid*sin_mphi)*dGenGrad_s - ((2*ig+m)*xMid*yMid*sin_mphi - m*yMid*yMid*cos_mphi )*dGenGrad_c ); **/
              dAy_dx += term*( (2*ig+m)*xMid*yMid*cos_mphi + m*yMid*yMid*sin_mphi )*dGenGrad_s;
              /** dAy_dy += term*( (((2*ig+m+1)*yMid*yMid + xMid*xMid)*cos_mphi - m*xMid*yMid*sin_mphi)*dGenGrad_s - (((2*ig+m+1)*yMid*yMid + xMid*xMid)*sin_mphi + m*xMid*yMid*cos_mphi )*dGenGrad_c ); **/
              dAy_dy += term*( ((2*ig+m+1)*yMid*yMid + xMid*xMid)*cos_mphi - m*xMid*yMid*sin_mphi )*dGenGrad_s;
              /** dAz_dx += (2*ig+m)*term*( ((2*ig+m)*xMid*sin_mphi - m*yMid*cos_mphi)*dGenGrad_c - ((2*ig+m)*xMid*cos_mphi + m*yMid*sin_mphi)*dGenGrad_s ); **/
              dAz_dx +=-(2*ig+m)*term*( (2*ig+m)*xMid*cos_mphi + m*yMid*sin_mphi )*GenGrad_s;
              /** dAz_dy += (2*ig+m)*term*( ((2*ig+m)*yMid*sin_mphi + m*xMid*cos_mphi)*dGenGrad_c - ((2*ig+m)*yMid*cos_mphi - m*xMid*sin_mphi)*dGenGrad_s ); **/
              dAz_dy +=-(2*ig+m)*term*( (2*ig+m)*yMid*cos_mphi - m*xMid*sin_mphi )*GenGrad_s;
            }
          }
          Ax = xMid*r*r*Ax;
          Ay = yMid*r*r*Ay;
          
          /** Update coordinates **/
          ux = 0.5*(px + pxLoop) - scaleA*Ax;
          uy = 0.5*(py + pyLoop) - scaleA*Ay;
          denom = (1.0 + delta)*(1.0 + delta) - ux*ux - uy*uy;
          denom = 1.0/sqrt(denom);
          xNext = x + step*bgg->zInterval*ux*denom;
          yNext = y + step*bgg->zInterval*uy*denom;
          pxNext = px + step*bgg->zInterval*scaleA*( (ux*dAx_dx + uy*dAy_dx)*denom + dAz_dx );
          pyNext = py + step*bgg->zInterval*scaleA*( (ux*dAx_dy + uy*dAy_dy)*denom + dAz_dy );
          iImpLoop++;
          if(iImpLoop>10) {
            //printf("HERE: daz_dx = %e\t By = %e\n", scaleA*dAz_dx, B[2]);
            break;
          }
        } while ( (fabs(xNext - xLoop) > epsImplConverge) || (fabs(yNext - yLoop) > epsImplConverge) ||
                 (fabs(pxNext - pxLoop) > epsImplConverge) || (fabs(pyNext - pyLoop) > epsImplConverge) );
        x = xNext;
        y = yNext;
        px = pxNext;
        py = pyNext;
        delta_s += step*bgg->zInterval*((1.0 + delta)*denom - 1.0);
        s += step*bgg->zInterval;

#ifdef DEBUG
        fprintf(fpdebug, "%le %le %le %le %le %le %le %le %i\n", 
                iz*bggData->dz, x, y, px, py, s, dAz_dy, dAz_dx, iImpLoop);
#endif

      }
 
      part[ip][0] = x;
      denom = (1.0 + delta)*(1.0 + delta) - px*px - py*py;
      denom = 1.0/sqrt(denom);
      part[ip][1] = px*denom;
      part[ip][2] = y;
      part[ip][3] = py*denom;
      part[ip][4] = s + delta_s;
    }
  }  else { 
    /* Element body */
    double B[3], p[3], dp[3], Bphi, Br;
    for (ip=0; ip<np; ip++) {
      x = part[ip][0];
      xp = part[ip][1];
      y = part[ip][2];
      yp = part[ip][3];
      s = part[ip][4];
      delta = part[ip][5];
      
      /* compute momenta (x, y, z) */
      denom = sqrt(1 + sqr(xp) + sqr(yp));
      p[2] = pCentral*(1+delta)/denom;
      p[0] = xp*p[2];
      p[1] = yp*p[2];
      gamma = sqrt(sqr(p[0]) + sqr(p[1]) + sqr(p[2]) + 1);
      
      /* Integrate through the magnet */
      for (iz=0; iz<bggData->nz; iz+=bgg->zInterval) {
        r = sqrt(sqr(x)+sqr(y));
        phi = atan2(y, x);
        
        /* Compute fields */
        Br = Bphi = B[2] = 0;
        for (im=0; im<bggData->nm; im++) {
          double mfact, term, sin_mphi, cos_mphi;
          m = bggData->m[im];
          if (bgg->mMaximum>0 && m>bgg->mMaximum)
            continue;
          mfact = dfactorial(m);
          sin_mphi = sin(m*phi);
          cos_mphi = cos(m*phi);
          for (ig=0; ig<igLimit; ig++) {
            term  = ipow(-1, ig)*mfact/(ipow(2, 2*ig)*factorial(ig)*factorial(ig+m))*ipow(r, 2*ig+m-1);
            B[2] += term*bggData->dCmns_dz[im][ig][iz]*r*sin_mphi;
            term *= bggData->Cmns[im][ig][iz];
            Br   += term*(2*ig+m)*sin_mphi;
            Bphi += m*term*cos_mphi;
          }
        }
        B[0] = (Br*cos(phi) - Bphi*sin(phi))*bgg->strength;
        B[1] = (Br*sin(phi) + Bphi*cos(phi))*bgg->strength;
        B[2] *= bgg->strength;
        
        /* Apply kicks */
        ds = step*bgg->zInterval*sqrt(1 + sqr(p[0]/p[2]) + sqr(p[1]/p[2]));
        dp[0] = -particleCharge*particleRelSign*ds/(particleMass*gamma*c_mks)*(p[1]*B[2] - p[2]*B[1]);
        dp[1] = -particleCharge*particleRelSign*ds/(particleMass*gamma*c_mks)*(p[2]*B[0] - p[0]*B[2]);
        dp[2] = -particleCharge*particleRelSign*ds/(particleMass*gamma*c_mks)*(p[0]*B[1] - p[1]*B[0]);
        
#ifdef DEBUG
        fprintf(fpdebug, "%.0f %le %le %le %le %le %le %le %le %le %le %le %le %le\n", 
                part[ip][6], ds, x, y, iz*bggData->dz, 
                B[0], B[1], B[2], 
                p[0], p[1], p[2], dp[0], dp[1], dp[2]);
#endif
        
        p[0] += dp[0];
        p[1] += dp[1];
        p[2] += dp[2];
        
        if (iz<izLast) {
          /* Drift forward */
          x += p[0]/p[2]*bggData->dz*bgg->zInterval;
          y += p[1]/p[2]*bggData->dz*bgg->zInterval;
          s += ds;
        }
      }
      if (iz<bggData->nz) {
        /* Drift forward */
        x += p[0]/p[2]*bggData->dz*(bggData->nz-1-(iz-1));
        y += p[1]/p[2]*bggData->dz*(bggData->nz-1-(iz-1));
        s += bggData->dz*(bggData->nz-1-iz)*sqrt(1+sqr(p[0]/p[2])+sqr(p[1]/p[2]));
      }
      
      part[ip][0] = x;
      part[ip][1] = p[0]/p[2];
      part[ip][2] = y;
      part[ip][3] = p[1]/p[2];
      part[ip][4] = s;
    }
  }


  /* Do misalignments */
  if (bgg->tilt)
    rotateBeamCoordinates(part, np, -bgg->tilt);
  if (bgg->dx || bgg->dy || bgg->dz)
    offsetBeamCoordinates(part, np, -bgg->dx, -bgg->dy, -bgg->dz);

  /* adjust for insertion length differing from field length */
  length = bgg->length; 
  fieldLength = bgg->fieldLength;
  if (length!=fieldLength)
    exactDrift(part, np, (length-fieldLength)/2);

#ifdef DEBUG  
  fflush(fpdebug);
#endif

  return np;
}

