/* Copyright 2016 by Michael Borland and Argonne National Laboratory,
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

  im = 0;
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
  long ip, ig, im, iz, m, igLimit, mMax;
  STORED_BGGEXP_DATA *bggData;
  double ds, x, y, xp, yp, delta, s, r, phi, denom;
  double B[3], p[3], dp[3], Bphi, Br, gamma, step,  length, fieldLength;
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
  if (bgg->nMaximum>0) {
    igLimit = bgg->nMaximum/2+1;
    if (igLimit>bggData->nGradients)
      igLimit = bggData->nGradients;
  }
    
  /* Element body */
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
    for (iz=0; iz<bggData->nz; iz++) {
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
      ds = step*sqrt(1 + sqr(p[0]/p[2]) + sqr(p[1]/p[2]));
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

      if (iz!=(bggData->nz-1)) {
        /* Drift forward */
        x += p[0]/p[2]*ds;
        y += p[1]/p[2]*ds;
        s += ds;
      }
    }
    part[ip][0] = x;
    part[ip][1] = p[0]/p[2];
    part[ip][2] = y;
    part[ip][3] = p[1]/p[2];
    part[ip][4] = s;
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

