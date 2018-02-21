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
  double zMin, zMax;   /* minimum and maximum z values */
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
  long *nstore;
  
  if (!fileHashTable)
    fileHashTable = hcreate(12);

  if (hfind(fileHashTable, filename, strlen(filename))) {
    printf("Using previously-stored BGGEXP data for filename %s\n", filename);
    fflush(stdout);
    im = *((long*)hstuff(fileHashTable));
    printf("Using BGGEXP table %ld for filename %s\n", im, filename);
    fflush(stdout);
    return im;
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
  storedBGGExpData[nBGGExpDataSets].zMin = DBL_MAX;
  storedBGGExpData[nBGGExpDataSets].zMax = -DBL_MAX;
  while ((readCode=SDDS_ReadPage(&SDDSin))>0) {
    if (!SDDS_GetParameter(&SDDSin, "m", &m) || m<1)
      bombElegantVA("Problem with value of m for page %ld of file %s for BGGEXP %s #%ld\n", readCode, filename, tcontext.elementName, tcontext.elementOccurrence);
    if (readCode==1) {
      long iz;
      double dz0, dz, *z, zMin, zMax;
      if ((nz = SDDS_RowCount(&SDDSin))<=1)
        bombElegantVA("Too few z values in file %s for BGGEXP %s #%ld\n", filename, tcontext.elementName, tcontext.elementOccurrence);
      if (!(z=SDDS_GetColumnInDoubles(&SDDSin, "z")))
        bombElegantVA("Problem reading column z from %s for BGGEXP %s #%ld\n", buffer, filename, tcontext.elementName, tcontext.elementOccurrence);
      find_min_max(&zMin, &zMax, z, nz);
      if (zMin<storedBGGExpData[nBGGExpDataSets].zMin)
        storedBGGExpData[nBGGExpDataSets].zMin = zMin;
      if (zMax>storedBGGExpData[nBGGExpDataSets].zMax)
        storedBGGExpData[nBGGExpDataSets].zMax = zMax;
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

  nstore = tmalloc(sizeof(*nstore));
  *nstore = nBGGExpDataSets;
  hadd(fileHashTable, filename, strlen(filename), (void*)nstore);

  printf("Done adding BGGEXP data from file %s\n", filename);
  fflush(stdout);

  return nBGGExpDataSets++;
}

long trackBGGExpansion(double **part, long np, BGGEXP *bgg, double pCentral, double **accepted, double *sigmaDelta2)
{
  long ip, ig, im, iz, irow, m, igLimit, izLast;
  STORED_BGGEXP_DATA *bggData;
  double ds, dz, x, y, xp, yp, delta, s, r, phi, denom;
  double gamma, step,  length;
  TRACKING_CONTEXT tcontext;
  double radCoef=0, isrCoef=0;
  double zMin, zMax, xVertex, zVertex, xEntry, zEntry, xExit, zExit;

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
  if (bgg->synchRad) {
    radCoef = ipow(particleCharge,4)/(6*PI*epsilon_o*ipow(c_mks,4)*ipow(particleMass,3));
    isrCoef = sqrt(55/(24*sqrt(3))*particleRadius*hbar_mks*ipow(particleCharge, 3))/sqr(particleMass*c_mks);
  }
  if (sigmaDelta2)
    *sigmaDelta2 = 0;

  getTrackingContext(&tcontext);

  if (!bgg->initialized) {
    //char *outputFile;
    bgg->initialized = 1;
    if (!(bgg->filename) || !strlen(bgg->filename)) {
      bombElegantVA("No filename given for BGGEXP %s #%ld\n", tcontext.elementName, tcontext.elementOccurrence);
    }
    bgg->dataIndex = addBGGExpData(bgg->filename);
#if !USE_MPI
    if (bgg->particleOutputFile && !bgg->SDDSpo) {
      bgg->SDDSpo = tmalloc(sizeof(*(bgg->SDDSpo)));
      bgg->particleOutputFile = compose_filename(bgg->particleOutputFile, tcontext.rootname);
      if (!SDDS_InitializeOutput(bgg->SDDSpo, SDDS_BINARY, 1, 
                                 NULL, NULL, bgg->particleOutputFile) ||
          0>SDDS_DefineParameter(bgg->SDDSpo, "SVNVersion", NULL, NULL, "SVN version number", NULL, SDDS_STRING, SVN_VERSION) ||
          !SDDS_DefineSimpleParameter(bgg->SDDSpo, "particleID", NULL, SDDS_LONG) ||
          !SDDS_DefineSimpleParameter(bgg->SDDSpo, "pCentral", "m$be$nc", SDDS_DOUBLE) ||
          (bgg->poIndex[0]=SDDS_DefineColumn(bgg->SDDSpo, "x", NULL, "m", NULL, NULL, SDDS_DOUBLE, 0 ))<0 ||
          (bgg->poIndex[1]=SDDS_DefineColumn(bgg->SDDSpo, "px", NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0))<0 ||
          (bgg->poIndex[2]=SDDS_DefineColumn(bgg->SDDSpo, "y", NULL, "m", NULL, NULL, SDDS_DOUBLE, 0 ))<0 ||
          (bgg->poIndex[3]=SDDS_DefineColumn(bgg->SDDSpo, "py", NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0))<0 ||
          (bgg->poIndex[4]=SDDS_DefineColumn(bgg->SDDSpo, "z", NULL, "m", NULL, NULL, SDDS_DOUBLE, 0))<0 ||
          (bgg->poIndex[5]=SDDS_DefineColumn(bgg->SDDSpo, "pz", NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0))<0 ||
          (bgg->poIndex[6]=SDDS_DefineColumn(bgg->SDDSpo, "Bx", NULL, "T", NULL, NULL, SDDS_DOUBLE, 0))<0 ||
          (bgg->poIndex[7]=SDDS_DefineColumn(bgg->SDDSpo, "By", NULL, "T", NULL, NULL, SDDS_DOUBLE, 0))<0 ||
          (bgg->poIndex[8]=SDDS_DefineColumn(bgg->SDDSpo, "Bz", NULL, "T", NULL, NULL, SDDS_DOUBLE, 0))<0 ||
          !SDDS_WriteLayout(bgg->SDDSpo)) {
        SDDS_SetError("Problem setting up particle output file for BGGEXP");
        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
      }
    }
#endif
  }
  bggData = storedBGGExpData+bgg->dataIndex;

  length = bgg->length; 
  
  zMin = bggData->zMin;
  zMax = bggData->zMax;
  
  if (bgg->isBend) {
    xVertex = bgg->xVertex;
    xEntry = bgg->xEntry;
    xExit = bgg->xExit;
    zVertex = bgg->zVertex;
    zEntry = bgg->zEntry;
    zExit = bgg->zExit;
  } else {
    xEntry = xExit = xVertex = 0;
    zVertex = (zMax+zMin)/2;
    zEntry = zVertex - length/2;
    zExit = zVertex + length/2;
    /*
    printf("fieldLength = %le, z:[%le, %le]\n", fieldLength, zMin, zMax);
    printf("vertex: (%le, %le)\n", zVertex, xVertex);
    printf("entry: (%le, %le)\n", zEntry, xEntry);
    printf("exit: (%le, %le)\n", zExit, xExit);
    */
  }

  step = (zMax-zMin)/(bggData->nz-1);
  /*
  if (fabs(step/bggData->dz-1)>1e-6) 
    bombElegantVA("Length mismatch for BGGEXP %s #%ld: %le vs %le\nlength=%le m, nz=%ld\n",
                  tcontext.elementName, tcontext.elementOccurrence, step, bggData->dz,
                  bgg->fieldLength, bggData->nz);
  */

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
    double px, py, pxNext, pyNext, pxLoop, pyLoop, ux=0, uy=0;
    double delta, r, phi, denom, scaleA, sin_phi, cos_phi;
    double Ax, dAx_dx, dAx_dy, Ay, dAy_dx, dAy_dy, dAz_dx, dAz_dy;
    double GenGrad_s, dGenGrad_s;
    double epsImplConverge = 1.e-14*1.0e-3;
    double Bx, By, Bz;

    double magnet_s;

    scaleA = -bgg->strength*particleCharge*particleRelSign/(pCentral*particleMass*c_mks);  /** [factor in parentheses of a = (q/p_0)*A] **/
    /* Element body */
    for (ip=0; ip<np; ip++) {
#if !USE_MPI
      if (bgg->SDDSpo) {
        if (!SDDS_StartPage(bgg->SDDSpo, bggData->nz+1) ||
            !SDDS_SetParameters(bgg->SDDSpo, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,
                                "particleID", (long)(part[ip][6]), "pCentral", pCentral, NULL)) {
          SDDS_SetError("Problem setting up particle output page for BGGEXP");
          SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
        }
      }
#endif

      /* Transform to canonical coordinates */
      delta = part[ip][5];
      denom = 1.0/sqrt(1.0 + part[ip][1]*part[ip][1] + part[ip][3]*part[ip][3]);
      x = part[ip][0];
      px = part[ip][1]*(1.0 + delta)*denom;
      y = part[ip][2];
      py = part[ip][3]*(1.0 + delta)*denom;
      s = part[ip][4];
      delta_s = 0.0;
      Bx = By = Bz = 0;

      /* Compute transformation from canonical accelerator coordinates -> magnet grid if bending magnet */
      denom = (1.0 + delta)*(1.0 + delta) - px*px - py*py;
      denom = 1.0/sqrt(denom);
      if (bgg->isBend) {
        phi = atan2(xVertex-xEntry, zVertex-zEntry);  /* angle of incoming x-plane w.r.t. magnet x */
        sin_phi = sin(phi);
        cos_phi = cos(phi);
        magnet_s = x*sin_phi/( cos_phi - sin_phi*px*denom );
      } else {
        phi = sin_phi = magnet_s = 0;
        cos_phi = 1;
      }

      x = (x + px*magnet_s*denom*denom)/cos_phi;
      y = (y + py*magnet_s*denom*denom);
      px = px*cos_phi + sin_phi/denom;
      delta_s = (1 + delta)*magnet_s*denom;

      x += xEntry + bgg->dxExpansion; /* Shift x according to entrance point & map shift */

      /* Drift backward from entrance plane to beginning of field map */
      denom = (1.0 + delta)*(1.0 + delta) - px*px - py*py;
      denom = 1.0/sqrt(denom);
      x -= (zEntry - zMin)*px*denom;
      y -= (zEntry - zMin)*py*denom;
      s -= (zEntry - zMin)*(1.0 + delta)*denom;

      /* Drift particles back 1/2 step from start of field map to prepare for implicit midpoint rule integration */
      denom = (1.0 + delta)*(1.0 + delta) - px*px - py*py;
      denom = 1.0/sqrt(denom);
      x -= 0.5*step*bgg->zInterval*px*denom;
      y -= 0.5*step*bgg->zInterval*py*denom;
      s -= 0.5*step*bgg->zInterval*(1.0 + delta)*denom;
     /* Integrate through the magnet */
      for (iz=irow=0; iz<bggData->nz; iz+=bgg->zInterval) {
#if !USE_MPI
        if (bgg->SDDSpo &&
            !SDDS_SetRowValues(bgg->SDDSpo, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, irow++,
                               bgg->poIndex[0], x,
                               bgg->poIndex[1], px*pCentral,
                               bgg->poIndex[2], y,
                               bgg->poIndex[3], py*pCentral,
                               bgg->poIndex[4], iz*bggData->dz,
                               bgg->poIndex[5], sqrt(sqr(pCentral*(1+delta))-(sqr(px)+sqr(py))*sqr(pCentral)),
                               bgg->poIndex[6], Bx,
                               bgg->poIndex[7], By,
                               bgg->poIndex[8], Bz,
                               -1)) {
          SDDS_SetError("Problem setting particle output data for BGGEXP");
          SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
        }
#endif
        r = sqrt(sqr(x)+sqr(y));
        phi = atan2(y, x);
        cos_phi = cos(phi);
        sin_phi = sin(phi);

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
            term  = ipow(-1, ig)*m_1fact*ipow(r, 2*ig+m-1)/(ipow(4, ig)*factorial(ig)*factorial(ig+m));
            dGenGrad_s = bggData->dCmns_dz[im][ig][iz];
            GenGrad_s = bggData->Cmns[im][ig][iz];
            /** Assume skew components Cmnc, dCmnc_dz = 0 **/
            /** Ax += term*( cos_mphi*GenGrad_s -  sin_mphi*GenGrad_c ); **/
            Ax += term*cos_mphi*dGenGrad_s;
            /** dAx_dx += term*( (((2*ig+m+1)*x*x + y*y)*cos_mphi + m*x*y*sin_mphi)*dGenGrad_s - (((2*ig+m+1)*x*x + y*y)*sin_mphi - m*x*y*cos_mphi )*dGenGrad_c ); **/
            dAx_dx += term*( ((2*ig+m)*cos_phi*x + r)*cos_mphi + m*cos_phi*y*sin_mphi )*dGenGrad_s;
            /** dAx_dy += term*( ((2*ig+m)*x*y*cos_mphi - m*x*x*sin_mphi)*dGenGrad_s - ((2*ig+m)*x*y*sin_mphi + m*x*x*cos_mphi )*dGenGrad_c ); **/
            dAx_dy += term*cos_phi*( (2*ig+m)*y*cos_mphi - m*x*sin_mphi )*dGenGrad_s;
            /** Ay += term*( cos_mphi*dGenGrad_s -  sin_mphi*dGenGrad_c ); **/
            Ay += term*cos_mphi*dGenGrad_s;
            /** dAy_dx += term*( ((2*ig+m)*x*y*cos_mphi + m*y*y*sin_mphi)*dGenGrad_s - ((2*ig+m)*x*y*sin_mphi - m*y*y*cos_mphi )*dGenGrad_c ); **/
            dAy_dx += term*sin_phi*( (2*ig+m)*x*cos_mphi + m*y*sin_mphi )*dGenGrad_s;
            /** dAy_dy += term*( (((2*ig+m+1)*y*y + x*x)*cos_mphi - m*x*y*sin_mphi)*dGenGrad_s - (((2*ig+m+1)*y*y + x*x)*sin_mphi + m*x*y*cos_mphi )*dGenGrad_c ); **/
            dAy_dy += term*( ((2*ig+m)*sin_phi*y + r)*cos_mphi - m*sin_phi*x*sin_mphi )*dGenGrad_s;
            /** dAz_dx += (2*ig+m)*term*( ((2*ig+m)*x*sin_mphi - m*y*cos_mphi)*dGenGrad_c - ((2*ig+m)*x*cos_mphi + m*y*sin_mphi)*dGenGrad_s ); **/
            dAz_dx +=-(2*ig+m)*term*( (2*ig+m)*cos_phi*cos_mphi + m*sin_phi*sin_mphi )*GenGrad_s;
            /** dAz_dy += (2*ig+m)*term*( ((2*ig+m)*y*sin_mphi + m*x*cos_mphi)*dGenGrad_c - ((2*ig+m)*y*cos_mphi - m*x*sin_mphi)*dGenGrad_s ); **/
            dAz_dy +=-(2*ig+m)*term*( (2*ig+m)*sin_phi*cos_mphi - m*cos_phi*sin_mphi )*GenGrad_s;
          }
        }
        Ax = x*r*Ax;
        Ay = y*r*Ay;
	dAz_dx = dAz_dx - bgg->By;
	dAz_dy = dAz_dy + bgg->Bx;

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
	  cos_phi = cos(phi);
 	  sin_phi = sin(phi);
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
              term  = ipow(-1, ig)*m_1fact*ipow(r, 2*ig+m-1)/(ipow(4, ig)*factorial(ig)*factorial(ig+m));
              dGenGrad_s = bggData->dCmns_dz[im][ig][iz];
              GenGrad_s = bggData->Cmns[im][ig][iz];
              /** Assume skew components Cmnc, dCmnc_dz = 0 **/
              /** Ax += term*( cos_mphi*dGenGrad_s -  sin_mphi*dGenGrad_c ); **/
              Ax += term*cos_mphi*dGenGrad_s;
              /** dAx_dx += term*( (((2*ig+m+1)*xMid*xMid + yMid*yMid)*cos_mphi + m*xMid*yMid*sin_mphi)*dGenGrad_s - (((2*ig+m+1)*xMid*xMid + yMid*yMid)*sin_mphi - m*xMid*yMid*cos_mphi )*dGenGrad_c ); **/
              dAx_dx += term*( ((2*ig+m)*cos_phi*xMid + r)*cos_mphi + m*cos_phi*yMid*sin_mphi )*dGenGrad_s;
              /** dAx_dy += term*( ((2*ig+m)*xMid*yMid*cos_mphi - m*xMid*xMid*sin_mphi)*dGenGrad_s - ((2*ig+m)*xMid*yMid*sin_mphi + m*xMid*xMid*cos_mphi )*dGenGrad_c ); **/
              dAx_dy += term*cos_phi*( (2*ig+m)*yMid*cos_mphi - m*xMid*sin_mphi )*dGenGrad_s;
              /** Ay += term*( cos_mphi*dGenGrad_s -  sin_mphi*dGenGrad_c ); **/
              Ay += term*cos_mphi*dGenGrad_s;
              /** dAy_dx += term*( ((2*ig+m)*xMid*yMid*cos_mphi + m*yMid*yMid*sin_mphi)*dGenGrad_s - ((2*ig+m)*xMid*yMid*sin_mphi - m*yMid*yMid*cos_mphi )*dGenGrad_c ); **/
              dAy_dx += term*sin_phi*( (2*ig+m)*xMid*cos_mphi + m*yMid*sin_mphi )*dGenGrad_s;
              /** dAy_dy += term*( (((2*ig+m+1)*yMid*yMid + xMid*xMid)*cos_mphi - m*xMid*yMid*sin_mphi)*dGenGrad_s - (((2*ig+m+1)*yMid*yMid + xMid*xMid)*sin_mphi + m*xMid*yMid*cos_mphi )*dGenGrad_c ); **/
              dAy_dy += term*( ((2*ig+m)*sin_phi*yMid + r)*cos_mphi - m*sin_phi*xMid*sin_mphi )*dGenGrad_s;
              /** dAz_dx += (2*ig+m)*term*( ((2*ig+m)*xMid*sin_mphi - m*yMid*cos_mphi)*dGenGrad_c - ((2*ig+m)*xMid*cos_mphi + m*yMid*sin_mphi)*dGenGrad_s ); **/
              dAz_dx +=-(2*ig+m)*term*( (2*ig+m)*cos_phi*cos_mphi + m*sin_phi*sin_mphi )*GenGrad_s;
              /** dAz_dy += (2*ig+m)*term*( ((2*ig+m)*yMid*sin_mphi + m*xMid*cos_mphi)*dGenGrad_c - ((2*ig+m)*yMid*cos_mphi - m*xMid*sin_mphi)*dGenGrad_s ); **/
              dAz_dy +=-(2*ig+m)*term*( (2*ig+m)*sin_phi*cos_mphi - m*cos_phi*sin_mphi )*GenGrad_s;
            }
          }
	  Ax = xMid*r*Ax;
	  Ay = yMid*r*Ay;
	  dAz_dx = dAz_dx - bgg->By;
	  dAz_dy = dAz_dy + bgg->Bx;
          
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

        if (bgg->synchRad || bgg->SDDSpo) {
	  ds = step*bgg->zInterval*(1.0 + delta)*denom;
	  /* compute Bx, By */
	  Bx = bgg->Bx;
          By = bgg->By;
          for (im=0; im<bggData->nm; im++) {
            double mfact, term, sin_mphi, cos_mphi;
            m = bggData->m[im];
            if (bgg->mMaximum>0 && m>bgg->mMaximum)
              continue;
            mfact = dfactorial(m);
            sin_mphi = sin(m*phi);
            cos_mphi = cos(m*phi);
            for (ig=0; ig<igLimit; ig++) {
              term  = ipow(-1, ig)*mfact*ipow(r, 2*ig+m-1)/(ipow(4, ig)*factorial(ig)*factorial(ig+m));
              GenGrad_s = bggData->Cmns[im][ig][iz];
              /** Assume skew components Cmnc = 0 **/
              Bx += term*( (2*ig+m)*cos_phi*sin_mphi - m*sin_phi*cos_mphi )*GenGrad_s;
              By += term*( (2*ig+m)*sin_phi*sin_mphi + m*cos_phi*cos_mphi )*GenGrad_s;
            }
          }
          Bx *= bgg->strength;
          By *= bgg->strength;
          Bz = bgg->strength*(dAy_dx - dAx_dy);
          if (bgg->synchRad) {
            double deltaTemp, B2, F;
            /* This is only valid for ultra-relatistic particles that radiate a small fraction of its energy in any step */
            /* It only uses Bx, By and ignores opening angle effects ~1/\gamma */
            B2 = sqr(Bx) + sqr(By);
            deltaTemp = delta - radCoef*pCentral*(1.0+delta)*B2*ds;
            F = isrCoef*pCentral*(1.0 + delta)*sqrt(ds)*pow(B2, 3./4.);
            if (bgg->isr && np!=1)
              deltaTemp += F*gauss_rn_lim(0.0, 1.0, srGaussianLimit, random_2);
            if (sigmaDelta2)
              *sigmaDelta2 += sqr(F);
            delta = deltaTemp;
          }
        }

        if (iz!=(bggData->nz))
          s += step*bgg->zInterval;

#ifdef DEBUG
        fprintf(fpdebug, "%le %le %le %le %le %le %le %le %i\n", 
                iz*bggData->dz, x, y, px, py, s, dAz_dy, dAz_dx, iImpLoop);
#endif

      }
      /* Drift particles back 1/2 step to end on last magnet grid point used */
      px = ux;  /* Assuming Ax -> 0 */
      py = uy;  /* Assuming Ay -> 0 */
      denom = (1.0 + delta)*(1.0 + delta) - px*px - py*py;
      denom = 1.0/sqrt(denom);
      x -= 0.5*step*bgg->zInterval*px*denom;
      y -= 0.5*step*bgg->zInterval*py*denom;
      s -= 0.5*step*bgg->zInterval*(1.0 + delta)*denom;
      denom = (1.0 + delta)*(1.0 + delta) - px*px - py*py;
      denom = 1.0/sqrt(denom);
      if (iz<bggData->nz) {
        /* Drift forward */
	x += bggData->dz*(bggData->nz-1-(iz-1))*px*denom;
	y += bggData->dz*(bggData->nz-1-(iz-1))*py*denom;
	s += bggData->dz*(bggData->nz-1-(iz-1))*(1.0 + delta)*denom;
      }
 
#if !USE_MPI
      if (bgg->SDDSpo) {
        if (!SDDS_SetRowValues(bgg->SDDSpo, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, irow++,
                               bgg->poIndex[0], x,
                               bgg->poIndex[1], px*pCentral,
                               bgg->poIndex[2], y,
                               bgg->poIndex[3], py*pCentral,
                               bgg->poIndex[4], (bggData->nz-1)*bggData->dz,
                               bgg->poIndex[5], sqrt(sqr(pCentral*(1+delta))-(sqr(px)+sqr(py))*sqr(pCentral)),
                               bgg->poIndex[6], Bx,
                               bgg->poIndex[7], By,
                               bgg->poIndex[8], Bz,
                               -1) ||
            !SDDS_WritePage(bgg->SDDSpo)) {
          SDDS_SetError("Problem setting particle output data for BGGEXP");
          SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
        }
      }
#endif

      /* Drift backward from end of field map to exit plane */
      denom = (1.0 + delta)*(1.0 + delta) - px*px - py*py;
      denom = 1.0/sqrt(denom);
      x -= (zMax - zExit)*px*denom;
      y -= (zMax - zExit)*py*denom;
      s -= (zMax - zExit)*(1.0 + delta)*denom;
      x -= xExit + bgg->dxExpansion;  /* Shift x according to exit point & map shift */

      /* Compute transformation from magnet grid -> accelerator coordinates if bending magnet */
      if (bgg->isBend) {
        phi = atan2(xVertex-xExit, zExit-zVertex);  /* angle of outgoing x-plane w.r.t. magnet x */
        sin_phi = sin(phi);
        cos_phi = cos(phi);
      } else {
        sin_phi = 0;
        cos_phi = 1;
      }
      denom = (1.0 + delta)*(1.0 + delta) - px*px - py*py;
      denom = 1.0/sqrt(denom);
      magnet_s = x*sin_phi/( cos_phi - sin_phi*px*denom );
      
      x = (x + px*magnet_s*denom*denom)/cos_phi;
      y = (y + py*magnet_s*denom*denom);
      px = px*cos_phi + sin_phi/denom;
      delta_s += (1 + delta)*magnet_s*denom;

      /* Update (non-canonical) particle coordinates */
      part[ip][0] = x;
      denom = (1.0 + delta)*(1.0 + delta) - px*px - py*py;
      denom = 1.0/sqrt(denom);
      part[ip][1] = px*denom;
      part[ip][2] = y;
      part[ip][3] = py*denom;
      part[ip][4] = s + delta_s;
      part[ip][5] = delta;
    }
  }  else { 
    /* Non-symplectic */
    double B[3], p[3], Bphi, Br, B2Max, pErr[3];
    double pOrig;
    double xpTemp, ypTemp, xpNew, ypNew;
    double xTemp, yTemp, xNew, yNew, preFactorDz, deltaTemp;
    double magnet_s;

    for (ip=0; ip<np; ip++) {
      B2Max = 0;
      x = part[ip][0];
      xp = part[ip][1];
      y = part[ip][2];
      yp = part[ip][3];
      s = part[ip][4];
      delta = part[ip][5];

      /* Compute transformation from accelerator coordinates */
      if (bgg->isBend) {
        phi = atan2(xVertex-xEntry, zVertex-zEntry);  /* angle of incoming x-plane w.r.t. magnet x */
        magnet_s = x*sin(phi)/( cos(phi) - xp*sin(phi) );
        x = (x + xp*magnet_s)/cos(phi);
        y =  y + yp*magnet_s;
        s += sqrt(1.0 + xp*xp + yp*yp)*magnet_s;
        yp = yp/(cos(phi) - xp*sin(phi));
        xp = (xp*cos(phi) + sin(phi))/(cos(phi) - xp*sin(phi));
      } else {
        phi = magnet_s = 0;
      }

      x += xEntry + bgg->dxExpansion; /* Shift x according to entrance point & map shift */

      /* Drift backward from entrance plane to beginning of field map */
      x -= (zEntry - zMin)*xp;
      y -= (zEntry - zMin)*yp;
      s -= (zEntry - zMin)*sqrt(1.0 + xp*xp + yp*yp);

      /* compute momenta (x, y, z) */
      denom = sqrt(1 + sqr(xp) + sqr(yp));
      pOrig = pCentral*(1+delta);
      p[2] = pOrig/denom;
      p[0] = xp*p[2];
      p[1] = yp*p[2];
      gamma = sqrt(sqr(p[0]) + sqr(p[1]) + sqr(p[2]) + 1);
      pErr[0] = pErr[1] = pErr[2] = 0;

#if !USE_MPI
      if (bgg->SDDSpo) {
        if (!SDDS_StartPage(bgg->SDDSpo, bggData->nz+1) ||
            !SDDS_SetParameters(bgg->SDDSpo, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,
                                "particleID", (long)(part[ip][6]), "pCentral", pCentral, NULL)) {
          SDDS_SetError("Problem setting up particle output page for BGGEXP");
          SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
        }
      }
#endif

      /* Integrate through the magnet */
      B[0] = B[1] = B[2] = 0;
      for (iz=irow=0; iz<bggData->nz-1; iz+=bgg->zInterval) {
	denom = sqrt(1 + sqr(xp) + sqr(yp));
	p[2] = pCentral*(1+delta)/denom;
	p[0] = xp*p[2];
	p[1] = yp*p[2];
	gamma = sqrt(sqr(p[0]) + sqr(p[1]) + sqr(p[2]) + 1);
#if !USE_MPI
        if (bgg->SDDSpo &&
            !SDDS_SetRowValues(bgg->SDDSpo, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, irow++,
                               bgg->poIndex[0], x,
                               bgg->poIndex[1], p[0],
                               bgg->poIndex[2], y,
                               bgg->poIndex[3], p[1],
                               bgg->poIndex[4], iz*bggData->dz,
                               bgg->poIndex[5], p[2],
                               bgg->poIndex[6], B[0],
                               bgg->poIndex[7], B[1],
                               bgg->poIndex[8], B[2],
                               -1)) {
          SDDS_SetError("Problem setting particle output data for BGGEXP");
          SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
        }
#endif
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
        B[0] = (bgg->Bx + (Br*cos(phi) - Bphi*sin(phi)))*bgg->strength;
        B[1] = (bgg->By + (Br*sin(phi) + Bphi*cos(phi)))*bgg->strength;
        B[2] *= bgg->strength;
        
	dz = bggData->dz*bgg->zInterval;
	preFactorDz = -dz*particleCharge*particleRelSign/(pCentral*particleMass*c_mks*(1.0+delta));
	preFactorDz =  preFactorDz*sqrt(1.0 + xp*xp + yp*yp);
	/* Apply prediction step */
	xTemp = x + dz*xp;
	yTemp = y + dz*yp;
	xpTemp = xp + preFactorDz*( (yp*B[2] - (1.0+xp*xp)*B[1]) + xp*yp*B[0] );
	ypTemp = yp + preFactorDz*( ((1.0+yp*yp)*B[0] - xp*B[2]) - xp*yp*B[1] );
	ds = dz*sqrt(1+sqr(xp)+sqr(yp));

        r = sqrt(sqr(xTemp)+sqr(yTemp));
        phi = atan2(yTemp, xTemp);

        /* Compute fields at next z location */
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
            B[2] += term*bggData->dCmns_dz[im][ig][iz+1]*r*sin_mphi;
            term *= bggData->Cmns[im][ig][iz+1];
            Br   += term*(2*ig+m)*sin_mphi;
            Bphi += m*term*cos_mphi;
          }
        }
        B[0] = (Br*cos(phi) - Bphi*sin(phi))*bgg->strength;
        B[1] = (Br*sin(phi) + Bphi*cos(phi))*bgg->strength;
        B[2] *= bgg->strength;

	dz = bggData->dz*bgg->zInterval;
	preFactorDz = -dz*particleCharge*particleRelSign/(pCentral*particleMass*c_mks*(1.0+delta));
	preFactorDz =  preFactorDz*sqrt(1.0 + xpTemp*xpTemp + ypTemp*ypTemp);
	/* Apply correction step */
	xNew = 0.5*(x + xTemp + dz*xpTemp);
	yNew = 0.5*(y + yTemp + dz*ypTemp);
	xpNew = 0.5*(xp + xpTemp + preFactorDz*( (ypTemp*B[2] - (1.0+xpTemp*xpTemp)*B[1]) + xpTemp*ypTemp*B[0] ));
	ypNew = 0.5*(yp + ypTemp + preFactorDz*( ((1.0+ypTemp*ypTemp)*B[0] - xpTemp*B[2]) - xpTemp*ypTemp*B[1] ));
	ds = 0.5*( ds + dz*sqrt(1+sqr(xpTemp)+sqr(ypTemp)) );

	x = xNew;
	y = yNew;
	xp = xpNew;
	yp = ypNew;
	s += ds;

	denom = sqrt(1 + sqr(xp) + sqr(yp));
	p[2] = pCentral*(1+delta)/denom;
	p[0] = xp*p[2];
	p[1] = yp*p[2];
        
#ifdef DEBUG
        fprintf(fpdebug, "%.0f %le %le %le %le %le %le %le %le %le %le\n", 
                part[ip][6], ds, x, y, iz*bggData->dz, 
                B[0], B[1], B[2], 
                p[0], p[1], p[2]);
#endif

        if (bgg->synchRad) {
          /* This is only valid for ultra-relatistic particles */
          double B2, F;
          B2 = sqr(B[0])+sqr(B[1]);
          if (B2>B2Max)
            B2Max = B2;
	  deltaTemp = delta - radCoef*pCentral*(1.0+delta)*B2*ds;
          F = isrCoef*pCentral*(1.0 + delta)*sqrt(ds)*pow(B2, 3./4.);
          if (bgg->isr && np!=1)
            deltaTemp += F*gauss_rn_lim(0.0, 1.0, srGaussianLimit, random_2);
          if (sigmaDelta2)
            *sigmaDelta2 += sqr(F);
	  delta = deltaTemp;
        }
      }
      if (iz<bggData->nz-1) {
        /* Drift forward */
        x += xp*bggData->dz*(bggData->nz-(iz-1));
        y += yp*bggData->dz*(bggData->nz-(iz-1));
        s += bggData->dz*(bggData->nz-(iz-1))*sqrt(1+sqr(xp)+sqr(yp));
      }
      
#if !USE_MPI
      if (bgg->SDDSpo) {
        if (!SDDS_SetRowValues(bgg->SDDSpo, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, irow++,
                               bgg->poIndex[0], x,
                               bgg->poIndex[1], p[0],
                               bgg->poIndex[2], y,
                               bgg->poIndex[3], p[1],
                               bgg->poIndex[4], (bggData->nz-1)*bggData->dz,
                               bgg->poIndex[5], p[2],
                               bgg->poIndex[6], B[0],
                               bgg->poIndex[7], B[1],
                               bgg->poIndex[8], B[2],
                               -1) ||
            !SDDS_WritePage(bgg->SDDSpo)) {
          SDDS_SetError("Problem setting particle output data for BGGEXP");
          SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
        }
      }
#endif

      /* Drift backward from end of field map to exit plane */
      x -= (zMax - zExit)*xp;
      y -= (zMax - zExit)*yp;
      s -= (zMax - zExit)*sqrt(1.0 + xp*xp + yp*yp);
      x -= xExit + bgg->dxExpansion;  /* Shift x according to exit point & map shift */
     
      /* Compute transformation from magnet grid -> accelerator coordinates if bending magnet */
      if (bgg->isBend) {
        phi = atan2(xVertex-xExit, zExit-zVertex);  /* angle of outgoing x-plane w.r.t. magnet x */
        magnet_s = x*sin(phi)/( cos(phi) - xp*sin(phi) );
        x = (x + xp*magnet_s)/cos(phi);
        y =  y + yp*magnet_s;
        s += sqrt(1.0 + xp*xp + yp*yp)*magnet_s;
        yp = yp/(cos(phi) - xp*sin(phi));
        xp = (xp*cos(phi) + sin(phi))/(cos(phi) - xp*sin(phi));
      } else {
        phi = magnet_s = 0;
      }

      part[ip][0] = x;
      part[ip][1] = xp; // p[0]/p[2];
      part[ip][2] = y;
      part[ip][3] = yp; // p[1]/p[2];
      part[ip][4] = s;
      /*
      printf("P: %le -> %le, change = %le, B2Max = %le\n",
             pOrig,
             sqrt(sqr(p[0])+sqr(p[1])+sqr(p[2])),
             sqrt(sqr(p[0])+sqr(p[1])+sqr(p[2]))-pOrig,
             B2Max);
      */
      if (bgg->synchRad) {
        double gamma0, beta0, gamma1, p1, beta1;
        gamma0 = sqrt(sqr(pCentral*(1+part[ip][5]))+1);
        beta0 = pCentral*(1+part[ip][5])/gamma0;
        p1 = sqrt(sqr(p[0])+sqr(p[1])+sqr(p[2]));
        gamma1 = sqrt(p1*p1+1);
        beta1 = p1/gamma1;
        part[ip][4] *= beta1/beta0;
        part[ip][5] = p1/pCentral - 1;
      }
    }
  }

  if (sigmaDelta2 && np)
    *sigmaDelta2 /= np;

  /* Do misalignments */
  if (bgg->tilt)
    rotateBeamCoordinates(part, np, -bgg->tilt);
  if (bgg->dx || bgg->dy || bgg->dz)
    offsetBeamCoordinates(part, np, -bgg->dx, -bgg->dy, -bgg->dz);

#ifdef DEBUG  
  fflush(fpdebug);
#endif

  return np;
}

