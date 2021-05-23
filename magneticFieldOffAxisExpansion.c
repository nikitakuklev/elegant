/* Copyright 2018 by Michael Borland, Ryan Lindberg, and Argonne National Laboratory,
 * all rights reserved.
 */
#include "mdb.h"
#include "SDDS.h"
#include "track.h"

typedef struct {
  long nz;             /* number of z points */
  double dz;           /* z spacing */
  double zMin, zMax;   /* minimum and maximum z values */
  long order;          /* order of derivative. 1 = quadrupole */
  long nDz1;           /* (number of z derivatives)+1 */
  double **fieldData;  /* data: Dz^m[Dx^n B] in T/m^(n+m), n=order */
} STORED_BOFFAXE_DATA;

static STORED_BOFFAXE_DATA *storedBOFFAXEData = NULL;
static long nBOFFAXEDataSets = 0;
static htab *fileHashTable = NULL;

void computeMagneticFieldFromOffAxisExpansion(double *B, double x, double y, long iz, long is, 
                                              BOFFAXE *boa, STORED_BOFFAXE_DATA *boaData);

#define BUFSIZE 16834

long addBOFFAXEData(char *filename, char *zColumn, char *dataColumn, long order)
{
  SDDS_DATASET SDDSin;
  TRACKING_CONTEXT tcontext;
  char buffer[BUFSIZE];
  char units1[BUFSIZE], units2[BUFSIZE];
  long readCode;
  long *nstore;
  
  if (!fileHashTable)
    fileHashTable = hcreate(12);

  if (hfind(fileHashTable, filename, strlen(filename))) {
    long istore;
    printf("Using previously-stored BOFFAXE data for filename %s\n", filename);
    fflush(stdout);
    istore = *((long*)hstuff(fileHashTable));
    printf("Using BOFFAXE table %ld for filename %s\n", istore, filename);
    fflush(stdout);
    return istore;
  }

  getTrackingContext(&tcontext);
  printf("Adding BOFFAXE data from file %s for element %s #%ld\n", filename, tcontext.elementName, tcontext.elementOccurrence);
  fflush(stdout);
  
  if (!SDDS_InitializeInputFromSearchPath(&SDDSin, filename))
    bombElegantVA("Unable to read file %s for BOFFAXE %s #%ld\n", filename, tcontext.elementName, tcontext.elementOccurrence);

  /* Check presence of z column */
  if (SDDS_CheckColumn(&SDDSin, zColumn, "m", SDDS_ANY_FLOATING_TYPE, stderr)!=SDDS_CHECK_OK) 
    bombElegantVA("Unable to find floating-point column \"%s\" with units \"m\" in file %s for BOFFAXE %s#%ld\n", 
                  zColumn, filename, tcontext.elementName, tcontext.elementOccurrence);

  /* Check presence of field column */
  if (order>1) {
    sprintf(units1, "T/m^%ld", order);
    sprintf(units2, "T/m$a%ld$n", order);
  } else {
    strcpy(units1, "T/m");
    strcpy(units2, "T/m");
  }

  if (SDDS_CheckColumn(&SDDSin, dataColumn, units1, SDDS_ANY_FLOATING_TYPE, NULL)!=SDDS_CHECK_OK &&
      SDDS_CheckColumn(&SDDSin, dataColumn, units2, SDDS_ANY_FLOATING_TYPE, NULL)!=SDDS_CHECK_OK) {
    if (order==1)
      bombElegantVA("Unable to find floating-point column \"%s\" with units \"%s\" in file %s for BOFFAXE %s#%ld\n", 
                    dataColumn, units1, filename, tcontext.elementName, tcontext.elementOccurrence);
    else 
      bombElegantVA("Unable to find floating-point column \"%s\" with units \"%s\" or \"%s\" in file %s for BOFFAXE %s#%ld\n", 
                    dataColumn, units1, units2, filename, tcontext.elementName, tcontext.elementOccurrence);
  }

  if (!(storedBOFFAXEData = SDDS_Realloc(storedBOFFAXEData, sizeof(*storedBOFFAXEData)*(nBOFFAXEDataSets+1))))
    bombElegantVA("Memory allocation error reading data from file %s for BOFFAXE %s #%ld\n", filename, tcontext.elementName, tcontext.elementOccurrence);

  storedBOFFAXEData[nBOFFAXEDataSets].nz = 0;
  storedBOFFAXEData[nBOFFAXEDataSets].dz = 0;
  storedBOFFAXEData[nBOFFAXEDataSets].zMin = DBL_MAX;
  storedBOFFAXEData[nBOFFAXEDataSets].zMax = -DBL_MAX;
  storedBOFFAXEData[nBOFFAXEDataSets].order = 0;
  storedBOFFAXEData[nBOFFAXEDataSets].fieldData = NULL;
  storedBOFFAXEData[nBOFFAXEDataSets].nDz1 = 1;
  
  /* find z-derivative columns, with names of form <dataColumn>Deriv, <dataColumn>Deriv2, <dataColumn>Deriv3,  ... */
  while (1) {
    if (storedBOFFAXEData[nBOFFAXEDataSets].nDz1==1)
      sprintf(buffer, "%sDeriv", dataColumn);
    else 
      sprintf(buffer, "%sDeriv%ld", dataColumn, storedBOFFAXEData[nBOFFAXEDataSets].nDz1);
    if (SDDS_CheckColumn(&SDDSin, buffer, NULL, SDDS_ANY_FLOATING_TYPE, NULL)==SDDS_CHECK_OK) {
      printf("Found %ldth derivative in column %s\n", storedBOFFAXEData[nBOFFAXEDataSets].nDz1, buffer);
      fflush(stdout);
      storedBOFFAXEData[nBOFFAXEDataSets].nDz1++;
    } else 
      break;
  } 
  storedBOFFAXEData[nBOFFAXEDataSets].fieldData 
    = tmalloc(sizeof(*(storedBOFFAXEData[nBOFFAXEDataSets].fieldData))*storedBOFFAXEData[nBOFFAXEDataSets].nDz1);

  while ((readCode=SDDS_ReadPage(&SDDSin))>0) {
    if (readCode==1) {
      long iz;
      double dz0, dz, *z, zMin, zMax;
      if ((storedBOFFAXEData[nBOFFAXEDataSets].nz = SDDS_RowCount(&SDDSin))<=1)
        bombElegantVA("Too few z values in file %s for BOFFAXE %s#%ld\n", 
                      filename, tcontext.elementName, tcontext.elementOccurrence);
      if (!(z=SDDS_GetColumnInDoubles(&SDDSin, zColumn)))
        bombElegantVA("Problem reading z column \"%s\" from file %s for BOFFAXE %s#%ld\n", 
                      zColumn, filename, tcontext.elementName, tcontext.elementOccurrence);
      find_min_max(&zMin, &zMax, z, storedBOFFAXEData[nBOFFAXEDataSets].nz);
      if (zMin<storedBOFFAXEData[nBOFFAXEDataSets].zMin)
        storedBOFFAXEData[nBOFFAXEDataSets].zMin = zMin;
      if (zMax>storedBOFFAXEData[nBOFFAXEDataSets].zMax)
        storedBOFFAXEData[nBOFFAXEDataSets].zMax = zMax;
      dz0 = z[1] - z[0];
      for (iz=1; iz<storedBOFFAXEData[nBOFFAXEDataSets].nz; iz++) {
        dz = z[iz] - z[iz-1];
        if (dz<=0)
          bombElegantVA("Data not monotonically increasing in z column from %s for BOFFAXE %s #%ld\n",
                        buffer, filename, tcontext.elementName, tcontext.elementOccurrence);
        if (fabs(dz0/dz-1)>1e-6)
          bombElegantVA("Data not uniformly increasing in z column from %s for BOFFAXE %s #%ld (%le vs %le)\n",
                        buffer, filename, tcontext.elementName, tcontext.elementOccurrence,
                        dz0, dz);
      }
      free(z);
      storedBOFFAXEData[nBOFFAXEDataSets].dz = dz0;
      if (!(storedBOFFAXEData[nBOFFAXEDataSets].fieldData[0] = SDDS_GetColumnInDoubles(&SDDSin, dataColumn)))
        bombElegantVA("Problem reading column \"%s\" from file %s for BOFFAXE %s#%ld\n", 
                      dataColumn, filename, tcontext.elementName, tcontext.elementOccurrence);
      for (iz=1; iz<storedBOFFAXEData[nBOFFAXEDataSets].nDz1; iz++) {
        if (iz==1)
          sprintf(buffer, "%sDeriv", dataColumn);
        else 
          sprintf(buffer, "%sDeriv%ld", dataColumn, iz);
      if (!(storedBOFFAXEData[nBOFFAXEDataSets].fieldData[iz] = SDDS_GetColumnInDoubles(&SDDSin, buffer)))
        bombElegantVA("Problem reading column \"%s\" from file %s for BOFFAXE %s#%ld\n", 
                      buffer, filename, tcontext.elementName, tcontext.elementOccurrence);
      }
    } else {
      bombElegantVA("Input file %s for BOFFAXE %s#%ld has more than one page", 
                    filename, tcontext.elementName, tcontext.elementOccurrence);
    }
  }
  SDDS_Terminate(&SDDSin);
  
  nstore = tmalloc(sizeof(*nstore));
  *nstore = nBOFFAXEDataSets;
  hadd(fileHashTable, filename, strlen(filename), (void*)nstore);

  printf("Done adding BOFFAXE data from file %s\n", filename);
  fflush(stdout);

  return nBOFFAXEDataSets++;
}

long trackMagneticFieldOffAxisExpansion(double **part, long np, BOFFAXE *boa, double pCentral, double **accepted, double *sigmaDelta2)
{
  long ip, iz, irow;
  STORED_BOFFAXE_DATA *boaData;
  double ds, dz, x, y, xp, yp, delta, s, denom;
  double length;
  TRACKING_CONTEXT tcontext;
  double radCoef=0, isrCoef=0;
  double zMin, zMax, zEntry, zExit;

  double B[3], p[3], B2Max, pErr[3];
  double pOrig;
  double xpTemp, ypTemp, xpNew, ypNew;
  double xTemp, yTemp, xNew, yNew, preFactorDz, deltaTemp;

#ifdef DEBUG
  static FILE *fpdebug = NULL;
  if (!fpdebug) {
    fpdebug = fopen("boaexp.deb", "w");
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
  if (boa->synchRad) {
    radCoef = ipow(particleCharge,4)/(6*PI*epsilon_o*ipow(c_mks,4)*ipow(particleMass,3));
    isrCoef = sqrt(55/(24*sqrt(3))*particleRadius*hbar_mks*ipow(particleCharge, 3))/sqr(particleMass*c_mks);
  }
  if (sigmaDelta2)
    *sigmaDelta2 = 0;

  getTrackingContext(&tcontext);

  if (!boa->initialized) {
    /* char *outputFile; */
    boa->initialized = 1;
    if (!(boa->filename) || !strlen(boa->filename))
      bombElegantVA("No filename given for BOFFAXE %s#%ld\n", tcontext.elementName, tcontext.elementOccurrence);
    if (!(boa->zColumn) || !strlen(boa->zColumn))
      bombElegantVA("No Z_COLUMN given for BOFFAXE %s#%ld\n", tcontext.elementName, tcontext.elementOccurrence);
    if (!(boa->fieldColumn) || !strlen(boa->fieldColumn))
      bombElegantVA("No FIELD_COLUMN given for BOFFAXE %s#%ld\n", tcontext.elementName, tcontext.elementOccurrence);
    if (boa->order<1)
      bombElegantVA("ORDER<1 given for BOFFAXE %s#%ld\n", tcontext.elementName, tcontext.elementOccurrence);
    if (boa->zInterval<1)
      bombElegantVA("Z_INTERVAL<1 given for BOFFAXE %s#%ld\n", tcontext.elementName, tcontext.elementOccurrence);
    boa->dataIndex = addBOFFAXEData(boa->filename, boa->zColumn, boa->fieldColumn, boa->order);
#if !USE_MPI
    if (boa->particleOutputFile && !boa->SDDSpo) {
      boa->SDDSpo = tmalloc(sizeof(*(boa->SDDSpo)));
      boa->particleOutputFile = compose_filename(boa->particleOutputFile, tcontext.rootname);
      if (!SDDS_InitializeOutputElegant(boa->SDDSpo, SDDS_BINARY, 1, 
                                 NULL, NULL, boa->particleOutputFile) ||
          0>SDDS_DefineParameter(boa->SDDSpo, "SVNVersion", NULL, NULL, "SVN version number", NULL, SDDS_STRING, SVN_VERSION) ||
          !SDDS_DefineSimpleParameter(boa->SDDSpo, "particleID", NULL, SDDS_LONG) ||
          !SDDS_DefineSimpleParameter(boa->SDDSpo, "pCentral", "m$be$nc", SDDS_DOUBLE) ||
          (boa->poIndex[0]=SDDS_DefineColumn(boa->SDDSpo, "x", NULL, "m", NULL, NULL, SDDS_DOUBLE, 0 ))<0 ||
          (boa->poIndex[1]=SDDS_DefineColumn(boa->SDDSpo, "px", NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0))<0 ||
          (boa->poIndex[2]=SDDS_DefineColumn(boa->SDDSpo, "y", NULL, "m", NULL, NULL, SDDS_DOUBLE, 0 ))<0 ||
          (boa->poIndex[3]=SDDS_DefineColumn(boa->SDDSpo, "py", NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0))<0 ||
          (boa->poIndex[4]=SDDS_DefineColumn(boa->SDDSpo, "z", NULL, "m", NULL, NULL, SDDS_DOUBLE, 0))<0 ||
          (boa->poIndex[5]=SDDS_DefineColumn(boa->SDDSpo, "pz", NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0))<0 ||
          (boa->poIndex[6]=SDDS_DefineColumn(boa->SDDSpo, "Bx", NULL, "T", NULL, NULL, SDDS_DOUBLE, 0))<0 ||
          (boa->poIndex[7]=SDDS_DefineColumn(boa->SDDSpo, "By", NULL, "T", NULL, NULL, SDDS_DOUBLE, 0))<0 ||
          (boa->poIndex[8]=SDDS_DefineColumn(boa->SDDSpo, "Bz", NULL, "T", NULL, NULL, SDDS_DOUBLE, 0))<0 ||
          !SDDS_WriteLayout(boa->SDDSpo)) {
        SDDS_SetError("Problem setting up particle output file for BOAEXP");
        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
      }
    }
#endif
  }
  boaData = storedBOFFAXEData+boa->dataIndex;
  if (boaData->nz%boa->zInterval!=0) 
      bombElegantVA("Z_INTERVAL (%ld) given for BOFFAXE %s#%ld doesn't evenly divide the number of data points (%ld) in the input file %s\n", 
                    boa->zInterval, tcontext.elementName, tcontext.elementOccurrence, boaData->nz);

  if (boa->fieldLength>0) {
    if (fabs(boa->fieldLength-(boaData->zMax-boaData->zMin))>1e-6*boaData->dz) 
      bombElegantVA("FIELD_LENGTH value %21.15e for BOFFAXE %s#%ld does not match z range %21.15e in file %s\n",
                    boa->fieldLength, tcontext.elementName, tcontext.elementOccurrence, boa->filename);
  }
  if ((length = boa->length)<0)
    bombElegantVA("Negative LENGTH value %21.15e for BOFFAXE %s#%ld is not permitted\n",
                  boa->length, tcontext.elementName, tcontext.elementOccurrence);
  
  zMin = -(boaData->zMax-boaData->zMin)/2;
  zMax = +(boaData->zMax-boaData->zMin)/2;
  
  zEntry = -length/2;
  zExit = length/2;

  /* Do misalignments */
  if (boa->dx || boa->dy || boa->dz)
    offsetBeamCoordinatesForMisalignment(part, np, boa->dx, boa->dy, boa->dz);
  if (boa->tilt)
    rotateBeamCoordinatesForMisalignment(part, np, boa->tilt);

  /* Non-symplectic integrator by R. Lindberg, taken from BGGEXP code */
  for (ip=0; ip<np; ip++) {
#if !USE_MPI
    if (boa->SDDSpo) {
      if (!SDDS_StartPage(boa->SDDSpo, boaData->nz+1) ||
          !SDDS_SetParameters(boa->SDDSpo, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,
                              "particleID", (long)(part[ip][6]), "pCentral", pCentral, NULL)) {
        SDDS_SetError("Problem setting up particle output page for BOFFAXE");
        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
      }
    }
#endif
    B2Max = 0;
    x = part[ip][0];
    xp = part[ip][1];
    y = part[ip][2];
    yp = part[ip][3];
    s = part[ip][4];
    delta = part[ip][5];

    /* Drift from entrance plane to beginning of field map */
    x -= (zEntry - zMin)*xp;
    y -= (zEntry - zMin)*yp;
    s -= (zEntry - zMin)*sqrt(1.0 + xp*xp + yp*yp);

    /* compute momenta (x, y, z) */
    denom = sqrt(1 + sqr(xp) + sqr(yp));
    pOrig = pCentral*(1+delta);
    p[2] = pOrig/denom;
    p[0] = xp*p[2];
    p[1] = yp*p[2];
    /* gamma = sqrt(sqr(p[0]) + sqr(p[1]) + sqr(p[2]) + 1); */
    pErr[0] = pErr[1] = pErr[2] = 0;


    /* Integrate through the magnet */
    B[0] = B[1] = B[2] = 0;
    dz = boaData->dz*boa->zInterval/boa->zSubdivisions;
    for (iz=irow=0; iz<boaData->nz-1; iz+=boa->zInterval) {
      long is;
      for (is=0; is<boa->zSubdivisions; is++) {
        denom = sqrt(1 + sqr(xp) + sqr(yp));
        p[2] = pCentral*(1+delta)/denom;
        p[0] = xp*p[2];
        p[1] = yp*p[2];
        /* gamma = sqrt(sqr(p[0]) + sqr(p[1]) + sqr(p[2]) + 1); */
      
        /* Compute fields */
        computeMagneticFieldFromOffAxisExpansion(B, x, y, iz, is, boa, boaData);
        
#if !USE_MPI
        if (is==0 && boa->SDDSpo &&
            !SDDS_SetRowValues(boa->SDDSpo, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, irow++,
                               boa->poIndex[0], x,
                               boa->poIndex[1], p[0]*pCentral,
                               boa->poIndex[2], y,
                               boa->poIndex[3], p[1]*pCentral,
                               boa->poIndex[4], boaData->zMin+iz*dz*boa->zSubdivisions + is*dz,
                               boa->poIndex[5], sqrt(sqr(pCentral*(1+delta))-(sqr(p[0])+sqr(p[1]))*sqr(pCentral)),
                               boa->poIndex[6], B[0],
                               boa->poIndex[7], B[1],
                               boa->poIndex[8], B[2],
                               -1)) {
          SDDS_SetError("Problem setting particle output data for BGGEXP");
          SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
        }
#endif

        preFactorDz = -dz*particleCharge*particleRelSign/(pCentral*particleMass*c_mks*(1.0+delta));
        preFactorDz =  preFactorDz*sqrt(1.0 + xp*xp + yp*yp);

        /* Apply prediction step */
        xTemp = x + dz*xp;
        yTemp = y + dz*yp;
        xpTemp = xp + preFactorDz*( (yp*B[2] - (1.0+xp*xp)*B[1]) + xp*yp*B[0] );
        ypTemp = yp + preFactorDz*( ((1.0+yp*yp)*B[0] - xp*B[2]) - xp*yp*B[1] );
        ds = dz*sqrt(1+sqr(xp)+sqr(yp));
      
        /* Compute fields at next z location */
        if (is==(boa->zSubdivisions-1))
          computeMagneticFieldFromOffAxisExpansion(B, x, y, iz+boa->zInterval, 0, boa, boaData);
        else
          computeMagneticFieldFromOffAxisExpansion(B, x, y, iz, is+1, boa, boaData);
          
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
        
        if (boa->synchRad) {
          /* This is only valid for ultra-relatistic particles */
          double B2, F;
          B2 = sqr(B[0])+sqr(B[1]);
          if (B2>B2Max)
            B2Max = B2;
          deltaTemp = delta - radCoef*pCentral*(1.0+delta)*B2*ds;
          F = isrCoef*pCentral*(1.0 + delta)*sqrt(ds)*pow(B2, 3./4.);
          if (boa->isr && np!=1)
            deltaTemp += F*gauss_rn_lim(0.0, 1.0, srGaussianLimit, random_2);
          if (sigmaDelta2)
            *sigmaDelta2 += sqr(F);
          delta = deltaTemp;
        }
      }
    }
    if (iz<(boaData->nz-1)) {
      /* Drift forward */
      x += xp*boaData->dz*(boaData->nz-(iz-1));
      y += yp*boaData->dz*(boaData->nz-(iz-1));
      s += boaData->dz*(boaData->nz-(iz-1))*sqrt(1+sqr(xp)+sqr(yp));
    }

    /* Drift backward from end of field map to exit plane */
    x -= (zMax - zExit)*xp;
    y -= (zMax - zExit)*yp;
    s -= (zMax - zExit)*sqrt(1.0 + xp*xp + yp*yp);
    
    part[ip][0] = x;
    part[ip][1] = xp; /*  p[0]/p[2]; */
    part[ip][2] = y;
    part[ip][3] = yp; /*  p[1]/p[2]; */
    part[ip][4] = s;

    if (boa->synchRad) {
      double gamma0, beta0, gamma1, p1, beta1;
      gamma0 = sqrt(sqr(pCentral*(1+part[ip][5]))+1);
      beta0 = pCentral*(1+part[ip][5])/gamma0;
      p1 = sqrt(sqr(p[0])+sqr(p[1])+sqr(p[2]));
      gamma1 = sqrt(p1*p1+1);
      beta1 = p1/gamma1;
      part[ip][4] *= beta1/beta0;
      part[ip][5] = p1/pCentral - 1;
    }

#if !USE_MPI
    if (boa->SDDSpo) {
      if (!SDDS_SetRowValues(boa->SDDSpo, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, irow++,
                             boa->poIndex[0], x,
                             boa->poIndex[1], p[0],
                             boa->poIndex[2], y,
                             boa->poIndex[3], p[1],
                             boa->poIndex[4], boaData->zMin+iz*dz*boa->zSubdivisions,
                             boa->poIndex[5], p[2],
                             boa->poIndex[6], B[0],
                             boa->poIndex[7], B[1],
                             boa->poIndex[8], B[2],
                             -1) ||
          !SDDS_WritePage(boa->SDDSpo)) {
        SDDS_SetError("Problem setting particle output data for BGGEXP");
        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
      }
    }
#endif
  }

  if (sigmaDelta2 && np)
    *sigmaDelta2 /= np;

  /* Do misalignments */
  if (boa->tilt)
    rotateBeamCoordinatesForMisalignment(part, np, -boa->tilt);
  if (boa->dx || boa->dy || boa->dz)
    offsetBeamCoordinatesForMisalignment(part, np, -boa->dx, -boa->dy, -boa->dz);

#ifdef DEBUG  
  fflush(fpdebug);
#endif

  return np;
}

void computeMagneticFieldFromOffAxisExpansion
(
 double *B, /* Bx, By, Bz output */
 double x,
 double y,
 long iz,
 long is,
 BOFFAXE *boa,
 STORED_BOFFAXE_DATA *boaData
 ) {
  long i;
  double delz = 0;
  double *fieldData;

  if (boa->order!=1 && boa->order!=2)
    bombElegant("Only ORDER=1 and ORDER=2 implemented for BOFFAXE", NULL);
  if (iz<0 || iz>=boaData->nz)
    B[0] = B[1] = B[2] = 0;
  else {
    fieldData = tmalloc(sizeof(*fieldData)*boaData->nDz1);
    if (is && boa->zSubdivisions>1) 
      delz = (1.0*is)/boa->zSubdivisions*boaData->dz;
    for (i=0; i<boaData->nDz1; i++) {
      long j;
      fieldData[i] = boaData->fieldData[i][iz];
      if (delz)
        for (j=i+1; j<boaData->nDz1; j++)
          fieldData[i] += boaData->fieldData[j][iz]*ipow(delz, j-i)/dfactorial(j-i);
    }
    if (boa->order==1) {
      /* quadrupole */
      B[0] = fieldData[0]*y;
      B[1] = fieldData[0]*x;
      B[2] = 0;
      /* Probably the dumbest possible way to do this... */
      if (boaData->nDz1>1 && (boa->expansionOrder==0 || boa->expansionOrder>1)) {
        B[2] = fieldData[1]*x*y;
        if (boaData->nDz1>2 && (boa->expansionOrder==0 || boa->expansionOrder>2)) {
          B[0] -= fieldData[2]*(x*x*y/4     + ipow(y,3)/12);
          B[1] -= fieldData[2]*(ipow(x,3)/12 + x*y*y/4);
          if (boaData->nDz1>3 && (boa->expansionOrder==0 || boa->expansionOrder>3)) {
            B[2] -= fieldData[3]*(ipow(x, 3)*y + x*ipow(y, 3))/12;
            if (boaData->nDz1>4 && (boa->expansionOrder==0 || boa->expansionOrder>4)) {
              B[0] += fieldData[4]*(ipow(x, 2)*ipow(y, 3)/36 + (5*ipow(x, 4)*y + ipow(y, 5))/720);
              B[1] += fieldData[4]*(ipow(x, 3)*ipow(y, 2)/36 + (ipow(x, 5) + 5*x*ipow(y, 4))/720);
              if (boaData->nDz1>5 && (boa->expansionOrder==0 || boa->expansionOrder>5)) {
                B[2] += fieldData[5]*(ipow(x*y, 3)/108 + (ipow(x, 5)*y + x*ipow(y, 5))/720);
                if (boaData->nDz1>6 && (boa->expansionOrder==0 || boa->expansionOrder>6)) {
                  B[0] -= fieldData[6]*(5*ipow(x, 4)*ipow(y, 3) + 3*ipow(x, 2)*ipow(y, 5))/4320;
                  B[1] -= fieldData[6]*(3*ipow(x, 5)*ipow(y, 2) + 5*ipow(x, 3)*ipow(y, 4))/4320;
                  if (boaData->nDz1>7 && (boa->expansionOrder==0 || boa->expansionOrder>7)) {
                    B[2] -= fieldData[7]*(ipow(x, 5)*ipow(y, 3) + ipow(x, 3)*ipow(y, 5))/4320;
                  }
                }
              }
            }
          }
        }
      }
    } else {
      /* sextupole */
      B[0] = fieldData[0]*x*y;
      B[1] = 0.5*fieldData[0]*(x*x-y*y);
      B[2] = 0;
      /* Probably the dumbest possible way to do this... */
      if (boaData->nDz1>1 && (boa->expansionOrder==0 || boa->expansionOrder>1)) {
        B[2] = fieldData[1]*(0.5*sqr(x)*y - ipow(y,3)/6);
        if (boaData->nDz1>2 && (boa->expansionOrder==0 || boa->expansionOrder>2)) {
          B[0] -= fieldData[2]*(x*ipow(y,3)/24 + ipow(x,3)*y/8);
          B[1] += fieldData[2]*(5*ipow(y,4)/96- ipow(x*y,2)/16 - ipow(x,4)/32);
          if (boaData->nDz1>3 && (boa->expansionOrder==0 || boa->expansionOrder>3)) {
            B[2] += fieldData[3]*(ipow(y,5)/96 - sqr(x)*ipow(y,3)/48 - ipow(x,4)*y/32);
            if (boaData->nDz1>4 && (boa->expansionOrder==0 || boa->expansionOrder>4)) {
              B[0] += fieldData[4]*(x*ipow(y,5)/1920 + ipow(x*y,3)/192 + 3*ipow(x,5)*y/640);
              B[1] += fieldData[4]*(-7*ipow(y,6)/3840 + ipow(x*sqr(y),2)/768 + ipow(sqr(x)*y,2)/256 + ipow(x,6)/1280);
              if (boaData->nDz1>5 && (boa->expansionOrder==0 || boa->expansionOrder>5)) {
                B[2] += fieldData[5]*(-ipow(y,7)/3840 + sqr(x)*ipow(y,5)/3840 + ipow(x,4)*ipow(y,3)/768 + ipow(x,6)*y/1280);
                if (boaData->nDz1>6 && (boa->expansionOrder==0 || boa->expansionOrder>6)) {
                  B[0] += fieldData[6]*(-ipow(x,3)*ipow(y,5)/11520 - ipow(x,5)*ipow(y,3)/5760 - ipow(x,7)*y/11520);
                  B[1] += fieldData[6]*(ipow(y,8)/30720 - ipow(x,6)*sqr(y)/11520 - ipow(x*y,4)/9216 - ipow(x,8)/92160);
                  if (boaData->nDz1>7 && (boa->expansionOrder==0 || boa->expansionOrder>7)) {
                    B[2] += fieldData[7]*(ipow(y,9)/276480 - ipow(x,4)*ipow(y,5)/46080 - ipow(x,6)*ipow(y,3)/34560 
                                          - ipow(x,8)*y/92160);
                  }
                }
              }
            }
          }
        }
      }
    }
    free(fieldData);
  }
  for (i=0; i<3; i++)
    B[i] *= boa->strength;
}
