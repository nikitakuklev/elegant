/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/


#include "mdb.h"
#include "madto.h"
#include <ctype.h>

static double bnEnge[3];
static double rhoEnge, lengthEnge, angleEnge, KgEnge;
double engeProfile(double z)
{
  return 1/(1 + exp(bnEnge[0] + bnEnge[1]*z + bnEnge[2]*sqr(z)));
}
double engeProfileTheta(double theta)
{
  return engeProfile(rhoEnge*sin(theta));
}

double engeEdgeFactor(double z)
{
  double ep;
  ep = engeProfile(z);
  return ep*(1-ep);
}


double engeOptimizationFunction(double *b, long *invalid)
{
  static double F[3];
  double temp1, temp2;
  
  bnEnge[0] = b[0];
  bnEnge[1] = b[1];
  bnEnge[2] = b[2];
  *invalid = 0;
  
  /* field goes to 1 at center of arc */
  F[0] = 1 - engeProfile(-rhoEnge*tan(angleEnge/2));

  /* effective length constraint */
  if (!gaussianQuadrature(engeProfileTheta, -angleEnge/2, 0.0, 100, 1e-8, &temp1)) {
    fprintf(stderr, "GQ #1 failed\n");
    *invalid = 1;
    return 0;
  }
  if (!gaussianQuadrature(engeProfile, 0.0, 10*lengthEnge, 1000, 1e-8, &temp2)) {
    fprintf(stderr, "GQ #2 failed\n");
    *invalid = 1;
    return 0;
  }
  F[1] = (lengthEnge/2 - (temp1 + temp2))/(lengthEnge/2);
  
  /* edge integral constraint */
  if (!gaussianQuadrature(engeEdgeFactor, -lengthEnge/2, 10*lengthEnge, 1000, 1e-8, &F[2])) {
    fprintf(stderr, "GQ #3 failed\n");
    *invalid = 1;
    return 0;
  }
  F[2] = (KgEnge - F[2])/KgEnge;
  
  return sqr(F[0]) + sqr(F[1]) + sqr(F[2]);
}


long computeEngeCoefficients(double *engeCoef, double rho, double length, double gap, double fint)
{
  double b[3], db[3], bMin[3], bMax[3];
  double result;
  static double lastData[4] = {-1, -1, -1, -1};
  static double lastResult[3];

  if (rho==lastData[0] && length==lastData[1] && gap==lastData[2] && fint==lastData[3]) {
    memcpy(engeCoef, lastResult, 3*sizeof(*engeCoef));
    return 1;
  }
  
  rhoEnge = rho;
  lengthEnge = length;
  angleEnge = length/rho;
  KgEnge = fint*gap;
  b[0] = -0.003183;
  b[1] = 1.911302/gap;
  b[2] = 0.0;
  db[0] = db[1] = db[2] = 1e-4;
  bMin[0] = bMin[1] = bMin[2] = -1e10;
  bMax[0] = bMax[1] = bMax[2] =  1e10;

  if (simplexMin(&result, b, db, bMin, bMax, NULL, 3, 1e-14, 1e-16, engeOptimizationFunction,
                 NULL, 500, 3, 12, 10.0, 10.0, 0)<0) {
    fprintf(stderr, "Problem finding enge coefficients, using defaults\n");
    engeCoef[0] = -0.003183;
    engeCoef[1] = 1.911302;
    engeCoef[2] = 0.0;
  }
  engeCoef[0] = b[0];
  engeCoef[1] = b[1]*gap;
  engeCoef[2] = b[2]*sqr(gap);

  lastData[0] = rho;
  lastData[1] = length;
  lastData[2] = gap;
  lastData[3] = fint;
  memcpy(lastResult, engeCoef, 3*sizeof(*lastResult));
  return 1;
}

void emitCosyDrift(FILE *fp, char *name, double length)
{
  fprintf(fp, "  PROCEDURE %s ;\n", name);
  fprintf(fp, "    DL %.15g ;\n", length);
  fprintf(fp, "  ENDPROCEDURE ;\n");
}

void emitCosyDipole(FILE *fp, char *name,
                    double length, double angle, 
                    double e1, double e2,
                    double h1, double h2,
                    double K1, double K2, double K3,
                    double hgap, double fint)
{
  double rho, engeCoef[3];
  long i, j;
  double sign = 1;
  
  if (angle<0)
    sign = -1;
  rho = length/(fabs(angle)+1e-300),
  fprintf(fp, "  PROCEDURE %s ;\n", name);
  fprintf(fp, "    VARIABLE W 1 4 4 ;\n");
  fprintf(fp, "    W(1,1) := 1 ;\n");
  fprintf(fp, "    W(2,1) := %.15g ;\n", K1*rho);
  fprintf(fp, "    W(3,1) := %.15g ;\n", K2*rho/2);
  fprintf(fp, "    W(4,1) := %.15g ;\n", K3*rho/6);
  for (i=1; i<=4; i++) 
    for (j=2; j<=4; j++)
      fprintf(fp, "    W(%ld, %ld) := 0 ;\n", i, j);

  if (angle<0) 
    fprintf(fp, "    CB ;\n");
  if (hgap!=0 && fint!=0) {
    computeEngeCoefficients(engeCoef, rho, length, 2*hgap, fint);
    fprintf(fp, "    FR 3 ;\n");
    fprintf(fp, "    FC 1 1 1 %.15g %.15g %.15g\n       0.0 0.0 0.0 ; \n",
            engeCoef[0], engeCoef[1], engeCoef[2]);
    fprintf(fp, "    FC 1 2 1 %.15g %.15g %.15g\n       0.0 0.0 0.0 ; \n",
            engeCoef[0], engeCoef[1], engeCoef[2]);
  } else {
    fprintf(fp, "    FR 1 ;\n");
  }
  fprintf(fp, "    MSS %.15g %.15g %.15g %.15g\n       %.15g %.15g %.15g W; \n",
          rho, 180/PI*fabs(angle), 
          hgap==0?1e-3:hgap, 
          sign*180/PI*e1, h1, 
          sign*180/PI*e2, h2);
  if (angle<0) 
    fprintf(fp, "    CB ;\n");
  fprintf(fp, "    FD ; \n");
  fprintf(fp, "  ENDPROCEDURE; \n");
}


void convert_to_cosy(char *outputfile, LINE_LIST *beamline, 
                     long cosyOrder, double pCentral,
                     double quadBoreRadius, double sextBoreRadius)
{
    ELEMENT_LIST *eptr;
    QUAD  *quad; KQUAD *kquad;
    SEXT  *sext; KSEXT *ksext;
    BEND  *bend;
    HCOR  *hcor;
    VCOR  *vcor;
    DRIFT *drift;
    CSBEND *csbend;
    CSRCSBEND *csrbend;
    CSRDRIFT *csrdrift;
    RFCA *rfca;
    FILE *fp;
    double BRho;
    long i;
    char *ptr;
#define NVAR 6
    char *varName[NVAR] = {
      "A", "B", "G", "R", "MU", "F"
      };
    char **nameUsed = NULL;
    long namesUsed = 0;

    fp = fopen_e(outputfile, "w", 0);

    fprintf(fp, "INCLUDE 'COSY';\n");
    fprintf(fp, "PROCEDURE RUN ;\n");
    for (i=0; i<NVAR; i++)
      fprintf(fp, "    VARIABLE %s 100 2 ; \n", varName[i]);
      
    BRho = pCentral*me_mks*c_mks/e_mks;

    /* emit procedures describing elements */
    eptr = &(beamline->elem);
    while (eptr) {
      while (ptr=strchr(eptr->name, ':'))
        *ptr = '_';
      for (i=0; i<namesUsed; i++) {
        if (strcmp(eptr->name, nameUsed[i])==0)
          break;
      }
      if (i!=namesUsed) {
        eptr = eptr->succ;
        continue;
      }
      nameUsed = trealloc(nameUsed, sizeof(*nameUsed)*(namesUsed+1));
      nameUsed[namesUsed++] = eptr->name;
      switch (eptr->type) {
      case T_QUAD:
        quad = (QUAD*)eptr->p_elem;
        fprintf(fp, "  PROCEDURE %s ;\n", eptr->name);
        fprintf(fp, "    { K1 = %e, R = %e, BRho = %e }\n",
                quad->k1, quadBoreRadius, BRho);
        fprintf(fp, "    MQ %.15g %.15g %.15g ; \n",
                quad->length, -quad->k1*BRho*quadBoreRadius, quadBoreRadius);
        fprintf(fp, "  ENDPROCEDURE; \n");
        break;
      case T_KQUAD:
        kquad = (KQUAD*)eptr->p_elem;
        fprintf(fp, "  PROCEDURE %s ;\n", eptr->name);
        fprintf(fp, "    { K1 = %e, R = %e, BRho = %e }\n",
                kquad->k1, quadBoreRadius, BRho);
        fprintf(fp, "    MQ %.15g %.15g %.15g ; \n",
                kquad->length, -kquad->k1*BRho*quadBoreRadius, quadBoreRadius);
        fprintf(fp, "  ENDPROCEDURE; \n");
        break;
      case T_SBEN:
        bend = (BEND*)eptr->p_elem;
        emitCosyDipole(fp, eptr->name, bend->length, bend->angle,
                       bend->e1, bend->e2, bend->h1, bend->h2,
                       bend->k1, bend->k2, 0.0,
                       bend->hgap, bend->fint);
        break;
      case T_CSBEND:
        csbend = (CSBEND*)eptr->p_elem;
        emitCosyDipole(fp, eptr->name, csbend->length, csbend->angle,
                       csbend->e1, csbend->e2, csbend->h1, csbend->h2,
                       csbend->k1, csbend->k2, csbend->k3,
                       csbend->hgap, csbend->fint);
        break;
      case T_CSRCSBEND:
        csrbend = (CSRCSBEND*)eptr->p_elem;
        emitCosyDipole(fp, eptr->name, csrbend->length, csrbend->angle,
                       csrbend->e1, csrbend->e2, csrbend->h1, csrbend->h2,
                       csrbend->k1, csrbend->k2, csrbend->k3,
                       csrbend->hgap, csrbend->fint);
        break;
      case T_DRIF:
        drift = (DRIFT*)eptr->p_elem;
        emitCosyDrift(fp, eptr->name, drift->length);
        break;
      case T_CSRDRIFT:
        csrdrift = (CSRDRIFT*)eptr->p_elem;
        emitCosyDrift(fp, eptr->name, csrdrift->length);
        break;
      case T_SEXT:
        sext = (SEXT*)eptr->p_elem;
        fprintf(fp, "  PROCEDURE %s ;\n", eptr->name);
        fprintf(fp, "    { K2 = %e, R = %e, BRho = %e }\n",
                sext->k2, sextBoreRadius, BRho);
        fprintf(fp, "    MSK %.15g %.15g %.15g ;\n",
                sext->length, sext->k2*BRho*sextBoreRadius, sextBoreRadius);
        fprintf(fp, "  ENDPROCEDURE; \n");
        break;
      case T_KSEXT:
        ksext = (KSEXT*)eptr->p_elem;
        fprintf(fp, "  PROCEDURE %s ;\n", eptr->name);
        fprintf(fp, "    { K2 = %e, R = %e, BRho = %e }\n",
                ksext->k2, sextBoreRadius, BRho);
        fprintf(fp, "    MSK %.15g %.15g %.15g ;\n",
                ksext->length, ksext->k2*BRho*sextBoreRadius, sextBoreRadius);
        fprintf(fp, "  ENDPROCEDURE; \n");
        break;
      case T_HCOR:
        hcor = (HCOR*)eptr->p_elem;
        emitCosyDrift(fp, eptr->name, hcor->length);
        break;
      case T_VCOR:
        vcor = (VCOR*)eptr->p_elem;
        emitCosyDrift(fp, eptr->name, vcor->length);
        break;
      case T_RFCA:
        rfca = (RFCA*)eptr->p_elem;
        fprintf(fp, "  PROCEDURE %s ;\n", eptr->name);
        fprintf(fp, "    VARIABLE V 1 2 2 ;\n");
        fprintf(fp, "    V(1,1) := %.15g ;\n", rfca->volt/1e3);
        fprintf(fp, "    V(1,2) := 0.0 ;\n");
        fprintf(fp, "    V(2,1) := 0.0 ;\n");
        fprintf(fp, "    V(2,2) := 0.0 ;\n");
        if (rfca->length)
          fprintf(fp, "    DL %.15g ;\n", rfca->length/2);
        fprintf(fp, "    RF V 1 %.15g %.15g 0.1 ;\n", rfca->freq, rfca->phase);
        if (rfca->length)
          fprintf(fp, "    DL %.15g ;\n", rfca->length/2);
        fprintf(fp, "  ENDPROCEDURE ;\n");
        break;
      default:
        emitCosyDrift(fp, eptr->name, 0.0);
        break;
      }
      eptr = eptr->succ;
    }

    /* emit beamline definition */
    fprintf(fp, "  PROCEDURE MACH ;");
    eptr = &(beamline->elem);
    i = 0;
    while (eptr) {
      if (i%4==0)
        fprintf(fp, "\n    %s ;", eptr->name);
      else
        fprintf(fp, " %s ;", eptr->name);
      i++;
      eptr = eptr->succ;
    }
    fprintf(fp, "\n  ENDPROCEDURE ;\n");
    
    fprintf(fp, "    OV %ld 2 1 ;\n", cosyOrder);
    fprintf(fp, "    RPE %.15g*PARA(1) ;\n", (sqrt(pCentral*pCentral+1)-1)*me_mks*sqr(c_mks)/e_mks/1e6);
    fprintf(fp, "    UM ; MACH; \n");
    fprintf(fp, "    TP MU; WRITE 7 ' DELTA-DEPENDENT TUNES: '   MU(1) MU(2) ;\n");
    fprintf(fp, "    GT MAP F MU A B G R;\n");
    fprintf(fp, "    WRITE 7 ' DELTA-DEPENDENT FIXED POINT ' F(1) F(2) F(3) F(4) ;\n");
    fprintf(fp, "    WRITE 7 ' DELTA-DEPENDENT ALPHAS ' A(1) A(2) ;\n");
    fprintf(fp, "    WRITE 7 ' DELTA-DEPENDENT BETAS  ' B(1) B(2) ;\n");
    fprintf(fp, "    WRITE 7 ' DELTA-DEPENDENT GAMMAS ' G(1) G(2) ;\n");
    fprintf(fp, "ENDPROCEDURE ;\n");
    fprintf(fp, "RUN ;\n");
    fprintf(fp, "END ;\n");
    }

