#include "mdb.h"
#include "SDDS.h"
#include "scan.h"

#define C_LIGHT_MKS 299792458.
#define E_CHARGE_PC 1.602176487e-7
#define E_MASS_MKS 9.1093837015e-31

typedef struct FRINGE_INT2
{
  double int1;
  double int2;
} FRINGE_INT2;
typedef struct FRINGE_INT3
{
  double int1;
  double int2;
  double int3;
} FRINGE_INT3;

void readInGGE(char *ggeFile, double **z, double **ggeD, double **ggeQ, double **ggeS, int *Nz);
void readInTrajectoryInput(char *input, double **x, double **xp, double *dz, int64_t *rows,
                           double *pCentral, int *Nsegments, double **zMagnetRef);

void LGBENDfringeCalc(double *z, double *ggeD, double *ggeQ, double *ggeS, double *xp,
                      double *x, int Nz, double invRigidity, double *zRef, int Nedges, char *output);

void CCBENDfringeCalc(double *z, double *ggeD, double *ggeQ, double *ggeS,
                      int Nz, double invRigidity, double zRef, double bendAngle, char *output, char *elementName);

//int refineNextMaximum(double *ggeD, int *zMaxInt, int edgeNum);

double integrateTrap(double *integrand, int startPt, int endPt, double dz);

void setStepFunction(double *heavside, double *maxD, int *zMaxInt, int *zEdgeInt, int edgeNum, double checkInt, double dz);

FRINGE_INT3 computeDipoleFringeInt(double *ggeD, double *heavside, double *z, double *zEdge,
                                   int *zMaxInt, int edgeNum, double invRigidity, double dz);

FRINGE_INT3 computeQuadrupoleFringeInt(double *ggeQ, double *heavside, double *z, double *zEdge,
                                       int *zMaxInt, int edgeNum, double invRigidity, double dz);

FRINGE_INT3 computeSextupoleFringeInt(double *ggeS, double *stepFuncS, double *z, double *zEdge,
                                      int *zMaxInt, int edgeNum, double invRigidity, double dz);

FRINGE_INT3 computeDipSextFringeInt(double *ggeD, double *stepFuncD, double *ggeS, double *stepFuncS,
                                    double *z, double *zEdge, int *zMaxInt, int edgeNum, double invRigidity, double dz);

#define SET_GGE 0
#define SET_CCBEND 1
#define SET_LGBEND 2
#define SET_PCENTRAL 3
#define SET_BENDANGLE 4
#define SET_ZREFFIELD 5
#define SET_TRAJECTORY 6
#define SET_ELEMENT_NAME 7
#define N_OPTIONS 8
char *option[N_OPTIONS] = {
  "gge",
  "ccbend",
  "lgbend",
  "pCentral",
  "bendAngle",
  "zRefField",
  "trajectory",
  "elementName",
};

#define USAGE "straightDipoleFringeCalc -gge=<filename> <outputfile>\n\
   -elementName=<string>\n\
   { -ccbend \n\
     -pCentral=<value> -bendAngle=<value> -zRefField=<value>\n\
   | -lgbend \n\
     -trajectory=<filename>\n\
   }\n\
\n\
-gge         Give generalized gradient expansion file, such as produced by computeRBGGE\n\
             or computeCBGGE.\n\
<outputFile> Name of file to which to write fringe parameters in the format expected by\n\
             elegant's load_parameters command\n\
-elementName Name of the element to which the data will be loaded in elegant.\n\
-ccbend      CCBEND mode. This requires three additioanl parameters.\n\
 -pCentral    Reference momentum\n\
 -bendAngle   Full bending angle in dipole\n\
 -zRefField   z location in gge file where reference magnetic field values will be taken\n\
-lgbend      LGBEND mode, which requires one additional parameter\n\
 -trajectory  SDDS file is from BGGEXP PARTICLE_OUTPUT_FILE feature, for a single\n\
              reference particle.\n\
               Parameters:\n\
                 pCentral          Reference momentum beta*gamma \n\
                 nMagnetSegments   Number of flat steps of non-zero field\n\
                 zReference#       z location in file where reference\n\
                                   magnetic field values will be taken,\n\
                                   1<= # <= nMagnetSegments\n\
               Columns: z, x, px, and pz where z is the same as that in the gge file\n\
                       (usually obtained from elegant tracking)"

int main(int argc, char **argv)
{
  SCANNED_ARG *scanned;
  long i_arg;

  double *z=NULL, *ggeD=NULL, *ggeQ=NULL, *ggeS=NULL, *xp, *x;
  double *zMagnetRef;

  double pCentral=DBL_MAX, invRigidity, bendAngle=DBL_MAX;

  int Nsegments, Nedges, ip, Nz;

  char *ggeFile=NULL, *output=NULL, *trajectoryFile = NULL, *elementName = NULL;
  int mode=0;

  zMagnetRef = malloc(sizeof(double) * 10);
  for (ip = 0; ip < 10; ip++)
    zMagnetRef[ip] = DBL_MAX;

  argc = scanargs(&scanned, argc, argv);
  if (argc < 2 || argc > (2 + N_OPTIONS)) {
    fprintf(stderr, "%s\n", USAGE);
    return (1);
  }
  for (i_arg = 1; i_arg < argc; i_arg++) {
    if (scanned[i_arg].arg_type == OPTION) {
      /* process options here */
      switch (match_string(scanned[i_arg].list[0], option, N_OPTIONS, 0)) {
      case SET_GGE:
        if (scanned[i_arg].n_items != 2) {
          fprintf(stderr, "invalid -gge syntax\n%s\n", USAGE);
          return (1);
        }
        ggeFile = scanned[i_arg].list[1];
        break;
      case SET_CCBEND:
        if (mode != 0)
          {
            fprintf(stderr, "too many modes given\n%s\n", USAGE);
            return (1);
          }
        mode = 1;
        break;
      case SET_LGBEND:
        if (mode != 0)
          {
            fprintf(stderr, "too many modes given\n%s\n", USAGE);
            return (1);
          }
        mode = 2;
        break;
      case SET_PCENTRAL:
        if (scanned[i_arg].n_items!=2 ||
            !sscanf(scanned[i_arg].list[1], "%lf", &pCentral))
          bomb("invalid -pCentral syntax", USAGE);
        break;
      case SET_BENDANGLE:
        if (scanned[i_arg].n_items!=2 ||
            !sscanf(scanned[i_arg].list[1], "%lf", &bendAngle))
          bomb("invalid -bendAngle syntax", USAGE);
        break;
      case SET_ZREFFIELD:
        if (scanned[i_arg].n_items!=2 ||
            !sscanf(scanned[i_arg].list[1], "%lf", &(zMagnetRef[0])))
          bomb("invalid -zRefField syntax", USAGE);
        break;
        break;
      case SET_TRAJECTORY:
        if (scanned[i_arg].n_items!=2)
          bomb("invalid -trajectory syntax", USAGE);
        trajectoryFile = scanned[i_arg].list[1];
        break;
      case SET_ELEMENT_NAME:
        if (scanned[i_arg].n_items!=2)
          bomb("invalid -elementName syntax", USAGE);
        elementName = scanned[i_arg].list[1];
        break;
      default:
        fprintf(stderr, "unknown option given\n%s\n", USAGE);
        return (1);
        break;
      }
    } else {
      if (!output)
        output = scanned[i_arg].list[0];
      else
        {
          fprintf(stderr, "too many filenames given\n%s\n", USAGE);
          return(1);
        }
    }
  }
  if (elementName == NULL)
    bomb("Must give -elementName option", NULL);
  if (ggeFile == NULL)
    {
      fprintf(stderr, "Must include -gge option\n");
      exit(1);      
    }  
  if (output == NULL)
    {
      fprintf(stderr, "No output file given\n");
      exit(1);      
    }  
  if (mode == 0)
    {
      fprintf(stderr, "Must include -ccbend or -lgbend option\n");
      exit(1);      
    }
  if (mode == 1)
    {
      if (pCentral == DBL_MAX)
        {
          fprintf(stderr, "Missing -pCentral option\n");
          exit(1);
        }
      if (bendAngle == DBL_MAX)
        {
          fprintf(stderr, "Missing -bendAngle option\n");
          exit(1);
        }
       if (zMagnetRef[0] == DBL_MAX)
        {
          fprintf(stderr, "Missing -zRefField option\n");
          exit(1);
        }
   }
  /* get GGE input */
  readInGGE(ggeFile, &z, &ggeD, &ggeQ, &ggeS, &Nz);
  /****** Note that ggeD = -CnmS0[m=1]; ggeQ = -2*CnmS0[m=2]; ggeS = -6*CnmS0[m=3] + 0.25*CnmS2[m=1] ******/

  if (mode == 1)
    {
      Nsegments = 1;
      Nedges = Nsegments + 1;
      zMagnetRef[0] -= z[0];
      for (ip = Nz - 1; ip >= 0; ip--)
        z[ip] -= z[0];
      /* add 0.01% of dz to eliminate rounding errors when converting to ints */
      //zMagnetRef[0] += 1.0e-4 * (z[1] - z[0]);
      invRigidity = 1.0e-12 * E_CHARGE_PC / (E_MASS_MKS * C_LIGHT_MKS * pCentral);
      CCBENDfringeCalc(z, ggeD, ggeQ, ggeS, Nz, invRigidity, zMagnetRef[0], bendAngle, output, elementName);
    }

  if (mode == 2)
    {
      int64_t trajRows, i;
      double dzTraj;
      if (!trajectoryFile || !strlen(trajectoryFile))
        bomb("trajectory file must be given in LGBEND mode", USAGE);

      readInTrajectoryInput(trajectoryFile, &x, &xp, &dzTraj, &trajRows, &pCentral, &Nsegments, &zMagnetRef);
      if (trajRows!=(Nz+1)) {
        fprintf(stderr, "Error: %ld rows in file %s but %ld rows in file %s\n",
                (long)trajRows, trajectoryFile, (long)Nz, ggeFile);
        exit(1);
      }
      for (i=0; i<Nz-1; i++)
        if (fabs(dzTraj-(z[i+1]-z[i]))>1e-6*dzTraj) {
          fprintf(stderr, "dz values are mismatched between the GGE and trajectory files\n");
          exit(1);
        }

      Nedges = Nsegments + 1;
      printf("Nedges = %i\n", Nedges);
      /* Internally set z[0] = 0 to simplify calculations */
      for (ip = 0; ip < Nedges - 1; ip++)
        zMagnetRef[ip] -= z[0];
      for (ip = Nz - 1; ip >= 0; ip--)
        z[ip] -= z[0];
      /* add 0.01% of dz to eliminate rounding errors when converting to ints */
      for (ip = 0; ip < Nedges - 1; ip++)
        zMagnetRef[ip] += 1.0e-4 * (z[1] - z[0]);
      invRigidity = 1.0e-12 * E_CHARGE_PC / (E_MASS_MKS * C_LIGHT_MKS * pCentral);
      LGBENDfringeCalc(z, ggeD, ggeQ, ggeS, xp, x, Nz, invRigidity, zMagnetRef, Nedges, output);
    }

  //fclose(outFP);
  free(z);
  free(ggeD);
  free(ggeQ);
  free(ggeS);
  exit(0);
}

void LGBENDfringeCalc(double *z, double *ggeD, double *ggeQ, double *ggeS, double *xp,
                      double *x, int Nz, double invRigidity, double *zRef, int Nedges, char *output)
{
  SDDS_DATASET SDDSout;
  double *stepFuncD, *stepFuncQ, *stepFuncS;
  double *zEdge, *maxD, *maxQ, *maxS;
  double *bendRad, *edgeAngle, *edgeX;

  FRINGE_INT3 dipFringeInt, sextFringeInt;
  FRINGE_INT3 quadFringeInt;

  double temp1, fringeInt, dz;

  int *zEdgeInt, *zMaxInt;
  int edgeNum, ip;

  dz = z[1] - z[0];
  /* allocate memory */
  stepFuncD = calloc(Nz, sizeof(double));
  stepFuncQ = calloc(Nz, sizeof(double));
  stepFuncS = calloc(Nz, sizeof(double));

  /* allocate memory */
  maxD = calloc(Nedges + 1, sizeof(double));
  maxQ = calloc(Nedges + 1, sizeof(double));
  maxS = calloc(Nedges + 1, sizeof(double));
  zEdge = calloc(Nedges, sizeof(double));
  bendRad = calloc(Nedges, sizeof(double));
  edgeAngle = calloc(Nedges, sizeof(double));
  edgeX = calloc(Nedges, sizeof(double));
  zEdgeInt = calloc(Nedges, sizeof(int));
  zMaxInt = calloc(Nedges + 1, sizeof(int));

  zMaxInt[0] = 0;
  for (edgeNum = 1; edgeNum < Nedges; edgeNum++)
    {
      zMaxInt[edgeNum] = (int)(zRef[edgeNum - 1] / dz);
      if ((zMaxInt[edgeNum] < 1) || (zMaxInt[edgeNum] > Nz - 2))
        {
          printf("Magnetic reference point is located outside the field map!\n");
          return;
        }
    }
  zMaxInt[Nedges] = Nz - 1;

  /* define multipole content for each body segment including outside magnet */
  for (edgeNum = 0; edgeNum < Nedges + 1; edgeNum++)
    {
      maxD[edgeNum] = ggeD[zMaxInt[edgeNum]];
      maxQ[edgeNum] = ggeQ[zMaxInt[edgeNum]];
      maxS[edgeNum] = ggeS[zMaxInt[edgeNum]];
      /* set multipole content outside magnet to small value */
      if ((edgeNum == 0) || (edgeNum == Nedges))
        {
          maxD[edgeNum] = 1.0e-10 * maxD[edgeNum];
          maxQ[edgeNum] = 1.0e-10 * maxQ[edgeNum];
          maxS[edgeNum] = 1.0e-10 * maxS[edgeNum];
        }
    }

  if (!SDDS_InitializeOutput(&SDDSout, SDDS_ASCII, 1, NULL, "Fringe integrals and body parameters for LGBEND", output))
    {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      exit(1);
    }
  if (!SDDS_DefineSimpleParameter(&SDDSout, "FringeSegmentNum", NULL, SDDS_LONG) ||
      !SDDS_DefineSimpleColumn(&SDDSout, "ElementName", NULL, SDDS_STRING) ||
      !SDDS_DefineSimpleColumn(&SDDSout, "ElementParameter", NULL, SDDS_STRING) ||
      !SDDS_DefineSimpleColumn(&SDDSout, "ParameterValue", NULL, SDDS_DOUBLE) ||
      !SDDS_WriteLayout(&SDDSout))
    {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      exit(1);
    }

  for (edgeNum = 0; edgeNum < Nedges; edgeNum++)
    {
      /* find hard edge according to integrated Bfield */
      fringeInt = integrateTrap(ggeD, zMaxInt[edgeNum], zMaxInt[edgeNum + 1], dz);
      zEdge[edgeNum] = (maxD[edgeNum + 1] * (z[zMaxInt[edgeNum + 1]] - z[zMaxInt[edgeNum]]) - fringeInt) / (maxD[edgeNum + 1] - maxD[edgeNum]);
      zEdge[edgeNum] += z[zMaxInt[edgeNum]];
      temp1 = (zEdge[edgeNum] - z[(int)(zEdge[edgeNum] / dz)]) / dz;
      /* angles at first and last edges set by input x' */
      if (edgeNum == 0)
        edgeAngle[edgeNum] = atan(xp[0]);
      else
        {
          if (edgeNum == Nedges - 1)
            edgeAngle[edgeNum] = atan(xp[Nz - 1]);
          else
            edgeAngle[edgeNum] = atan((1.0 - temp1) * xp[(int)(zEdge[edgeNum] / dz)] + temp1 * xp[(int)(zEdge[edgeNum] / dz) + 1]);
          /* interior angles set by linear interpolation to edge */
        }
      /* horizontal coordinate set by linear interpolation to edge */
      edgeX[edgeNum] = (1.0 - temp1) * x[(int)(zEdge[edgeNum] / dz)] + temp1 * x[(int)(zEdge[edgeNum] / dz) + 1];

      bendRad[edgeNum] = invRigidity * maxD[edgeNum + 1];
      bendRad[edgeNum] = 1.0 / bendRad[edgeNum];

      /* find step function field that gives same integrated field */
      zEdgeInt[edgeNum] = (int)(zEdge[edgeNum] / dz);
      setStepFunction(stepFuncD, maxD, zMaxInt, zEdgeInt, edgeNum, fringeInt, dz);

      for (ip = zMaxInt[edgeNum]; ip < zEdgeInt[edgeNum]; ip++)
        stepFuncQ[ip] = maxQ[edgeNum];
      stepFuncQ[zEdgeInt[edgeNum]] = ((maxQ[edgeNum] - maxQ[edgeNum + 1]) * stepFuncD[zEdgeInt[edgeNum]] - maxQ[edgeNum] * maxD[edgeNum + 1] + maxQ[edgeNum + 1] * maxD[edgeNum]) / (maxD[edgeNum] - maxD[edgeNum + 1]);
      for (ip = zEdgeInt[edgeNum] + 1; ip <= zMaxInt[edgeNum + 1]; ip++)
        stepFuncQ[ip] = maxQ[edgeNum + 1];

      for (ip = zMaxInt[edgeNum]; ip < zEdgeInt[edgeNum]; ip++)
        stepFuncS[ip] = maxS[edgeNum];
      stepFuncS[zEdgeInt[edgeNum]] = ((maxS[edgeNum] - maxS[edgeNum + 1]) * stepFuncD[zEdgeInt[edgeNum]] - maxS[edgeNum] * maxD[edgeNum + 1] + maxS[edgeNum + 1] * maxD[edgeNum]) / (maxD[edgeNum] - maxD[edgeNum + 1]);
      for (ip = zEdgeInt[edgeNum] + 1; ip <= zMaxInt[edgeNum + 1]; ip++)
        stepFuncS[ip] = maxS[edgeNum + 1];

      dipFringeInt = computeDipoleFringeInt(ggeD, stepFuncD, z, zEdge, zMaxInt, edgeNum, invRigidity, dz);
      quadFringeInt = computeQuadrupoleFringeInt(ggeQ, stepFuncQ, z, zEdge, zMaxInt, edgeNum, invRigidity, dz);
      sextFringeInt = computeSextupoleFringeInt(ggeS, stepFuncS, z, zEdge, zMaxInt, edgeNum, invRigidity, dz);

      if (!SDDS_StartPage(&SDDSout, 12) ||
          !SDDS_SetParameters(&SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "FringeSegmentNum", edgeNum + 1, NULL) ||
          !SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, 0,  "ElementName", "LGBEND", "ElementParameter", "intK0",  \
                              "ParameterValue", dipFringeInt.int2, NULL) ||
          !SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, 1,  "ElementName", "LGBEND", "ElementParameter", "intK2",  \
                              "ParameterValue", dipFringeInt.int3, NULL) ||
          !SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, 2,  "ElementName", "LGBEND", "ElementParameter", "intK4",  \
                              "ParameterValue", sextFringeInt.int3, NULL) ||
          !SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, 3,  "ElementName", "LGBEND", "ElementParameter", "intK5",  \
                              "ParameterValue", sextFringeInt.int2, NULL) ||
          !SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, 4,  "ElementName", "LGBEND", "ElementParameter", "intK6",  \
                              "ParameterValue", sextFringeInt.int1, NULL) ||
          !SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, 5,  "ElementName", "LGBEND", "ElementParameter", "intKIK", \
                              "ParameterValue", quadFringeInt.int1, NULL) ||
          !SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, 6,  "ElementName", "LGBEND", "ElementParameter", "intKII", \
                              "ParameterValue", quadFringeInt.int2, NULL) ||
          !SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, 7,  "ElementName", "LGBEND", "ElementParameter", "LONGIT_L",      \
                              "ParameterValue", fabs(fringeInt / maxD[edgeNum + 1]), NULL) ||
          !SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, 8,  "ElementName", "LGBEND", "ElementParameter", "EDGE_ANGLE",  \
                              "ParameterValue", edgeAngle[edgeNum], NULL) ||
          !SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, 9,  "ElementName", "LGBEND", "ElementParameter", "EDGE_X",     \
                              "ParameterValue", edgeX[edgeNum], NULL) ||
          !SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, 10, "ElementName", "LGBEND", "ElementParameter", "K1",     \
                              "ParameterValue", invRigidity * maxQ[edgeNum + 1], NULL) ||
          !SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, 11, "ElementName", "LGBEND", "ElementParameter", "K2",     \
                              "ParameterValue", invRigidity * (maxS[edgeNum + 1] - 0.25 * maxD[edgeNum + 1]), NULL) ||
          !SDDS_WritePage(&SDDSout))
        {
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
          exit(1);
        }
    }
  if (!SDDS_Terminate(&SDDSout))
    {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      exit(1);
    }
  free(stepFuncD);
  free(stepFuncQ);
  free(stepFuncS);
  free(maxD);
  free(maxQ);
  free(maxS);
  free(zEdge);
  free(bendRad);
  free(edgeAngle);
  free(edgeX);
  free(zEdgeInt);
  free(zMaxInt);
}

void CCBENDfringeCalc(double *z, double *ggeD, double *ggeQ, double *ggeS, int Nz, double invRigidity, 
                      double zRef, double bendAngle, char *output, char *elementName)
{
  SDDS_DATASET SDDSout;
  double *stepFuncD, *stepFuncQ, *stepFuncS;
  double *zEdge, *maxD, *maxQ, *maxS;

  FRINGE_INT3 dipFringeInt, sextFringeInt;
  FRINGE_INT3 quadFringeInt, dipSextFringeInt;

  double arcLength, fringeInt, dz = z[1] - z[0];

  int *zEdgeInt, *zMaxInt;
  int ip, edgeNum, Nedges = 2, iRow;

  /* allocate memory */
  stepFuncD = calloc(Nz, sizeof(double));
  stepFuncQ = calloc(Nz, sizeof(double));
  stepFuncS = calloc(Nz, sizeof(double));

  maxD = calloc(Nedges + 1, sizeof(double));
  maxQ = calloc(Nedges + 1, sizeof(double));
  maxS = calloc(Nedges + 1, sizeof(double));
  zEdge = calloc(Nedges, sizeof(double));

  zEdgeInt = calloc(Nedges, sizeof(int));
  zMaxInt = calloc(Nedges + 1, sizeof(int));

  /* set position of the reference points */
  zMaxInt[0] = 0;
  zMaxInt[1] = (int)(zRef / dz);
  if ((zMaxInt[1] < 1) || (zMaxInt[1] > Nz - 2))
    {
      printf("Magnetic reference point is located outside the field map!\n");
      return;
    }
  zMaxInt[2] = Nz - 1;

  /* set multipole content outside magnet to small value */
  maxD[0] = 1.0e-10 * ggeD[zMaxInt[0]];
  maxQ[0] = 1.0e-10 * ggeQ[zMaxInt[0]];
  maxS[0] = 1.0e-10 * ggeS[zMaxInt[0]];
  maxD[2] = 1.0e-10 * ggeD[zMaxInt[2]];
  maxQ[2] = 1.0e-10 * ggeQ[zMaxInt[2]];
  maxS[2] = 1.0e-10 * ggeS[zMaxInt[2]];

  /* set multipole content in body */
  maxD[1] = ggeD[zMaxInt[1]];
  maxQ[1] = ggeQ[zMaxInt[1]];
  maxS[1] = ggeS[zMaxInt[1]];

  if (!SDDS_InitializeOutput(&SDDSout, SDDS_ASCII, 1, NULL, "Fringe integrals and body parameters for CCBEND", output))
    {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      exit(1);
    }
  if (!SDDS_DefineSimpleColumn(&SDDSout, "ElementName", NULL, SDDS_STRING) ||
      !SDDS_DefineSimpleColumn(&SDDSout, "ElementType", NULL, SDDS_STRING) ||
      !SDDS_DefineSimpleColumn(&SDDSout, "ElementParameter", NULL, SDDS_STRING) ||
      !SDDS_DefineSimpleColumn(&SDDSout, "ParameterValue", NULL, SDDS_DOUBLE) ||
      !SDDS_WriteLayout(&SDDSout) ||
      !SDDS_StartPage(&SDDSout, 12*Nedges) ) 
    {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      exit(1);
    }

  for (edgeNum = iRow = 0; edgeNum < Nedges; edgeNum++)
    {
      char nameBuffer[256];
      /* find hard edge according to integrated Bfield */
      fringeInt = integrateTrap(ggeD, zMaxInt[edgeNum], zMaxInt[edgeNum + 1], dz);
      zEdge[edgeNum] = (maxD[edgeNum + 1] * (z[zMaxInt[edgeNum + 1]] - z[zMaxInt[edgeNum]]) - fringeInt) / (maxD[edgeNum + 1] - maxD[edgeNum]);
      zEdge[edgeNum] += z[zMaxInt[edgeNum]];

      /* find step function field that gives same integrated field */
      zEdgeInt[edgeNum] = (int)(zEdge[edgeNum] / dz);
      setStepFunction(stepFuncD, maxD, zMaxInt, zEdgeInt, edgeNum, fringeInt, dz);

      for (ip = zMaxInt[edgeNum]; ip < zEdgeInt[edgeNum]; ip++)
        stepFuncQ[ip] = maxQ[edgeNum];
      stepFuncQ[zEdgeInt[edgeNum]] = ((maxQ[edgeNum] - maxQ[edgeNum + 1]) * stepFuncD[zEdgeInt[edgeNum]] - maxQ[edgeNum] * maxD[edgeNum + 1] + maxQ[edgeNum + 1] * maxD[edgeNum]) / (maxD[edgeNum] - maxD[edgeNum + 1]);
      for (ip = zEdgeInt[edgeNum] + 1; ip <= zMaxInt[edgeNum + 1]; ip++)
        stepFuncQ[ip] = maxQ[edgeNum + 1];

      for (ip = zMaxInt[edgeNum]; ip < zEdgeInt[edgeNum]; ip++)
        stepFuncS[ip] = maxS[edgeNum];
      stepFuncS[zEdgeInt[edgeNum]] = ((maxS[edgeNum] - maxS[edgeNum + 1]) * stepFuncD[zEdgeInt[edgeNum]] - maxS[edgeNum] * maxD[edgeNum + 1] + maxS[edgeNum + 1] * maxD[edgeNum]) / (maxD[edgeNum] - maxD[edgeNum + 1]);
      for (ip = zEdgeInt[edgeNum] + 1; ip <= zMaxInt[edgeNum + 1]; ip++)
        stepFuncS[ip] = maxS[edgeNum + 1];

      dipFringeInt = computeDipoleFringeInt(ggeD, stepFuncD, z, zEdge, zMaxInt, edgeNum, invRigidity, dz);
      quadFringeInt = computeQuadrupoleFringeInt(ggeQ, stepFuncQ, z, zEdge, zMaxInt, edgeNum, invRigidity, dz);
      sextFringeInt = computeSextupoleFringeInt(ggeS, stepFuncS, z, zEdge, zMaxInt, edgeNum, invRigidity, dz);

      dipSextFringeInt = computeDipSextFringeInt(ggeD, stepFuncD, ggeS, stepFuncS, z, zEdge, zMaxInt, edgeNum, invRigidity, dz);

      // This needs to be looked at more...
      if(edgeNum == 1)
	arcLength = 0.5*bendAngle*(zEdge[1]-zEdge[0])/sin(0.5*bendAngle);
      else
	arcLength = bendAngle / (invRigidity * maxD[1]);

      if (!SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, iRow++,
                             "ElementName", elementName, "ElementType", "CCBEND", 
                             "ElementParameter", snprintf(nameBuffer, 256, "FRINGE%dK0", edgeNum+1)?nameBuffer:"?",
                              "ParameterValue", dipFringeInt.int2, NULL) ||
          !SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, iRow++,
                             "ElementName", elementName, "ElementType", "CCBEND", 
                             "ElementParameter", snprintf(nameBuffer, 256, "FRINGE%dK2", edgeNum+1)?nameBuffer:"?",
                              "ParameterValue", dipFringeInt.int3, NULL) ||
          !SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, iRow++,
                             "ElementName", elementName, "ElementType", "CCBEND", 
                             "ElementParameter", snprintf(nameBuffer, 256, "FRINGE%dK4", edgeNum+1)?nameBuffer:"?",
                              "ParameterValue", sextFringeInt.int3, NULL) ||
          !SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, iRow++,
                             "ElementName", elementName, "ElementType", "CCBEND", 
                             "ElementParameter", snprintf(nameBuffer, 256, "FRINGE%dK5", edgeNum+1)?nameBuffer:"?",
                              "ParameterValue", sextFringeInt.int2, NULL) ||
          !SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, iRow++,
                             "ElementName", elementName, "ElementType", "CCBEND", 
                             "ElementParameter", snprintf(nameBuffer, 256, "FRINGE%dK6", edgeNum+1)?nameBuffer:"?",
                              "ParameterValue", sextFringeInt.int1, NULL) ||
          !SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, iRow++,
                             "ElementName", elementName, "ElementType", "CCBEND", 
                             "ElementParameter", snprintf(nameBuffer, 256, "FRINGE%dK7", edgeNum+1)?nameBuffer:"?",
                              "ParameterValue", dipSextFringeInt.int1+dipSextFringeInt.int2, NULL) ||
          !SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, iRow++,
                             "ElementName", elementName, "ElementType", "CCBEND", 
                             "ElementParameter", snprintf(nameBuffer, 256, "FRINGE%dI0", edgeNum+1)?nameBuffer:"?",
                              "ParameterValue", quadFringeInt.int1, NULL) ||
          !SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, iRow++,
                             "ElementName", elementName, "ElementType", "CCBEND", 
                             "ElementParameter", snprintf(nameBuffer, 256, "FRINGE%dI1", edgeNum+1)?nameBuffer:"?",
                             "ParameterValue", quadFringeInt.int2, NULL))
        {
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
          exit(1);
        }
      if (edgeNum==0) 
        {
          if (!SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, iRow++,
                                 "ElementName", elementName, "ElementType", "CCBEND", 
                                 "ElementParameter", "L", 
                                 "ParameterValue", arcLength, NULL) ||
              !SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, iRow++,
                                 "ElementName", elementName, "ElementType", "CCBEND", 
                                 "ElementParameter", "ANGLE",
                                 "ParameterValue", bendAngle, NULL) ||
              !SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, iRow++,
                                 "ElementName", elementName, "ElementType", "CCBEND", 
                                 "ElementParameter", "K1", 
                                 "ParameterValue", invRigidity * maxQ[edgeNum + 1], NULL) ||
              !SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, iRow++,
                                 "ElementName", elementName, "ElementType", "CCBEND", 
                                 "ElementParameter", "K2", 
                                 "ParameterValue", invRigidity * (maxS[edgeNum + 1] - 0.25 * maxD[edgeNum + 1]), NULL))
            {
              SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
              exit(1);
            }
        }
    }
  if (!SDDS_WritePage(&SDDSout) || !SDDS_Terminate(&SDDSout))
    {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      exit(1);
    }
  free(stepFuncD);
  free(stepFuncQ);
  free(stepFuncS);
  free(maxD);
  free(maxQ);
  free(maxS);
  free(zEdge);
  free(zEdgeInt);
  free(zMaxInt);
}

void readInGGE(char *ggeFile, double **z, double **ggeD, double **ggeQ, double **ggeS, int *Nz)
{
  SDDS_DATASET SDDSin;
  int64_t i, rows;
  double *zValues, *CnmS0Values, *CnmS2Values;

  if (!SDDS_InitializeInput(&SDDSin, ggeFile))
    {
      fprintf(stderr, "Unable to read file %s\n", ggeFile);
      exit(1);
    }
  if (SDDS_CheckColumn(&SDDSin, "z", "m", SDDS_ANY_FLOATING_TYPE, stderr)!=SDDS_CHECK_OK)
    {
      fprintf(stderr, "Unable to find floating-point column \"z\" with units \"m\" in file %s\n", ggeFile);
      exit(1);
    }
  if (SDDS_CheckColumn(&SDDSin, "CnmS0", NULL, SDDS_ANY_FLOATING_TYPE, stderr)!=SDDS_CHECK_OK)
    {
      fprintf(stderr, "Unable to find floating-point column \"CnmS0\" in file %s\n", ggeFile);
      exit(1);
    }
  if (SDDS_CheckColumn(&SDDSin, "CnmS2", NULL, SDDS_ANY_FLOATING_TYPE, stderr)!=SDDS_CHECK_OK)
    {
      fprintf(stderr, "Unable to find floating-point column \"CnmS2\" in file %s\n", ggeFile);
      exit(1);
    }

  //Read first page
  if (SDDS_ReadPage(&SDDSin) <= 0)
    {
      fprintf(stderr, "Unable read the first page of %s\n", ggeFile);
      exit(1);
    }
  if ((rows = SDDS_RowCount(&SDDSin))<=1)
    {
      fprintf(stderr, "Too few z values in file %s\n", ggeFile);
      exit(1);
    }
  if (!(zValues=SDDS_GetColumnInDoubles(&SDDSin, "z")))
    {
      fprintf(stderr, "Problem reading column z from %s\n", ggeFile);
      exit(1);
    }
  if (!(CnmS0Values=SDDS_GetColumnInDoubles(&SDDSin, "CnmS0")))
    {
      fprintf(stderr, "Problem reading column CnmS0 from %s\n", ggeFile);
      exit(1);
    }
  if (!(CnmS2Values=SDDS_GetColumnInDoubles(&SDDSin, "CnmS2")))
    {
      fprintf(stderr, "Problem reading column CnmS2 from %s\n", ggeFile);
      exit(1);
    }
  *Nz = rows;
  *z = zValues;
  *ggeD = CnmS0Values;
  *ggeS = CnmS2Values;
  for (i = 0; i < rows; i++)
    {
      (*ggeD)[i] *= -1;
      (*ggeS)[i] *= .25;
    }

  //Read second page
  if (SDDS_ReadPage(&SDDSin) <= 0)
    {
      fprintf(stderr, "Unable read the second page of %s\n", ggeFile);
      exit(1);
    }
  if ((rows = SDDS_RowCount(&SDDSin)) != *Nz)
    {
      fprintf(stderr, "Unequal row count between pages %s\n", ggeFile);
      exit(1);
    }
  if (!(zValues=SDDS_GetColumnInDoubles(&SDDSin, "z")))
    {
      fprintf(stderr, "Problem reading column z from %s\n", ggeFile);
      exit(1);
    }
  if (!(CnmS0Values=SDDS_GetColumnInDoubles(&SDDSin, "CnmS0")))
    {
      fprintf(stderr, "Problem reading column CnmS0 from %s\n", ggeFile);
      exit(1);
    }
  *ggeQ = CnmS0Values;
  for (i = 0; i < rows; i++)
    {
      if ((*z)[i] != zValues[i])
        {
          fprintf(stderr, "z values differ between pages %s\n", ggeFile);
          exit(1);
        }
      (*ggeQ)[i] *= -2;
    }
  free(zValues);

  //Read third page
  if (SDDS_ReadPage(&SDDSin) <= 0)
    {
      fprintf(stderr, "Unable read the third page of %s\n", ggeFile);
      exit(1);
    }
  if ((rows = SDDS_RowCount(&SDDSin)) != *Nz)
    {
      fprintf(stderr, "Unequal row count between pages %s\n", ggeFile);
      exit(1);
    }
  if (!(zValues=SDDS_GetColumnInDoubles(&SDDSin, "z")))
    {
      fprintf(stderr, "Problem reading column z from %s\n", ggeFile);
      exit(1);
    }
  if (!(CnmS0Values=SDDS_GetColumnInDoubles(&SDDSin, "CnmS0")))
    {
      fprintf(stderr, "Problem reading column CnmS0 from %s\n", ggeFile);
      exit(1);
    }
  if (rows != *Nz)
    {
      fprintf(stderr, "Unequal row count between pages %s\n", ggeFile);
      exit(1);
    }
  for (i = 0; i < rows; i++)
    {
      if ((*z)[i] != zValues[i])
        {
          fprintf(stderr, "z values differ between pages %s\n", ggeFile);
          exit(1);
        }
      (*ggeS)[i] -= 6.0 * CnmS0Values[i];
    }
  free(zValues);
  free(CnmS0Values);
  if (!SDDS_Terminate(&SDDSin))
    {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      exit(1);
    }

  return;
}

void readInTrajectoryInput(char *input, double **x, double **xp, double *dz, int64_t *rows,
  double *pCentral, int *Nsegments, double **zMagnetRef)
{
  SDDS_DATASET SDDSin;
  int64_t i;
  int32_t segments;
  double *px, *pz, *z, dz0;
  char name[30];
  int ip;

  if (!SDDS_InitializeInput(&SDDSin, input))
    {
      fprintf(stderr, "Unable to read file %s\n", input);
      exit(1);
    }
  if (SDDS_CheckParameter(&SDDSin, "pCentral", NULL, SDDS_ANY_FLOATING_TYPE, stderr)!=SDDS_CHECK_OK)
    {
      fprintf(stderr, "Unable to find floating-point parameter \"pCentral\" in file %s\n", input);
      exit(1);
    }
  if (SDDS_CheckParameter(&SDDSin, "nMagnetSegments", NULL, SDDS_ANY_INTEGER_TYPE, stderr)!=SDDS_CHECK_OK)
    {
      fprintf(stderr, "Unable to find integer paramenter \"nMagnetSegments\" in file %s\n", input);
      exit(1);
    }
  if (SDDS_CheckColumn(&SDDSin, "x", "m", SDDS_ANY_FLOATING_TYPE, stderr)!=SDDS_CHECK_OK)
    {
      fprintf(stderr, "Unable to find floating-point column \"x\" with units \"m\" in file %s\n", input);
      exit(1);
    }
  if (SDDS_CheckColumn(&SDDSin, "z", "m", SDDS_ANY_FLOATING_TYPE, stderr)!=SDDS_CHECK_OK)
    {
      fprintf(stderr, "Unable to find floating-point column \"x\" with units \"m\" in file %s\n", input);
      exit(1);
    }
  if (SDDS_CheckColumn(&SDDSin, "px", NULL, SDDS_ANY_FLOATING_TYPE, stderr)!=SDDS_CHECK_OK)
    {
      fprintf(stderr, "Unable to find floating-point column \"px\" in file %s\n", input);
      exit(1);
    }
  if (SDDS_CheckColumn(&SDDSin, "pz", NULL, SDDS_ANY_FLOATING_TYPE, stderr)!=SDDS_CHECK_OK)
    {
      fprintf(stderr, "Unable to find floating-point column \"pz\" in file %s\n", input);
      exit(1);
    }
  if (SDDS_ReadPage(&SDDSin) <= 0)
    {
      fprintf(stderr, "Unable read the first page of %s\n", input);
      exit(1);
    }
  if ((*rows = SDDS_RowCount(&SDDSin))<=1)
    {
      fprintf(stderr, "Too few rows in file %s\n", input);
      exit(1);
    }
  if (!SDDS_GetParameterAsDouble(&SDDSin, "pCentral", pCentral))
    {
      fprintf(stderr, "Problem reading parameter pCentral from %s\n", input);
      exit(1);
    }
  if (!SDDS_GetParameterAsLong(&SDDSin, "nMagnetSegments", &segments))
    {
      fprintf(stderr, "Problem reading parameter nMagnetSegments from %s\n", input);
      exit(1);
    }
  *Nsegments = segments;
  for (ip = 1; ip <= *Nsegments; ip++)
    {
      sprintf(name, "zReference%d", ip);
      if (!SDDS_GetParameterAsDouble(&SDDSin, name, &(*zMagnetRef)[ip - 1])) 
        {
          fprintf(stderr, "Problem reading parameter %s from %s\n", name, input);
          exit(1);
        }
    }
  if (!(*x=SDDS_GetColumnInDoubles(&SDDSin, "x")))
    {
      fprintf(stderr, "Problem reading column x from %s\n", input);
      exit(1);
    }
  if (!(z=SDDS_GetColumnInDoubles(&SDDSin, "z")))
    {
      fprintf(stderr, "Problem reading column z from %s\n", input);
      exit(1);
    }
  if ((*dz = (z[*rows-1] - z[0])/(*rows-1))<=0)
    bomb("dz is non-positive in trajectory file", NULL);
  for (i=0; i<(*rows-1); i++)
    if (fabs((dz0 = z[i+1]-z[i])/(*dz)-1)>1e-6) 
      bomb("fractional dz variation of more than 1e-6 seen in trajectory file", NULL);
  if (!(px=SDDS_GetColumnInDoubles(&SDDSin, "px")))
    {
      fprintf(stderr, "Problem reading column px from %s\n", input);
      exit(1);
    }
  if (!(pz=SDDS_GetColumnInDoubles(&SDDSin, "pz")))
    {
      fprintf(stderr, "Problem reading column pz from %s\n", input);
      exit(1);
    }
  if (!(*xp = malloc(sizeof(**xp)*(*rows)))) {
    fprintf(stderr, "Problem allocating memory for xp\n");
  }
  for (i=0; i<*rows; i++) 
    (*xp)[i] = px[i]/pz[i];
  free(px);
  free(pz);
  if (!SDDS_Terminate(&SDDSin))
    {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      exit(1);
    }
}

/*
int refineNextMaximum(double *ggeD, int *zMaxInt, int edgeNum)
{
  double max;
  int search, searchRet, count = 0;

  search = zMaxInt[edgeNum];
  max = ggeD[search];
  search++;
  do
    {
      if (fabs(ggeD[search]) > max)
        {
          max = fabs(ggeD[search]);
          searchRet = search;
          count = 0;
        }
      else
        count++; // prevent local mimima due to noise 
      search++;
    }
  while (count < 10);

  return (searchRet);
}
*/

/*** integrates integrand using the trapezoidal rule from startPt to endPt ***/
double integrateTrap(double *integrand, int startPt, int endPt, double dz)
{
  double integral = 0.0;
  int ip;

  for (ip = startPt; ip < endPt; ip++)
    integral += 0.5 * (integrand[ip] + integrand[ip + 1]);
  integral *= dz;

  return (integral);
}

/*** sets heaviside step function such that the integrated bending field between  ***/
/*** zMaxInt[edgeNum] and zEdgeInt[edgeNum] matches that of field map = checkInt  ***/
void setStepFunction(double *heavside, double *maxD, int *zMaxInt, int *zEdgeInt, int edgeNum, double checkInt, double dz)
{
  double integral;

  int ip, test = -1;

  while (test < 0)
    {
      /* set constant value in segment */
      for (ip = zMaxInt[edgeNum]; ip < zEdgeInt[edgeNum]; ip++)
        heavside[ip] = maxD[edgeNum];
      /* set endpoint to zero  */
      heavside[zEdgeInt[edgeNum]] = 0.0;
      /* set constant value in next segment */
      for (ip = zEdgeInt[edgeNum] + 1; ip <= zMaxInt[edgeNum + 1]; ip++)
        heavside[ip] = maxD[edgeNum + 1];
      integral = integrateTrap(heavside, zMaxInt[edgeNum], zMaxInt[edgeNum + 1], dz);
      /* adjust endpoint in present segment to match integral */
      heavside[zEdgeInt[edgeNum]] = (checkInt - integral) / dz;

      if (maxD[edgeNum] > maxD[edgeNum + 1])
        { /* Assumes field in next segment is smaller */
          if ((maxD[edgeNum] > heavside[zEdgeInt[edgeNum]]) && (heavside[zEdgeInt[edgeNum]] > maxD[edgeNum + 1]))
            test = 1; /* endpoint correction is bewteen values in 2 segments -> success */
          else        /* Move edge location by one step and return to top */
            if (heavside[zEdgeInt[edgeNum]] > maxD[edgeNum])
              zEdgeInt[edgeNum]++;
            else
              zEdgeInt[edgeNum]--;
        }
      else
        { /* Assumes field in next segment is larger */
          if ((maxD[edgeNum + 1] > heavside[zEdgeInt[edgeNum]]) && (heavside[zEdgeInt[edgeNum]] > maxD[edgeNum]))
            test = 1; /* endpoint correction is bewteen values in 2 segments -> success */
          else        /* Move edge location by one step and return to top */
            if (heavside[zEdgeInt[edgeNum]] > maxD[edgeNum + 1])
              zEdgeInt[edgeNum]--;
            else
              zEdgeInt[edgeNum]++;
        }
    }
}

/*** Calculates dipole fringe field contributions ***/
FRINGE_INT3 computeDipoleFringeInt(double *ggeD, double *heavside, double *z, double *zEdge,
                                   int *zMaxInt, int edgeNum, double invRigidity, double dz)
{
  FRINGE_INT3 dipoleFringeInt;
  double intK0 = 0.0;
  double intK1 = 0.0;
  double intK2 = 0.0;
  double sumHeav = heavside[zMaxInt[edgeNum]] + heavside[zMaxInt[edgeNum + 1]];
  double prodHeav = heavside[zMaxInt[edgeNum]] * heavside[zMaxInt[edgeNum + 1]];

  int ip;

  for (ip = zMaxInt[edgeNum]; ip < zMaxInt[edgeNum + 1]; ip++)
    {
      intK0 += 0.5 * ((z[ip] - zEdge[edgeNum]) * (heavside[ip] - ggeD[ip]) + (z[ip + 1] - zEdge[edgeNum]) * (heavside[ip + 1] - ggeD[ip + 1]));
      intK1 += 0.5 * ((ggeD[ip] - heavside[ip]) + (ggeD[ip + 1] - heavside[ip + 1]));
      intK2 += 0.5 * ((sumHeav - ggeD[ip]) * ggeD[ip] + (sumHeav - ggeD[ip + 1]) * ggeD[ip + 1] - 2.0 * prodHeav);
    }
  dipoleFringeInt.int1 = dz * invRigidity * intK1;
  dipoleFringeInt.int2 = dz * invRigidity * intK0;
  dipoleFringeInt.int3 = dz * invRigidity * invRigidity * intK2;

  return (dipoleFringeInt);
}

/*** Calculates quadrupole fringe field contributions ***/
FRINGE_INT3 computeQuadrupoleFringeInt(double *ggeQ, double *stepFuncQ, double *z, double *zEdge,
                                       int *zMaxInt, int edgeNum, double invRigidity, double dz)
{
  FRINGE_INT3 quadFringeInt;
  double intQ0 = 0.0;
  double intQ1 = 0.0;
  double intQ2 = 0.0;

  int ip;

  for (ip = zMaxInt[edgeNum]; ip < zMaxInt[edgeNum + 1]; ip++)
    {
      intQ0 += 0.5 * ((ggeQ[ip] - stepFuncQ[ip]) + (ggeQ[ip + 1] - stepFuncQ[ip + 1]));
      intQ1 += 0.5 * ((z[ip] - zEdge[edgeNum]) * (ggeQ[ip] - stepFuncQ[ip]) + (z[ip + 1] - zEdge[edgeNum]) * (ggeQ[ip + 1] - stepFuncQ[ip + 1]));

      intQ2 += 0.5*( (z[ip]-zEdge[edgeNum])*(z[ip]-zEdge[edgeNum])*(ggeQ[ip] - stepFuncQ[ip])
		     + (z[ip+1]-zEdge[edgeNum])*(z[ip+1]-zEdge[edgeNum])*(ggeQ[ip+1] - stepFuncQ[ip+1]) );
    }
  quadFringeInt.int1 = dz * invRigidity * intQ0;
  quadFringeInt.int2 = dz * invRigidity * intQ1;
  quadFringeInt.int3 = dz * invRigidity * intQ2;

  return (quadFringeInt);
}

/*** Calculates sextupole/curvature ringe field contributions ***/
FRINGE_INT3 computeSextupoleFringeInt(double *ggeS, double *stepFuncS, double *z, double *zEdge,
                                      int *zMaxInt, int edgeNum, double invRigidity, double dz)
{
  FRINGE_INT3 sextFringeInt;
  double intK4 = 0.0;
  double intK5 = 0.0;
  double intK6 = 0.0;

  int ip;

  for (ip = zMaxInt[edgeNum]; ip < zMaxInt[edgeNum + 1]; ip++)
    {
      intK4 += 0.5 * ((z[ip] - zEdge[edgeNum]) * (z[ip] - zEdge[edgeNum]) * (ggeS[ip] - stepFuncS[ip]) + (z[ip + 1] - zEdge[edgeNum]) * (z[ip + 1] - zEdge[edgeNum]) * (ggeS[ip + 1] - stepFuncS[ip + 1]));
      intK5 += 0.5 * ((z[ip] - zEdge[edgeNum]) * (ggeS[ip] - stepFuncS[ip]) + (z[ip + 1] - zEdge[edgeNum]) * (ggeS[ip + 1] - stepFuncS[ip + 1]));
      intK6 += 0.5 * ((ggeS[ip] - stepFuncS[ip]) + (ggeS[ip + 1] - stepFuncS[ip + 1]));
    }
  sextFringeInt.int1 = dz * invRigidity * intK6;
  sextFringeInt.int2 = dz * invRigidity * intK5;
  sextFringeInt.int3 = dz * invRigidity * intK4;

  return (sextFringeInt);
}

FRINGE_INT3 computeDipSextFringeInt(double *ggeD, double *stepFuncD, double *ggeS, double *stepFuncS,
                                    double *z, double *zEdge, int *zMaxInt, int edgeNum, double invRigidity, double dz)
{
  FRINGE_INT3 dipSextFringeInt;

  double intI5 = 0.0;
  double temp6a, temp6b, intI6 = 0.0;
  double temp7, intI7 = 0.0;

  char intFileName[20];
  FILE *outFP;

  int ip, id;

  for(ip=zMaxInt[edgeNum]; ip<zMaxInt[edgeNum+1]; ip++) {
    intI5 += 0.5*invRigidity*(z[ip]-zEdge[edgeNum])*(z[ip]-zEdge[edgeNum])*(ggeS[ip] - stepFuncS[ip])*ggeD[ip];

    temp6a = 0.5*dz*invRigidity*(ggeD[ip] - stepFuncD[ip]);
    temp6b = 0.5*dz*invRigidity*(ggeS[ip] - stepFuncS[ip]);

    temp7 = 0.5*dz*invRigidity*(ggeS[ip] - stepFuncS[ip]);
    for(id=zMaxInt[edgeNum]; id<ip-1; id++) {
      intI6 += 0.5*(z[ip]-z[id])*( temp6a*(ggeS[id] - stepFuncS[id]) + temp6b*(ggeD[id] - stepFuncD[id]) );
      intI6 += 0.5*(z[ip]-z[id+1])*( temp6a*(ggeS[id+1] - stepFuncS[id+1]) + temp6b*(ggeD[id+1] - stepFuncD[id+1]) );

      intI7 += 0.5*temp7*(z[ip]-z[id])*ggeD[id];
      intI7 += 0.5*temp7*(z[ip]-z[id+1])*ggeD[id+1];

    }
    intI5 += 0.5*invRigidity*(z[ip+1]-zEdge[edgeNum])*(z[ip+1]-zEdge[edgeNum])*(ggeS[ip+1] - stepFuncS[ip+1])*ggeD[ip+1];

    temp6a = 0.5*dz*invRigidity*(ggeD[ip+1] - stepFuncD[ip+1]);
    temp6b = 0.5*dz*invRigidity*(ggeS[ip+1] - stepFuncS[ip+1]);

    temp7 = 0.5*dz*invRigidity*(ggeS[ip+1] - stepFuncS[ip+1]);
    for(id=zMaxInt[edgeNum]; id<ip; id++) {
      intI6 += 0.5*(z[ip+1]-z[id])*( temp6a*(ggeS[id] - stepFuncS[id]) + temp6b*(ggeD[id] - stepFuncD[id]) );
      intI6 += 0.5*(z[ip+1]-z[id+1])*( temp6a*(ggeS[id+1] - stepFuncS[id+1]) + temp6b*(ggeD[id+1] - stepFuncD[id+1]) );

      intI7 += 0.5*temp7*(z[ip+1]-z[id])*ggeD[id];
      intI7 += 0.5*temp7*(z[ip+1]-z[id+1])*ggeD[id+1];
    }
  }
  dipSextFringeInt.int1 = 0.5*dz*invRigidity*intI5;
  dipSextFringeInt.int2 = 0.5*dz*invRigidity*intI6;
  dipSextFringeInt.int3 = dz*invRigidity*intI7;

  return(dipSextFringeInt);
}
