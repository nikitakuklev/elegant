/*************************************************************************\
* Copyright (c) 2017 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2017 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* program: sddsmatchmoments.c
 * purpose: transform particle coordinates to match 6D beam moments computed by
 *          moments_output command of elegant.
 *
 * Michael Borland, 2017
 *
 */
#include "mdb.h"
#include "scan.h"
#include "SDDS.h"

#define SET_PIPE 0
#define SET_EXCLUDE 1
#define SET_MOMENTS_FILE 2
#define N_OPTIONS 3

char *option[N_OPTIONS] = {
  "pipe", "exclude", "momentsfile"
} ;

char *USAGE1="sddsmatchmoments [-pipe=[input][,output]] [<SDDSinputfile>] [<SDDSoutputfile>]\n\
  [-momentsFile=<filename>] [-exclude=[{x|y|z}][,centroids]]\n\n";
char *USAGE2="The input file must have columns x, xp, y, yp, and p; for example, an\n\
elegant beam output file is acceptable. \n\n\
-exclude        Exclude one or more planes from the transformation.\n\
-momentsFile    Provide the name of the elegant moments_output file.\n\n\
Program by Michael Borland.  ("__DATE__")\n";

long check_sdds_beam_column(SDDS_TABLE *SDDS_table, char *name, char *units);
void readMomentsFile(char *filename, double Mm[6][6], double Cm[6]);
void transformCoordinates(double *x, double *xp, double *y, double *yp, double *t, double *p, long np,
                          unsigned long flags, double desiredMoment[6][6], double desiredCentroid[6]);
void findTransformationMatrix4(double sigma[4][4], double desiredSigma[4][4], double M[4][4]);

#define EXCLUDE_X 0
#define EXCLUDE_Y 1 
#define EXCLUDE_Z 2
#define EXCLUDE_CENTROIDS 3
#define N_EXCLUDE 4
char *excludeOption[N_EXCLUDE] = {"x", "y", "z", "centroids"};
#define FL_EXCLUDE_X         0x01
#define FL_EXCLUDE_Y         0x02
#define FL_EXCLUDE_Z         0x04
#define FL_EXCLUDE_CENTROIDS 0x08

int main(int argc, char **argv)
{
  SDDS_DATASET SDDSin, SDDSout;
  char *inputfile, *outputfile, *momentsFile;
  long i_arg, rows, readCode;
  SCANNED_ARG *s_arg;
  unsigned long pipeFlags, excludeFlags;
  double *x=NULL, *xp=NULL, *y=NULL, *yp=NULL, *p=NULL, *t=NULL;
  double moment[6][6], centroid[6];

  SDDS_RegisterProgramName(argv[0]);
  argc = scanargs(&s_arg, argc, argv);
  if (argc<2) {
    fprintf(stderr, "%s%s\n", USAGE1, USAGE2);
    return(1);
  }

  inputfile = outputfile = momentsFile = NULL;
  pipeFlags = excludeFlags = 0;
  
  for (i_arg=1; i_arg<argc; i_arg++) {
    if (s_arg[i_arg].arg_type==OPTION) {
      switch (match_string(s_arg[i_arg].list[0], option, N_OPTIONS, 0)) {
      case SET_PIPE:
        if (!processPipeOption(s_arg[i_arg].list+1, s_arg[i_arg].n_items-1, &pipeFlags))
          SDDS_Bomb("invalid -pipe syntax");
        break;
      case SET_EXCLUDE:
        excludeFlags = 0;
        s_arg[i_arg].n_items--;
        if (!scanItemList(&excludeFlags, s_arg[i_arg].list+1, &s_arg[i_arg].n_items, 0, 
                          "x", -1, NULL, 0, FL_EXCLUDE_X,
                          "y", -1, NULL, 0, FL_EXCLUDE_Y,
                          "z", -1, NULL, 0, FL_EXCLUDE_Z,
                          "centroids", -1, NULL, 0, FL_EXCLUDE_CENTROIDS,
                          NULL)) 
          SDDS_Bomb("Invalid -exclude syntax/values");
        if (((excludeFlags&FL_EXCLUDE_X?1:0)+(excludeFlags&FL_EXCLUDE_Y?1:0)+(excludeFlags&FL_EXCLUDE_Z?1:0))>1)
          SDDS_Bomb("Exclude only one of x, y, or z");
        break;
      case SET_MOMENTS_FILE:
        if (s_arg[i_arg].n_items!=2)
          SDDS_Bomb("Invalid -momentsFile syntax/values");
        momentsFile = s_arg[i_arg].list[1];
        break;
      default:
        fprintf(stdout, "error: unknown switch: %s\n", s_arg[i_arg].list[0]);
        fflush(stdout);
        exit(1);
        break;
      }
    }
    else {
      if (inputfile==NULL)
        inputfile = s_arg[i_arg].list[0];
      else if (outputfile==NULL)
        outputfile = s_arg[i_arg].list[0];
      else
        SDDS_Bomb("too many filenames");
    }
  }
  
  processFilenames("sddsmatchmoments", &inputfile, &outputfile, pipeFlags, 0, NULL);

  readMomentsFile(momentsFile, moment, centroid);

  if (!SDDS_InitializeInput(&SDDSin, inputfile))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);

  /* check the input file for valid data */
  if (!check_sdds_beam_column(&SDDSin, "x", "m") || !check_sdds_beam_column(&SDDSin, "y", "m") ||
      !check_sdds_beam_column(&SDDSin, "xp", NULL) || !check_sdds_beam_column(&SDDSin, "yp", NULL) ||
      (!check_sdds_beam_column(&SDDSin, "p", "m$be$nc") && !check_sdds_beam_column(&SDDSin, "p", NULL))) {
    fprintf(stderr, 
            "sddsmatchtwiss: one or more data quantities (x, xp, y, yp, p) have the wrong units or are not present in %s", 
            inputfile);
    exit(1);
  }

  if (!SDDS_InitializeCopy(&SDDSout, &SDDSin, outputfile, "w") ||
      !SDDS_SaveLayout(&SDDSout) || !SDDS_WriteLayout(&SDDSout))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);

  rows = 0;
  while ((readCode=SDDS_ReadPage(&SDDSin))>0) {
    if ((rows=SDDS_RowCount(&SDDSin))==0) {
      if (!SDDS_CopyPage(&SDDSin, &SDDSout) || !SDDS_WritePage(&SDDSout))
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      continue;
    }
    
    if (!(x = SDDS_GetColumnInDoubles(&SDDSin, "x")) ||
        !(xp = SDDS_GetColumnInDoubles(&SDDSin, "xp")) ||
        !(y = SDDS_GetColumnInDoubles(&SDDSin, "y")) ||
        !(yp = SDDS_GetColumnInDoubles(&SDDSin, "yp")) ||
        !(t = SDDS_GetColumnInDoubles(&SDDSin, "t")) ||
        !(p = SDDS_GetColumnInDoubles(&SDDSin, "p")))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);

    /* Transform coordinates */
    transformCoordinates(x, xp, y, yp, t, p, rows, excludeFlags, moment, centroid);

    if (!SDDS_CopyPage(&SDDSout, &SDDSin) ||
        !(SDDS_SetColumnFromDoubles(&SDDSout, SDDS_SET_BY_NAME, x, rows, "x")) ||
        !(SDDS_SetColumnFromDoubles(&SDDSout, SDDS_SET_BY_NAME, y, rows, "y")) ||
        !(SDDS_SetColumnFromDoubles(&SDDSout, SDDS_SET_BY_NAME, xp, rows, "xp")) ||
        !(SDDS_SetColumnFromDoubles(&SDDSout, SDDS_SET_BY_NAME, yp, rows, "yp")) ||
        !(SDDS_SetColumnFromDoubles(&SDDSout, SDDS_SET_BY_NAME, t, rows, "t")) ||
        !(SDDS_SetColumnFromDoubles(&SDDSout, SDDS_SET_BY_NAME, p, rows, "p")) ||
        !SDDS_WritePage(&SDDSout))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    free(x);
    free(xp);
    free(y);
    free(yp);
    free(t);
    free(p);
  }
  if (!SDDS_Terminate(&SDDSout) || !SDDS_Terminate(&SDDSin))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  if (readCode==0)
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  free_scanargs(&s_arg,argc);
  return 0;
}

long check_sdds_beam_column(SDDS_TABLE *SDDS_table, char *name, char *units)
{
  char *units1;
  if (SDDS_GetColumnIndex(SDDS_table, name)<0)
    return(0);
  if (SDDS_GetColumnInformation(SDDS_table, "units", &units1, SDDS_GET_BY_NAME, name)!=SDDS_STRING) {
    SDDS_SetError("units field of column has wrong data type!");
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
  }
  if (!units) {
    if (!units1)
      return(1);
    if (SDDS_StringIsBlank(units1)) {
      free(units1);
      return(1);
    }
    return(0);
  }
  if (!units1)
    return(0);
  if (strcmp(units, units1)==0) {
    free(units1);
    return(1);
  }
  free(units1);
  return(0);
}

void readMomentsFile(char *filename, double moment[6][6], double centroid[6])
{
  SDDS_DATASET SDDSin;
  long i, j, rows;
  double *data;
  char s[100];

  if (!SDDS_InitializeInput(&SDDSin, filename))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  if (SDDS_ReadPage(&SDDSin)<=0 || (rows=SDDS_RowCount(&SDDSin))<=0)
    SDDS_Bomb("Problem reading moments file. Is it empty?");

  for (i=0; i<6; i++) {
    /* Get s[i] */
    sprintf(s, "s%ld", i+1);
    if (!(data=SDDS_GetColumnInDoubles(&SDDSin, s))) 
      fprintf(stderr, "Error: problem reading %s from moments file. Check existence and type.\n", s);
    moment[i][i] = sqr(data[rows-1]);
    free(data);

    /* Get c[i] */
    sprintf(s, "c%ld", i+1);
    if (!(data=SDDS_GetColumnInDoubles(&SDDSin, s))) 
      fprintf(stderr, "Error: problem reading %s from moments file. Check existence and type.\n", s);
    centroid[i] = data[rows-1];
    free(data);

    for (j=i+1; j<6; j++) {
      /* Get s[i][j] */
      sprintf(s, "s%ld%ld", i+1, j+1);
      if (!(data=SDDS_GetColumnInDoubles(&SDDSin, s))) 
        fprintf(stderr, "Error: problem reading %s from moments file. Check existence and type.\n", s);
      moment[i][j] = moment[j][i] = data[rows-1];
      free(data);
    }
  }

  if (!SDDS_Terminate(&SDDSin))
    SDDS_Bomb("Problem terminating moments file.");
}

void transformCoordinates(double *x, double *xp, double *y, double *yp, double *t, double *p, long np,
                          unsigned long flags, double desiredMoment[6][6], double desiredCentroid[6])
{
  long i, j, ip;
  double pAve = 0;

  if (flags&(FL_EXCLUDE_X+FL_EXCLUDE_Y+FL_EXCLUDE_Z)) {
    /* Do a 4x4 transformation */
    double **coord, sigma[4][4], dMoment[4][4], dCentroid[4], centroid[4], M[4][4];

    if (!(flags&FL_EXCLUDE_Z)) {
      /* transform (t, p) to (s, delta) */
      for (ip=0; ip<np; ip++)
        pAve += p[ip];
      pAve /= np;
      for (ip=0; ip<np; ip++) {
        double beta, gamma;
        gamma = sqrt(sqr(p[ip])-1);
        beta = p[ip]/gamma;
        t[ip] /= beta*c_mks;
        p[ip] = (p[ip]-pAve)/pAve;
      }
    }

    /* Copy data */
    coord = (double**)zarray_2d(sizeof(**coord), np, 4);
    if (flags&FL_EXCLUDE_X) {
      for (ip=0; ip<np; ip++) {
        coord[ip][0] = y[ip];
        coord[ip][1] = yp[ip];
        coord[ip][3] = t[ip];
        coord[ip][4] = p[ip];
      }
      for (i=0; i<4; i++) {
        dCentroid[i] = desiredCentroid[i+2];
        for (j=0; j<4; j++)
          dMoment[i][j] = desiredMoment[i+2][j+2];
      }
    } else if (flags&FL_EXCLUDE_Y) {
      long ii, jj;
      for (ip=0; ip<np; ip++) {
        coord[ip][0] = x[ip];
        coord[ip][1] = xp[ip];
        coord[ip][3] = t[ip];
        coord[ip][4] = p[ip];
      }
      for (i=0; i<4; i++) {
        ii = i<2 ? i : i+2;
        dCentroid[i] = desiredCentroid[ii];
        for (j=0; j<4; j++) {
          jj = j<2 ? j : j+2;
          dMoment[i][j] = desiredMoment[ii][jj];
        }
      }
    } else {
      for (ip=0; ip<np; ip++) {
        coord[ip][0] = x[ip];
        coord[ip][1] = xp[ip];
        coord[ip][2] = y[ip];
        coord[ip][3] = yp[ip];
      }
    }

    /* Compute centroids of input beam */
    for (i=0; i<4; i++) {
      centroid[i] = 0;
      for (ip=0; ip<np; ip++)
        centroid[i] += coord[ip][i];
      centroid[i] /= np;
    }

    /* Remove centroids */
    for (ip=0; ip<np; ip++) {
      for (i=0; i<4; i++)
        coord[ip][i] -= centroid[i];
    }
     
    /* Compute sigma matrix */
    for (i=0; i<4; i++) {
      for (j=0; j<=i; j++) {
        sigma[i][j] = 0;
        for (ip=0; ip<np; ip++) 
          sigma[i][j] += coord[ip][i]*coord[ip][j];
        sigma[j][i] = sigma[i][j] = sigma[i][j]/np;
      }
    }

    /* Find transform to get desired sigma matrix */
    findTransformationMatrix4(sigma, dMoment, M);
    
    /* Transform individual coordinates */
    for (ip=0; ip<np; ip++) {
      double copy[4];
      memcpy(&copy[0], &coord[ip][0], sizeof(double)*4);
      for (i=0; i<4; i++) {
        coord[ip][i] = 0;
        for (j=0; j<=i; j++)
          coord[ip][i] += M[i][j]*copy[j];
      }
    }

    if (!(flags&FL_EXCLUDE_CENTROIDS)) {
      /* Set centroids */
      for (ip=0; ip<np; ip++)
        for (i=0; i<4; i++)
          coord[ip][i] += dCentroid[i];
    } else {
      /* Restore centroids */
      for (ip=0; ip<np; ip++)
        for (i=0; i<4; i++)
          coord[ip][i] += centroid[i];
    }
    
    /* Copy data back to input arrays */
    if (flags&FL_EXCLUDE_X) {
      for (ip=0; ip<np; ip++) {
        y[ip]  = coord[ip][0];
        yp[ip] = coord[ip][1];
        t[ip]  = coord[ip][3];
        p[ip]  = coord[ip][4];
      }
    } else if (flags&FL_EXCLUDE_Y) {
      for (ip=0; ip<np; ip++) {
        x[ip]  = coord[ip][0];
        xp[ip] = coord[ip][1];
        t[ip]  = coord[ip][3];
        p[ip]  = coord[ip][4];
      }
    } else {
      for (ip=0; ip<np; ip++) {
        x[ip]  = coord[ip][0];
        xp[ip] = coord[ip][1];
        y[ip]  = coord[ip][2];
        yp[ip] = coord[ip][3];
      }
    }

    if (!(flags&FL_EXCLUDE_Z)) {
      /* transform (s, delta) to (t, p)*/
      for (ip=0; ip<np; ip++) {
        double beta, gamma;
        p[ip] = pAve*(1 + p[ip]);
        gamma = sqrt(sqr(p[ip])-1);
        beta = p[ip]/gamma;
        t[ip] = t[ip]*beta*c_mks;
      }
    }

    /* Free working memory */
    free_zarray_2d((void**)coord, np, 4);
  } else {
    SDDS_Bomb("Can't do 6x6 transformation yet");
  }

}

#include "gsl/gsl_linalg.h"
#include "gsl/gsl_blas.h"
#include "gsl/gsl_eigen.h"
#include "gsl/gsl_vector.h"

void findTransformationMatrix4(double sigma[4][4], double desiredSigma[4][4], double M[4][4])
{
  gsl_matrix *M1, *M2, *M3;
  long i, j;
  gsl_eigen_symm_workspace *ws;
  gsl_vector *vec;

  M1 = M2 = M3 = NULL;
  if (!(M1 = gsl_matrix_alloc(4, 4)) ||
      !(M2 = gsl_matrix_alloc(4, 4)) ||
      !(M3 = gsl_matrix_alloc(4, 4)))
    SDDS_Bomb("gsl_matrix_alloc failed");

  for (i=0; i<4; i++) {
    fprintf(stderr, "sigma[%ld]: ", i);
    for (j=0; j<4; j++)
      fprintf(stderr, "%le  ", sigma[i][j]);
    fprintf(stderr, "\n");
  }

  for (i=0; i<4; i++)
    for (j=0; j<4; j++)
      gsl_matrix_set(M1, i, j, sigma[i][j]);
  ws = gsl_eigen_symm_alloc(4);
  vec = gsl_vector_alloc(4);
  gsl_eigen_symm(M1, vec, ws);

  for (i=0; i<4; i++)
    for (j=0; j<4; j++)
      gsl_matrix_set(M1, i, j, sigma[i][j]);
  for (i=0; i<4; i++)
    fprintf(stderr,"eval[%ld] = %le\n", i, gsl_vector_get(vec, i));

  if (gsl_linalg_cholesky_decomp(M1))
    SDDS_Bomb("gsl_linalg_cholesky_decomp failed (4x4 case). Does sigma matrix have diagonal blocks of zeros?");
  /* Zero out the upper triangular part, which gsl_linalg_cholesky_decomp oddly doesn't touch */
  for (i=0; i<4; i++)
    for (j=i+1; j<4; j++)
      gsl_matrix_set(M1, i, j, 0.0);
  if (gsl_linalg_cholesky_invert(M1))
    SDDS_Bomb("gsl_linalg_cholesky_invert failed (4x4 case). Does sigma matrix have diagonal blocks of zeros?");

  for (i=0; i<4; i++)
    for (j=0; j<4; j++)
      gsl_matrix_set(M2, i, j, desiredSigma[i][j]);

  if (gsl_linalg_cholesky_decomp(M2))
    SDDS_Bomb("gsl_linalg_cholesky_decomp failed (4x4 case). Does sigma matrix have diagonal blocks of zeros?");

  /* Zero out the upper triangular part, which gsl_linalg_cholesky_decomp oddly doesn't touch */
  for (i=0; i<4; i++)
    for (j=i+1; j<4; j++)
      gsl_matrix_set(M2, i, j, 0.0);

  /* Compute the product M2*M1 */
  if (gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, M2, M1, 0.0, M3))
    SDDS_Bomb("gsl_blas_dgemm (4x4 case).");

  for (i=0; i<4; i++)
    for (j=0; j<4; j++)
      M[i][j] = gsl_matrix_get(M3, i, j);

  gsl_matrix_free(M1);
  gsl_matrix_free(M2);
  gsl_matrix_free(M3);
}
