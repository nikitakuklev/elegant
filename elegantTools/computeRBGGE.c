#include "mdb.h"
#include "SDDS.h"
#include "scan.h"
#include "fftpackC.h"
//#define DEBUG 1
//#define OLDFFT 1

#define TWOPI 6.28318530717958647692528676656

typedef struct COMPLEX
{
  double re;
  double im;
} COMPLEX;
COMPLEX fourierCoeffIntegralTrap(COMPLEX *Bint, double *x, double lambdaN, double dx, double xMax, int32_t Nx);
COMPLEX fourierCoeffIntegralSimp(COMPLEX *Bint, double *x, double lambdaN, double dx, double xMax, int32_t Nx);
COMPLEX fourierCoeffIntegralLinInterp(COMPLEX *Bint, double *x, double lambdaN, double dx, double xMax, int32_t Nx);
COMPLEX fourierCoeffIntegralLinInterpSkew(COMPLEX *Bint, double *x, double lambdaN, double dx, double xMax, int32_t Nx);
COMPLEX calcGGtopbottomA(COMPLEX *beta, double k, double *lambda, double yMax, int32_t grad_r, int32_t Ncoeff);
COMPLEX calcGGtopbottomB(COMPLEX *beta, double k, double *lambda, double yMax, int32_t grad_r, int32_t Ncoeff);
COMPLEX calcGGrightA(COMPLEX *beta, double k, double *tau, double xMax, int32_t grad_r, int32_t Ncoeff);
COMPLEX calcGGrightB(COMPLEX *beta, double k, double *tau, double xMax, int32_t grad_r, int32_t Ncoeff);
COMPLEX calcGGleftA(COMPLEX *beta, double k, double *tau, double xMax, int32_t grad_r, int32_t Ncoeff);
COMPLEX calcGGleftB(COMPLEX *beta, double k, double *tau, double xMax, int32_t grad_r, int32_t Ncoeff);
COMPLEX calcGGtopbottomk0A(COMPLEX *beta, double *lambda, double yMax, int32_t grad_r, int32_t Ncoeff);
COMPLEX calcGGtopbottomk0B(COMPLEX *beta, double *lambda, double yMax, int32_t grad_r, int32_t Ncoeff);
COMPLEX calcGGrightk0(COMPLEX *beta, double *tau, double xMax, int32_t grad_r, int32_t Ncoeff);
COMPLEX calcGGleftk0(COMPLEX *beta, double *tau, double xMax, int32_t grad_r, int32_t Ncoeff);
COMPLEX calcGGallSidesBz(COMPLEX *beta, double k, double *lambda, double yMax, int32_t Ncoeff);
void FFT(COMPLEX *field, int32_t isign, int32_t npts);
unsigned long IntCeilingPowerOf2(unsigned long i);

int computeGGderiv(char *topFile, char *bottomFile, char *leftFile, char *rightFile, char *outputFile, long derivatives, long multipoles, long fundamental);
int computeGGcos(char *topFile, char *bottomFile, char *leftFile, char *rightFile, char *outputFile, long derivatives, long multipoles, long fundamental);

int ReadInputFiles(long BzMode, char *topFile, char *bottomFile, char *leftFile, char *rightFile,
                   int32_t *Nx, int32_t *Ny, int32_t *Nfft,
                   double *dx, double *dy, double *dz,
                   COMPLEX ***ByTop, COMPLEX ***ByBottom, COMPLEX ***BxRight, COMPLEX ***BxLeft);

#define SET_TOP_BY 0
#define SET_BOTTOM_BY 1
#define SET_LEFT_BX 2
#define SET_RIGHT_BX 3
#define SET_NORMAL 4
#define SET_SKEW 5
#define SET_DERIVATIVES 6
#define SET_MULTIPOLES 7
#define SET_FUNDAMENTAL 8
#define N_OPTIONS 9

char *option[N_OPTIONS] = {
  "top", "bottom", "left", "right", "normal", "skew", "derivatives", "multipoles", "fundamental"};

#define USAGE "computeRBGGE -top=<filename> -bottom=<filename> -left=<filename> -right=<filename>\n\
             -normal=<output> [-skew=<output>] [-derivatives=<number>] [-multipoles=<number>] [-fundamental=<number>]\n\
-top         (x, y, z, Bx, By, Bz) map for top plane (y=constant, y>0).\n\
-bottom      (x, y, z, Bx, By, Bz) map for bottom plane (y=constant, y<0).\n\
-right       (x, y, z, Bx, By, Bz) map for right plane (x=constant, x<0).\n\
-left        (x, y, z, Bx, By, Bz) map for left plane (x=constant, x>0).\n\
-normal      Output file for normal-component generalized gradients.\n\
-skew        Output file for skew-component generalized gradients.\n\
-derivatives Number of derivatives vs z desired in output. Default: 7\n\
-multipoles  Number of multipoles desired in output. Default: 8\n\
-fundamental Fundamental multipole of sequence. 0=none (default), 1=dipole, 2=quadrupole, etc.\n\n\
Rectangular Boundary Generalized Gradient Expansion by Ryan Lindberg, Robert Soliday, and Michael Borland."

int main(int argc, char **argv)
{
  SCANNED_ARG *scanned;
  long i_arg;
  long multipoles = 8, derivatives = 7, fundamental=0;
  char *topFile = NULL, *bottomFile = NULL, *leftFile = NULL, *rightFile = NULL;
  char *normalOutputFile = NULL, *skewOutputFile = NULL;

#ifdef DEBUG
#ifdef OLDFFT
  fprintf(stderr, "Using old FFT\n");
#else
  fprintf(stderr, "Using new FFT\n");
#endif
#endif

  argc = scanargs(&scanned, argc, argv);
  if (argc < 2 || argc > (2 + N_OPTIONS))
    {
      fprintf(stderr, "%s\n", USAGE);
      return (1);
    }
  for (i_arg = 1; i_arg < argc; i_arg++)
    {
      if (scanned[i_arg].arg_type == OPTION)
        {
          /* process options here */
          switch (match_string(scanned[i_arg].list[0], option, N_OPTIONS, 0))
            {
            case SET_TOP_BY:
              if (scanned[i_arg].n_items != 2)
                {
                  fprintf(stderr, "invalid -top syntax\n%s\n", USAGE);
                  return (1);
                }
              topFile = scanned[i_arg].list[1];
              break;
            case SET_BOTTOM_BY:
              if (scanned[i_arg].n_items != 2)
                {
                  fprintf(stderr, "invalid -bottom syntax\n%s\n", USAGE);
                  return (1);
                }
              bottomFile = scanned[i_arg].list[1];
              break;
            case SET_LEFT_BX:
              if (scanned[i_arg].n_items != 2)
                {
                  fprintf(stderr, "invalid -left syntax\n%s\n", USAGE);
                  return (1);
                }
              leftFile = scanned[i_arg].list[1];
              break;
            case SET_RIGHT_BX:
              if (scanned[i_arg].n_items != 2)
                {
                  fprintf(stderr, "invalid -right syntax\n%s\n", USAGE);
                  return (1);
                }
              rightFile = scanned[i_arg].list[1];
              break;
            case SET_NORMAL:
              if (scanned[i_arg].n_items != 2)
                {
                  fprintf(stderr, "invalid -normal syntax\n%s\n", USAGE);
                  return (1);
                }
              normalOutputFile = scanned[i_arg].list[1];
              break;
            case SET_SKEW:
              if (scanned[i_arg].n_items != 2)
                {
                  fprintf(stderr, "invalid -skew syntax\n%s\n", USAGE);
                  return (1);
                }
              skewOutputFile = scanned[i_arg].list[1];
              break;
            case SET_DERIVATIVES:
              if (scanned[i_arg].n_items != 2 ||
                  sscanf(scanned[i_arg].list[1], "%ld", &derivatives) != 1 ||
                  derivatives <= 0)
                {
                  fprintf(stderr, "invalid -derivatives syntax\n%s\n", USAGE);
                  return (1);
                }
              break;
            case SET_MULTIPOLES:
              if (scanned[i_arg].n_items != 2 ||
                  sscanf(scanned[i_arg].list[1], "%ld", &multipoles) != 1 ||
                  multipoles <= 0)
                {
                  fprintf(stderr, "invalid -multipoles syntax\n%s\n", USAGE);
                  return (1);
                }
              break;
            case SET_FUNDAMENTAL:
              if (scanned[i_arg].n_items != 2 ||
                  sscanf(scanned[i_arg].list[1], "%ld", &fundamental) != 1 ||
                  fundamental < 0)
                {
                  fprintf(stderr, "invalid -fundamental syntax\n%s\n", USAGE);
                  return (1);
                }
              break;
            default:
              fprintf(stderr, "unknown option given\n%s\n", USAGE);
              return (1);
              break;
            }
        }
      else
        {
          fprintf(stderr, "too many files listed\n%s\n", USAGE);
          return (1);
        }
    }
  if ((topFile == NULL) || (bottomFile == NULL) || (leftFile == NULL) || (rightFile == NULL))
    {
      fprintf(stderr, "%s\n", USAGE);
      return (1);
    }
#ifdef DEBUG
    fprintf(stderr, "topFile=%s\n", topFile);
    fprintf(stderr, "bottomFile=%s\n", bottomFile);
    fprintf(stderr, "leftFile=%s\n", leftFile);
    fprintf(stderr, "rightFile=%s\n", rightFile);
    fprintf(stderr, "derivatives=%ld\n", derivatives);
    fprintf(stderr, "multipoles=%ld\n", multipoles);
    fprintf(stderr, "normalOutputFile=%s\n", normalOutputFile);
    fprintf(stderr, "skewOutputFile=%s\n", skewOutputFile);
#endif
  if (normalOutputFile != NULL)
    {
      if (computeGGderiv(topFile, bottomFile, leftFile, rightFile, normalOutputFile, derivatives, multipoles, fundamental))
        return 1;
    }
  if (skewOutputFile != NULL)
    {
      computeGGcos(topFile, bottomFile, leftFile, rightFile, skewOutputFile, derivatives, multipoles, fundamental);
    }
  return (0);
}

int ReadInputFiles(long BzMode, char *topFile, char *bottomFile, char *leftFile, char *rightFile,
                   int32_t *Nx, int32_t *Ny, int32_t *Nfft,
                   double *dx, double *dy, double *dz,
                   COMPLEX ***ByTop, COMPLEX ***ByBottom, COMPLEX ***BxRight, COMPLEX ***BxLeft)
{
  SDDS_DATASET SDDSInput;
  double *cvalues, *xvalues, *yvalues, *zvalues;
  int32_t rows, ik, ix, iy, n;
  int32_t tmpNx, tmpNy, tmpNfft;
  double tmpdx, tmpdy, tmpdz;
  char name[20];

  /* Read in By from the top face */
  if (SDDS_InitializeInput(&SDDSInput, topFile) != 1)
    {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      return (1);
    }
  if (BzMode)
    {
      sprintf(name, "Bz");
    }
  else
    {
      sprintf(name, "By");
    }
  if ((SDDS_CheckColumn(&SDDSInput, name, NULL, SDDS_ANY_NUMERIC_TYPE, stderr) != SDDS_CHECK_OKAY) ||
      (SDDS_CheckColumn(&SDDSInput, "x", NULL, SDDS_ANY_NUMERIC_TYPE, stderr) != SDDS_CHECK_OKAY) ||
      (SDDS_CheckColumn(&SDDSInput, "z", NULL, SDDS_ANY_NUMERIC_TYPE, stderr) != SDDS_CHECK_OKAY))
    {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      return (1);
    }
  if (SDDS_ReadPage(&SDDSInput) != 1)
    {
      fprintf(stderr, "Unable to read SDDS page\n");
      return (1);
    }
  rows = SDDS_RowCount(&SDDSInput);
  cvalues = SDDS_GetColumnInDoubles(&SDDSInput, name);
  if (cvalues == NULL)
    {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      return (1);
    }
  xvalues = SDDS_GetColumnInDoubles(&SDDSInput, "x");
  if (xvalues == NULL)
    {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      return (1);
    }
  zvalues = SDDS_GetColumnInDoubles(&SDDSInput, "z");
  if (zvalues == NULL)
    {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      return (1);
    }
  if (SDDS_Terminate(&SDDSInput) != 1)
    {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      return (1);
    }
  if ((zvalues[0] != zvalues[1]) || (xvalues[0] == xvalues[1])) {
      fprintf(stderr, "Error: Sort the input files with sddssort -col=z -col=y -col=x\n");
      return (1);
  }
  for (ix = 1; ix < rows; ix++)
    {
      if (zvalues[ix-1] != zvalues[ix])
        {
          *Nx = ix;
          *Nfft = rows / ix;
          *dx = (xvalues[ix-1] - xvalues[0]) / (ix - 1);
          *dz = (zvalues[rows-1] - zvalues[0]) / (*Nfft - 1);
          break;
        }
    }
  if (rows != *Nx * *Nfft)
    {
      fprintf(stderr, "Unexpected row count\n");
      return (1);
    }
#ifdef DEBUG
  fprintf(stderr, "Top file %s, Nx=%ld, Nfft=%ld, dx=%le, dz=%le\n",
          topFile, (long)*Nx, (long)*Nfft, *dx, *dz);
#endif

  *ByTop = calloc(*Nx, sizeof(COMPLEX *));
  for (ix = 0; ix < *Nx; ix++)
    {
      (*ByTop)[ix] = calloc(*Nfft, sizeof(COMPLEX));
    }
  n = 0;
  for (ik = 0; ik < *Nfft; ik++)
    {
      for (ix = 0; ix < *Nx; ix++)
        {
          (*ByTop)[ix][ik].re = cvalues[n];
          n++;
        }
    }
  free(cvalues);
  free(xvalues);
  free(zvalues);
  tmpNx = *Nx;
  tmpNfft = *Nfft;
  tmpdx = *dx;
  tmpdz = *dz;
  /* Read in By from the bottom face */
  if (SDDS_InitializeInput(&SDDSInput, bottomFile) != 1)
    {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      return (1);
    }
  if (BzMode)
    {
      sprintf(name, "Bz");
    }
  else
    {
      sprintf(name, "By");
    }
  if ((SDDS_CheckColumn(&SDDSInput, name, NULL, SDDS_ANY_NUMERIC_TYPE, stderr) != SDDS_CHECK_OKAY) ||
      (SDDS_CheckColumn(&SDDSInput, "x", NULL, SDDS_ANY_NUMERIC_TYPE, stderr) != SDDS_CHECK_OKAY) ||
      (SDDS_CheckColumn(&SDDSInput, "z", NULL, SDDS_ANY_NUMERIC_TYPE, stderr) != SDDS_CHECK_OKAY))
    {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      return (1);
    }
  if (SDDS_ReadPage(&SDDSInput) != 1)
    {
      fprintf(stderr, "Unable to read SDDS page\n");
      return (1);
    }
  rows = SDDS_RowCount(&SDDSInput);
  cvalues = SDDS_GetColumnInDoubles(&SDDSInput, name);
  if (cvalues == NULL)
    {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      return (1);
    }
  xvalues = SDDS_GetColumnInDoubles(&SDDSInput, "x");
  if (xvalues == NULL)
    {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      return (1);
    }
  zvalues = SDDS_GetColumnInDoubles(&SDDSInput, "z");
  if (zvalues == NULL)
    {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      return (1);
    }
  if (SDDS_Terminate(&SDDSInput) != 1)
    {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      return (1);
    }
  if ((zvalues[0] != zvalues[1]) || (xvalues[0] == xvalues[1])) {
      fprintf(stderr, "Error: Sort the input files with sddssort -col=z -col=y -col=x\n");
      return (1);
  }
  for (ix = 1; ix < rows; ix++)
    {
      if (zvalues[ix-1] != zvalues[ix])
        {
          *Nx = ix;
          *Nfft = rows / ix;
          *dx = (xvalues[ix-1] - xvalues[0]) / (ix - 1);
          *dz = (zvalues[rows-1] - zvalues[0]) / (*Nfft - 1);
          break;
        }
    }
  if (rows != *Nx * *Nfft)
    {
      fprintf(stderr, "Unexpected row count\n");
      return (1);
    }
#ifdef DEBUG
  fprintf(stderr, "Bottom file %s, Nx=%ld, Nfft=%ld, dx=%le, dz=%le\n",
          bottomFile, (long)*Nx, (long)*Nfft, *dx, *dz);
#endif

  if (tmpNx != *Nx)
    {
      fprintf(stderr, "Nx values differ in the input files\n");
      return (1);
    }
  if (tmpNfft != *Nfft)
    {
      fprintf(stderr, "Nfft values differ in the input files\n");
      return (1);
    }
  if ((tmpdx + 1e-9 < *dx) || (tmpdx - 1e-9 > *dx))
    {
      fprintf(stderr, "dx values differ in the input files\n");
      return (1);
    }
  if ((tmpdz + 1e-9 < *dz) || (tmpdz - 1e-9 > *dz))
    {
      fprintf(stderr, "dz values differ in the input files\n");
      return (1);
    }

  *ByBottom = calloc(*Nx, sizeof(COMPLEX *));
  for (ix = 0; ix < *Nx; ix++)
    {
      (*ByBottom)[ix] = calloc(*Nfft, sizeof(COMPLEX));
    }
  n = 0;
  for (ik = 0; ik < *Nfft; ik++)
    {
      for (ix = 0; ix < *Nx; ix++)
        {
          (*ByBottom)[ix][ik].re = cvalues[n];
          n++;
        }
    }
  free(cvalues);
  free(xvalues);
  free(zvalues);

  /* Read in Bx from the left face */
  if (SDDS_InitializeInput(&SDDSInput, leftFile) != 1)
    {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      return (1);
    }
  if (BzMode)
    {
      sprintf(name, "Bz");
    }
  else
    {
      sprintf(name, "Bx");
    }
  if ((SDDS_CheckColumn(&SDDSInput, name, NULL, SDDS_ANY_NUMERIC_TYPE, stderr) != SDDS_CHECK_OKAY) ||
      (SDDS_CheckColumn(&SDDSInput, "y", NULL, SDDS_ANY_NUMERIC_TYPE, stderr) != SDDS_CHECK_OKAY) ||
      (SDDS_CheckColumn(&SDDSInput, "z", NULL, SDDS_ANY_NUMERIC_TYPE, stderr) != SDDS_CHECK_OKAY))
    {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      return (1);
    }
  if (SDDS_ReadPage(&SDDSInput) != 1)
    {
      fprintf(stderr, "Unable to read SDDS page\n");
      return (1);
    }
  rows = SDDS_RowCount(&SDDSInput);
  cvalues = SDDS_GetColumnInDoubles(&SDDSInput, name);
  if (cvalues == NULL)
    {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      return (1);
    }
  yvalues = SDDS_GetColumnInDoubles(&SDDSInput, "y");
  if (yvalues == NULL)
    {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      return (1);
    }
  zvalues = SDDS_GetColumnInDoubles(&SDDSInput, "z");
  if (zvalues == NULL)
    {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      return (1);
    }
  if (SDDS_Terminate(&SDDSInput) != 1)
    {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      return (1);
    }
  if ((zvalues[0] != zvalues[1]) || (yvalues[0] == yvalues[1])) {
      fprintf(stderr, "Error: Sort the input files with sddssort -col=z -col=y -col=x\n");
      return (1);
  }
  for (iy = 1; iy < rows; iy++)
    {
      if (zvalues[iy-1] != zvalues[iy])
        {
          *Ny = iy;
          *Nfft = rows / iy;
          *dy = (yvalues[iy-1] - yvalues[0]) / (iy - 1);
          *dz = (zvalues[rows-1] - zvalues[0]) / (*Nfft - 1);
          break;
        }
    }
  if (rows != *Ny * *Nfft)
    {
      fprintf(stderr, "Unexpected row count\n");
      return (1);
    }

  if (tmpNfft != *Nfft)
    {
      fprintf(stderr, "Nfft values differ in the input files\n");
      return (1);
    }
  if ((tmpdz + 1e-9 < *dz) || (tmpdz - 1e-9 > *dz))
    {
      fprintf(stderr, "dz values differ in the input files\n");
      return (1);
    }
#ifdef DEBUG
  fprintf(stderr, "Left file %s, Ny=%ld, Nfft=%ld, dx=%le, dz=%le\n",
          leftFile, (long)*Ny, (long)*Nfft, *dx, *dz);
#endif

  *BxLeft = calloc(*Ny, sizeof(COMPLEX *));
  for (iy = 0; iy < *Ny; iy++)
    {
      (*BxLeft)[iy] = calloc(*Nfft, sizeof(COMPLEX));
    }
  n = 0;
  for (ik = 0; ik < *Nfft; ik++)
    {
      for (iy = 0; iy < *Ny; iy++)
        {
          (*BxLeft)[iy][ik].re = cvalues[n];
          n++;
        }
    }
  free(cvalues);
  free(yvalues);
  free(zvalues);
  tmpNy = *Ny;
  tmpdy = *dy;

  /* Read in Bx from the right face */
  if (SDDS_InitializeInput(&SDDSInput, rightFile) != 1)
    {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      return (1);
    }
  if (BzMode)
    {
      sprintf(name, "Bz");
    }
  else
    {
      sprintf(name, "Bx");
    }
  if ((SDDS_CheckColumn(&SDDSInput, name, NULL, SDDS_ANY_NUMERIC_TYPE, stderr) != SDDS_CHECK_OKAY) ||
      (SDDS_CheckColumn(&SDDSInput, "y", NULL, SDDS_ANY_NUMERIC_TYPE, stderr) != SDDS_CHECK_OKAY) ||
      (SDDS_CheckColumn(&SDDSInput, "z", NULL, SDDS_ANY_NUMERIC_TYPE, stderr) != SDDS_CHECK_OKAY))
    {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      return (1);
    }
  if (SDDS_ReadPage(&SDDSInput) != 1)
    {
      fprintf(stderr, "Unable to read SDDS page\n");
      return (1);
    }
  rows = SDDS_RowCount(&SDDSInput);
  cvalues = SDDS_GetColumnInDoubles(&SDDSInput, name);
  if (cvalues == NULL)
    {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      return (1);
    }
  yvalues = SDDS_GetColumnInDoubles(&SDDSInput, "y");
  if (xvalues == NULL)
    {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      return (1);
    }
  zvalues = SDDS_GetColumnInDoubles(&SDDSInput, "z");
  if (zvalues == NULL)
    {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      return (1);
    }
  if (SDDS_Terminate(&SDDSInput) != 1)
    {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      return (1);
    }
  if ((zvalues[0] != zvalues[1]) || (yvalues[0] == yvalues[1])) {
      fprintf(stderr, "Error: Sort the input files with sddssort -col=z -col=y -col=x\n");
      return (1);
  }
  for (iy = 1; iy < rows; iy++)
    {
      if (zvalues[iy-1] != zvalues[iy])
        {
          *Ny = iy;
          *Nfft = rows / iy;
          *dy = (yvalues[iy-1] - yvalues[0]) / (iy - 1);
          *dz = (zvalues[rows-1] - zvalues[0]) / (*Nfft - 1);
          break;
        }
    }
  if (rows != *Ny * *Nfft)
    {
      fprintf(stderr, "Unexpected row count\n");
      return (1);
    }
#ifdef DEBUG
  fprintf(stderr, "Right file %s, Ny=%ld, Nfft=%ld, dx=%le, dz=%le\n",
          rightFile, (long)*Ny, (long)*Nfft, *dx, *dz);
#endif

  if (tmpNy != *Ny)
    {
      fprintf(stderr, "Ny values differ in the input files\n");
      return (1);
    }
  if (tmpNfft != *Nfft)
    {
      fprintf(stderr, "Nfft values differ in the input files\n");
      return (1);
    }
  if ((tmpdy + 1e-9 < *dy) || (tmpdx - 1e-9 > *dy))
    {
      fprintf(stderr, "dx values differ in the input files\n");
      return (1);
    }
  if ((tmpdz + 1e-9 < *dz) || (tmpdz - 1e-9 > *dz))
    {
      fprintf(stderr, "dz values differ in the input files\n");
      return (1);
    }

  *BxRight = calloc(*Ny, sizeof(COMPLEX *));
  for (iy = 0; iy < *Ny; iy++)
    {
      (*BxRight)[iy] = calloc(*Nfft, sizeof(COMPLEX));
    }
  n = 0;
  for (ik = 0; ik < *Nfft; ik++)
    {
      for (iy = 0; iy < *Ny; iy++)
        {
          (*BxRight)[iy][ik].re = cvalues[n];
          n++;
        }
    }
  free(cvalues);
  free(yvalues);
  free(zvalues);

  return (0);
}

int computeGGderiv(char *topFile, char *bottomFile, char *leftFile, char *rightFile, char *outputFile, long derivatives, long multipoles, long fundamental)
{
  COMPLEX **ByTop, **ByBottom, **BxRight, **BxLeft;

  COMPLEX **betaTop, **betaBottom, **betaRight, **betaLeft;
  COMPLEX **genGradr_k;
  COMPLEX **derivGG;
  COMPLEX *Bint;
  COMPLEX genGradT, genGradB, genGradR, genGradL;

  double *lambda, *tau, *k, *x, *y;

  double xMax, yMax;
  double dx, dy, dz, dk, invNfft;

  int32_t n, ir, ix, Nx, iy, Ny, ik, Nfft, Nz;

  int32_t Ngrad = 8;
  int32_t Nderiv = 7;
  int32_t Ncoeff = 40;

  SDDS_DATASET SDDSOutput;
  char name[20];

  Ngrad = multipoles;
  Nderiv = 2 * derivatives - 1;

  if (ReadInputFiles(0, topFile, bottomFile, leftFile, rightFile, &Nx, &Ny, &Nfft, &dx, &dy, &dz, &ByTop, &ByBottom, &BxRight, &BxLeft) != 0)
    {
      return (1);
    }
  Nz = Nfft;
  
  dk = TWOPI / (dz * (double)Nfft);
  for (ix = 0; ix < Nx; ix++)
    {
      FFT(ByTop[ix], -1, Nfft);
      FFT(ByBottom[ix], -1, Nfft);
    }
  for (iy = 0; iy < Ny; iy++)
    {
      FFT(BxRight[iy], -1, Nfft);
      FFT(BxLeft[iy], -1, Nfft);
    }

  xMax = dx * 0.5 * (double)(Nx - 1);
  yMax = dy * 0.5 * (double)(Ny - 1);
#ifdef DEBUG
  fprintf(stderr, "xMax = %e, yMax=%e\n", xMax, yMax);
#endif
  x = calloc(Nx, sizeof(double));
  for (ix = 0; ix < Nx; ix++)
    x[ix] = -xMax + dx * (double)ix;
  y = calloc(Ny, sizeof(double));
  for (iy = 0; iy < Ny; iy++)
    y[iy] = -yMax + dy * (double)iy;
  lambda = calloc(Ncoeff, sizeof(double));
  tau = calloc(Ncoeff, sizeof(double));
  for (n = 0; n < Ncoeff; n++)
    {
      lambda[n] = PI * (double)n / (2.0 * xMax);
      tau[n] = PI * (double)n / (2.0 * yMax);
    }
  betaTop = calloc(Nfft, sizeof(COMPLEX *));
  betaBottom = calloc(Nfft, sizeof(COMPLEX *));
  betaRight = calloc(Nfft, sizeof(COMPLEX *));
  betaLeft = calloc(Nfft, sizeof(COMPLEX *));
  for (ik = 0; ik < Nfft; ik++)
    {
      betaTop[ik] = calloc(Ncoeff, sizeof(COMPLEX));
      betaBottom[ik] = calloc(Ncoeff, sizeof(COMPLEX));
      betaRight[ik] = calloc(Ncoeff, sizeof(COMPLEX));
      betaLeft[ik] = calloc(Ncoeff, sizeof(COMPLEX));
    }

  Bint = calloc(Nx, sizeof(COMPLEX));
  for (ik = 0; ik < Nfft; ik++)
    {
      for (ix = 0; ix < Nx; ix++)
        {
          Bint[ix].re = ByTop[ix][ik].re;
          Bint[ix].im = ByTop[ix][ik].im;
        }
      betaTop[ik][0] = fourierCoeffIntegralTrap(Bint, x, lambda[0], dx, xMax, Nx);
      betaTop[ik][0].re = 0.5 * betaTop[ik][0].re;
      betaTop[ik][0].im = 0.5 * betaTop[ik][0].im;
      for (n = 1; n < Ncoeff; n++)
        betaTop[ik][n] = fourierCoeffIntegralLinInterp(Bint, x, lambda[n], dx, xMax, Nx);

      for (ix = 0; ix < Nx; ix++)
        {
          Bint[ix].re = -ByBottom[ix][ik].re;
          Bint[ix].im = -ByBottom[ix][ik].im;
        }
      betaBottom[ik][0] = fourierCoeffIntegralTrap(Bint, x, lambda[0], dx, xMax, Nx);
      betaBottom[ik][0].re = 0.5 * betaBottom[ik][0].re;
      betaBottom[ik][0].im = 0.5 * betaBottom[ik][0].im;
      for (n = 1; n < Ncoeff; n++)
        betaBottom[ik][n] = fourierCoeffIntegralLinInterp(Bint, x, lambda[n], dx, xMax, Nx);
    }
  free(Bint);

  Bint = calloc(Ny, sizeof(COMPLEX));
  for (ik = 0; ik < Nfft; ik++)
    {
      for (iy = 0; iy < Ny; iy++)
        {
          Bint[iy].re = BxRight[iy][ik].re;
          Bint[iy].im = BxRight[iy][ik].im;
        }
      betaRight[ik][0] = fourierCoeffIntegralSimp(Bint, y, tau[0], dy, yMax, Ny);
      betaRight[ik][0].re = 0.5 * betaRight[ik][0].re;
      betaRight[ik][0].im = 0.5 * betaRight[ik][0].im;
      for (n = 1; n < Ncoeff; n++)
        betaRight[ik][n] = fourierCoeffIntegralLinInterp(Bint, y, tau[n], dy, yMax, Ny);

      for (iy = 0; iy < Ny; iy++)
        {
          Bint[iy].re = -BxLeft[iy][ik].re;
          Bint[iy].im = -BxLeft[iy][ik].im;
        }
      betaLeft[ik][0] = fourierCoeffIntegralSimp(Bint, y, tau[0], dy, yMax, Ny);
      betaLeft[ik][0].re = 0.5 * betaLeft[ik][0].re;
      betaLeft[ik][0].im = 0.5 * betaLeft[ik][0].im;
      for (n = 1; n < Ncoeff; n++)
        betaLeft[ik][n] = fourierCoeffIntegralLinInterp(Bint, y, tau[n], dy, yMax, Ny);
    }
  free(Bint);

  k = calloc(Nfft, sizeof(double));
  for (ik = 0; ik < Nfft / 2; ik++)
    k[ik] = dk * (double)ik;
  for (ik = Nfft / 2; ik < Nfft; ik++)
    k[ik] = -dk * (double)(Nfft - ik);
  genGradr_k = calloc(Ngrad, sizeof(COMPLEX *));
  for (ir = 0; ir < Ngrad; ir++)
    genGradr_k[ir] = calloc(Nfft, sizeof(COMPLEX));
  for (ir = 0; ir < Ngrad; ir++)
    {
      long ir1;
      if (fundamental)
        ir1 = fundamental*(2*ir+1);
      else
        ir1 = ir + 1;
      ik = 0;
      /* top and bottom need care for k->0 */
      genGradT = calcGGtopbottomk0A(betaTop[ik], lambda, yMax, ir1, Ncoeff);
      genGradB = calcGGtopbottomk0A(betaBottom[ik], lambda, yMax, ir1, Ncoeff);

      genGradR = calcGGrightA(betaRight[ik], k[ik], tau, xMax, ir1, Ncoeff);
      genGradL = calcGGleftA(betaLeft[ik], k[ik], tau, xMax, ir1, Ncoeff);

      genGradr_k[ir][ik].re = genGradT.re - genGradB.re + genGradR.re + genGradL.re;
      genGradr_k[ir][ik].im = genGradT.im - genGradB.im + genGradR.im + genGradL.im;
      for (ik = 1; ik < Nfft; ik++)
        {
          genGradT = calcGGtopbottomA(betaTop[ik], k[ik], lambda, yMax, ir1, Ncoeff);
          genGradB = calcGGtopbottomA(betaBottom[ik], k[ik], lambda, yMax, ir1, Ncoeff);
          genGradR = calcGGrightA(betaRight[ik], k[ik], tau, xMax, ir1, Ncoeff);
          genGradL = calcGGleftA(betaLeft[ik], k[ik], tau, xMax, ir1, Ncoeff);

          genGradr_k[ir][ik].re = genGradT.re - genGradB.re + genGradR.re + genGradL.re;
          genGradr_k[ir][ik].im = genGradT.im - genGradB.im + genGradR.im + genGradL.im;
        }
    }

  invNfft = 1.0 / (double)Nfft;
  derivGG = calloc(Nderiv, sizeof(COMPLEX *));
  for (n = 0; n < Nderiv; n++)
    derivGG[n] = calloc(Nfft, sizeof(COMPLEX));

#ifdef DEBUG
  fprintf(stderr, "Printing results...\n");
#endif

  if (SDDS_InitializeOutput(&SDDSOutput, SDDS_BINARY, 1, NULL, "computeRBGGE normal output", outputFile) != 1)
    {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      return (1);
    }
  if ((SDDS_DefineSimpleParameter(&SDDSOutput, "m", NULL, SDDS_LONG) != 1) ||
      (SDDS_DefineSimpleColumn(&SDDSOutput, "z", NULL, SDDS_DOUBLE) != 1))
    {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      return (1);
    }
  for (n = 0; n < Nderiv; n+=2)
    {
      sprintf(name, "CnmS%" PRId32, n);
      if (SDDS_DefineSimpleColumn(&SDDSOutput, name, NULL, SDDS_DOUBLE) != 1)
        {
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
          return (1);
        }
    }
  for (n = 0; n < Nderiv; n+=2)
    {
      sprintf(name, "dCnmS%" PRId32 "/dz", n);
      if (SDDS_DefineSimpleColumn(&SDDSOutput, name, NULL, SDDS_DOUBLE) != 1)
        {
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
          return (1);
        }
    }
  if (SDDS_WriteLayout(&SDDSOutput) != 1)
    {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      return (1);
    }
  for (ir = 0; ir < Ngrad; ir++)
    {
      /* Take derivatives */
      for (ik = 0; ik < Nfft; ik++)
        {
          derivGG[0][ik].re = -k[ik] * genGradr_k[ir][ik].im;
          derivGG[0][ik].im = k[ik] * genGradr_k[ir][ik].re;
        }
      for (n = 1; n < Nderiv; n++)
        for (ik = 0; ik < Nfft; ik++)
          {
            derivGG[n][ik].re = -k[ik] * derivGG[n - 1][ik].im;
            derivGG[n][ik].im = k[ik] * derivGG[n - 1][ik].re;
          }
      FFT(genGradr_k[ir], 1, Nfft);
      for (n = 0; n < Nderiv; n++)
        FFT(derivGG[n], 1, Nfft);
      if (SDDS_StartPage(&SDDSOutput, Nz) != 1)
        {
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
          return (1);
        }
      if (SDDS_SetParameters(&SDDSOutput, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE, "m", 
                             (int32_t)(fundamental?fundamental*(2*ir+1):ir), NULL) != 1)
        {
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
          return (1);
        }
      for (ik = 0; ik < Nz; ik++)
        {
          if (SDDS_SetRowValues(&SDDSOutput, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE, ik,
                                "z", dz * (double)ik,
                                "CnmS0", genGradr_k[ir][ik].re * invNfft,
                                NULL) != 1)
            {
              SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
              return (1);
            }
          for (n = 2; n < Nderiv; n+=2)
            {
              sprintf(name, "CnmS%" PRId32, n);
               if (SDDS_SetRowValues(&SDDSOutput, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE, ik,
                                    name, derivGG[n-1][ik].re * invNfft,
                                    NULL) != 1)
                {
                  SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
                  return (1);
                }
            }
          for (n = 0; n < Nderiv; n+=2)
            {
              sprintf(name, "dCnmS%" PRId32 "/dz", n);
               if (SDDS_SetRowValues(&SDDSOutput, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE, ik,
                                    name, derivGG[n][ik].re * invNfft,
                                    NULL) != 1)
                {
                  SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
                  return (1);
                }
            }
        }
      if (SDDS_WritePage(&SDDSOutput) != 1)
        {
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
          return (1);
        }
    }
  if (SDDS_Terminate(&SDDSOutput) != 1)
    {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      return (1);
    }

  return (0);
}

int computeGGcos(char *topFile, char *bottomFile, char *leftFile, char *rightFile, char *outputFile, long derivatives, long multipoles,
                 long fundamental)
{
  COMPLEX **ByTop, **ByBottom, **BxRight, **BxLeft;
  COMPLEX **BzTop, **BzBottom, **BzRight, **BzLeft;

  COMPLEX **betaTop, **betaBottom, **betaRight, **betaLeft;
  COMPLEX **genGradr_k;
  COMPLEX **derivGG;
  COMPLEX *Bint;
  COMPLEX genGradT, genGradB, genGradR, genGradL;

  double *lambda, *tau, *k, *x, *y;

  double xMax, yMax;
  double dx, dy, dz, dk, invNfft;

  int32_t n, ir, ix, Nx, iy, Ny, ik, Nfft, Nz;

  int32_t Ngrad = 8;
  int32_t Nderiv = 7;
  int32_t Ncoeff = 40;

  SDDS_DATASET SDDSOutput;

  char name[20];

  Ngrad = multipoles;
  Nderiv = 2 * derivatives - 1;

  if (ReadInputFiles(0, topFile, bottomFile, leftFile, rightFile, &Nx, &Ny, &Nfft, &dx, &dy, &dz, &ByTop, &ByBottom, &BxRight, &BxLeft) != 0)
    {
      return (1);
    }
  Nz = Nfft;

  dk = TWOPI / (dz * (double)Nfft);
  for (ix = 0; ix < Nx; ix++)
    {
      FFT(ByTop[ix], -1, Nfft);
      FFT(ByBottom[ix], -1, Nfft);
    }
  for (iy = 0; iy < Ny; iy++)
    {
      FFT(BxRight[iy], -1, Nfft);
      FFT(BxLeft[iy], -1, Nfft);
    }

  xMax = dx * 0.5 * (double)(Nx - 1);
  yMax = dy * 0.5 * (double)(Ny - 1);
  printf("xMax = %e, yMax=%e\n", xMax, yMax);
  x = calloc(Nx, sizeof(double));
  for (ix = 0; ix < Nx; ix++)
    x[ix] = -xMax + dx * (double)ix;
  y = calloc(Ny, sizeof(double));
  for (iy = 0; iy < Ny; iy++)
    y[iy] = -yMax + dy * (double)iy;
  lambda = calloc(Ncoeff, sizeof(double));
  tau = calloc(Ncoeff, sizeof(double));
  for (n = 0; n < Ncoeff; n++)
    {
      lambda[n] = PI * (double)n / (2.0 * xMax);
      tau[n] = PI * (double)n / (2.0 * yMax);
    }
  betaTop = calloc(Nfft, sizeof(COMPLEX *));
  betaBottom = calloc(Nfft, sizeof(COMPLEX *));
  betaRight = calloc(Nfft, sizeof(COMPLEX *));
  betaLeft = calloc(Nfft, sizeof(COMPLEX *));
  for (ik = 0; ik < Nfft; ik++)
    {
      betaTop[ik] = calloc(Ncoeff, sizeof(COMPLEX));
      betaBottom[ik] = calloc(Ncoeff, sizeof(COMPLEX));
      betaRight[ik] = calloc(Ncoeff, sizeof(COMPLEX));
      betaLeft[ik] = calloc(Ncoeff, sizeof(COMPLEX));
    }

  Bint = calloc(Nx, sizeof(COMPLEX));
  for (ik = 0; ik < Nfft; ik++)
    {
      for (ix = 0; ix < Nx; ix++)
        {
          Bint[ix].re = ByTop[ix][ik].re;
          Bint[ix].im = ByTop[ix][ik].im;
        }
      betaTop[ik][0] = fourierCoeffIntegralTrap(Bint, x, lambda[0], dx, xMax, Nx);
      betaTop[ik][0].re = 0.5 * betaTop[ik][0].re;
      betaTop[ik][0].im = 0.5 * betaTop[ik][0].im;
      for (n = 1; n < Ncoeff; n++)
        betaTop[ik][n] = fourierCoeffIntegralLinInterp(Bint, x, lambda[n], dx, xMax, Nx);

      for (ix = 0; ix < Nx; ix++)
        {
          Bint[ix].re = -ByBottom[ix][ik].re;
          Bint[ix].im = -ByBottom[ix][ik].im;
        }
      betaBottom[ik][0] = fourierCoeffIntegralTrap(Bint, x, lambda[0], dx, xMax, Nx);
      betaBottom[ik][0].re = 0.5 * betaBottom[ik][0].re;
      betaBottom[ik][0].im = 0.5 * betaBottom[ik][0].im;
      for (n = 1; n < Ncoeff; n++)
        betaBottom[ik][n] = fourierCoeffIntegralLinInterp(Bint, x, lambda[n], dx, xMax, Nx);
    }
  free(Bint);

  Bint = calloc(Ny, sizeof(COMPLEX));
  for (ik = 0; ik < Nfft; ik++)
    {
      for (iy = 0; iy < Ny; iy++)
        {
          Bint[iy].re = BxRight[iy][ik].re;
          Bint[iy].im = BxRight[iy][ik].im;
        }
      betaRight[ik][0] = fourierCoeffIntegralTrap(Bint, y, tau[0], dy, yMax, Ny);
      betaRight[ik][0].re = 0.5 * betaRight[ik][0].re;
      betaRight[ik][0].im = 0.5 * betaRight[ik][0].im;
      for (n = 1; n < Ncoeff; n++)
        betaRight[ik][n] = fourierCoeffIntegralLinInterp(Bint, y, tau[n], dy, yMax, Ny);

      for (iy = 0; iy < Ny; iy++)
        {
          Bint[iy].re = -BxLeft[iy][ik].re;
          Bint[iy].im = -BxLeft[iy][ik].im;
        }
      betaLeft[ik][0] = fourierCoeffIntegralTrap(Bint, y, tau[0], dy, yMax, Ny);
      betaLeft[ik][0].re = 0.5 * betaLeft[ik][0].re;
      betaLeft[ik][0].im = 0.5 * betaLeft[ik][0].im;
      for (n = 1; n < Ncoeff; n++)
        betaLeft[ik][n] = fourierCoeffIntegralLinInterp(Bint, y, tau[n], dy, yMax, Ny);
    }
  free(Bint);

  k = calloc(Nfft, sizeof(double));
  for (ik = 0; ik < Nfft / 2; ik++)
    k[ik] = dk * (double)ik;
  for (ik = Nfft / 2; ik < Nfft; ik++)
    k[ik] = -dk * (double)(Nfft - ik);
  genGradr_k = calloc(Ngrad, sizeof(COMPLEX *));
  for (ir = 0; ir < Ngrad; ir++)
    genGradr_k[ir] = calloc(Nfft, sizeof(COMPLEX));

  /* calculate skew gradients C_ir for ir != 0 */
  for (ir = 1; ir < Ngrad; ir++)
    {
      ik = 0;
      /* top and bottom need care for k->0 */
      genGradT = calcGGtopbottomk0B(betaTop[ik], lambda, yMax, ir, Ncoeff);
      genGradB = calcGGtopbottomk0B(betaBottom[ik], lambda, yMax, ir, Ncoeff);

      genGradR = calcGGrightk0(betaRight[ik], tau, xMax, ir, Ncoeff);
      genGradL = calcGGleftk0(betaLeft[ik], tau, xMax, ir, Ncoeff);

      genGradr_k[ir][ik].re = genGradT.re + genGradB.re + genGradR.re + genGradL.re;
      genGradr_k[ir][ik].im = genGradT.im + genGradB.im + genGradR.im + genGradL.im;
      for (ik = 1; ik < Nfft; ik++)
        {
          genGradT = calcGGtopbottomB(betaTop[ik], k[ik], lambda, yMax, ir, Ncoeff);
          genGradB = calcGGtopbottomB(betaBottom[ik], k[ik], lambda, yMax, ir, Ncoeff);
          genGradR = calcGGrightB(betaRight[ik], k[ik], tau, xMax, ir, Ncoeff);
          genGradL = calcGGleftB(betaLeft[ik], k[ik], tau, xMax, ir, Ncoeff);

          genGradr_k[ir][ik].re = genGradT.re + genGradB.re + genGradR.re + genGradL.re;
          genGradr_k[ir][ik].im = genGradT.im + genGradB.im + genGradR.im + genGradL.im;
        }
    }

  /* Get Bz on the boundary */
  if (ReadInputFiles(1, topFile, bottomFile, leftFile, rightFile, &Nx, &Ny, &Nfft, &dx, &dy, &dz, &BzTop, &BzBottom, &BzRight, &BzLeft) != 0)
    {
      return (1);
    }
  Nz = Nfft;

  for (ix = 0; ix < Nx; ix++)
    {
      FFT(BzTop[ix], -1, Nfft);
      FFT(BzBottom[ix], -1, Nfft);
    }
  for (iy = 0; iy < Ny; iy++)
    {
      FFT(BzRight[iy], -1, Nfft);
      FFT(BzLeft[iy], -1, Nfft);
    }

  /* calculate Fourier coefficients associated with Bz to get C_0 */
  Bint = calloc(Nx, sizeof(COMPLEX));
  for (ik = 0; ik < Nfft; ik++)
    {
      for (ix = 0; ix < Nx; ix++)
        {
          Bint[ix].re = BzTop[ix][ik].re;
          Bint[ix].im = BzTop[ix][ik].im;
        }
      betaTop[ik][0].re = 0.0;
      betaTop[ik][0].im = 0.0;
      for (n = 1; n < Ncoeff; n++)
        betaTop[ik][n] = fourierCoeffIntegralLinInterpSkew(Bint, x, lambda[n], dx, xMax, Nx);

      for (ix = 0; ix < Nx; ix++)
        {
          Bint[ix].re = BzBottom[ix][ik].re;
          Bint[ix].im = BzBottom[ix][ik].im;
        }
      betaBottom[ik][0].re = 0.0;
      betaBottom[ik][0].im = 0.0;
      for (n = 1; n < Ncoeff; n++)
        betaBottom[ik][n] = fourierCoeffIntegralLinInterpSkew(Bint, x, lambda[n], dx, xMax, Nx);
    }
  free(Bint);

  Bint = calloc(Ny, sizeof(COMPLEX));
  for (ik = 0; ik < Nfft; ik++)
    {
      for (iy = 0; iy < Ny; iy++)
        {
          Bint[iy].re = BzRight[iy][ik].re;
          Bint[iy].im = BzRight[iy][ik].im;
        }
      betaRight[ik][0].re = 0.0;
      betaRight[ik][0].im = 0.0;
      for (n = 1; n < Ncoeff; n++)
        betaRight[ik][n] = fourierCoeffIntegralLinInterpSkew(Bint, y, tau[n], dy, yMax, Ny);

      for (iy = 0; iy < Ny; iy++)
        {
          Bint[iy].re = BzLeft[iy][ik].re;
          Bint[iy].im = BzLeft[iy][ik].im;
        }
      betaLeft[ik][0].re = 0.0;
      betaLeft[ik][0].im = 0.0;
      for (n = 1; n < Ncoeff; n++)
        betaLeft[ik][n] = fourierCoeffIntegralLinInterpSkew(Bint, y, tau[n], dy, yMax, Ny);
    }
  free(Bint);

  for (ik = 0; ik < Nfft; ik++)
    {
      genGradT = calcGGallSidesBz(betaTop[ik], k[ik], lambda, yMax, Ncoeff);
      genGradB = calcGGallSidesBz(betaBottom[ik], k[ik], lambda, yMax, Ncoeff);
      genGradR = calcGGallSidesBz(betaRight[ik], k[ik], tau, xMax, Ncoeff);
      genGradL = calcGGallSidesBz(betaLeft[ik], k[ik], tau, xMax, Ncoeff);
      genGradr_k[0][ik].re = genGradT.re + genGradB.re + genGradR.re + genGradL.re;
      genGradr_k[0][ik].im = genGradT.im + genGradB.im + genGradR.im + genGradL.im;
    }

  invNfft = 1.0 / (double)Nfft;
  derivGG = calloc(Nderiv, sizeof(COMPLEX *));
  for (n = 0; n < Nderiv; n++)
    derivGG[n] = calloc(Nfft, sizeof(COMPLEX));

  printf("Printing results...\n");

  if (SDDS_InitializeOutput(&SDDSOutput, SDDS_BINARY, 1, NULL, "computeRBGGE skew output", outputFile) != 1)
    {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      return (1);
    }
  if ((SDDS_DefineSimpleParameter(&SDDSOutput, "m", NULL, SDDS_LONG) != 1) ||
      (SDDS_DefineSimpleColumn(&SDDSOutput, "z", NULL, SDDS_DOUBLE) != 1))
    {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      return (1);
    }
  
  for (n = 0; n < Nderiv; n+=2)
    {
      sprintf(name, "CnmC%" PRId32, n);
      if (SDDS_DefineSimpleColumn(&SDDSOutput, name, NULL, SDDS_DOUBLE) != 1)
        {
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
          return (1);
        }
    }
  for (n = 0; n < Nderiv; n+=2)
    {
      sprintf(name, "dCnmC%" PRId32 "/dz", n);
      if (SDDS_DefineSimpleColumn(&SDDSOutput, name, NULL, SDDS_DOUBLE) != 1)
        {
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
          return (1);
        }
    }
  
  if (SDDS_WriteLayout(&SDDSOutput) != 1)
    {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      return (1);
    }
  for (ir = 0; ir < Ngrad; ir++)
    {
      /* Take derivatives */
      if (ir == 0)
        for (ik = 0; ik < Nfft; ik++)
          {
            derivGG[0][ik].re = genGradr_k[ir][ik].re;
            derivGG[0][ik].im = genGradr_k[ir][ik].im;
            genGradr_k[ir][ik].re = 0.0;
            genGradr_k[ir][ik].im = 0.0;
          }
      else
        for (ik = 0; ik < Nfft; ik++)
          {
            derivGG[0][ik].re = -k[ik] * genGradr_k[ir][ik].im;
            derivGG[0][ik].im = k[ik] * genGradr_k[ir][ik].re;
          }
      for (n = 1; n < Nderiv; n++)
        for (ik = 0; ik < Nfft; ik++)
          {
            derivGG[n][ik].re = -k[ik] * derivGG[n - 1][ik].im;
            derivGG[n][ik].im = k[ik] * derivGG[n - 1][ik].re;
          }
      FFT(genGradr_k[ir], 1, Nfft);
      for (n = 0; n < Nderiv; n++)
        FFT(derivGG[n], 1, Nfft);

      if (SDDS_StartPage(&SDDSOutput, Nz) != 1)
        {
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
          return (1);
        }
      if (SDDS_SetParameters(&SDDSOutput, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE, "m", 
                             (fundamental?fundamental*(2*ir+1):ir), NULL) != 1)
        {
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
          return (1);
        }
      for (ik = 0; ik < Nz; ik++)
        {
          if (ir == 0)
            {
              if (SDDS_SetRowValues(&SDDSOutput, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE, ik,
                                    "z", dz * (double)ik,
                                    "CnmC0", 0.0,
                                    NULL) != 1)
                {
                  SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
                  return (1);
                }
            }
          else
            {
              if (SDDS_SetRowValues(&SDDSOutput, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE, ik,
                                    "z", dz * (double)ik,
                                    "CnmC0", genGradr_k[ir][ik].re * invNfft,
                                    NULL) != 1)
                {
                  SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
                  return (1);
                }
            }
          for (n = 2; n < Nderiv; n+=2)
            {
              sprintf(name, "CnmC%" PRId32, n);
              if (SDDS_SetRowValues(&SDDSOutput, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE, ik,
                                    name, derivGG[n-1][ik].re * invNfft,
                                    NULL) != 1)
                {
                  SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
                  return (1);
                }
            }
          for (n = 0; n < Nderiv; n+=2)
            {
              sprintf(name, "dCnmC%" PRId32 "/dz", n);
              if (SDDS_SetRowValues(&SDDSOutput, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE, ik,
                                    name, derivGG[n][ik].re * invNfft,
                                    NULL) != 1)
                {
                  SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
                  return (1);
                }
            }
 
        }
      if (SDDS_WritePage(&SDDSOutput) != 1)
        {
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
          return (1);
        }
    }
  if (SDDS_Terminate(&SDDSOutput) != 1)
    {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      return (1);
    }

  return (0);
}

COMPLEX fourierCoeffIntegralTrap(COMPLEX *Bint, double *x, double lambdaN, double dx, double xMax, int32_t Nx)
{
  COMPLEX integral;

  int32_t ix, NxM = Nx - 1;

  integral.re = 0.5 * Bint[0].re * cos((x[0] + xMax) * lambdaN);
  integral.im = 0.5 * Bint[0].im * cos((x[0] + xMax) * lambdaN);
  for (ix = 1; ix < NxM; ix++)
    {
      integral.re += Bint[ix].re * cos((x[ix] + xMax) * lambdaN);
      integral.im += Bint[ix].im * cos((x[ix] + xMax) * lambdaN);
    }
  integral.re += 0.5 * Bint[NxM].re * cos((x[NxM] + xMax) * lambdaN);
  integral.im += 0.5 * Bint[NxM].im * cos((x[NxM] + xMax) * lambdaN);
  integral.re = dx * integral.re / xMax;
  integral.im = dx * integral.im / xMax;

  return (integral);
}

COMPLEX fourierCoeffIntegralSimp(COMPLEX *Bint, double *x, double lambdaN, double dx, double xMax, int32_t Nx)
{
  COMPLEX integral;

  int32_t ix;

  integral.re = Bint[0].re * cos((x[0] + xMax) * lambdaN);
  integral.im = Bint[0].im * cos((x[0] + xMax) * lambdaN);
  for (ix = 1; ix < Nx - 2; ix++)
    {
      integral.re += 4.0 * Bint[ix].re * cos((x[ix] + xMax) * lambdaN);
      integral.im += 4.0 * Bint[ix].im * cos((x[ix] + xMax) * lambdaN);
      ix++;
      integral.re += 2.0 * Bint[ix].re * cos((x[ix] + xMax) * lambdaN);
      integral.im += 2.0 * Bint[ix].im * cos((x[ix] + xMax) * lambdaN);
    }
  integral.re += 4.0 * Bint[ix].re * cos((x[ix] + xMax) * lambdaN);
  integral.im += 4.0 * Bint[ix].im * cos((x[ix] + xMax) * lambdaN);
  ix++;
  integral.re += Bint[ix].re * cos((x[ix] + xMax) * lambdaN);
  integral.im += Bint[ix].im * cos((x[ix] + xMax) * lambdaN);
  integral.re = dx * integral.re / (3.0 * xMax);
  integral.im = dx * integral.im / (3.0 * xMax);

  return (integral);
}

COMPLEX fourierCoeffIntegralLinInterp(COMPLEX *Bint, double *x, double lambdaN, double dx, double xMax, int32_t Nx)
{
  COMPLEX integral;

  double slope, yint, x0, x1;

  int32_t ix;

  integral.re = 0.0;
  integral.im = 0.0;
  for (ix = 0; ix < Nx - 1; ix++)
    {
      x0 = x[ix] + xMax;
      x1 = x[ix + 1] + xMax;
      slope = (Bint[ix + 1].re - Bint[ix].re) / dx;
      yint = (x[ix + 1] * Bint[ix].re - x[ix] * Bint[ix + 1].re) / dx;
      integral.re += (yint / lambdaN) * (sin(x1 * lambdaN) - sin(x0 * lambdaN));
      integral.re += (slope / lambdaN) * (x[ix + 1] * sin(x1 * lambdaN) - x[ix] * sin(x0 * lambdaN) + (cos(x1 * lambdaN) - cos(x0 * lambdaN)) / lambdaN);
      slope = (Bint[ix + 1].im - Bint[ix].im) / dx;
      yint = (x[ix + 1] * Bint[ix].im - x[ix] * Bint[ix + 1].im) / dx;
      integral.im += (yint / lambdaN) * (sin(x1 * lambdaN) - sin(x0 * lambdaN));
      integral.im += (slope / lambdaN) * (x[ix + 1] * sin(x1 * lambdaN) - x[ix] * sin(x0 * lambdaN) + (cos(x1 * lambdaN) - cos(x0 * lambdaN)) / lambdaN);
    }
  integral.re = integral.re / xMax;
  integral.im = integral.im / xMax;

  return (integral);
}

COMPLEX fourierCoeffIntegralLinInterpSkew(COMPLEX *Bint, double *x, double lambdaN, double dx, double xMax, int32_t Nx)
{
  COMPLEX integral;

  double slope, yint, x0, x1;

  int32_t ix;

  integral.re = 0.0;
  integral.im = 0.0;
  for (ix = 0; ix < Nx - 1; ix++)
    {
      x0 = x[ix] + xMax;
      x1 = x[ix + 1] + xMax;
      slope = (Bint[ix + 1].re - Bint[ix].re) / dx;
      yint = (x[ix + 1] * Bint[ix].re - x[ix] * Bint[ix + 1].re) / dx;
      integral.re -= (yint / lambdaN) * (cos(x1 * lambdaN) - cos(x0 * lambdaN));
      integral.re += (slope / lambdaN) * (-x[ix + 1] * cos(x1 * lambdaN) + x[ix] * cos(x0 * lambdaN) + (sin(x1 * lambdaN) - sin(x0 * lambdaN)) / lambdaN);
      slope = (Bint[ix + 1].im - Bint[ix].im) / dx;
      yint = (x[ix + 1] * Bint[ix].im - x[ix] * Bint[ix + 1].im) / dx;
      integral.im -= (yint / lambdaN) * (cos(x1 * lambdaN) - cos(x0 * lambdaN));
      integral.im += (slope / lambdaN) * (-x[ix + 1] * cos(x1 * lambdaN) + x[ix] * cos(x0 * lambdaN) + (sin(x1 * lambdaN) - sin(x0 * lambdaN)) / lambdaN);
    }
  integral.re = integral.re / xMax;
  integral.im = integral.im / xMax;

  return (integral);
}

COMPLEX calcGGtopbottomA(COMPLEX *beta, double k, double *lambda, double yMax, int32_t grad_r, int32_t Ncoeff)
{
  COMPLEX genGrad;
  double klam, splus, sminus;
  double rcoef, factor;
  int32_t theta0r, deltan, sign, power;

  int32_t n, i;

  rcoef = 1.0;
  power = grad_r / 2;
  for (i = 0; i < power; i++)
    rcoef = -1.0 * rcoef; /* (-1)^{floor[r/2]} */
  for (i = grad_r; i >= 1; i--)
    rcoef *= 0.5 / (double)i; /* 1/(r!2^r) */
  if (grad_r % 2 == 0)
    theta0r = 0;
  else
    theta0r = 1;
  genGrad.re = 0.0;
  genGrad.im = 0.0;
  for (n = 0; n < Ncoeff; n++)
    {
      if (n % 2 == 0)
        deltan = theta0r;
      else
        deltan = 1 - theta0r;
      sign = 1;
      power = n / 2;
      for (i = 0; i < power; i++)
        sign = -1 * sign;
      klam = sqrt(k * k + lambda[n] * lambda[n]);
      splus = lambda[n] + klam;
      sminus = lambda[n] - klam;
      for (i = 1; i < grad_r; i++)
        {
          splus = splus * (lambda[n] + klam);   /* s_{+}^r */
          sminus = sminus * (lambda[n] - klam); /* s_{-}^r */
        }
      factor = (double)(sign * deltan) * rcoef / (klam * 2.0 * cosh(yMax * klam));
      genGrad.re += (splus - sminus) * factor * beta[n].re;
      genGrad.im += (splus - sminus) * factor * beta[n].im;
    }

  return (genGrad);
}

COMPLEX calcGGtopbottomB(COMPLEX *beta, double k, double *lambda, double yMax, int32_t grad_r, int32_t Ncoeff)
{
  COMPLEX genGrad;
  double klam, splus, sminus;
  double rcoef, factor;
  int32_t sign, power;

  int32_t n, i;

  rcoef = 1.0;
  for (i = grad_r; i >= 1; i--)
    rcoef *= 0.5 / (double)i; /* 1/(r!2^r) */

  genGrad.re = 0.0;
  genGrad.im = 0.0;
  for (n = 0; n < Ncoeff; n++)
    {
      if ((grad_r + n) % 2 == 0)
        {
          sign = -1;
          power = (grad_r + n) / 2;
          for (i = 1; i < power; i++)
            sign = -1 * sign; /* (-1)^{(ir+n)/2} */
        }
      else
        sign = 0;
      klam = sqrt(k * k + lambda[n] * lambda[n]);
      splus = lambda[n] + klam;
      sminus = lambda[n] - klam;
      for (i = 1; i < grad_r; i++)
        {
          splus = splus * (lambda[n] + klam);   /* s_{+}^r */
          sminus = sminus * (lambda[n] - klam); /* s_{-}^r */
        }
      factor = (double)sign * rcoef / (klam * 2.0 * sinh(yMax * klam));
      genGrad.re += (splus + sminus) * factor * beta[n].re;
      genGrad.im += (splus + sminus) * factor * beta[n].im;
    }

  return (genGrad);
}

COMPLEX calcGGrightA(COMPLEX *beta, double k, double *tau, double xMax, int32_t grad_r, int32_t Ncoeff)
{
  COMPLEX genGrad;
  double ktau, qplus, qminus;
  double rcoef, factor;
  int32_t theta0r, sign, power;

  int32_t n, i;

  rcoef = 1.0;
  for (i = grad_r; i >= 1; i--)
    rcoef *= 0.5 / (double)i; /* 1/(r!2^r) */
  if (grad_r % 2 == 0)
    theta0r = 0;
  else
    theta0r = 1;
  genGrad.re = 0.0;
  genGrad.im = 0.0;
  for (n = 1; n < Ncoeff; n += 2)
    {
      sign = 1;
      power = n / 2;
      for (i = 0; i < power; i++)
        sign = -1 * sign; /* (-1)^{floor[n/2]} */
      ktau = sqrt(k * k + tau[n] * tau[n]);
      qplus = -(tau[n] + ktau);
      qminus = tau[n] - ktau;
      for (i = 1; i < grad_r; i++)
        {
          qplus = -qplus * (tau[n] + ktau);  /* (-s_{+})^r */
          qminus = qminus * (tau[n] - ktau); /* s_{-}^r */
        }
      factor = rcoef * (-(double)(sign * theta0r) / (ktau * 2.0 * sinh(xMax * ktau)) + (double)(sign * (1 - theta0r)) / (ktau * 2.0 * cosh(xMax * ktau)));
      genGrad.re += (qminus - qplus) * factor * beta[n].re;
      genGrad.im += (qminus - qplus) * factor * beta[n].im;
    }

  return (genGrad);
}

COMPLEX calcGGrightB(COMPLEX *beta, double k, double *tau, double xMax, int32_t grad_r, int32_t Ncoeff)
{
  COMPLEX genGrad;
  double ktau, qplus, qminus;
  double rcoef, factor;
  int32_t theta0r, sign, power;

  int32_t n, i;

  rcoef = 1.0;
  for (i = grad_r; i >= 1; i--)
    rcoef *= 0.5 / (double)i; /* 1/(r!2^r) */
  if (grad_r % 2 == 0)
    theta0r = 0;
  else
    theta0r = 1;
  genGrad.re = 0.0;
  genGrad.im = 0.0;
  for (n = 0; n < Ncoeff; n += 2)
    {
      sign = 1;
      power = n / 2;
      for (i = 0; i < power; i++)
        sign = -1 * sign; /* (-1)^{floor[n/2]} */
      ktau = sqrt(k * k + tau[n] * tau[n]);
      qplus = -(tau[n] + ktau);
      qminus = tau[n] - ktau;
      for (i = 1; i < grad_r; i++)
        {
          qplus = -qplus * (tau[n] + ktau);  /* (-s_{+})^r */
          qminus = qminus * (tau[n] - ktau); /* s_{-}^r */
        }
      factor = rcoef * ((double)(sign * (1 - theta0r)) / (ktau * 2.0 * sinh(xMax * ktau)) - (double)(sign * theta0r) / (ktau * 2.0 * cosh(xMax * ktau)));
      genGrad.re += (qminus + qplus) * factor * beta[n].re;
      genGrad.im += (qminus + qplus) * factor * beta[n].im;
    }

  return (genGrad);
}

COMPLEX calcGGleftA(COMPLEX *beta, double k, double *tau, double xMax, int32_t grad_r, int32_t Ncoeff)
{
  COMPLEX genGrad;
  double ktau, qplus, qminus;
  double rcoef, factor;
  int32_t theta0r, sign, power;

  int32_t n, i;

  rcoef = 1.0;
  for (i = grad_r; i >= 1; i--)
    rcoef *= 0.5 / (double)i; /* 1/(r!2^r) */
  if (grad_r % 2 == 0)
    theta0r = 0;
  else
    theta0r = 1;
  genGrad.re = 0.0;
  genGrad.im = 0.0;
  for (n = 1; n < Ncoeff; n += 2)
    {
      sign = 1;
      power = n / 2;
      for (i = 0; i < power; i++)
        sign = -1 * sign; /* (-1)^{floor[n/2]} */
      ktau = sqrt(k * k + tau[n] * tau[n]);
      qplus = -(tau[n] + ktau);
      qminus = tau[n] - ktau;
      for (i = 1; i < grad_r; i++)
        {
          qplus = -qplus * (tau[n] + ktau);  /* (-s_{+})^r */
          qminus = qminus * (tau[n] - ktau); /* s_{-}^r */
        }
      factor = rcoef * (-(double)(sign * theta0r) / (ktau * 2.0 * sinh(xMax * ktau)) - (double)(sign * (1 - theta0r)) / (ktau * 2.0 * cosh(xMax * ktau)));
      genGrad.re += (qminus - qplus) * factor * beta[n].re;
      genGrad.im += (qminus - qplus) * factor * beta[n].im;
    }

  return (genGrad);
}

COMPLEX calcGGleftB(COMPLEX *beta, double k, double *tau, double xMax, int32_t grad_r, int32_t Ncoeff)
{
  COMPLEX genGrad;
  double ktau, qplus, qminus;
  double rcoef, factor;
  int32_t theta0r, sign, power;

  int32_t n, i;

  rcoef = 1.0;
  for (i = grad_r; i >= 1; i--)
    rcoef *= 0.5 / (double)i; /* 1/(r!2^r) */
  if (grad_r % 2 == 0)
    theta0r = 0;
  else
    theta0r = 1;
  genGrad.re = 0.0;
  genGrad.im = 0.0;
  for (n = 0; n < Ncoeff; n += 2)
    {
      sign = 1;
      power = n / 2;
      for (i = 0; i < power; i++)
        sign = -1 * sign; /* (-1)^{floor[n/2]} */
      ktau = sqrt(k * k + tau[n] * tau[n]);
      qplus = -(tau[n] + ktau);
      qminus = tau[n] - ktau;
      for (i = 1; i < grad_r; i++)
        {
          qplus = -qplus * (tau[n] + ktau);  /* (-s_{+})^r */
          qminus = qminus * (tau[n] - ktau); /* s_{-}^r */
        }
      factor = rcoef * ((double)(sign * (1 - theta0r)) / (ktau * 2.0 * sinh(xMax * ktau)) + (double)(sign * theta0r) / (ktau * 2.0 * cosh(xMax * ktau)));
      genGrad.re += (qminus + qplus) * factor * beta[n].re;
      genGrad.im += (qminus + qplus) * factor * beta[n].im;
    }

  return (genGrad);
}

COMPLEX calcGGtopbottomk0A(COMPLEX *beta, double *lambda, double yMax, int32_t grad_r, int32_t Ncoeff)
{
  COMPLEX genGrad;
  double splus;
  double rcoef, factor;
  int32_t theta0r, deltan, sign, power;

  int32_t n, i;

  rcoef = 1.0;
  power = grad_r / 2;
  for (i = 0; i < power; i++)
    rcoef = -1.0 * rcoef; /* (-1)^{floor[r/2]} */
  for (i = grad_r; i >= 1; i--)
    rcoef *= 0.5 / (double)i; /* 1/(r!2^r) */
  if (grad_r % 2 == 0)
    theta0r = 0;
  else
    theta0r = 1;
  genGrad.re = 0.0;
  genGrad.im = 0.0;
  for (n = 0; n < Ncoeff; n++)
    {
      if (n % 2 == 0)
        deltan = theta0r;
      else
        deltan = 1 - theta0r;
      sign = 1;
      power = n / 2;
      for (i = 0; i < power; i++)
        sign = -1 * sign; /* (-1)^{floor[n/2]} */
      splus = 1.0;
      for (i = 1; i < grad_r; i++)
        splus = splus * 2.0 * lambda[n]; /* (2*lambda)^(grad_r-1) */
      factor = (double)(sign * deltan) * rcoef * splus / (cosh(yMax * lambda[n]));
      genGrad.re += factor * beta[n].re;
      genGrad.im += factor * beta[n].im;
    }

  return (genGrad);
}

COMPLEX calcGGtopbottomk0B(COMPLEX *beta, double *lambda, double yMax, int32_t grad_r, int32_t Ncoeff)
{
  COMPLEX genGrad;
  double splus;
  double rcoef, factor;
  int32_t sign, power;

  int32_t n, i;

  rcoef = 1.0;
  for (i = grad_r; i >= 1; i--)
    rcoef *= 0.5 / (double)i; /* 1/(r!2^r) */
  //printf("1/rcoef = %e\n", 1.0/rcoef);
  if (grad_r == 2)
    {
      genGrad.re = -(rcoef / yMax) * beta[0].re;
      genGrad.im = -(rcoef / yMax) * beta[0].im;
    }
  else
    {
      genGrad.re = 0.0;
      genGrad.im = 0.0;
    }
  for (n = 1; n < Ncoeff; n++)
    {
      if ((grad_r + n) % 2 == 0)
        {
          sign = -1;
          power = (grad_r + n) / 2;
          for (i = 1; i < power; i++)
            sign = -1 * sign; /* (-1)^{(ir+n)/2} */
        }
      else
        sign = 0;
      splus = 1.0;
      for (i = 1; i < grad_r; i++)
        splus = splus * 2.0 * lambda[n]; /* (2*lambda)^(grad_r-1) */
      factor = (double)sign * rcoef * splus / (sinh(yMax * lambda[n]));
      genGrad.re += factor * beta[n].re;
      genGrad.im += factor * beta[n].im;
    }

  return (genGrad);
}

COMPLEX calcGGrightk0(COMPLEX *beta, double *tau, double xMax, int32_t grad_r, int32_t Ncoeff)
{
  COMPLEX genGrad;
  double qplus;
  double rcoef, factor;
  int32_t theta0r, sign, power;

  int32_t n, i;

  rcoef = 1.0;
  for (i = grad_r; i >= 1; i--)
    rcoef *= 0.5 / (double)i; /* 1/(r!2^r) */
  if (grad_r % 2 == 0)
    theta0r = 0;
  else
    theta0r = 1;
  switch (grad_r)
    {
    case 1:
      genGrad.re = rcoef * beta[0].re;
      genGrad.im = rcoef * beta[0].im;
      break;
    case 2:
      genGrad.re = (rcoef * 2.0 / xMax) * beta[0].re;
      genGrad.im = (rcoef * 2.0 / xMax) * beta[0].im;
      break;
    default:
      genGrad.re = 0.0;
      genGrad.im = 0.0;
    }
  for (n = 2; n < Ncoeff; n += 2)
    {
      sign = 1;
      power = n / 2;
      for (i = 0; i < power; i++)
        sign = -1 * sign; /* (-1)^{floor[n/2]} */
      qplus = -1.0;
      for (i = 1; i < grad_r; i++)
        qplus = -qplus * 2.0 * tau[n]; /* (-s_{+})^{r-1} */
      factor = rcoef * ((double)(sign * (1 - theta0r)) / sinh(xMax * tau[n]) - (double)(sign * theta0r) / cosh(xMax * tau[n]));

      genGrad.re += qplus * factor * beta[n].re;
      genGrad.im += qplus * factor * beta[n].im;
    }
  return (genGrad);
}

COMPLEX calcGGleftk0(COMPLEX *beta, double *tau, double xMax, int32_t grad_r, int32_t Ncoeff)
{
  COMPLEX genGrad;
  double qplus;
  double rcoef, factor;
  int32_t theta0r, sign, power;

  int32_t n, i;

  rcoef = 1.0;
  for (i = grad_r; i >= 1; i--)
    rcoef *= 0.5 / (double)i; /* 1/(r!2^r) */
  if (grad_r % 2 == 0)
    theta0r = 0;
  else
    theta0r = 1;
  switch (grad_r)
    {
    case 1:
      genGrad.re = -rcoef * beta[0].re;
      genGrad.im = -rcoef * beta[0].im;
      break;
    case 2:
      genGrad.re = (rcoef * 2.0 / xMax) * beta[0].re;
      genGrad.im = (rcoef * 2.0 / xMax) * beta[0].im;
      break;
    default:
      genGrad.re = 0.0;
      genGrad.im = 0.0;
    }
  for (n = 2; n < Ncoeff; n += 2)
    {
      sign = 1;
      power = n / 2;
      for (i = 0; i < power; i++)
        sign = -1 * sign; /* (-1)^{floor[n/2]} */
      qplus = -1.0;
      for (i = 1; i < grad_r; i++)
        qplus = -qplus * 2.0 * tau[n]; /* (-s_{+})^{r-1} */
      factor = rcoef * ((double)(sign * (1 - theta0r)) / sinh(xMax * tau[n]) + (double)(sign * theta0r) / cosh(xMax * tau[n]));
      genGrad.re += qplus * factor * beta[n].re;
      genGrad.im += qplus * factor * beta[n].im;
    }

  return (genGrad);
}

COMPLEX calcGGallSidesBz(COMPLEX *beta, double k, double *lambda, double yMax, int32_t Ncoeff)
{
  COMPLEX genGrad;
  double klambda, factor;
  int32_t sign, power;

  int32_t n, i;

  genGrad.re = 0.0;
  genGrad.im = 0.0;
  for (n = 1; n < Ncoeff; n += 2)
    {
      sign = 1;
      power = n / 2;
      for (i = 0; i < power; i++)
        sign = -1 * sign; /* (-1)^{floor[n/2]} */
      klambda = sqrt(k * k + lambda[n] * lambda[n]);
      factor = (double)sign / (2.0 * cosh(yMax * klambda));
      genGrad.re += factor * beta[n].re;
      genGrad.im += factor * beta[n].im;
    }

  return (genGrad);
}

#ifndef OLDFFT
void FFT(COMPLEX *field, int32_t isign, int32_t npts)
{
  double *real_imag;
  int32_t i;
#ifdef DEBUG
  static FILE *fpfft = NULL;
  if (!fpfft) {
    fpfft = fopen("computeRBGGE.fft", "w");
    fprintf(fpfft, "SDDS1\n");
    fprintf(fpfft, "&column name=i type=long &end\n");
    fprintf(fpfft, "&column name=Real type=double &end\n");
    fprintf(fpfft, "&column name=Imag type=double &end\n");
    fprintf(fpfft, "&data mode=ascii &end\n");
  }
#endif

  real_imag = tmalloc(sizeof(double) * (2 * npts + 2));
  for (i = 0; i < npts; i++)
    {
      real_imag[2 * i] = field[i].re;
      real_imag[2 * i + 1] = field[i].im;
    }
  if (isign == -1)
    {
      complexFFT(real_imag, npts, 0);
      for (i = 0; i < npts; i++)
        {
          field[i].re = npts*real_imag[2 * i];
          field[i].im = npts*real_imag[2 * i + 1];
        }
    }
  else
    {
      complexFFT(real_imag, npts, INVERSE_FFT);
      for (i = 0; i < npts; i++)
        {
          field[i].re = real_imag[2 * i];
          field[i].im = real_imag[2 * i + 1];
        }
    }
#ifdef DEBUG
  fprintf(fpfft, "%" PRId32 "\n", npts);
  for (i=0; i<npts; i++)
    fprintf(fpfft, "%" PRId32 " %le %le\n", i, field[i].re, field[i].im);
#endif
  free(real_imag);
}
#else
void FFT(COMPLEX *field, int32_t isign, int32_t npts)
{
  unsigned long mmax, m, hn, j, istep, i;
  double wtemp, wr, wpr, wpi, wi, theta, flt_isign;
  /* must be double to preserve accuracy */
  COMPLEX tempz;

  COMPLEX *fft;
  double tempr, tempi;

#ifdef DEBUG
  static FILE *fpfft = NULL;
  if (!fpfft) {
    fpfft = fopen("computeRBGGE.fft", "w");
    fprintf(fpfft, "SDDS1\n");
    fprintf(fpfft, "&column name=i type=long &end\n");
    fprintf(fpfft, "&column name=Real type=double &end\n");
    fprintf(fpfft, "&column name=Imag type=double &end\n");
    fprintf(fpfft, "&data mode=ascii &end\n");
  }
#endif

  fft = calloc(npts, sizeof(COMPLEX));

  hn = npts / 2;

  /* copy over */
  for (i = 0; i < npts; i++)
    fft[i] = field[i];

  /* Fourier Transform isign = -1, inverse isign = 1 */

  flt_isign = (double)(isign);
  /* first, sort into bit-reversed order */
  for (j = 0, i = 0; i < npts; i++) /* increment in regular order */
    {
      if (j > i) /* swap i and j = bit-reversal(i) once, if distinct */
        {
          tempz = fft[i];
          fft[i] = fft[j];
          fft[j] = tempz;
        }
      for (m = hn; ((j >= m) && (m >= 1)); m >>= 1)
        /* find bit-reversal of next i */
        {
          j -= m;
        }
      j += m;
    }

  /* next, apply Danielson-Lanczos algorithm */
  for (mmax = 1; (npts > mmax); mmax = istep)
    /* loop through log_base2(N) times */
    {
      istep = mmax << 1; /* = 2*mmax */
      /* initialize trig functions */
      theta = (flt_isign) * (PI / ((double)(mmax)));
      wtemp = sin(0.5 * theta);
      wpr = -2.0 * wtemp * wtemp;
      wpi = sin(theta);
      wr = 1.0;
      wi = 0.0;
      for (m = 0; m < mmax; m += 1)
        {
          for (i = m; i < npts; i += istep)
            {
              j = i + mmax;
              tempr = wr * fft[j].re - wi * fft[j].im;
              tempi = wr * fft[j].im + wi * fft[j].re;
              fft[j].re = fft[i].re - tempr;
              fft[j].im = fft[i].im - tempi;
              fft[i].re += tempr;
              fft[i].im += tempi;
            }
          /* update trig functions via recurrence relations */
          wr = wr + (wtemp = wr) * wpr - wi * wpi;
          wi = wi + wi * wpr + wtemp * wpi;
        }
    }

  /* copy */
  for (i = 0; i < npts; i++)
    {
      field[i] = fft[i];

    }
  free(fft);

#ifdef DEBUG
  fprintf(fpfft, "%ld\n", npts);
  for (i=0; i<npts; i++)
    fprintf(fpfft, "%ld %le %le\n", i, field[i].re, field[i].im);
#endif
}
#endif

unsigned long IntCeilingPowerOf2(unsigned long i)
{
  /* returns smallest non-negative x = 2^n  such that i <= x */

  unsigned long x;

  for (x = 1; (i > x); x <<= 1)
    {
    }
  return x;
}
