
#include "SDDS.h"
#include "scan.h"
#include "constants.h"
#include "mdb.h"
#include <gsl/gsl_sf_airy.h>
#include <gsl/gsl_complex_math.h>
#include <math.h>
  
#define CLO_HEIGHT 0
#define CLO_FREQUENCY_LIMIT 1
#define CLO_NUMBER 2
#define CLO_RADIUS 3
#define CLO_PIPE 4
#define N_OPTIONS 5

char *option[N_OPTIONS]= {
  "height", "frequencyLimit", "n", "radius", "pipe"
  };

char *USAGE = "csrImpedance <outputFile> -height=<valueInMeters> -radius=<valueInMeters> \
         -frequencyLimit=maximum=<valueInHz>,minimum=<ValueInHz> -n=<interger> -pipe[=out]";

typedef gsl_complex fcomplex; 
const double Z0 = 376.730313461770606;   /* mu_o * c_mks  */
const double C5=1.004524;

fcomplex sum_F0(double height, double k, double radius);
void csr_impedance(double height, double radius, double fmin, double fmax, long N, SDDS_DATASET *SDDS_out);
fcomplex Z_asymtotic(double height, double radius, double f);
void SetupOutput(char *filename, SDDS_DATASET *SDDS_out);


int main (int argc, char **argv)
{
  SDDS_DATASET SDDS_out;
  char *outputFile=NULL;
  SCANNED_ARG *s_arg;
  double height=-1, radius=-1, fmin=0, fmax=-1;
  long N=0, i_arg;
  unsigned long pipeFlags=0, dummyFlags=0;
    
  argc = scanargs(&s_arg, argc, argv);
  if (argc<2) 
    bomb(NULL, USAGE);
  for (i_arg=1; i_arg<argc; i_arg++) {
    if (s_arg[i_arg].arg_type==OPTION) {
      delete_chars(s_arg[i_arg].list[0], "_");
      switch (match_string(s_arg[i_arg].list[0], option, N_OPTIONS, 0)) {
      case CLO_HEIGHT:
        if (s_arg[i_arg].n_items!=2 || !get_double(&height, s_arg[i_arg].list[1]))
          SDDS_Bomb("Invalid -weight syntax provided");
        break;
      case CLO_RADIUS:
        if (s_arg[i_arg].n_items!=2 || !get_double(&radius, s_arg[i_arg].list[1]))
          SDDS_Bomb("Invalid -radius syntax provided");
        break;
      case CLO_NUMBER:
        if (s_arg[i_arg].n_items!=2 || !get_long(&N, s_arg[i_arg].list[1]))
          SDDS_Bomb("Invalid -radius syntax provided");
        break;
      case CLO_FREQUENCY_LIMIT:
        if (s_arg[i_arg].n_items<2)
          SDDS_Bomb("Invalid -frequencyLimit syntax provided"); 
        s_arg[i_arg].n_items -=1;
        if (!scanItemList(&dummyFlags, s_arg[i_arg].list+1, &s_arg[i_arg].n_items, 0,
                          "maximumFrequency", SDDS_DOUBLE, &fmax, 1, 0,
                          "minimumFrequency", SDDS_DOUBLE, &fmin, 1, 0, NULL))
          SDDS_Bomb("invalid -errorTemplate syntax");
        s_arg[i_arg].n_items +=1; 
        break;
      case CLO_PIPE:
        if (!processPipeOption(s_arg[i_arg].list+1, s_arg[i_arg].n_items-1, &pipeFlags))
          SDDS_Bomb("invalid -pipe syntax");
        if (pipeFlags!=USE_STDOUT)
          SDDS_Bomb("only -pipe=out syntax is valid!");
        break;
      }
    } else if (!outputFile)
      outputFile = s_arg[i_arg].list[0];
    else
      SDDS_Bomb("Too many files provided.");
  }
  if (!outputFile && !pipeFlags)
    SDDS_Bomb("output file is not provided!");
  if (!N)
    SDDS_Bomb("The number N is not provided.");
  if (fmin<0 || fmax<0)
    SDDS_Bomb("Invalid minimum or maximum frequency, should not be less than 0.");
  if (height<0)
    SDDS_Bomb("The height is not provided.");
  if (radius<0)
    SDDS_Bomb("The radius is not provided.");
 
  SetupOutput(outputFile, &SDDS_out);
  csr_impedance (height, radius, fmin, fmax, N, &SDDS_out); 
  if (!SDDS_Terminate(&SDDS_out))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  free_scanargs(&s_arg, argc); 
  return 0;
}

void csr_impedance (double height, double radius, double fmin, double fmax, long N, SDDS_DATASET *SDDS_out) 
{
  long Np, i; 
  double df, f, k;
  fcomplex *Z, *Z_LLimit;
  Np = pow(2,N)+1;
  df=(fmax-fmin)/(Np-1);
  
  Z=(fcomplex *)calloc(sizeof(*Z), Np);
  Z_LLimit=(fcomplex *)calloc(sizeof(*Z), Np);
  if (!SDDS_StartPage(SDDS_out, Np) ||
      !SDDS_SetParameters(SDDS_out, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,
                          "Height", height, "Radius", radius, "MaximumFrequency", fmax, 
                          "MinimumFrequency", fmin, "N", N, NULL))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors); 
  for (i=0; i<Np; i++) 
  {
    f = fmin+df*i;
    k = 2*PI*f/c_mks;
    
    if (fabs(k)>1e-6) {
      Z[i] = gsl_complex_mul_real(sum_F0(height,k,radius),2*PI/height*pow(2.0/k/radius,1.0/3)*Z0); 
      Z_LLimit[i] = Z_asymtotic (height, radius, k);
    }
    else
      Z[i] = gsl_complex_rect(0.0, 0.0);

    /* set it to zero if the number is too small to make sddsplot draw pictures properly */ 
    if (fabs(GSL_REAL(Z[i]))<1e-32)
      Z[i] = gsl_complex_rect (0.0, GSL_IMAG(Z[i]));
    if (fabs(GSL_IMAG(Z[i]))<1e-32)
      Z[i] = gsl_complex_rect (GSL_REAL(Z[i]), 0.0);
    if (fabs(GSL_REAL(Z_LLimit[i]))<1e-32)
      Z_LLimit[i] = gsl_complex_rect (0.0, GSL_IMAG(Z_LLimit[i]));
    if (fabs(GSL_IMAG(Z_LLimit[i]))<1e-32)
      Z_LLimit[i] = gsl_complex_rect (GSL_REAL(Z_LLimit[i]), 0.0);
    if (!SDDS_SetRowValues(SDDS_out,  SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, i, 
                           "f", f, "k", k, "ZReal", GSL_REAL(Z[i])*2*PI*radius, "ZImag", GSL_IMAG(Z[i])*2*PI*radius*-1, NULL))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  if (!SDDS_WritePage(SDDS_out))
     SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  free(Z_LLimit); 
  free(Z);
}

fcomplex sum_F0 (double h, double k, double r)
{
  int p, pMax=50;
  double beta, beta2, tmp=PI/h*pow(r/2/k/k, 1.0/3); 
  fcomplex sum=gsl_complex_rect(0.0,0.0), tmp1, tmp2, tmp3;
  double ai,bi,aip,bip;

  for (p=0; p<pMax; p++) {
    beta = (2*p+1)*tmp;
    beta2 = beta*beta;
    if (beta2<100.0) {
      ai = gsl_sf_airy_Ai (beta2, GSL_PREC_DOUBLE);
      bi = gsl_sf_airy_Bi (beta2, GSL_PREC_DOUBLE);
      aip = gsl_sf_airy_Ai_deriv (beta2, GSL_PREC_DOUBLE);
      bip = gsl_sf_airy_Bi_deriv (beta2, GSL_PREC_DOUBLE);
      /*   printf("%e\t%16.6e %16.6e %16.6e %16.6e\n",beta2,ai,bi,aip,bip); */
      tmp1=gsl_complex_rect(aip,-bip);
      tmp2=gsl_complex_rect(ai, -bi);
      tmp3=  gsl_complex_add(gsl_complex_mul_real(tmp1,aip),gsl_complex_mul_real(tmp2,beta2*ai));
      sum = gsl_complex_add (sum,tmp3); 
      /* printf ("%e\t%e\t%e\n",
         beta2, GSL_REAL(tmp3), GSL_IMAG(tmp3));    */
      /*   printf ("aip*bip=%lf, beta2*ai*bi=%lf\n", aip*bip, beta2*ai*bi); */
    }
    else
      break;
  } 
  return sum; 
}

fcomplex Z_asymtotic  (double h, double r, double k)
{
  double k2h3=pow(k,2.0)*pow(h,3.0), PI3=pow(PI,3.0); 
  return gsl_complex_mul_real( gsl_complex_rect(exp(-2*PI3/r/3/k2h3), -1.5*C5*pow(k2h3/PI3/r, 2.0)),PI*Z0/k/h/h);
}

void SetupOutput (char *filename, SDDS_DATASET *SDDS_out) {
  if (!SDDS_InitializeOutput(SDDS_out, SDDS_BINARY, 1, NULL, NULL, filename))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  if (SDDS_DefineParameter(SDDS_out, "Height", "height", "m", "The height of the chamber", NULL, SDDS_DOUBLE, NULL)<0 ||
      SDDS_DefineParameter(SDDS_out, "Radius", "radius", "m", "Radius", NULL, SDDS_DOUBLE, NULL)<0 ||
      SDDS_DefineParameter(SDDS_out, "MaximumFrequency", NULL, "Hz", NULL, NULL, SDDS_DOUBLE, NULL)<0 ||
      SDDS_DefineParameter(SDDS_out, "MinimumFrequency", NULL, "Hz", NULL, NULL, SDDS_DOUBLE, NULL)<0 ||
      SDDS_DefineParameter(SDDS_out, "N", "N", NULL, "The power of two", NULL, SDDS_LONG, NULL)<0 ||
      SDDS_DefineColumn(SDDS_out, "f", "f", "Hz", "frequency", NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineColumn(SDDS_out, "k", "k", NULL, "wave number", NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineColumn(SDDS_out, "ZReal", "Re[Z]", "Ohms", "the real part of the impedance", NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineColumn(SDDS_out, "ZImag", "Im[Z]", "Ohms", "the imaginary part of the impedance", NULL, SDDS_DOUBLE, 0)<0 ) 
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  if (!SDDS_WriteLayout(SDDS_out))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
}
