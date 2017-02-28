/*************************************************************************\
* Copyright (c) 2015 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2015 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

#include "mdb.h"
#include "track.h"
void bombElegantVA(char *template, ...);


int main(int argc, char **argv) 
{
  char *input, *output, *filterFile;
  SDDS_DATASET SDDSin, SDDSout;
  IIRFILTER filterBank[4];
  long i, rows, nFilters;
  double *x, *y;
  
  if (argc<4) 
    bomb("usage: iirFilter <input> <output> <filterFile>", NULL);
  input = argv[1];
  output = argv[2];
  filterFile = argv[3];

  for (i=0; i<4; i++)
    memset(&filterBank[i], 0, sizeof(filterBank[i]));

  if (!SDDS_InitializeInput(&SDDSin, input) ||
      SDDS_ReadPage(&SDDSin)!=1 ||
      (rows=SDDS_RowCount(&SDDSin))<=0 ||
      !(x=SDDS_GetColumnInDoubles(&SDDSin, "x")) ||
      !SDDS_Terminate(&SDDSin))
    SDDS_Bomb("Problem with input file");

  if (!SDDS_InitializeOutput(&SDDSout, SDDS_BINARY, 1, NULL, NULL, output) ||
      SDDS_DefineColumn(&SDDSout, "y", NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      !SDDS_WriteLayout(&SDDSout) ||
      !SDDS_StartPage(&SDDSout, rows))
    SDDS_Bomb("Problem with output file");

  nFilters = readIIRFilter(filterBank, 4, filterFile);
  y = tmalloc(sizeof(*y)*rows);
  for (i=0; i<rows; i++) 
    y[i] = applyIIRFilter(filterBank, nFilters, x[i]);
  if (!SDDS_SetColumn(&SDDSout, SDDS_SET_BY_NAME, y, rows, "y") || !SDDS_WritePage(&SDDSout) || !SDDS_Terminate(&SDDSout))
    SDDS_Bomb("Problem with output file (2)");

  free(x);
  free(y);
  
  freeIIRFilterMemory(filterBank, nFilters);

  return 0;
}

#if 0
void bombElegantVA(char *template, ...) 
{
  char *p;
  char c, *s;
  int i;
  long j;
  va_list argp;
  char buffer[256];
  float f;
  double d;
  
  va_start(argp, template);
  p = template;
  while (*p) {
    if (*p=='%') {
      switch (*++p) {
      case 'l':
        switch (*++p) {
        case 'd':
          j = va_arg(argp, long int);
          printf("%ld", j);
          break;
        case 'e':
          d =  va_arg(argp, double);
          printf("%21.15le", f);
          break;
        case 'f':
          d =  va_arg(argp, double);
          printf("%lf", f);
          break;
        case 'g':
          d =  va_arg(argp, double);
          printf("%21.15lg", f);
          break;
        default:
          printf("%%l%c", *p);
          break;
        }
        break;
      case 'c':
        c = va_arg(argp, int);
        putchar(c);
        break;
      case 'd':
        i = va_arg(argp, int);
        printf("%d", i);
        break;
      case 's':
        s = va_arg(argp, char *);
        fputs(s, stdout);
        break;
      case 'e':
        d = va_arg(argp, double);
        printf("%21.15e", d);
        break;
      case 'f':
        d = va_arg(argp, double);
        printf("%f", d);
        break;
      case 'g':
        d = va_arg(argp, double);
        printf("%21.15g", d);
        break;
      default:
        printf("%%%c", *p);
        break;
      }
    }
    else {
      putchar(*p);
    }
    p++;
  }
  exit(1);
}

  
void bombElegantTemplate(char *template, char *a1, char *a2, char *usage)
{
  fputs("error: ", stdout);
  if (a1 && a2)
    printf(template, a1, a2);
  else if (a1) 
    printf(template, a1);
  else
    printf(template);
  fputc('\n', stdout);
  if (usage)
    printf("usage: %s\n", usage);
  exit(1);
}


void bombElegant(char *error, char *usage)
{
    bombElegantTemplate(error, NULL, NULL, usage);
}

int readIIRFilter(IIRFILTER *filterBank, long maxFilters, char *inputFile)
{
  SDDS_DATASET SDDSin;
  int i, code;
  
  if (!inputFile)
    return 0;
  if (maxFilters<=0 || filterBank==NULL)
    bombElegantVA("Error: readIIRFilter called with invalid arguments for file %s\n", inputFile);
  if (!SDDS_InitializeInputFromSearchPath(&SDDSin, inputFile)) 
    bombElegantVA("Error: Problem with IIR filter file %s\n", inputFile);
  if (SDDS_CheckColumn(&SDDSin, (char*)"denominator", NULL, SDDS_ANY_FLOATING_TYPE, stdout)!=SDDS_CHECK_OK ||
      SDDS_CheckColumn(&SDDSin, (char*)"numerator", NULL, SDDS_ANY_FLOATING_TYPE, stdout)!=SDDS_CHECK_OK)
    bombElegantVA("Error: Problem with existence/type of columns in IIR filter file %s\n", inputFile);

  for (i=0; i<maxFilters; i++) {
    filterBank[i].nTerms = filterBank[i].iBuffer = 0;
    if (filterBank[i].an) 
      free(filterBank[i].an);
    if (filterBank[i].bn) 
      free(filterBank[i].bn);
    if (filterBank[i].xn)
      free(filterBank[i].xn);
    if (filterBank[i].yn)
      free(filterBank[i].yn);
    filterBank[i].an = filterBank[i].bn = filterBank[i].xn = filterBank[i].yn = NULL;
  }
    
  i = 0;
  while ((code=SDDS_ReadPage(&SDDSin))>0 && i<maxFilters) {
    filterBank[i].nTerms = SDDS_RowCount(&SDDSin);
    
    if (!(filterBank[i].an = SDDS_GetColumnInDoubles(&SDDSin, (char*)"denominator")) ||
        !(filterBank[i].bn = SDDS_GetColumnInDoubles(&SDDSin, (char*)"numerator")))
      bombElegantVA("Error: Problem retrieving column data for IIR filter file %s, page %ld\n", inputFile, code);

    if (!(filterBank[i].xn=calloc(filterBank[i].nTerms, sizeof(*(filterBank[i].xn)))) ||
        !(filterBank[i].yn=calloc(filterBank[i].nTerms, sizeof(*(filterBank[i].yn)))))
      bombElegantVA("Error: Memory allocation problem for IIR filter file %s, page %ld, %ld terms\n", inputFile, code, filterBank[i].nTerms);
    
    i++;
  }
  if (code==0)
    bombElegantVA("Error: IIR filter file %s has invalid page(s)\n", inputFile);
  if (i>=maxFilters)
    bombElegantVA("Error: IIR filter file %s has too many pages (%ld maximum allowed)\n", inputFile, maxFilters);

  if (!SDDS_Terminate(&SDDSin))
    bombElegantVA("Error: Problem terminating IIR filter file %s\n", inputFile);

  return i;
}


double applyIIRFilter(IIRFILTER *filterBank, long nFilters, double x)
{
  double y, output;
  long i, j, ib, d;
  
  output = 0;

  for (i=0; i<nFilters; i++) {
    if (filterBank[i].nTerms<=0)
      continue;
    if ((--(filterBank[i].iBuffer))<0)
      filterBank[i].iBuffer = filterBank[i].nTerms-1;

    ib = filterBank[i].iBuffer;
    filterBank[i].xn[ib] = x;
    filterBank[i].yn[ib] = 0; /* for convenience in doing sums */
    
    y = 0;
    d = 0;
    for (j=ib; j<filterBank[i].nTerms; j++, d++) {
      y += filterBank[i].xn[j]*filterBank[i].bn[d];
      y -= filterBank[i].yn[j]*filterBank[i].an[d];
    }
    for (j=0; j<ib; j++, d++) {
      if (d>=filterBank[i].nTerms)
        bombElegant("Problem (1) with circular buffer handling for IIR filter", NULL);
      y += filterBank[i].xn[j]*filterBank[i].bn[d];
      y -= filterBank[i].yn[j]*filterBank[i].an[d];
    }
    if (d!=filterBank[i].nTerms) {
      bombElegant("Problem (2) with circular buffer handling for IIR filter", NULL);
    }
    filterBank[i].yn[ib] = y;
    output += y;
  }

  return output;
}

void freeIIRFilterMemory(IIRFILTER *filterBank, long nFilters)
{
  long i;
  for (i=0; i<nFilters; i++) {
    if (filterBank[i].an) free(filterBank[i].an);
    if (filterBank[i].bn) free(filterBank[i].bn);
    if (filterBank[i].xn) free(filterBank[i].xn);
    if (filterBank[i].yn) free(filterBank[i].yn);
    filterBank[i].an = filterBank[i].bn = filterBank[i].xn = filterBank[i].yn = NULL;
    filterBank[i].iBuffer = filterBank[i].nTerms = 0;
  }
}

#endif
