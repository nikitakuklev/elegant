/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

#include <stdio.h>
#include "mdb.h"
#include "SDDS.h"
#include "constants.h"
#include "chbook.h"

book1 *chbook1(char *vName, char *units, double xmin, double xmax, long xbins)
{
  book1 *book;
  book = (book1 *)malloc(sizeof(*book));

  book->vname = vName;
  book->units = units;

  book->xmin = xmin;
  book->xmax = xmax;
  book->dx = (xmax-xmin)/(double)xbins;
  book->xbins = xbins;  
  book->length = xbins;
  book->count = 0;
  book->total = 0;

  book->value = calloc(sizeof(*book->value), book->length);

  return book;
}

void chfill1(book1 *bName, double x, double Frequency)
{
  long index;

  index = (long)((x - bName->xmin)/bName->dx);
  if (index < 0) index =0;
  if (index > bName->xbins-1) index = bName->xbins-1;
  
  bName->value[index] += Frequency;
  bName->total += Frequency;
  bName->count++;

  return;
}

void free_hbook1 (book1 *x)
{
  free(x->value);
  free(x);
  return;
}

book2 *chbook2(char *xName, char *yName, char *xunits, char *yunits,
               double xmin, double xmax, double ymin, double ymax, long xbins, long ybins)
{
  book2 *book;

  book = (book2 *)malloc(sizeof(*book));

  book->xname = xName;
  book->yname = yName;
  book->xunits = xunits;
  book->yunits = yunits;

  book->xmin = xmin;
  book->xmax = xmax;
  book->dx = (xmax-xmin)/(double)xbins;
  book->xbins = xbins;  
  book->ymin = ymin;
  book->ymax = ymax;
  book->dy = (ymax-ymin)/(double)ybins;
  book->ybins = ybins;  

  book->length = xbins*ybins;
  book->count = 0;
  book->total = 0;

  book->value = calloc(sizeof(*book->value), book->length);

  return book;
}

void chfill2(book2 *bName, double x, double y, double Frequency)
{
  long index[2], i;

  index[0] = (long)((x - bName->xmin)/bName->dx);
  if (index[0] < 0) index[0] =0;
  if (index[0] > (bName->xbins-1)) index[0] = bName->xbins-1;
  
  index[1] = (long)((y - bName->ymin)/bName->dy);
  if (index[1] < 0) index[1] =0;
  if (index[1] > (bName->ybins-1)) index[1] = bName->ybins-1;
  
  i = index[0]*bName->xbins + index[1];
  bName->value[i] += Frequency;
  bName->total += Frequency;
  bName->count++;

  return;
}

void free_hbook2 (book2 *x)
{
  free(x->value);
  free(x);
  return;
}

ntuple *chbookn(char **vName, char **units, long ND, double *xmin, double *xmax, long *xbins, long offset)
{
  ntuple *book;
  long i;

  book = (ntuple *)malloc(sizeof(*book));

  book->nD = ND;

  book->vname = calloc(sizeof(*book->vname), ND);
  book->units = calloc(sizeof(*book->units), ND);
  book->xmin = calloc(sizeof(*book->xmin), ND);
  book->xmax = calloc(sizeof(*book->xmax), ND);
  book->dx = calloc(sizeof(*book->dx), ND);
  book->xbins = calloc(sizeof(*book->xbins), ND);
  
  book->count = 0;
  book->total = 0;
  book->length = 1;
  for (i=0; i<ND; i++) {
    book->vname[i] = vName[i+offset];
    book->units[i] = units[i+offset];
    book->length *= xbins[i+offset];
    book->xbins[i] = xbins[i+offset];
    book->xmin[i] = xmin[i+offset];
    book->xmax[i] = xmax[i+offset];
    book->dx[i] = (xmax[i+offset]-xmin[i+offset])/(double)xbins[i+offset];
  }

  book->value = calloc(sizeof(*book->value), book->length);

  return book;
}

void chfilln(ntuple *bName, double *x, double Frequency, long offset)
{
  long index;
  long i, temp;

  index = 0;
  for (i=0; i<bName->nD; i++) {
    temp = (long) ((x[i+offset] - bName->xmin[i])/bName->dx[i]);
    if (temp < 0) temp = 0;
    if (temp > bName->xbins[i]-1) temp = bName->xbins[i]-1;
    if(i>0) index = index*bName->xbins[i-1]+temp;
    if(i==0) index = temp;
  }
  
  bName->value[index] += Frequency;
  bName->total += Frequency;
  bName->count++;
  return;
}

void free_hbookn (ntuple *x)
{
  free(x->vname);
  free(x->units);
  free(x->xmin);
  free(x->xmax);
  free(x->dx);
  free(x->xbins);
  free(x->value);
  free(x);
  return;
}

book1m *chbook1m(char **vName, char **units, double *xmin, double *xmax, long *xbins, long column_number)
{
  book1m *book;
  long i, j;
  char **ptr0, *buffer;

  book = (book1m *)malloc(sizeof(*book));

  book->nD = book->bins = 0;
  for (i=0; i<column_number; i++) {
    if (!xbins[i]) continue;
    book->nD++;
    if (book->bins < xbins[i])
      book->bins = xbins[i];
  }

  book->vname = calloc(sizeof(*book->vname), book->nD);
  book->units = calloc(sizeof(*book->units), book->nD);
  book->xmin = calloc(sizeof(*book->xmin), book->nD);
  book->xmax = calloc(sizeof(*book->xmax), book->nD);
  book->dx = calloc(sizeof(*book->dx), book->nD);
  
  book->count = 0;
  book->total = 0;
  book->length = book->bins;
  for (i=j=0; i<column_number; i++) {
    if (!xbins[i]) continue;
    book->vname[j] = vName[i];
    book->units[j] = units[i];
    book->xmin[j] = xmin[i];
    book->xmax[j] = xmax[i];
    book->dx[j] = (xmax[i]-xmin[i])/(double)book->bins;
    j++;
  }

  ptr0 = (char**)tmalloc((unsigned)(sizeof(*ptr0)*book->length));
  buffer =  (char*)tmalloc((unsigned)(sizeof(*buffer)*sizeof(double)*book->length*book->nD));
  for (i=0; i<book->length; i++)
    ptr0[i] = buffer+i*sizeof(double)*book->nD;
  book->value = (double**)ptr0;
  for (i=0; i<book->length; i++)
    for (j=0; j<book->nD; j++)
      book->value[i][j] = 0;

  return book;
}

void chfill1m(book1m *bName, double *x, double Frequency, long *xbins, long column_number){
  long index, i, j;

  for (i=j=0; i<column_number; i++) {
    if (!xbins[i]) continue;    
    index = (long)((x[i] - bName->xmin[j])/bName->dx[j]);
    if (index < 0) index =0;
    if (index > bName->bins-1) index = bName->bins-1;
    bName->value[index][j] += Frequency;
    j++;
  } 
  bName->total += Frequency;
  bName->count++;
  return; 
}

void free_hbook1m(book1m *x){
  free(x->vname);
  free(x->units);
  free(x->xmin);
  free(x->xmax);
  free(x->dx);
  free(*(x->value));
  free(x->value);
  free(x);
  return;
}

void chprint1(book1 *bName, char *filename, char *description, EXTERNAL_PARA *parameter_definition,
              void **sdds_value, long n_parameters, long normalize, long verbosity, long append)
{
  SDDS_DATASET outPage;
  double value;
  long i, last_index, index;
  char buffer[100];

  last_index = -1;
  if (!append) {
    /* Open file for writting */
    if (verbosity)
      fprintf( stdout, "Opening \"%s\" for writing...\n", filename);

    if (!SDDS_InitializeOutput(&outPage, SDDS_BINARY, 1, 
                               description, description, filename))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    for (i=0; i<n_parameters; i++) {
      if (!SDDS_ProcessParameterString(&outPage, parameter_definition[i].text, 0) ||
          (index=SDDS_GetParameterIndex(&outPage, parameter_definition[i].name)<0)) {
        fprintf(stdout, "Unable to define SDDS parameter for chprint1--name was:\n%s\n",
                parameter_definition[i].name);
        fflush(stdout);
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
        exit(1);
      }
      if (index!=(last_index+1))
        fprintf(stdout, "\7\7\7WARNING: parameter indices for SDDS file %s are not sequential--this will probably cause unexpected results\n", filename);
      fflush(stdout);
      last_index = index;
    }

    if (0>SDDS_DefineParameter(&outPage, "VariableName", NULL, bName->units, 
                               NULL, NULL, SDDS_STRING, NULL) ||
        0>SDDS_DefineParameter(&outPage, "VariableInterval", NULL, NULL, 
                               NULL, NULL, SDDS_DOUBLE, NULL) ||
        0>SDDS_DefineParameter(&outPage, "VariableMinimum", NULL, NULL, 
                               NULL, NULL, SDDS_DOUBLE, NULL) ||
        0>SDDS_DefineParameter(&outPage, "VariableMaximum", NULL, NULL, 
                               NULL, NULL, SDDS_DOUBLE, NULL) ||
        0>SDDS_DefineParameter(&outPage, "VariableDimension", NULL, NULL, 
                               NULL, NULL, SDDS_LONG, NULL) ||
        0>SDDS_DefineParameter(&outPage, "Total_count", NULL, NULL, 
                               NULL, NULL, SDDS_LONG, NULL))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);

    sprintf(buffer, "%sFrequency", bName->vname);
    if (0>SDDS_DefineColumn(&outPage, bName->vname, NULL, bName->units,
                            NULL, NULL, SDDS_DOUBLE, 0) ||
        0>SDDS_DefineColumn(&outPage, buffer, NULL, NULL, 
                            NULL, NULL, SDDS_DOUBLE, 0))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);

    if (!SDDS_WriteLayout(&outPage) )
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  } else {
    if (!SDDS_InitializeAppend(&outPage, filename))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }

  /* Write to output file */
  if (0>SDDS_StartPage(&outPage, bName->length) ||
      !SDDS_SetParameters(&outPage, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,
                          "VariableName", bName->vname,
                          "VariableInterval", bName->dx,
                          "VariableMinimum", bName->xmin,
                          "VariableMaximum", bName->xmax,
                          "VariableDimension", bName->length,
                          "Total_count", bName->count, NULL))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);

  for (i=0; i<n_parameters; i++) {
    if (!SDDS_SetParameters(&outPage, SDDS_SET_BY_INDEX|SDDS_PASS_BY_REFERENCE, 
                            i, sdds_value[i], -1))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }

  for (i=0; i<bName->length; i++) {
    value = ((double)i+0.5)*bName->dx+bName->xmin;
    if (normalize) {
      if (!SDDS_SetRowValues(&outPage, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, i,
                             0, value, 1, bName->value[i]/bName->total, -1))
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    } else {
      if (!SDDS_SetRowValues(&outPage, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, i,
                             0, value, 1, bName->value[i], -1))
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
  }
  if (!SDDS_WritePage(&outPage))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  if (!SDDS_Terminate(&outPage))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  return;
}

void chprint2(book2 *bName, char *filename, char *description, EXTERNAL_PARA *parameter_definition,
              void **sdds_value, long n_parameters, long normalize, long verbosity, long append)
{
  SDDS_DATASET outPage;
  char buffer[100];
  long i, last_index, index;

  last_index = -1;
  if (!append) {
    /* Open file for writting */
    if (verbosity)
      fprintf( stdout, "Opening \"%s\" for writing...\n", filename);
    if (!SDDS_InitializeOutput(&outPage, SDDS_BINARY, 1, 
                               description, description, filename))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);

    for (i=0; i<n_parameters; i++) {
      if (!SDDS_ProcessParameterString(&outPage, parameter_definition[i].text, 0) ||
          (index=SDDS_GetParameterIndex(&outPage, parameter_definition[i].name))<0) {
        fprintf(stdout, "Unable to define SDDS parameter for chprint2--name was:\n%s\n",
                parameter_definition[i].name);
        fflush(stdout);
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
        exit(1);
      }
      if (index!=(last_index+1))
        fprintf(stdout, "\7\7\7WARNING: parameter indices for SDDS file %s are not sequential--this will probably cause unexpected results\n", filename);
      fflush(stdout);
      last_index = index;
    }

    if (0>SDDS_DefineParameter(&outPage, "Variable1Name", NULL, bName->xunits, 
                               NULL, NULL, SDDS_STRING, NULL) ||
        0>SDDS_DefineParameter(&outPage, "Variable1Interval", NULL, NULL, 
                               NULL, NULL, SDDS_DOUBLE, NULL) ||
        0>SDDS_DefineParameter(&outPage, "Variable1Minimum", NULL, NULL, 
                               NULL, NULL, SDDS_DOUBLE, NULL) ||
        0>SDDS_DefineParameter(&outPage, "Variable1Maximum", NULL, NULL, 
                               NULL, NULL, SDDS_DOUBLE, NULL) ||
        0>SDDS_DefineParameter(&outPage, "Variable1Dimension", NULL, NULL, 
                               NULL, NULL, SDDS_LONG, NULL) ||
        0>SDDS_DefineParameter(&outPage, "Variable2Name", NULL, bName->yunits, 
                               NULL, NULL, SDDS_STRING, NULL) ||
        0>SDDS_DefineParameter(&outPage, "Variable2Interval", NULL, NULL, 
                               NULL, NULL, SDDS_DOUBLE, NULL) ||
        0>SDDS_DefineParameter(&outPage, "Variable2Minimum", NULL, NULL, 
                               NULL, NULL, SDDS_DOUBLE, NULL) ||
        0>SDDS_DefineParameter(&outPage, "Variable2Maximum", NULL, NULL, 
                               NULL, NULL, SDDS_DOUBLE, NULL) ||
        0>SDDS_DefineParameter(&outPage, "Variable2Dimension", NULL, NULL, 
                               NULL, NULL, SDDS_LONG, NULL) ||
        0>SDDS_DefineParameter(&outPage, "Total_count", NULL, NULL, 
                               NULL, NULL, SDDS_LONG, NULL))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    
    sprintf(buffer, "%s-%sFrequency", bName->xname, bName->yname);
    if (0>SDDS_DefineColumn(&outPage, buffer, NULL, NULL, 
                            NULL, NULL, SDDS_DOUBLE, 0))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);

    if (!SDDS_WriteLayout(&outPage))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  } else {
    if (!SDDS_InitializeAppend(&outPage, filename))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }

  /* Write to output file */
  if (0>SDDS_StartPage(&outPage, bName->length) ||
      !SDDS_SetParameters(&outPage, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,
                          "Variable1Name", bName->xname,
                          "Variable2Name", bName->yname,
                          "Variable1Interval", bName->dx,
                          "Variable1Minimum", bName->xmin,
                          "Variable1Maximum", bName->xmax,
                          "Variable1Dimension", bName->xbins,
                          "Variable2Interval", bName->dy,
                          "Variable2Minimum", bName->ymin,
                          "Variable2Maximum", bName->ymax,
                          "Variable2Dimension", bName->ybins,
                          "Total_count", bName->count, NULL))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  
  for (i=0; i<n_parameters; i++) {
    if (!SDDS_SetParameters(&outPage, SDDS_SET_BY_INDEX|SDDS_PASS_BY_REFERENCE, 
                            i, sdds_value[i], -1))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }

  for (i=0; i<bName->length; i++) {
    if (normalize) {
      if (!SDDS_SetRowValues(&outPage, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, i,
                             0, bName->value[i]/bName->total, -1))
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    } else {
      if (!SDDS_SetRowValues(&outPage, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, i,
                             0, bName->value[i], -1))
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
  }

  if (!SDDS_WritePage(&outPage))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  if (!SDDS_Terminate(&outPage))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);

  return;
}

void chprintn(ntuple *bName, char *filename, char *description, EXTERNAL_PARA *parameter_definition,
              void **sdds_value, long n_parameters, long normalize, long verbosity, long append)
{
  SDDS_DATASET outPage;
  char buffer[100], dimensionString[4];
  long i, Index[5], last_index, index;
 
  last_index = -1;
  if (!append) {
    /* Open file for writting */
    if (verbosity)
      fprintf( stdout, "Opening \"%s\" for writing...\n", filename);
    if (!SDDS_InitializeOutput(&outPage, SDDS_BINARY, 1, 
                               description, description, filename))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);

    for (i=0; i<n_parameters; i++) {
      if (!SDDS_ProcessParameterString(&outPage, parameter_definition[i].text, 0) ||
          (index=SDDS_GetParameterIndex(&outPage, parameter_definition[i].name))<0) {
        fprintf(stdout, "Unable to define SDDS parameter for chprintn--name was:\n%s\n",
                parameter_definition[i].name);
        fflush(stdout);
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
        exit(1);
      }
      if (index!=(last_index+1))
        fprintf(stdout, "\7\7\7WARNING: parameter indices for SDDS file %s are not sequential--this will probably cause unexpected results\n", filename);
      fflush(stdout);
      last_index = index;
    }

    for (i=0; i< bName->nD; i++) {
      sprintf(dimensionString, "%02ld", i);
      
      sprintf(buffer, "Variable%sName", dimensionString);
      if (0>SDDS_DefineParameter(&outPage, buffer, NULL, bName->units[i], 
                                 NULL, NULL, SDDS_STRING, NULL)) 
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);

      sprintf(buffer, "Variable%sMin", dimensionString);
      if (0>SDDS_DefineParameter(&outPage, buffer, NULL, NULL, 
                                 NULL, NULL, SDDS_DOUBLE, NULL)) 
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);

      sprintf(buffer, "Variable%sMax", dimensionString);
      if (0>SDDS_DefineParameter(&outPage, buffer, NULL, NULL, 
                                 NULL, NULL, SDDS_DOUBLE, NULL)) 
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);

      sprintf(buffer, "Variable%sInterval", dimensionString);
      if (0>SDDS_DefineParameter(&outPage, buffer, NULL, NULL, 
                                 NULL, NULL, SDDS_DOUBLE, NULL)) 
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);

      sprintf(buffer, "Variable%sDimension", dimensionString);
      if (0>SDDS_DefineParameter(&outPage, buffer, NULL, NULL, 
                                 NULL, NULL, SDDS_LONG, NULL)) 
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }

    if (0>SDDS_DefineParameter(&outPage, "Total_count", NULL, NULL, 
                               NULL, NULL, SDDS_LONG, NULL))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    if (0>SDDS_DefineColumn(&outPage, "Index", NULL, NULL, 
                            NULL, NULL, SDDS_LONG, 0) ||
        0>SDDS_DefineColumn(&outPage, "Frequency", NULL, NULL, 
                            NULL, NULL, SDDS_DOUBLE, 0))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    
    if (!SDDS_WriteLayout(&outPage))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  } else {
    if (!SDDS_InitializeAppend(&outPage, filename))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }

  /* Write to output file */
  if (0>SDDS_StartPage(&outPage, bName->length)) 
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);

  for (i=0; i<n_parameters; i++) {
    if (!SDDS_SetParameters(&outPage, SDDS_SET_BY_INDEX|SDDS_PASS_BY_REFERENCE, 
                            i, sdds_value[i], -1))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }

  Index[4] = n_parameters-1;
  for (i=0; i< bName->nD; i++) {
    Index[0] = Index[4]+1; Index[1]=Index[0]+1; Index[2]=Index[1]+1;
    Index[3]=Index[2]+1; Index[4]=Index[3]+1;
    if (!SDDS_SetParameters(&outPage, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE,
                            Index[0], bName->vname[i],
                            Index[1], bName->xmin[i],
                            Index[2], bName->xmax[i],
                            Index[3], bName->dx[i],
                            Index[4], bName->xbins[i], -1))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  if (!SDDS_SetParameters(&outPage, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,
                          "Total_count", bName->count, NULL))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  for (i=0; i<bName->length; i++) {
    if (normalize) {
      if (!SDDS_SetRowValues(&outPage, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, i,
                             0, i, 1, bName->value[i]/bName->total, -1))
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    } else {
      if (!SDDS_SetRowValues(&outPage, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, i,
                             0, i, 1, bName->value[i], -1))
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
  }
  if (!SDDS_WritePage(&outPage))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  if (!SDDS_Terminate(&outPage))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);

  return;
}

void chprint1m(book1m *bName, char *filename, char *description, EXTERNAL_PARA *parameter_definition, 
               void **sdds_value, long n_parameters, long normalize, long verbosity, long append)
{
  SDDS_DATASET outPage;
  char name[100], units[100], freq[100];
  long i, j, Index[2], last_index, index;
  double *value1, *value2;

  last_index = -1;
  if (!append) {
    /* Open file for writting */
    if (verbosity)
      fprintf( stdout, "Opening \"%s\" for writing...\n", filename);
    if (!SDDS_InitializeOutput(&outPage, SDDS_BINARY, 1, 
                               description, description, filename))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);

    for (i=0; i<n_parameters; i++) {
      if (!SDDS_ProcessParameterString(&outPage, parameter_definition[i].text, 0) ||
          (index=SDDS_GetParameterIndex(&outPage, parameter_definition[i].name))<0) {
        fprintf(stdout, "Unable to define SDDS parameter for chprint1m--name was:\n%s\n",
                parameter_definition[i].name);
        fflush(stdout);
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
        exit(1);
      }
      if (index!=(last_index+1))
        fprintf(stdout, "\7\7\7WARNING: parameter indices for SDDS file %s are not sequential--this will probably cause unexpected results\n", filename);
      fflush(stdout);
      last_index = index;
    }

    for (i=0; i< bName->nD; i++) {
      sprintf(name, "%s", bName->vname[i]);
      sprintf(units, "%s", bName->units[i]);
      sprintf(freq, "%sFrequency", bName->vname[i]);
      if (0>SDDS_DefineColumn(&outPage, name, NULL, units,
                               NULL, NULL, SDDS_DOUBLE, 0) ||
          0>SDDS_DefineColumn(&outPage, freq, NULL, NULL, 
                              NULL, NULL, SDDS_DOUBLE, 0))
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
    if (!SDDS_WriteLayout(&outPage) )
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  } else {
    if (!SDDS_InitializeAppend(&outPage, filename))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }

  /* Write to output file */
  if (0>SDDS_StartPage(&outPage, bName->length))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  for (i=0; i<n_parameters; i++) {
    if (!SDDS_SetParameters(&outPage, SDDS_SET_BY_INDEX|SDDS_PASS_BY_REFERENCE, 
                            i, sdds_value[i], -1))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }

  Index[1] = -1;
  value1 = calloc(sizeof(*value1), bName->length);
  value2 = calloc(sizeof(*value2), bName->length);
  for (i=0; i<bName->nD; i++){
    sprintf(name, "%s", bName->vname[i]);
    sprintf(freq, "%sFrequency", bName->vname[i]);
    for (j=0; j<bName->length; j++){
      value1[j] = ((double)j+0.5)*bName->dx[i]+bName->xmin[i];
      value2[j] = bName->value[j][i];
      if (normalize) 
        value2[j] /= bName->total;      
    }
    Index[0] = Index[1]+1; Index[1]=Index[0]+1;
    for (j=0; j<bName->length; j++){
      if (!SDDS_SetRowValues(&outPage, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, 
                             j, name, value1[j], freq, value2[j], NULL))
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);      
    }
    /*
    if (!SDDS_SetColumn(&outPage, SDDS_SET_BY_INDEX, value1, bName->length, Index[0]) ||
        !SDDS_SetColumn(&outPage, SDDS_SET_BY_INDEX, value2, bName->length, Index[1]))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);      
    */
  }
  if (!SDDS_WritePage(&outPage))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  if (!SDDS_Terminate(&outPage))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);

  free(value1); free(value2);
  return;
}

