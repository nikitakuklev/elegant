/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* prototypes for chbook.c */

typedef struct {
  char *vname, *units;
  double xmin, xmax, dx;
  long xbins;
  double *value, total;
  long length, count;
} book1;

typedef struct {
  char *xname, *yname, *xunits, *yunits;
  double xmin, xmax, dx;
  long xbins;
  double ymin, ymax, dy;
  long ybins;
  double *value, total;
  long length, count;
} book2;

typedef struct {
  char **vname, **units;
  long nD;
  double *xmin, *xmax, *dx;
  long *xbins;
  double *value, total;
  long length, count;
} ntuple;

typedef struct {
  char **vname, **units;
  double *xmin, *xmax, *dx;
  long nD;
  long bins;
  double **value, total;
  long length, count;
} book1m;

book1 *chbook1(char *vName, char *units, double xmin, double xmax, long xbins);
void chfill1(book1 *bName, double x, double Frequency);
void free_hbook1(book1 *x);

book2 *chbook2(char *xName, char *yName, char *xunits, char *yunits,
               double xmin, double xmax, double ymin, double ymax, long xbins, long ybins);
void chfill2(book2 *bName, double x, double y, double Frequency);
void free_hbook2(book2 *x);

ntuple *chbookn(char **vName, char **units, long NDimension, double *xmin, double *xmax, long *xbins, long offset);
void chfilln(ntuple *bName, double *x, double Frequency, long offset);
void free_hbookn(ntuple *x);

book1m *chbook1m(char **vName, char **units, double *xmin, double *xmax, long *xbins, long column_number);
void chfill1m(book1m *bName, double *x, double Frequency, long *xbins, long column_number);
void free_hbook1m(book1m *x);

typedef struct {
    char *name, *text; 
} EXTERNAL_PARA;
void chprint1(book1 *bName, char *filename, char *description, EXTERNAL_PARA *parameter_definition,
              void **sdds_value, long n_parameters, long normalize, long verbosity, long append);
void chprint2(book2 *bName, char *filename, char *description, EXTERNAL_PARA *parameter_definition,
              void **sdds_value, long n_parameters, long normalize, long verbosity, long append);
void chprintn(ntuple *bName, char *filename, char *description, EXTERNAL_PARA *parameter_definition,
              void **sdds_value, long n_parameters, long normalize, long verbosity, long append);
void chprint1m(book1m *bName, char *filename, char *description, EXTERNAL_PARA *parameter_definition, 
               void **sdds_value, long n_parameters, long normalize, long verbosity, long append);
