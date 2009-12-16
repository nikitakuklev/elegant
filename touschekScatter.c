/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

#include "mdb.h"
#include "track.h"
#include "match_string.h"
#include "touschekScatter.h"
#include "SDDS.h"
#include "constants.h"

#define DEBUG 0
#define MAXREGION 30
#define STEPS 100

/* global variable */
static TSCATTER_SPEC *tsSpec = NULL;
static long distIn = 0;
static SDDS_DATASET tranPage, longPage;
static double tranMin[4], tranMax[4], tranInterval[4], tranDelta[4], *tweight;
static double longMin[2], longMax[2], longInterval[2], longDelta[2], *lweight;
static long tranDimen[4], longDimen[2];

/* Initialization of the simulation */
int init_distInput();
void init_TSPEC (RUN *run);
TSCATTER *initTSCATTER (ELEMENT_LIST *eptr);
/* Calculate local and integrated Touschek scattering rate */
int TouschekRate(LINE_LIST *beamline);
void FIntegral(double tm, double B1, double B2, double *F, char *name);
double Fvalue (double t, double tm, double b1, double b2);
/* Monte Carlo simulation of Touschek scattering */
void TouschekDistribution (RUN *run, LINE_LIST *beamline);
void TouschekDistributionE(RUN *run, LINE_LIST *beamline);

void selectPart0(TSCATTER *tsptr, double *p1, double *p2, 
                double *dens1, double *dens2, double *ran1);
void selectPart1(TSCATTER *tsptr, double *p1, double *p2, 
                double *dens1, double *dens2, double *ran1);
double findDens(double *p);

void bunch2cm(double *p1, double *p2, double *q, double *beta, double *gamma);
void eulertrans(double *v0, double theta, double phi, double *v1, double *v);
void cm2bunch(double *p1, double *p2, double *q, double *beta, double *gamma);
double moeller(double beta0, double theta);
void pickPart(double *weight, long *index, long start, long end, 
              long *iTotal, double *wTotal, double weight_limit, double weight_ave);

char *compose_filename1(char *template, char *root_name, long index);
int str_inp(const char *string1, const char *string2);
void print_hbook (book1 *x, book1 *y, book1 *s, book1 *xp, book1 *yp, book1 *dp,
                  TSCATTER *tsptr, char *filename, char *description, int verbosity, double factor);

void setupTouschekEffect(NAMELIST_TEXT *nltext, RUN *run, LINE_LIST *beamline) 
{
  ELEMENT_LIST *eptr;
  int flag;

  /* Check if there is TScatter element along beamline. */
  eptr = &(beamline->elem);
  flag = 0;
  while (eptr) {
    if (eptr->type == T_TSCATTER) {
      flag = 1;
      break;
    }
    eptr = eptr->succ; 
  }
  if(!flag) 
    bomb("No TSCATTER element along beamline", NULL);                

   /* process the namelist text */
  set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
  set_print_namelist_flags(0);
  process_namelist(&touschek_scatter, nltext);
  if (echoNamelists) print_namelist(stdout, &touschek_scatter);

  if (!charge)    
    bomb("charge has to be given", NULL);
  if (frequency<1)
    bomb("frequency has to >=1", NULL);
  if (!emit_x && !emit_nx)
    bomb("bunch emittance-x has to be given", NULL);
  if (!emit_nx)
    emit_nx = sqrt(1.+ sqr(run->p_central))*emit_x;
  if (!emit_y && !emit_ny)
    bomb("bunch emittance-y has to be given", NULL);
  if (!emit_ny)
    emit_ny = sqrt(1.+ sqr(run->p_central))*emit_y;
  if (!sigma_dp)
    bomb("energy spread has to be given", NULL);
  if (!sigma_s)
    bomb("bunch length has to be given", NULL);
  if (!delta && !delta_mev)
    bomb("delta and delta_mev can not be zero at the same time", NULL);
  if (!delta_mev)
    delta_mev = (sqrt(sqr((delta+1.)*run->p_central)+1.)-sqrt(sqr(run->p_central)+1.))*me_mev;
  if (longDist || tranDist) {
    distIn = 1;
    if (!longDist || !tranDist)
      bomb("both tranDist and longDist need to be given", NULL);
    init_distInput();    
  }
  
  /* Initial Touschek calculation then calculate the Twiss function. */
  init_TSPEC (run);
  run_twiss_output(run, beamline, NULL, -1);
#if USE_MPI
  MPI_Barrier(MPI_COMM_WORLD); 
#endif
  /* Calculate Piwinski's scattering rate. */
  if (TouschekRate(beamline))
    bomb("Touschek scattering rate calculation error", NULL);
  /* Generate scattered particles at TScatter. 
     And track the scattered particles down to beamline. */
#if USE_MPI
  MPI_Barrier(MPI_COMM_WORLD); 
#endif
  if (estimate) 
    TouschekDistributionE(run, beamline);
  else
    TouschekDistribution(run, beamline);
    
#if USE_MPI
  MPI_Barrier(MPI_COMM_WORLD); 
#endif
  return;
}

int TouschekRate(LINE_LIST *beamline)
{
  double NP;
  double tm, B1, B2, F, rate, IntR, IntLength;
  ELEMENT_LIST *eptr;

  double betagamma, gamma; 
  double deltap0, sp2, sp4, beta2, betagamma2, sh2;
  double sx2, sxb2, dx2, dx_2;
  double sy2, syb2, dy2, dy_2;
  double a0, c0, c1, c2, c3;
  TWISS *twiss0;

  NP = tsSpec->charge/e_mks;
  a0 = sqr(re_mks)*c_mks*NP*NP/(8*sqrt(PI)*sigma_s);

  IntR = 0.;
  IntLength = 0.;
  rate = 0.;
  eptr = &(beamline->elem);
  while (eptr) {
    if (eptr->type == T_TSCATTER) {
      ((TSCATTER*)eptr->p_elem)->AveR = IntR/IntLength;
      ((TSCATTER*)eptr->p_elem)->p_rate = rate;
      ((TSCATTER*)eptr->p_elem)->total_scatter = IntR / c_mks * tsSpec->frequency;
      IntR = 0.;
      IntLength = 0.;
    }
    if(!(entity_description[eptr->type].flags&HAS_LENGTH) ||
       (eptr->pred && !(((DRIFT*)eptr->p_elem)->length)) ||
       !eptr->succ) {
      eptr = eptr->succ;
      continue;
    }
    betagamma = eptr->Pref_output;
    gamma = sqrt(sqr(betagamma)+1.);
    beta2 = sqr(betagamma)/sqr(gamma);
    betagamma2 = sqr(betagamma);

    deltap0 = sqrt(sqr(gamma+tsSpec->delta_mev/me_mev)-1.)-betagamma;
    tm = beta2*sqr(deltap0/betagamma);
    sp2 = sqr(tsSpec->delta_p0/betagamma);
    sp4 = sqr(sp2);

    twiss0 = eptr->twiss;
    sxb2 = twiss0->betax*tsSpec->emitN[0]/gamma;
    dx2  = sqr(twiss0->etax);
    sx2  = sxb2 + dx2*sp2;
    dx_2 = sqr(twiss0->alphax * twiss0->etax + twiss0->betax * twiss0->etapx);

    syb2 = twiss0->betay*tsSpec->emitN[1]/gamma;
    dy2  = sqr(twiss0->etay);
    sy2  = syb2 + dy2*sp2;
    dy_2 = sqr(twiss0->alphay * twiss0->etay + twiss0->betay * twiss0->etapy);

    sh2 = 1/(1/sp2+(dx2+dx_2)/sxb2+(dy2+dy_2)/syb2);
    c0 = sh2/(sqr(beta2*tsSpec->emitN[0]*tsSpec->emitN[1])*sp2);
    c1 = sqr(twiss0->betax)/(2*betagamma2*sxb2);
    c2 = sqr(twiss0->betay)/(2*betagamma2*syb2);
    c3 = sx2*sy2-sp4*dx2*dy2;
    B1 = c1*(1-sh2*dx_2/sxb2)+c2*(1-sh2*dy_2/syb2);
    B2 = sqr(B1)-c0*c3;   	  
    if (B2<0) {
      fprintf(stdout, "B2^2<0 at \"%s\" occurence %ld", eptr->name, eptr->occurence);
      exit(1);
    }
    B2=sqrt(B2);   	  

    FIntegral(tm, B1, B2, &F, eptr->name);

    rate = a0*sqrt(c0)*F/gamma/gamma;
    IntR += rate * ((DRIFT*)eptr->p_elem)->length;
    IntLength +=  ((DRIFT*)eptr->p_elem)->length;
    eptr = eptr->succ; 
  }

  return(0);
}

void FIntegral(double tm, double b1, double b2, double *F, char *name) 
{
  long maxRegion = MAXREGION;
  long steps = STEPS; /* number of integration steps per decade */
  long i, j;
  double test = 1e-5, simpsonCoeff[2] = {2.,4.};
  double h, tau0, tau, intF;
  double cof, sum, fun;
  /* Using simpson's rule to do the integral. 
     split integral into region with "steps" steps per region.
     integration intervals be [1,3], [3,9], [9,27], and so on. 
     i is the index over these intervals. j is the index over each step.
  */
  tau0 = tm;
  intF = 0.0;
  for (i=0; i<maxRegion; i++) {
    if (i==0) {
      h = tau0 * 2. / steps;
    } else {
      h = tau0 * 3. / steps;
    }
    sum = 0.0;
    for (j=0; j<=steps; j++) {
      tau = tau0 + h*j;
      cof = simpsonCoeff[j%2]; 
      if (j==0 || j==steps) 
        cof = 1.;
      fun = Fvalue(tau, tm, b1, b2);
      sum += cof*fun;
    }
    tau0 = tau;
    sum = (sum/3.0)*h;
    intF += sum;
    if (FABS(sum/intF)<test) 
      break;
    if ( i== maxRegion) 
      fprintf( stdout, "**Warning** Integral did not converge till tau= %g.\n",tau);
  }

  *F = intF;
  return;
}

double Fvalue (double t, double tm, double b1, double b2)
{
  double c0, c1, c2, result;

  c0 = (sqr(2.+1./t)*(t/tm/(1.+t)-1.)+1.-sqrt(tm*(1.+t)/t)-(4.+1./t)/(2.*t)*log(t/tm/(1.+t)))*sqrt(t/(1.+t));
  c1 = exp(-b1*t);
  c2 = dbesi0(b2*t);
  result = c0 * c1 * c2;
  /* If overflow/underflow use approximate equation for modified bessel function. */
  if (isnan(result) || result>FLT_MAX) {
    result=c0*exp(b2*t-b1*t)/sqrt(PIx2*b2*t);
  } 
  return result;
}

char *compose_filename1(char *template, char *root_name, long index)
{
  char *root, *ext, *filename, *inx;
  long i, divp;

  divp = str_inp(template, ".");
  if (divp<0)
    bomb("touschek_scatter: The filename has to be x.x format. x is 1 to anynumber of characters", NULL);
  if (str_in(template, "%s")) {
    root = root_name;
  }
  else {
    root = tmalloc(sizeof(char)*divp);
    for (i=0; i<divp; i++) {
      *(root+i) = *(template+i);
    }
    *(root+i) = '\0';    
  }

  ext =  tmalloc(sizeof(char)*(strlen(template)-divp+1));
  for(i=0; i<strlen(template)-divp; i++) {
    *(ext+i) = *(template+i+divp);
  }
  *(ext+i) = '\0';

  inx = tmalloc(sizeof(char)*5);
  sprintf(inx, "%04ld", index);

  filename =  tmalloc(sizeof(char)*(strlen(root)+strlen(inx)+strlen(ext)+1));
  strcat(filename,root);
  strcat(filename,inx);
  strcat(filename,ext);
  return(filename);
}

int str_inp(const char *string1, const char *string2)
{
  int i, j, flag;

  if (strlen(string1)<strlen(string2))
    return (-1);

  i=j=0;
  for (i=0;i<strlen(string1); i++) {
    if (*(string1+i) == *string2) {
      if (strlen(string2)+i>strlen(string1))
        return (-1);
      for (j=0; j<strlen(string2); j++) {
        if (*(string1+i+j) != *(string2+j)) {
          flag = 0;
          break;
        }
        flag=1;
      }
      if(flag)
        return (i);
    }
  }

  return (-1);
}

void print_hbook (book1 *x, book1 *y, book1 *s, book1 *xp, book1 *yp, book1 *dp,
                  TSCATTER *tsptr, char *filename, char *description, int verbosity, double factor)
{
  SDDS_DATASET outPage;
  long i;
  double *v1, *v2, *v3, *v4, *v5, *v6;

  v1 = calloc(sizeof(*v1), x->length);
  v2 = calloc(sizeof(*v2), y->length);
  v3 = calloc(sizeof(*v3), s->length);
  v4 = calloc(sizeof(*v4), xp->length);
  v5 = calloc(sizeof(*v5), yp->length);
  v6 = calloc(sizeof(*v6), dp->length);
  for(i=0; i<x->length; i++) {
    v1[i] = ((double)i+0.5)*x->dx + x->xmin;
    v2[i] = ((double)i+0.5)*y->dx + y->xmin;
    v3[i] = ((double)i+0.5)*s->dx + s->xmin;
    v4[i] = ((double)i+0.5)*xp->dx+ xp->xmin;
    v5[i] = ((double)i+0.5)*yp->dx+ yp->xmin;
    v6[i] = ((double)i+0.5)*dp->dx+ dp->xmin;
    x->value[i] *= factor;
    y->value[i] *= factor;
    s->value[i] *= factor;
    xp->value[i] *= factor;
    yp->value[i] *= factor;
    dp->value[i] *= factor;
  }

  /* Open file for writting */
  if (verbosity)
    fprintf( stdout, "Opening \"%s\" for writing...\n", filename);

  if (!SDDS_InitializeOutput(&outPage, SDDS_BINARY, 1, 
                             description, description, filename))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);

  if (0>SDDS_DefineParameter(&outPage, "Element_Name", NULL, NULL, 
                             NULL, NULL, SDDS_STRING, NULL) ||
      0>SDDS_DefineParameter(&outPage, "Element_Location", NULL, NULL, 
                             NULL, NULL, SDDS_DOUBLE, NULL) ||
      0>SDDS_DefineParameter(&outPage, "Piwinski_AveR", NULL, NULL, 
                             NULL, NULL, SDDS_DOUBLE, NULL) ||
      0>SDDS_DefineParameter(&outPage, "Piwinski_Rate", NULL, NULL, 
                             NULL, NULL, SDDS_DOUBLE, NULL) ||
      0>SDDS_DefineParameter(&outPage, "MC_Rate", NULL, NULL, 
                             NULL, NULL, SDDS_DOUBLE, NULL) ||
      0>SDDS_DefineParameter(&outPage, "Ignored_Rate", NULL, NULL, 
                             NULL, NULL, SDDS_DOUBLE, NULL) ||
      0>SDDS_DefineParameter(&outPage, "nbins", NULL, NULL, 
                             NULL, NULL, SDDS_LONG, NULL) ||
      0>SDDS_DefineParameter(&outPage, "Total_count", NULL, NULL, 
                             NULL, NULL, SDDS_LONG, NULL) ||
      0>SDDS_DefineParameter(&outPage, "VariableName_1", NULL, NULL, 
                             NULL, NULL, SDDS_STRING, NULL) ||
      0>SDDS_DefineParameter(&outPage, "Interval_1", NULL, NULL, 
                             NULL, NULL, SDDS_DOUBLE, NULL) ||
      0>SDDS_DefineParameter(&outPage, "Minimum_1", NULL, NULL, 
                             NULL, NULL, SDDS_DOUBLE, NULL) ||
      0>SDDS_DefineParameter(&outPage, "VariableName_2", NULL, NULL, 
                             NULL, NULL, SDDS_STRING, NULL) ||
      0>SDDS_DefineParameter(&outPage, "Interval_2", NULL, NULL, 
                             NULL, NULL, SDDS_DOUBLE, NULL) ||
      0>SDDS_DefineParameter(&outPage, "Minimum_2", NULL, NULL, 
                             NULL, NULL, SDDS_DOUBLE, NULL) ||
      0>SDDS_DefineParameter(&outPage, "VariableName_3", NULL, NULL, 
                             NULL, NULL, SDDS_STRING, NULL) ||
      0>SDDS_DefineParameter(&outPage, "Interval_3", NULL, NULL, 
                             NULL, NULL, SDDS_DOUBLE, NULL) ||
      0>SDDS_DefineParameter(&outPage, "Minimum_3", NULL, NULL, 
                             NULL, NULL, SDDS_DOUBLE, NULL) ||
      0>SDDS_DefineParameter(&outPage, "VariableName_4", NULL, NULL, 
                             NULL, NULL, SDDS_STRING, NULL) ||
      0>SDDS_DefineParameter(&outPage, "Interval_4", NULL, NULL, 
                             NULL, NULL, SDDS_DOUBLE, NULL) ||
      0>SDDS_DefineParameter(&outPage, "Minimum_4", NULL, NULL, 
                             NULL, NULL, SDDS_DOUBLE, NULL) ||
      0>SDDS_DefineParameter(&outPage, "VariableName_5", NULL, NULL, 
                             NULL, NULL, SDDS_STRING, NULL) ||
      0>SDDS_DefineParameter(&outPage, "Interval_5", NULL, NULL, 
                             NULL, NULL, SDDS_DOUBLE, NULL) ||
      0>SDDS_DefineParameter(&outPage, "Minimum_5", NULL, NULL, 
                             NULL, NULL, SDDS_DOUBLE, NULL) ||
      0>SDDS_DefineParameter(&outPage, "VariableName_6", NULL, NULL, 
                             NULL, NULL, SDDS_STRING, NULL) ||
      0>SDDS_DefineParameter(&outPage, "Interval_6", NULL, NULL, 
                             NULL, NULL, SDDS_DOUBLE, NULL) ||
      0>SDDS_DefineParameter(&outPage, "Minimum_6", NULL, NULL, 
                             NULL, NULL, SDDS_DOUBLE, NULL))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);

  if (0>SDDS_DefineColumn(&outPage, "x", NULL, NULL, 
                          NULL, NULL, SDDS_DOUBLE, 0) ||
      0>SDDS_DefineColumn(&outPage, "x_count", NULL, NULL, 
                          NULL, NULL, SDDS_DOUBLE, 0) ||
      0>SDDS_DefineColumn(&outPage, "y", NULL, NULL, 
                          NULL, NULL, SDDS_DOUBLE, 0) ||
      0>SDDS_DefineColumn(&outPage, "y_count", NULL, NULL, 
                          NULL, NULL, SDDS_DOUBLE, 0) ||
      0>SDDS_DefineColumn(&outPage, "s", NULL, NULL, 
                          NULL, NULL, SDDS_DOUBLE, 0) ||
      0>SDDS_DefineColumn(&outPage, "s_count", NULL, NULL, 
                          NULL, NULL, SDDS_DOUBLE, 0) ||
      0>SDDS_DefineColumn(&outPage, "xp", NULL, NULL, 
                          NULL, NULL, SDDS_DOUBLE, 0) ||
      0>SDDS_DefineColumn(&outPage, "xp_count", NULL, NULL, 
                          NULL, NULL, SDDS_DOUBLE, 0) ||
      0>SDDS_DefineColumn(&outPage, "yp", NULL, NULL, 
                          NULL, NULL, SDDS_DOUBLE, 0) ||
      0>SDDS_DefineColumn(&outPage, "yp_count", NULL, NULL, 
                          NULL, NULL, SDDS_DOUBLE, 0) ||
      0>SDDS_DefineColumn(&outPage, "dp/p", NULL, NULL, 
                          NULL, NULL, SDDS_DOUBLE, 0) ||
      0>SDDS_DefineColumn(&outPage, "dp/p_count", NULL, NULL, 
                          NULL, NULL, SDDS_DOUBLE, 0))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);

  if (!SDDS_WriteLayout(&outPage) )
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);

  /* Write to output file */
  if (0>SDDS_StartPage(&outPage, x->length) ||
      !SDDS_SetParameters(&outPage, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,
                          "Element_Name", tsptr->name,
                          "Element_Location", tsptr->s, 
                          "VariableName_1", x->vname,
                          "Interval_1", x->dx,
                          "Minimum_1", x->xmin,
                          "VariableName_2", y->vname,
                          "Interval_2", y->dx,
                          "Minimum_2", y->xmin,
                          "VariableName_3", s->vname,
                          "Interval_3", s->dx,
                          "Minimum_3", s->xmin,
                          "VariableName_4", xp->vname,
                          "Interval_4", xp->dx,
                          "Minimum_4", xp->xmin,
                          "VariableName_5", yp->vname,
                          "Interval_5", yp->dx,
                          "Minimum_5", yp->xmin,
                          "VariableName_6", dp->vname,
                          "Interval_6", dp->dx,
                          "Minimum_6", dp->xmin,
                          "Piwinski_AveR", tsptr->AveR,
                          "Piwinski_Rate", tsptr->p_rate,
                          "MC_Rate", tsptr->s_rate,
                          "Ignored_Rate", tsptr->i_rate,
                          "nbins", x->length,
                          "Total_count", x->count, NULL) ||
      !SDDS_SetColumn(&outPage, SDDS_SET_BY_NAME, v1,        x->length,  "x") ||
      !SDDS_SetColumn(&outPage, SDDS_SET_BY_NAME, x->value,  x->length,  "x_count") ||
      !SDDS_SetColumn(&outPage, SDDS_SET_BY_NAME, v2,        y->length,  "y") ||
      !SDDS_SetColumn(&outPage, SDDS_SET_BY_NAME, y->value,  y->length,  "y_count") ||
      !SDDS_SetColumn(&outPage, SDDS_SET_BY_NAME, v3,        s->length,  "s") ||
      !SDDS_SetColumn(&outPage, SDDS_SET_BY_NAME, s->value,  s->length,  "s_count") ||
      !SDDS_SetColumn(&outPage, SDDS_SET_BY_NAME, v4,        xp->length, "xp") ||
      !SDDS_SetColumn(&outPage, SDDS_SET_BY_NAME, xp->value, xp->length, "xp_count") ||
      !SDDS_SetColumn(&outPage, SDDS_SET_BY_NAME, v5,        yp->length, "yp") ||
      !SDDS_SetColumn(&outPage, SDDS_SET_BY_NAME, yp->value, yp->length, "yp_count") ||
      !SDDS_SetColumn(&outPage, SDDS_SET_BY_NAME, v6,        dp->length, "dp/p") ||
      !SDDS_SetColumn(&outPage, SDDS_SET_BY_NAME, dp->value, dp->length, "dp/p_count") ||
      !SDDS_WritePage(&outPage))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);

  if (!SDDS_Terminate(&outPage))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);

  return;
}

void TouschekDistribution(RUN *run, LINE_LIST *beamline)
{
  long i, j, total_event, n_left, eleCount=0;
  ELEMENT_LIST *eptr;
  TSCATTER *tsptr;
  double p1[6], p2[6], dens1, dens2, tdens;
  double theta, phi, qa[3], qb[3], beta[3], qabs, gamma;
  double beta0, cross, temp;
  static SDDS_TABLE SDDS_bunch, SDDS_loss;
  BEAM  Beam0, Beam, *beam0, *beam;
  double xrangeP[6], xrangeM[6];
  /*
  book1 *x, *y, *s, *xp, *yp, *dp;
  book1 *x0, *y0, *s0, *xp0, *yp0, *dp0;
  */
  book1m *iniBook, *outBook;
  char *Name[6]={"x","xp","y","yp","dt","dp"};
  char *Units[6]={"m","","m","","s",""};
  long bookBins[6];
  book1 *lossDis;
  double *weight;
  double ran1[11];
  long *index, iTotal, sTotal;
  double weight_limit, weight_ave, wTotal;
 
  eptr = &(beamline->elem);
  beam0 = &Beam0;
  beam = &Beam;
  sTotal = (int)beamline->revolution_length;

  while (eptr) {
    if (eptr->type == T_TSCATTER) {
      eleCount++;
      if(eleCount < i_start) {
        eptr = eptr->succ; 
        continue;
      }
      if(eleCount > i_end)
        break;
#if USE_MPI
      MPI_Barrier(MPI_COMM_WORLD); 
      if (myid==0) {
#endif
        tsptr = initTSCATTER (eptr);
        weight = (double*)malloc(sizeof(double)*n_simulated);
        beam0->particle = (double**)czarray_2d(sizeof(double), n_simulated, 7);
        beam0->original = (double**)czarray_2d(sizeof(double), n_simulated, 7);
        beam0->accepted = NULL;
        beam0->n_original = beam0->n_to_track = beam0->n_particle = n_simulated;
        beam0->n_accepted = beam0->n_saved = 0;
        beam0->p0_original = beam0->p0 = tsptr->betagamma;
        beam0->bunchFrequency = 0.;
        beam0->lostOnPass = tmalloc(sizeof(*(beam0->lostOnPass))*beam0->n_to_track);
        if (!distIn) {
          xrangeM[0] = -(xrangeP[0] = 0.5*tsSpec->range[0]*tsptr->sigx);
          xrangeM[1] = -(xrangeP[1] = 0.5*tsSpec->range[1]*tsptr->sigy);
          xrangeM[2] = -(xrangeP[2] = 0.5*tsSpec->range[2]*tsptr->sigz);
          xrangeM[3] = -(xrangeP[3] = 0.5*tsSpec->range[0]
                         *sqrt(tsSpec->emitN[0]/tsptr->gamma/tsptr->twiss[0][1])*(1.+fabs(tsptr->twiss[0][0])));
          xrangeM[4] = -(xrangeP[4] = 0.5*tsSpec->range[1]
                         *sqrt(tsSpec->emitN[1]/tsptr->gamma/tsptr->twiss[1][1])*(1.+fabs(tsptr->twiss[1][0])));
	
          xrangeM[0] = -(xrangeP[0] += fabs(0.5*tsSpec->range[2]*tsSpec->delta_p0/tsptr->betagamma*tsptr->disp[0][0]));
          xrangeM[1] = -(xrangeP[1] += fabs(0.5*tsSpec->range[2]*tsSpec->delta_p0/tsptr->betagamma*tsptr->disp[1][0]));
          xrangeM[3] = -(xrangeP[3] += fabs(0.5*tsSpec->range[2]*tsSpec->delta_p0/tsptr->betagamma*tsptr->disp[0][1]));
          xrangeM[4] = -(xrangeP[4] += fabs(0.5*tsSpec->range[2]*tsSpec->delta_p0/tsptr->betagamma*tsptr->disp[1][1]));
          xrangeM[5] = -(xrangeP[5] = 0.5*tsSpec->range[2]*tsSpec->delta_p0/tsptr->betagamma);
        } else {
          double dpmax = fabs(longMax[1])>fabs(longMin[1])? fabs(longMax[1]) : fabs(longMin[1]); 
          xrangeP[0] = tranMax[0] + fabs(dpmax*tsptr->disp[0][0]);
          xrangeM[0] = tranMin[0] - fabs(dpmax*tsptr->disp[0][0]);
          xrangeP[1] = tranMax[2] + fabs(dpmax*tsptr->disp[1][0]);
          xrangeM[1] = tranMin[2] - fabs(dpmax*tsptr->disp[1][0]);
          xrangeP[3] = tranMax[1] + fabs(dpmax*tsptr->disp[0][1]);
          xrangeM[3] = tranMin[1] - fabs(dpmax*tsptr->disp[0][1]);
          xrangeP[4] = tranMax[3] + fabs(dpmax*tsptr->disp[1][1]);
          xrangeM[4] = tranMin[3] - fabs(dpmax*tsptr->disp[1][1]);
          xrangeP[2] = longMax[0];
          xrangeM[2] = longMin[0];
          xrangeP[5] = longMax[1];
          xrangeM[5] = longMin[1];
        }
        if (initial) {
          tdens = 0.;
          tsptr->iniFile = compose_filename1(initial, run->rootname, eleCount);
          for (i=0; i<6; i++)
            bookBins[i] = tsSpec->nbins;
          iniBook = chbook1m(Name, Units, xrangeM, xrangeP, bookBins, 6);
        }
        if (distribution) {
          tsptr->disFile = compose_filename1(distribution, run->rootname, eleCount);
          for (i=0; i<6; i++)
            bookBins[i] = tsSpec->nbins;
          xrangeM[5] = -0.1;
          xrangeP[5] = 0.1;
          outBook = chbook1m(Name, Units, xrangeM, xrangeP, bookBins, 6);
        }
        if (output) {
          tsptr->outFile = compose_filename1(output, run->rootname, eleCount);
          lossDis = chbook1("los_distribution", "m", 0, sTotal, sTotal);
        }
        if (loss) {
          tsptr->losFile = compose_filename1(loss, run->rootname, eleCount);
          SDDS_BeamScatterLossSetup(&SDDS_loss, tsptr->losFile, SDDS_BINARY, 1, 
                                    "lost particle coordinates", run->runfile,
                                    run->lattice, "touschek_scatter");
        }
        if (bunch) {
          tsptr->bunFile = compose_filename1(bunch, run->rootname, eleCount);
          SDDS_BeamScatterSetup(&SDDS_bunch, tsptr->bunFile, SDDS_BINARY, 1, 
                                "scattered-beam phase space", run->runfile,
                                run->lattice, "touschek_scatter");
        }
        i = 0; j=0; total_event=0;
        while(1) {
          if(total_event>=n_simulated)
            break;
          /* Select the 11 random number then mix them. Use elegant run_setup seed */
          for (j=0; j<11; j++) {
            ran1[j] = random_1_elegant(1);
          }
          randomizeOrder((char*)ran1, sizeof(ran1[0]), 11, 0, random_4);
          
          total_event++;
          if (!distIn)
            selectPart0(tsptr, p1, p2, &dens1, &dens2, ran1);
          else
            selectPart1(tsptr, p1, p2, &dens1, &dens2, ran1);
          if (!dens1 || !dens2)
            continue;

          if (initial) {
            tdens += dens1 + dens2;
            chfill1m(iniBook, p1, dens1, bookBins, 6);
            chfill1m(iniBook, p1, dens2, bookBins, 6);
          }
          /* This is very important. Change from slop to MeV */
          for(j=3; j<5; j++) {
            p1[j] *= tsptr->pCentral_mev;
            p2[j] *= tsptr->pCentral_mev;
          }
          p1[5] = (p1[5]+1)*tsptr->pCentral_mev;
          p2[5] = (p2[5]+1)*tsptr->pCentral_mev;
            
          bunch2cm(p1,p2,qa,beta,&gamma);
          
          theta = (ran1[9]*0.9999+0.00005)*PI;
          phi = ran1[10]*PI;
	  
          temp = dens1*dens2*sin(theta);
          eulertrans(qa,theta,phi,qb,&qabs);
          cm2bunch(p1,p2,qb,beta,&gamma);
          p1[5] -= tsptr->pCentral_mev;
          p2[5] -= tsptr->pCentral_mev;
	  
          if(fabs(p1[5])>tsSpec->delta_mev || fabs(p2[5])>tsSpec->delta_mev) {
            beta0=qabs/sqrt(qabs*qabs+me_mev*me_mev);
            cross = moeller(beta0,theta);
            temp *= cross*beta0/gamma/gamma;
	    
            if(fabs(p1[5])>tsSpec->delta_mev) {
              tsptr->totalWeight += temp;
              p1[3] /= tsptr->pCentral_mev;
              p1[4] /= tsptr->pCentral_mev;
              p1[5] /= tsptr->pCentral_mev;
              if (distribution) {
                chfill1m(outBook, p1, temp, bookBins, 6);
              }
              tsptr->simuCount++;
              
              beam0->particle[i][0] = p1[0];
              beam0->particle[i][1] = p1[3];
              beam0->particle[i][2] = p1[1];
              beam0->particle[i][3] = p1[4];
              beam0->particle[i][4] = p1[2];
              beam0->particle[i][5] = p1[5];
              beam0->particle[i][6] = i+1;
              weight[i] = temp;
              i++;
            }
              
            if(i>=n_simulated)
              break;
	    
            if(fabs(p2[5])>tsSpec->delta_mev) {
              tsptr->totalWeight += temp;
              p2[3] /= tsptr->pCentral_mev;
              p2[4] /= tsptr->pCentral_mev;
              p2[5] /= tsptr->pCentral_mev;
              if (distribution) {
                chfill1m(outBook, p2, temp, bookBins, 6);
              }
              tsptr->simuCount++;
                
              beam0->particle[i][0] = p2[0];
              beam0->particle[i][1] = p2[3];
              beam0->particle[i][2] = p2[1];
              beam0->particle[i][3] = p2[4];
              beam0->particle[i][4] = p2[2];
              beam0->particle[i][5] = p2[5];
              beam0->particle[i][6] = i+1;
              weight[i] = temp;
              i++;
            }
          }
        }
        if (total_event/tsptr->simuCount > 20) {
          if (distribution_cutoff[0]<5 || distribution_cutoff[1]<5 ) 
            fprintf(stdout, "waring: Scattering rate is low, please use 5 sigma beam for better simulation.\n");
          else
            fprintf(stdout, "waring: Scattering rate is very low, please ignore the rate from Monte Carlo simulation. Use Piwinski's rate only\n"); 
        }
        tsptr->factor = tsptr->factor / (double)(total_event);
        tsptr->s_rate = tsptr->totalWeight * tsptr->factor;
        tsptr->i_rate = tsptr->s_rate*tsSpec->ignoredPortion;
	
        /* Pick tracking particles from the simulated scattered particles */
        index = (long*)malloc(sizeof(long)*tsptr->simuCount);
        for (i=0; i<tsptr->simuCount; i++) index[i]=i;
        if (tsSpec->ignoredPortion <=1e-6) {
          iTotal = tsptr->simuCount; 
          wTotal = tsptr->totalWeight;
        } else {
          iTotal = 0;
          wTotal =0.;
          weight_limit = tsptr->totalWeight*(1-tsSpec->ignoredPortion);
          weight_ave = tsptr->totalWeight/tsptr->simuCount;
          pickPart(weight, index, 0, tsptr->simuCount,  
                   &iTotal, &wTotal, weight_limit, weight_ave);
        }
        
        beam->particle = (double**)czarray_2d(sizeof(double), iTotal, 7);
        beam->original = (double**)czarray_2d(sizeof(double), iTotal, 7);
        beam->accepted = NULL;
        beam->n_original = beam->n_to_track = beam->n_particle = iTotal;
        beam->n_accepted = beam->n_saved = 0;
        beam->p0_original = beam->p0 = tsptr->betagamma;
        beam->bunchFrequency = 0.;
        beam->lostOnPass = tmalloc(sizeof(*(beam->lostOnPass))*beam->n_to_track);
	  
        for (i=0; i<iTotal; i++) {
          beam->original[i][0] = beam->particle[i][0] = beam0->particle[index[i]][0];
          beam->original[i][1] = beam->particle[i][1] = beam0->particle[index[i]][1];
          beam->original[i][2] = beam->particle[i][2] = beam0->particle[index[i]][2];
          beam->original[i][3] = beam->particle[i][3] = beam0->particle[index[i]][3];
          beam->original[i][4] = beam->particle[i][4] = beam0->particle[index[i]][4];
          beam->original[i][5] = beam->particle[i][5] = beam0->particle[index[i]][5];
          beam->original[i][6] = beam->particle[i][6] = i+1;
          weight[i] *= tsptr->factor;
        }
        if (distIn)
          tsptr->total_scatter *= tsptr->s_rate/tsptr->p_rate;

        if (bunch)
          dump_scattered_particles(&SDDS_bunch, beam->particle, (long)iTotal,
                                   weight, tsptr);
        if (distribution) {
          /*
          print_hbook(x, y, s, xp, yp, dp, tsptr, tsptr->disFile, 
                      "Distribution of Scattered particles", 1, x->count*tsptr->factor/tsptr->s_rate);
          */
          free_hbook1m(outBook);
        }
        if (initial) {
          /*
          print_hbook(x0, y0, s0, xp0, yp0, dp0, tsptr, tsptr->iniFile, 
                      "Distribution of simulated particles", 1, x0->count/tdens);
          */
          free_hbook1m(iniBook);
        }
#if USE_MPI
      }
      MPI_Barrier(MPI_COMM_WORLD); 
#endif
      if (do_track) {
        n_left = do_tracking(beam, NULL, (long)iTotal, NULL, beamline, 
                             &beam->p0, NULL, NULL, NULL, NULL, run, 1,
                             0, 1, 0, NULL,
                             NULL, NULL, beam->lostOnPass, eptr);
#if USE_MPI
        MPI_Barrier(MPI_COMM_WORLD); 
        if (myid==0) {
#endif
          if (loss) {
            dump_scattered_loss_particles(&SDDS_loss, beam->particle+n_left, beam->original,  
                                          beam->lostOnPass+n_left, beam->n_to_track-n_left, weight, tsptr);
            if (!SDDS_Terminate(&SDDS_loss)) {
              SDDS_SetError("Problem terminating 'losses' file (finish_output)");
              SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
            }
          }
          if (output) {
            for (i=0; i< beam->n_to_track-n_left; i++) {
              j = (beam->particle+n_left)[i][6]-1;
              chfill1(lossDis, (beam->particle+n_left)[i][4], weight[j]*tsptr->total_scatter/tsptr->s_rate);
            }
            chprint1(lossDis, tsptr->outFile, "Beam loss distribution in particles/s/m", NULL,
                     NULL, 0, 0, verbosity, 0);
            free_hbook1(lossDis);
          }
#if USE_MPI
        }
#endif
      }
#if USE_MPI
      if (myid==0) {
#endif
        free_beamdata(beam);
        free_beamdata(beam0);
        free(weight);
#if USE_MPI
      }
#endif
    }
    eptr = eptr->succ; 
  }
  return;
}

void TouschekDistributionE(RUN *run, LINE_LIST *beamline)
{
  ELEMENT_LIST *eptr;
  TSCATTER *tsptr;
  static SDDS_TABLE SDDS_bunch, SDDS_loss;
  BEAM  Beam0, *beam0;
  book1 *lossDis;
  long i, n_left, eleCount=0, sTotal;
  double *weight, **fiducial_part;
 
  eptr = &(beamline->elem);
  beam0 = &Beam0;
  sTotal = (int)beamline->revolution_length;

  while (eptr) {
    if (eptr->type == T_TSCATTER) {
      eleCount++;
      if(eleCount < i_start) {
        eptr = eptr->succ; 
        continue;
      }
      if(eleCount > i_end)
        break;
#if USE_MPI
      MPI_Barrier(MPI_COMM_WORLD); 
      if (myid==0) {
#endif
        weight = (double*)malloc(sizeof(double)*n_simulated);
        beam0->particle = (double**)czarray_2d(sizeof(double), n_simulated, 7);
        beam0->original = (double**)czarray_2d(sizeof(double), n_simulated, 7);
        beam0->accepted = NULL;
        beam0->n_original = beam0->n_to_track = beam0->n_particle = n_simulated;
        beam0->n_accepted = beam0->n_saved = 0;
        beam0->p0_original = beam0->p0 = tsptr->betagamma;
        beam0->bunchFrequency = 0.;
        beam0->lostOnPass = tmalloc(sizeof(*(beam0->lostOnPass))*beam0->n_to_track);
        if (output) {
          tsptr->outFile = compose_filename1(output, run->rootname, eleCount);
          lossDis = chbook1("los_distribution", "m", 0, sTotal, sTotal);
        }
        if (loss) {
          tsptr->losFile = compose_filename1(loss, run->rootname, eleCount);
          SDDS_BeamScatterLossSetup(&SDDS_loss, tsptr->losFile, SDDS_BINARY, 1, 
                                    "lost particle coordinates", run->runfile,
                                    run->lattice, "touschek_scatter");
        }
        if (bunch) {
          tsptr->bunFile = compose_filename1(bunch, run->rootname, eleCount);
          SDDS_BeamScatterSetup(&SDDS_bunch, tsptr->bunFile, SDDS_BINARY, 1, 
                                "scattered-beam phase space", run->runfile,
                                run->lattice, "touschek_scatter");
        }
        for (i=0; i<n_simulated; i++) {
          beam0->original[i][0] =beam0->particle[i][0] = 0;
          beam0->original[i][1] =beam0->particle[i][1] = 0;
          beam0->original[i][2] =beam0->particle[i][2] = 0;
          beam0->original[i][3] =beam0->particle[i][3] = 0;
          beam0->original[i][4] =beam0->particle[i][4] = 0;
          if (tsSpec->delta_mev >= tsptr->pCentral_mev)
            beam0->original[i][5] = beam0->particle[i][5] = ((double)i/(n_simulated-1)-0.5)*2.;
          else
            beam0->original[i][5] = beam0->particle[i][5] = ((double)i/(n_simulated-1)-0.5)*2.*tsSpec->delta_mev/tsptr->pCentral_mev;
          beam0->original[i][6] =beam0->particle[i][6] = i+1;
          weight[i] = 1.;
        }
        if (bunch)
          dump_scattered_particles(&SDDS_bunch, beam0->particle, (long)n_simulated,
                                   weight, tsptr);
#if USE_MPI
      }
      MPI_Barrier(MPI_COMM_WORLD); 
#endif
      if (do_track) {
        delete_phase_references();
        reset_special_elements(beamline, 1);
        reset_driftCSR();
        fiducial_part=(double**)czarray_2d(sizeof(double), 1, 7); 
        for (fiducial_part[0][6]=1, i=0; i<6;i++)
          fiducial_part[0][i] = 0.0;
        n_left = do_tracking(NULL, fiducial_part, 1, NULL, beamline, 
                             &beam0->p0, NULL, NULL, NULL, NULL, run, 1,
                             513, 1, 0, NULL,
                             NULL, NULL, NULL, eptr);
        beam0->p0 = beam0->p0_original;
        n_left = do_tracking(beam0, NULL, (long)n_simulated, NULL, beamline, 
                             &beam0->p0, NULL, NULL, NULL, NULL, run, 2,
                             0, 1, 0, NULL,
                             NULL, NULL, beam0->lostOnPass, eptr);
        free_czarray_2d((void**)fiducial_part, 1, 7);
#if USE_MPI
        MPI_Barrier(MPI_COMM_WORLD); 
        if (myid==0) {
#endif
          if (loss) {
            dump_scattered_loss_particles(&SDDS_loss, beam0->particle+n_left, beam0->original,  
                                          beam0->lostOnPass+n_left, beam0->n_to_track-n_left, weight, tsptr);
            if (!SDDS_Terminate(&SDDS_loss)) {
              SDDS_SetError("Problem terminating 'losses' file (finish_output)");
              SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
            }
          }
          if (output) {
            for (i=0; i< beam0->n_to_track-n_left; i++) {
              chfill1(lossDis, (beam0->particle+n_left)[i][4], 1.);
            }
            chprint1(lossDis, tsptr->outFile, "Beam loss distribution in particles/s/m", NULL,
                     NULL, 0, 0, verbosity, 0);
            free_hbook1(lossDis);
          }
#if USE_MPI
        }
#endif
      }
#if USE_MPI
      if (myid==0) {
#endif
        free_beamdata(beam0);
        free(weight);
#if USE_MPI
      }
#endif
    }
    eptr = eptr->succ; 
  }
  return;
}

/* Initialize beam parameter at each Scatter element */
TSCATTER *initTSCATTER (ELEMENT_LIST *eptr)
{
  TSCATTER *tsptr;
  double temp;

  tsptr = ((TSCATTER*)eptr->p_elem);

  tsptr->betagamma = eptr->Pref_output;
  tsptr->gamma = sqrt(sqr(tsptr->betagamma)+1.);
  tsptr->pCentral_mev = tsptr->gamma * me_mev;
  tsptr->s_rate = tsptr->i_rate = tsptr->totalWeight = 0;
  tsptr->simuCount = 0;
  tsptr->name = eptr->name;
  tsptr->s = eptr->end_pos;

  tsptr->twiss[0][0] = eptr->twiss->alphax;
  tsptr->twiss[0][1] = eptr->twiss->betax;
  tsptr->twiss[0][2] = (1.+ sqr(eptr->twiss->alphax))/eptr->twiss->betax;
  tsptr->twiss[1][0] = eptr->twiss->alphay;
  tsptr->twiss[1][1] = eptr->twiss->betay;
  tsptr->twiss[1][2] = (1.+ sqr(eptr->twiss->alphay))/eptr->twiss->betay;
  tsptr->twiss[2][0] = 0.0;
  tsptr->twiss[2][1] = tsSpec->sigz*eptr->Pref_output/tsSpec->delta_p0;
  tsptr->twiss[2][2] = 1./tsptr->twiss[2][1];
  tsptr->disp[0][0] = eptr->twiss->etax;
  tsptr->disp[0][1] = eptr->twiss->etapx;
  tsptr->disp[1][0] = eptr->twiss->etay;
  tsptr->disp[1][1] = eptr->twiss->etapy;

  tsptr->sigx = sqrt(tsSpec->emitN[0]/tsptr->gamma*tsptr->twiss[0][1]);
  tsptr->sigy = sqrt(tsSpec->emitN[1]/tsptr->gamma*tsptr->twiss[1][1]);
  tsptr->sigz = tsSpec->sigz;
  tsptr->sigxyz = tsptr->sigx * tsptr->sigy * tsptr->sigz;

  temp = sqr(tsSpec->charge/e_mks)*sqr(PI)*sqr(re_mks)*c_mks/4.;
  if (distIn) {
    tsptr->factor = temp*(tranMax[0]-tranMin[0])*sqr(tranMax[1]-tranMin[1])
      *(tranMax[2]-tranMin[2])*sqr(tranMax[3]-tranMin[3])
      *(longMax[0]-longMin[0])*sqr(longMax[1]-longMin[1]);
  } else {
    tsptr->factor = temp*pow(tsSpec->range[0],3.0)*pow(tsSpec->range[1],3.0)*pow(tsSpec->range[2],3.0)
      /pow(2*PI, 6.0)/tsptr->sigxyz;
  }

  return (tsptr);
}

/* Initialize Touschek Scattering Calculation */
void init_TSPEC (RUN *run)
{
  if (!(tsSpec = SDDS_Malloc(sizeof(*tsSpec))))
    bomb("memory allocation failure at setup Touscheck scatter", NULL);                

  tsSpec->nbins = nbins;
  tsSpec->charge = charge;
  tsSpec->frequency = frequency;
  tsSpec->ignoredPortion = ignored_portion;
  tsSpec->emitN[0] = emit_nx;
  tsSpec->emitN[1] = emit_ny;
  tsSpec->emitN[2] = sigma_s * sigma_dp * sqrt(sqr(run->p_central)+1.);
  tsSpec->range[0] = 2 * distribution_cutoff[0];
  tsSpec->range[1] = 2 * distribution_cutoff[1];
  tsSpec->range[2] = 2 * distribution_cutoff[2];
  tsSpec->sigz = sigma_s;
  tsSpec->delta_p0 = sigma_dp * run->p_central;
  tsSpec->delta_mev = delta_mev;

  return;
}

/* Initialize distribution function input */
int init_distInput ()
{
  long i, nrow;
  double sum;

  if (verbosity)
    fprintf( stdout, "Opening \"%s\" for checking presence of parameters.\n", tranDist);
  if (!SDDS_InitializeInput(&tranPage, tranDist))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  SDDS_ReadPage(&tranPage);
  if (!SDDS_GetParameters(&tranPage,
                          "Variable00Min", &tranMin[0],
                          "Variable01Min", &tranMin[1],
                          "Variable02Min", &tranMin[2],
                          "Variable03Min", &tranMin[3],
                          "Variable00Max", &tranMax[0],
                          "Variable01Max", &tranMax[1],
                          "Variable02Max", &tranMax[2],
                          "Variable03Max", &tranMax[3],
                          "Variable00Interval", &tranInterval[0],
                          "Variable01Interval", &tranInterval[1],
                          "Variable02Interval", &tranInterval[2],
                          "Variable03Interval", &tranInterval[3],
                          "Variable00Dimension", &tranDimen[0],
                          "Variable01Dimension", &tranDimen[1],
                          "Variable02Dimension", &tranDimen[2],
                          "Variable03Dimension", &tranDimen[3],
                          NULL) )
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  for (i=0; i<4; i++)
    tranDelta[i] = tranMax[i] - tranMin[i];
  nrow = SDDS_RowCount(&tranPage);
  tweight = (double*)malloc(sizeof(double)*nrow);
  tweight = SDDS_GetColumnInDoubles(&tranPage, "weight");
  for (sum=0., i=0; i<nrow; i++)
     sum += tweight[i];
  for (i=0; i<nrow; i++)
    tweight[i] /= sum;

  if (verbosity)
    fprintf( stdout, "Opening \"%s\" for checking presence of parameters.\n", longDist);
  if (!SDDS_InitializeInput(&longPage, longDist))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  SDDS_ReadPage(&longPage);
  if (!SDDS_GetParameters(&longPage,
                          "Variable00Min", &longMin[0],
                          "Variable01Min", &longMin[1],
                          "Variable00Max", &longMax[0],
                          "Variable01Max", &longMax[1],
                          "Variable00Interval", &longInterval[0],
                          "Variable01Interval", &longInterval[1],
                          "Variable00Dimension", &longDimen[0],
                          "Variable01Dimension", &longDimen[1],
                          NULL) )
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  longMin[0] *= c_mks;
  longMax[0] *= c_mks;
  longInterval[0] *= c_mks;
  for (i=0; i<2; i++)
    longDelta[i] = longMax[i] - longMin[i];

  nrow = SDDS_RowCount(&longPage);
  lweight = (double*)malloc(sizeof(double)*nrow);
  lweight = SDDS_GetColumnInDoubles(&longPage, "weight");
  for (sum=0., i=0; i<nrow; i++)
     sum += lweight[i];
  for (i=0; i<nrow; i++)
    lweight[i] /= sum;

  if (!SDDS_Terminate(&tranPage) || !SDDS_Terminate(&longPage))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);

  return (0);
}

/************************************************************************************\
 * Modified from S. Khan's code.                                                    *
 *      select two electrons (labelled a and b) for Touschek scattering             *
 * p1[x,y,s,xp,yp,dp/p]                                                             *
\************************************************************************************/

void selectPart0(TSCATTER *tsptr, double *p1, double *p2, 
                double *dens1, double *dens2, double *ran1)
{
  int i;
  double U[3], V1[3], V2[3], densa[3], densb[3];

  /* Select random particle's parameter in normal phase space */
  for (i=0; i<3; i++) {
    U[i]  = (ran1[i]-0.5) * tsSpec->range[i]*sqrt(tsSpec->emitN[i]/tsptr->gamma);
    V1[i] = (ran1[i+3]-0.5) * tsSpec->range[i]*sqrt(tsSpec->emitN[i]/tsptr->gamma);
    V2[i] = (ran1[i+6]-0.5) * tsSpec->range[i]*sqrt(tsSpec->emitN[i]/tsptr->gamma);
    densa[i] = exp(-0.5*(U[i]*U[i]+V1[i]*V1[i])/tsSpec->emitN[i]*tsptr->gamma);
    if (i==2)
      densb[i] = exp(-0.5*(U[i]*U[i]+V2[i]*V2[i])/tsSpec->emitN[i]*tsptr->gamma);
  }
  /* Change particle's parameter from normal phase space to real phase space */
  for (i=0; i<3; i++) {
    p1[i] = p2[i] = sqrt(tsptr->twiss[i][1])*U[i];
    p1[i+3] = (V1[i] - tsptr->twiss[i][0]*U[i])/sqrt(tsptr->twiss[i][1]);
    p2[i+3] = (V2[i] - tsptr->twiss[i][0]*U[i])/sqrt(tsptr->twiss[i][1]);
  }
  /* Dispersion correction */
  p1[0] = p1[0] + p1[5]*tsptr->disp[0][0];
  p1[1] = p1[1] + p1[5]*tsptr->disp[1][0];
  p1[3] = p1[3] + p1[5]*tsptr->disp[0][1];
  p1[4] = p1[4] + p1[5]*tsptr->disp[1][1];

  p2[0] = p1[0] - p2[5]*tsptr->disp[0][0];
  p2[1] = p1[1] - p2[5]*tsptr->disp[1][0];
  U[0] = p2[0]/sqrt(tsptr->twiss[0][1]);
  U[1] = p2[1]/sqrt(tsptr->twiss[1][1]);
  p2[3] = (V2[0] - tsptr->twiss[0][0]*U[0])/sqrt(tsptr->twiss[0][1]);
  p2[4] = (V2[1] - tsptr->twiss[1][0]*U[1])/sqrt(tsptr->twiss[1][1]);
  densb[0] = exp(-0.5*(U[0]*U[0]+V2[0]*V2[0])/tsSpec->emitN[0]*tsptr->gamma);
  densb[1] = exp(-0.5*(U[1]*U[1]+V2[1]*V2[1])/tsSpec->emitN[1]*tsptr->gamma);

  p2[0] = p1[0];
  p2[1] = p1[1];
  p2[3] = p2[3] + p2[5]*tsptr->disp[0][1];
  p2[4] = p2[4] + p2[5]*tsptr->disp[1][1];

  *dens1 = densa[0] * densa[1] * densa[2]; 
  *dens2 = densb[0] * densb[1] * densb[2]; 

  return;
}

void selectPart1(TSCATTER *tsptr, double *p1, double *p2, 
                double *dens1, double *dens2, double *ran1)
{
  /* select {z,deltaP} for p1 and p2 (use 3 random number) */
  p2[2] = p1[2] = ran1[0]*longDelta[0]+longMin[0];
  p1[5] = ran1[1]*longDelta[1]+longMin[1];
  p2[5] = ran1[2]*longDelta[1]+longMin[1];
  
  /* select {x,xp} {y,yp} for p1 and p2 (use 6 random number) */
  p1[0] = ran1[3]*tranDelta[0]+tranMin[0];
  p1[1] = ran1[4]*tranDelta[2]+tranMin[2];
  p1[3] = ran1[5]*tranDelta[1]+tranMin[1];
  p2[3] = ran1[6]*tranDelta[1]+tranMin[1];
  p1[4] = ran1[7]*tranDelta[3]+tranMin[3];
  p2[4] = ran1[8]*tranDelta[3]+tranMin[3];
  if (!(*dens1 = findDens(p1)))
    return;

  /* Dispersion correction */
  p1[0] = p1[0] + p1[5]*tsptr->disp[0][0];
  p1[1] = p1[1] + p1[5]*tsptr->disp[1][0];
  p1[3] = p1[3] + p1[5]*tsptr->disp[0][1];
  p1[4] = p1[4] + p1[5]*tsptr->disp[1][1];

  p2[0] = p1[0] - p2[5]*tsptr->disp[0][0];
  p2[1] = p1[1] - p2[5]*tsptr->disp[1][0];
  if (!(*dens2 = findDens(p2)))
    return;

  p2[0] = p1[0];
  p2[1] = p1[1];
  p2[3] = p2[3] + p2[5]*tsptr->disp[0][1];
  p2[4] = p2[4] + p2[5]*tsptr->disp[1][1];

  return;
}

double findDens(double *p)
{
  double p0[6];
  double x0[16], x1[8], x2[4], x3[2], x4;
  double t0[4], t1[2], t2;
  long i, j, k, l, index1, index2;
  long Ix[4], It[2];
  double Cx[4][2], Ct[2][2];

  for (i=0; i<3; i++) {
    p0[i*2] = p[i];
    p0[i*2+1] = p[i+3];
  }

  /* do longitudinal first */
  for (i=0; i<2; i++) {
    Ct[i][0] = (p0[i+4]-longMin[i])/longInterval[i]+0.5;
    It[i] = (long)Ct[i][0];
    Ct[i][1] = Ct[i][0] - It[i];
    Ct[i][0] = 1. - Ct[i][1];
    It[i]--;
  } 

  index1 = 0;
  for (i=0; i<2; i++) {
    for (j=0; j<2; j++) {
      if ((It[0]+i)<0 || (It[0]+i)==longDimen[0] ||
          (It[1]+j)<0 || (It[1]+j)==longDimen[1])
        t0[index1]=0.;
      else {
        index2 = (It[1]+j)+(It[0]+i)*longDimen[1];
        t0[index1] = lweight[index2];
      }
      index1++;
    }
  }

  if (!(t0[0] || t0[1] || t0[2] || t0[3]))
    return (0);

  for (i=0; i<2; i++) {
    t1[i] = (Ct[1][0]*t0[i*2] + Ct[1][1]*t0[i*2+1])/longInterval[1];
  }
  t2 = (Ct[0][0]*t1[0] + Ct[0][1]*t1[1])/longInterval[0];

  /* do transverse part */
  for (i=0; i<4; i++) {
    Cx[i][0] = (p0[i]-tranMin[i])/tranInterval[i]+0.5;
    Ix[i] = (long)Cx[i][0];
    Cx[i][1] = Cx[i][0] - Ix[i];
    Cx[i][0] = 1. - Cx[i][1];
    Ix[i]--;
  } 

  index1 = 0;
  for (i=0; i<2; i++) {
    for (j=0; j<2; j++) {
      for (k=0; k<2; k++) {
        for (l=0; l<2; l++) {
          if ((Ix[0]+i)<0 || (Ix[0]+i)==tranDimen[0] ||
              (Ix[1]+j)<0 || (Ix[1]+j)==tranDimen[1] ||
              (Ix[2]+k)<0 || (Ix[2]+k)==tranDimen[2] ||
              (Ix[3]+l)<0 || (Ix[3]+l)==tranDimen[3])
            x0[index1] = 0.;
          else {
            index2 = (Ix[3]+l)+(Ix[2]+k)*tranDimen[3]
              +(Ix[1]+j)*tranDimen[3]*tranDimen[2]
              +(Ix[0]+i)*tranDimen[3]*tranDimen[2]*tranDimen[1];
            x0[index1] = tweight[index2];
          }
          index1++;
        }
      }
    }
  }

  for (i=0; i<8; i++) {
    x1[i] = (Cx[3][0]*x0[i*2] + Cx[3][1]*x0[i*2+1])/tranInterval[3];
  }
  for (i=0; i<4; i++) {
    x2[i] = (Cx[2][0]*x1[i*2] + Cx[2][1]*x1[i*2+1])/tranInterval[2];
  }
  for (i=0; i<2; i++) {
    x3[i] = (Cx[1][0]*x2[i*2] + Cx[1][1]*x2[i*2+1])/tranInterval[1];
  }
  x4 = (Cx[0][0]*x3[0] + Cx[0][1]*x3[1])/tranInterval[0]; 

  return (x4*t2);
}

/* This function transfer particle's momentum from lab system to c.o.m system */
/* p1[3]=p1x*c, p1[4]=p1y*c, p1[5]=p1z*c in MeV */
/* p2[3]=p2x*c, p2[4]=p2y*c, p2[5]=p1z*c in MeV */
void bunch2cm(double *p1, double *p2, double *q, double *beta, double *gamma)
{
  double pp1, pp2, e1, e2, ee;
  int i;
  double bb, betap1, factor;

  pp1=0.0;
  pp2=0.0;
  for(i=3; i<6; i++) {
    pp1 = pp1 + sqr(p1[i]);
    pp2 = pp2 + sqr(p2[i]);
  }
  e1=sqrt(me_mev*me_mev+pp1);
  e2=sqrt(me_mev*me_mev+pp2);
  ee=e1+e2;

  betap1=0.0;
  bb=0.0;
  for(i=0; i<3; i++) {
    beta[i]=(p1[i+3]+p2[i+3])/ee;
    betap1=betap1+beta[i]*p1[i+3];
    bb=bb+beta[i]*beta[i];
  }

  *gamma = 1./sqrt(1.-bb);
  factor = ((*gamma)-1.)*betap1/bb;

  for(i=0; i<3; i++) {
    q[i]=p1[i+3]+factor*beta[i]-(*gamma)*e1*beta[i];
  }

  return;
}

/* Rotate scattered p in c.o.m system */
void eulertrans(double *v0, double theta, double phi, double *v1, double *v)
{
  double th, ph, s1, s2, c1, c2;
  double x0, y0, z0; 

  *v=sqrt(v0[0]*v0[0]+v0[1]*v0[1]+v0[2]*v0[2]);
  th=acos(v0[2]/(*v));
  ph=atan2(v0[1],v0[0]);

  s1=sin(th);
  s2=sin(ph);
  c1=cos(th);
  c2=cos(ph);

  x0=cos(theta);
  y0=sin(theta)*cos(phi);
  z0=sin(theta)*sin(phi);

  v1[0] = (*v) * (s1*c2*x0 - s2*y0 - c1*c2*z0);
  v1[1] = (*v) * (s1*s2*x0 + c2*y0 - c1*s2*z0);
  v1[2] = (*v) * (   c1*x0         +    s1*z0);

  return;
}

void cm2bunch(double *p1, double *p2, double *q, double *beta, double *gamma)
{
  int i;
  double pq, e, betaq, bb, factor;

  pq=0.0;
  for(i=0; i<3; i++) {
    pq = pq + q[i]*q[i];
  }

  e=sqrt(me_mev*me_mev+pq);

  betaq=0.0;
  bb=0.0;
  for(i=0; i<3; i++) {
    betaq = betaq + beta[i]*q[i];
    bb = bb + beta[i]*beta[i];
  }

  factor=((*gamma)-1)*betaq/bb;
  for(i=0; i<3; i++) {
    p1[i+3]= q[i] + (*gamma)*beta[i]*e + factor*beta[i];
    p2[i+3]=-q[i] + (*gamma)*beta[i]*e - factor*beta[i];
  }

  return;
}

double moeller(double beta0, double theta)
{
  double cross; 
  double beta2, st2;

  beta2=beta0*beta0;
  st2=sqr(sin(theta));

  cross = (1.-beta2)*(sqr(1.+1./beta2)*(4./st2/st2-3./st2)+1.+4./st2);
 
  return cross;
}

void pickPart(double *weight, long *index, long start, long end, 
              long *iTotal, double *wTotal, double weight_limit, double weight_ave)
{
  long i, i1, i2, N;
  double w1, w2;
  long *index1, *index2;
  double *weight1, *weight2;

  i1=i2=0;
  w1=w2=0.;
  N = end-start;
  if(N<3) return;  /* scattered particles normally appear in pair */
  index2 = (long*)malloc(sizeof(long)*N);
  weight2 = (double*)malloc(sizeof(double)*N);
  index1 = (long*)malloc(sizeof(long)*N);
  weight1 = (double*)malloc(sizeof(double)*N);
  
  for (i=start; i<end; i++) {
    if (weight[i] > weight_ave) {
      weight2[i2] = weight[i];
      index2[i2++] = index[i];
      w2 += weight[i];
    } else {
      weight1[i1] = weight[i];
      index1[i1++] = index[i];
      w1 += weight[i];
    }
  }
  if ((w2+ (*wTotal)) > weight_limit) {
    weight_ave = w2/(double)i2;
    for (i=0; i<i2; i++) {
      index[start+i]=index2[i];
      weight[start+i]=weight2[i];
    }
    free(weight1);
    free(index1);
    free(weight2);
    free(index2);
    pickPart(weight, index, start, start+i2,
             iTotal, wTotal, weight_limit, weight_ave);
    return;
  }

  *iTotal += i2;
  *wTotal += w2;
  weight_ave = w1/(double)i1;
  for (i=0; i<i2; i++) {
    index[start+i]=index2[i];
    weight[start+i]=weight2[i];
  }
  for (i=0; i<i1; i++) {
    index[start+i2+i] = index1[i];
    weight[start+i2+i] = weight1[i];
  }
  free(weight1);
  free(index1);
  free(weight2);
  free(index2);
  pickPart(weight, index, i2+start, end,
           iTotal, wTotal, weight_limit, weight_ave);
  return;
}
