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
#include "chbook.h"

#define DEBUG 0
#define NDIV 10000;

static TSCATTER_SPEC *tsSpec = NULL;
void init_TSPEC ();

int TouschekRate(LINE_LIST *beamline);
void FIntegral(double tm, double B1, double B2, double *F, char *name);
double Fvalue (double t, double tm, double b1, double b2);

void TouschekDistirbution(RUN *run, LINE_LIST *beamline);
TSCATTER *initTSCATTER (ELEMENT_LIST *eptr);

void selectPart(long iseed, TSCATTER *tsptr, double *p1, double *p2, 
                double *dens1, double *dens2);
void bunch2cm(double *p1, double *p2, double *q, double *beta);
void eulertrans(double *v0, double theta, double phi, double *v1, double *v);
void cm2bunch(double *p1, double *p2, double *q, double *beta);
double moeller(double vc, double theta, double phi);

char *compose_filename1(char *template, char *root_name, long index);
int str_inp(const char *string1, const char *string2);
void fill_hbook (book1 *x, book1 *y, book1 *s, book1 *xp, book1 *yp, 
                 book1 *dp, double *p, double weight);
void print_hbook (book1 *x, book1 *y, book1 *s, book1 *xp, book1 *yp, book1 *dp,
                  TSCATTER *tsptr, char *filename, char *description, int verbosity);

void setupTouschekEffect(NAMELIST_TEXT *nltext, RUN *run, LINE_LIST *beamline) 
{
  long i;
  ELEMENT_LIST *eptr;
  int flag;

  /* Check if there is TScatter element along beamline. Then calculate Twiss. */
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

  run_twiss_output(run, beamline, NULL, -1);
 
  /* process the namelist text */
  set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
  set_print_namelist_flags(0);
  process_namelist(&touschek_scatter, nltext);
  print_namelist(stdout, &touschek_scatter);

  if (!charge)    
    bomb("charge has to be given", NULL);
  charge /= 1e9;
  if (frequency<1)
    bomb("frequency has to >=1", NULL);
  if (!delta)
    bomb("energy aperture has to be given", NULL);

  if (!p_central_mev)
    bomb("beam energy has to be given", NULL);
  if (!emittance[0] || !emittance[1])
    bomb("bunch emittance has to be given", NULL);
  emittance[0] /= 1.e9;
  emittance[1] /= 1.e9;
  if (!sigma_dp)
    bomb("energy spread has to be given", NULL);
  if (!sigma_s)
    bomb("bunch length has to be given", NULL);
  sigma_s /= 1.e3;

  init_TSPEC ();

  /* Calculate integrated scattering rate at TScatter. */
  if (TouschekRate(beamline))
    bomb("Touschek scattering rate calculation error", NULL);
  /* Generate scatter particles at TScatter. 
     And track through scattered particles down to beamline. */
  TouschekDistirbution(run, beamline);

  return;
}

int TouschekRate(LINE_LIST *beamline)
{
  double NP;
  double tm, B1, B2, F, rate, IntR;
  ELEMENT_LIST *eptr;

  double pCentral, gamma; 
  double sp2, sp4, beta2, betagamma2, sh2;
  double sx2, sxb2, dx2, dx_2;
  double sy2, syb2, dy2, dy_2;
  double a0, c0, c1, c2, c3;
  TWISS *twiss0;

  NP = charge/e_mks;
  gamma = p_central_mev/me_mev;
  pCentral = sqrt(sqr(gamma)-1);

  sp2 = sqr(sigma_dp);
  sp4 = sqr(sp2);
  beta2 = sqr(pCentral)/sqr(gamma);
  betagamma2 = sqr(pCentral);
  a0 = sqr(re_mks)*NP/(4*sqrt(PI)*sqr(gamma)*sigma_s);
  tm = beta2*delta*delta;

  IntR = 0.;
  eptr = &(beamline->elem);
  while (eptr) {
    if (eptr->type == T_TSCATTER) {
      ((TSCATTER*)eptr->p_elem)->IntR = IntR;
      ((TSCATTER*)eptr->p_elem)->p_rate = rate*c_mks;
      IntR = 0.;
    }
    if(!(entity_description[eptr->type].flags&HAS_LENGTH) ||
       (eptr->pred && !(((DRIFT*)eptr->p_elem)->length)) ||
       !eptr->succ) {
      eptr = eptr->succ;
      continue;
    }

    twiss0 = eptr->twiss;
    sxb2 = twiss0->betax*emittance[0];
    dx2  = sqr(twiss0->etax);
    sx2  = sxb2 + dx2*sp2;
    dx_2 = sqr(twiss0->alphax * twiss0->etax + twiss0->betax * twiss0->etapx);

    syb2 = twiss0->betay*emittance[1];
    dy2  = sqr(twiss0->etay);
    sy2  = syb2 + dy2*sp2;
    dy_2 = sqr(twiss0->alphay * twiss0->etay + twiss0->betay * twiss0->etapy);

    sh2 = 1/(1/sp2+(dx2+dx_2)/sxb2+(dy2+dy_2)/syb2);
    c0 = sqrt(sh2)/(sigma_dp*betagamma2*emittance[0]*emittance[1]);
    c1 = sqr(twiss0->betax)/(2*betagamma2*sxb2);
    c2 = sqr(twiss0->betay)/(2*betagamma2*syb2);
    c3 = sx2*sy2-sp4*dx2*dy2;
    B1 = c1*(1-sh2*dx_2/sxb2)+c2*(1-sh2*dy_2/syb2);
    B2 = sqr(B1)-sqr(c0)*c3;   	  
    if (B2<0) {
      fprintf(stdout, "B2^2<0 at \"%s\" occurence %ld", eptr->name, eptr->occurence);
      exit(1);
    }
    B2=sqrt(B2);   	  

    FIntegral(tm, B1, B2, &F, eptr->name);

    rate = a0*c0*F*NP;
    IntR += rate * frequency * ((DRIFT*)eptr->p_elem)->length;

    eptr = eptr->succ; 
  }

  return(0);
}

void FIntegral(double tm, double b1, double b2, double *F, char *name) 
{
  double f0, f1, sum;
  double HPI, step;
  long converge=0;
  double t1, k0, k1, km;
  
  HPI = PI/2;
  step = HPI/NDIV;
  k0 = km = atan(sqrt(tm));
  f0 = Fvalue(tm, tm, b1, b2);
  sum = 0;
  f1 = 0;
  
  while (!converge) {
    k1 = k0 + step;
    t1 = sqr(tan(k1));
    if (isnan(t1)) {
      fprintf(stdout, "Integration failed at Element %s", name);
      break;
    }
    
    f1 = Fvalue(t1, tm, b1, b2);

    if (f1>0 && f1*(HPI-k1)<1e-3*sum)
      converge =1;
    
    sum +=(f1+f0)/2*step;
    k0 = k1;
    f0 = f1;    
  }

  *F = sum;

  return;
}

double Fvalue (double t, double tm, double b1, double b2)
{
  double c0, c1, c2, result;

  c0 = (sqr(2*t+1)*(t/tm/(1+t)-1)/t+t-sqrt(t*tm*(1+t))-(2+1/(2*t))*log(t/tm/(1+t)))*sqrt(1+t);
  c1 = exp(-b1*t);
  c2 = dbesi0(b2*t);
  result = c0 * c1 * c2;
  /* If overflow/underflow use approximate equation for modified bessel function. */
  if (isnan(result) || result>MAXFLOAT) {
    result=c0*exp(b2*t-b1*t)/sqrt(PIx2*b2*t);
  } 
  return result;
}

char *compose_filename1(char *template, char *root_name, long index)
{
  char *root, *ext, *filename, *inx;
  long i, divp;

  divp = str_inp(template, ".");
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

  ext =  tmalloc(sizeof(char)*(strlen(template)-divp));
  for(i=0; i<strlen(template)-divp; i++) {
    *(ext+i) = *(template+i+divp);
  }
  *(ext+i) = '\0';

  inx = tmalloc(sizeof(char)*4);
  sprintf(inx, "%03ld", index);

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

void fill_hbook (book1 *x, book1 *y, book1 *s, book1 *xp, book1 *yp, 
                 book1 *dp, double *p, double weight)
{
  chfill1(x, p[0], weight);
  chfill1(y, p[1], weight);
  chfill1(s, p[2], weight);
  chfill1(xp, p[3], weight);
  chfill1(yp, p[4], weight);
  chfill1(dp, p[5], weight);

  return;
}

void print_hbook (book1 *x, book1 *y, book1 *s, book1 *xp, book1 *yp, book1 *dp,
                  TSCATTER *tsptr, char *filename, char *description, int verbosity)
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
  }

  /* Open file for writting */
  if (verbosity)
    fprintf( stdout, "Opening \"%s\" for writing...\n", filename);

  if (!SDDS_InitializeOutput(&outPage, SDDS_BINARY, 1, 
                             description, description, filename))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);

  if (0>SDDS_DefineParameter(&outPage, "VariableName_1", NULL, NULL, 
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
                             NULL, NULL, SDDS_DOUBLE, NULL) ||
      0>SDDS_DefineParameter(&outPage, "Piwinski_IntR", NULL, NULL, 
                             NULL, NULL, SDDS_DOUBLE, NULL) ||
      0>SDDS_DefineParameter(&outPage, "Simulated_Rate", NULL, NULL, 
                             NULL, NULL, SDDS_DOUBLE, NULL) ||
      0>SDDS_DefineParameter(&outPage, "nbins", NULL, NULL, 
                             NULL, NULL, SDDS_LONG, NULL) ||
      0>SDDS_DefineParameter(&outPage, "Total_count", NULL, NULL, 
                             NULL, NULL, SDDS_LONG, NULL))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);

  if (0>SDDS_DefineColumn(&outPage, "value1", NULL, NULL, 
                          NULL, NULL, SDDS_DOUBLE, 0) ||
      0>SDDS_DefineColumn(&outPage, "weight_1", NULL, NULL, 
                          NULL, NULL, SDDS_DOUBLE, 0) ||
      0>SDDS_DefineColumn(&outPage, "value2", NULL, NULL, 
                          NULL, NULL, SDDS_DOUBLE, 0) ||
      0>SDDS_DefineColumn(&outPage, "weight_2", NULL, NULL, 
                          NULL, NULL, SDDS_DOUBLE, 0) ||
      0>SDDS_DefineColumn(&outPage, "value3", NULL, NULL, 
                          NULL, NULL, SDDS_DOUBLE, 0) ||
      0>SDDS_DefineColumn(&outPage, "weight_3", NULL, NULL, 
                          NULL, NULL, SDDS_DOUBLE, 0) ||
      0>SDDS_DefineColumn(&outPage, "value4", NULL, NULL, 
                          NULL, NULL, SDDS_DOUBLE, 0) ||
      0>SDDS_DefineColumn(&outPage, "weight_4", NULL, NULL, 
                          NULL, NULL, SDDS_DOUBLE, 0) ||
      0>SDDS_DefineColumn(&outPage, "value5", NULL, NULL, 
                          NULL, NULL, SDDS_DOUBLE, 0) ||
      0>SDDS_DefineColumn(&outPage, "weight_5", NULL, NULL, 
                          NULL, NULL, SDDS_DOUBLE, 0) ||
      0>SDDS_DefineColumn(&outPage, "value6", NULL, NULL, 
                          NULL, NULL, SDDS_DOUBLE, 0) ||
      0>SDDS_DefineColumn(&outPage, "weight_6", NULL, NULL, 
                          NULL, NULL, SDDS_DOUBLE, 0))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);

  if (!SDDS_WriteLayout(&outPage) )
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);

  /* Write to output file */
  if (0>SDDS_StartPage(&outPage, x->length) ||
      !SDDS_SetParameters(&outPage, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,
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
                          "Piwinski_IntR", tsptr->IntR,
                          "Simulated_Rate", tsptr->s_rate,
                          "nbins", x->length,
                          "Total_count", x->count, NULL) ||
      !SDDS_SetColumn(&outPage, SDDS_SET_BY_NAME, v1,        x->length,  "value1") ||
      !SDDS_SetColumn(&outPage, SDDS_SET_BY_NAME, x->value,  x->length,  "weight_1") ||
      !SDDS_SetColumn(&outPage, SDDS_SET_BY_NAME, v2,        y->length,  "value2") ||
      !SDDS_SetColumn(&outPage, SDDS_SET_BY_NAME, y->value,  y->length,  "weight_2") ||
      !SDDS_SetColumn(&outPage, SDDS_SET_BY_NAME, v3,        s->length,  "value3") ||
      !SDDS_SetColumn(&outPage, SDDS_SET_BY_NAME, s->value,  s->length,  "weight_3") ||
      !SDDS_SetColumn(&outPage, SDDS_SET_BY_NAME, v4,        xp->length, "value4") ||
      !SDDS_SetColumn(&outPage, SDDS_SET_BY_NAME, xp->value, xp->length, "weight_4") ||
      !SDDS_SetColumn(&outPage, SDDS_SET_BY_NAME, v5,        yp->length, "value5") ||
      !SDDS_SetColumn(&outPage, SDDS_SET_BY_NAME, yp->value, yp->length, "weight_5") ||
      !SDDS_SetColumn(&outPage, SDDS_SET_BY_NAME, v6,        dp->length, "value6") ||
      !SDDS_SetColumn(&outPage, SDDS_SET_BY_NAME, dp->value, dp->length, "weight_6") ||
      !SDDS_WritePage(&outPage))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);

  if (!SDDS_Terminate(&outPage))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);

  return;
}

void TouschekDistirbution(RUN *run, LINE_LIST *beamline)
{
  long i, j, n_left;
  ELEMENT_LIST *eptr;
  TSCATTER *tsptr;
  double p1[6], p2[6], dens1, dens2;
  double theta, phi, phi05, qa[3], qb[3], beta[3], qabs;
  double energy0, vc, cross, temp;
  static SDDS_TABLE SDDS_bunch, SDDS_loss;
  BEAM  beam1, *beam;
  double xrange[6];
  book1 *x, *y, *s, *xp, *yp, *dp;
  double *weight;

  eptr = &(beamline->elem);
  beam = &beam1;
  weight = (double*)malloc(sizeof(double)*n_particles_per_bunch);
  while (eptr) {
    if (eptr->type == T_TSCATTER) {
      beam->particle = (double**)czarray_2d(sizeof(double), n_particles_per_bunch, 7);
      beam->original = beam->accepted = NULL;
      beam->n_original = beam->n_to_track = beam->n_particle = n_particles_per_bunch;
      beam->n_accepted = beam->n_saved = 0;
      beam->p0_original = beam->p0 = tsSpec->pCentral;
      beam->bunchFrequency = 0.;
      beam->lostOnPass = tmalloc(sizeof(*(beam->lostOnPass))*beam->n_to_track);

      tsptr = initTSCATTER (eptr);
      tsptr->disFile = compose_filename1(distribution, run->rootname, eptr->occurence);
      tsptr->losFile = compose_filename1(loss, run->rootname, eptr->occurence);
      tsptr->bunFile = compose_filename1(bunch, run->rootname, eptr->occurence);
      SDDS_BeamScatterSetup(&SDDS_bunch, tsptr->bunFile, SDDS_BINARY, 1, 
                            "scattered-beam phase space", run->runfile,
                            run->lattice, "touschek_scatter");
      SDDS_BeamLossSetup(&SDDS_loss, tsptr->losFile, SDDS_BINARY, 1, 
                         "lost particle coordinates", run->runfile,
                         run->lattice, "touschek_scatter");

      xrange[0] = tsptr->range[0]/2.;
      xrange[1] = tsptr->range[1]/2.;
      xrange[2] = tsptr->range[2]/2.;
      xrange[3] = (tsptr->range[3]+fabs(tsptr->twiss[0][0]/tsptr->twiss[0][1])*tsptr->range[0])/2.;
      xrange[4] = (tsptr->range[4]+fabs(tsptr->twiss[1][0]/tsptr->twiss[1][1])*tsptr->range[1])/2.;
      xrange[5] = tsptr->range[5]/2.;
      xrange[0] = xrange[0]+fabs(tsSpec->range[2]*tsSpec->sigma_p*tsptr->disp[0][0]);
      xrange[3] = xrange[3]+fabs(tsSpec->range[2]*tsSpec->sigma_p*tsptr->disp[0][1]);
      x  = chbook1("x", -xrange[0], xrange[0], tsSpec->nbins);
      y  = chbook1("y", -xrange[1], xrange[1], tsSpec->nbins);
      s  = chbook1("s", -xrange[2], xrange[2], tsSpec->nbins);
      xp = chbook1("xp", -xrange[3], xrange[3], tsSpec->nbins);
      yp = chbook1("yp", -xrange[4], xrange[4], tsSpec->nbins);
      dp = chbook1("dp", -0.07, 0.07, tsSpec->nbins);

      i = 0;
      while(1) {
        if(i>=n_particles_per_bunch)
          break;

        selectPart(tsSpec->seed, tsptr, p1, p2, &dens1, &dens2);

        /* This is very important. Change from slop to MeV */
        for(j=3; j<6; j++) {
          p1[j] *= tsSpec->p0;
          p2[j] *= tsSpec->p0;
        }

        bunch2cm(p1,p2,qa,beta);
        theta = (random_1(tsSpec->seed)*0.996+0.002)*PI;
        temp = dens1*dens2*sin(theta);
        phi = random_1(tsSpec->seed)*PI;
        eulertrans(qa,theta,phi,qb,&qabs);
        cm2bunch(p1,p2,qb,beta);

        if(fabs(p1[5])>tsSpec->dp0 || fabs(p2[5])>tsSpec->dp0) {
          energy0 = sqrt(qabs*qabs+me_mev*me_mev);
          vc=qabs/energy0;
          phi05=phi+PIo2;
          cross = moeller(vc,theta,phi05);
          temp *= cross*vc;

          if(fabs(p1[5])>tsSpec->dp0) {
            tsptr->totalWeight += temp;
            p1[3] /= tsSpec->gamma;
            p1[4] /= tsSpec->gamma;
            p1[5] /= energy0;
            fill_hbook (x, y, s, xp, yp, dp, p1, temp);
            tsptr->simuCount++;

            if (temp > tsSpec->weight_limit) {
              tsptr->ignorWeight += temp;
              beam->particle[i][0] = p1[0];
              beam->particle[i][1] = p1[3];
              beam->particle[i][2] = p1[1];
              beam->particle[i][3] = p1[4];
              beam->particle[i][4] = p1[2]/tsSpec->gamma;
              beam->particle[i][5] = p1[5];
              weight[i] = temp;
              beam->particle[i][6] = ++i;
            }
          }

          if(i>=n_particles_per_bunch)
            break;

          if(fabs(p2[5])>tsSpec->dp0) {
            tsptr->totalWeight += temp;
            p2[3] /= tsSpec->gamma;
            p2[4] /= tsSpec->gamma;
            p2[5] /= energy0;
            fill_hbook (x, y, s, xp, yp, dp, p2, temp);
            tsptr->simuCount++;
           
            if (temp > tsSpec->weight_limit) {
              tsptr->ignorWeight += temp;
              beam->particle[i][0] = p2[0];
              beam->particle[i][1] = p2[3];
              beam->particle[i][2] = p2[1];
              beam->particle[i][3] = p2[4];
              beam->particle[i][4] = p2[2]/tsSpec->gamma;
              beam->particle[i][5] = p2[5];
              weight[i] = temp;
              beam->particle[i][6] = ++i;
            }
          }
        }
        if (tsptr->simuCount > (long)1e7) break;
      }

      tsptr->ignorWeight = tsptr->totalWeight - tsptr->ignorWeight;
      if (tsptr->ignorWeight/tsptr->totalWeight < 1e-2) {
        if (tsptr->simuCount < (long)1e6)
          fprintf(stderr, "weight_limit is too low, try to increase it\n");
        if (tsptr->simuCount > (long)1e7)
          fprintf(stderr, "you should use less particles to increase tracking speed\n");
      }
      if (tsptr->ignorWeight/tsptr->totalWeight > 1e-2) {
        if (tsptr->simuCount < (long)1e6)
          fprintf(stderr, "you should use more particles for accurate results\n");
        if (tsptr->simuCount > (long)1e7)
          fprintf(stderr, "weight_limit is too high, try to lower it\n");
      }

      tsptr->s_rate = tsptr->totalWeight * tsptr->factor / (double)(tsptr->simuCount);
      dump_scattered_particles(&SDDS_bunch, beam->particle, (long)n_particles_per_bunch, 
                               weight, tsptr, tsSpec);
      print_hbook(x, y, s, xp, yp, dp, tsptr, tsptr->disFile, 
                  "Distribution of Scattered particles", 1);
      n_left = do_tracking(beam, NULL, (long)n_particles_per_bunch, NULL, beamline, 
                         &tsSpec->pCentral, NULL, NULL, NULL, NULL, run, 1,
                         0, 1, 0, NULL,
                         NULL, NULL, beam->lostOnPass, eptr);
      /*
      printf("nleft=%ld, tp0=%g, rp0=%g\n", n_left, tsSpec->pCentral, run->p_central);
      printf("betax=%g, alfax=%g, etax=%g, etaxp=%g\n", 
             tsptr->twiss[0][1], tsptr->twiss[0][0],tsptr->disp[0][0],tsptr->disp[0][1]);
      printf("betax=%g, alfax=%g, etax=%g, etaxp=%g\n", 
             tsptr->twiss[1][1], tsptr->twiss[1][0],tsptr->disp[1][0],tsptr->disp[1][1]);
      */
      dump_lost_particles(&SDDS_loss, beam->particle+n_left, beam->lostOnPass+n_left,
                          beam->n_to_track-n_left, 1);
      if (!SDDS_Terminate(&SDDS_loss)) {
        SDDS_SetError("Problem terminating 'losses' file (finish_output)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      }
      free_beamdata(beam);
    }
    eptr = eptr->succ; 
  }
  free(weight);
  return;
}

TSCATTER *initTSCATTER (ELEMENT_LIST *eptr)
{
  TSCATTER *tsptr;

  tsptr = ((TSCATTER*)eptr->p_elem);

  tsptr->twiss[0][0] = eptr->twiss->alphax;
  tsptr->twiss[0][1] = eptr->twiss->betax;
  tsptr->twiss[0][2] = (1.+ sqr(eptr->twiss->alphax))/eptr->twiss->betax;
  tsptr->twiss[1][0] = eptr->twiss->alphay;
  tsptr->twiss[1][1] = eptr->twiss->betay;
  tsptr->twiss[1][2] = (1.+ sqr(eptr->twiss->alphay))/eptr->twiss->betay;
  tsptr->twiss[2][0] = 0.0;
  tsptr->twiss[2][1] = sqr(tsSpec->gamma/tsSpec->sigma_p);
  tsptr->twiss[2][2] = sqr(1./tsSpec->sz0);
  tsptr->disp[0][0] = eptr->twiss->etax;
  tsptr->disp[0][1] = eptr->twiss->etapx;
  tsptr->disp[1][0] = eptr->twiss->etay;
  tsptr->disp[1][1] = eptr->twiss->etapy;

  tsptr->sigx = sqrt(tsSpec->emit[0]*tsptr->twiss[0][1]); 
  tsptr->sigy = sqrt(tsSpec->emit[1]*tsptr->twiss[1][1]);
  tsptr->sigz = tsSpec->sz0;
  tsptr->sigxyz = tsptr->sigx * tsptr->sigy * tsptr->sigz;

  tsptr->range[0] = 2.*tsSpec->range[0]*tsptr->sigx;
  tsptr->range[1] = 2.*tsSpec->range[1]*tsptr->sigy;
  tsptr->range[2] = 2.*tsSpec->range[2]*tsptr->sigz;
  tsptr->range[3] = 2.*tsSpec->range[0]*sqrt(tsSpec->emit[0]/tsptr->twiss[0][1]);
  tsptr->range[4] = 2.*tsSpec->range[1]*sqrt(tsSpec->emit[1]/tsptr->twiss[1][1]);
  tsptr->range[5] = 2.*tsSpec->range[2]*tsSpec->sigma_p/tsSpec->gamma;

  tsptr->factor = pow(PI, 2.0)
    *8.*pow(tsSpec->range[0], 3.0)*pow(tsSpec->range[1], 3.0)*pow(tsSpec->range[2], 3.0)/pow(PI, 6.0)
    *pow(tsSpec->charge/e_mks, 2.0)
    *c_mks/tsSpec->gamma
    *pow(re_mks, 2.0)/4./tsptr->sigxyz;

  tsptr->s_rate = tsptr->totalWeight = tsptr->ignorWeight = 0.;
  tsptr->simuCount = 0;

  return (tsptr);
}

void init_TSPEC ()
{
  if (!(tsSpec = SDDS_Malloc(sizeof(*tsSpec))))
    bomb("memory allocation failure at setup Touscheck scatter", NULL);                

  tsSpec->seed = seed;
  tsSpec->ebeam = p_central_mev;
  tsSpec->delta = delta;
  tsSpec->sigma_p = sigma_dp;
  tsSpec->sigz = sigma_s;
  tsSpec->charge = charge;
  tsSpec->emit[0] = emittance[0];
  tsSpec->emit[1] = emittance[1];
  tsSpec->emit[2] = 1.0;

  tsSpec->gamma = tsSpec->ebeam/me_mev;
  tsSpec->p0 = sqrt(sqr(tsSpec->ebeam) - me_mev*me_mev);
  tsSpec->pCentral = tsSpec->p0/me_mev;
  tsSpec->dp0 = delta*tsSpec->p0/tsSpec->gamma;
  tsSpec->sz0  = tsSpec->sigz*tsSpec->gamma;
  tsSpec->range[0] = distribution_cutoff[0];
  tsSpec->range[1] = distribution_cutoff[1];
  tsSpec->range[2] = distribution_cutoff[2];
  tsSpec->nbins = nbins;
  tsSpec->weight_limit = weight_limit;

  return;
}

/************************************************************************************\
 * Modified from S. Khan's code.                                                    *
 *      select two electrons (labelled a and b) for Touschek scattering             *
\************************************************************************************/

void selectPart(long iseed, TSCATTER *tsptr, double *p1, double *p2, 
                double *dens1, double *dens2)
{
  long iseed0 = 13578931;
  int i,j;
  double densa[3], densb[3];

  if (iseed > 0) 
    iseed0 = iseed;
  for (i=0; i<3; i++) {
    p1[i] = p2[i] = (random_1(iseed0)-0.5) * tsptr->range[i];
    j = i+3;
    p1[j] = (random_2(iseed0)-0.5) * tsptr->range[j]
      -tsptr->twiss[i][0]/tsptr->twiss[i][1]*p1[i];
    densa[i]=exp(-0.5*(tsptr->twiss[i][2]*p1[i]*p1[i]
                       +tsptr->twiss[i][0]*p1[i]*p1[j]*2.
                       +tsptr->twiss[i][1]*p1[j]*p1[j])/tsSpec->emit[i]);
    if(i!=0) {
      p2[j] = (random_2(iseed0)-0.5) * tsptr->range[j]
        -tsptr->twiss[i][0]/tsptr->twiss[i][1]*p2[i];
      densb[i]=exp(-0.5*(tsptr->twiss[i][2]*p2[i]*p2[i]
                         +tsptr->twiss[i][0]*p2[i]*p2[j]*2.
                         +tsptr->twiss[i][1]*p2[j]*p2[j])/tsSpec->emit[i]);
    }
  }
  /* Dispersion correction */
  p1[0] = p1[0] + p1[5]*tsptr->disp[0][0]*tsSpec->gamma;
  p1[3] = p1[3] + p1[5]*tsptr->disp[0][1]*tsSpec->gamma;

  p2[0] = p1[0] - p2[5]*tsptr->disp[0][0]*tsSpec->gamma;
  p2[3] = (random_2(iseed0)-0.5) * tsptr->range[3]
    -tsptr->twiss[0][0]/tsptr->twiss[0][1]*p2[0];

  densb[0]=exp(-0.5*(tsptr->twiss[0][2]*p2[0]*p2[0]
                     +tsptr->twiss[0][0]*p2[0]*p2[3]*2.
                     +tsptr->twiss[0][1]*p2[3]*p2[3])/tsSpec->emit[0]);
  p2[0] = p1[0];
  p2[3] = p2[3] + p2[5]*tsptr->disp[0][1]*tsSpec->gamma;

  *dens1 = densa[0] * densa[1] * densa[2]; 
  *dens2 = densb[0] * densb[1] * densb[2]; 

  return;
}

void bunch2cm(double *p1, double *p2, double *q, double *beta)
{
  double pp1, pp2, e1, e2, ee;
  int i;
  double betap1, bb, gamma, gg1, factor1;

  pp1=0.0;
  pp2=0.0;
  for(i=0; i<3; i++) {
    pp1 = pp1 + sqr(p1[i+3]);
    pp2 = pp2 + sqr(p2[i+3]);
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

  gamma=1./sqrt(1.-bb);
  gg1=gamma/(gamma+1);
  factor1=gamma*(gg1*betap1-e1);

  for(i=0; i<3; i++) {
    q[i]=p1[i+3]+factor1*beta[i];
  }

  return;
}

void eulertrans(double *v0, double theta, double phi, double *v1, double *v)
{
  double th, ph, s1, s2, c1, c2;
  double x0, y0, z0; 

  *v=sqrt(v0[0]*v0[0]+v0[1]*v0[1]+v0[2]*v0[2]);
  th=acos(v0[2]/(*v));
  ph=PIo2+atan2(v0[1],v0[0]);

  s1=sin(th);
  s2=sin(ph);
  c1=cos(th);
  c2=cos(ph);

  x0=sin(theta)*cos(phi);
  y0=sin(theta)*sin(phi);
  z0=cos(theta);

  v1[0]=*v*(c2*x0-c1*s2*y0+s1*s2*z0);
  v1[1]=*v*(s2*x0+c1*c2*y0-s1*c2*z0);
  v1[2]=*v*(         s1*y0   +c1*z0);

  return;
}

void cm2bunch(double *p1, double *p2, double *q, double *beta)
{
  int i;
  double qq, e, betaq, bb, gamma, gg1, factor1, factor2;

  qq=0.0;
  for(i=0; i<3; i++) {
    qq = qq + q[i]*q[i];
  }

  e=sqrt(me_mev*me_mev+qq);

  betaq=0.0;
  bb=0.0;
  for(i=0; i<3; i++) {
    betaq = betaq + beta[i]*q[i];
    bb = bb + beta[i]*beta[i];
  }
  gamma=1./sqrt(1.-bb);
  gg1=gamma/(gamma+1);
  factor1=gamma*( gg1*betaq+e);
  factor2=gamma*(-gg1*betaq+e);

  for(i=0; i<3; i++) {
    p1[i+3]= q[i]+factor1*beta[i];
    p2[i+3]=-q[i]+factor2*beta[i];
  }

  return;
}

double moeller(double vc, double theta, double phi)
{
  double cross; 
  double vc2, st2, x;

  vc2=vc*vc;
  st2=sqr(sin(theta));
  x=1./vc2;

  cross = (1-vc2)*(sqr(x+1)*(4./st2/st2-3./st2)+1+4./st2);
 
  return cross;
}

