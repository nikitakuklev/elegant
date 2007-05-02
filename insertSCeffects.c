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
#include "complex.h"

#define DEBUG 0

typedef struct {
  char *name, *type, *exclude;
} SC_SPEC;

static SC_SPEC *scSpec = NULL;
static long No_scSpec = 0;
static SPACE_CHARGE *sc = NULL;

void linearSCKick(double *coord, ELEMENT_LIST *eptr, double *center);
int  nonlinearSCKick(double *coord, ELEMENT_LIST *eptr, double *center, 
                     double sigmax, double sigmay, double *kick);

long getSCMULTSpecCount() 
{
  return (No_scSpec);
}
char *getSCMULTName()
{
  return (sc->name);
}

void addSCSpec(char *name, char *type, char *exclude)
{
  if (!(scSpec 
	= SDDS_Realloc(scSpec,
		       sizeof(*scSpec)*(No_scSpec+1))))
    bomb("memory allocation failure", NULL);
  scSpec[No_scSpec].name = NULL;
  scSpec[No_scSpec].type = NULL;
  scSpec[No_scSpec].exclude = NULL;
  if ((name &&
       !SDDS_CopyString(&scSpec[No_scSpec].name, name)) ||
      (type &&
       !SDDS_CopyString(&scSpec[No_scSpec].type, type)) ||
      (exclude &&
       !SDDS_CopyString(&scSpec[No_scSpec].exclude, exclude)))
    bomb("memory allocation failure", NULL);
  
  No_scSpec++;
}

void clearSCSpecs() 
{
  while (No_scSpec--) {
    if (scSpec[No_scSpec].name)
      free(scSpec[No_scSpec].name);
    if (scSpec[No_scSpec].type)
      free(scSpec[No_scSpec].type);
    if (scSpec[No_scSpec].exclude)
      free(scSpec[No_scSpec].exclude);
  }
  free(scSpec);
  scSpec = NULL;
}

long insertSCMULT(char *name, long type, long *occurrence) 
{
  long i;
  for (i=0; i<No_scSpec; i++) {
    if (scSpec[i].exclude && wild_match(name, scSpec[i].exclude))
      continue;
    if (scSpec[i].name && !wild_match(name, scSpec[i].name))
      continue;
    if (scSpec[i].type && !wild_match(entity_name[type], scSpec[i].type))
      continue;
    (*occurrence)++;
    break;
  }

  if (*occurrence < sc->nskip || sc->nskip==0)
    return(0);

  *occurrence = 0;
  return(1);
}

#include "insertSCeffects.h"
void setupSCEffect(NAMELIST_TEXT *nltext, RUN *run, LINE_LIST *beamline) 
{
  long i;
  
  if (!No_scSpec && !(sc = SDDS_Realloc(sc, sizeof(*sc))))
    bomb("memory allocation failure", NULL);                

  /* process the namelist text */
  set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
  set_print_namelist_flags(0);
  process_namelist(&insert_sceffects, nltext);
  print_namelist(stdout, &insert_sceffects);

  if (clear) {
    clearSCSpecs();
    if (!name && !type)
      return;
  }
  if (disable)
    return;
  
  if (!name || !strlen(name))
    bomb("no name given", NULL);
  str_toupper(name);
  if (has_wildcards(name) && strchr(name, '-'))
    name = expand_ranges(name);
  if (type) {
    str_toupper(type);
    if (has_wildcards(type) && strchr(type, '-'))
      type = expand_ranges(type);
    for (i=0; i<N_TYPES; i++)
      if (wild_match(entity_name[i], type))
	break;
    if (i==N_TYPES) {
      fprintf(stderr, "type pattern %s does not match any known type", type);
      exit(1);
    }
  }
  if (exclude) {
    str_toupper(exclude);
    if (has_wildcards(exclude) && strchr(exclude, '-'))
      exclude = expand_ranges(exclude);
  }
  
  addSCSpec(name, type, exclude);

  sc->name = "SCElem";
  if (element_prefix) 
    sc->name = element_prefix;

  sc->nskip = 0;
  sc->nonlinear = 0; 
  sc->horizontal = sc->vertical = sc->longitudinal =0;
  if (skip)
    sc->nskip = skip;

  if (vertical)
    sc->vertical = vertical;

  if (horizontal)
    sc->horizontal = horizontal;

  if (longitudinal)
    sc->longitudinal =longitudinal;

  if (nonlinear)
    sc->nonlinear = nonlinear;
}

/* track through space charge element */
void trackThroughSCMULT(double **part, long np, ELEMENT_LIST *eptr)
{
  long i;
  double *coord;
  double kx, ky, sx;
  double center[3], kick[2];
  double sigmax, sigmay;
  int flag;

  if (!np)  
    return;
  /* compute bunch center */
  for(i=center[0]=center[1]=center[2]=0; i<np; i++) {
    coord = part[i];
    center[0] += coord[0];
    center[1] += coord[2];
    center[2] += coord[4];
  }
  center[0] /= np;
  center[1] /= np;
  center[2] /= np;
 
  /* apply kick to particles */
  if (!nonlinear) {
    for(i=0; i<np; i++) {
      coord = part[i];
      linearSCKick(coord, eptr, center);
    }
  }
  else {
    sigmax = computeRmsCoordinate(part, 0, np);
    sigmay = computeRmsCoordinate(part, 2, np);
    for(i=0; i<np; i++) {
      coord = part[i];

      if ((fabs(coord[0]-center[0])<sigmax) && (fabs(coord[2]-center[1]) < sigmay)) {
        linearSCKick(coord, eptr, center);
        continue;
      }

      if (sigmax/sigmay>0.99 && sigmax/sigmay < 1.01) {
        sx = 0.99 * sigmay;
        flag = nonlinearSCKick(coord, eptr, center, sx, sigmay, kick); 
        if(!flag) {
          linearSCKick(coord, eptr, center);
          continue;
        }
        kx = kick[0];
        ky = kick[1];
        sx = 1.01 * sigmay;
        flag = nonlinearSCKick(coord, eptr, center, sx, sigmay, kick); 
        if(!flag) {
          linearSCKick(coord, eptr, center);
          continue;
        }
        kx += kick[0];
        ky += kick[1];
        coord[1] += kx / 2.0;
        coord[3] += ky / 2.0;
      }
      else {
        flag = nonlinearSCKick(coord, eptr, center, sigmax, sigmay, kick); 
        if(!flag) {
          linearSCKick(coord, eptr, center);
          continue;
        }
        coord[1] += kick[0];
        coord[3] += kick[1];
      }
    }
  }

  sc->dmux=sc->dmuy=0.0;       /* reset space charge strength */
}

void linearSCKick(double *coord, ELEMENT_LIST *eptr, double *center)
{
  double k0, kx, ky;
  k0 = sc->c1 * exp(-sqr(coord[4]-center[2])/sqr(sc->sigmaz)/2.0);
  if (sc->horizontal) {
    kx = k0 * sc->dmux / eptr->twiss->betax;	/* From dmux to KL */
    coord[1] += kx*(coord[0]-center[0]);
  }
  if (sc->vertical) {
    ky = k0 * sc->dmuy / eptr->twiss->betay;	/* From dmuy to KL */
    coord[3] += ky*(coord[2]-center[1]);
  }
}

int nonlinearSCKick(double *coord, ELEMENT_LIST *eptr, double *center, 
                     double sigmax, double sigmay, double *kick)
{
  double k0, kx, ky, sqs;
  COMPLEX wa, wb, w1, w2, w;
  double temp;
  long flag;

  k0 = sc->c1 * exp(-sqr(coord[4]-center[2])/sqr(sc->sigmaz)/2.0) * sqrt(PI/2.0);

  sqs = sqrt(fabs(sqr(sigmax)-sqr(sigmay))*2.0);
  kx = k0 * sc->dmux * sigmax * sqrt(sigmax+sigmay) / sqrt(fabs(sigmax-sigmay)) / eptr->twiss->betax;
  ky = k0 * sc->dmuy * sigmay * sqrt(sigmax+sigmay) / sqrt(fabs(sigmax-sigmay)) / eptr->twiss->betay;

  w1.r = (coord[0]-center[0])/sqs;
  w1.i = (coord[2]-center[1])/sqs;
  w2.r = (coord[0]-center[0])*sigmay/sigmax/sqs;
  w2.i = (coord[2]-center[1])*sigmax/sigmay/sqs;
  temp = exp((-sqr(coord[0]-center[0])/sqr(sigmax)-sqr(coord[2]-center[1])/sqr(sigmay))/2.0);

  wofz(&w1.r, &w1.i, &wa.r, &wa.i, &flag);
  if(!flag) return(0);
  wofz(&w2.r, &w2.i, &wb.r, &wb.i, &flag);
  if(!flag) return(0);

  w.r = wa.r - temp*wb.r;
  w.i = wa.i - temp*wb.i;
  kick[0] = kx * w.i;
  kick[1] = ky * w.r;

  return(1);
}

void initializeSCMULT(ELEMENT_LIST *eptr, double **part, long np, double Po, long i_pass )
{
  static CHARGE *charge;
	
  if (!eptr->twiss)
    bomb("Twiss parameters must be calculated before SC tracking.", NULL);
		
  if (i_pass==0) {
    while(eptr) {
      if (eptr->type==T_CHARGE) {
        charge = (CHARGE*)eptr->p_elem;
        break;
      }
      eptr =eptr->succ;
    }
    if (charge==NULL) 
      bomb("No charge element is given.", NULL);
	
  }
  sc->sigmax = computeRmsCoordinate(part, 0, np);
  sc->sigmay = computeRmsCoordinate(part, 2, np);
  sc->sigmaz = computeRmsCoordinate(part, 4, np);
  sc->c0 = sqrt(2.0/PI) * re_mks * charge->charge / e_mks;
  sc->c1 = sc->c0/sqr(Po)/sqrt(sqr(Po)+1.0)/sc->sigmaz;
  /*       printf("c0=%.6g, c1=%.6g, sz=%.6g\n\n", sc->c0, sc->c1, sc->sigmaz); */

  sc->dmux=sc->dmuy=0.0;
  sc->length=0.0;
}

void accumulateSCMULT(double **part, long np, ELEMENT_LIST *eptr)
{
  TWISS *twiss0;
  double dmux, dmuy, temp;
  double length;
	
  twiss0 = (eptr->pred)->twiss;
  temp = sc->sigmax + sc->sigmay;
  dmux = twiss0->betax / sc->sigmax / temp;
  dmuy = twiss0->betay / sc->sigmay / temp;
  sc->sigmax = computeRmsCoordinate(part, 0, np);
  sc->sigmay = computeRmsCoordinate(part, 2, np);
  twiss0 = eptr->twiss;
  temp = sc->sigmax + sc->sigmay;
  dmux += twiss0->betax / sc->sigmax / temp;
  dmuy += twiss0->betay / sc->sigmay / temp;
	
  length = ((DRIFT*)eptr->p_elem)->length;
  sc->dmux += dmux * length /2.0;
  sc->dmuy += dmuy * length /2.0;
}

double computeRmsCoordinate(double **coord, long i1, long np)
{
  double vrms, x, xc;
  long i;
  
  if (!np)
    return(0.0);

  /* compute centroids */
  for (i=xc=0; i<np; i++) {
    xc  += coord[i][i1];
  }
  xc  /= np;

  for (i=vrms=0; i<np; i++) {
    vrms += sqr(x  = coord[i][i1]-xc );
  }
  return(sqrt(vrms/np));
}

