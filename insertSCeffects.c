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

#define DEBUG 0

typedef struct {
  char *name, *type, *exclude;
} SC_SPEC;

static SC_SPEC *scSpec = NULL;
static long scSpecs = 0;
static SPACE_CHARGE *sc = NULL;

double get_rms(double **coord, long i1, long np);

long getSCMULTSpecCount() 
{
	return (scSpecs);
}
char *getSCMULTName()
{
	return (sc->name);
}

void addSCSpec(char *name, char *type, char *exclude)
{
  if (!(scSpec 
	= SDDS_Realloc(scSpec,
		       sizeof(*scSpec)*(scSpecs+1))))
    bomb("memory allocation failure", NULL);
  scSpec[scSpecs].name = NULL;
  scSpec[scSpecs].type = NULL;
  scSpec[scSpecs].exclude = NULL;
  if ((name &&
       !SDDS_CopyString(&scSpec[scSpecs].name, name)) ||
      (type &&
       !SDDS_CopyString(&scSpec[scSpecs].type, type)) ||
      (exclude &&
       !SDDS_CopyString(&scSpec[scSpecs].exclude, exclude)))
    bomb("memory allocation failure", NULL);
  
  scSpecs++;
}

void clearSCSpecs() 
{
  while (scSpecs--) {
    if (scSpec[scSpecs].name)
      free(scSpec[scSpecs].name);
    if (scSpec[scSpecs].type)
      free(scSpec[scSpecs].type);
    if (scSpec[scSpecs].exclude)
      free(scSpec[scSpecs].exclude);
  }
  free(scSpec);
  scSpec = NULL;
}

long insertSCMULT(char *name, long type, long *occurance) 
{
  long i;
  for (i=0; i<scSpecs; i++) {
    if (scSpec[i].exclude && wild_match(name, scSpec[i].exclude))
      continue;
    if (scSpec[i].name && !wild_match(name, scSpec[i].name))
      continue;
    if (scSpec[i].type && !wild_match(entity_name[type], scSpec[i].type))
      continue;
    (*occurance)++;
    break;
  }

  if(*occurance < sc->nskip || sc->nskip==0)
  	  	return(0);

  *occurance = 0;
  return(1);
}

#include "insertSCeffects.h"
void setupSCEffect(NAMELIST_TEXT *nltext, RUN *run, LINE_LIST *beamline) 
{
  long i;
  
  if (!scSpecs && !(sc = SDDS_Realloc(sc, sizeof(*sc))))
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

  sc->name = NULL;
  if(element_prefix) 
  	sc->name = element_prefix;

  sc->nskip = sc->nonlinear = 0; 
  sc->horizontal = sc->vertical = sc->longitudinal =0;
  if(skip)
  	sc->nskip = skip;

  if(vertical)
  	sc->vertical = vertical;

  if(horizontal)
  	sc->horizontal = horizontal;

  if(longitudinal)
  	sc->longitudinal =longitudinal;

  if(nonlinear)
  	sc->nonlinear = nonlinear;
}

/* track through space charge element */
void trackThroughSCMULT(double **part, long np, ELEMENT_LIST *eptr)
{
	long i;
	double *coord;
	double k0, kx, ky;
	double x0, y0, z0;
	static long j;
	
	if (!np)  
	    return;

	/* compute bunch center */
	for(i=x0=y0=z0=0; i<np; i++) {
		coord = part[i];
		x0 += coord[0];
		y0 += coord[2];
		z0 += coord[4];
	}
	x0 /= np;
	y0 /= np;
	z0 /= np;

	/* apply kick to particles */
	for(i=0; i<np; i++) {
		coord = part[i];
		k0 = exp(-sqr(coord[4]-z0)/sqr(sc->sigmaz)/2.0);
		if(sc->horizontal) {
			kx = sc->c1 * k0 * (4*PI * sc->dmux / eptr->twiss->betax);	/* From dmux to KL */
			coord[1] += kx*(coord[0]-x0);
		}
		if(sc->vertical) {
			ky = sc->c1 * k0 * (4*PI * sc->dmuy / eptr->twiss->betay);	/* From dmuy to KL */
			coord[3] += ky*(coord[2]-y0);
/*			if(i==0) printf("j=%d, dmuy=%.6g, z=%.6g, k0=%.6g, ky=%.6g, betay=%.6g, c0=%.6g, c1=%.6g\n", ++j, sc->dmuy, coord[4]-z0, k0, ky, eptr->twiss->betay, sc->c0, sc->c1); */
		}
	}
	sc->dmux=sc->dmuy=0.0;												/* reset space charge strength */
}

void initializeSCMULT(ELEMENT_LIST *eptr, double **part, long np, double Po, long i_pass )
{
	CHARGE *charge;
	
	if(!eptr->twiss)
		bomb("Twiss parameters must be calculated before SC tracking.", NULL);
		
	if(i_pass==0) {
		while(eptr) {
			if(eptr->type==T_CHARGE) {
				charge = (CHARGE*)eptr->p_elem;
				break;
			}
			eptr =eptr->succ;
		}
		if(charge==NULL) 
			bomb("No charge element is given.", NULL);
	
		sc->c0 = re_mks * charge->charge / e_mks / pow(2*PI, 3./2.);
		sc->sigmax = get_rms(part, 0, np);
		sc->sigmay = get_rms(part, 2, np);
		sc->sigmaz = get_rms(part, 4, np);
	}
	sc->c1 = sc->c0/pow(Po, 3.0)/sc->sigmaz;
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

	sc->sigmax = get_rms(part, 0, np);
	sc->sigmay = get_rms(part, 2, np);
	twiss0 = eptr->twiss;
	temp = sc->sigmax + sc->sigmay;
	dmux += twiss0->betax / sc->sigmax / temp;
	dmuy += twiss0->betay / sc->sigmay / temp;
	
	length = ((DRIFT*)eptr->p_elem)->length;
	sc->dmux += dmux * length /2.0;
	sc->dmuy += dmuy * length /2.0;
}

double get_rms(double **coord, long i1, long np)
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

