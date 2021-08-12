/*************************************************************************\
 * Copyright (c) 2003 The University of Chicago, as Operator of Argonne
 * National Laboratory.
 * Copyright (c) 2003 The Regents of the University of California, as
 * Operator of Los Alamos National Laboratory.
 * This file is distributed subject to a Software License Agreement found
 * in the file LICENSE that is included with this distribution. 
 \*************************************************************************/

/*  file: cooler.cc 
 *  based on: tfeedback.c
 *
 */


#if defined(SOLARIS) && !defined(__GNUC__)
#include <sunmath.h>
#endif

#include <complex>
#include "mdb.h"
#include "track.h"
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>
//#include <algorithm>

#define cms c_mks
#define twopi PIx2

// Turn of MPI implementation for this element because the incoherent effects required *all* particle coordinates
// on one processor. For that reason, CPICKUP and CKICKER are set to UNIPROCESSOR mode, so that all the work
// is done by the master processor.

#undef USE_MPI

void coolerPickup(CPICKUP *cpickup, double **part0, long np0, long pass, double Po, long idSlotsPerBunch)
{

  /*Initialize some variables */
  double sum;
  long i;
  double *time0 = NULL;           /* array to record arrival time of each particle */
  long *ibParticle = NULL;        /* array to record which bucket each particle is in */
  long **ipBucket = NULL;         /* array to record particle indices in part0 array for all particles in each bucket */
  long *npBucket = NULL;          /* array to record how many particles are in each bucket */
  long nBuckets=0;

  // Allocate space for particle coordinates 
  cpickup->long_coords = (double *)malloc(sizeof(double)*np0);
  cpickup->horz_coords = (double *)malloc(sizeof(double)*np0);
  cpickup->vert_coords = (double *)malloc(sizeof(double)*np0);
  

  // MPI Required vars 
#if USE_MPI
  long npTotal;
  double sumTotal;
  MPI_Status mpiStatus;

  // this element does nothing in single particle mode (e.g., trajectory, orbit, ..) 
  if (notSinglePart==0)
    return;
#endif

  // runs only between user specified passes
  if ((cpickup->startPass>0 && pass<cpickup->startPass) || (cpickup->endPass>0 && pass>cpickup->endPass))
    return;

  // makes sure element is initialized 
  if (cpickup->initialized==0)
    initializeCoolerPickup(cpickup);

  // assign buckets normal
  if (isSlave || !notSinglePart) 
    index_bunch_assignments(part0, np0, cpickup->bunchedBeamMode?idSlotsPerBunch:0, Po, &time0, &ibParticle, &ipBucket, &npBucket, &nBuckets, -1);

  // assign buckets mpi
#if USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
  if (myid==0)
    MPI_Recv(&nBuckets, 1, MPI_LONG, 1, 1, MPI_COMM_WORLD, &mpiStatus);
  else if (myid==1)
    MPI_Send(&nBuckets, 1, MPI_LONG, 0, 1, MPI_COMM_WORLD);
#endif


  // memory management for bucket assignments
  if (cpickup->updateInterval>1 && pass%cpickup->updateInterval!=0) {
    if (isSlave || !notSinglePart) 
      free_bunch_index_memory(time0, ibParticle, ipBucket, npBucket, nBuckets);
    return;
  }



  if (cpickup->nBunches==0 || cpickup->nBunches!=nBuckets) {
    if (cpickup->nBunches!=nBuckets) {
      printf("Number of bunches has changed, re-initializing cooler pickup\n");
      fflush(stdout);
    }
    cpickup->nBunches = nBuckets;
    cpickup->filterOutput = (double*)SDDS_Realloc(cpickup->filterOutput, sizeof(*cpickup->filterOutput)*nBuckets);
    for (i=0; i<nBuckets; i++) {
      cpickup->filterOutput[i] = 0;
    }
    cpickup->pass0 = pass;
  }



  // Important stuff happens here
  //    filterOutput takes average time of arrival
  //    part_coords takes each particles individual time of arival
  
  sum = 0; // Used to average all particle arival times in the pickup
  for(i=0; i<np0; i++){
    sum += part0[i][4];
    cpickup->long_coords[i] = part0[i][4]; // store the particle time in the coord array to be used in kicker 
    cpickup->horz_coords[i] = part0[i][0] + cpickup->dx;   
    cpickup->vert_coords[i] = part0[i][2] + cpickup->dy; 
  }
  

    
  


#if USE_MPI
  if (myid==0)
    np0 = 0;
  MPI_Allreduce(&sum, &sumTotal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&np0,  &npTotal, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  cpickup->filterOutput[0] = sumTotal/npTotal;
#else
  cpickup->filterOutput[0] = sum/np0;
#endif


  //memory managment
  if (isSlave || !notSinglePart) 
    free_bunch_index_memory(time0, ibParticle, ipBucket, npBucket, nBuckets);
}




// Pickup initialization function{
void initializeCoolerPickup(CPICKUP *cpickup)
{
  if (cpickup->ID==NULL || !strlen(cpickup->ID))
    bombElegant("you must give an ID string for CPICKUP", NULL);

  cpickup->nBunches = 0;

  if (cpickup->updateInterval<1)
    cpickup->updateInterval = 1;

  cpickup->initialized = 1;
}


struct E_params {double x; double y; double gamma; double lambda;};

double Ex (double theta, void * params){
  //unpacking parameters
  struct E_params * p = (struct E_params *) params;
  double x = (p->x);
  double y = (p->y);
  double gamma = (p->gamma);
  double lambda = (p->lambda);
	
  // define some parameters for calc
  double k0 = twopi/lambda;
  double rho = sqrt(x*x + y*y);
  double phi = atan2(y,x); 
	
  double J0_term = gsl_sf_bessel_Jn(0, rho*k0*theta/(1 + pow(gamma*theta,2)) );  
  double J2_term = gsl_sf_bessel_Jn(2, rho*k0*theta/(1 + pow(gamma*theta,2)) ) * pow(gamma*theta,2) * cos(2*phi); 
  double I0_term = theta / pow(1 + pow(gamma*theta,2) ,4);

  return (J0_term + J2_term) * I0_term;
}

struct indexed_coord {int id; double t; double x; double y;};

int compare_indexed_coord(const void *arg1, const void *arg2){
  indexed_coord const *lhs = static_cast<indexed_coord const*>(arg1);
  indexed_coord const *rhs = static_cast<indexed_coord const*>(arg2);
  return (lhs->t < rhs->t) ? -1 :  ((rhs->t < lhs->t) ? 1 : 0);
}

void coolerKicker(CKICKER *ckicker, double **part0, long np0, LINE_LIST *beamline, long pass, long nPasses, char *rootname, double Po, long idSlotsPerBunch)
{
  long i,j;
  double *time0 = NULL;     /* array to record arrival time of each particle */
  long *ibParticle = NULL;  /* array to record which bucket each particle is in */
  long **ipBucket = NULL;		/* array to record particle indices in part0 array for all particles in each bucket */
  long *npBucket = NULL;    /* array to record how many particles are in each bucket */
  long iBucket, nBuckets=0;
  long updateInterval;


  double Ex0;
  if (ckicker->transverseMode!=0) {
    // integration variables	
    double error;
    size_t neval;

    struct E_params params = {0, 0, Po, ckicker->lambda_rad}; // {x, y, gamma, lambda}

    gsl_function F;
    F.function = &Ex;
    F.params = &params;

    // (gsl_function F, a, b, epsabs, epsrel, result, error, neval)
    gsl_integration_qng(&F, 0, ckicker->angle_rad, 0, 1e-7, &Ex0, &error, &neval);
  }


#if USE_MPI
  MPI_Status mpiStatus;
  double sumTotal;
  long npTotal;

  if (notSinglePart==0)
    /* this element does nothing in single particle mode (e.g., trajectory, orbit, ..) */
    return;
#endif

  if ((ckicker->startPass>0 && pass<ckicker->startPass) || (ckicker->endPass>0 && pass>ckicker->endPass))
    return;

  if (isSlave || !notSinglePart) 
    index_bunch_assignments(part0, np0, ckicker->bunchedBeamMode?idSlotsPerBunch:0, Po, &time0, &ibParticle, &ipBucket, &npBucket, &nBuckets, -1);

#if USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
  if (myid==0)
    MPI_Recv(&nBuckets, 1, MPI_LONG, 1, 1, MPI_COMM_WORLD, &mpiStatus);
  else if (myid==1)
    MPI_Send(&nBuckets, 1, MPI_LONG, 0, 1, MPI_COMM_WORLD);
#endif

  if (ckicker->initialized==0) {
    initializeCoolerKicker(ckicker, beamline, nPasses*nBuckets, rootname);
  }


  if (ckicker->startPass>0 && ckicker->startPass!=ckicker->pickup->startPass)
    bombElegantVA((char*)"CKICKER linked to CPICKUP with different START_PASS value (%ld vs %ld).", 
                  ckicker->startPass, ckicker->pickup->startPass);
  if (ckicker->endPass>0 && ckicker->endPass!=ckicker->pickup->endPass)
    bombElegantVA((char*)"CKICKER linked to CPICKUP with different END_PASS value (%ld vs %ld).", 
                  ckicker->endPass, ckicker->pickup->endPass);

  if ((updateInterval =  ckicker->pickup->updateInterval*ckicker->updateInterval)<=0) 
    bombElegantVA((char*)"CKICKER and CPICKUP with ID=%s have UPDATE_INTERVAL product of %d", ckicker->ID, updateInterval);
  if (pass%updateInterval!=0) {
    if (isSlave || !notSinglePart) 
      free_bunch_index_memory(time0, ibParticle, ipBucket, npBucket, nBuckets);
    return;
  }


  if (ckicker->nBunches==0 || ckicker->nBunches!=nBuckets) {
    if (ckicker->nBunches!=nBuckets) {
      printf("Number of bunches has changed, re-initializing cooler driver.\n");
      fflush(stdout);
    }
    ckicker->nBunches = nBuckets;
    ckicker->pass0 = pass;
  } 

  if (ckicker->nBunches!=ckicker->pickup->nBunches)
    bombElegant("mismatch in number of buckets between CKICKER and CPICKUP", NULL);

  
  
  // ================================= //
  // This is where the kick is applied //
  // ================================= //
  double sum, nom_kick;
  double s_avg, s_i, x_i, y_i;
  
  indexed_coord *incoherent_ar;
  incoherent_ar = (indexed_coord *)malloc(sizeof(indexed_coord)*np0);

  sum = 0;
  for(i=0; i<np0; i++){
    sum += part0[i][4];
    
    if(ckicker->incoherentMode!=0){
      incoherent_ar[i].id = i;
      incoherent_ar[i].x = ckicker->pickup->horz_coords[i];
      incoherent_ar[i].y = ckicker->pickup->vert_coords[i];
      incoherent_ar[i].t = ckicker->pickup->long_coords[i];
    }  
  }
  
  if(ckicker->incoherentMode!=0){
    qsort(incoherent_ar, np0, sizeof(indexed_coord),compare_indexed_coord ); 
  }


  // s_avg equals average arrival time at kicker - average time at pickup
#if USE_MPI
  if (myid==0)
    np0 = 0;
  MPI_Allreduce(&sum, &sumTotal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&np0,  &npTotal, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  s_avg = sumTotal/npTotal - ckicker->pickup->filterOutput[0];
#else
  s_avg = sum/np0 - ckicker->pickup->filterOutput[0];
#endif
  
  for(i=0; i<np0; i++){
    long pid = i;
    if (ckicker->incoherentMode!=0){
      long pid = incoherent_ar[i].id;
      
      for(j=i+1; j<np0-i; j++){   
        //if(abs(incoherent_ar[j].t-incoherent_ar[i].t) < ckicker->Nu*ckicker->lambda_rad){  
        if(abs(incoherent_ar[j].t-incoherent_ar[i].t) < 11*ckicker->lambda_rad){ 
          double incoherent_phase = abs(incoherent_ar[j].t-incoherent_ar[i].t)*twopi/ckicker->lambda_rad;
          double incoherent_strength = 1 - abs(incoherent_phase / (twopi*11));
          double kick = -(ckicker->strength * incoherent_strength) * sin(incoherent_phase + ckicker->phase*twopi); 
          
          part0[pid][5] += kick;
        }else{
          break;
        }   
      }
    }
      
    //particle coords
    s_i = part0[pid][4] - ckicker->pickup->long_coords[pid];
    x_i = part0[pid][0] - ckicker->magnification * ckicker->pickup->horz_coords[pid];
    y_i = part0[pid][2] - ckicker->magnification * ckicker->pickup->vert_coords[pid];   

    // on-axis kick = sin(difference in path length from avg + phase)
    nom_kick = -(ckicker->strength) * sin((s_avg - s_i) * (twopi/(ckicker->lambda_rad)) + ckicker->phase*twopi);     
    
		
    
    if (ckicker->transverseMode!=0) {
      double Exi, error;
      size_t neval;
      struct E_params params = {x_i, y_i, Po, ckicker->lambda_rad}; // {x, y, gamma, lambda}
			
      gsl_function F;
      F.function = &Ex;
      F.params = &params;

      gsl_integration_qng(&F, 0, ckicker->angle_rad/Po, 0, 1e-7, &Exi, &error, &neval);
      part0[pid][5] += nom_kick * Exi/Ex0;
    }else{
      part0[pid][5] += nom_kick;
    }   
  }//*/
 

#if USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  if (isSlave || !notSinglePart) 
    free_bunch_index_memory(time0, ibParticle, ipBucket, npBucket, nBuckets);

}

// Kicker initialization function

void initializeCoolerKicker(CKICKER *ckicker, LINE_LIST *beamline, long nPasses, char *rootname)
{
  ELEMENT_LIST *eptr;
  long pickupFound = 0;

  if (ckicker->ID==NULL || !strlen(ckicker->ID))
    bombElegant("you must give an ID string for CKICKER", NULL);

  eptr = beamline->elem;
  while (eptr) {
    if (eptr->type==T_CPICKUP && strcmp(ckicker->ID, ((CPICKUP*)eptr->p_elem)->ID)==0) {
      pickupFound = 1;
      ckicker->pickup = ((CPICKUP*)eptr->p_elem);
      break;
    }
    eptr = eptr->succ;
  }
  if (!pickupFound) 
    bombElegant("pickup not found for CKICKER", NULL);

  ckicker->nBunches = 0;

  if (ckicker->updateInterval<1)
    ckicker->updateInterval = 1;
  ckicker->initialized = CKICKER_MAIN_INIT;
}



