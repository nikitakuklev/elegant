/*************************************************************************\
* Copyright (c) 2017 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2017 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: gasScattering.c
 * purpose: Do tracking to find slope aperture starting from the end of each element.
 *
 * Michael Borland, 2017
 */

#include "mdb.h"
#include "track.h"
#include "gasScattering.h"

static SDDS_DATASET SDDSsa;
static long fireOnPass = 1;

#if USE_MPI
void gatherLostParticles(double ***lostParticles, long *nLost, long nSurvived, long n_processors, int myid);

static long nElements;
static ELEMENT_LIST **elementArray = NULL;
static double betax0, betay0;
static double alphax0, alphay0;
static long iElementName, iElementOccurence, iElementType, ixpOrig, iypOrig, iupOrig, ivpOrig, isOrig,
    ixLost, iyLost, ideltaLost, isLost;

static void slopeOffsetFunction(double **coord, long np, long pass, long i_elem, long n_elem, ELEMENT_LIST *eptr, double *pCentral)
{
  long ix, iy, id, ie, ip, particleID;
  MALIGN mal;
#ifdef DEBUG
  if (i_elem==0)
    fprintf(fpd, "processor %d has %ld particles on pass %ld\n",
            myid, np, pass);
#endif
  if (pass==fireOnPass) {
    for (ie=0; ie<nElements; ie++) {
      if (eptr==elementArray[ie])
        break;
    }
    if (ie==nElements) return;
    elementArray[ie] = eptr;
    mal.dxp = mal.dyp = 0;
    mal.dz = mal.dt = mal.de = mal.dx = mal.dy = 0;
    mal.startPID = mal.endPID = -1;
    for (ip=0; ip<np; ip++) {
      if ((particleID = coord[ip][6])<0) {
#ifdef DEBUG
        fprintf(fpd, "particleID = %ld, excluded\n", particleID);
#endif
        continue;
      }
      id = particleID%(nx*ny);
      if ((particleID-id)/(nx*ny)!=ie) {
#ifdef DEBUG
        fprintf(fpd, "ie = %ld, computed ie is %ld, excluded\n", ie, (particleID-id)/nDelta);
#endif
        continue;
      }
      if (id>nx*ny)
        bombElegant("invalid id value (>nx*ny)", NULL);
      ix = id%ny;
      iy = id/ny;
      mal.dxp = ix*(xpmax-xpmin)/(nx-1.0) + xpmin;
      mal.dyp = iy*(ypmax-ypmin)/(ny-1.0) + ypmin;
      if (twiss_scaling) {
        mal.dxp *= sqrt(betax0/elementArray[ie]->twiss->betax);
        mal.dyp *= sqrt(betay0/elementArray[ie]->twiss->betay);
      }
      offset_beam(coord+ip, 1, &mal, *pCentral);
#ifdef DEBUG
      if (fpd)
        fprintf(fpd, "Imparted error (%le, %le) to particle %ld (ie=%ld, pId=%ld) on processor %d\n", 
                mal.xp, mal.yp, ip, ie, (long)coord[ip][6], myid);
#endif
    }
#ifdef DEBUG
    if (fpd) {
      fprintf(fpd, "Slope kick done on processor %d\n", myid);
      fflush(fpd);
    }
#endif
  }
}
#endif

void setupGasScattering(
                                 NAMELIST_TEXT *nltext,
                                 RUN *run,
                                 VARY *control,
                                 long twissFlag
                                 )
{
  char description[200];
  
#if !USE_MPI
  bombElegant("gas_scattering command is not available in serial elegant.", NULL);
#endif

#if USE_MPI
  /* process namelist input */
  set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
  set_print_namelist_flags(0);
  if (processNamelist(&gas_scattering, nltext)==NAMELIST_ERROR)
    bombElegant(NULL, NULL);
  if (echoNamelists) print_namelist(stdout, &gas_scattering);
  if (twiss_scaling && !twissFlag)
    bombElegant("You must compute twiss parameters if twiss_scaling is desired", NULL);

  if (run->concat_order!=0)
    bombElegant("at present, gas_scattering is incompatible with concatenation", NULL);
  
  /* check for data errors */
  if (!output)
    bombElegant("no output filename specified", NULL);
  if (xpmin >= xpmax)
    bombElegant("xpmin >= xpmax",  NULL);
  if (ypmin >= ypmax)
    bombElegant("ypmin >= ypmax",  NULL);
  if (s_start>=s_end)
    bombElegant("s_start >= s_end", NULL);
  if (include_name_pattern && has_wildcards(include_name_pattern) && strchr(include_name_pattern, '-'))
    include_name_pattern = expand_ranges(include_name_pattern);
  if (include_type_pattern && has_wildcards(include_type_pattern) && strchr(include_type_pattern, '-'))
    include_type_pattern = expand_ranges(include_type_pattern);
  
  nElements = 0;

  output = compose_filename(output, run->rootname);
  sprintf(description, "Slope aperture search");

  if (myid==0) {
    if (!SDDS_InitializeOutput(&SDDSsa, SDDS_BINARY, 1, description, "momentum aperture",  output)) {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      exitElegant(1);
    }
    
    if ((iElementName=SDDS_DefineColumn(&SDDSsa, "ElementName", NULL, NULL, NULL, NULL, SDDS_STRING, 0))<0 ||
        (iElementType=SDDS_DefineColumn(&SDDSsa, "ElementType", NULL, NULL, NULL, NULL, SDDS_STRING, 0))<0 ||
        (iElementOccurence=SDDS_DefineColumn(&SDDSsa, "ElementOccurence", NULL, NULL, NULL, NULL, SDDS_LONG, 0))<0 ||
        (isOrig=SDDS_DefineColumn(&SDDSsa, "s", NULL, "m", NULL, NULL, SDDS_DOUBLE, 0))<0 ||
        (ixpOrig=SDDS_DefineColumn(&SDDSsa, "xp", NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0))<0 ||
        (iypOrig=SDDS_DefineColumn(&SDDSsa, "yp", NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0))<0 ||
        (iupOrig=SDDS_DefineColumn(&SDDSsa, "up", NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0))<0 ||
        (ivpOrig=SDDS_DefineColumn(&SDDSsa, "vp", NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0))<0 ||
        (ixLost=SDDS_DefineColumn(&SDDSsa, "xLost", NULL, "m", NULL, NULL, SDDS_DOUBLE, 0))<0 ||
        (iyLost=SDDS_DefineColumn(&SDDSsa, "yLost", NULL, "m", NULL, NULL, SDDS_DOUBLE, 0))<0 ||
        (ideltaLost=SDDS_DefineColumn(&SDDSsa, "deltaLost", NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0))<0 ||
        (isLost=SDDS_DefineColumn(&SDDSsa, "sLost", NULL, "m", NULL, NULL, SDDS_DOUBLE, 0))<0 ||
        SDDS_DefineParameter(&SDDSsa, "Step", NULL, NULL, NULL, NULL, SDDS_LONG, NULL)<0 ||
        SDDS_DefineParameter(&SDDSsa, "SVNVersion", NULL, NULL, "SVN version number", NULL, SDDS_STRING, SVN_VERSION)<0) {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      exitElegant(1);
    }
   
    if (!SDDS_SaveLayout(&SDDSsa) || !SDDS_WriteLayout(&SDDSsa)){
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      exitElegant(1);
    }
  }
#endif
}

void finishGasScattering()
{
#if USE_MPI
  if (SDDS_IsActive(&SDDSsa) && !SDDS_Terminate(&SDDSsa)) {
    SDDS_SetError("Problem terminating SDDS output (finishGasScattering)");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
#endif
}

long runGasScattering(
                              RUN *run,
                              VARY *control,
                              ERRORVAL *errcon,
                              LINE_LIST *beamline,
                              double *startingCoord
                              )
{    
#if USE_MPI
  double **coord;
  long nTotal, id, ip, ie, ix, iy, nLost, nLeft, nElem, nEachProcessor, code, iRow;
  long nWorkingProcessors = n_processors - 1;
  double **lostParticles;
  double pCentral;
  ELEMENT_LIST *elem, *elem0;
#if MPI_DEBUG
  short mpiDebug = 1;
#else
  short mpiDebug = 0;
#endif
#endif

#if !USE_MPI
  bombElegant("Can't do gas scattering in serial elegant.", NULL);
#endif

#if USE_MPI
  if (myid==0 || mpiDebug) {
    printf("Started multi-particle gas scattering algorithm\n");
    fflush(stdout);
  }
  
  elem = &(beamline->elem);
  if (twiss_scaling) {
    betax0 = elem->twiss->betax;
    alphax0 = elem->twiss->alphax;
    betay0 = elem->twiss->betay;
    alphay0 = elem->twiss->alphay;
  } else {
    betax0 = betay0 = 1;
    alphax0 = alphay0 = 0;
  }

  /* determine how many elements will be tracked */
  elem0 = NULL;
  nElem = 0;
  elementArray = NULL;
  while (elem && elem->end_pos<s_end) {
    if (elem->end_pos>=s_start && elem->end_pos<=s_end &&
        (!include_name_pattern || wild_match(elem->name, include_name_pattern)) &&
        (!include_type_pattern || wild_match(entity_name[elem->type], include_type_pattern)) ) {
      if (!elem0)
	elem0 = elem;
      nElem++;
      elementArray = SDDS_Realloc(elementArray, sizeof(*elementArray)*(nElem+10));
      elementArray[nElem-1] = elem;
    }
    elem = elem->succ;
  }
  nElements = nElem;
  if (nElem==0) 
    SDDS_Bomb("no elements found for gas scattering computation");
  if (verbosity>0 && (myid==0 || mpiDebug)) {
    printf("%ld elements selected for tracking\n", nElements);
    fflush(stdout);
  }
    
  nElem = nElements;
  nTotal = nElem*nx*ny;
  if (nTotal%nWorkingProcessors!=0) {
    printf("Warning: The number of working processors (%ld) does not evenly divide into the number of particles (nx=%ld, ny=%ld, nElem=%ld)\n",
           nWorkingProcessors, nx, ny, nElem);
    fflush(stdout);
    nEachProcessor =  (nTotal/nWorkingProcessors)+1;
  } else {
    nEachProcessor = nTotal/nWorkingProcessors;
  }
  
  if (myid==0 || mpiDebug) {
    printf("nTotal = %ld, nWorkingProcessors = %ld, nx = %ld, ny = %ld, nElements = %ld, nEachProcessor = %ld\n",
           nTotal, nWorkingProcessors, nx, ny, nElements, nEachProcessor);
    fflush(stdout);
  }
  
  if (myid!=0) {
    /* allocate and initialize array for tracking */
    coord = (double**)czarray_2d(sizeof(**coord), nEachProcessor, COORDINATES_PER_PARTICLE);
    lostParticles = (double**)czarray_2d(sizeof(**lostParticles), nEachProcessor, COORDINATES_PER_PARTICLE+1);	 
  } else {
    coord = (double**)czarray_2d(sizeof(**coord), 1, COORDINATES_PER_PARTICLE);
    lostParticles = (double**)czarray_2d(sizeof(**lostParticles), nEachProcessor*nWorkingProcessors, COORDINATES_PER_PARTICLE);
  }

  if (control->n_passes==1)
    fireOnPass = 0;
  else
    fireOnPass = 1;
  
  setTrackingOmniWedgeFunction(NULL);
  if (startingCoord)
    memcpy(coord[0], startingCoord, sizeof(double)*6);
  else
    memset(coord[0], 0, sizeof(**coord)*6);
  coord[0][6] = 1;
  pCentral = run->p_central;
  if (verbosity>1) {
    printf("Tracking fiducial particle\n");
    fflush(stdout);
  }
  delete_phase_references();
  reset_special_elements(beamline, RESET_INCLUDE_ALL&~RESET_INCLUDE_RANDOM);
  code = do_tracking(NULL, coord, 1, NULL, beamline, &pCentral, 
                     NULL, NULL, NULL, NULL, run, control->i_step, 
                     FIRST_BEAM_IS_FIDUCIAL+(verbosity>4?0:SILENT_RUNNING)+INHIBIT_FILE_OUTPUT, 1, 0, NULL, NULL, NULL, NULL, NULL);
  if (!code) {
    if (myid==0)
      printf("Fiducial particle lost. Don't know what to do.\n");
    exitElegant(1);
  }
  if (myid==0 || mpiDebug) {
    printf("Fiducial particle tracked.\n");
    fflush(stdout);
  }
  
  MPI_Barrier(MPI_COMM_WORLD);
  if (myid==0 || mpiDebug) {
    printf("Fiducalization completed\n");
    fflush(stdout);
  }

  nLost = 0;
  nLeft = 0;
  if (myid!=0) {
    nLeft = nEachProcessor;
    for (ip=0; ip<nLeft; ip++) 
      lostParticles[ip][6] = -2;
    for (ip=0; ip<nLeft; ip++) {
      if (startingCoord)
        memcpy(coord[ip], startingCoord, sizeof(**coord)*6);
      else
        memset(coord[ip], 0, sizeof(**coord)*6);
      coord[ip][6] = (myid-1)*nEachProcessor + ip;
      if (coord[ip][6]>=nTotal) {
        /* Don't track more buffer particles than needed */
        coord[ip][6] = -1;
        nLeft = ip+1;
      }
    }
    setTrackingOmniWedgeFunction(slopeOffsetFunction); 
    if (verbosity>1 && mpiDebug) {
      printf("Tracking\n");
      fflush(stdout);
    }
    nLost = nLeft;
    nLeft = do_tracking(NULL, coord, nLeft, NULL, beamline, &pCentral, 
                        NULL, NULL, NULL, NULL, run, control->i_step, 
                        FIDUCIAL_BEAM_SEEN+FIRST_BEAM_IS_FIDUCIAL+SILENT_RUNNING+INHIBIT_FILE_OUTPUT, 
                        control->n_passes, 0, NULL, NULL, NULL, lostParticles, NULL);
    nLost -= nLeft;
    printf("Done tracking nLeft = %ld, nLost = %ld\n", nLeft, nLost);
    fflush(stdout);
    if (verbosity>9) {
      for (ip=0; ip<nLost; ip++)
        printf("ip = %ld, particleID = %ld\n", ip, (long)lostParticles[ip][6]);
      fflush(stdout);
    }
    setTrackingOmniWedgeFunction(NULL);
  }

  if (verbosity>3 && (myid==0 || mpiDebug)) {
    printf("Waiting on barrier after tracking\n");
    fflush(stdout);
  }    

  MPI_Barrier(MPI_COMM_WORLD);

  nLeft = 0;
  gatherLostParticles(&lostParticles, &nLost, nLeft, n_processors, myid);
  if (myid==0 || mpiDebug) {
    printf("Lost-particle gather done, nLost = %ld\n", nLost); 
    for (ip=0; ip<nLost; ip++) {
      printf("ip=%ld, particleID=%ld\n", ip, (long)lostParticles[ip][6]);
    }
    fflush(stdout);
  }
  
  if (myid==0) {
    if (!SDDS_StartPage(&SDDSsa, nLost) ||
        !SDDS_SetParameters(&SDDSsa, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "Step",
                            control->i_step, NULL)) {
      SDDS_SetError("Problem writing SDDS table (doSlopeApertureSearch)");
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
    
    for (ip=iRow=0; ip<nLost; ip++) {
      if (verbosity>5) {
        printf("Processing ip=%ld, particleID=%ld\n", ip, (long)lostParticles[ip][6]);
        fflush(stdout);
      }
      if (lostParticles[ip][6]<0) {
        if (lostParticles[ip][6]==-2) {
          long j;
          printf("Problem with lost-particle accounting\n");
          for (j=0; j<COORDINATES_PER_PARTICLE+1; j++)
            printf("coord[%ld] = %le\n", j, lostParticles[ip][j]);
          bombElegant("problem with lost particle accounting!", NULL);
        }
        /* buffer particle, ignore */
        continue;
      }
      
      /* Figure out (ix, iy, ie) */
      particleID = lostParticles[ip][6];
      id = particleID%(nx*ny);
      ie = (particleID-id)/(nx*ny);
      ix = id%ny;
      iy = id/ny;
      if (verbosity>5) {
        printf("id = %ld, ie = %ld, ix = %ld, iy = %ld\n", 
               id, ie, ix, iy);
        fflush(stdout);
      }
      if (!SDDS_SetRowValues(&SDDSsa, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, iRow++,
                             iElementName, elementArray[ie]->name,
                             iElementType, entity_name[elementArray[ie]->type],
                             iElementOccurence, elementArray[ie]->occurence, 
                             isOrig, elementArray[ie]->end_pos,
                             ixpOrig, 
                             (ix*(xpmax-xpmin)/(nx-1.0) + xpmin)*(twiss_scaling?sqrt(betax0/elementArray[ie]->twiss->betax):1),
                             iypOrig,
                             (iy*(ypmax-ypmin)/(ny-1.0) + ypmin)*(twiss_scaling?sqrt(betay0/elementArray[ie]->twiss->betay):1),
                             ixLost, lostParticles[ip][0],
                             iyLost, lostParticles[ip][2],
                             ideltaLost, lostParticles[ip][5],
                             isLost, lostParticles[ip][4], 
                             -1)) {
        SDDS_SetError("Problem setting row values in SDDS table (doSlopeApertureSearch)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      }
      if (verbosity>5) {
        printf("Set data for row %ld: %s, %s, %ld, %le, %le, %le, %le %le, %le, %le\n",
               iRow-1,
               elementArray[ie]->name, entity_name[elementArray[ie]->type], elementArray[ie]->occurence,
               elementArray[ie]->end_pos, ix*(xpmax-xpmin)/(nx-1.0) + xpmin, iy*(ypmax-ypmin)/(ny-1.0) + ypmin,
               lostParticles[ip][0], lostParticles[ip][2], lostParticles[ip][5], lostParticles[ip][4]);
        fflush(stdout);
      }
    }
    if (verbosity>5) {
      printf("About to write page\n");
      fflush(stdout);
    }
    if (!SDDS_WritePage(&SDDSsa)) {
      SDDS_SetError("Problem writing SDDS table (doSlopeApertureSearch)");
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
    if (!inhibitFileSync)
      SDDS_DoFSync(&SDDSsa);
    free_czarray_2d((void**)coord, 1, COORDINATES_PER_PARTICLE);
    free_czarray_2d((void**)lostParticles, nEachProcessor*nWorkingProcessors, COORDINATES_PER_PARTICLE+1);	 
  }
  else {
    free_czarray_2d((void**)coord, nEachProcessor, COORDINATES_PER_PARTICLE);
    free_czarray_2d((void**)lostParticles, nEachProcessor, COORDINATES_PER_PARTICLE+1);	 
  }

  if (verbosity>3 && (myid==0 || mpiDebug)) {
    printf("Waiting on barrier after file operations\n");
    fflush(stdout);
  }
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  return 1;
}


