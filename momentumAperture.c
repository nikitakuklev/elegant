/*************************************************************************\
* Copyright (c) 2006 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2006 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: momentumAperture.c
 * purpose: Do tracking to find momentum aperture starting from the end of each element.
 * Ref: M. Belgrounne et al. PAC203, 896.
 *
 * Michael Borland, 2006
 */

#include "mdb.h"
#include "track.h"
#include "momentumAperture.h"

static SDDS_DATASET SDDSma;
static double momentumOffsetValue = 0;
static long fireOnPass = 1;
static double **turnByTurnCoord = NULL;
static long turnsStored = 0;

#include "fftpackC.h"
long determineTunesFromTrackingData(double *tune, double **turnByTurnCoord, long turns, double delta);
long multiparticleLocalMomentumAcceptance(RUN *run, VARY *control, ERRORVAL *errcon, LINE_LIST *beamline, double *startingCoord);
#if USE_MPI
void gatherLostParticles(double ***lostParticles, long *nLost, long nSurvived, long n_processors, int myid);
#endif

static void momentumOffsetFunction(double **coord, long np, long pass, double *pCentral)
{
  MALIGN mal;
  
  if (pass==fireOnPass) {
    mal.dxp = mal.dyp = mal.dz = mal.dt = mal.de = 0;
    mal.dx = x_initial;
    mal.dy = y_initial;
    mal.dp = momentumOffsetValue;
    mal.startPID = mal.endPID = -1;
    offset_beam(coord, np, &mal, *pCentral);
    turnsStored = 0;
  }
  if (turnByTurnCoord) {
    /* N.B.: the arrays are ordered differently than normal, e.g.,
     * turnByTurnCoord[0] is the data for x, not particle/turn 0
     */
    turnByTurnCoord[0][turnsStored] = coord[0][0]; 
    turnByTurnCoord[1][turnsStored] = coord[0][1]; 
    turnByTurnCoord[2][turnsStored] = coord[0][2]; 
    turnByTurnCoord[3][turnsStored] = coord[0][3]; 
    turnByTurnCoord[4][turnsStored] = coord[0][5]; 
    turnsStored++;
  }
}

#ifdef DEBUG
FILE *fpd = NULL;
#endif
static long nElements;
static ELEMENT_LIST **elementArray = NULL;
#if USE_MPI
static double deltaStep;
static long nDelta;

static void momentumOffsetFunctionOmni(double **coord, long np, long pass, long i_elem, long n_elem, ELEMENT_LIST *eptr, double *pCentral)
{
  long id, ie, ip, particleID;
  MALIGN mal;
#ifdef DEBUG
#if USE_MPI
  if (i_elem==0)
    fprintf(fpd, "processor %d has %ld particles on pass %ld\n",
            myid, np, pass);
#endif
#endif
  if (pass==fireOnPass) {
    for (ie=0; ie<nElements; ie++) {
      if (eptr==elementArray[ie])
        break;
    }
    if (ie==nElements) return;
    elementArray[ie] = eptr;
    mal.dxp = mal.dyp = mal.dz = mal.dt = mal.de = 0;    
    mal.dx = x_initial;
    mal.dy = y_initial;
    mal.startPID = mal.endPID = -1;
    for (ip=0; ip<np; ip++) {
      if ((particleID = coord[ip][6])<0) {
#ifdef DEBUG
        fprintf(fpd, "particleID = %ld, excluded\n", particleID);
#endif
        continue;
      }
      id = particleID%nDelta;
      if ((particleID-id)/nDelta!=ie) {
#ifdef DEBUG
        fprintf(fpd, "ie = %ld, computed ie is %ld, excluded\n", ie, (particleID-id)/nDelta);
#endif
        continue;
      }
      if (id>nDelta) 
        bombElegant("invalid id value (>nDelta)", NULL);
      mal.dp = delta_negative_limit + id*deltaStep;
      offset_beam(coord+ip, 1, &mal, *pCentral);
#ifdef DEBUG
#if USE_MPI
      if (fpd)
        fprintf(fpd, "Imparted error %le to particle %ld (ie=%ld, pId=%ld) on processor %d\n", mal.dp, ip, ie, (long)coord[ip][6], myid);
#endif
#endif
    }
#ifdef DEBUG
#if USE_MPI
    if (fpd) {
      fprintf(fpd, "Momentum kick done on processor %d\n", myid);
      fflush(fpd);
    }
#endif
#endif
  }
}
#endif

void setupMomentumApertureSearch(
                                 NAMELIST_TEXT *nltext,
                                 RUN *run,
                                 VARY *control
                                 )
{
  char description[200];
  
  /* process namelist input */
  set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
  set_print_namelist_flags(0);
  if (processNamelist(&momentum_aperture, nltext)==NAMELIST_ERROR)
    bombElegant(NULL, NULL);
  if (echoNamelists) print_namelist(stdout, &momentum_aperture);

  if (run->concat_order!=0)
    bombElegant("at present, momentum_aperture is incompatible with concatenation", NULL);
  
  /* check for data errors */
  if (!output)
    bombElegant("no output filename specified", NULL);
  if (delta_negative_limit>=0)
    bombElegant("delta_negative_limit >= 0", NULL);
  if (delta_positive_limit<=0) 
    bombElegant("delta_positive_limit <= 0", NULL);
  if (delta_negative_start>0)
    bombElegant("delta_negative_start > 0", NULL);
  if (delta_positive_start<0) 
    bombElegant("delta_positive_start < 0", NULL);
  if (delta_step_size<=0)
    bombElegant("delta_step_size <= 0", NULL);
  if (fabs(delta_negative_limit-delta_negative_start)<=delta_step_size/2)
    bombElegant("|delta_negative_limit-delta_negative_start| <= delta_step_size/2", NULL);
  if (delta_positive_limit-delta_positive_start<=delta_step_size/2)
    bombElegant("delta_positive_limit-delta_positive_start <= delta_step_size/2", NULL);
  if (splits<0)
    bombElegant("splits < 0", NULL);
  if (s_start>=s_end)
    bombElegant("s_start >= s_end", NULL);
  if (include_name_pattern && has_wildcards(include_name_pattern) && strchr(include_name_pattern, '-'))
    include_name_pattern = expand_ranges(include_name_pattern);
  if (include_type_pattern && has_wildcards(include_type_pattern) && strchr(include_type_pattern, '-'))
    include_type_pattern = expand_ranges(include_type_pattern);
  if (skip_elements<0)
    bombElegant("skip_elements < 0", NULL);
  if (process_elements<=0)
    bombElegant("process_elements <= 0", NULL);
  
  nElements = 0;
  if (output_mode==2) {
    /* Massively parallel algorithm using acceptance feature */
    if (include_name_pattern && has_wildcards(include_name_pattern))
      bombElegant("Can't have wildcards in name pattern for output_mode=2", NULL);
    if (include_type_pattern)
      bombElegant("Can't used include_type_pattern for output_mode=2", NULL);
  }
  
  output = compose_filename(output, run->rootname);
  sprintf(description, "Momentum aperture search");
#if SDDS_MPI_IO
  if (output_mode!=2) {
    SDDS_MPI_Setup(&SDDSma, 1, n_processors, myid, MPI_COMM_WORLD, 1);
    if (!SDDS_Parallel_InitializeOutput(&SDDSma, description, "momentum aperture",  output)) {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      exitElegant(1);
    }
  } else if (myid==0) {
    /* master only writes */
    if (!SDDS_InitializeOutput(&SDDSma, SDDS_BINARY, 1, description, "momentum aperture",  output)) {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      exitElegant(1);
    }
  }
#else
  if (!SDDS_InitializeOutput(&SDDSma, SDDS_BINARY, 1, description, "momentum aperture",  output)) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    exitElegant(1);
  }
#endif

  switch (output_mode) {
  case 1:
    if (SDDS_DefineColumn(&SDDSma, "ElementName", NULL, NULL, NULL, NULL, SDDS_STRING, 0)<0 ||
        SDDS_DefineColumn(&SDDSma, "s", NULL, "m", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
        SDDS_DefineColumn(&SDDSma, "ElementType", NULL, NULL, NULL, NULL, SDDS_STRING, 0)<0 ||
        SDDS_DefineColumn(&SDDSma, "ElementOccurence", NULL, NULL, NULL, NULL, SDDS_LONG, 0)<0 ||
        SDDS_DefineColumn(&SDDSma, "direction", NULL, NULL, NULL, NULL, SDDS_SHORT, 0)<0 ||
        SDDS_DefineColumn(&SDDSma, "deltaFound", NULL, NULL, NULL, NULL, SDDS_SHORT, 0)<0 ||
        SDDS_DefineColumn(&SDDSma, "delta", "$gd$R$blimit$n", NULL, NULL, NULL, SDDS_DOUBLE, 0)<0 ||
        SDDS_DefineColumn(&SDDSma, "lostOnPass", NULL, NULL, NULL, NULL, SDDS_LONG, 0)<0 ||
        SDDS_DefineColumn(&SDDSma, "sLost", NULL, "m", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
        SDDS_DefineColumn(&SDDSma, "xLost", NULL, "m", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
        SDDS_DefineColumn(&SDDSma, "yLost", NULL, "m", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
        SDDS_DefineColumn(&SDDSma, "deltaLost", NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0)<0 ||
        SDDS_DefineColumn(&SDDSma, "nuxLost", NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0)<0 ||
        SDDS_DefineColumn(&SDDSma, "nuyLost", NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0)<0 ||
        SDDS_DefineParameter(&SDDSma, "Step", NULL, NULL, NULL, NULL, SDDS_LONG, NULL)<0 ||
        SDDS_DefineParameter(&SDDSma, "SVNVersion", NULL, NULL, "SVN version number", NULL, SDDS_STRING, SVN_VERSION)<0) {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      exitElegant(1);
    }
    break;
  case 2:
#if USE_MPI
/*
    if (skip_elements>0)
      bombElegant("skip_elements can't be non-zero for output_mode=2.", NULL);
    if (s_start>0)
      bombElegant("s_start can't be non-zero for output_mode=2.", NULL);
 */
    if (myid==0) {
      if (SDDS_DefineColumn(&SDDSma, "ElementName", NULL, NULL, NULL, NULL, SDDS_STRING, 0)<0 ||
          SDDS_DefineColumn(&SDDSma, "s", NULL, "m", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
          SDDS_DefineColumn(&SDDSma, "ElementType", NULL, NULL, NULL, NULL, SDDS_STRING, 0)<0 ||
          SDDS_DefineColumn(&SDDSma, "ElementOccurence", NULL, NULL, NULL, NULL, SDDS_LONG, 0)<0 ||
          SDDS_DefineColumn(&SDDSma, "deltaPositive", "$gd$R$bpos$n", NULL, NULL, NULL, SDDS_DOUBLE, 0)<0 ||
          SDDS_DefineColumn(&SDDSma, "deltaNegative", "$gd$R$bneg$n", NULL, NULL, NULL, SDDS_DOUBLE, 0)<0 ||
          SDDS_DefineParameter(&SDDSma, "Step", NULL, NULL, NULL, NULL, SDDS_LONG, NULL)<0 ||
          SDDS_DefineParameter(&SDDSma, "SVNVersion", NULL, NULL, "SVN version number", NULL, SDDS_STRING, SVN_VERSION)<0) {
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
        exitElegant(1);
      }
    }
#else
    bombElegant("output_mode=2 not available in serial elegant.",  NULL);
#endif
    break;
  default:
    if (SDDS_DefineColumn(&SDDSma, "ElementName", NULL, NULL, NULL, NULL, SDDS_STRING, 0)<0 ||
        SDDS_DefineColumn(&SDDSma, "s", NULL, "m", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
        SDDS_DefineColumn(&SDDSma, "ElementType", NULL, NULL, NULL, NULL, SDDS_STRING, 0)<0 ||
        SDDS_DefineColumn(&SDDSma, "ElementOccurence", NULL, NULL, NULL, NULL, SDDS_LONG, 0)<0 ||
        SDDS_DefineColumn(&SDDSma, "deltaPositiveFound", NULL, NULL, NULL, NULL, SDDS_SHORT, 0)<0 ||
        SDDS_DefineColumn(&SDDSma, "deltaPositive", "$gd$R$bpos$n", NULL, NULL, NULL, SDDS_DOUBLE, 0)<0 ||
        SDDS_DefineColumn(&SDDSma, "lostOnPassPositive", NULL, NULL, NULL, NULL, SDDS_LONG, 0)<0 ||
        SDDS_DefineColumn(&SDDSma, "sLostPositive", NULL, "m", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
        SDDS_DefineColumn(&SDDSma, "xLostPositive", NULL, "m", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
        SDDS_DefineColumn(&SDDSma, "yLostPositive", NULL, "m", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
        SDDS_DefineColumn(&SDDSma, "deltaLostPositive", NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0)<0 ||
        SDDS_DefineColumn(&SDDSma, "nuxLostPositive", NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0)<0 ||
        SDDS_DefineColumn(&SDDSma, "nuyLostPositive", NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0)<0 ||
        SDDS_DefineColumn(&SDDSma, "deltaNegativeFound", NULL, NULL, NULL, NULL, SDDS_SHORT, 0)<0 ||
        SDDS_DefineColumn(&SDDSma, "deltaNegative", "$gd$R$bneg$n", NULL, NULL, NULL, SDDS_DOUBLE, 0)<0 ||
        SDDS_DefineColumn(&SDDSma, "lostOnPassNegative", NULL, NULL, NULL, NULL, SDDS_LONG, 0)<0 ||
        SDDS_DefineColumn(&SDDSma, "sLostNegative", NULL, "m", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
        SDDS_DefineColumn(&SDDSma, "xLostNegative", NULL, "m", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
        SDDS_DefineColumn(&SDDSma, "yLostNegative", NULL, "m", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
        SDDS_DefineColumn(&SDDSma, "deltaLostNegative", NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0)<0 ||
        SDDS_DefineColumn(&SDDSma, "nuxLostNegative", NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0)<0 ||
        SDDS_DefineColumn(&SDDSma, "nuyLostNegative", NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0)<0 ||
        SDDS_DefineParameter(&SDDSma, "Step", NULL, NULL, NULL, NULL, SDDS_LONG, NULL)<0 ||
        SDDS_DefineParameter(&SDDSma, "SVNVersion", NULL, NULL, "SVN version number", NULL, SDDS_STRING, SVN_VERSION)<0) {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      exitElegant(1);
    }
    break;
  }
  
  
#if !SDDS_MPI_IO
  /* In the version with parallel IO, the layout will be written later */
  if(!SDDS_SaveLayout(&SDDSma) || !SDDS_WriteLayout(&SDDSma)){
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    exitElegant(1);
  }
#else 
  if (output_mode==2 && myid==0) {
    /* In output_mode==2, the master does the IO */
    if(!SDDS_SaveLayout(&SDDSma) || !SDDS_WriteLayout(&SDDSma)){
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      exitElegant(1);
    }
  }
#endif  
}

void finishMomentumApertureSearch()
{
  if (SDDS_IsActive(&SDDSma) && !SDDS_Terminate(&SDDSma)) {
    SDDS_SetError("Problem terminating SDDS output (finishMomentumApertureSearch)");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
}

long doMomentumApertureSearch(
                              RUN *run,
                              VARY *control,
                              ERRORVAL *errcon,
                              LINE_LIST *beamline,
                              double *startingCoord
                              )
{    
  double **coord;
  long nElem, iElem;
  double deltaInterval, pCentral, deltaStart;
  ELEMENT_LIST *elem, *elem0;
  long side;
  double **lostParticles;	
  short **loserFound, *direction=NULL, **survivorFound;
  int32_t **lostOnPass, *ElementOccurence;
  double deltaLimit1[2], deltaLimit, **deltaWhenLost, delta;
  double deltaStart1[2];
  double **xLost, **yLost, **deltaSurvived, **sLost, deltaLost;
  double **xTuneSurvived, **yTuneSurvived;
  double *sStart;
  double nominalTune[2], tune[2];
  char **ElementName, **ElementType;
  long code, outputRow, jobCounter;
  long processElements, skipElements, deltaSign, split, slot;
  char s[1000];
  unsigned long fiducial_flag_save;
#if USE_MPI
#if defined(DEBUG)
  FILE *fpdMpi = NULL;
  sprintf(s, "%s-debug-%02d", output, myid);
  if (!(fpdMpi = fopen(s, "w")))
    fpdMpi = NULL;
  fprintf(fpdMpi, "ID = %d\n", myid);
  fflush(fpdMpi);
#endif
  notSinglePart = 0;
#else
#if defined(DEBUG)
  FILE *fpdeb = NULL;

  sprintf(s, "%s-Debug", output);
  fprintf(stderr, "Debug file %s\n", s);
  
  if (!(fpdeb = fopen(s, "w")))
    bombElegant("unable to open debug file for momentum aperture scan", NULL);
  fprintf(fpdeb, "SDDS1\n");
  fprintf(fpdeb, "&parameter name=ElementName type=string &end\n");
  fprintf(fpdeb, "&parameter name=Side type=string &end\n");
  fprintf(fpdeb, "&column name=delta type=double &end\n");
  fprintf(fpdeb, "&column name=Lost type=short &end\n");
  fprintf(fpdeb, "&column name=Loss Pass type=short &end\n");
  fprintf(fpdeb, "&data mode=ascii no_row_counts=1 &end\n");
  fflush(fpdeb);
#endif
#endif

  /* determine how many elements will be tracked */
  elem = &(beamline->elem);
  elem0 = NULL;
  nElem = 0;
  elementArray = NULL;
  skipElements = skip_elements;
  while (elem && nElem<process_elements && elem->end_pos<s_end) {
    if (elem->end_pos>=s_start && elem->end_pos<=s_end &&
        (!include_name_pattern || wild_match(elem->name, include_name_pattern)) &&
        (!include_type_pattern || wild_match(entity_name[elem->type], include_type_pattern)) ) {
      if (skipElements>0) {
        skipElements --;
        elem = elem->succ;
        continue;
      }
      if (!elem0)
	elem0 = elem;
      nElem++;
      elementArray = SDDS_Realloc(elementArray, sizeof(*elementArray)*(nElem+10));
      elementArray[nElem-1] = elem;
    }
    elem = elem->succ;
  }
  nElements = nElem;

#if !USE_MPI
  if (nElem==0) 
    SDDS_Bomb("no elements found between s_start and s_end for momentum aperture computation");
#else
  switch (output_mode) {
  case 1:
    if ((2*nElem<n_processors) && (myid==0)) {
      printf("Warning: The number of elements should be larger than the number of processors to avoid wasting resources.\nThe effective number of elements is %ld. The number of processors is %d.\n", 
             2*nElem, n_processors);
    }    
    break;
  case 2:
    break;
  default:
    if (nElem<n_processors && (myid==0)) {
      printf("Warning: The number of elements should be larger than the number of processors to avoid wasting resources.\nThe number of elements is %ld. The number of processors is %d.\n",
             nElem, n_processors);
    }    
    break;
  }
  
  if (verbosity) {
    verbosity = 0;
    if (myid == 0)
      printf ("Warning: In parallel version, limited intermediate information will be provided\n");
  }
#endif

  if (output_mode==2) {
    printf("Branching multiple particle mode\n");
    multiparticleLocalMomentumAcceptance(run, control, errcon, beamline, startingCoord);
    printf("Returning from LMA search main routine.\n");
    fflush(stdout);
    return 1;
  }
  printf("Running in search mode\n");
  
  /* allocate arrays for tracking */
  coord = (double**)czarray_2d(sizeof(**coord), 1, 7);

  /* allocate arrays for storing data for negative and positive momentum limits for each element */
  lostOnPass = (int32_t**)czarray_2d(sizeof(**lostOnPass), (output_mode?1:2), (output_mode?2:1)*nElem);
  loserFound = (short**)czarray_2d(sizeof(**loserFound), (output_mode?1:2), (output_mode?2:1)*nElem);
  survivorFound = (short**)czarray_2d(sizeof(**survivorFound), (output_mode?1:2), (output_mode?2:1)*nElem);
  deltaSurvived = (double**)czarray_2d(sizeof(**deltaSurvived), (output_mode?1:2), (output_mode?2:1)*nElem);
  xTuneSurvived = (double**)czarray_2d(sizeof(**xTuneSurvived), (output_mode?1:2), (output_mode?2:1)*nElem);
  yTuneSurvived = (double**)czarray_2d(sizeof(**yTuneSurvived), (output_mode?1:2), (output_mode?2:1)*nElem);
  xLost = (double**)czarray_2d(sizeof(**xLost), (output_mode?1:2), (output_mode?2:1)*nElem);
  yLost = (double**)czarray_2d(sizeof(**yLost), (output_mode?1:2), (output_mode?2:1)*nElem);
  deltaWhenLost = (double**)czarray_2d(sizeof(**deltaWhenLost), (output_mode?1:2), (output_mode?2:1)*nElem);
  sLost = (double**)czarray_2d(sizeof(**sLost), (output_mode?1:2), (output_mode?2:1)*nElem);
  sStart = (double*)tmalloc(sizeof(*sStart)*(output_mode?2:1)*nElem);
  ElementName = (char**)tmalloc(sizeof(*ElementName)*(output_mode?2:1)*nElem);
  ElementType = (char**)tmalloc(sizeof(*ElementType)*(output_mode?2:1)*nElem);
  ElementOccurence = (int32_t*)tmalloc(sizeof(*ElementOccurence)*(output_mode?2:1)*nElem);
  if (output_mode)
    direction = (short*)tmalloc(sizeof(*direction)*2*nElem);
  deltaLimit1[0] = delta_negative_limit;
  deltaLimit1[1] = delta_positive_limit;
  deltaStart1[0] = delta_negative_start;
  deltaStart1[1] = delta_positive_start;

  if (control->n_passes==1)
    fireOnPass = 0;
  else
    fireOnPass = 1;
  turnByTurnCoord = (double**)czarray_2d(sizeof(double), 5, control->n_passes);
  
  /* need to do this because do_tracking() in principle may realloc this pointer */
  lostParticles = (double**)czarray_2d(sizeof(double),1, 8);	 
 
  elem = elem0;
  iElem = 0;
  processElements = process_elements;
  skipElements = skip_elements;
  
  /* Prevent do_tracking() from recognizing these flags. Instead, we'll control behavior directly */
  fiducial_flag_save = beamline->fiducial_flag;
  
  beamline->fiducial_flag = 0;
  
  if (fiducialize || forbid_resonance_crossing) {
    long i;
    if (startingCoord)
      memcpy(coord[0], startingCoord, sizeof(double)*6);
    else
      memset(coord[0], 0, sizeof(**coord)*6);
    coord[0][6] = 1;
    pCentral = run->p_central;
    if (verbosity>1) 
      printf("Tracking fiducial particle\n");
    delete_phase_references();
    reset_special_elements(beamline, RESET_INCLUDE_ALL&~RESET_INCLUDE_RANDOM);
    code = do_tracking(NULL, coord, 1, NULL, beamline, &pCentral, 
                       NULL, NULL, NULL, NULL, run, control->i_step, 
                       FIRST_BEAM_IS_FIDUCIAL+(verbosity>1?0:SILENT_RUNNING)+INHIBIT_FILE_OUTPUT, 1, 0, NULL, NULL, NULL, lostParticles, NULL);
    if (!code) {
      printf("Fiducial particle lost. Don't know what to do.\n");
      exitElegant(1);
    }
    printf("Fiducialization completed\n");
    if (verbosity>2)
      for (i=0; i<6; i++) 
        printf("%le%s", coord[0][i], i==5?"\n":", ");
  }

  outputRow = -1;
  jobCounter = -1;

#if USE_MPI
  verbosity = 0;  
#endif

  while (elem && processElements>0) {
#ifdef DEBUG
    printf("checking element %s#%ld\n", elem->name, elem->occurence); fflush(stdout);
#endif
    if ((!include_name_pattern || wild_match(elem->name, include_name_pattern)) &&
        (!include_type_pattern || wild_match(entity_name[elem->type], include_type_pattern))) {
      if (elem->end_pos>s_end) 
        break;
#ifdef DEBUG
      printf("including element %s#%ld\n", elem->name, elem->occurence); fflush(stdout);
#endif
      if (output_mode==0) {
#if USE_MPI
        jobCounter++;
        if (myid != jobCounter%n_processors) {
          elem = elem->succ;
          continue;
        }
#endif
        outputRow++;
      }
#if !USE_MPI
      if (verbosity>0) {
        printf("Searching for energy aperture for %s #%ld at s=%em\n", elem->name, elem->occurence, elem->end_pos);
        fflush(stdout);
      }
#endif
      for (side=0; side<2; side++) {
        if (output_mode==1) {
#if USE_MPI
          jobCounter++;
          if (myid!=jobCounter%n_processors)
            continue;
#endif
          outputRow++;
          slot = 0;
          direction[outputRow] = (side==0?-1:1);
        } else
          slot = side;

#if USE_MPI
        if (myid==0) {
          if (output_mode==1) {
            sprintf(s, "About %.3g%% done: ", (jobCounter*50.0)/nElem);
            report_stats(stdout, s);
            fflush(stdout);
          } else if (side==0) {
            sprintf(s, "About %.3g%% done: ", (jobCounter*100.0)/nElem);
            report_stats(stdout, s);
            fflush(stdout);
          }
        }
#else
	if (output_mode==1) {
	  sprintf(s, "About %.3g%% done: ", (outputRow*50.0)/nElem);
          report_stats(stdout, s);
          fflush(stdout);
	} else if (side==0) {
	  sprintf(s, "About %.3g%% done: ", (outputRow*100.0)/nElem);
          report_stats(stdout, s);
          fflush(stdout);
        }
        
#endif

        ElementName[outputRow] = elem->name;
        ElementType[outputRow] = entity_name[elem->type];
        ElementOccurence[outputRow] = elem->occurence;
        sStart[outputRow] = elem->end_pos;
        deltaStart = deltaStart1[side];
        deltaSign = side==0 ? -1 : 1;
        lostOnPass[slot][outputRow] = -1;
        loserFound[slot][outputRow] = survivorFound[slot][outputRow] = 0;
        xLost[slot][outputRow] = yLost[slot][outputRow] = 
          deltaWhenLost[slot][outputRow] = sLost[slot][outputRow] = 
            deltaSurvived[slot][outputRow] =  xTuneSurvived[slot][outputRow] =
              yTuneSurvived[slot][outputRow] = 0;
        deltaLost = deltaSign*DBL_MAX/2;
        deltaInterval = delta_step_size*deltaSign;
        
        if (verbosity>1) {
          printf(" Searching for %s side from %e toward %e with interval %e\n", side==0?"negative":"positive",
                  deltaStart, deltaLimit1[side], delta_step_size);
          fflush(stdout);
        }

        if (forbid_resonance_crossing) {
          momentumOffsetValue = 0;
          setTrackingWedgeFunction(momentumOffsetFunction, 
                                   elem->succ?elem->succ:elem0); 
          if (startingCoord)
            memcpy(coord[0], startingCoord, sizeof(double)*6);
          else
            memset(coord[0], 0, sizeof(**coord)*6);
          coord[0][6] = 1;
          pCentral = run->p_central;
          code = do_tracking(NULL, coord, 1, NULL, beamline, &pCentral, 
                             NULL, NULL, NULL, NULL, run, control->i_step, 
                             (fiducialize?FIDUCIAL_BEAM_SEEN+FIRST_BEAM_IS_FIDUCIAL:0)+(verbosity>4?0:SILENT_RUNNING)+INHIBIT_FILE_OUTPUT,
	                     control->n_passes, 0, NULL, NULL, NULL, lostParticles, NULL);
          if (!code || !determineTunesFromTrackingData(nominalTune, turnByTurnCoord, turnsStored, 0.0)) {
            printf("Fiducial particle tune is undefined.\n");
            exitElegant(1);
          }
          if (verbosity>3)
            printf("  Nominal tunes: %e, %e\n", nominalTune[0], nominalTune[1]);
        }
        
        deltaLimit = deltaLimit1[side];
        for (split=0; split<=splits; split++) {
          delta = deltaStart;
          
#if defined(DEBUG)
#if !USE_MPI
          fprintf(fpdeb, "%s\n%s\n",
                  elem->name, side==0?"negative":"positive");
          fflush(fpdeb);
#else
          fprintf(fpdMpi, "%s\n%s\n",
                  elem->name, side==0?"negative":"positive");
          fflush(fpdMpi);
#endif
#endif

          while (fabs(delta) <= fabs(deltaLimit)) {
            setTrackingWedgeFunction(momentumOffsetFunction, 
                                     elem->succ?elem->succ:elem0); 
            momentumOffsetValue = delta;
            if (startingCoord)
              memcpy(coord[0], startingCoord, sizeof(double)*6);
            else
              memset(coord[0], 0, sizeof(**coord)*6);
            coord[0][6] = 1;
            pCentral = run->p_central;
            if (verbosity>3) {
              printf("  Tracking with delta0 = %e (%e, %e, %e, %e, %e, %e), pCentral=%e\n", 
                      delta, coord[0][0], coord[0][1], coord[0][2], coord[0][3], coord[0][4], coord[0][5],
                      pCentral);
              fflush(stdout);
            }
            lostParticles[0][7] = -1;
            if (!fiducialize) {
              delete_phase_references();
              reset_special_elements(beamline, RESET_INCLUDE_ALL&~RESET_INCLUDE_RANDOM);
            }
            code = do_tracking(NULL, coord, 1, NULL, beamline, &pCentral, 
                               NULL, NULL, NULL, NULL, run, control->i_step, 
                               (fiducialize?FIDUCIAL_BEAM_SEEN:0)+FIRST_BEAM_IS_FIDUCIAL+(verbosity>4?0:SILENT_RUNNING)+(allow_watch_file_output?0:INHIBIT_FILE_OUTPUT), 
			       control->n_passes, 0, NULL, NULL, NULL, lostParticles, NULL);
            if (code && turnsStored>2) {
              if (!determineTunesFromTrackingData(tune, turnByTurnCoord, turnsStored, delta)) {
                if (forbid_resonance_crossing) {
                  if (verbosity>3)
                    printf("   Resonance crossing detected (no tunes).  Particle lost\n");
                  code = 0; /* lost */
                }
              } else {
                if (verbosity>3) 
                  printf("   Tunes: %e, %e\n", tune[0], tune[1]);
                if (forbid_resonance_crossing &&
                    ( (((long)(2*tune[0])) - ((long)(2*nominalTune[0])))!=0 ||
                     (((long)(2*tune[1])) - ((long)(2*nominalTune[1])))!=0) ) {
                  /* crossed integer or half integer */
                  if (verbosity>3)
                    printf("   Resonance crossing detected (%e, %e -> %e, %e).  Particle lost\n",
                            nominalTune[0], nominalTune[1], tune[0], tune[1]);
                  code = 0;
                }
              }
            } else 
              tune[0] = tune[1] = -1;
            if (!code) {
              /* particle lost */
              if (verbosity>3) {
                long i;
                printf("  Particle lost with delta0 = %e at s = %e\n", delta, coord[0][4]);
                if (verbosity>4)
                  for (i=0; i<6; i++)
                    printf("   coord[%ld] = %e\n", i, coord[0][i]);
                fflush(stdout);
              }
              lostOnPass[slot][outputRow] = lostParticles[0][7];
              xLost[slot][outputRow] = coord[0][0];
              yLost[slot][outputRow] = coord[0][2];
              sLost[slot][outputRow] = coord[0][4];
              deltaLost = delta;
              deltaWhenLost[slot][outputRow] = (coord[0][5]-pCentral)/pCentral;
              loserFound[slot][outputRow] = 1;
              break;
            } else {
              if (verbosity>2)
                printf("  Particle survived with delta0 = %e\n", delta);
              if (verbosity>3)
                printf("     Final coordinates: %le, %le, %le, %le, %le, %le\n",
                        coord[0][0], coord[0][1], coord[0][2], coord[0][3], coord[0][4], coord[0][5]);
              deltaSurvived[slot][outputRow] = delta;
              survivorFound[slot][outputRow] = 1;
              xTuneSurvived[slot][outputRow] = tune[0];
              yTuneSurvived[slot][outputRow] = tune[1];
            }
            delta += deltaInterval;
          } /* delta search */
          if (split==0) {
            if (!survivorFound[slot][outputRow]) {
	      if (!soft_failure) {
		printf("Error: No survivor found for initial scan for  %s #%ld at s=%em\n", elem->name, elem->occurence, elem->end_pos);
		exit(1);
	      }
              printf("Warning: No survivor found for initial scan for  %s #%ld at s=%em\n", elem->name, elem->occurence, elem->end_pos);
	      deltaSurvived[slot][outputRow] = 0;
	      survivorFound[slot][outputRow] = 1;
	      split = splits;
            }
            if (!loserFound[slot][outputRow]) {
              if (!soft_failure) {
                printf("Error: No loss found for initial scan for  %s #%ld at s=%em\n", elem->name, elem->occurence, elem->end_pos);
                exitElegant(1);
              }
	      loserFound[slot][outputRow] = 1;
	      split = splits;
            }
          }
          deltaStart = deltaSurvived[slot][outputRow] - steps_back*deltaInterval;
          deltaInterval /= split_step_divisor;
          deltaStart += deltaInterval;
          deltaLimit = deltaLost;
          if ((deltaStart<0 && deltaSign==1) || (deltaStart>0 && deltaSign==-1))
            deltaStart = 0;
        } /* split loop */
#if defined(DEBUG)
#if !USE_MPI
        fprintf(fpdeb, "\n");
        fflush(fpdeb);
#else
        fprintf(fpdMpi, "\n");
        fflush(fpdMpi);
#endif
#endif
        if (verbosity>0) {
          printf("Energy aperture for %s #%ld at s=%em is %e\n", elem->name, elem->occurence, elem->end_pos,
                  deltaSurvived[slot][outputRow]);
          fflush(stdout);
        }
      } /* side loop */
      
      processElements --;
    } /* element loop */
    elem = elem->succ;
  } 

  outputRow++;

#if USE_MPI
#if defined(DEBUG)
  fprintf(fpdMpi, "Waiting on barrier after main loop.\n");
  fflush(fpdMpi);
#endif
  MPI_Barrier(MPI_COMM_WORLD);
#endif

#if SDDS_MPI_IO
  /* Open file here for parallel IO */
  if (!SDDS_LayoutWritten(&SDDSma)) {
    if (!SDDS_MPI_File_Open(SDDSma.MPI_dataset, SDDSma.layout.filename, SDDS_MPI_WRITE_ONLY)) 
      SDDS_MPI_BOMB("SDDS_MPI_File_Open failed.", &SDDSma.MPI_dataset->MPI_file);
    if (!SDDS_MPI_WriteLayout(&SDDSma))  
      SDDS_MPI_BOMB("SDDS_MPI_WriteLayout failed.", &SDDSma.MPI_dataset->MPI_file);
  }
#endif                                                                                       
  if (!SDDS_StartPage(&SDDSma, outputRow) ||
      !SDDS_SetParameters(&SDDSma, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "Step",
                          control->i_step, NULL)) {
    SDDS_SetError("Problem writing SDDS table (doMomentumApertureSearch)");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }

  if ((output_mode==0 && 
       (!SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, ElementName, outputRow, "ElementName") ||
        !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, sStart, outputRow, "s") ||
        !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, ElementType, outputRow, "ElementType") ||
        !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, ElementOccurence, outputRow, "ElementOccurence") ||
        !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, loserFound[1], outputRow, "deltaPositiveFound") ||
        !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, deltaSurvived[1], outputRow, "deltaPositive") ||
        !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, lostOnPass[1], outputRow, "lostOnPassPositive") ||
        !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, sLost[1], outputRow, "sLostPositive") ||
        !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, xLost[1], outputRow, "xLostPositive") ||
        !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, yLost[1], outputRow, "yLostPositive") ||
        !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, deltaWhenLost[1], outputRow, "deltaLostPositive") ||
        !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, loserFound[0], outputRow, "deltaNegativeFound") ||
        !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, deltaSurvived[0], outputRow, "deltaNegative") ||
        !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, lostOnPass[0], outputRow, "lostOnPassNegative") ||
        !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, sLost[0], outputRow, "sLostNegative") ||
        !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, xLost[0], outputRow, "xLostNegative") ||
        !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, yLost[0], outputRow, "yLostNegative") ||
        !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, deltaWhenLost[0], outputRow, "deltaLostNegative") ||
        !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, xTuneSurvived[0], outputRow, "nuxLostNegative") ||
        !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, xTuneSurvived[1], outputRow, "nuxLostPositive") ||
        !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, yTuneSurvived[0], outputRow, "nuyLostNegative") ||
        !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, yTuneSurvived[1], outputRow, "nuyLostPositive"))) ||
      (output_mode==1 && 
       (!SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, ElementName, outputRow, "ElementName") ||
        !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, sStart, outputRow, "s") ||
        !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, ElementType, outputRow, "ElementType") ||
        !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, ElementOccurence, outputRow, "ElementOccurence") ||
        !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, direction, outputRow, "direction") ||
        !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, loserFound[0], outputRow, "deltaFound") ||
        !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, deltaSurvived[0], outputRow, "delta") ||
        !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, lostOnPass[0], outputRow, "lostOnPass") ||
        !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, sLost[0], outputRow, "sLost") ||
        !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, xLost[0], outputRow, "xLost") ||
        !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, yLost[0], outputRow, "yLost") ||
        !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, deltaWhenLost[0], outputRow, "deltaLost") ||
        !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, xTuneSurvived[0], outputRow, "nuxLost") ||
        !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, yTuneSurvived[0], outputRow, "nuyLost")))) {
    SDDS_SetError("Problem writing SDDS table (doMomentumApertureSearch)");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
        
#if !SDDS_MPI_IO
  if (!SDDS_WritePage(&SDDSma)) {
#else
  if (!SDDS_MPI_WritePage(&SDDSma)) {
#endif
    SDDS_SetError("Problem writing SDDS table (doMomentumApertureSearch)");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  if (!inhibitFileSync)
    SDDS_DoFSync(&SDDSma);
      

  free_czarray_2d((void**)coord, 1, 7);
  free_czarray_2d((void**)lostOnPass, (output_mode?1:2), (output_mode?2:1)*nElem);
  free_czarray_2d((void**)loserFound, (output_mode?1:2), (output_mode?2:1)*nElem);
  free_czarray_2d((void**)survivorFound, (output_mode?1:2), (output_mode?2:1)*nElem);
  free_czarray_2d((void**)deltaSurvived, (output_mode?1:2), (output_mode?2:1)*nElem);
  free_czarray_2d((void**)xLost, (output_mode?1:2), (output_mode?2:1)*nElem);
  free_czarray_2d((void**)yLost, (output_mode?1:2), (output_mode?2:1)*nElem);
  free_czarray_2d((void**)deltaWhenLost, (output_mode?1:2), (output_mode?2:1)*nElem);
  free_czarray_2d((void**)sLost, (output_mode?1:2), (output_mode?2:1)*nElem);
  free_czarray_2d((void**)xTuneSurvived, (output_mode?1:2), (output_mode?2:1)*nElem);
  free_czarray_2d((void**)yTuneSurvived, (output_mode?1:2), (output_mode?2:1)*nElem);
  free_czarray_2d((void**)turnByTurnCoord, 5, control->n_passes);
  turnByTurnCoord = NULL;
  
  free(sStart);
  free(ElementName);
  free(ElementType);
  free(ElementOccurence);

  beamline->fiducial_flag = fiducial_flag_save;
  if (fiducialize) {
    delete_phase_references();
    reset_special_elements(beamline, RESET_INCLUDE_ALL&~RESET_INCLUDE_RANDOM);
  }  

  return 1;
}

long multiparticleLocalMomentumAcceptance(
                              RUN *run,
                              VARY *control,
                              ERRORVAL *errcon,
                              LINE_LIST *beamline,
                              double *startingCoord
                              )
{    
#if USE_MPI
  double **coord;
  long nTotal, ip, ie, idelta, nLeft, nLost, nElem, nEachProcessor, code;
  double pCentral, delta;
  double *sStart, **deltaSurvived, **deltaLost;
  char **ElementName, **ElementType;
  int32_t *ElementOccurence;
  short **loserFound;
  long n_working_processors = n_processors - 1;
  double **lostParticles;

#ifdef DEBUG
  double round = 0.5;
  long id;
  ELEMENT_LIST *elem0;
  char s[1000];
  char logFile[100];
  FILE *fpd;
  sprintf(logFile, "mmap-%03d.log", myid);
  fpd = fopen(logFile, "w");
  fprintf(fpd, "Started LMA algorithm for processor %d\n", myid);
#endif
  if (myid==0) {
    printf("Started multi-particle  LMA algorithm\n");
    fflush(stdout);
  }
  
  nElem = nElements;
  nDelta = (delta_positive_limit-delta_negative_limit)/delta_step_size+1;
  deltaStep = (delta_positive_limit-delta_negative_limit)/(nDelta-1);
  nTotal = nElem*nDelta;
  if (nTotal%n_working_processors!=0) {
    printf("Warning: The number of working processors (%ld) does not evenly divide into the number of particles (nDelta=%ld, nElem=%ld)\n",
            n_working_processors, nDelta, nElem);
    nEachProcessor =  (nTotal/n_working_processors)+1;
  } else {
    nEachProcessor = nTotal/n_working_processors;
  }
  
#ifdef DEBUG
  fprintf(fpd, "nTotal = %ld, n_working_processors = %ld, nDelta = %ld, deltaStep = %le, nElements = %ld, nEachProcessor = %ld\n",
          nTotal, n_working_processors, nDelta, deltaStep, nElements, nEachProcessor);
#endif
  if (myid==0) {
    printf("nTotal = %ld, n_working_processors = %ld, nDelta = %ld, deltaStep = %le, nElements = %ld, nEachProcessor = %ld\n",
            nTotal, n_working_processors, nDelta, deltaStep, nElements, nEachProcessor);
    fflush(stdout);
  }
  
  /* allocate and initialize array for tracking */
  if (myid==0) {
    coord = (double**)czarray_2d(sizeof(**coord), nEachProcessor*n_working_processors, 7);
    lostParticles = NULL;
  } else {
    coord = (double**)czarray_2d(sizeof(**coord), nEachProcessor, 7);
    lostParticles = (double**)czarray_2d(sizeof(double), nEachProcessor, 8);	 
  }
  
#ifdef DEBUG
  fprintf(fpd, "Arrays allocated\n");
  fflush(fpd);
#endif

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
  if (myid==0) {
    printf("Fiducial particle tracked.\n");
    fflush(stdout);
  }
  
#ifdef DEBUG
  fprintf(fpd, "Fiducialized\n"); fflush(fpd);
#endif
  MPI_Barrier(MPI_COMM_WORLD);
  printf("Fiducalization completed\n");
  fflush(stdout);

  if (myid!=0) {
    for (ip=0; ip<nEachProcessor; ip++) 
      lostParticles[ip][6] = -2;
    for (ip=0; ip<nEachProcessor; ip++) {
      if (startingCoord)
        memcpy(coord[ip], startingCoord, sizeof(**coord)*6);
      else
        memset(coord[ip], 0, sizeof(**coord)*6);
      coord[ip][6] = (myid-1)*nEachProcessor + ip;
      if (coord[ip][6]>=nTotal) {
        /* Don't track more buffer particles than needed */
        coord[ip][6] = -1;
        nEachProcessor = ip+1;
      }
    }
    setTrackingOmniWedgeFunction(momentumOffsetFunctionOmni); 
#ifdef DEBUG
    fprintf(fpd, "Tracking\n");
    fflush(fpd);
#endif
    printf("Tracking\n");
    fflush(stdout);
    nLeft = do_tracking(NULL, coord, nEachProcessor, NULL, beamline, &pCentral, 
                        NULL, NULL, NULL, NULL, run, control->i_step, 
                        FIDUCIAL_BEAM_SEEN+FIRST_BEAM_IS_FIDUCIAL+SILENT_RUNNING+INHIBIT_FILE_OUTPUT, 
                        control->n_passes, 0, NULL, NULL, NULL, lostParticles, NULL);
    nLost = nEachProcessor - nLeft;
#ifdef DEBUG
    fprintf(fpd, "Done tracking nLeft = %ld, nLost = %ld\n", nLeft, nLost);
    fflush(fpd);
#endif
    printf("Done tracking nLeft = %ld, nLost = %ld\n", nLeft, nLost);
    fflush(stdout);
    setTrackingOmniWedgeFunction(NULL);
  } else
    nLeft = 0;

  MPI_Barrier(MPI_COMM_WORLD);

  /* gather lost particle data to master */
  gatherLostParticles(&lostParticles, &nLost, nLeft, n_processors, myid);
  printf("Lost-particle gather done\n"); fflush(stdout);

  if (myid==0) {
    long slot;
    if (!SDDS_StartPage(&SDDSma, nElem) ||
        !SDDS_SetParameters(&SDDSma, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "Step",
                            control->i_step, NULL)) {
      SDDS_SetError("Problem writing SDDS table (doMomentumApertureSearch)");
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
    
    /* allocate arrays for storing data for negative and positive momentum limits for each element */
    deltaLost = (double**)czarray_2d(sizeof(**deltaLost), 2, nElem);
    for (ie=0; ie<nElem; ie++) {
      deltaLost[0][ie] = -DBL_MAX;
      deltaLost[1][ie] = DBL_MAX;
    }
    deltaSurvived = (double**)czarray_2d(sizeof(**deltaSurvived), 2, nElem);
    sStart = (double*)tmalloc(sizeof(*sStart)*nElem);
    ElementName = (char**)tmalloc(sizeof(*ElementName)*nElem);
    ElementType = (char**)tmalloc(sizeof(*ElementType)*nElem);
    ElementOccurence = (int32_t*)tmalloc(sizeof(*ElementOccurence)*nElem);
    loserFound = (short**)czarray_2d(sizeof(**loserFound), 2, nElem);
    
    for (ip=0; ip<nLost; ip++) {
      if (lostParticles[ip][6]<0) {
        if (lostParticles[ip][6]==-2) {
          bombElegant("problem with lost particle accounting!", NULL);
        }
        /* buffer particle, ignore */
        continue;
      }
      idelta = ((long)lostParticles[ip][6])%nDelta;
      ie = (lostParticles[ip][6]-idelta)/nDelta;
      delta = delta_negative_limit + idelta * deltaStep;
      if (delta>=0) {
        slot = 1;
        if (!loserFound[slot][ie] || delta<deltaLost[slot][ie]) {
          loserFound[slot][ie] = 1;
          deltaLost[slot][ie] = delta;
        }
      } else {
        slot = 0;
        if (!loserFound[slot][ie] || delta>deltaLost[slot][ie]) {
          loserFound[slot][ie] = 1;
          deltaLost[slot][ie] = delta;
        }
      }
      sStart[ie] = elementArray[ie]->end_pos;
      ElementName[ie] = elementArray[ie]->name;
      ElementType[ie] = entity_name[elementArray[ie]->type];
      ElementOccurence[ie] = elementArray[ie]->occurence;
    }
    for (ie=0; ie<nElem; ie++) {
      deltaSurvived[0][ie] = deltaLost[0][ie] + deltaStep;
      deltaSurvived[1][ie] = deltaLost[1][ie] - deltaStep;
    }
    if (!SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, ElementName, nElem, "ElementName") ||
        !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, sStart, nElem, "s") ||
        !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, ElementType, nElem, "ElementType") ||
        !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, ElementOccurence, nElem, "ElementOccurence") ||
        !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, deltaSurvived[1], nElem, "deltaPositive") ||
        !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, deltaSurvived[0], nElem, "deltaNegative")) {
      SDDS_SetError("Problem writing SDDS table (doMomentumApertureSearch)");
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
    if (!SDDS_WritePage(&SDDSma)) {
      SDDS_SetError("Problem writing SDDS table (doMomentumApertureSearch)");
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
    if (!inhibitFileSync)
      SDDS_DoFSync(&SDDSma);

    free(ElementType);
    free(ElementName);
    free(ElementOccurence);
    free(elementArray);
    free(sStart);
    free_czarray_2d((void**)deltaLost, 2, nElem);
    free_czarray_2d((void**)deltaSurvived, 2, nElem);
    /* free_czarray_2d((void**)coord, nEachProcessor*n_working_processors, 7); */
  } else {
    free_czarray_2d((void**)coord, nEachProcessor, 7);
  }
  
  MPI_Barrier(MPI_COMM_WORLD);
#else
  bombElegant("This LMA method is only available in Pelegant.", NULL);
#endif
  return 1;
}


long determineTunesFromTrackingData(double *tune, double **turnByTurnCoord, long turns, double delta)
{
  double amplitude[4], frequency[4], phase[4], dummy;

#if defined(DEBUG)
  long i;
  static FILE *fpd = NULL;
  if (!fpd) {
    fpd = fopen("tbt.sdds", "w");
    fprintf(fpd, "SDDS1\n&column name=Pass type=long &end\n");
    fprintf(fpd, "&column name=x type=double units=m &end\n");
    fprintf(fpd, "&column name=xp type=double &end\n");
    fprintf(fpd, "&column name=y type=double units=m &end\n");
    fprintf(fpd, "&column name=yp type=double &end\n");
    fprintf(fpd, "&column name=delta type=double &end\n");
    fprintf(fpd, "&parameter name=delta0 type=double &end\n");
    fprintf(fpd, "&data mode=ascii &end\n");
  }
  fprintf(fpd, "%ld\n", turns);
  for (i=0; i<turns; i++) {
    fprintf(fpd, "%ld %e %e %e %e %e\n", i,
            turnByTurnCoord[0][i],
            turnByTurnCoord[1][i],
            turnByTurnCoord[2][i],
            turnByTurnCoord[3][i],
            turnByTurnCoord[4][i]);
  }
#endif

  if (PerformNAFF(&frequency[0], &amplitude[0], &phase[0], 
		  &dummy, 0.0, 1.0, turnByTurnCoord[0], turns, 
		  NAFF_MAX_FREQUENCIES|NAFF_FREQ_CYCLE_LIMIT|NAFF_FREQ_ACCURACY_LIMIT,
		  0.0, 1, 200, 1e-12, 0, 0)!=1) {
    printf("Warning: NAFF failed for tune analysis from tracking (x).\n");
    return 0;
  }

  if (PerformNAFF(&frequency[1], &amplitude[1], &phase[1], 
		  &dummy, 0.0, 1.0, turnByTurnCoord[1], turns, 
		  NAFF_MAX_FREQUENCIES|NAFF_FREQ_CYCLE_LIMIT|NAFF_FREQ_ACCURACY_LIMIT,
		  0.0, 1, 200, 1e-12,0, 0)!=1) {
    printf("Warning: NAFF failed for tune analysis from tracking (xp).\n");
    return 0;
  }

  if (PerformNAFF(&frequency[2], &amplitude[2], &phase[2], 
		  &dummy, 0.0, 1.0, turnByTurnCoord[2], turns,
		  NAFF_MAX_FREQUENCIES|NAFF_FREQ_CYCLE_LIMIT|NAFF_FREQ_ACCURACY_LIMIT,
		  0.0, 1, 200, 1e-12, 0, 0)!=1) {
    printf("Warning: NAFF failed for tune analysis from tracking (y).\n");
    return 0;
  }

  if (PerformNAFF(&frequency[3], &amplitude[3], &phase[3], 
		  &dummy, 0.0, 1.0, turnByTurnCoord[3], turns,
		  NAFF_MAX_FREQUENCIES|NAFF_FREQ_CYCLE_LIMIT|NAFF_FREQ_ACCURACY_LIMIT,
		  0.0, 1, 200, 1e-12, 0, 0)!=1) {
    printf("Warning: NAFF failed for tune analysis from tracking (yp).\n");
    return 0;
  }

  tune[0] = adjustTuneHalfPlane(frequency[0], phase[0], phase[1]);
  tune[1] = adjustTuneHalfPlane(frequency[2], phase[2], phase[3]);
  return 1;
}

#if USE_MPI
void gatherLostParticles(double ***lostParticles, long *nLost, long nSurvived, long n_processors, int myid) 
{
  long work_processors = n_processors-1;
  int root = 0, i, nItems, displs;
  int my_nLost, *nLostCounts, nLost_total; 
 
  MPI_Status status;
  nLostCounts = malloc(sizeof(int) * n_processors);

  if (myid==0) {
    my_nLost = 0;
  }
  else {
    my_nLost = *nLost;
  }
  
  MPI_Gather(&my_nLost, 1, MPI_INT, nLostCounts, 1, MPI_INT, root, MPI_COMM_WORLD);
  MPI_Reduce(&my_nLost, &nLost_total, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

  if (myid==0) {
    /* set up the displacement array and the number of elements that are received from each processor */ 
    nLostCounts[0] = 0;
    displs = 0;
    *lostParticles = (double**)resize_czarray_2d((void**)(*lostParticles), sizeof(double), nLost_total, COORDINATES_PER_PARTICLE+1);
    for (i=1; i<=work_processors; i++) {
      /* gather information for lost particles */
      displs = displs+nLostCounts[i-1];
      nItems = nLostCounts[i]*(COORDINATES_PER_PARTICLE+1);
      MPI_Recv (&(*lostParticles)[displs][0], nItems, MPI_DOUBLE, i, 102, MPI_COMM_WORLD, &status);
    } 
  } else {
    /* send information for lost particles */
    MPI_Send (&(*lostParticles)[nSurvived][0], my_nLost*(COORDINATES_PER_PARTICLE+1), MPI_DOUBLE, root, 102, MPI_COMM_WORLD);  
  }

  if (myid==0) {
    *nLost = nLost_total;
  }
}
#endif
