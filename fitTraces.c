/* Copyright 1997 by Michael Borland and Argonne National Laboratory,
 * all rights reserved.
 */
/*
 * $Log: not supported by cvs2svn $
 */
/* file: fitTraces.c 
 *
 * Michael Borland, 1997
 */
#include "mdb.h"
#include "track.h"
#include "fitTraces.h"

typedef struct {
  long BPMs, traces;        /* number of BPMs and number of traces */
  char **BPMName;           /* BPM names */
  double **x, **y;          /* measured coordinates at BPMs */
  ELEMENT_LIST **element;   /* pointer to the element structure */
  long *latticeIndex;       /* index of BPM in the lattice */
  double **startingCoord;   /* starting coordinates of each trace, varied in fitting */
} FIT_TRACE_DATA;

typedef struct {
  long parameters;
  char **elementName, **parameterName;
  double *delta;
  double **paramData; /* will point to the actual location used to store the parameter value in
                         the element structure */
  ELEMENT_LIST **target;
  double *definedValue;
} FIT_TRACE_PARAMETERS ;

typedef struct {
  SDDS_TABLE SDDStable;
  long **traceDataIndex, *paramDataIndex;
} FIT_OUTPUT_DATA;


FIT_TRACE_PARAMETERS *fit_traces_readFitParametersFile(char *dataFile, LINE_LIST *beamline);
FIT_TRACE_DATA *fit_traces_readTraceDataFile(char *dataFile, LINE_LIST *beamline);
void find_trajectory_bpm_readouts(double *xReadout, double *yReadout, 
                                  ELEMENT_LIST **bpmElement, long *BPMIndex, long BPMs, 
                                  LINE_LIST *beamline, RUN *run,
                                  TRAJECTORY *trajBuffer, double *startingCoordinate, double momentum);
void fit_traces_findDerivatives(FIT_TRACE_DATA *traceData, FIT_TRACE_PARAMETERS *fitParam,
                                MATRIX *D, LINE_LIST *beamline, RUN *run);
double fit_trace_findReadbackErrors(MATRIX *readbackError, FIT_TRACE_DATA *traceData,
                                    LINE_LIST *beamline, RUN *run);
FIT_OUTPUT_DATA *fit_trace_setUpOutputFile(char *filename,
                                           FIT_TRACE_PARAMETERS *fitParam,
                                           FIT_TRACE_DATA *traceData, long iterations);
void fit_trace_saveParamValues(double *buffer, FIT_TRACE_DATA *traceData,
                                   FIT_TRACE_PARAMETERS *fitParam);
void fit_trace_restoreParamValues(double *buffer, FIT_TRACE_DATA *traceData,
                                  FIT_TRACE_PARAMETERS *fitParam, RUN *run);
double fit_trace_takeStep(MATRIX *D, MATRIX *Dt, MATRIX *DtD, MATRIX *DtDInv, MATRIX *DtDInvDt, 
                          MATRIX *readbackVector, MATRIX *paramVector,FIT_TRACE_DATA *traceData, 
                          FIT_TRACE_PARAMETERS *fitParam, LINE_LIST *beamline, RUN *run,
                          double convergenceFactor);
void fit_trace_setRowValues(FIT_OUTPUT_DATA *outputData, long iteration, 
                            double rmsError, double lastRmsError, double pass,
                            double convergenceFactor, FIT_TRACE_DATA *traceData, 
                            FIT_TRACE_PARAMETERS *fitParam);

void do_fit_trace_data(NAMELIST_TEXT *nltext, RUN *run, LINE_LIST *beamline)
{
  double result;
  long i, j, iteration, subIteration, parameters, readbacks;
  long iCoord, iBPM, iTrace, iUserParam, offset;
  double rmsError, lastRmsError;
  FIT_TRACE_PARAMETERS *fitParam;
  FIT_TRACE_DATA *traceData;
  MATRIX *D, *Dt, *DtD, *DtDInv, *DtDInvDt, *readbackVector, *paramVector;
  double *lastParameterValues, *startParameterValues;
  FIT_OUTPUT_DATA *outputData;
  double rmsErrorDelta[3];
  long iMin, iMax, pass=0;
  
  /* process namelist text */
  process_namelist(&fit_traces, nltext);
  print_namelist(stdout, &fit_traces);

  if (!trace_data_file || !fexists(trace_data_file)) {
    fprintf(stderr, "fit_traces: trace_data_file file not given or not found\n");
    exit(1);
  }
  if (!fit_parameters_file || !fexists(fit_parameters_file)) {
    fprintf(stderr, "fit_traces: fit_parameters_file file not given or not found\n");
    exit(1);
  }
  if (iterations<1) {
    fprintf(stderr, "fit_traces: iterations<1\n");
    exit(1);
  }
  if (!output) {
    fprintf(stderr, "fit_traces: output file not given\n");
    exit(1);
  }
  
  /* trace_data_file file contains column data giving
   * x readings, y readings, and BPM name.
   * Each page is a separate trace.
   */
  traceData = fit_traces_readTraceDataFile(trace_data_file, beamline);
  if (!traceData->traces) {
    fprintf(stderr, "fit_traces: no trace data in %s\n", trace_data_file);
    exit(1);
  }
  fprintf(stderr, "%ld traces with %ld BPMs:\n", traceData->traces, traceData->BPMs);
  for (iTrace=0; iTrace<traceData->traces; iTrace++) {
    fprintf(stderr, "Trace %ld: ", iTrace);
    for (iBPM=0; iBPM<traceData->BPMs; iBPM++) {
      fprintf(stderr, "(%10.3le, %10.3le)  ",
              traceData->x[iTrace][iBPM],
              traceData->y[iTrace][iBPM]);
    }
    fprintf(stderr, "\n");
  }
  
  /* fit_parameters_file file contains column data giving ElementName, ElementParameter
   */
  fitParam = fit_traces_readFitParametersFile(fit_parameters_file, beamline);

  outputData = fit_trace_setUpOutputFile(output, fitParam, traceData, iterations);
  
  parameters = 4*traceData->traces+fitParam->parameters;
  readbacks  = 2*traceData->BPMs*traceData->traces;
  lastParameterValues = tmalloc(sizeof(*lastParameterValues)*parameters);
  startParameterValues = tmalloc(sizeof(*startParameterValues)*parameters);
  
  /* check that there are enough traces for the given number of parameters */
  if (2*traceData->BPMs*traceData->traces<(4*traceData->traces+fitParam->parameters)) {
    fprintf(stderr, "fit_traces: too few traces for given number of fit parameters.\n");
    exit(1);
  }
  
  /* allocate matrices to solve the problem:
   * readbackVector = D*paramVector 
   * where D is the matrix of derivatives of readbacks wrt parameters 
   * least squares solution is
   * paramVector = Inv(Tr(D)*D)*Tr(D)*readbackVector
   * D is organized as follows:
   *   D(i, j) = derivative of ith readout w.r.t. jth parameter where
   *   i = 2*t*B + 2*b + c with
   *       t = trace index [0, T-1], b = BPM index [0, B-1], c = x/y coordinate index [0, 1]
   *   j = p+4*t on [0, 4T-1] with
   *         p = phase-space coord index, [0, 3]
   *         t = trace index [0, T-1]
   *     = u+4T on [4T, 4T+U-1] with
   *         u = user-parameter index, [0, U-1]
   * readbackVector is organized as follows:
   *   rV(i,0) = ith readout, where i is as above
   * parameterVector is organized as follows:
   *   pV(j,0) = jth parameter, where j is as above
   */
  m_alloc(&D, readbacks, parameters);
  m_alloc(&Dt, parameters, readbacks);
  m_alloc(&DtD, parameters, parameters);
  m_alloc(&DtDInv, parameters, parameters);
  m_alloc(&DtDInvDt, parameters, readbacks);
  m_alloc(&readbackVector, readbacks, 1);    /* really the vector of readback errors */
  m_alloc(&paramVector, parameters, 1);      /* really the vector of parameter deltas */

  lastRmsError = 0;
  for (iteration=0; iteration<iterations; iteration++) {
    fprintf(stderr, "Doing iteration %ld  ", iteration);
    subIteration = sub_iterations;

    fit_trace_setRowValues(outputData, iteration, rmsError, lastRmsError, pass, 
                           convergenceFactor, traceData, fitParam);

    lastRmsError = fit_trace_findReadbackErrors(readbackVector, traceData, beamline, run);

    fit_trace_saveParamValues(startParameterValues, traceData, fitParam);
    rmsErrorDelta[0] 
      = fit_trace_takeStep(D, Dt, DtD, DtDInv, DtDInvDt, readbackVector, paramVector,
                           traceData, fitParam, beamline, run,
                           convergenceFactor)
        - lastRmsError;
    fit_trace_restoreParamValues(startParameterValues, traceData, fitParam, run);
    fit_trace_findReadbackErrors(readbackVector, traceData, beamline, run);

    rmsErrorDelta[1] 
      = fit_trace_takeStep(D, Dt, DtD, DtDInv, DtDInvDt, readbackVector, paramVector,
                           traceData, fitParam, beamline, run,
                           convergenceFactor+convergenceFactorDelta) 
        - lastRmsError;
    fit_trace_restoreParamValues(startParameterValues, traceData, fitParam, run);
    fit_trace_findReadbackErrors(readbackVector, traceData, beamline, run);

    rmsErrorDelta[2] 
      = fit_trace_takeStep(D, Dt, DtD, DtDInv, DtDInvDt, readbackVector, paramVector,
                           traceData, fitParam, beamline, run,
                           convergenceFactor-convergenceFactorDelta) 
        - lastRmsError;
    fit_trace_restoreParamValues(startParameterValues, traceData, fitParam, run);
    fit_trace_findReadbackErrors(readbackVector, traceData, beamline, run);

    index_min_max(&iMin, &iMax, rmsErrorDelta, 3);
    if (iMin==1) {
      convergenceFactor += convergenceFactorDelta;
    } else if (iMin==2) {
      convergenceFactor -= convergenceFactorDelta;
    }
    fprintf(stderr, " C=%lf ", convergenceFactor);
    
    do {
      /* fit_trace_saveParamValues(lastParameterValues, traceData, fitParam); */
      lastRmsError = rmsError;
      rmsError = fit_trace_takeStep(D, Dt, DtD, DtDInv, DtDInvDt, readbackVector, paramVector,
                                    traceData, fitParam, beamline, run,
                                    convergenceFactor);
      if (rmsError<target)
        break;
      if (lastRmsError<rmsError) {
        convergenceFactor -= convergenceFactorBackoff;
        /* 
        fit_trace_restoreParamValues(lastParameterValues, traceData, fitParam, run);
        fit_trace_findReadbackErrors(readbackVector, traceData, beamline, run);
        */
      }
      if (convergenceFactor<convergenceFactorMin)
        convergenceFactor = convergenceFactorMin;
      else if (convergenceFactor>convergenceFactorMax)
        convergenceFactor = convergenceFactorMax;
      pass++;
    } while (--subIteration > 0);
    
    rmsError = fit_trace_findReadbackErrors(readbackVector, traceData, beamline, run);
    fprintf(stderr, "rms error is %le  C=%le\n", rmsError, convergenceFactor);
    if (!lastRmsError)
      lastRmsError = rmsError;

    if (rmsError<target)
      break;
  }
  fit_trace_setRowValues(outputData, iteration, rmsError, lastRmsError, pass, 
                         convergenceFactor, traceData, fitParam);
  if (!SDDS_Terminate(&outputData->SDDStable)) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
}

double fit_trace_findReadbackErrors
  (
   MATRIX *readbackError,
   FIT_TRACE_DATA *traceData,
   LINE_LIST *beamline,
   RUN *run
   )
{
  long iTrace, iBPM, iCoord, row;
  static double *x=NULL, *y=NULL;
  static long lastBPMs = 0, lastElements = 0;
  double startingCoord[4], p, sum;
  static TRAJECTORY *trajectory = NULL;
  
  if (!lastBPMs || lastBPMs!=traceData->BPMs) {
    lastBPMs = traceData->BPMs;
    if (!(x=SDDS_Realloc(x, sizeof(*x)*lastBPMs)) ||
        !(y=SDDS_Realloc(y, sizeof(*y)*lastBPMs))) {
      fprintf(stderr, "Memory allocation failure in fit_trace_findReadbackErrors (1)\n");
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      exit(1);
    }
  }
  if (!lastElements || lastElements!=beamline->n_elems) {
    if (!(trajectory = SDDS_Realloc(trajectory, sizeof(*trajectory)*beamline->n_elems))) {
      fprintf(stderr, "Memory allocation failure in fit_trace_findReadbackErrors (2)\n");
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      exit(1);
    }
    lastElements = beamline->n_elems;
  }

  /* track each trace to predict position */
  p = sqrt(sqr(run->ideal_gamma)+1);
  sum = 0;
  for (iTrace=0; iTrace<traceData->traces; iTrace++) {
    for (iCoord=0; iCoord<4; iCoord++)
      startingCoord[iCoord] = traceData->startingCoord[iTrace][iCoord];
    find_trajectory_bpm_readouts(x, y, traceData->element,
                                 traceData->latticeIndex, traceData->BPMs,
                                 beamline, run, trajectory, startingCoord, p);
    for (iBPM=0; iBPM<traceData->BPMs; iBPM++) {
      row = 2*iTrace*traceData->BPMs + 2*iBPM ;
      readbackError->a[row  ][0] = x[iBPM]-traceData->x[iTrace][iBPM];
      readbackError->a[row+1][0] = y[iBPM]-traceData->y[iTrace][iBPM];
      sum += sqr(readbackError->a[row  ][0])+sqr(readbackError->a[row+1][0]);
    }
  }
  return sqrt(sum);
}
 
void fit_traces_findDerivatives
  (
   FIT_TRACE_DATA *traceData,
   FIT_TRACE_PARAMETERS *fitParam,
   MATRIX *D,
   LINE_LIST *beamline,
   RUN *run
   )
{
  long iTrace, iBPM, iUserParam, iCoord, row, column;
  static double *x0=NULL, *y0=NULL, *x=NULL, *y=NULL;
  static long lastBPMs = 0, lastElements = 0;
  double startingCoord[4], refValue, p;
  static TRAJECTORY *trajectory = NULL;
  
  if (!lastBPMs || lastBPMs!=traceData->BPMs) {
    lastBPMs = traceData->BPMs;
    if (!(x0=SDDS_Realloc(x0, sizeof(*x0)*lastBPMs)) ||
        !(y0=SDDS_Realloc(y0, sizeof(*y0)*lastBPMs)) ||
        !(x=SDDS_Realloc(x, sizeof(*x)*lastBPMs)) ||
        !(y=SDDS_Realloc(y, sizeof(*y)*lastBPMs))) {
      fprintf(stderr, "Memory allocation failure in fit_traces_findDerivatives\n");
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      exit(1);
    }
  }
  if (!lastElements || lastElements!=beamline->n_elems) {
    if (!(trajectory = SDDS_Realloc(trajectory, sizeof(*trajectory)*beamline->n_elems))) {
      fprintf(stderr, "Memory allocation failure in fit_traces_findDerivatives\n");
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      exit(1);
    }
    lastElements = beamline->n_elems;
  }

  for (iCoord=0; iCoord<D->n; iCoord++) {
    for (iBPM=0; iBPM<D->m; iBPM++) {
      D->a[iCoord][iBPM] = 0;
    }
  }

  /* two types of parameters: 
     1. initial trajectory coordinates 
     2. element parameters declared by user
     */
  
  /* compute and store derivative data for initial trajectory coordintes */
  p = sqrt(sqr(run->ideal_gamma)+1);
  for (iTrace=0; iTrace<traceData->traces; iTrace++) {
    /* find reference traj for derivatives */
    for (iCoord=0; iCoord<4; iCoord++)
      startingCoord[iCoord] = traceData->startingCoord[iTrace][iCoord];
    find_trajectory_bpm_readouts(x0, y0, traceData->element,
                                 traceData->latticeIndex, traceData->BPMs,
                                 beamline, run, trajectory, startingCoord, p);
    for (iCoord=0; iCoord<4; iCoord++) {
      refValue = startingCoord[iCoord];
      startingCoord[iCoord] += 1e-6;
      find_trajectory_bpm_readouts(x, y, traceData->element,
                                   traceData->latticeIndex, traceData->BPMs,
                                   beamline, run, trajectory, startingCoord, p);
      startingCoord[iCoord] = refValue;
      column = iCoord+4*iTrace;
      for (iBPM=0; iBPM<traceData->BPMs; iBPM++) {
        row = 2*iTrace*traceData->BPMs+2*iBPM;
        D->a[row  ][column] = (x[iBPM]-x0[iBPM])/1e-6;
        D->a[row+1][column] = (y[iBPM]-y0[iBPM])/1e-6;
      }
    }
  }
  
  for (iTrace=0; iTrace<traceData->traces; iTrace++) {
    for (iCoord=0; iCoord<4; iCoord++)
      startingCoord[iCoord] = traceData->startingCoord[iTrace][iCoord];
    find_trajectory_bpm_readouts(x0, y0, traceData->element,
                                 traceData->latticeIndex, traceData->BPMs,
                                 beamline, run, trajectory, startingCoord, p);
    for (iUserParam=0; iUserParam<fitParam->parameters; iUserParam++) {
      refValue = *(fitParam->paramData[iUserParam]);
      *(fitParam->paramData[iUserParam]) += fitParam->delta[iUserParam];
      if (fitParam->target[iUserParam]->matrix) {
        free_matrices(fitParam->target[iUserParam]->matrix);
        free(fitParam->target[iUserParam]->matrix);
      }
      compute_matrix(fitParam->target[iUserParam], run, NULL);
      find_trajectory_bpm_readouts(x, y, traceData->element,
                                   traceData->latticeIndex, traceData->BPMs,
                                   beamline, run, trajectory, startingCoord, p);
      *(fitParam->paramData[iUserParam]) = refValue;
      column = 4*traceData->traces + iUserParam;
      for (iBPM=0; iBPM<traceData->BPMs; iBPM++) {
        row = 2*iTrace*traceData->BPMs+2*iBPM;
        D->a[row  ][column] = (x[iBPM]-x0[iBPM])/1e-6;
        D->a[row+1][column] = (y[iBPM]-y0[iBPM])/1e-6;
      }
    }
  }

}



void find_trajectory_bpm_readouts
  (
   double *xReadout, double *yReadout,  /* arrays in which to return x and y readouts */
   ELEMENT_LIST **bpmElement,           /* pointers to BPM elements in the beamline */
   long *latticeIndex,                  /* index in the lattice of BPMs */
   long BPMs, 
   LINE_LIST *beamline, RUN *run,
   TRAJECTORY *trajBuffer, double *startingCoordinate, double momentum)
{
  static double **particle = NULL;
  long tracking_flags = TEST_PARTICLES, nPart = 1;
  long iBPM, i;
  
  if (!particle) {
    particle = (double**)zarray_2d(sizeof(**particle), 1, 7);
    particle[0][4] = particle[0][5] = 0;
  }
  for (i=0; i<4; i++)
    particle[0][i] = startingCoordinate[i];
  
  if (!do_tracking(particle, &nPart, NULL, beamline, &momentum,
                   (double**)NULL, (BEAM_SUMS**)NULL, (long*)NULL,
                   trajBuffer, run, 0, tracking_flags, 1)) {
    fprintf(stderr, "Error tracking particle to find trajectory at BPMs.\n");
    exit(1);
  }
  for (iBPM=0; iBPM<BPMs; iBPM++) {
    xReadout[iBPM] 
      = computeMonitorReading(bpmElement[iBPM], 0, 
                              trajBuffer[latticeIndex[iBPM]].centroid[0],
                              trajBuffer[latticeIndex[iBPM]].centroid[2]);
    yReadout[iBPM] 
      = computeMonitorReading(bpmElement[iBPM], 1, 
                              trajBuffer[latticeIndex[iBPM]].centroid[0],
                              trajBuffer[latticeIndex[iBPM]].centroid[2]);
/*
    fprintf(stderr, "%ld: (%le, %le) -> (%le, %le)\n", iBPM,
            trajBuffer[latticeIndex[iBPM]].centroid[0],
            trajBuffer[latticeIndex[iBPM]].centroid[2],
            xReadout[iBPM], yReadout[iBPM]);
*/
  }

  
}

FIT_TRACE_PARAMETERS *fit_traces_readFitParametersFile
  (
   char *dataFile,
   LINE_LIST *beamline
   )
{
  SDDS_TABLE SDDSin;
  long i;
  FIT_TRACE_PARAMETERS *ftp;
  long parameterIndex, elementType;
  
  if (!SDDS_InitializeInput(&SDDSin, dataFile)) {
    fprintf(stderr, "Error: couldn't read file %s\n", dataFile);
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    exit(1);
  }
  if (SDDS_CHECK_OKAY!=SDDS_CheckColumn(&SDDSin, "ElementName", NULL, SDDS_STRING, stderr) ||
      SDDS_CHECK_OKAY!=SDDS_CheckColumn(&SDDSin, "ElementParameter", NULL, SDDS_STRING, stderr) ||
      SDDS_CHECK_OKAY!=SDDS_CheckColumn(&SDDSin, "Delta", NULL, SDDS_DOUBLE, stderr)) {
    fprintf(stderr, "Problem with column(s) in file %s\n", dataFile);
    exit(1);
  }
  if (SDDS_ReadPage(&SDDSin)<=0) {
    fprintf(stderr, "Problem reading data from file %s\n", dataFile);
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    exit(1);
  }
  if (!(ftp=malloc(sizeof(*ftp)))) {
    fprintf(stderr, "Error: memory allocation failure (fit_traces_readFitParametersFile)\n");
    exit(1);
  }
  if (!(ftp->parameters=SDDS_CountRowsOfInterest(&SDDSin))) {
    fprintf(stderr, "No data in file %s\n", dataFile);
    exit(1);
  }
  if (!(ftp->paramData=malloc(sizeof(*ftp->paramData)*ftp->parameters))) {
    fprintf(stderr, "Error: memory allocation failure (fit_traces_readFitParametersFile)\n");
    exit(1);
  }
  if (!(ftp->elementName = SDDS_GetColumn(&SDDSin, "ElementName")) ||
      !(ftp->parameterName = SDDS_GetColumn(&SDDSin, "ElementParameter")) || 
      !(ftp->delta = SDDS_GetColumn(&SDDSin, "Delta"))) {
    fprintf(stderr, "Problem reading data from file %s\n", dataFile);
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    exit(1);
  }
  if (!(ftp->target = malloc(sizeof(*ftp->target)*ftp->parameters)) ||
      !(ftp->definedValue = malloc(sizeof(*ftp->definedValue)*ftp->parameters))) {
    fprintf(stderr, "Error: memory allocation problem reading parameters file\n");
    exit(1);
  }
  for (i=0; i<ftp->parameters; i++) {
    if (!(ftp->target[i]=find_element(ftp->elementName[i], NULL, &(beamline->elem)))) {
      fprintf(stderr, "Error: element %s not found in beamline\n", ftp->elementName[i]);
    }
    elementType = ftp->target[i]->type;
    if ((parameterIndex=confirm_parameter(ftp->parameterName[i], elementType))<0) {
      fprintf(stderr, "Error: element %s does not have a parameter called %s\n", 
              ftp->parameterName[i]);
      exit(1);
    }
    if (entity_description[elementType].parameter[parameterIndex].type!=IS_DOUBLE) {
      fprintf(stderr, "Error: parameter %s of element %s is not a double value\n",
              ftp->parameterName[i], ftp->elementName[i]);
      exit(1);
    }
    ftp->paramData[i]
      = ((double*)(ftp->target[i]->p_elem +
                   entity_description[elementType].parameter[parameterIndex].offset));
    ftp->definedValue[i] = *(ftp->paramData[i]);
  }
  
  if (SDDS_ReadPage(&SDDSin)>1)
    fprintf(stderr, "Warning: file %s has multiple pages---only the first is used.\n", dataFile);
  SDDS_Terminate(&SDDSin);
  fprintf(stderr, "%ld fit parameters defined.\n", ftp->parameters);
  
  return ftp;
}


FIT_TRACE_DATA *fit_traces_readTraceDataFile
  (
   char *dataFile, 
   LINE_LIST *beamline
   ) 
{
  SDDS_TABLE SDDSin;
  long maxTraces, iTrace, indexFirst, iTraceFirst;
  long iBPMFirst, iBPM;
  FIT_TRACE_DATA *trace;
  char **BPMName;

  if (!SDDS_InitializeInput(&SDDSin, dataFile)) {
    fprintf(stderr, "Error: couldn't read file %s\n", dataFile);
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    exit(1);
  }
  if (SDDS_CHECK_OKAY!=SDDS_CheckColumn(&SDDSin, "x", "m", SDDS_ANY_NUMERIC_TYPE, stderr) ||
      SDDS_CHECK_OKAY!=SDDS_CheckColumn(&SDDSin, "y", "m", SDDS_ANY_NUMERIC_TYPE, stderr) ||
      SDDS_CHECK_OKAY!=SDDS_CheckColumn(&SDDSin, "BPMName", NULL, SDDS_STRING, stderr)) {
    fprintf(stderr, "Problem with column(s) x, y, or BPMName in file %s\n",
            dataFile);
    exit(1);
  }
  
  if (!(trace = malloc(sizeof(*trace)))) {
    fprintf(stderr, "Error trying to allocate space for traces.\n");
    exit(1);
  }
  maxTraces = 0;
  iTrace = 0;
  trace->x = trace->y = NULL;
  while (SDDS_ReadPage(&SDDSin)>0) {
    if (maxTraces<=iTrace) {
      if (!(trace->x = SDDS_Realloc(trace->x, (maxTraces+10)*sizeof(*trace->x))) ||
          !(trace->y = SDDS_Realloc(trace->y, (maxTraces+10)*sizeof(*trace->y))) ) {
        fprintf(stderr, "Error trying to allocate space for more traces.\n");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
        exit(1);        
      }
      maxTraces += 10;
    }
    if (!iTrace) {
      if (!(trace->BPMs = SDDS_CountRowsOfInterest(&SDDSin))) {
        fprintf(stderr, "No traces on first page of trace file\n");
        exit(1);
      }
    } else {
      if (trace->BPMs != SDDS_CountRowsOfInterest(&SDDSin)) {
        fprintf(stderr, "Fit traces have different numbers of data points.");
        exit(1);
      }
    }
    if (!(trace->x[iTrace]=SDDS_GetColumnInDoubles(&SDDSin, "x")) || 
        !(trace->y[iTrace]=SDDS_GetColumnInDoubles(&SDDSin, "y"))) {
      fprintf(stderr, "Error trying to read x or y values for trace\n");
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      exit(1);
    }
    if (!(BPMName=SDDS_GetColumn(&SDDSin, "BPMName"))) {
      fprintf(stderr, "Error trying to read BPM names for trace\n");
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      exit(1);
    }
    if (iTrace==0) {
      trace->BPMName = BPMName;
    } else {
      /* require the same BPM names on every page */
      for (iBPM=0; iBPM<trace->BPMs; iBPM++) {
        if (strcmp(trace->BPMName[iBPM], BPMName[iBPM])) {
          fprintf(stderr, "Fit traces have mismatched BPM names---all pages must have the names in the same order.\n");
          exit(1);
        }
      }
    }
    iTrace++;
  }

  trace->traces = iTrace;
  trace->startingCoord = (double**)zarray_2d(sizeof(**trace->startingCoord), trace->traces, 4);
  
  indexFirst = LONG_MAX;
  iBPMFirst = 0;
  if (!(trace->latticeIndex=malloc(sizeof(*trace->latticeIndex)*trace->BPMs)) || \
      !(trace->element=malloc(sizeof(*trace->element)*trace->BPMs))) {
    fprintf(stderr, "Memory allocation failure storing trace data.\n");
    exit(1);
  }
  for (iBPM=0; iBPM<trace->BPMs; iBPM++) {
      if (!(trace->element[iBPM]=find_element_index
            (trace->BPMName[iBPM], NULL, &(beamline->elem), trace->latticeIndex+iBPM))
          || trace->element[iBPM]->type!=T_MONI) {
        fprintf(stderr, "Element %s not found or not of type MONI\n", trace->BPMName[iBPM]);
        exit(1);
      }
      if (trace->latticeIndex[iBPM]<indexFirst) {
        indexFirst = trace->latticeIndex[iBPM];
        iBPMFirst = iBPM;
      }
    }

  for (iTrace=0; iTrace<trace->traces; iTrace++) {
    trace->startingCoord[iTrace][0] = trace->x[iTrace][iBPMFirst];
    trace->startingCoord[iTrace][2] = trace->y[iTrace][iBPMFirst];
    trace->startingCoord[iTrace][1] =
      trace->startingCoord[iTrace][3] = 0;
  }
  
  fprintf(stdout, "%ld traces found with %ld BPMs.\n", trace->traces, trace->BPMs);
  fflush(stdout);
  return trace;
}

FIT_OUTPUT_DATA *fit_trace_setUpOutputFile(char *filename,
                                           FIT_TRACE_PARAMETERS *fitParam,
                                           FIT_TRACE_DATA *traceData, 
                                           long iterations)
{
  long iTrace, iCoord, iUserParam;
  FIT_OUTPUT_DATA *outputData;
  char name[128];
  char *coordName[4] = {
    "x", "xp", "y", "yp",
  };
  char *coordUnits[4] = {
    "m", "rad", "m", "rad", 
  };
  
  if (!(outputData = malloc(sizeof(*outputData))) ||
      !(outputData->traceDataIndex 
        = (long**)zarray_2d(sizeof(**outputData->traceDataIndex), 
                            traceData->traces, 4)) ||
      !(outputData->paramDataIndex
        = malloc(sizeof(*outputData->paramDataIndex)*fitParam->parameters))) {
    fprintf(stderr, "memory allocation failure in fit_trace_setUpOutputFile");
    exit(1);
  }
  
  if (!SDDS_InitializeOutput(&outputData->SDDStable, SDDS_BINARY, 0, NULL, NULL,
                             filename) ||
      0>SDDS_DefineColumn(&outputData->SDDStable, "Iteration", NULL, NULL, 
                          "Major iteration number",
                          NULL, SDDS_LONG, 0) ||
      0>SDDS_DefineColumn(&outputData->SDDStable, "Passes", NULL, NULL, 
                          "Total number of iteration passes", 
                          NULL, SDDS_LONG, 0) ||
      0>SDDS_DefineColumn(&outputData->SDDStable, "RMSError", NULL, "m", 
                          "RMS error of predicted BPM readouts",
                          NULL, SDDS_DOUBLE, 0) ||
      0>SDDS_DefineColumn(&outputData->SDDStable, "RMSErrorImprovement", NULL, NULL, 
                          "Fractional change in RMS error",
                          NULL, SDDS_DOUBLE, 0) || 
      0>SDDS_DefineColumn(&outputData->SDDStable, "ConvergenceFactor", NULL, NULL, 
                          NULL, NULL, SDDS_DOUBLE, 0)) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  for (iTrace=0; iTrace<traceData->traces; iTrace++) {
    for (iCoord=0; iCoord<4; iCoord++) {
      sprintf(name, "Trace%ldInitial%s", iTrace, coordName[iCoord]);
      if (0>(outputData->traceDataIndex[iTrace][iCoord]=
             SDDS_DefineColumn(&outputData->SDDStable, name, NULL, coordUnits[iCoord],
                               NULL, NULL, SDDS_DOUBLE, 0))) {
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      }
    }
  }
  for (iUserParam=0; iUserParam<fitParam->parameters; iUserParam++) {
    sprintf(name, "%s.%s", fitParam->elementName[iUserParam], 
            fitParam->parameterName[iUserParam]);
    if (0>(outputData->paramDataIndex[iUserParam]=
           SDDS_DefineColumn(&outputData->SDDStable, name, NULL, NULL,
                             NULL, NULL, SDDS_DOUBLE, 0))) {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
  }
  if (!SDDS_WriteLayout(&outputData->SDDStable) 
      || !SDDS_StartPage(&outputData->SDDStable, iterations+1)) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  return outputData;
}


void fit_trace_saveParamValues
  (
   double *buffer, 
   FIT_TRACE_DATA *traceData,
   FIT_TRACE_PARAMETERS *fitParam
   )
{
  long iTrace, iCoord, iUserParam, offset;
  for (iTrace=0; iTrace<traceData->traces; iTrace++) 
    for (iCoord=0; iCoord<4; iCoord++) 
      buffer[iTrace*4+iCoord] = traceData->startingCoord[iTrace][iCoord];
  offset = traceData->traces*4;
  for (iUserParam=0; iUserParam<fitParam->parameters; iUserParam++) 
    buffer[iUserParam+offset] = *(fitParam->paramData[iUserParam]);
}

void fit_trace_restoreParamValues
  (
   double *buffer,
   FIT_TRACE_DATA *traceData,
   FIT_TRACE_PARAMETERS *fitParam,
   RUN *run
   )
{
  long iTrace, iCoord, iUserParam, offset;
  
  /* assert new trajectory starting values */
  for (iTrace=0; iTrace<traceData->traces; iTrace++) {
    for (iCoord=0; iCoord<4; iCoord++) {
      traceData->startingCoord[iTrace][iCoord] = buffer[iTrace*4+iCoord];
    }
  }

  /* assert new element parameter values and recompute matrices */
  offset = traceData->traces*4;
  for (iUserParam=0; iUserParam<fitParam->parameters; iUserParam++) {
    *(fitParam->paramData[iUserParam]) = buffer[offset+iUserParam];
    if (fitParam->target[iUserParam]->matrix) {
      free_matrices(fitParam->target[iUserParam]->matrix);
      free(fitParam->target[iUserParam]->matrix);
    }
    compute_matrix(fitParam->target[iUserParam], run, NULL);
  }
}

double fit_trace_takeStep
  (
   MATRIX *D, 
   MATRIX *Dt, 
   MATRIX *DtD, 
   MATRIX *DtDInv, 
   MATRIX *DtDInvDt, 
   MATRIX *readbackVector, 
   MATRIX *paramVector,
   FIT_TRACE_DATA *traceData, 
   FIT_TRACE_PARAMETERS *fitParam, 
   LINE_LIST *beamline, 
   RUN *run,
   double convergenceFactor
   )
{
  double rmsError;
  long iTrace, iCoord, offset, iUserParam;

  /* find derivatives with respect to all parameters, both the user-declared ones
   * and the starting trajectory values 
   */
  fit_traces_findDerivatives(traceData, fitParam, D, beamline, run);
  
  /* solve for the parameter error values.  Parameters must have these
   * values subtracted off in order to remove the readback errors
   */
  if (!m_trans(Dt, D)) 
    m_error("Tranposing D");
  if (!m_mult(DtD, Dt, D))
    m_error("Multiplying Dt.D");
  if (!m_invert(DtDInv, DtD))
    m_error("Inverting Dt.D");
  if (!m_mult(DtDInvDt, DtDInv, Dt))
    m_error("Multiplying Inv(Dt.D) and Dt");
  if (!m_mult(paramVector, DtDInvDt, readbackVector))
    m_error("Multipling DtDInvDt and readbackVector");
  if (convergenceFactor>1 &&
      !m_scmul(paramVector, paramVector, convergenceFactor))
    m_error("Multipling paramVector by convergenceFactor");

  /* assert new trajectory starting values */
  for (iTrace=0; iTrace<traceData->traces; iTrace++) {
    for (iCoord=0; iCoord<4; iCoord++) {
      traceData->startingCoord[iTrace][iCoord] -= 
        paramVector->a[iTrace*4+iCoord][0];
    }
  }

  /* assert new element parameter values and recompute matrices */
  offset = traceData->traces*4;
  for (iUserParam=0; iUserParam<fitParam->parameters; iUserParam++) {
    *(fitParam->paramData[iUserParam]) -= 
      paramVector->a[offset+iUserParam][0];
    if (fitParam->target[iUserParam]->matrix) {
      free_matrices(fitParam->target[iUserParam]->matrix);
      free(fitParam->target[iUserParam]->matrix);
    }
    compute_matrix(fitParam->target[iUserParam], run, NULL);
  }

  /* compute RMS error and return it */
  rmsError = fit_trace_findReadbackErrors(readbackVector, traceData, beamline, run);
  return rmsError;
}


void fit_trace_setRowValues
  (
   FIT_OUTPUT_DATA *outputData, 
   long iteration, 
   double rmsError, 
   double lastRmsError, 
   double pass,
   double convergenceFactor, 
   FIT_TRACE_DATA *traceData, 
   FIT_TRACE_PARAMETERS *fitParam
   ) 
{
  long iTrace, iCoord, iUserParam;
  
  if (!SDDS_SetRowValues(&outputData->SDDStable, 
                         SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,
                         iteration, 
                         "Iteration", iteration, 
                         "RMSError", rmsError,
                         "Passes", pass,
                         "RMSErrorImprovement", (rmsError-lastRmsError)/lastRmsError, 
                         "ConvergenceFactor", convergenceFactor,
                         NULL)) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  for (iTrace=0; iTrace<traceData->traces; iTrace++) {
    for (iCoord=0; iCoord<4; iCoord++) {
      if (!SDDS_SetRowValues(&outputData->SDDStable,
                             SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE,
                             iteration, 
                             outputData->traceDataIndex[iTrace][iCoord],
                             traceData->startingCoord[iTrace][iCoord],
                             -1)) {
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      }
    }
  }
  for (iUserParam=0; iUserParam<fitParam->parameters; iUserParam++) {
    if (!SDDS_SetRowValues(&outputData->SDDStable,
                           SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE,
                           iteration,
                           outputData->paramDataIndex[iUserParam],
                           *(fitParam->paramData[iUserParam]), 
                           -1)) {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
  }
  if (!SDDS_UpdatePage(&outputData->SDDStable, 0)) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
}

