/*************************************************************************\
* Copyright (c) 2010 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2010 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

#include "mdb.h"
#include "track.h"
#include "modulate.h"
#ifdef HAVE_GPU
#include <gpu_base.h>
#include <gpu_simple_rfca.h>
#endif

long loadModulationTable(double **t, double **value, char *file, char *timeColumn, char *amplitudeColumn);

void addModulationElements(MODULATION_DATA *modData, NAMELIST_TEXT *nltext, LINE_LIST *beamline, RUN *run)
{
  long n_items, n_added, firstIndexInGroup;
  ELEMENT_LIST *context;
  double sMin = -DBL_MAX, sMax = DBL_MAX;
  double *tData, *AData;
  long nData = 0;
  
  /* process namelist text */
  set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
  set_print_namelist_flags(0);
  if (processNamelist(&modulate_elements, nltext)==NAMELIST_ERROR)
    bombElegant(NULL, NULL);
  if (name==NULL) {
    if (!type)
      bombElegant("element name missing in modulate_elements namelist", NULL);
    SDDS_CopyString(&name, "*");
  }
  if (!expression && !(filename && time_column && amplitude_column))
    bombElegant("either expression or filename, time_column, and amplitude_column must all be given", NULL);
  if (expression && filename)
    bombElegant("only one of expression and filename may be given", NULL);
  if (item==NULL)
    bombElegant("item name missing in modulate_elements namelist", NULL);
  if (echoNamelists) print_namelist(stdout, &modulate_elements);

  if (filename) {
    /* Read data file */
    if ((nData = loadModulationTable(&tData, &AData, filename, time_column, amplitude_column))<=2)
      bombElegant("too few items in modulation table", NULL);
  } else
    nData = 0;
  
  n_added = 0;
  n_items = modData->nItems;
  context = NULL;
  if (after && strlen(after)) {
    if (!(context=find_element(after, &context, &(beamline->elem)))) {
      fprintf(stdout, "Element %s not found in beamline.\n", after);
      exitElegant(1);
    }
    sMin = context->end_pos;
    if (find_element(after, &context, &(beamline->elem))) {
      fprintf(stdout, "Element %s found in beamline more than once.\n", after);
      exitElegant(1);
    }
    fprintf(stdout, "%s found at s = %le m\n", after, sMin);
    fflush(stdout);
  }
  context = NULL;
  if (before && strlen(before)) {
    if (!(context=find_element(before, &context, &(beamline->elem)))) {
      fprintf(stdout, "Element %s not found in beamline.\n", before);
      exitElegant(1);
    }
    sMax = context->end_pos;
    if (find_element(before, &context, &(beamline->elem))) {
      fprintf(stdout, "Element %s found in beamline more than once.\n", after);
      exitElegant(1);
    }
    fprintf(stdout, "%s found at s = %le m\n", before, sMax);
    fflush(stdout);
  }
  if (after && before && sMin>sMax) {
    fprintf(stdout, "Element %s is not upstream of %s!\n",
            before, after);
    exitElegant(1);
  }
  if (type && has_wildcards(type) && strchr(type, '-'))
    type = expand_ranges(type);
  if (has_wildcards(name)) {
    if (strchr(name, '-'))
      name = expand_ranges(name);
    str_toupper(name);
  }
  
  firstIndexInGroup = -1;
  while ((context=wfind_element(name, &context, &(beamline->elem)))) {
    if (type && !wild_match(entity_name[context->type], type))
      continue;
    if ((sMin>=0 && context->end_pos<sMin) ||
        (sMax>=0 && context->end_pos>sMax) ||
        (s_start>=0 && context->end_pos<s_start) ||
        (s_end>=0 && context->end_pos>s_end) ||
        (start_occurence && context->occurence<start_occurence) ||
        (end_occurence && context->occurence>end_occurence) )
      continue;
    fprintf(stdout, "Adding modulation for %s#%ld at s=%le\n", context->name, context->occurence, context->end_pos);
    
    modData->element          = SDDS_Realloc(modData->element, sizeof(*modData->element)*(n_items+1));
    modData->expression       = SDDS_Realloc(modData->expression, sizeof(*modData->expression)*(n_items+1));
    modData->parameterNumber  = SDDS_Realloc(modData->parameterNumber, sizeof(*modData->parameterNumber)*(n_items+1));
    modData->flags            = SDDS_Realloc(modData->flags, sizeof(*modData->flags)*(n_items+1));
    modData->verboseThreshold = SDDS_Realloc(modData->verboseThreshold, sizeof(*modData->verboseThreshold)*(n_items+1));
    modData->lastVerboseValue = SDDS_Realloc(modData->lastVerboseValue, sizeof(*modData->lastVerboseValue)*(n_items+1));
    modData->unperturbedValue = SDDS_Realloc(modData->unperturbedValue, sizeof(*modData->unperturbedValue)*(n_items+1));
    modData->nData            = SDDS_Realloc(modData->nData, sizeof(*modData->nData)*(n_items+1));
    modData->dataIndex        = SDDS_Realloc(modData->dataIndex, sizeof(*modData->dataIndex)*(n_items+1));
    modData->timeData         = SDDS_Realloc(modData->timeData, sizeof(*modData->timeData)*(n_items+1));
    modData->modulationData   = SDDS_Realloc(modData->modulationData, sizeof(*modData->modulationData)*(n_items+1));
    modData->record           = SDDS_Realloc(modData->record, sizeof(*modData->record)*(n_items+1));
    modData->flushRecord      = SDDS_Realloc(modData->flushRecord, sizeof(*modData->flushRecord)*(n_items+1));
    modData->fpRecord         = SDDS_Realloc(modData->fpRecord, sizeof(*modData->fpRecord)*(n_items+1));

    modData->element[n_items] = context;
    modData->flags[n_items] = (multiplicative?MULTIPLICATIVE_MOD:0) + (differential?DIFFERENTIAL_MOD:0) 
      + (verbose?VERBOSE_MOD:0) + (refresh_matrix?REFRESH_MATRIX_MOD:0);
    modData->verboseThreshold[n_items] = verbose_threshold;
    modData->timeData[n_items] = modData->modulationData[n_items] = NULL;
    modData->expression[n_items] = NULL;
    modData->fpRecord[n_items] = NULL;
    modData->nData[n_items] = 0;
    modData->flushRecord[n_items] = flush_record;

    if (filename) {
      if ((modData->dataIndex[n_items] = firstIndexInGroup)==-1) {
        modData->timeData[n_items]      = tData;
        modData->modulationData[n_items] = AData;
        modData->nData[n_items] = nData;
      }
    } else
      cp_str(&modData->expression[n_items], expression);

    if ((modData->parameterNumber[n_items] = confirm_parameter(item, context->type))<0) {
      fprintf(stdout, "error: cannot modulate %s---no such parameter for %s (wildcard name: %s)\n",item, context->name, name);
      fflush(stdout);
      exitElegant(1);
    }

    modData->lastVerboseValue[n_items] = modData->unperturbedValue[n_items] 
      = parameter_value(context->name, context->type, modData->parameterNumber[n_items], beamline);

    if (modData->unperturbedValue[n_items]==0 && modData->flags[n_items]&MULTIPLICATIVE_MOD) {
      fprintf(stdout, "***\7\7\7 warning: you've specified multiplicative modulation for %s.%s, but the unperturbed value is zero.\nThis may be an error.\n", 
              context->name, item);
      fflush(stdout);
    }

    if (record 
#if USE_MPI
        && myid==0
#endif
        ) {
      modData->record[n_items] = compose_filename(record, run->rootname);
      record = NULL;
      if (!(modData->fpRecord[n_items] = fopen(modData->record[n_items], "w")))
        SDDS_Bomb("problem setting up  modulation record file");
      fprintf(modData->fpRecord[n_items], "SDDS1\n&column name=t, units=s, type=double &end\n");
      fprintf(modData->fpRecord[n_items], "&column name=Pass, type=long &end\n");
      fprintf(modData->fpRecord[n_items], "&column name=Amplitude, type=double &end\n");
      fprintf(modData->fpRecord[n_items], "&column name=OriginalValue, type=double &end\n");
      fprintf(modData->fpRecord[n_items], "&column name=NewValue, type=double &end\n");
      fprintf(modData->fpRecord[n_items], "&data mode=ascii, no_row_counts=1 &end\n");
    }
    
    modData->nItems = ++n_items;
    n_added++;
    if (firstIndexInGroup==-1)
      firstIndexInGroup = n_items-1;
  }
    
  if (!n_added) {
    fprintf(stdout, "error: no match given modulation\n");
    fflush(stdout);
    exitElegant(1);
  }
}

long loadModulationTable(double **t, double **value, char *file, char *timeColumn, char *amplitudeColumn)
{
  SDDS_TABLE SDDS_table;
  long i, count=0;

  if (!t || !value || !file || !timeColumn || !amplitudeColumn)
    bombElegant("NULL value pointer passed (loadModulationTable)", NULL);
  
  if (!SDDS_InitializeInputFromSearchPath(&SDDS_table, file) || 
      SDDS_ReadTable(&SDDS_table)!=1 ||
      !(count=SDDS_CountRowsOfInterest(&SDDS_table)) ||
      !(*value = SDDS_GetColumnInDoubles(&SDDS_table, amplitudeColumn)) ||
      !(*t = SDDS_GetColumnInDoubles(&SDDS_table, timeColumn)) ||
      !SDDS_Terminate(&SDDS_table)) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    bombElegant("Unable to read required data from modulation file.", NULL);
  }
  for (i=1; i<count; i++)
    if ((*t)[i-1]>=(*t)[i])
      bombElegant("time values not monotonically increasing for modulation data", NULL);
  return(count);
}

long applyElementModulations(MODULATION_DATA *modData, double pCentral, double **coord, long np, RUN *run, long iPass)
{
  long iMod, code, matricesUpdated, jMod;
  short modulationValid = 0;
  double modulation, value, t, lastValue;
  long type, param;
  char *p_elem;
 
  if (modData->nItems<=0)
    return 0;

#ifdef HAVE_GPU
  if (getGpuBase()->elementOnGpu)
    t = gpu_findFiducialTime(np,0, 0, pCentral, FID_MODE_TMEAN);
  else 
#endif
    t = findFiducialTime(coord, np, 0, 0, pCentral, FID_MODE_TMEAN);
  matricesUpdated = 0;


  for (iMod=0; iMod<modData->nItems; iMod++) {
    type = modData->element[iMod]->type;
    param = modData->parameterNumber[iMod];
    p_elem = (char*)(modData->element[iMod]->p_elem);
    value = modData->unperturbedValue[iMod];
    modulation = 0;
    
    if (!modData->expression[iMod]) {
      jMod = iMod;
      if (modData->dataIndex[iMod]!=-1)
        jMod = modData->dataIndex[iMod];
      code = 1;
      if (t<=modData->timeData[jMod][0]) {
        modulation = modData->modulationData[jMod][0];
        fprintf(stderr, "Warning: interpolation at t=%21.15le below modulation table range [%21.15le, %21.15le] for element %s, parameter %s\n",
                t, modData->timeData[jMod][0], modData->timeData[jMod][modData->nData[jMod]-1], modData->element[jMod]->name,
                entity_description[type].parameter[param].name);
      } else if (t>=modData->timeData[jMod][modData->nData[jMod]-1]) {
        modulation = modData->modulationData[jMod][modData->nData[jMod]-1];
        fprintf(stderr, "Warning: interpolation at t=%21.15le above modulation table range [%21.15le, %21.15le] for element %s, parameter %s\n",
                t, modData->timeData[jMod][0], modData->timeData[jMod][modData->nData[jMod]-1], modData->element[jMod]->name,
                entity_description[type].parameter[param].name);
      } else
        modulation = interp(modData->modulationData[jMod], modData->timeData[jMod], modData->nData[jMod], t, 0, 1, &code);
      if (code==0) {
        fprintf(stderr, "Error: interpolation failed for t=%21.15le for element %s, parameter %s\n",
                t, modData->element[jMod]->name, entity_description[type].parameter[param].name);
        exitElegant(1);
      }
      modulationValid = 1;
    } else {
      push_num(t);
      modulation = rpn(modData->expression[iMod]);
      rpn_clear();
    }

    if (modData->flags[iMod]&DIFFERENTIAL_MOD) {
      if (modData->flags[iMod]&MULTIPLICATIVE_MOD)
        value = (1+modulation)*value;
      else
        value = value + modulation;
    } else {
      if (modData->flags[iMod]&MULTIPLICATIVE_MOD)
        value = value*modulation;
      else
        value = modulation;
    }

    switch (entity_description[type].parameter[param].type)  {
    case IS_DOUBLE:
      lastValue = modData->lastVerboseValue[iMod];
      *((double*)(p_elem+entity_description[type].parameter[param].offset)) = value;
      if (modData->flags[iMod]&VERBOSE_MOD && fabs(value-lastValue)>modData->verboseThreshold[iMod]*(fabs(value)+fabs(lastValue))/2) {
        fprintf(stdout, "Modulation value for element %s#%ld, parameter %s changed to %21.15le at t = %21.15le (originally %21.15le)\n",
                modData->element[iMod]->name, modData->element[iMod]->occurence,
                entity_description[type].parameter[param].name, value, t, modData->unperturbedValue[iMod]);
	modData->lastVerboseValue[iMod] = value;
      }
      break;
    case IS_LONG:
      lastValue = modData->lastVerboseValue[iMod];
      value = (long)(value+0.5);
      *((long*)(p_elem+entity_description[type].parameter[param].offset)) = value;
      if (modData->flags[iMod]&VERBOSE_MOD && value!=lastValue) {
        fprintf(stdout, "Modulation value for element %s#%ld, parameter %s changed to %ld at t = %21.15le (originally %ld)\n",
                modData->element[iMod]->name, modData->element[iMod]->occurence,
                entity_description[type].parameter[param].name, (long)(value+0.5), t, (long)(modData->unperturbedValue[iMod]));
	modData->lastVerboseValue[iMod] = value;
      }
      break;
    default:
      break;
    }

    if (modData->fpRecord[iMod] 
#if USE_MPI
        && myid==0
#endif
        ) {
      fprintf(modData->fpRecord[iMod], "%21.15le %ld %21.15le %21.15le %21.15le\n",
              t, iPass, modulation, modData->unperturbedValue[iMod], value);
      if (modData->flushRecord[iMod]>0 && (iPass%modData->flushRecord[iMod])==0)
	fflush(modData->fpRecord[iMod]);
    }
    
    if (entity_description[type].flags&HAS_MATRIX && 
        entity_description[type].parameter[param].flags&PARAM_CHANGES_MATRIX &&
        ((modData->flags[iMod]&REFRESH_MATRIX_MOD) || (entity_description[type].flags&MATRIX_TRACKING))) {
      /* update the matrix */
      if (modData->element[iMod]->matrix) {
        free_matrices(modData->element[iMod]->matrix);
        tfree(modData->element[iMod]->matrix);
        modData->element[iMod]->matrix = NULL;
      }
      compute_matrix(modData->element[iMod], run, NULL);
      matricesUpdated ++;
    }
  }

  return matricesUpdated;
}



      
