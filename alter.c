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

static long printingEnabled = 1;

typedef struct {
  char *name, *item, *type, *exclude;
  double value;
  char *string_value;
  long differential, multiplicative, alter_at_each_step, alter_before_load_parameters, verbose, allow_missing_elements, allow_missing_parameters;
  long start_occurence, end_occurence, occurence_step;
  double s_start, s_end;
  char *after, *before;
} ALTER_SPEC;

static ALTER_SPEC *alterSpec = NULL;
static long alterSpecs = 0;

void reset_alter_specifications() {
  long i;
  for (i=0; i<alterSpecs; i++) {
    if (alterSpec[i].name) free(alterSpec[i].name);
    if (alterSpec[i].item) free(alterSpec[i].item);
    if (alterSpec[i].type) free(alterSpec[i].type);
    if (alterSpec[i].exclude) free(alterSpec[i].exclude);
    if (alterSpec[i].string_value) free(alterSpec[i].string_value);
    if (alterSpec[i].after) free(alterSpec[i].after);
    if (alterSpec[i].before) free(alterSpec[i].before);
  }
  free(alterSpec);
  alterSpec = NULL;
  alterSpecs = 0;
}

void setup_alter_element(NAMELIST_TEXT *nltext, RUN *run, LINE_LIST *beamline)
{
#include "alter.h"
    
#if !USE_MPI
    printingEnabled = 1;
#else
    printingEnabled = myid==1?1:0;
#endif

    /* process the namelist text */
    set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
    set_print_namelist_flags(0);
    if (processNamelist(&alter_elements, nltext)==NAMELIST_ERROR)
      bombElegant(NULL, NULL);
    if (echoNamelists) print_namelist(stdout, &alter_elements);

    if (!alter_at_each_step && alter_before_load_parameters) {
      printf("N.B.: alter_before_load_parameters has no effect unless alter_at_each_step=1\n");
      fflush(stdout);
    }
    
    if (!name || !strlen(name))
      bombElegant("no name given", NULL);
    if (has_wildcards(name) && strchr(name, '-'))
      name = expand_ranges(name);
    if (!item || !strlen(item))
      bombElegant("no item given", NULL);
    if (multiplicative) {
      if (!differential)
	/* convert to fractional change */
	value = value-1;
      differential = 0;
    }
    if ((((s_start>=0 && s_end>=0) ? 1 : 0) +
         ((start_occurence!=0 && end_occurence!=0) ? 1 : 0 ) +
         ((after || before) ? 1 : 0 ))>1)
      bombElegant("can't combine start_occurence/end_occurence, s_start/s_end, and after/before---use one method only", NULL);
    if (start_occurence>end_occurence) 
      bombElegant("start_occurence > end_occurence", NULL);
    if ((start_occurence!=0 && end_occurence!=0) && occurence_step<=0)
      bombElegant("occurence_step<=0", NULL);

    if (after || before) {
      ELEMENT_LIST *context;
      context = NULL;
      s_start = -DBL_MAX;
      s_end = DBL_MAX;
      if (after && strlen(after)) {
        if (!(context=find_element(after, &context, &(beamline->elem)))) {
	  if (printingEnabled)
	    fprintf(stdout, "Element %s not found in beamline.\n", after);
          exitElegant(1);
        }
        s_start = context->end_pos;
        if (find_element(after, &context, &(beamline->elem))) {
	  if (printingEnabled)
	    fprintf(stdout, "Element %s found in beamline more than once.\n", after);
          exitElegant(1);
        }
	if (printingEnabled) {
	  fprintf(stdout, "%s found at s = %le m\n", after, s_start);
	  fflush(stdout);
	}
      }
      context = NULL;
      if (before && strlen(before)) {
        if (!(context=find_element(before, &context, &(beamline->elem)))) {
	  if (printingEnabled)
	    fprintf(stdout, "Element %s not found in beamline.\n", before);
          exitElegant(1);
        }
        s_end = context->end_pos;
        if (find_element(before, &context, &(beamline->elem))) {
	  if (printingEnabled)
	    fprintf(stdout, "Element %s found in beamline more than once.\n", before);
          exitElegant(1);
        }
	if (printingEnabled) {
	  fprintf(stdout, "%s found at s = %le m\n", before, s_end);
	  fflush(stdout);
	}
      }
      if (s_start>s_end) 
        bombElegant("'after' element follows 'before' element!", NULL);
    }
    if (s_start>s_end)
      bombElegant("s_start > s_end", NULL);
    if (type) {
      long i;
      str_toupper(type);
      if (has_wildcards(type) && strchr(type, '-'))
        type = expand_ranges(type);
      for (i=0; i<N_TYPES; i++)
	if (wild_match(entity_name[i], type))
	  break;
      if (i==N_TYPES)
	bombElegant("type pattern does not match any known type", NULL);
    }    
    if (exclude && has_wildcards(exclude) && strchr(exclude, '-'))
      exclude = expand_ranges(exclude);

    alterSpec = SDDS_Realloc(alterSpec, sizeof(*alterSpec)*(alterSpecs+1));
    alterSpec[alterSpecs].name = 
      alterSpec[alterSpecs].item = 
        alterSpec[alterSpecs].type = 
          alterSpec[alterSpecs].exclude = 
            alterSpec[alterSpecs].string_value = 
              alterSpec[alterSpecs].after = 
                alterSpec[alterSpecs].before = NULL;
    cp_str(&alterSpec[alterSpecs].name, name);
    cp_str(&alterSpec[alterSpecs].item, item);
    if (type)
      cp_str(&alterSpec[alterSpecs].type, type);
    if (exclude)
      cp_str(&alterSpec[alterSpecs].exclude, exclude);
    if (string_value)
      cp_str(&alterSpec[alterSpecs].string_value, string_value);
    if (after)
      cp_str(&alterSpec[alterSpecs].after, after);
    if (before)
      cp_str(&alterSpec[alterSpecs].before, before);
    alterSpec[alterSpecs].value = value;
    alterSpec[alterSpecs].differential = differential;
    alterSpec[alterSpecs].multiplicative = multiplicative;
    alterSpec[alterSpecs].alter_at_each_step = alter_at_each_step;
    alterSpec[alterSpecs].alter_before_load_parameters = alter_before_load_parameters;
    alterSpec[alterSpecs].verbose = verbose;
    alterSpec[alterSpecs].allow_missing_elements = allow_missing_elements;
    alterSpec[alterSpecs].allow_missing_parameters = allow_missing_parameters;
    alterSpec[alterSpecs].start_occurence = start_occurence;
    alterSpec[alterSpecs].end_occurence = end_occurence;
    alterSpec[alterSpecs].occurence_step = occurence_step;
    alterSpec[alterSpecs].s_start = s_start;
    alterSpec[alterSpecs].s_end = s_end;
    alterSpecs++;
    
    do_alter_elements(run, beamline, 0, 0);
}

void do_alter_elements(RUN *run, LINE_LIST *beamline, short before_load_parameters, short per_step) 
{
    long thisType, lastType, iParam=0, nMatches;
    ELEMENT_LIST *context, *eptr;
    char *p_elem;
    char *p_elem0;
    char **changedDefinedParameter = NULL;
    long nChangedDefinedParameter = 0;
    long warningCountDown = 25;
    long i;
    
    for (i=0; i<alterSpecs; i++) {
      if (alterSpec[i].alter_at_each_step!=per_step)
        continue;
      if (alterSpec[i].alter_at_each_step && alterSpec[i].alter_before_load_parameters!=before_load_parameters)
        continue;
      context = NULL;
      lastType = -1;
      nMatches = 0;
      while ((eptr=wfind_element(alterSpec[i].name, &context, &(beamline->elem)))) {
        if (alterSpec[i].exclude && strlen(alterSpec[i].exclude) && wild_match(eptr->name, alterSpec[i].exclude))
          continue;
        if (alterSpec[i].start_occurence!=0 && alterSpec[i].end_occurence!=0) {
          if (eptr->occurence<alterSpec[i].start_occurence || eptr->occurence>alterSpec[i].end_occurence ||
              (eptr->occurence-alterSpec[i].start_occurence)%alterSpec[i].occurence_step!=0)
            continue;
        }
        if (alterSpec[i].s_start>=0 && alterSpec[i].s_end>=0 &&
            (eptr->end_pos<alterSpec[i].s_start || eptr->end_pos>alterSpec[i].s_end))
          continue;
        if (alterSpec[i].type && !wild_match(entity_name[context->type], alterSpec[i].type))
          continue;
        if ((thisType = eptr->type)!=lastType) {
          lastType = thisType;
          iParam = confirm_parameter(alterSpec[i].item, thisType);
        }
        if (iParam<0) {
          if (printingEnabled && warningCountDown>0) {
            fprintf(stderr, "%s: element %s does not have parameter %s\n", 
                    alterSpec[i].allow_missing_parameters?"Warning":"Error",
                    eptr->name, alterSpec[i].item);
            if (--warningCountDown==0 && alterSpec[i].allow_missing_parameters)
              fprintf(stderr, "*** Further messages suppressed!\n");
          }
          if (!alterSpec[i].allow_missing_parameters)
            exitElegant(1);
          continue;
        }
        nMatches++;
        p_elem = eptr->p_elem;
        p_elem0 = eptr->p_elem0;
        switch (entity_description[eptr->type].parameter[iParam].type) {
        case IS_DOUBLE:
          if (alterSpec[i].verbose && printingEnabled)
            fprintf(stdout, "Changing %s#%ld.%s from %21.15e to ",
                    eptr->name, eptr->occurence,
                    entity_description[eptr->type].parameter[iParam].name, 
                    *((double*)(p_elem+entity_description[eptr->type].parameter[iParam].offset)));
          /* this step could be very inefficient */
          if ((nMatches==1 || has_wildcards(alterSpec[i].name)) &&
              (!nChangedDefinedParameter ||
               match_string(eptr->name, changedDefinedParameter, 
                            nChangedDefinedParameter, EXACT_MATCH)==-1)) {
            change_defined_parameter(eptr->name, iParam, thisType, alterSpec[i].value, NULL, 
                                     alterSpec[i].differential?LOAD_FLAG_DIFFERENTIAL:
                                     (alterSpec[i].multiplicative?LOAD_FLAG_FRACTIONAL:LOAD_FLAG_ABSOLUTE));
            if (!(changedDefinedParameter=SDDS_Realloc(changedDefinedParameter,
                                                       sizeof(*changedDefinedParameter)*
                                                       (nChangedDefinedParameter+1)))) {
              bombElegant("memory allocation failure (alter_elements)", NULL);
            }
            changedDefinedParameter[nChangedDefinedParameter++] = eptr->name;
          }
          if (alterSpec[i].differential)
            *((double*)(p_elem+entity_description[eptr->type].parameter[iParam].offset)) += alterSpec[i].value;
          else if (alterSpec[i].multiplicative)
            *((double*)(p_elem+entity_description[eptr->type].parameter[iParam].offset)) *= 1+alterSpec[i].value;
          else
            *((double*)(p_elem+entity_description[eptr->type].parameter[iParam].offset)) = alterSpec[i].value;
          *((double*)(p_elem0+entity_description[eptr->type].parameter[iParam].offset)) = 
            *((double*)(p_elem+entity_description[eptr->type].parameter[iParam].offset)) ;
          if (alterSpec[i].verbose && printingEnabled) {
            fprintf(stdout, "%21.15e\n",
                    *((double*)(p_elem+entity_description[eptr->type].parameter[iParam].offset)));
            fflush(stdout);
          }
          break;
        case IS_LONG:
          if (alterSpec[i].verbose && printingEnabled)
            fprintf(stdout, "Changing %s#%ld.%s from %ld to ",
                    eptr->name, eptr->occurence,
                    entity_description[eptr->type].parameter[iParam].name, 
                    *((long*)(p_elem+entity_description[eptr->type].parameter[iParam].offset)));
          /* this step could be very inefficient */
          if ((nMatches==1 || has_wildcards(alterSpec[i].name)) &&
              (!nChangedDefinedParameter ||
               match_string(eptr->name, changedDefinedParameter, 
                            nChangedDefinedParameter, EXACT_MATCH)==-1)) {
            change_defined_parameter(eptr->name, iParam, thisType, alterSpec[i].value, NULL, 
                                     alterSpec[i].differential?LOAD_FLAG_DIFFERENTIAL:
                                     (alterSpec[i].multiplicative?LOAD_FLAG_FRACTIONAL:LOAD_FLAG_ABSOLUTE));
            if (!(changedDefinedParameter=SDDS_Realloc(changedDefinedParameter,
                                                       sizeof(*changedDefinedParameter)*
                                                       (nChangedDefinedParameter+1)))) {
              bombElegant("memory allocation failure (alter_elements)", NULL);
            }
            changedDefinedParameter[nChangedDefinedParameter++] = eptr->name;
          }
          if (alterSpec[i].differential)
            *((long*)(p_elem+entity_description[eptr->type].parameter[iParam].offset)) += 
              nearestInteger(alterSpec[i].value);
          else if (alterSpec[i].multiplicative)
            *((long*)(p_elem+entity_description[eptr->type].parameter[iParam].offset)) *= 1+alterSpec[i].value;
          else
            *((long*)(p_elem+entity_description[eptr->type].parameter[iParam].offset)) =
              nearestInteger(alterSpec[i].value);
          *((long*)(p_elem0+entity_description[eptr->type].parameter[iParam].offset)) = 
            *((long*)(p_elem+entity_description[eptr->type].parameter[iParam].offset)) ;
          if (alterSpec[i].verbose && printingEnabled) {
            fprintf(stdout, "%ld\n",
                    *((long*)(p_elem+entity_description[eptr->type].parameter[iParam].offset)));
            fflush(stdout);
          }
          break;
        case IS_STRING:
          if (alterSpec[i].string_value==NULL) {
            if (printingEnabled)
              fprintf(stderr, "Error: string_value is NULL for alter_elements, but parameter %s of %s is a string parameter\n",
                      entity_description[eptr->type].parameter[iParam].name, eptr->name);
            exitElegant(1);
          }
          /* unfortunately, can't free the existing pointer as I can't be sure that it isn't
           * pointing to static memory
           */
          if (alterSpec[i].verbose && printingEnabled)
            fprintf(stdout, "Changing %s#%ld.%s from %s to ",
                    eptr->name, eptr->occurence,
                    entity_description[eptr->type].parameter[iParam].name, 
                    *((char**)(p_elem+entity_description[eptr->type].parameter[iParam].offset)));
          cp_str((char**)(p_elem+entity_description[eptr->type].parameter[iParam].offset),
                 alterSpec[i].string_value);
          cp_str((char**)(p_elem0+entity_description[eptr->type].parameter[iParam].offset),
                 alterSpec[i].string_value);
          /* this step could be very inefficient */
          if ((nMatches==1 || has_wildcards(alterSpec[i].name)) &&
              (!nChangedDefinedParameter ||
               match_string(eptr->name, changedDefinedParameter, 
                            nChangedDefinedParameter, EXACT_MATCH)==-1)) {
            change_defined_parameter(eptr->name, iParam, thisType, alterSpec[i].value, alterSpec[i].string_value,
                                     alterSpec[i].differential?LOAD_FLAG_DIFFERENTIAL:
                                     (alterSpec[i].multiplicative?LOAD_FLAG_FRACTIONAL:LOAD_FLAG_ABSOLUTE));
            if (!(changedDefinedParameter=SDDS_Realloc(changedDefinedParameter,
                                                       sizeof(*changedDefinedParameter)*
                                                       (nChangedDefinedParameter+1)))) {
              bombElegant("memory allocation failure (alter_elements)", NULL);
            }
            changedDefinedParameter[nChangedDefinedParameter++] = eptr->name;
          }
          if (alterSpec[i].verbose && printingEnabled) {
            fprintf(stdout, "%s\n",
                    *((char**)(p_elem+entity_description[eptr->type].parameter[iParam].offset)));
            fflush(stdout);
          }
          break;
        }
        eptr->flags |= 
          PARAMETERS_ARE_PERTURBED |
            ((entity_description[eptr->type].parameter[iParam].flags&PARAM_CHANGES_MATRIX)?VMATRIX_IS_PERTURBED:0);
        if ((eptr->flags&PARAMETERS_ARE_PERTURBED) && (entity_description[eptr->type].flags&HAS_MATRIX) && eptr->matrix) {
          free_matrices(eptr->matrix);
          free(eptr->matrix);
          eptr->matrix = NULL;
        }
      }
      if (nMatches==0)  {
        if (alterSpec[i].allow_missing_elements) {
          if (printingEnabled)
            fprintf(stdout, "Warning: no matches for %s\n", alterSpec[i].name);
        } else {
          if (printingEnabled)
            fprintf(stdout, "Error: no matches for %s\n", alterSpec[i].name);
          exitElegant(1);
        }
      } else 
        compute_end_positions(beamline);
      if (changedDefinedParameter) {
        free(changedDefinedParameter);
        changedDefinedParameter = NULL;
        nChangedDefinedParameter = 0;
      }
    }
  }

