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
#include "alter.h"
#include "match_string.h"

#define DEBUG 0

void do_alter_element(NAMELIST_TEXT *nltext, RUN *run, LINE_LIST *beamline)
{
    long thisType, lastType, iParam=0, nMatches;
    ELEMENT_LIST *context, *eptr;
    char *p_elem;
    char **changedDefinedParameter = NULL;
    long nChangedDefinedParameter = 0;

    /* process the namelist text */
    set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
    set_print_namelist_flags(0);
    process_namelist(&alter_elements, nltext);
    print_namelist(stdout, &alter_elements);

    if (!name || !strlen(name))
      bomb("no name given", NULL);
    if (has_wildcards(name) && strchr(name, '-'))
      name = expand_ranges(name);
    if (!item || !strlen(item))
      bomb("no item given", NULL);
    if (multiplicative)
      /* convert to fractional change */
      value = value-1;
    if (multiplicative && differential)
      bomb("can't combine multiplicative and differential modes", NULL);
    if (type) {
      long i;
      str_toupper(type);
      if (has_wildcards(type) && strchr(type, '-'))
        type = expand_ranges(type);
      for (i=0; i<N_TYPES; i++)
	if (wild_match(entity_name[i], type))
	  break;
      if (i==N_TYPES)
	bomb("type pattern does not match any known type", NULL);
    }    
    if (exclude && has_wildcards(exclude) && strchr(exclude, '-'))
      exclude = expand_ranges(exclude);
      
    context = NULL;
    lastType = -1;
    nMatches = 0;
    while ((eptr=wfind_element(name, &context, &(beamline->elem)))) {
      if (exclude && strlen(exclude) && wild_match(eptr->name, exclude))
        continue;
      if (type && !wild_match(entity_name[context->type], type))
        continue;
      if ((thisType = eptr->type)!=lastType) {
        lastType = thisType;
        if ((iParam=confirm_parameter(item, thisType))<0) {
          fprintf(stderr, "%s: element %s does not have parameter %s\n", 
                  allow_missing_parameters?"Warning":"Error",
                  eptr->name, item);
          if (!allow_missing_parameters)
            exit(1);
        }
      }
      nMatches++;
      p_elem = eptr->p_elem;
      switch (entity_description[eptr->type].parameter[iParam].type) {
      case IS_DOUBLE:
        if (verbose)
          fprintf(stdout, "Changing %s.%s from %21.15e to ",
                  eptr->name, 
                  entity_description[eptr->type].parameter[iParam].name, 
                  *((double*)(p_elem+entity_description[eptr->type].parameter[iParam].offset)));
        /* this step could be very inefficient */
        if ((nMatches==1 || has_wildcards(name)) &&
	    (!nChangedDefinedParameter ||
	    match_string(eptr->name, changedDefinedParameter, 
			 nChangedDefinedParameter, EXACT_MATCH)==-1)) {
          change_defined_parameter(eptr->name, iParam, thisType, value, NULL, 
                                   differential?LOAD_FLAG_DIFFERENTIAL:
                                   (multiplicative?LOAD_FLAG_FRACTIONAL:LOAD_FLAG_ABSOLUTE));
	  if (!(changedDefinedParameter=SDDS_Realloc(changedDefinedParameter,
						    sizeof(*changedDefinedParameter)*
						     (nChangedDefinedParameter+1)))) {
	    bomb("memory allocation failure (alter_elements)", NULL);
	  }
	  changedDefinedParameter[nChangedDefinedParameter++] = eptr->name;
	}
        if (differential)
          *((double*)(p_elem+entity_description[eptr->type].parameter[iParam].offset)) += value;
        else if (multiplicative)
          *((double*)(p_elem+entity_description[eptr->type].parameter[iParam].offset)) *= 1+value;
        else
          *((double*)(p_elem+entity_description[eptr->type].parameter[iParam].offset)) = value;
        if (verbose) {
          fprintf(stdout, "%21.15e\n",
                  *((double*)(p_elem+entity_description[eptr->type].parameter[iParam].offset)));
            fflush(stdout);
        }
        break;
      case IS_LONG:
        if (verbose)
          fprintf(stdout, "Changing %s.%s from %ld to ",
                  eptr->name, 
                  entity_description[eptr->type].parameter[iParam].name, 
                  *((long*)(p_elem+entity_description[eptr->type].parameter[iParam].offset)));
        /* this step could be very inefficient */
        if ((nMatches==1 || has_wildcards(name)) &&
	    (!nChangedDefinedParameter ||
	    match_string(eptr->name, changedDefinedParameter, 
			 nChangedDefinedParameter, EXACT_MATCH)==-1)) {
          change_defined_parameter(eptr->name, iParam, thisType, value, NULL, 
                                   differential?LOAD_FLAG_DIFFERENTIAL:
                                   (multiplicative?LOAD_FLAG_FRACTIONAL:LOAD_FLAG_ABSOLUTE));
	  if (!(changedDefinedParameter=SDDS_Realloc(changedDefinedParameter,
						    sizeof(*changedDefinedParameter)*
						     (nChangedDefinedParameter+1)))) {
	    bomb("memory allocation failure (alter_elements)", NULL);
	  }
	  changedDefinedParameter[nChangedDefinedParameter++] = eptr->name;
	}
        if (differential)
          *((long*)(p_elem+entity_description[eptr->type].parameter[iParam].offset)) += 
            nearestInteger(value);
        else if (multiplicative)
          *((long*)(p_elem+entity_description[eptr->type].parameter[iParam].offset)) *= 1+value;
        else
          *((long*)(p_elem+entity_description[eptr->type].parameter[iParam].offset)) =
            nearestInteger(value);
        if (verbose) {
          fprintf(stdout, "%ld\n",
                  *((long*)(p_elem+entity_description[eptr->type].parameter[iParam].offset)));
          fflush(stdout);
        }
        break;
      case IS_STRING:
        if (string_value==NULL) {
          fprintf(stderr, "Error: string_value is NULL for alter_elements, but parameter %s of %s is a string parameter\n",
                  entity_description[eptr->type].parameter[iParam].name, eptr->name);
          exit(1);
        }
        /* unfortunately, can't free the existing pointer as I can't be sure that it isn't
         * pointing to static memory
         */
        if (verbose)
          fprintf(stdout, "Changing %s.%s from %s to ",
                  eptr->name, 
                  entity_description[eptr->type].parameter[iParam].name, 
                  *((char**)(p_elem+entity_description[eptr->type].parameter[iParam].offset)));
        cp_str((char**)(p_elem+entity_description[eptr->type].parameter[iParam].offset),
               string_value);
        /* this step could be very inefficient */
        if ((nMatches==1 || has_wildcards(name)) &&
	    (!nChangedDefinedParameter ||
	    match_string(eptr->name, changedDefinedParameter, 
			 nChangedDefinedParameter, EXACT_MATCH)==-1)) {
          change_defined_parameter(eptr->name, iParam, thisType, value, NULL, 
                                   differential?LOAD_FLAG_DIFFERENTIAL:
                                   (multiplicative?LOAD_FLAG_FRACTIONAL:LOAD_FLAG_ABSOLUTE));
	  if (!(changedDefinedParameter=SDDS_Realloc(changedDefinedParameter,
						    sizeof(*changedDefinedParameter)*
						     (nChangedDefinedParameter+1)))) {
	    bomb("memory allocation failure (alter_elements)", NULL);
	  }
	  changedDefinedParameter[nChangedDefinedParameter++] = eptr->name;
	}
        if (verbose) {
          fprintf(stdout, "%s\n",
                  *((char**)(p_elem+entity_description[eptr->type].parameter[iParam].offset)));
          fflush(stdout);
        }
        break;
      }
      eptr->flags |= 
        PARAMETERS_ARE_PERTURBED |
          ((entity_description[eptr->type].parameter[iParam].flags&PARAM_CHANGES_MATRIX)?VMATRIX_IS_PERTURBED:0);
    }
    if (nMatches==0)  {
      if (allow_missing_elements)
        fprintf(stdout, "Warning: no matches for %s\n", name);
      else {
        fprintf(stdout, "Error: no matches for %s\n", name);
        exit(1);
      }
    }
    if (changedDefinedParameter)
      free(changedDefinedParameter);
  }


