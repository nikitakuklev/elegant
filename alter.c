/* Copyright 1999 by Michael Borland and Argonne National Laboratory,
 * all rights reserved.
 */
#include "mdb.h"
#include "track.h"
#include "alter.h"
#include "match_string.h"

#define DEBUG 0

long do_alter_element(NAMELIST_TEXT *nltext, RUN *run, LINE_LIST *beamline)
{
    long i, index, thisType, lastType, iParam, nMatches;
    ELEMENT_LIST *context, *eptr;
    char *p_elem, *sdata;
    
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
      str_toupper(type);
      if (has_wildcards(type) && strchr(type, '-'))
        type = expand_ranges(type);
    }
    if (exclude && has_wildcards(exclude) && strchr(exclude, '-'))
      exclude = expand_ranges(exclude);
      
    context = NULL;
    lastType = -1;
    nMatches = 0;
    while (eptr=wfind_element(name, &context, &(beamline->elem))) {
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
        change_defined_parameter(eptr->name, iParam, thisType, value, NULL, 
                                 differential?LOAD_FLAG_DIFFERENTIAL:(multiplicative?LOAD_FLAG_FRACTIONAL:LOAD_FLAG_ABSOLUTE));
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
        change_defined_parameter(eptr->name, iParam, thisType, value, NULL, 
                                 differential?LOAD_FLAG_DIFFERENTIAL:(multiplicative?LOAD_FLAG_FRACTIONAL:LOAD_FLAG_ABSOLUTE));
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
        change_defined_parameter(eptr->name, iParam, thisType, 0, string_value, 
                                 differential?LOAD_FLAG_DIFFERENTIAL:(multiplicative?LOAD_FLAG_FRACTIONAL:LOAD_FLAG_ABSOLUTE));
        if (verbose) {
          fprintf(stdout, "%s\n",
                  *((char**)(p_elem+entity_description[eptr->type].parameter[iParam].offset)));
          fflush(stdout);
        }
        break;
      }
      eptr->flags |= 
        PARAMETERS_ARE_PERTURBED |
          (entity_description[eptr->type].parameter[iParam].changes_matrix?VMATRIX_IS_PERTURBED:0);
    }
    if (nMatches==0)
      fprintf(stdout, "Warning: no matches for %s\n", name);
}

