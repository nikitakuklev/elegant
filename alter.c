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
    long i, index, type, lastType, iParam, nMatches;
    ELEMENT_LIST *context, *eptr;
    char *p_elem;
    
    /* process the namelist text */
    set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
    set_print_namelist_flags(0);
    process_namelist(&alter_element, nltext);
    print_namelist(stdout, &alter_element);

    if (!strlen(name))
      bomb("no name given", NULL);
    if (!strlen(item))
      bomb("no item given", NULL);
    if (multiplicative)
      /* convert to fractional change */
      value = value-1;
    if (multiplicative && differential)
      bomb("can't combine multiplicative and differential modes", NULL);
    
    context = NULL;
    lastType = -1;
    nMatches = 0;
    while (eptr=wfind_element(name, &context, &(beamline->elem))) {
      if (exclude && strlen(exclude) && wild_match(eptr->name, exclude))
        continue;
      if ((type = eptr->type)!=lastType) {
        lastType = type;
        if ((iParam=confirm_parameter(item, type))<0) {
          fprintf(stderr, "%s: element %s does not have parameter %s\n", 
                  allow_missing_parameters?"Warning":"Error",
                  eptr->name, item);
          if (!allow_missing_parameters)
            exit(1);
        }
      }
      nMatches++;
      /* this step could be very inefficient */
      change_defined_parameter(eptr->name, iParam, type, value, NULL, 
                               differential?LOAD_FLAG_DIFFERENTIAL:(multiplicative?LOAD_FLAG_FRACTIONAL:LOAD_FLAG_ABSOLUTE));
      p_elem = eptr->p_elem;
      switch (entity_description[eptr->type].parameter[iParam].type) {
      case IS_DOUBLE:
        if (verbose)
          fprintf(stdout, "Changing %s.%s from %21.15e to ",
                  eptr->name, 
                  entity_description[eptr->type].parameter[iParam].name, 
                  *((double*)(p_elem+entity_description[eptr->type].parameter[iParam].offset)));
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
        if (differential)
          *((long*)(p_elem+entity_description[eptr->type].parameter[iParam].offset)) += value;
        else if (multiplicative)
          *((long*)(p_elem+entity_description[eptr->type].parameter[iParam].offset)) *= 1+value;
        else
          *((long*)(p_elem+entity_description[eptr->type].parameter[iParam].offset)) = value;
        if (verbose) {
          fprintf(stdout, "%ld\n",
                  *((long*)(p_elem+entity_description[eptr->type].parameter[iParam].offset)));
          fflush(stdout);
        }
        break;
      default:
        fprintf(stderr, "Can't alter non-numeric parameter %s for element %s\n", item, eptr->name);
        break;
      }
      eptr->flags |= 
        PARAMETERS_ARE_PERTURBED |
          (entity_description[eptr->type].parameter[iParam].changes_matrix?VMATRIX_IS_PERTURBED:0);
    }
    if (nMatches==0)
      fprintf(stdout, "Warning: no matches for %s\n", name);
}

