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
#include "insert_elements.h"

typedef struct {
  char *name, *type, *exclude, *elemDef;
  long nskip, add_end;
} ADD_ELEM;

static ADD_ELEM addElem;
static long add_elem_flag = 0;

long getAddElemFlag() 
{
  return (add_elem_flag);
}

char *getElemDefinition()
{
  return (addElem.elemDef);
} 

long getAddEndFlag()
{
  return (addElem.add_end);
} 

void do_insert_elements(NAMELIST_TEXT *nltext, RUN *run, LINE_LIST *beamline) 
{
  long i;

  /* process the namelist text */
  set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
  set_print_namelist_flags(0);
  process_namelist(&insert_elements, nltext);
  print_namelist(stdout, &insert_elements);

  if (disable)
    return;
  if (!skip && !add_at_end)
    bomb("skip and add_at_end can not be zero at the same time", NULL);

  add_elem_flag = 0;
  
  if ((!name || !strlen(name)) && !type)
    bomb("name or type needs to be given", NULL);
  if (!element_def || !strlen(element_def))
    bomb("element's definition is not given", NULL);
 
  if (name) {
    str_toupper(name);
    if (has_wildcards(name) && strchr(name, '-'))
      name = expand_ranges(name);
  }
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
  str_toupper(element_def);

  add_elem_flag = 1;
  addElem.nskip = skip;
  addElem.add_end =0;
  if (add_at_end)
    addElem.add_end = add_at_end;
  addElem.name = name;
  addElem.type = type;
  addElem.exclude = exclude;
  addElem.elemDef = element_def;

  delete_spaces(addElem.elemDef);
  beamline = get_beamline(NULL, beamline->name, run->p_central, 0);

  return;
}

long insertElem(char *name, long type, long *occurrence) 
{
  if (addElem.exclude && wild_match(name, addElem.exclude))
    return(0);
  if (addElem.name && !wild_match(name, addElem.name))
    return(0);
  if (addElem.type && !wild_match(entity_name[type], addElem.type))
    return(0);
  (*occurrence)++;

  if (*occurrence < addElem.nskip || addElem.nskip==0)
    return(0);

  *occurrence = 0;
  return(1);
}

