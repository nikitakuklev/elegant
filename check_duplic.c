/* Copyright 1994 by Michael Borland and Argonne National Laboratory,
 * all rights reserved.
 */
/* file: check_duplic.c
 * contents: check_duplic_elem(), check_duplic_line()
 * purpose: check for duplicate definitions of elements or beamlines.
 *          Subroutines should be called after each element/line is added
 *          since they only check to see that the most recently added
 *          item is not a duplication of another.
 *
 * Michael Borland, 1988 
 */
#include "mdb.h"
#include "track.h"


void check_duplic_elem(
                       ELEMENT_LIST **elem0,          /* root of element list */
                       ELEMENT_LIST **new_elem,       /* newly added element  */
                       long n_elems
                       )
{
  char *new_name;
  long i, comparison, hi, lo, mid;
  ELEMENT_LIST *elem, *insertionPoint, *elast;
  ELEMENT_LIST **elemArray = NULL;

  if (!elem0 || !*elem0)
    bomb("root pointer of element list is null (check_duplic_elem)", NULL);
  if (!new_elem || !*new_elem)
    bomb("new element pointer is null (check_duplic_elem)", NULL);
  
#ifdef DEBUG    
  fprintf(stderr, "%ld elements in list\n",
         n_elems);
#endif
  if (n_elems>=1) {
    elemArray = trealloc(elemArray, sizeof(*elemArray)*(n_elems+10));
    elem = *elem0;
    for (i=0; i<n_elems; i++) {
      elemArray[i] = elem;
      elem = elem->succ;
    }
    
    new_name = (*new_elem)->name;
    elem = *elem0;
    insertionPoint = NULL;

    lo = 0; 
    hi = n_elems-1;
    if ((comparison=strcmp(new_name, elemArray[lo]->name))==0) {
      fprintf(stderr, "error: duplicate definitions of element %s\n",
              new_name);
      exit(1);
    }
    else if (comparison<0)
      insertionPoint = elemArray[lo];
    if (!insertionPoint) {
      if ((comparison=strcmp(new_name, elemArray[hi]->name))==0) {
        fprintf(stderr, "error: duplicate definitions of element %s\n",
                new_name);
        exit(1);
      }
      else if (comparison<0) {
        do {
          mid = (lo+hi)/2;
          if ((comparison=strcmp(new_name, elemArray[mid]->name))==0) {
            fprintf(stderr, "error: duplicate definitions of element %s\n",
                    new_name);
            exit(1);
          }
          else if (comparison<0)
            hi = mid;
          else
            lo = mid;
        } while ((hi-lo)>1);
        if (strcmp(new_name, elemArray[hi]->name)<0)
          insertionPoint = elemArray[hi];
      }
    }
    if (insertionPoint) {
#ifdef DEBUG
      fprintf(stderr, "inserting %s before %s\n", new_name, insertionPoint->name);
#endif
      elast = (*new_elem)->pred;
      if ((*new_elem)->pred)
        (*new_elem)->pred->succ = NULL;
      (*new_elem)->pred = insertionPoint->pred;
      if (insertionPoint->pred)
        insertionPoint->pred->succ = (*new_elem);
      (*new_elem)->succ = insertionPoint;
      insertionPoint->pred = (*new_elem);
      if (insertionPoint == *elem0)
        *elem0 = (*new_elem);
      *new_elem = elast;
    }
  }

#ifdef DEBUG
  elem  = *elem0;
  while (elem && elem->name) {
    fprintf(stderr, "  %s\n", elem->name);
    elem = elem->succ;
  }
  fflush(stderr);
#endif

}

void check_duplic_line(
                       LINE_LIST *line,          /* root of line list */
                       LINE_LIST *new_line,      /* newly added line */
                       long n_lines
                       )
{
  char *new_name;
  long i;

  new_name = new_line->name;
  
  n_lines--;
  for (i=0; i<n_lines; i++) {
    if (strcmp(new_name, line->name)==0) {
      *line->name = 0;
      fprintf(stderr, "warning: duplicate definitions of line %s\n", 
             new_name);
    }
    line = line->succ;
  }
}


