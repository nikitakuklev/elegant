/* Copyright 1994 by Michael Borland and Argonne National Laboratory,
 * all rights reserved.
 */
/* file    : extend_list.c
 * contents: extend_line_list(), extend_elem_list()
 *           Routines to extend linked lists for program mtt.c
 *
 * Michael Borland, 1987.
 */
#include "mdb.h"
#include "track.h"

void extend_line_list(LINE_LIST **lptr)
{
    log_entry("extend_line_list");
    (*lptr)->succ = tmalloc((unsigned)sizeof(**lptr));
    ((*lptr)->succ)->pred = *lptr;
    *lptr                 = (*lptr)->succ;
    (*lptr)->succ         = NULL;
    (*lptr)->links        = NULL;
    log_exit("extend_line_list");
    }

void extend_elem_list(ELEMENT_LIST **eptr)
{
    log_entry("extend_elem_list");
    (*eptr)->succ = tmalloc((unsigned)sizeof(**eptr));
    ((*eptr)->succ)->pred = *eptr;
    *eptr                 = (*eptr)->succ;
    (*eptr)->succ         = NULL;
    (*eptr)->matrix       = NULL;
    (*eptr)->type   = (*eptr)->flags = (*eptr)->end_pos = 0;
    (*eptr)->p_elem = (*eptr)->name  = NULL;
    log_exit("extend_elem_list");
    }

