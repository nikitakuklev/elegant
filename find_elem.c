/* Copyright 1994 by Michael Borland and Argonne National Laboratory,
 * all rights reserved.
 */
/* file: find_elem.c
 *
 * Michael Borland, 1989-94
 */
#include "mdb.h"
#include "track.h"

ELEMENT_LIST *find_element(char *elem_name,  ELEMENT_LIST **context,  ELEMENT_LIST *elem)
{
    ELEMENT_LIST *eptr;

    log_entry("find_element");
    if (!elem_name)
        bomb("elem_name is NULL (find_element)", NULL);
    if (!elem)
        bomb("elem is NULL (find_element)", NULL);

    if (*context==NULL)
        eptr = elem;
    else
        eptr = *context = (*context)->succ;
    while (eptr && strcmp(eptr->name, elem_name)!=0)
        eptr = eptr->succ;
    log_exit("find_element");
    return(*context=eptr);
    }

ELEMENT_LIST *wfind_element(char *elem_name,  ELEMENT_LIST **context,  ELEMENT_LIST *elem)
{
    ELEMENT_LIST *eptr;

    log_entry("wfind_element");
    if (!elem_name)
        bomb("elem_name is NULL (wfind_element)", NULL);
    if (!elem)
        bomb("elem is NULL (wfind_element)", NULL);

    if (*context==NULL)
        eptr = elem;
    else
        eptr = *context = (*context)->succ;
    while (eptr && !wild_match(eptr->name, elem_name))
        eptr = eptr->succ;

    log_exit("wfind_element");
    return(*context=eptr);
    }

long confirm_parameter(char *item_name, long type)
{
    long i;
    PARAMETER *param;

    log_entry("confirm_parameter");

    if (!item_name)
        bomb("item_name is NULL (confirm_parameter)", NULL);
    
    param = entity_description[type].parameter;
    for (i=0; i<entity_description[type].n_params; i++) 
        if (strcmp(item_name, param[i].name)==0) {
            log_exit("confirm_parameter");
            return(i);
            }
    log_exit("confirm_parameter");
    return(-1);
    }

