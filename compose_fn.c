/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: compose_fn.c
 *
 * Michael Borland, 1994.
 */
#include "mdb.h"
#include "track.h"

char *compose_filename(char *template, char *root_name)
{
    char *ptr;
    
    if (str_in(template, "%s")) {
        ptr = tmalloc(sizeof(char)*(strlen(template)+strlen(root_name)+1));
        sprintf(ptr, template, root_name);
        return(ptr);
        }
    else
        return(template);
    }
