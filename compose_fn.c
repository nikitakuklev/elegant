/* Copyright 1994 by Michael Borland and Argonne National Laboratory,
 * all rights reserved.
 */
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
