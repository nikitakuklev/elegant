/* Copyright 1994 by Michael Borland and Argonne National Laboratory,
 * all rights reserved.
 */
#include <stdio.h>

void link_date(void)
{
    fprintf(stdout, "Link date: %s %s\n", __DATE__, __TIME__);
    fflush(stdout);
    }

    
