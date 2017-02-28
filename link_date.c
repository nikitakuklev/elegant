/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

#include <stdio.h>

void link_date(void)
{
#ifdef SVN_VERSION
    printf("Link date: %s %s, SVN revision: %s\n", __DATE__, __TIME__, SVN_VERSION);
#else
    printf("Link date: %s %s\n", __DATE__, __TIME__);
#endif
    fflush(stdout);
    }

    
