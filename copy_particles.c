/* Copyright 1994 by Michael Borland and Argonne National Laboratory,
 * all rights reserved.
 */
/* routine: copy_particles()
 * purpose: copy particle data from one array to another
 *
 * Michael Borland, 1989
 */
#include "mdb.h"
#include "track.h"

void copy_particles(
    double **copy,
    double **original,
    long n_particles
    )
{
    register long ip, ic;
    double *cptr, *optr;

    log_entry("copy_particles");
    if (!copy)
        bomb("can't copy particles--target is NULL pointer (copy_particles)", NULL);
    if (!original)
        bomb("can't copy particles--source is NULL pointer (copy_particles)", NULL);

    for (ip=n_particles-1; ip>=0; ip--) {
        cptr = *copy++;
        optr = *original++;
        if (!cptr)
            bomb("element of target array is NULL pointer (copy_particles)", NULL);
        if (!optr)
            bomb("element of source array is NULL pointer (copy_particles)", NULL);
        for (ic=6; ic>=0; ic--)
            *cptr++ = *optr++;
        }
    log_exit("copy_particles");
    }
 
