/* Copyright 1994 by Michael Borland and Argonne National Laboratory,
 * all rights reserved.
 */
/* routine: advance_values()
 * purpose: sequence an array of values in a systematic way such that
 *          a specified n-dimensional grid is covered.
 *
 * Michael Borland, 1988.
 */
#include "mdb.h"
#include "track.h"

long advance_values1(double *value, long n_values, long *value_index, double *initial, double *step, 
                            double **enumerated_value, long *counter, long *max_count, long *flags, long n_indices)
{
    long i, counter_changed;

    if ((counter_changed=advance_counter(counter, max_count, n_indices))<0)
        return(-1);

    for (i=0; i<n_values; i++)  {
        if (enumerated_value[i])
            value[i] = enumerated_value[i][counter[value_index[i]]];
        else {
            if (!(flags[i]&VARY_GEOMETRIC))
                value[i] = initial[i] + counter[value_index[i]]*step[i];
            else
                value[i] = initial[i]*ipow(step[i], counter[value_index[i]]);
            }
        }

    return(counter_changed);
    }



