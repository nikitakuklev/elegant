/* Copyright 1994 by Michael Borland and Argonne National Laboratory,
 * all rights reserved.
 */
/* file: compute_centroids()
 * contents: compute_centroids(), compute_sigmas(), zero_beam_sums(),
 *           accumulate_beam_sums()
 *
 * Michael Borland, 1989
 */


/* routine: compute_centroids()
 * purpose: compute centroids for (x, x', y, y', s, deltap/p)
 *
 * Michael Borland, 1989
 */
#include "mdb.h"
#include "track.h"

void compute_centroids(
    double *centroid,
    double **coordinates,
    long n_part
    )
{
    long i_part, i_coord;
    double sum[6], *part;

    log_entry("compute_centroids");
    for (i_coord=0; i_coord<6; i_coord++)
        sum[i_coord] = 0;

    for (i_part=0; i_part<n_part; i_part++) {
        part = coordinates[i_part];
        for (i_coord=0; i_coord<6; i_coord++) 
            sum[i_coord] += part[i_coord];
        }
    if (n_part)
        for (i_coord=0; i_coord<6; i_coord++)
            centroid[i_coord] = sum[i_coord]/n_part;
    log_exit("compute_centroids");
    }

/* routine: compute_sigmas()
 * purpose: compute sigmas for (x, x', y, y', s, deltap/p)
 *
 * Michael Borland, 1989
 */

void compute_sigmas(
    double *sigma,
    double *centroid,
    double **coordinates,
    long n_part
    )
{
    long i_part, i_coord;
    double sum2[6], *part, value;

    log_entry("compute_sigmas");
    for (i_coord=0; i_coord<6; i_coord++)
        sum2[i_coord] = 0;

    for (i_part=0; i_part<n_part; i_part++) {
        part = coordinates[i_part];
        for (i_coord=0; i_coord<6; i_coord++) 
            sum2[i_coord] += sqr(part[i_coord]-centroid[i_coord]);
        }
    if (n_part)
        for (i_coord=0; i_coord<6; i_coord++)
            if ((value=sum2[i_coord])>0)
                sigma[i_coord] = sqrt(value/n_part);
            else
                sigma[i_coord] = 0;

    log_exit("compute_sigmas");
    }

void zero_beam_sums(
    BEAM_SUMS *sums,
    long n
    )
{
    long i, j;
    log_entry("zero_beam_sums");
    for (i=0; i<n; i++) {
        for (j=0; j<4; j++)
            sums[i].maxabs[j] = 0;
        for (j=0; j<6; j++)
            sums[i].sum[j] = 0;
        for (j=0; j<21; j++)
            sums[i].sum2[j] = 0;
        sums[i].n_part = sums[i].z = sums[i].p0 = 0;
        }
    log_exit("zero_beam_sums");
    }

void accumulate_beam_sums(
    BEAM_SUMS *sums,
    double **coords,
    long n_part,
    double p_central
    )
{
    double value;
    double *part;
    long i_part, i, j, index;
    double *sum, *sum2, *maxabs;

    log_entry("accumulate_beam_sums");
    sum  = sums->sum;
    sum2 = sums->sum2;
    maxabs = sums->maxabs;
    if (!sums->n_part)
        sums->p0 = p_central;
    sums->n_part += n_part;
    for (i_part=n_part-1; i_part>=0; i_part--) {
        part = *coords++;
        for (i=0; i<4; i++)
            if ((value=FABS(part[i]))>maxabs[i])
                maxabs[i] = value;
        for (i=index=0; i<6; i++) {
            sum[i] += (value=part[i]);
            for (j=i; j<6; j++)
                sum2[index++] += value*part[j];
            }
        }
    log_exit("accumulate_beam_sums");
    }

void copy_beam_sums(
    BEAM_SUMS *target,
    BEAM_SUMS *source
    )
{
    long i, j, index;
    double *sum_s, *sum2_s, *maxabs_s;
    double *sum_t, *sum2_t, *maxabs_t;

    log_entry("copy_beam_sums");
    sum_s    = source->sum;
    sum2_s   = source->sum2;
    maxabs_s = source->maxabs;
    sum_t    = target->sum;
    sum2_t   = target->sum2;
    maxabs_t = target->maxabs;

    for (i=0; i<4; i++)
        maxabs_t[i] = maxabs_s[i];
    for (i=index=0; i<6; i++) {
        sum_t[i] = sum_s[i];
        for (j=i; j<6; j++, index++)
            sum2_t[index] = sum2_s[index];
        }
    target->n_part = source->n_part;
    target->z      = source->z;
    target->p0     = source->p0;
    log_exit("copy_beam_sums");
    }

