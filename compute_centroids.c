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

    for (i_coord=0; i_coord<6; i_coord++)
        sum[i_coord] = centroid[i_coord] = 0;

    for (i_part=0; i_part<n_part; i_part++) {
        part = coordinates[i_part];
        for (i_coord=0; i_coord<6; i_coord++) 
            sum[i_coord] += part[i_coord];
        }
    if (n_part)
        for (i_coord=0; i_coord<6; i_coord++)
            centroid[i_coord] = sum[i_coord]/n_part;

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

    }

void zero_beam_sums(
                    BEAM_SUMS *sums,
                    long n
                    )
{
  long i, j, k;
  for (i=0; i<n; i++) {
    for (j=0; j<4; j++)
      sums[i].maxabs[j] = 0;
    for (j=0; j<6; j++)
      sums[i].sum[j] = 0;
    for (j=0; j<6; j++)
      for (k=j; k<6; k++)
        sums[i].sum2[j][k] = 0;
#if DO_NORMEMIT_SUMS
    for (j=0; j<4; j++)
      sums[i].psum[j] = 0;
    for (j=0; j<4; j++)
      for (k=j; k<4; k++)
        sums[i].psum2[j][k] = 0;
#endif
    sums[i].n_part = sums[i].z = sums[i].p0 = 0;
  }
}

void accumulate_beam_sums(
                          BEAM_SUMS *sums,
                          double **coords,
                          long n_part,
                          double p_central
                          )
{
  double value;
  double *part, *part0;
  long i_part, i, j;
  double *sum, *maxabs;
  double beta0;
  
  sum  = sums->sum;
#if DO_NORMEMIT_SUMS
  psum = sums->psum;
#endif
  maxabs = sums->maxabs;
  if (!sums->n_part)
    sums->p0 = p_central;
  beta0 = sums->p0/sqrt(sqr(sums->p0)+1);
  sums->n_part += n_part;
  part0 = coords[0];

  for (i_part=n_part-1; i_part>=0; i_part--) {
    part = *coords++;
    for (i=0; i<4; i++)
      if ((value=FABS(part[i]))>maxabs[i])
        maxabs[i] = value;
    for (i=0; i<6; i++) {
      sum[i] += (value=part[i]);
      for (j=i; j<6; j++) {
        sums->sum2[i][j] += value*part[j];
      }
    }

#if DO_NORMEMIT_SUMS
    ds = part[5]-part0[5];
    pz = p_central*(1+part[5])/sqrt(1+sqr(part[1])+sqr(part[3]));
    for (i=0; i<4; i++)
      save[i] = part[i];
    ds = 0;
    part[0] -= ds*part[1];
    part[1] *= pz;
    part[2] -= ds*part[3];
    part[3] *= pz;
    for (i=0; i<4; i++) {
      psum[i] += (value=part[i]);
      for (j=i; j<4; j++)
        sums->psum2[i][j] += value*part[j];
    }
    for (i=0; i<4; i++)
      part[i] = save[i];
#endif
  }
}

void copy_beam_sums(
    BEAM_SUMS *target,
    BEAM_SUMS *source
    )
{
  long i, j;
  double *sum_s, *maxabs_s;
  double *sum_t, *maxabs_t;
#if DO_NORMEMIT_SUMS
  double *psum_s, *psum_t;
#endif
  sum_s    = source->sum;
  maxabs_s = source->maxabs;
  
  sum_t    = target->sum;
  maxabs_t = target->maxabs;

#if DO_NORMEMIT_SUMS
  psum_t   = target->psum;
  psum_s   = source->psum;
#endif

  for (i=0; i<4; i++) {
    maxabs_t[i] = maxabs_s[i];
  }
#if DO_NORMEMIT_SUMS
  for (i=0; i<2; i++)
    psum_t[i] = psum_s[i];
#endif
  for (i=0; i<6; i++) {
    sum_t[i] = sum_s[i];
    for (j=i; j<6; j++)
      target->sum2[i][j] = source->sum2[i][j];
#if DO_NORMEMIT_SUMS
    if (i<4) {
      psum_t[i] = psum_s[i];
      for (j=i; j<4; j++)
        target->psum2[i][j] = source->psum2[i][j];
    }
#endif
  }
  target->n_part = source->n_part;
  target->z      = source->z;
  target->p0     = source->p0;
}

