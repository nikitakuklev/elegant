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
    double *emit,
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
    if (emit)
      for (i_coord=0; i_coord<3;  i_coord++)
        emit[i_coord] = 0;
    
    for (i_part=0; i_part<n_part; i_part++) {
        part = coordinates[i_part];
        for (i_coord=0; i_coord<6; i_coord++) 
            sum2[i_coord] += sqr(part[i_coord]-centroid[i_coord]);
        }
    if (n_part) {
        for (i_coord=0; i_coord<6; i_coord++)
            if ((value=sum2[i_coord])>0)
                sigma[i_coord] = sqrt(value/n_part);
            else
                sigma[i_coord] = 0;
        if (emit) {
          for (i_coord=0; i_coord<6; i_coord+=2) {
            double sum12;
            sum12 = 0;
            for (i_part=0; i_part<n_part; i_part++) {
              part = coordinates[i_part];
              sum12 += (part[i_coord]-centroid[i_coord])*(part[i_coord+1]-centroid[i_coord+1]);
            }
            if ((emit[i_coord/2] = sqr(sigma[i_coord]*sigma[i_coord+1]) - sqr(sum12/n_part))>0)
              emit[i_coord/2] = sqrt(emit[i_coord/2]);
            else
              emit[i_coord/2] = 0;
          }
        }
      }
    }

void zero_beam_sums(
                    BEAM_SUMS *sums,
                    long n
                    )
{
  long i, j, k;
  for (i=0; i<n; i++) {
    for (j=0; j<6; j++)
      sums[i].maxabs[j] = 0;
    for (j=0; j<6; j++)
      sums[i].centroid[j] = 0;
    for (j=0; j<6; j++)
      for (k=j; k<6; k++)
        sums[i].sigma[j][k] = 0;
    sums[i].n_part = sums[i].z = sums[i].p0 = 0;
  }
}

void accumulate_beam_sums(
                          BEAM_SUMS *sums,
                          double **coord,
                          long n_part,
                          double p_central
                          )
{
  long i_part, i, j;
  double centroid[6];
  double Sij, value;
  
  if (!sums->n_part)
    sums->p0 = p_central;

  if (n_part) {
    /* maximum amplitudes */
    for (i=0; i<6; i++) {
      if (i==4)
        continue;  /* done below */
      for (i_part=0; i_part<n_part; i_part++) {
        if ((value=fabs(coord[i_part][i]))>sums->maxabs[i])
          sums->maxabs[i] = value;
      }
    }
    /* compute centroids for present beam and add in to existing centroid data */
    for (i=0; i<6; i++) {
      for (centroid[i]=i_part=0; i_part<n_part; i_part++)
        centroid[i] += coord[i_part][i];
      sums->centroid[i] = (sums->centroid[i]*sums->n_part+centroid[i])/(sums->n_part+n_part);
      centroid[i] /= n_part;
    }
    i = 4;
    for (i_part=0; i_part<n_part; i_part++)
      if ((value=fabs(coord[i_part][i]-centroid[i]))>sums->maxabs[i])
        sums->maxabs[i] = value;
    /* compute Sigma[i][j] for present beam and add to existing data */
    for (i=0; i<6; i++)
      for (j=i; j<6; j++) {
        for (Sij=i_part=0; i_part<n_part; i_part++)
          Sij += (coord[i_part][i]-centroid[i])*(coord[i_part][j]-centroid[j]);
        sums->sigma[i][j] = (sums->sigma[i][j]*sums->n_part + Sij)/(sums->n_part+n_part);
      }
  }
  
  sums->n_part += n_part;
}

void copy_beam_sums(
    BEAM_SUMS *target,
    BEAM_SUMS *source
    )
{
  long i, j;

  for (i=0; i<6; i++) {
    target->maxabs[i] = source->maxabs[i];
  }
  for (i=0; i<6; i++) {
    target->centroid[i] = source->centroid[i];
    for (j=i; j<6; j++)
      target->sigma[i][j] = source->sigma[i][j];
  }
  target->n_part = source->n_part;
  target->z      = source->z;
  target->p0     = source->p0;
}

