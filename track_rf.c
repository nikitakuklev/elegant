/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* routine: track_through_rf_deflector()
 * purpose: track particles through an RF deflector
 * 
 * Michael Borland, 1989
 */
#include "track.h"
#include "mdb.h"

void track_through_rf_deflector(
    double **final, 
    RFDF *rf_param,
    double **initial,
    long n_particles,
    double pc_central
    )
{
    double t_first;     /* time when first particle crosses cavity center */
    double t_part;      /* time at which a particle enters cavity */
    double Estrength;    /* |e.V.L/nSections|/(gap.m.c^2.k) */
    double x, xp, y, yp;
    double beta, dp_r, px, py, pz, beta_z, pc;
    double omega, Ephase;
    double cos_tilt, sin_tilt;
    double length, k, z;
    long ip, is;

    log_entry("track_through_rf_deflector");

    if (!rf_param->initialized) {
        rf_param->initialized = 1;
        for (ip=rf_param->t_first_particle=0; ip<n_particles; ip++) {
            pc = pc_central*(1+initial[ip][5]);
            beta = pc/sqrt(1+sqr(pc));
            beta_z = beta/(1 + sqr(initial[ip][1])+sqr(initial[ip][3]));
            rf_param->t_first_particle += (initial[ip][4]/beta + 
                                      rf_param->length/(2*beta_z))/c_mks;
            }
        if (n_particles)
            rf_param->t_first_particle /= n_particles;
        if (rf_param->n_kicks<1)
            rf_param->n_kicks = 10;
        rf_param->initialized = 1;
        if (rf_param->frequency==0) 
            bomb("RFDF cannot have frequency=0", NULL);
        }

    omega = 2*PI*rf_param->frequency;
    t_first = rf_param->t_first_particle;
    cos_tilt = cos(rf_param->tilt);
    sin_tilt = sin(rf_param->tilt);
    length = rf_param->length/rf_param->n_kicks;
    Estrength = (e_mks*rf_param->voltage/rf_param->n_kicks)/(me_mks*sqr(c_mks));
    k = omega/c_mks;
    Ephase = rf_param->phase*PI/180.0 + k*rf_param->length/2 + omega*(rf_param->time_offset-t_first);
    
    for (ip=0; ip<n_particles; ip++) {
        x  = initial[ip][0];
        xp = initial[ip][1];
        y  = initial[ip][2];
        yp = initial[ip][3];
        pc = pc_central*(1+initial[ip][5]);
        pz = pc/sqrt(1+sqr(xp)+sqr(yp));
        px = xp*pz;
        py = yp*pz;
        beta = pc/sqrt(1+sqr(pc));
        beta_z = pz/pc;
        t_part = initial[ip][4]/(c_mks*beta);
        z = 0;
        for (is=0; is<=rf_param->n_kicks; is++) {
          if (is==0 || is==rf_param->n_kicks) {
            /* first half-drift and last half-drift */
            t_part += (length/(2*c_mks*beta_z)); 
            x += xp*length/2;
            y += yp*length/2;
            z += length/2;
            if (is==rf_param->n_kicks)
              break;
          } else {
            t_part += (length/(c_mks*beta_z)); 
            x += xp*length;
            y += yp*length;
            z += length;
          }
          dp_r = Estrength*cos(t_part*omega - k*z + Ephase);
          xp = (px += dp_r*cos_tilt)/pz;
          yp = (py += dp_r*sin_tilt)/pz;
          beta_z = pz/pc;
          pc = sqrt(sqr(px)+sqr(py)+sqr(pz));
        }
        beta = pc/sqrt(1+sqr(pc));
        final[ip][0] = x;
        final[ip][1] = xp;
        final[ip][2] = y;
        final[ip][3] = yp;
        final[ip][4] = t_part*c_mks*beta;
        final[ip][5] = (pc-pc_central)/pc_central;
        final[ip][6] = initial[ip][6];
        }
    log_exit("track_through_rf_deflector");
    }

