/* Copyright 1994 by Michael Borland and Argonne National Laboratory,
 * all rights reserved.
 */
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
    double Epsi;         /* omega*(t_part-t_first)+phase      */
    double Bpsi;         /* omega*(t_part-t_first)+phase-PI/2 */
    double Estrength;    /* |e.V|/(gap.m.c^2.k) */
    double Bstrength;    /* -|eB|/(m.c.k) */
    double x, xp, y, yp, dxp, dyp;
    double beta, dbeta_r, dr, beta_z, gamma, pc;
    double omega, k, Ephase, kl;
    double cos_tilt, sin_tilt;
    double Estrengthp, cos_Epsi;
    double Bstrengthp, cos_Bpsi;
    double Bkick;
    double length, dt_part;
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
        if (rf_param->n_sections<1)
            rf_param->n_sections = 10;
        rf_param->initialized = 1;
        if (rf_param->gap==0 || rf_param->frequency==0) 
            bomb("RFDF cannot have gap=0 or frequency=0", NULL);
        }

    omega = 2*PI*rf_param->frequency;
    t_first = rf_param->t_first_particle;
    cos_tilt = cos(rf_param->tilt);
    sin_tilt = sin(rf_param->tilt);
    kl = (length=rf_param->length/rf_param->n_sections)*(k = omega/c_mks) ;
    Estrength = (e_mks*rf_param->voltage)/(me_mks*sqr(c_mks)*k*rf_param->gap);
    Bstrength = -e_mks*rf_param->B_field/(me_mks*c_mks*k);
    Ephase = rf_param->phase*PI/180.0 + omega*(rf_param->time_offset-t_first);
    Bkick  = -e_mks*rf_param->B_field/(me_mks); 

    for (ip=0; ip<n_particles; ip++) {
        x  = initial[ip][0] - rf_param->dx;
        xp = initial[ip][1];
        y  = initial[ip][2] - rf_param->dy;
        yp = initial[ip][3];
        pc = pc_central*(1+initial[ip][5]);
        beta = pc/(gamma=sqrt(1+sqr(pc)));
        beta_z = beta/sqrt(1+sqr(xp)+sqr(yp));
        t_part = initial[ip][4]/(c_mks*beta);
        for (is=0; is<rf_param->n_sections; is++) {
            Epsi    = t_part*omega + Ephase;
            Bpsi    = Epsi - PIo2;
            dbeta_r = (Estrengthp=Estrength/gamma)*
                        ((cos_Epsi=cos(Epsi)) - cos(kl/beta_z+Epsi))
                    + (Bstrengthp=Bstrength*beta_z/gamma)*
                        ((cos_Bpsi=cos(Bpsi)) - cos(kl/beta_z+Bpsi));
            dr      = Estrengthp/beta_z*
                       (length*cos_Epsi - (sin(kl/beta_z+Epsi)-sin(Epsi))*beta_z/k)
                    + Bstrengthp/beta_z*
                       (length*cos_Bpsi - (sin(kl/beta_z+Bpsi)-sin(Bpsi))*beta_z/k);
            t_part += (dt_part=length/(c_mks*beta_z)); 
            dxp = dbeta_r*cos_tilt/beta_z;
            dyp = dbeta_r*sin_tilt/beta_z;
            beta_z += 
                (Bkick*(xp*cos_tilt+yp*sin_tilt)*beta_z*dt_part/gamma 
                        - sqr(gamma)*pow3(beta_z)*(xp*dxp+yp*dyp))/
                                                    (1+sqr(beta_z*gamma));
            x      += xp*length + dr*cos_tilt;
            y      += yp*length + dr*sin_tilt;
            xp     += dxp;
            yp     += dyp;
            beta   = beta_z*sqrt(1+sqr(xp)+sqr(yp));
            if (beta<1)
                pc = beta/sqrt(1-sqr(beta));
            gamma   = sqrt(1+sqr(pc));
            }
        final[ip][0] = x + rf_param->dx;
        final[ip][1] = xp;
        final[ip][2] = y + rf_param->dy;
        final[ip][3] = yp;
        final[ip][4] = t_part*c_mks*beta;
        final[ip][5] = (pc-pc_central)/pc_central;
        final[ip][6] = initial[ip][6];
        }
    log_exit("track_through_rf_deflector");
    }

