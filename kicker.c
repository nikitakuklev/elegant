/* Copyright 1994 by Michael Borland and Argonne National Laboratory,
 * all rights reserved.
 */
/* file: kicker.c
 * contents: track_through_kicker()
 * 
 * Michael Borland, 1991
 */
#include "mdb.h"
#include "track.h"
#include "table.h"


void track_through_kicker(
    double **part, long np, KICKER *kicker, double p_central, long pass, long default_order
    )
{
    long i, n, ip;
    double time, time_offset, angle, t0, *coord, sum_amp, amplitude, dx, ds;
    double x, xp, y, yp, dp, s, cos_tilt, sin_tilt, curv;
    double theta_i, alpha_i, alpha_f;
    long xmem, ymem;

    if (np<=0)
        return;

    log_entry("track_through_kicker");

    if (!kicker->t_wf)
        set_up_kicker(kicker);

    if (kicker->fire_on_pass>pass) {
        drift_beam(part, np, kicker->length, default_order); 
        log_exit("track_through_kicker");
        return;
        }

    if (!kicker->t_wf || !kicker->amp_wf || !kicker->n_wf ||
        kicker->tmin>=kicker->tmax)        
        bomb("no (valid) waveform data for kicker", NULL);

    if (kicker->phase_reference==0) 
        kicker->phase_reference = unused_phase_reference();

    switch (get_phase_reference(&time_offset, kicker->phase_reference)) {
        case REF_PHASE_RETURNED:
            break;
        case REF_PHASE_NOT_SET:
        case REF_PHASE_NONEXISTENT:
            if (!kicker->fiducial_seen) {
                /* set reference phase so that the center of this bunch goes through
                 * at the desired phase.
                 */
                t0 = 0;
                for (ip=t0=0; ip<np; ip++)
                    t0 += part[ip][4]/(c_mks*beta_from_delta(p_central, part[ip][5]));
                t0 /= np;
                kicker->t_fiducial = t0;
                kicker->fiducial_seen = 1;
                }
            set_phase_reference(kicker->phase_reference, time_offset = kicker->t_fiducial);
            break;
        default:
            bomb("unknown return value from get_phase_reference()", NULL);
            break;
        }

    if (kicker->length<=0) {
        log_exit("track_through_kicker");
        return;
        }

    time_offset += kicker->time_offset;

    cos_tilt = cos(kicker->tilt);
    sin_tilt = sin(kicker->tilt);

    if (kicker->spatial_dependence) {
        xmem = rpn_create_mem("x");
        ymem = rpn_create_mem("y");
        }

    sum_amp = 0;
    for (ip=0; ip<np; ip++) {
        angle  = kicker->angle;
        time = part[ip][4]/(c_mks*beta_from_delta(p_central, part[ip][5])) - time_offset;

        if (time>kicker->tmax || time<kicker->tmin) {
            if (kicker->periodic) {
                n = (time-kicker->tmin)/(kicker->tmax-kicker->tmin);
                time -= n*(kicker->tmax-kicker->tmin);
                }
            }
  
        amplitude = 0;
        if (time>kicker->tmax || time<kicker->tmin)
            angle *= 0;
        else {
            /* should to a binary search here ! */
            for (i=1; i<kicker->n_wf; i++)
                if (kicker->t_wf[i]>=time)
                    break;
            if (i==kicker->n_wf) {
                fprintf(stderr, "error: waveform interpolation problem in track_through_kicker()\n");
                fprintf(stderr, "particle time is %21.15e\nwaveform limits are %21.15e to %21.15e\n",
                        time, kicker->tmin, kicker->tmax);
                exit(1);
                }
            i--;
            angle *= (amplitude=INTERPOLATE(kicker->amp_wf[i], kicker->amp_wf[i+1],
                               kicker->t_wf[i],   kicker->t_wf[i+1], time));
            }
        sum_amp += amplitude;

        coord = part[ip];

        if (!angle) {
            coord[0] += kicker->length*coord[1];
            coord[2] += kicker->length*coord[3];
            coord[4] += kicker->length*sqrt(1+sqr(coord[1])+sqr(coord[3]));
            continue;
            }

        x  =  coord[0]*cos_tilt + coord[2]*sin_tilt;
        y  = -coord[0]*sin_tilt + coord[2]*cos_tilt;
        xp =  coord[1]*cos_tilt + coord[3]*sin_tilt;
        yp = -coord[1]*sin_tilt + coord[3]*cos_tilt;
        s  = coord[4];
        dp = coord[5];

        if (kicker->spatial_dependence) {
            rpn_store(x, xmem);
            rpn_store(y, ymem);
            angle *= rpn(kicker->spatial_dependence);
            if (rpn_check_error()) exit(1);
            rpn_clear();
            }

        curv = sin(-angle)/kicker->length/(1+dp);

        theta_i = atan(xp);
        alpha_i = -theta_i;
        alpha_f = asin(kicker->length*curv + sin(alpha_i));
        if (fabs(curv*kicker->length)<1e-12) {
            x += (dx=kicker->length*xp);
            s += (ds=fabs(kicker->length/cos(alpha_i)));
            }
        else {
            x += (dx=(cos(alpha_f)  - cos(alpha_i))/curv);
            s += (ds=fabs((alpha_f-alpha_i)/curv));
            }
        xp = tan(theta_i - (alpha_f - alpha_i));
        y += kicker->length*yp;
        coord[0] = x*cos_tilt - y*sin_tilt;
        coord[2] = x*sin_tilt + y*cos_tilt;
        coord[1] = xp*cos_tilt - yp*sin_tilt;
        coord[3] = xp*sin_tilt + yp*cos_tilt;
        coord[4] = s;
        }
/*
    if (np)
        fprintf(stderr, "average kicker amplitude = %f\n", sum_amp/np);
 */

    log_exit("track_through_kicker");
    }

void set_up_kicker(KICKER *kicker)
{
    TABLE data;
    long i;

    log_entry("set_up_kicker");

    if (!kicker->waveform)
        bomb("no waveform filename given for kicker", NULL);

    if (!get_table(&data, kicker->waveform, 1, 0))
        bomb("unable to read waveform for kicker", NULL);

    if (data.n_data<=1)
        bomb("kicker waveform contains less than 2 points", NULL);

    kicker->t_wf   = data.c1;
    kicker->amp_wf = data.c2;
    kicker->n_wf   = data.n_data;
    for (i=0; i<kicker->n_wf-1; i++)
        if (kicker->t_wf[i]>kicker->t_wf[i+1])
            bomb("time values not monotonically increasing in kicker waveform", NULL);
    kicker->tmin = kicker->t_wf[0];
    kicker->tmax = kicker->t_wf[kicker->n_wf-1];
    tfree(data.xlab); tfree(data.ylab); tfree(data.title); tfree(data.topline);
    data.xlab = data.ylab = data.title = data.topline = NULL;
    log_exit("set_up_kicker");
    }
