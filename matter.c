/* Copyright 1994 by Michael Borland and Argonne National Laboratory,
 * all rights reserved.
 */
/* file: matter.c
 * contents: track_through_matter()
 *
 * Michael Borland, 1993
 */
#include "mdb.h"
#include "track.h"

#define SQRT_3 1.7320508075688772

void track_through_matter(
    double **part, long np, MATTER *matter, double Po
    )
{
    long ip;
    double L, Nrad, *coord, theta_rms, beta, P, E;
    double z1, z2, dx, dy, ds, t;

    log_entry("track_through_matter");

    if ((L=matter->length)==0)
        return;

    if (matter->Xo==0)
        Nrad = theta_rms = 0;
    else {
        Nrad = matter->length/matter->Xo;
        beta = Po/sqrt(sqr(Po)+1);
        theta_rms = 13.6/me_mev/Po/sqr(beta)*sqrt(Nrad)*(1+0.038*log(Nrad));
        }

    for (ip=0; ip<np; ip++) {
        coord = part[ip];
        if (Nrad) {
            if (!matter->elastic) {
                P = (1+coord[5])*Po;
                E = sqrt(sqr(P)+1)-1;
                beta = P/E;
                t = coord[4]/beta;
                }
            z1 = gauss_rn(0, random_2);
            z2 = gauss_rn(0, random_2);
            coord[0] += (dx=(z1/SQRT_3 + z2)*L*theta_rms/2 + L*coord[1]);
            coord[1] += z2*theta_rms;
            z1 = gauss_rn(0, random_2);
            z2 = gauss_rn(0, random_2);
            coord[2] += (dy=(z1/SQRT_3 + z2)*L*theta_rms/2 + L*coord[3]);
            coord[3] += z2*theta_rms;
            ds = sqrt(sqr(L)+sqr(dx)+sqr(dy));
            if (!matter->elastic) {
                E *= exp(-ds/matter->Xo);
                P = sqrt(sqr(E+1)-1);
                coord[5] = (P-Po)/Po;
                beta = P/E;
                coord[4] = t*beta+ds;
                }
            else
                coord[4] += ds;
            }
        else {
            coord[0] += L*coord[1];
            coord[2] += L*coord[3];
            coord[4] += L*sqrt(1+sqr(coord[1])+sqr(coord[3]));
            }
        }

    log_exit("track_through_matter");
    }
