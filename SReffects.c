/* Copyright 1995 by Michael Borland and Argonne National Laboratory,
 * all rights reserved.
 */
/* file: SReffects.c
 * contents: track_SReffects()
 *
 * Michael Borland, 1995
 */
#include "mdb.h"
#include "track.h"

void track_SReffects(double **coord, long np, SREFFECTS *SReffects, double Po, 
                             TWISS *twiss)
{
    long ip;
    double Fx, Fy, Fdelta, Ddelta, P, t;
    double gamma2Ratio, gammaRatio;
    double Srxp, Sryp, Srdelta, rx, ry, rdelta, Sxp, Syp;
    double xpEta, ypEta, ex, ey;
    double *part, beta;
    static long first = 1;

    gamma2Ratio = (sqr(Po)+1)/(sqr(SReffects->pRef)+1);
    gammaRatio = sqrt(gamma2Ratio);

    if (SReffects->DdeltaRef>0)
        bomb("DdeltaRef>0 in track_SReffects", NULL);
    /* compute P/Po change per turn due to SR losses */
    Ddelta = SReffects->DdeltaRef*gammaRatio*gamma2Ratio;

    /* damping rates less RF contributions */
    Fx = 1 + (SReffects->Jx - 1)*Ddelta;
    Fy = 1 + (SReffects->Jy - 1)*Ddelta;
    Fdelta = 1 + SReffects->Jdelta*Ddelta;

    /* RMS values for random added components */
    if (twiss->betax<=0 || twiss->betay<=0)
        bomb("Twiss parameters invalid in track_SReffects", NULL);
    Srdelta = sqrt(1-sqr(Fdelta))*SReffects->SdeltaRef*gammaRatio;
    Srxp    = sqrt((1-sqr(1+SReffects->Jx*Ddelta))*SReffects->exRef*gamma2Ratio*(1+sqr(twiss->alphax))/twiss->betax);
    Sryp    = sqrt((1-sqr(1+SReffects->Jy*Ddelta))*SReffects->eyRef*gamma2Ratio*(1+sqr(twiss->alphay))/twiss->betay);

    if (first) {
        first = 0;
        fprintf(stderr, "Damping/QE constants:\nFx = %e, Fy = %e, Fd = %e\nSrxp = %e, Sryp = %e, Srd = %e\n",
               Fx, Fy, Fdelta, Srxp, Sryp, Srdelta);
        fprintf(stderr, "Twiss parameters: betax = %e, etapx = %e, betay = %e\n",
               twiss->betax, twiss->etapx, twiss->betay);
        }
    

    for (ip=0; ip<np; ip++) {
        part     = coord[ip];
        xpEta    = part[5]*twiss->etapx;
        part[1]  = (part[1] - xpEta)*Fx + Srxp*gauss_rn(0, random_2) + xpEta;
        ypEta    = part[5]*twiss->etapy;
        part[3]  = (part[3] - ypEta)*Fy + Sryp*gauss_rn(0, random_2) + ypEta;
        P = (1+part[5])*Po;
        beta = P/sqrt(sqr(P)+1);
        t = part[4]/beta;
        part[5]  = Ddelta + part[5]*Fdelta + Srdelta*gauss_rn(0, random_2);
        P = (1+part[5])*Po;
        beta = P/sqrt(sqr(P)+1);
        part[4] = t*beta;
        }
    }


