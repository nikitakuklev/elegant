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
                             TWISS *twiss, RADIATION_INTEGRALS *radIntegrals)
{
    long ip;
    double Fx, Fy, Fdelta, Ddelta, P, t;
    double gamma2Ratio, gammaRatio;
    double Srxp, Sryp, Srdelta;
    double xpEta, ypEta;
    double *part, beta;
    static long first = 1;

    if (SReffects->pRef==0) {
      if (!radIntegrals) {
        bomb("Problem with SREFFECTS element: pRef=0 but no radiation integrals computed.  Use the twiss_output command to compute these.", NULL);
      }
      /* take data from radiation integrals */
      SReffects->pRef = Po;
      SReffects->Jx = radIntegrals->Jx;
      SReffects->Jx = radIntegrals->Jx;
      SReffects->Jy = radIntegrals->Jy;
      SReffects->Jdelta = radIntegrals->Jdelta;
      SReffects->exRef = radIntegrals->ex0/(1+SReffects->coupling);
      SReffects->eyRef = SReffects->exRef*SReffects->coupling;
      SReffects->SdeltaRef = radIntegrals->sigmadelta;
      SReffects->DdeltaRef = -radIntegrals->Uo/me_mev/Po;
/*
      fprintf(stdout, "SReffects set up:\nPref=%g, Jx/y/delta=%g/%g/%g\nex/ey=%g/%g  Sdelta=%g  Ddelta=%g\n",
              SReffects->pRef, SReffects->Jx, SReffects->Jy, SReffects->Jdelta,
              SReffects->exRef, SReffects->eyRef, SReffects->SdeltaRef, SReffects->DdeltaRef);
      fflush(stdout);
*/
    }

    gamma2Ratio = (sqr(Po)+1)/(sqr(SReffects->pRef)+1);
    gammaRatio = sqrt(gamma2Ratio);

    if (SReffects->DdeltaRef>0)
        bomb("DdeltaRef>0 in track_SReffects", NULL);
    /* compute P/Po change per turn due to SR losses at the present momentum */
    Ddelta = SReffects->DdeltaRef*gammaRatio*gamma2Ratio;

    /* RMS values for random added components */
    if (twiss->betax<=0 || twiss->betay<=0)
        bomb("Twiss parameters invalid in track_SReffects", NULL);
    /* damping rates less RF contributions */
    Fx = 1 + (SReffects->Jx - 1)*Ddelta;
    Fy = 1 + (SReffects->Jy - 1)*Ddelta;
    Fdelta = 1 + SReffects->Jdelta*Ddelta;
    if (SReffects->qExcite) {
      Srdelta = sqrt(1-sqr(Fdelta))*SReffects->SdeltaRef*gammaRatio;
      Srxp    = sqrt((1-sqr(1+SReffects->Jx*Ddelta))*SReffects->exRef*gamma2Ratio*(1+sqr(twiss->alphax))/twiss->betax);
      Sryp    = sqrt((1-sqr(1+SReffects->Jy*Ddelta))*SReffects->eyRef*gamma2Ratio*(1+sqr(twiss->alphay))/twiss->betay);
    } else {
      Srdelta = Srxp = Sryp = 0;
    }

    /* damping rates less RF contributions */
    if (!SReffects->damp)
      Fx = Fy = Fdelta = 1;
    if (!SReffects->loss)
      Ddelta = 0;

    if (SReffects->fraction!=1) {
      /* scale with fraction of effect */
      Fx = 1+(Fx-1)*SReffects->fraction;
      Fy = 1+(Fy-1)*SReffects->fraction;
      Fdelta = 1+(Fdelta-1)*SReffects->fraction;
      Ddelta *= SReffects->fraction;
      Srdelta *= SReffects->fraction;
      Srxp *= SReffects->fraction;
      Sryp *= SReffects->fraction;
    }

    if (first) {
        first = 0;
/*
        fprintf(stdout, "Damping/QE constants:\nFx = %e, Fy = %e, Fd = %e\nSrxp = %e, Sryp = %e, Srd = %e\n",
               Fx, Fy, Fdelta, Srxp, Sryp, Srdelta);
        fflush(stdout);
        fprintf(stdout, "Twiss parameters: betax = %e, etapx = %e, betay = %e\n",
               twiss->betax, twiss->etapx, twiss->betay);
        fflush(stdout);
*/
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


