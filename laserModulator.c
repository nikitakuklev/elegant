/*************************************************************************\
* Copyright (c) 2003 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2003 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: laserModulator.c 
 * Based on a MATLAB script by P. Emma.
 *
 * Michael Borland, 2004
 */

/*
 $Log: not supported by cvs2svn $
 Revision 1.1  2004/03/19 17:27:35  borland
 First version in repository.  May not compile and has had no testing!

*/
 
#include "mdb.h"
#include "track.h"
#include "complex.h"
#include "matlib.h"

/* Comment's from P. Emma's code:
%   Function to laser-modulate e- beam in an undulator.  The undulator is modeled
%   as a series of alternating sign constant bend-radius dipoles, rather than a
%   sinusoidal field (for calc. convenience).  Therefore the K-value for resonance
%   is slightly different.
%
%   INPUTS:     
%               lam:    Laser wavelength [m]
%               phi0:   Laser phase [rad]
%               Nper:   Number of undulator periods (must be even number) [ ]
%               w0:     Laser transverse beam size at waist [m]
%               P0:     Laser peak power [W]
*/

void trackLaserModulator(double **part, long nPart, LSRMDLTR *lsrMdltr, double Po)
{
  double c      = 2.99792458E8;          /* light speed [m/s] */
  double mc2    = 510.99906E-6;          /* e- rest mass [GeV] */
  double Z0     = 120*PI;                /* free-space impedance [ohms] */
  long Nper, harm, np, i, j, pm, jN;
  double gam0, E0, xl0, yl0, ku, lam, x, y, z, t, r2, phi0;
  double per, K, theta0, Lb, ds, R;
  double xpk, k, w, Ef0, ZR, tr, zr, w0, P0, *sr, *theta, *C, *S;
  VMATRIX *Rp, *Rm, *R1;
  COMPLEX Q, Efx, ctmp1, ctmp2, ctmp3;
  COMPLEX Efz, Bfx, Bfy, Bfz;
  double Bf[3], Ef[3], v[3], F[3], dkik[3];
  double tmp1, tmp2, tmp3, v0, vs;
  
  /* force number of periods to be even */
  Nper = 2*(lsrMdltr->undulatorPeriods/2);
  per = lsrMdltr->length/Nper;
  lam = lsrMdltr->laserWavelength;
  w0 = lsrMdltr->laserW0;
  P0 = lsrMdltr->laserPeakPower;
  
  gam0 = sqrt(sqr(Po)+1);
  harm   = 1;                      /*  % laser wavelength slips per und. period (normally = 1) */
  v0     = c*sqrt(1-1/sqr(gam0));  /*  % e- nominal velocity at bend center [m/s] */
  E0     = gam0*mc2;            /*  % e- energy [GeV] */
  xl0    = 0;                   /*  % laser hor. position at bend center [m] */
  yl0    = 0;                   /*  % laser ver. position at bend center [m] */
  ku     = 2*PI/per;           /*  % und wave number [1/m] */
  /* % und param. (for const. radius bend und.) [ ] */
  K      = sqrt(3*(2*sqr(gam0)*harm*lam/per-1)); 
  np     = 10*harm;                         /*  % number of slices per half-period */
  theta0 = K/gam0;                          /*  % bend angle over one-half dipole (quarter-period) */
  Lb     = per/2*(1 + 1./6.*sqr(K/gam0));     /*  % length of arc over an und. half-period [m] */
  ds     = Lb/np;                           /*  % step size along und [m] */ 

  /* % 6x6 linear transport matrix through +bend (half-per) */
  Rp = bend_matrix(ds, 2*theta0/np, theta0/np, theta0/np, 0, 0, 0, 0, 0, 0, 0, 
                   0, 0, 1, 1, 0, 0);
  /* % 6x6 linear transport matrix through -bend (half-per) */
  Rm = bend_matrix(ds, -2*theta0/np, -theta0/np, -theta0/np, 0, 0, 0, 0, 0, 0, 0, 
                   0, 0, 1, 1, 0, 0);
  /* % 6x6 linear transport matrix through half-bend (1/4-per) */
  R1 = bend_matrix(Lb/2, theta0, theta0/2, theta0/2, 0, 0, 0, 0, 0, 0, 0, 
                   0, 0, 1, 1, 0, 0);

  R = Lb/2/theta0;  /* bending radius */
  sr = tmalloc(sizeof(*sr)*(np+1));
  theta  = tmalloc(sizeof(*theta)*(np+1));
  C = tmalloc(sizeof(*C)*(np+1));
  S = tmalloc(sizeof(*S)*(np+1));
  for (i=0; i<=np; i++) {
    sr[i] = i*ds - Lb/2;
    theta[i] = sr[i]/R;
    C[i] = cos(theta[i]);  
    S[i] = sin(theta[i]);  
  }
  
  xpk    = -R*(1-cos(theta0));   /*  % peak x pos. at center of half-per (<0) [m] */
  k      = 2*PI/lam;             /*  % laser wave number [1/m] */
  w      = c*k;                  /*  % laser angular frequency [rad/sec] */
  Ef0    = 1e-9*sqrt(Z0*c/4/PI)*sqrt(16*P0/c/sqr(w0));         /*  % peak E-field [GV/m] */
  ZR     = k/2*sqr(w0);                                        /*  % Raleigh range [m] */
  
  /* % transport bunch through initial und termination */
  track_particles(part, R1, part, nPart);
  
  zr = -per/2*Nper;        /* % starting pos. of z (at center of 1st bend, origin at und center) [m] */
  tr = -Lb*(Nper+0.5)/v0;  /* % start time at start of 1st bend (origin at und center) [sec] */
  pm = 1;                  /* % start with -x at bend center (pm only +1 or -1) */
  for (jN=0; jN<(2*Nper+1); jN++) {    /* % loop over each period */
    fprintf(stderr, "Doing period %ld\n", jN);
    for (j=0; j<np; j++) {            /* % loop within each period */
      t = tr + j*ds/v0;       /* % time advances for each ds step [sec] */
      for (i=0; i<nPart; i++) {         /* loop over all particles */
        x = pm*(xpk + R*(1-C[j])) 
          + part[i][0]*C[j] 
            + part[i][4]*pm*S[j] - xl0;  /* % e- x-coordinates wrt laser (origin at bend center) [m] */
        y = part[i][2] - yl0;            /* % e- y-coordinates wrt laser [m] */
        z = zr + R*S[j] - part[i][0]*pm*S[j] 
          + part[i][4]*C[j];  /* % longitudinal e- z-coordinates wrt bend center [m] */
        r2 = sqr(x) + sqr(y);   /* % radial distance from e- to laser at each s location [m] */
        /* % Alex Chao's complex-Q [m]                */
        /* Q = 1/(z -i*ZR) */
        Q = cdiv(cassign(1, 0), cassign(z, -ZR));
        /* % complex x-E-field [GV/m] */
        /* Efx  = Ef0*exp(-i*w*t+i*k*z+i*phi0+i*k/2*r2.*Q)./(1+i*z/ZR)  */
        ctmp1 = cadd(cassign(0, -w*t+k*z+phi0), cmul(cassign(0, k/2*r2), Q));
        ctmp2 = cdiv(cexp(ctmp1), cassign(1, z/ZR));
        Efx = cmulr(ctmp2, Ef0);
        /* Ef   = [Efx, zeros(size(x)), -Efx.*Q.*x];             % 3D E-field vector due to laser [GV/m] */
        Efz = cmulr(cmul(Efx, Q), -x);
        Ef[0] = Efx.r;
        Ef[1] = 0;
        Ef[2] = Efz.r;
        /* Bf   = real([-Efx.*Q.^2.*x.*y, Efx.*(Q.^2.*x.^2-i*Q/k+1), -Efx.*Q.*y]); */
        Bfx = cmulr(cmul(Efx, cipowr(Q, 2)), -x*y);
        ctmp1 = cmulr(cmul(Q, cassign(0, 1)), -1./k);
        ctmp2 = cadd(cadd(cmulr(cipowr(Q, 2), sqr(x)), ctmp1), cassign(1, 0));
        Bfy = cmul(Efx, ctmp2);
        Bfz = cmulr(cmul(Efx, Q), -y);
        Bf[0] = Bfx.r;
        Bf[1] = Bfy.r;
        Bf[2] = Bfz.r;
        
        vs = c_mks*sqrt(1-1/sqr(gam0)/sqr(1+part[i][5]));       /*  % velocity in s-direction (along arc) [m/s] */
        /* % 3D velocity vector in x,y,z coordinates [m/s] */
        v[0] = vs*pm*S[j] + vs*part[i][1]*C[j];
        v[1] = vs*part[i][3];
        v[2] = vs*C[j];

        F[0] = Ef[0] + (v[1]*Bf[2] - v[2]*Bf[1]);
        F[1] = Ef[1] + (v[2]*Bf[0] - v[0]*Bf[2]);
        F[2] = Ef[2] + (v[0]*Bf[1] - v[1]*Bf[0]);
        /* % calculate kick in x', y', and dE/E at each j [ ] */
        dkik[0] = F[0]*ds/E0;
        dkik[1] = F[1]*ds/E0;
        dkik[2] = F[2]*ds/E0;
        
        /* % add x-kick, rotated to bunch coordinates, to each particle [ ] */
        part[i][1] = part[i][1] + dkik[0]*C[j] - dkik[2]*pm*S[j];   
        /* % add y-kick to each particle [ ] */
        part[i][3] = part[i][3] + dkik[1];
        /* % add z-kick, rotated to bunch coordinates, to each particle [ ]     */
        part[i][5] = part[i][5] + dkik[0]*pm*S[j] + dkik[2]*C[j];   
      }
      if (pm > 0) {
        track_particles(part, Rm, part, nPart);
      } else {
        track_particles(part, Rp, part, nPart);
      }
    }
    pm = -pm;    /* % toggle +1 -> -1 for which side of wiggle we are on */
    zr = zr + per/2;
    tr = tr + Lb/v0;
  }
  track_particles(part, R1, part, nPart);
}
