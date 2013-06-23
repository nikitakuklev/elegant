/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: matter.c
 * contents: track_through_matter()
 *
 * Michael Borland, 1993
 */
#include "mdb.h"
#include "track.h"

#define SQRT_3 (1.7320508075688772)
#define AMU (1.6605e-27)
#define SQR_PI (PI*PI)
#define ALPHA (1./137.036)

#define BS_Y0 (1e-8)

double radiationLength(long Z, double A, double rho);
double solveBrehmsstrahlungCDF(double F);

long track_through_matter(
                          double **part, long np, MATTER *matter, double Po, double **accepted, double z0
                          )
{
  long ip;
  double L, Nrad, *coord, theta_rms=0, beta, P, gamma=0.0;
  double z1, z2, dx, dy, ds, t=0.0, dGammaFactor;
  double K1, K2=0.0, sigmaTotal, probScatter=0.0, dgamma;
  double Xo, probBSScatter = 0, probERScatter=0;
  long nScatters=0, i_top, isLost;
  long sections, sections0=1, impulseMode;
  double L1, prob, probBS, probER;  
  long multipleScattering = 0;
  
  log_entry("track_through_matter");

  if (particleIsElectron==0)
    bombElegant("MATTER element doesn't work for particles other than electrons", NULL);
  
  if (matter->length!=0) {
    L = matter->length;
    impulseMode = 0;
  } else if (matter->lEffective!=0) {
    L = matter->lEffective;
    impulseMode = 1;
  }
  else 
    return np;

  if (matter->energyDecay && (matter->nuclearBrehmsstrahlung || matter->electronRecoil))
    bombElegant("ENERGY_DECAY=1 and NUCLEAR_BREHMSSTRAHLUNG=1 or ELECTRON_RECOIL=1 options to MATTER/SCATTER element are mutually exclusive", NULL);

  beta = Po/sqrt(sqr(Po)+1);
  if (matter->Xo==0) {
    if (matter->Z<1 || matter->A<1 || matter->rho==0)
      bombElegant("XO=0 but Z, A, or rho invalid for MATTER element", NULL);
    Xo = radiationLength(matter->Z, matter->A, matter->rho);
    /* printf("Computed radiation length for Z=%ld, A=%le, rho=%le is %le m\n",
           matter->Z, matter->A, matter->rho, Xo);
           */
  } else 
    Xo = matter->Xo;
  
  Nrad = matter->length/Xo;
  dGammaFactor = 1-exp(-Nrad);
  prob = probBS = probER = 0;
  if (Nrad<1e-3 || matter->nuclearBrehmsstrahlung || matter->electronRecoil) {
    if (matter->Z<1 || matter->A<1 || matter->rho<=0)
      bombElegant("MATTER element is too thin---provide Z, A, and rho for single-scattering calculation.", NULL);
    K1 = 4*matter->Z*(matter->Z+1)*sqr(particleRadius/(beta*Po));
    K2 = sqr(pow(matter->Z, 1./3.)*ALPHA/Po);
    sigmaTotal = K1*pow(PI, 3)/(sqr(K2)+K2*SQR_PI);
    probScatter = matter->rho/(AMU*matter->A)*matter->length*sigmaTotal;
    /* printf("K1=%le, K2=%le, mean expected number of scatters is %le\n", K1, K2, probScatter); */
    probBSScatter = 0;
    if (matter->nuclearBrehmsstrahlung) {
      probBSScatter = 4*matter->length/(3*Xo)*(-log(BS_Y0)-(1-BS_Y0)+3./8.*(1-BS_Y0*BS_Y0));
    }
    if (matter->electronRecoil) {
      probERScatter = matter->length*matter->rho/(AMU*matter->A)*PIx2*matter->Z*sqr(re_mks)/Po*(1/BS_Y0-1);
    }
    sections0 = probScatter/matter->pLimit+1;
    Nrad /= sections0;
    multipleScattering = 0;    
    L1 = L/sections0;
    prob = probScatter/sections0;
    probBS = probBSScatter/sections0;
    probER = probERScatter/sections0;
    printf("Sections=%ld, L1 = %le, probIS = %le, probBS = %le, probER = %le\n", sections0, L1, prob, probBS, probER);
  } else {
    multipleScattering = 1;
    theta_rms = 13.6/particleMassMV/Po/sqr(beta)*sqrt(Nrad)*(1+0.038*log(Nrad));
  }
  
  i_top = np-1;
  if (impulseMode)
    L = L1 = 0;
  for (ip=0; ip<=i_top; ip++) {
    coord = part[ip];
    isLost = 0;
    if (Nrad) {
      if (matter->energyDecay || matter->nuclearBrehmsstrahlung) {
        P = (1+coord[5])*Po;
        gamma = sqrt(sqr(P)+1);
        beta = P/gamma;
        t = coord[4]/beta;
      }
      if (multipleScattering) {
        /* use the multiple scattering formula */
        z1 = gauss_rn(0, random_2);
        z2 = gauss_rn(0, random_2);
        coord[0] += (dx=(z1/SQRT_3 + z2)*L*theta_rms/2 + L*coord[1]);
        coord[1] += z2*theta_rms;
        z1 = gauss_rn(0, random_2);
        z2 = gauss_rn(0, random_2);
        coord[2] += (dy=(z1/SQRT_3 + z2)*L*theta_rms/2 + L*coord[3]);
        coord[3] += z2*theta_rms;
        ds = sqrt(sqr(L)+sqr(dx)+sqr(dy));
      } else {
        /* model scattering using the cross section */
        double F, theta, phi, zs, dxp, dyp;
        long is;
        ds = dgamma = 0;
        sections = sections0;
        for (is=0; is<sections0 && !isLost; is++) {
          if (random_2(1)<prob) {
            nScatters ++;
            /* single-scattering computation */
            /* scatter occurs at location 0<=zs<=L */
            zs = L1*random_2(1);
            /* pick a value for CDF and get corresponding angle */
            F = random_2(1);
            theta = sqrt((1-F)*K2*SQR_PI/(K2+F*SQR_PI));
            phi = random_2(1)*PIx2;
            dxp = theta*sin(phi);
            dyp = theta*cos(phi);
            /* advance to location of scattering event */
            ds += zs*sqrt(1+sqr(coord[1])+sqr(coord[3]));
            /* scatter */
            coord[1] += dxp;
            coord[3] += dyp;
            /* advance to end of slice */
            coord[0] += dxp*(L1-zs);
            coord[2] += dyp*(L1-zs);
            ds += (L1-zs)*sqrt(1+sqr(coord[1])+sqr(coord[3]));
          } else {
            ds += L1*sqrt(1+sqr(coord[1])+sqr(coord[3]));
            coord[0] += coord[1]*L1;
            coord[2] += coord[3]*L1;
          }
          if (probBS!=0 && random_2(1)<probBS)
            gamma -= gamma*solveBrehmsstrahlungCDF(random_2(1));
          if (probER!=0 && random_2(1)<probER)
            gamma -= BS_Y0/(1-random_2(1)*(1-BS_Y0));
          if (gamma<=1) {
            isLost = 1;
            break;
          }
        }
      }
      if (!isLost) {
        if (probBSScatter) {
          P = sqrt(sqr(gamma)-1);
          coord[5] = (P-Po)/Po;
          beta = P/gamma;
          coord[4] = t*beta+ds;
        } else if (matter->energyDecay) {
          dgamma = gamma*dGammaFactor;
          if (matter->energyStraggle) {
            double dgamma1;
            /* very simple-minded estimate: StDev(dE) = Mean(dE)/2 */
            while ((dgamma1 = dgamma*(1+0.5*gauss_rn(0, random_2)))<0)
              ;
            dgamma = dgamma1;
          }
          gamma -= dgamma;
          if (gamma<=1) 
            isLost = 1;
          else {
            P = sqrt(sqr(gamma)-1);
            coord[5] = (P-Po)/Po;
            beta = P/gamma;
            coord[4] = t*beta+ds;
          }
        }
        else
          coord[4] += ds;
      }
      if (isLost) {
        swapParticles(part[ip], part[i_top]);
        if (accepted)
          swapParticles(accepted[ip], accepted[i_top]);
        part[i_top][4] = z0+ds;
        part[i_top][5] = 0;
        i_top --;
        ip --;
      }
    }
    else {
      coord[0] += L*coord[1];
      coord[2] += L*coord[3];
      coord[4] += L*sqrt(1+sqr(coord[1])+sqr(coord[3]));
    }
  }
  
  log_exit("track_through_matter");
  return (i_top+1);
}


double inelasticGasScattering(double Z, double gamma, double nL, double P)
{
  double C1, C2, St;
  
  C1 = 16*sqr(re_mks*Z)/(3*137)*log(183/pow(Z, 1./3.))*nL;
  C2 = 16*sqr(re_mks)*Z/(3*137)*nL;
  St = P;
  
  return -(-40*C1 + 81*C2 - 40*C2*log((5*gamma)/2.) 
           + sqrt(160*C2*(25*C1 - 35*C2 + 40*St + 25*C2*log((5*gamma)/2.)) 
                  + ipow(40*C1 - 81*C2 + 40*C2*log((5*gamma)/2.),2)))/(80.*C2);
}


double radiationLength(long Z, double A, double rho)
/* Returns radiation length for electrons in m for given Z and A (in AMUs). See PhysRevD.86.010001, page 329 */ 
{
  double fZ;
  double alpha, a2;
  double Lrad, Lradp;
  
  alpha = e_mks*e_mks/(4*PI*epsilon_o*hbar_mks*c_mks);
  a2 = sqr(alpha*Z);
  fZ = a2*(1/(1+a2) + 0.20206 + a2*(-0.0369 + a2*(0.0083 - a2*0.002)));
  switch (Z) {
  case 1:
    Lrad = 5.31;
    Lradp = 6.144;
    break;
  case 2:
    Lrad = 4.79;
    Lradp = 5.621;
    break;
  case 3:
    Lrad = 4.74;
    Lradp = 5.805;
    break;
  case 4:
    Lrad = 4.71;
    Lradp = 5.924;
    break;
  default:
    Lrad = log(184.15*pow(Z, -1./3.));
    Lradp = log(1194*pow(Z, -2./3.));
    break;
  }    
  /* factor is to convert to meters */
  return 1e-3/(4*alpha*sqr(re_mks)*NAvogadro/A*(Z*Z*(Lrad - fZ) + Z*Lradp))/rho;
}

double solveBrehmsstrahlungCDF(double F)
/* Solve F == G(y)/G(1) where G(y)=(ln(y/y0) - (y-y0) + 3/8*(y^2-y0^2)
 */
{
  static double *lnyTable = NULL, *FTable = NULL;
  static double dy, y0=BS_Y0;
  static long beenWarned = 0;
  double y;
  long nf = 1000, i, code;
  
  if (!FTable) {
    /* make a table of F(y) with nf points */
    double y1;
    lnyTable = tmalloc(sizeof(*lnyTable)*nf);
    FTable = tmalloc(sizeof(*FTable)*nf);
    dy = (1-y0)/(nf-1.);
    for (i=0; i<nf; i++) {
      y1 = y0 + i*dy;
      lnyTable[i] = log(y1/y0);
      FTable[i] = log(y1/y0)-(y1-y0)+3./8.*(sqr(y1)-sqr(y0));
    }
    for (i=0; i<nf; i++)
      FTable[i] /= FTable[nf-1];
  }

  y = interp(lnyTable, FTable, nf, F, 0, 1, &code);
  if (code==1) {
    y = y0*exp(y);
  } else {
    if (!beenWarned)  {
      beenWarned = 1; 
      printf("*** Warning: interpolation problem for brehmsstrahlung.\n");
    }
    y = y0;
  }
  return y;
}

