#include <stdio.h>
#include "mdb.h"
#include "constants.h"
#include "zibs.h"

/*
  converted from ZAP
  some coordinate convention:
  Zap uses coordinate x, z, l for our x, y, z.
  So I had to change several variable names to be consistent
  with our coordinate system. Here is the complete list to preserve my
  sanity:
  ZAP        here
  ----------------
  ipr     -> printFlag
  ns      -> steps
  istep      not used
  itit1      not used
  const   -> const
  dencon  -> dencon
  pnb     -> particles
  coulog  -> coulombLog
  ro      -> re_mks
  beta    -> beta (v/c)
  gm      -> gamma
  clite   -> c_mks
  epsx    -> emitx
  epsz    -> emity
  sigp    -> sigmaDelta
  sigl    -> sigmaz
  epszx   -> emityxRatio
  epsxz   -> emitxyRatio
  imax    -> elements (number of elements)
  dels    -> deltas
  ss      -> ss (magnet position)
  btx     -> betax pointer
  btz     -> betay pointer
  betrxz  -> betaxyRatio
  betrzx  -> betayxRatio
  dqcon   -> tuneShiftConst (space charge effect)
  dqscx   -> SCtuneShiftx
  dqscz   -> SCtuneShifty
  dqscmx  -> SCtuneShiftMax
  epsx0   -> emitx0 (natural emittance)
  epsz0   -> emity0
  tavx    -> tauxAve
  tavy    -> tauyAve
  tavz    -> tauzAve
  iz      -> i (element counter)
  imx     -> elements
  apx     -> alphax pointer
  ap      -> etax pointer
  app     -> etaxp pointer
  betax   -> btx (local variable)
  betaz   -> bty
  alx     -> alx
  eta     -> eta
  etap    -> etap
  am      -> atomicNumber
  cl      -> cz (constants used in scattering rate)
  cz      -> cy
  conl    -> conz
  conz    -> cony
  zintx   -> zintx
  zintz   -> zinty
  zintl   -> zintz
  ccz     -> cy
  iiz     -> k  (counter within steps loop)
  al      -> al (used to store powers of 10, whatever!)
  bl      -> bl 
  maxdec  -> maxDecades
  imax1   -> not used (this is part of a sleazy way zap uses the last
             element of an array to store the average values -- I couldn't
             recognize what was going on in this subroutine until
             I read the input routines and found out what the last element
             of the array was for.)
  imax2   -> not used (same stuff as above)
  fjoe    -> weight (used in summing contribution from each lattice point)
  gsr     -> transSRdampRate
  gsrl    -> longSRdampRate
*/

/* prototypes */
void IBSGrowthRates (double gamma, double emitx, double emity,
                     double sigmaDelta, double sigmaz,
                     double particles,
                     double emitx0, double sigmaDelta0, 
                     double transSRdampRate, double longSRdampRate,
                     double coupling,
                     double *s, double *betax, double *alphax, double *betay, 
                     double *alphay, double *etax, double *etaxp, long elements, 
                     long superperiods, long verbosity,
                     double *xGrowthRate, double *yGrowthRate, double *zGrowthRate) {
/*
     Subroutine uses simpson's rule integration
     to calculate bjorken/mtingwa integrals (eqn. 3.4)
     particle accelerators 13, 115(1983)

     Equivalent expressions are found in conte/martini
     particle accelerators 17, 1(1985)

     Integrals are broken into decades to optimize speed

     Use weighted average of lattice parameters in regions of non-zero eta
         to predict average lifetimes (in last line of integral outputs)
     The weighted lifetimes are based on an average value for script-h
         and are corrected to reflect the entire ring size

*/
  double test = 1e-5;
  double SCtuneShiftLimit = 0.25;
  double simpsonCoeff[2] = {2.,4.};
  long i, j, k;
  double beta, EMeV, coulombLogReturn, dencon;
  double emityxRatio, emitxyRatio, betaxyRatio, betayxRatio;
  double sumx, sumy, delta;
  double constant, tuneShiftConst, SCtuneShiftx, SCtuneShifty, SCtuneShiftMax, particlesSC;
  double tauxAve, tauyAve, tauzAve;
  double betaxAve, betayAve;
  double betx, alx, bety, eta, etap, phi, c1, c2, c3, cx, cy, cz, r1, a, b;
  double cscale, testLog, checkLog, cprime, conz, cony;
  double zintx, zinty, zintz;
  double ccy, td1, td2, tz1, tz2, ty1, ty2, tx1, tx2;
  double h, aloop, term, func, polyx, polyy, polyz, sumz;
  double alam, cof, f, coff, tmpx, tmpy, tmpz;
  double txi, tyi, tzi, weight, taux, tauy, tauz;
  double epsCheck, epzCheck;  

#define STEPS 10
#define MAXDECADES 30
  long maxDecades = MAXDECADES;
  long steps = STEPS; /* number of integration steps per decade */
  double al[MAXDECADES+2], bl[MAXDECADES+2], sqrt_al[MAXDECADES+2];
  double alam2d[MAXDECADES+2][STEPS+2], sqrtAlam2d[MAXDECADES+2][STEPS+2];
  long nsteps, noWarning;
    
  EMeV = sqrt(sqr(gamma) + 1) * me_mev;
  betaxAve = 0.0;
  betayAve = 0.0;
  for (i=1; i<elements; i++) {
    delta = s[i] - s[i-1];
    betaxAve += delta * (betax[i] + betax[i-1])/ 2.;
    betayAve += delta * (betay[i] + betay[i-1])/ 2.;
  }
  betaxAve /= s[elements-1];
  betayAve /= s[elements-1];
  beta = sqrt(1 - 1 / sqr( gamma ));
  noWarning = 1;
  coulombLogReturn = coulombLog(gamma, emitx, emity, betaxAve, betayAve, sigmaz, particles, noWarning);
  if (verbosity>3)
    fprintf( stdout, "Coulomb log: %g.\n", coulombLogReturn);
  constant = particles * coulombLogReturn * sqr(re_mks) * c_mks /
    (8 * PI * pow(beta,3) * pow(gamma,4) * emitx *  emity * sigmaDelta * sigmaz);
  dencon = particles/ (8 * sigmaz * sqrt(pow(PI,3) * emitx * emity));
  emityxRatio = sqrt( emity/ emitx );
  emitxyRatio = 1.0 / emityxRatio;
  sumx = 0.0;
  sumy = 0.0;
  
  if (elements < 2)
    bomb(NULL,"There are fewer than two elements in the twiss function arrays.\n");
  
  for (i=1; i<elements; i++) {
    delta = s[i] - s[i-1];
    betaxyRatio = sqrt(betax[i-1]/betay[i-1]);
    betayxRatio = 1.0/ betaxyRatio;
    sumx += 0.5 * betaxyRatio * delta/ (1 + emitxyRatio * betaxyRatio);
    sumy += 0.5 * betayxRatio * delta/ (1 + emityxRatio * betayxRatio);
    betaxyRatio = sqrt(betax[i]/betay[i]);
    betayxRatio = 1.0/ betaxyRatio;
    sumx += 0.5 * betaxyRatio * delta/ (1 + emitxyRatio * betaxyRatio);
    sumy += 0.5 * betayxRatio * delta/ (1 + emityxRatio * betayxRatio);
  }
  /* formula simplified from zap since we use the superperiods variable */
  tuneShiftConst = re_mks * superperiods * dencon/ pow(gamma,3);
  SCtuneShiftx = tuneShiftConst * sumx;
  SCtuneShifty = tuneShiftConst * sumy;
  SCtuneShiftMax = MAX(SCtuneShiftx, SCtuneShifty);
  if (verbosity>3)
    fprintf( stdout, "Space charge tune shifts:\n x: %g y: %g.\n", 
            SCtuneShiftx, SCtuneShifty);
  if (SCtuneShiftMax > SCtuneShiftLimit ) {
    particlesSC = SCtuneShiftLimit * particles/ SCtuneShiftMax;
    if (verbosity > 0)
      fprintf( stdout, "Warning: Space charge tune shift is excessive - should use fewer "
              "than %12.4g particles per bunch.\n", particlesSC);
  }

  tauxAve = 0.0;
  tauyAve = 0.0;
  tauzAve = 0.0;
  *xGrowthRate= 0.0; 
  *yGrowthRate= 0.0; 
  *zGrowthRate= 0.0; 

  steps = steps + steps%2;
  al[0] = 0;
  for (j=0; j<maxDecades; j++) {
    bl[j] = pow(10,j);
    al[j+1] = bl[j];
    sqrt_al[j] = sqrt(al[j]);
    h = (bl[j]-al[j])/steps;
    for (k=1; k<=steps; k++) {
      alam2d[j][k] = al[j]+k*h;
      sqrtAlam2d[j][k] = sqrt(alam2d[j][k]);
    }
  }
  
  for( i=0; i<elements; i++) {
    /* Does a better job at trapeze integration than original zap 
     */
    if (i==0) {
      weight = FABS(s[1]-s[0])/ 2.0/ s[elements-1];
    } else if (i==(elements-1) ) {
      weight = FABS(s[elements-1]-s[elements-2])/ 2.0/ s[elements-1];
    } else {
      weight = FABS(s[i+1]-s[i-1])/ 2.0/ s[elements-1];
    }
    /* original zap weighting 
    if (i==0)
      weight = 0.0;
    else
      weight = FABS(s[i]-s[i-1])/ s[elements-1];
      */
    if (!weight)
      continue;
    
    betx = betax[i];
    alx = alphax[i];
    bety = betay[i];
    eta = etax[i];
    etap = etaxp[i];
    
    phi=etap+(alx*eta/betx);
    c1=sqr(gamma*eta)/(emitx*betx);
    c3=betx/emitx;
    c2=c3*sqr(gamma*phi);
    cx=c1+c2;
    cz=sqr(gamma/sigmaDelta);
    cy=bety/emity;
    r1=3./cy;
    a=cx+cz;
    b=(c3+cy)*(c1+cz)+cy*c2;
    
    /*
      C=C3*CZ*(C1+CL)
      
      Define cprime=c*cscale to try to keep the value
      small enough for the vax in single precision.
      Test log(c) to see if it needs scaling
      */
    
    cscale = 1.0;
    testLog = 33.0;
    checkLog = log10(c3)+log10(cy)+log10(c1+cz);
    if( checkLog > testLog)
      cscale = pow(10.0,(testLog - checkLog));
    cprime = c3*cy*cscale*(c1+cz);
    
    conz = constant*cz;
    cony = constant*cy;

    if (verbosity>3) {
      fprintf( stdout, "constant= %12.5g  conz=%12.5g  cony=%12.5g"
              "   a=%12.5g b=%12.5g\n", constant, conz, cony, a, b);
    }
    
    /*     split integral into decades, with "steps" steps per decade
     */

    zintx = 0.0;
    zinty = 0.0;
    zintz = 0.0;
    
    /*     Constants for integration loop
           to keep the numbers reasonable, the numerator is
           scaled by 1/cprime and the denominator by 1/cprime**2
           the extra factor of cprime is accounted for after integrating
           */    
    /*     ccz=c**(-2./3.) */
    ccy = pow(cprime,(-2./3.));
    td1=(a+c3)*ccy;
    td2=1./(sqrt(ccy)*cscale*cy);
    tz1=(2.*a-cy-c3)/cprime;
    tz2=(b-2.*c3*cy)/cprime;
    ty1=(-a-c3+2.*cy)/cprime;
    /*     ty2=(b-r1*c+cy*c3)/cprime */
    ty2=(b+cy*c3)/cprime-r1/cscale;
    tx1=(2.*a*(cx-c3)-cy*cx-c3*(cy-cz-2.*c3-6.*c2))/cprime;
    tx2=(c3+cx)*((b+c3*cy)/cprime)-6./cscale+3.*c3*cy*(cz/cprime);
    al[0]= 0.0;
    for( j=0; j<maxDecades; j++ ) {
      /* 
        It seems that al and bl are integration limits. The
        integration intervals seem to be [0,10], [10,100], [100,1000],
        and so on. j is the index over these intervals.
        */
      /* 
        bl[j] = pow(10,j); 
         al[j+1] = bl[j];
      */
      h = (bl[j]-al[j])/ steps;
      aloop = al[j];

      /* Evaluate simpson's rule summation for one interval
         the integrand is calculated in the loop itself
         */
      term = sqrt((cy+aloop)*ccy) * 
        sqrt(aloop*ccy*aloop+td1*aloop+td2);
      func = sqrt_al[j]/term/term/term;
      polyz = tz1*aloop+tz2;
      polyx = tx1*aloop+tx2;
      polyy = ty1*aloop+ty2;
      /* First point in simpson's integration
       */
      sumz = func*polyz;
      sumx = func*polyx;
      sumy = func*polyy;
      /* split decade into "steps" steps (even number, usually 10) 
         and do simpsons rule. There should be 11 points, ten
         of which are treated below.
         */
      for( k=1; k<=steps; k++) {
        alam = alam2d[j][k];
/* second point (k=1) should have coefficient of 4 */
        cof = simpsonCoeff[k%2]; 
        term = sqrt((cy+alam)*ccy) * sqrt(alam*ccy*alam+td1*alam+td2);
        f = sqrtAlam2d[j][k]/term/term/term;
        coff = cof*f;
        polyz = tz1*alam+tz2;
        polyx = tx1*alam+tx2;
        polyy = ty1*alam+ty2;
        sumz += coff*polyz;
        sumx += coff*polyx;
        sumy += coff*polyy;
      }
      /* This is to compensate for the "2" coefficient as
         the coefficient for the last point should be 1.
         */
      sumz -= f*polyz;
      sumx -= f*polyx;
      sumy -= f*polyy;
      tmpz = (sumz/3.0)*h;
      tmpx = (sumx/3.0)*h;
      tmpy = (sumy/3.0)*h;
      zintz += tmpz;
      zintx += tmpx;
      zinty += tmpy;
      /*
        Test to see if integral has converged
        */
      if( FABS(tmpz/zintz)<test && FABS(tmpx/zintx)<test &&
         FABS(tmpy/zinty)<test ) break;
      if (j == maxDecades) 
        fprintf( stdout, "**Warning** Integral did not converge in %ld decades.\n",maxDecades);
    }
    nsteps = steps * j;

    /* divide answers by cprime to account for scaling */
    txi = constant * (zintx/cprime);
    tyi = cony * (zinty/cprime);
    tzi = conz * (zintz/cprime);

    *xGrowthRate += txi*weight;
    *yGrowthRate += tyi*weight;
    *zGrowthRate += tzi*weight;
  }
  taux = 1.0/ *xGrowthRate;
  tauy = 1.0/ *yGrowthRate;
  tauz = 1.0/ *zGrowthRate;
  
  /* the quantum excitation term in the denominator is
       SR  0    1
      g   e   -----
       x   x  1 + k 
    */
  if (transSRdampRate != 0)
    epsCheck = (transSRdampRate - *xGrowthRate/ (1 + coupling)) * emitx/ (transSRdampRate * emitx0/ (1 + coupling));
  if (longSRdampRate != 0)
    epzCheck = (longSRdampRate - *zGrowthRate) * sqr(sigmaDelta)/ (longSRdampRate * sqr(sigmaDelta0));
  
  /* output average values 
   */
  if (verbosity > 1) {
    fprintf( stdout, "(Weighted) average rates (1/sec): longitudinal= %15.6g"
            "   horizontal= %15.6g   vertical= %15.6g\n", *zGrowthRate, *xGrowthRate, *yGrowthRate);
/*    fprintf( stdout, "(Weighted) average lifetimes (sec): longitudinal= %15.6g"
            "   horizontal= %15.6g   vertical= %15.6g\n", tauz, taux, tauy);
*/
  }
  
  /*
    OUTPUT (GSR-ACUPL*GIBS(EPSHAT))*EPSHAT/(GSR*EPS0X) AND
    (GSRL-GLIBS(EPLHAT))*EPLHAT/(GSRL*EPS0L) RATIOS
    */
  if (verbosity>1 && transSRdampRate!=0.0) {
    fprintf( stdout, "coupling = %8.4f\n(transSRdampRate - IBSGrowthRate(emitxTrial)/(1+coupling)) * emitxTrial/(transSRdampRate*emitx0) = %15.6g\n", coupling, epsCheck);
    fprintf( stdout, "(longSRdampRate - IBSLongGrowthRate(emitzTrial)) * emitzTrial/(longSRdampRate*emitz0)    = %15.6g\n", epzCheck);
  }
  return;
}

double coulombLog (double gamma, double emitx, double emity,
                   double betaxAve, double betayAve, double sigz, double particles,
                   long noWarning) {
  double EMeV, transverseEnergy, tempeV, sigmaxcm, sigmaycm, sigmazcm;
  double volume, density, charge, debyeLength, rmax, rmin, rminClassical, rminQuantum;
  double value;
  long debug = 0;
  
  EMeV = sqrt(sqr(gamma) + 1) * me_mev;
  /*
    Calculate transverse temperature as 2*p*x'
    i.e., assume the transverse energy is temperature/2
    */
  transverseEnergy = 0.5e6 * (gamma * EMeV - me_mev) * (emitx/ betaxAve);
  tempeV = 2 * transverseEnergy;
  /*
    calculate beam volume to get density (in cm**-3) 
    */
  sigmaxcm = 100. * sqrt(emitx*betaxAve);
  sigmaycm = 100. * sqrt(emity*betayAve);
  sigmazcm = 100. * sigz;
  volume = 8.0 * sqrt(pow(PI,3)) * sigmaxcm * sigmaycm * sigmazcm;
  density = particles/ volume;
  
  charge = 1;
  debyeLength = 743.4 * sqrt(tempeV/ density)/ charge;
  rmax = MIN( sigmaxcm, debyeLength);
  if(!noWarning && debyeLength < 2.0 * sigmaxcm )
    fprintf( stdout, "Warning: The beam density probably corresponds to an unreasonably high space charge tune shift in this case. (debyeLength < 2.0 * sigmaxcm)\n");

  /*
    Calculate rmin as larger of classical distance of closest approach
    or quantum mechanical diffraction limit from nuclear radius
    */
  rminClassical = 1.44e-7 * sqr(charge)/ tempeV;
  rminQuantum = 1.9732858e-11/ (2.0 * sqrt(2.0e-6 * transverseEnergy * me_mev));
  rmin = MAX(rminClassical, rminQuantum);
  
  value = log(rmax/ rmin);
  if( !noWarning && (value < 0.0) ) 
    fprintf( stdout, "Warning: Colomb logarithm is less than zero.\n");
  if( debug ) {
    fprintf( stdout, "Coulomb outputs (in cm):\nrminClassical = %10.3g\nrminQuantum = %10.3g\ndebyeLength = %10.3g\nsigmaxcm = %10.3g\n", rminClassical, rminQuantum, debyeLength, sigmaxcm);
  }
  return value;
}
