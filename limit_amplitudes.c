/* Copyright 1994 by Michael Borland and Argonne National Laboratory,
 * all rights reserved.
 */
#include <ctype.h>
#include "track.h"
#include "mdb.h"

/* routine: rectangular_collimator()
 * purpose: eliminate particles that hit the walls of a non-zero length
 *          rectangular collimator
 *
 * Michael Borland, 1989
 */

#ifdef VAX_VMS
#define isnan(x) 0
#define isinf(x) 0
#endif

long rectangular_collimator(
    double **initial, RCOL *rcol, long np, double **accepted, double z,
    double Po
    )
{
    double length, *ini;
    long ip, itop, is_out;
    double xsize, ysize;
    double x_center, y_center;
    double x1, y1, zx, zy;

    log_entry("rectangular_collimator");

    xsize  = rcol->x_max;
    ysize  = rcol->y_max;
    x_center = rcol->dx;
    y_center = rcol->dy;

    if (!(xsize || ysize)) {
        log_exit("rectangular_collimator");
        return(np);
        }

    itop = np-1;
    for (ip=0; ip<np; ip++) {
        ini = initial[ip];
        if ((xsize && fabs(ini[0]-x_center) > xsize) ||
            (ysize && fabs(ini[2]-y_center) > ysize) ||
            isinf(ini[0]) || isinf(ini[2]) ||
            isnan(ini[0]) || isnan(ini[2]) ) {
            SWAP_PTR(initial[ip], initial[itop]);
            if (accepted)
                SWAP_PTR(accepted[ip], accepted[itop]);
            initial[itop][4] = z; /* record position of particle loss */
            initial[itop][5] = Po*(1+initial[itop][5]);
            --itop;
            --ip;
            --np;
            }
        }
    if (np==0 || (length=rcol->length)<=0) {
        log_exit("rectangular_collimator");
        return(np);
        }

    itop = np-1;
    for (ip=0; ip<np; ip++) {
        ini = initial[ip];
        x1 = ini[0] + length*ini[1];
        y1 = ini[2] + length*ini[3];
        is_out = 0;
        if (xsize && fabs(x1-x_center)>xsize)
            is_out = 1;
        if (ysize && fabs(y1-y_center)>ysize)
            is_out += 2;
        if (isinf(x1) || isinf(y1) || isnan(x1) || isnan(y1) )
            is_out += 4;
        if (is_out&4) {
            ini[4] = z+length;
            ini[0] = x1;
            ini[2] = y1;
            ini[5] = Po*(1+ini[5]);
            SWAP_PTR(initial[ip], initial[itop]);
            if (accepted)
                SWAP_PTR(accepted[ip], accepted[itop]);
            --itop;
            --ip;
            --np;
            }
        else if (is_out) {
            zx = zy = DBL_MAX;
            if (is_out&1 && ini[1]!=0)
                zx = (SIGN(ini[1])*xsize-(ini[0]-x_center))/ini[1];
            if (is_out&2 && ini[3]!=0)
                zx = (SIGN(ini[3])*ysize-(ini[2]-y_center))/ini[3];
            if (zx<zy) {
                ini[0] += ini[1]*zx;
                ini[2] += ini[3]*zx;
                ini[4] = z+zx; 
                }
            else {
                ini[0] += ini[1]*zy;
                ini[2] += ini[3]*zy;
                ini[4] = z+zy;
                }
            ini[5] = Po*(1+ini[5]);
            SWAP_PTR(initial[ip], initial[itop]);
            if (accepted)
                SWAP_PTR(accepted[ip], accepted[itop]);
            --itop;
            --ip;
            --np;
            }
        else {
            ini[4] += length*sqrt(1+sqr(ini[1])+sqr(ini[3]));
            ini[0] = x1;
            ini[2] = y1;
            }
        }
    log_exit("rectangular_collimator");
    return(np);
    }


/* routine: limit_amplitudes()
 * purpose: eliminate particles with (x,y) larger than given values
 *
 * Michael Borland, 1989
 */

long limit_amplitudes(
    double **coord, double xmax, double ymax, long np, double **accepted,
    double z, double Po, long extrapolate_z)
{
    long ip, itop, is_out;
    double *part;
    double dz, dzx, dzy;

    log_entry("limit_amplitudes");

    if (!(xmax || ymax)) {
        log_exit("limit_amplitudes");
        return(np);
        }

    itop = np-1;

    for (ip=0; ip<np; ip++) {
        part = coord[ip];
        is_out = 0;
        if (xmax && fabs(part[0])>xmax)
            is_out += 1;
        if (ymax && fabs(part[2])>ymax)
            is_out += 2;
        if (isinf(part[0]) || isinf(part[2]) || isnan(part[0]) || isnan(part[2]) )
            is_out += 4;
        dz = 0;
        if (is_out && !(is_out&4) && extrapolate_z) {
            /* find the actual position of loss, assuming a drift preceded with 
             * the same aperture 
             */
            dzx = dzy = -DBL_MAX;
            if (is_out&1 && part[1]!=0)
                dzx = (part[0]-SIGN(part[1])*xmax)/part[1];
            if (is_out&2 && part[3]!=0)
                dzy = (part[2]-SIGN(part[3])*ymax)/part[3];
            if (dzx>dzy) 
                dz = -dzx;
            else
                dz = -dzy;
            if (dz==-DBL_MAX)
                dz = 0;
            part[0] += dz*part[1];
            part[2] += dz*part[3];
            }
        if (is_out) {
            SWAP_PTR(coord[ip], coord[itop]);
            coord[itop][4] = z+dz;  /* record position of loss */
            coord[itop][5] = Po*(1+coord[itop][5]);
            if (accepted)
                SWAP_PTR(accepted[ip], accepted[itop]);
            --itop;
            --ip;
            np--;
            }
        }
    log_exit("limit_amplitudes");
    return(np);
    }

            
/* routine: elliptical_collimator()
 * purpose: eliminate particles that hit the walls of a non-zero length
 *          elliptical collimator
 *
 * Michael Borland, 1989
 */

long elliptical_collimator(
    double **initial, ECOL *ecol, long np, double **accepted, double z,
    double Po)
{
    double length, *ini;
    long ip, itop;
    double a2, b2;
    double dx, dy;

    log_entry("elliptical_collimator");

    a2 = sqr(ecol->x_max);
    b2 = sqr(ecol->y_max);
    dx = ecol->dx;
    dy = ecol->dy;

    if (!a2 || !b2) {
        log_exit("elliptical_collimator");
        return(np);
        }

    itop = np-1;
    for (ip=0; ip<np; ip++) {
        ini = initial[ip];
        if ((sqr(ini[0]-dx)/a2 + sqr(ini[2]-dy)/b2)>1 ||
            isinf(ini[0]) || isinf(ini[2]) ||
            isnan(ini[0]) || isnan(ini[2]) ) {
            SWAP_PTR(initial[ip], initial[itop]);
            initial[itop][4] = z;
            initial[itop][5] = sqrt(sqr(Po*(1+initial[itop][5]))+1);
            if (accepted)
                SWAP_PTR(accepted[ip], accepted[itop]);
            --itop;
            --ip;
            np--;
            }
        }
    if (np==0 || (length=ecol->length)<=0) {
        log_exit("elliptical_collimator");
        return(np);
        }
    z += length;

    itop = np-1;
    for (ip=0; ip<np; ip++) {
        ini = initial[ip];
        ini[0] += ini[1]*length;
        ini[2] += ini[3]*length;
        if ((sqr(ini[0]-dx)/a2 + sqr(ini[2]-dy)/b2)>1 ||
            isnan(ini[0]-dx) || isnan(ini[2]-dy)) {
            SWAP_PTR(initial[ip], initial[itop]);
            initial[itop][4] = z;
            initial[itop][5] = sqrt(sqr(Po*(1+initial[itop][5]))+1);
            if (accepted)
                SWAP_PTR(accepted[ip], accepted[itop]);
            --itop;
            --ip;
            np--;
            }
        else
            ini[4] += length*sqrt(1+sqr(ini[1])+sqr(ini[3]));
        }
    log_exit("elliptical_collimator");
    return(np);
    }


/* routine: elimit_amplitudes()
 * purpose: eliminate particles outside an ellipse with given semi-major
 *          and semi-minor axes.
 *
 * Michael Borland, 1989
 */

long elimit_amplitudes(
    double **coord,
    double xmax,       /* half-axis in x direction */ 
    double ymax,       /* half-axis in y direction */
    long np,
    double **accepted,
    double z,
    double Po,
    long extrapolate_z
    )
{
    long ip, itop;
    double *part;
    double a2, b2, c1, c2, c0, dz, det;

    log_entry("elimit_amplitudes");

    if (!(xmax || ymax)) {
        log_exit("elimit_amplitudes");
        return(np);
        }

    a2 = sqr(xmax);
    b2 = sqr(ymax);
    itop = np-1;

    for (ip=0; ip<np; ip++) {
        part = coord[ip];
        if (isinf(part[0]) || isinf(part[2]) || isnan(part[0]) || isnan(part[2]) ) {
            SWAP_PTR(coord[ip], coord[itop]);
            coord[itop][4] = z;
            coord[itop][5] = sqrt(sqr(Po*(1+coord[itop][5]))+1);
            if (accepted)
                SWAP_PTR(accepted[ip], accepted[itop]);
            --itop;
            --ip;
            np--;
            continue;
            }
        if ((sqr(part[0])/a2 + sqr(part[2])/b2) > 1) {
            c0 = sqr(part[0])/a2 + sqr(part[2])/b2 - 1;
            c1 = 2*(part[0]*part[1]/a2 + part[2]*part[3]/b2);
            c2 = sqr(part[1])/a2 + sqr(part[3])/b2;
            det = sqr(c1)-4*c0*c2;
            dz = 0;
            if (z>0 && extrapolate_z && c2 && det>=0) {
                if ((dz = (-c1+sqrt(det))/(2*c2))>0)
                    dz = (-c1-sqrt(det))/(2*c2);
                if ((z+dz)<0) 
                    dz = -z;
                part[0] += dz*part[1];
                part[2] += dz*part[3];
                }
            SWAP_PTR(coord[ip], coord[itop]);
            coord[itop][4] = z+dz;
            coord[itop][5] = sqrt(sqr(Po*(1+coord[itop][5]))+1);
            if (accepted)
                SWAP_PTR(accepted[ip], accepted[itop]);
            --itop;
            --ip;
            np--;
            }
        }
    log_exit("elimit_amplitudes");
    return(np);
    }

#define SQRT_3 1.7320508075688772
            
long beam_scraper(
    double **initial, SCRAPER *scraper, long np, double **accepted, double z,
    double Po
    )
{
    double length, *ini;
    long do_x, do_y, ip, itop;
    double limit;

    log_entry("beam_scraper");

    if (scraper->insert_from) {
        switch (toupper(scraper->insert_from[1])) {
          case 'Y': case 'V':
            scraper->direction = 1;
            break;
          case 'X': case 'H':
            scraper->direction = 0;
            break;
          default:
            fprintf(stdout, "Error: invalid scraper insert_from parameter: %s\n",
                    scraper->insert_from);
            fflush(stdout);
            bomb("scraper insert_from axis letter is not one of x, h, y, or v", NULL);
            break;
            }
        if (scraper->insert_from[0]=='-')
            scraper->direction += 2;
        scraper->insert_from = NULL;
        }
    if (scraper->direction<0 || scraper->direction>3)
      return np;
    
    if (scraper->direction==0 || scraper->direction==2) {
        do_x = scraper->direction==0 ? 1 : -1;
        do_y = 0;
        limit = scraper->position*do_x;
        }
    else {
        do_x = 0;
        do_y = scraper->direction==1 ? 1 : -1;
        limit = scraper->position*do_y;
        }    

    if (scraper->length && scraper->Xo) {
        /* scraper has material properties that scatter beam and
         * absorb energy--method is the same as track_through_matter()
         */
        double Nrad, theta_rms, beta, P, E;
        double x1, y1, z1, z2, dx, dy, ds, t;
        
        beta = Po/sqrt(sqr(Po)+1);
        Nrad = scraper->length/scraper->Xo;
        theta_rms = 13.6/(me_mev*Po*sqr(beta))*sqrt(Nrad)*(1+0.038*log(Nrad));
        itop = np-1;
        for (ip=0; ip<np; ip++) {
            ini = initial[ip];
            x1  = ini[0] + ini[1]*scraper->length;
            y1  = ini[2] + ini[3]*scraper->length;
            if ((do_x && (do_x*(ini[0]-scraper->dx)>limit || do_x*(x1-scraper->dx)>limit)) ||
                (do_y && (do_y*(ini[2]-scraper->dy)>limit || do_y*(y1-scraper->dy)>limit)) ) {
                /* scatter and/or absorb energy */
                if (!scraper->elastic) {
                    P = (1+ini[5])*Po;
                    E = sqrt(sqr(P)+1)-1;
                    beta = P/E;
                    t = ini[4]/beta;
                    }
                z1 = gauss_rn(0, random_2);
                z2 = gauss_rn(0, random_2);
                ini[0] += (dx=(z1/SQRT_3 + z2)*scraper->length*theta_rms/2 + scraper->length*ini[1]);
                ini[1] += z2*theta_rms;
                z1 = gauss_rn(0, random_2);
                z2 = gauss_rn(0, random_2);
                ini[2] += (dy=(z1/SQRT_3 + z2)*scraper->length*theta_rms/2 + scraper->length*ini[3]);
                ini[3] += z2*theta_rms;
                ds = sqrt(sqr(scraper->length)+sqr(dx)+sqr(dy));
                if (!scraper->elastic) {
                    E *= exp(-ds/scraper->Xo);
                    P = sqrt(sqr(E+1)-1);
                    ini[5] = (P-Po)/Po;
                    beta = P/E;
                    ini[4] = t*beta+ds;
                    }
                else
                    ini[4] += ds;
                }
            else {
                ini[0] = x1;
                ini[2] = y1;
                ini[4] += scraper->length*sqrt(1+sqr(ini[1])+sqr(ini[3]));
                }
            }
        log_exit("beam_scraper");
        return(np);
        }

    itop = np-1;
    for (ip=0; ip<np; ip++) {
        ini = initial[ip];
        if ((do_x && do_x*(ini[0]-scraper->dx) > limit) ||
            (do_y && do_y*(ini[2]-scraper->dy) > limit) ) {
            SWAP_PTR(initial[ip], initial[itop]);
            if (accepted)
                SWAP_PTR(accepted[ip], accepted[itop]);
            initial[itop][4] = z; /* record position of particle loss */
            initial[itop][5] = Po*(1+initial[itop][5]);
            --itop;
            --ip;
            --np;
            }
        }
    if (np==0 || (length=scraper->length)<=0) {
        log_exit("beam_scraper");
        return(np);
        }

    z += length;
    itop = np-1;
    for (ip=0; ip<np; ip++) {
        ini = initial[ip];
        ini[0] += length*ini[1];
        ini[2] += length*ini[3];
        if ((do_x && do_x*(ini[0]-scraper->dx) > limit) ||
            (do_y && do_y*(ini[2]-scraper->dy) > limit) ) {
            SWAP_PTR(initial[ip], initial[itop]);
            initial[itop][4] = z; /* record position of particle loss */
            initial[itop][5] = Po*(1+initial[itop][5]);
            if (accepted)
                SWAP_PTR(accepted[ip], accepted[itop]);
            --itop;
            --ip;
            --np;
            }
        else
            ini[4] += length*sqrt(1+sqr(ini[1])+sqr(ini[3]));
        }
    log_exit("beam_scraper");
    return(np);
    }

long track_through_pfilter(
    double **initial, PFILTER *pfilter, long np, double **accepted, double z,
    double Po
    )
{
  long ip, itop;
  static double *deltaBuffer=NULL;
  static long maxBuffer = 0;
  
  itop = np-1;

  if ((pfilter->lowerFraction || pfilter->upperFraction) && !pfilter->limitsFixed) {
    double level[2]={-1,-1}, limit[2], upper[2];
    long count = 0, i;
    if (maxBuffer<np &&
        !(deltaBuffer=SDDS_Realloc(deltaBuffer, sizeof(*deltaBuffer)*(maxBuffer=np))))
      SDDS_Bomb("memory allocation failure");
    for (ip=0; ip<np; ip++)
      deltaBuffer[ip] = initial[ip][5];
    /* eliminate lowest lowerfraction of particles and highest
       upperfraction */
    if (pfilter->lowerFraction>0 && pfilter->lowerFraction<1) {
      upper[count] = 0;
      level[count++] = pfilter->lowerFraction*100;
    }
    if (pfilter->upperFraction>0 && pfilter->upperFraction<1) {
      upper[count] = 1;
      level[count++] = 100-pfilter->upperFraction*100;
    }
    compute_percentiles(limit, level, count, deltaBuffer, np);
    itop = np-1;
    for (i=0; i<2; i++) {
      if (level[i]<0)
        break;
      if (upper[i]) {
        pfilter->pUpper = (1+limit[i])*Po;
        pfilter->hasUpper = 1;
      }
      else {
        pfilter->pLower = (1+limit[i])*Po;
        pfilter->hasLower = 1;
      }
      if (!pfilter->fixPLimits) {
        /* filter in next block so there are no discrepancies due to
         * small numerical differences
         */
        for (ip=0; ip<=itop; ip++) {
          if ((upper[i] && initial[ip][5]>limit[i]) ||
              (!upper[i] && initial[ip][5]<limit[i])) {
            SWAP_PTR(initial[ip], initial[itop]);
            initial[itop][4] = z;  /* record position of particle loss */
            initial[itop][5] = Po*(1+initial[itop][5]);  /* momentum at loss */
            if (accepted)
              SWAP_PTR(accepted[ip], accepted[itop]);
            --itop;
            --ip;
          }
        }
      }
    }
    if (pfilter->fixPLimits)
      pfilter->limitsFixed = 1;
  }
  
  if (pfilter->limitsFixed) {
    double p;
    for (ip=0; ip<=itop; ip++) {
      p = (1+initial[ip][5])*Po;
      if ((pfilter->hasUpper && p>pfilter->pUpper) ||
          (pfilter->hasLower && p<pfilter->pLower)) {
        SWAP_PTR(initial[ip], initial[itop]);
        initial[itop][4] = z;  /* record position of particle loss */
        initial[itop][5] = p;  /* momentum at loss */
          if (accepted)
            SWAP_PTR(accepted[ip], accepted[itop]);
        --itop;
        --ip;
      }
    }
  }
  
  if (pfilter->deltaLimit<0) {
#if defined(MINIMIZE_MEMORY)
    free(deltaBuffer);
    deltaBuffer = NULL;
    maxBuffer= 0;
#endif
    return itop+1;
  }
  for (ip=0; ip<=itop; ip++) {
    if (fabs(initial[ip][5])<pfilter->deltaLimit)
      continue;
    SWAP_PTR(initial[ip], initial[itop]);
    initial[itop][4] = z;  /* record position of particle loss */
    initial[itop][5] = Po*(1+initial[itop][5]);  /* momentum at loss */
    if (accepted)
      SWAP_PTR(accepted[ip], accepted[itop]);
    --itop;
    --ip;
  }
#if defined(MINIMIZE_MEMORY)
  free(deltaBuffer);
  deltaBuffer = NULL;
  maxBuffer= 0;
#endif
  return itop+1;
}

