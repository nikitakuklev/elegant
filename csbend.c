/* Copyright 1994 by Michael Borland and Argonne National Laboratory,
 * all rights reserved.
 */
/* file: csbend.c
 * contents:  track_through_canonical_sbend()
 *
 *
 * Michael Borland, 1991, 1992.
 */
#include "mdb.h"
#include "track.h"

void integrate_csbend_ord2(double *Qf, double *Qi, double s, long n, double rho0, double p0);
void integrate_csbend_ord4(double *Qf, double *Qi, double s, long n, double rho0, double p0);
void exactDrift(double **part, long np, double length);
void convertFromCSBendCoords(double **part, long np, double rho0, 
			     double cos_ttilt, double sin_ttilt, long ctMode);
void convertToCSBendCoords(double **part, long np, double rho0, 
			     double cos_ttilt, double sin_ttilt, long ctMode);
long applyLowPassFilter(double *histogram, long bins, double dx, double start, double end);
long correctDistribution(double *array, long npoints, double desiredSum);

static double Fy_0, Fy_x, Fy_x2, Fy_x3, Fy_x4;
static double Fy_y2, Fy_x_y2, Fy_x2_y2;
static double Fy_y4;

static double Fx_y, Fx_x_y, Fx_x2_y, Fx_x3_y;
static double Fx_y3, Fx_x_y3;

static double rho0, rho_actual, rad_coef, isrConstant;

static long particle_lost;
static double s_lost;

#if !defined(PARALLEL)
/* to avoid problems with HP parallel compiler */
extern unsigned long multipoleKicksDone ;
#endif

long track_through_csbend(double **part, long n_part, CSBEND *csbend, double p_error, double Po, double **accepted,
                          double z_start)
{
  double nh, betah2, gammah3, deltah4;
  double h, h2, h3;
  long i_part, i_top;
  double rho, s, Fx, Fy;
  double x, xp, y, yp, dp, y2, dp0;
  double n, beta, gamma, delta, fse, dp_prime;
  double tilt, etilt, cos_ttilt, sin_ttilt, ttilt;
  double *coord;
  double angle, e1, e2, Kg;
  double psi1, psi2, he1, he2;
  double Qi[6], Qf[6];
  double dcoord_etilt[6];
  double dxi, dyi, dzi;
  double dxf, dyf, dzf;
  double delta_xp;
  double e1_kick_limit, e2_kick_limit;
  
  if (!csbend)
    bomb("null CSBEND pointer (track_through_csbend)", NULL);

  if (csbend->angle==0) {
    exactDrift(part, n_part, csbend->length);
    return n_part;
  }

  if (csbend->integration_order!=2 && csbend->integration_order!=4)
    bomb("CSBEND integration_order is invalid--must be either 2 or 4", NULL);

  if (csbend->use_bn) {
    rho0 = csbend->length/csbend->angle;
    csbend->k1_internal = csbend->b1/rho0;
    csbend->k2_internal = csbend->b2/rho0;
    csbend->k3_internal = csbend->b3/rho0;
    csbend->k4_internal = csbend->b4/rho0;
  } else {
    csbend->k1_internal = csbend->k1;
    csbend->k2_internal = csbend->k2;
    csbend->k3_internal = csbend->k3;
    csbend->k4_internal = csbend->k4;
  }

  he1 = csbend->h1;
  he2 = csbend->h2;
  if (csbend->angle<0) {
    angle = -csbend->angle;
    e1    = -csbend->e1;
    e2    = -csbend->e2;
    etilt = csbend->etilt;
    tilt  = csbend->tilt + PI;      /* work in rotated system */
    rho0  = -csbend->length/angle;  /* temporarily keep the sign */
    n     = -sqr(rho0)*csbend->k1_internal;
    beta  = 0.5*csbend->k2_internal*pow3(rho0);
    gamma = csbend->k3_internal*pow4(rho0)/6.;
    delta = csbend->k4_internal*pow5(rho0)/24.;
    /* this is always a postive value now */
    rho0  = -rho0;
  }
  else {
    angle = csbend->angle;
    e1    = csbend->e1;
    e2    = csbend->e2;
    etilt = csbend->etilt;
    tilt  = csbend->tilt;
    rho0  = csbend->length/angle;
    n     = -sqr(rho0)*csbend->k1_internal;
    beta  = 0.5*csbend->k2_internal*pow3(rho0);
    gamma = csbend->k3_internal*pow4(rho0)/6.;
    delta = csbend->k4_internal*pow5(rho0)/24.;
  }

  fse = csbend->fse;
  h2 = sqr(h=1./rho0);
  h3 = h*h2;
  nh = n*h;
  betah2 = beta*h2;
  gammah3 = gamma*h3;
  deltah4 = delta*h2*h2;
  if (fse>-1)
    rho_actual = 1/((1+fse)*h);
  else
    rho_actual = 1e16/h;

  e1_kick_limit = csbend->edge1_kick_limit;
  e2_kick_limit = csbend->edge2_kick_limit;
  if (csbend->kick_limit_scaling) {
    e1_kick_limit *= rho0/rho_actual;
    e2_kick_limit *= rho0/rho_actual;
  }
  if (e1_kick_limit>0 || e2_kick_limit>0)
    fprintf(stdout, "rho0=%e  rho_a=%e fse=%e e1_kick_limit=%e e2_kick_limit=%e\n",
            rho0, rho_actual, csbend->fse, e1_kick_limit, e2_kick_limit);
    fflush(stdout);
  
  /* angles for fringe-field effects */
  Kg   = 2*csbend->hgap*csbend->fint;
  psi1 = Kg/rho_actual/cos(e1)*(1+sqr(sin(e1)));
  psi2 = Kg/rho_actual/cos(e2)*(1+sqr(sin(e2)));

  /* rad_coef is d((P-Po)/Po)/ds for the on-axis, on-momentum particle, where po is the momentum of
   * the central particle.
   */
  if (csbend->synch_rad)
    rad_coef = sqr(e_mks)*pow3(Po)*sqr(1+fse)/(6*PI*epsilon_o*sqr(c_mks)*me_mks*sqr(rho0));
  else
    rad_coef = 0;
  /* isrConstant is the RMS increase in dP/P per meter due to incoherent SR.  */
  if (csbend->isr)
    isrConstant = re_mks*sqrt(55.0/(24*sqrt(3))*pow5(Po)*
                              137.0359895/pow3(fabs(rho_actual)));
  else
    isrConstant = 0;
                            
  Fy_0  = 1;
  Fy_x  = -nh;
  Fy_x2 = Fy_x3 = Fy_x4 = Fy_y2 = Fy_x_y2 = Fy_x2_y2 = Fy_y4 = 0;
  if (csbend->nonlinear) {
    Fy_x2    = betah2;
    Fy_x3    = gammah3;
    Fy_y2    = (h*nh - 2*betah2)/2;
    Fy_x_y2  =  - (2*h*betah2 + nh*h2 + 6*gammah3)/2;
    Fy_x4    = deltah4;
    Fy_x2_y2 =  - (3*h*gammah3 - h3*nh - 2*h2*betah2 + 12*deltah4)/2;
    Fy_y4    = (12*h*gammah3 - h3*nh - 2*h2*betah2 + 24*deltah4)/24;
  }

  Fx_y    =  - nh;
  Fx_x_y  = Fx_x2_y = Fx_x3_y = Fx_y3 = Fx_x_y3 = 0;
  if (csbend->nonlinear) {
    Fx_x_y  = 2*betah2;
    Fx_x2_y = 3*gammah3;
    Fx_y3   =  - (2*h*betah2 + nh*h2 + 6*gammah3)/6;
    Fx_x3_y = 4*deltah4;
    Fx_x_y3 =  - (3*h*gammah3 - h3*nh - 2*h2*betah2 + 12*deltah4)/3;
  }

  ttilt = tilt + etilt;
  if (ttilt==0) {
    cos_ttilt = 1;
    sin_ttilt = 0;
  }
  else if (fabs(fabs(ttilt)-PI)<1e-12) {
    cos_ttilt = -1;
    sin_ttilt = 0;
  }
  else if (fabs(ttilt-PIo2)<1e-12) {
    cos_ttilt = 0;
    sin_ttilt = 1;
  }
  else if (fabs(ttilt+PIo2)<1e-12) {
    cos_ttilt = 0;
    sin_ttilt = -1;
  }
  else {
    cos_ttilt = cos(ttilt);
    sin_ttilt = sin(ttilt);
  }

  if (etilt) {
    /* compute final offsets due to error-tilt of the magnet */
    /* see pages 90-93 of notebook 1 about this */
    double q1a, q2a, q3a;
    double q1b, q2b, q3b;
    double qp1, qp2, qp3; 
    double dz, tan_alpha, k;

    q1a = (1-cos(angle))*rho0*(cos(etilt)-1);
    q2a = 0;
    q3a = (1-cos(angle))*rho0*sin(etilt);
    qp1 = sin(angle)*cos(etilt);
    qp2 = cos(angle);
    k = sqrt(sqr(qp1)+sqr(qp2));
    qp1 /= k;
    qp2 /= k;
    qp3 = sin(angle)*sin(etilt)/k;
    tan_alpha = 1./tan(angle)/cos(etilt);
    q1b = q1a*tan_alpha/(tan(angle)+tan_alpha);
    q2b = -q1b*tan(angle);
    dz  = sqrt(sqr(q1b-q1a)+sqr(q2b-q2a));
    q3b = q3a + qp3*dz;

    dcoord_etilt[0] = sqrt(sqr(q1b) + sqr(q2b));
    dcoord_etilt[1] = tan(atan(tan_alpha)-(PIo2-angle));
    dcoord_etilt[2] = q3b;
    dcoord_etilt[3] = qp3;
    dcoord_etilt[4] = dz*sqrt(1+sqr(qp3));
    dcoord_etilt[5] = 0;
#ifdef DEBUG
    fprintf(stdout, "pre-tilt offsets due to ETILT=%le:  %le %le %le %le %le\n",
            etilt, dcoord_etilt[0], dcoord_etilt[1], dcoord_etilt[2],
            dcoord_etilt[3], dcoord_etilt[4]);
    fflush(stdout);
#endif

    /* rotate by tilt to get into same frame as bend equations. */
    rotate_coordinates(dcoord_etilt, tilt);
#ifdef DEBUG
    fprintf(stdout, "offsets due to ETILT=%le:  %le %le %le %le %le\n",
            etilt, dcoord_etilt[0], dcoord_etilt[1], dcoord_etilt[2],
            dcoord_etilt[3], dcoord_etilt[4]);
    fflush(stdout);
#endif
  }
  else
    fill_double_array(dcoord_etilt, 6L, 0.0);

  dxi = -csbend->dx;
  dzi =  csbend->dz;
  dyi = -csbend->dy;
  
  /* must use the original angle here because the translation is done after
   * the final rotation back
   */
  dxf =  csbend->dx*cos(csbend->angle) + csbend->dz*sin(csbend->angle);
  dzf =  csbend->dx*sin(csbend->angle) - csbend->dz*cos(csbend->angle);
  dyf = csbend->dy;

  i_top = n_part-1;
#if !defined(PARALLEL)
  multipoleKicksDone += n_part*csbend->n_kicks*(csbend->integration_order==4?4:1);
#endif

  for (i_part=0; i_part<=i_top; i_part++) {
    if (!part) {
      fprintf(stdout, "error: null particle array found (working on particle %ld) (track_through_csbend)\n", i_part);
      fflush(stdout);
      abort();
    }
    if (!(coord = part[i_part])) {
      fprintf(stdout, "error: null coordinate pointer for particle %ld (track_through_csbend)\n", i_part);
      fflush(stdout);
      abort();
    }
    if (accepted && !accepted[i_part]) {
      fprintf(stdout, "error: null accepted particle pointer for particle %ld (track_through_csbend)\n", i_part);
      fflush(stdout);
      abort();
    }

    coord[4] += dzi*sqrt(1 + sqr(coord[1]) + sqr(coord[3]));
    coord[0]  = coord[0] + dxi + dzi*coord[1];
    coord[2]  = coord[2] + dyi + dzi*coord[3];

    x  =  coord[0]*cos_ttilt + coord[2]*sin_ttilt;
    y  = -coord[0]*sin_ttilt + coord[2]*cos_ttilt;
    xp =  coord[1]*cos_ttilt + coord[3]*sin_ttilt;
    yp = -coord[1]*sin_ttilt + coord[3]*cos_ttilt;
    s  = coord[4];
    dp = dp0 = coord[5];

    if (csbend->edge1_effects) {
      rho = (1+dp)*rho_actual;
      if (csbend->edge_order<2) {
        /* apply edge focusing */
        delta_xp = tan(e1)/rho*x;
        if (e1_kick_limit>0 && fabs(delta_xp)>e1_kick_limit)
          delta_xp = SIGN(delta_xp)*e1_kick_limit;
        xp += delta_xp;
        yp -= tan(e1-psi1/(1+dp))/rho*y;
      } else
        apply_edge_effects(&x, &xp, &y, &yp, rho, n, e1, he1, psi1*(1+dp), -1);
    }

    /* transform to curvilinear coordinates */
    xp *= (1+x/rho0);
    yp *= (1+x/rho0);

    /* load input coordinates into arrays */
    Qi[0] = x;  Qi[1] = xp;  Qi[2] = y;  Qi[3] = yp;  Qi[4] = 0;  Qi[5] = dp;

    if (csbend->edge1_effects && e1!=0 && rad_coef) {
      /* pre-adjust dp/p to anticipate error made by integrating over entire sector */
      y2 = y*y;
      Fx = (Fx_y + (Fx_x_y + (Fx_x2_y + Fx_x3_y*x)*x)*x + (Fx_y3 + Fx_x_y3*x)*y2)*y;
      Fy = Fy_0 + (Fy_x + (Fy_x2 + (Fy_x3 + Fy_x4*x)*x)*x)*x + (Fy_y2 + (Fy_x_y2 + Fy_x2_y2*x)*x + Fy_y4*y2)*y2;
      dp_prime = -rad_coef*(sqr(Fx)+sqr(Fy))*sqr(1+dp)*sqrt(sqr(1+x/rho0)+sqr(xp)+sqr(yp));
      Qi[5] -= dp_prime*x*tan(e1);
    }

    particle_lost = 0;
    if (!particle_lost) {
      if (csbend->integration_order==4)
        integrate_csbend_ord4(Qf, Qi, csbend->length, csbend->n_kicks, rho0, Po);
      else
        integrate_csbend_ord2(Qf, Qi, csbend->length, csbend->n_kicks, rho0, Po);
    }

    if (particle_lost) {
      if (!part[i_top]) {
        fprintf(stdout, "error: couldn't swap particles %ld and %ld--latter is null pointer (track_through_csbend)\n",
                i_part, i_top);
        fflush(stdout);
        abort();
      }
      SWAP_PTR(part[i_part], part[i_top]);
      if (accepted) {
        if (!accepted[i_top]) {
          fprintf(stdout, 
                  "error: couldn't swap acceptance data for particles %ld and %ld--latter is null pointer (track_through_csbend)\n",
                  i_part, i_top);
          fflush(stdout);
          abort();
        }
        SWAP_PTR(accepted[i_part], accepted[i_top]);
      }
      part[i_top][4] = z_start + s_lost;
      part[i_top][5] = Po*(1+part[i_top][5]);
      i_top--;
      i_part--;
      continue;
    }

    if (csbend->edge2_effects && e2!=0 && rad_coef) {
      /* post-adjust dp/p to correct error made by integrating over entire sector */
      x = Qf[0];
      xp = Qf[1];
      y = Qf[2];
      yp = Qf[3];
      dp = Qf[5];
      y2 = y*y;
      Fx = (Fx_y + (Fx_x_y + (Fx_x2_y + Fx_x3_y*x)*x)*x + (Fx_y3 + Fx_x_y3*x)*y2)*y;
      Fy = Fy_0 + (Fy_x + (Fy_x2 + (Fy_x3 + Fy_x4*x)*x)*x)*x + (Fy_y2 + (Fy_x_y2 + Fy_x2_y2*x)*x + Fy_y4*y2)*y2;
      dp_prime = -rad_coef*(sqr(Fx)+sqr(Fy))*sqr(1+dp)*sqrt(sqr(1+x/rho0)+sqr(xp)+sqr(yp));
      Qf[5] -= dp_prime*x*tan(e2);
    }

    /* get final coordinates */
    if (rad_coef) {
      double p0, p1;
      double beta0, beta1;
      /* fix previous distance information to reflect new velocity--since distance
       * is really time-of-flight at the current velocity 
       */
      p0 = Po*(1+dp0);
      beta0 = p0/sqrt(sqr(p0)+1);
      p1 = Po*(1+Qf[5]);
      beta1 = p1/sqrt(sqr(p1)+1);
      s = beta1*s/beta0 + Qf[4];
    }
    else
      s += Qf[4];
    x = Qf[0];  xp = Qf[1];  y = Qf[2];  yp = Qf[3];  dp = Qf[5];

    /* transform to cartesian coordinates */
    xp /= (1+x/rho0);
    yp /= (1+x/rho0);

    if (csbend->edge2_effects) {
      /* apply edge focusing */
      rho = (1+dp)*rho_actual;
      if (csbend->edge_order<2) {
        delta_xp = tan(e2)/rho*x;
        if (e2_kick_limit>0 && fabs(delta_xp)>e2_kick_limit)
          delta_xp = SIGN(delta_xp)*e2_kick_limit;
        xp += delta_xp;
        yp -= tan(e2-psi2/(1+dp))/rho*y;
      } else
        apply_edge_effects(&x, &xp, &y, &yp, rho, n, e2, he2, psi2*(1+dp), 1);
    }
    
    coord[0] =  x*cos_ttilt -  y*sin_ttilt + dcoord_etilt[0];
    coord[2] =  x*sin_ttilt +  y*cos_ttilt + dcoord_etilt[2];
    coord[1] = xp*cos_ttilt - yp*sin_ttilt + dcoord_etilt[1];
    coord[3] = xp*sin_ttilt + yp*cos_ttilt + dcoord_etilt[3];
    coord[4] = s;
    coord[5] = dp;

    coord[0] += dxf + dzf*coord[1];
    coord[2] += dyf + dzf*coord[3];
    coord[4] += dzf*sqrt(1+ sqr(coord[1]) + sqr(coord[3]));
  }

  return(i_top+1);
}

void integrate_csbend_ord2(double *Qf, double *Qi, double s, long n, double rho0, double p0)
{
  long i;
  double factor, f, phi, ds, dsh, dp, dist;
  double Fx, Fy, x, y, y2;
  double sine, cosi, tang;
  double sin_phi, cos_phi;
  double xp, yp, denom;

#define X0 Qi[0]
#define XP0 Qi[1]
#define Y0 Qi[2]
#define YP0 Qi[3]
#define S0 Qi[4]
#define DPoP0 Qi[5]

#define X Qf[0]
#define QX Qf[1]
#define Y Qf[2]
#define QY Qf[3]
#define S Qf[4]
#define DPoP Qf[5]

  if (!Qf)
    bomb("NULL final coordinates pointer ()", NULL);
  if (!Qi)
    bomb("NULL initial coordinates pointer (integrate_csbend_ord2)", NULL);
  if (n<1)
    bomb("invalid number of steps (integrate_csbend_ord2)", NULL);

  /* calculate canonical momenta (scaled to central momentum) */
  dp = DPoP0;
  f = (1+dp)/sqrt(sqr(1+X0/rho0) + sqr(XP0) + sqr(YP0));
  QX = XP0*f;
  QY = YP0*f;

  X = X0;
  Y = Y0;
  S = S0;
  DPoP = DPoP0;

  ds = s/n;
  dsh = ds/2;
  dist = 0;

  for (i=0; i<n; i++) {
    if (i==0) {
      /* do half-length drift */
      if ((f=sqr(1+DPoP)-sqr(QY))<=0) {
        particle_lost = 1;
        s_lost = dist;
        return;
      }
      f = sqrt(f);
      if (fabs(QX/f)>1) {
        particle_lost = 1;
        s_lost = dist;
        return;
      }
      phi = asin(sin_phi=QX/f);
      sine = sin(dsh/rho0+phi);
      if ((cosi = cos(dsh/rho0+phi))==0) {
        particle_lost = 1;
        s_lost = dist;
        return;
      }
      tang = sine/cosi;
      cos_phi = cos(phi);
      QX = f*sine;
      Y += QY*(factor=(rho0+X)*cos_phi/f*(tang-sin_phi/cos_phi));
      dist += factor*(1+DPoP);
      f = cos_phi/cosi;
      X  = rho0*(f-1) + f*X;
    }

    /* calculate the scaled fields */
    x = X;
    y = Y;
    y2   = y*y;
    Fx = (Fx_y + (Fx_x_y + (Fx_x2_y + Fx_x3_y*x)*x)*x + (Fx_y3 + Fx_x_y3*x)*y2)*y;
    Fy = Fy_0 + (Fy_x + (Fy_x2 + (Fy_x3 + Fy_x4*x)*x)*x)*x + (Fy_y2 + (Fy_x_y2 + Fy_x2_y2*x)*x + Fy_y4*y2)*y2;

    /* do kicks */
    QX += -ds*(1+X/rho0)*Fy/rho_actual;
    QY += ds*(1+X/rho0)*Fx/rho_actual;
    if (rad_coef || isrConstant) {
      denom = sqrt(sqr(1+DPoP)-sqr(QX)-sqr(QY));
      xp = QX/denom;
      yp = QY/denom;
      QX /= (1+DPoP);
      QY /= (1+DPoP);
      DPoP -= rad_coef*sqr(1+DPoP)*(sqr(Fx)+sqr(Fy))*sqrt(sqr(1+X/rho0)+sqr(xp)+sqr(yp))*ds 
        + isrConstant*sqrt(ds)*gauss_rn_lim(0.0, 1.0, 3.0, random_2);
      QX *= (1+DPoP);
      QY *= (1+DPoP);
    }
    
    if (i==n-1) {
      /* do half-length drift */
      if ((f=sqr(1+DPoP)-sqr(QY))<=0) {
        particle_lost = 1;
        s_lost = dist;
        return;
      }
      f = sqrt(f);
      if (fabs(QX/f)>1) {
        particle_lost = 1;
        s_lost = dist;
        return;
      }
      phi = asin(sin_phi=QX/f);
      sine = sin(dsh/rho0+phi);
      if ((cosi = cos(dsh/rho0+phi))==0) {
        particle_lost = 1;
        s_lost = dist;
        return;
      }
      tang = sine/cosi;
      cos_phi = cos(phi);
      QX = f*sine;
      Y += QY*(factor=(rho0+X)*cos_phi/f*(tang-sin_phi/cos_phi));
      dist += factor*(1+DPoP);
      f = cos_phi/cosi;
      X  = rho0*(f-1) + f*X;
    }
    else {
      /* do full-length drift */
      if ((f=sqr(1+DPoP)-sqr(QY))<=0) {
        particle_lost = 1;
        s_lost = dist;
        return;
      }
      f = sqrt(f);
      if (fabs(QX/f)>1) {
        particle_lost = 1;
        s_lost = dist;
        return;
      }
      phi = asin(sin_phi=QX/f);
      sine = sin(ds/rho0+phi);
      if ((cosi = cos(ds/rho0+phi))==0) {
        particle_lost = 1;
        s_lost = dist;
        return;
      }
      tang = sine/cosi;
      cos_phi = cos(phi);
      QX = f*sine;
      Y += QY*(factor=(rho0+X)*cos_phi/f*(tang-sin_phi/cos_phi));
      dist += factor*(1+DPoP);
      f = cos_phi/cosi;
      X  = rho0*(f-1) + f*X;
    }
  }

  /* convert back to slopes */
  f = (1+X/rho0)/sqrt(sqr(1+DPoP)-sqr(QX)-sqr(QY));
  Qf[1] *= f;
  Qf[3] *= f;
  Qf[4] += dist;
}


void integrate_csbend_ord4(double *Qf, double *Qi, double s, long n, double rho0, double p0)
{
  long i;
  double factor, f, phi, ds, dsh, dp, dist;
  double Fx, Fy, x, y, y2;
  double sine, cosi, tang;
  double sin_phi, cos_phi;
  double isrFactor;
  
#define X0 Qi[0]
#define XP0 Qi[1]
#define Y0 Qi[2]
#define YP0 Qi[3]
#define S0 Qi[4]
#define DPoP0 Qi[5]

#define X Qf[0]
#define QX Qf[1]
#define Y Qf[2]
#define QY Qf[3]
#define S Qf[4]
#define DPoP Qf[5]

  /* BETA is 2^(1/3) */
#define BETA 1.25992104989487316477

  if (!Qf)
    bomb("NULL final coordinates pointer ()", NULL);
  if (!Qi)
    bomb("NULL initial coordinates pointer (integrate_csbend_ord4)", NULL);
  if (n<1)
    bomb("invalid number of steps (integrate_csbend_ord4)", NULL);

  /* calculate canonical momenta (scaled to central momentum) */
  dp = DPoP0;
  f = (1+dp)/sqrt(sqr(1+X0/rho0) + sqr(XP0) + sqr(YP0));
  QX = XP0*f;
  QY = YP0*f;

  X = X0;
  Y = Y0;
  S = S0;
  DPoP = DPoP0;

  dist = 0;

  s /= n;
  isrFactor = isrConstant*sqrt(s/3);
  for (i=0; i<n; i++) {
    
    /* do first drift */
    dsh = s/2/(2-BETA);
    if ((f=sqr(1+DPoP)-sqr(QY))<=0) {
      particle_lost = 1;
      s_lost = dist;
      return;
    }
    f = sqrt(f);
    if (fabs(QX/f)>1) {
      particle_lost = 1;
      s_lost = dist;
      return;
    }
    phi = asin(sin_phi=QX/f);
    sine = sin(dsh/rho0+phi);
    if ((cosi = cos(dsh/rho0+phi))==0) {
      particle_lost = 1;
      s_lost = dist;
      return;
    }
    tang = sine/cosi;
    cos_phi = cos(phi);
    QX = f*sine;
    Y += QY*(factor=(rho0+X)*cos_phi/f*(tang-sin_phi/cos_phi));
    dist += factor*(1+DPoP);
    f = cos_phi/cosi;
    X  = rho0*(f-1) + f*X;
    
    /* do first kick */
    ds = s/(2-BETA);
    /* -- calculate the scaled fields */
    x = X;
    y = Y;
    y2   = y*y;
    Fx = (Fx_y + (Fx_x_y + (Fx_x2_y + Fx_x3_y*x)*x)*x + (Fx_y3 + Fx_x_y3*x)*y2)*y;
    Fy = Fy_0 + (Fy_x + (Fy_x2 + (Fy_x3 + Fy_x4*x)*x)*x)*x + (Fy_y2 + (Fy_x_y2 + Fy_x2_y2*x)*x + Fy_y4*y2)*y2;
    /* --do kicks */
    QX += -ds*(1+X/rho0)*Fy/rho_actual;
    QY += ds*(1+X/rho0)*Fx/rho_actual;
    if (rad_coef || isrConstant) {
      QX /= (1+DPoP);
      QY /= (1+DPoP);
      DPoP -= rad_coef*sqr(1+DPoP)*(sqr(Fx)+sqr(Fy))*(1+X/rho0)*ds 
        + isrFactor*gauss_rn_lim(0.0, 1.0, 3.0, random_2);
      QX *= (1+DPoP);
      QY *= (1+DPoP);
    }

    /* do second drift */
    dsh = s*(1-BETA)/(2-BETA)/2;
    if ((f=sqr(1+DPoP)-sqr(QY))<=0) {
      particle_lost = 1;
      s_lost = dist;
      return;
    }
    f = sqrt(f);
    if (fabs(QX/f)>1) {
      particle_lost = 1;
      s_lost = dist;
      return;
    }
    phi = asin(sin_phi=QX/f);
    sine = sin(dsh/rho0+phi);
    if ((cosi = cos(dsh/rho0+phi))==0) {
      particle_lost = 1;
      s_lost = dist;
      return;
    }
    tang = sine/cosi;
    cos_phi = cos(phi);
    QX = f*sine;
    Y += QY*(factor=(rho0+X)*cos_phi/f*(tang-sin_phi/cos_phi));
    dist += factor*(1+DPoP);
    f = cos_phi/cosi;
    X  = rho0*(f-1) + f*X;
    
    /* do second kick */
    ds = -s*BETA/(2-BETA);
    /* -- calculate the scaled fields */
    x = X;
    y = Y;
    y2   = y*y;
    Fx = (Fx_y + (Fx_x_y + (Fx_x2_y + Fx_x3_y*x)*x)*x + (Fx_y3 + Fx_x_y3*x)*y2)*y;
    Fy = Fy_0 + (Fy_x + (Fy_x2 + (Fy_x3 + Fy_x4*x)*x)*x)*x + (Fy_y2 + (Fy_x_y2 + Fy_x2_y2*x)*x + Fy_y4*y2)*y2;
    /* --do kicks */
    QX += -ds*(1+X/rho0)*Fy/rho_actual;
    QY += ds*(1+X/rho0)*Fx/rho_actual;
    if (rad_coef || isrConstant) {
      QX /= (1+DPoP);
      QY /= (1+DPoP);
      DPoP -= rad_coef*sqr(1+DPoP)*(sqr(Fx)+sqr(Fy))*(1+X/rho0)*ds +
        isrFactor*gauss_rn_lim(0.0, 1.0, 3.0, random_2);
      QX *= (1+DPoP);
      QY *= (1+DPoP);
    }
    
    /* do third drift */
    dsh = s*(1-BETA)/(2-BETA)/2;
    if ((f=sqr(1+DPoP)-sqr(QY))<=0) {
      particle_lost = 1;
      s_lost = dist;
      return;
    }
    f = sqrt(f);
    if (fabs(QX/f)>1) {
      particle_lost = 1;
      return;
    }
    phi = asin(sin_phi=QX/f);
    sine = sin(dsh/rho0+phi);
    if ((cosi = cos(dsh/rho0+phi))==0) {
      particle_lost = 1;
      s_lost = dist;
      return;
    }
    tang = sine/cosi;
    cos_phi = cos(phi);
    QX = f*sine;
    Y += QY*(factor=(rho0+X)*cos_phi/f*(tang-sin_phi/cos_phi));
    dist += factor*(1+DPoP);
    f = cos_phi/cosi;
    X  = rho0*(f-1) + f*X;
    
    /* do third kick */
    ds = s/(2-BETA);
    /* -- calculate the scaled fields */
    x = X;
    y = Y;
    y2   = y*y;
    Fx = (Fx_y + (Fx_x_y + (Fx_x2_y + Fx_x3_y*x)*x)*x + (Fx_y3 + Fx_x_y3*x)*y2)*y;
    Fy = Fy_0 + (Fy_x + (Fy_x2 + (Fy_x3 + Fy_x4*x)*x)*x)*x + (Fy_y2 + (Fy_x_y2 + Fy_x2_y2*x)*x + Fy_y4*y2)*y2;
    /* --do kicks */
    QX += -ds*(1+X/rho0)*Fy/rho_actual;
    QY += ds*(1+X/rho0)*Fx/rho_actual;
    if (rad_coef || isrConstant) {
      QX /= (1+DPoP);
      QY /= (1+DPoP);
      DPoP -= rad_coef*sqr(1+DPoP)*(sqr(Fx)+sqr(Fy))*(1+X/rho0)*ds +
        isrFactor*gauss_rn_lim(0.0, 1.0, 3.0, random_2);
      QX *= (1+DPoP);
      QY *= (1+DPoP);
    }
    
    /* do fourth drift */
    dsh = s/2/(2-BETA);
    if ((f=sqr(1+DPoP)-sqr(QY))<=0) {
      particle_lost = 1;
      s_lost = dist;
      return;
    }
    f = sqrt(f);
    if (fabs(QX/f)>1) {
      particle_lost = 1;
      s_lost = dist;
      return;
    }
    phi = asin(sin_phi=QX/f);
    sine = sin(dsh/rho0+phi);
    if ((cosi = cos(dsh/rho0+phi))==0) {
      particle_lost = 1;
      s_lost = dist;
      return;
    }
    tang = sine/cosi;
    cos_phi = cos(phi);
    QX = f*sine;
    Y += QY*(factor=(rho0+X)*cos_phi/f*(tang-sin_phi/cos_phi));
    dist += factor*(1+DPoP);
    f = cos_phi/cosi;
    X  = rho0*(f-1) + f*X;
  }

  /* convert back to slopes */
  f = (1+X/rho0)/sqrt(sqr(1+DPoP)-sqr(QX)-sqr(QY));
  Qf[1] *= f;
  Qf[3] *= f;
  Qf[4] += dist;
}

typedef struct {
  unsigned long lastMode;
#define CSRDRIFT_STUPAKOV          0x0001UL
#define CSRDRIFT_SALDIN54          0x0002UL
#define CSRDRIFT_OVERTAKINGLENGTH  0x0004UL
#define CSRDRIFT_ATTENUATIONLENGTH 0x0008UL
#define CSRDRIFT_SPREAD            0x0010UL
  long bins, valid;
  double dctBin, s0, ds0, zLast, z0;
  double S11, S12, S22;
  double *dGamma;
  double rho, bendingAngle, Po, perc68BunchLength, perc90BunchLength, peakToPeakWavelength, rmsBunchLength;
  /* for Saldin eq 54 (NIM A 398 (1997) 373-394) mode: */
  FILE *fpSaldin;
  long nSaldin;
  double *FdNorm;   /* Saldin's Fd(sh, x)/Fd(sh, 0), sh = bunch-length*gamma^3/rho */
  double *xSaldin;  /* distance from end of bend */
  double lastFdNorm; /* last value obtained from interpolation */
  /* for Stupakov mode */
  long SGOrder, SGHalfWidth, SGDerivHalfWidth, SGDerivOrder;
  double binRangeFactor;
  double GSConstant, MPCharge;
  char *StupakovOutput;
  SDDS_DATASET SDDS_Stupakov;
  long StupakovFileActive, StupakovOutputInterval;
  long trapazoidIntegration;
  double highFrequencyCutoff0, highFrequencyCutoff1;
} CSR_LAST_WAKE;
CSR_LAST_WAKE csrWake;

#define DERBENEV_CRITERION_DISABLE 0
#define DERBENEV_CRITERION_EVAL 1
#define DERBENEV_CRITERION_ENFORCE 2
#define N_DERBENEV_CRITERION_OPTIONS 3
static char *derbenevCriterionOption[N_DERBENEV_CRITERION_OPTIONS] = {
  "disable", "evaluate", "enforce"};

long track_through_csbendCSR(double **part, long n_part, CSRCSBEND *csbend, double p_error, 
                             double Po, double **accepted, double z_start, double z_end,
                             CHARGE *charge, char *rootname)
{
  double nh, betah2, gammah3, deltah4;
  double h, h2, h3, he1, he2;
  static long csrWarning = 0;
  static double *beta0=NULL, *ctHist=NULL, *ctHistDeriv=NULL;
  static double *dGamma=NULL, *T1=NULL, *T2=NULL, *denom=NULL;
  static long maxParticles = 0, maxBins = 0 ;
  static char *particleLost=NULL;
  double x, xp, y, yp, p1, beta1, p0;
  double ctLower, ctUpper, dct, slippageLength, phiBend, slippageLength13;
  long diSlippage, diSlippage4;
  long nBins, nBinned;
  long i_part, i_top, kick;
  double rho=0.0, Fx, Fy, y2;
  double n, beta, gamma, delta, fse, dp_prime;
  double tilt, etilt, cos_ttilt, sin_ttilt, ttilt;
  double *coord;
  double angle, e1, e2, Kg;
  double psi1, psi2;
  double Qi[6], Qf[6];
  double dcoord_etilt[6];
  double dxi, dyi, dzi;
  double dxf, dyf, dzf;
  double delta_xp;
  double macroParticleCharge, CSRConstant;
  long iBin, iBinBehind;
  long csrInhibit = 0;
  double derbenevRatio = 0;
  TRACKING_CONTEXT tContext;

  reset_driftCSR();
  getTrackingContext(&tContext);

  if (!csbend)
    bomb("null CSBEND pointer (track_through_csbend)", NULL);

  if (csbend->angle==0) {
    exactDrift(part, n_part, csbend->length);
    return n_part;
  }

  if (csbend->integration_order!=2 && csbend->integration_order!=4)
    bomb("CSBEND integration_order is invalid--must be either 2 or 4", NULL);

  macroParticleCharge = 0;
  if (charge) {
    macroParticleCharge = charge->macroParticleCharge;
  } else if (csbend->bins && !csrWarning) {
    fprintf(stdout, "Warning: you asked for CSR on CSBEND but didn't give a CHARGE element\n");
    fflush(stdout);
    csrWarning = 1;
  }
  
  if ((nBins=csbend->bins)<2)
    bomb("Less than 2 bins for CSR!", NULL);

  if (csbend->SGDerivHalfWidth<=0)
    csbend->SGDerivHalfWidth = csbend->SGHalfWidth;
  if (csbend->SGDerivHalfWidth<=0)
    csbend->SGDerivHalfWidth = 1;

  if (csbend->SGDerivOrder<=0)
    csbend->SGDerivOrder = csbend->SGOrder;
  if (csbend->SGDerivOrder<=0)
    csbend->SGDerivOrder = 1;
  
  if (n_part>maxParticles &&
      (!(beta0=SDDS_Realloc(beta0, sizeof(*beta0)*(maxParticles=n_part))) ||
       !(particleLost=SDDS_Realloc(particleLost, sizeof(*particleLost)*n_part))))
    bomb("Memory allocation failure (track_through_csbendCSR)", NULL);

  if (csbend->use_bn) {
    rho0 = csbend->length/csbend->angle;
    csbend->k1_internal = csbend->b1/rho0;
    csbend->k2_internal = csbend->b2/rho0;
    csbend->k3_internal = csbend->b3/rho0;
    csbend->k4_internal = csbend->b4/rho0;
  } else {
    csbend->k1_internal = csbend->k1;
    csbend->k2_internal = csbend->k2;
    csbend->k3_internal = csbend->k3;
    csbend->k4_internal = csbend->k4;
  }

  he1 = csbend->h1;
  he2 = csbend->h2;
  if (csbend->angle<0) {
    angle = -csbend->angle;
    e1    = -csbend->e1;
    e2    = -csbend->e2;
    etilt = csbend->etilt;
    tilt  = csbend->tilt + PI;
    rho0  = -csbend->length/angle;
    n     = -sqr(rho0)*csbend->k1_internal;
    beta  = 0.5*csbend->k2_internal*pow3(rho0);
    gamma = csbend->k3_internal*pow4(rho0)/6.;
    delta = csbend->k4_internal*pow5(rho0)/24.;
    rho0  = -rho0;
  }
  else {
    angle = csbend->angle;
    e1    = csbend->e1;
    e2    = csbend->e2;
    etilt = csbend->etilt;
    tilt  = csbend->tilt;
    rho0  = csbend->length/angle;
    n     = -sqr(rho0)*csbend->k1_internal;
    beta  = 0.5*csbend->k2_internal*pow3(rho0);
    gamma = csbend->k3_internal*pow4(rho0)/6.;
    delta = csbend->k4_internal*pow5(rho0)/24.;
  }

  fse = csbend->fse;
  h2 = sqr(h=1./rho0);
  h3 = h*h2;
  nh = n*h;
  betah2 = beta*h2;
  gammah3 = gamma*h3;
  deltah4 = delta*h2*h2;
  if (fse>-1)
    rho_actual = 1/((1+fse)*h);
  else
    rho_actual = 1e16/h;

  /* angles for fringe-field effects */
  Kg   = 2*csbend->hgap*csbend->fint;
  psi1 = Kg/rho_actual/cos(e1)*(1+sqr(sin(e1)));
  psi2 = Kg/rho_actual/cos(e2)*(1+sqr(sin(e2)));

  /* rad_coef is d((P-Po)/Po)/ds for the on-axis, on-momentum particle, where po is the momentum of
   * the central particle.
   */
  if (csbend->synch_rad)
    rad_coef = sqr(e_mks)*pow3(Po)*sqr(1+fse)/(6*PI*epsilon_o*sqr(c_mks)*me_mks*sqr(rho0));
  else
    rad_coef = 0;
  /* isrConstant is the RMS increase in dP/P per meter due to incoherent SR.  */
  if (csbend->isr) 
    isrConstant = re_mks*sqrt(55.0/(24*sqrt(3))*pow5(Po)*
                              137.0359895/pow3(fabs(rho_actual)));
  else
    isrConstant = 0;
  Fy_0  = 1;
  Fy_x  = -nh;
  Fy_x2 = Fy_x3 = Fy_x4 = Fy_y2 = Fy_x_y2 = Fy_x2_y2 = Fy_y4 = 0;
  if (csbend->nonlinear) {
    Fy_x2    = betah2;
    Fy_x3    = gammah3;
    Fy_y2    = (h*nh - 2*betah2)/2;
    Fy_x_y2  =  - (2*h*betah2 + nh*h2 + 6*gammah3)/2;
    Fy_x4    = deltah4;
    Fy_x2_y2 =  - (3*h*gammah3 - h3*nh - 2*h2*betah2 + 12*deltah4)/2;
    Fy_y4    = (12*h*gammah3 - h3*nh - 2*h2*betah2 + 24*deltah4)/24;
  }

  Fx_y    =  - nh;
  Fx_x_y  = Fx_x2_y = Fx_x3_y = Fx_y3 = Fx_x_y3 = 0;
  if (csbend->nonlinear) {
    Fx_x_y  = 2*betah2;
    Fx_x2_y = 3*gammah3;
    Fx_y3   =  - (2*h*betah2 + nh*h2 + 6*gammah3)/6;
    Fx_x3_y = 4*deltah4;
    Fx_x_y3 =  - (3*h*gammah3 - h3*nh - 2*h2*betah2 + 12*deltah4)/3;
  }

  ttilt = tilt + etilt;
  if (ttilt==0) {
    cos_ttilt = 1;
    sin_ttilt = 0;
  }
  else if (fabs(fabs(ttilt)-PI)<1e-12) {
    cos_ttilt = -1;
    sin_ttilt = 0;
  }
  else if (fabs(ttilt-PIo2)<1e-12) {
    cos_ttilt = 0;
    sin_ttilt = 1;
  }
  else if (fabs(ttilt+PIo2)<1e-12) {
    cos_ttilt = 0;
    sin_ttilt = -1;
  }
  else {
    cos_ttilt = cos(ttilt);
    sin_ttilt = sin(ttilt);
  }


  if (etilt) {
    /* compute final offsets due to error-tilt of the magnet */
    /* see pages 90-93 of notebook 1 about this */
    double q1a, q2a, q3a;
    double q1b, q2b, q3b;
    double qp1, qp2, qp3; 
    double dz, tan_alpha, k;

    q1a = (1-cos(angle))*rho0*(cos(etilt)-1);
    q2a = 0;
    q3a = (1-cos(angle))*rho0*sin(etilt);
    qp1 = sin(angle)*cos(etilt);
    qp2 = cos(angle);
    k = sqrt(sqr(qp1)+sqr(qp2));
    qp1 /= k;
    qp2 /= k;
    qp3 = sin(angle)*sin(etilt)/k;
    tan_alpha = 1./tan(angle)/cos(etilt);
    q1b = q1a*tan_alpha/(tan(angle)+tan_alpha);
    q2b = -q1b*tan(angle);
    dz  = sqrt(sqr(q1b-q1a)+sqr(q2b-q2a));
    q3b = q3a + qp3*dz;

    dcoord_etilt[0] = sqrt(sqr(q1b) + sqr(q2b));
    dcoord_etilt[1] = tan(atan(tan_alpha)-(PIo2-angle));
    dcoord_etilt[2] = q3b;
    dcoord_etilt[3] = qp3;
    dcoord_etilt[4] = dz*sqrt(1+sqr(qp3));
    dcoord_etilt[5] = 0;

    /* rotate by tilt to get into same frame as bend equations. */
    rotate_coordinates(dcoord_etilt, tilt);
  }
  else
    fill_double_array(dcoord_etilt, 6L, 0.0);

  dxi = -csbend->dx;
  dzi =  csbend->dz;
  dyi = -csbend->dy;

  /* must use the original angle here because the translation is done after
   * the final rotation back
   */
  dxf =  csbend->dx*cos(csbend->angle) + csbend->dz*sin(csbend->angle);
  dzf =  csbend->dx*sin(csbend->angle) - csbend->dz*cos(csbend->angle);
  dyf = csbend->dy;

  if (csbend->particleOutputFile && !csbend->particleFileActive) {
    /* set up SDDS output file for particle coordinates inside bend */
    csbend->particleFileActive = 1;
    csbend->particleOutputFile = compose_filename(csbend->particleOutputFile, rootname);
    if (!SDDS_InitializeOutput(&csbend->SDDSpart, SDDS_BINARY, 1, 
                               NULL, NULL, csbend->particleOutputFile) ||
        !SDDS_DefineSimpleParameter(&csbend->SDDSpart, "Pass", NULL, SDDS_LONG) ||
        !SDDS_DefineSimpleParameter(&csbend->SDDSpart, "Kick", NULL, SDDS_LONG) ||
        !SDDS_DefineSimpleParameter(&csbend->SDDSpart, "pCentral", "m$be$nc", SDDS_DOUBLE) ||
        !SDDS_DefineSimpleParameter(&csbend->SDDSpart, "Angle", NULL, SDDS_DOUBLE) ||
        (csbend->xIndex=SDDS_DefineColumn(&csbend->SDDSpart, "x", NULL, "m", 
                                                NULL, NULL, SDDS_DOUBLE, 0 ))<0 ||
        (csbend->xpIndex=SDDS_DefineColumn(&csbend->SDDSpart, "xp", NULL, NULL, 
                                                 NULL, NULL, SDDS_DOUBLE, 0))<0 ||
        (csbend->tIndex=SDDS_DefineColumn(&csbend->SDDSpart, "t", NULL, "s", 
                                                NULL, NULL, SDDS_DOUBLE, 0))<0 ||
        (csbend->pIndex=SDDS_DefineColumn(&csbend->SDDSpart, "p", NULL, "m$be$nc", 
                                                NULL, NULL, SDDS_DOUBLE, 0))<0 ||
        !SDDS_WriteLayout(&csbend->SDDSpart)) {
      SDDS_SetError("Problem setting up particle output file for CSR");
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
    }
  }
  
  if (csbend->histogramFile && !csbend->wakeFileActive) {
    /* set up SDDS output file for CSR monitoring */
    csbend->wakeFileActive = 1;
    csbend->histogramFile = compose_filename(csbend->histogramFile, rootname);
    if (!SDDS_InitializeOutput(&csbend->SDDSout, SDDS_BINARY, 1, NULL, NULL, csbend->histogramFile) ||
        !SDDS_DefineSimpleParameter(&csbend->SDDSout, "Pass", NULL, SDDS_LONG) ||
        !SDDS_DefineSimpleParameter(&csbend->SDDSout, "Kick", NULL, SDDS_LONG) ||
        !SDDS_DefineSimpleParameter(&csbend->SDDSout, "pCentral", "m$be$nc", SDDS_DOUBLE) ||
        !SDDS_DefineSimpleParameter(&csbend->SDDSout, "Angle", NULL, SDDS_DOUBLE) ||
        !SDDS_DefineSimpleParameter(&csbend->SDDSout, "SlippageLength", "m", SDDS_DOUBLE) ||
        !SDDS_DefineSimpleParameter(&csbend->SDDSout, "TotalBunchLength", "m", SDDS_DOUBLE) ||
        !SDDS_DefineSimpleParameter(&csbend->SDDSout, "BinSize", "m", SDDS_DOUBLE) ||
        !SDDS_DefineSimpleParameter(&csbend->SDDSout, "dsKick", "m", SDDS_DOUBLE) ||
        !SDDS_DefineSimpleParameter(&csbend->SDDSout, "DerbenevRatio", NULL, SDDS_DOUBLE) ||
        !SDDS_DefineSimpleColumn(&csbend->SDDSout, "s", "m", SDDS_DOUBLE) ||
        !SDDS_DefineSimpleColumn(&csbend->SDDSout, "LinearDensity", "C/s", SDDS_DOUBLE) ||
        !SDDS_DefineSimpleColumn(&csbend->SDDSout, "LinearDensityDeriv", "C/s$a2$n", SDDS_DOUBLE) ||
        !SDDS_DefineSimpleColumn(&csbend->SDDSout, "DeltaGamma", NULL, SDDS_DOUBLE) ||
        !SDDS_DefineSimpleColumn(&csbend->SDDSout, "GammaDeriv", "1/m", SDDS_DOUBLE) ||
        !SDDS_DefineSimpleColumn(&csbend->SDDSout, "DeltaGammaT1", NULL, SDDS_DOUBLE) ||
        !SDDS_DefineSimpleColumn(&csbend->SDDSout, "DeltaGammaT2", NULL, SDDS_DOUBLE) ||
        !SDDS_WriteLayout(&csbend->SDDSout)) {
      SDDS_SetError("Problem setting up wake output file for CSR");
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
    }
  }

  /*  prepare arrays for CSR integrals */
  nBins = csbend->bins;
  if (!(ctHist=SDDS_Realloc(ctHist, sizeof(*ctHist)*nBins)) ||
      !(ctHistDeriv=SDDS_Realloc(ctHistDeriv, sizeof(*ctHistDeriv)*nBins)) ||
      !(denom=SDDS_Realloc(denom, sizeof(*denom)*nBins)) ||
      !(T1=SDDS_Realloc(T1, sizeof(*T1)*nBins)) ||
      !(T2=SDDS_Realloc(T2, sizeof(*T2)*nBins)) ||
      !(dGamma=SDDS_Realloc(dGamma, sizeof(*dGamma)*nBins)))
      bomb("memory allocation failure (track_through_csbendCSR)", NULL);

  csrWake.dGamma = dGamma;
  csrWake.bins = nBins;
  csrWake.ds0 = csbend->length/csbend->n_kicks;
  csrWake.zLast = csrWake.z0 = z_end;
  csrWake.highFrequencyCutoff0 = csbend->highFrequencyCutoff0;
  csrWake.highFrequencyCutoff1 = csbend->highFrequencyCutoff1;
  
#if !defined(PARALLEL)  
  multipoleKicksDone += n_part*csbend->n_kicks*(csbend->integration_order==4?4:1);
#endif

  /* check particle data, transform coordinates, and handle edge effects */
  for (i_part=0; i_part<n_part; i_part++) {
    if (!part) {
      fprintf(stdout, "error: null particle array found (working on particle %ld) (track_through_csbend)\n", i_part);
      fflush(stdout);
      abort();
    }
    if (!(coord = part[i_part])) {
      fprintf(stdout, "error: null coordinate pointer for particle %ld (track_through_csbend)\n", i_part);
      fflush(stdout);
      abort();
    }
    if (accepted && !accepted[i_part]) {
      fprintf(stdout, "error: null accepted particle pointer for particle %ld (track_through_csbend)\n", i_part);
      fflush(stdout);
      abort();
    }

    /* adjust for element offsets */
    coord[4] += dzi*sqrt(1 + sqr(coord[1]) + sqr(coord[3]));
    coord[0]  = coord[0] + dxi + dzi*coord[1];
    coord[2]  = coord[2] + dyi + dzi*coord[3];

    /* perform tilt transformations and save some data */
    x  =  coord[0]*cos_ttilt + coord[2]*sin_ttilt;
    y  = -coord[0]*sin_ttilt + coord[2]*cos_ttilt;
    coord[0] = x;
    coord[2] = y;
    xp =  coord[1]*cos_ttilt + coord[3]*sin_ttilt;
    yp = -coord[1]*sin_ttilt + coord[3]*cos_ttilt;
    coord[1] = xp;
    coord[3] = yp;
    p0 = Po*(1+coord[5]);
    beta0[i_part] = p0/sqrt(p0*p0+1);
    coord[4] /= beta0[i_part];
    particleLost[i_part] = 0;

#undef X
#undef Y
#define X coord[0]
#define Y coord[2]
#define XP coord[1]
#define YP coord[3]
#define CT coord[4]
#define DP coord[5]
    if (csbend->edge1_effects) {
      /* apply edge focusing */
      rho = (1+DP)*rho_actual;
      if (csbend->edge_order<2) {
        delta_xp = tan(e1)/rho*X;
        XP += delta_xp;
        YP -= tan(e1-psi1/(1+DP))/rho*Y;
      }
      else
        apply_edge_effects(&X, &XP, &Y, &YP, rho, n, e1, he1, psi1*(1+DP), -1);
    }

    /* transform to curvilinear coordinates */
    XP *= (1+X/rho0);
    YP *= (1+X/rho0);

    if (csbend->edge1_effects && e1!=0 && rad_coef) {
      /* pre-adjust dp/p to anticipate error made by integrating over entire sector */
      y2 = Y*Y;
      Fx = (Fx_y + (Fx_x_y + (Fx_x2_y + Fx_x3_y*X)*X)*X
            + (Fx_y3 + Fx_x_y3*X)*y2)*Y;
      Fy = Fy_0 + (Fy_x + (Fy_x2 + (Fy_x3 + Fy_x4*X)*X)*X)*X
        + (Fy_y2 + (Fy_x_y2 + Fy_x2_y2*X)*X + Fy_y4*y2)*y2;
      dp_prime = -rad_coef*(sqr(Fx)+sqr(Fy))*sqr(1+DP)*
        sqrt(sqr(1+X/rho0)+sqr(XP)+sqr(YP));
      DP -= dp_prime*X*tan(e1);
    }
  }

  if (csbend->csr)
    CSRConstant = 2*macroParticleCharge*e_mks/pow(3*rho0*rho0, 1./3.)/(4*PI*epsilon_o*me_mks*sqr(c_mks));
  else
    CSRConstant = 0;
  /* Now do the body of the sector dipole */
  for (kick=phiBend=0; kick<csbend->n_kicks; kick++) {
    for (i_part=0; i_part<n_part; i_part++) {
      coord = part[i_part];
      if (particleLost[i_part])
        continue;
      
      /* load input coordinates into arrays */
      Qi[0] = X;
      Qi[1] = XP;
      Qi[2] = Y;
      Qi[3] = YP;
      Qi[4] = 0;  
      Qi[5] = DP;

      if (csbend->integration_order==4)
        integrate_csbend_ord4(Qf, Qi, csbend->length/csbend->n_kicks, 1, rho0, Po);
      else
        integrate_csbend_ord2(Qf, Qi, csbend->length/csbend->n_kicks, 1, rho0, Po);
      particleLost[i_part] = particle_lost;
      
      /* retrieve coordinates from arrays */
      X  = Qf[0];  
      XP = Qf[1];  
      Y  = Qf[2];  
      YP = Qf[3];  
      DP = Qf[5];
      if (rad_coef) {
        /* convert additional distance traveled to ct using mean velocity */
        p1 = Po*(1+DP);
        beta1 = p1/sqrt(p1*p1+1);
        CT += Qf[4]*2/(beta0[i_part]+beta1);
        beta0[i_part] = beta1;
      } else {
        CT += Qf[4]/beta0[i_part];  
      }
    }

    if (n_part>1 && csbend->derbenevCriterionMode) {
      /* evaluate Derbenev criterion: sigma_x/sigma_z << (sigma_z/R)^(1/3) */
      long code;
      double Sz, Sx;
      switch (code=match_string(csbend->derbenevCriterionMode, derbenevCriterionOption, N_DERBENEV_CRITERION_OPTIONS, 0)) {
      case DERBENEV_CRITERION_DISABLE:
        break;
      case DERBENEV_CRITERION_EVAL:
      case DERBENEV_CRITERION_ENFORCE:
        rms_emittance(part, 4, 5, n_part, &Sz, NULL, NULL);
        Sz = sqrt(Sz);
        rms_emittance(part, 0, 1, n_part, &Sx, NULL, NULL);
        Sx = sqrt(Sx);
        derbenevRatio = (Sx/Sz)/pow(rho/Sz, 1./3.);
        if (derbenevRatio>0.1) {
          if (code==DERBENEV_CRITERION_EVAL)
            fprintf(stderr, "Warning: Using 1-D CSR formalism but Derbenev criterion not satisfied (%le > 0.1).\n",
                    derbenevRatio);
          else {
            csrInhibit = 1;
            fprintf(stderr, "Warning: Derbenev criterion not satisfied (%le > 0.1)---not applying CSR\n",
                    derbenevRatio);
          }
        }
        break;
      default:
        fprintf(stderr, "Error: invalid value for DERBENEV_CRITERION_MODE. Give 'disable', 'evaluate', or 'enforce'\n");
        exit(1);
        break;
      }
    }
    if (n_part>1 && !csrInhibit) {
      /* compute CSR potential function */
      if (kick==0 || !csbend->binOnce) {
        /* - first make a density histogram */
        ctLower = ctUpper = dct = 0;
        if ((nBinned = 
             binParticleCoordinate(&ctHist, &maxBins,
                                   &ctLower, &ctUpper, &dct, &nBins, 
                                   csbend->binRangeFactor<1.1?1.1:csbend->binRangeFactor, 
                                   part, n_part, 4))!=n_part) {
          fprintf(stdout, "Only %ld of %ld particles binned for CSR\n", nBinned, n_part);
          fflush(stdout);
        }
        
        /* - smooth the histogram, normalize to get linear density, and 
           copy in preparation for taking derivative
           */
        if (csbend->highFrequencyCutoff0>0) {
          long nz;
          nz = applyLowPassFilter(ctHist, nBins, dct, csbend->highFrequencyCutoff0, csbend->highFrequencyCutoff1);
          if (nz) {
            fprintf(stdout, "Warning: high pass filter resulted in negative values in %ld bins\n",
                    nz);
            fflush(stdout);
          }
        }
        if (csbend->SGHalfWidth>0) {
          SavitzyGolaySmooth(ctHist, nBins, csbend->SGOrder, csbend->SGHalfWidth, csbend->SGHalfWidth,  0);
          correctDistribution(ctHist, nBins, 1.0*nBinned);
        }
        for (iBin=0; iBin<nBins; iBin++) {
          denom[iBin] = pow(dct*iBin, 1./3.);
          ctHistDeriv[iBin] = (ctHist[iBin] /= dct);
        }
        /* - compute derivative with smoothing.  The deriv is w.r.t. index number and
         * I won't scale it now as it will just fall out in the integral 
         */
        SavitzyGolaySmooth(ctHistDeriv, nBins, csbend->SGDerivOrder, 
                           csbend->SGDerivHalfWidth, csbend->SGDerivHalfWidth, 1);
      } else {
        ctLower += rho0*angle/csbend->n_kicks;
        ctUpper += rho0*angle/csbend->n_kicks;
      }
      
      
      phiBend += angle/csbend->n_kicks;
      slippageLength = rho0*ipow(phiBend, 3)/24.0;
      slippageLength13 = pow(slippageLength, 1./3.);
      diSlippage = slippageLength/dct;
      diSlippage4 = 4*slippageLength/dct;
      for (iBin=0; iBin<nBins; iBin++) {
        double term1, term2;
        long count;
        T1[iBin] = T2[iBin] = 0;
        term1 = term2 = 0;
        if (CSRConstant) {
          if (csbend->steadyState) {
            if (!csbend->trapazoidIntegration) {
              for (iBinBehind=iBin+1; iBinBehind<nBins; iBinBehind++)
                T1[iBin] += ctHistDeriv[iBinBehind]/denom[iBinBehind-iBin];
            }
            else {
              if ((iBinBehind=iBin+1)<nBins)
                term1 = ctHistDeriv[iBinBehind]/denom[iBinBehind-iBin];
              for (count=0, iBinBehind=iBin+1; iBinBehind<nBins; iBinBehind++, count++)
                T1[iBin] += (term2=ctHistDeriv[iBinBehind]/denom[iBinBehind-iBin]);
              if ((iBin+1)<nBins)
                T1[iBin] += 0.3*sqr(denom[1])*(2*ctHistDeriv[iBin+1]+3*ctHistDeriv[iBin])/dct;
              if (count>1)
                T1[iBin] -= (term1+term2)/2;
            }
          } else {
            if (!csbend->trapazoidIntegration) {
              for (iBinBehind=iBin+1; iBinBehind<=(iBin+diSlippage) && iBinBehind<nBins; iBinBehind++)
                T1[iBin] += ctHistDeriv[iBinBehind]/denom[iBinBehind-iBin];
            }
            else {
              if ((iBinBehind = iBin+1)<nBins && iBinBehind<=(iBin+diSlippage))
                term1 = ctHistDeriv[iBinBehind]/denom[iBinBehind-iBin]/2;
              for (count=0, iBinBehind = iBin+1; iBinBehind<=(iBin+diSlippage) && iBinBehind<nBins; 
                   count++, iBinBehind++)
                T1[iBin] += (term2=ctHistDeriv[iBinBehind]/denom[iBinBehind-iBin]);
              if (diSlippage>0 && (iBin+1)<nBins)
                T1[iBin] += 0.3*sqr(denom[1])*(2*ctHistDeriv[iBin+1]+3*ctHistDeriv[iBin])/dct;
              if (count>1)
                T1[iBin] -= (term1+term2)/2;
            }
            if ((iBin+diSlippage)<nBins)
              T2[iBin] += ctHist[iBin+diSlippage];
            if ((iBin+diSlippage4)<nBins)
              T2[iBin] -= ctHist[iBin+diSlippage4];
          }
          /* there is no negative sign here because my derivative is w.r.t. -s
             in notation of Saldin, et. al. */
          T1[iBin] *= CSRConstant*csbend->length/csbend->n_kicks; 
          /* keep the negative sign on this term, which has no derivative */
          T2[iBin] *= -CSRConstant*csbend->length/csbend->n_kicks/slippageLength13;
        }
        dGamma[iBin] = T1[iBin]+T2[iBin];
      }

      if (CSRConstant) {
        for (i_part=0; i_part<n_part; i_part++) {
          long nBins1;
          nBins1 = nBins-1;
          coord = part[i_part];
          if (!particleLost[i_part]) {
            double f;
            /* apply CSR kick */
            iBin = (f=(CT-ctLower)/dct);
            f -= iBin;
            if (iBin>=0 && iBin<nBins1)
              DP += ((1-f)*dGamma[iBin]+f*dGamma[iBin+1])/Po;
          }
        }
      }
      
      if (csbend->particleFileActive && kick%csbend->particleOutputInterval==0) {
        long ip;
        /* dump particle data at this location */
        if (!SDDS_StartPage(&csbend->SDDSpart, n_part) ||
            !SDDS_SetParameters(&csbend->SDDSpart, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, 
                                "Pass", -1, "Kick", kick, "pCentral", Po, "Angle", phiBend, 
                                NULL))
          SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
        for (ip=0; ip<n_part; ip++) {
          if (!SDDS_SetRowValues(&csbend->SDDSpart, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE,
                                 ip, 
                                 csbend->xIndex, part[ip][0],
                                 csbend->xpIndex, part[ip][1],
                                 csbend->tIndex, part[ip][4],
                                 csbend->pIndex, Po*(1+part[ip][5]),
                                 -1)) 
            SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
        }
        if (!SDDS_WritePage(&csbend->SDDSpart))
          SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
      }

      if (tContext.sliceAnalysis && tContext.sliceAnalysis->active &&
	  kick!=(csbend->n_kicks-1) &&
	  (csbend->sliceAnalysisInterval==0 ||
	   kick%csbend->sliceAnalysisInterval==0)) {
	convertFromCSBendCoords(part, n_part, rho0, cos_ttilt, sin_ttilt, 1);
	performSliceAnalysisOutput(tContext.sliceAnalysis, part, n_part, 
				   0, tContext.step, Po, 
				   macroParticleCharge*n_part,
				   tContext.elementName, 
				   z_start + (kick*(z_end-z_start))/(csbend->n_kicks-1),
				   1);
	convertToCSBendCoords(part, n_part, rho0, cos_ttilt, sin_ttilt, 1);
      }

      if (csbend->wakeFileActive && 
          ((!csbend->outputLastWakeOnly && kick%csbend->outputInterval==0) ||
           (csbend->outputLastWakeOnly && kick==(csbend->n_kicks-1)))) {
        /* scale the linear density and its derivative to get C/s and C/s^2 
         * ctHist is already normalized to dct, but ctHistDeriv requires an additional factor
         */
        for (iBin=0; iBin<nBins; iBin++) {
          ctHist[iBin] *= macroParticleCharge*c_mks;
          ctHistDeriv[iBin] *= macroParticleCharge*sqr(c_mks)/dct;
        }
        if (!SDDS_StartPage(&csbend->SDDSout, nBins) ||
            !SDDS_SetColumn(&csbend->SDDSout, SDDS_SET_BY_NAME, dGamma, nBins, "DeltaGamma") ||
            !SDDS_SetColumn(&csbend->SDDSout, SDDS_SET_BY_NAME, T1, nBins, "DeltaGammaT1") ||
            !SDDS_SetColumn(&csbend->SDDSout, SDDS_SET_BY_NAME, T2, nBins, "DeltaGammaT2") ||
            !SDDS_SetColumn(&csbend->SDDSout, SDDS_SET_BY_NAME, ctHist, nBins, "LinearDensity") ||
            !SDDS_SetColumn(&csbend->SDDSout, SDDS_SET_BY_NAME, ctHistDeriv, nBins, "LinearDensityDeriv") ||
            !SDDS_SetParameters(&csbend->SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, 
                                "Pass", -1, "Kick", kick, "dsKick", csbend->length/csbend->n_kicks,
                                "pCentral", Po, "Angle", phiBend, "SlippageLength", slippageLength,
                                "TotalBunchLength", ctUpper-ctLower,
                                "BinSize", dct, 
                                "DerbenevRatio", derbenevRatio, NULL))
          SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
        if (csbend->binOnce) {
          /* fix these arrays so they can be used again */
          ctHist[iBin] /= macroParticleCharge*c_mks;
          ctHistDeriv[iBin] /= macroParticleCharge*sqr(c_mks)/dct;
        }
        /* use T1 array to output s and T2 to output dGamma/ds */
        for (iBin=0; iBin<nBins; iBin++) {
          T1[iBin] = ctLower-(ctLower+ctUpper)/2.0+dct*(iBin+0.5);
          T2[iBin] = dGamma[iBin]/(csbend->length/csbend->n_kicks);
        }
        if (!SDDS_SetColumn(&csbend->SDDSout, SDDS_SET_BY_NAME, T1, nBins, "s") ||
            !SDDS_SetColumn(&csbend->SDDSout, SDDS_SET_BY_NAME, T2, nBins, "GammaDeriv") ||
            !SDDS_WritePage(&csbend->SDDSout))
          SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
      }
    }
  }

  if (!csbend->binOnce && n_part>1 && !csrInhibit && !csbend->csrBlock) {
    /* prepare some data for use by CSRDRIFT element */
    csrWake.dctBin = dct;
    ctLower = ctUpper = dct = 0;
    if ((nBinned = 
         binParticleCoordinate(&ctHist, &maxBins,
                               &ctLower, &ctUpper, &dct, &nBins, 
                               csbend->binRangeFactor<1.1?1.1:csbend->binRangeFactor, 
                               part, n_part, 4))!=n_part) {
      fprintf(stdout, "Only %ld of %ld particles binned for CSR\n", nBinned, n_part);
      fflush(stdout);
    }
    csrWake.s0 = ctLower + dzf;
  } else {
    ctLower = ctUpper = dct = 0;
    csrWake.dctBin = dct;
    csrWake.s0 = ctLower + dzf;
  }
  
  /* remove lost particles, handle edge effects, and transform coordinates */    
  i_top = n_part-1;
  for (i_part=0; i_part<=i_top; i_part++) {
    coord = part[i_part];
    if (csbend->edge2_effects && e2!=0 && rad_coef) {
      /* post-adjust dp/p to correct error made by integrating over entire sector */
      y2 = Y*Y;
      Fx = (Fx_y + (Fx_x_y + (Fx_x2_y + Fx_x3_y*X)*X)*X
            + (Fx_y3 + Fx_x_y3*X)*y2)*Y;
      Fy = Fy_0 + 
        (Fy_x + (Fy_x2 + (Fy_x3 + Fy_x4*X)*X)*X)*X 
          + (Fy_y2 + (Fy_x_y2 + Fy_x2_y2*X)*X + Fy_y4*y2)*y2;
      dp_prime = -rad_coef*(sqr(Fx)+sqr(Fy))*sqr(1+DP)*
        sqrt(sqr(1+X/rho0)+sqr(XP)+sqr(YP));
      DP -= dp_prime*X*tan(e2);
    }

    /* convert CT to distance traveled at final velocity */
    p1 = Po*(1+DP);
    beta1 = p1/sqrt(sqr(p1)+1);
    coord[4] = CT*beta1;
    
    /* transform to cartesian coordinates */
    XP /= (1+X/rho0);
    YP /= (1+X/rho0);

    if (particleLost[i_part] || p1<=0) {
      if (!part[i_top]) {
        fprintf(stdout, "error: couldn't swap particles %ld and %ld--latter is null pointer (track_through_csbend)\n",
                i_part, i_top);
        fflush(stdout);
        abort();
      }
      SWAP_PTR(part[i_part], part[i_top]);
      if (accepted) {
        if (!accepted[i_top]) {
          fprintf(stdout, 
                  "error: couldn't swap acceptance data for particles %ld and %ld--latter is null pointer (track_through_csbend)\n",
                  i_part, i_top);
          fflush(stdout);
          abort();
        }
        SWAP_PTR(accepted[i_part], accepted[i_top]);
      }
      part[i_top][4] = z_start + s_lost;
      part[i_top][5] = Po*(1+part[i_top][5]);
      i_top--;
      i_part--;
      continue;
    }

    if (csbend->edge2_effects) {
      /* apply edge focusing */
      rho = (1+DP)*rho_actual;
      if (csbend->edge_order<2) {
        delta_xp = tan(e2)/rho*X;
        XP += delta_xp;
        YP -= tan(e2-psi2/(1+DP))/rho*Y;
      }
      else 
        apply_edge_effects(&X, &XP, &Y, &YP, rho, n, e2, he2, psi2*(1+DP), 1);
    }

    coord = part[i_part];
    x  =  X*cos_ttilt -  Y*sin_ttilt + dcoord_etilt[0];
    y  =  X*sin_ttilt +  Y*cos_ttilt + dcoord_etilt[2];
    xp = XP*cos_ttilt - YP*sin_ttilt + dcoord_etilt[1];
    yp = XP*sin_ttilt + YP*cos_ttilt + dcoord_etilt[3];
    X  = x;
    Y  = y;
    XP = xp;
    YP = yp;
    coord[0] += dxf + dzf*coord[1];
    coord[2] += dyf + dzf*coord[3];
    coord[4] += dzf*sqrt(1+ sqr(coord[1]) + sqr(coord[3]));
  }

  if (n_part>1 && !csbend->csrBlock) {
    /* prepare more data for CSRDRIFT */
    long imin, imax;
    double S55;
    
    rms_emittance(part, 0, 1, i_top+1, &csrWake.S11, &csrWake.S12, &csrWake.S22);
    rms_emittance(part, 4, 5, i_top+1, &S55, NULL, NULL);
    csrWake.rmsBunchLength = sqrt(S55);
    csrWake.perc68BunchLength = approximateBeamWidth(0.6826, part, i_top+1, 4)/2;
    csrWake.perc90BunchLength = approximateBeamWidth(0.9, part, i_top+1, 4)/2;
#ifdef DEBUG
      fprintf(stderr, "rms bunch length = %le, percentile bunch length (68, 90) = %le, %le\n",
              csrWake.rmsBunchLength, csrWake.perc68BunchLength, csrWake.perc90BunchLength);
#endif
    if (macroParticleCharge) {
      index_min_max(&imin, &imax, csrWake.dGamma, csrWake.bins);
      csrWake.peakToPeakWavelength = 2*fabs(1.0*imax-imin)*dct;
    } else {
      csrWake.peakToPeakWavelength = csrWake.perc68BunchLength;
    }

    csrWake.valid = 1;
    csrWake.rho = rho_actual;
    csrWake.bendingAngle = fabs(angle);
    csrWake.Po = Po;
    csrWake.SGOrder = csbend->SGOrder;
    csrWake.SGDerivOrder = csbend->SGDerivOrder;
    csrWake.SGHalfWidth = csbend->SGHalfWidth;
    csrWake.SGDerivHalfWidth = csbend->SGDerivHalfWidth;
    csrWake.GSConstant = CSRConstant*pow(3*rho0*rho0, 1./3.)/2;  /* used for G. Stupakov's drift formulae */
    csrWake.MPCharge = macroParticleCharge;
    csrWake.binRangeFactor = csbend->binRangeFactor;
    csrWake.trapazoidIntegration = csbend->trapazoidIntegration;
  }
  
#if defined(MINIMIZE_MEMORY)
  /* leave dGamma out of this because that memory is used by CSRDRIFT */
  free(beta0);
  free(ctHist);
  free(ctHistDeriv);
  free(T1);
  free(T2);
  free(denom);
  free(particleLost);
  beta0 = ctHist = ctHistDeriv = T1 = T2 = denom = NULL;
  particleLost  = NULL;
  maxBins = maxParticles = 0;
#endif

  return(i_top+1);
}

long binParticleCoordinate(double **hist, long *maxBins,
                           double *lower, double *upper, double *binSize, long *bins,
                           double expansionFactor,
                           double **particleCoord, long nParticles, long coordinateIndex)
{
  long iBin, iParticle, nBinned;
  double value;
  
  if (*binSize<=0 && *bins<1)
    return -1;
  if (*binSize>0 && *bins>1)
    return -2;
  if (*lower==*upper) {
    /* find range of points */
    *upper = -(*lower = DBL_MAX);
    for (iParticle=0; iParticle<nParticles; iParticle++) {
      value = particleCoord[iParticle][coordinateIndex];
      if (value<*lower)
        *lower = value;
      if (value>*upper)
        *upper = value;
    }
    if (expansionFactor>1) {
      double center, range;
      center = (*lower+*upper)/2;
      range = (*upper-*lower)*expansionFactor;
      *lower = center-range/2;
      *upper = center+range/2;
    }
  }
  
  if (*binSize>0)
    /* bin size given, so determine the number of bins */
    *bins = (*upper-*lower)/(*binSize);
  *binSize = (*upper-*lower)/(*bins);

  /* realloc if necessary */
  if (*bins>*maxBins &&
      !(*hist=SDDS_Realloc(*hist, sizeof(**hist)*(*maxBins=*bins))))
    bomb("Memory allocation failure (binParticleCoordinate)", NULL);
    
  for (iBin=0; iBin<*bins; iBin++)
    (*hist)[iBin] = 0;
  nBinned = 0;
  for (iParticle=nBinned=0; iParticle<nParticles; iParticle++) {
    /* the coordinate of the bin center is (iBin+0.5)*(*binSize) + *lower */
    iBin = (particleCoord[iParticle][coordinateIndex] - *lower)/(*binSize);
    if (iBin<0 || iBin>(*bins-1))
        continue;
    (*hist)[iBin] += 1;
    nBinned++;
  }
  return nBinned;
}

void computeSaldinFdNorm(double **FdNorm, double **x, long *n, double sMax, long ns,
                         double Po, double radius, double angle, double dx, char *normMode);
long track_through_driftCSR_Stupakov(double **part, long np, CSRDRIFT *csrDrift, 
                                     double Po, double **accepted, double zStart, char *rootname);

long track_through_driftCSR(double **part, long np, CSRDRIFT *csrDrift, 
                            double Po, double **accepted, double zStart, char *rootname)
{
  long iPart, iKick, iBin, binned, nKicks, iSpreadMode=0;
  double *coord, p, beta, dz, ct0=0.0, factor, dz0, dzFirst;
  double ctmin, ctmax, spreadFactor, dct;
  double zTravel, attenuationLength, thetaRad=0.0, sigmaZ, overtakingLength, criticalWavelength, wavelength=0.0;
  static char *spreadMode[3] = {"full", "simple", "radiation-only"};
  static char *wavelengthMode[3] = {"sigmaz", "bunchlength", "peak-to-peak"};
  static char *bunchlengthMode[3] = {"rms", "68-percentile", "90-percentile"};
  unsigned long mode;
  static long warned = 0;
  long nBins1;
  
  if (np<=1 || !csrWake.valid || !csrDrift->csr) {
    exactDrift(part, np, csrDrift->length);
    return np;
  }
  nBins1 = csrWake.bins - 1;

  mode = 
    (csrDrift->spread?CSRDRIFT_SPREAD:0) +
      (csrDrift->useOvertakingLength?CSRDRIFT_OVERTAKINGLENGTH:0) +
        (csrDrift->useSaldin54?CSRDRIFT_SALDIN54:0) +
          (csrDrift->attenuationLength>0?CSRDRIFT_ATTENUATIONLENGTH:0) +
            (csrDrift->useStupakov?CSRDRIFT_STUPAKOV:0) ;
  if (bitsSet(mode)>1) {
    fprintf(stdout, "Error: Too many modes set for CSRDRIFT.\n");
    exit(1);
  }
  if (csrWake.lastMode && csrWake.lastMode!=mode) {
    fprintf(stdout, "Error: CSRDRIFT mode changed between dipoles. Pick one mode following each dipole.\n");
    exit(1);
  }
  csrWake.lastMode = mode;
  
  if (mode&CSRDRIFT_STUPAKOV)
    return track_through_driftCSR_Stupakov(part, np, csrDrift, Po, accepted, zStart, rootname);

  if (!warned) {
    fprintf(stdout, "Warning: USE_STUPAKOV=1 is recommended for CSRDRIFT elements.\n");
    fprintf(stdout, "This is the most physical model available at this time in elegant.\n");
    warned = 1;
  }
  
  dct = csrWake.dctBin;
  if (csrDrift->dz>0) {
    if ((nKicks = csrDrift->length/csrDrift->dz)<1)
      nKicks = 1;
  } else 
    nKicks = csrDrift->nKicks;
  if (nKicks<=0)
    bomb("nKicks=0 in CSR drift.", NULL);
  dz = (dz0=csrDrift->length/nKicks)/2;
  
  sigmaZ = 0;
  switch (match_string(csrDrift->bunchlengthMode, bunchlengthMode, 3, 0)) {
  case 0:
    sigmaZ = csrWake.rmsBunchLength;
    break;
  case 1:
    sigmaZ = csrWake.perc68BunchLength;
    break;
  case 2:
    sigmaZ = csrWake.perc90BunchLength;
    break;
  default:
    bomb("invalid bunchlength_mode for CSRDRIFT.  Use rms or percentile.", NULL);
  }
  
  overtakingLength = pow(24*sigmaZ*csrWake.rho*csrWake.rho, 1./3.);

  if (mode&CSRDRIFT_OVERTAKINGLENGTH)
    attenuationLength = overtakingLength*csrDrift->overtakingLengthMultiplier;
  else
    attenuationLength = csrDrift->attenuationLength;
  
  if (mode&CSRDRIFT_SPREAD) {
    iSpreadMode = 0;
    if (csrDrift->spreadMode && 
        (iSpreadMode=match_string(csrDrift->spreadMode, spreadMode, 3, 0))<0)
      bomb("invalid spread_mode for CSR DRIFT.  Use full, simple, or radiation-only", NULL);
    switch (match_string(csrDrift->wavelengthMode, wavelengthMode, 3, 0)) {
    case 0:
    case 1:
      /* bunch length */
      wavelength = sigmaZ;
      break;
    case 2:
      /* peak-to-peak */
      wavelength = csrWake.peakToPeakWavelength;
      break;
    default:
      bomb("invalid wavelength_mode for CSR DRIFT.  Use sigmaz or peak-to-peak", NULL);
      break;
    }
    criticalWavelength = 4.19/ipow(csrWake.Po, 3)*csrWake.rho;
    thetaRad = 0.5463e-3/(csrWake.Po*0.511e-3)/pow(criticalWavelength/wavelength, 1./3.);
  }

  if (mode&CSRDRIFT_SALDIN54) {
    if (csrWake.FdNorm==NULL) {
      if (csrDrift->nSaldin54Points<20) 
        csrDrift->nSaldin54Points = 20;
      computeSaldinFdNorm(&csrWake.FdNorm, &csrWake.xSaldin, &csrWake.nSaldin,
                          2*sigmaZ, csrDrift->nSaldin54Points, csrWake.Po, csrWake.rho, csrWake.bendingAngle, dz,
                          csrDrift->normMode);
      if (csrDrift->Saldin54Output)  {
        long ix;
        if (!csrDrift->fpSaldin) {
          csrDrift->Saldin54Output = compose_filename(csrDrift->Saldin54Output, rootname);
          csrDrift->fpSaldin = fopen(csrDrift->Saldin54Output, "w");
          fprintf(csrDrift->fpSaldin, "SDDS1\n&column name=z, type=double &end\n&column name=Factor, type=double &end\n");
          fprintf(csrDrift->fpSaldin, "&data mode=ascii no_row_counts=1 &end\n");
        } else
          fprintf(csrDrift->fpSaldin, "\n");
        for (ix=0; ix<csrWake.nSaldin; ix++) 
          fprintf(csrDrift->fpSaldin, "%le %le\n", csrWake.xSaldin[ix], csrWake.FdNorm[ix]);
        fflush(csrDrift->fpSaldin);
      }
    }
  }

  dzFirst = zStart - csrWake.zLast;
  zTravel = zStart-csrWake.z0;  /* total distance traveled by radiation to reach this point */
#ifdef DEBUG
  fprintf(stdout, "CSR in drift:\n");
  fprintf(stdout, "zStart = %21.15le, zLast = %21.15le, zTravel = %21.15le\n", zStart, csrWake.zLast,
          zTravel);
  fprintf(stdout, "dzFirst = %21.15e, s0 = %21.15e\n", dzFirst, csrWake.s0);
#endif

  for (iKick=0; iKick<nKicks; iKick++) {
    /* first drift is dz=dz0/2, others are dz0 */
    if (iKick==1)
      dz = dz0;
    zTravel += dz;

    ctmin = DBL_MAX;
    ctmax = -DBL_MAX;

    /* propagate particles forward, converting s to c*t=s/beta */
    for (iPart=0; iPart<np; iPart++) {
      coord = part[iPart];
      coord[0] += coord[1]*dz;
      coord[2] += coord[3]*dz;
      p = Po*(1+coord[5]);
      beta = p/sqrt(p*p+1);
      coord[4] = (coord[4]+dz*sqrt(1+sqr(coord[1])+sqr(coord[3])))/beta;
#ifdef DEBUG
      if (coord[4]>ctmax)
        ctmax = coord[4];
      if (coord[4]<ctmin)
        ctmin = coord[4];
#endif
    }

    factor = 1;
    if (csrWake.dGamma) {
      /* propagate wake forward */
      csrWake.s0 += dz+dzFirst;   /* accumulates position of back end of the radiation pulse */
      ct0 = csrWake.s0;
      
      if (attenuationLength>0) {
        /* attenuate wake */
        if ((factor = exp(-(dz+dzFirst)/attenuationLength))<1) {
          for (iBin=0; iBin<csrWake.bins; iBin++)
            csrWake.dGamma[iBin] *= factor;
        }
      }
      /* factor to account for difference in drift lengths here and in
       * csrcsbend integration.  Use dz0 here because that is the
       * length integrated by each kick.  Add dzFirst to account for any
       * length we may have missed due to intervening non-drift elements.
       */
      factor = (dz0+dzFirst)/csrWake.ds0;
    }
    if (mode&CSRDRIFT_SPREAD) {
      /* compute loss of on-axis field due to spread of beam using a simple-minded
       * computation of beam sizes */
      switch (iSpreadMode) {
      case 0:  /* full */
        factor *= (spreadFactor =
                   sqrt(csrWake.S11/(csrWake.S11 + 
                                     2*zTravel*csrWake.S12 + 
                                     zTravel*zTravel*(sqr(thetaRad)+csrWake.S22))));
        break;
      case 1: /* simple */
        factor *= (spreadFactor =
                   sqrt(csrWake.S11/(csrWake.S11 + zTravel*zTravel*(sqr(thetaRad)+csrWake.S22))));
        break;
      case 2: /* radiation only */
        factor *= (spreadFactor =
                   sqrt(csrWake.S11/(csrWake.S11 + sqr(zTravel*thetaRad))));
        break;
      default:
        bomb("invalid spread code---programming error!", NULL);
        break;
      }
    }
    
    if (mode&CSRDRIFT_SALDIN54) {
      long code;
      double f0 = 0;
      if (zTravel<=csrWake.xSaldin[csrWake.nSaldin-1]) 
        factor *= (f0=interp(csrWake.FdNorm, csrWake.xSaldin, csrWake.nSaldin, zTravel, 0, 1, &code));
      else 
        factor = 0;
      csrWake.lastFdNorm = f0;
#ifdef DEBUG
      fprintf(csrWake.fpSaldin, "%le %le\n", zTravel, f0);
      fflush(csrWake.fpSaldin);
#endif
      if (!code) {
        fprintf(stderr, "Warning: interpolation failure for Saldin eq. 54\n");
        fprintf(stderr, "zTravel = %le,  csrWake available up to %le\n",
                zTravel, csrWake.xSaldin[csrWake.nSaldin-1]);
        factor = 0;
      }
    }
    
    dzFirst = 0;

    /* apply kick to each particle and convert back to normal coordinates */
    for (iPart=binned=0; iPart<np; iPart++) {
      coord = part[iPart];
      if (csrWake.dGamma) {
        double f;
        iBin = (f=(coord[4]-ct0)/dct);
        f -= iBin;
        if (iBin>=0 && iBin<nBins1) {
          coord[5] += ((1-f)*csrWake.dGamma[iBin]+f*csrWake.dGamma[iBin+1])/Po*factor;
          binned ++;
        }
      }
      p = (1+coord[5])*Po;
      beta = p/sqrt(p*p+1);
      coord[4] = beta*coord[4];
    }
    if (csrWake.dGamma && np!=binned) {
      fprintf(stdout, "only %ld of %ld particles binned for CSR drift\n",
              binned, np);
#ifdef DEBUG
      fprintf(stdout, "beam ct min, max = %21.15e, %21.15e\n",
              ctmin, ctmax);
      fprintf(stdout, "wake ct0 = %21.15e, ct1 = %21.15e\n",
              ct0, ct0+csrWake.dctBin*csrWake.bins);
#endif
      fflush(stdout);
    }
  }
  /* do final drift of dz0/2 */
  dz = dz0/2;
  for (iPart=0; iPart<np; iPart++) {
    coord = part[iPart];
    coord[0] += coord[1]*dz;
    coord[2] += coord[3]*dz;
    coord[4] += dz*sqrt(1+sqr(coord[1])+sqr(coord[3]));
  }    

  csrWake.zLast = zStart+csrDrift->length;
  
  if (csrWake.dGamma) {
    /* propagate wake forward */
    csrWake.s0 += dz;
    ct0 = csrWake.s0;
    
    if (attenuationLength>0) {
      /* attenuate wake */
      if ((factor = exp(-dz/attenuationLength))<1) {
        for (iBin=0; iBin<csrWake.bins; iBin++)
            csrWake.dGamma[iBin] *= factor;
        }
    }
  }

  return np;
}

/* this should be called before starting to track a beamline to make sure that
 * CSR drift elements upstream of all CSRBEND elements get treated like ordinary
 * drift spaces. */

long reset_driftCSR()
{
  csrWake.lastMode = 0;
  if (csrWake.valid && csrWake.FdNorm) {
    fprintf(stdout, "Last value of normalization factor for CSR wake was %le\n",
            csrWake.lastFdNorm);
  }
  csrWake.valid = csrWake.bins = 0;
  csrWake.dctBin = csrWake.s0 = csrWake.ds0 = csrWake.zLast =
    csrWake.z0 = csrWake.S11 = csrWake.S12 = csrWake.S22 = 0;
  csrWake.dGamma = NULL;
  csrWake.nSaldin = 0;
  if (csrWake.FdNorm) {
    free(csrWake.FdNorm);
    free(csrWake.xSaldin);
    csrWake.FdNorm = csrWake.xSaldin = NULL;
  }
  if (csrWake.StupakovFileActive) {
    if (!SDDS_Terminate(&csrWake.SDDS_Stupakov))
      bomb("problem terminating data file for Stupakov output from CSRDRIFT", NULL);
    csrWake.StupakovFileActive = 0;
  }
  return 1;
}

double SolveForPsiSaldin54(double xh, double sh);
double Saldin5354Factor(double xh, double sh, double phihm, double xhLowerLimit);

void computeSaldinFdNorm(double **FdNorm, double **x, long *n, double sMax, long ns,
                         double Po, double radius, double bendingAngle, double dx, 
                         char *normMode)
{
  double xEnd, sh, beta, gamma, xh, dx0;
  long ix, is;
  double phihs, phihm, xhLowerLimit, xUpperLimit, s, f, fx;
  double t1, t2, f0, fmax;
  char *allowedNormMode[2] = {"first", "peak"};

  gamma = sqrt(sqr(Po)+1);
  beta = Po/gamma;

  if ((xEnd = sMax/(1-beta))>1000 || isnan(xEnd) || isinf(xEnd)) {
    fprintf(stderr, "Warning: the extent of the CSR drift wake decay was limited at 1km\n");
    xEnd = 1000;
  }

  *n = 100;
  dx0 = xEnd/(100*(*n));
  if (dx<dx0) {
    *n = xEnd/(100*dx);
    if (*n>100000) {
      *n = 100000;
      fprintf(stderr, "Note: the CSR drift wake decay table size hit the limit of 100k points\n");
    }
  } else 
    dx = dx0;
  fx = pow(xEnd/dx, 1./(*n));

  if (!(*FdNorm = calloc(sizeof(**FdNorm), (*n))) ||
      !(*x = malloc(sizeof(**x)*(*n))))
    bomb("memory allocation failure (computeSaldinFdNorm)", NULL);

  for (ix=0; ix<*n; ix++)
    (*x)[ix] = ix==0 ? 0 : ipow(fx, ix-1)*dx;
  for (is=0; is<ns; is++) {
    /* don't use s=0 as it is singular */
    s = (is+1.0)*sMax/ns;
    sh = s*ipow(gamma, 3)/radius;
    phihm = bendingAngle*gamma;
    t1 = 12*sh;
    t2 = sqrt(64+144*sh*sh);
    phihs = pow(t1+t2, 1./3.) - pow(-t1 + t2, 1./3.);
    xhLowerLimit = -1;
    if (phihs>phihm)
      xhLowerLimit = sh - phihm - ipow(phihm, 3)/6 + sqrt(sqr(ipow(phihm, 3)-6*sh) + 9*ipow(phihm, 4))/6;
    xUpperLimit = 0.999*s/(1-beta);
    for (ix=0; ix<*n; ix++) {
      if ((*x)[ix]>=xUpperLimit)
        break;
      xh = (*x)[ix]*gamma/radius;
      (*FdNorm)[ix] += Saldin5354Factor(xh, sh, phihm, xhLowerLimit);
    }
  }

  /* average over s */
  for (ix=0; ix<*n; ix++)
    (*FdNorm)[ix] /= ns;
  
  /* get the first nonzero and also the maximum value of Fd */
  for (ix=f0=fmax=0; ix<*n; ix++) {
    f= (*FdNorm)[ix];
    if (f0==0 && f>0)
      f0 = f;
    if (fmax<f)
      fmax = f;
  }
  if (fmax>f0/0.99) {
    fprintf(stderr, "Warning: possible problem with SALDIN54 drift mode: too few (%ld) points. Max/start-1 is %le\n",
            ns,
            fmax/f0-1);
  }
  switch (match_string(normMode, allowedNormMode, 2, 0)) {
  case 0:
    /* first */
    f = f0;
    break;
  case 1:
    /* peak */
    f = fmax;
    break;
  default:
    fprintf(stderr, "Error: unknown Saldin-54 normalization mode: %s\n", normMode);
    exit(1);
    break;
  }
  if (f)
    for (ix=0; ix<*n; ix++)
      (*FdNorm)[ix] /= f;
  else
    for (ix=0; ix<*n; ix++)
      (*FdNorm)[ix] = 0;
}

double SolveForPsiSaldin54(double xh, double sh)
{
  double s_sum, s_diff2, bestSol;
  double solList[4] = {-1, -1, -1, -1};
  long nSols=0, sol;

  s_sum = (-2*xh - sqrt(-8 + 4*pow(xh,2) - 
                         (4*pow(2,0.3333333333333333)*(-1 + pow(xh,2)))/
                         pow(2 + 9*pow(sh,2) - 3*pow(xh,2) - 6*sh*pow(xh,3) + 
                               3*pow(xh,4) + sqrt(4*pow(-1 + pow(xh,2),3) + 
                                                    pow(2 + 9*pow(sh,2) - 3*pow(xh,2) - 6*sh*pow(xh,3) + 
                                                          3*pow(xh,4),2)),0.3333333333333333) + 
                         2*pow(2,0.6666666666666666)*
                         pow(2 + 9*pow(sh,2) - 3*pow(xh,2) - 6*sh*pow(xh,3) + 
                               3*pow(xh,4) + sqrt(4*pow(-1 + pow(xh,2),3) + 
                                                    pow(2 + 9*pow(sh,2) - 3*pow(xh,2) - 6*sh*pow(xh,3) + 
                                                          3*pow(xh,4),2)),0.3333333333333333)))/2.;
  if (!isnan(s_sum)) {
    s_diff2 = (-16 + 8*pow(xh,2) + (4*pow(2,0.3333333333333333)*(-1 + pow(xh,2)))/
              pow(2 + 9*pow(sh,2) - 3*pow(xh,2) - 6*sh*pow(xh,3) + 
                    3*pow(xh,4) + sqrt(4*pow(-1 + pow(xh,2),3) + 
                                         pow(2 + 9*pow(sh,2) - 3*pow(xh,2) - 6*sh*pow(xh,3) + 
                                               3*pow(xh,4),2)),0.3333333333333333) - 
              2*pow(2,0.6666666666666666)*
              pow(2 + 9*pow(sh,2) - 3*pow(xh,2) - 6*sh*pow(xh,3) + 
                    3*pow(xh,4) + sqrt(4*pow(-1 + pow(xh,2),3) + 
                                         pow(2 + 9*pow(sh,2) - 3*pow(xh,2) - 6*sh*pow(xh,3) + 
                                               3*pow(xh,4),2)),0.3333333333333333) + 
              (16*(-3*sh + pow(xh,3)))/
              sqrt(-8 + 4*pow(xh,2) - 
                   (4*pow(2,0.3333333333333333)*(-1 + pow(xh,2)))/
                   pow(2 + 9*pow(sh,2) - 3*pow(xh,2) - 6*sh*pow(xh,3) + 
            3*pow(xh,4) + sqrt(4*pow(-1 + pow(xh,2),3) + 
                                 pow(2 + 9*pow(sh,2) - 3*pow(xh,2) - 6*sh*pow(xh,3) + 
                                       3*pow(xh,4),2)),0.3333333333333333) + 
                   2*pow(2,0.6666666666666666)*
                   pow(2 + 9*pow(sh,2) - 3*pow(xh,2) - 6*sh*pow(xh,3) + 
                         3*pow(xh,4) + sqrt(4*pow(-1 + pow(xh,2),3) + 
                                              pow(2 + 9*pow(sh,2) - 3*pow(xh,2) - 6*sh*pow(xh,3) + 
                                                    3*pow(xh,4),2)),0.3333333333333333)))/4.;
    if (s_diff2>=0) {
      solList[0] = s_sum+sqrt(s_diff2);
      solList[1] = s_sum+sqrt(s_diff2);
      nSols = 2;
    }
  }
  
  s_sum =    (-2*xh + sqrt(-8 + 4*pow(xh,2) - 
                            (4*pow(2,0.3333333333333333)*(-1 + pow(xh,2)))/
                            pow(2 + 9*pow(sh,2) - 3*pow(xh,2) - 6*sh*pow(xh,3) + 
                                  3*pow(xh,4) + sqrt(4*pow(-1 + pow(xh,2),3) + 
                                                       pow(2 + 9*pow(sh,2) - 3*pow(xh,2) - 6*sh*pow(xh,3) + 
                                                             3*pow(xh,4),2)),0.3333333333333333) + 
                            2*pow(2,0.6666666666666666)*
                            pow(2 + 9*pow(sh,2) - 3*pow(xh,2) - 6*sh*pow(xh,3) + 
                                  3*pow(xh,4) + sqrt(4*pow(-1 + pow(xh,2),3) + 
                                                       pow(2 + 9*pow(sh,2) - 3*pow(xh,2) - 6*sh*pow(xh,3) + 
                                                             3*pow(xh,4),2)),0.3333333333333333)))/2.;
  if (!isnan(s_sum)) {
    s_diff2 = (-16 + 8*pow(xh,2) + (4*pow(2,0.3333333333333333)*(-1 + pow(xh,2)))/
              pow(2 + 9*pow(sh,2) - 3*pow(xh,2) - 6*sh*pow(xh,3) + 
                    3*pow(xh,4) + sqrt(4*pow(-1 + pow(xh,2),3) + 
                                         pow(2 + 9*pow(sh,2) - 3*pow(xh,2) - 6*sh*pow(xh,3) + 
                                               3*pow(xh,4),2)),0.3333333333333333) - 
              2*pow(2,0.6666666666666666)*
              pow(2 + 9*pow(sh,2) - 3*pow(xh,2) - 6*sh*pow(xh,3) + 
                    3*pow(xh,4) + sqrt(4*pow(-1 + pow(xh,2),3) + 
                                         pow(2 + 9*pow(sh,2) - 3*pow(xh,2) - 6*sh*pow(xh,3) + 
                                               3*pow(xh,4),2)),0.3333333333333333) - 
              (16*(-3*sh + pow(xh,3)))/
              sqrt(-8 + 4*pow(xh,2) - 
                   (4*pow(2,0.3333333333333333)*(-1 + pow(xh,2)))/
                   pow(2 + 9*pow(sh,2) - 3*pow(xh,2) - 6*sh*pow(xh,3) + 
                         3*pow(xh,4) + sqrt(4*pow(-1 + pow(xh,2),3) + 
                                              pow(2 + 9*pow(sh,2) - 3*pow(xh,2) - 6*sh*pow(xh,3) + 
                                                    3*pow(xh,4),2)),0.3333333333333333) + 
                   2*pow(2,0.6666666666666666)*
                   pow(2 + 9*pow(sh,2) - 3*pow(xh,2) - 6*sh*pow(xh,3) + 
                         3*pow(xh,4) + sqrt(4*pow(-1 + pow(xh,2),3) + 
                                              pow(2 + 9*pow(sh,2) - 3*pow(xh,2) - 6*sh*pow(xh,3) + 
                                                    3*pow(xh,4),2)),0.3333333333333333)))/4.;
  
    if (s_diff2>=0) {
      solList[nSols] = s_sum+sqrt(s_diff2);
      solList[nSols+1] = s_sum-sqrt(s_diff2);
      nSols += 2;
    }
  }
  bestSol = solList[0];
  for (sol=0; sol<nSols; sol++) {
    if (solList[sol]>bestSol) {
      bestSol = solList[sol];
    }     
  }
  return bestSol;
}

double Saldin5354Factor(double xh, double sh, double phihm, double xhLowerLimit)
{
  double t1, t2, f, psi, psi2;
  if (xh<xhLowerLimit) {
    /* use Saldin 53 */
    t1 = (ipow(phihm, 3) + 3*xh*sqr(phihm) - 6*sh);
    t2 = 3*(phihm+2*xh);
    f = 2/(phihm+2*xh)*(1 + (t1 + t2)/sqrt(t1*t1+sqr(phihm*t2))) - 1/sh;
  } else {
    if ((psi = SolveForPsiSaldin54(xh, sh))>=0) {
      psi2 = psi*psi;
      f =  4*(2*xh*(psi2+1)+psi*(psi2+2))/
        (4*xh*xh*(psi2+1)+4*xh*psi*(psi2+2)+psi2*(psi2+4)) - 1/sh;
    } else
      return 0;
  }
  if (isnan(f) || isinf(f))
    f = 0;
  return f;
}

void exactDrift(double **part, long np, double length)
{
  long i;
  double *coord;
  for (i=0; i<np; i++) {
    coord = part[i];
    coord[0] += coord[1]*length;
    coord[2] += coord[3]*length;
    coord[4] += length*sqrt(1+sqr(coord[1])+sqr(coord[3]));
  }
}


double SolveForPhiStupakov(double x, double ds, double phim);
void DumpStupakovOutput(char *filename, SDDS_DATASET *SDDSout, long *active,
                        double zTravel, double *ctHist, double *ctHistDeriv,
                        double *dGamma, long nBins, double dct, 
                        double MPCharge, double dz,
                        long nCaseC, long nCaseD1,long nCaseD2,
                        double x, double dsMax, double phi0, double phi1) ;

static double SolveForPhiStupakovDiffSum = 0;
static long SolveForPhiStupakovDiffCount = 0;

long track_through_driftCSR_Stupakov(double **part, long np, CSRDRIFT *csrDrift, 
                            double Po, double **accepted, double zStart, char *rootname)
{
  long iPart, iKick, iBin, binned, nKicks;
  long nCaseC, nCaseD1, nCaseD2;
  double ctLower, ctUpper, ds;
  long nBins, maxBins, nBinned, diBin;
  double *coord, p, beta, dz, factor, dz0, dzFirst;
  double zTravel, dct, zOutput;
  double *ctHist=NULL, *ctHistDeriv=NULL, *phiSoln=NULL;
  double length;
  long nBins1;
  double dsMax, x;
  TRACKING_CONTEXT tContext;

  getTrackingContext(&tContext);

  SolveForPhiStupakovDiffCount = 0;
  SolveForPhiStupakovDiffSum = 0;
  
  length = csrDrift->length;
  if (zStart!=csrWake.zLast) {
    length += (dzFirst = zStart-csrWake.zLast);
    /* propagate beam back so we can tranverse the missing length including CSR
     */
    for (iPart=0; iPart<np; iPart++) {
      coord = part[iPart];
      coord[0] -= dzFirst*coord[1];
      coord[2] -= dzFirst*coord[3];
      coord[4] -= dzFirst*sqrt(1+sqr(coord[1])+sqr(coord[3]));
    }
    zStart = csrWake.zLast;
  }
  zOutput = zStart;  /* absolute coordinate used for output of data vs z or s */
  
  if (csrDrift->dz>0) {
    if ((nKicks = length/csrDrift->dz+0.5)<1)
      nKicks = 1;
  } else 
    nKicks = csrDrift->nKicks;
  if (nKicks<=0)
    bomb("nKicks=0 in CSR drift.", NULL);
  dz = (dz0=length/nKicks)/2;
  
  zTravel = zStart-csrWake.z0;  /* total distance traveled by radiation to reach this point */

  maxBins = nBins = csrWake.bins;
  nBins1 = nBins-1;
  if (!(ctHist=SDDS_Malloc(sizeof(*ctHist)*nBins)) ||
      !(ctHistDeriv=SDDS_Malloc(sizeof(*ctHistDeriv)*nBins)) ||
      !(phiSoln=SDDS_Malloc(sizeof(*phiSoln)*nBins)))
    bomb("memory allocation failure (track_through_driftCSR)", NULL);
  
  for (iKick=0; iKick<nKicks; iKick++) {
    /* first drift is dz=dz0/2, others are dz0 */
    if (iKick==1)
      dz = dz0;
    zTravel += dz;
    zOutput += dz;
    
    x = zTravel/csrWake.rho;
    dsMax = csrWake.rho/24*pow(csrWake.bendingAngle, 3)
      *(csrWake.bendingAngle+4*x)/(csrWake.bendingAngle+x);
    /* propagate particles forward, converting s to c*t=s/beta */
    for (iPart=0; iPart<np; iPart++) {
      coord = part[iPart];
      coord[0] += coord[1]*dz;
      coord[2] += coord[3]*dz;
      p = Po*(1+coord[5]);
      beta = p/sqrt(p*p+1);
      coord[4] = (coord[4]+dz*sqrt(1+sqr(coord[1])+sqr(coord[3])))/beta;
    }

    /* bin the particle distribution */
    ctLower = ctUpper = dct = 0;
    if ((nBinned = 
         binParticleCoordinate(&ctHist, &maxBins,
                               &ctLower, &ctUpper, &dct, &nBins, 
                               csrWake.binRangeFactor<1.1?1.1:csrWake.binRangeFactor,
                               part, np, 4))!=np) {
      fprintf(stdout, "Only %ld of %ld particles binned for CSRDRIFT\n", nBinned, np);
      fflush(stdout);
    }
      
    /* - smooth the histogram, normalize to get linear density, and 
       copy in preparation for taking derivative
       */
    if (csrWake.highFrequencyCutoff0>0) {
      long nz;
      nz = applyLowPassFilter(ctHist, nBins, dct, csrWake.highFrequencyCutoff0, csrWake.highFrequencyCutoff1);
      if (nz) {
        fprintf(stdout, "Warning: high pass filter resulted in negative values in %ld bins\n",
                nz);
        fflush(stdout);
      }
    }
    if (csrWake.SGHalfWidth>0) {
      SavitzkyGolaySmooth(ctHist, nBins, csrWake.SGOrder, csrWake.SGHalfWidth, csrWake.SGHalfWidth,  0);
      correctDistribution(ctHist, nBins, 1.0*nBinned);
    }
    for (iBin=0; iBin<nBins; iBin++)
      ctHistDeriv[iBin] = (ctHist[iBin] /= dct);
    /* - compute derivative with smoothing.  The deriv is w.r.t. index number and
     * I won't scale it now as it will just fall out in the integral 
     */
    SavitzkyGolaySmooth(ctHistDeriv, nBins, csrWake.SGDerivOrder, 
                       csrWake.SGDerivHalfWidth, csrWake.SGDerivHalfWidth, 1);

    /* Case C */ 
    nCaseC = 0;
    nCaseD1 = 0;
    nCaseD2 = 0;
    for (iBin=0; iBin<nBins; iBin++) {
      double f;
      ds = csrWake.rho/6*sqr(csrWake.bendingAngle)*(csrWake.bendingAngle + 3*x);
      diBin = ds/dct;
      if (iBin+diBin<nBins) {
        f = -1/(csrWake.bendingAngle+2*x); 
        csrWake.dGamma[iBin] = f*ctHist[iBin+diBin];
        nCaseC++;
      } else
        csrWake.dGamma[iBin] = 0;
    }
    /* Case D */
    for (iBin=0; iBin<nBins; iBin++) {
      phiSoln[iBin] = -1;
      if ((ds = iBin*dct)>dsMax)
        break;
      phiSoln[iBin] = SolveForPhiStupakov(x, iBin*dct/csrWake.rho, csrWake.bendingAngle);
    }
    for (iBin=0; iBin<nBins; iBin++) {
      long jBin, first, count;
      double term1, term2;
      diBin = dsMax/dct;
      if (iBin+diBin<nBins) {
        nCaseD1 ++;
        csrWake.dGamma[iBin] += ctHist[iBin+diBin]/(csrWake.bendingAngle+2*x);
      }
      first = 1;
      count = 0;
      for (jBin=iBin; jBin<nBins; jBin++) {
        double phi, denom;
        if ((phi = phiSoln[jBin-iBin])>=0) {
          /* I put in a negative sign here because my s is opposite in direction to 
           * Saldin et al. and Stupakov, so my derivative has the opposite sign.
           * Note lack of ds factor here as I use the same one in my unnormalized derivative.
           */
          if (phi>0) {
            /* ^^^ If I test phi+2*x here, I get noisy, unphysical results very close
             * to the dipole exit 
             */
            term2 = ctHistDeriv[jBin]/(phi+2*x);
            csrWake.dGamma[iBin] -= term2;
            if (first) {
              term1 = term2;
              first = 0;
            }
            count++;
            nCaseD2++;
          }
        } else
          break;
      }
      if (count>1 && csrWake.trapazoidIntegration)
        /* trapazoid rule correction for ends */
        csrWake.dGamma[iBin] += (term1+term2)/2;
    }
    /* the minus sign adjusts for Stupakov using wake<0 to indicate energy gain
     */
    factor = -4/csrWake.rho*csrWake.GSConstant*dz0;
    for (iBin=0; iBin<nBins; iBin++)
      csrWake.dGamma[iBin] *= factor;

    if ((csrDrift->StupakovOutput || csrWake.StupakovFileActive) && 
        (csrDrift->StupakovOutputInterval<2 || iKick%csrDrift->StupakovOutputInterval==0)) {
      double x, dsMax, phi0, phi1;
      if (!csrWake.StupakovFileActive) {
        if (!SDDS_CopyString(&csrWake.StupakovOutput, csrDrift->StupakovOutput))
          bomb("string copying problem preparing Stupakov output for CSRDRIFT", NULL);
        csrWake.StupakovOutput = compose_filename(csrWake.StupakovOutput, rootname);
      }
      x = zTravel/csrWake.rho;
      dsMax = csrWake.rho/24*pow(csrWake.bendingAngle, 3)
        *(csrWake.bendingAngle+4*x)/(csrWake.bendingAngle+x);
      phi0 = SolveForPhiStupakov(x, 0.0, csrWake.bendingAngle);
      phi1 = SolveForPhiStupakov(x, dsMax/csrWake.rho*0.999, csrWake.bendingAngle);
      
      /* note that the contents of ctHist and ctHistDeriv are corrupted by this operation */
      DumpStupakovOutput(csrWake.StupakovOutput, &csrWake.SDDS_Stupakov, 
                         &csrWake.StupakovFileActive, zTravel,
                         ctHist, ctHistDeriv, csrWake.dGamma, nBins, dct, csrWake.MPCharge,
                         dz0, nCaseC, nCaseD1, nCaseD2,
                         x, dsMax/csrWake.rho, phi0, phi1);
    }
    
    /* apply kick to each particle and convert back to normal coordinates */
    for (iPart=binned=0; iPart<np; iPart++) {
      double f;
      coord = part[iPart];
      iBin = (f=(coord[4]-ctLower)/dct);
      f -= iBin;
      if (iBin>=0 && iBin<nBins1) {
        coord[5] += ((1-f)*csrWake.dGamma[iBin] + f*csrWake.dGamma[iBin+1])/Po;
        binned ++;
      }
      p = (1+coord[5])*Po;
      beta = p/sqrt(p*p+1);
      coord[4] = beta*coord[4];
    }
    if (tContext.sliceAnalysis && tContext.sliceAnalysis->active &&
	(csrDrift->sliceAnalysisInterval==0 ||
	 iKick%csrDrift->sliceAnalysisInterval==0))  
	performSliceAnalysisOutput(tContext.sliceAnalysis, part, np, 
				   0, tContext.step, Po, 
				   csrWake.MPCharge*np,
				   tContext.elementName, 
				   zOutput, 0);
    if (np!=binned) {
      fprintf(stdout, "only %ld of %ld particles binned for CSR drift\n",
              binned, np);
      fflush(stdout);
    }
  }
  
  /* do final drift of dz0/2 */
  dz = dz0/2;
  for (iPart=0; iPart<np; iPart++) {
    coord = part[iPart];
    coord[0] += coord[1]*dz;
    coord[2] += coord[3]*dz;
    coord[4] += dz*sqrt(1+sqr(coord[1])+sqr(coord[3]));
  }    

  csrWake.zLast = zStart + length;
  free(ctHist);
  free(ctHistDeriv);
  free(phiSoln);
#if DEBUG
  if (SolveForPhiStupakovDiffCount)
    fprintf(stdout, "Phi solution accuracy for %ld solutions: %le\n",
            SolveForPhiStupakovDiffCount, SolveForPhiStupakovDiffSum/SolveForPhiStupakovDiffCount);
#endif
  return np;
}

static double SolveForPhiStupakov_x, SolveForPhiStupakov_4x;

double SolveForPhiStupakovFn(double phi)
{
  return phi*phi*phi*(phi+SolveForPhiStupakov_4x)/(phi+SolveForPhiStupakov_x);
}

/* solve for phi:  ds=phi^3/24*(phi+4*x)/(phi+x), where ds = (s-s')/rho */

double SolveForPhiStupakov(double x, double ds, double phim)
{
  double phi;
  static double phiLast = -1;
  
  if (ds<0)
    return -1;
  if (ds==0)
    return 0;
  
  ds *= 24;
  SolveForPhiStupakov_x = x;
  SolveForPhiStupakov_4x = 4*x;

  if (phiLast==-1)
    phiLast = phim/2;

  /* try phim first */
  if (fabs(ds-SolveForPhiStupakovFn(phim))<ds/1e4) {
    phiLast = phim;
    return phim;
  }
  
  /* try a solution with Newton's method */
  phi = zeroNewton(SolveForPhiStupakovFn, ds, phiLast, phim/1000, 3, ds/1e4);
  if (phi<0 || phi>phim || fabs(ds - SolveForPhiStupakovFn(phi))>ds/1e4) 
    /* try a more plodding method */
    phi = zeroInterp(SolveForPhiStupakovFn, ds, 0, phim*1.01, phim/100, ds/1e4);
  if (phi<0 || phi>phim)
    return -1;
  phiLast = phi;
  SolveForPhiStupakovDiffCount ++;
  SolveForPhiStupakovDiffSum += fabs(ds - SolveForPhiStupakovFn(phi));
  return phi;
}


/* this procedure destroys the contents of ctHist and ctHistDeriv ! */

void DumpStupakovOutput(char *filename, SDDS_DATASET *SDDSout, long *active,
                        double zTravel, double *ctHist, double *ctHistDeriv,
                        double *dGamma, long nBins, double dct, 
                        double MPCharge, double dz,
                        long nCaseC, long nCaseD1, long nCaseD2,
                        double x, double dsMax, double phi0, double phi1) 
{
  long i;
  if (!*active) {
    if (!SDDS_InitializeOutput(SDDSout, SDDS_BINARY, 1, NULL, NULL, filename) ||
        !SDDS_DefineSimpleParameter(SDDSout, "z", "m", SDDS_DOUBLE) ||
        !SDDS_DefineSimpleParameter(SDDSout, "CaseC", "#", SDDS_LONG) ||
        !SDDS_DefineSimpleParameter(SDDSout, "CaseD1", "#", SDDS_LONG) ||
        !SDDS_DefineSimpleParameter(SDDSout, "CaseD2", "#", SDDS_LONG) ||
        !SDDS_DefineSimpleParameter(SDDSout, "x", NULL, SDDS_DOUBLE) ||
        !SDDS_DefineSimpleParameter(SDDSout, "dsMax", NULL, SDDS_DOUBLE) ||
        !SDDS_DefineSimpleParameter(SDDSout, "phi0", NULL, SDDS_DOUBLE) ||
        !SDDS_DefineSimpleParameter(SDDSout, "phi1", NULL, SDDS_DOUBLE) ||
        !SDDS_DefineSimpleColumn(SDDSout, "s", "m", SDDS_DOUBLE) ||
        !SDDS_DefineSimpleColumn(SDDSout, "LinearDensity", "C/s", SDDS_DOUBLE) ||
        !SDDS_DefineSimpleColumn(SDDSout, "LinearDensityDeriv", "C/s$a2$n", SDDS_DOUBLE) ||
        !SDDS_DefineSimpleColumn(SDDSout, "DeltaGamma", NULL, SDDS_DOUBLE) ||
        !SDDS_DefineSimpleColumn(SDDSout, "GammaDeriv", "1/m", SDDS_DOUBLE) ||
        !SDDS_WriteLayout(SDDSout)) {
      SDDS_SetError("Problem setting up output file for CSRDRIFT (Stupakov mode)");
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
    }
    *active = 1;
  }
  for (i=0; i<nBins; i++) {
    ctHist[i] *= MPCharge*c_mks;
    ctHistDeriv[i] *= MPCharge*sqr(c_mks)/dct;
  }
  if (!SDDS_StartPage(SDDSout, nBins) ||
      !SDDS_SetColumn(SDDSout, SDDS_SET_BY_NAME, dGamma, nBins, "DeltaGamma")  ||
      !SDDS_SetColumn(SDDSout, SDDS_SET_BY_NAME, ctHist, nBins, "LinearDensity")  ||
      !SDDS_SetColumn(SDDSout, SDDS_SET_BY_NAME, ctHistDeriv, nBins, "LinearDensityDeriv")) {
    SDDS_SetError("Problem writing to output file for CSRDRIFT (Stupakov mode)");
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
  }
  /* use ctHist array for output of s and ctHistDeriv for dGamma/ds */
  for (i=0; i<nBins; i++) {
    ctHist[i] = dct*(i+0.5-nBins/2);
    ctHistDeriv[i] = dGamma[i]/dz;
  }
  if (!SDDS_SetColumn(SDDSout, SDDS_SET_BY_NAME, ctHist, nBins, "s") ||
      !SDDS_SetColumn(SDDSout, SDDS_SET_BY_NAME, ctHistDeriv, nBins, "GammaDeriv") ||
      !SDDS_SetParameters(SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, 
                          "z", zTravel, "CaseC", nCaseC,
                          "CaseD1", nCaseD1, "CaseD2", nCaseD2, 
                          "x", x, "dsMax", dsMax, "phi0", phi0, "phi1", phi1,
                          NULL) ||
      !SDDS_WritePage(SDDSout)) {
    SDDS_SetError("Problem writing to output file for CSRDRIFT (Stupakov mode)");
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
  }
    
}


void apply_edge_effects(
                        double *x, double *xp, double *y, double *yp, 
                        double rho, double n, double beta, double he, double psi, long which_edge
                        )
{
  double h, tan_beta, tan2_beta, sec_beta, sec2_beta, h2;
  double R21, R43;
  double T111, T133, T211, T441, T331, T221, T233, T243, T431, T432;
  double x0, xp0, y0, yp0;

  h = 1/rho;
  R21 = h*(tan_beta=tan(beta));
  R43 = -h*tan(beta-psi);

  h2 = sqr(h);
  T111 = which_edge*h/2*(tan2_beta=sqr(tan_beta));
  T133 = -which_edge*h/2*(sec2_beta=sqr(sec_beta=1./cos(beta)));
  T211 = which_edge==-1?
    -n*h2*tan_beta:
    -h2*(n+tan2_beta/2)*tan_beta;
  T441 = -(T331 = T221 = -which_edge*h*tan2_beta);
  T233 =  which_edge==-1?
    h2*(n+.5+tan2_beta)*tan_beta:
    h2*(n-tan2_beta/2)*tan_beta;
  T243 = which_edge*h*tan2_beta;
  T431 = h2*(2*n+(which_edge==1?sec2_beta:0))*tan_beta;
  T432 = which_edge*h*sec2_beta;
  if (he!=0) {
    double term;
    term = h/2*he*sec2_beta*sec_beta;
    T211 += term;
    T233 -= term;
    T431 -= 2*term;
  }

  x0 = *x;  xp0 = *xp;  y0 = *y;  yp0 = *yp;
  *x  = x0  + T111*sqr(x0) + T133*sqr(y0);
  *xp = xp0 + R21*x0 + T211*sqr(x0) + T221*x0*xp0 + T233*sqr(y0) + T243*y0*yp0;
  *y  = y0  + T331*x0*y0;
  *yp = yp0 + R43*y0 + T441*yp0*x0 + T431*x0*y0 + T432*xp0*y0;
}

/* this is used solely to convert coordinates inside the element for
 * the purpose of generating output.  It ignores misalignments.
 */

void convertFromCSBendCoords(double **part, long np, double rho0, 
			     double cos_ttilt, double sin_ttilt, 
			     long ctMode)
{
  long ip;
  double x, y, xp, yp, *coord;

  for (ip=0; ip<np; ip++) {
    coord = part[ip];

    XP /= (1+X/rho0);
    YP /= (1+X/rho0);
    
    x  =  X*cos_ttilt -  Y*sin_ttilt;
    y  =  X*sin_ttilt +  Y*cos_ttilt;
    xp = XP*cos_ttilt - YP*sin_ttilt;
    yp = XP*sin_ttilt + YP*cos_ttilt;

    X  = x;
    Y  = y;
    XP = xp;
    YP = yp;

    if (ctMode)
      coord[4] /= c_mks;
  }
}


/* this is used solely to undo the transformation done by 
 * convertFromCSBendCoords
 */

void convertToCSBendCoords(double **part, long np, double rho0, 
			     double cos_ttilt, double sin_ttilt, long ctMode)
{
  long ip;
  double x, y, xp, yp, *coord;

  for (ip=0; ip<np; ip++) {
    coord = part[ip];

    x  =   X*cos_ttilt +  Y*sin_ttilt;
    y  =  -X*sin_ttilt +  Y*cos_ttilt;
    xp =  XP*cos_ttilt + YP*sin_ttilt;
    yp = -XP*sin_ttilt + YP*cos_ttilt;

    X  = x;
    Y  = y;
    XP = xp;
    YP = yp;

    XP *= (1+X/rho0);
    YP *= (1+X/rho0);

    if (ctMode)
      coord[4] *= c_mks;
  }
}

#include "fftpackC.h"
long applyLowPassFilter(double *histogram, long bins, double dx, double start, double end)
{
  long i, i1, i2;
  double fraction, dfraction, sum;
  double *realimag, dfrequency, length;
  long frequencies;

  if (!(realimag = (double*)malloc(sizeof(*realimag)*(bins+2))))
    SDDS_Bomb("allocation failure");

  if (end<start)
    end = start;

  frequencies = bins/2 + 1;
  length = dx*(bins-1);
  dfrequency = 1.0/length;
  realFFT2(realimag, histogram, bins, 0);

  i1 = start*frequencies;
  if (i1<0) 
    i1=0;
  if (i1>frequencies-1)
    i1 = frequencies-1;
  
  i2 = end*frequencies;
  if (i2<0) 
    i2=0;
  if (i2>frequencies-1)
    i2 = frequencies-1;
  
  dfraction = i1==i2? 0 : 1./(i2-i1);
  fraction = 1;
  for (i=i1; i<=i2; i++) {
    realimag[2*i  ] *= fraction;
    realimag[2*i+1] *= fraction;
    if ((fraction -= dfraction)>1)
      fraction = 0;
  }
  for (; i<frequencies; i++) {
    realimag[2*i  ] = 0;
    realimag[2*i+1] = 0;
  }

  realFFT2(realimag, realimag, bins, INVERSE_FFT);

  /* copy data to input buffer.
   * normalize to keep the sum constant
   * don't allow negative values 
   */
  for (i=sum=0; i<bins; i++) {
    sum += histogram[i];
    histogram[i] = realimag[i];
  }
  free(realimag);
  return correctDistribution(histogram, bins, sum);
}

long correctDistribution(double *array, long npoints, double desiredSum)
{
  double sum, factor;
  long nz, i;
  for (i=nz=sum=0; i<npoints; i++) {
    if (array[i]<0) {
      nz ++;
      array[i] = 0;
    }
    sum += array[i];
  }
  if (!sum)
    return nz;
  factor = desiredSum/sum;
  for (i=0; i<npoints; i++)
    array[i] *= factor;
  return nz;
}


