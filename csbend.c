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

static double Fy_0, Fy_x, Fy_x2, Fy_x3, Fy_x4;
static double Fy_y2, Fy_x_y2, Fy_x2_y2;
static double Fy_y4;

static double Fx_y, Fx_x_y, Fx_x2_y, Fx_x3_y;
static double Fx_y3, Fx_x_y3;

static double rho0, rho_actual, rad_coef;

static long particle_lost;
static double s_lost;
extern unsigned long multipoleKicksDone ;

long binParticleCoordinate(double **hist, long *maxBins,
                           double *lower, double *upper, double *binSize, long *bins,
                           double expansionFactor,
                           double **particleCoord, long nParticles, long coordinateIndex);

long track_through_csbend(double **part, long n_part, CSBEND *csbend, double p_error, double Po, double **accepted,
                          double z_start)
{
  static double nh, betah2, gammah3, deltah4;
  static double h, h2, h3;
  long i_part, i_top;
  double rho, s, Fx, Fy;
  double x, xp, y, yp, dp, y2, dp0;
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
  double e1_kick_limit, e2_kick_limit;
  
  if (!csbend)
    bomb("null CSBEND pointer (track_through_csbend)", NULL);

  if (csbend->angle==0)
    bomb("angle = 0 for csbend", NULL);

  if (csbend->integration_order!=2 && csbend->integration_order!=4)
    bomb("CSBEND integration_order is invalid--must be either 2 or 4", NULL);

  if (csbend->use_bn) {
    rho0 = csbend->length/csbend->angle;
    csbend->k2 = csbend->b2/rho0;
    csbend->k3 = csbend->b3/rho0;
    csbend->k4 = csbend->b4/rho0;
  }
  
  if (csbend->angle<0) {
    angle = -csbend->angle;
    e1    = -csbend->e1;
    e2    = -csbend->e2;
    etilt = csbend->etilt;
    tilt  = csbend->tilt + PI;
    rho0  = -csbend->length/angle;
    n     = -sqr(rho0)*csbend->k1;
    beta  = 0.5*csbend->k2*pow3(rho0);
    gamma = csbend->k3*pow4(rho0)/6.;
    delta = csbend->k4*pow5(rho0)/24.;
    rho0  = -rho0;
  }
  else {
    angle = csbend->angle;
    e1    = csbend->e1;
    e2    = csbend->e2;
    etilt = csbend->etilt;
    tilt  = csbend->tilt;
    rho0  = csbend->length/angle;
    n     = -sqr(rho0)*csbend->k1;
    beta  = 0.5*csbend->k2*pow3(rho0);
    gamma = csbend->k3*pow4(rho0)/6.;
    delta = csbend->k4*pow5(rho0)/24.;
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

  dxf =  csbend->dx*cos(angle) + (-csbend->dz)*sin(angle);
  dzf = -csbend->dx*sin(angle) + (-csbend->dz)*cos(angle);
  dyf = csbend->dy;

  i_top = n_part-1;
  multipoleKicksDone += n_part*csbend->n_kicks*(csbend->integration_order==4?4:1);

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
      /* apply edge focusing */
      rho = (1+dp)*rho_actual;
      delta_xp = tan(e1)/rho*x;
      if (e1_kick_limit>0 && fabs(delta_xp)>e1_kick_limit)
        delta_xp = SIGN(delta_xp)*e1_kick_limit;
      xp += delta_xp;
      yp -= tan(e1-psi1/(1+dp))/rho*y;
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
      delta_xp = tan(e2)/rho*x;
      if (e2_kick_limit>0 && fabs(delta_xp)>e2_kick_limit)
        delta_xp = SIGN(delta_xp)*e2_kick_limit;
      xp += delta_xp;
      yp -= tan(e2-psi2/(1+dp))/rho*y;
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
    if (rad_coef) {
      denom = sqrt(sqr(1+DPoP)-sqr(QX)-sqr(QY));
      xp = QX/denom;
      yp = QY/denom;
      QX /= (1+DPoP);
      QY /= (1+DPoP);
      DPoP -= rad_coef*sqr(1+DPoP)*(sqr(Fx)+sqr(Fy))*sqrt(sqr(1+X/rho0)+sqr(xp)+sqr(yp))*ds;
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
    if (rad_coef) {
      QX /= (1+DPoP);
      QY /= (1+DPoP);
      DPoP -= rad_coef*sqr(1+DPoP)*(sqr(Fx)+sqr(Fy))*(1+X/rho0)*ds;
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
    if (rad_coef) {
      QX /= (1+DPoP);
      QY /= (1+DPoP);
      DPoP -= rad_coef*sqr(1+DPoP)*(sqr(Fx)+sqr(Fy))*(1+X/rho0)*ds;
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
    if (rad_coef) {
      QX /= (1+DPoP);
      QY /= (1+DPoP);
      DPoP -= rad_coef*sqr(1+DPoP)*(sqr(Fx)+sqr(Fy))*(1+X/rho0)*ds;
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
  long bins, new;
  double dctBin, s0, ds0, zLast, z0;
  double thetaRad;
  double S11, S12, S22;
  double *dGamma;
} CSR_LAST_WAKE;
CSR_LAST_WAKE csrWake = {
  0, 0, 
  0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0,
  NULL};

long track_through_csbendCSR(double **part, long n_part, CSRCSBEND *csbend, double p_error, 
                             double Po, double **accepted, double z_start, double z_end,
                             CHARGE *charge, char *rootname)
{
  static double nh, betah2, gammah3, deltah4;
  static double h, h2, h3;
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
  double rho, Fx, Fy, y2;
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
  double wavelength, criticalWavelength;
  
  if (!csbend)
    bomb("null CSBEND pointer (track_through_csbend)", NULL);

  if (csbend->angle==0)
    bomb("angle = 0 for csbend", NULL);

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
  if (csbend->SGDerivOrder<=0)
    csbend->SGDerivOrder = csbend->SGOrder;
  
  if (n_part>maxParticles &&
      (!(beta0=SDDS_Realloc(beta0, sizeof(*beta0)*(maxParticles=n_part))) ||
       !(particleLost=SDDS_Realloc(particleLost, sizeof(*particleLost)*n_part))))
    bomb("Memory allocation failure (track_through_csbendCSR)", NULL);

  if (csbend->use_bn) {
    rho0 = csbend->length/csbend->angle;
    csbend->k2 = csbend->b2/rho0;
    csbend->k3 = csbend->b3/rho0;
    csbend->k4 = csbend->b4/rho0;
  }

  if (csbend->angle<0) {
    angle = -csbend->angle;
    e1    = -csbend->e1;
    e2    = -csbend->e2;
    etilt = csbend->etilt;
    tilt  = csbend->tilt + PI;
    rho0  = -csbend->length/angle;
    n     = -sqr(rho0)*csbend->k1;
    beta  = 0.5*csbend->k2*pow3(rho0);
    gamma = csbend->k3*pow4(rho0)/6.;
    delta = csbend->k4*pow5(rho0)/24.;
    rho0  = -rho0;
  }
  else {
    angle = csbend->angle;
    e1    = csbend->e1;
    e2    = csbend->e2;
    etilt = csbend->etilt;
    tilt  = csbend->tilt;
    rho0  = csbend->length/angle;
    n     = -sqr(rho0)*csbend->k1;
    beta  = 0.5*csbend->k2*pow3(rho0);
    gamma = csbend->k3*pow4(rho0)/6.;
    delta = csbend->k4*pow5(rho0)/24.;
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

  dxf =  csbend->dx*cos(angle) + (-csbend->dz)*sin(angle);
  dzf = -csbend->dx*sin(angle) + (-csbend->dz)*cos(angle);
  dyf = csbend->dy;

  if (csbend->histogramFile && !csbend->fileActive) {
    /* set up SDDS output file for CSR monitoring */
    csbend->fileActive = 1;
    csbend->histogramFile = compose_filename(csbend->histogramFile, rootname);
    if (!SDDS_InitializeOutput(&csbend->SDDSout, SDDS_BINARY, 1, NULL, NULL, csbend->histogramFile) ||
        !SDDS_DefineSimpleParameter(&csbend->SDDSout, "Pass", NULL, SDDS_LONG) ||
        !SDDS_DefineSimpleParameter(&csbend->SDDSout, "Kick", NULL, SDDS_LONG) ||
        !SDDS_DefineSimpleParameter(&csbend->SDDSout, "pCentral", "m$be$nc", SDDS_DOUBLE) ||
        !SDDS_DefineSimpleParameter(&csbend->SDDSout, "Angle", NULL, SDDS_DOUBLE) ||
        !SDDS_DefineSimpleParameter(&csbend->SDDSout, "SlippageLength", "m", SDDS_DOUBLE) ||
        !SDDS_DefineSimpleParameter(&csbend->SDDSout, "TotalBunchLength", "m", SDDS_DOUBLE) ||
        !SDDS_DefineSimpleParameter(&csbend->SDDSout, "BinSize", "m", SDDS_DOUBLE) ||
        !SDDS_DefineSimpleColumn(&csbend->SDDSout, "s", "m", SDDS_DOUBLE) ||
        !SDDS_DefineSimpleColumn(&csbend->SDDSout, "LinearDensity", "C/s", SDDS_DOUBLE) ||
        !SDDS_DefineSimpleColumn(&csbend->SDDSout, "LinearDensityDeriv", "C/s$a2$n", SDDS_DOUBLE) ||
        !SDDS_DefineSimpleColumn(&csbend->SDDSout, "DeltaGamma", NULL, SDDS_DOUBLE) ||
        !SDDS_DefineSimpleColumn(&csbend->SDDSout, "DeltaGammaT1", NULL, SDDS_DOUBLE) ||
        !SDDS_DefineSimpleColumn(&csbend->SDDSout, "DeltaGammaT2", NULL, SDDS_DOUBLE) ||
        !SDDS_WriteLayout(&csbend->SDDSout)) {
      SDDS_SetError("Problem setting up output file for CSR");
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
  
  multipoleKicksDone += n_part*csbend->n_kicks*(csbend->integration_order==4?4:1);
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
      delta_xp = tan(e1)/rho*X;
      XP += delta_xp;
      YP -= tan(e1-psi1/(1+DP))/rho*Y;
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

  CSRConstant = 2*macroParticleCharge*e_mks/pow(3*rho0*rho0, 1./3.)/(4*PI*epsilon_o*me_mks*sqr(c_mks));
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
    
    /* compute CSR potential function */
    /* - first make a density histogram */
    ctLower = ctUpper = dct = 0;
    if ((nBinned = 
         binParticleCoordinate(&ctHist, &maxBins,
                               &ctLower, &ctUpper, &dct, &nBins, 
                               1.2, part, n_part, 4))!=n_part) {
      fprintf(stdout, "Only %ld of %ld particles binned for CSR\n", nBinned, n_part);
      fflush(stdout);
    }
    
    /* - smooth the histogram, normalize to get linear density, and 
       copy in preparation for taking derivative
       */
    SavitzyGolaySmooth(ctHist, nBins, csbend->SGOrder, csbend->SGHalfWidth, csbend->SGHalfWidth,  0);
    for (iBin=0; iBin<nBins; iBin++) {
      denom[iBin] = pow(dct*iBin, 1./3.);
      ctHistDeriv[iBin] = (ctHist[iBin] /= dct);
    }
    /* - compute derivative with smoothing.  The deriv is w.r.t. index number and
     * I won't scale it now as it will just fall out in the integral 
     */
    SavitzyGolaySmooth(ctHistDeriv, nBins, csbend->SGDerivOrder, 
                       csbend->SGDerivHalfWidth, csbend->SGDerivHalfWidth, 1);

    phiBend += angle/csbend->n_kicks;
    slippageLength = rho0*ipow(phiBend, 3)/24.0;
    slippageLength13 = pow(slippageLength, 1./3.);
    diSlippage = slippageLength/dct;
    diSlippage4 = 4*slippageLength/dct;
    for (iBin=0; iBin<nBins; iBin++) {
      T1[iBin] = T2[iBin] = 0;
      if (CSRConstant) {
        if (csbend->steadyState) {
          for (iBinBehind=iBin+1; iBinBehind<nBins; iBinBehind++)
            T1[iBin] += ctHistDeriv[iBinBehind]/denom[iBinBehind-iBin];
        } else {
          for (iBinBehind=iBin+1; iBinBehind<=(iBin+diSlippage) && iBinBehind<nBins; iBinBehind++)
            T1[iBin] += ctHistDeriv[iBinBehind]/denom[iBinBehind-iBin];
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
    
    for (i_part=0; i_part<n_part; i_part++) {
      coord = part[i_part];
      if (!particleLost[i_part]) {
        /* apply CSR kick */
        iBin = (CT-ctLower)/dct;
        if (iBin>=0 && iBin<nBins)
          DP += dGamma[iBin]/Po;
      }
    }

    if (csbend->fileActive && kick%csbend->outputInterval==0) {
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
                              "Pass", -1, "Kick", kick, 
                              "pCentral", Po, "Angle", phiBend, "SlippageLength", slippageLength,
                              "TotalBunchLength", ctUpper-ctLower,
                              "BinSize", dct, NULL))
        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
      /* use T1 array to output s */
      for (iBin=0; iBin<nBins; iBin++)
        T1[iBin] = ctLower-(ctLower+ctUpper)/2.0+dct*(iBin+0.5);
      if (!SDDS_SetColumn(&csbend->SDDSout, SDDS_SET_BY_NAME, T1, nBins, "s") ||
          !SDDS_WritePage(&csbend->SDDSout))
        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
    }
  }
  
  /* prepare soem data for use by CSRDRIFT element */
  csrWake.dctBin = dct;
  ctLower = ctUpper = dct = 0;
  if ((nBinned = 
       binParticleCoordinate(&ctHist, &maxBins,
                             &ctLower, &ctUpper, &dct, &nBins, 
                             1.2, part, n_part, 4))!=n_part) {
    fprintf(stdout, "Only %ld of %ld particles binned for CSR\n", nBinned, n_part);
    fflush(stdout);
  }
  csrWake.s0 = ctLower;
  
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

    if (particleLost[i_part]) {
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
      delta_xp = tan(e2)/rho*X;
      XP += delta_xp;
      YP -= tan(e2-psi2/(1+DP))/rho*Y;
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

  /* prepare more data for CSRDRIFT */
  rms_emittance(part, 0, 1, i_top+1, &csrWake.S11, &csrWake.S12, &csrWake.S22);

  /* compute angular spread of radiation, using Wiedemann's formula */
  /* wavelength ~ (sigma z)/2 */
  wavelength = beam_width(0.6826, part, i_top+1, 4)/2;
  criticalWavelength = 4.19/ipow(Po, 3)*rho_actual;
  csrWake.thetaRad = 0.5463e-3/(Po*0.511e-3)/pow(criticalWavelength/wavelength, 1./3.);

#ifdef DEBUG
  fprintf(stdout, "wavelength = %le, critWL = %le, thetaRad = %le\n",
          wavelength, criticalWavelength, csrWake.thetaRad);
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
    *upper = -(*lower = DBL_MAX);
    for (iParticle=0; iParticle<nParticles; iParticle++) {
      value = particleCoord[iParticle][coordinateIndex];
      if (value<*lower)
        *lower = value;
      if (value>*upper)
        *upper = value;
    }
    if (expansionFactor!=1) {
      double center, range;
      center = (*lower+*upper)/2;
      range = (*upper-*lower)*expansionFactor;
      *lower = center-range/2;
      *upper = center+range/2;
    }
  }
  
  if (*binSize>0)
    *bins = (*upper-*lower)/(*binSize);
  *binSize = (*upper-*lower)/(*bins);
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


long track_through_driftCSR(double **part, long np, CSRDRIFT *csrDrift, 
                            double Po, double **accepted, double zStart)
{
  long iPart, iKick, iBin, binned;
  double *coord, t, p, beta, dz, ct0, factor, dz0, dzFirst;
  double ctmin, ctmax, spreadFactor;
  double zTravel;
  
  dz = (dz0=csrDrift->length/csrDrift->nKicks)/2;
  dzFirst = zStart - csrWake.zLast;
#ifdef DEBUG
  fprintf(stdout, "dzFirst = %21.15e\n", dzFirst);
  fprintf(stdout, "ct0 = %21.15e\n", csrWake.s0+dz+dzFirst);
#endif

  ctmin = DBL_MAX;
  ctmax = -DBL_MAX;

  
  zTravel = zStart-csrWake.z0;  /* total distance traveled by radiation to reach this point */
  for (iKick=0; iKick<csrDrift->nKicks; iKick++) {
    /* first drift is dz=dz0/2, others are dz0 */
    if (iKick==1)
      dz = dz0;
    zTravel += dz;
    
    /* propagate particles forward, converting s to c*t=s/beta */
    for (iPart=0; iPart<np; iPart++) {
      coord = part[iPart];
      coord[0] += coord[1]*dz;
      coord[2] += coord[3]*dz;
      p = Po*(1+coord[5]);
      beta = p/sqrt(p*p+1);
      coord[4] = (coord[4]+dz*sqrt(1+sqr(coord[1])+sqr(coord[3])))/beta;
      if (iKick==0) {
        if (coord[4]>ctmax)
          ctmax = coord[4];
        if (coord[4]<ctmin)
          ctmin = coord[4];
      }
    }
#ifdef DEBUG
    if (iKick==0) {
      fprintf(stdout, "ctmin, max = %21.15e, %21.15e\n",
              ctmin, ctmax);
      fflush(stdout);
    }
#endif

    factor = 1;
    if (csrWake.dGamma) {
      /* propagate wake forward */
      csrWake.s0 += dz+dzFirst;   /* accumulates position of back end of the radiation pulse */
      ct0 = csrWake.s0;
      
      if (csrDrift->attenuationLength>0) {
        /* attenuate wake */
        if ((factor = exp(-(dz+dzFirst)/csrDrift->attenuationLength))<1) {
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
    if (csrDrift->spread) {
      /* compute loss of on-axis field due to spread of beam using a simple-minded
       * computation of beam sizes */
      factor *= (spreadFactor =
                 sqrt(csrWake.S11/(csrWake.S11 + 
                                   2*zTravel*csrWake.S12 + 
                                   zTravel*zTravel*(sqr(csrWake.thetaRad)+csrWake.S22))));
#ifdef DEBUG
      fprintf(stdout, "Spread factor is %le\n", spreadFactor);
      fprintf(stdout, "S11 = %le  S12 = %le  S22 = %le\nthetaRad = %le  z = %le\n",
              csrWake.S11, csrWake.S12, csrWake.S22, csrWake.thetaRad, zTravel);
#endif
    }
    
    dzFirst = 0;

    /* apply kick to each particle and convert back to normal coordinates */
    for (iPart=binned=0; iPart<np; iPart++) {
      coord = part[iPart];
      if (csrWake.dGamma) {
        iBin = (coord[4]-ct0)/csrWake.dctBin;
        if (iBin>=0 && iBin<csrWake.bins) {
          coord[5] += csrWake.dGamma[iBin]/Po*factor;
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

  csrWake.zLast += csrDrift->length;

  if (csrWake.dGamma) {
    /* propagate wake forward */
    csrWake.s0 += dz;
    ct0 = csrWake.s0;
    
    if (csrDrift->attenuationLength>0) {
      /* attenuate wake */
      if ((factor = exp(-dz/csrDrift->attenuationLength))<1) {
        for (iBin=0; iBin<csrWake.bins; iBin++)
            csrWake.dGamma[iBin] *= factor;
      }
    }
  }
  return np;
}
