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

    log_entry("track_through_csbend");

    if (!csbend)
        bomb("null CSBEND pointer (track_through_csbend)", NULL);

    if (csbend->angle==0)
        bomb("angle = 0 for csbend", NULL);

    if (csbend->integration_order!=2 && csbend->integration_order!=4)
        bomb("CSBEND integration_order is invalid--must be either 2 or 4", NULL);

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
    rho_actual = 1/((1+fse)*h);

    /* angles for fringe-field effects */
    Kg   = 2*csbend->hgap*csbend->fint;
    psi1 = Kg/rho0/cos(e1)*(1+sqr(sin(e1)));
    psi2 = Kg/rho0/cos(e2)*(1+sqr(sin(e2)));

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

    log_entry("track_through_csbend.1");

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
        printf("pre-tilt offsets due to ETILT=%le:  %le %le %le %le %le\n",
            etilt, dcoord_etilt[0], dcoord_etilt[1], dcoord_etilt[2],
             dcoord_etilt[3], dcoord_etilt[4]);
#endif

        /* rotate by tilt to get into same frame as bend equations. */
        rotate_coordinates(dcoord_etilt, tilt);
#ifdef DEBUG
        printf("offsets due to ETILT=%le:  %le %le %le %le %le\n",
            etilt, dcoord_etilt[0], dcoord_etilt[1], dcoord_etilt[2],
             dcoord_etilt[3], dcoord_etilt[4]);
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

    log_exit("track_through_csbend.1");

    log_entry("track_through_csbend.2");
    i_top = n_part-1;

    for (i_part=0; i_part<=i_top; i_part++) {
        if (!part) {
            fprintf(stderr, "error: null particle array found (working on particle %ld) (track_through_csbend)\n", i_part);
            abort();
            }
        if (!(coord = part[i_part])) {
            fprintf(stderr, "error: null coordinate pointer for particle %ld (track_through_csbend)\n", i_part);
            abort();
            }
        if (accepted && !accepted[i_part]) {
            fprintf(stderr, "error: null accepted particle pointer for particle %ld (track_through_csbend)\n", i_part);
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
            rho = (1+dp)*rho0;
            xp += tan(e1)/rho*x;
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
                fprintf(stderr, "error: couldn't swap particles %ld and %ld--latter is null pointer (track_through_csbend)\n",
                    i_part, i_top);
                abort();
                }
            SWAP_PTR(part[i_part], part[i_top]);
            if (accepted) {
                if (!accepted[i_top]) {
                    fprintf(stderr, 
                        "error: couldn't swap acceptance data for particles %ld and %ld--latter is null pointer (track_through_csbend)\n",
                        i_part, i_top);
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
            rho = (1+dp)*rho0;
            xp += tan(e2)/rho*x;
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
    log_exit("track_through_csbend.2");

    log_exit("track_through_csbend");
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

