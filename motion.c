/* Copyright 1994 by Michael Borland and Argonne National Laboratory,
 * all rights reserved.
 */
/* routine: motion()
 * purpose: integrate equations of motion for a particle in a 
 *          static or time-dependent field.  This module is appropriate
 *          for elements with no coordinate system curvature.  lorentzX.c
 *          is better if there is coordinate system curvature.
 *
 * The equations of motion are
 *
 *      '         _                                                  _
 *     P  =  -eE (X, tau)/(m.c.omega) - e/(m.gamma.omega).eps   P B (X, tau)
 *      i       i                                            ijk j k
 *
 *      '
 *     X  =  P /gamma
 *      i     i
 *
 *       _   ____        _     _ 
 * where P = beta.gamma, X = k.x, and derivatives are with respect to
 * tau=omega.time.   E has units of V/m and B has units of T.
 *
 * For resonant fields, omega is the angular frequency and k=omega/c.
 * For ramped fields, omega is 1/ramp_time and k=omega/c.
 * For static fields, I arbitrarily take omega=2*PI c/L, where L is the total length.
 *
 * Michael Borland, 1989
 */
#include "mdb.h"
#include "match_string.h"
#include "track.h"
#include "matlib.h"

void (*set_up_derivatives(void *field, long field_type, double *kscale, 
        double z_start, double *Z_end, double *tau_start, double **part, long n_part, 
        double Po, double *X_limit, double *Y_limit, double *accuracy,
        long *n_steps))(double *, double *, double);
void derivatives_tmcf_mode(
    double *qp,       /* derivatives w.r.t. tau */
    double *q,        /* X,Y,Z,Px,Py,Pz */
    double tau     
    );
void derivatives_rftmEz0(
    double *qp,       /* derivatives w.r.t. tau */
    double *q,        /* X,Y,Z,Px,Py,Pz */
    double tau     
    );
void derivatives_ce_plates(
    double *qp,       /* derivatives w.r.t. tau */
    double *q,        /* X,Y,Z,Px,Py,Pz */
    double tau 
    );
void derivatives_tw_plates(
    double *qp,       /* derivatives w.r.t. tau */
    double *q,        /* X,Y,Z,Px,Py,Pz */
    double tau 
    );
void derivatives_tw_linac(
    double *qp,       /* derivatives w.r.t. tau */
    double *q,        /* X,Y,Z,Px,Py,Pz */
    double tau 
    );
void derivatives_twmta(
    double *qp,       /* derivatives w.r.t. tau */
    double *q,        /* X,Y,Z,Px,Py,Pz */
    double tau 
    );
void derivatives_mapSolenoid(
    double *qp,       /* derivatives w.r.t. tau */
    double *q,        /* X,Y,Z,Px,Py,Pz */
    double tau     
    );

void setupRftmEz0FromFile(RFTMEZ0 *rftmEz0, double frequency, double length);
void setupRftmEz0SolenoidFromFile(RFTMEZ0 *rftmEz0, double length, double k);
void setupMapSolenoidFromFile(MAP_SOLENOID *mapSol, double length);
double exit_function(double *qp, double *q, double phase);
double *select_fiducial(double **part, long n_part, char *mode);
void select_integrator(char *desired_method);
void input_impulse_tw_linac(double *P, double *q);
void output_impulse_tw_linac(double *P, double *q);
void input_impulse_twmta(double *P, double *q);
void output_impulse_twmta(double *P, double *q);

static void (*derivatives)(double *xqp, double *xq, double xt);
static void (*input_impulse)(double *P, double *q);
static void (*output_impulse)(double *P, double *q);
static long (*integrator)();

long analyzeSpacing(double *z, long nz, double *dzReturn, FILE *fpError);

#define OUT_OF_APERTURE (-1)

/* end of element, in scaled coordinates */
double Z_end;

/* offsets to be added to particle position to get coordinates relative
 * to element transverse center.
 */
double X_offset, Y_offset;

/* position of centers of X and Y apertures relative to straight-ahead
 * centerline
 */
double X_aperture_center, Y_aperture_center;

/* maximum coordinates relative to aperture center
 * that a particle can have without being considered lost
 */
double X_limit, Y_limit;

/* coordinates reached by particle, relative to
 * aperture center of element, when lost.
 */
double X_out, Y_out, Z_out;
long limit_hit, radial_limit, derivCalls=0;

/* flag to indicate that element changes the central momentum */
static long change_p0;

/* rotation matrices for tilt (roll), yaw, pitch */
static MATRIX *partRot, *fieldRot;
/* geometrical center of the element used for rotation */
double X_center, Y_center, Z_center;

#if defined(DEBUG)
static long starting_integration;
static double initial_phase=0, final_phase=0;
#endif

long motion(
    double **part,
    long n_part,
    void *field,
    long field_type,
    double *pCentral,
    double *dgamma,
    double *dP,
    double **accepted,
    double z_start
    )
{
    double tau, tau_limit, P_central;
    double tau_start;
    long n_eq = 6;
    static double q[6], dqdt[6];         /* X,Y,Z,Px,Py,Pz */
    static long misses[6], accmode[6];
    static double accuracy[6], tiny[6];
    double tolerance, hmax, hrec;
    double *coord;
    long i_part, i_top, rk_return;
    double kscale, *P, Po, gamma;
    static double gamma0, P0[3];
    long n_steps;
    double tol_factor, end_factor;
    double X, Y;

    if (!n_part)
        return(0);
    change_p0 = 0;
    P_central = *pCentral;
    
    log_entry("motion");

    X_aperture_center = Y_aperture_center = 0;
    X_limit = Y_limit = 0;
    X_offset = Y_offset = 0;
    radial_limit = 0;
    input_impulse = output_impulse = NULL;
    derivatives = set_up_derivatives(field, field_type, &kscale, z_start, &Z_end, 
                                     &tau_start, part, n_part, P_central, 
                                     &X_limit, &Y_limit, &tolerance, &n_steps);
    if (n_steps<2)
        bomb("too few steps for integration", NULL);

    P = q+3;
    i_top = n_part-1;
    *dgamma = dP[0] = dP[1] = dP[2] = 0;

    fill_long_array(accmode, 6, 3);
    fill_double_array(tiny , 6, 1e-16);
    hrec = 0;

    for (i_part=0; i_part<=i_top; i_part++) {
        tol_factor = 1;
        end_factor = 1.5;
        coord = part[i_part];
        q[0] = coord[0]*kscale;
        q[1] = coord[2]*kscale;
        if ((X_limit && fabs(q[0]-X_aperture_center)>X_limit) || isnan(q[0]) || isinf(q[0]) ||
            (Y_limit && fabs(q[1]-Y_aperture_center)>Y_limit) || isnan(q[1]) || isinf(q[1])) {
            SWAP_PTR(part[i_part], part[i_top]);
            part[i_top][4] = z_start;   /* record position of loss */
            part[i_top][5] = P_central*(1+part[i_top][5]);
            if (accepted)
                SWAP_PTR(accepted[i_part], accepted[i_top]);
            i_top--;
            i_part--;
            }
        else {
            fill_double_array(accuracy, 3, tolerance*kscale);
            fill_double_array(accuracy+3, 3, tolerance);
            fill_long_array(misses, 6, 0);
            q[2] = 0.0;
            Po   = P_central*(1+coord[5]);
            gamma0 = gamma = sqrt(1+sqr(Po));
            P0[2] = P[2] = Po/sqrt(1+sqr(coord[1])+sqr(coord[3]));
            P0[0] = P[0] = coord[1]*P[2];
            P0[1] = P[1] = coord[3]*P[2];            
            X_out = Y_out = Z_out = 0;
            tau_limit = end_factor*Z_end*gamma/P[2] +
                      (tau = tau_start + coord[4]*kscale*gamma/Po);
            hmax = (tau_limit - tau)/n_steps;
            if (integrator==rk_odeint3_na || integrator==mmid_odeint3_na)
                hrec = (tau_limit-tau)/n_steps;
            if (hrec<=0)
                hrec = hmax;
            limit_hit = 0;
#if defined(DEBUG)
            starting_integration = 1;
            initial_phase = final_phase = 0;
#endif
            if (input_impulse)
                input_impulse(P, q);
            derivCalls = 0;
            switch (rk_return = (*integrator)(q, derivatives, n_eq, accuracy, accmode, tiny, misses,
                  &tau, tau_limit, sqr(tolerance)*(tau_limit-tau), hrec, hmax, &hrec, exit_function, 
                  sqr(tolerance)*kscale)) {
                case DIFFEQ_ZERO_STEPSIZE:
                case DIFFEQ_CANT_TAKE_STEP:
                case DIFFEQ_OUTSIDE_INTERVAL:
                case DIFFEQ_XI_GT_XF:
                case DIFFEQ_EXIT_COND_FAILED:
                    fprintf(stdout, "Integration failure: %s\n", diffeq_result_description(rk_return));
                    fflush(stdout);
/*                    for (i=0; i<6; i++) 
                        fprintf(stdout, "%11.4e  ", coord[i]);
                        fflush(stdout);
                    fprintf(stdout, "\ngamma = %.4e,  P=(%.4e, %.4e, %.4e)\n",
                        gamma0, P0[0], P0[1], P0[2]);
                    fflush(stdout);
                    fprintf(stdout, "tolerance = %e     end_factor = %e\n",
                        tolerance, end_factor);
                    fflush(stdout);
*/
                    SWAP_PTR(part[i_part], part[i_top]);
                    part[i_top][4] = z_start;
                    part[i_top][5] = sqrt(sqr(P_central*(1+part[i_top][5]))+1);
                    if (accepted)
                        SWAP_PTR(accepted[i_part], accepted[i_top]);
                    i_top--;
                    i_part--;
                    break;
                default:
                    (*derivatives)(dqdt, q, tau);
                    Po      = sqrt(sqr(P[0])+sqr(P[1])+sqr(P[2]));
                    gamma   = sqrt(1+sqr(Po));
                    *dgamma += gamma-gamma0;
                    dP[0]   += P[0] - P0[0];
                    dP[1]   += P[1] - P0[1];
                    dP[2]   += P[2] - P0[2];
#if defined(DEBUG)
                    fprintf(stdout, "initial, final phase: %21.15e, %21.15e\n", initial_phase, final_phase);
                    fflush(stdout);
                    fprintf(stdout, "deriv calls: %ld\n", derivCalls);
                    fflush(stdout);
                    fprintf(stdout, "Exit: tau=%le, ef=%le\n", tau, exit_function(dqdt, q, tau));
                    fflush(stdout);
                    fprintf(stdout, "Pout = %21.15e  Zout = %21.15le\n", 
                            Po, q[2]/kscale);
                    fflush(stdout);
#endif
                    if (!limit_hit) {
                        if (isnan(q[0]) || isinf(q[0]) || isnan(q[1]) || isinf(q[1])) {
                            X_out = Y_out = DBL_MAX;
                            Z_out = q[2];
                            limit_hit = 1;
                            }
                        else if (radial_limit) {
                            if (X_limit) {
                                X = fabs(q[0]-X_aperture_center);
                                Y = fabs(q[1]-Y_aperture_center);
                                if (sqr(X_limit)<(sqr(X)+sqr(Y))) {
                                    X_out = X;
                                    Y_out = Y;
                                    Z_out = q[2];
                                    limit_hit = 1;
                                    }
                                }
                            }
                        else {
                            if ((X_limit && (X=fabs(q[0]-X_aperture_center))>X_limit) ||
                                (Y_limit && (Y=fabs(q[1]-Y_aperture_center))>Y_limit) ) {
                                X_out = X;
                                Y_out = Y;
                                Z_out = q[2];
                                limit_hit = 1;
                                }
                            }
                        }
                    if (limit_hit || P[2]<=0) {
                         SWAP_PTR(part[i_part], part[i_top]);
                         /* record position of loss */
                         part[i_top][0] = (X_out + X_aperture_center)/kscale;
                         part[i_top][2] = (Y_out + Y_aperture_center)/kscale;
                         part[i_top][4] = Z_out/kscale + z_start;
                         part[i_top][5] = Po;
                         if (accepted)
                             SWAP_PTR(accepted[i_part], accepted[i_top]);
                         i_top--;
                         i_part--;
                         continue;
                         }
                    if (output_impulse)
                        output_impulse(P, q);
                    coord[0] = q[0]/kscale;
                    coord[2] = q[1]/kscale;
                    coord[1] = P[0]/P[2];
                    coord[3] = P[1]/P[2];
                    coord[5] = (Po-P_central)/P_central;
                    coord[4] = (tau-tau_start)*Po/(kscale*gamma);
                    break;
                }
            }
        }
    if (change_p0)
      do_match_energy(part, n_part, pCentral, 0);
    
    log_exit("motion");
    return(i_top+1);
    }
 
 
void *field_global;
 
void (*set_up_derivatives(
                          void *field,
                          long field_type,
                          double *kscale,
                          double z_start,
                          double *Z_end,
                          double *tau_start,
                          double **part,
                          long n_part,
                          double P_central,
                          double *X_limit,
                          double *Y_limit,
                          double *accuracy,
                          long *n_steps
                          ))(double *, double *, double)
{

  TMCF_MODE *tmcf;
  CE_PLATES *cep;
  TW_PLATES *twp;
  TW_LINAC *twla;
  TWMTA *twmta;
  RFTMEZ0 *rftmEz0;
  MAP_SOLENOID *mapSol;
  
  double *fiducial, gamma, Po, Pz, gamma_w, omega;
  double Escale, Bscale;

  field_global = field;

  Escale = (Bscale = e_mks/me_mks)/c_mks;
  change_p0 = 0;
  
  switch (field_type) {
  case T_MAPSOLENOID:
    mapSol = field;
    *kscale = 1;
    if (!mapSol->initialized)
      setupMapSolenoidFromFile(mapSol, mapSol->length);
    *tau_start = 0;
    *Z_end = mapSol->length;
    X_offset = -(X_center = X_aperture_center = mapSol->dx);
    Y_offset = -(Y_center = Y_aperture_center = mapSol->dy);
    Z_center = mapSol->length/2;
    /* assume 1000m aperture for now */     
    *X_limit = *Y_limit = 1e3;
    *accuracy = mapSol->accuracy;
    *n_steps = mapSol->n_steps;
    mapSol->zMap0 = mapSol->length/2-mapSol->lUniform/2;
    mapSol->zMap1 = mapSol->length/2+mapSol->lUniform/2;
    mapSol->BxUniformScaled = mapSol->BxUniform*e_mks/(me_mks*c_mks);
    mapSol->ByUniformScaled = mapSol->ByUniform*e_mks/(me_mks*c_mks);
    select_integrator(mapSol->method);
    setupRotate3Matrix((void**)&partRot, -mapSol->eTilt, -mapSol->eYaw, -mapSol->ePitch);
    setupRotate3Matrix((void**)&fieldRot, mapSol->eTilt, mapSol->eYaw, mapSol->ePitch);
    return(derivatives_mapSolenoid);
    break;
  case T_RFTMEZ0:
    rftmEz0 = field;
    *kscale = (omega=PIx2*rftmEz0->frequency)/c_mks;
    if (rftmEz0->change_p0)
      change_p0 = 1;
    if (!rftmEz0->Ez) {
      setupRftmEz0FromFile(rftmEz0, rftmEz0->frequency, rftmEz0->length);
      setupRftmEz0SolenoidFromFile(rftmEz0, rftmEz0->length, *kscale);
      rftmEz0->initialized = 1;
    }
    if (!rftmEz0->fiducial_part) {
      /* This is the fiducial particle--the phase offset is set so 
       * that its phase when it reaches the cavity center will be 
       * approximately rftmEz0->phase + omega*rftmEz0->time_offset.
       * Several other quantities are calculated here also.
       */
      rftmEz0->k = *kscale;
      rftmEz0->fiducial_part = fiducial = select_fiducial(part, n_part, rftmEz0->fiducial);
      if (fiducial) {
        Po = P_central*(1+fiducial[5]);
        gamma = sqrt(1+sqr(Po));
        Pz = Po*sqrt(1-sqr(fiducial[0])-sqr(fiducial[2]));
        rftmEz0->phase0 =
          get_reference_phase((long)(rftmEz0->phase_reference+.5),
                              -rftmEz0->k*gamma*(rftmEz0->length/(2*Pz) + fiducial[4]/Po));
      }
      else {
        /* phase to v=c */
        rftmEz0->phase0 = get_reference_phase((long)(rftmEz0->phase_reference+.5),
                                           -rftmEz0->k*rftmEz0->length/2 + z_start);
        rftmEz0->fiducial_part = part[0];
      }
    }
    *tau_start = rftmEz0->phase + PIx2*rftmEz0->frequency*rftmEz0->time_offset +
      rftmEz0->phase0;
    *Z_end = rftmEz0->length*(*kscale);
    X_offset = -(X_center = X_aperture_center = rftmEz0->k*rftmEz0->dx);
    Y_offset = -(Y_center = Y_aperture_center = rftmEz0->k*rftmEz0->dy);
    Z_center = *Z_end/2;
    /* assume 1000m aperture for now */     
    *X_limit = rftmEz0->k*1e3;
    *Y_limit = rftmEz0->k*1e3;
    *accuracy = rftmEz0->accuracy;
    *n_steps = rftmEz0->n_steps;
    select_integrator(rftmEz0->method);
    setupRotate3Matrix((void**)&partRot, rftmEz0->eTilt, rftmEz0->eYaw, rftmEz0->ePitch);
    setupRotate3Matrix((void**)&fieldRot, rftmEz0->eTilt, rftmEz0->eYaw, rftmEz0->ePitch);
    return(derivatives_rftmEz0);
    break;
  case T_TMCF:
    tmcf = field;
    *kscale = (omega=PIx2*tmcf->frequency)/c_mks;
    if (!tmcf->fiducial_part) {
      /* This is the fiducial particle--the phase offset is set so 
       * that its phase when it reaches the cavity center will be 
       * approximately tmcf->phase + omega*tmcf->time_offset.
       * Several other quantities are calculated here also.
       */
      tmcf->k = (omega=PIx2*tmcf->frequency)/c_mks;
      *kscale = tmcf->k;
      tmcf->fiducial_part = fiducial = select_fiducial(part, n_part, tmcf->fiducial);
      if (fiducial) {
        Po = P_central*(1+fiducial[5]);
        gamma = sqrt(1+sqr(Po));
        Pz = Po*sqrt(1-sqr(fiducial[0])-sqr(fiducial[2]));
        tmcf->phase0 =
          get_reference_phase((long)(tmcf->phase_reference+.5),
                              -tmcf->k*gamma*(tmcf->length/(2*Pz) + fiducial[4]/Po));
      }
      else {
        /* phase to v=c */
        tmcf->phase0 = get_reference_phase((long)(tmcf->phase_reference+.5),
                                           -tmcf->k*tmcf->length/2 + z_start);
        tmcf->fiducial_part = part[0];
      }
      tmcf->Er   *= e_mks/(me_mks*c_mks*omega);
      tmcf->Ez   *= e_mks/(me_mks*c_mks*omega);
      tmcf->Bphi *= e_mks/(me_mks*omega);
    }
    *tau_start = tmcf->phase + PIx2*tmcf->frequency*tmcf->time_offset +
      tmcf->phase0;
    *Z_end = tmcf->length*(*kscale);
    X_offset = tmcf->k*tmcf->radial_offset*cos(tmcf->tilt);
    Y_offset = tmcf->k*tmcf->radial_offset*sin(tmcf->tilt);
    *X_limit = tmcf->k*tmcf->x_max;
    *Y_limit = tmcf->k*tmcf->y_max;
    X_aperture_center = tmcf->k*tmcf->dx;
    Y_aperture_center = tmcf->k*tmcf->dy;
    *accuracy = tmcf->accuracy;
    *n_steps = tmcf->n_steps;
    select_integrator(tmcf->method);
    return(derivatives_tmcf_mode);
  case T_CEPL:
    cep = field;
    *kscale = 1./(c_mks*cep->ramp_time);
    if (!cep->fiducial_part) {
      /* This is the fiducial particle--the tau offset tau0 is set 
       * so that its tau when it reaches the cavity center will be 
       * approximately 0. Several other quantities are calculated 
       * here also.
       */
      *kscale = cep->k = 1./(c_mks*cep->ramp_time);
      cep->fiducial_part = fiducial = select_fiducial(part, n_part, cep->fiducial);
      if (fiducial) {
        Po = P_central*(1+fiducial[5]);
        gamma = sqrt(1+sqr(Po));
        Pz = Po*sqrt(1-sqr(fiducial[0])-sqr(fiducial[2]));
        cep->tau0 = get_reference_phase(cep->phase_reference,
                                        -cep->k*gamma*(cep->length/(2*Pz) + fiducial[4]/Po ) );
      }
      else {
        /* phase to v=c */
        cep->tau0 = get_reference_phase((long)(cep->phase_reference+.5),
                                        -cep->k*cep->length/2 + z_start);
        cep->fiducial_part = part[0];
      }
    }
    cep->E_scaled = e_mks*cep->voltage*cep->ramp_time/
      (cep->gap*me_mks*c_mks);
    cep->E_static = e_mks*cep->static_voltage*cep->ramp_time/
      (cep->gap*me_mks*c_mks);
    cep->cos_tilt = cos(cep->tilt);
    cep->sin_tilt = sin(cep->tilt);
    *tau_start = cep->time_offset/cep->ramp_time + cep->tau0;
    *Z_end = cep->length*(*kscale);
    *accuracy = cep->accuracy;
    *n_steps = cep->n_steps;
    X_offset = Y_offset = 0;
    *X_limit = cep->k*cep->x_max;
    *Y_limit = cep->k*cep->y_max;
    X_aperture_center = cep->k*cep->dx;
    Y_aperture_center = cep->k*cep->dy;
    select_integrator(cep->method);
    return(derivatives_ce_plates);
  case T_TWPL:
    twp = field;
    *kscale = 1./(c_mks*twp->ramp_time);
    if (!twp->fiducial_part) {
      /* This is the fiducial particle--the tau offset tau0 is set 
       * so that its tau when it reaches the cavity center will be 
       * approximately omega*time_offset. Several other quantities are calculated 
       * here also.
       */
      *kscale = twp->k = 1./(c_mks*twp->ramp_time);
      twp->fiducial_part = fiducial = select_fiducial(part, n_part, twp->fiducial);
      if (fiducial) {
        Po = P_central*(1+fiducial[5]);
        gamma = sqrt(1+sqr(Po));
        Pz = Po*sqrt(1-sqr(fiducial[0])-sqr(fiducial[2]));
        twp->tau0 = get_reference_phase(twp->phase_reference,
                                        -twp->k*gamma*(twp->length/(2*Pz) + fiducial[4]/Po) );
      }
      else {
        /* phase to v=c */
        twp->tau0 = get_reference_phase((long)(twp->phase_reference+.5),
                                        -twp->k*twp->length/2 + z_start);
        twp->fiducial_part = part[0];
      }
    }
    twp->E_scaled = e_mks*twp->voltage*twp->ramp_time/
      (twp->gap*me_mks*c_mks);
    twp->E_static = e_mks*twp->static_voltage*twp->ramp_time/
      (twp->gap*me_mks*c_mks);
    twp->cos_tilt = cos(twp->tilt);
    twp->sin_tilt = sin(twp->tilt);
    *tau_start = twp->time_offset/twp->ramp_time + twp->tau0;
    *Z_end = twp->length*(*kscale);
    *accuracy = twp->accuracy;
    *n_steps = twp->n_steps;
    X_offset = Y_offset = 0;
    *X_limit = twp->k*twp->x_max;
    *Y_limit = twp->k*twp->y_max;
    X_aperture_center = twp->k*twp->dx;
    Y_aperture_center = twp->k*twp->dy;
#if defined(DEBUG)
    fprintf(stdout, "TWPL parameters:\n");
    fflush(stdout);
    fprintf(stdout, "l=%le acc=%le x_max=%le y_max=%le\n",
            twp->length, twp->accuracy, twp->x_max, twp->y_max);
    fflush(stdout);
    fprintf(stdout, "dx=%le dy=%le method=%s ramp_time=%le phiref=%ld\n",
            twp->dx, twp->dy, twp->method, twp->ramp_time,
            twp->phase_reference);
    fflush(stdout);
#endif
    select_integrator(twp->method);
    return(derivatives_tw_plates);
  case T_TWLA:
    twla = (TW_LINAC*)field;
    radial_limit = 1;
    *kscale = (omega=twla->frequency*PIx2)/c_mks;
    twla->alphaS = twla->alpha/(*kscale);
    if (twla->change_p0)
      change_p0 = 1;
    if (!twla->fiducial_part) {
      /* This is the fiducial particle--the phase offset phase0 is 
       * set so that its initial phase will be twla->phase + 
       * omega*twla->time_offset. Several other quantities are 
       * calculated here also.
       */
      twla->kz = *kscale/twla->beta_wave;
      twla->fiducial_part = fiducial = select_fiducial(part, n_part, twla->fiducial);
      if (fiducial) {
        Po = P_central*(1+fiducial[5]);
        gamma = sqrt(1+sqr(Po));
        /* find the phase offset, phase0 = -omega*t_fiducial , or else equivalent value
         * for phase reference */
        twla->phase0 = get_reference_phase(twla->phase_reference,
                                           -omega/c_mks*gamma*fiducial[4]/Po );
      }
      else {
        /* use v=c */
        twla->phase0 = get_reference_phase(twla->phase_reference, -omega/c_mks*z_start);
        twla->fiducial_part = part[0];
      }
    }
    Escale /= omega;
    Bscale /= omega;
    twla->ErS    = -twla->Ez*twla->kz/2 * Escale * twla->focussing / (*kscale);
    twla->BphiS  = -twla->Ez*omega/(2*sqr(c_mks)) * Bscale * twla->focussing / (*kscale);
    twla->EzS    = twla->Ez * Escale;
    twla->BsolS  = twla->B_solenoid * Bscale;
    /* calculate initial tau value, less omega*t: 
     *    tau_start = omega*(t_offset-t_fid)+phase 
     *              = omega*t_offset + phase + phase0
     */
    *tau_start = twla->phase + omega*twla->time_offset + twla->phase0;
    *Z_end = twla->length*(*kscale);
    X_offset = X_aperture_center = (*kscale)*twla->dx;
    Y_offset = Y_aperture_center = (*kscale)*twla->dy;
    *X_limit = (*kscale)*twla->x_max;
    *Y_limit = (*kscale)*twla->y_max;
    *accuracy = twla->accuracy;
    *n_steps = twla->n_steps;
    select_integrator(twla->method);
    derivatives    = derivatives_tw_linac;
    input_impulse  = input_impulse_tw_linac;
    output_impulse = output_impulse_tw_linac;
    return(derivatives_tw_linac);
  case T_TWMTA:
    twmta = field;
    radial_limit = 0;
    *kscale = (omega=twmta->frequency*PIx2)/c_mks;
    twmta->alphaS = twmta->alpha/(*kscale);
    if (!twmta->fiducial_part) {
      /* This is the fiducial particle--the phase offset phase0 is 
       * set so that its initial phase will be twmta->phase + 
       * omega*twmta->time_offset. Several other quantities are 
       * calculated here also.
       */
      twmta->kz = (*kscale)/twmta->beta_wave;
      gamma_w = 1/sqrt(1-sqr(twmta->beta_wave));
      twmta->ky = sqrt(sqr(twmta->kx)+sqr(twmta->kz/gamma_w));
      twmta->Kx = twmta->kx/(*kscale);
      twmta->Ky = twmta->ky/(*kscale);
      twmta->Kz = twmta->kz/(*kscale);
      /*
        fprintf(stdout, "TWMTA kx, ky, kz = %e, %e, %e\n", twmta->kx, twmta->ky, twmta->kz);
        fflush(stdout);
        fprintf(stdout, "      Ex, Ey, Ez = %e, %e, %e\n", twmta->ExS, twmta->EyS, twmta->EzS);
        fflush(stdout);
        fprintf(stdout, "      Bx, By, Bz = %e, %e\n", twmta->BxS, twmta->ByS);
        fflush(stdout);
        */
      twmta->fiducial_part = fiducial = select_fiducial(part, n_part, twmta->fiducial);
      if (fiducial) {
        /* find the phase offset, phase0 = -omega*t_fiducial , or else equivalent value
         * for phase reference */
        Po = P_central*(1+fiducial[5]);
        gamma = sqrt(1+sqr(Po));
        twmta->phase0 = get_reference_phase(twmta->phase_reference,
                                            -omega/c_mks*gamma*fiducial[4]/Po );
      }
      else {
        /* use v=c */
        twmta->phase0 = get_reference_phase(twmta->phase_reference, -omega/c_mks*z_start);
        fiducial = part[0];
      }
    }
    Escale /= omega;
    Bscale /= omega;
    twmta->EzS = twmta->Ez * Escale;
    twmta->ExS = twmta->EzS*twmta->kx*twmta->kz/(sqr(twmta->ky)-sqr(twmta->kx));
    twmta->EyS = -twmta->ExS*twmta->ky/twmta->kx;
    twmta->BxS = twmta->Ez/omega*Bscale*twmta->ky*
      (sqr(twmta->kx)-sqr(twmta->ky)+sqr(twmta->kz))/(sqr(twmta->ky)-sqr(twmta->kx));
    twmta->ByS = twmta->BxS*twmta->kx/twmta->ky;
    twmta->BsolS = twmta->Bsol*Bscale;
    /* calculate initial tau value, less omega*t: 
     *    tau_start = omega*(-t_fid)+phase 
     *              = phase + phase0
     */
    *tau_start = twmta->phase + twmta->phase0;
    *Z_end = twmta->length*(*kscale);
    X_offset = X_aperture_center = (*kscale)*twmta->dx;
    Y_offset = Y_aperture_center = (*kscale)*twmta->dy;
    *X_limit = (*kscale)*twmta->x_max;
    *Y_limit = (*kscale)*twmta->y_max;
    *accuracy = twmta->accuracy;
    *n_steps = twmta->n_steps;
    select_integrator(twmta->method);
    input_impulse  = input_impulse_twmta;
    output_impulse = output_impulse_twmta;
    return(derivatives_twmta);
  default:
    bomb("invalid mode type (set_up_derivatives", NULL);
    break;
  }
  exit(1);
}

void derivatives_rftmEz0(
    double *qp,       /* derivatives w.r.t. phase */
    double *q,        /* X,Y,Z,Px,Py,Pz */
    double tau 
    )
{
  double gamma, *P, *Pp;
  double E[3], BOverGamma[3], XYZ[3];
  RFTMEZ0 *rftmEz0;
  long iz, ir;
  double R, Zoffset, BphiOverRG, ErOverR, BrOverRG;
  double X, Y, Z, sinPhase, cosPhase;
  double B1, B2;
  
  derivCalls++;
  P  = q+3;
  gamma = sqrt(sqr(P[0])+sqr(P[1])+sqr(P[2])+1);
  
  /* X' = Px/gamma, etc. */
  qp[0] = P[0]/gamma;
  qp[1] = P[1]/gamma;
  qp[2] = P[2]/gamma;
  
  /* find coordinates relative to element center */
  XYZ[0] = q[0]-X_center;
  XYZ[1] = q[1]-Y_center;
  XYZ[2] = q[2]-Z_center;
  /* perform coordinate rotation into the element frame */
  rotate3(XYZ, partRot);
  X = XYZ[0];
  Y = XYZ[1];
  Z = XYZ[2] + Z_center;  /* Z is expressed relative to element start */
  R = sqrt(sqr(X)+sqr(Y));

  /* get scaled RF fields */
  rftmEz0 = field_global;
  iz = q[2]/rftmEz0->dZ;
  if (iz<0 || iz>=rftmEz0->nz) {
    E[0] = E[1] = E[2] = 0;
    BOverGamma[0] = BOverGamma[1] = BOverGamma[2] = 0;
  } else {
    if (iz==(rftmEz0->nz-1))
      iz -= 1;
    Zoffset = Z-iz*rftmEz0->dZ;
    sinPhase = sin(tau);
    cosPhase = cos(tau);
    E[2] = (rftmEz0->Ez[iz] + Zoffset*rftmEz0->dEzdZ[iz])*sinPhase*rftmEz0->Ez_peak;
    ErOverR = -(rftmEz0->dEzdZ[iz] + 
           (rftmEz0->dEzdZ[iz+1]-rftmEz0->dEzdZ[iz])/rftmEz0->dZ*Zoffset)/2*sinPhase*rftmEz0->Ez_peak;
    E[0] = ErOverR*X;
    E[1] = ErOverR*Y;
    BphiOverRG = E[2]/2/gamma*cosPhase;
    BOverGamma[0] = -BphiOverRG*Y;
    BOverGamma[1] = BphiOverRG*X;
    BOverGamma[2] = 0;
  }

  /* add scaled solenoid fields */
  if (rftmEz0->BzSol) {
    iz = (Z - rftmEz0->Z0Sol)/rftmEz0->dZSol;
    Zoffset = Z - (iz*rftmEz0->dZSol + rftmEz0->Z0Sol);
    ir = R/rftmEz0->dRSol;
    
    if (iz>=0 && iz<rftmEz0->nzSol && ir<rftmEz0->nrSol) {
      if (iz==rftmEz0->nzSol-1)
        iz -= 1;
      if (ir==rftmEz0->nrSol-1)
        ir -= 1;
      BrOverRG = 0;
      if (R>0) {
        /* compute Br/R */
        B1 = rftmEz0->BrSol[ir][iz] +
          (rftmEz0->BrSol[ir+1][iz] - rftmEz0->BrSol[ir][iz])/rftmEz0->dRSol*(R-ir*rftmEz0->dRSol);
        B2 = rftmEz0->BrSol[ir][iz+1] +
          (rftmEz0->BrSol[ir+1][iz+1] - rftmEz0->BrSol[ir][iz+1])/rftmEz0->dRSol*(R-ir*rftmEz0->dRSol);
        BrOverRG = (B1 + (B2-B1)*Zoffset/rftmEz0->dZSol)/gamma/R*rftmEz0->solenoidFactor;
        BOverGamma[0] += X*BrOverRG;
        BOverGamma[1] += Y*BrOverRG;
      }
      /* compute Bz/Gamma */
      B1 = rftmEz0->BzSol[ir][iz] +
        (rftmEz0->BzSol[ir+1][iz] - rftmEz0->BzSol[ir][iz])/rftmEz0->dRSol*(R-ir*rftmEz0->dRSol);
      B2 = rftmEz0->BzSol[ir][iz+1] +
        (rftmEz0->BzSol[ir+1][iz+1] - rftmEz0->BzSol[ir][iz+1])/rftmEz0->dRSol*(R-ir*rftmEz0->dRSol);
      BOverGamma[2] = 
        (B1 + (B2-B1)*Zoffset/rftmEz0->dZSol)/gamma*rftmEz0->solenoidFactor;
    }
  }
  /* E and BOverGama contain field components in the element frame for
   * the present particle position (which was expressed in the element frame
   * also).  Rotate these components back to the lab frame. 
   */
  rotate3(E, fieldRot);
  rotate3(BOverGamma, fieldRot);
  
  /* (Px,Py,Pz)' = (Ex,Ey,Ez) + (Px,Py,Pz)x(Bx,By,Bz)/gamma */
  Pp = qp+3;
  Pp[0] = -(E[0] + (P[1]*BOverGamma[2]-P[2]*BOverGamma[1]));
  Pp[1] = -(E[1] + (P[2]*BOverGamma[0]-P[0]*BOverGamma[2]));
  Pp[2] = -(E[2] + (P[0]*BOverGamma[1]-P[1]*BOverGamma[0]));
}

void derivatives_mapSolenoid(
    double *qp,       /* derivatives w.r.t. phase */
    double *q,        /* X,Y,Z,Px,Py,Pz */
    double tau 
    )
{
  double gamma, *P, *Pp;
  double BOverGamma[3];
  double XYZ[3];
  MAP_SOLENOID *mapSol;
  long iz, ir;
  double R, Zoffset, BphiOverRG, ErOverR, BrOverRG;
  double X, Y, Z, sinPhase, cosPhase;
  double B1, B2;
  
  derivCalls++;
  P  = q+3;
  gamma = sqrt(sqr(P[0])+sqr(P[1])+sqr(P[2])+1);
  
  /* X' = Px/gamma, etc. */
  qp[0] = P[0]/gamma;
  qp[1] = P[1]/gamma;
  qp[2] = P[2]/gamma;

  /* find coordinates relative to element center */
  XYZ[0] = q[0]-X_center;
  XYZ[1] = q[1]-Y_center;
  XYZ[2] = q[2]-Z_center;
  /* perform coordinate rotation */
  rotate3(XYZ, partRot);
  X = XYZ[0];
  Y = XYZ[1];
  Z = XYZ[2] + Z_center;
  R = sqrt(sqr(X)+sqr(Y));

  mapSol = field_global;
  iz = Z/mapSol->dz;
  Zoffset = Z - iz*mapSol->dz;
  ir = R/mapSol->dr;
  
  BOverGamma[0] = BOverGamma[1] = BOverGamma[2] = 0;
  if (iz>=0 && iz<mapSol->nz && ir<mapSol->nr) {
    if (iz==mapSol->nz-1)
        iz -= 1;
    if (ir==mapSol->nr-1)
      ir -= 1;
    BrOverRG = 0;
    if (R>0) {
      /* compute Br/R */
      B1 = mapSol->Br[ir][iz] +
        (mapSol->Br[ir+1][iz] - mapSol->Br[ir][iz])/mapSol->dr*(R-ir*mapSol->dr);
      B2 = mapSol->Br[ir][iz+1] +
        (mapSol->Br[ir+1][iz+1] - mapSol->Br[ir][iz+1])/mapSol->dr*(R-ir*mapSol->dr);
      BrOverRG = (B1 + (B2-B1)*Zoffset/mapSol->dz)/gamma/R*mapSol->factor;
      BOverGamma[0] += X*BrOverRG;
      BOverGamma[1] += Y*BrOverRG;
    }
    /* compute Bz/Gamma */
    B1 = mapSol->Bz[ir][iz] +
      (mapSol->Bz[ir+1][iz] - mapSol->Bz[ir][iz])/mapSol->dr*(R-ir*mapSol->dr);
    B2 = mapSol->Bz[ir][iz+1] +
      (mapSol->Bz[ir+1][iz+1] - mapSol->Bz[ir][iz+1])/mapSol->dr*(R-ir*mapSol->dr);
    BOverGamma[2] = (B1 + (B2-B1)*Zoffset/mapSol->dz)/gamma*mapSol->factor;
  }

  if ((mapSol->BxUniform || mapSol->ByUniform) &&
      (Z>mapSol->zMap0 && Z<mapSol->zMap1)) {
    /* add components from uniform field, if any */
    BOverGamma[0] += mapSol->BxUniformScaled/gamma;
    BOverGamma[1] += mapSol->ByUniformScaled/gamma;
  }
  /* BOverGama contains field components in the element frame for
   * the present particle position (which was expressed in the element frame
   * also).  Rotate these components back to the lab frame. 
   */
  rotate3(BOverGamma, fieldRot);
  
  /* (Px,Py,Pz)' = (Ex,Ey,Ez) + (Px,Py,Pz)x(Bx,By,Bz)/gamma */
  Pp = qp+3;
  Pp[0] = P[1]*BOverGamma[2]-P[2]*BOverGamma[1];
  Pp[1] = P[2]*BOverGamma[0]-P[0]*BOverGamma[2];
  Pp[2] = P[0]*BOverGamma[1]-P[1]*BOverGamma[0]; 
  /* Since the field components and momenta are in the lab frame, 
   * I don't need to perform a rotation of the forces. 
   */
}

void derivatives_tmcf_mode(
    double *qp,       /* derivatives w.r.t. tau */
    double *q,        /* X,Y,Z,Px,Py,Pz */
    double tau     
    )
{
    double gamma, *P, *Pp;
    double phase, sin_phase, cos_phase;
    double sin_tilt, cos_tilt;
    double X, Y, R, Z;
    static double E[3], B[3];
    TMCF_MODE *tmcf;

    derivCalls++;
    /* gamma = sqrt(P^2+1) */
    P  = q+3;
    gamma = sqrt(sqr(P[0])+sqr(P[1])+sqr(P[2])+1);

    /* X' = Px/gamma, etc. */
    qp[0] = P[0]/gamma;
    qp[1] = P[1]/gamma;
    qp[2] = P[2]/gamma;

    /* Calculate scaled RF fields E and B, for which                */
    /* (Px,Py,Pz)' = (Ex,Ey,Ez) + (Px,Py,Pz)x(Bx,By,Bz)/gamma       */
    Pp = qp+3;
    tmcf = field_global;
    if ((Z=q[2])<=Z_end && Z>=0) {
        /* The electric fields are multiplied by sin(w.t) and the
         * magnetic field by cos(w.t).  This is consistent with the
         * data dumped by SF02 and SHY. 
         */
        E[2] = tmcf->Ez*(sin_phase=sin(phase=(tau)));
        Y = q[1] + Y_offset;
        X = q[0] + X_offset;
        R = sqrt(sqr(X)+sqr(Y));
        E[0] = tmcf->Er*(cos_tilt=(R?X/R:0))*sin_phase;
        E[1] = tmcf->Er*(sin_tilt=(R?Y/R:0))*sin_phase;
        B[2] = 0;
        /* B is scaled by gamma here */
        B[1] =  tmcf->Bphi*(cos_phase=cos(phase))/gamma;
        B[0] = -B[1]*sin_tilt;
        B[1] *= cos_tilt;

        /* The sign of the electron is taken care of here. */
        Pp[0] = -(E[0] - P[2]*B[1]            );
        Pp[1] = -(E[1] + P[2]*B[0]            );
        Pp[2] = -(E[2] + (P[0]*B[1]-P[1]*B[0]));
        }
    else
        Pp[0] = Pp[1] = Pp[2] = 0;

    }

void derivatives_ce_plates(
    double *qp,       /* derivatives w.r.t. tau */
    double *q,        /* X,Y,Z,Px,Py,Pz */
    double tau 
    )
{
    register double gamma, *P, *Pp, t_ramp;
    double E, z;
    CE_PLATES *cep;

    derivCalls++;
    /* gamma = sqrt(P^2+1) */
    P  = q+3;
    gamma = sqrt(sqr(P[0])+sqr(P[1])+sqr(P[2])+1);

    /* X' = Px/gamma, etc. */
    qp[0] = P[0]/gamma;
    qp[1] = P[1]/gamma;
    qp[2] = P[2]/gamma;

    /* get scaled RF fields E and B, for which                */
    /* (Px,Py,Pz)' = (Ex,Ey,Ez) + (Px,Py,Pz)x(Bx,By,Bz)/gamma */
    Pp = qp+3;
    cep = field_global;
    if ((z=q[2])<=Z_end && z>=0) {
        /* Electrons are deflected to positive coordinates by a negative
         * electric field.
         */
        if ((t_ramp = tau)<1 && t_ramp>0)
            E = cep->E_scaled*t_ramp;
        else if (t_ramp<=0)
            E = 0;
        else
            E = cep->E_scaled;
        E += cep->E_static;

        /* The sign of the electron charge is taken care of here. */
        Pp[0] = -E*cep->cos_tilt;
        Pp[1] = -E*cep->sin_tilt;
        Pp[2] = 0;
        }
    else
        Pp[0] = Pp[1] = Pp[2] = 0;

    }

void derivatives_tw_plates(
    double *qp,       /* derivatives w.r.t. tau */
    double *q,        /* X,Y,Z,Px,Py,Pz */
    double tau 
    )
{
    double gamma, *P, *Pp, t_ramp;
    double E, z, B;
    TW_PLATES *twp;

    derivCalls++;
    /* gamma = sqrt(P^2+1) */
    P  = q+3;
    gamma = sqrt(sqr(P[0])+sqr(P[1])+sqr(P[2])+1);

    /* X' = Px/gamma, etc. */
    qp[0] = P[0]/gamma;
    qp[1] = P[1]/gamma;
    qp[2] = P[2]/gamma;

    /* get scaled RF fields E and B, for which                */
    /* (Px,Py,Pz)' = (Ex,Ey,Ez) + (Px,Py,Pz)x(Bx,By,Bz)/gamma */
    Pp = qp+3;
    twp = field_global;
    if ((z=q[2])<=Z_end && z>=0) {
        /* Electrons are deflected to positive coordinates by a negative
         * electric field.  The magnetic field deflects in the same direction.
         */
        if ((t_ramp = tau+z-Z_end)<1 && t_ramp>0)
            E = twp->E_scaled*t_ramp;
        else if (t_ramp<0)
            E = 0;
        else
            E = twp->E_scaled;
        E += twp->E_static;
        B  = -E/gamma;
        Pp[0] = -(E-P[2]*B)*twp->cos_tilt;
        Pp[1] = -(E-P[2]*B)*twp->sin_tilt;
        Pp[2] = -B*(P[0]*twp->cos_tilt+P[1]*B*twp->sin_tilt);
        }
    else
        Pp[0] = Pp[1] = Pp[2] = 0;

    }


void derivatives_tw_linac(
    double *qp,       /* derivatives w.r.t. tau */
    double *q,        /* X,Y,Z,Px,Py,Pz */
    double tau 
    )
{
    register double gamma, *P, *Pp;
    static double E[3], B[3];
    double phase, X, Y, Z, cos_phase, sin_phase, droop;
    TW_LINAC *twla;

    derivCalls++;
    /* gamma = sqrt(P^2+1) */
    P  = q+3;
    gamma = sqrt(sqr(P[0])+sqr(P[1])+sqr(P[2])+1);

    /* X' = Px/gamma, etc. */
    qp[0] = P[0]/gamma;
    qp[1] = P[1]/gamma;
    qp[2] = P[2]/gamma; 

    /* get scaled RF fields E and B, for which                */
    /* (Px,Py,Pz)' = (Ex,Ey,Ez) + (Px,Py,Pz)x(Bx,By,Bz)/gamma */
    Pp = qp+3;
    twla = field_global;
    droop = 1;
    if ((Z=q[2])<=Z_end && Z>=0) {
        E[2] = twla->EzS*(sin_phase=sin(phase=Z/twla->beta_wave - tau));
#if defined(DEBUG)
        if (starting_integration) {
            starting_integration = 0;
            initial_phase = phase;
            }
        else
            final_phase = phase;
#endif
        Y = q[1] + Y_offset;
        X = q[0] + X_offset;
        E[0]  = twla->ErS*(cos_phase=cos(phase));
        E[1]  = E[0]*Y;
        E[0] *= X;
        /* B is scaled by gamma here */
        B[2]  = twla->BsolS/gamma;
        B[1]  = twla->BphiS*cos_phase/gamma;
        B[0]  = -B[1]*Y;
        B[1] *= X;

        /* The sign of the electron is taken care of here. */
        if (twla->alphaS)
            droop = exp(-twla->alphaS*Z);
        Pp[0] = -(E[0] + P[1]*B[2] - P[2]*B[1])*droop;
        Pp[1] = -(E[1] + P[2]*B[0] - P[0]*B[2])*droop;
        Pp[2] = -(E[2] + P[0]*B[1] - P[1]*B[0])*droop;
        }
    else
        Pp[0] = Pp[1] = Pp[2] = 0;

    }

void input_impulse_tw_linac(double *P, double *q)
{
    register double gamma;
    static double B[3];
    TW_LINAC *twla;

    /* gamma = sqrt(P^2+1) */
    gamma = sqrt(sqr(P[0])+sqr(P[1])+sqr(P[2])+1);

    
    twla = (TW_LINAC*)field_global;
    if (twla->BsolS==0)
        return;
    
    B[0] = -twla->BsolS*q[0]/(2*gamma);
    B[1] = -twla->BsolS*q[1]/(2*gamma);
    P[0] +=  P[2]*B[1];
    P[1] += -P[2]*B[0];
    P[2] += -P[0]*B[1]+P[1]*B[0];
    }

void output_impulse_tw_linac(double *P, double *q)
{
    register double gamma;
    static double B[3];
    TW_LINAC *twla;

    /* gamma = sqrt(P^2+1) */
    gamma = sqrt(sqr(P[0])+sqr(P[1])+sqr(P[2])+1);

    
    twla = (TW_LINAC*)field_global;
    if (twla->BsolS==0)
        return;
    
    B[0] = twla->BsolS*q[0]/(2*gamma);
    B[1] = twla->BsolS*q[1]/(2*gamma);
    P[0] +=  P[2]*B[1];
    P[1] += -P[2]*B[0];
    P[2] += -P[0]*B[1]+P[1]*B[0];
    }

void derivatives_twmta(
    double *qp,       /* derivatives w.r.t. tau */
    double *q,        /* X,Y,Z,Px,Py,Pz */
    double tau 
    )
{
    register double gamma, *P, *Pp;
    static double E[3], B[3];
    double phase, X, Y, Z, droop;
    double sin_phase, cos_phase, cosh_KyY, sinh_KyY, cos_KxX, sin_KxX;
    TWMTA *twmta;

    derivCalls++;
    /* gamma = sqrt(P^2+1) */
    P  = q+3;
    gamma = sqrt(sqr(P[0])+sqr(P[1])+sqr(P[2])+1);

    /* X' = Px/gamma, etc. */
    qp[0] = P[0]/gamma;
    qp[1] = P[1]/gamma;
    qp[2] = P[2]/gamma;

    /* get scaled RF fields E and B, for which                */
    /* (Px,Py,Pz)' = (Ex,Ey,Ez) + (Px,Py,Pz)x(Bx,By,Bz)/gamma */
    Pp = qp+3;
    twmta = field_global;
    droop = 1;
    if ((Z=q[2])<=Z_end && Z>=0) {
        Y = q[1] + Y_offset;
        X = q[0] + X_offset;
        sin_phase = sin(phase=Z/twmta->beta_wave - tau);
        cos_phase = cos(phase);
#if defined(DEBUG)
        if (starting_integration) {
            starting_integration = 0;
            initial_phase = phase;
            }
        else
            final_phase = phase;
#endif
        cosh_KyY = cosh(twmta->Ky*Y);
        sinh_KyY = sinh(twmta->Ky*Y);
        cos_KxX  = cos(twmta->Kx*X);
        sin_KxX  = sin(twmta->Kx*X);
        E[0]  = twmta->ExS*sin_KxX*cosh_KyY*cos_phase;
        E[1]  = twmta->EyS*cos_KxX*sinh_KyY*cos_phase;
        E[2]  = twmta->EzS*cos_KxX*cosh_KyY*sin_phase;
        /* B is scaled by gamma here */
        B[0]  = twmta->BxS/gamma*cos_KxX*sinh_KyY*cos_phase;
        B[1]  = twmta->ByS/gamma*sin_KxX*cosh_KyY*cos_phase;
        B[2]  = twmta->BsolS/gamma;

        if (twmta->alphaS)
            droop = exp(-twmta->alphaS*Z);

        /* The sign of the electron is taken care of here. */
        Pp[0] = -(E[0] + P[1]*B[2] - P[2]*B[1])*droop;
        Pp[1] = -(E[1] + P[2]*B[0] - P[0]*B[2])*droop;
        Pp[2] = -(E[2] + P[0]*B[1] - P[1]*B[0])*droop;
        }
    else
        Pp[0] = Pp[1] = Pp[2] = 0;

    }

void input_impulse_twmta(double *P, double *q)
{
    register double gamma;
    static double B[3];
    TWMTA *twmta;

    /* gamma = sqrt(P^2+1) */
    gamma = sqrt(sqr(P[0])+sqr(P[1])+sqr(P[2])+1);

    
    twmta = (TWMTA*)field_global;
    if (twmta->BsolS==0)
        return;
    
    B[0] = -twmta->BsolS*q[0]/(2*gamma);
    B[1] = -twmta->BsolS*q[1]/(2*gamma);
    P[0] +=  P[2]*B[1];
    P[1] += -P[2]*B[0];
    P[2] += -P[0]*B[1]+P[1]*B[0];
    }

void output_impulse_twmta(double *P, double *q)
{
    register double gamma;
    static double B[3];
    TWMTA *twmta;

    /* gamma = sqrt(P^2+1) */
    gamma = sqrt(sqr(P[0])+sqr(P[1])+sqr(P[2])+1);

    
    twmta = (TWMTA*)field_global;
    if (twmta->BsolS==0)
        return;
    
    B[0] = twmta->BsolS*q[0]/(2*gamma);
    B[1] = twmta->BsolS*q[1]/(2*gamma);
    P[0] +=  P[2]*B[1];
    P[1] += -P[2]*B[0];
    P[2] += -P[0]*B[1]+P[1]*B[0];
    }

double exit_function(
    double *qp,
    double *q,
    double phase
    )
{
    double X, Y, dZ;

    if ((dZ=Z_end-q[2])>=0 && !limit_hit) {
        if (radial_limit) {
            if (X_limit) {
                X = fabs(q[0]-X_aperture_center);
                Y = fabs(q[1]-Y_aperture_center);
                if (sqr(X_limit)<(sqr(X)+sqr(Y))) {
                    X_out = X;
                    Y_out = Y;
                    Z_out = q[2];
                    limit_hit = 1;
                    }
                }
            }
        else {
            if ((X_limit && (X=fabs(q[0]-X_aperture_center))>X_limit) ||
                (Y_limit && (Y=fabs(q[1]-Y_aperture_center))>Y_limit) ) {
                X_out = X;
                Y_out = Y;
                Z_out = q[2];
                limit_hit = 1;
                }
            }
        }
    if (q[5]<=0)  /* Pz */
        return(0.0);
    return(dZ);
    }

#define FID_MEDIAN 0
#define FID_MINIMUM 1
#define FID_MAXIMUM 2
#define FID_AVERAGE 3
#define FID_FIRST 4
#define FID_LIGHT 5
#define N_KNOWN_MODES 6
static char *known_mode[N_KNOWN_MODES] = {
    "median", "minimum", "maximum", "average", "first", "light"
    } ;

double *select_fiducial(double **part, long n_part, char *var_mode_in)
{
  long i;
  long i_best, i_var, fid_mode;
  double value, best_value, sum;
  char *var, *mode, *var_mode, *ptr;

  log_entry("select_fiducial");

  fid_mode = FID_AVERAGE;
  i_var = -1;    /* t, p average */
  cp_str(&var_mode, var_mode_in==NULL?"average":var_mode_in);

  if (n_part<=1) {
    log_exit("select_fiducial");
    return(part[0]);
  }

  if (strlen(var_mode)!=0) {
    if ((ptr=strchr(var_mode, ',')) && (var=get_token(var_mode))) {
      if (*var!='t' && *var!='p')
        bomb("no known variable listed for fiducial--must be 't' or 'p'\n", NULL);
      if (*var=='t')
        i_var = 4;
      else
        i_var = 5;
    } 
    if ((mode=get_token(var_mode))) {
      if ((fid_mode=match_string(mode, known_mode, N_KNOWN_MODES, 0))<0) {
        fputs("error: no known mode listed for fiducialization--must be one of:\n", stdout);
        for (i=0; i<N_KNOWN_MODES; i++)
          fprintf(stdout, "    %s\n", known_mode[i]);
          fflush(stdout);
        exit(1);
      }
    }
  }
  if (i_var==-1 && fid_mode!=FID_AVERAGE) {
    fprintf(stdout, "unless you use average mode for fiducialization, you must specify t or p coordinate.\n");
    fflush(stdout);
    exit(1);
  }
  if (i_var==-1) {
    double deltaAve, tAve;
    static double coord[6];
    for (i=deltaAve=tAve=0; i<n_part; i++) {
      deltaAve += part[i][5];
      tAve += part[i][4];
    }
    deltaAve /= n_part;
    tAve /= n_part;
    coord[0] = coord[1] = coord[2] = coord[3] = 0;
    coord[4] = tAve;
    coord[5] = deltaAve;
    return coord;
  }
  else {
    switch (fid_mode) {
    case FID_MEDIAN:
      if ((i_best = find_median_of_row(&best_value, part, i_var, n_part))<0)
        bomb("error: computation of median failed (select_fiducial)", NULL);
      break;
    case FID_MINIMUM:
      i_best = 0;
      best_value = DBL_MAX;
      for (i=0; i<n_part; i++) {
        if (best_value>(value=part[i][i_var])) {
          best_value = value;
          i_best = i;
        }
      }
      break;
    case FID_MAXIMUM:
      i_best = 0;
      best_value = -DBL_MAX;
      for (i=0; i<n_part; i++) {
        if (best_value<(value=part[i][i_var])) {
          best_value = value;
          i_best = i;
        }
      }
      break;
    case FID_AVERAGE:
      for (i=sum=0; i<n_part; i++)
        sum += part[i][i_var];
      sum /= n_part;
      best_value = DBL_MAX;
      for (i=0; i<n_part; i++) {
        if (best_value>(value=fabs(sum-part[i][i_var]))) {
          best_value = value;
          i_best = i;
        }
      }
      break;
    case FID_FIRST:
      i_best = 0;
      break;
    case FID_LIGHT:
      /* special case--return NULL pointer to indicate that phasing is to v=c */
      log_exit("select_fiducial");
      return(NULL);
    default:
      bomb("apparent programming error in select_fiducial", NULL);
      break;
    }
    if (i_best<0 || i_best>n_part)
      bomb("fiducial particle selection returned invalid value!", NULL);
  }
  
  log_exit("select_fiducial");
  return(part[i_best]);
}


#define RUNGE_KUTTA 0
#define BULIRSCH_STOER 1
#define NA_RUNGE_KUTTA 2
#define MODIFIED_MIDPOINT 3
#define N_METHODS 4
static char *method[N_METHODS] = {
    "runge-kutta", "bulirsch-stoer", "non-adaptive runge-kutta", "modified midpoint"
    } ;

void select_integrator(char *desired_method)
{
    long i;

    log_entry("select_integrator");

#if defined(DEBUG)
    fprintf(stdout, "select_integrator called with pointer %lx\n", (long)(desired_method));
    fflush(stdout);
    fprintf(stdout, "this translates into string %s\n", desired_method);
    fflush(stdout);
#endif
    switch (match_string(desired_method, method, N_METHODS, 0)) {
      case RUNGE_KUTTA:
        fprintf(stdout, "Warning: adaptive integrator chosen.  \"non-adaptive runge-kutta\" is recommended.\n");
        fflush(stdout);
        integrator = rk_odeint3;
        break;
      case BULIRSCH_STOER:
        fprintf(stdout, "Warning: adaptive integrator chosen.  \"non-adaptive runge-kutta\" is recommended.\n");
        fflush(stdout);
        integrator = bs_odeint3;
        break;
      case NA_RUNGE_KUTTA:
        integrator = rk_odeint3_na;
        break;
      case MODIFIED_MIDPOINT:
        integrator = mmid_odeint3_na;
        break;
      default:
        fprintf(stdout, "error: unknown integration method %s requested.\n", desired_method);
        fflush(stdout);
        fprintf(stdout, "Available methods are:\n");
        fflush(stdout);
        for (i=0; i<N_METHODS; i++)
            fprintf(stdout, "    %s\n", method[i]);
            fflush(stdout);
        exit(1);
        break;
        }
    log_exit("select_integrator");
    }

void setupRftmEz0FromFile(RFTMEZ0 *rftmEz0, double frequency, double length)
{
  SDDS_DATASET SDDSin;
  long i;
  double *z;
  double k, omega, Ez, EzMax;
  
  if (!rftmEz0)
    bomb("setupRftmEz0FromFile: NULL pointer passed.", NULL);

  if (rftmEz0->initialized)
    return;
  
  /* verify input parameters */
  if (!rftmEz0->inputFile || !rftmEz0->zColumn || !rftmEz0->EzColumn)
    bomb("RFTMEZ0 must have INPUTFILE, ZCOLUMN, and EZCOLUMN specified.", NULL);
  if (rftmEz0->radial_order!=1)
    bomb("RFTMEZ0 restricted to radial_order=1 at present", NULL);

  /* read (z, Ez) data from the file  */
  if (!SDDS_InitializeInputFromSearchPath(&SDDSin, rftmEz0->inputFile) || SDDS_ReadPage(&SDDSin)!=1) {
    fprintf(stderr, "Error: unable to open or read RFTMEZ0 file %s\n", 
            rftmEz0->inputFile);
    exit(1);
  }
  if ((rftmEz0->nz=SDDS_RowCount(&SDDSin))<0 || rftmEz0->nz<2) {
    fprintf(stderr, "Error: no data or insufficient data in RFTMEZ0 file %s\n", 
            rftmEz0->inputFile);
    exit(1);
  }
  if (SDDS_CheckColumn(&SDDSin, rftmEz0->zColumn, "m", SDDS_ANY_FLOATING_TYPE,
                       stderr)!=SDDS_CHECK_OK) {
    fprintf(stderr, "Error: problem with column %s in RFTMEZ0 file %s.  Check existence, type, and units.\n",
            rftmEz0->zColumn, rftmEz0->inputFile);
    exit(1);
  }
  if (SDDS_CheckColumn(&SDDSin, rftmEz0->EzColumn, NULL, SDDS_ANY_FLOATING_TYPE,
                       stderr)!=SDDS_CHECK_OK) {
    fprintf(stderr, "Error: problem with column %s in RFTMEZ0 file %s.  Check existence and type.\n",
            rftmEz0->zColumn, rftmEz0->inputFile);
    exit(1);
  }
  if (!(z=SDDS_GetColumnInDoubles(&SDDSin, rftmEz0->zColumn)) ||
      !(rftmEz0->Ez=SDDS_GetColumnInDoubles(&SDDSin, rftmEz0->EzColumn))) {
    fprintf(stderr, "Error: problem retrieving RFTMEZ0 data from SDDS structure for file %s\n",
            rftmEz0->inputFile);
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  rftmEz0->z0 = z[0];

  if (!SDDS_Terminate(&SDDSin)) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  if (!(rftmEz0->dEzdZ=SDDS_Malloc(sizeof(*(rftmEz0->dEzdZ))*rftmEz0->nz))) {
    fprintf(stderr, "Error: memory allocation failure setting up RFTMEZ0 file %s\n",
            rftmEz0->inputFile);
    exit(1);
  }

  /* compute dEz/dz for use in computing Er(z) */
  if (!analyzeSpacing(z, rftmEz0->nz, &rftmEz0->dz, stderr)) {
    fprintf(stderr, "Error: problem with z points from RFTMEZ0 file %s\n",
            rftmEz0->inputFile);
    exit(1);
  }
  
  if (fabs((length-(z[rftmEz0->nz-1]-z[0]))/rftmEz0->dz)>1e-6) {
    fprintf(stderr, "Error: declared length and length from fields from RFTMEZ0 file %s do not agree to 1/10^6 tolerance\n",
            rftmEz0->inputFile);
    exit(1);
  }
  free(z);

  /* find maximum Ez */
  EzMax = 0;
  for (i=0; i<rftmEz0->nz; i++) {
    if ((Ez=fabs(rftmEz0->Ez[i]))>EzMax)
      EzMax = Ez;
  }
  
  /* perform scaling and normalization (scaling by Ez_peak occurs when derivatives are computed) */
  k = (PIx2*frequency)/c_mks;
  rftmEz0->dZ = rftmEz0->dz*k;
  fflush(stdout);
  if (EzMax)
    for (i=0; i<rftmEz0->nz; i++)
      rftmEz0->Ez[i] *= e_mks/(me_mks*c_mks*PIx2*frequency)/EzMax;
  else
    for (i=0; i<rftmEz0->nz; i++)
      rftmEz0->Ez[i] = 0;

  /* take derivative except at endpoints */
  for (i=1; i<rftmEz0->nz-1; i++)
    rftmEz0->dEzdZ[i] = (rftmEz0->Ez[i+1]-rftmEz0->Ez[i-1])/(2*rftmEz0->dZ);

  /* for endpoints, assume field is zero outside region of data */
  rftmEz0->dEzdZ[0] = rftmEz0->Ez[1]/(2*rftmEz0->dZ);
  rftmEz0->dEzdZ[rftmEz0->nz-1] = rftmEz0->Ez[rftmEz0->nz-1]/(2*rftmEz0->dZ);

}

void setupRftmEz0SolenoidFromFile(RFTMEZ0 *rftmEz0, double length, double k) 
{
  SDDS_DATASET SDDSin;
  double *z, *r, dz, dr, *rTemp, z0;
  long page, ir, iz;
  
  if (rftmEz0->initialized)
    return;
  
  rftmEz0->BrSol = rftmEz0->BzSol = NULL;
  
  /* check for solenoid input file */
  if (!rftmEz0->solenoidFile)
    return;
  
  if (!rftmEz0->solenoid_zColumn || !rftmEz0->solenoid_rColumn ||
      !rftmEz0->solenoidBzColumn || !rftmEz0->solenoidBrColumn) 
    SDDS_Bomb("missing column name for solenoid for RFTMEZ0 element");

  if (!SDDS_InitializeInputFromSearchPath(&SDDSin, rftmEz0->solenoidFile)) {
    fprintf(stderr, "Error: unable to open or read RFTMEZ0 solenoid file %s\n", 
            rftmEz0->solenoidFile);
    exit(1);
  }
  if (SDDS_CheckColumn(&SDDSin, rftmEz0->solenoid_zColumn, "m", SDDS_ANY_FLOATING_TYPE,
                       stderr)!=SDDS_CHECK_OK) {
    fprintf(stderr, "Error: problem with column %s in RFTMEZ0 solenoid file %s.  Check existence, type, and units.\n",
            rftmEz0->solenoid_zColumn, rftmEz0->solenoidFile);
    exit(1);
  }
  if (SDDS_CheckColumn(&SDDSin, rftmEz0->solenoid_rColumn, "m", SDDS_ANY_FLOATING_TYPE,
                       stderr)!=SDDS_CHECK_OK) {
    fprintf(stderr, "Error: problem with column %s in RFTMEZ0 solenoid file %s.  Check existence, type, and units.\n",
            rftmEz0->solenoid_rColumn, rftmEz0->solenoidFile);
    exit(1);
  }

  if (SDDS_CheckColumn(&SDDSin, rftmEz0->solenoidBzColumn, "T", SDDS_ANY_FLOATING_TYPE,
                       stderr)!=SDDS_CHECK_OK) {
    fprintf(stderr, "Error: problem with column %s in RFTMEZ0 solenoid file %s.  Check existence, type, and units.\n",
            rftmEz0->solenoidBzColumn, rftmEz0->solenoidFile);
    exit(1);
  }
  if (SDDS_CheckColumn(&SDDSin, rftmEz0->solenoidBrColumn, "T", SDDS_ANY_FLOATING_TYPE,
                       stderr)!=SDDS_CHECK_OK) {
    fprintf(stderr, "Error: problem with column %s in RFTMEZ0 solenoid file %s.  Check existence, type, and units.\n",
            rftmEz0->solenoidBzColumn, rftmEz0->solenoidFile);
    exit(1);
  }

  z = NULL;
  r = NULL;
  while ((page=SDDS_ReadPage(&SDDSin))>0) {
    if (page==1) {
      if ((rftmEz0->nzSol=SDDS_RowCount(&SDDSin))<0 || rftmEz0->nzSol<2) {
        fprintf(stderr, "Error: no data or insufficient data in RFTMEZ0 solenoid file %s\n", 
                rftmEz0->solenoidFile);
        exit(1);
      }
      if (!(z=SDDS_GetColumnInDoubles(&SDDSin, rftmEz0->solenoid_zColumn))) {
        fprintf(stderr, "Error: problem getting z data from RFTMEZ0 solenoid file %s\n",
                rftmEz0->solenoidFile);
        exit(1);
      }
      if (!analyzeSpacing(z, rftmEz0->nzSol, &dz, stderr)) {
        fprintf(stderr, "Problem with z spacing of solenoid data from RFTMEZ0 file %s (page %ld)\n",
                rftmEz0->solenoidFile, page);
        exit(1);
      }
      z0 = z[0];
      if (fabs((length-(z[rftmEz0->nzSol-1]-z[0]))/dz)>1e-6) {
        fprintf(stderr, "Error: declared length and length from fields from RFTMEZ0 solenoid file %s do not agree to 1/10^6 tolerance\n",
                rftmEz0->solenoidFile);
        
        exit(1);
      }
      free(z);
    } else {
      if (rftmEz0->nzSol!=SDDS_RowCount(&SDDSin)) {
        fprintf(stderr, "Error: page %ld of RFTMEZ0 file %s has only %ld rows (%ld expected)\n",
                page, rftmEz0->solenoidFile, SDDS_RowCount(&SDDSin), rftmEz0->nzSol);
        exit(1);
      }
    }

    if (!(rTemp = SDDS_GetColumnInDoubles(&SDDSin, rftmEz0->solenoid_rColumn))) {
      fprintf(stderr, "Error: problem getting r data from RFTMEZ0 solenoid file %s\n",
              rftmEz0->solenoidFile);
      exit(1);
    }
    if (!(r = SDDS_Realloc(r, sizeof(*r)*page)))
      SDDS_Bomb("memory allocation failure (setupRftmEz0SolenoidFromFile)");
    r[page-1] = rTemp[0];
    free(rTemp);
    
    if (!(rftmEz0->BzSol = SDDS_Realloc(rftmEz0->BzSol, sizeof(*rftmEz0->BzSol)*page)) ||
        !(rftmEz0->BrSol = SDDS_Realloc(rftmEz0->BrSol, sizeof(*rftmEz0->BrSol)*page)) )
      SDDS_Bomb("memory allocation failure (setupRftmEz0SolenoidFromFile)");
    if (!(rftmEz0->BzSol[page-1] = 
          SDDS_GetColumnInDoubles(&SDDSin, rftmEz0->solenoidBzColumn)) ||
        !(rftmEz0->BrSol[page-1] = 
          SDDS_GetColumnInDoubles(&SDDSin, rftmEz0->solenoidBrColumn)) )  {
       fprintf(stderr, "Error: problem getting field data from RFTMEZ0 solenoid file %s\n",
               rftmEz0->solenoidFile);
       exit(1);
    }
    rftmEz0->nrSol = page;
  }

  SDDS_Terminate(&SDDSin);
  
  if (!analyzeSpacing(r, rftmEz0->nrSol, &dr, stderr)) {
    fprintf(stderr, "Problem with r spacing of solenoid data from RFTMEZ0 file %s\nr values are:\n",
            rftmEz0->solenoidFile);
    for (ir=0; ir<rftmEz0->nrSol; ir++)
      fprintf(stderr, "%le\n", r[ir]);
    exit(1);
  }
  free(r);
  
  rftmEz0->dRSol = k*dr;
  rftmEz0->dZSol = k*dz;
  rftmEz0->Z0Sol = k*(z0-rftmEz0->z0);

  /* perform scaling */
  for (ir=0; ir<rftmEz0->nrSol; ir++) 
    for (iz=0; iz<rftmEz0->nzSol; iz++)  {
      rftmEz0->BrSol[ir][iz] *= e_mks/(me_mks*k*c_mks);
      rftmEz0->BzSol[ir][iz] *= e_mks/(me_mks*k*c_mks);
    }
}

long analyzeSpacing(double *z, long nz, double *dzReturn, FILE *fpError) 
{
  double dzMin, dzMax, dz;
  long i;
  
  if (nz<3) {
    if (fpError)
      fprintf(fpError, "Problem with point spacing: less than 3 points\n");
    return 0;
  }
  
  dzMin = dzMax = z[1]-z[0];
  for (i=2; i<nz; i++) {
    if ((dz=z[i]-z[i-1])<dzMin)
      dzMin = dz;
    if (dz>dzMax)
      dzMax = dz;
  }
  if (dzMin<=0 || dzMax<=0) {
    if (fpError)
      fprintf(fpError, "Error: points are not monotonically increasing\n");
    return 0;
  }
  if (fabs(1-dzMin/dzMax)>0.0001) {
    if (fpError)
      fprintf(fpError, "Error: spacing of points is not uniform (0.01% tolerance)\n");
    return 0;
  }
  *dzReturn = (z[nz-1]-z[0])/(nz-1.0);
  return 1;
}

void setupMapSolenoidFromFile(MAP_SOLENOID *mapSol, double length)
{
  SDDS_DATASET SDDSin;
  double *z, *r, dz, dr, *rTemp, z0;
  long page, ir, iz;
  
  if (mapSol->initialized)
    return;
  mapSol->initialized = 1;
  
  mapSol->Br = mapSol->Bz = NULL;
  
  /* check for solenoid input file */
  if (!mapSol->inputFile)
    return;
  
  if (!mapSol->zColumn || !mapSol->rColumn ||
      !mapSol->BzColumn || !mapSol->BrColumn) 
    SDDS_Bomb("missing column name for solenoid for MAPSOLENOID element");

  if (!SDDS_InitializeInputFromSearchPath(&SDDSin, mapSol->inputFile)) {
    fprintf(stderr, "Error: unable to open or read MAPSOLENOID file %s\n", 
            mapSol->inputFile);
    exit(1);
  }
  if (SDDS_CheckColumn(&SDDSin, mapSol->zColumn, "m", SDDS_ANY_FLOATING_TYPE,
                       stderr)!=SDDS_CHECK_OK) {
    fprintf(stderr, "Error: problem with column %s in MAPSOLENOID file %s.  Check existence, type, and units.\n",
            mapSol->zColumn, mapSol->inputFile);
    exit(1);
  }
  if (SDDS_CheckColumn(&SDDSin, mapSol->rColumn, "m", SDDS_ANY_FLOATING_TYPE,
                       stderr)!=SDDS_CHECK_OK) {
    fprintf(stderr, "Error: problem with column %s in MAPSOLENOID file %s.  Check existence, type, and units.\n",
            mapSol->rColumn, mapSol->inputFile);
    exit(1);
  }

  if (SDDS_CheckColumn(&SDDSin, mapSol->BzColumn, "T", SDDS_ANY_FLOATING_TYPE,
                       stderr)!=SDDS_CHECK_OK) {
    fprintf(stderr, "Error: problem with column %s in MAPSOLENOID file %s.  Check existence, type, and units.\n",
            mapSol->BzColumn, mapSol->inputFile);
    exit(1);
  }
  if (SDDS_CheckColumn(&SDDSin, mapSol->BrColumn, "T", SDDS_ANY_FLOATING_TYPE,
                       stderr)!=SDDS_CHECK_OK) {
    fprintf(stderr, "Error: problem with column %s in MAPSOLENOID file %s.  Check existence, type, and units.\n",
            mapSol->BzColumn, mapSol->inputFile);
    exit(1);
  }

  z = NULL;
  r = NULL;
  while ((page=SDDS_ReadPage(&SDDSin))>0) {
    if (page==1) {
      if ((mapSol->nz=SDDS_RowCount(&SDDSin))<0 || mapSol->nz<2) {
        fprintf(stderr, "Error: no data or insufficient data in MAPSOLENOID file %s\n", 
                mapSol->inputFile);
        exit(1);
      }
      if (!(z=SDDS_GetColumnInDoubles(&SDDSin, mapSol->zColumn))) {
        fprintf(stderr, "Error: problem getting z data from MAPSOLENOID file %s\n",
                mapSol->inputFile);
        exit(1);
      }
      if (!analyzeSpacing(z, mapSol->nz, &dz, stderr)) {
        fprintf(stderr, "Problem with z spacing of solenoid data from MAPSOLENOID file %s (page %ld)\n",
                mapSol->inputFile, page);
        exit(1);
      }
      if (fabs((length-(z[mapSol->nz-1]-z[0]))/dz)>1e-6) {
        fprintf(stderr, "Error: declared length and length from fields from MAP_SOLENOID file %s do not agree to 1/10^6 tolerance\n",
                mapSol->inputFile);
        
        exit(1);
      }
      free(z);
    } else {
      if (mapSol->nz!=SDDS_RowCount(&SDDSin)) {
        fprintf(stderr, "Error: page %ld of MAPSOLENOID file %s has only %ld rows (%ld expected)\n",
                page, mapSol->inputFile, SDDS_RowCount(&SDDSin), mapSol->nz);
        exit(1);
      }
    }

    if (!(rTemp = SDDS_GetColumnInDoubles(&SDDSin, mapSol->rColumn))) {
      fprintf(stderr, "Error: problem getting r data from MAPSOLENOID file %s\n",
              mapSol->inputFile);
      exit(1);
    }
    if (!(r = SDDS_Realloc(r, sizeof(*r)*page)))
      SDDS_Bomb("memory allocation failure (setupmapSolSolenoidFromFile)");
    r[page-1] = rTemp[0];
    free(rTemp);
    
    if (!(mapSol->Bz = SDDS_Realloc(mapSol->Bz, sizeof(*mapSol->Bz)*page)) ||
        !(mapSol->Br = SDDS_Realloc(mapSol->Br, sizeof(*mapSol->Br)*page)) )
      SDDS_Bomb("memory allocation failure (setupmapSolSolenoidFromFile)");
    if (!(mapSol->Bz[page-1] = 
          SDDS_GetColumnInDoubles(&SDDSin, mapSol->BzColumn)) ||
        !(mapSol->Br[page-1] = 
          SDDS_GetColumnInDoubles(&SDDSin, mapSol->BrColumn)) )  {
       fprintf(stderr, "Error: problem getting field data from MAPSOLENOID file %s\n",
               mapSol->inputFile);
       exit(1);
    }
    mapSol->nr = page;
  }

  SDDS_Terminate(&SDDSin);
  
  if (!analyzeSpacing(r, mapSol->nr, &dr, stderr)) {
    fprintf(stderr, "Problem with r spacing of solenoid data from MAPSOLENOID file %s\nr values are:\n",
            mapSol->inputFile);
    for (ir=0; ir<mapSol->nr; ir++)
      fprintf(stderr, "%le\n", r[ir]);
    exit(1);
  }
  free(r);
  
  /* perform scaling */
  for (ir=0; ir<mapSol->nr; ir++) 
    for (iz=0; iz<mapSol->nz; iz++)  {
      mapSol->Br[ir][iz] *= e_mks/(me_mks*c_mks);
      mapSol->Bz[ir][iz] *= e_mks/(me_mks*c_mks);
    }

  mapSol->dr = dr;
  mapSol->dz = dz;
}

