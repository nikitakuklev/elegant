/* Copyright 1994 by Michael Borland and Argonne National Laboratory,
 * all rights reserved.
 */
/* file: simple_rfca.c
 * contents: simple_rf_cavity()
 *
 * Michael Borland, 1991
 */
#include "mdb.h"
#include "track.h"

static char *fiducialModeChoice[4] = {
    "light", "tmean", "first", "pmaximum",
    };

unsigned long parseFiducialMode(char *modeString)
{
    long code;

    if (!modeString)
        return FID_MODE_TMEAN;
    switch (code=match_string(modeString, fiducialModeChoice, 4, 0)) {
      case 0:
      case 1:
      case 2:
      case 3:
        return FID_MODE_LIGHT<<code;
      default:
        return 0;
        }
    }

double findFiducialTime(double **part, long np, double s0, double sOffset,
                        double p0, unsigned long mode)
{
    if (mode&FID_MODE_LIGHT)
        return (s0+sOffset)/c_mks;
    if (!np || mode&FID_MODE_FIRST)
        return (part[0][4]+sOffset)/(c_mks*beta_from_delta(p0, np?part[0][5]:0.0));
    if (mode&FID_MODE_PMAX) {
        long ibest, i;
        double best;
        best = part[0][5];
        ibest = 0;
        for (i=1; i<np; i++)
            if (best<part[i][5]) {
                best = part[i][5];
                ibest = i;
                }
        return (part[ibest][4]+sOffset)/(c_mks*beta_from_delta(p0, part[ibest][5]));
        }
    if (mode&FID_MODE_TMEAN) {
        double tsum;
        long ip;
        
        for (ip=tsum=0; ip<np; ip++) {
            tsum += (part[ip][4]+sOffset)/(c_mks*beta_from_delta(p0, part[ip][5]));
            }
        return tsum/np;
        }
    bomb("invalid fiducial mode in findFiducialTime", NULL);
	return(1);
    }


long simple_rf_cavity(
    double **part, long np, RFCA *rfca, double **accepted, double *P_central, double zEnd
    )
{
    long ip, same_dgamma;
    double timeOffset;
    double P, gamma, dgamma, phase, length, volt, To;
    double *coord, t, t0, omega, beta_i, tau, dt;
    static long been_warned = 0;
#ifdef DEBUG
    static FILE *fplog = NULL;
#endif

    log_entry("simple_rf_cavity");

    if (!been_warned) {        
        if (rfca->freq<1e3 && rfca->freq)  {
            fprintf(stderr, "\7\7\7warning: your RFCA frequency is less than 1kHz--this may be an error\n");
            been_warned = 1;
            }
        if (fabs(rfca->volt)<100 && rfca->volt) {
            fprintf(stderr, "\7\7\7warning: your RFCA voltage is less than 100V--this may be an error\n");
            been_warned = 1;
            }
        if (been_warned) {
            fprintf(stderr, "units of parameters for RFCA are as follows:\n");
            print_dictionary_entry(stderr, T_RFCA, 0);
            }
        }

    if (!part)
        bomb("NULL particle data pointer (simple_rf_cavity)", NULL);
    for (ip=0; ip<np; ip++)
        if (!part[ip]) {
            fprintf(stderr, "NULL pointer for particle %ld (simple_rf_cavity)\n", ip);
            abort();
            }
    if (!rfca)
        bomb("NULL rfca pointer (simple_rf_cavity)", NULL);

#ifdef DEBUG
    if (!fplog) {
        fplog = fopen_e("rfca.debug", "w", 0);
        fprintf(fplog, "11\nnp\nlength\nvoltage\nphase\nQ\nt0\nt\n$gD$rp\nR$b12$n\nR$b22$n\n\n\n1000\n");
        }
    fprintf(fplog, "%ld  %le  %le  %le  %le ",
        np, rfca->length, rfca->volt, rfca->phase, rfca->Q);
#endif

    if (np<=0) {
        log_exit("simple_rf_cavity");
        return(np);
        }

    if (rfca->change_t && rfca->Q)
        bomb("incompatible RF cavity parameters: change_t!=0 && Q!=0", NULL);

    length = rfca->length;

    if (rfca->volt==0) {
        if (rfca->length) {
            for (ip=0; ip<np; ip++) {
                coord = part[ip];
                coord[0] += coord[1]*length;
                coord[2] += coord[3]*length;
                coord[4] += length;
                }
            }
        log_exit("simple_rf_cavity");
        return(np);
        }

    omega = PIx2*rfca->freq;
    volt  = rfca->volt/(1e6*me_mev);
    if (omega)
        tau = rfca->Q/omega;
    else
        tau = 0;

    if (rfca->phase_reference==0) 
        rfca->phase_reference = unused_phase_reference();

    switch (get_phase_reference(&phase, rfca->phase_reference)) {
        case REF_PHASE_RETURNED:
            break;
        case REF_PHASE_NOT_SET:
        case REF_PHASE_NONEXISTENT:
            if (!rfca->fiducial_seen) {
                unsigned long mode;
                if (!(mode = parseFiducialMode(rfca->fiducial)))
                    bomb("invalid fiducial mode for RFCA element", NULL);
                t0 = findFiducialTime(part, np, zEnd-length, length/2, *P_central, mode);
                rfca->phase_fiducial = -omega*t0;
                rfca->fiducial_seen = 1;
                }
            set_phase_reference(rfca->phase_reference, phase=rfca->phase_fiducial);
            break;
        default:
            bomb("unknown return value from get_phase_reference()", NULL);
            break;
        }

    if (omega) {
        t0 = -rfca->phase_fiducial/omega;
        To = PIx2/omega;
        }
    else
        t0 = To = 0;
    phase += rfca->phase*PI/180.0;

    same_dgamma = 0;
    if (omega==0 && tau==0) {
        dgamma = volt*sin(phase);
        same_dgamma = 1;
        }

    timeOffset = 0;
    if (omega && rfca->change_t) {
      coord = part[0];
      P     = *P_central*(1+coord[5]);
      beta_i = P/(gamma=sqrt(sqr(P)+1));
      t     = coord[4]/(c_mks*beta_i);
      if (omega!=0 && t>(0.9*To) && rfca->change_t)
        timeOffset = ((long)(t/To+0.5))*To;
    }
    
    for (ip=0; ip<np; ip++) {
        coord = part[ip];
        if (length) {
          /* apply initial drift */
          coord[0] += coord[1]*length/2;
          coord[2] += coord[3]*length/2;
          coord[4] += length/2*sqrt(1+sqr(coord[1])+sqr(coord[3]));
        }
        
        /* compute energy kick */
        P     = *P_central*(1+coord[5]);
        beta_i = P/(gamma=sqrt(sqr(P)+1));
        t     = coord[4]/(c_mks*beta_i)-timeOffset;
        if ((dt = t-t0)<0)
            dt = 0;
        if  (!same_dgamma)
            dgamma = volt*sin(omega*t+phase)*(tau?sqrt(1-exp(-dt/tau)):1);

        /* apply energy kick */
        add_to_particle_energy(coord, t, *P_central, dgamma);

        if (length) {
          /* apply final drift */
          coord[0] += coord[1]*length/2;
          coord[2] += coord[3]*length/2;
          coord[4] += length/2.0*sqrt(1+sqr(coord[1])+sqr(coord[3]));
        }
      }

    if (rfca->change_p0)
        do_match_energy(part, np, P_central, 0);

    log_exit("simple_rf_cavity");
    return(np);
    }

void add_to_particle_energy(double *coord, double timeOfFlight, double Po, double dgamma)
{
  double gamma, gamma1, PRatio, P, P1, Pz1, Pz;

  P = Po*(1+coord[5]);                    /* old momentum */
  gamma1 = (gamma=sqrt(P*P+1)) + dgamma;  /* new gamma */
  if (gamma1<=1)
    gamma1 = 1+1e-7;
  P1 = sqrt(gamma1*gamma1-1);             /* new momentum */
  coord[5] = (P1-Po)/Po;                  

  /* adjust s for the new particle velocity */
  coord[4] = timeOfFlight*c_mks*P1/gamma1;

  /* adjust slopes so that Px and Py are conserved */
  Pz = P/sqrt(1+sqr(coord[1])+sqr(coord[3]));
  Pz1 = sqrt(Pz*Pz + gamma1*gamma1 - gamma*gamma);
  PRatio = Pz/Pz1;
  coord[1] *= PRatio;
  coord[3] *= PRatio;
}
