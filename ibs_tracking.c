/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* Intra-beam scattering element */

#include "mdb.h"
#include "track.h"
#include "zibs.h"
#include "SDDS.h"
#include "constants.h"

#define DEBUG 0

void setup_track_IBS(ELEMENT_LIST *element);
void inflateEmittance(double **coord, double Po, 
		      long offset, long np, double factor);
void SDDS_IBScatterSetup(SDDS_TABLE *SDDS_table, char *filename, long mode, long lines_per_row, char *contents, 
                         char *command_file, char *lattice_file, char *caller, long isRing);
void dump_IBScatter(SDDS_TABLE *SDDS_table, IBSCATTER *IBS, long pass);
void reset_IBS_output(ELEMENT_LIST *element);

#if DEBUG
static FILE *fpdeb = NULL;
#endif

void track_IBS(double **coord, long np, IBSCATTER *IBS, double Po, 
               ELEMENT_LIST *element, CHARGE *charge, long i_pass, RUN *run)
{
  double randomNumber;
  ONE_PLANE_PARAMETERS longitBeamParam, xBeamParam, yBeamParam;
  double beta0, beta1, p, vz;
  double RNSigma[3], RNSigmaCheck[3]={0,0,0};
  long ipart, icoord, ihcoord;
  static SDDS_TABLE outPage;
  static long isInit=0, doOut=0;

  if (!(IBS->s))
    setup_track_IBS(element);
#if DEBUG
  if (!fpdeb) {
    if (!(fpdeb=fopen("ibs.sdds", "w")))
      bomb("can't open ibs.sdds", NULL);
    fprintf(fpdeb, "SDDS1\n&column name=Pass type=long &end\n");
    fprintf(fpdeb, "&column name=xRate type=float &end\n");
    fprintf(fpdeb, "&column name=zRate type=float &end\n");
    fprintf(fpdeb, "&column name=xRateEmit type=float &end\n");
    fprintf(fpdeb, "&column name=zRateEmit type=float &end\n");
    fprintf(fpdeb, "&column name=xRN type=float &end\n");
    fprintf(fpdeb, "&column name=xRNCheck type=float &end\n");
    fprintf(fpdeb, "&column name=zRN type=float &end\n");
    fprintf(fpdeb, "&column name=zRNCheck type=float &end\n");
    fprintf(fpdeb, "&column name=emitx type=float &end\n");
    fprintf(fpdeb, "&column name=emitz type=float &end\n");
    fprintf(fpdeb, "&column name=sigmaz type=float &end\n");
    fprintf(fpdeb, "&column name=sigmaDelta type=float &end\n");
    fprintf(fpdeb, "&column name=charge type=float &end\n");
    fprintf(fpdeb, "&data mode=ascii no_row_counts=1 &end\n");
    fflush(fpdeb);
  }
#endif

  compute_transverse_parameters(&xBeamParam, coord, np, 0);
  IBS->emitx = IBS->emitx0 = xBeamParam.emittance;
  compute_transverse_parameters(&yBeamParam, coord, np, 2);
  IBS->emity = IBS->emity0 = yBeamParam.emittance;
  compute_longitudinal_parameters(&longitBeamParam, coord, np, Po);
  IBS->emitl = IBS->emitl0 = longitBeamParam.emittance;  /* units are seconds */
  vz = IBS->revolutionLength/IBS->dT;
  IBS->sigmaz = IBS->sigmaz0 = vz*sqrt(longitBeamParam.s11);
  IBS->sigmaDelta = IBS->sigmaDelta0 = sqrt(longitBeamParam.s22);

  if (charge)
    IBS->charge = charge->macroParticleCharge*np;   
  IBSRate (fabs(IBS->charge/e_mks), 
           IBS->coupling, IBS->elements, 1, IBS->verbosity, IBS->isRing,
           IBS->emitx0, IBS->emity0, IBS->sigmaDelta0, IBS->sigmaz0,
           IBS->s, IBS->pCentral, IBS->betax, IBS->alphax, IBS->betay,
           IBS->alphay, IBS->etax, IBS->etaxp, IBS->etay, IBS->etayp,
           IBS->xRateVsS, IBS->yRateVsS, IBS->zRateVsS,
           &IBS->xGrowthRate, &IBS->yGrowthRate, &IBS->zGrowthRate, 1);

  RNSigma[0] = RNSigma[1] = RNSigma[2] = 0;
  if (IBS->isRing && IBS->factor>0)
    IBS->dT *= IBS->factor;
  if (!IBS->smooth) {
    if (IBS->do_x && IBS->xGrowthRate>0)
      RNSigma[0] = sqrt(sqr(1 + IBS->dT * IBS->xGrowthRate)-1)*sqrt(xBeamParam.s22);
    if (IBS->do_y && IBS->yGrowthRate>0)
      RNSigma[1] = sqrt(sqr(1 + IBS->dT * IBS->yGrowthRate)-1)*sqrt(yBeamParam.s22);
    if (IBS->do_z && IBS->zGrowthRate>0)
      RNSigma[2] = sqrt(sqr(1 + IBS->dT * IBS->zGrowthRate)-1)*IBS->sigmaDelta0;
    for (icoord=1, ihcoord=0; icoord<6; icoord+=2, ihcoord++) {
      if (RNSigma[ihcoord]) {
        RNSigmaCheck[ihcoord] = 0;
        if (icoord!=5) {
          for (ipart=0; ipart<np; ipart++) {
            randomNumber = gauss_rn_lim(0.0, RNSigma[ihcoord], 3.0, random_2);
            coord[ipart][icoord] += randomNumber;
            RNSigmaCheck[ihcoord] += sqr(randomNumber);
          }
        } else {
          for (ipart=0; ipart<np; ipart++) {
            randomNumber = gauss_rn_lim(0.0, RNSigma[ihcoord], 3.0, random_2);
            p = Po*(1+coord[ipart][5]);
            beta0 = p/sqrt(p*p+1);
            coord[ipart][icoord] += randomNumber;
            p = Po*(1+coord[ipart][5]);
            beta1 = p/sqrt(p*p+1);
            coord[ipart][4] *= beta1/beta0;
            RNSigmaCheck[ihcoord] += sqr(randomNumber);
          }
        }
        RNSigmaCheck[ihcoord] = sqrt(RNSigmaCheck[ihcoord]/np);
      }
    }
  } else {
    /* inflate each emittance by the prescribed factor */
    inflateEmittance(coord, Po, 0, np, (1+IBS->dT*IBS->xGrowthRate));
    inflateEmittance(coord, Po, 2, np, (1+IBS->dT*IBS->yGrowthRate));
    inflateEmittance(coord, Po, 4, np, (1+IBS->dT*IBS->zGrowthRate));
  }

#if DEBUG
  fprintf(fpdeb, "%ld %le %le %le %le %le %le %le %le %le %le %le %le %le\n",
	  i_pass, 
	  IBS->xGrowthRate, IBS->zGrowthRate, 
	  IBS->xGrowthRate*IBS->emitx0, IBS->zGrowthRate*IBS->emitl0*vz,
	  RNSigma[0], RNSigmaCheck[0], RNSigma[2], RNSigmaCheck[2],
	  IBS->emitx0, IBS->emitl0*vz, IBS->sigmaz0, IBS->sigmaDelta0,
	  charge?
	  fabs(charge->macroParticleCharge*np):
	  fabs(IBS->charge));
  fflush(fpdeb);
#endif
  /* update beam emittance information after IBS scatter for IBSCATTER */
  if (!IBS->isRing) {
    compute_transverse_parameters(&xBeamParam, coord, np, 0);
    IBS->emitx = xBeamParam.emittance;
    compute_transverse_parameters(&yBeamParam, coord, np, 2);
    IBS->emity = yBeamParam.emittance;
    compute_longitudinal_parameters(&longitBeamParam, coord, np, Po);
    IBS->emitl = longitBeamParam.emittance;  /* units are seconds */
    IBS->sigmaz = vz*sqrt(longitBeamParam.s11);
    IBS->sigmaDelta = sqrt(longitBeamParam.s22);
  }

  if(IBS->filename) {
    if (!isInit) {
      SDDS_IBScatterSetup(&outPage, IBS->filename, SDDS_BINARY, 1, "IBS scatter growth rate output", 
                          run->runfile, run->lattice, "ibs_tracking", IBS->isRing);
      isInit = 1;
    }

    if ((int)i_pass/IBS->interval > doOut) {
      doOut++;
      reset_IBS_output(element);
   }
    if (IBS->output) {
      dump_IBScatter(&outPage, IBS, i_pass);
      IBS->output = 0;
    }
  }

  return;
}

void inflateEmittance(double **coord, double Po, 
		      long offset, long np, double factor)
{
  long i;
  double factorSqrt;
  if (!np)
    return;
  factorSqrt = sqrt(factor);
  if (offset==4) {
    /* longitudinal */
    double tc, dpc, *time, beta, p;
    time = tmalloc(sizeof(*time)*np);
    for (i=dpc=tc=0; i<np; i++) {
      dpc += coord[i][5];
      p = Po*(1+coord[i][5]);
      beta = p/sqrt(p*p+1);
      time[i] = coord[i][4]/beta;
      tc += time[i];
    }
    tc /= np;
    dpc /= np;
    for (i=0; i<np; i++) {
      time[i] = (time[i]-tc)*factorSqrt+tc;
      coord[i][5] = (coord[i][5]-dpc)*factorSqrt+dpc;
      p = Po*(1+coord[i][5]);
      beta = p/sqrt(p*p+1);
      coord[i][4] = time[i]*beta;
    }
    free(time);
  } else {
    double c[2];
    for (i=c[0]=c[1]=0; i<np; i++) {
      c[0] += coord[i][offset+0];
      c[1] += coord[i][offset+1];
    }
    c[0] /= np;
    c[1] /= np;
    for (i=0; i<np; i++) {
      coord[i][offset+0] = (coord[i][offset+0]-c[0])*factorSqrt+c[0];
      coord[i][offset+1] = (coord[i][offset+1]-c[1])*factorSqrt+c[1];
    }
  }
}

void reset_IBS_output(ELEMENT_LIST *element)
{
  ELEMENT_LIST *element0;
  IBSCATTER *IBS;

  element0 = element;
  while (element) {
    if (element->type==T_IBSCATTER) {
      IBS = (IBSCATTER*)element->p_elem; 
      IBS->output = 1;
    }
    element = element->succ;
  }
  element = element0;
}
/* Set twiss parameter arrays etc. */
void setup_track_IBS(ELEMENT_LIST *element)
{
  long count, nElements, i, isRing = 0;
  double startRingPos, finalPos;
  double s0, s1, dt, p0, gamma;
  ELEMENT_LIST *element0, *elementStartRing;
  IBSCATTER *IBS;
  
  if (!element->twiss) 
    bomb("Twiss parameters must be calculated for IBS tracking.", NULL);

  /* Find out start point of ring */
  element0 = elementStartRing = element;
  startRingPos = 0;

  count = 0;
  while (element) {
    if (element->type==T_RECIRC) {
      startRingPos = element->end_pos;
      elementStartRing = element->succ;
      break;
    }
    count++;
    element = element->succ;
  }
  if (elementStartRing!=element0) {
    element = elementStartRing;
  } else {
    element = element0;
    count = 0;
  }

  element0 = elementStartRing = element;
  nElements =0;
  s0 = s1 = startRingPos;
  dt =0;
  while (element) {
    s0 = s1;
    s1 = element->end_pos;
    if (s1 > s0 && element->pred) {
      p0 = (element->Pref_output + element->pred->Pref_output)/2.;
      gamma = sqrt(p0*p0+1);
      dt += (s1-s0)*gamma/p0/c_mks;
    }

    if (element->type==T_IBSCATTER) {
      IBS = (IBSCATTER*)element->p_elem; 
      IBS->revolutionLength = element->end_pos - startRingPos;
      IBS->dT = dt;
      dt = 0.;
      startRingPos = element->end_pos;
      IBS->elements = nElements;
      IBS->offset = count;
      IBS->output = 1;
      count = count + nElements +1;
      isRing = IBS->isRing;
      if (!(IBS->name = SDDS_Realloc(IBS->name, sizeof(*(IBS->name))*nElements)) ||
          !(IBS->s = SDDS_Realloc(IBS->s, sizeof(*(IBS->s))*nElements)) ||
          !(IBS->pCentral = SDDS_Realloc(IBS->pCentral, sizeof(*(IBS->pCentral))*nElements)) ||
          !(IBS->betax = SDDS_Realloc(IBS->betax, sizeof(*(IBS->betax))*nElements)) ||
          !(IBS->alphax = SDDS_Realloc(IBS->alphax, sizeof(*(IBS->alphax))*nElements)) ||
          !(IBS->betay = SDDS_Realloc(IBS->betay, sizeof(*(IBS->betay))*nElements)) ||
          !(IBS->alphay = SDDS_Realloc(IBS->alphay, sizeof(*(IBS->alphay))*nElements)) ||
          !(IBS->etax = SDDS_Realloc(IBS->etax, sizeof(*(IBS->etax))*nElements)) ||
          !(IBS->etaxp = SDDS_Realloc(IBS->etaxp, sizeof(*(IBS->etaxp))*nElements)) ||
          !(IBS->etay = SDDS_Realloc(IBS->etay, sizeof(*(IBS->etay))*nElements)) ||
          !(IBS->etayp = SDDS_Realloc(IBS->etayp, sizeof(*(IBS->etayp))*nElements)) ||
          !(IBS->xRateVsS = SDDS_Realloc(IBS->xRateVsS, sizeof(*(IBS->xRateVsS))*nElements)) ||
          !(IBS->yRateVsS = SDDS_Realloc(IBS->yRateVsS, sizeof(*(IBS->yRateVsS))*nElements)) ||
          !(IBS->zRateVsS = SDDS_Realloc(IBS->zRateVsS, sizeof(*(IBS->zRateVsS))*nElements)))
        bomb("memory allocation failure in setup_track_IBS", NULL);
      for (i=0; i<IBS->elements; i++) {
        cp_str(&IBS->name[i], elementStartRing->name);
        IBS->s[i] = elementStartRing->end_pos;
        IBS->pCentral[i] = elementStartRing->Pref_output;
        IBS->betax[i] = elementStartRing->twiss->betax;
        IBS->alphax[i] = elementStartRing->twiss->alphax;
        IBS->betay[i] = elementStartRing->twiss->betay;
        IBS->alphay[i] = elementStartRing->twiss->alphay;
        IBS->etax[i] = elementStartRing->twiss->etax;
        IBS->etaxp[i] = elementStartRing->twiss->etapx;
        IBS->etay[i] = elementStartRing->twiss->etay;
        IBS->etayp[i] = elementStartRing->twiss->etapy;
        elementStartRing = elementStartRing->succ;
      }
      nElements = -1;
      elementStartRing = element->succ;
    }
    nElements ++;
    finalPos = element->end_pos;
    element = element->succ;
  }
  if (isRing)
    if (finalPos != IBS->s[IBS->elements-1])
      bomb("You must have IBSCATTER at the end of the RING", NULL);
  element = element0;
}

#define IBSCATTER_RING_PARAMETERS 17
#define IBSCATTER_LINAC_PARAMETERS 23
static SDDS_DEFINITION ibscatter_print_parameter[IBSCATTER_LINAC_PARAMETERS] = {
  {"Charge", "&parameter name=Charge, type=double, units=\"nC\", description=\"Bunch charge in nC\" &end"},
  {"Particles", "&parameter name=Particles, type=double, description=\"Number of particles in bunch\" &end"},
  {"s", "&parameter name=s, type=double, units=\"m\", description=\"IBScatter element location in beamline\" &end"},
  {"Pass", "&parameter name=Pass, type=long, description=\"Pass number\" &end"},
  {"StartPoint", "&parameter name=StartPoint, type=long, description=\"IStart point in the beamline\" &end"},
  {"Elements", "&parameter name=Elements, type=long, description=\"Number of elements to integrate\" &end"},
  {"ds", "&parameter name=ds, type=double, units=\"m\", description=\"Integrated length for IBS rate\" &end"},  
  {"dt", "&parameter name=dt, type=double, units=\"m\", description=\"Integrated time for IBS scattering\" &end"},  
  {"xGrowthRate", "&parameter name=xGrowthRate, symbol=\"g$bIBS,x$n\", units=\"1/s\", type=double, description=\"Accumulated IBS emittance growth rate in the horizontal plane\" &end"},
  {"yGrowthRate", "&parameter name=yGrowthRate, symbol=\"g$bIBS,y$n\", units=\"1/s\", type=double, description=\"Accumulated IBS emittance growth rate in the vertical plane\" &end"},
  {"zGrowthRate", "&parameter name=zGrowthRate, symbol=\"g$bIBS,z$n\", units=\"1/s\", type=double, description=\"Accumulated IBS emittance growth rate in the longitudinal plane\" &end"},
  {"enx0", "&parameter name=enx0, symbol=\"$gge$r$bx$n\", units=\"m$be$nc $gp$rm\", type=double, description=\"Normalized initial horizontal emittance\" &end"},
  {"eny0", "&parameter name=eny0, symbol=\"$gge$r$by$n\", units=\"m$be$nc $gp$rm\", type=double, description=\"Normalized initial vertical emittance\" &end"},
  {"emitx0", "&parameter name=emitx0, symbol=\"$ge$r$bx,Input$n\", units=\"$gp$rm\", type=double, description=\"Initial horizontal emittance\" &end"},
  {"emity0", "&parameter name=emity0, symbol=\"$ge$r$by,Input$n\", units=\"$gp$rm\", type=double, description=\"Initial vertical emittance\" &end"},
  {"sigmaDelta0", "&parameter name=sigmaDelta0, symbol=\"$gs$r$bd,Input$n\", type=double, description=\"Initial momentum spread\" &end"},
  {"sigmaz0", "&parameter name=sigmaz0, symbol=\"$gs$r$bz,Input$n\", units=m, type=double, description=\"Initial bunch length\" &end"},
  {"enx", "&parameter name=enx, symbol=\"$gge$r$bx$n\", units=\"m$be$nc $gp$rm\", type=double, description=\"Normalized horizontal emittance with IBS\" &end"},
  {"eny", "&parameter name=eny, symbol=\"$gge$r$by$n\", units=\"m$be$nc $gp$rm\", type=double, description=\"Normalized vertical emittance with IBS\" &end"},
  {"emitx", "&parameter name=emitx, symbol=\"$ge$r$bx$n\", units=\"$gp$rm\", type=double, description=\"Horizontal emittance with IBS\" &end"},
  {"emity", "&parameter name=emity, symbol=\"$ge$r$by$n\", units=\"$gp$rm\", type=double, description=\"Vertical emittance with IBS\" &end"},
  {"sigmaDelta", "&parameter name=sigmaDelta, symbol=\"$gs$r$bd$n\", type=double, description=\"Momentum spread with IBS\" &end"},
  {"sigmaz", "&parameter name=sigmaz, symbol=\"$gs$r$bz$n\", units=m, type=double, description=\"Bunch length with IBS\" &end"},
};
#define IBSCATTER_COLUMNS 5
static SDDS_DEFINITION ibscatter_print_column[IBSCATTER_COLUMNS] = {
    {"ElementName", "&column name=ElementName, type=string, description=\"Element name\", format_string=%10s &end"},
    {"s", "&column name=s, units=m, type=double, description=\"Distance\" &end"},
    {"dIBSRatex", "&column name=dIBSRatex, units=\"1/(m s)\", type=double, description=\"Horizontal IBS Emittance Growth Rate\" &end"},
    {"dIBSRatey", "&column name=dIBSRatey, units=\"1/(m s)\", type=double, description=\"Vertical IBS Emittance Growth Rate\" &end"},
    {"dIBSRatez", "&column name=dIBSRatez, units=\"1/(m s)\", type=double, description=\"Longitudinal IBS Emittance Growth Rate\" &end"},
};

void SDDS_IBScatterSetup(SDDS_TABLE *SDDS_table, char *filename, long mode, long lines_per_row, char *contents, 
                         char *command_file, char *lattice_file, char *caller, long isRing)
{
    log_entry("SDDS_IBScatterSetup");

    if (isRing)
      SDDS_ElegantOutputSetup(SDDS_table, filename, mode, lines_per_row, contents, command_file, lattice_file,
                            ibscatter_print_parameter, IBSCATTER_RING_PARAMETERS, ibscatter_print_column, IBSCATTER_COLUMNS,
                            caller, SDDS_EOS_NEWFILE|SDDS_EOS_COMPLETE);
    else
      SDDS_ElegantOutputSetup(SDDS_table, filename, mode, lines_per_row, contents, command_file, lattice_file,
                            ibscatter_print_parameter, IBSCATTER_LINAC_PARAMETERS, ibscatter_print_column, IBSCATTER_COLUMNS,
                            caller, SDDS_EOS_NEWFILE|SDDS_EOS_COMPLETE);

    log_exit("SDDS_IBScatterSetup");
}

void dump_IBScatter(SDDS_TABLE *SDDS_table, IBSCATTER *IBS, long pass)
{
  long i;
  double gamma;

  log_entry("dump_IBScatter");

  if (!IBS->elements)
    return;

  if (!SDDS_StartTable(SDDS_table, IBS->elements)) {
    SDDS_SetError("Problem starting SDDS table (dump_IBScatter)");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }

  for (i=0; i<IBS->elements; i++) {
    if (!SDDS_SetRowValues(SDDS_table, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, i,
                           0, IBS->name[i], 1, IBS->s[i], 
                           2, IBS->xRateVsS[i], 3, IBS->yRateVsS[i], 4, IBS->zRateVsS[i], -1)) {
      SDDS_SetError("Problem setting SDDS row values (dump_IBScatter)");
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
  }
  gamma = sqrt(ipow(IBS->pCentral[IBS->elements-1], 2)+1.);
  if ((!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "Particles", fabs(IBS->charge/e_mks), NULL))||
      (!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "Charge", IBS->charge*1e9, NULL))||
      (!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "s", IBS->s[IBS->elements-1], NULL))||
      (!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "Pass", pass, NULL))||
      (!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "StartPoint", IBS->offset, NULL))||
      (!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "Elements", IBS->elements, NULL))||
      (!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "ds", IBS->revolutionLength, NULL))||
      (!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "dt", IBS->dT, NULL))||
      (!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "xGrowthRate", IBS->xGrowthRate, NULL))||
      (!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "yGrowthRate", IBS->yGrowthRate, NULL))||
      (!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "zGrowthRate", IBS->zGrowthRate, NULL))||
      (!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "enx0", IBS->emitx0*gamma, NULL))||
      (!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "eny0", IBS->emity0*gamma, NULL))||
      (!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "emitx0", IBS->emitx0, NULL))||
      (!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "emity0", IBS->emity0, NULL))||
      (!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "sigmaDelta0", IBS->sigmaDelta0, NULL))||
      (!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "sigmaz0", IBS->sigmaz0, NULL))){
    SDDS_SetError("Problem setting SDDS parameters (dump_IBScatter)");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }

  if (!IBS->isRing) {
    if ((!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "emitx", IBS->emitx, NULL))||
        (!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "emity", IBS->emity, NULL))||
        (!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "enx", IBS->emitx*gamma, NULL))||
        (!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "eny", IBS->emity*gamma, NULL))||
        (!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "sigmaDelta", IBS->sigmaDelta, NULL))||
        (!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "sigmaz", IBS->sigmaz, NULL))){
      SDDS_SetError("Problem setting SDDS parameters (dump_IBScatter)");
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
  }

  if (!SDDS_WriteTable(SDDS_table)) {
    SDDS_SetError("Problem writing SDDS table (dump_IBScatter)");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  } 

  SDDS_UpdatePage(SDDS_table, 0);
       
  SDDS_DoFSync(SDDS_table);
  log_exit("dump_IBScatter");    
}
