/* Copyright 1999 by Michael Borland and Argonne National Laboratory,
 * all rights reserved.
 */
/* Intra-beam scattering element */

#include "mdb.h"
#include "track.h"
#include "zibs.h"

#define DEBUG 0

void setup_track_IBS(IBSCATTER *IBS, ELEMENT_LIST *element);

#define DEBUG 1

#if DEBUG
static FILE *fpdeb = NULL;
static long pass = 0;
#endif

void track_IBS(double **coord, long np, IBSCATTER *IBS, double Po, 
               ELEMENT_LIST *element, RADIATION_INTEGRALS *radIntegrals0,
               CHARGE *charge)
{
  double gamma, emitx, emity, sigmaDelta, sigmaz;
  double xGrowthRate, yGrowthRate, zGrowthRate;
  double sigmax, sigmay, sigmat, emitl, dT;
  RADIATION_INTEGRALS radIntegrals;
  ONE_PLANE_PARAMETERS longitBeamParam, xBeamParam, yBeamParam;
  double vz;
  double RNSigma[3];
  long ipart, icoord, ihcoord;
  
  if (!(IBS->s))
    setup_track_IBS(IBS, element);
#if DEBUG
  if (!fpdeb) {
    if (!(fpdeb=fopen("ibs.sdds", "w")))
      bomb("can't open ibs.sdds", NULL);
    fprintf(fpdeb, "SDDS1\n&column name=Pass type=long &end\n");
    fprintf(fpdeb, "&column name=xRate type=double &end\n");
    fprintf(fpdeb, "&column name=yRate type=double &end\n");
    fprintf(fpdeb, "&column name=zRate type=double &end\n");
    fprintf(fpdeb, "&column name=xRN type=double &end\n");
    fprintf(fpdeb, "&column name=yRN type=double &end\n");
    fprintf(fpdeb, "&column name=zRN type=double &end\n");
    fprintf(fpdeb, "&column name=emitx type=double &end\n");
    fprintf(fpdeb, "&column name=emity type=double &end\n");
    fprintf(fpdeb, "&column name=emitz type=double &end\n");
    fprintf(fpdeb, "&column name=sigmaz type=double &end\n");
    fprintf(fpdeb, "&column name=sigmaDelta type=double &end\n");
    fprintf(fpdeb, "&column name=charge type=double &end\n");
    fprintf(fpdeb, "&data mode=ascii no_row_counts=1 &end\n");
    fflush(fpdeb);
  }
#endif

  memcpy((void*)&radIntegrals, radIntegrals0, sizeof(radIntegrals));
  computeRadiationIntegrals(&radIntegrals, Po, IBS->revolutionLength);
  compute_transverse_parameters(&xBeamParam, coord, np, 0);
  emitx = xBeamParam.emittance;
  compute_transverse_parameters(&yBeamParam, coord, np, 2);
  emity = yBeamParam.emittance;
  compute_longitudinal_parameters(&longitBeamParam, coord, np, Po);
  emitl = longitBeamParam.emittance;  /* units are seconds */
  gamma = sqrt(sqr(Po)+1);
  vz = c_mks*Po/gamma;
  sigmax = sqrt(xBeamParam.s11);
  sigmay = sqrt(yBeamParam.s11);
  sigmaz = vz*(sigmat = sqrt(longitBeamParam.s11));
  sigmaDelta = sqrt(longitBeamParam.s22);
  IBSGrowthRates(gamma, 
                 emitx, emity, sigmaDelta, 
                 sigmaz,
                 charge?fabs(charge->macroParticleCharge*np/e_mks):fabs(IBS->charge/e_mks),
                 radIntegrals.ex0, radIntegrals.sigmadelta,
                 2.0/radIntegrals.taux, 2.0/radIntegrals.taudelta,
                 IBS->coupling, 
                 IBS->s, IBS->betax, IBS->alphax, IBS->betay, IBS->alphay,
                 IBS->etax, IBS->etaxp, IBS->elements, 1, 0, 
                 &xGrowthRate, &yGrowthRate, &zGrowthRate
                 );
  RNSigma[0] = RNSigma[1] = RNSigma[2] = 0;
  dT = IBS->revolutionLength/vz;
  if (xGrowthRate>0)
    RNSigma[0] = sqrt((sqr(1+dT*xGrowthRate)-1))*sqrt(xBeamParam.s22);
  if (yGrowthRate>0)
    RNSigma[1] = sqrt((sqr(1+dT*yGrowthRate)-1))*sqrt(yBeamParam.s22);
  if (zGrowthRate>0)
    RNSigma[2] = sqrt((sqr(1+dT*zGrowthRate)-1))*sigmaDelta;
#if DEBUG
  fprintf(fpdeb, "%ld %le %le %le %le %le %le %le %le %le %le %le %le\n",
          pass++, xGrowthRate, yGrowthRate, zGrowthRate, 
	  RNSigma[0], RNSigma[1], RNSigma[2],
	  emitx, emity, sigmaz*sigmaDelta, sigmaz, sigmaDelta,
	  charge?
	  fabs(charge->macroParticleCharge*np):
	  fabs(IBS->charge));
  fflush(fpdeb);
#endif

  for (icoord=1, ihcoord=0; icoord<6; icoord+=2, ihcoord++) {
    if (RNSigma[ihcoord])
      for (ipart=0; ipart<np; ipart++) {
        coord[ipart][icoord] += RNSigma[ihcoord]*gauss_rn(0, random_2);
      }
  }
}


/* Set twiss parameter arrays etc. */
void setup_track_IBS(IBSCATTER *IBS, ELEMENT_LIST *element)
{
  long iElement, nElements, offset, count;
  double startRingPos, endRingPos;
  ELEMENT_LIST *element0, *elementStartRing;
  
  if (!element->twiss) 
    bomb("Twiss parameters must be calculated for IBS tracking.", NULL);

  element0 = elementStartRing = element;
  offset = count = 0;
  startRingPos = endRingPos = 0;
  while (element) {
    if (element->type==T_RECIRC) {
      offset = count;
      startRingPos = element->end_pos;
      elementStartRing = element->succ;
    }
    count ++;
    endRingPos = element->end_pos;
    element = element->succ;
  }
  if (elementStartRing!=element) {
    element = elementStartRing;
  } else 
    element = element0;
  IBS->revolutionLength = endRingPos-startRingPos;
  
  if (!(IBS->elements = nElements = count-offset))
    bomb("No elements in ring for IBS tracking.", NULL);

  if (!(IBS->s = SDDS_Realloc(IBS->s, sizeof(*(IBS->s))*nElements)) ||
      !(IBS->betax = SDDS_Realloc(IBS->betax, sizeof(*(IBS->betax))*nElements)) ||
      !(IBS->alphax = SDDS_Realloc(IBS->alphax, sizeof(*(IBS->alphax))*nElements)) ||
      !(IBS->betay = SDDS_Realloc(IBS->betay, sizeof(*(IBS->betay))*nElements)) ||
      !(IBS->alphay = SDDS_Realloc(IBS->alphay, sizeof(*(IBS->alphay))*nElements)) ||
      !(IBS->etax = SDDS_Realloc(IBS->etax, sizeof(*(IBS->etax))*nElements)) ||
      !(IBS->etaxp = SDDS_Realloc(IBS->etaxp, sizeof(*(IBS->etaxp))*nElements)))
    bomb("memory allocation failure in setup_track_IBS", NULL);

  iElement = 0;
  while (element) {
    IBS->s[iElement] = element->end_pos;
    IBS->betax[iElement] = element->twiss->betax;
    IBS->alphax[iElement] = element->twiss->alphax;
    IBS->betay[iElement] = element->twiss->betay;
    IBS->alphay[iElement] = element->twiss->alphay;
    IBS->etax[iElement] = element->twiss->etax;
    IBS->etaxp[iElement] = element->twiss->etapx;
    iElement++;
    element = element->succ;
  }
}


