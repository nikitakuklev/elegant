/* Copyright 1999 by Michael Borland and Argonne National Laboratory,
 * all rights reserved.
 */
/*
 * $Log: not supported by cvs2svn $
 * Revision 1.7  2000/05/10 01:52:15  borland
 * Now correctly compute matrices for alpha magnet and stray field in the presence
 * of changes in central momentum.
 * SASE FEL now gives slice average plus non-slice computation when slices are
 * requested.
 *
 * Revision 1.6  2000/04/21 20:52:47  soliday
 * Added include fdlibm.h for Bessel function with Borland C.
 *
 * Revision 1.5  2000/04/20 20:22:35  borland
 * Added ability to do computations for slices.
 *
 * Revision 1.4  2000/01/25 19:49:16  borland
 * Removed unnecessary array and inserted free statement for another.
 * Now uses average momentum rather than central momentum.
 *
 * Revision 1.3  1999/10/12 21:50:00  borland
 * All printouts now go to the stdout rather than stderr.  fflush statements,
 * some unnecessary, were added in a mostly automated fashion.
 *
 * Revision 1.2  1999/08/05 15:40:23  soliday
 * Added WIN32 and Linux support
 *
 * Revision 1.1  1999/07/01 19:19:57  borland
 * First versions in repository.
 *
 */
#include "mdb.h"
#include "track.h"
#include "sasefel.h"
#if defined(__BORLANDC__)
#include <fdlibm.h>
#endif

long DefineSASEParameters(SASEFEL_OUTPUT *sasefelOutput, long slice);

void setupSASEFELAtEnd(NAMELIST_TEXT *nltext, RUN *run, OUTPUT_FILES *output_data)
{
  SASEFEL_OUTPUT *sasefelOutput;
  SDDS_DATASET *SDDSout;

  /* process namelist text */
  process_namelist(&sasefel, nltext);
  print_namelist(stdout, &sasefel);

  if (beta<0)
    bomb("beta < 0", NULL);
  if (undulator_K<=0)
    bomb("undulator_K <= 0", NULL);
  if (undulator_period<=0)
    bomb("undulator_period <= 0", NULL);
  if (n_slices<0 || slice_fraction<0 ||
      (n_slices==0 && slice_fraction>0) ||
      (n_slices>0 && slice_fraction<=0) ||
      n_slices*slice_fraction>1)
    bomb("invalid slice parameters", NULL);
  sasefelOutput = &(output_data->sasefel);
  if (sasefelOutput->active && sasefelOutput->filename) {
    SDDS_Terminate(&(sasefelOutput->SDDSout));
    SDDS_ClearErrors();
    if (sasefelOutput->betaToUse) free(sasefelOutput->betaToUse);
    if (sasefelOutput->charge) free(sasefelOutput->charge);
    if (sasefelOutput->pCentral) free(sasefelOutput->pCentral);
    if (sasefelOutput->rmsBunchLength) free(sasefelOutput->rmsBunchLength);
    if (sasefelOutput->Sdelta) free(sasefelOutput->Sdelta);
    if (sasefelOutput->emit) free(sasefelOutput->emit);
    if (sasefelOutput->lightWavelength) free(sasefelOutput->lightWavelength);
    if (sasefelOutput->saturationLength) free(sasefelOutput->saturationLength);
    if (sasefelOutput->gainLength) free(sasefelOutput->gainLength);
    if (sasefelOutput->noisePower) free(sasefelOutput->noisePower);
    if (sasefelOutput->saturationPower) free(sasefelOutput->saturationPower);
    if (sasefelOutput->PierceParameter) free(sasefelOutput->PierceParameter);
    if (sasefelOutput->etaDiffraction) free(sasefelOutput->etaDiffraction);
    if (sasefelOutput->etaEmittance) free(sasefelOutput->etaEmittance);
    if (sasefelOutput->etaEnergySpread) free(sasefelOutput->etaEnergySpread);
    if (sasefelOutput->betaToUseIndex) free(sasefelOutput->betaToUseIndex);
    if (sasefelOutput->chargeIndex) free(sasefelOutput->chargeIndex);
    if (sasefelOutput->pCentralIndex) free(sasefelOutput->pCentralIndex);
    if (sasefelOutput->rmsBunchLengthIndex) free(sasefelOutput->rmsBunchLengthIndex);
    if (sasefelOutput->SdeltaIndex) free(sasefelOutput->SdeltaIndex);
    if (sasefelOutput->emitIndex) free(sasefelOutput->emitIndex);
    if (sasefelOutput->lightWavelengthIndex) free(sasefelOutput->lightWavelengthIndex);
    if (sasefelOutput->saturationLengthIndex) free(sasefelOutput->saturationLengthIndex);
    if (sasefelOutput->gainLengthIndex) free(sasefelOutput->gainLengthIndex);
    if (sasefelOutput->noisePowerIndex) free(sasefelOutput->noisePowerIndex);
    if (sasefelOutput->saturationPowerIndex) free(sasefelOutput->saturationPowerIndex);
    if (sasefelOutput->PierceParameterIndex) free(sasefelOutput->PierceParameterIndex);
    if (sasefelOutput->etaDiffractionIndex) free(sasefelOutput->etaDiffractionIndex);
    if (sasefelOutput->etaEmittanceIndex) free(sasefelOutput->etaEmittanceIndex);
    if (sasefelOutput->etaEnergySpreadIndex) free(sasefelOutput->etaEnergySpreadIndex);
    if (sasefelOutput->sliceFound) free(sasefelOutput->sliceFound);
  }
  
  sasefelOutput->active = 1;
  sasefelOutput->beta = beta;
  sasefelOutput->undulatorK = undulator_K;
  sasefelOutput->undulatorPeriod = undulator_period;
  sasefelOutput->nSlices = n_slices;
  sasefelOutput->sliceFraction = slice_fraction;

  if (!(sasefelOutput->betaToUse = malloc(sizeof(*(sasefelOutput->betaToUse))*(n_slices+2))) ||
      !(sasefelOutput->charge = malloc(sizeof(*(sasefelOutput->charge))*(n_slices+2))) ||
      !(sasefelOutput->pCentral = malloc(sizeof(*(sasefelOutput->pCentral))*(n_slices+2))) ||
      !(sasefelOutput->rmsBunchLength = malloc(sizeof(*(sasefelOutput->rmsBunchLength))*(n_slices+2))) ||
      !(sasefelOutput->Sdelta = malloc(sizeof(*(sasefelOutput->Sdelta))*(n_slices+2))) ||
      !(sasefelOutput->emit = malloc(sizeof(*(sasefelOutput->emit))*(n_slices+2))) ||
      !(sasefelOutput->lightWavelength = malloc(sizeof(*(sasefelOutput->lightWavelength))*(n_slices+2))) ||
      !(sasefelOutput->saturationLength = malloc(sizeof(*(sasefelOutput->saturationLength))*(n_slices+2))) ||
      !(sasefelOutput->gainLength = malloc(sizeof(*(sasefelOutput->gainLength))*(n_slices+2))) ||
      !(sasefelOutput->noisePower = malloc(sizeof(*(sasefelOutput->noisePower))*(n_slices+2))) ||
      !(sasefelOutput->saturationPower = malloc(sizeof(*(sasefelOutput->saturationPower))*(n_slices+2))) ||
      !(sasefelOutput->PierceParameter = malloc(sizeof(*(sasefelOutput->PierceParameter))*(n_slices+2))) ||
      !(sasefelOutput->etaDiffraction = malloc(sizeof(*(sasefelOutput->etaDiffraction))*(n_slices+2))) ||
      !(sasefelOutput->etaEmittance = malloc(sizeof(*(sasefelOutput->etaEmittance))*(n_slices+2))) ||
      !(sasefelOutput->etaEnergySpread = malloc(sizeof(*(sasefelOutput->etaEnergySpread))*(n_slices+2)))) 
    bomb("memory allocation failure (setupSASEFELAtEnd)", NULL);

  if (!(sasefelOutput->betaToUseIndex = malloc(sizeof(*(sasefelOutput->betaToUseIndex))*(n_slices+2))) ||
      !(sasefelOutput->chargeIndex = malloc(sizeof(*(sasefelOutput->chargeIndex))*(n_slices+2))) ||
      !(sasefelOutput->pCentralIndex = malloc(sizeof(*(sasefelOutput->pCentralIndex))*(n_slices+2))) ||
      !(sasefelOutput->rmsBunchLengthIndex = malloc(sizeof(*(sasefelOutput->rmsBunchLengthIndex))*(n_slices+2))) ||
      !(sasefelOutput->SdeltaIndex = malloc(sizeof(*(sasefelOutput->SdeltaIndex))*(n_slices+2))) ||
      !(sasefelOutput->emitIndex = malloc(sizeof(*(sasefelOutput->emitIndex))*(n_slices+2))) ||
      !(sasefelOutput->lightWavelengthIndex = malloc(sizeof(*(sasefelOutput->lightWavelengthIndex))*(n_slices+2))) ||
      !(sasefelOutput->saturationLengthIndex = malloc(sizeof(*(sasefelOutput->saturationLengthIndex))*(n_slices+2))) ||
      !(sasefelOutput->gainLengthIndex = malloc(sizeof(*(sasefelOutput->gainLengthIndex))*(n_slices+2))) ||
      !(sasefelOutput->noisePowerIndex = malloc(sizeof(*(sasefelOutput->noisePowerIndex))*(n_slices+2))) ||
      !(sasefelOutput->saturationPowerIndex = malloc(sizeof(*(sasefelOutput->saturationPowerIndex))*(n_slices+2))) ||
      !(sasefelOutput->PierceParameterIndex = malloc(sizeof(*(sasefelOutput->PierceParameterIndex))*(n_slices+2))) ||
      !(sasefelOutput->etaDiffractionIndex = malloc(sizeof(*(sasefelOutput->etaDiffractionIndex))*(n_slices+2))) ||
      !(sasefelOutput->etaEmittanceIndex = malloc(sizeof(*(sasefelOutput->etaEmittanceIndex))*(n_slices+2))) ||
      !(sasefelOutput->etaEnergySpreadIndex = malloc(sizeof(*(sasefelOutput->etaEnergySpreadIndex))*(n_slices+2))) ||
      !(sasefelOutput->sliceFound = malloc(sizeof(*(sasefelOutput->sliceFound))*(n_slices+2)))) 
    bomb("memory allocation failure (setupSASEFELAtEnd)", NULL);

  if (output) {

    sasefelOutput->filename = compose_filename(output, run->rootname);
    SDDSout = &(sasefelOutput->SDDSout);
    if (!SDDS_InitializeOutput(SDDSout, SDDS_BINARY, 0, NULL, NULL, sasefelOutput->filename) ||
        !SDDS_DefineSimpleParameter(SDDSout, "Step", NULL, SDDS_LONG) ||
        !SDDS_DefineSimpleParameter(SDDSout, "undulatorK", NULL, SDDS_DOUBLE) ||
        !SDDS_DefineSimpleParameter(SDDSout, "undulatorPeriod", "m", SDDS_DOUBLE)) {
      fprintf(stdout, "Unable define SDDS parameter for file %s\n", sasefelOutput->filename);
      fflush(stdout);
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
    if (n_slices<=1) {
      if (!DefineSASEParameters(sasefelOutput, 0)) {
        fprintf(stdout, "Unable define SDDS parameters for file %s\n", sasefelOutput->filename);
        fflush(stdout);
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      }
    } else {
      long slice;
      /* slice 0 is the nominal (no slicing)
       * slice N+1 is the average over the slices 
       */
      for (slice=0; slice<=n_slices+1; slice++) {
        if (!DefineSASEParameters(sasefelOutput, slice)) {
          fprintf(stdout, "Unable define SDDS parameters for file %s\n", sasefelOutput->filename);
          fflush(stdout);
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
      }
    }
    if (!SDDS_WriteLayout(SDDSout))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
}

long DefineSASEParameters(SASEFEL_OUTPUT *sasefelOutput, long slice)
{
  SDDS_DATASET *SDDSout;
  char buffer[100], sliceNumString[20];
  
  SDDSout = &(sasefelOutput->SDDSout);
  if (slice && sasefelOutput->nSlices>1) {
    if (slice<=sasefelOutput->nSlices) {
      sprintf(sliceNumString, "Slice%02ld", slice);
    } else {
      /* "slice" N+1 is the average over all slices */
      sprintf(sliceNumString, "Ave");
    }
  }
  else
    sliceNumString[0] = 0;

  sprintf(buffer, "beta%s", sliceNumString);
  if ((sasefelOutput->betaToUseIndex[slice] = 
       SDDS_DefineParameter(SDDSout, buffer, NULL, "m", NULL, NULL, SDDS_DOUBLE, NULL))<0)
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  sprintf(buffer, "charge%s", sliceNumString);
  if ((sasefelOutput->chargeIndex[slice] = 
       SDDS_DefineParameter(SDDSout, buffer, NULL, "C", NULL, NULL, SDDS_DOUBLE, NULL))<0)
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  sprintf(buffer, "rmsBunchLength%s", sliceNumString);
  if ((sasefelOutput->rmsBunchLengthIndex[slice] = 
       SDDS_DefineParameter(SDDSout, buffer, NULL, "s", NULL, NULL, SDDS_DOUBLE, NULL))<0)
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  sprintf(buffer, "Sdelta%s", sliceNumString);
  if ((sasefelOutput->SdeltaIndex[slice] = 
       SDDS_DefineParameter(SDDSout, buffer, NULL, NULL, NULL, NULL, SDDS_DOUBLE, NULL))<0)
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  sprintf(buffer, "emit%s", sliceNumString);
  if ((sasefelOutput->emitIndex[slice] = 
       SDDS_DefineParameter(SDDSout, buffer, NULL, "m", NULL, NULL, SDDS_DOUBLE, NULL))<0)
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  sprintf(buffer, "pCentral%s", sliceNumString);
  if ((sasefelOutput->pCentralIndex[slice] = 
       SDDS_DefineParameter(SDDSout, buffer, NULL, "m$be$nc", NULL, NULL, SDDS_DOUBLE, NULL))<0)
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  sprintf(buffer, "lightWavelength%s", sliceNumString);
  if ((sasefelOutput->lightWavelengthIndex[slice] = 
       SDDS_DefineParameter(SDDSout, buffer, NULL, "m", NULL, NULL, SDDS_DOUBLE, NULL))<0)
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  sprintf(buffer, "gainLength%s", sliceNumString);
  if ((sasefelOutput->gainLengthIndex[slice] = 
       SDDS_DefineParameter(SDDSout, buffer, NULL, "m", NULL, NULL, SDDS_DOUBLE, NULL))<0)
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  sprintf(buffer, "noisePower%s", sliceNumString);
  if ((sasefelOutput->noisePowerIndex[slice] = 
       SDDS_DefineParameter(SDDSout, buffer, NULL, "W", NULL, NULL, SDDS_DOUBLE, NULL))<0)
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  sprintf(buffer, "saturationPower%s", sliceNumString);
  if ((sasefelOutput->saturationPowerIndex[slice] = 
       SDDS_DefineParameter(SDDSout, buffer, NULL, "W", NULL, NULL, SDDS_DOUBLE, NULL))<0)
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  sprintf(buffer, "saturationLength%s", sliceNumString);
  if ((sasefelOutput->saturationLengthIndex[slice] = 
       SDDS_DefineParameter(SDDSout, buffer, NULL, "m", NULL, NULL, SDDS_DOUBLE, NULL))<0)
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  sprintf(buffer, "PierceParameter%s", sliceNumString);
  if ((sasefelOutput->PierceParameterIndex[slice] = 
       SDDS_DefineParameter(SDDSout, buffer, NULL, NULL, NULL, NULL, SDDS_DOUBLE, NULL))<0)
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  sprintf(buffer, "etaDiffraction%s", sliceNumString);
  if ((sasefelOutput->etaDiffractionIndex[slice] = 
       SDDS_DefineParameter(SDDSout, buffer, NULL, NULL, NULL, NULL, SDDS_DOUBLE, NULL))<0)
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  sprintf(buffer, "etaEmittance%s", sliceNumString);
  if ((sasefelOutput->etaEmittanceIndex[slice] = 
       SDDS_DefineParameter(SDDSout, buffer, NULL, NULL, NULL, NULL, SDDS_DOUBLE, NULL))<0)
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  sprintf(buffer, "etaEnergySpread%s", sliceNumString);
  if ((sasefelOutput->etaEnergySpreadIndex[slice] = 
       SDDS_DefineParameter(SDDSout, buffer, NULL, NULL, NULL, NULL, SDDS_DOUBLE, NULL))<0)
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
}


void doSASEFELAtEndOutput(SASEFEL_OUTPUT *sasefelOutput, long step)
{
  SDDS_DATASET *SDDSout;
  long slice, i;
  char buffer[100];
  
  if (!sasefelOutput || !sasefelOutput->active) 
    SDDS_Bomb("doSASEFELAtEndOutput called without proper setup!");
  if (!sasefelOutput->filename)
    return;
  SDDSout = &(sasefelOutput->SDDSout);
  if (!SDDS_StartPage(SDDSout, 0) ||
      !SDDS_SetParameters(SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,
                          "Step", step,
                          "undulatorK", sasefelOutput->undulatorK,
                          "undulatorPeriod", sasefelOutput->undulatorPeriod,
                          NULL)) {
    fprintf(stdout, "Unable write data to file %s\n", sasefelOutput->filename);
    fflush(stdout);
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  for (slice=0; slice<=sasefelOutput->nSlices+1; slice++) {
    if (!SDDS_SetParameters(SDDSout, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE,
                            sasefelOutput->betaToUseIndex[slice], sasefelOutput->betaToUse[slice], 
                            sasefelOutput->chargeIndex[slice], sasefelOutput->charge[slice], 
                            sasefelOutput->pCentralIndex[slice], sasefelOutput->pCentral[slice], 
                            sasefelOutput->rmsBunchLengthIndex[slice], sasefelOutput->rmsBunchLength[slice], 
                            sasefelOutput->SdeltaIndex[slice], sasefelOutput->Sdelta[slice], 
                            sasefelOutput->emitIndex[slice], sasefelOutput->emit[slice], 
                            sasefelOutput->lightWavelengthIndex[slice], sasefelOutput->lightWavelength[slice], 
                            sasefelOutput->saturationLengthIndex[slice], sasefelOutput->saturationLength[slice], 
                            sasefelOutput->gainLengthIndex[slice], sasefelOutput->gainLength[slice], 
                            sasefelOutput->noisePowerIndex[slice], sasefelOutput->noisePower[slice], 
                            sasefelOutput->saturationPowerIndex[slice], sasefelOutput->saturationPower[slice], 
                            sasefelOutput->PierceParameterIndex[slice], sasefelOutput->PierceParameter[slice], 
                            sasefelOutput->etaDiffractionIndex[slice], sasefelOutput->etaDiffraction[slice], 
                            sasefelOutput->etaEmittanceIndex[slice], sasefelOutput->etaEmittance[slice], 
                            sasefelOutput->etaEnergySpreadIndex[slice], sasefelOutput->etaEnergySpread[slice], 
                            -1))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  if (!SDDS_WritePage(SDDSout))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
}


void computeSASEFELAtEnd(SASEFEL_OUTPUT *sasefelOutput, double **particle, long particles, 
                         double Po, double charge)
{
  double emitx, emity;
  double *time, deltaAve, deltaRMS, tRMS;
  double S11, S12, S22, S33, S34, S44;
  long i, slice, nSlices;
  double xLimit[2], percentLevel[2];
  long count, slicesFound=0, j;
  double aveCoord[6], rmsCoord[6];
  
  if (!particles) {
    fprintf(stdout, "no particles left---can't compute FEL parameters");
    fflush(stdout);
    /* fill in some dummy values */
    for (slice=0; slice<sasefelOutput->nSlices+1; slice++) {
      sasefelOutput->betaToUse[slice] = sasefelOutput->charge[slice] = DBL_MAX;
      sasefelOutput->pCentral[slice] = sasefelOutput->rmsBunchLength[slice] = DBL_MAX;
      sasefelOutput->Sdelta[slice] = sasefelOutput->emit[slice] = DBL_MAX;
      sasefelOutput->lightWavelength[slice] = sasefelOutput->saturationLength[slice] = DBL_MAX;
      sasefelOutput->gainLength[slice] = sasefelOutput->noisePower[slice] = DBL_MAX;
      sasefelOutput->saturationPower[slice] = sasefelOutput->PierceParameter[slice] = DBL_MAX;
      sasefelOutput->etaDiffraction[slice] = sasefelOutput->etaEmittance[slice] = DBL_MAX;
      sasefelOutput->etaEnergySpread[slice] = DBL_MAX;
    }
    return;
  }
  if (!(time=malloc(sizeof(*time)*particles)))
    SDDS_Bomb("memory allocation failure (computeSASEFELAtEnd)");
  computeTimeCoordinates(time, Po, particle, particles);

  /* compute normal values (over entire beam) */
  if (sasefelOutput->nSlices==0) {
    percentLevel[0] = 10;
    percentLevel[1] = 90;
  } else {
    percentLevel[0] = 50 - 100*sasefelOutput->sliceFraction/2.0;
    percentLevel[1] = 50 + 100*sasefelOutput->sliceFraction/2.0;
  }
  
  /* find center of energy distribution */
  for (i=deltaAve=0; i<particles; i++)
    deltaAve += particle[i][5];
  deltaAve /= particles;
  /* compute rms energy spread */
  for (i=deltaRMS=0; i<particles; i++)
    deltaRMS += sqr(particle[i][5]-deltaAve);
  sasefelOutput->Sdelta[0] = deltaRMS = sqrt(deltaRMS/particles);
  
  /* compute rms-equivalent time value so that Q/(sqrt(2*PI)*tRMS) is
   * a good estimate of peak current.  I use the 10% and 90% points of
   * the distribution to compute peak current, then get equivalent tRMS.
   */
  compute_percentiles(xLimit, percentLevel, 2, time, particles);
  sasefelOutput->rmsBunchLength[0] = tRMS 
    = (xLimit[1] - xLimit[0])/(0.8*sqrt(2*PI));
  
  emitx = rms_emittance(particle, 0, 1, particles, &S11, NULL, NULL);
  emity = rms_emittance(particle, 2, 3, particles, &S33, NULL, NULL);
  sasefelOutput->emit[0] = sqrt(emitx*emity);
  sasefelOutput->pCentral[0] = Po*(1+deltaAve);
  sasefelOutput->charge[0] = charge;  
  
  if (sasefelOutput->beta==0)
    sasefelOutput->betaToUse[0] = sqrt(S11*S33/(emitx*emity));
  else
    sasefelOutput->betaToUse[0] = sasefelOutput->beta;
  
  ComputeSASEFELParameters(&sasefelOutput->lightWavelength[0], &sasefelOutput->saturationLength[0], 
                           &sasefelOutput->gainLength[0],
                           &sasefelOutput->noisePower[0], &sasefelOutput->saturationPower[0], 
                           &sasefelOutput->PierceParameter[0],
                           &sasefelOutput->etaDiffraction[0], &sasefelOutput->etaEmittance[0], 
                           &sasefelOutput->etaEnergySpread[0],
                           charge, tRMS, sasefelOutput->undulatorPeriod,
                           sasefelOutput->undulatorK, 
                           sasefelOutput->betaToUse[0], sasefelOutput->emit[0], 
                           sasefelOutput->Sdelta[0], sasefelOutput->pCentral[0],
                           1);
  if (sasefelOutput->nSlices>1) {
    /* compute values for each slice, plus average */
    nSlices = sasefelOutput->nSlices;
    
    sasefelOutput->betaToUse[nSlices+1] = sasefelOutput->charge[nSlices+1] = 0;
    sasefelOutput->pCentral[nSlices+1] = sasefelOutput->rmsBunchLength[nSlices+1] = 0;
    sasefelOutput->Sdelta[nSlices+1] = sasefelOutput->emit[nSlices+1] = 0;
    sasefelOutput->lightWavelength[nSlices+1] = sasefelOutput->saturationLength[nSlices+1] = 0;
    sasefelOutput->gainLength[nSlices+1] = sasefelOutput->noisePower[nSlices+1] = 0;
    sasefelOutput->saturationPower[nSlices+1] = sasefelOutput->PierceParameter[nSlices+1] = 0;
    sasefelOutput->etaDiffraction[nSlices+1] = sasefelOutput->etaEmittance[nSlices+1] = 0;
    sasefelOutput->etaEnergySpread[nSlices+1] = 0;
  
    for (slice=1; slice<=sasefelOutput->nSlices; slice++) {
      /* find boundaries of slice in time */
      percentLevel[0] = 100*(0.5-sasefelOutput->nSlices*sasefelOutput->sliceFraction/2.0 + 
                             (slice-1)*sasefelOutput->sliceFraction);
      if (percentLevel[0]<0)
        percentLevel[0] = 0;
      percentLevel[1] = percentLevel[0] + 100*sasefelOutput->sliceFraction;
      if (percentLevel[1]>100)
        percentLevel[1] = 100;
      
      /* compute rms-equivalent time value so that Q/(sqrt(2*PI)*tRMS) is
       * the average current in the slice
       */
      compute_percentiles(xLimit, percentLevel, 2, time, particles);
      sasefelOutput->rmsBunchLength[slice] = (xLimit[1] - xLimit[0])/sqrt(2*PI);
      
      /* find center of energy distribution */
      for (j=0; j<6; j++)
        aveCoord[j] = 0;
      for (i=count=0; i<particles; i++) {
        if (time[i]>xLimit[0] && time[i]<xLimit[1]) {
          count++;
          for (j=0; j<6; j++)
            aveCoord[j] += particle[i][j];
        }
      }
      if (count<2) {
        /* fill in some dummy values */
        sasefelOutput->sliceFound[slice] = 0;
        sasefelOutput->betaToUse[slice] = sasefelOutput->charge[slice] = DBL_MAX;
        sasefelOutput->pCentral[slice] = sasefelOutput->rmsBunchLength[slice] = DBL_MAX;
        sasefelOutput->Sdelta[slice] = sasefelOutput->emit[slice] = DBL_MAX;
        sasefelOutput->lightWavelength[slice] = sasefelOutput->saturationLength[slice] = DBL_MAX;
        sasefelOutput->gainLength[slice] = sasefelOutput->noisePower[slice] = DBL_MAX;
        sasefelOutput->saturationPower[slice] = sasefelOutput->PierceParameter[slice] = DBL_MAX;
        sasefelOutput->etaDiffraction[slice] = sasefelOutput->etaEmittance[slice] = DBL_MAX;
        sasefelOutput->etaEnergySpread[slice] = DBL_MAX;
      }
      slicesFound++;
      sasefelOutput->sliceFound[slice] = 1;
      for (j=0; j<6; j++) {
        aveCoord[j] /= count;
        rmsCoord[j] = 0;
      }
      S12 = S34 = 0;
      
      /* compute rms energy spread and transverse moments */
      for (i=deltaRMS=0; i<particles; i++)
        if (time[i]>xLimit[0] && time[i]<xLimit[1]) {
          for (j=0; j<6; j++)
            rmsCoord[j] += sqr(particle[i][j]-aveCoord[j]);
          S12 += (particle[i][0]-aveCoord[0])*(particle[i][1]-aveCoord[1]);
          S34 += (particle[i][2]-aveCoord[2])*(particle[i][3]-aveCoord[3]);
        }
      for (j=0; j<6; j++)
        rmsCoord[j] = sqrt(rmsCoord[j]/count);
      S12 /= count;
      S34 /= count;
      S11 = sqr(rmsCoord[0]);
      S22 = sqr(rmsCoord[1]);
      S33 = sqr(rmsCoord[2]);
      S44 = sqr(rmsCoord[3]);
      
      sasefelOutput->Sdelta[slice] = rmsCoord[5];
      emitx = sqrt(S11*S22-sqr(S12));
      emity = sqrt(S33*S44-sqr(S34));
      sasefelOutput->emit[slice] = sqrt(emitx*emity);
      sasefelOutput->pCentral[slice] = Po*(1+aveCoord[5]);
      sasefelOutput->charge[slice] = charge*sasefelOutput->sliceFraction;  

      if (sasefelOutput->beta==0)
        sasefelOutput->betaToUse[slice] = sqrt(S11*S33/(emitx*emity));
      else
        sasefelOutput->betaToUse[slice] = sasefelOutput->beta;
      
      ComputeSASEFELParameters(&sasefelOutput->lightWavelength[slice], 
                               &sasefelOutput->saturationLength[slice], 
                               &sasefelOutput->gainLength[slice],
                               &sasefelOutput->noisePower[slice],
                               &sasefelOutput->saturationPower[slice], 
                               &sasefelOutput->PierceParameter[slice],
                               &sasefelOutput->etaDiffraction[slice],
                               &sasefelOutput->etaEmittance[slice], 
                               &sasefelOutput->etaEnergySpread[slice],
                               sasefelOutput->charge[slice], 
                               sasefelOutput->rmsBunchLength[slice], 
                               sasefelOutput->undulatorPeriod,
                               sasefelOutput->undulatorK, 
                               sasefelOutput->betaToUse[slice], sasefelOutput->emit[slice], 
                               sasefelOutput->Sdelta[slice], sasefelOutput->pCentral[slice],
                               1);
      sasefelOutput->lightWavelength[nSlices+1] += sasefelOutput->lightWavelength[slice];
      sasefelOutput->saturationLength[nSlices+1] += sasefelOutput->saturationLength[slice];
      sasefelOutput->gainLength[nSlices+1] += sasefelOutput->gainLength[slice];
      sasefelOutput->noisePower[nSlices+1] += sasefelOutput->noisePower[slice];
      sasefelOutput->saturationPower[nSlices+1] += sasefelOutput->saturationPower[slice];
      sasefelOutput->PierceParameter[nSlices+1] += sasefelOutput->PierceParameter[slice];
      sasefelOutput->etaDiffraction[nSlices+1] += sasefelOutput->etaDiffraction[slice];
      sasefelOutput->etaEmittance[nSlices+1] += sasefelOutput->etaEmittance[slice];
      sasefelOutput->etaEnergySpread[nSlices+1] += sasefelOutput->etaEnergySpread[slice];
      sasefelOutput->betaToUse[nSlices+1] += sasefelOutput->betaToUse[slice];
      sasefelOutput->emit[nSlices+1] += sasefelOutput->emit[slice];
      sasefelOutput->Sdelta[nSlices+1] += sasefelOutput->Sdelta[slice];
      sasefelOutput->pCentral[nSlices+1] += sasefelOutput->pCentral[slice];
      sasefelOutput->charge[nSlices+1] += sasefelOutput->charge[slice];
      sasefelOutput->rmsBunchLength[nSlices+1] += sasefelOutput->rmsBunchLength[slice];
    }
    
    if (!slicesFound)
      bomb("No valid slices found for SASE FEL computation.", NULL);

    sasefelOutput->lightWavelength[nSlices+1] /= slicesFound;
    sasefelOutput->saturationLength[nSlices+1] /= slicesFound;
    sasefelOutput->saturationPower[nSlices+1] /= slicesFound;
    sasefelOutput->gainLength[nSlices+1] /= slicesFound;
    sasefelOutput->noisePower[nSlices+1] /= slicesFound;
    sasefelOutput->PierceParameter[nSlices+1] /= slicesFound;
    sasefelOutput->etaDiffraction[nSlices+1] /= slicesFound;
    sasefelOutput->etaEmittance[nSlices+1] /= slicesFound;
    sasefelOutput->etaEnergySpread[nSlices+1] /= slicesFound;
    sasefelOutput->betaToUse[nSlices+1] /= slicesFound;
    sasefelOutput->emit[nSlices+1] /= slicesFound;
    sasefelOutput->Sdelta[nSlices+1] /= slicesFound;
    sasefelOutput->pCentral[nSlices+1] /= slicesFound;
    sasefelOutput->charge[nSlices+1] /= slicesFound;
    sasefelOutput->rmsBunchLength[nSlices+1] /= slicesFound;
            
  }
  
  free(time);
}

void storeSASEFELAtEndInRPN(SASEFEL_OUTPUT *sasefelOutput)
{
  rpn_store(sasefelOutput->lightWavelength[0],
            rpn_create_mem("SASE.lightWavelength"));
  rpn_store(sasefelOutput->gainLength[0],
            rpn_create_mem("SASE.gainLength"));
  rpn_store(sasefelOutput->saturationPower[0],
            rpn_create_mem("SASE.saturationPower"));
  rpn_store(sasefelOutput->saturationLength[0],
            rpn_create_mem("SASE.saturationLength"));
}
