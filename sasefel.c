/* Copyright 1999 by Michael Borland and Argonne National Laboratory,
 * all rights reserved.
 */
/*
 * $Log: not supported by cvs2svn $
 */
#include "mdb.h"
#include "track.h"
#include "sasefel.h"

void setupSASEFELAtEnd(NAMELIST_TEXT *nltext, RUN *run, OUTPUT_FILES *output_data)
{
  SASEFEL_OUTPUT *sasefelOutput;
  SDDS_DATASET *SDDSout;
  
  /* process namelist text */
  process_namelist(&sasefel, nltext);
  print_namelist(stderr, &sasefel);

  if (beta<0)
    bomb("beta < 0", NULL);
  if (undulator_K<=0)
    bomb("undulator_K <= 0", NULL);
  if (undulator_period<=0)
    bomb("undulator_period <= 0", NULL);
  
  sasefelOutput = &(output_data->sasefel);
  if (sasefelOutput->active && sasefelOutput->filename) {
    SDDS_Terminate(&(sasefelOutput->SDDSout));
    SDDS_ClearErrors();
  }
  
  sasefelOutput->active = 1;
  sasefelOutput->beta = beta;
  sasefelOutput->undulatorK = undulator_K;
  sasefelOutput->undulatorPeriod = undulator_period;
  if (output) {
    sasefelOutput->filename = compose_filename(output, run->rootname);
    SDDSout = &(sasefelOutput->SDDSout);
    if (!SDDS_InitializeOutput(SDDSout, SDDS_BINARY, 0, NULL, NULL, sasefelOutput->filename) ||
        !SDDS_DefineSimpleParameter(SDDSout, "Step", NULL, SDDS_LONG) ||
        !SDDS_DefineSimpleParameter(SDDSout, "undulatorK", NULL, SDDS_DOUBLE) ||
        !SDDS_DefineSimpleParameter(SDDSout, "undulatorPeriod", "m", SDDS_DOUBLE) ||
        !SDDS_DefineSimpleParameter(SDDSout, "beta", NULL, SDDS_DOUBLE) ||
        !SDDS_DefineSimpleParameter(SDDSout, "charge", "C", SDDS_DOUBLE) ||
        !SDDS_DefineSimpleParameter(SDDSout, "rmsBunchLength", "s", SDDS_DOUBLE) ||
        !SDDS_DefineSimpleParameter(SDDSout, "Sdelta", NULL, SDDS_DOUBLE) ||
        !SDDS_DefineSimpleParameter(SDDSout, "emit", NULL, SDDS_DOUBLE) ||
        !SDDS_DefineSimpleParameter(SDDSout, "pCentral", "m$be$nc", SDDS_DOUBLE) ||
        !SDDS_DefineSimpleParameter(SDDSout, "lightWavelength", "m", SDDS_DOUBLE) ||
        !SDDS_DefineSimpleParameter(SDDSout, "saturationLength", "m", SDDS_DOUBLE) ||
        !SDDS_DefineSimpleParameter(SDDSout, "gainLength", "m", SDDS_DOUBLE) ||
        !SDDS_DefineSimpleParameter(SDDSout, "noisePower", "W", SDDS_DOUBLE) ||
        !SDDS_DefineSimpleParameter(SDDSout, "saturationPower", "W", SDDS_DOUBLE) ||
        !SDDS_DefineSimpleParameter(SDDSout, "PierceParameter", NULL, SDDS_DOUBLE) ||
        !SDDS_DefineSimpleParameter(SDDSout, "etaDiffraction", NULL, SDDS_DOUBLE) ||
        !SDDS_DefineSimpleParameter(SDDSout, "etaEmittance", NULL, SDDS_DOUBLE) ||
        !SDDS_DefineSimpleParameter(SDDSout, "etaEnergySpread", NULL, SDDS_DOUBLE) ||
        !SDDS_WriteLayout(SDDSout)) {
      fprintf(stderr, "Unable define SDDS parameter for file %s\n", sasefelOutput->filename);
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    }
  }
}

void doSASEFELAtEndOutput(SASEFEL_OUTPUT *sasefelOutput, long step)
{
  SDDS_DATASET *SDDSout;
  
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
                          "beta", sasefelOutput->betaToUse,
                          "charge", sasefelOutput->charge,
                          "rmsBunchLength", sasefelOutput->rmsBunchLength,
                          "Sdelta", sasefelOutput->Sdelta,
                          "emit", sasefelOutput->emit,
                          "pCentral", sasefelOutput->pCentral,
                          "lightWavelength", sasefelOutput->lightWavelength,
                          "saturationLength", sasefelOutput->saturationLength,
                          "gainLength", sasefelOutput->gainLength,
                          "noisePower", sasefelOutput->noisePower,
                          "saturationPower", sasefelOutput->saturationPower,
                          "PierceParameter", sasefelOutput->PierceParameter,
                          "etaDiffraction", sasefelOutput->etaDiffraction,
                          "etaEmittance", sasefelOutput->etaEmittance,
                          "etaEnergySpread", sasefelOutput->etaEnergySpread,
                          NULL) ||
      !SDDS_WritePage(SDDSout)) {
    fprintf(stderr, "Unable write data to file %s\n", sasefelOutput->filename);
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
  }
}


void computeSASEFELAtEnd(SASEFEL_OUTPUT *sasefelOutput, double **particle, long particles, 
                         double Po, double charge)
{
  double emitx, emity, bunchLength, Sdelta;
  double *time, *delta, tAve, deltaAve, deltaRMS, tRMS, S11, S33;
  long i;
  double xLimit[2], percentLevel[2] = {10, 90}, tSpread;
  
  if (!(time=SDDS_Malloc(sizeof(*time)*particles)))
    SDDS_Bomb("memory allocation failure (computeSASEFELAtEnd)");
  if (!(delta=SDDS_Malloc(sizeof(*delta)*particles)))
    SDDS_Bomb("memory allocation failure (computeSASEFELAtEnd)");
  if (!particles) {
    fprintf(stderr, "no particles left---can't compute FEL parameters");
    sasefelOutput->lightWavelength = sasefelOutput->saturationLength =
        sasefelOutput->gainLength = sasefelOutput->noisePower =
            sasefelOutput->saturationPower = sasefelOutput->PierceParameter = 
              sasefelOutput->etaDiffraction = sasefelOutput->etaEmittance =
                sasefelOutput->etaEnergySpread = DBL_MAX;
    return;
  }
  
  /* find center of energy distribution */
  for (i=deltaAve=0; i<particles; i++)
    deltaAve += (delta[i] = particle[i][5]);
  deltaAve /= particles;
  /* compute rms energy spread */
  for (i=deltaRMS=0; i<particles; i++)
    deltaRMS += sqr(particle[i][5]-deltaAve);
  sasefelOutput->Sdelta = deltaRMS = sqrt(deltaRMS/particles);

  /* compute rms-equivalent time value so that Q/(sqrt(2*PI)*tRMS) is
   * a good estimate of peak current.  I use the 10% and 90% points of
   * the distribution to compute peak current, then get equivalent tRMS.
   */
  computeTimeCoordinates(time, Po, particle, particles);
  compute_percentiles(xLimit, percentLevel, 2, time, particles);
  sasefelOutput->rmsBunchLength = tRMS 
    = (xLimit[1] - xLimit[0])/(0.8*sqrt(2*PI));

  emitx = rms_emittance(particle, 0, 1, particles, &S11, NULL, NULL);
  emity = rms_emittance(particle, 2, 3, particles, &S33, NULL, NULL);
  sasefelOutput->emit = sqrt(emitx*emity);
  sasefelOutput->pCentral = Po;
  sasefelOutput->charge = charge;  

  if (sasefelOutput->beta==0)
    sasefelOutput->betaToUse = sqrt(S11*S33/(emitx*emity));
  else
    sasefelOutput->betaToUse = sasefelOutput->beta;

  ComputeSASEFELParameters(&sasefelOutput->lightWavelength, &sasefelOutput->saturationLength, 
                           &sasefelOutput->gainLength,
                           &sasefelOutput->noisePower, &sasefelOutput->saturationPower, 
                           &sasefelOutput->PierceParameter,
                           &sasefelOutput->etaDiffraction, &sasefelOutput->etaEmittance, 
                           &sasefelOutput->etaEnergySpread,
                           charge, tRMS, sasefelOutput->undulatorPeriod,
                           sasefelOutput->undulatorK, 
                           sasefelOutput->betaToUse, sasefelOutput->emit, 
                           sasefelOutput->Sdelta, sasefelOutput->pCentral=Po,
                           1);
}

void storeSASEFELAtEndInRPN(SASEFEL_OUTPUT *sasefelOutput)
{
  rpn_store(sasefelOutput->lightWavelength,
            rpn_create_mem("SASE.lightWavelength"));
  rpn_store(sasefelOutput->gainLength,
            rpn_create_mem("SASE.gainLength"));
  rpn_store(sasefelOutput->saturationPower,
            rpn_create_mem("SASE.saturationPower"));
  rpn_store(sasefelOutput->saturationLength,
            rpn_create_mem("SASE.saturationLength"));
}
