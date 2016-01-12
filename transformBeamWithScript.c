/*************************************************************************\
 * Copyright (c) 2002 The University of Chicago, as Operator of Argonne
 * National Laboratory.
 * Copyright (c) 2002 The Regents of the University of California, as
 * Operator of Los Alamos National Laboratory.
 * This file is distributed subject to a Software License Agreement found
 * in the file LICENSE that is included with this distribution. 
\*************************************************************************/

#include "mdb.h"
#include "mdbsun.h"
#include "track.h"
#ifdef USE_GSL
#include "gsl/gsl_poly.h"
#endif

long transformBeamWithScript(SCRIPT *script, double pCentral, CHARGE *charge, 
                             BEAM *beam, double **part, long np, long *nLost,
                             char *mainRootname, long iPass, long driftOrder)
{
#if USE_MPI
  return transformBeamWithScript_p(script, pCentral, charge, beam, part, np,  nLost, mainRootname, iPass, driftOrder);
#else
  return transformBeamWithScript_s(script, pCentral, charge, beam, part, np,  nLost, mainRootname, iPass, driftOrder);
#endif
}

long transformBeamWithScript_s(SCRIPT *script, double pCentral, CHARGE *charge, 
                             BEAM *beam, double **part, long np, long *nLost,
                             char *mainRootname, long iPass, long driftOrder)
{
  char *rootname=NULL, *input, *output=NULL;
  char *cmdBuffer0, *cmdBuffer1=NULL;
  SDDS_DATASET SDDSout, SDDSin;
  double *data = NULL;
  char *dataname[6] = {"x","xp","y","yp","t","p"};
  long i, j, npNew, nameLength, doDrift;
  char passString[20];
  long k, lostIndex;

  doDrift = 0;
  if (script->onPass>=0) {
    if (script->onPass!=iPass)
      doDrift = 1;
  } else {
    if (script->startPass>=0 && script->startPass>iPass)
      doDrift = 1;
    if (script->endPass>=0 && iPass>script->endPass)
      doDrift = 1;
    if (script->passInterval>0) {
      if (script->startPass<0)
        bombElegant("Problem with script element: START_PASS<0 but PASS_INTERVAL>0", NULL);
      if ((iPass-script->startPass)%script->passInterval!=0)
        doDrift = 1;
    }
  }

  if (doDrift) {
    drift_beam(part, np, script->length, driftOrder);
    return np;
  }
  
  if (!script->rootname || !strlen(script->rootname)) {
    /* generate random rootname */
    if (!(rootname = tmpname(NULL)))
      bombElegant("problem generating temporary filename for script", NULL);
  } else 
    rootname = compose_filename(script->rootname, mainRootname);
  if (!rootname)
    bombElegant("problem generating temporary rootname for script", NULL);
  nameLength = (script->directory?strlen(script->directory):0) + \
    strlen(rootname) + strlen(script->inputExtension) +
    strlen(script->outputExtension) + 4;
  if (!(input = malloc(sizeof(*input)*nameLength)) ||
      !(output = malloc(sizeof(*output)*nameLength)))
    bombElegant("problem generating temporary filename for script", NULL);

  /* prepare command */
  if (script->directory && strlen(script->directory)) {
#if defined(_WIN32)
    sprintf(input, "%s\\%s.%s", script->directory, rootname, script->inputExtension);
    sprintf(output, "%s\\%s.%s", script->directory, rootname, script->outputExtension);
#else
    sprintf(input, "%s/%s.%s", script->directory, rootname, script->inputExtension);
    sprintf(output, "%s/%s.%s", script->directory, rootname, script->outputExtension);
#endif
  } else {
    sprintf(input, "%s.%s", rootname, script->inputExtension);
    sprintf(output, "%s.%s", rootname, script->outputExtension);
  }
  if (rootname!=script->rootname)
    free(rootname);

 sprintf(passString, "%ld", iPass);

  if (!(cmdBuffer0=malloc(sizeof(char)*(strlen(script->command)+10*strlen(input)+10*strlen(output)+strlen(passString)))) ||
      !(cmdBuffer1=malloc(sizeof(char)*(strlen(script->command)+10*strlen(input)+10*strlen(output)+strlen(passString)))))
    bombElegant("memory allocation failure making command buffer for script", NULL);
  replaceString(cmdBuffer0, script->command, "%i", input, 9, 0);
  replaceString(cmdBuffer1, cmdBuffer0, "%o", output, 9, 0);
 
  replaceString(cmdBuffer0, cmdBuffer1, "%p", passString, 9, 0);
  strcpy_ss(cmdBuffer1, cmdBuffer0);
  
  /* substitute numerical parameters */
  for (i=0; i<10; i++) {
    long count = 0;
    char tag[10], value[25], *ptr;
    sprintf(tag, "%%np%ld", i);
    ptr = cmdBuffer1;
    while ((ptr=strstr(ptr, tag))) {
      count ++;
      ptr += 3;
    }
    if (!count) continue;
    sprintf(value, "%21.15e", script->NP[i]);
    if (!(cmdBuffer0 = SDDS_Realloc(cmdBuffer0, sizeof(*cmdBuffer1)*(strlen(cmdBuffer1)+count*25+1))) ||
        !(cmdBuffer1 = SDDS_Realloc(cmdBuffer1, sizeof(*cmdBuffer1)*(strlen(cmdBuffer1)+count*25+1))))
      SDDS_Bomb("memory allocation failure");
    replaceString(cmdBuffer0, cmdBuffer1, tag, value, count, 0);
    strcpy_ss(cmdBuffer1, cmdBuffer0);
  }
  /* substitute string parameters */
  for (i=0; i<10; i++) {
    long count = 0;
    char tag[10], *ptr;
    if (!script->SP[i] || strlen(script->SP[i])==0)
      continue;
    sprintf(tag, "%%sp%ld", i);
    ptr = cmdBuffer1;
    while ((ptr=strstr(ptr, tag))) {
      count ++;
      ptr += 3;
    }
    if (!count) continue;
    if (!(cmdBuffer0 = 
          SDDS_Realloc(cmdBuffer0, sizeof(*cmdBuffer1)*(strlen(cmdBuffer1)+count*strlen(script->SP[i])+1))) ||
        !(cmdBuffer1 = 
          SDDS_Realloc(cmdBuffer1, sizeof(*cmdBuffer1)*(strlen(cmdBuffer1)+count*strlen(script->SP[i])+1))))
      SDDS_Bomb("memory allocation failure");
    replaceString(cmdBuffer0, cmdBuffer1, tag, script->SP[i], count, 0);
    strcpy_ss(cmdBuffer1, cmdBuffer0);
  }
  interpret_escaped_quotes(cmdBuffer1);
  
  if (script->verbosity>0) {
    fprintf(stdout, "%s\n", cmdBuffer1);
    fflush(stdout);
  }
 
  /* dump the data to script input file */
  SDDS_ForceInactive(&SDDSout);
  SDDS_PhaseSpaceSetup(&SDDSout, input, SDDS_BINARY, 1, "script input", 
                       "unknown", "unknown",
                       "transformBeamWithScript");
  dump_phase_space(&SDDSout, part, np, 0, pCentral, charge?charge->macroParticleCharge*np                    :0.0, beam?beam->id_slots_per_bunch:0);

  if (!SDDS_Terminate(&SDDSout))
    SDDS_Bomb("problem terminating script input file");

#if defined(CONDOR_COMPILE)
  _condor_ckpt_disable();
#endif

  /* run the script */
  if (script->useCsh)
    executeCshCommand(cmdBuffer1);
  else 
    system(cmdBuffer1);

#if defined(CONDOR_COMPILE)
  _condor_ckpt_enable();
#endif

  if (script->verbosity>0) {
    fprintf(stdout, "Command completed\n");
    fflush(stdout);
  }

  /* read the data from script output file */
  if (!fexists(output)) 
    SDDS_Bomb("unable to find script output file");
  if (!SDDS_InitializeInput(&SDDSin, output)) {
    SDDS_SetError("Unable to read script output file");
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
  }

  if (!check_sdds_column(&SDDSin, "x", "m") ||
      !check_sdds_column(&SDDSin, "y", "m") ||
      !check_sdds_column(&SDDSin, "xp", NULL) ||
      !check_sdds_column(&SDDSin, "yp", NULL) ||
      !check_sdds_column(&SDDSin, "p", "m$be$nc") ||
      !check_sdds_column(&SDDSin, "t", "s")) {
    if (!check_sdds_column(&SDDSin, "p", "m$be$nc") &&
        check_sdds_column(&SDDSin, "p", NULL)) {
      fprintf(stdout, "Warning: p has no units in script output file.  Expected m$be$nc\n");
      fflush(stdout);
    } else {
      fprintf(stdout, 
              "necessary data quantities (x, x', y, y', t, p) have the wrong units or are not present in script output");
      fflush(stdout);
      exitElegant(1);
    }
  }
 
  if (SDDS_ReadPage(&SDDSin)!=1) {
    SDDS_SetError("Unable to read script output file");
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
  }
  npNew = SDDS_RowCount(&SDDSin);

  if (script->verbosity>0) {
    fprintf(stdout, "%ld particles in script output file (was %ld)\n", npNew, np);
    fflush(stdout);
  }

  if (!npNew) {
    return 0;
  }
  if (npNew>np)
    if (script->noNewParticles)
      bombElegant("The number of particles increased after the SCRIPT element without seting the correct flag.\n Please set the correct NO_NEW_PARTICLES flag for the SCRIPT element!", NULL);

  if (!script->useParticleID) {
    if ((!script->noNewParticles || (npNew!=np)) && (iPass==0)){
      fprintf (stdout, "Warning: There is no particle ID available to find which particles are lost. The particles lost in the SCRIPT element will not be recorded!\n");
    }
    if ((!script->noNewParticles) && (npNew < np)) {
      /* Do particle accounting only */
      beam->n_to_track = npNew+*nLost;
      /* move lost particles into the upper part of the arrays */
      if (beam->lost)
	for (i=0; i<np-npNew; i++) {
	  swapParticles(beam->lost[npNew+i], beam->lost[beam->n_to_track+i]);
	}
    }
  }    
      
  if (npNew>np) {
    /* We may have to resize the arrays in the BEAM structure */
    
    fprintf(stdout, "Increasing number of particles from %ld (%ld active) to %ld (%ld active)\n",
            np+(nLost?*nLost:0), np, npNew+(nLost?*nLost:0), npNew);
    fflush(stdout);
    
    if (!beam) {
      fprintf(stderr, "Error: script element increased the number of particles from %ld to %ld\n.",
              np, npNew);
      fprintf(stderr, "This happened (apparently) during a pre-tracking stage, which isn't allowed\n");
      exitElegant(1);
    }
    /* Check that previous particle counts are correct */
    if ((np+(nLost?*nLost:0))!=beam->n_to_track) {
      fprintf(stderr, "Particle accounting problem in SCRIPT element:\n");
      fprintf(stderr, "np = %ld, *nLost = %ld, beam->n_to_track = %ld, beam->n_particle=%ld\n",
              np, (nLost?*nLost:0), beam->n_to_track, beam->n_particle);
      fprintf(stderr, "This could happen if the particleID is not unique.\n");
      exitElegant(1);
    }

    if ((npNew+(nLost?*nLost:0)) > beam->n_particle) {
      if (beam->original==beam->particle) {
        /* This means, oddly enough, that the particle array and original array are the same because the
         * separate original array wasn't needed.  n_original gives the size of both arrays (including
         * live and lost particles).  To avoid confusion, we'll copy the data to a new array before
         * doing anything else, even though it means the original array is not used for anything and
         * contains a useless frozen copy of the present beam.
         * Use n_original since that's the size of the array, including lost particles. 
         */
        beam->particle = (double**)czarray_2d(sizeof(double), beam->n_original, 7);
        copy_particles(beam->particle, beam->original, beam->n_original);
      }
      /* resize the particle array, leaving space for the lost particle data at the top */
      if (!(beam->particle = (double**)resize_czarray_2d((void**)beam->particle,sizeof(double), npNew+(nLost?*nLost:0), 7)) ||
          !(beam->lost = (double**)resize_czarray_2d((void**)beam->lost,sizeof(double), npNew+(nLost?*nLost:0), 8))) {
        fprintf(stderr, "Memory allocation failure increasing particle array size to %ld\n",
                npNew+(nLost?*nLost:0));
        exitElegant(1);
      }
      beam->n_particle = npNew+(nLost?*nLost:0);
      /* move lost particles into the upper part of the arrays */
      if (beam->lost)
	for (i=(nLost?*nLost:0)-1; i>=0; i--) {
	  swapParticles(beam->lost[np+i], beam->lost[npNew+i]);
      }
    }
    
    if (beam->accepted)  {
      /* this data is invalid when particles are added */
      free_czarray_2d((void**)beam->accepted, np+(nLost?*nLost:0), 7);
      beam->accepted = NULL;
    }
    beam->n_to_track = npNew+(nLost?*nLost:0);
    fprintf(stdout, "beam->n_particle = %ld, beam->n_to_track = %ld\n",
	    beam->n_particle, beam->n_to_track);
   
    part = beam->particle;
  }

  for (i=0; i<6; i++) {
    if (!(data = SDDS_GetColumnInDoubles(&SDDSin, dataname[i]))) {
      SDDS_SetError("Unable to read script output file");
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
    }
    for (j=0; j<npNew; j++) {
      part[j][i] = data[j];
    }
    free(data);
    data = NULL;
  }

  if (script->useParticleID ) {
    if (!(data = SDDS_GetColumnInDoubles(&SDDSin, "particleID"))) {
      SDDS_SetError("Unable to read particleID from script output file. Please set USE_PARTICLE_ID=0 if this is desired.\n");
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
    }
    if (!script->noNewParticles) {
      if (iPass==0)
        fprintf (stdout, "Warning: New particles are added in the SCRIPT element. The particles lost in the SCRIPT element will not be recorded!\n");
      for (j=0; j<npNew; j++)
        part[j][6] = data[j];
    } else {
      if (npNew<np) {
	/* Find out which particles are lost */
	lostIndex = npNew;
	for (i=0; i<np; i++) {
	  for (j=0; j<npNew; j++) {
	    if (part[i][6]==data[j])
	      break; /* this particle survived */
	  }
	  if (j==npNew) { /* record the lost particle */
	    /* Copy the lost particles in the SCRIPT element into the upper part of the arrays. */
	    if (nLost)
              (*nLost)++;
	    if (beam->lost) {
	      for (k=0; k<7; k++) {
		beam->lost[lostIndex][k] = part[i][k];
	      }
	      beam->lost[lostIndex][7] = (double) iPass; 
	      lostIndex++;
	    }
	  }
	}
      }
    }
    if (data)
      free(data);
  }

  /* assign new particle IDs if there are new particles */
  if (npNew>np && !script->useParticleID) {
    fprintf(stdout, "Changing particle ID for new particles to make them sequential\n");
    for (j=0; j<npNew; j++)
      part[j][6] = j+1;
  }

  if (charge) {
    double totalCharge, oldMacroParticleCharge;
    oldMacroParticleCharge = charge->macroParticleCharge;
    if (!SDDS_GetParameterAsDouble(&SDDSin, "Charge", &totalCharge)) {
      SDDS_SetError("Unable to read Charge parameter from script output file");
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
    }
    charge->charge = totalCharge;
    charge->macroParticleCharge = 0;
    if (npNew)
      charge->macroParticleCharge = totalCharge/npNew;
    if (oldMacroParticleCharge!=0 && fabs((charge->macroParticleCharge-oldMacroParticleCharge)/oldMacroParticleCharge)>1e-8) {
      printf("*** Warning: macro-particle charge changed after SCRIPT element, from %le to %le. This may indicate a problem.\n",
             oldMacroParticleCharge, charge->macroParticleCharge);
      printf("             Please ensure that the Charge parameter is set correctly in the output file from your script.\n");
    }
 }

  if (SDDS_ReadPage(&SDDSin)!=-1)
    SDDS_Bomb("Script output file has multiple pages");
  if (!SDDS_Terminate(&SDDSin))
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);

  if (script->verbosity) {
    fprintf(stdout, "Done processing particle input file from script\n");
    fflush(stdout);
  }

  /* convert (t, p) data to (s, delta) */
  for (j=0; j<npNew; j++) {
    double p, beta;
    p = part[j][5];
    part[j][5] = (p-pCentral)/pCentral;
    beta = p/sqrt(sqr(p)+1);
    part[j][4] *= beta*c_mks;
  }
  if (!script->useParticleID) {
    fprintf(stdout, "Assigning fresh IDs to all particles.\n");
    for (j=0; j<npNew; j++) {
      part[j][6] = j+1;
    }
  }
  
  if (!script->keepFiles) {
    /* delete the input and output files */
    remove(input);
    remove(output);
  }

  /* clean up */
  free(cmdBuffer0);
  free(cmdBuffer1);

  return npNew;
}

long transformBeamWithScript_p(SCRIPT *script, double pCentral, CHARGE *charge, 
                             BEAM *beam, double **part, long np, long *nLost,
                             char *mainRootname, long iPass, long driftOrder)
{
  char *rootname=NULL, *input, *output=NULL;
  char *cmdBuffer0, *cmdBuffer1=NULL;
  SDDS_DATASET SDDSout, SDDSin;
  double *data = NULL;
  char *dataname[6] = {"x","xp","y","yp","t","p"};
  long i, j, npNew, nameLength, doDrift;
  char passString[20];
#if !USE_MPI
  long k, lostIndex;
#else
  long npTotal=0, rootnameLength;
#endif

  doDrift = 0;
  if (script->onPass>=0) {
    if (script->onPass!=iPass)
      doDrift = 1;
  } else {
    if (script->startPass>=0 && script->startPass>iPass)
      doDrift = 1;
    if (script->endPass>=0 && iPass>script->endPass)
      doDrift = 1;
    if (script->passInterval>0) {
      if (script->startPass<0)
        bombElegant("Problem with script element: START_PASS<0 but PASS_INTERVAL>0", NULL);
      if ((iPass-script->startPass)%script->passInterval!=0)
        doDrift = 1;
    }
  }

  if (doDrift) {
    drift_beam(part, np, script->length, driftOrder);
    return np;
  }
  
#if USE_MPI
  if (notSinglePart) { 
    MPI_Allreduce (&np, &(beam->n_to_track_total), 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
#if MPI_DEBUG
    printf("n_to_track_total = %ld\n", beam->n_to_track_total);
#endif
    if (!beam->n_to_track_total)
       return 0;
  }	
#endif

  if (!script->rootname || !strlen(script->rootname)) {
    /* generate random rootname */
    if (isMaster && !(rootname = tmpname(NULL)))
      bombElegant("problem generating temporary filename for script", NULL);
#if SDDS_MPI_IO
    if (isMaster)
      rootnameLength = strlen(rootname)+1;
    MPI_Bcast(&rootnameLength, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    if (isSlave) /* Master and slave could have different rootname length if use C-shell, which calls tmpname one more time on master for execution.*/
      rootname = malloc(rootnameLength); 
    /* As different processors will have different process names, we need
     make sure they have the same file name for parallel I/O */
    MPI_Bcast(rootname, rootnameLength, MPI_CHAR, 0, MPI_COMM_WORLD);
#if MPI_DEBUG
    printf("rootname = %s\n", rootname);
#endif
#endif
  } else 
    rootname = compose_filename(script->rootname, mainRootname);
  if (!rootname)
    bombElegant("problem generating temporary rootname for script", NULL);
  nameLength = (script->directory?strlen(script->directory):0) + \
    strlen(rootname) + strlen(script->inputExtension) +
    strlen(script->outputExtension) + 4;
  if (!(input = malloc(sizeof(*input)*nameLength)) ||
      !(output = malloc(sizeof(*output)*nameLength)))
    bombElegant("problem generating temporary filename for script", NULL);

  if (0) {
    long pidMin, pidMax;
    pidMin = LONG_MAX;
    pidMax = LONG_MIN;
    for (j=0; j<np; j++) {
      if (pidMin>part[j][6])
        pidMin = part[j][6];
      if (pidMax<part[j][6])
        pidMax = part[j][6];
    }
    fprintf(stdout, "before, pid limits: %ld, %ld\n", pidMin, pidMax);
  }
  
  /* prepare command */
  if (script->directory && strlen(script->directory)) {
#if defined(_WIN32)
    sprintf(input, "%s\\%s.%s", script->directory, rootname, script->inputExtension);
    sprintf(output, "%s\\%s.%s", script->directory, rootname, script->outputExtension);
#else
    sprintf(input, "%s/%s.%s", script->directory, rootname, script->inputExtension);
    sprintf(output, "%s/%s.%s", script->directory, rootname, script->outputExtension);
#endif
  } else {
    sprintf(input, "%s.%s", rootname, script->inputExtension);
    sprintf(output, "%s.%s", rootname, script->outputExtension);
  }
  if (rootname!=script->rootname)
    free(rootname);

 sprintf(passString, "%ld", iPass);

  if (!(cmdBuffer0=malloc(sizeof(char)*(strlen(script->command)+10*strlen(input)+10*strlen(output)+strlen(passString)))) ||
      !(cmdBuffer1=malloc(sizeof(char)*(strlen(script->command)+10*strlen(input)+10*strlen(output)+strlen(passString)))))
    bombElegant("memory allocation failure making command buffer for script", NULL);
  replaceString(cmdBuffer0, script->command, "%i", input, 9, 0);
  replaceString(cmdBuffer1, cmdBuffer0, "%o", output, 9, 0);
 
  replaceString(cmdBuffer0, cmdBuffer1, "%p", passString, 9, 0);
  strcpy_ss(cmdBuffer1, cmdBuffer0);
  
  /* substitute numerical parameters */
  for (i=0; i<10; i++) {
    long count = 0;
    char tag[10], value[25], *ptr;
    sprintf(tag, "%%np%ld", i);
    ptr = cmdBuffer1;
    while ((ptr=strstr(ptr, tag))) {
      count ++;
      ptr += 3;
    }
    if (!count) continue;
    sprintf(value, "%21.15e", script->NP[i]);
    if (!(cmdBuffer0 = SDDS_Realloc(cmdBuffer0, sizeof(*cmdBuffer1)*(strlen(cmdBuffer1)+count*25+1))) ||
        !(cmdBuffer1 = SDDS_Realloc(cmdBuffer1, sizeof(*cmdBuffer1)*(strlen(cmdBuffer1)+count*25+1))))
      SDDS_Bomb("memory allocation failure");
    replaceString(cmdBuffer0, cmdBuffer1, tag, value, count, 0);
    strcpy_ss(cmdBuffer1, cmdBuffer0);
  }
  /* substitute string parameters */
  for (i=0; i<10; i++) {
    long count = 0;
    char tag[10], *ptr;
    if (!script->SP[i] || strlen(script->SP[i])==0)
      continue;
    sprintf(tag, "%%sp%ld", i);
    ptr = cmdBuffer1;
    while ((ptr=strstr(ptr, tag))) {
      count ++;
      ptr += 3;
    }
    if (!count) continue;
    if (!(cmdBuffer0 = 
          SDDS_Realloc(cmdBuffer0, sizeof(*cmdBuffer1)*(strlen(cmdBuffer1)+count*strlen(script->SP[i])+1))) ||
        !(cmdBuffer1 = 
          SDDS_Realloc(cmdBuffer1, sizeof(*cmdBuffer1)*(strlen(cmdBuffer1)+count*strlen(script->SP[i])+1))))
      SDDS_Bomb("memory allocation failure");
    replaceString(cmdBuffer0, cmdBuffer1, tag, script->SP[i], count, 0);
    strcpy_ss(cmdBuffer1, cmdBuffer0);
  }
  interpret_escaped_quotes(cmdBuffer1);
  
  if (script->verbosity>0) {
    fprintf(stdout, "%s\n", cmdBuffer1);
    fflush(stdout);
  }
 
  /* dump the data to script input file */
#if USE_MPI
  if (notSinglePart || (!notSinglePart&&isMaster))
#endif
  {
    SDDS_ForceInactive(&SDDSout);
    SDDS_PhaseSpaceSetup(&SDDSout, input, SDDS_BINARY, 1, "script input", 
			 "unknown", "unknown",
			 "transformBeamWithScript");
#if USE_MPI
    dump_phase_space(&SDDSout, part, np, 0, pCentral, charge?charge->macroParticleCharge*beam->n_to_track_total:0.0, beam?beam->id_slots_per_bunch:0);
#else
    dump_phase_space(&SDDSout, part, np, 0, pCentral, charge?charge->macroParticleCharge*np                    :0.0, beam?beam->id_slots_per_bunch:0);
#endif

    if (!SDDS_Terminate(&SDDSout))
      SDDS_Bomb("problem terminating script input file");
  }
#if defined(CONDOR_COMPILE)
  _condor_ckpt_disable();
#endif

  /* run the script */
  if (isMaster) /* This will be done on the master */
  {  
    if (script->useCsh)
      executeCshCommand(cmdBuffer1);
    else 
      system(cmdBuffer1);
  }
#if defined(CONDOR_COMPILE)
  _condor_ckpt_enable();
#endif

  if (script->verbosity>0) {
    fprintf(stdout, "Command completed\n");
    fflush(stdout);
  }

  /* read the data from script output file */
#if SDDS_MPI_IO
#ifdef MPI_DEBUG
  printf("Waiting on MPI barrier, notSinglePart = %ld\n", notSinglePart);
  fflush(stdout);
#endif

  MPI_Barrier(MPI_COMM_WORLD);

  if (notSinglePart) {
    if (!fexists(output)) 
      SDDS_Bomb("unable to find script output file");
    SDDSin.parallel_io = 1;
    /* set up parallel IO information */      
    SDDS_MPI_Setup(&SDDSin, 1, n_processors, myid, MPI_COMM_WORLD, 0);
    if (!SDDS_MPI_InitializeInput(&SDDSin, output)) {
      SDDS_SetError("Unable to read script output file");
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
    }
  } else if (isMaster){
    if (!fexists(output)) 
      SDDS_Bomb("unable to find script output file");
    if (!SDDS_InitializeInput(&SDDSin, output)) {
      SDDS_SetError("Unable to read script output file");
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
    }
  } else 
    bombElegant("Invalid branch for SCRIPT element. Seek expert help.", NULL);
#else  
  if (!fexists(output)) 
    SDDS_Bomb("unable to find script output file");
  if (!SDDS_InitializeInput(&SDDSin, output)) {
    SDDS_SetError("Unable to read script output file");
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
  }
#endif

#if USE_MPI
  if (notSinglePart || (!notSinglePart&&isMaster))
#endif
  if (!check_sdds_column(&SDDSin, "x", "m") ||
      !check_sdds_column(&SDDSin, "y", "m") ||
      !check_sdds_column(&SDDSin, "xp", NULL) ||
      !check_sdds_column(&SDDSin, "yp", NULL) ||
      !check_sdds_column(&SDDSin, "p", "m$be$nc") ||
      !check_sdds_column(&SDDSin, "t", "s")) {
    if (!check_sdds_column(&SDDSin, "p", "m$be$nc") &&
        check_sdds_column(&SDDSin, "p", NULL)) {
      fprintf(stdout, "Warning: p has no units in script output file.  Expected m$be$nc\n");
      fflush(stdout);
    } else {
      fprintf(stdout, 
              "necessary data quantities (x, x', y, y', t, p) have the wrong units or are not present in script output");
      fflush(stdout);
      exitElegant(1);
    }
  }
 
#if !SDDS_MPI_IO  
  if (SDDS_ReadPage(&SDDSin)!=1) {
    SDDS_SetError("Unable to read script output file");
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
  }
#else
  if (notSinglePart)
    SDDS_MPI_ReadPage(&SDDSin);
  else if (isMaster) 
    SDDS_ReadPage(&SDDSin);
#ifdef MPI_DEBUG
  printf("Finished parallel read of SCRIPT output file\n");
  fflush(stdout);
#endif
  if (notSinglePart || (!notSinglePart&&isMaster))
#endif
  npNew = SDDS_RowCount(&SDDSin);

#if !USE_MPI
  if (script->verbosity>0) {
    fprintf(stdout, "%ld particles in script output file (was %ld)\n", npNew, np);
    fflush(stdout);
  }

  if (!npNew) {
    return 0;
  }
  if (npNew>np)
#else
  if (!notSinglePart) 
    MPI_Bcast(&npNew, 1, MPI_LONG, 0, MPI_COMM_WORLD);
  else {
    if (beam) {
      npTotal = beam->n_to_track_total;
      beam->n_to_track_total = SDDS_MPI_TotalRowCount(&SDDSin);
      if (script->verbosity>0) {
        fprintf(stdout, "%ld particles in script output file (was %ld)\n", beam->n_to_track_total, npTotal);
        fflush(stdout);
      }	
    }
  }
  if ((!notSinglePart && (npNew>np))||(notSinglePart && (beam->n_to_track_total>npTotal)))
#endif
  if (script->noNewParticles)
    bombElegant("The number of particles increased after the SCRIPT element without seting the correct flag.\n Please set the correct NO_NEW_PARTICLES flag for the SCRIPT element!", NULL);

  if (!script->useParticleID) {
    if ((!script->noNewParticles || (npNew!=np)) && (iPass==0)){
      fprintf (stdout, "Warning: There is no particle ID available to find which particles are lost. The particles lost in the SCRIPT element will not be recorded!\n");
    }
    if ((!script->noNewParticles) && (npNew < np)) {
      /* Do particle accounting only */
      beam->n_to_track = npNew+*nLost;
      /* move lost particles into the upper part of the arrays */
      if (beam->lost)
	for (i=0; i<np-npNew; i++) {
	  swapParticles(beam->lost[npNew+i], beam->lost[beam->n_to_track+i]);
	}
    }
  }    
      
  if (npNew>np) {
    /* We may have to resize the arrays in the BEAM structure */
    
    fprintf(stdout, "Increasing number of particles from %ld (%ld active) to %ld (%ld active)\n",
            np+(nLost?*nLost:0), np, npNew+(nLost?*nLost:0), npNew);
    fflush(stdout);
    
    if (!beam) {
      fprintf(stderr, "Error: script element increased the number of particles from %ld to %ld\n.",
              np, npNew);
      fprintf(stderr, "This happened (apparently) during a pre-tracking stage, which isn't allowed\n");
      exitElegant(1);
    }
    /* Check that previous particle counts are correct */
    if ((np+(nLost?*nLost:0))!=beam->n_to_track) {
      fprintf(stderr, "Particle accounting problem in SCRIPT element:\n");
      fprintf(stderr, "np = %ld, *nLost = %ld, beam->n_to_track = %ld, beam->n_particle=%ld\n",
              np, (nLost?*nLost:0), beam->n_to_track, beam->n_particle);
      fprintf(stderr, "This could happen if the particleID is not unique.\n");
      exitElegant(1);
    }
#ifdef MPI_DEBUG
    printf("Script check 1\n");
    fflush(stdout);
#endif

    if ((npNew+(nLost?*nLost:0)) > beam->n_particle) {
      if (beam->original==beam->particle) {
        /* This means, oddly enough, that the particle array and original array are the same because the
         * separate original array wasn't needed.  n_original gives the size of both arrays (including
         * live and lost particles).  To avoid confusion, we'll copy the data to a new array before
         * doing anything else, even though it means the original array is not used for anything and
         * contains a useless frozen copy of the present beam.
         * Use n_original since that's the size of the array, including lost particles. 
         */
        beam->particle = (double**)czarray_2d(sizeof(double), beam->n_original, 7);
        copy_particles(beam->particle, beam->original, beam->n_original);
      }
      /* resize the particle array, leaving space for the lost particle data at the top */
      if (!(beam->particle = (double**)resize_czarray_2d((void**)beam->particle,sizeof(double), npNew+(nLost?*nLost:0), 7)) ||
          !(beam->lost = (double**)resize_czarray_2d((void**)beam->lost,sizeof(double), npNew+(nLost?*nLost:0), 8))) {
        fprintf(stderr, "Memory allocation failure increasing particle array size to %ld\n",
                npNew+(nLost?*nLost:0));
        exitElegant(1);
      }
      beam->n_particle = npNew+(nLost?*nLost:0);
#if !USE_MPI
      /* move lost particles into the upper part of the arrays */
      if (beam->lost)
	for (i=(nLost?*nLost:0)-1; i>=0; i--) {
	  swapParticles(beam->lost[np+i], beam->lost[npNew+i]);
      }
#endif
    }
    
    if (beam->accepted)  {
      /* this data is invalid when particles are added */
      free_czarray_2d((void**)beam->accepted, np+(nLost?*nLost:0), 7);
      beam->accepted = NULL;
    }
    beam->n_to_track = npNew+(nLost?*nLost:0);
    fprintf(stdout, "beam->n_particle = %ld, beam->n_to_track = %ld\n",
	    beam->n_particle, beam->n_to_track);
   
    part = beam->particle;
  }

#if USE_MPI
#ifdef MPI_DEBUG
    printf("Script check 2\n");
    fflush(stdout);
#endif
  /* Particles could be redistributed, move lost particles into the upper part of the arrays */
  if ((np != npNew) && beam && nLost && beam->lost)
    for (i=0; i<=(nLost?*nLost:0)-1; i++) {
      swapParticles(beam->lost[np+i], beam->lost[npNew+i]);
    }
  if ((isSlave && notSinglePart) || (isMaster && !notSinglePart)) 
#endif
  for (i=0; i<6; i++) {
    if (!(data = SDDS_GetColumnInDoubles(&SDDSin, dataname[i]))) {
      SDDS_SetError("Unable to read script output file");
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
    }
    for (j=0; j<npNew; j++) {
      part[j][i] = data[j];
    }
    free(data);
    data = NULL;
  }

#if USE_MPI
#ifdef MPI_DEBUG
    printf("Script check 3\n");
    fflush(stdout);
#endif
  if (!notSinglePart)
    MPI_Bcast(part[0], 7*npNew, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
  if (script->useParticleID ) {
#if USE_MPI
#ifdef MPI_DEBUG
    printf("Script check 4\n");
    fflush(stdout);
#endif
    if ((isSlave && notSinglePart) || (isMaster && !notSinglePart))
#endif
    if (!(data = SDDS_GetColumnInDoubles(&SDDSin, "particleID"))) {
      SDDS_SetError("Unable to read particleID from script output file. Please set USE_PARTICLE_ID=0 if this is desired.\n");
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
    }
    if (!script->noNewParticles) {
      if (iPass==0)
        fprintf (stdout, "Warning: New particles are added in the SCRIPT element. The particles lost in the SCRIPT element will not be recorded!\n");
      for (j=0; j<npNew; j++)
        part[j][6] = data[j];
    } else {
#if !USE_MPI
      if (npNew<np) {
	/* Find out which particles are lost */
	lostIndex = npNew;
	for (i=0; i<np; i++) {
	  for (j=0; j<npNew; j++) {
	    if (part[i][6]==data[j])
	      break; /* this particle survived */
	  }
	  if (j==npNew) { /* record the lost particle */
	    /* Copy the lost particles in the SCRIPT element into the upper part of the arrays. */
	    if (nLost)
              (*nLost)++;
	    if (beam->lost) {
	      for (k=0; k<7; k++) {
		beam->lost[lostIndex][k] = part[i][k];
	      }
	      beam->lost[lostIndex][7] = (double) iPass; 
	      lostIndex++;
	    }
	  }
	}
      }
#else
#ifdef MPI_DEBUG
      printf("Script check 5\n");
      fflush(stdout);
#endif
      /* Even though particle ID is available, we will not record lost particle coordinates due to particle
	 redistribution in Pelegant. The particles move from one processor to another make it very complicated 
	 to find out which particles are lost. While it should be not hard to find the lost particles within 
	 the script if there is such an interest. */
      if (beam && (beam->n_to_track_total<npTotal) && (iPass==0))
	fprintf (stdout, "Warning: Lost particle coordinates in the SCRIPT element will not be recorded.\n");
#endif
    }
    if (data)
      free(data);
#ifdef MPI_DEBUG
      printf("Script check 5.1\n");
      fflush(stdout);
#endif
  }

  if (0) {
    long pidMin, pidMax;
    pidMin = LONG_MAX;
    pidMax = LONG_MIN;
    for (j=0; j<npNew; j++) {
      if (pidMin>part[j][6])
        pidMin = part[j][6];
      if (pidMax<part[j][6])
        pidMax = part[j][6];
    }
    fprintf(stdout, "after, np=%ld, pid limits: %ld, %ld\n", npNew, pidMin, pidMax);
  }
  
  /* assign new particle IDs if there are new particles */
#if !USE_MPI
  if (npNew>np && !script->useParticleID) {
    fprintf(stdout, "Changing particle ID for new particles to make them sequential\n");
    for (j=0; j<npNew; j++)
      part[j][6] = j+1;
  }
#else
#ifdef MPI_DEBUG
    printf("Script check 6\n");
    fflush(stdout);
#endif
  if (isSlave && notSinglePart) {
    if (beam && (beam->n_to_track_total>npTotal) && !script->useParticleID) {
      long sum=0, tmp, my_offset=0, *offset = tmalloc(n_processors*sizeof(*offset));
      MPI_Allgather (&npNew, 1, MPI_LONG, offset, 1, MPI_LONG, workers);
      tmp = offset[0];
      for (i=1; i<n_processors; i++) {
        sum += tmp;
        tmp = offset[i];
        offset[i] = sum;
      }
      offset[0] = 0;
      my_offset = offset[myid-1];
      tfree(offset);
      for (j=0; j<npNew; j++){
	part[j][6] = j+1+my_offset;
      }
    }
  }
#endif

#if USE_MPI
#ifdef MPI_DEBUG
    printf("Script check 7\n");
    fflush(stdout);
#endif
  if ((!notSinglePart&&isMaster) || notSinglePart)
#endif
  if (charge) {
    double totalCharge, oldMacroParticleCharge;
    oldMacroParticleCharge = charge->macroParticleCharge;
    if (!SDDS_GetParameterAsDouble(&SDDSin, "Charge", &totalCharge)) {
      SDDS_SetError("Unable to read Charge parameter from script output file");
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
    }
    charge->charge = totalCharge;
    charge->macroParticleCharge = 0;
#if USE_MPI
#ifdef MPI_DEBUG
    printf("Script check 8\n");
    fflush(stdout);
#endif
    if (!notSinglePart) {
      if (npNew)
	charge->macroParticleCharge = totalCharge/npNew;
    } else
      if (beam->n_to_track_total)
	charge->macroParticleCharge = totalCharge/beam->n_to_track_total;
#else
    if (npNew)
      charge->macroParticleCharge = totalCharge/npNew;
#endif
    if (oldMacroParticleCharge!=0 && fabs((charge->macroParticleCharge-oldMacroParticleCharge)/oldMacroParticleCharge)>1e-8) {
      printf("*** Warning: macro-particle charge changed after SCRIPT element, from %le to %le. This may indicate a problem.\n",
             oldMacroParticleCharge, charge->macroParticleCharge);
      printf("             Please ensure that the Charge parameter is set correctly in the output file from your script.\n");
    }
 }

#if USE_MPI
#ifdef MPI_DEBUG
    printf("Script check 9\n");
    fflush(stdout);
#endif
  if (charge && notSinglePart) {
    MPI_Bcast(&(charge->charge), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&charge->macroParticleCharge, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }
  if (notSinglePart || (isMaster && !notSinglePart)) 
#endif
  {
  if (SDDS_ReadPage(&SDDSin)!=-1)
    SDDS_Bomb("Script output file has multiple pages");
  if (!SDDS_Terminate(&SDDSin))
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
  }
  if (script->verbosity) {
    fprintf(stdout, "Done processing particle input file from script\n");
    fflush(stdout);
  }

  /* convert (t, p) data to (s, delta) */
  for (j=0; j<npNew; j++) {
    double p, beta;
    p = part[j][5];
    part[j][5] = (p-pCentral)/pCentral;
    beta = p/sqrt(sqr(p)+1);
    part[j][4] *= beta*c_mks;
  }
  if (!script->useParticleID) {
    fprintf(stdout, "Assigning fresh IDs to all particles.\n");
    for (j=0; j<npNew; j++) {
      part[j][6] = j+1;
    }
  }
  
#if USE_MPI
#ifdef MPI_DEBUG
    printf("Script check 10\n");
    fflush(stdout);
#endif
  if (isMaster)
#endif
  if (!script->keepFiles) {
    /* delete the input and output files */
    remove(input);
    remove(output);
  }

  /* clean up */
  free(cmdBuffer0);
  free(cmdBuffer1);

#if USE_MPI  
#ifdef MPI_DEBUG
    printf("Script check 11, npNew = %ld\n", npNew);
    fflush(stdout);
#endif
#endif

  return npNew;
}

