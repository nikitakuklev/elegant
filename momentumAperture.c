/*************************************************************************\
* Copyright (c) 2006 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2006 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: momentumAperture.c
 * purpose: Do tracking to find momentum aperture starting from the end of each element.
 * Ref: M. Belgrounne et al. PAC203, 896.
 *
 * Michael Borland, 2006
 */

#define DEBUG 1


#include "mdb.h"
#include "track.h"
#include "momentumAperture.h"

static SDDS_DATASET SDDSma;
static double momentumOffsetValue = 0;

#define WEDGE_METHOD 1

static void momentumOffsetFunction(double **coord, long np, long pass)
{
  long i;
  if (pass==1)
    for (i=0; i<np; i++) {
      coord[i][0] += x_initial;
      coord[i][2] += y_initial;
      coord[i][5] += momentumOffsetValue;
    }
}

void setupMomentumApertureSearch(
                                 NAMELIST_TEXT *nltext,
                                 RUN *run,
                                 VARY *control
                                 )
{
  char description[200];
  
  /* process namelist input */
  set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
  set_print_namelist_flags(0);
  process_namelist(&momentum_aperture, nltext);
  print_namelist(stdout, &momentum_aperture);

  /* check for data errors */
  if (!output)
    bomb("no output filename specified", NULL);
  if (delta_points<2)
    bomb("delta_points should be at least 2", NULL);
  if (delta_negative_limit>=0)
    bomb("delta_negative_limit should be negative", NULL);
  if (delta_positive_limit<=0)
    bomb("delta_positive_limit should be positive", NULL);
  if (s_start>=s_end)
    bomb("s_start >= s_end", NULL);
  if (include_name_pattern && has_wildcards(include_name_pattern) && strchr(include_name_pattern, '-'))
    include_name_pattern = expand_ranges(include_name_pattern);
  if (n_splits<0)
    bomb("n_splits < 0", NULL);
  
  output = compose_filename(output, run->rootname);
  sprintf(description, "Momentum aperture search");
  if (!SDDS_InitializeOutput(&SDDSma, SDDS_BINARY, 1, description, "momentum aperture",  output)) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    exit(1);
  }
  if (SDDS_DefineColumn(&SDDSma, "ElementName", NULL, NULL, NULL, NULL, SDDS_STRING, 0)<0 ||
      SDDS_DefineColumn(&SDDSma, "s", NULL, "m", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineColumn(&SDDSma, "deltaPositiveFound", NULL, NULL, NULL, NULL, SDDS_LONG, 0)<0 ||
      SDDS_DefineColumn(&SDDSma, "deltaPositive", "$gd$R$bpos$n", NULL, NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineColumn(&SDDSma, "lostOnPassPositive", NULL, NULL, NULL, NULL, SDDS_LONG, 0)<0 ||
      SDDS_DefineColumn(&SDDSma, "sLostPositive", NULL, "m", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineColumn(&SDDSma, "xLostPositive", NULL, "m", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineColumn(&SDDSma, "yLostPositive", NULL, "m", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineColumn(&SDDSma, "deltaLostPositive", NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineColumn(&SDDSma, "deltaNegativeFound", NULL, NULL, NULL, NULL, SDDS_LONG, 0)<0 ||
      SDDS_DefineColumn(&SDDSma, "deltaNegative", "$gd$R$bneg$n", NULL, NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineColumn(&SDDSma, "lostOnPassNegative", NULL, NULL, NULL, NULL, SDDS_LONG, 0)<0 ||
      SDDS_DefineColumn(&SDDSma, "sLostNegative", NULL, "m", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineColumn(&SDDSma, "xLostNegative", NULL, "m", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineColumn(&SDDSma, "yLostNegative", NULL, "m", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineColumn(&SDDSma, "deltaLostNegative", NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0)<0 ||
#if defined(DEBUG)
      SDDS_DefineColumn(&SDDSma, "xco", NULL, "m", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineColumn(&SDDSma, "xpco", NULL, "", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineColumn(&SDDSma, "yco", NULL, "m", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineColumn(&SDDSma, "ypco", NULL, "", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
#endif
      !SDDS_SaveLayout(&SDDSma) || !SDDS_WriteLayout(&SDDSma)) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    exit(1);
  }
  
}

void finishMomentumApertureSearch()
{
  if (SDDS_IsActive(&SDDSma) && !SDDS_Terminate(&SDDSma)) {
    SDDS_SetError("Problem terminating SDDS output (finishMomentumApertureSearch)");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
}


long doMomentumApertureSearch(
                              RUN *run,
                              VARY *control,
                              ERRORVAL *errcon,
                              LINE_LIST *beamline,
                              double *startingCoord
                              )
{    
  double **coord;
  long ip, nElem, iElem;
  double deltaInterval, pCentral;
  ELEMENT_LIST *elem, *elem0;
  long **lostOnPass, **loserFound, **survivorFound, *lostOnPass0, side;
  double deltaStart[2], **deltaSurvived, delta, delta0;
  double **xLost, **yLost, **deltaLost, **sLost;
  double *sStart, length;
  char **ElementName;
  long points, splitsLeft, iTraj, survivorJustSeen, code;
  TRAJECTORY *traj;
#if defined(DEBUG)
  double *xco, *xpco, *yco, *ypco;
  FILE *fpdeb = NULL;
  char s[1000];

  sprintf(s, "%s-Debug", output);
  fprintf(stderr, "Debug file %s\n", s);
  
  if (!(fpdeb = fopen(s, "w")))
    bomb("unable to open debug file for momentum aperture scan", NULL);
  fprintf(fpdeb, "SDDS1\n");
  fprintf(fpdeb, "&parameter name=ElementName type=string &end\n");
  fprintf(fpdeb, "&parameter name=Side type=string &end\n");
  fprintf(fpdeb, "&column name=delta type=double &end\n");
  fprintf(fpdeb, "&column name=Lost type=short &end\n");
  fprintf(fpdeb, "&column name=LossPass type=short &end\n");
  fprintf(fpdeb, "&data mode=ascii no_row_counts=1 &end\n");
  fflush(fpdeb);
#endif

  /* determine how many elements will be tracked */
  elem = &(beamline->elem);
  elem0 = NULL;
  nElem = 0;
  while (elem) {
    length = 0;
    if (entity_description[elem->type].flags&HAS_LENGTH)
      length = ((DRIFT*)(elem->p_elem))->length;
    if ((elem->end_pos-length)>=s_start && (elem->end_pos-length)<=s_end &&
        (!include_name_pattern || wild_match(elem->name, include_name_pattern))) {
      if (!elem0)
	elem0 = elem;
      nElem++;
    } else if ((elem->end_pos-length)>s_end)
      break;
    elem = elem->succ;
  }
  if (nElem==0) 
    SDDS_Bomb("no elements found between s_start and s_end for momentum aperture computation");

  /* allocate arrays for tracking */
  coord = (double**)czarray_2d(sizeof(**coord), 1, 7);

  /* allocate arrays for storing data for negative and positive momentum limits for each element */
  lostOnPass = (long**)czarray_2d(sizeof(**lostOnPass), 2, nElem);
  loserFound = (long**)czarray_2d(sizeof(**loserFound), 2, nElem);
  survivorFound = (long**)czarray_2d(sizeof(**survivorFound), 2, nElem);
  deltaSurvived = (double**)czarray_2d(sizeof(**deltaSurvived), 2, nElem);
  xLost = (double**)czarray_2d(sizeof(**xLost), 2, nElem);
  yLost = (double**)czarray_2d(sizeof(**yLost), 2, nElem);
  deltaLost = (double**)czarray_2d(sizeof(**deltaLost), 2, nElem);
  sLost = (double**)czarray_2d(sizeof(**sLost), 2, nElem);
  sStart = (double*)tmalloc(sizeof(*sStart)*nElem);
  ElementName = (char**)tmalloc(sizeof(*ElementName)*nElem);
#if defined(DEBUG)
  xco = (double*)tmalloc(sizeof(*xco)*nElem);
  xpco = (double*)tmalloc(sizeof(*xpco)*nElem);
  yco = (double*)tmalloc(sizeof(*yco)*nElem);
  ypco = (double*)tmalloc(sizeof(*ypco)*nElem);
#endif

  /* start the output page */
  if (!SDDS_StartTable(&SDDSma, nElem)) {
    SDDS_SetError("Unable to start SDDS table (doMomentumApertureSearch)");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }

  deltaStart[0] = delta_negative_limit;
  deltaStart[1] = delta_positive_limit;

  /* need to do this because do_tracking() in principle may realloc this pointer */
  lostOnPass0 = tmalloc(sizeof(*lostOnPass0)*1);
  
  traj = NULL;
  if (startingCoord) {
    /* track to get the closed orbit vs s */
    traj = tmalloc(sizeof(*traj)*(beamline->n_elems+1));
    memcpy(coord[0], startingCoord, sizeof(**coord)*6);
    memcpy(traj[0].centroid, startingCoord, sizeof(**coord)*6);
    traj[0].elem = NULL;
    delete_phase_references();
    reset_special_elements(beamline, 1);
    pCentral = run->p_central;
    if (!do_tracking(NULL, coord, 1, NULL, beamline, &pCentral, 
                     NULL, NULL, NULL, traj+1, run, control->i_step, 
                     SILENT_RUNNING, 1, 0, NULL, NULL, NULL, lostOnPass0,
                     elem0)) {
      fprintf(stdout, "Error: particle lost while determining central trajectory!\n");
      return 0;
    }
  }
  

  elem = elem0;
  iElem = 0;
  iTraj = -1;

  while (elem) {
    iTraj++;
    if (!include_name_pattern || wild_match(elem->name, include_name_pattern)) {
      length = 0;
      if (entity_description[elem->type].flags&HAS_LENGTH)
        length = ((DRIFT*)(elem->p_elem))->length;
      if ((elem->end_pos-length)>s_end) 
        break;
      if (verbosity>0) {
        fprintf(stdout, "Searching for energy aperture for %s #%ld at s=%em\n", elem->name, elem->occurence, elem->end_pos-length);
        fflush(stdout);
      }
      ElementName[iElem] = elem->name;
      sStart[iElem] = elem->end_pos-length;
      if (startingCoord && iTraj && traj[iTraj].elem!=elem->pred) {
        fprintf(stdout, "Fatal error: book-keeping error for central trajectory.\n");
        fprintf(stdout, "iTraj = %ld, traj.elem = %p, traj.elem->pred = %p, traj.elem->succ = %p, elem = %p\n",
                iTraj, traj[iTraj].elem, 
                traj[iTraj].elem->pred, 
                traj[iTraj].elem->succ, elem);
        exit(1);
      }
      for (side=0; side<2; side++) {
#if defined(DEBUG)
        fprintf(fpdeb, "%s\n%s\n",
                elem->name, side==0?"negative":"positive");
        fflush(fpdeb);
#endif
        survivorJustSeen = 0;
        deltaInterval = deltaStart[side]/(delta_points-1);
        lostOnPass[side][iElem] = -1;
        loserFound[side][iElem] = survivorFound[side][iElem] = 0;
        xLost[side][iElem] = yLost[side][iElem] = 
          deltaLost[side][iElem] = sLost[side][iElem] = 
            deltaSurvived[side][iElem] =  0;
        if (verbosity>1) {
          fprintf(stdout, " Searching for %s side from 0 to %e with interval %e\n", side==0?"negative":"positive",
                  deltaStart[side], deltaInterval);
          fflush(stdout);
        }
        points = delta_points;
        splitsLeft = n_splits;
        delta0 = 0;
        if (scan_from_outside) {
          delta0 = deltaStart[side];
          deltaInterval *= -1;
        }
        do {
          for (ip=0; ip<points; ip++) {
            delta = delta0 + ip*deltaInterval;
            delete_phase_references();
            reset_special_elements(beamline, 1);
            pCentral = run->p_central;
#if WEDGE_METHOD
            setTrackingWedgeFunction(momentumOffsetFunction, elem); 
            momentumOffsetValue = delta;
            if (startingCoord)
              memcpy(coord[0], startingCoord, sizeof(double)*6);
            else
              memset(coord[0], 0, sizeof(**coord)*6);
#else
            if (startingCoord)
              memcpy(coord[0], traj[iTraj].centroid, sizeof(double)*6);
            else
              memset(coord[0], 0, sizeof(**coord)*6);
            coord[0][5] += delta;
            coord[0][0] += x_initial;
            coord[0][2] += y_initial;
#endif
            coord[0][6] = 1;
#if defined(DEBUG)
            xco[iElem] = coord[0][0];
            xpco[iElem] = coord[0][1];
            yco[iElem] = coord[0][2];
            ypco[iElem] = coord[0][3];
#endif
            if (verbosity>3) {
              fprintf(stdout, "  Tracking with delta0 = %e (%e, %e, %e, %e, %e, %e)\n", 
                      delta, coord[0][0], coord[0][1], coord[0][2], coord[0][3], coord[0][4], coord[0][5]);
              fflush(stdout);
            }
            lostOnPass0[0] = -1;
            code = do_tracking(NULL, coord, 1, NULL, beamline, &pCentral, 
                               NULL, NULL, NULL, NULL, run, control->i_step, 
                               SILENT_RUNNING, control->n_passes, 0, NULL, NULL, NULL, lostOnPass0,
                               WEDGE_METHOD?NULL:elem);
#if defined(DEBUG)
            fprintf(fpdeb, "%21.15e %hd %hd\n", delta, code?(short)0:(short)1, (short)lostOnPass0[0]);
#endif
            if (code) {
              /* particle survived */
              if (verbosity>2) {
                fprintf(stdout, "  Particle survived with delta0 = %e, sf = %e m\n", delta, coord[0][4]);
                fflush(stdout);
              }
              if (!survivorFound[side][iElem] || fabs(delta)>fabs(deltaSurvived[side][iElem]))  {
                if (verbosity>4)
                  fprintf(stdout, "  New survivor-of-record for this scan\n");
                deltaSurvived[side][iElem] = delta;
              }
              survivorFound[side][iElem] = 1;
              if (scan_from_outside) {
                if (survivorJustSeen)
                  break;
                survivorJustSeen = 1;
              }
            } else {
              /* particle lost */
              if (verbosity>3) {
                long i;
                fprintf(stdout, "  Particle lost with delta0 = %e\n", delta);
                if (verbosity>4)
                  for (i=0; i<6; i++)
                    fprintf(stdout, "   coord[%ld] = %e\n", i, coord[0][i]);
                fflush(stdout);
              }
              lostOnPass[side][iElem] = lostOnPass0[0];
              xLost[side][iElem] = coord[0][0];
              yLost[side][iElem] = coord[0][2];
              sLost[side][iElem] = coord[0][4];
              deltaLost[side][iElem] = (coord[0][5]-pCentral)/pCentral;
              loserFound[side][iElem] = 1;
              if (!scan_from_outside)
                break;
              survivorJustSeen = 0;
            }
          }
          points = 1;
          deltaInterval /= 2;
          if (survivorFound[side][iElem]) 
            delta0 = deltaSurvived[side][iElem] + (scan_from_outside?-1:1)*deltaInterval;
          else
            break;
        } while (--splitsLeft >= 0);
#if defined(DEBUG)
        fprintf(fpdeb, "\n");
        fflush(fpdeb);
#endif
      }
      iElem++;
    }
    elem = elem->succ;
  }    

  if (!SDDS_StartPage(&SDDSma, nElem) ||
      !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, ElementName, nElem, "ElementName") ||
      !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, sStart, nElem, "s") ||
      !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, loserFound[1], nElem, "deltaPositiveFound") ||
      !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, deltaSurvived[1], nElem, "deltaPositive") ||
      !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, lostOnPass[1], nElem, "lostOnPassPositive") ||
      !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, sLost[1], nElem, "sLostPositive") ||
      !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, xLost[1], nElem, "xLostPositive") ||
      !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, yLost[1], nElem, "yLostPositive") ||
      !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, deltaLost[1], nElem, "deltaLostPositive") ||
      !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, loserFound[0], nElem, "deltaNegativeFound") ||
      !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, deltaSurvived[0], nElem, "deltaNegative") ||
      !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, lostOnPass[0], nElem, "lostOnPassNegative") ||
      !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, sLost[0], nElem, "sLostNegative") ||
      !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, xLost[0], nElem, "xLostNegative") ||
      !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, yLost[0], nElem, "yLostNegative") ||
      !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, deltaLost[0], nElem, "deltaLostNegative") ||
#if defined(DEBUG)
      !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, xco, nElem, "xco") ||
      !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, xpco, nElem, "xpco") ||
      !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, yco, nElem, "yco") ||
      !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, ypco, nElem, "ypco") ||
#endif
      !SDDS_WritePage(&SDDSma)) {
    SDDS_SetError("Problem writing SDDS table (doMomentumApertureSearch)");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  SDDS_DoFSync(&SDDSma);
      

  free_czarray_2d((void**)coord, 1, 7);
  free_czarray_2d((void**)lostOnPass, 2, nElem);
  free_czarray_2d((void**)loserFound, 2, nElem);
  free_czarray_2d((void**)survivorFound, 2, nElem);
  free_czarray_2d((void**)deltaSurvived, 2, nElem);
  free_czarray_2d((void**)xLost, 2, nElem);
  free_czarray_2d((void**)yLost, 2, nElem);
  free_czarray_2d((void**)deltaLost, 2, nElem);
  free_czarray_2d((void**)sLost, 2, nElem);
  free(sStart);
  free(ElementName);
  free(lostOnPass0);
  if (traj)
    free(traj);
  
  return 1;
}


