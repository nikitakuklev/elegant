/* Copyright 1994 by Michael Borland and Argonne National Laboratory,
 * all rights reserved.
 */
/* file: amplif.c
 * purpose: compute amplification factors 
 *
 * Michael Borland, 1992
 */
#include "mdb.h"
#include "track.h"
#include "match_string.h"
#include "correctDefs.h"

void compute_amplification_factors(
                                   NAMELIST_TEXT *nltext,
                                   RUN *run,
                                   CORRECTION *correct,
                                   long do_closed_orbit,
                                   LINE_LIST *beamline
                                   )
{
#include "amplif.h"
  ELEMENT_LIST *eptr, *emax, *emaxu, *emaxc;
  long type_code, iparam, iplane, n_part, i, j, nsum, n_kicks;
  double start[7], p, max_pos, max_kick, max_Ac, max_Au, rms_pos, rms_kick, original_value;
  static double **one_part = NULL;
  static TRAJECTORY *traj = NULL, *trajc = NULL;
  static long n_traj = 0;
  static FILE *fpout = NULL, *fpcof = NULL, *fpuof = NULL, *fpkf = NULL;
  double *kick, *Cij;
  static char name_header[256], unit_header[256], printf_string[256], unit_pos[32], unit_kick[32], *Ai_unit, *Cij_unit;
  static char s[256], description[256];
  double *Ac_vs_z, *Au_vs_z;
  CORMON_DATA *CM;
  STEERING_LIST *SL;

  log_entry("compute_amplification_factors");

  if (!one_part)
    one_part = (double**)zarray_2d(sizeof(**one_part), 1, 7);
  Au_vs_z = tmalloc(sizeof(*Au_vs_z)*(n_traj=beamline->n_elems+1));
  Ac_vs_z = tmalloc(sizeof(*Ac_vs_z)*(n_traj=beamline->n_elems+1));

  /* process namelist input */
  set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
  set_print_namelist_flags(0);
  process_namelist(&amplification_factors, nltext);
  print_namelist(stderr, &amplification_factors);

  if (plane)
    str_tolower(plane);
  if (plane && plane[0]=='v')
    iplane = 2;
  else
    iplane = 0;
  n_kicks = 0;
  if (correct->mode!=-1) {
    CM = iplane==0?correct->CMx:correct->CMy;
    SL = iplane==0?&correct->SLx:&correct->SLy;
    Cij = tmalloc(sizeof(*Cij)*(n_kicks=CM->ncor));
  }
  else
    Cij = NULL;

  if (type) {
    str_toupper(type);
    if ((type_code=match_string(type, entity_name, N_TYPES, EXACT_MATCH))<0) {
      fprintf(stderr, "warning: no exact match for type %s\n", type);
      if ((type_code=match_string(type, entity_name, N_TYPES, DCL_STYLE_MATCH|RETURN_FIRST_MATCH))<0)
        bomb("unknown type name", NULL);
      fprintf(stderr, "assuming that you meant type %s\n", entity_name[type_code]);
      cp_str(&type, entity_name[type_code]);
    }
  }
  else
    type_code = -1;
  if (name)
    str_toupper(name);

  if (!type && !name)
    bomb("at least one of type and name must be defined", NULL);
  if (!item)
    bomb("item must be given", NULL);
  str_toupper(item);
  if (type_code!=-1 && !confirm_parameter(item, type_code))
    bomb("specified element type does not have specified item as a parameter", NULL);

  if (fpout)
    fclose(fpout);
  if (output)
    fpout = fopen_e(output=compose_filename(output, run->rootname), "w", 0);
  else
    fpout = NULL;
  if (corrected_orbit_function)
    fpcof = fopen_e(corrected_orbit_function = compose_filename(corrected_orbit_function, run->rootname), "w", 
                    FOPEN_SAVE_IF_EXISTS);
  if (uncorrected_orbit_function)
    fpuof = fopen_e(uncorrected_orbit_function = compose_filename(uncorrected_orbit_function, run->rootname), "w", 
                    FOPEN_SAVE_IF_EXISTS);
  if (kick_function)
    fpkf = fopen_e(kick_function = compose_filename(kick_function, run->rootname), "w", FOPEN_SAVE_IF_EXISTS);

  if (number_to_do>0 && maximum_z>0)
    bomb("only one of number_to_do and maximum_z may be positive", NULL);

  if (type_code==-1) {
    eptr = NULL;
    iparam = -1;
    while ((eptr = wfind_element(name, &eptr, &(beamline->elem)))) {
      if ((iparam=confirm_parameter(item, eptr->type))>=0)
        break;
    }
    if (iparam<0 || !eptr)
      bomb("no element exists of given name with given item as a parameter", NULL);
  }
  else {
    eptr = &(beamline->elem);
    while (eptr) {
      if (eptr->type==type_code) {
        if ((iparam=confirm_parameter(item, eptr->type))>=0)
          break;
        bomb("elements of the given type do not have the given item as a parameter", NULL);
      }
      eptr = eptr->succ;
    }
    if (!eptr)
      bomb("no items of the given type exist in the beamline", NULL);
  }

  /* prepare labels */
  name_header[0] = unit_header[0] = printf_string[0] = 0;
  sprintf(unit_pos, "m/%s", entity_description[eptr->type].parameter[iparam].unit);
  cp_str(&Ai_unit, unit_pos);
  add_to_standard_headers(name_header, unit_header, printf_string,  "rms pos.", unit_pos, "%13.6le", 15);
  add_to_standard_headers(name_header, unit_header, printf_string,  "max pos.", unit_pos, "%13.6le", 15);
  add_to_standard_headers(name_header, unit_header, printf_string,     "z_max", "m", "%13.6le", 15);
  if (correct->mode!=-1) {
    sprintf(unit_kick, "rad/%s", entity_description[eptr->type].parameter[iparam].unit);
    cp_str(&Cij_unit, unit_kick);
    add_to_standard_headers(name_header, unit_header, printf_string,  "rms kick", unit_kick, "%13.6le", 15);
    add_to_standard_headers(name_header, unit_header, printf_string,  "max kick", unit_kick, "%13.6le", 15);
  }
  add_to_standard_headers(name_header, unit_header, printf_string,    "z_elem", "m", "%13.6le", 15);
  add_to_standard_headers(name_header, unit_header, printf_string,   "element",  "", "%10s#%ld", 10);

  if (correct->mode!=-1) {
    fputs("Amplification factors with correction for", stderr);
    if (fpuof) {
      if (correct->n_xy_cycles!=1)
        bomb("n_xy_cycles!=1--can't provide uncorrected amplification function\nuse separate run or set n_xy_cycles=1", NULL);
    }
  }
  else {
    fputs("Amplification factors", stderr);
    if (fpcof)
      bomb("can't compute corrected-orbit amplification function if you don't do correction", NULL);
    if (fpkf)
      bomb("can't compute kick function if you don't do orbit correction", NULL);
  }

  description[0] = 0;
  if (number_to_do>0) {
    sprintf(description, "First %ld elements", number_to_do);
  }
  else if (maximum_z>0) {
    sprintf(description, "Elements with z<%em", maximum_z);
  }
  else {
    sprintf(description, "All elements");
  }
  if (type_code!=-1) {
    sprintf(s, " of type %s", type);
    strcat(description, s);
  }
  if (name) {
    sprintf(s, " named %s", name);
    strcat(description, s);
  }

  sprintf(s, ", when %s is changed (by %g %s)",
          item, change, entity_description[eptr->type].parameter[iparam].unit);
  strcat(description, s);

  fputs(description, stderr);

  if (fpout) {
    fprintf(fpout, "SDDS1\n&description text=\"Amplification factors for beamline %s from %s\" &end\n",
            beamline->name, run->lattice);
    fprintf(fpout, "&parameter name=WithCorrection, type=string, fixed_value=%s &end\n",
            correct->mode!=-1?"yes":"no");
    fprintf(fpout, "&parameter name=GroupDescription, type=string, fixed_value=\"%s\" &end\n",
            description);
    fprintf(fpout, "&column name=%sResponseRMS, units=%s, type=double &end\n", iplane?"x":"y", unit_pos);
    fprintf(fpout, "&column name=%sResponseMaximum, units=%s, type=double &end\n", iplane?"x":"y", unit_pos);
    fprintf(fpout, "&column name=sMaximum, units=m, type=double &end\n");
    fprintf(fpout, "&column name=KickRMS, units=%s, type=double &end\n", unit_kick);
    fprintf(fpout, "&column name=KickMaximum, units=%s, type=double &end\n", unit_kick);
    fprintf(fpout, "&column name=sActuator, units=m, type=double &end\n");
    fprintf(fpout, "&column name=FullElementName, type=string &end\n");
    fprintf(fpout, "&data mode=ascii, no_row_counts=1 &end\n");
  }
  if (fpuof) {
    fprintf(fpuof, "SDDS1\n&description text=\"Uncorrected amplification functions for beamline %s from %s\" &end\n",
            beamline->name, run->lattice);
    fprintf(fpuof, "&parameter name=GroupDescription, type=string, fixed_value=\"%s\" &end\n",
            description);
    fprintf(fpuof, "&parameter name=Actuator, type=string &end\n");
    fprintf(fpuof, "&column name=s, units=m, type=double &end\n");
    fprintf(fpuof, "&column name=%sResponse, units=%s, type=double &end\n", iplane==0?"x":"y", unit_pos);
    fprintf(fpuof, "&column name=ElementName, type=string &end\n&column name=ElementOccurence, type=long &end\n");
    fprintf(fpuof, "&data mode=ascii &end\n");
  }
  if (fpcof) {
    fprintf(fpcof, "SDDS1\n&description text=\"Corrected amplification functions for beamline %s from %s\" &end\n",
            beamline->name, run->lattice);
    fprintf(fpcof, "&parameter name=GroupDescription, type=string, fixed_value=\"%s\" &end\n",
            description);
    fprintf(fpcof, "&parameter name=Actuator, type=string &end\n");
    fprintf(fpcof, "&column name=s, units=m, type=double &end\n");
    fprintf(fpcof, "&column name=%sResponse, units=%s, type=double &end\n", iplane==0?"x":"y", unit_pos);
    fprintf(fpcof, "&column name=ElementName, type=string &end\n&column name=ElementOccurence, type=long &end\n");
    fprintf(fpcof, "&data mode=ascii &end\n");
  }
  if (fpkf) {
    fprintf(fpkf, "SDDS1\n&description text=\"Correction kick amplification functions for beamline %s from %s\" &end\n",
            beamline->name, run->lattice);
    fprintf(fpkf, "&parameter name=GroupDescription, type=string, fixed_value=\"%s\" &end\n",
            description);
    fprintf(fpkf, "&column name=s, units=m, type=double &end\n");
    fprintf(fpkf, "&column name=C, units=%s, type=double &end\n", Cij_unit);
    fprintf(fpkf, "&column name=ElementName, type=string &end\n&column name=ElementOccurence, type=long &end\n");
    fprintf(fpkf, "&data mode=ascii, no_row_counts=1 &end\n");
  }

  fputs(name_header, stderr);
  fputc('\n', stderr);
  fputs(unit_header, stderr);
  fputc('\n', stderr);
  fflush(stderr);

  if (!name)
    eptr = &(beamline->elem);
  else
    eptr = NULL;
  while (number_to_do && ((name && (eptr=wfind_element(name, &eptr, &(beamline->elem)))) || (!name && eptr)) &&
         !(maximum_z && eptr->end_pos>maximum_z)) {
    if (type_code!=-1 && (eptr->type!=type_code)) {
      if (!name)
        eptr = eptr->succ;
      continue;
    }
    if ((iparam=confirm_parameter(item, eptr->type))==-1)
      continue;
    number_to_do--;
    fprintf(stderr, "\nWorking on element %s#%ld at z=%em\n", eptr->name, eptr->occurence, eptr->end_pos);
    if (entity_description[eptr->type].parameter[iparam].type!=IS_DOUBLE)
      bomb("item is not floating-point type", NULL);

    original_value = *((double*)(eptr->p_elem+entity_description[eptr->type].parameter[iparam].offset));
    *((double*)(eptr->p_elem+entity_description[eptr->type].parameter[iparam].offset)) += change;
    if (entity_description[eptr->type].parameter[iparam].changes_matrix) {
      if (eptr->matrix) {
        free_matrices(eptr->matrix);
        tfree(eptr->matrix);
        eptr->matrix = NULL;
      }
      compute_matrix(eptr, run, NULL);
    }

    if (correct->mode!=-1) {
      if (!do_correction(correct, run, beamline, NULL, NULL, 0, 1))
        bomb("correction failed", NULL);
      traj  = correct->traj[0];
      trajc = correct->traj[2];
    }
    else if (do_closed_orbit) {
      run_closed_orbit(run, beamline, start, NULL, 1);
      traj = beamline->closed_orbit;
    }
    else {
      if (n_traj && n_traj<(beamline->n_elems+1)) {
        if (traj)
          tfree(traj);
        traj = NULL;
        n_traj = beamline->n_elems+1;
      }
      if (!traj) 
        traj = tmalloc(sizeof(*traj)*(beamline->n_elems+1));
      n_part = 1;
      p = sqrt(sqr(run->ideal_gamma)-1);
      fill_double_array(*one_part, 7, 0.0);
      if (!do_tracking(one_part, &n_part, NULL, beamline, &p, (double**)NULL, (BEAM_SUMS**)NULL, (long*)NULL,
                       traj, run, 0, TEST_PARTICLES+TIME_DEPENDENCE_OFF, 1, NULL))
        bomb("tracking failed for test particle (compute_amplification_factors)", NULL);
    }

    *((double*)(eptr->p_elem+entity_description[eptr->type].parameter[iparam].offset)) = original_value;
    if (entity_description[eptr->type].parameter[iparam].changes_matrix) {
      if (eptr->matrix) {
        free_matrices(eptr->matrix);
        tfree(eptr->matrix);
        eptr->matrix = NULL;
      }
      compute_matrix(eptr, run, NULL);
    }

    max_pos = -DBL_MAX;
    emax = NULL;
    for (i=nsum=rms_pos=0; i<=beamline->n_elems; i++) {
      if (!traj[i].n_part)
        continue;
      Au_vs_z[i] += sqr(traj[i].centroid[iplane]);
      if (trajc) {
        if (fabs(trajc[i].centroid[iplane])>max_pos) {
          max_pos = fabs(trajc[i].centroid[iplane]);
          emax = trajc[i].elem;
        }
        rms_pos += sqr(trajc[i].centroid[iplane]);
        Ac_vs_z[i] += sqr(trajc[i].centroid[iplane]);
      }
      else {
        if (fabs(traj[i].centroid[iplane])>max_pos) {
          max_pos = fabs(traj[i].centroid[iplane]);
          emax = traj[i].elem;
        }
        rms_pos += sqr(traj[i].centroid[iplane]);
      }
      nsum++;
    }
    if (nsum)
      rms_pos = sqrt(rms_pos/nsum);
    if (fpuof) {
      fprintf(fpuof, "%s#%ld\n%ld\n", eptr->name, eptr->occurence, n_traj-2);
      for (i=1; i<n_traj-1; i++)
        if (traj[i].elem)
          fprintf(fpuof, "%e %e %s %ld\n", traj[i].elem->end_pos, traj[i].centroid[iplane]/change, 
                  traj[i].elem->name, traj[i].elem->occurence);
        else
          fprintf(fpuof, "0 0 ? 0\n");
    }
    
    if (fpcof) {
      fprintf(fpcof, "%s#%ld\n%ld\n", eptr->name, eptr->occurence, n_traj-2);
      for (i=1; i<n_traj-1; i++) 
        if (trajc[i].elem)
          fprintf(fpcof, "%e %e %s %ld\n", trajc[i].elem->end_pos, trajc[i].centroid[iplane]/change, 
                  trajc[i].elem->name, trajc[i].elem->occurence);
        else
          fprintf(fpcof, "0 0 ? 0\n");
    }
    if (correct->mode==-1 && emax) {
      fprintf(stderr, printf_string,
              rms_pos/change, max_pos/change, emax->end_pos, eptr->end_pos, eptr->name, eptr->occurence);
      fputc('\n', stderr);
      if (fpout) {
        fprintf(fpout, printf_string,
                rms_pos/change, max_pos/change, emax->end_pos, eptr->end_pos, eptr->name, eptr->occurence);
        fputc('\n', fpout);
      }
    }
    else {
      kick = CM->kick[correct->n_iterations];
      max_kick = -DBL_MAX;
      for (i=rms_kick=0; i<n_kicks; i++) {
        Cij[i] += sqr(kick[i]);
        if (fabs(kick[i])>max_kick)
          max_kick = fabs(kick[i]);
        rms_kick += sqr(kick[i]);
      }
      if (n_kicks)
        rms_kick = sqrt(rms_kick/n_kicks);
      fprintf(stderr, printf_string,
              rms_pos/change, max_pos/change, emax->end_pos, rms_kick/change, 
              max_kick/change, eptr->end_pos, eptr->name, eptr->occurence);
      fputc('\n', stderr);
      if (fpout) {
        fprintf(fpout, printf_string,
                rms_pos/change, max_pos/change, emax->end_pos, rms_kick/change, 
                max_kick/change, eptr->end_pos, eptr->name, eptr->occurence);
        fputc('\n', fpout);
      }
    }

    if (!name)
      eptr = eptr->succ;
  }

  if (fpuof) 
    fprintf(fpuof, "ResponseRMS\n%ld\n", n_traj-2);
  if (fpcof) 
    fprintf(fpcof, "ResponseRMS\n%ld\n", n_traj-2);
  if (nsum) {
    max_Ac = 0;
    max_Au = 0;
    emaxc = NULL;
    emaxu = NULL;
    for (i=1; i<n_traj-1; i++) {
      if (!traj[i].elem) {
        if (fpuof)
          fprintf(fpuof, "0 0 ? 0\n");
        if (fpcof)
          fprintf(fpcof, "0 0 ? 0\n");
        continue;
      }
      Au_vs_z[i] = sqrt(Au_vs_z[i])/change;
      Ac_vs_z[i] = sqrt(Ac_vs_z[i])/change;
      if (fpuof && i<n_traj-1)
        fprintf(fpuof, "%e %e %s %ld\n", traj[i].elem->end_pos, Au_vs_z[i],
                traj[i].elem->name, traj[i].elem->occurence);
      if (fpcof && i<n_traj-1)
        fprintf(fpcof, "%e %e %s %ld\n", traj[i].elem->end_pos, Ac_vs_z[i], 
                traj[i].elem->name, traj[i].elem->occurence);
      if (Ac_vs_z[i]>max_Ac) {
        max_Ac = Ac_vs_z[i];
        emaxc = trajc[i].elem;
      }
      if (Au_vs_z[i]>max_Au) {
        max_Au = Au_vs_z[i];
        emaxu = traj[i].elem;
      }
    }
    fputc('\n', stderr);
    if (emaxc && correct->mode!=-1)
      fprintf(stderr, "maximum corrected-orbit amplification function is %e %s at %s at z=%em\n",
              max_Ac, Ai_unit, emaxc->name, emaxc->end_pos);
    if (emaxu)
      fprintf(stderr, "maximum orbit amplification function is %e %s at %s at z=%em\n",
              max_Au, Ai_unit, emax->name, emax->end_pos);
  }
  if (fpcof)
    fclose(fpcof);
  if (fpuof)
    fclose(fpuof);
  tfree(Au_vs_z);
  tfree(Ac_vs_z);

  if (n_kicks && fpkf) {
    for (i=max_kick=0; i<n_kicks; i++) {
      Cij[i] = sqrt(Cij[i])/change;
      if (max_kick<Cij[i]) {
        max_kick=Cij[i];
        j = i;
      }
      fprintf(fpkf, "%e %e %s %ld\n", CM->ucorr[i]->end_pos, Cij[i],
              CM->ucorr[i]->name, CM->ucorr[i]->occurence);
    }
    fclose(fpkf);
    if (max_kick)
      fprintf(stderr, "maximum kick amplification factor is %e %s from %s#%ld.%s at %e m\n",
              max_kick, Cij_unit, CM->ucorr[j]->name, CM->ucorr[j]->occurence,
              SL->corr_param[CM->sl_index[j]], CM->ucorr[j]->end_pos);
  }
  if (Cij)
    tfree(Cij);

  log_exit("compute_amplification_factors");
}
