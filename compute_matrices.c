/* Copyright 1994 by Michael Borland and Argonne National Laboratory,
 * all rights reserved.
 */
/* contents: compute_matrices(), drift_matrix(), sextupole_matrix() 
 *           plus more...
 * Michael Borland, 1989.
 */
#include "track.h"
#include "mdb.h"

#define DEBUG 0

long determine_bend_flags(ELEMENT_LIST *eptr, long edge1_effects, long edge2_effects);

VMATRIX *full_matrix(ELEMENT_LIST *elem, RUN *run, long order) 
{
    VMATRIX *M1, *M2, *tmp, *Md;
    ELEMENT_LIST *member;
    long i;
    double Pref_input;

    if (!elem) {
        fputs("error: NULL element pointer passed to full_matrix", stderr);
        abort();
        }
    
#ifdef WATCH_MEMORY
    fprintf(stderr, "start full_matrix: CPU: %6.2lf  PF: %6ld  MEM: %6ld\n",
           cpu_time()/100.0, page_faults(), memory_count());
#endif
    
    return accumulate_matrices(elem, run, NULL, order, 1);
    }

VMATRIX *accumulate_matrices(ELEMENT_LIST *elem, RUN *run, VMATRIX *M0, long order, long full_matrix_only)
{
  VMATRIX *M1, *M2, *tmp;
  ELEMENT_LIST *member;
  double Pref_input;
  long i;
  
  if (!elem) {
    fputs("error: NULL element pointer passed to accumulate_matrices", stderr);
    abort();
  }
  
  initialize_matrices(M1=tmalloc(sizeof(*M1)), order);
  initialize_matrices(M2=tmalloc(sizeof(*M2)), order);
  for (i=0; i<6; i++)
    M1->R[i][i] = M2->R[i][i] = 1;
  if (M0)
    copy_matrices1(M1, M0);
  
  member = elem;
  while (member) {
    if (member->type<0 || member->type>=N_TYPES) {
      fprintf(stderr, "error: bad element type %ld (accumulate_matrices)\n", member->type);
      fprintf(stderr, "element name is %s and end position is %em\n", 
             (member->name?member->name:"{null}"), member->end_pos);
      abort();
    }
    if (member->pred)
      Pref_input = member->pred->Pref_output;
    else
      Pref_input = member->Pref_input;
    if (!member->matrix || Pref_input!=member->Pref_input)
      compute_matrix(member, run, NULL);
    if ((entity_description[member->type].flags&HAS_MATRIX) && !member->matrix) {
      fprintf(stderr, "programming error: matrix not computed for element %s\n",
              member->name);
      abort();
    }
    if (member->matrix) {
      concat_matrices(M2, member->matrix, M1);
      tmp = M2;
      M2  = M1;
      M1  = tmp;
      if (!full_matrix_only) {
        if (member->accumMatrix)
          free_matrices(member->accumMatrix);
        else 
          member->accumMatrix = tmalloc(sizeof(*(member->accumMatrix)));
        copy_matrices(member->accumMatrix, M2);
      }
    }
    member = member->succ;
  }
  if (M2) {
    free_matrices(M2); tfree(M2); M2 = NULL;
  }
  return M1;
}

VMATRIX *append_full_matrix(ELEMENT_LIST *elem, RUN *run, VMATRIX *M0, long order) 
{
  return accumulate_matrices(elem, run, M0, order, 1);
#if 0  
  VMATRIX *M1, *M2, *tmp, *Md;
    ELEMENT_LIST *member;
    double Pref_input;

    log_entry("append_full_matrix");
    if (!elem) {
        fputs("error: NULL element pointer passed to append_full_matrix", stderr);
        abort();
        }
    if (!M0) {
        fputs("error: NULL initial matrix pointer passed to append_full_matrix", stderr);
        abort();
        }
    
    initialize_matrices(M1=tmalloc(sizeof(*M1)), order);
    initialize_matrices(M2=tmalloc(sizeof(*M2)), order);
    copy_matrices1(M1, M0);
    
    member = elem;
    while (member) {
        if (member->type<0 || member->type>=N_TYPES) {
            fprintf(stderr, "error: bad element type %ld (full_matrix)\n", member->type);
            fprintf(stderr, "element name is %s and end position is %em\n", 
                   (member->name?member->name:"{null}"), member->end_pos);
            abort();
            }
        if (member->pred)
            Pref_input = member->pred->Pref_output;
        else
            Pref_input = member->Pref_input;
        if (!member->matrix || Pref_input!=member->Pref_input)
            compute_matrix(member, run, NULL);
        if ((entity_description[member->type].flags&HAS_MATRIX) && !member->matrix) {
            fprintf(stderr, "programming error: matrix not computed for element %s\n",
                    member->name);
            abort();
            }
        if (member->matrix) {
            concat_matrices(M2, member->matrix, M1);
            tmp = M2;
            M2  = M1;
            M1  = tmp;
            }
        member = member->succ;
        }
    if (M2) {
        free_matrices(M2); tfree(M2); M2 = NULL;
        }
    log_exit("append_full_matrix");
    return(M1);
#endif
    }

long fill_in_matrices(
                      ELEMENT_LIST *elem,
                      RUN *run
                      )
{
    ELEMENT_LIST *member;
    long n_elements;
    
    log_entry("fill_in_matrices");
    
    n_elements = 0;
    member = elem;
    while (member) {
        if (member->type<0 || member->type>=N_TYPES) {
            fprintf(stderr, "error: bad element type %ld (fill_in_matrices)\n", member->type);
            fprintf(stderr, "element name is %s and end position is %em\n", 
                   (member->name?member->name:"{null}"), member->end_pos);
            abort();
            }
        if ((member->matrix==NULL || (member->pred && member->pred->Pref_output!=member->Pref_input)) &&
            entity_description[member->type].flags&HAS_MATRIX) {
          compute_matrix(member, run, NULL);
          n_elements++;
        }
        member = member->succ;
        }
    log_exit("fill_in_matrices");    
    return(n_elements);
    }

long calculate_matrices(
                        LINE_LIST *line,            /* Beamline to calculate matrices for. */
                        RUN *run
                        )
{
    ELEMENT_LIST *member;
    long n_elements;
    
    log_entry("calculate_matrices");
    
    n_elements = 0;
    member = &(line->elem);
    line->elem_recirc = NULL;
    line->i_recirc    = 0;
    while (member) {
        if (!member->matrix || (member->pred && member->pred->Pref_output!=member->Pref_input))
            compute_matrix(member, run, NULL);
        if (member->type==T_RECIRC) {
            line->elem_recirc = member;
            line->i_recirc    = n_elements;
            }
        n_elements++;
        member = member->succ;
        }
    log_exit("calculate_matrices");    
    return(n_elements);
    }

VMATRIX *drift_matrix(double length, long order)
{
    VMATRIX *M;
    double *C, **R, ***T;
    
    log_entry("drift_matrix");
    
    M = tmalloc(sizeof(*M));
    M->order = (order>2?2:order);
    initialize_matrices(M, M->order);
    R = M->R;
    C = M->C;
    T = M->T;
    
    C[4] = length;
    R[0][0] = R[1][1] = R[2][2] = R[3][3] = R[4][4] = R[5][5] = 1;
    R[0][1] = R[2][3] = length;
    
    if (order>1)
        T[4][1][1] = T[4][3][3] = length/2;
    
    log_exit("drift_matrix");
    return(M);
    }

VMATRIX *sextupole_matrix(double K2, double length, long maximum_order, double tilt, double fse)
{
    VMATRIX *M;
    double *C, **R, ***T;
    double temp;
    
    log_entry("sextupole_matrix");
    
    K2 *= (1+fse);

    M = tmalloc(sizeof(*M));
    initialize_matrices(M, M->order=MIN(2,maximum_order));
    R = M->R;
    C = M->C;
    
    R[0][0] = R[1][1] = R[2][2] = R[3][3] = R[4][4] = R[5][5] = 1;
    C[4] = R[0][1] = R[2][3] = length;
    if (M->order>=2) {
        T = M->T;
        temp = K2*length/2;   /* temp = ks^2*l */
        T[1][0][0] = -(T[1][2][2] = temp);
        T[3][2][0] = 2*temp;
        temp *= length;       /* temp = ks^2*l^2 */
        T[0][0][0] = -temp/2;
        T[0][2][2] = temp/2;
        T[1][1][0] = -temp;
        T[1][3][2] = temp;
        T[2][2][0] =  T[3][3][0] = T[3][2][1] = temp;
        temp *= length;       /* temp = ks^2*l^3 */
        T[0][1][0] = -temp/3.;
        T[1][1][1] = -temp/3.;
        T[0][3][2] = T[2][3][0] = T[2][2][1] = temp/3.;
        T[1][3][3] = temp/3.;
        T[3][3][1] = 2*temp/3;
        temp *= length;       /* temp = ks^2*l^4 */
        T[0][1][1] = -temp/12;
        T[0][3][3] =  temp/12;
        T[2][3][1] =  temp/6;
        /* path length terms--same as for drift */
        T[4][1][1] = T[4][3][3] = length/2;
        }
    tilt_matrices(M, tilt);
    log_exit("sextupole_matrix");
    return(M);
    }


VMATRIX *solenoid_matrix(double length, double ks, long max_order)
{
    VMATRIX *M;
    double *C, **R, ***T;
    double cos_ksl, sin_ksl, ksl, CS, S2, temp;
    
    log_entry("solenoid_matrix");
    
    if (ks==0 || length==0)
        return(drift_matrix(length, max_order));
    
    /* defined Ks = -B/(B.rho) as in MAD, which is different from TRANSPORT Ks = -B/(2*B.rho) */
    ks /= 2;
    ksl = ks*length;

    M = tmalloc(sizeof(*M));
    M->order = MIN(max_order, 2);
    initialize_matrices(M, M->order);
    R = M->R;
    C = M->C;
    T = M->T;
    
    cos_ksl = cos(ksl);
    sin_ksl = sin(ksl);
    C[4] = length;
    
    R[0][0] = R[1][1] = R[2][2] = R[3][3] = sqr(cos_ksl);
    R[2][0] = R[3][1] = -(R[0][2] = R[1][3] = CS = sin_ksl*cos_ksl);
    R[4][4] = R[5][5] = 1;
    
    R[0][1] = R[2][3] = CS/ks;
    R[1][0] = R[3][2] = -ks*CS;
    
    S2 = sqr(sin_ksl);
    R[1][2] = -(R[3][0] = ks*S2);
    R[2][1] = -(R[0][3] = S2/ks);
    
    if (M->order==2) {
        double sin_2ksl, cos_2ksl;
        sin_2ksl = sin(2*ksl);
        cos_2ksl = cos(2*ksl);
        
        temp = ksl*sin_2ksl;
        T[0][5][0] = T[1][5][1] = T[2][5][2] = T[3][5][3] = temp;
        
        T[0][5][1] = T[2][5][3] = sin_2ksl/(2*ks) - length*cos_2ksl;
        
        temp = -ksl*cos_2ksl;
        T[0][5][2] = T[1][5][3] = temp;
        T[3][5][1] = T[2][5][0] = -temp;
        
        T[2][5][1] = -(T[0][5][3] = (1.0 - cos_2ksl)/(2*ks) - length*sin_2ksl);
        T[1][5][0] = T[3][5][2] = 0.5*ks*(2*ksl*cos_2ksl + sin_2ksl);
        T[1][5][2] = T[3][5][0] = 0.5*ks*(1.0 - cos_2ksl + 2*ksl*sin_2ksl);
        T[3][5][0] = -T[1][5][2];
        
        T[4][1][1] = T[4][3][3] = length/2;
        }
    log_exit("solenoid_matrix");
    return(M);
    }


VMATRIX *compute_matrix(
                        ELEMENT_LIST *elem,
                        RUN *run,
                        VMATRIX *Mspace        /* currently ignored--intended to allow memory to be passed to put matrix in */
                        )
{
    QUAD *quad; BEND *bend; SEXT *sext; HCOR *hcor; HVCOR *hvcor;
    VCOR *vcor; ALPH *alph; WATCH *watch; DRIFT *drift;
    SOLE *sole; ROTATE *rot; QFRING *qfring;
    MONI *moni; HMON *hmon; VMON *vmon; 
    KSEXT *ksext; KSBEND *ksbend; KQUAD *kquad; NIBEND *nibend; NISEPT *nisept;
    SAMPLE *sample; STRAY *stray; CSBEND *csbend; RFCA *rfca; ENERGY *energy;
    MATTER *matter; MALIGN *malign; MATR *matr; MODRF *modrf;
    long bend_flags;
    
    log_entry("compute_matrix");
    
    if (elem->pred)
        elem->Pref_input = elem->pred->Pref_output;
    else 
        elem->Pref_input = run->p_central;
    elem->Pref_output = elem->Pref_input;

    elem->matrix = NULL;
    
    switch (elem->type) {
      case T_DRIF:
        drift = (DRIFT*)elem->p_elem;
        elem->matrix = drift_matrix(drift->length, (drift->order?drift->order:run->default_order));
        break;
      case T_RBEN: case T_SBEN:
        bend = (BEND*)elem->p_elem;
        bend_flags = determine_bend_flags(elem, bend->edge1_effects, bend->edge2_effects);
        elem->matrix = 
            bend_matrix(bend->length, bend->angle, bend->e1, bend->e2, bend->k1, bend->k2, bend->tilt, bend->fint, 
                        bend->hgap*2, bend->fse, bend->etilt,
                        bend->order?bend->order:run->default_order, bend->edge_order, bend_flags,
                        bend->TRANSPORT);
        if (bend->dx || bend->dy || bend->dz) {
            if (bend->tilt)
                bomb("can't misalign tilted bending magnet--sorry.", NULL);
            misalign_matrix(elem->matrix, bend->dx, bend->dy, bend->dz, bend->angle);
            }
        break;
      case T_QUAD:
        quad = (QUAD*)elem->p_elem;
        elem->matrix = quadrupole_matrix(quad->k1, quad->length, 
                                         quad->order?quad->order:run->default_order, quad->tilt, quad->ffringe,
                                         quad->fse);
        if (quad->dx || quad->dy || quad->dz)
            misalign_matrix(elem->matrix, quad->dx, quad->dy, quad->dz, 0.0);
        break;
      case T_SEXT:
        sext = (SEXT*)elem->p_elem;
        elem->matrix = sextupole_matrix(sext->k2, sext->length, 
                                        sext->order?sext->order:run->default_order, sext->tilt,
                                        sext->fse);
        if (sext->dx || sext->dy || sext->dz)
            misalign_matrix(elem->matrix, sext->dx, sext->dy, sext->dz, 0.0);
        break;
      case T_ALPH:
        alph = (ALPH*)elem->p_elem;
#if DEBUG
        print_elem(stderr, elem);
        fprintf(stderr, "part = %ld, threshold = %e, order = %ld\nxmax = %e, xs1 = %e, xs2 = %e\n",
               alph->part, alph->threshold, alph->order, alph->xs1, alph->xs2);
        fprintf(stderr, "dp1 = %e, dp2 = %e, dx = %e, dy = %e, gradient = %e\n",
               alph->dp1, alph->dp2, alph->dx, alph->dy, alph->gradient);
        print_elem(stderr, elem);
#endif
        if (alph->xmax)
            alph->gradient = run->p_central*sqr(ALPHA_CONST/alph->xmax);
        else
            bomb("supply xmax for alpha magnet", NULL);
        elem->matrix = alpha_magnet_matrix(alph->gradient, sqrt(sqr(run->p_central)+1),
                                           (alph->order?alph->order:run->default_order), alph->part);
        if ((alph->xs1 || alph->xs2 || alph->dp1!=-1 || alph->dp2!=1 ||
             alph->xPuck!=-1 || alph->widthPuck!=0) && alph->part==0)
            bomb("alpha-magnet scraper not supported for full magnet", NULL);
        if (alph->dx || alph->dy || alph->dz)
            misalign_matrix(elem->matrix, alph->dx, alph->dy, alph->dz, 0.0);
        break;  
      case T_ROTATE:
        rot = (ROTATE*)elem->p_elem;
        elem->matrix = rotation_matrix(rot->tilt);
        break;
      case T_MONI:
        moni = (MONI*)elem->p_elem;
        elem->matrix = drift_matrix(moni->length, moni->order?moni->order:run->default_order);
        break;
      case T_HMON:
        hmon = (HMON*)elem->p_elem;
        elem->matrix = drift_matrix(hmon->length, hmon->order?hmon->order:run->default_order);
        break;
      case T_VMON:
        vmon = (VMON*)elem->p_elem;
        elem->matrix = drift_matrix(vmon->length, vmon->order?vmon->order:run->default_order);
        break;
      case T_SOLE:
        sole = (SOLE*) elem->p_elem;
        if (sole->B && !sole->ks)
          sole->ks = -sole->B*e_mks/(me_mks*c_mks*elem->Pref_input);
        elem->matrix = solenoid_matrix(sole->length, sole->ks,
                                       sole->order?sole->order:run->default_order);
        if (sole->dx || sole->dy || sole->dz)
            misalign_matrix(elem->matrix, sole->dx, sole->dy, sole->dz, 0.0);
        break;
      case T_VCOR:
        vcor = (VCOR*) elem->p_elem;
        elem->matrix = corrector_matrix(vcor->length, vcor->kick, vcor->tilt+PI/2,
                                        vcor->b2, vcor->calibration,
                                        vcor->edge_effects, vcor->order?vcor->order:run->default_order);
        break;
      case T_HCOR:
        hcor = (HCOR*) elem->p_elem;
        elem->matrix = corrector_matrix(hcor->length, hcor->kick, hcor->tilt, 
                                        hcor->b2, hcor->calibration,
                                        hcor->edge_effects, hcor->order?hcor->order:run->default_order);
        break;
      case T_MATR:
        matr = (MATR*)elem->p_elem;
        if (!matr->matrix_read) 
            bomb("matrix not read in for MATR element", NULL);
        elem->matrix = &(((MATR*)elem->p_elem)->M);
        break;
      case T_WATCH:
        watch = (WATCH*)elem->p_elem;
        if (!watch->initialized) {
            elem->matrix = NULL;
            set_up_watch_point(watch, run);
            }
        break;
      case T_QFRING:
        qfring = (QFRING*)elem->p_elem;
        elem->matrix = qfringe_matrix(qfring->k1, qfring->length, 
                                      qfring->tilt, qfring->direction, qfring->order?qfring->order:run->default_order,
                                      qfring->fse);
        if (qfring->dx || qfring->dy || qfring->dz)
            misalign_matrix(elem->matrix, qfring->dx, qfring->dy, qfring->dz, 0.0);
        break;
      case T_MALIGN: 
        malign = (MALIGN*)elem->p_elem;
        if (malign->on_pass==-1)
            elem->matrix = misalignment_matrix((MALIGN*)elem->p_elem, run->default_order);
        else
            elem->matrix = drift_matrix(0.0, run->default_order);
        break;
      case T_KSBEND:
        ksbend = (KSBEND*)elem->p_elem;
        if (ksbend->n_kicks<1)
            bomb("n_kicks must be > 0 for KSBEND element", NULL);
        ksbend->flags = determine_bend_flags(elem, ksbend->edge1_effects, ksbend->edge2_effects);
        elem->matrix = 
            bend_matrix(ksbend->length, ksbend->angle, ksbend->e1, ksbend->e2, ksbend->k1, 
                        ksbend->k2, ksbend->tilt, ksbend->fint, 
                        ksbend->hgap*2, ksbend->fse, ksbend->etilt,
                        (run->default_order?run->default_order:1), ksbend->edge_order, ksbend->flags,
                        ksbend->TRANSPORT);
        if (ksbend->dx || ksbend->dy || ksbend->dz) {
            if (ksbend->tilt)
                bomb("can't misalign tilted bending magnet", NULL);
            misalign_matrix(elem->matrix, ksbend->dx, ksbend->dy, ksbend->dz, ksbend->angle);
            }
        break;
      case T_KQUAD:
        kquad = (KQUAD*)elem->p_elem;
        if (kquad->bore)
            kquad->k1 = kquad->B/kquad->bore*(e_mks/(me_mks*c_mks*elem->Pref_input));
        if (kquad->n_kicks<1)
            bomb("n_kicks must by > 0 for KQUAD element", NULL);
        elem->matrix = quadrupole_matrix(kquad->k1, kquad->length, 
                                         (run->default_order?run->default_order:1), kquad->tilt, 0.0,
                                         kquad->fse);
        if (kquad->dx || kquad->dy || kquad->dz)
            misalign_matrix(elem->matrix, kquad->dx, kquad->dy, kquad->dz, 0.0);
        readErrorMultipoleData(&(kquad->systematicMultipoleData),
                                  kquad->systematic_multipoles);
        readErrorMultipoleData(&(kquad->randomMultipoleData),
                                  kquad->random_multipoles);
        break;
      case T_KSEXT:
        ksext = (KSEXT*)elem->p_elem;
        if (ksext->bore)
            ksext->k2 = 2*ksext->B/sqr(ksext->bore)*(e_mks/(me_mks*c_mks*elem->Pref_input));
        if (ksext->n_kicks<1)
            bomb("n_kicks must by > 0 for KSEXT element", NULL);
        elem->matrix = sextupole_matrix(ksext->k2, ksext->length, 
                                        (run->default_order?run->default_order:2), ksext->tilt,
                                        ksext->fse);
        if (ksext->dx || ksext->dy || ksext->dz)
            misalign_matrix(elem->matrix, ksext->dx, ksext->dy, ksext->dz, 0.0);
        readErrorMultipoleData(&(ksext->systematicMultipoleData),
                                  ksext->systematic_multipoles);
        readErrorMultipoleData(&(ksext->randomMultipoleData),
                                  ksext->random_multipoles);
        break;
      case T_MAGNIFY:
        elem->matrix = magnification_matrix((MAGNIFY*)elem->p_elem);
        break;
      case T_SAMPLE:
        sample = (SAMPLE*)elem->p_elem;
        if (sample->interval<=0)
            bomb("sample interval invalid", NULL);
        if (sample->fraction>1 || sample->fraction<=0)
            bomb("sample fraction invalid", NULL);
        break;
      case T_HVCOR:
        hvcor = (HVCOR*) elem->p_elem;
        elem->matrix = hvcorrector_matrix(hvcor->length, hvcor->xkick, hvcor->ykick, 
                                          hvcor->tilt, hvcor->b2, hvcor->xcalibration, hvcor->ycalibration, 
                                          hvcor->edge_effects, hvcor->order?hvcor->order:run->default_order);
        break;
      case T_NIBEND:
        nibend = (NIBEND*)elem->p_elem;
        bend_flags = determine_bend_flags(elem, 1, 1);
        elem->matrix = 
            bend_matrix(nibend->length, nibend->angle, nibend->e1, nibend->e2, 0.0,
                        0.0, nibend->tilt, nibend->fint, nibend->hgap*2,
                        nibend->fse, nibend->etilt,
                        (run->default_order?run->default_order:1), 0L, bend_flags, 0);
        break;
      case T_NISEPT:
        nisept = (NISEPT*)elem->p_elem;
        bend_flags = determine_bend_flags(elem, 1, 1);
        elem->matrix = 
            bend_matrix(nisept->length, nisept->angle, nisept->e1, nisept->angle-nisept->e1, 0.0,
                        0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                        (run->default_order?run->default_order:1), 0, bend_flags, 0);
        break;
      case T_STRAY:
        stray = (STRAY*)elem->p_elem;
        elem->matrix = stray_field_matrix(stray->length, &stray->lBx, &stray->gBx, elem->end_theta,
                                          (stray->order?stray->order:run->default_order),
                                          run->p_central);
        break;
      case T_CSBEND:
        csbend = (CSBEND*)elem->p_elem;
        if (csbend->n_kicks<1)
            bomb("n_kicks must be > 0 for CSBEND element", NULL);
        csbend->flags = determine_bend_flags(elem, csbend->edge1_effects, csbend->edge2_effects);
        elem->matrix = 
            bend_matrix(csbend->length, csbend->angle, csbend->e1, csbend->e2, csbend->k1, 
                        csbend->k2, csbend->tilt, csbend->fint, 
                        csbend->hgap*2, csbend->fse, csbend->etilt,
                        (run->default_order?run->default_order:1), 1L, csbend->flags, 0);
        if (csbend->dx || csbend->dy || csbend->dz) {
            if (csbend->tilt)
                bomb("can't misalign tilted bending magnet", NULL);
            misalign_matrix(elem->matrix, csbend->dx, csbend->dy, csbend->dz, csbend->angle);
            }
        break;
      case T_RFCA: 
        rfca = (RFCA*)elem->p_elem;
        elem->matrix = rf_cavity_matrix(rfca->length, rfca->volt, rfca->freq, rfca->phase, 
                                        &elem->Pref_output, run->default_order?run->default_order:1);
        if (!rfca->change_p0) {
            elem->matrix->C[5] = (elem->Pref_output-elem->Pref_input)/elem->Pref_input;
            elem->Pref_output = elem->Pref_input;
            }
        break;
      case T_MODRF: 
        modrf = (MODRF*)elem->p_elem;
        elem->matrix = rf_cavity_matrix(modrf->length, modrf->volt, modrf->freq, modrf->phase, 
                                        &elem->Pref_output, run->default_order?run->default_order:1);
        elem->matrix->C[5] = (elem->Pref_output-elem->Pref_input)/elem->Pref_input;
        elem->Pref_output = elem->Pref_input;
        break;
      case T_ENERGY:
        energy = (ENERGY*)elem->p_elem;
        if (energy->central_energy || energy->central_momentum) {
            if (energy->central_energy) 
                elem->Pref_output = sqrt(sqr(energy->central_energy)-1);
            else
                elem->Pref_output = energy->central_momentum;
            }
        break;
      case T_MATTER:
        matter = (MATTER*)elem->p_elem;
        elem->matrix = drift_matrix(matter->length, run->default_order);
        break;
      case T_KPOLY: case T_RFDF:  case T_RFTM:  case T_RMDF:  case T_TMCF: case T_CEPL:  case T_TWPL:  case T_TWLA:  
      case T_TWMTA: case T_RCOL:  case T_PEPPOT: case T_MAXAMP: case T_ECOL: case T_TRCOUNT: 
      case T_RECIRC: case T_SCRAPER: case T_CENTER: case T_MULT: case T_SCATTER: case T_RAMPRF: case T_RAMPP: 
      case T_KICKER: case T_RFMODE:
      default:
        if (entity_description[elem->type].flags&HAS_LENGTH)
            elem->matrix = drift_matrix(*((double*)elem->p_elem), run->default_order);
        if ((entity_description[elem->type].flags&HAS_MATRIX) && !elem->matrix) {
            fprintf(stderr, "error: failed to compute matrix for %s, which should have a matrix.\n",
                    elem->name);
            abort();
            }
        break;
        }
    
    log_exit("compute_matrix");
    return(elem->matrix);
    }

char *watch_mode[N_WATCH_MODES] = {
    "coordinates", "parameters", "centroids", "fft",
    } ;
char *fft_window_name[N_FFT_WINDOWS] = {
    "hanning", "parzen", "welch", "uniform",
    } ;

void set_up_watch_point(WATCH *watch, RUN *run)
{
    char *mode, *qualifier;

    log_entry("set_up_watch_point");
    
    if (watch->interval<=0 || watch->fraction<=0)
        bomb("interval or fraction is non-positive for WATCH element", NULL);
    if (!watch->mode || watch->mode[0]==0)
        bomb("mode must be given for WATCH element", NULL);
    mode = watch->mode;
    if (qualifier=strchr(mode, ' '))
         *qualifier++ = 0;
    if ((watch->mode_code=match_string(watch->mode, watch_mode, N_WATCH_MODES, 0))<0)
        bomb("unknown watch mode", NULL);
    if (watch->mode_code!=WATCH_COORDINATES)
      watch->start_pass = 0;
    if (watch->label && str_in(watch->label, "%s")) {
        char *buffer;
        buffer = tmalloc(sizeof(*buffer)*(strlen(watch->label)+strlen(run->rootname)+1));
        sprintf(buffer, watch->label, run->rootname);
        free(watch->label);
        watch->label = buffer;
        }
    watch->filename = compose_filename(watch->filename, run->rootname);
    SDDS_WatchPointSetup(watch, SDDS_BINARY, 1, run->runfile, run->lattice, "set_up_watch_point", qualifier);
    watch->initialized = 1;
    watch->count = 0;
    log_exit("set_up_watch_point");
    }

VMATRIX *magnification_matrix(MAGNIFY *magnif)
{
    VMATRIX *M;
    
    log_entry("magnification_matrix");
    
    M = tmalloc(sizeof(*M));
    initialize_matrices(M, M->order=1);
    M->R[0][0] = magnif->mx;
    M->R[1][1] = magnif->mxp;
    M->R[2][2] = magnif->my;
    M->R[3][3] = magnif->myp;
    M->R[4][4] = magnif->ms;
    M->R[5][5] = magnif->mdp;
    
    log_exit("magnification_matrix");
    return(M);
    }


void reset_special_elements(LINE_LIST *beamline)
{
    ELEMENT_LIST *eptr;
    NIBEND *nibend; NISEPT *nisept;
    
    log_entry("reset_special_elements");
    
    eptr = &(beamline->elem);
    while (eptr) {
        switch (eptr->type) {
          case T_NIBEND:
            nibend = (NIBEND*)eptr->p_elem;
            nibend->zeta_offset = 0;
            break;
          case T_NISEPT:
            nisept = (NISEPT*)eptr->p_elem;
            nisept->fse_opt = 0;
            break;
          case T_TMCF:
            ((TMCF_MODE*)eptr->p_elem)->fiducial_part = NULL;
            break;
          case T_CEPL:
            ((CE_PLATES*)eptr->p_elem)->fiducial_part = NULL;
            break;
          case T_TWPL:
            ((TW_PLATES*)eptr->p_elem)->fiducial_part = NULL;
            break;
          case T_TWLA:
            ((TW_LINAC*)eptr->p_elem)->fiducial_part = NULL;
            break;
          case T_TWMTA:
            ((TWMTA*)eptr->p_elem)->fiducial_part = NULL;
            break;
          case T_RFCA:
            ((RFCA*)eptr->p_elem)->fiducial_seen = 0;
            break;
          case T_MODRF:
            ((MODRF*)eptr->p_elem)->fiducial_seen = 0;
            break;
          case T_RFMODE:
            ((RFMODE*)eptr->p_elem)->initialized = 0;
            if (((RFMODE*)eptr->p_elem)->fprec)
                fclose(((RFMODE*)eptr->p_elem)->fprec);
            break;
          case T_TRFMODE:
            ((TRFMODE*)eptr->p_elem)->initialized = 0;
            if (((TRFMODE*)eptr->p_elem)->fprec)
                fclose(((TRFMODE*)eptr->p_elem)->fprec);
            break;
          case T_RAMPRF:
            ((RAMPRF*)eptr->p_elem)->Ts = 0;
            ((RAMPRF*)eptr->p_elem)->fiducial_seen = 0;
            break;
          case T_KQUAD:
            ((KQUAD*)eptr->p_elem)->randomMultipoleData.randomized = 0;
            break;
          case T_KSEXT:
            ((KSEXT*)eptr->p_elem)->randomMultipoleData.randomized = 0;
            break;
          default:
            break;
            }
        eptr = eptr->succ;
        }
    log_exit("reset_special_elements");
    }

VMATRIX *stray_field_matrix(double length, double *lB, double *gB, double theta, long order, double p_central)
{
    /* factor in equation theta = CTHETA*B(T)*L(m)/(beta*gamma): */
#define CTHETA 5.8667907921396181e+02
    VMATRIX *M;
    double xkick, ykick, Bx, By;
    double rho;
    
    log_entry("stray_field_matrix");
    
    Bx = lB[0] + cos(theta)*gB[0] + sin(theta)*gB[3];
    By = lB[1] + gB[1];
    if (By) {
        rho = p_central/(CTHETA*By);
        xkick = -asin(length/rho);
        }
    else
        xkick = 0;
    if (Bx) {
        rho = p_central/(CTHETA*Bx);
        ykick = asin(length/rho);
        }
    else
        ykick = 0;
    
    M = hvcorrector_matrix(length, xkick, ykick, 0.0, 0.0, 1.0, 1.0, 0, order);
    
    log_exit("stray_field_matrix");
    return(M);
    }


VMATRIX *rf_cavity_matrix(double length, double voltage, double frequency, double phase, double *P_central, long order)
{
    VMATRIX *M;
    double *C, **R, dP, gamma, dgamma;
    double cos_phase, sin_phase;

    log_entry("rf_cavity_matrix");

    M = tmalloc(sizeof(*M));
    M->order = 1;
    initialize_matrices(M, M->order);
    R = M->R;
    C = M->C;
    
    if (*P_central<=0) {
        fprintf(stderr, "error: P_central = %g\n", *P_central);
        abort();
        }

    C[4] = length;
    voltage /= 1e6;  /* convert to MV */
    sin_phase = sin(PI*phase/180.0);
    cos_phase = cos(PI*phase/180.0);
    dgamma = voltage/me_mev*sin_phase;
    gamma = sqrt(sqr(*P_central)+1);
    dP    = sqrt(sqr(gamma+dgamma)-1) - *P_central;

    R[0][0] = R[2][2] = R[4][4] = 1;
    R[1][1] = R[3][3] = R[5][5] = 1/(1+dP/(*P_central));
    if (fabs(dP/(*P_central))>1e-14)
        R[0][1] = R[2][3] = length*(*P_central)/dP*log(1 + dP/(*P_central));
    else
        R[0][1] = R[2][3] = length;
    R[5][4] = (voltage/me_mev)*cos_phase/(gamma + dgamma*sin_phase)*(PIx2*frequency/c_mks);

    *P_central += dP;

    log_exit("rf_cavity_matrix");
    return(M);
    }

