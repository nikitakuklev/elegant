/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* contents: compute_matrices(), drift_matrix(), sextupole_matrix() 
 *           plus more...
 * Michael Borland, 1989.
 */
#include "track.h"
#include "mdb.h"
#include "matlib.h"

#define DEBUG 0

long determine_bend_flags(ELEMENT_LIST *eptr, long edge1_effects, long edge2_effects);
VMATRIX *matrixFromExplicitMatrix(EMATRIX *emat, long order);

VMATRIX *full_matrix(ELEMENT_LIST *elem, RUN *run, long order) 
{

    if (!elem) {
        fputs("error: NULL element pointer passed to full_matrix", stdout);
        abort();
        }
    
#ifdef WATCH_MEMORY
    fprintf(stdout, "start full_matrix: CPU: %6.2lf  PF: %6ld  MEM: %6ld\n",
           cpu_time()/100.0, page_faults(), memory_count());
    fflush(stdout);
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
    fputs("error: NULL element pointer passed to accumulate_matrices", stdout);
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
      fprintf(stdout, "error: bad element type %ld (accumulate_matrices)\n", member->type);
      fflush(stdout);
      fprintf(stdout, "element name is %s and end position is %em\n", 
             (member->name?member->name:"{null}"), member->end_pos);
      fflush(stdout);
      abort();
    }
    if (member->pred)
      Pref_input = member->pred->Pref_output;
    else
      Pref_input = member->Pref_input;
    if (!member->matrix || Pref_input!=member->Pref_input)
      compute_matrix(member, run, NULL);
    if ((entity_description[member->type].flags&HAS_MATRIX) && !member->matrix) {
      fprintf(stdout, "programming error: matrix not computed for element %s\n",
              member->name);
      fflush(stdout);
      abort();
    }
    if (member->matrix) {
      concat_matrices(M2, member->matrix, M1,
                      entity_description[member->type].flags&HAS_RF_MATRIX?
                      CONCAT_EXCLUDE_S0:0);
      tmp = M2;
      M2  = M1;
      M1  = tmp;
      if (!full_matrix_only) {
        if (member->accumMatrix)
          free_matrices(member->accumMatrix);
        else 
          member->accumMatrix = tmalloc(sizeof(*(member->accumMatrix)));
        copy_matrices(member->accumMatrix, M1);
      }
    } else if (!full_matrix_only) {
      if (member->accumMatrix)
        free_matrices(member->accumMatrix);
      else 
        member->accumMatrix = tmalloc(sizeof(*(member->accumMatrix)));
      copy_matrices(member->accumMatrix, M1);
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
        fputs("error: NULL element pointer passed to append_full_matrix", stdout);
        abort();
        }
    if (!M0) {
        fputs("error: NULL initial matrix pointer passed to append_full_matrix", stdout);
        abort();
        }
    
    initialize_matrices(M1=tmalloc(sizeof(*M1)), order);
    initialize_matrices(M2=tmalloc(sizeof(*M2)), order);
    copy_matrices1(M1, M0);
    
    member = elem;
    while (member) {
        if (member->type<0 || member->type>=N_TYPES) {
            fprintf(stdout, "error: bad element type %ld (full_matrix)\n", member->type);
            fflush(stdout);
            fprintf(stdout, "element name is %s and end position is %em\n", 
                   (member->name?member->name:"{null}"), member->end_pos);
            fflush(stdout);
            abort();
            }
        if (member->pred)
            Pref_input = member->pred->Pref_output;
        else
            Pref_input = member->Pref_input;
        if (!member->matrix || Pref_input!=member->Pref_input)
            compute_matrix(member, run, NULL);
        if ((entity_description[member->type].flags&HAS_MATRIX) && !member->matrix) {
            fprintf(stdout, "programming error: matrix not computed for element %s\n",
                    member->name);
            fflush(stdout);
            abort();
            }
        if (member->matrix) {
          concat_matrices(M2, member->matrix, M1,
                          entity_description[member->type].flags&HAS_RF_MATRIX?
                        CONCAT_EXCLUDE_S0:0);
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
            fprintf(stdout, "error: bad element type %ld (fill_in_matrices)\n", member->type);
            fflush(stdout);
            fprintf(stdout, "element name is %s and end position is %em\n", 
                   (member->name?member->name:"{null}"), member->end_pos);
            fflush(stdout);
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

VMATRIX *wiggler_matrix(double length, double radius, long order)
{
    VMATRIX *M;
    double **R, *C;
    double kl;

    M = tmalloc(sizeof(*M));
    M->order = 1;
    initialize_matrices(M, M->order);
    R = M->R;
    C = M->C;

    R[0][0] = R[1][1] = R[2][2] = R[3][3] = R[4][4] = R[5][5] = 1;

    if (length) {
      C[4] = length;
      R[0][1] = length;
      kl = length/(SQRT2*fabs(radius));
      R[2][2] = R[3][3] = cos(kl);
      R[2][3] = sin(kl)/(kl/length);
      R[3][2] = -(kl/length)*sin(kl);
    }

    return(M);
    }

VMATRIX *sextupole_matrix(double K2, double length, long maximum_order, double tilt, double fse)
{
    VMATRIX *M;
    double *C, **R, ***T, ****U;
    double temp;
    
    log_entry("sextupole_matrix");
    
    K2 *= (1+fse);

    M = tmalloc(sizeof(*M));
    initialize_matrices(M, M->order=MIN(3,maximum_order));
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

      if (M->order >= 3) {
        U = M->Q;
        U[0][0][0][0] = ipow(K2,2)*ipow(length,4)/48.0 ;
        U[0][1][0][0] = ipow(K2,2)*ipow(length,5)/48.0 ;
        U[0][1][1][0] = ipow(K2,2)*ipow(length,6)/144.0 ;
        U[0][1][1][1] = ipow(K2,2)*ipow(length,7)/1008.0 ;
        U[0][2][2][0] = ipow(K2,2)*ipow(length,4)/48.0 ;
        U[0][2][2][1] = -ipow(K2,2)*ipow(length,5)/240.0 ;
        U[0][3][2][0] = ipow(K2,2)*ipow(length,5)/40.0 ;
        U[0][3][2][1] = ipow(K2,2)*ipow(length,6)/360.0 ;
        U[0][3][3][0] = ipow(K2,2)*ipow(length,6)/240.0 ;
        U[0][3][3][1] = ipow(K2,2)*ipow(length,7)/1008.0 ;
        U[0][5][0][0] = K2*ipow(length,2)/4.0 ;
        U[0][5][1][0] = K2*ipow(length,3)/6.0 ;
        U[0][5][1][1] = K2*ipow(length,4)/24.0 ;
        U[0][5][2][2] = -K2*ipow(length,2)/4.0 ;
        U[0][5][3][2] = -K2*ipow(length,3)/6.0 ;
        U[0][5][3][3] = -K2*ipow(length,4)/24.0 ;
        U[1][0][0][0] = ipow(K2,2)*ipow(length,3)/12.0 ;
        U[1][1][0][0] = 5.0*ipow(K2,2)*ipow(length,4)/48.0 ;
        U[1][1][1][0] = ipow(K2,2)*ipow(length,5)/24.0 ;
        U[1][1][1][1] = ipow(K2,2)*ipow(length,6)/144.0 ;
        U[1][2][2][0] = ipow(K2,2)*ipow(length,3)/12.0 ;
        U[1][2][2][1] = -ipow(K2,2)*ipow(length,4)/48.0 ;
        U[1][3][2][0] = ipow(K2,2)*ipow(length,4)/8.0 ;
        U[1][3][2][1] = ipow(K2,2)*ipow(length,5)/60.0 ;
        U[1][3][3][0] = ipow(K2,2)*ipow(length,5)/40.0 ;
        U[1][3][3][1] = ipow(K2,2)*ipow(length,6)/144.0 ;
        U[1][5][0][0] = K2*length/2.0 ;
        U[1][5][1][0] = K2*ipow(length,2)/2.0 ;
        U[1][5][1][1] = K2*ipow(length,3)/6.0 ;
        U[1][5][2][2] = -K2*length/2.0 ;
        U[1][5][3][2] = -K2*ipow(length,2)/2.0 ;
        U[1][5][3][3] = -K2*ipow(length,3)/6.0 ;
        U[2][2][0][0] = ipow(K2,2)*ipow(length,4)/48.0 ;
        U[2][2][1][0] = ipow(K2,2)*ipow(length,5)/40.0 ;
        U[2][2][1][1] = ipow(K2,2)*ipow(length,6)/240.0 ;
        U[2][2][2][2] = ipow(K2,2)*ipow(length,4)/48.0 ;
        U[2][3][0][0] = -ipow(K2,2)*ipow(length,5)/240.0 ;
        U[2][3][1][0] = ipow(K2,2)*ipow(length,6)/360.0 ;
        U[2][3][1][1] = ipow(K2,2)*ipow(length,7)/1008.0 ;
        U[2][3][2][2] = ipow(K2,2)*ipow(length,5)/48.0 ;
        U[2][3][3][2] = ipow(K2,2)*ipow(length,6)/144.0 ;
        U[2][3][3][3] = ipow(K2,2)*ipow(length,7)/1008.0 ;
        U[2][5][2][0] = -K2*ipow(length,2)/2.0 ;
        U[2][5][2][1] = -K2*ipow(length,3)/6.0 ;
        U[2][5][3][0] = -K2*ipow(length,3)/6.0 ;
        U[2][5][3][1] = -K2*ipow(length,4)/12.0 ;
        U[3][2][0][0] = ipow(K2,2)*ipow(length,3)/12.0 ;
        U[3][2][1][0] = ipow(K2,2)*ipow(length,4)/8.0 ;
        U[3][2][1][1] = ipow(K2,2)*ipow(length,5)/40.0 ;
        U[3][2][2][2] = ipow(K2,2)*ipow(length,3)/12.0 ;
        U[3][3][0][0] = -ipow(K2,2)*ipow(length,4)/48.0 ;
        U[3][3][1][0] = ipow(K2,2)*ipow(length,5)/60.0 ;
        U[3][3][1][1] = ipow(K2,2)*ipow(length,6)/144.0 ;
        U[3][3][2][2] = 5.0*ipow(K2,2)*ipow(length,4)/48.0 ;
        U[3][3][3][2] = ipow(K2,2)*ipow(length,5)/24.0 ;
        U[3][3][3][3] = ipow(K2,2)*ipow(length,6)/144.0 ;
        U[3][5][2][0] = -K2*length ;
        U[3][5][2][1] = -K2*ipow(length,2)/2.0 ;
        U[3][5][3][0] = -K2*ipow(length,2)/2.0 ;
        U[3][5][3][1] = -K2*ipow(length,3)/3.0 ;
        U[4][1][0][0] = -K2*ipow(length,2)/4.0 ;
        U[4][1][1][0] = -K2*ipow(length,3)/6.0 ;
        U[4][1][1][1] = -K2*ipow(length,4)/24.0 ;
        U[4][2][2][1] = K2*ipow(length,2)/4.0 ;
        U[4][3][2][0] = K2*ipow(length,2)/2.0 ;
        U[4][3][2][1] = K2*ipow(length,3)/3.0 ;
        U[4][3][3][0] = K2*ipow(length,3)/6.0 ;
        U[4][3][3][1] = K2*ipow(length,4)/8.0 ;
      }
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
        /* These terms per P. Emma */
        T[4][0][0] = T[4][2][2] = sqr(ks)*length/2;
        T[4][3][0] = -(T[4][2][1] = ks*length);
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
    VCOR *vcor; ALPH *alph; DRIFT *drift;
    SOLE *sole; ROTATE *rot; QFRING *qfring;
    MONI *moni; HMON *hmon; VMON *vmon; 
    KSEXT *ksext; KSBEND *ksbend; KQUAD *kquad; NIBEND *nibend; NISEPT *nisept;
    SAMPLE *sample; STRAY *stray; CSBEND *csbend; RFCA *rfca; ENERGY *energy;
    RFCW *rfcw; 
    MATTER *matter; MALIGN *malign; MATR *matr; MODRF *modrf;
    CSRCSBEND *csrcsbend;
    CSRDRIFT *csrdrift;
    LSCDRIFT *lscdrift;
    WIGGLER *wiggler;
    double ks, wigglerRadius, Pref_output;
    VARY rcContext;
    long fiducialize;

    getRunControlContext(&rcContext);
    fiducialize = 1;
    if (rcContext.ready) {
      if ((rcContext.fiducial_flag&FIRST_BEAM_IS_FIDUCIAL) &&
	  rcContext.i_step!=0)
	fiducialize = 0;
    }
    if (elem->pred)
        elem->Pref_input = elem->pred->Pref_output;
    else 
        elem->Pref_input = run->p_central;

    /* Pref_output is the assumed output value of Pref */
    Pref_output = elem->Pref_input;
    if (!fiducialize)
      /* Assume we've already fiducialized and use the previous Pref output value */
      Pref_output = elem->Pref_output;
    /* This variable is used to pass the input momentum to some elements and 
     * get the output momentum back.  It will be reset to the local variable
     * Pref_output if fiducializing 
     */
    elem->Pref_output = elem->Pref_input;

    elem->matrix = NULL;
    
    switch (elem->type) {
      case T_DRIF:
        drift = (DRIFT*)elem->p_elem;
        elem->matrix = drift_matrix(drift->length, (drift->order?drift->order:run->default_order));
        break;
      case T_WIGGLER:
	wiggler = (WIGGLER*)elem->p_elem;
	if (wiggler->K) {
	  long poles;
	  double period;
	  poles = 2*(wiggler->poles/2)+1;
	  period = 2*(wiggler->length/poles);
	  wiggler->radiusInternal = sqrt(sqr(elem->Pref_input)+1)*period/(PIx2*wiggler->K);
	} else
	  wiggler->radiusInternal = wiggler->radius;
	if (wiggler->radiusInternal==0) {
	  fprintf(stderr, "Error: wiggler radius is zero\n");
	  fprintf(stderr, "Parameters are length=%e, poles=%ld, radius=%e, K=%e\n",
		  wiggler->length, wiggler->poles, wiggler->radiusInternal, wiggler->K);
	}
        elem->matrix = wiggler_matrix(wiggler->length, wiggler->radiusInternal,
				      run->default_order);
	break;
      case T_SCRIPT:
        elem->matrix = drift_matrix(((SCRIPT*)elem->p_elem)->length, 
				    run->default_order);
	break;
      case T_RBEN: case T_SBEN:
        bend = (BEND*)elem->p_elem;
        bend->edgeFlags = determine_bend_flags(elem, bend->edge1_effects, bend->edge2_effects);
        if (bend->use_bn) {
          bend->k1_internal = bend->b1/(bend->length/bend->angle);
          bend->k2_internal = bend->b2/(bend->length/bend->angle);
        }
        else {
          bend->k1_internal = bend->k1;
          bend->k2_internal = bend->k2;
        }
        elem->matrix = 
            bend_matrix(bend->length, bend->angle, bend->e1, bend->e2, 
                        bend->h1, bend->h2, 
                        bend->k1_internal, bend->k2_internal, 
                        bend->tilt, bend->fint, 
                        bend->hgap*2, bend->fse, bend->etilt,
                        bend->order?bend->order:run->default_order, bend->edge_order, 
			bend->edgeFlags, bend->TRANSPORT);
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
                                         quad->fse, quad->xkick, quad->ykick, quad->fringeType); 
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
        print_elem(stdout, elem);
        fprintf(stdout, "part = %ld, threshold = %e, order = %ld\nxmax = %e, xs1 = %e, xs2 = %e\n",
               alph->part, alph->threshold, alph->order, alph->xs1, alph->xs2);
        fflush(stdout);
        fprintf(stdout, "dp1 = %e, dp2 = %e, dx = %e, dy = %e, gradient = %e\n",
               alph->dp1, alph->dp2, alph->dx, alph->dy, alph->gradient);
        fflush(stdout);
        print_elem(stdout, elem);
#endif
        if (alph->xmax)
            alph->gradient = elem->Pref_input*sqr(ALPHA_CONST/alph->xmax);
        else
            bomb("supply xmax for alpha magnet", NULL);
        elem->matrix = alpha_magnet_matrix(alph->gradient, sqrt(sqr(run->p_central)+1),
                                           (alph->order?alph->order:run->default_order), alph->part);
        if ((alph->xs1 || alph->xs2 || alph->dp1!=-1 || alph->dp2!=1 ||
             alph->xPuck!=-1 || alph->widthPuck!=0) && alph->part==0)
            bomb("alpha-magnet scraper not supported for full magnet", NULL);
        if (alph->tilt)
            tilt_matrices(elem->matrix, alph->tilt);
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
        ks = sole->ks;
        if (sole->B && !ks)
          ks = -sole->B*e_mks/(me_mks*c_mks*elem->Pref_input);
        elem->matrix = solenoid_matrix(sole->length, ks,
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
        break;
      case T_HISTOGRAM:
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
        if (malign->on_pass==-1 || malign->forceModifyMatrix)
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
            bend_matrix(ksbend->length, ksbend->angle, ksbend->e1, ksbend->e2, 
                        ksbend->h1, ksbend->h2, ksbend->k1, 
                        ksbend->k2, ksbend->tilt, ksbend->fint, 
                        ksbend->hgap*2, ksbend->fse, ksbend->etilt,
                        ksbend->nonlinear?2:(run->default_order?run->default_order:1),
                        ksbend->edge_order, ksbend->flags,
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
                                         kquad->fse, kquad->xkick, kquad->ykick, NULL);
        if (kquad->dx || kquad->dy || kquad->dz)
            misalign_matrix(elem->matrix, kquad->dx, kquad->dy, kquad->dz, 0.0);
        readErrorMultipoleData(&(kquad->systematicMultipoleData),
                                  kquad->systematic_multipoles, 0);
        readErrorMultipoleData(&(kquad->randomMultipoleData),
                                  kquad->random_multipoles, 0);
        readErrorMultipoleData(&(kquad->steeringMultipoleData),
                                  kquad->steering_multipoles, 1);
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
                                  ksext->systematic_multipoles, 0);
        readErrorMultipoleData(&(ksext->randomMultipoleData),
                                  ksext->random_multipoles, 0);
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
        nibend->edgeFlags = determine_bend_flags(elem, 1, 1);
        elem->matrix = 
            bend_matrix(nibend->length, nibend->angle, nibend->e1, nibend->e2, 
                        0.0, 0.0, 0.0, 0.0, nibend->tilt, nibend->fint, nibend->hgap*2,
                        nibend->fse, nibend->etilt,
                        (run->default_order?run->default_order:1), 0L, nibend->edgeFlags, 0);
        break;
      case T_NISEPT:
        nisept = (NISEPT*)elem->p_elem;
        nisept->edgeFlags = determine_bend_flags(elem, 1, 1);
        elem->matrix = 
            bend_matrix(nisept->length, nisept->angle, nisept->e1, nisept->angle-nisept->e1, 
                        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                        (run->default_order?run->default_order:1), 0, nisept->edgeFlags, 0);
        break;
      case T_STRAY:
        stray = (STRAY*)elem->p_elem;
        elem->matrix = 
          stray_field_matrix(stray->length, &stray->lBx, &stray->gBx, elem->end_theta,
                                          (stray->order?stray->order:run->default_order),
                                          elem->Pref_input, stray->WiInitialized?stray->Wi:NULL);
        break;
      case T_CSBEND:
        csbend = (CSBEND*)elem->p_elem;
        if (csbend->n_kicks<1)
            bomb("n_kicks must be > 0 for CSBEND element", NULL);
        csbend->edgeFlags = determine_bend_flags(elem, csbend->edge1_effects, csbend->edge2_effects);
        if (csbend->use_bn) {
          csbend->k1_internal = csbend->b1/(csbend->length/csbend->angle);
          csbend->k2_internal = csbend->b2/(csbend->length/csbend->angle);
          csbend->k3_internal = csbend->b3/(csbend->length/csbend->angle);
          csbend->k4_internal = csbend->b4/(csbend->length/csbend->angle);
        } else {
          csbend->k1_internal = csbend->k1;
          csbend->k2_internal = csbend->k2;
          csbend->k3_internal = csbend->k3;
          csbend->k4_internal = csbend->k4;
        }
        elem->matrix = 
            bend_matrix(csbend->length, csbend->angle, csbend->e1, csbend->e2, 
                        csbend->h1, csbend->h2, csbend->k1_internal, 
                        csbend->k2_internal, csbend->tilt, csbend->fint, 
                        csbend->hgap*2, csbend->fse, csbend->etilt,
                        csbend->nonlinear?2:(run->default_order?run->default_order:1),
                        csbend->edge_order, csbend->edgeFlags, 0);
        if (csbend->dx || csbend->dy || csbend->dz) {
          if (csbend->tilt)
              bomb("can't misalign tilted bending magnet", NULL);
          misalign_matrix(elem->matrix, csbend->dx, csbend->dy, csbend->dz, csbend->angle);
        }
        break;
      case T_CSRCSBEND:
        csrcsbend = (CSRCSBEND*)elem->p_elem;
        if (csrcsbend->n_kicks<1)
            bomb("n_kicks must be > 0 for CSRCSBEND element", NULL);
        csrcsbend->edgeFlags = determine_bend_flags(elem, csrcsbend->edge1_effects, csrcsbend->edge2_effects);
        if (csrcsbend->use_bn) {
          csrcsbend->k1_internal = csrcsbend->b1/(csrcsbend->length/csrcsbend->angle);
          csrcsbend->k2_internal = csrcsbend->b2/(csrcsbend->length/csrcsbend->angle);
          csrcsbend->k3_internal = csrcsbend->b3/(csrcsbend->length/csrcsbend->angle);
          csrcsbend->k4_internal = csrcsbend->b4/(csrcsbend->length/csrcsbend->angle);
        } else {
          csrcsbend->k1_internal = csrcsbend->k1;
          csrcsbend->k2_internal = csrcsbend->k2;
          csrcsbend->k3_internal = csrcsbend->k3;
          csrcsbend->k4_internal = csrcsbend->k4;
        }
        elem->matrix = 
            bend_matrix(csrcsbend->length, csrcsbend->angle, csrcsbend->e1, csrcsbend->e2, 
                        csrcsbend->h1, csrcsbend->h2, csrcsbend->k1_internal, 
                        csrcsbend->k2_internal, csrcsbend->tilt, csrcsbend->fint, 
                        csrcsbend->hgap*2, csrcsbend->fse, csrcsbend->etilt,
                        csrcsbend->nonlinear?2:(run->default_order?run->default_order:1),
                        csrcsbend->edge_order, csrcsbend->edgeFlags, 0);
        if (csrcsbend->dx || csrcsbend->dy || csrcsbend->dz) {
            if (csrcsbend->tilt)
                bomb("can't misalign tilted bending magnet", NULL);
            misalign_matrix(elem->matrix, csrcsbend->dx, csrcsbend->dy, csrcsbend->dz, csrcsbend->angle);
            }
        break;
      case T_RFCA: 
        rfca = (RFCA*)elem->p_elem;
        elem->matrix = rf_cavity_matrix(rfca->length, rfca->volt, rfca->freq, rfca->phase, 
                                        &elem->Pref_output, run->default_order?run->default_order:1,
                                        rfca->end1Focus, rfca->end2Focus,
                                        rfca->bodyFocusModel, 
                                        fiducialize*(rfca->change_p0 || run->always_change_p0),
					Pref_output);
        if (rfca->dx || rfca->dy)
          misalign_matrix(elem->matrix, rfca->dx, rfca->dy, 0.0, 0.0);
        break;
      case T_RFCW: 
        rfcw = (RFCW*)elem->p_elem;
        elem->matrix = rf_cavity_matrix(rfcw->length, rfcw->volt, rfcw->freq, rfcw->phase, 
                                        &elem->Pref_output, run->default_order?run->default_order:1,
                                        rfcw->end1Focus, rfcw->end2Focus,
                                        rfcw->bodyFocusModel, 
                                        fiducialize*(rfca->change_p0 || run->always_change_p0),
					Pref_output);
        if (rfcw->dx || rfcw->dy)
          misalign_matrix(elem->matrix, rfcw->dx, rfcw->dy, 0.0, 0.0);
        break;
      case T_MODRF: 
        modrf = (MODRF*)elem->p_elem;
        elem->matrix = rf_cavity_matrix(modrf->length, modrf->volt, modrf->freq, modrf->phase, 
                                        &elem->Pref_output, run->default_order?run->default_order:1,
                                        0, 0, NULL,
                                        fiducialize*run->always_change_p0, Pref_output);
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
      case T_CSRDRIFT:
        csrdrift = (CSRDRIFT*)elem->p_elem;
        if ((csrdrift->useOvertakingLength?1:0)+
            (csrdrift->spread?1:0)+(csrdrift->attenuationLength>0?1:0)>1) 
          bomb("Give one and only one of SPREAD, ATTENUATION_LENGTH, or USE_OVERTAKING_LENGTH for CSRDRIFT", NULL);
        elem->matrix = drift_matrix(csrdrift->length, run->default_order);
        break;
      case T_LSCDRIFT:
        lscdrift = (LSCDRIFT*)elem->p_elem;
        elem->matrix = drift_matrix(lscdrift->length, run->default_order);
        break;
      case T_TWISSELEMENT:
        elem->matrix = twissTransformMatrix((TWISSELEMENT*)elem->p_elem, NULL);
        break;
      case T_REFLECT:
        elem->matrix = tmalloc(sizeof(*(elem->matrix)));
        initialize_matrices(elem->matrix, elem->matrix->order = 1);
        elem->matrix->R[0][0] = elem->matrix->R[2][2] = 
          elem->matrix->R[4][4] = elem->matrix->R[5][5] = 1;
        elem->matrix->R[1][1] = elem->matrix->R[3][3] = -1;
        break;
      case T_LTHINLENS:
        elem->matrix = lightThinLensMatrix((LTHINLENS*)elem->p_elem);
        break;
      case T_LMIRROR:
        elem->matrix = lightMirrorMatrix((LMIRROR*)elem->p_elem);
        break;
      case T_EMATRIX:
        elem->matrix = matrixFromExplicitMatrix((EMATRIX*)elem->p_elem, run->default_order);
        break;
      case T_KPOLY: case T_RFDF:  case T_RFTMEZ0:  case T_RMDF:  case T_TMCF: case T_CEPL:  
      case T_TWPL:  case T_TWLA:  
      case T_TWMTA: case T_RCOL:  case T_PEPPOT: case T_MAXAMP: 
      case T_ECOL: case T_TRCOUNT: 
      case T_RECIRC: case T_SCRAPER: case T_CENTER: case T_MULT: 
      case T_SCATTER: case T_RAMPRF: case T_RAMPP: 
      case T_KICKER: case T_RFMODE: case T_REMCOR: case T_MAPSOLENOID:
      case T_DSCATTER: case T_LSRMDLTR: 
      default:
        if (entity_description[elem->type].flags&HAS_LENGTH)
            elem->matrix = drift_matrix(*((double*)elem->p_elem), run->default_order);
        if ((entity_description[elem->type].flags&HAS_MATRIX) && !elem->matrix) {
            fprintf(stdout, "error: failed to compute matrix for %s, which should have a matrix.\n",
                    elem->name);
            fflush(stdout);
            abort();
            }
        break;
        }

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
    if ((qualifier=strchr(mode, ' ')))
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
    watch->flushSample = -1;
    log_exit("set_up_watch_point");
    }

void set_up_histogram(HISTOGRAM *histogram, RUN *run)
{
  if (histogram->interval<=0)
    bomb("interval is non-positive for HISTOGRAM element", NULL);
  if (histogram->bins<=2)
    bomb("number of bins is less than 2 for HISTOGRAM element", NULL);
  if (!histogram->xData && !histogram->yData && !histogram->longitData)
    bomb("no data selected for HISTOGRAM element", NULL);
  if (histogram->binSizeFactor<=0)
    bomb("bin_size_factor is non-positive for HISTOGRAM element", NULL);
  
  histogram->filename = compose_filename(histogram->filename, run->rootname);
  
  SDDS_HistogramSetup(histogram, SDDS_BINARY, 1, run->runfile, run->lattice, "set_up_histogram");
  histogram->initialized = 1;
  histogram->count = 0;
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


void reset_special_elements(LINE_LIST *beamline, long includeRF)
{
    ELEMENT_LIST *eptr;
    NIBEND *nibend; NISEPT *nisept;
    
    log_entry("reset_special_elements");

    eptr = &(beamline->elem);
    while (eptr) {
        switch (eptr->type) {
          case T_KICKER:
            ((KICKER*)eptr->p_elem)->fiducial_seen = 0;
            break;
          case T_NIBEND:
            nibend = (NIBEND*)eptr->p_elem;
            nibend->initialized = 0;
            break;
          case T_NISEPT:
            nisept = (NISEPT*)eptr->p_elem;
            nisept->fse_opt = 0;
            break;
          case T_TMCF:
            if (includeRF) {
              ((TMCF_MODE*)eptr->p_elem)->fiducial_part = NULL;
            }
            break;
          case T_CEPL:
            ((CE_PLATES*)eptr->p_elem)->fiducial_part = NULL;
            break;
          case T_TWPL:
            ((TW_PLATES*)eptr->p_elem)->fiducial_part = NULL;
            break;
          case T_TWLA:
            if (includeRF) {
              ((TW_LINAC*)eptr->p_elem)->fiducial_part = NULL;
            }
            break;
          case T_TWMTA:
            if (includeRF) {
              ((TWMTA*)eptr->p_elem)->fiducial_part = NULL;
            }
            break;
          case T_RFCA:
            if (includeRF) {
              ((RFCA*)eptr->p_elem)->fiducial_seen = 0;
            }
            break;
          case T_RFCW:
            if (includeRF) {
              ((RFCW*)eptr->p_elem)->rfca.fiducial_seen = 0;
            }
            break;
          case T_MODRF:
            if (includeRF) {
              ((MODRF*)eptr->p_elem)->fiducial_seen = 0;
            }
            break;
          case T_RFMODE:
            if (includeRF) {
              ((RFMODE*)eptr->p_elem)->initialized = 0;
              if (SDDS_IsActive(&((RFMODE*)eptr->p_elem)->SDDSrec))
                SDDS_Terminate(&((RFMODE*)eptr->p_elem)->SDDSrec);
            }
            break;
          case T_FRFMODE:
            if (includeRF)
              ((FRFMODE*)eptr->p_elem)->initialized = 0;
            break;
          case T_TRFMODE:
            if (includeRF) {
              ((TRFMODE*)eptr->p_elem)->initialized = 0;
              if (((TRFMODE*)eptr->p_elem)->fprec)
                fclose(((TRFMODE*)eptr->p_elem)->fprec);
            }
            break;
          case T_FTRFMODE:
            if (includeRF)
              ((FTRFMODE*)eptr->p_elem)->initialized = 0;
            break;
          case T_RAMPRF:
            if (includeRF) {
              ((RAMPRF*)eptr->p_elem)->Ts = 0;
              ((RAMPRF*)eptr->p_elem)->fiducial_seen = 0;
            }
            break;
          case T_KQUAD:
            ((KQUAD*)eptr->p_elem)->randomMultipoleData.randomized = 0;
            break;
          case T_KSEXT:
            ((KSEXT*)eptr->p_elem)->randomMultipoleData.randomized = 0;
            break;
          case T_MATR:
            if (includeRF)
              ((MATR*)eptr->p_elem)->fiducialSeen = 0;
            break;
          case T_EMATRIX:
            if (includeRF)
              ((EMATRIX*)eptr->p_elem)->fiducialSeen = 0;
            break;
          default:
            break;
            }
        eptr = eptr->succ;
        }
    log_exit("reset_special_elements");
    }

VMATRIX *stray_field_matrix(double length, double *lB, double *gB, double theta, long order, double p_central, 
                            void *Wi)
{
    /* factor in equation theta = CTHETA*B(T)*L(m)/(beta*gamma): */
#define CTHETA 5.8667907921396181e+02
    VMATRIX *M;
    double xkick, ykick, Bx, By;
    double rho;
#ifdef DEBUG_STRAY
    static FILE *fp=NULL;
    if (!fp) {
      fp = fopen_e("stray.erl", "w", 0);
      fprintf(fp, "SDDS1\n&column name=xKick type=double &end\n");
      fprintf(fp, "&column name=yKick type=double &end\n");
      fprintf(fp, "&column name=p type=double &end\n");
      fprintf(fp, "&column name=By type=double &end\n");
      fprintf(fp, "&column name=Bx type=double &end\n");
      fprintf(fp, "&data mode=ascii no_row_counts=1 &end\n");
    }
#endif

    if (Wi) {
      /* BLocalTotal = BLocal + Wi*GBlobal */
      MATRIX *Bg, *Bl;
      m_alloc(&Bg, 3, 1);
      m_alloc(&Bl, 3, 1);
      Bg->a[0][0] = gB[0];
      Bg->a[1][0] = gB[1];
      Bg->a[2][0] = gB[2];
      m_mult(Bl, (MATRIX*)Wi, Bg);
      Bx = Bl->a[0][0] + lB[0];
      By = Bl->a[1][0] + lB[1];
#ifdef DEBUG_STRAY
      m_show(Bg, "%10.3e", "Bg: \n", stdout);
      m_show(Wi, "%10.3e", "Wi: \n", stdout);
      m_show(Bl, "%10.3e", "Wi*Bg: \n", stdout);
      fprintf(stdout, "Bx = %e, By = %e\n", Bx, By);
#endif
      m_free(&Bg);
      m_free(&Bl);
    } else {
      if (gB[0] || gB[1] || gB[2]) 
        bomb("to use global stray fields, you must do a floor coordinate computation", NULL);
      Bx = lB[0];
      By = lB[1];
    }
    if (By) {
        rho = p_central/(CTHETA*By);
        xkick = asin(length/rho);
        }
    else
        xkick = 0;
    if (Bx) {
        rho = p_central/(CTHETA*Bx);
        ykick = -asin(length/rho);
        }
    else
        ykick = 0;
#ifdef DEBUG_STRAY
    fprintf(fp, "%le %le %le %le %le\n", xkick, ykick, p_central, By, Bx);
    fflush(fp);
#endif

    M = hvcorrector_matrix(length, xkick, ykick, 0.0, 0.0, 1.0, 1.0, 0, order);

    return(M);
    }


VMATRIX *rf_cavity_matrix(double length, double voltage, double frequency, double phase, 
                          double *P_central, long order, long end1Focus, long end2Focus,
                          char *bodyFocusModel, long change_p0, double Preference)
{
    VMATRIX *M, *Medge, *Mtot, *tmp;
    double *C, **R, dP, gamma, dgamma, dgammaMax;
    double cos_phase, sin_phase;
    double inverseF[2] = {0,0};
    long end, useSRSModel;

    M = tmalloc(sizeof(*M));
    M->order = 1;
    initialize_matrices(M, M->order);
    R = M->R;
    C = M->C;
    
    if (*P_central<=0) {
        fprintf(stdout, "error: P_central = %g\n", *P_central);
        fflush(stdout);
        abort();
        }

    C[4] = length;
    voltage /= 1e6;  /* convert to MV */
    sin_phase = sin(PI*phase/180.0);
    cos_phase = cos(PI*phase/180.0);
    dgamma = (dgammaMax=voltage/me_mev)*sin_phase;
    gamma = sqrt(sqr(*P_central)+1);
    dP    = sqrt(sqr(gamma+dgamma)-1) - *P_central;

    useSRSModel = 0;
    if (bodyFocusModel) {
      char *modelName[2] = { "none", "srs" };
      switch (match_string(bodyFocusModel, modelName, 2, 0)) {
      case 0:
        break;
      case 1:
        useSRSModel = 1;
        break;
      default:
        fprintf(stderr, "Error: bodyFocusModel=%s not understood for RFCA\n", bodyFocusModel);
        exit(1);
        break;
      }
    }

    if (!useSRSModel) {
      R[0][0] = R[2][2] = R[4][4] = 1;
      R[1][1] = R[3][3] = R[5][5] = 1/(1+dP/(*P_central));
      if (fabs(dP/(*P_central))>1e-14)
        R[0][1] = R[2][3] = length*(*P_central)/dP*log(1 + dP/(*P_central));
      else
        R[0][1] = R[2][3] = length;
    } else {
      /* note that Rosenzweig and Serafini use gamma in places
       * where they should probably use momentum, but I'll keep
       * their expressions for now.
       */
      double alpha, sin_alpha, gammaf;
      gammaf = gamma+dgamma;
      if (fabs(sin_phase)>1e-6)
        alpha = log(gammaf/gamma)/(2*SQRT2*sin_phase);
      else
        alpha = dgammaMax/gamma/(2*SQRT2);
      R[0][0] = R[2][2] = cos(alpha);
      R[1][1] = R[3][3] = R[0][0]*gamma/gammaf;
      R[0][1] = R[2][3] = 2*SQRT2*gamma*length/dgammaMax*(sin_alpha=sin(alpha));
      R[1][0] = R[3][2] = -sin_alpha*dgammaMax/(length*gammaf*2*SQRT2);
    }
    
    R[5][4] = (voltage/me_mev)*cos_phase/(gamma + dgamma)*(PIx2*frequency/c_mks);

    if (length && (end1Focus || end2Focus)) {
      if (end1Focus) {
        inverseF[0] = dgamma/(2*length*gamma);
        }
      if (end2Focus)
        inverseF[1] = -dgamma/(2*length*(gamma+dgamma));
      Medge = tmalloc(sizeof(*Medge));
      initialize_matrices(Medge, Medge->order = 1);
      Medge->R[0][0] = Medge->R[1][1] = Medge->R[2][2] = Medge->R[3][3] = 
        Medge->R[4][4] = Medge->R[5][5] = 1;
      Mtot = tmalloc(sizeof(*Mtot));
      initialize_matrices(Mtot, Mtot->order = 1);
      for (end=0; end<2; end++) {
        if (inverseF[end]) {
          Medge->R[1][0] = Medge->R[3][2] = -inverseF[end];
        } else
          continue;
        if (end==0) {
          concat_matrices(Mtot, M, Medge, CONCAT_EXCLUDE_S0);
        } else {
          concat_matrices(Mtot, Medge, M, CONCAT_EXCLUDE_S0);
        }
        tmp = Mtot;
        Mtot = M;
        M = tmp;
      }
      free_matrices(Medge);
      tfree(Medge);
      free_matrices(Mtot);
      tfree(Mtot);
    }

    if (change_p0)
      Preference = *P_central + dP;

    if (!change_p0) {
      /* The matrix implicitly included acceleration. 
       * Must change the reference momentum of the matrix to exclude this.
       */
      long i;
      Mtot = tmalloc(sizeof(*Mtot));
      Medge = tmalloc(sizeof(*Medge));
      initialize_matrices(Mtot, Mtot->order = 1);
      initialize_matrices(Medge, Medge->order = 1);
      Medge->R[0][0] = Medge->R[1][1] = Medge->R[2][2] = Medge->R[3][3] = 
        Medge->R[4][4] = Medge->R[5][5] = 1;
      Medge->C[5] = (*P_central+dP-Preference)/Preference;
      concat_matrices(Mtot, Medge, M, CONCAT_EXCLUDE_S0);
      tmp = Mtot;
      Mtot = M;
      M = tmp;
      free_matrices(Medge);
      free_matrices(Mtot);
      tfree(Medge);
      tfree(Mtot);
    }

    *P_central = Preference;

    return(M);
}

VMATRIX *twissTransformMatrix(TWISSELEMENT *twissWanted,
                              TWISS *twissInput)
{
  VMATRIX *M;
  double beta1, beta2, alpha1, alpha2;
  
  M = tmalloc(sizeof(*M));
  initialize_matrices(M, M->order=1);
  M->R[0][0] = M->R[1][1] = M->R[2][2] = M->R[3][3] = 
    M->R[4][4] = M->R[5][5] = 1;
  if (twissInput==NULL)
    return M;

  beta1 = twissInput->betax;
  beta2 = twissWanted->betax;
  alpha1 = twissInput->alphax;
  alpha2 = twissWanted->alphax;
  M->R[0][0] = beta2/sqrt(beta1*beta2);
  M->R[1][0] = (alpha1-alpha2)/sqrt(beta1*beta2);
  M->R[1][1] = beta1/sqrt(beta1*beta2);
  
  beta1 = twissInput->betay;
  beta2 = twissWanted->betay;
  alpha1 = twissInput->alphay;
  alpha2 = twissWanted->alphay;
  M->R[2][2] = beta2/sqrt(beta1*beta2);
  M->R[3][2] = (alpha1-alpha2)/sqrt(beta1*beta2);
  M->R[3][3] = beta1/sqrt(beta1*beta2);
  
  return M;
}

/* thin lens for light optics */
VMATRIX *lightThinLensMatrix(LTHINLENS *ltl)
{
  VMATRIX  *M;
  double *C, **R;
  long i;
  
  M = tmalloc(sizeof(*M));
  M->order = 1;
  initialize_matrices(M, M->order);
  R = M->R;
  C = M->C;

  for (i=0; i<6; i++) {
    C[i] = 0;
    R[i][i] = 1;
  }
  if (ltl->fx)
    R[1][0] = -1/ltl->fx;
  if (ltl->fy)
    R[3][2] = -1/ltl->fy;
  if (ltl->tilt)
    tilt_matrices(M, ltl->tilt);
  if (ltl->pitch)
    pitch_matrices(M, ltl->pitch);
  if (ltl->yaw)
    yaw_matrices(M, ltl->yaw);
  if (ltl->dx || ltl->dy || ltl->dz) 
    misalign_matrix(M, ltl->dx, ltl->dy, ltl->dz, 0.0);
  return M;
}

/* mirror for light optics */
VMATRIX *lightMirrorMatrix(LMIRROR *lm)
{
  VMATRIX  *M;
  double *C, **R;
  long i;
  
  M = tmalloc(sizeof(*M));
  M->order = 1;
  initialize_matrices(M, M->order);
  R = M->R;
  C = M->C;

  for (i=0; i<6; i++) {
    C[i] = 0;
    R[i][i] = 1;
  }
  if (lm->Rx)
    R[1][0] = -2/(lm->Rx*cos(lm->theta));
  if (lm->Ry)
    R[3][2] = -2/(lm->Ry/cos(lm->theta));
  if (lm->tilt)
    tilt_matrices(M, lm->tilt);
  if (lm->pitch)
    pitch_matrices(M, lm->pitch);
  if (lm->yaw)
    yaw_matrices(M, lm->yaw);
  if (lm->dx || lm->dy || lm->dz) 
    misalign_matrix(M, lm->dx, lm->dy, lm->dz, 2*PI-2*lm->theta);
  return M;
}

/* explicit matrix input from user */
VMATRIX *matrixFromExplicitMatrix(EMATRIX *emat, long order)
{
  long i, j, k;
  VMATRIX  *M;
  double *C, **R, ***T;

  if (emat->order!=0)
    order = emat->order;
  M = tmalloc(sizeof(*M));
  initialize_matrices(M, M->order=order);
  C = M->C;
  R = M->R;
  T = M->T;
  
  for (i=0; i<6; i++) 
    C[i] = emat->C[i];

  if (order>=1)
    for (i=0; i<6; i++) 
      for (j=0; j<6; j++) 
        R[i][j] = emat->R[i][j];

  if (order>=2)
    for (i=0; i<6; i++) 
      for (j=0; j<6; j++) 
        for (k=0; k<=j; k++)
          T[i][j][k] = emat->T[i][j][k];

  return M;
}


