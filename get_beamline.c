/* Copyright 1994 by Michael Borland and Argonne National Laboratory,
 * all rights reserved.
 */
/* routine: get_beamline()
 * purpose: read a mad-format lattice and return a pointer to a linked
 *          list for the beamline specified. 
 * 
 *	    It is assumed that the MAD-style input file is in
 *          the usual units of meters, radians, etc.  The 
 *          output has the same units.
 *
 * Michael Borland, 1989
 */
#include "mdb.h"
#include "track.h"
#include <ctype.h>
#include "match_string.h"

void show_elem(ELEMENT_LIST *eptr, long type);
void process_rename_request(char *s, char **name, long n_names);

static ELEMENT_LIST *elem;   /* elem: root of linked-list of ELEM structures */
static LINE_LIST *line;      /* line: root of linked-list of LINE structures */

void lfree(void *ptr) 
{
    static FILE *fpl = NULL;
    if (!fpl)
        fpl = fopen_e("free.log", "w", 0);
    fprintf(fpl, "%x freed\n", ptr);
    tfree(ptr);
    }

#define MAX_LINE_LENGTH 16384 
#define MAX_FILE_NESTING 10
LINE_LIST *get_beamline(char *madfile, char *use_beamline, double p_central, long echo)
{
  long type, i, i_elem;
  long occurence, iMad;
  static ELEMENT_LIST *eptr, *eptr1;
  static LINE_LIST *lptr;
  static long n_elems, n_lines;
  FILE *fp_mad[MAX_FILE_NESTING];
  char s[MAX_LINE_LENGTH], *cfgets(), *ptr, t[MAX_LINE_LENGTH];

  log_entry("get_beamline");

  if (madfile) {
#ifdef DEBUG
    fprintf(stdout, "reading from file %s\n", madfile);
    fflush(stdout);
#endif
    fp_mad[0] = fopen_e(madfile, "r", 0);
    iMad = 0;
    
    elem = tmalloc(sizeof(*elem));
    line = tmalloc(sizeof(*line));
    
    elem->pred = elem->succ = NULL;
    elem->name = NULL;
    eptr      = elem;
    n_elems   = 0;        /* number of physical elements in linked-list */
    line->pred = line->succ = NULL;
    line->name = NULL;
    lptr      = line;        
    n_lines   = 0;        /* number of line definitions in linked-list  */
    
    /* assemble linked-list of simple elements and a separate linked-list
       of fully expanded line definitions */
    while (iMad>=0) {
      while (cfgets(s, MAX_LINE_LENGTH, fp_mad[iMad])) {
        if (echo) {
          fprintf(stdout, "%s\n", s);
          fflush(stdout);
        }
        if (s[0]=='#' && strncmp(s, "#INCLUDE:", strlen("#INCLUDE:"))==0) {
          if (++iMad==MAX_FILE_NESTING) 
            bomb("files nested too deeply", NULL);
          ptr = get_token(s+strlen("#INCLUDE:"));
          if (echo) {
            fprintf(stdout, "reading file %s\n", ptr);
          }
          fp_mad[iMad] = fopen_e(ptr, "r", 0);
          continue;
        }
        strcpy(t, s);
        if ((type = tell_type(s, elem))==T_NODEF) {
          if (!is_blank(s))
            fprintf(stdout, "warning: no recognized statement on line: %s\n", t);
          fflush(stdout);
          continue;
        }
#ifdef DEBUG
        fprintf(stdout, "type code = %ld\n", type);
        fflush(stdout);
#endif
        if (type==T_RENAME)
          process_rename_request(s, entity_name, N_TYPES);
        else if (type==T_TITLE) {
          fgets(s, MAX_LINE_LENGTH, fp_mad[iMad]);
          compressString(s, " ");
          if (s[i=strlen(s)-1]==' ')
            s[i] = 0;
        }    
        else if (type==T_USE || type==T_RETURN) 
          break;
        else if (type==T_LINE) {
#ifdef DEBUG
          fprintf(stdout, "current element list is:\n");
          fflush(stdout);
          print_elem_list(stdout, elem);
#endif
          fill_line(line, n_lines, elem, n_elems, s);
          check_duplic_line(line, lptr, n_lines+1);
#ifdef DEBUG 
          fprintf(stdout, "\n****** expanded line %s:\n", lptr->name);
          fflush(stdout);
          print_line(stdout, lptr); 
#endif
          extend_line_list(&lptr);
          n_lines++;
        }
        else {
          if (type==T_ECOPY) {
#ifdef DEBUG
            fprintf(stdout, "copying existing element\n");
            fflush(stdout);
#endif
            strcpy(s, t);
            copy_named_element(eptr, s, elem);
          }
          else {
#ifdef DEBUG
            fprintf(stdout, "creating new element\n");
            fflush(stdout);
#endif
            fill_elem(eptr, s, type, fp_mad[iMad]); 
          }
          check_duplic_elem(&elem, &eptr, n_elems);
#ifdef DEBUG
          print_elem(stdout, elem);
#endif
          extend_elem_list(&eptr);
          n_elems++;
        }
      }
      fclose(fp_mad[iMad]);
      iMad--;
    }
    if (n_elems==0 || n_lines==0) {
      fprintf(stdout, "Insufficient (recognizable) data in file.\n");
      fflush(stdout);
      exit(1);
    }
    /* since the lists were being extended before it was known that
       the was another object to put in them, must eliminate references
       to the most recently added nodes. 
       */
    (eptr->pred)->succ = NULL;
    (lptr->pred)->succ = NULL;
    eptr = eptr->pred;
    lptr = lptr->pred;
  }
  else {
    s[0] = 0;
    type = T_NODEF;
  }
  

  if (type!=T_USE && use_beamline==NULL) {
    if (n_lines==0)
      bomb("no beam-line defined\n", NULL);
    fprintf(stdout, "no USE statement--will use line %s\n", lptr->name);
    fflush(stdout);
  }
  else {
    if (!use_beamline) {
      if ((ptr=get_token(s))==NULL) 
        bomb("no line named in USE statement", NULL);
    }
    else
      ptr = str_toupper(use_beamline);
    lptr = line;
    while (lptr) {
      if (strcmp(lptr->name, ptr)==0) 
        break;
      lptr = lptr->succ;
    }
    if (lptr==NULL) {
      fprintf(stdout, "no definition of beam-line %s\n", ptr);
      fflush(stdout);
      exit(1);
    }
  }

  /* these really aren't necessary, since I clear memory upon allocation */
  lptr->elem_recirc = lptr->elem_twiss = lptr->elast = NULL;
  lptr->twiss0 = NULL;
  lptr->matrix = NULL;

  /* go through and give occurence numbers to each element */
  eptr = &(lptr->elem);
  while (eptr) {
    eptr->occurence = 0;
    lptr->elast = eptr;
    eptr = eptr->succ;
  }

  eptr = &(lptr->elem);
  lptr->flags = 0;
  while (eptr) {
    eptr->Pref_input = eptr->Pref_output = p_central;
    if (eptr->occurence==0) {
      /* this is the first occurence of this element--go through and find any others */
      eptr->occurence = occurence = 1;
      eptr1 = eptr->succ;
      while (eptr1) {
        if (strcmp(eptr->name, eptr1->name)==0)
          eptr1->occurence = ++occurence;
        eptr1 = eptr1->succ;
      }
    }
    if (eptr->type==T_SREFFECTS)
      lptr->flags |= BEAMLINE_TWISS_WANTED;
    eptr = eptr->succ;
  }

  compute_end_positions(lptr);
  
  return(lptr);
}

double compute_end_positions(LINE_LIST *lptr) 
{
    double z, l, theta, z_recirc;
    static ELEMENT_LIST *eptr, *eptr1;
    long i_elem;
    
    /* use length data to establish z coordinates at end of each element */
    /* also check for duplicate recirculation elements and set occurence numbers to 0 */
    eptr = &(lptr->elem);
    z = z_recirc = 0;
    theta = 0;
    i_elem = 0;
    do {
        if (!(entity_description[eptr->type].flags&HAS_LENGTH))
            l = 0;
        else
            l = (*((double*)eptr->p_elem));
        if (eptr->type==T_SBEN || eptr->type==T_RBEN)
            theta += ((BEND*)eptr->p_elem)->angle;
        else if (eptr->type==T_KSBEND)
            theta += ((KSBEND*)eptr->p_elem)->angle;
        else if (eptr->type==T_NIBEND)
            theta += ((NIBEND*)eptr->p_elem)->angle;
        else if (eptr->type==T_NISEPT)
            theta += ((NISEPT*)eptr->p_elem)->angle;
        else if (eptr->type==T_CSBEND)
            theta += ((CSBEND*)eptr->p_elem)->angle;
        else if (eptr->type==T_RECIRC) {
            if (lptr->elem_recirc)
                bomb("multiple recirculation (RECIRC) elements in beamline--this doesn't make sense", NULL);
            lptr->elem_recirc = eptr;
            lptr->i_recirc = i_elem;
            z_recirc = z;
            }
        if (l<0) {
            fprintf(stdout, "warning(1): element %s has negative length = %e\n", eptr->name, l);
            fflush(stdout);
          }
        eptr->end_pos = z + l;
        eptr->end_theta = theta ;
        z = eptr->end_pos;
        i_elem++;
        } while ((eptr=eptr->succ));

    lptr->revolution_length = z - z_recirc;
    return lptr->revolution_length;
    }

void show_elem(ELEMENT_LIST *eptr, long type)
{
    long j;
    char *ptr;
    PARAMETER *parameter;

    log_entry("show_elem");

    parameter = entity_description[type].parameter;
    fprintf(stdout,  "%s %s at z=%em:\n", 
        entity_name[type], eptr->name, eptr->end_pos);
    fflush(stdout);
    for (j=0; j<entity_description[type].n_params; j++) {
        switch (parameter[j].type) {
            case IS_DOUBLE:
                fprintf(stdout,  "    %s = %.16e with offset %ld\n", 
                    parameter[j].name, 
                    *(double*)(eptr->p_elem+parameter[j].offset),
                    parameter[j].offset);
                fflush(stdout);
                break;
            case IS_LONG:
                fprintf(stdout,  "    %s = %ld with offset %ld\n", 
                    parameter[j].name, 
                    *(long *)(eptr->p_elem+parameter[j].offset),
                    parameter[j].offset);
                fflush(stdout);
                break;
            case IS_STRING:
                if ((ptr = *(char**)(eptr->p_elem+parameter[j].offset))) {
                    fprintf(stdout,  "    %s = %e\n", parameter[j].name, ptr);
                    fflush(stdout);
                  }
                else {
                    fprintf(stdout,  "    %s = %e\n", parameter[j].name, ptr);
                    fflush(stdout);
                  }
                break;
            }
        }
    log_exit("show_elem");
    }

void free_elements(ELEMENT_LIST *elemlist)
{
    ELEMENT_LIST *eptr;

    log_entry("free_elements");

    if (elemlist) {
        eptr= elemlist;
        }
    else  {
        eptr = elem;
        elem = NULL;
#ifdef DEBUG
        fprintf(stdout, "freeing elements in main list\n");
        fflush(stdout);
#endif
        }
    while (eptr) {
#ifdef DEBUG
        fprintf(stdout, "freeing memory for element %s of type %s\n", 
                eptr->name?eptr->name:"NULL",
                (eptr->type>=0 && eptr->type<N_TYPES)?entity_name[eptr->type]:"NULL");
        fflush(stdout);
#endif
        if (eptr->type==T_WATCH) {
            WATCH *wptr;
            if ((wptr = (WATCH*)eptr->p_elem)) {
                if (wptr->initialized && !SDDS_Terminate(&wptr->SDDS_table)) {
                    SDDS_SetError("Problem terminate watch-point SDDS file (free_elements)");
                    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
                    }
                }
            }
#ifdef DEBUG
        fprintf(stdout, "pointers: p_elem = %x   name = %x   matrix = %x\n",
            eptr->p_elem, eptr->name, eptr->matrix);
        fflush(stdout);
        fprintf(stdout, "          pred = %x     succ = %x  \n",
            eptr->pred, eptr->succ);
        fflush(stdout);
#endif
        tfree(eptr->p_elem);
        eptr->p_elem = NULL;
        tfree(eptr->name);
        eptr->name = NULL;
        if (entity_description[eptr->type].flags&HAS_MATRIX && eptr->matrix) {
            free_matrices(eptr->matrix);
            free(eptr->matrix);
            eptr->matrix = NULL;
            }
        if (eptr->accumMatrix) {
          free_matrices(eptr->accumMatrix);
          free(eptr->accumMatrix);
          eptr->accumMatrix = NULL;
        }
        if (eptr->succ) {
            eptr = eptr->succ;
            free(eptr->pred);
            eptr->pred = NULL;
            }
        else {
            free(eptr);
            break;
            }
        }
    log_exit("free_elements");
    }

void free_beamlines(LINE_LIST *beamline)
{
    LINE_LIST *lptr;

    log_entry("free_beamlines");

    if (beamline) {
        lptr = beamline;
        }
    else {
        lptr = line;
        line = NULL;
#ifdef DEBUG
        fprintf(stdout, "freeing main beamline list\n", NULL);
        fflush(stdout);
#endif
        }
    while (lptr) {
#ifdef DEBUG
        fprintf(stdout, "*************************\nfreeing memory for beamline %s with %ld elements\n",
            lptr->name?lptr->name:"NULL", lptr->n_elems);
        fflush(stdout);
        fprintf(stdout, "pointers:   name = %x    succ = %x   pred = %x\n",
            lptr->name, lptr->succ, lptr->pred);
        fflush(stdout);
#endif
        if (lptr->definition) {
            tfree(lptr->definition);
            lptr->definition = NULL;
            }
        if (lptr->name) {
            tfree(lptr->name);
            lptr->name = NULL;
            }
        if (lptr->n_elems) {
            free_elements((lptr->elem).succ);
            /* should free name etc. for lptr->elem also */
            lptr->n_elems = 0;
            lptr->flags = 0;
            }
        if (lptr->succ) {
            lptr = lptr->succ;
            tfree(lptr->pred);
            lptr->pred = NULL;
            }
        else {
            tfree(lptr);
            break;
            }
        }
    log_exit("free_beamlines");
    }

void delete_matrix_data(LINE_LIST *beamline)
{
    LINE_LIST *lptr;
    ELEMENT_LIST *eptr;

    log_entry("delete_matrix_data");

    if (beamline) {
        lptr = beamline;
        }
    else {
        lptr = line;
        }
    while (lptr) {
        if (lptr->n_elems) {
            eptr = &(lptr->elem);
            while (eptr) {
                if (entity_description[eptr->type].flags&HAS_MATRIX && eptr->matrix) {
                    free_matrices(eptr->matrix);
                    tfree(eptr->matrix);
                    eptr->matrix = NULL;
                    }
                if (eptr->accumMatrix) {
                  free_matrices(eptr->accumMatrix);
                  tfree(eptr->accumMatrix);
                  eptr->accumMatrix = NULL;
                }
                eptr = eptr->succ;
                }
            lptr->n_elems = 0;
            lptr->flags = 0;
            }
        lptr = lptr->succ;
        }
    log_exit("delete_matrix_data");
    }


/* routine: do_save_lattice()
 * purpose: save the element and beamline definitions to a file
 * 
 * Michael Borland, 1991
 */
#include "save_lattice.h"

void do_save_lattice(NAMELIST_TEXT *nltext, RUN *run, LINE_LIST *beamline)
{
    FILE *fp;
    ELEMENT_LIST *eptr;
    LINE_LIST *lptr;
    long j;
    double dvalue;
    long lvalue;
    char *ptr;
    PARAMETER *parameter;
    char s[16384], t[1024], name[1024];

    log_entry("do_save_lattice");

    /* process the namelist text */
    set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
    set_print_namelist_flags(0);
    process_namelist(&save_lattice, nltext);
    print_namelist(stdout, &save_lattice);

    /* check for valid data */
    if (filename==NULL)
        bomb("no filename given to save lattice to", NULL);
    if (str_in(filename, "%s"))
        filename = compose_filename(filename, run->rootname);
    fp = fopen_e(filename, "w", FOPEN_INFORM_OF_OPEN);

    eptr = elem;
    while (eptr) {
        parameter = entity_description[eptr->type].parameter;
        if (strpbrk(eptr->name, ":.,/-_+ "))
            sprintf(name, "\"%s\"", eptr->name);
        else
            strcpy(name, eptr->name);
        sprintf(s, "%s: %s,", name, entity_name[eptr->type]);
        for (j=0; j<entity_description[eptr->type].n_params; j++) {
            switch (parameter[j].type) {
                case IS_DOUBLE:
                    dvalue = *(double*)(eptr->p_elem+parameter[j].offset);
                    if (dvalue!=parameter[j].number) {
                        /* value is not the default, so add to output */
                        sprintf(t, "%s=%.16g", parameter[j].name, dvalue);
                        strcat(s, t);
                        if (j!=entity_description[eptr->type].n_params-1)
                            strcat(s, ",");
                        }
                    break;
                case IS_LONG:
                    lvalue = *(long *)(eptr->p_elem+parameter[j].offset);
                    if (lvalue!=parameter[j].integer) {
                        /* value is not the default, so add to output */
                        sprintf(t, "%s=%ld", parameter[j].name, lvalue);
                        strcat(s, t);
                        if (j!=entity_description[eptr->type].n_params-1)
                            strcat(s, ",");
                        }
                    break;
                case IS_STRING:
                    if ((ptr = *(char**)(eptr->p_elem+parameter[j].offset))) {
                        sprintf(t, "%s=\"%s\"", parameter[j].name, ptr);
                        strcat(s, t);
                        if (j!=entity_description[eptr->type].n_params-1)
                            strcat(s, ",");
                        }                    
                    break;
                }
            }
        if (s[j=strlen(s)-1]==',')
            s[j] = 0;
        print_with_continuation(fp, s, 79);
        eptr = eptr->succ;
        }
    lptr = line;
    while (lptr) {
        print_with_continuation(fp, lptr->definition, 79);
        lptr = lptr->succ;
        }

    if (beamline && beamline->name)
        fprintf(fp, "USE,%s\n", beamline->name);

    fprintf(fp, "RETURN\n");
    fclose(fp);

    log_exit("do_save_lattice");
    }

void print_with_continuation(FILE *fp, char *s, long endcol)
{
    char c, *ptr;
    long l;

    while ((l=strlen(s))) {
        if (l>endcol) {
          ptr = s+endcol-2;
          while (ptr!=s && *ptr!=',')
            ptr--;
          if (ptr==s)
            c = *(ptr = s + endcol - 1);
          else {
            ptr++;
            c = *ptr;
          }
          *ptr = 0;
          fputs(s, fp);
          fputs("&\n", fp);
          s = ptr;
          *ptr = c;
        }
        else {
          fputs(s, fp);
          fputc('\n', fp);
          log_exit("print_with_continuation");
          return;
        }
      }
  }

void change_defined_parameter_values(char **elem_name, long *param_number, long *type, double *value, long n_elems)
{
    ELEMENT_LIST *eptr;
    char *p_elem;
    long i_elem, elem_type, data_type, param;

    log_entry("change_defined_parameter_values");

    for (i_elem=0; i_elem<n_elems; i_elem++) {
        eptr = NULL;
        elem_type = type[i_elem];
        param     = param_number[i_elem];
        data_type = entity_description[elem_type].parameter[param].type;
        while (find_element(elem_name[i_elem], &eptr, elem)) {
            p_elem = eptr->p_elem;
            switch (data_type) {
                case IS_DOUBLE:
                    *((double*)(p_elem+entity_description[elem_type].parameter[param].offset)) = value[i_elem];
#if DEBUG
                    fprintf(stdout, "   changing parameter %s of %s #%ld to %e\n",
                           entity_description[elem_type].parameter[param].name,
                           eptr->name, eptr->occurence,
                           *((double*)(p_elem+entity_description[elem_type].parameter[param].offset)));
                    fflush(stdout);
#endif
                    break;
                case IS_LONG:
                    *((long*)(p_elem+entity_description[elem_type].parameter[param].offset)) = 
                      nearestInteger(value[i_elem]);
#if DEBUG
                    fprintf(stdout, "   changing parameter %s of %s #%ld to %ld\n",
                           entity_description[elem_type].parameter[param].name,
                           eptr->name, eptr->occurence,
                           *((long*)(p_elem+entity_description[elem_type].parameter[param].offset)));
                    fflush(stdout);
#endif
                    break;
                case IS_STRING:
                default:
                    bomb("unknown/invalid variable quantity", NULL);
                    exit(1);
                }
            }
        }
    log_exit("change_defined_parameter_values");
    }

void change_defined_parameter(char *elem_name, long param, long elem_type, 
                              double value, char *valueString, unsigned long mode)
{
  ELEMENT_LIST *eptr;
  char *p_elem;
  long data_type;

  log_entry("change_defined_parameter");

  eptr = NULL;
  data_type = entity_description[elem_type].parameter[param].type;
  if (mode&LOAD_FLAG_IGNORE)
    return;
  while (find_element(elem_name, &eptr, elem)) {
    p_elem = eptr->p_elem;
    switch (data_type) {
    case IS_DOUBLE:
      if (valueString) {
        if (!sscanf(valueString, "%lf", &value)) {
          fprintf(stdout, "Error (change_defined_parameter): unable to scan double from \"%s\"\n", valueString);
          fflush(stdout);
          exit(1);
        }
      }
      if (mode&LOAD_FLAG_VERBOSE)
        fprintf(stdout, "Changing definition (mode %s) %s.%s from %e to ", 
                (mode&LOAD_FLAG_ABSOLUTE)?"absolute":
                ((mode&LOAD_FLAG_DIFFERENTIAL)?"differential":
                 (mode&LOAD_FLAG_FRACTIONAL)?"fractional":"unknown"),
                elem_name, entity_description[elem_type].parameter[param].name,
                *((double*)(p_elem+entity_description[elem_type].parameter[param].offset)));
        fflush(stdout);
      if (mode&LOAD_FLAG_ABSOLUTE)
        *((double*)(p_elem+entity_description[elem_type].parameter[param].offset)) = value;
      else if (mode&LOAD_FLAG_DIFFERENTIAL)
        *((double*)(p_elem+entity_description[elem_type].parameter[param].offset)) += value;
      else if (mode&LOAD_FLAG_FRACTIONAL)
        *((double*)(p_elem+entity_description[elem_type].parameter[param].offset)) *= 1+value;
      if (mode&LOAD_FLAG_VERBOSE)
        fprintf(stdout, "%e\n", 
                *((double*)(p_elem+entity_description[elem_type].parameter[param].offset)));
        fflush(stdout);
      break;
    case IS_LONG:
      if (valueString) {
        if (!sscanf(valueString, "%lf", &value)) {
          fprintf(stdout, "Error (change_defined_parameter): unable to scan double from \"%s\"\n", valueString);
          fflush(stdout);
          exit(1);
        }
      }
      if (mode&LOAD_FLAG_VERBOSE)
        fprintf(stdout, "Changing definition (mode %s) %s.%s from %ld to ",
                (mode&LOAD_FLAG_ABSOLUTE)?"absolute":
                ((mode&LOAD_FLAG_DIFFERENTIAL)?"differential":
                 (mode&LOAD_FLAG_FRACTIONAL)?"fractional":"unknown"),
                elem_name, entity_description[elem_type].parameter[param].name,
                *((long*)(p_elem+entity_description[elem_type].parameter[param].offset)));
        fflush(stdout);
      if (mode&LOAD_FLAG_ABSOLUTE)
        *((long*)(p_elem+entity_description[elem_type].parameter[param].offset)) = 
          nearestInteger(value);
      else if (mode&LOAD_FLAG_DIFFERENTIAL)
        *((long*)(p_elem+entity_description[elem_type].parameter[param].offset)) += 
          nearestInteger(value);
      else if (mode&LOAD_FLAG_FRACTIONAL)
        *((long*)(p_elem+entity_description[elem_type].parameter[param].offset)) *= 1+value;
      if (mode&LOAD_FLAG_VERBOSE)
        fprintf(stdout, "%ld\n",
                *((long*)(p_elem+entity_description[elem_type].parameter[param].offset)));
        fflush(stdout);
      break;
    case IS_STRING:
      if (mode&LOAD_FLAG_VERBOSE)
        fprintf(stdout, "Changing definition %s.%s from %s to %s\n",
                elem_name, entity_description[elem_type].parameter[param].name,
                *((char**)(p_elem+entity_description[elem_type].parameter[param].offset)),
                valueString);
        fflush(stdout);
      if (!SDDS_CopyString(((char**)(p_elem+entity_description[elem_type].parameter[param].offset)), 
                           valueString)) {
        fprintf(stdout, "Error (change_defined_parameter): unable to copy string parameter value\n");
        fflush(stdout);
        exit(1);
      }
      break;
    default:
      bomb("unknown/invalid variable quantity", NULL);
      exit(1);
    }
  }
  log_exit("change_defined_parameter");
}

void process_rename_request(char *s, char **name, long n_names)
{
    long i;
    char *old, *new, *ptr;

    log_entry("process_rename_request");
    if (!(ptr=strchr(s, '='))) 
        bomb("invalid syntax for RENAME", NULL);
    *ptr++ = 0;
    old = s;
    str_toupper(trim_spaces(old = s));
    while (*old==',' || *old==' ')
        old++;
    str_toupper(trim_spaces(new = ptr));
    if (match_string(new, name, n_names, EXACT_MATCH)>=0) {
        fprintf(stdout, "error: can't rename to name %s--already exists\n", new);
        fflush(stdout);
        exit(1);
        }
    if ((i=match_string(old, name, n_names, EXACT_MATCH))<0) {
        fprintf(stdout, "error: can't rename %s to %s--%s not recognized\n", old, new, old);
        fflush(stdout);
        exit(1);
        }
    fprintf(stdout, "warning: element %s now known as %s\n", old, new);
    fflush(stdout);
    cp_str(name+i, new);
    log_exit("process_rename_request");
    }

    
long nearestInteger(double value)
{
  if (value<0)
    return -1*((long)(-value+0.5));
  return (long)(value+0.5);
}

