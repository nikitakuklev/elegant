/* Copyright 1994 by Michael Borland and Argonne National Laboratory,
 * all rights reserved.
 */
/* routine: print_line()
 * purpose: print a beam-line with magnet parameters
 * 
 * Michael Borland, 1989
 */
#include "track.h"
#include "mdb.h"

void print_line(FILE *fp, LINE_LIST *lptr)
{
    ELEMENT_LIST *eptr;
    long i, j;
    char *ptr;
    PARAMETER *parameter;

    log_entry("print_line");

    eptr = &(lptr->elem);
    fprintf(fp, "line name: %s\n", lptr->name);
    fprintf(fp, "line has %ld elements\n", lptr->n_elems);
    
    for (i=0; i<lptr->n_elems; i++) {
        parameter = entity_description[eptr->type].parameter;
        fprintf(fp, "%s %s at z=%em, theta=%e, part of %s:\n", 
            entity_name[eptr->type], eptr->name, eptr->end_pos,
            eptr->end_theta, (eptr->part_of?eptr->part_of:"?"));
        for (j=0; j<entity_description[eptr->type].n_params; j++) {
            switch (parameter[j].type) {
                case IS_DOUBLE:
                    fprintf(fp, "    %s = %.15g %s\n", 
                        parameter[j].name, 
                        *(double*)(eptr->p_elem+parameter[j].offset),
                        parameter[j].unit);
                    break;
                case IS_LONG:
                    fprintf(fp, "    %s = %ld %s\n", 
                        parameter[j].name, 
                        *(long *)(eptr->p_elem+parameter[j].offset),
                        parameter[j].unit);
                    break;
                case IS_STRING:
                    if ((ptr = *(char**)(eptr->p_elem+parameter[j].offset)))
                        fprintf(fp, "    %s = \"%s\"\n", parameter[j].name, ptr);
                    else
                        fprintf(fp, "    %s = \"\"\n", parameter[j].name);
                    break;
                }
            }
        eptr = eptr->succ;
        if (eptr==NULL && (i+1)!=lptr->n_elems) {
            fputs("line list ends prematurely", stderr);
            exit(1);
            }
        }
    log_exit("print_line");
    }

void print_elem(FILE *fp, ELEMENT_LIST *eptr)
{
    char *ptr;
    PARAMETER *parameter;
    long j;

    log_entry("print_elem");

    parameter = entity_description[eptr->type].parameter;
    fprintf(fp, "%s %s at z=%em, theta=%e:\n", 
        entity_name[eptr->type], eptr->name, eptr->end_pos, eptr->end_theta);
    for (j=0; j<entity_description[eptr->type].n_params; j++) {
        switch (parameter[j].type) {
            case IS_DOUBLE:
                fprintf(fp, "    %s = %.15g %s\n", 
                    parameter[j].name, 
                    *(double*)(eptr->p_elem+parameter[j].offset),
                    parameter[j].unit);
                break;
            case IS_LONG:
                fprintf(fp, "    %s = %ld %s\n", 
                    parameter[j].name, 
                    *(long *)(eptr->p_elem+parameter[j].offset),
                    parameter[j].unit);
                break;
            case IS_STRING:
                if ((ptr = *(char**)(eptr->p_elem+parameter[j].offset)))
                    fprintf(fp, "    %s = \"%s\"\n", parameter[j].name, ptr);
                else
                    fprintf(fp, "    %s = \"\"\n", parameter[j].name);
                break;
            }
        }
    log_exit("print_elem");
    }

void print_elem_list(FILE *fp, ELEMENT_LIST *eptr)
{
    long j;
    char *ptr;
    PARAMETER *parameter;

    log_entry("print_elem_list");

    while (eptr) {
        parameter = entity_description[eptr->type].parameter;
        fprintf(fp, "%s %s at z=%em, theta=%e:\n", 
            entity_name[eptr->type], eptr->name, eptr->end_pos, eptr->end_theta);
        if (eptr->type==T_MATR) {
            fprintf(fp, "    ORDER = %ld\n", eptr->matrix->order);
            }
        else {
            for (j=0; j<entity_description[eptr->type].n_params; j++) {
                switch (parameter[j].type) {
                    case IS_DOUBLE:
                        fprintf(fp, "    %s = %.15g %s\n", 
                            parameter[j].name, 
                            *(double*)(eptr->p_elem+parameter[j].offset),
                            parameter[j].unit);
                        break;
                    case IS_LONG:
                        fprintf(fp, "    %s = %ld %s\n", 
                            parameter[j].name, 
                            *(long *)(eptr->p_elem+parameter[j].offset),
                            parameter[j].unit);
                        break;
                    case IS_STRING:
                        if ((ptr = *(char**)(eptr->p_elem+parameter[j].offset)))
                            fprintf(fp, "    %s = \"%s\"\n", parameter[j].name, ptr);
                        else
                            fprintf(fp, "    %s = %e\n", parameter[j].name, ptr);
                        break;
                    }
                }
            }
        eptr = eptr->succ;
        }
    log_exit("print_elem_list");
    }

void print_elem_names(FILE *fp, ELEMENT_LIST *eptr)
{

    while (eptr) {
      if (eptr->name)
        fprintf(fp, "%s\n", eptr->name);
      eptr = eptr->succ;
    }
  }

