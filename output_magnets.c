/* Copyright 1994 by Michael Borland and Argonne National Laboratory,
 * all rights reserved.
 */
/* routine: output_magnets
 * purpose: write mpl format data to give plot of magnets.
 *
 * Michael Borland, 1988, 1991
 */
#include "mdb.h"
#include "track.h"

void output_magnets(char *filename, char *line_name, LINE_LIST *beamline)
{
    ELEMENT_LIST *eptr;
    QUAD  *qptr; BEND  *bptr; 
    KQUAD *kqptr; KSBEND *kbptr; CSBEND *cbptr;
    long n_points, iPhase;
    double start, end, total_length, dz, value;
    FILE *fpm;

    log_entry("output_magnets");

    total_length = 0;
    eptr = &(beamline->elem);
    n_points = 0;
    while (eptr) {
        switch (eptr->type) {
            case T_QUAD: case T_KQUAD:
                n_points += 5;
                break;
            case T_RBEN:  case T_SBEN: case T_KSBEND: case T_CSBEND:
                n_points += 6;
                break;
            case T_SEXT: case T_KSEXT:
                n_points += 5;
                break;
            case T_HCOR:
            case T_VCOR:
            case T_HVCOR:
                n_points += 3;
                break;
            case T_DRIF:
                n_points += 1;
                break;
            case T_HMON:
            case T_VMON:
                n_points += 4;
                break;
            case T_MONI:
                n_points += 5;
                break;
            case T_MULT:
                n_points += 7;
                break;
            case T_MARK:    /* zero-length drift */
                break;
            case T_RFCA: case T_TWLA: case T_RAMPRF:
            case T_MODRF:
                n_points += 9;
                break;
            default:
                if (entity_description[eptr->type].flags&HAS_LENGTH)
                    n_points += 1;
                else
                    fprintf(stderr, "warning: element %s, entity type (%s) will not appear in magnet output.\n",
                        eptr->name, entity_name[eptr->type]);
                break;
            }
        total_length = eptr->end_pos;
        eptr = eptr->succ;
        }

    fpm = fopen_e(filename, "w", FOPEN_SAVE_IF_EXISTS);

    start = end = 0;
    n_points += 1;
    fprintf(fpm, "SDDS1\n&description text=\"magnet layout for beamline %s\" &end\n", line_name);
    fprintf(fpm, "&column name=ElementName, type=string &end\n");
    fprintf(fpm, "&column name=s, units=m, type=double &end\n&column name=Profile, type=double &end\n");
    fprintf(fpm, "&data mode=ascii &end\n%ld\n", n_points);

    eptr = &(beamline->elem);
    fprintf(fpm, "_BEGIN_ 0 0\n");
    while (eptr) {
         switch (eptr->type) {
            case T_QUAD:
                qptr = (QUAD*)eptr->p_elem;
                fprintf(fpm, "%s %e  %d\n", eptr->name, start, SIGN(qptr->k1));
                end = start+qptr->length;
                fprintf(fpm, "%s %e  %d\n%s %e  0\n%s %e 0 %\n%s %e 0\n", eptr->name, 
                        end, SIGN(qptr->k1), eptr->name, end, eptr->name, start, eptr->name, end);
                start = end;
                break;
            case T_KQUAD:
                kqptr = (KQUAD*)eptr->p_elem;
                if (kqptr->bore)
                    value = kqptr->B;
                else
                    value = kqptr->k1;
                fprintf(fpm, "%s %e  %d\n", eptr->name, start, SIGN(kqptr->k1));
                end = start+kqptr->length;
                fprintf(fpm, "%s %e  %d\n%s %e  0\n%s %e 0 %\n%s %e 0\n", eptr->name, 
                        end, SIGN(kqptr->k1), eptr->name, end, eptr->name, start, eptr->name, end);
                start = end;
                break;
            case T_RBEN:  case T_SBEN:
                bptr = (BEND*)eptr->p_elem;
                end  = start+bptr->length;
                if (bptr->angle>0)
                    fprintf(fpm, 
                            "%s %e .33333333\n%s %e .33333333\n%s %e  0\n%s %e  0\n%s %e  0\n%s %e 0\n",
                            eptr->name, start, eptr->name, end, eptr->name, end, eptr->name, start, 
                            eptr->name, start, eptr->name, end);
                else if (bptr->angle<0)
                    fprintf(fpm, 
                            "%s %e -.33333333\n%s %e -.33333333\n%s %e  0\n%s %e  0\n%s %e  0\n%s %e 0\n",
                            eptr->name, start, eptr->name, end, eptr->name, end, eptr->name, start, 
                            eptr->name, start, eptr->name, end);
                start = end;
                break;
            case T_KSBEND:
                kbptr = (KSBEND*)eptr->p_elem;
                end  = start+kbptr->length;
                if (kbptr->angle>0)
                    fprintf(fpm, 
                            "%s %e .33333333\n%s %e .33333333\n%s %e  0\n%s %e  0\n%s %e  0\n%s %e 0\n",
                            eptr->name, start, eptr->name, end, eptr->name, end, eptr->name, start, 
                            eptr->name, start, eptr->name, end);
                else if (kbptr->angle<0)
                    fprintf(fpm, 
                            "%s %e -.33333333\n%s %e -.33333333\n%s %e  0\n%s %e  0\n%s %e  0\n%s %e 0\n",
                            eptr->name, start, eptr->name, end, eptr->name, end, eptr->name, start, 
                            eptr->name, start, eptr->name, end);
                start = end;
                break;
            case T_SEXT:
                end = start+((SEXT*)eptr->p_elem)->length;
                fprintf(fpm, "%s %e  .5\n%s %e .5\n%s %e 0\n%s %e 0\n%s %e 0\n", 
                        eptr->name, start, eptr->name, end, eptr->name, end, eptr->name, start, 
                        eptr->name, end);
                start = end;
                break;
            case T_KSEXT:
                end = start+((KSEXT*)eptr->p_elem)->length;
                fprintf(fpm, "%s %e  .5\n%s %e .5\n%s %e 0\n%s %e 0\n%s %e 0\n", 
                        eptr->name, start, eptr->name, end, eptr->name, end, eptr->name, start, 
                        eptr->name, end);
                start = end;
                break;
            case T_HCOR:
                end    = start+((HCOR*)eptr->p_elem)->length;
                fprintf(fpm, "%s %e .25\n%s %e .25\n%s %e 0\n", 
                        eptr->name, start, eptr->name, end, eptr->name, end);
                start = end;
                break;
            case T_VCOR:
                end    = start+((VCOR*)eptr->p_elem)->length;
                fprintf(fpm, "%s %e -.25\n%s %e -.25\n%s %e 0\n", 
                        eptr->name, start, eptr->name, end, eptr->name, end);
                start = end;
                break;
            case T_HVCOR:
                end    = start+((HVCOR*)eptr->p_elem)->length;
                fprintf(fpm, "%s %e .25\n%s %e -.25\n%s %e 0\n",
                        eptr->name, start, eptr->name, end, eptr->name, end);
                start = end;
                break;
            case T_DRIF:
                start = (end = start+((DRIFT*)eptr->p_elem)->length);
                fprintf(fpm, "%s %e  0\n", eptr->name, end);
                break;
            case T_HMON:
                dz = ((HMON*)eptr->p_elem)->length/2;
                fprintf(fpm, "%s %e 0.125\n%s %e 0\n%s %e 0\n%s %e 0\n",
                        eptr->name, start+dz, eptr->name, start+2*dz, 
                        eptr->name, start+dz, eptr->name, start+2*dz);
                start += 2*dz;
                break;
            case T_VMON:
                dz = ((VMON*)eptr->p_elem)->length/2;
                fprintf(fpm, "%s %e -0.125\n%s %e 0\n%s %e 0\n%s %e 0\n",
                        eptr->name, start+dz, eptr->name, start+2*dz, 
                        eptr->name, start+dz, eptr->name, start+2*dz);
                start += 2*dz;
                break;
            case T_MONI:
                dz = ((MONI*)eptr->p_elem)->length/2;
                fprintf(fpm, "%s %e 0.125\n%s %e 0\n%s %e -0.125\n%s %e 0\n%s %e 0\n",
                        eptr->name, start+dz, eptr->name, start+2*dz, eptr->name, start+dz, 
                        eptr->name, start, eptr->name, start+2*dz);
                start += 2*dz;
                break;
            case T_MULT:
                dz = ((MULT*)eptr->p_elem)->length/3;
                fprintf(fpm, "%s %e 0.6666\n%s %e 0.6666\n%s %e 0\n%s %e -0.6666\n%s %e -0.6666\n%s %e 0\n%s %e 0\n",
                        eptr->name, start+dz, eptr->name, start+2*dz, eptr->name, start+3*dz,
                        eptr->name, start+2*dz, eptr->name, start+dz, eptr->name, start, eptr->name, start+3*dz);
                start += 3*dz;
                break;
            case T_MARK:    /* zero-length drift */
                break;
            case T_CSBEND:
                cbptr = (CSBEND*)eptr->p_elem;
                end  = start+cbptr->length;
                if (cbptr->angle>0)
                    fprintf(fpm, 
                        "%s %e .33333333\n%s %e .33333333\n%s %e  0\n%s %e  0\n%s %e  0\n%s %e 0\n",
                            eptr->name, start, eptr->name, end, eptr->name, end, 
                            eptr->name, start, eptr->name, start, eptr->name, end);
                else if (cbptr->angle<0) 
                    fprintf(fpm, 
                        "%s %e -.33333333\n%s %e -.33333333\n%s %e  0\n%s %e  0\n%s %e  0\n%s %e 0\n",
                            eptr->name, start, eptr->name, end, eptr->name, end, 
                            eptr->name, start, eptr->name, start, eptr->name, end);
                start = end;
                break;
            case T_RFCA: case T_TWLA: case T_RAMPRF:
            case T_MODRF:
                dz = ((DRIFT*)eptr->p_elem)->length;
                dz /= 8;
                for (iPhase=0; iPhase<9; iPhase++) {
                  fprintf(fpm, "%s %e %e\n",
                          eptr->name, start+dz*iPhase,
                          0.5*sin((iPhase/8.0)*PIx2));
                }
                start += dz*8;
                break;
            default:
                if (entity_description[eptr->type].flags&HAS_LENGTH) {
                    dz = ((DRIFT*)eptr->p_elem)->length;
                    fprintf(fpm, "%s %e 0\n", eptr->name, start+=dz);
                    }
                break;
            }
        eptr = eptr->succ;
        }
    log_exit("output_magnets");
    fclose(fpm);
    }

