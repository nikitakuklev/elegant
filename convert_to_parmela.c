/* Copyright 1994 by Michael Borland and Argonne National Laboratory,
 * all rights reserved.
 */
/* routine: convert_to_parmela
 * purpose: convert parsed mad lattice to parmela format (writes output
 *          to a file)
 * 
 * Michael Borland, 1990.
 */
#include "mdb.h"
#include "madto.h"

#define PC_CONST (33.3563603408)
#define RAD_TO_DEG(x) ( (int) ((x)*180./PI) )

void convert_to_parmela(char *outputfile, LINE_LIST *beamline, long flip_k, double angle_tolerance, 
                        char *header_file, char *ender_file,
                        double quad_aperture, double sext_aperture, double pc)
{
    FILE *fp;
    ELEMENT_LIST *eptr;
    double rho, K_cen;
    QUAD *quad;
    SEXT *sext;
    BEND *bend;
    HCOR *hcor;
    DRIFT *drift;
    char s[100];
    int quad_sign;
    FILE *fpi;

    fp = fopen_e(outputfile, "w", 0);

    quad_sign = flip_k?-1:1;
    
    if (header_file) {
        fpi = fopen_e(header_file, "r", 1);
        while (fgets(s, 256, fpi))
            fputs(s, fp);
        fclose(fpi);
        }

    fprintf(fp, "run 0 2856 0 %le\n", K_cen=0.511*(sqrt(sqr(pc/0.511)+1)-1) );
    eptr = &(beamline->elem);
    pc *= PC_CONST;
    while (eptr) {
        switch (eptr->type) {
          case T_QUAD:
            quad = (QUAD*)eptr->p_elem;
            if (quad->tilt!=0)
                fprintf(fp, "rotate 0 %le 0 %le\n",
                        quad_aperture, quad->tilt*180/PI);
            fprintf(fp, "quad %le %le 0 %le\n",
                    quad->length*100, quad_aperture*100, 
                    quad->k1*pc*quad_sign/100);
            if (quad->tilt!=0)
                fprintf(fp, "rotate 0 %le 0 %le\n",
                        quad_aperture, -quad->tilt*180/PI);
            break;
          case T_SBEN: case T_RBEN:
            bend  = (BEND*)eptr->p_elem;
            if (eptr->type==T_RBEN) {
                rho = bend->length/(2*sin(bend->angle/2));
                bend->length = rho*bend->angle;
                }
            else {
                rho = bend->length/bend->angle;
                }
            if (bend->tilt!=0)
                fprintf(fp, "rotate 1e6 0 %le\n", bend->tilt*180/PI);
            fprintf(fp, "bend %le 1e6 0 %le %le %le %le\n",
                    100*bend->length, K_cen, bend->angle*180/PI,
                    bend->e1*180/PI, bend->e2*180/PI);
            if (bend->tilt!=0)
                fprintf(fp, "rotate 1e6 0 %le\n", -bend->tilt*180/PI);
            break;
          case T_DRIF:
            drift = (DRIFT*)eptr->p_elem;
            fprintf(fp, "drift %le 1e6 0\n", 100*drift->length);
            break;
          case T_VCOR:
          case T_HCOR:
            printf("HCOR/VCOR not supported for PARMELA--drift used\n");
            hcor = (HCOR*)eptr->p_elem;
            fprintf(fp, "drift %le 1e6 0\n", 100*hcor->length);
            break;
          case T_SEXT:
            sext = (SEXT*)eptr->p_elem;
            printf("SEXTUPOLE not supported for PARMELA--drift assumed\n");
            fprintf(fp, "drift %le 1e6 0\n", 100*sext->length);
            break;
          default:
            printf("unknown element type--element '%s' ignored\n",
                   eptr->name);
            break;
            }
        eptr = eptr->succ;
        }
        
    if (ender_file) {
        fpi = fopen_e(ender_file, "r", 1);
        while (fgets(s, 256, fpi))
            fputs(s, fp);
        fclose(fpi);
        }

    }

