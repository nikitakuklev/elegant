/* Copyright 1994 by Michael Borland and Argonne National Laboratory,
 * all rights reserved.
 */
/* file: floor.c
 * purpose: computation/output of floor coordinates
 *
 * Michael Borland, 1993
 */
#include "mdb.h"
#include "track.h"

static SDDS_TABLE SDDS_floor;

#define IC_S 0
#define IC_X 1
#define IC_Z 2
#define IC_THETA 3
#define IC_ELEMENT 4
#define IC_OCCURENCE 5
#define IC_TYPE 6
#define N_COLUMNS 7
static SDDS_DEFINITION column_definition[N_COLUMNS] = {
    {"s", "&column name=s, units=m, type=double, description=\"Distance\" &end"},
    {"X", "&column name=X, units=m, type=double, description=\"Transverse survey coordinate\" &end"},
    {"Z", "&column name=Z, units=m, type=double, description=\"Longitudinal survey coordinate\" &end"},
    {"theta", "&column name=theta, symbol=\"$gQ$r\", units=radians, type=double, description=\"Survey angle\" &end"},
    {"ElementName", "&column name=ElementName, type=string, description=\"Element name\", format_string=%10s &end"},
    {"ElementOccurence", 
         "&column name=ElementOccurence, type=long, description=\"Occurence of element\", format_string=%6ld &end"},
    {"ElementType", "&column name=ElementType, type=string, description=\"Element-type name\", format_string=%10s &end"},
    } ;

void output_floor_coordinates(NAMELIST_TEXT *nltext, RUN *run, LINE_LIST *beamline)
{
    double *data, length;
    ELEMENT_LIST *elem, *last_elem;
    double dX, dZ, dtheta, rho, angle;
    long is_bend, is_magnet, n_points, row_index;
    BEND *bend; KSBEND *ksbend; CSBEND *csbend; MALIGN *malign; CSRCSBEND *csrbend;
    char label[200];
#include "floor.h"

    log_entry("output_floor_coordinates");

    /* process namelist input */
    set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
    set_print_namelist_flags(0);
    process_namelist(&floor_coordinates, nltext);
    print_namelist(stdout, &floor_coordinates);

    if (magnet_centers && vertices_only)
        bomb("you can simultaneously request magnet centers and vertices only output", NULL);
    if (filename)
        filename = compose_filename(filename, run->rootname);
    else
        bomb("filename must be given for floor coordinates", NULL);
    
    SDDS_ElegantOutputSetup(&SDDS_floor, filename, SDDS_BINARY, 1, "floor coordinates", 
                            run->runfile, run->lattice, NULL, 0,
                            column_definition, N_COLUMNS, "floor coordinates", 
                            SDDS_EOS_NEWFILE|SDDS_EOS_COMPLETE);

    n_points = beamline->n_elems+1;
    if (vertices_only)
        n_points = 2;
    last_elem = NULL;
    if (include_vertices || vertices_only) {
        elem = &(beamline->elem);
        while (elem) {
            switch (elem->type) {
              case T_RBEN: 
              case T_SBEN: 
              case T_KSBEND: 
              case T_CSBEND:
              case T_CSRCSBEND:
                n_points ++;
                break;
              default:
                break;
                }
            last_elem = elem;
            elem = elem->succ;
            }
        }

    if (!SDDS_StartTable(&SDDS_floor, 2*n_points)) {
        SDDS_SetError("Unable to start SDDS table (output_floor_coordinates)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }

    row_index = 0;
    if (!SDDS_SetRowValues(&SDDS_floor, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, row_index++,
                           IC_S, (double)0.0, IC_X, X0, IC_Z, Z0, IC_THETA, theta0,
                           IC_ELEMENT, "_BEG_", IC_OCCURENCE, (long)1, IC_TYPE, "MARK", -1)) {
        SDDS_SetError("Unable to set SDDS row (output_floor_coordinates)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
        
    elem = &(beamline->elem);
    data = tmalloc(sizeof(*data)*4);
    while (elem) {
        dX = dZ = dtheta = length = 0;
        is_bend = is_magnet = 0;
        switch (elem->type) {
          case T_RBEN: case T_SBEN:
            is_bend = 1;
            bend = (BEND*)elem->p_elem;
            dtheta = -(angle=bend->angle);
            rho = (length=bend->length)/bend->angle;
            break;
          case T_KSBEND:
            is_bend = 1;
            ksbend = (KSBEND*)elem->p_elem;
            dtheta = -(angle=ksbend->angle);
            rho = (length=ksbend->length)/ksbend->angle;
            break;
          case T_CSBEND:
            is_bend = 1;
            csbend = (CSBEND*)elem->p_elem;
            dtheta = -(angle=csbend->angle);
            rho = (length=csbend->length)/csbend->angle;
            break;
          case T_CSRCSBEND:
            is_bend = 1;
            csrbend = (CSRCSBEND*)elem->p_elem;
            dtheta = -(angle=csrbend->angle);
            rho = (length=csrbend->length)/csrbend->angle;
            break;
          case T_MALIGN:
            malign = (MALIGN*)elem->p_elem;
            dX = malign->dx;
            dZ = malign->dz;
            dtheta = atan(malign->dxp);
            break;
          default:
            if (entity_description[elem->type].flags&HAS_LENGTH)
                length= dZ = *((double*)(elem->p_elem));
            if (entity_description[elem->type].flags&IS_MAGNET)
                is_magnet = 1;
            break;
            }
        if (is_bend) {
            if (include_vertices || vertices_only || magnet_centers) {
                /* output data for the vertex point */
                if (angle && !isnan(rho))
                  data[0] += (dZ=rho*tan(angle/2));
                else
                  data[0] += (dZ=length/2);
                dX = 0;
                rotate_xy(&dX, &dZ, theta0);
                data[1] = X0+dX;
                data[2] = Z0+dZ;
                data[3] = theta0;
                sprintf(label, "%s-VP", elem->name);
                if (!SDDS_SetRowValues(&SDDS_floor, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, row_index++,
                                       IC_S, data[0], IC_X, data[1], IC_Z, data[2], IC_THETA, data[3],
                                       IC_ELEMENT, label, IC_OCCURENCE, elem->occurence, IC_TYPE, "VERTEX-POINT", -1)) {
                    SDDS_SetError("Unable to set SDDS row (output_floor_coordinates)");
                    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
                    }
                }
            /* calculate offsets for end of bend magnet */
            if (angle && !isnan(rho)) {
              dX = rho*(cos(angle)-1);
              dZ = rho*sin(angle);
            } else {
              dX = 0;
              dZ = length;
            }              
          }
        if (is_magnet && magnet_centers) {
            dX /= 2;
            dZ /= 2;
            sprintf(label, "%s-C", elem->name);
            }
        else
            strcpy(label, elem->name);

        rotate_xy(&dX, &dZ, theta0);
        theta0 += dtheta;
        while (theta0>PIx2)
            theta0 -= PIx2;
        while (theta0<-PIx2)
            theta0 += PIx2;
        X0 += dX;
        Z0 += dZ;

        if (!vertices_only || (!last_elem || elem==last_elem)) {
            if (!is_bend && magnet_centers)
                data[0] = elem->end_pos - length/2;
            else
                data[0] = elem->end_pos;
            data[1] = X0;
            data[2] = Z0;
            data[3] = theta0;
            if (!SDDS_SetRowValues(&SDDS_floor, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, row_index++,
                                   IC_S, data[0], IC_X, data[1], IC_Z, data[2], IC_THETA, data[3],
                                   IC_ELEMENT, label, IC_OCCURENCE, elem->occurence, IC_TYPE, 
                                   entity_name[elem->type], -1)) {
                SDDS_SetError("Unable to set SDDS row (output_floor_coordinates)");
                SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
                }
            }

        if (is_magnet && magnet_centers) {
            X0 += dX;
            Z0 += dZ; 
            }

        elem = elem->succ;
        }
    if (!SDDS_WriteTable(&SDDS_floor)) {
        SDDS_SetError("Unable to write floor coordinate data (output_floor_coordinates)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
    if (!SDDS_Terminate(&SDDS_floor)) {
        SDDS_SetError("Unable to terminate SDDS file (output_floor_coordinates)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
    log_exit("output_floor_coordinates");
    }


void final_floor_coordinates(LINE_LIST *beamline, double *X, double *Z, double *Theta)
{
  double *data, length;
  ELEMENT_LIST *elem, *last_elem;
  double dX, dZ, dtheta, rho, angle;
  double X0, Z0, theta0;
  long is_bend, is_magnet;
  BEND *bend; KSBEND *ksbend; CSBEND *csbend; MALIGN *malign;
  CSRCSBEND *csrbend;

  last_elem = NULL;
  elem = &(beamline->elem);
  X0 = Z0 = theta0 = 0;
  while (elem) {
    dX = dZ = dtheta = length = 0;
    is_bend = is_magnet = 0;
    switch (elem->type) {
    case T_RBEN: case T_SBEN:
      is_bend = 1;
      bend = (BEND*)elem->p_elem;
      dtheta = -(angle=bend->angle);
      rho = (length=bend->length)/bend->angle;
      break;
    case T_KSBEND:
      is_bend = 1;
      ksbend = (KSBEND*)elem->p_elem;
      dtheta = -(angle=ksbend->angle);
      rho = (length=ksbend->length)/ksbend->angle;
      break;
    case T_CSBEND:
      is_bend = 1;
      csbend = (CSBEND*)elem->p_elem;
      dtheta = -(angle=csbend->angle);
      rho = (length=csbend->length)/csbend->angle;
      break;
    case T_CSRCSBEND:
      is_bend = 1;
      csrbend = (CSRCSBEND*)elem->p_elem;
      dtheta = -(angle=csrbend->angle);
      rho = (length=csrbend->length)/csrbend->angle;
      break;
    case T_MALIGN:
      malign = (MALIGN*)elem->p_elem;
      dX = malign->dx;
      dZ = malign->dz;
      dtheta = atan(malign->dxp);
      break;
    default:
      if (entity_description[elem->type].flags&HAS_LENGTH)
        length= dZ = *((double*)(elem->p_elem));
      if (entity_description[elem->type].flags&IS_MAGNET)
        is_magnet = 1;
      break;
    }
    if (is_bend) {
      /* calculate offsets for end of bend magnet */
      if (angle && !isnan(rho)) {
        dX = rho*(cos(angle)-1);
        dZ = rho*sin(angle);
      } else {
        dX = 0;
        dZ = length;
      }
    }

    rotate_xy(&dX, &dZ, theta0);
    theta0 += dtheta;
    while (theta0>PIx2)
      theta0 -= PIx2;
    while (theta0<-PIx2)
      theta0 += PIx2;
    X0 += dX;
    Z0 += dZ;

    elem = elem->succ;
  }
  *X = X0;
  *Z = Z0;
  *Theta = theta0;
}

