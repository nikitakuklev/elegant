/* Copyright 1994 by Michael Borland and Argonne National Laboratory,
 * all rights reserved.
 */
/* routine: tilt_matrices()
 * purpose: rotate matrices through a given angle 
 *
 * Michael Borland, 1989
 */
#include "track.h"
#include "mdb.h"

void tilt_matrices0(VMATRIX *M, double tilt)
{
    VMATRIX *rot, *irot, *Mr;

    log_entry("tilt_matrices0");

    rot  = rotation_matrix(tilt);
    irot = rotation_matrix(-tilt);
    Mr = tmalloc(sizeof(*Mr));
    initialize_matrices(Mr, M->order);

    concat_matrices(Mr,     M, rot, 0);
    concat_matrices(M  , irot, Mr , 0);
    free_matrices(rot); tfree(rot); rot = NULL;
    free_matrices(irot); tfree(irot); irot = NULL;
    free_matrices(Mr); tfree(Mr); Mr = NULL;
    log_exit("tilt_matrices0");
    }

void tilt_matrices(VMATRIX *M, double tilt)
{
    static VMATRIX Rot, IRot;
    static long initialized=0;
    VMATRIX Mr;
    double sin_tilt, cos_tilt;


    log_entry("tilt_matrices");
    
    if (tilt==0) {
        log_exit("tilt_matrices");
        return;
        }

    if (!initialized) {
        initialized = 1;
        initialize_matrices(&Rot, 1);
        initialize_matrices(&IRot, 1);
        }

    initialize_matrices(&Mr, M->order);
    
    if (fabs(tilt-PI/2)<1e-12) {
        sin_tilt = 1;
        cos_tilt = 0;
        }
    else if (fabs(tilt+PI/2)<1e-12) {
        sin_tilt = -1;
        cos_tilt = 0;
        }
    else if (fabs(tilt-PI)<1e-12 || fabs(tilt+PI)<1e-12) {
        sin_tilt = 0;
        cos_tilt = -1;
        }
    else {
        sin_tilt = sin(tilt);
        cos_tilt = cos(tilt);
        }

    /* Rotation matrix for (x, y) */
    IRot.R[0][0] =   Rot.R[0][0] =  cos_tilt ;
    IRot.R[0][2] = -(Rot.R[0][2] =  sin_tilt);
    IRot.R[2][0] = -(Rot.R[2][0] = -sin_tilt);
    IRot.R[2][2] =   Rot.R[2][2] =  cos_tilt ;

    /* Rotation matrix for (x', y') */
    IRot.R[1][1] =   Rot.R[1][1] =  cos_tilt ;
    IRot.R[1][3] = -(Rot.R[1][3] =  sin_tilt);
    IRot.R[3][1] = -(Rot.R[3][1] = -sin_tilt);
    IRot.R[3][3] =   Rot.R[3][3] =  cos_tilt ;

    IRot.R[4][4] = IRot.R[5][5] = Rot.R[4][4]  = Rot.R[5][5]  = 1;

    concat_matrices(&Mr,     M, &Rot, 0);
    concat_matrices(M  , &IRot, &Mr , 0);

    free_matrices(&Mr); 

    log_exit("tilt_matrices");
    }

/* routine: rotation_matrix()
 * purpose: return R matrix for rotation of beamline about given angle
 *
 * Michael Borland, 1989
 */

VMATRIX *rotation_matrix(double tilt)
{
    VMATRIX *Rot;
    double sin_tilt, cos_tilt;

    log_entry("rotation_matrix");

    if (fabs(tilt-PI/2)<1e-12) {
        sin_tilt = 1;
        cos_tilt = 0;
        }
    else if (fabs(tilt-PI)<1e-12) {
        sin_tilt = 0;
        cos_tilt = -1;
        }
    else {
        sin_tilt = sin(tilt);
        cos_tilt = cos(tilt);
        }

    Rot = tmalloc(sizeof(*Rot));
    initialize_matrices(Rot, 1);

    /* Rotation matrix for (x, y) */
    Rot->R[0][0] =  cos_tilt ;
    Rot->R[0][2] =  sin_tilt ;
    Rot->R[2][0] = -sin_tilt ;
    Rot->R[2][2] =  cos_tilt ;

    /* Rotation matrix for (x', y') */
    Rot->R[1][1] =  cos_tilt ;
    Rot->R[1][3] =  sin_tilt ;
    Rot->R[3][1] = -sin_tilt ;
    Rot->R[3][3] =  cos_tilt ;

    Rot->R[4][4] = Rot->R[5][5] = 1;

    log_exit("rotation_matrix");
    return(Rot);
    }

void rotate_coordinates(double *coord, double angle)
{
    static double x, xp, y, yp;
    static double sin_a, cos_a;

    if (!angle)
        return;
    sin_a = sin(angle);
    cos_a = cos(angle);
    x = coord[0]; xp = coord[1]; y = coord[2]; yp = coord[3];
    coord[0] =   x*cos_a + y*sin_a;
    coord[2] =  -x*sin_a + y*cos_a;
    coord[1] =  xp*cos_a + yp*sin_a;
    coord[3] = -xp*sin_a + yp*cos_a;
    }

