#ifndef _GPU_FUNCS_H_
#define _GPU_FUNCS_H_
/* collection of gpu functions from various files common to multiple elements
 */

/* from tilt_matrices.c */
void gpu_rotateBeamCoordinates(long np, double angle);

/* from malign_mat.c */
void gpu_offsetBeamCoordinates(long np, double dx, double dy, double dz);

#ifdef __cplusplus
extern "C" {
#endif

/* from do_tracking.c */
void gpu_center_beam(CENTER *center, long np, long iPass, double p0);

/* from csbend.c (template conflict with other calls in gpu_csbend.cu) */
void gpu_addCorrectorRadiationKick(long np, ELEMENT_LIST *elem, long type,
           double Po, double *sigmaDelta2, long disableISR);

/* from do_tracking.c */
void gpu_offset_beam(long nToTrack, MALIGN *offset, double P_central);

void gpu_do_match_energy(long np, double* P_central, long change_beam);

void gpu_matr_element_tracking(VMATRIX *M, MATR *matr,
                               long np, double z);

void gpu_ematrix_element_tracking(VMATRIX *M, EMATRIX *matr,
                                  long np, double z, double *P_central);

void gpu_collect_trajectory_data(double* centroid, long np);

void gpu_set_central_momentum(long np, double  P_new, double *P_central);

#ifdef __cplusplus
}
#endif

#endif /* _GPU_FUNCS_H_ */
