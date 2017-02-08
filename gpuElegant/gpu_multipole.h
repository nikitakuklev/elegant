#ifndef _GPU_MULTIPOLE_H_
#define _GPU_MULTIPOLE_H_

#include <gpu_track.h>
#include <gpu_base.h>

#ifdef __cplusplus
extern "C" {
#endif

long gpu_multipole_tracking2(long n_part, ELEMENT_LIST *elem,
       double p_error, double Po, double **accepted, double z_start,
       MAXAMP *maxamp, APERTURE_DATA *apFileData,
       double *sigmaDelta2);

#ifdef __cplusplus
}
#endif

#endif /* _GPU_MULTIPOLE_H_ */
