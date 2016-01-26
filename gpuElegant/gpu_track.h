#ifndef _GPUTRACK_H_
#define _GPUTRACK_H_

#ifdef  __cplusplus
#include <complex> // Needed to avoid conflicts w/ track
#endif

#ifndef COORDINATES_PER_PARTICLE
#include <mdb.h>
#include <track.h>

#undef pow3
#define pow3(x) x*x*x
#undef pow4
#define pow4(x) x*x*x*x
#undef pow5
#define pow5(x) x*x*x*x*x
#undef pow6
#define pow6(x) x*x*x*x*x*x

#ifdef max
#undef max
#endif

#ifdef min
#undef min
#endif

#ifdef m_alloc
#undef m_alloc
#endif

#endif /* #ifndef COORDINATES_PER_PARTICLE */

#endif
