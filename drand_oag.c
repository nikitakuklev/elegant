
/* This routine is copied from the drand.c under the SDDS directory
   for parallelizing the random number generator. */

#include "mdb.h"
#include "track.h"
#include "f2c.h"

extern double dlaran_(integer *seed); 

double random_1_elegant(long iseed)
{

    static short initialized = 0;
    static integer seed[4] = {0,0,0,0};

    if (!initialized || iseed<0) {
        if (iseed<0)
          iseed = -iseed;
#if (!USE_MPI)
        random_2(-(iseed+2));
#else
        random_2(-(iseed+2*(myid+4)));
#endif
        random_3(-(iseed+4));
        random_4(-(iseed+6));
        seed[3] = ((iseed & 4095)/2)*2+1;
        seed[2] = (iseed >>= 12) & 4095;
        seed[1] = (iseed >>= 12) & 4095;
        seed[0] = (iseed >>= 12) & 4095;
        initialized = 1;
        }
    if (!initialized)
        bomb("random_1_elegant not properly initialized", NULL);

    return dlaran_(seed);

}
