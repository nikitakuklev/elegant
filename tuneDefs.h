#include "matlib.h"

/* structure for tune correction information */
typedef struct {
    double tunex, tuney;    /* desired tunes */
    char **name;            /* names of quadrupole families */
    long n_families;        /* number of families */
    long n_iterations;      /* number of times to repeat correction */
    double gain;            /* gain for correction */
    MATRIX *T;              /* Nfx2 matrix to give quadrupole strength changes to change 
                               chromaticities by given amount */
    MATRIX *dK1;           /* Nfx1 matrix of quadrupole strength changes */
    MATRIX *dtune;         /* 2x1 matrix of desired tune changes */
    } TUNE_CORRECTION;


/* prototypes for tune.c */
void setup_tune_correction(NAMELIST_TEXT *nltext, RUN *run, LINE_LIST *beamline, TUNE_CORRECTION *tune);
void do_tune_correction(TUNE_CORRECTION *tune, RUN *run, LINE_LIST *beamline, double *clorb, long step, long last_iteration);

