/* prototypes for zibs.c */
/* ZAP-based intra-beam scattering routines by L. Emery */

double coulombLog (double gamma, double emitx, double emity,
                   double betaxAve, double betayAve, double sigz, double particles, 
                   long noWarning);

void IBSGrowthRates (double gamma, double emitx, double emity,
                     double sigmaDelta, double sigmaz,
                     double particles,
                     double emitx0, double sigmaDelta0,
                     double transSRdampRate, double longSRdampRate,
                     double coupling,
                     double *s, double *betax, double *alphax, double *betay, 
                     double *alphay, double *etax, double *etaxp, long elements, 
                     long superperiods, long verbosity,
                     double *xGrowthRate, double *yGrowthRate, double *zGrowthRate);

