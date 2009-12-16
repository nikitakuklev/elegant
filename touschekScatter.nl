/* file: touschekScatter.nl
 * purpose: namelist for simulating Touschek scatter effects
 * 
 * Aimin Xiao, 2007
 */
#include "namelist.h"

#namelist touschek_scatter static
        long nbins = 100;
        double charge = 0;
        double frequency = 1;
        double emit_x = 0;
        double emit_nx = 0;
        double emit_y = 0;
        double emit_ny = 0;
        double sigma_dp = 0;
        double sigma_s = 0;
        double distribution_cutoff[3] = {3, 3, 3};
        double delta = 0;
        double delta_mev = 0;
        STRING bunch = NULL;
        STRING loss = NULL;
        STRING distribution = NULL;
        STRING initial = NULL;
        STRING output = NULL;
        STRING longDist = NULL;
        STRING tranDist = NULL;
        long n_simulated = 5E6;
        double ignored_portion = 0.05;
        long i_start = 0;
        long i_end = 1;
        long do_track = 0;
	     long estimate = 0;
        long verbosity = 0;
#end
/*
charge in C.
frequency in Hz
delta in % "energy aperture of interested"
emittance in m.rad
sigma_dp is the realative energy spread, sigma_s in m.
*/
