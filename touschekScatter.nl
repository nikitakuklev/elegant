/* file: touschekScatter.nl
 * purpose: namelist for simulating Touschek scatter effects
 * 
 * Aimin Xiao, 2007
 */
#include "namelist.h"

#namelist touschek_scatter static
        long seed = 123456789;
        long nbins = 100;
        double charge = 0;
        double frequency = 1;
        double delta = 0;
        double p_central_mev = 0.0;
        double emittance[2]  = {0, 0};
        double sigma_dp = 0.0;
        double sigma_s = 0.0;
        double distribution_cutoff[3] = {3, 3, 3};
        long save_initial_coordinates = 0;
        STRING bunch = "%s.bun";
        STRING loss = "%s.los";
        STRING distribution = "%s.dis";
        long n_particles_per_bunch = 1000;
        long total_scattered_particles = 10000;
        double weight_limit = 1.0;
        long verbosity = 0;
#end
/*
charge in nC.
frequency in Hz
delta in % "energy aperture of interested"
emittance in nm
sigma_dp in real, sigma_s in mm.
*/
