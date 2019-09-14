/* file: ionEffects.nl
 * purpose: setup command for ion simulation
 * 
 * M. Borland, J. Calvey, 2017
 */
#include "namelist.h"

#namelist ion_effects static
          STRING pressure_profile = NULL;
          STRING ion_properties = NULL;
          STRING beam_output = NULL;
          long beam_output_all_locations = 0;
          STRING ion_density_output = NULL;
          long ion_output_all_locations = 1;
          long ion_species_output = 0;
          STRING field_calculation_method = "gaussian";
	  double distribution_fit_target = 0.1;
	  double distribution_fit_tolerance = 1e-3;
	  long distribution_fit_evaluations = 100;
	  long distribution_fit_passes = 1;
	  long distribution_fit_restarts = 1;
	  STRING fit_residual_type = NULL;
          long macro_ions = 0;
          long symmetrize = 0;
          long generation_interval = 1;
          long multiple_ionization_interval = 100;
          double ion_span[2] = {0, 0};
	  double ion_bin_divisor[2] = {10.0, 10.0};
	  double ion_range_multiplier[2] = {2.0, 2.0};
	  double ion_sigma_limit_multiplier[2] = {0, 0};
	  STRING ion_histogram_output = NULL;
	  double ion_histogram_output_s_start = 0;
	  double ion_histogram_output_s_end = 10;
	  long ion_histogram_output_interval = 1000;
	  long ion_histogram_min_output_bins = 200;
	  long ion_histogram_max_bins = 10000;
          long verbosity = 0;
#end

