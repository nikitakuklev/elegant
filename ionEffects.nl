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
          long grid_points_x = 0;
          long grid_points_y = 0;
          STRING field_calculation_method = NULL;
          long macro_ions = 0;
          long generation_interval = 1;
          long multiple_ionization_interval = 100;
          double x_span = 0;
          double y_span = 0;
          long verbosity = 0;
#end

