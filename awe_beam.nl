/* file: awe_beam.nl
 * contents: namelist awe beam input
 * 
 * Michael Borland, 1989
 */
#include "namelist.h"

#namelist awe_beam
    STRING input = NULL;
    STRING input_type = "elegant";
    STRING dump_selection_string = NULL;
    long n_particles_per_ring = 1;
    long one_random_bunch = 0;
    long reuse_bunch = 0;
    long prebunched = 0;
    long sample_interval = 1;
    long n_dumps_to_skip = 0;
    long center_transversely = 0;
    long center_arrival_time = 0;
    double sample_fraction = 1;
    double p_lower = 0.0;
    double p_upper = 0.0;
#end

