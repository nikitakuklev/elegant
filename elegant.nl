/* file: track_nl.nl
 * contents: namelist and other stuff for track main module
 * 
 * Michael Borland, 1989
 */
#include "namelist.h"

#namelist run_setup static
    STRING lattice = NULL;
    STRING use_beamline = NULL;
    STRING rootname = NULL;
    STRING output = NULL;
    STRING centroid = NULL;
    STRING sigma = NULL;
    STRING final = NULL;
    STRING acceptance = NULL;
    STRING losses = NULL;
    STRING magnets = NULL;
    STRING semaphore_file = NULL;
    STRING parameters = NULL;
    long combine_bunch_statistics = 0;
    long wrap_around = 1;
    long default_order = 2;
    long concat_order = 0;
    long print_statistics = 0;
    long random_number_seed = 987654321;
    long correction_iterations = 1;
    long echo_lattice = 0;
    double p_central = 0.0;
    long always_change_p0 = 0;
    STRING expand_for = NULL;
    long tracking_updates = 1;
#end

#namelist track static
    long center_on_orbit = 0;
    long center_momentum_also = 1;
    long soft_failure = 1;
    long use_linear_chromatic_matrix = 0;
    long longitudinal_ring_only = 0;
#end

#namelist print_dictionary static
    STRING filename = NULL;
    long latex_form = 0;
#end
