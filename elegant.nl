/* file: track_nl.nl
 * contents: namelist and other stuff for track main module
 * 
 * Michael Borland, 1989
 */
#include "namelist.h"

#namelist global_settings static
     long inhibit_fsync = 0;
     long allow_overwrite = 1;
     long echo_namelists = 1;
     long mpi_randomization_mode = 3;
     long exact_normalized_emittance = 0;
     double SR_gaussian_limit = 3.0;
     long inhibit_seed_permutation = 0;
     STRING log_file = NULL;
     STRING error_log_file = NULL;
     long share_tracking_based_matrices = 1;
     long tracking_based_matrices_store_limit = 5000;
     long parallel_tracking_based_matrices = 1;
     long mpi_io_force_file_sync = 0;
     long mpi_io_read_buffer_size = 0;
     long mpi_io_write_buffer_size = 0;
     long usleep_mpi_io_kludge = 0;
     double tracking_matrix_step_factor = 1;
     double tracking_matrix_points = 9;
     double tracking_matrix_step_size[6] = {5e-5, 5e-5, 5e-5, 5e-5, 5e-5, 5e-5};
     long warning_limit = 10;
     short malign_method = 0;
     double slope_limit = SLOPE_LIMIT;
     double coord_limit = COORD_LIMIT;
#end

#namelist change_start static,struct
          STRING element_name = NULL;
          long ring_mode = 0;
          long element_occurence = 1;
          long delta_position = 0;
#end

#namelist change_end static,struct
          STRING element_name = NULL;
          long element_occurence = 1;
          long delta_position = 0;
#end

#namelist run_setup static
    STRING lattice = NULL;
    STRING use_beamline = NULL;
    STRING rootname = NULL;
    STRING output = NULL;
    STRING centroid = NULL;
    STRING bpm_centroid = NULL;
    STRING sigma = NULL;
    STRING final = NULL;
    STRING acceptance = NULL;
    STRING losses = NULL;
    long losses_include_global_coordinates = 0;
    double losses_s_limit[2] = {-DBL_MAX, DBL_MAX};
    STRING magnets = NULL;
    STRING semaphore_file = NULL;
    STRING parameters = NULL;
    long suppress_parameter_defaults = 0;
    STRING rfc_reference_output = NULL;
    long combine_bunch_statistics = 0;
    long wrap_around = 1;
    long final_pass = 0;
    long default_order = 2;
    long concat_order = 0;
    long print_statistics = 0;
    long show_element_timing = 0;
    long monitor_memory_usage = 0;
    long random_number_seed = 987654321;
    long correction_iterations = 1;
    long echo_lattice = 0;
    double p_central = 0.0;
    double p_central_mev = 0.0;
    long always_change_p0 = 0;
    long load_balancing_on = 0;
    long random_sequence_No = 1;
    STRING expand_for = NULL;
    long tracking_updates = 1;
    STRING search_path = NULL;
    long element_divisions = 0;
    long back_tracking = 0;
#end

#namelist change_particle,struct
    STRING name = "electron";
    double mass_ratio = 0;
    double charge_ratio = 0;
#end

#namelist track static
    long center_on_orbit = 0;
    long center_momentum_also = 1;
    long offset_by_orbit = 0;
    long offset_momentum_also = 1;
    long soft_failure = 1;
    long use_linear_chromatic_matrix = 0;
    long longitudinal_ring_only = 0;
    long ibs_only = 0;
    long stop_tracking_particle_limit = -1;
    long check_beam_structure = 0;
    STRING interrupt_file = "%s.interrupt";
#end

#namelist print_dictionary static
    STRING filename = NULL;
    long latex_form = 0;
    long SDDS_form = 0;
#end

#namelist semaphores
        STRING started = "%s.started";
        STRING done = "%s.done";
        STRING failed = "%s.failed";
#end

#namelist include_commands,struct
	  STRING filename = NULL;
	  long disable = 0;
#end
