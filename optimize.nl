/* file: optimize.nl
 * purpose: namelists for optimizing beamline
 * 
 * Michael Borland, 1991
 */
#include "namelist.h"

static char *optimize_mode[N_OPTIM_MODES] = {
    "minimize", "maximize"
    } ;

static char *optimize_method[N_OPTIM_METHODS] = {
    "simplex", "grid", "sample", "powell", "randomsample", "randomwalk",
    } ;

#namelist optimization_term static
    STRING term = NULL;
    double weight = 1.0;
#end

#namelist optimization_setup static
    STRING equation = NULL;
    STRING mode = "minimize";
    STRING method = "simplex";
    double tolerance = -0.01;
    double target = -DBL_MAX;
    long soft_failure = 1;
    long n_passes = 2;
    long n_evaluations = 500;
    long n_restarts = 0;
    double restart_worst_term_factor = 1;
    long matrix_order = 1;
    STRING log_file = NULL;
    long verbose = 1;
    long output_sparsing_factor = 1;
    long balance_terms = 0;
    double simplex_divisor = 3;
    double simplex_pass_range_factor = 1;
#end

#namelist optimization_variable static
    STRING name = NULL;
    STRING item = NULL;
    double lower_limit = 0;
    double upper_limit = 0;
    double step_size = 1;
#end

#namelist optimization_constraint static
    STRING quantity = NULL;
    double lower = 0;
    double upper = 0;
#end

#namelist optimize static
     long summarize_setup = 0;
#end
