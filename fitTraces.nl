/* file: fitTraces.nl
 * purpose: namelists for fitting multiple traces through a beamline
 * 
 * Michael Borland, 1997
 */
/*
 * $Log: not supported by cvs2svn $
 * Revision 1.2  1997/10/20 14:57:11  borland
 * Improved trace fitting and related routines.  Added output of traces
 * after fitting.  Fixed some output-related bugs.
 *
 * Revision 1.1  1997/08/13 20:03:45  borland
 * First version.
 *
 */

#include "namelist.h"

#namelist fit_traces static
        STRING trace_data_file = NULL;
        STRING fit_parameters_file = NULL;
        long iterations = 10;
        long sub_iterations = 10;
        double convergence_factor = 1;
        double convergence_factor_divisor = 2;
        double convergence_factor_multiplier = 1.5;
        long convergence_increase_steps = 5;
        double convergence_factor_min = 1.0;
        double convergence_factor_max = 2.0;
        double position_change_limit = 0;
        double slope_change_limit = 0;
        long trace_sub_iterations = 500;
        double trace_convergence_factor = 0.9;
        double bpm_calibration_error_limit = 0.1;
        STRING fit_output_file = NULL;
        STRING trace_output_file = NULL;
        double target = 1e-12;
#end

