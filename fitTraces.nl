/* file: fitTraces.nl
 * purpose: namelists for fitting multiple traces through a beamline
 * 
 * Michael Borland, 1997
 */
/*
 * $Log: not supported by cvs2svn $
 */

#include "namelist.h"

#namelist fit_traces static
        STRING trace_data_file = NULL;
        STRING fit_parameters_file = NULL;
        long iterations = 10;
        long sub_iterations = 10;
        double convergenceFactor = 1;
        double convergenceFactorDelta = 0.005;
        double convergenceFactorBackoff = 0.5;
        double convergenceFactorMin = 1.0;
        double convergenceFactorMax = 2.0;
        STRING output = NULL;
        double target = 1e-12;
#end

