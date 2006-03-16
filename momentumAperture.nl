/* file: momentumAperture.nl
 * contents: namelist for momentum aperture vs s
 * 
 * Michael Borland, 2006
 */
#include "namelist.h"

#namelist momentum_aperture static
    STRING output = NULL;
    double x_initial = 0;
    double y_initial = 0;
    double delta_negative_start = 0.0;
    double delta_negative_limit = -0.03;
    double delta_positive_start = 0.0;
    double delta_positive_limit = 0.03;
    long delta_points = 21;
    long n_splits = 3;
    long two_steps_back = 1;
    double s_start = 0;
    double s_end = DBL_MAX;
    STRING include_name_pattern = NULL;
    long verbosity = 0;        
#end


