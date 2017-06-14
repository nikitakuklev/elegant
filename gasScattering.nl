/* file: gasScattering.nl
 * contents: namelist for gas scattering simulation
 * 
 * Michael Borland, 2017
 */
#include "namelist.h"

#namelist gas_scattering static
    STRING output = NULL;
    STRING log_file = NULL;
    double xpmin = -0.01;
    double xpmax = 0.01;
    long nx = 11;
    double ypmin = 0.0;
    double ypmax = 0.01;
    long ny = 11;
    long twiss_scaling = 0;
    double s_start = 0;
    double s_end = DBL_MAX;
    STRING include_name_pattern = NULL;
    STRING include_type_pattern = NULL;
    long verbosity = 1;        
    long soft_failure = 0;
    long allow_watch_file_output = 0;
#end


