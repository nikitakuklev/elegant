/* file: analyze.nl
 * purpose: namelist definition for transport analysis
 *
 * M.Borland, 1992 
 */
#include "namelist.h"

#namelist analyze_map static
    STRING output = NULL;
    STRING printout = NULL;
    double delta_x = 5e-5;
    double delta_xp = 5e-5;
    double delta_y = 5e-5;
    double delta_yp = 5e-5;
    double delta_s  = 5e-5;
    double delta_dp = 5e-5;
    long center_on_orbit = 0;
    long verbosity = 0;
    long n_points = 2;  /* backward compatibility. ignored */
    long printout_order = 2;
#end

