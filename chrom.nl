/* file: chrom.nl
 * contents: namelist for chromaticity correction
 * 
 * Michael Borland, 1989
 */
#include "namelist.h"

#namelist chromaticity
    STRING sextupoles = NULL;
    double dnux_dp = 0;
    double dnuy_dp = 0;
    double sextupole_tweek = 1e-3;
    long n_iterations = 1;
    STRING strength_log = NULL;
    long change_defined_values = 0;
    double strength_limit = 0;
#end


