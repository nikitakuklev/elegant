/* file: load_parameters.nl
 * purpose: namelist for loading parameters from external file
 * 
 * Michael Borland, 1993
 */
#include "namelist.h"

#namelist load_parameters
        STRING filename = NULL;
        long change_defined_values = 0;
        long clear_settings = 0;
        long allow_missing_elements = 0;
        long allow_missing_parameters = 0;
        long force_occurence_data = 0;
        long verbose = 0;
#end

