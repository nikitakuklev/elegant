/* file: tune.nl
 * contents: namelist for tune correction
 * 
 * Michael Borland, 1989
 */
#include "namelist.h"

#namelist correct_tunes
    STRING quadrupoles = NULL;
    double tune_x = -1;
    double tune_y = -1;
    long n_iterations = 1;
    STRING strength_log = NULL;
    long change_defined_values = 0;
#end


