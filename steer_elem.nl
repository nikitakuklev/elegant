/* file: steer_elem.nl 
 * purpose: namelist definition for adding correction elements
 *
 * Michael Borland, 1991
 */
#include "namelist.h"

#namelist steering_element static
    STRING name = NULL;
    STRING item = NULL;
    STRING plane = "h";
    double tweek = 1e-3;
    double limit = 0;
#end

