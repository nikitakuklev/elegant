/* file: alter.nl
 * purpose: namelist for altering element properties
 * 
 * Michael Borland, 1999
 */
#include "namelist.h"

#namelist alter_element static
        STRING name = NULL;
        STRING item = NULL;
        STRING exclude = NULL;
        double value = 0;
        long differential = 0;
        long multiplicative = 0;
        long verbose = 0;
        long allow_missing_parameters = 0;
#end


