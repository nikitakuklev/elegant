/* file: frequencyMap.nl
 * purpose: namelist definition for frequency map analysis
 *
 * M.Borland, 1992 
 */
#include "namelist.h"

#namelist frequency_map static
    STRING output = NULL;
    double xmin =  0.0;
    double xmax =  0.1;
    double ymin =  0.0;
    double ymax =  0.1;
    long nx = 21;
    long ny = 11;
    long verbosity = 1;
#end

