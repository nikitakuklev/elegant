/* file: insert_elements.nl
 * purpose: namelist for insert elements into beamline
 * 
 * A.Xiao, 2007
 */
#include "namelist.h"

#namelist insert_elements static
        STRING name = NULL;
        STRING type = NULL;
        STRING exclude = NULL;
        long skip = 0;
        long disable = 0;
        long add_at_end = 0;
        STRING element_def = NULL;
        long verbose = 1;
#end


