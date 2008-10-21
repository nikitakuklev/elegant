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
        long skip = 1;
        long disable = 0;
        long add_at_end = 0;
        STRING element_def = NULL;
        long verbose = 1;
        long total_occurence = 0;
        long occurence[10]={0,0,0,0,0,0,0,0,0,0};
#end


