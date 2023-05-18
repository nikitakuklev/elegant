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
        long insert_before = 0;
        long add_at_start = 0;
        long add_at_end = 0;
        double s_start = -1;
        double s_end = -1;
        STRING element_def = NULL;
        long verbose = 0;
        long total_occurrences = 0;
        long occurrence[100]={0};
        long allow_no_matches = 0;
        long allow_no_insertions = 0;
#end


