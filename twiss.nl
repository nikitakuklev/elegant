/* file: twiss.nl
 * contents: namelist for twiss_output
 * 
 * Michael Borland, 1989
 */
#include "namelist.h"

#namelist tune_shift_with_amplitude,struct
    long turns = 1000;
    double x0 = 1e-10;
    double y0 = 1e-10;
    double x1 = 1e-6;
    double y1 = 1e-6;
    long use_concatenation = 1;
#end

#namelist twiss_output
    STRING filename = NULL;
    long matched = 1;
    long output_at_each_step = 0;
    long output_before_tune_correction = 0;
    long final_values_only = 0;
    long statistics = 0;
    long radiation_integrals = 0;
    double beta_x = 1;
    double alpha_x = 0;
    double eta_x = 0;
    double etap_x = 0;
    double beta_y = 1;
    double alpha_y = 0;
    double eta_y = 0;
    double etap_y = 0;
    STRING reference_file = NULL;
    STRING reference_element = NULL;
    long reference_element_occurrence = 0;
    long concat_order = 3;
    long higher_order_chromaticity = 0;
#end

