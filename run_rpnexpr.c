/* Copyright 1994 by Michael Borland and Argonne National Laboratory,
 * all rights reserved.
 */
/* file: run_rpn.c
 *
 * Michael Borland, 1991
 */
#include "mdb.h"
#include "track.h"
#include "run_rpnexpr.h"

void run_rpn_expression(NAMELIST_TEXT *nltext)
{
/*    expression = NULL; 
 */

    log_entry("run_rpn_expression");
   
    /* process the namelist text */
    set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
    set_print_namelist_flags(0);
    process_namelist(&rpn_expression, nltext);
    print_namelist(stdout, &rpn_expression);

    if (expression) {
        rpn(expression);
        if (rpn_check_error()) exit(1);
        }
    log_exit("run_rpn_expression");
    }
