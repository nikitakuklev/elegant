/* Copyright 1994 by Michael Borland and Argonne National Laboratory,
 * all rights reserved.
 */
/* file: optimize.c
 *
 * Michael Borland, 1991
 */
#include "mdb.h"
#include "track.h"
#include "optimize.h"
#include "match_string.h"

char *gen_pcode(char *s);

#define DEBUG 0

void do_optimization_setup(OPTIMIZATION_DATA *optimization_data, NAMELIST_TEXT *nltext, RUN *run, LINE_LIST *beamline)
{
    long i;

    log_entry("do_optimization_setup");

    /* process the namelist text */
    set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
    set_print_namelist_flags(0);
    process_namelist(&optimization_setup, nltext);
    print_namelist(stderr, &optimization_setup);

    /* check validity of input values, and copy into structure */
    if ((optimization_data->equation=equation)==NULL || strlen(equation)==0)
        bomb("equation is invalid in optimize namelist", NULL);
    if ((optimization_data->mode=match_string(mode, optimize_mode, N_OPTIM_MODES, EXACT_MATCH))<0)
        bomb("unknown optimization mode", NULL);
    if ((optimization_data->method=match_string(method, optimize_method, N_OPTIM_METHODS, EXACT_MATCH))<0)
        bomb("unknown optimization method", NULL);
    if ((optimization_data->tolerance=tolerance)==0)
        bomb("tolerance == 0", NULL);
    if ((optimization_data->n_passes=n_passes)<=0)
        bomb("n_passes <= 0", NULL);
    if ((optimization_data->n_evaluations=n_evaluations)<=0)
        bomb("n_evaluations <= 0", NULL);
    optimization_data->soft_failure = soft_failure;
    if (log_file) {
        if (str_in(log_file, "%s"))
            log_file = compose_filename(log_file, run->rootname);
        if (strcmp(log_file, "/dev/tty")==0 || strcmp(log_file, "tt:")==0)
            optimization_data->fp_log = stderr;
        else if ((optimization_data->fp_log=fopen_e(log_file, "w", FOPEN_RETURN_ON_ERROR|FOPEN_SAVE_IF_EXISTS))==NULL)
            bomb("unable to open log file", NULL);
        }
    if (optimization_data->mode==OPTIM_MODE_MAXIMUM && target!=-DBL_MAX)
      target = -target;
    optimization_data->target = target;

    /* reset flags for elements that may have been varied previously */
    if (optimization_data->variables.n_variables)
        set_element_flags(beamline, optimization_data->variables.element, NULL, NULL, NULL, optimization_data->variables.n_variables,
                PARAMETERS_ARE_STATIC, 0, 1, 0);

    /* initialize other elements of the structure */
    optimization_data->new_data_read = 0;
    optimization_data->variables.n_variables = 0;
    optimization_data->covariables.n_covariables = 0;
    optimization_data->constraints.n_constraints = 0;
    log_exit("do_optimization_setup");
    }

void add_optimization_variable(OPTIMIZATION_DATA *optimization_data, NAMELIST_TEXT *nltext, RUN *run, LINE_LIST *beamline)
{
    long n_variables;
    OPTIM_VARIABLES *variables;
    ELEMENT_LIST *context;
    /* these are used to append a dummy name to the variables list for use with final parameters output: */
    static char *extra_name = "optimized";
    static char *extra_unit = "";

    log_entry("add_optimization_variable");

    if ((n_variables = optimization_data->variables.n_variables)==0) {
        if (optimization_data->new_data_read)
            bomb("improper sequencing of variation and tracking", NULL);
        optimization_data->new_data_read = 1;
        }

    variables = &(optimization_data->variables);
    variables->element = trealloc(variables->element, sizeof(*variables->element)*(n_variables+2));
    variables->item = trealloc(variables->item, sizeof(*variables->item)*(n_variables+2));
    variables->lower_limit = trealloc(variables->lower_limit, sizeof(*variables->lower_limit)*(n_variables+2));
    variables->upper_limit = trealloc(variables->upper_limit, sizeof(*variables->upper_limit)*(n_variables+2));
    variables->step = trealloc(variables->step, sizeof(*variables->step)*(n_variables+2));
    variables->orig_step = trealloc(variables->orig_step, sizeof(*variables->orig_step)*(n_variables+2));
    variables->varied_quan_name = trealloc(variables->varied_quan_name, sizeof(*variables->varied_quan_name)*(n_variables+2));
    variables->varied_quan_unit = trealloc(variables->varied_quan_unit, sizeof(*variables->varied_quan_unit)*(n_variables+2));
    variables->varied_type = trealloc(variables->varied_type, sizeof(*variables->varied_type)*(n_variables+2));
    variables->varied_param = trealloc(variables->varied_param, sizeof(*variables->varied_param)*(n_variables+2));
    variables->varied_quan_value = trealloc(variables->varied_quan_value, sizeof(*variables->varied_quan_value)*(n_variables+2));
    variables->initial_value = trealloc(variables->initial_value, sizeof(*variables->initial_value)*(n_variables+2));
    variables->memory_number = trealloc(variables->memory_number, sizeof(*variables->memory_number)*(n_variables+2));

    /* process namelist text */
    /* can't use automatic defaults, because of DBL_MAX being a nonconstant object */
    set_namelist_processing_flags(0);
    set_print_namelist_flags(0);
    name = item = NULL;
    step_size = 1;
    lower_limit = -(upper_limit = DBL_MAX);
    process_namelist(&optimization_variable, nltext);
    print_namelist(stderr, &optimization_variable);

    /* check for valid input */
    if (name==NULL)
        bomb("element name missing in optimization_variable namelist", NULL);
    str_toupper(name);
    context = NULL;
    if (!find_element(name, &context, &(beamline->elem))) {
        fprintf(stderr, "error: cannot vary element %s--not in beamline\n", name);
        exit(1);
        }
    cp_str(&variables->element[n_variables], name);
    variables->varied_type[n_variables] = context->type;
    if (item==NULL)
        bomb("item name missing in optimization_variable list", NULL);
    str_toupper(item);
    if ((variables->varied_param[n_variables] = confirm_parameter(item, context->type))<0) {
        fprintf(stderr, "error: cannot vary %s--no such parameter for %s\n",item, name);
        exit(1);
        }
    cp_str(&variables->item[n_variables], item);
    cp_str(&variables->varied_quan_unit[n_variables], 
        entity_description[context->type].parameter[variables->varied_param[n_variables]].unit);
    if (!get_parameter_value(variables->varied_quan_value+n_variables, name, variables->varied_param[n_variables],
            context->type, beamline))
        bomb("unable to get initial value for parameter", NULL);
    if (lower_limit>=upper_limit)
        bomb("lower_limit >= upper_limit", NULL);

    variables->initial_value[n_variables] = variables->varied_quan_value[n_variables];
    variables->orig_step[n_variables]     = step_size;
    variables->varied_quan_name[n_variables]  = tmalloc(sizeof(char)*(strlen(name)+strlen(item)+3));
    sprintf(variables->varied_quan_name[n_variables], "%s.%s", name, item);
    variables->lower_limit[n_variables] = lower_limit;
    variables->upper_limit[n_variables] = upper_limit;
    rpn_store(variables->initial_value[n_variables], 
              variables->memory_number[n_variables] = rpn_create_mem(variables->varied_quan_name[n_variables]));

    variables->varied_quan_name[n_variables+1] = extra_name;
    variables->varied_quan_unit[n_variables+1] = extra_unit;

    optimization_data->variables.n_variables += 1;
    log_exit("add_optimization_variable");
    }

void add_optimization_covariable(OPTIMIZATION_DATA *optimization_data, NAMELIST_TEXT *nltext, RUN *run, LINE_LIST *beamline)
{
#include "optim_covariable.h"
    long n_covariables;
    OPTIM_COVARIABLES *covariables;
    ELEMENT_LIST *context;

    log_entry("add_optimization_covariable");

    n_covariables = optimization_data->covariables.n_covariables;

    covariables = &(optimization_data->covariables);
    covariables->element = trealloc(covariables->element, sizeof(*covariables->element)*(n_covariables+1));
    covariables->item = trealloc(covariables->item, sizeof(*covariables->item)*(n_covariables+1));
    covariables->equation = trealloc(covariables->equation, sizeof(*covariables->equation)*(n_covariables+1));
    covariables->pcode = trealloc(covariables->pcode, sizeof(*covariables->pcode)*(n_covariables+1));
    covariables->varied_quan_name = trealloc(covariables->varied_quan_name, sizeof(*covariables->varied_quan_name)*(n_covariables+1));
    covariables->varied_quan_unit = trealloc(covariables->varied_quan_unit, sizeof(*covariables->varied_quan_unit)*(n_covariables+1));
    covariables->varied_type = trealloc(covariables->varied_type, sizeof(*covariables->varied_type)*(n_covariables+1));
    covariables->varied_param = trealloc(covariables->varied_param, sizeof(*covariables->varied_param)*(n_covariables+1));
    covariables->varied_quan_value = trealloc(covariables->varied_quan_value, sizeof(*covariables->varied_quan_value)*(n_covariables+1));
    covariables->memory_number = trealloc(covariables->memory_number, sizeof(*covariables->memory_number)*(n_covariables+1));


/*    name = item = equation = NULL; 
 */
    
    /* process namelist text */
    set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
    set_print_namelist_flags(0);
    process_namelist(&optimization_covariable, nltext);
    print_namelist(stderr, &optimization_covariable);

    /* check for valid input */
    if (name==NULL)
        bomb("element name missing in optimization_variable namelist", NULL);
    str_toupper(name);
    context = NULL;
    if (!find_element(name, &context, &(beamline->elem))) {
        fprintf(stderr, "error: cannot vary element %s--not in beamline\n", name);
        exit(1);
        }
    cp_str(&covariables->element[n_covariables], name);
    covariables->varied_type[n_covariables] = context->type;
    if (item==NULL)
        bomb("item name missing in optimization_variable list", NULL);
    str_toupper(item);
    if ((covariables->varied_param[n_covariables] = confirm_parameter(item, context->type))<0) {
        fprintf(stderr, "error: cannot vary %s--no such parameter for %s\n",item, name);
        exit(1);
        }
    cp_str(&covariables->item[n_covariables], item);
    cp_str(&covariables->varied_quan_unit[n_covariables], 
        entity_description[context->type].parameter[covariables->varied_param[n_covariables]].unit);
    if (!get_parameter_value(covariables->varied_quan_value+n_covariables, name, covariables->varied_param[n_covariables],
            context->type, beamline))
        bomb("unable to get initial value for parameter", NULL);
    covariables->varied_quan_name[n_covariables]  = tmalloc(sizeof(char)*(strlen(name)+strlen(item)+3));
    sprintf(covariables->varied_quan_name[n_covariables], "%s.%s", name, item);
    cp_str(&covariables->equation[n_covariables], equation);
    rpn_store(covariables->varied_quan_value[n_covariables] = rpn(equation),
        covariables->memory_number[n_covariables] = rpn_create_mem(covariables->varied_quan_name[n_covariables]) );
    if (rpn_check_error()) exit(1);
    fprintf(stderr, "Initial value of %s is %e %s\n", covariables->varied_quan_name[n_covariables], covariables->varied_quan_value[n_covariables],
        covariables->varied_quan_unit[n_covariables]);
    covariables->pcode[n_covariables] = gen_pcode(covariables->equation[n_covariables]);
    optimization_data->covariables.n_covariables += 1;
    log_exit("add_optimization_covariable");
    }

void add_optimization_constraint(OPTIMIZATION_DATA *optimization_data, NAMELIST_TEXT *nltext, RUN *run, LINE_LIST *beamline)
{
    long n_constraints;
    OPTIM_CONSTRAINTS *constraints;

    log_entry("add_optimization_constraint");

    constraints = &(optimization_data->constraints);
    n_constraints = constraints->n_constraints;

    constraints->quantity   = trealloc(constraints->quantity, sizeof(*constraints->quantity)*(n_constraints+1));
    constraints->index      = trealloc(constraints->index, sizeof(*constraints->index)*(n_constraints+1));
    constraints->lower      = trealloc(constraints->lower, sizeof(*constraints->lower)*(n_constraints+1));
    constraints->upper      = trealloc(constraints->upper, sizeof(*constraints->upper)*(n_constraints+1));

/*
    quantity = NULL;
    lower = upper = 0;
 */
    
    /* process namelist text */
    set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
    set_print_namelist_flags(0);
    process_namelist(&optimization_constraint, nltext);
    print_namelist(stderr, &optimization_constraint);

    /* check for valid input */
    if (quantity==NULL)
        bomb("quantity name missing in optimization_constraint namelist", NULL);
    constraints->quantity[n_constraints] = quantity;
    if ((constraints->index[n_constraints]=get_final_property_index(quantity))<0)
        bomb("no match for quantity to be constrained", NULL);
    if (lower>upper)
        bomb("lower > upper for constraint", NULL);
    constraints->lower[n_constraints] = lower;
    constraints->upper[n_constraints] = upper;

    constraints->n_constraints += 1;
    log_exit("add_optimization_constraint");
    }


void summarize_optimization_setup(OPTIMIZATION_DATA *optimization_data)
{
    OPTIM_VARIABLES *variables;
    OPTIM_CONSTRAINTS *constraints;
    OPTIM_COVARIABLES *covariables;
    long i;

    log_entry("summarize_optimization_setup");

    constraints = &(optimization_data->constraints);
    variables = &(optimization_data->variables);
    covariables = &(optimization_data->covariables);

    if (optimization_data->fp_log) {
        fprintf(optimization_data->fp_log, "\nOptimization to be performed using method '%s' in mode '%s' with tolerance %e.\n", 
            optimize_method[optimization_data->method], optimize_mode[optimization_data->mode], optimization_data->tolerance);
        fprintf(optimization_data->fp_log, "    As many as %ld function evaluations will be performed on each of %ld passes.\n",
            optimization_data->n_evaluations, optimization_data->n_passes);

        fprintf(optimization_data->fp_log, "    The quantity to be optimized is defined by the equation\n        %s\n    subject to %ld constraints%c\n",
            optimization_data->equation, constraints->n_constraints,
            constraints->n_constraints>0?':':'.');
        for (i=0; i<constraints->n_constraints; i++)
            fprintf(optimization_data->fp_log, "        %13.6e <= %10s <= %13.6e\n",
                constraints->lower[i], constraints->quantity[i], constraints->upper[i]);

        fprintf(optimization_data->fp_log, "    The following variables will be used in the optimization:\n");
        fprintf(optimization_data->fp_log, "    name      initial value   lower limit    upper limit     step size\n");
        fprintf(optimization_data->fp_log, "--------------------------------------------------------------------------\n");
        for (i=0; i<variables->n_variables; i++) {
            fprintf(optimization_data->fp_log, "%12s  %13.6e", variables->varied_quan_name[i], variables->initial_value[i]);
            if (variables->lower_limit[i]==variables->upper_limit[i])
                fprintf(optimization_data->fp_log, "        -              -      ");
            else
                fprintf(optimization_data->fp_log, "  %13.6e  %13.6e", variables->lower_limit[i], variables->upper_limit[i]);
            if (variables->orig_step[i]==0)
                fprintf(optimization_data->fp_log, "        -\n");
            else
                fprintf(optimization_data->fp_log, "  %13.6e\n", variables->orig_step[i]);
            }
        fputc('\n', optimization_data->fp_log);
        if (covariables->n_covariables) {
            fprintf(optimization_data->fp_log, "    The following covariables will be used in the optimization:\n");
            fprintf(optimization_data->fp_log, "    name      equation\n");
            fprintf(optimization_data->fp_log, "--------------------------------------------------------------------------\n");
            for (i=0; i<covariables->n_covariables; i++)
                fprintf(optimization_data->fp_log, "%12s  %s\n", covariables->varied_quan_name[i], covariables->equation[i]);
            fputc('\n', optimization_data->fp_log);
            }
        fputc('\n', optimization_data->fp_log);
        fflush(optimization_data->fp_log);
        }
    if (!log_file || optimization_data->fp_log!=stderr) {
        fprintf(stderr, "\nOptimization to be performed using method '%s' in mode '%s' with tolerance %e.\n", 
            optimize_method[optimization_data->method], optimize_mode[optimization_data->mode], optimization_data->tolerance);
        fprintf(stderr, "    As many as %ld function evaluations will be performed on each of %ld passes.\n",
            optimization_data->n_evaluations, optimization_data->n_passes);
         fprintf(stderr, "    The quantity to be optimized is defined by the equation\n        %s\n    subject to %ld constraints%c\n",
            optimization_data->equation, constraints->n_constraints,
            constraints->n_constraints>0?':':'.');
        for (i=0; i<constraints->n_constraints; i++)
            fprintf(stderr, "        %13.6e <= %10s <= %13.6e\n",
                constraints->lower[i], constraints->quantity[i], constraints->upper[i]);
        fprintf(stderr, "    The following variables will be used in the optimization:\n");
        fprintf(stderr, "    name      initial value   lower limit    upper limit     step size\n");
        fprintf(stderr, "--------------------------------------------------------------------------\n");
        for (i=0; i<variables->n_variables; i++) {
            fprintf(stderr, "%12s  %13.6e", variables->varied_quan_name[i], variables->initial_value[i]);
            if (variables->lower_limit[i]==variables->upper_limit[i])
                fprintf(stderr, "        -              -      ");
            else
                fprintf(stderr, "  %13.6e  %13.6e", variables->lower_limit[i], variables->upper_limit[i]);
            if (variables->orig_step[i]==0)
                fprintf(stderr, "        -\n");
            else
                fprintf(stderr, "  %13.6e\n", variables->orig_step[i]);
            }
        fputc('\n', stderr);
        if (covariables->n_covariables) {
            fprintf(stderr, "    The following covariables will be used in the optimization:\n");
            fprintf(stderr, "    name      equation\n");
            fprintf(stderr, "--------------------------------------------------------------------------\n");
            for (i=0; i<covariables->n_covariables; i++)
                fprintf(stderr, "%12s  %s\n", covariables->varied_quan_name[i], covariables->equation[i]);
            fputc('\n', stderr);
            }
        fputc('\n', stderr);
        }
    log_exit("summarize_optimization_setup");
    }

/* variables needed to do tracking for optimization--this data has to be held in global
 * variables like this since the optimization function has a simple calling syntax
 */
static RUN *run;
static VARY *control;
static ERROR *error;
static BEAM *beam;
static OPTIMIZATION_DATA *optimization_data;
static OUTPUT_FILES *output;
static LINE_LIST *beamline;
static long beam_type_code, n_evaluations_made, n_passes_made;
static double *final_property_value;
static long final_property_values;

static long optim_func_flags;

void do_optimize(NAMELIST_TEXT *nltext, RUN *run1, VARY *control1, ERROR *error1, LINE_LIST *beamline1, 
            BEAM *beam1, OUTPUT_FILES *output1, OPTIMIZATION_DATA *optimization_data1, long beam_type1)
{
    static long optUDFcount = 0;
    double optimization_function(double *values, long *invalid);
    void optimization_report(double result, double *value, long pass, long n_evals, long n_dim);
    OPTIM_VARIABLES *variables;
    OPTIM_COVARIABLES *covariables;
    OPTIM_CONSTRAINTS *constraints;
    double result;
    long i;

    log_entry("do_optimize");

    run               = run1;
    control           = control1;
    error             = error1;
    beamline          = beamline1;
    beam              = beam1;
    output            = output1,
    optimization_data = optimization_data1;
    beam_type_code   = beam_type1;
    n_evaluations_made = 0;
    n_passes_made      = 0;
    
    variables = &(optimization_data->variables);
    covariables = &(optimization_data->covariables);
    constraints = &(optimization_data->constraints);

    if (variables->n_variables==0)
        bomb("no variables specified for optimization", NULL);

    /* set the end-of-optimization hidden variable to 0 */
    variables->varied_quan_value[variables->n_variables] = 0;

    /* process namelist text */
    set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
    set_print_namelist_flags(0);
    process_namelist(&optimize, nltext);
    print_namelist(stderr, &optimize);

    if (summarize_setup)
        summarize_optimization_setup(optimization_data);

    for (i=0; i<variables->n_variables; i++) {
        variables->step[i] = variables->orig_step[i];
        if (variables->step[i]==0) {
            if (variables->lower_limit[i]==variables->upper_limit[i])
                fprintf(stderr, "Note: step size for %s set to %e.\n", variables->varied_quan_name[i], variables->step[i] = 1);
            else
                fprintf(stderr, "Note: step size for %s set to %e.\n", variables->varied_quan_name[i],
                    variables->step[i] = (variables->upper_limit[i]-variables->lower_limit[i])/10);
            }
        }
    
    final_property_values = count_final_properties();
    final_property_value = tmalloc(sizeof(*final_property_value)*final_property_values);

    if (!output->sums_vs_z) {
        output->sums_vs_z = tmalloc(sizeof(*output->sums_vs_z));
        output->n_z_points = 0;
        }

    variables->varied_quan_value[i] = 0;   /* end-of-optimization indicator */
    control->i_step = control->n_steps = 1;
/*    control->n_indices = 0; */
/*    control->i_vary = 1; */
    optim_func_flags = FINAL_SUMS_ONLY + INHIBIT_FILE_OUTPUT + SILENT_RUNNING;

    if (!optimization_data->UDFcreated) {
        char UDFname[100];
        sprintf(UDFname, "optUDF%ld", optUDFcount++);
        create_udf(UDFname, optimization_data->equation);
        cp_str(&optimization_data->UDFname, UDFname);
        optimization_data->UDFcreated = 1;
        }

    switch (optimization_data->method) {
        case OPTIM_METHOD_SIMPLEX:
            fputs("Starting simplex optimization.", stderr);
            if (simplexMin(&result, variables->varied_quan_value, variables->step, variables->lower_limit,
                           variables->upper_limit, variables->n_variables, optimization_data->target, 
                           optimization_data->tolerance, optimization_function, optimization_report,
                           optimization_data->n_evaluations, optimization_data->n_passes)<0) {
                if (result>optimization_data->tolerance) {
                    if (!optimization_data->soft_failure)
                        bomb("optimization unsuccessful--aborting", NULL);
                    else
                        fputs("warning: optimization unsuccessful--continuing", stderr);
                    }
                else
                    fputs("warning: maximum number of passes reached in simplex optimization", stderr);
                }
            break;
        case OPTIM_METHOD_GRID:
            fputs("Starting grid-search optimization.", stderr);
            if (!grid_search_min(&result, variables->varied_quan_value, variables->lower_limit, variables->upper_limit, variables->step,
                    variables->n_variables, optimization_function)) {
                if (!optimization_data->soft_failure)
                    bomb("optimization unsuccessful--aborting", NULL);
                else 
                    fputs("warning: optimization unsuccessful--continuing", stderr);
                }
            break;
        case OPTIM_METHOD_SAMPLE:
            fputs("Starting grid-sample optimization.", stderr);
            if (!grid_sample_min(&result, variables->varied_quan_value, variables->lower_limit, variables->upper_limit, variables->step,
                    variables->n_variables, optimization_function, optimization_data->n_evaluations*1.0)) {
                if (!optimization_data->soft_failure)
                    bomb("optimization unsuccessful--aborting", NULL);
                else
                    fputs("warning: optimization unsuccessful--continuing", stderr);
                }
            break;
        default:
            bomb("unknown optimization method code (do_optimize())", NULL);
            break;
        }

    /* evaluate once more at the optimimum point to get all parameters right and to get additional output */
    optim_func_flags = 0;
    variables->varied_quan_value[variables->n_variables] = 1;   /* indicates end-of-optimization */
    result = optimization_function(variables->varied_quan_value, &i);

    /* change values in element definitions so that new lattice can be saved */
    change_defined_parameter_values(variables->element, variables->varied_param, variables->varied_type,
            variables->varied_quan_value, variables->n_variables);
    if (control->n_elements_to_vary)
        change_defined_parameter_values(control->element, control->varied_param, control->varied_type, 
            control->varied_quan_value, control->n_elements_to_vary);
    if (covariables->n_covariables)
        change_defined_parameter_values(covariables->element, covariables->varied_param, covariables->varied_type,
            covariables->varied_quan_value, covariables->n_covariables);

    if (optimization_data->fp_log) {
        fprintf(optimization_data->fp_log, "Optimization results:\n    '%s' has value %.15g\n", optimization_data->equation, 
                optimization_data->mode==OPTIM_MODE_MAXIMUM?-result:result);
        fprintf(optimization_data->fp_log, "    A total of %ld function evaluations were made.\n", n_evaluations_made);
        if (constraints->n_constraints) {
            fprintf(optimization_data->fp_log, "Constraints:\n");
            for (i=0; i<constraints->n_constraints; i++)
                fprintf(optimization_data->fp_log, "%10s: %23.15e\n", constraints->quantity[i], 
                        final_property_value[constraints->index[i]]);
            }
        fprintf(optimization_data->fp_log, "Optimum values of variables and changes from initial values:\n");
        for (i=0; i<variables->n_variables; i++)
            fprintf(optimization_data->fp_log, "%10s: %23.15e  %23.15e\n", variables->varied_quan_name[i], 
                    variables->varied_quan_value[i], variables->varied_quan_value[i]-variables->initial_value[i]);
        for (i=0; i<covariables->n_covariables; i++)
            fprintf(optimization_data->fp_log, "%10s: %23.15e\n", covariables->varied_quan_name[i], covariables->varied_quan_value[i]);
        fflush(optimization_data->fp_log);
        }
    if (!log_file || optimization_data->fp_log!=stderr) {
        fprintf(stderr, "Optimization results:\n    '%s' has value %.15g\n", optimization_data->equation, 
                optimization_data->mode==OPTIM_MODE_MAXIMUM?-result:result);
        fprintf(stderr, "    A total of %ld function evaluations were made.\n", n_evaluations_made);
        if (constraints->n_constraints) {
            fprintf(stderr, "Constraints:\n");
            for (i=0; i<constraints->n_constraints; i++)
                fprintf(stderr, "%10s: %23.15e\n", constraints->quantity[i], final_property_value[constraints->index[i]]);
            }
        fprintf(stderr, "Optimum values of variables and changes from initial values:\n");
        for (i=0; i<variables->n_variables; i++)
            fprintf(stderr, "%10s: %23.15e  %23.15e\n", variables->varied_quan_name[i], 
                    variables->varied_quan_value[i], variables->varied_quan_value[i]-variables->initial_value[i]);
        for (i=0; i<covariables->n_covariables; i++)
            fprintf(stderr, "%10s: %23.15e\n", covariables->varied_quan_name[i], covariables->varied_quan_value[i]);
        }
    for (i=0; i<variables->n_variables; i++)
        variables->initial_value[i] = variables->varied_quan_value[i];
    log_exit("do_optimize");
    }

/* Next three lines from elegant.c: */
#define SET_AWE_BEAM     5
#define SET_BUNCHED_BEAM 6
#define SET_SDDS_BEAM   33

static char *twiss_name[22] = {
    "betax", "alphax", "nux", "etax", "etapx", 
    "betay", "alphay", "nuy", "etay", "etapy",
    "max.betax", "max.etax", "max.etapx", 
    "max.betay", "max.etay", "max.etapy",
    "min.betax", "min.etax", "min.etapx", 
    "min.betay", "min.etay", "min.etapy",
    };
static long twiss_mem[22] = {
    -1, -1, -1, -1, -1, 
    -1, -1, -1, -1, -1, 
    -1, -1, -1,
    -1, -1, -1, 
    -1, -1, -1,
    -1, -1, -1, 
    };
static char *radint_name[8] = {
    "ex0", "Sdelta0",
    "Jx", "Jy", "Jdelta",
    "taux", "tauy", "taudelta",
  } ;
static long radint_mem[8] = {
  -1, -1, 
  -1, -1, -1,
  -1, -1, -1,
} ;

double optimization_function(double *value, long *invalid)
{
    double rpn(char *expression);
    OPTIM_VARIABLES *variables;
    OPTIM_CONSTRAINTS *constraints;
    OPTIM_COVARIABLES *covariables;
    double conval, result;
    long i;
    VMATRIX *M;
    TWISS twiss_ave, twiss_min, twiss_max;
    
    log_entry("optimization_function");

    n_evaluations_made++;

    variables = &(optimization_data->variables);
    constraints = &(optimization_data->constraints);
    covariables = &(optimization_data->covariables);

    /* assert variable values and store in rpn memories */
    delete_phase_references();
    reset_special_elements(beamline);

    assert_parameter_values(variables->element, variables->varied_param, variables->varied_type,
        value, variables->n_variables, beamline);
    for (i=0; i<variables->n_variables; i++)
        rpn_store(value[i], variables->memory_number[i]);
    /* set element flags to indicate variation of parameters */
    set_element_flags(beamline, variables->element, NULL, variables->varied_type, variables->varied_param, 
            variables->n_variables, PARAMETERS_ARE_VARIED, VMATRIX_IS_VARIED, 0, 0);

    if (covariables->n_covariables) {
        /* calculate values of covariables and assert these as well */
        for (i=0; i<covariables->n_covariables; i++) {
            rpn_store(covariables->varied_quan_value[i]=rpn(covariables->pcode[i]), covariables->memory_number[i]);
            if (rpn_check_error()) exit(1);
            }
        assert_parameter_values(covariables->element, covariables->varied_param, covariables->varied_type,
            covariables->varied_quan_value, covariables->n_covariables, beamline);
        /* set element flags to indicate variation of parameters */
        set_element_flags(beamline, covariables->element, NULL, covariables->varied_type, covariables->varied_param, 
                covariables->n_covariables, PARAMETERS_ARE_VARIED, VMATRIX_IS_VARIED, 0, 0);
        }

    if (optimization_data->fp_log) {
        fprintf(optimization_data->fp_log, "new variable values for evaluation %ld of pass %ld:\n", n_evaluations_made, n_passes_made+1);
        fflush(optimization_data->fp_log);
        for (i=0; i<variables->n_variables; i++)
            fprintf(optimization_data->fp_log, "    %10s: %23.15e\n", variables->varied_quan_name[i], value[i]);
        fflush(optimization_data->fp_log);
        if (covariables->n_covariables) {
            fprintf(optimization_data->fp_log, "new covariable values:\n");
            for (i=0; i<covariables->n_covariables; i++)
                fprintf(optimization_data->fp_log, "    %10s: %23.15e\n", covariables->varied_quan_name[i], covariables->varied_quan_value[i]);
            }
        fflush(optimization_data->fp_log);
        }

    /* compute matrices for perturbed elements */
    i = compute_changed_matrices(beamline, run);
    if (i) {
        beamline->flags &= ~BEAMLINE_CONCAT_CURRENT;
        beamline->flags &= ~BEAMLINE_TWISS_CURRENT;
        beamline->flags &= ~BEAMLINE_RADINT_CURRENT;
        }
    if (i && beamline->matrix) {
        free_matrices(beamline->matrix);
        free(beamline->matrix);
        beamline->matrix = NULL;
        }
    
    if (optimization_data->fp_log) {
        fprintf(optimization_data->fp_log, "%ld matrices (re)computed\n", i);
        fflush(optimization_data->fp_log);
        }

    /* generate initial beam distribution and track it */
    switch (beam_type_code) {
      case SET_AWE_BEAM:
        bomb("beam type code of SET_AWE_BEAM in optimization_function--this shouldn't happen", NULL);
        break;
      case SET_BUNCHED_BEAM:
        new_bunched_beam(beam, run, control, output, 0);
        break;
      case SET_SDDS_BEAM:
        new_sdds_beam(beam, run, control, output, 0);
        break;
      default:
        bomb("unknown beam type code in optimization", NULL);
        break;
        }

    control->i_step++;       /* to prevent automatic regeneration of beam */
    zero_beam_sums(output->sums_vs_z, output->n_z_points+1);
    if (beamline->flags&BEAMLINE_TWISS_WANTED) {
        if (twiss_mem[0]==-1) {
            for (i=0; i<22; i++)
                twiss_mem[i] = rpn_create_mem(twiss_name[i]);
            }
        /* get twiss mode and (beta, alpha, eta, etap) for both planes */
        update_twiss_parameters(run, beamline);
        /* store twiss parameters for last element */
        for (i=0; i<5; i++) {
            rpn_store(*((&beamline->elast->twiss->betax)+i)/(i==2?PIx2:1), twiss_mem[i]);
            rpn_store(*((&beamline->elast->twiss->betay)+i)/(i==2?PIx2:1), twiss_mem[i+5]);
            }
        /* store statistics */
        compute_twiss_statistics(beamline, &twiss_ave, &twiss_min, &twiss_max);
        rpn_store(twiss_max.betax, twiss_mem[10]);
        rpn_store(twiss_max.etax,  twiss_mem[11]);
        rpn_store(twiss_max.etapx, twiss_mem[12]);
        rpn_store(twiss_max.betay, twiss_mem[13]);
        rpn_store(twiss_max.etay,  twiss_mem[14]);
        rpn_store(twiss_max.etapy, twiss_mem[15]);
        rpn_store(twiss_min.betax, twiss_mem[16]);
        rpn_store(twiss_min.etax,  twiss_mem[17]);
        rpn_store(twiss_min.etapx, twiss_mem[18]);
        rpn_store(twiss_min.betay, twiss_mem[19]);
        rpn_store(twiss_min.etay,  twiss_mem[20]);
        rpn_store(twiss_min.etapy, twiss_mem[21]);
        }
    if (beamline->flags&BEAMLINE_RADINT_WANTED) {
      if (radint_mem[0]==-1) {
        for (i=0; i<8; i++)
          radint_mem[i] = rpn_create_mem(radint_name[i]);
      }
      /* radiation integrals already updated by update_twiss_parameters above
         which is guaranteed to be called
         */
      rpn_store(beamline->radIntegrals.ex0, radint_mem[0]);
      rpn_store(beamline->radIntegrals.sigmadelta, radint_mem[1]);
      rpn_store(beamline->radIntegrals.Jx, radint_mem[2]);
      rpn_store(beamline->radIntegrals.Jy, radint_mem[3]);
      rpn_store(beamline->radIntegrals.Jdelta, radint_mem[4]);
      rpn_store(beamline->radIntegrals.taux, radint_mem[5]);
      rpn_store(beamline->radIntegrals.tauy, radint_mem[6]);
      rpn_store(beamline->radIntegrals.taudelta, radint_mem[7]);
    }
    
    for (i=0; i<variables->n_variables; i++)
      variables->varied_quan_value[i] = value[i];

    track_beam(run, control, error, variables, beamline, beam, output, optim_func_flags);

    /* compute final parameters and store in rpn memories */
    if (!output->sums_vs_z)
        bomb("sums_vs_z element of output structure is NULL--programming error (optimization_function)", NULL);
    if ((i=compute_final_properties(final_property_value, output->sums_vs_z+output->n_z_points, beam->n_to_track, beam->p0, 
                                    M=full_matrix(&(beamline->elem), run, 1), beam->particle, control->i_step))
        != final_property_values) {
        fprintf(stderr, "error: compute_final_properties computed %ld quantities when %ld were expected (optimization_function)\n",
                i, final_property_values);
        abort();
        }
    rpn_store_final_properties(final_property_value, final_property_values);
    free_matrices(M); free(M); M = NULL;

    /* check constraints */
    *invalid = 0;
    if (optimization_data->fp_log && constraints->n_constraints) {
        fprintf(optimization_data->fp_log, "    Constraints:\n");
        fflush(optimization_data->fp_log);
        }
    for (i=0; i<constraints->n_constraints; i++) {
        if (optimization_data->fp_log)
            fprintf(optimization_data->fp_log, "    %10s: %23.15e", constraints->quantity[i], 
                    final_property_value[constraints->index[i]]);
        if ((conval=final_property_value[constraints->index[i]])<constraints->lower[i] ||
                conval>constraints->upper[i]) {
            *invalid = 1;
            if (optimization_data->fp_log) {
                fprintf(optimization_data->fp_log, " ---- invalid\n\n");
                fflush(optimization_data->fp_log);
                }
            log_exit("optimization_function");
            return(0.0);
            }
        if (optimization_data->fp_log) {
            fputc('\n', optimization_data->fp_log);
            fflush(optimization_data->fp_log);
            }
        }

    /* compute and return optimized quantity */
    rpn_clear();    /* clear rpn stack */
    result = rpn(optimization_data->UDFname);
    if (isnan(result) || isinf(result)) {
      result = sqrt(DBL_MAX);
    }
    
    if (rpn_check_error())
        exit(1);

    if (optimization_data->fp_log) {
        fprintf(optimization_data->fp_log, "equation evaluates to %23.15e\n\n", result);
        fflush(optimization_data->fp_log);
        }

    if (optimization_data->mode==OPTIM_MODE_MAXIMUM)
        result *= -1;
    log_exit("optimization_function");
    
    return(result);
    }

void optimization_report(double result, double *value, long pass, long n_evals, long n_dim)
{
    OPTIM_VARIABLES *variables;
    OPTIM_COVARIABLES *covariables;
    OPTIM_CONSTRAINTS *constraints;
    long i;

    log_entry("optimization_report");

    if (!optimization_data->fp_log)  {
        log_exit("optimization_report");
        return;
        }

    variables = &(optimization_data->variables);
    covariables = &(optimization_data->covariables);
    constraints = &(optimization_data->constraints);

    fprintf(optimization_data->fp_log, "Optimization pass %ld completed:\n    '%s' has value %23.15e\n", 
        pass, optimization_data->equation, optimization_data->mode==OPTIM_MODE_MAXIMUM?-result:result);
    n_passes_made = pass;

    if (constraints->n_constraints) {
        fprintf(optimization_data->fp_log, "    Constraints:\n");
        for (i=0; i<constraints->n_constraints; i++)
            fprintf(optimization_data->fp_log, "    %10s: %23.15e\n", constraints->quantity[i], 
                    final_property_value[constraints->index[i]]);
        }

    fprintf(optimization_data->fp_log, "    Values of variables:\n");
    for (i=0; i<n_dim; i++)
        fprintf(optimization_data->fp_log, "    %10s: %23.15e\n", variables->varied_quan_name[i], value[i]);
    fflush(optimization_data->fp_log);
    log_exit("optimization_report");
    }
