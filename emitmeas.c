/* Copyright 1994 by Michael Borland and Argonne National Laboratory,
 * all rights reserved.
 */
/* program: emitmeas
 * purpose: process elegant runs to determine emittance of beam.
 *          Also processes experimental data using elegant-computed
 *          matrices.
 * Michael Borland, 1989, 1990, 1991, 1993, 1994.
 */
/*
 * $Log: not supported by cvs2svn $
 * Revision 1.2  1997/09/05 18:16:12  emery
 * Fixed a non-initialized memory error
 * which occurs when the uncertainty column isn't specified
 * in sigma_x_file
 *
 * Revision 1.1  1997/08/28 19:51:44  borland
 * First version in repository.
 *
 */
#include "mdb.h"
#include "matlib.h"
#include "match_string.h"
#include "SDDS.h"
#include "scan.h"
#include "table.h"
#include <time.h>
#include <memory.h>

#define SET_ERROR_LEVEL 0
#define SET_N_ERROR_SETS 1
#define SET_SEED 2
#define SET_AUTONICE 3
#define SET_FILTER 4
#define SET_DEVIATION_LIMIT 5
#define SET_USE_WIDTHS 6
#define SET_SIGMAX_FILE 7
#define SET_SIGMAY_FILE 8
#define SET_X_FIT_FILE 9
#define SET_Y_FIT_FILE 10
#define SET_VARIABLE_NAME 11
#define SET_RESOLUTION 12
#define SET_LIMIT_MODE 13
#define SET_IGNORE_PLANE 14
#define SET_UNCERTAINTY_FRAC 15
#define SET_MINIMUM_UNCERT 16
#define SET_ADD_RESOLUTION 17
#define SET_FIND_UNCERT 18
#define SET_FIXED_UNCERT 19
#define SET_EMITTANCE_RESULTS 20
#define SET_CONSTANT_WEIGHTING 21
#define SET_VERBOSITY 22
#define N_OPTIONS 23

char *option[N_OPTIONS] = {
    "error_level", "n_error_sets", "seed", "autonice",
    "filter", "deviation_limit", "use_widths",
    "sigma_x_file", "sigma_y_file", "x_fit_output", "y_fit_output",
    "variable_name", "resolution", "limit_mode", "ignore_plane",
    "uncertainty_fraction", "minimum_uncertainty", "add_resolution",
    "find_uncertainties", "fixed_uncertainty",
    "emittance_results", "constant_weighting", "verbosity"
    } ;

#define USAGE "emitmeas inputfile [outputfile]\n\
 [-error_level=value_in_mm[{`gaussian',n_sigmas | `uniform'}]]\n\
 [-limit_mode={`resolution' | `zero'}[{,`reject'}]\n\
 [-n_error_sets=number] [-constant_weighting]\n\
 [-seed=integer] [-filter=quantity_name,lower,upper{,`notch'}]\n\
 [-deviation_limit=level1_in_mm{,level2...}] [-use_widths]\n\
 [-sigma_x_file=filename,variable-name,sigma-name[,uncertainty-name]]\n\
 [-sigma_y_file=filename,variable-name,sigma-name[,uncertainty-name]]\n\
 [-x_fit_output=filename] [-y_fit_output=filename]\n\
 [-variable_name=name] [-resolution=x_resolution_mm,y_resolution_mm]\n\
 [-ignore_plane={`x' | `y'}] [-uncertainty_fraction=x_value,y_value]\n\
 [-minimum_uncertainty=x_value_mm,y_value_mm]] [-add_resolution] \n\
 [-find_uncertainties] [-fixed_uncertainty=x_value_mm,y_value_mm]\n\
 [-emittance_results=filename] [-verbosity=level]\n\n\
Program by Michael Borland. (This is version 16, March 1994.)"

#define N_HELP_LINES 19
static char *additional_help[N_HELP_LINES] = {
USAGE,
"\nThis program computes the emittance from numerical or experimental",
"(K1, sigma) data.",
"inputfile is an elegant \"final\" parameters output file (SDDS format),",
"which contains the R matrix and simulated beam sigmas as a function",
"of some variable or variables.",
"Optional experimental data is supplied via -sigma_?_file, and is in",
"SDDS format, with columns (K1 (1/m^2), width (m)) .",
"-error_level and -n_error_sets allow the addition of uniform random",
"    measurement errors to the sigmas.  -limit_mode allows you to specify",
"    what to do if a modified sigma is below the resolution.",
"-deviation_limit allows exclusion of data that falls to far from the fit.",
"-use_widths asks the program to use half-widths from elegant rather than",
"     sigmas.  These widths are supposed to be equal to sigma for a guassian",
"     beam, but differ from x_rms or y_rms for non-gaussian beams.",
"-resolution allows specification of half-width measurement resolution,",
"    which is added/subtracted in quadrature from the sigma or width values.",
"-variable_name allows specification of the quantity that was varied in the",
"    simulation/experiment.  If not given, the first varied quadrupole is chosen."
    } ; 

#define N_INCREMENT 10

double solve_normal_form(MATRIX *F, MATRIX *sF, MATRIX *P, MATRIX *M, MATRIX *C, double *s2_fit);
double solve_normal_form_opt(MATRIX *F, MATRIX *sF, MATRIX *P, MATRIX *M, MATRIX *C, double dev_limit,
        int *n_used, double *s2_fit);
void get_sigma_data(double *sigma, double *uncert, double *variable_data,
    int n_configs, char *file, char *variable_name, char *sigma_name, char *uncert_name,
    int do_filter, double filter_lower, double filter_upper, FILE *fpo);
void print_fit(char *filename, double *variable_data, char *variable_name,
    double *sigma2_fit, double resol, char *sigma_name, int n_pts);
double propagate_errors_for_emittance(double **Sigma, double **Covar);
void set_up_covariance_matrix(MATRIX *K, double *sigma, double *uncert, int n_configs, int equal_weights);
double estimate_uncertainty(double *uncert, MATRIX *S, MATRIX *sS, MATRIX *R, MATRIX *s2, MATRIX *K, double dev_limit, 
        int n_configs, double uncert_min, double *fit_sig2_return);
int make_tweeked_data_set(MATRIX *s2, double *sigma, double error_level, double error_sigmas, int error_type_code, 
    int n_configs, double resol, int reject_at_limit, double limit, int *n_at_resol);

static char *x_width_symbol = "BS$bx$n (m)";
static char *y_width_symbol = "BS$by$n (m)";

#define GAUSSIAN_ERRORS 0
#define UNIFORM_ERRORS  1
#define N_ERROR_TYPES   2
char *error_type[N_ERROR_TYPES] = {
    "gaussian", "uniform"
    } ;

#define LIMIT_AT_RESOLUTION 0
#define LIMIT_AT_ZERO       1
#define N_LIMIT_OPTIONS     2
char *limit_option[N_LIMIT_OPTIONS] = {
    "resolution", "zero"
    } ;

#define MAX_N_TRIES 100

main(
    int argc,
    char **argv
    )
{
    FILE *fpo, *fpe;
    SDDS_TABLE main_input;
    double *R11;        /* R11 matrix element for ith configuration */
    double *R12;        /* R12 matrix element for ith configuration */
    double *data;       /* utility pointer */
    double *sigmax, *uncertx;               /* sigma in x plane for ith configuration */
    double *R33, *R34, *sigmay, *uncerty;   /* similar data for y plane */
    int n_configs, i_config, max_n_configs, new_configs;
    MATRIX *Rx, *Ry, *s2x, *s2y;
    MATRIX *Sx, *Sy, *sSx, *sSy, *Kx, *Ky;
    double S11_sum, S11_sum2, S12_sum, S12_sum2, S22_sum, S22_sum2;
    double S33_sum, S33_sum2, S34_sum, S34_sum2, S44_sum, S44_sum2;
    double betax, alphax, betax_sum, betax_sum2, alphax_sum, alphax_sum2;
    double betay, alphay, betay_sum, betay_sum2, alphay_sum, alphay_sum2;
    int i_R11, i_R12, i_R33, i_R34, i_sx, i_sy, i_variable;
    double emitx, emity;
    SCANNED_ARG *scanned;
    int i_arg, i;
    char *input, *output, *x_fit_file, *y_fit_file;
    char *variable_name;
    double *variable_data, *x_fit_sig2, *y_fit_sig2;
    int n_error_sets, seed, i_error, error_type_code;
    double error_sigmas;
    double emitx_max, emitx_min;
    double emity_max, emity_min;
    double emitx_sum, emity_sum;
    double emitx2_sum, emity2_sum;
    double emitx_sd, emity_sd;
    double error_level, x_error_level, y_error_level;
    double md_x, md_y, *dev_limit, md;
    int i_dev, n_dev_limits;
    int nx_used, nx_used_sum;
    int ny_used, ny_used_sum;
    int n_filters, *i_filter, *notch_filter, do_filter_variable;
    char **filter_quan;
    double *filter_lower, *filter_upper, variable_filter_lower, variable_filter_upper;
    int n_good_fits_x, n_good_fits_y;
    int use_widths;
    char *sigma_x_file, *sigma_x_vname, *sigma_x_sname, *sigma_x_uname;
    char *sigma_y_file, *sigma_y_vname, *sigma_y_sname, *sigma_y_uname;
    double x_resol, y_resol;
    double x_limit, y_limit;
    int limit_code, reject_at_limit;
    int n_xresol, n_yresol, ignore_x, ignore_y;
    double x_uncert_frac, y_uncert_frac;
    double x_uncert_min, y_uncert_min;
    double x_fixed_uncert, y_fixed_uncert;
    int equal_weights_x_fit, equal_weights_y_fit, constant_weighting;
    int add_resolution, find_uncert, verbosity;

    char *x_width_name = "Sx", *y_width_name = "Sy";

    argc = scanargs(&scanned, argc, argv);
    if (argc<2 || argc>(2+N_OPTIONS)) {
        for (i=0; i<N_HELP_LINES; i++)
            puts(additional_help[i]);
        exit(1);
        }

    input = output = NULL;
    dev_limit = NULL;
    n_error_sets = 1;
    error_level = n_dev_limits = 0;
    seed = -1;
    filter_quan = NULL;
    use_widths = 0;
    sigma_x_file = sigma_y_file = NULL;
    sigma_x_vname = sigma_y_vname = NULL;
    sigma_x_sname = sigma_y_sname = NULL;
    sigma_x_uname = sigma_y_uname = NULL;
    x_fit_file = y_fit_file = NULL;
    variable_name = NULL;
    x_resol = y_resol = 0;
    n_filters = 0;
    filter_lower = filter_upper = NULL;
    filter_quan  = NULL;
    i_filter = notch_filter = NULL;
    error_type_code = UNIFORM_ERRORS;
    error_sigmas = 1;
    limit_code = LIMIT_AT_ZERO;
    reject_at_limit = 0;
    ignore_x = ignore_y = 0;
    x_uncert_frac = y_uncert_frac = 0;
    x_uncert_min = y_uncert_min = 0;
    add_resolution = find_uncert = constant_weighting = 0;
    x_fixed_uncert = y_fixed_uncert = verbosity = 0;
    fpe = NULL;
    
    for (i_arg=1; i_arg<argc; i_arg++) {
        if (scanned[i_arg].arg_type==OPTION) {
            switch (match_string(scanned[i_arg].list[0], option,
                                 N_OPTIONS, 0)) {
                /* process options here */
              case SET_ERROR_LEVEL:
                if (scanned[i_arg].n_items==2) {
                    if (!sscanf(scanned[i_arg].list[1], "%lf", &error_level))
                        bomb("invalid -error_level syntax", USAGE);
                    }
                else if (scanned[i_arg].n_items>=3 && scanned[i_arg].n_items<=4) {
                    if (!sscanf(scanned[i_arg].list[1], "%lf", &error_level))
                        bomb("invalid -error_level syntax", USAGE);
                    if ((error_type_code=match_string(scanned[i_arg].list[2], error_type, N_ERROR_TYPES, 0))<0)
                        bomb("unknown error type", USAGE);
                    if (scanned[i_arg].n_items==4) {
                        if (error_type_code!=GAUSSIAN_ERRORS || 
                            !sscanf(scanned[i_arg].list[3], "%lf", &error_sigmas) ||
                            error_sigmas<0)
                            bomb("invalid -error_level syntax", USAGE);
                        }
                    }
                else
                    bomb("invalid -error_level syntax", USAGE);
                error_level /= 1e3;
                break;
              case SET_N_ERROR_SETS:
                if (scanned[i_arg].n_items!=2 ||
                    !sscanf(scanned[i_arg].list[1], "%d", &n_error_sets) ||
                    n_error_sets<0)
                    bomb("invalid -n_error_sets syntax", USAGE);
                break;
              case SET_SEED:
                if (scanned[i_arg].n_items!=2 ||
                    !sscanf(scanned[i_arg].list[1], "%d", &seed) ||
                    seed<0)
                    bomb("invalid -seed syntax", USAGE);
                break;
              case SET_AUTONICE:
                puts("warning: autonice not supported");
                break;
              case SET_FILTER:
                n_filters += 1;
                filter_lower = trealloc(filter_lower, sizeof(*filter_lower)*n_filters);
                filter_upper = trealloc(filter_upper, sizeof(*filter_upper)*n_filters);
                filter_quan  = trealloc(filter_quan , sizeof(*filter_quan)*n_filters);
                i_filter     = trealloc(i_filter    , sizeof(*i_filter)*n_filters);
                notch_filter = trealloc(notch_filter, sizeof(*notch_filter)*n_filters);
                notch_filter[n_filters-1] = 0;
                if (!(scanned[i_arg].n_items==4  || scanned[i_arg].n_items==5) ||
                    !(filter_quan[n_filters-1]=scanned[i_arg].list[1]) ||
                    !sscanf(scanned[i_arg].list[2], "%lf", filter_lower+n_filters-1) ||
                    !sscanf(scanned[i_arg].list[3], "%lf", filter_upper+n_filters-1) ||
                    filter_lower[n_filters-1]>=filter_upper[n_filters-1] ||
                    (scanned[i_arg].n_items==5 && !(notch_filter[n_filters-1]=(scanned[i_arg].list[4][0]=='n'))))
                    bomb("invalid -filter syntax", USAGE);
                str_toupper(filter_quan[n_filters-1]);
                break;
              case SET_DEVIATION_LIMIT:
                if ((n_dev_limits=scanned[i_arg].n_items-1)<1)
                    bomb("invalid -deviation_limit syntax", USAGE);
                dev_limit = tmalloc(sizeof(double)*n_dev_limits);
                for (i=0; i<n_dev_limits; i++) {
                    if (!sscanf(scanned[i_arg].list[i+1], "%lf", dev_limit+i) ||
                        dev_limit[i]==0) 
                        bomb("invalid -deviation_limit syntax", USAGE);
                    if (dev_limit[i]>0)
                        dev_limit[i] /= 1e3;
                    }
                break;
              case SET_USE_WIDTHS:
                cp_str(&x_width_name, "Wx");
                cp_str(&y_width_name, "Wy");
                use_widths = 1;
                break;
              case SET_SIGMAX_FILE:
                if (scanned[i_arg].n_items!=4 && scanned[i_arg].n_items!=5)
                    bomb("invalid -sigma_x_file syntax", USAGE);
                sigma_x_file = scanned[i_arg].list[1];
                sigma_x_vname = scanned[i_arg].list[2];
                sigma_x_sname = scanned[i_arg].list[3];
                if (scanned[i_arg].n_items==5)
                    sigma_x_uname = scanned[i_arg].list[4];
                break;
              case SET_SIGMAY_FILE:
                if (scanned[i_arg].n_items!=4 && scanned[i_arg].n_items!=5)
                    bomb("invalid -sigma_y_file syntax", USAGE);
                sigma_y_file = scanned[i_arg].list[1];
                sigma_y_vname = scanned[i_arg].list[2];
                sigma_y_sname = scanned[i_arg].list[3];
                if (scanned[i_arg].n_items==5)
                    sigma_y_uname = scanned[i_arg].list[4];
                break;
              case SET_X_FIT_FILE:
                if (scanned[i_arg].n_items!=2 ||
                    !(x_fit_file=scanned[i_arg].list[1]) )
                    bomb("invalid -x_fit_file syntax", USAGE);
                break;
              case SET_Y_FIT_FILE:
                if (scanned[i_arg].n_items!=2 ||
                    !(y_fit_file=scanned[i_arg].list[1]) )
                    bomb("invalid -y_fit_file syntax", USAGE);
                break;
              case SET_VARIABLE_NAME:
                if (scanned[i_arg].n_items!=2)
                    bomb("invalid -variable_name syntax", USAGE);
                variable_name = scanned[i_arg].list[1];
                break;
              case SET_RESOLUTION:
                if (scanned[i_arg].n_items!=3 ||
                    !sscanf(scanned[i_arg].list[1], "%lf", &x_resol) ||
                    !sscanf(scanned[i_arg].list[2], "%lf", &y_resol) ||
                    x_resol<0 || y_resol<0) 
                    bomb("invalid -resolution syntax", USAGE);
                x_resol /= 1e3;
                y_resol /= 1e3;
                break;
              case SET_LIMIT_MODE:
                if (scanned[i_arg].n_items<2 || scanned[i_arg].n_items>4 ||
                    (limit_code=match_string(scanned[i_arg].list[1], limit_option,
                                             N_LIMIT_OPTIONS, 0))<0)
                    bomb("invalid -limit_mode syntax", USAGE);
                if (scanned[i_arg].n_items==3) {
                    if (scanned[i_arg].list[2][0]=='r')
                        reject_at_limit = 1;
                    else
                        bomb("invalid -limit_mode syntax", USAGE);
                    }
                break;
              case SET_IGNORE_PLANE:
                if (scanned[i_arg].n_items!=2 ||
                    !((ignore_x=(scanned[i_arg].list[1][0]=='x')) ||
                      (ignore_y=(scanned[i_arg].list[1][0]=='y')) ) )
                    bomb("invalid -ignore_plane syntax", USAGE);
                break;
              case SET_UNCERTAINTY_FRAC:
                if (scanned[i_arg].n_items!=3 ||
                    !sscanf(scanned[i_arg].list[1], "%lf", &x_uncert_frac) ||
                    x_uncert_frac<=0 ||
                    !sscanf(scanned[i_arg].list[1], "%lf", &y_uncert_frac) ||
                    y_uncert_frac<=0)
                    bomb("invalid -uncertainty_fraction syntax", USAGE);
                break;
              case SET_MINIMUM_UNCERT:
                if (scanned[i_arg].n_items!=3 ||
                    !sscanf(scanned[i_arg].list[1], "%lf", &x_uncert_min) ||
                    x_uncert_min<=0 ||
                    !sscanf(scanned[i_arg].list[1], "%lf", &y_uncert_min) ||
                    y_uncert_min<=0)
                    bomb("invalid -minimum_uncertainty syntax", USAGE);
                x_uncert_min /= 1e3;
                y_uncert_min /= 1e3;
                break;
              case SET_ADD_RESOLUTION:
                add_resolution = 1;
                break;
              case SET_FIND_UNCERT:
                find_uncert = 1;
                break;
              case SET_FIXED_UNCERT:
                if (scanned[i_arg].n_items!=3 ||
                    !sscanf(scanned[i_arg].list[1], "%lf", &x_fixed_uncert) ||
                    !sscanf(scanned[i_arg].list[2], "%lf", &y_fixed_uncert) ||
                    x_fixed_uncert<0 || y_fixed_uncert<0) 
                    bomb("invalid -fixed_uncertainty syntax", USAGE);
                x_fixed_uncert /= 1e3;
                y_fixed_uncert /= 1e3;
                break;
              case SET_EMITTANCE_RESULTS:
                if (scanned[i_arg].n_items!=2 ||
                    !(fpe=fopen_e(scanned[i_arg].list[1], "w", 0)))
                    bomb("invalid -emittance_results syntax", USAGE);
                break;
              case SET_CONSTANT_WEIGHTING:
                constant_weighting = 1;
                break;
              case SET_VERBOSITY:
                if (scanned[i_arg].n_items!=2 ||
                    !sscanf(scanned[i_arg].list[1], "%d", &verbosity) ||
                    verbosity<0)
                    bomb("invalid -verbosity syntax", USAGE);
                break;
              default:
                bomb("unknown option given", USAGE);
                break;
                }
            }
        else {
            if (!input)
                input = scanned[i_arg].list[0];
            else if (!output)
                output = scanned[i_arg].list[0];
            else
                bomb("too many filenames given", USAGE);
            }
        }

    if (!output)
        fpo = stderr;
    else
        fpo = fopen_e(output, "w", 0);

    if (ignore_x && ignore_y)
        bomb("can't ignore both x and y planes!", USAGE);

    if (verbosity>0)
        fprintf(fpo, "input file: %s\n", input);
    if (n_error_sets && error_level) {
        if (verbosity>0) {
            if (error_type_code==UNIFORM_ERRORS) {
                if (error_level>0)
                    fprintf(fpo, "%d error sets will be used with uniform error level of +/- %le mm\n",
                            n_error_sets, error_level*1e3);
                else 
                    fprintf(fpo, "%d error sets will be used with uniform error level to be deterined from initial fit\n",
                            n_error_sets);
                }
            else {
                if (error_level>0) 
                    fprintf(fpo, "%d error sets will be used with %.3lf-sigma gaussian distribution with sigma of %le mm\n",
                            n_error_sets, error_sigmas, error_level*1e3);
                else
                    fprintf(fpo, "%d error sets will be used with %.3lf-sigma gaussian distribution with sigma to be determined from initial fit\n",
                            n_error_sets, error_sigmas);
                }
            fprintf(fpo, "beam-sizes are limited at %s\n", limit_option[limit_code]);
            if (reject_at_limit)
                fprintf(fpo, "error set will be rejected if any beam-size is below the limit\n");
            }
        x_limit = (limit_code==LIMIT_AT_ZERO?0.0:x_resol);
        y_limit = (limit_code==LIMIT_AT_ZERO?0.0:y_resol);
        }
    else if (error_level) {
        if (verbosity>0) {
            if (error_type_code==UNIFORM_ERRORS) {
                if (error_level>0)
                    fprintf(fpo, "1 error set will be used with uniform error level of +/- %le mm\n",
                            error_level*1e3);
                else
                    fprintf(fpo, "1 error set will be used with uniform error level to be determined from initial fit\n");
                }
            else {
                if (error_level>0) 
                    fprintf(fpo, "1 error set will be used with %.3lf-sigma gaussian distribution with sigma of %le mm\n",
                            error_sigmas, error_level*1e3);
                else
                fprintf(fpo, "1 error set will be used with %.3lf-sigma gaussian distribution with sigma to be determined from initial fit\n",
                        error_sigmas);
                }
            fprintf(fpo, "beam-sizes are limited at %s\n", limit_option[limit_code]);
            if (reject_at_limit)
                fprintf(fpo, "error set will be rejected if any beam-size is below the limit\n");
            }
        x_limit = (limit_code==LIMIT_AT_ZERO?0.0:x_resol);
        y_limit = (limit_code==LIMIT_AT_ZERO?0.0:y_resol);
        }

    if (fpe)
        fprintf(fpe, "$ge$r$bx$n ($gp$r-mm-mr)\n$ge$r$by$n ($gp$r-mm-mr)\n\n\n%-10d\n",
                n_error_sets);

    if (verbosity>0) {
        if (use_widths)
            fprintf(fpo, "using beam-sizes computed from widths\n");
        else
            fprintf(fpo, "using beam-sizes computed by averaging over coordinates\n");
        fprintf(fpo, "horizontal/vertical resolution: %le/%le mm\n",
                x_resol*1e3, y_resol*1e3);
        if (n_filters) {
            fprintf(fpo, "the following filters will be applied to the data:\n");
            for (i=0; i<n_filters; i++) 
                if (!notch_filter[i])
                    fprintf(fpo, "%s: [%le, %le]\n", filter_quan[i], filter_lower[i], filter_upper[i]);
                else
                    fprintf(fpo, "%s: <%le  >%le\n", filter_quan[i], filter_lower[i], filter_upper[i]);
            fprintf(fpo, "N.B.: only the first non-notch filter is applied to data in sigma files\n");
            }
        }

    if (!SDDS_InitializeInput(&main_input, input))
        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);

    i_sx  = SDDS_GetColumnIndex(&main_input, x_width_name);
    i_sy  = SDDS_GetColumnIndex(&main_input, y_width_name);
    i_R11 = SDDS_GetColumnIndex(&main_input, "R11");
    i_R12 = SDDS_GetColumnIndex(&main_input, "R12");
    i_R33 = SDDS_GetColumnIndex(&main_input, "R33");
    i_R34 = SDDS_GetColumnIndex(&main_input, "R34");
    if (i_R11<0 || i_R12<0 || i_R33<0 || i_R34<0 || i_sx<0 || i_sy<0)
        bomb("input file has wrong format or missing data", NULL);

    if (variable_name) {
        if ((i_variable=SDDS_GetColumnIndex(&main_input, variable_name))<0)
            bomb("no match for variable for fit output/filtering", NULL);
        }
    else {
        char **column_name;
        long column_names;
        if (!(column_name = SDDS_GetColumnNames(&main_input, &column_names)))
            SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
        for (i_variable=0; i_variable<column_names; i_variable++)
            if (column_name[i_variable][0]=='Q')
                break;
        if (i_variable==column_names)
            bomb("you did not specify -variable_name, and there is no obvious choice in the input data", NULL);
        variable_name = column_name[i_variable];
        if (verbosity>0) 
            fprintf(fpo, "will use %s as independent variable\n", variable_name);
        }

    do_filter_variable = 0;
    for (i=0; i<n_filters; i++) {
        if ((i_filter[i]=SDDS_GetColumnIndex(&main_input, filter_quan[i]))<0)
            bomb("filter quantity not found in input file", NULL);
        if (i_filter[i]==i_variable && !notch_filter[i] && !do_filter_variable) {
            do_filter_variable = 1;
            variable_filter_lower = filter_lower[i];
            variable_filter_upper = filter_upper[i];
            }
        }

    i_config = n_configs = 0;
    max_n_configs = N_INCREMENT;
    R11 = tmalloc(sizeof(double)*max_n_configs);
    R12 = tmalloc(sizeof(double)*max_n_configs);
    R33 = tmalloc(sizeof(double)*max_n_configs);
    R34 = tmalloc(sizeof(double)*max_n_configs);
    sigmax = tmalloc(sizeof(double)*max_n_configs);
    sigmay = tmalloc(sizeof(double)*max_n_configs);
    uncertx = tmalloc(sizeof(double)*max_n_configs);
    uncerty = tmalloc(sizeof(double)*max_n_configs);
    variable_data = tmalloc(sizeof(double)*max_n_configs);
    x_fit_sig2 = tmalloc(sizeof(double)*max_n_configs);
    y_fit_sig2 = tmalloc(sizeof(double)*max_n_configs);

    while (SDDS_ReadTable(&main_input)>0) {
        new_configs = SDDS_CountRowsOfInterest(&main_input);
        if ((i_config+new_configs)>=max_n_configs) {
            max_n_configs += new_configs+N_INCREMENT;
            R11 = trealloc(R11, sizeof(double)*max_n_configs);
            R12 = trealloc(R12, sizeof(double)*max_n_configs);
            R33 = trealloc(R33, sizeof(double)*max_n_configs);
            R34 = trealloc(R34, sizeof(double)*max_n_configs);
            sigmax = trealloc(sigmax, sizeof(double)*max_n_configs);
            sigmay = trealloc(sigmay, sizeof(double)*max_n_configs);
            uncertx = trealloc(uncertx, sizeof(double)*max_n_configs);
            uncerty = trealloc(uncerty, sizeof(double)*max_n_configs);
            variable_data = trealloc(variable_data, sizeof(double)*max_n_configs);
            x_fit_sig2 = trealloc(x_fit_sig2, sizeof(double)*max_n_configs);
            y_fit_sig2 = trealloc(y_fit_sig2, sizeof(double)*max_n_configs);
            }
        if (!(data = SDDS_GetColumn(&main_input, "R11")))
            SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
        memcpy(R11+i_config, data, sizeof(*data)*new_configs);
        free(data);
        if (!(data = SDDS_GetColumn(&main_input, "R12")))
            SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
        memcpy(R12+i_config, data, sizeof(*data)*new_configs);
        free(data);
        if (!(data = SDDS_GetColumn(&main_input, "R33")))
            SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
        memcpy(R33+i_config, data, sizeof(*data)*new_configs);
        free(data);
        if (!(data = SDDS_GetColumn(&main_input, "R34")))
            SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
        memcpy(R34+i_config, data, sizeof(*data)*new_configs);
        free(data);

        if (!(data = SDDS_GetColumn(&main_input, x_width_name)))
            SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
        memcpy(sigmax+i_config, data, sizeof(*data)*new_configs);
        free(data);

        if (!(data = SDDS_GetColumn(&main_input, y_width_name)))
            SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
        memcpy(sigmay+i_config, data, sizeof(*data)*new_configs);
        free(data);
 
        if (!(data = SDDS_GetColumn(&main_input, variable_name)))
            SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
        memcpy(variable_data+i_config, data, sizeof(*data)*new_configs);
        free(data);

        n_configs = new_configs+i_config;
        for ( ; i_config<n_configs; i_config++) {
            uncertx[i_config] = uncerty[i_config] = 0;
            if (x_fit_file)
                x_fit_sig2[i_config] = 0;
            if (y_fit_file) 
                y_fit_sig2[i_config] = 0;
            }
        }
    if (!SDDS_Terminate(&main_input))
        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
    if (n_configs<4)
        bomb("too few configurations in data", NULL);
    if (verbosity>0)
        printf("%ld configurations read from primary input file\n", n_configs);

    if (sigma_x_file)
        get_sigma_data(sigmax, uncertx, variable_data, n_configs, sigma_x_file, sigma_x_vname,
                       sigma_x_sname, sigma_x_uname, do_filter_variable, 
                       variable_filter_lower, variable_filter_upper, fpo);
    if (sigma_y_file)
        get_sigma_data(sigmay, uncerty, variable_data, n_configs, sigma_y_file, sigma_y_vname,
                       sigma_y_sname, sigma_y_uname, do_filter_variable, 
                       variable_filter_lower, variable_filter_upper, fpo);

    if (!ignore_x) {
        if (verbosity>0)
            fprintf(fpo, "** conditions of analysis for x plane:\n");
        if (x_resol && verbosity>0) {
            if (add_resolution)
                fprintf(fpo, "    resolution will be added in quadrature to beam sizes (simulates measurement)\n");
            else 
                fprintf(fpo, "    resolution will be subtracted in quadrature from beam sizes\n");
            }
        for (i_config=0; i_config<n_configs; i_config++) {
            if (add_resolution && x_resol)  
                sigmax[i_config] = sqrt(sqr(sigmax[i_config]) + sqr(x_resol));
            }
        if (x_resol && add_resolution)
            x_resol = 0;
        if (x_uncert_frac && x_fixed_uncert)
            bomb("can't specify uncertainty fraction and fixed uncertainty together", NULL);
        if (constant_weighting) {
            if (verbosity>0) 
                fprintf(fpo, "    will use constant-weighted fit\n");
            for (i_config=0; i_config<n_configs; i_config++) {
                if (sigmax[i_config]==0)
                    bomb("can't do constant weighting--there is a sigma that is zero", NULL);
                uncertx[i_config] = 1;
                }
            }
        else if ((x_uncert_frac || x_fixed_uncert) && sigmax[0]!=0) {
            equal_weights_x_fit = 0;
            if (verbosity>0) {
                if (x_uncert_frac)
                    fprintf(fpo, "    will use fractional uncertainties for x.  Fraction = %le\n", x_uncert_frac);
                else
                    fprintf(fpo, "    will use fixed uncertainty of %le mm for x.\n", x_fixed_uncert*1e3);
                }
            for (i_config=0; i_config<n_configs; i_config++)
                uncertx[i_config] = (x_uncert_frac ? sigmax[i_config]*x_uncert_frac : x_fixed_uncert);
            }
        else if (uncertx[0]==0 || sigmax[0]==0) {
            equal_weights_x_fit = 1;
            if (verbosity>0) 
                fprintf(fpo, "    will use equal-weights fitting in x\n");
            for (i_config=0; i_config<n_configs; i_config++)
                uncertx[i_config] = 1;
            }
        else {
            equal_weights_x_fit = 0;
            if (verbosity>0) 
                fprintf(fpo, "    will use provided uncertainties for x\n");
            }
        if (x_uncert_min) {
            if (verbosity>0) 
                fprintf(fpo, "    minimum x uncertainty to be used is %le mm\n", x_uncert_min*1e3);
            for (i_config=0; i_config<n_configs; i_config++)
                if (uncertx[i_config] < x_uncert_min)
                    uncertx[i_config] = x_uncert_min;
            }                    
        }
    if (!ignore_y) {
        if (verbosity>0)
            fprintf(fpo, "** conditions of analysis for y plane: \n");
        if (y_resol && verbosity>0) {
            if (add_resolution)
                fprintf(fpo, "    resolution will be added in quadrature to beam sizes (simulates measurement)\n");
            else 
                fprintf(fpo, "    resolution will be subtracted in quadrature from beam sizes\n");
            }
        for (i_config=0; i_config<n_configs; i_config++) {
            if (add_resolution && y_resol)  
                sigmay[i_config] = sqrt(sqr(sigmay[i_config]) + sqr(y_resol));
            }
        if (y_resol && add_resolution)
            y_resol = 0;
        if (y_uncert_frac && y_fixed_uncert)
            bomb("can't specify both uncertainty fraction and fixed uncertainty", NULL);
        if (constant_weighting) {
            if (verbosity>0) 
                fprintf(fpo, "    will use constant-weighted fit\n");
            for (i_config=0; i_config<n_configs; i_config++) {
                if (sigmay[i_config]==0)
                    bomb("can't do constant weighting--there is a sigma that is zero", NULL);
                uncerty[i_config] = 1;
                }
            }
        else if ((y_uncert_frac || y_fixed_uncert) && sigmay[0]!=0) {
            equal_weights_y_fit = 0;
            if (verbosity>0) {
                if (y_uncert_frac) 
                    fprintf(fpo, "    will use fractional uncertainties for y.  Fraction = %le\n", y_uncert_frac);
                else
                    fprintf(fpo, "    will used fixed uncertainty of %le mm for y.\n", y_fixed_uncert*1e3);
                }
            for (i_config=0; i_config<n_configs; i_config++) {
                uncerty[i_config] = (y_uncert_frac ? sigmay[i_config]*y_uncert_frac : y_fixed_uncert);
                }
            }
        else if (uncerty[0]==0 || sigmay[0]==0) {
            equal_weights_y_fit = 1;
            if (verbosity>0) 
                fprintf(fpo, "    will use equal-weights fitting in y\n");
            for (i_config=0; i_config<n_configs; i_config++)
                uncerty[i_config] = 1;
            }
        else {
            equal_weights_y_fit = 0;
            if (verbosity>0) 
                fprintf(fpo, "    will use provided uncertainties for y\n");
            }
        if (y_uncert_min) {
            if (verbosity>0) 
                fprintf(fpo, "    minimum y uncertainty to be used is %le mm\n", y_uncert_min*1e3);
            for (i_config=0; i_config<n_configs; i_config++)
                if (uncerty[i_config] < y_uncert_min)
                    uncerty[i_config] = y_uncert_min;
            }                    
        }


    if (!ignore_x && verbosity>1) {
        fprintf(fpo, "x data: \n");
        if (!constant_weighting && !equal_weights_x_fit)
            for (i_config=0; i_config<n_configs; i_config++)
                fprintf(fpo, "%le %le  +/- %le\n", 
                        variable_data[i_config], sigmax[i_config], uncertx[i_config]);
        else
            for (i_config=0; i_config<n_configs; i_config++)
                fprintf(fpo, "%le %le \n", 
                        variable_data[i_config], sigmax[i_config]);
        }
    if (!ignore_y && verbosity>1) {
        fprintf(fpo, "y data: \n");
        if (!constant_weighting && !equal_weights_y_fit) 
            for (i_config=0; i_config<n_configs; i_config++)
                fprintf(fpo, "%le %le  +/- %le\n", 
                        variable_data[i_config], sigmay[i_config], uncerty[i_config]);
        else
            for (i_config=0; i_config<n_configs; i_config++)
                fprintf(fpo, "%le %le \n", 
                        variable_data[i_config], sigmay[i_config]);
        }
    if (verbosity>0) 
        fputs("\n", fpo);

    m_alloc(&Rx,  n_configs, 3);
    m_alloc(&Ry,  n_configs, 3);
    m_alloc(&s2x, n_configs, 1);
    m_alloc(&s2y, n_configs, 1);
    m_alloc(&Sx,  3, 1);
    m_alloc(&Sy,  3, 1);
    m_alloc(&sSx,  3, 3);
    m_alloc(&sSy,  3, 3);
    m_alloc(&Kx,  n_configs, n_configs);  m_zero(Kx);
    m_alloc(&Ky,  n_configs, n_configs);  m_zero(Ky);

    if (seed<0) {
        /* generate seed from system clock */
        seed = (int)time((time_t)0);
        }
    random_1(-seed);

    if (n_error_sets<=0)
        n_error_sets = 1;
    if (n_dev_limits<=0) {
        dev_limit = tmalloc(sizeof(double));
        n_dev_limits = 1;
        dev_limit[0] = 0;
        }

    for (i_config=0; i_config<n_configs; i_config++) {
        Rx->a[i_config][0]  =  sqr(R11[i_config]);
        Rx->a[i_config][1]  =  2*R11[i_config]*R12[i_config];
        Rx->a[i_config][2]  =  sqr(R12[i_config]);
        Ry->a[i_config][0]  =  sqr(R33[i_config]);
        Ry->a[i_config][1]  =  2*R33[i_config]*R34[i_config];
        Ry->a[i_config][2]  =  sqr(R34[i_config]);
        }

    x_error_level = y_error_level = 0;

    for (i_dev=0; i_dev<n_dev_limits; i_dev++) {
        if (dev_limit[i_dev]>0)
            fprintf(fpo, "\n\n*** maximum acceptable deviation from fit: %le mm\n",
                1e3*dev_limit[i_dev]);
        else if (dev_limit[i_dev]<0)
            fprintf(fpo, "\n\n*** maximum acceptable deviation from fit: %le * mean-deviation\n",
                -dev_limit[i_dev]);
        else
            dev_limit[i_dev] = 0;
        emitx_sum = emitx2_sum = 0;
        emity_sum = emity2_sum = 0;
        emitx_max = -(emitx_min = HUGE);
        emity_max = -(emity_min = HUGE);
        S11_sum = S11_sum2 = S12_sum = S12_sum2 = S22_sum = S22_sum2 = 0;
        S33_sum = S33_sum2 = S34_sum = S34_sum2 = S44_sum = S44_sum2 = 0;
        betax_sum = betax_sum2 = alphax_sum = alphax_sum2 = 0;
        betay_sum = betay_sum2 = alphay_sum = alphay_sum2 = 0;
        md_x = md_y = nx_used_sum = ny_used_sum = 0;
        n_good_fits_x = n_good_fits_y = 0;
        n_xresol = n_yresol = 0; 
        for (i_error=0; i_error<n_error_sets; i_error++) {
            for (i_config=0; i_config<n_configs; i_config++) {
                s2x->a[i_config][0] =  sqr( sigmax[i_config] ) - (add_resolution?0:sqr(x_resol));
                s2y->a[i_config][0] =  sqr( sigmay[i_config] ) - (add_resolution?0:sqr(y_resol));
                }

            if (error_level<0) {
                if (!ignore_x) {
                    if (x_error_level==0) {
                        /* must do initial fit to find error level */
                        set_up_covariance_matrix(Kx, sigmax, uncertx, n_configs, equal_weights_x_fit);
                        x_error_level = estimate_uncertainty(uncertx, Sx, sSx, Rx, s2x, Kx, dev_limit[i_dev], 
                                                n_configs, x_uncert_min, x_fit_sig2);
                        fprintf(fpo, "\n*** using error level of %le mm for x\n", x_error_level*1e3);
                        if (x_fit_file) {
                            print_fit(x_fit_file, variable_data, variable_name, x_fit_sig2, x_resol, x_width_symbol, n_configs);
                            x_fit_file = NULL;
                            }
                        }
                    if (!make_tweeked_data_set(s2x, sigmax, x_error_level, error_sigmas, error_type_code, n_configs, 
                                x_resol, reject_at_limit, x_limit, &n_xresol))  {
                        fprintf(fpo, "fatal error: failed to get acceptable error set\n");
                        exit(1);
                        }
                    }
                if (!ignore_y) {
                    if (y_error_level==0) {
                        /* must do initial fit to find error level */
                        set_up_covariance_matrix(Ky, sigmay, uncerty, n_configs, equal_weights_y_fit);
                        y_error_level = estimate_uncertainty(uncerty, Sy, sSy, Ry, s2y, Ky, dev_limit[i_dev], 
                                                n_configs,  y_uncert_min, y_fit_sig2);
                        fprintf(fpo, "\n*** using error level of %le mm for y\n", y_error_level*1e3);
                        if (y_fit_file) {
                            print_fit(y_fit_file, variable_data, variable_name, y_fit_sig2, y_resol, y_width_symbol, n_configs);
                            y_fit_file = NULL;
                            }
                        }
                    if (!make_tweeked_data_set(s2y, sigmay, y_error_level, error_sigmas, error_type_code, n_configs, 
                                y_resol, reject_at_limit, y_limit, &n_yresol))  {
                        fprintf(fpo, "fatal error: failed to get acceptable error set\n");
                        exit(1);
                        }
                    }
                }
            else if (error_level>0) {
                if (!ignore_x) {
                    if (!make_tweeked_data_set(s2x, sigmax, error_level, error_sigmas, error_type_code, n_configs, 
                                x_resol, reject_at_limit, x_limit, &n_xresol))  {
                        fprintf(fpo, "fatal error: failed to get acceptable error set\n");
                        exit(1);
                        }
                    }
                if (!ignore_y) {
                    if (!make_tweeked_data_set(s2y, sigmay, error_level, error_sigmas, error_type_code, n_configs, 
                                y_resol, reject_at_limit, y_limit, &n_yresol)) {
                        fprintf(fpo, "fatal error: failed to get acceptable error set\n");
                        exit(1);
                        }
                    }
                }

            if (!ignore_x) 
                set_up_covariance_matrix(Kx, sigmax, uncertx, n_configs, equal_weights_x_fit);
            if (!ignore_y)
                set_up_covariance_matrix(Ky, sigmay, uncerty, n_configs, equal_weights_y_fit);

            if (!ignore_x) {
                if (find_uncert) {
                    md = estimate_uncertainty(uncertx, Sx, sSx, Rx, s2x, Kx, dev_limit[i_dev], n_configs, x_uncert_min, NULL);
                    fprintf(fpo, "*** using uncertainty %le mm for x\n", md*1e3);
                    }
                md_x += solve_normal_form_opt(Sx, sSx, Rx, s2x, Kx, dev_limit[i_dev], &nx_used, x_fit_sig2);
                if (nx_used && (emitx = Sx->a[0][0]*Sx->a[2][0]-sqr(Sx->a[1][0]))>0) {
                    emitx = sqrt(emitx);
                    S11_sum += Sx->a[0][0]; S11_sum2 += sqr(Sx->a[0][0]);
                    S12_sum += Sx->a[1][0]; S12_sum2 += sqr(Sx->a[1][0]);
                    S22_sum += Sx->a[2][0]; S22_sum2 += sqr(Sx->a[2][0]);
                    betax = Sx->a[0][0]/emitx; 
                    betax_sum += betax;  betax_sum2  += sqr(betax);
                    alphax = -Sx->a[1][0]/emitx;
                    alphax_sum += alphax; alphax_sum2 += sqr(alphax);
                    emitx_sum  += emitx;
                    emitx2_sum += sqr(emitx);
                    if (emitx_max<emitx)
                        emitx_max = emitx;
                    if (emitx_min>emitx)
                        emitx_min = emitx;
                    nx_used_sum += nx_used;
                    n_good_fits_x++;
                    }
                }

            if (!ignore_y) {       
                if (find_uncert) {
                    md = estimate_uncertainty(uncerty, Sy, sSy, Ry, s2y, Ky, dev_limit[i_dev], n_configs, y_uncert_min, NULL);
                    fprintf(fpo, "*** using uncertainty %le mm for y\n", md*1e3);
                    }
                md_y += solve_normal_form_opt(Sy, sSy, Ry, s2y, Ky, dev_limit[i_dev], &ny_used, y_fit_sig2);
                if (ny_used && (emity = Sy->a[0][0]*Sy->a[2][0]-sqr(Sy->a[1][0]))>0) {
                    emity = sqrt(emity);
                    S33_sum += Sy->a[0][0]; S33_sum2 += sqr(Sy->a[0][0]);
                    S34_sum += Sy->a[1][0]; S34_sum2 += sqr(Sy->a[1][0]);
                    S44_sum += Sy->a[2][0]; S44_sum2 += sqr(Sy->a[2][0]);
                    betay = Sy->a[0][0]/emity; 
                    betay_sum += betay;  betay_sum2  += sqr(betay);
                    alphay = -Sy->a[1][0]/emity;
                    alphay_sum += alphay; alphay_sum2 += sqr(alphay);
                    emity_sum  += emity;
                    emity2_sum += sqr(emity);
                    if (emity_max<emity)
                        emity_max = emity;
                    if (emity_min>emity)
                        emity_min = emity;
                    ny_used_sum += ny_used;
                    n_good_fits_y++;
                    }
                }

            if (fpe)
                fprintf(fpe, "%le %le\n", emitx*1e6, emity*1e6);
            }
    
        if (n_error_sets>1) {
            if (!ignore_x) {
                if (n_xresol)
                    fprintf(fpo, "\nwarning: on average %.2lf of %d x points were at resolution limit\n", 
                            (1.0*n_xresol)/n_error_sets, n_configs);
                if (n_good_fits_x) {
                    fprintf(fpo, "\nx-plane emittance: %le +/- %le mm-mrad\n",
                        emitx=1e6*emitx_sum/n_good_fits_x, 
                        emitx_sd=1e6*sqrt(emitx2_sum/n_good_fits_x-sqr(emitx_sum/n_good_fits_x)) );
                    fprintf(fpo, "    min/max emittance: %le / %le mm-mrad\n",
                       1e6*emitx_min, 1e6*emitx_max);
                    fprintf(fpo, "    average weighted rms deviation for x plane: %le mm\n", 
                        1e3*md_x/n_good_fits_x);
                    fprintf(fpo, "     sigma-matrix: S11 = %le +/- %le mm^2\n",
                        1e6*S11_sum/n_good_fits_x,
                        1e6*sqrt(S11_sum2/n_good_fits_x-sqr(S11_sum/n_good_fits_x)) );
                    fprintf(fpo, "                   S12 = %le +/- %le mm-mrad\n",
                        1e6*S12_sum/n_good_fits_x,
                        1e6*sqrt(S12_sum2/n_good_fits_x-sqr(S12_sum/n_good_fits_x)) );
                    fprintf(fpo, "                   S22 = %le +/- %le mrad^2\n",
                        1e6*S22_sum/n_good_fits_x,
                        1e6*sqrt(S22_sum2/n_good_fits_x-sqr(S22_sum/n_good_fits_x)) );
                    fprintf(fpo, "                 betax = %le +/- %le m\n",
                            betax_sum/n_good_fits_x,
                            sqrt(betax_sum2/n_good_fits_x-sqr(betax_sum/n_good_fits_x)) );
                    fprintf(fpo, "                alphax = %le +/- %le m\n",
                            alphax_sum/n_good_fits_x,
                            sqrt(alphax_sum2/n_good_fits_x-sqr(alphax_sum/n_good_fits_x)) );
                    fprintf(fpo, "    average number of points used for each of %d fits: %.1lf\n",
                        n_good_fits_x, (double)nx_used_sum/n_good_fits_x);
                    }
                else
                    fprintf(fpo, "\n*** all fits failed for x plane\n");
                }
            if (!ignore_y) {
                if (n_yresol)
                    fprintf(fpo, "\nwarning: on average %.2lf of %d y points were at resolution limit\n", 
                            (1.0*n_yresol)/n_error_sets, n_configs);
                if (n_good_fits_y) {
                    fprintf(fpo, "\ny-plane emittance: %le +/- %le mm-mrad\n",
                        emity=1e6*emity_sum/n_good_fits_y, 
                        emity_sd=1e6*sqrt(emity2_sum/n_good_fits_y-sqr(emity_sum/n_good_fits_y)) );
                    fprintf(fpo, "    average weighted rms deviation for y plane: %le mm\n", 
                        1e3*md_y/n_good_fits_y);
                    fprintf(fpo, "    min/max emittance: %le / %le mm-mrad\n",
                       1e6*emity_min, 1e6*emity_max);
                    fprintf(fpo, "     sigma-matrix: S33 = %le +/- %le mm^2\n",
                        1e6*S33_sum/n_good_fits_y,
                        1e6*sqrt(S33_sum2/n_good_fits_y-sqr(S33_sum/n_good_fits_y)) );
                    fprintf(fpo, "                   S34 = %le +/- %le mm-mrad\n",
                        1e6*S34_sum/n_good_fits_y,
                        1e6*sqrt(S34_sum2/n_good_fits_y-sqr(S34_sum/n_good_fits_y)) );
                    fprintf(fpo, "                   S44 = %le +/- %le mrad^2\n",
                        1e6*S44_sum/n_good_fits_y,
                        1e6*sqrt(S44_sum2/n_good_fits_y-sqr(S44_sum/n_good_fits_y)) );
                    fprintf(fpo, "                 betay = %le +/- %le m\n",
                            betay_sum/n_good_fits_y,
                            sqrt(betay_sum2/n_good_fits_y-sqr(betay_sum/n_good_fits_y)) );
                    fprintf(fpo, "                alphay = %le +/- %le m\n",
                            alphay_sum/n_good_fits_y,
                            sqrt(alphay_sum2/n_good_fits_y-sqr(alphay_sum/n_good_fits_y)) );
                    fprintf(fpo, "    average number of points used for each of %d fits: %.1lf\n",
                        n_good_fits_y, (double)ny_used_sum/n_good_fits_y);
                    }
                else
                    fprintf(fpo, "\n*** all fits failed for y plane\n");
                }
            }
        else {
            if (!ignore_x) {
                if (n_xresol)
                    fprintf(fpo, "\nwarning: %d of %d x points were at resolution limit\n", 
                            n_xresol, n_configs);
                if (n_good_fits_x) {
                    double S11, S12, S22, sS11, sS12, sS22, s_emitx;
                    S11 = 1e6*Sx->a[0][0]; sS11 = 1e6*sqrt(sSx->a[0][0]);
                    S12 = 1e6*Sx->a[1][0]; sS12 = 1e6*sqrt(sSx->a[1][1]);
                    S22 = 1e6*Sx->a[2][0]; sS22 = 1e6*sqrt(sSx->a[2][2]);
                    if (emitx = sqrt(S11*S22 - sqr(S12))) 
                        s_emitx = 1e6*propagate_errors_for_emittance(Sx->a, sSx->a);
                    else
                        s_emitx = 0;
                    if (!equal_weights_x_fit && !constant_weighting) {
                        fprintf(fpo, "\nx-plane emittance: %le +/- %le mm-mrad\n",  emitx, s_emitx);
                        fprintf(fpo, "    weighted rms deviation for x plane: %le mm\n", 
                                1e3*md_x);
                        fprintf(fpo, "     sigma-matrix: S11 = %le +/- %le mm^2\n", S11, sS11);
                        fprintf(fpo, "                   S12 = %le +/- %le mm-mrad\n", S12, sS12);
                        fprintf(fpo, "                   S22 = %le +/- %le mrad^2\n",  S22, sS22);
                        fprintf(fpo, "     Twiss-parameters:  beta = %le +/- %le m\n", 
                                        S11/emitx, (S11/emitx)*sqrt(sqr(s_emitx/emitx)+sqr(sS11/S11)));
                        fprintf(fpo, "                       alpha = %le +/- %le\n",
                                       -S12/emitx, (S12/emitx)*sqrt(sqr(s_emitx/emitx)+sqr(sS12/S12)));
                        fprintf(fpo, "                       gamma = %le +/- %le 1/m\n",
                                        S22/emitx, (S22/emitx)*sqrt(sqr(s_emitx/emitx)+sqr(sS22/S22)));
                        fprintf(fpo, "    number of points used for fit: %d\n", nx_used);
                        }
                    else {
                        fprintf(fpo, "\nx-plane emittance: %le mm-mrad\n",  emitx);
                        fprintf(fpo, "    weighted rms deviation for x plane: %le mm\n", 
                                1e3*md_x);
                        fprintf(fpo, "     sigma-matrix: S11 = %le mm^2\n", S11);
                        fprintf(fpo, "                   S12 = %le mm-mrad\n", S12);
                        fprintf(fpo, "                   S22 = %le mrad^2\n",  S22);
                        fprintf(fpo, "     Twiss-parameters:  beta = %le m\n", S11/emitx);
                        fprintf(fpo, "                       alpha = %le \n", -S12/emitx);
                        fprintf(fpo, "                       gamma = %le 1/m\n", S22/emitx);
                        fprintf(fpo, "    number of points used for fit: %d\n", nx_used);
                        }
                    }
                else
                    fprintf(fpo, "\n*** fit failed for x plane\n");
                }
            if (!ignore_y) {
                if (n_yresol)
                    fprintf(fpo, "\nwarning: %d of %d y points were at resolution limit\n", 
                            n_yresol, n_configs);
                if (n_good_fits_y) {
                    double S33, S34, S44, sS33, sS34, sS44, s_emity;
                    S33 = 1e6*Sy->a[0][0]; sS33 = 1e6*sqrt(sSy->a[0][0]);
                    S34 = 1e6*Sy->a[1][0]; sS34 = 1e6*sqrt(sSy->a[1][1]);
                    S44 = 1e6*Sy->a[2][0]; sS44 = 1e6*sqrt(sSy->a[2][2]);
                    if (emity = sqrt(S33*S44 - sqr(S34))) 
                        s_emity = 1e6*propagate_errors_for_emittance(Sy->a, sSy->a);
                    else
                        s_emity = 0;
                    if (!equal_weights_y_fit && !constant_weighting) {
                        fprintf(fpo, "\ny-plane emittance: %le +/- %le mm-mrad\n",  emity, s_emity);
                        fprintf(fpo, "    weighted rms deviation for y plane: %le mm\n", 
                                1e3*md_y);
                        fprintf(fpo, "     sigma-matrix: S33 = %le +/- %le mm^2\n", S33, sS33);
                        fprintf(fpo, "                   S34 = %le +/- %le mm-mrad\n", S34, sS34);
                        fprintf(fpo, "                   S44 = %le +/- %le mrad^2\n",  S44, sS44);
                        fprintf(fpo, "     Twiss-parameters:  beta = %le +/- %le m\n", 
                                        S33/emity, (S33/emity)*sqrt(sqr(s_emity/emity)+sqr(sS33/S33)));
                        fprintf(fpo, "                       alpha = %le +/- %le\n",
                                       -S34/emity, (S34/emity)*sqrt(sqr(s_emity/emity)+sqr(sS34/S34)));
                        fprintf(fpo, "                       gamma = %le +/- %le 1/m\n",
                                        S44/emity, (S44/emity)*sqrt(sqr(s_emity/emity)+sqr(sS44/S44)));
                        fprintf(fpo, "    number of points used for fit: %d\n", ny_used);
                        }
                    else {
                        fprintf(fpo, "\ny-plane emittance: %le mm-mrad\n",  emity);
                        fprintf(fpo, "    weighted rms deviation for y plane: %le mm\n", 
                                1e3*md_y);
                        fprintf(fpo, "     sigma-matrix: S33 = %le mm^2\n", S33);
                        fprintf(fpo, "                   S34 = %le mm-mrad\n", S34);
                        fprintf(fpo, "                   S44 = %le mrad^2\n",  S44);
                        fprintf(fpo, "     Twiss-parameters:  beta = %le m\n", S33/emity);
                        fprintf(fpo, "                       alpha = %le \n", -S34/emity);
                        fprintf(fpo, "                       gamma = %le 1/m\n", S44/emity);
                        fprintf(fpo, "    number of points used for fit: %d\n", ny_used);
                        }
                    }
                else
                    fprintf(fpo, "\n*** fit failed for y plane\n");
                }
            }
        }
    fclose(fpo);

    if (x_fit_file)
        print_fit(x_fit_file, variable_data, variable_name, x_fit_sig2, x_resol, x_width_symbol, n_configs);

    if (y_fit_file)
        print_fit(y_fit_file, variable_data, variable_name, y_fit_sig2, y_resol, y_width_symbol, n_configs);

    }

void print_fit(
    char *filename, double *variable_data, char *variable_name,
    double *sigma2_fit, double resol, char *sigma_name,
    int n_pts
    )
{
    int i, n_good;
    FILE *fp;


    for (i=n_good=0; i<n_pts; i++)
        if (sigma2_fit[i]>=0)
            n_good++;

    fp = fopen_e(filename, "w", 0);
    fprintf(fp, "%s\n%s\nEMITMEAS fit output\n\n%d\n", variable_name, sigma_name,
        n_good);

    for (i=0; i<n_pts; i++) 
        if (sigma2_fit[i]>=0)
            fprintf(fp, "%le %le\n", variable_data[i], sqrt(sigma2_fit[i]+sqr(resol)));

    fclose(fp);
    }

double solve_normal_form_opt(
    MATRIX *F,       /* Mx1 matrix of Fit coefficients (returned) */
    MATRIX *sF,      /* MxM matrix of errors in Fit coefficients (returned) */
    MATRIX *P,       /* NxM matrix of Parameters of fit (provided) */
    MATRIX *M,       /* Nx1 column matrix of Measured quantities (provided) */
                     /* M = P.F is the equation being solved */
    MATRIX *K,      /* NxN inverse covariance matrix for Measured quantities.
                        K[i][j] = delta(i,j)/uncert[i]^2 */
    double dev_limit,/* limit on deviation for any point used in final fit */
    int *n_used,    /* number of points used in final fit */
    double *s2_fit  /* sigma^2 from fit (returned) */
    )
{
    int i, j, i_good;
    double rms_error, error;
    double *s2_fit2;
    int *index, n_good, *good;
    MATRIX *Pc, *Mc, *Kc;

    rms_error = solve_normal_form(F, sF, P, M, K, s2_fit);

    if (dev_limit==0) {
        *n_used = M->n;
        return(rms_error);
        }

    good = tmalloc(sizeof(*good)*M->n);
    if (dev_limit<0)
        dev_limit = -dev_limit*rms_error;
    for (i=n_good=0; i<M->n; i++) {
        error = fabs(sqrt(s2_fit[i]) - sqrt(M->a[i][0]));
        if (dev_limit>error) {
            n_good++;
            good[i] = 1;
            }
        else {
            s2_fit[i] = -1;  /* used to mark excluded points for caller */
            good[i] = 0;
            }
        }
    if (n_good==0) {
        *n_used = 0;
        return(0.0);
        }

    m_alloc(&Pc, n_good, P->m);
    m_alloc(&Mc, n_good, 1);
    m_alloc(&Kc, n_good, n_good);
    m_zero(Kc);
    index = tmalloc(sizeof(int)*n_good);
    s2_fit2 = tmalloc(sizeof(*s2_fit2)*n_good);

    for (i=i_good=0; i<M->n; i++) {
        if (good[i]) {
            index[i_good] = i;
            for (j=0; j<P->m; j++)
                Pc->a[i_good][j] = P->a[i][j];
            for (j=0; j<Mc->m; j++)
                Mc->a[i_good][j] = M->a[i][j];
            Kc->a[i_good][i_good] = K->a[i][i];
            i_good++;
            }
        }            


    *n_used = n_good;
    rms_error = solve_normal_form(F, sF, Pc, Mc, Kc, s2_fit2);
    for (i=0; i<n_good; i++)
        s2_fit[index[i]] = s2_fit2[i];

    free(good);
    free(s2_fit2);
    free(index);
    m_free(&Pc);
    m_free(&Mc);
    m_free(&Kc);

    return(rms_error);
    }

double solve_normal_form(
    MATRIX *F,       /* Mx1 matrix of Fit coefficients (returned) */
    MATRIX *sF,       /* MxM covariance matrix for Fit coefficients (returned) */
    MATRIX *P,       /* NxM matrix of Parameters of fit (provided) */
    MATRIX *M,       /* Nx1 column matrix of Measured sigma^2 (provided) */
    MATRIX *K,       /* NxN inverse covariance matrix for Measured quantities (provided) .
                        K[i][j] = delta(i,j)/uncert[i]^2 */
    /* M = P.F is the equation being solved                          */
    /* the solution is F = inverse(transpose(P) K P) tranpose(P) K M = T.M */
    /*           and  sF = T.K.transpose(T) */
    double *s2_fit      /* sigma^2 from fit (returned) */
    )
{
    int n, m, i;
    MATRIX *Pt, *Pt_K, *Pt_K_P, *Inv_Pt_K_P, *Inv_PtKP_PtK;
    MATRIX *Mp, *Tt, *TC, *C;
    double error, rms_error;

    n = M->n;        /* n is the number of experimental measurements */
    m = P->m;        /* m is 3 */
    m_alloc(&Pt , m, n);
    m_alloc(&Pt_K, m, n);
    m_alloc(&Pt_K_P, m, m);
    m_alloc(&Inv_Pt_K_P, m, m);
    m_alloc(&Inv_PtKP_PtK, m, n);
    m_alloc(&Mp, M->n, M->m);
    m_alloc(&Tt, n, m);
    m_alloc(&TC, m, n);
    m_alloc(&C, K->n, K->m);

    /* find the fit */
    if (!m_trans(Pt, P))
        bomb("matrix error--call was: m_trans(Pt, P)", NULL);
    if (!m_mult(Pt_K, Pt, K))
        bomb("matrix error--call was: m_mult(Pt_K, Pt, K)", NULL);
    if (!m_mult(Pt_K_P, Pt_K, P))
        bomb("matrix error--call was: m_mult(Pt_K_P, Pt_K, P)", NULL);
    if (!m_invert(Inv_Pt_K_P, Pt_K_P))
        bomb("matrix error--call was: m_invert(Inv_Pt_K_P, Pt_K_P)", NULL);
    if (!m_mult(Inv_PtKP_PtK, Inv_Pt_K_P, Pt_K))
        bomb("matrix error--call was: m_mult(Inv_PtKP_PtK, Inv_Pt_K_P, Pt_K)", NULL);
    if (!m_mult(F, Inv_PtKP_PtK, M))
        bomb("matrix error--call was: m_mult(F, Inv_PtKP_PtK, M)", NULL);
    m_zero(sF);
    if (m_invert(C, K)) {
        if (!m_trans(Tt, Inv_PtKP_PtK))
            bomb("matrix error--call was: m_trans(Tt, Inv_PtKP_PtK)", NULL);
        if (!m_mult(TC, Inv_PtKP_PtK, C))
            bomb("matrix error--call was: m_mult(TC, Inv_PtKP_PtK, C)", NULL);
        if (!m_mult(sF, TC, Tt))
            bomb("matrix error--call was: m_mult(sF, TC, Tt)", NULL);
        }

    /* evaluate the fit */
    if (!m_mult(Mp, P, F))
        bomb("matrix error--call was: m_mult(Mp, P, F)", NULL);
    for (i=rms_error=0; i<Mp->n; i++) {
        if ((s2_fit[i] = Mp->a[i][0])<0)
            bomb("bad fit--negative sigma^2!\n", NULL);
        error = sqrt(Mp->a[i][0]) - sqrt(M->a[i][0]);
        rms_error += sqr(error);
        }
    if (Mp->n>3)
        rms_error = sqrt(rms_error/(Mp->n-3));

    m_free(&Pt);
    m_free(&Pt_K);
    m_free(&Pt_K_P);
    m_free(&Inv_Pt_K_P);
    m_free(&Inv_PtKP_PtK);
    m_free(&Mp);
    m_free(&Tt);
    m_free(&TC);
    m_free(&C);

    return(rms_error);
    }

void get_sigma_data(
    double *sigma, 
    double *uncert,
    double *variable_data,   /* should match first column of input data */
    int n_configs,          /* number of data points expected */
    char *file,
    char *variable_name, char *sigma_name, char *uncert_name,
    int do_filter, double filter_lower, double filter_upper,
    FILE *fpo
    )
{
    SDDS_TABLE table;
    int i, j, n_data;
    double *vdata, *sdata, *udata;
    
    udata = NULL;
    if (!SDDS_InitializeInput(&table, file) || !SDDS_ReadTable(&table))
        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
    if ((n_data=SDDS_CountRowsOfInterest(&table))<n_configs) {
        printf("too little data in file %s--only %ld points present\n", file, n_data);
        exit(1);
        }
    if (!(vdata=SDDS_GetColumn(&table, variable_name)) || !(sdata=SDDS_GetColumn(&table, sigma_name)) ||
        (uncert_name && !(udata=SDDS_GetColumn(&table, uncert_name))))
        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
    
    for (i=j=0; i<n_configs && j<n_data; i++, j++) {
        if (do_filter) {
            if (vdata[j]<filter_lower || vdata[j]>filter_upper) {
                i--;
                continue;
                }
            }
        sigma[i] = sdata[j];
        if (uncert_name)
            uncert[i] = udata[j];
        else
            uncert[i] = 0;
        }
    if (i!=n_configs) {
        printf("insufficient data in file %s", file);
        exit(1);
        }
    free(sdata);
    free(vdata);
    if (udata) 
      free(udata);
    fprintf(fpo, "%d sigmas taken from file %s.\n", n_configs, file);
    }


double propagate_errors_for_emittance(double **Sigma, double **Covar)
{
    double S11, S12, S22, emit;
    double sum;

    S11 = Sigma[0][0];  S12 = Sigma[1][0];  S22 = Sigma[2][0];
    if (!(emit = sqrt(S11*S22-sqr(S12))))
        return(-1.0);
    sum = 0;
    /* diagonal terms */
    sum += (sqr(S22)*Covar[0][0]/4 + sqr(S11)*Covar[2][2]/4 +
                sqr(S12)*Covar[1][1])/sqr(emit);
    /* off-diagonal terms */
    sum += 2*(-S12*S22/2*Covar[0][1] + S11*S22/4*Covar[0][2] + 
              -S11*S12/2*Covar[1][2])/sqr(emit);
    return(sqrt(sum));
    }



int make_tweeked_data_set(MATRIX *s2, double *sigma, double error_level, double error_sigmas, int error_type_code, 
    int n_configs, double resol, int reject_at_limit, double limit, int *n_at_resol)
{
    int i_config, n_tries;
    double tweek;

    for (i_config=n_tries=0; i_config<n_configs && n_tries<MAX_N_TRIES; i_config++) {
        tweek = (error_type_code==UNIFORM_ERRORS?
                    2*error_level*(random_1(1)-.5):
                    gauss_rn_lim(0.0, error_level, error_sigmas, random_1));
        s2->a[i_config][0] = sqr( sigma[i_config] + tweek ) - sqr(resol);
        if (s2->a[i_config][0]<sqr(limit)) {
            if (reject_at_limit) {
                n_tries++;
                i_config--;
                continue;
                }
            else {
                *n_at_resol += 1;
                s2->a[i_config][0] = sqr(limit);
                }
            }
        }
    if (n_tries==MAX_N_TRIES)
        return(0);
    return(1);
    }

void set_up_covariance_matrix(MATRIX *K, double *sigma, double *uncert, int n_configs, int equal_weights)
/* actually, K is the inverse of the covariance matrix */
{
    int i_config;

    for (i_config=0; i_config<n_configs; i_config++) {
        if (sigma[i_config]==0 || uncert[i_config]==0 || equal_weights)
            K->a[i_config][i_config] = 1;
        else
            K->a[i_config][i_config] = 1./sqr(2*sigma[i_config]*uncert[i_config]);
        }
    }

double estimate_uncertainty(double *uncert, MATRIX *S, MATRIX *sS, MATRIX *R, MATRIX *s2, 
                MATRIX *K, double dev_limit, int n_configs, double uncert_min, double *fit_sig2_return)
{
    double md, *fit_sig2;
    int n_used, i_config;

    if (fit_sig2_return)
        fit_sig2 = fit_sig2_return;
    else
        fit_sig2 = tmalloc(sizeof(*fit_sig2)*n_configs);

    /* find initial fit with supplied covariance matrix */
    md = solve_normal_form_opt(S, sS, R, s2, K, 0.0, &n_used, fit_sig2);
    if (!md || !n_used)
        bomb("unable to find initial fit (1)", NULL);

    /* calculate new covariance matrix */
    for (i_config=0; i_config<n_configs; i_config++)
        K->a[i_config][i_config] = 1./sqr(2*md)/s2->a[i_config][0];

    /* do second fit, excluding points that lie to far out */
    md = solve_normal_form_opt(S, sS, R, s2, K, dev_limit, &n_used, fit_sig2);
    if (!md || !n_used)
        bomb("unable to find initial fit (2)", NULL);

    /* calculate new covariance matrix */
    if (uncert_min && md<uncert_min)
        md = uncert_min;
    for (i_config=0; i_config<n_configs; i_config++) {
        K->a[i_config][i_config] = 1./sqr(2*md)/s2->a[i_config][0];
        uncert[i_config] = md;
        }

    if (!fit_sig2_return)
        free(fit_sig2);

    return(md);
    }
