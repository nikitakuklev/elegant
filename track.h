/* file: track.h
 *       contains definitions used to manipulate mad-format lattices,
 *       and to do tracking through such lattices.
 *       Lengths and multipoles are expressed in meters^n, angles in
 *       radians.
 * The external arrays are in track_data.c
 *
 * Michael Borland, 1987
 */
#include <stdio.h>
#include "namelist.h"
#include "matlib.h"
#include "SDDS.h"
#include "rpn.h"

#define malloc_verify(n) 1

/* Variable-order transport matrix structure */

typedef struct {
    double *C, **R, ***T, ****Q;
    long order;
    } VMATRIX;

/* structure for storing Twiss parameters */
typedef struct {
    double betax, alphax, phix, etax, etapx, apx;
    double betay, alphay, phiy, etay, etapy, apy;
    double Cx, Cy;
    } TWISS;

/* structure for accumulating beam moments */

typedef struct {
    double sum[6];       /* sums for centroids */
    double sum2[21];     /* sums for off-diagonal first moments */
    double maxabs[4];    /* maximum values for transverse coordinates */
    long n_part;         /* number of particles */
    double z;            /* z location */
    double p0;           /* reference momentum (beta*gamma) */
    } BEAM_SUMS;

typedef struct {
    double c1, c2;          /* centroids for two coordinate planes */
    double min1, min2;      /* minimum value for each plane */
    double max1, max2;      /* maximum value for each plane */
    double s11, s12, s22;   /* <(xi-<xi>)*(xj-<xj>)> */
    double emittance;       /* sqrt(s11*s22-s22^2) */
    double S1, S2;          /* sqrt(<xi^2>) */
    } ONE_PLANE_PARAMETERS;

/* Node structure for linked-list of element definitions: */

typedef struct element_list {
    double end_pos, end_theta;
    char *name;
    char *definition_text;
    char *p_elem;       /* pointer to the element structure */
    long type;
    long occurence;     /* greater than 1, if assigned */
    long flags;
#define PARAMETERS_ARE_STATIC    0
#define PARAMETERS_ARE_VARIED    1
#define PARAMETERS_ARE_PERTURBED 2
#define VMATRIX_IS_VARIED 4
#define VMATRIX_IS_PERTURBED 8
    double Pref_input, Pref_output;
    VMATRIX *matrix;      /* matrix of this element */
    TWISS *twiss;
    VMATRIX *accumMatrix; /* accumulated matrix to the end of this element */
    char *part_of;     /* name of lowest-level line that this element is part of */
    struct element_list *pred, *succ;
    } ELEMENT_LIST;

typedef struct {
    double centroid[6];
    long n_part;
    ELEMENT_LIST *elem;
    } TRAJECTORY;

/* structure to store data on links between elements */
typedef struct {
    char **target_name;            /* names of target elements */
    ELEMENT_LIST ***target_elem;   /* arrays of pointers to target element pointers */
    long *n_targets;               /* number of targets with given name */
    char **item;                   /* names of items to be changed */
    double *initial_value;         /* initial value of the parameter */
    long *target_param;            /* parameter (item) type code */
    char **source_name;            /* names of source elements, parallel to target_name */
    long *source_position;         /* before, after, etc. */
    long *flags;
#define STATIC_LINK 1
#define DYNAMIC_LINK 2
#define POST_CORRECTION_LINK 4
#define TURN_BY_TURN_LINK 8
#define LINK_ELEMENT_DEFINITION 16
    ELEMENT_LIST ***source_elem;   /* arrays of pointers to source element pointers, parallel to target_elem */
    long *n_parameters;            /* numbers of parameters for each source */
    char **equation;               /* rpn equations for items in terms of parameters of source */
    long n_links;
    } ELEMENT_LINKS;

/* Node structure for linked-list of beamline definitions: */

typedef struct line_list {
    char *name;
    char *definition;
    ELEMENT_LIST elem;     /* linked list of elements that make up this beamline */
    long n_elems, ncat_elems;
    ELEMENT_LIST ecat;     /* linked list of concatenated elements that are equivalent to the beamline */
    long i_recirc;                /* refers to element index in elem list */
    ELEMENT_LIST *elem_recirc;    /* pointer to element in elem list */
    long icat_recirc;             /* refers to element index in ecat list */
    ELEMENT_LIST *ecat_recirc;    /* pointer to element in ecat list */
    TWISS *twiss0;         /* initial Twiss parameters */
    ELEMENT_LIST *elem_twiss;  /* pointer to element for which twiss0 holds entering Twiss parameters.
                              Usually &elem or elem_recirc. */
    ELEMENT_LIST *elast;    /* pointer to last element &elem list */
    double tune[2];          /* x and y tunes from start of elem_twiss to end of line */
    double chromaticity[2];  /* dNUx/dp and dNUy/dp */
    double acceptance[4];    /* in pi-meter-radians for x and y, plus z locations of limits (4 doubles in all) */
    char *acc_limit_name[2];  /* names of elements at which acceptance is limited, in x and y */
    TRAJECTORY *closed_orbit;  /* closed orbit, if previously calculated, starting at recirc element */
    VMATRIX *matrix;       /* matrix from start of elem_twiss to end of line */
    char *part_of;         /* name of lowest-level line that this line is part of */
    ELEMENT_LINKS *links;   /* pointer to element links for this beamline */
    struct line_list *pred, *succ;
    double revolution_length;
    long flags;
/* flags to indicate status of operations for beamline
 * X_CURRENT : operation is current
 * X_DONE    : operation has been previously done, but may not be current
 */
#define BEAMLINE_CONCAT_CURRENT 0x00000001
#define BEAMLINE_CONCAT_DONE    0x00000002
#define BEAMLINE_TWISS_CURRENT  0x00000004
#define BEAMLINE_TWISS_DONE     0x00000008
#define BEAMLINE_TWISS_WANTED   0x00000010
    } LINE_LIST;

/* structure for passing information on run conditions */

typedef struct {
    double ideal_gamma, p_central;
    long default_order, concat_order, print_statistics;
    long combine_bunch_statistics, wrap_around, tracking_updates; 
    char *runfile, *lattice, *acceptance, *centroid, *sigma, 
         *final, *output, *rootname, *losses;
    } RUN;

/* structure containing information for variation of parameters */
typedef struct {
    long at_start;               /* indicates that present state is start of variation */
    long n_indices;              /* number of indices for variation */
    long *index, *index_limit;   /* current indices, index limits */
    long n_elements_to_vary;    
    long *element_index;         /* index used for each element */
    char **element;              /* names of elements being varied, e.g., "Q1" */
    char **item;                 /* names of item to vary for each element, e.g., "K1" */
    double *initial, *final;     /* rangse to vary over */
    double *step;                /* step sizes */
    double **enumerated_value;   /* list of values to take on, if enumerated list given */
    char **varied_quan_name;     /* e.g., "Q1[K1]" */
    char **varied_quan_unit;
    long *varied_type;           /* type code of varied element, e.g., T_QUAD */
    double *varied_quan_value;   /* current value */
    long *varied_param;          /* parameter number of varied value */
    long *flags;                 /* flag bits for variation */
#define VARY_GEOMETRIC 1
    long i_vary;                 /* running count of number of vary steps taken */
    long i_step;
    long n_steps;                /* number of error sets/bunches levels */
    double bunch_frequency;      /* bunch interval, if timing is varied */
    long n_passes;               /* number of times to go through beamline */
    long new_data_read;          /* new data has been read for variation of elements */
    LINE_LIST *cell;             /* cell to be varied along with main beamline */
    } VARY;
void check_VARY_structure(VARY *_control, char *caller);

/* structure containing information for random errors */
typedef struct {
    long n_items;
    char **name;                 /* names of elements being varied, e.g., "Q1" */
    char **item;                 /* name of item to vary for each element, e.g., "K1" */
    char **quan_name;            /* full name of perturbed quantity, e.g., dQ1[K1] */
    char **quan_unit;
    double *error_level;         /* e.g., sigma for gaussian errors */
    double *error_cutoff;        /* e.g., 3 sigma */
    long *error_type;
    long *elem_type;             /* type code of element, e.g., T_QUAD */
    long *param_number;          /* parameter number of varied value */
    long *flags;                 /* flag bits follow: */      
#define BIND_ERRORS 1
#define FRACTIONAL_ERRORS 2
#define ANTIBIND_ERRORS 4
#define BIND_ERRORS_MASK (BIND_ERRORS|ANTIBIND_ERRORS)
#define PRE_CORRECTION 8
#define POST_CORRECTION 16
#define NAME_IS_LINE 32
#define NONADDITIVE_ERRORS 64
    long *bind_number;           /* how many consecutive elements to bind */
    double *unperturbed_value;   /* current value without errors */
    double *error_value;         /* current error value */
    FILE *fp_log;                /* file to log error values to */
    long new_data_read;          /* new data has been read for control of tracking */
    } ERROR;

/* see correction.c for additional explanation of the next three structures */

typedef struct {
    /* arrays for information on individual correcting elements */
    char **corr_name;                      /* names of groups of correcting elements */
    long *corr_type;                       /* type numbers */
    char **corr_param;                     /* parameter names */
    long *param_offset;                    /* offset of correcting parameter in element structure */
    long *param_index;                     /* index of correcting parameter in entity description */
    double *corr_tweek;                    /* tweek values--amount to change parameter by to get dkick/dparam */
    double *corr_limit;                    /* limiting absolute value of the parameter */
    long n_corr_types;
    } STEERING_LIST;

typedef struct {
    /* information on useful monitors and correctors */
    long *mon_index;                       /* index of monitor in trajectory array */
    long nmon, ncor;                       /* numbers of monitors and correctors */
    ELEMENT_LIST **umoni, **ucorr;         /* arrays of pointers to monitor and corrector elements */
    double *kick_coef;                     /* dkick/dparam (==1 for hkick, vkick, and hvkick elements) */
    long *sl_index;                        /* index of steering list entry for each corrector */

    /* arrays for holding corrector information for output */
    double **kick, **posi;
    /* copies of input specifications for correction */
    double corr_fraction, corr_accuracy, corr_limit, bpm_noise, default_tweek, bpm_noise_cutoff;
    long fixed_length, bpm_noise_distribution;
    /* correction matrix plus working matrices.  dK = T*Qo is the vector of corrector kicks */
    /* Cij = dX(monitor i)/dK(corrector j) */
    MATRIX *T, *Qo, *dK, *C; 
    /* information about last correction */
    long n_cycles_done, inverse_computed;
    } CORMON_DATA;

typedef struct {
    long mode;
#define TRAJECTORY_CORRECTION 0
#define ORBIT_CORRECTION 1
    long method, verbose, track_before_and_after, n_iterations, n_xy_cycles;
    long prezero_correctors, start_from_centroid, use_actual_beam, response_only;
    double clorb_accuracy;
    double clorb_iterations;
    STEERING_LIST SLx, SLy;
    CORMON_DATA *CMx, *CMy;
    TRAJECTORY **traj;
    } CORRECTION ;
    
/* structures containing information for optimization */

typedef struct {
    char **element;                       /* names of element being varied, e.g., "Q1" */
    char **item;                          /* names of item to vary, e.g., "K1" */
    double *lower_limit, *upper_limit;    /* ranges to allow variation over */
    double *step;                         /* step size */
    double *orig_step;                    /* original step size */
    char **varied_quan_name;              /* e.g., "Q1[K1]" */
    char **varied_quan_unit;
    long *varied_type;                    /* type codes of varied element, e.g., T_QUAD */
    double *varied_quan_value;            /* current values */
    long *varied_param;                   /* parameter numbers of varied parameters */
    double *initial_value;                /* initial values of varied parameters */
    long *memory_number;                  /* rpn memory numbers of varied parameters */
    long n_variables;
    } OPTIM_VARIABLES ;

typedef struct {
    char **element;                       /* names of elements being varied, e.g., "Q1" */
    char **item;                          /* names of items to vary, e.g., "K1" */
    char **equation;                      /* rpn expressions for values of the variables */
    char **pcode;                         /* Pseudo-code versions of the equations */
    char **varied_quan_name;              /* e.g., "Q1[K1]" */
    char **varied_quan_unit;
    long *varied_type;                    /* type codes of varied element, e.g., T_QUAD */
    double *varied_quan_value;            /* current values */
    long *varied_param;                   /* parameter numbers of varied parameters */
    long *memory_number;                  /* rpn memory numbers of varied parameters */
    long n_covariables;
    } OPTIM_COVARIABLES ;

typedef struct {
    char **quantity;                    /* any quantity compiled by compute_final_parameters */
    long *index;                        /* indices in array compile by compute_final_parameters */
    double *lower, *upper;              /* upper and lower limits of allowed ranges */
    long n_constraints;
    } OPTIM_CONSTRAINTS ;

#define OPTIM_MODE_MINIMUM     0
#define OPTIM_MODE_MAXIMUM     1
#define N_OPTIM_MODES          2

#define OPTIM_METHOD_SIMPLEX   0
#define OPTIM_METHOD_GRID      1
#define OPTIM_METHOD_SAMPLE    2
#define N_OPTIM_METHODS        3

typedef struct {
    long mode, method;
    double tolerance, target;
    long i_pass, n_passes, n_evaluations;
    long soft_failure, UDFcreated;
    FILE *fp_log;
    char *equation;              /* rpn equation for thing to optimize */    
    char *UDFname;
    OPTIM_VARIABLES variables;
    OPTIM_CONSTRAINTS constraints;    
    OPTIM_COVARIABLES covariables;
    long update_periodic_twiss_parameters;    /* flag: user must request this */
    long new_data_read;          /* new data has been read for optimization */
    } OPTIMIZATION_DATA;

/* structure to store particle coordinates */
typedef struct {
    double **original;      /* original particle data */
    long n_original;        /* number of particles, and also the length of all arrays */
    long n_saved;           /* number of particles saved in original array, if this is being done */
    double p0_original;     /* initial central momentum */
    double **particle;      /* current/final coordinates */
    long n_to_track;        /* initial number of particles being tracked.  Often equal to n_original, but not always. */
    long p0;                /* current/final central momentum */
    double **accepted;      /* coordinates of accepted particles, with loss info on lost particles */
    long n_accepted;        /* final number of particles being tracked. */
    } BEAM;

/* structure to hold information for output files specified in run_setup namelist */

typedef struct {
    SDDS_TABLE SDDS_output, SDDS_accept, SDDS_centroid, SDDS_sigma, SDDS_final, SDDS_losses;
    long output_initialized, accept_initialized, centroid_initialized, sigma_initialized,
         final_initialized, losses_initialized;
    BEAM_SUMS *sums_vs_z;
    long n_z_points;
    } OUTPUT_FILES;

/* structure for chromaticity correction information */
typedef struct {
    double chromx, chromy;    /* desired chromaticities */
    double strengthLimit;     /* maximum absolute value of strength */
    char **name;              /* names of sextupole families */
    long n_families;          /* number of families */
    long n_iterations;        /* number of times to repeat correction */
    MATRIX *T;                /* Nfx2 matrix to give sextupole strength changes to change 
                                 chromaticities by given amount */
    MATRIX *dK2;              /* Nfx1 matrix of sextupole strength changes */
    MATRIX *dchrom;           /* 2x1 matrix of desired chromaticity changes */
    } CHROM_CORRECTION;

/* structure for tune correction information */
typedef struct {
    double tunex, tuney;    /* desired tunes */
    char **name;            /* names of quadrupole families */
    long n_families;        /* number of families */
    long n_iterations;      /* number of times to repeat correction */
    MATRIX *T;              /* Nfx2 matrix to give quadrupole strength changes to change 
                               chromaticities by given amount */
    MATRIX *dK1;           /* Nfx1 matrix of quadrupole strength changes */
    MATRIX *dtune;         /* 2x1 matrix of desired tune changes */
    } TUNE_CORRECTION;


/* data arrays for awe dumps, found in dump_particlesX.c */
#define N_BEAM_QUANTITIES 9
extern char *beam_quan[N_BEAM_QUANTITIES];
extern char *beam_unit[N_BEAM_QUANTITIES];
#define N_FINAL_QUANTITIES 76
extern char *final_quan[N_FINAL_QUANTITIES];
extern char *final_unit[N_FINAL_QUANTITIES];

/* entity type codes */
#define T_ECOPY -32767 
#define T_RENAME -5
#define T_RETURN -4
#define T_TITLE -3
#define T_USE   -2
#define T_NODEF -1
#define N_MADCOMS 5
#define T_LINE  0
#define T_QUAD  1
#define T_SBEN  2
#define T_RBEN  3
#define T_DRIF  4
#define T_SEXT  5
#define T_OCT   6
#define T_MULT  7
#define T_SOLE  8
#define T_HCOR  9
#define T_VCOR  10
#define T_RFCA  11
#define T_ELSE  12
#define T_HMON  13
#define T_VMON  14
#define T_MONI  15
#define T_RCOL  16
#define T_ECOL  17
#define T_MARK  18
/* Explicit matrix input: */
#define T_MATR  19
/* Alpha magnet, or one half thereof: */
#define T_ALPH  20
/* TM RF DeFlector with spatially constant fields--semi-analytical solution: */
#define T_RFDF  21
/* TM-mode RF cavity using fields from SUPERFISH and NAG integrator: */
#define T_RFTM  22
/* RaMped, spatially-constant electric DeFlector--semi-analytical solution: */
#define T_RMDF  23
/* TM-mode RF cavity with spatially Constant Fields, using NAG integrator: */
#define T_TMCF  24
/* Ramped, spatially-Constant Electric deflector PLates, using NAG integrator: */
#define T_CEPL  25
#define T_WATCH  26
/* Traveling-Wave deflector PLates (TEM fields), using NAG integrator: */
#define T_TWPL 27
#define T_MALIGN 28
/* Traveling-Wave Linear Accelerator, using NAG and first space harmonic only*/
#define T_TWLA 29
/* Pepper-pot plate */
#define T_PEPPOT 30
/* Energy matching */
#define T_ENERGY 31
/* Maximum amplitudes (after elements with matrices) */
#define T_MAXAMP 32
/* Coordinate rotation */
#define T_ROTATE 33
/* Define point from which transmission is to be calculated */
#define T_TRCOUNT 34
/* Define point from which recirculation is to begin */
#define T_RECIRC 35
/* Quadrupole triangular fringe-field element */
#define T_QFRING 36
/* Scraper insertable from one side */
#define T_SCRAPER 37
/* Trajectory correction */
#define T_CENTER 38
/* kicker magnet */
#define T_KICKER 39
/* kick sextupole */
#define T_KSEXT 40
/* kick sector bend */
#define T_KSBEND 41
/* kick quadrupole */
#define T_KQUAD 42
#define T_MAGNIFY 43
#define T_SAMPLE 44
#define T_HVCOR 45
#define T_SCATTER 46
#define T_NIBEND 47
#define T_KPOLY 48
#define T_NISEPT 49
#define T_RAMPRF 50
#define T_RAMPP  51
#define T_STRAY 52
#define T_CSBEND 53
#define T_TWMTA 54
#define T_MATTER 55
#define T_RFMODE 56
#define T_TRFMODE 57
#define T_ZLONGIT 58
#define T_SREFFECTS 59
#define T_MODRF 60
#define T_BMAPXY 61
#define N_TYPES 62

extern char *entity_name[N_TYPES];
extern char *madcom_name[N_MADCOMS];
extern char *entity_text[N_TYPES];

/* number of parameters for physical elements
 * a zero indicates an unsupported element
 */
#define N_QUAD_PARAMS 9
#define N_BEND_PARAMS 21
#define N_DRIFT_PARAMS 2
#define N_SEXT_PARAMS 8
#define N_OCTU_PARAMS 0
#define N_MULT_PARAMS 12
#define N_SOLE_PARAMS 7
#define N_HCOR_PARAMS 8
#define N_VCOR_PARAMS 8
#define N_RFCA_PARAMS 9
#define N_ELSE_PARAMS 0
#define N_HMON_PARAMS 8
#define N_VMON_PARAMS 8
#define N_MONI_PARAMS 10
#define N_RCOL_PARAMS 5
#define N_ECOL_PARAMS 5
#define N_MARK_PARAMS 1
#define N_MATR_PARAMS 3
#define N_ALPH_PARAMS 10
#define N_RFDF_PARAMS 12
#define N_RFTM_PARAMS 15
#define N_RMDF_PARAMS 10
#define N_TMCF_PARAMS 18
#define N_CEPL_PARAMS 16
#define N_TWPL_PARAMS 16
#define N_WATCH_PARAMS 5
#define N_MALIGN_PARAMS 9
#define N_TWLA_PARAMS 18
#define N_PEPPOT_PARAMS 6
#define N_ENERGY_PARAMS 4
#define N_MAXAMP_PARAMS 3
#define N_ROTATE_PARAMS 1
#define N_TRCOUNT_PARAMS 1
#define N_RECIRC_PARAMS 1
#define N_QFRING_PARAMS 9
#define N_SCRAPER_PARAMS 8
#define N_CENTER_PARAMS 4
#define N_KICKER_PARAMS 9
#define N_KSEXT_PARAMS 11
#define N_KSBEND_PARAMS 27
#define N_KQUAD_PARAMS 11
#define N_MAGNIFY_PARAMS 6
#define N_SAMPLE_PARAMS 2
#define N_HVCOR_PARAMS 10
#define N_SCATTER_PARAMS 5
#define N_NIBEND_PARAMS 18
#define N_KPOLY_PARAMS 8
#define N_RAMPRF_PARAMS 9
#define N_RAMPP_PARAMS 1
#define N_NISEPT_PARAMS 9
#define N_STRAY_PARAMS 7
#define N_CSBEND_PARAMS 24
#define N_MATTER_PARAMS 3
#define N_RFMODE_PARAMS 15
#define N_TRFMODE_PARAMS 10
#define N_TWMTA_PARAMS 17
#define N_ZLONGIT_PARAMS 13
#define N_MODRF_PARAMS 13
#define N_SREFFECTS_PARAMS 8
#define N_BMAPXY_PARAMS 5

typedef struct {
    char *name;            /* parameter name */
    char *unit;            /* parameter unit */
    long type;              /* parameter data type */
    long changes_matrix;    /* logical value */
    long offset;           /* offset position in structure */
    /* one of the following will have a meaningful value--I don't use a union,
     * since these cannot be initialized (except the first member)
     */
    char *string;
    double number;
    long integer;
    } PARAMETER;

/* maximum slope and coordinate allowed for particles in certain routines
 * (for KSBENDs, KQUADs, KSEXTs, CSBENDs, and MULT elements)
 */
#define SLOPE_LIMIT 1.0L
#define COORD_LIMIT 10.0L

#define IS_DOUBLE 1
#define IS_LONG 2
#define IS_STRING 3
#define DEFAULT_FREQUENCY 2856e6
#define DEFAULT_GAP 0.01
#define DEFAULT_FINT 0.5
#define DEFAULT_N_SECTIONS 10
#define DEFAULT_ACCURACY 1e-4
#define DEFAULT_INTEG_METHOD "runge-kutta"
#define DEFAULT_FIDUCIAL_MODE "t,median"
#define DEFAULT_RAMP_TIME 1e-9
#define DEFAULT_RADIAL_OFFSET 1.0
#define DEFAULT_BETA_WAVE 1.0
#define DEFAULT_THRESHOLD 1e-12
#define DEFAULT_NIBEND_TYPE "linear"
#define DEFAULT_N_KICKS 4

/* bit definitions for flags word in ELEMENT_DESCRIPTION */
#define HAS_MATRIX      0x00000001UL
#define HAS_LENGTH      0x00000002UL
#define DONT_CONCAT     0x00000004UL
#define OFFSETS_CHECKED 0x00000008UL
#define IS_MAGNET       0x00000010UL
#define MATRIX_TRACKING 0x00000020UL
typedef struct {
    long n_params;
    unsigned long flags;
    long structure_size;      /* in bytes */
    PARAMETER *parameter;
    } ELEMENT_DESCRIPTION;

extern ELEMENT_DESCRIPTION entity_description[N_TYPES];

/* names and storage structure for quadrupole physical parameters */
extern PARAMETER quad_param[N_QUAD_PARAMS];

typedef struct {
    double length, k1, tilt, ffringe;
    double dx, dy, dz, fse;
    long order;
    } QUAD;

/* names and storage structure for bending magnet physical parameters */
extern PARAMETER bend_param[N_BEND_PARAMS];

typedef struct {
    double length, angle, k1, e1, e2, tilt;
    double k2, h1, h2, hgap, fint;
    double dx, dy, dz;
    double fse;     /* Fractional Strength Error */
    double etilt;   /* error tilt angle */
    long edge1_effects, edge2_effects;
    long order, edge_order, TRANSPORT;
    } BEND;

/* names and storage structure for drift length physical parameters */
extern PARAMETER drift_param[N_DRIFT_PARAMS];

typedef struct {
    double length;
    long order;
    } DRIFT;

/* names and storage structure for sextupole physical parameters */
extern PARAMETER sext_param[N_SEXT_PARAMS];

typedef struct {
    double length, k2, tilt;
    double dx, dy, dz, fse;
    long order;
    } SEXT;

/* names and storage structure for solenoid */
extern PARAMETER sole_param[N_SOLE_PARAMS] ;
   
typedef struct {
    double length, ks, B;
    double dx, dy, dz;
    long order;
    } SOLE;

/* names and storage structure for arbitrary multipole */
extern PARAMETER mult_param[N_MULT_PARAMS];
   
typedef struct {
    double length, KnL, tilt, bore, BnL;
    double dx, dy, dz, factor;
    long order, n_kicks, synch_rad;
    } MULT;

/* names and storage structure for horizontal corrector physical parameters */
extern PARAMETER hcor_param[N_HCOR_PARAMS] ;
   
typedef struct {
    double length, kick, tilt, b2, calibration;
    long edge_effects, order, steering;
    } HCOR;

/* names and storage structure for vertical corrector physical parameters */
extern PARAMETER vcor_param[N_VCOR_PARAMS] ;
   
typedef struct {
    double length, kick, tilt, b2, calibration;
    long edge_effects, order, steering;
    } VCOR;

/* names and storage structure for RF cavity physical parameters */
extern PARAMETER rfca_param[N_RFCA_PARAMS] ;
   
typedef struct {
    double length, volt, phase, freq, Q;
    long phase_reference, change_p0, change_t;
    char *fiducial;
    /* for internal use only: */
    long fiducial_seen;
    double phase_fiducial;
    } RFCA;


/* names and storage structure for modulated RF cavity physical parameters */
extern PARAMETER modrf_param[N_MODRF_PARAMS] ;
   
typedef struct {
    double length, volt, phase, freq, Q;
    long phase_reference;
    double amMag, amPhase, amFreq, pmMag, pmPhase, pmFreq;
    char *fiducial;
    /* for internal use only: */
    long fiducial_seen;
    double phase_fiducial; /* -omega0*t0 */
    } MODRF;

/* names and storage structure for horizontal monitor physical parameters */
extern PARAMETER hmon_param[N_HMON_PARAMS] ;
   
typedef struct {
    double length, dx, dy, weight, tilt, calibration;
    long order;
    char *readout;   /* rpn equation for x readout as function of x and y */
    } HMON;

/* names and storage structure for vertical monitor physical parameters */
extern PARAMETER vmon_param[N_VMON_PARAMS] ;
   
typedef struct {
    double length, dx, dy, weight, tilt, calibration;
    long order;
    char *readout;   /* rpn equation for y readout as function of x and y */
    } VMON;

/* names and storage structure for two-plane monitor physical parameters */
extern PARAMETER moni_param[N_MONI_PARAMS] ;
   
typedef struct {
    double length, dx, dy, weight, tilt, xcalibration, ycalibration;
    long order;
    char *x_readout, *y_readout; /* rpn equations for x and y readouts as function of actual x and y */
    } MONI;

/* names and storage structure for rectangular collimator physical parameters */
extern PARAMETER rcol_param[N_RCOL_PARAMS] ;
   
typedef struct {
    double length, x_max, y_max, dx, dy;
    } RCOL;

/* names and storage structure for elliptical collimator physical parameters */
extern PARAMETER ecol_param[N_ECOL_PARAMS] ;
   
typedef struct {
    double length, x_max, y_max, dx, dy;
    } ECOL;

/* storage structure for marker */

typedef struct {
    long fitpoint;
    /* values for internal use: */
    long init_flags;       /* 1==twiss_mem initialized, 2===centroid_mem initialized */
    long *twiss_mem;       /* betax, alphax, NUx, etax, etaxp, betay, ... */
    long *centroid_mem;    /* (x, xp, y, yp, s, dp, Pcen, n) */
    long *sigma_mem;       /* (x, xp, y, yp, s, dp) */
    } MARK;

/* storage structure for alpha magnet */
extern PARAMETER alph_param[N_ALPH_PARAMS] ;

/* x_max(m) = ALPHA_CONST*sqrt(beta*gamma/grad_B(G/cm)) */
#define ALPHA_CONST 0.750498604674380

typedef struct {
    double xmax;        /* 75.05*sqrt(beta*gamma/gradient) in meters */
    double xs1, xs2;    /* for momentum filtration */
    double dp1, dp2;    /* for momentum filtration */
    double dx, dy, dz;
    long part;
    long order;
    double gradient;   
    } ALPH;

/* names and storage structure for RF deflector cavity done with
 * semi-analytical method
 */
extern PARAMETER rfdf_param[N_RFDF_PARAMS] ;
   
typedef struct {
    double length, phase, tilt, frequency, voltage, gap;
    double time_offset;             /* equivalent to phase */
    double B_field;
    long n_sections, phase_reference;
    double dx, dy;
    /* for internal use only */
    double t_first_particle;        /* not to be set by user! */
    long   initialized;             /* ditto */
    } RFDF;

/* TM-mode RF-cavity using NAG integrator and output from SUPERFISH
 */
extern PARAMETER rftm_param[N_RFTM_PARAMS] ;

typedef struct {
    double length, frequency, phase, Ez_peak, time_offset;
    double radial_offset, tilt;
    double accuracy;
    long phase_reference, n_steps;
    double dx, dy;
    char *filename, *method, *fiducial;
    /* variables for internal use only: */
    double *fiducial_part;
    double phase0;
    double k;
    double **Er, **Bphi, **Ez;
    long nr, nz;
    double dr, dz;
    double r_min, z_min;
    } TM_MODE;

/* names and storage structure for ramped deflector plates using 
 * semi-analytical method.
 */
extern PARAMETER rmdf_param[N_RMDF_PARAMS] ;
   
typedef struct {
    double length, tilt, ramp_time, voltage, gap;
    double time_offset;             
    long n_sections, phase_reference;
    double dx, dy;
    double t_first_particle;        /* not to be set by user! */
    long   initialized;             /* ditto */
    } RMDF;

/* Constant-field TM deflector cavity.  Er, Ez, and Bphi are all constant over
 * the length of the cavity.  Uses NAG integrator.
 */
extern PARAMETER tmcf_param[N_TMCF_PARAMS];

typedef struct {
    double length, frequency, phase, time_offset;
    double radial_offset, tilt;
    double Er, Bphi, Ez;
    double accuracy;
    double x_max, y_max;
    double dx, dy;
    long phase_reference, n_steps;
    char *method, *fiducial;
    /* variables for internal use only: */
    double *fiducial_part;            
    double phase0;        /* phase at which fiducial particle reaches center */
    double k;             /* omega/c */
    } TMCF_MODE;

/* Constant electric field deflector plates using NAG integrator.
 */
extern PARAMETER cepl_param[N_CEPL_PARAMS];

typedef struct {
    /* variables set by the user (assigned values in compute_matrices): */
    double length, ramp_time, time_offset;
    double voltage, gap, static_voltage, tilt;
    double accuracy;
    double x_max, y_max, dx, dy;
    long phase_reference, n_steps;
    char *method, *fiducial;
    /* variables for internal use only: */
    double *fiducial_part;
    double tau0;        /* t/ramp_time at which fiducial particle reaches center */
    double E_scaled;    /* e.voltage.ramp_time/(gap.m.c) */
    double E_static;    /* e.static_voltage.ramp_time/(gamp.m.c) */
    double k;           /* 1/(c.ramp_time) */
    double sin_tilt, cos_tilt;
    } CE_PLATES;

/* storage structure for watch points  */
extern PARAMETER watch_param[N_WATCH_PARAMS];

#define WATCH_COORDINATES 0
#define WATCH_PARAMETERS 1
#define WATCH_CENTROIDS 2
#define WATCH_FFT 3
#define N_WATCH_MODES 4
extern char *watch_mode[N_WATCH_MODES];
#define FFT_HANNING 0
#define FFT_PARZEN 1
#define FFT_WELCH 2
#define FFT_UNIFORM 3
#define N_FFT_WINDOWS 4
extern char *fft_window_name[N_FFT_WINDOWS];

typedef struct {
    double fraction;
    long interval;
    char *filename, *label, *mode;
    /* internal variables for SDDS output */
    long initialized, count, mode_code, window_code;
    SDDS_TABLE SDDS_table;
    } WATCH;

/* Traveling wave (TEM) deflector plates using NAG integrator.
 */
extern PARAMETER twpl_param[N_TWPL_PARAMS];

typedef struct {
    /* variables set by the user (assigned values in compute_matrices): */
    double length, ramp_time, time_offset;
    double voltage, gap, static_voltage, tilt;
    double accuracy;
    double x_max, y_max;
    double dx, dy;
    long phase_reference, n_steps;
    char *method, *fiducial;
    /* variables for internal use only: */
    double *fiducial_part;
    double tau0;        /* t/ramp_time at which fiducial particle reaches center */
    double E_scaled;    /* e.voltage.ramp_time/(gap.m.c) */
    double E_static;    /* e.static_voltage.ramp_time/(gamp.m.c) */
    double k;           /* 1/(c.ramp_time) */
    double sin_tilt, cos_tilt;
    } TW_PLATES;

/* names and storage structure for misalignment physical parameters */
extern PARAMETER malign_param[N_MALIGN_PARAMS] ;

typedef struct {
    double dxp, dyp, dx, dy, dz, dt, dp, de;
    long on_pass;
    } MALIGN;

/* Traveling-Wave Linear Accelerator, using NAG and first space harmonic 
 */
extern PARAMETER twla_param[N_TWLA_PARAMS];

typedef struct {
    /* variables set by the user (assigned values in compute_matrices): */
    double length, frequency, phase, time_offset;
    double Ez, B_solenoid, accuracy;
    double x_max, y_max;
    double dx, dy, beta_wave, alpha;
    long phase_reference, n_steps, focussing;
    char *method, *fiducial; 
    /* variables for internal use only: */
    double *fiducial_part;            
    double EzS, BsolS;
    double ErS, BphiS;    /* calculated from Ez */
    double phase0;        /* phase at which fiducial particle reaches center of
                             first cell */
    double kz;            /* omega/(c*beta_wave) */
    double alphaS;        /* alpha/k */
    } TW_LINAC;

/* storage structure for pepper-pot plates */
extern PARAMETER peppot_param[N_PEPPOT_PARAMS];

typedef struct {
    double length, radii, transmission, tilt;
    double theta_rms;
    long n_holes;
    /* internal variables */
    double *x, *y;
    } PEPPOT;

/* storage structure for energy matching */
extern PARAMETER energy_param[N_ENERGY_PARAMS];

typedef struct {
    double central_energy;
    double central_momentum;
    long match_beamline;
    long match_particles;
    } ENERGY;

/* storage structure for amplitude limits */
extern PARAMETER maxamp_param[N_MAXAMP_PARAMS];

typedef struct {
    double x_max, y_max;
    long elliptical;
    } MAXAMP;

/* storage structure for beam rotation */
extern PARAMETER rotate_param[N_ROTATE_PARAMS];

typedef struct {
    double tilt;
    } ROTATE;

/* storage structure for transmission count */
extern PARAMETER trcount_param[N_TRCOUNT_PARAMS];

typedef struct {
    long dummy;
    } TRCOUNT;

/* storage structure for recirculation point */
extern PARAMETER recirc_param[N_RECIRC_PARAMS];

typedef struct {
    long i_recirc_element;
    } RECIRC;

/* names and storage structure for quadrupole fringe field parameters */
extern PARAMETER qfring_param[N_QFRING_PARAMS];

typedef struct {
    double length, k1, tilt;
    double dx, dy, dz, fse;
    long direction, order;
    } QFRING;

/* names and storage structure for beam-scraper parameters */
extern PARAMETER scraper_param[N_SCRAPER_PARAMS];

typedef struct {
    double length, position;
    double dx, dy;
    double Xo;                  /* radiation length--if zero, all particles are absorbed */
    char *insert_from;          /* one of +x, -x, +y, -y --- replaces direction */
    long elastic;               /* selects (in)elastic scattering, if Xo is nonzero */
    long direction;             /* obsolete--used internally, however */
/* scraper 'direction' indicates the side from which the scraper comes in:
 *                       1
 *                       ^ y
 *                       |
 *                       |
 *          0    x  <----+-----  2      particles travel into page
 *                       |
 *                       | 
 *                       3
 */
    } SCRAPER;

/* names and storage structure for beam-centering parameters */
extern PARAMETER center_param[N_CENTER_PARAMS];

typedef struct {
    long x, xp, y, yp;
    } CENTER;

/* names and storage structure for time-dependent kicker */
extern PARAMETER kicker_param[N_KICKER_PARAMS];

typedef struct {
    double length, angle, tilt, time_offset;
    long periodic, phase_reference, fire_on_pass;
    char *waveform, *spatial_dependence;
    /* for internal use only: */
    double *t_wf, *amp_wf;         /* amplitude vs time for waveform */
    double tmin, tmax;
    long n_wf;                     /* number of points in waveform */
    long fiducial_seen;
    double t_fiducial;
    FILE *fpdebug;
    } KICKER;

/* names and storage structure for kick sextupole physical parameters */
extern PARAMETER sext_param[N_SEXT_PARAMS];

typedef struct {
    double length, k2, tilt, bore, B;
    double dx, dy, dz, fse;
    long n_kicks, synch_rad;
    } KSEXT;

/* names and storage structure for symplectic bending magnet physical parameters */
extern PARAMETER ksbend_param[N_KSBEND_PARAMS];

typedef struct {
    double length, angle, k1, k2, k3, k4, e1, e2, tilt;
    double h1, h2, hgap, fint;
    double dx, dy, dz;
    double fse;     /* Fractional Strength Error */
    double etilt;   /* error tilt angle */
    long n_kicks, nonlinear, synch_rad;
    long edge1_effects, edge2_effects, edge_order, paraxial, TRANSPORT;
    char *method; 
    /* for internal use only */
    long flags;
    } KSBEND;

/* names and storage structure for kick quadrupole physical parameters */
extern PARAMETER kquad_param[N_KQUAD_PARAMS];

typedef struct {
    double length, k1, tilt, bore, B;
    double dx, dy, dz, fse;
    long n_kicks, synch_rad;
    } KQUAD;

/* names and storage structure for magnifier physical parameters */
typedef struct {
    double mx, mxp, my, myp, ms, mdp;
    } MAGNIFY;

/* names and storage structure for beam sample physical parameters */
typedef struct {
    double fraction;
    long interval;
    } SAMPLE;

/* names and storage structure for horizontal/vertical corrector physical parameters */
extern PARAMETER hvcor_param[N_HVCOR_PARAMS] ;
   
typedef struct {
    double length, xkick, ykick, tilt, b2, xcalibration, ycalibration;
    long edge_effects, order, steering;
    } HVCOR;

/* names and storage structure for matrix input from a file */
extern PARAMETER matr_param[N_MATR_PARAMS] ;
   
typedef struct {
    double length;
    char *filename;
    long order;
    /* for internal use only */
    long matrix_read;
    VMATRIX M;
    } MATR;

/* names and storage structure for scattering element physical parameters */
typedef struct {
    double x, xp, y, yp, dp;
    } SCATTER;

/* names and storage structure for polynomial kick terms */
extern PARAMETER kpoly_param[N_KPOLY_PARAMS];

typedef struct {
    double coefficient, tilt, dx, dy, dz, factor;
    long order;
    char *plane;
    /* for internal use only: */
    long yplane; 
    } KPOLY;

/* names and storage structure for numerically integrated bending magnet physical parameters */
extern PARAMETER nibend_param[N_NIBEND_PARAMS];

typedef struct {
    double length, angle, e1, e2, tilt;
    double dx, dy, dz;
    double fint, hgap;     /* used to calculate flen */
    double fp1, fp2;       /* fringe-field parameters 1 and 2 */
    double fse;            /* Fractional Strength Error */
    double etilt;          /* error tilt angle */
    double accuracy;
    char *model, *method;
    long synch_rad;
    /* for internal use only: */
    double flen;            /* distance from iron edge to end of fringe field */
    double rho0;            /* central bending radius */ 
    double zeta_offset;     /* offset of fringe field perpendicular to pole face necessary
                             * to get correct bend angle */
    double last_zeta_offset; 
    double x_correction;    /* used to fix spurious central trajectory offsets */
    long negative_angle;    /* used to keep track of need to invert signs before and after integration */
    } NIBEND;

/* names and storage structure for numerically integrated septum magnet physical parameters */
extern PARAMETER nisept_param[N_NISEPT_PARAMS];

typedef struct {
    double length;          /* straight-line length */
    double angle;           /* desired bending angle */
    double e1;
    double b1;              /* cartesian gradient is Bo*b1 */
    double q1_ref;          /* reference coordinate for Bo, gradient */
    double flen;            /* fringe-field length */
    double accuracy;
    char *method;
    char *model;            /* linear only for this element */
    /* for internal use only: */
    double e2;              /* e2 = angle - e1 */
    double rho_ideal;       /* length/angle */
    double rho0;            /* bending radius along q1=q1_ref */ 
    double fse_opt;         /* optimum fse to get right bending angle */
    double last_fse_opt;    
    double q1_offset;
    long negative_angle;
    } NISEPT;

/* names and storage structure for ramped RF cavity physical parameters */
extern PARAMETER ramprf_param[N_RAMPRF_PARAMS] ;
   
typedef struct {
    double length, volt, phase, freq;
    long phase_reference;        /* only meaningful if no frequency ramping */
    char *vwaveform, *pwaveform, *fwaveform, *fiducial;
    /* for internal use only: */
    double Ts;                 /* accumulated time-of-flight of central particle */
    double *t_Vf, *Vfactor;    /* (time, V/volt) pairs */
    double *t_dP, *dPhase;     /* (time, delta phase) pairs */
    double *t_ff, *ffactor;    /* (time, frequency/freq) pairs */
    long n_Vpts, n_Ppts, n_fpts;
    double phase_fiducial;
    long fiducial_seen;
    } RAMPRF;

/* names and storage structure for momentum ramp */
extern PARAMETER rampp_param[N_RAMPP_PARAMS] ;
   
typedef struct {
    char *waveform;
    /* for internal use only: */
    double Po;                 /* set to central momentum for first beam seen */
    double *t_Pf, *Pfactor;    /* (time, P/Po) pairs */
    long n_pts;
    } RAMPP;

/* names and storage structure for stray fields */
extern PARAMETER stray_param[N_STRAY_PARAMS] ;
typedef struct {
    double length;
    double lBx, lBy;         /* local field, in Tesla */
    double gBx, gBy, gBz;    /* global field, in Tesla */
    long order;
    } STRAY;

/* names and storage structure for canonically-integrated bending magnet physical parameters */
extern PARAMETER csbend_param[N_CSBEND_PARAMS];

typedef struct {
    double length, angle, k1, k2, k3, k4, e1, e2, tilt;
    double h1, h2, hgap, fint;
    double dx, dy, dz;
    double fse;     /* Fractional Strength Error */
    double etilt;   /* error tilt angle */
    long n_kicks, nonlinear, synch_rad;
    long edge1_effects, edge2_effects;
    long integration_order;
    /* for internal use only: */
    long flags;
    } CSBEND;

/* names and storage structure for traveling-wave muffin-tin accelerator */
extern PARAMETER twmta_param[N_TWMTA_PARAMS];

typedef struct {
    /* variables set by the user (assigned values in compute_matrices): */
    double length, frequency, phase;
    double Ez, accuracy;
    double x_max, y_max;
    double dx, dy, kx;
    double beta_wave, Bsol, alpha;
    long phase_reference, n_steps;
    char *method, *fiducial; 
    /* variables for internal use only: */
    double *fiducial_part;            
    double ky;
    double ExS, EyS, EzS;
    double BxS, ByS, BsolS;
    double phase0;       /* phase at which fiducial particle reaches center of
                            first cell */
    double kz;           /* omega/(c*beta_wave) */
    double Kx, Ky, Kz;   /* kx/ko etc. */
    double alphaS;       /* alpha/k */
    } TWMTA;

/* names and storage structure for matter physical parameters */
extern PARAMETER matter_param[N_MATTER_PARAMS];

typedef struct {
    double length;
    double Xo;       /* radiation length */
    long elastic;
    } MATTER;

/* names and storage structure for RF mode physical parameters */
extern PARAMETER rfmode_param[N_RFMODE_PARAMS];

typedef struct {
    double Ra, Q, freq;        /* shunt impedance, Q, frequency */
    double charge;             /* total initial charge */
    double initial_V;          /* initial voltage */
    double initial_phase;      /* phase of initial voltage */
    double beta;               /* the cavity beta (default is 0) */
    double bin_size;           /* size of charge bins */
    long n_bins;               /* number of charge bins */
    long preload;              /* preload with steady-state voltage for point bunch */
    double preload_factor;     /* factor to multiply preload voltage by--usually 1 */
    long rigid_until_pass;     /* beam is "rigid" until this pass */
    long sample_interval;      /* sample interval for record file */
    char *record;              /* name of file to record (t, V) in */
    long single_pass;          /* controls accumulation of voltage from turn-to-turn */
    /* for internal use: */
    double mp_charge;          /* charge per macroparticle */
    long initialized;          /* indicates that beam has been seen */
    double V;                  /* magnitude of voltage */
    double Vr, Vi;             /* real, imaginary components of voltage phasor at t=tlast */
    double last_t;             /* time at which last particle was seen */
    double last_phase;         /* phase at t=last_t */
    FILE *fprec;               /* pointer to file for recording (t, Vr) */
    } RFMODE;

/* names and storage structure for transverse RF mode physical parameters */
extern PARAMETER trfmode_param[N_TRFMODE_PARAMS];

typedef struct {
    double Ra, Q, freq;        /* shunt impedance, Q, frequency */
    double charge;             /* total initial charge */
    double beta;               /* the cavity beta (default is 0) */
    double bin_size;           /* size of charge bins */
    long n_bins;               /* number of charge bins */
    long sample_interval;      /* sample interval for record file */
    char *record;              /* name of file to record (t, V) in */
    long single_pass;          /* controls accumulation of voltage from turn-to-turn */
    /* for internal use: */
    double mp_charge;          /* charge per macroparticle */
    long initialized;          /* indicates that beam has been seen */
    double Vx;                 /* magnitude of voltage */
    double Vxr, Vxi;           /* real, imaginary components of voltage phasor at t=tlast */
    double Vy;                 /* magnitude of voltage */
    double Vyr, Vyi;           /* real, imaginary components of voltage phasor at t=tlast */
    double last_t;             /* time at which last particle was seen */
    double last_xphase;        /* phase at t=last_t */
    double last_yphase;        /* phase at t=last_t */
    FILE *fprec;               /* pointer to file for recording (t, Vr) */
    } TRFMODE;

/* names and storage structure for longitudinal impedance physical parameters */
extern PARAMETER zlongit_param[N_ZLONGIT_PARAMS];

typedef struct {
    double charge;             /* total initial charge */
    long broad_band;           /* flag */
    double Ra, Q, freq;        /* shunt impedance, Q, frequency */
    char *Zreal, *Zimag;       /* impedance vs frequency files */
    double bin_size;           /* size of charge bins */
    long n_bins;               /* number of charge bins--must be 2^n */
    char *wakes;               /* name of file to save wake potentials in */
    long wake_interval;        /* interval (in turns) between outupt of wakes */
    long area_weight;          /* flag to turn on area-weighting */
    long interpolate;          /* flag to turn on interpolation */
    /* for internal use: */
    long initialized;          /* indicates that files are loaded */
    double *Z;                 /* n_Z (Re Z, Im Z) pairs */
    /* variables for SDDS output of wakes */
    SDDS_TABLE SDDS_wake;
    long SDDS_wake_initialized;
    } ZLONGIT;

/* names and storage structure for SR effects */
extern PARAMETER sreffects_param[N_SREFFECTS_PARAMS];

typedef struct {
    double Jx, Jy, Jdelta, exRef, eyRef, SdeltaRef, DdeltaRef, pRef;
    } SREFFECTS;

/* names and storage structure for SR effects */
extern PARAMETER bmapxy_param[N_BMAPXY_PARAMS];

typedef struct {
  double length, strength, accuracy;
  char *method, *filename;
  /* these are set by the program when the file is read */
  long points, nx, ny;
  double *Fx, *Fy;
  double xmin, xmax, dx;
  double ymin, ymax, dy;
} BMAPXY;

/* macros for bending magnets */ 
#define SAME_BEND_PRECEDES 1 
#define SAME_BEND_FOLLOWS 2 
#define BEND_EDGE1_EFFECTS 4 
#define BEND_EDGE2_EFFECTS 8 
#define BEND_EDGE_EFFECTS (BEND_EDGE1_EFFECTS+BEND_EDGE2_EFFECTS)

#define IS_BEND(type) ((type)==T_SBEN || (type)==T_RBEN)

/* flags for run_awe_beam and run_bunched_beam */
#define TRACK_PREVIOUS_BUNCH 1

/* flags for do_tracking/track_beam flag word */
#define FINAL_SUMS_ONLY 1
#define TEST_PARTICLES 2
#define BEGIN_AT_RECIRC 4
#define TEST_PARTICLE_LOSSES 8
#define SILENT_RUNNING 16
#define TIME_DEPENDENCE_OFF 32
#define INHIBIT_FILE_OUTPUT 64

/* return values for get_reference_phase and check_reference_phase */
#define REF_PHASE_RETURNED 1
#define REF_PHASE_NOT_SET  2
#define REF_PHASE_NONEXISTENT 3

/* definitions for use by bunched beam routines */

#define GAUSSIAN_BEAM 0
#define HARD_EDGE_BEAM 1
#define UNIFORM_ELLIPSE 2
#define SHELL_BEAM 3
#define DYNAP_BEAM 4
#define LINE_BEAM 5
#define N_BEAM_TYPES 6

extern char *beam_type[N_BEAM_TYPES];

typedef struct {
    double emit, beta, alpha, eta, etap;
    double cutoff;
    long beam_type;
    double cent_posi, cent_slope;
    } TRANSVERSE;

typedef struct {
    double sigma_dp, sigma_s, dp_s_coupling;
    double cutoff;
    long beam_type;
    double cent_s, cent_dp;
    } LONGITUDINAL;

void zero_centroid(double **particle, long n_particles, long coord);
long generate_bunch(double **particle, long n_particles, TRANSVERSE *x_plane,  TRANSVERSE *y_plane,
    LONGITUDINAL *longit, long *enforce_rms_params, long limit_invar, long symmetrize, 
    long elliptical_symmetry, double Po);
void set_beam_centroids(double **particle, long offset, long n_particles, double cent_posi, 
    double cent_slope);

/* prototypes for alpha_matrix.c: */
extern VMATRIX *alpha_magnet_matrix(double gradient, double xgamma, long maximum_order,
    long part);
extern long alpha_magnet_tracking(double **particle, VMATRIX *M, ALPH *alpha, long n_part,
    double **accepted, double P_central, double z);
 
/* prototypes for awe_beam13.c: */
extern long get_particles(double ***particle, char **input, long n_input, long one_dump, long n_skip);
extern void setup_awe_beam(BEAM *beam, NAMELIST_TEXT *nltext, RUN *run, VARY *control,
    ERROR *errcon, OPTIM_VARIABLES *optim, OUTPUT_FILES *output, LINE_LIST *beamline, long n_elements);
extern long run_awe_beam(RUN *run, VARY *control, ERROR *errcon,
    LINE_LIST *beamline, long n_elements, BEAM *beam, OUTPUT_FILES *output, long flags);
extern long new_awe_beam(BEAM *beam, RUN *run, VARY *control, OUTPUT_FILES *output, long flags);
extern void finish_awe_beam(OUTPUT_FILES *output, RUN *run, VARY *control, ERROR *errcon,
    LINE_LIST *beamline, long n_elements, BEAM *beam);
extern void adjust_arrival_time_data(double **coord, long np, double Po);
 
/* prototypes for bend_matrix6.c: */
extern VMATRIX *bend_matrix(double length, double angle, double ea1, double ea2,             
    double k1, double k2, double tilt, double fint, double gap, double fse, double etilt,
    long order, long edge_order, long flags, long TRANSPORT);
extern VMATRIX *edge_matrix(double beta, double h, double n, long which_edge,             
    double gK, long order, long all_terms, long TRANSPORT);
extern VMATRIX *corrector_matrix(double length, double kick, double tilt, double b2, double calibration,
    long do_edges, long max_order);
extern VMATRIX *hvcorrector_matrix(double length, double xkick, double ykick, double tilt, double b2,
    double xcalibration, double ycalibration, long do_edges, long max_order);
extern VMATRIX *sbend_matrix(double t0, double h, double ha, double n,         
    double beta, double xgamma, long order);
 
/* prototypes for bunched_beam12.c: */
extern void setup_bunched_beam(BEAM *beam, NAMELIST_TEXT *nltext, RUN *run, VARY *control,
    ERROR *errcon, OPTIM_VARIABLES *optim, OUTPUT_FILES *output, LINE_LIST *beamline, long n_elements);
extern long new_bunched_beam(BEAM *beam, RUN *run, VARY *control, OUTPUT_FILES *output, long flags);
extern long run_bunched_beam(RUN *run, VARY *control, ERROR *errcon, OPTIM_VARIABLES *optim, LINE_LIST *beamline, long n_elements,
    BEAM *beam, OUTPUT_FILES *output, long flags);
extern void finish_bunched_beam(OUTPUT_FILES *output, RUN *run, VARY *control, ERROR *errcon, OPTIM_VARIABLES *optim,
                                LINE_LIST *beamline, long n_elements, BEAM *beam);
extern char *brief_number(double x, char *buffer);

extern long track_beam(RUN *run, VARY *control, ERROR *errcon, OPTIM_VARIABLES *optim,
    LINE_LIST *beamline, BEAM *beam, OUTPUT_FILES *output, long flags);
extern void finish_output(OUTPUT_FILES *output, RUN *run, VARY *control, ERROR *errcon, OPTIM_VARIABLES *optim,
                          LINE_LIST *beamline, long n_elements, BEAM *beam);
extern void setup_output(OUTPUT_FILES *output, RUN *run, VARY *control, ERROR *errcon, OPTIM_VARIABLES *optim,
    LINE_LIST *beamline);

/* prototypes for cfgets.c: */
extern char *cfgets(char *s, long n, FILE *fpin);
extern void delete_spaces(char *s);
extern void str_to_upper_quotes(char *s);
 
/* prototypes for check_duplic.c: */
extern void check_duplic_elem(ELEMENT_LIST **elem, ELEMENT_LIST **new_elem, long n_elems);
extern void check_duplic_line(LINE_LIST *line, LINE_LIST *new_line, long n_lines);
 
/* prototypes for compute_centroids.c: */
extern void compute_centroids(double *centroid, double **coordinates, long n_part);
extern void compute_sigmas(double *sigma, double *centroid, double **coordinates, long n_part);
extern void zero_beam_sums(BEAM_SUMS *sums, long n);
extern void accumulate_beam_sums(BEAM_SUMS *sums, double **coords, long n_part, double p_central);
extern void copy_beam_sums(BEAM_SUMS *target, BEAM_SUMS *source);
 
/* prototypes for compute_matrices13.c: */
extern VMATRIX *full_matrix(ELEMENT_LIST *elem, RUN *run, long order);
extern VMATRIX *append_full_matrix(ELEMENT_LIST *elem, RUN *run, VMATRIX *M0, long order);
VMATRIX *accumulate_matrices(ELEMENT_LIST *elem, RUN *run, VMATRIX *M0, long order, long full_matrix_only);
extern long fill_in_matrices(ELEMENT_LIST *elem, RUN *run);
extern long calculate_matrices(LINE_LIST *line, RUN *run);
extern VMATRIX *drift_matrix(double length, long order);
extern VMATRIX *sextupole_matrix(double K2, double length, long maximum_order, double tilt, double fse);
extern VMATRIX *solenoid_matrix(double length, double ks, long max_order);
extern VMATRIX *compute_matrix(ELEMENT_LIST *elem, RUN *run, VMATRIX *Mspace);
extern void set_up_watch_point(WATCH *watch, RUN *run);
extern VMATRIX *magnification_matrix(MAGNIFY *magnif);
extern void reset_special_elements(LINE_LIST *beamline);
extern VMATRIX *stray_field_matrix(double length, double *lB, double *gB, double theta, long order, double p_central);
extern VMATRIX *rf_cavity_matrix(double length, double voltage, double frequency, double phase, double *P_central, long order);

/* prototypes for concat_beamline2.c: */
extern void copy_matrices1(VMATRIX *M1,  VMATRIX *M0);
extern void free_elements1(ELEMENT_LIST *elemlist);
extern void concatenate_beamline(LINE_LIST *beamline, RUN *run);
 
/* prototypes for concat_mat.c: */
extern void concat_matrices(VMATRIX *M2, VMATRIX *M1, VMATRIX *M0);
 
/* prototypes for copy_particles.c: */
extern void copy_particles(double **copy, double **original, long n_particles);
 
/* prototypes for correct.c: */
extern void correction_setup(CORRECTION *_correct, NAMELIST_TEXT *nltext, RUN *run, LINE_LIST *beamline);
double computeMonitorReading(ELEMENT_LIST *elem, long coord, double x, double y, 
                             unsigned long flags);
#define COMPUTEMONITORREADING_TILT_0 0x0001UL
#define COMPUTEMONITORREADING_CAL_1  0x0002UL
void setMonitorCalibration(ELEMENT_LIST *elem, double calib, long coord);
double getMonitorCalibration(ELEMENT_LIST *elem, long coord);

extern long do_correction(CORRECTION *correct, RUN *run, LINE_LIST *beamline, double *starting_coords, 
        BEAM *beam, long sim_step);
extern long find_closed_orbit(TRAJECTORY *clorb, double clorb_acc, long clorb_iter, LINE_LIST *beamline, VMATRIX *M, 
    RUN *run, double dp, long start_from_recirc, long fixed_length, double *starting_point, double iter_fraction);
extern void add_steering_element(CORRECTION *correct, LINE_LIST *beamline, RUN *run, NAMELIST_TEXT *nltext);
extern void rotate_xy(double *x, double *y, double angle);
extern void compute_trajcor_matrices(CORMON_DATA *CM, STEERING_LIST *SL, long coord, RUN *run, LINE_LIST *beamline, long find_only, long invert);
extern void compute_orbcor_matrices(CORMON_DATA *CM, STEERING_LIST *SL, long coord, RUN *run, LINE_LIST *beamline, long find_only, long invert, long fixed_length);

extern void setup_corrector_output(char *filename, RUN *run);
extern void dump_corrector_data(CORMON_DATA *CM, STEERING_LIST *SL, long index, char *plane, long step);
extern void setup_cormon_stats(char *filename, RUN *run);
extern void dump_cormon_stats(long verbose, long plane, double **kick, long n_kicks, 
    double **position, long n_positions, double *Cdp, long n_iterations, long cycle,
    long final_cycle, long step);
extern void setup_orb_traj_output(char *filename, char *mode, RUN *run);
extern void dump_orb_traj(TRAJECTORY *traj, long n_elems, char *description, long step);

extern void setup_correction_matrix_output(NAMELIST_TEXT *nltext, RUN *run, LINE_LIST *beamline, CORRECTION *correct,
                                    long *do_response);
extern void run_response_output(RUN *run, LINE_LIST *beamline, CORRECTION *correct, long tune_corrected);
extern void finish_response_output(void);

/* prototypes for counter.c: */
extern long advance_values1(double *value, long n_values, long *value_index, double *initial, double *step, 
                            double **enumerated_value, long *counter, long *max_count, long *flags, long n_indices);

/* prototypes for do_tracking13.c: */
extern double beta_from_delta(double p, double delta);
extern long do_tracking(double **coord, long *n_original, long *effort, LINE_LIST *beamline, double *P_central,    
    double **accepted, BEAM_SUMS **sums_vs_z, long *n_z_points, TRAJECTORY *traj_vs_z, RUN *run, long step,
    long flags, long n_passes);
extern void do_element_misalignment(ELEMENT_LIST *elem, double **coord, long n, long mode);
extern void offset_beam(double **coord, long n_to_track, MALIGN *offset, double P_central);
extern void do_match_energy(double **coord, long np, double *P_central, long change_beam);
extern void set_central_energy(double **coord, long np, double new_energy, double *P_central);
extern void set_central_momentum(double **coord, long np, double  P_new, double *P_central);
extern void center_beam(double **part, CENTER *center, long np);
void drift_beam(double **part, long np, double length, long order);
void scatter(double **part, long np, double Po, SCATTER *scatter);
void store_fitpoint_twiss_parameters(MARK *fpt, char *name, long occurence, TWISS *twiss);
void store_fitpoint_beam_parameters(MARK *fpt, char *name, long occurence, double **coord, long np, double Po);

extern void track_through_kicker(double **part, long np, KICKER *kicker, double p_central, long pass,
      long order);
extern long simple_rf_cavity(double **part, long np, RFCA *rfca, double **accepted, double *P_central,
                             double zEnd);
extern long modulated_rf_cavity(double **part, long np, MODRF *modrf, double P_central, double zEnd);
extern void set_up_kicker(KICKER *kicker);
void add_to_particle_energy(double *coord, double timeOfFlight, double Po, double dgamma);


#define FID_MODE_LIGHT   0x001UL
#define FID_MODE_TMEAN   0x002UL
#define FID_MODE_FIRST   0x004UL
#define FID_MODE_PMAX    0x008UL
double findFiducialTime(double **part, long np, double s0, double sOffset,
                        double p0, unsigned long mode);
extern unsigned long parseFiducialMode(char *mode);

/* prototypes for final_props.c */
extern void SDDS_FinalOutputSetup(SDDS_TABLE *SDDS_table, char *filename, long mode, long lines_per_row,
                           char *contents, char *command_file, char *lattice_file, 
                           char **varied_quantity_name, char **varied_quantity_unit, long varied_quantities,
                           char **error_element_name, char **error_element_unit, long error_elements,
                           char **optimization_quantity_name, char **optimization_quantity_unit, long optimization_quantities,
                           char *caller);
extern void dump_final_properties(SDDS_TABLE *SDDS_table, BEAM_SUMS *sums,
     double *varied_quan, char *first_varied_quan_name, long n_varied_quan,
     double *perturbed_quan, char *first_perturbed_quan_name, long n_perturbed_quan,
     double *optim_quan, char *first_optim_quan_name, long n_optim_quan,
     long step, double **particle, long n_original, double p_central, VMATRIX *M);
extern long compute_final_properties
    (double *data, BEAM_SUMS *sums, long n_original, double p_central, VMATRIX *M, double **coord, long step);
extern void rpn_store_final_properties(double *value, long number);
extern long get_final_property_index(char *name);
extern long count_final_properties();

extern double beam_width(double fraction, double **coord, long n_part, long sort_coord);
extern double rms_emittance(double **coord, long i1, long i2, long n);
extern double rms_longitudinal_emittance(double **coord, long n, double Po);
extern double rms_norm_emittance(double **coord, long i1, long i2, long ip, long n, double Po);
extern void compute_longitudinal_parameters(ONE_PLANE_PARAMETERS *bp, double **coord, long n, double Po);

/* prototypes for matrix_output.c: */
void simplify_units(char *buffer, char **numer, long n_numer, char **denom, long n_denom);
void run_matrix_output(RUN *run, LINE_LIST *beamline);
void setup_matrix_output(NAMELIST_TEXT *nltext, RUN *run, LINE_LIST *beamline);

/* prototypes for twiss.c: */
VMATRIX *compute_periodic_twiss(double *betax, double *alphax, double *etax, double *etaxp,
    double *phix, double *betay, double *alphay, double *etay, double *etayp, double *phiy,
    ELEMENT_LIST *elem, double *clorb, RUN *run);
void propagate_twiss_parameters(TWISS *twiss0, double *tune, ELEMENT_LIST *elem, long plane, RUN *run, double *traj);
long get_twiss_mode(long *mode, double *x_twiss, double *y_twiss);
void compute_twiss_parameters(RUN *run, LINE_LIST *beamline, double *starting_coord, long matched, 
    double beta_x, double alpha_x, double eta_x, double etap_x, 
    double beta_y, double alpha_y, double eta_y, double etap_y);
void dump_twiss_parameters(TWISS *twiss0, ELEMENT_LIST *elem, long n_elem, 
    double *tune, double *chromaticity, double *acceptance, double alphac,
    long final_values_only, long tune_corrected, RUN *run);
void setup_twiss_output(NAMELIST_TEXT *nltext, RUN *run, LINE_LIST *beamline, long *do_twiss_output);
void run_twiss_output(RUN *run, LINE_LIST *beamline, double *starting_coord, long tune_corrected);
void finish_twiss_output(void);
void copy_doubles(double *source, double *target, long n);

/* prototypes for elegant.c: */
extern char *compose_filename(char *template, char *root_name);
extern double find_beam_p_central(char *input);
void center_beam_on_coords(double **particle, long n_part, double *coord, long center_momentum_also);
void link_date(void);
void check_heap();
void do_print_dictionary(char *filename);
void print_dictionary_entry(FILE *fp, long type);

/* prototypes for error.c: */
extern void error_setup(ERROR *errcon, NAMELIST_TEXT *nltext, RUN *run_cond, LINE_LIST *beamline);
extern void add_error_element(ERROR *errcon, NAMELIST_TEXT *nltext, LINE_LIST *beamline);
extern double parameter_value(char *pname, long elem_type, long param, LINE_LIST *beamline);
extern double perturbation(double xamplitude, double xcutoff, long xerror_type);
 
/* prototypes for extend_list.c: */
extern void extend_line_list(LINE_LIST **lptr);
extern void extend_elem_list(ELEMENT_LIST **eptr);
 
/* prototypes for get_beamline5.c: */
extern void show_elem(ELEMENT_LIST *eptr, long type);
extern LINE_LIST *get_beamline(char *madfile, char *use_beamline, double p_central);
extern void show_elem(ELEMENT_LIST *eptr, long type);
extern void free_elements(ELEMENT_LIST *elemlist);
extern void free_beamlines(LINE_LIST *beamline);
extern void do_save_lattice(NAMELIST_TEXT *nl, RUN *run, LINE_LIST *beamline);
void print_with_continuation(FILE *fp, char *s, long endcol);
void change_defined_parameter_values(char **elem_name, long *param_number, long *type, double *value, long n_elems);
void change_defined_parameter(char *elem_name, long param_number, long type, double value, char *valueString);

/* prototypes for limit_amplitudes4.c: */
extern long rectangular_collimator(double **initial, RCOL *rcol, long np, double **accepted, double z, double P_central);
extern long limit_amplitudes(double **coord, double xmax, double ymax, long np, double **accepted, double z, double P_central,
                                   long extrapolate_z);
extern long elliptical_collimator(double **initial, ECOL *ecol, long np, double **accepted, double z, double P_central);
extern long elimit_amplitudes(double **coord, double xmax, double ymax, long np, double **accepted, double z,
    double P_central, long extrapolate_z);
extern long beam_scraper(double **initial, SCRAPER *scraper, long np, double **accepted, double z,
    double P_central);
 
/* prototypes for kick_sbend.c: */
long track_through_kick_sbend(double **part, long n_part, KSBEND *ksbend, double p_error, double Po,
    double **accepted, double z_start);
void bend_edge_kicks(double *x, double *xp, double *y, double *yp, double rho, double n, double beta, 
    double psi,  long which_edge);

/* prototypes for mad_parse4.c: */
extern long is_simple(char *s);
extern void fill_line(LINE_LIST *line, long nl, ELEMENT_LIST *elem, long ne, char *s);
extern ELEMENT_LIST *expand_line(ELEMENT_LIST *leptr, LINE_LIST *lptr,
    char *s, LINE_LIST *line, long nl, ELEMENT_LIST *elem, long ne, char *part_of);
extern long is_simple(char *s);
extern void fill_elem(ELEMENT_LIST *eptr, char *s, long type, FILE *fp_input);
extern long expand_phys(ELEMENT_LIST *leptr, char *entity, ELEMENT_LIST *elem_list,     
    long ne, LINE_LIST *line_list, long nl, long reverse, long multiplier, char *part_of);
extern void copy_element(ELEMENT_LIST *e1, ELEMENT_LIST *e2, long reverse);
void copy_named_element(ELEMENT_LIST *eptr, char *s, ELEMENT_LIST *elem);
extern long copy_line(ELEMENT_LIST *e1, ELEMENT_LIST *e2, long ne, long reverse, char *part_of);
extern long tell_type(char *s, ELEMENT_LIST *elem);
extern char *get_param_name(char *s);
extern char *find_param(char *s, char *param);
extern void unknown_parameter(char *parameter, char *element, char *type_name, char *caller);
extern void parse_element(char *p_elem, PARAMETER *parameter, long n_params,
    char *string, ELEMENT_LIST *eptr, char *type_name);
extern void parse_pepper_pot(PEPPOT *peppot, FILE *fp, char *name);
 
/* prototypes for malign_mat.c: */
extern void misalign_matrix(VMATRIX *M, double dx, double dy, double dz, double bend_angle);
extern VMATRIX *misalignment_matrix(MALIGN *malign, long order);
extern void offset_matrix(VMATRIX *M, double dx, double dxp, double dy, double dyp);

/* prototypes for matrix7.c: */
extern void print_matrices(FILE *fp, char *string, VMATRIX *M);
extern void initialize_matrices(VMATRIX *M, long order);
extern void null_matrices(VMATRIX *M);
extern void track_particles(double **final, VMATRIX *M, double  **initial, long n_part);
extern void free_matrices(VMATRIX *M);
extern void free_nonlinear_matrices(VMATRIX *M);
extern void set_matrix_pointers(double **C, double ***R, double ****T, double *****Q, VMATRIX *M);
extern long read_matrices(VMATRIX *M, FILE *fp);
extern void filter_matrices(VMATRIX *M, double threshold);
extern void random_matrices(VMATRIX *M, double C0, double R0, double T0, double Q0);
extern void copy_matrices(VMATRIX *M1, VMATRIX *M0);
extern long check_matrix(VMATRIX *M, char *comment);
 
/* prototypes for motion4.c: */
extern long motion(double **part, long n_part, void *field, long field_type, double P_central, double *dgamma,
    double *dP, double **accepted, double z_start);
 
/* prototypes for multipole.c: */
extern long multipole_tracking(double **particle, long n_part, MULT *multipole, double p_error, double Po, double **accepted, double z_start);
extern long multipole_tracking2(double **particle, long n_part, ELEMENT_LIST *elem, double p_error, double Po, double **accepted, double z_start);

/* prototypes for output_magnets.c: */
extern void output_magnets(char *filename, char *line_name, LINE_LIST *beamline);
 
/* prototypes for pepper_pot2.c: */
extern long pepper_pot_plate(double **initial, PEPPOT *peppot, long np, double **accepted);
 
/* prototypes for phase_reference.c: */
extern long get_phase_reference(double *phase, long phase_ref_number);
extern long set_phase_reference(long phase_ref_number, double phase);
extern void delete_phase_references(void);
extern long unused_phase_reference(void);
extern double get_reference_phase(long phase_ref, double phase0);

/* prototypes for print_line2.c: */
extern void print_line(FILE *fp, LINE_LIST *lptr);
extern void print_elem_list(FILE *fp, ELEMENT_LIST *eptr);
 
/* prototypes for quad_matrix3.c: */
extern VMATRIX *quadrupole_matrix(double K1, double l, long maximum_order, double tilt, double ffringe, double fse);
extern VMATRIX *quad_fringe(double l, double ko, long order, long reverse, double fse);
extern void qfringe_R_matrix(double *R11, double *R21, double *R12, double *R22, double dk_dz, double l);
extern void qfringe_T_matrix(double *T116, double *T126, double *T216, double *T226,
    double *T511, double *T512, double *T522, double dk_dz, double l, long reverse);
extern VMATRIX *qfringe_matrix(double K1, double l, double tilt, long direction, long order, double fse);

/* prototypes for tilt_matrices.c: */
extern void tilt_matrices0(VMATRIX *M, double tilt);
extern void tilt_matrices(VMATRIX *M, double tilt);
extern VMATRIX *rotation_matrix(double tilt);
extern void rotate_coordinates(double *coord, double angle);

/* prototypes for track_ramp.c: */
extern void track_through_ramped_deflector(double **final, RMDF *ramp_param, double **initial, long n_particles, double pc_central);
 
/* prototypes for track_rf2.c: */
extern void track_through_rf_deflector(double **final, RFDF *rf_param, double **initial, long n_particles, double pc_central);
 
/* prototypes for vary4.c: */
extern void vary_setup(VARY *_control, NAMELIST_TEXT *nltext, RUN *run, LINE_LIST *beamline);
extern void add_varied_element(VARY *_control, NAMELIST_TEXT *nltext, RUN *run, LINE_LIST *beamline);
extern long vary_beamline(VARY *_control, ERROR *errcon, RUN *run, LINE_LIST *beamline);
extern long perturb_beamline(VARY *_control, ERROR *errcon, RUN *run, LINE_LIST *beamline);
extern ELEMENT_LIST *find_element(char *elem_name,  ELEMENT_LIST **context, ELEMENT_LIST *elem);
extern ELEMENT_LIST *wfind_element(char *elem_name,  ELEMENT_LIST **context, ELEMENT_LIST *elem);
ELEMENT_LIST *find_element_index(char *elem_name,  ELEMENT_LIST **context,  ELEMENT_LIST *elem, long *index);
extern long confirm_parameter(char *item_name, long type);
extern void set_element_flags(LINE_LIST *beamline, char **elem_name, long *elem_perturb_flags,
    long *type, long *param, long n_elems,
    long pflag, long mflag, long overwrite, long permit_flags);
extern void assert_parameter_values(char **elem_name, long *param_number, long *type, double *value, long n_elems,
    LINE_LIST *beamline);
long get_parameter_value(double *value, char *elem_name, long param_number, long type, LINE_LIST *beamline);
extern void assert_perturbations(char **elem_name, long *param_number, long *type, long n_elems,
    double *amplitude, double *cutoff, long *error_type, double *perturb, long *elem_perturb_flags,
    long *bind_number, FILE *fp_log, long step, LINE_LIST *beamline, long permit_flags);
extern long compute_changed_matrices(LINE_LIST *beamline, RUN *run);

/* prototypes for routines in optimize.c */

void do_optimization_setup(OPTIMIZATION_DATA *_optimize, NAMELIST_TEXT *nltext, RUN *run, LINE_LIST *beamline);
void add_optimization_variable(OPTIMIZATION_DATA *_optimize, NAMELIST_TEXT *nltext, RUN *run, LINE_LIST *beamline);
void add_optimization_constraint(OPTIMIZATION_DATA *_optimize, NAMELIST_TEXT *nltext, RUN *run, LINE_LIST *beamline);
void summarize_optimization_setup(OPTIMIZATION_DATA *_optimize);
void do_optimize(NAMELIST_TEXT *nltext, RUN *run1, VARY *control1, ERROR *error1, LINE_LIST *beamline1, 
            BEAM *beam1, OUTPUT_FILES *output1, OPTIMIZATION_DATA *optimization_data1, long beam_type1);
void add_optimization_covariable(OPTIMIZATION_DATA *_optimize, NAMELIST_TEXT *nltext, RUN *run, LINE_LIST *beamline);

/* prototype for sample.c */
long sample_particles(double **initial, SAMPLE *samp, long np, double **accepted, double z, double p0);


/* prototype for run_rpnexpr.c */
void run_rpn_expression(NAMELIST_TEXT *nltext);

/* prototypes for trace.c */
void process_trace_request(NAMELIST_TEXT *nltext);
void log_entry(char *routine);
void log_exit(char *routine);

/* flag word for trace mode */
extern long trace_mode;
#define TRACE_ENTRY 1
#define TRACE_HEAP_VERIFY 2
#define TRACE_MEMORY_LEVEL 4
#define TRACEBACK_ON 8

/* global particle ID counter */
extern long particleID;

/* prototypes for chrom.c */
void setup_chromaticity_correction(NAMELIST_TEXT *nltext, RUN *run, LINE_LIST *beamline, CHROM_CORRECTION *chrom);
void do_chromaticity_correction(CHROM_CORRECTION *chrom, RUN *run, LINE_LIST *beamline, double *clorb,
        long step, long last_iteration);
void computeChromaticities(double *chromx, double *chromy, TWISS *twiss, VMATRIX *M);

/* prototypes for lorentz.c */
long  lorentz(double **part, long n_part, void *field, long field_type, double P_central, double **accepted);
void lorentz_report(void);

/* prototypes for kick_poly.c */
long polynomial_kicks(double **particle, long n_part, KPOLY *kpoly, double p_error, double Po,
    double **accepted, double z_start);

/* prototypes for ramp_p.c */
long ramp_momentum(double **coord, long np, RAMPP *rampp, double *P_central, long pass);

/* prototypes for ramped_rfca.c */
long ramped_rf_cavity(double **part, long np, RAMPRF *ramprf, double P_central, 
                      double L_central, double z_cavity, long pass);

/* prototypes for closed_orbit.c */
extern void dump_closed_orbit(TRAJECTORY *traj, long n_elems, long step);
void finish_clorb_output(void);
void run_closed_orbit(RUN *run, LINE_LIST *beamline, double *starting_coord, BEAM *beam, long do_output);
void setup_closed_orbit(NAMELIST_TEXT *nltext, RUN *run, LINE_LIST *beamline);

/* prototypes for aperture_search.c */
void setup_aperture_search(NAMELIST_TEXT *nltext, RUN *run, VARY *control);
long do_aperture_search(RUN *run, VARY *control, ERROR *errcon, LINE_LIST *beamline);
long do_aperture_search_mp(RUN *run, VARY *control, ERROR *errcon, LINE_LIST *beamline);
long do_aperture_search_sp(RUN *run, VARY *control, ERROR *errcon, LINE_LIST *beamline);
void finish_aperture_search(RUN *run, VARY *control, ERROR *errcon, LINE_LIST *beamline);

/* prototypes for analyze.c */
void setup_transport_analysis(NAMELIST_TEXT *nltext, RUN *run, VARY *control, ERROR *errcon);
void do_transport_analysis(RUN *run, VARY *control, ERROR *errcon, LINE_LIST *beamline, double *orbit);
void finish_transport_analysis(RUN *run, VARY *control, ERROR *errcon, LINE_LIST *beamline);

/* prototypes for tune.c */
void setup_tune_correction(NAMELIST_TEXT *nltext, RUN *run, LINE_LIST *beamline, TUNE_CORRECTION *tune);
void do_tune_correction(TUNE_CORRECTION *tune, RUN *run, LINE_LIST *beamline, double *clorb, long step, long last_iteration);

/* prototypes for link_elements.c */
void element_link_control(ELEMENT_LINKS *links, NAMELIST_TEXT *nltext, RUN *run_cond, LINE_LIST *beamline);
void add_element_links(ELEMENT_LINKS *links, NAMELIST_TEXT *nltext, LINE_LIST *beamline);
long assert_element_links(ELEMENT_LINKS *links, RUN *run_cond, LINE_LIST *beamline, long flags);
void reset_element_links(ELEMENT_LINKS *links, RUN *run_cond, LINE_LIST *beamline);

void compute_amplification_factors(NAMELIST_TEXT *nltext, RUN *run, CORRECTION *correct,
    long closed_orbit, LINE_LIST *beamline);

void track_through_matter(double **part, long np, MATTER *matter, double Po);

void track_through_rfmode(double **part, long np, RFMODE *rfmode, double Po,
    char *element_name, double element_z, long pass, long n_passes);
void set_up_rfmode(RFMODE *rfmode, char *element_name, double element_z, long n_passes, RUN *run, long n_particles,
                   double Po, double Lo);

void track_through_trfmode(double **part, long np, TRFMODE *trfmode, double Po,
    char *element_name, double element_z, long pass, long n_passes);
void set_up_trfmode(TRFMODE *rfmode, char *element_name, double element_z, long n_passes, RUN *run, long n_particles);

void track_through_zlongit(double **part, long np, ZLONGIT *zlongit, double Po, RUN *run, long i_pass);
void set_up_zlongit(ZLONGIT *zlongit, RUN *run);

void track_SReffects(double **coord, long n, SREFFECTS *SReffects, double Po, 
                             TWISS *twiss);

long track_through_csbend(double **part, long n_part, CSBEND *csbend, double p_error, double Po, double **accepted,
    double z_start);

void output_floor_coordinates(NAMELIST_TEXT *nltext, RUN *run, LINE_LIST *beamline);

void setup_load_parameters(NAMELIST_TEXT *nltext, RUN *run, LINE_LIST *beamline);
long do_load_parameters(LINE_LIST *beamline, long change_definitions);
#define NO_LOAD_PARAMETERS 0
#define PARAMETERS_LOADED 1
#define PARAMETERS_ENDED 2
void finish_load_parameters();

typedef struct {
    char *name, *text; 
    } SDDS_DEFINITION;
#define SDDS_EOS_NEWFILE 1
#define SDDS_EOS_COMPLETE 2
extern void SDDS_ElegantOutputSetup(SDDS_TABLE *SDDS_table, char *filename, long mode, long lines_per_row,
                             char *contents, char *command_file, char *lattice_file, SDDS_DEFINITION *parameter_definition,
                             long n_parameters, SDDS_DEFINITION *column_definition, long n_columns,
                             char *caller, long flags);
extern void SDDS_PhaseSpaceSetup(SDDS_TABLE *SDDS_table, char *filename, long mode, long lines_per_row, char *contents,
                          char *command_file, char *lattice_file, char *caller);
extern void SDDS_BeamLossSetup(SDDS_TABLE *SDDS_table, char *filename, long mode, long lines_per_row, char *contents, 
                          char *command_file, char *lattice_file, char *caller);
extern void SDDS_SigmaMatrixSetup(SDDS_TABLE *SDDS_table, char *filename, long mode, long lines_per_row,
                           char *command_file, char *lattice_file, char *caller);
extern void SDDS_WatchPointSetup(WATCH *waatch, long mode, long lines_per_row,
                          char *command_file, char *lattice_file, char *caller, char *qualifier);
extern void dump_watch_particles(WATCH *watch, long step, long pass, double **particle, long particles, double Po,
                                 double length);
extern void dump_watch_parameters(WATCH *watch, long step, long pass, long n_passes, double **particle, long particles, 
                           long original_particles,  double Po);
extern void dump_watch_FFT(WATCH *watch, long step, long pass, long n_passes, double **particle, long particles,
                           long original_particles,  double Po);
extern void do_watch_FFT(double **data, long n_data, long slot, long window_code);
extern void dump_lost_particles(SDDS_TABLE *SDDS_table, double **particle, long particles, long step);
extern void dump_centroid(SDDS_TABLE *SDDS_table, BEAM_SUMS *sums, LINE_LIST *beamline, long n_elements, long bunch,
                          double p_central);
