/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: final_props.c
 * purpose: routines used to compute and output values for
 *          'final' file.
 * M. Borland, 1993.
 */

#include "mdb.h"
#include "mdbsun.h"
#include "track.h"
#include "matlib.h"

#define ANALYSIS_BINS 10000

static double tmp_safe_sqrt;
#define SAFE_SQRT(x) ((tmp_safe_sqrt=(x))<0?0.0:sqrt(tmp_safe_sqrt))

#define FINAL_PROPERTY_PARAMETERS (96+9)
#define FINAL_PROPERTY_LONG_PARAMETERS 5
#define F_SIGMA_OFFSET 0
#define F_SIGMA_QUANS 7
#define F_CENTROID_OFFSET F_SIGMA_OFFSET+F_SIGMA_QUANS
#define F_CENTROID_QUANS 7
#define F_SIGMAT_OFFSET F_CENTROID_OFFSET+F_CENTROID_QUANS
#define F_SIGMAT_QUANS 15
#define F_T_OFFSET F_SIGMAT_OFFSET+F_SIGMAT_QUANS
#define F_T_QUANS 5
#define F_EMIT_OFFSET F_T_OFFSET+F_T_QUANS
#define F_EMIT_QUANS 5
#define F_NEMIT_OFFSET F_EMIT_OFFSET+F_EMIT_QUANS
#define F_NEMIT_QUANS 4
#define F_WIDTH_OFFSET F_NEMIT_OFFSET+F_NEMIT_QUANS
#define F_WIDTH_QUANS 16
#define F_PERC_OFFSET F_WIDTH_OFFSET+F_WIDTH_QUANS
#define F_PERC_QUANS 9
#define F_RMAT_OFFSET F_PERC_OFFSET+F_PERC_QUANS
#define F_RMAT_QUANS 31
#define F_STATS_OFFSET F_RMAT_OFFSET+F_RMAT_QUANS
#define F_STATS_QUANS 5
#define F_N_OFFSET F_STATS_OFFSET+F_STATS_QUANS
#define F_N_QUANS 1
#if (F_N_QUANS+F_N_OFFSET)!=FINAL_PROPERTY_PARAMETERS
#error "FINAL_PROPERTY_PARAMETERS is inconsistent with parameter offsets"
#endif

static SDDS_DEFINITION final_property_parameter[FINAL_PROPERTY_PARAMETERS] = {
/* beginning of type=double parameters */
    {"Sx",    "&parameter name=Sx, symbol=\"$gs$r$bx$n\", units=m, type=double &end"},
    {"Sxp",    "&parameter name=Sxp, symbol=\"$gs$r$bx'$n\", type=double &end"},
    {"Sy",    "&parameter name=Sy, symbol=\"$gs$r$by$n\", units=m, type=double &end"},
    {"Syp",    "&parameter name=Syp, symbol=\"$gs$r$by'$n\", type=double &end"},
    {"Ss",    "&parameter name=Ss, units=m, symbol=\"$gs$r$bs$n\", type=double &end"},
    {"Sdelta",    "&parameter name=Sdelta, symbol=\"$gs$bd$n$r\", type=double &end"},
    {"St", "&parameter name=St, symbol=\"$gs$r$bt$n\", type=double, units=s, &end"},
    {"Cx", "&parameter name=Cx, symbol=\"<x>\", units=m, type=double &end"},
    {"Cxp", "&parameter name=Cxp, symbol=\"<x'>\", type=double &end"},
    {"Cy", "&parameter name=Cy, symbol=\"<y>\", units=m, type=double &end"},
    {"Cyp", "&parameter name=Cyp, symbol=\"<y'>\", type=double &end"},
    {"Cs", "&parameter name=Cs, symbol=\"<s>\", units=m, type=double &end"},
    {"Cdelta", "&parameter name=Cdelta, symbol=\"<$gd$r>\", type=double &end"},
    {"Ct", "&parameter name=Ct, symbol=\"<t>\", units=s, type=double &end"},
    {"s12",    "&parameter name=s12, symbol=\"$gs$r$b12$n\", units=m, type=double &end"},
    {"s13",    "&parameter name=s13, symbol=\"$gs$r$b13$n\", units=\"m$a2$n\", type=double &end"},
    {"s14",    "&parameter name=s14, symbol=\"$gs$r$b14$n\", units=m, type=double &end"},
    {"s15",    "&parameter name=s15, symbol=\"$gs$r$b15$n\", units=\"m^a2$n\", type=double &end"},
    {"s16",    "&parameter name=s16, symbol=\"$gs$r$b16$n\", units=m, type=double &end"},
    {"s23",    "&parameter name=s23, symbol=\"$gs$r$b23$n\", units=m, type=double &end"},
    {"s24",    "&parameter name=s24, symbol=\"$gs$r$b24$n\", type=double &end"},
    {"s25",    "&parameter name=s25, symbol=\"$gs$r$b25$n\", units=m, type=double &end"},
    {"s26",    "&parameter name=s26, symbol=\"$gs$r$b26$n\", type=double &end"},
    {"s34",    "&parameter name=s34, symbol=\"$gs$r$b34$n\", units=m, type=double &end"},
    {"s35",    "&parameter name=s35, symbol=\"$gs$r$b35$n\", units=\"m^a2$n\", type=double &end"},
    {"s36",    "&parameter name=s36, symbol=\"$gs$r$b36$n\", units=m, type=double &end"},
    {"s45",    "&parameter name=s45, symbol=\"$gs$r$b45$n\", units=m, type=double &end"},
    {"s46",    "&parameter name=s46, symbol=\"$gs$r$b46$n\", type=double &end"},
    {"s56",    "&parameter name=s56, symbol=\"$gs$r$b56$n\", units=m, type=double &end"},
    {"Transmission", "&parameter name=Transmission, description=Transmission, type=double &end"},
    {"pCentral", "&parameter name=pCentral, symbol=\"p$bcen$n\", units=\"m$be$nc\", type=double &end"},
    {"pAverage", "&parameter name=pAverage, symbol=\"p$bave$n\", units=\"m$be$nc\", type=double &end"},
    {"KAverage",  "&parameter name=KAverage, symbol=\"K$bave$n\", units=MeV, type=double &end"},
    {"Charge", "&parameter name=Charge, units=C, type=double &end"},
    {"ex", "&parameter name=ex, symbol=\"$ge$r$bx$n\", units=m, type=double &end"},
    {"ey", "&parameter name=ey, symbol=\"$ge$r$by$n\", units=m, type=double &end"},
    {"ecx", "&parameter name=ecx, symbol=\"$ge$r$bx,c$n\", units=m, type=double &end"},
    {"ecy", "&parameter name=ecy, symbol=\"$ge$r$by,c$n\", units=m, type=double &end"},
    {"el", "&parameter name=el, symbol=\"$ge$r$bl$n\", units=s, type=double &end"},
    {"enx", "&parameter name=enx, symbol=\"$ge$r$bx,n$n\", type=double, units=m &end"},
    {"eny", "&parameter name=eny, symbol=\"$ge$r$by,n$n\", type=double, units=m &end"},
    {"ecnx", "&parameter name=ecnx, symbol=\"$ge$r$bx,cn$n\", type=double, units=m &end"},
    {"ecny", "&parameter name=ecny, symbol=\"$ge$r$by,cn$n\", type=double, units=m &end"},
    {"Wx", "&parameter name=Wx, type=double, units=m, symbol=\"W$bx$n\" &end"},
    {"Wy", "&parameter name=Wy, type=double, units=m, symbol=\"W$by$n\" &end"},
    {"Dt", "&parameter name=Dt, type=double, units=s, symbol=\"$gD$rt\" &end"},
    {"Ddelta", "&parameter name=Ddelta, type=double, symbol=\"$gDd$r\" &end"},
    {"Dt50", "&parameter name=Dt50, type=double, units=s, symbol=\"$gD$rt$b50$n\", &end"},
    {"Dt60", "&parameter name=Dt60, type=double, units=s, symbol=\"$gD$rt$b60$n\", &end"},
    {"Dt70", "&parameter name=Dt70, type=double, units=s, symbol=\"$gD$rt$b70$n\", &end"},
    {"Dt80", "&parameter name=Dt80, type=double, units=s, symbol=\"$gD$rt$b80$n\", &end"},
    {"Dt90", "&parameter name=Dt90, type=double, units=s, symbol=\"$gD$rt$b90$n\", &end"},
    {"Dt95", "&parameter name=Dt95, type=double, units=s, symbol=\"$gD$rt$b95$n\", &end"},
    {"Ddelta50", "&parameter name=Ddelta50, type=double, symbol=\"$gDd$r$b50$n\", &end"},
    {"Ddelta60", "&parameter name=Ddelta60, type=double, symbol=\"$gDd$r$b60$n\", &end"},
    {"Ddelta70", "&parameter name=Ddelta70, type=double, symbol=\"$gDd$r$b70$n\", &end"},
    {"Ddelta80", "&parameter name=Ddelta80, type=double, symbol=\"$gDd$r$b80$n\", &end"},
    {"Ddelta90", "&parameter name=Ddelta90, type=double, symbol=\"$gDd$r$b90$n\", &end"},
    {"Ddelta95", "&parameter name=Ddelta95, type=double, symbol=\"$gDd$r$b95$n\", &end"},
    {"DtPerc10", "&parameter name=DtPerc10, type=double, units=s, &end"},
    {"DtPerc20", "&parameter name=DtPerc20, type=double, units=s, &end"},
    {"DtPerc30", "&parameter name=DtPerc30, type=double, units=s, &end"},
    {"DtPerc40", "&parameter name=DtPerc40, type=double, units=s, &end"},
    {"DtPerc50", "&parameter name=DtPerc50, type=double, units=s, &end"},
    {"DtPerc60", "&parameter name=DtPerc60, type=double, units=s, &end"},
    {"DtPerc70", "&parameter name=DtPerc70, type=double, units=s, &end"},
    {"DtPerc80", "&parameter name=DtPerc80, type=double, units=s, &end"},
    {"DtPerc90", "&parameter name=DtPerc90, type=double, units=s, &end"},
    {"R11", "&parameter name=R11, type=double, symbol=\"R$b11$n\" &end"},
    {"R12", "&parameter name=R12, type=double, units=m, symbol=\"R$b12$n\" &end"},
    {"R13", "&parameter name=R13, type=double, symbol=\"R$b13$n\" &end"},
    {"R14", "&parameter name=R14, type=double, units=m, symbol=\"R$b14$n\" &end"},
    {"R15", "&parameter name=R15, type=double, units=m, symbol=\"R$b15$n\" &end"},
    {"R16", "&parameter name=R16, type=double, units=m, symbol=\"R$b16$n\" &end"},
    {"R21", "&parameter name=R21, type=double, units=m, symbol=\"R$b21$n\" &end"},
    {"R22", "&parameter name=R22, type=double, symbol=\"R$b22$n\" &end"},
    {"R23", "&parameter name=R23, type=double, units=1/m, symbol=\"R$b23$n\" &end"},
    {"R24", "&parameter name=R24, type=double, symbol=\"R$b24$n\" &end"},
    {"R25", "&parameter name=R25, type=double, symbol=\"R$b25$n\" &end"},
    {"R26", "&parameter name=R26, type=double, symbol=\"R$b26$n\" &end"},
    {"R31", "&parameter name=R31, type=double, symbol=\"R$b31$n\" &end"},
    {"R32", "&parameter name=R32, type=double, units=m, symbol=\"R$b32$n\" &end"},
    {"R33", "&parameter name=R33, type=double, symbol=\"R$b33$n\" &end"},
    {"R34", "&parameter name=R34, type=double, units=m, symbol=\"R$b34$n\" &end"},
    {"R35", "&parameter name=R35, type=double, units=m, symbol=\"R$b35$n\" &end"},
    {"R36", "&parameter name=R36, type=double, units=m, symbol=\"R$b36$n\" &end"},
    {"R41", "&parameter name=R41, type=double, units=1/m, symbol=\"R$b41$n\" &end"},
    {"R42", "&parameter name=R42, type=double, symbol=\"R$b42$n\" &end"},
    {"R43", "&parameter name=R43, type=double, units=1/m, symbol=\"R$b43$n\" &end"},
    {"R44", "&parameter name=R44, type=double, symbol=\"R$b44$n\" &end"},
    {"R45", "&parameter name=R45, type=double, symbol=\"R$b45$n\" &end"},
    {"R46", "&parameter name=R46, type=double, symbol=\"R$b46$n\" &end"},
    {"R51", "&parameter name=R51, type=double, symbol=\"R$b51$n\" &end"},
    {"R52", "&parameter name=R52, type=double, units=m, symbol=\"R$b52$n\" &end"},
    {"R53", "&parameter name=R53, type=double, symbol=\"R$b53$n\" &end"},
    {"R54", "&parameter name=R54, type=double, units=m, symbol=\"R$b54$n\" &end"},
    {"R55", "&parameter name=R55, type=double, symbol=\"R$b55$n\" &end"},
    {"R56", "&parameter name=R56, type=double, units=m, symbol=\"R$b56$n\" &end"},
    {"detR", "&parameter name=detR, type=double, symbol=\"det R\" &end"},
    {"CPU", "&parameter name=CPU, type=double, units=s &end"},
/* beginning of type=long parameters */
    {"MEM", "&parameter name=MEM, type=long, units=pages &end"},
    {"PF", "&parameter name=PF, type=long, units=pages &end"},
    {"Step",  "&parameter name=Step, type=long &end"},
    {"Steps",  "&parameter name=Steps, type=long &end"},
    {"Particles", "&parameter name=Particles, description=\"Number of particles\", type=long &end"},
    } ;


void SDDS_FinalOutputSetup(SDDS_TABLE *SDDS_table, char *filename, long mode, long lines_per_row,
                           char *contents, char *command_file, char *lattice_file, 
                           char **varied_quantity_name, char **varied_quantity_unit, long varied_quantities,
                           char **error_element_name, char **error_element_unit, long error_elements,
                           long *error_element_index, long *error_element_duplicates,
                           char **optimization_quantity_name, char **optimization_quantity_unit, 
                           long optimization_quantities,
                           char *caller)
{
  long i, duplicates=0;
  SDDS_ElegantOutputSetup(SDDS_table, filename, mode, lines_per_row, contents, command_file, 
                          lattice_file, final_property_parameter, FINAL_PROPERTY_PARAMETERS, NULL, 0,
                          caller, SDDS_EOS_NEWFILE);
  if ((varied_quantities &&
       !SDDS_DefineSimpleParameters(SDDS_table, varied_quantities, varied_quantity_name, varied_quantity_unit, SDDS_DOUBLE)) ||
      (optimization_quantities &&
       !SDDS_DefineSimpleParameters(SDDS_table, optimization_quantities, optimization_quantity_name, optimization_quantity_unit,
                                    SDDS_DOUBLE))) {
    fprintf(stdout, "Problem defining extra SDDS parameters in file %s (%s)\n", filename, caller);
    fflush(stdout);
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    exit(1);
  }
  for (i=0; i<error_elements; i++) {
    if (!SDDS_DefineSimpleParameter(SDDS_table, error_element_name[i], 
                                    error_element_unit[i], SDDS_DOUBLE)) {
      /* determine if the error is just because of a duplication */
      long j;
      for (j=0; j<i; j++)
        if (strcmp(error_element_name[i], error_element_name[j])==0)
          break;
      if (i==j) {
        fprintf(stdout, "Problem defining extra SDDS parameter %s in file %s (%s)\n", error_element_name[i],
                filename, caller);
        fflush(stdout);
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
        exit(1);
      }
      duplicates++;
    }
    if ((error_element_index[i] = SDDS_GetParameterIndex(SDDS_table, error_element_name[i]))<0) {
      fprintf(stdout, "Problem defining extra SDDS parameter %s in file %s (%s): couldn't retrieve index\n",
               error_element_name[i],
              filename, caller);
      fflush(stdout);
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      exit(1);
    }
  }
  *error_element_duplicates = duplicates;
  
  if (!SDDS_WriteLayout(SDDS_table)) {
    fprintf(stdout, "Unable to write SDDS layout for file %s (%s)\n", filename, caller);
    fflush(stdout);
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    exit(1);
  }
}

void dump_final_properties
    (SDDS_TABLE *SDDS_table, BEAM_SUMS *sums,
     double *varied_quan, char *first_varied_quan_name, long n_varied_quan,
     long totalSteps,
     double *perturbed_quan, long *perturbed_quan_index, 
     long perturbed_quan_duplicates, long n_perturbed_quan,
     double *optim_quan, char *first_optim_quan_name, long n_optim_quan,
     long step, double **particle, long n_original, double p_central, VMATRIX *M,
     double charge)
{
    long n_computed, n_properties;
    double *computed_properties;
    long index, i;

    log_entry("dump_final_properites");

    if (!SDDS_table)
        bomb("NULL SDDS_TABLE pointer (dump_final_properties)", NULL);
    if (!sums)
        bomb("NULL beam sums pointer (dump_final_properites)", NULL);
    if (n_varied_quan && (!varied_quan || !first_varied_quan_name))
        bomb("Unexpected NULL pointer for varied quanitity values/names (dump_final_properties)", NULL);
    if (n_perturbed_quan && (!perturbed_quan || !perturbed_quan_index))
        bomb("Unexpected NULL pointer for perturbed quantity values/names (dump_final_properties)", NULL);
    if (n_optim_quan && (!optim_quan || !first_optim_quan_name)) 
        bomb("Unexpected NULL pointer for optimization quantity values/names (dump_final_properties)", NULL);
    if (!M)
        bomb("NULL matrix pointer (dump_final_properties)", NULL);
    if (!particle)
        bomb("NULL particle coordinates pointer (dump_final_properties)", NULL);

    if (!SDDS_StartTable(SDDS_table, 0)) {
        SDDS_SetError("Problem starting SDDS table (dump_final_properties)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }

    if ((n_properties=SDDS_ParameterCount(SDDS_table)) !=
        (FINAL_PROPERTY_PARAMETERS+n_varied_quan+n_perturbed_quan+n_optim_quan-perturbed_quan_duplicates)) {
        fprintf(stdout, "error: the number of parameters (%ld) defined for the SDDS table for the final properties file is not equal to the number of quantities (%ld) for which information is provided (dump_final_properties)\n",
                n_properties, 
                FINAL_PROPERTY_PARAMETERS+n_varied_quan+n_perturbed_quan+n_optim_quan-perturbed_quan_duplicates);
        fflush(stdout);
        abort();
        }
    computed_properties = tmalloc(sizeof(*computed_properties)*n_properties);

    if ((n_computed=compute_final_properties
                       (computed_properties, sums, n_original, p_central, M, particle, step,
                        totalSteps, charge))!=
        (n_properties-(n_varied_quan+n_perturbed_quan+n_optim_quan-perturbed_quan_duplicates))) {
        fprintf(stdout, "error: compute_final_properties computed %ld quantities--%ld expected. (dump_final_properties)",
            n_computed, n_properties-(n_varied_quan+n_perturbed_quan+n_optim_quan-perturbed_quan_duplicates));
        fflush(stdout);
        abort();
        }

    if ((index=SDDS_GetParameterIndex(SDDS_table, "MEM"))<0) {
        SDDS_SetError("Problem getting SDDS index of Step parameter (dump_final_properties)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
    for (i=0; i<FINAL_PROPERTY_LONG_PARAMETERS; i++)
        if (!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE,
                                i+index, (long)computed_properties[i+index], -1)) {
            SDDS_SetError("Problem setting SDDS parameter values (dump_final_properties)");
            SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
            }
    if ((index=SDDS_GetParameterIndex(SDDS_table, "Sx"))<0) {
        SDDS_SetError("Problem getting SDDS index of Sx parameter (dump_final_properties)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
    for (i=0; i<FINAL_PROPERTY_PARAMETERS-FINAL_PROPERTY_LONG_PARAMETERS; i++)
        if (!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, 
                               i+index, computed_properties[i+index], -1)) {
            SDDS_SetError("Problem setting SDDS parameter values for computed properties (dump_final_properties)");
            SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
            }
    free(computed_properties);
    
    if (first_varied_quan_name) {
        if ((index=SDDS_GetParameterIndex(SDDS_table, first_varied_quan_name))<0) {
            SDDS_SetError("Problem getting SDDS index of first varied quantity parameter (dump_final_properties)");
            SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
            }
        for (i=0; i<n_varied_quan; i++)
            if (!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, 
                                   i+index, varied_quan[i], -1)) {
                SDDS_SetError("Problem setting SDDS parameter values for varied quantities (dump_final_properties)");
                SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
                }
        }

    if (perturbed_quan_index) {
      for (i=0; i<n_perturbed_quan; i++)
        if (!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE,
                                perturbed_quan_index[i], (double)0.0, -1)) {
          SDDS_SetError("Problem setting SDDS parameter values for perturbed quantities (dump_final_properties)");
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
      for (i=0; i<n_perturbed_quan; i++) {
        double value;
        if (!SDDS_GetParameterByIndex(SDDS_table, perturbed_quan_index[i], &value) ||
            !SDDS_SetParameters(SDDS_table, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE,
                                perturbed_quan_index[i], perturbed_quan[i]+value, -1)) {
          SDDS_SetError("Problem setting SDDS parameter values for perturbed quantities (dump_final_properties)");
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
      }
    }

    if (first_optim_quan_name) {
        if ((index=SDDS_GetParameterIndex(SDDS_table, first_optim_quan_name))<0) {
            SDDS_SetError("Problem getting SDDS index of first optimization quantity parameter (dump_final_properties)");
            SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
            }
        for (i=0; i<n_optim_quan; i++)
            if (!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE,
                                   i+index, optim_quan[i], -1)) {
                SDDS_SetError("Problem setting SDDS parameter values for optimization quantities (dump_final_properties)");
                SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
                }
        }
    
    if (!SDDS_WriteTable(SDDS_table)) {
        SDDS_SetError("Problem writing SDDS data for final properties (dump_final_properties)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
    if (!SDDS_DoFSync(SDDS_table))
      fprintf(stdout, "Warning: problem fsync'ing final properties output file\n");

    log_exit("dump_final_properties");
    }

/* routine: compute_final_properties
 * purpose: compute sigmas, centroids, number of particles, transmission,
 *          plus Rij (i<=4, j<=6) at end of beamline
 *
 */

long compute_final_properties
  (double *data, BEAM_SUMS *sums, long n_original, double p_central, VMATRIX *M, double **coord, 
   long step, long steps, double charge)
{
  register long i, j;
  long i_data, index, offset;
  double dp_min, dp_max, Ddp;
  double p_sum, gamma_sum, p, sum, tc, tmin, tmax, dt, t, pAverage;
  double **R, centroid[6];
  MATRIX Rmat;
  static double *tData = NULL, *deltaData = NULL;
  static long percDataMax = 0;
  double percLevel[12] = {25, 20, 15, 10, 5, 2.5, 75, 80, 85, 90, 95, 97.5};
  double tPosition[12], deltaPosition[12];
  double percLevel2[9] = {10,20,30,40,50,60,70,80,90};
  double tPosition2[9];
  log_entry("compute_final_properties");

  if (!data)
    bomb("return data array is null (compute_final_properties)", NULL);
  if (!sums)
    bomb("beam sums element is null (compute_final_properties)", NULL);
  if (!M || !M->C || !M->R)
    bomb("invalid/null transport map (compute_final_properties)", NULL);
  if (!coord)
    bomb("particle coordinate array is null (compute_final_properties)", NULL);
  
  /* compute centroids and sigmas */
  if (sums->n_part) {
    for (i=0; i<6; i++) 
      centroid[i] = data[i+F_CENTROID_OFFSET] = sums->centroid[i];
    for (i=0; i<6; i++)
      data[i+F_SIGMA_OFFSET   ] = sqrt(sums->sigma[i][i]);
    offset = F_SIGMAT_OFFSET;
    index = 0;
    /* sigma matrix elements sij */
    for (i=0; i<6; i++) {
      /* skip the diagonal element */
      for (j=i+1; j<6; j++) 
        data[offset++] = sums->sigma[i][j];
    }
    /* time centroid, sigma, and delta */
    tmax = dp_max = -(tmin=dp_min=DBL_MAX);
    if ((!tData || sums->n_part>percDataMax) &&
        (!(tData = malloc(sizeof(*tData)*(percDataMax=sums->n_part))) ||
         !(deltaData = malloc(sizeof(*deltaData)*(percDataMax=sums->n_part)))))
      bomb("memory allocation failure (compute_final_properties)", NULL);
    for (i=sum=0; i<sums->n_part; i++) {
      if (!coord[i]) {
        fprintf(stdout, "coordinate element for particle %ld is null (compute_final_properties)\n", i);
        fflush(stdout);
        abort();
      }
      if (coord[i][5]>dp_max)
        dp_max = coord[i][5];
      if (coord[i][5]<dp_min)
        dp_min = coord[i][5];
      p = p_central*(1+coord[i][5]);
      deltaData[i] = coord[i][5];
      sum += ( t = tData[i] = coord[i][4]/(p/sqrt(sqr(p)+1)*c_mks) );
      if (t<tmin)
        tmin = t;
      if (t>tmax)
        tmax = t;
    }
    dt = tmax-tmin;
    Ddp = dp_max - dp_min;
    data[6+F_CENTROID_OFFSET] = (tc = sum/sums->n_part);
    for (i=sum=0; i<sums->n_part; i++)
      sum += sqr( tData[i] - tc);
    data[6+F_SIGMA_OFFSET] = sqrt(sum/sums->n_part);
    /* results of these calls used below */
    approximate_percentiles(tPosition, percLevel, 12, tData, sums->n_part, ANALYSIS_BINS);
    approximate_percentiles(tPosition2, percLevel2, 9, tData, sums->n_part, ANALYSIS_BINS);
    approximate_percentiles(deltaPosition, percLevel, 12, deltaData, sums->n_part, ANALYSIS_BINS);
  }
  else {
    for (i=0; i<7; i++) 
      data[i+F_CENTROID_OFFSET] = data[i+F_SIGMA_OFFSET] = 0;
    dt = Ddp = tmin = tmax = 0;
    for (i=0; i<12; i++)
      tPosition[i] = deltaPosition[i] = 0;
  }

  /* transmission */
  if (n_original)
    data[F_T_OFFSET] = ((double)sums->n_part)/n_original;
  else
    data[F_T_OFFSET] = 0;
  /* lattice momentum */
  data[F_T_OFFSET+1] = p_central;

  /* compute average momentum and kinetic energy */
  p_sum = gamma_sum = 0;
  pAverage = p_central;
  for (i=0; i<sums->n_part; i++) {
    p_sum     += (p = (1+coord[i][5])*p_central);
    gamma_sum += sqrt(sqr(p)+1);
  }
  if (sums->n_part) {
    pAverage = data[F_T_OFFSET+2] = p_sum/sums->n_part;
    data[F_T_OFFSET+3] = (gamma_sum/sums->n_part-1)*me_mev;
  }
  else
    data[F_T_OFFSET+2] = data[F_T_OFFSET+3] = 0;
  /* beam charge */
  data[F_T_OFFSET+4] = charge;
  
  /* compute "sigma" from width of particle distributions for x and y */
  if (coord && sums->n_part>3) {
    
    data[F_WIDTH_OFFSET] = approximateBeamWidth(0.6826F, coord, sums->n_part, 0L)/2.;
    data[F_WIDTH_OFFSET+1] = approximateBeamWidth(0.6826F, coord, sums->n_part, 2L)/2.;
    data[F_WIDTH_OFFSET+2] = dt;
    data[F_WIDTH_OFFSET+3] = Ddp;
    for (i=0; i<6; i++) {
      data[F_WIDTH_OFFSET+ 4+i] =     tPosition[6+i]-    tPosition[0+i];
      data[F_WIDTH_OFFSET+10+i] = deltaPosition[6+i]-deltaPosition[0+i];
    }
    for (i=0; i<9; i++)
      data[F_PERC_OFFSET+i] = tPosition2[i]-tmin;
  }
  else {
    for (i=0; i<10; i++)
      data[F_WIDTH_OFFSET] = 0;
    for (i=0; i<9; i++)
      data[F_PERC_OFFSET+i] = tPosition2[i]-tmin;
  }

  /* compute emittances */
  data[F_EMIT_OFFSET]   = rms_emittance(coord, 0, 1, sums->n_part, NULL, NULL, NULL);
  data[F_EMIT_OFFSET+1] = rms_emittance(coord, 2, 3, sums->n_part, NULL, NULL, NULL);
  /* corrected transverse emittances */
  if (sums->sigma[5][5]) {
    data[F_EMIT_OFFSET+2] = SAFE_SQRT(sqr(data[F_EMIT_OFFSET]) - 
                                      (sqr(sums->sigma[0][5])*sums->sigma[1][1] -
                                       2*sums->sigma[0][1]*sums->sigma[0][5]*sums->sigma[1][5] +
                                       sqr(sums->sigma[1][5])*sums->sigma[0][0])/sums->sigma[5][5]);
    data[F_EMIT_OFFSET+3] = SAFE_SQRT(sqr(data[F_EMIT_OFFSET+1]) - 
                                      (sqr(sums->sigma[2][5])*sums->sigma[3][3] -
                                       2*sums->sigma[2][3]*sums->sigma[2][5]*sums->sigma[3][5] +
                                       sqr(sums->sigma[3][5])*sums->sigma[2][2])/sums->sigma[5][5]);
  } else 
    data[F_EMIT_OFFSET+2] = data[F_EMIT_OFFSET+3] = 0;

  data[F_EMIT_OFFSET+4] = rms_longitudinal_emittance(coord, sums->n_part, p_central);
  
  /* compute normalized emittances */
  for (i=0; i<4; i++)
    data[F_NEMIT_OFFSET+i]   = pAverage*data[F_EMIT_OFFSET+i];

  R = M->R;
  i_data = F_RMAT_OFFSET;
  for (i=0; i<5; i++) 
    for (j=0; j<6; j++)
      data[i_data++] = R[i][j];
  Rmat.a = R;
  Rmat.n = Rmat.m = 6;
  data[i_data] = m_det(&Rmat);
  
  /* run time statistics */
  i_data = F_STATS_OFFSET;
#if defined(UNIX) || defined(VAX_VMS)
  data[i_data++] = cpu_time()/100.0;
  data[i_data++] = memory_count();
  data[i_data++] = page_faults();
#else
  data[i_data++] = 0;
  data[i_data++] = 0;
  data[i_data++] = 0;
#endif
  data[i_data++] = step;
  data[i_data++] = steps;
  
  /* number of particles */
  data[i_data=F_N_OFFSET] = sums->n_part;

  log_exit("compute_final_properties");
  return(i_data+1);
}

double beam_width(double fraction, double **coord, long n_part, 
    long sort_coord)
{
  static long i_median, i_lo, i_hi;
  static double dx0, dx1;

  log_entry("beam_width");

  /* sort data in ascending order by the coordinate for which the
   * beam width is to be determined */
  set_up_row_sort(sort_coord, 6L, sizeof(**coord), double_cmpasc);
  qsort((void*)coord, n_part, sizeof(*coord), row_compare);
  for (i_median=1; i_median<n_part; i_median++) {
    if (coord[i_median][sort_coord]<coord[i_median-1][sort_coord])
      bomb("sort failure in beam_width()", NULL);
  }

  /* find indices of particles that are at +/- fraction/2 from the median */
  i_median = n_part/2;
  if ((i_lo = i_median - fraction/2.*n_part)<0) {
    fprintf(stdout, "warning: i_lo < 0 in beam_width\ni_median = %ld, n_part = %ld, fraction = %e\n",
            i_lo, n_part, fraction);
    fflush(stdout);
  }
  if ((i_hi = i_median + fraction/2.*n_part)>=n_part) {
    fprintf(stdout, "warning: i_hi >= n_part in beam_width!\ni_median = %ld, n_part = %ld, fraction = %e\n",
            i_hi, n_part, fraction);
    fflush(stdout);
  }

  if (i_lo!=0 && i_hi!=(n_part-1)) {
    dx0 = coord[i_hi][sort_coord]  -coord[i_lo][sort_coord];
    dx1 = coord[i_hi+1][sort_coord]-coord[i_lo-1][sort_coord];
    log_exit("beam_width");    
    return(INTERPOLATE(dx0, dx1, ((double)i_hi-i_lo+1)/n_part, ((double)i_hi-i_lo+3)/n_part, fraction));
  }
  else {
    log_exit("beam_width");    
    return(coord[i_hi][sort_coord]-coord[i_lo][sort_coord]);
  }
}

double rms_emittance(double **coord, long i1, long i2, long n, 
                     double *s11Return, double *s12Return, double *s22Return)
{
  double s11, s12, s22, x, xp;
  double xc, xpc;
  long i;
  
  if (!n)
    return(0.0);

  /* compute centroids */
  for (i=xc=xpc=0; i<n; i++) {
    xc  += coord[i][i1];
    xpc += coord[i][i2];
  }
  xc  /= n;
  xpc /= n;

  for (i=s11=s12=s22=0; i<n; i++) {
    s11 += sqr(x  = coord[i][i1]-xc );
    s22 += sqr(xp = coord[i][i2]-xpc);
    s12 += x*xp;
  }
  if (s11Return)
    *s11Return = s11/n;
  if (s22Return)
    *s22Return = s22/n;
  if (s12Return)
    *s12Return = s12/n;
  return(SAFE_SQRT(s11*s22-sqr(s12))/n);
}

double rms_norm_emittance(double **coord, long i1, long i2, long ip, long n, double Po)
{
    return Po*rms_emittance(coord, i1, i2, n, NULL, NULL, NULL);
#if 0
    double s11, s12, s22;
    double x, px, xc, pxc;
    long i;
    
    if (!n)
        return(0.0);

    /* compute centroids */
    for (i=xc=pxc=0; i<n; i++) {
        xc  += coord[i][i1];
        pxc += coord[i][i2]*(1+coord[i][ip])/sqrt(1+sqr(coord[i][1])+sqr(coord[i][3]));
        }
    xc  /= n;
    pxc /= n;

    for (i=s11=s12=s22=0; i<n; i++) {
        s11 += sqr(  x = coord[i][i1] - xc );
        s22 += sqr( px = coord[i][i2]*(1+coord[i][ip])/sqrt(1+sqr(coord[i][1])+sqr(coord[i][3])) - pxc );
        s12 += x*px;
        }

    return(Po*SAFE_SQRT(s11*s22-sqr(s12))/n);
#endif
  }

double rms_longitudinal_emittance(double **coord, long n, double Po)
{
    double s11, s12, s22, dt, ddp;
    double tc, dpc, beta, P;
    long i;
    static double *time = NULL;
    static long max_n = 0;

    if (!n)
        return(0.0);

    if (n>max_n)
        time = trealloc(time, sizeof(*time)*(max_n=n));

    log_entry("rms_logitudinal_emittance");

    /* compute centroids */
    for (i=tc=dpc=0; i<n; i++) {
        P = Po*(1+coord[i][5]);
        beta = P/sqrt(P*P+1);
        time[i] = coord[i][4]/(beta*c_mks);
        tc  += time[i];
        dpc += coord[i][5];
        }
    tc  /= n;
    dpc /= n;

    for (i=s11=s12=s22=0; i<n; i++) {
        s11 += sqr(dt  =  time[i]    - tc);
        s22 += sqr(ddp = coord[i][5] - dpc);
        s12 += dt*ddp;
        }

    log_exit("rms_longitudinal_emittance");
    return(SAFE_SQRT(s11*s22-sqr(s12))/n);
    }

void compute_longitudinal_parameters(ONE_PLANE_PARAMETERS *bp, double **coord, long n, double Po)
{
    double s11, s12, s22, S1, S2, dt, ddp;
    double tc, dpc, beta, P;
    double tmin, tmax, dpmin, dpmax;
    long i;
    static double *time = NULL;
    static long max_n = 0;

    log_entry("compute_longitudinal_parameters");

    if (!n)
        return;
    if (!bp) {
        fprintf(stdout, "NULL ONE_PLANE_PARAMETERS pointer passed to compute_longitudinal_parameters\n");
        fflush(stdout);
        abort();
        }

    if (n>max_n)
        time = trealloc(time, sizeof(*time)*(max_n=n));

    /* compute centroids */
    tmax  = -(tmin  = DBL_MAX);
    dpmax = -(dpmin = DBL_MAX);
    for (i=tc=dpc=0; i<n; i++) {
        P = Po*(1+coord[i][5]);
        beta = P/sqrt(P*P+1);
        time[i] = coord[i][4]/(beta*c_mks);
        tc  += time[i];
        dpc += coord[i][5];
        if (coord[i][5]>dpmax)
            dpmax = coord[i][5];
        if (coord[i][5]<dpmin)
            dpmin = coord[i][5];
        if (time[i]>tmax)
            tmax = time[i];
        if (time[i]<tmin)
            tmin = time[i];
        }
    bp->c1 = (tc  /= n);
    bp->c2 = (dpc /= n);
    bp->min1 = tmin;
    bp->min2 = dpmin;
    bp->max1 = tmax;
    bp->max2 = dpmax;

    for (i=S1=S2=s11=s12=s22=0; i<n; i++) {
        s11 += sqr(dt  =  time[i]    - tc);
        s22 += sqr(ddp = coord[i][5] - dpc);
        s12 += dt*ddp;
        S1 += sqr(time[i]);
        S2 += sqr(coord[i][5]);
        }
    bp->s11 = (s11 /= n);
    bp->s12 = (s12 /= n);
    bp->s22 = (s22 /= n);
    bp->S1  = sqrt(S1/n);
    bp->S2  = sqrt(S2/n);
    bp->emittance = SAFE_SQRT(s11*s22-sqr(s12));

    log_exit("compute_longitudinal_properties");
    }

void compute_transverse_parameters(ONE_PLANE_PARAMETERS *bp, double **coord, long n, long plane)
{
  double s11, s12, s22, S1, S2, dx, dxp;
  double xc, xpc;
  double xmin, xmax, xpmin, xpmax;
  long i, offset;

  if (!n)
    return;
  if (!bp) {
    fprintf(stdout, "NULL ONE_PLANE_PARAMETERS pointer passed to compute_transverse_parameters\n");
    fflush(stdout);
    abort();
  }
  offset = plane?2:0;
  
  /* compute centroids */
  xmax  = -(xmin  = DBL_MAX);
  xpmax = -(xpmin = DBL_MAX);
  for (i=xc=xpc=0; i<n; i++) {
    xc  += coord[i][offset+0];
    xpc += coord[i][offset+1];
    if (coord[i][offset]>xmax)
      xmax = coord[i][offset];
    if (coord[i][offset]<xmin)
      xmin = coord[i][offset];
    if (coord[i][offset+1]>xpmax)
      xpmax = coord[i][offset+1];
    if (coord[i][offset+1]<xpmin)
      xpmin = coord[i][offset+1];
  }
  
  bp->c1 = (xc  /= n);
  bp->c2 = (xpc /= n);
  bp->min1 = xmin;
  bp->min2 = xpmin;
  bp->max1 = xmax;
  bp->max2 = xpmax;

  for (i=S1=S2=s11=s12=s22=0; i<n; i++) {
    s11 += sqr(dx = coord[i][offset] - xc);
    s22 += sqr(dxp = coord[i][offset+1] - xpc);
    s12 += dx*dxp;
    S1 += sqr(coord[i][offset]);
    S2 += sqr(coord[i][offset+1]);
  }
  bp->s11 = (s11 /= n);
  bp->s12 = (s12 /= n);
  bp->s22 = (s22 /= n);
  bp->S1  = sqrt(S1/n);
  bp->S2  = sqrt(S2/n);
  bp->emittance = SAFE_SQRT(s11*s22-sqr(s12));
}

void rpn_store_final_properties(double *value, long number)
{
    static long *memory_number = NULL;
    long i;
    log_entry("rpn_store_final_parameters");
    if (number!=FINAL_PROPERTY_PARAMETERS) {
        fprintf(stdout, "error: number of values (%ld) being stored != FINAL_PROPERTY_PARAMETERS (rpn_store_final_parameters)\n",
                number);
        fflush(stdout);
        abort();
        }
    if (!memory_number) {
        memory_number = tmalloc(sizeof(*memory_number)*number);
        for (i=0; i<number; i++) 
            memory_number[i] = rpn_create_mem(final_property_parameter[i].name);
        }
    for (i=0; i<number; i++)
        rpn_store(value[i], memory_number[i]);
    log_exit("rpn_store_final_parameters");
    }

long get_final_property_index(char *name)
{
    long i;
    for (i=0; i<FINAL_PROPERTY_PARAMETERS; i++)
        if (strcmp(name, final_property_parameter[i].name)==0)
            return(i);
    return(-1);
    }

long count_final_properties()
{
    return(FINAL_PROPERTY_PARAMETERS);
    }

double approximateBeamWidth(double fraction, double **part, long nPart, long iCoord)
{
  double *hist, *cdf;
  long maxBins=ANALYSIS_BINS, bins=ANALYSIS_BINS, i50, iLo, iHi, i;
  double xMin, xMax, dx;
  xMin = xMax = dx = 0;
  
  /* make histogram of the coordinate */
  hist = tmalloc(sizeof(*hist)*bins);
  binParticleCoordinate(&hist, &maxBins, &xMin, &xMax, &dx, &bins, 
                        1.01, part, nPart, iCoord);

  /* sum histogram to get CDF */
  cdf = hist;
  for (i=1; i<bins; i++)
    cdf[i] += cdf[i-1];
  /* normalize CDF and find 50% point */
  i50 = bins/2;
  for (i=0; i<bins; i++) {
    cdf[i] /= cdf[bins-1];
    if (cdf[i]<0.50)
      i50 = i;
  }
  /* find locations containing half the indicated area around the 50% point */
  iLo = iHi = i50;
  for (i=i50; i<bins; i++) {
    if ((cdf[i]-0.5)<fraction/2)
      iHi = i;
    else 
      break;
  }
  for (i=i50; i>=0; i--) {
    if ((0.5-cdf[i])<fraction/2)
      iLo = i;
    else break;
  }
  free(hist);
  return (iHi-iLo)*dx;
}

