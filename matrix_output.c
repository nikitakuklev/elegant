/* Copyright 1994 by Michael Borland and Argonne National Laboratory,
 * all rights reserved.
 */
/* file: matrix_output()
 * purpose: handle user requests for matrix output
 *
 * Michael Borland, 1991
 */
#include "mdb.h"
#include "track.h"
#include "match_string.h"

static long output_now = 0;
static long n_outputs = 0;
static FILE **fp_printout;      /* for printout */
static long *print_full_only = NULL;
static long *SDDS_order = NULL;
static long *print_order = NULL;
static char *unit[6] = {"m", "rad", "m", "rad", "m", "1"};
static char **start_name = NULL;
static long *start_occurence = NULL;
static char **SDDS_match = NULL;

static SDDS_TABLE *SDDS_matrix = NULL;
static long *SDDS_matrix_initialized = NULL;
static long *SDDS_matrix_count = NULL;

#define IC_S 0
#define IC_ELEMENT 1
#define IC_OCCURENCE 2
#define IC_TYPE 3
#define N_COLUMNS 4
static SDDS_DEFINITION column_definition[N_COLUMNS] = {
    {"s", "&column name=s, units=m, type=double, description=\"Distance\" &end"},
    {"ElementName", "&column name=ElementName, type=string, description=\"Element name\", format_string=%10s &end"},
    {"ElementOccurence", 
         "&column name=ElementOccurence, type=long, description=\"Occurence of element\", format_string=%6ld &end"},
    {"ElementType", "&column name=ElementType, type=string, description=\"Element-type name\", format_string=%10s &end"},
    } ;

#define IP_STEP 0
#define N_PARAMETERS 1
static SDDS_DEFINITION parameter_definition[N_PARAMETERS] = {
    {"Step", "&parameter name=Step, type=long, description=\"Simulation step\" &end"},
    } ;


void SDDS_set_matrices(SDDS_TABLE *SDDS_table, VMATRIX *M, long order,
    ELEMENT_LIST *elem, long i_element, long n_elements);

void setup_matrix_output(
    NAMELIST_TEXT *nltext,
    RUN *run, 
    LINE_LIST *beamline
    )
{
#include "matrix_output.h"
    long i, j, k, l;
    char *denom[4], *numer[4];
    long n_denom, n_numer;
    char s[100], t[100];
    char buffer[SDDS_MAXLINE];

    log_entry("setup_matrix_output");

    /* process namelist input */
    set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
    set_print_namelist_flags(0);
    process_namelist(&matrix_output, nltext);
    print_namelist(stderr, &matrix_output);

    /* check for validity of namelist inputs */
    if (printout==NULL && SDDS_output==NULL)
        bomb("no output filenames given in namelist matrix_output", NULL);
    if (printout && !(printout_order>0 && printout_order<4))
        bomb("printout_order is invalid", NULL);
    if (SDDS_output && !(SDDS_output_order>0 && SDDS_output_order<4))
        bomb("SDDS_output_order is invalid", NULL);
    printout   = compose_filename(printout, run->rootname);
    SDDS_output = compose_filename(SDDS_output, run->rootname);

    fp_printout = trealloc(fp_printout, sizeof(*fp_printout)*(n_outputs+1));
    print_full_only = trealloc(print_full_only, sizeof(*print_full_only)*(n_outputs+1));
    print_order= trealloc(print_order, sizeof(*print_order)*(n_outputs+1));
    start_name = trealloc(start_name, sizeof(*start_name)*(n_outputs+1));
    start_occurence = trealloc(start_occurence, sizeof(*start_occurence)*(n_outputs+1));
    SDDS_match = trealloc(SDDS_match, sizeof(*SDDS_match)*(n_outputs+1));
    SDDS_order= trealloc(SDDS_order, sizeof(*SDDS_order)*(n_outputs+1));
    SDDS_matrix= trealloc(SDDS_matrix, sizeof(*SDDS_matrix)*(n_outputs+1));
    SDDS_matrix_initialized= trealloc(SDDS_matrix_initialized, sizeof(*SDDS_matrix_initialized)*(n_outputs+1));
    SDDS_matrix_count= trealloc(SDDS_matrix_count, sizeof(*SDDS_matrix_count)*(n_outputs+1));

    if (start_from)
        cp_str(start_name+n_outputs, start_from);
    else 
        start_name[n_outputs] = NULL;
    start_occurence[n_outputs] = start_from_occurence;
    print_order[n_outputs] = printout?printout_order:0;
    print_full_only[n_outputs] = full_matrix_only;

    SDDS_order[n_outputs]   = SDDS_output?SDDS_output_order:0;
    if (SDDS_output_match)
        cp_str(SDDS_match+n_outputs, SDDS_output_match);
    else
        SDDS_match[n_outputs] = NULL;
    SDDS_matrix_initialized[n_outputs] = SDDS_matrix_count[n_outputs] = 0;
    SDDS_ZeroMemory(SDDS_matrix+n_outputs, sizeof(*SDDS_matrix));

    if (printout) {
        fp_printout[n_outputs] = fopen_e(printout, "w", 0);
        }
    else
        fp_printout[n_outputs] = NULL;

    if (SDDS_output) {
        SDDS_ElegantOutputSetup(SDDS_matrix+n_outputs, SDDS_output, SDDS_BINARY, 1, "matrix", 
                                run->runfile, run->lattice, parameter_definition, N_PARAMETERS,
                                column_definition, N_COLUMNS, "setup_matrix_output", 
                                SDDS_EOS_NEWFILE);

        for (i=0; i<6; i++) {
            sprintf(buffer, "&column name=C%ld, symbol=\"C$b%ld$n\", type=double ", i+1, i+1);
                    if (SDDS_StringIsBlank(unit[i]))
                        strcpy(t, " &end");
                    else
                        sprintf(t, "units=%s &end", unit[i]);
                    strcat(buffer, t);
                    if (!SDDS_ProcessColumnString(SDDS_matrix+n_outputs, buffer, 0)) {
                        SDDS_SetError("Problem defining SDDS matrix output Rij columns (setup_matrix_output)");
                        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
                        }

            }
        if (SDDS_output_order>=1) {
            for (i=0; i<6; i++) {
                for (j=0; j<6; j++) {
                    sprintf(buffer, "&column name=R%ld%ld, symbol=\"R$b%ld%ld$n\", type=double ", i+1, j+1, i+1, j+1);
                    if (i==j)
                        strcpy(t, " &end");
                    else
                        sprintf(t, "units=%s/%s &end", unit[i], unit[j]);
                    strcat(buffer, t);
                    if (!SDDS_ProcessColumnString(SDDS_matrix+n_outputs, buffer, 0)) {
                        SDDS_SetError("Problem defining SDDS matrix output Rij columns (setup_matrix_output)");
                        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
                        }
                    }
                }
            }
        if (SDDS_output_order>=2) {
            n_numer = 1;
            n_denom = 2;
            for (i=0; i<6; i++) {
                for (j=0; j<6; j++) {
                    numer[0] = unit[i];
                    denom[0] = unit[j];
                    for (k=0; k<=j; k++) {
                        sprintf(buffer, "&column name=T%ld%ld%ld, symbol=\"T$b%ld%ld%ld$n\", type=double ",
                                i+1, j+1, k+1, i+1, j+1, k+1);
                        numer[0] = unit[i];
                        denom[0] = unit[j];
                        denom[1] = unit[k];
                        simplify_units(s, numer, n_numer, denom, n_denom);
                        if (SDDS_StringIsBlank(s))
                            sprintf(t, " &end");
                        else
                            sprintf(t, "units=%s &end", s);
                        strcat(buffer, t);
                        if (!SDDS_ProcessColumnString(SDDS_matrix+n_outputs, buffer, 0)) {
                            SDDS_SetError("Problem defining SDDS matrix output Tijk columns (setup_matrix_output)");
                            SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
                            }
                        }
                    }
                }
            }

        if (SDDS_output_order>=3) {
            n_numer = 1;
            n_denom = 3;
            for (i=0; i<6; i++) {
                for (j=0; j<6; j++) {
                    for (k=0; k<=j; k++) {
                        for (l=0; l<=k; l++) {
                            sprintf(buffer, "&column name=U%ld%ld%ld%ld, symbol=\"U$b%ld%ld%ld%ld$n\", type=double ",
                                    i+1, j+1, k+1, l+1, i+1, j+1, k+1, l+1);
                            numer[0] = unit[i];
                            denom[0] = unit[j];
                            denom[1] = unit[k];
                            denom[2] = unit[l];
                            simplify_units(t, numer, n_numer, denom, n_denom);
                            if (SDDS_StringIsBlank(s))
                                sprintf(t, " &end");
                            else
                                sprintf(t, "units=%s &end", s);
                            strcat(buffer, t);
                            if (!SDDS_ProcessColumnString(SDDS_matrix+n_outputs, buffer, 0)) {
                                SDDS_SetError("Problem defining SDDS matrix output Uijk columns (setup_matrix_output)");
                                SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
                                }
                            }
                        }
                    }
                }
            }

        if (!SDDS_WriteLayout(SDDS_matrix+n_outputs)) {
            SDDS_SetError("Unable to write SDDS layout for matrix output");
            SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
            }
        SDDS_matrix_initialized[n_outputs] = 1;
        }
    else
        SDDS_matrix_initialized[n_outputs] = 0;

    if (!output_at_each_step) {
        /* user wants output now */
        output_now = n_outputs;
        n_outputs++;
        run_matrix_output(run, beamline);
        n_outputs--;
        if (SDDS_matrix_initialized[n_outputs]) {
            if (!SDDS_Terminate(SDDS_matrix+n_outputs)) {
                SDDS_SetError("Problem terminating SDDS output (setup_matrix_output)");
                SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
                }
            SDDS_matrix_initialized[n_outputs] = 0;
            }
        if (fp_printout[n_outputs])
            fclose(fp_printout[n_outputs]);
        fp_printout[n_outputs] = NULL;
        }
    else {
        /* otherwise, the information is saved, so advance array counter */
        n_outputs++;
        }
    output_now = 0;

    log_exit("setup_matrix_output");
    }            

void run_matrix_output(
    RUN *run, 
    LINE_LIST *beamline
    )
{
    ELEMENT_LIST *member, *first_member;
    long i, n_elem_no_matrix, n_elements, sfo;
    long i_SDDS_output, n_SDDS_output;
    long i_output, output_order;
    VMATRIX *M1, *M2, *tmp;
    char s[256];
    double z0;

    if (n_outputs==0)
        return;

    log_entry("run_matrix_output");

    calculate_matrices(beamline, run);

    for (i_output=output_now; i_output<n_outputs; i_output++) {
        output_order= MAX(print_order[i_output], SDDS_order[i_output]);
        initialize_matrices(M1=tmalloc(sizeof(*M1)), output_order);
        initialize_matrices(M2=tmalloc(sizeof(*M2)), output_order);
        for (i=0; i<6; i++)
            M1->R[i][i] = 1;
        
        n_elements = n_elem_no_matrix = 0;
        first_member = member = &(beamline->elem);
        z0 = 0;
        sfo = start_occurence[i_output];
        if (start_name[i_output]!=NULL) {
            member = NULL;
            while (sfo--)
                if (!(first_member=find_element(start_name[i_output], &member, &(beamline->elem))))
                    bomb("can't find specified occurence of given element for matrix output", NULL);
            fprintf(stderr, "starting matrix output from element %s at z=%e m\n", first_member->name, first_member->end_pos);
            member = first_member;
            if (first_member->pred)
                z0 = first_member->pred->end_pos;
            }
        beamline->elem_recirc = NULL;
        beamline->i_recirc    = 0;
        while (member) {
#if DEBUG
            fprintf(stderr, "working on matrix for %s\n", member->name);
            print_elem(stderr, member);
#endif
            if (!member->matrix)
                compute_matrix(member, run, NULL);
            n_elements++;
            member = member->succ;
            }
        if (SDDS_matrix_initialized[i_output]) {
            member = first_member;
            n_SDDS_output = 1;        /* for the beginning element */
            while (member) {
                if (entity_description[member->type].flags&HAS_MATRIX) {
                    if (!SDDS_match[i_output] || wild_match(member->name, SDDS_match[i_output]))
                        n_SDDS_output++;
                    }
                member = member->succ;
                }
            }
        if (fp_printout[i_output])
            print_line(fp_printout[i_output], beamline);
        if (SDDS_matrix_initialized[i_output]) {
            ELEMENT_LIST start_elem;
            start_elem.end_pos = 0;
            start_elem.name = "_BEG_";
            start_elem.type = T_MARK;
            start_elem.occurence = 1;
            SDDS_set_matrices(SDDS_matrix+i_output, M1, SDDS_order[i_output], &start_elem, 0, n_SDDS_output);
            }
        i_SDDS_output = 1;
        member = first_member;
        while (member) {
            if (entity_description[member->type].flags&HAS_MATRIX) {
                if (fp_printout[i_output] && !print_full_only[i_output]) {
#ifdef DEBUG
                    fprintf(stderr, "doing printout for matrix of %s %s\n", 
                           entity_name[member->type], member->name);
#endif
                    if (member->matrix->order>print_order[i_output]) {
                        SWAP_LONG(member->matrix->order, print_order[i_output]);
                        print_matrices(fp_printout[i_output], member->name, member->matrix);
                        SWAP_LONG(member->matrix->order, print_order[i_output]);
                        }
                    else
                        print_matrices(fp_printout[i_output], member->name, member->matrix);
                    }
#ifdef DEBUG
                fprintf(stderr, "concatenating matrix of %s %s\n", 
                       entity_name[member->type], member->name);
#endif
                concat_matrices(M2, member->matrix, M1);
                tmp = M2;
                M2  = M1;
                M1  = tmp;
                if (fp_printout[i_output] && !print_full_only[i_output]) {
                  if (M1->order > print_order[i_output]) {
                    SWAP_LONG(M1->order, print_order[i_output]);
                    print_matrices(fp_printout[i_output], "Concantenated matrix after last element", M1);
                    SWAP_LONG(M1->order, print_order[i_output]);
                  } else {
                    print_matrices(fp_printout[i_output], "Concantenated matrix after last element", M1);
                  }
                }
                if (SDDS_matrix_initialized[i_output] && 
                    (!SDDS_match[i_output] || wild_match(member->name, SDDS_match[i_output]))) {
                    SDDS_set_matrices(SDDS_matrix+i_output, M1, SDDS_order[i_output], member,
                                    i_SDDS_output++, n_SDDS_output);
                    }
                }
            else {
                n_elem_no_matrix++;
                }
            member = member->succ;
            }
        if (SDDS_matrix_initialized[i_output]) {
            if (!SDDS_WriteTable(SDDS_matrix+i_output)) {
                SDDS_SetError("Unable to write matrix data (run_matrix_output)");
                SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
                }
            if (!SDDS_EraseData(SDDS_matrix+i_output)) {
                SDDS_SetError("Unable to erase matrix data (run_matrix_output)");
                SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
                }
            if (n_elem_no_matrix)
                fprintf(stderr, "warning: %ld elements had no matrix\n", n_elem_no_matrix);
            }
        if (fp_printout[i_output]) {
            if (n_elem_no_matrix)
                sprintf(s, "full matrix---WARNING: %ld ELEMENTS HAD NO MATRIX",
                        n_elem_no_matrix);
            else
                sprintf(s, "full matrix");
            SWAP_LONG(M1->order, print_order[i_output]);
            print_matrices(fp_printout[i_output], s, M1);
            SWAP_LONG(M1->order, print_order[i_output]);
            }
      }
    log_exit("run_matrix_output");
    }


void simplify_units(char *buffer, char **numer, long n_numer, char **denom, long n_denom)
{
    long in, id, i, j;
    long *mult_numer, *mult_denom;
    char t[256];

    log_entry("simplify_units");

#ifdef DEBUG
    fprintf(stderr, "*** simplify_units called\nnumerators: ");
    for (in=0; in<n_numer; in++)
        fprintf(stderr, "%s ", numer[in]);
    fprintf(stderr, "\ndenominators: ");
    for (id=0; id<n_denom; id++)
        fprintf(stderr, "%s ", denom[id]);
    fputc('\n', stderr);
#endif

    for (in=0; n_numer>0 && in<n_numer; in++) {
        for (id=0; n_denom>0 && id<n_denom && n_numer>0 && in<n_numer; id++) {
            if (strcmp(numer[in], denom[id])==0) {
                for (i=in; i<n_numer-1; i++)
                    numer[i] = numer[i+1];
                for (i=id; i<n_denom-1; i++)
                    denom[i] = denom[i+1];
                n_numer--; in--;
                n_denom--; id--;
                }
            }
        }
    for (id=0; n_denom>0 && id<n_denom; id++) {
        if (strcmp(denom[id], "1")==0) {
            for (i=id; i<n_denom-1; i++)
                denom[i] = denom[i+1];
            n_denom--; id--;
            }
        }

#ifdef DEBUG
    fprintf(stderr, "common factors divided out\nnumerators: ");
    for (in=0; in<n_numer; in++)
        fprintf(stderr, "%s ", numer[in]);
    fprintf(stderr, "\ndenominators: ");
    for (id=0; id<n_denom; id++)
        fprintf(stderr, "%s ", denom[id]);
    fputc('\n', stderr);
#endif

    mult_numer = tmalloc(sizeof(*mult_numer)*(n_numer?n_numer:1));
    for (in=0; in<n_numer; in++) {
        mult_numer[in] = 1;
        for (i=in+1; i<n_numer; i++) {
            if (strcmp(numer[in], numer[i])==0) {
                mult_numer[in]++;
                for (j=i; j<n_numer-1; j++)
                    numer[j] = numer[j+1];
                n_numer--;
                i--;
                }
            }
        }

    mult_denom = tmalloc(sizeof(*mult_denom)*(n_denom?n_denom:1));
    for (id=0; id<n_denom; id++) {
        mult_denom[id] = 1;
        for (i=id+1; i<n_denom; i++) {
            if (strcmp(denom[id], denom[i])==0) {
                mult_denom[id]++;
                for (j=i; j<n_denom-1; j++)
                    denom[j] = denom[j+1];
                n_denom--;
                i--;
                }
            }
        }

#ifdef DEBUG
    fprintf(stderr, "multiplicities established\nnumerators: ");
    for (in=0; in<n_numer; in++)
        fprintf(stderr, "%s^%ld ", numer[in], mult_numer[in]);
    fprintf(stderr, "\ndenominators: ");
    for (id=0; id<n_denom; id++)
        fprintf(stderr, "%s^%ld ", denom[id], mult_denom[id]);
    fputc('\n', stderr);
#endif

    if (n_numer==0 && n_denom==0)
        strcpy(buffer, "");
    else {
        buffer[0] = 0;
        if (n_numer<1)
            strcpy(buffer, "1");
        else {
            for (in=0; in<n_numer; in++) {
                if (mult_numer[in]==1)
                    sprintf(t, "%s%s", numer[in], (in==n_numer-1?"":"-"));
                else
                    sprintf(t, "%s$a%ld$n%s", numer[in], mult_numer[in], (in==n_numer-1?"":"-"));
                strcat(buffer, t);
                }
            }
        if (n_denom==1)
            strcat(buffer, "/");
        else if (n_denom>1)
            strcat(buffer, "/(");
        if (n_denom>=1) {
            for (id=0; id<n_denom; id++) {
                if (mult_denom[id]==1)
                    sprintf(t, "%s%s", denom[id], (id==n_denom-1?"":"-"));
                else
                    sprintf(t, "%s$a%ld$n%s", denom[id], mult_denom[id], (id==n_denom-1?"":"-"));
                strcat(buffer, t);
                }
            if (n_denom>1)
                strcat(buffer, ")");
            }
        }
    tfree(mult_numer); mult_numer = NULL;
    tfree(mult_denom); mult_denom = NULL;
    log_exit("simplify_units");
    }        
       
void finish_matrix_output()
{
    long i_output;
    for (i_output=0; i_output<n_outputs; i_output++) {
        if (fp_printout[i_output])
            fclose(fp_printout[i_output]);
        if (SDDS_matrix_initialized[i_output] && !SDDS_Terminate(SDDS_matrix+i_output)) {
            SDDS_SetError("Problem terminating SDDS matrix output (finish_matrix_output)");
            SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
            }
        SDDS_matrix_initialized[i_output] = 0;
        fp_printout[i_output] = NULL;
        }
    n_outputs = 0;
    }

void SDDS_set_matrices(SDDS_TABLE *SDDS_table, VMATRIX *M, long order,
    ELEMENT_LIST *elem, long i_element, long n_elements
    )
{
    register long i, j, k, l, index;

    log_entry("SDDS_set_matrices");

    if (!M || !(M->C) || !(M->R) || (order>1 && !(M->T)) || (order>2 && !(M->Q)))
        bomb("bad matrix passed to SDDS_set_matrices()", NULL);

    if (i_element==0) {
        if (!SDDS_StartTable(SDDS_table, n_elements)) {
            SDDS_SetError("Problem starting SDDS table for matrix output (SDDS_set_matrices)");
            SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
            }
        }
            
    if (!SDDS_SetRowValues(SDDS_table, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, i_element, 
                           IC_S, elem->end_pos, IC_ELEMENT, elem->name, IC_OCCURENCE, elem->occurence,
                           IC_TYPE, entity_name[elem->type], -1)) {
        SDDS_SetError("Problem setting row values for SDDS matrix output (SDDS_set_matrices)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }

    index = IC_TYPE+1;
    for (i=0; i<6; i++, index)
        if (!SDDS_SetRowValues(SDDS_table, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, i_element, 
                               index++, M->C[i], -1)) {
            SDDS_SetError("Problem setting row values for SDDS matrix output (SDDS_set_matrices)");
            SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
            }

    for (i=0; i<6; i++)
        for (j=0; j<6; j++)
            if (!SDDS_SetRowValues(SDDS_table, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, i_element, 
                                   index++, M->R[i][j], -1)) {
                SDDS_SetError("Problem setting row values for SDDS matrix output (SDDS_set_matrices)");
                SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
                }
 
    if (order>1)
        for (i=0; i<6; i++)
            for (j=0; j<6; j++)
                for (k=0; k<=j; k++)
                    if (!SDDS_SetRowValues(SDDS_table, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, i_element,
                                           index++, M->T[i][j][k], -1)) {
                        SDDS_SetError("Problem setting row values for SDDS matrix output (SDDS_set_matrices)");
                        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
                        }

    if (order>2)
        for (i=0; i<6; i++)
            for (j=0; j<6; j++)
                for (k=0; k<=j; k++)
                    for (l=0; l<=k; l++)
                        if (!SDDS_SetRowValues(SDDS_table, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, i_element,
                                               index++, M->Q[i][j][k][l], -1)) {
                            SDDS_SetError("Problem setting row values for SDDS matrix output (SDDS_set_matrices)");
                            SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
                            }


    log_exit("SDDS_set_matrices");
    }

