/* Copyright 1994 by Michael Borland and Argonne National Laboratory,
 * all rights reserved.
 */
/* routine: get_phase_reference()
 * purpose: store/retrieve reference phases for time-dependent elements
 *
 * Michael Borland, 1989
 */
#include "track.h"
#include "mdb.h"

/* flag bit values */
#define FL_REF_PHASE_SET 1

struct phase_reference {
    long ref_number;
    long flags;
    double phase;
    } *reference = NULL;
static long n_references=0;


long get_phase_reference(
    double *phase,
    long phase_ref_number
    )
{
    long i;
#ifdef DEBUG    
    printf("get_phase_reference(%le, %ld)\n", *phase, phase_ref_number);
#endif

    log_entry("get_phase_reference");

    if (phase_ref_number==0) {
        log_exit("get_phase_reference");
        return(0);
        }

    for (i=0; i<n_references; i++) {
        if (reference[i].ref_number==phase_ref_number) {
            if (reference[i].flags&FL_REF_PHASE_SET) {
                *phase = reference[i].phase;
#ifdef DEBUG
                printf("returning REF_PHASE_RETURNED\n");
#endif
                log_exit("get_phase_reference");
                return(REF_PHASE_RETURNED);
                }
#ifdef DEBUG
                printf("returning REF_PHASE_NOT_SET\n");
#endif
            log_exit("get_phase_reference");
            return(REF_PHASE_NOT_SET);
            }
        }
#ifdef DEBUG
                printf("returning REF_PHASE_NONEXISTENT\n");
#endif
    log_exit("get_phase_reference");
    return(REF_PHASE_NONEXISTENT);
    }

long set_phase_reference(
    long phase_ref_number,  /* number of the phase reference group */
    double phase            /* phase to assert for fiducial particle */
    )
{
    long i;

    log_entry("set_phase_reference");
    
#ifdef DEBUG
    printf("set_phase_reference(%ld, %le)\n", phase_ref_number, phase);
#endif

    if (phase_ref_number==0) {
        log_exit("set_phase_reference");
        return(0);
        }

    for (i=0; i<n_references; i++) {
        if (reference[i].ref_number==phase_ref_number) {
            reference[i].phase = phase;
            reference[i].flags = FL_REF_PHASE_SET;
#ifdef DEBUG
            printf("existing phase reference set\n");
#endif
            log_exit("set_phase_reference");
            return(1);
            }
        }
    if (phase_ref_number>LONG_MAX/2)
        bomb("please use a small integer for the phase_reference number", NULL);
    reference = trealloc(reference, sizeof(*reference)*(++n_references));
    reference[i].ref_number = phase_ref_number;
    reference[i].phase = phase;
    reference[i].flags = FL_REF_PHASE_SET;
#ifdef DEBUG
    printf("new phase reference set\n");
#endif
    log_exit("set_phase_reference");
    return(1);
    }

void delete_phase_references(void)
{
    long i;

    log_entry("delete_phase_references");

#ifdef DEBUG
    printf("phase references deleted\n");
#endif
    for (i=0; i<n_references; i++)
        reference[i].flags = 0;
    log_exit("delete_phase_references");
    }

long unused_phase_reference()
{
    static long big=LONG_MAX;

    log_entry("unused_phase_reference");

    reference = trealloc(reference, sizeof(*reference)*(n_references+1));
    reference[n_references].ref_number = big;
    reference[n_references].phase = 0;
    reference[n_references].flags = 0;
    n_references++;
    log_exit("unused_phase_reference");
    return(big--);
    }

double get_reference_phase(long phase_ref, double phase0)
    /* routine to maintain emulate old get_reference_phase(),
     * which is obsolete and should be phased out (pun intended)
     */
{
    double phase;

    log_entry("get_reference_phase");
#ifdef DEBUG
    printf("obsolete routine get_reference_phase called\n");
#endif
    switch (get_phase_reference(&phase, phase_ref)) {
        case REF_PHASE_RETURNED:
            log_exit("get_reference_phase");
            return(phase);
            break;
        case REF_PHASE_NOT_SET:
        case REF_PHASE_NONEXISTENT:
        default:
            set_phase_reference(phase_ref, phase0);
            log_exit("get_reference_phase");
            return(phase0);
            break;
        }
    }

