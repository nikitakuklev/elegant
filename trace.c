/* Copyright 1994 by Michael Borland and Argonne National Laboratory,
 * all rights reserved.
 */
/* file: trace.c
 * purpose: routines for tracing program calls
 * method: a file is maintained that contains the calling sequence.  It is
 *         updated by file positioning using a stack of file positions.
 *         In addition, a list of strings is maintained so that a trace-back
 *         can be given in event of a crash.
 *
 * M.Borland, 1992 
 */
#include "mdb.h"
#include "mdbsun.h"
#include "track.h"
#include "trace.h"
#include <signal.h>

#define MAX_LENGTH 78

static long trace_level = 0;
static long max_level = 0;
static char **routine_name = NULL;
static FILE *fp = NULL;
static FILE *fpmem = NULL;
static long *file_pos = NULL;
static char blank[MAX_LENGTH];
static long in_trace_routine = 0;
static long memory_level = 0;

void process_trace_request(NAMELIST_TEXT *nltext)
{
    long i;

    if (fp)
        fclose(fp);
    fp = NULL;

    /* process the namelist text */
    set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
    set_print_namelist_flags(0);
    process_namelist(&trace, nltext);
    print_namelist(stderr, &trace);

    if (record_allocation && filename)
        keep_alloc_record(filename);

#if defined(VAX_VMS)
    fputs("warning: program trace is not supported on this system", stderr);
    return;
#endif

    trace_mode = 0;

    if (!filename && trace_on)
        bomb("provide filename for program trace", NULL);
    if (memory_log)
        fpmem = fopen_e(memory_log, "w", 0);
    else
        fpmem = NULL;
    trace_mode += (memory_log?TRACE_MEMORY_LEVEL:0);
    trace_mode += (traceback_on?TRACEBACK_ON:0);
    if (trace_on) {
        fp = fopen_e(filename, "w", FOPEN_SAVE_IF_EXISTS);
        for (i=0; i<MAX_LENGTH; i++)
            blank[i] = ' ';
        trace_mode += TRACE_ENTRY + (heap_verify_depth?TRACE_HEAP_VERIFY:0);
        printf("trace activated into file %s\n", filename);
        }
    }

void log_entry(char *routine)
{
    long len;

#if defined(VAX_VMS)
    return;
#endif

    if (!trace_mode)
        return;

    in_trace_routine = 1;

    if (trace_level<0) {
        fprintf(stderr, "error: trace level is negative (log_entry)\n");
        fprintf(stderr, "calling routine is %s\n", routine);
        exit(1);
        }

    if (trace_mode&TRACE_MEMORY_LEVEL)
        memory_level = memory_count();

    if (trace_mode&TRACE_ENTRY) {
        /* keep track of calls in a file and in an internal list */
        if (max_level<=trace_level) {
            file_pos = trealloc(file_pos, sizeof(*file_pos)*(max_level=trace_level+1));
            routine_name  = trealloc(routine_name, sizeof(*routine_name)*(max_level));
            file_pos[trace_level] = ftell(fp);
            routine_name[trace_level]  = tmalloc(sizeof(**routine_name)*(MAX_LENGTH+1));
            }
        else {
            fseek(fp, file_pos[trace_level], 0);
            }
        strncpy(routine_name[trace_level], routine, MAX_LENGTH);
        if ((len=strlen(routine))>MAX_LENGTH)
            len = MAX_LENGTH;
        fwrite(routine, sizeof(*routine), len=strlen(routine), fp);
        if (len<MAX_LENGTH)
            fwrite(blank, sizeof(*blank), MAX_LENGTH-len, fp);
        fputc('\n', fp);
        fflush(fp);
        }
    if (trace_mode&TRACEBACK_ON) {
        /* keep track of calls internally */
        if (max_level<=trace_level) {
            routine_name  = trealloc(routine_name, sizeof(*routine_name)*(max_level=trace_level+1));
            routine_name[trace_level]  = tmalloc(sizeof(**routine_name)*(MAX_LENGTH+1));
            }
        strncpy(routine_name[trace_level], routine, MAX_LENGTH);
        }
    trace_level++;
    in_trace_routine = 0;
    }

void log_exit(char *routine)
{
    long memlev;

#if defined(VAX_VMS)
    return;
#endif
    if (!trace_mode)
        return;

    in_trace_routine = 2;
    if (trace_mode&TRACE_MEMORY_LEVEL) {
        memlev = memory_count();
        if (memlev!=memory_level) {
            fprintf(fpmem, "memory changed to %ld inside %s\n", memlev, 
                    routine_name[trace_level-1]?routine_name[trace_level-1]:"{NULL}");
            fflush(fpmem);
            }
        memory_level = memlev;
        }
    if (trace_mode&TRACE_ENTRY) {
        if (trace_level<=0) {
            fprintf(stderr, "error: trace level is nonpositive (log_exit)\n");
            fprintf(stderr, "calling routine is %s\n", routine);
            exit(1);
            }
        fseek(fp, file_pos[trace_level-1], 0);
        fwrite(blank, sizeof(*blank), MAX_LENGTH, fp);
        fflush(fp);
        }
    trace_level--;
    in_trace_routine = 0;
    }

void traceback_handler(int sig)
{
    long i;
    switch (sig) {
#if !defined(_WIN32)
        case SIGHUP: fprintf(stderr, "\nTerminated by SIGHUP\n"); break;
        case SIGQUIT: fprintf(stderr, "\nTerminated by SIGQUIT\n"); break;
        case SIGTRAP: fprintf(stderr, "\nTerminated by SIGTRAP\n"); break;
        case SIGBUS: 
            fprintf(stderr, "\nTerminated by SIGBUS"); 
            break;
#endif
        case SIGINT: fprintf(stderr, "\nTerminated by SIGINT\n"); break;
        case SIGABRT: fprintf(stderr, "\nTerminated by SIGABRT\n"); break;
        case SIGILL: fprintf(stderr, "\nTerminated by SIGILL\n"); break;
        case SIGFPE: 
            fprintf(stderr, "\nTerminated by SIGFPE"); 
            break;
        case SIGSEGV: 
            fprintf(stderr, "\nTerminated by SIGSEGV"); 
            break;
        default:      fprintf(stderr, "\nTerminated by unknown signal\n"); break;
        }   
    fprintf(stderr, "Program trace-back:\n");
    for (i=0; i<trace_level; i++) {
        fputs(routine_name[i], stderr);
        fputc('\n', stderr);
        }
    if (in_trace_routine==1)
        fprintf(stderr, "log_entry\n");
    else if (in_trace_routine==2)
        fprintf(stderr, "log_exit\n");
    fflush(stderr);
    fflush(stderr);    /* to force flushing of output sent to stderr by other parts of the code */
    exit(1);
    }
