/* Copyright 1994 by Michael Borland and Argonne National Laboratory,
 * all rights reserved.
 */
/* file: subprocess.c
 *
 * Michael Borland, 1994
 */
#include <stdio.h>
#if defined(_WIN32)
#include <process.h>
#else
#include <unistd.h>
#endif
#include "mdb.h"
#include "track.h"
#include "subprocess.h"
#include <signal.h>

#if defined(_WIN32)
#define popen(command, mode) _popen(command, mode)
#endif

static FILE *fp = NULL;
static int pid;

/* dummy signal handler for use with sigpause */
void subprocess_sigusr1()
{
}

void run_subprocess(NAMELIST_TEXT *nltext, RUN *run)
{
  static char buffer[1024];
  char *ptr, *ptr0;

  log_entry("run_subprocess");

  /* process the namelist text */
  set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
  set_print_namelist_flags(0);
  process_namelist(&subprocess, nltext);
  print_namelist(stdout, &subprocess);

  if (command) {
    buffer[0] = 0;
    ptr0 = command;
    while ((ptr=strstr(ptr0, "%s"))) {
      if (ptr!=command && *(ptr-1)=='%') {
        *(ptr-1) = 0;
        strcat(buffer, ptr0);
        strcat(buffer, "%s");
        ptr += 2;
        ptr0 = ptr;
      }
      else {
        if (!run || !run->rootname)
          bomb("rootname must be initialized prior to subprocess execution if \%s substitution is used", NULL);
        *ptr = 0;
        ptr += 2;
        strcat(buffer, ptr0);
        strcat(buffer, run->rootname);
        ptr0 = ptr;
      }
    }
    strcat(buffer, ptr0);
    fprintf(stdout, "%s\n", buffer);
    fflush(stdout);
    executeCshCommand(buffer);
  }

  log_exit("run_subprocess");
}

void executeCshCommand(char *cmd)
{
#if !defined(_WIN32)
  signal(SIGUSR1, subprocess_sigusr1);
#endif

  if (!fp) {
    /* open a pipe and start csh */
    fp = popen("csh", "w");
    pid = getpid();
  }
  
  fprintf(fp, "%s\nkill -USR1 %d\n", cmd, pid);
  fflush(fp);

#if !defined(_WIN32)
  /* pause until SIGUSR1 is received */
  sigpause(SIGUSR1);

  /* back to default behavior for sigusr1 */
  signal(SIGUSR1, SIG_DFL);
#endif
}

