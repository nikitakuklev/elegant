/* $Log: not supported by cvs2svn $
 */
/* program: replaceText
 * purpose: replace strings in a file.
 * a simplified variant of the replace program.
 * 
 * Michael Borland, 2000
 */
#include "mdb.h"
#include "scan.h"

#define SET_ORIGINAL 0
#define SET_REPLACEMENT 1
#define SET_VERBOSE 2
#define SET_STRINGS 3
#define SET_STDIO 4
#define SET_FILTER 5
#define N_OPTIONS 6

char *option[N_OPTIONS] = {
  "original", "replacement", "verbose", "strings", "stdio", "filter",
} ;

#define USAGE "replaceText {input [output] | -stdio | -strings=string[,string...]}\n\
-original=string[,string...] -replacement=string[,string...]\n\
[-verbose] [-filter=string[,string...]]\n\n\
-stdio        Take input from standard input and deliver output to standard output.\n\
-strings      Perform replacement on the listed strings.\n\
-original     List of strings to replace.\n\
-replacement  Parallel list of replacements for original strings.\n\
-verbose      Report the number of replacements made.\n\
-filter       Perform replacement only on lines containing one or more of the\n\
              strings in the filter list.\n\n\
Program by Michael Borland.  (This is version 1, April 2000)."

int replace_string1(char *t, char *s, char *orig, char *repl);
long passesFilters(char *s, char **filter, long filters);

main(
    int argc,
    char **argv
    )
{
  char **orig, **repl, *input, *output, *ptr, **filter;
  char **string_list;
  int n_pairs, i_pair, filters;
  char s[1024], t[1024];
  FILE *fpi, *fpo;
  SCANNED_ARG *s_arg;
  int i_arg, *count, n_strings, use_stdio;
  char *input_fn, *output_fn;
  int tmp_file, verbose;

  argc = scanargs(&s_arg, argc, argv);
  if (argc<2 || argc>(3+N_OPTIONS)) 
    bomb("too few or too many arguments", USAGE);

  input = output = NULL;
  orig = repl = filter =  NULL;
  n_pairs = verbose = filters = 0;
  string_list = NULL;
  use_stdio = tmp_file = 0;

  for (i_arg=1; i_arg<argc; i_arg++) {
    if (s_arg[i_arg].arg_type==OPTION) {
      switch (match_string(s_arg[i_arg].list[0], option, N_OPTIONS, 0)) {
      case SET_ORIGINAL:
        if ((i_pair=s_arg[i_arg].n_items-1)<1 ||
            (n_pairs && i_pair!=n_pairs) ||
            !(orig=s_arg[i_arg].list+1))
          bomb("invalid -original syntax", USAGE);
        n_pairs = i_pair;
        break;
      case SET_REPLACEMENT:
        if ((i_pair=s_arg[i_arg].n_items-1)<1 ||
            (n_pairs && i_pair!=n_pairs) ||
            !(repl=s_arg[i_arg].list+1))
          bomb("invalid -replacement syntax", USAGE);
        n_pairs = i_pair;
        break;
      case SET_VERBOSE:
        verbose = 1;
        break;
      case SET_STRINGS:
        if (s_arg[i_arg].n_items<2)
          bomb("invalid -strings syntax", USAGE);
        string_list = s_arg[i_arg].list+1;
        n_strings = s_arg[i_arg].n_items-1;
        break;
      case SET_STDIO:
        use_stdio = 1;
        break;
      case SET_FILTER:
        if (s_arg[i_arg].n_items<2)
          bomb("invalid -filter syntax", USAGE);
        filter = s_arg[i_arg].list+1;
        filters = s_arg[i_arg].n_items-1;
        break;
      default:
        bomb("unknown option given", USAGE);
        break;
      }
    }
    else {
      if (!input)
        fpi = fopen_e(input=s_arg[i_arg].list[0], "r", 0);
      else if (!output)
        fpo = fopen_e(output=s_arg[i_arg].list[0], "w", FOPEN_SAVE_IF_EXISTS);
      else
        bomb("too many filenames listed", USAGE);
    }
  }                
  
  if (!(repl && orig))
    bomb("insufficient information supplied", USAGE);
  if (!input && !string_list && !use_stdio)
    bomb("unsufficient information supplied", USAGE);

  count = tmalloc(sizeof(int)*n_pairs);
  for (i_pair=0; i_pair<n_pairs; i_pair++) {
    interpret_escapes(orig[i_pair]);
    interpret_escapes(repl[i_pair]);
    count[i_pair] = 0;
  }

  if (string_list) {
    int i_string;
    for (i_string=0; i_string<n_strings; i_string++) {
      strcpy(s, string_list[i_string]);
      for (i_pair=0; i_pair<n_pairs; i_pair++) {
        count[i_pair] += replace_string1(t, s, orig[i_pair], repl[i_pair]);
        strcpy(s, t);
      }
      puts(t);
    }
    exit(0);
  }

  if (input) {
    tmp_file = 0;
#ifndef VAX_VMS
    if (!output) {
      output = tmpname(NULL);
      tmp_file = 1;
      fpo = fopen_e(output, "w", 0);
    }
#else
    if (!output)
      fpo = fopen_e(output=input, "w", 0);
#endif
    
    clean_filename(input);
    cp_str(&input_fn, input);
    if (ptr=strchr(input, '.'))
      *ptr = 0;
    
    if (output) {
      cp_str(&output_fn, output);
      clean_filename(output);
      if (ptr=strchr(output, '.')) 
        *ptr = 0;
    }
    else 
      output_fn = input_fn;
  }
  else {
    fpi = stdin;
    fpo = stdout;
    verbose = 0;
  }
  
  while (fgets(s, 1024, fpi)) {
    if (!filters || passesFilters(s, filter, filters))
      for (i_pair=0; i_pair<n_pairs; i_pair++) {
        count[i_pair] += replace_string1(t, s, orig[i_pair], repl[i_pair]);
        strcpy(s, t);
      }
    if (fputs(s, fpo)==EOF) {
      puts("error writing to file--program aborted"); 
      exit(1);
    }
  }

  if (verbose)
    for (i_pair=0; i_pair<n_pairs; i_pair++) 
      printf("%d replacements of \"%s\" by \"%s\"\n",
             count[i_pair], orig[i_pair], repl[i_pair]);

  if (tmp_file) {
    sprintf(s, "%s~", input_fn);
    rename(input_fn, s);
    rename(output_fn, input_fn);
  }
  exit(0);
}

int replace_string1(char *t, char *s, char *orig, char *repl)
{
  char *ptr0, *ptr1;
  int count, lorig;
  char temp;

  ptr0 = s; 
  t[0] = 0;
  count = 0;
  lorig = strlen(orig);

  while ( ((ptr1=ptr0) && orig[0]=='^' && (orig[1]==0 || strncmp(ptr0, orig+1, strlen(orig+1))==0)) || 
         (orig[0]!='^' && (ptr1=str_in(ptr0, orig))) || 
         orig[lorig-1]=='$') {
    if (orig[lorig-1]=='$') {
      if (((long)strlen(ptr0))<=(lorig-1))
        break;
      if (ptr1=ptr0+strlen(ptr0)-(lorig-1)) {
        if (strncmp(ptr1, orig, lorig-1)!=0)
          break;
      }
      else if (lorig!=1)
        break;
    }
    if (orig[0]=='^') {
      strcat(t, repl);
      strcat(t, ptr0+lorig-1);
      return(1);
    }
    if (orig[lorig-1]=='$') {
      *ptr1 = 0;
      strcat(t, ptr0);
      strcat(t, repl);
      return(1);
    }
    count++;
    *ptr1 = 0;
    strcat(t, ptr0);
    strcat(t, repl);
    ptr0 = ptr1+lorig;
  }
  if (strlen(ptr0));
  strcat(t, ptr0);
  return(count);
}

long passesFilters(char *s, char **filter, long filters)
{
  while (filters--) {
    if (!strstr(s, *filter++))
      return 0;
  }
  return 1;
}
