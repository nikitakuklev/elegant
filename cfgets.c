/* Copyright 1994 by Michael Borland and Argonne National Laboratory,
 * all rights reserved.
 */
/* prototypes for this file are in prot.out */
/* routine: cfgets()
 * purpose: read from a file, skipping MAD-style comments.
 *          Returns the next non-comment line, including continuations 
 *          (marked with an & in the last column).  The routine also 
 *          deletes all spaces from the line, as a kludge for the parsing 
 *          routines.
 * 
 * Michael Borland, 1987
 */
#include "mdb.h"
#include "track.h"
#include <ctype.h>
#if defined(_WIN32)
#include <stdlib.h>
#else
#if defined(UNIX) && defined(GNU_C)
long toupper(char c);
long tolower(char c);
#endif
#endif

void delete_spaces(char *s);
void str_to_upper_quotes(char *s);
char *cfgets1(char *s, long n, FILE *fpin);

char *cfgets(char *s, long n, FILE *fpin)
{
  s[0] = 0;
  if (!cfgets1(s, n, fpin))
    return NULL;
  str_to_upper_quotes(s);
  return s;
}

char *cfgets1(char *s, long n, FILE *fpin)
{
  register long l;
  
  while (fgets(s, n, fpin)) {
    if (s[0]=='!')
      continue;
    chop_nl(s);
    delete_spaces(s);
    l = strlen(s);
    while (l!=0 && s[l]<27)
      l--;
    if (s[l]=='&') {
      s[l] = 0;
      cfgets(s+l, n-l, fpin);
      return(s);
    }
    else {
      return(s);
    }
  }
  return NULL;
}

void delete_spaces(char *s)
{
    char *ptr, *ptr0;
    ptr0 = s;
    while (*s) {
        if (*s=='"' && (ptr0==s || *(s-1)!='\\')) {
            s++;
            while (*s && (*s!='"' || *(s-1)=='\\'))
                s++;
            if (*s=='"' && *(s-1)!='\\')
                s++;
            }
        else if (*s==' ' || *s=='\011') {
            ptr = s++;
            while (*s==' ' || *s=='\011')
                s++;
            strcpy(ptr, s);
            s = ptr;
            }
        else
            s++;
        }
    }

void str_to_upper_quotes(char *s)
{
  char *ptr0;
  ptr0 = s;
  if (!s)
    return;
  while (*s) {
    if (*s=='"' && (ptr0==s || *(s-1)!='\\')) {
      while (*++s && (*s!='"' || *(s-1)=='\\'))
        ;
    }
    else 
      *s = toupper(*s);
    s++;
  }
}
