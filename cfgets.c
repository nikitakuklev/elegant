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

#if defined(UNIX) && defined(GNU_C)
long toupper(char c);
long tolower(char c);
#endif

void delete_spaces(char *s);
void str_to_upper_quotes(char *s);

char *cfgets(char *s, long n, FILE *fpin)
{
    register long l;
    static long level = 0;

    level++;
    while (fgets(s, n, fpin)) {
        chop_nl(s);
        delete_spaces(s);
        if (s[0]!='!') {
            l = strlen(s);
            while (l!=0 && s[l]<27)
                l--;
            if (s[l]=='&') {
                cfgets(s+l, n-l, fpin);
                if (level==1) 
                    str_to_upper_quotes(s);
                level--;
                return(s);
                }
            else {
                if (level==1) 
                    str_to_upper_quotes(s);
                level--;
                return(s);
                }
            }
        }
    return(NULL);
    }

void delete_spaces(char *s)
{
    char *ptr;
    while (*s) {
        if (*s=='"') {
            s++;
            while (*s && *s!='"')
                s++;
            if (*s=='"')
                s++;
            }
        else if (*s==' ') {
            ptr = s++;
            while (*s==' ')
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
    while (*s) {
        if (*s=='"') {
            while (*++s && *s!='"')
                ;
            }
        else 
            *s = toupper(*s);
        s++;
        }
    }
