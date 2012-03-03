/*
 --------------------------------------------------------------------------- #
 Center for Environmental Genomics
 Copyright (C) 2009-2012 University of Washington.
 
 Authors:
 Vaughn Iverson
 vsi@uw.edu
 --------------------------------------------------------------------------- #
 This file is part of SEAStAR.
 
 SEAStAR is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 SEAStAR is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with SEAStAR.  If not, see <http://www.gnu.org/licenses/>.
 --------------------------------------------------------------------------- #
 
 =============================================================================
 Name        : seastar_shared.h
 Description : Library of functions used in common by multiple SEAStAR 
               subprojects 
 =============================================================================
*/

#ifndef _STDIO_H_
#include <stdio.h>
#endif

#ifndef _STDLIB_H_
#include <stdlib.h>
#endif

#ifndef _STRING_H_
#include <string.h>
#endif

#ifndef	_CTYPE_H_
#include <ctype.h>
#endif

#ifndef __MATH_H__
#include <math.h>
#endif

#ifndef ZLIB_H
#include <zlib.h>
#endif

#ifndef UTSTRING_H
#include <utstring.h>
#endif

// Build a default --version string for all tools when SS_BUILD_VERSION isn't passed to the compiler by the build env.

#ifndef SS_BUILD_VERSION
#define SS_BUILD_VERSION "No version information available. Executable built on: " __DATE__ ", " __TIME__ 
#endif

// Common macros

#define ss_trunc_utstring(s, n)                          \
do {                                                     \
(s)->d[(s)->i-(n)] = '\0';                               \
(s)->i = (s)->i - (n);                                   \
} while(0)

#define ss_trim_utstring(s, n)                           \
do {                                                     \
(s)->d[(n)] = '\0';                                      \
(s)->i = (n);                                            \
} while(0)

#define ss_strcat_utstring(s,b)                          \
do {                                                     \
size_t ss_ut_len = strlen(b);                            \
utstring_reserve(s,ss_ut_len);                           \
strcat(&(s)->d[(s)->i], b);                              \
(s)->i += ss_ut_len;                                     \
(s)->d[(s)->i]='\0';                                     \
} while(0)

#define ss_rev_utstring(s)                               \
do {                                                     \
char ss_ut_c;                                            \
for (int i = 0; i < ((s)->i)>>1; i++) {                  \
ss_ut_c = (s)->d[i];                                     \
(s)->d[i]=(s)->d[(s)->i - i - 1];                        \
(s)->d[(s)->i - i - 1]= ss_ut_c;                         \
}                                                        \
} while(0)
        
/////////////////////////
// Typedefs 
/////////////////////////

typedef struct ss_rand_struct {
    
    int idx;
    unsigned int rbuf[624];
    
} ss_rand_inst;

/////////////////////////
// Prototypes 
/////////////////////////

ss_rand_inst *ss_rseed(unsigned int seed);

double ss_rand(ss_rand_inst *r_ptr);

char *ss_get_utstring(FILE *f, UT_string *s);

char *ss_get_cat_utstring(FILE *f, UT_string *s);

char *ss_gzget_utstring(gzFile f, UT_string *s);

long unsigned int ss_stream_writer(UT_string *fn, int pipe_fd, int zipped);

long unsigned int ss_stream_reader(UT_string *fn, int pipe_fd);

