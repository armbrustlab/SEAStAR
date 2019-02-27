/*
 --------------------------------------------------------------------------- #
 Center for Environmental Genomics
 Copyright (C) 2009-2013 University of Washington.

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

/*
   NDEBUG is required to avoid linker errors in XCode related to a known
   interaction between OpenMP and the __builtin_expect hint used by assert.h
   This may be safely disabled for debugging purposes when using GCC on any
   platform.
*/
#define NDEBUG 1

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

// LUT for rev comp transliteration
static const char ss_rc_code[] =
"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXTVGHXXCDXXMXKNXXXYSAXBWXRXXXXXXXtvghXXcdXXmXknXXXysaXbwxrXXXXXX";
//                                                                ABCD  GH  K MN   RST VWXY       abcd  gh  k mn   rst vwxy

/////////////////////////
// Common macros
/////////////////////////

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
utstring_bincpy(s,b,strlen(b))

#define ss_rev_utstring(s)                               \
do {                                                     \
char ss_ut_c;                                            \
for (int i = 0; i < ((s)->i)>>1; i++) {                  \
ss_ut_c = (s)->d[i];                                     \
(s)->d[i]=(s)->d[(s)->i - i - 1];                        \
(s)->d[(s)->i - i - 1]= ss_ut_c;                         \
}                                                        \
} while(0)

#define ss_rev_comp_utstring(s)                          \
do {                                                     \
char ss_ut_c;                                            \
for (int i = 0; i < ((s)->i)>>1; i++) {                  \
ss_ut_c = (s)->d[i];                                     \
(s)->d[i] = ss_rc_code[(s)->d[(s)->i - i - 1]];          \
(s)->d[(s)->i - i - 1] = ss_rc_code[ss_ut_c];            \
}                                                        \
} while(0)

/////////////////////////
// JSON stdout macros
/////////////////////////

#define JSON_BEG_OBJ printf("{");
#define JSON_STR_PROP(key,value) printf("\"%s\":\"%s\"",key,value);
#define JSON_FLT_PROP(key,value) printf("\"%s\":%f",key,value);
#define JSON_FLT15_PROP(key,value) printf("\"%s\":%.15f",key,value);
#define JSON_INT_PROP(key,value) printf("\"%s\":%d",key,value);
#define JSON_LIT_PROP(key,value) printf("\"%s\":%s",key,value);
#define JSON_SUB_PROP(key) printf("\"%s\" : ",key);
#define JSON_END_OBJ printf("}");
#define JSON_BEG_ARR printf("[");
#define JSON_STR_ELEM(value) printf("\"%s\"",value);
#define JSON_FLT_ELEM(value) printf("%f",value);
#define JSON_FLT1_ELEM(value) printf("%.1f",value);
#define JSON_INT_ELEM(value) printf("%d",value);
#define JSON_LIT_ELEM(value) printf("%s",value);
#define JSON_END_ARR printf("]");
#define JSON_COMMA printf(",");

#define JSON_N_STR_ARRAY_OBJ(key,vals,cont,comma) if (vals->count) {  \
        JSON_SUB_PROP(key); \
        JSON_BEG_ARR; \
        JSON_STR_ELEM(vals->cont[0]); \
        for (int i=1; i<vals->count; i++) { \
            JSON_COMMA; \
            JSON_STR_ELEM(vals->cont[i]); \
        } \
        JSON_END_ARR; \
        if (comma) { JSON_COMMA; } \
      }

#define JSON_N_FLT_ARRAY_OBJ(key,array,num,comma) if (array) {  \
        JSON_SUB_PROP(key); \
        JSON_BEG_ARR; \
        JSON_FLT1_ELEM(array[1]); \
        for (int i=2; i<=num; i++) { \
            JSON_COMMA; \
            JSON_FLT1_ELEM(array[i]); \
        } \
        JSON_END_ARR; \
        if (comma) { JSON_COMMA; } \
      }

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

int ss_is_colorspace_fastq(char *fn);

int ss_is_fastq(char *fn);

int ss_check_fastq_files(char *r1, char *r2, char *s1, char *s2);

gzFile ss_get_gzFile(char *fn, char *mode);
