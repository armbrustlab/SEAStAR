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
 Name        : seastar_shared.c
 Description : Library of functions used in common by multiple SEAStAR 
               subprojects 
 =============================================================================
*/


/////////////////////////
// Function prototypes
/////////////////////////

#include <seastar_shared.h>

#define UTSTRING_DEFAULT_LEN 10

// Helper Macros

#define ss_utstring_freebuflen(s) ((unsigned)((s)->n - (s)->i))
#define ss_utstring_setlen(s) ((s)->i = strlen(s->d))

/////////////////////////////////////////////////////////////////
// 
// ss_stream_writer
//
// This function opens a text file (gzipped or not) and shoves
// reads from an open file handle (e.g. a pipe) into this file.
// This is useful for easily implementing a "writer thread".
//

long unsigned int ss_stream_writer(UT_string *fn, int pipe_fd, int zipped) {
    
    long unsigned int cnt = 0;
    
    UT_string *data;
    utstring_new(data);
    
    FILE *file = NULL;
    
    FILE *pipe_out = fdopen(pipe_fd,"r"); 
    
    if (!(utstring_len(fn))) {
        fclose(pipe_out);
        return(0);
    }
    
    if (ss_get_utstring(pipe_out, data)) {
        
        if (zipped) {  
            
            // Add a trailing .gz if not present
            if ((utstring_len(fn) < 4) || (strcmp(utstring_body(fn)+utstring_len(fn)-3, ".gz"))) {
                ss_strcat_utstring(fn, ".gz");
            }
            
            if (!(file = gzopen(utstring_body(fn), "w"))) {
                fprintf(stderr, "Error! Could not open filename for output: %s\n", utstring_body(fn));
                exit(EXIT_FAILURE);
            }
            
            do {
                gzputs(file, utstring_body(data));
                cnt++;
            } while (ss_get_utstring(pipe_out, data));
            
            gzclose(file); 
            
        } else {
            
            // Trim off a trailing .gz if present
            if ((utstring_len(fn) >= 4) && !(strcmp(utstring_body(fn)+utstring_len(fn)-3, ".gz"))) {
                ss_trunc_utstring(fn, 3);
            }
            
            if (!(file = fopen(utstring_body(fn), "w"))) {
                fprintf(stderr, "Error! Could not open filename for output: %s\n", utstring_body(fn));
                exit(EXIT_FAILURE);
            }
            
            do {
                fputs(utstring_body(data), file);
                cnt++;
            } while (ss_get_utstring(pipe_out, data));
            
            fclose(file); 
        }
    }
    
    fclose(pipe_out);
    
    utstring_free(data);
    
    return(cnt);
}
 
/////////////////////////////////////////////////////////////////
// 
// ss_stream_reader
//
// This function opens a text file (gzipped or not) and shoves
// the resulting reads into another file handle (eg. a pipe)
// This is useful for easily implementing a "reader thread".
//

long unsigned int ss_stream_reader(UT_string *fn, int pipe_fd) {
    
    long unsigned int cnt = 0;
    
    UT_string *data;
    utstring_new(data);
    
    FILE *file = NULL;
    
    FILE *pipe_in = fdopen(pipe_fd,"w");
    
    if (!(utstring_len(fn))) {
        fclose(pipe_in);
        return(0);
    }
    
    // Try to open the file
    if (!(file = gzopen(utstring_body(fn), "r"))) {
        utstring_printf(fn, ".gz");
        if (!(file = gzopen(utstring_body(fn), "r"))) {
            fclose(pipe_in);
            return(0);
        }
    } 
    
    while (ss_gzget_utstring(file, data)) {
        fputs(utstring_body(data), pipe_in);
        cnt++;
    }
    
    fclose(pipe_in);
    gzclose(file);
    utstring_free(data);
    
    return(cnt);
}

/////////////////////////////////////////////////////////////////
// 
// ss_get_utstring
//
// This function abstracts the fgets line reading function for
// use with UT_string type variable length strings.
//

char *ss_get_utstring(FILE *f, UT_string *s) {
    
    char *retval = NULL;
    
    utstring_clear(s);
    
    if (ss_utstring_freebuflen(s) <= 1) {
        utstring_reserve(s, UTSTRING_DEFAULT_LEN);
    }
    
    retval = fgets(utstring_body(s), ss_utstring_freebuflen(s), f);
    ss_utstring_setlen(s);
    
    while (retval && (utstring_body(s)[utstring_len(s)-1] != '\n')) {
        utstring_reserve(s, utstring_len(s)); // Double the alloced buf
        retval = fgets(utstring_body(s)+utstring_len(s), ss_utstring_freebuflen(s), f);
        ss_utstring_setlen(s);
    }
    
    return(utstring_len(s) ? utstring_body(s) : NULL);
}

/////////////////////////////////////////////////////////////////
// 
// ss_get_cat_utstring
//
// This function abstracts the fgets line reading function for
// use with UT_string type variable length strings.
//

char *ss_get_cat_utstring(FILE *f, UT_string *s) {
    
    char *retval = NULL;
    
    if (ss_utstring_freebuflen(s) <= 1) {
        utstring_reserve(s, UTSTRING_DEFAULT_LEN);
    }

    retval = fgets(utstring_body(s)+utstring_len(s), ss_utstring_freebuflen(s), f);
    ss_utstring_setlen(s);
    
    while (retval && (utstring_body(s)[utstring_len(s)-1] != '\n')) {
        utstring_reserve(s, utstring_len(s)); // Double the alloced buf
        retval = fgets(utstring_body(s)+utstring_len(s), ss_utstring_freebuflen(s), f);
        ss_utstring_setlen(s);
    }
    
    return(utstring_len(s) ? utstring_body(s) : NULL);
}

/////////////////////////////////////////////////////////////////
// 
// ss_gzget_utstring
//
// This function abstracts the gzgets line reading function for
// use with UT_string type variable length strings.
//

char *ss_gzget_utstring(gzFile f, UT_string *s) {
        
    char *retval = NULL;
    
    utstring_clear(s);
    
    if (ss_utstring_freebuflen(s) <= 1) {
        utstring_reserve(s, UTSTRING_DEFAULT_LEN);
    }
    
    retval = gzgets(f, utstring_body(s), ss_utstring_freebuflen(s));
    ss_utstring_setlen(s);
    
    while (retval && (utstring_body(s)[utstring_len(s)-1] != '\n')) {
        utstring_reserve(s, utstring_len(s)); // Double the alloced buf
        retval = gzgets(f, utstring_body(s)+utstring_len(s), ss_utstring_freebuflen(s));
        ss_utstring_setlen(s);
    }
    
    return(utstring_len(s) ? utstring_body(s) : NULL);
}

/////////////////////////////////////////////////////////////////
// A portable pseudo random number generator
// This is a re-entrant implementation of the Mersenne Twister
// algortihm. 
//
// To use: 
// 1) declare a context pointer variable of type ss_rand_inst
// 2) call ss_rseed with a starting seed value, assign the return to the pointer
// 3) call ss_rand with the conext pointer. It returns a double between 0.0-1.0
// 4) repeat. When finished, call free on the context pointer.
//

ss_rand_inst *ss_rseed(unsigned int seed) {
    
    ss_rand_inst *r_ptr = NULL;
    
    // Malloc an instance structure for this generator
    if (!(r_ptr = malloc(sizeof(ss_rand_inst)))) {  
        return (NULL);        
    }
    
    r_ptr->idx = 0;
    
    r_ptr->rbuf[0] = seed;
    
    // Initialize the generator state
    for (int i = 1; i < 624; i++) {
        r_ptr->rbuf[i] = (0x6c078965 * (r_ptr->rbuf[i-1] ^ ((r_ptr->rbuf[i-1])>>30)) + i) & 0xffffffff;
    }
    
    return (r_ptr);
}

//////////////////////////////////////////////////

double ss_rand(ss_rand_inst *r_ptr) {
    
    unsigned int y = 0;
    
    if (r_ptr->idx == 0) { // Generate new untempered numbers once all have been used
        
        for (int i = 0; i < 624; i++) {
            
            y = (r_ptr->rbuf[i] & 0x8000000) + (r_ptr->rbuf[(i+1) % 624] & 0x7FFFFFFF);
            r_ptr->rbuf[i] = r_ptr->rbuf[(i+397) % 624] ^ (y >> 1);
            
            if ((y % 2) != 0) {
                r_ptr->rbuf[i] = r_ptr->rbuf[i] ^ 0x9908b0df;   
            }
        }
    }
    
    // generate a tempered pseudo random number from the current state
    y = r_ptr->rbuf[r_ptr->idx];
    y = y ^ (y >> 11);
    y = y ^ ((y << 7) & 0x9d2c5680);
    y = y ^ ((y << 15) & 0xefc60000);
    y = y ^ (y >> 18);
    
    r_ptr->idx = (r_ptr->idx + 1) % 624;
    
    return((double) y / (unsigned int) 0xffffffff);
}

/////////////////////////////////////////////////////////////////
// 
// ss_is_colorspace_fastq
//
// This function checks if a file is in colorspace fastq format.
// It checks for the file name as given and a gzipped form with 
// ".gz" suffix.
//

int ss_is_colorspace_fastq(char *fn) {
    int colorspace = 0;
    gzFile fq_hnd = NULL;
    UT_string *str;
    UT_string *file_name;
    utstring_new(str);
    utstring_new(file_name);
    
    utstring_printf(file_name, fn);
    
    if (ss_is_fastq(utstring_body(file_name))) {
        if (fq_hnd = ss_get_gzFile(utstring_body(file_name), "r")) {
            ss_gzget_utstring(fq_hnd, str);
            
            // Detect colorspace read alignments by the naming convention that preserves the
            // primer base and first color in the header line between '@' and '+'.
            // e.g.  @TA+49|lambda:1231_1682_1381
            //        ^^
            // This checks syntax and protects against buffer overflow below.
            if (strspn(utstring_body(str), "@+ACGTNacgtn") >= 3) {
                if (utstring_body(str)[3] == '+') {
                    colorspace = 1;
                }
            }
            
            // Close file
            gzclose(fq_hnd);
        }
    }
    
    utstring_free(str);
    utstring_free(file_name);
    
    return colorspace;
}

/////////////////////////////////////////////////////////////////
// 
// ss_is_fastq
//
// This function checks if a file is in valid fastq format.
// It checks for the file name as given and a gzipped form with 
// ".gz" suffix.
//

int ss_is_fastq(char *fn) {
    int fastq = 0;
    gzFile fq_hnd = NULL;
    UT_string *str;
    UT_string *file_name;
    utstring_new(str);
    utstring_new(file_name);
    
    utstring_printf(file_name, fn);
    
    if (fq_hnd = ss_get_gzFile(utstring_body(file_name), "r")) {
        ss_gzget_utstring(fq_hnd, str);
        
        if (utstring_body(str)[0] == '@') {
            ss_gzget_utstring(fq_hnd, str);
            ss_gzget_utstring(fq_hnd, str);
            if (utstring_body(str)[0] == '+') {
                fastq = 1;
            }
        }
        // Close file
        gzclose(fq_hnd);
    }
    
    utstring_free(str);
    utstring_free(file_name);
    
    return fastq;
}

/////////////////////////////////////////////////////////////////
// 
// ss_check_fastq_files
//
// This function checks if files are all colorspace or all nucleotide.
// If data is mixed colorspace and nucleotide, exit with an error.
// If data is all colorspace return 1.  If all nucleotide return 0.
// Files are checked with names as given and with a ".gz" suffix.
//

int ss_check_fastq_files(char *r1, char *r2, char *s1, char *s2) {
    // Check if input fastq is colorspace
    // If some files are colorspace and some are basespace, throw an error
    int fcount = 0;
    int cscount = 0;
    fcount += ss_is_fastq(r1);
    fcount += ss_is_fastq(r2);
    fcount += ss_is_fastq(s1);
    fcount += ss_is_fastq(s2);
    cscount += (ss_is_fastq(r1) && ss_is_colorspace_fastq(r1));
    cscount += (ss_is_fastq(r2) && ss_is_colorspace_fastq(r2));
    cscount += (ss_is_fastq(s1) && ss_is_colorspace_fastq(s1));
    cscount += (ss_is_fastq(s2) && ss_is_colorspace_fastq(s2));

    if (cscount && (cscount != fcount)) {
        fprintf(stderr, "Error: Mixed colorspace and basespace FASTQ files detected\n");
        exit(EXIT_FAILURE);
    }
    return cscount ? 1 : 0;  // 1 if colorspace, 0 if nucleotide
}

/////////////////////////////////////////////////////////////////
//
// ss_get_gzFile
//
// This function tries to open the file with name as given and failing that
// tries to open the file with ".gz" extension.  It uses gzopen and returns a 
// gzFile.  If file opening fails a NULL pointer is returned.
//

gzFile ss_get_gzFile(char *fn, char *mode) {
    gzFile fq_hnd = NULL;
    UT_string *file_name;
    utstring_new(file_name);
    
    utstring_printf(file_name, fn);
    if (!(fq_hnd = gzopen(utstring_body(file_name), mode))) {
        utstring_printf(file_name, ".gz");
        fq_hnd = gzopen(utstring_body(file_name), mode);
    }
    return fq_hnd;
}
