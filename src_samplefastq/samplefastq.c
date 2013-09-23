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
 Name        : samplefastq.c
 Description : Read FASTQ format file(s) and randomly samplea the reads.
 =============================================================================
*/

#include <seastar_shared.h>

#include <argtable2.h>
#ifdef _OPENMP
#include <omp.h>
#endif

/////////////////////////
// Macros 
/////////////////////////

/////////////////////////
// Function prototypes
/////////////////////////

unsigned long int fq_stream_counter(UT_string *fq_fn);

unsigned long int fq_stream_sampler(UT_string *fq_fn, int pipe_fd, double sample_frac, unsigned int seed, unsigned long int *org);

/////////////////////////////////////////////////////////////////
// Entry point
int main(int argc, char *argv[]) {
    
#ifndef _OPENMP
    fprintf(stderr, "\nERROR: Program built with compiler lacking OpenMP support.\n");
    fprintf(stderr, "See SEAStAR README file for information about suitable compilers.\n");
    exit(EXIT_FAILURE);
#endif    
    
    ///////////////////////////
    // Variable declarations
    ///////////////////////////
    
    // Input filenames
    UT_string *in_read1_fq_fn, *in_read2_fq_fn, *in_single1_fq_fn, *in_single2_fq_fn;
    utstring_new(in_read1_fq_fn);
    utstring_new(in_read2_fq_fn);
    utstring_new(in_single1_fq_fn);
    utstring_new(in_single2_fq_fn);
    
    // Output filenames
    UT_string *out_read1_fq_fn, *out_read2_fq_fn, *out_single1_fq_fn, *out_single2_fq_fn;
    utstring_new(out_read1_fq_fn);
    utstring_new(out_read2_fq_fn);
    utstring_new(out_single1_fq_fn);
    utstring_new(out_single2_fq_fn);
    
    // Flags
    int singles_flag = 0;
    unsigned int seed_val = 0;
    double read_frac = 0.0;
    
    // Read counters
    unsigned long int mp_org = 0, R1_org = 0, R2_org = 0, singlet1_org = 0, singlet2_org = 0;
    unsigned long int mp_cnt = 0, R1_cnt = 0, R2_cnt = 0, singlet1_cnt = 0, singlet2_cnt = 0;

    // All done with variable declarations!!
    
    ///////////////////////////////////
    // Command line argtable settings
    ///////////////////////////////////
    
    struct arg_int *seed = arg_int0(NULL, "seed", "<int>", "Seed value for the pseudo-random number generator. [12345]");
    struct arg_dbl *sample = arg_dbl1("f", "frac_sample", "<float> or <int>", "Fraction (or count) of input reads (singlets+mate-pairs) to sample.");
    struct arg_rem *sample_r = arg_rem(NULL, "values < 1.0 are interpreted as a fraction of reads, >= 1.0 is converted to an integer count of reads.");
    struct arg_lit *gzip = arg_lit0("z", "gzip", "Output converted files in gzip compressed format. [NULL]");
    struct arg_file *input = arg_file1(NULL, NULL, "<in_prefix>", "Input file prefix: (e.g. <in_prefix>_single.fastq [<in_prefix>_read1.fastq <in_prefix>_read2.fastq]) ");
    struct arg_file *output = arg_file1(NULL, NULL, "<out_prefix>", "Output file prefix: (e.g. <out_prefix>_single.fastq [<out_prefix>_read1.fastq <out_prefix>_read2.fastq]) ");
    struct arg_lit *version = arg_lit0(NULL,"version","Print the build version and exit.");  
    struct arg_lit *h = arg_lit0("h", "help", "Request help.");
    struct arg_end *end = arg_end(20);
    
    void *argtable[] = {h,version,seed,sample,sample_r,gzip,input,output,end};
    int arg_errors = 0;
    
    ////////////////////////////////////////////////////////////////////////
    // Handle command line processing (via argtable2 library) 
    ////////////////////////////////////////////////////////////////////////
    
    arg_errors = arg_parse(argc, argv, argtable);
    
	if (version->count) {
		fprintf(stderr, "%s version: %s\n", argv[0], SS_BUILD_VERSION);
		exit(EXIT_SUCCESS);
    }
	
    if (h->count) {
        arg_print_syntaxv(stderr, argtable, "\n\n");
        arg_print_glossary(stderr, argtable, "%-25s %s\n");
        fprintf(stderr, "\nInput prefix is the part of the filename before the final\n");
        fprintf(stderr, "underscore \"_\" (eg. not including \"_read1.fq\").\n");
        fprintf(stderr, "The output prefix will be used to name the output files.\n");
        fprintf(stderr, "Mate-paired read files are automatically used if present.\n");
        fprintf(stderr, "\nNote: Input and output files may be gzipped.\n");
        fprintf(stderr, "\nAlso note: Using --frac_sample with a read count to sample \n");
        fprintf(stderr, "will result in reading the inputs twice. Use a fraction if possible.\n");        
        exit(EXIT_FAILURE);
    }    
    
    if (arg_errors) { 
        arg_print_errors(stderr, end, "samplefastq");
        arg_print_syntaxv(stderr, argtable, "\n");
        exit(EXIT_FAILURE);
    }
    
    if (seed->count) {
        seed_val = seed->ival[0];
    } else {
        seed_val = 12345;
    }
    
    // Check for null string prefixes
    if (!(strlen(input->filename[0]) && strlen(output->filename[0]))) {
        fprintf(stderr, "Error: NULL prefix strings are not permitted.\n");
        exit(EXIT_FAILURE);
    }

    // Construct input filenames
    utstring_printf(in_read1_fq_fn, "%s.read1.fastq", input->filename[0]);
    utstring_printf(in_read2_fq_fn, "%s.read2.fastq", input->filename[0]);

    utstring_printf(in_single1_fq_fn, "%s.single.fastq", input->filename[0]);
    
    FILE *in_single1_file = NULL;
    
    // Try to open a singlet fastq file
    if (!(in_single1_file = gzopen(utstring_body(in_single1_fq_fn), "r"))) {
        utstring_printf(in_single1_fq_fn, ".gz");
        
        if (!(in_single1_file = gzopen(utstring_body(in_single1_fq_fn), "r"))) {
            utstring_clear(in_single1_fq_fn);
            utstring_printf(in_single1_fq_fn, "%s.single1.fastq", input->filename[0]);
            utstring_printf(in_single2_fq_fn, "%s.single2.fastq", input->filename[0]);
            singles_flag = 1;
        } 
    }

    if (in_single1_file) {
        gzclose(in_single1_file);
    }
    
    // Output filenames

    utstring_printf(out_read1_fq_fn, "%s.read1.fastq", output->filename[0]);
    utstring_printf(out_read2_fq_fn, "%s.read2.fastq", output->filename[0]);
    
    if (singles_flag == 1) {
        utstring_printf(out_single1_fq_fn, "%s.single1.fastq", output->filename[0]);
        utstring_printf(out_single2_fq_fn, "%s.single2.fastq", output->filename[0]);
    } else {
        utstring_printf(out_single1_fq_fn, "%s.single.fastq", output->filename[0]);
    }
    
    /////////////////////////////////////////////////////////////////////////////////////////
    // Begin processing!

#ifdef _OPENMP    
    omp_set_num_threads(8);
#endif    
    
    if (sample->dval[0] > 1.0) {    // Calculate the fraction of reads            

#pragma omp parallel sections default(shared)       
        {
            
#pragma omp section 
            {   // Mate-pair counter
                mp_org = fq_stream_counter(in_read1_fq_fn);
            }
        
#pragma omp section 
            {   // Singlet counter
                singlet1_org = fq_stream_counter(in_single1_fq_fn);
            }
            
#pragma omp section 
            {   // Singlet file #2 counter
                singlet2_org = fq_stream_counter(in_single2_fq_fn);
            }            
            
        }    
        read_frac = sample->dval[0]/(double)(singlet1_org + singlet2_org + mp_org);
    } else {
        read_frac = sample->dval[0];
    }

    int r1_pipe[2];
    int r2_pipe[2];
    int s1_pipe[2];
    int s2_pipe[2];
    pipe(r1_pipe);
    pipe(r2_pipe);
    pipe(s1_pipe);
    pipe(s2_pipe);
    
#pragma omp parallel sections default(shared)       
    {
        
#pragma omp section 
        {   // Read1 reader
            R1_cnt = fq_stream_sampler(in_read1_fq_fn, r1_pipe[1], read_frac, seed_val, &R1_org);
        }
        
#pragma omp section 
        {   // Read1 writer
            ss_stream_writer(out_read1_fq_fn, r1_pipe[0], gzip->count);
        }
        
#pragma omp section 
        {   // Read2 reader
            R2_cnt = fq_stream_sampler(in_read2_fq_fn, r2_pipe[1], read_frac, seed_val, &R2_org);
        }
        
#pragma omp section 
        {   // Read2 writer
            ss_stream_writer(out_read2_fq_fn, r2_pipe[0], gzip->count);
        }
        
#pragma omp section 
        {   // Single1 reader
            singlet1_cnt = fq_stream_sampler(in_single1_fq_fn, s1_pipe[1], read_frac, seed_val+1, &singlet1_org);
        }
        
#pragma omp section 
        {   // Single1 writer
            ss_stream_writer(out_single1_fq_fn, s1_pipe[0], gzip->count);
        }
        
#pragma omp section 
        {   // Single2 reader
            singlet2_cnt = fq_stream_sampler(in_single2_fq_fn, s2_pipe[1], read_frac, seed_val-1, &singlet2_org);
        }
        
#pragma omp section 
        {   // Single2 writer
            ss_stream_writer(out_single2_fq_fn, s2_pipe[0], gzip->count);
        }
    } 
    
    utstring_free(in_read1_fq_fn);
    utstring_free(in_read2_fq_fn);
    utstring_free(in_single1_fq_fn);
    utstring_free(in_single2_fq_fn);
    utstring_free(out_read1_fq_fn);
    utstring_free(out_read2_fq_fn);
    utstring_free(out_single1_fq_fn);
    utstring_free(out_single2_fq_fn);
    
    if ((R1_org + R2_org) && !(singlet1_cnt + singlet2_cnt)) {
        fprintf(stderr, "\nWarning! read1/read2 files were processed, but no corresponding input singlets were found.\n");
    } 
    
    if ((R1_cnt != R2_cnt) || (R1_org != R2_org)) {
        
        fprintf(stderr, "\nERROR! read1 and read2 fastq files did not contain an equal number of reads!!! %lu %lu\n", R1_org, R2_org);
        exit(EXIT_FAILURE);
        
    } else if (!(mp_org+singlet1_org+singlet2_org)) {

        fprintf(stderr, "Error! No reads found in input files, or input(s) not found.\n");
        exit(EXIT_FAILURE);
        
    } else {
   
        mp_cnt = R1_cnt;
        mp_org = R1_org;
        
        printf("\nMatepairs: Before: %lu, After: %lu\n", mp_org, mp_cnt);
        printf("Singlets: Before: %lu %lu After: %lu %lu\n", singlet1_org, singlet2_org, singlet1_cnt, singlet2_cnt);
        printf("Total Reads Processed: %lu, Reads retained: %lu\n", (2*mp_org)+singlet1_org+singlet2_org, (2*mp_cnt)+singlet1_cnt+singlet2_cnt);
    
        exit(EXIT_SUCCESS);
    }
}

/////////////////////////////////////////////////////////////////
// 
// fq_stream_counter
//
// This function takes an open FASTQ file (gzipped or not) and counts
// the resulting reads.
// This is useful for easily implementing a "read counter thread".
//

unsigned long int fq_stream_counter(UT_string *fq_fn) {

    long unsigned int cnt = 0;
    
    UT_string *data;
    utstring_new(data);
    
    FILE *fq_file = NULL;
    
    if (!(utstring_len(fq_fn))) {
        return(0);
    }
    
    // Try to open the read1 fastq file
    if (!(fq_file = gzopen(utstring_body(fq_fn), "r"))) {
        utstring_printf(fq_fn, ".gz");
        if (!(fq_file = gzopen(utstring_body(fq_fn), "r"))) {
            return(0);    // Not a mate-paired run
        }
    } 
   
    while (ss_gzget_utstring(fq_file, data)) {
        cnt++;
    }

    gzclose(fq_file);
    
    utstring_free(data);
    
    return(cnt/4);
}

/////////////////////////////////////////////////////////////////
// 
// fq_stream_sampler
//
// This function takes an open FASTQ file (gzipped or not) and 
// randomly samples the resulting reads. 
// This is useful for easily implementing a "read sampler thread".
//

unsigned long int fq_stream_sampler(UT_string *fq_fn, int pipe_fd, double sample_frac, unsigned int seed, unsigned long int *org) {
    
    UT_string *data;
    utstring_new(data);
    ss_rand_inst *rnd_state = ss_rseed(seed);
    unsigned long int cnt = 0;
    
    FILE *fq_file = NULL;

    FILE *pipe_in = fdopen(pipe_fd, "w");
    
    if (!(utstring_len(fq_fn))) {
        fclose(pipe_in);
        return(0);
    }
    
    // Make sure count starts at zero
    *org = 0;
    
    // Try to open the read1 fastq file
    if (!(fq_file = gzopen(utstring_body(fq_fn), "r"))) {
        utstring_printf(fq_fn, ".gz");
        if (!(fq_file = gzopen(utstring_body(fq_fn), "r"))) {
            fclose(pipe_in);
            return(0);    // Not a mate-paired run
        }
    } 
    
    while (ss_gzget_utstring(fq_file, data)) {

        (*org)++;
        
        if (ss_rand(rnd_state) <= sample_frac) {
            fputs(utstring_body(data), pipe_in);
            ss_gzget_utstring(fq_file, data);
            fputs(utstring_body(data), pipe_in);
            ss_gzget_utstring(fq_file, data);
            fputs(utstring_body(data), pipe_in);
            ss_gzget_utstring(fq_file, data);
            fputs(utstring_body(data), pipe_in);
            cnt++;
        } else {
            ss_gzget_utstring(fq_file, data);
            ss_gzget_utstring(fq_file, data);
            ss_gzget_utstring(fq_file, data);
        }
    }

    fclose(pipe_in);
    
    gzclose(fq_file);
    
    utstring_free(data);
    free(rnd_state);
    
    return(cnt);
}
