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
 Name        : fastq_nodup
 Description : Read FASTQ format file(s) and reject reads that appear to be
               PCR over-amplification duplicates 
 =============================================================================
*/

#include <seastar_shared.h>

#include <argtable2.h>
#include <uthash.h>
#include <utlist.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/////////////////////////
// Macros 
/////////////////////////

// This macro converts a string representation of read data (i.e. "ACGTN") into 32-bits of 2 bit binary coded numeric data. 
#define CONVERT(s) (c[(int)s[15]]<<30 |  \
                    c[(int)s[14]]<<28 |  \
                    c[(int)s[13]]<<26 |  \
                    c[(int)s[12]]<<24 |  \
                    c[(int)s[11]]<<22 |  \
                    c[(int)s[10]]<<20 |  \
                    c[(int)s[9]]<<18 |   \
                    c[(int)s[8]]<<16 |   \
                    c[(int)s[7]]<<14 |   \
                    c[(int)s[6]]<<12 |   \
                    c[(int)s[5]]<<10 |   \
                    c[(int)s[4]]<<8 |    \
                    c[(int)s[3]]<<6 |    \
                    c[(int)s[2]]<<4 |    \
                    c[(int)s[1]]<<2 |    \
                    c[(int)s[0]])

/////////////////////////
// Typedefs 
/////////////////////////

typedef struct read_list_st {
    unsigned int prefix; 
    unsigned int quality;
    unsigned long int best;
    unsigned int rejected;
    struct read_list_st *next;    // Note, this must be called "next", for utlist.h
} read_list_entry;

typedef struct reject_list_st {
    unsigned long int id;
    UT_hash_handle hh;
} reject_list_entry;

/////////////////////////
// Globals 
/////////////////////////

// LUT for transliteration, note that 'N' becomes 'A'
static const char c[] = 
"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\x00N\x01NNN\x02NNNNNN\x00NNNNN\x03NNNNNNNNNNNN\x00N\x01NNN\x02NNNNNN\x00NNNNN\x03NNNNNNNNNNN";

/////////////////////////
// Function prototypes
/////////////////////////

unsigned long int reject_handler(int reject_pipe_fd, reject_list_entry **rejects);

unsigned long int fq_stream_prefixer(UT_string *fq_fn, int pipe_fd, int pipe2_fd);

unsigned long int match_rejecter(int r1_pipe, int r2_pipe, int reject_pipe, read_list_entry **r_list, int read_len, int slop, unsigned int *eLUT, unsigned int verbose);

unsigned long int fq_stream_rejecter(UT_string *fq_fn, int pipe_fd, UT_string *out_prefix, int no_pre, reject_list_entry *rejects, unsigned long int *rej_cnt);

unsigned long int fq_stream_singlet_rejecter(UT_string *fq_fn, int pipe_fd, UT_string *out_prefix, int no_pre, int read_len, read_list_entry **r_list1, read_list_entry **r_list2, int seed, unsigned long int *rej_cnt);

unsigned long int mismatch_rejecter(int reject_pipe_fd, read_list_entry **r_list, unsigned long int r_start,unsigned long int r_stop, int read_len, int slop, int r_slop, unsigned int *eLUT, unsigned int verbose);

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
    
    reject_list_entry *rejects = NULL;
    
    unsigned long int num_rejects = 0, read1_rejects = 0, read2_rejects = 0, singlet1_rejects = 0, singlet2_rejects = 0;
    
    // Input filenames
    UT_string *in_read1_fq_fn, *in_read2_fq_fn, *in_single1_fq_fn, *in_single2_fq_fn;
    utstring_new(in_read1_fq_fn);
    utstring_new(in_read2_fq_fn);
    utstring_new(in_single1_fq_fn);
    utstring_new(in_single2_fq_fn);
    
    // Output filenames
    UT_string *out_read1_fq_fn, *out_read2_fq_fn, *out_single1_fq_fn, *out_single2_fq_fn, *out_mates_fq_fn;
    utstring_new(out_read1_fq_fn);
    utstring_new(out_read2_fq_fn);
    utstring_new(out_single1_fq_fn);
    utstring_new(out_single2_fq_fn);
    utstring_new(out_mates_fq_fn);
    
    // Read name prefix
    UT_string *out_read_prefix;
    utstring_new(out_read_prefix);
    
    // Parameters
    int read_len = 0, r_slop = 0, slop = 0; 

    // Flags
    // Read counters
    long unsigned int mp_cnt = 0, singlet1_cnt = 0, singlet2_cnt = 0, R1_org = 0, R2_org = 0; 
    int singles_flag = 0;

    // Hash table struct vars
    read_list_entry **r_list1 = NULL;
    read_list_entry **r_list2 = NULL;
    long unsigned int read_list_length = 0;
    
    // Seed for the random number generator
    int rnd_seed = 12345;
    
    // Error LUT
    unsigned int ErrLUT[0x10000];
    
    ////////////////////////////////////////////////////////////////////////
    // All done with variable declarations!!
    
    ///////////////////////////////////
    // Command line argtable settings
    ///////////////////////////////////
    
    struct arg_lit *gzip = arg_lit0("z", "gzip", "Output converted files in gzip compressed format. [NULL]");
    struct arg_int *len = arg_int0("l","index_len","<u>","Match length of both mates [14]");
    struct arg_int *mismatch = arg_int0("d","index_err","<u>","Mismatches within indexed mate [1]");
    struct arg_int *m_mismatch = arg_int0("e","match_err","<u>","Mismatches within unindexed mate [3]");
    struct arg_int *seed = arg_int0(NULL,"seed","<n>","Seed used by random number generator [12345]");
    struct arg_str *pre_read_id = arg_str0(NULL, "prefix", "<string>", "Prefix to add to read identifiers. [out_prefix]");
    struct arg_lit *no_pre = arg_lit0(NULL, "no_prefix", "Do not change the read names in any way. [NULL]");
    struct arg_file *input = arg_file1(NULL, NULL, "<in_prefix>", "Input file prefix: (e.g. <in_prefix>_single.fastq [<in_prefix>_read1.fastq <in_prefix>_read2.fastq]) ");
    struct arg_file *output = arg_file1(NULL, NULL, "<out_prefix>", "Output file prefix: (e.g. <out_prefix>_single.fastq [<out_prefix>_read1.fastq <out_prefix>_read2.fastq]) ");
    struct arg_lit *v = arg_lit0("v", "verbose", "Print lots of information about run.");
    struct arg_lit *version = arg_lit0(NULL,"version","Print the build version and exit.");    
    struct arg_lit *h = arg_lit0("h", "help", "Request help.");
    struct arg_end *end = arg_end(20);
    
    void *argtable[] = {h,version,v,gzip,pre_read_id,no_pre,len,mismatch,m_mismatch,seed,input,output,end};

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
        fprintf(stderr, "\nInput and output \"prefixes\" are the part of the filename before:\n");
        fprintf(stderr, "_single.fastq _read1.fastq _read2.fastq  Mate-paired files\n");
        fprintf(stderr, "(_read1 _read2) are required. A singlet read file is automatically used if present.\n");
        
        exit(EXIT_FAILURE);
    }    
    
    if (arg_errors) { 
        arg_print_errors(stderr, end, "fastqnodup");
        arg_print_syntaxv(stderr, argtable, "\n");
        exit(EXIT_FAILURE);
    }
    
    // Validate read len
    // len should not be set greater than 16 because code in this program
    // assumes that the 2 bit encoded read prefixes can be represented with a 
    // 32 bit unsigned int.  If this size limit is relaxed then care must be
    // taken to change affected 32 unsigned ints to 62 bit unsigned ints.
    if (len->count) {
        read_len = 2*len->ival[0];
        if ((read_len <= 0) || (read_len > 32)) {
            fprintf(stderr, "index_len must be between 1 and 16, inclusive.\n");
            exit(EXIT_FAILURE);
        }        
    } else {
        read_len = 14*2;
    }
    read_list_length = (long unsigned int)1<<read_len;
    
    // Validate read err
    if (mismatch->count) {
        r_slop = mismatch->ival[0];
        if ((r_slop < 0) || (r_slop > 2)) {
            fprintf(stderr, "index_err must be between 0 and 2 inclusive\n");
            exit(EXIT_FAILURE);
        }        
    } else {
        r_slop = 1;
    }
    
    // Validate mate err
    if (m_mismatch->count) {
        slop = m_mismatch->ival[0];
        if (slop < 0) {
            fprintf(stderr, "match_err must be >= 0\n");
            exit(EXIT_FAILURE);
        }        
    } else {
        slop = 3;
    }

    // Seed the pseudo-random number generator (ss_rand)
    
    if (seed->count) {    // negative seed for ss_rand
        rnd_seed = (seed->ival[0] > 0) ? -seed->ival[0] : seed->ival[0];
    }    
    
    // Allocate read index table
    if (!(r_list1 = calloc(read_list_length, sizeof(read_list_entry *)))) {
        fprintf(stderr, "Calloc failed! r_hash allocation\n");
        exit(EXIT_FAILURE);
    } 

    if (!(r_list2 = calloc(read_list_length, sizeof(read_list_entry *)))) {
        fprintf(stderr, "Calloc failed! r_hash allocation\n");
        exit(EXIT_FAILURE);
    } 

    // Construct input filenames
    if (input->count) {
        utstring_printf(in_read1_fq_fn, "%s.read1.fastq", input->filename[0]);
        utstring_printf(in_read2_fq_fn, "%s.read2.fastq", input->filename[0]);
        utstring_printf(in_single1_fq_fn, "%s.single.fastq", input->filename[0]);
    } else {
        fprintf(stderr, "in_prefix not specified\n");
    }
    
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
    
    // Verbose output
    if (v->count) {    
        if (singles_flag) {
            printf("Input files: %s %s %s %s\n", utstring_body(in_read1_fq_fn), utstring_body(in_read2_fq_fn), utstring_body(in_single1_fq_fn), utstring_body(in_single2_fq_fn));
            printf("Output files: %s %s %s %s\n", utstring_body(out_read1_fq_fn), utstring_body(out_read2_fq_fn), utstring_body(out_single1_fq_fn), utstring_body(out_single2_fq_fn));
        } else {
            printf("Input files: %s %s %s\n", utstring_body(in_read1_fq_fn), utstring_body(in_read2_fq_fn), utstring_body(in_single1_fq_fn));
            printf("Output files: %s %s %s\n", utstring_body(out_read1_fq_fn), utstring_body(out_read2_fq_fn), utstring_body(out_single1_fq_fn));
        }
        printf("Read len: %d\n", read_len/2);
        printf("Index mismatches: %d  Mate mismatches: %d \n", r_slop, slop);
    }
    
    if (pre_read_id->count) {
        
        if (no_pre->count) {
            fprintf(stderr, "Error: Both --prefix and --no_prefix were specified.\n");
            exit(EXIT_FAILURE);
        }
        
        if (! strlen(pre_read_id->sval[0])) {
            fprintf(stderr, "Read ID prefix may not be zero length.\n");
            exit(EXIT_FAILURE);
        } 
        
        if (strchr(pre_read_id->sval[0], ':') || strchr(pre_read_id->sval[0], '|') || strchr(pre_read_id->sval[0], '+') || strchr(pre_read_id->sval[0], '/')) {
            fprintf(stderr, "Read ID prefix '%s' may not contain the characters ':', '|', '+' or '/'.\n", pre_read_id->sval[0]);
            exit(EXIT_FAILURE);
        }
        
        // Build default read ID prefix
        ss_strcat_utstring(out_read_prefix, pre_read_id->sval[0]);
        
    } else {
        
        if (!no_pre->count) {
            if (strchr(output->filename[0], ':') || strchr(output->filename[0], '|') || strchr(output->filename[0], '+') || strchr(output->filename[0], '/')) {
                fprintf(stderr, "Read ID prefix '%s' (from output prefix) may not contain the characters ':', '|', '+' or '/'.\n", output->filename[0]);
                fprintf(stderr, "Hint: Use the --prefix or --no_prefix parameter if the output file prefix contains path information.\n");
                exit(EXIT_FAILURE);
            }        
            
            // Build default read ID prefix
            ss_strcat_utstring(out_read_prefix, output->filename[0]);
        }
    }
    
    // Check for null string prefixes
    if (!(strlen(input->filename[0]) && strlen(output->filename[0]))) {
        fprintf(stderr, "Error: NULL prefix strings are not permitted.\n");
        exit(EXIT_FAILURE);
    }
    
    // Begin processing!
    
#ifdef _OPENMP    
    omp_set_num_threads(9);
#endif    
    
    int r1_pipe[2];
    int r2_pipe[2];
    pipe(r1_pipe);
    pipe(r2_pipe);
    int r1_pipe2[2];
    int r2_pipe2[2];
    pipe(r1_pipe2);
    pipe(r2_pipe2);

    int reject_pipe[2];
    pipe(reject_pipe);

    int rej_dup = dup(reject_pipe[1]);

#pragma omp parallel for default(shared)       
    for (int x = 0; x <= 0xffff; x++) {
        unsigned int err = ((x>>2)&0x3333)+(x&0x3333);
        err = ((err>>8)+err)&0xff;
        err = (err>>4)+(err&0xf);
        ErrLUT[x] = err;
    }    
 
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // The following parallel code builds a hash table of duplicate mate-pairs of reads to reject.   
    // These threads only examine matching duplicate index reads
    
#pragma omp parallel sections default(shared)       
    {
        
#pragma omp section 
        {   // Read1 reader
            fq_stream_prefixer(in_read1_fq_fn, r1_pipe[1], r1_pipe2[1]);
        }
        
#pragma omp section 
        {   // Read2 reader
            fq_stream_prefixer(in_read2_fq_fn, r2_pipe[1], r2_pipe2[1]);
        }        

#pragma omp section 
        {   // Read1 processor
            match_rejecter(r1_pipe[0], r2_pipe[0], reject_pipe[1], r_list1, read_len, slop, ErrLUT, v->count);
        }
        
#pragma omp section 
        {   // Read2 processor
            match_rejecter(r2_pipe2[0], r1_pipe2[0], rej_dup, r_list2, read_len, slop, ErrLUT, v->count);
        } 
        
#pragma omp section 
        {   // Reject handler
            num_rejects = reject_handler(reject_pipe[0], &rejects);
        }
    }
   
    // Only do this analysis if the r_slop (-d paramter) is one or more
    if (r_slop > 0) {
        
        pipe(reject_pipe);
        int rej_p[8];
        rej_p[0] = reject_pipe[1];
        for (int x = 1; x < 8; x++) {
            rej_p[x] = dup(rej_p[x-1]);
        }
        
        /////////////////////////////////////////////////////////////////////////////////////////////////
        // The following parallel code extends the hash table of duplicate mate-pairs of reads to reject.
        // These threads examine mis-matching duplicate index reads
        
        
#pragma omp parallel sections default(shared)       
        {
            
#pragma omp section 
            {   // Read1 reader
                mismatch_rejecter(rej_p[0], r_list1, (0*(1L<<(read_len-2))), (1*(1L<<(read_len-2)))-1, read_len, slop, r_slop, ErrLUT, v->count);
            }
            
#pragma omp section 
            {   // Read1 reader
                mismatch_rejecter(rej_p[1], r_list1, (1*(1L<<(read_len-2))), (2*(1L<<(read_len-2)))-1, read_len, slop, r_slop, ErrLUT, v->count);
            }
            
#pragma omp section 
            {   // Read1 reader
                mismatch_rejecter(rej_p[2], r_list1, (2*(1L<<(read_len-2))), (3*(1L<<(read_len-2)))-1, read_len, slop, r_slop, ErrLUT, v->count);
            }
            
#pragma omp section 
            {   // Read1 reader
                mismatch_rejecter(rej_p[3], r_list1, (3*(1L<<(read_len-2))), (4*(1L<<(read_len-2)))-1, read_len, slop, r_slop, ErrLUT, v->count);
            }
            
#pragma omp section 
            {   // Read1 reader
                mismatch_rejecter(rej_p[4], r_list2, (0*(1L<<(read_len-2))), (1*(1L<<(read_len-2)))-1, read_len, slop, r_slop, ErrLUT, v->count);
            }
            
#pragma omp section 
            {   // Read1 reader
                mismatch_rejecter(rej_p[5], r_list2, (1*(1L<<(read_len-2))), (2*(1L<<(read_len-2)))-1, read_len, slop, r_slop, ErrLUT, v->count);
            }
            
#pragma omp section 
            {   // Read1 reader
                mismatch_rejecter(rej_p[6], r_list2, (2*(1L<<(read_len-2))), (3*(1L<<(read_len-2)))-1, read_len, slop, r_slop, ErrLUT, v->count);
            }
            
#pragma omp section 
            {   // Read1 reader
                mismatch_rejecter(rej_p[7], r_list2, (3*(1L<<(read_len-2))), (4*(1L<<(read_len-2)))-1, read_len, slop, r_slop, ErrLUT, v->count);
            }
            
#pragma omp section 
            {   // Reject handler
                num_rejects = reject_handler(reject_pipe[0], &rejects);
            }
        }
    }
    
    int s1_pipe[2];
    int s2_pipe[2];
    pipe(r1_pipe);
    pipe(r2_pipe);
    pipe(s1_pipe);
    pipe(s2_pipe);
   
    /////////////////////////////////////////////////////////////////////////////////////////////////
    // The following parallel code uses the hash table of rejected duplicate mate-pairs to re-read
    // the input files and write outputs without the rejected reads.     
    
#pragma omp parallel sections default(shared)       
    {
        
#pragma omp section 
        {   // Read1 reader
            R1_org = fq_stream_rejecter(in_read1_fq_fn, r1_pipe[1], out_read_prefix, no_pre->count, rejects, &read1_rejects);
        }

#pragma omp section 
        {   // Read1 writer
            ss_stream_writer(out_read1_fq_fn, r1_pipe[0], gzip->count);
        }
        
#pragma omp section 
        {   // Read2 reader
            R2_org = fq_stream_rejecter(in_read2_fq_fn, r2_pipe[1], out_read_prefix, no_pre->count, rejects, &read2_rejects);
        }
        
#pragma omp section 
        {   // Read2 writer
            ss_stream_writer(out_read2_fq_fn, r2_pipe[0], gzip->count);
        }
        
#pragma omp section 
        {   // Single1 reader
            singlet1_cnt = fq_stream_singlet_rejecter(in_single1_fq_fn, s1_pipe[1], out_read_prefix, no_pre->count, read_len, r_list1, r_list2, rnd_seed, &singlet1_rejects);
        }
        
#pragma omp section 
        {   // Single1 writer
            ss_stream_writer(out_single1_fq_fn, s1_pipe[0], gzip->count);
        }
        
#pragma omp section 
        {   // Single2 reader
            singlet2_cnt = fq_stream_singlet_rejecter(in_single2_fq_fn, s2_pipe[1], out_read_prefix, no_pre->count, read_len, r_list1, r_list2, rnd_seed+1, &singlet2_rejects);
        }
        
#pragma omp section 
        {   // Single2 writer
            ss_stream_writer(out_single2_fq_fn, s2_pipe[0], gzip->count);
        }
    }    
    
    if (R1_org != R2_org) {
        fprintf(stderr, "\nERROR! read1 and read2 fastq files did not contain an equal number of reads. %lu %lu\n", R1_org, R2_org);
        exit(EXIT_FAILURE);
    } 
    
    if (!(R1_org+singlet1_cnt+singlet2_cnt)) {
        
        fprintf(stderr, "ERROR! No reads found in input files, or input(s) not found.\n");
        exit(EXIT_FAILURE);
        
    }
    
    if (mp_cnt && !(singlet1_cnt + singlet2_cnt)) {
        fprintf(stderr, "\nWarning! read1/read2 files were processed, but no corresponding input singlets were found.\n");
    } 
    
    num_rejects = read1_rejects;
    mp_cnt = R1_org;

    printf("\nMatepairs: %lu  Rejected: %lu (%0.1f%%)\n", mp_cnt, num_rejects, 100.0*num_rejects/(double)mp_cnt);
    printf("Singlets: %lu %lu  Rejected: %lu %lu (%0.1f%%)\n", singlet1_cnt, singlet2_cnt, 
           singlet1_rejects, singlet2_rejects, 100.0*(singlet1_rejects+singlet2_rejects)/(double)(singlet1_cnt+singlet2_cnt));
    
    // Free up all of the data structures.

#pragma omp parallel for default(shared)       
    for (int x = 0; x < read_list_length; x++) {
        read_list_entry *entry = NULL;
        
        while (r_list1[x]) {
            entry = r_list1[x];
            LL_DELETE(r_list1[x],r_list1[x]);
            free(entry);
        }
        
        while (r_list2[x]) {
            entry = r_list2[x];
            LL_DELETE(r_list2[x],r_list2[x]);
            free(entry);
        }
    }
    
    free(r_list1);
    free(r_list2);
    
    exit(EXIT_SUCCESS);
}

/////////////////////////////////////////////////////////////////
// 
// fq_stream_prefixer
//
// This function takes a FASTQ filename (gzipped or not) and 
// generates bitpacked prefixes and quality sums for the reads.

unsigned long int fq_stream_prefixer(UT_string *fq_fn, int pipe_fd, int pipe2_fd) {

    UT_string *head_data;
    utstring_new(head_data);
    UT_string *seq_data;
    utstring_new(seq_data);
    UT_string *extra_data;
    utstring_new(extra_data);
    UT_string *qual_data;
    utstring_new(qual_data);
    
    unsigned long int cnt = 0;
    
    unsigned int pq[2];
    
    char *q = NULL;
    
    FILE *fq_file = NULL;
    
    FILE *pipe_in = fdopen(pipe_fd, "wb");
    FILE *pipe2_in = fdopen(pipe2_fd, "wb");

    if (!(utstring_len(fq_fn))) {
        fclose(pipe_in);
        fclose(pipe2_in);
        return(0);
    }
    
    // Try to open the fastq file
    if (!(fq_file = gzopen(utstring_body(fq_fn), "r"))) {
        utstring_printf(fq_fn, ".gz");
        if (!(fq_file = gzopen(utstring_body(fq_fn), "r"))) {
            fclose(pipe_in);
            fclose(pipe2_in);
            return(0);    
        }
    } 
    
    while (ss_gzget_utstring(fq_file, head_data)) {
        
        ss_gzget_utstring(fq_file, seq_data);
        ss_gzget_utstring(fq_file, extra_data);
        ss_gzget_utstring(fq_file, qual_data);
                
        pq[0] = CONVERT(utstring_body(seq_data));
        
        pq[1] = 0;
        for (q = utstring_body(qual_data); *q; q++) {
            pq[1] += (unsigned int) (*q - 33);
        }
        
        fwrite(pq, sizeof(unsigned int), 2, pipe_in);
        fwrite(pq, sizeof(unsigned int), 2, pipe2_in);
        
        cnt++;
    }
    
    fclose(pipe_in);
    fclose(pipe2_in);
    
    gzclose(fq_file);
    
    utstring_free(head_data);
    utstring_free(seq_data);
    utstring_free(extra_data);
    utstring_free(qual_data);
    
    return(cnt);
}


//////////////////////////////////////////////////////////////////////////////////////////////
// 
// match_rejecter
//
// This function takes streams of bitpacked prefixes and quality sums for the reads,
// and generates a reject list for perfect matches of the index read (with slop for the mate).

unsigned long int match_rejecter(int r1_pipe_fd, int r2_pipe_fd, int reject_pipe_fd, read_list_entry **r_list, int read_len, int slop, unsigned int *eLUT, unsigned int verbose) {
    
    FILE *r1_in = fdopen(r1_pipe_fd, "rb");
    FILE *r2_in = fdopen(r2_pipe_fd, "rb");
    FILE *rej_out = fdopen(reject_pipe_fd, "wb");
    
    read_list_entry **entry1 = NULL, *entry2 = NULL, *r2p = NULL;
    
    unsigned int pq_r1[2], pq_r2[2];
    
    unsigned int qual = 0, read1_prefix = 0, read2_prefix = 0;
    
    unsigned long int mp_cnt = 0;
    // may be too large for unsigned int so used unsigned long first
    unsigned long int read_mask_tmp = ((long unsigned int)1<<read_len) - 1;
    unsigned int read_mask = read_mask_tmp;
    
    while (fread(pq_r1, sizeof(unsigned int), 2, r1_in)) {
        fread(pq_r2, sizeof(unsigned int), 2, r2_in);
        
        // Prepared read prefixes and quality summary info
        read1_prefix = pq_r1[0] & read_mask;
        read2_prefix = pq_r2[0] & read_mask;
        qual = pq_r1[1] + pq_r2[1];
        
        
        // Look to see if a read1 with this prefix exists
        entry1 = &r_list[read1_prefix];
        entry2 = NULL;
        
        unsigned int err = 0;
        
        // Search for matching read2 entries for the exact read1 prefix 
        LL_FOREACH(*entry1,r2p) {
            err = 0;
            unsigned int xor = read2_prefix^r2p->prefix;
            
            xor = ((xor>>1)|xor)&0x55555555;
            err = eLUT[((xor>>16)+xor)&0xffff];
            
            // If the read2 prefix is within slop mismatches, make it entry2
            if (err <= slop) {
                entry2 = r2p;
                break;
            }
        }
        
        if (!entry2) {    // If no matching read2 prefix was found
            
            // This hasn't changed in this function
            // entry1 = &r_list[read1_prefix];
            
            if (!(entry2 = malloc(sizeof(read_list_entry)))) {  
                fprintf(stderr, "malloc failed: read2 hash entry\n");
                exit(EXIT_FAILURE);        
            }
            
            // Generate a new mate for the read1 prefix
            entry2->prefix = read2_prefix;
            entry2->quality = 0;
            entry2->best = 0;
            entry2->rejected = -1;    // Will be incremented below
            LL_PREPEND(*entry1, entry2);
        }
        
        entry2->rejected += 1;
        
        // If the quality of this read is better...
        if (qual > entry2->quality) {   // For new entry2 will always be true
            
            if (entry2->rejected) { // For new entry2 will never be true 
                fwrite(&(entry2->best), sizeof(unsigned long int), 1, rej_out);
                
                if (verbose) {
                    fprintf(stderr, "REJECT: %lu by %lu (0,%d)\n", entry2->best, mp_cnt, err);
                }
            }
            
            entry2->best = mp_cnt;
            entry2->quality = qual;
            entry2->prefix = read2_prefix;    // Replace it, since this is better qual
            
        } else {
            fwrite(&mp_cnt, sizeof(unsigned long int), 1, rej_out);
            
            if (verbose) {
                fprintf(stderr, "REJECT: %lu by %lu (0,%d) \n", mp_cnt, entry2->best, err);
            }
        }
        
        mp_cnt++;    // Inc the counter
        
    }
    
    fclose(r1_in);
    fclose(r2_in);
    fclose(rej_out);
    
    return(mp_cnt);
}


/////////////////////////////////////////////////////////////////
// 
// mismatch_rejecter
//
// This function takes the data structure built by fq_stream_rejector
// and looks for up to two mismatch duplicates to remove.
//

unsigned long int mismatch_rejecter(int reject_pipe_fd, read_list_entry **r_list, unsigned long int r_start, unsigned long int r_stop, int read_len, int slop, int r_slop, unsigned int *eLUT, unsigned int verbose) {
    
    FILE *rej_out = fdopen(reject_pipe_fd, "wb");
    
    read_list_entry **entry1 = NULL, *entry2 = NULL, *r2p = NULL, **read1 = NULL, *read2 = NULL;
    
    unsigned long int r1p = 0, read1_prefix = 0;
    
    unsigned long int mp_cnt = 0;
    
    if (r_slop == 0) {
        fprintf(stderr, "ASSERT: Function mismatch_rejecter should not be called with r_slop == 0.\n");
        exit(EXIT_FAILURE);
    }   
    
    for (read1_prefix = r_start; read1_prefix <= r_stop; read1_prefix++) {
        
        // Lookup the entry for this prefix
        read1 = &r_list[read1_prefix];
        
        // Search for matching read2 entries for read1
        LL_FOREACH(*read1,read2) {
            
            unsigned int err = 0;
            unsigned int miss = 1;
            
            entry2 = NULL;
            
            // Mismatch 1, try all 3*read_len possibilities
            for (int i = 0; i < read_len; i+=2) {
                unsigned int rptmp = (read1_prefix & ~(3 << i));
                
                unsigned int same_a = (read1_prefix >> i) & 3;
                for (int a = 0; a < 4; a++) {
                    if (a != same_a) {    // Don't repeat the exact match!
                        r1p = rptmp | (a << i);
                        entry1 = &r_list[r1p];
                        LL_FOREACH(*entry1,r2p) {
                            err = 0;
                            unsigned int xor = read2->prefix^r2p->prefix;
                            
                            xor = ((xor>>1)|xor)&0x55555555;
                            err = eLUT[((xor>>16)+xor)&0xffff];
                            
                            if (err <= slop) {
                                entry2 = r2p;
                                goto entry2_found;
                            }
                        }
                    }
                }
            }
            
            miss = 2;
            
            // If no match has been found, then try 2 mismatches
            if (r_slop > 1) {
                // Mismatch 2, try all 3*3*readlen*(read_len-1) possibilities
                for (int i = 0; i < read_len-2; i+=2) {
                    for (int j = i+2; j < read_len; j+=2) {
                        unsigned int rptmp = (read1_prefix & ~((3 << i) | (3 << j)));
                        unsigned int same_a = (read1_prefix >> i) & 3;
                        unsigned int same_b = (read1_prefix >> j) & 3;
                        
                        for (int x = 0; x < 16; x++) {
                            unsigned int a = x&3;
                            unsigned int b = x>>2;
                            if ((a != same_a) && (b != same_b)) {
                                r1p = rptmp | (a << i) | (b << j);
                                entry1 = &r_list[r1p];
                                LL_FOREACH(*entry1, r2p) {
                                    
                                    err = 0;
                                    unsigned int xor = read2->prefix^r2p->prefix;
                                    
                                    xor = ((xor>>1)|xor)&0x55555555;
                                    err = eLUT[((xor>>16)+xor)&0xffff];
                                    
                                    if (err <= slop) {
                                        entry2 = r2p;
                                        goto entry2_found;
                                    }
                                }
                            }
                        }    
                    } 
                }
            }        
            
        entry2_found:
            
            if (entry2) {
                if (read2->quality > entry2->quality) {
                    fwrite(&(entry2->best), sizeof(unsigned long int), 1, rej_out);
                    read2->rejected += 1;
                    if (verbose) {
                        fprintf(stderr, "REJECT: %lu by %lu (%d,%d) \n", entry2->best, read2->best, miss, err);
                    }
                    
                } else if (read2->quality < entry2->quality) {
                    fwrite(&(read2->best), sizeof(unsigned long int), 1, rej_out);
                    entry2->rejected += 1;
                    if (verbose) {
                        fprintf(stderr, "REJECT: %lu by %lu (%d,%d)\n", read2->best, entry2->best, miss, err);
                    }
                    
                } else if (read2->best < entry2->best) {
                    fwrite(&(entry2->best), sizeof(unsigned long int), 1, rej_out);
                    read2->rejected += 1;                
                    if (verbose) {
                        fprintf(stderr, "REJECT: %lu by %lu (%d,%d) TIE!\n", entry2->best, read2->best, miss, err);
                    }
                    
                } else if (read2->best > entry2->best) {
                    fwrite(&(read2->best), sizeof(unsigned long int), 1, rej_out);
                    entry2->rejected += 1;
                    if (verbose) {
                        fprintf(stderr, "REJECT: %lu by %lu (%d,%d) TIE!\n", read2->best, entry2->best, miss, err);
                    }
                    
                } else {  // Should never get here!
                    fprintf(stderr, "ASSERT: Duplicate read numbers encountered in reject code!\n");
                    exit(EXIT_FAILURE);    
                }
            }
            
            mp_cnt++;    // Inc the counter
            
        }
        
    }    
    
    fclose(rej_out);
    
    return(mp_cnt);
}


/////////////////////////////////////////////////////////////////
// 
// fq_stream_rejecter
//
// This function takes an open FASTQ file (gzipped or not) and 
// selects/rejects the resulting reads. 
//

unsigned long int fq_stream_rejecter(UT_string *fq_fn, int pipe_fd, UT_string *out_prefix, int no_pre, reject_list_entry *rejects, unsigned long int *rej_cnt) {
    
    UT_string *new_head_data;
    utstring_new(new_head_data);
    UT_string *head_data;
    utstring_new(head_data);
    UT_string *seq_data;
    utstring_new(seq_data);
    UT_string *extra_data;
    utstring_new(extra_data);
    UT_string *qual_data;
    utstring_new(qual_data);
    
    unsigned long int mp_cnt = 0;
    
    reject_list_entry *rej = NULL;
    
    char *start = NULL; 
    char *end = NULL;
    
    FILE *fq_file = NULL;
    
    FILE *pipe_in = fdopen(pipe_fd, "w");
    
    if (!(utstring_len(fq_fn))) {
        fclose(pipe_in);
        return(0);
    }
    
    // Try to open the read1 fastq file
    if (!(fq_file = gzopen(utstring_body(fq_fn), "r"))) {
        utstring_printf(fq_fn, ".gz");
        if (!(fq_file = gzopen(utstring_body(fq_fn), "r"))) {
            fclose(pipe_in);
            return(0);    
        }
    } 

    while (ss_gzget_utstring(fq_file, head_data)) {
        
        ss_gzget_utstring(fq_file, seq_data);
        ss_gzget_utstring(fq_file, extra_data);
        ss_gzget_utstring(fq_file, qual_data);
        
        HASH_FIND(hh, rejects, &mp_cnt, sizeof(unsigned long int), rej);
        
        if (!rej) {    // If not rejected, then write 'em out!
            
            // Fixup the read name
            utstring_clear(new_head_data);
            
            end = strchr(utstring_body(head_data), ':');
            
            if (!end || no_pre) {
                utstring_concat(new_head_data, head_data);
            } else {
                end++;
                if ((start = strchr(utstring_body(head_data), '|'))) {
                    *start = '\0';
                    utstring_printf(new_head_data, "@%s|%s:%s",utstring_body(head_data)+1,utstring_body(out_prefix),end);
                } else {
                    utstring_printf(new_head_data, "@%.2s+%s:%s",utstring_body(head_data)+1,utstring_body(out_prefix),end);
                } 
            }
            
            fputs(utstring_body(new_head_data), pipe_in);
            fputs(utstring_body(seq_data), pipe_in);
            fputs(utstring_body(extra_data), pipe_in);
            fputs(utstring_body(qual_data), pipe_in);
            
        } else {
            (*rej_cnt)++;
        }

        mp_cnt++;
    }
        
    fclose(pipe_in);
    
    gzclose(fq_file);
    
    utstring_free(new_head_data);
    utstring_free(head_data);
    utstring_free(seq_data);
    utstring_free(extra_data);
    utstring_free(qual_data);
    
    return(mp_cnt);
}


/////////////////////////////////////////////////////////////////
// 
// fq_stream_singlet_rejecter
//
// This function takes an open FASTQ file (gzipped or not) and 
// selects/rejects the resulting reads. 
//

unsigned long int fq_stream_singlet_rejecter(UT_string *fq_fn, int pipe_fd, UT_string *out_prefix, int no_pre, int read_len, read_list_entry **r_list1, read_list_entry **r_list2, int seed, unsigned long int *rej_cnt) {
    
    UT_string *new_head_data;
    utstring_new(new_head_data);
    UT_string *head_data;
    utstring_new(head_data);
    UT_string *seq_data;
    utstring_new(seq_data);
    UT_string *extra_data;
    utstring_new(extra_data);
    UT_string *qual_data;
    utstring_new(qual_data);
    
    unsigned long int cnt = 0;
    // may be too large for unsigned int so used unsigned long first
    unsigned long int read_mask_tmp = ((long unsigned int)1<<read_len) - 1;
    unsigned int read_mask = read_mask_tmp;
    
    unsigned int read_prefix = 0; 

    int rej = 0;
    
    read_list_entry **entry1 = NULL, *entry2 = NULL;
    
    char *start = NULL; 
    char *end = NULL;
    
    ss_rand_inst *rnd_state = NULL;
    rnd_state = ss_rseed(seed);
    
    FILE *fq_file = NULL;
    
    FILE *pipe_in = fdopen(pipe_fd, "w");
    
    if (!(utstring_len(fq_fn))) {
        fclose(pipe_in);
        return(0);
    }
    
    // Try to open the read1 fastq file
    if (!(fq_file = gzopen(utstring_body(fq_fn), "r"))) {
        utstring_printf(fq_fn, ".gz");
        if (!(fq_file = gzopen(utstring_body(fq_fn), "r"))) {
            fclose(pipe_in);
            return(0);    
        }
    } 
    
    while (ss_gzget_utstring(fq_file, head_data)) {
        
        ss_gzget_utstring(fq_file, seq_data);
        ss_gzget_utstring(fq_file, extra_data);
        ss_gzget_utstring(fq_file, qual_data);
        

        read_prefix = CONVERT(utstring_body(seq_data)) & read_mask;
        // Look to see if a read1 with this prefix exists
        entry1 = &r_list1[read_prefix];
        if (!*entry1) {
            entry1 = &r_list2[read_prefix];
        }
    
        rej = 0;
        
        if (*entry1) {
            // Figure out the reject rate for mates of this read1 prefix
            unsigned int rejc = 0, cnt = 0;
            LL_FOREACH(*entry1,entry2) {
                cnt += entry2->rejected + 1;
                rejc += entry2->rejected;
            }
            
            // Choose retain or reject randomly, biased by mate reject rate
            rej = (ss_rand(rnd_state) <= (double)rejc/(double)cnt);
            
            *rej_cnt += rej;
        } 
        
        if (!rej) {    // If not rejected, then write 'em out!
            
            // Fixup the read name
            utstring_clear(new_head_data);
            
            end = strchr(utstring_body(head_data), ':');
            
            if (!end || no_pre) {
                utstring_concat(new_head_data, head_data);
            } else {
                end++;
                if ((start = strchr(utstring_body(head_data), '|'))) {
                    *start = '\0';
                    utstring_printf(new_head_data, "@%s|%s:%s",utstring_body(head_data)+1,utstring_body(out_prefix),end);
                } else {
                    utstring_printf(new_head_data, "@%.2s+%s:%s",utstring_body(head_data)+1,utstring_body(out_prefix),end);
                } 
            }
            
            fputs(utstring_body(new_head_data), pipe_in);
            fputs(utstring_body(seq_data), pipe_in);
            fputs(utstring_body(extra_data), pipe_in);
            fputs(utstring_body(qual_data), pipe_in);
        }
        cnt++;
    }
    
    fclose(pipe_in);
    
    gzclose(fq_file);
    
    utstring_free(new_head_data);
    utstring_free(head_data);
    utstring_free(seq_data);
    utstring_free(extra_data);
    utstring_free(qual_data);
    
    return(cnt);
}

/////////////////////////////////////////////////////
// reject_handler
//

unsigned long int reject_handler(int reject_pipe_fd, reject_list_entry **rejects) {

    unsigned long int cnt = 0;
    
    reject_list_entry *rej = NULL;
    
    FILE *rej_pipe = fdopen(reject_pipe_fd,"rb");
    
    long unsigned int num = 0;
    
    while (fread(&num, sizeof(long unsigned int), 1, rej_pipe)) {
        
        HASH_FIND(hh, *rejects, &num, sizeof(unsigned long int), rej); 
        if (!rej) {
            if (!(rej = malloc(sizeof(reject_list_entry)))) {  
                fprintf(stderr, "malloc failed: reject hash entry\n");
                exit(EXIT_FAILURE);        
            }
            rej->id = num;
            HASH_ADD(hh, *rejects, id, sizeof(unsigned long int), rej);
            cnt++;
        }
        
    }
    
    fclose(rej_pipe);

    return cnt;
}
