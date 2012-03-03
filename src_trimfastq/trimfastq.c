/*
 --------------------------------------------------------------------------- #
 Center for Environmental Genomics
 Copyright (C) 2009-2012 University of Washington.
 
 Authors:
 Vaughn Iverson
 vsi@uw.edu
 
 Chris Berthiaume
 chrisbee@uw.edu
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
 Name        : trimfastq
 Description : Read FASTQ format file(s) and trim/reject reads on quality
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
// Globals 
/////////////////////////

// Probability parmameter
double err_prob = 1.0;

// Entropy parameter
double entropy_cutoff = -1.0;

// Strict entropy
int strict_ent = 0;

// Trimming Method parameter [0=>cs, 1=>nt]
int method_flag = 0;

// Minimum output read length
int min_len = 0;
int min_mate_len = 0;
int fix_len = 0;

/////////////////////////
// Function prototypes
/////////////////////////

int trimmer(char *qualdata);

double entropy_calc(char *seqdata, int len);

unsigned long int fq_stream_trimmer(UT_string *fq_fn, int pipe_fd, UT_string *out_prefix, int no_pre, int len_pre, unsigned long int *comp_cnt, unsigned long int *org, char split);

/////////////////////////////////////////////////////////////////
// Entry point

int main(int argc, char *argv[]) {
    
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
    UT_string *out_read1_fq_fn, *out_read2_fq_fn, *out_single1_fq_fn, *out_single2_fq_fn, *out_mates_fq_fn;
    utstring_new(out_read1_fq_fn);
    utstring_new(out_read2_fq_fn);
    utstring_new(out_single1_fq_fn);
    utstring_new(out_single2_fq_fn);
    utstring_new(out_mates_fq_fn);
    
    // Read name prefix
    UT_string *out_read_prefix;
    utstring_new(out_read_prefix);
    
    // Flags
    int singles_flag = 0;   // 1 when two output singles files being written
    int num_input_singles_files = 0;
    
    // Read counters
    
    unsigned long int mp_org = 0, R1_org = 0, R2_org = 0, singlet1_org = 0, singlet2_org = 0;
    unsigned long int mp_cnt = 0, R1_cnt = 0, R2_cnt = 0, singlet1_cnt = 0, singlet2_cnt = 0, s1_cnt = 0, s2_cnt = 0;
    unsigned long int comp_r1 = 0, comp_r2 = 0, comp_s1 = 0, comp_s2 = 0;
    unsigned long int read1_singlet_cnt = 0, read2_singlet_cnt = 0;

    ////////////////////////////////////////////////////////////////////////
    // All done with variable declarations!!
    
    ///////////////////////////////////
    // Command line argtable settings
    ///////////////////////////////////
    
    struct arg_lit *gzip = arg_lit0("z", "gzip", "Output converted files in gzip compressed format. [NULL]");
    struct arg_lit *inv_singles = arg_lit0("v", "invert_singles", "Causes singles output to be the inverse of the input. 2->1 or 1->2 [NULL]");
    struct arg_lit *num_singles = arg_lit0("s", "singles", "Write two singlet files, one for each mate-paired input file. [NULL]");
    struct arg_rem *sing_rem = arg_rem(NULL, "Note! -v is only valid when there are input singlet reads. -s is only valid when there are NO input singlet reads.");    
    struct arg_lit *colorspace = arg_lit0("c", "color_space", "Use color-space errors to trim and reject reads [default: use nucleotide space errors]");
    struct arg_str *pre_read_id = arg_str0(NULL, "prefix", "<string>", "Prefix to add to read identifiers. [out_prefix]");
    struct arg_lit *no_pre = arg_lit0(NULL, "no_prefix", "Do not change the read names in any way. [NULL]");
    struct arg_lit *pre_read_len = arg_lit0(NULL, "add_len", "Add the final trimmed length value to the read id prefix. [length not added]");
    struct arg_dbl *prob = arg_dbl0("p","correct_prob","<d>","Probability that output reads are correct in nucleotide space (or colorspace with -c) [0.5]");
    struct arg_int *fixed_len = arg_int0("f","fixed_len","<u>","Trim all reads to a fixed length, still filtering on quality [no fixed length]");
    struct arg_int *len = arg_int0("l","min_read_len","<u>","Minimum length of a singlet or longest-mate in nucleotides [24]");
    struct arg_int *mate_len = arg_int0("m","min_mate_len","<u>","Minimum length of the shortest mate in nucleotides [min_read_len]");
    struct arg_dbl *entropy = arg_dbl0("e","entropy_filter","<d>","Remove reads with per position information below given value (in bits per dinucleotide) [No filter]");
    struct arg_lit *entropy_strict = arg_lit0(NULL, "entropy_strict", "Reject reads for low entropy overall, not just the retained part after trimming [NULL]");
    struct arg_lit *mates = arg_lit0(NULL, "mates_file", "Produce a Velvet compatible mate-paired output file (e.g. <out_prefix>_mates.fastq) with read2 mates reversed. [none]");
    struct arg_file *input = arg_file1(NULL, NULL, "<in_prefix>", "Input file prefix: (e.g. <in_prefix>_single.fastq [<in_prefix>_read1.fastq <in_prefix>_read2.fastq]) ");
    struct arg_file *output = arg_file1(NULL, NULL, "<out_prefix>", "Output file prefix: (e.g. <out_prefix>_single.fastq [<out_prefix>_read1.fastq <out_prefix>_read2.fastq]) ");
    struct arg_lit *version = arg_lit0(NULL,"version","Print the build version and exit."); 
    struct arg_lit *h = arg_lit0("h", "help", "Request help.");
    struct arg_end *end = arg_end(20);
    
    void *argtable[] = {h,version,colorspace,gzip,inv_singles,num_singles,sing_rem,prob,len,mate_len,fixed_len,pre_read_id,pre_read_len,no_pre,entropy,entropy_strict,mates,input,output,end};
    int arg_errors = 0;
        
    ////////////////////////////////////////////////////////////////////////
    // Handle command line processing (via argtable2 library) 
    ////////////////////////////////////////////////////////////////////////
    
    arg_errors = arg_parse(argc, argv, argtable);

	if (version->count) {
		fprintf(stderr, "%s version: %s\n", argv[0], SS_BUILD_VERSION);
		exit(EXIT_FAILURE);
    }    
	
    if (h->count) {
        arg_print_syntaxv(stderr, argtable, "\n\n");
        arg_print_glossary(stderr, argtable, "%-25s %s\n");
        fprintf(stderr, "\nInput and output \"prefixes\" are the part of the filename before:\n");
        fprintf(stderr, "_single.fastq [_read1.fastq _read2.fastq]  A singlets (single) file\n");
        fprintf(stderr, "is required.  Mate-paired read files are automatically used if present.\n");
        fprintf(stderr, "Three fastq output files only produced for mate-paired inputs.\n");
        // fprintf(stderr, "\nNote! Input and output fastq files may be gzipped.\n");
        
        exit(EXIT_FAILURE);
    }    
    
    if (arg_errors) { 
        arg_print_errors(stderr, end, "trimfastq");
        arg_print_syntaxv(stderr, argtable, "\n");
        exit(EXIT_FAILURE);
    }
    
    // Validate method
    if (colorspace->count) {
        method_flag = 0;    // cs space
    } else {
        method_flag = 1;    // nt space
    }    

    // Validate entropy
    if (entropy->count) {
        entropy_cutoff = entropy->dval[0];
        if ((entropy_cutoff < 0.0) || (entropy_cutoff > 4.0)) {
            fprintf(stderr, "entropy_filter must be [0.0 - 4.0] \n");
            exit(EXIT_FAILURE);
        }
        strict_ent = entropy_strict->count;
    } else {
        if (entropy_strict->count) {
            fprintf(stderr, "Error: --entropy_strict requires --entropy_filter.\n");
            exit(EXIT_FAILURE);
        } 
        entropy_cutoff = -1.0;
    }    
    
    // Validate error_prob
    if (prob->count) {
        err_prob = prob->dval[0];
        if ((err_prob < 0.0) || (err_prob > 1.0)) {
            fprintf(stderr, "error_prob must be 0.0 - 1.0 inclusive\n");
            exit(EXIT_FAILURE);
        }
    } else {
        err_prob = 0.5;
    }    

    // Validate min read len
    if (len->count) {
        min_len = len->ival[0];
        if (min_len <= 0) {
            fprintf(stderr, "min_read_len must be > 0\n");
            exit(EXIT_FAILURE);
        }        
    } else {
        min_len = 24;
    }

    // Validate min mate len
    if (mate_len->count) {
        min_mate_len = mate_len->ival[0];
        if (min_mate_len <= 0) {
            fprintf(stderr, "min_mate_len must be > 0\n");
            exit(EXIT_FAILURE);
        }        
        if (min_mate_len > min_len) {
            fprintf(stderr, "min_mate_len must be <= min_len\n");
            exit(EXIT_FAILURE);
        }        
    } else {
        min_mate_len = min_len;
    }
    
    if (fixed_len->count) {
        fix_len = min_mate_len = min_len = fixed_len->ival[0];
        if ((mate_len->count) || (len->count)) {
            fprintf(stderr, "fixed_len cannot be used with min_read_len or min_mate_len\n");
            exit(EXIT_FAILURE);
        }
        if (fix_len <= 0) {
            fprintf(stderr, "fixed_len must be > 0\n");
            exit(EXIT_FAILURE);
        }
    } else {
        fix_len = 0;
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
                fprintf(stderr, "Hint: Use the --prefix parameter if the output file prefix contains path information.\n");
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
    
    // Construct input filenames
    utstring_printf(in_read1_fq_fn, "%s.read1.fastq", input->filename[0]);
    utstring_printf(in_read2_fq_fn, "%s.read2.fastq", input->filename[0]);
    
    utstring_printf(in_single1_fq_fn, "%s.single.fastq", input->filename[0]);
    
    FILE *in_single1_file = NULL;
    
    num_input_singles_files = 1;
    
    // Try to open a singlet fastq file
    if (!(in_single1_file = gzopen(utstring_body(in_single1_fq_fn), "r"))) {
        utstring_printf(in_single1_fq_fn, ".gz");
        
        if (!(in_single1_file = gzopen(utstring_body(in_single1_fq_fn), "r"))) {
            utstring_clear(in_single1_fq_fn);
            utstring_printf(in_single1_fq_fn, "%s.single1.fastq", input->filename[0]);
            utstring_printf(in_single2_fq_fn, "%s.single2.fastq", input->filename[0]);
            num_input_singles_files = 2;
            
            if ((in_single1_file = gzopen(utstring_body(in_single1_fq_fn), "r")) || (in_single1_file = gzopen(utstring_body(in_single2_fq_fn), "r"))) {
                singles_flag = 1;   // Two singlet outputs
            } else {
                utstring_printf(in_single1_fq_fn, ".gz");
                utstring_printf(in_single2_fq_fn, ".gz");
                
                if ((in_single1_file = gzopen(utstring_body(in_single1_fq_fn), "r")) || (in_single1_file = gzopen(utstring_body(in_single2_fq_fn), "r"))) {
                    singles_flag = 1;   // Two singlet outputs
                } else {
                    singles_flag = num_singles->count;  // Number of singlet outputs set by -s parm
                    if (inv_singles->count) {
                        fprintf(stderr, "Error: Invalid option -v, No input singlet file(s) found. Use -s to select multiple output singlet files.\n");
                        exit(EXIT_FAILURE);
                    }
                }
            }
        } 
    }
    
    if (in_single1_file) {
        gzclose(in_single1_file);
        if (num_singles->count) {
            fprintf(stderr, "Error: Invalid option -s, Input singlet file(s) found, use -v to change the number of output singlet files.\n");
            exit(EXIT_FAILURE);
        }
    }

    // singles->count inverts the current singles file input scheme
    singles_flag = (singles_flag ^ inv_singles->count);
    
    // Output filenames
    
    utstring_printf(out_read1_fq_fn, "%s.read1.fastq", output->filename[0]);
    utstring_printf(out_read2_fq_fn, "%s.read2.fastq", output->filename[0]);
    
    if (singles_flag == 1) {
        utstring_printf(out_single1_fq_fn, "%s.single1.fastq", output->filename[0]);
        utstring_printf(out_single2_fq_fn, "%s.single2.fastq", output->filename[0]);
    } else {
        utstring_printf(out_single1_fq_fn, "%s.single.fastq", output->filename[0]);
    }

    if (mates->count) {
        utstring_printf(out_mates_fq_fn, "%s.mates.fastq", output->filename[0]);
    }
    
    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Begin processing!
    
#ifdef _OPENMP    
    omp_set_num_threads(10);
#endif    
    
    // This is the value of a non-valid pipe descriptor    
#define NO_PIPE 0
    
    int r1_pipe[2];
    int r2_pipe[2];
    int s1_pipe[2];
    int s2_pipe[2];
    pipe(r1_pipe);
    pipe(r2_pipe);
    pipe(s1_pipe);
    pipe(s2_pipe);
    
    int r1_out_pipe[2];
    int r2_out_pipe[2];
    int mates_out_pipe[2];
    int s1_out_pipe[2];
    int s2_out_pipe[2];
    pipe(r1_out_pipe);
    pipe(r2_out_pipe);
    pipe(mates_out_pipe);
    pipe(s1_out_pipe);
    pipe(s2_out_pipe);
    
    
#pragma omp parallel sections default(shared)       
    {
        
#pragma omp section 
        {   // Read1 reader
            fq_stream_trimmer(in_read1_fq_fn, r1_pipe[1], out_read_prefix, no_pre->count, pre_read_len->count, &comp_r1, &R1_org, '\0');
        }
        
#pragma omp section 
        {   // Read1 writer
            R1_cnt = ss_stream_writer(out_read1_fq_fn, r1_out_pipe[0], gzip->count) / 4;
        }
        
#pragma omp section 
        {   // Read2 reader
            fq_stream_trimmer(in_read2_fq_fn, r2_pipe[1], out_read_prefix, no_pre->count, pre_read_len->count, &comp_r2, &R2_org, '\0');
        }
        
#pragma omp section 
        {   // Read2 writer
            R2_cnt = ss_stream_writer(out_read2_fq_fn, r2_out_pipe[0], gzip->count) / 4;
        }
        
#pragma omp section 
        {   // Single1 reader
            
            // When there is only one input singles file, but two output singles files, then supply which mate to use for this stream in the split parameter
            if ((singles_flag) && (num_input_singles_files == 1)) {
                singlet1_cnt = fq_stream_trimmer(in_single1_fq_fn, s1_pipe[1], out_read_prefix, no_pre->count, pre_read_len->count, &comp_s1, &singlet1_org, '1');
            } else {
                singlet1_cnt = fq_stream_trimmer(in_single1_fq_fn, s1_pipe[1], out_read_prefix, no_pre->count, pre_read_len->count, &comp_s1, &singlet1_org, '\0');
            }
        }
        
#pragma omp section 
        {   // Single1 writer
            s1_cnt = ss_stream_writer(out_single1_fq_fn, s1_out_pipe[0], gzip->count) / 4;
        }
        
#pragma omp section 
        {   // Single2 reader
            
            // When there is only one input singles file, but two output singles files, then supply which mate to use for this stream in the split parameter
            if ((singles_flag) && (num_input_singles_files == 1)) {
                singlet2_cnt = fq_stream_trimmer(in_single1_fq_fn, s2_pipe[1], out_read_prefix, no_pre->count, pre_read_len->count, &comp_s2, &singlet2_org, '2');
            } else {
                singlet2_cnt = fq_stream_trimmer(in_single2_fq_fn, s2_pipe[1], out_read_prefix, no_pre->count, pre_read_len->count, &comp_s2, &singlet2_org, '\0');
            }
        }
        
#pragma omp section 
        {   // Single2 writer
            s2_cnt = ss_stream_writer(out_single2_fq_fn, s2_out_pipe[0], gzip->count) / 4;
        }

#pragma omp section 
        {   // Velvet mates writer
            ss_stream_writer(out_mates_fq_fn, mates_out_pipe[0], gzip->count);
        }
        
#pragma omp section 
        {   // Dispatcher

            
            // Allocate data buffer strings
            
            UT_string *r1_data;
            utstring_new(r1_data);
            UT_string *r2_data;
            utstring_new(r2_data);
            UT_string *s1_data;
            utstring_new(s1_data);
            UT_string *s2_data;
            utstring_new(s2_data);
            
            UT_string *mate_data;
            utstring_new(mate_data);

            UT_string *rev_data;
            utstring_new(rev_data);

            // Read header scanf pattern string:
            // Parses out the three bead coordinates from the first FASTQ "header" line for a read
            // e.g.  @TA+lambda:1231_1682_1381/1
            //                  xxxx yyyy zzzz    where xxxx, yyyy and zzzz are the bead coordinates
            
            const char *head_scanner = "@%*[^:]:%hu_%hu_%hu"; 
            
            // Encoded header 64-bit values 
            // Each solid read coordiante is xxxx_yyyy_zzzz, the code below packs these three
            // (< 16-bit) unsigned values into a 64-bit unsigned integer as:
            // 0000xxxxyyyyzzzz in hex notation
            long unsigned int single1_head_int = 0, single2_head_int = 0, read1_head_int = 0, read2_head_int = 0, last_head_int = 0;
            
            // Pointers to short int parts of the above 64-bit ints (used by sscanf below)
            short unsigned int *single1_1 = (short unsigned int *)(&single1_head_int)+2;
            short unsigned int *single1_2 = (short unsigned int *)(&single1_head_int)+1;
            short unsigned int *single1_3 = (short unsigned int *)(&single1_head_int);
            short unsigned int *single2_1 = (short unsigned int *)(&single2_head_int)+2;
            short unsigned int *single2_2 = (short unsigned int *)(&single2_head_int)+1;
            short unsigned int *single2_3 = (short unsigned int *)(&single2_head_int);
            short unsigned int *read1_1 = (short unsigned int *)(&read1_head_int)+2;
            short unsigned int *read1_2 = (short unsigned int *)(&read1_head_int)+1;
            short unsigned int *read1_3 = (short unsigned int *)(&read1_head_int);
            short unsigned int *read2_1 = (short unsigned int *)(&read2_head_int)+2;
            short unsigned int *read2_2 = (short unsigned int *)(&read2_head_int)+1;
            short unsigned int *read2_3 = (short unsigned int *)(&read2_head_int);

            // Pipes
            FILE *r1_in = fdopen(r1_pipe[0],"r"); 
            FILE *r2_in = fdopen(r2_pipe[0],"r");                 
            
            FILE *s1_in = fdopen(s1_pipe[0],"r");
            FILE *s2_in = fdopen(s2_pipe[0],"r");             

            FILE *mates_out = fdopen(mates_out_pipe[1],"w"); 
            
            FILE *r1_out = fdopen(r1_out_pipe[1],"w"); 
            FILE *r2_out = fdopen(r2_out_pipe[1],"w");                 
            
            FILE *s1_out = fdopen(s1_out_pipe[1],"w");
            FILE *s2_out = fdopen(s2_out_pipe[1],"w");             
            
            if (!singles_flag) {
                fclose(s2_out);
                s2_out = s1_out;
            }
            
            // Prime flags for loop below
            int single1_hungry = 0;
            int single2_hungry = 0;
            int read1_hungry = 0;
            int read2_hungry = 0;
           
            // Prime buffers
            if (ss_get_utstring(r1_in, r1_data)) {
                sscanf(utstring_body(r1_data), head_scanner, read1_1, read1_2, read1_3);
                if ((mates->count) && (read1_head_int >= read2_head_int)) {
                    utstring_clear(mate_data);
                    utstring_concat(mate_data, r1_data);    // head
                    ss_get_utstring(r1_in, rev_data);       // seq
                    utstring_concat(r1_data, rev_data);
                    ss_trunc_utstring(rev_data, 1);  // remove newline
                    ss_rev_utstring(rev_data);
                    ss_strcat_utstring(rev_data, "\n"); // add newline back
                    utstring_concat(mate_data, rev_data);
                    ss_get_utstring(r1_in, rev_data);       // extra
                    utstring_concat(r1_data, rev_data);
                    utstring_concat(mate_data, rev_data);
                    ss_get_utstring(r1_in, rev_data);       // qual
                    utstring_concat(r1_data, rev_data);
                    ss_trunc_utstring(rev_data, 1);  // remove newline
                    ss_rev_utstring(rev_data);
                    ss_strcat_utstring(rev_data, "\n"); // add newline back
                    utstring_concat(mate_data, rev_data);
                } else {
                    ss_get_cat_utstring(r1_in, r1_data);
                    ss_get_cat_utstring(r1_in, r1_data);
                    ss_get_cat_utstring(r1_in, r1_data);
                }
            } else {
                read1_head_int = 0xffffffffffffffff;
            }
            
            if (ss_get_utstring(r2_in, r2_data)) {
                sscanf(utstring_body(r2_data), head_scanner, read2_1, read2_2, read2_3);
                    ss_get_cat_utstring(r2_in, r2_data);
                    ss_get_cat_utstring(r2_in, r2_data);
                    ss_get_cat_utstring(r2_in, r2_data);
            } else {
                read2_head_int = 0xffffffffffffffff;
            }

            if (ss_get_utstring(s1_in, s1_data)) {
                sscanf(utstring_body(s1_data), head_scanner, single1_1, single1_2, single1_3);
                ss_get_cat_utstring(s1_in, s1_data);
                ss_get_cat_utstring(s1_in, s1_data);
                ss_get_cat_utstring(s1_in, s1_data);
            } else {
                single1_head_int = 0xffffffffffffffff;
            }
            
            if (ss_get_utstring(s2_in, s2_data)) {
                sscanf(utstring_body(s2_data), head_scanner, single2_1, single2_2, single2_3);
                ss_get_cat_utstring(s2_in, s2_data);
                ss_get_cat_utstring(s2_in, s2_data);
                ss_get_cat_utstring(s2_in, s2_data);
            } else {
                single2_head_int = 0xffffffffffffffff;
            }            

            while ((read1_head_int & read2_head_int & single1_head_int & single2_head_int) != 0xffffffffffffffff) {
                                
                if ((read1_head_int <= read2_head_int) && (read1_head_int < single1_head_int) && (read1_head_int < single2_head_int)) {
                    
                    last_head_int = read1_head_int;
                    
                    if (read1_head_int == read2_head_int) {  // Write a trimmed mate-pair
                        
                        fputs(utstring_body(r1_data), r1_out);
                        fputs(utstring_body(r2_data), r2_out);
                        
                        if (mates->count) {
                            fputs(utstring_body(r2_data), mates_out);
                            fputs(utstring_body(mate_data), mates_out);
                        }

                        read1_hungry = read2_hungry = 1;

                    } else {    // write read1 singlet
                    
                        read1_singlet_cnt++;
                        fputs(utstring_body(r1_data), s1_out);
                        read1_hungry = 1;
                    
                    }
                    
                } else if ((read2_head_int < single1_head_int) && (read2_head_int < single2_head_int)) {
                    
                    last_head_int = read2_head_int;
                    read2_singlet_cnt++;
                    fputs(utstring_body(r2_data), s2_out);
                    read2_hungry = 1;
                    
                } else if (single1_head_int < single2_head_int) {

                    last_head_int = single1_head_int;
                    fputs(utstring_body(s1_data), s1_out);
                    single1_hungry = 1;
                    
                } else { // single2 must be the smallest
                    
                    last_head_int = single2_head_int;
                    fputs(utstring_body(s2_data), s2_out);
                    single2_hungry = 1;
                    
                }
            
                // Deal with hungry input buffers
                
                if (read1_hungry) {
                    if (ss_get_utstring(r1_in, r1_data)) {
                        sscanf(utstring_body(r1_data), head_scanner, read1_1, read1_2, read1_3);
                        
                        if (read1_head_int <= last_head_int) {
                            fprintf(stderr, "Error: Out of order read detected! %s\n", utstring_body(r1_data));
                            exit(EXIT_FAILURE);
                        }
                        if ((mates->count) && (read1_head_int >= read2_head_int)) {
                            utstring_clear(mate_data);
                            utstring_concat(mate_data, r1_data);    // head
                            ss_get_utstring(r1_in, rev_data);       // seq
                            utstring_concat(r1_data, rev_data);
                            ss_trunc_utstring(rev_data, 1);  // remove newline
                            ss_rev_utstring(rev_data);
                            ss_strcat_utstring(rev_data, "\n"); // add newline back
                            utstring_concat(mate_data, rev_data);
                            ss_get_utstring(r1_in, rev_data);       // extra
                            utstring_concat(r1_data, rev_data);
                            utstring_concat(mate_data, rev_data);
                            ss_get_utstring(r1_in, rev_data);       // qual
                            utstring_concat(r1_data, rev_data);
                            ss_trunc_utstring(rev_data, 1);  // remove newline
                            ss_rev_utstring(rev_data);
                            ss_strcat_utstring(rev_data, "\n"); // add newline back
                            utstring_concat(mate_data, rev_data);
                        } else {
                            ss_get_cat_utstring(r1_in, r1_data);
                            ss_get_cat_utstring(r1_in, r1_data);
                            ss_get_cat_utstring(r1_in, r1_data);
                        }
                    } else {
                        read1_head_int = 0xffffffffffffffff;
                    }                
                    read1_hungry = 0;
                }

                if (read2_hungry) {
                    if (ss_get_utstring(r2_in, r2_data)) {
                        sscanf(utstring_body(r2_data), head_scanner, read2_1, read2_2, read2_3);
                        
                        if (read2_head_int <= last_head_int) {
                            fprintf(stderr, "Error: Out of order read detected! %s\n", utstring_body(r2_data));
                            exit(EXIT_FAILURE);
                        }
                        ss_get_cat_utstring(r2_in, r2_data);
                        ss_get_cat_utstring(r2_in, r2_data);
                        ss_get_cat_utstring(r2_in, r2_data);
                    } else {
                        read2_head_int = 0xffffffffffffffff;
                    }
                    read2_hungry = 0;
                }

                if (single1_hungry) {
                    if (ss_get_utstring(s1_in, s1_data)) {
                        sscanf(utstring_body(s1_data), head_scanner, single1_1, single1_2, single1_3);
                        
                        if (single1_head_int <= last_head_int) {
                            fprintf(stderr, "Error: Out of order read detected! %s\n", utstring_body(s1_data));
                            exit(EXIT_FAILURE);
                        }
                        
                        ss_get_cat_utstring(s1_in, s1_data);
                        ss_get_cat_utstring(s1_in, s1_data);
                        ss_get_cat_utstring(s1_in, s1_data);
                    } else {
                        single1_head_int = 0xffffffffffffffff;
                    }
                    single1_hungry = 0;
                }
                
                if (single2_hungry) {
                    if (ss_get_utstring(s2_in, s2_data)) {
                        sscanf(utstring_body(s2_data), head_scanner, single2_1, single2_2, single2_3);
                        
                        if (single2_head_int <= last_head_int) {
                            fprintf(stderr, "Error: Out of order read detected! %s\n", utstring_body(s2_data));
                            exit(EXIT_FAILURE);
                        }
                        
                        ss_get_cat_utstring(s2_in, s2_data);
                        ss_get_cat_utstring(s2_in, s2_data);
                        ss_get_cat_utstring(s2_in, s2_data);
                    } else {
                        single2_head_int = 0xffffffffffffffff;
                    }
                    single2_hungry = 0;
                }
            }
                        
            fclose(r1_in); 
            fclose(r2_in);                 
            
            fclose(s1_in);
            fclose(s2_in);             
            
            fclose(mates_out); 
            
            fclose(r1_out); 
            fclose(r2_out);                 
            
            fclose(s1_out);
            
            if (singles_flag) {
                fclose(s2_out);
            }
           
            // Free buffers
            utstring_free(r1_data);
            utstring_free(r2_data);
            utstring_free(s1_data);
            utstring_free(s2_data);
            utstring_free(mate_data);
            utstring_free(rev_data);

        }        
    } 

    if (R1_org != R2_org) {
        
        fprintf(stderr, "\nERROR! read1 and read2 fastq files did not contain an equal number of reads. %lu %lu\n", R1_org, R2_org);
        exit(EXIT_FAILURE);
        
    } 

    if (!(R1_org+singlet1_org+singlet2_org)) {
    
        fprintf(stderr, "ERROR! No reads found in input files, or input(s) not found.\n");
        exit(EXIT_FAILURE);
    
    }
    
    if ((R1_org + R2_org) && !(singlet1_cnt + singlet2_cnt)) {
        fprintf(stderr, "\nWarning! read1/read2 files were processed, but no corresponding input singlets were found.\n");
    } 
    
    if (entropy->count) {
        fprintf(stderr, "\nLow complexity reads discarded: Read1: %lu, Read2: %lu, Singlets: %lu %lu\n", comp_r1, comp_r2, comp_s1, comp_s2);
    }

    mp_org = R1_org;
    mp_cnt = R1_cnt;
    
    fprintf(stderr, "\nMatepairs: Before: %lu, After: %lu\n", mp_org, mp_cnt);
    fprintf(stderr, "Singlets: Before: %lu %lu After: %lu %lu\n", singlet1_org, singlet2_org, s1_cnt, s2_cnt);
    fprintf(stderr, "Read1 singlets: %lu, Read2 singlets: %lu, Original singlets: %lu %lu\n", read1_singlet_cnt, read2_singlet_cnt, singlet1_cnt, singlet2_cnt);
    fprintf(stderr, "Total Reads Processed: %lu, Reads retained: %lu\n", 2*mp_org+singlet1_org+singlet2_org, 2*mp_cnt+s1_cnt+s2_cnt);
    
    utstring_free(in_read1_fq_fn);
    utstring_free(in_read2_fq_fn);
    utstring_free(in_single1_fq_fn);
    utstring_free(in_single2_fq_fn);
    utstring_free(out_read1_fq_fn);
    utstring_free(out_read2_fq_fn);
    utstring_free(out_single1_fq_fn);
    utstring_free(out_single2_fq_fn);
    utstring_free(out_mates_fq_fn);
    
    utstring_free(out_read_prefix);
    
    exit(EXIT_SUCCESS);
}

////////////////////////////////////////////////////////////////////////////
// Function trimmer
// Parm: qualdata -- pointer to buffer of fastq encoded Phred quality scores
// Returns: int -- Position to trim to
// Uses global parm values: method_flag, err_prob 

int trimmer(char *qualdata) {
    
    char *q = qualdata;
    double prob = 1.0;
    double psave = 0.0, p2 = 0.0;
    int diff = 0;
    
    if (method_flag) {  // nt space => joint cs space error probability
        psave = pow(10.0, (((double)-(*q - 33))/10.0));
        q++;
        do {
            p2 = psave;
            psave = pow(10.0, (((double)-(*q - 33))/10.0));
            prob *= (1.0 - (p2 * psave));
            q++;
        } while ((prob >= err_prob) && *q);
        
    } else {  // cs space cumulative error prob
        do {
            prob *= (1.0 - pow(10.0, (((double)-(*q - 33))/10.0)));
            q++;
        } while ((prob >= err_prob) && *q);
    }
    
    diff = q - qualdata - 1;

    // If fixed_length reads requested, then truncate to fix_len 
    if (fix_len && (diff >= fix_len)) {
        diff = fix_len;
    }
    
    return diff;
}

////////////////////////////////////////////////////////////////////////////
// Function entropy_calc
// Parm: seqdata -- pointer to buffer of fastq encoded sequence data
//       len -- length of read to consider, NOTE! overridden by global strict_ent
// Returns: int -- (as boolean) keep or reject read
// Uses global parm value: entropy_cutoff 

double entropy_calc(char *seqdata, int len) {
    
    // LUT for color to index
    static const int LUT[128] = {    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, \
        -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, \
        -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, \
        -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, \
        -1, 0,-1, 1,-1,-1,-1, 2,-1,-1,-1,-1,-1,-1, 4,-1, \
        -1,-1,-1,-1, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, \
        -1, 0,-1, 1,-1,-1,-1, 2,-1,-1,-1,-1,-1,-1, 4,-1, \
        -1,-1,-1,-1, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1  };
    
    char *s = seqdata;
    double bits = 0.0;
    double l2 = log(2.0);
    
    len = (strict_ent) ? strlen(seqdata) : len;
    
    int l = len;
    
    // One count for each color pair [ACGTN][ACGTN]
    int counts[5][5] = {{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}};
    
    // Loop through the colorspace values
    for (int x = 0; x < (len-1); x++) {
        counts[LUT[(int)(*s)]][LUT[(int)(*(s+1))]]++;
        s++;
    }
    
    // Discount the dinucleotides with N's from length
    for (int x = 0; x < 5; x++) {
        l -= counts[x][4];
        l -= counts[4][x];
    }
    
    l += counts[4][4];    // Don't double subtract NN.
    
    // calculate information
    for (int x = 0; x < 4; x++) {
        for (int y = 0; y < 4; y++) {
            bits += (counts[x][y]) ? (double) counts[x][y] * (log((double) l / (double) counts[x][y]) / l2) : 0.0;
        }
    }
    
    // Calc bits/nt (based on di-nucleotide info)
    bits = bits / (double) len;    // Use full length, to penalize Ns
    
    return (bits);
}


/////////////////////////////////////////////////////////////////
// 
// fq_stream_trimmer
//
// This function takes an open FASTQ file (gzipped or not) and 
// trims the resulting reads. 
// This is useful for easily implementing a "read trimmer thread".
//

unsigned long int fq_stream_trimmer(UT_string *fq_fn, int pipe_fd, UT_string *out_prefix, int no_pre, int len_pre, unsigned long int *comp_cnt, unsigned long int *org, char split) {
    
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
    
    char *start = NULL; 
    char *end = NULL;
    char *suffix = NULL;
    
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
    
    int x = 0;
    
    while (ss_gzget_utstring(fq_file, head_data)) {
        
        ss_gzget_utstring(fq_file, seq_data);
        ss_gzget_utstring(fq_file, extra_data);
        ss_gzget_utstring(fq_file, qual_data);
        
        if (!split || ((suffix = strchr(utstring_body(head_data), '/')) && (suffix[1] == split))) {

            (*org)++;
            
            if ((x = trimmer(utstring_body(qual_data))) >= min_len) {  // Keep at least some of read
                
                // Reject read if complexity is too low
                if ((entropy_cutoff < 0.0) || (entropy_calc(utstring_body(seq_data), x) >= entropy_cutoff)) {                    
                    
                    // Truncate sequence
                    
                    ss_trim_utstring(seq_data, x);
                    ss_strcat_utstring(seq_data, "\n");
                    ss_trim_utstring(qual_data, x);
                    ss_strcat_utstring(qual_data, "\n");
                    
                    // Fixup the read name
                    utstring_clear(new_head_data);
                    
                    end = strchr(utstring_body(head_data), ':');
                    
                    if (no_pre) {
                        if ((start = strchr(utstring_body(head_data), '|'))) {
                            start++;
                        } else {
                            start = utstring_body(head_data) + 4;
                        }
                        *end = '\0';
                    } else {
                        start = utstring_body(out_prefix);
                    }
                    
                    end++;
                    
                    if (len_pre) {
                        utstring_printf(new_head_data, "@%.2s+%u|%s:%s",utstring_body(head_data)+1,x,start,end);
                    } else {
                        utstring_printf(new_head_data, "@%.2s+%s:%s",utstring_body(head_data)+1,start,end);
                    }
                    
                    fputs(utstring_body(new_head_data), pipe_in);
                    fputs(utstring_body(seq_data), pipe_in);
                    fputs(utstring_body(extra_data), pipe_in);
                    fputs(utstring_body(qual_data), pipe_in);
                    cnt++;
                    
                } else {
                    (*comp_cnt)++;
                }
            }
        }
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
