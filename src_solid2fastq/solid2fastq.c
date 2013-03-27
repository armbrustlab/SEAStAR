/*
 --------------------------------------------------------------------------- #
 Center for Environmental Genomics
 Copyright (C) 2009-2012 University of Washington.
 
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
 
 Authors:
 Vaughn Iverson
 vsi@uw.edu
 
 Chris Berthiaume
 chrisbee@uw.edu
 --------------------------------------------------------------------------- #
 
 ============================================================================
 Name        : solid2fastq.c
 Description : Convert ABI SOLiD output files to FASTQ format
 ============================================================================
*/

#include <seastar_shared.h>
#include <argtable2.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#define READ_SUFFIX_R3 "/1"
#define READ_SUFFIX_F3 "/2"

//////////////////////////////////////////////
// Function prototypes, defined below main()
//////////////////////////////////////////////

int convert_one_read_cs(UT_string *prefix, UT_string *suffix, UT_string *head, UT_string *rdata, UT_string *qdata, UT_string *output);

int convert_one_read_nt(UT_string *prefix, UT_string *suffix, UT_string *head, UT_string *seq, UT_string *sep, UT_string *qual, UT_string *output);

long unsigned int convert_one_stream(UT_string *prefix, UT_string *suffix, int no_suffix, const char *infile_prefix, int bc_flag, const char *bc_id, int out_pipe);

long unsigned int open_input_files(UT_string *prefix, UT_string *suffix, const char *infile_prefix, int bc_flag, const char *bc_id, FILE **csfasta_file, FILE **qual_file, FILE **fastq_file);

/////////////////////////////////////////////////////////////////
// Entry point
/////////////////////////////////////////////////////////////////

int main(int argc, char *argv[]) {
    
    ///////////////////////////
    // Variable declarations
    ///////////////////////////
    
    // Output filenames
    UT_string *read1_fq_fn, *read2_fq_fn, *single1_fq_fn, *single2_fq_fn;
    utstring_new(read1_fq_fn);
    utstring_new(read2_fq_fn);
    utstring_new(single1_fq_fn);
    utstring_new(single2_fq_fn);
    
    // Read name prefix/suffix
    UT_string *read_prefix;
    utstring_new(read_prefix);
    UT_string *read_suffix_R3, *read_suffix_F3;
    utstring_new(read_suffix_R3);
    utstring_new(read_suffix_F3);
    
    // Pipes used for inter-thread communication
    int R_pipe[2];
    int F_pipe[2];
    int r1_pipe[2];
    int r2_pipe[2];
    int s1_pipe[2];
    int s2_pipe[2];
    
    // Flags
    int R3_hungry = 0, F3_hungry = 0;
    
    // Read counters
    unsigned long int mp_cnt = 0, R3_singlet_cnt = 0, F3_singlet_cnt = 0;
    unsigned long int R3_cnt = 0, F3_cnt = 0;
    
    // All done with variable declarations!!
    
    ///////////////////////////////////
    // Command line argtable settings
    ///////////////////////////////////
    
    struct arg_str *pre_read_id = arg_str0(NULL, "prefix", "<string>", "Prefix to add to read identifiers. [out_prefix]");
    struct arg_str *bc = arg_str0("b", "bc", "<BC>", "Barcode string in filename. Only needed for barcoded runs. [NULL]");
    struct arg_lit *no_pre = arg_lit0("n", "no_prefix", "Do not change the read names in any way. [NULL]");
    struct arg_lit *no_suf = arg_lit0("x", "no_suffix", "Suppress /1 or /2 suffix additions to read IDs. [NULL]");
    struct arg_lit *gzip = arg_lit0("z", "gzip", "Output converted files in gzip compressed format. [NULL]");
    struct arg_lit *singles = arg_lit0("s", "singles", "Write two singlet files, one for each mate-paired input file. [NULL]");
    struct arg_file *input = arg_file1(NULL, NULL, "<in_prefix>", "Input file prefix: (e.g. <in_prefix>_F3.csfasta <in_prefix>_F3_QV.qual [<in_prefix>_R3.csfasta <in_prefix>_R3_QV.qual] or <in_prefix>_<bc>_F3.fastq <in_prefix>_<bc>_R3.fastq)");
    struct arg_file *output = arg_file1(NULL, NULL, "<out_prefix>", "Output file prefix: (e.g. <out_prefix>.single.fastq [<out_prefix>.read1.fastq <out_prefix>.read2.fastq]) ");
	struct arg_lit *version = arg_lit0(NULL,"version","Print the build version and exit.");  
    struct arg_lit *h = arg_lit0("h", "help", "Request help.");
    struct arg_end *end = arg_end(20);
    
    void *argtable[] = {h,version,pre_read_id,no_pre,no_suf,gzip,bc,singles,input,output,end};
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
        fprintf(stderr, "\nInput prefix is the part of the filename shown in the\n");
        fprintf(stderr, "string on the \"Title\" line of the CSFASTA format inputs.\n");
        fprintf(stderr, "The output prefix will be used to name the output files and,\n");
        fprintf(stderr, "by default, to prefix the output read names.\n");
        fprintf(stderr, "Mate-paired read files are automatically used if present.\n");
        fprintf(stderr, "Three fastq output files only produced for mate-paired inputs.\n");
        fprintf(stderr, "\nNote: Input and output files may be gzipped.\n");
        exit(EXIT_FAILURE);
    }    
    
    if (arg_errors) { 
        arg_print_errors(stderr, end, "solid2fastq");
        arg_print_syntaxv(stderr, argtable, "\n");
        exit(EXIT_FAILURE);
    }
        
    if (pre_read_id->count) {
        
        if (no_pre->count) {
            fprintf(stderr, "Error: Both --prefix and --no_prefix were specified.\n");
            exit(EXIT_FAILURE);
        }
        
        if (!strlen(pre_read_id->sval[0])) {
            fprintf(stderr, "Read ID prefix may not be zero length.\n");
            exit(EXIT_FAILURE);
        } 
        
        if (strchr(pre_read_id->sval[0], ':') || strchr(pre_read_id->sval[0], '|') || strchr(pre_read_id->sval[0], '+') || strchr(pre_read_id->sval[0], '/')) {
            fprintf(stderr, "Read ID prefix '%s' may not contain the characters ':', '|', '+' or '/'.\n", pre_read_id->sval[0]);
            exit(EXIT_FAILURE);
        }
        
        // Build default read ID prefix
        utstring_printf(read_prefix, "%s:", pre_read_id->sval[0]);
        
    } else {
        
        if (no_pre->count) {
            utstring_printf(read_prefix, ":", output->filename[0]);
        } else {
            
            if (strchr(output->filename[0], ':') || strchr(output->filename[0], '|') || strchr(output->filename[0], '+') || strchr(output->filename[0], '/')) {
                fprintf(stderr, "Read ID prefix '%s' (from output prefix) may not contain the characters ':', '|', '+' or '/'.\n", output->filename[0]);
                fprintf(stderr, "Hint: Use the --prefix parameter if the output file prefix contains path information.\n");
                exit(EXIT_FAILURE);
            }        
            
            // Build default read ID prefix
            utstring_printf(read_prefix, "%s:", output->filename[0]);
        }
    }

    // Set up suffixes
    utstring_printf(read_suffix_R3, READ_SUFFIX_R3);
    utstring_printf(read_suffix_F3, READ_SUFFIX_F3);
    
    // Check for null string filename prefixes
    if (!(strlen(input->filename[0]) && strlen(output->filename[0]))) {
        fprintf(stderr, "Error: NULL prefix strings are not permitted.\n");
        exit(EXIT_FAILURE);
    }
    
    // Construct filenames
    // Output filenames
    
    utstring_printf(read1_fq_fn, "%s.read1.fastq", output->filename[0]);
    utstring_printf(read2_fq_fn, "%s.read2.fastq", output->filename[0]);
    
    if (singles->count) {
        utstring_printf(single1_fq_fn, "%s.single1.fastq", output->filename[0]);
        utstring_printf(single2_fq_fn, "%s.single2.fastq", output->filename[0]);
    } else {
        utstring_printf(single1_fq_fn, "%s.single.fastq", output->filename[0]);
    }
    
    ///////////////////////////////////
    // Begin processing!
    ///////////////////////////////////
    
#ifdef _OPENMP    
    omp_set_num_threads(7);
#endif    
    
    pipe(R_pipe);
    pipe(F_pipe);
    pipe(r1_pipe);
    pipe(r2_pipe);
    pipe(s1_pipe);
    pipe(s2_pipe);
    
#pragma omp parallel sections default(shared)       
    {
        
#pragma omp section 
        {   // R3_Reader / Converter
            R3_cnt = convert_one_stream(read_prefix, read_suffix_R3, no_suf->count, input->filename[0], bc->count, bc->sval[0], R_pipe[1]);
        }
        
#pragma omp section
        {   // F3_reader / Converter
            F3_cnt = convert_one_stream(read_prefix, read_suffix_F3, no_suf->count, input->filename[0], bc->count, bc->sval[0], F_pipe[1]);
        }
        
#pragma omp section 
        {   // read2_writer
            ss_stream_writer(read2_fq_fn, r2_pipe[0], gzip->count);
        }
        
#pragma omp section 
        {   // read1_writer
            ss_stream_writer(read1_fq_fn, r1_pipe[0], gzip->count);
        }

#pragma omp section 
        {   // single2_writer
            ss_stream_writer(single2_fq_fn, s2_pipe[0], gzip->count);
        }
        
#pragma omp section 
        {   // single1_writer
            ss_stream_writer(single1_fq_fn, s1_pipe[0], gzip->count);
        }
        
#pragma omp section
        {   // Demuxer 
            
            // Read header scanf pattern string:
            // Parses out the three bead coordinates from the first FASTQ "header" line for a read
            // e.g.  @lambda:TA+1231_1682_1381/1
            //                  xxxx yyyy zzzz    where xxxx, yyyy and zzzz are the bead coordinates
            
            const char *head_scanner = "@%*[^:]:%hu_%hu_%hu"; 
            
            // Encoded header 64-bit values 
            // Each solid read coordiante is xxxx_yyyy_zzzz, the code below packs these three
            // (< 16-bit) unsigned values into a 64-bit unsigned integer as:
            // 0000xxxxyyyyzzzz in hex notation

            long unsigned int R3_head_int = 0, F3_head_int = 0;
            
            // Pointers to short int parts of the above 64-bit ints (used by sscanf below)
            short unsigned int *R3_1 = (short unsigned int *)(&R3_head_int)+2;
            short unsigned int *R3_2 = (short unsigned int *)(&R3_head_int)+1;
            short unsigned int *R3_3 = (short unsigned int *)(&R3_head_int);
            short unsigned int *F3_1 = (short unsigned int *)(&F3_head_int)+2;
            short unsigned int *F3_2 = (short unsigned int *)(&F3_head_int)+1;
            short unsigned int *F3_3 = (short unsigned int *)(&F3_head_int);
            
            UT_string *R_data, *F_data;
            utstring_new(R_data);
            utstring_new(F_data);
            
            FILE *R_in = fdopen(R_pipe[0],"r");                 
            FILE *F_in = fdopen(F_pipe[0],"r");                 

            FILE *r1_out = fdopen(r1_pipe[1],"w"); 
            FILE *r2_out = fdopen(r2_pipe[1],"w");                 

            FILE *s1_out = fdopen(s1_pipe[1],"w");
            FILE *s2_out = fdopen(s2_pipe[1],"w");             

            if (!singles->count) {
                fclose(s2_out);
                s2_out = s1_out;
            }
            
            // Prime flags for loop below
            R3_hungry = 0;
            F3_hungry = 0;
            
            // Prime buffers
            if (ss_get_utstring(R_in, R_data)) {
                sscanf(utstring_body(R_data), head_scanner, R3_1, R3_2, R3_3);
                ss_get_cat_utstring(R_in, R_data);
                ss_get_cat_utstring(R_in, R_data);
                ss_get_cat_utstring(R_in, R_data);
            } else {
                R3_head_int = 0xffffffffffffffff;
            }
            
            if (ss_get_utstring(F_in, F_data)) {
                sscanf(utstring_body(F_data), head_scanner, F3_1, F3_2, F3_3);
                ss_get_cat_utstring(F_in, F_data);
                ss_get_cat_utstring(F_in, F_data);
                ss_get_cat_utstring(F_in, F_data);
            } else {
                F3_head_int = 0xffffffffffffffff;
            }
            
            while ((F3_head_int != 0xffffffffffffffff) || (R3_head_int != 0xffffffffffffffff)) {   // Figure out which read(s) to write where
                
                if (R3_head_int == F3_head_int) {  // Matched reads, write both
                    fputs(utstring_body(F_data), r2_out);
                    fputs(utstring_body(R_data), r1_out);
                    R3_hungry = 1;
                    F3_hungry = 1;
                    mp_cnt++;
                } else if (R3_head_int > F3_head_int) { // write F3 singlet
                    fputs(utstring_body(F_data), s2_out);
                    F3_hungry = 1;
                    F3_singlet_cnt++;
                } else {    // write R3 singlet
                    fputs(utstring_body(R_data), s1_out);
                    R3_hungry = 1;
                    R3_singlet_cnt++;
                }  
                
                if (R3_hungry) {
                    
                    if (ss_get_utstring(R_in, R_data)) {
                        sscanf(utstring_body(R_data), head_scanner, R3_1, R3_2, R3_3);
                        ss_get_cat_utstring(R_in, R_data);
                        ss_get_cat_utstring(R_in, R_data);
                        ss_get_cat_utstring(R_in, R_data);
                    } else {
                        R3_head_int = 0xffffffffffffffff;
                    }
                    
                    R3_hungry = 0;
                }
                
                if (F3_hungry) {
                    
                    if (ss_get_utstring(F_in, F_data)) {
                        sscanf(utstring_body(F_data), head_scanner, F3_1, F3_2, F3_3);
                        ss_get_cat_utstring(F_in, F_data);
                        ss_get_cat_utstring(F_in, F_data);
                        ss_get_cat_utstring(F_in, F_data);
                    } else {
                        F3_head_int = 0xffffffffffffffff;
                    }
                    
                    F3_hungry = 0;
                }
             }
            
            // Close pipes
            fclose(R_in);
            fclose(F_in);
            fclose(r1_out);
            fclose(r2_out);
            fclose(s1_out);
            if (singles->count) {
                fclose(s2_out);
            }
            
            // Free UT_strings
            utstring_free(R_data);
            utstring_free(F_data);
        }            
    }

    // Free variable length strings

    utstring_free(read1_fq_fn);
    utstring_free(read2_fq_fn);
    utstring_free(single1_fq_fn);
    utstring_free(single2_fq_fn);
    utstring_free(read_prefix);
    
    if ((R3_cnt) && (!F3_cnt)) { // Only R3 Input?  =  Suspicious
        fprintf(stderr, "Warning! No F3 reads processed. Check for missing or empty F3 csfasta/qual files.\n");
    }
    
    if (F3_singlet_cnt+R3_singlet_cnt+2*mp_cnt) {
        printf("\nMatepairs: %lu, Matepaired Reads: %lu\n", mp_cnt, mp_cnt*2);
        printf("F3 singlets: %lu, R3 singlets: %lu, Total singlets: %lu\n", F3_singlet_cnt, R3_singlet_cnt, F3_singlet_cnt+R3_singlet_cnt);
        printf("Total Reads Processed: %lu\n", F3_singlet_cnt+R3_singlet_cnt+2*mp_cnt);
        exit(EXIT_SUCCESS);
    } else {
        fprintf(stderr, "Error! No reads found in input files, or input(s) not found.\n");
        exit(EXIT_FAILURE);
    }
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

long unsigned int convert_one_stream(UT_string *prefix, UT_string *suffix, int no_suffix_flag, const char *infile_prefix, int bc_flag, const char *bc_id, int out_pipe)
{

    long unsigned int cnt = 0;
    
    UT_string *csfasta_head, *qual_head, *fastq_head;
    utstring_new(csfasta_head);
    utstring_new(qual_head);
    utstring_new(fastq_head);
    
    UT_string *csfasta_data, *qual_data, *fastq_seq, *fastq_sep, *fastq_qual;
    utstring_new(csfasta_data);
    utstring_new(qual_data);
    utstring_new(fastq_seq);
    utstring_new(fastq_sep);
    utstring_new(fastq_qual);
    
    UT_string *actual_suffix;
    utstring_new(actual_suffix);
    
    UT_string *out_buf;
    utstring_new(out_buf);
    
    UT_string *csfasta_fn, *qual_fn;
    utstring_new(csfasta_fn);
    utstring_new(qual_fn);

    // Input filename base
    UT_string *infile_base;
    utstring_new(infile_base);

    FILE *csfasta_file = NULL;
    FILE *qual_file = NULL;
    FILE *fastq_file = NULL;
    FILE *out_file = NULL;
    
    open_input_files(prefix, suffix, infile_prefix, bc_flag, bc_id, &csfasta_file, &qual_file, &fastq_file);
    out_file = fdopen(out_pipe, "w");
    
    if (no_suffix_flag) {
        utstring_printf(actual_suffix, "");
    } else {
        utstring_printf(actual_suffix, utstring_body(suffix));
    }
    
    // At least one of sequence and qual files found
    if (fastq_file) {
        // Read initial data from input files, eating '#' comment lines from tops of files.
        ss_gzget_utstring(fastq_file, fastq_head);
        while (utstring_body(fastq_head)[0] == '#') { ss_gzget_utstring(fastq_file, fastq_head); }
        do {
            
            ss_gzget_utstring(fastq_file, fastq_seq);
            ss_gzget_utstring(fastq_file, fastq_sep);
            ss_gzget_utstring(fastq_file, fastq_qual);
            
            convert_one_read_nt(prefix, actual_suffix, fastq_head, fastq_seq, fastq_sep, fastq_qual, out_buf);
            
            fputs(utstring_body(out_buf), out_file);
            
            ss_gzget_utstring(fastq_file, fastq_head);
            
            cnt++;
            
        } while (utstring_len(fastq_head));
    } else if (csfasta_file && qual_file) {
        
        // Read initial data from input files, eating '#' comment lines from tops of files.
        ss_gzget_utstring(csfasta_file, csfasta_head);
        while (utstring_body(csfasta_head)[0] == '#') { ss_gzget_utstring(csfasta_file, csfasta_head); }
        ss_gzget_utstring(qual_file, qual_head);
        while (utstring_body(qual_head)[0] == '#') { ss_gzget_utstring(qual_file, qual_head); }
    
        do {

            ss_gzget_utstring(csfasta_file, csfasta_data);
            ss_gzget_utstring(qual_file, qual_data);

            convert_one_read_cs(prefix, actual_suffix, csfasta_head, csfasta_data, qual_data, out_buf);
        
            fputs(utstring_body(out_buf), out_file);
        
            ss_gzget_utstring(csfasta_file, csfasta_head);
            ss_gzget_utstring(qual_file, qual_head);
        
            cnt++;
        
        } while (utstring_len(csfasta_head));
    }
    
    fclose(out_file);
    if (csfasta_file) gzclose(csfasta_file);
    if (qual_file) gzclose(qual_file);
    if (fastq_file) gzclose(fastq_file);
    utstring_free(csfasta_head);
    utstring_free(qual_head);
    utstring_free(fastq_head);
    utstring_free(csfasta_data);
    utstring_free(qual_data);
    utstring_free(fastq_seq);
    utstring_free(fastq_sep);
    utstring_free(fastq_qual);
    utstring_free(out_buf);
    utstring_free(infile_base);
    utstring_free(actual_suffix);
    
    return (cnt);
}

//////////////////////////////////////////////////////
// Convert a csfasta/qual sequence record to SEAStAR comptabile fastq format.
// Colorspace numbers are converted to pseudo-basespace equivalents and quality
// is encoded as ASCII characters following the Sanger fastq convention.  The
// input read should have a header formatted as @nnn_nnn_nnn_SUFFIX.  The
// header is converted to @XX+PREFIX:nnn_nnn_nnn/{1,2} where /1 or /2 is passed
// in as suffix and XX are the first two characters of the sequence.  The
// first two characters of the sequence and quality are removed from output.
int convert_one_read_cs(UT_string *prefix, UT_string *suffix, UT_string *head, UT_string *rdata, UT_string *qdata, UT_string *output) {

    // LUT for transliteration
    static const char col_convert[] = 
    "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNACGTNNNNNNNNNNNNNANCNNNGNNNNNNNNNNNNTNNNNNNNNNNNNANCNNNGNNNNNNNNNNNNTNNNNNNNNNNN";    
    
    static const char qual_convert[] = 
    "!!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";    
    
    char *ptr = NULL;
    char *tokptr = NULL;
    char *save_ptr = NULL;
    
    // Convert the data string
    
    ptr = utstring_body(rdata)+1;  // +1 skips primer base
    
    while (*ptr != '\n') {
        *ptr = col_convert[(int)*ptr];
        ptr++;
    }
    
    // Detect reads named like: >nnn_nnn_nnn_SUFFIX
    //                         And remove:  ^^^^^^^                 
    
    int uc = 0;
    for (char *c = utstring_body(head); *c; c++) {
        if (*c == '_' && ++uc == 3) {
            ss_trim_utstring(head,c-utstring_body(head));  // Trunc at the 3rd '_'
            break;
        }
    }

    utstring_clear(output);  // Reset the output buffer
    
    utstring_printf(output, "@%.2s+%s%s%s\n%s+\n",
                    utstring_body(rdata),
                    utstring_body(prefix),
                    utstring_body(head)+1, // +1 skips '>'
                    utstring_body(suffix),
                    utstring_body(rdata)+2 // +2 skips base and primer color
                    );
    
    // Note, we reuse the rdata buffer now to build the FASTQ quality chars
    ptr = utstring_body(rdata)+2;  // Get the buffer ready for qualities
    
    strtok_r(utstring_body(qdata), " ", &save_ptr);    // Tokenize on spaces
    tokptr = strtok_r(NULL, " ", &save_ptr);    // Eat the first quality value
    
    do {
        *ptr++ = qual_convert[atoi(tokptr) + 1];
    } while((tokptr = strtok_r(NULL, " \n", &save_ptr)));
    
    ss_strcat_utstring(output, (utstring_body(rdata)+2));

    return(1);
}

//////////////////////////////////////////////////////
// Convert a fastq nucleotide record to SEAStAR comptabile fastq format.  The
// input read should have a header formatted as @nnn_nnn_nnn_SUFFIX.  The
// header is converted to @PREFIX:nnn_nnn_nnn/{1,2} where /1 or /2 is passed in
// as suffix.
int convert_one_read_nt(UT_string *prefix, UT_string *suffix, UT_string *head, UT_string *seq, UT_string *sep, UT_string *qual, UT_string *output) {
    // Detect reads named like: @nnn_nnn_nnn_SUFFIX
    //                         And remove:  ^^^^^^^
    
    int uc = 0;
    for (char *c = utstring_body(head); *c; c++) {
        if (*c == '_' && ++uc == 3) {
            ss_trim_utstring(head,c-utstring_body(head));  // Trunc at the 3rd '_'
            break;
        }
    }
    
    utstring_clear(output);  // Reset the output buffer
    
    utstring_printf(output, "@%s%s%s\n%s%s%s",
                    utstring_body(prefix),
                    utstring_body(head)+1, // +1 skips '@'
                    utstring_body(suffix),
                    utstring_body(seq),
                    utstring_body(sep),
                    utstring_body(qual)
                    );
    
    return(1);
}

//////////////////////////////////////////////////////
// Open files matching prefix, read type (F3 or R3/F5), and barcode.  Matching
// fastq file names are attempted first, then csfasta/qual files.
//
// Will identify the following file name patterns
// /1 variants
// <in_prefix>_<BC>_F3.fastq
// <in_prefix>_F3[.csfasta|_QV.qual|.QV.qual]
// <in_prefix>_<BC>_F3[.csfasta|.QV.qual]
// <in_prefix>_F3_[<BC>.csfasta|_QV_<BC>.qual]
//
// /2 variants
// <in_prefix>_<BC>_R3.fastq
// <in_prefix>_[R3|F5-BC|F5-P2|F5-DNA|F5-RNA][.csfasta|_QV.qual|.QV.qual]
// <in_prefix>_<BC>_[R3|F5-BC|F5-P2|F5-DNA|F5-RNA][.csfasta|.QV.qual]
// <in_prefix>_[R3|F5-BC|F5-P2|F5-DNA|F5-RNA]_[<BC>.csfasta|_QV_<BC>.qual]
long unsigned int open_input_files(UT_string *prefix, UT_string *suffix, const char *infile_prefix, int bc_flag, const char *bc_id, FILE **csfasta_file, FILE **qual_file, FILE **fastq_file)
{
    UT_string *csfasta_fn, *qual_fn, *fastq_fn;
    utstring_new(csfasta_fn);
    utstring_new(qual_fn);
    utstring_new(fastq_fn);
    
    // Input filename base
    UT_string *infile_base;
    utstring_new(infile_base);
    
    // Construct F3 and R3 file names, try to open files
    // F3
    if (!strcmp(utstring_body(suffix), READ_SUFFIX_F3)) {
        // Check for input nucleotide fastq from SOLiD first.
        // It's OK to assume a bc_id string is always applicable here because
        // SOLiD nucleotide FASTQ files extracted from XSQ files will always
        // have a barcode, even if that barcode is only "default" for
        // non-barcoded runs.
        utstring_printf(fastq_fn, "%s_%s_F3.fastq", infile_prefix, bc_id);
        
        // A barcode wasn't specified and/or there was no matching FASTQ file
        if (!(bc_flag && (*fastq_file = ss_get_gzFile(utstring_body(fastq_fn), "r")))) {
            
            // Different versions of SOLiD instruments differ in both barcode
            // string placement within the native file name and in use of _QV or
            // .QV to designate quality files.  Handle all variations in the
            // following loops.
            for (int dot_QV = 0; dot_QV < 2; dot_QV++) {
                for (int bc_before_suffix = 0; bc_before_suffix < 2; bc_before_suffix++) {
                    utstring_clear(infile_base);
                    utstring_clear(csfasta_fn);
                    utstring_clear(qual_fn);
                    
                    if (bc_flag && bc_before_suffix) {
                        utstring_printf(infile_base, "%s", infile_prefix);
                    } else {
                        utstring_printf(infile_base, "%s_F3", infile_prefix);
                    }
                    if (bc_flag) {
                        if (bc_before_suffix) {
                            utstring_printf(csfasta_fn, "%s_%s_F3.csfasta", utstring_body(infile_base), bc_id);
                            utstring_printf(qual_fn, "%s_%s_F3.QV.qual", utstring_body(infile_base), bc_id);
                        } else {
                            utstring_printf(csfasta_fn, "%s_%s.csfasta", utstring_body(infile_base), bc_id);
                            utstring_printf(qual_fn, "%s_QV_%s.qual", utstring_body(infile_base), bc_id);
                        }
                    } else {
                        utstring_printf(csfasta_fn, "%s.csfasta", utstring_body(infile_base));
                        if (dot_QV) {
                            utstring_printf(qual_fn, "%s.QV.qual", utstring_body(infile_base));
                        } else {
                            utstring_printf(qual_fn, "%s_QV.qual", utstring_body(infile_base));
                        }
                    }
                    
                    // Try to open the csfasta file
                    if (!(*csfasta_file = ss_get_gzFile(utstring_body(csfasta_fn), "r"))) {
                        continue;  // Found correct read1 csfasta file name
                    }
                    
                    // Try to open the qual file
                    if (*qual_file = ss_get_gzFile(utstring_body(qual_fn), "r")) {
                        break;  // found correct read1 qual file name
                    }
                }
                if (*qual_file) break;  // found correct read1 file names
            }
        }
    } else {
        // R3
#define TRIES 5
        // Check for input nucleotide fastq from SOLiD first.
        // It's OK to assume a bc_id string is always applicable here because
        // SOLiD nucleotide FASTQ files extracted from XSQ files will always
        // have a barcode, even if that barcode is only "default" for
        // non-barcoded runs.
        utstring_printf(fastq_fn, "%s_%s_R3.fastq", infile_prefix, bc_id);  // assuming this is always R3, but should add F5-etc... support later
        
        // A barcode wasn't specified and/or there was no matching FASTQ file
        if (!(bc_flag && (*fastq_file = ss_get_gzFile(utstring_body(fastq_fn), "r")))) {
            
            // Try the various possible versions of /2 file suffix
            char bc_suf[TRIES][7] = { "R3", "F5-BC", "F5-P2", "F5-DNA", "F5-RNA" };
            for (int try = 0; try < TRIES; try++) {
                
                // Different versions of SOLiD instruments differ in both barcode
                // string placement within the native file name and in use of _QV or
                // .QV to designate quality files.  Handle all variations in the
                // following loops.
                for (int dot_QV = 0; dot_QV < 2; dot_QV++) {
                    for (int bc_before_suffix = 0; bc_before_suffix < 2; bc_before_suffix++) {
                        utstring_clear(infile_base);
                        utstring_clear(csfasta_fn);
                        utstring_clear(qual_fn);
                        if (bc_flag && bc_before_suffix) {
                            utstring_printf(infile_base, "%s", infile_prefix);
                        } else {
                            utstring_printf(infile_base, "%s_%s", infile_prefix, bc_suf[try]);
                        }
                        
                        if (bc_flag) {
                            if (bc_before_suffix) {
                                utstring_printf(csfasta_fn, "%s_%s_%s.csfasta", utstring_body(infile_base), bc_id, bc_suf[try]);
                                utstring_printf(qual_fn, "%s_%s_%s.QV.qual", utstring_body(infile_base), bc_id, bc_suf[try]);
                            } else {
                                utstring_printf(csfasta_fn, "%s_%s.csfasta", utstring_body(infile_base), bc_id);
                                utstring_printf(qual_fn, "%s_QV_%s.qual", utstring_body(infile_base), bc_id);
                            }
                        } else {
                            utstring_printf(csfasta_fn, "%s.csfasta", utstring_body(infile_base));
                            if (dot_QV) {
                                utstring_printf(qual_fn, "%s.QV.qual", utstring_body(infile_base));
                            } else {
                                utstring_printf(qual_fn, "%s_QV.qual", utstring_body(infile_base));
                            }
                        }
                        
                        // Try to open the csfasta file
                        if (!(*csfasta_file = ss_get_gzFile(utstring_body(csfasta_fn), "r"))) {
                            continue;  // Found correct read2 csfasta file name
                        }
                        
                        // Try to open the qual file
                        if (*qual_file = ss_get_gzFile(utstring_body(qual_fn), "r")) {
                            break; // found correct files
                        }
                    }
                    if (*qual_file) break;  // found correct files
                }
                if (*qual_file) break; // found correct files
            }
        }
    }
    
    utstring_free(infile_base);
    utstring_free(csfasta_fn);
    utstring_free(qual_fn);
    utstring_free(fastq_fn);
}
