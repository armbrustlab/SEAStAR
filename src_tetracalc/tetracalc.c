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
 
 ============================================================================
 Name        : tetracalc
 Description : Reads fasta files and calculates the "Tetra" sequence 
               correlation between these sequences.
 For more info see: 
 ============================================================================
*/

#include <seastar_shared.h>
#include "tune_tabs.h"
#include <assert.h>
#include <argtable2.h>
#include <utlist.h>

// This maintains a version of the program that builds and runs correctly
// when OpenMP support is disabled in--or absent from--the compiler.
#ifdef _OPENMP
#include <omp.h>
#endif

/////////////////////////
// Macros 
/////////////////////////


/////////////////////////
// Typedefs 
/////////////////////////

// This is the "per sequence" structure, corresponding to each sequence in a multi-FASTA
typedef struct seq_info_st {
    UT_string *seq_name;
    UT_string *seq;
    unsigned int start;
    unsigned int stop;
    struct seq_info_st *next;    // Note, this must be called "next", for utlist.h
} seq_info;

// This is the "per scaffold" structure, with multi-FASTA semantics determined by the '-f' commandline parm
typedef struct seq_stat_st {
    seq_info *seqs;
    UT_string *scaffold_name;
    unsigned int tot_len;
    unsigned int num_seqs;
    unsigned int tetra[4][4][4][4];
    unsigned int tri_f[4][4][4];
    unsigned int tri_b[4][4][4];
    unsigned int di_c[4][4];
    unsigned int nuc[4];
    long unsigned int index;
    double Z[256];
    double Z3[64];
    struct pair_stat_st *best;
    struct seq_info_st *next;    // Note, this must be called "next", for utlist.h
} seq_entry;

// This struct holds the results of sequence correlation analysis for each pair of sequences 
typedef struct pair_stat_st {
    seq_entry *x_seq;
    seq_entry *y_seq;
    double cor;
    //double slope;
    //double intercept;
    double cor3;
    //double slope3;
    //double intercept3;
    double cor_thresh;
    double sec_thresh;
    double tri_thresh;
    // double tetra_thresh;
} stat_entry;

/////////////////////////
// Function prototypes
/////////////////////////

// All functions defined below main() in this file

seq_info *read_fasta_seq(gzFile *fn);        // Pointer to open file descriptor

void cleanup(char *seq);    // Cleans up sequence symbols 

seq_entry *new_seq_entry();    // seq_entry create / destroy

void free_seq_entry(seq_entry *seq);

void reg_corr(stat_entry *s_ptr); // Perform regression analysis

void calc_Z_stats(seq_entry *seq); // Calculate Teeling Z-stats

// Calculate whether this correlation represents the best one seen
// for either of the sequences.

int new_best(stat_entry *this_cor, unsigned int second_length);

// int new_best_tuned(stat_entry *this_cor);

int seq_entry_length_cmp(const void *ptr2, const void *ptr1);

double *interp_tune_table(double value, float *row_tab, float *val_tab, int tab_rows, int tab_cols);

double calc_thresh(stat_entry *cor, float *col_tab, double *thresh_tab, int tab_cols);

double pass_correlation(stat_entry *pair, unsigned int second_length);

stat_entry *max_correlation(stat_entry *pair1, stat_entry *pair2, unsigned int second_length);

int limit_thresh(double *thresh_tab, int tab_cols, double limit);

/////////////////////////////////////////////////////////////////
// Entry point
/////////////////////////////////////////////////////////////////

int main (int argc, char *argv[]) {
    
#ifndef _OPENMP
    fprintf(stderr, "\nERROR: Program built with compiler lacking OpenMP support.\n");
    fprintf(stderr, "See SEAStAR README file for information about suitable compilers.\n");
    exit(EXIT_FAILURE);
#endif    
    
    ///////////////////////////
    // Variable declarations
    /////////////////////////// 

    // Input file handles
    gzFile *f_hnd = NULL;    
    
    // Sequence table
    seq_entry *seq_list = NULL, *seq = NULL;
    
    // Sequences
    seq_info *seq_ptr = NULL;
    
    // Array of sequence pointers, make OpenMP possible!
    seq_entry **seq_array = NULL, **seq_pptr = NULL;
    long unsigned int num_seq = 0;
    
    // LUTs for calculating tetra array index
    unsigned int seq_LUT[256];
    for (int x = 0; x < 256; x++) {
        seq_LUT[x] = 100;
    }
    seq_LUT['A'] = 0; seq_LUT['C'] = 1; seq_LUT['G'] = 2; seq_LUT['T'] = 3;
    const int max_tetra_sum = 4*seq_LUT['T'];    

    // Correlation matrix
    stat_entry *cor_entries = NULL;
    stat_entry **cor_matrix = NULL;
    
    // Parameters with default values
    unsigned int length_cutoff = 5000;
    unsigned int second_length = 0;
    double cor_thresh = 0.9;
    double sec_thresh = 0;
    double tri_thresh = 0;
    
    unsigned long int all_seq_len = 0;
    
    double *tetra_merge_tab = NULL;
    double *tri_merge_tab = NULL;
    double *tri_cutoff_tab = NULL;
    double *tetra_cutoff_tab = NULL;

    ///////////////////////////////////
    // Command line argtable settings
    ///////////////////////////////////
    
    struct arg_rem *sep1 = arg_rem(NULL, "\n==== Essential options ====\n");
    struct arg_lit *h = arg_lit0("h","help","Request help.");
    struct arg_lit *version = arg_lit0(NULL,"version","Print the build version and exit.");
    struct arg_dbl *merge_tar = arg_dbl0(NULL,"merge_tar","<n>","Target fraction of within genome merging [0.95]");
    struct arg_int *min_len = arg_int0("m","min_len","<n>","Minimum sequence length for consideration [5000]");
    struct arg_int *num_threads = arg_int0(NULL,"num_threads","<n>","Number of threads to use on multicore systems. [NUM_PROCESSORS]");
    struct arg_rem *sep2 = arg_rem(NULL, "\n==== Advanced / Experimental options ====\n");
    struct arg_lit *fixed = arg_lit0("f","fixed","Use fixed correlation thresholds instead of the merge_tar tuned thresholds [FALSE].");
    struct arg_dbl *cor_parm = arg_dbl0("t","cor_thresh","<n>","Minimum tetra-nuc correlation (R value) for a link to be reported [0.9]");
    struct arg_dbl *sec_parm = arg_dbl0("s","sec_thresh","<n>","Secondary minimum tetra-nuc correlation  [cor_thresh]");
    struct arg_dbl *tri_parm = arg_dbl0("r","tri_thresh","<n>","Minimum tri-nuc correlation  [cor_thresh]");
    struct arg_int *sec_len = arg_int0("c","sec_len","<n>","Maximum sequence length for consideration of tri-nuc [min_len]");    
    struct arg_lit *sf = arg_lit0(NULL,"seq_files","Treat the sequences in each input file as a single scaffold.");
    struct arg_lit *greedy = arg_lit0("g","greedy","Greedy recalculation when merged seqs were the best correlation for other seqs. [FALSE]");
    struct arg_int *chunk = arg_int0(NULL,"chunk","<n>","Break input sequences into chunks of N bases for analysis [FALSE]");    
    struct arg_file *zstats = arg_file0("z","z_stats","<file>","Output the Z-statistics for each sequence sequence, as a table to <file>");
    struct arg_lit *verbose = arg_lit0(NULL,"verbose","Output verbose runtime statistics to stderr [FALSE].");
    struct arg_rem *sep3 = arg_rem(NULL, "\n==== Input files ====\n");    
    struct arg_file *f = arg_filen(NULL,NULL,"<file>",1,10000,"Input file(s), nucleotides in FASTA format. STDIN with '-'. gzip supported");
    struct arg_end *end = arg_end(20);
    void *argtable[] = {sep1,h,version,merge_tar,min_len,num_threads,sep2,fixed,cor_parm,sec_parm,tri_parm,sec_len,sf,greedy,chunk,zstats,verbose,sep3,f,end};
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
        fprintf(stderr, "\ntetracalc is a tool that uses tetra- (and tri-) nucleotide usage statistics to cluster\n");
        fprintf(stderr, "input sequences into taxonomically correlated (genomic) \"bins\".\n");
        fprintf(stderr, "\nUsage: tetracalc [options] <file1> [<file2>]... > <outfile.json>\n");
        fprintf(stderr, "\n[options] : ");
        
        arg_print_syntaxv(stderr, argtable, "\n");
        arg_print_glossary(stderr, argtable, "%-25s %s\n");
		fprintf(stderr, "\n");
        exit(EXIT_FAILURE);
    }
    
    if (arg_errors) { 
        arg_print_errors(stderr, end, "tetracalc");
        arg_print_syntaxv(stderr, argtable, "\n");
        exit(EXIT_FAILURE);
    }

    if (min_len->count) {
        if (min_len->ival[0] > 0) {
            length_cutoff = min_len->ival[0];
        } else {
            fprintf(stderr, "Parameter min_len (-m) must be > 0\n");
            exit(EXIT_FAILURE);
        }
    } else { // Default value.
        length_cutoff = 5000;
    }
    
    if (sec_len->count) {
        if (!fixed->count) {
            fprintf(stderr, "Parameter sec_len (-c) may only be used with --fixed\n");
            exit(EXIT_FAILURE);
        }
        if (sec_len->ival[0] > 0) {
            second_length = sec_len->ival[0];
        } else {
            fprintf(stderr, "Parameter sec_len (-c) must be > 0\n");
            exit(EXIT_FAILURE);
        }
    } else { // Default value.
        second_length = length_cutoff;
    }
    
    if (cor_parm->count) {
        if (!fixed->count) {
            fprintf(stderr, "Parameter cor_thresh (-t) may only be used with --fixed\n");
            exit(EXIT_FAILURE);
        }
        if ((cor_parm->dval[0] >= 0.0) && (cor_parm->dval[0] <= 1.0)) {
            cor_thresh = cor_parm->dval[0];
        } else {
            fprintf(stderr, "Parameter cor_thresh (-t) must be [0.0 - 1.0]\n");
            exit(EXIT_FAILURE);
        }
    } else { // Default value.
        cor_thresh = 0.9;
    }
    
    if (sec_parm->count) {
        if (!fixed->count) {
            fprintf(stderr, "Parameter sec_thresh (-s) may only be used with --fixed\n");
            exit(EXIT_FAILURE);
        }
        if ((sec_parm->dval[0] >= 0.0) && (sec_parm->dval[0] <= 1.0)) {
            sec_thresh = sec_parm->dval[0];
        } else {
            fprintf(stderr, "Parameter sec_thresh (-s) must be [0.0 - 1.0]\n");
            exit(EXIT_FAILURE);
        }
    } else { // Default value.
        sec_thresh = cor_thresh;
    }

    if (tri_parm->count) {
        if (!fixed->count) {
            fprintf(stderr, "Parameter tri_thresh (-r) may only be used with --fixed\n");
            exit(EXIT_FAILURE);
        }
        if ((tri_parm->dval[0] >= 0.0) && (tri_parm->dval[0] <= 1.0)) {
            tri_thresh = tri_parm->dval[0];
        } else {
            fprintf(stderr, "Parameter tri_thresh (-r) must be [0.0 - 1.0]\n");
            exit(EXIT_FAILURE);
        }
    } else { // Default value.
        tri_thresh = cor_thresh;
    }

    if (!fixed->count) {
        double merge_target = 0.95;
        if (merge_tar->count) {
            merge_target = merge_tar->dval[0];
        }
        if ((merge_target >= 0.0) && (merge_target <= 1.0)) {
            tetra_merge_tab = interp_tune_table(merge_target, tune_row_vals, tetra_table_vals, num_tune_rows, num_tune_cols);
            tri_merge_tab = interp_tune_table(merge_target, tune_row_vals, tri_table_vals, num_tune_rows, num_tune_cols);
            tri_cutoff_tab = interp_tune_table(0.99, tune_row_vals, tetra_table_vals, num_tune_rows, num_tune_cols);
            tetra_cutoff_tab = interp_tune_table(0.99, tune_row_vals, tri_table_vals, num_tune_rows, num_tune_cols);
            second_length = 0;
        } else {
            fprintf(stderr, "Parameter --merge_tar must be [0.0 - 1.0]\n");
            exit(EXIT_FAILURE);
        }

        // Limit the permissible adaptive thresholds to some minimum values.

        if (cor_parm->count) {
            limit_thresh(tetra_merge_tab, num_tune_cols, cor_thresh);
        }
        
        
        if (tri_parm->count) {
            limit_thresh(tri_merge_tab, num_tune_cols, tri_thresh);
        }
        
        
        if (sec_parm->count) {
            limit_thresh(tri_cutoff_tab, num_tune_cols, sec_thresh);
        }
    }

    if (chunk->count) {
        if (length_cutoff > chunk->ival[0]) {
            fprintf(stderr, "min_len (-m) is greater than chunk size (--chunk).\n");
            exit(EXIT_FAILURE);
        }
    }
    
#ifdef _OPENMP
    if (num_threads->count) {
        omp_set_num_threads(num_threads->ival[0]);
    } else {
        omp_set_num_threads(omp_get_num_procs());
    }
    //    omp_set_num_threads(1);
#endif
    
    ////////////////////////////////////////////////////////////////////////
    // Read in sequence files
    ////////////////////////////////////////////////////////////////////////
    
    if (verbose->count) { 
        fprintf(stderr,"Num files: %d\n", f->count);
    }
    
    for (unsigned int x = 0; x < f->count; x++) {
        
        if (strcmp(f->filename[x], "-")) {
            
            // Try to open the fasta file
            if (!(f_hnd = gzopen(f->filename[x], "r"))) {
                int errorno = 0;
                const char *errstr = gzerror(f_hnd, &errorno);
                fprintf(stderr, "%s -- FASTA file not found!\n Error %d, %s", f->filename[0], errorno, errstr);
                exit(EXIT_FAILURE);
            }
            
        } else {
            // Try to open the fasta file
            if (!(f_hnd = gzdopen(STDIN_FILENO, "r"))) {
                int errorno = 0;
                const char *errstr = gzerror(f_hnd, &errorno);
                fprintf(stderr, "%s -- Error reading FASTA file from stdin (via '-')\n Error %d, %s", f->filename[0], errorno, errstr);
                exit(EXIT_FAILURE);
            }
        }
        
        seq_info *sp = NULL;
        
        if (sf->count) {    // Combine sequences in a multi-FASTA file into a single scaffold

            unsigned int ns = 0;
            unsigned int tl = 0;
            
            seq_ptr = read_fasta_seq(f_hnd);

            // Calc seq stats for whole file.
            LL_FOREACH(seq_ptr,sp) {
                ns++;
                unsigned int len = utstring_len(sp->seq);
                tl += len;
                sp->start = 0;
                sp->stop = len - 1;
                cleanup(utstring_body(sp->seq));  // Uppercase, [^ACGT] --> N
            }
            
            // Enforce the minimum seq length for analysis
            if (tl >= length_cutoff) {
                seq = new_seq_entry();
                seq->seqs = seq_ptr;
                seq->tot_len = tl;
                seq->num_seqs = ns;
                
                // Lop off any directory path in the filename
                char *tmp = strrchr((char *) f->filename[x], '/');
                ss_strcat_utstring(seq->scaffold_name, (tmp) ? tmp+1 : (char *) f->filename[x]);
                LL_PREPEND(seq_list, seq);
                num_seq++;
            } else {
                // Free up all of the sequences
                seq = new_seq_entry();
                seq->seqs = seq_ptr;
                free_seq_entry(seq);
            }
            
        } else { // !(sf->count) // Consider sequences in a multi-FASTA file as independent sequences
            
            seq_ptr = read_fasta_seq(f_hnd);
            
            while (seq_ptr) {
                
                if (utstring_len(seq_ptr->seq) >= length_cutoff) {
                    
                    cleanup(utstring_body(seq_ptr->seq));  // Uppercase, [^ACGT] --> N
                    
                    if (chunk->count) {
                        
                        UT_string *org_seq_name = seq_ptr->seq_name;
                        
                        // Go through all of the loaded sequences and split them up
                        for (int i = 0; i+chunk->ival[0] <= utstring_len(seq_ptr->seq); i+=chunk->ival[0]) {
                            seq = new_seq_entry();
                            seq->tot_len = chunk->ival[0];
                            seq->num_seqs = 1;
              
                            if (i) {
                                if (!(seq->seqs = calloc(1, sizeof(seq_info)))) {
                                    fprintf(stderr, "calloc failed: chunk seq_info creation\n");
                                    exit(EXIT_FAILURE);
                                }
                                seq->seqs->seq = seq_ptr->seq;
                                seq->seqs->next = NULL;
                            } else {
                                seq->seqs = seq_ptr;
                            }
                            
                            seq->seqs->start = i;
                            seq->seqs->stop = i + chunk->ival[0] - 1;

                            utstring_new(seq->seqs->seq_name);
                            utstring_concat(seq->seqs->seq_name, org_seq_name);
                            utstring_printf(seq->seqs->seq_name,".%u.%u",seq->seqs->start,seq->seqs->stop);
                            utstring_concat(seq->scaffold_name, seq->seqs->seq_name);
                            
                            LL_PREPEND(seq_list, seq);
                            num_seq++;
                        }

                        utstring_free(org_seq_name);
                        sp = seq_ptr->next;
                        seq_ptr->next = NULL;
                        seq_ptr = sp;
                        
                    } else {
                        seq = new_seq_entry();
                        seq->seqs = seq_ptr;
                        seq_ptr = seq_ptr->next;
                        seq->seqs->next = NULL;
                        
                        seq->num_seqs = 1;
                        seq->tot_len = utstring_len(seq->seqs->seq);
                        utstring_concat(seq->scaffold_name, seq->seqs->seq_name);
                        seq->seqs->start = 0;
                        seq->seqs->stop = seq->tot_len - 1;

                        LL_PREPEND(seq_list, seq);
                        num_seq++;
                    }
                } else {
                    // Free up this sequence and move on to the next one
                    seq = new_seq_entry();
                    seq->seqs = seq_ptr;
                    seq_ptr = seq_ptr->next;
                    seq->seqs->next = NULL;
                    free_seq_entry(seq);
                }
            }
        }
        
        // Close fasta file
        gzclose(f_hnd);
        
    }
    
    ////////////////////////////////////////////////////////////////////////
    // Prepare sequences
    ////////////////////////////////////////////////////////////////////////
    
    if (!(seq_array = calloc(num_seq, sizeof(seq_entry *)))) {
        fprintf(stderr, "calloc failed: sequence stats creation\n");
        exit(EXIT_FAILURE);        
    }
    
    seq_pptr = seq_array;
    
    LL_FOREACH(seq_list, seq) {
        *seq_pptr++ = seq;
        all_seq_len += seq->tot_len;
    }
    
    ////////////////////////////////////////////////////////////////////////
    // Process sequences for tetranucleotide Z stats
    ////////////////////////////////////////////////////////////////////////
    
    if (verbose->count) { 
        fprintf(stderr, "Processing %lu sequences\n", num_seq);
        fprintf(stderr, "Total length %lu\n", all_seq_len);
    }
    
    if (num_seq < 2) {
        fprintf(stderr, "ERROR: There must be two or more sequences above the minimum length to cluster\n");
        exit(EXIT_FAILURE);
    }
    
#pragma omp parallel for schedule(static) default(none) private(seq,seq_ptr) shared(num_seq,seq_array,stderr,seq_LUT)
    
    for (long unsigned int x = 0; x < num_seq; x++) {

        seq = seq_array[x];
        seq->index = x;
        
//        fprintf(stderr,"Sequence %d %s\n", x, seq->seq_name);            
        
        double expect, front_tri, back_tri, center_di, var;
        unsigned int c1, c2, c3, c4, r1, r2, r3, r4;
        
        LL_FOREACH(seq->seqs, seq_ptr) {
            
            char *c_ptr = utstring_body(seq_ptr->seq) + seq_ptr->start;
            char *c_stop = utstring_body(seq_ptr->seq) + seq_ptr->stop - 2;

//            fprintf(stderr,"Sequence %s\n", seq_ptr->seq_name);    
            
            // Walk the string and accumulate tetra frequencies
            while (c_ptr != c_stop) {
                
                // Indices from the forward bases of the tetranuc
                c1 = seq_LUT[(int)*c_ptr]; 
                c2 = seq_LUT[(int)*(c_ptr+1)];    
                c3 = seq_LUT[(int)*(c_ptr+2)]; 
                c4 = seq_LUT[(int)*(c_ptr+3)];
                
//                fprintf(stderr,"%d%d%d%d\n", c1,c2,c3,c4);
                
                c_ptr++;
                
                // Test for 'N's in this tetra. If so, skip it.
                if ((c1+c2+c3+c4) > max_tetra_sum) {
                    continue;
                }
                
                // Reverse complement tetra indices
                r1 = 3-c4;    
                r2 = 3-c3;    
                r3 = 3-c2;    
                r4 = 3-c1;
                
                // calc tetranuc frequencies forward strand
                seq->tetra[c1][c2][c3][c4]++;
                
                // calc tetranuc frequencies reverse strand
                seq->tetra[r1][r2][r3][r4]++;
                
                // calc front trinuc frequencies forward strand
                seq->tri_f[c1][c2][c3]++;
                
                // calc front trinuc frequencies reverse strand
                seq->tri_f[r1][r2][r3]++;
                
                // calc back trinuc frequencies forward strand
                seq->tri_b[c2][c3][c4]++;
                
                // calc back trinuc frequencies reverse strand
                seq->tri_b[r2][r3][r4]++;
                
                // calc center dinuc frequencies forward strand
                seq->di_c[c2][c3]++;
                
                // calc center dinuc frequencies reverse strand
                seq->di_c[r2][r3]++;
                
                // nuc frequences for center of "front" trinucs
                seq->nuc[c2]++;
                seq->nuc[r2]++;

            }
        }
        
        calc_Z_stats(seq);  // Calculate the sequence Z-Stats
    } 

    /////////////////////////////////////////////////////////////////////////////////////////
    // Calculate an N X N trangular matrix of correlation coeffs from N^2 Z stat regressions
    /////////////////////////////////////////////////////////////////////////////////////////

    long unsigned int num_regressions = ((num_seq * num_seq) - num_seq)/2;
    
    if (verbose->count) { 
        fprintf(stderr,"Processing %lu regressions\n", num_regressions);
    }
    
    if (!(cor_matrix = calloc(num_regressions, sizeof(stat_entry *)))) {
        fprintf(stderr, "calloc failed: correlation matrix.  %lu * %lu byte structs\n", num_regressions, sizeof(stat_entry *));
        exit(EXIT_FAILURE);        
    }
    
    if (!(cor_entries = calloc(num_regressions, sizeof(stat_entry)))) {
        fprintf(stderr, "calloc failed: correlation entries.  %lu * %lu byte structs\n", num_regressions, sizeof(stat_entry));
        exit(EXIT_FAILURE);
    }
    
    long unsigned int matches = 0;
    
    //#pragma omp parallel for schedule(guided) default(none) shared(cor_entries,cor_matrix,num_seq,seq_array,cor_cutoff,tri_cutoff,sec_cutoff,second_length) reduction(+:matches)
#pragma omp parallel for schedule(static) default(none) shared(fixed,tetra_merge_tab,tri_merge_tab,tri_cutoff_tab,tetra_cutoff_tab,num_tune_cols,tune_col_vals,cor_entries,cor_matrix,num_seq,seq_array,cor_thresh,tri_thresh,sec_thresh,second_length) reduction(+:matches)
    for (long unsigned int y = 1; y < num_seq; y++) {
        for (long unsigned int x = 0; x < y; x++) {
            // Pascal's triangle indexing
            long unsigned int idx = ((y*(y-1))/2) + x;
            
            stat_entry *this_cor = &cor_entries[idx];
            cor_matrix[idx] = this_cor;
            
            this_cor->x_seq = seq_array[x];
            this_cor->y_seq = seq_array[y];

            reg_corr(this_cor);

            if (fixed->count) {
                this_cor->cor_thresh = cor_thresh;
                this_cor->sec_thresh = sec_thresh;
                this_cor->tri_thresh = tri_thresh;
                // this_cor->tetra_thresh = sec_thresh;
            } else {
                this_cor->cor_thresh = calc_thresh(this_cor, tune_col_vals, tetra_merge_tab, num_tune_cols);
                this_cor->sec_thresh = calc_thresh(this_cor, tune_col_vals, tri_cutoff_tab, num_tune_cols);
                this_cor->tri_thresh = calc_thresh(this_cor, tune_col_vals, tri_merge_tab, num_tune_cols);
                // this_cor->tetra_thresh = calc_thresh(this_cor, tune_col_vals, tetra_cutoff_tab, num_tune_cols);
            }
            
            // Determine if this the current top correlation match for each sequence

            matches += new_best(this_cor, second_length);
        }
    }

    // Write diagnostic information on the Z-statistics

    if (zstats->count) {
        
        FILE *file = fopen(zstats->filename[0], "w");
        
        if (!file) {
        	    fprintf(stderr, "ERROR opening Z-statistics file for writing: %s", zstats->filename[0]);
                exit(EXIT_FAILURE);
        } else {
            for (long unsigned int y = 0; y < num_seq; y++) {
                fprintf(file,"%s\t",utstring_body(seq_array[y]->scaffold_name)+1);     // clip off leading '>'
            }
            fprintf(stderr,"\n");
            for (long unsigned int x = 0; x < 256; x++) {
                for (long unsigned int y = 0; y < num_seq; y++) {
                    fprintf(file,"%f\t",seq_array[y]->Z[x]);    
                }
                fprintf(file,"\n");
            }
            fclose(file);
        }
    }
    
    if (verbose->count) { 
        fprintf(stderr,"\nMatching pairs: %lu\n", matches);
    }
    
    stat_entry *max_pair = NULL;
 
    do {
    
        num_regressions = ((num_seq * num_seq) - num_seq)/2;
        max_pair = NULL;
        
#pragma omp parallel default(none) shared(merge_tar,tune_col_vals,tetra_merge_tab,tri_cutoff_tab,tri_merge_tab,num_tune_cols,stderr,num_seq,seq_array,cor_matrix,num_regressions,second_length,max_pair)
        {
            stat_entry *local_max_pair = NULL;
            
            // Find the current maximum scoring match
#pragma omp for schedule(static)
            for (long unsigned int idx = 0; idx < num_seq; idx++) {
                
                //            if (seq_array[idx]->best) {
                //                fprintf(stderr,"\n&&&&&& Checking Idx: %lu %f %f\n",idx,seq_array[idx]->best->cor,seq_array[idx]->best->cor3);
                //                fprintf(stderr,"Best IDs: %lu %lu\n",seq_array[idx]->best->x_seq->index,seq_array[idx]->best->y_seq->index);
                //                fprintf(stderr,"Best IDs: %s\n\n%s\n",utstring_body(seq_array[idx]->best->x_seq->scaffold_name),utstring_body(seq_array[idx]->best->y_seq->scaffold_name));
                //            }
                
                // Test code
                //if ((seq_array[idx]->best) && ((seq_array[idx]->best->x_seq->seqs == NULL) || (seq_array[idx]->best->y_seq->seqs == NULL))) {
                //    fprintf(stderr,"\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
                //    fprintf(stderr,"Error in corr: %lu between %lu and %lu, lengths: %u %u \n",idx,seq_array[idx]->best->x_seq->index,seq_array[idx]->best->y_seq->index,seq_array[idx]->best->x_seq->tot_len,seq_array[idx]->best->y_seq->tot_len);
                //    fprintf(stderr,"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n\n");
                //}
                
                // Early out if this sequence doesn't have a best cor
                if (!seq_array[idx]->best) {
                    continue;
                }
                
                if (max_correlation(seq_array[idx]->best, local_max_pair, second_length) == seq_array[idx]->best) {
                    local_max_pair = seq_array[idx]->best;
                }
            }
#pragma omp critical
            max_pair = max_correlation(max_pair, local_max_pair, second_length);
        }
        
        // If this is true, merge these sequences and rinse/repeat
        if (max_pair) {
            
            if (verbose->count) {
                
                /* // Exhaustive nearest-neighbor check code
                num_regressions = ((num_seq * num_seq) - num_seq)/2;
                
                // If everything is working above, this shouldn't produce any output!!!
                    for (long unsigned int idx = 0; idx < num_regressions; idx++) {
                        if (!(max_correlation(max_pair, cor_matrix[idx], second_length) == max_pair)) {
                            
                            // General stats to stderr
                            fprintf(stderr,"@@@@@@@@@ BAD BEST COR @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
                            fprintf(stderr,"%s\n\n", utstring_body(cor_matrix[idx]->x_seq->scaffold_name));
                            fprintf(stderr,"%s\n", utstring_body(cor_matrix[idx]->y_seq->scaffold_name));
                            fprintf(stderr,"%d\t%d\n", cor_matrix[idx]->x_seq->tot_len, cor_matrix[idx]->y_seq->tot_len);
                            fprintf(stderr,"R\tSlope\tIntercept\n");
                            fprintf(stderr,"Tetra: \t%f\t%f\t%f\n", cor_matrix[idx]->cor, cor_matrix[idx]->cor_thresh, cor_matrix[idx]->sec_thresh);
                            fprintf(stderr,"Tri: \t%f\t%f\t%f\n", cor_matrix[idx]->cor3, cor_matrix[idx]->tri_thresh, cor_matrix[idx]->tetra_thresh);
                        }
                }
                
                */
                 
                // General stats to stderr
                fprintf(stderr,"************************************************************************** %lu\n", num_seq);
                fprintf(stderr,"%s %d\n\n", utstring_body(max_pair->x_seq->scaffold_name), max_pair->x_seq->index);
                fprintf(stderr,"%s %d\n", utstring_body(max_pair->y_seq->scaffold_name), max_pair->y_seq->index);
                fprintf(stderr,"%d\t%d\n", max_pair->x_seq->tot_len, max_pair->y_seq->tot_len);
                fprintf(stderr,"R\tSlope\tIntercept\n");
                fprintf(stderr,"Tetra: \t%f\t%f\t%f\n", max_pair->cor, max_pair->cor_thresh, max_pair->sec_thresh);
                // fprintf(stderr,"Tri: \t%f\t%f\t%f\n", max_pair->cor3, max_pair->tri_thresh, max_pair->tetra_thresh);
                fprintf(stderr,"Tri: \t%f\t%f\n", max_pair->cor3, max_pair->tri_thresh);
            }
                
            // Initialize the temp new seq struct
            seq = new_seq_entry();
            memcpy(seq, max_pair->y_seq, sizeof(seq_entry));

            seq_entry *seq2 = max_pair->x_seq;
            
            // Add seqs to the new seq list
            LL_CONCAT(seq->seqs, seq2->seqs);
            
            seq2->seqs = NULL;
            max_pair->y_seq->seqs = NULL;
            max_pair->y_seq->scaffold_name = NULL;
            
            // Make a new scaffold name
            ss_strcat_utstring(seq->scaffold_name, "\n");
            utstring_concat(seq->scaffold_name, seq2->scaffold_name);
            
            seq->tot_len += seq2->tot_len;
            seq->num_seqs += seq2->num_seqs;
            seq->best = NULL;
            seq->index = -1;

            for (unsigned int i = 0; i < 256; i++) {
                // Casting below eliminates need for 3 extra enclosed loops
                (((unsigned int *)seq->tetra)[i]) += (((unsigned int *)seq2->tetra)[i]);
                if (i < 64) {
                    (((unsigned int *)seq->tri_f)[i]) += (((unsigned int *)seq2->tri_f)[i]);
                    (((unsigned int *)seq->tri_b)[i]) += (((unsigned int *)seq2->tri_b)[i]);
                    if (i < 16) {
                        (((unsigned int *)seq->di_c)[i]) += (((unsigned int *)seq2->di_c)[i]);
                        if (i < 4) {
                            seq->nuc[i] += seq2->nuc[i];
                        }
                    }
                }
            }
            
            calc_Z_stats(seq);   // Recalculate the sequence Z-stats
            
            // Now replace the two merged seqs with just this one.
            
            long unsigned int old_idx1 = max_pair->x_seq->index;
            long unsigned int old_idx2 = max_pair->y_seq->index;
            
            // If old_idx1 is the last sequence, then swap them to make the code below simpler.
            if (old_idx1 == num_seq - 1) {
                long unsigned int tmp = old_idx2;
                old_idx2 = old_idx1;
                old_idx1 = tmp;
            }

            free_seq_entry(seq_array[old_idx1]);
            free_seq_entry(seq_array[old_idx2]);

            seq_array[old_idx1] = seq;
            seq->index = old_idx1;

            // Detect case where neither of the merged sequences is the last one in the list
            if (old_idx2 != num_seq - 1) {
                seq_array[old_idx2] = seq_array[num_seq - 1];
                seq_array[old_idx2]->index = old_idx2;
                
                // Move the precalculated entries in the triangular correlation table to the correct new positions
                long unsigned int y = num_seq-1;
                long unsigned int y_pre = (y*(y-1))/2;  // precalculate the 'y' part of the old correlation matrix row
                long unsigned int idx = -1;
                for (long unsigned int x = 0; x < y; x++) {
                    
                    // Pascal's triangle indexing
                    if (x < old_idx2) {
                        idx = ((old_idx2*(old_idx2-1))/2) + x;
                    } else if (x > old_idx2) { // x > old_idx2
                        idx = ((x*(x-1))/2) + old_idx2;
                        // Keep the x and y sequence pointers consistent.
                        cor_matrix[y_pre + x]->y_seq = seq_array[x];
                        cor_matrix[y_pre + x]->x_seq = seq_array[old_idx2];
                    } else continue;
                    
                    // Swap these entries.
                    stat_entry *tmp_ptr = cor_matrix[y_pre + x];
                    cor_matrix[y_pre + x] = cor_matrix[idx];
                    cor_matrix[idx] = tmp_ptr;
                    
                }
                // Free the unused entries??
                // memset(&cor_matrix[y_pre],0,(num_seq-1)*sizeof(stat_entry));
            }
            
            num_seq--;
            
            // Intelligently recalc invalid entries in the correlation table

            int c1 = 0;
            int c2 = 0;
            int c3 = 0;
            int c4 = 0;
            int c5 = 0;
            
            // #pragma omp parallel for schedule(guided) default(none) shared(greedy,stderr,cor_matrix,num_seq,seq_array,old_idx1,global_max_cor,global_max_cor3,max_pair,cor_cutoff,sec_cutoff,tri_cutoff,second_length) reduction(+:c1,c2,c3,c4,c5)
#pragma omp parallel for schedule(static) default(none) shared(verbose,fixed,tetra_merge_tab,tri_merge_tab,tri_cutoff_tab,tetra_cutoff_tab,num_tune_cols,tune_col_vals,greedy,stderr,cor_matrix,num_seq,seq_array,old_idx1,max_pair,cor_thresh,sec_thresh,tri_thresh,second_length) reduction(+:c1,c2,c3,c4,c5)
            for (long unsigned int z = 0; z < num_seq; z++) {
                long unsigned int x = 0;
                long unsigned int y = 0;
                
                if (z < old_idx1) {
                    x = z;
                    y = old_idx1;
                } else if (z > old_idx1) {
                    y = z;
                    x = old_idx1;
                } else continue;
                
                // Pascal's triangle indexing
                long unsigned int idx = ((y*(y-1))/2) + x;
                
                stat_entry *this_cor = cor_matrix[idx];
                stat_entry *other_cor = cor_matrix[((num_seq*(num_seq-1))/2) + z];
                stat_entry *prev_best = NULL;
                
                int merged_best = 0;
                //double prev_cor = -1.0;
                //double prev_cor3 = -1.0;
                //double prev_cor_thresh = 0.0;
                
                if (seq_array[z]->best == this_cor) {  // If sequence z's previous best was the merged seq old_idx1
                    prev_best = other_cor;
                    memcpy(prev_best, this_cor, sizeof(stat_entry));
                    // prev_cor = this_cor->cor;
                    // prev_cor3 = this_cor->cor3;
                    // prev_cor_thresh = this_cor->cor_thresh;
                    merged_best = 1;
                } else if (seq_array[z]->best == other_cor) {  // If sequence z's previous best was the merged seq old_idx2 (now moved to edge of matrix)
                    prev_best = other_cor;
                    // prev_cor = other_cor->cor;
                    // prev_cor3 = other_cor->cor3;
                    // prev_cor_thresh = other_cor->cor_thresh;
                    merged_best = 2;
                }
                
                this_cor->x_seq = seq_array[x];
                this_cor->y_seq = seq_array[y];
                
                reg_corr(this_cor);
                
                if (fixed->count) {
                    this_cor->cor_thresh = cor_thresh;
                    this_cor->sec_thresh = sec_thresh;
                    this_cor->tri_thresh = tri_thresh;
                    // this_cor->tetra_thresh = sec_thresh;
                } else {
                    this_cor->cor_thresh = calc_thresh(this_cor, tune_col_vals, tetra_merge_tab, num_tune_cols);
                    this_cor->sec_thresh = calc_thresh(this_cor, tune_col_vals, tri_cutoff_tab, num_tune_cols);
                    this_cor->tri_thresh = calc_thresh(this_cor, tune_col_vals, tri_merge_tab, num_tune_cols);
                    // this_cor->tetra_thresh = calc_thresh(this_cor, tune_col_vals, tetra_cutoff_tab, num_tune_cols);
                }
                
                // Determine if this the current top correlation match for each sequence

                // The case below...
                // 1)  If the seq best was one of the merged sequences. and this correlation isn't higher (ie it got worse, unlikely...), then a new high for that sequence needs to be found from among all other correlations.
                // 1a) *** ---> Greedy Shortcut... Just go with the merged seqs as the best anyway... Done.  Evaluate this case!

                // Previous best match was with one of the merged sequences...
                
                if (merged_best) {
                    c1++;
                    // Recalculate the new best correlation match for this sequence using all other sequence pairs.  (This should be the exception...)
                    
                    if ((greedy->count) || (max_correlation(this_cor, prev_best, second_length) == this_cor)) {
                            seq_array[z]->best = NULL;
                            new_best(this_cor, second_length);
                            c2++;
                        } else {
                            c3++;
                            seq_array[z]->best = NULL;
                            
                            for (long unsigned int zidx = 0; zidx < num_seq; zidx++) {
                                long unsigned int xx = 0;
                                long unsigned int yy = 0;
                                
                                if (zidx < z) {
                                    xx = zidx;
                                    yy = z;
                                } else if (zidx > z) {
                                    yy = zidx;
                                    xx = z;
                                } else continue;
                                
                                // Pascal's triangle indexing
                                long unsigned int iidx = ((yy*(yy-1))/2) + xx;
                                new_best(cor_matrix[iidx], second_length);
                            }
                            
                            if (seq_array[z]->best != this_cor) {
                                c4++;
                                if (verbose->count) {
                                    if (seq_array[z]->best) {
                                        fprintf(stderr,"^^^^^^^^^^^^^^^^^^^^^ Full Recalc case: %lu %f %f %e\n",z,this_cor->cor,seq_array[z]->best->cor,seq_array[z]->best->cor-this_cor->cor);
                                    } else {
                                        fprintf(stderr,"^^^^^^^^^^^^^^^^^^^^^ No best found! %lu %f %f \n",z,this_cor->cor,this_cor->cor3);
                                    }
                                }
                            }
                        }
                } else {
                    c5++;
                    new_best(this_cor, second_length);
                }
            }
            if (verbose->count) {
                fprintf(stderr,"++++++++++++++ Merged best: %d Greedy: %d Recalc: %d Better: %d Best not merged: %d %.2f%% %.2f%% %.2f%%\n",c1,c2,c3,c4,c5, 100.0*c1/(c1+c5), 100.0*c2/(c2+c3),100.0*c4/c3);
            }
        }
    } while (max_pair);
    
    
    ////////////////////////////////////////////////////////////////////////
    // Output merged results!
    
    num_regressions = ((num_seq * num_seq) - num_seq)/2;

    // If everything is working above, this shouldn't produce any output!!!
    if (verbose->count) {
        for (long unsigned int idx = 0; idx < num_regressions; idx++) {
            if (pass_correlation(cor_matrix[idx], second_length)) {
                
                // General stats to stderr
                fprintf(stderr,"------------------------------------------------------------\n");
                fprintf(stderr,"%s\n\n", utstring_body(cor_matrix[idx]->x_seq->scaffold_name));
                fprintf(stderr,"%s\n", utstring_body(cor_matrix[idx]->y_seq->scaffold_name));
                fprintf(stderr,"%d\t%d\n", cor_matrix[idx]->x_seq->tot_len, cor_matrix[idx]->y_seq->tot_len);
                fprintf(stderr,"R\tSlope\tIntercept\n");
                fprintf(stderr,"Tetra: \t%f\t%f\t%f\n", cor_matrix[idx]->cor, cor_matrix[idx]->cor_thresh, cor_matrix[idx]->sec_thresh);
                // fprintf(stderr,"Tri: \t%f\t%f\t%f\n", cor_matrix[idx]->cor3, cor_matrix[idx]->tri_thresh, cor_matrix[idx]->tetra_thresh);
                fprintf(stderr,"Tri: \t%f\t%f\n", cor_matrix[idx]->cor3, cor_matrix[idx]->tri_thresh);
            }
        }
    }

    if (verbose->count) { 
        fprintf(stderr,"\nNumber of unmerged sequences remaining: %d\n", num_seq);
    }
    
#pragma omp parallel for schedule(static) default(none) shared(stderr,cor_entries,cor_matrix,num_seq,seq_array,second_length)

    for (long unsigned int y = 1; y < num_seq; y++) {
        for (long unsigned int x = 0; x < y; x++) {
            // Pascal's triangle indexing
            long unsigned int idx = ((y*(y-1))/2) + x;
            
            stat_entry * this_cor = cor_matrix[idx];
            
            this_cor->x_seq = seq_array[x];
            this_cor->y_seq = seq_array[y];
            
            // Determine if this is the current top correlation match for each sequence
            
            // if (!new_best(this_cor, -1.0, -1.0, -1.0, 0)) {
            //     fprintf(stderr, "\nNonmatching %lu %lu\n", y, x);
            // }
            
        }
    }

    if (verbose->count) {
        fprintf(stderr, "\nRemaining sequences:\n");
    }
    
    qsort(seq_array, num_seq, sizeof(seq_entry *), seq_entry_length_cmp);
    
    if (verbose->count) {
        fprintf(stderr,"****RESULTS****\n");
        if (num_seq == 1) {
            fprintf(stderr,"--- %lu %u %lu\t\n%s\n", 0, seq_array[0]->tot_len, all_seq_len, utstring_body(seq_array[0]->scaffold_name));
        } else {
            for (long unsigned int idx = 0; idx < num_seq; idx++) {
                if (seq_array[idx]->best) {
                fprintf(stderr,"--- %lu %u %.3f %.3f %u %.3f %.3f %.3f\t\n%s\n", idx, seq_array[idx]->tot_len, seq_array[idx]->best->cor, seq_array[idx]->best->cor3, (seq_array[idx]->best->x_seq->index == idx) ? seq_array[idx]->best->y_seq->index : seq_array[idx]->best->x_seq->index, seq_array[idx]->best->cor_thresh, seq_array[idx]->best->sec_thresh, seq_array[idx]->best->tri_thresh, utstring_body(seq_array[idx]->scaffold_name));
                } else {
                    fprintf(stderr,"--- %lu %u\n%s\n", idx, seq_array[idx]->tot_len, utstring_body(seq_array[idx]->scaffold_name));
                }
            }
        }
    }    

    // Write output in JSON format on STDOUT
    
    JSON_BEG_OBJ;
    JSON_SUB_PROP("clusters");
    JSON_BEG_ARR;
    int prev_flag = 0;
    for (long unsigned int idx = 0; idx < num_seq; idx++) {
        if (seq_array[idx]->seqs->next) {
            if (prev_flag) {
                JSON_COMMA;
            } else {
               prev_flag = 1;
            }
            JSON_BEG_ARR;
            seq_info *sp = NULL;
            LL_FOREACH(seq_array[idx]->seqs, sp) {
                JSON_STR_ELEM(utstring_body(sp->seq_name)+1);
                if (sp->next) {
                    JSON_COMMA;
                }
            }
            JSON_END_ARR;
        }
    }
    JSON_END_ARR;
    JSON_END_OBJ;
    printf("\n");
    
    ////////////////////////////////////////////////////////////////////////
    // Free sequences
    ////////////////////////////////////////////////////////////////////////
    
    free(cor_matrix);
    
    free(cor_entries);
    
    for (long unsigned int x = 0; x < num_seq; x++) {
        free_seq_entry(seq_array[x]); 
    } 
    
    free(seq_array);
    
    if (!fixed->count) {
        free(tetra_merge_tab);
        free(tri_cutoff_tab);
        free(tri_merge_tab);
    }
    
    exit(EXIT_SUCCESS);
}

////////////////////////////////////////////////////////////////////////
// Function: reg_corr - Regression/Correlation analysis of two Z arrays 
// Returns: a stats structure 

void reg_corr(stat_entry *s_ptr) {
    
    double *xZ = s_ptr->x_seq->Z;
    double *yZ = s_ptr->y_seq->Z;

    double sum_x = 0.0;
    double sum_y = 0.0;
    double sum_xy = 0.0;
    double sum_x2 = 0.0;
    double sum_y2 = 0.0;
    double N = 256.0;

    for (int z = 0; z < 256; z++) {
        double zX = *xZ++;
        double zY = *yZ++;
        
        sum_x += zX;
        sum_x2 += zX * zX;
        sum_xy += zX * zY;
        sum_y2 += zY * zY;
        sum_y += zY;
    }
    
    // Calculate regression/correlation parameters
    double sum_x_2 = sum_x * sum_x;
    double sum_y_2 = sum_y * sum_y;
    
    // s_ptr->slope = ((N * sum_xy) - (sum_x * sum_y)) / ((N * sum_x2) - (sum_x_2));
    // s_ptr->intercept = (sum_y / N) - (s_ptr->slope * (sum_x / N));
    
    double slope = ((N * sum_xy) - (sum_x * sum_y)) / ((N * sum_x2) - (sum_x_2));
    double intercept = (sum_y / N) - (slope * (sum_x / N));
    
//    s_ptr->cor = s_ptr->slope * sqrt(((N * sum_x2) - sum_x_2) / ((N * sum_y2) - sum_y_2));
    s_ptr->cor = slope * sqrt(((N * sum_x2) - sum_x_2) / ((N * sum_y2) - sum_y_2));
    
    // Now do trinucleotides too
    
    xZ = s_ptr->x_seq->Z3;
    yZ = s_ptr->y_seq->Z3;
    
    sum_x = 0.0;
    sum_y = 0.0;
    sum_xy = 0.0;
    sum_x2 = 0.0;
    sum_y2 = 0.0;
    N = 64.0;
    
    for (int z = 0; z < 64; z++) {
        double zX = *xZ++;
        double zY = *yZ++;
        
        sum_x += zX;
        sum_x2 += zX * zX;
        sum_xy += zX * zY;
        sum_y2 += zY * zY;
        sum_y += zY;
    }
    
    // Calculate regression/correlation parameters
    sum_x_2 = sum_x * sum_x;
    sum_y_2 = sum_y * sum_y;
    
    // s_ptr->slope3 = ((N * sum_xy) - (sum_x * sum_y)) / ((N * sum_x2) - (sum_x_2));
    // s_ptr->intercept3 = (sum_y / N) - (s_ptr->slope3 * (sum_x / N));
    
    // s_ptr->cor3 = s_ptr->slope3 * sqrt(((N * sum_x2) - sum_x_2) / ((N * sum_y2) - sum_y_2));

    double slope3 = ((N * sum_xy) - (sum_x * sum_y)) / ((N * sum_x2) - (sum_x_2));
    double intercept3 = (sum_y / N) - (slope3 * (sum_x / N));
    
    s_ptr->cor3 = slope3 * sqrt(((N * sum_x2) - sum_x_2) / ((N * sum_y2) - sum_y_2));

    return;
}

////////////////////////////////////////////////////////////////////////
// Function: new_seq_entry 
// Returns: an alloc'ed and zeroed seq_entry struct

seq_entry *new_seq_entry() {
    seq_entry *seq = NULL;
    
    if (!(seq = calloc(1, sizeof(seq_entry)))) {  
        fprintf(stderr, "calloc failed: sequence stats creation\n");
        exit(EXIT_FAILURE);        
    }
    
    utstring_new(seq->scaffold_name);
    return(seq);
}
               
////////////////////////////////////////////////////////////////////////
// Function: free_seq_entry
// Free sequence struct and alloced strings within
// Returns: void

void free_seq_entry(seq_entry *seq) {
    
    seq_info *s_ptr = seq->seqs;
    seq_info *s_ptr2 = NULL;

    // Free up sequence info
    while (s_ptr) {
        s_ptr2 = s_ptr;
        s_ptr = s_ptr2->next;
        utstring_free(s_ptr2->seq_name);
        // Only free this for chunk 0.  Other chunks share the seq string.
        if (s_ptr2->start == 0) {
            utstring_free(s_ptr2->seq);
        }
        free(s_ptr2);
    };
    
    if (seq->scaffold_name)
        utstring_free(seq->scaffold_name);
    
    free(seq);
}
               
////////////////////////////////////////////////////////////////////////
// Function: cleanup 
// Returns cleaned-up input with non-[atgcATGC] turned to N,
// and all characters are uppercased in output
// Note, that the counts array is not initialized, so this function
// can be used to accumulate residue counts.  It must be zero
// pre-initialized before the first time it is used with cleanup() 
//  

void cleanup(char *sequence) {
    
    // LUT for transliteration
    static const char clean[] = 
    "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNANCNNNGNNNNNNNNNNNNTNNNNNNNNNNNNANCNNNGNNNNNNNNNNNNTNNNNNNNNNNN";
    
    // Walk through the input buffer, counting base occurrances and 
    // cleaning up the characters using the "clean" LUT.  Anything not
    // in [ATGCatgc] becomes 'N'.  [atgc] become [ATGC].
    
    while (*sequence) {
        *sequence = clean[(int) *sequence]; 
        sequence++;
    }
    
    return;
}

////////////////////////////////////////////////////////////////////////
// Function: read_fasta_seq 
// Takes an open filepointer, and pointers to sequence, header buffer,
// and coding position buffer pointers.
// 
// Returns allocated and filled buffers in each of the pointer locations.
// Return value is 0 when EOF or FASTA format error occurs
//  

seq_info *read_fasta_seq(gzFile *fn)    // Pointer to open file descriptor
{
    seq_info *seqs = NULL;
    seq_info *seq_ptr = NULL;
    z_off_t fpos = 0;        // File position, used to read next sequence

    UT_string *buf = NULL;        // Temporary input buffer for file input
    utstring_new(buf);
    
    do {
        
        if (!(seq_ptr =  calloc(1, sizeof(seq_info)))) {
            fprintf(stderr, "calloc failed: sequence info creation\n");
            exit(EXIT_FAILURE);
        }
        
        utstring_new(seq_ptr->seq_name);
        utstring_new(seq_ptr->seq);
        
        // Now read sequence data, first char of each line must be alphabetic.
        // Stop at EOF or line with a non-alpha prefix

        // If the next header isn't already in the read buffer
        if (utstring_body(buf)[0] != '>') {
            // Eat any lines before header, then grab header
            while ((ss_gzget_utstring(fn, seq_ptr->seq_name)) && (utstring_body(seq_ptr->seq_name)[0] != '>'));
        } else {
            utstring_concat(seq_ptr->seq_name, buf);    // Just use the header we already have
        }
        
        ss_trim_utstring(seq_ptr->seq_name, strcspn(utstring_body(seq_ptr->seq_name), "\n\r"));
        
        // Now loop to read in file data
        while ((ss_gzget_utstring(fn, buf)) && isalpha(utstring_body(buf)[0])) {
            // trim the string at first whitespace
            ss_trim_utstring(buf, strcspn(utstring_body(buf), " \t\n\r"));
            utstring_concat(seq_ptr->seq, buf);
        }
        
        // If not the end of file, restore the file pointer by one line
        if (utstring_len(seq_ptr->seq)) {
            LL_PREPEND(seqs,seq_ptr);
        }
        
    } while (utstring_len(seq_ptr->seq));
        
    utstring_free(seq_ptr->seq_name);
    utstring_free(seq_ptr->seq);
    free(seq_ptr);
        
    utstring_free(buf);

    return(seqs);  //  Done with this file
}

///////////////////////////////////////////////////////////////////
// Calculate the "Z-statistics" for one sequence
///////////////////////////////////////////////////////////////////

void calc_Z_stats(seq_entry *seq) {

    double expect, front_tri, back_tri, center_di, var;
    int c1, c2, c3, c4;
    
    for (c1 = 0; c1 < 4; c1++) {
        for (c2 = 0; c2 < 4; c2++) {
            for (c3 = 0; c3 < 4; c3++) {
                for (c4 = 0; c4 < 4; c4++) {
                    // Equation from: Teeling, et al.
                    // Application of tetranucleotide frequencies for the assignment of genomic fragments
                    // Environmental Microbiology (2004) 6(9), 938–947
                    front_tri = seq->tri_f[c1][c2][c3];
                    back_tri = seq->tri_b[c2][c3][c4];
                    center_di = seq->di_c[c2][c3];
                    expect = (front_tri * back_tri) / center_di;
                    var = expect * (((center_di - front_tri) * (center_di - back_tri)) / (center_di * center_di));
                    seq->Z[c1*64+c2*16+c3*4+c4] = ((double) seq->tetra[c1][c2][c3][c4] - expect) / sqrt(var);
                    if (isnan(seq->Z[c1*64+c2*16+c3*4+c4])) {
                        seq->Z[c1*64+c2*16+c3*4+c4] = 0.0;
                    }
                }
                // Similar calculation as above, but using trinucleotides instead, hence "Z3"
                expect = ((double) seq->di_c[c1][c2] * (double) seq->di_c[c2][c3]) / (double) seq->nuc[c2];
                var = expect * ((((double) seq->nuc[c2] - (double) seq->di_c[c1][c2]) * ((double) seq->nuc[c2] - (double) seq->di_c[c2][c3])) / ((double) seq->nuc[c2] * (double) seq->nuc[c2]));
                seq->Z3[c1*16+c2*4+c3] = ((double) seq->tri_f[c1][c2][c3] - expect) / sqrt(var);
                if (isnan(seq->Z3[c1*16+c2*4+c3])) {
                    seq->Z3[c1*16+c2*4+c3] = 0.0;
                }
            }
        }
    }

    return;
}

///////////////////////////////////////////////////////////////////
// Calculate which of two possible merges is the best
///////////////////////////////////////////////////////////////////

int new_best(stat_entry *this_cor, unsigned int second_length) {
    
    // Determine if this the current top correlation match for each sequence
    
    int matches = 0;

    if (max_correlation(this_cor, this_cor->x_seq->best, second_length) == this_cor) {
#pragma omp critical (new_best)
        this_cor->x_seq->best = max_correlation(this_cor, this_cor->x_seq->best, second_length);
    }
        
    if (max_correlation(this_cor, this_cor->y_seq->best, second_length) == this_cor) {
#pragma omp critical (new_best)
        this_cor->y_seq->best = max_correlation(this_cor, this_cor->y_seq->best, second_length);
    }
    
    return matches;
}

/////////////////////////////////////////////////////////////////////////
// Sort cmp function for qsorting output sequence clusters by length
/////////////////////////////////////////////////////////////////////////

int seq_entry_length_cmp(const void *ptr2, const void *ptr1) {
    
    int retval = 0;
    int sc = (*((seq_entry **) ptr1))->tot_len - (*((seq_entry **) ptr2))->tot_len;
    retval = (sc > 0) ? 1 : -1;
    
    return retval;
}

/////////////////////////////////////////////////////////////////////////
// Interpolate a merge-confidence specific lookup table from the tetra
// or tri nucleotide specific
/////////////////////////////////////////////////////////////////////////

double *interp_tune_table(double value, float *row_tab, float *val_tab, int tab_rows, int tab_cols) {
    
    double *out_tab = NULL;
    
    if (!(out_tab =  calloc(tab_cols, sizeof(double)))) {
        fprintf(stderr, "calloc failed: interpolation table\n");
        exit(EXIT_FAILURE);
    }
    
    for (int c = 0; c < tab_cols; c++) {
        out_tab[c] = -1.0;
        for (int r = 1; r < tab_rows; r++) {
            float v1 = *(val_tab+r*tab_cols+c);
            float v0 = *(val_tab+(r-1)*tab_cols+c);
            if ((v1 <= value) && (v0 >= value)) {
                out_tab[c] = (double) row_tab[r-1] + ((row_tab[r] - row_tab[r-1]) * ((v0 - value) / (v0 - v1)));
                break;
            }
        }
        if (out_tab[c] == -1.0)
            out_tab[c] = row_tab[tab_rows-1];
    }
    
    return out_tab;
}

/////////////////////////////////////////////////////////////////////////
// Place a floor on the thresholds used by smaller sized sequences.
/////////////////////////////////////////////////////////////////////////

int limit_thresh(double *thresh_tab, int tab_cols, double limit) {

    int limited = 0;
    for (int c = 0; c < tab_cols; c++) {
        if (thresh_tab[c] < limit) {
            thresh_tab[c] = limit;
            limited++;
        }
    }
    return limited;
}

/////////////////////////////////////////////////////////////////////////
// Calculate the threshold for a given scaffold size from an interp
// table calculated by interp_tune_table (above).
/////////////////////////////////////////////////////////////////////////

double calc_thresh(stat_entry *cor, float *col_tab, double *thresh_tab, int tab_cols) {
    
    double output = 0.0;
    
    //double min = 0;
    //double max = 0;
    
    unsigned int size = (cor->x_seq->tot_len > cor->y_seq->tot_len) ? cor->y_seq->tot_len : cor->x_seq->tot_len;
    
    //if (cor->x_seq->tot_len > cor->y_seq->tot_len) {
    //    min = log2(cor->y_seq->tot_len);
    //    max = log2(cor->x_seq->tot_len);
    //} else {
    //    max = log2(cor->y_seq->tot_len);
    //    min = log2(cor->x_seq->tot_len);
    //}
    
    //unsigned int size = exp2(min + (max-min)/4.0);
    
    // Entries in the col
    double log_min = col_tab[0];
    double log_max = col_tab[tab_cols-1];
    double log_step = (log_max - log_min) / (tab_cols-1.0);
    double log_size = log2(size);
    
    if (log_size < log_min) {
        output = thresh_tab[0];
    } else if (log_size >= log_max) {
        output = thresh_tab[tab_cols-1];
    } else { // Interpolate in log space
        int pos = (int) trunc((log_size - log_min) / log_step);
        output = thresh_tab[pos] + ((thresh_tab[pos+1] - thresh_tab[pos]) * (log_size - col_tab[pos])/log_step);
    }
    
    return output;
}

/////////////////////////////////////////////////////////////////////////
// Checks to see if a given correlation pair passes either of its
// merge threshold conditions (tetranuc or trinuc)
/////////////////////////////////////////////////////////////////////////
double pass_correlation(stat_entry *pair, unsigned int second_length) {
    
    double pass = 0.0;
    
    if ((pair) && (pair->cor3 >= pair->tri_thresh)) {
        if (pair->cor >= pair->cor_thresh) {
            pass = pair->cor;
        } else if ((pair->cor >= pair->sec_thresh) &&
                   ((!second_length) ||
                    (pair->x_seq->tot_len <= second_length) ||
                    (pair->y_seq->tot_len <= second_length)))
        {
            pass = -pair->cor3;
            //pass = -pair->cor;
        }
    }
    
    return pass;
}

/////////////////////////////////////////////////////////////////////////
// Returns a pointer to whichever of the passed in correlation pairs
// has the greater score.  Returns NULL if neither has a passing correlation
/////////////////////////////////////////////////////////////////////////
stat_entry *max_correlation(stat_entry *pair1, stat_entry *pair2, unsigned int second_length) {
    
    stat_entry *max = NULL;
    
    double pass_p1 = pass_correlation(pair1, second_length);
    double pass_p2 = pass_correlation(pair2, second_length);
    
    if (pass_p1 > 0.0) {
        if (pass_p2 > pass_p1) {
            max = pair2;
        } else {
            max = pair1;
        }
    } else if (pass_p2 > 0.0) {
        max = pair2;
    } else if (pass_p1 < 0.0) {
        if (pass_p2 < pass_p1) {
            max = pair2;
        } else {
            max = pair1;
        }
    } else if (pass_p2 < 0.0) {
        max = pair2;
    }
    
    return max;
}
