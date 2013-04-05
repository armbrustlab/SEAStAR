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
 Name        : ref_select
 Description : Reads SAM format alignment files and a referece seq
 dataset and selects the most likely reference sequences, plus
 a bunch of other stuff...
 ============================================================================
*/

#include <seastar_shared.h>

#include <assert.h>
#include <argtable2.h>
#include <uthash.h>
#include <utlist.h>

// This maintains a version of the program that builds and runs correctly
// when OpenMP support is disabled in--or absent from--the compiler.
#ifdef _OPENMP
#include <omp.h>
#define UNTHREAD 0  // You can unthread reads separately from OpenMP here
#else
//#error "No OpenMP !"
#define UNTHREAD 1     // This always needs to be = 1
#endif

/////////////////////////
// Macros 
/////////////////////////

// Hash FOREACH macro.  Note, this will not work for structs which can
// be present in multiple hashes!!! (ie. hh is not their only UT_hash_handle)
#define ss_HASH_FOREACH(head,el) for(el=head;el;el=el->hh.next)

#define ss_HASH_ADD_UTSTR_PTR(head,utstr_field,add)                                \
HASH_ADD_KEYPTR(hh,head,utstring_body((add)->utstr_field),utstring_len((add)->utstr_field),add)

#define ss_HASH_ADD_UTSTR_STRUCT(head,utstr_field,add)                             \
HASH_ADD_KEYPTR(hh,head,utstring_body(&(add)->utstr_field),utstring_len(&(add)->utstr_field),add)

/////////////////////////
// Typedefs 
/////////////////////////

typedef struct pairing_problem_st {
    unsigned int start;
    unsigned int end;
    struct pairing_problem_st *next;  // Note, this must be called "next", for utlist.h
    char type;
} pairing_problem;

typedef struct read_list_st {
    unsigned int r_id;
    int seq_pos;
    unsigned int map_count;
    struct seq_st *seq_ptr;
    UT_hash_handle hh;
    struct read_list_st *prev;    // Note, this must be called "prev", for utlist.h
    struct read_list_st *next;    // Note, this must be called "next", for utlist.h
} read_list_entry;

// This structure is used to build the ref seq catalog table
typedef struct ref_seq_st {
    unsigned int ref_seq_len;
    unsigned int uncovered;
    float ref_bit_score;
    float ref_coverage; 
    float ref_gc_percent;
    double ref_abundance;
    UT_hash_handle hh;
    UT_string seq_name;
    UT_string seq_desc;
    struct seq_st *sequences;
} ref_seq_entry;

// This structure is used to build the seq catalog hash table
typedef struct seq_st {
    unsigned int seq_len;
    unsigned int offset;    // Holds the offset from this chunk to the beginning of sequence
    int seq_num;
    float bit_score;
    float abs_bit_score;    // Only used when length normalized bitscoring is in use. Used to accurately recalculate bit_score
    float org_bit_score;
    float coverage;
    float gc_percent;
    float mp_insert_mean;
    float mp_insert_stdev;
    unsigned int mp_good_pairs;
    double abundance;
    ref_seq_entry *ref_seq;
    UT_hash_handle hh;
    read_list_entry *reads;
    UT_string seq_id;
    float **recon_array;    // Pointer to array of arrays of floats (coverage, and sequence errors)
    struct seq_st *prev;    // Note, this must be called "prev", for utlist.h
    struct seq_st *next;    // Note, this must be called "next", for utlist.h
    struct seq_list_st *shared;
    struct seq_list_st *forward;
    struct seq_list_st *backward;
    int num_shared;
    int num_forward;
    int num_backward;
    float uncovered;
    float adj_coverage_mean;
    float adj_coverage_stddev;
    float adj_coverage_min;
    float adj_coverage_max;
    pairing_problem *mp_issues;
} seq_entry;

// Read name to id mapping
typedef struct read_id_st {
    UT_string read_name;
    int read_id;
    UT_hash_handle hh;                
} read_name_id;

// The above struct is necessary because the struct below is not hashable due
// to the fact that it is alloced/indexed in an array that may need to be realloced
// as it grows, which would mess up all of the pointers of a hash.

// This structure is used to build a table of reads
typedef struct read_st {
    float bit_info;
    unsigned short int seq_count;
    unsigned short int read_len;
    unsigned short int edit_dist;
    unsigned short int second_edit_dist;
    unsigned int mate_id;
    read_name_id *read_name_ptr;
    read_list_entry *read_chain;
    read_list_entry *second_read_chain;
    char read_type;  // C or N for Color or Nucleotide
} read_entry;

// Exclude seq table entry
typedef struct exclude_st {
    UT_string seq_id;
    UT_hash_handle hh;                
} exclude_entry;

// Detail seq table entry
typedef struct detail_st {
    UT_string seq_name;
    UT_hash_handle hh;                
} detail_entry;

// seq hash entry
typedef struct seq_hash_st {
    float shared_bitscore;
    int seq_num;
    int num_reads;
    int pos_sum;
    UT_hash_handle hh;                
} seq_hash_entry;

// seq list entry
typedef struct seq_list_st {
    float shared_bitscore;
    int seq_num;
    int mean_pos;
    seq_entry *seq_ptr;
} seq_list_entry;

/////////////////////////
// Globals 
/////////////////////////
// NONE!

/////////////////////////
// Function prototypes
/////////////////////////

// CMP functions for qsort
int read_entry_count_cmp(const void *ptr1, const void *ptr2);

// This is for debugging VSI
//int read_entry_id_cmp(const void *ptr1, const void *ptr2);

int catalog_entry_bitscore_cmp(const void *ptr2, const void *ptr1); 
int selected_list_bitscore_cmp(const void *ptr2, const void *ptr1);
int seq_entry_bitscore_cmp(const void *ptr2, const void *ptr1); 

long unsigned int sam_stream_reader(UT_string *fn, int pipe_fd);
long unsigned int sam_header_stream_reader(UT_string *fn, int pipe_fd);

int create_catalog_entry(char *seq_id_str, int seq_len, char *ref_seq_id_str,
                         char *desc_str, exclude_entry **exclude_table,
                         detail_entry **detail_table, seq_entry **seq_table,
                         ref_seq_entry **catalog, int invert_ex,
                         int sep, int det_all, int split_len);

/////////////////////////////////////////////////////////////////
// Entry point
/////////////////////////////////////////////////////////////////

int main(const int argc, char *argv[]) {
    
    ///////////////////////////
    // Variable declarations
    ///////////////////////////    
    
    // Weight table used by sequence reconstruction algorithm
    float weight[2][16] = {
        {
            1.0e6, 1.0, 1.0, 1.25,
            1.0, 1.25, 1.25, 1.25,
            1.0, 1.25, 1.25, 1.25,
            1.25, 1.25, 1.25, 1.25
        },
        {
            1.0e6, 0.0, 0.0, 0.334,
            0.0, 0.334, 0.334, 0.75,
            0.0, 0.334, 0.334, 0.75,
            0.334, 0.75, 0.75, 0.99
        }
    };     
    
    // LUT for transliteration
    static const char code[] = 
    "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\x01N\x02NNN\x04NNNNNN\x00NNNNN\x08NNNNNNNNNNNN\x01N\x02NNN\x04NNNNNN\x00NNNNN\x08NNNNNNNNNNN";
    
    // Build table for colorspace decoder ring
    char c2s[128][128];
    memset(c2s, 0, 128*128);
    c2s['N']['A'] = 'N'; c2s['N']['C'] = 'N'; c2s['N']['G'] = 'N'; c2s['N']['T'] = 'N'; c2s['N']['N'] = 'N';
    c2s['A']['A'] = 'A'; c2s['A']['C'] = 'C'; c2s['A']['G'] = 'G'; c2s['A']['T'] = 'T';
    c2s['C']['A'] = 'C'; c2s['C']['C'] = 'A'; c2s['C']['G'] = 'T'; c2s['C']['T'] = 'G';
    c2s['G']['A'] = 'G'; c2s['G']['C'] = 'T'; c2s['G']['G'] = 'A'; c2s['G']['T'] = 'C';
    c2s['T']['A'] = 'T'; c2s['T']['C'] = 'G'; c2s['T']['G'] = 'C'; c2s['T']['T'] = 'A';    
    
    // Create ambiguity bitcode matching table
    char match[16][16];
    memset(match, 0, 16*16);
    match[1][1] = 1; match[2][1] = 2; match[4][1] = 4; match[8][1] = 8;
    match[1][2] = 2; match[2][2] = 1; match[4][2] = 8; match[8][2] = 4;
    match[1][4] = 4; match[2][4] = 8; match[4][4] = 1; match[8][4] = 2;
    match[1][8] = 8; match[2][8] = 4; match[4][8] = 2; match[8][8] = 1;
    
    // Calculate remainder of table
    for (int i = 1; i < 16; i++) {
        for (int j = 1; j < 16; j++) {
            if (!match[i][j]) {
                for (int k = 0; k < 4; k++) {
                    for (int l = 0; l < 4; l++) {
                        match[i][j] |= match[i&(1<<k)][j&(1<<l)];
                    }
                }
            }
        }
    }
    
    // Nucleotide code complement table
    int comp[16] =   { 0, 8, 4, 0, 2, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0 };
    
    // GC content calc table
    int gc_tab[16] = { 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }; 
    
    // Ambiguity code lookup table
    char cc[] = "XACMGRSVTWYHKDBN";
    
    // Used in postprocessing
    int max_sel_seq_len = 0;  
    
    // Misc working char/str variables
    char *ptr = NULL, *ptr2 = NULL, *tokptr = NULL;
    
    UT_string *str, *str2, *str3, *str4;
    UT_string *name_str;
    UT_string *id_str;
    
    utstring_new(str);
    utstring_new(str2);
    utstring_new(str3);
    utstring_new(str4);
    utstring_new(name_str);
    utstring_new(id_str);
    
    char *s_id = NULL;
    
    // Used to build filenames for outputting selected reads
    UT_string *file_name;
    utstring_new(file_name);    
    
    UT_string *read_out_fn;
    utstring_new(read_out_fn);
    
    // Controls where "second chance" (best edit distance + 1) read mappings
    // are used in abundance estimation when the best mapping taxa have been
    // eliminated (eg. via bit score threshold)
    
    float bit_thresh = 0.00;
    float bit_frac = 0.00;
    
    float share_lim = 500.0;    // Bit score limited for mate-pair analysis output
    float mate_lim = 250.0;
    
    int split_len = 0;
    unsigned int mp_cutoff_len = 0;
    
    int sim_frac_skip_reads = 0;
    int read_cnt_limit = 50000;
    int skip_flag = 0;
    
    // Hash tables
    ref_seq_entry *catalog = NULL, *entry = NULL, **cat_array = NULL, **cat_ptr = NULL;
    seq_entry *seq_table = NULL, *seq = NULL, **seq_array = NULL, **seq_ptr = NULL, *selected = NULL, *s_ptr = NULL;
    exclude_entry *exclude_table = NULL, *ex_ptr = NULL;
    detail_entry *detail_table = NULL, *det_ptr = NULL;
    int det_all = 0;
    seq_hash_entry *shared_seq = NULL;
    
    // Read tables
    
    unsigned long int read_info_len = 0, read_info_inc = 0;
    read_info_len = read_info_inc = 125000;
    read_entry *read_table = NULL, *r_ptr = NULL, *r_ptr2 = NULL;    
    if (!(read_table = calloc(read_info_inc, sizeof(read_entry)))) {  
        fprintf(stderr, "ERROR: calloc failed: read_table array creation\n");
        exit(EXIT_FAILURE);        
    }
    
    read_list_entry *read = NULL, *r2 = NULL, *r3 = NULL;    // Working read entry pointer
    read_name_id *skipped = NULL;    // Reads that have been skipped using --sim_frac.
    read_name_id *read_names = NULL;    // Mapping from name to internal ID
    read_name_id *read_name_ptr = NULL;
    
    // Counters
    unsigned int read_id = 0, new_read_id = 1;    // Used so an array of reads can be used
    
    // Working variables
    unsigned int num = 0;        // Number of refs mapped
    int ed = 0;    // min and current read edit distance
    float num_seq = 0.0;        // Number of ref seqs
    unsigned int n_taxa = 0;
    unsigned int n_seq = 0;
    unsigned int num_reads = 0;
    int position = 0;
    int rnd_seed = 12345;
    ss_rand_inst *rnd_state = NULL;
    
    int cov_window = 0;
    int colorspace_recon = 0;    // Set to one if FASTQ input files are SOLiD data
    
    ///////////////////////////////////
    // Command line argtable settings
    ///////////////////////////////////
  
    struct arg_rem *sep1 = arg_rem(NULL, "\n==== Essential options (see also \"Input files\" below) ====\n");

    struct arg_lit *h = arg_lit0("h","help","Request help.");
 	struct arg_lit *version = arg_lit0("v","version","Print the build version and exit.");
    struct arg_lit *v = arg_lit0(NULL,"verbose","Print detailed status to stderr during run. [FALSE]");
    struct arg_str *detail = arg_strn("d","detail","<seq_id>",0,100,"Only produce per-base statistics and reconstructed seqeunce for listed seq ids. [ALL SEQS]");
    struct arg_file *detail_file = arg_file0(NULL,"detail_file","<detail_file>","Filename of file containing sequence ids to treat as in --detail [NULL]");
    struct arg_lit *seq_recon = arg_lit0("q","recon_seq","Reconstruct sequences from input reads [NULL]");
    struct arg_file *seq_ref = arg_file0(NULL,"ref","<ref_file>","Filename of file containing reference sequences which will be output in place of reconstructed sequence [NULL]");
    struct arg_lit *mp_analysis = arg_lit0("m","mate_pair","Perform mate-pair analysis [NULL]");
    struct arg_int *split = arg_int0(NULL,"split","<n>","Split reference sequences into n-base pieces for analysis. [no splitting]");
    struct arg_lit *sep = arg_lit0(NULL,"separate_strands","Top and bottom stands of reference sequences are considered separately [FALSE]");
    struct arg_rem *sep_r = arg_rem(NULL, "=== NOTE: --separate_strands is exclusive of mate-pairing (-m), sequence reconstruction (-q) and splitting (--split)");
    struct arg_lit *rollup = arg_lit0(NULL,"rollup","Output reference taxa level summary stats [FALSE]");
    struct arg_int *seed = arg_int0(NULL,"seed","<n>","Seed used by random number generator [12345]");
    struct arg_int *num_threads = arg_int0(NULL,"num_threads","<n>","Number of threads to use on multicore systems. [NUM_PROCESSORS]");
    
    struct arg_rem *sep2 = arg_rem(NULL, "\n==== Options for bitscore calculation and reference selection ====\n");

    struct arg_dbl *t = arg_dbl0("t","bit_thresh","<n>","Minimum bitscore value for a ref sequence to be selected [0.0]");
    struct arg_dbl *frac = arg_dbl0("f","bit_fraction","<n>","Minimum bitscore value, as a fraction of the top scoring sequence [0.0]");
    struct arg_rem *t_frac = arg_rem(NULL, "=== NOTE: if both -f and -t are used, the highest resulting threshold will be used.");
    struct arg_int *l = arg_int0("l","read_map_limit","<n>","Maximum number of ref mappings for a read to be considered [50000]");
    struct arg_lit *s = arg_lit0("s","second_chance_reads","Allow reads to be mapped to edit dist +1 mappings when their best ref is eliminated. [FALSE]");
    struct arg_lit *r = arg_lit0(NULL,"relax_read_sharing","Allow reads to be mapped to all sequences, regardless of edit distance. [FALSE]");
    struct arg_lit *nt = arg_lit0(NULL,"all_taxa","Use count of all taxa in catalog to calc bitscores. [Count taxa from read mappings]");
    struct arg_lit *a = arg_lit0("a","absolute_bitscores","Use absolute bitscores to select ref seqs. Default is length normalized scores (bits/nt).");
    struct arg_lit *no_rand = arg_lit0(NULL,"no_rand","Assign duplicate reads within a reference to the first position seen. [randomly assign]");
    struct arg_dbl *sim_frac = arg_dbl0(NULL,"sim_frac","<n>","Fraction of aligned reads to use. For use with simulated readsets. [1.0]");
    struct arg_str *exclude = arg_strn("e","exclude","<seq_id>",0,100,"Exclude these sequences from all analysis [NULL]");
    struct arg_file *exclude_file = arg_file0(NULL,"exclude_file","<exclude_file>","Filename of file containing sequence ids to exclude from analysis [NULL]");
    struct arg_lit *invert_ex = arg_lit0(NULL,"invert_exclude","--exclude sequences are inverted: exclude all sequences except those in the exclude file. [NULL]");
    
    struct arg_rem *sep3 = arg_rem(NULL, "\n==== Options for coverage calculations ====\n");

    struct arg_lit *detailed = arg_lit0(NULL,"per_base","Per base statistics will be generated and output for detail sequences [summaries only]");
    struct arg_int *w = arg_int0("w","cov_window","<n>","Assumed read length, used for per position coverage maps [actual read len when available or 49]");
    struct arg_lit *detect_dups = arg_lit0(NULL,"detect_dups","Detect likely collapsed duplications. Assumes uniform coverage (e.g. single genome). [FALSE]");

    struct arg_rem *sep4 = arg_rem(NULL, "\n==== Options for sequence reconstruction and read filtering ====\n");
    
    struct arg_dbl *ambig_tol = arg_dbl0(NULL,"ambig_tol","<n>","Tolerance for permitting ambiguity codes in reconstruction 0.0-1.0 (1 = no ambig) [0.2]");
    struct arg_str *read_out = arg_str0(NULL,"read_output","<prefix>","Enable matching read output, to files with this prefix. [NULL]");
    struct arg_lit *read_out_nomap = arg_lit0(NULL,"output_nomatch","Modifies --read_output to write non-matching reads. [NULL]");
    struct arg_lit *read_out_gz = arg_lit0(NULL,"read_output_gzip","Modifies --read_output to write gzipped output files [FALSE]");

    struct arg_rem *sep5 = arg_rem(NULL, "\n==== Options for generating mate-pairing statistics ====\n");
    
    struct arg_dbl *mp_share_lim = arg_dbl0(NULL,"mp_share_lim","<n>","Minimum bitscore for shared sequence reporting [500.0]");
    struct arg_dbl *mp_mate_lim = arg_dbl0(NULL,"mp_mate_lim","<n>","Minimum bitscore for sequence mate-pair link reporting [250.0]");
    struct arg_lit *mp_strict = arg_lit0(NULL,"mp_strict","No mate-pair links between sequences for pairs that both map within a single sequence. [FALSE]");
    struct arg_lit *mp_inserts = arg_lit0(NULL,"mp_inserts","Perform insert size estimation using all sequences. [FALSE]");
    struct arg_lit *mp_circ = arg_lit0(NULL,"mp_circular","Allow circular self-linking mate-pairs joining ends of a single sequence. [FALSE]");
    struct arg_int *mp_cutoff = arg_int0(NULL,"mp_cutoff","<n>","Per pair insert size cutoff for inclusion in per base statistics. [0 (no filter)]");
    
    struct arg_rem *sep6 = arg_rem(NULL, "\n==== Input files ====\n");
    
    struct arg_file *read_files = arg_filen("r","read_file","<read_file>",0,100,"Input filename(s) for FASTQ read files [NULL]");
    struct arg_file *rev_read_files = arg_filen(NULL,"rev_read_file","<read_file>",0,100,"Input filename(s) for FASTQ read files for opposite strand mates [NULL]");
    struct arg_lit *f_bwa_samse = arg_lit0(NULL,"old_bwa_samse","Use the old BWA-style \"samse\" alignment format instead of SAM. [FALSE]");
    struct arg_file *c = arg_file0("c","catalog","<catalog_file>","Filename of reference catlog file [Required for --old_bwa_samse, optional otherwise]");
    struct arg_file *fr = arg_filen(NULL,"rev_align","<sam_file>",0,500,"Input files from opposite strand mate alignments (e.g. bwa samse -n)");
    struct arg_file *f = arg_filen(NULL,NULL,"<sam_file>",1,500,"SAM format input files from alignments (e.g. from: bwa samse -n)");
    struct arg_end *end = arg_end(20);

    void *argtable[] = {sep1, h, version, v, detail, detail_file, seq_recon, seq_ref, mp_analysis, split, sep, sep_r, rollup, seed, num_threads, sep2, t, frac, t_frac, l, s, r,
        nt, a, no_rand, sim_frac, exclude, exclude_file, invert_ex, sep3, detailed, w, detect_dups, sep4, ambig_tol, read_out, read_out_nomap, read_out_gz, sep5, mp_share_lim, mp_mate_lim, mp_strict,
        mp_inserts, mp_circ, mp_cutoff, sep6, read_files, rev_read_files, f_bwa_samse, c, fr, f, end};
    
    int arg_errors = 0;
    
    // All done with variable declarations!!
    
    ////////////////////////////////////////////////////////////////////////
    // Handle command line processing (via argtable2 library) 
    ////////////////////////////////////////////////////////////////////////
    
    arg_errors = arg_parse(argc, argv, argtable);
    
	if (version->count) {
		fprintf(stderr, "%s version: %s\n", argv[0], SS_BUILD_VERSION);
		exit(EXIT_SUCCESS);
    }
	
    if (h->count) {
        fprintf(stderr,"\nref_select is a utility for calculating a wide variety of useful\n");
        fprintf(stderr,"information from alignments of next-generation sequence reads to\n");
        fprintf(stderr,"large reference sequence datasets.\n\n"); 
        
        fprintf(stderr,"For example, it can be used to build a meta-genomic assembly graph\n");
        fprintf(stderr,"from paired-reads aligned back to de novo assembled contigs. In this\n");
        fprintf(stderr,"case, for SOLiD colorspace contigs, it will also reconstruct the\n");
        fprintf(stderr,"nucleotide versions of the contigs from the aligned reads.\n\n"); 
        
        fprintf(stderr,"It can also be used to precisely select (and accurately estimate the\n");
        fprintf(stderr,"relative abundances of) sequences from an alignment with a reference\n");
        fprintf(stderr,"library (e.g. metagenomic reads to 16S rDNA database, or mRNA\n"); 
        fprintf(stderr,"transcript reads to a database of genes from one or more genomes.\n\n");
      
        fprintf(stderr,"ref_select processes SAM format alignment files into a JSON format\n");
        fprintf(stderr,"\"sequence graph\" file ready to be processed by the graph_ops tool.\n");
        fprintf(stderr,"Optionally, ref_select will also read FASTQ files and filter out reads\n");
        fprintf(stderr,"(and mates) which are or aren't present in the provided SAM alignment.\n\n");
        
        fprintf(stderr,"The JSON output file is structured as a graph of \"nodes\" (selected\n");
        fprintf(stderr,"sequences) and \"edges\" (connections between nodes based on read-pairing)\n");
        fprintf(stderr,"A variety of useful information for the nodes and edges is included in\n");
        fprintf(stderr,"the output JSON file, which can be accessed using the graph_ops tool or\n");
        fprintf(stderr,"custom code written in any language with a JSON library.\n");
        
        fprintf(stderr, "\nUsage: ref_select [options] <sam_file1> [<sam_file2>]... > <outfile.json>\n");
        fprintf(stderr, "\n[options] : ");
        
        arg_print_syntaxv(stderr, argtable, "\n");
        arg_print_glossary(stderr, argtable, "%-30s %s\n");
        fprintf(stderr,"\n"); 
        exit(EXIT_SUCCESS);
    }
    
    if (arg_errors) { 
        fprintf(stderr, "ERROR: Invalid arguments:\n");
        arg_print_errors(stderr, end, "ref_select");
        arg_print_syntaxv(stderr, argtable, "\n");
        exit(EXIT_FAILURE);
    }
    
    if (t->count) {
        bit_thresh = t->dval[0];
    } else {
        bit_thresh = 0.0;
    }
    
    if (frac->count) {
        bit_frac = frac->dval[0];
    } else {
        bit_frac = 0.0;
    }
    
    if (l->count) {
        read_cnt_limit = l->ival[0];
    } else {
        read_cnt_limit = 50000;
    }
    
    if (split->count) {
        split_len = split->ival[0];
    } else {
        split_len = 0;
    }
    
    if (sep->count && mp_analysis->count) {       
        fprintf(stderr, "ERROR: Separate strand analysis (--separate_strands) and mate-pair analysis (-m) are mutually exclusive.\n");
        exit(EXIT_FAILURE);    
    }

    if (sep->count && seq_recon->count) {
        fprintf(stderr, "ERROR: Separate strand analysis (--separate_strands) and sequence reconstruction (-q) are mutually exclusive.\n");
        exit(EXIT_FAILURE);    
    }    
    
    if (sep->count && split->count) {
        fprintf(stderr, "ERROR: Separate strand analysis (--separate_strands) and reference splitting (--split) are mutually exclusive.\n");
        exit(EXIT_FAILURE);    
    }  
    
    if (sim_frac->count) {
        if ((sim_frac->dval[0] >= 1.0) || (sim_frac->dval[0] <= 0.0)) {
            fprintf(stderr, "ERROR: Simulated random fraction of input reads (--sim_frac) must be between 0 and 1 exclusive.\n");
            exit(EXIT_FAILURE);            
        }
    }
        
    if (ambig_tol->count) {       
        if (!seq_recon->count) {
            fprintf(stderr, "ERROR: Ambiguity code tolerance (-ambig_tol) requires sequence reconstruction (with -q)\n");
            exit(EXIT_FAILURE);            
        } else {
            // Update scoring matrix
            for (int x = 1; x < 16; x++) {
                weight[1][x] *= (0.5 + (2.5 * ambig_tol->dval[0]));
                weight[1][x] = (weight[1][x] > 1.0) ? 1.0 : weight[1][x];
            }
        }
    }
    
    if (read_out->count) {
        if (!read_files->count) {
            fprintf(stderr, "ERROR: Read output (--read_out) requested but no input reads provided (with -r)\n");
            exit(EXIT_FAILURE);            
        }
        if (strlen(read_out->sval[0]) == 0) {
            fprintf(stderr, "ERROR: Read output (--read_out). Zero length prefix is not permitted.\n");
            exit(EXIT_FAILURE);            
        }
    }
    
    if (read_out_nomap->count && !read_out->count) {
        fprintf(stderr, "ERROR: Non-matching Read output (--read_out_nomatch) requested but output prefix not specified (with --read_output)\n");
        exit(EXIT_FAILURE);            
    }
    
    if (read_out_gz->count && !read_out->count) {
        fprintf(stderr, "ERROR: Gzipped Read output (--read_output_gzip) requested but output prefix not specified (with --read_output)\n");
        exit(EXIT_FAILURE);            
    }
    
    if (seq_recon->count && seq_ref->count) {
        fprintf(stderr, "ERROR: Sequence reconstruction (-q) and sequence from external reference file (--ref) are mutually exclusive.\n");
        exit(EXIT_FAILURE);
    }
    
    if (seq_recon->count) {
        if (!read_files->count) {
            fprintf(stderr, "ERROR: Sequence reconstruction (-q) requested but no input reads provided (with -r)\n");       
            exit(EXIT_FAILURE);
        } else {
            // Verify all of the input files, check to see if any are colorspace
            colorspace_recon = 0;
            int nucleotide_recon = 0;
            char *file_name = NULL;
            
            for (int x = 0; x < read_files->count + rev_read_files->count; x++) {
                
                if (x < read_files->count) {
                    file_name = (char *) read_files->filename[x];
                } else {
                    file_name = (char *) rev_read_files->filename[x - read_files->count];
                }
                
                // Check if this is a valid fastq file
                if (!ss_is_fastq(file_name)) {
                    fprintf(stderr, "ERROR: %s is not a valid FASTQ file!\n", file_name);
                    exit(EXIT_FAILURE);
                }
                
                // Check if this file is a SOLiD colorspace encoded FASTQ
                if (ss_is_colorspace_fastq(file_name)) {
                    colorspace_recon = 1;
                } else {
                    nucleotide_recon = 1;
                }
            }
        }
    }
    
    if (s->count && r->count) {
        fprintf(stderr, "ERROR: Second chance reads (-s) and relaxed mapping (--relax_read_sharing) are mutually exclusive.\n");
        exit(EXIT_FAILURE);            
    } 
    
    if (w->count) {
        cov_window = w->ival[0];
    } else {
        cov_window = 49;
    }
    
    if (mp_share_lim->count) {       
        if (!mp_analysis->count) {
            fprintf(stderr, "ERROR: mp_share_lim may only be used when mate pair analysis (-m) is enabled.\n");
            exit(EXIT_FAILURE);            
        } else {
            if (mp_share_lim->dval[0] < 0.0) {
                fprintf(stderr, "ERROR: mp_share_lim must be >= 0.0.\n");
                exit(EXIT_FAILURE);            
            } else {
                share_lim = mp_share_lim->dval[0];
            }
        }
    }
    
    if (mp_mate_lim->count) {       
        if (!mp_analysis->count) {
            fprintf(stderr, "ERROR: mp_mate_lim may only be used when mate pair analysis (-m) is enabled.\n");
            exit(EXIT_FAILURE);            
        } else {
            if (mp_mate_lim->dval[0] < 0.0) {
                fprintf(stderr, "ERROR: mp_mate_lim must be >= 0.0.\n");
                exit(EXIT_FAILURE);            
            } else {
                mate_lim = mp_mate_lim->dval[0];
            }
        }
    }
    
    if ((!mp_inserts->count) && mp_circ->count) {
        fprintf(stderr, "ERROR: --mp_circular requires --mp_inserts.\n");
        exit(EXIT_FAILURE);            
    } 
    
    if (mp_cutoff->count) {
        if (!mp_inserts->count) {
            fprintf(stderr, "ERROR: --mp_cutoff requires --mp_inserts.\n");
            exit(EXIT_FAILURE);            
        } 

        if (mp_cutoff->ival[0] < 0) {
             fprintf(stderr, "ERROR: --mp_cutoff must be >= 0\n");
             exit(EXIT_FAILURE);  
        } else {
             mp_cutoff_len = mp_cutoff->ival[0];
        }
    } 
    
    if ((detect_dups->count) && !(mp_analysis->count)) {
        fprintf(stderr, "ERROR: --detect_dups requires --mp_analysis\n");
        exit(EXIT_FAILURE);
    }
    
    if (seq_ref->count) {
        gzFile ref_hnd;
        if (!(ref_hnd = gzopen(seq_ref->filename[0], "r"))) {
            int errorno = 0;
            const char *errstr = gzerror(ref_hnd, &errorno);
            fprintf(stderr, "ERROR: %s -- Reference file not found!\n Error %d, %s\n", seq_ref->filename[0], errorno, errstr);
            
            exit(EXIT_FAILURE);
        }
        gzclose(ref_hnd);
    }
    
    if (f_bwa_samse->count && !(c->count)) {
        fprintf(stderr, "ERROR: --old_bwa_samse requires use of a catalog file (-c)\n\n");
        exit(EXIT_FAILURE);
    }
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    // End argtable parameter validation
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    ///////////////////////////////////////////////////
    // Seed the pseudo-random number generator
    ///////////////////////////////////////////////////
    
    if (seed->count) {    // negative seed for ss_rand
        rnd_seed = (seed->ival[0] > 0) ? -seed->ival[0] : seed->ival[0];
    }
    
    rnd_state = ss_rseed(rnd_seed);
    
    ///////////////////////////////////////////////////
    // Set the appropriate number of threads
    // Note, 2x is to ensure at least 2 threads, and
    // to provide adequate opportunity for parallelism
    // this should be carefully tuned.
    
#ifdef _OPENMP
    if (num_threads->count) {
        if (num_threads->ival[0] < 2) {
            omp_set_num_threads(2);
            fprintf(stderr, "WARNING: --num_threads must be at least 2. Proceeding with 2 threads.\n");
        } else {
            omp_set_num_threads(num_threads->ival[0]);
        }
    } else {
        omp_set_num_threads(omp_get_num_procs());
    }
#endif
    
    ///////////////////////////////////////////////////
    // Read exclude file into hash table keyed by id
    ///////////////////////////////////////////////////

    // Handle explicit excludes on the command line
    
    if (exclude->count) {
        for (int i = 0; i < exclude->count; i++) {
            if (!(ex_ptr = calloc(1, sizeof(exclude_entry)))) {
                fprintf(stderr, "ERROR: calloc failed: exclude entry creation\n");
                exit(EXIT_FAILURE);
            }

            ss_strcat_utstring(&ex_ptr->seq_id, exclude->sval[i]);
            ss_HASH_ADD_UTSTR_STRUCT(exclude_table, seq_id, ex_ptr);
        }
    }

    // The exclude file is a file with the following single field (per line):
    // 1) Sequence ID (unique)
    
    if (exclude_file->count) {
        gzFile ex_hnd = NULL;
        // Try to open the exclude file
        if (!(ex_hnd = gzopen(exclude_file->filename[0], "r"))) {
            int errorno = 0;
            const char *errstr = gzerror(ex_hnd, &errorno);
            fprintf(stderr, "ERROR: %s -- Exclude file not found!\n Error %d, %s\n", exclude_file->filename[0], errorno, errstr);
            
            exit(EXIT_FAILURE);
        }
        
        while (ss_gzget_utstring(ex_hnd, str)) {
            
            if (!(ex_ptr = calloc(1, sizeof(exclude_entry)))) {  
                fprintf(stderr, "ERROR: calloc failed: exclude entry creation\n");
                exit(EXIT_FAILURE);        
            }
            
            ss_trunc_utstring(str, 1);
            utstring_concat(&ex_ptr->seq_id, str);
            ss_HASH_ADD_UTSTR_STRUCT(exclude_table, seq_id, ex_ptr);
            
        }
        
        // Close exclude file
        gzclose(ex_hnd);
    }
    
    // Process the detail output command line parms
    
    if (detail->count) {
        for (int i = 0; i < detail->count; i++) {
            if (!(det_ptr = calloc(1, sizeof(detail_entry)))) {
                fprintf(stderr, "ERROR: calloc failed: detail entry creation\n");
                exit(EXIT_FAILURE);
            }
            
            ss_strcat_utstring(&det_ptr->seq_name, detail->sval[i]);
            ss_HASH_ADD_UTSTR_STRUCT(detail_table, seq_name, det_ptr);
        }
    }
    
    if (detail_file->count) {
        
        gzFile sel_hnd = NULL;
        // Try to open the exclude file
        if (!(sel_hnd = gzopen(detail_file->filename[0], "r"))) {
            int errorno = 0;
            const char *errstr = gzerror(sel_hnd, &errorno);
            fprintf(stderr, "ERROR: %s -- Select file not found!\n Error %d, %s\n", detail_file->filename[0], errorno, errstr);
            
            exit(EXIT_FAILURE);
        }
        
        while (ss_gzget_utstring(sel_hnd, str)) {
            
            if (!(det_ptr = calloc(1, sizeof(detail_entry)))) {
                fprintf(stderr, "ERROR: calloc failed: detail entry creation\n");
                exit(EXIT_FAILURE);
            }
            
            ss_trunc_utstring(str, 1);
            
            utstring_concat(&det_ptr->seq_name, str);
            ss_HASH_ADD_UTSTR_STRUCT(detail_table, seq_name, det_ptr);
        }
        
        // Close exclude file
        gzclose(sel_hnd);
    }
    
    // If neither detail option is used, select all sequences by default
    if (!(detail->count || detail_file->count)) {
        det_all = 1;
    }
    
    ///////////////////////////////////////////////////
    // Read catalog file into hash table keyed by id
    ///////////////////////////////////////////////////
    
    // The catalog file is a tab separated file with the following fields:
    // 1) Sequence ID (unique)
    // 2) Sequence Length
    // 3) Reference Sequence Name (not necessarily unique)
    // 4) Sequence Description (used for output only)
    
    utstring_clear(str);
    
    // Pipe used for inter-thread communication
    int pipe1[2];

    if (c->count) {
        
        pipe(pipe1);
        
#pragma omp parallel sections default(shared)
        {
        
#pragma omp section
        {   // catalog file reader
            UT_string *catalog_file;
            utstring_new(catalog_file);
            ss_strcat_utstring(catalog_file, c->filename[0]);
            ss_stream_reader(catalog_file, pipe1[1]);
            utstring_free(catalog_file);
        }
        
#pragma omp section
        {   // catalog file datastruct builder
            
            FILE *R_in = fdopen(pipe1[0],"r");
            
            while (ss_get_utstring(R_in, str)) {

                char *save_ptr = NULL;
                char *seq_id_str = strtok_r(utstring_body(str), "\t", &save_ptr);        // Tokenize on tabs;
                char *len_str = strtok_r(NULL, "\t", &save_ptr);
                char *ref_seq_id_str = strtok_r(NULL, "\t", &save_ptr);
                char *desc_str = strtok_r(NULL, "\t", &save_ptr);
                
                if (desc_str) {  // Remove the '/n' from the end of the desc string.
                    desc_str[strlen(desc_str)-1] = '\0';
                }
                
                int seq_len = atoi(len_str);
                
                if (seq_len <= 0) {
                    fprintf(stderr, "ERROR: Catalog file: Sequence length must be greater than 0.\n");
                    exit(EXIT_FAILURE);
                }
                
                create_catalog_entry(seq_id_str, seq_len, ref_seq_id_str, desc_str, &exclude_table, &detail_table, &seq_table, &catalog, invert_ex->count, sep->count, det_all, split_len);
                
            }
            
            fclose(R_in);
            
        }  // End omp section
        }   // End omp parallel sections
        
    } else {  // Read the header of each SAM file and build the catalog automatically
        
        for (int x = 0; x < (f->count + fr->count); x++) {
            
            utstring_clear(file_name);
            
            // Figure out which parameter structure to get the next filename from
            if (x < f->count) {
                ss_strcat_utstring(file_name, f->filename[x]);
            } else {
                ss_strcat_utstring(file_name, fr->filename[x-f->count]);
            }
            
            pipe(pipe1);
            
#pragma omp parallel sections default(shared)
            {
            
#pragma omp section
            {   // alignment file reader
                unsigned long int cnt = 0;
                cnt = sam_header_stream_reader(file_name, pipe1[1]);
                
                if (!cnt) {
                    fprintf(stderr, "ERROR: Alignment file: %s is empty or does not exist.\n", utstring_body(file_name));
                    exit(EXIT_FAILURE);
                }
            }
            
#pragma omp section 
            {
                // Line begins with '@', that is, this is a header line. e.g.
                // @SQ     SN:NODE_78      LN:6232
                FILE *R_in = fdopen(pipe1[0],"r");
            
                int cnt = 0;
            
                while (ss_get_utstring(R_in, str)) {
                    
                    cnt++;
                    
                    char *input_buf = utstring_body(str);
                    
                    char *save_ptr = NULL;
                    
                    char *tok_ptr = strtok_r(input_buf, "\t", &save_ptr);        // Tokenize on tabs;
                    char *len_ptr = NULL;
                    char *seq_id_ptr = NULL;
                    char *desc_ptr = NULL;
                    
                    if (!strncmp(tok_ptr, "@SQ", 3)) {  // This is a SeQuence header line
                        while ((tok_ptr = strtok_r(NULL, "\t", &save_ptr))) {
                            
                            if (!seq_id_ptr && (!strncmp(tok_ptr, "SN:", 3))) {    // this is the sequence ID
                                seq_id_ptr = tok_ptr + 3;
                            } else if (!len_ptr && (!strncmp(tok_ptr, "LN:", 3))) { // this is the sequence length
                                len_ptr = tok_ptr + 3;
                            } else if (!len_ptr && (!strncmp(tok_ptr, "SP:", 3))) { // this is the optional species name (description)
                                desc_ptr = tok_ptr + 3;
                            }
                        }
                        
                        if (!desc_ptr) {
                            desc_ptr = seq_id_ptr;
                        }
                        
                        int seq_len = atoi(len_ptr);
                        
                        if (seq_len <= 0) {
                            fprintf(stderr, "ERROR: SAM file format error %s: Line %lu  %s\nSequence length (LN:) must be greater than 0.\n",utstring_body(file_name), cnt, input_buf);
                            exit(EXIT_FAILURE);
                        }
                        
                        if (!strlen(seq_id_ptr)) {
                            fprintf(stderr, "ERROR: SAM file format error %s: Line %lu  %s\nInvalid sequence name string (SN:)\n",utstring_body(file_name), cnt, input_buf);
                            exit(EXIT_FAILURE);
                        }
                        
                        HASH_FIND_STR(catalog, seq_id_ptr, entry);

                        if (entry) {
                            if ((entry->sequences) && (entry->sequences[0].seq_len != seq_len)) {
                                fprintf(stderr, "ERROR: SAM error %s: Line %lu  %s\nInvalid sequence length (LN:) %d, previous SAM file has inconsistent length: %d\n",utstring_body(file_name), cnt, seq_len, entry->sequences[0].seq_len);
                                exit(EXIT_FAILURE);
                            }
                        } else {
                            create_catalog_entry(seq_id_ptr, seq_len, seq_id_ptr, desc_ptr, &exclude_table, &detail_table, &seq_table, &catalog, invert_ex->count, sep->count, det_all, split_len);
                        }
                    }
                }
            }
            
            }
        }
    }
    
    // Set up vars containing the number of ref_seqs and seqs
    // in the catalog and output some stats.
    
    n_taxa = (unsigned int) HASH_COUNT(catalog);
    
    num_seq = (float) HASH_COUNT(seq_table);
    n_seq = HASH_COUNT(seq_table);
    
    if (v->count) {
        fprintf(stderr, "INFO: Num ref_seqs: %u\n", n_taxa);
        fprintf(stderr, "INFO: Num seqs: %.0f %u\n", num_seq, n_seq);
    }
    
    /////////////////////////////////////////////////////////////////
    // cat_array and seq_array are arrays of pointers to the elements
    // of the (now unchanging) hash tables of ref_sequences and
    // sequences. These are used in the code below to use OMP to
    // parallelize for-loops over these structures.
    
    // Alloc and populate an array of pointers to the catalog entries
    if (!(cat_array = calloc(n_taxa, sizeof(ref_seq_entry *)))) {
        fprintf(stderr, "ERROR: calloc failed: catalog array creation\n");
        exit(EXIT_FAILURE);
    }
    
    cat_ptr = cat_array;
    
    ss_HASH_FOREACH(catalog,entry) {
        *cat_ptr++ = entry;
    }
    
    // Alloc and populate an array of pointers to the sequence entires
    if (!(seq_array = calloc(n_seq, sizeof(seq_entry *)))) {
        fprintf(stderr, "ERROR: calloc failed: sequence array creation\n");
        exit(EXIT_FAILURE);
    }
    
    seq_ptr = seq_array;
    
    ss_HASH_FOREACH(seq_table,seq) {
        seq->seq_num = (seq_ptr - seq_array);
        *seq_ptr++ = seq;
    }
    
    /////////////////////////////////////////////////////////////////////
    // Now read in read mappings from the reference alignment file(s)
    /////////////////////////////////////////////////////////////////////
    
    int limit_skipped = 0;
    
    unsigned short int r_len = 0;
    
    int read_rev = 0; // Flag to flip the sign of the aligned position of a read
    
    for (int x = 0; x < (f->count + fr->count); x++) {
        
        utstring_clear(file_name);
        
        // Figure out which parameter structure to get the next filename from
        if (x < f->count) {
            ss_strcat_utstring(file_name, f->filename[x]);
            read_rev = 1;
        } else {
            ss_strcat_utstring(file_name, fr->filename[x-f->count]);
            read_rev = -1; // Flag to flip the sign of the aligned position of these reads
        }
        
        // This is a bit of a hack to get mate-pairs / paired-ends working
        // Because both mates in a pair  have the same read name, we add a prefix so
        // that the mates can be descriminated when reading the alignment files
        
        char read_order[3];  // read_type strings are prefixes for the read names
        char mate_order[2];  
        
        if (strstr(utstring_body(file_name), "read1")) {
            strcpy(read_order, "1~");
            strcpy(mate_order, "2");
        } else if (strstr(utstring_body(file_name), "read2")) {
            strcpy(read_order, "2~");
            strcpy(mate_order, "1");
        } else {
            strcpy(read_order, "S~");
            *mate_order = '\0';
        }
        
        utstring_clear(str);
        
        pipe(pipe1);
        
#pragma omp parallel sections default(shared)       
        {
            
#pragma omp section 
            {   // alignment file reader
                unsigned long int cnt = 0;
                if (f_bwa_samse->count) {
                    cnt = ss_stream_reader(file_name, pipe1[1]);
                } else {
                    cnt = sam_stream_reader(file_name, pipe1[1]);
                }
                
                if (!cnt) {
                    fprintf(stderr, "ERROR: Alignment file: %s is empty or does not exist.\n", utstring_body(file_name));
                    exit(EXIT_FAILURE);
                }
            }
            
#pragma omp section 
            {   // alignment file datastruct builder
                
                FILE *R_in = fdopen(pipe1[0],"r");      
                char *save_ptr = NULL;
                char file_read_type = '\0';
                
                while (ss_get_utstring(R_in, str)) {
                    
                    if (utstring_body(str)[0] == '>') {    // New read
                        
                        skip_flag = 0;
                        strtok_r(utstring_body(str), " ", &save_ptr);  // Tokenize on spaces
                        num = atoi(strtok_r(NULL, " ", &save_ptr));    // Tokenize on spaces
                        
                        // Detect colorspace read alignments by the naming convention that preserves the
                        // primer base and first color in the header line between '>' and '+'.
                        // e.g.  >TA+49|lambda:1231_1682_1381
                        //        ^^
                        // Determine if this alignment file contains colorspace or nucleotide space reads
                        if (!file_read_type) {
                            file_read_type = 'N';
                            // This checks syntax and protects against buffer overflow below.
                            if (strspn(utstring_body(str), ">+ACGTNacgtn") >= 3) {
                                if (utstring_body(str)[3] == '+') {
                                    file_read_type = 'C';
                                }
                            }
                        }
                        
                        if (num && (num <= read_cnt_limit)) {    // Check to see if this read is over read_cnt_limit
                            
                            // Look for the separator that indicates an embedded read_length
                            // in the read_identifier
                            tokptr = strchr(utstring_body(str), '|');
                            
                            // If found, then parse out the length and reposition the pointer 
                            // on the remainder of the identifier
                            if (tokptr) {

                                *tokptr++ = '\0';
                                
                                r_len = atoi(utstring_body(str)+((utstring_body(str)[3] == '+') ? 4 : 1));    // Skip the '>..+' or '>'

                                if (!r_len) {   // If the parse failed
                                    fprintf(stderr, "ERROR: Alignment file: read %s has an invalid readlength.\n", utstring_body(str));
                                    exit(EXIT_FAILURE);
                                }
                                
                                if (!w->count) {
                                    if (new_read_id == 1) {    
                                        cov_window = r_len;
                                    } else {
                                        cov_window = (r_len > cov_window) ? r_len : cov_window; 
                                    }
                                }
                                
                            } else {  // Else, use the default length
                                tokptr = utstring_body(str)+((utstring_body(str)[3] == '+') ? 4 : 1);            // Skip the '>..+' or '>'
                                r_len = cov_window;        // Default
                            }
                            
                            utstring_clear(name_str);

                            utstring_printf(name_str, "%s%s", read_order, tokptr);
                            
                            HASH_FIND_STR(read_names, utstring_body(name_str), read_name_ptr);
                            
                            if (read_name_ptr) {
                                
                                read_id = read_name_ptr->read_id;
                                
                            } else {
                                
                                if (sim_frac->count) {    // If we are skipping simulated reads
                                    
                                    HASH_FIND_STR(skipped, utstring_body(name_str), read_name_ptr);    // Check to see if this read has already been skipped
                                    
                                    if (read_name_ptr || (ss_rand(rnd_state) > sim_frac->dval[0])) {  // Simulate only a random fraction of reads being used
                                        skip_flag = 1;    // then don't use this read
                                        
                                        if (!read_name_ptr) {    // If this read hasn't already been seen
                                            sim_frac_skip_reads++;
                                            // Add this read to the skipped hash so that we can see later if it has already been rejected
                                            if (!(read_name_ptr = calloc(1, sizeof(read_name_id)))) {  
                                                fprintf(stderr, "ERROR: calloc failed: read skip entry creation\n");
                                                exit(EXIT_FAILURE);        
                                            }
                                            utstring_concat(&read_name_ptr->read_name, name_str);
                                            ss_HASH_ADD_UTSTR_STRUCT(skipped, read_name, read_name_ptr);
                                        }
                                    } else {
                                        read_id = new_read_id;
                                    }
                                } else {
                                    read_id = new_read_id;
                                }
                            }
                            
                        } else {
                            
                            // Begin skipping read mappings
                            // because there are more than read_cnt_limit of them
                            if (num > read_cnt_limit) {
                                limit_skipped++;
                            }
                            skip_flag = 1;
                        }
                        
                        r_ptr = NULL;
                        
                    } else if (!skip_flag) {
                        
                        s_id = strtok_r(utstring_body(str), "\t", &save_ptr);    // Tokenize on tabs
                        
                        // read_rev flips the sign of the position for reads 
                        // from input files indicated as "opposite strand" mates 
                        position = read_rev * atoi(strtok_r(NULL, "\t", &save_ptr));    
                        
                        ed = atoi(strtok_r(NULL, "\t", &save_ptr));
                        ed = (r->count) ? 0 : ed;  // Relaxed edit distance, all ed = 0
                        
                        // Look to see if this is an excluded sequence
                        HASH_FIND_STR(exclude_table, s_id, ex_ptr);
                        
                        // If not excluded, then process it!
                        if ((!invert_ex->count && !ex_ptr) || (invert_ex->count && ex_ptr)) {
                            
                            if (!(r_ptr) && read_id == new_read_id) {    // This is the first ref for this read
                                
                                // Reallocate the read_info array if it has grown too big
                                if (read_id >= read_info_len) {
                                    read_info_len += read_info_inc;
                                    if (!(read_table = realloc(read_table, read_info_len*sizeof(read_entry)))) {  
                                        fprintf(stderr, "ERROR: realloc failed: read_table array extension\n");
                                        exit(EXIT_FAILURE);        
                                    }
                                }
                                r_ptr = &read_table[read_id];
                                r_ptr->seq_count = 0;
                                r_ptr->read_len = r_len;
                                r_ptr->bit_info = 0.0;
                                r_ptr->read_chain = NULL;
                                r_ptr->edit_dist = ed;    
                                r_ptr->second_read_chain = NULL;
                                r_ptr->second_edit_dist = ed + 2;
                                r_ptr->mate_id = 0;
                                r_ptr->read_type = file_read_type;
                                
                                if (!(read_name_ptr = calloc(1, sizeof(read_name_id)))) {  
                                    fprintf(stderr, "ERROR: calloc failed: read name mapping creation\n");
                                    exit(EXIT_FAILURE);        
                                }
                                read_name_ptr->read_id = read_id;
                                r_ptr->read_name_ptr = read_name_ptr;
                                
                                utstring_concat(&read_name_ptr->read_name, name_str);
                                ss_HASH_ADD_UTSTR_STRUCT(read_names, read_name, read_name_ptr);
                                
                                // Set the mate pointers if found
                                
                                if (!r_ptr->mate_id && *mate_order)    { // Not singlet and no mate yet
                                    
                                    utstring_body(name_str)[0] = *mate_order;        // Change name to mate name
                                    
                                    HASH_FIND_STR(read_names, utstring_body(name_str), read_name_ptr);    // Try to find the mate
                                    if (read_name_ptr) {    // If found
                                        r_ptr->mate_id = read_name_ptr->read_id;    // Set the mate_ptr to the mate
                                        read_table[read_name_ptr->read_id].mate_id = read_id;  // Set the pointer coming back, too
                                    }
                                }
                                
                                new_read_id++;    // Committed to using this read_id in the array
                            }
                            
                            r_ptr = &read_table[read_id];
                            
                            if (ed < r_ptr->edit_dist) {    // If this is a better match than has been seen
                                
                                // We need to go back and undo a lot of stuff
                                
                                // If there has already been a read ed reshuffle
                                // There is no "third best" so all of these reads need to go!
                                
                                // remove this read from all of the refseqs it had mapped to
                                DL_FOREACH(r_ptr->read_chain,read) {
                                    HASH_DEL(read->seq_ptr->reads,read);
                                }
                                
                                if (s->count) { // If second chance reads enabled
                                    
                                    while (r_ptr->second_read_chain) {
                                        r2 = r_ptr->second_read_chain;
                                        DL_DELETE(r_ptr->second_read_chain,r_ptr->second_read_chain);
                                        free(r2); 
                                    }
                                    
                                    // If the edit distance isn't within two
                                    if (ed < r_ptr->edit_dist - 2) {
                                        // Kill all of these mappings, not good enough
                                        while (r_ptr->read_chain) {
                                            r2 = r_ptr->read_chain;
                                            DL_DELETE(r_ptr->read_chain,r_ptr->read_chain);
                                            free(r2);
                                        }
                                    } else {    // Otherwise, make this the new second string...
                                        r_ptr->second_read_chain = r_ptr->read_chain;
                                        r_ptr->second_edit_dist = r_ptr->edit_dist;
                                    }
                                } else {
                                    while (r_ptr->read_chain) {
                                        r2 = r_ptr->read_chain;
                                        DL_DELETE(r_ptr->read_chain,r_ptr->read_chain);
                                        free(r2);
                                    }
                                }
                                r_ptr->read_chain = NULL;    // Start fresh!
                                r_ptr->seq_count = 0;
                                r_ptr->edit_dist = ed;
                            }
                            
                            // Use position to determine the correct prefix to use
                            
                            utstring_clear(id_str);
                            
                            if (split_len) {
                                utstring_printf(id_str, "%s|%u", s_id, ((abs(position)-1) / split_len));
                            } else if (sep->count && (position < 0)) {
                                utstring_printf(id_str, "%s|R", s_id);
                            } else {
                                ss_strcat_utstring(id_str, s_id);
                            }
                            
                            HASH_FIND_STR(seq_table, utstring_body(id_str), seq);
                            
                            if (!seq && split_len) { 
                                
                                ss_strcat_utstring(id_str, "$");
                                HASH_FIND_STR(seq_table, utstring_body(id_str), seq);
                                
                                if (!seq) { 
                                    
                                    // Do over with chunk ID-1, since the last chunk includes a remainder
                                    utstring_clear(id_str);
                                    utstring_printf(id_str, "%s|%u$", s_id, ((abs(position)-1) / split_len)-1); 
                                    HASH_FIND_STR(seq_table, utstring_body(id_str), seq);
                                    
                                } 
                            }
                            
                            if (!seq) {    // Still not found!!
                                fprintf(stderr, "ERROR: Read (%s) mapped to seq: %s (at pos: %d) -- Sequence not found in catalog\n", utstring_body(&r_ptr->read_name_ptr->read_name)+2, utstring_body(id_str), position);
                                exit(EXIT_FAILURE);    
                            }
                            
                            if (ed == r_ptr->edit_dist) {    // If this is among the best matches
                                
                                HASH_FIND_INT(seq->reads, &read_id, read);
                                
                                // NOTE!  Read only added if it is not already found in this Reference Sequence
                                // Reads that map to more than one place in a reference sequence are not counted multiple times!
                                // Perhaps this should be an option!                            
                                
                                if (!read) {
                                    
                                    r_ptr->seq_count++;
                                    
                                    // Add this read_id to the read list for this ref
                                    if (!(read = calloc(1, sizeof(read_list_entry)))) {  
                                        fprintf(stderr, "ERROR: calloc failed: read list entry creation\n");
                                        exit(EXIT_FAILURE);        
                                    }
                                    read->r_id = read_id;
                                    
                                    read->seq_pos = (position < 0) ? (position + seq->offset) : (position - seq->offset);
                                    
                                    read->seq_ptr = seq;
                                    read->map_count = 1;
                                    HASH_ADD_INT(seq->reads, r_id, read);
                                    DL_PREPEND(r_ptr->read_chain, read);
                                    
                                } else {
                                    
                                    read->map_count++;    // Keep track of how many per reference
                                    // Decide whether to randomly replace the read mapping with this new one.
                                    if ((!no_rand->count) && (ss_rand(rnd_state) <= (1.0/(double)read->map_count))) {
                                        // Replace the read mapping with this new one
                                        
                                        read->seq_pos = (position < 0) ? (position + seq->offset) : (position - seq->offset);
                                        
                                        read->seq_ptr = seq;
                                    }                                
                                }
                                
                                //                            assert(abs(read->seq_pos) <= seq->seq_len);
                                
                            } else if (s->count && (ed <= r_ptr->second_edit_dist)) {    // If second chance reads enabled and this is among the second best matches
                                
                                HASH_FIND_INT(seq->reads, &read_id, read);
                                
                                if (!read) {
                                    
                                    // If the edit distance is better than +2, tighten up the criteria to +1
                                    if (ed <= r_ptr->second_edit_dist) {
                                        r_ptr->second_edit_dist = ed;
                                    }
                                    
                                    // Add this read_id to the second chance read list for this ref
                                    if (!(read = calloc(1, sizeof(read_list_entry)))) {  
                                        fprintf(stderr, "ERROR: calloc failed: read list entry creation\n");
                                        exit(EXIT_FAILURE);        
                                    }
                                    read->r_id = read_id;
                                    
                                    read->seq_pos = (position < 0) ? (position + seq->offset) : (position - seq->offset);
                                    
                                    read->seq_ptr = seq;
                                    read->map_count = 1;
                                    // Note, 2nd chance, NOT added to seq read lists
                                    // HASH_ADD_INT(seq->ref_seq->reads, r_id, read);
                                    DL_PREPEND(r_ptr->second_read_chain, read);
                                    
                                } else {
                                    
                                    read->map_count++;    // Keep track of how many per reference
                                    // Decide whether to randomly replace the read mapping with this new one.
                                    if ((!no_rand->count) && (ss_rand(rnd_state) <= (1.0/(double)read->map_count))) {
                                        // Replace the read mapping with this new one
                                        
                                        read->seq_pos = (position < 0) ? (position + seq->offset) : (position - seq->offset);
                                        
                                        read->seq_ptr = seq;
                                    }
                                }
                            }
                        }
                    }
                }
                fclose(R_in);
            }
        }
        
        num_reads = new_read_id-1;
    }
    
    if (v->count) {
         fprintf(stderr, "INFO: Total number of reads used: %d\n", num_reads);
    
         if (sim_frac->count) {  // Simulate only a random fraction of reads being used
             fprintf(stderr, "INFO: Reads skipped because of simulated fraction (--sim_frac): %d of %d (%f %f)\n", sim_frac_skip_reads, sim_frac_skip_reads+num_reads, num_reads/(double)(sim_frac_skip_reads+num_reads), sim_frac->dval[0]);
         }
    
         fprintf(stderr, "INFO: Reads skipped because mapped to more than %d ref seqs: %d\n", read_cnt_limit, limit_skipped);
    }  
    
    //////////////////////////////////////////////////////////////////
    // All done with file reads, time to start processing the data!
    //////////////////////////////////////////////////////////////////
    
    // Calculate the actual number of sequences in the mapped dataset
    unsigned int seq_tot = 0;
#pragma omp parallel for reduction(+ : seq_tot) schedule(static)
    for (int z = 0; z < n_seq; z++) {
        if (seq_array[z]->reads) {
            seq_tot++;
        } else {
            // Done with this ref seq if it has no reads.
            seq_array[z] = NULL;
        }
    }        
    
    if (v->count) {
         fprintf(stderr, "INFO: Actual number of ref sequences in readset: %d\n", seq_tot);
    }
    
    if (nt->count) {
        num_seq = (float) n_seq;
    } else {
        num_seq = (float) seq_tot;        
    }
    
    if (num_seq == 1.0) {
        fprintf(stderr, "WARNING: num_seq == 1, setting to 2\n");
        num_seq = 2.0;
    } else if (num_seq == 0.0) {
        fprintf(stderr, "ERROR: num_seq == 0 because %d reads aligned with references in the input.\n", num_reads);
        fprintf(stderr, "Check your alignment tool and reference sequences.\n");
        exit(EXIT_FAILURE);
    }
    
    if (v->count) {
         fprintf(stderr, "INFO: Number of taxa used for scoring: %f\n", num_seq);
    }
    
    float total_bits = 0.0;
    double l2 = log(2.0);
    
    // Calculate the read information scores from counts
#pragma omp parallel for schedule(static) reduction(+ : total_bits)
    for (int x=1; x <= read_id; x++) {
        if (read_table[x].seq_count) {
            read_table[x].bit_info = log(num_seq / (double) read_table[x].seq_count) / l2; 
        } else {
            read_table[x].bit_info = 0.0;
        }
        
        total_bits += read_table[x].bit_info;
    }
    
    if (v->count) {
         fprintf(stderr, "INFO: Total bits: %f\n", total_bits);    
    }
    
    // Calculate the sequence scores
    double max_val = 0.0;
    int max_seq = 0;
    double my_max_val = 0.0;
    int my_max_seq = -1;
    
#pragma omp parallel private(seq,read,my_max_seq,my_max_val) 
    {
        my_max_val = 0.0;
        my_max_seq = -1;
#pragma omp for schedule(static)
        for (int z = 0; z < n_seq; z++) {
            
            if ((seq = seq_array[z])) {
                
                seq->bit_score = 0.0;
                float bs = 0.0;
                
                ss_HASH_FOREACH(seq->reads,read) {
                    bs += read_table[read->r_id].bit_info;
                }
                
                if (a->count) { // if absolute bitscores in use
                    seq->bit_score = bs;
                } else {
                    seq->abs_bit_score = bs;
                    seq->bit_score = (double) seq->abs_bit_score / (double) seq->seq_len;
                }
                
                seq->org_bit_score = seq->bit_score;
                
                if (seq->bit_score > my_max_val) {
                    my_max_val = seq->bit_score;
                    my_max_seq = z;
                }
            }
        }
        
        if (my_max_val > max_val) {    // Avoid the critical section if possible
#pragma omp critical 
            if (my_max_val > max_val) {
                max_val = my_max_val;
                max_seq = my_max_seq;
            }
        }
    }
    
    if (v->count) {
         fprintf(stderr, "INFO: Finished calculating initial sequence scores: %s : %f\n", utstring_body(&seq_array[max_seq]->seq_id), max_val);
    }
    
    /////////////////////////////////////////////////////////////////////////////////////
    // Build mate-pair / common sequence stats when requested to do mp analysis

    if (mp_analysis->count) {    
      
        // fprintf(stderr, "seq_hash_entry: %lu, seq_list_entry: %lu\n",sizeof(seq_hash_entry), sizeof(seq_list_entry));
        
#pragma omp parallel default(none) private(seq,read,r2,r_ptr,r_ptr2,r3,shared_seq) shared(n_seq,seq_array,read_table,stderr,mp_strict)
        {
            seq_hash_entry *seq_hashes = NULL, *fwd_seq_hashes = NULL, *bwd_seq_hashes = NULL, *free_seq_hashes = NULL, *last_seq_hash = NULL;
            
#pragma omp for schedule(dynamic,500)
            for (int z = 0; z < n_seq; z++) {    // Walk through all sequences
                
                if ((seq = seq_array[z])) {    // IF this sequence has any reads
                    
                    seq->shared = seq->forward = seq->backward = NULL;
                    seq_hashes = fwd_seq_hashes = bwd_seq_hashes = NULL;
                    
                    ss_HASH_FOREACH(seq->reads,read) {            // Go through all of its reads
                        
                        r_ptr = &read_table[read->r_id];
                        
                        DL_FOREACH(r_ptr->read_chain, r2) {    // For each mapping of that read to some sequence
                            
                            if (r2->seq_ptr != seq) {        // If that sequence isn't *this* sequence
                                
                                // Look up to see if we've already seen a read shared between these two
                                HASH_FIND_INT(seq_hashes,&r2->seq_ptr->seq_num,shared_seq);    
                                
                                // If not, make a new shared read struct
                                if (!shared_seq) {
                                    
                                    if (!free_seq_hashes) {    // If there are no more structs free, make one!
                                        if (!(free_seq_hashes = malloc(sizeof(seq_hash_entry)))) {
                                            fprintf(stderr, "ERROR: malloc failed: seq_hash_entry list\n");
                                            exit(EXIT_FAILURE);    
                                        }
                                        free_seq_hashes->hh.next = NULL;
                                    }
                                    
                                    shared_seq = free_seq_hashes;
                                    free_seq_hashes = free_seq_hashes->hh.next;
                                    memset(shared_seq, 0, sizeof(seq_hash_entry));
                                    shared_seq->seq_num = r2->seq_ptr->seq_num;
                                    HASH_ADD_INT(seq_hashes,seq_num,shared_seq);
                                } 
                                
                                // Update the info shared in common
                                shared_seq->shared_bitscore += r_ptr->bit_info;
                                shared_seq->num_reads++;
                                // Position at center of read -- Negative if sequences are on opposite strands
                                shared_seq->pos_sum += (abs(read->seq_pos) + (r_ptr->read_len >> 1)) * ((((long int) read->seq_pos * (long int) r2->seq_pos) < 0) ? -1 : 1);
                            }
                        }
                        
                        // Now do mate pairs too
                        int mate = read_table[read->r_id].mate_id;
                        
                        if (mate) {
                            
                            // Go through all mappings for this mate and look to see if any of them are to the current sequence
                            unsigned int self = 0;
                            
                            r_ptr2 = &read_table[mate]; 
                            
                            // IF mp_strict is TRUE, then never use mate-pairs that both map within a single reference sequence
                            // (the current seq, or any other) to also link two separate reference sequences. This prevents 
                            // duplicated reference sequences from improperly/unnecessarily "cross-linking" to each other.
                            
                            if (mp_strict->count) {
                                
                                DL_FOREACH(r_ptr2->read_chain, r2) {    // For each mapping of the mate to some sequence
                                    if (r2->seq_ptr == seq) {        // If that sequence *is* this sequence
                                        self++;                        // Then this mate pair is self linked
                                        break;
                                    } else {
                                        // Alternatively:
                                        // If the sequence the mate maps to also has the original read, then it is self linked there
                                        DL_FOREACH(r_ptr->read_chain, r3) {    // For each mapping of that read to some sequence
                                            if (r3->seq_ptr == r2->seq_ptr) {        // If that sequence *is* this sequence
                                                self++;                        // Then this mate pair is self linked
                                                break;
                                            }    
                                        }
                                    }
                                }
                            }                            
                            
                            DL_FOREACH(r_ptr2->read_chain, r2) {    // For each mapping of that read to some sequence
                                
                                // Reads that map to the zero position don't count (probably an aligner bug)
                                // If mp_strict is TRUE, self will be non-zero for all reads with a mate that maps to
                                // some shared reference sequence. In that case, only proceed when the current sequence is 
                                // one of those self-linked seqeunces.
                                
                                if ((read->seq_pos && (!self || r2->seq_ptr == seq))) {
                                    
                                    // If this is a "forward" mate
                                    if (((read->seq_pos < 0) && (utstring_body(&r_ptr2->read_name_ptr->read_name)[0] == '1')) || 
                                        ((read->seq_pos > 0) && (utstring_body(&r_ptr2->read_name_ptr->read_name)[0] == '2'))) {
                                        
                                        // Look up to see if we've already seen a read shared between these two
                                        HASH_FIND_INT(fwd_seq_hashes,&r2->seq_ptr->seq_num,shared_seq);    
                                        
                                        // If not, make a new shared read struct
                                        if (!shared_seq) {
                                            
                                            if (!free_seq_hashes) {    // If there are no more structs free, make one!
                                                if (!(free_seq_hashes = malloc(sizeof(seq_hash_entry)))) {
                                                    fprintf(stderr, "ERROR: malloc failed: seq_hash_entry list\n");
                                                    exit(EXIT_FAILURE);    
                                                }
                                                free_seq_hashes->hh.next = NULL;
                                            }
                                            
                                            shared_seq = free_seq_hashes;
                                            free_seq_hashes = free_seq_hashes->hh.next;
                                            memset(shared_seq, 0, sizeof(seq_hash_entry));
                                            shared_seq->seq_num = r2->seq_ptr->seq_num;
                                            HASH_ADD_INT(fwd_seq_hashes,seq_num,shared_seq);
                                        } 
                                        
                                        // Update the info shared in common: minimum bitscore of the mapping of the two mates
                                        shared_seq->shared_bitscore += (r_ptr2->bit_info < read_table[read->r_id].bit_info) ? r_ptr2->bit_info : read_table[read->r_id].bit_info;
                                        shared_seq->num_reads++;
                                        // Position at center of read -- Negative if sequences are on opposite strands  NOTE! SOLiD Specific!  Same stand Same Orientation!
                                        shared_seq->pos_sum += (abs(read->seq_pos) + (r_ptr2->read_len >> 1)) * ((((long int) read->seq_pos * (long int) r2->seq_pos) < 0) ? -1 : 1);
                                        
                                    } else { // Else this is a "backward" mate
                                        
                                        // Look up to see if we've already seen a read shared between these two
                                        HASH_FIND_INT(bwd_seq_hashes,&r2->seq_ptr->seq_num,shared_seq);    
                                        
                                        // If not, make a new shared read struct
                                        if (!shared_seq) {
                                            
                                            if (!free_seq_hashes) {    // If there are no more structs free, make one!
                                                if (!(free_seq_hashes = malloc(sizeof(seq_hash_entry)))) {
                                                    fprintf(stderr, "ERROR: malloc failed: seq_hash_entry list\n");
                                                    exit(EXIT_FAILURE);    
                                                }
                                                free_seq_hashes->hh.next = NULL;
                                            }
                                            
                                            shared_seq = free_seq_hashes;
                                            free_seq_hashes = free_seq_hashes->hh.next;
                                            memset(shared_seq, 0, sizeof(seq_hash_entry));
                                            shared_seq->seq_num = r2->seq_ptr->seq_num;
                                            HASH_ADD_INT(bwd_seq_hashes,seq_num,shared_seq);
                                        } 
                                        
                                        // Update the info shared in common: minimum bitscore of the mapping of the two mates
                                        shared_seq->shared_bitscore += (r_ptr2->bit_info < read_table[read->r_id].bit_info) ? r_ptr2->bit_info : read_table[read->r_id].bit_info;;
                                        shared_seq->num_reads++;
                                        // Position at center of read -- Negative if sequences are on opposite strands  NOTE! SOLiD Specific!  Same stand Same Orientation!
                                        shared_seq->pos_sum += (abs(read->seq_pos) + (r_ptr2->read_len >> 1)) * ((((long int) read->seq_pos * (long int) r2->seq_pos) < 0) ? -1 : 1);
                                    }
                                }
                            }
                        }
                    }
                    
                    // Once all reads have been processed for this sequence, find the best alternative sequences
                    
                    if (seq_hashes) {
                        
                        int cnt = 0;
                        
                        if (!(seq->shared = calloc(HASH_COUNT(seq_hashes), sizeof(seq_list_entry)))) {
                            fprintf(stderr, "ERROR: calloc failed: seq_list_entry list\n");
                            exit(EXIT_FAILURE);                            
                        }
                        
                        // Make a sorted singly linked list of shared sequences
                        ss_HASH_FOREACH(seq_hashes, shared_seq) {
                            seq->shared[cnt].shared_bitscore = shared_seq->shared_bitscore; 
                            seq->shared[cnt].mean_pos = shared_seq->pos_sum / shared_seq->num_reads; 
                            seq->shared[cnt].seq_num = shared_seq->seq_num;
                            seq->shared[cnt].seq_ptr = seq_array[shared_seq->seq_num];
                            cnt++;
                            last_seq_hash = shared_seq;
                        }
                        seq->num_shared = cnt;
                        
                        qsort(seq->shared, cnt, sizeof(seq_list_entry), seq_entry_bitscore_cmp);
                        
                        // Add hash entries back to the free pool
                        last_seq_hash->hh.next = free_seq_hashes;
                        free_seq_hashes = seq_hashes;
                    } else {
                        seq->num_shared = 0;
                        seq->shared = NULL;
                    }
                    
                    // Again for forward mates
                    if (fwd_seq_hashes) {
                        
                        int cnt = 0;
                        
                        if (!(seq->forward = calloc(HASH_COUNT(fwd_seq_hashes), sizeof(seq_list_entry)))) {
                            fprintf(stderr, "ERROR: calloc failed: seq_list_entry list\n");
                            exit(EXIT_FAILURE);                            
                        }
                        
                        // Make a sorted singly linked list of shared sequences
                        ss_HASH_FOREACH(fwd_seq_hashes, shared_seq) {
                            seq->forward[cnt].shared_bitscore = shared_seq->shared_bitscore; 
                            seq->forward[cnt].mean_pos = shared_seq->pos_sum / shared_seq->num_reads; 
                            seq->forward[cnt].seq_num = shared_seq->seq_num;
                            seq->forward[cnt].seq_ptr = seq_array[shared_seq->seq_num];
                            cnt++;
                            last_seq_hash = shared_seq;
                        }
                        seq->num_forward = cnt;
                        
                        qsort(seq->forward, cnt, sizeof(seq_list_entry), seq_entry_bitscore_cmp);
                        
                        // Add hash entries back to the free pool
                        last_seq_hash->hh.next = free_seq_hashes;
                        free_seq_hashes = fwd_seq_hashes;
                    } else {
                        seq->num_forward = 0;
                        seq->forward = NULL;
                    }
                    
                    // And again for backward mates
                    if (bwd_seq_hashes) {
                        
                        int cnt = 0;
                        
                        if (!(seq->backward = calloc(HASH_COUNT(bwd_seq_hashes), sizeof(seq_list_entry)))) {
                            fprintf(stderr, "ERROR: calloc failed: seq_list_entry list\n");
                            exit(EXIT_FAILURE);                            
                        }
                        
                        // Make a sorted singly linked list of shared sequences
                        ss_HASH_FOREACH(bwd_seq_hashes, shared_seq) {
                            seq->backward[cnt].shared_bitscore = shared_seq->shared_bitscore; 
                            seq->backward[cnt].mean_pos = shared_seq->pos_sum / shared_seq->num_reads; 
                            seq->backward[cnt].seq_num = shared_seq->seq_num;
                            seq->backward[cnt].seq_ptr = seq_array[shared_seq->seq_num];
                            cnt++;
                            last_seq_hash = shared_seq;
                        }
                        seq->num_backward = cnt;
                        
                        qsort(seq->backward, cnt, sizeof(seq_list_entry), seq_entry_bitscore_cmp);
                        
                        // Add hash entries back to the free pool
                        last_seq_hash->hh.next = free_seq_hashes;
                        free_seq_hashes = bwd_seq_hashes;
                    } else {
                        seq->num_backward = 0;
                        seq->backward = NULL;
                    }
                }
            }
            
            // Free all of the used hash structs
            while (free_seq_hashes) {
                last_seq_hash = free_seq_hashes->hh.next;
                free(free_seq_hashes);    // Free up the pool of hashes used by this thread
                free_seq_hashes = last_seq_hash;
            }
        }    
    }
    
    if (v->count) {    
        fprintf(stderr, "INFO: Pre-selection\n");
    }
    
    //////////////////////////////////////////////////////////////////////////////
    // Capture the maximum sequence score and recalculate the bit threshold
    // if the -f (bit_fraction) parameter is in use.
    
    if (frac->count) {
        double temp;
        
        temp = seq_array[max_seq]->bit_score * bit_frac;
        
        bit_thresh = (temp > bit_thresh) ? temp : bit_thresh;
        
        if (v->count) {    
            fprintf(stderr, "INFO: Ref sequence selection threshold: %f\n", bit_thresh);
        }
    }
    
    // Now, remove reads from taxa (and remove taxa with no reads)
    unsigned int loop_count = 1;
    
    do {
        seq_entry *prev_seq = seq_array[max_seq];
        int prev = max_seq;
        
        if (v->count) {
            fprintf(stderr, "INFO: Selected: %u  %s : %f : %d\n", loop_count, utstring_body(&seq_array[max_seq]->seq_id), max_val, max_seq);
        }
        
        loop_count++;
        
        max_val = -1.0e50;
        max_seq = -1;
        
        // Start with the reads of the previous top scoring sequence
        // Walk through them, and remove all references to those reads in other ref seqs
        // Adjusting the reference bitscores to reflect each removed read
        
        if (a->count) { // If using absolute bitscores
            
            ss_HASH_FOREACH(prev_seq->reads, read) {
                
                r_ptr = &read_table[read->r_id];
                
                DL_FOREACH(r_ptr->read_chain, r2) {
                    if (r2 != read) {
                        seq = r2->seq_ptr;
                        HASH_DEL(seq->reads, r2);
                        seq->bit_score -= r_ptr->bit_info;
                    }
                }
            }
            
        } else { // !a->count = Length normalized bitscores
            
            ss_HASH_FOREACH(prev_seq->reads, read) {
                
                r_ptr = &read_table[read->r_id];
                
                DL_FOREACH(r_ptr->read_chain, r2) {
                    if (r2 != read) {
                        seq = r2->seq_ptr;
                        HASH_DEL(seq->reads, r2);
                        
                        seq->abs_bit_score -= r_ptr->bit_info;
                        seq->bit_score = (double) seq->abs_bit_score / (double) seq->seq_len;
                    }
                }
            }
        }
        
        HASH_DEL(seq_table, prev_seq);
        ss_HASH_ADD_UTSTR_STRUCT(selected, seq_id, prev_seq);
        seq_array[prev] = NULL;            
        
        // Now, calculate the new max score
        
#pragma omp parallel private(seq,my_max_seq,my_max_val) 
        {
            my_max_val = -1.0e50;
            my_max_seq = -1;
            
#pragma omp for schedule(static)
            for (int z = 0; z < n_seq; z++) {
                if ((seq = seq_array[z])) {
                    if (!seq->reads) {
                        // Blow away entries with no reads
                        seq_array[z] = NULL;
                    } else if (seq->bit_score > my_max_val) {
                        my_max_val = seq->bit_score;
                        my_max_seq = z;
                    }
                }
            }
            
            if (my_max_val >= max_val) {    // Avoid the critical section if possible
#pragma omp critical 
                if ((my_max_val > max_val) || ((my_max_val > max_val) && (my_max_seq < max_seq))) {
                    max_val = my_max_val;
                    max_seq = my_max_seq;
                }
            }
        }
        
    } while ((max_seq >= 0) && (max_val >= bit_thresh));
    
    if (v->count) {    
        fprintf(stderr, "INFO: Alloc Coverage / Reconstruction buffers\n");
    }
    
    // create an array of selected sequence pointers so OpenMP
    // can parallelize the reconstruction code across sequences
    
    int num_sel = 0 ;
    num_sel = HASH_COUNT(selected);
    seq_entry **sel_list = NULL, **sel_ptr = NULL;
    
    if (!(sel_list = calloc(num_sel, sizeof(seq_entry **)))) {
        fprintf(stderr, "ERROR: calloc failed: selection array creation\n");
        exit(EXIT_FAILURE);    
    }
    
    // memset(sel_list, 0, num_sel * sizeof(seq_entry **));
    sel_ptr = sel_list;
    
    // Populate the array of selected sequence pointers
    ss_HASH_FOREACH(selected,seq) {
        *sel_ptr++ = seq;
    }
    
    /////////////////////////////////////////////////////////////////////////////////////////
    // If a detailed coverage report and/or sequence reconstruction has been requested
    // Allocate per sequence, per position arrays to hold read/coverage data
    
    int recon_array_size = 0;
    
    if (detailed->count || seq_recon->count || mp_inserts->count || seq_ref->count) {
        
        /***********************************************************************************
         
         The seq->recon_array buffer is allocated for all sequences when det_all == 1 
         (the default), or for only those sequences present in the detail_table (populated
         with the --detail or --select_file command line options.
         
         This buffer is structured as a 2D array of floats:   
         
            seq->recon_array[M][N] 
         
         where N is the seq->seq_len+2 (to allow for summation and other operations on
             each row in the sequence recon dynamic programming algorithms.)
         
         and where M (recon_array_size) is a variable number of rows based on the following 
         cases:
         
         (M = 1)
         
         seq->recon_array[0][0..N] ---> Per-base read coverage (always present for selected seqs)

         (M += 1) -- When reference sequence is externally read in from a fasta file, (mutually exclusive of
         sequence reconstruction, below) that sequence is available in seq->recon_array[1] as:
         
         (char *) &(seq->recon_array[1][0]) ---> ~ 1/4 of the bytes of the float array cast to char
         with the remainder zero filled.
         
         (M += 15) -- These rows are allocated when --recon_seq is used:
         
         seq->recon_array[1..15][0..N] ---> Nucleotide dynamic programming tables
                                            indexed by: 2=A,3=C, (... for MGRSVTWYHKDBN)

         (M += 15) -- These rows are allocated when --recon_seq is used and there are colorspace reads:
         
         seq->recon_array[16..30][0..N] ---> Colorspace dynamic programming tables
         indexed by the colorspace equivalent of:
         
         NOTE! AFTER SEQUENCE RECONSTRUCTION the final sequence string is available in
         (char *) &(seq->recon_array[1][0]), but in this case, this over-writes the first
         row of the dynamic programming table, which is no longer needed at this point.

         (M += 2) -- These last two rows are allocated when the --mate_pair option is used
         
         seq->recon_array[M-2][0..N]  ---> Per-base mate-pair insert size
         seq->recon_array[M-1][0..N]  ---> Per-base physical matepair coverage
         
         
         ***********************************************************************************/

        recon_array_size = 1;
        
        if (mp_inserts->count) {
            // MP Insert stats are always in the last two arrays (recon_array_size - 1 and recon_array_size - 2)
            recon_array_size = 3;
        }

        if (seq_recon->count) {
            // If there are colorspace reads then these arrays need an extra
            // 15 elements to keep track of colorspace errors
            recon_array_size += (colorspace_recon ? 30 : 15);
        }
        
        if (seq_ref->count) {
            // Use sequence from an external fasta reference file rather than reconstruction
            recon_array_size += 1;
        }
        
        ss_HASH_FOREACH(selected,seq) {
            
            if (!det_all) {
                HASH_FIND_STR(detail_table, utstring_body(&seq->ref_seq->seq_name), det_ptr);                
            }
            
            // Only alloc for selected detail sequences (to save memory) 
            // unless det_all is true, then do everything.
            if (det_all || det_ptr) {
                
                // Keep track of the maximum selected sequence length
                if (seq->seq_len > max_sel_seq_len) {
                    max_sel_seq_len = seq->seq_len;
                }
                
                // calloc a coverage array with an element for each sequence position
                if (!((seq->recon_array = calloc(recon_array_size, sizeof(float **))) && 
                      (seq->recon_array[0] = calloc(recon_array_size * (seq->seq_len+2), sizeof(float))))) {
                    
                    fprintf(stderr, "ERROR: callocs failed: coverage array creation\n");
                    exit(EXIT_FAILURE);        
                }
                
                for (int y = 1; y < recon_array_size; y++) {
                    // Set the pointer to each 1/array_size of the array alloced above.
                    seq->recon_array[y] = seq->recon_array[0] + y * (seq->seq_len+2);
                }
            }
        }
    }    
    
    // Abundance estimation
    // For each read in the sorted table
    // Go through and weed out the read mappings that are to nonselected
    // ref sequences, or those that score below the cutoff
    
    if (v->count) {
        fprintf(stderr, "INFO: Done with ref select.  Now estimating abundance.\n");
    }
    
    int second_chance_reads = 0;
    
#pragma omp parallel for default(none) shared(read_table,selected,bit_thresh,num_reads,s) private(r_ptr,read,seq,r2) reduction(+ : second_chance_reads) schedule(static)
    for (int read_id = 1; read_id <= num_reads; read_id++) {
        r_ptr = &read_table[read_id];
        r_ptr->seq_count = 0;
        read = r_ptr->read_chain;
        
        // Go through each mapping for this read
        // This loop doesn't use a DL_FOREACH(r_ptr->read_chain,read) construct
        // because it potentially deletes the entry for read, which DL_FOREACH
        // doesn't handle well.
        do {
            HASH_FIND_STR(selected, utstring_body(&read->seq_ptr->seq_id), seq);
            
            if (seq && (seq->bit_score >= bit_thresh)) {
                r_ptr->seq_count++;
                read = read->next;
            } else {
                r2 = read;
                read = read->next;
                DL_DELETE(r_ptr->read_chain, r2);
            }
        } while (read);    
        
        // If this read has no mappings left
        if (s->count && !(r_ptr->read_chain) && (r_ptr->second_read_chain)) {
            
            // And there are best_edit_distance + 1 reads available, use them
            read = r_ptr->second_read_chain;
            
            // Go through each second best mapping for this read
            
            // This loop doesn't use a DL_FOREACH(r_ptr->second_read_chain,read) construct
            // because it potentially deletes the entry for read, which DL_FOREACH
            // doesn't handle well.
            
            do {
                HASH_FIND_STR(selected, utstring_body(&read->seq_ptr->seq_id), seq);
                if (seq && (seq->bit_score >= bit_thresh)) {
                    r_ptr->seq_count++;
                    read = read->next;
                } else {
                    r2 = read;
                    read = read->next;
                    DL_DELETE(r_ptr->second_read_chain, r2);
                }
            } while (read);    
            
            r_ptr->read_chain = r_ptr->second_read_chain;
            r_ptr->second_read_chain = NULL;
            r_ptr->edit_dist = r_ptr->second_edit_dist;
            if (r_ptr->seq_count) {
                second_chance_reads++;
            }
        }
    } // End for
    
    if (v->count) {
        fprintf(stderr, "INFO: Done reallocating read mappings to remaining taxa.\n");
    }
    
    ////////////////////////////////////////////////////////////////////////////////////
    // This code estimates the mate-pair insertion size for all reference sequences
    
    if (mp_inserts->count) {
        
        if (v->count) {    
            fprintf(stderr, "INFO: Mate-pair insert estimation\n");
        }
        
#pragma omp parallel default(none) private(seq,read,r2,r_ptr,det_ptr) shared(mp_cutoff_len,recon_array_size,det_all,v,num_sel,sel_list,read_table,stderr,mp_strict,mp_circ)
        {
            
#pragma omp for schedule(dynamic,50)
            for (int z = 0; z < num_sel; z++) {    // Walk through all sequences
                
                if ((seq = sel_list[z])) {    // IF this sequence has any reads
                    
                    if (seq->recon_array) {    // If this sequence is chosen for detailed analysis
                        
                        double sum_dist = 0.0;
                        int mate_cnt = 0;
                        int circ_cnt = 0;
                        double circ_dist = 0.0;
                        double mean = 0.0;
                        double sse = 0.0;
                        double stdev = 0.0;
                        
                        ss_HASH_FOREACH(seq->reads, read) { 
                            // Go through all of the "read1" reads for this reference
                            if ((read->map_count == 1) && (read_table[read->r_id].mate_id) && (utstring_body(&read_table[read->r_id].read_name_ptr->read_name)[0] == '1')) {
                                r_ptr = &read_table[read_table[read->r_id].mate_id];
                                
                                DL_FOREACH(r_ptr->read_chain, r2) {
                                    // If this read maps uniquely, and it maps to this sequence, and it is on the same strand
                                    if ((r2->map_count == 1) && (r2->seq_ptr == read->seq_ptr) && 
                                        ((read->seq_pos >= 0 && r2->seq_pos >= 0) || (read->seq_pos < 0 && r2->seq_pos < 0))) {
                                        
                                        int start = abs(read->seq_pos);
                                        int end = abs(r2->seq_pos);
                                        
                                        //assert((start == read->seq_pos) || (start == -read->seq_pos));
                                        //assert((end == r2->seq_pos) || (end == -r2->seq_pos));
                                        
                                        if (start > end) {    // Swap so end is always greater
                                            int tmp = end;
                                            end = start;
                                            start = tmp;
                                            r_ptr = &read_table[read->r_id];
                                        }
                                        
                                        // "Normal pairing" top or bottom strand
                                        if (read->seq_pos < r2->seq_pos) {
                                            
                                            // Add the length of whichever read is "end" so that the full length of that read is included
                                            end += ((read->seq_pos > 0) ? r_ptr->read_len : read_table[read->r_id].read_len) - 1;
                                            
                                            // Clip to sequence length
                                            end = (end > seq->seq_len) ? seq->seq_len : end; 
                                            
                                            if (r_ptr->read_type == 'C') {
                                                start--;
                                            }
                                            
                                            if ((1 + end - start <= mp_cutoff_len) || (!mp_cutoff_len)) {    // Optionally filter on valid insert sizes
          
                                                mate_cnt++;
                                                sum_dist += end - start;
                                                
                                                for (int x = start; x < end; x++) {
                                                    seq->recon_array[recon_array_size-1][x] += 1.0;
                                                    seq->recon_array[recon_array_size-2][x] += 1 + end - start;
                                                }
                                                
                                            }
                                            
                                        } else if (mp_circ->count) {  // Not "normal pairing", mates are reversed on strand, or connect off ends of circular sequence.
                                            
                                            start += ((read->seq_pos > 0) ? r_ptr->read_len : read_table[read->r_id].read_len) - 1;
                                            
                                            // Clip end to sequence length
                                            end = (end > seq->seq_len) ? seq->seq_len : end; 
                                            start = (start > seq->seq_len) ? seq->seq_len : start; 
                                            
                                            if (r_ptr->read_type == 'C') {
                                                end--;
                                            }
                                            
                                            if ((1 + start + (seq->seq_len - end) < mp_cutoff_len) || (!mp_cutoff_len)) {    // Optionally filter on valid insert sizes

                                                mate_cnt++;
                                                circ_cnt++;
                                                // Modulo seq_len arithmetic ahead!
                                                sum_dist += start + (seq->seq_len - end);
                                                circ_dist += start + (seq->seq_len - end);
                                                
                                                for (int x = 0; x < start; x++) {
                                                    seq->recon_array[recon_array_size-1][x] += 1.0;
                                                    seq->recon_array[recon_array_size-2][x] += 1 + start + (seq->seq_len - end);
                                                    
                                                } 
                                                
                                                for (int x = end; x < seq->seq_len; x++) {
                                                    seq->recon_array[recon_array_size-1][x] += 1.0;
                                                    seq->recon_array[recon_array_size-2][x] += 1 + start + (seq->seq_len - end);
                                                } 
                                                
                                            }
                                        }
                                        
                                        break;
                                        
                                    }
                                }
                            }
                            
                        }
                        
                        if (mate_cnt) {
                            
                            mean = sum_dist/mate_cnt;
                            
                            ss_HASH_FOREACH(seq->reads, read) { // Go through all of the read1 reads for this reference
                                
                                if ((read->map_count == 1) && (read_table[read->r_id].mate_id) && (utstring_body(&read_table[read->r_id].read_name_ptr->read_name)[0] == '1')) {
                                    r_ptr = &read_table[read_table[read->r_id].mate_id];
                                    DL_FOREACH(r_ptr->read_chain, r2) {
                                        if ((r2->map_count == 1) && (r2->seq_ptr == read->seq_ptr) && ((read->seq_pos * r2->seq_pos) > 0)) {
                                            int diff = 0;
                                            if ((read->seq_pos < r2->seq_pos) && 
                                                (((diff = (abs(read->seq_pos - r2->seq_pos) + 
                                                           ((read->seq_pos > 0) ? r_ptr->read_len : read_table[read->r_id].read_len) - 1)) < mp_cutoff_len) || 
                                                            (!mp_cutoff_len))) {
                                                    
                                                    sse += pow(mean - diff, 2);
                                                
                                                } else if ((mp_circ->count) && 
                                                           (((diff = (seq->seq_len - abs(read->seq_pos - r2->seq_pos) + 
                                                                     ((read->seq_pos > 0) ? r_ptr->read_len : read_table[read->r_id].read_len) - 1)) < mp_cutoff_len) || 
                                                                      (!mp_cutoff_len))) { // Modulo seq_len arithmetic
                                                               
                                                    sse += pow(mean - diff, 2);
                                                               
                                            }
                                            break;
                                        }
                                    }
                                    stdev = sqrt(sse/mate_cnt);
                                }
                            }
                        }
                        
#pragma omp critical                        
                        if (mate_cnt) {
                            seq->mp_insert_mean = mean;
                            seq->mp_insert_stdev = stdev;
                            seq->mp_good_pairs = mate_cnt;
                            
                            if (v->count) {
                                fprintf(stderr, "INFO: %s (len: %d, pairs: %d) estimated insert size: %.1f +/- %.1f\n", utstring_body(&seq->seq_id), seq->seq_len, mate_cnt, mean, stdev);
                                if (circ_cnt) {
                                    fprintf(stderr, "INFO: Circular (pairs: %d) estimated insert size: %.1f\n", circ_cnt, circ_dist/circ_cnt);    
                                }
                            }
                            
                        } else {
                            if (v->count) {    
                                fprintf(stderr, "INFO: %s (len: %d) No mate pairs!\n", utstring_body(&seq->seq_id), seq->seq_len);
                            }
                        }
                    }    
                }    
            }    
        }  
        if (v->count) {
             fprintf(stderr, "INFO: End Mate_pair analysis\n");
        }
        
    }
    
    ////////////////////////////////////////////////////////////////////////////////////
    // Load FASTQ read info if sequence reconstruction is requested
    // This must be done before read table is qsorted below...

    if (seq_recon->count || read_out->count) {

        if (v->count) {
            fprintf(stderr, "INFO: Start FASTQ processing\n");
        }
        
        for (int fx = 0; fx < read_files->count + rev_read_files->count; fx++) {
            
            utstring_clear(file_name);
            
            if (fx < read_files->count) {
                ss_strcat_utstring(file_name, read_files->filename[fx]);
                read_rev = 1; // Flag to flip the sign of the aligned position of a read
            } else {
                ss_strcat_utstring(file_name, rev_read_files->filename[fx - read_files->count]);
                read_rev = -1; // Flag to flip the sign of the aligned position of a read
            }        
            
            if (v->count) {
                fprintf(stderr, "INFO: Reading FASTQ file: %s ...\n", utstring_body(file_name));
            }
            
            // This is a bit of a hack to get mate-pairs / paired-ends working
            // Because both mates in a pair have the same read name, we add a prefix so
            // that the mates can be descriminated when reading the alignment files
            char read_order[3];  // read_order strings are prefixes for the read names
            
            if (strstr(utstring_body(file_name), "read1")) {
                strcpy(read_order, "1~");
            } else if (strstr(utstring_body(file_name), "read2")) {
                strcpy(read_order, "2~");
            } else {
                strcpy(read_order, "S~");
            }
            
            if (read_out->count) {
                // Build a new output filename, use prefix (including any path info in it...)
                
                utstring_clear(read_out_fn);
                
                ss_strcat_utstring(read_out_fn, read_out->sval[0]);
                
                // Read to the last '/' and use that as the output filename suffix
                char *slash_pos = strrchr(utstring_body(file_name),'/');
                
                if (slash_pos) {
                    slash_pos++;
                } else {
                    slash_pos = utstring_body(file_name);
                }
                
                ss_strcat_utstring(read_out_fn, slash_pos);
            }
            
            utstring_clear(str);
            
            // Pipe used for inter-thread communication
            int pipe2[2];
            pipe(pipe1);
            pipe(pipe2);
            
#pragma omp parallel sections default(shared)       
            {
                
#pragma omp section 
                {   // FASTQ file reader
                    ss_stream_reader(file_name, pipe1[1]);
                }
                
#pragma omp section 
                {   // FASTQ file writer
                    if (read_out->count) {
                        ss_stream_writer(read_out_fn, pipe2[0], read_out_gz->count);
                    } 
                }                
                
#pragma omp section 
                {   // FASTQ Read handler
                    
                    FILE *R_in = fdopen(pipe1[0], "r");
                    
                    FILE *R_out = NULL;
                    if (read_out->count) {
                        R_out = fdopen(pipe2[1], "w");       
                    }
                    
                    while (ss_get_utstring(R_in, str4)) {
                        
                        // Skip the length field, if present
                        
                        // Look for the separator that indicates an embedded read_length
                        // in the read_identifier
                        
                        tokptr = NULL;
                        tokptr = strchr(utstring_body(str4), '|');
                        
                        // If found, then parse out the length and reposition the pointer 
                        // on the remainder of the identifier
                        
                        utstring_clear(name_str);
                        ss_strcat_utstring(name_str, read_order);

                        if (tokptr) {
                            ss_strcat_utstring(name_str, tokptr+1);
                        } else if (utstring_len(str4) >= 4 && *(utstring_body(str4)+3) == '+') {
                            // Colorspace read, use the default length (skipping over "@..+")
                            ss_strcat_utstring(name_str, utstring_body(str4)+4);
                        } else {  // Else, only skip @
                            ss_strcat_utstring(name_str, utstring_body(str4)+1);
                        }
                        
                        if (utstring_body(name_str)[utstring_len(name_str) - 3] == '/') {
                            ss_trunc_utstring(name_str, 3);
                        }
                        
                        // Find the info for this read
                        HASH_FIND_STR(read_names, utstring_body(name_str), read_name_ptr);
                        
                        ss_get_utstring(R_in, str);     // Sequence line
                        ss_get_utstring(R_in, str2);    // Extra comment line
                        ss_get_utstring(R_in, str3);    // Quality line
                        
                        // If read output enabled and this read maps to some selected sequence
                        if (read_out->count) {
                            
                            int write_it = 0;
                            
                            // If this is a singlet
                            if (*read_order == 'S') {
                                // If this read maps
                                if (read_name_ptr && (read_table[read_name_ptr->read_id].read_chain)) {
                                    if (!read_out_nomap->count) {    // Output matching singlet reads
                                        write_it = 1;
                                    }
                                } else if (read_out_nomap->count) {    // Output non-matching singlet reads
                                    write_it = 1;
                                }
                            } else {  // Paired read
                                
                                if (read_name_ptr && (read_table[read_name_ptr->read_id].read_chain)) {
                                    if (!read_out_nomap->count) {    // Output matching singlet reads
                                        write_it = 1;
                                    }
                                } else {    // Not selected
                                    // Go look to see if the mate is selected
                                    utstring_body(name_str)[0] = (*read_order == '1') ? '2' : '1';
                                    // Find out if the mate mapped
                                    read_name_id *read_name_ptr2 = NULL;
                                    HASH_FIND_STR(read_names, utstring_body(name_str), read_name_ptr2);
                                    utstring_body(name_str)[0] = *read_order;  // Set it back!
                                    
                                    // If this read's mate mapped to selected sequence, then output this read
                                    if (read_name_ptr2 && (read_table[read_name_ptr2->read_id].read_chain)) {
                                        if (!read_out_nomap->count) {    // Output non-matching reads with matching mates
                                            write_it = 1;                                    }
                                    } else if (read_out_nomap->count) {  // Output non-matching reads with non-matching mates
                                        write_it = 1;
                                    }
                                }
                            }
                            
                            if (write_it) {
                                fputs(utstring_body(str4), R_out);    // Name line
                                fputs(utstring_body(str), R_out);     // Sequence line
                                fputs(utstring_body(str2), R_out);    // Extra comment line
                                fputs(utstring_body(str3), R_out);    // Quality line
                            }
                        }

                        if (read_name_ptr && seq_recon->count) {    // If found
                         
                            // The value of color flag is used to put nucleotide data in
                            // the correct table when there is colorspace data being used
                            
                            // When colorspace reads in use, all color data goes in upper table (+15)

                            // Detect this case and use it to provide colorspace reconstruction 

                            char primer_base = '\0';
                            char first_color = '\0';
                            char *fastq_header_ptr = utstring_body(str4);
                            
                            // Get pointer to read_struct
                            r_ptr = &read_table[read_name_ptr->read_id];

                            if (r_ptr->read_type == 'C') {
                                if (*(fastq_header_ptr+3) != '+') {    // Also check for '+' in "@NN+ ..."
                                                                       // If not present, then something is wrong
                                    fprintf(stderr, "ERROR: Colorspace read: %d %d %d %s %s is missing primer+first color.\n", r_ptr->read_len,
                                            seq->seq_len, read->seq_pos, utstring_body(&seq->seq_id),
                                            utstring_body(&r_ptr->read_name_ptr->read_name));
                                    exit(EXIT_FAILURE);
                                } else {
                                    primer_base = *(fastq_header_ptr+1);
                                    first_color = *(fastq_header_ptr+2);
                                }
                            }
                            
                            // Set the read length from the sequence length
                            int len = r_ptr->read_len = utstring_len(str)-1;
                            
                            // For each read-sequence mapping
                            DL_FOREACH(r_ptr->read_chain,read) {
                                
                                seq = read->seq_ptr;
                                
                                if (seq->recon_array) {
                                    
                                    int inc = ((read_rev * read->seq_pos) < 0) ? -1 : 1;    // Negative pos = - strand
                                    int pos = abs(read->seq_pos);
                                    int z = 0;
                                    
                                    // Sometimes reads map at/off the ends, don't use those...
                                    if ((pos >= 0) && (pos + len - 1 < seq->seq_len)) {
                                        
                                        // Working pointer into sequence string
                                        ptr = (inc == 1) ? utstring_body(str) : utstring_body(str) + len - 1;

                                        // If this read is colorspace...
                                        if (r_ptr->read_type == 'C') {
                                        
                                            for (int x = pos; x < pos + len; x++) {
                                           
                                                // Convert nucleotide/color to bitcode
                                                z = code[(int)*ptr];
                                                // Move to next position on sequence
                                                ptr += inc;
                                            
                                                // Fill out ambiguity error table for this nucleotide / color
                                                for (int y = 1; y < 16; y++) {
                                                    seq->recon_array[y+15][x] += weight[((z & y) != 0)][y];
                                                }
                                            }
                                        
                                            // This is the nt code calculated from the colorspace read prefix
                                            z = code[(int)c2s[(int)primer_base][(int)first_color]];
                                            
                                            // If minus strand, then compliment the nt and put it in the correct position
                                            if (inc < 0) {
                                                z = comp[z];
                                                pos = pos + len;    // This is correct, not -1!  see example below...
                                                
                                                //   12345    Coord
                                                //   A      Nt goes at aligned position
                                                // T30110    Forward cs read
                                                //   AACAA    Sequence
                                                //   01100T Rev cs Read
                                                //       T  <- This is the case above, nt goes at aligned pos + len
                                                //   12345  Coord
                                            }
                                            
                                            // Fill out ambiguity error table for this nucleotide
                                            for (int y = 1; y < 16; y++) {
                                                seq->recon_array[y][pos] += weight[((z & y) != 0)][y];
                                            }
                                            
                                        } else {  // This is a nucleotide read

                                            for (int x = pos+1; x <= pos + len; x++) {
                                                
                                                // Convert nucleotide/color to bitcode
                                                z = code[(int)*ptr];
                                                // Rev comp when necessary
                                                if (inc < 0) {
                                                    z = comp[z];
                                                }
                                                // Move to next position on sequence
                                                ptr += inc;
                                                
                                                // Fill out ambiguity error table for this nucleotide / color
                                                for (int y = 1; y < 16; y++) {
                                                    seq->recon_array[y][x] += weight[((z & y) != 0)][y];
                                                }
                                            }
                                        }
                                        
                                    } else {
                                        if (v->count) {
                                             fprintf(stderr, "WARNING: Skipped improperly aligning read: %d %d %d %s %s\n", r_ptr->read_len, seq->seq_len, read->seq_pos, utstring_body(&seq->seq_id), utstring_body(&r_ptr->read_name_ptr->read_name));
                                        }
                                    }
                                }        
                            }
                        }             
                    }
                    fclose(R_in);
                    
                    if (R_out) {
                        fclose(R_out);       
                    }
                    
                }
            }
        }
    }
    
    if (v->count) {
         fprintf(stderr, "INFO: Abundance Estimation \n");
    }
    
    ////////////////////////////////////////////////////////////////////////////////////
    // Abundance estimation code
    
    // Sort the read table on the remaining seq_count calced above.
    qsort(read_table+1, num_reads, sizeof(read_entry), read_entry_count_cmp);
    
    if (v->count) {
        fprintf(stderr, "INFO: Done Sorting the read table by information. Num Reads: %d\n", num_reads);
    }    
    
    // Now, estimate abundance...
    
    double total_abund = 0.0;
    int remaining_reads = 0;
    
    // This can't be parallelized because it calculates cumulatively on the 
    // readlist sorted in increasing order of number of refs matched
    for (int read_id = 1; read_id <= num_reads; read_id++) {
        
        double tmp = 0.0;
        r_ptr = &read_table[read_id];
        
        // Special case for unique reads
        if (r_ptr->seq_count == 1) {
            
            remaining_reads++;
            read = r_ptr->read_chain;
            seq = read->seq_ptr;
            
            double add = 1.0;
            tmp = add/((double)seq->seq_len);
            total_abund += tmp;
            seq->abundance += tmp;

            seq->coverage += tmp * r_ptr->read_len;
            
            int start = abs(read->seq_pos);
            
            // Do not account for reads that are mapped off sequence ends
            if ((seq->recon_array) && ((start > 0) && (start + r_ptr->read_len - 1 <= seq->seq_len))) {
                
                // start += (start == 0);  
                int end = start + r_ptr->read_len;
                // end = (end > seq->seq_len) ? seq->seq_len : end;
                
                if (r_ptr->read_type == 'C') {
                    if (read->seq_pos > 0) {
                        start--;
                    } else {
                        end++;
                    }
                }
                for (int x = start; x < end; x++) {
                    seq->recon_array[0][x] += add;
                } 
            }
        }    
        // Non unique reads    
        else if (r_ptr->seq_count > 1) {
            
            remaining_reads++;
            // Go through each mapping for this read
            unsigned int n = r_ptr->seq_count;
            double scale = 0.0;
            unsigned int map_scale = 0;
            // Go through the reads to calculate their total current
            // abundance estimates.  This will be a scaling factor.
            
            DL_FOREACH(r_ptr->read_chain,read) {
                seq = read->seq_ptr; 
                if (seq->abundance == 0.0) {
                    n--;
                    map_scale += read->map_count;    // This is used for scale when no abundances
                } else {
                    scale += seq->abundance*(double)read->map_count;
                }                
            }
            
            if (n >= 1) { // If at least one ref has positive abundance
                DL_FOREACH(r_ptr->read_chain,read) {
                    seq = read->seq_ptr;
                    if (seq->abundance > 0.0) {
                        double add = ((seq->abundance*(double)read->map_count)/scale);
                        tmp = add/((double)seq->seq_len);
                        total_abund += tmp;
                        seq->abundance += tmp;
                        
                        seq->coverage += tmp * r_ptr->read_len;
                        
                        int start = abs(read->seq_pos);
                        
                        // Do not account for reads that are mapped off sequence ends
                        if ((seq->recon_array) && ((start > 0) && (start + r_ptr->read_len - 1 <= seq->seq_len))) {
                            
                            // start += (start == 0);  
                            int end = start + r_ptr->read_len;
                            // end = (end > seq->seq_len) ? seq->seq_len : end;
                            
                            if (r_ptr->read_type == 'C') {
                                if (read->seq_pos > 0) {
                                    start--;
                                } else {
                                    end++;
                                }
                            }                            
                            for (int x = start; x < end; x++) {
                                seq->recon_array[0][x] += add;
                            } 
                        }
                    }                 
                }
            } else {    // If none of these refs have abundance, split the read by map count
                
                scale = (double) map_scale;
                
                DL_FOREACH(r_ptr->read_chain,read) {
                    seq = read->seq_ptr;
                    double add = ((double)read->map_count/scale);
                    tmp = add/((double)seq->seq_len);
                    total_abund += tmp;
                    seq->abundance += tmp;
                    
                    seq->coverage += tmp * r_ptr->read_len;
                    
                    int start = abs(read->seq_pos);
                    
                    // Do not account for reads that are mapped off sequence ends
                    if ((seq->recon_array) && ((start > 0) && (start + r_ptr->read_len - 1 <= seq->seq_len))) {
                        
                        // start += (start == 0);  
                        int end = start + r_ptr->read_len;
                        // end = (end > seq->seq_len) ? seq->seq_len : end;
                        
                        if (r_ptr->read_type == 'C') {
                            if (read->seq_pos > 0) {
                                start--;
                            } else {
                                end++;
                            }
                        }                        
                        for (int x = start; x < end; x++) {
                            seq->recon_array[0][x] += add;
                        } 
                    }
                }
            }
        }
    }
    
    if (v->count) {
        fprintf(stderr, "INFO: Done estimating abundance, time to output results.\n");
    }
    
    ////////////////////////////////////////////////////////////////////////////
    // Output reconstructed sequence if requested

    if (seq_recon->count) {    
      
        if (v->count) {
            fprintf(stderr, "INFO: Processing reconstruction arrays.\n");
        }
        
        max_sel_seq_len += 2;    // Extra 2 rows of working space.
        
#pragma omp parallel default(shared)
        {
            
            // Allocate working structures
            
            float *dynamic[16];    // Dynamic programming score table pointer array
            char *bases[16];    // Called nucleotide table pointer array
            
            // calloc a scoring array with 16 elements for each sequence position
            if (!(dynamic[0] = calloc(16 * max_sel_seq_len, sizeof(float)))) {
                fprintf(stderr, "ERROR: calloc failed: dynamic array creation\n");
                exit(EXIT_FAILURE);        
            }
            
            // calloc a base array with 16 elements for each sequence position
            if (!(bases[0] = calloc(16 * max_sel_seq_len, sizeof(char)))) {
                fprintf(stderr, "ERROR: calloc failed: base array creation\n");
                exit(EXIT_FAILURE);        
            }
            
            if (v->count) {
                fprintf(stderr, "INFO: Working arrays alloced...\n");
            }
            
            // Initialize each sequence coverage array
            for (int y = 0; y < 16; y++) {
                // Set the pointer to each 1/array_size of the array alloced above.
                dynamic[y] = dynamic[0] + y * max_sel_seq_len;
                bases[y] = bases[0] + y * max_sel_seq_len;
            }
            
#pragma omp for private(ptr, seq)
            for (int i = 0; i < num_sel; i++) {
                
                seq = sel_list[i];
                
                if (seq->recon_array) {
                    
                    float min; 
                    float score;
                    int selected_nt;
                    int nt;
                    
                    // Clear the working storage
                    memset(dynamic[0], 0, (16 * max_sel_seq_len * sizeof(float)));
                    memset(bases[0], 0, (16 * max_sel_seq_len * sizeof(char)));    
                    
                    // For each position
                    for (int x = seq->seq_len; x >= 1; x--) {
                        
                        min = 1.0e20;
                        // For each prior nucleotide possibility
                        for (int y = 1; y < 16; y++) {
                            
                            if (colorspace_recon) {
                            
                                min = 1.0e20;
                                
                                // For each nuc/color possibility
                                for (int z = 15; z > 0; z--) {
                                    
                                    // calc the next nucleotide for this nuc/color combo (j,k)
                                    nt = match[y][z];
                                    
                                    // calc the dynamic error sum of this Nuc + Color combo
                                    score = dynamic[nt][x + 1] + seq->recon_array[z+15][x];
                                    
                                    // Keep track of which next nt value minimizes the error
                                    if (score <= min) { // <= because lower z colors that are = in score have less ambiguity
                                        min = score;
                                        selected_nt = nt;
                                    }
                                }
                                
                                // Account for the nucleotide score at this position
                                min += seq->recon_array[y][x];
                                
                                // For each color possibility
                                for (int z = 1; z < y; z++) {
                                    // Check to see if there is a less ambiguous call that is at least as good
                                    if (((z | y) == y) && (seq->recon_array[z+15][x] <= min)) {
                                        selected_nt = bases[z][x];    // Get out of ambiguity land, if possible
                                    }
                                }
                            } else {
                                // For each nucleotide possibility
                                if (y == 1) {  // Nothing in here depends on y, so just do it the first time.
                                    for (int nt = 15; nt > 0; nt--) {
                                        
                                        // calc the dynamic error sum of this Nuc + Nuc combo
                                        score = dynamic[nt][x + 1] + seq->recon_array[nt][x];
                                        
                                        // Keep track of which next nt value minimizes the error
                                        if (score <= min) { // <= because lower z colors that are = in score have less ambiguity
                                            min = score;
                                            selected_nt = nt;
                                        }
                                    }
                                }
                            }
                            
                            // Update the dynamic error sum for this position as the best color error + the base error of this nucleotide choice
                            dynamic[y][x] = min;
                            //  Remember the nucleotide choice for the next position
                            bases[y][x] = selected_nt;
                        }
                    }
                    
                    // Setup to begin outputting nucleotide sequence
                    min = dynamic[1][1];
                    int best = 1;
                    
                    int gc_sum = 0;
                    int gc_cnt = 0;
                    
                    // This table will be reused to write the sequence string
                    // Except the coverage table: recon_array[0][0]...
                    ptr = (char *) &(seq->recon_array[1][0]);
                    memset(ptr, 0, 2*seq->seq_len);
                    
                    // Find the initial nucleotide that minimizes the overall dynamic error sum
                    for (int z = 2; z < 16; z++) {
                        if (min > dynamic[z][1]) {
                            min = dynamic[z][1];
                            best = z;
                        }
                    }
                                      
                    if (seq->recon_array[0][1] != 0.0) {
                        *ptr++ = cc[best];
                        gc_cnt += (comp[best] != 0);
                        gc_sum += gc_tab[best];
                    } else {
                        *ptr++ = 'N';    
                    }
                    
                    // Loop through all of the bases
                    for (int x = 1; x < seq->seq_len; x++) {
                        
                        // The next base is looked up based on the current base and position
                        best = bases[best][x];
                        
                        // If not zero coverage, output the nucleotide base
                        if ((seq->recon_array[0][x] != 0.0) && (seq->recon_array[15+(colorspace_recon*15)][x] != 0.0)) {
                            *ptr++ = cc[best];
                            gc_cnt += (comp[best] != 0);
                            gc_sum += gc_tab[best];
                        } else {
                            *ptr++ = 'N';
                        }
                    }

                    *ptr++ = '\0';
                    
                    if (gc_cnt != 0) {
                        seq->gc_percent = (float) gc_sum / (float) gc_cnt;
                    } else {
                        seq->gc_percent = 0.0;
                    }
                }
            }
            
            // Free the dynamic structres
            free(dynamic[0]);
            free(bases[0]);
        }
    }    
    
    ////////////////////////////////////////////////////////////////////////////
    // Search for problems and dump all other results
    
    double total_fractional_abundance = 0.0;
    int num_selected = 0;
    
    // Calculate summary stats
    
    if (selected && selected->bit_score >= bit_thresh) {
        
#pragma omp parallel for default(none) private(seq) shared(num_sel,sel_list,total_abund,cov_window,mp_inserts,recon_array_size) reduction(+ : num_selected, total_fractional_abundance) 
        for (int z = 0; z < num_sel; z++) {    // Walk through all sequences
            
            if (!(seq = sel_list[z])) continue;
            
            total_fractional_abundance += seq->abundance/total_abund;
            num_selected++;
            
            // Only do all of this if this entry's sequences have recon_arrays.
            if (seq->recon_array) {
                double max = 0.0;
                double min = 1000000000.0;
                double total = 0.0;
                double mean = 0.0;
                double stddev = 0.0;
                int uncovered = 0;
                
                if (seq->seq_len > 2*cov_window) {
                    
                    for (int x = 1; x <= seq->seq_len; x++) {
                        double val = seq->recon_array[0][x];
                        
                        // Don't let the sequence ends skew the statistics
                        // Reads won't map at the ends...
                        if ((x >= cov_window) && (x < seq->seq_len-cov_window)) {
                            total += val;
                            max = (max < val) ? val : max;
                            min = (min > val) ? val : min;
                        }
                        uncovered += (val == 0);
                    }
                    
                    if (mp_inserts->count) {
                        
                        for (int x = 1; x <= seq->seq_len; x++) {
                            double val = seq->recon_array[recon_array_size-1][x];
                            
                            // Prenormalize the per position insert size
                            
                            seq->recon_array[recon_array_size-2][x] = val ? (seq->recon_array[recon_array_size-2][x] / val) : 0.0;
                            
                        }
                    }
                    
                    // Don't let the sequence ends skew the statistics
                    // Reads won't map at the ends...
                    
                    mean = total / (double) (seq->seq_len - 2*cov_window);
                    
                    for (int x = cov_window; x < seq->seq_len-cov_window; x++) {
                        seq->recon_array[0][x];
                        stddev += pow((seq->recon_array[0][x] - mean),2.0);
                    }
                    stddev = sqrt(stddev/(double)(seq->seq_len - 2*cov_window));
                    seq->uncovered = (double)uncovered / seq->seq_len;
                    seq->adj_coverage_mean = mean;
                    seq->adj_coverage_stddev = stddev;
                    seq->adj_coverage_min = min;
                    seq->adj_coverage_max = max;
                    seq->ref_seq->uncovered += uncovered;
                    
                }
            }
        }
    } else {
        fprintf(stderr, "\nWARNING: No taxa met the bit threshold limit!\n");
    }

    
    ////////////////////////////////////////////////////////////////////////////
    // Search for assembly problems within contig stats
    
    double mean_coverage = 0.0;
    float cov_threshold = 0.0;
    
    if (mp_analysis->count) {
        
        unsigned int total_length = 0;
        
        // Now, write all of the nodes up front
        if (selected && selected->bit_score >= bit_thresh) {
            
            // Calculate the mean coverage of all selected contigs
            if (detect_dups->count) {
             
#pragma omp parallel for default(none) private(seq) shared(num_sel,sel_list) reduction(+ : mean_coverage, total_length)
                for (int i = 0; i < num_sel; i++) {
                    seq = sel_list[i];
                    mean_coverage += seq->coverage * seq->seq_len;
                    total_length += seq->seq_len;
                }
                
                // Complete the average
                mean_coverage /= (double) total_length;
                cov_threshold = mean_coverage * 2.0;  // 200% of the mean coverage of all contigs

                if (v->count) {
                   fprintf(stderr, "\nINFO: Mean coverage (used for dup detection) = %.2f\nINFO: Duplication coverage threshold = %.2f\n", mean_coverage, cov_threshold);
                }
            }

#pragma omp parallel for default(none) private(seq) shared(cov_threshold,num_sel,sel_list,detect_dups,mean_coverage,stderr,mp_inserts,recon_array_size)
            for (int i = 0; i < num_sel; i++) {
                seq = sel_list[i];
                if (seq->recon_array) {
                    float *read_cov = seq->recon_array[0];
                    
                    // Look for collapsed duplicated regions
                    // This code assumes that "normal" coverage of contigs in this dataset is expected to
                    // be relatively constant.  That is, these contigs aren't from a metagenome or a mixed culture.
                    
                    if (detect_dups->count) {
                        
                        // Now find problem regions
                        float max_seen = 0.0;
                        float current_max = 0.0;
                        unsigned int current_start = 1;
                        unsigned int max_start = 0;
                        unsigned int max_end = 0;
                        float next_val = 0.0;
                        
                        for (int x = 1; x <= seq->seq_len; x++) {
                            next_val = ((read_cov[x] > cov_threshold) ? 1.0 : -1.0);
                            
                            // Is this a continuation of the current segmemt?
                            if (next_val + current_max > next_val) {
                                
                                current_max += next_val;
                                
                            } else {  // Or is it a reset to a possible new segment?
                                
                                current_max = next_val;
                                current_start = x;
                                
                                if (max_seen > 60) {  // Was the previous segment a significant region?
                                    pairing_problem *pp = malloc(sizeof(pairing_problem));
                                    if (!pp) {
                                        fprintf(stderr, "ERROR: malloc failure: pairing_problem\n");
                                        exit(EXIT_FAILURE);
                                    }
                                    pp->start = max_start;
                                    pp->end = max_end;
                                    pp->type = 'D';
                                    LL_PREPEND(seq->mp_issues, pp);
                                    x = max_end + 1;  // Start looking at the end of the last segment
                                    max_seen = 0.0;
                                    max_start = max_end = 0;
                                    
                                }
                            }
                            // Update the max segment
                            if (current_max > max_seen) {
                                max_seen = current_max;
                                max_start = current_start;
                                max_end = x;
                            }
                        }
                        // Get the last region if there was one
                        if (max_seen > 60) {  // Was this a substantial region?
                            pairing_problem *pp = malloc(sizeof(pairing_problem));
                            if (!pp) {
                                fprintf(stderr, "ERROR: malloc failure: pairing_problem\n");
                                exit(EXIT_FAILURE);
                            }
                            pp->start = max_start;
                            pp->end = max_end;
                            pp->type = 'D';
                            LL_PREPEND(seq->mp_issues, pp);
                        }
                    }
                    
                    // Look for misassemblies by inspecting the mate-pair analysis, if used
                    
                    if (mp_inserts->count && seq->mp_good_pairs) {
                        
                        float *physical_cov = seq->recon_array[recon_array_size-1];
                        
                        // This detects any regions of "physical coverage" where there is a dip below the moving average,
                        // which then recovers well above the moving average
                        
                        float prev = 0.0;
                        const unsigned int step = 75;
                        unsigned int front_max_pos = seq->seq_len;
                        unsigned int back_max_pos = 0;
                        float phys_cov_threshold = 0.0;
                        
                        // Find the top "level" of physical coverage at the beginning of the contig
                        for (int x = 1; x <= seq->seq_len; x += step) {
                            if ((physical_cov[x] < prev) && (physical_cov[x] > read_cov[x])) {
                                front_max_pos = x-step;
                                break;
                            } else {
                                prev = physical_cov[x];
                            }
                        }
                        
                        prev = 0.0;
                        // Find the top "level" of physical coverage at the end of the contig
                        for (int x = 1+(step * (seq->seq_len / step)); x >= 1; x -= step) {
                            if ((physical_cov[x] < prev) && (physical_cov[x] > read_cov[x])) {
                                back_max_pos = x+step;
                                break;
                            } else {
                                prev = physical_cov[x];
                            }
                        }
                        
                        if (front_max_pos < back_max_pos) {  // If these are equal, then there are no breaks, or this contig is too short
                            
                            float front_max_val = physical_cov[front_max_pos];
                            float back_max_val = physical_cov[back_max_pos];
                            phys_cov_threshold = (front_max_val + back_max_val) * 0.25;  // 50% of the mean of front and back maxima
                            
                            fprintf(stderr, "DIAG: %s pairing min/max %f %d %f %d thresh: %f\n",utstring_body(&seq->seq_id),front_max_val,front_max_pos,back_max_val,back_max_pos,(front_max_val + back_max_val) * 0.25);
                            
                            // Now find problem regions
                            float max_seen = 0.0;
                            float current_max = 0.0;
                            unsigned int current_start = front_max_pos+1;
                            unsigned int max_start = 0;
                            unsigned int max_end = 0;
                            float next_val = 0.0;
                            
                            for (int x = front_max_pos+1; x < back_max_pos; x++) {
                                next_val = ((physical_cov[x] <= phys_cov_threshold) ? 1.0 : -10.0);
                                
                                // Is this a continuation of the current segmemt?
                                if (next_val + current_max > next_val) {
                                    
                                    current_max += next_val;
                                    
                                } else {  // Or is it a reset to a possible new segment?
                                    
                                    current_max = next_val;
                                    current_start = x;
                                    
                                    if (max_seen > 60) {  // Was the previous segment a significant region?
                                        pairing_problem *pp = malloc(sizeof(pairing_problem));
                                        if (!pp) {
                                            fprintf(stderr, "ERROR: malloc failure: pairing_problem\n");
                                            exit(EXIT_FAILURE);
                                        }
                                        pp->start = max_start;
                                        pp->end = max_end;
                                        pp->type = 'F';
                                        LL_PREPEND(seq->mp_issues, pp);
                                        max_seen = 0.0;
                                        x = max_end + 1;  // Start looking at the end of the last segment
                                        max_start = max_end = 0;
                                    }
                                }
                                // Update the max segment
                                if (current_max > max_seen) {
                                    max_seen = current_max;
                                    max_start = current_start;
                                    max_end = x;
                                }
                            }
                            // Get the last region if there was one
                            if (max_seen > 60) {  // Was this a substantial region?
                                pairing_problem *pp = malloc(sizeof(pairing_problem));
                                if (!pp) {
                                    fprintf(stderr, "ERROR: malloc failure: pairing_problem\n");
                                    exit(EXIT_FAILURE);
                                }
                                pp->start = max_start;
                                pp->end = max_end;
                                pp->type = 'F';
                                LL_PREPEND(seq->mp_issues, pp);
                            }
                        }
                        
                        // This detects "uncovered ends", which are misassemblies that are relatively short end sequences
                        
                        if ((phys_cov_threshold > 10.0) && (seq->mp_insert_mean > 0)) {
                            for (int x = 1; x <= seq->seq_len; x++) {
                                // IF no mate-pairs span this location, but there is coverage here...
                                if ((physical_cov[x] == 0.0) && (read_cov[x] != 0.0)) {
                                    int start = x;
                                    // Then count the number of positions in a row where this is true
                                    while ((physical_cov[x] == 0.0) && (read_cov[x] != 0.0) && (x <= seq->seq_len)) {
                                        x++;
                                    }
                                    // If the number of positions is more than a threshold (currently 60) then note this....
                                    if ((x <= seq->seq_len) && (x-start >= 60)) {
                                        pairing_problem *pp = malloc(sizeof(pairing_problem));
                                        if (!pp) {
                                            fprintf(stderr, "ERROR: malloc failure: pairing_problem\n");
                                            exit(EXIT_FAILURE);
                                        }
                                        pp->start = start;
                                        pp->end = x-1;
                                        pp->type = 'Z';
                                        LL_PREPEND(seq->mp_issues, pp);
                                    }
                                }
                            }
                        }
                        
                        // This detects small inserts, which are misassemblies that improperly insert sequences.
                        // Look for runs of identical physical cov above some threshold.
                        // But only in sequences significantly longer than the estimated MP insert size.
                        
                        if ((phys_cov_threshold > 10.0) && (seq->seq_len > 2.5*seq->mp_insert_mean)) {
                            float prev = 0.0;
                            for (int x = 1; x <= seq->seq_len; x++) {
                                // If there is identical mp coverage, and positive read coverage here...
                                if ((physical_cov[x] == prev) && (read_cov[x] != 0.0)) {
                                    int start = x;
                                    // Then count the number of positions in a row where this is true
                                    while ((physical_cov[x] == prev) && (read_cov[x] != 0.0) && (x <= seq->seq_len)) {
                                        x++;
                                    }
                                    // If the number of positions is more than a threshold (currently 100) then note this....
                                    if ((x <= seq->seq_len) && (x-start >= 100)) {
                                        pairing_problem *pp = malloc(sizeof(pairing_problem));
                                        if (!pp) {
                                            fprintf(stderr, "ERROR: malloc failure: pairing_problem\n");
                                            exit(EXIT_FAILURE);
                                        } 
                                        pp->start = start;
                                        pp->end = x-1;
                                        pp->type = 'I';
                                        LL_PREPEND(seq->mp_issues, pp);
                                    }
                                }    
                                prev = seq->recon_array[recon_array_size-1][x];
                            }
                        }
                    }
                }
            }
        }
    }
    
    if (seq_ref->count) {
        
        pipe(pipe1);

#pragma omp parallel sections default(shared)       
        {
            
            
#pragma omp section 
            {   // reference fasta file reader
                UT_string *ref_file;
                utstring_new(ref_file);
                ss_strcat_utstring(ref_file, seq_ref->filename[0]);
                ss_stream_reader(ref_file, pipe1[1]);
                utstring_free(ref_file);
            }
            
#pragma omp section 
            {   // reference file parser, store nucleotide sequence for selected sequences
                
                FILE *ref_in = fdopen(pipe1[0],"r");
                int new_bases = 0;
                int bases_added = 0;
                int gc_sum = 0;
                int gc_cnt = 0;
                int base_index;                
                seq = NULL;
                int num_sel = HASH_COUNT(selected);
                int num_added = 0;
                UT_string *line;
                utstring_new(line);
                
                if (v->count) {
                    fprintf(stderr, "INFO: Loading sequence from external reference fasta file\n");
                }
                
                while (ss_get_utstring(ref_in, line)) {
                    if (*(utstring_body(line)) == '>') {
                        if (seq) {  // Store GC for previous sequence
                            if (gc_sum != 0) {
                                seq->gc_percent = (float) gc_sum / (float) gc_cnt;
                            } else {
                                seq->gc_percent = 0.0;
                            }
                            *ptr2 = '\0';
                            seq = NULL;
                        }
                        
                        char *save_ptr = NULL;
                        tokptr = strtok_r(utstring_body(line), "\t \n", &save_ptr);  // Tokenize on whitespace
                        tokptr++;  // skip '>'
                        
                        // Look to see if sequence is selected and marked to have
                        // its nucleotide sequence output
                        HASH_FIND_STR(selected, tokptr, seq);
                        if (seq) {
                            num_added++;
                            ptr2 = (char *) &(seq->recon_array[1][0]);  // Array to hold reference sequence
                            gc_sum = 0;
                            gc_cnt = 0;
                            bases_added = 0;
                        }
                    } else if (seq) {
                        ptr = utstring_body(line);  // Next part of reference sequence to add from external fasta file
                        new_bases = strlen(ptr) - 1;  // number of new bases without newline
                        if ((bases_added + new_bases) > seq->seq_len) {
                            fprintf(stderr, "ERROR: Sequence in external reference fasta file for %s longer than listed in catalog file.  %i bp > %i bp\n", utstring_body(&(seq->seq_id)), bases_added+new_bases, seq->seq_len);
                            exit(EXIT_FAILURE);
                        }
                        for (int i = 0; i < new_bases; i++) {
                            base_index = code[*ptr];
                            gc_cnt += (comp[base_index] != 0);
                            gc_sum += gc_tab[base_index];
                            *ptr2++ = *ptr++;
                        }
                        bases_added += new_bases;
                    }
                }
                fclose(ref_in);
                if (seq) {
                    if (gc_sum != 0) {  // store GC for previous sequence
                        seq->gc_percent = (float) gc_sum / (float) gc_cnt;
                    } else {
                        seq->gc_percent = 0.0;
                    }
                    *ptr2 = '\0';
                }
                
                utstring_free(line);
                
                if (v->count) {
                    fprintf(stderr, "INFO: Loaded nucleotide sequence from external reference fasta file for %i/%i selected sequences\n", num_added, num_sel);
                }
            }  // End omp section 
        }   // End omp parallel sections
    }
    
    if (v->count) {
        fprintf(stderr, "\nINFO: Total abundance = %.12f\n", total_abund);
        fprintf(stderr, "\nINFO: Num taxa selected = %d out of %d\n", num_selected, seq_tot);
        fprintf(stderr, "\nINFO: Total abundance = %f%%\n", total_fractional_abundance*100.0);
        fprintf(stderr, "INFO: Reads accounted for: %d out of %d (%.2f%%)\n", remaining_reads, num_reads, 100.0*(double)remaining_reads/(double)num_reads);
        fprintf(stderr, "INFO: Second chance reads %d\n", second_chance_reads);
    }
    
    JSON_BEG_OBJ;
    JSON_STR_PROP("SEASTAR_tool", argv[0]);  JSON_COMMA;
    JSON_STR_PROP("SEASTAR_version", SS_BUILD_VERSION);  JSON_COMMA;
    JSON_SUB_PROP("runtime_parameters");
    JSON_BEG_OBJ;

    JSON_LIT_PROP("all_taxa", nt->count ? "true" : "false");  JSON_COMMA;
    JSON_FLT_PROP("bit_thresh",bit_thresh);  JSON_COMMA;
    JSON_FLT_PROP("bit_fraction",bit_frac);  JSON_COMMA;
    JSON_INT_PROP("read_map_limit", read_cnt_limit);  JSON_COMMA;
    JSON_LIT_PROP("second_chance_reads", s->count ? "true" : "false");  JSON_COMMA;
    JSON_LIT_PROP("relax_read_sharing", r->count ? "true" : "false");  JSON_COMMA;
    JSON_LIT_PROP("absolute_bitscores", a->count ? "true" : "false");  JSON_COMMA;

    JSON_INT_PROP("cov_window", cov_window);  JSON_COMMA;

    if (detail->count) {
        JSON_N_STR_ARRAY_OBJ("detail_list",detail,sval);
    }
    if (detail_file->count) { JSON_STR_PROP("detail_file",detail_file->filename[0]); JSON_COMMA; }

    if (seq_recon->count) { 
        JSON_FLT_PROP("ambig_tol",ambig_tol->count ? ambig_tol->dval[0] : 0.2);  JSON_COMMA;
    }
    
    if (seq_ref->count) {
        JSON_STR_PROP("ref", seq_ref->filename[0]);  JSON_COMMA;
    }
    
    JSON_N_STR_ARRAY_OBJ("read_files",read_files,filename);
    JSON_N_STR_ARRAY_OBJ("rev_read_files",rev_read_files,filename);
    if (read_out->count) { 
        JSON_STR_PROP("read_out",read_out->sval[0]); JSON_COMMA;
        JSON_LIT_PROP("output_nonmatch", read_out_nomap->count ? "true" : "false");  JSON_COMMA;
        JSON_LIT_PROP("read_output_gzip", read_out_gz->count ? "true" : "false");  JSON_COMMA;
    }
    if (exclude->count) {
        JSON_N_STR_ARRAY_OBJ("exclude_list",exclude,sval);
    }
    if (exclude_file->count) { JSON_STR_PROP("exclude_file",exclude_file->filename[0]); JSON_COMMA; }
    JSON_LIT_PROP("invert_ex",invert_ex->count ? "true" : "false"); JSON_COMMA;
    JSON_STR_PROP("catalog",c->filename[0]);  JSON_COMMA;
    if (split->count) { JSON_INT_PROP("split", split->ival[0]);  JSON_COMMA; }
    JSON_LIT_PROP("separate_strands", sep->count ? "true" : "false");  JSON_COMMA;
    if (mp_analysis->count) { 
        JSON_FLT_PROP("mp_share_lim", share_lim);  JSON_COMMA;
        JSON_FLT_PROP("mp_mate_lim", mate_lim);  JSON_COMMA;
        JSON_LIT_PROP("mp_strict", mp_strict->count ? "true" : "false");  JSON_COMMA;
        JSON_LIT_PROP("mp_inserts", mp_inserts->count ? "true" : "false");  JSON_COMMA;
        JSON_LIT_PROP("mp_circ", mp_circ->count ? "true" : "false");  JSON_COMMA;
        JSON_LIT_PROP("detect_dups", detect_dups->count ? "true" : "false");  JSON_COMMA;
    }

    if (sim_frac->count) { JSON_FLT_PROP("sim_frac", sim_frac->dval[0]);  JSON_COMMA; }

    JSON_LIT_PROP("old_bwa_samse", f_bwa_samse->count ? "true" : "false");  JSON_COMMA;
    JSON_INT_PROP("seed", seed->count ? seed->ival : 12345);  JSON_COMMA;
    JSON_LIT_PROP("no_rand", no_rand->count ? "true" : "false");  JSON_COMMA;
    JSON_N_STR_ARRAY_OBJ("rev_align_files",fr,filename);
    JSON_N_STR_ARRAY_OBJ("align_files",f,filename);
    JSON_LIT_PROP("rollup", rollup->count ? "true" : "false");

    JSON_END_OBJ;

    if (selected && selected->bit_score >= bit_thresh) {
        int seen = 0;        
        JSON_COMMA;
        JSON_SUB_PROP("nodes");
        JSON_BEG_OBJ;
        
        qsort(sel_list, num_sel, sizeof(seq_entry *), selected_list_bitscore_cmp); 
        
        for (int z = 0; z < num_sel; z++) {    // Walk through all sequences
            
            if (!(seq = sel_list[z]) || (seq->abundance == 0.0)) continue;
            
            if (seen) JSON_COMMA;
            seen = 1;
            
            JSON_SUB_PROP(utstring_body(&seq->seq_id));
            JSON_BEG_OBJ; 
            JSON_FLT_PROP("bits", seq->bit_score);  JSON_COMMA;
            JSON_FLT_PROP("rd_cnt", seq->abundance*seq->seq_len);  JSON_COMMA;
            JSON_FLT_PROP("int_cov", seq->abundance);  JSON_COMMA;
            JSON_FLT15_PROP("rel_ab", seq->abundance/total_abund);  JSON_COMMA;
            JSON_FLT_PROP("cov", seq->coverage);  JSON_COMMA;
            JSON_FLT_PROP("rd_len", seq->coverage/seq->abundance);  JSON_COMMA;
            JSON_INT_PROP("seq_len", seq->seq_len);  JSON_COMMA;
            if (seq->gc_percent != 0.0) {
                JSON_FLT_PROP("pct_gc", 100.0*seq->gc_percent);  JSON_COMMA;
            }
            JSON_FLT_PROP("pct_uncov", 100.0*seq->uncovered);  JSON_COMMA;
            if (mp_analysis->count) {
                JSON_INT_PROP("mp_pairs", seq->mp_good_pairs);  JSON_COMMA;  
                JSON_INT_PROP("mp_sh", seq->num_shared);  JSON_COMMA; 
                JSON_INT_PROP("mp_fwd", seq->num_forward);  JSON_COMMA;
                JSON_INT_PROP("mp_bwd", seq->num_backward);  JSON_COMMA;
                if (seq->mp_issues) {
                    JSON_SUB_PROP("contig_problems");
                    pairing_problem *pp = NULL;
                    JSON_BEG_ARR;
                    LL_FOREACH(seq->mp_issues, pp) {
                        if (pp != seq->mp_issues) JSON_COMMA;
                        JSON_BEG_OBJ;
                        JSON_INT_PROP("start", pp->start);  JSON_COMMA;
                        JSON_INT_PROP("end", pp->end);  JSON_COMMA;
                        switch (pp->type) {
                            case 'F' :
                                JSON_STR_PROP("type", "Physical coverage break");
                                break;
                            case 'D' :
                                JSON_STR_PROP("type", "Collapsed duplication");
                                break;
                            case 'Z' :
                                JSON_STR_PROP("type", "Uncovered end");        
                                break;
                            case 'I' :
                                JSON_STR_PROP("type", "Small insert");        
                                break;
                            default :
                                JSON_STR_PROP("type", "Unknown!");        
                        }
                        JSON_END_OBJ;
                    }
                    JSON_END_ARR;  JSON_COMMA;    
                }
                if (seq->seq_len > 2*cov_window) {
                    if (mp_inserts->count) {
                        JSON_FLT_PROP("mp_ins_mean", seq->mp_insert_mean);  JSON_COMMA;
                        JSON_FLT_PROP("mp_ins_stdev", seq->mp_insert_stdev);  JSON_COMMA;
                    }
                    JSON_FLT_PROP("adj_cov", seq->adj_coverage_mean);  JSON_COMMA;
                    if (seq->recon_array) {
                        JSON_FLT_PROP("adj_cov_stddev", seq->adj_coverage_stddev);  JSON_COMMA;
                        JSON_FLT_PROP("adj_cov_min", seq->adj_coverage_min);  JSON_COMMA;
                        JSON_FLT_PROP("adj_cov_max", seq->adj_coverage_max);  JSON_COMMA;
                    }
                } else {
                    JSON_LIT_PROP("short_seq", "true");  JSON_COMMA;
                }
            }
            if (seq->recon_array && (detailed->count || seq->mp_issues)) {
                JSON_N_FLT_ARRAY_OBJ("per_nt_cov",seq->recon_array[0],seq->seq_len);  
                if (mp_inserts->count) {        
                    JSON_N_FLT_ARRAY_OBJ("per_nt_phys_cov",seq->recon_array[recon_array_size-1],seq->seq_len);  
                    JSON_N_FLT_ARRAY_OBJ("per_nt_mp_ins",seq->recon_array[recon_array_size-2],seq->seq_len);  
                }
            }
            JSON_STR_PROP("name", utstring_body(&seq->ref_seq->seq_name));  JSON_COMMA;
            JSON_STR_PROP("desc", utstring_body(&seq->ref_seq->seq_desc));
            if ((seq_recon->count || seq_ref->count) && seq->recon_array) {
                JSON_COMMA;
                JSON_STR_PROP("recon_seq", (char *) &(seq->recon_array[1][0]));
            }
            JSON_END_OBJ;
        }
        
        JSON_END_OBJ;
        
    }

    // Output edge data when mp analysis is enabled

    if (mp_analysis->count) {
        JSON_COMMA;
        JSON_SUB_PROP("shared_seq_edges");
        JSON_BEG_ARR;
        int seen = 0;
        for (int z = 0; z < num_sel; z++) {    // Walk through all sequences
            if (!(seq = sel_list[z]) || (seq->abundance == 0.0)) continue;
            
            // Shared sequence "edges"
            for (int x = 0; x < seq->num_shared; x++) {
                
                if (seq->shared[x].shared_bitscore < share_lim) {
                    break;
                }
                
                if (seq->shared[x].seq_ptr->abundance == 0.0)
                    continue;
                
                HASH_FIND_STR(selected, utstring_body(&seq->shared[x].seq_ptr->seq_id), s_ptr);
                
                if (s_ptr && s_ptr->num_shared) {
                    int y;
                    // Find the entry index (y) for this share in the other sequence's shared dataset
                    for (y = 0; (y < s_ptr->num_shared) && (seq != s_ptr->shared[y].seq_ptr); y++);
                    
                    assert(y < s_ptr->num_shared);
                    
                    if (seen) JSON_COMMA;
                    seen = 1;
                    JSON_BEG_OBJ;
                    JSON_STR_PROP("n1", utstring_body(&seq->seq_id));  JSON_COMMA;
                    JSON_STR_PROP("n2", utstring_body(&seq->shared[x].seq_ptr->seq_id));  JSON_COMMA;
                    JSON_FLT_PROP("bits", seq->shared[x].shared_bitscore);  JSON_COMMA;
                    JSON_INT_PROP("p1", seq->shared[x].mean_pos);  JSON_COMMA;
                    JSON_INT_PROP("p2", s_ptr->shared[y].mean_pos);
                    JSON_END_OBJ;
                }
            }
            seq->num_shared = 0;    // Don't make this edge again
        }
        JSON_END_ARR;
        
        JSON_COMMA;
        JSON_SUB_PROP("internal_edges");
        
        // This is to restore the original read_table for debugging
        // qsort(read_table+1, num_reads, sizeof(read_entry), read_entry_id_cmp);
        
        JSON_BEG_ARR;
        seen = 0;
        for (int z = 0; z < num_sel; z++) {    // Walk through all sequences
            if (!(seq = sel_list[z]) || (seq->abundance == 0.0)) continue;
            
            // Backward edges
            for (int x = 0; x < seq->num_backward; x++) {
                
                if (seq != seq->backward[x].seq_ptr) {
                    continue;   // Only handling self-links here 
                }
                
                if (seq->backward[x].shared_bitscore < mate_lim) {
                    break;
                }
                
                if (seq->num_forward) {
                    int y;
                    
                    // Find the entry index (y) for this backward in the other sequence's dataset
                    if (1) {
                        // if (seq->backward[x].mean_pos > 0) {  // Only output internal links from entirely "good" seqs
                        for (y = 0; (y < seq->num_forward) && (seq != seq->forward[y].seq_ptr); y++);
                        
                        if (seq == seq->forward[y].seq_ptr) {  // Self links...
                            if (seen) JSON_COMMA;
                            seen = 1; 
                            JSON_BEG_OBJ;
                            JSON_STR_PROP("n1", utstring_body(&seq->seq_id));  JSON_COMMA;
                            JSON_STR_PROP("n2", utstring_body(&seq->seq_id));  JSON_COMMA;
                            JSON_FLT_PROP("bits", seq->backward[x].shared_bitscore);  JSON_COMMA;
                            JSON_INT_PROP("p1", seq->forward[y].mean_pos);  JSON_COMMA;
                            JSON_INT_PROP("p2", seq->backward[x].mean_pos);
                            JSON_END_OBJ;
                        // } else {
                            // Something is wrong!
                           // fprintf(stderr, "\nWARNING: Sequence %s has internal pairing problems!\n", utstring_body(&seq->seq_id));
                        }
                    }
                } else {
                    // Something is wrong!
                    fprintf(stderr, "\nWARNING: Sequence %s has internal pairing problems!\n", utstring_body(&seq->seq_id));
                }
            }
        }
        JSON_END_ARR;
        
        JSON_COMMA;
        JSON_SUB_PROP("edges");
        JSON_BEG_ARR;
        seen = 0;
        for (int z = 0; z < num_sel; z++) {    // Walk through all sequences
            if (!(seq = sel_list[z]) || (seq->abundance == 0.0)) continue;
            
            // Backward edges
            for (int x = 0; x < seq->num_backward; x++) {
                
                if (seq->backward[x].shared_bitscore < mate_lim) {
                    break;
                }
                
                if ((seq == seq->backward[x].seq_ptr) || (seq->backward[x].seq_ptr->abundance == 0.0)) {
                    continue;   // Self-links have already been handled
                }
                
                HASH_FIND_STR(selected, utstring_body(&seq->backward[x].seq_ptr->seq_id), s_ptr);
                
                if (s_ptr && s_ptr->num_backward) {
                    int y;
                    
                    if (seen) JSON_COMMA;
                    seen = 1; 
                    JSON_BEG_OBJ;
                    
                    // Find the entry index (y) for this backward in the other sequence's dataset
                    if (seq->backward[x].mean_pos > 0) {
                        for (y = 0; (y < s_ptr->num_forward) && (seq != s_ptr->forward[y].seq_ptr); y++);    // Note reversal below to make BF into FB
                        
                        JSON_STR_PROP("n1", utstring_body(&seq->backward[x].seq_ptr->seq_id));  JSON_COMMA;
                        JSON_STR_PROP("n2", utstring_body(&seq->seq_id));  JSON_COMMA;
                        JSON_STR_PROP("dir", "FB");  JSON_COMMA;
                        JSON_FLT_PROP("bits", seq->backward[x].shared_bitscore);  JSON_COMMA;
                        JSON_INT_PROP("p1", s_ptr->forward[y].mean_pos);  JSON_COMMA;
                        JSON_INT_PROP("p2", seq->backward[x].mean_pos);
                    } else {
                        for (y = 0; (y < s_ptr->num_backward) && (seq != s_ptr->backward[y].seq_ptr); y++);
                        
                        assert(y < s_ptr->num_backward);
                        
                        JSON_STR_PROP("n1", utstring_body(&seq->seq_id));  JSON_COMMA;
                        JSON_STR_PROP("n2", utstring_body(&seq->backward[x].seq_ptr->seq_id));  JSON_COMMA;
                        JSON_STR_PROP("dir", "BB");  JSON_COMMA;
                        JSON_FLT_PROP("bits", seq->backward[x].shared_bitscore);  JSON_COMMA;
                        JSON_INT_PROP("p1", seq->backward[x].mean_pos);  JSON_COMMA;
                        JSON_INT_PROP("p2", s_ptr->backward[y].mean_pos);
                    }
                    JSON_END_OBJ;
                }
            }
            
            // Forward edges
            for (int x = 0; x < seq->num_forward; x++) {
                
                if (seq->forward[x].shared_bitscore < mate_lim) {
                    break;
                }
                
                if ((seq == seq->forward[x].seq_ptr) || (seq->forward[x].seq_ptr->abundance == 0.0)) {
                    continue;   // Self-links have already been handled
                }
                
                HASH_FIND_STR(selected, utstring_body(&seq->forward[x].seq_ptr->seq_id), s_ptr);
                
                if (s_ptr && s_ptr->num_forward && (seq != seq->forward[x].seq_ptr)) {
                    
                    int y;
                    
                    if (seen) JSON_COMMA;
                    seen = 1; 
                    JSON_BEG_OBJ;
                    
                    // Find the entry index (y) for this share in the other sequence's dataset
                    if (seq->forward[x].mean_pos > 0) {
                        // Self links are handled in the backward section above...
                        for (y = 0; (y < s_ptr->num_backward) && (seq != s_ptr->backward[y].seq_ptr); y++);
                        
                        assert(y < s_ptr->num_backward);
                        
                        JSON_STR_PROP("n1", utstring_body(&seq->seq_id));  JSON_COMMA;
                        JSON_STR_PROP("n2", utstring_body(&seq->forward[x].seq_ptr->seq_id));  JSON_COMMA;
                        JSON_STR_PROP("dir", "FB");  JSON_COMMA;
                        JSON_FLT_PROP("bits", seq->forward[x].shared_bitscore);  JSON_COMMA;
                        JSON_INT_PROP("p1", seq->forward[x].mean_pos);  JSON_COMMA;
                        JSON_INT_PROP("p2", s_ptr->backward[y].mean_pos);
                        
                    } else {
                        for (y = 0; (y < s_ptr->num_forward) && (seq != s_ptr->forward[y].seq_ptr); y++);
                        
                        assert(y < s_ptr->num_forward);
                        
                        JSON_STR_PROP("n1", utstring_body(&seq->seq_id));  JSON_COMMA;
                        JSON_STR_PROP("n2", utstring_body(&seq->forward[x].seq_ptr->seq_id));  JSON_COMMA;
                        JSON_STR_PROP("dir", "FF");  JSON_COMMA;
                        JSON_FLT_PROP("bits", seq->forward[x].shared_bitscore);  JSON_COMMA;
                        JSON_INT_PROP("p1", seq->forward[x].mean_pos);  JSON_COMMA;
                        JSON_INT_PROP("p2", s_ptr->forward[y].mean_pos);
                    } 
                    JSON_END_OBJ;
                }
            }
            seq->num_forward = seq->num_backward = 0;    // Don't make these edges again
        }
        JSON_END_ARR;
    }

    // Output rollup stats when splitting is on and/or the roll-up flag is set

    if (split_len || rollup->count) {
        
        JSON_COMMA;
        JSON_SUB_PROP("rollup_stats");
        JSON_BEG_OBJ;
        
        total_abund = 0.0;
    #pragma omp parallel for default(none) shared(cat_array, n_taxa) private(entry, seq) reduction(+ : total_abund) schedule(static)
        for (int z = 0; z < n_taxa; z++) {
            entry = cat_array[z];
            entry->ref_coverage = 0.0;
            entry->ref_abundance = 0.0;            
            entry->ref_bit_score = 0.0;
            entry->ref_gc_percent = 0.0;
            DL_FOREACH(entry->sequences, seq) {
                entry->ref_coverage += seq->coverage * (double) seq->seq_len;
                entry->ref_abundance += seq->abundance * (double) seq->seq_len;
                entry->ref_bit_score += seq->bit_score * (double) seq->seq_len;
                entry->ref_gc_percent += seq->gc_percent * (double) seq->seq_len;
            }
            entry->ref_coverage /= (double) entry->ref_seq_len;
            entry->ref_abundance /= (double) entry->ref_seq_len;
            total_abund += entry->ref_abundance;
            entry->ref_bit_score /= (double) entry->ref_seq_len;
            entry->ref_gc_percent /= (double) entry->ref_seq_len;
        }
        
        if (v->count) {
            fprintf(stderr, "\nINFO: Total abundance (again) = %.12f\n", total_abund);
        }
        // Sort by bit_score
        
        qsort(cat_array, n_taxa, sizeof(ref_seq_entry *), catalog_entry_bitscore_cmp);
        
        for (int z = 0; (z < n_taxa) && (cat_array[z]->ref_bit_score > bit_thresh); z++) {
            
            entry = cat_array[z];
            
            if (z) JSON_COMMA;
            JSON_SUB_PROP(utstring_body(&entry->seq_name));
            JSON_BEG_OBJ;
            JSON_SUB_PROP("node_list");
            JSON_BEG_ARR;
            DL_FOREACH(entry->sequences, seq) {
                if (entry->sequences != seq) JSON_COMMA;
                JSON_STR_ELEM(utstring_body(&seq->seq_id));    
            }
            JSON_END_ARR;  JSON_COMMA;
            
            JSON_FLT_PROP("bits", entry->ref_bit_score);  JSON_COMMA;
            JSON_FLT_PROP("rd_count", entry->ref_abundance*entry->ref_seq_len);  JSON_COMMA;
            JSON_FLT_PROP("int_cov", entry->ref_abundance);  JSON_COMMA;
            JSON_FLT_PROP("rel_ab", entry->ref_abundance/total_abund);  JSON_COMMA;
            JSON_FLT_PROP("cov", entry->ref_coverage);  JSON_COMMA;
            JSON_FLT_PROP("rd_len", entry->ref_coverage/entry->ref_abundance);  JSON_COMMA;
            JSON_INT_PROP("seq_len", entry->ref_seq_len);  JSON_COMMA;
            if (entry->ref_gc_percent != 0.0) {
                JSON_FLT_PROP("pct_gc", 100.0*entry->ref_gc_percent);  JSON_COMMA;
            }
            JSON_FLT_PROP("pct_uncov", 100.0*entry->uncovered/entry->ref_seq_len);  JSON_COMMA;
            JSON_STR_PROP("desc", utstring_body(&entry->seq_desc));
            JSON_END_OBJ;
        }
        JSON_END_OBJ;
    }

    JSON_COMMA;
    JSON_SUB_PROP("run_stats");
    JSON_BEG_OBJ;
    JSON_INT_PROP("num_selected_nodes", num_selected);  JSON_COMMA;
    JSON_INT_PROP("original_nodes", seq_tot);  JSON_COMMA;
    JSON_FLT_PROP("total_node_abundance", total_fractional_abundance);  JSON_COMMA;
    if (detect_dups->count) {
        JSON_FLT_PROP("mean_coverage", mean_coverage);  JSON_COMMA;
        JSON_FLT_PROP("coverage_duplication_threshold", cov_threshold);  JSON_COMMA;
    }
    JSON_INT_PROP("reads_assigned_to_selected", remaining_reads);  JSON_COMMA;
    JSON_INT_PROP("original_reads_aligned", num_reads);  JSON_COMMA;
    JSON_INT_PROP("second_chance_assigned_reads", second_chance_reads);
    JSON_END_OBJ;

    JSON_END_OBJ;

    printf("\n");
        
    free(sel_list);    
    
    utstring_free(str);
    utstring_free(str2);
    utstring_free(str3);
    utstring_free(str4);
    utstring_free(name_str);
    utstring_free(id_str);
    utstring_free(read_out_fn);
    utstring_free(file_name);
    
    exit(EXIT_SUCCESS);
}

///////////////////////////////////////////////////////////////////
// Integer comparison function for qsort / bsearch
///////////////////////////////////////////////////////////////////

int read_entry_count_cmp(const void *ptr1, const void *ptr2) {
    int sc = ((read_entry *)ptr1)->seq_count - ((read_entry *)ptr2)->seq_count;
    if (!sc) {
        // This ensures the sort stable when seq_count is equal.  Sort by original read_id
        sc = ((read_entry *)ptr1)->read_name_ptr->read_id - ((read_entry *)ptr2)->read_name_ptr->read_id;
    }
    return sc;
}

// This is for debugging only 
// int read_entry_id_cmp(const void *ptr1, const void *ptr2) {
//     int sc = ((read_entry *)ptr1)->read_name_ptr->read_id - ((read_entry *)ptr2)->read_name_ptr->read_id;
//    return sc;
// }

int catalog_entry_bitscore_cmp(const void *ptr2, const void *ptr1) {
    
    int retval = 0;
    float sc = (*((ref_seq_entry **) ptr1))->ref_bit_score - (*((ref_seq_entry **) ptr2))->ref_bit_score;
    
    if (sc == 0.0) {
        // This ensures the sort is stable when seq_count is equal.  Secondary sort by (unique) sequence name
        retval = strcmp(utstring_body(&(*((ref_seq_entry **) ptr1))->seq_name), utstring_body(&(*((ref_seq_entry **) ptr2))->seq_name));
    } else {
        retval = (sc > 0.0) ? 1 : -1;
    }
    
    return retval;
}

int selected_list_bitscore_cmp(const void *ptr2, const void *ptr1) {
    
    int retval = 0;
    float sc = (*((seq_entry **) ptr1))->bit_score - (*((seq_entry **) ptr2))->bit_score;
    
    if (sc == 0.0) {
        // This ensures the sort is stable when seq_count is equal.  Secondary sort by (unique) sequence name
        retval = strcmp(utstring_body(&(*((seq_entry **) ptr1))->seq_id), utstring_body(&(*((seq_entry **) ptr2))->seq_id));
    } else {
        retval = (sc > 0.0) ? 1 : -1;
    }
    
    return retval;
}            

int seq_entry_bitscore_cmp(const void *ptr2, const void *ptr1) {
    
    int retval = 0;
    float sc = ((seq_list_entry *) ptr1)->shared_bitscore - ((seq_list_entry *) ptr2)->shared_bitscore;
    
    if (sc == 0.0) {
        // This ensures the sort is stable when seq_count is equal.  Sort by original read_id
        retval = ((seq_list_entry *) ptr1)->seq_num - ((seq_list_entry *) ptr2)->seq_num;
    } else {
        retval = (sc > 0.0) ? 1 : -1;
    }
    
    return retval;
}

/////////////////////////////////////////////////////////////////
// 
// sam_stream_reader
//
// This function opens a SAM format text file (gzipped or not) and 
// converts each line to a simpler format mirroring the original
// BWA "samse" output format.
//
// The resulting data is written to another file handle (eg. a pipe)
// This is useful for easily implementing a SAM "reader thread".
//

long unsigned int sam_stream_reader(UT_string *fn, int pipe_fd) {
    
    long unsigned int cnt = 1;
    
    char *ptr = NULL;
    char *tok_ptr = NULL;
    char *ref_ptr = NULL;
    char *pos_ptr = NULL;
    char *bwa_ed_ptr = NULL;
    char *ed_ptr = NULL;
    char *xa_ptr = NULL;
    
    int bitfield = 0;
    int ref_cnt = 0;
    int mapped_cnt = 0;
    int field_cnt = 0;
    
    char *save_ptr = NULL;
    
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
        
        char *input_buf = utstring_body(data);
        unsigned int length = utstring_len(data);
        
        if (input_buf[0] != '@') {    
            
            input_buf[--length] = '\0'; // Get rid of the \n
            
            // Split on spaces and tabs
            strtok_r(input_buf, "\t", &save_ptr);    // input_buf is now the read_id only
            
            if (!(tok_ptr = strtok_r(NULL, "\t", &save_ptr))) {      // get the second field
                fprintf(stderr, "ERROR: SAM File format error! %s  Line: %lu  %s\n", utstring_body(fn), cnt, input_buf);
                fprintf(stderr, "Missing bitfield (col 2)\n"); 
                exit(EXIT_FAILURE);    
            }    
            
            bitfield = atoi(tok_ptr);
            
            if (!(bitfield & 4)) {    // Continue only if this is a mapped read (4 bit is not set)
                
                mapped_cnt++;
                ref_cnt = 1;
                bwa_ed_ptr = NULL;
                ed_ptr = NULL;
                xa_ptr = NULL;
                
                if (!(ref_ptr = strtok_r(NULL, "\t", &save_ptr))) {  // Get the reference ID
                    fprintf(stderr, "ERROR: SAM File format error! %s  Line: %lu  %s\n", utstring_body(fn), cnt, input_buf);
                    fprintf(stderr, "Missing reference ID (col 3)\n");
                    exit(EXIT_FAILURE);    
                }
                
                if (!(pos_ptr = strtok_r(NULL, "\t", &save_ptr))) {    // Get the reference position
                    fprintf(stderr, "ERROR: SAM File format error! %s  Line: %lu  %s\n", utstring_body(fn), cnt, input_buf);
                    fprintf(stderr, "Missing reference position (col 4)\n");
                    exit(EXIT_FAILURE);    
                }
                
                field_cnt = 4;
                
                // Now loop through tokens looking for "XM:i:", "NM:i:", "XA:Z:"
                while ((tok_ptr = strtok_r(NULL, "\t", &save_ptr))) {
                    field_cnt++; 
                    if ((field_cnt < 12) || ((*tok_ptr != 'X') && (*tok_ptr != 'N'))) continue;    // cheap early out
                    
                    if (!bwa_ed_ptr && (!strncmp(tok_ptr, "XM:i:", 5))) {    // this is the BWA edit distance for the first match
                        bwa_ed_ptr = tok_ptr + 5; 
                    } else if (!ed_ptr && (!strncmp(tok_ptr, "NM:i:", 5))) {    // this is edit distance we use for non-BWA aligners
                        ed_ptr = tok_ptr + 5; 
                    } else if (!xa_ptr && (!strncmp(tok_ptr, "XA:Z:", 5))) { // this is the match info for all other matches (BWA)
                        xa_ptr = tok_ptr + 5; 
                        
                        // Count semicolons to figure out how many more references there are
                        for (ptr = xa_ptr; *ptr; ptr++) {
                            ref_cnt += (*ptr == ';');
                        }
                    }
                }
                
                if (! (ed_ptr || bwa_ed_ptr)) {
                    fprintf(stderr, "ERROR: SAM File format error! %s  Line: %lu  %s\n", utstring_body(fn), cnt, input_buf);
                    fprintf(stderr, "Missing edit distance (XM:i:N or NM:i:N)\n");
                    //                    fprintf(stderr, "%s\n",utstring_body(keep));
                    exit(EXIT_FAILURE);    
                }
                
                fprintf(pipe_in,">%s %d %d\n", input_buf, ref_cnt, ref_cnt);
                
                if (bwa_ed_ptr) {
                    ed_ptr = bwa_ed_ptr;
                }
                
                fprintf(pipe_in,"%s\t%c%s\t%s\n", ref_ptr, (bitfield&16)?'-':'+', pos_ptr, ed_ptr);    // 16 bit in bitfield is '-' strand when set
                
                if (xa_ptr) {
                    
                    // Count semicolons in the XA string to determine
                    // This string needs to be split on semicolons
                    ref_ptr = strtok_r(xa_ptr, ",", &save_ptr);
                    
                    do {
                        pos_ptr = strtok_r(NULL, ",", &save_ptr); // position, with sign
                        ptr = strtok_r(NULL, ",", &save_ptr);     // Skip
                        ed_ptr = strtok_r(NULL, ";", &save_ptr);  // edit distance
                        
                        if (!(pos_ptr && ed_ptr)) {
                            fprintf(stderr, "ERROR: SAM File format error! %s  Line: %lu  %s\n", utstring_body(fn), cnt, input_buf);
                            fprintf(stderr, "Missing position or edit distance (XA:Z:A,B,C,D;...)\n");
                            exit(EXIT_FAILURE);    
                        }
                        
                        fprintf(pipe_in,"%s\t%s\t%s\n", ref_ptr, pos_ptr, ed_ptr);    
                        
                    } while ((ref_ptr = strtok_r(NULL, ",", &save_ptr)));
                }
            }
        } 
        cnt++;
    }
    
    fclose(pipe_in);
    gzclose(file);
    utstring_free(data);
    
    return(cnt);
}

/*************************************************************************
 create_catalog_entry
**************************************************************************/

int create_catalog_entry(char *seq_id_str, int seq_len, char *ref_seq_id_str,
                         char *desc_str, exclude_entry **exclude_table,
                         detail_entry **detail_table, seq_entry **seq_table,
                         ref_seq_entry **catalog, int invert_ex,
                         int sep, int det_all, int split_len) {

    int include_seq = 1;
    
    // int exclude = (*exclude_table != NULL);
    
    seq_entry *seq = NULL;
    
    exclude_entry *ex_ptr = NULL;
    
    ref_seq_entry *entry = NULL;
    
    if (*exclude_table) {
        // Look to see if this is an excluded sequence
        HASH_FIND_STR(*exclude_table, seq_id_str, ex_ptr);
        include_seq = ((!invert_ex && !ex_ptr) || (invert_ex && ex_ptr));
    }
    
    // If not excluded, then add it!
    if (include_seq) {
        // If this is a whole line
        // Make a new seq cat struct
        
        // If no splitting, split len is seq_len
        int local_split = split_len ? split_len : seq_len;
        
        // For each split "chunk" of the reference seq, create a sequence entry
        // Note! This can be due to splitting, or stranded analysis, which are
        // mututally exclusive.
        
        int last_chunk = 0;
        
        if (split_len) {
            last_chunk = ((seq_len / local_split) - 1 + (seq_len < local_split));
        } else if (sep) {
            last_chunk = 1;
        }
        
        // Alloc as many records as there are "chunks"
        if (!(seq = calloc((size_t)last_chunk+1, sizeof(seq_entry)))) {
            fprintf(stderr, "ERROR: calloc failed: sequence catalog creation\n");
            exit(EXIT_FAILURE);
        }
        
        for (unsigned int chunk = 0; chunk <= last_chunk; chunk++) {
            
            ss_strcat_utstring(&seq->seq_id, seq_id_str);
            
            if (split_len) { // Only add suffix when splitting is active
                if (chunk == last_chunk) {
                    // Add $ to end of the chunk id indicating last chunk
                    utstring_printf(&seq->seq_id, "|%u$", chunk);
                } else {
                    utstring_printf(&seq->seq_id, "|%u", chunk);
                }
            } else if (sep && chunk) {
                ss_strcat_utstring(&seq->seq_id, "|R");
            }
            
            // If this is the last chunk, calculate the true sequence length
            
            if ((chunk == last_chunk) || (sep)) {
                // If last chunk or separate strand analysis
                if (seq_len <= local_split) {
                    seq->seq_len = seq_len;     // Last chunk is only chunk
                    seq->offset = 0;
                } else {
                    seq->seq_len = local_split + (seq_len % local_split); // Else, last chunk + remainder
                    seq->offset = chunk * local_split;
                }
            } else {
                seq->seq_len = local_split; // Else, full chunk
                seq->offset = chunk * local_split;
            }
            
            seq->bit_score = -1.0;    // Initialize to invalid score
            seq->abs_bit_score = -1.0;    // Initialize to invalid score
            seq->org_bit_score = -1.0;    // Initialize to invalid score
            
            ss_HASH_ADD_UTSTR_STRUCT(*seq_table, seq_id, seq);
            
            // Search to see if a refseq matching this seq is already alloc'ed
            HASH_FIND_STR(*catalog, ref_seq_id_str, entry);
            
            if (entry) {  // Already seen ref sequence
                
                seq->ref_seq = entry;
                entry->ref_seq_len += seq->seq_len;
                DL_PREPEND(entry->sequences, seq);
                
            } else {    // ref seq is new
                
                // Allocate and populate a ref seq catalog record
                if (!(entry = calloc(1, sizeof(ref_seq_entry)))) {
                    fprintf(stderr, "ERROR: calloc failed: refseq catalog creation\n");
                    exit(EXIT_FAILURE);
                }
                
                seq->ref_seq = entry;
                entry->ref_seq_len = seq->seq_len;
                DL_PREPEND(entry->sequences, seq);
                
                ss_strcat_utstring(&entry->seq_name, ref_seq_id_str);
                
                if (desc_str) {
                    ss_strcat_utstring(&entry->seq_desc, desc_str);
                } else {
                    utstring_init(&entry->seq_desc);
                }
                
                // Add record to the catalog hash table
                ss_HASH_ADD_UTSTR_STRUCT(*catalog, seq_name, entry);
            }
            
            seq++;  // Move to next seq record in the calloced array.
        }
    }
    
    return (0);
}
        

/////////////////////////////////////////////////////////////////
//
// sam_header_stream_reader
//
// This function opens a sam file (gzipped or not) and shoves
// the resulting header lines into another file handle (eg. a pipe)
// This function stops reading as soon as the first non-header line
// (doesn't begin with '@') is encountered.
//

long unsigned int sam_header_stream_reader(UT_string *fn, int pipe_fd) {
    
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
    
    while (ss_gzget_utstring(file, data) && utstring_body(data)[0] == '@') {
        fputs(utstring_body(data), pipe_in);
        cnt++;
    }
    
    fclose(pipe_in);
    gzclose(file);
    utstring_free(data);
    
    return(cnt);
}

