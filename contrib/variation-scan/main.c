#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <time.h>
#include <ctype.h>
#include <pwd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/times.h>

#define BASE_STR_LEN 10
#define MAX_HOSTNAME 256
#define ALPHABET_SIZE 93
#define GetIndex(c) ((int)c - 33)
#define GEN 1
#define STR 2
#define TRI 3
#define VAR 4
#define STRLIST 5
#define SITE 6
#define SCAN 7
#define VSCAN 8
#define RNG 9


typedef struct _stdvar{
  int  id;
  void *mem;
  struct _stdvar *next;
} memstd;

typedef struct _node{
  char  *info;
  int   leaf;
  struct _node *children[ALPHABET_SIZE];
} TRIE;


typedef struct _string {
  char   *buffer;
  size_t length; //Total length of string
  size_t size;   //Current size of string
} string;

typedef struct _variant {
  string    *chromosome;
  string         *start;
  string           *end;
  char        strand[2];
  string            *id;
  string            *SO;
  string     *reference;
  string       *alleles;
  string          *freq;
  struct _variant *prev;
  struct _variant *next;
} variant;

typedef struct _range {
  int start[2]         ;
  int end[2]           ;
  int length[2]        ;
  int left_flank[2]    ;
  int right_flank[2]   ;
  variant *var_info ;
  struct _range *prev;
  struct _range *next;
} range;

typedef struct _stringlist {
  string            *element;
  struct _stringlist   *next;
} stringlist;

typedef struct _site {
  string   *sequence;
  float      weight;
  float        pval;
} site;

typedef struct _scan {
  int           offset;
  site              *D;
  site              *R;
  struct _scan   *next;
} scan;

typedef struct _varscan{
  variant *variation;
  scan    *scan_info;
  scan    *bestD;
  scan    *bestR;
  struct _varscan  *next;
} varscan;

typedef struct _lth {
  float score;
  float wdiff;
  float pratio;
} lth;

typedef struct _uth {
  float pval;
} uth;

typedef struct _threshold {
  struct _lth lower;
  struct _uth upper;
} threshold;

int CheckConvertMatrixOutput(char *out_dir, char *matrix_filename);
char *isStringListSame(char *stringtosplit, char separator);
threshold *sethreshold(threshold *cutoff, char *type, char *value);
void printHeaderSingleVariants(string *varscanFile, char*input, char *output, string *CMD);
void printHeaderHaplotypes(string *varscanFile, char*input, char *output, string *CMD);
void ScanHaplosequences(string *mscanquick_file, string *varscanFile, char *matrix_name, int matrix_size, threshold *cutoff);
void ScanSingleVariants(string *mscanquick_file, string *varscanFile, char *matrix_name, int matrix_size, threshold *cutoff);
void processLocus(varscan *locus, char *matrix_name, threshold *cutoff,FILE *fout);
void processHaplotypes(varscan *locus, range *intersect, int is_indel,int nb_sites1, int nb_sites2, int matrix_size ,char *matrix_name, threshold *cutoff, FILE *fout);
void compareAlleles(varscan *firstvar,varscan *secndvar,site *firstsite,site *secndsite,int first_offset,int second_offset,char *strand,threshold *cutoff,FILE *fout,char *matrix_name);
void compareHapAlleles(range *intersect,int is_indel,int matrix_size,varscan *firstvar,varscan *secndvar,site *firstsite,site *secndsite,int first_offset,int second_offset,char *strand,threshold *cutoff,FILE *fout,char *matrix_name);
void printVarScan (FILE *fout,char *matrix_name,varscan *firstvar, variant *bestvariant, variant *worstvariant,site *bestsite, site *worstsite, float bestpval, float worstpval,int bestoffset,int worstoffset, char *strand);
void varscanfree(varscan *delete);
varscan *varscanew(void);
varscan *varscanewToList(memstd **List);
varscan *varscanadd(varscan *group);
void varscanend(varscan *group);
void sitefree(site *delete);
site *sitenew(void);
site *sitenewToList(memstd **List);
site *sitefill(site *element, char *seq, char *weight, char *pval);
void scanfree(scan *delete);
scan *scanew(void);
scan *scanewToList(memstd **List);
scan *scanadd(scan *group);
scan *scanfill(scan *element,char *offset_scan, char *sequence_D, char *weight_D, char *pval_D,char *sequence_R, char *weight_R, char *pval_R );
void scanend(scan *group);
string *GetVariantHapIndex(string *offset_and_length,char *alleles,char *sequence);
string *GetVariantIndex(string *offset_and_length,char *sequence);
char *splitstr(char *str_to_split, char sep);
char *SplitOffsetFromTotalVars(char *variants_info);
//void CreateFastaFromHaplosequences(string *line ,char **token, unsigned long int *nb_variation, unsigned long int top_variation, int *nb_seq, FILE *fh_varsequence,FILE *fh_fasta_sequence);
//void CreateFastaFromSingleVariants(string *line, char **token, unsigned long int *nb_variation, unsigned long int top_variation, int *nb_seq, FILE *fh_varsequence, FILE *fh_fasta_sequence);
void CreateFastaFromVarseqHaplotypes(string *varsequence, string *fasta_sequence , int matrix_size, unsigned long int *nb_variation, unsigned long int *top_variation, unsigned long int *nb_seq);
void CreateFastaFromVarseqVariants(string *varsequence, string *fasta_sequence , int matrix_size, unsigned long int *nb_variation, unsigned long int *top_variation, unsigned long int *nb_seq);
void CreateFastaFromFastaHaplotypes(string *prev_fasta, string *new_fasta , int matrix_size, unsigned long int *nb_variation, unsigned long int *top_variation, unsigned long int *nb_seq);
void CreateFastaFromFastaVariants(string *prev_fasta, string *new_fasta , int matrix_size, unsigned long int *nb_variation, unsigned long int *top_variation, unsigned long int *nb_seq);
TRIE *ListInputMatrixFormats(void);
string *GetProgramPath(string *program_path, char *program_name, int die_on_error, stringlist *preferred_path);
void argToList(stringlist *list, char *value);
void strlistfree(stringlist *delete);
void strlistend(stringlist *group);
stringlist *strlistnew(void);
stringlist *strlistnewToList(memstd **List);
stringlist *strlistadd(stringlist *group);
void rangefree(range *delete);
range *rangenew(void);
range *rangenewToList(memstd **List);
range *rangeadd(range *group);
range *rangefill(range *element,int start1,int start2,int end1,int end2,int length1,int length2,int left_flank1,int left_flank2,int right_flank1,int right_flank2);
void rangend(range *group);
void varfree(variant *delete);
variant *varnew(void);
variant *varnewToList(memstd **List);
variant *varadd(variant *group);
variant *varfill(variant *element,char *chrInfo, char *startInfo, char *endInfo, char *strandInfo, char *idInfo, char *SOInfo, char *refInfo, char *allelesInfo, char *freqInfo);
void varend(variant *group);
void cgiMessage(char *message_type,char *color,char *fmt,va_list ap);
void cgiWarning(char *fmt,va_list ap);
void RsatInfo(char *fmt, ...);
void RsatWarning(char *fmt, ...);
void cgiError(char *fmt,va_list ap);
void RsatFatalError(char *fmt, ...) ;
void *_malloc(size_t size,char *func);
void *_realloc(void *ptr,size_t size,char *func);
void TrieInsert(TRIE *trie,char *key,char *value);
char* TrieSearch(TRIE *trie,char *key);
TRIE *TrieStart(void);
TRIE *TrieNew(void);
void TrieEnd(TRIE *trie);
void strfree(string *delete);
char *strccat(string *destn,char *fmt, ...);
char *strcopy(string *destn,char *orign);
char *strfmt(string *strn, char *fmt, ...);
string *strnew(void);
string *strnewToList(memstd **List);
string *_strmalloc(size_t size,char *func);
string *strlimt(string *totest);
char **initokadd(string *line,char **token,int numtok);
char **getokens(int numtok);
char *stralloc(string *resize,int size);
void *_MemTrackMalloc(size_t size,memstd **list, char *func);
void rstdlist(memstd *list);
memstd *MemTrackElmt(void);
memstd *MemTrackNew(void);
memstd *MemTrackAdd(void *ptr, int type, memstd **list);
memstd *MemTrackUpd(void *old, void *new, memstd *list);
memstd *ListStdUpd(void *lastptr,void *currptr,memstd *list);
memstd *rstdelem(void *ptr,memstd *list);
void rlist(memstd *list);
memstd *relem(void *ptr,memstd *list);
void InitMemLists(void);
void ReadProperties(void);
time_t InitRSAT(char *program,char *cmd);
int UpdateExecTimeLogFile(time_t start_time,time_t done_time);
int CheckValOpt(char *value[]);
int doit(char *command,int dry,int die_on_error,int verbose,int batch,char *job_prefix,FILE *log_handle,FILE *err_handle,char *cluster_queue);
string *AlphaDate(string *strtime,time_t *start_time);
time_t StartScript(char *program,char *cmd);
void ReportExecutionTime(time_t start_time);
FILE *OpenInputFile(FILE *filehandle,char *filename);
FILE *OpenOutputFile(FILE *filehandle,char *filename);
FILE *OpenAppendFile(FILE *filehandle,char *filename);
int CheckOutDir(string *output_dir,mode_t Umask, mode_t Chmod);
string *SplitFileName(string *token[],char *filename);
string *Get_pub_temp(string *public_temp_dir);
string *Get_temp_dir(string *tmp_dir);
string *make_temp_file(string *tmp_file,char *tmp_dir,char *tmp_prefix,int add_date,int make_dir, int protect);
string *limlinetok(string *line, char **token,int numtok);
string *Get_contigs_file(string *contigs_file,string *genome_dir);
string *Get_contig_file(string *contig_file,string *genome_dir);
TRIE *Get_file_seq_name(string *genome_dir);
string *Get_species_dir_from_supported_file(string *species_dir,char *species,char *assembly,char *release,char *supported_file);
string *Get_data_dir(string *data_dir);
string *Get_supported_file(string *supported_file);
string *Get_genomes_dir(string *genomes_dir);
string *Get_species_dir(string *species_dir,char *species, char *assembly, char *release, char *species_suffix);
string *Get_genome_dir(string *genome_dir,char *species, char *assembly, char *release, char *species_suffix);
string *Get_variation_dir(string *var_dir,char *species, char *assembly, char *release, char *species_suffix);
char *Get_sequence(char *sequence_file);
int switch_strand(char *sequence);

memstd *RsatMemTracker = NULL;

string *RSAT         = NULL;
string *HTML         = NULL;
string *LOGS         = NULL;
string *WWW_TMP      = NULL;
string *counter_file = NULL;

string *BIN           = NULL;
string *LIB           = NULL;
string *PYTHON        = NULL;
string *SCRIPTS       = NULL;

string *date                   = NULL;
string *HOSTNAME               = NULL;

unsigned int *PID              = NULL;
string *LOGIN                  = NULL;

string *CMD                    = NULL;
string *PROGRAM                = NULL;

string *log_file               = NULL;
string *exec_time_log_file     = NULL;
string *start_time_log_file    = NULL;
string *web_attacks_log_file   = NULL;
string *denied_access_log_file = NULL;

time_t start_time;
int verbose = 2;

void help(){
  printf(
  "\n"
  "NAME\n"
  "    variation-scan\n"
  "\n"
  "VERSION\n"
  "    2.0.7\n"
  "\n"
  "DESCRIPTION\n"
  "    Scan variant sequences with position specific scoring matrices (PSSM)\n"
  "    and report variations that affect the binding score, in order to predict\n"
  "    regulatory variants.\n"
  "\n"
  "AUTHORS\n"
  "    Alejandra Medina Rivera <amedina@liigh.unam.mx>\n"
  "    Walter Santana Garcia <wsantana@lcg.unam.mx>\n"
  "    Jacques van Helden <Jacques.van-Helden@univ-amu.fr>\n"
  "\n"
  "CATEGORY\n"
  "    util\n"
  "\n"
  "USAGE\n"
  "     variation-scan [-i sequence_file] -m matrix_file -bg backgournd_file [-calc_distrib] [-o outputfile] [-v #] [...]\n"
  "\n"
  "  Example\n"
  "INPUT FORMAT\n"
  "  Sequence file\n"
  "    variation-scan takes as input a variation file in the format produced by\n"
  "    retrieve-variation-seq. For details about this format, see\n"
  "    retrieve-variation-seq output format.\n"
  "\n"
  "  Matrix file\n"
  "    A list of matrix in transfanc/tab format\n"
  "\n"
  "  Background file\n"
  "    oligo-analysis format\n"
  "\n"
  "OUTPUT FORMAT\n"
  "    A tab delimited file with the following column content.\n"
  "\n"
  "    1. matrix_ac\n"
  "       Name of the matrix (generally the transcription factor name)\n"
  "\n"
  "    2. var_id\n"
  "       ID of the variation\n"
  "\n"
  "    3. var_class\n"
  "       Variation type, according to SNP Ontology (SO) nomenclature\n"
  "\n"
  "    4. var_coord\n"
  "       Coordinates of the variation\n"
  "\n"
  "    5. best_w\n"
  "       Best weigth for the putative site\n"
  "\n"
  "    6. worst_w\n"
  "       Worst weigth for the putative site\n"
  "\n"
  "    7. w_diff\n"
  "       Difference between best and worst weigth\n"
  "\n"
  "    8. best_pval\n"
  "       P_value of the best putative site\n"
  "\n"
  "    9. worst_pval\n"
  "       P_value of the worst putative site\n"
  "\n"
  "   10.pval_ratio\n"
  "       Ratio between worst and best pval ( pval_ratio = worst_pval/best_pval )\n"
  "\n"
  "   11. best_allele\n"
  "       Allele in the best putative site\n"
  "\n"
  "   12. worst_allele\n"
  "       Allele in the worst putative site\n"
  "\n"
  "   13. best_offest\n"
  "       Offset of the best putative site\n"
  "\n"
  "   14. worst_offset\n"
  "       Offset of the worst putative site\n"
  "\n"
  "   15. min_offset_diff\n"
  "       Difference minimal between best and worst putative site\n"
  "\n"
  "   16. best_strand\n"
  "       Strand of the best putative site\n"
  "\n"
  "   17. worst_strand\n"
  "       Strand of the worst putative site\n"
  "\n"
  "   18. str_change\n"
  "       Indicate if strand have change between the offset of min_offset_diff\n"
  "\n"
  "   19. best_seq\n"
  "       Sequence of the worst putative site\n"
  "\n"
  "   20. worst_seq\n"
  "       Sequence of the worst putative site\n"
  "\n"
/*  "   21. is_ref_better\n"
  "       Flag if the reference allele is the allele with the best score\n"
  "\n"*/
  "   21. minor_alle_freq\n"
  "       Minor allele frequency\n"
  "\n"
  "SEE ALSO\n"
  "  download-ensembl-genome\n"
  "    retrieve-variation-seq uses the sequences downloaded from Ensembl using\n"
  "    the tool download-ensembl-genome.\n"
  "\n"
  "  download-ensembl-variations\n"
  "    retrieve-variation-seq uses variation coordinates downloaded from\n"
  "    Ensembl using the tool download-ensembl-variations.\n"
  "\n"
  "  variation-scan\n"
  "    Scan variation sequences with one or several position-specific scoring\n"
  "    matrices.\n"
  "\n"
  "WISH LIST\n"
  "OPTIONS\n"
  "    -v #\n"
  "        Level of verbosity (detail in the warning messages during execution)\n"
  "\n"
  "    -h  Display full help message\n"
  "\n"
  "    -help\n"
  "        Same as -h\n"
  "\n"
  "    -i inputfile\n"
  "        Variation file RSAT format, i.e. varBed. For more details, see\n"
  "        <retrieve-variation-seq>\n"
  "\n"
  "    -m matrixfile\n"
  "        Matrix file transfac/tab format\n"
  "\n"
  "    -top_matrices #\n"
  "        Only work with the # top matrix\n"
  "\n"
  "    -bg Background file\n"
  "         oligo-analysis file format\n"
  "\n"
  "    -m_format matrix_format\n"
  "        Matrix file format\n"
  "\n"
  "        Supported formats: alignace, assembly, cb, clustal, cluster-buster,\n"
  "                            consensus, sequences, feature, footprintdb, gibbs,\n"
  "                            infogibbs, info-gibbs, jaspar, homer, mscan, meme,\n"
  "                            meme_block, motifsampler, stamp, stamp-transfac, tab,\n"
  "                            tf, transfac, cis-bp, uniprobe, yeastract and encode.\n"
  "\n"
  "    -mml #\n"
  "        Length of the longest Matrix, this values has to be consistent with\n"
  "        the one used io for retrieving the variant sequences (see\n"
  "        <retrieve-variation-seq>).\n"
  "\n"
/*  "    -top_matrix #\n"
  "        Only work with the # top matrix\n"
  "\n"*/
  "    -top_variation #\n"
  "        Only work with the # top variation\n"
  "\n"
  "    -lth [score | w_diff | pval_ratio] #\n"
  "        Only return variations with type score > #\n"
  "\n"
  "    -uth [pval] #\n"
/*  "    -html #\n"
  "        Convert the tab-delimited file into an HTML file, which facilitates\n"
  "        the inspection of the results with a Web browser. The HTML file has\n"
  "        the same name as the output file, but the extension (.tab, .txt) is\n"
  "        replaced by the .html extension\n"*/
  "\n"
  "    -calc_distrib\n"
  "        Calculate and save distribution of matrices\n"
  "\n"
  "    -distrib_dir path/to/dir\n"
  "        Directory to store the distribution files. Mandatory if\n"
  "        -calc_distrib is being used.\n"
  "\n"
/*  "    -distrib_list #\n"
  "        Name of the file containing the list of matrix distrib file name\n"
  "\n"
  "        NOTE: This file must be in the same directory as the distrib file\n"
  "\n"
  "    -only_biggest\n"
  "        Only return the biggest difference of score between two alleles of a\n"
  "        variation regarthless of the window, this option is usefull for\n"
  "        insertions and deletions\n"
  "\n"*/
  "    -o outputfile\n"
  "        The output file is in fasta format.\n"
  "\n"
  "        If no output file is specified, the standard output is used. This\n"
  "        allows to use the command within a pipe.\n"
  "\n"
  );
  rlist(RsatMemTracker);
  exit(0);
}


int main(int argc, char *argv[]){
  /////////////////////////////////////////////////
  // Declare variables
  /////////////////////////////////////////////////
  //Files
  FILE *fout            = stdout;
  FILE *fin             = stdin;
  //FILE *fh_popen        = NULL;
  //Numeric
  //int no_offset      =  1;
  int nb_matrix      =  0;
  unsigned long int nb_variation   =  0;
  unsigned long int nb_seq         =  0;
  char *end_nb = NULL;
  //int nb_var         =  0;
  //int output_lines   =  0;
  //int flank_len      = 29; // Default -mml value
  //float pval_limit   =  1;
  int i              =  0;
  //int retstatus      =  0;

  //Lower thresholds
  //float score      = 1;
  //float w_diff     = 1;
  //float pval_ratio = 1;
  //Upper thresholds
  //float pval       = 1;
  //NOTE. I set to NAN due to lack of information
  //for lowest possible score.
  threshold cutoff  = {.upper.pval   = 0.0,
                       .lower.score  = NAN,
                       .lower.wdiff  = 0.0,
                       .lower.pratio = 0.0 };

  //Structs
  //struct stat dir_exists;
  TRIE *trieMatrixFmts      = NULL;
  string *program_full_path = NULL;
  string *mscanquick_cmd    = NULL;
  string *deletetmps_cmd    = NULL;
  string *out_dir           = NULL;
  string *varsequence       = NULL;
  string *variationscan_file = NULL;

  /////////////////////////////////
  // Input variables
  /////////////////////////////////
  //Files
  char *input                  = NULL;
  char *output                 = NULL;
  stringlist *matrix_files     = NULL;
  string *distrib_dir          = NULL;
  char *distrib_list           = NULL;
  char *bg                     = NULL;
  string *matrix_format        = NULL;
  char *type                   = NULL;

  //Flags
  int only_biggest     =    0;
  int calc_distrib     =    0;
  int debug            =    0;

  //Numeric
  //int top_matrix       =   -1;
  unsigned long int top_variation    =    0;
  int top_matrices     =    0;
  int mml              =   29;

  //Thresholds
  //float lth            =    0;
  //float uth            =    0;

  /////////////////////////////////////////////////
  // Read Arguments
  /////////////////////////////////////////////////

  //Initialize RsatMemTracker for memory tracking
  //it also add contents to RSAT and Initialize CMD,
  //PROGRAM, and PID.
  InitMemLists();

  //Start Creating CMD string.
  strcopy(CMD,"variation-scan");
  strcopy(PROGRAM,"variation-scan");

  //Create list of matrix files
  matrix_files  = strlistnewToList(&RsatMemTracker);
  matrix_format = strnewToList(&RsatMemTracker);
  distrib_dir   = strnewToList(&RsatMemTracker);

  strcopy(distrib_dir,"");
  strcopy(matrix_format,"transfac");
  strcopy(matrix_files->element,"");

  if (argc <= 1) help();

  for (i = 1; i < argc; i++) {
    if ( strcmp(argv[i],"-h") == 0 || strcmp(argv[i],"--help") == 0) help();
    else if ( strcmp(argv[i],"-i") == 0 && CheckValOpt(argv+i)) {
      strccat(CMD," %s %s",argv[i],argv[i+1]);
      input = argv[++i];
    } else if ( strcmp(argv[i],"-o") == 0 && CheckValOpt(argv+i)) {
      strccat(CMD," %s %s",argv[i],argv[i+1]);
      output = argv[++i];
    } else if ( strcmp(argv[i],"-m") == 0 && CheckValOpt(argv+i) ) {
      strccat(CMD," %s %s",argv[i],argv[i+1]);
      argToList(matrix_files, argv[i+1]);
      i++;
    } else if ( strcmp(argv[i],"-top_matrices") == 0 && CheckValOpt(argv+i) ) {
      strccat(CMD," %s %s",argv[i],argv[i+1]);
    	top_matrices = atoi(argv[++i]);
    } else if ( strcmp(argv[i],"-bg") == 0 && CheckValOpt(argv+i) ) {
      strccat(CMD," %s %s",argv[i],argv[i+1]);
      bg = argv[++i];
    } else if ( strcmp(argv[i],"-m_format") == 0 && CheckValOpt(argv+i) ) {
      strccat(CMD," %s %s",argv[i],argv[i+1]);
      strcopy(matrix_format,argv[++i]);
      //matrix_format = argv[++i];
    } else if (strcmp(argv[i],"-mml") == 0 && CheckValOpt(argv+i) ) {
      strccat(CMD," %s %s",argv[i],argv[i+1]);
      mml = atoi(argv[++i]);
    } else if ( strcmp(argv[i],"-top_variation") == 0 && CheckValOpt(argv+i) ) {
      strccat(CMD," %s %s %s",argv[i],argv[i+1]);
      top_variation = strtoul(argv[++i], &end_nb, 10);
      //top_variation = atoi(argv[++i]);
    } else if ( strcmp(argv[i],"-lth") == 0 && CheckValOpt(argv+i) ) {
      strccat(CMD," %s %s %s",argv[i],argv[i+1],argv[i+2]);
      type =  argv[++i];
      sethreshold(&cutoff,type,argv[++i]);
      //lth  =  atof(argv[++i]);
    } else if ( strcmp(argv[i],"-uth") == 0 && CheckValOpt(argv+i) ) {
      strccat(CMD," %s %s %s",argv[i],argv[i+1],argv[i+2]);
      type =  argv[++i];
      sethreshold(&cutoff,type,argv[++i]);
      //uth  =  atof(argv[++i]);
    } else if ( strcmp(argv[i],"-calc_distrib") == 0 ) {
      strccat(CMD," %s",argv[i]);
      calc_distrib = 1;
    } else if ( strcmp(argv[i],"-distrib_dir") == 0 && CheckValOpt(argv+i) ) {
      strccat(CMD," %s %s %s",argv[i],argv[i+1],argv[i+2]);
      strcopy(distrib_dir, argv[++i]);
    } else if ( strcmp(argv[i],"-distrib_list") == 0 && CheckValOpt(argv+i) ) {
      strccat(CMD," %s %s",argv[i],argv[i+1]);
      distrib_list  = argv[++i];
    } else if ( strcmp(argv[i],"-only_biggest") == 0 && CheckValOpt(argv+i) ) {
      strccat(CMD," %s",argv[i]);
      only_biggest = 1;
    } else if ( strcmp(argv[i],"-v") == 0 && CheckValOpt(argv+i) ) {
      strccat(CMD," %s %s",argv[i],argv[i+1]);
      verbose = atoi(argv[++i]);
    } else if ( strcmp(argv[i],"-debug") == 0 ) {
      debug   = 1;
    } else {
      RsatFatalError("Invalid option", argv[i], NULL);
    }

  }


  /////////////////////////////////////////////////
  //Initialize script
  /////////////////////////////////////////////////
  start_time = StartScript(PROGRAM->buffer,CMD->buffer);

  //Allocate memory for strings
  program_full_path = strnewToList(&RsatMemTracker);
  mscanquick_cmd    = strnewToList(&RsatMemTracker);
  out_dir           = strnewToList(&RsatMemTracker);
  deletetmps_cmd    = strnewToList(&RsatMemTracker);
  varsequence       = strnewToList(&RsatMemTracker);
  variationscan_file = strnewToList(&RsatMemTracker);


  //Get supported matrix formats as a TRIE
  trieMatrixFmts = ListInputMatrixFormats();

  /////////////////////////////////////////////////
  // Validate Arguments
  /////////////////////////////////////////////////
  if (input){
    strcopy(varsequence, input);
    fin = OpenInputFile(fin,input);
  } else {
    strfmt(varsequence, "");
  }
  if (output){
    strcopy(variationscan_file, output);
    fout = OpenOutputFile(fout,output);
  } else {
    strfmt(variationscan_file, "");
  }
  if(!bg)
    RsatFatalError("No Background file was supplied. Use -bg option to specify it", NULL);
  if( strcmp(matrix_files->element->buffer, "") == 0  )
    RsatFatalError("No Matrix file was supplied. Use -m option to specify it", NULL);
  if( TrieSearch(trieMatrixFmts,matrix_format->buffer) == NULL )
    RsatFatalError("Invalid input format for matrix. Please see list of supported formats with --help", NULL);
  if(  distrib_list &&  access(distrib_list,F_OK) == -1)
    RsatFatalError("Distribution list", distrib_list, "does not exists", NULL);
  if(calc_distrib && (strcmp(distrib_dir->buffer,"") == 0))
    RsatFatalError("-calc_distrib option must be used together with -distrib_dir. Please supply a directory with it", NULL);

  /*if( !type ){
    if( lth != 0){
      if (strcmp(type,"score") == 0) {
        score      = lth;
      } else if (strcmp(type,"w_diff") == 0) {
        w_diff     = lth;
      } else if (strcmp(type,"pval_ratio") == 0) {
        pval_ratio = lth;
      } else {
        RsatFatalError("This is not a valid type -lth option",type,".Please recheck.",NULL);
      }
    } else if (uth != 0){
      if (strcmp(type,"score") == 0) {
        pval = uth;
      } else {
        RsatFatalError("This is not a valid type -uth option",type,".Please recheck.",NULL);

      }
    }
  }*/

  /////////////////////////////////////////////////
  // Check if matrix-scan-quick is available
  /////////////////////////////////////////////////
  GetProgramPath(program_full_path ,"matrix-scan-quick", 0, NULL);
  //strfmt(mscanquick_cmd, "%s -h", program_full_path->buffer);
  if(strcmp(program_full_path->buffer, "") == 0) {
    RsatFatalError("Unable to find matrix-scan-quick. It is not installed. Please recheck",NULL);
    //if ( (fh_popen = popen(mscanquick_cmd->buffer,"r")) == NULL ) RsatFatalError("Unable to popen in main()",NULL);
    //retstatus = pclose(fh_popen);
    //RsatFatalError("Unable to find matrix-scan-quick. It is not installed. Please recheck",NULL);
  }

  /////////////////////////////////////////////////
  // Calculate distribution of scores
  /////////////////////////////////////////////////
  RsatInfo("Calculating distributions", NULL);

  //Create both tmp directories
  make_temp_file(out_dir,"","variation_scan",1,1,0);

  if ( strcmp(distrib_dir->buffer,"") == 0 )
    strcopy(distrib_dir, out_dir->buffer);

  if ( CheckOutDir(out_dir,0,0755) == 0 )
    RsatFatalError("Unable to create directory",out_dir->buffer,"in main()",NULL);

  if ( CheckOutDir(distrib_dir,0,0755) == 0 )
    RsatFatalError("Unable to create directory",distrib_dir->buffer,"in main()",NULL);

  //Declare variables
  FILE *fh_matrix_list            =   NULL;
  FILE *fh_distrib_list           =   NULL;
  char   **token                  =   NULL;
  char   *index                   =   NULL;
  string *dir_and_file[2]         = {NULL};
  string *matrixprefix[2]         = {NULL};
  string *line                    =   NULL;
  string *ls_cmd                  =   NULL;
  string *bgprefix                =   NULL;
  string *out_distrib_list        =   NULL;
  string *out_distrib_list_sort   =   NULL;
  string *input_matrix_list       =   NULL;
  string *matrix_distrib_cmd      =   NULL;
  string *sort_matrix_cmd         =   NULL;
  string *convert_matrix_cmd      =   NULL;
  string *convert_bg_cmd          =   NULL;
  string *bg_inclusive            =   NULL;
  stringlist *matrixID            =   NULL;
  stringlist *matrixfileprefix    =   NULL;


  //Allocate memory for variables
  dir_and_file[0]         = strnewToList(&RsatMemTracker);
  dir_and_file[1]         = strnewToList(&RsatMemTracker);
  matrixprefix[0]         = strnewToList(&RsatMemTracker);
  matrixprefix[1]         = strnewToList(&RsatMemTracker);
  bgprefix                = strnewToList(&RsatMemTracker);
  ls_cmd                  = strnewToList(&RsatMemTracker);
  convert_matrix_cmd      = strnewToList(&RsatMemTracker);
  sort_matrix_cmd         = strnewToList(&RsatMemTracker);
  convert_bg_cmd          = strnewToList(&RsatMemTracker);
  out_distrib_list        = strnewToList(&RsatMemTracker);
  out_distrib_list_sort   = strnewToList(&RsatMemTracker);
  input_matrix_list       = strnewToList(&RsatMemTracker);
  matrix_distrib_cmd      = strnewToList(&RsatMemTracker);
  bg_inclusive            = strnewToList(&RsatMemTracker);
  line                    = strnewToList(&RsatMemTracker);
  token                   = getokens(4);


  SplitFileName(dir_and_file,bg);
  index = strrchr(dir_and_file[1]->buffer, '.');
  if (index) {
    *index = '\0';
    dir_and_file[1]->size = strlen(dir_and_file[1]->buffer) + 1;
  }
  //Copy file bg prefix, create file name for distrib list and distrib list sort
  strcopy(bgprefix , dir_and_file[1]->buffer);
  strfmt(out_distrib_list_sort, "%s/%s_list.sort.tab", distrib_dir->buffer, dir_and_file[1]->buffer);
  strccat(dir_and_file[1], "_list.tab");


  //Prepare output distribution list, open it and print header to it
  strfmt(out_distrib_list , "%s/%s", distrib_dir->buffer, dir_and_file[1]->buffer );
  fh_distrib_list = OpenOutputFile(fh_distrib_list, out_distrib_list->buffer  );
  fprintf(fh_distrib_list, "#MATRIX_PREFIX\tMATRIX_ID\tMATRIX_FILE\tDISTRIB_FILE\tDB\tBG_PREFIX\n" );

  //Prepare stringlists: matrixID and matrixfileprefix
  matrixID          = strlistnewToList(&RsatMemTracker);
  matrixfileprefix  = strlistnewToList(&RsatMemTracker);
  strcopy(matrixID->element,"");
  strcopy(matrixfileprefix->element,"");

  for (stringlist *curr_matrix = matrix_files; curr_matrix != NULL ; curr_matrix = curr_matrix->next) {

    //Remove directory from file path and delete '.' extension
    SplitFileName(matrixprefix,curr_matrix->element->buffer);
    index = strrchr(matrixprefix[1]->buffer, '.');
    if (index) {
      *index = '\0';
      matrixprefix[1]->size = strlen(matrixprefix[1]->buffer) + 1;
    }
    //printf("This is matrixprefix[1] %s\n", matrixprefix[1]->buffer);

    //Convert-matrixes; Split matrix files into individual tab-formatted files
    strfmt(convert_matrix_cmd,"%s/perl-scripts/convert-matrix -v 1 -from %s -to tab -split -i %s -o %s/%s",
    RSAT->buffer, matrix_format->buffer, curr_matrix->element->buffer, out_dir->buffer, matrixprefix[1]->buffer);
    if(verbose >= 6) RsatInfo("\nConverting matrices from", matrix_format->buffer, "to tab format", convert_matrix_cmd->buffer, NULL);
    doit(convert_matrix_cmd->buffer,0,1,0,0,NULL,NULL,NULL,NULL);

    //Test if convert-matrix created successfully the matrix files in tab format
    if ( CheckConvertMatrixOutput(out_dir->buffer, matrixprefix[1]->buffer) == 0 ){
      RsatFatalError("Matrix file could not be converted to tab format. Please check the input matrix format is consistent with the -m_format option.", NULL);
    }

    //Prepare current input matrix list for reading
    strfmt(input_matrix_list, "%s/%s_matrix_list.tab", out_dir->buffer, matrixprefix[1]->buffer);
    //Open file
    fh_matrix_list  = OpenInputFile(fh_matrix_list,   input_matrix_list->buffer );

    //Print header to output
    i = 0;
    token[0] = line->buffer;

    while ( fread( (line->buffer + line->size),1,1,fh_matrix_list) == 1 ) {
      //If '\t' is found assign next char address as the next token
      if (line->buffer[line->size] == '\t') {
        token[++i] = line->buffer + line->size + 1;
        line->buffer[line->size] = '\0';
      }

      //If end of line is found,proceed to process line
      if (line->buffer[line->size] == '\n') {

        //Reinitialize counters
        line->buffer[line->size] = '\0';
        line->size = 0;
        i = 0;

        //Filters for comments and more
        if(token[0][0] == '#'){initokadd(line,token,4);continue;}
        if(token[0][0] == ';'){initokadd(line,token,4);continue;}
        if(token[0][0] == ' '){initokadd(line,token,4);continue;}
        if(token[0][0] == '\t'){initokadd(line,token,4);continue;}
        if(token[0][0] == '\0'){initokadd(line,token,4);continue;}
        //Check if nb of top matrixes has been reached
        if(top_matrices && !(nb_matrix < top_matrices)) break;
        nb_matrix++;

        //WSG I removed this piece of code, no need. NOTE Look for the variables rm
        //Add content to list of matrixes
        //argToList(matrixID, token[1]);
        //argToList(matrixfileprefix, matrixprefix[1]->buffer);

        //Print content to distribution list
        fprintf(fh_distrib_list, "%s\t%s\t%s/%s_%s.tab\t%s_%s_%s.distrib\t.\t%s\n",
                matrixprefix[1]->buffer,
                token[1],
                out_dir->buffer, matrixprefix[1]->buffer, token[1],
                bgprefix->buffer, matrixprefix[1]->buffer, token[1],
                bgprefix->buffer );

        //If successful continue to the next start of line for reading
        initokadd(line,token,4);
        continue;
      }

      //Resize line string if limit has reached
      limlinetok(line,token,4);

      //Keep track of size count for each char read
      line->size++;

    }

    //Close input matrix list fh
    fclose(fh_matrix_list);
  }

  //Close output distrib list fh
  fclose(fh_distrib_list);

  //Remove tmp variables
  RsatMemTracker = relem((void*)token, RsatMemTracker);

  if(verbose >= 12) RsatInfo("Distribution list :", out_distrib_list->buffer, NULL);

  /////////////////////////////////////////////////////////////////////////
  // Sort matrix distrib list file by matrix length
  /////////////////////////////////////////////////////////////////////////
  strfmt(sort_matrix_cmd,"\npython %s/contrib/variation-scan/sort-matrix.py -i %s -o %s",
          RSAT->buffer, out_distrib_list->buffer, out_distrib_list_sort->buffer);
  if(verbose >= 6)  RsatInfo("\nSorting matrix list by size :", sort_matrix_cmd->buffer,"\n", NULL);
  doit(sort_matrix_cmd->buffer,0,1,0,0,NULL,NULL,NULL,NULL);
  if(verbose >= 12) RsatInfo("Sorted distribution list :", out_distrib_list_sort->buffer, NULL);


  /////////////////////////////////////////////////////////////////////////
  // Convert background model to inclusive format
  /////////////////////////////////////////////////////////////////////////

  strfmt(bg_inclusive, "%s/vscan_bg_inclusive_file.inclusive", out_dir->buffer);
  strfmt(convert_bg_cmd, "%s/convert-background-model -i %s -from oligos -to inclusive "
                         "-o %s -bg_pseudo 0.01", SCRIPTS->buffer, bg, bg_inclusive->buffer);
  if(verbose >= 6) RsatInfo("\nConverting bg model to inclusive", convert_bg_cmd->buffer, NULL);
  doit(convert_bg_cmd->buffer,0,1,0,0,NULL,NULL,NULL,NULL);

  /////////////////////////////////////////////////////////////////////////
  // Detect number of fields for haplotype or single variant analysis
  /////////////////////////////////////////////////////////////////////////
  //Declare variables
  FILE *fh_varsequence    = NULL;
  //FILE *fh_fasta_sequence = NULL;

  //string *fasta_sequence  = NULL;
  int num_tokens          =    0;
  int haplotype           =    0;

  //Allocate memory for variables
  //fasta_sequence = strnewToList(&RsatMemTracker);
  token          = getokens(11);

  //Prepare variables for buffer reading
  i = 0;
  line->size = 0;
  token[0] = line->buffer;

  // Open stream if needed
  if( strcmp(varsequence->buffer, "") != 0){
    fh_varsequence    = OpenInputFile(fh_varsequence, varsequence->buffer );
  } else {
    fh_varsequence = fin;
  }
  //fh_fasta_sequence = OpenOutputFile(fh_fasta_sequence, fasta_sequence->buffer );

  //if(verbose >= 6) RsatInfo("Creating fasta sequence file", fasta_sequence->buffer, NULL);

  while ( fread( (line->buffer + line->size),1,1,fh_varsequence) == 1 ) {
    //If '\t' is found assign next char address as the next token
    if (line->buffer[line->size] == '\t') {
      token[++i] = line->buffer + line->size + 1;
      line->buffer[line->size] = '\0';
    }
    //If end of line is found,proceed to process line
    if (line->buffer[line->size] == '\n') {
      //Save total number of tokens
      num_tokens = i;
      //Reinitialize counters
      line->buffer[line->size] = '\0';
      line->size = 0;
      i = 0;

      //Filters for comments and more
      if(token[0][0] == ';') {initokadd(line,token,11);continue;}
      if(token[0][0] == ' ') {initokadd(line,token,11);continue;}
      if(token[0][0] == '\t'){initokadd(line,token,11);continue;}
      if(token[0][0] == '\0'){initokadd(line,token,11);continue;}
      if(token[0][0] == '#') {
        break;
      }
    }
    //Resize line string if limit has reached
    limlinetok(line,token,11);

    //Keep track of size count for each char read
    line->size++;
  }

  /////////////////////////////////////////////////
  // Print Header
  /////////////////////////////////////////////////
  //printf("\n\nHERE1!!!!!!!!!!!!!!!!!!! \n\n" );
  //fclose(fout);
  //printf("\n\nHERE2!!!!!!!!!!!!!!!!!!! \n\n" );

  //Test for number of fields i.e. which kind of sequence file is it reading
  //(Haplotypes or single variants) and create fasta.
  //Add +1 due token index starting on 0.
  num_tokens += 1;
  if(num_tokens == 11){
    haplotype = 1;
    if(verbose >= 6) RsatInfo("Detected haplotype variants.",NULL);
    printHeaderHaplotypes(variationscan_file,   input, output, CMD);

    //CreateFastaFromHaplosequences(line, token, &nb_variation, top_variation, &nb_seq, fh_varsequence, fh_fasta_sequence);
  } else if (num_tokens == 10){
    haplotype = 0;
    if(verbose >= 6) RsatInfo("Detected single variants.",NULL);
    printHeaderSingleVariants(variationscan_file, input, output, CMD);
    //CreateFastaFromSingleVariants(line, token, &nb_variation, top_variation, &nb_seq, fh_varsequence, fh_fasta_sequence);
  } else {
    RsatFatalError("This file does not contain the correct number of fields",NULL);
  }

  //Close input variation sequence and output fasta fh
  if( strcmp(varsequence->buffer, "") != 0) {
    fclose(fh_varsequence);
  }

  //fclose(fh_fasta_sequence);

  //Remove tmp variables
  RsatMemTracker = relem((void*)token, RsatMemTracker);



  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  ////                       CODE-SECTION-IN-WORK                      ////

  //Declare variables
  FILE *fh_distrib_list_sort  =   NULL;
  string *mscanquick_file     =   NULL;
  char *matrix_prefix         =   NULL;
  char *matrix_id             =   NULL;
  char *matrix_file           =   NULL;
  int   matrix_size           =      0;
  char *distrib_file          =   NULL;
  char *db                    =   NULL;
  char *bg_prefix             =   NULL;

  int curr_matrix             =      0;
  int prev_matrix_size        =      0;

  string *prev_fasta  = NULL;
  string *new_fasta   = NULL;
  //Allocate memory for variables
  token           = getokens(7);
  mscanquick_file = strnewToList(&RsatMemTracker);
  prev_fasta      = strnewToList(&RsatMemTracker);
  new_fasta       = strnewToList(&RsatMemTracker);

  //Prepare variables for buffer reading
  i = 0;
  line->size = 0;
  token[0] = line->buffer;

  //Open file containing all matrices files sorted by decreasing length
  //and the distribution name for each matrix. It will be used to generate
  //fasta sequences each time shorter for matrix-scan-quick.
  fh_distrib_list_sort = OpenInputFile(fh_distrib_list_sort, out_distrib_list_sort->buffer);

  //////////////////////////////////////////////////////////////////////////////////
  // Get distributions for each matrix,create fastas and execute matrix-scan-quick
  /////////////////////////////////////////////////////////////////////////////////
  // #MATRIX_PREFIX\tMATRIX_ID\tMATRIX_FILE\tMATRIX_SIZE\tDISTRIB_FILE\tDB\tBG_PREFIX
  while ( fread( (line->buffer + line->size),1,1,fh_distrib_list_sort) == 1 ) {
    //If '\t' is found assign next char address as the next token
    if (line->buffer[line->size] == '\t') {
      token[++i] = line->buffer + line->size + 1;
      line->buffer[line->size] = '\0';
    }
    //If end of line is found,proceed to process line
    if (line->buffer[line->size] == '\n') {
      //Save total number of tokens
      //num_tokens = i;
      //Reinitialize counters
      line->buffer[line->size] = '\0';
      line->size = 0;
      i = 0;

      //Skip commented line
      if(token[0][0] == '#') {initokadd(line,token,7);continue;}

      //Assign tokens to new variable names for an easier code reading
      matrix_prefix         =   token[0];
      matrix_id             =   token[1];
      matrix_file           =   token[2];
      matrix_size           =   atoi(token[3]);
      distrib_file          =   token[4];
      db                    =   token[5];
      bg_prefix             =   token[6];

      //Counter for number of current matrix
      curr_matrix++;

      //////////////////////////////////////////////////////////////////////////
      ///////////////////    Create fasta sequences     ////////////////////////
      /////////////////////////////////////////////////////////////////////////

      //Check if the matrix size of the previous matrix is the same as the current
      //The same fasta sequence file will be analyzed if this is the case.
      if (matrix_size != prev_matrix_size) {
        //For the first matrix fastas are created directly from varSeq files
        if(curr_matrix == 1) {
          //printf("HELLO\n\n\n");
          //Open input varsequence file and ouput fasta file
          // Verbose of fasta sequence creation
          strfmt(prev_fasta, "%s/%d_vscan_fasta.fa",out_dir->buffer,matrix_size);
          if(verbose >= 6) RsatInfo("Creating fasta sequence file ", prev_fasta->buffer, NULL);
          //Create fasta files from haplotype varSeq
          if (haplotype){
            CreateFastaFromVarseqHaplotypes(varsequence,prev_fasta , matrix_size, &nb_variation, &top_variation, &nb_seq);
          //Create fasta files from single variants varSeq
          } else {
            CreateFastaFromVarseqVariants(varsequence, prev_fasta , matrix_size, &nb_variation, &top_variation, &nb_seq);
          }

        //For other matrices different from the 1st one fastas are created
        //directly from the previous fasta sequences. By decresing the length
        //of flanking sequences
        } else {
          //Reinitialize counters
          nb_seq = 0;
          nb_variation = 0;

          // Verbose of fasta sequence creation
          strfmt(new_fasta, "%s/%d_vscan_fasta.fa",out_dir->buffer,matrix_size);
          if(verbose >= 6) RsatInfo("Creating fasta sequence file ", prev_fasta->buffer, NULL);
          //Create fasta files from haplotype varSeq
          if (haplotype){
            CreateFastaFromVarseqHaplotypes(varsequence,new_fasta , matrix_size, &nb_variation, &top_variation, &nb_seq);
            //CreateFastaFromFastaHaplotypes(prev_fasta,new_fasta , matrix_size, &nb_variation, &top_variation, &nb_seq);

          //Create fasta files from single variants varSeq
          } else {
            CreateFastaFromFastaVariants(prev_fasta, new_fasta , matrix_size, &nb_variation, &top_variation, &nb_seq);

          }

          ///////////////////////////////////////////////////////////////
          //Remove tmp prev_fasta sequence
          if (!debug) {
            strfmt(deletetmps_cmd,"rm %s",prev_fasta->buffer);
            doit(deletetmps_cmd->buffer,0,0,0,0,NULL,NULL,NULL,NULL);
          }
          //UPDATE new fasta to be previous fasta for future loops
          strcopy(prev_fasta, new_fasta->buffer);
        }

      }

      ///////////////////////////////////////////////////////////////////////////
      ///////    Create matrix distributions (matrix-distrib)       ////////////
      //////////////////////////////////////////////////////////////////////////

      //Prepare commandmd for matrix-distrib
      strfmt(matrix_distrib_cmd,"%s/matrix-distrib -m %s -matrix_format tab -decimals 1 "
      "-bgfile %s -bg_pseudo 0.01 -bg_format oligos -pseudo 1 -o %s/%s",
      SCRIPTS->buffer,
      matrix_file,
      bg,
      distrib_dir->buffer, distrib_file);

      //Execute matrix-distrib command
      if(verbose >= 6) RsatInfo("\nCalculating p-val distribution for matrix",
                                 matrix_id, matrix_distrib_cmd->buffer,"\n", NULL);
      doit(matrix_distrib_cmd->buffer,0,1,0,0,NULL,NULL,NULL,NULL);

      ///////////////////////////////////////////////////////////////////////////
      ///////        Scan fasta sequences (matrix-scan-quick)       ////////////
      //////////////////////////////////////////////////////////////////////////

      ////////////////////////////////////////////
      //Prepare cmd for matrix-scan-quick
      strfmt(mscanquick_file,"%s/%s_%s_%s.mscan",
                              out_dir->buffer, bg_prefix, matrix_prefix, matrix_id);
      strfmt(mscanquick_cmd, "%s -i %s -m %s -pseudo 1 -decimals 1 -2str -origin start -bgfile %s "
            "-name %s -t 1 -distrib %s/%s -o %s",
            program_full_path->buffer,
            prev_fasta->buffer,
            matrix_file,
            bg_inclusive->buffer,
            matrix_id,
            distrib_dir->buffer, distrib_file,
            mscanquick_file->buffer);


      //Execute matrix-scan-quick
      if(verbose >= 6) RsatInfo("Running matrix-scan-quick", mscanquick_cmd->buffer,"\n", NULL);
      doit(mscanquick_cmd->buffer,0,1,0,0,NULL,NULL,NULL,NULL);


      ///////////////////////////////////////////////////////////////////////////
      ///////          Analyze scans from matrix-scan-quick)        ////////////
      //////////////////////////////////////////////////////////////////////////

      //Analyze scans for haplotypes
      if (haplotype){
        if(verbose >= 6) RsatInfo("Detected haplotype variants. Analyzing variations for matrix", matrix_id, NULL);
        ScanHaplosequences(mscanquick_file, variationscan_file, matrix_id, matrix_size, &cutoff);
      //Analyze scans for single variants
      } else {
        if(verbose >= 6) RsatInfo("Detected single variants. Analyzing variations for matrix", matrix_id, NULL);
        ScanSingleVariants(mscanquick_file, variationscan_file, matrix_id, matrix_size, &cutoff);
      }

      ///////////////////////////////////////////////////////////////////////////
      ////////            Remove tmp matrix-scan-quick files             ///////
      //////////////////////////////////////////////////////////////////////////
      if (!debug) {
        strfmt(deletetmps_cmd,"rm %s",mscanquick_file->buffer);
        doit(deletetmps_cmd->buffer,0,0,0,0,NULL,NULL,NULL,NULL);
      }

      ///////////////////////////////////////////////////////////////////////////
      ////////    Update previous matrix size value with the current      ///////
      ///////////////////////////////////////////////////////////////////////////
      prev_matrix_size = matrix_size;

      //If successful continue to the next start of line for reading
      initokadd(line,token,7);continue;
    }
    //Resize line string if limit has reached
    limlinetok(line,token,7); // QUESTION: SHOULD token here be == 7?

    //Keep track of size count for each char read
    line->size++;
  }



  ////                             E N D                               ////
  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////



  //Update execution log files
  ReportExecutionTime(start_time);
  //Free all objects
  rlist(RsatMemTracker);
  //Close file handlers
  fclose(fin);
  //fclose(fout);

  return 0;
}
/*Counts the number of files created after convert-matrix has been executed
  in variation-scan.  Returns the number of files by executing the command
  ls -lhtr out_dir/matrix_filname_*.tab | wc -l and popen */
int CheckConvertMatrixOutput(char *out_dir, char *matrix_filename){
  //Declare variables
  FILE   *fh_popen = NULL;
  string *cmd_countFiles = NULL;
  string *popen_out = NULL;

  int nb_files = 0;

  //Allocate memory
  cmd_countFiles   = strnewToList(&RsatMemTracker);
  popen_out        = strnewToList(&RsatMemTracker);

  //Create cmd
  strfmt(cmd_countFiles, "ls -lhtr %s/%s_*.tab | wc -l ",out_dir,matrix_filename);

  //Execute cmd in cmdline and write obtained FILENAME from mktemp to tmp_file
  if ( (fh_popen = popen(cmd_countFiles->buffer,"r")) == NULL ){
    RsatFatalError("Unable to popen in CheckConvertMatrixOutput()",NULL);
  }
  popen_out->size = 0;
  while ( (popen_out->buffer[popen_out->size] = fgetc(fh_popen)) != '\n' ) {
    //printf("This is the out of popen : %s \n",tmp_file->buffer );
    popen_out->size++;
    strlimt(popen_out);
  }

  //Add '\0'-end to string and +1 to size
  popen_out->buffer[popen_out->size] = '\0';
  popen_out->size++;
  //Close filehandler
  pclose(fh_popen);

  //Convert number of files string to char
  //Substract -1 because a log file
  //*_matrix_list.tab is always created
  nb_files = atoi(popen_out->buffer) - 1;

  if (verbose >= 12) RsatInfo("CheckConvertMatrixOutput() result", popen_out->buffer,NULL);

  //Remove tmp variables
  RsatMemTracker = relem( (void*)popen_out,RsatMemTracker );
  RsatMemTracker = relem( (void*)cmd_countFiles,RsatMemTracker );

  return nb_files;
}

void printHeaderSingleVariants(string *varscanFile, char*input, char *output, string *CMD){
  //Declare variables
  FILE *fh_varscan_output = stdout;

  //Open filehandler
  if(output) {
    fh_varscan_output = OpenOutputFile(fh_varscan_output, varscanFile->buffer);
  }

  //Print command-line and input/output files
  fprintf(fh_varscan_output,"; %s\n",CMD->buffer);

  //Test if input file IS stdin
  if(input) {
    fprintf(fh_varscan_output,"; Input file\n; \t%s\n",input);
  }

  //Test if output file IS stdout
  if(output) {
    fprintf(fh_varscan_output,"; Output file\n; \t%s\n",output);
  }

  //Print header
  fprintf(fh_varscan_output,
          ";\tcolumn content\n"
          ";\t1\tmatrix_ac\tAccession number of the positions-pecific scoring matrix\n"
          ";\t2\tvar_id\tID of the variation\n"
          ";\t3\tvar_class\tVariation type, according to SNP Ontology (SO) nomenclature\n"
          ";\t4\tvar_coord\tCoordinates of the variation\n"
          ";\t5\tbest_w\tBest weight for the putative site\n"
          ";\t6\tworst_w\tWorst weight for the putative site\n"
          ";\t7\tw_diff\tDifference between best and worst weight\n"
          ";\t8\tbest_pval\tP_value of the best putative site\n"
          ";\t9\tworst_pval\tP_value of the worst putative site\n"
          ";\t10\tpval_ratio\tRatio between worst and best pval ( pval_ratio = worst_pval/best_pval )\n"
          ";\t11\tbest_variant\tVariant in the best putative site\n"
          ";\t12\tworst_variant\tVariant in the worst putative site\n"
          ";\t13\tbest_offest\tOffset of the best putative site\n"
          ";\t14\tworst_offset\tOffset of the worst putative site\n"
          ";\t15\tmin_offset_diff\tDifference minimal between best and worst putative site\n"
          ";\t16\tbest_strand\tStrand of the best putative site\n"
          ";\t17\tworst_strand\tStrand of the worst putative site\n"
          ";\t18\tstr_change\tIndicate if strand has changed between the offset of min_offset_diff\n"
          ";\t19\tbest_seq\tSequence of the worst putative site\n"
          ";\t20\tworst_seq\tSequence of the worst putative site\n"
          ";\t21\tminor_alle_freq\tMinor allele frequency\n"
          "#ac_motif\tvar_id\tvar_class\tvar_coord\tbest_w\tworst_w\tw_diff\tbest_pval\tworst_pval\t"
          "pval_ratio\tbest_variant\tworst_variant\tbest_offset\tworst_offset\tmin_offset_diff\tbest_strand\t"
          "worst_strand\tstr_change\tbest_seq\tworst_seq\tminor_allele_freq\n");


  //Close fh
  if(output) {
    fclose(fh_varscan_output);
  }

  return;
}


void printHeaderHaplotypes(string *varscanFile, char*input, char *output, string *CMD){
  //Declare variables
  FILE *fh_varscan_output = stdout;

  //Open filehandler
  if(output) {
    fh_varscan_output = OpenOutputFile(fh_varscan_output, varscanFile->buffer);
  }

  //Print command-line and input/output files
  fprintf(fh_varscan_output,"; %s\n",CMD->buffer);

  //Test if input file IS stdin
  if(input) {
    fprintf(fh_varscan_output,"; Input file\n; \t%s\n",input);
  }

  //Test if output file IS stdout
  if(output) {
    fprintf(fh_varscan_output,"; Output file\n; \t%s\n",output);
  }

  //Print header
  fprintf(fh_varscan_output,
          ";\tcolumn content\n"
          ";\t1\tmatrix_ac\tAccession number of the positions-pecific scoring matrix\n"
          ";\t2\thap_id\tID of the haplotype where the variants come from. \n"
          ";\t3\tvar_id\tID of the variation\n"
          ";\t4\tvar_class\tVariation type, according to SNP Ontology (SO) nomenclature\n"
          ";\t5\tvar_coord\tCoordinates of the variation\n"
          ";\t6\tbest_w\tBest weight for the putative site\n"
          ";\t7\tworst_w\tWorst weight for the putative site\n"
          ";\t8\tw_diff\tDifference between best and worst weight\n"
          ";\t9\tbest_pval\tP_value of the best putative site\n"
          ";\t10\tworst_pval\tP_value of the worst putative site\n"
          ";\t11\tpval_ratio\tRatio between worst and best pval ( pval_ratio = worst_pval/best_pval )\n"
          ";\t12\tbest_variant\tVariant in the best putative site\n"
          ";\t13\tworst_variant\tVariant in the worst putative site\n"
          ";\t14\tbest_offest\tOffset of the best putative site\n"
          ";\t15\tworst_offset\tOffset of the worst putative site\n"
          ";\t16\tmin_offset_diff\tDifference minimal between best and worst putative site\n"
          ";\t17\tbest_strand\tStrand of the best putative site\n"
          ";\t18\tworst_strand\tStrand of the worst putative site\n"
          ";\t19\tstr_change\tIndicate if strand has change between the offset of min_offset_diff\n"
          ";\t20\tbest_seq\tSequence of the worst putative site\n"
          ";\t21\tworst_seq\tSequence of the worst putative site\n"
          ";\t22\tminor_alle_freq\tMinor allele frequency\n"
          "#ac_motif\thap_id\tvar_id\tvar_class\tvar_coord\tbest_w\tworst_w\tw_diff\tbest_pval\tworst_pval\t"
          "pval_ratio\tbest_variant\tworst_variant\tbest_offset\tworst_offset\tmin_offset_diff\tbest_strand\t"
          "worst_strand\tstr_change\tbest_seq\tworst_seq\tminor_allele_freq\n");

  //Close fh
  if(output) {
    fclose(fh_varscan_output);
  }

  return;
}

void ScanHaplosequences(string *mscanquick_file, string *varscanFile, char *matrix_name, int matrix_size, threshold *cutoff) {

  //Declare variables
  FILE *fh_mscanquick_input = NULL;
  FILE *varscan_file        = stdout;
  //FILE *fh_varscan_output   = NULL;

  int  i = 0;
  int  j = 0;

  int nb_vars           = 0;
  int var_offset        = 0;
  int var_length        = 0;
  int real_start_offset = 0;
  int real_end_offset   = 0;
  int eval_start_offset = 0;
  int eval_end_offset   = 0;

  int isCoordDiff       = 0;
  int isAllelDiff       = 0;
  int matrix_length     = 0;

  //char *offset[3]    ={NULL};
  char **token        = NULL;
  char *character     = NULL;
  char **varfield     = NULL;
  char **varcoord     = NULL;
  char **offset_info  = NULL;

  site *curr_site     = NULL;
  scan *curr_scan     = NULL;
  scan *prev_scan     = NULL;
  varscan  *locus        = NULL;
  varscan  *curr_varscan = NULL;
  range    *intersect    = NULL;
  range    *curr_intersect = NULL;

  string *curr_group      = NULL;
  string *line            = NULL;
  string *prev_seq_nb     = NULL;
  string *curr_seq_nb     = NULL;

  string *str_allele1     = NULL;
  string *str_allele2     = NULL;

  char *total_vars        = NULL;
  char *offset_list       = NULL;
  char *offset_element    = NULL;
  char *str_offset        = NULL;
  char *str_length        = NULL;
  char *var_coord         = NULL;
  char *id                = NULL;
  char *so_term           = NULL;
  char *allele1           = NULL;
  char *allele2           = NULL;
  char *tmp_str_allele1   = NULL;
  char *tmp_str_allele2   = NULL;

  char *allele_freq       = NULL;

  int total                 = 0;
  int is_indel              = 0;
  int multiple_variants     = 0;
  int offset                = 0;
  int length                = 0;
  int start                 = 0;
  int end                   = 0;
  int prev_nb_allele_sites  = 0;
  int curr_nb_allele_sites  = 0;

  //Allocate memory for variables
  varfield    = getokens(9);
  varcoord    = getokens(4);
  offset_info = getokens(3);
  token       = getokens(9);

  curr_group   = strnewToList(&RsatMemTracker);
  line         = strnewToList(&RsatMemTracker);
  prev_seq_nb  = strnewToList(&RsatMemTracker);
  curr_seq_nb  = strnewToList(&RsatMemTracker);
  str_allele1  = strnewToList(&RsatMemTracker);
  str_allele2  = strnewToList(&RsatMemTracker);


  //Prepare variables for buffer reading
  line->size = 0;
  token[0] = line->buffer;


  //Initialize strings
  strcopy(curr_group,"");

  //Open input file
  fh_mscanquick_input = OpenInputFile(fh_mscanquick_input, mscanquick_file->buffer);

  //Open output file
  if(strcmp(varscanFile->buffer,"") != 0) {
    varscan_file = OpenAppendFile(varscan_file, varscanFile->buffer);
  }

  while ( fread( (line->buffer + line->size),1,1,fh_mscanquick_input) == 1 ) {
    //If '\t' is found assign next char address as the next token
    if (line->buffer[line->size] == '\t') {
      token[++i] = line->buffer + line->size + 1;
      line->buffer[line->size] = '\0';
    }

    //If end of line is found,proceed to process line
    if (line->buffer[line->size] == '\n') {
      //Reinitialize counters
      line->buffer[line->size] = '\0';
      line->size = 0;
      i = 0;

      //Filters for comments
      if(token[0][0] == '#'){initokadd(line,token,9);continue;}

      /////////////////////////////////////////////////////////////
      // Split the seq_id identifier from matrix scan quick
      // i.e. the first column.
      varfield[0] = token[0];
      for(j = 0; j < 9 - 1; j++) {
        varfield[j+1] = splitstr(varfield[j],';');
      }

      /////////////////////////////////////////////////////////////
      // Split the variant coordinate, separated by '_'
      // i.e. the second varfield.
      varcoord[0] = varfield[1];
      for(j = 0; j < 4 - 1; j++){
        varcoord[j+1] = splitstr( varcoord[j] , '_');
      }
      // NOTE IMPORTANT!! NEED to come back for this chunk!
      //Test if it is a deletion or insertion and move +1 the start
      //because ggctgtgcGCcggctccc the first letter in sequence is the same,
      //        ggctgtgcGcggctccc  and it is not necessary in the comparison
      //if((strcmp(varfield[6],"deletion")  == 0) ||
      //  (strcmp(varfield[6],"insertion") == 0)) eval_start_offset++ ;

      /////////////////////////////////////////////////////////////////////
      //Add and Fetch information to site
      curr_site = sitenew();
      sitefill(curr_site,token[6],token[7],token[8]);

      /////////////////////////////////////////////////////////////////////
      // FIRST group creation of all
      if ( !locus ) {
        ///////////////////////////////////////////////
        // Create varscan object and add it to memlist
        locus = varscanewToList(&RsatMemTracker);
        curr_varscan = locus;
        curr_scan = locus->scan_info;
        prev_scan = curr_scan;

        /////////////////////////////////////////////////////////
        // Split total number of variants from offset and length
        offset_list = SplitOffsetFromTotalVars(varfield[8]);
        total_vars  = varfield[8];
        total       = atoi(total_vars);

        multiple_variants = total > 1 ? 1 : 0;
        ////////////////////////////////////////////////
        // Process only multiple variants for sequences
        if (multiple_variants) {
          //Save list of alleles for hap1 and hap2
          strcopy(str_allele1, varfield[3]);
          strcopy(str_allele2, varfield[2]);

          ////////////////////////////////////////////////////////
          //NOTE WSG. IMPORTANT. I need to come back and fix this
          //and assess effectively how to save the buffer strings
          //Save the real start of strings buffer for
          tmp_str_allele1 = str_allele1->buffer;
          tmp_str_allele2 = str_allele2->buffer;

          /////////////////////////////////////////////
          // Split tokens for variant info and offsets
          //           and create ranges
          //Create ranges
          for (int index = 0; index < total; index++) {
            //Allocate new range
            if(!intersect){
              intersect = rangenewToList(&RsatMemTracker);
              curr_intersect = intersect;
            } else {
              curr_intersect = rangeadd( intersect );
            }


            /////////////////////////////////////////////
            // Split tokens for variant info and offsets

            //Split variant offset from list
            offset_element  = splitstr(offset_list,'|');
            str_offset = offset_list;
            str_length = splitstr(offset_list,'_');

            //Convert char* to int for offset and length
            offset  = atoi(str_offset) + 1;// NOTE: IMPORTANT +1
            length  = atoi(str_length); //due to change in 0-index to 1-index

            //Get real start and end of the sequences
            start = offset - matrix_size + 1;
            end   = offset + length + matrix_size - 1;

            //Split variant coordinates from list
             var_coord = splitstr(varfield[4], ',');
            //Split variant ID from list
             id = splitstr(varfield[5], ',');
            //Split variant SO term from list
             so_term = splitstr(varfield[6], ',');
             //Split variant SO term from list
             allele1 = splitstr(tmp_str_allele1, ',');
            //Split variant Alleles from list
             allele2 = splitstr(tmp_str_allele2, ',');
            //Split variant Allele freq from list
             allele_freq = splitstr(varfield[7], ',');

             //NOTE IMPORTANT WSG. I change this flag state
             //in order to check the presence of an insertion
             //or deletion. so_term = varfield[6]
             if(strcmp(varfield[6],"SNV") != 0 && strcmp(varfield[6],"substitution") != 0){
               is_indel = 1;
             }


             //Fill range
             rangefill(curr_intersect,offset,0, offset + length - 1, 0, length ,0,start,0, end,0);
             varfill(curr_intersect->var_info, "",varfield[4] ,"","+",varfield[5],varfield[6],tmp_str_allele1,tmp_str_allele2,varfield[7]) ;

            //Update tokens first position for further splitting
            offset_list    = offset_element;
            varfield[4]    = var_coord;
            varfield[5]    = id;
            varfield[6]    = so_term;
            tmp_str_allele1    = allele1;
            tmp_str_allele2    = allele2;
            varfield[7]    = allele_freq;
          }
          }
          ///////////////////////////////////////////////////////////
          // Add variation information and the first site
          varfill(curr_varscan->variation,
                  varcoord[0],
                  varcoord[1],
                  varcoord[2],
                  varcoord[3],
                  varfield[5],
                  varfield[6],
                  varfield[4],
                  varfield[2],
                  varfield[7]);
          //curr_scan->offset = -matrix_size + atoi(token[4]) - real_start_offset + 1; //Added a +1 in order to go from [-matrix_length,0]
          curr_scan->offset = atoi(token[4]);
          curr_scan->D = curr_site;

          //Initialize best strand alleles
          curr_varscan->bestD = curr_scan;
          curr_varscan->bestR = curr_scan;
          //Update variant/haplotype index
          strcopy(curr_group,varfield[0]);
          //Update counter of sites
          curr_nb_allele_sites = 1;
          //Continue to next line
          initokadd(line,token,9);
          continue;

      }

      /////////////////////////////////////////////////////////////////////
      //  UPDATE Flags of Diffrent loci and Different Alleles
      isCoordDiff =  (strcmp(varfield[0],curr_group->buffer) == 0) ? 0 : 1;
      isAllelDiff =  (strcmp(varfield[2],curr_varscan->variation->alleles->buffer)  == 0) ? 0 : 1;
      /////////////////////////////////////////////////////////////////////
      // Update previous group of alleles number
      strcopy(curr_group,varfield[0]);

      ///////////////////////////////////////////////////////////////
      ///////////////////////////////////////////////////////////////
      // Add another allele to list
      if(!isCoordDiff && isAllelDiff) {
        //Start counter of new allele
        prev_nb_allele_sites = curr_nb_allele_sites;
        curr_nb_allele_sites = 1;

        //Add new VARSCAN element for a different allele
        curr_varscan = varscanadd(locus);
        curr_scan = curr_varscan->scan_info;
        prev_scan = curr_scan;

        //////////////////////////////////////////////////////////////////
        // Add variation information for the new allele and the first site
        varfill(curr_varscan->variation,
                varcoord[0],
                varcoord[1],
                varcoord[2],
                varcoord[3],
                varfield[5],
                varfield[6],
                varfield[4],
                varfield[2],
                varfield[7]);
        //curr_scan->offset = -matrix_size + atoi(token[4]) - real_start_offset + 1; //Added a +1 in order to go from [-matrix_length,0]
        curr_scan->offset = atoi(token[4]);
        curr_scan->D = curr_site;

        //Initialize best strand alleles
        curr_varscan->bestD = curr_scan;
        curr_varscan->bestR = curr_scan;

        //Continue to next line
        initokadd(line,token,9);
        continue;

      ////////////////////////////////////////////////////////////////
      ///////////////////////////////////////////////////////////////
      // Process the whole locus due to change in genomic coordinates
      } else if(isCoordDiff) {
        ////////////////////////////////////////////////////////////////
        //Check if the locus to process has a different variants or
        //not by testing the existance of the intersect object
        if( !intersect ){
          /////////////////////////////////////////////
          //  Process Haplotype
          processHaplotypes(locus, NULL, is_indel,prev_nb_allele_sites, curr_nb_allele_sites, matrix_size ,matrix_name, cutoff,varscan_file);
          intersect = NULL;
          is_indel = 0;
        ////////////////////////////////////////////////////////////////
        //If there is no intersect object process the locus as single
        //variants
        } else {
          processHaplotypes(locus, intersect, is_indel,prev_nb_allele_sites, curr_nb_allele_sites, matrix_size ,matrix_name, cutoff,varscan_file);
          intersect = NULL;
          is_indel = 0;
          //processLocus(locus, matrix_name, cutoff,varscan_file);
        }
        ////////////////////////////////////////////////////////////////
        // Create new objects for the current line, i.e. new locus
        ///////////////////////////////////////////////
        // Create varscan object and add it to memlist
        locus = varscanewToList(&RsatMemTracker);
        curr_varscan = locus;
        curr_scan = locus->scan_info;
        prev_scan = curr_scan;

        /////////////////////////////////////////////////////////
        // Split total number of variants from offset and length
        //printf("\n This is line %s\n", token[0]);
        //printf("This is varfield[8], previous to SplitOffsetFromTotalVars() %s\n", varfield[8]);

        offset_list = SplitOffsetFromTotalVars(varfield[8]);
        total_vars  = varfield[8];
        total       = atoi(total_vars);
        //printf("\n This is offset_list %s\n", offset_list );

        multiple_variants = total > 1 ? 1 : 0;
        ////////////////////////////////////////////////
        // Process only multiple variants for sequences
        if (multiple_variants) {
          //Save list of alleles for hap1 and hap2
          strcopy(str_allele1, varfield[3]);
          strcopy(str_allele2, varfield[2]);

          ////////////////////////////////////////////////////////
          //NOTE WSG. IMPORTANT. I need to come back and fix this
          //and assess effectively how to save the buffer strings
          //Save the real start of strings buffer for
          tmp_str_allele1 = str_allele1->buffer;
          tmp_str_allele2 = str_allele2->buffer;

          /////////////////////////////////////////////
          // Split tokens for variant info and offsets
          //           and create ranges
          //Create ranges
          for (int index = 0; index < total; index++) {
            //Allocate new range
            if(!intersect){
              intersect = rangenewToList(&RsatMemTracker);
              curr_intersect = intersect;
            } else {
              curr_intersect = rangeadd( intersect );
            }


            /////////////////////////////////////////////
            // Split tokens for variant info and offsets

            //Split variant offset from list
            offset_element  = splitstr(offset_list,'|');
            str_offset = offset_list;
            str_length = splitstr(offset_list,'_');

            //Convert char* to int for offset and length
            offset  = atoi(str_offset) + 1; // NOTE: IMPORTANT +1
            length  = atoi(str_length);

            //Get real start and end of the sequences
            start = offset - matrix_size + 1;
            end   = offset + length + matrix_size - 1;

            //Split variant coordinates from list
             var_coord = splitstr(varfield[4], ',');
            //Split variant ID from list
             id = splitstr(varfield[5], ',');
            //Split variant SO term from list
             so_term = splitstr(varfield[6], ',');
             //Split variant SO term from list
             allele1 = splitstr(tmp_str_allele1, ',');
            //Split variant Alleles from list
             allele2 = splitstr(tmp_str_allele2, ',');
            //Split variant Allele freq from list
             allele_freq = splitstr(varfield[7], ',');

             //NOTE IMPORTANT WSG. I change this flag state
             //in order to check the presence of an insertion
             //or deletion. so_term = varfield[6]
             if(strcmp(varfield[6],"SNV") != 0 && strcmp(varfield[6],"substitution") != 0){
               is_indel = 1;
             }
             //Fill range
             rangefill(curr_intersect,offset,0, offset + length - 1, 0, length ,0,start,0, end,0);
             varfill(curr_intersect->var_info, "",varfield[4] ,"","+",varfield[5],varfield[6],tmp_str_allele1,tmp_str_allele2,varfield[7]) ;

            //Update tokens first position for further splitting
            offset_list      = offset_element;
            varfield[4]      = var_coord;
            varfield[5]      = id;
            varfield[6]      = so_term;
            tmp_str_allele1  = allele1;
            tmp_str_allele2  = allele2;
            varfield[7]      = allele_freq;
          }
          }
          ///////////////////////////////////////////////////////////
          // Add variation information and the first site
          varfill(curr_varscan->variation,
                  varcoord[0],
                  varcoord[1],
                  varcoord[2],
                  varcoord[3],
                  varfield[5],
                  varfield[6],
                  varfield[4],
                  varfield[2],
                  varfield[7]);
          //curr_scan->offset = -matrix_size + atoi(token[4]) - real_start_offset + 1; //Added a +1 in order to go from [-matrix_length,0]
          curr_scan->offset = atoi(token[4]);
          curr_scan->D = curr_site;

          //Initialize best strand alleles
          curr_varscan->bestD = curr_scan;
          curr_varscan->bestR = curr_scan;
          //Update variant/haplotype index
          strcopy(curr_group,varfield[0]);

          //Update number of allele sites
          curr_nb_allele_sites = 1;
          //Continue to next line
          initokadd(line,token,9);
          continue;

      }

      //////////////////////////////////////////////////////////////
      //Add information to current scan offset
      if (strcmp(token[3],"D") == 0) {
        //D scan
        curr_scan = scanadd(curr_scan);
        //curr_scan->offset = -matrix_length + atoi(token[4]) - real_start_offset + 1; //Added a +1 in order to go from [-matrix_length,0]
        curr_scan->offset = atoi(token[4]);
        curr_scan->D = curr_site;
      } else {
        //R scan
        curr_scan->R = curr_site;
        if (curr_varscan->bestD->D->weight < curr_scan->D->weight) curr_varscan->bestD = curr_scan;
        if (curr_varscan->bestR->R->weight < curr_scan->R->weight) curr_varscan->bestR = curr_scan;
        prev_scan = curr_scan;

      }
      //Update number of allele sites
      curr_nb_allele_sites++;

      //If successful continue to the next start of line for reading
      initokadd(line,token,9);
      continue;

    }

    //Resize line string if limit has reached
    limlinetok(line,token,9);

    //Keep track of size count for each char read
    line->size++;

  }
  ///////////////////////////
  //Process remaining varscan
  if ( ! intersect  ) {
    processHaplotypes(locus, NULL , is_indel, prev_nb_allele_sites, curr_nb_allele_sites,matrix_size ,matrix_name, cutoff,varscan_file);
    intersect = NULL;
    //processHaplotypes(locus, intersect ,matrix_name, cutoff,fout);
  } else {
    processHaplotypes(locus, intersect, is_indel, prev_nb_allele_sites, curr_nb_allele_sites,matrix_size ,matrix_name, cutoff,varscan_file);
    intersect = NULL;
    //processLocus(locus, matrix_name, cutoff, varscan_file);
  }


  //Remove tmp variables
  RsatMemTracker = relem((void*)str_allele2, RsatMemTracker);
  RsatMemTracker = relem((void*)str_allele1, RsatMemTracker);
  RsatMemTracker = relem((void*)curr_seq_nb, RsatMemTracker);
  RsatMemTracker = relem((void*)prev_seq_nb, RsatMemTracker);
  RsatMemTracker = relem((void*)line, RsatMemTracker);
  RsatMemTracker = relem((void*)curr_group, RsatMemTracker);
  RsatMemTracker = relem((void*)token, RsatMemTracker);
  RsatMemTracker = relem((void*)offset_info, RsatMemTracker);
  RsatMemTracker = relem((void*)varcoord, RsatMemTracker);
  RsatMemTracker = relem((void*)varfield, RsatMemTracker);



  //Close input matrix-scan-quick fh
  fclose(fh_mscanquick_input);

  //Close output fh if NOT stdout
  if(strcmp(varscanFile->buffer,"") != 0) fclose(varscan_file);

  return;
}

void processHaplotypes(varscan *locus, range *intersect,int is_indel ,int nb_sites1, int nb_sites2, int matrix_size ,char *matrix_name, threshold *cutoff, FILE *fout){
  //Declare variables
  varscan *firstvar  = locus;
  varscan *secndvar  = locus->next;
  scan    *firstscan = firstvar->scan_info;
  scan    *secndscan = secndvar->scan_info;

  // Verbose
  if(verbose >= 13) RsatInfo("Processing locus '",
                             locus->variation->chromosome->buffer,
                             locus->variation->start->buffer,
                             locus->variation->end->buffer,
                             locus->variation->strand,
                             locus->variation->id->buffer,
                             "' at", matrix_name, NULL);

  //If the number of sites added in both alleles are the same, then all are SNV
  //NOTE IMPORTANT Perhaps not? An unlikely combination of indels can give rise to this?
  if ( is_indel ) {
    //Get bestD site for each scan
    firstscan = firstvar->bestD;
    secndscan = secndvar->bestD;
    //Compare strand D
    compareHapAlleles(intersect,is_indel,matrix_size,firstvar,secndvar,firstscan->D,secndscan->D,firstscan->offset,secndscan->offset,"D",cutoff, fout, matrix_name);
    //Compare strand R
    firstscan = firstvar->bestR;
    secndscan = secndvar->bestR;
    compareHapAlleles(intersect,is_indel,matrix_size,firstvar,secndvar,firstscan->R,secndscan->R,firstscan->offset,secndscan->offset,"R",cutoff, fout, matrix_name);

  //Process indels in the haplotypes
  } else {
    //Iterate through every variant in order to make pairwise comparisons
    while (firstscan != NULL) {
      //Compare strand D
      compareHapAlleles(intersect,is_indel,matrix_size,firstvar,secndvar,firstscan->D,secndscan->D,firstscan->offset,secndscan->offset,"D",cutoff, fout, matrix_name);
      //Compare strand R
      compareHapAlleles(intersect,is_indel,matrix_size,firstvar,secndvar,firstscan->R,secndscan->R,firstscan->offset,secndscan->offset,"R",cutoff, fout, matrix_name);

      //Go to next offset
      firstscan = firstscan->next;
      secndscan = secndscan->next;
    }

  }
  //Remove tmp varaibles from locus
  RsatMemTracker = relem((void*)locus,     RsatMemTracker);
  if(intersect) RsatMemTracker = relem((void*)intersect, RsatMemTracker);
  //printf("ALMOST EXIT!\n" );
  return;
}

void compareHapAlleles(range *intersect,int is_indel,int matrix_size,varscan *firstvar,varscan *secndvar,site *firstsite,site *secndsite,int first_offset,int second_offset,char *strand,threshold *cutoff,FILE *fout,char *matrix_name){
  //IMPORTANT
  //NOTE FOR ME: I NEED TO KNOW WHAT TO DO TO THE OFFSETS
  //printf("\n\nENTEREEED compareHapAlleles\n" );
  //printf("1st: %d 2nd: %d\n", first_offset, second_offset);
  //Declare variables
  site *bestsite  = NULL;
  site *worstsite = NULL;

  string *final_varcoord = NULL;
  string *final_id       = NULL;
  string *final_soterm   = NULL;
  string *final_allele1  = NULL;
  string *final_allele2  = NULL;
  string *final_allefreq = NULL;
  string *bestallele     = NULL;
  string *worstallele    = NULL;

  string *final2_varcoord = NULL;
  string *final2_id       = NULL;
  string *final2_soterm   = NULL;
  string *final2_allele1  = NULL;
  string *final2_allele2  = NULL;
  string *final2_allefreq = NULL;

  string *var_coord       = NULL;
  string *id              = NULL;
  string *so_term         = NULL;
  string *allele_freq     = NULL;
  string *final_offset    = NULL;

  variant *bestvariant  = NULL;
  variant *worstvariant = NULL;

  range *intersect_tmp = intersect;

  float bestpval  = 0.0;
  float worstpval = 0.0;

  int bestoffset  = 0;
  int worstoffset = 0;
  int j           = 0;

  //ALlocate memory for variables
  final_varcoord = strnewToList(&RsatMemTracker);
  final_id       = strnewToList(&RsatMemTracker);
  final_soterm   = strnewToList(&RsatMemTracker);
  final_allele1  = strnewToList(&RsatMemTracker);
  final_allele2  = strnewToList(&RsatMemTracker);
  final_allefreq = strnewToList(&RsatMemTracker);
  final_offset   = strnewToList(&RsatMemTracker);

  //Print DBEUG
  //printf("INSIDE compareAlleles This is firstscan BEST %s: %d %s %.2f %.2e  \n",strand,first_offset,firstsite->sequence->buffer,firstsite->weight,firstsite->pval );
  //printf("INSIDE compareAlleles This is secndscan BEST %s: %d %s %.2f %.2e  \n",strand,second_offset,secndsite->sequence->buffer,secndsite->weight,secndsite->pval );

  //Test if this is a single variant
  if(intersect != NULL) {
    //printf("FOUND single variant\n");
    //////////////////////////////////////////////////////////////////////////
    // Get intersection with ranges
    /// Create haplogroup
    while(intersect_tmp != NULL ) {
      //rintf("Starting range comparisons\n");
      //printf("left_flank[0]: %d and right_flank[0]: %d and first_offset: %d \n", intersect_tmp->left_flank[0],intersect_tmp->right_flank[0], first_offset);
      if(first_offset >= intersect_tmp->left_flank[0] && first_offset <= intersect_tmp->end[0] ) {
      if(j == 0){
      strcopy(final_varcoord, intersect_tmp->var_info->start->buffer);
      strcopy(final_id, intersect_tmp->var_info->id->buffer);
      strcopy(final_soterm, intersect_tmp->var_info->SO->buffer);
      strcopy(final_allele2, intersect_tmp->var_info->reference->buffer);
      strcopy(final_allele1, intersect_tmp->var_info->alleles->buffer);
      strcopy(final_allefreq, intersect_tmp->var_info->freq->buffer);
      strfmt(final_offset, "%d",  intersect_tmp->start[0] - (first_offset + matrix_size - 1) ); // Works for SNPS

      } else {
      strccat(final_varcoord,",%s", intersect_tmp->var_info->start->buffer);
      strccat(final_id,",%s", intersect_tmp->var_info->id->buffer);
      strccat(final_soterm,",%s", intersect_tmp->var_info->SO->buffer);
      strccat(final_allele2,",%s", intersect_tmp->var_info->reference->buffer);
      strccat(final_allele1,",%s", intersect_tmp->var_info->alleles->buffer);
      strccat(final_allefreq,",%s", intersect_tmp->var_info->freq->buffer);
      strccat(final_offset, ",%d",  intersect_tmp->start[0] - (first_offset + matrix_size - 1) ); // Works for SNPS

      }
      //Update counter
      j++;
    }
    // Go to next range element
    intersect_tmp = intersect_tmp->next;
  }
    //Update counter
    j = 0;
    //If an insertion/deletion exists in one of the haplotypes, create a second group of var info
    if( is_indel ){
      //printf("Indel found\n");
      intersect_tmp = intersect;
      //ALlocate memory for variables
      final2_varcoord = strnewToList(&RsatMemTracker);
      final2_id       = strnewToList(&RsatMemTracker);
      final2_soterm   = strnewToList(&RsatMemTracker);
      final2_allele1  = strnewToList(&RsatMemTracker);
      final2_allele2  = strnewToList(&RsatMemTracker);
      final2_allefreq = strnewToList(&RsatMemTracker);
      // Create another group for the second offset
      //Second offset
      while(intersect_tmp != NULL ){
        //printf("Starting range comparisons for the second indel\n");

        if(second_offset >= intersect_tmp->left_flank[1] && second_offset <= intersect_tmp->end[1]) {
          if(j == 0){
            //printf("Starting range comparisons for the second indel j == 0\n");

            strcopy(final2_varcoord, intersect_tmp->var_info->start->buffer);
            strcopy(final2_id, intersect_tmp->var_info->id->buffer);
            strcopy(final2_soterm, intersect_tmp->var_info->SO->buffer);
            strcopy(final2_allele2, intersect_tmp->var_info->reference->buffer);
            strcopy(final2_allele1, intersect_tmp->var_info->alleles->buffer);
            strcopy(final2_allefreq, intersect_tmp->var_info->freq->buffer);

          } else {
            //printf("Starting range comparisons for the second indel j > 0\n");

            strccat(final2_varcoord,",%s", intersect_tmp->var_info->start->buffer);
            strccat(final2_id,",%s", intersect_tmp->var_info->id->buffer);
            strccat(final2_soterm,",%s", intersect_tmp->var_info->SO->buffer);
            strccat(final2_allele2,",%s", intersect_tmp->var_info->reference->buffer);
            strccat(final2_allele1,",%s", intersect_tmp->var_info->alleles->buffer);
            strccat(final2_allefreq,",%s", intersect_tmp->var_info->freq->buffer);
          }
          //Update counter
          j++;
        }
        // Go to next range element
        intersect_tmp = intersect_tmp->next;
    }

    }
  } else {
    //printf("\nENTERED default of single variants\n" );
    strcopy(final_allele1,firstvar->variation->alleles->buffer);
    strcopy(final_allele2,secndvar->variation->alleles->buffer);
    strcopy(final_allefreq,firstvar->variation->freq->buffer);
    strcopy(final_id,firstvar->variation->id->buffer);
    strcopy(final_soterm,firstvar->variation->SO->buffer);
    strcopy(final_varcoord,firstvar->variation->reference->buffer); //This is varfield[4], i.e. var coordinate
    strfmt(final_offset, "%d", - (first_offset + matrix_size - 1) + matrix_size  ); // Works for SNPS


    //Adjust names

  }
  ////////////////////////////////////////////
  //Perform comparisons of BOTH sites
  ////////////////////////////////////////////
  //printf("final_allele1 %s\n",final_allele1->buffer );
  //printf("final_allele2 %s\n",final_allele2->buffer );
  //printf("final_allefreq %s\n",final_allefreq->buffer );
  //printf("final_id %s\n",final_id->buffer );
  //printf("final_soterm %s\n",final_soterm->buffer );
  //printf("final_varcoord %s\n",final_varcoord->buffer );
  //printf("Performing sites comparisons\n");
  var_coord   = final_varcoord;
  allele_freq = final_allefreq;
  so_term     = final_soterm;
  id          = final_id;

  if (firstsite->weight > secndsite->weight) {
    bestsite    = firstsite;
    bestvariant = firstvar->variation;
    bestallele  = final_allele1;
    bestoffset  = first_offset;
    //printf("\n");

  } else {
    bestsite    = secndsite;
    bestvariant = secndvar->variation;
    bestallele  = final_allele2;
    bestoffset  = second_offset;

    //For haplotypes with indels
    if(is_indel){
      bestallele  = final2_allele2;
      var_coord   = final2_varcoord;
      allele_freq = final2_allefreq;
      so_term     = final2_soterm;
      id          = final2_id;
    }

  }

  if (firstsite->weight < secndsite->weight) {
    worstsite    = firstsite;
    worstvariant = firstvar->variation;
    worstallele  = final_allele1;
    worstoffset  = first_offset;
    //For haplotypes with indels
    if(is_indel){
       worstallele = final2_allele2;
       strccat(var_coord  ,";%s", final_varcoord->buffer);
       strccat(allele_freq,";%s", final_allefreq->buffer);
       strccat(so_term    ,";%s", final_soterm->buffer);
       strccat(id         ,";%s", final_id->buffer);
     }

  } else {
    worstsite    = secndsite;
    worstvariant = secndvar->variation;
    worstallele  = final_allele2;
    worstoffset  = second_offset;

    //For haplotypes with indels
    if(is_indel){
       worstallele = final2_allele2;
       strccat(var_coord  ,";%s", final2_varcoord->buffer);
       strccat(allele_freq,";%s", final2_allefreq->buffer);
       strccat(so_term    ,";%s", final2_soterm->buffer);
       strccat(id         ,";%s", final2_id->buffer);

    }

  }

  bestpval = (firstsite->pval < secndsite->pval) ? firstsite->pval : secndsite->pval;
  worstpval = (firstsite->pval > secndsite->pval) ? firstsite->pval : secndsite->pval;
  //Filter line output for each threshold
  if( cutoff->upper.pval   && (cutoff->upper.pval   < bestpval) ) return;
  if( (!isnan(cutoff->lower.score)) && (cutoff->lower.score  > bestsite->weight) ) return;
  if( cutoff->lower.wdiff  && (cutoff->lower.wdiff  > bestsite->weight - worstsite->weight) ) return;
  if( cutoff->lower.pratio && (cutoff->lower.pratio > worstpval/bestpval) ) return;




  ///////////////////////////////////////////////////////////////////////////////
  // Prepare to print
  //Declare variables
  int offset_diff = 0;
  int smallest    = 0;
  int biggest     = 0;

  //Get biggest and smallest offset
  if(bestoffset > worstoffset){
    smallest = bestoffset;
    biggest  = worstoffset;
  } else {
    smallest = worstoffset;
    biggest  = bestoffset;
  }
  offset_diff = smallest - biggest;
  //printf("About to print\n");
  //printf("These are the pair of sequences %s %s\n", bestsite->sequence->buffer, worstsite->sequence->buffer );
  fprintf(fout, "%s\t%s:%s-%s_%s\t%s\t%s\t%s\t%.2f\t%.2f\t%.2f\t%.2e\t%.2e\t%.2f\t%s\t%s\t%s\t%s\t%d\t%s\t%s\t0\t%s\t%s\t%s\n",
          matrix_name,
          firstvar->variation->chromosome->buffer, firstvar->variation->start->buffer, firstvar->variation->end->buffer, firstvar->variation->strand,
          id->buffer,
          so_term->buffer,
          var_coord->buffer,
          bestsite->weight,
          worstsite->weight,
          bestsite->weight - worstsite->weight,
          bestsite->pval,
          worstsite->pval,
          worstpval/bestpval,
          bestallele->buffer, //
          worstallele->buffer,//
          final_offset->buffer,//
          final_offset->buffer,//
          offset_diff,//
          strand,
          strand,
          bestsite->sequence->buffer,
          worstsite->sequence->buffer,
          allele_freq->buffer);

          //PRINTED
          //printf("%s\t%s:%s-%s_%s\t%s\t%s\t%s\t%.2f\t%.2f\t%.2f\t%.2e\t%.2e\t%.2f\t%s\t%s\t%d\t%d\t%d\t%s\t%s\t0\t%s\t%s\t%s\n",
          //        matrix_name,
          //        firstvar->variation->chromosome->buffer, firstvar->variation->start->buffer, firstvar->variation->end->buffer, firstvar->variation->strand,
          //        id->buffer,
          //        so_term->buffer,
          //        var_coord->buffer,
          //        bestsite->weight,
          //        worstsite->weight,
          //        bestsite->weight - worstsite->weight,
          //        bestsite->pval,
          //        worstsite->pval,
          //        worstpval/bestpval,
          //        bestallele->buffer, //
          //        worstallele->buffer,//
          //        bestoffset,//
          //        worstoffset,//
          //        offset_diff,//
          //        strand,
          //        strand,
          //        bestsite->sequence->buffer,
          //        worstsite->sequence->buffer,
          //        allele_freq->buffer);
  //Remove tmp variables
  if(final2_varcoord) RsatMemTracker = relem((void*)final2_varcoord, RsatMemTracker);
  if(final2_id)       RsatMemTracker = relem((void*)final2_id, RsatMemTracker);
  if(final2_soterm)   RsatMemTracker = relem((void*)final2_soterm, RsatMemTracker);
  if(final2_allele1)  RsatMemTracker = relem((void*)final2_allele1, RsatMemTracker);
  if(final2_allele2)  RsatMemTracker = relem((void*)final2_allele2, RsatMemTracker);
  if(final2_allefreq) RsatMemTracker = relem((void*)final2_allefreq, RsatMemTracker);

  RsatMemTracker = relem((void*)final_offset, RsatMemTracker);
  RsatMemTracker = relem((void*)final_allefreq, RsatMemTracker);
  RsatMemTracker = relem((void*)final_allele2, RsatMemTracker);
  RsatMemTracker = relem((void*)final_allele1, RsatMemTracker);
  RsatMemTracker = relem((void*)final_soterm, RsatMemTracker);
  RsatMemTracker = relem((void*)final_id, RsatMemTracker);
  RsatMemTracker = relem((void*)final_varcoord, RsatMemTracker);

  return;
}


void ScanSingleVariants(string *mscanquick_file, string *varscanFile, char *matrix_name, int matrix_size, threshold *cutoff){
//void ScanSingleVariants(string *line, char **token, char *matrix_name, string * mscanquick_file, threshold *cutoff, FILE *fout) {
  //Declare variables
  //  printf("\n\nIT ENTERED SCANSINGLEVARIANTS\n" );
  FILE *fh_mscanquick_input = NULL;
  FILE *varscan_file        = stdout;
  int  i = 0;
  int  j = 0;

  int nb_vars           = 0;
  int var_offset        = 0;
  int var_length        = 0;
  int real_start_offset = 0;
  int real_end_offset   = 0;
  int eval_start_offset = 0;
  int eval_end_offset   = 0;

  int isCoordDiff       = 0;
  int isAllelDiff       = 0;
  int matrix_length     = 0;

  //char *offset[3]    ={NULL};
  char **token        = NULL;
  char *character     = NULL;
  char **varfield     = NULL;
  char **varcoord     = NULL;
  char **offset_info  = NULL;

  site *curr_site     = NULL;
  scan *curr_scan     = NULL;
  scan *prev_scan     = NULL;
  varscan  *locus        = NULL;
  varscan  *curr_varscan = NULL;

  string *curr_group      = NULL;
  string *line            = NULL;

  //Allocate memory for variables
  token       = getokens(9);
  varfield    = getokens(8);
  varcoord    = getokens(4);
  offset_info = getokens(3);

  curr_group = strnewToList(&RsatMemTracker);
  line       = strnewToList(&RsatMemTracker);

  //printf("\n\nIT ENTERED SCANSINGLEVARIANTS 2\n" );
  //Prepare variables for buffer reading
  line->size = 0;
  token[0] = line->buffer;



  //printf("\n\nIT ENTERED SCANSINGLEVARIANTS 3\n" );

  //Initialize strings
  strcopy(curr_group,"");
  //printf("\n\nIT ENTERED SCANSINGLEVARIANTS 4\n" );

  ////////////////////////////////////////////////////////////
  //Open input file
  fh_mscanquick_input = OpenInputFile(fh_mscanquick_input, mscanquick_file->buffer);
  //Open output file
  if(strcmp(varscanFile->buffer,"") != 0) {
    varscan_file = OpenAppendFile(varscan_file, varscanFile->buffer);
  }




  while ( fread( (line->buffer + line->size),1,1,fh_mscanquick_input) == 1 ) {
    //If '\t' is found assign next char address as the next token
    if (line->buffer[line->size] == '\t') {
      token[++i] = line->buffer + line->size + 1;
      line->buffer[line->size] = '\0';
    }

    //If end of line is found,proceed to process line
    if (line->buffer[line->size] == '\n') {
      //Reinitialize counters
      line->buffer[line->size] = '\0';
      line->size = 0;
      i = 0;

      //Filters for comments
      if(token[0][0] == '#'){initokadd(line,token,9);continue;}

      //Split 1st field ';' delimited
      character   = token[0];
      varfield[0] = token[0];
      //varcoord[0] = token[0];
      for(int index = 0; character[index] != '\0'; index++) {
        //Split the genomic coordinate in varcoord[0]
        //i.e. chr_start_end_strand
        if(character[index] == '_'){
          character[index] = '\0';
          varcoord[0] = varfield[1];
          varcoord[++j]    = character + index + 1;
        }
        //Split the 8 ';'-separated fields of token[0]
        //i.e. nb;coordinate;id;Allele2;Allele1;SO;Freq;totalvars|offset-length
        if(character[index] == ';') {
          character[index] = '\0';
          varfield[++i] = character + index + 1;
          //if(i == 2) {
            //NOTE WSG. Assess the possibility to include
            //the *Test code here in order to diminish
            //computational cost for splitting the same
            //variant information string i.e. token[0]
          //}
        }
        //Split the last ';'-separated field of token[0]
        //i.e.totalvars|offset-length[|offset-length]{1,n}
        if(character[index] == '|') {
          character[index] = '\0';
          offset_info[0] = varfield[7];
          offset_info[1] = character + index + 1;
          for (index += 1; character[index] != '\0'; index++) {
            if(character[index] == '_') {
              character[index] = '\0';
              offset_info[2] = character + index + 1;
              break;
            }
          }
          break;
        }
      }
      //Reinitialize counters
      i = 0;
      j = 0;
      //Get variant offset information
      //printf("This is offset_info[0]: %s",offset_info[0]);
      //printf("This is offset_info[1]: %s",offset_info[1]);
      //printf("This is offset_info[2]: %s\n",offset_info[2]);

      nb_vars    = atoi(offset_info[0]);
      var_offset = atoi(offset_info[1]) + 1;
      var_length = atoi(offset_info[2]);

      //Get matrix length
      if(!matrix_length) matrix_length = atoi(token[5]) - atoi(token[4]) + 1;
      //Get real offset start and end
      real_start_offset = var_offset - matrix_length + 1;
      real_end_offset   = var_offset + var_length - 1; //NOTE WSG IMPORTANT! This is crucial for deletions and insertions
      //Get eval offset start and end
      eval_start_offset = real_start_offset;
      eval_end_offset   = real_end_offset;
      //Test if it is a deletion or insertion and move +1 the start
      //because ggctgtgcGCcggctccc the first letter in sequence is the same,
      //        ggctgtgcGcggctccc  and it is not necessary in the comparison
      if((strcmp(varfield[5],"deletion")  == 0) ||
        (strcmp(varfield[5],"insertion") == 0)) eval_start_offset++ ;

      //If offset is out of boundaries, skip the line
      if(atoi(token[4]) < eval_start_offset || atoi(token[4]) > eval_end_offset)
        {initokadd(line,token,9);continue;}

      //Add and Fetch information to site
      curr_site = sitenew();
      sitefill(curr_site,token[6],token[7],token[8]);

      if(!locus){
        locus = varscanewToList(&RsatMemTracker);
        curr_varscan = locus;
        curr_scan = locus->scan_info;
        prev_scan = curr_scan;

        /*printf("This is varcoord[0] : %s\t",varcoord[0]);
        printf("This is varcoord[1] : %s\t",varcoord[1]);
        printf("This is varcoord[2] : %s\t",varcoord[2]);
        printf("This is varcoord[3] : %s\t",varcoord[3]);

        printf("This is varfield[0] : %s\t",varfield[0]);
        printf("This is varfield[1] : %s\t",varfield[1]);
        printf("This is varfield[2] : %s\t",varfield[2]);
        printf("This is varfield[3] : %s\t",varfield[3]);
        printf("This is varfield[4] : %s\t",varfield[4]);
        printf("This is varfield[5] : %s\n",varfield[5]);*/



        varfill(curr_varscan->variation,varcoord[0],varcoord[1],varcoord[2],varcoord[3],varfield[4],varfield[5],"",varfield[2],varfield[6]);
        curr_scan->offset = -matrix_length + atoi(token[4]) - real_start_offset + 1; //Added a +1 in order to go from [-matrix_length,0]
        curr_scan->D = curr_site;

        //Initialize best strand alleles
        curr_varscan->bestD = curr_scan;
        curr_varscan->bestR = curr_scan;

        strcopy(curr_group,varfield[0]);
        initokadd(line,token,9);
        continue;
      }
      //*Test for variant information
      /*isCoordDiff = ((strcmp(varcoord[0],curr_varscan->variation->chromosome->buffer) == 0 ) &&
                     (strcmp(varcoord[1],curr_varscan->variation->start->buffer) == 0) &&
                     (strcmp(varcoord[2],curr_varscan->variation->end->buffer)  == 0)) ? 0 : 1;*/
      //isAllelDiff =  (strcmp(varfield[1],curr_varscan->variation->alleles->buffer)  == 0) ? 0 : 1;
      isCoordDiff =  (strcmp(varfield[0],curr_group->buffer) == 0) ? 0 : 1;
      isAllelDiff =  (strcmp(varfield[2],curr_varscan->variation->alleles->buffer)  == 0) ? 0 : 1;
      //Update previous group of alleles number
      strcopy(curr_group,varfield[0]);
      //Add another allele to list
      if(!isCoordDiff && isAllelDiff) {
        curr_varscan = varscanadd(locus);
        curr_scan = curr_varscan->scan_info;
        prev_scan = curr_scan;


        varfill(curr_varscan->variation,varcoord[0],varcoord[1],varcoord[2],varcoord[3],varfield[4],varfield[5],"",varfield[2],varfield[6]);
        curr_scan->offset = -matrix_length + atoi(token[4]) - real_start_offset + 1; //Added a +1 in order to go from [-matrix_length,0]
        curr_scan->D = curr_site;

        //Initialize best strand alleles
        curr_varscan->bestD = curr_scan;
        curr_varscan->bestR = curr_scan;

        initokadd(line,token,9);
        continue;

      //Process the whole locus due to change in
      //genomic coordinates
      } else if(isCoordDiff) {
        processLocus(locus, matrix_name, cutoff,varscan_file);
        locus = varscanewToList(&RsatMemTracker);
        curr_varscan = locus;
        curr_scan = locus->scan_info;
        prev_scan = curr_scan;

        varfill(curr_varscan->variation,varcoord[0],varcoord[1],varcoord[2],varcoord[3],varfield[4],varfield[5],"",varfield[2],varfield[6]);
        curr_scan->offset = -matrix_length + atoi(token[4]) - real_start_offset + 1; //Added a +1 in order to go from [-matrix_length,0]
        curr_scan->D = curr_site;

        //Initialize best strand alleles
        curr_varscan->bestD = curr_scan;
        curr_varscan->bestR = curr_scan;

        initokadd(line,token,9);
        continue;
      }

      //Add information to current scan offset
      if (strcmp(token[3],"D") == 0) {
        /*printf("curr_varscan->scan_info %p\n", (void*)(curr_varscan->scan_info));
        printf("curr_varscan->scan_info->offset %d \n", curr_varscan->scan_info->offset);
        printf("curr_varscan->scan_info->D->sequence %s \n", curr_varscan->scan_info->D->sequence->buffer );
        printf("curr_varscan->scan_info->R->sequence %s \n", curr_varscan->scan_info->R->sequence->buffer );*/
        curr_scan = scanadd(curr_scan);
        curr_scan->offset = -matrix_length + atoi(token[4]) - real_start_offset + 1; //Added a +1 in order to go from [-matrix_length,0]
        curr_scan->D = curr_site;
      } else {
        curr_scan->R = curr_site;
        if (curr_varscan->bestD->D->weight < curr_scan->D->weight) curr_varscan->bestD = curr_scan;
        if (curr_varscan->bestR->R->weight < curr_scan->R->weight) curr_varscan->bestR = curr_scan;
        prev_scan = curr_scan;
        //printf("curr_scan->D. offset :  %d weight: %.2f pval: %.2e sequence: %s \n", curr_scan->offset, curr_scan->D->weight , curr_scan->D->pval,curr_scan->D->sequence->buffer);
        //printf("curr_scan->R. offset :  %d weight: %.2f pval: %.2e sequence: %s \n", curr_scan->offset, curr_scan->R->weight , curr_scan->R->pval,curr_scan->R->sequence->buffer);
        //printf("curr_varscan->bestD->D->weight. offset :  %d weight: %.2f pval: %.2e sequence: %s \n", curr_varscan->bestD->offset, curr_varscan->bestD->D->weight , curr_varscan->bestD->D->pval, curr_varscan->bestD->D->sequence->buffer);
        //printf("curr_varscan->bestR->R->weight. offset :  %d weight: %.2f pval: %.2e sequence: %s \n\n", curr_varscan->bestR->offset, curr_varscan->bestR->R->weight , curr_varscan->bestR->R->pval, curr_varscan->bestD->D->sequence->buffer);

      }

      //If successful continue to the next start of line for reading
      initokadd(line,token,9);
      continue;
    }

    //Resize line string if limit has reached
    limlinetok(line,token,9);

    //Keep track of size count for each char read
    line->size++;

  }
  //Process remaining varscan
  processLocus(locus, matrix_name, cutoff, varscan_file);

  //Remove temporal variables IMPORTANT
  RsatMemTracker = relem((void*)line, RsatMemTracker);

  //Close input matrix-scan-quick fh
  fclose(fh_mscanquick_input);

  //Close output fh if NOT stdout
  if(strcmp(varscanFile->buffer,"") != 0) fclose(varscan_file);
  return;
}

void processLocus(varscan *locus, char *matrix_name, threshold *cutoff, FILE *fout){
  //Tmp strings
  string *SO_tmp1 = NULL;
  string *SO_tmp2 = NULL;
  char *isSO_tmp1 = NULL;
  char *isSO_tmp2 = NULL;

  // Verbose
  if(verbose >= 13) RsatInfo("Processing locus ",
                             locus->variation->chromosome->buffer,
                             locus->variation->start->buffer,
                             locus->variation->end->buffer,
                             locus->variation->strand,
                             locus->variation->id->buffer,
                             "at", matrix_name, NULL);

  //Allocate memory for variables
   SO_tmp1 = strnewToList(&RsatMemTracker);
   SO_tmp2 = strnewToList(&RsatMemTracker);

  //Iterate through every variant in order to make pairwise comparisons
	for(varscan *firstvar =locus; firstvar != NULL; firstvar = firstvar->next){
		for(varscan *secndvar = firstvar->next; secndvar != NULL; secndvar = secndvar->next){

      //NOTE. IMPORTANT WSG. Temporal patch when multiallelic sites are found. The input SO from varSEq
      //The definitive solution is on retrieve-varseq for everythin that is not an haplotype
      //is a list of SNV. However this is failling on strcmp. So I decided temporaly to condense this list and test it each time
      //1       14541   14542   +       rs1045833       SNV,SNV A       C       NA      agaagtccccgccccagctgtgtggcctcCgccagccttccgctccttgaagctggtct
      //1       14541   14542   +       rs1045833       SNV,SNV A       T       NA      agaagtccccgccccagctgtgtggcctcTgccagccttccgctccttgaagctggtct
      //1       14541   14542   +       rs1045833       SNV,SNV A       A       NA      agaagtccccgccccagctgtgtggcctcAgccagccttccgctccttgaagctggtct

      //Copy SO from variants to tmp string
      strcopy(SO_tmp1,firstvar->variation->SO->buffer);
      strcopy(SO_tmp2,secndvar->variation->SO->buffer);

      isSO_tmp1 = isStringListSame(SO_tmp1->buffer, ',');
      isSO_tmp2 = isStringListSame(SO_tmp2->buffer, ',');
      if(isSO_tmp1) strcopy(firstvar->variation->SO,isSO_tmp1);
      if(isSO_tmp2) strcopy(secndvar->variation->SO,isSO_tmp2);


      //printf("This is the first variation %s %s AND this is the second variation %s %s \n",
      //firstvar->variation->id->buffer,firstvar->variation->alleles->buffer,
      //secndvar->variation->id->buffer,secndvar->variation->alleles->buffer);

      //Process comparisons by coordinates if both alleles are SNVs
      if((strcmp(firstvar->variation->SO->buffer,"SNV") == 0) && (strcmp(secndvar->variation->SO->buffer,"SNV") == 0)){
        scan *firstscan = firstvar->scan_info;
        scan *secndscan = secndvar->scan_info;


        //printf("This is the firstscan offset %d AND this is the secondscan offset %d\n",firstscan->offset,secndscan->offset);

        while (firstscan != NULL) {
          //Compare strand D
          compareAlleles(firstvar,secndvar,firstscan->D,secndscan->D,firstscan->offset,secndscan->offset,"D",cutoff, fout, matrix_name);
          //Compare strand R
          compareAlleles(firstvar,secndvar,firstscan->R,secndscan->R,firstscan->offset,secndscan->offset,"R",cutoff, fout, matrix_name);
          //Continue to next offset
          firstscan = firstscan->next;
          secndscan = secndscan->next;
        }
      //Process comaprisons by getting the best site of each allele if an INDEL is found
      } else {
        scan *firstscan = firstvar->bestD;
        scan *secndscan = secndvar->bestD;
        //printf("\nThis is firstscan BEST D: %d %s %.2f %.2e  \n",firstscan->offset,firstscan->D->sequence->buffer,firstscan->D->weight,firstscan->D->pval );
        //printf("This is secndscan BEST D: %d %s %.2f %.2e  \n",secndscan->offset,secndscan->D->sequence->buffer,secndscan->D->weight,secndscan->D->pval );



        //Compare strand D
        compareAlleles(firstvar,secndvar,firstscan->D,secndscan->D,firstscan->offset,secndscan->offset,"D",cutoff, fout, matrix_name);
        //Compare strand R
        firstscan = firstvar->bestR;
        secndscan = secndvar->bestR;
        //printf("This is firstscan BEST R: %d %s %.2f %.2e  \n",firstscan->offset,firstscan->R->sequence->buffer,firstscan->R->weight,firstscan->R->pval );
        //printf("This is secndscan BEST R: %d %s %.2f %.2e  \n",secndscan->offset,secndscan->R->sequence->buffer,secndscan->R->weight,secndscan->R->pval );

        compareAlleles(firstvar,secndvar,firstscan->R,secndscan->R,firstscan->offset,secndscan->offset,"R",cutoff, fout, matrix_name);
      }
		}
	}

  //Remove tmp variables from locus
  RsatMemTracker = relem((void*)SO_tmp1, RsatMemTracker);
  RsatMemTracker = relem((void*)SO_tmp2, RsatMemTracker);
  RsatMemTracker = relem((void*)locus, RsatMemTracker);

  return;
}

void compareAlleles(varscan *firstvar,varscan *secndvar,site *firstsite,site *secndsite,int first_offset,int second_offset,char *strand,threshold *cutoff,FILE *fout,char *matrix_name){
  //Declare variables
  site *bestsite  = NULL;
  site *worstsite = NULL;

  variant *bestvariant  = NULL;
  variant *worstvariant = NULL;

  float bestpval  = 0.0;
  float worstpval = 0.0;

  int bestoffset  = 0;
  int worstoffset = 0;
  //Print DBEUG
  //printf("INSIDE compareAlleles This is firstscan BEST %s: %d %s %.2f %.2e  \n",strand,first_offset,firstsite->sequence->buffer,firstsite->weight,firstsite->pval );
  //printf("INSIDE compareAlleles This is secndscan BEST %s: %d %s %.2f %.2e  \n",strand,second_offset,secndsite->sequence->buffer,secndsite->weight,secndsite->pval );

  //Perform comparisons
  if (firstsite->weight > secndsite->weight) {
    bestsite    = firstsite;
    bestvariant = firstvar->variation;
    bestoffset  = first_offset;
  } else {
    bestsite    = secndsite;
    bestvariant = secndvar->variation;
    bestoffset  = second_offset;

  }
  if (firstsite->weight < secndsite->weight) {
    worstsite    = firstsite;
    worstvariant = firstvar->variation;
    worstoffset  = first_offset;

  } else {
    worstsite    = secndsite;
    worstvariant = secndvar->variation;
    worstoffset  = second_offset;

  }

  bestpval = (firstsite->pval < secndsite->pval) ? firstsite->pval : secndsite->pval;
  worstpval = (firstsite->pval > secndsite->pval) ? firstsite->pval : secndsite->pval;
  //Filter line output for each threshold
  if( cutoff->upper.pval   && (cutoff->upper.pval   < bestpval) ) return;
  if( (!isnan(cutoff->lower.score)) && (cutoff->lower.score  > bestsite->weight) ) return;
  if( cutoff->lower.wdiff  && (cutoff->lower.wdiff  > bestsite->weight - worstsite->weight) ) return;
  if( cutoff->lower.pratio && (cutoff->lower.pratio > worstpval/bestpval) ) return;

  //Print line to output
  printVarScan(fout,matrix_name,firstvar,bestvariant,worstvariant,bestsite,worstsite,bestpval,worstpval,bestoffset,worstoffset,strand);

  return;
}

void printVarScan (FILE *fout,char *matrix_name,varscan *firstvar, variant *bestvariant, variant *worstvariant,site *bestsite, site *worstsite, float bestpval, float worstpval,int bestoffset,int worstoffset, char *strand) {
  //Declare variables
  int offset_diff = 0;
  int smallest    = 0;
  int biggest     = 0;

  //Get biggest and smallest offset
  if(bestoffset > worstoffset){
    smallest = bestoffset;
    biggest  = worstoffset;
  } else {
    smallest = worstoffset;
    biggest  = bestoffset;
  }
  offset_diff = smallest - biggest;
  //printf("These are the pair of sequences %s %s\n", bestsite->sequence->buffer, worstsite->sequence->buffer );
  fprintf(fout, "%s\t%s\t%s\t%s:%s-%s_%s\t%.2f\t%.2f\t%.2f\t%.2e\t%.2e\t%.2f\t%s\t%s\t%d\t%d\t%d\t%s\t%s\t0\t%s\t%s\t%s\n",
          matrix_name,
          firstvar->variation->id->buffer,
          firstvar->variation->SO->buffer,
          firstvar->variation->chromosome->buffer, firstvar->variation->start->buffer, firstvar->variation->end->buffer, firstvar->variation->strand,
          bestsite->weight,
          worstsite->weight,
          bestsite->weight - worstsite->weight,
          bestsite->pval,
          worstsite->pval,
          worstpval/bestpval,
          bestvariant->alleles->buffer,
          worstvariant->alleles->buffer,
          bestoffset,
          worstoffset,
          offset_diff,
          strand,
          strand,
          bestsite->sequence->buffer,
          worstsite->sequence->buffer,
          firstvar->variation->freq->buffer);
  return;
}

//NOTE IMPORTANT WSG PATCH
char *isStringListSame(char *stringtosplit, char separator){
  //Declare variables
  int status = 0;
  char *curr_token = NULL;
  stringlist *list = NULL;
  stringlist *tmp  = NULL;
  stringlist *prev = NULL;
  //Allocate memory for variables
  list = strlistnewToList(&RsatMemTracker);

  //Initialize variables
  curr_token = stringtosplit;
  tmp        = list;
  prev       = tmp;
  strcopy(tmp->element,"");
  //printf("This is stringtosplit %s\n",stringtosplit);
  //Split string
  for (int i = 0; stringtosplit[i] != '\0'; i++) {
    if (stringtosplit[i] == separator) {
      stringtosplit[i] = '\0';
      if (strcmp(tmp->element->buffer,"") == 0) {
        strcopy(tmp->element, curr_token);
        curr_token = stringtosplit + i + 1;
      } else {
        tmp = strlistadd(list);
        strcopy(tmp->element, curr_token);
        curr_token = stringtosplit + i + 1;
        //Compare if both are equals
        if(strcmp(prev->element->buffer, tmp->element->buffer) != 0){
          status = 1;
          break;
        }
        prev = tmp;
      }
    }
  }
  if(status){
    //Remove tmp variables from locus
    RsatMemTracker = relem((void*)list, RsatMemTracker);
    return NULL;
  } else {
    //Last element
    tmp = strlistadd(list);
    strcopy(tmp->element, curr_token);
    //Compare if both are equals
    if(strcmp(prev->element->buffer, tmp->element->buffer) != 0){
      //Remove tmp variables from locus
      RsatMemTracker = relem((void*)list, RsatMemTracker);
      return NULL;
    }
  }
  //Remove tmp variables from locus
  RsatMemTracker = relem((void*)list, RsatMemTracker);
  return stringtosplit;
}

threshold *sethreshold(threshold *cutoff, char *type, char *value) {
  if (strcmp(type,"score") == 0) {
    cutoff->lower.score  = atof(value);
  } else if (strcmp(type,"w_diff") == 0) {
    cutoff->lower.wdiff  = atof(value);
  } else if (strcmp(type,"pval_ratio") == 0) {
    cutoff->lower.pratio = atof(value);
  } else if (strcmp(type,"pval") == 0) {
    cutoff->upper.pval   = atof(value);
  } else {
    RsatFatalError("This is not a valid type -lth/-uth option",type,".Please recheck.",NULL);
  }
  return cutoff;
}
/*Not sure if this is right, varscan is a list with the first node as variant*/
void varscanfree(varscan *delete){
  if(!delete) return;
  scanend(delete->scan_info);
  varfree(delete->variation);
  free(delete);
  return;
}

varscan *varscanew(void){
  //Declare variables
  varscan *new = NULL;
  new = (varscan*)_malloc(sizeof(variant),"varscanew");

  //Initialize attributes
  new->variation = varnew();
  new->scan_info = scanew();
  new->bestD     = new->scan_info;
  new->bestR     = new->scan_info;
  new->next      = NULL;

  return new;
}

varscan *varscanewToList(memstd **List){
  //Declare variables
  varscan *new = NULL;
  //Create element and add it to list
  new = varscanew();
  MemTrackAdd((void*)new, VSCAN, List);

  return new;
}

varscan *varscanadd(varscan *group){
  //Declare variables
  varscan *tmp = NULL;
  varscan *new = NULL;

  new = varscanew();
  //Create tmp
  tmp = group;
  //Attach variant to group
  do {
    if (tmp->next == NULL) {
      tmp->next = new;
      break;
    }
    tmp = tmp->next;
  } while(tmp != NULL );

  if(verbose >= 14) RsatInfo("A variation scan was added to group successfully",NULL);

  return new;
}

/*varscan *varscanfill(){
	//Initialize attributes
  strcopy(element->chromosome,    chrInfo);
  strcopy(element->start     ,  startInfo);
  strcopy(element->end       ,    endInfo);
  if( strandInfo[0] == '-' ) strcpy(element->strand, "-");
  strcopy(element->id       ,      idInfo);
  strcopy(element->SO       ,      SOInfo);
  strcopy(element->reference,     refInfo);
  strcopy(element->alleles  , allelesInfo);
  strcopy(element->freq     ,    freqInfo);

  return element;
}*/

void varscanend(varscan *group){
  varscan *element = NULL;
  if(!group) return;

  do {
    element = group;
    //Pass to next list element
    group = group->next;
    //Release current element
    varscanfree(element);
  } while(group != NULL);

  if(verbose >= 14) RsatInfo("Variation-scan group has been successfully deleted.",NULL);

  return;
}


void sitefree(site *delete){
  strfree(delete->sequence);
  free(delete);
  return;
}

site *sitenew(void){
  //Declare variables
  site *new = NULL;
  new = (site *)_malloc(sizeof(site),"sitenew");

  //Initialize attributes
  new->sequence = strnew();
  new->weight   = 0.0;
  new->pval     = 0.0;

  return new;
}

site *sitenewToList(memstd **List){
  site *new = NULL;
  new = sitenew();
  MemTrackAdd((void*)new,SITE, List);

  return new;
}

site *sitefill(site *element, char *seq, char *weight, char *pval){
  //NOTE WSG IMPORTANT. This is important to assess!
  //Declare variables
  /*int decimals = 0;

  //Prevent decimals form being longer than 2 after the dot
  for (int i = 0; weight[i] != '\0'; i++) {
    if(weight[i] == '.') {
      decimals = 0;
      continue;
    }
    decimals++;
    if(decimals >= 2){
      weight[i + 1] = '\0';
      break;
    }
  }
  //Prevent decimals form being longer than 2 after the dot
  for (int i = 0; pval[i] != '\0'; i++) {
    if(pval[i] == '.') {
      decimals = 0;
      continue;
    }
    decimals++;
    if(decimals >= 10){
      pval[i + 1] = '\0';
      break;
    }
  }*/
  //Initialize attributes
  strcopy(element->sequence,seq);
  element->weight = atof(weight);
  element->pval   = atof(pval);
  //NOTE.WSG. Both weight and pval are truncated for the first 2
  //decimal points
  //element->weight = floorf(atof(weight) * 100)/100;
  //element->pval   = floorf(atof(pval) * 100)/100;

  return element;
}

void scanfree(scan *delete){
  if(!delete) return;
  sitefree(delete->D);
  sitefree(delete->R);
  free(delete);
  return;
}

scan *scanew(void){
  //Declare variables and allocate memory
  scan *new = NULL;
  new = (scan *)_malloc(sizeof(scan),"scanew");

  //Initialize attributes
  new->offset = 0;
  new->D      = NULL;//sitenew();
  new->R      = NULL;//sitenew();
  new->next   = NULL;

  return new;
}

scan *scanewToList(memstd **List){
  scan *new = NULL;
  new = scanew();
  MemTrackAdd((void*)new, SCAN, List);
  return new;
}

scan *scanadd(scan *group){
  //Generate new element
  scan *tmp = NULL;
  scan *new = NULL;
  new = scanew();

  //Create tmp
  tmp = group;
  //Attach variant to group
  do {
    if (tmp->next == NULL) {
      tmp->next = new;
      break;
    }
    tmp = tmp->next;
  } while(tmp != NULL );

  if(verbose >= 14) RsatInfo("Scan was added to group successfully",NULL);

  return new;
}

scan *scanfill(scan *element,char *offset_scan, char *sequence_D, char *weight_D, char *pval_D,char *sequence_R, char *weight_R, char *pval_R ){
  //Initialize attributes
  element->offset = atoi(offset_scan);
  sitefill(element->D,sequence_D,weight_D,pval_D);
  sitefill(element->R,sequence_R,weight_R,pval_R);

  return element;
}

void scanend(scan *group){
  scan *element = NULL;
  if(!group) return;

  do {
    element = group;
    //Pass to next list element
    group = group->next;
    //Release current element
    scanfree(element);
  } while(group != NULL);

  if(verbose >= 14) RsatInfo("Scan group has been successfully deleted.",NULL);

  return;
}

void CreateFastaFromVarseqHaplotypes(string *varsequence, string *fasta_sequence , int matrix_size, unsigned long int *nb_variation, unsigned long int *top_variation, unsigned long int *nb_seq){
  //printf("\n\n AT LEAST ENTERED? with fasta to print: %s\n\n",fasta_sequence->buffer);
  //Declare variables
  string *line                    = NULL;
  string *offset_and_length_tmp   = NULL;
  string *offset_list_hap1        = NULL;
  //string *offset_and_length_final = NULL;
  string *chrom                   = NULL;
  string *chr_start                   = NULL;
  string *chr_end                     = NULL;
  string *sequence1               = NULL;
  string *sequence2               = NULL;
  string *str_allele1             = NULL;
  string *str_allele2             = NULL;
  string *new_offset_hap1         = NULL;
  string *new_offset              = NULL;
  string *final_varcoord          = NULL;
  string *final_id                = NULL;
  string *final_soterm            = NULL;
  string *final_allele1           = NULL;
  string *final_allele2           = NULL;
  string *final_allefreq          = NULL;
  string *haplo_seq1              = NULL;
  string *haplo_seq2              = NULL;


  variant *HaploGroup             = NULL;


  FILE *fh_varsequence      = NULL;
  FILE *fh_fasta_sequence   = NULL;

  char **token          = NULL;
  char *str_offset      = NULL;
  char *str_length      = NULL;
  char *str_length_hap1 = NULL;
  char *sequence        = NULL;
  char *total_vars      = NULL;
  char *offset_list     = NULL;
  char *offset_element_hap1 = NULL;
  char *offset_element  = NULL;
  char *var_coord       = NULL;
  char *id              = NULL;
  char *so_term         = NULL;
  char *tmp_str_allele1 = NULL;
  char *tmp_str_allele2 = NULL;
  char *tmp_str_offset_list_hap1 = NULL;
  char *allele1         = NULL;
  char *allele2         = NULL;
  char *allele_freq     = NULL;

  int i                 = 0;
  int total             = 0;
  int offset_hap1       = 0;
  int offset            = 0;
  int length_hap1       = 0;
  int length            = 0;
  int start_hap1        = 0;
  int start             = 0;
  int end_hap1          = 0;
  int end               = 0;
  int multiple_variants = 0;
  int isAllelDiff       = 0;
  int isLineOdd         = 0;
  int skipHapGroup      = 0;
  int forbidden_char    = 0;


  //Allocate memory for variables
  line                  = strnewToList(&RsatMemTracker);
  offset_and_length_tmp = strnewToList(&RsatMemTracker);
  offset_list_hap1      = strnewToList(&RsatMemTracker);
  chrom                 = strnewToList(&RsatMemTracker);
  chr_start             = strnewToList(&RsatMemTracker);
  chr_end               = strnewToList(&RsatMemTracker);
  sequence1             = strnewToList(&RsatMemTracker);
  sequence2             = strnewToList(&RsatMemTracker);
  str_allele1           = strnewToList(&RsatMemTracker);
  str_allele2           = strnewToList(&RsatMemTracker);
  new_offset_hap1       = strnewToList(&RsatMemTracker);
  new_offset            = strnewToList(&RsatMemTracker);
  final_varcoord        = strnewToList(&RsatMemTracker);
  final_id              = strnewToList(&RsatMemTracker);
  final_soterm          = strnewToList(&RsatMemTracker);
  final_allele1         = strnewToList(&RsatMemTracker);
  final_allele2         = strnewToList(&RsatMemTracker);
  final_allefreq        = strnewToList(&RsatMemTracker);
  haplo_seq1            = strnewToList(&RsatMemTracker);
  haplo_seq2            = strnewToList(&RsatMemTracker);
  token                 = getokens(11);

  //Initialize values
  strfmt(chrom,"");
  strfmt(chr_start,"");
  strfmt(chr_end,"");
  strfmt(offset_and_length_tmp,"");
  strfmt(offset_list_hap1, "");

  ///////////////////////////////////////////////////////
  //Open filehandlers

  //Test if input file IS stdin
  if(strcmp(varsequence->buffer,"") != 0) {
    fh_varsequence    = OpenInputFile(fh_varsequence, varsequence->buffer);
  } else {
    fh_varsequence    = stdin;
  }

  fh_fasta_sequence = OpenOutputFile(fh_fasta_sequence, fasta_sequence->buffer);



  //Prepare variables for buffer reading
  line->size = 0;
  token[0] = line->buffer;

  while ( fread( (line->buffer + line->size),1,1,fh_varsequence) == 1 ) {
    // Test for lines with forbidden characters and save flag
    if(line->buffer[line->size] == ';' ||
       line->buffer[line->size] == '|' ||
       line->buffer[line->size] == '_')
       {
        forbidden_char = 1;

       }

    //If '\t' is found assign next char address as the next token
    if (line->buffer[line->size] == '\t') {
      token[++i] = line->buffer + line->size + 1;
      line->buffer[line->size] = '\0';
      int temporal = i - 1;
      //printf("HELLO %s\n",token[temporal] );
    }

    //If end of line is found,proceed to process line
    if (line->buffer[line->size] == '\n') {
      //Reinitialize counters
      line->buffer[line->size] = '\0';
      line->size = 0;
      //printf("HELLO 10 %s\n",token[i] );
      i = 0;

      //Filters for comments and more
      if(token[0][0] == '#'){forbidden_char = 0;initokadd(line,token,11);continue;}
      if(token[0][0] == ';'){forbidden_char = 0;initokadd(line,token,11);continue;}
      if(token[0][0] == ' '){forbidden_char = 0;initokadd(line,token,11);continue;}
      if(token[0][0] == '\t'){forbidden_char = 0;initokadd(line,token,11);continue;}
      if(token[0][0] == '\0'){forbidden_char = 0;initokadd(line,token,11);continue;}

      //Raise error if forbidden characters where found
      if (forbidden_char) {
        RsatFatalError("A forbidden character ';','|' or '_' was found at variant",
                        token[0],token[1],token[2],token[3],token[4],token[5],token[6],token[7],NULL);
      }

      //Add +1 to sequence counter
      (*nb_seq)++;

      //Test if this is a locus that should be Skipped, i.e. homozygous
      if(skipHapGroup){skipHapGroup = 0;continue;}
      ///////////////////////////////////////////////////
      //Get offset and length of all variants in sequence
      // to assess if more than 1 variant is present

      //printf("This is offset_and_length_tmp before GetVariantIndex:%s\n",offset_and_length_tmp->buffer );
      //printf("THIS IS LINE: %s\n", line->buffer);
      //printf("This is     token[10]         before GetVariantIndex:%s\n",token[10] );

      //Initialize offset_and_length_tmp
      strfmt(offset_and_length_tmp,"");
      GetVariantHapIndex(offset_and_length_tmp,token[8],token[10]);
      //printf("This is offset_and_length_tmp after  GetVariantIndex:%s\n",offset_and_length_tmp->buffer );
      //printf("This is      token[10]        after  GetVariantIndex:%s\n",token[10] );

      //Split total number of variants from offset and length
      offset_list = SplitOffsetFromTotalVars(offset_and_length_tmp->buffer);
      total_vars  = offset_and_length_tmp->buffer;
      total = atoi(total_vars);

      multiple_variants = total > 1 ? 1 : 0;
      //Flags for condition testing
      isAllelDiff = (strcmp(token[7],token[8]) == 0) ? 1 : 0;
      isLineOdd  = (*nb_seq % 2 != 0) ? 1 : 0;

      //If both alleles are the same for the first Haplotype, skip this line
      //and its next haplotype line.
      if(isAllelDiff && isLineOdd){skipHapGroup = 1;continue;}
      //printf("This is the number of total vars: %d\n", total );
      //printf("This is the  offset_length of vars: %s\n\n", offset_list );

      //////////////////////////////////////////////////////////////////////
      // Process only an odd line number when found in multiple variants seq
      if( *nb_seq % 2 != 0  && multiple_variants) {
        //printf("\n\n ENTERED CHUNK 1?\n\n");
        //Save sequence from haplotype 1
        strcopy(sequence1,token[10]);
        //Save list of alleles for hap1 and hap2
        strcopy(str_allele1, token[7]);
        strcopy(str_allele2, token[8]);
        //Save offset list for haplotype 1
        strcopy(offset_list_hap1, offset_list);

        /////////////////////////////////////////////////
        //Only non-overlapping variants will be split
        //if()
        //Allocate memory for variants
        //HaploGroup = varnewToList(&RsatMemTracker);

        //If successful continue to the next start of line for reading
        initokadd(line,token,11);
        continue;
      //////////////////////////////////////////////////////////////////////
      // Process only an even line number when found in multiple variants seq
      } else if(*nb_seq % 2 == 0  && multiple_variants) {
        //printf("\n\n ENTERED CHUNK 2?\n\n");
        ////////////////////////////////////////////////////////
        //NOTE WSG. IMPORTANT. I need to come back and fix this
        //and assess effectively how to save the buffer strings
        //Save the real start of strings buffer for
        tmp_str_allele1          = str_allele1->buffer;
        tmp_str_allele2          = str_allele2->buffer;
        tmp_str_offset_list_hap1 = offset_list_hap1->buffer;

        //Save sequence from haplotype 2
        strcopy(sequence2,token[10]);

        //Start creating ranges
        range *intersect[total];
        int hap_start           = 0;
        int hap_end             = 0;
        int curr_var            = 0;
        //printf("\n\n PASSED THE STRING MANIPULTATION CHUNK 2?\n\n");
        //NOTE IMPORTANT TEMPORAL. Perhaps a linked list to avoid this problem?
        for (int k = 0; k < total; k++) {
          intersect[k] = NULL;
        }
        //printf("\n\n PASSED THE INTERSECT[K] INIT CHUNK 2?\n\n");

        //Create ranges
        for (curr_var = 0; curr_var < total; curr_var++) {
          //Allocate new range
          intersect[curr_var] = rangenew();
          //printf("\n\n PASSED THE RANGENEW CHUNK 2?\n\n");
          /////////////////////////////////////////////
          // Split tokens for variant info and offsets

          //Split variant offset from Haplotype 1 list
          offset_element_hap1  = splitstr(tmp_str_offset_list_hap1,'|');
          //str_offset = offset_list;
          str_length_hap1 = splitstr(tmp_str_offset_list_hap1,'_');
          //Split variant offset from Haplotype 2 list
          offset_element  = splitstr(offset_list,'|');
          str_offset = offset_list;
          str_length = splitstr(offset_list,'_');

          //Convert char* to int for offset and length for Haplotype 1
          offset_hap1  = atoi(tmp_str_offset_list_hap1);
          length_hap1  = atoi(str_length_hap1);
          //Convert char* to int for offset and length for Haplotype 2
          offset  = atoi(str_offset);
          length  = atoi(str_length);

          //Get real start and end of the sequences for Haplotype 2
          start_hap1 = offset_hap1 - matrix_size + 1;
          end_hap1   = offset_hap1 + length_hap1 + matrix_size - 1;
          //Get real start and end of the sequences for Haplotype 2
          start = offset - matrix_size + 1;
          end   = offset + length + matrix_size - 1;

          //printf("token[4] %s token[5] %s token[6] %s token[7] %s token[8] %s token[9] %s\n",token[4], token[5],token[6],token[7],token[8],token[9] );
          //Split variant coordinates from list
           var_coord = splitstr(token[4], ',');
          //Split variant ID from list
           id = splitstr(token[5], ',');
          //Split variant SO term from list
           so_term = splitstr(token[6], ',');
           //Split variant SO term from list
           allele1 = splitstr(tmp_str_allele1, ',');
          //Split variant Alleles from list
           allele2 = splitstr(tmp_str_allele2, ',');
          //Split variant Allele freq from list
           allele_freq = splitstr(token[9], ',');
           //printf("\n\n PASSED THE SPLISTR CHUNK 2?\n\n");

           //Fill range
           rangefill(intersect[curr_var], offset_hap1,offset,
                                          offset_hap1 + length_hap1 - 1, offset + length - 1,
                                          length_hap1, length,
                                          start_hap1, start,
                                          end_hap1, end);
           varfill(intersect[curr_var]->var_info, "",token[4] ,"","+",token[5],token[6],tmp_str_allele1,tmp_str_allele2,token[9]) ;
           //printf("\n\n PASSED THE RANGEFILL AND VARFILL CHUNK 2?\n\n");

          //Update tokens first position for further splitting
          tmp_str_offset_list_hap1 = offset_element_hap1;
          offset_list        = offset_element;
          token[4]           = var_coord;
          token[5]           = id;
          token[6]           = so_term;
          tmp_str_allele1    = allele1;
          tmp_str_allele2    = allele2;
          token[9]           = allele_freq;
          //printf("\n\n PASSED THE REASSIGNMENT CHUNK 2?\n\n");


          /////////////////////////////////////////////////////////////
          //Only test for overlaps when the 2nd variant has been found
          if(curr_var > 0) {
            //printf("\n\n ENTERED CHUNK 4?\n\n");

            //for(hap_end = curr_var - 1; hap_end != hap_start; hap_end--){
            //  if( intersect[curr_var]->start < intersect[hap_end]->end) break;
            //}
            ////-----------------IMPORTANT DELETE-------------
            //printf("start curr_var %d  intersect curr_var-1-end %d matrix size%d\n",
             //intersect[curr_var]->start[0], intersect[curr_var - 1]->end[0], matrix_size);
             //+1 REALLY IMPORTANT
            if(intersect[curr_var]->start[0] - intersect[curr_var - 1]->end[0] +1 > matrix_size ){

              //printf("\n\n ENTERED CHUNK 5?\n\n");

            //If all the previous variants were tested and they do not
            //overlap with the currrent variant split them into different
            //haplotype groups
            //Get new total number of variants for this haplo group
            int total_variants  = curr_var - hap_start;
              ///Create haplogroup
              for (int j = hap_start; j < curr_var; j++) {
                //printf("\n\n ENTERED CHUNK 6?\n\n");

                if(j == hap_start) {
                  strfmt(new_offset_hap1, "%d|%d_%d", total_variants,
                         intersect[j]->start[0] - intersect[hap_start]->left_flank[0],
                         intersect[j]->end[0]   - intersect[j]->start[0] + 1 );
                  strfmt(new_offset, "%d|%d_%d", total_variants,
                         intersect[j]->start[1] - intersect[hap_start]->left_flank[1],
                         intersect[j]->end[1]   - intersect[j]->start[1] + 1 );
                  strcopy(final_varcoord, intersect[j]->var_info->start->buffer);
                  strcopy(final_id, intersect[j]->var_info->id->buffer);
                  strcopy(final_soterm, intersect[j]->var_info->SO->buffer);
                  strcopy(final_allele1, intersect[j]->var_info->reference->buffer);
                  strcopy(final_allele2, intersect[j]->var_info->alleles->buffer);
                  strcopy(final_allefreq, intersect[j]->var_info->freq->buffer);
                } else {
                  strccat(new_offset_hap1, "|%d_%d",
                         intersect[j]->start[0] - intersect[hap_start]->left_flank[0],
                         intersect[j]->end[0]   - intersect[j]->start[0] + 1 );
                  strccat(new_offset, "|%d_%d",
                         intersect[j]->start[1] - intersect[hap_start]->left_flank[1],
                         intersect[j]->end[1]   - intersect[j]->start[1] + 1 );
                  strccat(final_varcoord,",%s", intersect[j]->var_info->start->buffer);
                  strccat(final_id,",%s", intersect[j]->var_info->id->buffer);
                  strccat(final_soterm,",%s", intersect[j]->var_info->SO->buffer);
                  strccat(final_allele1,",%s", intersect[j]->var_info->reference->buffer);
                  strccat(final_allele2,",%s", intersect[j]->var_info->alleles->buffer);
                  strccat(final_allefreq,",%s", intersect[j]->var_info->freq->buffer);
                }
              }
              //Prepare sequences to print
              strcopy(haplo_seq1,sequence1->buffer);
              strcopy(haplo_seq2,sequence2->buffer);
              //-----DELETE-------
              //printf("This is the haplo_seq1 BEFORE %s\n", haplo_seq1->buffer );
              //printf("This is the haplo_seq2 BEFORE %s\n", haplo_seq2->buffer );

              haplo_seq1->buffer[intersect[curr_var - 1]->right_flank[0]] = '\0';
              haplo_seq2->buffer[intersect[curr_var - 1]->right_flank[1]] = '\0';

              //printf("This is the haplo_seq1 AFTER %s\n", haplo_seq1->buffer );
              //printf("This is the haplo_seq2 AFTER %s\n", haplo_seq2->buffer );

              //////////////////////////////////////////////////////
              //Skip the current haplotype if this is an homozygous
              if(strcmp(final_allele1->buffer,final_allele2->buffer) == 0){hap_start = curr_var;continue;}

              //printf("\n\n PASSED THE COPY CHUNK 5?\n\n");
              //Print fasta header and sequence 1
              fprintf(fh_fasta_sequence, ">%lu;%s_%s_%s_%s;%s;%s;%s;%s;%s;%s;%s\n",
              *nb_variation,
              token[0], token[1], token[2], token[3],
              final_allele2->buffer,
              final_allele1->buffer,
              final_varcoord->buffer,
              final_id->buffer,
              final_soterm->buffer,
              final_allefreq->buffer,
              new_offset_hap1->buffer);
              fprintf(fh_fasta_sequence, "%s\n", haplo_seq1->buffer + intersect[hap_start]->left_flank[0]);

              //Print fasta header and sequence 2
              fprintf(fh_fasta_sequence, ">%lu;%s_%s_%s_%s;%s;%s;%s;%s;%s;%s;%s\n",
              *nb_variation,
              token[0], token[1], token[2], token[3],
              final_allele1->buffer,
              final_allele1->buffer,
              final_varcoord->buffer,
              final_id->buffer,
              final_soterm->buffer,
              final_allefreq->buffer,
              new_offset->buffer);
              fprintf(fh_fasta_sequence, "%s\n", haplo_seq2->buffer + intersect[hap_start]->left_flank[1]);

              //printf("\n\n PASSED THE PRINTING CHUNK 5?\n\n");
              //printf( ">%lu;%s_%s_%s_%s;%s;%s;%s;%s;%s;%s;%s\n",
              //*nb_variation,
              //token[0], token[1], token[2], token[3],
              //final_allele2->buffer,
              //final_allele1->buffer,
              //final_varcoord->buffer,
              //final_id->buffer,
              //final_soterm->buffer,
              //final_allefreq->buffer,
              //new_offset_hap1->buffer);
              //printf( "%s\n", haplo_seq1->buffer + intersect[hap_start]->left_flank[0]);

              //Print fasta header and sequence 2
              //printf( ">%lu;%s_%s_%s_%s;%s;%s;%s;%s;%s;%s;%s\n",
              //*nb_variation,
              //token[0], token[1], token[2], token[3],
              //final_allele1->buffer,
              //final_allele1->buffer,
              //final_varcoord->buffer,
              //final_id->buffer,
              //final_soterm->buffer,
              //final_allefreq->buffer,
              //new_offset->buffer);
              //printf( "%s\n", haplo_seq2->buffer + intersect[hap_start]->left_flank[1]);
              //Add +1 to number of variations
              (*nb_variation)++;

              //Update hap_start
              hap_start = curr_var;

            }



          }

        }
        //////////////////////////////////////////////////////////////
        // Print the remaining haplogroup
        int total_variants  = curr_var - hap_start;
          ///Create haplogroup
          //IMPORTANT I changed the 'j <=' for 'j<' because curr_var Updates
          //in the last loop of the previous for WATCH OUT
          for (int j = hap_start; j < curr_var; j++) {
            //printf("\n\n PASSED 1 part THE -REMAINING- PRINTING CHUNK 5?\n\n");

            if(j == hap_start) {
              //printf("HAP_START 1 \n" );
              strfmt(new_offset_hap1, "%d|%d_%d", total_variants,
                     intersect[j]->start[0] - intersect[hap_start]->left_flank[0],
                     intersect[j]->end[0]   - intersect[j]->start[0] + 1 );
              strfmt(new_offset, "%d|%d_%d", total_variants,
                     intersect[j]->start[1] - intersect[hap_start]->left_flank[1],
                     intersect[j]->end[1]   - intersect[j]->start[1] + 1 );
              strcopy(final_varcoord, intersect[j]->var_info->start->buffer);
              strcopy(final_id, intersect[j]->var_info->id->buffer);
              strcopy(final_soterm, intersect[j]->var_info->SO->buffer);
              strcopy(final_allele1, intersect[j]->var_info->reference->buffer);
              strcopy(final_allele2, intersect[j]->var_info->alleles->buffer);
              strcopy(final_allefreq, intersect[j]->var_info->freq->buffer);
            } else {
              //printf("HAP_START 2 \n" );
              strccat(new_offset_hap1, "|%d_%d",
                     intersect[j]->start[0] - intersect[hap_start]->left_flank[0],
                     intersect[j]->end[0]   - intersect[j]->start[0] + 1 );
              strccat(new_offset, "|%d_%d",
                     intersect[j]->start[1] - intersect[hap_start]->left_flank[1],
                     intersect[j]->end[1]   - intersect[j]->start[1] + 1 );
              strccat(final_varcoord,",%s", intersect[j]->var_info->start->buffer);
              strccat(final_id,",%s", intersect[j]->var_info->id->buffer);
              strccat(final_soterm,",%s", intersect[j]->var_info->SO->buffer);
              strccat(final_allele1,",%s", intersect[j]->var_info->reference->buffer);
              strccat(final_allele2,",%s", intersect[j]->var_info->alleles->buffer);
              strccat(final_allefreq,",%s", intersect[j]->var_info->freq->buffer);
            }
          }
          //printf("HAP_START 3 \n" );
          //Skip homozygous haplotypes
          if(strcmp(final_allele1->buffer,final_allele2->buffer) != 0){

          //Prepare sequences to print
          strcopy(haplo_seq1,sequence1->buffer);
          strcopy(haplo_seq2,sequence2->buffer);
          haplo_seq1->buffer[intersect[curr_var - 1]->right_flank[0] ] = '\0';
          haplo_seq2->buffer[intersect[curr_var - 1]->right_flank[1] ] = '\0';
          //printf("\n\n PASSED 2 part THE -REMAINING- PRINTING CHUNK 5?\n\n");
          //Print fasta header and sequence 1
          fprintf(fh_fasta_sequence, ">%lu;%s_%s_%s_%s;%s;%s;%s;%s;%s;%s;%s\n",
          *nb_variation,
          token[0], token[1], token[2], token[3],
          final_allele2->buffer,
          final_allele1->buffer,
          final_varcoord->buffer,
          final_id->buffer,
          final_soterm->buffer,
          final_allefreq->buffer,
          new_offset_hap1->buffer);
          fprintf(fh_fasta_sequence, "%s\n", haplo_seq1->buffer + intersect[hap_start]->left_flank[0]);

          //Print fasta header and sequence 2
          fprintf(fh_fasta_sequence, ">%lu;%s_%s_%s_%s;%s;%s;%s;%s;%s;%s;%s\n",
          *nb_variation,
          token[0], token[1], token[2], token[3],
          final_allele1->buffer,
          final_allele1->buffer,
          final_varcoord->buffer,
          final_id->buffer,
          final_soterm->buffer,
          final_allefreq->buffer,
          new_offset->buffer);
          fprintf(fh_fasta_sequence, "%s\n", haplo_seq2->buffer + intersect[hap_start]->left_flank[1]);
          //printf("\n\n PASSED THE -REMAINING- PRINTING CHUNK 5?\n\n");
          //Print fasta header and sequence 1
          //printf(">%lu;%s_%s_%s_%s;%s;%s;%s;%s;%s;%s;%s\n",
          //*nb_variation,
          //token[0], token[1], token[2], token[3],
          //final_allele2->buffer,
          //final_allele1->buffer,
          //final_varcoord->buffer,
          //final_id->buffer,
          //final_soterm->buffer,
          //final_allefreq->buffer,
          //new_offset_hap1->buffer);
          //printf("%s\n", haplo_seq1->buffer + intersect[hap_start]->left_flank[0]);

          //Print fasta header and sequence 2
          //printf(">%lu;%s_%s_%s_%s;%s;%s;%s;%s;%s;%s;%s\n",
          //*nb_variation,
          //token[0], token[1], token[2], token[3],
          //final_allele1->buffer,
          //final_allele1->buffer,
          //final_varcoord->buffer,
          //final_id->buffer,
          //final_soterm->buffer,
          //final_allefreq->buffer,
          //new_offset->buffer);
          //printf("%s\n", haplo_seq2->buffer + intersect[hap_start]->left_flank[1]);
          //Add +1 to number of variations QUESTION CAREFUL?
          (*nb_variation)++;
        }
        //Remove tmp variables
        RsatMemTracker = relem( (void*)HaploGroup, RsatMemTracker );
        for (curr_var = 0; curr_var < total; curr_var++) {
          //Allocate new range
          rangefree(intersect[curr_var]);
        }
        //If successful continue to the next start of line for reading
        initokadd(line,token,11);
        continue;
      }

      //printf("\n\n IT WENT FOR THE DEFAULT CHUNK ?\n\n");

      //////////////////////////////////////////////////////////////////////
      // Start looking for variant intersections according to their positions
      //in the sequence.
      //////////////////////////////////////////////////////////////////////
      // If only one variant is found process sequence as normal

      //Split offset_and_length
      for (int j = 0; offset_list[j] != '\0'; j++) {
        //Split offset from length
        if( offset_list[j] == '_'){
          offset_list[j] = '\0';
          str_length = offset_list + j + 1;
          break;
        }
      }
    //printf("\n\n PASSED THE DEFAULT SPLIT 1?\n\n");
      //Change variable name for consistency with code below
      str_offset = offset_list;
      //Convert to integers total|offset_length
      //total   = atoi(offset_and_length_tmp->buffer);
      offset  = atoi(str_offset);
      length  = atoi(str_length);

      //Get real start and end of the sequences
      start = offset - matrix_size + 1;
      end   = offset + length + matrix_size - 1;

      offset = offset - start;
      //printf("\n\n PASSED THE DEFAULT CALCULATION OF COORDINATES 1?\n\n");

      //TODO CAREFULLY ASSESS !!!!!!!!!
      token[10][end] = '\0';
      sequence = token[10] + start ;
      //sequence[end+1] = '\0';
      //printf("\n\n PASSED THE DEFAULT SEQUENCE REARREGEMENT?\n\n");
      //printf(">%lu;%s_%s_%s_%s;%s;%s;%s;%s;%s;%s;%d|%d_%d\n",
      //*nb_variation,token[0], token[1], token[2], token[3], token[8], token[7],token[4], token[5], token[6],token[9],
      //total,offset,length);
      //printf("%s\n", sequence);

      //Print fasta header and sequence
      fprintf(fh_fasta_sequence, ">%lu;%s_%s_%s_%s;%s;%s;%s;%s;%s;%s;%d|%d_%d\n",
      *nb_variation,token[0], token[1], token[2], token[3], token[8], token[7],token[4], token[5], token[6],token[9],
      total,offset,length);
      fprintf(fh_fasta_sequence, "%s\n", sequence);
      //TODO  UNTIL HERE !!!!!!!!!!!!!!!!!!

      //If locus is different add +1 to variation counter
      //if( !((strcmp(chrom->buffer, token[0]) == 0) &&
      //      (strcmp(start->buffer, token[1]) == 0) &&
      //      (strcmp(end->buffer, token[2]) == 0))) {
      //NOTE IMPORTANT I need to move this strcmp to the first lines and use
      //a flag for this and also a flag if they are the same in and even line number

      if (strcmp(token[7],token[8]) == 0 ){
            (*nb_variation)++;
            strcopy(chrom,token[0]);
            strcopy(chr_start,token[1]);
            strcopy(chr_end,token[2]);
          }
      //Check if nb of top variations has been reached
      if(*top_variation && !(*nb_variation < *top_variation)) break;


      //If successful continue to the next start of line for reading
      initokadd(line,token,11);
      continue;
    }

    //Resize line string if limit has reached
    limlinetok(line,token,11);

    //Keep track of size count for each char read
    line->size++;

  }

  //Close filehandlers
  //fclose(fh_varsequence);
  fclose(fh_fasta_sequence);
  //Test if input file IS stdin and Close
  if(strcmp(varsequence->buffer,"") != 0) {
    fclose(fh_varsequence);
  }

  //Remove tmp variables
  RsatMemTracker = relem( (void*)token, RsatMemTracker );
  RsatMemTracker = relem( (void*)haplo_seq2     , RsatMemTracker );
  RsatMemTracker = relem( (void*)haplo_seq1     , RsatMemTracker );
  RsatMemTracker = relem( (void*)final_allefreq , RsatMemTracker );
  RsatMemTracker = relem( (void*)final_allele2  , RsatMemTracker );
  RsatMemTracker = relem( (void*)final_allele1  , RsatMemTracker );
  RsatMemTracker = relem( (void*)final_soterm   , RsatMemTracker );
  RsatMemTracker = relem( (void*)final_id       , RsatMemTracker );
  RsatMemTracker = relem( (void*)final_varcoord , RsatMemTracker );
  RsatMemTracker = relem( (void*)new_offset     , RsatMemTracker );
  RsatMemTracker = relem((void*)new_offset_hap1, RsatMemTracker);
  RsatMemTracker = relem((void*)str_allele2, RsatMemTracker);
  RsatMemTracker = relem((void*)str_allele1, RsatMemTracker);
  RsatMemTracker = relem((void*)sequence2, RsatMemTracker);
  RsatMemTracker = relem((void*)sequence1, RsatMemTracker);
  RsatMemTracker = relem((void*)chr_end, RsatMemTracker);
  RsatMemTracker = relem((void*)chr_start, RsatMemTracker);
  RsatMemTracker = relem((void*)chrom, RsatMemTracker);
  RsatMemTracker = relem((void*)offset_list_hap1, RsatMemTracker);
  RsatMemTracker = relem((void*)offset_and_length_tmp, RsatMemTracker);
  RsatMemTracker = relem( (void*)line, RsatMemTracker );

  return;
}

void CreateFastaFromVarseqVariants(string *varsequence, string *fasta_sequence , int matrix_size, unsigned long int *nb_variation, unsigned long int *top_variation, unsigned long int *nb_seq){
  //Declare variables
  string *line                    = NULL;
  string *offset_and_length_tmp   = NULL;
  //string *offset_and_length_final = NULL;
  string *chrom                   = NULL;
  string *chr_start                   = NULL;
  string *chr_end                     = NULL;

  FILE *fh_varsequence      = NULL;
  FILE *fh_fasta_sequence   = NULL;

  char **token     = NULL;
  char *str_offset = NULL;
  char *str_length = NULL;
  char *sequence   = NULL;

  int i = 0;
  int total  = 0;
  int offset = 0;
  int length = 0;
  int start  = 0;
  int end    = 0;
  int forbidden_char = 0;


  //Allocate memory for variables
  line              = strnewToList(&RsatMemTracker);
  offset_and_length_tmp = strnewToList(&RsatMemTracker);
  chrom             = strnewToList(&RsatMemTracker);
  chr_start             = strnewToList(&RsatMemTracker);
  chr_end               = strnewToList(&RsatMemTracker);
  token             = getokens(10);

  //Initialize values
  strfmt(chrom,"");
  strfmt(chr_start,"");
  strfmt(chr_end,"");


  ///////////////////////////////////////////////////////
  //Open filehandlers

  //Test if input file IS stdin
  if(strcmp(varsequence->buffer,"") != 0) {
    fh_varsequence    = OpenInputFile(fh_varsequence, varsequence->buffer);
  } else {
    fh_varsequence    = stdin;
  }

  //Open filehandlers
  fh_fasta_sequence = OpenOutputFile(fh_fasta_sequence, fasta_sequence->buffer);

  //Prepare variables for buffer reading
  line->size = 0;
  token[0] = line->buffer;

  while ( fread( (line->buffer + line->size),1,1,fh_varsequence) == 1 ) {
    // Test for lines with forbidden characters and save flag
    if(line->buffer[line->size] == ';' ||
       line->buffer[line->size] == '|' ||
       line->buffer[line->size] == '_')
       {
        forbidden_char = 1;

       }

    //If '\t' is found assign next char address as the next token
    if (line->buffer[line->size] == '\t') {
      token[++i] = line->buffer + line->size + 1;
      line->buffer[line->size] = '\0';
    }

    //If end of line is found,proceed to process line
    if (line->buffer[line->size] == '\n') {
      //Reinitialize counters
      line->buffer[line->size] = '\0';
      line->size = 0;
      i = 0;


      //Filters for comments and more
      if(token[0][0] == '#'){forbidden_char = 0;initokadd(line,token,10);continue;}
      if(token[0][0] == ';'){forbidden_char = 0;initokadd(line,token,10);continue;}
      if(token[0][0] == ' '){forbidden_char = 0;initokadd(line,token,10);continue;}
      if(token[0][0] == '\t'){forbidden_char = 0;initokadd(line,token,10);continue;}
      if(token[0][0] == '\0'){forbidden_char = 0;initokadd(line,token,10);continue;}

      //Add +1 to sequence counter
      (*nb_seq)++;

      //Raise error if forbidden characters where found
      if (forbidden_char) {
        RsatFatalError("A forbidden character ';','|' or '_' was found at variant",
                        token[0],token[1],token[2],token[3],token[4],token[5],token[6],token[7],NULL);
      }

      //Get offset and length of all variants in sequence
      GetVariantIndex(offset_and_length_tmp,token[9]);
      //strcopy(offset_and_length_final, offset_and_length_tmp->buffer);
      //Split offset_and_length
      for (int j = 0; offset_and_length_tmp->buffer[j] != '\0'; j++) {
        //Split the number of variants from offset and length
        if( offset_and_length_tmp->buffer[j] == '|'){
          offset_and_length_tmp->buffer[j] = '\0';
          str_offset = offset_and_length_tmp->buffer + j + 1;
        }
        //Split offset from length
        if( offset_and_length_tmp->buffer[j] == '_'){
          offset_and_length_tmp->buffer[j] = '\0';
          str_length = offset_and_length_tmp->buffer + j + 1;
          break;
        }
      }

      //Convert to integers total|offset_length
      total   = atoi(offset_and_length_tmp->buffer);
      offset  = atoi(str_offset);
      length  = atoi(str_length);

      //Get real start and end of the sequences
      start = offset - matrix_size + 1;
      end   = offset + length + matrix_size - 1;

      offset = offset - start;

      //TODO CAREFULLY ASSESS !!!!!!!!!
      //int subscript = end +1;
      token[9][end] = '\0';
      sequence = token[9] + start ;
      //printf("\n\nTHIS IS is %s and %c\n",sequence, token[9][end] );
      //sequence[subscript] = '\0';

      //Print fasta header and sequence
      fprintf(fh_fasta_sequence, ">%lu;%s_%s_%s_%s;%s;%s;%s;%s;%s;%d|%d_%d\n",
      *nb_variation,token[0], token[1], token[2], token[3], token[7], token[6],token[4], token[5], token[8],
      total,offset,length);
      fprintf(fh_fasta_sequence, "%s\n", sequence);
      //TODO  UNTIL HERE !!!!!!!!!!!!!!!!!!

      //If locus is different add +1 to variation counter
      /*if( !((strcmp(chrom->buffer, token[0]) == 0) &&
            (strcmp(start->buffer, token[1]) == 0) &&
            (strcmp(end->buffer, token[2]) == 0))) {*/
      if (strcmp(token[6],token[7]) == 0 ){
            (*nb_variation)++;
            strcopy(chrom,token[0]);
            strcopy(chr_start,token[1]);
            strcopy(chr_end,token[2]);
          }
      //Check if nb of top variations has been reached
      //printf("top_variation %lu nb_variation %lu \n", *top_variation, *nb_variation);
      if(*top_variation && !(*nb_variation < *top_variation))break;

      //If successful continue to the next start of line for reading
      initokadd(line,token,10);
      continue;
    }

    //Resize line string if limit has reached
    limlinetok(line,token,10);

    //Keep track of size count for each char read
    line->size++;

  }

  //Close filehandlers
  //fclose(fh_varsequence);
  fclose(fh_fasta_sequence);
  //Test if input file IS stdin and Close
  if(strcmp(varsequence->buffer,"") != 0) {
    fclose(fh_varsequence);
  }

  //Remove tmp variables
  RsatMemTracker = relem( (void*)token, RsatMemTracker );
  RsatMemTracker = relem((void*)chr_end, RsatMemTracker);
  RsatMemTracker = relem((void*)chr_start, RsatMemTracker);
  RsatMemTracker = relem((void*)chrom, RsatMemTracker);
  RsatMemTracker = relem((void*)offset_and_length_tmp, RsatMemTracker);
  RsatMemTracker = relem( (void*)line, RsatMemTracker );

  return;
}

void CreateFastaFromFastaHaplotypes(string *prev_fasta, string *new_fasta , int matrix_size, unsigned long int *nb_variation, unsigned long int *top_variation, unsigned long int *nb_seq){
  /*printf("\n\n AT LEAST ENTERED? with fasta to print: %s\n\n",fasta_sequence->buffer);
  //Declare variables
  string *line                    = NULL;
  string *offset_and_length_tmp   = NULL;
  string *offset_list_hap1        = NULL;
  //string *offset_and_length_final = NULL;
  string *chrom                   = NULL;
  string *chr_start                   = NULL;
  string *chr_end                     = NULL;
  string *sequence1               = NULL;
  string *sequence2               = NULL;
  string *str_allele1             = NULL;
  string *str_allele2             = NULL;
  string *new_offset_hap1         = NULL;
  string *new_offset              = NULL;
  string *final_varcoord          = NULL;
  string *final_id                = NULL;
  string *final_soterm            = NULL;
  string *final_allele1           = NULL;
  string *final_allele2           = NULL;
  string *final_allefreq          = NULL;
  string *haplo_seq1              = NULL;
  string *haplo_seq2              = NULL;


  variant *HaploGroup             = NULL;


  FILE *fh_varsequence      = NULL;
  FILE *fh_fasta_sequence   = NULL;

  char **token          = NULL;
  char *str_offset      = NULL;
  char *str_length      = NULL;
  char *str_length_hap1 = NULL;
  char *sequence        = NULL;
  char *total_vars      = NULL;
  char *offset_list     = NULL;
  char *offset_element_hap1 = NULL;
  char *offset_element  = NULL;
  char *var_coord       = NULL;
  char *id              = NULL;
  char *so_term         = NULL;
  char *tmp_str_allele1 = NULL;
  char *tmp_str_allele2 = NULL;
  char *tmp_str_offset_list_hap1 = NULL;
  char *allele1         = NULL;
  char *allele2         = NULL;
  char *allele_freq     = NULL;

  int i                 = 0;
  int total             = 0;
  int offset_hap1       = 0;
  int offset            = 0;
  int length_hap1       = 0;
  int length            = 0;
  int start_hap1        = 0;
  int start             = 0;
  int end_hap1          = 0;
  int end               = 0;
  int multiple_variants = 0;
  int isAllelDiff       = 0;
  int isLineOdd         = 0;
  int skipHapGroup      = 0;
  unsigned long int nb_line = 0;


  //Allocate memory for variables
  line                  = strnewToList(&RsatMemTracker);
  offset_and_length_tmp = strnewToList(&RsatMemTracker);
  offset_list_hap1      = strnewToList(&RsatMemTracker);
  chrom                 = strnewToList(&RsatMemTracker);
  chr_start             = strnewToList(&RsatMemTracker);
  chr_end               = strnewToList(&RsatMemTracker);
  sequence1             = strnewToList(&RsatMemTracker);
  sequence2             = strnewToList(&RsatMemTracker);
  str_allele1           = strnewToList(&RsatMemTracker);
  str_allele2           = strnewToList(&RsatMemTracker);
  new_offset_hap1       = strnewToList(&RsatMemTracker);
  new_offset            = strnewToList(&RsatMemTracker);
  final_varcoord        = strnewToList(&RsatMemTracker);
  final_id              = strnewToList(&RsatMemTracker);
  final_soterm          = strnewToList(&RsatMemTracker);
  final_allele1         = strnewToList(&RsatMemTracker);
  final_allele2         = strnewToList(&RsatMemTracker);
  final_allefreq        = strnewToList(&RsatMemTracker);
  haplo_seq1            = strnewToList(&RsatMemTracker);
  haplo_seq2            = strnewToList(&RsatMemTracker);
  token                 = getokens(9);

  //Initialize values
  strfmt(chrom,"");
  strfmt(chr_start,"");
  strfmt(chr_end,"");
  strfmt(offset_and_length_tmp,"");
  strfmt(offset_list_hap1, "");

  //Open filehandlers
  fh_varsequence    = OpenInputFile(fh_varsequence, input->buffer);
  fh_fasta_sequence = OpenOutputFile(fh_fasta_sequence, fasta_sequence->buffer);

  //Prepare variables for buffer reading
  line->size = 0;
  token[0] = line->buffer;

  while ( fread( (line->buffer + line->size),1,1,fh_varsequence) == 1 ) {
    //If '\t' is found assign next char address as the next token
    if (line->buffer[line->size] == ';') {
      token[++i] = line->buffer + line->size + 1;
      line->buffer[line->size] = '\0';
      int temporal = i - 1;
      printf("HELLO %s\n",token[temporal] );
    }

    //If end of line is found,proceed to process line
    if (line->buffer[line->size] == '\n') {
      //Reinitialize counters
      line->buffer[line->size] = '\0';
      line->size = 0;
      printf("HELLO 10 %s\n",token[i] );
      i = 0;

      //Filters for comments and more

      //Add +1 to number of line
      nb_line++;
      //Add +1 to sequence counter
      (*nb_seq)++;

      //Test if this is a locus that should be Skipped, i.e. homozygous
      //if(skipHapGroup){skipHapGroup--;continue;}
      ///////////////////////////////////////////////////
      //Get offset and length of all variants in sequence
      // to assess if more than 1 variant is present

      printf("This is offset_and_length_tmp before GetVariantIndex:%s\n",offset_and_length_tmp->buffer );
      printf("THIS IS LINE: %s\n", line->buffer);
      printf("This is     token[10]         before GetVariantIndex:%s\n",token[8] );

      //Initialize offset_and_length_tmp
      strfmt(offset_and_length_tmp,"");
      GetVariantIndex(offset_and_length_tmp,token[8]);
      printf("This is offset_and_length_tmp after  GetVariantIndex:%s\n",offset_and_length_tmp->buffer );
      printf("This is      token[10]        after  GetVariantIndex:%s\n",token[10] );

      //Split total number of variants from offset and length
      offset_list = SplitOffsetFromTotalVars(offset_and_length_tmp->buffer);
      total_vars  = offset_and_length_tmp->buffer;
      total = atoi(total_vars);

      multiple_variants = total > 1 ? 1 : 0;
      //Flags for condition testing
      //isAllelDiff = (strcmp(token[7],token[8]) == 0) ? 1 : 0;
      //isLineOdd  = (*nb_seq % 2 != 0) ? 1 : 0;

      //If both alleles are the same for the first Haplotype, skip this line
      //and its next haplotype line.
      ///if(isAllelDiff && isLineOdd){skipHapGroup = 4;continue;}
      printf("This is the number of total vars: %d\n", total );
      printf("This is the  offset_length of vars: %s\n\n", offset_list );

      //////////////////////////////////////////////////////////////////////
      // Process only an odd line number when found in multiple variants seq
      if( *nb_line % 2 != 0  && multiple_variants) {
        printf("\n\n ENTERED CHUNK 1?\n\n");
        //Save sequence from haplotype 1
        strcopy(sequence1,token[10]);
        //Save list of alleles for hap1 and hap2
        strcopy(str_allele1, token[7]);
        strcopy(str_allele2, token[8]);
        //Save offset list for haplotype 1
        strcopy(offset_list_hap1, offset_list);

        /////////////////////////////////////////////////
        //Only non-overlapping variants will be split
        //if()
        //Allocate memory for variants
        //HaploGroup = varnewToList(&RsatMemTracker);

        //If successful continue to the next start of line for reading
        initokadd(line,token,9);
        continue;
      //////////////////////////////////////////////////////////////////////
      // Process only an even line number when found in multiple variants seq
      } else if(*nb_line % 2 == 0  && multiple_variants) {
        printf("\n\n ENTERED CHUNK 2?\n\n");
        ////////////////////////////////////////////////////////
        //NOTE WSG. IMPORTANT. I need to come back and fix this
        //and assess effectively how to save the buffer strings
        //Save the real start of strings buffer for
        tmp_str_allele1          = str_allele1->buffer;
        tmp_str_allele2          = str_allele2->buffer;
        tmp_str_offset_list_hap1 = offset_list_hap1->buffer;

        //Save sequence from haplotype 2
        strcopy(sequence2,token[10]);

        //Start creating ranges
        range *intersect[total];
        int hap_start           = 0;
        int hap_end             = 0;
        int curr_var            = 0;
        printf("\n\n PASSED THE STRING MANIPULTATION CHUNK 2?\n\n");
        //NOTE IMPORTANT TEMPORAL. Perhaps a linked list to avoid this problem?
        for (int k = 0; k < total; k++) {
          intersect[k] = NULL;
        }
        printf("\n\n PASSED THE INTERSECT[K] INIT CHUNK 2?\n\n");

        //Create ranges
        for (curr_var = 0; curr_var < total; curr_var++) {
          //Allocate new range
          intersect[curr_var] = rangenew();
          printf("\n\n PASSED THE RANGENEW CHUNK 2?\n\n");
          /////////////////////////////////////////////
          // Split tokens for variant info and offsets

          //Split variant offset from Haplotype 1 list
          offset_element_hap1  = splitstr(tmp_str_offset_list_hap1,'|');
          //str_offset = offset_list;
          str_length_hap1 = splitstr(tmp_str_offset_list_hap1,'_');
          //Split variant offset from Haplotype 2 list
          offset_element  = splitstr(offset_list,'|');
          str_offset = offset_list;
          str_length = splitstr(offset_list,'_');

          //Convert char* to int for offset and length for Haplotype 1
          offset_hap1  = atoi(tmp_str_offset_list_hap1);
          length_hap1  = atoi(str_length_hap1);
          //Convert char* to int for offset and length for Haplotype 2
          offset  = atoi(str_offset);
          length  = atoi(str_length);

          //Get real start and end of the sequences for Haplotype 2
          start_hap1 = offset_hap1 - matrix_size + 1;
          end_hap1   = offset_hap1 + length_hap1 + matrix_size - 1;
          //Get real start and end of the sequences for Haplotype 2
          start = offset - matrix_size + 1;
          end   = offset + length + matrix_size - 1;

          printf("token[4] %s token[5] %s token[6] %s token[7] %s token[8] %s token[9] %s\n",token[4], token[5],token[6],token[7],token[8],token[9] );
          //Split variant coordinates from list
           var_coord = splitstr(token[4], ',');
          //Split variant ID from list
           id = splitstr(token[5], ',');
          //Split variant SO term from list
           so_term = splitstr(token[6], ',');
           //Split variant SO term from list
           allele1 = splitstr(tmp_str_allele1, ',');
          //Split variant Alleles from list
           allele2 = splitstr(tmp_str_allele2, ',');
          //Split variant Allele freq from list
           allele_freq = splitstr(token[9], ',');
           printf("\n\n PASSED THE SPLISTR CHUNK 2?\n\n");

           //Fill range
           rangefill(intersect[curr_var], offset_hap1,offset,
                                          offset_hap1 + length_hap1 - 1, offset + length - 1,
                                          length_hap1, length,
                                          start_hap1, start,
                                          end_hap1, end);
           varfill(intersect[curr_var]->var_info, "",token[4] ,"","+",token[5],token[6],tmp_str_allele1,tmp_str_allele2,token[9]) ;
           printf("\n\n PASSED THE RANGEFILL AND VARFILL CHUNK 2?\n\n");

          //Update tokens first position for further splitting
          tmp_str_offset_list_hap1 = offset_element_hap1;
          offset_list        = offset_element;
          token[4]           = var_coord;
          token[5]           = id;
          token[6]           = so_term;
          tmp_str_allele1    = allele1;
          tmp_str_allele2    = allele2;
          token[9]           = allele_freq;
          printf("\n\n PASSED THE REASSIGNMENT CHUNK 2?\n\n");


          /////////////////////////////////////////////////////////////
          //Only test for overlaps when the 2nd variant has been found
          if(curr_var > 0) {
            printf("\n\n ENTERED CHUNK 4?\n\n");

            //for(hap_end = curr_var - 1; hap_end != hap_start; hap_end--){
            //  if( intersect[curr_var]->start < intersect[hap_end]->end) break;
            //}
            ////-----------------IMPORTANT DELETE-------------
            printf("start curr_var %d  intersect curr_var-1-end %d matrix size%d\n",
             intersect[curr_var]->start[0], intersect[curr_var - 1]->end[0], matrix_size);
             //+1 REALLY IMPORTANT
            if(intersect[curr_var]->start[0] - intersect[curr_var - 1]->end[0] +1 > matrix_size ){

              printf("\n\n ENTERED CHUNK 5?\n\n");

            //If all the previous variants were tested and they do not
            //overlap with the currrent variant split them into different
            //haplotype groups
            //Get new total number of variants for this haplo group
            int total_variants  = curr_var - hap_start;
              ///Create haplogroup
              for (int j = hap_start; j < curr_var; j++) {
                printf("\n\n ENTERED CHUNK 6?\n\n");

                if(j == hap_start) {
                  strfmt(new_offset_hap1, "%d|%d_%d", total_variants,
                         intersect[j]->start[0] - intersect[hap_start]->left_flank[0],
                         intersect[j]->end[0]   - intersect[j]->start[0] + 1 );
                  strfmt(new_offset, "%d|%d_%d", total_variants,
                         intersect[j]->start[1] - intersect[hap_start]->left_flank[1],
                         intersect[j]->end[1]   - intersect[j]->start[1] + 1 );
                  strcopy(final_varcoord, intersect[j]->var_info->start->buffer);
                  strcopy(final_id, intersect[j]->var_info->id->buffer);
                  strcopy(final_soterm, intersect[j]->var_info->SO->buffer);
                  strcopy(final_allele1, intersect[j]->var_info->reference->buffer);
                  strcopy(final_allele2, intersect[j]->var_info->alleles->buffer);
                  strcopy(final_allefreq, intersect[j]->var_info->freq->buffer);
                } else {
                  strccat(new_offset_hap1, "|%d_%d",
                         intersect[j]->start[0] - intersect[hap_start]->left_flank[0],
                         intersect[j]->end[0]   - intersect[j]->start[0] + 1 );
                  strccat(new_offset, "|%d_%d",
                         intersect[j]->start[1] - intersect[hap_start]->left_flank[1],
                         intersect[j]->end[1]   - intersect[j]->start[1] + 1 );
                  strccat(final_varcoord,",%s", intersect[j]->var_info->start->buffer);
                  strccat(final_id,",%s", intersect[j]->var_info->id->buffer);
                  strccat(final_soterm,",%s", intersect[j]->var_info->SO->buffer);
                  strccat(final_allele1,",%s", intersect[j]->var_info->reference->buffer);
                  strccat(final_allele2,",%s", intersect[j]->var_info->alleles->buffer);
                  strccat(final_allefreq,",%s", intersect[j]->var_info->freq->buffer);
                }
              }
              //Prepare sequences to print
              strcopy(haplo_seq1,sequence1->buffer);
              strcopy(haplo_seq2,sequence2->buffer);
              //-----DELETE-------
              printf("This is the haplo_seq1 BEFORE %s\n", haplo_seq1->buffer );
              printf("This is the haplo_seq2 BEFORE %s\n", haplo_seq2->buffer );

              haplo_seq1->buffer[intersect[curr_var - 1]->right_flank[0]] = '\0';
              haplo_seq2->buffer[intersect[curr_var - 1]->right_flank[1]] = '\0';

              printf("This is the haplo_seq1 AFTER %s\n", haplo_seq1->buffer );
              printf("This is the haplo_seq2 AFTER %s\n", haplo_seq2->buffer );

              //////////////////////////////////////////////////////
              //Skip the current haplotype if this is an homozygous
              if(strcmp(final_allele1->buffer,final_allele2->buffer) == 0){hap_start = curr_var;continue;}

              printf("\n\n PASSED THE COPY CHUNK 5?\n\n");
              //Print fasta header and sequence 1
              fprintf(fh_fasta_sequence, ">%lu;%s_%s_%s_%s;%s;%s;%s;%s;%s;%s;%s\n",
              *nb_variation,
              token[0], token[1], token[2], token[3],
              final_allele2->buffer,
              final_allele1->buffer,
              final_varcoord->buffer,
              final_id->buffer,
              final_soterm->buffer,
              final_allefreq->buffer,
              new_offset_hap1->buffer);
              fprintf(fh_fasta_sequence, "%s\n", haplo_seq1->buffer + intersect[hap_start]->left_flank[0]);

              //Print fasta header and sequence 2
              fprintf(fh_fasta_sequence, ">%lu;%s_%s_%s_%s;%s;%s;%s;%s;%s;%s;%s\n",
              *nb_variation,
              token[0], token[1], token[2], token[3],
              final_allele1->buffer,
              final_allele1->buffer,
              final_varcoord->buffer,
              final_id->buffer,
              final_soterm->buffer,
              final_allefreq->buffer,
              new_offset->buffer);
              fprintf(fh_fasta_sequence, "%s\n", haplo_seq2->buffer + intersect[hap_start]->left_flank[1]);

              printf("\n\n PASSED THE PRINTING CHUNK 5?\n\n");
              printf( ">%lu;%s_%s_%s_%s;%s;%s;%s;%s;%s;%s;%s\n",
              *nb_variation,
              token[0], token[1], token[2], token[3],
              final_allele2->buffer,
              final_allele1->buffer,
              final_varcoord->buffer,
              final_id->buffer,
              final_soterm->buffer,
              final_allefreq->buffer,
              new_offset_hap1->buffer);
              printf( "%s\n", haplo_seq1->buffer + intersect[hap_start]->left_flank[0]);

              //Print fasta header and sequence 2
              printf( ">%lu;%s_%s_%s_%s;%s;%s;%s;%s;%s;%s;%s\n",
              *nb_variation,
              token[0], token[1], token[2], token[3],
              final_allele1->buffer,
              final_allele1->buffer,
              final_varcoord->buffer,
              final_id->buffer,
              final_soterm->buffer,
              final_allefreq->buffer,
              new_offset->buffer);
              printf( "%s\n", haplo_seq2->buffer + intersect[hap_start]->left_flank[1]);
              //Add +1 to number of variations
              (*nb_variation)++;

              //Update hap_start
              hap_start = curr_var;

            }



          }

        }
        //////////////////////////////////////////////////////////////
        // Print the remaining haplogroup
        int total_variants  = curr_var - hap_start;
          ///Create haplogroup
          //IMPORTANT I changed the 'j <=' for 'j<' because curr_var Updates
          //in the last loop of the previous for WATCH OUT
          for (int j = hap_start; j < curr_var; j++) {
            printf("\n\n PASSED 1 part THE -REMAINING- PRINTING CHUNK 5?\n\n");

            if(j == hap_start) {
              printf("HAP_START 1 \n" );
              strfmt(new_offset_hap1, "%d|%d_%d", total_variants,
                     intersect[j]->start[0] - intersect[hap_start]->left_flank[0],
                     intersect[j]->end[0]   - intersect[j]->start[0] + 1 );
              strfmt(new_offset, "%d|%d_%d", total_variants,
                     intersect[j]->start[1] - intersect[hap_start]->left_flank[1],
                     intersect[j]->end[1]   - intersect[j]->start[1] + 1 );
              strcopy(final_varcoord, intersect[j]->var_info->start->buffer);
              strcopy(final_id, intersect[j]->var_info->id->buffer);
              strcopy(final_soterm, intersect[j]->var_info->SO->buffer);
              strcopy(final_allele1, intersect[j]->var_info->reference->buffer);
              strcopy(final_allele2, intersect[j]->var_info->alleles->buffer);
              strcopy(final_allefreq, intersect[j]->var_info->freq->buffer);
            } else {
              printf("HAP_START 2 \n" );
              strccat(new_offset_hap1, "|%d_%d",
                     intersect[j]->start[0] - intersect[hap_start]->left_flank[0],
                     intersect[j]->end[0]   - intersect[j]->start[0] + 1 );
              strccat(new_offset, "|%d_%d",
                     intersect[j]->start[1] - intersect[hap_start]->left_flank[1],
                     intersect[j]->end[1]   - intersect[j]->start[1] + 1 );
              strccat(final_varcoord,",%s", intersect[j]->var_info->start->buffer);
              strccat(final_id,",%s", intersect[j]->var_info->id->buffer);
              strccat(final_soterm,",%s", intersect[j]->var_info->SO->buffer);
              strccat(final_allele1,",%s", intersect[j]->var_info->reference->buffer);
              strccat(final_allele2,",%s", intersect[j]->var_info->alleles->buffer);
              strccat(final_allefreq,",%s", intersect[j]->var_info->freq->buffer);
            }
          }
          printf("HAP_START 3 \n" );
          //Skip homozygous haplotypes
          if(strcmp(final_allele1->buffer,final_allele2->buffer) != 0){

          //Prepare sequences to print
          strcopy(haplo_seq1,sequence1->buffer);
          strcopy(haplo_seq2,sequence2->buffer);
          haplo_seq1->buffer[intersect[curr_var - 1]->right_flank[0] ] = '\0';
          haplo_seq2->buffer[intersect[curr_var - 1]->right_flank[1] ] = '\0';
          printf("\n\n PASSED 2 part THE -REMAINING- PRINTING CHUNK 5?\n\n");
          //Print fasta header and sequence 1
          fprintf(fh_fasta_sequence, ">%lu;%s_%s_%s_%s;%s;%s;%s;%s;%s;%s;%s\n",
          *nb_variation,
          token[0], token[1], token[2], token[3],
          final_allele2->buffer,
          final_allele1->buffer,
          final_varcoord->buffer,
          final_id->buffer,
          final_soterm->buffer,
          final_allefreq->buffer,
          new_offset_hap1->buffer);
          fprintf(fh_fasta_sequence, "%s\n", haplo_seq1->buffer + intersect[hap_start]->left_flank[0]);

          //Print fasta header and sequence 2
          fprintf(fh_fasta_sequence, ">%lu;%s_%s_%s_%s;%s;%s;%s;%s;%s;%s;%s\n",
          *nb_variation,
          token[0], token[1], token[2], token[3],
          final_allele1->buffer,
          final_allele1->buffer,
          final_varcoord->buffer,
          final_id->buffer,
          final_soterm->buffer,
          final_allefreq->buffer,
          new_offset->buffer);
          fprintf(fh_fasta_sequence, "%s\n", haplo_seq2->buffer + intersect[hap_start]->left_flank[1]);
          printf("\n\n PASSED THE -REMAINING- PRINTING CHUNK 5?\n\n");
          //Print fasta header and sequence 1
          printf(">%lu;%s_%s_%s_%s;%s;%s;%s;%s;%s;%s;%s\n",
          *nb_variation,
          token[0], token[1], token[2], token[3],
          final_allele2->buffer,
          final_allele1->buffer,
          final_varcoord->buffer,
          final_id->buffer,
          final_soterm->buffer,
          final_allefreq->buffer,
          new_offset_hap1->buffer);
          printf("%s\n", haplo_seq1->buffer + intersect[hap_start]->left_flank[0]);

          //Print fasta header and sequence 2
          printf(">%lu;%s_%s_%s_%s;%s;%s;%s;%s;%s;%s;%s\n",
          *nb_variation,
          token[0], token[1], token[2], token[3],
          final_allele1->buffer,
          final_allele1->buffer,
          final_varcoord->buffer,
          final_id->buffer,
          final_soterm->buffer,
          final_allefreq->buffer,
          new_offset->buffer);
          printf("%s\n", haplo_seq2->buffer + intersect[hap_start]->left_flank[1]);
          //Add +1 to number of variations QUESTION CAREFUL?
          (*nb_variation)++;
        }
        //Remove tmp variables
        RsatMemTracker = relem( (void*)HaploGroup, RsatMemTracker );
        for (curr_var = 0; curr_var < total; curr_var++) {
          //Allocate new range
          rangefree(intersect[curr_var]);
        }
        //If successful continue to the next start of line for reading
        initokadd(line,token,9);
        continue;
      }

      printf("\n\n IT WENT FOR THE DEFAULT CHUNK ?\n\n");

      //////////////////////////////////////////////////////////////////////
      // Start looking for variant intersections according to their positions
      //in the sequence.
      //////////////////////////////////////////////////////////////////////
      // If only one variant is found process sequence as normal

      //Split offset_and_length
      for (int j = 0; offset_list[j] != '\0'; j++) {
        //Split offset from length
        if( offset_list[j] == '_'){
          offset_list[j] = '\0';
          str_length = offset_list + j + 1;
          break;
        }
      }
    printf("\n\n PASSED THE DEFAULT SPLIT 1?\n\n");
      //Change variable name for consistency with code below
      str_offset = offset_list;
      //Convert to integers total|offset_length
      //total   = atoi(offset_and_length_tmp->buffer);
      offset  = atoi(str_offset);
      length  = atoi(str_length);

      //Get real start and end of the sequences
      start = offset - matrix_size + 1;
      end   = offset + length + matrix_size - 1;

      offset = offset - start;
      printf("\n\n PASSED THE DEFAULT CALCULATION OF COORDINATES 1?\n\n");

      token[10][end] = '\0';
      sequence = token[10] + start ;
      //sequence[end+1] = '\0';
      printf("\n\n PASSED THE DEFAULT SEQUENCE REARREGEMENT?\n\n");
      printf(">%lu;%s_%s_%s_%s;%s;%s;%s;%s;%s;%s;%d|%d_%d\n",
      *nb_variation,token[0], token[1], token[2], token[3], token[8], token[7],token[4], token[5], token[6],token[9],
      total,offset,length);
      printf("%s\n", sequence);

      //Print fasta header and sequence
      fprintf(fh_fasta_sequence, ">%lu;%s_%s_%s_%s;%s;%s;%s;%s;%s;%s;%d|%d_%d\n",
      *nb_variation,token[0], token[1], token[2], token[3], token[8], token[7],token[4], token[5], token[6],token[9],
      total,offset,length);
      fprintf(fh_fasta_sequence, "%s\n", sequence);

      //If locus is different add +1 to variation counter

      if (strcmp(token[7],token[8]) == 0 ){
            //(*nb_variation)++;
            strcopy(chrom,token[0]);
            strcopy(chr_start,token[1]);
            strcopy(chr_end,token[2]);
          }
      //Check if nb of top variations has been reached
      //if(*top_variation && !(*nb_variation < *top_variation)) break;


      //If successful continue to the next start of line for reading
      initokadd(line,token,9);
      continue;
    }

    //Resize line string if limit has reached
    limlinetok(line,token,9);

    //Keep track of size count for each char read
    line->size++;

  }

  //Close filehandlers
  fclose(fh_varsequence);
  fclose(fh_fasta_sequence);

  //Remove tmp variables
  RsatMemTracker = relem( (void*)token, RsatMemTracker );
  RsatMemTracker = relem( (void*)haplo_seq2     , RsatMemTracker );
  RsatMemTracker = relem( (void*)haplo_seq1     , RsatMemTracker );
  RsatMemTracker = relem( (void*)final_allefreq , RsatMemTracker );
  RsatMemTracker = relem( (void*)final_allele2  , RsatMemTracker );
  RsatMemTracker = relem( (void*)final_allele1  , RsatMemTracker );
  RsatMemTracker = relem( (void*)final_soterm   , RsatMemTracker );
  RsatMemTracker = relem( (void*)final_id       , RsatMemTracker );
  RsatMemTracker = relem( (void*)final_varcoord , RsatMemTracker );
  RsatMemTracker = relem( (void*)new_offset     , RsatMemTracker );
  RsatMemTracker = relem((void*)new_offset_hap1, RsatMemTracker);
  RsatMemTracker = relem((void*)str_allele2, RsatMemTracker);
  RsatMemTracker = relem((void*)str_allele1, RsatMemTracker);
  RsatMemTracker = relem((void*)sequence2, RsatMemTracker);
  RsatMemTracker = relem((void*)sequence1, RsatMemTracker);
  RsatMemTracker = relem((void*)chr_end, RsatMemTracker);
  RsatMemTracker = relem((void*)chr_start, RsatMemTracker);
  RsatMemTracker = relem((void*)chrom, RsatMemTracker);
  RsatMemTracker = relem((void*)offset_list_hap1, RsatMemTracker);
  RsatMemTracker = relem((void*)offset_and_length_tmp, RsatMemTracker);
  RsatMemTracker = relem( (void*)line, RsatMemTracker );

  */
  return;
}

void CreateFastaFromFastaVariants(string *prev_fasta, string *new_fasta , int matrix_size, unsigned long int *nb_variation, unsigned long int *top_variation, unsigned long int *nb_seq){
  //Declare variables
  string *line                    = NULL;
  string *offset_and_length_tmp   = NULL;
  //string *offset_and_length_final = NULL;
  string *chrom                   = NULL;
  string *chr_start                   = NULL;
  string *chr_end                     = NULL;

  FILE *fh_varsequence      = NULL;
  FILE *fh_fasta_sequence   = NULL;

  char **token     = NULL;
  char *str_offset = NULL;
  char *str_length = NULL;
  char *sequence   = NULL;

  int i = 0;
  int total  = 0;
  int offset = 0;
  int length = 0;
  int start  = 0;
  int end    = 0;

  unsigned long int nb_line = 0;
  //Allocate memory for variables
  line              = strnewToList(&RsatMemTracker);
  offset_and_length_tmp = strnewToList(&RsatMemTracker);
  chrom             = strnewToList(&RsatMemTracker);
  chr_start             = strnewToList(&RsatMemTracker);
  chr_end               = strnewToList(&RsatMemTracker);
  token             = getokens(8);

  //Initialize values
  strfmt(chrom,"");
  strfmt(chr_start,"");
  strfmt(chr_end,"");

  //Open filehandlers
  fh_varsequence    = OpenInputFile(fh_varsequence, prev_fasta->buffer);
  fh_fasta_sequence = OpenOutputFile(fh_fasta_sequence, new_fasta->buffer);

  //Prepare variables for buffer reading
  line->size = 0;
  token[0] = line->buffer;

  while ( fread( (line->buffer + line->size),1,1,fh_varsequence) == 1 ) {
    //If '\t' is found assign next char address as the next token
    if (line->buffer[line->size] == ';') {
      token[++i] = line->buffer + line->size + 1;
      line->buffer[line->size] = '\0';
    }

    //If end of line is found,proceed to process line
    if (line->buffer[line->size] == '\n') {
      //Reinitialize counters
      line->buffer[line->size] = '\0';
      line->size = 0;
      i = 0;

      //Filters for comments and more
      //if(token[0][0] == '#'){initokadd(line,token,8);continue;}
      //if(token[0][0] == ';'){initokadd(line,token,8);continue;}
      //if(token[0][0] == ' '){initokadd(line,token,8);continue;}
      //if(token[0][0] == '\t'){initokadd(line,token,8);continue;}
      //if(token[0][0] == '\0'){initokadd(line,token,8);continue;}
      //Add +1 to number of line
      nb_line++;
      //Add +1 to sequence counter
      //(*nb_seq)++;

      //Process header of fasta
      if(nb_line % 2 == 1){
        //Get offset and length of all variants in sequence
        //GetVariantIndex(offset_and_length_tmp,token[9]);
        //strcopy(offset_and_length_final, offset_and_length_tmp->buffer);
        //Split offset_and_length
        for (int j = 0; token[7][j] != '\0'; j++) {
          //Split the number of variants from offset and length
          if( token[7][j] == '|'){
            token[7][j] = '\0';
            str_offset = token[7] + j + 1;
          }
          //Split offset from length
          if( token[7][j] == '_'){
            token[7][j] = '\0';
            str_length = token[7] + j + 1;
            break;
          }
        }

        //Convert to integers total|offset_length
        total   = atoi(token[7]);
        offset  = atoi(str_offset);
        length  = atoi(str_length);

        //Get real start and end of the sequences
        start = offset - matrix_size + 1;
        end   = offset + length + matrix_size - 1;

        offset = offset - start;

        //Print fasta header and sequence
        fprintf(fh_fasta_sequence, "%s;%s;%s;%s;%s;%s;%s;%d|%d_%d\n",
        token[0], token[1], token[2], token[3], token[4], token[5],token[6],total,offset,length);
      } else {
        //int susbscript = end + 1;
        //TODO CAREFULLY ASSESS !!!!!!!!!
        token[0][end] = '\0';
        sequence = token[0] + start ;
        //sequence[susbscript] = '\0';
        fprintf(fh_fasta_sequence, "%s\n", sequence);
      }
      //TODO  UNTIL HERE !!!!!!!!!!!!!!!!!!

      //If locus is different add +1 to variation counter
      /*if( !((strcmp(chrom->buffer, token[0]) == 0) &&
            (strcmp(start->buffer, token[1]) == 0) &&
            (strcmp(end->buffer, token[2]) == 0))) {*/
      /*if (strcmp(token[6],token[7]) == 0 ){
            (*nb_variation)++;
            strcopy(chrom,token[0]);
            strcopy(start,token[1]);
            strcopy(end,token[2]);
          }*/
      //Check if nb of top variations has been reached
      //if(*top_variation && !(*nb_variation < *top_variation)) break;

      //If successful continue to the next start of line for reading
      initokadd(line,token,8);
      continue;
    }

    //Resize line string if limit has reached
    limlinetok(line,token,8);

    //Keep track of size count for each char read
    line->size++;

  }

  //Close filehandlers
  fclose(fh_varsequence);
  fclose(fh_fasta_sequence);

  //Remove tmp variables
  RsatMemTracker = relem( (void*)token, RsatMemTracker );
  RsatMemTracker = relem((void*)chr_end, RsatMemTracker);
  RsatMemTracker = relem((void*)chr_start, RsatMemTracker);
  RsatMemTracker = relem((void*)chrom, RsatMemTracker);
  RsatMemTracker = relem((void*)offset_and_length_tmp, RsatMemTracker);
  RsatMemTracker = relem( (void*)line, RsatMemTracker );

  return;
}

void DEPRECATED_CreateFastaFromHaplosequences(string *line, char **token, unsigned long int *nb_variation, unsigned long int *top_variation, int *nb_seq, FILE *fh_varsequence, FILE *fh_fasta_sequence) {
  //Declare variables
  int i = 0;
  string *offset_and_length = NULL;
  string *chrom  = NULL;
  string *start  = NULL;
  string *end    = NULL;
  //Prepare variables for buffer reading
  line->size = 0;
  token[0] = line->buffer;
  //Allocate memory for variables
  offset_and_length = strnewToList(&RsatMemTracker);
  chrom     = strnewToList(&RsatMemTracker);
  start     = strnewToList(&RsatMemTracker);
  end       = strnewToList(&RsatMemTracker);
  //Initialize values
  strfmt(chrom,"");
  strfmt(start,"");
  strfmt(end,"");

  while ( fread( (line->buffer + line->size),1,1,fh_varsequence) == 1 ) {
    //If '\t' is found assign next char address as the next token
    if (line->buffer[line->size] == '\t') {
      token[++i] = line->buffer + line->size + 1;
      line->buffer[line->size] = '\0';
    }

    //If end of line is found,proceed to process line
    if (line->buffer[line->size] == '\n') {
      //Reinitialize counters
      line->buffer[line->size] = '\0';
      line->size = 0;
      i = 0;
      //If locus is different add +1 to variaiton counter
      if( !((strcmp(chrom->buffer, token[0]) == 0) &&
            (strcmp(start->buffer, token[1]) == 0) &&
            (strcmp(end->buffer, token[2]) == 0)) ) {
              (*nb_variation)++;
              strcopy(chrom,token[0]);
              strcopy(start,token[1]);
              strcopy(end,token[2]);
            }
      //Check if nb of top variations has been reached
      if(*top_variation && !(*nb_variation < *top_variation)) break;
      //Add +1 to sequence counter
      (*nb_seq)++;
      //Get offset and length of all variants in sequence
      GetVariantIndex(offset_and_length,token[10]);

      //Print fasta header and sequence
      fprintf(fh_fasta_sequence, ">%s:%s_%s_%s;%s;%s;%s;%s;%s;%s;%s\n",
      token[0], token[1], token[2], token[3],token[5], token[4],token[7], token[8], token[6], token[9],offset_and_length->buffer);
      fprintf(fh_fasta_sequence, "%s\n", token[9]);
      //If successful continue to the next start of line for reading
      initokadd(line,token,11);
      continue;
    }

    //Resize line string if limit has reached
    limlinetok(line,token,11);

    //Keep track of size count for each char read
    line->size++;

  }
  //Remove tmp variables
  RsatMemTracker = relem((void*)end, RsatMemTracker);
  RsatMemTracker = relem((void*)start, RsatMemTracker);
  RsatMemTracker = relem((void*)chrom, RsatMemTracker);
  RsatMemTracker = relem((void*)offset_and_length, RsatMemTracker);

  return;
}

void DEPRECATED_CreateFastaFromSingleVariants(string *line, char **token, unsigned long int *nb_variation, unsigned long int *top_variation, int *nb_seq, FILE *fh_varsequence, FILE *fh_fasta_sequence) {
  //Declare variables
  int i = 0;
  string *offset_and_length = NULL;
  string *chrom  = NULL;
  string *start  = NULL;
  string *end    = NULL;
  //Prepare variables for buffer reading
  line->size = 0;
  token[0] = line->buffer;
  //Allocate memory for variables
  offset_and_length = strnewToList(&RsatMemTracker);
  chrom     = strnewToList(&RsatMemTracker);
  start     = strnewToList(&RsatMemTracker);
  end       = strnewToList(&RsatMemTracker);
  //Initialize values
  strfmt(chrom,"");
  strfmt(start,"");
  strfmt(end,"");
  //(*nb_variation)++;

  while ( fread( (line->buffer + line->size),1,1,fh_varsequence) == 1 ) {
    //If '\t' is found assign next char address as the next token
    if (line->buffer[line->size] == '\t') {
      token[++i] = line->buffer + line->size + 1;
      line->buffer[line->size] = '\0';
    }

    //If end of line is found,proceed to process line
    if (line->buffer[line->size] == '\n') {
      //Reinitialize counters
      line->buffer[line->size] = '\0';
      line->size = 0;
      i = 0;


      //Add +1 to sequence counter
      (*nb_seq)++;
      //Get offset and length of all variants in sequence
      GetVariantIndex(offset_and_length,token[9]);

      //Print fasta header and sequence
      fprintf(fh_fasta_sequence, ">%lu;%s_%s_%s_%s;%s;%s;%s;%s;%s;%s\n",
      *nb_variation,token[0], token[1], token[2], token[3], token[7], token[6],token[4], token[5], token[8],offset_and_length->buffer);
      fprintf(fh_fasta_sequence, "%s\n", token[9]);

      //If locus is different add +1 to variation counter
      /*if( !((strcmp(chrom->buffer, token[0]) == 0) &&
            (strcmp(start->buffer, token[1]) == 0) &&
            (strcmp(end->buffer, token[2]) == 0))) {*/
      if (strcmp(token[6],token[7]) == 0 ){
            (*nb_variation)++;
            strcopy(chrom,token[0]);
            strcopy(start,token[1]);
            strcopy(end,token[2]);
          }
      //Check if nb of top variations has been reached
      if(*top_variation && !(*nb_variation < *top_variation)) break;

      //If successful continue to the next start of line for reading
      initokadd(line,token,11);
      continue;
    }

    //Resize line string if limit has reached
    limlinetok(line,token,11);

    //Keep track of size count for each char read
    line->size++;

  }
  //Remove tmp variables
  RsatMemTracker = relem((void*)end, RsatMemTracker);
  RsatMemTracker = relem((void*)start, RsatMemTracker);
  RsatMemTracker = relem((void*)chrom, RsatMemTracker);
  RsatMemTracker = relem((void*)offset_and_length, RsatMemTracker);

  return;
}

/* Returns a string in the forma:
total|offset_length[,offset_length]{1,n}
for all the variants in the haplotype.
*/
string *GetVariantHapIndex(string *offset_and_length,char *alleles,char *sequence){
  //Declare variables
  int i               = 0;
  int offset          = 0;
  int var_length      = 0;
  int upcase_stat     = 0;
  int nb_variants     = 0;

  char *allele        = NULL;

  string *tmp         = NULL;
  string *tmp_alleles = NULL;
  variant *variants   = NULL;
  variant *curr_var   = NULL;

  //Allocate memory for variable
  tmp          = strnewToList(&RsatMemTracker);
  tmp_alleles  = strnewToList(&RsatMemTracker);
  variants     = varnewToList(&RsatMemTracker);


  //Format string
  strfmt(tmp,"");
  strfmt(tmp_alleles,"");
  strfmt(offset_and_length, "");

  //Copy alleles to tmp string variables
  strcopy(tmp_alleles, alleles);

  //Initialize first allele token
  allele = tmp_alleles->buffer;
  //printf("This is alleles %s and sequence %s\n", tmp_alleles->buffer, sequence);
  ///////////////////////////////////////////////
  //Split the list of alleles
  for (i = 0; tmp_alleles->buffer[i] != '\0'; i++){

    //Split allele tokens when ',' found
    if(tmp_alleles->buffer[i] == ','){

      //Add +1 to number of variants
      nb_variants++;
      //Insert null
      tmp_alleles->buffer[i] = '\0';
      //Convert size of variant to string
      strfmt(tmp, "%d", strlen(allele));
      //printf("This is %d)size of variant %s with allele %s\n",nb_variants,tmp->buffer, allele );
      //var_length = 0;
      //Add information directly to the first variant
      //NOTE:
      //variant->start will be used to store variant length and
      //variant->alleles will be used to store allele information
      if(nb_variants == 1){
        varfill(variants,"",tmp->buffer,"","","","","",allele,"" );
        curr_var = variants;
      } else {
        curr_var = varadd(variants);
        varfill(curr_var,"",tmp->buffer,"","","","","",allele,"" );
      }

      //Update start of next token
      allele = tmp_alleles->buffer + i + 1;
    }

    //Add +1 to size of the current variant
    var_length++;

  }
  //printf("Number of total vars %d\n",nb_variants );
  ////////////////////////////////////////////
  // Add the remaining variant
  //Add +1 to number of variants
  nb_variants++;
  //printf("Number of total vars %d\n",nb_variants );

  // Convert size of variant to string
  strfmt(tmp, "%d", strlen(allele));
  //Test if no list of variants was found
  //i.e. only a single variation.
  //printf("This is %d)size of variant %s with allele %s\n",nb_variants,tmp->buffer, allele );

  //Single variant processing
  if(curr_var == NULL){
    varfill(variants,"",tmp->buffer,"","","","","",allele,"" );
  //Multiple variants were found
  } else {
    curr_var = varadd(variants);
    varfill(curr_var,"",tmp->buffer,"","","","","",allele,"" );
  }

  //Reinitialize curr_var to the first element in variants
  curr_var = variants;
  //Reinitialize tmp string
  strfmt(tmp,"");

  ////////////////////////////////////////////
  // Start iterating through the sequence
  for (i = 0; sequence[i] != '\0'; i++) {
    upcase_stat = isupper(sequence[i]);
    if( upcase_stat ) {
      //Assign counter value to offset
      offset = i;

      //Update forward the index counter,
      //i.e. add the length of the variant
      // but consider -1 of the i++
      i += atoi(curr_var->start->buffer) - 1;

      //Append offset and length of variant
      strccat(tmp,"|%d_%s", offset, curr_var->start->buffer);

      //Check if this is the last variant of list
      if(curr_var->next == NULL){
        break;
      //Go to next variant
      } else {
        curr_var = curr_var->next;
      }

    }

  }
  //printf("Number of total vars %d\n",nb_variants );
  strfmt(offset_and_length,"%d%s",nb_variants,tmp->buffer);

  //Remove tmp variables
  RsatMemTracker = relem((void*)variants,RsatMemTracker);
  RsatMemTracker = relem((void*)tmp_alleles,RsatMemTracker);
  RsatMemTracker = relem((void*)tmp,RsatMemTracker);

  return offset_and_length;
}

string *GetVariantIndex(string *offset_and_length,char *sequence){
  //Declare variables
  int i           = 0;
  int offset      = 0;
  int var_length  = 0;
  int upcase_stat = 0;
  int nb_variants = 0;
  string *tmp     = NULL;

  //Allocate memory for variable
  tmp = strnewToList(&RsatMemTracker);
  //Format string
  strfmt(tmp,"");
  strfmt(offset_and_length, "");

  //Test if each letter on sequence is uppercase
  for (i = 0; sequence[i] != '\0'; i++) {
    upcase_stat = isupper(sequence[i]);
    if( upcase_stat ) {
      var_length++;
      if(! (var_length > 1)) offset  =  i;
    }
    if (!upcase_stat && var_length != 0) {
      nb_variants++;
      if(nb_variants > 1) strccat(tmp,"|");
      strccat(tmp,"%d_%d", offset, var_length);
      offset     = 0;
      var_length = 0;
    }
  }
  strccat(offset_and_length,"%d|%s",nb_variants,tmp->buffer);

  //Remove tmp variables
  RsatMemTracker = relem((void*)tmp,RsatMemTracker);
  return offset_and_length;
}

/*Splits the output string of GetVariantIndex(),
i.e. total|offset_length[,offset_length]{1,n}.
It modifies the variants_info string by inserting
a '/0' character where the first '|' is found and
returns a pointer to the total of variatns and to
the start of the offset list.
 */
char *SplitOffsetFromTotalVars(char *variants_info){
  //Initialize variables
  char *str_offset = NULL;

  //Split offset_and_length
  for (int j = 0; variants_info[j] != '\0'; j++) {
    //Split the number of variants from offset and length
    if( variants_info[j] == '|'){
      variants_info[j] = '\0';
      str_offset = variants_info + j + 1;
      break;
    }

  }

  return str_offset;

}
/**/

/*Split a char array into tokens separated by a single
  character. Returns a pointer to the start of the next
  token*/
char *splitstr(char *str_to_split, char sep){

  for (int i = 0; str_to_split[i] != '\0' ; i++) {
    if(str_to_split[i] == sep ){
      str_to_split[i] = '\0';
      return str_to_split + i + 1;

    }
  }

  return NULL;
}
///////////////////////////////////////
///////////////////////////////////////
///////////////////////////////////////
///////////////////////////////////////
/*List available input/output matrix formats*/
TRIE *ListInputMatrixFormats(void){
  //Declare variables
  TRIE *trieMtxfmts = NULL;

  //Start trie
  trieMtxfmts = TrieStart();

  //Insert key-value to trie
  TrieInsert(trieMtxfmts,"alignace","1");
  TrieInsert(trieMtxfmts,"assembly","1");
  TrieInsert(trieMtxfmts,"cb","1");
  TrieInsert(trieMtxfmts,"clustal","1");
  TrieInsert(trieMtxfmts,"cluster-buster","1");
  TrieInsert(trieMtxfmts,"consensus","1");
  TrieInsert(trieMtxfmts,"sequences","1");
  TrieInsert(trieMtxfmts,"feature","1");
  TrieInsert(trieMtxfmts,"footprintdb","1");
  TrieInsert(trieMtxfmts,"gibbs","1");
  TrieInsert(trieMtxfmts,"infogibbs","1");
  TrieInsert(trieMtxfmts,"info-gibbs","1");
  TrieInsert(trieMtxfmts,"jaspar","1");
  TrieInsert(trieMtxfmts,"homer","1");
  TrieInsert(trieMtxfmts,"mscan","1");
  TrieInsert(trieMtxfmts,"meme","1");
  TrieInsert(trieMtxfmts,"meme_block","1");
  TrieInsert(trieMtxfmts,"motifsampler","1");
  TrieInsert(trieMtxfmts,"stamp","1");
  TrieInsert(trieMtxfmts,"stamp-transfac","1");
  TrieInsert(trieMtxfmts,"tab","1");
  TrieInsert(trieMtxfmts,"tf","1");
  TrieInsert(trieMtxfmts,"transfac","1");
  TrieInsert(trieMtxfmts,"cis-bp","1");
  TrieInsert(trieMtxfmts,"uniprobe","1");
  TrieInsert(trieMtxfmts,"yeastract","1");
  TrieInsert(trieMtxfmts,"encode","1");

  return trieMtxfmts;
}

void argToList(stringlist *list, char *value){
  //Declare variables
  stringlist *tmp = NULL;
  //Add information to string list
  if( strcmp(list->element->buffer,"") == 0 ) {
     strcopy(list->element, value);
  } else {
     tmp = strlistadd(list);
     strcopy(tmp->element, value);
  }

  return;
}

string *GetProgramPath(string *program_path, char *program_name, int die_on_error, stringlist *preferred_path){
  //Declare variables
  //int path_found = 0;
  string *possible_path = NULL;
  string *cmd_which     = NULL;
  FILE *fh_popen        = NULL;

  //Check for Global variants in order to search their directories
  if(strcmp(BIN->buffer    , "") == 0) RsatFatalError("Global variable BIN is not defined"     , NULL);
  if(strcmp(PYTHON->buffer , "") == 0) RsatFatalError("Global variable PYTHON is not defined"  , NULL);
  if(strcmp(SCRIPTS->buffer, "") == 0) RsatFatalError("Global variable SCRIPTS is not defined" , NULL);

  //Allcoate memory for variables
  possible_path = strnewToList(&RsatMemTracker);
  cmd_which     = strnewToList(&RsatMemTracker);

  //Initialize program_path as empty string
  strfmt(program_path,"");

  //If a list of directories has been supplied test if program exists there
  if (preferred_path != NULL) {
    while ( preferred_path != NULL ) {
      strfmt(possible_path, "%s/%s", preferred_path->element->buffer, program_name);
      if ( access(possible_path->buffer,F_OK) != -1 ) { //NOTE DO NOT FORGET
        strcopy(program_path,possible_path->buffer);
        break;
      }
      preferred_path = preferred_path->next;
    }
  }

  //Check if program exists on common RSAT directories: BIN, PYTHON and SCRIPTS
  if(strcmp(program_path->buffer, "") == 0) {
    strfmt(possible_path,"%s/%s", BIN->buffer,     program_name);
    if ( (access(possible_path->buffer,F_OK) != -1) && (access(possible_path->buffer,X_OK) != -1) ) {
      strcopy(program_path,possible_path->buffer);
      //Remove tmp allocated variables
      RsatMemTracker = relem((void*)possible_path, RsatMemTracker);
      RsatMemTracker = relem((void*)cmd_which    , RsatMemTracker);
      if (verbose >= 4) RsatInfo("GetProgramPath() path found", program_path->buffer, ".", NULL);
      return program_path;
    }
    strfmt(possible_path,"%s/%s", PYTHON->buffer,  program_name);
    if ( (access(possible_path->buffer,F_OK) != -1) && (access(possible_path->buffer,X_OK) != -1) ) {
      strcopy(program_path,possible_path->buffer);
      //Remove tmp allocated variables
      RsatMemTracker = relem((void*)possible_path, RsatMemTracker);
      RsatMemTracker = relem((void*)cmd_which    , RsatMemTracker);
      if (verbose >= 4) RsatInfo("GetProgramPath() path found", program_path->buffer, ".", NULL);
      return program_path;
    }
    strfmt(possible_path,"%s/%s", SCRIPTS->buffer, program_name);
    if ( (access(possible_path->buffer,F_OK) != -1) && (access(possible_path->buffer,X_OK) != -1) ) {
      strcopy(program_path,possible_path->buffer);
      //Remove tmp allocated variables
      RsatMemTracker = relem((void*)possible_path, RsatMemTracker);
      RsatMemTracker = relem((void*)cmd_which    , RsatMemTracker);
      if (verbose >= 4) RsatInfo("GetProgramPath() path found", program_path->buffer, ".", NULL);
      return program_path;
    }
  }
  //If the path has not ben found yet, find the program anywhere in the user path
  if( strcmp(program_path->buffer,"") == 0 ) {
    strfmt(cmd_which, "which %s", program_name);
    //Execute cmd in cmdline and write obtained output to program_path
    if ( (fh_popen = popen(cmd_which->buffer,"r")) == NULL ) RsatFatalError("Unable to popen in GetProgramPath()",NULL);
    program_path->size = 0;
    while ( (program_path->buffer[program_path->size] = fgetc(fh_popen)) != '\n' ) {
      program_path->size++;
      strlimt(program_path);
    }
    //Add '\0'-end to string and +1 to size
    program_path->buffer[program_path->size] = '\0';
    program_path->size++;
    //Close filehandler
    pclose(fh_popen);
  }
  //Check if the program path has been found
  if (strcmp(program_path->buffer, "") == 0) {
    if (die_on_error) {
      RsatFatalError("The program", program_name, "is not found in your path.", NULL);
    } else {
      //Remove tmp allocated variables
      RsatMemTracker = relem((void*)possible_path, RsatMemTracker);
      RsatMemTracker = relem((void*)cmd_which    , RsatMemTracker);
      RsatWarning("The program", program_name, "is not found in your path.", NULL);
      return NULL;
    }
  }
  //Check if the program can be executed
  if( access(program_path->buffer,X_OK) != -1 ){
    if (die_on_error) {
      RsatFatalError("The program", program_name, "cannot be run.", NULL);
    } else {
      //Remove tmp allocated variables
      RsatMemTracker = relem((void*)possible_path, RsatMemTracker);
      RsatMemTracker = relem((void*)cmd_which    , RsatMemTracker);
      RsatWarning("The program", program_name, "cannot be run.", NULL);
      return NULL;
    }
  }
  //Remove tmp allocated variables
  RsatMemTracker = relem((void*)possible_path, RsatMemTracker);
  RsatMemTracker = relem((void*)cmd_which    , RsatMemTracker);

  if (verbose >= 4) RsatInfo("GetProgramPath() path found", program_path->buffer, ".", NULL);

  return program_path;
}

void strlistfree(stringlist *delete){
  if(!delete) return;
  strfree(delete->element);
  free(delete);
  return;
}

void strlistend(stringlist *group){
  stringlist *element = NULL;
  if(!group) return;

  do {
    element = group;
    //Pass to next list element
    group = group->next;
    //Release current element
    strlistfree(element);
  } while(group != NULL);

  if(verbose >= 14) RsatInfo("String list has been successfully deleted.", NULL);

  return;
}

stringlist *strlistnew(void){
  //Declare variable and allocate memory
  stringlist *new = NULL;
  new = (stringlist *)_malloc(sizeof(stringlist),"strlistnew");

  //Initialize attributes
  new->element = strnew();
  new->next    = NULL;

  return new;
}

stringlist *strlistnewToList(memstd **List){
  stringlist *new = NULL;
  new = strlistnew();
  MemTrackAdd((void*)new, STRLIST, List);
  return new;
}

//Returns the new added element, i.e. last element
stringlist *strlistadd(stringlist *group){
  //Generate new element
  stringlist *tmp = NULL;
  stringlist *new = NULL;
  new = strlistnew();

  //Create tmp
  tmp = group;
  //Attach variant to group
  do {
    if (tmp->next == NULL) {
      tmp->next = new;
      break;
    }
    tmp = tmp->next;
  } while(tmp != NULL );

  if(verbose >= 14) RsatInfo("String was added to list successfully",NULL);

  return new;
}

///////RANGES
void rangefree(range *delete){
  if(!delete) return;
  if(delete->var_info) varfree(delete->var_info);
  //if(delete->prev) rangefree(delete->prev);
  //if(delete->next) rangefree(delete->next);
  free(delete);
  return;
}

range *rangenew(void){
  //Declare variable and allocate memory
  range *new = NULL;
  new = (range *)_malloc(sizeof(range),"rangenew");

  //Initialize attributes
  new->start[0]       =        0;
  new->start[1]       =        0;
  new->end[0]         =        0;
  new->end[1]         =        0;
  new->length[0]      =        0;
  new->length[1]      =        0;
  new->left_flank[0]  =        0;
  new->left_flank[1]  =        0;
  new->right_flank[0] =        0;
  new->right_flank[1] =        0;
  new->var_info       = varnew();
  new->prev           =     NULL;
  new->next           =     NULL;

  return new;
}

range *rangenewToList(memstd **List){
  range *new = NULL;
  new = rangenew();
  MemTrackAdd((void*)new, RNG, List);
  return new;
}

range *rangeadd(range *group){
  //Generate new element
  range *tmp = NULL;
  range *new = NULL;
  new = rangenew();

  //Create tmp
  tmp = group;
  //Attach variant to group
  do {
    if (tmp->next == NULL) {
      tmp->next = new;
      new->prev = tmp;
      break;
    }
    tmp = tmp->next;
  } while(tmp != NULL );

  if(verbose >= 14) RsatInfo("Range was added to group successfully",NULL);

  return new;
}

range *rangefill(range *element,int start1,int start2,int end1,int end2,int length1,int length2,int left_flank1,int left_flank2,int right_flank1,int right_flank2){
  //Initialize attributes
  element->start[0]       = start1;
  element->start[1]       = start2;
  element->end[0]         = end1;
  element->end[1]         = end2;
  element->length[0]      = length1;
  element->length[1]      = length2;
  element->left_flank[0]  = left_flank1;
  element->left_flank[1]  = left_flank2;
  element->right_flank[0] = right_flank1;
  element->right_flank[1] = right_flank2;

  return element;
}

void rangend(range *group){
  range *element = NULL;
  if(!group) return;

  do {
    element = group;
    //Pass to next list element
    group = group->next;
    //Release current element
    rangefree(element);
  } while(group != NULL);

  if(verbose >= 14) RsatInfo("Range group has been successfully deleted.",NULL);

  return;
}

///////VARIANTS
void varfree(variant *delete){
  if(!delete) return;
  strfree(delete->chromosome);
  strfree(delete->start     );
  strfree(delete->end       );
  strfree(delete->id        );
  strfree(delete->SO        );
  strfree(delete->reference );
  strfree(delete->alleles   );
  strfree(delete->freq      );
  free(delete);
  return;
}

variant *varnew(void){
  //Declare variable and allocate memory
  variant *new = NULL;
  new = (variant *)_malloc(sizeof(variant),"varnew");

  //Initialize attributes
  new->chromosome = strnew();
  new->start      = strnew();
  new->end        = strnew();
  new->strand[0]   =  '+';
  new->strand[1]   = '\0';
  new->id         = strnew();
  new->SO         = strnew();
  new->reference  = strnew();
  new->alleles    = strnew();
  new->freq       = strnew();
  new->prev       = NULL;
  new->next       = NULL;

  return new;
}

variant *varnewToList(memstd **List){
  variant *new = NULL;
  new = varnew();
  MemTrackAdd((void*)new, VAR, List);
  return new;
}

//Returns the new added element, i.e. last element
variant *varadd(variant *group){
  //Generate new element
  variant *tmp = NULL;
  variant *new = NULL;
  new = varnew();

  //Create tmp
  tmp = group;
  //Attach variant to group
  do {
    if (tmp->next == NULL) {
      tmp->next = new;
      new->prev = tmp;
      break;
    }
    tmp = tmp->next;
  } while(tmp != NULL );

  if(verbose >= 14) RsatInfo("Variant was added to group successfully",NULL);

  return new;
}

variant *varfill(variant *element,char *chrInfo, char *startInfo, char *endInfo, char *strandInfo, char *idInfo, char *SOInfo, char *refInfo, char *allelesInfo, char *freqInfo){
  //Initialize attributes
  strcopy(element->chromosome,    chrInfo);
  strcopy(element->start     ,  startInfo);
  strcopy(element->end       ,    endInfo);
  if( strandInfo[0] == '-' ) strcpy(element->strand, "-");
  strcopy(element->id       ,      idInfo);
  strcopy(element->SO       ,      SOInfo);
  strcopy(element->reference,     refInfo);
  strcopy(element->alleles  , allelesInfo);
  strcopy(element->freq     ,    freqInfo);

  return element;
}

void varend(variant *group){
  variant *element = NULL;
  if(!group) return;

  do {
    element = group;
    //Pass to next list element
    group = group->next;
    //Release current element
    varfree(element);
  } while(group != NULL);

  if(verbose >= 14) RsatInfo("Variants group has been successfully deleted.",NULL);

  return;
}

void cgiMessage(char *message_type,char *color,char *fmt,va_list ap){
  if(!message_type) message_type = "Information";
  if(!color) color = "#006600";

  fprintf(stderr, "<blockquote class='%s'>\n<font color='%s'><b> %s : </b> ", message_type,color,message_type);
  for (char *str = fmt ; str != NULL ; str = va_arg(ap,char*)) {
    fprintf(stderr,"%s<br>\n",str);
  }
  fprintf(stderr, "</font>\n</blockquote> <br><hr size=3>\n");

  //va_end(ap);

  return;
}

void cgiWarning(char *fmt,va_list ap){
  //strcat(message,"<P>");
  cgiMessage("Warning","#FFAA00",fmt,ap);
  return;
}

void RsatInfo(char *fmt, ...) {
  va_list ap;
  va_start(ap,fmt);

  if(getenv("RSA_OUTPUT_CONTEXT") && strcmp(getenv("RSA_OUTPUT_CONTEXT"),"cgi") == 0){
    cgiMessage(NULL,NULL,fmt,ap);
  }else{
    fprintf(stderr,";INFO\t");
    for (char *str = fmt ; str != NULL ; str = va_arg(ap,char*)) {
      fprintf(stderr," %s",str);
    }
    fprintf(stderr,"\n");
  }

  va_end(ap);
  return;
}

void RsatWarning(char *fmt, ...) {
  va_list ap;
  va_start(ap,fmt);

  if(getenv("RSA_OUTPUT_CONTEXT") && strcmp(getenv("RSA_OUTPUT_CONTEXT"),"cgi") == 0){
    cgiWarning(fmt,ap);
  }else{
    fprintf(stderr,";WARNING\t");
    for (char *str = fmt ; str != NULL ; str = va_arg(ap,char*)) {
      fprintf(stderr," %s",str);
    }
    fprintf(stderr,"\n");
  }

  va_end(ap);
  return;
}

void cgiError(char *fmt,va_list ap){
  char *site = NULL;
  char *admi = NULL;
  char hostname[MAX_HOSTNAME];

  if(gethostname(hostname,MAX_HOSTNAME) != 0) strcpy(hostname,"User");
  //printf("Hello1 %s\n",hostname );

  site = getenv("rsat_site");
  admi = getenv("rsat_server_admin");

  if (!site) site = "rsat_site";
  if (!admi) admi = "rsat_admin";

  fprintf(stderr, "<blockquote class='%s'>\n<font color='%s'><b> %s : </b> ", "Error","#DD0000","Error");
  fprintf(stderr,"Error ocurred on RSAT site: %s; host server: %s; admin: %s <br>\n",site,hostname,admi);
  fprintf(stderr, "</font>\n</blockquote> <br><hr size=3>\n");

  cgiMessage("Error","#DD0000",fmt,ap);

  fprintf(stderr,"</body></html>\n");

  //Remove all memory pile tracers
  //rlist(RsatMemTracker);

  //exit(0);

  return;
}

void RsatFatalError(char *fmt, ...) {
  char *context = NULL;

  va_list ap;
  va_start(ap,fmt);
  context = getenv("RSA_OUTPUT_CONTEXT") ? getenv("RSA_OUTPUT_CONTEXT") : "screen";
  if (setenv("RSA_ERROR","1",1) == -1) {
    fprintf(stderr, "Error:\n\tUnable to set RSA_ERROR=1 in RsatFatalError()\n");
  }
  if(strcmp(context,"cgi") == 0) {
     cgiError(fmt,ap);
  } else {
    fprintf(stderr,"Error\n\t");
    for (char *str = fmt ; str != NULL ; str = va_arg(ap,char*)) {
      fprintf(stderr,"%s ",str);
    }
    fprintf(stderr,"\n");
  }

  va_end(ap);

  //Remove all memory pile tracers
  rlist(RsatMemTracker);

  exit(0);

  return;
}

void *_malloc(size_t size,char *func){
  void *pointer = NULL;
  char bytes[BASE_STR_LEN];

  sprintf(bytes,"%zd",size);
  if(verbose >= 14) RsatInfo("Allocating",bytes,"bytes at",func,NULL);

  pointer = malloc(size);

  if(pointer == NULL){
    RsatFatalError("Unable to allocate",bytes,"bytes at",func,NULL);
  }

  return pointer;
}

void *_realloc(void *ptr,size_t size,char *func){
  void *pointer = NULL;
  char bytes[BASE_STR_LEN];

  sprintf(bytes,"%zd",size);
  if(verbose >= 14) RsatInfo("Reallocating",bytes,"bytes at",func,NULL);

  pointer = realloc(ptr,size);

  if(pointer == NULL){
    RsatFatalError("Unable to allocate",bytes,"bytes at",func,NULL);
  }
  return pointer;
}

TRIE *TrieStart(void){
  //Declare variables
  TRIE *node = NULL;
  //Allocate memory for trie
  //TODO.WSG Update node = (TRIE *)_stdmalloc(sizeof(TRIE),"TrieStart");
  node = (TRIE *)_malloc(sizeof(TRIE),"TrieStart");
  MemTrackAdd((void*)node, TRI, &RsatMemTracker);

  //node = (TRIE *)malloc(sizeof(TRIE));
  //if(!(node)) RsatFatalError("Unable to allocate memory for a node",NULL);

  //Initialize node
  node->leaf = 0;
  node->info = NULL;
  for (int i = 0; i < ALPHABET_SIZE; i++) {
    node->children[i] = NULL;
  }

  return node;
}

/*Creates a trie node, with all children and info initialized with NULL.
  Leaf attribute is set to 0. It does not receive parameters, and a pointer
  to the node is returned. On failure a FatalError is raised.*/
TRIE *TrieNew(void){
  //Declare variables
  TRIE *node = NULL;
  //Allocate memory for trie
  node = (TRIE *)_malloc(sizeof(TRIE),"TrieNew");
  //node = (TRIE *)malloc(sizeof(TRIE));
  //if(!(node)) RsatFatalError("Unable to allocate memory for a node",NULL);

  //Initialize node
  node->leaf = 0;
  node->info = NULL;
  for (int i = 0; i < ALPHABET_SIZE; i++) {
    node->children[i] = NULL;
  }

  return node;
}

/*Inserts a pair of key-value to TRIE. It receives a TRIE pointer
  (where the key-value is going to be inserted), a char pointer
  to key name and a char pointer to the value pair. If the key
  already exists or the char letter is not supported on the
  current alphabet a warning is raised. On failure a FatalError
  is raised*/
void TrieInsert(TRIE *trie,char *key,char *value){
  //Declare variables
  TRIE *tmp = trie;
  int value_length = strlen(value);
  int charIndex = 0;

  //Iterate through every key letter and insert them to trie
  for (int i = 0; key[i] != '\0' ; i++) {
    //Convert key letter to alphabet index
    charIndex = GetIndex(key[i]);
    //Test if index is a valid alphabet index
    if(charIndex < 0 || charIndex > 93){
      RsatWarning("Unvalid character for this Alphabet, key not inserted",NULL);
      return;
    }
    //If letter has no entry at trie add it
    if (! (tmp->children[charIndex]) ){
      tmp->children[charIndex] = TrieNew();
    }
    //Pass to next trie letter
    tmp = tmp->children[charIndex];
  }
  //Test if last trie letter is a leaf, i.e. if key already exists
  if (tmp->leaf) {
    RsatWarning("This key already,exists,cannot duplicate entries",NULL); //Erase structure?
    return;
  //Change leaf status to identify letter as end of key string
  } else {
    tmp->leaf = 1;
  }
  //Allocate enough memory for value and add it to tmp->info
  tmp->info = (char *)_malloc( sizeof(char) * (value_length+1) , "TrieInsert" );
  //if( (tmp->info = (char *)malloc(sizeof(char) * (value_length+1))) == NULL){
  //  RsatFatalError("Unable to allocate memory for value on leaf",NULL);
  //}
  //Copy value content to leaf trie
  strcpy(tmp->info,value);

  return;
}

/*Searches a key in a trie. It returns a char pointer to the value
  content if the key was found or a NULL otherwise*/
char *TrieSearch(TRIE *trie,char *key){
  //Declare variables
  TRIE *tmp = trie;
  int charIndex = 0;

  //Iterate through every key letter
  for (int i = 0; key[i] != '\0' ; i++) {
    charIndex = GetIndex(key[i]);
    //If a null entry is found, i.e. a key letter is missing
    //return null
    if( !(tmp->children[charIndex]) ){
      return NULL;
    }

    tmp = tmp->children[charIndex];
  }
  //Return pointer to value if last key letter on trie is a leaf
  //i.e. is a legitimate trie
  return (tmp->leaf)? tmp->info : NULL;
}

/*Removes entire allocated memory for a given trie. It does not
  return any value*/
void TrieEnd(TRIE *trie){
  //If trie node info is not empty,i.e. has a value, free it
  if(trie->info != NULL) free(trie->info);
  //Iterate through every non-null children node recursively
  for (int i = 0; i < ALPHABET_SIZE; i++) {
    if(trie->children[i] != NULL){
      TrieEnd(trie->children[i]);
    }
  }
  //Free node
  free(trie);
  return;
}

//STR

string *strnew(void){
  string *new = NULL;
  //NOTE.(2017-04-10) By doing this, _strmalloc
  //is completely isolated. Should I remove it?
  //new = _strmalloc(sizeof(string),"strnew");
  new = (string *)_malloc(sizeof(string),"strnew");
  new->buffer = (char *)_malloc(sizeof(char) * BASE_STR_LEN,"strnew");
  new->length = BASE_STR_LEN;
  new->size   = 0;

  return new;
}

/*NOTE WSG(2017-08-10). Updates the content of the list by directly
  accessing to the pointer contents of list. Be CAREFULL for thread
  safeties options. This is why I use **list and not *list. */
string *strnewToList(memstd **List){
  string *new = NULL;
  new = strnew();
  MemTrackAdd((void*)new, STR, List);
  //ListStrAdd(new,List);
  return new;
}

char *strccat(string *destn,char *fmt, ...){
  int nchar;
  va_list ap;

  //Safety chek that this a '\0'-terminated string.
  if(destn->buffer[destn->size - 1] != '\0' || fmt == NULL) {
    printf("\n\n\nUnable to concatenate correctly\n\n" );
    return NULL;
  }

  //Format string with va_list where buffer == '\0'
  va_start(ap,fmt);
  nchar = vsnprintf(destn->buffer + destn->size - 1,destn->length - (destn->size - 1),fmt,ap);
  va_end(ap);

  //If nchar was not successfull,retry,now with enough space in buffer
  if(nchar + 1 > destn->length - (destn->size - 1)){
    destn->buffer = stralloc(destn,destn->size + nchar);
    va_start(ap,fmt);
    nchar = vsnprintf(destn->buffer + destn->size - 1,destn->length - (destn->size - 1),fmt,ap);
    va_end(ap);
  }

  //Update size
  destn->size += nchar;

  return destn->buffer;
}
/*Copies a '\0'-termianted char array to a string struct*/
char *strcopy(string *destn,char *orign){
  int i;

  if(orign == NULL) return NULL;
  destn->size = 0;

  for (i = 0; orign[i] != '\0'; i++) {

    if(destn->size >= destn->length){
      destn->buffer = stralloc(destn,destn->size + BASE_STR_LEN/2);
    }
    destn->buffer[i]= orign[i];
    destn->size++;
  }

  if(destn->size >= destn->length){
    destn->buffer = stralloc(destn,destn->size + BASE_STR_LEN/2);
  }

  destn->buffer[i] = '\0';
  destn->size++;

  return destn->buffer;
}

char *strfmt(string *strn, char *fmt, ...){
  int nchar;
  va_list ap;

  //Safety chek that this at least and empty string, i.e.'\0'-terminated.
  if(fmt == NULL) return NULL;

  //Format string with va_list where buffer == '\0'
  va_start(ap,fmt);
  strn->size = 0;
  nchar = vsnprintf(strn->buffer,strn->length,fmt,ap);
  va_end(ap);
  //If nchar was not successfull,retry,now with enough space in buffer, add +1 for \0
  if(nchar + 1 > strn->length){
    strn->buffer = stralloc(strn,nchar + 1);
    va_start(ap,fmt);
    nchar = vsnprintf(strn->buffer,strn->length,fmt,ap) ;
    va_end(ap);
  }
  //Update size
  strn->size = nchar + 1;

  return strn->buffer;
}

void strfree(string *delete){
  if(!delete) return;
  free(delete->buffer);
  free(delete);
  return;
}

string *_strmalloc(size_t size,char *func){
  string *pointer = NULL;
  pointer = (string *)_malloc(size,func);
  //NOTE(2017-04-10).I passed this line to a more
  //generic function.QUESTION. Without this line,
  //perhaps I can suppress this function?
  //ListStrAdd(pointer,&ListStr);
  return pointer;
}
/*Test if a string buffer has reached its max length. If true,
  it reallocs new memory by an adding factor of +BASE_STR_LEN/2.*/
string *strlimt(string *totest){
  if (totest->size >= totest->length ) {
   totest->buffer = stralloc(totest, (int)(totest->length + BASE_STR_LEN/2));
 }
 return totest;
}

/*Initializes all token to first address,i.e. token[0]*/
char **initokadd(string *line,char **token,int numtok){
  for (int i = 0; i < numtok; i++) {
    token[i] = token[0];
  }
  return token;
}


//Allocates memory for an array of tokens i.e. char*,
//It returns a pointer for the start of the array.
char **getokens(int numtok){
  char **tokens = NULL;

  //Allocate memory for token array
  tokens = (char**)_MemTrackMalloc(sizeof(char*) * numtok, &RsatMemTracker, "getokens");
  //Initialize all tokens in NULL
  for (int i = 0; i < numtok; i++) {
    tokens[i] = NULL;
  }

  return tokens;
}

/*It reallocs the buffer of the given string with a new
  given size.The new total length of the string is the
  size parameter passed to the function.It returns the
  pointer to the new strings buffer*/
char *stralloc(string *resize,int size){
  resize->buffer = (char *)_realloc(resize->buffer,sizeof(char) * size,"stralloc");
  resize->length = size;
  return resize->buffer;
}

void *_MemTrackMalloc(size_t size, memstd **list, char *func){
  void *pointer = NULL;
  pointer = _malloc(size,func);
  MemTrackAdd(pointer, GEN, list);
  return pointer;
}

/*Creates an element entry for a backtracking memory list.*/
memstd *MemTrackElmt(void){
  //Request memory
  memstd *list = NULL;
  list = (memstd *)_malloc(sizeof(memstd),"MemTrackElmt");
  //Initialize element
  list->id   = 0;
  list->mem  = NULL;
  list->next = NULL;
  return list;
}

/*This function creates a pile for backtracking allocated
  memory.Everything but the id is Initialized to NULL*/
memstd *MemTrackNew(void){
  //Request memory
  memstd *list = NULL;
  list = (memstd *)_malloc(sizeof(memstd),"MemTrackNew");
  //Initialize element
  list->id   = 0;
  list->mem  = NULL;
  list->next = NULL;

  return list;
}

/*NOTE WSG(2017-08-10). Updates the content of the list by directly
  accessing to the pointer contents of list. Be CAREFULL for thread
  safety options for a future paralelization. This is why I use
  **list and not *list. */
/*This function creates adds an entry to the backtracking memory
  allocated pile. A malloc pointer is needed to add it to new->mem */
memstd *MemTrackAdd(void *ptr, int type, memstd **list){
  //Generate new element
  memstd *new = NULL;
  new = MemTrackElmt();
  //Add node information
  new->next = *list;
  new->mem  = ptr;
  new->id   = type;
  //Add node to list
  *list = new;
  //NOTE WSG. Remove printf("ListStdAdd ADDED; pointer: %p, id : %d, mem : %p ,next->id: %d\n", new, new->id, new->mem, new->next->id );
  if(verbose >= 14) RsatInfo("Pointer element was added to list",NULL);

  return *list;
}

memstd *MemTrackUpd(void *old, void *new, memstd *list){
  memstd *start = NULL;

  //Keep registry of the first struct element of the
  //pile in order to return it at functions return
  start = list;

  do {
    if(list->mem == old) {
      list->mem = new;
      if(verbose >= 14) RsatInfo("Found memory address in MemTrackUpd. Updating.",NULL);
      //Return the first element of the pile
      return start;

    }
    //Pass to the next struct element to test
    list  = list->next;
  } while(list != NULL);

  if(verbose >= 14) RsatInfo("Unable to update memory address. Not found.",NULL);

  //Return the first element of the pile
  return start;
}


memstd *ListStdUpd(void *lastptr,void *currptr,memstd *list){
  while (list->next != NULL) {
    if(list->mem == lastptr){
      list->mem = currptr;
      return list;
    }
    list = list->next;
  }
  return NULL;
}

void rlist(memstd *list){
  memstd *entry = NULL;
  if(!list) return;

  //TODO WSG. Remove
  /*memstd *listdbg = NULL;
  listdbg = list;
  do {
      printf("This is rstDlist %p list \n", listdbg);
      printf("This is rstDlist %d id \n", listdbg->id);
      printf("This is rstDlist %p mem \n", listdbg->mem);
      printf("This is rstDlist %p next \n\n", listdbg->next);

      listdbg = listdbg->next;
  } while(listdbg != NULL);
  printf("1.-END of PRINTING inside rstDlist\n" );*/

  do {
    //TODO WSG. Remove printf("1.-START of DELETING inside rstDlist\n" );
    entry = list;

    //Release contained pointer
    if (list->id == GEN) {
      free(list->mem);
    } else if (list->id == STR) {
      strfree( (string*)list->mem );
    } else if (list->id == TRI) {
      TrieEnd( (TRIE*)list->mem );
    } else if (list->id == VAR) {
      varend( (variant*)list->mem );
    } else if (list->id == STRLIST) {
      strlistend((stringlist*)list->mem);
    } else if (list->id == SITE){
      sitefree((site*)list->mem);
    } else if (list->id == SCAN){
      scanend((scan*)list->mem);
    } else if (list->id == VSCAN){
      varscanend((varscan*)list->mem);
    } else if (list->id == RNG){
      rangend((range*)list->mem);
    } else {
      RsatInfo("This is not a valid memory object to free", NULL);
    }
    //Pass to next list element
    list = list->next;
    //Release current element
    free(entry);

  } while(list != NULL);
  //TODO WSG. Remove printf("2.-End of DELETING inside rstDlist\n" );

  if(verbose >= 14) RsatInfo("memstd list has been successfully deleted.",NULL);
  return;
}

memstd *relem(void *ptr,memstd *list){
  memstd *start = NULL;
  memstd *entry = NULL;

  //NOTE. WSG Remove printf("\nrstDelem ; This is ListStd pointer: %p, id : %d, mem : %p ,next->id: %d\n", list, list->id, list->mem, list->next->id );
  //NOTE WSG.Until I find a safe way to print pointer without a priori
  //knowledge of its length I will not print it.
  //char memory_address[BASE_STR_LEN];
  //sprintf(memory_address,"%p",ptr);

  //Keep registry of the first struct element of the
  //pile in order to return it at functions return
  start = list;

  do {
    if(list->mem == ptr){
      //if(verbose >= 14) RsatInfo("Found memory address",memory_address,". Deleting.",NULL);
      //NOTE. WSG printf("rstDelem FOUND; pointer: %p, id : %d, mem : %p ,next->id: %d\n", list, list->id, list->mem, list->next->id );
      if(verbose >= 14) RsatInfo("Found memory address in relem. Deleting.",NULL);
      if(start == list){
        //New start is the penultimate element of pile
        start = list->next;
      } else {
        //Link elements of pile, prev and next struct
        //of the desired element to delete
        entry->next = list->next;
      }
      //Release memstd elem and contained pointer
      if (list->id == GEN) {
        free(list->mem);
      } else if (list->id == STR) {
        strfree( (string*)list->mem );
      } else if (list->id == TRI) {
        TrieEnd( (TRIE*)list->mem );
      } else if (list->id == VAR) {
        varend( (variant*)list->mem );
      } else if (list->id == STRLIST){
        strlistend((stringlist*)list->mem);
      } else if (list->id == SITE){
        sitefree((site*)list->mem);
      } else if (list->id == SCAN){
        scanend((scan*)list->mem);
      } else if (list->id == VSCAN){
        varscanend((varscan*)list->mem);
      } else if (list->id == RNG){
        rangend((range*)list->mem);
      } else {
        RsatInfo("This is not a valid memory object to free", NULL);
      }
      free(list);
      //Return the first element of the pile
      return start;

    }
    //Keep registry of the previous already tested struct
    //in order to relink if found in-between pile.
    entry = list;
    //Pass to the next struct element to test
    list  = list->next;
  } while(list != NULL);

  if(verbose >= 14) RsatInfo("Unable to delete memory address. Not found.",NULL);

  //Return the first element of the pile
  return start;
}

/*
  NOTE.WSG (2017-14-06). I should come back and recheck a better way to initialize
  both RsatMemTracker. When it is created with MemTrackNew() the first element at
  "mem" is NULL so I ask for memory manually and asigned it to the first element.
  The first element is RSAT environmental string variable, and the second
  is getpid(),i.e.unsigned int (ALWAYS!).

  NOTE.WSG (2017-14-06). For some reason valgrind is not aware that fh are being closed
  at exit(). However this is not a problem at all. Just double-check.
*/
void InitMemLists(void){
  //Create list for memory collection
  RsatMemTracker = MemTrackNew();

  //Add content to RSAT variable! i.e. get "RSAT".This is the
  //first element in RsatMemTracker.
  RSAT = strnew();
  RsatMemTracker->id  = STR;
  RsatMemTracker->mem = (void*)RSAT;
  if( !strcopy(RSAT,getenv("RSAT")) ) RsatFatalError("$RSAT environmental variable is not defined,please use 'source [PATHTORSAT]/rsat/RSAT_config.bashrc'",NULL);

  //Add content to PID variable!
  PID          = (unsigned int *)_MemTrackMalloc(sizeof(unsigned int), &RsatMemTracker,"InitMemLists");
  //PID          = (unsigned int *)_malloc(sizeof(unsigned int),"InitRSAT");
  //ListStd->mem = (void*)PID;
  *PID         = getpid();

  //Reserve memory for variables PROGRAM and CMD
  CMD     = strnewToList(&RsatMemTracker);
  PROGRAM = strnewToList(&RsatMemTracker);

  return;
}

void ReadProperties(void){
  //Check if $RSAT environmental variable is contained in string-object RSAT, else return
  if (RSAT != NULL && strcmp(RSAT->buffer,"") == 0){
    RsatWarning("$RSAT environmental variable is not defined,please use 'source [PATHTORSAT]/rsat/RSAT_config.bashrc",NULL);
    return;
  }
  //Declare variables
  FILE *fh_RSAT_configProps = NULL;
  string *RSAT_configProps  = NULL;
  string *line              = NULL;

  int j = 0;
  char **token = NULL;

  //Allocate memory for variables
  RSAT_configProps  = strnewToList(&RsatMemTracker);
  line              = strnewToList(&RsatMemTracker);
  token             = getokens(2);

  //Get RSAT_config.props file path and test for existance
  strfmt(RSAT_configProps,"%s/RSAT_config.props",RSAT->buffer);
  fh_RSAT_configProps = OpenInputFile(fh_RSAT_configProps,RSAT_configProps->buffer);
  if(verbose >= 5) RsatInfo("Reading property file",RSAT_configProps->buffer,NULL);

  token[0] = line->buffer;

  while ( fread((line->buffer + line->size),1,1,fh_RSAT_configProps) == 1) {
    //If '=',i.e. separator, is found report next address as token[1]
    if (line->buffer[line->size] == '=') {
      token[++j] = line->buffer + line->size + 1;
      line->buffer[line->size] = '\0';
    }

    //If end of line is found,proceed to export key-value pair for
    //environmental variables i.e. token[0] and token[1] respectively.
    if (line->buffer[line->size] == '\n') {

      //Reinitialize counters
      line->buffer[line->size] = '\0';
      line->size = 0;
      j = 0;

      //Filters for comments and more
      if(token[0][0] == '#'){initokadd(line,token,2);continue;}
      if(token[0][0] == ';'){initokadd(line,token,2);continue;}
      if(token[0][0] == ' '){initokadd(line,token,2);continue;}
      if(token[0][0] == '\t'){initokadd(line,token,2);continue;}
      if(token[0][0] == '\0'){initokadd(line,token,2);continue;}

      //Export this pair of key-value as
      //environmental variables.
      if(verbose >= 12) RsatInfo("Setting environmental variable",token[0],"=",token[1],NULL);
      if (setenv(token[0],token[1],1) == -1){
        RsatFatalError("Unable to set ",token[0],"= ",token[1]," in InitRSAT()",NULL);
      }
      //If successful continue to the next start of line for reading
      initokadd(line,token,2);
      continue;
    }

    //Resize line string if limit has reached
    limlinetok(line,token,2);

    //Keep track of size count for each char read
    line->size++;

    //If current size of buffer is about to reach its maximum length
    //realloc current buffer into a much bigger one.
    //if(line->length - line->size == 1){ //NOTE.WSG. Should I constraint this to equal '0' as a more effective way?
    //  stralloc( line,(int)(line->length + BASE_STR_LEN/4) );
      //printf("This is token[0]: %p ;token[1] :%p; length: %d \n",token[0],token[1],(int)(token[1] - token[0])); //NOTE.Assess.
    //  token[1] = line->buffer + (token[1] - token[0] - 1);
    //  token[0] = line->buffer;
    //}

  }

  //Close filehandler
  fclose(fh_RSAT_configProps);

  //Remove tmp allocated variables
  RsatMemTracker = relem( (void*)line,RsatMemTracker );
  RsatMemTracker = relem( (void*)RSAT_configProps,RsatMemTracker );
  RsatMemTracker = relem( (void*)token,RsatMemTracker );

  return;
}

time_t InitRSAT(char *program,char *cmd){
  struct passwd *pwd;
  struct tm *info;

  time_t start_time;

  char   *login      = NULL;
  char   *rsat_site  = NULL;

  string *rsat_tmp   = NULL;
  string *ip_address = NULL;

  //Set umask
  umask(0022);

  //ReadProperties from RSAT_config.props
  ReadProperties();

  //Get start_time for this initialization
  if( time(&start_time) == -1 ) RsatFatalError("Time is not available in StartScript()", NULL);
  info = localtime(&start_time);

  //Create tmp-string variables
  rsat_tmp   = strnewToList(&RsatMemTracker);
  ip_address = strnewToList(&RsatMemTracker);

  //Check if environmental variable RSAT_WWW is defined
  if(!(getenv("rsat_www")) || (strcmp(getenv("rsat_www"),"auto") == 0)){
    if( !strcopy(ip_address,getenv("HTTP_HOST")) ) strcopy(ip_address,"localhost");
    if(verbose >= 4) RsatInfo("Host IP address",ip_address->buffer,NULL);
    strfmt(rsat_tmp,"http://%s/rsat/",ip_address->buffer);
    if (setenv("rsat_www",rsat_tmp->buffer,1) == -1){
      RsatFatalError("Unable to set rsat_www=",rsat_tmp->buffer,"in InitRSAT()",NULL);
    }
    if(setenv("rsat_ws",rsat_tmp->buffer,1) == -1){
      RsatFatalError("Unable to set rsat_ws=",rsat_tmp->buffer,"in InitRSAT()",NULL);
    }
    if(verbose >= 4) RsatInfo("RSAT_WWW",getenv("rsat_www"),NULL);
  }

  //Directories variables
  BIN     = strnewToList(&RsatMemTracker);
  LIB     = strnewToList(&RsatMemTracker);
  PYTHON  = strnewToList(&RsatMemTracker);
  SCRIPTS = strnewToList(&RsatMemTracker);

  strfmt(BIN,"%s/bin",RSAT->buffer);
  strfmt(LIB,"%s/ext_lib",RSAT->buffer);
  strfmt(PYTHON,"%s/python-scripts",RSAT->buffer);
  strfmt(SCRIPTS,"%s/perl-scripts",RSAT->buffer);

  //Get memory for RSAT variables
  HTML         = strnewToList(&RsatMemTracker);
  LOGS         = strnewToList(&RsatMemTracker);
  WWW_TMP      = strnewToList(&RsatMemTracker);
  counter_file = strnewToList(&RsatMemTracker);

  date    = strnewToList(&RsatMemTracker);
  HOSTNAME= strnewToList(&RsatMemTracker);

  LOGIN   = strnewToList(&RsatMemTracker);

  //NOTE. Delete as soon as InitMemLists() works
  //CMD     = strnewToList();
  //PROGRAM = strnewToList();

  log_file               = strnewToList(&RsatMemTracker);
  exec_time_log_file     = strnewToList(&RsatMemTracker);
  start_time_log_file    = strnewToList(&RsatMemTracker);
  web_attacks_log_file   = strnewToList(&RsatMemTracker);
  denied_access_log_file = strnewToList(&RsatMemTracker);

  //Add content to RSAT variables
  if( !strfmt(WWW_TMP,"%s/tmp",getenv("rsat_www")) ) strcopy(WWW_TMP,"");
  strfmt(HTML,"%s/public_html",RSAT->buffer);
  strfmt(LOGS,"%s/logs",RSAT->buffer);
  strfmt(counter_file,"%s/count-file",LOGS->buffer);

  //CHANGED WSG Start date
  AlphaDate(date,&start_time);
  //date->size = strlen(date->buffer) + 1;

  //Hostname size sticks to standard SUSv2(256),which is bigger
  //than the 64 POSIX standards specifics.
  stralloc(HOSTNAME,MAX_HOSTNAME);
  if(gethostname(HOSTNAME->buffer,MAX_HOSTNAME) != 0) strcpy(HOSTNAME->buffer,"User");
  HOSTNAME->size = strlen(HOSTNAME->buffer) + 1;

  //Get Login
  login = getlogin();
  if(!login){                                  //If failed,getpwuid()
    pwd = getpwuid(getuid());
    login = pwd->pw_name;
    strcopy(LOGIN,(!login) ? "Kilroy" : login);//If failed again,use "Kilroy"
  } else {
    strcopy(LOGIN,login);
  }
  //Associated program variables
  //strcopy(CMD,cmd);//NOTE.This is redundant in main
  //strcopy(PROGRAM,program);//NOTE.This is redundant in main

  //Associated variables to rsat_site
  if((rsat_site = getenv("rsat_site"))){
    strfmt(log_file              ,"%s/log-file_%s_%04d_%02d.txt",LOGS->buffer,rsat_site,1900 + info->tm_year,1 + info->tm_mon);
    strfmt(exec_time_log_file    ,"%s/exec_time_log_%s_%04d_%02d.txt",LOGS->buffer,rsat_site,1900 + info->tm_year,1 + info->tm_mon);
    strfmt(start_time_log_file   ,"%s/start_time_log_%s_%04d_%02d.txt",LOGS->buffer,rsat_site,1900 + info->tm_year,1 + info->tm_mon);
    strfmt(web_attacks_log_file  ,"%s/web_attacks_log_%s_%04d_%02d.txt",LOGS->buffer,rsat_site,1900 + info->tm_year,1 + info->tm_mon);
    strfmt(denied_access_log_file,"%s/denied_access_log_%s_%04d_%02d.txt",LOGS->buffer,rsat_site,1900 + info->tm_year,1 + info->tm_mon);
  }else{
    strcopy(log_file              ,"");
    strcopy(exec_time_log_file    ,"");
    strcopy(start_time_log_file   ,"");
    strcopy(web_attacks_log_file  ,"");
    strcopy(denied_access_log_file,"");
  }

  //Remove tmp strings RSAT_configProps and line
  RsatMemTracker = relem( (void*)ip_address,RsatMemTracker );
  RsatMemTracker = relem( (void*)rsat_tmp  ,RsatMemTracker );

  //TODO WSG.Remove DBG print
  //printMemstr(RsatMemTracker);

  return start_time;
}

int UpdateExecTimeLogFile(time_t start_time,time_t done_time){
  FILE *stime_log_file_fh = NULL;

  char *remote_addr = NULL;
  string *end_tstr  = NULL;

  //Allocate memory for variables
  end_tstr = strnewToList(&RsatMemTracker);

  //Get end process date
  AlphaDate(end_tstr,&done_time);

  //Get remote address i.e. Web queries
  remote_addr = getenv("REMOTE_ADDR");
  if(!remote_addr) remote_addr = "";

  //Write header of the exect time log file if file does not exist
  if (access(exec_time_log_file->buffer,F_OK) == -1) {
    if(verbose >= 4) RsatInfo("Creating execution time log file",exec_time_log_file->buffer,NULL);
    stime_log_file_fh = OpenOutputFile(stime_log_file_fh,exec_time_log_file->buffer);
    fprintf(stime_log_file_fh,"#start_date.time\tdone_date.time\tseconds\tPID\thostname\tusername\tscript_name\tcommand\tremote_addr\n");
    fclose(stime_log_file_fh);
    chmod(exec_time_log_file->buffer,0666);
  }

  //Update process information to log file
  if(verbose >= 4) RsatInfo("Updating execution time log file",exec_time_log_file->buffer,NULL);
  stime_log_file_fh = OpenAppendFile(stime_log_file_fh,exec_time_log_file->buffer);
  fprintf(stime_log_file_fh, "%s\t%s\t%.2f\t%d\t%s\t%s\t%s\t%s\t%s\n",date->buffer,end_tstr->buffer,difftime(done_time,start_time),*PID,HOSTNAME->buffer,LOGIN->buffer,PROGRAM->buffer,CMD->buffer,remote_addr);
  fclose(stime_log_file_fh);
  chmod(exec_time_log_file->buffer,0666);

  //Remove tmp string end_tstr
  RsatMemTracker = relem( (void*)end_tstr,RsatMemTracker );

  return 1;
}

int CheckValOpt(char *value[]) {
  if (strncmp(value[1],"-",1) == 0) {
    RsatFatalError("Missing value for option",value[0],NULL);
  }
  return 1;
}
/*
  NOTE.WSG (2017-14-06). As soon as I debug the rest of the code. I should return to this function
  erase printfs and make sure that errors are being reported accurately.
*/
int doit(char *command,int dry,int die_on_error,int verbose,int batch,char *job_prefix,FILE *log_handle,FILE *err_handle,char *cluster_queue){
  int error;
  //printf("1.-Hello inside error on doit()\n" );

  string *error_message = NULL;

  if(log_handle){
    fprintf(log_handle,"\n\n%s\n",command);
  }

  if(!command){
    RsatWarning("doit() was called with empty command.",NULL);
    return 0;
  }
  //printf("1.-Hello inside error on doit()\n" );
  error = system(command);
  //printf("2.-Hello inside error on doit()\n" );

  if (error != 0){
    error_message = strnewToList(&RsatMemTracker);
    if (error == -1){ strfmt(error_message,"Could not execute the command\n\t%s",command);
    } else {
      strfmt(error_message,"Error\t%d\tocurred during execution of the command:\n\t%s",error,command);
    }


    //char *error_message   = NULL;
    //if( (error_message    = (char *)malloc(sizeof(char) * BASE_STR_LEN)) == NULL) RsatFatalError("Unable to allocate memory in make_temp_file()",NULL);

    //if(error == -1) sprintf(error_message,"Could not execute the command\n\t%s",command);
    //else{
    //  sprintf(error_message,"Error\t%d\tocurred during execution of the command:\n\t%s",error,command);
    //}

    if(err_handle) fprintf(err_handle, "\n\n%s\n",error_message->buffer);

    if(die_on_error){ RsatFatalError(error_message->buffer,NULL);
    } else {
      RsatWarning(error_message->buffer,NULL);
      RsatMemTracker = relem( (void*)error_message,RsatMemTracker );
      //free(error_message);
      return 0;
    }
  }
  return 1;
}

string *AlphaDate(string *strtime,time_t *start_time){
  struct tm *info;
  info = localtime(start_time);

  strfmt(strtime,"%02d-%02d-%02d.%02d%02d%02d",1900 + info->tm_year,1 + info->tm_mon,info->tm_mday,info->tm_hour,info->tm_min,info->tm_sec);
  return strtime;
}

time_t StartScript(char *program,char *cmd){
  time_t start_time;

  FILE *stime_log_file_fh = NULL;

  char *remote_addr = NULL;

  //Initialize RSAT, e.g. environmental variables
  start_time = InitRSAT(program,cmd);

  //If log-file specified,report stats to it
  if(getenv("start_time") && (strcmp(getenv("start_time"),"1") ==0 )){
    //Program name
    if (!program) program = "undefined";

    //Get remote address i.e. Web queries
    remote_addr = getenv("REMOTE_ADDR");
    if(!remote_addr) remote_addr = "";

    //Write header of the exect time log file if file does not exist
    if (access(start_time_log_file->buffer,F_OK) == -1) {
      if(verbose >= 4) RsatInfo("Creating start script log file",start_time_log_file->buffer,NULL);
      stime_log_file_fh = OpenOutputFile(stime_log_file_fh,start_time_log_file->buffer);
      fprintf(stime_log_file_fh,"#start_date.time\thostname\tPID\tusername\tscript_name\tcommand\tremote_addr\n");
      fclose(stime_log_file_fh);
      chmod(start_time_log_file->buffer,0666);
    }

    //Update process information to log file
    if(verbose >= 4) RsatInfo("Updating start script log file",start_time_log_file->buffer,NULL);
    stime_log_file_fh = OpenAppendFile(stime_log_file_fh,start_time_log_file->buffer);
    fprintf(stime_log_file_fh, "%s\t%s\t%d\t%s\t%s\t%s\t%s\n",date->buffer,HOSTNAME->buffer,*PID,LOGIN->buffer,PROGRAM->buffer,CMD->buffer,remote_addr);
    fclose(stime_log_file_fh);
    chmod(start_time_log_file->buffer,0666);

  }

  return start_time;
}

void ReportExecutionTime(time_t start_time){
  //Declare variables
  time_t end_time;
  struct tms process_time;

  string *end_tstr = NULL;

  //Allocate memory for variables
  end_tstr = strnewToList(&RsatMemTracker);

  if( time(&end_time) == -1 ) RsatFatalError("Time is not available in ReportExecutionTime()", NULL);

  //AlphaDate(start_tstr,&start_time);
  AlphaDate(end_tstr,&end_time);

  if(times(&process_time) == -1){
    RsatWarning("Unable to retrieve time information",NULL);
    return;
  }

  //Report execution time
  printf(
  "; Host name\t%s\n"
  "; Job started\t%s\n"
  "; Job done\t%s\n"
  "; Seconds\t%.2f\n"
  ";\tuser\t%.2f\n"
  ";\tsystem\t%.2f\n"
  ";\tcuser\t%.2f\n"
  ";\tcsystem\t%.2f\n",
  HOSTNAME->buffer,
  date->buffer,
  end_tstr->buffer,
  difftime(end_time,start_time),
  (float)process_time.tms_utime/(float)sysconf (_SC_CLK_TCK),
  (float)process_time.tms_stime/(float)sysconf (_SC_CLK_TCK),
  (float)process_time.tms_cutime/(float)sysconf (_SC_CLK_TCK),
  (float)process_time.tms_cstime/(float)sysconf (_SC_CLK_TCK));

  //Update server script
  UpdateExecTimeLogFile(start_time,end_time);

  //Remove temporal strings
  RsatMemTracker = relem( (void*)end_tstr,RsatMemTracker );

  return;
}
/*
  NOTE.WSG (2017-14-06). I should add the uncompress and compress
  forms for I/O operations.
*/
FILE *OpenInputFile(FILE *filehandle,char *filename){
  if(verbose >= 10) RsatInfo("Opening input stream",filename,NULL);
  if (access(filename,F_OK) != -1) {
    if((filehandle = fopen(filename,"r")) == NULL){
      RsatFatalError("Fail to open file",filename,"in OpenInputFile()",NULL);
    }
  }else{
    RsatFatalError("Fail to open file,",filename,"does not exist",NULL);
  }
  return filehandle;
}

FILE *OpenOutputFile(FILE *filehandle,char *filename){
  if(verbose >= 10) RsatInfo("Opening output stream",filename,NULL);
    if((filehandle = fopen(filename,"w")) == NULL){
      RsatFatalError("Fail to open file",filename,"in OpenOutputFile()",NULL);
    }

  return filehandle;
}

FILE *OpenAppendFile(FILE *filehandle,char *filename){
  if(verbose >= 10) RsatInfo("Opening output stream",filename,NULL);
    if((filehandle = fopen(filename,"a")) == NULL){
      RsatFatalError("Fail to open file",filename,"in OpenAppendFile()",NULL);
    }

  return filehandle;
}

//TODO WSG. Update from here and below
/*Test if the passed directory PATH already exists, if no, create it with the appropiate
  masks and permissions. A string containing the directory PATH, a Umask and
  a Chmod code is needed. The return value is 0 on failure or 1 on success. */
int CheckOutDir(string *output_dir,mode_t Umask, mode_t Chmod){
  //Declare variables
  struct stat output_exists;
  string *cmd = NULL;
  //char cmd[BASE_STR_LEN];

  //Allocate memory for variables
  cmd = strnewToList(&RsatMemTracker);

  //Validate input
  if(!output_dir) return 0;

  //Check if directory already exists
  if ( stat(output_dir->buffer,  &output_exists) == 0 && S_ISDIR(output_exists.st_mode) ){
    if(verbose >= 12) RsatWarning("Directory",output_dir->buffer,"already exists",NULL);
    RsatMemTracker = relem( (void*)cmd, RsatMemTracker );
    return 1;
  }

  //Set default values for umask and chmod if 0's have been passed as parameter
  if(!Umask) Umask = 0002;
  if(!Chmod) Chmod = 0777;

  //Specify a mask and create new directory
  umask(Umask);
  if ( mkdir(output_dir->buffer,Chmod) == -1 ){
    strfmt(cmd,"mkdir -p %s",output_dir->buffer);
    //sprintf(cmd,"mkdir -p %s",output_dir);
    if (system(cmd->buffer) == -1) return 0;
  }

  //Specify permissions for directory
  chmod(output_dir->buffer,Chmod);

  //Remove tmp strings
  RsatMemTracker = relem( (void*)cmd, RsatMemTracker );

  //Test if directory has been created satisfactorily
  if( !(stat(output_dir->buffer,  &output_exists) == 0 && S_ISDIR(output_exists.st_mode)) ) return 0;
  if(verbose >= 12) RsatInfo("Directory",output_dir->buffer,"was created successfully!",NULL);

  return 1;
}

//NOTE. Missing substitution from /+ to /
/*Splits a given filename in PATH and filename. An array of 2 elements
  is needed to save each part and also the filename is required.*/
string *SplitFileName(string *token[],char *filename){
  //Declare variables
  char *split_index    = NULL;
  string *tmp_filename = NULL;

  //Validate input
  if(filename == NULL) return NULL;

  //Allocate memory for variables
  tmp_filename = strnewToList(&RsatMemTracker);

  //Initialize tokens with empty string
  strcopy(token[0],"");
  strcopy(token[1],"");

  //Copy filename to a string
  strcopy(tmp_filename,filename);

  //Find position of '/' char at filename
  split_index = strrchr(tmp_filename->buffer,'/');

  //If index has been found
  if(split_index != NULL){
    //Insert a '\0' at split index in order to split the
    //content of the tmp filename
    *split_index = '\0';
    //Copy both parts to separate strings
    strcopy(token[0],tmp_filename->buffer);
    strcopy(token[1],split_index + 1);
  }else{
    strcopy(token[1],tmp_filename->buffer);
  }

  //Remove tmp strings
  RsatMemTracker = relem( (void*)tmp_filename,RsatMemTracker );

  return token[0];
}

/*Gets absolute PATH to $RSAT/public_html/tmp. A string is needed to write down
  the resulting PATH. It returns a pointer to the string object*/
string *Get_pub_temp(string *public_temp_dir){
  strfmt(public_temp_dir,"%s/public_html/tmp",RSAT->buffer);
  return public_temp_dir;
}

/*Gets absolute path to RSAT tmp directory. A string is needed to write down
  the resulting PATH. It returns a pointer to the string object. If fails
  a FatalError is raised.*/
string *Get_temp_dir(string *tmp_dir){
  char user[]    = "temp_user";
  char *username = NULL;
  char *home     = NULL;
  struct passwd *pwd;
  struct tm *info;
  time_t myTime;

  //Get time_t
  if( time(&myTime) == -1 ) RsatFatalError("Time is not available in Get_temp_dir()", NULL);
  info = localtime(&myTime);

  //Get real userID
  pwd = getpwuid(getuid());
  username = (pwd == NULL) ? user : pwd->pw_name;

  //Build the tmp_dir PATH
  if ( getenv("RSA_OUTPUT_CONTEXT")  ){
    if( strcmp(getenv("RSA_OUTPUT_CONTEXT"),"cgi") == 0 || strcmp(getenv("RSA_OUTPUT_CONTEXT"),"RSATWS") == 0){
      Get_pub_temp(tmp_dir);
      strccat(tmp_dir,"/%s",username);
    }
  } else {
    home = getenv("HOME");
    if(!home) RsatFatalError("Unable to find environmental variable $HOME", NULL);
    strfmt(tmp_dir,"%s/.rsat_tmp_dir",home);
  }
  //Build directory name using current date and append it to PATH
  strccat(tmp_dir,"/%04d/%02d/%02d",1900 + info->tm_year,1 + info->tm_mon,info->tm_mday);

  if (verbose >= 12) RsatInfo("Get_temp_dir() result", tmp_dir->buffer,NULL);

  return tmp_dir;
}

/*    tmp_file: string to write down the resulting tmp file name.
      tmp_dir: if NULL, the default RSAT temporal dir is used.
      tmp_prefix: prefix for the file name.
      add_date (value 0 or 1): if 1, the date is added to the suffix.
      make_dir (value 0 or 1): if 1, create a temporary directory rather than temporary file.
      protect  (value 0 or 1): if 1, create a index.html for cgi purpose.*/
string *make_temp_file(string *tmp_file,char *tmp_dir,char *tmp_prefix,int add_date,int make_dir, int protect){
  //Declare variables
  FILE *fh_popen = NULL;

  string *dir_and_file[2] = {NULL};

  string *real_tmp_dir    = NULL;
  string *cmd_mktmp       = NULL;

  //char *real_tmp_dir    = NULL;
  //char *dir_and_file[2] = {NULL};
  //char cmd_mktmp[BASE_STR_LEN];
  //if( (real_tmp_dir    = (char *)malloc(sizeof(char) * BASE_STR_LEN)) == NULL) RsatFatalError("Unable to allocate memory in make_temp_file()",NULL);
  //if( (dir_and_file[0] = (char *)malloc(sizeof(char) * BASE_STR_LEN)) == NULL) RsatFatalError("Unable to allocate memory in make_temp_file()",NULL);
  //if( (dir_and_file[1] = (char *)malloc(sizeof(char) * BASE_STR_LEN)) == NULL) RsatFatalError("Unable to allocate memory in make_temp_file()",NULL);

  //Allocate memory for variables
  //dir_and_file = getokens(2);

  real_tmp_dir    = strnewToList(&RsatMemTracker);
  cmd_mktmp       = strnewToList(&RsatMemTracker);
  dir_and_file[0] = strnewToList(&RsatMemTracker);
  dir_and_file[1] = strnewToList(&RsatMemTracker);

  //If tmp_dir has been passed without NULL,copy it to real_tmp_dir
  //else initialize real_tmp_dir with NULL, IN ORDER to use the
  //default RSAT temporal directory
  if(tmp_dir) {
    strcopy(real_tmp_dir,tmp_dir);
  } else {
    strcopy(real_tmp_dir,"");
    //real_tmp_dir[0] = '\0';
  }

  //If tmp_file has been passed without argument,i.e. no string to write
  //result,delete tmp allocated variables,then exit function and return NULL
  if (!tmp_file) {
    RsatMemTracker = relem( (void*)dir_and_file[1], RsatMemTracker );
    RsatMemTracker = relem( (void*)dir_and_file[0], RsatMemTracker );
    RsatMemTracker = relem( (void*)cmd_mktmp, RsatMemTracker );
    RsatMemTracker = relem( (void*)real_tmp_dir, RsatMemTracker );

    //free(dir_and_file[0]);
    //free(dir_and_file[1]);
    //free(real_tmp_dir);

    return NULL;
  }

  //If a prefix has been passed make sure there is no path by splitting filename
  if (tmp_prefix) {
    //Split tmp_prefix in PATH and FILENAME
    SplitFileName(dir_and_file,tmp_prefix);

    //Test if a dir PATH is already written at real_tmp_dir(passed to function) or has been
    //found by splitting tmp_prefix
    if (strcmp(real_tmp_dir->buffer,"") != 0  && strcmp(dir_and_file[0]->buffer,"") != 0){
      //Append found PATH dir from splitted FILENAME to function's passed directory
      strccat(real_tmp_dir,"/%s",dir_and_file[0]->buffer);
    //If ONLY a dir PATH from splitted FILENAME is available,copy it to real_tmp_dir
    } else if(strcmp(dir_and_file[0]->buffer,"") != 0){
      strcopy(real_tmp_dir,dir_and_file[0]->buffer);
    }
  //If no prefix has been passed,write down a default 'tmp' name to resulting string
  } else {
    strcopy(dir_and_file[1],"tmp");
  }

  //If no dir PATH is contained yet on real_tmp_dir,because no dir PATH was passed and
  //no dir PATH was found on splitted tmp_prefix, retrieve default RSAT tmp dir PATH and
  //copy it to real_tmp_dir
  if (strcmp(real_tmp_dir->buffer,"") == 0){
    Get_temp_dir(real_tmp_dir);
  }

  //Check if directory exists or has been created successfully
  if ( CheckOutDir(real_tmp_dir,0,0755) == 0 ) RsatFatalError("Unable to create",real_tmp_dir->buffer,"in make_temp_file()",NULL);

  //If this parameter has been passed as 1 to function, check if cgi environment is found and
  //create and index.html with the legend 'Access forbidden'
  if(protect){
    if(getenv("RSA_OUTPUT_CONTEXT") && strcmp(getenv("RSA_OUTPUT_CONTEXT"),"cgi") == 0){
      //Declare variables
      string *index_file = NULL;
      FILE *fh_index = NULL;

      //Allocate memory for variables
      index_file = strnewToList(&RsatMemTracker);

      //Create index.html file on tmp directory
      strfmt(index_file,"%s/index.html",tmp_dir);
      //sprintf(index_file,"%s/index.html",tmp_dir);

      //Check if index.html exists and write content to it if it does not exists
      if(access(index_file->buffer,F_OK) != -1){
        fh_index = OpenOutputFile(fh_index,index_file->buffer);
        fprintf(fh_index, "<html>\n<b>Access forbidden</b>\n</html>\n");
        fclose(fh_index);
      }
      //Remove tmp strings
      RsatMemTracker = relem( (void*)index_file,RsatMemTracker );
    }
  }

  //If this parameter has been passed as 1 to function, add date and time to FILENAME
  if(add_date){
    //Declare variables
    struct tm *info;
    time_t myTime;

    //Get date/time and append it to FILENAME
    if( time(&myTime) == -1 ) RsatFatalError("Time is not available in Get_temp_dir()", NULL);
    info = localtime(&myTime);
    strccat(dir_and_file[1],"_%02d-%02d-%02d.%02d%02d%02d",1900 + info->tm_year,1 + info->tm_mon,info->tm_mday,info->tm_hour,info->tm_min,info->tm_sec);
    //sprintf(dir_and_file[1],"%s_%02d-%02d-%02d.%02d%02d%02d",dir_and_file[1],1900 + info->tm_year,1 + info->tm_mon,info->tm_mday,info->tm_hour,info->tm_min,info->tm_sec);
  }
  //Check whether a directoy or a file is going to be created
  if(make_dir){
    strfmt(cmd_mktmp,"mktemp -u %s/%s_XXXXXX",real_tmp_dir->buffer,dir_and_file[1]->buffer);
  } else {
    strfmt(cmd_mktmp,"mktemp -u -d %s/%s_XXXXXX",real_tmp_dir->buffer,dir_and_file[1]->buffer);
  }


  //Execute cmd in cmdline and write obtained FILENAME from mktemp to tmp_file
  if ( (fh_popen = popen(cmd_mktmp->buffer,"r")) == NULL ) RsatFatalError("Unable to popen in make_temp_file()",NULL);
  tmp_file->size = 0;
  while ( (tmp_file->buffer[tmp_file->size] = fgetc(fh_popen)) != '\n' ) {
    //printf("This is the out of popen : %s \n",tmp_file->buffer );
    tmp_file->size++;
    strlimt(tmp_file);
  }
  //Add '\0'-end to string and +1 to size
  tmp_file->buffer[tmp_file->size] = '\0';
  tmp_file->size++;
  //Close filehandler
  pclose(fh_popen);

  //Remove tmp variables
  RsatMemTracker = relem( (void*)dir_and_file[1],RsatMemTracker );
  RsatMemTracker = relem( (void*)dir_and_file[0],RsatMemTracker );
  RsatMemTracker = relem( (void*)cmd_mktmp,RsatMemTracker );
  RsatMemTracker = relem( (void*)real_tmp_dir,RsatMemTracker );

  if (verbose >= 12) RsatInfo("make_temp_file() result", tmp_file->buffer,NULL);
  return tmp_file;
}

//NOTE. All tokens must be NULL when they havent filled yet.
/*Test if a line has reached its max length. If true,
  it reallocs new memory by an adding factor of +BASE_STR_LEN/2,
  and also remaps the token pos pointers to new line buffer*/
string *limlinetok(string *line, char **token,int numtok){
  int *tokenpos = NULL;
  int i = 0;

  //Test if current line buffer is at max length
  if (line->size == line->length - 2 ) {
    //NOTE WSG. Remove printf("This is inside-inside line->size == %d and line->length == %d\n", (int)line->size, (int)line->length );
    //Find memory to allocate a temporal array for distance between
    //tokens and buffer memory pointer.NOTE WSG(2017-20-10) Do not
    //forget to add a stdnewToList() for a more cleaner way.
    tokenpos = (int*)_MemTrackMalloc(sizeof(int) * numtok, &RsatMemTracker, "limlinetok");
    //Fill temporal array with distance between tokens
    // and buffer memory pointer
    for (i = 0; i < numtok; i++) {
      //If NULL is found in tokens to remap, break. This
      //means there are no further tokens in array to remap
      if (token[i] == NULL) break;
      tokenpos[i] = (int)( (line->buffer - token[i]) * -1 ); //Times *-1 to get +number
    }

    //Realloc buffer with more space
    line->buffer = stralloc(line, (int)(line->length + BASE_STR_LEN/2));
    //Remap tokens to their new buffer location
    for (i = 0; i < numtok; i++) {
      if (token[i] == NULL) break;
      token[i] = line->buffer + tokenpos[i];//Perhaps a casting is needed.
    }

    //Free memory for temporal array
    RsatMemTracker = relem( (void*)tokenpos, RsatMemTracker );

  }

  return line;
}

/*Retrieves the absolute path of the contigs.txt file on a given valid
  genome directory. A string is needed to write down the path as 1st
  argument. As second argument it takes a valid genome directory and
  appends the /contigs.txt to the genome directory. It returns a pointer
  to the given string. This function always succeeds.
  The current tab file format is:
  id   accession       version type    length  description
  11      chromosome:GRCh37:11:1:135006516:1      GRCh37  chromosome      135006516       chromosome 11*/
string *Get_contigs_file(string *contigs_file,string *genome_dir){
  strfmt(contigs_file,"%s/contigs.txt",genome_dir->buffer);
  if(verbose >= 5) RsatInfo("Get_contigs_file() result", contigs_file->buffer,NULL);

  return contigs_file;
}

/*Retrieves the absolute path of the contig.tab file on a given valid
  genome directory. A string is needed to write down the path as 1st
  argument. As second argument it takes a valid genome directory and
  appends the /contig.tab to the genome directory. It returns a pointer
  to the given string. This function always succeeds.
  The current txt file format is:
  file  accession
  chromosome_GRCh37_11_1_135006516_1.raw  chromosome:GRCh37:11:1:135006516:1*/
string *Get_contig_file(string *contig_file,string *genome_dir){
  strfmt(contig_file,"%s/contig.tab",genome_dir->buffer);
  if(verbose >= 5) RsatInfo("Get_contig_file() result", contig_file->buffer,NULL);

  return contig_file;
}

/*It builds a TRIE search struct, with chromosomes as keys and filenames
  as values. It takes as parameter a valid genome directory installed on
  $RSAT/data/genomes. The trie is built by parsing contigs.txt, which
  contains (id->accession) and contig.tab, which contains (rawFile<-accession).
  Both files are contained inside the passed directory. Finally, it returns
  a pointer to the current trie. If an error occurs a Fatal Error is raised.*/
TRIE *Get_file_seq_name(string *genome_dir){
  //Declare variables
  FILE *fh_contig  = NULL;
  FILE *fh_files   = NULL;

  char **token   = NULL;
  //char *token[3] = {NULL};

  string *line         = NULL;
  string *contig_file  = NULL;
  string *contigs_file = NULL;
  string *chromos_file = NULL;

  //char chromos_file[BASE_STR_LEN];
  //int i = 0;
  int j = 0;

  TRIE *trieContig = NULL;
  TRIE *trieChromo = NULL;

  //Allocate memory for variables
  token = getokens(3);

  line         = strnewToList(&RsatMemTracker);
  contig_file  = strnewToList(&RsatMemTracker);
  contigs_file = strnewToList(&RsatMemTracker);
  chromos_file = strnewToList(&RsatMemTracker);

  //Get absolute path for contig.tab and contigs.txt
  Get_contig_file(contig_file,genome_dir);
  Get_contigs_file(contigs_file,genome_dir);

  //Open files
  fh_contig  = OpenInputFile(fh_contig,contig_file->buffer);
  fh_files   = OpenInputFile(fh_files,contigs_file->buffer);

  //Start both search-tries,the 1st maps (id->accession)
  //the 2nd maps (id->rawFile).
  trieContig = TrieStart();
  trieChromo = TrieStart();

  //1st token is start of line.
  token[0] = line->buffer;

  //Parse contig.tab file (id->accession)
  while ( fread( (line->buffer + line->size), 1, 1, fh_contig ) == 1 ) {
    //Get tokens by substituting '\t' for '\0' to obtain key-value
    //pair for trie, i.e. token[0] and token[1] respectively.
    if (line->buffer[line->size] == '\t' && j < 2) {
      token[++j] = line->buffer + line->size + 1;
         line->buffer[line->size] = '\0';
    }

    //Insert Key-Value pair at the end of line
    if (line->buffer[line->size] == '\n') {
      //Reinitialize counters
      j = 0;
      line->size = 0;
      //If a comment is found skip line
      if (token[0][0] == '-') continue;
      //Insert key-value to trie
      TrieInsert(trieContig,token[1],token[0]);
      continue;
    }
    //Resize line string if limit has reached
    limlinetok(line,token,3);

    //Add +1 to current line size
    line->size++;
  }

  //Parse contigs.txt file (id->accession)
  while ( fread( (line->buffer + line->size), 1, 1, fh_files ) == 1 ) {
    //Get tokens by substituting '\t' for '\0' to obtain key-value
    //pair for trie, i.e. token[0] and token[1] respectively.
    if (line->buffer[line->size] ==  '\t' && j < 2) {
      token[++j] = line->buffer + line->size + 1;
         line->buffer[line->size] = '\0';
    }

    //Insert Key-Value pair at the end of line
    if (line->buffer[line->size] == '\n') {
      line->buffer[line->size] = '\0';
      //Reinitialize counters
      line->size = 0;
      j = 0;
      //Absolute path of the file
      strfmt(chromos_file,"%s/%s",genome_dir->buffer,token[0]);
      //sprintf(chromos_file,"%s/%s",genome_dir,token[0]);
      //Test for file access, if one fails raise a FatalError
      if(access(chromos_file->buffer,F_OK) == -1) {
        TrieEnd(trieContig);
        TrieEnd(trieChromo);
        fclose(fh_contig);
        fclose(fh_files);
        RsatFatalError("File",chromos_file->buffer,"does not exist",NULL);
      }
      //Insert key-value to trie
      TrieInsert(trieChromo,TrieSearch(trieContig,token[1]),token[0]);
      continue;
    }
    //Resize line string if limit has reached
    limlinetok(line,token,3);

    //Add +1 to current line size
    line->size++;
  }

  //Remove trie
  //TrieEnd(trieContig);

  //Close file handlers
  fclose(fh_contig);
  fclose(fh_files);

  //Remove allocate variables
  RsatMemTracker = relem( (void*)trieContig  ,RsatMemTracker );
  RsatMemTracker = relem( (void*)line        ,RsatMemTracker );
  RsatMemTracker = relem( (void*)contig_file ,RsatMemTracker );
  RsatMemTracker = relem( (void*)contigs_file,RsatMemTracker );
  RsatMemTracker = relem( (void*)chromos_file,RsatMemTracker );
  RsatMemTracker = relem( (void*)token       ,RsatMemTracker );

  return trieChromo;
}

/*Parses the $RSAT/public_html/data/supported_organisms_ensembl.tab file,which contains
  the installed organism from ensembl at the current computer/server. A string is needed
  to write down the path. Finally, it returns the pointer to the string,in case of failure
  a Fatal Error is raised. The current tab file format is:
    #id	species	assembly_version	db	ensembl_version	update_date	species_directory
  Where desired tokens are:
    Species: token[1],assembly: token[2],Release: token[4],Directory: token[6]
  NOTE WSG(2017-10-10).If this file format changes I will not be able to parse
  it anymore,PLEASE review carefully with JvH and AMR */
string *Get_species_dir_from_supported_file(string *species_dir,char *species,char *assembly,char *release,char *supported_file){
  string *line = NULL;

  char *rsat_idx          = NULL;
  char **token            = NULL;
  //char *token[7]          = {NULL};

  int   j                 = 0;

  FILE *fh_supportedFile  = NULL;

  //Allocate memory for array of tokens
  token = getokens(7);

  //Allocate memory for line string
  line = strnewToList(&RsatMemTracker);

  fh_supportedFile = OpenInputFile(fh_supportedFile,supported_file);
  token[j] = line->buffer;

  //TODO WSG. Remove
  //printf("\n ---INSIDE Get_species_dir_from_supported_file\n");
  //printMemstr(RsatMemTracker);

  //!feof(fh_supportedFile) for this while?
  while ( fread((line->buffer + line->size),1,1,fh_supportedFile) == 1) {
    //For each tab found write a '\0' in order to create token
    if (line->buffer[line->size] == '\t') {
      token[++j] = line->buffer + line->size + 1;
      line->buffer[line->size] = '\0';
    }

    //When end-of-line is found,process the tokens and also write a '\0' instead of '\n'
    if (line->buffer[line->size] == '\n') {
      line->buffer[line->size] = '\0';
      line->size = 0;
      j = 0;
      token[1][0] = toupper(token[1][0]);

      //Check if release and assembly where passed as query options in order to compare properly
      if (release && assembly) {
        if (strcmp(token[1],species) == 0 && strcmp(token[2],assembly) == 0 && strcmp(token[4],release) == 0) {
          //Test if species_directory field is empty
          if (token[6] != '\0'){
            if( (rsat_idx = strchr(token[6],'}')) != NULL ){
              strfmt( species_dir, "%s%s", RSAT->buffer, (rsat_idx + 1) );
            } else {
              RsatFatalError("Get_species_dir_from_supported_file() could not identify species",
                              species,"from release",release,
                              "in the organism table,species_dir field was empty\n",NULL);
            }
          } else {
            RsatFatalError("Get_species_dir_from_supported_file() could not identify species",
                            species,"from release",release,
                            "in the organism table,species_dir field was empty\n",NULL);
          }
          RsatMemTracker = relem( (void*)token,RsatMemTracker );
          RsatMemTracker = relem( (void*)line ,RsatMemTracker );
          fclose(fh_supportedFile);
          return species_dir;
        }
      //Check if ONLY release was passed as query option in order to compare properly
      } else if (release) {

        if (strcmp(token[1],species) == 0 && strcmp(token[4],release) == 0) {
          //Test if species_directory field is empty
          if (token[6] != '\0') {
            if( (rsat_idx = strchr(token[6],'}')) != NULL ) {
              strfmt( species_dir, "%s%s", RSAT->buffer, (rsat_idx + 1) );
            } else {
              RsatFatalError("Get_species_dir_from_supported_file() could not identify species",
                              species,"from release",release,
                              "in the organism table,species_dir field was not correct\n",NULL);
            }
          } else {
            RsatFatalError("Get_species_dir_from_supported_file() could not identify species",
                            species,"from release",release,
                            "in the organism table,species_dir field was empty\n",NULL);
          }
          RsatMemTracker = relem( (void*)token,RsatMemTracker );
          RsatMemTracker = relem( (void*)line ,RsatMemTracker );
          fclose(fh_supportedFile);
          return species_dir;
        }
      }
      initokadd(line,token,7);
      continue;

    }
    //Resize line string if limit has reached
    limlinetok(line,token,7);
    //Add +1 to current line size
    line->size++;

  }
  RsatMemTracker = relem( (void*)token,RsatMemTracker );
  RsatMemTracker = relem( (void*)line ,RsatMemTracker );
  fclose(fh_supportedFile);
  RsatFatalError("Get_species_dir_from_supported_file() could not identify species",
                  species,"from release",release,
                  "in the organism table.\n\tPlease refer to install-ensembl-genome\n",NULL);
  return NULL; //NOTE:NULL or species_dir?
}
/*Retrieves the $RSAT/public_html/data absolute path. A string is needed to write
  down the path. Finally, it returns the pointer to the string.*/
string *Get_data_dir(string *data_dir){
  strfmt(data_dir,"%s/public_html/data",RSAT->buffer);
  //strcpy(data_dir,getenv("RSAT"));
  //strcat(data_dir,"/public_html/data");
  if(verbose >= 5) RsatInfo("Get_data_dir() result", data_dir->buffer,NULL);
  return data_dir;
}

/*Retrieves the supported_organisms_ensembl.tab absolute file path,which contains
  the installed organism from ensembl at the current computer/server. A string is
  needed to write down the path. Finally, it returns the pointer to the string */
string *Get_supported_file(string *supported_file){
  //It doesn't support a predefined PATH to attach names
  //TODO WSG. Remove
  //printf("\n ---INSIDE Get_supported_file\n");
  //printMemstr(RsatMemTracker);

  supported_file = Get_data_dir(supported_file);

  //TODO WSG. Remove
  //printf("\n---1 INSIDE Get_supported_file_after_Get_data_dir\n");
  //printMemstr(RsatMemTracker);

  strccat(supported_file,"/supported_organisms_ensembl.tab");
  if(verbose >= 5) RsatInfo("Get_supported_file() result", supported_file->buffer,NULL);

  //strcat(supported_file,"/supported_organisms_ensembl.tab");
  //TODO WSG. Remove
  //printf("\n---2 INSIDE Get_supported_file_after_strccat\n");
  //printMemstr(RsatMemTracker);

  return supported_file;
}

/*Returns the absolute path for the genomes directories at $RSAT/public_html/data*/
string *Get_genomes_dir(string *genomes_dir){
  Get_data_dir(genomes_dir);
  strccat(genomes_dir,"/genomes");
  if(verbose >= 5) RsatInfo("Get_genomes_dir() result", genomes_dir->buffer,NULL);

  return genomes_dir;
}

/*Retrieves the species directory path for the queried species,assembly, release and/or species suffix.
  A string is needed to write down the path. Finally, it returns the pointer to the string. If one
  of the dependant function fails it will raise a Fatal Error. */
string *Get_species_dir(string *species_dir,char *species, char *assembly, char *release, char *species_suffix){
  string *supported_file = NULL;
  supported_file = strnewToList(&RsatMemTracker);

  //Will modify the species string, i.e. uppercase of first letter.
  //char supported_file[BASE_STR_LEN];
  species[0] = toupper(species[0]);
  //TODO WSG. Remove
  //printf("\n---INSIDE Get_species_dir\n");
  //printMemstr(RsatMemTracker);

  supported_file = Get_supported_file(supported_file);

  //Parse and query the $RSAT/public_html/data/supported_organisms_ensembl.tab file
  species_dir = Get_species_dir_from_supported_file(species_dir,species,assembly,release,supported_file->buffer);
  if(verbose >= 5) RsatInfo("Get_species_dir() result", species_dir->buffer,NULL);

  //Remove tmp allocated variable
  RsatMemTracker = relem( (void*)supported_file,RsatMemTracker );

  return species_dir;
}

/*Retrieves genome directory path for the queried species, assembly, release and/or species suffix.
  A string is needed to write down the path. Finally, it returns the pointer to the string*/
string *Get_genome_dir(string *genome_dir,char *species, char *assembly, char *release, char *species_suffix){
  //TODO WSG. Remove
  //printf("\n ---INSIDE  Get_species_dir\n");
  //printMemstr(RsatMemTracker);

  genome_dir = Get_species_dir(genome_dir,species,assembly,release,species_suffix);
  strccat(genome_dir,"/genome");
  if(verbose >= 5) RsatInfo("Get_genome_dir() result",genome_dir->buffer,NULL);
  return genome_dir;
}

/*Retrieves variations directory path for the queried species, assembly, release and/or species suffix.
  A string is needed to write down the path. Finally, it returns the pointer to the string*/
string *Get_variation_dir(string *var_dir,char *species, char *assembly, char *release, char *species_suffix){
  var_dir = Get_species_dir(var_dir,species,assembly,release,species_suffix);
  strccat(var_dir,"/variations");
  if(verbose >= 5) RsatInfo("Get_variation_dir() result",var_dir->buffer,NULL);

  return var_dir;
}

//NOTE. WSG(2017-21-04). I use this for .raw files to load it in a variable
/*Loads an entire content of a file into a memory variable, a char* for the name of the file
  is needed. If function succeeds a pointer to the variable holding the file content is
  returned, on failure a FatalError is raised.*/
char *Get_sequence(char *sequence_file){
  FILE *fh_sequence = NULL;
  char *sequence    = NULL;
  long long read_bytes;
  struct stat file;

  stat(sequence_file,&file);
  fh_sequence = OpenInputFile(fh_sequence,sequence_file);

  if((sequence = (char *)malloc(sizeof(char) * (file.st_size) )) == NULL) RsatFatalError("Unable to allocate memory for Get_sequence()",NULL);
  if ( (read_bytes = fread(sequence,sizeof(char),file.st_size,fh_sequence)) !=  file.st_size ) RsatFatalError("Unable to read from stream %lld objects in Get_sequence(), read %lld",file.st_size,read_bytes,NULL);

  fclose(fh_sequence);
  //printf("Success for file : %s\n", sequence_file);
  return sequence;
}

int switch_strand(char *sequence){
  int i;
  for (i=0;sequence[i] == '\0';i++){
    switch (sequence[i]){
      case 'A':
          sequence[i] = 'T';
          break;

      case 'C':
          sequence[i] = 'G';
          break;

      case 'G':
          sequence[i] = 'C';
          break;

      case 'T':
          sequence[i] = 'A';
          break;

      case 'N':
          break;
      case ',':
          break;
      default:
          RsatWarning("This allele sequence is not supported");
          return 0;
    }
  }
  return 1;
}
