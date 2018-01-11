#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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
  int          offset;
  site              D;
  site              R;
  struct _scan  *next;
} scan;

typedef struct _varscan {
  variant     *variation;
  scan        *scan_info;
  struct _varscan  *next;
} varscan;

string *GetProgramPath(string *program_path, char *program_name, int die_on_error, stringlist *preferred_path);
void argToList(stringlist *list, char **value);
void strlistfree(stringlist *delete);
void strlistend(stringlist *group);
stringlist *strlistnew(void);
stringlist *strlistnewToList(memstd **List);
stringlist *strlistadd(stringlist *group);
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
int CheckOutDir(string *output_dir,int Umask, int Chmod);
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

void printHeader(int PhasedFile, FILE *fh_outputFile);
void processHaplotypes(int mml, variant *firstVar, variant *HaploGroup, variant *lastVar, variant *printVar, char *sequence, FILE *fout,
                       string *varCoords, string *IDs, string *SOs, string *alleleFreqs, string *Haplotype1, string *Haplotype2,
                       string *Haplotype1Sequence, string *Haplotype2Sequence );
int CheckOutOfIndex( char *chr,char *start, char *end, int mml, long long maxsize );

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
  "    retrieve-variation-seq\n"
  "\n"
  "VERSION\n"
  "    2.0\n"
  "\n"
  "DESCRIPTION\n"
  "    Given a set of IDs for polymorphic variations, retrieve the\n"
  "    corresponding variants and their flanking sequences, in order to scan\n"
  "    them wiht the tool variation-scan.\n"
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
  "     retrieve-variation-seq -species species_name (-release # | -assembly assembly)  \\\n"
  "       [-i #inputfile] [-format variation_format] \\\n"
  "       [-col ID_column] [-mml #] [-o outputfile] [-v #] [...]\n"
  "\n"
    "Example\n"
    "    Get variation sequence of Homo_sapiens from a bed file\n"
    "\n"
    "      retrieve-variation-seq -v 2 \\\n"
    "        -species Homo_sapiens -release 84 -format bed -assembly GCRh37 \\\n"
    "        -i $RSAT/public_html/demo_files/sample_regions_for_variations_hg19.bed \\\n"
    "        -mml 30 \\\n"
    "        -o variations.varSeq\n"
    "\n"
  "INPUT FORMAT\n"
    "Genomic coordinate file\n"
    "  The option -i allows to specify a genomic coordinate file in bed format.\n"
    "  The program only takes into account the 3 first columns of the bed file,\n"
    "  which specify the genomic coordinates.\n"
"\n"
    "  Note (from Jacques van Helden): the UCSC genome browser adopts a\n"
    "  somewhat inconsistent convention for start and end coordinates: the\n"
    "  start position is zero-based (first nucleotide of a chromosome/scaffold\n"
    "  has coordinate 0), but the end position is considered not included in\n"
    "  the selection. This is equivalent to have a zero-based coordinate for\n"
    "  the start, and a 1-base coordinate for the end.\n"
"\n"
    "Example of bed file\n"
    "   chr1   3473041 3473370\n"
    "   chr1   4380371 4380650\n"
    "   chr1   4845581 4845781\n"
    "   chr1   4845801 4846260\n"
"\n"
    "  The definition of the BED format is provided on the UCSC Genome Browser\n"
    "  web site (http://genome.ucsc.edu/FAQ/FAQformat#format1).\n"
"\n"
    "  This program only takes into account the 3 first columns, which specify\n"
    "  the genomic coordinates.\n"
"\n"
    "  1. chrom\n"
    "      The name of the chromosome (e.g. chr3, chrY, chr2_random) or\n"
    "      scaffold (e.g. scaffold10671).\n"
"\n"
    "  2. chromStart\n"
    "      The starting position of the feature in the chromosome or scaffold.\n"
    "      For RSAT programs, the first base in a chromosome is numbered 1\n"
    "      (this differs from the UCSC-specific zero-based notation for the\n"
    "      start).\n"
"\n"
    "      Note from Jacques van Helden: the UCSC genome browser adopts a\n"
    "      somewhat inconsistent convention for start and end coordinates: the\n"
    "      start position is zero-based (first nucleotide of a\n"
    "      chromosome/scaffold has coordinate 0), and the end position is\n"
    "      considered not included in the selection. This is equivalent to have\n"
    "      a zero-based coordinate for the start, and a 1-base coordinate for\n"
    "      the end. We find this representation completely counter-intuitive,\n"
    "      and we herefore decided to adopt a 'normal' convention, where:\n"
"\n"
    "      start and end position represent the first and last positions\n"
    "      included in the region of interest.\n"
    "      start and end positions are provided in one-based notation (first\n"
    "      base of a chromosome or contig has coordinate 1).\n"
"\n"
    "  3. chromEnd\n"
    "      The ending position of the feature in the chromosome or scaffold.\n"
"\n"
    "Variation file\n"
    "  See download-ensembl-variation output format.\n"
"\n"
    "Variation ID list\n"
    "  A tab delimited file with id of variation in column.\n"
"\n"
  "OUTPUT FORMAT\n"
"\n"
  "CASE 1 . In case phase information had been found by a supplied varBed\n"
  "         A tab delimited file with the following column content.\n"
"\n"
  "    1. chrom\n"
  "       Chromosome name of Haplotype\n"
"\n"
  "    2. start\n"
  "       Start position of Haplotype\n"
"\n"
  "    3. end\n"
  "       End position of Haplotype\n"
"\n"
  "    4. strand\n"
  "       Strand of Haplotype\n"
"\n"
  "    5. var_coord\n"
  "       Comma-separated list of individual variation coordinates\n"
"\n"
  "    6. id\n"
  "       Comma-separated identifier list of the individual variations\n"
"\n"
  "    7. SO term\n"
  "       Comma-separated sequence ontology (SO) term describing each variation type\n"
"\n"
  "    8. hap1\n"
  "       Comma-separated list from the individual alleles comprising Haplotype 1\n"
"\n"
  "    9. hap\n"
  "       Comma-separated list from the individual alleles comprising the current Haplotype\n"
"\n"
  "    10. Allele \n"
  "       Comma-separated list from the individual allele frequencies\n"
"\n"
  "    11. seq\n"
  "       Sequence of the current Haplotype\n"
"\n"
"\n"
  "CASE 2 . In case phase information had been found by a supplied varBed\n"
  "         A tab delimited file with the following column content.\n"
"\n"
  "    1. chrom\n"
  "        The name of the chromosome (e.g. 1, X, 8...)\n"
"\n"
  "    2. start\n"
  "        The starting position of the variation in the chromosome\n"
"\n"
  "    3. end\n"
  "        The ending position of the variation in the chromosome\n"
"\n"
  "    4. strand\n"
  "        The strand of the variation in the chromosome\n"
"\n"
  "    5. id\n"
  "        ID of the variation\n"
"\n"
  "    8. soterm\n"
  "        SO Term of the the variation(s)\n"
"\n"
  "    7. ref_var\n"
  "        Allele of the variation in the reference sequence\n"
"\n"
  "    8. alt_var\n"
  "        Allele of the variation in the sequence\n"
"\n"
  "    9. allele_freq\n"
  "        Allele frequency\n"
"\n"
  "   10. seq\n"
  "        Sequence of the current variant, flanked by a user-specified neighbouring region\n"
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
  "    -species species_name\n"
  "        Species name. This name must correspond to the species of the\n"
  "        variation/bed/id file if provided.\n"
"\n"
  "    -species_suffix\n"
  "        Species name. This name must correspond to the species of the\n"
  "        variation/bed/id file if provided.\n"
"\n"
  "    -release #\n"
  "        The version of ensembl database (e.g. 84).\n"
"\n"
  "        Note: each Ensembl release contains a specific assembly for each\n"
  "        species. When the option -release is used, the option -assembly\n"
  "        should thus in principle not be used.\n"
"\n"
  "    -assembly #\n"
  "        Assembly (e.g. GRCh37 for the assembly 37 of the Human genome).\n"
"\n"
  "        Note: genome assemblies can cover several successive ensemble\n"
  "        versions. In case of ambiguity, the latest corresponding ensembl\n"
  "        version is used.\n"
"\n"
  "    -i input_file\n"
  "        Input File.\n"
"\n"
  "        The input file specifies a list of query variations. Each row\n"
  "        corresponds to one query.\n"
"\n"
  "        The variations can be provided in various formats (see option\n"
  "        -format below).\n"
"\n"
  "    -format variation_format\n"
  "        Format of the input file\n"
"\n"
  "        Supported formats:\n"
"\n"
  "        varBed\n"
  "            Format of variation files used by all RSAT scripts.\n"
"\n"
  "        id  tab-delimited file with all variation IDs in a given column,\n"
  "            which can be specified by the option -col.\n"
"\n"
  "        bed General format for the description of genomic features (see\n"
  "            https://genome.ucsc.edu/FAQ/FAQformat.html#format1).\n"
"\n"
  "    -mml #\n"
  "        Length of the longest Matrix.\n"
"\n"
  "        The program will adapt the length of the flanks to be extracted on\n"
  "        each side of the variants, in order to be able to align the longest\n"
  "        matrix on both sides.\n"
"\n"
  "        The flanking size on each size will be\n"
"\n"
  "            n = mml -1\n"
"\n"
  "        For example, if the longest matrix of the database contains 30\n"
  "        columns, the flanking size will be 29 base pairs, which is\n"
  "        sufficient to align the SNP at all the positions of the matrix.\n"
"\n"
  "    -col #\n"
  "        Column containing the variation IDs with the input format 'id'.\n"
"\n"
  "        Default : 1\n"
"\n"
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
  //Initialize RsatMemTracker for memory tracking
  //it also add contents to RSAT and Initialize CMD,
  //PROGRAM, and PID.
  InitMemLists();

  //Start Creating CMD string.
  strcopy(CMD,"retrieve-variation-seq");
  strcopy(PROGRAM,"retrieve-variation-seq");

  /////////////////////////////////////////////////
  // Declare variables
  /////////////////////////////////////////////////
  struct stat dir_exists;
  struct stat file;

  FILE *fh_stdin        = NULL;
  TRIE *trieChrom       = NULL;
  FILE *fout            = stdout;
  FILE *fin             = stdin;

  char *input          = NULL;
  char *output         = NULL;
  char *species        = NULL;
  char *species_suffix = NULL;
  char *release    = NULL;
  char *assembly   = NULL;
  char *format     = NULL;
  char *col        = NULL;
  char *seq_search = NULL;
  char *sequence   = NULL;

  char character;

  string *line                = NULL;
  string *genome_dir          = NULL;
  string *variant_dir         = NULL;
  string *outfile_stdin       = NULL;
  string *outfile_var_info    = NULL;
  string *outfile_var_sorted  = NULL;
  string *get_variations_cmd  = NULL;
  string *sort_variations_cmd = NULL;
  string *deletetmps_cmd      = NULL;

  string *curr_chr           = NULL;
  string *seq_file           = NULL;

  long long left_flank       = 0;
  long long right_flank      = 0;
  long long sequence_maxsize = 0;

  int mml       = 29;
  int phased    =  0;
  int inputVars =  0;
  int i;
  int j;
  int k;

  //Default value for column number
  col = (char*)_MemTrackMalloc(sizeof(char) * 2, &RsatMemTracker, "main");
  strcpy(col,"1");

  /////////////////////////////////////////////////
  // Read Arguments
  /////////////////////////////////////////////////
  if (argc <= 1) help();

  for (i = 1; i < argc; i++) {
    if (strcmp(argv[i],"-h") == 0 || strcmp(argv[i],"--help") == 0) help();
    else if (strcmp(argv[i],"-i") == 0 && CheckValOpt(argv+i)) {
      strccat(CMD," %s %s",argv[i],argv[i+1]);
      input = argv[++i];
    } else if (strcmp(argv[i],"-o") == 0 && CheckValOpt(argv+i)) {
      strccat(CMD," %s %s",argv[i],argv[i+1]);
      output = argv[++i];
    } else if (strcmp(argv[i],"-species") == 0 && CheckValOpt(argv+i)) {
      strccat(CMD," %s %s",argv[i],argv[i+1]);
      species = argv[++i];
    } else if (strcmp(argv[i],"-species_suffix") == 0 && CheckValOpt(argv+i)) {
      strccat(CMD," %s %s",argv[i],argv[i+1]);
      species_suffix = argv[++i];
    } else if (strcmp(argv[i],"-release") == 0 && CheckValOpt(argv+i)) {
      strccat(CMD," %s %s",argv[i],argv[i+1]);
      release = argv[++i];
    } else if (strcmp(argv[i],"-assembly") == 0 && CheckValOpt(argv+i)) {
      strccat(CMD," %s %s",argv[i],argv[i+1]);
      assembly = argv[++i];
    } else if (strcmp(argv[i],"-format") == 0 && CheckValOpt(argv+i)) {
      strccat(CMD," %s %s",argv[i],argv[i+1]);
      format = argv[++i];
    } else if (strcmp(argv[i],"-mml") == 0 && CheckValOpt(argv+i)) {
      strccat(CMD," %s %s",argv[i],argv[i+1]);
      mml = atoi(argv[++i]);
    } else if (strcmp(argv[i],"-col") == 0 && CheckValOpt(argv+i)) {
      strccat(CMD," %s %s",argv[i],argv[i+1]);
      RsatMemTracker = relem( (void*)col, RsatMemTracker);
      col = argv[++i];
    } else if (strcmp(argv[i],"-v") == 0 && CheckValOpt(argv+i)) {
      strccat(CMD," %s %s",argv[i],argv[i+1]);
      verbose = atoi(argv[++i]);
    } else {
      RsatFatalError("Invalid option", argv[i], NULL);
    }
  }

  /////////////////////////////////////////////////
  //Initialize script
  /////////////////////////////////////////////////
  start_time = StartScript(PROGRAM->buffer,CMD->buffer);

  //Allocate memory for strings
  line                = strnewToList(&RsatMemTracker);
  genome_dir          = strnewToList(&RsatMemTracker);
  variant_dir         = strnewToList(&RsatMemTracker);
  outfile_var_sorted  = strnewToList(&RsatMemTracker);
  get_variations_cmd  = strnewToList(&RsatMemTracker);
  sort_variations_cmd = strnewToList(&RsatMemTracker);
  deletetmps_cmd      = strnewToList(&RsatMemTracker);
  /////////////////////////////////////////////////
  // Validate Arguments
  /////////////////////////////////////////////////
  if (!species) RsatFatalError("No species specified. Use -species",NULL);
  if (!(assembly || release)) RsatFatalError("No assembly and ensembl version specified. Use at least one of these options: -release -assembly",NULL);
  if (!format) RsatFatalError("No input format specified. Use -format",NULL);

  // Retrieve and check genome and variation directories
  Get_genome_dir(genome_dir,species,assembly,release,species_suffix);
  Get_variation_dir(variant_dir,species,assembly,release,species_suffix);

  //QUESTION WSG(2017-06-18). stat() needs to have execute permissions on all path folders
  if ( !(stat(genome_dir->buffer,  &dir_exists) == 0 && S_ISDIR(dir_exists.st_mode)) ) RsatFatalError("Genome directory" ,  genome_dir->buffer, "does not exists or granted permissions were not properly set. Use download-ensembl-variation before retrieve-variation-seq or check for access/execution permissions.",NULL);
  if ( !(stat(variant_dir->buffer, &dir_exists) == 0 && S_ISDIR(dir_exists.st_mode)) ) RsatFatalError("Variation directory",variant_dir->buffer,"does not exists or granted permissions were not properly set. Use download-ensembl-variation before retrieve-variation-seq or check for access/execution permissions.",NULL);

  //Check if chromosome files are not missing
  trieChrom = Get_file_seq_name(genome_dir);

  /////////////////////////////////////////////////
  // Print first part of header
  /////////////////////////////////////////////////

  //Open outfile
  if (output) fout = OpenOutputFile(fout,output);
  //Print first part of header
  printf("; %s\n",CMD->buffer);

  fprintf(fout,"; %s\n",CMD->buffer);
  if(input)  fprintf(fout,"; Input  file\n; \t%s\n",input);
  if(output) fprintf(fout,"; Output file\n; \t%s\n",output);


  //Check for stdin and write it down to a file at tmp in $RSAT
  if(!input){
    //Allocate memory for string
    outfile_stdin = strnewToList(&RsatMemTracker);

    //Create name ok tmp file
    make_temp_file(outfile_stdin,"","variation_stdin_rsat",1,0,0);
    strccat(outfile_stdin,".%s",format);

    //Write stdin to file
    fh_stdin = OpenOutputFile(fh_stdin,outfile_stdin->buffer);
    while ( fread(&character,1,1,fin) == 1){
      if ( fwrite(&character,1,1,fh_stdin) != 1)
        RsatFatalError("Unable to read properly from stdin in main()",NULL);
    }
    fclose(fh_stdin);

    //Input is now file at .tmp in $RSAT
    input = outfile_stdin->buffer;
  }

  //Initialize cmd if format eq to bed or id
  if (strcmp(format,"bed") == 0 || strcmp(format,"id") == 0) {
    //Allocate memory for string
    outfile_var_info   = strnewToList(&RsatMemTracker);

    strfmt(get_variations_cmd,"%s/perl-scripts/variation-info -i %s -species %s",RSAT->buffer,input,species);
    if(release)  strccat(get_variations_cmd," -release %s",release);
    if(assembly) strccat(get_variations_cmd," -assembly %s",assembly);
  }

  /////////////////////////////////////////////////
  // Retrieve variations from varBed format
  /////////////////////////////////////////////////
  if (strcmp(format,"varBed") == 0) {
    fin = OpenInputFile(fin,input);
    while( fread( (line->buffer + line->size) , 1, 1, fin) == 1 ) {
      //Substitute '\s' for '\0'
      if(line->buffer[line->size] == ' '){
        line->buffer[line->size] = '\0';
      }
      //Test for end of line
      if(line->buffer[line->size] == '\n'){
        //Reinitialize parameters for line reading
        line->buffer[line->size] = '\0';

        //Skip lines
        if(line->buffer[0] == '#') break;
        if(line->buffer[0] == ';'){
          //This parse depends heavily on the varBed header file
          if( strcmp( line->buffer + 2,"Phased" ) == 0){
            phased = (strcmp( line->buffer + line->size - 4, "True" ) == 0) ? 1 : 0;
            if(phased && verbose >= 12) RsatInfo("Found Status Phased = True",NULL);
            break;
          }
        }
        //Reinitialize counter
        line->size = 0;
        continue;
      }
      //Add +1 to current line size
      line->size++;
      //Resize line string if limit has reached
      strlimt(line);
    }
    //Clean line size string and close fh
    line->size = 0;
    fclose(fin);

  }
  /////////////////////////////////////////////////
  // Retrieve variations from id format
  /////////////////////////////////////////////////
  else if (strcmp(format,"id") == 0) {
    //Create temporary file
    make_temp_file(outfile_var_info,"","variation_rsat_fromIDs", 1,0,0);
    strccat(outfile_var_info,".varBed");

    //get_variations command completed for id
    strccat(get_variations_cmd," -format id -col %s -o %s",col,outfile_var_info->buffer);

    doit(get_variations_cmd->buffer,0,0,1,0,NULL,NULL,NULL,NULL);
    input = outfile_var_info->buffer;

  }
  /////////////////////////////////////////////////
  // Retrieve variations from BED format
  /////////////////////////////////////////////////
  else if (strcmp(format,"bed") == 0){
    //Create temprary file
    make_temp_file(outfile_var_info,"","variation_rsat_from_bed",1,0,0);
    strccat(outfile_var_info,".varBed");

    //get_variations command completed for bed file
    strccat(get_variations_cmd," -format bed -o %s",outfile_var_info->buffer);

    //doit
    doit(get_variations_cmd->buffer,0,0,1,0,NULL,NULL,NULL,NULL);
    input = outfile_var_info->buffer;
  }
  //If format file is not recognized,raise FatalError
  else {
    RsatFatalError("Format",format,"is not a valid format. Please use any of these: varBed,id or bed",NULL);
  }

  //Sort file
  make_temp_file(outfile_var_sorted,"","variation_sorted_rsat",1,0,0);
  strccat(outfile_var_sorted,".varBed");
  strfmt(sort_variations_cmd,"grep -v '^;\\|^#' %s | sort -k1,1 -k2,2n > %s",input, outfile_var_sorted->buffer);
  doit(sort_variations_cmd->buffer,0,0,1,0,NULL,NULL,NULL,NULL);
  input = outfile_var_sorted->buffer;

  /////////////////////////////////////////////////
  // Print last part of header
  /////////////////////////////////////////////////
  printHeader(phased,fout);

  /////////////////////////////////////////////////
  // Retrieve sequences
  /////////////////////////////////////////////////
  char **token      = NULL;
  char *alt_allele  = NULL;

  //Allocate memory for variables
  token = getokens(12);
  curr_chr = strnewToList(&RsatMemTracker);
  seq_file = strnewToList(&RsatMemTracker);

  //Initialize variables
  i  = 0;
  j  = 0;
  inputVars = 0;
  strcopy(curr_chr,"");
  fin = OpenInputFile(fin,input);
  token[0] = line->buffer;
  if (phased) {

    //Declare variables
    string  *varCoords            = NULL;
    string  *IDs                  = NULL;
    string  *SOs                  = NULL;
    string  *alleleFreqs          = NULL;
    string  *Haplotype1           = NULL;
    string  *Haplotype2           = NULL;
    string  *Haplotype1Sequence   = NULL;
    string  *Haplotype2Sequence   = NULL;

    variant *HaploGroup           = NULL;
    variant *firstVar             = NULL;
    variant *lastVar              = NULL;
    variant *printVar             = NULL;

    int isChrDiff                 =    0;

    //Allocate memory for variables
    varCoords           = strnewToList(&RsatMemTracker);
    IDs                 = strnewToList(&RsatMemTracker);
    SOs                 = strnewToList(&RsatMemTracker);
    alleleFreqs         = strnewToList(&RsatMemTracker);
    Haplotype1          = strnewToList(&RsatMemTracker);
    Haplotype2          = strnewToList(&RsatMemTracker);
    Haplotype1Sequence  = strnewToList(&RsatMemTracker);
    Haplotype2Sequence  = strnewToList(&RsatMemTracker);


    while( fread( (line->buffer + line->size) , 1, 1, fin) == 1 ) {

      if (line->buffer[line->size] == '\t' && j < 12) {
        token[++j] = line->buffer + line->size + 1;
           line->buffer[line->size] = '\0';
      }

      if(line->buffer[line->size] == '\n'){
        //Reinitialize parameters for line reading
        line->buffer[line->size] = '\0';
        line->size = 0;
        j = 0;

        ///////////////
        //Skip lines
        if(line->buffer[line->size] == '#')  continue;
        if(line->buffer[line->size] == ';')  continue;
        if(line->buffer[line->size] == '\0') continue;

        inputVars++;
        ///////////////
        //Load fasta
        if(strncmp(token[0],"chr",3) == 0) token[0] = token[0] + 3;
        //Test if new chromsome is different from previous and load file
        isChrDiff = (strcmp(token[0],curr_chr->buffer) != 0) ? 1 : 0;
        if( isChrDiff ){
          seq_search = TrieSearch(trieChrom,token[0]);
          if(seq_search == NULL){
            RsatWarning("Unable to locate file for this chr",token[0],"at",genome_dir->buffer,".Skipping line.",NULL);
            token[0] = line->buffer;
            continue;
          }
          //Process remaining variants from last chromosome
          while ( HaploGroup != NULL )   {
            //printf("LINE1 BEFORE \n");
            //printf("This is firstVar,HaploGroup,lastVar and printVar : %p, %p, %p and %p\n", (void*)firstVar, (void*)HaploGroup, (void*)lastVar, (void*)printVar );
            //Process Haplotypes
            processHaplotypes(mml, firstVar, HaploGroup, printVar, printVar, sequence, fout,
                              varCoords, IDs, SOs, alleleFreqs, Haplotype1, Haplotype2,
                              Haplotype1Sequence, Haplotype2Sequence);
            //printf("This is firstVar,HaploGroup,lastVar and printVar : %p, %p, %p and %p\n", (void*)firstVar, (void*)HaploGroup, (void*)lastVar, (void*)printVar );
            //printf("LINE1 AFTER\n");
           //Continue to next variant at center
           HaploGroup = HaploGroup->next;
           //Free unnecessary variants
           /*printf("This is RsatMemTracker pointer %p\n", (void*)RsatMemTracker );
           printf("This is RsatMemTracker mem %p\n", RsatMemTracker->mem );
           printf("This is RsatMemTracker id %d\n",RsatMemTracker->id );
           printf("This is RsatMemTracker var start %s-%s\n",((variant*)RsatMemTracker->mem)->start->buffer, ((variant*)RsatMemTracker->mem)->end->buffer);*/
           if (HaploGroup == NULL) RsatMemTracker= relem((void*)firstVar,RsatMemTracker);

         }
         //Load new chromosome
         strfmt(seq_file,"%s/%s",genome_dir->buffer,seq_search);
         if(sequence) free(sequence); //NOTE WSG.Pending to update
         sequence = Get_sequence(seq_file->buffer);
         stat(seq_file->buffer,&file);
         sequence_maxsize = (long long)file.st_size;
         strcopy(curr_chr,token[0]);
        }
        if( CheckOutOfIndex( token[0],token[1], token[2], mml, sequence_maxsize ) != 1 ) continue;
        //Test if alleles are in '-' strand and convert them to '+'
        if (token[3][0] == '-'){
          if(switch_strand(token[6]) == 0){
            RsatWarning("This is not a valid allele at",token[0],token[1],token[2],"Skipped.",NULL);
            continue;
          }
        } else if (token[3][0] != '+') {
          RsatWarning("Strand information does not match any know annotation.Skipped.",NULL);
          continue;
        }
        //Process HaploGroups
        if ( isChrDiff ) {
          HaploGroup          = varnewToList(&RsatMemTracker);
          varfill(HaploGroup, token[0], token[1], token[2], token[3], token[4], token[7], token[5], token[6], token[9]);
          firstVar = HaploGroup;
          /*printf("FirstVar This is Variant pointer %p\n", (void*)firstVar);
          printf("FirstVar This is Variant ID %s\n", firstVar->id->buffer);
          printf("FirstVar This is Variant Chr %s\n", firstVar->chromosome->buffer);
          printf("FirstVar This is Variant start %s\n", firstVar->start->buffer);
          printf("FirstVar This is Variant end %s\n", firstVar->end->buffer);
          printf("FirstVar This is Variant strand %s\n", firstVar->strand);
          printf("FirstVar This is Variant id %s\n", firstVar->id->buffer);
          printf("FirstVar This is Variant so %s\n", firstVar->SO->buffer);
          printf("FirstVar This is Variant ref %s\n", firstVar->reference->buffer);
          printf("FirstVar This is Variant all %s\n", firstVar->alleles->buffer);
          printf("FirstVar This is Variant freq %s\n", firstVar->freq->buffer);
          printf("FirstVar This is Variant prev %p\n", firstVar->prev);
          printf("FirstVar This is Variant next %p\n\n", firstVar->next);*/
          continue;
        } else {
          lastVar = varadd(HaploGroup);
          varfill(lastVar,token[0], token[1], token[2], token[3], token[4], token[7], token[5], token[6], token[9]);
          /*printf("LastVar This is Variant pointer %p\n", (void*)lastVar);
          printf("LastVar This is Variant ID %s\n", lastVar->id->buffer);
          printf("LastVar This is Variant Chr %s\n", lastVar->chromosome->buffer);
          printf("LastVar This is Variant start %s\n", lastVar->start->buffer);
          printf("LastVar This is Variant end %s\n", lastVar->end->buffer);
          printf("LastVar This is Variant strand %s\n", lastVar->strand);
          printf("LastVar This is Variant id %s\n", lastVar->id->buffer);
          printf("LastVar This is Variant so %s\n", lastVar->SO->buffer);
          printf("LastVar This is Variant ref %s\n", lastVar->reference->buffer);
          printf("LastVar This is Variant all %s\n", lastVar->alleles->buffer);
          printf("LastVar This is Variant freq %s\n", lastVar->freq->buffer);
          printf("LastVar This is Variant prev %p\n", lastVar->prev);
          printf("LastVar This is Variant next %p\n\n", lastVar->next);*/
          if( !isChrDiff && (atoi(lastVar->start->buffer) < atoi(lastVar->prev->end->buffer))  ) {
            RsatFatalError("End is bigger than Start, this is not a valid coordinate.",NULL);
          }
          if ( isChrDiff || ( mml < (atoi(lastVar->start->buffer) - atoi(HaploGroup->end->buffer)) ) ) {
            while (  ( mml < (atoi(lastVar->start->buffer) - atoi(HaploGroup->end->buffer)) ) )   {
              //printf("LINE2 BEFORE \n");
              //printf("This is firstVar,HaploGroup,lastVar and printVar : %p, %p, %p and %p\n", (void*)firstVar, (void*)HaploGroup, (void*)lastVar, (void*)printVar );
              //Process Haplotypes
              processHaplotypes(mml, firstVar, HaploGroup, lastVar, printVar, sequence, fout,
                                varCoords, IDs, SOs, alleleFreqs, Haplotype1, Haplotype2,
                                Haplotype1Sequence, Haplotype2Sequence);
              //printf("This is firstVar,HaploGroup,lastVar and printVar : %p, %p, %p and %p\n", (void*)firstVar, (void*)HaploGroup, (void*)lastVar, (void*)printVar );
              //printf("LINE2 AFTER \n");
              //Continue to next variant at center
              HaploGroup = HaploGroup->next;
              //Free unnecessary variants
              variant *tmp = firstVar;
              if ( HaploGroup == lastVar ) {
                do {
                  firstVar = firstVar->next;
                  varfree(firstVar->prev);
                } while(firstVar != HaploGroup);
                HaploGroup->prev = NULL;
                MemTrackUpd((void*)tmp, (void *)HaploGroup, RsatMemTracker);
                break;
              }

            }

          } else {
            continue;
          }
        }
        continue;

      }

      //Resize line string if limit has reached
      limlinetok(line,token,12);

      //Add +1 to current line size
      line->size++;
    }
    //Process last variants in list
    while ( HaploGroup != NULL )   {
      //printf("LINE3 BEFORE\n");
      //Process Haplotypes
      processHaplotypes(mml, firstVar, HaploGroup, printVar, printVar, sequence, fout,
                        varCoords, IDs, SOs, alleleFreqs, Haplotype1, Haplotype2,
                        Haplotype1Sequence, Haplotype2Sequence);
      //printf("LINE3 AFTER\n");
     //Continue to next variant at center
     HaploGroup = HaploGroup->next;
   }

  } else {

    while( fread( (line->buffer + line->size) , 1, 1, fin) == 1 ) {

      if (line->buffer[line->size] == '\t' && j < 12) {
        token[++j] = line->buffer + line->size + 1;
           line->buffer[line->size] = '\0';
      }

      if(line->buffer[line->size] == '\n'){
        //Reinitialize parameters for line reading
        line->buffer[line->size] = '\0';
        line->size = 0;
        j = 0;

        ///////////////
        //Skip lines
        if(line->buffer[line->size] == '#')  continue;
        if(line->buffer[line->size] == ';')  continue;
        if(line->buffer[line->size] == '\0') continue;

        ///////////////
        //Load fasta
        if(strncmp(token[0],"chr",3) == 0) token[0] = token[0] + 3;
        //Test if new chromsome is different from previous and load file
        if(strcmp(token[0],curr_chr->buffer) != 0){
          seq_search = TrieSearch(trieChrom,token[0]);
          if(seq_search == NULL){
            RsatWarning("Unable to locate file for this chr",token[0],"at",genome_dir->buffer,".Skipping line.",NULL);
            token[0] = line->buffer;
            continue;
          }
          strfmt(seq_file,"%s/%s",genome_dir->buffer,seq_search);
          if(sequence) free(sequence); //NOTE WSG.Pending to update
          sequence = Get_sequence(seq_file->buffer);
          stat(seq_file->buffer,&file);
          sequence_maxsize = (long long)file.st_size;
          strcopy(curr_chr,token[0]);
        }
        if( CheckOutOfIndex( token[0],token[1], token[2], mml, sequence_maxsize ) != 1 ) continue;
        ////////////////////////////////////////
        //Retrieve sequences for each variant
        k = 0;
        left_flank  = atoi(token[1]) - mml -1;
        right_flank = atoi(token[2]);

        //Test if alleles are in '-' strand and convert them to '+'
        if (token[3][0] == '-'){
          if(switch_strand(token[6]) == 0){
            RsatWarning("This is not a valid allele at",token[0],token[1],token[2],"Skipped.",NULL);
            continue;
          }
        } else if (token[3][0] != '+') {
          RsatWarning("Strand information does not match any know annotation.Skipped.",NULL);
          continue;
        }
        //ALT alleles
        alt_allele = token[6];
        do {
          for (i = 0; alt_allele[i] != '\0' ; i++) {
            if ( alt_allele[i] == ',') {
              alt_allele[i] = '\0';
              break;
            }
          }
          fprintf(fout,"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t",token[0],token[1],token[2],token[3],token[4],token[7],token[5],alt_allele,token[9]);
          for (k = left_flank; k < left_flank + mml; k++) {
            fprintf(fout,"%c",tolower(sequence[k]));
          }
          fprintf(fout,"%s", alt_allele);
          for (k = right_flank; k < right_flank + mml; k++) {
            fprintf(fout,"%c",tolower(sequence[k]));
          }
          fprintf(fout,"\n");
          alt_allele = alt_allele + i + 1;
        } while( alt_allele  != token[7] );

        //REF allele
        fprintf(fout,"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t",token[0],token[1],token[2],token[3],token[4],token[7],token[5],token[5],token[9]);
        for (k = left_flank; k < left_flank + mml; k++) {
          fprintf(fout,"%c",tolower(sequence[k]));
        }
        fprintf(fout,"%s", token[5]);
        for (k = right_flank; k < right_flank + mml; k++) {
          fprintf(fout,"%c",tolower(sequence[k]));
        }
        fprintf(fout,"\n");
        continue;

      }

      //Resize line string if limit has reached
      limlinetok(line,token,12);

      //Add +1 to current line size
      line->size++;
    }

  }


  //Free raw sequence memory
  free(sequence);
  //TrieEnd(trieChrom);

  //Close file handlers
  fclose(fin);
  fclose(fout);

  //Remove tmp sorted file
  strfmt(deletetmps_cmd,"rm %s",outfile_var_sorted->buffer);
  doit(deletetmps_cmd->buffer,0,0,0,0,NULL,NULL,NULL,NULL);

  //Update execution log files
  ReportExecutionTime(start_time);

  //Free all objects
  rlist(RsatMemTracker);

    //////////////////////
    ///TODO WSG. Remove!!!
    //fclose(fin);
    //fclose(fout);
    //rlist(RsatMemTracker);
    //printf("The end!\n");
    exit(0);

    ///
    //////////////////////
  return 0;
}

///////////////////////////////////////
///////////////////////////////////////
///////////////////////////////////////
///////////////////////////////////////
void argToList(stringlist *list, char **value){
  //Declare variables
  stringlist *tmp = NULL;

  //Add information to string list
  if( strcmp(list->element->buffer,"") == 0 ) {
     strcopy(list->element, value[1]);
  } else {
     tmp = strlistadd(list);
     strcopy(tmp->element, value[1]);
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
    fclose(fh_popen);
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
    if (error == -1) strfmt(error_message,"Could not execute the command\n\t%s",command);
    else{
      strfmt(error_message,"Error\t%d\tocurred during execution of the command:\n\t%s",error,command);
    }


    //char *error_message   = NULL;
    //if( (error_message    = (char *)malloc(sizeof(char) * BASE_STR_LEN)) == NULL) RsatFatalError("Unable to allocate memory in make_temp_file()",NULL);

    //if(error == -1) sprintf(error_message,"Could not execute the command\n\t%s",command);
    //else{
    //  sprintf(error_message,"Error\t%d\tocurred during execution of the command:\n\t%s",error,command);
    //}

    if(err_handle) fprintf(err_handle, "\n\n%s\n",error_message->buffer);

    if(die_on_error) RsatFatalError(error_message->buffer,NULL);
    else{
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
int CheckOutDir(string *output_dir,int Umask, int Chmod){
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
  if ( CheckOutDir(real_tmp_dir,0,755) == 0 ) RsatFatalError("Unable to create",real_tmp_dir->buffer,"in make_temp_file()",NULL);

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
    tmp_file->size++;
    strlimt(tmp_file);
  }
  //Add '\0'-end to string and +1 to size
  tmp_file->buffer[tmp_file->size] = '\0';
  tmp_file->size++;
  //Close filehandler
  fclose(fh_popen);

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
        printf("1HOLA!!!!\n" );
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
    //Check if ONLY assembly was passed as query option in order to compare properly
    } else if (assembly) {

      if (strcmp(token[1],species) == 0 && strcmp(token[2],assembly) == 0) {
        //Test if species_directory field is empty
        if (token[6] != '\0') {
          if( (rsat_idx = strchr(token[6],'}')) != NULL ) {
            strfmt( species_dir, "%s%s", RSAT->buffer, (rsat_idx + 1) );
          } else {
            RsatFatalError("Get_species_dir_from_supported_file() could not identify species",
                            species,"from assembly",assembly,
                            "in the organism table,species_dir field was not correct\n",NULL);
          }
        } else {
          RsatFatalError("Get_species_dir_from_supported_file() could not identify species",
                          species,"from assembly",assembly,
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


///////////////////////////////////////
// Exclusive functions program
void printHeader(int PhasedFile, FILE *fh_outputFile){
  if (PhasedFile) {
    fprintf(fh_outputFile,
            "; Column contents\n"
            "; chr          - Chromosome of Haplotype\n"
            "; start        - Start position of Haplotype\n"
            "; end          - End position of Haplotype\n"
            "; strand       - Strand of Haplotype\n"
            "; var_coord    - Comma-separated list of individual variation coordinates\n"
            "; id           - Comma-separated identifier list of the individual variations\n"
            "; soterm       - Comma-separated sequence ontology (SO) term describing each variation type\n"
            "; hap1         - Comma-separated list from the individual alleles comprising Haplotype 1\n"
            "; hap          - Comma-separated list from the individual alleles comprising the current Haplotype\n"
            "; allele_freq  - Comma-separated list from the individual allele frequencies\n"
            "; seq          - Sequence of the current Haplotype\n"
            "#chr\tstart\tend\tstrand\tvar_coord\tid\tsoterm\thap1\thap\tallele_freq\tseq\n");
  } else {
    fprintf(fh_outputFile,
            "; Column contents\n"
            "; chr          - Chromosome\n"
            "; start        - Start position\n"
            "; end          - End position\n"
            "; strand       - Strand\n"
            "; id           - Identifier\n"
            "; soterm       - Sequence ontology (SO) term describing variation type\n"
            "; ref_var      - Reference variant\n"
            "; alt_var      - Alternative variant\n"
            "; allele_freq  - Frequency of the current allele\n"
            "; seq          - Sequence of the current variant, flanked by a user-specified neighbouring region\n"
            "#chr\tstart\tend\tstrand\tid\tsoterm\tref_var\talt_var\tallele_freq\tseq\n");
  }
  return;
}


void processHaplotypes(int mml, variant *firstVar, variant *HaploGroup, variant *lastVar, variant *printVar, char *sequence, FILE *fout,
                       string *varCoords, string *IDs, string *SOs, string *alleleFreqs, string *Haplotype1, string *Haplotype2,
                       string *Haplotype1Sequence, string *Haplotype2Sequence ) {
  //Declare variables
  long long i = 0;
  long long left_flank     =          0;
  long long right_flank    =          0;
  long long tmp_left_flank = left_flank;

  //NOTE WSG. I suppressed -1 because start is 0-based and 1-terminated
  left_flank     = atoi(HaploGroup->start->buffer) - mml;
  right_flank    = atoi(HaploGroup->end->buffer)   + mml;
  tmp_left_flank = left_flank;

  for (printVar = firstVar; printVar != lastVar; printVar = printVar->next) {
    if (printVar == firstVar) {
      strfmt(varCoords   , "%s:%s-%s", printVar->chromosome->buffer,printVar->start->buffer,printVar->end->buffer);
      strfmt(IDs         , "%s"      , printVar->id->buffer         );
      strfmt(SOs         , "%s"      , printVar->SO->buffer         );
      strfmt(alleleFreqs , "%s"      , printVar->freq->buffer       );
      strfmt(Haplotype1  , "%s"      , printVar->reference->buffer  );
      strfmt(Haplotype2  , "%s"      , printVar->alleles->buffer    );
      Haplotype1Sequence->size = 0;
      Haplotype2Sequence->size = 0;
        //printf("1.-Did it succeeded?\n" );
        //printf("This is the printVar, firstVar, HaploGroup and lastVar :%p, %p, %p and %p\n", (void*)printVar,(void*)firstVar,(void*)HaploGroup,(void*)lastVar);
        //printf("This is Haplotype1Sequence %s\n", Haplotype1Sequence->buffer );
        //printf("This is Haplotype2Sequence %s\n\n", Haplotype2Sequence->buffer );


    } else {
      strccat(varCoords  , ",%s:%s-%s", printVar->chromosome->buffer,printVar->start->buffer,printVar->end->buffer);
      strccat(IDs        , ",%s"      , printVar->id->buffer        );
      strccat(SOs        , ",%s"      , printVar->SO->buffer        );
      strccat(alleleFreqs, ",%s"      , printVar->freq->buffer      );
      strccat(Haplotype1 , ",%s"      , printVar->reference->buffer );
      strccat(Haplotype2 , ",%s"      , printVar->alleles->buffer   );
    }
    for ( i = tmp_left_flank; i < atoi(printVar->start->buffer); i++ ) {
        //Write nucleotide to strings
        /*printf("2.-Did it succeeded?\n" );
        printf("This is the firstVar, HaploGroup and lastVar : %p, %p and %p\n", (void*)firstVar,(void*)HaploGroup,(void*)lastVar);
        printf("This is char %c on sequence[%lld]\n\n", sequence[i], i);*/
      Haplotype1Sequence->buffer[Haplotype1Sequence->size] = tolower(sequence[i]);
      Haplotype2Sequence->buffer[Haplotype2Sequence->size] = tolower(sequence[i]);
      //Update strings sizes
      Haplotype1Sequence->size++;
      Haplotype2Sequence->size++;
      //printf("This is current size and length: %d && %d and last char %c\n", (int)Haplotype2Sequence->size, (int)Haplotype2Sequence->length,Haplotype2Sequence->buffer[Haplotype2Sequence->size - 1]);
      //Test if strings have reached their limits
      strlimt(Haplotype1Sequence);
      strlimt(Haplotype2Sequence);
    }
    Haplotype1Sequence->buffer[Haplotype1Sequence->size] = '\0';
    Haplotype2Sequence->buffer[Haplotype2Sequence->size] = '\0';
    //Update strings sizes
    Haplotype1Sequence->size++;
    Haplotype2Sequence->size++;
    //Test if strings have reached their limits
    strlimt(Haplotype1Sequence);
    strlimt(Haplotype2Sequence);
    //Add Hap1/Hap2 alleles to their respective sequences
    strccat(Haplotype1Sequence, "%s", printVar->reference->buffer );
    strccat(Haplotype2Sequence, "%s", printVar->alleles->buffer   );

    //Update strings sizes to overwrite the trailing '\0' character from
    //the strccat operation
    Haplotype1Sequence->size--;
    Haplotype2Sequence->size--;
    /*printf("This is Variant pointer %p\n", (void*)printVar);
    printf("This is Variant ID %s\n", printVar->id->buffer);
    printf("This is Variant Chr %s\n", printVar->chromosome->buffer);
    printf("This is Variant start %s\n", printVar->start->buffer);
    printf("This is Variant end %s\n", printVar->end->buffer);
    printf("This is Variant strand %s\n", printVar->strand);
    printf("This is Variant id %s\n", printVar->id->buffer);
    printf("This is Variant so %s\n", printVar->SO->buffer);
    printf("This is Variant ref %s\n", printVar->reference->buffer);
    printf("This is Variant all %s\n", printVar->alleles->buffer);
    printf("This is Variant freq %s\n", printVar->freq->buffer);
    printf("This is Variant prev %p\n", printVar->prev);
    printf("This is Variant next %p\n\n", printVar->next);*/

    tmp_left_flank = atoi(printVar->end->buffer);
    //printf("\nCHECKPOINT 1!\n");
  }

  for ( i = tmp_left_flank; i < right_flank; i++ ) {
    //Write nucleotide to strings
    Haplotype1Sequence->buffer[Haplotype1Sequence->size] = tolower(sequence[i]);
    Haplotype2Sequence->buffer[Haplotype2Sequence->size] = tolower(sequence[i]);
    //Update strings sizes
    Haplotype1Sequence->size++;
    Haplotype2Sequence->size++;
    //Test if strings have reached their limits
    strlimt(Haplotype1Sequence);
    strlimt(Haplotype2Sequence);
  }
  Haplotype1Sequence->buffer[Haplotype1Sequence->size] = '\0';
  Haplotype2Sequence->buffer[Haplotype2Sequence->size] = '\0';
  //Update strings sizes
  Haplotype1Sequence->size++;
  Haplotype2Sequence->size++;
  //Print content of both Haolotype lines
  fprintf( fout, "%s\t%lld\t%lld\t+\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
                 HaploGroup->chromosome->buffer,
                 left_flank,
                 right_flank,
                 varCoords->buffer,
                 IDs->buffer,
                 SOs->buffer,
                 Haplotype1->buffer,
                 Haplotype2->buffer,
                 alleleFreqs->buffer,
                 Haplotype2Sequence->buffer);
   fprintf( fout, "%s\t%lld\t%lld\t+\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
                  HaploGroup->chromosome->buffer,
                  left_flank,
                  right_flank,
                  varCoords->buffer,
                  IDs->buffer,
                  SOs->buffer,
                  Haplotype1->buffer,
                  Haplotype1->buffer,
                  alleleFreqs->buffer,
                  Haplotype1Sequence->buffer);

     //printf("\nIt entered!\n");
     //printf("This is Variant pointer %p\n", (void*)lastVar);
     //printf("This is Variant ID %s\n", lastVar->id->buffer);
     //printf("This is Variant Chr %s\n", lastVar->chromosome->buffer);
     //printf("This is Variant start %s\n", lastVar->start->buffer);
     //printf("This is Variant end %s\n", lastVar->end->buffer);
     //printf("This is Variant strand %s\n", lastVar->strand);
     //printf("This is Variant id %s\n", lastVar->id->buffer);
     //printf("This is Variant so %s\n", lastVar->SO->buffer);
     //printf("This is Variant ref %s\n", lastVar->reference->buffer);
     //printf("This is Variant all %s\n", lastVar->alleles->buffer);
     //printf("This is Variant freq %s\n", lastVar->freq->buffer);
     //printf("This is Variant prev %p\n", lastVar->prev);
     //printf("This is Variant next %p\n", lastVar->next);

     /*printf("HaploGroup This is Variant pointer %p\n", (void*)HaploGroup);
     printf("HaploGroup This is Variant ID %s\n", HaploGroup->id->buffer);
     printf("HaploGroup This is Variant Chr %s\n", HaploGroup->chromosome->buffer);
     printf("HaploGroup This is Variant start %s\n", HaploGroup->start->buffer);
     printf("HaploGroup This is Variant end %s\n", HaploGroup->end->buffer);
     printf("HaploGroup This is Variant strand %s\n", HaploGroup->strand);
     printf("HaploGroup This is Variant id %s\n", HaploGroup->id->buffer);
     printf("HaploGroup This is Variant so %s\n", HaploGroup->SO->buffer);
     printf("HaploGroup This is Variant ref %s\n", HaploGroup->reference->buffer);
     printf("HaploGroup This is Variant all %s\n", HaploGroup->alleles->buffer);
     printf("HaploGroup This is Variant freq %s\n", HaploGroup->freq->buffer);
     printf("HaploGroup This is Variant prev %p\n", HaploGroup->prev);
     printf("HaploGroup This is Variant next %p\n\n", HaploGroup->next);*/
  return;

}

int CheckOutOfIndex( char *chr,char *start, char *end, int mml, long long maxsize ){
  //Declare variables
  long long num_start = 0;
  long long num_end   = 0;

  //Convert to value
  num_start = atoi(start);
  num_end   = atoi(end);

  if( num_start < 0  ){
    RsatWarning("Unable to process current variant at", chr, start, end,"The start is out of range for the sequence",NULL);
    return 0;
  } else if ( num_end > (maxsize -1) ) {
    RsatWarning("Unable to process current variant at", chr, start, end,"The end is out of range for the sequence",NULL);
    return 0;
  } else if ( num_end   + mml >  (maxsize -1) ) {
    RsatWarning("Unable to process current variant at", chr, start, end,"The start sequence is out of range for the sequence and will overflow.",NULL);
    return 0;
  } else if ( num_start - mml < 0 ) {
    RsatWarning("Unable to process current variant at", chr, start, end,"The end is out of range for the sequence and will overflow.",NULL);
    return 0;
  }
  return 1;
}
