## Import functions to run R
import os
import itertools
from snakemake.utils import R

################################################################
## Environment variables
RSAT = os.environ['RSAT']

configfile: os.path.join(RSAT, "snakefiles", "workflows", "kmer_sequence_classification.json")


################################################################
ORG=config["genome"]["organism"]
#ORG="Brachypodium_distachyon"

RES_DIR = os.path.join(config["dir"]["results"] , ORG)
#RES_DIR = os.path.join("results", ORG)
SEQ_DIR = os.path.join(RES_DIR, "sequences")
SEQ_TYPE="upstream_TSS"
SEQ_LEN=2000
SEQ_PREFIX = "Brachypodium_distachyon_gene_flanking_region_" + SEQ_TYPE + "_" + str(SEQ_LEN) + "_unique"
SEQ_FILE = os.path.join(SEQ_DIR, SEQ_PREFIX + ".fasta")

POS_OL = 6 ## Oligonucleotide length for position-analysis
POS_CI = config["positions"]["window_size"] ## Window width for position-analysis
POS_STRAND = config["positions"]["strand"]
POS_SUFFIX = "positions_" + str(POS_OL) + "nt_ci" + str(POS_CI) + POS_STRAND
POS_PREFIX = SEQ_PREFIX + "_" + POS_SUFFIX
#POS_DIR = os.path.join(RES_DIR, "k-mers", "positions", POS_PREFIX)
POS_DIR = SEQ_DIR
POS_FILE = os.path.join(POS_DIR, POS_PREFIX + ".tab")

print("SEQ_DIR\t" + SEQ_DIR)
print("SEQ_FILE\t" + SEQ_FILE)
print("POS_SUFFIX\t" + POS_SUFFIX)
print("POS_DIR\t" + POS_DIR)
print("POS_PREFIX\t" + POS_PREFIX)
print("POS_FILE\t" + POS_FILE)

################################################################################################################################
################################################################################################################################

################################################################
## Define rules
rule positions:
    input: "{any_file}.fasta"
    output: "{any_file}_" +  POS_SUFFIX + ".tab"
    params: verbosity = config["verbosity"], \
            window_size = str(POS_CI), \
            origin=config["positions"]["origin"], \
            offset = config["positions"]["offset"], \
            k = str(POS_OL), \
            strand = POS_STRAND, \
            return_fields="occ,exp_occ,freq_per_window,freq_per_word,chi,sig,rank,graphs,html",
#            return_fields="distrib,occ,coverage,index,exp_occ,freq_per_window,freq_per_word,chi,sig,rank,graphs,clusters,matrices,html",
            opt=""
    shell: "position-analysis -v {params.verbosity} -i {input} -ci {params.window_size} -l {params.k} {params.strand} -last 200 -return {params.return_fields} -o {output}"

#########
## All
rule all:
     input: POS_FILE

# ##################################################################
# ## Fetch the corresponding sequences of agiven refrences genome
# ## using the coordinates (BED) file
# rule fetch_sequences_from_BED:
#     """  Given a coordinates (BED) file and a reference genome fetch and retrieve the corresponding sequences in fasta format.
#          Program: RSAT fetch-sequences
#     """
#     input:
#         BED_DIR + "/{promoters}.bed"
#     output:
#         SEQUENCES_DIR + "/{promoters}.fasta"
#     message:
#         " Fetching sequences from {wildcards.promoters}.bed in the reference genome {REFERENCE_GENOME}."
#     params:
#         ref_genome = REFERENCE_GENOME, \
#         v = "1", \
#         head_format = "galaxy"
#     shell:
#         'fetch-sequences -v {params.v} -i {input} -genome {params.ref_genome} -header_format {params.head_format} -o {output}' 
