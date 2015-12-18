"""Description: This workflow collects regulatory sequences for a set
of input genes (e.g. co-expression cluster, regulons, ...) and runs
various motif discovery approaches as well as functional enrichment
analysis (Gene Ontology processes).  In parallel, a random control is
performed by running exactly the same approaches on random gene
clusters of the same size as the original cluster.  


Authors: Bruno Contreras Moreira, Claire Rioualen, Jaime Castro
Mondragon, Jacques van Helden 

Tested with some maize clusters from
  http://www.ncbi.nlm.nih.gov/pubmed/25918418

Usage: 
    snakemake -p  -c "qsub {params.qsub}" -j 12 \
        -s ${RSAT}/snakefiles/workflows/gene-cluster-motifs.py \
        --configfile configfile \
        [targets]

Flowcharts:
    snakemake -p \
        -s ${RSAT}/snakefiles/workflows/gene-cluster-motifs.py \
        --configfile configfile \
        --force flowcharts

The configuration file must be specified with the snakemake opeion
--configfile. This requires a snakemake version >= 3.2.

"""

#================================================================#
#                        Imports                                 #
#================================================================#

from snakemake.utils import R
import os
import sys
import time
import datetime
import pandas as pd
import numpy as np

## Define verbosity
if not "verbosity" in config.keys():
    verbosity = 1

#================================================================#
#                         Includes                               #
#================================================================#

## Define the path to RSAT snakemake libraries
RSAT = os.environ["RSAT"]
RULES = os.path.join(RSAT, "snakefiles/rules")

## Inclute the rules to be used for this workflow
include: os.path.join(RULES, "retrieve_seq.rules")
include: os.path.join(RULES, "util.py")

#================================================================#
#                      Data & wildcards                          #
#================================================================#

#================================================================#
#           Check input and output files and directories         #
#================================================================#


################################################################
## Output directory
if not "outdir" in config.keys():
    sys.exit("Output directory must be defined in config file (option outdir).")
OUTDIR = config["outdir"]
if not os.path.exists(OUTDIR):
    time_warn("Creating output directory\t" + OUTDIR)
    os.makedirs(OUTDIR)

################################################################
## Input file describing the gene cluster composition. 
## This must be a tab-delimited file containing at least two columns:
##   - gene name / ID
##   - cluster name
## Additional columns are ignored.
if not "cluster_file" in config.keys():
    sys.exit("Output directory must be defined in config file (option cluster_file).")
## Read the count table for the current sample
cluster_desc = pd.read_csv(config["cluster_file"], sep="\t", header=None, index_col=False, names=["gene", "cluster"])

## Extract cluster names
cluster_names = np.unique(cluster_desc["cluster"])
time_warn("Cluster names: " + "; ".join(cluster_names))
for cluster in cluster_names:
    cluster_dir = os.path.join(OUTDIR,cluster)
    if not os.path.exists(cluster_dir):
        time_warn("Cluster directory\t" + cluster_dir)
        os.makedirs(cluster_dir)


## Extract the list of cluster names from the second column of the
## cluster table.


#================================================================#
#                         Workflow                               #
#================================================================#


#SPECIES=zea_mays
#ORTH_SPECIES=sorghum_bicolor

#OUTDIR="results/sequences"
#ORTH_SPECIES_FULL="Sorghum_bicolor.Sorbi1.29"
#ORTHIDENT=70
#RNDSAMPLES=50 # how many random replicates to be generated per input gene cluster
#MEMEMKORD=2   # for MEME bg


# SEQ=expand(os.path.join(OUTDIR, "{gene_cluster}", "{gene_cluster}.raw.rm.fna") ,gene_cluster = GENE_CLUSTERS )


## TEMPORARY
ONE_CLUSTER_PATH="~/marseille/protocolos/regulons/Zea_mays/E2F/regulonE2F"
ONE_CLUSTER_PATH="/home/rsat/test/gene-cluster-motifs_test/regulons/Zea_mays/E2F/regulonE2F"
SEQ=ONE_CLUSTER_PATH + config["retrieve_seq"]["suffix"] + ".fna"
rule all:
    input: SEQ

rule list_param:
    run:
        print("RSAT\t" + RSAT)
        print("RULES\t" + RULES)
        print("OUTDIR\t" + OUTDIR)
        print("SEQ\t" + SEQ)


################################################################

## Obtain upstream sequences of all gene_clusters, convert consensus and get orths
## uses iupac2meme from MEME suite
# .PHONY: get-seqs $(GENE_CLUSTERS)
# get-seqs: $(GENE_CLUSTER) 
# $(GENE_CLUSTER): list_param 
# 	@echo Retrieving upstream sequences of $@
# 	@retrieve-seq -org ${SPECIES} -feattype gene -from ${FROM} -to ${TO} \
# 		-noorf -i $@/gene_cluster$@.txt -label id -o $@/gene_cluster$@.raw.fna
# 	@retrieve-seq -org ${SPECIES} -feattype gene -from ${FROM} -to ${TO} \
# 		-noorf -rm -i $@/gene_cluster$@.txt -label id -o $@/gene_cluster$@.rm.fna
# 	@echo
# 	@echo Converting cognate consensus of $@ 
# 	@perl -ne 'print `iupac2meme $$_`' $@/consensus$@.txt > $@/consensus$@.meme
# 	@convert-matrix -i $@/consensus$@.meme -from meme -to tf -o $@/consensus$@.tf
# 	@echo
# 	@echo Retrieving orthologues of $@
# 	@get-orthologs-compara -i $@/gene_cluster$@.txt -ref_org ${ORTH_SPECIES} \
# 		-ident_query ${ORTHIDENT} -o $@/orthologues.txt
# 	@retrieve-seq -org ${ORTH_SPECIES_FULL} -feattype gene -from ${FROM} -to ${TO} \
# 		-noorf -i $@/orthologues.txt -label id -o $@/orthologues.raw.fna
# 	@retrieve-seq -org ${ORTH_SPECIES_FULL} -feattype gene -from ${FROM} -to ${TO} \
# 		-noorf -rm -i $@/orthologues.txt -label id -o $@/orthologues.rm.fna
# 	@cat $@/gene_cluster$@.raw.fna $@/orthologues.raw.fna > \
# 		$@/gene_cluster$@.orthologues.raw.fna
# 	@echo
# 	@echo Creating ${RNDSAMPLES} random clusters based on $@
# 	@for r in `seq 1 ${RNDSAMPLES}`; do \
# 		echo Replicate $$r; \
# 		NSEQS=`wc -l $@/gene_cluster$@.txt`; \
# 		RND_CLUSTER=$@/random$@-$$r.txt; \
#         RNDRAW=$@/random$@-$$r.raw.fna; \
#         RNDMSK=$@/random$@-$$r.rm.fna; \
#         random-genes -n $$NSEQS} -org ${SPECIES} -feattype gene -g 1 -o $$RND_CLUSTER; \
#         retrieve-seq -org ${SPECIES} -feattype gene -from ${FROM} -to ${TO} \
#         	-noorf -i $$RND_CLUSTER -label id -o $$RNDRAW; \
#         retrieve-seq -org ${SPECIES} -feattype gene -from ${FROM} -to ${TO} \
#         	-noorf -rm -i $$RND_CLUSTER -label id -o $$RNDMSK; \
# 	done
# ################################################################


# ## Data import & merging.

# IMPORT = expand(OUTDIR + "{samples}/{samples}.fastq", samples=SAMPLE_IDS) 

# ## Graphics & reports
# GRAPHICS = expand(OUTDIR + "dag.pdf")
# REPORT = expand(OUTDIR + "report.html")

# #----------------------------------------------------------------#
# # Quality control
# #----------------------------------------------------------------#

# RAW_QC = expand(OUTDIR + "{samples}/{samples}_fastqc/", samples=SAMPLE_IDS)
# RAW_READNB = expand(OUTDIR + "{samples}/{samples}_fastq_readnb.txt", samples=SAMPLE_IDS)

# #----------------------------------------------------------------#
# # Alignment<
# #----------------------------------------------------------------#

# ## to avoid duplicates, fasta sequence should be moved to {genome} directly...
# BWA_INDEX = expand(config["dir"]["genomes"] + "{genome}/BWAIndex/{genome}.fa.bwt", genome=GENOME)
# BOWTIE2_INDEX = expand(config["dir"]["genomes"] + "{genome}/Bowtie2Index/{genome}.fa.1.bt2", genome=GENOME)

# MAPPING = expand(OUTDIR + "{alignment}.sam", alignment=ALIGNMENT)

# # Sorted and converted reads (bam, bed)
# SORTED_MAPPED_READS_BWA = expand(OUTDIR + "{alignment}_sorted_pos.bam", alignment=ALIGNMENT)
# BAM_READNB = expand(OUTDIR + "{alignment}_sorted_pos_bam_readnb.txt", alignment=ALIGNMENT)
# SORTED_READS_BED = expand(OUTDIR + "{alignment}_sorted_pos.bed", alignment=ALIGNMENT)
# BED_FEAT_COUNT = expand(OUTDIR + "{alignment}_sorted_pos_bed_nb.txt", alignment=ALIGNMENT)

# # ----------------------------------------------------------------
# # Peak-calling
# # ----------------------------------------------------------------

# PEAKS = expand(OUTDIR + "{peakcalling}.bed", peakcalling=PEAKCALLING)

# # ----------------------------------------------------------------
# # Peak analysis
# # ----------------------------------------------------------------

# GET_FASTA = expand(OUTDIR + "{peakcalling}.fasta", peakcalling=PEAKCALLING)
# PURGE_PEAKS = expand(OUTDIR + "{peakcalling}_purged.fasta", peakcalling=PEAKCALLING)
# PEAKS_LENGTH = expand(OUTDIR + "{peakcalling}_purged_length.png", peakcalling=PEAKCALLING)
# PEAK_MOTIFS = expand(OUTDIR + "{motifs}_peak-motifs_synthesis.html", motifs=MOTIFS)

# ## Oligo analysis # ! missing f* input exception
# OLIGO = config['oligo_analysis']['count_oligo'].split()
# OLIGO_ANALYSIS = expand(OUTDIR + "{peakcalling}_purged_oligo{oligo}.txt", peakcalling=PEAKCALLING, oligo=OLIGO)

# #================================================================#
# #                        Rule all                                #
# #================================================================#

# rule all: 
# 	"""
# 	Run all the required analyses
# 	"""
# 	input: IMPORT, GRAPHICS, PEAK_MOTIFS#RAW_QC, BWA_INDEX, MAPPING, PEAKS, 
# 	params: qsub=config["qsub"]
# 	shell: "echo Job done    `date '+%Y-%m-%d %H:%M'`"

# ##================================================================#
# ##                          Report                                #
# ##================================================================#

# #NOW = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")

# ##rule report:
# ##    """
# ##    Generate a report with the list of datasets + summary of the results.
# ##    """
# ## see Scerevisiae report

# ##----------------------------------------------------------------#
# ## Build the report (including DAG and rulegraph flowcharts).
# #from snakemake.utils import report

# ## Bulleted list of samples for the report
# #SAMPLE_IDS_OL=report_numbered_list(SAMPLE_IDS)
# #RAW_READS_OL=report_numbered_list(IMPORT)
# #RAW_QC_OL=report_numbered_list(RAW_QC)

# #MAPPING_OL=report_numbered_list(MAPPING)

# #PEAKFILES_OL=report_numbered_list(PEAKS)

# ##	input: GRAPHICS, IMPORT, TRIMMED_READS_SICKLE, TRIMMED_QC, RAW_QC, MAPPED_READS_BWA, RAW_READNB, BAM_READNB, BED_READNB, PEAKS_MACS2, FETCH_MACS2_PEAKS, PURGE_MACS2_PEAKS #redundant for flowcharts

# #NOW = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")

# #rule report:
# #    """
# #    Generate a report with the list of datasets + summary of the results.
# #    """
# #    input:  dag=config["dir"]["reports"] + "dag.pdf", \
# #            dag_png=config["dir"]["reports"] + "dag.png", \
# #            rulegraph=config["dir"]["reports"] + "rule.pdf", \
# #            rulegraph_png=config["dir"]["reports"] + "rule.png"
# #    output: html=config["dir"]["reports"] + "report.html"
# #    run:
# #        report("""
# #        ===========================================
# #        ChIP-seq analysis - - - P.aeruginosa ParABS
# #        ===========================================
# #        
# #        :Date:                 {NOW}
# #        :Project:              P aeruginosa
# #        :Analysis workflow:    Claire Rioualen
# #        
# #        Contents
# #        ========
# #        
# #        - `Flowcharts`_
# #        - `Datasets`_
# #             - `Samples`_
# #             - `Raw reads`_
# #             - `Mapping`_
# #             - `Peaks`_
# #             - `QC reports`_

# #        -----------------------------------------------------

# #        Flowcharts
# #        ==========

# #        - Sample treatment: dag_
# #        - Workflow: rulegraph_

# #        .. image:: rulegraph.png

# #        -----------------------------------------------------

# #        Datasets
# #        ========
# #        
# #        Samples
# #        -------

# #        {SAMPLE_IDS_OL} 

# #        Raw reads 
# #        ---------

# #        {RAW_READS_OL}

# #        Mapping
# #        -------

# #        {MAPPING_OL}

# #        Peaks
# #        -----

# #        {PEAKFILES_OL}

# #        QC reports
# #        ----------

# #        {RAW_QC_OL}

# #        -----------------------------------------------------

# #        """, output.html, metadata="Claire Rioualen (claire.rioualen@inserm.fr)", **input)
