################################################################
## Calculate statistics on genome sizes, genics/intergenic proportions, protein sizes, ...

include ${RSAT}/makefiles/util.mk
MAKEFILE=${RSAT}/makefiles/genome_statistics.mk

ORGANISMS=`supported-organisms -format full | cut -f 1 | grep -v ";"`

#ORG=Mycoplasma_genitalium
#ORG=Acinetobacter_sp_ADP1
ORG=Saccharomyces_cerevisiae
GENOME_DIR=${RSAT}/data/genomes/${ORG}/genome
GENOME=${GENOME_DIR}/contigs.txt
GENOME_SIZE=`sequence-lengths -i ${GENOME} -format filelist -sum`
genome_size:
	@echo ${ORG} ${GENOME_SIZE}

FEATURES=${GENOME_DIR}/feature.tab

FEATURE_NB=`awk '$$2=="CDS"' ${FEATURES} | wc -l`
count_cds:
	@echo ${ORG} ${FEATURE_NB}

list_organisms:
	@echo ${ORGANISMS}

################################################################
## genic versus intergenic fraction for each genome
coding_or_not: coding_or_not_all_orgs coding_or_not_cat


## Calculate genome statistics for one organism
STAT_DIR=data/genome_stats
STAT_FILE_ONE_ORG=${STAT_DIR}/${ORG}_stats.tab
coding_or_not_one_org:
	@Mkdir -p ${STAT_DIR}
	@(cd ${STAT_DIR}; coding-or-not -org ${ORG} -return stats)
	@echo ${STAT_FILE_ONE_ORG}

## Calculate genome statistics for all organisms
STAT_FILE=${STAT_DIR}/stats_all_organisms
coding_or_not_all_orgs:
	${MAKE} iterate_organisms ORG_TASK=coding_or_not_one_org

## Merge all genome statistics in a single file
coding_or_not_cat:
	grep "^;" ${STAT_FILE_ONE_ORG} > ${STAT_FILE}.tab
	cat ${STAT_DIR}/*_stats.tab | grep -v "^;"  >> ${STAT_FILE}.tab
	text-to-html -i ${STAT_FILE}.tab -o ${STAT_FILE}.html -font variable -chunk 1000
	@echo ${STAT_FILE}.tab

################################################################
## Protein sizes

## Calculate protein size statistics for one organism
PROT_DIR=${RSAT}/data/genomes/${ORG}/genome
PROT_SEQ=${PROT_DIR}/${ORG}_aa.fasta
aa_lengths_one_org:
	@mkdir -p ${PROT_STAT_DIR}
	cat ${PROT_SEQ} | sequence-lengths -format fasta | cut -f 2 > ${LEN_FILE}.txt
	classfreq -ci 100 -v  -i ${LEN_FILE}.txt -o ${LEN_FILE}_distrib.tab
	stats -i ${LEN_FILE}.txt -o ${LEN_FILE}_stats.tab -table -add organism_name ${ORG}
	@echo ${LEN_FILE}_stats.tab

## Calculate protein size statistics for all organisms
aa_lengths_all_orgs:
	${MAKE} iterate_organisms ORG_TASK=aa_lengths_one_org

## Merge all protein sizes in a single file
CDS_FILE=${GENOME_DIR}/cds.tab
PROT_STAT_DIR=data/protein_stats
LEN_FILE=${PROT_STAT_DIR}/${ORG}_protein_lengths
PROT_STAT_FILE=${PROT_STAT_DIR}/protein_stats_all_organisms
aa_lengths_cat:
	grep "^;" ${LEN_FILE}_stats.tab > ${PROT_STAT_FILE}.tab
	cat ${PROT_STAT_DIR}/*_stats.tab | grep -v "^;"  >> ${PROT_STAT_FILE}.tab
	text-to-html -i ${PROT_STAT_FILE}.tab -o ${PROT_STAT_FILE}.html -font variable -chunk 1000
	@echo ${PROT_STAT_FILE}.tab

################################################################
## Perform all tasks
NB_ORGS=`${MAKE} list_organisms | wc -w`
all:
	@echo "Calculating genome and protein statistics for ${NB_ORGS} organisms"
	@${MAKE} aa_lengths_all_orgs
	@${MAKE} aa_lengths_cat
	@${MAKE} coding_or_not_all_orgs
	@${MAKE} coding_or_not_cat
