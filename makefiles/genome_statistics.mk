################################################################
## Calculate statistics on genome sizes

include ${RSAT}/makefiles/util.mk

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

coding_or_not_cat:
	grep "^;" ${STAT_FILE_ONE_ORG} > ${STAT_FILE}.tab
	cat ${STAT_DIR}/*_stats.tab | grep -v "^;"  >> ${STAT_FILE}.tab
	text-to-html -i ${STAT_FILE}.tab -o ${STAT_FILE}.html -font variable -chunk 1000
	@echo ${STAT_FILE}.tab

STAT_FILE=${STAT_DIR}/stats_all_organisms
coding_or_not_all_orgs:
	${MAKE} iterate_organisms ORG_TASK=coding_or_not_one_org

STAT_DIR=data/genome_stats
STAT_FILE_ONE_ORG=${STAT_DIR}/${ORG}_stats.tab
coding_or_not_one_org:
	@Mkdir -p ${STAT_DIR}
	@(cd ${STAT_DIR}; coding-or-not -org ${ORG} -return stats)
	@echo ${STAT_FILE_ONE_ORG}

################################################################
## Protein sizes
CDS_FILE=${GENOME_DIR}/cds.tab
PROT_STAT_DIR=data/protein_stats
LEN_FILE=${PROT_STAT_DIR}/${ORG}_protein_lengths
PROT_STAT_FILE=${PROT_STAT_DIR}/protein_stats_all_organisms
aa_lengths_cat:
	grep "^;" ${LEN_FILE}_stats.tab > ${PROT_STAT_FILE}.tab
	cat ${PROT_STAT_DIR}/*_stats.tab | grep -v "^;"  >> ${PROT_STAT_FILE}.tab
	text-to-html -i ${PROT_STAT_FILE}.tab -o ${PROT_STAT_FILE}.html -font variable -chunk 1000
	@echo ${PROT_STAT_FILE}.tab

aa_lengths_all_orgs:
	${MAKE} iterate_organisms ORG_TASK=aa_lengths_one_org


NCBI_DIR=${HOME}/downloads/ftp.ncbi.nih.gov/genomes
EU_PROT_DIR=${NCBI_DIR}/${ORG}
BACT_PROT_DIR=${NCBI_DIR}/Bacteria/${ORG}
aa_lengths_one_org:
	@if [ -d "${BACT_PROT_DIR}" ] ; then	\
		echo "${ORG}	prokaryote ${BACT_PROT_DIR}" ; \
		${MAKE} aa_lengths_one_org_2 PROT_DIR=${BACT_PROT_DIR} ; \
	fi
	@if [ -d "${EU_PROT_DIR}" ] ; then						 \
		echo "${ORG}	eukaryote ${EU_PROT_DIR}"; \
		${MAKE} aa_lengths_one_org_2 PROT_DIR='${EU_PROT_DIR}' ; \
	fi

PROT_DIR=${NCBI_DIR}/${ORG}
aa_lengths_one_org_2:
	@mkdir -p ${PROT_STAT_DIR}
#	grep -v '^--' ${CDS_FILE} | awk -F\t '{print $$15}' | sequence-lengths -format multi | cut -f 2 > ${LEN_FILE}.txt
	cat ${PROT_DIR}/*.faa | sequence-lengths -format fasta | cut -f 2 > ${LEN_FILE}.txt
	classfreq -ci 100 -v  -i ${LEN_FILE}.txt -o ${LEN_FILE}_distrib.tab
	stats -i ${LEN_FILE}.txt -o ${LEN_FILE}_stats.tab -table -add organism_name ${ORG}
	@echo ${LEN_FILE}_stats.tab