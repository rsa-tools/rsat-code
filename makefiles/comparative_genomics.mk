###############################################################
# Comparative genomics with Saccharomyces downstream sequences

include ${RSAT}/makefiles/util.mk
MAKEFILE=${RSAT}/makefiles/comparative_genomics.mk

################################################################
## Variables
V=1
SIDE=3
REF_ORG=Saccharomyces_cerevisiae

################################################################
## Perform the complete analysis
ORG=S_kudriavzevii
#ORG=S_bayanus
OTHER_ORGANISMS=S_paradoxus S_mikatae S_bayanus S_kluyveri S_castellii S_kudriavzevii
ORGANISMS=S_cerevisiae ${OTHER_ORGANISMS}

################################################################
## Iterators

## iterate over organisms
ORG_TASK=utr2wc_one_org

## Original files from SGD distribution
SEQ_LEN=500
GENOME_DIR=${RSAT}/downloads/genome-ftp.stanford.edu/pub/yeast/data_download/sequence/fungal_genomes/${ORG}
UTR_FASTA=${GENOME_DIR}/${SOURCE}/orf_dna/utr${SIDE}_${SEQ_LEN}.fasta.gz

################################################################
## Collect gene descriptions
DESCR_FILE=${ORF_DIR}/S_cerevisiae_ORF_descriptions.tab
descriptions:
#	cat ${RSAT}/data/genomes/${REF_ORG}/genome/feature.tab | grep CDS | cut -f 1,3,8 > ${DECR_FILE}
	gene-info -i ${ORF_FILE} -o ${DESCR_FILE} -org ${REF_ORG}


################################################################
## Convert UTR sequences to wcconsensus format
## (convenient: one line per sequence)
SOURCE=WashU
UTR_DIR=data/utr${SIDE}
SEQ_DIR=${UTR_DIR}/sequences
UTR_WC=${SEQ_DIR}/${ORG}_${SOURCE}_utr${SIDE}_${SEQ_LEN}.wc.gz
STRICTO=								\
	${SEQ_DIR}/S_cerevisiae*_utr${SIDE}_${SEQ_LEN}.wc.gz		\
	${SEQ_DIR}/S_paradoxus*_utr${SIDE}_${SEQ_LEN}.wc.gz		\
	${SEQ_DIR}/S_mikatae*_utr${SIDE}_${SEQ_LEN}.wc.gz		\
	${SEQ_DIR}/S_kudriavzevii*_utr${SIDE}_${SEQ_LEN}.wc.gz		\
	${SEQ_DIR}/S_bayanus*_utr${SIDE}_${SEQ_LEN}.wc.gz	

LATO=									\
	${STRICTO}							\
	${SEQ_DIR}/S_kluyveri*_utr${SIDE}_${SEQ_LEN}.wc.gz		\
	${SEQ_DIR}/S_castellii*_utr${SIDE}_${SEQ_LEN}.wc.gz	

UTR_FILES=${LATO}
PROXIMITY=lato

list_utr_files:
	@echo ${UTR_FILES}

utr2wc_one_org:
	@echo "Converting UTR to format wc for organism ${ORG}	${UTR_WC}"
	@mkdir -p ${SEQ_DIR}
	gunzip -c ${UTR_FASTA}						\
		| perl -pe 's|^\>\S+_Contig|\>Contig|'			\
		| perl -pe 's/^\>(\S+)\s/\>${ORG}_${SOURCE}_$$1_/'	\
		| perl -pe 's/^\>S_/\>/'				\
		| convert-seq -from fasta -to wc -lw 0			\
		> ${UTR_WC}
#	gzip -f ${UTR_WC}

################################################################
## Collect non-coding sequences from Saccharomyces cerevisiae
SCE_UTR_FILE=${SEQ_DIR}/S_cerevisiae_SGD_utr${SIDE}_${SEQ_LEN}.wc.gz
SIDE_NAME=downstream
MAKE_UP=${MAKE}   SIDE=5 SIDE_NAME=upstream SEQ_LEN=2000 FROM=-2000 TO=-1 SUB_SEQ_LEN=800 SUB_FROM=-800 SUB_TO=-1 
MAKE_DOWN=${MAKE} SIDE=3 SIDE_NAME=downstream SEQ_LEN=500 FROM=1 TO=500 SUB_SEQ_LEN=200 SUB_FROM=1 SUB_TO=200

## Execute one task (UP_TASK=...) with all the parameters for upstream sequence
up:
	${MAKE_UP} ${UP_TASK}

up_Sce:
	${MAKE_UP} utr_Sce 

down_Sce:
	${MAKE_DOWN} utr_Sce 

FROM=1
TO=${SEQ_LEN}
RETRIEVE_CMD=	retrieve-seq -org ${REF_ORG}			\
		-i ${ORF_FILE}						\
		-label gene						\
		-type ${SIDE_NAME}					\
		-from ${FROM} -to ${TO} -format fasta			\
		| perl -pe 's/^>(\S+)\s+/>cerevisiae_SGD_$$1_/'	\
		| convert-seq -from fasta -to wc -lw 0	\
		> ${SCE_UTR_FILE}
utr_Sce:
	@echo ${RETRIEVE_CMD}
	${RETRIEVE_CMD}
#	gzip -f ${SCE_UTR_FILE}
	@echo "S_cerevisiae downstream sequences saved in ${SCE_UTR_FILE}"

################################################################
## Convert downstream sequences for one alternate yeast species
utr2wc:
	@${MAKE} utr2wc_one_org ORG=S_bayanus SOURCE=WashU
	@${MAKE} utr2wc_one_org ORG=S_castellii SOURCE=WashU
	@${MAKE} utr2wc_one_org ORG=S_kluyveri SOURCE=WashU
	@${MAKE} utr2wc_one_org ORG=S_kudriavzevii SOURCE=WashU
	@${MAKE} utr2wc_one_org ORG=S_mikatae SOURCE=WashU

	@${MAKE} utr2wc_one_org ORG=S_bayanus SOURCE=MIT
	@${MAKE} utr2wc_one_org ORG=S_mikatae SOURCE=MIT
	@${MAKE} utr2wc_one_org ORG=S_paradoxus SOURCE=MIT

################################################################
## Select subsequences of a specific length
sub_seq_all_org:
	@${MAKE} iterate_organisms ORG_TASK=sub_seq_one_org
	@${MAKE} iterate_organisms ORG_TASK=purge_seq_one_org

SUB_FROM=1
SUB_TO=200
SUB_SEQ_LEN=200
SUB_UTR_WC_PREFIX=${SEQ_DIR}/${ORG}_utr${SIDE}_${SUB_SEQ_LEN}
SUB_UTR_WC=${SUB_UTR_WC_PREFIX}.wc.gz
SUB_UTR_WC_PURGED=${SUB_UTR_WC_PREFIX}_purged.wc.gz
BIG_UTR_WC=`ls -1 ${SEQ_DIR}/${ORG}_*_utr${SIDE}_${SEQ_LEN}.wc.gz`
sub_seq_one_org:
	@echo "sub-sequence	${ORG}	${BIG_UTR_WC}	${SUB_FROM}	${SUB_TO}	${SUB_UTR_WC}"
	sub-sequence -format wc -i ${BIG_UTR_WC} -from ${SUB_FROM} -to ${SUB_TO} -o ${SUB_UTR_WC}

purge_seq_one_org:
	@echo "Purging sequences for genome	${SUB_UTR_WC_PURGED}"
	purge-sequence -1str -format wc			\
		-ml ${PURGE_ML} -mis ${PURGE_MIS}	\
		-i ${SUB_UTR_WC}			\
		| convert-seq -from fasta -to wc	\
		> ${SUB_UTR_WC_PURGED}
	@echo "Purged sequences for genome	${SUB_UTR_WC_PURGED}"
	@echo

################################################################
## Merge redundant genomes (sequenced by both MIT and WashU
merge_sources:
	@mkdir -p ${REDUNDANT_DIR}
	@for org in S_bayanus S_mikatae; do				\
		${MAKE} discard_redundant_files ORG_TO_MERGE=$${org};	\
		${MAKE} merge_sources_one_org ORG_TO_MERGE=$${org};	\
		${MAKE} count_merged ORG_TO_MERGE=$${org};		\
	done

################################################################
## Merge the sequences from MIT and WashU for a given organism
ORG_TO_MERGE=S_bayanus
MERGED_FILE=${SEQ_DIR}/${ORG_TO_MERGE}_WASHU-MIT_utr${SIDE}_${SEQ_LEN}.wc.gz
MIT_FILE=${ORG_TO_MERGE}_MIT_utr${SIDE}_${SEQ_LEN}.wc.gz		
WASHU_FILE=${ORG_TO_MERGE}_WashU_utr${SIDE}_${SEQ_LEN}.wc.gz 
MERGE_ML=200
REDUNDANT_DIR=${SEQ_DIR}/redundant
discard_redundant_files:
	@echo "Discarding redundant files for organism ${ORG_TO_MERGE}	utr${SIDE}	${SEQ_LEN}"
	mv -f ${SEQ_DIR}/${MIT_FILE} ${REDUNDANT_DIR}
	mv -f ${SEQ_DIR}/${WASHU_FILE} ${REDUNDANT_DIR}

merge_sources_one_org:
	@echo "Merging sources for organism ${ORG_TO_MERGE}	utr${SIDE}	${SEQ_LEN}"
	cat ${REDUNDANT_DIR}/${MIT_FILE} ${REDUNDANT_DIR}/${WASHU_FILE}	\
		| purge-sequence -1str -format wc -mis 0 -ml ${SEQ_LEN}	\
		| convert-seq -from fasta -to wc			\
		| grep -vE '\\\n{${SEQ_LEN}}'				\
		> ${MERGED_FILE}

count_merged:
	wc ${REDUNDANT_DIR}/${MIT_FILE} ${REDUNDANT_DIR}/${WASHU_FILE} ${MERGED_FILE}

################################################################
## Iterate over S_cerevisiae ORFs
ORF_LIST=`cat ${ORF_FILE}`
ORF_DIR=data/orfs
ORF_FILE=${ORF_DIR}/S_cerevisiae_orfs.txt
#ORF_FILE=${ORF_DIR}/to_do.txt
#ORF_FILE=${ORF_DIR}/done.txt
ORF_TASK=ortho_one_orf
iterate_orfs:
	@for orf in ${ORF_LIST}; do				\
		${MAKE} ${ORF_TASK} ORF=$${orf};	\
	done

################################################################
## Get the list of S_cerevisiae ORFs from the locus_tag attributes
orf_list:
	grep -v '^--' ${RSAT}/data/genomes/${REF_ORG}/genome/cds_locus_tag.tab \
		| cut -f 2 | grep -v '^Q' > ${ORF_FILE}

all:
	${MAKE_DOWN} analyze_all_orfs PROXIMITY=stricto UTR_FILES="${STRICTO}"
	${MAKE_UP} analyze_all_orfs PROXIMITY=stricto UTR_FILES="${STRICTO}"
	${MAKE_DOWN} analyze_all_orfs PROXIMITY=lato UTR_FILES="${LATO}"
	${MAKE_UP} analyze_all_orfs PROXIMITY=lato UTR_FILES="${LATO}"
	@${MAKE} index

analyze_all_orfs:
	@${MAKE} open_index
	@${MAKE} iterate_orfs ORF_TASK=analyze_one_orf
	@${MAKE} close_index

#ALL_TASKS= collect_one_orf purge_one_orf align_one_orf tree_one_orf profile_one_orf oligos_one_orf clean_one_orf
ALL_TASKS= collect_one_orf purge_one_orf align_one_orf profile_one_orf index_one_orf
analyze_one_orf: 
	@${MAKE} -k ${ALL_TASKS}

################################################################
## Collect the ortholog sequences for one specific ORF in a separate
## file
ORF=YBL037W

#ORF=YAL003W
#ORF=YBR020W
#ORF=YER150W
#ORF=YAL064W-B
#ORF=YAL067C
#ORF=YJR048W
ORTHO_PREFIX=${ORF}_utr${SIDE}_${SEQ_LEN}_ortho
ORTHO_DIR=${UTR_DIR}/orthologs/${PROXIMITY}
COLLECT_DIR=${ORTHO_DIR}/unaligned
COLLECT_FILE=${COLLECT_DIR}/${ORTHO_PREFIX}.fasta
collect_one_orf:
	@mkdir -p ${COLLECT_DIR}
	zgrep -i ${ORF} ${UTR_FILES}					\
		| perl -pe 's|${UTR_DIR}/||'				\
		| perl -pe 's|_utr${SIDE}_${SEQ_LEN}.wc.gz||'		\
		| perl -pe 's|.*:||'					\
		| convert-seq -from wc -to fasta -o ${COLLECT_FILE}
	@echo "Collected orthologs for ORF ${ORF}	${COLLECT_FILE}"

################################################################
## Purge ortholog sequences to remove almost completely duplicated
## sequences (paralogs + genomes sequences 2 times)
PURGED_DIR=${ORTHO_DIR}/purged
PURGED_FILE=${PURGED_DIR}/${ORTHO_PREFIX}_purged.fasta
PURGE_MIS=3
PURGE_ML=100
purge_one_orf:
	@mkdir -p ${PURGED_DIR}
	purge-sequence -1str -ml ${PURGE_ML} -mis ${PURGE_MIS}	\
		-i ${COLLECT_FILE}				\
		| convert-seq -from fasta -to wc		\
		| grep -vE '\\n+\\'				\
		| convert-seq -from wc -to fasta		\
		> ${PURGED_FILE}
	@echo "Purged orthologs for ORF ${ORF}	${PURGED_FILE}"

################################################################
## Align one ortholog group
ALIGN_DIR=${ORTHO_DIR}/aligned
ALIGN_FILE=${ALIGN_DIR}/${ORTHO_PREFIX}.aln
DND_FILE_ORI=${PURGED_DIR}/${ORTHO_PREFIX}_purged.dnd
DND_FILE=${ALIGN_DIR}/${ORTHO_PREFIX}.dnd
clustal_one_orf: align_one_orf tree_one_orf

## Multiple alignment
align_one_orf:
	@mkdir -p ${ALIGN_DIR}
	clustalw -TYPE=DNA -ALIGN=t		\
		-INFILE=${PURGED_FILE}		\
		-OUTFILE=${ALIGN_FILE}
	@mv ${DND_FILE_ORI} ${DND_FILE}
	@echo "Aligned orthologs for orf ${ORF}	${ALIGN_FILE}"

## Neighbour-Joining tree with bootstrap
TREE_DIR=${ORTHO_DIR}/NJtrees
TREE_FILE_TMP=${ALIGN_DIR}/${ORTHO_PREFIX}.phb
TREE_FILE=${TREE_DIR}/${ORTHO_PREFIX}.phb
tree_one_orf:
	@mkdir -p ${TREE_DIR}
	clustalw -TYPE=DNA -BOOTSTRAP=100	\
		-INFILE=${ALIGN_FILE}
	mv ${TREE_FILE_TMP} ${TREE_DIR}
	@echo "Calculated NJ tree for orf ${ORF}	${TREE_FILE}"

################################################################
## Position-analysis
POSITION_DIR=${UTR_DIR}/genome-wise_oligos/positions
POSITION_PREFIX=${OL}nt_ci${CI}${STR}${NOOV}
POSITION_FILE=${POSITION_DIR}/${ORG}_${POSITION_PREFIX}_positions.tab
POS_ORIGIN=0
CI=20
POSITION_OPTIONS=-v 1 -l ${OL} ${STR} ${NOOV}	\
	-ci ${CI} -nogrouprc -origin ${POS_ORIGIN} -sort	\
	-return chi,distrib,exp,rank		\
	-nofilter				\
	-i ${SUB_UTR_WC_PURGED} -format wc	\
	-o ${POSITION_FILE}

POSITION_CMD=position-analysis ${POSITION_OPTIONS}
positions:
	@${MAKE} iterate_oligo_lengths OLIGO_TASK=positions_one_len

positions_one_len:
	@${MAKE} iterate_organisms ORG_TASK=positions_one_genome_one_len
	@${MAKE} positions_compare

positions_compare:
	compare-scores -sc 4 -o ${POSITION_DIR}/${POSITION_PREFIX}_comparison_chi2.tab -files ${POSITION_DIR}/*_${POSITION_PREFIX}_positions.tab
	compare-scores -sc 3 -o ${POSITION_DIR}/${POSITION_PREFIX}_comparison_occ.tab -files ${POSITION_DIR}/*_${POSITION_PREFIX}_positions.tab

positions_one_genome_one_len:
	@echo ${POSITION_CMD}
	@${POSITION_CMD}

################################################################
## Oligo-analysis
OL=6
STR=-1str
NOOV=-noov
OLIGO_RETURN_FIELDS=occ,freq,mseq,proba,rank,zscore
OLIGO_PREFIX=oligos_${OL}nt${STR}${NOOV}
OLIGO_GENERIC_OPTIONS=				\
	-sort -v ${V}				\
	-return ${OLIGO_RETURN_FIELDS}		\
	-l ${OL} ${STR} ${NOOV}
OLIGO_CMD=oligo-analysis			\
	${OLIGO_GENERIC_OPTIONS}		\
	${OLIGO_OPTIONS}

## ORF-wise oligo analysis
ORF_OLIGO_DIR=${ORTHO_DIR}/ORF-wise_oligos
THOSIG=0
THMSSIG=0
ORF_OLIGO_PREFIX=${ORTHO_PREFIX}_${OLIGO_PREFIX}_osig${THOSIG}_msig${THMSSIG}
ORF_OLIGO_FILE=${ORF_OLIGO_DIR}/${ORF_OLIGO_PREFIX}
ORF_OLIGO_OPTIONS=				\
		 -thosig ${THOSIG}		\
		-thmssig ${THMSSIG}		\
		-org ${REF_ORG} -bg intergenic	\
		-i ${PURGED_FILE}		\
		-o ${ORF_OLIGO_FILE}
oligos_one_orf:
	@mkdir -p ${ORF_OLIGO_DIR}
	@${MAKE} oligos OLIGO_OPTIONS="${ORF_OLIGO_OPTIONS}"

################################################################
## Genome-wise oligo-analysis
oligos_all_genomes:
	@${MAKE} iterate_organisms ORG_TASK=oligos_one_genome NOOV=-ovlp
	@${MAKE} iterate_organisms ORG_TASK=oligos_one_genome NOOV=-noov

compare_oligos_all_genomes:
	@echo "Comparing each pair of genomes"
	@for org1 in ${ORGANISMS}; do							\
		for org2 in ${ORGANISMS}; do						\
			${MAKE} compare_oligos_two_genomes ORG1=$${org1} ORG2=$${org2};	\
		done ;									\
	done

compare_oligos_two_genomes:
	@${MAKE} compare_oligos_two_genomes_one_sc SC=3 NOOV=-noov
	@${MAKE} compare_oligos_two_genomes_one_sc SC=9 NOOV=-noov
	@${MAKE} compare_oligos_two_genomes_one_sc SC=10 NOOV=-noov
	@${MAKE} compare_oligos_two_genomes_one_sc SC=3 NOOV=-ovlp
	@${MAKE} compare_oligos_two_genomes_one_sc SC=9 NOOV=-ovlp

## Compare oligonucleotide composition between genomes
ORG1=S_cerevisiae
ORG2=S_kluyveri
SC=9
OLIGOS_ORG1=${OLIGO_MARKOV_DIR}/${ORG1}_${OLIGO_PREFIX}_mkv${MKV}.gz
OLIGOS_ORG2=${OLIGO_MARKOV_DIR}/${ORG2}_${OLIGO_PREFIX}_mkv${MKV}.gz
GENOME_OLIGO_COMPA_DIR=${OLIGO_MARKOV_DIR}/comparisons
OLIGOS_COMPA_PREFIX=${GENOME_OLIGO_COMPA_DIR}/${ORG1}_vs_${ORG2}_${OLIGO_PREFIX}_mkv${MKV}_sc${SC}
IMG_FORMAT=jpg
compare_oligos_two_genomes_one_sc:
	@echo "Comparing oligos	${ORG1}	${ORG2}	${SC}	${NOOV}"
	@mkdir -p ${GENOME_OLIGO_COMPA_DIR}
	compare-scores -sc ${SC}				\
		-i ${OLIGOS_ORG1} -i ${OLIGOS_ORG2}	\
		-o ${OLIGOS_COMPA_PREFIX}.tab
	XYgraph -xcol 2 -ycol 3				\
		-pointsize 3				\
		-title '${OLIGO_PREFIX}'		\
		-xleg1 '${ORG1}'			\
		-yleg1 '${ORG2}'			\
		-format ${IMG_FORMAT}			\
		-i ${OLIGOS_COMPA_PREFIX}.tab		\
		-o ${OLIGOS_COMPA_PREFIX}.${IMG_FORMAT}	\
		-htmap > ${OLIGOS_COMPA_PREFIX}.html

oligos_one_genome:
	@${MAKE} iterate_oligo_lengths OLIGO_TASK=oligo_freq_one_genome_one_length
	@${MAKE} markov_one_genome

## Count oligonucleotide frequencies in the full set of sequences of one organism
OLIGO_FREQ_DIR=${UTR_DIR}/genome-wise_oligos/frequencies
OLIGO_FREQ_PREFIX=${ORG}_${OLIGO_PREFIX}_freq.tab.gz
OLIGO_FREQ_FILE=${OLIGO_FREQ_DIR}/${OLIGO_FREQ_PREFIX}
OLIGO_LENGTHS= 1 2 3 4 5 6 7 8
OLIGO_FREQ_OPTIONS=				\
	-i ${SUB_UTR_WC_PURGED} -format wc	\
	-o ${OLIGO_FREQ_FILE}
oligo_freq_one_genome_one_length:
	@mkdir -p ${OLIGO_FREQ_DIR}
	@echo
	@echo "Analyzing oligo frequencies for genome	${ORG}	${OLIGO_FREQ_FILE}"
	@${MAKE} oligos OLIGO_RETURN_FIELDS="occ,freq,rank" OLIGO_OPTIONS="${OLIGO_FREQ_OPTIONS}"
	@echo "Done oligos for genome	${ORG}	${GENOME_MARKOV_FILE}"

## Detect over-represented oligos in the full set of sequences for one organism
MKV=3
OLIGO_MARKOV_DIR=${UTR_DIR}/genome-wise_oligos/markov
OLIGO_MARKOV_PREFIX=${ORG}_${OLIGO_PREFIX}_mkv${MKV}.tab.gz
OLIGO_MARKOV_FILE=${OLIGO_MARKOV_DIR}/${OLIGO_MARKOV_PREFIX}
MARKOV_OPTIONS=					\
	-markov ${MKV}				\
	-i ${SUB_UTR_WC_PURGED} -format wc	\
	-o ${OLIGO_MARKOV_FILE}
markov_one_genome:
	@mkdir -p ${OLIGO_MARKOV_DIR}
	@echo
	@echo "Analyzing markov for genome	${ORG}	${OLIGO_MARKOV_FILE}"
	@${MAKE} oligos OLIGO_OPTIONS="${MARKOV_OPTIONS}"
	@echo "Done markov for genome	${ORG}	${OLIGO_MARKOV_FILE}"

oligos:
	@${MAKE} my_command MY_COMMAND="${OLIGO_CMD}"

################################################################
## Convert a clustal file into a profile matrix
PROFILE_DIR=${ORTHO_DIR}/profiles
PROFILE_FILE=${PROFILE_DIR}/${ORTHO_PREFIX}.matrix
profile_one_orf:
	@mkdir -p ${PROFILE_DIR}
	convert-matrix -v ${V} -format clustal -pseudo 0 -decimal 2	\
		-margins						\
		-return profile,parameters				\
	 	-i ${ALIGN_FILE}					\
		-o ${PROFILE_FILE} 
	@echo "Created profile for orf ${ORF}	${PROFILE_FILE}"

################################################################
## Remove the unaligned sequence files for one ORF
clean_one_orf:
	rm -f ${COLLECT_FILE} ${PURGED_FILE}


################################################################
## Addd one ORF to the HTML index
INDEX_FILE=index.html
index: open_index index_all_orfs close_index

open_index:
	@echo "Opening index file ${INDEX_FILE}"
	@echo "<html><body><table border=1 cellpadding=3 cellspacing=3>" > ${INDEX_FILE}

	@echo "<th>Sensu lato<br>5'</th>" >> ${INDEX_FILE}
#	@echo "<th></th>" >> ${INDEX_FILE}

	@echo "<th>Sensu stricto<br>5'</th>" >> ${INDEX_FILE}
#	@echo "<th></th>" >> ${INDEX_FILE}

	@echo "<th>ORF</th>" >> ${INDEX_FILE}
#	@echo "<th></th>" >> ${INDEX_FILE}

	@echo "<th>Sensu stricto<br>3'</th>" >> ${INDEX_FILE}
#	@echo "<th></th>" >> ${INDEX_FILE}

	@echo "<th>Sensu lato<br>3'</th>" >> ${INDEX_FILE}
#	@echo "<th></th>" >> ${INDEX_FILE}

	@echo "<th>Description</th>" >> ${INDEX_FILE}

close_index:
	@echo "Closing index file ${INDEX_FILE}"
	@echo "</table></body></html>" >> ${INDEX_FILE}

index_all_orfs:
	@${MAKE} iterate_orfs ORF_TASK=index_one_orf

ORF_DESCR=`grep -E "[;\[]${ORF};" ${DESCR_FILE}`
description_one_orf:
	@echo ${DESCR_FILE}	${ORF}	${ORF_DESCR}

index_one_orf:
	@echo "Indexing ORF ${ORF}"
	@echo "<tr>" >> ${INDEX_FILE}

	@${MAKE} index_one_orf_one_group PROXIMITY=lato UTR_FILES='${LATO}' SIDE=5 SEQ_LEN=2000 SIDE_NAME=utr5 CELL_BG='#CCCCFF'
	@${MAKE} index_one_orf_one_group PROXIMITY=stricto UTR_FILES='${STRICTO}' SIDE=5 SEQ_LEN=2000 SIDE_NAME=utr5 CELL_BG='#AAAAFF'
	@echo "<th BGCOLOR='#FFDDDD'><tt>${ORF}</tt></th>" >> ${INDEX_FILE}
#	@echo "<td></td>" >> ${INDEX_FILE}
	@${MAKE} index_one_orf_one_group PROXIMITY=stricto UTR_FILES='${STRICTO}' SIDE=3 SEQ_LEN=500 SIDE_NAME=utr3 CELL_BG='#AAFFAA'
	@${MAKE} index_one_orf_one_group PROXIMITY=lato UTR_FILES='${LATO}' SIDE=3 SEQ_LEN=500 SIDE_NAME=utr3 CELL_BG='#CCFFCC'

	@echo "<td WIDTH=400 BGCOLOR='#FFDDDD'>${ORF_DESCR}</td>" >> ${INDEX_FILE}
	@echo "</tr>" >> ${INDEX_FILE}
	@echo "" >> ${INDEX_FILE}

index_one_orf_one_group:
	@echo "<td bgcolor=${CELL_BG}>" >> ${INDEX_FILE}
	@echo "<font size=-1>" >> ${INDEX_FILE}
	${MAKE} link_one_file FILE_TO_LINK=${COLLECT_FILE} LINK_TEXT=raw
	${MAKE} link_one_file FILE_TO_LINK=${PURGED_FILE} LINK_TEXT=purged
	${MAKE} link_one_file FILE_TO_LINK=${ALIGN_FILE} LINK_TEXT=align
	${MAKE} link_one_file FILE_TO_LINK=${TREE_FILE} LINK_TEXT=tree
	${MAKE} link_one_file FILE_TO_LINK=${PROFILE_FILE} LINK_TEXT=profile
	${MAKE} link_one_file FILE_TO_LINK=${ORF_OLIGO_FILE} LINK_TEXT=oligos
	@echo "</font></td>" >> ${INDEX_FILE}
#	@echo "<td></td>" >> ${INDEX_FILE}

link_one_file:
	if [ -f "${FILE_TO_LINK}" ] ; then						\
		echo "[<a href=${FILE_TO_LINK}>${LINK_TEXT}</a>]"  >> ${INDEX_FILE};	\
	fi
	if [ -f "${FILE_TO_LINK}.gz" ] ; then						\
		echo "[<a href=${FILE_TO_LINK}.gz>${LINK_TEXT}.gz</a>]"  >> ${INDEX_FILE};	\
	fi

################################################################
## Synchronize the results to the server
SERVER=rsat.ulb.ac.be
SERVER_DIR=rsa-tools/data/comparative_genomics
BIGRE=jvanheld@${SERVER}:${SERVER_DIR}
TO_SYNC=data
TO_SYNC_DIR=`dirname ${TO_SYNC}`
RSYNC_CMD=${RSYNC} ${OPT} ${TO_SYNC} ${BIGRE}/${TO_SYNC_DIR}/
to_bigre:
	@echo ${RSYNC_CMD}
	@${RSYNC_CMD}


