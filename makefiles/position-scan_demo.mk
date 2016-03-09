################################################################
## Demo for the RSAT tool matrix-clustering
##
## Authors: Jaime Castro-Mondragon
## Date: 2016


include ${RSAT}/makefiles/util.mk
MAKEFILE=${RSAT}/makefiles/position-scan_demo.mk

################################################################
## Parameters for the analysis
MATRIX_FORMAT=tf
SEQ_FORMAT=fasta
MKV_ORDER=1
V=2
BIN=25
PVAL=1e-3
MARKOV_ORDER=1
OPT=

## Define a set of demo files
JUND_PREFIX=Jun_Chip_seq

## Choose a particular demo set
PROFILE_PREFIX=${JUND_PREFIX}

## Define file locations based on the chosen demo set
INPUT_FILES_DIR=${RSAT}/public_html/demo_files
MATRIX_FILE=${INPUT_FILES_DIR}/${PROFILE_PREFIX}_matrices.tf
SEQUENCE_FILE=${INPUT_FILES_DIR}/${PROFILE_PREFIX}_sequences.fasta

list_param:
	@echo "PROFILE_PREFIX		${PROFILE_PREFIX}"
	@echo "INPUT_FILES_DIR		${INPUT_FILES_DIR}"
	@echo "MATRIX_FILE		${MATRIX_FILE}"
	@echo "SEQUENCE_FILE		${SEQUENCE_FILE}"
#	@echo "		${}"




################################################################
## Run position-scan on one demo set (the particular cases will be
## specified below)
TITLE='position-scan result'
POSITION_PROFILE_BASENAME=${PROFILE_PREFIX}_bin_size_${BIN}_pval${PVAL}_mkv_${MARKOV_ORDER}
POSITION_PROFILE_DIR=results/position_profile_results/${PROFILE_PREFIX}/${BIN}_nt_bin/pval${PVAL}/mkv_${MARKOV_ORDER}
POSITION_PROFILE_FILE_PREFIX=${POSITION_PROFILE_DIR}/${POSITION_PROFILE_BASENAME}
PROFILE_CMD=position-scan -v ${V} \
		-matrix ${MATRIX_FILE} \
                -matrix_format ${MATRIX_FORMAT} \
		-title '${TITLE}' \
		-seq ${SEQUENCE_FILE} \
                -seq_format ${SEQ_FORMAT} \
                -bginput \
                -markov ${MARKOV_ORDER} \
                -bin ${BIN} \
		-pval ${PVAL} ${OPT}\
		-o ${POSITION_PROFILE_FILE_PREFIX}

_profiles:
	@echo
	@echo "Running position-scan	${PROFILE_PREFIX}	${OPT}"
	${MAKE} my_command MY_COMMAND="${PROFILE_CMD}"
	@echo "		${PROFILE_CMD}"
	@echo "		${POSITION_PROFILE_FILE_PREFIX}_report.html"


## Cluster motifs resulting from 12 independent analysis of peak-motifs (Chen data set). 
HOCOMOCO_MATRICES=${RSAT}/public_html/motif_databases/HOCOMOCO/HOCOMOCO_NonRedundant_2015-11-23_Human_Ncor0.8_cor0.85.tf
profile_Jun_ChIPseq_peaks:
	@echo
	@echo "Running position-scan with Hocomoco Human motifs on ChIP-seq peaks of JUND"
	${MAKE} _profiles PROFILE_PREFIX=${JUND_PREFIX}  MATRIX_FILE=${HOCOMOCO_MATRICES}\
		TITLE='Hocomoco motifs in JunD peaks'
