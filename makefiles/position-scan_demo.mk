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
RANDOM_SEQ=
RANDOM_MOTIFS=

## Define a set of demo files
JUND_PREFIX=Jun_Chip_seq
JUND_RANDOM_SEQ_PREFIX=Jun_Chip_seq_random_sequences


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
		-pval ${PVAL} \
                ${OPT} \
                ${RANDOM_SEQ} \
                ${RANDOM_MOTIFS} \
		-o ${POSITION_PROFILE_FILE_PREFIX}

_profiles:
	@echo
	@echo "Running position-scan	${PROFILE_PREFIX}	${OPT}"
	${MAKE} my_command MY_COMMAND="${PROFILE_CMD}"
	@echo "		${PROFILE_CMD}"
	@echo "		${POSITION_PROFILE_FILE_PREFIX}_report.html"


################################################################
## Profiles of a set of PSSMs in a set of JUND ChIP-seq peaks
## taken from Encode project
HOCOMOCO_MATRICES=${RSAT}/public_html/demo_files/HOCOMOCO_Human_Ncor0.8_cor0.85_JUND_Demo.tf
profile_Jun_ChIPseq_peaks:
	@echo
	@echo "Running position-scan with Hocomoco Human motifs on ChIP-seq peaks of JUND"
	${MAKE} _profiles PROFILE_PREFIX=${JUND_PREFIX}  MATRIX_FILE=${HOCOMOCO_MATRICES} BIN=31\
                TITLE='Hocomoco motifs in JunD peaks'

################################################################
## Profiles of a set of PSSMs in a set of JUND ChIP-seq peaks
## taken from Encode project
## NEGATIVE CONTROL: the original sequences are randomly permuted
HOCOMOCO_MATRICES=${RSAT}/public_html/demo_files/HOCOMOCO_Human_Ncor0.8_cor0.85_JUND_Demo.tf
profile_Jun_ChIPseq_peaks_random_sequences:
	@echo
	@echo "Running position-scan with Hocomoco Human motifs on ChIP-seq peaks of JUND"
	@echo "Negative Control: suffled sequences. The nucleotide interdependence is lost"
	${MAKE} _profiles PROFILE_PREFIX=${JUND_PREFIX} RANDOM_MOTIFS='-rand_seq'\
               POSITION_PROFILE_BASENAME='Jun_Chip_seq_random_sequences_bin_size_${BIN}_pval${PVAL}_mkv_${MARKOV_ORDER}'\
               MATRIX_FILE=${HOCOMOCO_MATRICES}\
               POSITION_PROFILE_DIR=${POSITION_PROFILE_DIR}'/RANDOM_SEQ'\
	       TITLE='Hocomoco motifs in JunD peaks; Random Sequences'
