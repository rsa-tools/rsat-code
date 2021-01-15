################################################################
## Tester for merge-matrices

include ${RSAT}/makefiles/util.mk
MAKEFILE=${RSAT}/makefiles/merge-matrices_demo.mk

## Verbosity
V=2

# ## Study case 1
# PEAK_SET=Sox2
# PEAK_SEQ=${RSAT}/public_html/demo_files/${PEAK_SET}_peaks.fasta
# SEQ_NAME=${PEAK_SET}
# SEQ_FILE=${PEAK_SEQ}

## Study case 2
SEQ_PREFIX=MET_up800-noorf
SEQ_FILE=${RSAT}/public_html/demo_files/${SEQ_PREFIX}.fasta
IN_MATRIX_PREFIX=CACGTn
OUT_MATRIX_PREFIX=CACGTk
KMER_MATRICES=${RSAT}/public_html/demo_files/matrices_from_kmers/${IN_MATRIX_PREFIX}.tf

## Run merge-matrices
CALC_MODE=sum
MERGED_DIR=results/merge-matrices_test/${SEQ_NAME}
FINAL_MATRICES=${MERGED_DIR}/${OUT_MATRIX_PREFIX}_${CALC_MODE}_merged.tf
ALL_MATRICES=`echo; ls -1 ${MERGED_DIR}/${OUT_MATRIX_PREFIX}_${CALC_MODE}_*.tf`
ALL_FEATURES=`echo; ls -1 ${MERGED_DIR}/${OUT_MATRIX_PREFIX}_${CALC_MODE}_*.ft`
ALL_SEQ=`echo; ls -1 ${MERGED_DIR}/${OUT_MATRIX_PREFIX}_${CALC_MODE}_*.fasta`

## Iterations for rescan-matrix
ITERATIONS=3
SCAN_PARAM=-uth pval 0.001

param:
	@echo "Parameters"
	@echo "	SEQ_PREFIX		${SEQ_PREFIX}"
	@echo "	SEQ_FILE		${SEQ_FILE}"
	@echo "	IN_MATRIX_PREFIX	${IN_MATRIX_PREFIX}"
	@echo "	OUT_MATRIX_PREFIX	${OUT_MATRIX_PREFIX}"
	@echo "	KMER_MATRICES		${KMER_MATRICES}"
	@echo "	CALC_MODE		${CALC_MODE}"
	@echo "	MERGED_DIR		${MERGED_DIR}"
	@echo "	ITERATIONS		${ITERATIONS}"
	@echo "	SCAN_PARAM		${SCAN_PARAM}"
	@echo "	OPT			${OPT}"
	@echo "	MERGE_COMMAND		${MERGE_COMMAND}"
	@echo "	FINAL_MATRICES		${FINAL_MATRICES}"
	@echo "	ALL_MATRICES		${ALL_MATRICES}"
	@echo "	ALL_FEATURES		${ALL_FEATURES}"
	@echo "	ALL_SEQ			${ALL_SEQ}"


MERGE_COMMAND=merge-matrices -v ${V}	\
	-i ${KMER_MATRICES}		\
	-in_format tf			\
	-out_format tf			\
	-calc ${CALC_MODE}		\
	-id ${OUT_MATRIX_PREFIX}	\
	-name ${OUT_MATRIX_PREFIX}	\
	-seq ${SEQ_FILE}		\
	-iterations ${ITERATIONS}	\
	-scan_param '${SCAN_PARAM}' 	\
	${OPT} -o ${FINAL_MATRICES}

merge_kmers:
	@echo "Merging matrices"
	@mkdir -p ${MERGED_DIR}
	${MERGE_COMMAND}
	@echo "	SEQ_PREFIX		${SEQ_PREFIX}"
	@echo "	SEQ_FILE		${SEQ_FILE}"
	@echo "	IN_MATRIX_PREFIX	${IN_MATRIX_PREFIX}"
	@echo "	KMER_MATRICES		${KMER_MATRICES}"
	@echo "	MERGED_DIR		${MERGED_DIR}"
	@echo "	OUT_MATRIX_PREFIX	${OUT_MATRIX_PREFIX}"
	@echo "	FINAL_MATRICES		${FINAL_MATRICES}"
	@echo "	ALL_MATRICES		${ALL_MATRICES}"
	@echo "	ALL_FEATURES		${ALL_FEATURES}"
	@echo "	ALL_SEQ			${ALL_SEQ}"
