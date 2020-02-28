################################################################
## Tester for rescan-matrix

include ${RSAT}/makefiles/util.mk
MAKEFILE=${RSAT}/makefiles/rescan-matrix_demo.mk

## Verbosity
V=2

# ## Study case 1
# PEAK_SET=Sox2
# PEAK_SEQ=${RSAT}/public_html/demo_files/${PEAK_SET}_peaks.fasta
# SEQ_NAME=${PEAK_SET}
# SEQ_FILE=${PEAK_SEQ}

## Study case 2
METSET=yeast_MET_promoters
MET_SEQ=${RSAT}/public_html/demo_files/MET_up800-noorf.fasta
MET_MATRIX_PREFIX=CACGTk
MET_MATRICES=${RSAT}/public_html/demo_files/matrices_from_kmers/CACGTk.tf
SEQ_NAME=${METSET}
SEQ_FILE=${MET_SEQ}
MATRICES=${MET_MATRICES}
MATRIX_PREFIX=${MET_MATRIX_PREFIX}

ITERATIONS=3

## Riun rescan-matrix
RESCAN_COMMAND=rescan-matrix -v ${V}	\
	-seq ${SEQ_FILE}		\
	-m ${MATRICES}			\
	-matrix_format tf		\
	-iterations ${ITERATIONS} 	\
	-o ${RESCAND_MATRICES}

RESCAND_DIR=results/merge-matrices_test/${SEQ_NAME}
RESCAND_MATRICES=${RESCAND_DIR}/${MATRIX_PREFIX}_rescanned-from-${SEQ_NAME}.tf
_rescan:
	@echo "Merging matrices"
	@mkdir -p ${RESCAND_DIR}
	${RESCAN_COMMAND}
	@echo "	SEQ_FILE	${SEQ_FILE}"
	@echo "	MATRICES	${MATRICES}"
	@echo "	RESCAND_DIR	${RESCAND_DIR}"
	@echo "	RESCAND_MATRICES	${RESCAND_MATRICES}"
	@echo "	SITES		${SITES}"

rescan_kmers:
	${MAKE} _rescan SEQ_NAME=${METSET} SEQ_FILE=${MET_SEQ} MATRICES=${MET_MATRICES} MATRIX_PREFIX=${MET_MATRIX_PREFIX}

