################################################################
## Quick test to compare the native matrix-scan with the quick version


include ${RSAT}/makefiles/util.mk
MAKEFILE=${RSAT}/makefiles/test_scan_vs_quick.mk

V=2
SEQ=${RSAT}/public_html/demo_files/matrix_scan_demo_sequences.fasta
MATRICES=${RSAT}/public_html/demo_files/matrix_scan_demo_matrices.tf
QUICK=-quick
RES_DIR=results
SITES=${RES_DIR}/Eve12${QUICK}_sites.ft
BG_FILE=${RSAT}/data/genomes/Drosophila_melanogaster/oligo-frequencies/1nt_upstream-noorf_Drosophila_melanogaster-ovlp-1str.freq
ORIGIN=genomic
scan:
	@mkdir -p ${RES_DIR}
	matrix-scan ${QUICK} -v ${V} \
		-i ${SEQ} -seq_format fasta \
		-m ${MATRICES} -matrix_format tf \
		-pseudo 1 -decimals 1 -2str \
		-origin ${ORIGIN}\
		-bgfile ${BG_FILE} \
		-bg_pseudo 0.01 -n score \
		-return sites -uth pval 1e-4 \
		-o ${SITES}
	@echo ${SITES}

quick:
	${MAKE} scan QUICK=-quick

slow:
	${MAKE} scan QUICK=''

all: quick slow


