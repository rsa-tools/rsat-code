################################################################
## tester for the program matrix-scan
include ${RSAT}/makefiles/util.mk
MAKEFILE=${RSAT}/makefiles/matrix-scan_test.mk

################################################################
## List parameters
V=2
SEQ=${RSAT}/public_html/demo_files/Dmelanogaster_eve_up5000.fasta
MATRIX_FILE=${RSAT}/public_html/demo_files/Dmelanogaster_segmentation_12matrices.tf
MATRIX_NAME_OPT=-matrix_name Kr -matrix_name bcd -matrix_name eve
RES_DIR=results/matrix-scan_test
SITES=${RES_DIR}/Dmelanogaster_eve_up5000_sites${SUFFIX}
SUFFIX=
list_param:
	@echo "matrix-scan test parameters"
	@echo "SEQ		${SEQ}"
	@echo "MATRIX_FILE	${MATRIX_FILE}"
	@echo "RES_DIR		${RES_DIR}"
	@echo "SITES		${SITES}.ft"
	@echo "SITES (quick)	${SITES}-quick.ft"
	@echo "SUFFIX		${SUFFIX}"
	@echo "MATRIX_NAME_OPT	${MATRIX_NAME_OPT}"
	@echo "OPT		${OPT}"
	@echo "MATRIX_SCAN_CMD	${MATRIX_SCAN_CMD}"

################################################################
## Scan Drosophila even-skipped promoter with matrices corresponding
## to 12 transcription factors involved in embryonic segmentation.
MATRIX_SCAN_CMD=matrix-scan -v ${V} \
		-matrix_format transfac -m ${MATRIX_FILE} \
		${MATRIX_NAME_OPT} \
		-pseudo 1 -decimals 1 -2str -bginput -markov 0 -bg_pseudo 0.01 \
		-seq_format fasta -n score -origin genomic \
		-return limits,sites,pval -lth score 1 -uth pval 1e-4 \
		-i ${SEQ} \
		${OPT} -o ${SITES}.ft
matrix_parameters:
	@echo
	@echo "Matrix parameters"
	convert-matrix -v ${V} -matrix_format transfac -i ${MATRIX_FILE} -from transfac -to tab -pseudo 1 -decimals 2 -return counts,frequencies,margins,parameters

scan_eve:
	@echo
	@echo "SEQ		${SEQ}"
	@echo "MATRIX_FILE	${MATRIX_FILE}"
	@echo "${DATE}	Scanning"
	@mkdir -p ${RES_DIR}
	${MATRIX_SCAN_CMD}
	@echo "${DATE}	Result file"
	@echo "	${SITES}.ft"

## Select a subset of matrices specified by their name
JASPAR=${RSAT}/public_html/motif_databases/JASPAR/jaspar_core_insects_2013-11.tf
scan_selected_names:
	${MAKE} scan_eve OPT='-matrix_name hb,eve' SUFFIX=_selected_names_hb_eve

## Select a subset of JASPAR matrices specified by their accession
JASPAR=${RSAT}/public_html/motif_databases/JASPAR/jaspar_core_insects_2013-11.tf
scan_jaspar_selected_acs:
	${MAKE} scan_eve MATRIX_FILE=${JASPAR}  OPT='-matrix_ac MA0049.1,MA0221.1' SUFFIX=_jaspar_selected_ids_hb_eve

## Select a subset of JASPAR matrices specified by their name
JASPAR=${RSAT}/public_html/motif_databases/JASPAR/jaspar_core_insects_2013-11.tf
scan_jaspar_selected_names:
	${MAKE} scan_eve MATRIX_FILE=${JASPAR}  OPT='-matrix_name hb,eve' SUFFIX=_jaspar_selected_names_hb_eve

################################################################
## Scan even-skipped promoter with the option -quick
scan_eve_quick:
	${MAKE} scan_eve OPT=-quick SUFFIX=-quick 

