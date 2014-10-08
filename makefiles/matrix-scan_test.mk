################################################################
## tester for the program matrix-scan
include ${RSAT}/makefiles/util.mk
MAKEFILE=${RSAT}/makefiles/matrix-scan_test.mk

################################################################
## Scan Drosophila even-skipped promoter with matrices corresponding
## to 12 transcription factors involved in embryonic segmentation.
V=2
SEQ=${RSAT}/public_html/demo_files/Dmelanogaster_eve_up5000.fasta
MATRIX_FILE=${RSAT}/public_html/demo_files/Dmelanogaster_segmentation_12matrices.tf
RES_DIR=results/matrix-scan_test
SITES=${RES_DIR}/Dmelanogaster_eve_up5000_sites${SUFFIX}.ft
SUFFIX=
scan_eve:
	@echo
	@echo "SEQ		${SEQ}"
	@echo "MATRIX_FILE	${MATRIX_FILE}"
	@echo "Scanning"
	@mkdir -p ${RES_DIR}
	matrix-scan -v ${V} \
		-matrix_format transfac -m ${MATRIX_FILE} \
		-pseudo 1 -decimals 1 -2str -bginput -markov 0 -bg_pseudo 0.01 \
		-seq_format fasta -n score -origin genomic \
		-return limits,sites,pval -lth score 1 -uth pval 1e-4 \
		-i ${SEQ} \
		${OPT} -o ${SITES}
	@echo "Result file"
	@echo "	${SITES}"

## Select a subset of matrices specified by their name
JASPAR=${RSAT}/public_html/data/motif_databases/JASPAR/jaspar_core_insects_2013-11.tf
scan_selected_names:
	${MAKE} scan_eve OPT='-matrix_name hb,eve' SUFFIX=_selected_names_hb_eve

## Select a subset of JASPAR matrices specified by their accession
JASPAR=${RSAT}/public_html/data/motif_databases/JASPAR/jaspar_core_insects_2013-11.tf
scan_jaspar_selected_acs:
	${MAKE} scan_eve MATRIX_FILE=${JASPAR}  OPT='-matrix_ac MA0049.1,MA0221.1' SUFFIX=_jaspar_selected_ids_hb_eve

## Select a subset of JASPAR matrices specified by their name
JASPAR=${RSAT}/public_html/data/motif_databases/JASPAR/jaspar_core_insects_2013-11.tf
scan_jaspar_selected_names:
	${MAKE} scan_eve MATRIX_FILE=${JASPAR}  OPT='-matrix_name hb,eve' SUFFIX=_jaspar_selected_names_hb_eve

################################################################
## Scan even-skipped promoter with the option -quick
scan_eve_quick:
	${MAKE} scan_eve OPT=-quick SUFFIX=-quick
