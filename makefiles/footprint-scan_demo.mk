################################################################
## Demonstration for the too footprint-scan

include ${RSAT}/makefiles/util.mk
MAKEFILE=${RSAT}/makefiles/footprint-scan_demo.mk

LEXA_MATRIX=${RSAT}/public_html/demo_files/LexA.2nt_upstream-noorf-ovlp-2str.20.tf
ORG=Escherichia_coli_K_12_substr__MG1655_uid57779
TAXON=Enterobacteriales
GENE=lexA
QUERY=-q ${GENE}
BATCH=
TASK=query_seq,filter_dyads,orthologs,ortho_seq,purge,dyads,maps,gene_index,index

list_param:
	@echo "Footprint-scan parameters"
	@echo "	LEXA_MATRIX	${LEXA_MATRIX}"
	@echo "	ORG		${ORG}"
	@echo "	TAXON		${TAXON}"
	@echo "	GENE		${GENE}"
	@echo "	QUERY		${QUERY}"
	@echo "	TASK		${TASK}"
	@echo "	SKIP		${SKIP}"
	@echo "	LAST		${LAST}"
	@echo "	BATCH		${BATCH}"
	@echo "	FP_SCAN_DIR	${FP_SCAN_DIR}"

## Generic command for footprint-scan (will be adapted in other
## targets by changing parameters).
FP_SCAN_DIR=results/footprint-scan_demo
_fp_scan:
	@mkdir -p ${FP_SCAN_DIR}
	@echo
	@echo "Running footprint-scan 	${ORG}	${TAXON}	${QUERY}"
	footprint-scan -v ${V} -nodie -synthesis -sep_genes \
		-org ${ORG} -taxon Enterobacteriales \
		-m ${LEXA_MATRIX} -matrix_format transfac \
		${QUERY} \
		-info_lines \
		-bgfile ${RSAT}/public_html/data/genomes/${ORG}/oligo-frequencies/1nt_upstream-noorf_${ORG}-ovlp-2str.freq.gz \
		-plot_format jpg -map_format jpg \
		-o ${FP_SCAN_DIR}


################################################################
## Run footprint-scan with selected gene
fp_scan_some_genes: list_param
	@${MAKE} _fp_scan QUERY='-q lexA -q recA -q uvrA -q uvrB'
