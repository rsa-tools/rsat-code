################################################################
## Demonstration for the tool footprint-scan

include ${RSAT}/makefiles/util.mk
MAKEFILE=${RSAT}/makefiles/footprint-scan_demo.mk

LEXA_MATRIX=${RSAT}/public_html/demo_files/LexA.2nt_upstream-noorf-ovlp-2str.20.tf
ORG=Escherichia_coli_K_12_substr__MG1655_uid57779
TAXON=Enterobacteriaceae
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
		-org ${ORG} -taxon Enterobacteriaceae \
		-m ${LEXA_MATRIX} -matrix_format transfac \
		${QUERY} \
		-info_lines \
		-bgfile ${RSAT}/public_html/data/genomes/${ORG}/oligo-frequencies/1nt_upstream-noorf_${ORG}-ovlp-2str.freq.gz \
		-plot_format jpg -map_format jpg \
		${BATCH} ${OPT} \
		-o ${FP_SCAN_DIR}


################################################################
## Run footprint-scan with selected gene
fp_scan_some_genes: list_param
	@${MAKE} _fp_scan QUERY='-q lexA -q recA -q uvrA -q uvrB'

################################################################
## Run footprint scan with all the genes of the genome of
## interest. This costs several hours of computation, it should better
## run on a cluster with the option -batch (see target
## fp_scan_all_genes_batch) or with a resricted number of genes
## (target fp_scan_gene_slice).
fp_scan_all_genes: list_param
	@${MAKE} _fp_scan  QUERY='-all_genes'

## Run the analysis of each gene of a genome in batch. This requires a
## PC cluster and a properly configured job manager (see cluster
## options in RSAT_config.props).
fp_scan_all_genes_batch: list_param
	@${MAKE} fp_scan_all_genes BATCH='-batch'

## Run footprint scan for a limited number of genes amon all the
## genes of the genome (for quick tests and debugging)
SKIP=0
LAST=25
fp_scan_gene_slice: list_param
	@${MAKE} fp_scan_all_genes OPT='-skip ${SKIP} -last ${LAST}'

