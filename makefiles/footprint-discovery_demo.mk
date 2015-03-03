################################################################
## Demonstration for the too footprint-discovery

include ${RSAT}/makefiles/util.mk
MAKEFILE=${RSAT}/makefiles/footprint-discovery_demo.mk

ORG=Escherichia_coli_K_12_substr__MG1655_uid57779
TAXON=Enterobacteriales
GENE=lexA
QUERY=-q ${GENE}
BATCH=
TASK=query_seq,filter_dyads,orthologs,ortho_seq,purge,dyads,maps,gene_index,index

list_param:
	@echo "Parameters"
	@echo "ORG		${ORG}"
	@echo "TAXON		${TAXON}"
	@echo "GENE		${GENE}"
	@echo "QUERY		${QUERY}"
	@echo "TASK		${TASK}"
	@echo "SKIP		${SKIP}"
	@echo "LAST		${LAST}"
	@echo "BATCH		${BATCH}"
	@echo "FP_DISCO_DIR	${FP_DISCO_DIR}"

## Run footprint-discovery for a single gene of interest.
FP_DISCO_DIR=results/footprint-discovery_demo
_fp_disco:
	@mkdir -p ${FP_DISCO_DIR}
	@echo
	@echo "Running footprint-discovery 	${GEBE}	${ORG}	${TAXON}"
	footprint-discovery  -v 1 -org ${ORG} -taxon ${TAXON} \
		${QUERY} \
		-lth occ 1 \
		-lth occ_sig 0 \
		-uth rank 50 \
		-return occ,proba,rank \
		-filter \
		-bg_model taxfreq \
		-task ${TASK} ${BATCH} ${OPT} \
		-o ${FP_DISCO_DIR}

fp_disco_one_gene: list_param
	@${MAKE} _fp_disco QUERY='-q ${GENE}'

fp_disco_all_genes: list_param
	@${MAKE} _fp_disco  QUERY='-all_genes'

SKIP=0
LAST=25
fp_disco_gene_slice: list_param
	@${MAKE} fp_disco_all_genes OPT='-skip ${SKIP} -last ${LAST}'

fp_disco_all_genes_batch: list_param
	@${MAKE} fp_disco_all_genes BATCH='-batch'
