################################################################
## Demonstration for the tool footprint-discovery

include ${RSAT}/makefiles/util.mk
MAKEFILE=${RSAT}/makefiles/footprint-discovery_demo.mk

#ORG=Escherichia_coli_K_12_substr__MG1655_uid57779
V=1
ORG=Escherichia_coli_GCF_000005845.2_ASM584v2
TAXON=Gammaproteobacteria
GENE=lexA
QUERY=-q ${GENE}
BATCH=
DIAMOND=-diamond
UNIQUE_OPT=
#UNIQUE_OPT=-unique_species
FP_DISCO_DIR=results/footprints${INFER_OPERONS}/${BG_MODEL}
BG_MODEL=taxfreq
SKIP=0
LAST=25
#INFER_OPERONS=-infer_operons
INFER_OPERONS=
#TASK=all
TASK=query_seq,filter_dyads,orthologs,ortho_seq,purge,dyads,maps,gene_index,index
#TASK=operons,query_seq,filter_dyads,orthologs,ortho_seq,purge,dyads,maps,gene_index,index
#TASK==bg_model,operons,query_seq,filter_dyads,orthologs,ortho_seq,purge,dyads,map_dyads,map_pssm,gene_index,index
#TASK==bg_model,operons,query_seq,filter_dyads,orthologs,ortho_seq,purge,dyads,map_dyads,map_pssm,network,gene_index,network_index,index

list_param:
	@echo "Parameters"
	@echo "	ORG		${ORG}"
	@echo "	TAXON		${TAXON}"
	@echo "	GENE		${GENE}"
	@echo "	GENE_FILE	${GENE_FILE}"
	@echo "	QUERY		${QUERY}"
	@echo "	TASK		${TASK}"
	@echo "	DIAMOND		${DIAMOND}"
	@echo "	BG_MODEL	${BG_MODEL}"
	@echo "	SKIP		${SKIP}"
	@echo "	LAST		${LAST}"
	@echo "	BATCH		${BATCH}"
	@echo "	FP_DISCO_DIR	${FP_DISCO_DIR}"
	@echo "	UNIQUE_OPT	${UNIQUE_OPT}"
	@echo "	INFER_OPERONS	${INFER_OPERONS}"
	@echo "	OPT		${OPT}"
	@echo "	FP_CMD		${FP_CMD}"

## Generic command for footprint-discovery (will be adapted in other
## targets by changing parameters).
FP_CMD=footprint-discovery -v ${V} ${DIAMOND} \
	-org ${ORG} \
	-taxon ${TAXON} \
	${QUERY} ${INFER_OPERONS} \
	${UNIQUE_OPT} \
	-sep_genes \
	-lth occ 1 \
	-lth occ_sig 0 \
	-uth rank 50 \
	-return occ,proba,rank \
	-filter \
	-bg_model ${BG_MODEL} \
	-task ${TASK} ${BATCH} ${OPT} \
	-o ${FP_DISCO_DIR}
_fp_disco:
	@mkdir -p ${FP_DISCO_DIR}
	@echo
	@echo "Running footprint-discovery 	${ORG}	${TAXON}	${QUERY}"
	@${FP_CMD}

################################################################
## Run footprint-discovery with one selected gene
fp_disco_one_gene: list_param
	@${MAKE} _fp_disco QUERY='-q ${GENE}'


################################################################
## Run footprint discovery with selected genes
##
## THIS FAILS HERE:  &InitQueryOutput(), footprint.lib.pl line389
SELECTED_GEENES=${YEAST_DEMO_GENES}
fp_disco_selected_genes:
	@echo "Selected ${ORG} genes	${GENE_FILE}"
	@${MAKE} _fp_disco QUERY='-genes ${GENE_FILE}' 

YEAST_DEMO_GENES=${RSAT}/public_html/demo_files/footprint-discovery_selected_yeast_genes.tab
fp_disco_selected_genes_yeast:
	@${MAKE} fp_disco_selected_genes \
		GENE_FILE=${YEAST_DEMO_GENES} \
		ORG=Saccharomyces_cerevisiae \
		TAXON=Saccharomycetales

################################################################
## Run footprint discovery with all the genes of the genome of
## interest. This costs several hours of computation, it should better
## run on a cluster with the option -batch (see target
## fp_disco_all_genes_batch) or with a resricted number of genes
## (target fp_disco_gene_slice).
fp_disco_all_genes: list_param
	@${MAKE} _fp_disco  QUERY='-all_genes'

################################################################
## (Re)generae summary index
index_all_genes: list_param
	@${MAKE} _fp_disco  QUERY='-all_genes' TASK=index,gene_index

## Run the analysis of each gene of a genome in batch. This requires a
## PC cluster and a properly configured job manager (see cluster
## options in RSAT_config.props).
fp_disco_all_genes_batch: list_param
	@${MAKE} fp_disco_all_genes BATCH='-batch'

## Run footprint discovery for a limited number of genes among all the
## genes of the genome (for quick tests and debugging)
fp_disco_gene_slice: list_param
	@${MAKE} fp_disco_all_genes OPT='-skip ${SKIP} -last ${LAST}'

