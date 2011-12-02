################################################################
## test for the pathway extraction tool

include ${RSAT}/makefiles/util.mk
MAKEFILE=${RSAT}/makefiles/pathway_extraction.mk

ORG=Escherichia_coli_K12

QUERY=met

################################################################
## Pathway extraction from a single group of genes
##
## As study case, we extract all E.coli genes whose name starts with
## 'met' (the methionine genes) and extract a sub-graph of the generic
## metabolic network (MetaCyc).

## Run the successive steps
one_group: genes gene_ids one_pathway

## Extract all E.coli genes whose name starts with ${QUERY}
genes:
	@echo "Getting genes for query ${QUERY}.*"
	gene-info -org ${ORG} -feattype CDS -full -q '^${QUERY}.*' -o ${QUERY}_genes.tab
	@echo
	@echo ${QUERY}_genes.tab

## Select the first column, containing gene Ids.
gene_ids:
	@echo "Selecting gene IDs"
	grep -v "^;" ${QUERY}_genes.tab | cut -f 1 > ${QUERY}_genes_IDs.txt
	@echo
	@echo `wc ${QUERY}_genes_IDs.txt` IDs
	@echo ${QUERY}_genes_IDs.txt

## Extract a pathway connecting at best the reactions catalyzed by these gene products
METAB_DB=MetaCyc
METAB_NETWORK=MetaCyc_directed_141
GRAPH_FILE=${RSAT}//public_html/data/metabolic_networks/networks/${METAB_DB}/${METAB_NETWORK}.txt
OUTDIR=results/${ORG}/${QUERY}
OUT_PREFIX=${OUTDIR}/${QUERY}_genes_IDs_Escherichia_coli_strain_K12_-83333_${METAB_NETWORK}_annot_pred_pathways
one_pathway:
	@echo "Extracting pathways from genes ${QUERY}_genes_IDs.txt"
	pathway-extractor -i ${QUERY}_genes_IDs.txt \
		-g ${GRAPH_FILE} \
		-ger ${RSAT}/data/metabolic_networks/GER_files/GPR_Uniprot_112011_${ORG}.tab \
		-o ${OUTDIR} \
		-t temp_dir
	@echo
	@echo ${OUT_PREFIX}*

## View the dot-formatted network
VIEWER_DOT=dotty
view_dot:
	${VIEWER_DOT} ${OUT_PREFIX}.dot

## View the  png-formatted image
VIEWER_PNG=open -a Preview
view_png:
	${VIEWER_PNG} ${OUT_PREFIX}.png
