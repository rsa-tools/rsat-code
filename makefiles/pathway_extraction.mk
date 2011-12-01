################################################################
## test for the pathway extraction tool

include ${RSAT}/makefiles/util.mk
MAKEFILE=${RSAT}/makefiles/pathway_extraction.mk

ORG=Escherichia_coli_K12

QUERY=met

## Extract all E.coli genes whose name starts with ${QUERY}
genes:
	gene-info -org ${ORG} -feattype CDS -full -q '^${QUERY}.*' -o ${QUERY}_genes.tab
	@echo ${QUERY}_genes.tab

## Select the first column, containing gene Ids.
gene_ids:
	grep -v "^;" ${QUERY}_genes.tab | cut -f 1 > ${QUERY}_genes_IDs.txt
	@echo ${QUERY}_genes_IDs.txt

## Extract a pathway connecting at best the reactions catalyzed by these gene products
METAB_DB=MetaCyc
METAB_NETWORK=MetaCyc_directed_141
NETWORK=${RSAT}//public_html/data/metabolic_networks/networks/${METAB_DB}/${METAB_NETWORK}.txt
#OUTDIR=results/${ORG}/${QUERY}
OUTDIR=${OUTDIR}/${QUERY}_${ORG}_${METAB_NETWORK}_pred_pathways.dot
extract_pathway:
	pathway-extractor -i ${QUERY}_genes_IDs.txt \
		-g ${NETWORK} \
		-ger ${RSAT}/data/metabolic_networks/GER_files/GPR_Uniprot_112011_${ORG}.tab \
		-o ${OUTDIR} \
		-t temp_dir


VIEWER=dotty
view_pathway:
	${VIEWER}