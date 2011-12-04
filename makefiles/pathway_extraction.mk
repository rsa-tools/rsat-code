################################################################
## test for the pathway extraction tool
include ${RSAT}/makefiles/util.mk
MAKEFILE=${RSAT}/makefiles/pathway_extraction.mk

ORG=Escherichia_coli_K12
QUERY=met
RES_DIR=results/${ORG}/${QUERY}

V=1

################################################################
## Get network data from the main RSAT server
get_data:
	@echo "Synchronizing data from the server to local RSAT"
	rsync -ruptvl -e ssh -z rsat@rsat.ulb.ac.be:rsa-tools/public_html/data/metabolic_networks ${RSAT}/public_html/data/
	@echo "Data synchronized"

################################################################
## Pathway extraction from a single group of genes
##
## As study case, we extract all E.coli genes whose name starts with
## 'met' (the methionine genes) and extract a sub-graph of the generic
## metabolic network (MetaCyc).

## Run the successive steps
one_query: genes gene_ids one_pathway

## Extract all E.coli genes whose name starts with ${QUERY}
GENES=${RES_DIR}/${QUERY}_genes.tab
genes:
	@mkdir -p ${RES_DIR}
	@echo "Getting genes for query ${QUERY}.*"
	gene-info -org ${ORG} -feattype CDS -full -q '^${QUERY}.*' -o ${GENES}
	@echo
	@echo ${GENES}

## Select the first column, containing gene Ids.
IDS=${RES_DIR}/${QUERY}_IDs.tab
gene_ids:
	@echo "Selecting gene IDs"
	grep -v "^;" ${QUERY}_genes.tab | cut -f 1 > ${IDS}
	@echo
	@echo `wc ${IDS}` IDs
	@echo ${IDS}

## Extract a pathway connecting at best the reactions catalyzed by these gene products
DB=MetaCyc
NETWORK=MetaCyc_directed_141
NETWORK_FILE=${RSAT}/public_html/data/metabolic_networks/networks/${DB}/${NETWORK}.txt
RES_PREFIX=${RES_DIR}/${QUERY}_IDs_Escherichia_coli_strain_K12_-83333_${NETWORK}_pred_pathways_annot
one_pathway:
	@echo "Extracting pathways from genes ${IDS}"
	pathway-extractor -v ${V} -i ${IDS} \
		-g ${NETWORK_FILE} \
		-ger ${RSAT}/data/metabolic_networks/GER_files/GER_Uniprot_${ORG}.tab \
		-o ${RES_DIR} \
		-t temp_dir
	@echo
	@echo ${RES_PREFIX}*

## View the dot-formatted network
VIEWER_DOT=dotty
view_dot:
	${VIEWER_DOT} ${RES_PREFIX}.dot

## View the  png-formatted image
VIEWER_PNG=open -a Preview
view_png:
	${VIEWER_PNG} ${RES_PREFIX}.png



################################################################
##
## EC co-occurrences
##
## Detect pairs of EC which co-occur frequently in proteins. 
PFAM_LINKS_DIR=/no_backup/MICROME_RESULTS/data/pfam_links/
UNIPROT_EC=${PFAM_LINKS_DIR}/uniprot2ec.tab
EC_NAMES=${PFAM_LINKS_DIR}/EC_names.tab
EC_EC_DIR=results/ec_vs_ec
EC_VS_EC=${EC_EC_DIR}/EC_vs_EC
ec_vs_ec:
## TO BE ADDED TO THE COMMAND BELOW
##	-rnames ${EC_NAMES}.tab
	@mkdir -p ${EC_EC_DIR}
	compare-classes -v ${V} -i ${UNIPROT_EC} \
		-triangle -distinct \
		-lth QR 1 \
		-return occ,freq,jac_sim,proba,rank \
		-lth sig 0 \
		-quick ${OPT} \
		-o ${EC_VS_EC}.tab	
	@echo ${EC_VS_EC}.tab
	@text-to-html -i ${EC_VS_EC}.tab \
		-o ${EC_VS_EC}.html
	@echo ${EC_VS_EC}.html


################################################################
## Establish the link between PFAM domains and EC numbers.
##
## We use a tricky way to establish this link by comparing two
## association tables:
## 1. Uniprot Accession - EC
## 2. Uniprot Accession - PFAM
##
## The link cannot be established unambigusously from the original
## Uniprot records, because a single protein can be associated to
## several domains (e.g. P00562 (AK2H_ECOLI): aspartokinase+homoserine
## dehydrogenase). Saddly enough the links from Uniprot to PFAM are
## documented at the level of the protein as a whole, and not at the
## level of the domain. Regions are documented separately, but the
## name associated to a region (e.g. Aspartokinase) is not the same as
## the one associated to the PFAM record (e.g. PF00696. AA_kinase).
##
## The trick used here is to link PFAM to EC numbers via Uniprot
## entries using compare-classes in order to compute the mutual
## coverage for each PFAM - EC association.
PFAM_NAMES=${PFAM_LINKS_DIR}/PFAM_names.tab
UNIPROT_PFAM=${PFAM_LINKS_DIR}/uniprot2pfam.tab
PFAM_EC_DIR=results/pfam_vs_ec
PFAM_VS_EC=${PFAM_EC_DIR}/PFAM_vs_EC
pfam_vs_ec:
## TO BE ADDED TO THE COMMAND BELOW
##	-names ${PFAM_NAMES}.tab -rnames ${EC_NAMES}.tab
	@mkdir -p ${PFAM_EC_DIR}
	compare-classes -v ${V} \
		-q ${UNIPROT_PFAM} \
		-r ${UNIPROT_EC} \
		-lth QR 1 \
		-return occ,freq,jac_sim,proba,rank \
		-quick ${OPT} \
		-o ${PFAM_VS_EC}.tab
	@echo ${PFAM_VS_EC}.tab
	@text-to-html -i ${PFAM_VS_EC}.tab \
		-o ${PFAM_VS_EC}.html
	@echo ${PFAM_VS_EC}.html

