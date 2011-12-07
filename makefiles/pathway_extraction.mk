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
EC_EC_DIR=results/ec_vs_ec${SUFFIX}
EC_VS_EC=${EC_EC_DIR}/EC_vs_EC
ec_vs_ec:
	@mkdir -p ${EC_EC_DIR}
	compare-classes -v ${V} -i ${UNIPROT_EC} \
		-rnames ${EC_NAMES} \
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
	@${MAKE} GRAPH=${EC_VS_EC} graph clusters_from_graph



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
##
## We temporarily impose a threshold on intersection (QR >= 5) and
## Jaccard similarity (sim_jac >= 10%) to reduce the number of
## time-consuming computations for the hypergeometric significance,
## but this could be released in a second time.
##
PFAM_NAMES=${PFAM_LINKS_DIR}/Pfam_names.tab
UNIPROT_PFAM=${PFAM_LINKS_DIR}/uniprot2pfam.tab
PFAM_EC_DIR=results/pfam_vs_ec${SUFFIX}
PFAM_VS_EC=${PFAM_EC_DIR}/PFAM_vs_EC
pfam_vs_ec:
	@mkdir -p ${PFAM_EC_DIR}
	compare-classes -v ${V} \
		-q ${UNIPROT_PFAM} -qnames ${PFAM_NAMES} \
		-r ${UNIPROT_EC} -rnames ${EC_NAMES} \
		-lth QR 5 -lth jac_sim 0.1 -lth sig 0 \
		-return occ,freq,jac_sim,proba,rank \
		-quick ${OPT} \
		-o ${PFAM_VS_EC}.tab
	@echo ${PFAM_VS_EC}.tab
	@text-to-html -i ${PFAM_VS_EC}.tab \
		-o ${PFAM_VS_EC}.html
	@echo ${PFAM_VS_EC}.html
	${MAKE} graph GRAPH=${PFAM_VS_EC}

it:
	@${MAKE} GRAPH=${PFAM_VS_EC} graph clusters_from_graph

################################################################
## quick test for ec_vs_ec or pfam_vs_ec
TOP_QUICK=10000
TEST_TASK=pfam_vs_ec
quick_test:
	${MAKE} V=2 OPT='-max_lines ${TOP_QUICK}' SUFFIX='_top${TOP_QUICK}'  ${TEST_TASK}


################################################################
## Convert association table into a graph
GRAPH=${EC_VS_EC}
graph: graph_gml graph_dot

graph_gml:
	convert-graph -i  ${GRAPH}.tab -from tab -to gml -scol 1 -tcol 2 -wcol 19 -ewidth -ecolors fire -min 0 -max 1 -o ${GRAPH}.gml
	@echo ${GRAPH}.gml

graph_dot:
	convert-graph -i  ${GRAPH}.tab -from tab -to dot -scol 1 -tcol 2 -wcol 19 -ewidth -ecolors fire -min 0 -max 1 -o ${GRAPH}.dot
	@echo ${GRAPH}.dot

################################################################
## Extract clusters from an association graph (EC versus EC, PFAM
## versus EC, ...)
clusters_from_graph:
	convert-graph -from tab -to tab -wcol 19 -scol 1 -tcol 2 -i  ${GRAPH}.tab -o ${GRAPH}_graph_for_mcl.tab
	@echo ${GRAPH}_graph_for_mcl.tab
	mcl ${GRAPH}_graph_for_mcl.tab -I 2.5 --abc -V all -o ${GRAPH}_mcl_clusters.mcl
	@echo ${GRAPH}_mcl_clusters.mcl
	convert-classes -from mcl -to tab -i ${GRAPH}_mcl_clusters.mcl -o ${GRAPH}_mcl_clusters.tab
	@echo ${GRAPH}_mcl_clusters.tab
	contingency-table -v 1 -margin -col1 2 -col2 1 -i ${GRAPH}_mcl_clusters.tab -o ${GRAPH}_mcl_cluster_sizes.tab
	@echo ${GRAPH}_mcl_cluster_sizes.tab
	classfreq -v 1 -col 2 -ci 1 -i ${GRAPH}_mcl_cluster_sizes.tab -o ${GRAPH}_mcl_cluster_size_distrib.tab 
	@echo ${GRAPH}_mcl_cluster_size_distrib.tab 
	XYgraph -format png -title1 'Cluster size distribution' -lines -xleg1 'Cluster size' -yleg1 'Number of clusters' -xmin 0 -xcol 2 -ycol 4 \
		-i ${GRAPH}_mcl_cluster_size_distrib.tab \
		-o ${GRAPH}_mcl_cluster_size_distrib.png
	@echo ${GRAPH}_mcl_cluster_size_distrib.png
	XYgraph -format png -title1 'Cluster size distribution' -lines -xleg1 'Cluster size' -yleg1 'Number of clusters' -xmin 0 -xcol 2 -ycol 4 -ylog 2 -xlog 2 \
		-i ${GRAPH}_mcl_cluster_size_distrib.tab \
		-o ${GRAPH}_mcl_cluster_size_distrib_xylog2.png
	@echo ${GRAPH}_mcl_cluster_size_distrib_xylog2.png
