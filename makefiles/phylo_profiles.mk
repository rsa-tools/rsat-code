################################################################
## Calculate phylogenetic profiles for a selected organism and taxon

include ${RSAT}/makefiles/util.mk
MAKEFILE=${RSAT}/makefiles/phylo_profiles.mk

#TAXON=Fungi
#ORG=Saccharomyces_cerevisiae

TAXON=Bacteria
#ORG=Escherichia_coli_K12
ORG=Escherichia_coli_K_12_substr__MG1655_uid57779
TAXON=Bacteria
DEPTH=5
#RES_DIR=results/phylo_profiles/${ORG}/${TAXON}
RES_ROOT=results
RES_DIR=${RES_ROOT}/phylo_profiles/${ORG}_${TAXON}_depth${DEPTH}
ORG_DIR=${RES_ROOT}/${ORG}


## Threshold on BLAST identity
TH_ID=30
## Threshold on BLAST alignment length
TH_LEN=50
## Threshold on BLAST expect (e-value)
TH_EXPECT=1e-10

list_param:
	@echo "ORG		${ORG}"
	@echo "TAXON		${TAXON}"
	@echo "DEPTH		${DEPTH}"
	@echo "TH_ID		${TH_ID}"
	@echo "TH_EXPECT	${TH_EXPECT}"
	@echo "TH_LEN		${TH_LEN}"
	@echo "RES_DIR		${RES_DIR}"


all_taxa:
	${MAKE} one_taxon TAXON=Fungi
	${MAKE} one_taxon TAXON=Archaea
	${MAKE} one_taxon TAXON=Bacteria
	${MAKE} one_taxon TAXON=Archaea
	${MAKE} merged_ortho
	${MAKE} merged_profiles

one_taxon: ortho genus_species profiles_prev profiles_sig profile_pairs sig_vs_MI
#one_taxon: genus_species profiles profiles_sig profile_pairs sig_vs_MI

################################################################
## Identify all the putative orthologs (BBH) 
ORTHO=${RES_DIR}/${ORG}_vs_${TAXON}_bbh_len${TH_LEN}_ident${TH_ID}_e${TH_EXPECT}
ortho:
	@echo
	@echo "Detecting all orthologs (BBH)	${ORG}	${TAXON}"
	@mkdir -p ${RES_DIR}
	@echo ${ORTHO}.tab
	get-orthologs -v 2 \
		-all \
		-org ${ORG} \
		-taxon ${TAXON} \
		-uth rank 1 -lth ali_len ${TH_LEN} -lth ident ${TH_ID} -uth e_value ${TH_EXPECT} \
		-return e_value,bit_sc,ident,ali_len \
		-o ${ORTHO}.tab
	@echo "	${ORTHO}.tab"


################################################################
## Generate a tab-delimited file with the genus (col1), species (col2)
## and both names (col3) for all the species found in the orthology
## table. 
GENUS_SPECIES=${RES_DIR}/${TAXON}_genus_species.tab
genus_species:
	@echo
	@echo "Generating a genus-species file"
	grep -v '^;' ${ORTHO}.tab | cut -f 2 | sort -u | perl -pe 's|_|\t|' | awk '{print $$1"\t"$$2"\t"$$1"_"$$2}' > ${GENUS_SPECIES}
	@echo ${GENUS_SPECIES}

################################################################
## Convert ortholog relationships into a table of phylogenetic profiles
PROFILES=${ORTHO}_profiles_evalue
profiles_prev:
	convert-classes -v 2 \
		-i ${ORTHO}.tab \
		-from tab -to profiles \
		-ccol 2 -mcol 3 -scol 4 -null "NA" \
		-o ${PROFILES}.tab
	@echo ${PROFILES}.tab

################################################################
## Convert ortholog relationships into a table of phylogenetic
## profiles, using sig=-log10(eval) as similarity metrics between
## orthologs
PROFILES_SIG=${ORTHO}_profiles_sig
profiles_sig:
	awk -F '\t' '$$1 !~ ";" && $$1 !~"#" {print $$2"\t"$$3"\t"(-log($$4)/log(10))}' ${ORTHO}.tab \
	| convert-classes -v 2 \
		-from tab -to profiles \
		-ccol 1 -mcol 2 -scol 3 -null "NA" \
		-o ${PROFILES_SIG}.tab
	@echo ${PROFILES_SIG}.tab


################################################################
## Generate a neetwork of co-occurrence (presence/absence)
PROFILE_PAIRS=${ORTHO}_profile_pairs
PRIMARY_NAMES=${ORG_DIR}/cds_primary_names.tab
INFINITE=300
profile_pairs:
	@mkdir -p ${ORG_DIR}
	grep primary ${CDS_NAMES} > ${PRIMARY_NAMES}
	compare-profiles -v 2 -i ${PROFILES}.tab \
		-distinct -return counts,jaccard,hyper,entropy \
		-na "NA" -inf ${INFINITE} -lth AB 2 -lth sig 0 \
		-names ${PRIMARY_NAMES} \
		-o ${PROFILE_PAIRS}.tab
	@echo ${PROFILE_PAIRS}.tab

################################################################
## Compare mutual information with the hypergeometric significance
IMG_FORMAT=pdf
sig_vs_MI:
	XYgraph  -format ${IMG_FORMAT} -xcol 27 -ycol 15,18 -legend \
		-i ${PROFILE_PAIRS}.tab \
		-xleg1 "I(A,B)" -yleg1 "hypergeometric significance" \
		-title1 "${SUFFIX}" \
		-o ${PROFILE_PAIRS}_sig_vs_MI.${IMG_FORMAT} 
	@echo ${PROFILE_PAIRS}_sig_vs_MI.${IMG_FORMAT} \


################################################################
## Merge the profiles from the trhee main taxa
MERGED_DIR=${ORG_DIR}/all
MERGED_ORTHO=${MERGED_DIR}/${ORG}_vs_all_bbh_len${TH_LEN}_ident${TH_ID}_e${TH_EXPECT}
TAXA=Fungi Bacteria Archaea
merged_ortho:
	@mkdir -p ${MERGED_DIR}
	@echo "" > ${MERGED_ORTHO}.tab
	@for taxon in ${TAXA}; do \
		${MAKE} _merge_one_taxon TAXON=$${taxon} ; \
	done
	@echo ${MERGED_ORTHO}.tab

_merge_one_taxon:
	@echo "	Concatening file	${ORTHO}.tab"
	grep -v '^;' ${ORTHO}.tab >> ${MERGED_ORTHO}.tab

merged_profiles:
	@${MAKE} TAXON=all genus_species 
	@${MAKE} TAXON=all profiles
	@${MAKE} TAXON=all profile_pairs


################################################################
################################################################
################################################################
## THIS IS A SECOND SCRIPT TO OBTAIN PHYLO PROFILES, I SHOULD CHECK
## THE REDUNDANT TASKS AND CLEAN
################################################################
################################################################


################################################################
## Generate phylogenetic profiles

## Iterate over organisms
ORGANISMS=Escherichia_coli_K_12_substr__MG1655_uid57779 Bacillus_subtilis_168_uid57675 Pseudomonas_putida_KT2440_uid57843 Salmonella_enterica_serovar_Typhimurium_LT2_uid57799

## Run all the tasks for one organism
one_org: select_species bbh profiles gene_names compa sig_vs_mi network_genes result_summary


################################################################
## Select a set of species in a way to reduce redundancy between
## closely related species. This is a bit tricky: we cut the organism
## tree at a given depth (e.g. 5) and select a single species of each
## taxon at this depth.
SPECIES=${RES_DIR}/selected_species_${TAXON}_depth${DEPTH}.tab
select_species:
	@echo "Selecting species	${TAXON}	depth=${DEPTH}"
	@mkdir -p ${RES_DIR}
	@supported-organisms -return ID,taxonomy | awk '$$1 == "${ORG}"' > ${SPECIES}
	@supported-organisms -return ID,taxonomy -taxon ${TAXON} | perl -pe 's|; |\t|g' \
		| cut -f 1-${DEPTH} | sort -k 2 -u | perl -pe 's|\t|; |g' | perl -pe 's|; |\t|' \
		>> ${SPECIES}
	@echo "species	`wc -l ${SPECIES}`"


################################################################
## Identify all the putative orthologs according to the criterion of
## bidirectional best hits (BBH)
V=2
SUFFIX=${ORG}_vs_${TAXON}_eval_${TH_EXPECT}_ident${TH_ID}_len${TH_LEN}
BBH=${RES_DIR}/bbh_${SUFFIX}
bbh:
	@echo
	@echo "Detecting all orthologs (BBH)	${ORG}	${TAXON}"
	@mkdir -p ${RES_DIR}
	get-orthologs -v ${V} \
		-i ${CDS} \
		-org ${ORG} \
		-org_list ${SPECIES} \
		-uth rank 1 -lth ali_len ${TH_LEN} -lth ident ${TH_ID} -uth e_value ${TH_EXPECT} \
		-return e_value,bit_sc,ident,ali_len \
		-o ${BBH}.tab
	@echo "BBH	${BBH}.tab"
#		-taxon ${TAXON}


################################################################
## Convert ortholog table into a profile table 
## with the IDs of the putative orthologs 
profiles: profiles_id profiles_boolean profiles_evalue

PROFILES=${RES_DIR}/profiles_${ORG}_vs_${TAXON}_eval_${TH_EXPECT}_ident${TH_ID}_len${TH_LEN}
profiles_id:
	@echo
	@echo "Computing phylogenetic profiles (gene IDs) from BBH"
	@wc -l ${BBH}.tab
	convert-classes -v 2 -i ${BBH}.tab \
		-from tab -to profiles \
		-ccol 2 -mcol 3 -scol 1 -null "<NA>" \
		| grep -v '^;' \
		> ${PROFILES}_ids.tab
	@echo "ID profiles	${PROFILES}_ids.tab"

## Convert ortholog table into a Boolean profile table 
profiles_boolean:
	@echo
	@echo "Computing Boolean phylogenetic profiles from BBH"
	@wc -l ${BBH}.tab
	convert-classes -v 2 -i ${BBH}.tab \
		-from tab -to profiles \
		-ccol 2 -mcol 3  -null 0 \
		| grep -v '^;' \
		> ${PROFILES}_boolean.tab
	add-gene-info -org ${ORG} -before -i ${PROFILES}_boolean.tab -info name > ${PROFILES}_boolean_names.tab
	@echo "Boolean profiles	${PROFILES}_boolean_names.tab"

## Convert ortholog table into a profile table with E-values 
profiles_evalue:
	@echo
	@echo "Computing E-value phylogenetic profiles from BBH"
	@wc -l ${BBH}.tab
	convert-classes -v 2 -i ${BBH}.tab \
		-from tab -to profiles \
		-ccol 2 -mcol 3  -scol 4 -null "NA" \
		| grep -v '^;' \
		> ${PROFILES}_Evalue.tab
	@echo "E-value profiles	${PROFILES}_boolean_names.tab"

# ## Convert the orthology into "classes", where each class (second
# ## column) corresponds to a gene from Saccharomyces cerevisiae, and
# ## indicates the set of genomes (first column) in which this gene is
# ## present.
# profile_classes:
# 	convert-classes -from tab -to tab -mcol 2 -ccol 3 -scol 5 \
# 		-i ${BBH}.tab \
# 		-o ${BBH}_classes.tab
#
#
# ## Compare profiles using the program compare-classes.
# compare-classes -v 3 \
# 	-i ${BBH}_classes.tab \
# 	-lth QR 1 -lth sig 0 -sort sig -sc 3 \
# 	-return occ,proba,dotprod,jac_sim,rank \
# 	-o ${BBH}_gene_pairs.tab


################################################################
## Pairwise comparisons between each gene pair 

################################################################
## Extract the primary name of each gene
CDS=${RSAT}/data/genomes/${ORG}/genome/cds.tab
CDS_NAMES=${RSAT}/data/genomes/${ORG}/genome/cds_names.tab
GENE_NAMES=${RES_DIR}/${ORG}_gene_names.tab
gene_names:
	@echo
	@echo "Getting gene names	${ORG}"
	@mkdir -p ${RES_DIR}
	grep -v '^--'  ${CDS} | cut -f 1 | add-gene-info -org ${ORG} -info name -o ${GENE_NAMES}
	@echo "	${GENE_NAMES}"

## Compare profiles using compare-classes
COMPA=${PROFILES}_compa
SIG_COL=`grep -P '^;\t\d+\tsig' ${COMPA}.tab | cut -f 2`
MIN_SPEC=5
compa:
	@echo
	@echo "Extracting co-occurrence network from phylogenetic profiles"
	@grep -v '^;' ${BBH}.tab | grep -v '^#' \
		| awk '{print $$2"\t"$$3"\t"$$4}' \
		| compare-classes -v ${V}  -i /dev/stdin -sc 3 \
		-return occ,freq,proba,entropy,jac_sim,rank -sort sig \
		-triangle -distinct \
		-lth sig 0 -lth Q ${MIN_SPEC} -lth R ${MIN_SPEC} -lth QR ${MIN_SPEC} \
		-rnames ${GENE_NAMES} -qnames ${GENE_NAMES} \
		-o ${COMPA}.tab
	@echo "	${COMPA}.tab"
	@text-to-html -i ${COMPA}.tab -o ${COMPA}.html
	@echo "	${COMPA}.html"
	@echo
	@echo "Generating network graph	SIG_COL=${SIG_COL}"
	@convert-graph -i ${COMPA}.tab -from tab -to gml -scol 3 -tcol 4 \
		-wcol ${SIG_COL} -ewidth -ecolors fire \
		-o ${COMPA}.gml 
	@echo "	${COMPA}.gml"


################################################################
## Compare mutual information with the hypergeometric significance
IMG_FORMAT=pdf
SIG_COL=`grep -P '^;\t\d+\tsig' ${COMPA}.tab | cut -f 2`
MI_COL=`grep -P '^;\t\d+\tI.Q.R' ${COMPA}.tab | cut -f 2`
sig_vs_mi:
	@echo "Comparing mutual information (${MI_COL}) to hypergeometric significance ${SIG_COL}"
	XYgraph  -format ${IMG_FORMAT} -xcol ${MI_COL} -ycol ${SIG_COL} \
		-i ${COMPA}.tab \
		-xleg1 "I(A,B)" \
		-yleg1 "hypergeometric significance" \
		-title1 "${SUFFIX}" \
		-o ${COMPA}_sig_vs_MI.${IMG_FORMAT} 
	@echo ${COMPA}_sig_vs_MI.${IMG_FORMAT} \

################################################################
## Extract the set of genes fond in the network, with their names and
## description. This file can be loaded in CyctoScape (for example)
## with the function "Import > Attributes from table".
network_genes:
	@echo
	@echo "Extracting description for network genes"
	@grep -v '^;' ${COMPA}.tab | grep -v '^#' | cut -f 1,2 | perl -pe 's|\t|\n|g' | sort -u > ${COMPA}_node_IDs.tab
	@wc -l ${COMPA}_node_IDs.tab
	add-gene-info -org ${ORG} -i ${COMPA}_node_IDs.tab \
		-info name,id,left,right,strand,descr | cut -f 2-10 \
		> ${COMPA}_node_descr.tab
	@echo "	${COMPA}_node_descr.tab"


################################################################
## Compute the distribution of node degrees to detect "hubs"
degree_distrib:
	@echo
	@echo "Computing degree distribution"
	graph-node-degree -v 1 -i ${COMPA}.gml -in_format gml -all -sort -o ${COMPA}_degrees.tab
	@echo "	${COMPA}_degrees.tab"
	classfreq -v 1 -i ${COMPA}_degrees.tab -col 4 -ci 1 -o ${COMPA}_degree_distrib.tab
	@echo "	${COMPA}_degree_distrib.tab"
	@${MAKE} _degree_distrib_graph LOG=-log
	@${MAKE} _degree_distrib_graph LOG=-ylog
	@${MAKE} _degree_distrib_graph LOG=''
	@${MAKE} filter_hubs MAX_DEGREE=50
	@${MAKE} filter_hubs MAX_DEGREE=20

_degree_distrib_graph:
	@XYgraph -i ${COMPA}_degree_distrib.tab \
		-xcol 1 -ycol 4,5,6 -xleg1 "degree" \
		-yleg1 "Number of genes" -title1 "${SUFFIX}" \
		-format ${IMG_FORMAT} \
		${LOG} -lines -pointsize 0 \
		-o ${COMPA}_degree_distrib${LOG}.${IMG_FORMAT}
	@echo "	${COMPA}_degree_distrib${LOG}.${IMG_FORMAT}"

################################################################
## Filter the graph by degree: suppress hubs above a given degree
## threshold.
MAX_DEGREE=50
HUB_LIST=${COMPA}_hubs_deg_gt${MAX_DEGREE}
filter_hubs:
	@echo 
	@echo "Filtering out hubs (degree > ${MAX_DEGREE}"
	grep -v '^;'  ${COMPA}_degrees.tab | grep -v '^#'| awk -F'\t' '$$4 > ${MAX_DEGREE}'  | cut -f 1 > ${HUB_LIST}.tab
	@echo "	${HUB_LIST}.tab"
	@grep -v -f ${HUB_LIST}.tab ${COMPA}.tab > ${COMPA}_maxdeg${MAX_DEGREE}.tab
	@echo "	${COMPA}_maxdeg${MAX_DEGREE}.tab"
	@convert-graph -i ${COMPA}_maxdeg${MAX_DEGREE}.tab -from tab -to gml -scol 3 -tcol 4 \
		-wcol ${SIG_COL} -ewidth -ecolors fire \
		-o ${COMPA}_maxdeg${MAX_DEGREE}.gml 
	@echo "	${COMPA}_maxdeg${MAX_DEGREE}.gml"

## Print a result summary
SUMMARY=${RES_DIR}/result_summary.tab
EDGE_NB=`grep -v '^;' ${COMPA}.tab | grep -v '^\#'| wc -l | awk '{print $$1}'`
NODE_NB=`grep -v '^;' ${COMPA}_node_IDs.tab | grep -v '^\#'| wc -l | awk '{print $$1}'`
SPECIES_NB=`grep -v '^;' ${SPECIES} | grep -v '^\#'| wc -l | awk '{print $$1}'`
result_summary:
	@echo
	@echo "Generating result summary"
	@echo "Result summary" > ${SUMMARY}
	@${MAKE} list_param >> ${SUMMARY}
	@echo "SPECIES_NB	${SPECIES_NB}" >> ${SUMMARY}
	@echo "EDGE_NB		${EDGE_NB}" >> ${SUMMARY}
	@echo "NODE_NB		${NODE_NB}" >> ${SUMMARY}
	@echo "BBH		${BBH}.tab" >> ${SUMMARY}
	@echo "NETWORK		${COMPA}.tab" >> ${SUMMARY}
	@echo "GRAPH		${COMPA}.gml" >> ${SUMMARY}
	@echo "NODE DESCR	${COMPA}_node_descr.tab" >> ${SUMMARY}
	@echo "	${SUMMARY}"
