################################################################
## Calculate phylogenetic profiles for a selected organism and taxon

include ${RSAT}/makefiles/util.mk
MAKEFILE=${RSAT}/makefiles/phylo_profiles.mk

#TAXON=Fungi
#ORG=Saccharomyces_cerevisiae

TAXON=Bacteria
#ORG=Escherichia_coli_K12
ORG=Escherichia_coli_K_12_substr__MG1655_uid57779


CDS=${RSAT}/data/genomes/${ORG}/genome/cds.tab
CDS_NAMES=${RSAT}/data/genomes/${ORG}/genome/cds_names.tab
TH_ID=30
TH_LEN=50
TH_2=1e-10

all_taxa:
	${MAKE} one_taxon TAXON=Fungi
	${MAKE} one_taxon TAXON=Archaea
	${MAKE} one_taxon TAXON=Bacteria
	${MAKE} one_taxon TAXON=Archaea
	${MAKE} merged_ortho
	${MAKE} merged_profiles

#one_taxon: ortho genus_species profiles profiles_sig profile_pairs sig_vs_MI
one_taxon: genus_species profiles profiles_sig profile_pairs sig_vs_MI

################################################################
## Identify all the putative orthologs (BBH) 
RES_DIR=results/phylo_profiles/${ORG}/${TAXON}
ORTHO=${RES_DIR}/${ORG}_vs_${TAXON}_bbh_len50_ident30_e1e-10
ortho:
	@mkdir -p ${RES_DIR}
	@echo ${ORTHO}.tab
	get-orthologs -v 2 \
		-all \
		-org ${ORG} \
		-taxon ${TAXON} \
		-uth rank 1 -lth ali_len 50 -lth ident 30 -uth e_value 1e-10 \
		-return e_value,bit_sc,ident,ali_len \
		-o ${ORTHO}.tab

GENUS_SPECIES=${RES_DIR}/${TAXON}_genus_species.tab
genus_species:
	grep -v '^;' ${ORTHO}.tab | cut -f 2 | sort -u | perl -pe 's|_|\t|' | awk '{print $$1"\t"$$2"\t"$$0}' > ${GENUS_SPECIES}
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
## Generate a neetwork of co-presence/absence
PROFILE_PAIRS=${ORTHO}_profile_pairs
profile_pairs:
	@grep primary ${CDS_NAMES} > results/${ORG}/cds_primary_names.tab
	compare-profiles -v 2 -i ${PROFILES}.tab \
		-distinct -return counts,jaccard,hyper,entropy \
		-na "NA" -inf Inf -lth AB 2 -lth sig 0 \
		-names results/${ORG}/cds_primary_names.tab \
		-o ${PROFILE_PAIRS}.tab
	@echo ${PROFILE_PAIRS}.tab

################################################################
## Compare mutual information with the hypergeometric significance
IMG_FORMAT=pdf
sig_vs_MI:
	XYgraph  -format ${IMG_FORMAT} -xcol 27 -ycol 15,18 -legend \
		-i ${PROFILE_PAIRS}.tab \
		-xleg1 "I(A,B)" -yleg1 "hypergeometric significance" \
		-title1 "${ORG}_vs_${TAXON}_bbh_len50_ident30_e1e-10" \
		-o ${PROFILE_PAIRS}_sig_vs_MI.${IMG_FORMAT} 
	@echo ${PROFILE_PAIRS}_sig_vs_MI.${IMG_FORMAT} \


################################################################
## Merge the profiles from the trhee main taxa
MERGED_DIR=results/${ORG}/all
MERGED_ORTHO=${MERGED_DIR}/${ORG}_vs_all_bbh_len50_ident30_e1e-10
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

all: gene_names select_species bbh profiles compa

################################################################
## Extract the primary name of each gene
GENE_NAMES=${RES_DIR}/gene_names.tab
CDS=${RSAT}/data/genomes/${ORG}/genome/cds.tab
gene_names:
	@mkdir -p ${RES_DIR}
	grep -v '^--'  ${CDS} | cut -f 1 | add-gene-info -org ${ORG} -info name -o ${GENE_NAMES}
	@echo ${GENE_NAMES}


################################################################
## Select a set of species in a way to reduce redundancy between
## closely related species. This is a bit tricky: we cut the organism
## tree at a given depth (e.g. 5) and select a single species of each
## taxon at this depth.
#ORG=Saccharomyces_cerevisiae
TAXON=Fungi
ORG=Escherichia_coli_K12
TAXON=Bacteria
DEPTH=4
RES_DIR=results/profiles/${ORG}/${TAXON}_depth${DEPTH}
SPECIES=${RES_DIR}/selected_species_${TAXON}_depth${DEPTH}.tab
select_species:
	@echo "Selecting species	${TAXON}	depth=${DEPTH}"
	@mkdir -p ${RES_DIR}
	@supported-organisms -return ID,taxonomy -taxon ${TAXON} | perl -pe 's|; |\t|g' \
		| cut -f 1-${DEPTH} | sort -k 2 -u | perl -pe 's|\t|; |g' | perl -pe 's|; |\t|' \
		> ${SPECIES}
	@echo "species	`wc -l ${SPECIES}`"


################################################################
## Identify all the putative orthologs according to the criterion of
## bidirectional best hits (BBH)
V=2
BBH=${RES_DIR}/bbh_${ORG}_vs_${TAXON}_eval_1e-10_ident30_len50
bbh:
	@mkdir -p ${RES_DIR}
	get-orthologs -v ${V} \
		-i ${CDS} \
		-org ${ORG} \
		-org_list ${SPECIES} \
		-uth rank 1 -lth ali_len 50 -lth ident 30 -uth e_value 1e-10 \
		-return e_value,bit_sc,ident,ali_len \
		-o ${BBH}.tab
	@echo "BBH	${BBH}.tab"
#		-taxon ${TAXON}


################################################################
## Convert ortholog table into a profile table 
## with the IDs of the putative orthologs 
profiles: profiles_id profiles_boolean profiles_evalue

PROFILES=${RES_DIR}/profiles_${ORG}_vs_${TAXON}_eval_1e-10_ident30_len50
profiles_id:
	convert-classes -v 2 -i ${BBH}.tab \
		-from tab -to profiles \
		-ccol 2 -mcol 3 -scol 1 -null "<NA>" \
		-o ${PROFILES}_ids.tab
	@echo "ID profiles	${PROFILES}_ids.tab"

## Convert ortholog table into a Boolean profile table 
profiles_boolean:
	convert-classes -v 2 -i ${BBH}.tab \
		-from tab -to profiles \
		-ccol 2 -mcol 3  -null 0 \
		-o ${PROFILES}_boolean.tab
	add-gene-info -org ${ORG} -before -i ${PROFILES}_boolean.tab -info name -o ${PROFILES}_boolean_names.tab
	@echo "Boolean profiles	${PROFILES}_boolean_names.tab"

## Convert ortholog table into a profile table with E-values 
profiles_evalue:
	convert-classes -v 2 -i ${BBH}.tab \
		-from tab -to profiles \
		-ccol 2 -mcol 3  -scol 4 -null "NA" \
		-o ${PROFILES}_Evalue.tab


################################################################
## Pairwise comparisons between each gene pair 
COMPA=${PROFILES}_compa
SIG_COL=`grep -P '^;\t\d+\tsig' results/profiles/Escherichia_coli_K12/Bacteria_depth5/profiles_Escherichia_coli_K12_vs_Bacteria_eval_1e-10_ident30_len50_compa.tab | cut -f 2`
MIN_SPEC=5
compa:
	grep -v '^;' ${BBH}.tab | grep -v '^#' \
		| awk '{print $$2"\t"$$3"\t"$$4}' \
		| compare-classes -v ${V}  -i /dev/stdin -sc 3 \
		-return occ,freq,proba,entropy,jac_sim,rank -sort sig \
		-triangle -distinct \
		-lth sig 0 -lth Q ${MIN_SPEC} -lth R ${MIN_SPEC} -lth QR ${MIN_SPEC} \
		-rnames ${GENE_NAMES} -qnames ${GENE_NAMES} \
		-o ${COMPA}.tab
	@echo ${COMPA}.tab
	@text-to-html -i ${COMPA}.tab -o ${COMPA}.html
	@echo ${COMPA}.html
	convert-graph -i ${COMPA}.tab -from tab -to gml -scol 3 -tcol 4 \
		-wcol ${SIG_COL} -ewidth -ecolors fire \
		-o ${COMPA}.gml 
	@echo ${COMPA}.gml 

