################################################################
## Calculate phylogenetic profiles for a selected organism and taxon

include ${RSAT}/makefiles/util.mk
MAKEFILE=makefiles/phylo_profiles.mk

TAXON=Fungi
ORG=Saccharomyces_cerevisiae
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

one_taxon:  genus_species profiles profiles_sig profile_pairs sig_vs_MI

################################################################
## Identify all the putative orthologs (BBH) 
RES_DIR=results/${ORG}/${TAXON}
ORTHO=${RES_DIR}/${ORG}_vs_${TAXON}_bbh_len50_ident30_e1e-10
ortho:
	@mkdir -p ${RES_DIR}
	@echo ${ORTHO}.tab
	get-orthologs -v 2 \
		-i ${CDS} \
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
profiles:
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

