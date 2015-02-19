################################################################
## Demonstrator for variation-scan and related programs.
include ${RSAT}/makefiles/util.mk
MAKEFILE=${RSAT}/makefiles/variation-scan_demo.mk

V=2
SPECIES=Homo_sapiens
ASSEMBLY=GRCh37
ENSEMBL_VERSION=75

## The folder does not depend on the ensembl version anymore.  The
## following option should be used only for specific purpose (when
## installing several ensembl versions of the same assembly).
#SPECIES_SUFFIX=ensembl${ENSEMBL_VERSION}
#SPECIES_SUFFIX_OPT=-species_suffix ${SPECIES_SUFFIX}
SPECIES_SUFFIX=
SPECIES_SUFFIX_OPT=

################################################################
## The test matrix comes from Ballester et al.
MATRIX=${RSAT}/public_html/demo_files/do798+do735_mmus_hnf6_liver.transfac

## Variants selected to illustate the typology of cases, including
## non-trivial cases with >2 variants.
DEMO_DIR=${RSAT}/public_html/demo_files/
VARIANTS=variation_demo_set

################################################################
## Count the number of variations per chromosome/contig for the selected organism
ORG=${SPECIES}_${ASSEMBLY}
VARIATION_DIR=${RSAT}/public_html/data/genomes/${ORG}/variations
variation_stats:
	@echo "Statistics about variations"
	@echo "	SPECIES		${SPECIES}"	
	@echo "	ASSEMBLY	${ASSEMBLY}"
	@echo "	ORG		${ORG}"	
	@echo "	VARIATION_DIR	${VARIATION_DIR}"
	@echo "Computing number of lines per contig (this can take time)"
	@(cd ${VARIATION_DIR}; wc -l *.varBed | sort -n)

all: convert_var 

################################################################
## Convert variations from VCF (variation X file) format into the
## format supported as input by RSAT retrieve-var.
RESULT_DIR=results/variation_scan_demo
VARIANT_FORMAT_IN=vcf
VARIANT_FORMAT_OUT=varBed
CONVERT_VAR_CMD=convert-variations \
	-i ${DEMO_DIR}/${VARIANTS}.${VARIANT_FORMAT_IN}  \
	-e_version ${ENSEMBL_VERSION} \
	-v ${V} -from ${VARIANT_FORMAT_IN} -to ${VARIANT_FORMAT_OUT} \
	-o ${RESULT_DIR}/${VARIANTS}.${VARIANT_FORMAT_OUT}
convert_var:
	@echo ""
	@echo "Converting variations from ${VARIANT_FORMAT_IN} to ${VARIANT_FORMAT_OUT}"
	@mkdir -p ${RESULT_DIR}
	@echo "${CONVERT_VAR_CMD}"
	@${CONVERT_VAR_CMD}
	@echo "Converted variation file"
	@echo "	${RESULT_DIR}/${VARIANTS}.${VARIANT_FORMAT_OUT}"

################################################################
## get variation information from IDs or from a bed file region

## Bed file from the demo directory
BED_VARIANTS=${DEMO_DIR}/Ballester_etal_elife_2014_module_beyondprimates_conserved_hg18_lift_to_hg19.bed
VAR_FROM_BED_OUT=${RESULT_DIR}/Ballester_etal_elife_2014_module_beyondprimates_conserved_hg18_lift_to_hg19
ID_VARIANTS=${DEMO_DIR}/variation_demo_set_MWeirauch_cell_2014_15SNPs_IDs.txt
VAR_FROM_ID_OUT=${RESULT_DIR}/variation_demo_set_MWeirauch_cell_2014_15SNPs

VAR_INFO_CMD=variation-info -v ${V}\
	-species ${SPECIES} \
	-e_version ${ENSEMBL_VERSION} \
	-a_version ${ASSEMBLY} \
	${SPECIES_SUFFIX_OPT} 

VAR_INFO_BED_CMD=${VAR_INFO_CMD} \
	-i ${BED_VARIANTS} \
	-format bed \
	-o ${VAR_FROM_BED_OUT}.varBed

VAR_INFO_ID_CMD=${VAR_INFO_CMD} \
	-i ${ID_VARIANTS} \
	-format id \
	-o ${VAR_FROM_ID_OUT}.varBed

get_var_from_bed:
	@echo "${VAR_INFO_BED_CMD}"
	@${VAR_INFO_BED_CMD}
	@echo "Out file"
	@echo "	${VAR_FROM_BED_OUT}"

get_var_from_ID:
	@echo "${VAR_INFO_ID_CMD}"
	@${VAR_INFO_ID_CMD}
	@echo "Out file"
	@echo "	${VAR_FROM_ID_OUT}"



################################################################
## Retrieve the sequences surrounding a set of input variations
RETRIEVE_VAR_CMD=retrieve-variation-seq  \
	-v ${V} \
	-species ${SPECIES} \
	-e_version ${ENSEMBL_VERSION} \
	-a_version ${ASSEMBLY} \
	${SPECIES_SUFFIX_OPT} 

RETRIEVE_VAR_CMD_VARBED=${RETRIEVE_VAR_CMD} \
	-i ${RESULT_DIR}/${VARIANTS}.varBed \
	-mml 30 -format varBed \
	-o ${RESULT_DIR}/${VARIANTS}_rsat_var.varSeq

retrieve_var_varbed:
	@echo "${RETRIEVE_VAR_CMD_VARBED}"
	@${RETRIEVE_VAR_CMD_VARBED}
	@echo "Out file"
	@echo "	${RESULT_DIR}/${VARIANTS}_rsat_var.varSeq"


RETRIEVE_VAR_CMD_BED=${RETRIEVE_VAR_CMD} \
	-i  ${BED_VARIANTS}\
	-mml 30 -format bed \
	-o ${VAR_FROM_BED_OUT}.varSeq

retrieve_var_bed:
	@echo "${RETRIEVE_VAR_CMD_BED}"
	@${RETRIEVE_VAR_CMD_BED}
	@echo "Out file"
	@echo "${VAR_FROM_BED_OUT}.varSeq"



RETRIEVE_VAR_CMD_ID=${RETRIEVE_VAR_CMD} \
	-i  ${ID_VARIANTS} \
	-mml 30 -format id \
	-o ${VAR_FROM_ID_OUT}.varSeq

retrieve_var_id:
	@echo "${RETRIEVE_VAR_CMD_ID}"
	@${RETRIEVE_VAR_CMD_ID}
	@echo "Out file"
	@echo "	${VAR_FROM_ID_OUT}.varSeq"



################################################################
## Scan selected variations with the matrix of interest
PVAL=0.1
PVAL_RATIO=2
BG_MODEL=public_html/demo_files/all_human_ENCODE_DNAse_mk1_bg.ol
VAR_SCAN_RES=${RESULT_DIR}/${VARIANTS}_rsat_var_scan_pval${PVAL}_pvalratio${PVAL_RATIO}
VAR_SCAN_CMD=variation-scan -v ${V} \
	-i ${RESULT_DIR}/${VARIANTS}_rsat_var.varSeq \
	-m ${MATRIX} -bg ${BG_MODEL} \
	-uth pval ${PVAL} \
	-lth pval_ratio ${PVAL_RATIO} \
	-o ${VAR_SCAN_RES}.tab
variation_scan:
	@echo "${VAR_SCAN_CMD}"
	@${VAR_SCAN_CMD}
	@echo "Output file"
	@echo "	${VAR_SCAN_RES}.tab"


