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
VARIANT_IDS=${DEMO_DIR}/variation_demo_set_MWeirauch_cell_2014_15SNPs_IDs.txt
VAR_FROM_ID_OUT=${RESULT_DIR}/variation_demo_set_MWeirauch_cell_2014_15SNPs

VAR_INFO_CMD=variation-info -v ${V}\
	-species ${SPECIES} \
	-e_version ${ENSEMBL_VERSION} \
	-a_version ${ASSEMBLY} \
	${SPECIES_SUFFIX_OPT} ${OPT}

VAR_INFO_BED_CMD=${VAR_INFO_CMD} \
	-i ${BED_VARIANTS} \
	-format bed \
	-o ${VAR_FROM_BED_OUT}.varBed

VAR_INFO_ID_CMD=${VAR_INFO_CMD} \
	-i ${VARIANT_IDS} \
	-format id \
	-o ${VAR_FROM_ID_OUT}.varBed

varinfo_from_bed_regions:
	@echo ""
	@echo "Getting variation information from genomic region (input bed file)."
	@echo "BED_VARIANTS	${BED_VARIANTS}"
	@echo "${DATE}	${VAR_INFO_BED_CMD}"
	@${VAR_INFO_BED_CMD}
	@echo "${DATE}	Collected variations from bed file";
	@echo "Output file: "
	@wc -l ${VAR_FROM_BED_OUT}.varBed

varinfo_from_ids:
	@echo ""
	@echo "Getting variation information from variant IDs"
	@echo "VARIANT_IDS	${VARIANT_IDS}"
	@echo "${VAR_INFO_ID_CMD}"
	@${VAR_INFO_ID_CMD}
	@echo "${DATE}	Collected variations from ID file";
	@echo "Output file: "
	@wc -l ${VAR_FROM_ID_OUT}.varBed



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
retrieve_varseq_from_varbed:
	@echo ""
	@echo "Retrieving variation sequences from variation info file"
	@echo "Input file	${RESULT_DIR}/${VARIANTS}.varBed"
	@echo "${RETRIEVE_VAR_CMD_VARBED}"
	@${RETRIEVE_VAR_CMD_VARBED}
	@echo "Out file"
	@echo "	${RESULT_DIR}/${VARIANTS}_rsat_var.varSeq"

## Still valid ? To check
RETRIEVE_VAR_CMD_BED=${RETRIEVE_VAR_CMD} \
	-i  ${BED_VARIANTS}\
	-mml 30 -format bed \
	-o ${VAR_FROM_BED_OUT}.varSeq
retrieve_var_bed:
	@echo "${RETRIEVE_VAR_CMD_BED}"
	@${RETRIEVE_VAR_CMD_BED}
	@echo "Out file"
	@echo "${VAR_FROM_BED_OUT}.varSeq"

## Still valid ? To check
RETRIEVE_VAR_CMD_ID=${RETRIEVE_VAR_CMD} \
	-i  ${VARIANT_IDS} \
	-mml 30 -format id \
	-o ${VAR_FROM_ID_OUT}.varSeq
retrieve_var_id:
	@echo "${RETRIEVE_VAR_CMD_ID}"
	@${RETRIEVE_VAR_CMD_ID}
	@echo "Out file"
	@echo "	${VAR_FROM_ID_OUT}.varSeq"


################################################################
## Scan selected variations with the matrix of interest
PVAL=0.001
PVAL_RATIO=10
BG_MODEL=public_html/demo_files/all_human_ENCODE_DNAse_mk1_bg.ol
VARSCAN_RES=${RESULT_DIR}/${VARIANTS}_rsat_var_scan_pval${PVAL}_pvalratio${PVAL_RATIO}
VARSCAN_CMD=variation-scan -v ${V} \
	-i ${RESULT_DIR}/${VARIANTS}_rsat_var.varSeq \
	-m ${MATRIX} -bg ${BG_MODEL} \
	-uth pval ${PVAL} \
	-lth pval_ratio ${PVAL_RATIO} \
	-o ${VARSCAN_RES}.tab
variation_scan:
	@echo "${VARSCAN_CMD}"
	@${VARSCAN_CMD}
	@echo "Output file"
	@echo "	${VARSCAN_RES}.tab"


## Scan variation with all motifs in JASPAR core vertebrate database
## (~200 motifs)
VARSCAN_RES_JASPAR=${RESULT_DIR}/${VARIANTS}_vs_JASPAR_rsat_var_scan_pval${PVAL}_pvalratio${PVAL_RATIO}
scan_variations_with_jaspar:
	@echo ""
	@echo "Scanning variations with all motifs from JASPAR core vertebrate"
	@${MAKE} variation_scan \
		MATRIX=${RSAT}/public_html/motif_databases/JASPAR/jaspar_core_vertebrates_2015_03.tf \
		VARSCAN_RES=${VARSCAN_RES_JASPAR}
	@echo 

## Scan variations from Weireauch et al. (Cell., 2014) with Jaspar core Vertebrates
WEIRAUCH_VARSEQ=public_html/demo_files/variation_demo_set_MWeirauch_cell_2014_15SNPs.var-seq
WEIRAUCH_JASPAR=${RESULT_DIR}/varscan_weirauch-snps_vs_JASPAR_pval${PVAL}_pvalratio${PVAL_RATIO}
varscan_weireauch_with_jaspar:
	@echo ""
	@echo "Scanning variations with all motifs from JASPAR core vertebrate"
	@variation-scan  -v ${V} \
		-m_format transfac \
		-m ${RSAT}/public_html/motif_databases/JASPAR/jaspar_core_vertebrates_2015_03.tf \
		-i ${WEIRAUCH_VARSEQ} \
		-bg ${RSAT}/public_html/data/genomes/Homo_sapiens_GRCh37/oligo-frequencies/3nt_upstream-noorf_Homo_sapiens_GRCh37-ovlp-1str.freq.gz \
		-lth score 1 \
		-lth w_diff 1 \
		-lth pval_ratio ${PVAL_RATIO} \
		-uth pval ${PVAL} \
		-o ${WEIRAUCH_JASPAR}.tab
	@echo "	${WEIRAUCH_JASPAR}.tab"
	@txt-to-html -i ${WEIRAUCH_JASPAR}.tab \
		-o ${WEIRAUCH_JASPAR}.html
	@echo "	${WEIRAUCH_JASPAR}.html"


## Scan variations from Weireauch et al. (Cell., 2014) with Cisbp core Vertebrates
WEIRAUCH_CISBP=${RESULT_DIR}/varscan_weirauch-snps_vs_cisBP_pval${PVAL}_pvalratio${PVAL_RATIO}
varscan_weireauch_with_cisbp:
	@echo ""
	@echo "Scanning variations with all motifs from CISBP core vertebrate"
	@variation-scan  -v ${V} \
		-m_format transfac \
		-m ${RSAT}/public_html/motif_databases/cisBP/cisBP_Homo_sapiens_2014-10.tf
		-i ${WEIRAUCH_VARSEQ} \
		-bg ${RSAT}/public_html/data/genomes/Homo_sapiens_GRCh37/oligo-frequencies/3nt_upstream-noorf_Homo_sapiens_GRCh37-ovlp-1str.freq.gz \
		-lth score 1 \
		-lth w_diff 1 \
		-lth pval_ratio ${PVAL_RATIO} \
		-uth pval ${PVAL} \
		-o ${WEIRAUCH_CISBP}.tab
	@echo "	${WEIRAUCH_CISBP}.tab"
	@txt-to-html -i ${WEIRAUCH_CISBP}.tab \
		-o ${WEIRAUCH_CISBP}.html
	@echo "	${WEIRAUCH_CISBP}.html"


################################################################
## Compare regulatory variations from various sources
## - variation-scan results of Weirauch variations scanned with
##   JASPAR matrices
## - is-rsnp results with Weirauch variations
## - HaploReg results with Weirauch variations
COMPA=${RESULT_DIR}/regvar_comparisons_weinrauch
WEINRAUCH_CISBP=${RESULT_DIR}/weirauch-snps_cisbp
WEINRAUCH_HAPLOREG=${RESULT_DIR}/weirauch-snps_haploreg
compare_regvar:
	compare-reg-var -v ${V} \
		-file weinrauch_jaspar ${WEIRAUCH_JASPAR}.tab \
		-file weinrauch_jaspar2 ${WEIRAUCH_JASPAR}.tab \
		-file weinrauch_jaspar3 ${WEIRAUCH_JASPAR}.tab \
		-o ${COMPA}.tab
	@echo "	${COMPA}.tab"
	text-to-html -i ${COMPA}.tab -o ${COMPA}.html
	@echo "	${COMPA}.html"

## TEMPORARY: WE DON't HAVE THE scisbp and haploreg results yet
## Yvon will collect them from the databases, and we will then convert
## their format to varscan format.
#		-file weinrauch_cisbp ${WEINRAUCH_CISBP}.tab \
#		-file weinrauch_halporeg ${WEINRAUCH_HAPLOREG}.tab \


################
## Test convert-varScan

isRSNP_FILE=${DEMO_DIR}/isRSNP_result.csv
CONV_VARSCAN_RES=${RESULT_DIR}/convert_varScan_results


CONVERT_VARSCAN_FROM_isRSNP_CMD=convert-varScan -i ${isRSNP_FILE} -from isRSNP -to varScan -o ${CONV_VARSCAN_RES}_from_isRSNP_to_varScan.tab

CONVERT_VARSCAN_FROM_varScan_CMD=convert-varScan -i ${VARSCAN_RES}.tab -to isRSNP -from varScan -o ${CONV_VARSCAN_RES}_from_varScan_to_isRSNP.tab

convert_varScan_from_isRSNP_to_varScan:
	@echo ${CONVERT_VARSCAN_FROM_isRSNP_CMD}
	@${CONVERT_VARSCAN_FROM_isRSNP_CMD}

convert_var_Scan_from_varScan_to_isRSNP:
	@echo ${CONVERT_VARSCAN_FROM_varScan_CMD}
	@${CONVERT_VARSCAN_FROM_varScan_CMD}



