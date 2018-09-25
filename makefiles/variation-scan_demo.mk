################################################################
## Demonstrator for variation-scan and related programs.
include ${RSAT}/makefiles/util.mk
MAKEFILE=${RSAT}/makefiles/variation-scan_demo.mk

V=2
SPECIES=Homo_sapiens
#ASSEMBLY=GRCh37
ASSEMBLY=GRCh38

## The folder does not depend on the ensembl version anymore.  The
## following option should be used only for specific purpose (when
## installing several ensembl versions of the same assembly).
#SPECIES_SUFFIX=ensembl${ENSEMBL_RELEASE}
#SPECIES_SUFFIX_OPT=-species_suffix ${SPECIES_SUFFIX}
SPECIES_SUFFIX=
SPECIES_SUFFIX_OPT=

################################################################
## The test matrix comes from Weirauch et al.
MATRIX=${RSAT}/public_html/demo_files/variation_demo_set_MWeirauch_cell_2014_15SNPs_TFs.tf

## Variants selected to illustate the typology of cases, including
## non-trivial cases with >2 variants.
DEMO_DIR=${RSAT}/public_html/demo_files/
#VARSEQ_DEMO_FILE=${DEMO_DIR}/${VARIANTS}.varseq
VARIANTS=variation_demo_set_MWeirauch_cell_2014_15SNPs
#VARIANTS=variation_demo_set
VARIANT_INDEL=variation_testfile
ORG=${SPECIES}_${ASSEMBLY}
VARIATION_DIR=${RSAT}/public_html/data/genomes/${ORG}/variations


################################################################
## List parameters
list_param:
	@echo
	@echo "variation-scan demo"
	@echo "	SPECIES		${SPECIES}"
	@echo "	ASSEMBLY	${ASSEMBLY}"
	@echo "	ORG		${ORG}"
	@echo "	VARIATION_DIR	${VARIATION_DIR}"
	@echo "	VARIANTS	${VARIANTS}"
	@echo "	DEMO_VARIANTS	${DEMO_VARIANTS}"
	@echo "	BED_VARIANTS	${BED_VARIANTS}"
	@echo "	VARIANT_IDS	${VARIANT_IDS}"

################################################################
## Count the number of variations per chromosome/contig for the selected organism
variation_stats:
	@echo "Statistics about variations"
	@echo "	SPECIES		${SPECIES}"
	@echo "	ASSEMBLY	${ASSEMBLY}"
	@echo "	ORG		${ORG}"
	@echo "	VARIATION_DIR	${VARIATION_DIR}"
	@echo "Computing number of lines per contig (this can take time)"
	@wc -l ${VARIATION_DIR}/*.varBed | sort -n


## Create result directory
mk_result_dir:
	@mkdir -p ${RESULT_DIR}

################################################################
## Convert variations from VCF (variation X file) format into the
## format supported as input by RSAT retrieve-var.
#RESULT_DIR=results/variation_scan_demo/${VARIANTS}
RESULT_DIR=results/variation_scan_demo
VARIANT_FORMAT_IN=vcf
DEMO_VARIANTS=${DEMO_DIR}/${VARIANT_INDEL}.${VARIANT_FORMAT_IN}
VARIANT_FORMAT_OUT=varBed
PHASED=
CONVERT_VAR_CMD=convert-variations \
	-i ${DEMO_VARIANTS}  \
	-v ${V} -from ${VARIANT_FORMAT_IN} -to ${VARIANT_FORMAT_OUT} ${PHASED}\
	-o ${RESULT_DIR}/${VARIANT_INDEL}${PHASED}.${VARIANT_FORMAT_OUT}
convert_var:
	@echo ""
	@echo "Converting variations from ${VARIANT_FORMAT_IN} to ${VARIANT_FORMAT_OUT}"
	@echo "${CONVERT_VAR_CMD}"
	@${CONVERT_VAR_CMD}
	@echo "Converted variation file"
	@echo "	${RESULT_DIR}/${VARIANT_INDEL}.${VARIANT_FORMAT_OUT}"

convert_var_phased:
	${MAKE} convert_var PHASED=-phased

################################################################
## get variation information from either variant IDs or bed file region.

## Bed file from the demo directory
BED_VARIANTS=${DEMO_DIR}/Ballester_etal_elife_2014_module_beyondprimates_conserved_hg18_lift_to_hg19.bed
VAR_FROM_BED_OUT=${RESULT_DIR}/Ballester_etal_elife_2014_module_beyondprimates_conserved_hg18_lift_to_hg19
VARIANT_IDS=${DEMO_DIR}/variation_demo_set_MWeirauch_cell_2014_15SNPs_IDs.txt
VAR_FROM_ID_OUT=${RESULT_DIR}/variation_demo_set_MWeirauch_cell_2014_15SNPs


## Generic parameters for variation-info command
VAR_INFO_CMD=variation-info -v ${V}\
	-species ${SPECIES} \
	-release ${ENSEMBL_RELEASE} \
	-assembly ${ASSEMBLY} \
	${SPECIES_SUFFIX_OPT} ${OPT}

## Get the variations that overlap a set of genomic regions specificed
## in a BED file. In this example we use a set of peaks from Ballester
## et al., 2010.
VAR_INFO_BED_CMD=${VAR_INFO_CMD} \
	-i ${BED_VARIANTS} \
	-format bed \
	-o ${VAR_FROM_BED_OUT}.varBed
varinfo_from_bed_regions: mk_result_dir
	@echo ""
	@echo "Getting variation information from genomic region (input bed file)."
	@echo "BED_VARIANTS	${BED_VARIANTS}"
	@echo "${DATE}	${VAR_INFO_BED_CMD}"
	@${VAR_INFO_BED_CMD}
	@echo "${DATE}	Collected variations from bed file";
	@echo "Output file: "
	@wc -l ${VAR_FROM_BED_OUT}.varBed


## Get variations from a list of user-specified IDs.
VAR_INFO_ID_CMD=${VAR_INFO_CMD} \
	-i ${VARIANT_IDS} \
	-format id \
	-o ${VAR_FROM_ID_OUT}.varBed
varinfo_from_ids: mk_result_dir
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

## The retrieve-var command can be used in various modalities:
##
## (1) provide a varBed file, which contains a description of one
## variation per row.
##
## (2) specifying genomic regions in a bed file. The program starts by
## identifying the variations that overlap the regions of the bed
## file, and then retrieve their sequences.
##
## (3) a list of variant identifiers provided in a text file (with one
## ID per row). The program then extracts the information about each
## specified variation (varBed info) and then their sequences.
##
## RETRIEVE_VAR_CMD is the common part of the command, we then
## specified the modality-specific parameters:
## RETRIEVE_VAR_CMD_VARBED and
## -release ${ENSEMBL_RELEASE}
RETRIEVE_VAR_CMD=retrieve-variation-seq  \
	-v ${V} \
	-species ${SPECIES} \
	-assembly ${ASSEMBLY} \
	${SPECIES_SUFFIX_OPT}

VARSEQ_DEMO=${DEMO_DIR}/${VARIANTS}.varseq
VARSEQ_OUT=${RESULT_DIR}/${VARIANTS}.varseq

RETRIEVE_VAR_CMD_VARBED=${RETRIEVE_VAR_CMD} \
	-i ${RESULT_DIR}/${VARIANTS}.varBed \
	-mml 30 -format varBed \
	-o ${VARSEQ_OUT}
retrieve_varseq_from_varBed: mk_result_dir
	@echo ""
	@echo "Retrieving variation sequences from variation info file"
	@echo "Input file	${RESULT_DIR}/${VARIANTS}.varBed"
	@echo "${RETRIEVE_VAR_CMD_VARBED}"
	@${RETRIEVE_VAR_CMD_VARBED}
	@echo "Out file"
	@echo "	${VARSEQ_OUT}"


## Retrieve sequences of the variations that overlap a set of genomic
## coordinates specified in a bed file.
RETRIEVE_VAR_CMD_BED=${RETRIEVE_VAR_CMD} \
	-i  ${BED_VARIANTS}\
	-mml 30 -format bed \
	-o ${VAR_FROM_BED_OUT}.varseq
retrieve_var_bed: mk_result_dir
	@echo "${RETRIEVE_VAR_CMD_BED}"
	@${RETRIEVE_VAR_CMD_BED}
	@echo "Out file"
	@echo "${VAR_FROM_BED_OUT}.varseq"
	@echo "Counting number of alleles per variation"
	@grep -v '^;' ${VAR_FROM_BED_OUT}.varseq | grep -v '^\#' \
		| cut -f 5 | sort | uniq -c | sort -k 1 -n \
		| classfreq -v 1 -ci 1 -col 1 \
		>  ${VAR_FROM_BED_OUT}_alleles_per_variant.tab
	@echo ${VAR_FROM_BED_OUT}_alleles_per_variant.tab

## Retrieve sequences of the variations specified by a list of IDs.
RETRIEVE_VAR_CMD_ID=${RETRIEVE_VAR_CMD} \
	-i  ${VARIANT_IDS} \
	-mml 30 -format id \
	-o ${VAR_FROM_ID_OUT}.varseq
retrieve_var_id: mk_result_dir
	@echo "${RETRIEVE_VAR_CMD_ID}"
	@${RETRIEVE_VAR_CMD_ID}
	@echo "Out file"
	@echo "	${VAR_FROM_ID_OUT}.varseq"


################################################################
## Scan selected variations with the matrix of interest
PVAL=0.001
PVAL_RATIO=10
BG_MODEL=public_html/demo_files/all_human_ENCODE_DNAse_mk1_bg.ol
VARSCAN_RES=${RESULT_DIR}/${VARIANTS}_scan_pval${PVAL}_pvalratio${PVAL_RATIO}
VARSCAN_CMD=variation-scan -v ${V} \
	-i ${VARSEQ_DEMO} \
	-m ${MATRIX} -bg ${BG_MODEL} \
	-uth pval ${PVAL} \
	-lth pval_ratio ${PVAL_RATIO} \
	-o ${VARSCAN_RES}.tab
variation_scan: mk_result_dir
	@echo "${VARSCAN_CMD}"
	@${VARSCAN_CMD}
	@echo "Output file"
	@echo "	${VARSCAN_RES}.tab"


## Scan variation with all motifs in JASPAR core vertebrate database
## (~200 motifs)
VARSCAN_RES_JASPAR=${RESULT_DIR}/${VARIANTS}_vs_JASPAR_rsat_var_scan_pval${PVAL}_pvalratio${PVAL_RATIO}
scan_variations_with_jaspar: mk_result_dir
	@echo ""
	@echo "Scanning variations with all motifs from JASPAR core vertebrate"
	@${MAKE} variation_scan \
		MATRIX=${RSAT}/public_html/motif_databases/JASPAR/jaspar_core_vertebrates_2015_03.tf \
		VARSCAN_RES=${VARSCAN_RES_JASPAR}
	@echo

## Scan variations from Weireauch et al. (Cell., 2014) with Jaspar core Vertebrates
JASPAR_CORE_VERTEBRATE=${RSAT}/public_html/motif_databases/JASPAR/jaspar_core_vertebrates_2015_03.tf
WEIRAUCH_VARSEQ=public_html/demo_files/variation_demo_set_MWeirauch_cell_2014_15SNPs.varseq
WEIRAUCH_JASPAR=${RESULT_DIR}/varscan_weirauch-snps_vs_JASPAR_pval${PVAL}_pvalratio${PVAL_RATIO}
varscan_weireauch_with_jaspar: mk_result_dir
	@echo ""
	@echo "Scanning variations with all motifs from JASPAR core vertebrate"
	@echo "	JASPAR_CORE_VERTEBRATE	${JASPAR_CORE_VERTEBRATE}"
	@echo "	WEIRAUCH_VARSEQ		${WEIRAUCH_VARSEQ}"
	@variation-scan  -v ${V} \
		-m_format transfac \
		-m ${JASPAR_CORE_VERTEBRATE} \
		-i ${WEIRAUCH_VARSEQ} \
		-bg ${RSAT}/public_html/data/genomes/${SPECIES}_${ASSEMBLY}/oligo-frequencies/3nt_upstream-noorf_${SPECIES}_${ASSEMBLY}-ovlp-1str.freq.gz \
		-lth score 1 \
		-lth w_diff 1 \
		-lth pval_ratio ${PVAL_RATIO} \
		-uth pval ${PVAL} \
		-o ${WEIRAUCH_JASPAR}.tab
	@echo "	${WEIRAUCH_JASPAR}.tab"
	@text-to-html -i ${WEIRAUCH_JASPAR}.tab \
		-o ${WEIRAUCH_JASPAR}.html
	@echo "	${WEIRAUCH_JASPAR}.html"


## Scan variations from Weireauch et al. (Cell., 2014) with Cisbp core Vertebrates
WEIRAUCH_CISBP=${RESULT_DIR}/varscan_weirauch-snps_vs_cisBP_pval${PVAL}_pvalratio${PVAL_RATIO}
varscan_weireauch_with_cisbp: mk_result_dir
	@echo ""
	@echo "Scanning variations with all motifs from CISBP core vertebrate"
	@variation-scan  -v ${V} \
		-m_format transfac \
		-m ${RSAT}/public_html/motif_databases/cisBP/cisBP_${SPECIES}_2014-10.tf
		-i ${WEIRAUCH_VARSEQ} \
		-bg ${RSAT}/public_html/data/genomes/${SPECIES}_${ASSEMBLY}/oligo-frequencies/3nt_upstream-noorf_${SPECIES}_${ASSEMBLY}-ovlp-1str.freq.gz \
		-lth score 1 \
		-lth w_diff 1 \
		-lth pval_ratio ${PVAL_RATIO} \
		-uth pval ${PVAL} \
		-o ${WEIRAUCH_CISBP}.tab
	@echo "	${WEIRAUCH_CISBP}.tab"
	@text-to-html -i ${WEIRAUCH_CISBP}.tab \
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
compare_regvar: mk_result_dir
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

convert_varScan_from_isRSNP_to_varScan: mk_result_dir
	@echo ${CONVERT_VARSCAN_FROM_isRSNP_CMD}
	@${CONVERT_VARSCAN_FROM_isRSNP_CMD}

convert_var_Scan_from_varScan_to_isRSNP: mk_result_dir
	@echo ${CONVERT_VARSCAN_FROM_varScan_CMD}
	@${CONVERT_VARSCAN_FROM_varScan_CMD}
