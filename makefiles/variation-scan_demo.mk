################################################################
## Demonstrator for variation-scan and related programs.

include ${RSAT}/makefiles/util.mk
MAKEFILE=${RSAT}/makefiles/variation-scan_demo.mk

SNP=rs2736191
MATRIX=${RSAT}/public_html/demo_files/do798+do735_mmus_hnf6_liver.transfac
VARIANTS=${RSAT}/public_html/demo_files/test_complex_tab

E_VERSION=72
V=3

CONVERT_VAR_CMD=convert-variations -i ${VARIANTS}.vcf  -e_version ${E_VERSION} -v ${V} -from vcf -to rsat-var -o ${VARIANTS}.rsat_var

convert_var:
	@echo "${CONVERT_VAR_CMD}"
	@${CONVERT_VAR_CMD}

ORG=Homo_sapiens
RETRIEVE_VAR_CMD=retrieve-variation-seq  -v ${V} -species ${ORG}  -e_version ${E_VERSION} -i ${VARIANTS}.rsat_var  -mml 30 -o ${VARIANTS}_rsat_var.seq -format rsat-var

retrieve_var:
	@echo "${RETRIEVE_VAR_CMD}"
	@${RETRIEVE_VAR_CMD}

PVAL=1e-3
PVAL_RATIO=10
BG_MODEL=public_html/demo_files/all_human_ENCODE_DNAse_mk1_bg.ol
VAR_SCAN_CMD=variation-scan -i ${VARIANTS}_rsat_var.seq -m ${MATRIX} -bg ${BG_MODEL} -uth pval ${PVAL} -lth pval_ratio ${PVAL_RATIO} -o ${VARIANTS}_rsat_var_scan_pval${PVAL}_pvalratio${PVAL_RATIO}.tab

variation_scan:
	@echo "${VAR_SCAN_CMD}"
	@${VAR_SCAN_CMD}


