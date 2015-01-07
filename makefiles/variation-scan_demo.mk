################################################################
## Demonstrator for variation-scan and related programs.

include ${RSAT}/makefiles/util.mk
MAKEFILE=${RSAT}/makefiles/variation-scan_demo.mk

################################################################
## The test matrix comes from Ballester et al.
MATRIX=${RSAT}/public_html/demo_files/do798+do735_mmus_hnf6_liver.transfac

## Variants selected to illustate the typology of cases, including
## non-trivial cases with >2 variants.
VARIANTS=${RSAT}/public_html/demo_files/variation_demo_set

E_VERSION=72
V=2
################################################################
## Convert variations from VCF (variation X file) format into the
## format supported as input by RSAT retrieve-var.
TO=rsat-var
CONVERT_VAR_CMD=convert-variations \
	-i ${VARIANTS}.vcf  \
	-e_version ${E_VERSION} \
	-v ${V} -from vcf -to ${TO} \
	-o ${VARIANTS}.${TO}
convert_var:
	@echo ""
	@echo "Converting variations from VCF to ${TO}"
	@echo "${CONVERT_VAR_CMD}"
	@${CONVERT_VAR_CMD}

################################################################
## Retrieve the sequences surrounding a set of input variations
ORG=Homo_sapiens
SPECIES_SUFFIX=ensembl72
A_VERSION=GRCh37
RETRIEVE_VAR_CMD=retrieve-variation-seq  \
	-v ${V} \
	-species ${ORG} \
	-e_version ${E_VERSION} \
	-a_version ${A_VERSION} \
	-species_suffix ${SPECIES_SUFFIX} \
	-i ${VARIANTS}.rsat-var \
	-mml 30 -format rsat-var \
	-o ${VARIANTS}_rsat_var.seq
retrieve_var:
	@echo "${RETRIEVE_VAR_CMD}"
	@${RETRIEVE_VAR_CMD}
	@echo "Out file" ${VARIANTS}_rsat_var.seq

################################################################
## Scan selected variations with the matrix of interest
PVAL=1
PVAL_RATIO=1
BG_MODEL=public_html/demo_files/all_human_ENCODE_DNAse_mk1_bg.ol
VAR_SCAN_CMD=variation-scan -v ${V} -i ${VARIANTS}_rsat_var.seq -m ${MATRIX} -bg ${BG_MODEL} -uth pval ${PVAL} -lth pval_ratio ${PVAL_RATIO} -o ${VARIANTS}_rsat_var_scan_pval${PVAL}_pvalratio${PVAL_RATIO}.tab
variation_scan:
	@echo "${VAR_SCAN_CMD}"
	@${VAR_SCAN_CMD}


