################################################################
## Tester for crer-scan

include ${RSAT}/makefiles/util.mk
MAKEFILE=${RSAT}/makefiles/crer-scan_test.mk

DEMO_DIR=${RSAT}/public_html/demo_files
SITES_DORSAL=${DEMO_DIR}/Drosophila_melanogaster_all_upstream3000-noorf_Dorsal_mkv2_pval0.001_sites.ft

## Choose python version. The script is compatible with both 2.7 and
## 3. Version 2.7 at least is required to use argsparse.
PYTHON=python3
#PYTHON=python2.7
RES_DIR=results/crer-scan_test/
CRER_DIR=${RES_DIR}/crers
DORSAL_PREFIX=Drosophila_melanogaster_all_upstream3000-noorf_Dorsal_mkv2_pval0.001_sites${IN_SUFFIX}
SITES=${SITES_DORSAL}
DORSAL_OUT=${CRER_DIR}/${DORSAL_PREFIX}${OUT_SUFFIX}_crer.ft
demo_dorsal:
	mkdir -p ${CRER_DIR}
	${PYTHON} python-scripts/crer_scan.py  -v ${V} \
		-i ${SITES_DORSAL} -s -in_format ft \
		-lth_crer_size 30 -uth_crer_size 500 -lth_crer_sites 2 -uth_crer_sites 1000 \
		-lth_crer_sites_distance 1 -uth_crer_sites_distance 1000 \
		-uth_crer_eval 0.0001 -uth_crer_pval 0.0001 \
		-lth_crer_sig 2.0 -uth_site_pval 0.0001 \
		-lth_score 0 -uth_score 1000 -uth_overlap 1 -return_limits_filtered -o ${DORSAL_OUT}
	@echo "	 ${DORSAL_OUT}"


################################################################
## Scan even-skipped promoter with the binding sites predicted from 12
## TFBM (matrices) of Drosophila segmentation genes.
LTH_SIG=0.1
EVE_PREFIX=Drosophila_melanogaster_eve_segmentation_sites_pval0.001${IN_SUFFIX}
EVE_UPSTREAM_SITES=${DEMO_DIR}/${EVE_PREFIX}.ft
EVE_OUT=${CRER_DIR}/${EVE_PREFIX}${OUT_SUFFIX}_crer.ft
demo_eve:
	@echo "Detecting CRERs for segmentation TF in even-skipped upstream sequence"
	@echo "Sites	${EVE_UPSTREAM_SITES}"
	mkdir -p ${CRER_DIR}
	${PYTHON} python-scripts/crer_scan.py -v ${V} \
		-i  ${EVE_UPSTREAM_SITES} \
		-s \
		-in_format ft \
		-lth_crer_size 100 \
		-uth_crer_size 500 \
		-lth_crer_sites 2 \
		-lth_crer_sites_distance 1 \
		-uth_crer_sites_distance 100 \
		-uth_site_pval 1e-3 \
		-lth_crer_sig ${LTH_SIG} \
		-uth_overlap 1 \
		-return_limits_filtered \
		-o ${EVE_OUT}
	@echo "	 ${EVE_OUT}"


## Run eve demo with a sites file that contains no comment lines
demo_eve_nocomments:
	@${MAKE} demo_eve IN_SUFFIX=_nocomments

## Check if the results with/without comments in the input file are the same
check_eve_nocomments:
	diff ${EVE_OUT} ${CRER_DIR}/${EVE_PREFIX}_nocomments_crer.ft

phython2.7_vs_3:
	@${MAKE} demo_eve PYTHON=python2.7 OUT_SUFFIX=python2.7
	@${MAKE} demo_eve PYTHON=python3 OUT_SUFFIX=python3


