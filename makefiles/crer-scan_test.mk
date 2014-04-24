################################################################
## Tester for crer-scan

include ${RSAT}/makefiles/util.mk
MAKEFILE=${RSAT}/makefiles/crer-scan_test.mk

SITES_DIR=data
SITES_DORSAL=${SITES_DIR}/Drosophila_melanogaster_all_upstream3000-noorf_Dorsal_mkv2_pval0.001_sites.ft

RES_DIR=crer-scan_test/results
SITES=${SITES_DORSAL}
demo_dorsal:
	mkdir -p ${RES_DIR}
	python3 python-scripts/crer_scan.py  \
		-i ${SITES_DORSAL} -s -in_format ft \
		-lth_crer_size 30 -uth_crer_size 500 -lth_crer_sites 2 -uth_crer_sites 1000 \
		-lth_crer_sites_distance 1 -uth_crer_sites_distance 1000 \
		-uth_crer_eval 0.0001 -uth_crer_pval 0.0001 \
		-lth_crer_sig 2.0 -uth_site_pval 0.0001 \
		-lth_score 0 -uth_score 1000 -uth_overlap 1 -return_limits_filtered -o None -v 2 


CRER_DIR=${RES_DIR}/crers
demo_eve:
	@echo "Detectiong CRERS	${SITES_DIR}/matrix-scan_2014-03_segmentation_DEMO1.ft"
	mkdir -p ${CRER_DIR}
	python3 python-scripts/crer_scan.py -v ${V} \
		-i ${RSAT}/public_html/demo_files/Drosophila_melanogaster_eve_segmentation_sites_pval0.001.ft \
		-s \
		-in_format ft \
		-lth_crer_size 30 -uth_crer_size 500 \
		-lth_crer_sites 2 -uth_crer_sites 1000 \
		-lth_crer_sites_distance 1 -uth_crer_sites_distance 1000 \
		-uth_site_pval 0.0001\
		-uth_crer_eval 0.0001 \
		-uth_crer_pval 0.0001 \
		-lth_crer_sig 2 \
		-lth_score 0 -uth_score 1000 \
		-uth_overlap 1 \
		-return_limits_filtered \
		-o ${CRER_DIR}/Drosophila_melanogaster_eve_segmentation_sitespval0.001_crer.ft
	@echo "	 ${CRER_DIR}/Drosophila_melanogaster_eve_segmentation_sitespval0.001_crer.ft"


