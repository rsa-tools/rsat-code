################################################################
## Run a set of diagnostic test to ensure that all the pieces of RSAT
## have been properly installed.

include ${RSAT}/makefiles/util.mk
MAKEFILE=makefiles/install_tests.mk

all: 	test_dir  \
	path \
	os \
	perl_modules \
	supported_organisms \
	r_version \
	purge_seq \
	crer_scan_python2 \
	crer_scan_python3 \
	matrix_clustering \
	ws_stub \
	ws_stub_test \
	ws_nostub_test \
	zip

TEST_DIR=./install_tests
test_dir:
	@mkdir -p ${TEST_DIR}

## RSAT path
PATH_FILE=${TEST_DIR}/RSAT_path.txt
path: test_dir
	@echo
	@echo "Checking RSAT_path"
	@uname -mrs > ${PATH_FILE}
	@echo "	${PATH_FILE}"

## Operating system properties
OS_FILE=${TEST_DIR}/operating_system.txt
os: test_dir
	@echo
	@echo "Checking operating system"
	@uname -mrs > ${OS_FILE}
	@echo "	${OS_FILE}"

## Check Perl modules
NB_MISSING_PERL=`wc -l ${TEST_DIR}/perl_modules_check_missing.txt`
perl_modules: test_dir
	@echo
	@echo "Checking Perl modules"
	@echo "	${TEST_DIR}/perl_modules_check_log.txt"
	@echo "	${TEST_DIR}/perl_modules_check_err.txt"
	@make -f ${RSAT}/makefiles/install_rsat.mk perl_modules_check \
		1> ${TEST_DIR}/perl_modules_check_log.txt \
		2> ${TEST_DIR}/perl_modules_check_err.txt
	@awk '($$1=="Fail") && ($$2 != "Object::InsideOut")' check_perl_modules_eval.txt > ${TEST_DIR}/perl_modules_check_missing.txt
	@echo "	Missing modules	${NB_MISSING_PERL}"

NB_SUPPORTED_ORGANISMS=`cat ${TEST_DIR}/supported_organisms.tsv | wc -l`
supported_organisms:
	@echo
	@echo "	${TEST_DIR}/supported_organisms.tsv"
	@echo "	${TEST_DIR}/supported_organisms_err.txt"
	@supported-organisms -return last_update,source,ID \
		-o ${TEST_DIR}/supported_organisms.tsv \
		2> ${TEST_DIR}/supported_organisms_err.txt
	@echo "Supported organisms	${NB_SUPPORTED_ORGANISMS}"

r_version:
	@echo
	@echo "R version"
	@echo "	${TEST_DIR}/R_version.txt"
	@echo "	${TEST_DIR}/R_version_err.txt"
	@R --version|grep '^R version' | awk '{print $3}' \
		1> ${TEST_DIR}/R_version.txt \
		2> ${TEST_DIR}/R_version_err.txt

ws_stub:
	@echo
	@echo "	${TEST_DIR}/ws_stub_log.txt"
	@echo "	${TEST_DIR}/ws_stub_msg.txt"
	@make -f ${RSAT}/makefiles/init_rsat.mk  ws_init ws_stub \
		1> ${TEST_DIR}/ws_stub_log.txt \
		2> ${TEST_DIR}/ws_stub_msg.txt


ws_stub_test:
	@echo
	@echo "	${TEST_DIR}/ws_stub_test_result.txt"
	@echo "	${TEST_DIR}/ws_stub_test_msg.txt"
	@make -f ${RSAT}/makefiles/init_rsat.mk  ws_init ws_stub_test \
		1> ${TEST_DIR}/ws_stub_test_result.txt \
		2> ${TEST_DIR}/ws_stub_test_msg.txt

ws_nostub_test:
	@echo
	@echo "	${TEST_DIR}/ws_nostub_test_result.txt"
	@echo "	${TEST_DIR}/ws_nostub_test_msg.txt"
	@make -f ${RSAT}/makefiles/init_rsat.mk  ws_init ws_nostub_test \
		1> ${TEST_DIR}/ws_nostub_test_result.txt \
		2> ${TEST_DIR}/ws_nostub_test_msg.txt

purge_seq:
	@echo
	@echo "Testing random-seq | purge-sequence"
	@echo "	${TEST_DIR}/rand_purged.fa"
	@echo "	${TEST_DIR}/rand_purged_err.txt"
	@random-seq -l 100 -n 10 | purge-sequence \
		1> ${TEST_DIR}/rand_purged.fa \
		2> ${TEST_DIR}/rand_purged_err.txt

crer_scan_python2:
	@echo
	@echo "crer-scan with python 2.7"
	@echo "	${TEST_DIR}/crer_scan_python2_log.txt"
	@echo "	${TEST_DIR}/crer_scan_python2_err.txt"
	@make  PYTHON=python2.7 -f makefiles/crer-scan_test.mk  demo_eve \
		1> ${TEST_DIR}/crer_scan_python2_log.txt \
		2> ${TEST_DIR}/crer_scan_python2_err.txt 
	@echo "	results/crer-scan_test//crers/Drosophila_melanogaster_eve_segmentation_sites_pval0.001_crer.ft"

crer_scan_python3:
	@echo
	@echo "crer-scan with python 3"
	@echo "	${TEST_DIR}/crer_scan_python3_log.txt"
	@echo "	${TEST_DIR}/crer_scan_python3_err.txt"
	@make  PYTHON=python3 -f makefiles/crer-scan_test.mk  demo_eve \
		1> ${TEST_DIR}/crer_scan_python3_log.txt \
		2> ${TEST_DIR}/crer_scan_python3_err.txt 
	@echo "	results/crer-scan_test//crers/Drosophila_melanogaster_eve_segmentation_sites_pval0.001_crer.ft"

matrix_clustering: test_dir
	@echo
	@echo "Running matrix-clustering demo"
	@echo "	${TEST_DIR}/matrix-clustering_log.txt"
	@echo "	${TEST_DIR}/matrix-clustering_err.txt"
	@make -f ${RSAT}/makefiles/matrix-clustering_demo.mk \
		cluster_peakmotifs_Oct4 \
		1> ${TEST_DIR}/matrix-clustering_log.txt \
		2> ${TEST_DIR}/matrix-clustering_err.txt

## Create a zip archive with the test results
ARCHIVE=install_tests_${TIME}.zip
zip:  test_dir
	@echo
	@echo "Creating archive with results and log files"
	@echo "	${ARCHIVE}"
	@zip -ry --quiet ${ARCHIVE} ${TEST_DIR}


