################################################################
## Run a set of diagnostic test to ensure that all the pieces of RSAT
## have been properly installed.

include ${RSAT}/makefiles/util.mk
MAKEFILE=makefiles/install_tests.mk


all: test_dir path os matrix_clustering zip

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


matrix_clustering:
	@echo
	@echo "Running matrix-clustering demo"
	@make -f ${RSAT}/makefiles/matrix-clustering_demo.mk \
		cluster_peakmotifs_Oct4 \
		1> ${TEST_DIR}/matrix-clustering_log.txt \
		2> ${TEST_DIR}/matrix-clustering_err.txt
	@echo "	${TEST_DIR}/matrix-clustering_log.txt"
	@echo "	${TEST_DIR}/matrix-clustering_err.txt"

## Create a zip archive with the test results
ARCHIVE=install_tests_${TIME}.zip
zip: 
	@echo
	@echo "Creating archive with results and log files"
	@zip -ry --quiet ${ARCHIVE} ${TEST_DIR}
	@echo "	${ARCHIVE}"
