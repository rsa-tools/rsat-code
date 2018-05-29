################################################################
## Testers for RSAT Web services

include ${RSAT}/makefiles/util.mk
MAKEFILE=makefiles/web_services_demo.mk

PYTHON=python

CLIENT_DIR=public_html/web_services/clients
SERVER=http://pedagogix-tagc.univ-mrs.fr/rsat
CLIENT_SCRIPT=${CLIENT_DIR}/python/matrix_scan_soap.py

param:
	@echo "Parameters"
	@echo "	SERVER		${SERVER}"
	@echo "	CLIENT_SCRIPT	${CLIENT_SCRIPT}"
	@echo "	CLIENT_DIR	${CLIENT_DIR}"

matrix_scan:
	${PYTHON} ${CLIENT_SCRIPT} \
		--sequence_file public_html/demo_files/Dmelanogaster_eve_up5000.fasta \
		--matrix_file public_html/demo_files/Dmelanogaster_segmentation_12matrices.tf \
		--uth_pval 1e-3 --server ${SERVER} \
		${OPT}
