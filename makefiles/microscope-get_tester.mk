################################################################
## Tests for peak-motifs

include ${RSAT}/makefiles/util.mk
MAKEFILE=${RSAT}/makefiles/microscope-get_tester.mk

V=2

DIR_DATA=data
DIR_RESULTS=results/microscope

PYTHON=python2.7

## Collect list of organisms supported in Microscope
ORG_FILE=${DIR_RESULTS}/organisms.tab
organisms:
	@mkdir -p ${DIR_RESULTS}
	@echo
	@echo "Collecting organisms from MICROSCOPE"
	microscope_get organisms
	@echo "${DATE} Organisms collected"
	@echo "	${ORG_FILE}"
