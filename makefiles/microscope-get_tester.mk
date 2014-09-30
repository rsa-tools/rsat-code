################################################################
## Tests for peak-motifs

#ORG=Escherichia_coli_K-12_DH10B
ORG=Escherichia_coli_K12_MG1655
#ORG=Bacillus_subtilis

################################################################
## Variables
V=1
MAKE=make -s -f ${MAKEFILE}
DATE=`date +%Y-%M-%d_%H:%M:%S`


################################################################
## List of targets
usage:
	@echo "usage: make [-OPT='options'] target"
	@echo "implemented targets"
	@perl -ne 'if (/^([a-z]\S+):/){ print "\t$$1\n";  }' ${MAKEFILE}
#include ${RSAT}/makefiles/util.mk
MAKEFILE=${RSAT}/makefiles/microscope-get_tester.mk

V=2

DIR_DATA=data
RES_DIR=results/microscope

PYTHON=python2.7

help:
	microscope_get --help

################################################################
## Collect list of organisms supported in Microscope
ORG_FILE=${RES_DIR}/organisms.tab
organisms:
	@mkdir -p ${RES_DIR}
	@echo
	@echo "Collecting supported organisms from MICROSCOPE"
	microscope_get organisms -o ${ORG_FILE}
	@echo "${DATE} Organisms collected"
	@echo "	${ORG_FILE}"


################################################################
## Get description of all reactions
REACT_DIR=${RES_DIR}/reactions
REACTIONS_ALL=${REACT_DIR}/reactions.tab
reactions:
	@mkdir -p ${REACT_DIR}
	@echo
	@echo "${DATE}	Collecting all reactions"
	microscope_get reactions -o ${REACTIONS_ALL}
	@echo "${DATE}	Reactions collected"
	@echo "	${REACTIONS_ALL}"


################################################################
## Collect gene-protein-reactio table for a given organism
ORG_DIR=${RES_DIR}/${ORG}
GPR=${ORG_DIR}/${ORG}_gpr.tab
gpr:
	@mkdir -p ${ORG_DIR}
	@echo
	@echo "${DATE}	Collecting GPR for organism ${ORG}"
	microscope_get gpr -org ${ORG} -o ${GPR}
	@echo "${DATE}	GPR collected"
	@echo "	${GPR}"

################################################################
## Get list of reactions for one organism
ORG_REACTIONS=${ORG_DIR}/${ORG}_reactions.tab
reactions_one_species:
	@mkdir -p ${ORG_DIR}
	@echo
	@echo "${DATE}	Collecting reactions for organism	${ORG}"
	microscope_get reactionList -org ${ORG} -o ${ORG_REACTIONS}
	@echo "${DATE}	Collected reactions for organism	${ORG}"
	@echo "	${ORG_REACTIONS}"

ORGANISMS=`grep -v '^\#' ${ORG_FILE}|cut -f 2 | xargs`
NB_ORGANISMS=`grep -v '^\#' ${ORG_FILE}|cut -f 2 | wc -l`
list_organisms:
	@echo "Organisms"
	@echo "${ORGANISMS}"

################################################################
## Iterate a task over all organisms
ORG_TASK=gpr reactions_one_species
all_organisms:
	@echo "Iterating over ${NB_ORGANISMS} organisms"
	@for org in ${ORGANISMS}; do \
		${MAKE} ${ORG_TASK} ORG=$$org; \
	done
