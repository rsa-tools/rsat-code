############################################################
#
# $Id: install_EnsemblGenomes.mk,v 1.3 2011/04/10 13:29:22 rsat Exp $
#
# Time-stamp: <>
#
############################################################

include ${RSAT}/makefiles/util.mk
MAKEFILE=${RSAT}/makefiles/install_ensembl_genomes.mk

DESCRIPTION="This script shows how to install genomes from Ensembl in RSAT

${MAKE} available_species
	Get the list of species available at Ensembl

${MAKE} available_species
"
descr:
	@echo
	@echo ${DESCRIPTION}
	@echo

param:
	@echo
	@echo "Default parameters"
	@echo "------------------"
	@echo "SPECIES	${SPECIES}"

available_species:


################################################################
## Return the list of organisms supported at ensembl.org
list_supported:
	supported-organisms-ensembl


## One genome
ORG=Pan_troglodytes_EnsEMBL
FROM=-2000

install_one_ensembl_genome:
	install-organism -v ${V} -org ${ORG} -task config -up_from ${FROM}
	install-organism -v ${V} -org ${ORG} -task allup,dyads,oligos,start_stop,upstream_freq	\
		-batch -rm
