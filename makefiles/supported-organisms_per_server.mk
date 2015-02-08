################################################################
## Check the supported organisms on the specific RSAT servers, taking
## in consideration their taxon specificity.

include ${RSAT}/makefiles/util.mk
MAKEFILE=makefiles/supported-organisms_per_server.mk

TAXON=Fungi
SERVER=rsat-tagc.univ-mrs.fr
SERVER_PATH=rsat
URL=http://${SERVER}/${SERVER_PATH}
RES_DIR=results/supported-organisms_per_server
SUPPORTED_FILE=${RES_DIR}/supported_${TAXON}_${SERVER}_${SERVER_PATH}.tab
supported_taxon:
	@echo "Supported ${TAXON} at ${URL}"
	@mkdir -p ${RES_DIR}
	@supported-organisms-server -taxon ${TAXON} -url ${URL} ${OPT} \
		-o ${SUPPORTED_FILE}
	@wc -l ${SUPPORTED_FILE}


fungi:
	@${MAKE} supported_taxon TAXON=Fungi SERVER=rsat-tagc.univ-mrs.fr SERVER_PATH=rsat

plants:
	@${MAKE} supported_taxon TAXON=Viridiplantae SERVER=floresta.eead.csic.es SERVER_PATH=rsat

