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
	@echo "Collecting supported ${TAXON} at ${URL}"
	@mkdir -p ${RES_DIR}
	@supported-organisms-server -taxon ${TAXON} -url ${URL} ${OPT} \
		-return last_update,source,ID,taxonomy \
		-o ${SUPPORTED_FILE}
	@wc -l ${SUPPORTED_FILE}


all: fungi plants prokaryotes protists metazoa

fungi:
	@${MAKE} supported_taxon TAXON=Fungi SERVER=rsat-tagc.univ-mrs.fr SERVER_PATH=rsat

plants:
	@${MAKE} supported_taxon TAXON=Viridiplantae SERVER=floresta.eead.csic.es SERVER_PATH=rsat

prokaryotes:
	@${MAKE} supported_taxon TAXON=Bacteria SERVER=embnet.ccg.unam.mx SERVER_PATH=rsa-tools
	@${MAKE} supported_taxon TAXON=Archaea SERVER=embnet.ccg.unam.mx SERVER_PATH=rsa-tools

EUKARYOTA_FILE=${RES_DIR}/supported_Eukaryota_rsat01.biologie.ens.fr_rsat.tab
PROTIST_FILE=${RES_DIR}/supported_Protists_rsat01.biologie.ens.fr_rsat.tab
protists:
	@${MAKE} supported_taxon TAXON=Eukaryota SERVER=rsat01.biologie.ens.fr SERVER_PATH=rsat
	@echo "Filtering protists"
	@awk '$$2=="ensemblgenomes"' ${EUKARYOTA_FILE} \
		| grep -v 'Fungi' \
		| grep -v 'Metazoa' \
		| grep -v 'Viridiplantae' \
		> ${PROTIST_FILE}
	@wc -l ${PROTIST_FILE}

metazoa:
	@${MAKE} supported_taxon TAXON=Metazoa SERVER=rsat.sb-roscoff.fr SERVER_PATH=
