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
	@echo ""
	@echo "Collecting supported ${TAXON} at ${URL}"
	@echo "	SERVER	${SERVER}"
	@echo "	TAXON	${TAXON}"
#	@echo "	OPT	${OPT}"
	@mkdir -p ${RES_DIR}
	@supported-organisms-server -taxon ${TAXON} -url ${URL} ${OPT} \
		-return last_update,source,ID,taxonomy \
		-o ${SUPPORTED_FILE}
	@wc -l ${SUPPORTED_FILE}

all: fungi plants prokaryotes protists metazoa teaching

fungi:
#	@${MAKE} supported_taxon TAXON=Fungi SERVER=rsat-tagc.univ-mrs.fr SERVER_PATH=rsat
	@${MAKE} supported_taxon TAXON=Metazoa SERVER=rsat.france-bioinformatique.fr SERVER_PATH=fungi

plants:
	@${MAKE} supported_taxon TAXON=Viridiplantae SERVER=rsat.eead.csic.es SERVER_PATH=plants
#	@${MAKE} supported_taxon TAXON=Viridiplantae SERVER=floresta.eead.csic.es SERVER_PATH=rsat

prokaryotes:
	@${MAKE} supported_taxon TAXON=Bacteria SERVER=embnet.ccg.unam.mx SERVER_PATH=rsat
	@${MAKE} supported_taxon TAXON=Archaea SERVER=embnet.ccg.unam.mx SERVER_PATH=rsat

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
	@${MAKE} supported_taxon TAXON=Metazoa SERVER=rsat.france-bioinformatique.fr SERVER_PATH=metazoa


teaching:
#	@${MAKE} supported_taxon TAXON=Organisms SERVER=pedagogix-tagc.univ-mrs.fr SERVER_PATH=rsat
	@${MAKE} supported_taxon TAXON=Metazoa SERVER=rsat.france-bioinformatique.fr SERVER_PATH=teaching

## Generic server with all organisms hosted by the Institut Francais de Bioinformatique
generic:
	@${MAKE} supported_taxon TAXON=Metazoa SERVER=rsat.france-bioinformatique.fr SERVER_PATH=rsat



