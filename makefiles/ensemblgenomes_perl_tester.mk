################################################################
## Examples of utilization for the Perl scripts establishing a
## connection to Ensembl and Ensembl Genomes databases.

include ${RSAT}/makefiles/util.mk
MAKEFILE=makefiles/ensemblgenomes_perl_tester.mk

DB=ensembl
#DB=ensemblgenomes
RESULT_DIR=results/ensemblgenomes
AVAILABLE_SPECIES=${RESULT_DIR}/available_species_${DB}_release${ENSEMBL_RELEASE}_${DAY}.txt
available_species:
	@mkdir -p ${RESULT_DIR}
	@echo "Collecting available species at ${DB}	release ${ENSEMBL_RELEASE}"
	install-ensembl-genome -v ${V} -db ${DB} -available_species \
		-version ${ENSEMBL_RELEASE} \
		-o ${AVAILABLE_SPECIES}
	@echo "	${AVAILABLE_SPECIES}"

AVAILABLE_PER_TAXID=${RESULT_DIR}/available_ensemblgenomes_${QUERY_TYPE}_${TAXID}_${DB}_release${ENSEMBL_RELEASE}_${DAY}.txt
QUERY_TYPE=taxid
TAXID=4751
available_per_taxid:
	@mkdir -p ${RESULT_DIR}
	supported-organisms-ensemblgenomes -v ${V} \
		-query_type taxid -q ${TAXID} \
		-o ${AVAILABLE_PER_TAXID}
	@echo "	${AVAILABLE_PER_TAXID}"

ORG=Saccharomyces_cerevisiae
install_one_species:
	install-ensembl-genome -v ${V} -db ${DB} -species ${ORG}

install_yeast:
	${MAKE} install_one_species DB=ensemblgenomes ORG=Saccharomyces_cerevisiae

install_coli:
	${MAKE} install_one_species DB=ensemblgenomes ORG=Escherichia_coli_str_k_12_substr_mg1655


AVAILABLE_GROUP=${RESULT_DIR}/available_species_${DB}_release${ENSEMBL_RELEASE}_${DAY}_${ORG_GROUP}.txt
select_one_group:
	@echo
	@echo "Selecting group	${ORG_GROUP}	from db ${DB}	${AVAILABLE_SPECIES}"
	grep -i ${ORG_GROUP} ${AVAILABLE_SPECIES} > ${AVAILABLE_GROUP}
	@echo "	${AVAILABLE_GROUP}"

ORG_GROUP=Fungi
install_one_group:
	@echo
	@echo "Installing group ${ORG_GROUP} from db ${DB}"
	install-ensembl-genome -v ${V} -db ${DB} -species_file ${AVAILABLE_GROUP} -nodie
	@echo "Installed group ${ORG_GROUP} from db ${DB}"

select_fungi:
	@${MAKE} select_one_group DB=ensemblgenomes ORG_GROUP=Fungi

install_fungi:
	@${MAKE} install_one_group DB=ensemblgenomes ORG_GROUP=Fungi
