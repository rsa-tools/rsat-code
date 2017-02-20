################################################################
## Examples of utilization for the Perl scripts establishing a
## connection to Ensembl and Ensembl Genomes databases.

include ${RSAT}/makefiles/util.mk
MAKEFILE=makefiles/install-ensembl-genomes_demo.mk

list_param:
	@echo
	@echo "Parameters"
	@echo "	ENSEMBL_RELEASE		${ENSEMBL_RELEASE}"
	@echo "	ENSEMBLGENOMES_RELEASE	${ENSEMBLGENOMES_RELEASE}"
	@echo "	DB			${DB}"
	@echo "	RESULT_DIR		${RESULT_DIR}"
	@echo "	AVAILABLE_SPECIES	${AVAILABLE_SPECIES}"
	@echo "	ORG			${ORG}"
	@echo "	QUERY_TYPE		${QUERY_TYPE}"

DB=ensembl
#DB=ensemblgenomes
RESULT_DIR=results/ensemblgenomes
AVAILABLE_SPECIES=${RESULT_DIR}/available_species_${DB}_release${ENSEMBL_RELEASE}_${DAY}.txt
available_species:
	@mkdir -p ${RESULT_DIR}
	@echo "Collecting available species at ${DB}	release ${ENSEMBL_RELEASE}"
	install-ensembl-genome -v ${V} -db ${DB} -available_species \
		-release ${ENSEMBL_RELEASE} \
		-o ${AVAILABLE_SPECIES}
	@echo "	${AVAILABLE_SPECIES}"

available_ensembl:
	@${MAKE} available_species DB=ensembl

available_ensemblgenomes:
	@${MAKE} available_species DB=ensemblgenomes

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

## Shuffle the lines of a species file, in order to avoid staying
## blocked in the hundreds of strains for the same species. This is
## useful to install bacteria.
shuffle_one_group:
	cat ${AVAILABLE_GROUP} bacteria_to_install.txt  |  perl -MList::Util=shuffle -e 'print shuffle(<STDIN>);' > ${AVAILABLE_GROUP_SHUFFLED}
	@echo "Shuffled ${AVAILABLE_GROUP}"
	@echo "	${AVAILABLE_GROUP_SHUFFLED}"

ORG_GROUP=Fungi
install_one_group:
	@echo
	@echo "Installing group ${ORG_GROUP} from db ${DB}"
	install-ensembl-genome -v ${V} -db ${DB} -species_file ${AVAILABLE_GROUP} -nodie
	@echo "Installed group ${ORG_GROUP} from db ${DB}"

install_one_group_shuffled:
	@${MAKE} install_one_group AVAILABLE_GROUP=${AVAILABLE_GROUP_SHUFFLED}

select_fungi:
	@${MAKE} select_one_group DB=ensemblgenomes ORG_GROUP=Fungi

install_fungi:
	@${MAKE} select_one_group DB=ensemblgenomes ORG_GROUP=Fungi
	@${MAKE} install_one_group DB=ensemblgenomes ORG_GROUP=Fungi

select_plants:
	@${MAKE} select_one_group DB=ensemblgenomes ORG_GROUP=Plants

install_plants:
	@${MAKE} select_one_group DB=ensemblgenomes ORG_GROUP=Plants
	@${MAKE} install_one_group DB=ensemblgenomes ORG_GROUP=Plants

select_metazoa:
	@${MAKE} select_one_group DB=ensemblgenomes ORG_GROUP=Metazoa

install_metazoa:
	@${MAKE} select_one_group DB=ensemblgenomes ORG_GROUP=Metazoa
	@${MAKE} install_one_group DB=ensemblgenomes ORG_GROUP=Metazoa


AVAILABLE_GROUP_SHUFFLED=${RESULT_DIR}/available_species_${DB}_release${ENSEMBL_RELEASE}_${DAY}_${ORG_GROUP}_shuffled.txt
install_bacteria_shuffled:
	@${MAKE} select_one_group DB=ensemblgenomes ORG_GROUP=Bacteria
	@${MAKE} shuffle_one_group DB=ensemblgenomes ORG_GROUP=Bacteria
	@${MAKE} install_one_group_shuffled DB=ensemblgenomes ORG_GROUP=Bacteria
