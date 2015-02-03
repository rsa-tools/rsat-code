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
