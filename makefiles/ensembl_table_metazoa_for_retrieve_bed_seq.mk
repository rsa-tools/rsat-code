## ensembl_table_metazoa_for_retrieve_bed_seq.mk makefile
## Usage: make -f ensembl_table_metazoa_for_retrieve_bed_seq.mk

include ${RSAT}/makefiles/util.mk

MAKEFILE=${RSAT}/makefiles/ensembl_table_metazoa_for_retrieve_bed_seq.mk

ORG=Homo_sapiens_GRCh37

CHROMOSOME_PREFIX=chr

## Command to create the equivance name table to idneitfy the raw files for metazoa genomes

CHROMOSOME_FILE_NAME=${RSAT}/data/genomes/${ORG}/genome/chromosome*.raw

CHROMOSOME_EQUIVALENCE_TABLE=${RSAT}/data/genomes/${ORG}/genome/common_to_ensemble_chromosome_names.tab

create_equivalence_table:
	ls -1 ${CHROMOSOME_FILE_NAME} | grep -v repeatmasked | perl -ane 's/.+\/// ; s/\.raw//;chomp; @line=split("_",$$_); print "${CHROMOSOME_PREFIX}".$$line[2]."\t".$$_."\n" ' > ${CHROMOSOME_EQUIVALENCE_TABLE}
	@echo "Equivalence table for" ${ORG} ${CHROMOSOME_EQUIVALENCE_TABLE}

create_equivalence_table_human37:
	${MAKE} create_equivalence_table ORG=Homo_sapiens_GRCh37

create_equivalence_table_human38:
	${MAKE} create_equivalence_table ORG=Homo_sapiens_GRCh38

create_equivalence_table_mouse38:
	${MAKE} create_equivalence_table ORG=Mus_musculus_GRCm38

create_equivalence_table_mouse37:
	${MAKE} create_equivalence_table ORG=Mus_musculus_NCBIM37




