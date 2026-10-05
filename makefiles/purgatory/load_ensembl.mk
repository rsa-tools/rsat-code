## ##############################################################
## Load a genome from ENSEMBL MYSQL distribution

include ${RSAT}/makefiles/util.mk
MAKEFILE=${RSAT}/makefiles/load_ensembl.mk

## ENSEMBL directory
ENSEMBL_DIR=${RSAT}/downloads/ftp.ensembl.org/pub/current_worm/data/mysql/caenorhabditis_elegans_core_27_130
CREATE_SCRIPT=`ls -1 ${ENSEMBL_DIR}/*.sql`

## MySQL parameters
DB=ensembl
LOGIN=ensembl
PASSWD=ensembl
MYSQL=mysql -u ${LOGIN} -D ${DB} -p${PASSWD}

## ##############################################################
## Creating database and tables

## This should be executed only once for ever
create_database:
	@echo "Creating ensembl database"
	mysql -u root -e "create database ensembl"
	mysql -u root -e "GRANT ALL ON ${DB}.* TO ${LOGIN}@'%' identified by '${PASSWD}'"

create_tables:
	@echo "Creating tables from script ${CREATE_SCRIPT}"
	cat ${CREATE_SCRIPT} | ${MYSQL} 

## ##############################################################
## Loading
TABLE=gene
TABLE_FILE=${ENSEMBL_DIR}/${TABLE}.txt.table
ALL_TABLES=`ls ${ENSEMBL_DIR}/*.txt.table.gz | perl -pe 's|${ENSEMBL_DIR}/||g' | perl -pe 's|.txt.table.gz||g' | xargs`
TABLES=${ALL_TABLES}
one_table_load:
	@echo "Loading table ${TABLE} from file ${TABLE_FILE}"
	@if [ -f "${TABLE_FILE}.gz" ] ; then						\
		gunzip -f ${TABLE_FILE}.gz ; \
	fi
	${MYSQL} -e "LOAD DATA LOCAL INFILE '${TABLE_FILE}' INTO TABLE ${TABLE};"
	gzip -f ${TABLE_FILE}


tables_load:
	@echo "Loading all tables ${TABLES}"
	@for table in ${TABLES}; do \
		${MAKE} one_table_load TABLE=$${table} ; \
	done

COUNT=`${MYSQL} -e "SELECT COUNT(*) FROM ${TABLE};" | grep -v '^+' | grep -iv 'count'`
one_table_count:
	@echo "${COUNT}	rows in table ${TABLE}"

tables_count:
	@echo "Loading all tables ${TABLES}"
	@for table in ${TABLES}; do \
		${MAKE} one_table_count TABLE=$${table} ; \
	done

## ##############################################################
## deletinf and dropping
one_table_delete:
	@echo "Deleting table ${TABLE} from file ${TABLE_FILE}"
	${MYSQL} -e "DELETE FROM ${TABLE};"

TABLES=`ls ${ENSEMBL_DIR}/*.txt.table.gz | perl -pe 's|${ENSEMBL_DIR}/||g' | perl -pe 's|.txt.table.gz||g' | xargs`
tables_delete:
	@echo "Deleting all tables ${TABLES}"
	@for table in ${TABLES}; do \
		${MAKE} one_table_delete TABLE=$${table} ; \
	done

drop_database:
	@echo "Dropping ensembl database"
	mysql -u root -e "drop database ensembl"
