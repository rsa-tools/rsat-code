################################################################
## Load GO (gene ontology) in a mysql database

include ${RSAT}/makefiles/util.mk

WHEN=now
MAKEFILE=${RSAT}/makefiles/go.mk
V=1

USER=go
PASS=go
DB=go
MYSQL=mysql -u ${USER} -D ${DB} -p${PASS}

## Create a 'go' database (as mysql root)
## And grant all permissions to the user 'go'
drop_db:
	mysql -u root -e "drop database ${DB};"

create_db:
	mysql -u root -e "create database ${DB};"
	mysql -u root -e "grant ALL on ${USER}.* to ${DB}@localhost identified by '${PASS}';"


################################################################
## Create an load tables from the TERMDB (GO terms)
TERMDB_DIR=~/go/go_200407-termdb-tables
TABLE=term
TABLE_SCRIPT=${TERMDB_DIR}/${TABLE}.sql
TABLES=`ls ${TERMDB_DIR}/*.sql | xargs basename | perl -pe 's|.sql||g'`
list_tables:
	@echo ${TABLES}

all_tables: list_tables
	@for t in ${TABLES}; do ${MAKE} one_table TABLE=$${t}; done

one_table: create_one_table load_one_table

create_one_table:
	@echo "Creating table	${TABLE}"
	${MYSQL} < ${TABLE_SCRIPT}

TABLE_DATA=${TERMDB_DIR}/${TABLE}.txt
load_one_table:
	@echo "Loading table	${TABLE}"
	${MYSQL} -e 'LOAD DATA LOCAL INFILE "${TABLE_DATA}" INTO TABLE ${TABLE};'


################################################################
## Extract GO term names and relationships for RSAT

TERM=${TERMDB_DIR}/term.txt
TERM_RSAT=${RSAT}/data/GO/term.tab
TERM2TERM=${TERMDB_DIR}/term2term.txt
TERM2TERM_RSAT=${RSAT}/data/GO/term2term.tab
go_for_rsat:
	@echo ";id	name" >  ${TERM_RSAT}
	cut -f 1,2 ${TERM} >> ${TERM_RSAT}
	@echo "term stored in ${TERM_RSAT}"

	@echo ";parent	child" >  ${TERM2TERM_RSAT}
	cut -f 3,4 ${TERM2TERM} >> ${TERM2TERM_RSAT}
	@echo "term2term stored in ${TERM2TERM_RSAT}"

