################################################################
## Compute some statistics about the transcription factor binding
## motifs stored in one of the TF databases supported on this RSAT
## instance.

include ${RSAT}/makefiles/util.mk
MAKEFILE=${RSAT}/makefiles/TFBM_DB_stats.mk

list_param:
	@echo
	@echo "Parameters"
	@echo "	DB_NAME		${DB_NAME}"
	@echo "	DB_PATH		${DB_PATH}"
	@echo "	DB_FORMAT	${DB_FORMAT}"
	@echo "	MATRICES	${MATRICES}"
	@echo "	LOG_BASE	${LOG_BASE}"
	@echo "	DB_STATS	${DB_STATS}"
	@echo "	DATABASES	${DATABASES}"

LOG_BASE=2
DB_NAME=jaspar_core_nonredundant_fungi
DB_TABLE=${RSAT}/public_html/motif_databases/db_matrix_files.tab
DB_PATH=`awk '$$1=="${DB_NAME}" {print $$3}' ${DB_TABLE}`
DB_FORMAT=`awk '$$1=="${DB_NAME}" {print $$2}' ${DB_TABLE}`
MATRICES=${RSAT}/public_html/motif_databases/${DB_PATH}
DB_STATS=${DB_NAME}_param_table_log${LOG_BASE}.tsv
param_table_one_db:
	@echo "Parameter table for TFBM database	${DB_NAME}"
	@convert-matrix -v ${V} -i ${MATRICES} -from ${DB_FORMAT} -to param_table -return parameters -base ${LOG_BASE} -o ${DB_STATS}
	@echo "	${DB_STATS}"

DATABASES=`grep -v '^;' ${DB_TABLE} | grep -v '^\#' | cut -f 1 | sort | xargs`
param_table_all_dbs:
	@for db in ${DATABASES}; do \
		${MAKE} DB_NAME=$${db} param_table_one_db; \
	done
