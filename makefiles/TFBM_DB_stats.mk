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

INFO_DISTRIB=${DB_NAME}_param_table_log${LOG_BASE}_info_distrib
IMG_FORMAT=pdf
info_distrib_one_db:
	@classfreq -i ${DB_STATS}  -col 20 -ci 1 -v 1 -o ${INFO_DISTRIB}.tsv
	@echo "	${INFO_DISTRIB}.tsv"
	@XYgraph -v 1 -i ${INFO_DISTRIB}.tsv \
		-xcol 1 -ycol 4,5,6  \
		-xsize 300 -ysize 200 \
		-col 20 -ci 1 -lines -format ${IMG_FORMAT} \
		-xleg1 "Information content per motif" \
		-yleg1 "Number of TF binding motifs" \
		-title "${DB_NAME}" \
		${OPT} -o ${INFO_DISTRIB}.${IMG_FORMAT}
	@echo "	${INFO_DISTRIB}.${IMG_FORMAT}"

DATABASES=`grep -v '^;' ${DB_TABLE} | grep -v '^\#' | cut -f 1 | sort | xargs`
DB_TASK= param_table_one_db info_distrib_one_db
param_table_all_dbs:
	@for db in ${DATABASES}; do \
		${MAKE} DB_NAME=$${db} ${DB_TASK}; \
	done
