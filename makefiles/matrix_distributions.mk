################################################################
## Compute matrix distributions for each matrix of a collection
## (Jaspar, Transfac).

include ${RSAT}/makefiles/util.mk
MAKEFILE=${RSAT}/makefiles/matrix_distributions.mk

DB=JASPAR
MATRIX_DB_FILE=${RSAT}/public_html/data/motif_databases/JASPAR/jaspar_core_vertebrates_2009_10.tf 
MATRIX_DIR=${RSAT}/public_html/data/motif_databases/JASPAR
SPLIT_DIR=${MATRIX_DIR}/separate_matrices_${DB}
SPLIT_PREFIX=${SPLIT_DIR}/${DB}
MATRIX_LIST=${SPLIT_PREFIX}_matrix_list.tab
## Split the matrix database into separate files
split_db:
	@echo
	@echo "Splitting matrix database into separate files"
	@mkdir -p ${SPLIT_DIR}
	convert-matrix -v 1 -from tf -to tf -split -i ${MATRIX_DB_FILE} -o ${SPLIT_PREFIX}
	@echo "${SPLIT_DIR}"
	@echo "${MATRIX_LIST}"

list_param:
	@echo "MATRIX_DB_FILE	${MATRIX_DB_FILE}"
	@echo "MATRIX_DIR	${MATRIX_DIR}"
	@echo "SPLIT_DIR	${SPLIT_DIR}"
	@echo "MATRIX_ID	${MATRIX_ID}"
	@echo "MATRIX_FILE	${MATRIX_FILE}"
	@echo "BG_OL		${BG_OL}"
	@echo "BG_PREFIX	${BG_PREFIX}"
	@echo "BG_FILE		${BG_FILE}"
	@echo "DISTRIB_DIR	${DISTRIB_DIR}"
	@echo "DISTRIB_FILE	${DISTRIB_FILE}"

BG_OL=3
ORG=Homo_sapiens_EnsEMBL
BG_PREFIX=${BG_OL}nt_upstream-noorf_${ORG}-ovlp-1str
BG_FILE=${RSAT}/data/genomes/${ORG}/oligo-frequencies/${BG_PREFIX}.freq
MATRIX_ID=`head -1 ${MATRIX_LIST} | cut -f 2`
MATRIX_FILE=${SPLIT_DIR}/${DB}_${MATRIX_ID}.tf
DISTRIB_DIR=${MATRIX_DIR}/pval_distributions_${DB}
DISTRIB_FILE=${DISTRIB_DIR}/${DB}_${MATRIX_ID}_distrib_bg_${BG_PREFIX}
one_matrix_distrib:
	@echo
	@echo "${DATE}	Computing distribution for matrix	${MATRIX_ID}"
	mkdir -p ${DISTRIB_DIR}
	matrix-distrib -v ${V} -m ${MATRIX_FILE}  -matrix_format tf -decimals 1 \
		-bgfile ${BG_FILE} -bg_pseudo 0.01 -bg_format oligos -pseudo 1 \
		-o ${DISTRIB_FILE}.tab
	@echo "	${DISTRIB_FILE}.tab"

MATRIX_IDS=`grep -v '^;' ${MATRIX_LIST} | cut -f 2 | xargs`
list_matrix_ids:
	@echo "Matrix IDs	${MATRIX_IDS}"

matrix_distrib_all:
	@for m in ${MATRIX_IDS}; do \
		${MAKE} one_matrix_distrib MATRIX_ID=$${m} ; \
	done

TIME_FILE=${MATRIX_DIR}/distrib_calc_time.txt
matrix_distrib_all_time: 
	time ${MAKE} matrix_distrib_all >& ${TIME_FILE}
