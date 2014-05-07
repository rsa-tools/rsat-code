################################################################
## Demonstration for compare-matrices, based on selected study cases

include ${RSAT}/makefiles/util.mk
MAKEFILE=${RSAT}/makefiles/compare-matrices_demo.mk
V=2

## Default parameters for comparison
MIN_W=5
MIN_WR=0.3
MIN_COR=0.75
MIN_NCOR=0.4

COMPA_SUFFIX=w${MIN_W}_wr${MIN_WR}_cor${MIN_COR}_Ncor${MIN_NCOR}

################################################################
## Case 1: peak-motifs result versus JASPAR database
PEAKMO_PREFIX=peak-motifs_result_Chen_Oct4
PEAKMO_MATRICES=${RSAT}/public_html/demo_files/${PEAKMO_PREFIX}_matrices.tf

JASPAR_PREFIX=jaspar_core_vertebrates_2013-11
JASPAR_MATRICES=${RSAT}/public_html/data/motif_databases/JASPAR/${JASPAR_PREFIX}.tf

PEAKMO_VS_JASPAR_DIR=results/peakmo_vs_jaspar
PEAKMO_VS_JASPAR=${PEAKMO_VS_JASPAR_DIR}/${PEAKMO_PREFIX}__vs__${JASPAR_PREFIX}_${COMPA_SUFFIX}
peakmo_vs_jaspar:
	@mkdir -p ${PEAKMO_VS_JASPAR_DIR}
	compare-matrices -v ${V} \
		-mode matches \
		-format1 transfac -file1 ${PEAKMO_MATRICES} \
		-format2 transfac -file2 ${JASPAR_MATRICES} \
		-mode matches \
		-DR \
		-uth offset_rank 1 \
		-lth w ${MIN_W} \
		-lth Wr ${MIN_WR} \
		-lth cor ${MIN_COR} \
		-lth Ncor ${MIN_NCOR} \
		-return matrix_name,matrix_id \
		-return cor,Ncor,logoDP,NIcor,NsEucl,SSD,NSW,match_rank \
		-return width,strand,offset,consensus,alignments_1ton \
		-sort Ncor \
		-o ${PEAKMO_VS_JASPAR}
	@echo ${PEAKMO_VS_JASPAR}

peakmo_vs_jaspar_param:
	@echo "PEAKMO_MATRICES	${PEAKMO_MATRICES}"
	@echo "JASPAR_MATRICES	${JASPAR_MATRICES}"
	@echo "Result	${PEAKMO_VS_JASPAR}_index.html"

################################################################
## Case 2: 
DB_PREFIX=${REGULONDB_PREFIX}
DB_DIR=${REGULONDB_DIR}
DB_MATRICES=${DB_DIR}/${DB_PREFIX}.tf
DB_COMPA_DIR=results/${DB_PREFIX}_vs_itself
DB_COMPA_RESULT=${DB_COMPA_DIR}/${DB_PREFIX}_vs_itself_${COMPA_SUFFIX}
DB_COMPA_CMD=compare-matrices -v ${V} \
		-mode matches -distinct \
		-format transfac -file ${DB_MATRICES} \
		-mode matches \
		-DR \
		-uth offset_rank 1 \
		-lth w ${MIN_W} \
		-lth Wr ${MIN_WR} \
		-lth cor ${MIN_COR} \
		-lth Ncor ${MIN_NCOR} \
		-return matrix_name,matrix_id \
		-return cor,Ncor,logoDP,NIcor,NsEucl,SSD,NSW,match_rank \
		-return width,strand,offset,consensus \
		-sort Ncor \
		-o ${DB_COMPA_RESULT}
TIME_FILE=time_${DB_PREFIX}_vs_itself.txt
DB_COMPA_EDGE_NB=`grep -v '^;' ${DB_COMPA_RESULT}.tab | grep -v '^\#' | wc -l`
db_vs_itself:
	@echo ""
	@echo "Comparing DB with itself	${DB_PREFIX}"
	@mkdir -p ${DB_COMPA_DIR}
	(time ${DB_COMPA_CMD}) >& time_${DB_PREFIX}_vs_itself.txt
	@echo "	${DB_COMPA_RESULT}_index.html"
	@echo "	${TIME_FILE}"
	@echo "${DB_PREFIX}_vs_itself	${DB_COMPA_EDGE_NB}	edges" > ${DB_PREFIX}_vs_itself_edges.txt
	@echo ${DB_PREFIX}_vs_itself_edges.txt

## Generate a graph of motif similarity (nodes = motifs, edges = similarity between two motifs)
db_vs_itself_graph:

## Display parameters for the matrix comparison
db_vs_itself_param:
	@echo "DB_PREFIX		${DB_PREFIX}"
	@echo "DB_DIR			${DB_DIR}"
	@echo "DB_MATRICES		${DB_MATRICES}"
	@echo "DB_COMPA_DIR		${DB_COMPA_DIR}"
	@echo "DB_COMPA_RESULT (tab)	${DB_COMPA_RESULT}.tab"
	@echo "DB_COMPA_RESULT (html)	${DB_COMPA_RESULT}.html"
	@echo "DB_COMPA_RESULT (index)	${DB_COMPA_RESULT}_index.html"
	@echo "NB OF LINKS		${DB_COMPA_EDGE_NB}"

## Permute all matrices of a database (for negative control)
DB_MATRICES_PERM=${DB_DIR}/${DB_PREFIX}_perm.tf
permute_db:
	@echo
	@echo "Permuting columns for ${DB_PREFIX}"
	permute-matrix -i ${DB_MATRICES} -in_format tf -out_format tf -o ${DB_MATRICES_PERM}
	@echo "	${DB_MATRICES_PERM}"

## RegulonDB
REGULONDB_PREFIX=regulonDB_2012-05
REGULONDB_DIR=${RSAT}/public_html/data/motif_databases/REGULONDB
REGULONDB_MATRICES=${REGULONDB_DIR}/${REGULONDB_PREFIX}.tf
regulondb_vs_itself:
	@${MAKE} db_vs_itself DB_PREFIX=${REGULONDB_PREFIX} DB_DIR=${REGULONDB_DIR}

permute_regulondb:
	@${MAKE} permute_db  DB_PREFIX=${REGULONDB_PREFIX} DB_DIR=${REGULONDB_DIR}

regulondb_vs_permuted:
	@${MAKE} db_vs_itself DB_PREFIX=${REGULONDB_PREFIX}_perm DB_DIR=${RSAT}/public_html/data/motif_databases/REGULONDB


## JASPAR core insects
JASPAR_GROUPS=all insects vertebrates nematods fungi urochordates plants
JASPAR_GROUP=insects
JASPAR_PREFIX=jaspar_core_${JASPAR_GROUP}_2013-11
JASPAR_DIR=${RSAT}/public_html/data/motif_databases/JASPAR
JASPAR_MATRICES=${JASPAR_DIR}/${JASPAR_PREFIX}.tf
jaspar_one_group_vs_itself:
	${MAKE} db_vs_itself DB_PREFIX=${JASPAR_PREFIX} DB_DIR=${JASPAR_DIR}

permute_jaspar_one_group:
	@${MAKE} permute_db  DB_PREFIX=${JASPAR_PREFIX} DB_DIR=${JASPAR_DIR}

jaspar_one_group_vs_permuted:
	@${MAKE} db_vs_itself DB_PREFIX=${JASPAR_PREFIX}_perm DB_DIR=${RSAT}/public_html/data/motif_databases/JASPAR

jaspar:
	@for g in ${JASPAR_GROUPS}; do \
		${MAKE} JASPAR_GROUP=$$g jaspar_one_group_vs_itself; \
		${MAKE} JASPAR_GROUP=$$g permute_jaspar_one_group; \
		${MAKE} JASPAR_GROUP=$$g jaspar_one_group_vs_permuted; \
	done
