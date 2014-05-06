################################################################
## Demonstration for compare-matrices, based on selected study cases

include ${RSAT}/makefiles/util.mk
MAKEFILE=${RSAT}/makefiles/compare-matrices_demo.mk
V=2

################################################################
## Case 1: peak-motifs result versus JASPAR database
PEAKMO_PREFIX=peak-motifs_result_Chen_Oct4
PEAKMO_MATRICES=${RSAT}/public_html/demo_files/${PEAKMO_PREFIX}_matrices.tf

JASPAR_PREFIX=jaspar_core_vertebrates_2013-11
JASPAR_MATRICES=${RSAT}/public_html/data/motif_databases/JASPAR/${JASPAR_PREFIX}.tf

PEAKMO_VS_JASPAR_DIR=results/peakmo_vs_jaspar
PEAKMO_VS_JASPAR=${PEAKMO_VS_JASPAR_DIR}/${PEAKMO_PREFIX}__vs__${JASPAR_PREFIX}
peakmo_vs_jaspar:
	@mkdir -p ${PEAKMO_VS_JASPAR_DIR}
	compare-matrices -v ${V} \
		-mode matches \
		-format1 transfac -file1 ${PEAKMO_MATRICES} \
		-format2 transfac -file2 ${JASPAR_MATRICES} \
		-mode matches \
		-DR \
		-uth offset_rank 1 \
		-lth w 5 \
		-lth Wr 0.3 \
		-lth cor 0.75 \
		-lth Ncor 0.4 \
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
DB_COMPA_RESULT=${DB_COMPA_DIR}/${DB_PREFIX}_vs_itself
db_vs_itself:
	@mkdir -p ${DB_COMPA_DIR}
	compare-matrices -v ${V} \
		-mode matches \
		-format transfac -file ${DB_MATRICES} \
		-mode matches \
		-DR \
		-uth offset_rank 1 \
		-lth w 5 \
		-lth Wr 0.3 \
		-lth cor 0.75 \
		-lth Ncor 0.4 \
		-return matrix_name,matrix_id \
		-return cor,Ncor,logoDP,NIcor,NsEucl,SSD,NSW,match_rank \
		-return width,strand,offset,consensus,alignments_1ton \
		-sort Ncor \
		-o ${DB_COMPA_RESULT}
	@echo ${DB_COMPA_RESULT}_index.html

db_vs_itself_param:
	@echo "DB_PREFIX	${DB_PREFIX}"
	@echo "DB_DIR		${DB_DIR}"
	@echo "DB_MATRICES	${DB_MATRICES}"
	@echo "DB_COMPA_DIR	${DB_COMPA_DIR}"
	@echo "DB_COMPA_RESULT	${DB_COMPA_RESULT}"


REGULONDB_PREFIX=regulonDB_2012-05
REGULONDB_DIR=${RSAT}/public_html/data/motif_databases/REGULONDB
REGULONDB_MATRICES=${REGULONDB_DIR}/${REGULONDB_PREFIX}.tf
regulondb_vs_itself:
	@${MAKE} db_vs_itself DB_PREFIX=regulonDB_2012-05 DB_DIR=${RSAT}/public_html/data/motif_databases/REGULONDB

