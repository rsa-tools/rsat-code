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
		-sort Ncor ${OPT} \
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
		-sort Ncor ${OPT} \
		-o ${DB_COMPA_RESULT}
#TIME_FILE=${DB_COMPA_DIR}/time_${DB_PREFIX}_vs_itself.txt
TIME_FILE=time_measurements/time_${DB_PREFIX}_vs_itself.txt
db_vs_itself:
	@echo ""
	@echo "Comparing DB with itself	${DB_PREFIX}"
	@mkdir -p ${DB_COMPA_DIR}
	(time -p ${DB_COMPA_CMD}) >& ${TIME_FILE}
	@echo "	${DB_COMPA_RESULT}_index.html"
	@echo "	${TIME_FILE}"
	@${MAKE} db_vs_itself_graph
	@${MAKE} graph_stats

STATS_FILE=${DB_COMPA_DIR}/${DB_PREFIX}_vs_itself_stats.tab
DB_MATRIX_NB=`grep '^AC' ${DB_MATRICES} | wc -l | perl -pe 's|^\s+||'`
DB_COMPA_EDGE_NB=`grep -v '^;' ${DB_COMPA_RESULT}.tab | grep -v '^\#' | wc -l | perl -pe 's|^\s+||'`
DB_COMPA_TIME=`grep '^user' ${TIME_FILE} | perl -pe 's|^user\s+||'`
graph_stats:
	@echo
	@echo "Collecting stats for ${DB_PREFIX}"
	@echo "DB_PREFIX	${DB_PREFIX}" > ${STATS_FILE}
	@echo "matrices  	${DB_MATRIX_NB}" >> ${STATS_FILE}
	@echo "graph_edges  	${DB_COMPA_EDGE_NB}" >> ${STATS_FILE}
	@echo "compa_time  	${DB_COMPA_TIME}" >> ${STATS_FILE}
	@echo "	${STATS_FILE}"

## Generate a graph of motif similarity (nodes = motifs, edges = similarity between two motifs)
WCOL=6
db_vs_itself_graph:
	@echo ""
	@echo "Generating motif similarity graph"
	perl -pe 's|\.\dnt\S+||g' ${DB_COMPA_RESULT}.tab \
		| convert-graph  -from tab -to gml \
		-scol 3 -tcol 4 -wcol ${WCOL} -undirected -ewidth -ecolors fire -min -1 -max 1 \
		> ${DB_COMPA_RESULT}.gml
	@echo "	${DB_COMPA_RESULT}.gml"
	@perl -pe 's|\.\dnt\S+||g' ${DB_COMPA_RESULT}.tab \
		| convert-graph -from tab -to dot \
		-scol 3 -tcol 4 -wcol ${WCOL} -undirected -ewidth -ecolors fire -min ${MIN_NCOR} -max 1 \
		| perl -pe 's|\.\dnt\S+||g' \
		> ${DB_COMPA_RESULT}.dot
	@echo "	${DB_COMPA_RESULT}.dot"
	@neato -Tdot ${DB_COMPA_RESULT}.dot > ${DB_COMPA_RESULT}_neato.dot
	@echo "	${DB_COMPA_RESULT}_neato.dot"
	@neato -Tpdf ${DB_COMPA_RESULT}.dot > ${DB_COMPA_RESULT}.pdf
	@echo "	${DB_COMPA_RESULT}.pdf"


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

permuted_regulondb_vs_itself:
	@${MAKE} db_vs_itself DB_PREFIX=${REGULONDB_PREFIX}_perm DB_DIR=${RSAT}/public_html/data/motif_databases/REGULONDB


## JASPAR core insects
JASPAR_GROUPS=all insects vertebrates nematodes fungi urochordates plants
JASPAR_GROUP=vertebrates
JASPAR_PREFIX=jaspar_core_${JASPAR_GROUP}_2013-11
JASPAR_DIR=${RSAT}/public_html/data/motif_databases/JASPAR
JASPAR_MATRICES=${JASPAR_DIR}/${JASPAR_PREFIX}.tf
jaspar_one_group_vs_itself:
	${MAKE} db_vs_itself DB_PREFIX=${JASPAR_PREFIX} DB_DIR=${JASPAR_DIR}

permute_jaspar_one_group:
	@${MAKE} permute_db DB_PREFIX=${JASPAR_PREFIX} DB_DIR=${JASPAR_DIR}

permuted_jaspar_one_group_vs_itself:
	@${MAKE} db_vs_itself DB_PREFIX=${JASPAR_PREFIX}_perm DB_DIR=${RSAT}/public_html/data/motif_databases/JASPAR

JASPAR_TASK=jaspar_one_group_vs_itself permute_jaspar_one_group permuted_jaspar_one_group_vs_itself
iterate_jaspar:
	@for g in ${JASPAR_GROUPS}; do \
		${MAKE} DB_PREFIX=jaspar_core_$${g}_2013-11 DB_DIR=${RSAT}/public_html/data/motif_databases/JASPAR JASPAR_GROUP=$$g ${JASPAR_TASK}; \
	done

iterate_jaspar_perm:
	@for g in ${JASPAR_GROUPS}; do \
		${MAKE} DB_PREFIX=jaspar_core_$${g}_2013-11_perm DB_DIR=${RSAT}/public_html/data/motif_databases/JASPAR JASPAR_GROUP=$$g ${JASPAR_TASK}; \
	done
