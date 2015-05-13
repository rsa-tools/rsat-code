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
## Compute cross-comparisons between a set of motifs discovered using
## peak-motifs.
#PEAKMO_PREFIX=peak-motifs_result_Chen_Oct4
#PEAKMO_MATRICES=${RSAT}/public_html/demo_files/${PEAKMO_PREFIX}_matrices.tf
PEAKMO_VS_PEAKMO_DIR=results/peakmo_vs_peakmo
PEAKMO_VS_PEAKMO=${PEAKMO_VS_PEAKMO_DIR}/${PEAKMO_PREFIX}__vs-itself_${COMPA_SUFFIX}
#SCORES=cor,Ncor,NcorS,logoDP,NIcor,NsEucl,SSD,NSW,match_rank,zscores
PLOT_FORMAT=pdf
peakmo_vs_peakmo:
	@mkdir -p ${PEAKMO_VS_PEAKMO_DIR}
	@echo 
	@echo "Comparing moti collection against itself, using z-scores"
	compare-matrices -v ${V} \
		-mode matches \
		-format1 transfac -file1 ${PEAKMO_MATRICES} \
		-format2 transfac -file2 ${PEAKMO_MATRICES} \
		-mode scores \
		-DR \
		-lth w ${MIN_W} \
		-return matrix_name,matrix_id \
		-return ${SCORES},zscores \
		-return width,strand,offset,consensus \
		-sort Ncor ${OPT} \
		-o ${PEAKMO_VS_PEAKMO}
	@echo "	${PEAKMO_VS_PEAKMO}"
	@echo
	@XYgraph -i  results/peakmo_vs_peakmo/peak-motifs_result_Chen_Oct4__vs-itself_w5_wr0.3_cor0.75_Ncor0.4.tab \
		-xcol 42 -xleg1 "match rank" \
		-ycol 32 -yleg1 "mean z-score" \
		-lines -hline 'red' 0 -r_plot \
		-format ${PLOT_FORMAT} \
		-o ${PEAKMO_VS_PEAKMO}_rank_vs_zscore.${PLOT_FORMAT}
	@echo "	${PEAKMO_VS_PEAKMO}_rank_vs_zscore.${PLOT_FORMAT}"

################################################################
## Case 1: peak-motifs result versus JASPAR database
PEAKMO_PREFIX=peak-motifs_result_Chen_Oct4
PEAKMO_MATRICES=${RSAT}/public_html/demo_files/${PEAKMO_PREFIX}_matrices.tf
JASPAR_VERSION=2015_03
JASPAR_PREFIX=jaspar_core_vertebrates_${JASPAR_VERSION}
JASPAR_MATRICES=${RSAT}/public_html/motif_databases/JASPAR/${JASPAR_PREFIX}.tf
PEAKMO_VS_JASPAR_DIR=results/peakmo_vs_jaspar
PEAKMO_VS_JASPAR=${PEAKMO_VS_JASPAR_DIR}/${PEAKMO_PREFIX}__vs__${JASPAR_PREFIX}_${COMPA_SUFFIX}
SCORES=cor,Ncor,NcorS,logoDP,NIcor,NsEucl,SSD,NSW,match_rank
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
		-return ${SCORES} \
		-return width,strand,offset,consensus,alignments_1ton \
		-sort Ncor ${OPT} \
		-o ${PEAKMO_VS_JASPAR}
	@echo ${PEAKMO_VS_JASPAR}


peakmo_vs_jaspar_param:
	@echo "PEAKMO_MATRICES	${PEAKMO_MATRICES}"
	@echo "JASPAR_MATRICES	${JASPAR_MATRICES}"
	@echo "Result	${PEAKMO_VS_JASPAR}_index.html"

################################################################
## Use compare-matrices-quick (less options than compare-matrices, but
## 100 times faster)
COMPA_SUFFIX_QUICK=w${MIN_W}_Ncor${MIN_NCOR}
PEAKMO_VS_JASPAR_QUICK_DIR=results/peakmo_vs_jaspar_quick
PEAKMO_VS_JASPAR_QUICK=${PEAKMO_VS_JASPAR_QUICK_DIR}/${PEAKMO_PREFIX}__vs__${JASPAR_PREFIX}_${COMPA_SUFFIX_QUICK}_quick
peakmo_vs_jaspar_quick:
	@mkdir -p ${PEAKMO_VS_JASPAR_QUICK_DIR}
	compare-matrices-quick -v ${V} \
	-file1 ${PEAKMO_MATRICES} \
	-file2 ${JASPAR_MATRICES} \
	-lth_w ${MIN_W} \
	-lth_ncor ${MIN_NCOR} \
	-o ${PEAKMO_VS_JASPAR_QUICK}
	@echo ${PEAKMO_VS_JASPAR_QUICK}


################################################################
## Case 2: compare all the motifs from a reference database
## (RegulonDB, Jaspar).
DB_PREFIX=${JASPAR_PREFIX}
DB_DIR=${JASPAR_DIR}
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
TIME_DIR=${DB_COMPA_DIR}
TIME_FILE=${DB_COMPA_DIR}/time_${DB_PREFIX}_vs_itself.txt
db_vs_itself:
	@echo ""
	@echo "Comparing DB with itself	${DB_PREFIX}"
	@mkdir -p ${DB_COMPA_DIR}
	@mkdir -p ${TIME_DIR}
	@echo "Time and log file	${TIME_FILE}"
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

################################################################
## Generate a graph of motif similarity (nodes = motifs, edges =
## similarity between two motifs)
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
	${MAKE} db_vs_itself_layout

## Choose a graph layout algorithm among those supported by GraphViz
## (dot, neato, fdp, sfdp, twopi, circo) The layered out graph is
## exported in dot format (which can be used to extract node
## coordinates) and to pdf (for visualisation).
LAYOUT=neato
db_vs_itself_layout:
	@${LAYOUT} -Tdot ${DB_COMPA_RESULT}.dot ${OPT} > ${DB_COMPA_RESULT}_${LAYOUT}.dot
	@echo "	${DB_COMPA_RESULT}_${LAYOUT}.dot"
	@${LAYOUT} -Tpdf ${DB_COMPA_RESULT}.dot ${OPT} > ${DB_COMPA_RESULT}_${LAYOUT}.pdf
	@echo "	${DB_COMPA_RESULT}_${LAYOUT}.pdf"


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

################################################################
## RegulonDB
REGULONDB_PREFIX=regulonDB_2014-04-11
REGULONDB_DIR=${RSAT}/public_html/motif_databases/REGULONDB
REGULONDB_MATRICES=${REGULONDB_DIR}/${REGULONDB_PREFIX}.tf
regulondb_vs_itself:
	@${MAKE} db_vs_itself DB_PREFIX=${REGULONDB_PREFIX} DB_DIR=${REGULONDB_DIR}

permute_regulondb:
	@${MAKE} permute_db  DB_PREFIX=${REGULONDB_PREFIX} DB_DIR=${REGULONDB_DIR}

permuted_regulondb_vs_itself:
	@${MAKE} db_vs_itself DB_PREFIX=${REGULONDB_PREFIX}_perm DB_DIR=${REGULONDB_DIR}

################################################################
## footprintDB (database compiled by Bruno Contreras)
FOOTPRINTDB_PREFIX=footprintDB.plants.motif
FOOTPRINTDB_DIR=${RSAT}/public_html/motif_databases/footprintDB
FOOTPRINTDB_MATRICES=${FOOTPRINTDB_DIR}/${FOOTPRINTDB_PREFIX}.tf
footprintdb_vs_itself:
	@${MAKE} db_vs_itself DB_PREFIX=${FOOTPRINTDB_PREFIX} DB_DIR=${FOOTPRINTDB_DIR}

permute_footprintdb:
	@${MAKE} permute_db  DB_PREFIX=${FOOTPRINTDB_PREFIX} DB_DIR=${FOOTPRINTDB_DIR}

permuted_footprintdb_vs_itself:
	@${MAKE} db_vs_itself DB_PREFIX=${FOOTPRINTDB_PREFIX}_perm DB_DIR=${FOOTPRINTDB_DIR}

################################################################
## JASPAR
JASPAR_GROUPS=all insects vertebrates nematodes fungi urochordates plants
JASPAR_GROUP=vertebrates
JASPAR_PREFIX=jaspar_core_${JASPAR_GROUP}_${JASPAR_VERSION}
JASPAR_DIR=${RSAT}/public_html/motif_databases/JASPAR
JASPAR_MATRICES=${JASPAR_DIR}/${JASPAR_PREFIX}.tf
jaspar_one_group_vs_itself:
	${MAKE} db_vs_itself DB_PREFIX=${JASPAR_PREFIX} DB_DIR=${JASPAR_DIR}

permute_jaspar_one_group:
	@${MAKE} permute_db DB_PREFIX=${JASPAR_PREFIX} DB_DIR=${JASPAR_DIR}

permuted_jaspar_one_group_vs_itself:
	@${MAKE} db_vs_itself DB_PREFIX=${JASPAR_PREFIX}_perm DB_DIR=${RSAT}/public_html/motif_databases/JASPAR

JASPAR_TASK=jaspar_one_group_vs_itself permute_jaspar_one_group permuted_jaspar_one_group_vs_itself
iterate_jaspar:
	@for g in ${JASPAR_GROUPS}; do \
		${MAKE} DB_PREFIX=jaspar_core_$${g}_${JASPAR_VERSION} DB_DIR=${RSAT}/public_html/motif_databases/JASPAR JASPAR_GROUP=$$g ${JASPAR_TASK}; \
	done

iterate_jaspar_perm:
	@for g in ${JASPAR_GROUPS}; do \
		${MAKE} DB_PREFIX=jaspar_core_$${g}_${JASPAR_VERSION}_perm DB_DIR=${RSAT}/public_html/motif_databases/JASPAR JASPAR_GROUP=$$g ${JASPAR_TASK}; \
	done



