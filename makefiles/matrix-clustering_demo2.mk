################################################################
## Demo for the RSAT tool matrix-clustering
##
## Authors: Jaime Castro & Jacques van Helden
## Date: Jan-April 2014


#include ${RSAT}/makefiles/util.mk
## We include compare-matrices_demo.mk to get the parameters (matrix files)
include ${RSAT}/makefiles/compare-matrices_demo.mk
MAKEFILE=${RSAT}/makefiles/matrix-clustering_demo2.mk

## Define a set of demo files
PEAKMO_PREFIX=peak-motifs_result_Chen_Oct4
FOOTPRINT_DISCO_PREFIX=footprint-discovery_LexA

## Choose a particular demo set
DEMO_PREFIX=${PEAKMO_PREFIX}

## Define file locations based on the chosen demo set
MATRIX_DIR=${RSAT}/public_html/demo_files/
MATRIX_FILE=${MATRIX_DIR}/${DEMO_PREFIX}_matrices.tf

## Verbosity
V=2

################################################################
## Compare all matrices from the input file, with specific parameters
## to ensure that all distances are computed.
COMPA_DIR=results/${DEMO_PREFIX}/pairwise_comparisons
COMPA_FILE=${COMPA_DIR}/${DEMO_PREFIX}_compa.tab
compa:
	@echo
	@echo "Motif comparisons"
	@mkdir -p ${COMPA_DIR}
	compare-matrices  -v ${V} \
		-file ${MATRIX_FILE} -format tf \
		-lth w 1 -lth cor -1 -lth Ncor -1 \
		-return Ncor,strand,offset,width,consensus \
		-o ${COMPA_FILE}
	@echo "	${COMPA_FILE}"

_cluster_old:
	cat cluster_motifs.R | \
		R  --slave --no-save --no-restore --no-environ \
		--args "infile='${COMPA_FILE}';outfile='results/peak-motifs_7nt_merged_oligos_positions_compa.json';score='Ncor'" \
		> cluster_log.txt

compa_peak_motifs:
	${MAKE} compa DEMO_PREFIX=${PEAKMO_PREFIX}

compa_footprint_discovery:
	${MAKE} compa DEMO_PREFIX=${FOOTPRINT_DISCO_PREFIX}

################################################################
## Run matrix-clusteringon one demo set (the particular cases will be
## specified below)
MIN_NCOR=0.3
CLUSTER_DIR=results/${DEMO_PREFIX}/motif_clusters_Ncor${MIN_NCOR}
CLUSTER_PREFIX=${COMPA_DIR}/${DEMO_PREFIX}_clustering
CLUSTER_CMD=matrix-clustering -v ${V} \
		-i ${MATRIX_FILE} -format tf \
		-lth Ncor ${MIN_NCOR} \
		-export newick -d3_base link -hclust_method average \
		-labels name,consensus ${OPT} \
		-o ${CLUSTER_PREFIX}
CLUSTER_TIME_FILE=${CLUSTER_PREFIX}_time_log.txt
cluster:
	@echo
	@echo "Running matrix-clustering	${DEMO_PREFIX}"
	@echo "	verbosity +time in file	${CLUSTER_TIME}"
	${CLUSTER_CMD}
#	(time ${CLUSTER_CMD}) >& ${CLUSTER_TIME_FILE}
	@echo "		${CLUSTER_PREFIX}_index.html"

## Cluster motifs resulting from peak-motifs (Chen Oct4 data set)
cluster_peakmo_no_threshold:
	@echo
	@echo "Running matrix-clustering on motifs discovered by peak-motifs (Oct 4 dataset from Chen 2008)"
	${MAKE} cluster DEMO_PREFIX=${PEAKMO_PREFIX} MIN_NCOR=0

## Cluster motifs resulting from peak-motifs (Chen Oct4 data set)
cluster_peakmo_threhsolds:
	@echo
	@echo "Running matrix-clustering on motifs discovered by peak-motifs (Oct 4 dataset from Chen 2008)"
	${MAKE} cluster DEMO_PREFIX=${PEAKMO_PREFIX} MIN_NCOR=0.3

## Cluster motifs resulting from footprint-discovery (LexA in Enterobacteriales)
cluster_footprints:
	@echo
	@echo "Running matrix-clustering on motifs discovered by footprint-discovery (query gene=LexA; taxon=Enterobacteriales)"
	${MAKE} cluster DEMO_PREFIX=${FOOTPRINT_DISCO_PREFIX}

## Cluster all motifs from RegulonDB
RDB_CLUSTER_DIR=results/regulondDB_clusters
RDB_CLUSTERS=${RDB_CLUSTER_DIR}/RDB_clusters
RDB_PREFIX=regulonDB_2012-05
RDB_MATRICES=${RSAT}/data/motif_databases/REGULONDB/${RDB_PREFIX}.tf
cluster_rdb:
	@echo "Clustering all matrices from RegulonDB"
	${MAKE} cluster DEMO_PREFIX=${RDB_PREFIX} MATRIX_FILE=${RDB_MATRICES}


cluster_jaspar_one_group:
	@echo "Clustering all matrices from JASPAR ${JASPAR_GROUP}"
	${MAKE} cluster DEMO_PREFIX=${JASPAR_PREFIX} MATRIX_FILE=${JASPAR_MATRICES}
