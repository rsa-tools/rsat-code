################################################################
## Demo for the RSAT tool matrix-clustering
##
## Authors: Jaime Castro & Jacques van Helden
## Date: Jan-April 2014


#include ${RSAT}/makefiles/util.mk
## We include compare-matrices_demo.mk to get the parameters (matrix files)
include ${RSAT}/makefiles/compare-matrices_demo.mk
MAKEFILE=${RSAT}/makefiles/matrix-clustering_demo2.mk

################################################################
## Parameters for the analysis
MIN_NCOR=0.4
MIN_COR=0.75
HCLUST_METHOD=average
#MIN_W=4
## Verbosity
V=2

## Define a set of demo files
PEAKMO_PREFIX=peak-motifs_result_Chen_Oct4
FOOTPRINT_DISCO_PREFIX=footprint-discovery_LexA
PEAKMO_NEG_CONTROL_PREFIX=peak-motifs_result_Chen_Oct4_permuted
OCT4_PREFIX=peak-motifs_Oct4


## Choose a particular demo set
DEMO_PREFIX=${PEAKMO_PREFIX}

## Define file locations based on the chosen demo set
MATRIX_DIR=${RSAT}/public_html/demo_files/
MATRIX_FILE=${MATRIX_DIR}/${DEMO_PREFIX}_matrices.tf


################################################################
## Compare all matrices from the input file, with specific parameters
## to ensure that all distances are computed.
COMPA_DIR=results/${DEMO_PREFIX}
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
CLUSTER_PREFIX=${DEMO_PREFIX}_hclust${HCLUST_METHOD}_Ncor${MIN_NCOR}_cor${MIN_COR}
CLUSTER_DIR=results/${DEMO_PREFIX}/hclust${HCLUST_METHOD}/Ncor${MIN_NCOR}_cor${MIN_COR}
CLUSTER_FILE_PREFIX=${CLUSTER_DIR}/${CLUSTER_PREFIX}
CLUSTER_CMD=matrix-clustering -v ${V} \
		-i ${MATRIX_FILE} -format tf \
		-lth Ncor ${MIN_NCOR} \
		-lth cor ${MIN_COR} \
		-cons \
		-export newick -d3_base file -hclust_method ${HCLUST_METHOD} \
		-label name,consensus ${OPT} \
		-o ${CLUSTER_FILE_PREFIX}
CLUSTER_TIME_FILE=${CLUSTER_FILE_PREFIX}_time_log.txt
cluster:
	@echo
	@echo "Running matrix-clustering	${DEMO_PREFIX}"
	@echo "	verbosity +time in file	${CLUSTER_TIME}"
	${CLUSTER_CMD}
#	(time ${CLUSTER_CMD}) >& ${CLUSTER_TIME_FILE}
	@echo "		${CLUSTER_FILE_PREFIX}_index.html"

## Cluster motifs resulting from peak-motifs (Chen Oct4 data set)
cluster_peakmo_no_threshold:
	@echo
	@echo "Running matrix-clustering on motifs discovered by peak-motifs (Oct 4 dataset from Chen 2008)"
	${MAKE} cluster DEMO_PREFIX=${PEAKMO_PREFIX} MIN_NCOR=-1 MIN_COR=-1

## Cluster motifs resulting from peak-motifs (Chen Oct4 data set)
cluster_peakmo_threhsolds:
	@echo
	@echo "Running matrix-clustering on motifs discovered by peak-motifs (Oct 4 dataset from Chen 2008)"
	${MAKE} cluster DEMO_PREFIX=${PEAKMO_PREFIX}

## Cluster permuted motifs resulting from peak-motifs (Chen Oct4 data set)
cluster_peakmo_neg_control:
	@echo
	@echo "Running matrix-clustering on permuted motifs discovered by peak-motifs (Oct 4 dataset from Chen 2008)"
	${MAKE} cluster DEMO_PREFIX=${PEAKMO_NEG_CONTROL_PREFIX}

## Cluster motifs resulting from peak-motifs (Chen Oct4 data set)
cluster_peakmo_Oct4_threhsolds:
	@echo
	@echo "Running matrix-clustering on motifs discovered by peak-motifs (Oct 4 dataset from Chen 2008)"
	${MAKE} cluster DEMO_PREFIX=${OCT4_PREFIX} MIN_NCOR=0.4 MIN_COR=0.7 
## We should add this option: OPT='-lth w 5'

## Cluster motifs resulting from footprint-discovery (LexA in Enterobacteriales)
cluster_footprints:
	@echo
	@echo "Running matrix-clustering on motifs discovered by footprint-discovery (query gene=LexA; taxon=Enterobacteriales)"
	${MAKE} cluster DEMO_PREFIX=${FOOTPRINT_DISCO_PREFIX}



## Cluster all motifs from RegulonDB
RDB_CLUSTER_DIR=results/regulondDB_clusters
RDB_CLUSTERS=${RDB_CLUSTER_DIR}/RDB_clusters
#RDB_PREFIX=regulonDB_2012-05
RDB_PREFIX=regulonDB_2014-04-11
RDB_MATRICES=${RSAT}/data/motif_databases/REGULONDB/${RDB_PREFIX}.tf
#RDB_PREFIX=regulonDB_2012-05_MOD2
#RDB_MATRICES=${RSAT}/data/motif_databases/REGULONDB/regulonDB_2012-05_MOD2.tf
cluster_rdb:
	@echo "Clustering all matrices from RegulonDB"
	${MAKE} cluster DEMO_PREFIX=${RDB_PREFIX} MATRIX_FILE=${RDB_MATRICES} MIN_NCOR=0.4


cluster_jaspar_one_group:
	@echo "Clustering all matrices from JASPAR ${JASPAR_GROUP}"
	${MAKE} cluster DEMO_PREFIX=${JASPAR_PREFIX} MATRIX_FILE=${JASPAR_MATRICES}
