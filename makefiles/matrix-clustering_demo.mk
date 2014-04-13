################################################################
## Demo for the RSAT tool matrix-clustering
##
## Authors: Jaime Castro & Jacques van Helden
## Date: Jan-April 2014


include ${RSAT}/makefiles/util.mk
MAKEFILE=${RSAT}/makefiles/matrix-clustering_demo.mk

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
COMPA_FILE=${COMPA_DIR}/peak-motifs_7nt_merged_oligos_positions_compa.tab
compa:
	@echo
	@echo "Motif comparisons"
	@mkdir -p ${COMPA_DIR}
	compare-matrices  -v ${V} \
		-file ${MATRIX_FILE} -format tf \
		-lth w 1 -lth cor -1 -lth Ncor -1 \
		-return Ncor,strand,offset,width,consensus,alignments_1ton \
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
## Demo 1: clustering between motifs discovered by peak-motifs
cluster_peakmo:
	@echo "Clustering of motifs discovered by peak-motifs"
	matrix-clustering -v ${V} \
		-i ${PEAKMO_MATRIX_FILE} -format tf \
		-export newick -d3_base file -hclust_method average \
		-labels name,consensus \
		-o results/peakmo_clustering/peakmo_example
	@echo "	results/peakmo_clustering/peakmo_example"

## Cluster all motifs from RegulonDB
RDB_CLUSTER_DIR=results/regulondDB_clusters
RDB_CLUSTERS=${RDB_CLUSTER_DIR}/RDB_clusters
cluster_rdb:
	@echo "Clustering all matrices from RegulonDB"
	matrix-clustering -v ${V} -i data/RDB_PSSMs.tf -format transfac -o ${RDB_CLUSTERS}
	@echo ${RDB_CLUSTERS}
