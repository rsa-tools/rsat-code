################################################################
## Demo for the RSAT tool matrix-clustering
##
## Authors: Jaime Castro & Jacques van Helden
## Date: Jan-April 2014


include ${RSAT}/makefiles/util.mk
MAKEFILE=${RSAT}/makefiles/matrix-clustering_demo.mk

## Define a set of demo files
PEAKMO_DEMO_MATRICES=${RSAT}/public_html/demo_files/matrix-clustering_demo_peak-motifs_matrices.tf
## We should add some other demo files (e.g. RegulonDB)

## Select one of the demo files for the test
MATRIX_FILE=${PEAKMO_DEMO_MATRICES}

V=2

################################################################
## Compare all matrices from the input file, with specific parameters
## to ensure that all distances are computed.
COMPA_DIR=results/pairwise_motif_comparisons
COMPA_FILE=${COMPA_DIR}/peak-motifs_7nt_merged_oligos_positions_compa.tab
compa_alone:
	@echo
	@echo "Motif comparisons"
	compare-matrices  -v 1 \
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
################################################################
## Demo 1: clustering between motifs discovered by peak-motifs
cluster_peakmo:
	@echo "Clustering of motifs discovered by peak-motifs"
	matrix-clustering -v ${V} \
		-i ${PEAKMO_DEMO_MATRICES} -format tf \
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
