################################################################
## Demo for the RSAT tool matrix-clustering
##
## Authors: Jaime Castro & Jacques van Helden
## Date: 2014 - 2015


## include ${RSAT}/makefiles/util.mk
## We include compare-matrices_demo.mk to get the parameters (matrix files)
include ${RSAT}/makefiles/compare-matrices_demo.mk
MAKEFILE=${RSAT}/makefiles/matrix-clustering_demo.mk

################################################################
## Parameters for the analysis
MIN_NCOR=0.4
MIN_COR=0.6
HCLUST_METHOD=average
RETURN_FIELDS=heatmap,align_consensus
MIN_W=5
V=2

## Define a set of demo files
PEAKMO_PREFIX=peak-motifs_result_Chen_Oct4
FOOTPRINT_DISCO_PREFIX=footprint-discovery_LexA
OCT4_PREFIX=peak-motifs_Oct4
ES_CELLS_PREFIX=ES_CELL_ANALYSIS

## Choose a particular demo set
MATRIX_PREFIX=${PEAKMO_PREFIX}

## Define file locations based on the chosen demo set
MATRIX_DIR=${RSAT}/public_html/demo_files
MATRIX_FILE=${MATRIX_DIR}/${MATRIX_PREFIX}_matrices.tf
FILE_TABLE=${MATRIX_DIR}/${MATRIX_PREFIX}_motif_set.tab
RANGE_TABLE=${MATRIX_DIR}/${MATRIX_PREFIX}_metric_ranges.txt

list_param:
	@echo "MATRIX_PREFIX		${MATRIX_PREFIX}"
	@echo "MATRIX_DIR		${MATRIX_DIR}"
	@echo "MATRIX_FILE		${MATRIX_FILE}"
	@echo "PERMUTED_PREFIX		${PERMUTED_PREFIX}"
	@echo "PERMUTED_DIR		${PERMUTED_DIR}"
	@echo "PERMUTED_MATRIX_FILE	${PERMUTED_MATRIX_FILE}"
	@echo "JASPAR_VERSION		${JASPAR_VERSION}"
	@echo "JASPAR_GROUPS		${JASPAR_GROUPS}"
	@echo "JASPAR_GROUP		${JASPAR_GROUP}"
	@echo "CISBP_GROUPS		${CISBP_GROUPS}"
	@echo "CISBP_GROUP		${CISBP_GROUP}"
#	@echo "		${}"




################################################################
## Run matrix-clustering on one demo set (the particular cases will be
## specified below)
TITLE='matrix-clustering result'
CLUSTER_PREFIX=${MATRIX_PREFIX}_hclust-${HCLUST_METHOD}_Ncor${MIN_NCOR}_cor${MIN_COR}
CLUSTER_DIR=results/matrix-clustering_results/${MATRIX_PREFIX}/${HCLUST_METHOD}_linkage/Ncor${MIN_NCOR}_cor${MIN_COR}
CLUSTER_FILE_PREFIX=${CLUSTER_DIR}/${CLUSTER_PREFIX}
CLUSTER_CMD=matrix-clustering -v ${V} \
		-i ${MATRIX_FILE} -matrix_format tf -motif_collection_name '${MATRIX_PREFIX}' \
		-title '${TITLE}' \
		-lth Ncor ${MIN_NCOR} \
		-lth cor ${MIN_COR} \
		-lth w ${MIN_W} \
		-hclust_method ${HCLUST_METHOD} \
		-label_in_tree name ${OPT} \
		-o ${CLUSTER_FILE_PREFIX} \
		-metric_build_tree ${METRIC_BUILD_TREE} \
		-return ${RETURN_FIELDS} 

CLUSTER_MULTI_SET_CMD=matrix-clustering -v ${V} \
		-file_table ${FILE_TABLE} \
		-title '${TITLE}' \
		-lth Ncor ${MIN_NCOR} \
		-lth cor ${MIN_COR} \
		-lth w ${MIN_W} \
		-hclust_method ${HCLUST_METHOD} \
		-label_in_tree name ${OPT} \
		-o ${CLUSTER_FILE_PREFIX} \
		-metric_build_tree ${METRIC_BUILD_TREE} \
		-return ${RETURN_FIELDS} 

NB_CLUSTER_CMD=matrix-clustering -v ${V} \
	        -i ${MATRIX_FILE} -matrix_format tf -motif_collection_name '${MATRIX_PREFIX}' \
		-range_th_table ${RANGE_TABLE} \
		-title '${TITLE}' \
		-hclust_method ${HCLUST_METHOD} \
		-label_in_tree name ${OPT} \
		-o ${CLUSTER_FILE_PREFIX} \
		-metric_build_tree ${METRIC_BUILD_TREE} \
		-return ${RETURN_FIELDS} 

_cluster:
	@echo
	@echo "Running matrix-clustering	${MATRIX_PREFIX}	${OPT}"
	${MAKE} my_command MY_COMMAND="${CLUSTER_CMD}"
	@echo "		${CLUSTER_CMD}"
	@echo "		${CLUSTER_FILE_PREFIX}_SUMMARY.html"

_cluster_multi:
	@echo
	@echo "Running matrix-clustering	${MATRIX_PREFIX}	${OPT}"
	${MAKE} my_command MY_COMMAND="${CLUSTER_MULTI_SET_CMD}"
	@echo "		${CLUSTER_MULTI_SET_CMD}"
	@echo "		${CLUSTER_FILE_PREFIX}_SUMMARY.html"

_nb_cluster:
	@echo
	@echo "Running matrix-clustering	${MATRIX_PREFIX}	${OPT}"
	${MAKE} my_command MY_COMMAND="${NB_CLUSTER_CMD}"
	@echo "\n${NB_CLUSTER_CMD}"

## Cluster motifs resulting from two independent analysis of peak-motifs (Chen data set) with Oct4 and Sox2 peaks. 
cluster_peakmotifs_ES_cells_analysis:
	@echo
	@echo "Running matrix-clustering on motifs discovered by peak-motifs (Oct4, Sox2 and Nanog dataset from Chen 2008)"
	${MAKE} _cluster_multi MATRIX_PREFIX=${ES_CELLS_PREFIX} \
		TITLE='Oct4-Sox2-Nanog motifs from peak motifs' \
		METRIC_BUILD_TREE=Ncor


## Cluster motifs resulting from peak-motifs (Chen Oct4 data set)
cluster_peakmotifs_Oct4:
	@echo
	@echo "Running matrix-clustering on motifs discovered by peak-motifs (Oct 4 dataset from Chen 2008)"
	${MAKE} _cluster MATRIX_PREFIX=${OCT4_PREFIX} \
		TITLE='Oct4 motifs peak motifs' \
		COLLECTION=${OCT4_PREFIX} \
		METRIC_BUILD_TREE=Ncor


## Cluster motifs resulting from peak-motifs (Chen Oct4 data set)
cluster_peakmotifs_Oct4_nb_clusters:
	@echo
	@echo "Running matrix-clustering on motifs discovered by peak-motifs (Oct 4 dataset from Chen 2008)"
	${MAKE} _nb_cluster MATRIX_PREFIX=${OCT4_PREFIX} \
		TITLE='Oct4 motifs peak motifs' \
		COLLECTION=${OCT4_PREFIX} \
		METRIC_BUILD_TREE=Ncor \
                RETURN_FIELDS=nb_clusters


cluster_peakmotifs_Oct4_roots_only:
	@echo
	@echo "Running matrix-clustering on motifs discovered by peak-motifs (Oct 4 dataset from Chen 2008)"
	${MAKE} _cluster MATRIX_PREFIX=${OCT4_PREFIX} \
		TITLE='Oct4 motifs peak motifs' \
		COLLECTION=${OCT4_PREFIX} \
		METRIC_BUILD_TREE=Ncor \
		RETURN_FIELDS=root_matrices

## Cluster motifs resulting from peak-motifs (Chen Oct4 data set),
## without any threshold
cluster_peakmotifs_Oct4_no_threshold:
	@echo
	@echo "Running matrix-clustering on motifs discovered by peak-motifs (Oct 4 dataset from Chen 2008)"
	@${MAKE} _cluster MATRIX_PREFIX=${PEAKMO_PREFIX} MIN_NCOR=-1 MIN_COR=-1 \
		TITLE='Peak-motifs results for Oct4 ChIP-seq peaks - no thresholds' \
		COLLECTION=${PEAKMO_PREFIX}

## Permutation test with peak-motifs (Chen Oct4 data set)
cluster_peakmotifs_Oct4_permute:
	@${MAKE} ${MATRIX_PREFIX}=${PEAKMO_PREFIX} permute_matrices cluster_permuted_matrices \
		TITLE='Permuted matrices from peak-motifs results with Oct4 ChIP-seq peaks' \
		COLLECTION=${PEAKMO_PREFIX}

## Rndomize input matrices by permuting their columns
PERMUTED_PREFIX=${MATRIX_PREFIX}_permuted
PERMUTED_DIR=results/matrix-clustering_results/${PERMUTED_PREFIX}
PERMUTED_MATRIX_FILE=${PERMUTED_DIR}/${PERMUTED_PREFIX}_matrices.tf
permute_matrices: list_param
	@echo
	@echo "Permuting matrices	${MATRIX_FILE}"
	@mkdir -p ${PERMUTED_DIR}
	@permute-matrix -i ${MATRIX_FILE} \
		-in_format transfac -out_format transfac \
		-o ${PERMUTED_MATRIX_FILE}
	@echo "	${PERMUTED_MATRIX_FILE}"

## Run clustering on permuted matrices
cluster_permuted_matrices:
	@echo
	@echo "Clustering permuted matrices	${PERMUTED_MATRIX_FILE}"
	${MAKE} _cluster MATRIX_PREFIX=${PERMUTED_PREFIX} MATRIX_FILE=${PERMUTED_MATRIX_FILE}

## Cluster motifs resulting from footprint-discovery (LexA in Enterobacteriales)
cluster_footprints:
	@echo
	@echo "Running matrix-clustering on motifs discovered by footprint-discovery (query gene=LexA; taxon=Enterobacteriales)"
	${MAKE} _cluster MATRIX_PREFIX=${FOOTPRINT_DISCO_PREFIX} \
		COLLECTION=${FOOTPRINT_DISCO_PREFIX} \
		TITLE='Matrices found in LexA orthologous promoters in Enterobacteriales' \
		COLLECTION='LexA_motifs' \
		METRIC_BUILD_TREE=Ncor \
		RETURN_FIELDS=align_consensus,heatmap


## Cluster all motifs from RegulonDB
RDB_CLUSTER_DIR=results/matrix-clustering_results/regulondDB_clusters
RDB_CLUSTERS=${RDB_CLUSTER_DIR}/RDB_clusters
RDB_PREFIX=regulonDB_2015-08-07
RDB_MATRICES=${RSAT}/public_html/motif_databases/REGULONDB/${RDB_PREFIX}.tf
cluster_regulondb:
	@echo
	@echo "Clustering all matrices from RegulonDB"
	${MAKE} _cluster MATRIX_PREFIX=${RDB_PREFIX} MATRIX_FILE=${RDB_MATRICES} \
		TITLE='RegulonDB database' \
		COLLECTION=${RDB_PREFIX} \
		METRIC_BUILD_TREE=Ncor \
		RETURN_FIELDS=align_consensus,heatmap

## Permutation test with RegulonDB
cluster_regulondb_permute:
	${MAKE} MATRIX_PREFIX=${RDB_PREFIX} MATRIX_FILE=${RDB_MATRICES} permute_matrices
	${MAKE} MATRIX_PREFIX=${RDB_PREFIX} cluster_permuted_matrices


################################################################
## Cluster all motifs from HOCOMOCO
HOCOMOCO_CLUSTER_DIR=results/matrix-clustering_results/HOCOMOCO_clusters
HOCOMOCO_CLUSTERS=${HOCOMOCO_CLUSTER_DIR}/HOCOMOCO_clusters
HOCOMOCO_PREFIX=HOCOMOCO_2015-08-07
HOCOMOCO_MATRICES=${RSAT}/public_html/motif_databases/HOCOMOCO/${HOCOMOCO_PREFIX}.tf
cluster_hocomoco:
	@echo
	@echo "Clustering all matrices from HOCOMOCO"
	${MAKE} _cluster MATRIX_PREFIX=${HOCOMOCO_PREFIX} MATRIX_FILE=${HOCOMOCO_MATRICES} \
		TITLE='HOCOMOCO motifs' \
		COLLECTION=${RDB_PREFIX} \
		METRIC_BUILD_TREE=Ncor \
		RETURN_FIELDS=align_consensus,heatmap


################################################################
## Cluster one jaspar group
JASPAR_GROUPS=nematodes fungi urochordates plants vertebrates insects all 
JASPAR_GROUP=vertebrates
JASPAR_VERSION=2015_03
JASPAR_PREFIX=jaspar_core_${JASPAR_GROUP}_${JASPAR_VERSION}
JASPAR_DIR=${RSAT}/public_html/motif_databases/JASPAR
JASPAR_MATRICES=${JASPAR_DIR}/${JASPAR_PREFIX}.tf
cluster_jaspar_all_groups:
	@for g in ${JASPAR_GROUPS}; do \
		${MAKE} cluster_jaspar_one_group JASPAR_GROUP=$${g} ; \
	done

cluster_jaspar_one_group:
	@echo
	@echo "Clustering all matrices from JASPAR ${JASPAR_GROUP}"
	${MAKE} _cluster MATRIX_PREFIX=${JASPAR_PREFIX} \
		MATRIX_FILE=${JASPAR_MATRICES} \
		MIN_COR=0.6 MIN_NCOR=0.4 \
		TITLE='Jaspar core ${JASPAR_GROUP} database' \
		COLLECTION=${JASPAR_PREFIX} \
		METRIC_BUILD_TREE=Ncor \

## Permutation test with RegulonDB
cluster_jaspar_one_group_permute:
	${MAKE} MATRIX_PREFIX=${JASPAR_PREFIX} permute_matrices cluster_permuted_matrices

#############################
## Cluster one cisBP group
## Default: Mus_musculus
CISBP_GROUPS=Homo_sapiens Mus_musculus
CISBP_GROUP=Mus_musculus
CISBP_PREFIX=cisBP_${CISBP_GROUP}_2014-10
CISBP_DIR=${RSAT}/public_html/motif_databases/cisBP
CISBP_MATRICES=${CISBP_DIR}/${CISBP_PREFIX}.tf
cluster_cisbp_one_group:
	@echo
	@echo "Clustering all matrices from cisBP ${CISBP_GROUP}"
	${MAKE} _cluster MATRIX_PREFIX=${CISBP_PREFIX} \
		MATRIX_FILE=${CISBP_MATRICES} \
		MIN_COR=0.6 MIN_NCOR=0.4
		TITLE='cisBP ${CISBP_GROUP} database' \
		COLLECTION=${CISBP_PREFIX}

################################################################
## Send some jobs to the queue to check if it works
enqueue_some_jobs:
	${MAKE}  cluster_peakmotifs_Oct4 WHEN=queue HCLUST_METHOD=single
	${MAKE}  cluster_peakmotifs_Oct4 WHEN=queue HCLUST_METHOD=average
	${MAKE}  cluster_peakmotifs_Oct4 WHEN=queue HCLUST_METHOD=complete
