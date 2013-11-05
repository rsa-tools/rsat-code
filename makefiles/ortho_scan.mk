## ##############################################################
## Scan upstream sequences of orthologous genes with a
## position-specific scoring matrix (PSSM)

include ${RSAT}/makefiles/ortho_disco.mk
MAKEFILE=${RSAT}/makefiles/ortho_scan.mk

ortho_disco_targets:
	${MAKE} usage MAKEFILE=${RSAT}/makefiles/ortho_disco.mk


################################################################
## Scan the promoters of one group of orhtologous genes with one
## position-specific scoring matrix
FACTOR=CRP
GENE=fruR
MATRIX_DIR=data/RELEASE/Sites_mtx
MATRIX=${MATRIX_DIR}/${FACTOR}.RegulonDB.mtx 
RESULT_DIR=results/matches/${FACTOR}
LTH_SCORE=5
STR=-2str
PSEUDO=1
MKV=1
SCAN_DIR=results/matches/${FACTOR}/per_gene
SCAN_FILE=${SCAN_DIR}/matches_TF_${FACTOR}_up_${GENE}${PURGE_SUFFIX}_${TAXON}${STR}_ps${PSEUDO}_mkv${MKV}
SCAN_RETURN=limits,sites,rank
SCAN_CMD=matrix-scan -v ${V} -i ${SEQ} \
		-m ${MATRIX} \
		-lth score ${LTH_SCORE} -return ${SCAN_RETURN} \
		${STR} -origin -0 -pseudo ${PSEUDO} \
		-bginput -markov ${MKV} -o ${SCAN_FILE}.tab; echo ${SCAN_FILE}.tab
one_scan:
	@mkdir -p ${SCAN_DIR}
	@${MAKE} my_command MY_COMMAND="${SCAN_CMD}"

MAP_CMD=feature-map -i ${SCAN_FILE}.tab \
		-scalebar -scalestep 25 -xsize 800 -dot \
		-title "${FACTOR} matches in ${TAXON} ${GENE} promoters" \
		-o ${SCAN_FILE}.png ; echo ${SCAN_FILE}.png
one_map:
	@${MAKE} my_command MY_COMMAND="${MAP_CMD}"


one_scan_map:
	@${MAKE} my_command MY_COMMAND="${SCAN_CMD}; ${MAP_CMD}"

all: orthologs upstream scan map
