## ##############################################################
## Scan upstream sequences of orthologous genes with a
## position-specific scoring matrix (PSSM)

include ${RSAT}/makefiles/ortho_disco.mk
MAKEFILE=${RSAT}/makefiles/ortho_scan.mk

ortho_disco_targets:
	${MAKE} usage MAKEFILE=${RSAT}/makefiles/ortho_disco.mk

FACTOR=LexA
GENE=LEXA
MATRIX_DIR=data/RELEASE/Sites_mtx
MATRIX=${MATRIX_DIR}/${FACTOR}.RegulonDB.mtx 
RESULT_DIR=results/matches/${FACTOR}
LTH_SCORE=5
STR=-2str
PSEUDO=1
MKV=1
SCAN_DIR=results/matches/${FACTOR}/per_gene
SCAN_FILE=${SCAN_DIR}/matches_${FACTOR}_${TAXON}${STR}_ps${PSEUDO}_mkv${MKV}
scan:
	@mkdir -p ${SCAN_DIR}
	matrix-scan -v ${V} -i ${SEQ} \
		-m ${MATRIX} \
		-lth score ${LTH_SCORE} -return limits,sites,rank \
		${STR} -origin -0 -pseudo ${PSEUDO} \
		-bginput -markov ${MKV} -o ${SCAN_FILE}.tab
	@echo ${SCAN_FILE}.tab

map:
	feature-map -i ${SCAN_FILE}.tab \
		-scalebar -scalestep 25 -xsize 800 -dot \
		-title "${FACTOR} matches in ${TAXON} ${GENE} promoters" \
		-o ${SCAN_FILE}.png -htmap
	@echo ${SCAN_FILE}.png

all: orthologs upstream scan map