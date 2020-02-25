################################################################
## Tester for matrix-from-patterns


include ${RSAT}/makefiles/util.mk
MAKEFILE=${RSAT}/makefiles/matrix-from-patterns_demo.mk

## Verbosity
V=2

## Study case 1
PEAK_SET=Sox2
PEAK_SEQ=${RSAT}/public_html/demo_files/${PEAK_SET}_peaks.fasta
SEQ_NAME=${PEAK_SET}
SEQ_FILE=${PEAK_SEQ}

## Study case 2
METSET=yeast_MET_promoters
MET_SEQ=${RSAT}/public_html/demo_files/MET_up800-noorf.fasta
SEQ_NAME=${METSET}
SEQ_FILE=${MET_SEQ}

################################################################1
## List parameters
param:
	@echo "Parameters"
	@echo "	PEAK_SET			${PEAK_SET}"
	@echo "	PEAK_SEQ		${PEAK_SEQ}"
	@echo "	METSET			${METSET}"
	@echo "	MET_SEQ			${MET_SEQ}"
	@echo "	SEQ_NAME		${SEQ_NAME}"
	@echo "	SEQ_FILE		${SEQ_FILE}"
	@echo "	OLIGO_PREFIX		${OLIGO_PREFIX}"
	@echo "	OLIGOS			${OLIGOS}"
	@echo "	ASSEMBLY		${ASSEMBLY}"
	@echo "	CLUSTERING_OPT		${CLUSTERING_OPT}"
	@echo "	CLUSTERING_PREFIX	${CLUSTERING_PREFIX}"
	@echo "	DIR		${PSSM_DIR}"
	@echo "	PSSM_PREFIX		${PSSM_PREFIX}"
	@echo "	SIG_MATRICES		${SIG_MATRICES}"
	@echo "	COUNT_MATRICES		${COUNT_MATRICES}"
	@echo "	LOGOS			${LOGOS}"

################################################################
## Discover over-reprsented k-mers
OLIGO_DIR=results/oligos/${SEQNAME}
MKV=1
OL=6
OLIGO_PREFIX=${SEQ_NAME}_oligos-2str-noov_${OL}nt_mkv${MKV}
OLIGOS=${OLIGO_DIR}/${OLIGO_PREFIX}
oligos:
	@echo ""
	@echo "Discovering over-represented oligonucleotides	${SEQ_NAME}"
	@mkdir -p ${OLIGO_DIR}
	time oligo-analysis  -v ${V} -quick -i ${SEQ_FILE} -sort \
		-lth occ_sig 0 -uth rank 100 -return occ,proba,rank \
		-2str -noov -seqtype dna -l ${OL} -markov ${MKV} -pseudo 0.01 \
		-o ${OLIGOS}.tab
	@echo "	${OLIGOS}.tab"

################################################################
## Assemble patterns
ASSEMBLY=${OLIGOS}.asmb
assembly:
	@echo
	@echo "Assemblink k-mers"
	time pattern-assembly  -v ${V} -i ${OLIGOS}.tab \
		-2str -maxfl 1 -subst 1 -max_asmb_width 20 -toppat 100 -max_asmb_size 50 -max_asmb_width 20 -max_asmb_nb 10 \
		-o ${ASSEMBLY}
	@echo "	${ASSEMBLY}"

################################################################
## Run matrix-from-patterns
CLUSTERING_OPT=counts
CLUSTERING_PREFIX=clustering-${CLUSTERING_OPT}
PSSM_DIR=results/${SEQ_NAME}_${CLUSTERING_PREFIX}
PSSM_PREFIX=${PSSM_DIR}/${OLIGO_PREFIX}_pssm
SIG_MATRICES=${PSSM_PREFIX}_sig_matrices.tf
COUNT_MATRICES=${PSSM_PREFIX}_count_matrices.tf
LOGOS=`ls ${PSSM_PREFIX}_*.png | grep -v _rc | xargs`
QUICK=-quick
matrices:
	@echo
	@echo "Running matrix-from-patterns"
	@mkdir -p ${PSSM_DIR}
	@echo "	PSSM_DIR		${PSSM_DIR}"
	time matrix-from-patterns -v ${V} \
		-sites \
		-seq ${SEQ_FILE} \
		-asmb ${ASSEMBLY} \
		-bginput \
		-toppat 100 -max_asmb_nb 10 -max_asmb_width 20 -subst 1 -prefix oligos_${OL}nt \
		-flanks 2 -collect_method matrix-scan${QUICK} -logo \
		-cluster ${CLUSTERING_OPT} ${OPT} \
		-o ${PSSM_PREFIX}
	@echo "	CLUSTERING_OPT		${CLUSTERING_OPT}"
	@echo "	PSSM_PREFIX		${PSSM_PREFIX}"
	@echo "	SIG_MATRICES		${SIG_MATRICES}"
	@echo "	COUNT_MATRICES		${COUNT_MATRICES}"
	@echo "	LOGOS			${LOGOS}"

# ################################################################
# ## Run matrix-from-patterns with matrix-clustering option in order to avoid redundancy between the motifs
# matrices_clustered:
# 	@echo
# 	@echo "Running matrix-from-patterns"
# 	matrix-from-patterns -v ${V}  -sites -seq ${PEAK_SEQ} \
# 		-pl ${OLIGOS}.tab \
# 		-bgfile ${BG_FILE} \
# 		-toppat 100 -max_asmb_nb 10 -max_asmb_width 20 -subst 1 -prefix oligos_${OL}nt \
# 		-flanks 2 -collect_method matrix-scan-quick -logo \
# 		-clustering \
# 		-o ${PSSM_PREFIX}
# 	@echo "	${PSSM_PREFIX}"

