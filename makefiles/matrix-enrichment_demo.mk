################################################################
## Usage: make -f matrix_enrichment_demo.mk

include ${RSAT}/makefiles/util.mk
MAKEFILE=${RSAT}/makefiles/matrix-enrichment_demo.mk

DATE=20160427
DEMO_FILE_DIR=${RSAT}/public_html/demo_files

V=2

## Outfile
MTXE_OUT=./results/matrix_enrichment/${DATE}
MTXE_CAP_OUT=${MTXE_OUT}/capstarr-seq


#######################
## CapStarr-seq demo ##
#######################

## 
MOTIFS_FILE=${DEMO_FILE_DIR}/Capstarr_selected_TFs.tf
MTX_FORMAT=transfac
PSEUDO=1
OPT=

## Bg model based in all humam promoter regions
BACKGROUND_MODEL=${DEMO_FILE_DIR}/All_hpromoterRegions_BG_model_mkv_2.oligos
BG_FORMAT=oligos

## Test sequence sets
HELA_SEQ=${DEMO_FILE_DIR}/CapStarrseq_Active_Prom_HELA_merge_IP.fasta
K562_SEQ=${DEMO_FILE_DIR}/CapStarrseq_Active_Prom_K562_merge_IP.fasta

## Control sequence set
INACTIVE_SEQ=${DEMO_FILE_DIR}/CapStarrseq_InactiveProm_FDR95_All_samples.fasta

## sequence options
HELA_SEQ_OPTIONS=-seq HELA ${HELA_SEQ}
K562_SEQ_OPTIONS=-seq K562 ${K562_SEQ}
INACTIVE_SEQ_OPTIONS=-seq inactive ${INACTIVE_SEQ}

SEQUENCES=${HELA_SEQ_OPTIONS} ${K562_SEQ_OPTIONS} ${INACTIVE_SEQ_OPTIONS}

TITLE_CAPSTARR="'TF enrichment in CapStarr-seq data'"


###################
## ES cells demo ##
###################

TITLE_ES_CELLS="'OCT4 motifs in OCT_SOX_NANOG sequences'"
MTXE_ES_CELLS=${MTXE_OUT}/OCT4_enrichment/OCT4_enrichment
OCT_MOTIFS=${DEMO_FILE_DIR}/MEME_ChIP_Oct4_matrices.tf

OCT_SEQ=${DEMO_FILE_DIR}/Oct4_peaks.fasta
NANOG_SEQ=${DEMO_FILE_DIR}/Nanog_peaks.fasta
SOX_SEQ=${DEMO_FILE_DIR}/Sox2_peaks.fasta

## sequence options
OCT_SEQ_OPTIONS=-seq Oct4 ${OCT_SEQ}
NANOG_SEQ_OPTIONS=-seq Nanog ${NANOG_SEQ}
SOX_SEQ_OPTIONS=-seq Sox2 ${SOX_SEQ}

ES_CELLS_SEQ=${OCT_SEQ_OPTIONS} ${NANOG_SEQ_OPTIONS} ${SOX_SEQ_OPTIONS}


## Command
MATRIX_ENRICHMENT_CDM=matrix-enrichment -v ${V} \
	-matrix ${MOTIFS_FILE} -matrix_format ${MTX_FORMAT} -pseudo ${PSEUDO} \
	-bgfile ${BACKGROUND_MODEL} -bg_format ${BG_FORMAT} \
	${SEQUENCES} \
	-title ${TITLE} ${OPT}\
	-o ${PREFIX}

## Enrichment with BG model based in one sequences
MATRIX_ENRICHMENT_BG_INPUT_CDM=matrix-enrichment -v ${V} \
	-matrix ${MOTIFS_FILE} -matrix_format ${MTX_FORMAT} -pseudo ${PSEUDO} \
	-bg_input all -markov 2 \
	${SEQUENCES} \
	-title ${TITLE} ${OPT} \
	-o ${PREFIX}


_enrichment:
	@echo
	@echo "Running 	${PREFIX}	${OPT}"
	${MAKE} my_command MY_COMMAND="${MATRIX_ENRICHMENT_CDM}"
	@echo "		${MATRIX_ENRICHMENT_CDM}"
	@echo "		${PREFIX}_report.html"

_bg_input:
	@echo
	@echo "Running 	${PREFIX}	${OPT}"
	${MAKE} my_command MY_COMMAND="${MATRIX_ENRICHMENT_BG_INPUT_CDM}"
	@echo "	${MATRIX_ENRICHMENT_BG_INPUT_CDM}"
	@echo "	${PREFIX}_report.html"

capstar_seq_demo:
	@echo
	@echo ${TITLE_CAPSTARR}
	@echo ${MATRIX_ENRICHMENT_CDM}
	${MAKE} _enrichment TITLE=${TITLE_CAPSTARR} PREFIX=${MTXE_CAP_OUT}

es_cells_demo:
	@echo
	@echo ${TITLE_ES_CELLS}
	${MAKE} _bg_input TITLE=${TITLE_ES_CELLS} PREFIX=${MTXE_ES_CELLS} MOTIFS_FILE=${OCT_MOTIFS} SEQUENCES="${ES_CELLS_SEQ}"
