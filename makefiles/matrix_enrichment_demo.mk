################################################################
## template.mk makefile
## Usage: make -f matrix_enrichment_demo.mk


include ${RSAT}/makefiles/util.mk

MAKEFILE=${RSAT}/makefiles/matrix_enrichment_demo.mk

DATE=20160218
DEMO_FILE_DIR=${RSAT}/public_html/demo_files

V=2
## CapStarr-seq demo 

## 
MOTIFS_FILE=${DEMO_FILE_DIR}/Capstarr_selected_TFs.tf
MTX_FORMAT=transfac
PSEUDO=1

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


TITLE_CAPSTARR='TF enrichment in CapStarr-seq data'

## Outfile
MTXE_OUT=./results/matrix_enrichment/${DATE}
MTXE_CAP_OUT=${MTXE_OUT}/capstarr-seq

## test command
MATRIX_ENRICHMENT_CDM=matrix-enrichment -v ${V} \
	-matrix ${MOTIFS_FILE} -matrix_format ${MTX_FORMAT} -pseudo ${PSEUDO} \
	-bgfile ${BACKGROUND_MODEL} -bg_format ${BG_FORMAT} \
	${SEQUENCES} \
	-title ${TITLE_CAPSTARR} \
	-o ${MTXE_CAP_OUT}

capstar_seq_demo:
	@echo ${TITLE_CAPSTARR}
	@echo ${MATRIX_ENRICHMENT_CDM}
	@${MATRIX_ENRICHMENT_CDM}
