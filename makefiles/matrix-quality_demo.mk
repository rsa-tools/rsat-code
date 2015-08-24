################################################################
## template.mk makefile
## Usage: make -f matrix-quality_demo.mk

include ${RSAT}/makefiles/util.mk

MAKEFILE=${RSAT}/makefiles/matrix-quality_demo.mk

DATE=20150428
DEMO_FILE_DIR=${RSAT}/public_html/demo_files
MOTIFS_FILE=${DEMO_FILE_DIR}/Ballester_etal_elife_2014_4TFs_motifs.tf

BACKGROUND_MODEL=${DEMO_FILE_DIR}/all_human_ENCODE_DNAse_mk1_bg.ol

HNF6_SING_SET=${DEMO_FILE_DIR}/Ballester_etal_elife_2014_hg18_hnf6_singletons.fa
HNF6_SET=-seq hnf6_sing ${HNF6_SING_SET}
HNF6_PLOTS=-plot hnf6_sing nwd  -plot hnf6_sing occ_proba

HNF4A_SING_SET=${DEMO_FILE_DIR}/Ballester_etal_elife_2014_hg18_hnf4a_singletons.fa
HNF4A_SET=-seq hnf4a_sing ${HNF4A_SING_SET}
HNF4A_PLOTS=-plot hnf4a_sing nwd -plot hnf4a_sing occ_proba

FOXA1_SING_SET=${DEMO_FILE_DIR}/Ballester_etal_elife_2014_hg18_foxa1_singletons.fa
FOXA1_SET=-seq foxa1_sing ${FOXA1_SING_SET}
FOXA1_PLOTS=-plot foxa1_sing nwd -plot foxa1_sing occ_proba

CEBPA_SING_SET=${DEMO_FILE_DIR}/Ballester_etal_elife_2014_hg18_cebpa_singletons.fa
CEBPA_SET=-seq cebpa_sing ${CEBPA_SING_SET}
CEBPA_PLOTS=-plot cebpa_sing nwd  -plot cebpa_sing occ_proba

ZOO_CHIP_SEQ_SETS=${HNF6_SET} ${HNF4A_SET} ${FOXA1_SET} ${CEBPA_SET}
ZOO_CHIP_SEQ_PLOTS=${HNF6_PLOTS} ${HNF4A_PLOTS} ${FOXA1_PLOTS} ${CEBPA_PLOTS}

TASKS=

TITLE_ZOO='Compare enrichment in reported TF singleton binding sites'
MTX_FORMAT=transfac
RPLOT=-r_plots 
V=2
PSEUDO=1
BG_PSEUDO=0.01

MTXQ_ZOO_CMD=matrix-quality -v ${V} -html_title ${TITLE_ZOO}  -ms ${MOTIFS_FILE} -matrix_format ${MTX_FORMAT} \
	-pseudo ${PSEUDO}  -seq_format fasta ${ZOO_CHIP_SEQ_SETS}  -bgfile ${BACKGROUND_MODEL} \
	-bg_format oligo-analysis -bg_pseudo ${BG_PSEUDO}  ${ZOO_CHIP_SEQ_PLOTS} -o ${MTXQ_ZOO_OUT} ${TASKS} ${RPLOT}


MTXQ_OUT=./results/matrix_quality/${DATE}
MTXQ_ZOO_OUT=${MTXQ_OUT}/zoo_chip_enrichment

zoo_chip_quality:
	@echo ${MTXQ_ZOO_CMD}
	@${MTXQ_ZOO_CMD}


QUALITY_LIST=${MTXQ_ZOO_OUT}_quality_prefix_files.txt
COMPARE_QUALITY=${MTXQ_OUT}/zoo_chip_enrichment/compare_quality
compare_zoo_chip_quality:
	@ls -1 ${MTXQ_OUT}/*/*_synthesis.html | perl -pe 's/_synthesis.+//'> ${QUALITY_LIST}
	compare-qualities -quality_list ${QUALITY_LIST} -cluster -o ${COMPARE_QUALITY}
