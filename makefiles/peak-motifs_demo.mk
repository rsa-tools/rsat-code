################################################################
## Tests for peak-motifs

include ${RSAT}/makefiles/util.mk
MAKEFILE=${RSAT}/makefiles/peak-motifs_demo.mk

V=2

DIR_DATA=data
DIR_PEAKMO=results/peak-motifs_demo/Oct4_Chen2008_sites_from_Jaspar


################################################################
##  List default parameters
list_params:
	@echo "Peak-motifs demo parameters"
	@echo "	MOTIF_PREFIX	${MOTIF_PREFIX}"
	@echo "	DIR_DATA	${DIR_DATA}"
	@echo "	DIR_PEAKMO	${DIR_PEAKMO}"
	@echo "	DISCO		${DISCO}"
	@echo "	MIN_OL		${MIN_OL}"
	@echo "	MAX_OL		${MAX_OL}"
	@echo "	PM_TASK		${PM_TASK}"
	@echo "	PM_CMD		${PM_CMD}"

################################################################
## Run peak-motifs on the peaks
MIN_OL=6
MAX_OL=7
DISCO=oligos,dyads,positions
PM_TASK=purge,seqlen,composition,disco,merge_motifs,split_motifs,motifs_vs_motifs,timelog,archive,motifs_vs_db,scan,small_summary,synthesis
MOTIF_PREFIX=Oct4_Chen2008_sites_from_Jaspar
PM_CMD=peak-motifs -v ${V} \
		-title ${MOTIF_PREFIX} \
		-i ${RSAT}/public_html/demo_files/peak-motifs_demo.fa \
		-markov auto \
		-disco ${DISCO} \
		-nmotifs 5 -minol ${MIN_OL} -maxol ${MAX_OL} \
		-no_merge_lengths -2str \
		-origin center \
		-source galaxy \
		-motif_db jaspar_core_vertebrates tf ${RSAT}/public_html/motif_databases/JASPAR/jaspar_core_vertebrates_2013-11.tf \
		-scan_markov 1 -source galaxy \
		-task ${PM_TASK} \
		-prefix peak-motifs \
		-noov -img_format png \
		-outdir ${DIR_PEAKMO}
peakmo_demo:
	@echo
	@echo "Running peak motifs demo	${MOTIF_PREFIX}"
	@mkdir -p ${DIR_PEAKMO}
	${PM_CMD}
	@echo "	${DIR_PEAKMO}"

