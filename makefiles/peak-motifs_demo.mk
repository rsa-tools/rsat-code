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
	@echo "DIR_PEAKMO	${DIR_PEAKMO}"

################################################################
## Run peak-motifs on the peaks
MIN_OL=6
MAX_OL=8
DISCO=oligos,dyads,positions
PM_TASK=purge,seqlen,composition,disco,collect_motifs,motifs_vs_motifs,timelog,archive,synthesis,small_summary,motifs_vs_db,scan
MOTIF_PREFIX=Oct4_Chen2008_sites_from_Jaspar
peakmo_Chen_Oct4:
	@echo
	@echo "Running peak motifs	${MOTIF_PREFIX}"
	@mkdir -p ${DIR_PEAKMO}
	peak-motifs -v ${V} \
		-title ${MOTIF_PREFIX} \
		-i ${RSAT}/public_html/demo_files/peak-motifs_demo.fa \
		-markov auto \
		-disco ${DISCO} \
		-nmotifs 5 -minol ${MIN_OL} -maxol ${MAX_OL} \
		-no_merge_lengths -2str \
		-origin center \
		-motif_db jaspar_core_vertebrates tf ${RSAT}/public_html/motif_databases/JASPAR/jaspar_core_vertebrates_2013-11.tf \
		-scan_markov 1 -source galaxy \
		-task ${PM_TASK} \
		-prefix peak-motifs \
		-noov -img_format png \
		-outdir ${DIR_PEAKMO}
	@echo "	${DIR_PEAKMO}"

