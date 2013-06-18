################################################################
## Tests for peak-motifs

include ${RSAT}/makefiles/util.mk
makefile=${RSAT}/makefiles/peak-motifs_tests.mk

V=2

DIR_DATA=data
DIR_RESULTS=results
DIR_PEAKS=${DIR_RESULTS}/peaks

DIR_READS=${DIR_DATA}/reads
FACTOR=
TEST=
CTRL=
TEST_READS=${DIR_READS}/${TEST}_reads.bed
CTRL_READS=${DIR_READS}/${CTRL}_reads.bed
GENOME=mm9

################################################################
##  List default parameters
list_params:
	@echo "DIR_DATA	${DIR_DATA}"
	@echo "DIR_RESULTS	${DIR_RESULTS}"
	@echo "DIR_READS	${DIR_READS}"
	@echo "DIR_PEAKS	${DIR_PEAKS}"
	@echo "GENOME		${GENOME}"
	@echo "FACTOR		${FACTOR}"
	@echo "TEST		${TEST}"
	@echo "CTRL		${CTRL}"
	@echo "PEAK_PREFIX	${PEAK_PREFIX}"

################################################################
## Peak calling with SWEMBL
SWEMBL_PEAKS_DIR=${DIR_PEAKS}/SWEMBL
SWEMBL_PREFIX=SWEMBL_${TEST}_vs_${CTRL}_peaks_R${SWEMBL_R}_nof
SWEMBL_PEAKS=${SWEMBL_PEAKS_DIR}/${SWEMBL_PREFIX}
SWEMBL_R=0.01
SWEBML_F=0
swembl:
	@echo
	@echo "Peak calling with SWEMBL"
	@echo "Log file	${SWEMBL_PEAKS}_log.txt"
	@mkdir -p ${SWEMBL_PEAKS_DIR}
	time (SWEMBL -i ${TEST_READS} -R ${SWEMBL_R} -B  -r ${CTRL_READS} -o ${SWEMBL_PEAKS}.bed) >& ${SWEMBL_PEAKS}_log.txt
	@echo "	${SWEMBL_PEAKS}.bed"
	@${MAKE} peak_len_distrib PEAKS=${SWEMBL_PEAKS}

################################################################
## Run SWEMBL with a series of values for the parameter R in order to
## test the impact fo R on peak numbers, and then the impact of peak
## number on motifs (enrichment, discovery).
SWEMBL_R_VALUES=0.2 0.17 0.15 0.12 0.1 0.05 0.02 0.01 0.005 0.002 0.001
SWEMBL_TASK=swembl
swembl_series:
	for r in ${SWEMBL_R_VALUES} ; do \
		${MAKE} SWEMBL_R=$${r} ${SWEMBL_TASK}; \
	done

################################################################
## Plot peak length distribution for one peak file
PEAK_PREFIX=${SWEMBL_PREFIX}
PEAKS=${SWEMBL_PEAKS}
peak_len_distrib:
	@echo "Computing peak length distribution ${PEAKS}"
	@echo " `grep -v '^#' ${PEAKS}.bed  | wc -l` peaks"
	grep -v '^#' ${PEAKS}.bed  | grep -v '^Region' | classfreq -v 1 -ci 20 -col 5 -o ${PEAKS}_len_distrib.tab
	@echo "	${PEAKS}_len_distrib.tab"
	XYgraph -i ${PEAKS}_len_distrib.tab \
		-lines -legend \
		-xcol 3 -xmin 0 -xmax 2000 -xgstep1 200 -xgstep2 50 -xsize 800 -xleg1 'Peak length' \
		-ycol 4,6 -ymin 0 -ysize 300 -yleg1 'Number of peaks' \
		-title1 '${PEAK_PREFIX}' \
		-title2 'Peak length distribution' \
		-o ${PEAKS}_len_distrib.png
	@echo "	${PEAKS}_len_distrib.png"

################################################################
## Fetch peak sequences 
EXTEND=0
FETCH_CHUNK=10000
PEAK_SEQ=${PEAKS}.fasta
fetch_one_peak_set:
	@echo
	@echo "Fetching peaks	${PEAKS}.bed"
	fetch-sequences -i ${PEAKS}.bed -header_format galaxy -extend ${EXTEND} -chunk ${FETCH_CHUNK} -genome ${GENOME} ${OPT} -o ${PEAK_SEQ}

## Fetch sequences for a series of SWEMBL peaks
fetch_swembl_peaks:
	@${MAKE} swembl_series SWEMBL_TASK=fetch_one_peak_set

################################################################
## Plot the number of peaks as a function of the value of parameter R
## for SWEMBL.
PEAK_NUMBERS=${SWEMBL_PEAKS_DIR}/peak_numbers_SWEMBL_${TEST}_vs_${CTRL}
PLOT_FORMAT=pdf 
swembl_peak_numbers:
	@echo
	@echo "Computing number of peaks returned by SWEMBL as a function of parametere R"
	@echo "#R.param	nb.peaks" > ${PEAK_NUMBERS}.tab
	${MAKE} swembl_series SWEMBL_TASK=_one_peak_number
	@echo "	${PEAK_NUMBERS}.tab"
	XYgraph -i ${PEAK_NUMBERS}.tab \
		-xcol 1 -xleg1 'SWEMBL parameter R' -xsize 400\
		-ycol 2 -yleg1 'Number of peaks' -ysize 250  \
		-lines \
		-title1 'Peak calling with SWEMBL' \
		-title2 'Dataset: HNF4 from Schmidt et al. (2010)' \
		-format ${PLOT_FORMAT} -o ${PEAK_NUMBERS}.${PLOT_FORMAT}
	@echo "	${PEAK_NUMBERS}.${PLOT_FORMAT}"

PEAK_NB=`grep -v '^\#' ${SWEMBL_PEAKS}.bed | wc -l  | perl -pe 's| +||'`
_one_peak_number:
	@echo "${SWEMBL_R}	${PEAK_NB}	${SWEMBL_PEAKS}.bed" >> ${PEAK_NUMBERS}.tab


################################################################
## Peak calling with MACS
##
## BEWARE: a single run with MACS costs 100 minutes on my Mac.
##
## Pour utiliser peaksplitter, il faut retourner le fichier "wig", ce
## qui permet à ce programme de trouver la forme des peaks et les
## couper.
##
## (Morgane) Les options à utiliser sont  --wig --single-profile.
## http://liulab.dfci.harvard.edu/MACS/README.html
##
## Voici les commandes de MACS, il faut bien sûr les lancer après
## avoir fait le mapping (et s'assurer au préalable de la qualité des
## reads, et les nettoyer si besoin, et connaître le type FASTQ
## utilisé pour adapter les bon paramètres pour le mapping):
##
##   MACS_MFOLD=5,30
##   MACS_MFOLD=10,30 if the dataset is of good quality (less permissive threshold)
##   MACS_PVAL=1e-5
##   ## bw: should be the average size of fragment before sonication default is 300
##   BW=300
## This value should be set to the sonication fragment size expected from wet experiment
##
## by default, the reads are normalized to the bigger dataset. if the
## control is larger than the experiment, add --to-small so that the
## treatment reads are not artificially changed
##
## ${MACS} -t ${READS_TREAT}.ENS.bed -c ${READS_CTL}.ENS.bed --format BED  --gsize ${GSIZE} --name "macs"  --mfold ${MACS_MFOLD} --pvalue ${MACS_PVAL} --bw ${BW} --wig --single-profile --verbose 2 --diag &> ${READS_TREAT}vs${READS_CTL}/MACS.out
##
## peak-splitter has to be installed separately
## it can be run within MACS with the option --call-subpeaks
##
## I still use it separately, so I can use filtered peaks (with additional control treatment and/or filtered by FDR), instead of the default MACS output
## BED is this peak file
## @${PEAK_SPLITTER} -f -p ${BED} -w ${WIG} -o ${OUTDIR} > ${OUTDIR}/peak_splitter.log
MACS=macs14
GENOME_SIZE=3000000000
TAG_SIZE=35
MACS_MFOLD=10,30
MACS_PVAL=1e-5
BW=300
MACS_UP_DIR=../../..
MACS_PEAKS_DIR=${DIR_PEAKS}/MACS
MACS_PREFIX=${TEST}_vs_${CTRL}_macs14_pval${MACS_PVAL}
macs:
	@echo
	@echo "Peak calling with MACS"
	@echo 
	@mkdir -p ${MACS_PEAKS_DIR}
	@echo "MACS directory	${MACS_PEAKS_DIR}"
	@echo "Log file	${MACS_PREFIX}_log.txt"
	time (cd ${MACS_PEAKS_DIR}; ${MACS} --treatment ${MACS_UP_DIR}/${TEST_READS} --control ${MACS_UP_DIR}/${CTRL_READS} --format BED  \
		--gsize ${GENOME_SIZE} --mfold ${MACS_MFOLD} --pvalue ${MACS_PVAL} \
		--bw ${BW} --wig --single-profile --verbose 2 --diag \
		--tsize ${TAG_SIZE} \
		--name ${MACS_PREFIX} > ${MACS_PREFIX}.out) &> ${MACS_PEAKS_DIR}/${MACS_PREFIX}_log.txt
	@echo "	${MACS_PREFIX}.out"
	R --vanilla --slave --file=${MACS_PREFIX}_model.r
	@${MAKE} sort_macs_peaks
	@${MAKE} fetch_macs_seqs


MACS_PEAKS=${MACS_PEAKS_DIR}/${MACS_PREFIX}_peaks
MACS_PEAKS_SORTED=${MACS_PEAKS}_sorted
MACS_SUMMITS=${MACS_PEAKS_DIR}/${MACS_PREFIX}_summits
MACS_SUMMITS_SORTED=${MACS_SUMMITS}_sorted
sort_macs_peaks:
	@echo
	@echo "Sorting MACS peaks and summits by decreasing score values"
	@sort -nr -k 5 ${MACS_PEAKS}.bed > ${MACS_PEAKS_SORTED}.bed
	@echo "	${MACS_PEAKS_SORTED}.bed"
	@sort -nr -k 5 ${MACS_SUMMITS}.bed > ${MACS_SUMMITS_SORTED}.bed
	@echo "	${MACS_SUMMITS_SORTED}.bed"

## Fetch peak sequences as deined by MACS + fixed-width peaks centered around MACS peak summits
fetch_macs_seqs:
	@${MAKE} PEAKS=${MACS_PEAKS_SORTED} EXTEND=0 fetch_one_peak_set
	@${MAKE} PEAKS=${MACS_SUMMITS_SORTED} EXTEND=100 fetch_one_peak_set

## Test various parameters for MACS
MACS_PVALUES=1e-5 1e-6 1e-7
#MACS_PVALUES=1e-8 1e-7 1e-6 1e-5
MACS_TASK=macs
macs_pval_series:
	@for pval in ${MACS_PVALUES}; do \
		${MAKE} ${MACS_TASK} MACS_PVAL=$${pval}; \
	done
#	${MAKE} macs MACS_MFOLD=5,30

################################################################
################################################################
####                                                        ####
####                    MOTIF ANALYSIS                      ####
####                                                        ####
################################################################
################################################################

################################################################
## Run peak-motifs on the peaks
MINOL=6
MAXOL=7
MIN_MKV=auto
MAX_MKV=auto
STR=-2str
NOOV=-noov
DISCO=oligos,dyads,positions
PM_TASK=purge,seqlen,composition,disco,merge_words,collect_motifs,motifs_vs_motifs,timelog,archive,synthesis,motifs_vs_db,scan
TOP_PEAKS=0
MOTIF_PREFIX=${PEAK_PREFIX}${SUMMITS}_top${TOP_PEAKS}
DIR_MOTIFS=${DIR_RESULTS}/motifs_${FACTOR}/${MOTIF_PREFIX}
GALAXY=-source galaxy
peakmo:
	@echo
	@echo "Running peak motifs	${MOTIF_PREFIX}"
	@mkdir -p ${DIR_MOTIFS}
	peak-motifs  -v ${V} -title "${MOTIF_PREFIX}" \
		-i ${PEAK_SEQ} \
		-minol ${MINOL} -maxol ${MAXOL} \
		-min_markov ${MIN_MKV} -max_markov ${MAX_MKV} \
		${STR} ${NOOV} \
		-disco ${DISCO} \
		-motif_db jaspar_core_vertebrates tf ${RSAT}/public_html/data/motif_databases/JASPAR/jaspar_core_vertebrates_2009_10.tf \
		-motif_db jaspar_pbm_mouse tf ${RSAT}/public_html/data/motif_databases/JASPAR/jaspar_pbm_mouse_2009_10.tf \
		${GALAXY} \
		-task ${PM_TASK} -prefix ${MOTIF_PREFIX} \
		-img_format png  \
		-top_peaks ${TOP_PEAKS} ${OPT} \
		-outdir ${DIR_MOTIFS}
	@echo "	${DIR_MOTIFS}"

## Run peak-motifs with MACS peaks
peakmo_macs:
	@${MAKE} peakmo  PEAK_PREFIX=${MACS_PREFIX} PEAKS=${MACS_PEAKS_SORTED}

peakmo_macs_summits:
	@${MAKE} peakmo  PEAK_PREFIX=${MACS_PREFIX} SUMMITS=_summits PEAKS=${MACS_SUMMITS_SORTED}

## Run peak-motifs with SWEMBL peaks
peakmo_swembl:
	@${MAKE} peakmo  PEAK_PREFIX=${SWEMBL_PREFIX} PEAKS=${SWEMBL_PEAKS}


## Run peak-motifs with increasing number of top peaks
TOP_VALUES=000100 000200 000300 000500 001000 002000 003000 005000 007000 010000 020000 030000 050000 070000 100000
TOP_TASK=peakmo_macs
peakmo_top_series:
	@for top in ${TOP_VALUES}; do \
		${MAKE} ${TOP_TASK} TOP_PEAKS=$${top}; \
	done


################################################################
## Run matrix quality to compare peak enrichment between SWEMBL and/or
## MACS peak sets.
REF_MATRICES=${DIR_DATA}/ref_matrices/${FACTOR}_matrices.tf
MATRICES=${REF_MATRICES}
ORG=Mus_Musculus_EnsEMBL
QUALITY_BG=${RSAT}/data/genomes/${ORG}/oligo-frequencies/2nt_upstream-noorf_${ORG}-noov-1str.freq.gz 
PEAK_CALLERS=MACS_SWEMBL
MATRIX_SOURCE=ref_matrices
QUALITY_PREFIX=${FACTOR}_${PEAK_CALLERS}_${MATRIX_SOURCE}
QUALITY_DIR=${DIR_RESULTS}/peak_quality/${QUALITY_PREFIX}
QUALITY_OUT=${QUALITY_DIR}/${QUALITY_PREFIX}
QUALITY_TASK=all
MACS_PEAK_SETS=-seq MACS_summits_P1e-7 ${DIR_PEAKS}/MACS/${TEST}_vs_${CTRL}_macs14_pval1e-7_peaks.fasta -perm MACS_summits_P1e-7 0 \
		-seq MACS_peaks_P1e-7 ${DIR_PEAKS}/MACS/${TEST}_vs_${CTRL}_macs14_pval1e-7_peaks.fasta -perm MACS_peaks_P1e-7 0 \
		-seq MACS_summits_P1e-6 ${DIR_PEAKS}/MACS/${TEST}_vs_${CTRL}_macs14_pval1e-6_peaks.fasta -perm MACS_summits_P1e-6 0 \
		-seq MACS_peaks_P1e-6 ${DIR_PEAKS}/MACS/${TEST}_vs_${CTRL}_macs14_pval1e-6_peaks.fasta -perm MACS_peaks_P1e-6 0 \
		-seq MACS_summits_P1e-5 ${DIR_PEAKS}/MACS/${TEST}_vs_${CTRL}_macs14_pval1e-5_peaks.fasta -perm MACS_summits_P1e-5 10 \
		-seq MACS_peaks_P1e-5 ${DIR_PEAKS}/MACS/${TEST}_vs_${CTRL}_macs14_pval1e-5_peaks.fasta -perm MACS_peaks_P1e-5 10

SWEMBL_PEAK_SETS=-seq SWEMBL_R0.1 ${SWEMBL_PEAKS_DIR}/SWEMBL_${TEST}_vs_${CTRL}_peaks_R0.1_nof.fasta -perm SWEMBL_R0.1 0 \
		-seq SWEMBL_R0.02 ${SWEMBL_PEAKS_DIR}/SWEMBL_${TEST}_vs_${CTRL}_peaks_R0.02_nof.fasta  -perm SWEMBL_R0.02 0 \
		-seq SWEMBL_R0.05 ${SWEMBL_PEAKS_DIR}/SWEMBL_${TEST}_vs_${CTRL}_peaks_R0.05_nof.fasta  -perm SWEMBL_R0.05 0 \
		-seq SWEMBL_R0.01 ${SWEMBL_PEAKS_DIR}/SWEMBL_${TEST}_vs_${CTRL}_peaks_R0.01_nof.fasta  -perm SWEMBL_R0.01 0 \
		-seq SWEMBL_R0.005 ${SWEMBL_PEAKS_DIR}/SWEMBL_${TEST}_vs_${CTRL}_peaks_R0.005_nof.fasta  -perm SWEMBL_R0.005 0 \
		-seq SWEMBL_R0.002 ${SWEMBL_PEAKS_DIR}/SWEMBL_${TEST}_vs_${CTRL}_peaks_R0.002_nof.fasta  -perm SWEMBL_R0.002 0 \
		-seq SWEMBL_R0.001 ${SWEMBL_PEAKS_DIR}/SWEMBL_${TEST}_vs_${CTRL}_peaks_R0.001_nof.fasta  -perm SWEMBL_R0.001 10
peak_quality:
	@mkdir -p ${QUALITY_DIR}
	matrix-quality -v ${V} \
		-ms ${MATRICES} -top 7 -matrix_format transfac \
		-seq_format fasta \
		${MACS_PEAK_SETS} \
		${SWEMBL_PEAK_SETS} \
		-pseudo 1 -img_format png,pdf \
		-bgfile ${QUALITY_BG} \
		-bg_format oligos -decimals 1 \
		-task ${QUALITY_TASK} ${OPT} \
		-o ${QUALITY_OUT}
	@echo ${QUALITY_OUT}

peak_quality_swembl_only:
	@${MAKE} peak_quality PEAK_CALLERS=SWEMBL MACS_PEAK_SETS=""

peak_quality_macs_only:
	@${MAKE} peak_quality PEAK_CALLERS=MACS SWEMBL_PEAK_SETS=""

################################################################
## Run matrix-quality with discovered matrices. For this , we use a
## 20-fold cross-validation because the matrices resulting from
## peak-motifs can contain several thousands of sites.
DISCO_MATRICES=${DIR_RESULTS}/motifs/${FACTOR}_discovered_matrices.tf
peak_quality_disco:
	${MAKE} MATRICES=${DISCO_MATRICES} MATRIX_SOURCE=discovered_matrices peak_quality OPT='-kfold 20'

