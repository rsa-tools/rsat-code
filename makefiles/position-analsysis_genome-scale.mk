################################################################
## Extract all the sequences surrounding a given genomic reference
## point (TSS, TTS, start codon, ...) for a given genome, and run
## position-analysis to discover k-mers characterized by specific
## positional profiles relative to this reference point.
################################################################

include ${RSAT}/makefiles/util.mk
MAKEFILE=makefiles/position-analsysis_genome-scale.mk

ORG=chlamydomonas_reinhardtii
FEATTYPE=gene
FROM=-500
TO=100
SEQTYPE=upstream
NOORF=-noorf
SEQ_PREFIX=${ORG}_${SEQTYPE}_${FEATTYPE}_${FROM}_${TO}${NOORF}
RESULT_DIR=${RSAT}/public_html/tmp/allseq_position-analysis/
SEQ_DIR=${RESULT_DIR}/${SEQ_PREFIX}_sequences
#RESULT_DIR=results/allseq_position-analysis/${SEQ_PREFIX}

WINDOW_SIZE=50

## Print the parameters
param:
	@echo "Parameters"
	@echo
	@echo "Sequences"
	@echo "	ORG		${ORG}"
	@echo "	FEATTYPE	${FEATTYPE}"
	@echo "	SEQTYPE		${SEQTYPE}"
	@echo "	FROM		${FROM}"
	@echo "	TO		${TO}"
	@echo "	NOORF		${NOORF}"
	@echo "	SEQ_PREFIX	${SEQ_PREFIX}"
	@echo "	RESULT_DIR	${RESULT_DIR}"
	@echo "	SEQ_DIR		${SEQ_DIR}"
	@echo "	ALLSEQ		${ALLSEQ}"
	@echo
	@echo "Position-analysis"
	@echo "	WINDOW_SIZE	${WINDOW_SIZE}"
	@echo "	KMER_SIZE	${KMER_SIZE}"
	@echo "	STRAND		${STRAND}"
	@echo "	KMER_OVERLAPS	${KMER_OVERLAPS}"
	@echo "	NB_CLUSTERS	${NB_CLUSTERS}"
	@echo "	POS_ORIGIN	${POS_ORIGIN}"
	@echo "	POS_OFFSET	${POS_OFFSET}"
	@echo "	NB_CLUSTERS	${NB_CLUSTERS}"
	@echo "	POS_PREFIX	${POS_PREFIX}"
	@echo "	POS_DIR		${POS_DIR}"
	@echo "	POS_FILE	${POS_FILE}"

## Collect all the sequences surrounding the user-selected reference point
ALLSEQ=${SEQ_DIR}/${SEQ_PREFIX}
allseq:
	@echo
	@echo "Collecting sequences	${SEQ_PREFIX}"
	@echo "	SEQ_DIR	${SEQ_DIR}"
	@mkdir -p ${SEQ_DIR}
	retrieve-seq -org ${ORG} -all -feattype ${FEATTYPE} -type ${SEQTYPE} \
		-from ${FROM} -to ${TO} -format fasta ${NOORF} \
		-o ${ALLSEQ}.fasta.gz
	@echo "	${ALL_SEQ}.fasta.gz"
	${MAKE} _allseq_len_distrib

## Compute the length distribution of the sequences
ALLSEQ_LEN=${ALLSEQ}_len_distrib_ci${WINDOW_SIZE}
_allseq_len_distrib:
	@echo
	@echo "Computing sequence length distribution"
	@sequence-lengths -i ${ALLSEQ}.fasta.gz | classfreq -v 1 -col 2 -ci ${WINDOW_SIZE} -o ${ALLSEQ_LEN}.tab
	@echo "	${ALLSEQ_LEN}.tab"
	@XYgraph -i ${ALLSEQ_LEN}.tab -xcol 3 -ycol 4,5,6 -legend \
		-title1 ${SEQ_PREFIX} -title2 'Length distribution' \
		-xleg 'Sequence length' -yleg '${FEATTYPE} nb' \
		-lines -xsize 600 -ysize 400 \
		-o ${ALLSEQ_LEN}.png
	@echo "	 ${ALLSEQ_LEN}.png"

## Detect positionally biased k-mers
KMER_SIZE=6
STRAND=-1str
KMER_OVERLAPS=-noov
POS_ORIGIN=end
NB_CLUSTERS=6
POS_OFFSET=-${TO}
POS_PREFIX=${SEQ_PREFIX}_${KMER_SIZE}nt_win${WINDOW_SIZE}${STRAND}${KMER_OVERLAPS}
POS_DIR=${RESULT_DIR}/${POS_PREFIX}
POS_FILE=${POS_DIR}/${POS_PREFIX}
POS_RETURN=distrib,occ,coverage,index,exp_occ,freq_per_window,freq_per_word,chi,sig,rank,graphs,matrices,clusters,html
POS_TASK=all
positions:
	@echo
	@echo "Running position-analysis"
	@mkdir -p ${POS_DIR}
	@echo "	POS_DIR	${POS_DIR}"
	position-analysis -v ${V} -i ${ALLSEQ}.fasta.gz -ci ${WINDOW_SIZE} -l ${KMER_SIZE} \
		${STRAND} ${KMER_OVERLAPS} \
		-origin ${POS_ORIGIN} -offset ${POS_OFFSET} \
		-sort -filter -return occ -task ${POS_TASK} \
		-lth occ 1 -lth sig 0 -clust_nb ${NB_CLUSTERS} -max_asmb_nb 5 \
		${OPT} -o ${POS_FILE}.tab
	@echo "	${POS_FILE}.tab"

# position-analysis -i $RSAT/public_html/tmp/_www/2015/10/14/tmp_sequence_2015-10-14.161508_Q2iLaN.fasta.purged -v 2 -sort -nofilter -lth_occ 1 -lth_sig 5 -clust_nb 2 -max_asmb_nb 2 -return chi,sig,rank,distrib,clusters,matrices,graphs,index -2str -noov -l 6 -ci 10 -origin center -o $RSAT/public_html/tmp/_www/2015/10/14/position-analysis_2015-10-14.161508_3u3Ab4/position-analysis.tab 