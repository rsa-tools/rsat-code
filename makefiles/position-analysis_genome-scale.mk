################################################################
## Extract all the sequences surrounding a given genomic reference
## point (TSS, TTS, start codon, ...) for a given genome, and run
## position-analysis to discover k-mers characterized by specific
## positional profiles relative to this reference point.
################################################################

include ${RSAT}/makefiles/util.mk
MAKEFILE=makefiles/position-analysis_genome-scale.mk

ORG=chlamydomonas_reinhardtii
FEATTYPE=gene
WINDOW_SIZE=50
FROM=-500
TO=110
MIN_POS=${FROM}
MAX_POS=-1
REPEAT_MASK=-rm
SEQTYPE=upstream
NOORF=-noorf
SEQ_PREFIX=${ORG}_${SEQTYPE}_${FEATTYPE}_${FROM}_${TO}${NOORF}${REPEAT_MASK}
#RESULT_DIR=${RSAT}/public_html/tmp/allseq_position-analysis/
RESULT_DIR=results/allseq_position-analysis/${SEQ_PREFIX}
SEQ_DIR=${RESULT_DIR}/${SEQ_PREFIX}_sequences


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
	@echo "	REPEAT_MASK	${REPEAT_MASK}"
	@echo "	SEQ_PREFIX	${SEQ_PREFIX}"
	@echo "	RESULT_DIR	${RESULT_DIR}"
	@echo "	SEQ_DIR		${SEQ_DIR}"
	@echo "	ALLSEQ		${ALLSEQ}"
	@echo
	@echo "Position-analysis"
	@echo "	WINDOW_SIZE		${WINDOW_SIZE}"
	@echo "	KMER_SIZE		${KMER_SIZE}"
	@echo "	STRAND			${STRAND}"
	@echo "	KMER_OVERLAP		${KMER_OVERLAP}"
	@echo "	NB_CLUSTERS		${NB_CLUSTERS}"
	@echo "	POS_ORIGIN		${POS_ORIGIN}"
	@echo "	POS_OFFSET		${POS_OFFSET}"
	@echo "	NB_CLUSTERS		${NB_CLUSTERS}"
	@echo "	MAX_KMERS_TO_CLUSTER	${MAX_KMERS_TO_CLUSTER}"
	@echo "	MAX_ASMB_PER_CLUSTER	${MAX_ASMB_PER_CLUSTER}"
	@echo "	MAX_ASMB_NB		${MAX_ASMB_NB}"
	@echo "	POS_SIG			${POS_SIG}"
	@echo "	POS_MAX_GRAPHS		${POS_MAX_GRAPH}"
	@echo "	POS_PREFIX		${POS_PREFIX}"
	@echo "	POS_DIR			${POS_DIR}"
	@echo "	POS_FILE		${POS_FILE}"

## Collect all the sequences surrounding the user-selected reference point
ALLSEQ=${SEQ_DIR}/${SEQ_PREFIX}
allseq:
	@echo
	@echo "Collecting sequences	${SEQ_PREFIX}"
	@echo "	SEQ_DIR	${SEQ_DIR}"
	@mkdir -p ${SEQ_DIR}
	retrieve-seq -org ${ORG} -all -feattype ${FEATTYPE} -type ${SEQTYPE} \
		-from ${FROM} -to ${TO} -format fasta ${NOORF} ${REPEAT_MASK} \
		-o ${ALLSEQ}.fasta
	@echo "	${ALL_SEQ}.fasta"
	${MAKE} _allseq_len_distrib

## Compute the length distribution of the sequences
ALLSEQ_LEN=${ALLSEQ}_len_distrib_ci${WINDOW_SIZE}
_allseq_len_distrib:
	@echo
	@echo "Computing sequence length distribution"
	@sequence-lengths -i ${ALLSEQ}.fasta | classfreq -v 1 -col 2 -ci ${WINDOW_SIZE} -o ${ALLSEQ_LEN}.tab
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
KMER_OVERLAP=-noov
POS_ORIGIN=end
NB_CLUSTERS=12
POS_OFFSET=-${TO}
POS_PREFIX=${SEQ_PREFIX}_${KMER_SIZE}nt_win${WINDOW_SIZE}${STRAND}${KMER_OVERLAP}_sig${POS_SIG}
POS_DIR=${RESULT_DIR}/${POS_PREFIX}
POS_MAX_GRAPHS=50
POS_FILE=${POS_DIR}/${POS_PREFIX}
POS_RETURN=distrib,occ,coverage,index,exp_occ,freq_per_window,freq_per_word,chi,sig,rank,graphs,matrices,clusters,html
POS_TASK=all
POS_SIG=5
MAX_KMERS_TO_CLUSTER=150
MAX_ASMB_PER_CLUSTER=5
MAX_ASMB_NB=30
positions_with_bug:
	@echo
	@echo "Running position-analysis"
	@mkdir -p ${POS_DIR}
	@echo "	POS_DIR	${POS_DIR}"
	position-analysis -v ${V} -i ${ALLSEQ}.fasta -ci ${WINDOW_SIZE} -l ${KMER_SIZE} \
		${STRAND} ${KMER_OVERLAP} \
		-origin ${POS_ORIGIN} -offset ${POS_OFFSET} \
		-sort -filter -return occ -task ${POS_TASK} \
		-lth occ 1 -lth sig ${POS_SIG} -clust_nb ${NB_CLUSTERS} -max_asmb_nb 5 \
		${OPT} -o ${POS_FILE}.tab
	@echo "	${POS_FILE}.tab"

# position-analysis -i $RSAT/public_html/tmp/_www/2015/10/14/tmp_sequence_2015-10-14.161508_Q2iLaN.fasta.purged -v 2 -sort -nofilter -lth_occ 1 -lth_sig 5 -clust_nb 2 -max_asmb_nb 2 -return chi,sig,rank,distrib,clusters,matrices,graphs,index -2str -noov -l 6 -ci 10 -origin center -o $RSAT/public_html/tmp/_www/2015/10/14/position-analysis_2015-10-14.161508_3u3Ab4/position-analysis.tab 
POS_RETURN=chi,sig,rank,distrib,clusters,matrices,graphs,index
positions:
	@echo
	@echo "Running position-analysis"
	@mkdir -p ${POS_DIR}
	@echo "	POS_DIR	${POS_DIR}"
	position-analysis -v 2 \
		-i ${ALLSEQ}.fasta \
		-sort \
		-filter \
		-lth_occ 1 \
		-lth_sig ${POS_SIG} \
		-max_graphs ${POS_MAX_GRAPHS} \
		-clust_nb ${NB_CLUSTERS} \
		-toppat ${MAX_KMERS_TO_CLUSTER} \
		-max_asmb_per_cluster ${MAX_ASMB_PER_CLUSTER} \
		-max_asmb_nb ${MAX_ASMB_NB} \
		-return ${POS_RETURN} \
		-minpos ${MIN_POS} \
		-maxpos ${MAX_POS} \
		${STRAND} \
		${KMER_OVERLAP} \
		-l ${KMER_SIZE} \
		-ci ${WINDOW_SIZE} \
		-origin ${POS_ORIGIN} -offset ${POS_OFFSET} \
		${OPT} \
		-o ${POS_FILE}.tab
	@echo "	${POS_FILE}.tab"

## Compare results obtained with/without the -rm option
masked_or_not:
	@make -f makefiles/position-analysis_genome-scale.mk  allseq positions REPEAT_MASK=-rm; 
	@make -f makefiles/position-analysis_genome-scale.mk  allseq positions REPEAT_MASK=""; 
	@compare-scores -sc1 5 -sc2 4 -i results/allseq_position-analysis/${ORG}_upstream_gene_-500_110-noorf/${ORG}_upstream_gene_-500_110-noorf_6nt_win50-1str-noov_sig5/${ORG}_upstream_gene_-500_110-noorf_6nt_win50-1str-noov_sig5.tab -i results/allseq_position-analysis/${ORG}_upstream_gene_-500_110-noorf-rm/${ORG}_upstream_gene_-500_110-noorf-rm_6nt_win50-1str-noov_sig5/${ORG}_upstream_gene_-500_110-noorf-rm_6nt_win50-1str-noov_sig5.tab |awk '{print $0"\t"$3-$2}' | sort -n -k 4


################################################################
## Cluser the PSSM returned by position-analysis
MATRICES=${POS_FILE}_pssm_count_matrices
CLUSTERED_MATRICES=${MATRICES}_clustered
cluster_matrices:
	@echo
	@echo "Clustering the matrices"
	matrix-clustering -v ${V} \
		-i ${MATRICES}.tf \
		 -max_matrices 300 -matrix_format transfac \
		-hclust_method average \
		-title ${POS_PREFIX} \
		-metric_build_tree 'Ncor' \
		-lth w 5 \
		-lth cor 0.6 \
		-lth Ncor 0.4 \
		-label_in_tree name \
		-return json,heatmap \
		-o ${MATRICES}_clustered
	@echo " ${MATRICES}_clustered"
	@echo " ${MATRICES}_clustered_SUMMARY.html"

################################################################
## Compute a table with k-mer occurrences per position window
## (columns) in each position window (rows).
OCC_PER_SEQ=${POS_DIR}/${POS_PREFIX}_occ_per_seq
occ_per_seq:
	@echo
	@echo "Computing occurrences per position for each sequence"
	@mkdir -p '${OCC_PER_SEQ}'
	@position-analysis -v 2 \
		-i ${ALLSEQ}.fasta \
		-return occ_per_seq \
		${STRAND} \
		${KMER_OVERLAP} \
		-l ${KMER_SIZE} \
		-ci ${WINDOW_SIZE} \
		-origin ${POS_ORIGIN} -offset ${POS_OFFSET} \
		-selected_kmers ${POS_FILE}.tab \
		-minpos ${MIN_POS} \
		-maxpos ${MAX_POS} \
		${OPT} \
		-o ${OCC_PER_SEQ}.tab
	@echo "	${OCC_PER_SEQ}.tab"



