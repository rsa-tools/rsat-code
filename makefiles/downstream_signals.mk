################################################################
## Motif discovery in all the gene downstream regions, in order to
## detect poly-adenylation and termination signals. This makefile
## reproduces and extends the approaches developed in van Helden et
## al. (2001), in particular the detection of positionally biased
## oligonucleotides (k-mers).
##
## van Helden,J., del Olmo,M. and Perez-Ortín,J.E. (2000) Statistical
## analysis of yeast genomic downstream sequences reveals putative
## polyadenylation signals. Nucleic Acids Res, 28, 1000-1010.
##
## Note: it seems that in EnsemblGenomes a large number of 3' UTR are
## exactly 100bp long, which suggetss that when the annotators ignore
## the actual length, the set it to 100bp by default. To avoid
## including too many coding sequences (and the resulting effect of
## aligned stop codons on oligo profiles) we set the upstream limit to
## 100.

include ${RSAT}/makefiles/util.mk


MAKEFILE=${RSAT}/makefiles/downstream_signals.mk

ORG=Saccharomyces_cerevisiae_R64-1-1
RES_DIR=results/downstream_signals
ORG_DIR=${RES_DIR}/${ORG}

dir:
	@mkdir -p ${ORG_DIR}

################################################################
## Collect all sequences around transcription termination site (TTS).
SEQ_PREFIX=${ORG}_downstream_${FEATTYPE}_from-${UP}_to${DOWN}
SEQ=${ORG_DIR}/${SEQ_PREFIX}.fasta.gz
FEATTYPE=mRNA
UP=100
DOWN=100
seq: dir
	@echo
	@echo "Retrieving sequences around ${FEATTYPE} ends for organism	${ORG}"
	retrieve-seq  -org ${ORG} -feattype ${FEATTYPE} -type downstream -format fasta \
		-label organism_name,id,name,sequence_type,current_from,current_to,ctg,orf_strand,reg_left,reg_right \
		-from -${UP} -to +${DOWN} -all \
		-o ${SEQ}
	@echo "	${SEQ}"


################################################################
## Run position-analysis to detect positionally biased
## oligonucleotides around the TTS.
OL=6
CI=20
NOOV=-noov
STR=-1str
POS_PREFIX=${SEQ_PREFIX}_${OL}nt${STR}${NOOV}_ci${CI}
POS_DIR=${ORG_DIR}/${POS_PREFIX}
POS_FILE=${POS_DIR}/${POS_PREFIX}
OPT_POS_TO_MATRICES=-return clusters,matrices -clust_nb 8 -max_asmb_nb 5
pos:
	@echo
	@echo "Position-analysis with ${OL}nt ${OPT}	${ORG}"
	@mkdir -p ${POS_DIR}
	position-analysis -i ${SEQ} -v ${V} -sort -nofilter -lth_occ 1 -max_graphs 100 \
		-return chi,sig,rank,distrib,exp_occ,graphs,index ${STR} ${NOOV} -l ${OL} -ci ${CI} \
		-origin start -offset ${UP} \
		${OPT} \
		-o ${POS_FILE}.tab
	@echo "	${POS_FILE}.tab"
	@echo "	${POS_FILE}_index.html"


################################################################
## Run position-analysis with a position-based Markov model: in each
## positional window, estimate the background k-mer frequencies (prior
## probabilities) from the composition of smaller k-mers in the same
## window.
##
## We take a Markov order of 2 to discard the impact of codons, which
## create strong signals at the transition between coding and
## non-coding sequences.
MKV=2
pos_mkv:
	${MAKE} pos POS_PREFIX=${SEQ_PREFIX}_${OL}nt${STR}${NOOV}_ci${CI}_mkv${MKV} OPT='-markov ${MKV}'

## Iterate over oligonucleotide lengths
OLIGO_LEN=1 2 3 4 5 6
TASK=pos
all_len_pos:
	@echo
	@echo "Iterating task	${TASK}	over oligo lengths ${OLIGO_LENGTHS}"
	@for ol in ${OLIGO_LEN}; do \
		${MAKE} ${TASK} OL=$${ol} ; \
	done

all_len_mkv:
	${MAKE} all_len_pos TASK=pos_mkv OLIGO_LEN="4 5 6" MKV=2
	${MAKE} all_len_pos TASK=pos_mkv OLIGO_LEN="3 4 5 6" MKV=1
	${MAKE} all_len_pos TASK=pos_mkv OLIGO_LEN="2 3 4 5 6" MKV=0

