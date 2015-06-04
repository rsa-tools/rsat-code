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

include ${RSAT}/makefiles/util.mk


MAKEFILE=${RSAT}/makefiles/downstream_signals.mk

ORG=Saccharomyces_cerevisiae_R64-1-1
RES_DIR=results/downstream_signals
ORG_DIR=${RES_DIR}/${ORG}

dir:
	@mkdir -p ${ORG_DIR}

################################################################
## Collect all sequences around transcription termination site (TTS).
TTS_SEQ=${ORG_DIR}/${ORG}_downstream_tts_from-${UP_TTS}_to${DOWN_TTS}.fasta.gz
UP_TTS=200
DOWN_TTS=100
tts_seq: dir
	@echo
	@echo "Retrieving sequences around the TTS for organism	${ORG}"
	retrieve-seq  -org ${ORG} -feattype mRNA -type downstream -format fasta \
		-label organism_name,id,name,sequence_type,current_from,current_to,ctg,orf_strand,reg_left,reg_right \
		-from -${UP_TTS} -to +${DOWN_TTS} -all \
		-o ${TTS_SEQ}
	@echo "	${TTS_SEQ}"


################################################################
## Run position-analysis to detect positionally biased
## oligonucleotides around the TTS.
OL=6
CI=20
NOOV=-noov
STR=-1str
POS_FILE=${ORG_DIR}/${ORG}_downstream_tts_from-${UP_TTS}_to${DOWN_TTS}_${OL}nt${STR}${NOOV}_ci${CI}
OPT_POS_TO_MATRICES=-return clusters,matrices -clust_nb 8 -max_asmb_nb 5 \
pos:
	@echo
	@echo "Position-analysis with ${OL}nt	${ORG}"
	position-analysis -i ${TTS_SEQ} -v ${V} -sort -nofilter -lth_occ 1 -max_graphs 100 \
		-return chi,sig,rank,distrib,exp_occ,graphs,index ${STR} ${NOOV} -l ${OL} -ci ${CI} \
		-origin start -offset ${UP_TTS} \
		${OPT} \
		-o ${POS_FILE}.tab
	@echo "	${POS_FILE}.tab"

OLIGO_LEN=1 2 3 4 5 6
pos_all_len:
	@echo
	@echo "Iterating over oligonucleotide lengths ${OLIGO_LENGTHS}"
	@for ol in ${OLIGO_LEN}; do \
		${MAKE} pos OL=$${ol} ; \
	done
