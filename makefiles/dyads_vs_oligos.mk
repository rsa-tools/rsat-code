## ##############################################################
## This scripts checks the consistency between olligo-analysis and
## dyad-analysis, by comparing the results obtained for
## hexanucleotides (oligo-analysis), and for pairs of trinucleotides
## with spacing 0 (dyad-analysis).
##
## This script has been useful for debugging prooblems with the
## treatment of N in dyad-analysis.

include ${RSAT}/makefiles/util.mk

all_tests:
	@for str in -1str -2str ; do  					\
		for noov in -noov -ovlp; do 				\
			${MAKE} one_test NOOV=$${noov} STR=$${str} ;	\
		done 							\
	done

one_test: oligos dyads all_compa

V=1 

## Sequence file
#SEQ_FILE=PURA_Escherichia_coli_K12_Gammaproteobacteria_up_purged.fasta
SEQ=ACEB
SEQ_FILE=${SEQ}.fasta

################################################################
## Common parmeters for oligo-analysis and dyad-analysis
## The option zeroocc is used to make sure that all patterns are
## counted, including those which are not found in the sequence set.
NOOV=-noov
STR=-2str
RETURN=-return occ,freq,proba,rank,zscore,ratio
ORG=Escherichia_coli_K12
BG=upstream-noorf
BG_OPT=-org ${ORG} -bg ${BG}
ZEROOCC=-zeroocc
OPTIONS=-v ${V} -i ${SEQ_FILE} ${STR} ${NOOV} ${RETURN} ${BG_OPT} ${ZEROOCC} -sort
SUFFIX=${STR}${NOOV}${ZEROOCC}_${BG}

################################################################
## Run oligo-analysis
OLIGOS=oligos_${SEQ}_${OL}nt${SUFFIX}
OL=6
oligos:
	oligo-analysis -l ${OL} ${OPTIONS} -o ${OLIGOS}.tab
	@echo ${OLIGOS}.tab

## oligo-analysis with a Markov chain model estimated from the input sequence
MKV=2
markov:
	${MAKE} oligos BG='mkv${MKV}' BG_OPT=' -markov ${MKV}'

################################################################
## Run dyad-analysis
DYADS=dyads_${SEQ}_${ML}nt_sp${SP}${SUFFIX}
DYADS_NOSPACING=${DYADS}_nospacing
ML=3
SP=0-20
dyads:
	dyad-analysis -l ${ML} -sp ${SP} -return monad_freq ${OPTIONS} -o ${DYADS}.tab 
	@echo ${DYADS}.tab
	grep 'n{0}' ${DYADS}.tab \
		| perl -pe 's/n\{0\}//g' \
		> ${DYADS_NOSPACING}.tab
	@echo ${DYADS_NOSPACING}.tab

################################################################
## Compare one column of the oligo-analysis and dyad-analysis result
## files
COMPA=compa${SP}_${SEQ}${SUFFIX}_sc${SC}
SC=3
LOG=
DIFF=diff_${COMPA}
one_compa:
	compare-scores -i ${OLIGOS}.tab -i ${DYADS_NOSPACING}.tab -sc ${SC} -o ${COMPA}.tab
	awk -F'\t' '$$2 != $$3 {print $$0"\t"($$3-$$2)}' ${COMPA}.tab > ${DIFF}
	@echo ${COMPA}.tab
	XYgraph -i ${COMPA}.tab \
		-xcol 2 -ycol 3 \
		-format jpg ${LOG} \
		-title1 '${SEQ_FILE}' \
		-title2 'score column ${SC}' \
		-xleg1 'oligos' \
		-yleg1 'dyads' \
		-o ${COMPA}.jpg
	@echo ${COMPA}.jpg

################################################################
## Compare all the comparable columns between oligo-analysis and dyad-analysis
all_compa:
	${MAKE} one_compa SC=3
	${MAKE} one_compa SC=4
	${MAKE} one_compa SC=5
	${MAKE} one_compa SC=6
	${MAKE} one_compa SC=7 LOG='-xlog 10 -ylog 10'
	${MAKE} one_compa SC=8 LOG='-xlog 10 -ylog 10'
	${MAKE} one_compa SC=9
	${MAKE} one_compa SC=10
	${MAKE} one_compa SC=11
	${MAKE} one_compa SC=13
	${MAKE} one_compa SC=14


################################################################
## Check the distribution of different statistics in the pattern
## discovery file.

## Generate random sequences
RAND_SEQ=rand_l${RAND_L}_r${RAND_R}
RAND_L=50
RAND_R=10
RAND_BG=
rand_seq:
	random-seq -l ${RAND_L} -r ${RAND_R} ${RAND_BG} -o ${RAND_SEQ}.fasta
	@echo ${RAND_SEQ}.fasta

## Test the distribution of score in the random sequences
RAND_PARAM=SEQ='${RAND_SEQ}' BG='monad' BG_OPT=''
rand_dyads:
	${MAKE} dyads ${RAND_PARAM}

rand_dyad_distrib:
	${MAKE} score_distrib ${RAND_PARAM}

DISCO_FILE=${DYADS}
score_distrib:
	grep -v "^;" ${DISCO_FILE}.tab | cut -f ${SC} | classfreq -v -o ${DISCO_FILE}_sc${SC}_distrib.tab