################################################################
## Test the quick option of dyad-analysis.  Check that the slow mode
## and the quick mode return exactly the same results.

include ${RSAT}/makefiles/util.mk
MAKEFILE=${RSAT}/makefiles/quick_test.mk
V=1

################################################################
## Generate a ranom sequence
L=100
N=100
QUICK_DIR=quick_test
SEQ_PREFIX=rand_L${L}_N${N}
SEQ_EXT=fa
SEQ_DIR=${QUICK_DIR}
SEQ_FILE=${SEQ_DIR}/${SEQ_PREFIX}.${SEQ_EXT}

dir:
		@mkdir -p ${QUICK_DIR}

rand: dir
	@echo
	@echo "Generating random sequence"
	@random-seq -l ${L} -n ${N} -o ${SEQ_FILE}
	@echo ${SEQ_FILE}

################################################################
## Analyze dyads
GROUPING=
COUNT_MODE=
STR=-1str
NOOV=-noov
SP=0-20
ML=3
DYAD_CMD=dyad-analysis -v ${V} -i ${SEQ_FILE} -sp ${SP} -l ${ML} ${STR} ${NOOV} ${GROUPING} -return ${RETURN} ${COUNT_MODE} -o ${DYADS}.tab ${OPT}
DYAD_SUFFIX=dyads_${ML}nt_sp${SP}${STR}${NOOV}${GROUPING}
DYADS=${QUICK_DIR}/${SEQ_PREFIX}_${DYAD_SUFFIX}${COUNT_MODE}
RETURN=occ,freq,ratio,zscore,proba,rank
dyads: dir
	@echo
	@echo "Running dyad-analysis	${COUNT_MODE}"
	time ${DYAD_CMD}
	@echo "	${DYADS}.tab"

dyads_quick:
	${MAKE} dyads COUNT_MODE=-quick 

OL=6
OLIGO_CMD=oligo-analysis -v ${V} -i ${SEQ_FILE} -l ${OL} ${STR} ${NOOV} ${GROUPING} -return ${RETURN} ${COUNT_MODE} -o ${OLIGOS}.tab ${OPT}
OLIGO_SUFFIX=oligos_${OL}nt_${STR}${NOOV}${GROUPING}
OLIGOS=${QUICK_DIR}/${SEQ_PREFIX}_${OLIGO_SUFFIX}${COUNT_MODE}
oligos: dir
	@echo
	@echo "Running oligo-analysis	${COUNT_MODE}"
	time ${OLIGO_CMD}
	@echo "	${OLIGOS}.tab"

oligos_quick:
	${MAKE} oligos COUNT_MODE=-quick 

################################################################
## Directly run count-words on the command line
count_words: dir
	@echo
	@echo "Running count-words (dyads)"
	time count-words -v 1 -i ${SEQ_FILE} -l ${ML} ${NOOV} ${STR} ${GROUPING} -sp ${SP} > ${DYADS}_cw.tab
	@wc -l ${DYADS}_cw.tab
	@echo
	@echo "Running count-words (oligos)"
	time count-words -v 1 -i ${SEQ_FILE} -l ${OL} ${NOOV} ${STR} ${GROUPING} > ${OLIGOS}_cw.tab
	@wc -l ${OLIGOS}_cw.tab

################################################################
## Run count-words with the input as STDIN
## Compare results with STDIN and option  -i.
##
## Usage: make -i -f $RSAT/makefiles/quick_test.mk count_words_stdin
##
## The option -i is required because the command diff returns an error
## on the STDER when there is a diff between two files (and there is
## always a diff, at least in the command line)s
count_words_stdin: count_words
	@echo
	@echo "Running count-words STDIN (dyads)"
	cat ${SEQ_FILE}  | count-words -v 1 -l ${ML} ${NOOV} ${STR} ${GROUPING} -sp ${SP} > ${DYADS}_cw_stdin.tab
	@echo "Line numbers for dyad files" 
	@wc -l ${DYADS}_cw.tab ${DYADS}_cw_stdin.tab
	@echo
	@echo "Running count-words STDIN (oligos)"
	time cat  ${SEQ_FILE} | count-words -v 1 -l ${OL} ${NOOV} ${STR} ${GROUPING} > ${OLIGOS}_cw_stdin.tab
	@echo "Line numbers for oligo files" 
	@wc -l ${OLIGOS}_cw.tab ${OLIGOS}_cw_stdin.tab
	@diff ${DYADS}_cw.tab ${DYADS}_cw_stdin.tab > ${DYADS}_cw_sdtin_diff.txt
	@diff ${OLIGOS}_cw.tab ${OLIGOS}_cw_stdin.tab > ${OLIGOS}_cw_sdtin_diff.txt

################################################################
## Compare the results of dyad-analysis in slow and quick mode
SC=3
COMPA=${QUICK_DIR}/${SEQ_PREFIX}_${DYAD_SUFFIX}_compa_sc${SC}
compare:
	@echo
	@echo "Comparing quick and slow modes"
	compare-scores -v 1 -basenames -ic 2 -sc ${SC} -i ${DYADS}.tab -i ${DYADS}-quick.tab | awk '{print $$0"\t"($$3-$$2)}'  > ${COMPA}.tab 
	@echo "	${COMPA}.tab"

all: dyads dyads_quick oligos oligos_quick count_words compare


################################################################
## A minimal test: small sequence, small oligo, a few spacings
SMALL_TASK=rand all
small:
	${MAKE} ${SMALL_TASK} L=3 N=10 ML=1 OL=2 SP=0-1 
