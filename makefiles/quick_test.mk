################################################################
## Test the quick option of dyad-analysis.  Check that the slow mode
## and the quick mode return exactly the same results.

include ${RSAT}/makefiles/util.mk
MAKEFILE=${RSAT}/makefiles/quick_test.mk
V=1

################################################################
## Generate a ranom sequence
SEQ_LEN=100
SEQ_NB=100
QUICK_DIR=quick_test
SEQ_PREFIX=rand_L${SEQ_LEN}_N${SEQ_NB}
SEQ_EXT=.fa
SEQ_DIR=${QUICK_DIR}
SEQ_FILE=${SEQ_DIR}/${SEQ_PREFIX}.${SEQ_EXT}

dir:
		@mkdir -p ${QUICK_DIR}

rand: dir
	@echo
	@echo "Generating random sequence"
	@random-seq -l ${SEQ_LEN} -n ${SEQ_NB} -o ${SEQ_FILE}
	@echo ${SEQ_FILE}

################################################################
## Analyze dyads
GROUPING=
COUNT_MODE=
STR=-1str
NOOV=-noov
SP=0-20
ML=3
DYAD_CMD=dyad-analysis -v ${V} -i ${SEQ_FILE} -sp ${SP} -l ${ML} ${STR} ${NOOV} ${GROUPING} -return ${RETURN} ${COUNT_MODE} -o ${DYADS}.tab
DYAD_SUFFIX=dyads_${ML}nt_sp${SP}${STR}${NOOV}${GROUPING}
DYADS=${QUICK_DIR}/${SEQ_PREFIX}_${DYAD_SUFFIX}${COUNT_MODE}
RETURN=occ
dyads: dir
	@echo
	@echo "Running dyad-analysis	${COUNT_MODE}"
	time ${DYAD_CMD}
	@echo "	${DYADS}.tab"

OL=6
OLIGO_CMD=oligo-analysis -v ${V} -i ${SEQ_FILE} -l ${OL} ${STR} ${NOOV} ${GROUPING} -return ${RETURN} ${COUNT_MODE} -o ${OLIGOS}.tab
OLIGO_SUFFIX=oligos_${OL}nt_${STR}${NOOV}${GROUPING}
OLIGOS=${QUICK_DIR}/${SEQ_PREFIX}_${OLIGO_SUFFIX}${COUNT_MODE}
oligos: dir
	@echo
	@echo "Running oligo-analysis	${COUNT_MODE}"
	time ${OLIGO_CMD}
	@echo "	${OLIGOS}.tab"

quick:
	${MAKE} dyads COUNT_MODE=-quick 
	${MAKE} oligos COUNT_MODE=-quick 


count_words: dir
	@echo
	@echo "Running count-words"
	time count-words -v 2 -i ${SEQ_FILE} -l ${ML} ${NOOV} ${STR} ${GROUPING} -sp ${SP} > ${DYADS}_cw.tab
	time count-words -v 2 -i ${SEQ_FILE} -l ${OL} ${NOOV} ${STR} ${GROUPING} > ${OLIGOS}_cw.tab

SC=3
COMPA=${QUICK_DIR}/${SEQ_PREFIX}_${DYAD_SUFFIX}_compa_sc${SC}
compare:
	@echo
	@echo "Comparing quick and slow modes"
	compare-scores -v 1 -basenames -ic 2 -sc ${SC} -i ${DYADS}.tab -i ${DYADS}-quick.tab | awk '{print $$0"\t"($$3-$$2)}'  > ${COMPA}.tab 
	@echo "	${COMPA}.tab"

all: dyads oligos quick count_words compare


################################################################
## A minimal test: small sequence, small oligo, a few spacings
small:
	${MAKE} all SEQ_LEN=3 SEQ_NB=10 ML=1 OL=2 SP=0-1 
