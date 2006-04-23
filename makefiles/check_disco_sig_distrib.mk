################################################################
## Check the distibution of significance returned by string-based
## pattern discovery programs, by analyzing random sequences or random
## selections of promoters.

include ${RSAT}/makefiles/util.mk


## Choice of the pattern discovery program
DISCO_CMD=${DYAD_CMD}
PATTERNS=${DYADS}
SUFFIX=${DYAD_SUFFIX}
DISCO_FILES=${DYAD_FILES}

## Model organism
ORG=Saccharomyces_cerevisiae

## Directory and sequences
WD=`pwd`
DIR=${RAND_DIR}
SEQ=${RAND_SEQ}
SEQ_FILE=${RAND_SEQ_FILE}
SEQ_PREFIX=${RAND_SEQ_PREFIX}


################################################################
## Apply pattern discovery algorithm to one sequence file

## dyad-analysis command
STR=-2str
NOOV=-noov
BG=monads
SP=0-20
ML=3
DYAD_SUFFIX=dyads_${ML}nt_sp${SP}${STR}${NOOV}_${BG}
DYADS=${DIR}/${SEQ}_${DYAD_SUFFIX}
DYAD_FILES=`ls ${RAND_DIR}/${SEQ_PREFIX}_test*_${DYAD_SUFFIX}.tab*`
V=1
DYAD_CMD=dyad-analysis -v ${V}\
	-i ${SEQ_FILE} \
	-sort -type any ${STR} ${NOOV} \
	-org ${ORG} -return occ,freq,rank,proba,monad_freq \
	-l ${ML} -spacing ${SP} -bg ${BG} \
	-o ${DYADS}.tab.gz \
	${OPT}

## Run pattern discovery on one sequence set
one_disco:
	@echo
	@echo ${DISCO_CMD}
	make my_command MY_COMMAND="${DISCO_CMD}"
	@echo ${PATTERNS}.tab

################################################################
## Generate one set of random sequence 
RAND_OL=6
L=200
N=20
BG_DESC=bg_${RAND_BG}_${RAND_OL}nt_${ORG}
RAND_DIR=${WD}/results/rand_seq_${BG_DESC}_n${N}_l${L}
RAND_SEQ_PREFIX=rand_L${L}_n${N}_${BG_DESC}
RAND_SEQ=${RAND_SEQ_PREFIX}_test${TEST}
RAND_SEQ_FILE=${DIR}/${SEQ}.fasta.gz
RAND_BG=equi
RAND_SEQ_CMD=mkdir -p ${RAND_DIR}; random-seq -org ${ORG} -bg ${RAND_BG} -ol ${RAND_OL} -l ${L} -r ${N} -o ${RAND_SEQ_FILE}
one_rand_seq: 
	@echo
	@echo "${RAND_SEQ_CMD}"
	${RAND_SEQ_CMD}
	@echo "${RAND_SEQ_FILE}"


ONE_TEST_CMD=${RAND_SEQ_CMD}; ${DISCO_CMD}
one_test: 
	@echo
#	@echo "${ONE_TEST_CMD}"
	@${MAKE} my_command MY_COMMAND="${ONE_TEST_CMD}"
	@echo "${RAND_SEQ_FILE}"
	@echo ${PATTERNS}.tab

################################################################
## Run a series of tests
TEST=1
TESTS_10=1 2 3 4 5 6 7 8 9 10
TESTS_100=0 1 2 3 4 5 6 7 8 9 \
	10 11 12 13 14 15 16 17 18 19 \
	 20 21 22 23 24 25 26 27 28 29 \
	 30 31 32 33 34 35 36 37 38 39 \
	 40 41 42 43 44 45 46 47 48 49 \
	 50 51 52 53 54 55 56 57 58 59 \
	 60 61 62 63 64 65 66 67 68 69 \
	 70 71 72 73 74 75 76 77 78 79 \
	 80 81 82 83 84 85 86 87 88 89 \
	 90 91 92 93 94 95 96 97 98 99
TESTS=${TESTS_100}
test_series:
	@for t in ${TESTS} ; do \
		date ; \
		${MAKE} one_test TEST=$${t}; \
	done


list_disco_files:
	@echo ${WD}
	@echo ${DIR}
	@echo ${DISCO_FILES}

################################################################
## Calculate score distribution and draw the graph
SC=7
SCORE_FILE=${DIR}/${SEQ_PREFIX}_${SUFFIX}_distrib
score_distrib: list_disco_files
	@echo
	@echo "Score distribution"
	zcat ${DISCO_FILES} | grep -v "^;" | cut -f ${SC} | classfreq -v -o ${SCORE_FILE}.tab ${OPT}
	@echo ${SCORE_FILE}.tab 
	XYgraph	-i ${SCORE_FILE}.tab \
		-lines \
		-xcol 3 -ycol 4,5,6 -legend \
		-xleg1 'score column ${SC}' \
		-yleg1 'Absolute frequency' -ylog \
		-title1 '${SEQ_PREFIX}' \
		-title2 '${SUFFIX}' \
		-o ${SCORE_FILE}.jpg

	@echo ${SCORE_FILE}.jpg

################################################################
## Remove random sequences
rm_all_rand_seq:
	@for t in ${TESTS} ; do \
		${MAKE} rm_one_rand_seq TEST=$${t}; \
	done

rm_one_rand_seq:
	@echo "deleting file	${RAND_SEQ_FILE}"
	@rm -f "${RAND_SEQ_FILE}"


################################################################
## Synchronize result files from merlin
## Exclude sequence and dyad files, only synchronize the distributions
from_merlin:
	rsync --exclude '*fasta*' --exclude '*.tab.gz' -ruptvl -e ssh merlin.scmbb.ulb.ac.be:test/dyad_sig_distrib/results .
