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
DIR=rand_seq
SEQ=${RAND_SEQ}
SEQ_FILE=${RAND_SEQ_FILE}
SEQ_PREFIX=${RAND_SEQ_PREFIX}

STR=-2str
NOOV=-noov
BG=monads
SP=0-20
ML=3
DYAD_SUFFIX=dyads_${ML}nt_sp${SP}${STR}${NOOV}_${BG}
DYADS=${DIR}/${SEQ}_${DYAD_SUFFIX}
DYAD_FILES=`ls ${DIR}/${SEQ_PREFIX}_test*_${DYAD_SUFFIX}.tab`
V=1
DYAD_CMD=dyad-analysis -v ${V}\
	-i ${SEQ_FILE} \
	-sort -type any ${STR} ${NOOV} \
	-org ${ORG} -return occ,freq,rank,proba,monad_freq \
	-l ${ML} -spacing ${SP} -bg ${BG} \
	-o ${DYADS}.tab \
	${OPT}

################################################################
## Apply pattern discovery algorithm to one sequence file
one_disco:
	${DISCO_CMD}
	@echo ${PATTERNS}.tab

RAND_OL=6
SEQ_LEN=1000
SEQ_NB=10
RAND_SEQ_PREFIX=rand_L${SEQ_LEN}_n${SEQ_NB}_bg_${RAND_OL}nt_${ORG}
RAND_SEQ=${RAND_SEQ_PREFIX}_test${TEST}
RAND_SEQ_FILE=${DIR}/${SEQ}.fasta
one_rand_seq: 
	@mkdir -p ${DIR}
	random-seq -org ${ORG} -bg upstream-noorf -ol ${RAND_OL} -l ${SEQ_LEN} -r ${SEQ_NB} -o ${RAND_SEQ_FILE}
	@echo "${RAND_SEQ_FILE}"

one_test: one_rand_seq one_disco

TEST=1
TESTS=1 2 3 4 5 6 7 8 9 10
test_series:
	@for t in ${TESTS} ; do \
		date ; \
		${MAKE} one_test TEST=$${t}; \
	done
	${MAKE} score_distrib

list_disco_files:
	@echo ${DISCO_FILES}

SC=9
SCORE_FILE=${DIR}/${SEQ_PREFIX}_${SUFFIX}_distrib
score_distrib: list_disco_files
	@echo
	@echo "Score distribution"
	cat ${DISCO_FILES} | grep -v "^;" | cut -f ${SC} | classfreq -v -o ${SCORE_FILE}.tab ${OPT}
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

problem:
	dyad-analysis  -i ${SEQ_FILE} -v 1 -sort -timeout 3600 -type any -2str -noov -org Saccharomyces_cerevisiae -lth occ 1 -lth occ_sig 0 -return occ -l 3 -spacing 0-0 -v 2
