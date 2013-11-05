################################################################
#
# Check random models by generating random sequences and counting
# oligonucleotide occurrences

MAKEFILE=${RSAT}/makefiles/random_models.mk
include ${RSAT}/makefiles/util.mk

DATE=`date +%Y%m%d_%H%M%S`

all: random_seq all_oligos

REPET=10000
SEQ_LEN=1000
DATA=data
SEQ_PREFIX=rand_L${SEQ_LEN}_R${REPET}_${MODEL_SUFFIX}
SEQ_FILE=${DATA}/${SEQ_PREFIX}.fta.gz
#MODEL=iid
random_seq:
	${MAKE} random_seq_${MODEL}

################################################################
# Generate a random sequence with a Bernouilli process and equiprobable nucleotides
random_seq_iid: 
	@mkdir -p ${DATA}
	random-seq -l ${SEQ_LEN} -r ${REPET} -o ${SEQ_FILE}

################################################################
# Generate a random sequence with a Markov chain process
MKV=3
MKV_OL=4
MODEL=markov
MODEL_SUFFIX=mkv${MKV}
ORG=Saccharomyces_cerevisiae
RANDOM_SEQ_CMD=random-seq -l ${SEQ_LEN} -r ${REPET} -o ${SEQ_FILE} -bg upstream -org ${ORG} -ol ${MKV_OL}
random_seq_markov:
	@mkdir -p ${DATA}
	@echo ${RANDOM_SEQ_CMD}
	${RANDOM_SEQ_CMD}

all_tests:
	@${MAKE} all_oligos NOOV=-noov
	@${MAKE} all MKV=0 MKV_OL=1
	@${MAKE} all MKV=5 MKV_OL=6
	@${MAKE} all_oligos 

################################################################
# Count oligonucleotide occurrences

OLIGO_LENGTHS=1 2 3 4 
#OLIGO_LENGTHS=1 2 3 4 5 6 7 8
all_oligos:
	@for ol in ${OLIGO_LENGTHS} ; do		\
		${MAKE} count_oligos OL=$${ol} ;	\
	done

OL=4
STR=-1str
NOOV=-ovlp
OLIGO_DIR=results/oligos
OLIGO_FILE=${OLIGO_DIR}/${SEQ_PREFIX}_${OL}nt${STR}${NOOV}.tab
OLIGO_CMD=oligo-analysis -v 1 -i ${SEQ_FILE}	\
		-l ${OL} -return occ -table	\
		${STR} ${NOOV}			\
		| perl -pe 's/^; seq/seq/'	\
		> ${OLIGO_FILE}
count_oligos:
	@mkdir -p ${OLIGO_DIR}
	@echo "${DATE}	${OLIGO_CMD}"
	${OLIGO_CMD}

