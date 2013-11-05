################################################################
## Calculate oligonucleotide distributions in the whole set of
## upstream sequences of a given organism

include ${RSAT}/makefiles/util.mk
MAKEFILE=${RSAT}/makefiles/calibrate_oligos.mk

STR=-2str
NOOV=-noov
UP_LEN=1000
OL=6

################################################################
## Oligo distributions in R random selections of N genes
R=1000
PWD=`pwd`
CALIB_DIR=${PWD}/results/${ORG}/rand_gene_selections/${OL}nt${STR}${NOOV}_N${N}_L${UP_LEN}_R${R}
CALIB_CMD=calibrate-oligos  -v 1 -r 1000 -sn ${N} -ol ${OL}			\
		-org ${ORG}							\
		-sl ${UP_LEN} -task all,clean_oligos -start 1 ${STR} ${NOOV}	\
		-outdir ${CALIB_DIR} 
calibrateN:
	@${MAKE} my_command MY_COMMAND="${CALIB_CMD}" JOB_PREFIX='calib.${ORG}.${OL}nt${STR}${NOOV}_N${N}_L${UP_LEN}_R${R}'
#	@${MAKE} multi MULTI_TASK=calibrate OPT='-bg calibN  -last 5'

################################################################
## Calculate oligo distributions for random selections of various sizes
N_VALUES=1 2 3 4 5 10 20 30 40 50 100 150 200
n_series:
	@for n in ${N_VALUES}; do \
		${MAKE} calibrateN N=$${n}; \
	done

## Tun the N series in a batch queue on the PC cluster
batch_n_series:
	${MAKE} n_series WHEN=queue