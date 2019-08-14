## Load site-specific options for this RSAT instance
include ${RSAT}/makefiles/util.mk

################################################################
## Variables
V=1
MAKEFILE=${RSAT}/makefiles/variation-info_benchmark.mk
MAKE=make -s -f ${MAKEFILE}
DATE=`date +%Y-%M-%d_%H:%M:%S`
DAY=`date +%Y%m%d`
TIME=`date +%Y%m%d_%H%M%S`
SSH_OPT = -e ssh 
RSYNC_OPT= -ruptvlz  ${SSH_OPT} 
RSYNC = rsync  ${RSYNC_OPT}
WGET=wget --passive-ftp -np -rNL

################################################################
## List of targets
usage:
	@echo "usage: make [-OPT='options'] target"
	@echo "implemented targets"
	@perl -ne 'if (/^([a-z]\S+):/){ print "\t$$1\n";  }' ${MAKEFILE}

V=2
N=2000
CHUNK=500
IN_PREFIX=diabetes_DA-LD
QUERIES=${IN_PREFIX}_SNP-IDs.txt
TOP_QUERIES=${IN_PREFIX}_top${N}_SNP-IDs.txt
OUT_PREFIX=`hostname`_${IN_PREFIX}_top${N}
VARINFO_PREFIX=${OUT_PREFIX}_chunk${CHUNK}_SNP-info
VARBED=${VARINFO_PREFIX}.varBed
VARINFO_OUT=${VARINFO_PREFIX}_out.txt
VARINFO_LOG=${VARINFO_PREFIX}_log.txt
VARSEQ=${OUT_PREFIX}.varSeq
VARSEQ_OUT=${OUT_PREFIX}_varseq_out.txt
VARSEQ_LOG=${OUT_PREFIX}_varseq_log.txt

################################################################
##LIst parameter values
param:
	@echo "Parameters"
	@echo "	N		${N}"
	@echo "	CHUNK		${CHUNK}"
	@echo "	QUERIES		${QUERIES}"
	@echo "	TOP_QUERIES	${TOP_QUERIES}"
	@echo "	VARBED		${VARBED}"
	@echo "	VARINFO_OUT	${VARINFO_OUT}"
	@echo "	VARINFO_LOG	${VARINFO_LOG}"
	@echo "	VARSEQ		${VARSEQ}"
	@echo "	VARSEQ_OUT	${VARSEQ_OUT}"
	@echo "	VARSEQ_LOG	${VARSEQ_LOG}"


################################################################
## Run one benchmark test for variation-info with a given number of
## query IDs and a given chunk size
VARINFO_CMD=(time variation-info -v ${V} \
	-species Homo_sapiens -assembly GRCh38  \
	-format id \
	-chunk_size ${CHUNK} \
	-i ${TOP_QUERIES} \
	-o ${VARBED} > ${VARINFO_OUT}) 2> ${VARINFO_LOG} &
varinfo_one_test: param
	@head -n ${N} ${QUERIES} > ${TOP_QUERIES};
	${VARINFO_CMD}

################################################################
## Run one benchmark test for variation-info with a given number of
## query IDs and a given chunk size
VARSEQ_CMD=(time retrieve-variation-seq -v ${V} \
	-species Homo_sapiens -assembly GRCh38  \
	-format varBed \
	-mml 30 \
	-i ${VARBED} \
	-o ${VARSEQ} > ${VARSEQ_OUT}) 2> ${VARSEQ_LOG} &
varseq_one_test: param
	${VARSEQ_CMD}

################################################################
## Test chunk values in a given range
# CHUNK_SERIES=0100 0200 0300 0400 0500 0600 0700 0800 0900 1000 1100 1200 1300 1400 1500 1600 1700 1800 1900 2000
CHUNK_SERIES=0100 0200 0300 0400 0500 0600 0700 0800 0900 1000 1100 1200 1300 1400 1500 1600 1800 2000
one_chunk_series:
	@for c in ${CHUNK_SERIES}; do \
		echo "	${TEST}	one_chunk_series	N=${N}	CHUNK=$${c}"; \
		${MAKE} ${TEST}_one_test CHUNK=$${c} ; \
	done

################################################################
## Test N values for a given chunk value
N_SERIES=01000 02000 03000 04000 05000 06000 08000 10000 12000 15000 18000
TEST=varinfo
one_n_series:
	@for n in ${N_SERIES}; do \
		echo "	${TEST}	one_n_series	N=$${n}	CHUNK=${CHUNK}"; \
		${MAKE} ${TEST}_one_test N=$${n} ; \
	done

################################################################
## Run tests with indicative values for the number of queries and chunk size
done_one_by_one:
	@${MAKE} TEST=varinfo one_chunk_series N=10 CHUNK_SERIES='01 02 05 10'
	@${MAKE} TEST=varinfo one_chunk_series N=100 CHUNK_SERIES='010 020 030 040 050 060 070 080 090 100'
	@${MAKE} TEST=varinfo one_chunk_series N=1000 CHUNK_SERIES='0100 0200 0300 0400 0500 0600 0700 0800 0900 1000'
	@${MAKE} TEST=varinfo one_chunk_series N=2000 CHUNK_SERIES='0100 0200 0300 0400 0500 0600 0700 0800 0900 1000 1100 1200 1300 1400 1500 1600 1800 2000'
	@${MAKE} TEST=varinfo one_chunk_series N=3000 CHUNK_SERIES='0100 0200 0300 0400 0500 0600 0700 0800 0900 1000 1100 1200 1300 1400 1500 1600 1800 2000'
	@${MAKE} TEST=varinfo one_chunk_series N=4000 CHUNK_SERIES='0100 0200 0300 0400 0500 0600 0700 0800 0900 1000 1100 1200 1300 1400 1500 1600 1800 2000'
	@${MAKE} TEST=varinfo one_chunk_series N=10000 CHUNK_SERIES='0100 0200 0300 0400 0500 0600 0700 0800 0900 1000 1100 1200 1300 1400 1500 1600 1800 2000'
	@${MAKE} TEST=varinfo one_chunk_series N=18000 CHUNK_SERIES='0100 0200 0300 0400 0500 0600 0700 0800 0900 1000 1100 1200 1300 1400 1500 1600 1800 2000'
	@${MAKE} TEST=varinfo one_chunk_series N=5000 CHUNK_SERIES='0100 0200 0300 0400 0500 0600 0700 0800 0900 1000 1100 1200 1300 1400 1500 1600 1800 2000'
	@${MAKE} TEST=varinfo one_chunk_series N=15000 CHUNK_SERIES='0100 0200 0300 0400 0500 0600 0700 0800 0900 1000 1100 1200 1300 1400 1500 1600 1800 2000'
	@${MAKE} TEST=varinfo one_n_series CHUNK=0500
	@${MAKE} TEST=varinfo one_n_series CHUNK=1000
	@${MAKE} TEST=varinfo one_n_series CHUNK=2000
	@${MAKE} TEST=varinfo one_n_series CHUNK=3000

next:
	@${MAKE}  TEST=varseq one_n_series CHUNK=0500

to_do_one_by_one:
	@${MAKE} TEST=varinfo one_chunk_series N=6000 CHUNK_SERIES='0100 0200 0300 0400 0500 0600 0700 0800 0900 1000 1100 1200 1300 1400 1500 1600 1800 2000'
	@${MAKE} TEST=varinfo varinfo_one_test N=5000 CHUNK=5000
	@${MAKE} TEST=varinfo varinfo_one_test N=10000 CHUNK=5000
	@${MAKE} TEST=varinfo varinfo_one_test N=15000 CHUNK=3000
	@${MAKE} TEST=varinfo varinfo_one_test N=15000 CHUNK=5000
	@${MAKE} TEST=varinfo varinfo_one_test N=18000 CHUNK=3000
	@${MAKE} TEST=varinfo varinfo_one_test N=18000 CHUNK=5000

################################################################
## Collect user time for the different tests
WHICH_TIME=user
get_time:
	grep ${WHICH_TIME} `hostname`_*${TEST}_log.txt
