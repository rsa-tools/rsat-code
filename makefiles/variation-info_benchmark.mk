## Load site-specific options for the cluster + other parameters
include ${RSAT}/RSAT_config.mk

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
N=5000
CHUNK=1000
IN_PREFIX=diabetes_DA-LD
QUERIES=${IN_PREFIX}_SNP-IDs.txt
TOP_QUERIES=${IN_PREFIX}_top${N}_SNP-IDs.txt
OUT_PREFIX=`hostname`_${IN_PREFIX}_top${N}_chunk${CHUNK}_SNP-info
VARBED=${OUT_PREFIX}.varBed
OUT=${OUT_PREFIX}_out.txt
LOG=${OUT_PREFIX}_log.txt

param:
	@echo "Parameters"
	@echo "	N		${N}"
	@echo "	CHUNK		${CHUNK}"
	@echo "	QUERIES		${QUERIES}"
	@echo "	TOP_QUERIES	${TOP_QUERIES}"
	@echo "	VARBED		${VARBED}"
	@echo "	OUT		${OUT}"
	@echo "	LOG		${LOG}"

one_test: param
	@head -n ${N} ${QUERIES} > ${TOP_QUERIES};
	(time variation-info -v ${V} \
		-species Homo_sapiens -assembly GRCh38  \
		-format id \
		-chunk_size ${CHUNK} \
		-i ${TOP_QUERIES} \
		-o ${VARBED} > ${OUT}) 2> ${LOG} &
#	@${MAKE} param

some_tests:
	@${MAKE} one_test N=10 CHUNK=01
	@${MAKE} one_test N=10 CHUNK=02
	@${MAKE} one_test N=10 CHUNK=05
	@${MAKE} one_test N=10 CHUNK=10

	@${MAKE} one_test N=100 CHUNK=010
	@${MAKE} one_test N=100 CHUNK=020
	@${MAKE} one_test N=100 CHUNK=050
	@${MAKE} one_test N=100 CHUNK=100

	@${MAKE} one_test N=1000 CHUNK=0100
	@${MAKE} one_test N=1000 CHUNK=0200
	@${MAKE} one_test N=1000 CHUNK=0500
	@${MAKE} one_test N=1000 CHUNK=1000

	@${MAKE} one_test N=5000 CHUNK=0100
	@${MAKE} one_test N=5000 CHUNK=0200
	@${MAKE} one_test N=5000 CHUNK=0500
	@${MAKE} one_test N=5000 CHUNK=0600
	@${MAKE} one_test N=5000 CHUNK=0700
	@${MAKE} one_test N=5000 CHUNK=0800
	@${MAKE} one_test N=5000 CHUNK=0900
	@${MAKE} one_test N=5000 CHUNK=1000
	@${MAKE} one_test N=5000 CHUNK=1100
	@${MAKE} one_test N=5000 CHUNK=1200
	@${MAKE} one_test N=5000 CHUNK=1300
	@${MAKE} one_test N=5000 CHUNK=1400
	@${MAKE} one_test N=5000 CHUNK=1500
	@${MAKE} one_test N=5000 CHUNK=2000
	@${MAKE} one_test N=5000 CHUNK=5000

	@${MAKE} one_test N=10000 CHUNK=1000
	@${MAKE} one_test N=10000 CHUNK=2000
	@${MAKE} one_test N=10000 CHUNK=5000

	@${MAKE} one_test N=15000 CHUNK=1000
	@${MAKE} one_test N=15000 CHUNK=2000
	@${MAKE} one_test N=15000 CHUNK=3000
	@${MAKE} one_test N=15000 CHUNK=5000


# CHUNK_SERIES=0100 0200 0300 0400 0500 0600 0700 0800 0900 1000 1100 1200 1300 1400 1500 1600 1700 1800 1900 2000
CHUNK_SERIES=0100 0200 0400 0600 0800 1000 1200 1400 1600 1800 2000
one_series:
	@for c in ${CHUNK_SERIES}; do \
		echo "	one_series	N=${N}	CHUNK=${CHUNK}"; \
		${MAKE} one_test CHUNK=$${c} ; \
	done

get_time:
	grep user `hostname`_*log.txt


one_series_old:
	@${MAKE} one_test N=5000 CHUNK=0100
	@${MAKE} one_test N=5000 CHUNK=0200
	@${MAKE} one_test N=5000 CHUNK=0300
	@${MAKE} one_test N=5000 CHUNK=0400
	@${MAKE} one_test N=5000 CHUNK=0500
	@${MAKE} one_test N=5000 CHUNK=0600
	@${MAKE} one_test N=5000 CHUNK=0700
	@${MAKE} one_test N=5000 CHUNK=0800
	@${MAKE} one_test N=5000 CHUNK=0900
	@${MAKE} one_test N=5000 CHUNK=1000
	@${MAKE} one_test N=5000 CHUNK=1100
	@${MAKE} one_test N=5000 CHUNK=1200
	@${MAKE} one_test N=5000 CHUNK=1300
	@${MAKE} one_test N=5000 CHUNK=1400
	@${MAKE} one_test N=5000 CHUNK=1500
	@${MAKE} one_test N=5000 CHUNK=1600
	@${MAKE} one_test N=5000 CHUNK=1700
	@${MAKE} one_test N=5000 CHUNK=1800
	@${MAKE} one_test N=5000 CHUNK=1900
	@${MAKE} one_test N=5000 CHUNK=2000
