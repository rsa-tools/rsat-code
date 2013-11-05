.Suffixes: _distrib.tab _poisson.tab _negbin.tab

################################################################
## makefile for getting info for the motif discovery competition

include ${RSAT}/makefiles/util.mk
include ${RSAT}/makefiles/upstream_calibrations.mk

MAKEFILE=${RSAT}/makefiles/motif_disco_competition_2003.mk

################################################################
#### General parameters

## level of verbosity
# V=1

# ################################################################
# ### commands
# MAKE=make -s -f ${MAKEFILE}

# SSH_OPT = -e ssh 
# RSYNC_OPT= -ruptvlz  ${SSH_OPT} 
# RSYNC = rsync  ${RSYNC_OPT}

# WGET=wget --passive-ftp -np -rNL

# ################################################################
# #### list of targets
# usage:
# 	@echo "usage: make [-OPT='options'] target"
# 	@echo "implemented targets"
# 	@perl -ne 'if (/^([a-z]\S+):/){ print "\t$$1\n";  }' ${MAKEFILE}

################################################################
#### Synchronization between different machines

#### rsync
TO_SYNC=results
#SERVER_LOCATION=rubens.ulb.ac.be:/rubens/dsk3/genomics/motif_discovery_competition/
USER=jvanheld
SERVER_LOCATION=${USER}@${SERVER}:${SERVER_DIR}
SERVER_DIR=motif_discovery_competition_2003/
#SERVER=merlin.ulb.ac.be
SERVER=164.15.109.32
EXCLUDE=--exclude '*~' --exclude oligos --exclude '*.wc' --exclude random_genes.tab --exclude '*.fasta' --exclude '*.fasta.gz' --exclude '*.wc.gz'
all_to_merlin:
	${MAKE} to_merlin TO_SYNC='*.html' 
	${MAKE} to_merlin TO_SYNC='makefiles' 

to_merlin:
	${MAKE} one_dir_to_server SERVER=164.15.109.32

one_dir_to_server:
	${RSYNC} ${EXCLUDE} ${TO_SYNC} ${SERVER_LOCATION}

from_server:
	${RSYNC} ${EXCLUDE} ${SERVER_LOCATION}${TO_SYNC} .

## Synchronize calibrations from merlin
from_merlin:
	${MAKE} from_server SERVER_DIR=motif_discovery_competition_2003/ 

# Temporary
#	${MAKE} from_server SERVER_DIR=./ TO_SYNC=results
#	${MAKE} from_server SERVER_DIR=makefiles/ TO_SYNC=results
#	${MAKE} from_server TO_SYNC='*.xls'
#	${MAKE} from_server TO_SYNC='*.tab'
#	${MAKE} from_server TO_SYNC=my_calibrate-oligos.R


from_liv:
	${MAKE} from_server SERVER=liv.bmc.uu.se ${SERVER_DIR}=motif_discovery_competition_2003/

################################################################
#### retrieve information from the competition web site 
INFO_SITE=http://www.cs.washington.edu/homes/tompa/competition/
get_info:
	${WGET} ${INFO_SITE}

################################################################
#### uncompress the data
DATA=data
DATA_ORI=www.cs.washington.edu/homes/tompa/competition/data.zip
uncompress_data:
	@mkdir -p ${DATA}
	unzip -d ${DATA} ${DATA_ORI}
	mv -f ${DATA}/data\[1\]/* ${DATA}
	rmdir ${DATA}/data\[1\] 
	(cd data ;					\
	mv -f S.cerevisiae Saccharomyces_cerevisiae;	\
	mv -f H.sapiens Homo_sapiens;			\
	mv -f M.musculus Mus_musculus;			\
	mv -f D.melanogaster Drosophila_melanogaster)

################################################################
#### parameters for analyzing one data set
ORG=Saccharomyces_cerevisiae
CURRENT_SET=yst01
SEQ_FILE=${DATA}/${ORG}/${CURRENT_SET}.fasta

## Directories
WORK_DIR=`pwd`
RES_DIR=${WORK_DIR}/results
ORG_DIR=${RES_DIR}/${ORG}
SEQ_LEN_DIR=${ORG_DIR}/sequence_lengths

################################################################
#### Create the result directory for the current set
current_dir:
	@mkdir -p ${ORG_DIR}

ORGANISMS=Saccharomyces_cerevisiae Drosophila_melanogaster Mus_musculus Homo_sapiens
ORG_TASK=calib_script

################################################################
#### iterate over all data sets for a given organism
ORG_SETS=`ls -1 ${DATA}/${ORG}/*.fasta | grep -v purged | perl -pe 's|${DATA}/${ORG}/||g'  | perl -pe 's|\.fasta||g'`
TASK=seq_len
iterate_sets:
	@echo "iterating over data sets for organism ${ORG}"
	@echo ${ORG_SETS}
	@for set in ${ORG_SETS}; do			\
		${MAKE} ${TASK} CURRENT_SET=$${set} ;	\
	done

################################################################
## Iterate some task over all oligonucleotide lengths
OL_TASK=background
OLIGO_LENGTHS=5 6 7 8
iterate_oligo_lengths:
	@echo "iterating over oligo lengths for organism ${ORG}"
	@echo ${OLIGO_LENGTHS}
	@for ol in ${OLIGO_LENGTHS}; do			\
		echo "oligo length $${ol}" ;		\
		${MAKE} ${OL_TASK} OLIGO_LEN=$${ol} ;	\
	done

################################################################
#### Calculate background model, depending on the size of the sequences in the current set
UPSTREAML_DIR=${ORG_DIR}/upstreamL_frequencies
background: seq_lengths_one_org all_up bg_oligos

################################################################
#### report sequence length for the current set
CURRENT_SEQ_LEN=`sequence-lengths -i ${SEQ_FILE} | cut -f 2 | sort -u`
SEQ_NB=`sequence-lengths -i ${SEQ_FILE} | wc -l | awk '{print $$1}'`
SEQ_LEN_FILE=${SEQ_LEN_DIR}/${CURRENT_SET}_lengths.txt
seq_len:  current_dir
	@mkdir -p ${SEQ_LEN_DIR}
	@echo "set ${CURRENT_SET} length ${CURRENT_SEQ_LEN}	${SEQ_NB} seqs	${SEQ_LEN_FILE}"
	@echo ${CURRENT_SEQ_LEN} > ${SEQ_LEN_FILE}

seq_lengths_one_org:
	${MAKE} iterate_sets TASK=seq_len

seq_lengths:
	${MAKE} iterate_organisms ORG_TASK=seq_lengths_one_org TASK=seq_len

################################################################
## Generate organism-specific makefiles for calibrating pattern
## frequencies in upstream sequences of the same sizes as in the test
## sets
MAKE_DIR=${WORK_DIR}/makefiles
CALIB_SCRIPT=${MAKE_DIR}/calibrate_${ORG}.mk
calib_scripts:
	${MAKE} iterate_organisms ORG_TASK=calib_scripts_one_org

calib_scripts_one_org:
	@mkdir -p ${MAKE_DIR}
	@echo "Generating calibration script for organism ${ORG}"
	@echo 'include $${RSAT}/makefiles/upstream_calibrations.mk' > ${CALIB_SCRIPT}
	@echo "ORG=${ORG}" >> ${CALIB_SCRIPT}
	@echo "calibrate:" >> ${CALIB_SCRIPT}
	@${MAKE} iterate_sets TASK=one_calib_script
	@echo "Script saved in file ${CALIB_SCRIPT}"

one_calib_script:
	@echo "	make -s -f $${RSAT}/makefiles/motif_disco_competition_2003.mk calibrate_oligos ORG=${ORG} SEQ_LEN=${CURRENT_SEQ_LEN} N=${SEQ_NB} STR=-2str STR=-2str NOOV=-noov" >> ${CALIB_SCRIPT}


################################################################
#### retrieve all upstream sequences of the same length as in the
#### current set
ALL_UP_FILE=${ORG_DIR}/${ORG}_allup${SEQ_LEN}.fasta.gz
ALL_UP_FILE_PURGED_UNCOMP=${ORG_DIR}/${ORG}_allup${SEQ_LEN}_purged.fasta
ALL_UP_FILE_PURGED=${ALL_UP_FILE_PURGED_UNCOMP}.gz
all_up:
	@echo "${ORG}	${CURRENT_SET}	Retrieving all upstream sequences for bakground model"
	@echo ${ALL_UP_FILE}
	@mkdir -p ${UPSTREAML_DIR}
	retrieve-seq -all -org ${ORG} -from -1 -to -${SEQ_LEN} -o ${ALL_UP_FILE}

all_up_purge:
	@echo "Purging all upstream sequences"
	@echo ${ALL_UP_FILE_PURGED}
	purge-sequence -i ${ALL_UP_FILE} -o ${ALL_UP_FILE_PURGED_UNCOMP}
	@echo "Compressing purged upstream sequences"
	gzip -f ${ALL_UP_FILE_PURGED_UNCOMP}

################################################################
#### calculate background oligo frequencies for the current set
BG_STR=-1str
BG_OLIGO_FILE=${UPSTREAML_DIR}/${ORG}_allup${SEQ_LEN}_${OLIGO_LEN}nt${BG_STR}${NOOV}_freq.tab
BG_OLIGO_FILE_PURGED=${UPSTREAML_DIR}/${ORG}_allup${SEQ_LEN}_${OLIGO_LEN}nt${BG_STR}${NOOV}_purged_freq.tab
OLIGO_LEN=6
bg_oligos:
	${MAKE} iterate_oligo_lengths OL_TASK=bg_oligos_one_ol_nopurge

## Calculate oligo frequencies in all upstream sequences
bg_oligos_one_ol: bg_oligos_one_ol_nopurge bg_oligos_one_ol_purge  bg_oligos_one_ol_purged_vs_not
	@echo "${ORG}	${CURRENT_SET}	Calculating background oligo frequencies"
	@echo ${BG_OLIGO_FILE}

## Calculate oligo frequencies in all upstream sequences, non purged
BG_OLIGOS_CMD=oligo-analysis ${NOOV} -i ${ALL_UP_FILE} -l ${OLIGO_LEN} -v ${V} ${BG_STR} -o ${BG_OLIGO_FILE} -return occ,freq
bg_oligos_one_ol_nopurge:
	${MAKE} my_command MY_COMMAND="${BG_OLIGOS_CMD}"

## Calculate oligo frequencies in all upstream sequences, purged
BG_OLIGOS_CMD_PURGE=oligo-analysis ${NOOV} -i ${ALL_UP_FILE_PURGED} -l ${OLIGO_LEN} -v ${V} ${BG_STR} -o ${BG_OLIGO_FILE_PURGED} -return occ,freq
bg_oligos_one_ol_purge:
	${MAKE} my_command MY_COMMAND="${BG_OLIGOS_CMD_PURGE}"

## Compare oligo frequencies between purged and not purged upstream sequences
bg_oligos_one_ol_purged_vs_not:
	compare-scores											\
		-i ${BG_OLIGO_FILE_PURGED}								\
		-i $ ${BG_OLIGO_FILE}  -sc 3								\
		-o ${UPSTREAML_DIR}/${ORG}_allup${SEQ_LEN}_${OLIGO_LEN}nt_${BG_STR}${NOOV}_purged_vs_not.tab
	XYgraph -xcol 2 -ycol 3										\
		-i ${UPSTREAML_DIR}/${ORG}_allup${SEQ_LEN}_${OLIGO_LEN}nt_$ (BG_STR)${NOOV}_purged_vs_not.tab	\
		-o ${UPSTREAML_DIR}/${ORG}_allup${SEQ_LEN}_${OLIGO_LEN}nt_$ (BG_STR)${NOOV}_purged_vs_not.jpg

#pattern_disco: oligos dyads

################################################################
#### Detect over-represented oligonucleotides
#OLIGO_DIR=${ORG_DIR}/oligos
#MODEL=-expfreq ${BG_OLIGO_FILE}
#MODEL_SUFFIX=allup
#OLIGO_FILE=${OLIGO_DIR}/oligos_${CURRENT_SET}_${OLIGO_LEN}nt${STR}${NOOV}_sig${THOSIG}_${MODEL_SUFFIX}
#OLIGO_FILE_HTML=${OLIGO_FILE}.html
#oligos:
#	@echo "${ORG}	${CURRENT_SET}	Detecting over-represented oligonucleotides	${OLIGO_FILE}"
#	@mkdir -p ${OLIGO_DIR}
#	oligo-analysis -v ${V}			\
#		-i ${SEQ_FILE}			\
#		${MODEL}			\
#		-l ${OLIGO_LEN} ${NOOV}		\
#		-return occ,freq,proba,rank	\
#		-sort -thosig ${THOSIG}		\
#		-seqtype dna			\
#		-o ${OLIGO_FILE}
#	text-to-html -i ${OLIGO_FILE} -o ${OLIGO_FILE_HTML}

# ################################################################
# #### synthezise the results
# INDEX_FILE=${ORG}_index.html
# index:
# 	@echo "<html>" > ${INDEX_FILE}
# 	@echo "<table border=1 cellpading=3>" > ${INDEX_FILE}
# 	${MAKE} iterate_sets TASK=index_one_row
# 	@echo "</table>" >> ${INDEX_FILE}
# 	@echo "</html>" >> ${INDEX_FILE}
# 	@echo "${ORG}	Index file generated	${INDEX_FILE}"

# index_one_row:
# 	@echo "<tr>" >> ${INDEX_FILE}
# 	@echo "<td>${CURRENT_SET}</td>" >> ${INDEX_FILE}
# 	@echo "<td><a href=${OLIGO_FILE_HTML}>oligos</a></td>" >> ${INDEX_FILE}
# 	@echo "</tr>" >> ${INDEX_FILE}

################################################################
#### Run multiple-family-analysis for one organism

#### Save the list of fasta files in a text file
FASTA_FILES= ls -1 ${DATA}/${ORG}/*.fasta | grep -v purged
SEQ_LIST_FILE=${ORG}_files.txt
list_fasta_files:
	${FASTA_FILES}
	${FASTA_FILES} > ${SEQ_LIST_FILE}

################################################################
#### Parameters for pattern discovery (oligo-analysis and
#### dyad-analysis)
STR=-2str
THOSIG=0
NOOV=-noov
#NOOV=-ovlp
#MULTI_TASK=purge,oligos,dyads,merge,slide,maps,synthesis,sql
MULTI_TASK=purge,oligos,oligo_maps,synthesis,sql
MULTI_BG=upstream
MULTI_EXP=-bg ${MULTI_BG}
PURGE=-purge
#PURGE=-nopurge
MULTI_DIR=${ORG_DIR}/multi/${MULTI_BG}_bg${PURGE}
MIN_OL=6
MAX_OL=8
MIN_SP=0
MAX_SP=20
SORT=score
MULTI_OPT=
SKIP=
MULTI_CMD=multiple-family-analysis -v ${V}				\
		-thosig ${THOSIG}					\
		${PURGE}						\
		-org ${ORG}						\
		-seq ${SEQ_LIST_FILE}					\
		-outdir ${MULTI_DIR}					\
		${STR}							\
		-minol ${MIN_OL} -maxol ${MAX_OL}			\
		-minsp ${MIN_SP} -maxsp ${MAX_SP}			\
		${MULTI_EXP}						\
		-sort ${SORT} -task ${MULTI_TASK}			\
		${NOOV} ${MULTI_OPT}					\
		${SKIP}							

## generic call for multiple-family-analysis
multi:
	@echo ${MULTI_CMD}
	${MAKE} my_command MY_COMMAND="${MULTI_CMD}"

## run multiple-family-analysis with default upstream calibration
## (same upstream length for all sets)
multi_upstream: 
	${MAKE} multi 

## For the time being, we run dyad-anaysis only with the default
## upstream background
DYAD_TASK=-task dyads,dyad_maps
multi_upstream_dyads: 
	${MAKE} multi MULTI_OPT='${DYAD_TASK}'

## run multiple-family-analysis with upstream frequencies calculated
## for each sequence length
multi_upstreamL:
	make iterate_oligo_lengths OL_TASK=multi_upstreamL_one_length

multi_upstreamL_synthesis:
	${MAKE} multi							\
		MULTI_DIR=${ORG_DIR}/multi/upstreamL_bg${PURGE}		\
		MULTI_EXP='-oligo_exp_freq ${BG_OLIGO_FILE}'		\
		MULTI_TASK=merge_oligos,oligo_maps,report,synthesis

multi_upstreamL_one_length:
	${MAKE} multi											\
		MULTI_DIR=${ORG_DIR}/multi/upstreamL_bg${PURGE}						\
		MULTI_EXP='-oligo_exp_freq ${BG_OLIGO_FILE}' MIN_OL=${OLIGO_LEN} MAX_OL=${OLIGO_LEN}

multi_calibN:
	${MAKE} multi MULTI_BG=calibN MULTI_OPT="-calib_dir ${CALIBN_DIR}" 

CALIBRATE_TASK=-task calibrate
CALIB1_DIR=${ORG_DIR}/calibrations_1gene
multi_calib1:
	${MAKE} multi MULTI_BG=calib1 MULTI_OPT="${CALIBRATE_TASK} -calib_dir ${CALIB1_DIR}"

#multi_calib1_all:
#	@for org in ${ORGANISMS} ; do							\
#		for ol in 5 6 7 8 ; do							\
#			${MAKE} multi_calib1 ORG=$${org} MIN_OL=$${ol} MAX_OL===$${ol};	\
#		done ;									\
#	done

################################################################
## Calculate the effect of the number of sequences on mean and variance
SEQ_NUMBER_SERIES=1 2 3 4 5 6 7 8 9 10 15 20 40 60 80 100

seq_nb_series:
	for nb in ${SEQ_NUMBER_SERIES} ; do		\
		${MAKE} calibrate_oligos N=$${nb} ;	\
	done

HUMAN_LENGTHS=500 1000 1500 2000 3000
seq_nb_series_human:
	@for len in ${HUMAN_LENGTHS}; do						\
	${MAKE} seq_nb_series ORG=Homo_sapiens SEQ_LEN=$${len} REPET=1000;	\
	done

YEAST_LENGTHS=500 1000
seq_nb_series_yeast:
	@for len in ${YEAST_LENGTHS}; do						\
		${MAKE} seq_nb_series ORG=Saccharomyces_cerevisiae SEQ_LEN=$${len} REPET=1000;	\
	done

FLY_LENGTHS=1500 2000 2500 3000
seq_nb_series_fly:
	@for len in ${FLY_LENGTHS}; do							\
		${MAKE} seq_nb_series ORG=Drosophila_melanogaster SEQ_LEN=$${len} REPET=1000 ;	\
	done

MOUSE_LENGTHS=500 1000 1500
seq_nb_series_mouse:
	@for len in ${MOUSE_LENGTHS}; do							\
		${MAKE} seq_nb_series ORG=Mus_musculus SEQ_LEN=$${len} REPET=1000 ;	\
	done

################################################################
## Effet of sequence length on mean and variance
LENGTH_SERIES= 100 200 300 500 1000 1500 2000 2500 3000
seq_length_series:
	@for len in ${LENGTH_SERIES} ; do			\
		${MAKE} calibrate_oligos SEQ_LEN=$${len} N=20 ;	\
	done

LENGTH_SERIES_YEAST=100 200 300 400 500 600 700 800 900 1000
seq_len_series_yeast:
	${MAKE} seq_length_series ORG=Saccharomyces_cerevisiae LENGTH_SERIES='${LENGTH_SERIES_YEAST}'

seq_len_series_human:
	${MAKE} seq_length_series ORG=Homo_sapiens

seq_len_series_mouse:
	${MAKE} seq_length_series ORG=Mus_musculus

seq_len_series_fly:
	${MAKE} seq_length_series ORG=Drosophila_melanogaster

################################################################
## Calculate oligonucleotide distributions of occurrences for random
## gene selections
RAND_DIR=rand_gene_selections
SEQ_LEN=500
CALIB_TASK=all,clean_oligos
START=1
REPET=1000
CALIBN_DIR=${ORG_DIR}/${RAND_DIR}
CURRENT_CALIBN_DIR=${CALIBN_DIR}/${OLIGO_LEN}nt${STR}${NOOV}_N${N}_L${SEQ_LEN}_R${REPET}
CALIBRATE_CMD=								\
	calibrate-oligos -v ${V}					\
		-r ${REPET} -sn ${N} -ol ${OLIGO_LEN} -sl ${SEQ_LEN}	\
		-task ${CALIB_TASK}					\
		-start ${START}						\
		${END}							\
		${STR} ${NOOV}						\
		-outdir ${CURRENT_CALIBN_DIR}					\
		-org ${ORG}

## Run the program immediately (WHEN=now) or submit it to a queue (WHEN=queue)
WHEN=now
N=5
calibrate_oligos:
	${MAKE} calibrate_oligos_${WHEN}

## Run the program immediately
calibrate_oligos_now:
	(umask 0002; ${CALIBRATE_CMD})

## Submit the program to a queue
JOB_DIR=jobs
JOB=`mktemp ${JOB_DIR}/job.XXXXXX`

calibrate_oligos_queue:
	@mkdir -p ${JOB_DIR}
	@for job in ${JOB} ; do						\
		echo "Job $${job}" ;					\
		echo "${CALIBRATE_CMD}" > $${job} ;		\
		qsub -m e -q rsa@merlin.ulb.ac.be -N $${job} -j oe	\
			-o $${job}.log $${job} ;	\
	done

## A quick test for calibrate_oligos
test: calibN_test all_up_test bg_oligos_test

calib1_test:
	${MAKE} calibrate_oligos ORG=Mycoplasma_genitalium N=1 SEQ_LEN=100 REPET=100 OLIGO_LEN=4

calibN_test:
	${MAKE} calibrate_oligos ORG=Mycoplasma_genitalium N=1 SEQ_LEN=100 REPET=100 OLIGO_LEN=4

all_up_test:
	${MAKE} all_up ORG=Mycoplasma_genitalium SEQ_LEN=100 MIN_OL=6 MAX_OL=4

bg_oligos_test:
	${MAKE} bg_oligos_one_ol_nopurge ORG=Mycoplasma_genitalium SEQ_LEN=100 OLIGO_LEN=4 V=3


ACCESS=results
give_access:
	find ${ACCESS} -type d -exec chmod 775 {} \;
	find ${ACCESS} -type f -exec chmod 664 {} \;


## ##############################################################
## Fit all previously calculated distributions with a Poisson and a
## negbin, respectively
## THIS IS OBSOLETE: calibrate-oligos now automatically exports the negbin and poisson files

DISTRIB_FILES=`find ${RES_DIR} -name '*_distrib.tab'`
GOOD_DISTRIB_FILES=`find ${RES_DIR} -name '*_distrib.tab' -exec wc {} \; | awk '$$1 >= 2000 {print $$4}'`
DISTRIB_LAW=negbin
FITTING_CMD=					\
	fit-distribution -v ${V}		\
	-distrib ${DISTRIB_LAW}			\
	-i ${DISTRIB_FILE}			\
	-o ${FITTING_FILE}

find_distrib_files:
	@echo ${DISTRIB_FILES}

good_distrib_files:
	@echo
	@echo "Good distrib files"
	@find ${RES_DIR} -name '*_distrib.tab' -exec wc {} \; | awk '$$1 >= 2000'

bad_distrib_files:
	@echo
	@echo "Bad distrib files"
	@find ${RES_DIR} -name '*_distrib.tab' -exec wc {} \; | awk '$$1 < 2000'

good_fittings:
	@echo  ${GOOD_DISTRIB_FILES}
	@for file in ${GOOD_DISTRIB_FILES} ; do			\
		${MAKE} one_fitting DISTRIB_FILE=$${file} ;	\
	done

check_distrib_files: good_distrib_files bad_distrib_files

all_fittings:
	${MAKE} all_fittings_${WHEN} DISTRIB_LAW=poisson
	${MAKE} all_fittings_${WHEN} DISTRIB_LAW=negbin


all_fittings_queue:
	@for infile in ${DISTRIB_FILES} ; do							\
		${MAKE} one_fitting_queue DISTRIB_FILE=${WORK_DIR}/$${infile} \
		FITTING_FILE=${WORK_DIR}/`echo $${infile} | perl -pe 's/distrib.tab/${DISTRIB_LAW}.tab/'` \
		JOB=`mktemp ${JOB_DIR}/job.XXXXXX`;		\
	done

all_fittings_now:
	@for infile in ${DISTRIB_FILES} ; do							\
		fit-distribution -v 1 -distrib ${DISTRIB_LAW} -i $${infile} -o ${WORK_DIR}/`echo $${infile} | perl -pe 's/distrib.tab/${DISTRIB_LAW}.tab/'`;				\
	done


## Submit the program to a queue
one_fitting_queue:
	@mkdir -p ${JOB_DIR}
	@for job in ${JOB} ; do						\
	echo "Job $${job}";						\
	echo "${FITTING_CMD}" > $${job};				\
	qsub -m e -q rsa@merlin.ulb.ac.be -N $${job} -j oe		\
		-o $${job}.log $${job};					\
	done

one_fitting_one_law:
	@echo "Fitting ${DISTRIB_LAW}	${FITTING_FILE}"
	@fit-distribution -v 1 -distrib ${DISTRIB_LAW} -i ${DISTRIB_FILE} -o ${FITTING_FILE} 
#	@head -50 ${DISTRIB_FILE}| fit-distribution -v 1 -distrib ${DISTRIB_LAW}  -o ${FITTING_FILE} 

_distrib.tab_poisson.tab:
	@fit-distribution -v 1 -distrib poisson -i $< -o $@

_distrib.tab_negbin.tab:
	@fit-distribution -v 1 -distrib negbin -i $< -o $@

###################################################################
# extract mean and variance from random distributions (stat files)
###################################################################

#LIST_STATS_FILES=`find ${ORG_DIR}/${RAND_DIR}/${OLIGO_LEN}nt${STR}${NOOV}_*_L${SEQ_LEN}_R${REPET} -name '*_stats.tab'`

list_stats_files:
	rm -f ${ORG}_${OLIGO_LEN}nt_${STR}${NOOV}_l${SEQ_LEN}_r${REPET}_files.tmp
	@for n in ${SEQ_NUMBER_SERIES};do				\
	find ${ORG_DIR}/${RAND_DIR}/${OLIGO_LEN}nt${STR}${NOOV}_N$${n}_L${SEQ_LEN}_R${REPET}/${ORG}_${OLIGO_LEN}nt_${STR}${NOOV}_n$${n}_l${SEQ_LEN}_r${REPET}_stats.tab >> ${ORG_DIR}/${RAND_DIR}/${ORG}_${OLIGO_LEN}nt_${STR}${NOOV}_l${SEQ_LEN}_r${REPET}_files.tmp;	\
	done

LIST_STATS_FILES=`more ${ORG_DIR}/${RAND_DIR}/${ORG}_${OLIGO_LEN}nt_${STR}${NOOV}_l${SEQ_LEN}_r${REPET}_files.tmp`

join_mean_var_N:
	${MAKE} list_stats_files
	compare-scores -sc 3 -files ${LIST_STATS_FILES} > \
	${ORG_DIR}/${RAND_DIR}/${OLIGO_LEN}nt${STR}${NOOV}_L${SEQ_LEN}_R${REPET}_allN_means.tab
	compare-scores -sc 4 -files ${LIST_STATS_FILES} > \
	${ORG_DIR}/${RAND_DIR}/${OLIGO_LEN}nt${STR}${NOOV}_L${SEQ_LEN}_R${REPET}_allN_vars.tab
	rm -f ${ORG_DIR}/${RAND_DIR}/${ORG}_${OLIGO_LEN}nt_${STR}${NOOV}_l${SEQ_LEN}_r${REPET}_files.tmp

join_all_mean_var_N:
	${MAKE} join_mean_var_N ORG=Saccharomyces_cerevisiae SEQ_LEN=500 REPET=1000
	${MAKE} join_mean_var_N ORG=Saccharomyces_cerevisiae SEQ_LEN=1000 REPET=1000
	${MAKE} join_mean_var_N ORG=Homo_sapiens SEQ_LEN=1000 REPET=1000
	${MAKE} join_mean_var_N ORG=Homo_sapiens SEQ_LEN=2000 REPET=1000
	${MAKE} join_mean_var_N ORG=Drosophila_melanogaster SEQ_LEN=1500 REPET=1000
	${MAKE} join_mean_var_N ORG=Drosophila_melanogaster SEQ_LEN=2000 REPET=1000
	${MAKE} join_mean_var_N ORG=Drosophila_melanogaster SEQ_LEN=2500 REPET=1000
	${MAKE} join_mean_var_N ORG=Drosophila_melanogaster SEQ_LEN=3000 REPET=1000
	${MAKE} join_mean_var_N ORG=Mus_musculus SEQ_LEN=500 REPET=1000
	${MAKE} join_mean_var_N ORG=Mus_musculus SEQ_LEN=1000 REPET=1000
	${MAKE} join_mean_var_N ORG=Mus_musculus SEQ_LEN=1500 REPET=1000
#	${MAKE} join_mean_var_N ORG=Saccharomyces_cerevisiae SEQ_LEN=500 REPET=10000
#	${MAKE} join_mean_var_N ORG=Saccharomyces_cerevisiae SEQ_LEN=1000 REPET=10000
#	${MAKE} join_mean_var_N ORG=Homo_sapiens SEQ_LEN=1000 REPET=10000
#	${MAKE} join_mean_var_N ORG=Homo_sapiens SEQ_LEN=2000 REPET=10000



BACKGROUND=calibN
OLD_BG_DIR=results/multi/${BACKGROUND}_bg
OLD_DIR=${OLD_BG_DIR}/${ORG}
NEW_DIR=results/${ORG}/multi/${BACKGROUND}_bg/

reorganize: rm_junk organize_all_bg mv_calib1 discard_oldies mv_seq_lengths mv_background_dirs

rm_junk:
	find . -name .DS_Store -exec rm {} \;

organize_all_bg:
	rm -rf results/multi/upstream_bg/Saccharomyces_cerevisiae_bk
	for bg in calibN calib1 upstream; do		\
		${MAKE} organize_one_bg BACKGROUND=$${bg};	\
	done

organize_one_bg:
	${MAKE} iterate_organisms ORG_TASK=organize_one_dir
	rmdir ${OLD_BG_DIR}

organize_one_dir:
	@echo "Old dir	${OLD_DIR}"
	@echo "New dir	${NEW_DIR}"
	@mkdir -p ${NEW_DIR}
	rsync -ruptvl ${OLD_DIR}/* ${NEW_DIR}/
	rm -rf ${OLD_DIR}

mv_calib1:
	${MAKE} iterate_organisms ORG_TASK=mv_calib1_one_org

mv_calib1_one_org:
	@mkdir -p results/${ORG}/calibrations_1gene/
	rsync -ruptvl results/${ORG}/multi/calib1_bg/calibrations/* results/${ORG}/calibrations_1gene/
	rm -rf results/${ORG}/multi/calib1_bg/calibrations

discard_oldies:
	mkdir -p old_results
	rsync -ruptvl results/multi old_results/
	rm -rf results/multi

mv_seq_lengths: mv_seq_lengths_mouse mv_seq_lengths_human mv_seq_lengths_fly mv_seq_lengths_yeast

mv_seq_lengths_mouse:
	mkdir -p results/Mus_musculus/sequence_lengths
	rsync -ruptvl results/Mus_musculus/mus*/*_lengths.txt results/Mus_musculus/sequence_lengths/
	rm -rf results/Mus_musculus/mus*

mv_seq_lengths_yeast:
	mkdir -p results/Saccharomyces_cerevisiae/sequence_lengths
	rsync -ruptv results/Saccharomyces_cerevisiae/yst*/*_lengths.txt results/Saccharomyces_cerevisiae/sequence_lengths/
	rm -rf results/Saccharomyces_cerevisiae/yst*

mv_seq_lengths_human:
	mkdir -p results/Homo_sapiens/sequence_lengths
	rsync -ruptv results/Homo_sapiens/hm*/*_lengths.txt results/Homo_sapiens/sequence_lengths/
	rm -rf results/Homo_sapiens/hm*

mv_seq_lengths_fly:
	mkdir -p results/Drosophila_melanogaster/sequence_lengths
	rsync -ruptv results/Drosophila_melanogaster/dm*/*_lengths.txt results/Drosophila_melanogaster/sequence_lengths/
	rm -rf results/Drosophila_melanogaster/dm*

mv_background_dirs:
	${MAKE} iterate_organisms ORG_TASK=mv_background_dir_one_org

mv_background_dir_one_org:
	mv results/${ORG}/background_frequencies results/${ORG}/upstreamL_frequencies 


################################################################
## Check the difference between different calibrations
## WARNING: there is an obvious problem with the files calib1: the
## frequency and variance is twice lower than with other calibration
## systems
COMPARISON_DIR=results/${ORG}/calib_comparisons
COMPARISON_PREFIX=${OLIGO_LEN}nt_${ORG}_L${SEQ_LEN}${STR}${NOOV}_R${REPET}_comparison
COMPARISON_FILE = ${COMPARISON_DIR}/${COMPARISON_PREFIX}
CALIB1_PREFIX=${OLIGO_LEN}nt_upstream_L${SEQ_LEN}_${ORG}${NOOV}${STR}
CALIB1_FILE=${CALIB1_DIR}/${CALIB1_PREFIX}_negbin.tab
CALIBN_PREFIX=${ORG}_${OLIGO_LEN}nt_${STR}${NOOV}_n1_l${SEQ_LEN}_r${REPET}
CALIBN_FILE= results/${ORG}/rand_gene_selections/${OLIGO_LEN}nt${STR}${NOOV}_N1_L${SEQ_LEN}_R${REPET}/${CALIBN_PREFIX}_negbin.tab
compare_calibrations:
	mkdir -p ${COMPARISON_DIR}
	compare-scores -sc1 5 -sc2 4			\
		-i ${CALIB1_FILE}			\
		-i ${CALIBN_FILE}			\
		-o ${COMPARISON_FILE}_variance.tab

	XYgraph -i ${COMPARISON_FILE}_variance.tab		\
		-size 600					\
		-xcol 2 -ycol 3					\
		-title1 "${COMPARISON_PREFIX}"			\
		-xleg1 "Single-gene calibration (variance)"	\
		-xleg2 "${CALIB1_PREFIX}"			\
		-yleg1 "Random selection, 1 gene (variance)"	\
		-yleg2 "${CALIBN_PREFIX}"			\
		-o ${COMPARISON_FILE}_variance.jpg

	compare-scores -sc1 3 -sc2 2 -sc3 4	\
		-i ${CALIB1_FILE}		\
		-i ${CALIBN_FILE}		\
		-i ${BG_OLIGO_FILE}		\
		-o ${COMPARISON_FILE}_mean.tab

	XYgraph -i ${COMPARISON_FILE}_mean.tab			\
		-size 600					\
		-xcol 2 -ycol 3					\
		-title1 "${COMPARISON_PREFIX}"			\
		-xleg1 "Single-gene calibration (mean)"		\
		-xleg2 "${CALIB1_PREFIX}"			\
		-yleg1 "Random selection, 1 gene (mean)"	\
		-yleg2 "${CALIBN_PREFIX}"			\
		-o ${COMPARISON_FILE}_mean.jpg


some_comparisons:
	make compare_calibrations ORG=Saccharomyces_cerevisiae SEQ_LEN=500 OLIGO_LEN=6
	make compare_calibrations ORG=Drosophila_melanogaster SEQ_LEN=2000 OLIGO_LEN=7

################################################################
## Compile report
report: generate_reports compile_reports

generate_reports:
	${MAKE} iterate_organisms		\
		ORG_TASK=multi_calibN		\
		MULTI_TASK=report,synthesis

compile_reports: 
	@echo "Generating report	${REPORT_DIR}"
	@${MAKE} report_headers
	@${MAKE} report_organisms
	@echo "Report generated"
	@echo ${RESULT_FILE}
	@echo ${PARAM_FILE}

REPORT_DIR=${RES_DIR}/reports
RESULT_FILE=${REPORT_DIR}/Simonis_vanHelden_results.txt
PARAM_FILE=${REPORT_DIR}/Simonis_vanHelden_parameters.txt
report_headers:
	@echo ">name of contact" > ${RESULT_FILE}
	@echo "Jacques van Helden" >> ${RESULT_FILE}
	@echo ">email" >> ${RESULT_FILE}
	@echo "jvanheld@bigre.ulb.ac.be" >> ${RESULT_FILE}
	@echo ">program name" >> ${RESULT_FILE}
	@echo "multiple-family-analyis" >> ${RESULT_FILE}
	@cp ${RESULT_FILE} ${PARAM_FILE}

report_organisms:
	${MAKE} iterate_organisms ORG_TASK=report_one_organism

ORG_REPORT_DIR=${ORG_DIR}/multi/calibN_bg-purge/mdc_report
ORG_RESULT_FILE=${ORG_REPORT_DIR}/${ORG}_files_bg_calibN-purge_6nt_8nt-noov-2str_sig0_results.txt
ORG_PARAM_FILE=${ORG_REPORT_DIR}/${ORG}_files_bg_calibN-purge_6nt_8nt-noov-2str_sig0_parameters.txt
report_one_organism:
	cat ${ORG_RESULT_FILE} >> ${RESULT_FILE}
	cat ${ORG_PARAM_FILE} >> ${PARAM_FILE}

################################################################
## Archive the reports
DATE=`date +%Y%m%d_%H%M%S`
ARCHIVE_DIR=archives/
archive:
	mkdir -p ${ARCHIVE_DIR}
	${MAKE} archive_report REPORT="${ARCHIVE_DIR}/report_${DATE}.zip"

archive_report:
	@echo ${REPORT}
	find results/*/multi -name '*_results.txt' -exec zip -ry ${REPORT} {} \;
	find results/*/multi -name '*_parameters.txt' -exec zip -ry ${REPORT} {} \;
	zip -ry -n '*~' ${REPORT} results/reports
	find results/*/multi -name '*_selection' -exec zip -ry ${REPORT} {} \;
	zip -ry ${REPORT} file_index.html
	@echo report archived in ${REPORT}

# oligo-analysis  -v 3 -l 7 -noov -2str -i calibrations/tmp_all_up_2000.fasta -return occ -distrib


################################################################
## Last minut task list
to_clean:
	rm -rf results/Mus_musculus/rand_gene_selections/6nt-2str-noov_N0_L_R10000
	rm -rf results/*/multi/calib1_bg
	rm -rf results/*/multi/calibN_bg
	rm -rf results/*/multi/upstream_bg
done:
	make iterate_organisms ORG_TASK=multi_calibN MIN_OL=6 MAX_OL=6
	make iterate_organisms ORG_TASK=multi_calibN MIN_OL=7 MAX_OL=7

running_on_liv:
	make -f makefiles/calibrate_Drosophila_melanogaster_sorted.mk calibrate OLIGO_LEN=5 >& calib_Drosophila_melanogaster_5nt_log.txt &

running_on_brol:
	make iterate_organisms ORG_TASK=multi_upstream_dyads MIN_OL=5 MAX_OL=8

queued_on_merlin:
	make -f makefiles/calibrate_Drosophila_melanogaster_sorted.mk calibrate OLIGO_LEN=5 >& calib_Drosophila_melanogaster_5nt_log.txt WHEN=queue &
	make -sk iterate_organisms ORG_TASK=multi_calibN MIN_OL=5 MAX_OL=5 WHEN=queue
	make -sk iterate_organisms ORG_TASK=multi_calibN MIN_OL=8 MAX_OL=8 WHEN=queue
	make iterate_organisms ORG_TASK=multi_upstream_dyads MIN_OL=5 MAX_OL=8 WHEN=queue
	make iterate_organisms ORG_TASK=multi_calib1 MIN_OL=6 MAX_OL=6 WHEN=queue
	make iterate_organisms ORG_TASK=multi_calib1 MIN_OL=5 MAX_OL=5 WHEN=queue
	make iterate_organisms ORG_TASK=multi_calib1 MIN_OL=7 MAX_OL=7 WHEN=queue
	make iterate_organisms ORG_TASK=multi_calib1 MIN_OL=8 MAX_OL=8 WHEN=queue

	make iterate_organisms ORG_TASK=multi_calibN MIN_OL=6 MAX_OL=6 WHEN=queue
	make iterate_organisms ORG_TASK=multi_calibN MIN_OL=7 MAX_OL=7 WHEN=queue
	make iterate_organisms ORG_TASK=multi_calibN MIN_OL=5 MAX_OL=5 WHEN=queue
	make iterate_organisms ORG_TASK=multi_calibN MIN_OL=8 MAX_OL=8 WHEN=queue

	make -s -f ${RSAT}/makefiles/motif_disco_competition_2003.mk calibrate_oligos ORG=Homo_sapiens SEQ_LEN=2000 N=35 STR=-2str STR=-2str NOOV=-noov REPET=100 WHEN=queue OLIGO_LEN=8

to_do:
	make iterate_organisms ORG_TASK=multi_calibN THOSIG=-1
	make iterate_organisms ORG_TASK=multi_calibN PURGE=-nopurge MIN_OL=8 MAX_OL=8 WHEN=queue
	make iterate_organisms ORG_TASK=multi_calibN PURGE=-nopurge MIN_OL=5 MAX_OL=8 WHEN=queue THOSIG=-1


## ##############################################################
## publish some directory on the web site

PUBLISH_SITE=jvanheld@www.bigre.ulb.ac.be:public_html/motif_discovery_competition_2003/
PUBLISH_DIR=evaluation
publish_eval:
	${RSYNC} --exclude '*~' ${PUBLISH_DIR} ${PUBLISH_SITE}

TO_PUBLISH=file_index.html
#PUBLISH_SERVER=jvanheld@rsat.ulb.ac.be:
PUBLISH_SITE=${PUBLISH_SERVER}${PUBLISH_DIR}
PUBLISH_DIR=/home/jvanheld/rsa-tools/public_html/data/motif_discovery_competition
publish:
	for pub in ${TO_PUBLISH}; do					\
		${MAKE} publish_one_item ITEM_TO_PUBLISH=$${pub} ;	\
	done

ITEM_TO_PUBLISH=results/Homo_sapiens/multi/calibN_bg-purge
ITEM_DIR=`dirname ${ITEM_TO_PUBLISH}`
publish_one_item:
	mkdir -p ${PUBLISH_DIR}/${ITEM_DIR}
	${RSYNC} ${ITEM_TO_PUBLISH} ${PUBLISH_SITE}/${ITEM_DIR}/
