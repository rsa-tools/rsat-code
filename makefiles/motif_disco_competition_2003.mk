.Suffixes: _distrib.tab _poisson.tab _negbin.tab

################################################################
## makefile for getting info for the motif discovery competition


################################################################
#### General parameters

## level of verbosity
V=1

################################################################
### commands
MAKEFILE=${RSAT}/makefiles/motif_disco_competition_2003.mk
MAKE=make -s -f ${MAKEFILE}

SSH_OPT = -e ssh 
RSYNC_OPT= -ruptvlz  ${SSH_OPT} 
RSYNC = rsync  ${RSYNC_OPT}

WGET=wget --passive-ftp -np -rNL

################################################################
#### list of targets
usage:
	@echo "usage: make [-OPT='options'] target"
	@echo "implemented targets"
	@perl -ne 'if (/^([a-z]\S+):/){ print "\t$$1\n";  }' ${MAKEFILE}

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
to_server:
	${MAKE} rsync_one_dir TO_SYNC=TASK_LIST.html
	${MAKE} rsync_one_dir TO_SYNC=results

one_dir_to_server:
	${RSYNC} ${EXCLUDE} ${TO_SYNC} ${SERVER_LOCATION}

update_from_server:
	${RSYNC} ${EXCLUDE} ${SERVER_LOCATION}${TO_SYNC} .

## Synchronize calibrations from merlin
from_merlin:
	${MAKE} update_from_server SERVER_DIR=motif_discovery_competition_2003/ TO_SYNC=results

# Temporary
#	${MAKE} update_from_server SERVER_DIR=./ TO_SYNC=results
#	${MAKE} update_from_server SERVER_DIR=makefiles/ TO_SYNC=results
#	${MAKE} update_from_server TO_SYNC='*.xls'
#	${MAKE} update_from_server TO_SYNC='*.tab'
#	${MAKE} update_from_server TO_SYNC=my_calibrate-oligos.R


from_liv:
	${MAKE} update_from_server SERVER=liv.bmc.uu.se ${SERVER_DIR}=motif_discovery_competition_2003/

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

################################################################
## Iterate over all organisms
ORGANISMS=Saccharomyces_cerevisiae Drosophila_melanogaster Mus_musculus Homo_sapiens
ORG_TASK=calib_script
iterate_organisms:
	@echo "Iterating task ${ORG_TASK} over organisms"
	@echo ${ORGANISMS}
	@for org in ${ORGANISMS} ; do			\
		${MAKE} ${ORG_TASK} ORG=$${org} ;	\
	done

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
background: seq_len all_up bg_oligos

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
	@echo "include ${RSAT}/makefiles/upstream_calibrations.mk" > ${CALIB_SCRIPT}
	@echo "ORG=${ORG}" >> ${CALIB_SCRIPT}
	@echo "calibrate:" >> ${CALIB_SCRIPT}
	@${MAKE} iterate_sets TASK=one_calib_script
	@echo "Script saved in file ${CALIB_SCRIPT}"

one_calib_script:
	@echo "	${MAKE} calibrate_oligos ORG=${ORG} SEQ_LEN=${CURRENT_SEQ_LEN} N=${SEQ_NB} STR=-2str STR=-2str NOOV=-noov" >> ${CALIB_SCRIPT}


################################################################
#### retrieve all upstream sequences of the same length as in the
#### current set
ALL_UP_FILE=${UPSTREAML_DIR}/${ORG}_allup${SEQ_LEN}.fasta.gz
ALL_UP_FILE_PURGED_UNCOMP=${UPSTREAML_DIR}/${ORG}_allup${SEQ_LEN}_purged.fasta
ALL_UP_FILE_PURGED=${ALL_UP_FILE_PURGED_UNCOMP}.gz
all_up:
	@echo "${ORG}	${CURRENT_SET}	Retrieving all upstream sequences for bakground model"
	@echo ${ALL_UP_FILE}
	@mkdir -p ${UPSTREAML_DIR}
	retrieve-seq -all -org ${ORG} -from -1 -to -${SEQ_LEN} -o ${ALL_UP_FILE}
	@echo "Purging all upstream sequences"
	@echo ${ALL_UP_FILE_PURGED}
	purge-sequence -i ${ALL_UP_FILE} -o ${ALL_UP_FILE_PURGED_UNCOMP}
	@echo "Compressing purged upstream sequences"
	gzip -f ${ALL_UP_FILE_PURGED_UNCOMP}

################################################################
#### calculate background oligo frequencies for the current set
BG_OLIGO_FILE=${UPSTREAML_DIR}/${ORG}_allup${SEQ_LEN}_${OLIGO_LEN}nt_1str${NOOV}_freq.tab
BG_OLIGO_FILE_PURGED=${UPSTREAML_DIR}/${ORG}_allup${SEQ_LEN}_${OLIGO_LEN}nt_1str${NOOV}_purged_freq.tab
OLIGO_LEN=6
bg_oligos:
	${MAKE} iterate_oligo_lengths OL_TASK=bg_oligos

bg_oligos_one_ol:
	@echo "${ORG}	${CURRENT_SET}	Calculating background oligo frequencies"
	@echo ${BG_OLIGO_FILE}
	oligo-analysis ${NOOV} -i ${ALL_UP_FILE} -l ${OLIGO_LEN} -v ${V} -1str -o ${BG_OLIGO_FILE} -return occ,freq
	oligo-analysis ${NOOV} -i ${ALL_UP_FILE_PURGED} -l ${OLIGO_LEN} -v ${V} -1str -o ${BG_OLIGO_FILE_PURGED} -return occ,freq
	compare-scores											\
		-i ${BG_OLIGO_FILE_PURGED}								\
		-i $ ${BG_OLIGO_FILE}  -sc 3								\
		-o ${UPSTREAML_DIR}/${ORG}_allup${SEQ_LEN}_${OLIGO_LEN}nt_1str${NOOV}_purged_vs_not.tab
	XYgraph -xcol 2 -ycol 3										\
		-i ${UPSTREAML_DIR}/${ORG}_allup${SEQ_LEN}_${OLIGO_LEN}nt_1str${NOOV}_purged_vs_not.tab	\
		-o ${UPSTREAML_DIR}/${ORG}_allup${SEQ_LEN}_${OLIGO_LEN}nt_1str${NOOV}_purged_vs_not.jpg

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

################################################################
#### synthezise the results
INDEX_FILE=${ORG}_index.html
index:
	@echo "<html>" > ${INDEX_FILE}
	@echo "<table border=1 cellpading=3>" > ${INDEX_FILE}
	${MAKE} iterate_sets TASK=index_one_row
	@echo "</table>" >> ${INDEX_FILE}
	@echo "</html>" >> ${INDEX_FILE}
	@echo "${ORG}	Index file generated	${INDEX_FILE}"

index_one_row:
	@echo "<tr>" >> ${INDEX_FILE}
	@echo "<td>${CURRENT_SET}</td>" >> ${INDEX_FILE}
	@echo "<td><a href=${OLIGO_FILE_HTML}>oligos</a></td>" >> ${INDEX_FILE}
	@echo "</tr>" >> ${INDEX_FILE}

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
MULTI_TASK=purge,oligos,maps,synthesis,sql
MULTI_BG=upstream
MULTI_EXP=-bg ${MULTI_BG}
PURGE=-purge
#PURGE=-nopurge
MULTI_DIR=${ORG_DIR}/multi/${MULTI_BG}_bg${PURGE}
MIN_OL=5
MAX_OL=8
MIN_SP=0
MAX_SP=20
SORT=score
MULTI_OPT=
MULTI_CMD=multiple-family-analysis -v ${V}				\
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
		-user jvanheld -password jvanheld -schema multifam


## generic call for multiple-family-analysis
multi:
	@echo ${MULTI_CMD}
	${MAKE} my_command MY_COMMAND="${MULTI_CMD}"

## run multiple-family-analysis with default upstream calibration
## (same lengths for all sets)
DYAD_TASK=-task dyads
multi_upstream: 
	${MAKE} multi MULTI_OPT='${DYAD_TASK}'

## run multiple-family-analysis with upstream frequencies calculated
## for each sequence length
#multi_upstreamL:
#	${MAKE} multi MULTI_DIR=${ORG_DIR}/multi/upstreamL_bg MULTI_EXP='-oligo_exp_freq ${BG_OLIGO_FILE}' MIN_OL=6 MAX_OL=6

multi_calibN:
	${MAKE} multi MULTI_BG=calibN MULTI_OPT="-calib_dir ${CALIBN_DIR}" MIN_OL=6 MAX_OL=6

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
calibrate_oligos_test:
	${MAKE} calibrate_oligos ORG=Mycoplasma_genitalium N=10 SEQ_LEN=200 STR=-1str NOOV=-ovlp R=10


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
	# ${MAKE} join_mean_var_N ORG=Saccharomyces_cerevisiae SEQ_LEN=500 REPET=10000
	# ${MAKE} join_mean_var_N ORG=Saccharomyces_cerevisiae SEQ_LEN=1000 REPET=10000
	# ${MAKE} join_mean_var_N ORG=Homo_sapiens SEQ_LEN=1000 REPET=10000
	# ${MAKE} join_mean_var_N ORG=Homo_sapiens SEQ_LEN=2000 REPET=10000

my_command:
	${MAKE} command_${WHEN}

command_queue:
	@mkdir -p ${JOB_DIR}
	@for job in ${JOB} ; do						\
		echo "Job $${job}" ;					\
		echo "${MY_COMMAND}" > $${job} ;		\
		qsub -m e -q rsa@merlin.ulb.ac.be -N $${job} -j oe	\
			-o $${job}.log $${job} ;	\
	done

command_now:
	${MY_COMMAND}


## ##############################################################
## publish some directory on the web site

PUBLISH_SITE=jvanheld@www.scmbb.ulb.ac.be:public_html/motif_discovery_competition_2003/
PUBLISH_DIR=evaluation
publish:
	${RSYNC} --exclude '*~' ${PUBLISH_DIR} ${PUBLISH_SITE}


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
