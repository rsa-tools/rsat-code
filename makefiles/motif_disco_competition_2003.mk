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
RSYNC_OPT= -ruptvl  ${SSH_OPT} 
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

#### publish
SYNC_DIR=results
#SERVER_LOCATION=rubens.ulb.ac.be:/rubens/dsk3/genomics/motif_discovery_competition/
SERVER_LOCATION=merlin.ulb.ac.be:motif_discovery_competition_2003/
EXCLUDE=--exclude '*~' --exclude oligos --exclude '*.wc' --exclude random_genes.tab
publish:
	${RSYNC} ${EXCLUDE} ${SYNC_DIR} ${SERVER_LOCATION}

update_from_server:
	${RSYNC} ${EXCLUDE} ${SERVER_LOCATION}${SYNC_DIR} .

## Synchronize calibrations from merlin
from_merlin:
	${MAKE} update_from_server SERVER_LOCATION=jvanheld@merlin.ulb.ac.be:

from_liv:
	${MAKE} update_from_server SERVER_LOCATION=jvanheld@liv.bmc.uu.se:/Users/jvanheld/motif_discovery_competition_2003/


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

RES_DIR=results
CURRENT_RES_DIR=${RES_DIR}/${ORG}/${CURRENT_SET}
################################################################
#### Create the result directory for the current set
current_dir:
	@mkdir -p ${CURRENT_RES_DIR}

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
#### Calculate background model, depending on the size of the sequences in the current set
EXP_FREQ_DIR=${RES_DIR}/${ORG}/background_frequencies
background: seq_len all_up bg_oligos 

################################################################
#### report sequence length for the current set
CURRENT_SEQ_LEN=`sequence-lengths -i ${SEQ_FILE} | cut -f 2 | sort -u`
SEQ_NB=`sequence-lengths -i ${SEQ_FILE} | wc -l | awk '{print $$1}'`
SEQ_LEN_FILE=${CURRENT_RES_DIR}/${CURRENT_SET}_lengths.txt
seq_len:  current_dir
	@echo "set ${CURRENT_SET} length ${CURRENT_SEQ_LEN}	${SEQ_NB} seqs	${SEQ_LEN_FILE}"
	@echo ${CURRENT_SEQ_LEN} > ${SEQ_LEN_FILE}

MAKE_DIR=makefiles
CALIB_SCRIPT=${MAKE_DIR}/calibrate_${ORG}.mk
calib_script:
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
ALL_UP_FILE=${EXP_FREQ_DIR}/${ORG}_allup${SEQ_LEN}.fasta.gz
ALL_UP_FILE_PURGED_UNCOMP=${EXP_FREQ_DIR}/${ORG}_allup${SEQ_LEN}_purged.fasta
ALL_UP_FILE_PURGED=${ALL_UP_FILE_PURGED_UNCOMP}.gz
all_up:
	@echo "${ORG}	${CURRENT_SET}	Retrieving all upstream sequences for bakground model"
	@mkdir -p ${EXP_FREQ_DIR}
	retrieve-seq -all -org ${ORG} -from -1 -to -${SEQ_LEN} -o ${ALL_UP_FILE}
	purge-sequence -i ${ALL_UP_FILE} -o ${ALL_UP_FILE_PURGED_UNCOMP}
	gzip -f ${ALL_UP_FILE_PURGED_UNCOMP}

################################################################
#### calculate background oligo frequencies for the current set
BG_OLIGO_FILE=${EXP_FREQ_DIR}/${ORG}_allup${SEQ_LEN}_${OLIGO_LEN}nt_1str${NOOV}_freq.tab
BG_OLIGO_FILE_PURGED=${EXP_FREQ_DIR}/${ORG}_allup${SEQ_LEN}_${OLIGO_LEN}nt_1str${NOOV}_purged_freq.tab
OLIGO_LEN=6
bg_oligos:
	@echo "${ORG}	${CURRENT_SET}	Calculating background oligo frequencies"
	oligo-analysis ${NOOV} -i ${ALL_UP_FILE} -l ${OLIGO_LEN} -v ${V} -1str -o ${BG_OLIGO_FILE} -return occ,freq
	oligo-analysis ${NOOV} -i ${ALL_UP_FILE_PURGED} -l ${OLIGO_LEN} -v ${V} -1str -o ${BG_OLIGO_FILE_PURGED} -return occ,freq
	compare-scores											\
		-i ${BG_OLIGO_FILE_PURGED}								\
		-i $ ${BG_OLIGO_FILE}  -sc 3								\
		-o ${EXP_FREQ_DIR}/${ORG}_allup${SEQ_LEN}_${OLIGO_LEN}nt_1str${NOOV}_purged_vs_not.tab
	XYgraph -xcol 2 -ycol 3										\
		-i ${EXP_FREQ_DIR}/${ORG}_allup${SEQ_LEN}_${OLIGO_LEN}nt_1str${NOOV}_purged_vs_not.tab	\
		-o ${EXP_FREQ_DIR}/${ORG}_allup${SEQ_LEN}_${OLIGO_LEN}nt_1str${NOOV}_purged_vs_not.jpg

################################################################
#### Parameters for pattern discovery (oligo-analysis and
#### dyad-analysis)
STR=-2str
THOSIG=0
NOOV=-noov
pattern_disco: oligos dyads

################################################################
#### Detect over-represented oligonucleotides
OLIGO_DIR=${CURRENT_RES_DIR}/oligos
OLIGO_FILE=${OLIGO_DIR}/oligos_${CURRENT_SET}_${OLIGO_LEN}nt${STR}${NOOV}_sig${THOSIG}
OLIGO_FILE_HTML=${OLIGO_DIR}/oligos_${CURRENT_SET}_${OLIGO_LEN}nt${STR}${NOOV}_sig${THOSIG}.html
oligos:
	@echo "${ORG}	${CURRENT_SET}	Detecting over-represented oligonucleotides	${OLIGO_FILE}"
	@mkdir -p ${OLIGO_DIR}
	oligo-analysis -v ${V}			\
		-i ${SEQ_FILE}			\
		-expfreq ${BG_OLIGO_FILE}	\
		-l ${OLIGO_LEN} ${NOOV}		\
		-return occ,freq,proba,rank	\
		-sort -thosig ${THOSIG}		\
		-seqtype dna			\
		-o ${OLIGO_FILE}
	text-to-html -i ${OLIGO_FILE} -o ${OLIGO_FILE_HTML}

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

#### Index fasta files
FASTA_FILES= ls -1 ${DATA}/${ORG}/*.fasta | grep -v purged
SEQ_LIST_FILE=${ORG}_files.txt
list_fasta_files:
	${FASTA_FILES}
	${FASTA_FILES} > ${SEQ_LIST_FILE}

MULTI_TASK=purge,oligos,dyads,maps,synthesis,sql
MULTI_DIR=${RES_DIR}/multi/${ORG}
MIN_OL=6
MAX_OL=6
NOOV=-noov
multi:
	@multiple-family-analysis -v ${V}				\
		-org ${ORG}						\
		-seq ${SEQ_LIST_FILE}					\
		-outdir ${MULTI_DIR}					\
		-2str -minol ${MIN_OL} -maxol ${MAX_OL}			\
		-bg upstream -sort name -task ${MULTI_TASK}		\
		${NOOV}							\
		-user jvanheld -password jvanheld -schema multifam


################################################################
## Calculate the effect of the number of sequences on mean and variance
SEQ_NUMBER_SERIES=1 2 3 4 5 6 7 8 9 10 15  20 30 40 50 60 70 80 90 100
seq_nb_series:
	fo*r nb in ${SEQ_NUMBER_SERIES} ; do			\
		${MAKE} calibrate_oligos N=$${nb};	\
	done

HUMAN_LENGTHS=500 1000 1500 2000 3000
seq_nb_series_human:
	${MAKE} seq_nb_series ORG=Homo_sapiens SEQ_LEN=1000 REPET=1000
	${MAKE} seq_nb_series ORG=Homo_sapiens SEQ_LEN=2000 REPET=1000

seq_nb_series_yeast:
	${MAKE} seq_nb_series ORG=Saccharomyces_cerevisiae SEQ_LEN=500 REPET=1000
	${MAKE} seq_nb_series ORG=Saccharomyces_cerevisiae SEQ_LEN=1000 REPET=1000

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
		${MAKE} calibrate_oligos SEQ_LEN=$${len} N=20;	\
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
SEQ_LEN=500
CALIB_TASK=all,clean_oligos
START=1
REPET=10000
CALIBRATE_CMD=							\
	calibrate-oligos.pl -v ${V}				\
		-r ${REPET} -sn ${N} -l ${OL} -sl ${SEQ_LEN}	\
		-task ${CALIB_TASK}				\
		-start ${START}					\
		${END}						\
		${STR} ${NOOV}					\
		-org ${ORG}

## Run the program immetiately (WHEN=now) or submit it to a queue (WHEN=queue)
WHEN=now
calibrate_oligos: calibrate_oligos_${WHEN}

## Run the program immediately
calibrate_oligos_now:
	${CALIBRATE_CMD}

## Submit the program to a queue
JOB_DIR=jobs
JOB=`mktemp job.XXXXXX`
calibrate_oligos_queue:
	@mkdir -p ${JOB_DIR}
	@for job in ${JOB} ; do						\
		echo "Job $${job}" ;					\
		echo "${CALIBRATE_CMD}" > ${JOB_DIR}/$${job} ;		\
		qsub -m e -q rsa@merlin.ulb.ac.be -N $${job} -j oe	\
			-o ${JOB_DIR}/$${job}.log ${JOB_DIR}/$${job} ;	\
	done

## A quick test for calibrate_oligos
calibrate_oligos_test:
	${MAKE} calibrate_oligos ORG=Mycoplasma_genitalium N=10 SEQ_LEN=200 STR=-1str NOOV=-ovlp R=10


give_access:
	dinf . -type d -exec chmod 775 {} \;
	dinf . -type f -exec chmod 664 {} \;


## ##############################################################
## Fit all previously calculated distributions with a Poisson and a
## negbin, respectively
DISTRIB_FILES=`find ${RES_DIR} -name '*_distrib.tab'` 
find_distrib_files:
	@echo ${DISTRIB_FILES}

check_distrib_files:
	@find ${RES_DIR} -name '*_distrib.tab' -exec wc {} \; | awk '$$1 < 2000'

all_fittings:
	@for file in ${DISTRIB_FILES} ; do			\
		${MAKE} one_fitting DISTRIB_FILE=$${file} ;	\
	done

one_fitting:
	${MAKE} one_fitting_one_law DISTRIB_LAW=poisson
	${MAKE} one_fitting_one_law DISTRIB_LAW=negbin

DISTRIB_LAW=negbin
DISTRIB_FILE=results/Saccharomyces_cerevisiae/rand_gene_selections/6nt-2str-noov_N10_L1000_R1000/Saccharomyces_cerevisiae_6nt_-2str-noov_n10_l1000_r1000_distrib.tab
FITTING_FILE=`echo ${DISTRIB_FILE} | perl -pe 's/distrib.tab/${DISTRIB_LAW}.tab/'` 
one_fitting_one_law:
	@echo "Fitting ${DISTRIB_LAW}	${FITTING_FILE}"
	@fit-distribution -v 1 -distrib ${DISTRIB_LAW} -i ${DISTRIB_FILE} -o ${FITTING_FILE} 
#	@head -50 ${DISTRIB_FILE}| fit-distribution -v 1 -distrib ${DISTRIB_LAW}  -o ${FITTING_FILE} 

