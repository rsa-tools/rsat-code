################################################################
## makefile for getting info for the motif discovery competition


################################################################
#### General parameters

## level of verbosity
V=1

################################################################
### commands
MAKEFILE=${RSAT}/makefiles/motif_disco_competition_2003.mk
MAKE=make -f ${MAKEFILE}

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
#### publish
SERVER_LOCATION=rubens.ulb.ac.be:/rubens/dsk3/genomics/motif_discovery_competition
publish:
	${RSYNC} --exclude '*~' . ${SERVER_LOCATION}

update_from_server:
	${RSYNC} --exclude '*~' ${SERVER_LOCATION}/'*' .

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
ORG_SETS=`ls -1 ${DATA}/${ORG}/*.fasta | perl -pe 's|${DATA}/${ORG}/||g'  | perl -pe 's|\.fasta||g'`
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
SEQ_LEN=`sequence-lengths -i ${SEQ_FILE} | cut -f 2 | sort -u`
SEQ_NB=`sequence-lengths -i ${SEQ_FILE} | wc -l | awk '{print $$1}'`
SEQ_LEN_FILE=${CURRENT_RES_DIR}/${CURRENT_SET}_lengths.txt
seq_len:  current_dir
	@echo "set ${CURRENT_SET} length ${SEQ_LEN}	${SEQ_NB} seqs	${SEQ_LEN_FILE}"
	@echo ${SEQ_LEN} > ${SEQ_LEN_FILE}

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
	@echo "	${MAKE} calibrate_oligos ORG=${ORG} SEQ_LEN=${SEQ_LEN} N=${SEQ_NB} STR=-2str STR=-2str NOOV=-noov" >> ${CALIB_SCRIPT}


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
FASTA_FILES= ls -1 ${DATA}/${ORG}/*.fasta
SEQ_LIST_FILE=${ORG}_files.txt
list_fasta_files:
	${FASTA_FILES}
	${FASTA_FILES} > ${SEQ_LIST_FILE}

MULTI_TASK=oligos
MULTI_DIR=${RES_DIR}/multi/${ORG}
MIN_OL=6
MAX_OL=6
NOOV=-noov
multi:
	@multiple-family-analysis -v ${V}					\
		-org ${ORG}						\
		-seq ${SEQ_LIST_FILE}					\
		-outdir ${MULTI_DIR}					\
		-2str -minol ${MIN_OL} -maxol ${MAX_OL}			\
		-bg upstream -sort name -task ${MULTI_TASK}		\
		${NOOV}							\
		-user jvanheld -password jvanheld -schema multifam

