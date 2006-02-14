
################################################################
#
# Check the calibration of upstream frequencies for some organisms
# - position analysis of all promoters

include ${RSAT}/makefiles/util.mk

MAKEFILE=${RSAT}/makefiles/upstream_calibrations.mk

#### programs
#NICE=nice -n 20
#MAKE=${NICE} make -s
DATE=`date +%Y%m%d_%H%M%S`

COMPRESS=.gz

### verbosity
V=1
STR=-1str
NOOV=-noov

### Default upstream sequence length
#UP_LEN=1000
UP_LEN=800

#### organisms
#ORG=Mycoplasma_genitalium
#ORG=Homo_sapiens
ORG=Saccharomyces_cerevisiae
FUNGI= 	Saccharomyces_cerevisiae		\
	Schizosaccharomyces_pombe
#	Saccharomyces_cerevisiae_no_mito

PROCARYOTES= \
	Mycoplasma_genitalium			\
	Escherichia_coli_K12			\
	Bacillus_subtilis			\
	Salmonella_typhimurium_LT2

OTHER_ORGS= \
	Homo_sapiens				\
	Drosophila_melanogaster			\
	Caenorhabditis_elegans			\
	Arabidopsis_thaliana			\
	Plasmodium_falciparum
#	Mus_musculus				\

ORGANISMS=${FUNGI} ${PROCARYOTES} ${OTHER_ORGANISMS}


#### directories
RESULT=`pwd`/results
ORG_DIR=${RESULT}/${ORG}
POS_DIR=${ORG_DIR}/position_analysis
OLIGO_DIR=${ORG_DIR}/oligo-frequencies
MS_MARKOV_DIR=${ORG_DIR}/MS_Markov_models
RAND_FAM_DIR=${ORG_DIR}/random_gene_families
RAND_MULTI_DIR=${RAND_FAM_DIR}/rand_r${R}_n${N}

SEQ_PREFIX=${ORG}_allup${UP_LEN}${NOORF}
SEQ_FILE_UNCOMPRESSED=${ORG_DIR}/${SEQ_PREFIX}.fta
SEQ_FILE=${SEQ_FILE_UNCOMPRESSED}${COMPRESS}

dirs:
	mkdir -p ${ORG_DIR}
	mkdir -p ${POS_DIR}
	mkdir -p ${OLIGO_DIR}
	mkdir -p ${MS_MARKOV_DIR}
	mkdir -p ${OLIGO_DISTRIB_DIR}
	mkdir -p ${RAND_MULTI_DIR}


################################################################
#### Perform all the analyses
####
#### With the -noov option, the double strand counting is required,
#### since the overlapping reverse complement exclude each other, the
#### double strand counts is not the sum of the single strand count of
#### direct and reverse complement.

all:
	${MAKE} iterate_organisms TASK=all_tasks

all_tasks:  
	${MAKE} upstream_calibrations NOORF=${NOORF}
	${MAKE} upstream_calibrations NOORF=''
	${MAKE} random_families 

upstream_calibrations: seqs position_analysis  markov oligo_tables clean_seqs

################################################################
####
#### Iterators
####
# iterate_organisms:
# 	@for org in ${ORGANISMS} ; do			\
# 		echo "treating organism $${org}" ;	\
# 		${MAKE} ${TASK} ORG=$${org} ;		\
# 	done

FAMILY_SIZES= 20 50 10 5 100
iterate_family_sizes:
	@for n in ${FAMILY_SIZES}; do		\
		${MAKE} ${TASK} N=$${n};	\
	done

SEQ_LENGTHS=1000 2000 3000 1500 2500 500
iterate_seq_lengths:
	@for len in ${SEQ_LENGTHS}; do		\
		${MAKE} ${TASK} SEQ_LEN=$${len};	\
	done


################################################################
#### retrieve upstream sequences
#NOORF=-noorf
NOORF=
RETRIEVE_CMD=retrieve-seq -org ${ORG} -all -o ${SEQ_FILE} ${NOORF} -from -1 -to -${UP_LEN}
seqs: dirs
	@echo ""
	@echo "Retrieving upstream sequences	${RETRIEVE_CMD}"
	${RETRIEVE_CMD}

clean_seqs:
	\rm -f ${SEQ_FILE}

################################################################
#### position analysis
SEQ_TYPE=allup${NOORF}
POS_OL=6
CI=100
POS_FILE=${POS_DIR}/${ORG}_${SEQ_TYPE}_distrib_${POS_OL}nt_ci${CI}${STR}${NOOV}${COMPRESS}
POSITION_ANALYSIS_CMD=							\
	position-analysis -v ${V} -i ${SEQ_FILE} -format fasta -ci ${CI}	\
		-l ${POS_OL} ${STR} ${NOOV}				\
		-return distrib,exp,chi -origin -0 -nofilter		\
		-o ${POS_FILE}

one_position_analysis:
	@echo ${DATE} ${POSITION_ANALYSIS_CMD}
	@mkdir -p ${POS_DIR}
	@${MAKE} my_command MY_COMMAND="${POSITION_ANALYSIS_CMD}"

POS_OLIGO_LENGTHS=1 2 3 4 5 6
position_analysis:
	@for ol in ${POS_OLIGO_LENGTHS}; do			\
		${MAKE} one_position_analysis POS_OL=$${ol};	\
	done	


################################################################
#### gene per gene analysis of oligonucleotide frequencies
OL=6
OLIGO_FILE=${OLIGO_DIR}/${ORG}_allup${NOORF}_${OL}nt_mkv${MKV}${STR}${NOOV}${COMPRESS}
MKV=3
OLIGO_ANALYSIS_CMD=						\
	oligo-analysis -v ${V} -i ${SEQ_FILE} -o ${OLIGO_FILE}		\
		${NOORF} ${STR} ${NOOV} -l ${OL} -markov ${MKV}	\
		-return occ,freq,proba,zscore,mseq,rank -sort
one_markov: dirs
	@echo ${DATE}	${OLIGO_ANALYSIS_CMD}
	@${MAKE} my_command MY_COMMAND="${OLIGO_ANALYSIS_CMD}"

markov_one_ol:
	@for m in ${MARKOV_ORDERS} ;do		\
		${MAKE} one_markov MKV=$${m};	\
	done 

markov:
	${MAKE} markov_one_ol OL=2 MARKOV_ORDERS='0' 
	${MAKE} markov_one_ol OL=3 MARKOV_ORDERS='0 1'
	${MAKE} markov_one_ol OL=4 MARKOV_ORDERS='0 1'
	${MAKE} markov_one_ol OL=5 MARKOV_ORDERS='0 1 2'
	${MAKE} markov_one_ol OL=6 MARKOV_ORDERS='0 1 2 3'
	${MAKE} markov_one_ol OL=7 MARKOV_ORDERS='0 1 2 3 4'
	${MAKE} markov_one_ol OL=8 MARKOV_ORDERS='0 1 2 3 4 5'


################################################################
## Calculate a table with the occurrences of each nucleotide (columns)
## in each upstream sequence (rows)

FROM=-1
TO=-250
OLIGO_TABLE=${OLIGO_DIR}/upstream_${FROM}_${TO}_per_gene_${OL}nt${STR}${NOOV}.tab${COMPRESS}

OLIGO_TABLE_CMD=sub-sequence -i ${SEQ_FILE} -from ${FROM} -to ${TO} | oligo-analysis -v ${V} -format fasta -table -l ${OL} ${STR} ${NOOV} -o ${OLIGO_TABLE}

one_oligo_table:
	mkdir -p ${OLIGO_DIR}
	@echo
	@echo ${DATE}
	@echo "${OLIGO_TABLE_CMD}"
	@${MAKE} my_command MY_COMMAND="${OLIGO_TABLE_CMD}"

oligo_windows:
	${MAKE} one_oligo_table FROM=-0001 TO=-0250
	${MAKE} one_oligo_table FROM=-0251 TO=-0500
	${MAKE} one_oligo_table FROM=-0501 TO=-0750
	${MAKE} one_oligo_table FROM=-0751 TO=-1000

	${MAKE} one_oligo_table FROM=-1001 TO=-1250
	${MAKE} one_oligo_table FROM=-1251 TO=-1500
	${MAKE} one_oligo_table FROM=-1501 TO=-1750
	${MAKE} one_oligo_table FROM=-1751 TO=-2000

	${MAKE} one_oligo_table FROM=-2001 TO=-2250
	${MAKE} one_oligo_table FROM=-2251 TO=-2500
	${MAKE} one_oligo_table FROM=-2501 TO=-2750
	${MAKE} one_oligo_table FROM=-2751 TO=-3000

	${MAKE} one_oligo_table FROM=-3001 TO=-3250
	${MAKE} one_oligo_table FROM=-3251 TO=-3500
	${MAKE} one_oligo_table FROM=-3501 TO=-3750
	${MAKE} one_oligo_table FROM=-3751 TO=-4000

	${MAKE} one_oligo_table FROM=-4001 TO=-4250
	${MAKE} one_oligo_table FROM=-4251 TO=-4500
	${MAKE} one_oligo_table FROM=-4501 TO=-4750
	${MAKE} one_oligo_table FROM=-4751 TO=-5000

OLIGO_LENGTHS=1 2 3 4 5 6
oligo_tables:
	@for ol in ${OLIGO_LENGTHS}; do			\
		${MAKE} oligo_windows OL=$${ol};	\
	done


################################################################
#
# Analyze random gene families to check the rate of false positive
#

random_families_all_organisms:
	${MAKE} iterate_organisms TASK=random_families

random_families:
	@${MAKE} iterate_family_sizes TASK=random_families_one_size
	@${MAKE} index_results
	${MAKE} rand_stats

random_families_one_size:
	@${MAKE} dirs
	@${MAKE} rand_fam
	@${MAKE} rand_multi

N=8
R=10000
RAND_FAM_FILE=${RAND_FAM_DIR}/rand_r${R}_n${N}.fam 
rand_fam:
	random-genes -org ${ORG} -n ${N} -r ${R} -o ${RAND_FAM_FILE}

#MAPS=,maps
MULTI_TASKS=upstream,oligos,synthesis,sql,clean${MAPS}
rand_multi: dirs
	multiple-family-analysis -v ${V} -org ${ORG}		\
		-i ${RAND_FAM_FILE}				\
		-outdir ${RAND_MULTI_DIR} -task ${MULTI_TASKS}



################################################################
#
# Index the results
#

INDEX_FILE=${RESULT}/index.html
index_results:
	@echo '<html><body><table>' > ${INDEX_FILE} 
	@${MAKE} iterate_organisms TASK=index_organism 
	@echo '</table></html></body>' >> ${INDEX_FILE} 

index_organism:
	@echo 'Indexing organism ${ORG}'
	@echo '<tr>' >> ${INDEX_FILE}
	@echo '<th><a href=${ORG}>${ORG}</a></th>' >> ${INDEX_FILE}
	@${MAKE} iterate_family_sizes TASK=index_rand_multi
	@echo '</tr>' >> ${INDEX_FILE}

SYNTH_TABLE=${ORG}/random_gene_families/rand_r${R}_n${N}/synthetic_tables/
index_rand_multi:
	@echo '<td><a href=${SYNTH_TABLE}>n${N} r${R}</a></td>' \
		>> ${INDEX_FILE}


################################################################
#
# Analyze the result of pattern discovery in random gene sets

rand_stats:
	@${MAKE} iterate_organisms TASK=sig_distrib
	@${MAKE} iterate_organisms TASK=sig_distrib_top
	@${MAKE} compare_sig_distrib TOP=_top
	@${MAKE} compare_sig_distrib TOP=''

################################################################
#### calculate the distribution of significance index

#### number of patterns per sig
PATTERN_TABLE=${RAND_MULTI_DIR}/sql_export/Pattern.tab
SIG_DISTRIB_DIR=${RAND_MULTI_DIR}/stats
SIG_DISTRIB_TOP=${SIG_DISTRIB_DIR}/sig_distrib_top
SIG_DISTRIB=${SIG_DISTRIB_DIR}/sig_distrib
sig_distrib:
	@echo 'Calculating sig distrib	${ORG}'
	@mkdir -p ${SIG_DISTRIB_DIR}
	@cat ${PATTERN_TABLE}			\
		| grep -v  '^--'		\
		| cut -f 11			\
		| classfreq -v -ci 0.5 -min 0	\
		> ${SIG_DISTRIB}.tab

#### number of families per sig
sig_distrib_top:
	@echo 'Calculating top sig distrib	${ORG}'
	@mkdir -p ${SIG_DISTRIB_DIR}
	@cat ${PATTERN_TABLE}			\
		| grep -v  '^--'		\
		| awk '$$12 == 1'		\
		| cut -f 11			\
		| classfreq -v -ci 0.5 -min 0	\
		> ${SIG_DISTRIB_TOP}.tab

################################################################
#### compare sig distributions between genomes
TMP=tmp
TOP=_top
COMPA_DIR=${RESULT}/comparisons/rand_R${R}_N${N}
COMPA_FILE=${COMPA_DIR}/sig_distrib_comparison${TOP}
compare_sig_distrib:
	${MAKE} compare_sig_distrib_table
	${MAKE} compare_sig_distrib_graph
	${MAKE} compare_sig_distrib_graph YLOG='-ylog'

compare_sig_distrib_table:
	@echo 'Comparing sig distributions	${TOP}	${COMPA_FILE}'
	@mkdir -p ${COMPA_DIR}
	@\rm -rf ${TMP}
	@for org in ${ORGANISMS} ; do		\
		${MAKE} tmp_copy ORG=$${org} ;	\
	done

	compare-scores -null 0 -sc 6 -files `ls -t ${TMP}/*`	\
		| perl -pe 's|_1||g'				\
		| perl -pe 's|${TMP}/||g'			\
		> ${TMP}/comparison.tab

	grep '^;' ${TMP}/comparison.tab		\
		> ${COMPA_FILE}.tab
	grep -v '^;' ${TMP}/comparison.tab	\
		| sort -n			\
		>> ${COMPA_FILE}.tab

	text-to-html -i ${COMPA_FILE}.tab -font variable	\
		| perl -pe 's/_/ /g'				\
		>  ${COMPA_FILE}.html
	@\rm -rf ${TMP}

compare_sig_distrib_graph:
	XYgraph -i ${COMPA_FILE}.tab				\
		-o ${COMPA_FILE}${YLOG}.gif				\
		-xcol 1 -ycol '2-10'				\
		${YLOG}						\
		-xsize 1000					\
		-title1 'sig distributions ${TOP}'		\
		-xleg1 'sig index [sig = -log(E_value)]'	\
		-yleg1 'number of patterns with sig >= x'	\
		-xgstep2 0.5 -ygstep2 10			\
		-lines -legend -header


tmp_copy:
	@mkdir -p ${TMP}
	cp  ${SIG_DISTRIB}${TOP}.tab ${TMP}/${ORG}


################################################################
## Calculate the distribution of occurrences for each oligonucleotide
## in all upstream sequencces
OLIGO_DISTRIB_DIR=${ORG_DIR}/oligo_distrib
OLIGO_DISTRIB_FILE=${OLIGO_DISTRIB_DIR}/${SEQ_PREFIX}_${OL}nt${STR}${NOOV}

OLIGO_DISTRIB_CMD= \
	oligo-analysis -v ${V} -i ${SEQ_FILE} ${STR} ${NOOV} -l ${OL} -distrib -o ${OLIGO_DISTRIB_FILE} ; \
	fit-distribution			\
		-v ${V} -distrib poisson \
		-i ${OLIGO_DISTRIB_FILE}		\
		-o ${OLIGO_DISTRIB_FILE}_poisson.tab ; \
	fit-distribution			\
		-v ${V} -distrib negbin \
		-i ${OLIGO_DISTRIB_FILE}		\
		-o ${OLIGO_DISTRIB_FILE}_negbin.tab


oligo_distrib_one_org:
	${MAKE} seqs
	${MAKE} iterate_oligo_lengths OLIGO_TASK=oligo_distrib_one_org_one_size

oligo_distrib_one_org_one_size: 
	@echo ""
	@echo "Analyzing oligo distributions for one organism	${ORG}	${OL}	${UP_LEN}	${NOORF}"
	@${MAKE} my_command MY_COMMAND="${OLIGO_DISTRIB_CMD}"
	@echo OLIGO_DISTRIB_FILE ${OLIGO_DISTRIB_FILE}
	@echo "Done oligo distributions for one organism	${ORG}	${OL}	${UP_LEN}	${NOORF}"


################################################################
## Organism-specific settings
oligo_distrib_all_organisms: oligo_distrib_fungi oligo_distrib_procaryotes oligo_distrib_other_orgs

GROUP=${FUNGI}
oligo_distrib_one_group:
	@echo "Analysing group	${GROUP}	${UP_LEN}"
	for org in ${GROUP}; do 			\
		${MAKE} NOORF='' seqs ORG=$${org} ; 		\
		${MAKE} NOORF='' oligo_distrib_one_org ORG=$${org} ; 	\
	done
	@echo "Done group	${GROUP}	${UP_LEN}"

oligo_distrib_fungi:
	${MAKE} oligo_distrib_one_group GROUP="${FUNGI}" UP_LEN=800

oligo_distrib_procaryotes:
	${MAKE} oligo_distrib_one_group GROUP="${PROCARYOTES}" UP_LEN=200
	${MAKE} oligo_distrib_one_group GROUP="${PROCARYOTES}" UP_LEN=400

oligo_distrib_other_orgs:
	${MAKE} oligo_distrib_one_group GROUP="${OTHER_ORGS}" UP_LEN=200
	${MAKE} oligo_distrib_one_group GROUP="${OTHER_ORGS}" UP_LEN=500
	${MAKE} oligo_distrib_one_group GROUP="${OTHER_ORGS}" UP_LEN=1000
	${MAKE} oligo_distrib_one_group GROUP="${OTHER_ORGS}" UP_LEN=2000

################################################################
## Create background models for MotifSampler
MARKOV_ORDERS=0 1 2 3 4 5
MS_markov_models:
	${MAKE} seqs COMPRESS=''
	@for m in ${MARKOV_ORDERS} ;do		\
		${MAKE} MS_one_markov_model MKV=$${m} COMPRESS='';	\
	done 
	rm -f ${SEQ_FILE_UNCOMPRESSED}

MS_MARKOV_FILE=${MS_MARKOV_DIR}/MS_mkv${MKV}_${SEQ_PREFIX}.txt
MS_one_markov_model: dirs
	CreateBackgroundModel -n ${ORG} -f ${SEQ_FILE} -o ${MKV} -b ${MS_MARKOV_FILE}

################################################################
## Drosophila analysis


################################################################
## Create background models for MEME

## TO BE DONE
