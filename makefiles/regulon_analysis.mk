################################################################
## Discovery patterns in yeast regulons

include ${RSAT}/makefiles/util.mk

WHEN=now
MAKEFILE=${RSAT}/makefiles/regulon_analysis.mk
V=1
WORK_DIR=`pwd`

################################################################
## queue tasks on merlin

generate_all_random:
	for family in ${ALL_RANDOM}; do			\
		${MAKE} generate_random RAND_FAM_FILE=$${family};			\
	done

#### beware, to change, purge-sequence doesnt use unique identifiers ####
#### temporary files used by it have the same name ####
upstream_regulons:
	for family in ${REGULONS}; do				\
		${MAKE} multi MULTI_TASK=upstream,purge FAM=${WORK_DIR}/$${family} WHEN=queue;	\
	done

upstream_random:
	for family in ${ALL_RANDOM}; do				\
		${MAKE} multi MULTI_TASK=upstream,purge FAM=${WORK_DIR}/$${family} WHEN=queue;	\
	done

MWIDTH_LIST = 6 8 10
MATRIX_TASKS = gibbs consensus meme
analyse_matrix:
	for mwidth in ${MWIDTH_LIST}; do \
		for task in ${MATRIX_TASKS}; do							\
			for family in ${ALL_RANDOM} ${REGULONS}; do				\
				${MAKE} multi MULTI_TASK=$${task} FAM=${WORK_DIR}/$${family} MATRIX_WIDTH=$${mwidth} WHEN=queue;	\
			done;									\
		done;										\
QUEUE_TASKS = oligos oligo_maps gibbs consensus MotifSampler dyads dyad_maps
queue_all:
#	for family in ${REGULONS} ${RAND_FAM}; do				\
#		${MAKE} multi MULTI_TASK=upstream,purge FAM=${WORK_DIR}/$${family} WHEN=now;	\
#	done
	for task in ${QUEUE_TASKS}; do							\
		for family in ${REGULONS} ${RAND_FAM}; do				\
			${MAKE} multi MULTI_TASK=$${task} FAM=$${family} WHEN=queue;	\
		done;									\
	done

RSAT_TASKS=oligos dyads
analyse_rsat:
	for task in ${RSAT_TASKS}; do			\
		for family in ${ALL_RANDOM} ${REGULONS}; do		\
			${MAKE} multi MULTI_TASK=$${task} FAM=${WORK_DIR}/$${family} WHEN=queue;	\
		done;														\
	done

MAX_OL_LIST=6 8
AFTER_TASKS=oligo_maps,merge_oligos,dyad_maps,synthesis,sql
after_queue:
	for family in ${ALL_RANDOM} ${REGULONS}; do					\
		for maxol in ${MAX_OL_LIST}; do \
			${MAKE} multi MULTI_TASK=${AFTER_TASKS} FAM=${WORK_DIR}/$${family} MAX_OL=$${maxol} WHEN=queue;	\
		done; \
	done

################################################################
## Data
################################################################

ALL_RANDOM=random_selection1.fam random_selection2.fam random_selection3.fam random_selection4.fam random_selection5.fam	\
random_selection6.fam random_selection7.fam random_selection8.fam random_selection9.fam random_selection10.fam

################################################################
### generate random gene selection
################################################################
## Generate random gene selections

REGULONS=regulons_TF_aMAZE
REGULONS_FILE=data/${REGULONS}.fam
#RAND_FAM=random_selection1
RAND_FAM=random_selection
RAND_FAM_FILE=data/${RAND_FAM}.fam
#REGULONS=regulons_TF_aMAZE
#RAND_FAM=random_selection
#REGULONS_FILE=data/${REGULONS}.fam
#RAND_FAM_FILE=data/${RAND_FAM}.fam
rand_fam:
	random-genes -fam ${REGULONS_FILE} -o ${RAND_FAM_FILE}	\
		 -features data/Feature_nomito.tab -org ${ORG}

generate_random: rand_fam

#	random-genes -fam ${REGULONS_FILE} -o data/${RAND_FAM_FILE} -features data/Feature_nomito.tab -org ${ORG};

################################################################
#### Parameters for pattern discovery (oligo-analysis and
#### dyad-analysis)
FAM=${REGULONS}
#FAM_FILE=data/${FAMS}
FAM_FILE=${WORK_DIR}/data/${FAM}.fam

MULTI_INPUT=-i ${FAM_FILE}
ORG=Saccharomyces_cerevisiae
STR=-2str
THOSIG=0
NOOV=-noov
MULTI_TASK=upstream,purge,oligos,merge_oligos,oligo_maps,dyads,dyad_maps,consensus,gibbs,MotifSampler,meme,synthesis,sql
NOORF=-noorf
MULTI_BG=upstream${NOORF}
MULTI_EXP=-bg ${MULTI_BG}
PURGE=-purge
RES_DIR=${WORK_DIR}/results
ORG_DIR=${RES_DIR}/${ORG}
MULTI_DIR=${ORG_DIR}/${FAM}/multi/${MULTI_BG}_bg${PURGE}
MIN_OL=6
MAX_OL=8
MIN_SP=0
MAX_SP=20
MATRIX_WIDTH=10
SPS=2
SORT=score
MULTI_OPT=
SKIP=
KNOWN=-known data/factor_site.tab
MULTI_CMD=multiple-family-analysis -v ${V}					\
		${MULTI_INPUT}							\
		-thosig ${THOSIG}						\
		${PURGE} ${NOORF}							\
		-org ${ORG}							\
		-outdir ${MULTI_DIR}						\
		${STR}								\
		-minol ${MIN_OL} -maxol ${MAX_OL}				\
		-minsp ${MIN_SP} -maxsp ${MAX_SP}				\
		-width ${MATRIX_WIDTH} -sps ${SPS}				\
		${MULTI_EXP}							\
		-sort ${SORT} -task ${MULTI_TASK}				\
		${NOOV} ${MULTI_OPT}						\
		${SKIP} ${KNOWN}							

list_directories:
	@echo working	${WORK_DIR}
	@echo results	${RES_DIR}
	@echo organism	${ORG_DIR}
	@echo multi	${MULTI_DIR}

## call multiple-family-analysis
multi:
	@echo ${MULTI_CMD}
	${MAKE} my_command MY_COMMAND="${MULTI_CMD}"


SERVER=merlin.bigre.ulb.ac.be
LOGIN=jvanheld
SERVER_DIR=research/seq_analysis/yeast/regul_fam/regulons_Nicolas_2004
SERVER_LOCATION=${LOGIN}@${SERVER}:${SERVER_DIR}
TO_SYNC=results
TARGET_DIR=.
TO_SERVER_CMD=${RSYNC} ${TO_SYNC} ${SERVER_LOCATION}
to_server:
	@echo ${TO_SERVER_CMD}
	${TO_SERVER_CMD}

FROM_SERVER_CMD=${RSYNC} ${SERVER_LOCATION}${TO_SYNC} ${TARGET_DIR}
from_server:
	@mkdir -p ${TARGET_DIR}
	@echo ${FROM_SERVER_CMD}
	${FROM_SERVER_CMD}

################################################################
# Synchronize results of discriminant analysis
#
DISCRIM_DIR=/home/nicolas/discriminant_analysis/gctfya_discr_an
get_discrim:
	${MAKE} from_server TO_SYNC=${DISCRIM_DIR} SERVER_DIR='' TARGET_DIR=results/discriminant_analysis/gctfya_discr_an

################################################################
# Compare results obtained with different pattern discovery programs
#
CRITERIA=ln.exp exp ln.Pval Pval unadjusted.information adjusted.information MAP model.map betaprior.map
compare:
	@${MAKE} compare_one_criterion PROGRAM=oligo-analysis CRITERION=max.score CI=0.5 COMPA_MIN=0  COMPA_SC=9 COMPA_MAX=auto

	@${MAKE} compare_one_criterion PROGRAM=gibbs CRITERION=MAP CI=5 COMPA_MIN=-60  COMPA_SC=9 COMPA_MAX=140
	@${MAKE} compare_one_criterion PROGRAM=gibbs CRITERION=model.map CI=5 COMPA_MIN=0 COMPA_SC=9 COMPA_MAX=auto
	@${MAKE} compare_one_criterion PROGRAM=gibbs CRITERION=betaprior.map CI=5 COMPA_MIN=-100 COMPA_SC=9 COMPA_MAX=auto

	@${MAKE} compare_one_criterion PROGRAM=consensus CRITERION=ln.exp CI=1 COMPA_MIN=-70 COMPA_SC=8 COMPA_MAX=1
	@${MAKE} compare_one_criterion PROGRAM=consensus CRITERION=ln.Pval CI=1 COMPA_MIN=-70 COMPA_SC=8 COMPA_MAX=1
	@${MAKE} compare_one_criterion PROGRAM=consensus CRITERION=adjusted.information CI=1 COMPA_MIN=0 COMPA_MAX=auto
	@${MAKE} compare_one_criterion PROGRAM=consensus CRITERION=unadjusted.information CI=1 COMPA_MIN=0 COMPA_MAX=auto



CRITERION=max.score
PROGRAM=oligo-analysis
CI=1
COMPA_DIR=results/comparisons
COMPA_MIN=0
COMPA_MAX=auto
COMPA_SC=9
XGSTEP2=1
REGULON_SYNTHESIS=${ORG_DIR}/${REGULONS}/multi/upstream_bg-purge/synthetic_tables/${REGULONS}_bg_upstream_up-800_-1_noorf-purge_6nt_8nt-noov-2str_sig0.html
RAND_FAM_SYNTHESIS=${ORG_DIR}/${RAND_FAM}/multi/upstream_bg-purge/synthetic_tables/${RAND_FAM}_bg_upstream_up-800_-1_noorf-purge_6nt_8nt-noov-2str_sig0.html
compare_one_criterion:
	@mkdir -p ${COMPA_DIR}
	@echo "comparing regulons and random selections for criterion ${CRITERION}"
	@echo ";" ${REGULONS} > ${COMPA_DIR}/${REGULONS}.${PROGRAM}.${CRITERION}
	@cat ${REGULON_SYNTHESIS} | grep ">${CRITERION}<"			\
		| perl -pe 's|<tr><td align=right>${CRITERION}</td><td>||g'	\
		|perl -pe 's|</td></tr>||'					\
		| sort -nr							\
		| classfreq -ci ${CI} -min ${COMPA_MIN} -max ${COMPA_MAX} -v	\
		>> ${COMPA_DIR}/${REGULONS}.${PROGRAM}.${CRITERION}

	@echo ";" ${RAND_FAM} > ${COMPA_DIR}/${RAND_FAM}.${PROGRAM}.${CRITERION}
	@cat ${RAND_FAM_SYNTHESIS}						\
		| grep ">${CRITERION}<"						\
		| perl -pe 's|<tr><td align=right>${CRITERION}</td><td>||g'	\
		| perl -pe 's|</td></tr>||'					\
		| sort -nr							\
		| classfreq -ci ${CI} -min ${COMPA_MIN} -max ${COMPA_MAX} -v	\
		>> ${COMPA_DIR}/${RAND_FAM}.${PROGRAM}.${CRITERION} 

	@compare-scores -null 0						\
		-i ${COMPA_DIR}/${REGULONS}.${PROGRAM}.${CRITERION}	\
		-i ${COMPA_DIR}/${RAND_FAM}.${PROGRAM}.${CRITERION}	\
		-ic 3 -sc ${COMPA_SC} -numeric				\
		| perl -pe 's|${COMPA_DIR}/||g'				\
		> ${COMPA_DIR}/compa.${PROGRAM}.${CRITERION}

	XYgraph											\
		-i ${COMPA_DIR}/compa.${PROGRAM}.${CRITERION}					\
		-o ${COMPA_DIR}/compa.${PROGRAM}.${CRITERION}.jpg				\
		-xmin ${COMPA_MIN}								\
		-xmax ${COMPA_MAX}								\
		-format jpg -xcol 1 -ycol 2,3							\
		-ymin 0 -ymax 1									\
		-ygstep1 0.1 -ygstep2 0.02							\
		-xgstep2 ${XGSTEP2}								\
		-title1 "${PROGRAM} ${CRITERION}" -xleg1 ${CRITERION}				\
		-title2 'Comparison between annotated regulons and random gene selections'	\
		-yleg1 "Number of families" -legend -lines -xsize 800 -ysize 400



	XYgraph											\
		-i ${COMPA_DIR}/compa.${PROGRAM}.${CRITERION}					\
		-o ${COMPA_DIR}/compa.${PROGRAM}.${CRITERION}_rand_vs_reg.jpg			\
		-xmin 0								\
		-xmax 1								\
		-format jpg -xcol 3 -ycol 2							\
		-ymin 0 -ymax 1									\
		-ygstep1 0.1 -ygstep2 0.02							\
		-xgstep1 0.1 -xgstep2 0.02								\
		-title1 "${PROGRAM} ${CRITERION}"						\
		-yleg1 "regulons"								\
		-title2 'Comparison between annotated regulons and random gene selections'	\
		-xleg1 "random selections"							\
		 -legend -lines -xsize 800 -ysize 400


##################### NICOLAS' SCRIPTS #####################################################################################
## Summary of the results : retrieve best scores, group random selections, calc classfreq and compare random and regulons ##
############################################################################################################################
CLUSTER_DIR=${ORG_DIR}/home/jvanheld/regulons
ALL_RANDOM_DIR=${ORG_DIR}/all_random
MATRIX_SERIE=6-8-10

summarize:
	${MAKE} get_all_best_scores
	${MAKE} group_all_random
	${MAKE} get_all_classfreq
	${MAKE} compare_all

############## retrieve best scores for each analysis ##############

get_all_best_scores:
	@echo "CLUSTER_DIR="${CLUSTER_DIR}
	for fam in ${REGULONS} ${ALL_RANDOM}; do			\
		for maxol in ${MAX_OL_LIST}; do \
			${MAKE} init_rsat_scores MULTI_DIR=${CLUSTER_DIR}/$${fam}/multi/${MULTI_BG}_bg${PURGE} MAX_OL=$${maxol}; \
			${MAKE} get_rsat_scores MULTI_DIR=${CLUSTER_DIR}/$${fam}/multi/${MULTI_BG}_bg${PURGE} \
			CURRENT_LIST="`cut -f 2 $${fam} | grep -v ^% | sort -u | perl -pe s/"\n"/" "/ | perl -pe s/' $$'//`" MAX_OL=$${maxol}; \
		done; \
		for mwidth in ${MWIDTH_LIST}; do		\
			${MAKE} init_matrix_scores MULTI_DIR=${CLUSTER_DIR}/$${fam}/multi/${MULTI_BG}_bg${PURGE} MATRIX_WIDTH=$${mwidth}; \
			${MAKE} get_matrix_scores MULTI_DIR=${CLUSTER_DIR}/$${fam}/multi/${MULTI_BG}_bg${PURGE} \
			CURRENT_LIST="`cut -f 2 $${fam} | grep -v ^% | sort -u | perl -pe s/"\n"/" "/ | perl -pe s/' $$'//`" MATRIX_WIDTH=$${mwidth}; \
		done; \
		${MAKE} init_matrix_scores MATRIX_WIDTH=${MATRIX_SERIE} MULTI_DIR=${CLUSTER_DIR}/$${fam}/multi/${MULTI_BG}_bg${PURGE};	\
		${MAKE} get_all_matrix_widths MULTI_DIR=${CLUSTER_DIR}/$${fam}/multi/${MULTI_BG}_bg${PURGE}	\
		CURRENT_LIST="`cut -f 2 $${fam} | grep -v ^% | sort -u | perl -pe s/"\n"/" "/ | perl -pe s/' $$'//`";	\
	done

get_all_matrix_widths:
	for famm in ${CURRENT_LIST}; do \
		${MAKE} group_matrix_scores NFAM=$${famm};	\
	done

group_matrix_scores:
	 echo ${NFAM}"	"`grep '^#' ${MULTI_DIR}/${NFAM}/MotifSampler_${NFAM}/${NFAM}_up800_noorf-s1-p0.2-n4-w*-x1-r5.sites | sort -bgu +13 | \
	 tail -1 | awk '{print $$14}'` >> ${MULTI_DIR}/synthetic_tables/motifsampler_up800_noorf-s1-p0.2-n4-w${MATRIX_SERIE}-x1-r5_best_LL.scores
	 echo ${NFAM}"	"`grep '^#' ${MULTI_DIR}/${NFAM}/MotifSampler_${NFAM}/${NFAM}_up800_noorf-s1-p0.2-n4-w*-x1-r5.sites | sort -bgu +11 | \
	 tail -1 | awk '{print $$12}'` >> ${MULTI_DIR}/synthetic_tables/motifsampler_up800_noorf-s1-p0.2-n4-w${MATRIX_SERIE}-x1-r5_best_IC.scores
	 echo ${NFAM}"	"`grep '^#' ${MULTI_DIR}/${NFAM}/MotifSampler_${NFAM}/${NFAM}_up800_noorf-s1-p0.2-n4-w*-x1-r5.sites | sort -bgu +9 | \
	 tail -1 | awk '{print $$10}'` >> ${MULTI_DIR}/synthetic_tables/motifsampler_up800_noorf-s1-p0.2-n4-w${MATRIX_SERIE}-x1-r5_best_CS.scores
	 echo ${NFAM}"	"`grep '^model map' ${MULTI_DIR}/${NFAM}/gibbs_${NFAM}/${NFAM}_up800_noorf-L*-n*-d-n | cut -f 4 -d ' ' | sort -gu | tail -1 | \
	 perl -pe s/';$$'//` >> ${MULTI_DIR}/synthetic_tables/gibbs_up800_noorf-L${MATRIX_SERIE}-d-n_best_MM.scores
	 echo ${NFAM}"	"`grep '^model map' ${MULTI_DIR}/${NFAM}/gibbs_${NFAM}/${NFAM}_up800_noorf-L*-n*-d-n | cut -f 8 -d ' ' | sort -gru | tail -1`	\
	 >> ${MULTI_DIR}/synthetic_tables/gibbs_up800_noorf-L${MATRIX_SERIE}-d-n_best_BP.scores
	 echo ${NFAM}"	"`grep '^ MAP' ${MULTI_DIR}/${NFAM}/gibbs_${NFAM}/${NFAM}_up800_noorf-L*-n*-d-n | awk '{print $$4}' | sort -gu | tail -1` \
	 >> ${MULTI_DIR}/synthetic_tables/gibbs_up800_noorf-L${MATRIX_SERIE}-d-n_best_MAP.scores
	 echo ${NFAM}"	"`grep ^MOTIF ${MULTI_DIR}/${NFAM}/meme_${NFAM}/${NFAM}_up800_noorf_tcm-w*-2str | awk '{print $$11}' | sort -gu | tail -1` \
	 >> ${MULTI_DIR}/synthetic_tables/meme_up800_noorf_tcm-w${MATRIX_SERIE}-2str_best_LLR.scores
	 echo ${NFAM}"	"`grep ^MOTIF ${MULTI_DIR}/${NFAM}/meme_${NFAM}/${NFAM}_up800_noorf_tcm-w*-2str | awk '{print log($$14)}' | sort -gru | tail -1` \
	 >> ${MULTI_DIR}/synthetic_tables/meme_up800_noorf_tcm-w${MATRIX_SERIE}-2str_best_Eval.scores
	 echo ${NFAM}"	"`grep '^unadjusted information' ${MULTI_DIR}/${NFAM}/consensus_${NFAM}/${NFAM}_up800_noorf-c2-L*-n* | cut -d ' ' -f 4 | \
	 sort -gu | tail -1` >> ${MULTI_DIR}/synthetic_tables/consensus_up800_noorf-c2-L${MATRIX_SERIE}_best_UI.scores
	 echo ${NFAM}"	"`grep '^sample size adjusted information' ${MULTI_DIR}/${NFAM}/consensus_${NFAM}/${NFAM}_up800_noorf-c2-L*-n* | cut -d ' ' -f 6 | \
	 sort -gu | tail -1` >> ${MULTI_DIR}/synthetic_tables/consensus_up800_noorf-c2-L${MATRIX_SERIE}_best_AI.scores
	 echo ${NFAM}"	"`grep '^ln(e-value)' ${MULTI_DIR}/${NFAM}/consensus_${NFAM}/${NFAM}_up800_noorf-c2-L*-n* | cut -d ' ' -f 3 | sort -gru | \
	 tail -1` >> ${MULTI_DIR}/synthetic_tables/consensus_up800_noorf-c2-L${MATRIX_SERIE}_best_lnEF.scores
	 echo ${NFAM}"	"`grep '^ln(p-value)' ${MULTI_DIR}/${NFAM}/consensus_${NFAM}/${NFAM}_up800_noorf-c2-L*-n* | cut -d ' ' -f 3 | sort -gru | \
	 tail -1` >> ${MULTI_DIR}/synthetic_tables/consensus_up800_noorf-c2-L${MATRIX_SERIE}_best_lnPval.scores

init_rsat_scores:
	@echo "selecting best scores in "${MULTI_DIR}" oligos "${MIN_OL}"_"${MAX_OL}
	@echo "Set	Oligo-analysis Sig" > ${MULTI_DIR}/synthetic_tables/oligos_up800_noorf_bg_upstream-noorf_${MIN_OL}-${MAX_OL}nt${STR}_sig${THOSIG}${NOOV}_best.scores
	@echo "Set	Dyad-analysis Sig" > ${MULTI_DIR}/synthetic_tables/dyads_up800_noorf_bg_upstream-noorf_l3_sp${MIN_SP}-${MAX_SP}${STR}_any_sig${THOSIG}${NOOV}_best.scores
	@echo "Set	Oligos+dyads Sig" > ${MULTI_DIR}/synthetic_tables/oligos_dyads_up800_noorf_bg_upstream-noorf_${MIN_OL}-${MAX_OL}nt_l3_sp${MIN_SP}-${MAX_SP}${STR}_sig${THOSIG}${NOOV}_best.scores
init_matrix_scores:
	@echo "selecting best scores in "${MULTI_DIR}" matrix width "${MATRIX_WIDTH}
	@echo "Set	Motif Sampler LL" > ${MULTI_DIR}/synthetic_tables/motifsampler_up800_noorf-s1-p0.2-n4-w${MATRIX_WIDTH}-x1-r5_best_LL.scores
	@echo "Set	Motif Sampler IC" > ${MULTI_DIR}/synthetic_tables/motifsampler_up800_noorf-s1-p0.2-n4-w${MATRIX_WIDTH}-x1-r5_best_IC.scores
	@echo "Set	Motif Sampler CS" > ${MULTI_DIR}/synthetic_tables/motifsampler_up800_noorf-s1-p0.2-n4-w${MATRIX_WIDTH}-x1-r5_best_CS.scores
	@echo "Set	Gibbs Sampler MM" > ${MULTI_DIR}/synthetic_tables/gibbs_up800_noorf-L${MATRIX_WIDTH}-d-n_best_MM.scores
	@echo "Set	Gibbs Sampler BP" > ${MULTI_DIR}/synthetic_tables/gibbs_up800_noorf-L${MATRIX_WIDTH}-d-n_best_BP.scores
	@echo "Set	Gibbs Sampler MAP" > ${MULTI_DIR}/synthetic_tables/gibbs_up800_noorf-L${MATRIX_WIDTH}-d-n_best_MAP.scores
	@echo "Set	MEME LLR" > ${MULTI_DIR}/synthetic_tables/meme_up800_noorf_tcm-w${MATRIX_WIDTH}-2str_best_LLR.scores
	@echo "Set	MEME lnEval" > ${MULTI_DIR}/synthetic_tables/meme_up800_noorf_tcm-w${MATRIX_WIDTH}-2str_best_Eval.scores
	@echo "Set	Consensus UI" > ${MULTI_DIR}/synthetic_tables/consensus_up800_noorf-c2-L${MATRIX_WIDTH}_best_UI.scores
	@echo "Set	Consensus AI" > ${MULTI_DIR}/synthetic_tables/consensus_up800_noorf-c2-L${MATRIX_WIDTH}_best_AI.scores
	@echo "Set	Consensus lnEF" > ${MULTI_DIR}/synthetic_tables/consensus_up800_noorf-c2-L${MATRIX_WIDTH}_best_lnEF.scores
	@echo "Set	Consensus lnPval" > ${MULTI_DIR}/synthetic_tables/consensus_up800_noorf-c2-L${MATRIX_WIDTH}_best_lnPval.scores


get_rsat_scores:
	for famm in ${CURRENT_LIST}; do				\
		${MAKE} select_rsat_best_scores NFAM=$${famm};		\
	done

get_matrix_scores:
	for famm in ${CURRENT_LIST}; do				\
		${MAKE} select_matrix_best_scores NFAM=$${famm};	\
	done

select_rsat_best_scores:
	 echo ${NFAM}"	"`grep -v '^;' ${MULTI_DIR}/${NFAM}/oligos_${NFAM}/${NFAM}_up800_noorf_oligos_bg_upstream-noorf_${MIN_OL}-${MAX_OL}nt${STR}_sig${THOSIG}${NOOV} | \
	 cut -f 8 | sort -gu | tail -1` >> ${MULTI_DIR}/synthetic_tables/oligos_up800_noorf_bg_upstream-noorf_${MIN_OL}-${MAX_OL}nt${STR}_sig${THOSIG}${NOOV}_best.scores
	 echo ${NFAM}"	"`grep -v '^;' ${MULTI_DIR}/${NFAM}/dyads_${NFAM}/${NFAM}_up800_noorf_dyads_bg_upstream-noorf_l3_sp${MIN_SP}-${MAX_SP}${STR}_any_sig${THOSIG}${NOOV} | \
	 cut -f 8 | perl -pe s/'   '// | sort -gu | tail -1` >> ${MULTI_DIR}/synthetic_tables/dyads_up800_noorf_bg_upstream-noorf_l3_sp${MIN_SP}-${MAX_SP}${STR}_any_sig${THOSIG}${NOOV}_best.scores
	 mkdir -p ${MULTI_DIR}/synthetic_tables/temp
	 echo ${NFAM}"	"`grep -v '^;' ${MULTI_DIR}/${NFAM}/oligos_${NFAM}/${NFAM}_up800_noorf_oligos_bg_upstream-noorf_${MIN_OL}-${MAX_OL}nt${STR}_sig${THOSIG}${NOOV} | cut -f 8 \
	 | sort -gu | tail -1` > ${MULTI_DIR}/synthetic_tables/temp/${NFAM}_oligos_dyads_best_sig.tmp
	 echo ${NFAM}"	"`grep -v '^;' ${MULTI_DIR}/${NFAM}/dyads_${NFAM}/${NFAM}_up800_noorf_dyads_bg_upstream-noorf_l3_sp${MIN_SP}-${MAX_SP}${STR}_any_sig${THOSIG}${NOOV} | cut -f 8 \
	 | perl -pe s/'   '// | sort -gu | tail -1` >> ${MULTI_DIR}/synthetic_tables/temp/${NFAM}_oligos_dyads_best_sig.tmp
	 echo ${NFAM}"	"`cut -f 2 ${MULTI_DIR}/synthetic_tables/temp/${NFAM}_oligos_dyads_best_sig.tmp | sort -gu | tail -1` >> \
	 ${MULTI_DIR}/synthetic_tables/oligos_dyads_up800_noorf_bg_upstream-noorf_${MIN_OL}-${MAX_OL}nt_l3_sp${MIN_SP}-${MAX_SP}${STR}_sig${THOSIG}${NOOV}_best.scores
	 rm ${MULTI_DIR}/synthetic_tables/temp/${NFAM}_oligos_dyads_best_sig.tmp

select_matrix_best_scores:

	 echo ${NFAM}"	"`grep '^#' ${MULTI_DIR}/${NFAM}/MotifSampler_${NFAM}/${NFAM}_up800_noorf-s1-p0.2-n4-w${MATRIX_WIDTH}-x1-r5.sites | sort -bgu +13 | \
	 tail -1 | awk '{print $$14}'` >> ${MULTI_DIR}/synthetic_tables/motifsampler_up800_noorf-s1-p0.2-n4-w${MATRIX_WIDTH}-x1-r5_best_LL.scores
	 echo ${NFAM}"	"`grep '^#' ${MULTI_DIR}/${NFAM}/MotifSampler_${NFAM}/${NFAM}_up800_noorf-s1-p0.2-n4-w${MATRIX_WIDTH}-x1-r5.sites | sort -bgu +11 | \
	 tail -1 | awk '{print $$12}'` >> ${MULTI_DIR}/synthetic_tables/motifsampler_up800_noorf-s1-p0.2-n4-w${MATRIX_WIDTH}-x1-r5_best_IC.scores
	 echo ${NFAM}"	"`grep '^#' ${MULTI_DIR}/${NFAM}/MotifSampler_${NFAM}/${NFAM}_up800_noorf-s1-p0.2-n4-w${MATRIX_WIDTH}-x1-r5.sites | sort -bgu +9 | \
	 tail -1 | awk '{print $$10}'` >> ${MULTI_DIR}/synthetic_tables/motifsampler_up800_noorf-s1-p0.2-n4-w${MATRIX_WIDTH}-x1-r5_best_CS.scores
	 echo ${NFAM}"	"`grep '^model map' ${MULTI_DIR}/${NFAM}/gibbs_${NFAM}/${NFAM}_up800_noorf-L${MATRIX_WIDTH}-n*-d-n | cut -f 4 -d ' ' | sort -gu | tail -1 | \
	 perl -pe s/';$$'//` >> ${MULTI_DIR}/synthetic_tables/gibbs_up800_noorf-L${MATRIX_WIDTH}-d-n_best_MM.scores
	 echo ${NFAM}"	"`grep '^model map' ${MULTI_DIR}/${NFAM}/gibbs_${NFAM}/${NFAM}_up800_noorf-L${MATRIX_WIDTH}-n*-d-n | cut -f 8 -d ' ' | sort -gru | tail -1`	\
	 >> ${MULTI_DIR}/synthetic_tables/gibbs_up800_noorf-L${MATRIX_WIDTH}-d-n_best_BP.scores
	 echo ${NFAM}"	"`grep '^ MAP' ${MULTI_DIR}/${NFAM}/gibbs_${NFAM}/${NFAM}_up800_noorf-L${MATRIX_WIDTH}-n*-d-n | awk '{print $$3}' | sort -gu | tail -1` \
	 >> ${MULTI_DIR}/synthetic_tables/gibbs_up800_noorf-L${MATRIX_WIDTH}-d-n_best_MAP.scores
	 echo ${NFAM}"	"`grep ^MOTIF ${MULTI_DIR}/${NFAM}/meme_${NFAM}/${NFAM}_up800_noorf_tcm-w${MATRIX_WIDTH}-2str | awk '{print $$11}' | sort -gu | tail -1` \
	 >> ${MULTI_DIR}/synthetic_tables/meme_up800_noorf_tcm-w${MATRIX_WIDTH}-2str_best_LLR.scores
	 echo ${NFAM}"	"`grep ^MOTIF ${MULTI_DIR}/${NFAM}/meme_${NFAM}/${NFAM}_up800_noorf_tcm-w${MATRIX_WIDTH}-2str | awk '{print log($$14)}' | sort -gru | tail -1` \
	 >> ${MULTI_DIR}/synthetic_tables/meme_up800_noorf_tcm-w${MATRIX_WIDTH}-2str_best_Eval.scores
	 echo ${NFAM}"	"`grep '^unadjusted information' ${MULTI_DIR}/${NFAM}/consensus_${NFAM}/${NFAM}_up800_noorf-c2-L${MATRIX_WIDTH}-n* | cut -d ' ' -f 4 | \
	 sort -gu | tail -1` >> ${MULTI_DIR}/synthetic_tables/consensus_up800_noorf-c2-L${MATRIX_WIDTH}_best_UI.scores
	 echo ${NFAM}"	"`grep '^sample size adjusted information' ${MULTI_DIR}/${NFAM}/consensus_${NFAM}/${NFAM}_up800_noorf-c2-L${MATRIX_WIDTH}-n* | cut -d ' ' -f 6 | \
	 sort -gu | tail -1` >> ${MULTI_DIR}/synthetic_tables/consensus_up800_noorf-c2-L${MATRIX_WIDTH}_best_AI.scores
	 echo ${NFAM}"	"`grep '^ln(e-value)' ${MULTI_DIR}/${NFAM}/consensus_${NFAM}/${NFAM}_up800_noorf-c2-L${MATRIX_WIDTH}-n* | cut -d ' ' -f 3 | sort -gru | \
	 tail -1` >> ${MULTI_DIR}/synthetic_tables/consensus_up800_noorf-c2-L${MATRIX_WIDTH}_best_lnEF.scores
	 echo ${NFAM}"	"`grep '^ln(p-value)' ${MULTI_DIR}/${NFAM}/consensus_${NFAM}/${NFAM}_up800_noorf-c2-L${MATRIX_WIDTH}-n* | cut -d ' ' -f 3 | sort -gru | \
	 tail -1` >> ${MULTI_DIR}/synthetic_tables/consensus_up800_noorf-c2-L${MATRIX_WIDTH}_best_lnPval.scores




######################################
###### group all random results ######
######################################
group_all_random:
	rm -rf ${ALL_RANDOM_DIR}
	for maxol in ${MAX_OL_LIST}; do \
		${MAKE} group_random_rsat MAX_OL=$${maxol}; \
	done
	for mwidth in ${MWIDTH_LIST}; do \
		${MAKE} group_random_matrix MATRIX_WIDTH=$${mwidth}; \
	done
	${MAKE} group_random_matrix MATRIX_WIDTH=${MATRIX_SERIE}

group_random_rsat:
	mkdir -p ${ALL_RANDOM_DIR}
	touch ${ALL_RANDOM_DIR}/oligos_up800_noorf_bg_upstream-noorf_${MIN_OL}-${MAX_OL}nt${STR}_sig${THOSIG}${NOOV}_all_random_best.scores
	touch ${ALL_RANDOM_DIR}/dyads_up800_noorf_bg_upstream-noorf_l3_sp${MIN_SP}-${MAX_SP}${STR}_any_sig${THOSIG}${NOOV}_all_random_best.scores
	touch ${ALL_RANDOM_DIR}/oligos_dyads_up800_noorf_bg_upstream-noorf_${MIN_OL}-${MAX_OL}nt_l3_sp${MIN_SP}-${MAX_SP}${STR}_sig${THOSIG}${NOOV}_all_random_best.scores
	for randfam in ${ALL_RANDOM}; do			\
		grep -v '^Set' ${CLUSTER_DIR}/$${randfam}/multi/${MULTI_BG}_bg${PURGE}/synthetic_tables/oligos_up800_noorf_bg_upstream-noorf_${MIN_OL}-${MAX_OL}nt${STR}_sig${THOSIG}${NOOV}_best.scores \
		>> ${ALL_RANDOM_DIR}/oligos_up800_noorf_bg_upstream-noorf_${MIN_OL}-${MAX_OL}nt${STR}_sig${THOSIG}${NOOV}_all_random_best.scores;	\
		grep -v '^Set' ${CLUSTER_DIR}/$${randfam}/multi/${MULTI_BG}_bg${PURGE}/synthetic_tables/dyads_up800_noorf_bg_upstream-noorf_l3_sp${MIN_SP}-${MAX_SP}${STR}_any_sig${THOSIG}${NOOV}_best.scores \
		>> ${ALL_RANDOM_DIR}/dyads_up800_noorf_bg_upstream-noorf_l3_sp${MIN_SP}-${MAX_SP}${STR}_any_sig${THOSIG}${NOOV}_all_random_best.scores;	\
		grep -v '^Set' ${CLUSTER_DIR}/$${randfam}/multi/${MULTI_BG}_bg${PURGE}/synthetic_tables/oligos_dyads_up800_noorf_bg_upstream-noorf_${MIN_OL}-${MAX_OL}nt_l3_sp${MIN_SP}-${MAX_SP}${STR}_sig${THOSIG}${NOOV}_best.scores \
		>> ${ALL_RANDOM_DIR}/oligos_dyads_up800_noorf_bg_upstream-noorf_${MIN_OL}-${MAX_OL}nt_l3_sp${MIN_SP}-${MAX_SP}${STR}_sig${THOSIG}${NOOV}_all_random_best.scores;	\
	done

group_random_matrix:
	touch ${ALL_RANDOM_DIR}/motifsampler_up800_noorf-s1-p0.2-n4-w${MATRIX_WIDTH}-x1-r5_all_random_best_LL.scores
	touch ${ALL_RANDOM_DIR}/motifsampler_up800_noorf-s1-p0.2-n4-w${MATRIX_WIDTH}-x1-r5_all_random_best_IC.scores
	touch ${ALL_RANDOM_DIR}/motifsampler_up800_noorf-s1-p0.2-n4-w${MATRIX_WIDTH}-x1-r5_all_random_best_CS.scores
	touch ${ALL_RANDOM_DIR}/gibbs_up800_noorf-L${MATRIX_WIDTH}-d-n_all_random_best_MM.scores
	touch ${ALL_RANDOM_DIR}/gibbs_up800_noorf-L${MATRIX_WIDTH}-d-n_all_random_best_BP.scores
	touch ${ALL_RANDOM_DIR}/gibbs_up800_noorf-L${MATRIX_WIDTH}-d-n_all_random_best_MAP.scores
	touch ${ALL_RANDOM_DIR}/meme_up800_noorf_tcm-w${MATRIX_WIDTH}-2str_all_random_best_LLR.scores
	touch ${ALL_RANDOM_DIR}/meme_up800_noorf_tcm-w${MATRIX_WIDTH}-2str_all_random_best_Eval.scores
	touch ${ALL_RANDOM_DIR}/consensus_up800_noorf-c2-L${MATRIX_WIDTH}_all_random_best_UI.scores
	touch ${ALL_RANDOM_DIR}/consensus_up800_noorf-c2-L${MATRIX_WIDTH}_all_random_best_AI.scores
	touch ${ALL_RANDOM_DIR}/consensus_up800_noorf-c2-L${MATRIX_WIDTH}_all_random_best_lnEF.scores
	touch ${ALL_RANDOM_DIR}/consensus_up800_noorf-c2-L${MATRIX_WIDTH}_all_random_best_lnPval.scores

	for randfam in ${ALL_RANDOM}; do			\
		grep -v '^Set' ${CLUSTER_DIR}/$${randfam}/multi/${MULTI_BG}_bg${PURGE}/synthetic_tables/motifsampler_up800_noorf-s1-p0.2-n4-w${MATRIX_WIDTH}-x1-r5_best_LL.scores \
		>> ${ALL_RANDOM_DIR}/motifsampler_up800_noorf-s1-p0.2-n4-w${MATRIX_WIDTH}-x1-r5_all_random_best_LL.scores;	\
		grep -v '^Set' ${CLUSTER_DIR}/$${randfam}/multi/${MULTI_BG}_bg${PURGE}/synthetic_tables/motifsampler_up800_noorf-s1-p0.2-n4-w${MATRIX_WIDTH}-x1-r5_best_IC.scores \
		>> ${ALL_RANDOM_DIR}/motifsampler_up800_noorf-s1-p0.2-n4-w${MATRIX_WIDTH}-x1-r5_all_random_best_IC.scores;	\
		grep -v '^Set' ${CLUSTER_DIR}/$${randfam}/multi/${MULTI_BG}_bg${PURGE}/synthetic_tables/motifsampler_up800_noorf-s1-p0.2-n4-w${MATRIX_WIDTH}-x1-r5_best_CS.scores \
		>> ${ALL_RANDOM_DIR}/motifsampler_up800_noorf-s1-p0.2-n4-w${MATRIX_WIDTH}-x1-r5_all_random_best_CS.scores;	\
		grep -v '^Set' ${CLUSTER_DIR}/$${randfam}/multi/${MULTI_BG}_bg${PURGE}/synthetic_tables/gibbs_up800_noorf-L${MATRIX_WIDTH}-d-n_best_MM.scores \
		>> ${ALL_RANDOM_DIR}/gibbs_up800_noorf-L${MATRIX_WIDTH}-d-n_all_random_best_MM.scores;	\
		grep -v '^Set' ${CLUSTER_DIR}/$${randfam}/multi/${MULTI_BG}_bg${PURGE}/synthetic_tables/gibbs_up800_noorf-L${MATRIX_WIDTH}-d-n_best_BP.scores \
		>> ${ALL_RANDOM_DIR}/gibbs_up800_noorf-L${MATRIX_WIDTH}-d-n_all_random_best_BP.scores;	\
		grep -v '^Set' ${CLUSTER_DIR}/$${randfam}/multi/${MULTI_BG}_bg${PURGE}/synthetic_tables/gibbs_up800_noorf-L${MATRIX_WIDTH}-d-n_best_MAP.scores \
		>> ${ALL_RANDOM_DIR}/gibbs_up800_noorf-L${MATRIX_WIDTH}-d-n_all_random_best_MAP.scores;	\
		grep -v '^Set' ${CLUSTER_DIR}/$${randfam}/multi/${MULTI_BG}_bg${PURGE}/synthetic_tables/meme_up800_noorf_tcm-w${MATRIX_WIDTH}-2str_best_LLR.scores \
		>> ${ALL_RANDOM_DIR}/meme_up800_noorf_tcm-w${MATRIX_WIDTH}-2str_all_random_best_LLR.scores;	\
		grep -v '^Set' ${CLUSTER_DIR}/$${randfam}/multi/${MULTI_BG}_bg${PURGE}/synthetic_tables/meme_up800_noorf_tcm-w${MATRIX_WIDTH}-2str_best_Eval.scores \
		>> ${ALL_RANDOM_DIR}/meme_up800_noorf_tcm-w${MATRIX_WIDTH}-2str_all_random_best_Eval.scores;	\
		grep -v '^Set' ${CLUSTER_DIR}/$${randfam}/multi/${MULTI_BG}_bg${PURGE}/synthetic_tables/consensus_up800_noorf-c2-L${MATRIX_WIDTH}_best_UI.scores \
		>> ${ALL_RANDOM_DIR}/consensus_up800_noorf-c2-L${MATRIX_WIDTH}_all_random_best_UI.scores;	\
		grep -v '^Set' ${CLUSTER_DIR}/$${randfam}/multi/${MULTI_BG}_bg${PURGE}/synthetic_tables/consensus_up800_noorf-c2-L${MATRIX_WIDTH}_best_AI.scores \
		>> ${ALL_RANDOM_DIR}/consensus_up800_noorf-c2-L${MATRIX_WIDTH}_all_random_best_AI.scores;	\
		grep -v '^Set' ${CLUSTER_DIR}/$${randfam}/multi/${MULTI_BG}_bg${PURGE}/synthetic_tables/consensus_up800_noorf-c2-L${MATRIX_WIDTH}_best_lnEF.scores \
		>> ${ALL_RANDOM_DIR}/consensus_up800_noorf-c2-L${MATRIX_WIDTH}_all_random_best_lnEF.scores;	\
		grep -v '^Set' ${CLUSTER_DIR}/$${randfam}/multi/${MULTI_BG}_bg${PURGE}/synthetic_tables/consensus_up800_noorf-c2-L${MATRIX_WIDTH}_best_lnPval.scores \
		>> ${ALL_RANDOM_DIR}/consensus_up800_noorf-c2-L${MATRIX_WIDTH}_all_random_best_lnPval.scores;	\
	done

#########################################################################################
###### calculate class frequencies and compare regulons with all random selections ######
#########################################################################################

COMP_DIR_NS=results/comp_ns
MULTI_REG_DIR=${ORG_DIR}/home/jvanheld/regulons/${REGULONS}/multi/${MULTI_BG}_bg${PURGE}/synthetic_tables

get_all_classfreq:
	mkdir -p ${COMP_DIR_NS}
	for maxol in ${MAX_OL_LIST}; do \
		${MAKE} get_rsat_classfreq MULTI_DIR=${ALL_RANDOM_DIR} CURSET=all_random MAX_OL=$${maxol} CUREXT=all_random_; \
		${MAKE} get_rsat_classfreq MULTI_DIR=${MULTI_REG_DIR} CURSET=${REGULONS} MAX_OL=$${maxol} CUREXT=""; \
	done
	for mwidth in ${MWIDTH_LIST}; do		\
		${MAKE} get_matrix_classfreq MULTI_DIR=${ALL_RANDOM_DIR} CURSET=all_random MATRIX_WIDTH=$${mwidth} CUREXT=all_random_; \
		${MAKE} get_matrix_classfreq MULTI_DIR=${MULTI_REG_DIR} CURSET=${REGULONS} MATRIX_WIDTH=$${mwidth} CUREXT=""; \
	done
	${MAKE} get_matrix_classfreq MULTI_DIR=${ALL_RANDOM_DIR} CURSET=all_random MATRIX_WIDTH=${MATRIX_SERIE} CUREXT=all_random_; \
	${MAKE} get_matrix_classfreq MULTI_DIR=${MULTI_REG_DIR} CURSET=${REGULONS} MATRIX_WIDTH=${MATRIX_SERIE} CUREXT=""; \

#### IMPORTANT NOTE ####
## RSAT AND Gibbs Sampler can return empty score (RSAT as been set to a threshold of 0 ###
## in these case, for oligo and dyads, empty has been set to 0 ##
## for Gibbs Beta prior, empty has been set to 0 ##
## for Gibbs Model Map, empty has been set to -1E10 ##
## for Gibbs MAP, empty has been set to 

get_rsat_classfreq:
	perl -pe s/'\t$$'/'\t0'/ ${MULTI_DIR}/oligos_up800_noorf_bg_upstream-noorf_${MIN_OL}-${MAX_OL}nt${STR}_sig${THOSIG}${NOOV}_${CUREXT}best.scores | cut -f 2 | \
	classfreq -ci 0.05 -min 0 | cut -f 3,9 > ${COMP_DIR_NS}/oligos_${MIN_OL}-${MAX_OL}nt_${CURSET}_classfreq.tab
	perl -pe s/'\t$$'/'\t0'/ ${MULTI_DIR}/dyads_up800_noorf_bg_upstream-noorf_l3_sp${MIN_SP}-${MAX_SP}${STR}_any_sig${THOSIG}${NOOV}_${CUREXT}best.scores | cut -f 2 | \
	classfreq -ci 0.05 -min 0 | cut -f 3,9 > ${COMP_DIR_NS}/dyads_sp${MIN_SP}-${MAX_SP}_${CURSET}_classfreq.tab
	perl -pe s/'\t$$'/'\t0'/ ${MULTI_DIR}/oligos_dyads_up800_noorf_bg_upstream-noorf_${MIN_OL}-${MAX_OL}nt_l3_sp${MIN_SP}-${MAX_SP}${STR}_sig${THOSIG}${NOOV}_${CUREXT}best.scores | \
	cut -f 2 | classfreq -ci 0.05 -min 0 | cut -f 3,9 > ${COMP_DIR_NS}/oligos_dyads_${MIN_OL}-${MAX_OL}nt_sp${MIN_SP}-${MAX_SP}_${CURSET}_classfreq.tab

get_matrix_classfreq:
	perl -pe s/'\t$$'/'\t0'/ ${MULTI_DIR}/motifsampler_up800_noorf-s1-p0.2-n4-w${MATRIX_WIDTH}-x1-r5_${CUREXT}best_LL.scores | perl -pe s/'\t-inf$$'/'\t0'/ | cut -f 2 | \
	classfreq -ci 1 -min 0 | cut -f 3,9 > ${COMP_DIR_NS}/motifsampler_w${MATRIX_WIDTH}_LL_${CURSET}_classfreq.tab
	perl -pe s/'\t$$'/'\t1'/ ${MULTI_DIR}/motifsampler_up800_noorf-s1-p0.2-n4-w${MATRIX_WIDTH}-x1-r5_${CUREXT}best_IC.scores | cut -f 2 | \
	classfreq -ci 0.01 -min 1 | cut -f 3,9 > ${COMP_DIR_NS}/motifsampler_w${MATRIX_WIDTH}_IC_${CURSET}_classfreq.tab
	perl -pe s/'\t$$'/'\t1'/ ${MULTI_DIR}/motifsampler_up800_noorf-s1-p0.2-n4-w${MATRIX_WIDTH}-x1-r5_${CUREXT}best_CS.scores | cut -f 2 | \
	classfreq -ci 0.01 -min 1 | cut -f 3,9 > ${COMP_DIR_NS}/motifsampler_w${MATRIX_WIDTH}_CS_${CURSET}_classfreq.tab
	perl -pe s/'\t$$'/'\t-1e-10'/ ${MULTI_DIR}/gibbs_up800_noorf-L${MATRIX_WIDTH}-d-n_${CUREXT}best_MM.scores | cut -f 2  | \
	classfreq -ci 2 -min -2 | cut -f 3,9 > ${COMP_DIR_NS}/gibbs_L${MATRIX_WIDTH}_MM_${CURSET}_classfreq.tab
	perl -pe s/'\t$$'/'\t0'/ ${MULTI_DIR}/gibbs_up800_noorf-L${MATRIX_WIDTH}-d-n_${CUREXT}best_BP.scores | cut -f 2 | \
	classfreq -ci 2 -min -400 | cut -f 3,9 > ${COMP_DIR_NS}/gibbs_L${MATRIX_WIDTH}_BP_${CURSET}_classfreq.tab
	perl -pe s/'\t$$'/'\t-300'/ ${MULTI_DIR}/gibbs_up800_noorf-L${MATRIX_WIDTH}-d-n_${CUREXT}best_MAP.scores | cut -f 2 | \
	classfreq -ci 2 -min -300 | cut -f 3,9 > ${COMP_DIR_NS}/gibbs_L${MATRIX_WIDTH}_MAP_${CURSET}_classfreq.tab
	perl -pe s/'\t$$'/'\t0'/ ${MULTI_DIR}/meme_up800_noorf_tcm-w${MATRIX_WIDTH}-2str_${CUREXT}best_LLR.scores | cut -f 2 | \
	classfreq -ci 2 -min 0 | cut -f 3,9 > ${COMP_DIR_NS}/meme_w${MATRIX_WIDTH}_LLR_${CURSET}_classfreq.tab
	perl -pe s/'\t$$'/'\t20'/ ${MULTI_DIR}/meme_up800_noorf_tcm-w${MATRIX_WIDTH}-2str_${CUREXT}best_Eval.scores | cut -f 2 | \
	classfreq -ci 0.1 -min -80 | cut -f 3,9 > ${COMP_DIR_NS}/meme_w${MATRIX_WIDTH}_Eval_${CURSET}_classfreq.tab
	cut -f 2 ${MULTI_DIR}/consensus_up800_noorf-c2-L${MATRIX_WIDTH}_${CUREXT}best_UI.scores | \
	classfreq -ci 0.01 -min 6 | cut -f 3,9 > ${COMP_DIR_NS}/consensus_L${MATRIX_WIDTH}_${CURSET}_UI_classfreq.tab
	cut -f 2 ${MULTI_DIR}/consensus_up800_noorf-c2-L${MATRIX_WIDTH}_${CUREXT}best_AI.scores | \
	classfreq -ci 0.01 -min 6 | cut -f 3,9 > ${COMP_DIR_NS}/consensus_L${MATRIX_WIDTH}_${CURSET}_AI_classfreq.tab
	cut -f 2 ${MULTI_DIR}/consensus_up800_noorf-c2-L${MATRIX_WIDTH}_${CUREXT}best_lnEF.scores | \
	classfreq -ci 0.1 -min -80 | cut -f 3,9 > ${COMP_DIR_NS}/consensus_L${MATRIX_WIDTH}_${CURSET}_lnEF_classfreq.tab
	cut -f 2 ${MULTI_DIR}/consensus_up800_noorf-c2-L${MATRIX_WIDTH}_${CUREXT}best_lnPval.scores | \
	classfreq -ci 1 -min -600 | cut -f 3,9 > ${COMP_DIR_NS}/consensus_L${MATRIX_WIDTH}_${CURSET}_lnPval_classfreq.tab

RAND_SET=all_random
COMP_ID=comparison

compare_all:
	for maxol in ${MAX_OL_LIST}; do \
		${MAKE} compare_rsat MAX_OL=$${maxol}; \
	done
	for mwidth in ${MWIDTH_LIST}; do		\
		${MAKE} compare_matrix MATRIX_WIDTH=$${mwidth}; \
	done
	${MAKE} compare_matrix MATRIX_WIDTH=${MATRIX_SERIE}

compare_rsat:
	@compare-scores -null 0						\
		-i ${COMP_DIR_NS}/oligos_${MIN_OL}-${MAX_OL}nt_${REGULONS}_classfreq.tab	\
		-i ${COMP_DIR_NS}/oligos_${MIN_OL}-${MAX_OL}nt_${RAND_SET}_classfreq.tab	\
		-numeric > ${COMP_DIR_NS}/oligos_${MIN_OL}-${MAX_OL}nt_${COMP_ID}.tab

	@compare-scores -null 0						\
		-i ${COMP_DIR_NS}/dyads_sp${MIN_SP}-${MAX_SP}_${REGULONS}_classfreq.tab	\
		-i ${COMP_DIR_NS}/dyads_sp${MIN_SP}-${MAX_SP}_${RAND_SET}_classfreq.tab \
		-numeric > ${COMP_DIR_NS}/dyads_sp${MIN_SP}-${MAX_SP}_${COMP_ID}.tab

	@compare-scores -null 0						\
		-i ${COMP_DIR_NS}/oligos_dyads_${MIN_OL}-${MAX_OL}nt_sp${MIN_SP}-${MAX_SP}_${REGULONS}_classfreq.tab \
		-i ${COMP_DIR_NS}/oligos_dyads_${MIN_OL}-${MAX_OL}nt_sp${MIN_SP}-${MAX_SP}_${RAND_SET}_classfreq.tab \
		-numeric > ${COMP_DIR_NS}/oligos_dyads_${MIN_OL}-${MAX_OL}nt_sp${MIN_SP}-${MAX_SP}_${COMP_ID}.tab

compare_matrix:

	@compare-scores -null 0						\
		-i ${COMP_DIR_NS}/motifsampler_w${MATRIX_WIDTH}_LL_${REGULONS}_classfreq.tab	\
		-i ${COMP_DIR_NS}/motifsampler_w${MATRIX_WIDTH}_LL_${RAND_SET}_classfreq.tab	\
		-numeric > ${COMP_DIR_NS}/motifsampler_w${MATRIX_WIDTH}_LL_${COMP_ID}.tab

	@compare-scores -null 0						\
		-i ${COMP_DIR_NS}/motifsampler_w${MATRIX_WIDTH}_IC_${REGULONS}_classfreq.tab	\
		-i ${COMP_DIR_NS}/motifsampler_w${MATRIX_WIDTH}_IC_${RAND_SET}_classfreq.tab	\
		-numeric > ${COMP_DIR_NS}/motifsampler_w${MATRIX_WIDTH}_IC_${COMP_ID}.tab

	@compare-scores -null 0						\
		-i ${COMP_DIR_NS}/motifsampler_w${MATRIX_WIDTH}_CS_${REGULONS}_classfreq.tab	\
		-i ${COMP_DIR_NS}/motifsampler_w${MATRIX_WIDTH}_CS_${RAND_SET}_classfreq.tab	\
		-numeric > ${COMP_DIR_NS}/motifsampler_w${MATRIX_WIDTH}_CS_${COMP_ID}.tab

	@compare-scores -null 0						\
		-i ${COMP_DIR_NS}/gibbs_L${MATRIX_WIDTH}_MM_${REGULONS}_classfreq.tab	\
		-i ${COMP_DIR_NS}/gibbs_L${MATRIX_WIDTH}_MM_${RAND_SET}_classfreq.tab	\
		-numeric > ${COMP_DIR_NS}/gibbs_L${MATRIX_WIDTH}_MM_${COMP_ID}.tab

	@compare-scores -null 0						\
		-i ${COMP_DIR_NS}/gibbs_L${MATRIX_WIDTH}_BP_${REGULONS}_classfreq.tab \
		-i ${COMP_DIR_NS}/gibbs_L${MATRIX_WIDTH}_BP_${RAND_SET}_classfreq.tab	\
		-numeric > ${COMP_DIR_NS}/gibbs_L${MATRIX_WIDTH}_BP_${COMP_ID}.tab

	@compare-scores -null 0						\
		-i ${COMP_DIR_NS}/gibbs_L${MATRIX_WIDTH}_MAP_${REGULONS}_classfreq.tab \
		-i ${COMP_DIR_NS}/gibbs_L${MATRIX_WIDTH}_MAP_${RAND_SET}_classfreq.tab \
		-numeric > ${COMP_DIR_NS}/gibbs_L${MATRIX_WIDTH}_MAP_${COMP_ID}.tab

	@compare-scores -null 0						\
		-i ${COMP_DIR_NS}/meme_w${MATRIX_WIDTH}_LLR_${REGULONS}_classfreq.tab \
		-i ${COMP_DIR_NS}/meme_w${MATRIX_WIDTH}_LLR_${RAND_SET}_classfreq.tab \
		-numeric > ${COMP_DIR_NS}/meme_w${MATRIX_WIDTH}_LLR_${COMP_ID}.tab

	@compare-scores -null 0						\
		-i ${COMP_DIR_NS}/meme_w${MATRIX_WIDTH}_Eval_${REGULONS}_classfreq.tab \
		-i ${COMP_DIR_NS}/meme_w${MATRIX_WIDTH}_Eval_${RAND_SET}_classfreq.tab \
		-numeric > ${COMP_DIR_NS}/meme_w${MATRIX_WIDTH}_Eval_${COMP_ID}.tab

	@compare-scores -null 0 						\
		-i ${COMP_DIR_NS}/consensus_L${MATRIX_WIDTH}_${REGULONS}_UI_classfreq.tab \
		-i ${COMP_DIR_NS}/consensus_L${MATRIX_WIDTH}_${RAND_SET}_UI_classfreq.tab \
		-numeric > ${COMP_DIR_NS}/consensus_L${MATRIX_WIDTH}_UI_${COMP_ID}.tab

	@compare-scores -null 0 						\
		-i ${COMP_DIR_NS}/consensus_L${MATRIX_WIDTH}_${REGULONS}_AI_classfreq.tab \
		-i ${COMP_DIR_NS}/consensus_L${MATRIX_WIDTH}_${RAND_SET}_AI_classfreq.tab \
		-numeric > ${COMP_DIR_NS}/consensus_L${MATRIX_WIDTH}_AI_${COMP_ID}.tab

	@compare-scores -null 0 						\
		-i ${COMP_DIR_NS}/consensus_L${MATRIX_WIDTH}_${REGULONS}_lnEF_classfreq.tab \
		-i ${COMP_DIR_NS}/consensus_L${MATRIX_WIDTH}_${RAND_SET}_lnEF_classfreq.tab \
		-numeric > ${COMP_DIR_NS}/consensus_L${MATRIX_WIDTH}_lnEF_${COMP_ID}.tab

	@compare-scores -null 0 						\
		-i ${COMP_DIR_NS}/consensus_L${MATRIX_WIDTH}_${REGULONS}_lnPval_classfreq.tab \
		-i ${COMP_DIR_NS}/consensus_L${MATRIX_WIDTH}_${RAND_SET}_lnPval_classfreq.tab \
		-numeric > ${COMP_DIR_NS}/consensus_L${MATRIX_WIDTH}_lnPval_${COMP_ID}.tab


################################################################
## Compare families with catalogs (GO)
CATALOG=GO
CATALOG_DIR=${RSAT}/public_html/data/genomes/${ORG}/catalogs/
CATALOG_FILE=${CATALOG_DIR}/${CATALOG}.tab
COMPA_DIR=${ORG_DIR}/catalog_comparison
COMPA_FILE=${COMPA_DIR}/${FAM}_vs_${CATALOG}
COMPARE_CMD=compare-classes -v ${V} -q ${FAM_FILE}	\
		-r ${CATALOG_FILE}			\
		-lth occ 2 -lth sig 0			\
		-sort sig				\
		-return occ,percent,proba,members	\
		-o ${COMPA_FILE}.tab			\
		-dot ${COMPA_FILE}.dot 

compare_one_catalog:
	@mkdir -p ${COMPA_DIR}
	@echo ${COMPARE_CMD}
	@${COMPARE_CMD}
	text-to-html -i ${COMPA_FILE}.tab -o ${COMPA_FILE}.html -font variable

self_compare:
	${MAKE} compare_one_catalog CATALOG=${FAM} CATALOG_FILE=${FAM_FILE}

CATALOGS=cellzome_complexes TF_targets MIPS GO ms_complexes 
compare_catalogs:
	@${MAKE} self_compare
	@for cat in ${CATALOGS}; do				\
		${MAKE} compare_one_catalog CATALOG=$${cat} ;	\
	done

