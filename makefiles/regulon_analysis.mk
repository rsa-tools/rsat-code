################################################################
## Discovery patterns in yeast regulons

include ${RSAT}/makefiles/util.mk

WHEN=now
MAKEFILE=${RSAT}/makefiles/regulon_analysis.mk

################################################################
## queue tasks on merlin
QUEUE_TASKS = oligos oligo_maps gibbs consensus MotifSampler dyads dyad_maps
queue_all:
	for family in ${REGULONS} ${RAND_FAM}; do				\
		${MAKE} multi MULTI_TASK=upstream FAM=$${family} WHEN=now;	\
	done
	for task in ${QUEUE_TASKS}; do							\
		for family in ${REGULONS} ${RAND_FAM}; do				\
			${MAKE} multi MULTI_TASK=$${task} FAM=$${family} WHEN=queue;	\
		done;									\
	done

after_queue:
	for family in ${REGULONS} ${RAND_FAM}; do					\
		${MAKE} multi MULTI_TASK=clean,synthesis,sql FAM=$${family} WHEN=now;	\
	done


################################################################
## Generate random gene selections
REGULONS=regulons_TF_aMAZE.fam
RAND_FAM=random_selection.fam
REGULONS_FILE=data/${REGULONS}
RAND_FAM_FILE=data/${RAND_FAM}
rand_fam:
	random-genes -fam ${REGULONS_FILE} -o ${RAND_FAM_FILE}	\
		 -features data/Feature_nomito.tab -org ${ORG}

################################################################
#### Parameters for pattern discovery (oligo-analysis and
#### dyad-analysis)
FAM=${REGULONS}
FAM_FILE=data/${FAMS}

MULTI_INPUT=-i ${FAM}
ORG=Saccharomyces_cerevisiae
STR=-2str
THOSIG=0
NOOV=-noov
MULTI_TASK=upstream,purge,oligos,merge_oligos,oligo_maps,dyads,dyad_maps,consensus,gibbs,MotifSampler,synthesis,sql
MULTI_BG=upstream
MULTI_EXP=-bg ${MULTI_BG}
PURGE=-purge
RES_DIR=`pwd`/results
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
MULTI_CMD=multiple-family-analysis -v ${V}					\
		${MULTI_INPUT}							\
		-thosig ${THOSIG}						\
		${PURGE}							\
		-org ${ORG}							\
		-outdir ${MULTI_DIR}						\
		${STR}								\
		-minol ${MIN_OL} -maxol ${MAX_OL}				\
		-minsp ${MIN_SP} -maxsp ${MAX_SP}				\
		-width ${MATRIX_WIDTH} -sps ${SPS}				\
		${MULTI_EXP}							\
		-sort ${SORT} -task ${MULTI_TASK}				\
		${NOOV} ${MULTI_OPT}						\
		${SKIP}							

## call multiple-family-analysis
multi:
	@echo ${MULTI_CMD}
	${MAKE} my_command MY_COMMAND="${MULTI_CMD}"


SERVER=merlin.ulb.ac.be
TO_SYNC=results
to_merlin:
	${RSYNC} ${TO_SYNC} jvanheld@${SERVER}:regulons/

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



