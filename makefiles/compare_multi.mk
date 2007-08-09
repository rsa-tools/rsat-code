################################################################
## Run multiple-family-analysis on a set of regulons (or other groups
## of suposedly coregulated genes) + a set of randomly selected
## promoters (negative control).
##
## Apply then compare-scores to compare the distributions of scores
## for the motifs obtained with the positive and negative data sets,
## respectively.


include ${RSAT}/makefiles/util.mk
MAKEFILE=${RSAT}/makefiles/compare_multi.mk
V=1


################################################################
## Compare score distributions for patterns (oligo-analysis or
## dyad-analysis) Number of regulons with at least one pattern
SQL_FILE=pattern
RESTRICT_COL=4
PATTERN_TYPE=oligo
SCORE_COL=11
SCORE_NAME=sig
#ORTHO_FAM_NB=80
#RAND_NB=400
ORTHO_FAM_NB=`grep -v '^\-\-' ${ORTHO_FAM} | sort -u | wc -l`
RAND_NB=${ORTHO_FAM_NB}
COMPA_DIR=results/program_comparison/${TAXON}/taxfreq/comparisons
COMPA=${COMPA_DIR}/${PATTERN_TYPE}_${SCORE_NAME}_distrib
DISTRIB_COL_REG=6
DISTRIB_COL_RAND=7
_compare_distrib:
	@echo ${ORTHO_FAM}
	@echo ${ORTHO_FAM_NB}
	@mkdir -p ${COMPA_DIR}
	compare-score-distrib -v 1 \
		-restrict ${RESTRICT_COL} ${PATTERN_TYPE} \
		-sc ${SCORE_COL} \
		-i ${MULTI_DIR}/sql_export/${SQL_FILE}.tab -file_factor 1 ${ORTHO_FAM_NB} \
		-i ${RAND_MULTI_DIR}/sql_export/${SQL_FILE}.tab  -file_factor 2 ${RAND_NB} \
		-null -inf -null inf -null na -null nan \
		-o ${COMPA}.tab
	@echo ${COMPA}.tab

## Graph with (inverse or direct) cumulative distributions (Y) as a
## function of the ${SCORE_NAME} score (X)
	XYgraph -i ${COMPA}.tab \
		-xcol 1 -ycol ${DISTRIB_COL_REG},${DISTRIB_COL_RAND} -lines \
		-title1 '${PATTERN_TYPE} ${SCORE_NAME}' \
		-title2 'Motifs per family' \
		-xsize 800 -ysize 500 \
		-xleg1 '${SCORE_NAME}' -yleg1 'Motifs per family' \
		-legend -ymin 0 \
		${DISTRIB_GRAPH_OPT} \
		-format ${IMG_FORMAT} -o ${COMPA}_Nicum.${IMG_FORMAT}
	@echo  ${COMPA}_Nicum.${IMG_FORMAT}

## Graph with inverse cumulative distributions (Y) as a function of
## the sig score (X) Logarithmic Y axis
	XYgraph -i ${COMPA}.tab \
		-xcol 1 -ycol  ${DISTRIB_COL_REG},${DISTRIB_COL_RAND} -lines \
		-title1 '${PATTERN_TYPE} ${SCORE_NAME}' \
		-title2 'Number of motifs per family' \
		-xsize 800 -ysize 500 \
		-xleg1 '${SCORE_NAME}' -yleg1 'Motifs per family' \
		-legend \
		-ylog 2 \
		${DISTRIB_GRAPH_OPT} \
		-format ${IMG_FORMAT} -o ${COMPA}_Nicum_log.${IMG_FORMAT}
	@echo  ${COMPA}_Nicum_log.${IMG_FORMAT}

## ROC-like Graph comparing inverse cumulative score distributions between regulons (Y) and random selections (X)
	XYgraph -i ${COMPA}.tab  -same_limits \
		-xcol ${DISTRIB_COL_RAND} -ycol  ${DISTRIB_COL_REG},${DISTRIB_COL_RAND} -lines \
		-pointsize 0 \
		-xsize 500 -ysize 500 \
		-title1 '${PATTERN_TYPE} ${SCORE_NAME}' \
		-title2 'Motifs per family' \
		-xleg1 "${RAND_NB} random selections" -yleg1 "${ORTHO_FAM_NB} regulons" \
		${ROC_GRAPH_OPT} \
		-format ${IMG_FORMAT} -o ${COMPA}_Nicum_roc.${IMG_FORMAT}
	@echo ${COMPA}_Nicum_roc.${IMG_FORMAT}

RANK_COL=12
_compare_distrib_rank1:
	compare-score-distrib -v 1 \
		-restrict ${RESTRICT_COL} ${PATTERN_TYPE} \
		-restrict ${RANK_COL} 1 -sc ${SCORE_COL} \
		-i ${MULTI_DIR}/sql_export/${SQL_FILE}.tab -file_factor 1 ${ORTHO_FAM_NB} \
		-i ${RAND_MULTI_DIR}/sql_export/${SQL_FILE}.tab  -file_factor 2 ${RAND_NB} \
		-null -inf -null inf -null na -null nan\
		-o ${COMPA}_rank1.tab
	@echo ${COMPA}_rank1.tab

## Graph comparing inverse cumulative distributions (Y) as a function of the sig score (X)
	XYgraph -i ${COMPA}_rank1.tab \
		-xcol 1 -ycol  ${DISTRIB_COL_REG},${DISTRIB_COL_RAND} -lines \
		-title1 '${PATTERN_TYPE} ${SCORE_NAME}' \
		-title2 'fraction of families with at least one motif' \
		-xsize 800 -ysize 500 \
		-xleg1 '${SCORE_NAME}' -yleg1 'fraction of families' \
		-legend \
		-ygstep1 0.1 -ygstep2 0.05 -ymin 0 -ymax 1 \
		${DISTRIB_GRAPH_OPT} \
		-format ${IMG_FORMAT} -o ${COMPA}_rank1_Nicum.${IMG_FORMAT}
	@echo  ${COMPA}_rank1_Nicum.${IMG_FORMAT}

## Graph comparing inverse cumulative distributions (Y) as a function of the sig score (X)
## Logarithmic Y axis
	XYgraph -i ${COMPA}_rank1.tab \
		-xcol 1 -ycol  ${DISTRIB_COL_REG},${DISTRIB_COL_RAND} -lines \
		-title1 '${PATTERN_TYPE} ${SCORE_NAME}' \
		-title2 'fraction of families with at least one motif' \
		-xsize 800 -ysize 500 \
		-xleg1 '${SCORE_NAME}' -yleg1 'fraction of families' \
		-legend \
		-ylog 2 -ymax 1 \
		${DISTRIB_GRAPH_OPT} \
		-format ${IMG_FORMAT} -o ${COMPA}_rank1_Nicum_log.${IMG_FORMAT}
	@echo  ${COMPA}_rank1_Nicum_log.${IMG_FORMAT}


## Graph comparing inverse cumulative score distributions between regulons (Y) and random selections (X)
	XYgraph -i ${COMPA}_rank1.tab \
		-xcol ${DISTRIB_COL_RAND} -ycol  ${DISTRIB_COL_REG},${DISTRIB_COL_RAND} -lines \
		-pointsize 0 \
		-xsize 500 -ysize 500 \
		-title1 '${PATTERN_TYPE} ${SCORE_NAME}' \
		-title2 'fraction of families with at least one motif' \
		-xleg1 "${RAND_NB} random selections" -yleg1 "${ORTHO_FAM_NB} ortholog groups" \
		-gstep1 0.1 -gstep2 0.05 \
		-min 0 -max 1 \
		${ROC_GRAPH_OPT} \
		-format ${IMG_FORMAT} -o ${COMPA}_rank1_Nicum_roc.${IMG_FORMAT}
	@echo ${COMPA}_rank1_Nicum_roc.${IMG_FORMAT}

compare_oligos:
	@${MAKE} _compare_distrib PATTERN_TYPE=oligo
	@${MAKE} _compare_distrib_rank1 PATTERN_TYPE=oligo

compare_dyads:
	@${MAKE} _compare_distrib PATTERN_TYPE=dyad
	@${MAKE} _compare_distrib_rank1 PATTERN_TYPE=dyad

################################################################
## Compare score distributions of MotifSampler
_compare_MS_one_score:
	${MAKE} _compare_distrib SQL_FILE=matrix \
		RESTRICT_COL=4 PATTERN_TYPE=MotifSampler 
	${MAKE} _compare_distrib_rank1 SQL_FILE=matrix \
		RESTRICT_COL=4 PATTERN_TYPE=MotifSampler \
		RANK_COL=6

compare_MS_IC:
	${MAKE} _compare_MS_one_score \
		SCORE_COL=36 \
		SCORE_NAME=InfoContent
compare_MS_LL:
	${MAKE} _compare_MS_one_score \
		SCORE_COL=37 \
		SCORE_NAME=LogLikelihood
compare_MS_CS:
	${MAKE} _compare_MS_one_score \
		SCORE_COL=38 \
		SCORE_NAME=ConsensusScore

compare_MS_info:
	${MAKE} _compare_MS_one_score \
		SCORE_COL=17 SCORE_NAME=total.information

compare_MS_infocol:
	${MAKE} _compare_MS_one_score \
		SCORE_COL=19 \
		SCORE_NAME=information.per.column

compare_MS: compare_MS_IC compare_MS_LL compare_MS_CS compare_MS_info compare_MS_infocol

################################################################
## Compare score distributions of MEME
_compare_meme_one_score:
	${MAKE} _compare_distrib SQL_FILE=matrix \
		RESTRICT_COL=4 PATTERN_TYPE=meme
	${MAKE} _compare_distrib_rank1 SQL_FILE=matrix \
		RESTRICT_COL=4 PATTERN_TYPE=meme \
		RANK_COL=6

compare_meme_llr:
	${MAKE} _compare_meme_one_score \
		SCORE_COL=34 \
		SCORE_NAME=llr \
		DISTRIB_COL_REG=4 DISTRIB_COL_RAND=5

compare_meme_Evalue:
	${MAKE} _compare_meme_one_score \
		SCORE_COL=35 \
		SCORE_NAME=E-value \
		DISTRIB_GRAPH_OPT=-xlog \
		DISTRIB_COL_REG=4 DISTRIB_COL_RAND=5

compare_meme_info:
	${MAKE} _compare_meme_one_score \
		SCORE_COL=17 SCORE_NAME=total.information

compare_meme_infocol:
	${MAKE} _compare_meme_one_score \
		SCORE_COL=19 \
		SCORE_NAME=information.per.column

compare_meme: compare_meme_Evalue compare_meme_llr compare_meme_info compare_meme_infocol


################################################################
## Compare score distributions of consensus
_compare_consensus_one_score:
	${MAKE} _compare_distrib SQL_FILE=matrix \
		RESTRICT_COL=4 PATTERN_TYPE=consensus
	${MAKE} _compare_distrib_rank1 SQL_FILE=matrix \
		RESTRICT_COL=4 PATTERN_TYPE=consensus \
		RANK_COL=6

compare_consensus_lnEvalue:
	${MAKE} _compare_consensus_one_score \
		SCORE_COL=26 SCORE_NAME=lnE-value \
		DISTRIB_COL_REG=4 DISTRIB_COL_RAND=5

compare_consensus_lnPvalue:
	${MAKE} _compare_consensus_one_score \
		SCORE_COL=24 SCORE_NAME=lnP-value \
		DISTRIB_COL_REG=4 DISTRIB_COL_RAND=5

compare_consensus_UI:
	${MAKE} _compare_consensus_one_score \
		SCORE_COL=28 SCORE_NAME=Unadjusted_info 

compare_consensus_AI:
	${MAKE} _compare_consensus_one_score \
		SCORE_COL=27 SCORE_NAME=Adjusted_info  

compare_consensus_info:
	${MAKE} _compare_consensus_one_score \
		SCORE_COL=17 SCORE_NAME=total.information

compare_consensus_infocol:
	${MAKE} _compare_consensus_one_score \
		SCORE_COL=19 SCORE_NAME=information.per.column

compare_consensus:  compare_consensus_lnPvalue compare_consensus_UI compare_consensus_AI compare_consensus_lnEvalue compare_consensus_info compare_consensus_infocol 

################################################################
## Compare score distributions of gibbs
_compare_gibbs_one_score:
	${MAKE} _compare_distrib SQL_FILE=matrix \
		RESTRICT_COL=4 PATTERN_TYPE=gibbs 
	${MAKE} _compare_distrib_rank1 SQL_FILE=matrix \
		RESTRICT_COL=4 PATTERN_TYPE=gibbs \
		RANK_COL=6 

compare_gibbs_MAP:
	${MAKE} _compare_gibbs_one_score \
		SCORE_COL=29 SCORE_NAME=MAP

compare_gibbs_MAP_per_site:
	${MAKE} _compare_gibbs_one_score \
		SCORE_COL=30 SCORE_NAME=MAP.per.site

compare_gibbs_betaprior_map:
	${MAKE} _compare_gibbs_one_score \
		SCORE_COL=31 SCORE_NAME=betaprior.MAP

compare_gibbs_model_map:
	${MAKE} _compare_gibbs_one_score \
		SCORE_COL=32 SCORE_NAME=model.MAP

compare_gibbs_info:
	${MAKE} _compare_gibbs_one_score \
		SCORE_COL=17 SCORE_NAME=total.information

compare_gibbs_infocol:
	${MAKE} _compare_gibbs_one_score \
		SCORE_COL=19 SCORE_NAME=information.per.column

compare_gibbs: compare_gibbs_MAP compare_gibbs_MAP_per_site compare_gibbs_info compare_gibbs_infocol

problems_gibbs: compare_gibbs_betaprior_map compare_gibbs_model_map

################################################################
## Compare score distributions of AlignACE
_compare_AlignACE_one_score:
	${MAKE} _compare_distrib SQL_FILE=matrix \
		RESTRICT_COL=4 PATTERN_TYPE=AlignACE 
	${MAKE} _compare_distrib_rank1 SQL_FILE=matrix \
		RESTRICT_COL=4 PATTERN_TYPE=AlignACE \
		RANK_COL=6 

compare_AlignACE_MAP:
	${MAKE} _compare_AlignACE_one_score \
		SCORE_COL=29 SCORE_NAME=MAP

compare_AlignACE_MAP_per_site:
	${MAKE} _compare_AlignACE_one_score \
		SCORE_COL=30 SCORE_NAME=MAP.per.site

compare_AlignACE_info:
	${MAKE} _compare_AlignACE_one_score \
		SCORE_COL=17 SCORE_NAME=total.information

compare_AlignACE_infocol:
	${MAKE} _compare_AlignACE_one_score \
		SCORE_COL=19 SCORE_NAME=information.per.column

compare_AlignACE: compare_AlignACE_MAP compare_AlignACE_MAP_per_site compare_AlignACE_info compare_AlignACE_infocol
