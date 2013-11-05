## CVS: added a target bg_model_meme

################################################################
## Evaluation of pattern discovery methods

MAKEFILE=${RSAT}/makefiles/pattern_disco_eval.mk

include ${RSAT}/makefiles/util.mk

ORG=Saccharomyces_cerevisiae
SITES=data/${ORG}/transfac_sites_${ORG}
IMG_FORMAT=jpg
site_map:
	feature-map -i ${SITES}.tab -o ${SITES}.${IMG_FORMAT} -format	\
		${IMG_FORMAT} -legend -scalebar -scalestep 50 -mspacing 2	\
		-label id -from -1000 -to 50 -mlen 1000

################################################################
## Generate MEME models
NOORF=-noorf
NOOV=-ovlp
MEME_BG_DIR=bg_models/meme/${ORG}
MEME_BG_FILE=`pwd`/${MEME_BG_DIR}/${ORG}_upstream-noorf${NOOV}-1str_6.freq
bg_model_meme:
	@echo "Generating background model for MEME	${ORG}	${MEME_BG_FILE}"
	@mkdir -p ${MEME_BG_DIR}
	@echo "# MEME background model ${ORG} upstream-noorf" > ${MEME_BG_FILE}
	for ol in 1 2 3 4 5 6; do ${MAKE} bg_model_meme_one_size OL=$${ol}; done
	@echo "Generated background model for MEME	${ORG}	${MEME_BG_FILE}"

OL=6
RSAT_BG_FILE=${RSAT}/data/genomes/${ORG}/oligo-frequencies/${OL}nt_upstream-noorf_${ORG}${NOOV}-1str.freq.gz
bg_model_meme_one_size:
	@echo ${RSAT_BG_FILE}
	@echo "# seq	frequency ${OL}nt_upstream-noorf${NOOV}-1str" >> ${MEME_BG_FILE}
	@gunzip -c ${RSAT_BG_FILE} | grep -v '^;' | cut -f 1,3 >> ${MEME_BG_FILE}

################################################################
## Generate MotifSampler model
NOORF=-noorf
NOOV=-ovlp
MOTIFSAMPLER_BG_DIR=bg_models/MotifSampler/${ORG}
MOTIFSAMPLER_BG_FILE=`pwd`/${MOTIFSAMPLER_BG_DIR}/${ORG}_upstream-noorf${NOOV}-1str_6.freq
ORDER=6
bg_model_MotifSampler:
	@echo "Generating background model for MotifSampler     ${ORG}  ${MOTIFSAMPLER_BG_FILE}"
	@mkdir -p ${MOTIFSAMPLER_BG_DIR}
	@echo "# MotifSampler background model ${ORG} upstream-noorf" > ${MOTIFSAMPLER_BG_FILE}
	${MAKE} CreateBackgroundModel -f -b ${MOTIFSAMPLER_BG_FILE} -o ${ORDER} -n ${ORG}
	@echo "Generated background model for MotifSampler      ${ORG}  ${MOTIFSAMPLER_BG_FILE}"

################################################################
## Run pattern discovery programs with the regulons
REGULON_FILE=gene_factor_${ORG}.tab
REGULONS=data/${ORG}/${REGULON_FILE}
REGULON_DIR=results/${ORG}/regulons/
MULTI_TASK=upstream,purge,oligos,oligo_maps,dyads,dyad_maps,synthesis,sql,clean,consensus,meme,gibbs
FAM=${REGULONS}
MULTI_DIR=${REGULON_DIR}
multi:
	@mkdir -p ${MULTI_DIR}
	multiple-family-analysis -v 1 -i ${FAM} -org ${ORG} -task	\
		${MULTI_TASK} -2str -sort score -outdir ${MULTI_DIR}	\
		${OPT} -user multifam -password multifam -login		\
		-MEME_bfile ${MEME_BG_FILE} \
		multifam -schema multifam ${OPT}

################################################################
## Load the results of multiple-family-aalysis into the database
SQL_DIR=${MULTI_DIR}/sql_export/sql_scripts/mysql
load:
	(cd ${SQL_DIR}; make load)

################################################################
## Select random genes for the negative control
RAND_DIR=results/${ORG}/random_genes
RAND_GENES=${RAND_DIR}/rand.fam
RAND_REPET=1
rand_genes:
	@mkdir -p ${RAND_DIR}
	random-genes -org ${ORG} -fam ${REGULONS} -o ${RAND_GENES} -r ${RAND_REPET}


################################################################
## Run pattern discovery programs with the random gene selections 
rand_multi:
	${MAKE} multi MULTI_DIR=${RAND_DIR} FAM=${RAND_GENES}

################################################################
## synchronize on mamaze
to_bigre:
	rsync -ruptvl -e ssh . jvanheld@mamaze.bigre.ulb.ac.be:pattern_disco_eval

## ##############################################################
## Score distributions

DATA_TYPE=random_genes
PROGRAM=meme
SCORE_COLUMN=2
SQL_FILE=results/${ORG}/${DATA_TYPE}/sql_export/matrix.tab
meme_score_distrib:
	${MAKE} score_distrib PROGRAM=meme SCORE_COLUMN=2 DATA_TYPE=random_genes
	${MAKE} score_distrib PROGRAM=meme SCORE_COLUMN=2 DATA_TYPE=regulons

consensus_score_distrib:
	${MAKE} score_distrib PROGRAM=consensus SCORE_COLUMN=13 DATA_TYPE=random_genes
	${MAKE} score_distrib PROGRAM=consensus SCORE_COLUMN=13 DATA_TYPE=regulons

## Calculate distibution of E-values obtained with MEME
SCORE_DIR=results/${ORG}/scores
SCORE_ENUM=${SCORE_DIR}/${ORG}_${PROGRAM}_${DATA_TYPE}_scores
SCORE_DISTRIB=${SCORE_ENUM}_distrib
ENUM_CMD=grep -v '^;' ${SQL_FILE} \
	| grep -v '^--' \
	| grep ${PROGRAM}  | cut -f 1,${SCORE_COLUMN} \
	| awk -F '\t' '{print $$1"\t"$$2"\t"(-log($$2)/log(10))}' >> ${SCORE_ENUM}.txt
score_distrib:
	@mkdir -p ${SCORE_DIR}
	@echo "; dataset	score	-log10(score)" >  ${SCORE_ENUM}.txt
	@echo ${ENUM_CMD}
	${ENUM_CMD}
	@echo "Scores	${PROGRAM}	${DATA_TYPE}	${SCORE_ENUM}.txt"
	@${MAKE} score_distrib_graph

## Draw frequency polygons with score distributions
TITLE= '${PROGRAM} scores with ${DATA_TYPE}'
score_distrib_graph:
	cut -f 3 ${SCORE_ENUM}.txt | classfreq  -v 1 -ci 1 -o ${SCORE_DISTRIB}.tab
	XYgraph -xcol 3 -ycol 4,5,6 -lines -xsize 800 -ysize 400 \
		-title1 ${TITLE} \
		-yleg1 'number of motifs' \
		-xleg1 'sig=-log10(E-value)' \
		-i ${SCORE_DISTRIB}.tab \
		-o ${SCORE_DISTRIB}.jpg
	@echo "Distrib	${PROGRAM}	${DATA_TYPE}	${SCORE_DISTRIB}.tab"
	@echo "Graph	${PROGRAM}	${DATA_TYPE}	${SCORE_DISTRIB}.jpg"
