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
## Run pattern discovery programs with the regulons
REGULONS=data/${ORG}/gene_factor_${ORG}.tab
REGULON_DIR=results/${ORG}/regulons/
MULTI_TASK=upstream,purge,oligos,oligo_maps,dyads,dyad_maps,synthesis,sql,clean,consensus,meme,gibbs
FAM=${REGULONS}
MULTI_DIR=${REGULON_DIR}
multi:
	@mkdir -p ${MULTI_DIR}
	multiple-family-analysis -v 1 -i ${FAM} -org ${ORG} -task	\
		${MULTI_TASK} -2str -sort score -outdir ${MULTI_DIR}	\
		${OPT} -user multifam -password multifam -login		\
		multifam -schema multifam

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
to_scmbb:
	rsync -ruptvl -e ssh . jvanheld@mamaze.scmbb.ulb.ac.be:pattern_disco_eval