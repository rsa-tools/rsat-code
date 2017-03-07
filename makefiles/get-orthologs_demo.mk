################################################################
## Demonstrator for get-orthologs
##
## We show here the different ways to tune the options of
## get-orthologs in order to obtain single blast hits, best hits or
## bidirectional best hits (BBH).

include ${RSAT}/makefiles/util.mk
MAKEFILE=${RSAT}/makefiles/get-orthologs_demo.mk

################################################################
## Specify the default parameters

## Gene of interest
GENE=lexA

## Query organism
ORG=Escherichia_coli_K_12_substr__MG1655_uid57779

## Reference taxon
TAXON=Enterobacteriaceae

## Taxonomic depth
DEPTH=5

## Output fields
OUTPUT_FIELDS=ident,ali_len,e_value,rank,s_rank

## Additional options can be passed to get-orthologs (by default we
## leave them empty)
OPT=

## Result directory
RES_DIR=results/get-orthologs_demo

################################################################
## Definition of the targets

## Create the result directory
res_dir:
	@mkdir -p ${RES_DIR}


## Lis parameters
list_param:
	@echo "	GENE		${GENE}"
	@echo "	ORG		${ORG}"
	@echo "	TAXON		${TAXON}"
	@echo "	DEPTH		${DEPTH}"
	@echo "	OUTPUT_FIELDS	${OUTPUT_FIELDS}"
	@echo "	RES_DIR		${RES_DIR}"


################################################################
## Genome blast
## THIS REQUIRES WRITE PERMISSIONS ON THE RSAT DIRECTORY
genome_blast_formatdb:
	@echo
	@echo "genome-blast -task formatdb	${ORG}	${TAXON}"
	genome-blast -v ${V} -q ${ORG} -dbtaxon ${TAXON} -reciprocal -task formatdb ${OPT}

genome_blast_blastall:
	@echo
	@echo "genome-blast -task formatdb	${ORG}	${TAXON}"
	genome-blast -v ${V} -q ${ORG} -dbtaxon ${TAXON} -reciprocal -task blastall ${OPT}

## Count the number of hits in a given file
HIT_FILE=${RES_DIR}/${GENE}_${TAXON}_depth${DEPTH}_BBH.tab
HITS=`grep -v '^;' ${HIT_FILE} | grep -v '^\#' | wc -l`
count_hits:
	@echo "	${HITS} hits	${RES_DIR}/${GENE}_${TAXON}_all_hits.tab"

all_tests: list_param all_hits filtered_hits best_hits bbh_manual bbh_auto bbh_depth


## Get all matching genes
all_hits: res_dir
	@echo
	@echo "Collecting all hits for ${ORG} gene ${GENE} in ${TAXON}"
	get-orthologs -v ${V} -org ${ORG} \
		-taxon ${TAXON} \
		-return ${OUTPUT_FIELDS} \
		-q ${GENE} ${OPT} \
		-o ${RES_DIR}/${GENE}_${TAXON}_all_hits.tab
	@${MAKE} count_hits HIT_FILE=${RES_DIR}/${GENE}_${TAXON}_all_hits.tab

## Filter hits on several criteria
filtered_hits: res_dir
	@echo
	@echo "Collecting filtered hits for ${ORG} gene ${GENE} in ${TAXON}"
	get-orthologs -v ${V} -org ${ORG} \
		-taxon ${TAXON} \
		-return ${OUTPUT_FIELDS} \
		-lth ident 30 -lth ali_len 50 -uth e_value 1e-5 \
		-q ${GENE} ${OPT} \
		-o ${RES_DIR}/${GENE}_${TAXON}_filtered_hits.tab
	@${MAKE} count_hits HIT_FILE=${RES_DIR}/${GENE}_${TAXON}_filtered_hits.tab


## Unidirectional best hits
best_hits: res_dir
	@echo
	@echo "Collecting unidirectional best hits (BH) for ${ORG} gene ${GENE} in ${TAXON}"
	get-orthologs -v ${V} -org ${ORG} \
		-taxon ${TAXON} \
		-return ${OUTPUT_FIELDS} \
		-lth ident 30 -lth ali_len 50 -uth e_value 1e-5 \
		-q ${GENE} -uth rank 1 ${OPT} \
		-o ${RES_DIR}/${GENE}_${TAXON}_unidirBH.tab
	@${MAKE} count_hits HIT_FILE=${RES_DIR}/${GENE}_${TAXON}_unidirBH.tab

## Bidirectional best hits (BBH)
bbh_manual: res_dir
	@echo
	@echo "Collecting bidirectional best hits (BBH) for ${ORG} gene ${GENE} in ${TAXON}"
	get-orthologs -v ${V} -org ${ORG} \
		-taxon ${TAXON} \
		-return ${OUTPUT_FIELDS} \
		-lth ident 30 -lth ali_len 50 -uth e_value 1e-5 \
		-q ${GENE} -uth rank 1 -uth s_rank ${OPT} \
		-o ${RES_DIR}/${GENE}_${TAXON}_BBH_manual.tab 
	@${MAKE} count_hits HIT_FILE=${RES_DIR}/${GENE}_${TAXON}_BBH_manual.tab

## Bidirectional best hits (BBH) in automatic mode
bbh_auto: res_dir
	@echo
	@echo "Collecting bidirectional best hits (BBH) for ${ORG} gene ${GENE} in ${TAXON}"
	get-orthologs -v ${V} -org ${ORG} \
		-taxon ${TAXON} \
		-return ${OUTPUT_FIELDS} \
		-type BBH \
		-q ${GENE} ${OPT} \
		-o ${RES_DIR}/${GENE}_${TAXON}_BBH_auto.tab 
	@${MAKE} count_hits HIT_FILE=${RES_DIR}/${GENE}_${TAXON}_BBH.tab

## Bidirectional best hits (BBH)
bbh_depth: res_dir
	@echo
	@echo "Getting bidirectional best hits (BBH) at depth ${DEPTH}"
	get-orthologs -v ${V} -org ${ORG} \
		-taxon ${TAXON} -depth ${DEPTH} \
		-return ${OUTPUT_FIELDS} \
		-lth ident 30 -lth ali_len 50 -uth e_value 1e-5 \
		-q ${GENE} -uth rank 1 ${OPT}  \
		-o ${RES_DIR}/${GENE}_${TAXON}_depth${DEPTH}_BBH.tab
	@${MAKE} count_hits HIT_FILE=${RES_DIR}/${GENE}_${TAXON}_depth${DEPTH}_BBH.tab

