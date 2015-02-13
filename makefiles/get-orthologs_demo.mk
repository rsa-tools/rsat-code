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

## Query organism
ORG=Escherichia_coli_K_12_substr__MG1655_uid57779

## Reference taxon
TAXON=Enterobacteriales

## Gene of interest
GENE=lexA

## Output fields
OUTPUT_FIELDS=ident,ali_len,e_value,rank

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
	@echo "	ORG		${ORG}"
	@echo "	GENE		${GENE}"
	@echo "	TAXON		${TAXON}"
	@echo "	OUTPUT_FIELDS	${OUTPUT_FIELDS}"
	@echo "	RES_DIR		${RES_DIR}"

## Get all matching genes
all_hits: res_dir
	@echo
	@echo "Collecting all hits for ${ORG} gene ${GENE} in ${TAXON}"
	get-orthologs -v 1 -org ${ORG} \
		-taxon ${TAXON} \
		-return ${OUTPUT_FIELDS} \
		-q ${GENE} ${OPT} \
		-o ${RES_DIR}/${GENE}_all_hits_${TAXON}.tab 
	@echo "	${RES_DIR}/${GENE}_all_hits_${TAXON}.tab "

## Filter hits on several criteria
filtered_hits:
	@echo
	@echo "Collecting filtered hits for ${ORG} gene ${GENE} in ${TAXON}"
	get-orthologs -v 1 -org ${ORG} \
		-taxon ${TAXON} \
		-return ${OUTPUT_FIELDS} \
		-lth ident 30 -lth ali_len 50 -uth e_value 1e-10 \
		-q ${GENE} ${OPT} \
		-o ${RES_DIR}/${GENE}_all_hits_${TAXON}.tab 
	@echo "	${RES_DIR}/${GENE}_filtered_hits_${TAXON}.tab "


## Unidirectional best hits
best_hits:
	@echo
	@echo "Getting unidirectional best hits"
	get-orthologs -v 1 -org ${ORG} \
		-taxon ${TAXON} \
		-return ${OUTPUT_FIELDS} \
		-lth ident 30 -lth ali_len 50 -uth e_value 1e-10 \
		-q ${GENE} -uth rank 1 ${OPT} \
		-o ${RES_DIR}/${GENE}_best_hits_${TAXON}.tab 
	@echo "	${RES_DIR}/${GENE}_best_hits_${TAXON}.tab "

## Bidirectional best hits (BBH)
bbh:
	@echo
	@echo "Getting bidirectional best hits (BBH)"
	get-orthologs -v 1 -org ${ORG} \
		-taxon ${TAXON} \
		-return ${OUTPUT_FIELDS} \
		-lth ident 30 -lth ali_len 50 -uth e_value 1e-10 \
		-q ${GENE} -uth rank 1 ${OPT} \
		-o ${RES_DIR}/${GENE}_best_hits_${TAXON}.tab 
	@echo "	${RES_DIR}/${GENE}_best_hits_${TAXON}.tab "

