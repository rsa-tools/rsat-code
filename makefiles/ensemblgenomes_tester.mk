################################################################
## Tester for the python script nsemblgenome;py

include ${RSAT}/makefiles/util.mk
MAKEFILE=${RSAT}/makefiles/ensemblgenomes_tester.mk
PYTHON=python2.7
#PYTHON=python3

DATE=`date +%Y%m%d`
V=2
SCRIPT=${RSAT}/python-scripts/ensembl_genomes
NULL=NA
DIVISION=EnsemblBacteria
DIVISION_OPT=--division ${DIVISION}
TIME=time

################################################################
## Get help messages

## Default help message displayed when the program is called without any argument
help:
	${PYTHON} ${SCRIPT}

## The option -h displays a short general help, and list supported
## tasks
help_script: 
	${PYTHON} ${SCRIPT} -h

## Help about the task retrieve_species
help_species:
	${PYTHON} ${SCRIPT} retrieve_species -h

## Help about the task retrieve_features
help_features:
	${PYTHON} ${SCRIPT} retrieve_features -h


## For developers: python documentation for the ensembl access functions.
doc_python_functions:
	${PYTHON} -m pydoc Scripts/EnsemblGenomes/ensembl_functions.py

################################################################
## Collect the description of all supported species
SPECIES_LIST=species_list.txt
RES_DIR=results
GENOMES_DIR=${RES_DIR}/ensembl_genomes
ORG_DESC=${GENOMES_DIR}/species_descriptions${DIVISION}_${DATE}.tab
species_list:
	@echo "Collecting species descriptions from EnsembLGenomes"
	@mkdir -p ${GENOMES_DIR}
	time ${PYTHON} ${SCRIPT} retrieve_species -v ${V} -f ${SPECIES_LIST} -o ${ORG_DESC} --null ${NULL} ${DIVISION_OPT}
	@echo "	description table	${ORG_DESC}"
	@echo "	species_list    	${GENOMES_DIR}/${SPECIES_LIST}"

################################################################
## Get features (genes, transcripts, proteins) from EnsemblGenomes for
## one species
ORG_TASK=get_features_one_species
#ORG_TASK=gene2go_one_species

#ORG=campylobacter_jejuni_subsp_jejuni_bh_01_0142
ORG=mycoplasma_genitalium_g37

PARSING_DEPTH=5
ORG_DIR=${GENOMES_DIR}/${ORG}
get_features_one_species:
	@echo "Collecting features for species ${ORG}"
	${PYTHON} ${SCRIPT}  retrieve_features -v ${V} --outputdir ${GENOMES_DIR}  -d ${PARSING_DEPTH} -org ${ORG}
	@echo "	${ORG_DIR}"

## Retrieve features for Mycoplasma genitalium (convenient for quick
## testing, because this is a very small bacterial genome)
mgen:
	${MAKE} ${ORG_TASK} ORG=mycoplasma_genitalium_g37

## Retrieve features for Bacillus subilis
bsub:
	${MAKE} ${ORG_TASK} ORG=bacillus_subtilis_subsp_subtilis_str_168

## Retrieve features for Escherichia coli
coli:
	${MAKE} ${ORG_TASK} ORG=escherichia_coli_str_k_12_substr_mg1655


## TO BE CHECKED LATER: Retrieve features for drosophila_ananassae,
## which poses poblem
#dana:
#	${MAKE} ${ORG_TASK} ORG=drosophila_ananassae

################################################################
## Get features (genes, transcripts, proteins) for a list of selected
## organisms specified in a text file (one organism per line).
ORG_FILE=data/genoscope_tests/organism_selection_genoscope.txt
get_features_species_list:
	@echo "Collecting features for organisms specified in file ${ORG_FILE}"
	${PYTHON} ${SCRIPT} retrieve_features -v ${V} --outputdir ${GENOMES_DIR}  -d ${PARSING_DEPTH} -f ${ORG_FILE}
	@echo "	${GENOMES_DIR}"


## NOTE: GO in separate makefile go_analysis.mk
