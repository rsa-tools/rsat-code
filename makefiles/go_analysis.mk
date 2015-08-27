################################################################
## Run GO enrichment analysis

include ${RSAT}/makefiles/util.mk
MAKEFILE=makefiles/go_analysis.mk
PYTHON=python2.7
#PYTHON=python -m cProfile -o profiling_result.txt
SCRIPT=${RSAT}/python-scripts/go_analysis.py


ORG=mycoplasma_genitalium_g37
ORGANISMS=

################################################################
## Help messages

## Generic help about GO functions
help:
	${PYTHON} ${SCRIPT} --help


## Get help about one particular task
TASK=download_go
help_task:
	${PYTHON} ${SCRIPT} ${TASK} --help

################################################################
## Get GO annotations for each gene
##

# ## Note: to get go annotation we first need to call the mgen target in
# ## order to retrieve protein to GO relations.


# install_goatools:
# 	git clone https://github.com/tanghaibao/goatools.git 
# go_all: getGOOboFile gene2GOAnnotation 
## For each GO term associated with a gene g
## We now need to declare also the relations between ancestor(t) and
## g.

################################################################
## First we need the full GO ontology (in obo format). 
GOOBO=gene_ontology_ext.obo
GOOBO_URL=http://www.geneontology.org/ontology/obo_format_1_2/${GOOBO}
GOSLIM=goslim_generic.obo
GOSLIM_URL=http://www.geneontology.org/GO_slims/${GOSLIM}
GO_DIR=results/GO
download_go:
	@echo
	@echo "Creating GO Directory	${GO_DIR}"
	@mkdir -p ${GO_DIR}
	${PYTHON} ${SCRIPT} download_go -o ${GO_DIR}/${GOOBO}

## Parse the content of the obo-formatted file
GO_REL=${GO_DIR}/GO_relations.tab
GO_DESC=${GO_DIR}/GO_description.tab
parse_go:
	@echo
	@echo "Parsing obo file	${GO_DIR}/${GOBO}"
	${PYTHON} ${SCRIPT} parse_go -f ${GO_DIR}/${GOOBO} --output_dir ${GO_DIR}
	@echo "GO term descriptions	${GO_DESC}"
	@echo "GO term relations	${GO_REL}"

## Parse the content of the obo-formatted file
get_annotations:
	@echo
	@echo "Getting gene annotations from Ensembl REST Web services"
	${PYTHON} ${SCRIPT} get_annotations -org ${ORG}

## Expand GO annotations, i.e. the gene-GO associations are
## transmitted from each class to all its ancestral classes
GO_ANNOT=annotations_table_${ORG}.tab
expand_annot:
	@echo
	@echo "Expading gene annotations from each class to its ancestor classes"
	${PYTHON} ${SCRIPT} expand -a ${GO_ANNOT} -d ${GO_DESC} -r ${GO_REL}

expand_org:
	@echo
	@echo "Expading gene annotations for organism ${ORG}"
	${PYTHON} ${SCRIPT} expand -org ${ORG}

ORG_DIR=results/ensembl_genomes/${ORG}
gene2go_one_species:
	@echo "Collecting gene-GO relationships for organism	${ORG}"
	@grep -w GO ${ORG_DIR}/proteins_xrefs.tab   | cut -f1,4 > ${ORG_DIR}/gene2GO.tab
	@echo "	${ORG_DIR}/gene2GO.tab"
	@echo
	@echo "Expanding GO annotations to ancestor classes for	${ORG}"
	@${PYTHON} ${SCRIPT} ancestor -i ${ORG_DIR}/gene2GO.tab -g ${GO_DIR}/gene_ontology_ext.obo -o ${ORG_DIR}/gene2GO_full.tab
	@echo "	${ORG_DIR}/gene2GO_full.tab"

