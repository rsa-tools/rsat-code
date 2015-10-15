################################################################
## Run GO enrichment analysis

include ${RSAT}/makefiles/util.mk
MAKEFILE=${RSAT}/makefiles/go_analysis.mk
PYTHON=python2.7
#PYTHON=python -m cProfile -o profiling_result.txt
SCRIPT=${RSAT}/python-scripts/go_analysis.py


#ORG=mycoplasma_genitalium_g37
ORG=bacillus_subtilis_subsp_subtilis_str_168
ORGANISMS=escherichia_coli_str_k_12_substr_mg1655 \
	saccharomyces_cerevisiae \
	drosophila_melanogaster \
	caenorhabditis_elegans \
	arabidopsis_thaliana \
	bacillus_subtilis_subsp_subtilis_str_168 \
	mycoplasma_genitalium_g37 \
	pseudomonas_aeruginosa_pao1_ve13

ORG_UCFIRST= `perl -e 'print ucfirst(${ORG});'`
ORG_DIR=${RSAT}/public_html/data/genomes/${ORG_UCFIRST}
GO_INSTALL_DIR=${ORG_DIR}/go_annotations

################################################################
## Help messages

## Generic help about GO functions
help:
	${PYTHON} ${SCRIPT} --help


list_param:
	@echo
	@echo "go_analysis.mk"
	@echo "	ORGANISMS		${ORGANISMS}"
	@echo "	ORG			${ORG}"
	@echo "	ORG_UCFIRST		${ORG_UCFIRST}"
	@echo "	GO_INSTALL_DIR		${GO_INSTALL_DIR}"
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

################################################################
## Collect organism-specific annotations

## Get GO annotations fo one specific organism and export the result
## in a tab-delimited file.
GO_DIR=results/GO_annotations/${ORG_UCFIRST}
get_annotations:
	@echo
	@echo "Getting GO annotations"
	@echo "	ORG	${ORG}"
	@echo "	GO_DIR	${GO_DIR}"
	@mkdir -p ${GO_DIR}
	${PYTHON} ${SCRIPT} get_annotations -org ${ORG} --output_dir ${GO_DIR}

## Expand GO annotations, i.e. the gene-GO associations are
## transmitted from each class to all its ancestral classes
GO_ANNOT=GO_annotations_${ORG}.tab
expand_annot:
	@echo
	@echo "Go annotation directory	${GO_DIR}"
	@mkdir -p ${GO_DIR}
	@echo "Expading gene annotations from each class to its ancestor classes"
	${PYTHON} ${SCRIPT} expand -a ${GO_DIR}/${GO_ANNOT} -d ${GO_DESC} -r ${GO_REL} --output_dir ${GO_DIR}

expand_org:
	@echo
	@echo "Expading gene annotations for organism ${ORG}"
	${PYTHON} ${SCRIPT} expand -org ${ORG}  --output_dir ${GO_DIR}

install_annot:
	@echo
	@echo "GO_INSTALL_DIR	${GO_INSTALL_DIR}"
	@mkdir -p ${GO_INSTALL_DIR}
	${MAKE} OPT='--output_dir ${GO_INSTALL_DIR}' get_annotations expand_annot

install_yeast:
	${MAKE} install_annot ORG=saccharomyces_cerevisiae

################################################################
## Collect info from cross-reference files (obsolete method)
ORG_DIR=results/ensembl_genomes/${ORG}
gene2go_one_species:
	@echo "Collecting gene-GO relationships for organism	${ORG}"
	@grep -w GO ${ORG_DIR}/proteins_xrefs.tab   | cut -f1,4 > ${ORG_DIR}/gene2GO.tab
	@echo "	${ORG_DIR}/gene2GO.tab"
	@echo
	@echo "Expanding GO annotations to ancestor classes for	${ORG}"
	@${PYTHON} ${SCRIPT} ancestor -i ${ORG_DIR}/gene2GO.tab -g ${GO_DIR}/gene_ontology_ext.obo -o ${ORG_DIR}/gene2GO_full.tab
	@echo "	${ORG_DIR}/gene2GO_full.tab"

################################################################
## Iterate task over all selected organisms
ORG_TASK=get_annotations expand_annot
all_organisms:
	@echo
	@echo "Iterating task over selected organisms"
	@echo "	${ORGANISMS}" 
	@for org in ${ORGANISMS}; do \
		${MAKE} ORG=$${org} ${ORG_TASK}; \
	done

