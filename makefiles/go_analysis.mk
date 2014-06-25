################################################################
## Run GO enrichment analysis

include ${RSAT}/makefiles/util.mk
MAKEFILE=${RSAT}/makefiles/go_analysis.mk
PYTHON=python2.7
SCRIPT=${RSAT}/python-scripts/go_analysis


#ORG=mycoplasma_genitalium_g37
ORG=escherichia_coli_str_k_12_substr_mg1655
RSAT_ORG=Escherichia_coli_K_12_substr__MG1655_uid57779

################################################################
## Get GO annotations for each gene
##

## To get go annotation we first need to call the mgen target
## in order to retrieve protein to GO relations.

install_goatools:
	git clone https://github.com/tanghaibao/goatools.git 

go_all: getGOOboFile gene2GOAnnotation 

## Collect the definitions for all the terms of the GO ontology (in
## obo format).
GOOBO=gene_ontology_ext.obo
GOOBO_URL=http://www.geneontology.org/ontology/obo_format_1_2/${GOOBO}
GOSLIM=goslim_generic.obo
GOSLIM_URL=http://www.geneontology.org/GO_slims/${GOSLIM}
GO_DIR=GO_data
download_go_terms:
	@echo
	@echo "Creating GO Directory	${GO_DIR}"
	@mkdir -p ${GO_DIR}
	@echo
	@echo "Downloading GO ontology	${GOOBO_URL}"
	${PYTHON} ${SCRIPT} download_go -o ${GO_DIR}/${GOOBO}
	@echo "	${GO_DIR}/${GOOBO}"

################################################################
## Collect gene info for one species, using EnsemblGenomes REST web
## services.

################################################################
## Collect gene - GO annotations for a selected organism
ORG_DIR=${GO_DIR}/${ORG}
gene2go_one_species:
	@echo "Collecting gene-GO relationships for organism	${ORG}"
	@mkdir -p ${ORG_DIR}
	@grep -w GO results/${ORG}/proteins_xrefs.tab   | cut -f1,4 > results/${ORG}/gene2GO.tab
	@echo "	results/${ORG}/gene2GO.tab"
	@echo
	@echo "Expanding GO annotations to ancestor classes for	${ORG}"
	@${PYTON} ${SCRIPT} ancestor -i results/${ORG}/gene2GO.tab -g ${GO_DIR}/gene_ontology_ext.obo -o results/${ORG}/gene2GO_full.tab
	@echo "	results/${ORG}/gene2GO_full.tab"

################################################################
## Download reference data sets from RegulonDB

download_regulondb: download_regulondb_tf2gene download_regulondb_operons

## Download gene-factor relationships from RegulonDB
RDB_TF_GENE=network_tf_gene.txt
RDB_TF_GENE_URL=http://regulondb.ccg.unam.mx/menu/download/datasets/files/${RDB_TF_GENE}
RDB_DIR=results/regulonDB
download_regulondb_tf2gene:
	@echo
	@echo "Downloading RegulonDB TF-gene network"
	@mkdir -p ${RDB_DIR}
	@wget --no-clobber --no-host-directories --directory-prefix ${RDB_DIR} ${OPT} ${RDB_TF_GENE_URL}
	@echo "	${RDB_DIR}/${RDB_TF_GENE}"

## Download operons from RegulonDB
RDB_OPERONS=OperonSet.txt
RDB_OPERONS_URL=http://regulondb.ccg.unam.mx/menu/download/datasets/files/${RDB_OPERONS}
download_regulondb_operons:
	@echo
	@echo "Downloading RegulonDB operons"
	@mkdir -p ${RDB_DIR}
	@wget --no-clobber --no-host-directories --directory-prefix ${RDB_DIR} ${OPT} ${RDB_OPERONS_URL}
	@echo "	${RDB_DIR}/${RDB_OPERONS}"

## Generate a TF-regulon network from RegulonDB
regulondb_tf_gene_network:
	@echo
	@echo "Preparing gene - TF relationships"
	@add-gene-info -org ${ORG} -i ${RDB_DIR}/${RDB_TF_GENE} -col 
