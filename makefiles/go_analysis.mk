################################################################
## Run GO enrichment analysis

include ${RSAT}/makefiles/util.mk
MAKEFILE=${RSAT}/makefiles/go_analysis.mk
PYTHON=python2.7
SCRIPT=${RSAT}/python-scripts/go_analysis


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
