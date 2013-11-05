## ##############################################################
## Analyze conserved motifs in upstream sequences of archaeal
## sequences

include ${RSAT}/makefiles/util.mk
MAKEFILE=makefiles/ortho_disco_validation.mk

GENE=LEXA
REF_ORG=Escherichia_coli_K12
TAXON=Gammaproteobacteria
MAIN_DIR=${HOME}/research/ortho_disco/regulondb_validation
ALL_GENES=${REGULONDB}

test:
	@echo "${ALL_GENES}"


################################################################
## Validation of the method : 
## - Discover patterns in upstream regions of all genes from RegulonDB
## - create a graph with genes as nodes and arcs representing co-dicovered
##   patterns
## - compare the discovered graph with the co-factor graph off regulonDB
##   (to be done)

## Select target genes from RegulonDB
REGULONDB_TABLE=data/gene_factor_Escherichia_coli_K12
REGULONDB_GENES=data/regulondb_genes.tab
regulondb_genes:
	@cut -f 1 ${REGULONDB_TABLE}.tab | sort -u > ${REGULONDB_GENES}
	@echo "${REGULONDB_GENES}"

## Select transcription factors from RegulonDB
REGULONDB_FACTORS=data/regulondb_factors.tab
regulondb_factors:
	@cut -f 2 ${REGULONDB_TABLE}.tab | sort -u  > ${REGULONDB_FACTORS}
	@echo "${REGULONDB_FACTORS}"

## Select genes and transcription factors from RegulonDB, both will be used
## for pattern discovery
REGULONDB_GF=data/regulondb_genes_and_factors.tab
regulondb_gf: regulondb_genes regulondb_factors
	cat ${REGULONDB_GENES} ${REGULONDB_FACTORS} | tr 'a-z' 'A-Z'  | sort -u > ${REGULONDB_GF}
	@echo "${REGULONDB_GF}"

## Run all the tasks for regulonDB analysis
REGULONDB=`cat ${REGULONDB_GF}| grep -v '^;' | xargs`
REGULON_TAXON=
regulondb_analysis:
	@echo "Analyzing genes from RegulonDB"
	@${MAKE} all_tasks_all_genes REF_ORG=Escherichia_coli_K12 TAXON=${REGULON_TAXON}  ALL_GENES="${REGULONDB}"
	@${MAKE} index_regulondb

## Index the results obtained with RegulonDB
index_regulondb:
	@${MAKE} index_results REF_ORG=Escherichia_coli_K12 TAXON=${REGULON_TAXON} ALL_GENES="${REGULONDB}"



## ##############################################################
## A small comparison for didactic purposes only
TEST_GENES=${COMPA_DIR}/test_genes.tab
TEST_DYAD_FILES=${COMPA_DIR}/test_dyad_files
test_genes:
	grep -v ';' data/regulondb_genes_and_factors.tab | grep -E "(PUR)|(PUT)|(MET)" > ${TEST_GENES}
	@echo "Test genes	${TEST_GENES}"
	grep -f ${TEST_GENES} ${DYAD_FILE_LIST} > ${TEST_DYAD_FILES}
	@echo "test dyad files	${TEST_DYAD_FILES}"
	${MAKE} dyad_profiles DYAD_FILE_LIST=${TEST_DYAD_FILES} COMPA_TABLE=${COMPA_DIR}/test_${REF_ORG}_${TAXON}_dyad_profiles.tab
	${MAKE} dyad_classes DYAD_FILE_LIST=${TEST_DYAD_FILES} COMPA_CLASSES=${COMPA_DIR}/test_${REF_ORG}_${TAXON}_dyad_classes.tab
	@grep '#' ${TEST_GENES} ${GENE_PAIRS}.tab > ${GENE_PAIRS}_test.tab
	@grep -f ${TEST_GENES} ${GENE_PAIRS}.tab | sort  >> ${GENE_PAIRS}_test.tab
	@echo ${GENE_PAIRS}_test.tab
	grep -f ${TEST_GENES} ${GENE_PAIRS}.tab | convert-graph -from tab -to dot -scol 1 -tcol 2 -wcol 8 -o ${GENE_PAIRS}_test.dot
	@echo ${GENE_PAIRS}_test.dot

## ##############################################################
## Compare the clustering results (clusters from Stephane Robin) with
## the regulons annotated in RegulonDB
REGULONDB_TABLE_UC=data/gene_factor_Escherichia_coli_K12_uppercase.tab 
CLUSTERS=${COMPA_DIR}/Erdos2_Param_Ward_Q12_clusters
clusters_vs_regulons:
	cat ${REGULONDB_TABLE} | tr a-z A-Z > ${REGULONDB_TABLE_UC}
	compare-classes -v 1 -r ${REGULONDB_TABLE_UC} -q ${CLUSTERS}.tab \
		-o ${CLUSTERS}_vs_regulondb_xtab.tab -matrix RandQ
	compare-classes -v 1 -r ${REGULONDB_TABLE_UC} -q ${CLUSTERS}.tab \
		-return jac_sim,members,occ,percent,proba,rank -sort sig \
		-o ${CLUSTERS}_vs_regulondb_stats.tab


FACTOR_GENE=data/factor_gene_Escherichia_coli_K12_regulonDB
GENE_GENE=data/gene_gene_Escherichia_coli_K12_regulonDB
regulondb_gene_pairs:
	awk -F '\t' '{print $$2"\t"$$1}' ${REGULONDB_TABLE_UC} | grep -v ";" | tr a-z A-Z > ${FACTOR_GENE}.tab
	@echo ${FACTOR_GENE}.tab
	compare-classes -v ${V} -i ${FACTOR_GENE}.tab -o ${GENE_GENE}.tab -return occ -lth QR 1 -sort names -distinct -triangle
	@echo ${GENE_GENE}.tab

	convert-graph -i  ${GENE_GENE}.tab -o  ${GENE_GENE}.gml -from tab -to gml
	@echo ${GENE_GENE}.gml
	convert-graph -i  ${GENE_GENE}.tab -o  ${GENE_GENE}.dot -from tab -to dot
	@echo ${GENE_GENE}.dot

compare_regulondb_gene_pairs:
	grep -v "^;" ${GENE_GENE}.tab | grep -v "^#" | awk '{print $$1"_"$$2"\tregulonDB"}'  >  ${GENE_GENE}_sets.tab
	grep -v "^;" ${GENE_PAIRS}.tab | grep -v "^#" | awk '{print $$1"_"$$2"\tpredicted"}' >>  ${GENE_GENE}_sets.tab
	compare-classes -i ${GENE_GENE}_sets.tab \
		 -v ${V} -return occ,proba \
		 -o ${GENE_GENE}_sets_compa.tab
	text-to-html -i ${GENE_GENE}_sets_compa.tab -o ${GENE_GENE}_sets_compa.html 

################################################################
## Draw a ggraph with RegulonDB data
regulondb_graph:
	convert-graph  -from tab -to gml -i ${FACTOR_GENE}.tab -o ${FACTOR_GENE}.gml
	@echo  ${FACTOR_GENE}.gml
	convert-graph  -from tab -to dot -i ${FACTOR_GENE}.tab -o ${FACTOR_GENE}.dot
	@echo  ${FACTOR_GENE}.dot


################################################################
## Synchronize results from BIGRE to my machine
from_bigre:
	(cd ..; rsync -ruptvlz -e ssh --exclude jobs --exclude *.fasta ${OPT} \
	jvanheld@merlin.bigre.ulb.ac.be:research/ortho_disco/regulondb_validation .)

to_bigre:
	(cd ..; rsync -ruptvlz -e ssh ${OPT} regulondb_validation jvanheld@merlin.bigre.ulb.ac.be:research/ortho_disco/)

