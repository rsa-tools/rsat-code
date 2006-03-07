## ##############################################################
## Discover conserved motifs in upstream sequences of orthologs for each gene
## of a genome, and infer a graph of co-regulation by comparing discovered
## motifs between each pair of genes

## THE FOLLOWING VARIABLES HAVE TO BE DEFINED
## GENE
## REF_ORG
## TAXON
## MAIN_DIR
##

## Default parameters
GENE=lexA
REF_ORG=Escherichia_coli_K12
TAXON=Gammaproteobacteria
MAIN_DIR=.

include ${RSAT}/makefiles/util.mk
MAKEFILE=${RSAT}/makefiles/ortho_disco.mk
MAKE=make -s -f ${MAKEFILE} REF_ORG=${REF_ORG} TAXON=${TAXON} GENE=${GENE} MAIN_DIR=${MAIN_DIR}

## Verbosity level
V=1

## List the parameters
list_parameters:
	@echo ""
	@echo "REF_ORG      ${REF_ORG}"
	@echo "TAXON        ${TAXON}"
	@echo "MAIN_DIR     ${MAIN_DIR}"
	@echo "RESULT_DIR   ${RESULT_DIR}"

## List of all genes to be analyzed. 
## By default, all the genes of the considered organism (REF_ORG)
ALL_GENES=` grep -v '^--' ${RSAT}/data/genomes/${REF_ORG}/genome/feature.tab | awk -F '\t' '$$2=="CDS"'| cut -f 3 | sort -u | xargs`
list_all_genes:
	@echo ""
	@echo "All genes"
	@echo "${ALL_GENES}"


################################################################
## Run all the tasks on all the genes
all_tasks_all_genes:
	@echo "All tasks applied on all genes	${REF_ORG}	${TAXON}"
	@${MAKE} list_all_genes 
	@${MAKE} list_parameters
	@for g in ${ALL_GENES} ; do \
		${MAKE} all_tasks REF_ORG=${REF_ORG} TAXON=${TAXON} MAIN_DIR=${MAIN_DIR} GENE=$${g} ; \
	done

################################################################
## Run all the tasks for a single gene
ALL_TASKS_CMD=${ORTHO_CMD} ; ${RETRIEVE_CMD} ; ${DYAD_CMD} ; ${ASSEMBLE_CMD}; ${MAP_CMD}
all_tasks:
	${MAKE} my_command MY_COMMAND="${ALL_TASKS_CMD}" JOB_PREFIX=${REF_ORG}_${TAXON}_${GENE}

################################################################
## Identify orthologs for a given gene (${GENE}) in the taxon of
## interest (${TAXON})
RESULT_DIR=${MAIN_DIR}/results/${REF_ORG}/${TAXON}
GENE_DIR=${RESULT_DIR}/${GENE}
PREFIX=${GENE_DIR}/${GENE}_${REF_ORG}_${TAXON}
ORTHOLOGS = ${PREFIX}_orthologs.tab
ORTHO_CMD=mkdir -p ${GENE_DIR} ; \
	get-orthologs -q ${GENE} -org ${REF_ORG} -taxon ${TAXON} -uth rank 1 -o ${ORTHOLOGS}
orthologs:
	@echo
	@echo "Getting orthologs	${GENE}		${REF_ORG}	${TAXON}"
	@@echo "${ORTHO_CMD}"
	@${ORTHO_CMD}
	@echo "${ORTHOLOGS}"

################################################################
## Retrieve upstream sequences of set of orthologous genes
SEQ=${PREFIX}_up.fasta
PURGED=${PREFIX}_up_purged.fasta
PURGE_ML=30
RETRIEVE_CMD=retrieve-seq-multigenome -i ${ORTHOLOGS} -o ${SEQ} -noorf ; purge-sequence -i ${SEQ} -o ${PURGED} -ml ${PURGE_ML} -mis 0 -2str
upstream:
	@echo "${RETRIEVE_CMD}"
	@${RETRIEVE_CMD}
	@echo ${SEQ}
	@echo ${PURGED}

################################################################
## Discover over-represented dyads in the promoters of the set of
## orthologous genes
DYADS=${PREFIX}_dyads.tab
DYAD_CMD=dyad-analysis -v ${V} -i ${PURGED} -sort -type any -2str -noov \
		-lth occ 1 -lth occ_sig 0 -return occ,rank,proba -l 3 -spacing 0-20 \
		-o ${DYADS}
dyads:
	@echo "${DYAD_CMD}"
	@${DYAD_CMD}
	@echo ${DYADS}

################################################################
## Count matches between discovered dyads and known sites
KNOWN_SITES=data/CRP_known_sites.tab
count_matches:
	count-matches -file1 ${DYADS} -file2 ${KNOWN_SITES}  -return id,match,weight,seq -lth weight 6 -2str

################################################################
## Assemble the over-represented dyad to form motif
ASSEMBLY=${PREFIX}_dyads.asmb
ASSEMBLE_CMD=pattern-assembly -v ${V} -i ${DYADS} -o ${ASSEMBLY} -subst 0 -2str -maxfl 1 -maxpat 50
assemble:
	@echo "${ASSEMBLE_CMD}"
	@${ASSEMBLE_CMD}
	@echo ${ASSEMBLY}

################################################################
## draw a feature-map with the instances of the over-represented dyads
MAP=${PREFIX}_dyads.png
MAP_CMD=dna-pattern -i ${SEQ} -pl ${DYADS} -limits -origin -0 \
		| convert-features -from dnapat -to ft \
		| feature-map -scorethick -legend -scalebar -scalestep 50 \
		 -format png -title '${GENE} ${REF_ORG} ${TAXON}' -o ${MAP}
map:
	@echo "${MAP_CMD}"
	@${MAP_CMD}
	@echo ${MAP} 

################################################################
## Index the results for the whole set of genes
INDEX_FILE=${MAIN_DIR}/results/${REF_ORG}/${TAXON}/index_${REF_ORG}_${TAXON}.html
index_results:
	@echo "Indexing results	${REF_ORG}	${TAXON}	${INDEX_FILE}"
	@echo "<html>" > ${INDEX_FILE}
	@echo "<body>" >> ${INDEX_FILE}
	@echo "<h1>${REF_ORG} ; ${TAXON}</h1>" >> ${INDEX_FILE}
	@echo "<table border=1 cellpadding=3 cellspacing=3 align=center>" >> ${INDEX_FILE}
#	${MAKE} index_one_result GENE=lexA
	@for g in ${ALL_GENES} ; do  \
		${MAKE} index_one_result GENE=$${g}  ; \
	done
	@echo "</table>" >> ${INDEX_FILE}
	@echo "</body></html>" >> ${INDEX_FILE}
	@echo "Results indexed	${INDEX_FILE}"

################################################################
## Add the analysis of a single gene to the index file
MAX_SIG=`grep -v '^;' ${DYADS} | cut -f 8 | sort -nr | head -1`
index_one_result:
	@echo "<tr>" >> ${INDEX_FILE}
	@echo "<td><a href=${GENE_DIR}>${GENE}</a></td>" >> ${INDEX_FILE}
	@echo "<td>${MAX_SIG}</td>" >> ${INDEX_FILE}
	@echo "<td><a href=${SEQ}>seq</a></td>" >> ${INDEX_FILE}
	@echo "<td><a href=${PURGED}>purged</a></td>" >> ${INDEX_FILE}
	@echo "<td><a href=${DYADS}>dyads</a></td>" >> ${INDEX_FILE}
	@echo "<td><a href=${ASSEMBLY}>assembly</a></td>" >> ${INDEX_FILE}
	@echo "<td><a href=${MAP}>map</a></td>" >> ${INDEX_FILE}
	@echo "</tr>" >> ${INDEX_FILE}


################################################################
## Compare dyads discovered in the different genes
COMPA_DIR=${RESULT_DIR}/comparisons
COMPA_TABLE=${COMPA_DIR}/${REF_ORG}_${TAXON}_dyad_profiles.tab
COMPA_CLASSES=${COMPA_DIR}/${REF_ORG}_${TAXON}_dyad_classes.tab
DYAD_FILE_LIST=${RESULT_DIR}/dyad_files.txt
dyad_file_list:
	@echo
	@echo "Generating the list of dyad files"
	(cd ${RESULT_DIR}; find . -name '*_${REF_ORG}_${TAXON}_dyads.tab'  > ${DYAD_FILE_LIST})
	@echo ${DYAD_FILE_LIST}
	@echo "	`cat ${DYAD_FILE_LIST} | wc -l `	dyad files"

dyad_profiles: 
	@echo
	@echo "Calculating dyad profiles"
	@mkdir -p ${COMPA_DIR}	
	(cd ${RESULT_DIR}; compare-scores -null "NA" -sc 8 \
		-suppress "\./" \
		-suppress "_dyads.tab" \
		-suppress "_${REF_ORG}_${TAXON}" \
		-filelist ${DYAD_FILE_LIST}  \
		| transpose-table \
		| perl -pe 's/\/\S+//' \
		> ${COMPA_TABLE})
	@echo ${COMPA_TABLE}
	@echo "	`grep -v '^;' ${COMPA_TABLE} | grep -v '^#' | wc -l`	profiles"

dyad_classes: 
	@echo
	@echo "Generating dyads/gene file"
	@mkdir -p ${COMPA_DIR}	
	(cd ${RESULT_DIR}; compare-scores -null "NA" -sc 8 \
		-format classes \
		-suppress "\./" \
		-suppress "_dyads.tab" \
		-suppress "_${REF_ORG}_${TAXON}" \
		-filelist ${DYAD_FILE_LIST}  \
		| perl -pe 's/\/\S+//' \
		| sort +1 \
		> ${COMPA_CLASSES})
	@echo ${COMPA_CLASSES}
	@echo "	`grep -v '^;' ${COMPA_CLASSES} | cut -f 2 | sort -u | wc -l`	genes"
	@echo "	`grep -v '^;' ${COMPA_CLASSES} | cut -f 1 | sort -u | wc -l`	dyads"

################################################################
## Compare profiles of dyad significance between each pair of genes
PROFILE_PAIRS=${COMPA_DIR}/${REF_ORG}_${TAXON}_profile_pairs
GENE_PAIRS=${COMPA_DIR}/${REF_ORG}_${TAXON}_gene_pairs
#GENE_PAIRS=boum
gene_pairs:
	@echo
	@echo "Calculating gene pairs"
	compare-classes -v ${V} -i ${COMPA_CLASSES} \
		-return occ,dotprod,jac_sim,proba \
		-sc 3 -lth RandQ 1 -distinct -triangle -sort dotprod \
		-o ${GENE_PAIRS}.tab
	@echo ${GENE_PAIRS}.tab
	@echo "	`grep -v ';' ${GENE_PAIRS}.tab | grep -v '^#' |  wc -l`	gene pairs"

profile_pairs:
	@echo
	@echo "Calculating profile pairs"
	compare-profiles -v ${V} -i ${COMPA_TABLE} -base 2 -distinct -return dotprod -lth AB 1 \
		-o ${PROFILE_PAIRS}.tab
	@echo ${PROFILE_PAIRS}.tab
	@echo "	`grep -v ';' ${PROFiLE_PAIRS}.tab | grep -v '^#' |  wc -l`	profile pairs"

MIN_SCORE=1
SCORE_COL=8
SCORE=dp
PAIR_GRAPH=${GENE_PAIRS}_${SCORE}${MIN_SCORE}
gene_pair_graphs:
	@echo
	@echo "Generating gene pair graphs"
	awk -F '\t' '$$${SCORE_COL} >= ${MIN_SCORE}' ${GENE_PAIRS}.tab \
		| grep -v '^;' \
		| convert-graph -from tab -scol 1 -tcol 3 -wcol ${SCORE_COL} -to gml \
		-o ${PAIR_GRAPH}.gml
	@echo ${PAIR_GRAPH}.dot
	awk -F '\t' '$$${SCORE_COL} >= ${MIN_SCORE}' ${GENE_PAIRS}.tab \
		| convert-graph -from tab -wcol ${SCORE_COL} -to dot \
		-o ${PAIR_GRAPH}.dot
	@echo ${PAIR_GRAPH}.gml

comparisons:
	@${MAKE} dyad_file_list
	@${MAKE} dyad_classes
	@${MAKE} gene_pairs
	@${MAKE} gene_pair_graphs MIN_SCORE=0
	@${MAKE} gene_pair_graphs MIN_SCORE=1
	@${MAKE} gene_pair_graphs MIN_SCORE=2
	@${MAKE} gene_pair_graphs MIN_SCORE=3
	@${MAKE} gene_pair_graphs MIN_SCORE=5
	@${MAKE} gene_pair_graphs MIN_SCORE=10
	@${MAKE} gene_pair_graphs MIN_SCORE=10
	@${MAKE} gene_pair_graphs MIN_SCORE=20
	@${MAKE} gene_pair_graphs MIN_SCORE=30
	@${MAKE} gene_pair_graphs MIN_SCORE=50


################################################################
## Validation of the method : 
## - Discover patterns in upstream regions of all genes from RegulonDB
## - create a graph with genes as nodes and arcs representing co-dicovered
##   patterns
## - compare the discovered graph with the co-factor graph off regulonDB
##   (to be done)

## Select target genes from RegulonDB
REGULONDB_TABLE=data/gene_factor_Escherichia_coli_K12.tab 
REGULONDB_GENES=data/regulondb_genes.tab
regulondb_genes:
	@cut -f 1 ${REGULONDB_TABLE} | sort -u > ${REGULONDB_GENES}
	@echo "${REGULONDB_GENES}"

## Select transcription factors from RegulonDB
REGULONDB_FACTORS=data/regulondb_factors.tab
regulondb_factors:
	@cut -f 2 ${REGULONDB_TABLE} | sort -u  > ${REGULONDB_FACTORS}
	@echo "${REGULONDB_FACTORS}"

## Select genes and transcription factors from RegulonDB, both will be used
## for pattern discovery
REGULONDB_GF=data/regulondb_genes_and_factors.tab
regulondb_gf: regulondb_genes regulondb_factors
	cat ${REGULONDB_GENES} ${REGULONDB_FACTORS} | tr 'a-z' 'A-Z'  | sort -u > ${REGULONDB_GF}
	@echo "${REGULONDB_GF}"

## Run all the tasks for reuglonDB analysis
REGULONDB=`cat ${REGULONDB_GF}| grep -v '^;' | xargs`
REGULON_TAXON=
regulondb_analysis:
	@echo "Analyzing genes from RegulonDB"
	@${MAKE} all_tasks_all_genes REF_ORG=Escherichia_coli_K12 TAXON=${REGULON_TAXON}  ALL_GENES="${REGULONDB}"
	@${MAKE} index_regulondb

## Index the results obtained with RegulonDB
index_regulondb:
	@${MAKE} index_results REF_ORG=Escherichia_coli_K12 TAXON=${REGULON_TAXON} ALL_GENES="${REGULONDB}"


