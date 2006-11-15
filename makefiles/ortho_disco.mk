## ##############################################################
## Discover conserved motifs in upstream sequences of orthologs for each gene
## of a genome, and infer a graph of co-regulation by comparing discovered
## motifs between each pair of genes

## THE FOLLOWING VARIABLES HAVE TO BE DEFINED FOR EACH NEW ANALYSIS
## GENE
## REF_ORG
## TAXON
## MAIN_DIR
##

## Default parameters
GENE=CRP
REF_ORG=Escherichia_coli_K12
TAXON=Gammaproteobacteria
MAIN_DIR=.

include ${RSAT}/makefiles/util.mk
MAKEFILE=${RSAT}/makefiles/ortho_disco.mk
MAKE=make -s -f ${MAKEFILE} REF_ORG=${REF_ORG} TAXON=${TAXON} GENE=${GENE} MAIN_DIR=${MAIN_DIR}

## Verbosity level
V=1

## ##############################################################
## List the parameters
RESULT_DIR=${MAIN_DIR}/results/${REF_ORG}/${TAXON}
SEQ_TAX_DIR=${MAIN_DIR}/data/sequences/ortho_promoters/${REF_ORG}/${TAXON}
SEQ_ORG_DIR=${MAIN_DIR}/data/sequences/gene_promoters/${REF_ORG}
GENE_DIR=${RESULT_DIR}/motifs/${GENE}
list_parameters:
	@echo ""
	@echo "REF_ORG      ${REF_ORG}"
	@echo "TAXON        ${TAXON}"
	@echo "MAIN_DIR     ${MAIN_DIR}"
	@echo "SEQ_ORG_DIR  ${SEQ_ORG_DIR}"
	@echo "SEQ_TAX_DIR  ${SEQ_TAX_DIR}"
	@echo "RESULT_DIR   ${RESULT_DIR}"
	@echo "GENE_DIR     ${GENE_DIR}"

################################################################
## List of all genes to be analyzed. 
## By default, all the genes of the considered organism (REF_ORG)
ALL_GENES=`grep -v '^--' ${RSAT}/data/genomes/${REF_ORG}/genome/cds_names.tab | grep primary | cut -f 2 | sort -u  | perl -pe 's/\n/ /g'`
list_all_genes:
	@echo ""
	@echo "All genes"
	@echo "${ALL_GENES}"


################################################################
## Run all the tasks on all the genes
all_tasks_all_genes:
	@echo "All tasks applied on all genes	${REF_ORG}	${TAXON}"
#	@${MAKE} list_all_genes ALL_GENES='${ALL_GENES}'
	@${MAKE} list_parameters
	@echo ""
	@for g in ${ALL_GENES} ; do \
		${MAKE} all_tasks REF_ORG=${REF_ORG} TAXON=${TAXON} MAIN_DIR=${MAIN_DIR} GENE=$${g} ; \
	done

## Retrieve ortholog clusters and upstream sequences for all the genes
all_sequences:
	@echo 
	@echo "Retrieving sequences for all genes	${REF_ORG}	${TAXON}"
	@for g in ${ALL_GENES} ; do \
		${MAKE} ortho_seq_tasks REF_ORG=${REF_ORG} TAXON=${TAXON} MAIN_DIR=${MAIN_DIR} GENE=$${g} ; \
	done


## Analyze dyads in all the genes (sequences must have been retrieved before)
all_dyads:
	@echo
	@echo "Analyzing dyads for all genes	${REF_ORG}	${TAXON}"
	@for g in ${ALL_GENES} ; do \
		${MAKE} dyad_tasks REF_ORG=${REF_ORG} TAXON=${TAXON} MAIN_DIR=${MAIN_DIR} GENE=$${g} ; \
	done

## Analyze dyads with organism-specific filter in all the genes (sequences must have been retrieved before)
all_dyads_filtered:
	@echo
	@echo "Analyzing dyads  with organism-specific filter for all genes	${REF_ORG}	${TAXON}"
	@for g in ${ALL_GENES} ; do \
		${MAKE} _dyad_tasks_filtered REF_ORG=${REF_ORG} TAXON=${TAXON} MAIN_DIR=${MAIN_DIR} GENE=$${g} ; \
	done
#	${MAKE} all_tasks_all_genes ALL_TASKS_CMD="${FILTER_DYADS_CMD}; ${DYAD_CMD} ; ${ASSEMBLE_CMD}; ${MAP_CMD}" \
#		JOB_PREFIX=${REF_ORG}_${TAXON}_${GENE} \


################################################################
## Apply one task to all the genes
GENE_TASK=match_known_sites
iterate_genes:
	@for g in ${ALL_GENES} ; do \
		${MAKE} ${GENE_TASK} REF_ORG=${REF_ORG} TAXON=${TAXON} MAIN_DIR=${MAIN_DIR} GENE=$${g} ; \
	done

################################################################
## Run all the tasks for a single gene
ALL_TASKS_CMD=${ORTHO_CMD} ; ${RETRIEVE_CMD} ; ${DYAD_CMD} ; ${ASSEMBLE_CMD}; ${MAP_CMD}
#ALL_TASKS_CMD=${FILTER_DYADS_CMD}
all_tasks:
	${MAKE} my_command MY_COMMAND="${ALL_TASKS_CMD}" JOB_PREFIX=${REF_ORG}_${TAXON}_${GENE}

## get orthologs + upstream sequences for a single gene
ortho_seq_tasks:
	${MAKE} my_command MY_COMMAND="${ORTHO_CMD} ; ${RETRIEVE_CMD}" JOB_PREFIX=${REF_ORG}_${TAXON}_${GENE}

## analyze dyads
dyad_tasks:
	${MAKE} my_command MY_COMMAND="${DYAD_CMD} ; ${ASSEMBLE_CMD}; ${MAP_CMD}" JOB_PREFIX=${REF_ORG}_${TAXON}_${GENE}

_dyad_tasks_filtered:
	${MAKE} my_command \
		MY_COMMAND="${FILTER_DYADS_CMD}; ${DYAD_CMD} ; ${ASSEMBLE_CMD}; ${MAP_CMD}" JOB_PREFIX=${REF_ORG}_${TAXON}_${GENE} \
		FILTER_SUFFIX=_filtered DYAD_OPT='-accept ${DYAD_FILTER}' 

################################################################
## Identify orthologs for a given gene (${GENE}) in the taxon of
## interest (${TAXON})
PREFIX=${GENE}_${REF_ORG}_${TAXON}
ORTHOLOGS = ${ORTHO_DIR}/${PREFIX}_orthologs.tab
ORTHO_DIR=${MAIN_DIR}/data/orthologs/${REF_ORG}/${TAXON}
ORTHO_CMD=mkdir -p ${ORTHO_DIR} ; \
	get-orthologs -q ${GENE} -org ${REF_ORG} -taxon ${TAXON} -uth rank 1 -o ${ORTHOLOGS} -return all
ORTHO_NB=`cat  ${ORTHOLOGS} | grep -v '^\#' | grep -v '^;' | wc -l `
orthologs:
	@echo
	@echo "Getting orthologs	${GENE}		${REF_ORG}	${TAXON}"
	@@echo "${ORTHO_CMD}"
	@${ORTHO_CMD}
	@echo "${ORTHOLOGS}"
	${MAKE} ortho_nb

ortho_nb:
	@echo "${ORTHO_NB} orthologs found for gene ${GENE} in ${TAXON}"

################################################################
## Retrieve upstream sequences of set of orthologous genes
NOORF=-noorf
#SEQ_TAX_DIR=${GENE_DIR}
SEQ=${SEQ_TAX_DIR}/${PREFIX}_up.fasta
PURGED=${SEQ_TAX_DIR}/${PREFIX}_up_purged.fasta
PURGE_ML=30
RETRIEVE_CMD=mkdir -p ${SEQ_TAX_DIR}; retrieve-seq-multigenome -feattype CDS,tRNA,rRNA -i ${ORTHOLOGS} -o ${SEQ} ${NOORF} ; purge-sequence -i ${SEQ} -o ${PURGED} -ml ${PURGE_ML} -mis 0 -2str -mask_short ${PURGE_ML}; gzip -f ${SEQ}; gzip -f ${PURGED}
upstream:
	@echo "${RETRIEVE_CMD}"
	@${RETRIEVE_CMD}
	@echo ${SEQ}
	@echo ${PURGED}

################################################################
## Discover over-represented dyads in the promoters of the set of
## orthologous genes
BG=monads
STR=-2str
NOOV=-noov
RETURN=occ,freq,proba,rank
SUFFIX=${STR}${NOOV}_${BG}_dyads${FILTER_SUFFIX}
DYADS=${GENE_DIR}/${PREFIX}${SUFFIX}
DYAD_CMD=mkdir -p ${GENE_DIR}; \
	dyad-analysis -v ${V} -i ${PURGED} -sort -type any ${STR} ${NOOV} \
		-lth occ 1 -lth occ_sig 0 -return ${RETURN} -l 3 -spacing 0-20 \
		-bg ${BG} -org ${REF_ORG} \
		-o ${DYADS}.tab ${DYAD_OPT} 
dyads:
	@echo "${DYAD_CMD}"
	@${DYAD_CMD}
	@echo ${DYADS}.tab

################################################################
## Run dyad analysis using the sequence of the reference organism as filter
FILTER_SEQ=${SEQ_ORG_DIR}/${GENE}_${REF_ORG}_up${NOORF}.fasta.gz
DYAD_FILTER=${GENE_DIR}/${GENE}_${REF_ORG}_dyad_filter
#DYADS_FILTERED=${GENE_DIR}/${PREFIX}${SUFFIX}_filtered
FILTER_DYADS_CMD= \
	mkdir -p ${SEQ_ORG_DIR} ; retrieve-seq -feattype CDS,tRNA,rRNA  -org ${REF_ORG} -q ${GENE} ${NOORF} -o ${FILTER_SEQ} ; echo 'Filter sequence	${FILTER_SEQ}'; \
	dyad-analysis -v 0 -i ${FILTER_SEQ} -type any ${STR} ${NOOV} -lth occ 1 -return occ -l 3 -spacing 0-20 -o ${DYAD_FILTER} ; echo 'Dyad filter	${DYAD_FILTER}'
filter_dyads:
	@echo
	@echo 'Filtering dyads	${GENE}	${TAXON}'
	@${FILTER_DYADS_CMD}

dyads_filtered:
	@${MAKE} filter_dyads
	@${MAKE} -s dyads FILTER_SUFFIX=_filtered DYAD_OPT='-accept ${DYAD_FILTER}' ; echo 'Filtered dyads	${DYADS}_filtered.tab' 
	@${MAKE} -s assemble FILTER_SUFFIX=_filtered ; echo 'Filtered dyad assembly	${DYADS}_filtered.asmb'
	@${MAKE} -s map FILTER_SUFFIX=_filtered ; echo 'Filtered dyad map	${DYADS}_filtered.png' 

################################################################
## Count matches between discovered dyads and known sites
KNOWN_SITES=data/sites_per_gene.tab
KNOWN_SITES_GENE=${GENE_DIR}/${GENE}_gene_known_sites.tab
MIN_W=6
KNOWN_SITE_MATCHES=${DYADS}_vs_known_sites_w${MIN_W}
COMPARE_PATTERNS_CMD= \
	awk '$$5 == "${GENE}" {print $$1"\t"$$3"_"$$2"_"$$4}' ${KNOWN_SITES} > ${KNOWN_SITES_GENE}.tab ; \
	perl -pe 's|[a-z]||g' ${KNOWN_SITES_GENE}.tab  > ${KNOWN_SITES_GENE}_noflanks.tab ; \
	compare-patterns -slide -v 1 -file2 ${DYADS}.tab -file1 ${KNOWN_SITES_GENE}.tab  -return id,match,strand,offset,weight,rel_w,Pval,Eval_p,sig_p,Eval_f,sig_f,seq -lth weight ${MIN_W} -2str -o ${KNOWN_SITE_MATCHES}.tab ; \
	compare-patterns -slide -v 1 -file2 ${DYADS}.tab -file1 ${KNOWN_SITES_GENE}_noflanks.tab  -return id,match,strand,offset,weight,rel_w,Pval,Eval_p,sig_p,Eval_f,sig_f,seq -lth weight ${MIN_W} -2str -o ${KNOWN_SITE_MATCHES}_noflanks.tab ; \
	compare-patterns -slide -v 1 -file1 ${DYADS}.tab -file2 ${KNOWN_SITES_GENE}.tab  -table weight -null "." -lth weight ${MIN_W} -2str -o ${KNOWN_SITE_MATCHES}_weight_table.tab ; \
	compare-patterns -slide -v 1 -file1 ${DYADS}.tab -file2 ${KNOWN_SITES_GENE}_noflanks.tab  -table weight -null "." -lth weight ${MIN_W} -2str -o ${KNOWN_SITE_MATCHES}_noflanks_weight_table.tab 
match_known_sites:
	@echo 
	@echo "Matching known sites for gene	${GENE}"
#	@echo "${COMPARE_PATTERNS_CMD}"
	${COMPARE_PATTERNS_CMD}
	@echo ${KNOWN_SITES_GENE}.tab
	@echo ${KNOWN_SITE_MATCHES}.tab
	@echo ${KNOWN_SITE_MATCHES}_weight_table.tab
	@echo ${KNOWN_SITES_GENE}_noflanks.tab
	@echo ${KNOWN_SITE_MATCHES}_noflanks.tab
	@echo ${KNOWN_SITE_MATCHES}_noflanks_weight_table.tab

################################################################
## Assemble the over-represented dyad to form motif
ASSEMBLY=${DYADS}.asmb
ASSEMBLE_CMD=pattern-assembly -v ${V} -i ${DYADS}.tab -o ${ASSEMBLY} -subst 0 ${STR} -maxfl 1 -maxpat 100
assemble:
	@echo "${ASSEMBLE_CMD}"
	@${ASSEMBLE_CMD}
	@echo ${ASSEMBLY}

################################################################
## draw a feature-map with the instances of the over-represented dyads
MAP=${DYADS}.png
MAP_CMD=dna-pattern -i ${SEQ} -pl ${DYADS}.tab -limits -origin -0 \
		| convert-features -from dnapat -to ft \
		| feature-map -scorethick -legend -scalebar -scalestep 50 \
		 -format png -title '${GENE} ${REF_ORG} ${TAXON}' -o ${MAP}
map:
	@echo "${MAP_CMD}"
	@${MAP_CMD}
	@echo ${MAP} 

################################################################
## Index the results for the whole set of genes
index_results: index_results_html index_results_tab 

#INDEX_FILE=${MAIN_DIR}/results/${REF_ORG}/${TAXON}/index_${REF_ORG}_${TAXON}${SUFFIX}.tab
#index_results_tab:
#	@echo "Indexing results	${REF_ORG}	${TAXON}	${INDEX_FILE}"
#	@echo "; ${REF_ORG} ; ${TAXON} ; ${SUFFIX}" > ${INDEX_FILE}
#	@for g in ${ALL_GENES} ; do  \
#		${MAKE} index_one_result_tab GENE=$${g}  ; \
#	done
#	@echo "Results indexed"
#	@echo ${INDEX_FILE}

################################################################
## Index the results in a "tab-delimted html file" (without the HTML
## table, which is heavy to load for browsers when there are 5000
## rows)
INDEX_FILE=${MAIN_DIR}/results/${REF_ORG}/${TAXON}/index_${REF_ORG}_${TAXON}${SUFFIX}.html
index_results_tab:
	@echo "Indexing results	${REF_ORG}	${TAXON}	${INDEX_FILE}"
	@echo "<html>" > ${INDEX_FILE}
	@echo "<body>" >> ${INDEX_FILE}
	@echo "<h1>${REF_ORG} ; ${TAXON} ; ${SUFFIX}</h1>" >> ${INDEX_FILE}
	@echo "<PRE>" >> ${INDEX_FILE}
	@for g in ${ALL_GENES} ; do  \
		${MAKE} index_one_result GENE=$${g}  ; \
	done
	@echo "</PRE>" >> ${INDEX_FILE}
	@echo "</body></html>" >> ${INDEX_FILE}
	@echo "Results indexed"
	@echo ${INDEX_FILE}

## Add the analysis of a single gene to the index file
SIG_COLUMN=9
MAX_SIG=`grep -v '^;' ${DYADS}.tab | cut -f ${SIG_COLUMN} | sort -nr | grep -v "Binary" | head -1`
MAX_SIG_FILTERED=`grep -v '^;' ${DYADS}_filtered.tab | cut -f ${SIG_COLUMN} | sort -nr | grep -v "Binary" | head -1`
ROOT_DIR=${RESULT_DIR}/
#MATCH_LINKS="	<a href=${KNOWN_SITE_MATCHES}.tab>matches</a>	<a href=${KNOWN_SITE_MATCHES}_weight_table.tab>match table</a>"
MATCH_LINKS=
index_one_result:
	@echo "Indexing result for gene	${GENE}"
	@echo "${GENE}	<a href=${SEQ}>seq</a>	<a href=${PURGED}>purged</a>	${MAX_SIG}	<a href=${DYADS}.tab>dyads</a>	<a href=${ASSEMBLY}>asmb</a>	<a href=${MAP}>map</a>	${MATCH_LINKS}" >> ${INDEX_FILE}


################################################################
## Index the results for the whole set of genes
## Index all results in a HTML table
## This is convenient to read, but heavy to load when there are 5000 genes
INDEX_TABLE=${MAIN_DIR}/results/${REF_ORG}/${TAXON}/index_table_${REF_ORG}_${TAXON}${SUFFIX}.html
index_results_html:
	@echo "Indexing results	${REF_ORG}	${TAXON}	${INDEX_TABLE}"
	@echo "<html>" > ${INDEX_TABLE}
	@echo "<body>" >> ${INDEX_TABLE}
	@echo "<h1>${REF_ORG} ; ${TAXON} ; ${SUFFIX}</h1>" >> ${INDEX_TABLE}
	@echo "<table border=1 cellpadding=3 cellspacing=3 align=center>" >> ${INDEX_TABLE}
#	${MAKE} index_one_result GENE=lexA
	@for g in ${ALL_GENES} ; do  \
		${MAKE} index_one_result_html GENE=$${g}  ; \
	done
	@echo "</table>" >> ${INDEX_TABLE}
	@echo "</body></html>" >> ${INDEX_TABLE}
	@echo "Results indexed	${INDEX_TABLE}"

################################################################
## Index one result in a HTML table
## This is convenient to read, but heavy to load when there are 5000 genes
index_one_result_html:
	@echo "Indexing result for gene	${GENE}"
	@echo "<tr>" >> ${INDEX_TABLE}
	@echo "<td><a href=${GENE_DIR}>${GENE}</a></td>" | perl -pe 's|${ROOT_DIR}||' >> ${INDEX_TABLE}
	echo "<td width=5></td>" >> ${INDEX_TABLE} ; \
	if [ -f "${DYADS}.tab" ] ; then \
		echo "<td>${MAX_SIG}</td>"  >> ${INDEX_TABLE} ; \
	else echo "<td></td>" >> ${INDEX_TABLE} ; \
	fi
	if [ -f "${DYADS}.tab" ] ; then \
		echo "<td><a href=${DYADS}.tab>dyads</a></td>"  | perl -pe 's|${ROOT_DIR}||' >> ${INDEX_TABLE} ; \
	else echo "<td></td>" >> ${INDEX_TABLE} ; \
	fi
	if [ -f "${ASSEMBLY}" ] ; then \
		echo "<td><a href=${ASSEMBLY}>assembly</a></td>"  | perl -pe 's|${ROOT_DIR}||' >> ${INDEX_TABLE} ; \
	else echo "<td></td>" >> ${INDEX_TABLE} ; \
	fi
	if [ -f "${MAP}" ] ; then \
		echo "<td><a href=${MAP}>map</a></td>"  | perl -pe 's|${ROOT_DIR}||' >> ${INDEX_TABLE} ; \
	else echo "<td></td>" >> ${INDEX_TABLE} ; \
	fi
	if [ -f "${KNOWN_SITE_MATCHES}.tab" ] ; then \
		echo "<td><a href=${KNOWN_SITE_MATCHES}.tab>matches</a></td>"  | perl -pe 's|${ROOT_DIR}||' >> ${INDEX_TABLE} ; \
	else echo "<td></td>" >> ${INDEX_TABLE} ; \
	fi
	if [ -f "${KNOWN_SITE_MATCHES}_weight_table.tab" ] ; then \
		echo "<td><a href=${KNOWN_SITE_MATCHES}_weight_table.tab>match table</a></td>"  | perl -pe 's|${ROOT_DIR}||' >> ${INDEX_TABLE} ; \
	else echo "<td></td>" >> ${INDEX_TABLE} ; \
	fi
	echo "<td width=5></td>" >> ${INDEX_TABLE} ; \
	if [ -f "${DYADS}_filtered.tab" ] ; then \
		echo "<td>${MAX_SIG_FILTERED}</td>"  >> ${INDEX_TABLE} ; \
	else echo "<td></td>" >> ${INDEX_TABLE} ; \
	fi
	if [ -f "${DYADS}_filtered.tab" ] ; then \
		echo "<td><a href=${DYADS}_filtered.tab>f.dyads</a></td>"  | perl -pe 's|${ROOT_DIR}||' >> ${INDEX_TABLE} ; \
	else echo "<td></td>" >> ${INDEX_TABLE} ; \
	fi
	if [ -f "${DYADS}_filtered.asmb" ] ; then \
		echo "<td><a href=${DYADS}_filtered.asmb>f.assembly</a></td>"  | perl -pe 's|${ROOT_DIR}||' >> ${INDEX_TABLE} ; \
	else echo "<td></td>" >> ${INDEX_TABLE} ; \
	fi
	if [ -f "${DYADS}_filtered.png" ] ; then \
		echo "<td><a href=${DYADS}_filtered.png>f.map</a></td>"  | perl -pe 's|${ROOT_DIR}||' >> ${INDEX_TABLE} ; \
	else echo "<td></td>" >> ${INDEX_TABLE} ; \
	fi
	echo "<td width=5></td>" >> ${INDEX_TABLE} ; \
	if [ -f "${ORTHOLOGS}" ] ; then \
		echo "<td><a href=${ORTHOLOGS}>${ORTHO_NB} ortho</a></td>"  | perl -pe 's|${ROOT_DIR}||' >> ${INDEX_TABLE} ; \
	else echo "<td></td>" >> ${INDEX_TABLE} ; \
	fi
	if [ -f "${SEQ}" ] ; then \
		echo "<td><a href=${SEQ}>seq</a></td>"  | perl -pe 's|${ROOT_DIR}||' >> ${INDEX_TABLE} ; \
	else echo "<td></td>" >> ${INDEX_TABLE} ; \
	fi
	if [ -f "${PURGED}" ] ; then \
		echo "<td><a href=${PURGED}>purged</a></td>"  | perl -pe 's|${ROOT_DIR}||' >> ${INDEX_TABLE} ; \
	else echo "<td></td>" >> ${INDEX_TABLE} ; \
	fi
	if [ -f "${SEQ}.gz" ] ; then \
		echo "<td><a href=${SEQ}.gz>seq</a></td>"  | perl -pe 's|${ROOT_DIR}||' >> ${INDEX_TABLE} ; \
	else echo "<td></td>" >> ${INDEX_TABLE} ; \
	fi
	if [ -f "${PURGED}.gz" ] ; then \
		echo "<td><a href=${PURGED}.gz>purged</a></td>"  | perl -pe 's|${ROOT_DIR}||' >> ${INDEX_TABLE} ; \
	else echo "<td></td>" >> ${INDEX_TABLE} ; \
	fi
	@echo "</tr>" >> ${INDEX_TABLE}


################################################################
## Compare dyads discovered in the different genes
COMPA_DIR=${RESULT_DIR}/gene_pair_network
DYAD_FILE_LIST=${RESULT_DIR}/files${SUFFIX}.txt
dyad_file_list:
	@echo
	@echo "Generating the list of dyad files	${REF_ORG}	${TAXON}"
	find ${RESULT_DIR} -name '*${SUFFIX}.tab' \
		> ${DYAD_FILE_LIST}
	@echo ${DYAD_FILE_LIST}
	@echo "	`cat ${DYAD_FILE_LIST} | wc -l `	dyad files"

## ##############################################################
## Generate dyad profiles: one row per gene, one column per dyad,
## cells indicate significance.
## Beware : the resulting matrix is very sparse (most values are NA),
## and takes A LOT of space (>100Mb for Pseudomonas only).  It is not
## necessary to use it anymore, the target "dyad_classes" gives
## similar results but is MUCH faster and takes MUCH less space
DYAD_PROFILES=${COMPA_DIR}/${REF_ORG}_${TAXON}${SUFFIX}_dyad_profiles.tab
dyad_profiles: dyad_file_list
	@echo
	@echo "Calculating dyad profiles	${REF_ORG}	${TAXON}"
	@mkdir -p ${COMPA_DIR}	
	compare-scores -null "NA" -sc ${SIG_COLUMN} \
		-suppress "${RESULT_DIR}/motifs/" \
		-suppress "_dyads.tab" \
		-suppress "_${REF_ORG}_${TAXON}" \
		-filelist ${DYAD_FILE_LIST}  \
		| transpose-table \
		| perl -pe 's/\/\S+//g' \
		> ${DYAD_PROFILES}
	@echo ${DYAD_PROFILES}
	@echo "	`grep -v '^;' ${DYAD_PROFILES} | grep -v '^#' | wc -l`	dyad profiles"


## ##############################################################
## Compare each pair of genes on the basis of their dyad profiles, and
## select pairs of genes profiles having a significant intersection
## between their dyad profiles.
## 
## This is computationally much heavier than gene_pairs, I just leave
## it for history and cross-validation
profile_pairs:
	@echo
	@echo "Calculating profile pairs	${REF_ORG}	${TAXON}"
	compare-profiles -v ${V} -i ${DYAD_PROFILES} -base 2 -distinct -return dotprod -lth AB 1 \
		-o ${PROFILE_PAIRS}.tab
	@echo ${PROFILE_PAIRS}.tab
	@echo "	`grep -v ';' ${PROFiLE_PAIRS}.tab | grep -v '^#' |  wc -l`	profile pairs"

## ##############################################################
## Collect the significant dyads for each gene and store the result in
## a single file with 3 columns:
## - dyad
## - gene
## - significance
DYAD_CLASSES=${COMPA_DIR}/${REF_ORG}_${TAXON}${SUFFIX}_dyad_classes.tab
DYAD_CLASS_CMD= mkdir -p ${COMPA_DIR}	 ; \
	compare-scores -null NA -sc ${SIG_COLUMN} \
		-format classes \
		-suppress '${RESULT_DIR}/motifs/' \
		-suppress '_dyads.tab' \
		-suppress '_${REF_ORG}_${TAXON}' \
		-filelist ${DYAD_FILE_LIST}  \
		| perl -pe 's/\/\S+//g' \
		| sort +1 \
		> ${DYAD_CLASSES}
dyad_classes: dyad_file_list
	@echo
	@echo "Generating dyads/gene file	${REF_ORG}	${TAXON}"
	${MAKE} my_command MY_COMMAND="${DYAD_CLASS_CMD}"
	@echo ${DYAD_CLASSES}
	@echo "	`grep -v '^;' ${DYAD_CLASSES} | cut -f 2 | sort -u | wc -l`	genes"
	@echo "	`grep -v '^;' ${DYAD_CLASSES} | cut -f 1 | sort -u | wc -l`	dyads"

## ##############################################################
## Compare discovered dyads between each pair of gene and return the
## results as a table with one row per pair of genes, and different
## significance statistics.
PROFILE_PAIRS=${COMPA_DIR}/${REF_ORG}_${TAXON}${SUFFIX}_profile_pairs
GENE_PAIRS=${COMPA_DIR}/${REF_ORG}_${TAXON}${SUFFIX}_gene_pairs
GENE_PAIR_RETURN=occ,dotprod,jac_sim,proba,entropy,rank
V2=3
GENE_PAIR_CMD= \
	compare-classes -v ${V2} -i ${DYAD_CLASSES} \
		-return ${GENE_PAIR_RETURN} \
		-sc 3 -lth QR 1 -distinct -triangle -sort dotprod \
		-o ${GENE_PAIRS}.tab ; \
	text-to-html -chunk 200 -i ${GENE_PAIRS}.tab -o  ${GENE_PAIRS}.html -font variable 
gene_pairs:
	@echo ; \
	@echo "Calculating gene pairs	${REF_ORG}	${TAXON}" ; \
	${MAKE} my_command MY_COMMAND="${GENE_PAIR_CMD}"
	echo ${GENE_PAIRS}.tab
	@echo ${GENE_PAIRS}.html 
	@echo "	`grep -v ';' ${GENE_PAIRS}.tab | grep -v '^#' |  wc -l`	gene pairs"


################################################################
## Convert the gene pairs into a graph (2 formats : .dot and .gml)
MIN_SCORE=1
SCORE_COL=15
SCORE=dp
PAIR_GRAPH=${GENE_PAIRS}_${SCORE}${MIN_SCORE}
PAIR_GRAPH_CMD=	grep -v '^;' ${GENE_PAIRS}.tab \
		| awk -F '\t' '$$${SCORE_COL} >= ${MIN_SCORE}'  \
		| convert-graph -from tab -scol 2 -tcol 3 -wcol ${SCORE_COL} -to dot \
		-o ${PAIR_GRAPH}.dot ; \
	grep -v '^;' ${GENE_PAIRS}.tab \
		| awk -F '\t' '$$${SCORE_COL} >= ${MIN_SCORE}' \
		| grep -v '^;' \
		| convert-graph -from tab -scol 2 -tcol 3 -wcol ${SCORE_COL} -to gml \
		-o ${PAIR_GRAPH}.gml

gene_pair_graph:
	@echo
	@echo "Generating gene pair graph	${REF_ORG}	${TAXON}"
	${MAKE} my_command MY_COMMAND="${PAIR_GRAPH_CMD}"
	@echo ${PAIR_GRAPH}.gml
	@echo ${PAIR_GRAPH}.dot

## ##############################################################
## This task combines dyad_classes, gene_pairs and pair_graph, in order tos
## end them in one shot to the cluster.
gene_network:
#	${MAKE} my_command MY_COMMAND="${DYAD_CLASS_CMD} ${GENE_PAIR_CMD}; ${PAIR_GRAPH_CMD}"
	${MAKE} my_command MY_COMMAND="${DYAD_CLASS_CMD} ${GENE_PAIR_CMD}; ${PAIR_GRAPH_CMD}"

################################################################
## Generate gene pair graphs with various levels of threshold
PAIR_GRAPH_SCORES=0 1 2 3 5 10 20 30 50
gene_pair_graph_series:
	@for s in ${PAIR_GRAPH_SCORES} ; do \
		${MAKE} gene_pair_graph MIN_SCORE=$${s} ; \
	done

## ##############################################################
## Same as gene_pairs but return the dyads characterizing each pair of
## gene (intersection, differences)
gene_pairs_members:
	${MAKE} gene_pairs GENE_PAIR_RETURN=occ,dotprod,jac_sim,proba,members,rank \
		GENE_PAIRS=${COMPA_DIR}/${REF_ORG}_${TAXON}${SUFFIX}_gene_pairs_dyads

################################################################
## Compare dot product, hypergeometric sig, adn jaccard distance
gene_pair_score_compa:
	@echo Comparing scores
	XYgraph -i ${GENE_PAIRS}.tab \
		-ylog 10 \
		-xlog 10 \
		-format png -xcol 8 -ycol 12 \
		-size 300 \
		-pointsize 4 \
		-xleg1 "dot product" \
		-yleg1 "hypergeometric sig = -log10(E-value)" \
		-o ${GENE_PAIRS}_dp_vs_sig.png
	XYgraph -i ${GENE_PAIRS}.tab \
		-ylog 10 \
		-xlog 10 \
		-format png -xcol 8 -ycol 7 \
		-size 300 \
		-pointsize 4 \
		-xleg1 "dot product" \
		-yleg1 "Jaccard distance" \
		-o ${GENE_PAIRS}_dp_vs_jaccard.png
	XYgraph -i ${GENE_PAIRS}.tab \
		-ylog 10 \
		-xlog 10 \
		-format png -xcol 12 -ycol 7 \
		-size 300 \
		-pointsize 4 \
		-xleg1 "hypergeometric sig = -log10(E-value)" \
		-yleg1 "Jaccard distance" \
		-o ${GENE_PAIRS}_sig_vs_jaccard.png
	@echo ${GENE_PAIRS}_sig_vs_jaccard.png

## ##############################################################
## Use MCL to extract clusters from the graph of pairwse relationships
## between genes showing similar over-represened motifs in their
## upstream sequences.
## MCL can be obtained at: http://micans.org/mcl/
#GRAPH=coli_gamma_common_motifs_dotprod
INFLATION=2.0
MCL_DIR=${RESULT_DIR}/clusters/mcl
MCL_FILE=${MCL_DIR}/${REF_ORG}_${TAXON}${SUFFIX}_sc${SCORE_COL}_mcl_I${INFLATION}
GRAPH_FILE=${GENE_PAIRS}_sc${SCORE_COL}.tab
mcl:
	@echo
	@mkdir -p ${MCL_DIR}
	grep -v '^;' ${GENE_PAIRS}.tab \
		| grep -v '^\#' \
		| awk '{print $$2"\t"$$3"\t"$$${SCORE_COL}}' \
		> ${GRAPH_FILE} ; echo ${GRAPH_FILE} 
	mcl ${GRAPH_FILE} --abc -I ${INFLATION} -o ${MCL_FILE}.mic >& mcl_log.txt 
	convert-classes -from mcl -to tab -i ${MCL_FILE}.mic -o ${MCL_FILE}.tab ; echo ${MCL_FILE}.tab 
	convert-graph -from tab -to gml -i ${MCL_FILE}.tab -o ${MCL_FILE}.gml ; echo ${MCL_FILE}.gml
	convert-graph -from tab -to dot -i ${MCL_FILE}.tab -o ${MCL_FILE}.dot ; echo ${MCL_FILE}.dot 


################################################################
## Test various values for the inflation parameter of mcl
INFLATION_VALUES=1.2 1.5 \
	2.0 2.5 \
	3.0 3.5 \
	4.0 4.5 \
	5.0 5.5
INFLATION_TASK=mcl
mcl_inflation_series:
	@for i in ${INFLATION_VALUES}; do \
		echo ; \
		echo "Inflation $${i}" ; \
		${MAKE} -s ${INFLATION_TASK} INFLATION=$${i}; \
	done

################################################################
## Compare the clustering result to RegulonDB (for E.coli only)
MCL_VS_REG=${MCL_FILE}__vs__regulonDB
mcl_vs_regulondb:
	cat ${REGULONDB_TABLE}.tab | tr a-z A-Z > ${REGULONDB_TABLE}_uc.tab
	compare-classes -v ${V} -r ${REGULONDB_TABLE}_uc.tab -q ${MCL_FILE}.tab -return occ,percent,proba,members,rank -lth QR 2 -sort sig -o ${MCL_VS_REG}.tab
	@echo ${MCL_VS_REG}.tab
	@text-to-html -i  ${MCL_VS_REG}.tab -o  ${MCL_VS_REG}.html -font variable -chunk 250
	@echo ${MCL_VS_REG}.html
	compare-classes -v ${V} -r ${REGULONDB_TABLE}_uc.tab -q ${MCL_FILE}.tab -return occ -matrix QR -o ${MCL_VS_REG}_conting.tab
	@echo ${MCL_VS_REG}_conting.tab


################################################################
## Dun all the post-discovery tasks
## - dyad classes
## - gene pairs
## - gene pair graphs
## - clustering on the gene pairs
comparisons: 
	@echo "This target is obsolete. Try	make gene_pair_network"

gene_pair_network:
	@${MAKE} dyad_file_list
	@${MAKE} dyad_classes
	@${MAKE} gene_pairs
	@${MAKE} gene_pair_graph_series
	@${MAKE} mcl_inflation_series


################################################################
## Motif discovery using MEME
MEME=meme
MEME_MOD=anr
MEME_NMOTIFS=5
MEME_EVT=1
MEME_MINW=6
MEME_MAXW=25
MEME_SUFFIX=_MEME_${MEME_MOD}_nmotifs${MEME_NMOTIFS}_evt${MEME_EVT}_minw${MEME_MINW}_maxw${MEME_MAXW}
MEME_FILE=${GENE_DIR}/${PREFIX}${MEME_SUFFIX}
#meme SREBP_I/SREBP_I_up2000_mrna_purged.fasta -revcomp -dna -mod anr -minw 6 -maxw 25 
MEME_CMD=gunzip -f ${PURGED}; \
	${MEME} ${PURGED} -dna -mod ${MEME_MOD} -nmotifs ${MEME_NMOTIFS} \
	-evt ${MEME_EVT} -revcomp -text \
	-minw ${MEME_MINW} -maxw ${MEME_MAXW} \
	${MEME_OPT} > ${MEME_FILE}; \
	gzip -f ${PURGED}
#dyad-analysis -v ${V} -i ${PURGED} -sort -type any ${STR} ${NOOV} \
#		-lth occ 1 -lth occ_sig 0 -return ${RETURN} -l 3 -spacing 0-20 \
#		-bg ${BG} -org ${REF_ORG} \
#		-o ${DYADS}.tab
meme:
	@echo
	@echo "Starting MEME for gene ${GENE}"
	${MAKE} my_command MY_COMMAND="${MEME_CMD}" JOB_PREFIX=meme_${REF_ORG}_${TAXON}_${GENE}
	@echo ${MEME_FILE}

## match MEME motifs with MAST
MAST=mast
MAST_FILE=${MEME_FILE}_mast
MAST_CMD=gunzip -f ${SEQ} ; \
	${MAST} ${MEME_FILE} -d ${SEQ} -mev ${MEME_EVT} -stdout > ${MAST_FILE}.html; \
	mast ${MEME_FILE} -d ${SEQ} -mev ${MEME_EVT} -text -stdout > ${MAST_FILE}.txt; \
	gzip -f ${SEQ}
mast:
	@echo
	@echo "Starting MAST for gene ${GENE}"
	${MAKE} my_command MY_COMMAND="${MAST_CMD}" JOB_PREFIX=mast_${REF_ORG}_${TAXON}_${GENE}
	@echo ${MAST_FILE}.html
	@echo ${MAST_FILE}.txt
