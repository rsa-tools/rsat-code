################################################################
## Solutions for the exercise "from sites to motifs" with bacterial
## regulons

include ${RSAT}/makefiles/util.mk
MAKEFILE=exo_from_sites_to_motifs.mk


# ## Create a directory for the practicals
# make_dir:
# 	cd $HOME
# 	mkdir -p bioinfostats_2014/MetJ_regulon_analysis
# 	cd  bioinfostats_2014/MetJ_regulon_analysis

# ## Get the full list of supported organisms
# supported:
# 	supported-organisms

## Restrict the list to Bacteria from the genus Escherichia
supported_escherichia:
	supported-organisms -taxon Escherichia

## Count the number of supported species  from the genus Escherichia
supported_escherichia_nb:
	supported-organisms -taxon Escherichia | wc -l

## Select the MG1655 strain (tricky way)
## Result: Escherichia_coli_K_12_substr__MG1655_uid57779
find_mg1655:
	supported-organisms -taxon Escherichia | grep MG1655

## Define an environment variable for readiability and facility when
## typing the subsequent commands.
ORG=Escherichia_coli_K_12_substr__MG1655_uid57779

## Retrieve-seq options to extract the non-coding sequences upstream of
## all the genes, with a maximum length of 400bp (positions relative
## to the CDS), and store them in a fasta file.
SEQ_DIR=results/sequences
SEQ_FILE=${SEQ_DIR}/allup_noorf_${ORG}.fasta
allup_seq:
	@echo
	@echo "Retrieving all upstream sequences from ${ORG}"
	@mkdir -p ${SEQ_DIR}
	retrieve-seq  -org ${ORG} -feattype CDS -type upstream -format fasta -label name -noorf  -imp_pos  -all -o ${SEQ_FILE}
	@echo "	${SEQ_FILE}"

## Statistics on sequence lengths
SEQ_LEN=${SEQ_DIR}/allup_noorf_${ORG}_len
LEN_CI=10
IMG_FORMAT=pdf
allup_len_distrib:
	@echo
	@echo "Computing length distribution"
	sequence-lengths -i ${SEQ_FILE} -o ${SEQ_LEN}.tab
	classfreq -v 1 -ci ${LEN_CI} -i ${SEQ_LEN}.tab -o ${SEQ_LEN}_distrib_ci${LEN_CI}.tab
	XYgraph -i ${SEQ_LEN}_distrib.tab -format ${IMG_FORMAT} \
		-lines -legend \
		-title1 "Distribution of upstream non-coding sequence lengths" \
		-title2 "${ORG}" \
		-xcol 3 -xsize 600 -xleg1 "upstream length" \
		-ycol 4,5,6 -ysize 400 -yleg1 "Number of genes" \
		-format ${IMG_FORMAT} \
		-o ${SEQ_LEN}_distrib_ci${LEN_CI}.${IMG_FORMAT}
	@echo "	${SEQ_LEN}.tab"
	@echo "	${SEQ_LEN}_distrib_ci${LEN_CI}.tab"
	@echo "	${SEQ_LEN}_distrib_ci${LEN_CI}.${IMG_FORMAT}"

## Optional: retrieve an organism-specific background model
bg_choose:
	choose-background-model -org ${ORG} -bg upstream-noorf -type oligo -l 3 -1str -noov

## Store the background file in an environment variable
BG_FILE=`choose-background-model -org ${ORG} -bg upstream-noorf -type oligo -l 3 -1str -ovlp`
bg_file:
	echo "	${BG_FILE}.gz"
	zless ${BG_FILE}.gz

bg_convert:
	convert-background-model -v 1 -i ${BG_FILE} -from oligos -to transitions

## Extract the first matrix from the mme file, and convert it to transfac format (more readable)
MEME_MATRICES=results/meme_motif/MetJ_motif_oops_meme
MATRIX=results/meme_motif/MetJ_motif_oops_meme_top1.tf
matrix_to_transfac:
	@echo
	@echo "Converting MEME result to transfac format"
	convert-matrix -i ${MEME_MATRICES}.txt -top 1 -from meme -to transfac -o ${MATRIX}
	@echo "	 ${MATRIX}"

REGULONDB_MATRICES=${RSAT}/public_html/data/motif_databases/REGULONDB/regulonDB_2012-05.tf 
MEME_MATRICES_VS_REGULONDB=${MEME_MATRICES}_vs_RegulonDB.tab
meme_vs_regulondb:
	compare-matrices -v 1 -format1 meme -file1 ${MEME_MATRICES}.txt \
		-file2 ${REGULONDB_MATRICES} -format2 tf \
		-strand DR \
		-return cor,Ncor,logoDP,NsEucl,NSW,match_rank,matrix_id,matrix_name,width,strand,offset,consensus,alignments_1ton \
		-lth cor 0.7 -lth Ncor 0.4 -uth match_rank 50 -o ${MEME_MATRICES_VS_REGULONDB}


## Convert names of regulated genes to identifiers, to avoid problems
## with synonyms
REFERENCE=regulated_genes
REFERENCE_GENES=data/annotated_genes_MetJ/${REFERENCE}_MetJ
ID_COL=1
_reference_ids:
	@echo
	@echo "Getting gene IDs for ${REFERENCE}"
	gene-info -i ${REFERENCE_GENES}.txt -org ${ORG} | grep -v '#' | sort -u > ${REFERENCE_GENES}_info.tab
	cut -f ${ID_COL} ${REFERENCE_GENES}_info.tab > ${REFERENCE_GENES}_ids.txt
	@echo "	${REFERENCE_GENES}_info.tab"
	@echo "	${REFERENCE_GENES}_ids.txt"

################################################################
## Compare predicted genes to reference genes This target is used with
## different values for the predicted (regexp, scan) and reference
## (genes, first operon gene).
COMMON_GENES=${PREDICTED_GENES}_inter_${REFERENCE}
_predicted_vs_known: _reference_ids
	comm -1 -2 ${REFERENCE_GENES}_ids.txt ${PREDICTED_GENES}_ids.txt \
		| add-linenb | \
		add-gene-info -org ${ORG} -info name,descr \
		-o ${COMMON_GENES}.tab
	@wc -l ${COMMON_GENES}.tab

_predicted_vs_first_operon_genes:
	@${MAKE} _predicted_vs_known REFERENCE=regulated_first_operon_genes


MEME_REGEXP=[CA]TGGA[CT][GA]TCT[AT]AAC[AGT]
REGEXP_RESULT=${SITES_DIR}/predicted_sites_allup_MetJ_meme_regexp
REGEXP_GENES=${REGEXP_RESULT}_genes
regexp_search:
	@echo
	@echo "Scanning promoters wih regular expression from the top meme motif"
	dna-pattern -v 1 -p '${MEME_REGEXP}' -i ${SEQ_FILE} -o ${REGEXP_RESULT}_sites.ft
	grep -v "^;" ${REGEXP_RESULT}_sites.ft | grep -v '^#' | cut -f 4 | sort -u \
		| add-gene-info -org ${ORG} -info id,descr -o ${REGEXP_GENES}_info.tab
	@cut -f 2 ${REGEXP_GENES}_info.tab | sort -u > ${REGEXP_GENES}_ids.txt
	@echo "	${REGEXP_RESULT}_sites.ft"
	@echo "`grep -v '^;' ${REGEXP_RESULT}_sites.ft | grep -v '#' | wc -l` sites found with regexp"
	@echo "	${REGEXP_GENES}_info.tab"
	@wc -l ${REGEXP_GENES}_ids.txt

regexp_vs_known:
	${MAKE} _predicted_vs_known REFERENCE=regulated_genes PREDICTED_GENES=${REGEXP_GENES}
	${MAKE} _predicted_vs_first_operon_genes REFERENCE=regulated_genes PREDICTED_GENES=${REGEXP_GENES}

################################################################
## Use matrix-scan to predict binding sites for the selected factor in
## all upstream sequences. I set the verbosity to 2 in order to get
## messages about progress of the task.
PVAL=0.001
MKV=2
SITES_DIR=results/predicted_sites
SITES=${SITES_DIR}/predicted_sites_allup_MetJ_mkv${MKV}_pval${PVAL}.ft
scan:
	@echo
	@echo "Scanning ${SEQ_FILE} with matrix ${MATRIX}"
	@mkdir -p ${SITES_DIR}
	matrix-scan -v 2 -matrix_format transfac -quick \
		-i ${SEQ_FILE} -seq_format fasta \
		-m ${MATRIX} -pseudo 1 -decimals 1 -n score \
		-2str -origin end \
		-bginput -markov ${MKV} -bg_pseudo 0.01 -return sites -return pval -uth pval ${PVAL}  \
		-o ${SITES}
	@echo "	${SITES}"

## Count the number of predicted sites witha  given p-value threshold
PVAL_THRESHOLDS=0.001 0.0001 0.00001 0.000001
count_sites:
	@for p in ${PVAL_THRESHOLDS}; do \
		${MAKE} _count_sites_one_threshold PVAL_THRESHOLD=$${p} ; \
	done

_count_sites_one_threshold:
	echo "p-value threshold	${PVAL_THRESHOLD}"
	grep -v '^;' ${SITES} | grep -v '^#' | awk -F'\t' '$$9 < ${PVAL_THRESHOLD}'| wc -l

## Count the number of predicted target genes
count_genes:
	for p in ${PVAL_THRESHOLDS}; do \
		${MAKE} _count_genes_one_threshold PVAL_THRESHOLD=$${p} ; \
	done

GENES_DIR=results/predicted_genes
PVAL_THRESHOLD=0.00001
SCAN_PREDICTED_GENES=${GENES_DIR}/predicted_genes_allup_MetJ_mkv${MKV}_pval${PVAL_THRESHOLD}
_count_genes_one_threshold:
	@mkdir -p ${GENES_DIR}
	grep -v '^;' ${SITES} | grep -v '^#' | awk -F'\t' '$$9 <= ${PVAL_THRESHOLD}' \
		| cut -f 1 | sort -u | add-linenb | add-gene-info -org ${ORG} -info id,descr \
		> ${SCAN_PREDICTED_GENES}.tab
	@cut -f 3 ${SCAN_PREDICTED_GENES}.tab | sort -u > ${SCAN_PREDICTED_GENES}_ids.txt
	@wc -l ${SCAN_PREDICTED_GENES}.tab

## Compare the list of predicted target genes to annotated regulated
## genes
compare_genes_all_thresholds:
	@echo
	@echo "Comparing predicted and ${REFERENCE} target genes"
	@for p in ${PVAL_THRESHOLDS}; do \
		${MAKE} _predicted_vs_known PVAL_THRESHOLD=$${p} PREDICTED_GENES=${SCAN_PREDICTED_GENES}; \
		${MAKE} _predicted_vs_first_operon_genes PVAL_THRESHOLD=$${p} PREDICTED_GENES=${SCAN_PREDICTED_GENES}; \
	done


################################################################
## Permute the columns of the matrix of interest, in order to run a
## negative control (measure the empirical rate of false positives).
PERM_MATRIX=results/meme_motif/MetJ_motif_oops_meme_top1_permuted.tf
permute_matrix:
	@echo
	@echo "Permuting matrix ${MATRIX}"
	permute-matrix -i ${MATRIX} \
		-in_format transfac  \
		-out_format transfac  \
		-o ${PERM_MATRIX}
	@echo "	${PERM_MATRIX}"

SITES_PERMMAT=${SITES_DIR}/predicted_sites_allup_MetJ_permuted_mkv${MKV}_pval${PVAL}.ft
scan_permuted_matrix:
	${MAKE} SITES=${SITES_PERMMAT} MATRIX=${PERM_MATRIX} scan

count_sites_permuted_matrix:
	${MAKE} SITES=${SITES_PERMMAT} MATRIX=${PERM_MATRIX} count_sites

################################################################
### Generate artificial random sequences to perform a second negative
### control (scan random sequences with the MetJ motif).
RAND_MKV=6
RAND_SEQ=${SEQ_DIR}/randseq_bg_upstream-noorf_${ORG}_mkv${RAND_MKV}.fasta
randseq:
	@echo
	@echo "Generating random sequences"
	random-seq -format fasta -type DNA -bg upstream-noorf \
		-org ${ORG} -markov ${RAND_MKV} \
		-i ${SEQ_FILE} -template_format fasta \
		-o ${RAND_SEQ}
	@echo "	${RAND_SEQ}"

SITES_RANDSEQ=${SITES_DIR}/predicted_sites_allup_MetJ_permuted_mkv${MKV}_pval${PVAL}.ft
scan_randseq:
	${MAKE} SITES=${SITES_RANDSEQ} SEQ_FILE=${RAND_SEQ} scan

count_sites_randseq:
	${MAKE} SITES=${SITES_RANDSEQ} SEQ_FILE=${RAND_SEQ} count_sites
