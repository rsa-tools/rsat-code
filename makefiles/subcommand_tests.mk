################################################################
## Test the rsat subcommands

targets:
	@echo "Targets"
	@echo "	targets			list targets of this makefile"
	@echo "	list_param		list parameters"
	@echo "	randseq			random-seq"
	@echo "	purgeseq		purge-sequence"
	@echo "	download_jaspar		download jaspar PSSM collection"
	@echo "	retrieve_matrix		retrieve-matrix"
	@echo "	convert_matrix		convert-matrix"
	@echo "	compare_matrices	compare-matrices"
	@echo "	download_peaks		download test peaks"
	@echo "	oligos			oligo-analysis"
	@echo "	positions		position-analysis"
	@echo "	assembly		pattern-assembly"
	@echo "	matrix_from_patterns	matrix-from-patterns"
	@echo "	create_background	create-background-model"
	@echo "	matrix_distrib		matrix-distrib"
	@echo "	matrix_quality		matrix-quality"
	@echo "	peakmo			peak-motifs"

################################################################
## List global parameters
V=1
RESULT_DIR=rsat_subcommand_results
list_param:
	@echo "rsat path	`which rsat`"
	@echo "	RESULT_DIR	${RESULT_DIR}"

################################################################
## Generate random sequences
RANDSEQ_DIR=${RESULT_DIR}/random-seq_result
RANDSEQ=${RANDSEQ_DIR}/randseq_l1000_n10.fasta
randseq:
	@echo "Testing random-seq"
	@mkdir -p ${RANDSEQ_DIR}
	@echo "	RANDSEQ_DIR	${RANDSEQ_DIR}"
	@echo "	Generating random sequences"
	rsat random-seq -l 1000 -n 10 -seed 123 -o ${RANDSEQ}
	@echo "	RANDSEQ		${RANDSEQ}"
#	md5sum ${RANDSEQ} > ${RANDSEQ}.md5
#	@echo "	md5sum		${RANDSEQ}.md5"

################################################################
## Download JASPAR non-redundant vertebrate collection
JASPAR_FILE=JASPAR2020_CORE_vertebrates_non-redundant_pfms.tf
JASPAR_URL=http://teaching.rsat.eu/motif_databases/JASPAR/Jaspar_2020/nonredundant/${JASPAR_FILE}
JASPAR_DIR=${RESULT_DIR}/data/jaspar
JASPAR=${JASPAR_DIR}/${JASPAR_FILE}
download_jaspar:
	@echo "Downloading JASPAR NR vertebrates"
	@mkdir -p ${JASPAR_DIR}
	@if [ -f ${JASPAR} ] ; \
	then echo "	JASPAR file already there"; \
	else wget --no-clobber ${JASPAR_URL} -O ${JASPAR}; \
	fi
	@echo "	JASPAR	${JASPAR}"


################################################################
## retrieve-matrix test
## Retrieve selected matrices from JASPAR
RETRIEVE_MATRIX_DIR=${RESULT_DIR}/retrieve-matrix_result
MATRIX_BASENAME=sox-oct_matrices
MATRICES=${RETRIEVE_MATRIX_DIR}/${MATRIX_BASENAME}.tf
retrieve_matrix: download_jaspar
	@echo "Testing retrieve-matrix"
	@mkdir -p ${RETRIEVE_MATRIX_DIR}
	rsat retrieve-matrix -v ${V} -i ${JASPAR} -id POU5F1 -id SOX2 -id Pou5f1::Sox2 \
		-o ${MATRICES}
	@echo "	MATRICES	${MATRICES}"


################################################################
## convert-matrix (including the generation of logos)
CONVERT_MATRIX_DIR=${RESULT_DIR}/convert-matrix_result
CONVERTED_MATRICES=${CONVERT_MATRIX_DIR}/${MATRIX_BASENAME}.tab
convert_matrix: download_jaspar
	@echo "Testing convert-matrix"
	@mkdir -p ${CONVERT_MATRIX_DIR}
	rsat convert-matrix -v ${V} -i ${MATRICES} -from transfac -to tab \
		-return weights,margins,logo -decimals 2 \
		-o ${CONVERTED_MATRICES}
	@echo "	CONVERTED_MATRICES	${CONVERTED_MATRICES}"

################################################################
## Compare matrices
MATRIX_COMPA_DIR=${RESULT_DIR}/compare-matrices_result
MATRIX_COMPA=${MATRIX_COMPA_DIR}/sox-oct_vs_itself
compare_matrices: retrieve_matrix
	@echo "Testing compare-matices"
	@mkdir -p ${MATRIX_COMPA_DIR}
	rsat compare-matrices -v ${V} -file ${MATRICES} -format transfac \
		-mode matches -distinct -strand DR \
		-lth cor 0.6 -lth Ncor 0.33 -uth match_rank 50 \
		-return cor,Ncor,logoDP,match_rank,matrix_id,matrix_name,width,consensus,alignments_1ton \
		-o ${MATRIX_COMPA}
	@echo "	MATRIX_COMPA_DIR	${MATRIX_COMPA_DIR}"
	@echo "	MATRIX_COMPA		${MATRIX_COMPA}_index.html"

################################################################
## matrix-clustering test
MATRIX_CLUSTERING_DIR=${RESULT_DIR}/matrix-clustering_result
matrix_clustering:
	@echo "Testing matrix-clustering"
	rsat matrix-clustering -h
	@echo "	MATRIX_CLUSTERING_DIR	${MATRIX_CLUSTERING_DIR}"


################################################################
## Download peak sequences for tests
PEAK_BASENAME=peak-motifs_demo
PEAK_FILE=${PEAK_BASENAME}.fa
PEAK_URL=http://teaching.rsat.eu//demo_files/${PEAK_FILE}
PEAKMO_DIR=${RESULT_DIR}/peak-motifs_result
PEAKS=${PEAKMO_DIR}/${PEAK_FILE}
download_peaks:
	@echo "	Downloading peak sequences from ${PEAK_URL}"
	@mkdir -p ${PEAKMO_DIR}
	@if [ -f ${PEAKS} ] ; \
	then echo "	Peak file already there"; \
	else wget --no-clobber ${PEAK_URL} -O ${PEAKS}; \
	fi
	@echo "	PEAKS	${PEAKS}"

################################################################
## Purge sequences to mask repeats
PURGESEQ_DIR=${RESULT_DIR}/purge-sequence_result
PURGED_PEAKS=${PURGESEQ_DIR}/${PEAK_BASENAME}_purged.fa
purgeseq: download_peaks
	@echo "Testing purge-seq"
	@mkdir -p ${PURGESEQ_DIR}
	@echo "	PURGESEQ_DIR	${PURGESEQ_DIR}"
	@echo "	Purging sequences"
	rsat purge-sequence -i ${PEAKS} -o ${PURGED_PEAKS}
	@echo "	PURGED_PEAKS		${PURGED_PEAKS}"

################################################################
## oligo-analysis test
OLIGO_DIR=${RESULT_DIR}/oligo-analysis_result
OLIGO_BASENAME=${PEAK_BASENAME}_6nt_2str_noov_sig0
OLIGO_PREFIX=${OLIGO_DIR}/${OLIGO_BASENAME}
OLIGOS=${OLIGO_PREFIX}.tsv
oligos: download_peaks
	@echo "Testing oligo-analysis"
	@mkdir -p ${OLIGO_DIR}
	rsat oligo-analysis -v ${V} -i ${PURGED_PEAKS} \
		-l 6 -2str -noov \
		-return occ,freq,proba,rank \
		-markov 4 -lth occ_sig 0 \
		-o ${OLIGOS}
	@echo "	OLIGO_DIR	${OLIGO_DIR}"
	@echo "	OLIGOS		${OLIGOS}"

################################################################
## position-analysis test
##
## Includes options to cluster k-mers based on their positional
## profiles, which reveals 2 clusters:
## 1. k-mers concentrated in the middle of the peaks
## 2. k-mers avoided in the middle of the peaks
##
POSITION_DIR=${RESULT_DIR}/position-analysis_result
POSITION_BASENAME=${POSITION_DIR}/${PEAK_BASENAME}_6nt_ci25
POSITIONS=${POSITION_BASENAME}.tsv
positions: download_peaks
	@echo "Testing position-analysis"
	@mkdir -p ${POSITION_DIR}
	rsat position-analysis -v ${V} -i ${PURGED_PEAKS} \
		-l 6 -2str -noov -ci 25 -lth_sig 0 -lth_occ 1 \
		-clust_nb 2 -max_asmb_per_cluster 2 \
		-origin center -maxpos 500 -minpos -500 \
		-sort -return chi,sig,rank,distrib,clusters,matrices,graphs,index \
		-o ${POSITIONS}
	@echo "	POSITION_DIR	${POSITION_DIR}"
	@echo "	POSITIONS	${POSITIONS}"
	@echo "	index		${POSITION_BASENAME}_index.html"
	@echo "	graph index	${POSITION_BASENAME}_graph_index.html"

################################################################
## Test pattern-assembly with oligos
ASSEMBLY=${OLIGO_PREFIX}_assembly.txt
assembly:
	@echo "Testing pattern-assembly"
	rsat pattern-assembly -v ${V} \
		-i ${OLIGOS} \
		-maxfl 1 -2str -subst 1 -match 5 \
		-o ${ASSEMBLY}
	@echo "	ASSEMBLY	${ASSEMBLY}"

################################################################
## Test matrix-from-patterns
OLIGO_MATRICES=${OLIGO_PREFIX}_pssm
matrix_from_patterns:
	@echo "Testing matrix-from-patterns"
	rsat matrix-from-patterns -seq ${PEAKS} \
		-pl ${OLIGOS} \
		-sc 9 -subst 1 -maxfl 1 -match 5 -2str \
		-sites -collect_method matrix-scan-quick -flanks 3 -logo  \
		-o ${OLIGO_MATRICES}
	@echo "	OLIGO_MATRICES		${OLIGO_MATRICES}"
	@echo "	PSSM (transfac format)	${OLIGO_MATRICES}_count_matrices.tf"
	@echo "	Sites (feature format)	${OLIGO_MATRICES}_sig_sites.ft"



################################################################
## Test create-background-model
BG_MKV=1
BG_DIR=${RESULT_DIR}/background-models
BG_FORMAT=oligos
BG_FILE=${BG_DIR}/${PEAK_BASENAME}_bg-model_markov${BG_MKV}_${BG_FORMAT}.tsv
create_background: download_peaks
	@echo "Testing create-background-model"
	@mkdir -p ${BG_DIR}
	rsat create-background-model -v 1 \
		-i ${PEAKS} \
		-markov ${BG_MKV} -out_format ${BG_FORMAT} \
		-o ${BG_FILE}
	@echo "	BG_DIR	${BG_DIR}"
	@echo "	BG_FILE	${BG_FILE}"

################################################################
## Test matrix-distrib
MATRIX_DISTRIB_DIR=${RESULT_DIR}/matrix-distrib_result
MATRIX_DISTRIB_PREFIX=${MATRIX_DISTRIB_DIR}/${MATRIX_BASENAME}_distrib
MATRIX_DISTRIB=${MATRIX_DISTRIB_PREFIX}.tsv
matrix_distrib:
	@echo "Testing matrix-distrib"
	@mkdir -p ${MATRIX_DISTRIB_DIR}
	rsat matrix-distrib -v 1 \
		-top 1 \
		-m ${MATRICES} -matrix_format transfac \
		-pseudo 1 -decimals 1 \
		-bg_format ${BG_FORMAT} \
		-bgfile ${BG_FILE} \
		-bg_pseudo 0.01 \
		-o ${MATRIX_DISTRIB}
	@echo "	MATRICES		${MATRICES}"
	@echo "	BG_FILE			${BG_FILE}"
	@echo "	MATRIX_DISTRIB_DIR	${MATRIX_DISTRIB_DIR}"
	@echo "	MATRIX_DISTRIB		${MATRIX_DISTRIB}"
	rsat XYgraph -i ${MATRIX_DISTRIB} \
		-title1 "Distribution of weights" \
		-title2 "Score probability" \
		-xcol 1 -ycol 2 -legend -lines -pointsize 1 \
		-xleg1 "Weight" \
		-yleg1 "Frequency" \
		-format pdf -r_plot \
		-o ${MATRIX_DISTRIB_PREFIX}_weigh-distrib.pdf
	@echo "	Weight distrib graph	 ${MATRIX_DISTRIB_PREFIX}_weigh-distrib.pdf"
	rsat XYgraph \
		-i ${MATRIX_DISTRIB} \
		-title1 'Distribution of weights  (log scale)' \
		-title2 'Score probability and P-value' \
		-xcol 1 \
		-ycol 2,4 \
		-legend \
		-lines \
		-pointsize 1 -ylog \
		-xleg1 'weight' -yleg1 'Frequency (log scale)' \
		-r_plot \
		-format pdf \
		-o ${MATRIX_DISTRIB_PREFIX}_weigh-distrib_ylog.pdf
	@echo "	Weight distrib graph (ylog) ${MATRIX_DISTRIB_PREFIX}_weigh-distrib_ylog.pdf"


################################################################
## Test matrix-quality
QUALITY_DIR=${RESULT_DIR}/matrix-quality_result
QUALITY_PREFIX=${QUALITY_DIR}/${OLIGO_BASENAME}_quality
matrix_quality:
	@echo "Testing matrix-quality"
	@mkdir -p ${QUALITY_DIR}
	rsat matrix-quality -v ${V} \
		-ms ${OLIGO_MATRICES}_count_matrices.tf \
		-matrix_format transfac \
		-seq peaks ${PEAKS} -seq_format fasta \
		-bgfile ${BG_FILE} -bg_format ${BG_FORMAT} \
		-kfold 10 -perm peaks 1 \
		-o ${QUALITY_PREFIX}
	@echo "	QUALITY_DIR	${QUALITY_DIR}"
	@echo "	index file	${QUALITY_PREFIX}_synthesis.html"

################################################################
## peak-motifs test
PEAKMO_TASK=purge,seqlen,composition,disco,merge_motifs,split_motifs,motifs_vs_motifs,timelog,synthesis,small_summary,scan,motifs_vs_db
peakmo: download_jaspar download_peaks
	@echo "Testing peak-motifs"
	@mkdir -p ${PEAKMO_DIR}
	@echo "	Running peak-motifs"
	rsat peak-motifs  -v 3 \
		-title Oct4_Chen2008_sites_from_Jaspar \
		-i ${PEAKS} \
		-markov auto \
		-disco oligos,positions \
		-nmotifs 5 -minol 6 -maxol 6 -no_merge_lengths -2str -origin center \
		-scan_markov 1 -source galaxy \
		-prefix peak-motifs -noov -img_format png \
		-outdir ${PEAKMO_DIR} \
		-motif_db JASPAR transfac ${JASPAR} \
		-task ${PEAKMO_TASK}  \
		&> ${PEAKMO_DIR}/peak-motifs_log.txt
	@echo "	Log file	${PEAKMO_DIR}/peak-motifs_log.txt"
	@echo "	Result page	${PEAKMO_DIR}/peak-motifs_synthesis.html"


