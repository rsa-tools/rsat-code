################################################################
## Test the rsat subcommands

targets:
	@echo "Targets"
	@echo "	targets			list targets of this makefile"
	@echo "	all			run all the targets (may take some time)"
	@echo "	list_param		list parameters"
	@echo "	randseq			random-seq"
	@echo "	purgeseq		purge-sequence"
	@echo "	download_jaspar		download jaspar PSSM collection"
	@echo "	sequence_lengths	sequence-lengths"
	@echo "	classfreq		classfreq"
	@echo "	xygraph			XYgraph"
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
	@echo "	download_organism	download-organism"
	@echo "	supported_local		supported-organisms"
	@echo "	supported_ensembl	supported-organisms-ensembl -db ensembl"
	@echo "	supported_ensemblg	supported-organisms-ensembl -db ensemblgenomes"
	@echo "	supported_ucsc		supported-organisms-ucsc"
	@echo "	gene_info		gene-info"
	@echo "	add_gene_info		add-gene-info"
	@echo "	retrieve_seq		retrieve-seq"
	@echo "	fetch_sequences		fetch-sequences"

################################################################
## List global parameters
V=1
RESULT_DIR=rsat_subcommand_results
list_param:
	@echo "rsat path	`which rsat`"
	@echo "	RESULT_DIR	${RESULT_DIR}"

################################################################
## Run all targets
all: targets list_param randseq purgeseq download_jaspar sequence_lengths classfreq xygraph retrieve_matrix convert_matrix compare_matrices download_peaks oligos positions assembly matrix_from_patterns create_background matrix_distrib matrix_quality peakmo download_organism supported_local supported_ensembl supported_ensemblg supported_ucsc gene_info add_gene_info retrieve_seq fetch_sequences


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
retrieve_matrix: download_jaspar create_background
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
## Check the lengths of the downloaded peaks
PEAK_LENGTHS=${PEAKS}_lengths.tsv
sequence_lengths: download_peaks
	@echo "Testing sequence-lengths"
	rsat sequence-lengths -i ${PEAKS} -o ${PEAK_LENGTHS}
	@echo "	PEAK_LENGTHS	${PEAK_LENGTHS}"

################################################################
## Compute the distribution of peak sequence lengths with classfreq
PEAK_LEN_DISTRIB=${PEAKS}_length_distrib.tsv
classfreq: sequence_lengths
	@echo "Testing classfreq"
	rsat classfreq -i ${PEAK_LENGTHS} -v ${V} -ci 10 -o ${PEAK_LEN_DISTRIB}
	@echo "	PEAK_LEN_DISTRIB	${PEAK_LEN_DISTRIB}"

################################################################
## XYgraph: plot the distribution of peak sequence lengths
PEAK_LEN_DISTRIB_GRAPH=${PEAKS}_length_distrib.png
xygraph: classfreq
	@echo "Testing XYgraph"
	rsat XYgraph -i ${PEAK_LEN_DISTRIB} \
		-title 'Peak length distribution' \
		-xsize 800 -ysize 400 -lines \
		-xcol 3 -ycol 4,5,6 -format png \
		-xleg1 "Peak length" -yleg1 "Frequencies" -legend -ymin 0 \
		-xgstep1 50 -xgstep2 10 \
		-ygstep1 200 -ygstep2 50 \
		-o ${PEAK_LEN_DISTRIB_GRAPH}
	@echo "	PEAK_LEN_DISTRIB_GRAPH	${PEAK_LEN_DISTRIB_GRAPH}"

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
	rsat create-background-model -v ${V} \
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
	rsat matrix-distrib -v ${V} \
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


################################################################
## Download an organism from an RSAT server
ORGANISM=Saccharomyces_cerevisiae
download_organism:
	@echo "Testing download-organism"
	@echo "	ORGANISM	${ORGANISM}"
	rsat download-organism -v ${V} -org ${ORGANISM}

################################################################
## List the organisms supported on the local RSAT instance
SUPPORTED_DIR=${RESULT_DIR}/supported-organisms-x_results
SUPPORTED_LOCAL=${SUPPORTED_DIR}/supported-organisms_local.tsv
supported_local:
	@echo "Testing supported-organisms"
	@mkdir -p ${SUPPORTED_DIR}
	@echo "	SUPPORTED_DIR	${SUPPORTED_DIR}"
	rsat supported-organisms -v ${V} -o ${SUPPORTED_LOCAL}
	@echo "	SUPPORTED_LOCAL	${SUPPORTED_LOCAL}"

################################################################
## List the organisms supported at ensembl (http://ensembl.org)
ENSEMBL_DB=ensembl
SUPPORTED_ENSEMBL=${SUPPORTED_DIR}/supported-organisms-ensembl_db-${ENSEMBL_DB}.tsv
supported_ensembl:
	@echo "Testing supported-organisms-ensembl"
	@mkdir -p ${SUPPORTED_DIR}
	@echo "	SUPPORTED_DIR	${SUPPORTED_DIR}"
	rsat supported-organisms-ensembl -v ${V} -db ensembl -o ${SUPPORTED_ENSEMBL}
	@echo "	ENSEMBL_DB		${ENSEMBL_DB}"
	@echo "	SUPPORTED_ENSEMBL	${SUPPORTED_ENSEMBL}"
	@echo "	Number of organisms	`grep -v ';' ${SUPPORTED_ENSEMBL} | wc -l`"

################################################################
## List the organisms supported at ensembl (http://ensembl.org)
SUPPORTED_ENSEMBLGENOMES=${SUPPORTED_DIR}/supported-organisms-ensembl_db-${ENSEMBL_DB}.tsv
supported_ensemblg:
	@echo "Testing supported-organisms-ensembl with ensemblgenomes"
	@mkdir -p ${SUPPORTED_DIR}
	@echo "	SUPPORTED_DIR	${SUPPORTED_DIR}"
	rsat supported-organisms-ensembl -v ${V} -db ensemblgenomes -o ${SUPPORTED_ENSEMBLGENOMES}
	@echo "	SUPPORTED_ENSEMBLGENOMES	${SUPPORTED_ENSEMBLGENOMES}"
	@echo "	Number of organisms	`grep -v ';' ${SUPPORTED_ENSEMBLGENOMES} | wc -l`"

################################################################
## List the organisms supported at ensembl (http://ensembl.org)
## NOTE: THIS SCRIPT IS STILL NOT FINISHED (JvH, 2020-01-04)
TAXID=4751
SUPPORTED_ENSEMBLGENOMES_TAXID=${SUPPORTED_DIR}/supported-organisms-ensemblgenomes_branch${TAXID}.tsv
supported_ensemblg_taxid:
	@echo "Testing supported-organisms-ensemblgenomes"
	@mkdir -p ${SUPPORTED_DIR}
	@echo "	SUPPORTED_DIR	${SUPPORTED_DIR}"
	rsat supported-organisms-ensemblgenomes -v ${V} \
		-query_type branch -q ${TAXID} \
		-o ${SUPPORTED_ENSEMBLGENOMES_TAXID}
	@echo "	TAXID				${TAXID}"
	@echo "	SUPPORTED_ENSEMBLGENOMES_TAXID	${SUPPORTED_ENSEMBLGENOMES_TAXID}"

################################################################
## List the organisms supported at UCSC genome browser (https://genome.ucsc.edu/)
SUPPORTED_UCSC=${SUPPORTED_DIR}/supported-organisms-ucsc.tsv
supported_ucsc:
	@echo "Testing supported-organisms-ucsc"
	@mkdir -p ${SUPPORTED_DIR}
	@echo "	SUPPORTED_DIR	${SUPPORTED_DIR}"
	rsat supported-organisms-ucsc -v ${V} -o ${SUPPORTED_UCSC}
	@echo "	SUPPORTED_UCSC	${SUPPORTED_UCSC}"
	@echo "	Number of organisms	`grep -v ';' ${SUPPORTED_UCSC} | wc -l`"

################################################################
## Get inforamtion about genes for an organism installed locally
GENE_INFO_DIR=${RESULT_DIR}/gene-info_result
FEATTYPE=gene
GENE_INFO=${GENE_INFO_DIR}/${ORGANISM}_MET_${FEATTYPE}_info.tsv
gene_info:
	@echo "Testing gene-info"
	@mkdir -p ${GENE_INFO_DIR}
	@echo "	GENE_INFO_DIR	${GENE_INFO_DIR}"
	rsat gene-info -org ${ORGANISM} -q 'MET\d+' -feattype ${FEATTYPE} -o ${GENE_INFO}
	@echo "	GENE_INFO	${GENE_INFO}"

################################################################
## Add user-specified information to a gene table
GENE_INFO_ADDED=${GENE_INFO_DIR}/${ORGANISM}_MET_${FEATTYPE}_info-added.tsv
add_gene_info: gene_info
	@echo "Testing add-gene-info"
	rsat add-gene-info -i ${GENE_INFO} \
		-org ${ORGANISM} \
		-info upstr_size,downstr_size,names \
		-col 1 -feattype ${FEATTYPE} \
		-o ${GENE_INFO_ADDED}
	@echo "	GENE_INFO_ADDED	${GENE_INFO_ADDED}"

################################################################
## Retrieve sequences from a locally installed organism
RETRIEVE_SEQ_DIR=${RESULT_DIR}/retrieve-seq_result
RETRIEVED_SEQ=${RETRIEVE_SEQ_DIR}/${ORGANISM}_MET_${FEATTYPE}_upstream-noorf.fasta
retrieve_seq: gene_info
	@echo "Testing retrieve-seq"
	@mkdir -p ${RETRIEVE_SEQ_DIR}
	@echo "	RETRIEVE_SEQ_DIR	RETRIEVE_${SEQ_DIR}"
	rsat retrieve-seq -org ${ORGANISM} \
		-i ${GENE_INFO} \
		-type upstream -noorf -feattype ${FEATTYPE} -label id,name \
		-o ${RETRIEVED_SEQ}
	@echo "	RETRIEVED_SEQ		${RETRIEVED_SEQ}"

################################################################
## Fetch sequences from UCSC genome browser
FETCHED_SEQ_DIR=${RESULT_DIR}/fetch-sequences_result
CEBPA_PEAKS_BASENAME=fetch-sequences_Schmidt_2011_mm9_CEBPA_SWEMBL_R0.12_702peaks
CEBPA_PEAKS_COORD_URL=http://metazoa.rsat.eu/demo_files/${CEBPA_PEAKS_BASENAME}.bed
CEBPA_PEAKS_BED=${FETCHED_SEQ_DIR}/${CEBPA_PEAKS_BASENAME}.bed
CEBPA_PEAKS_SEQ=${FETCHED_SEQ_DIR}/${CEBPA_PEAKS_BASENAME}.fasta
fetch_sequences:
	@echo "Testing fetch-sequences"
	@mkdir -p ${FETCHED_SEQ_DIR}
	@echo "	FETCHED_SEQ_DIR	${FETCHED_SEQ_DIR}"
	@echo "	CEBPA_PEAKS_COORD_URL	${CEBPA_PEAKS_COORD_URL}"
	wget --no-clobber ${CEBPA_PEAKS_COORD_URL} -O ${CEBPA_PEAKS_BED}
	@echo "	CEBPA_PEAKS_BED	${CEBPA_PEAKS_BED}"
	rsat fetch-sequences -i ${CEBPA_PEAKS_BED} -genome mm9 -o ${CEBPA_PEAKS_SEQ}
	@echo "	CEBPA_PEAKS_SEQ	${CEBPA_PEAKS_SEQ}"

