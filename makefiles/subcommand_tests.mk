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
MATRIX_COMPA_DIR=${RESULT_DIR}/compare-matrices_results
MATRIX_COMPA=${MATRIX_COMPA_DIR}/sox-oct_vs_itself
compare_matrices: retrieve_matrix
	@echo "Testing compare-matices"
	@mkdir -p ${MATRIX_COMPA_DIR}
	rsat compare-matrices -v ${V} -file ${MATRICES} -format transfac \
		-mode matches -distinct -strand DR \
		-lth cor 0.7 -lth Ncor 0.4 -uth match_rank 50 \
		-return cor,Ncor,logoDP,match_rank,matrix_id,matrix_name,width,consensus,alignments_1ton \
		-o ${MATRIX_COMPA}
	@echo "	MATRIX_COMPA_DIR	${MATRIX_COMPA_DIR}"
	@echo "	MATRIX_COMPA		${MATRIX_COMPA}_index.html"

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
## peak-motifs test
PEAKMO_TASK=purge,seqlen,composition,disco,merge_motifs,split_motifs,motifs_vs_motifs,timelog,synthesis,small_summary,scan
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


