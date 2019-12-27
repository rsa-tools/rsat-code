################################################################
## Test the rsat subcommands

targets:
	@echo "Targets"
	@echo "	targets		list targets of this makefile"
	@echo "	list_param	list parameters"
	@echo "	randseq		random-seq test"
	@echo "	peakmo		peak-motifs test"

RESULT_DIR=rsat_subcommand_results
list_param:
	@echo "rsat path	`which rsat`"
	@echo "	RESULT_DIR	${RESULT_DIR}"

RANDSEQ_DIR=${RESULT_DIR}/random-seq_result
RANDSEQ=${RANDSEQ_DIR}/randseq_l1000_n10.fasta
randseq:
	@echo "Testing random-seq"
	@mkdir -p ${RANDSEQ_DIR}
	@echo "	RANDSEQ_DIR	${RANDSEQ_DIR}"
	@echo "	Generating random sequences"
	rsat random-seq -l 1000 -n 10 -seed 123 -o ${RANDSEQ}
	@echo "	RANDSEQ		${RANDSEQ}"


PEAKMO_DIR=${RESULT_DIR}/peak-motifs_result
JASPAR_URL=http://teaching.rsat.eu/motif_databases/JASPAR/Jaspar_2020/nonredundant/JASPAR2020_CORE_vertebrates_non-redundant_pfms.tf
peakmo:
	@echo "Testing peak-motifs"
	@echo "	Downloading peak sequences"
	wget --no-clobber http://teaching.rsat.eu//demo_files/peak-motifs_demo.fa
	@echo "	Downloading JASPAR NR vertebrates"
	wget --no-clobber ${JASPAR_URL}
	@echo "	Running peak-motifs"
	rsat peak-motifs  -v 3 \
		-title Oct4_Chen2008_sites_from_Jaspar \
		-i peak-motifs_demo.fa \
		-markov auto \
		-disco oligos,positions \
		-nmotifs 5 -minol 6 -maxol 6 -no_merge_lengths -2str -origin center \
		-scan_markov 1 -source galaxy \
		-prefix peak-motifs -noov -img_format png \
		-outdir ${PEAKMO_DIR} \
		-task purge,seqlen,composition,disco,merge_motifs,split_motifs,motifs_vs_motifs,timelog,synthesis,small_summary,scan  \
		&> peak-motifs_log.txt
	@echo "	Log file	peak-motifs_log.txt"
	@echo "	Result page	peak-motifs_result/peak-motifs_synthesis.html"


