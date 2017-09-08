################################################################
## Tester for matrix-from-patterns


include ${RSAT}/makefiles/util.mk
MAKEFILE=${RSAT}/makefiles/matrix-from-patterns_demo.mk

PEAKSET=Sox2
PEAK_SEQ=${RSAT}/public_html/demo_files/${PEAKSET}_peaks.fasta
V=2

################################################################
## Discover over-reprsented k-mers
OLIGO_DIR=results/oligos/${PEAKSET}
MKV=4
OL=6
OLIGOS=${OLIGO_DIR}/peak-motifs_oligos-2str-noov_${OL}nt_mkv${MKV}
oligos:
	@echo ""
	@echo "Discovering over-represented oligonucleotides"
	@mkdir -p ${OLIGO_DIR}
	oligo-analysis  -v ${V} -quick -i ${PEAK_SEQ} -sort -lth ratio 1 \
		-lth occ_sig 0 -uth rank 100 -return occ,proba,rank \
		-2str -noov -seqtype dna -l ${OL} -markov ${MKV} -pseudo 0.01 \
		-o ${OLIGOS}.tab
	@echo "	${OLIGOS}.tab"

assembly:
	@echo
	@echo "Assemblink k-mers"
	pattern-assembly  -v ${V} -i ${OLIGOS}.tab \
		-2str -maxfl 1 -subst 1 -max_asmb_width 20 -toppat 100 -max_asmb_size 50 -max_asmb_width 20 -max_asmb_nb 10 \
		-o ${OLIGOS}.asmb
	@echo "	${OLIGOS}.asmb"

################################################################
## Run matrix-from-patterns
matrices_all:
	@echo
	@echo "Running matrix-from-patterns"
	matrix-from-patterns -v ${V}  -sites -seq ${PEAK_SEQ} \
		-pl ${OLIGOS}.tab \
		-bgfile ${BG_FILE} \
		-toppat 100 -max_asmb_nb 10 -max_asmb_width 20 -subst 1 -prefix oligos_${OL}nt \
		-flanks 2 -collect_method matrix-scan-quick -logo \
		-no_clustering \
		-o ${OLIGOS}_pssm
	@echo "	${OLIGOS}_pssm"

################################################################
## Run matrix-from-patterns with matrix-clustering option in order to avoid redundancy between the motifs
matrices_clustered:
	@echo
	@echo "Running matrix-from-patterns"
	matrix-from-patterns -v ${V}  -sites -seq ${PEAK_SEQ} \
		-pl ${OLIGOS}.tab \
		-bgfile ${BG_FILE} \
		-toppat 100 -max_asmb_nb 10 -max_asmb_width 20 -subst 1 -prefix oligos_${OL}nt \
		-flanks 2 -collect_method matrix-scan-quick -logo \
		-clustering \
		-o ${OLIGOS}_pssm
	@echo "	${OLIGOS}_pssm"

