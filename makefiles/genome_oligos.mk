################################################################
## Detect over- and under-represented oligonucleotides in full genomes
## and in different sequence types (coding, upstream, intergenic, ...).

include ${RSAT}/makefiles/util.mk
MAKEFILE=${RSAT}/makefiles/genome_oligos.mk

all: seqs seq_len_distrib all_seq_all_markov clean_seqs


## Retrieve sequences of different types
##ORG=Neisseria_meningitidis_Z2491
genome_segments:
	@echo ""
	@echo "Retrieving sequences"
	mkdir -p ${SEQ_DIR}
	(cd ${SEQ_DIR}; coding-or-not -org ${ORG} -return seq,ncs,orf,cs,div,conv,tandem)
	convert-seq -i ${RSAT}/data/genomes/${ORG}/genome/contigs.txt -from filelist -to wc -lw 0 -o ${SEQ_DIR}/${ORG}_genome.wc
	gzip -f ${SEQ_DIR}/*.wc

NOORF=-noorf
ORG=Saccharomyces_cerevisiae
SEQ_DIR=data/${ORG}/sequences
SEQ_TYPE=upstream
SEQ_PREFIX=${ORG}_${SEQ_TYPE}
SEQ_FILE=${SEQ_DIR}/${SEQ_PREFIX}.wc.gz
seqs:
	@echo "Retrieving all sequences	${ORG}	${SEQ_TYPE}"
	@mkdir -p ${SEQ_DIR}
	retrieve-seq -all -org ${ORG} ${NOORF} -type ${SEQ_TYPE} -format wc -lw 0 -nocomment -label id,name -o ${SEQ_FILE}
	@echo ${SEQ_FILE}

## Delete sequence files
clean_seqs:
	@echo ""
	@echo "Cleaning sequences"
	rm -rf ${SEQ_DIR}

################################################################
## Calculate sequence length distribution for one sequence type
LEN_DIR=data/${ORG}/sequence_lengths/
LEN_FILE=${LEN_DIR}/${SEQ_PREFIX}_len.tab
one_seq_type_len_distrib:
	@mkdir -p ${LEN_DIR}
	sequence-lengths -i ${SEQ_FILE} -format wc -o ${LEN_FILE}
	@echo ${LEN_FILE}
	cut -f 2  ${LEN_FILE} | classfreq -v 1 -min 0 -o ${LEN_DISTRIB}.tab
	@echo ${LEN_DISTRIB}.tab

################################################################
## Plot th histogram of a sequence length distribution
LEN_DISTRIB=${LEN_DIR}/${SEQ_PREFIX}_len_distrib
one_seq_type_len_histo:
	XYgraph -i ${LEN_DISTRIB}.tab \
		-lines -xlog 2 \
		-o ${LEN_DISTRIB}.jpg \
		-xcol 3 -ycol 4,5,6 \
		-title1 '${ORG} ${SEQ_TYPE}' \
		-title2 'sequence length distribution' \
		-xlab 'length' \
		-ylab 'number of sequences'
	@echo ${LEN_DISTRIB}.jpg

################################################################
## Calculate sequence length distributions for all sequence types
seq_len_distrib:
	@echo ""
	@echo "Calculating sequence lengths"
	@${MAKE} one_seq_type_len_distrib SEQ_TYPE=genome
	@for type in intergenic_segments gene_segments ORF_sequences intergenic_segments_tandem intergenic_segments_convergent intergenic_segments_divergent,upstream,downstream; do \
		${MAKE} one_seq_type_len_distrib SEQ_TYPE=$${type} ; \
		${MAKE} one_seq_type_len_histo SEQ_TYPE=$${type} ; \
	done

################################################################
## Detect ove- and under-represented oligonucleotides

## A single analysis
#SEQ_TYPE=genome
MKV_DIR=results/${ORG}/oligos
OL=6
STR=-2str
MKV=3
NOOV=-noov
MKV_FILE=${MKV_DIR}/${SEQ_TYPE}_${OL}nt${STR}${NOOV}_mkv${MKV}
V=1
OLIGO_CMD=oligo-analysis -v ${V} -i ${SEQ_FILE} -format wc -two_tail \
		-l ${OL} -markov ${MKV} ${STR} ${NOOV} -lth occ_sig 0 \
		-return occ,freq,ratio,proba,rank,zscore,like -sort  \
		-o ${MKV_FILE}
one_seq_type_one_markov_one_str:
	@mkdir -p ${MKV_DIR}
	@echo ""
	@echo ${OLIGO_CMD}
	@${OLIGO_CMD}
	@echo ${MKV_FILE}

## Iterate over strands (single or both)
one_seq_type_one_markov:
#	for str in '-2str' '-1str' ; do 
	for str in '-1str' ; do \
		${MAKE} one_seq_type_one_markov_one_str STR=$${str} ; \
	done

## Iterate over selected oligo lengths and markov orders
one_seq_type_all_markov:
	@${MAKE} one_seq_type_one_markov OL=1 MKV=0
	@${MAKE} one_seq_type_one_markov OL=2 MKV=0
	@${MAKE} one_seq_type_one_markov OL=3 MKV=1
	@${MAKE} one_seq_type_one_markov OL=4 MKV=2
	@${MAKE} one_seq_type_one_markov OL=5 MKV=3

## all markov models for hexamers (just because I like them)
	@${MAKE} one_seq_type_one_markov OL=6 MKV=0
	@${MAKE} one_seq_type_one_markov OL=6 MKV=1
	@${MAKE} one_seq_type_one_markov OL=6 MKV=2
	@${MAKE} one_seq_type_one_markov OL=6 MKV=3
	@${MAKE} one_seq_type_one_markov OL=6 MKV=4

	@${MAKE} one_seq_type_one_markov OL=7 MKV=5
	@${MAKE} one_seq_type_one_markov OL=8 MKV=6
	@${MAKE} one_seq_type_one_markov OL=9 MKV=7
	@${MAKE} one_seq_type_one_markov OL=10 MKV=8


## Iterate over all sequence types
all_seq_all_markov:
	@for type in genome intergenic_segments ORF_sequences; do \
		${MAKE} one_seq_type_all_markov SEQ_TYPE=$${type}; \
	done
