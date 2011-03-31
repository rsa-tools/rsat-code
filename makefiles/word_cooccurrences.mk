################################################################
##
## Count joint distribution of two words in randomly generated sequences
##
## Implementation : Maud Vidick and Jacques van Helden
##
## March 2011.

include ${RSAT}/makefiles/util.mk
MAKEFILE=${RSAT}/makefiles/word_cooccurrences.mk


ORG=Saccharomyces_cerevisiae
#ORG=Mus_musculus_EnsEMBL
L=10000
N=10000
MKV=2
BG=upstream-noorf

## Those parameters are automatically derived from the previous ones
PREFIX=${ORG}_${BG}_m${MKV}_L${L}_N${N}
RES_DIR=results/${PREFIX}

################################################################
## Define the targets


## Run all the required tasks in teh appropriate order
all: randseq count ct heatmap compress
	@echo "Fininshed by job"

compress:
	gzip ${SEQ}

uncompress:
	gunzip ${SEQ}

## Generate a set of random sequences
SEQ=${RES_DIR}/${PREFIX}_randseq.fasta
randseq:
	@echo
	@echo "Generating random sequences"
	@mkdir -p ${RES_DIR}
	random-seq -l ${L} -n ${N} -format fasta \
		-bg ${BG} -org ${ORG} -markov ${MKV} -type dna \
		-o ${SEQ}
	@echo ${SEQ}

## Count word occurrences
STR=-1str
PATTERNS= patterns.txt
COUNTS=${RES_DIR}/${PREFIX}_counts
V=1
count:
	@echo
	@echo "Counting word occurrences"
	dna-pattern -v ${V} -pl ${PATTERNS} -i ${SEQ} -format fasta \
		-ovlp ${STR} -return table -return rowsum -o colsum \
		-o ${COUNTS}.tab
	@echo ${COUNTS}.tab

## Compute the contingency table
ct:
	@echo
	@echo "Computing contingency table"
	contingency-table -v 1 -i ${COUNTS}.tab -col1 3 -col2 4 -density -sort num_labels \
		 -o ${COUNTS}_xtab.tab
	@echo ${COUNTS}_xtab.tab

## Draw the heat map
heatmap:
	@echo
	@echo "Drawing heat map"
	draw-heatmap -rownames -i ${COUNTS}_xtab.tab -o ${COUNTS}_xtab.png
	@echo ${COUNTS}_xtab.png