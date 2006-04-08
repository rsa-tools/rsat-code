## ##############################################################
## This scripts checks the consistency between olligo-analysis and
## dyad-analysis, by comparing the results obtained for
## hexanucleotides (oligo-analysis), and for pairs of trinucleotides
## with spacing 0 (dyad-analysis).
##
## This script has been useful for debugging prooblems with the
## treatment of N in dyad-analysis.

include ${RSAT}/makefiles/util.mk

all: oligos dyads all_compa

V=1 
#SEQ_FILE=PURA_Escherichia_coli_K12_Gammaproteobacteria_up_purged.fasta
SEQ_FILE=ACEB.fasta

################################################################
## Common parmeters for oligo-analysis and dyad-analysis
NOOV=-noov
STR=-2str
RETURN=-return occ,freq,proba,rank,zscore,ratio
ORG=Escherichia_coli_K12
BG=-bg upstream-noorf -org ${ORG}
OPTIONS=-v ${V} -i ${SEQ_FILE} ${STR} ${NOOV} ${BG} ${RETURN} -zeroocc -sort
SUFFIX=${STR}${NOOV}

################################################################
## Run oligo-analysis
OLIGOS=oligos${SUFFIX}.tab
oligos:
	oligo-analysis -l 6 ${OPTIONS} -o ${OLIGOS}
	@echo ${OLIGOS}

################################################################
## Run dyad-analysis
DYADS=dyads_sp${SP}${SUFFIX}.tab
DYADS_NOSPACING=dyads_sp${SP}_${SUFFIX}_nospacing.tab
SP=0-1
dyads:
	dyad-analysis -l 3 -sp ${SP} ${OPTIONS} -o ${DYADS} 
	@echo ${DYADS}
	grep 'n{0}' ${DYADS} \
		| perl -pe 's/n\{0\}//g' \
		> ${DYADS_NOSPACING}
	@echo ${DYADS_NOSPACING}

################################################################
## Compare one column of the oligo-analysis and dyad-analysis result
## files
COMPA=compa${SP}_${SUFFIX}_sc${SC}
SC=3
LOG=
DIFF=diff_${COMPA}
one_compa:
	compare-scores -i ${OLIGOS} -i ${DYADS_NOSPACING} -sc ${SC} -o ${COMPA}.tab
	awk -F'\t' '$$2 != $$3 {print $$0"\t"($$3-$$2)}' ${COMPA}.tab > ${DIFF}
	@echo ${COMPA}.tab
	XYgraph -i ${COMPA}.tab \
		-xcol 2 -ycol 3 \
		-format jpg ${LOG} \
		-title1 '${SEQ_FILE}' \
		-title2 'score column ${SC}' \
		-xleg1 'oligos' \
		-yleg1 'dyads' \
		-o ${COMPA}.jpg
	@echo ${COMPA}.jpg

################################################################
## Compare all the comparable columns between oligo-analysis and dyad-analysis
all_compa:
	${MAKE} one_compa SC=3
	${MAKE} one_compa SC=4
	${MAKE} one_compa SC=5
	${MAKE} one_compa SC=6
	${MAKE} one_compa SC=7 LOG='-xlog 10 -ylog 10'
	${MAKE} one_compa SC=8 LOG='-xlog 10 -ylog 10'
	${MAKE} one_compa SC=9
	${MAKE} one_compa SC=10
	${MAKE} one_compa SC=11
	${MAKE} one_compa SC=13
	${MAKE} one_compa SC=14
