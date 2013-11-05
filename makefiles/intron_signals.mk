################################################################
##
## Predict potential splicing signals by detecting positionally biased words
## in the 3' or 5' ends of first intron sequences (Homo sapiens).

include ${RSAT}/makefiles/util.mk
MAKEFILE=makefile

SEQ_TYPE=firstintron
ORG=Homo_sapiens

################################################################
## Get all intron sequences from Ensembl.  Beware, this can take several hours
## because sequences have to be fetched one by one.
SEQ_DIR=data/sequences
SEQ=${SEQ_DIR}/all_${SEQ_TYPE}_${ORG}${RM}
RM=-rm
get_sequences:
	@echo
	@echo "Retrieving all sequences from Ensembl"
	@mkdir -p ${SEQ_DIR}
	retrieve-ensembl-seq.pl -v 2 \
		-org ${ORG} \
		-all -feattype intron -firstintron \
		-mask coding ${RM} \
		-o ${SEQ}.fasta

################################################################
## Compute sequence length distributions
seq_len_distrib:
	sequence-lengths -i ${SEQ}.fasta \
		| classfreq -v -col 2 -ci 100 \
		-o ${SEQ}_len_distrib.tab
	@echo ${SEQ}_len_distrib.tab
	XYgraph -i ${SEQ}_len_distrib.tab \
		-xcol 3 -ycol 7,8,9 \
		-xsize 800 -ysize 400 \
		-ymin 0 -ymax 1 -legend \
		-xleg1 'first intron size (log scale)' \
		-yleg1 'frequency' \
		-xlog 10 -lines \
		-o ${SEQ}_len_distrib.png
	@echo ${SEQ}_len_distrib.png

################################################################
## Detect positionally biased words
ORI=end
FROM=-200
CI=1
TO=-1
OL=3
POS_DIR=results/${ORG}/${SEQ_TYPE}/${OL}nt
POS=${POS_DIR}/${OL}nt_ci${CI}_all_${SEQ_TYPE}_${ORI}_from${FROM}_to${TO}
POS_CMD=sub-sequence -i ${SEQ}.fasta  -origin ${ORI} \
		-from ${FROM} -to ${TO} \
		| position-analysis -v 2  \
		-origin ${ORI} \
		-1str -noov -sort \
		-return chi,distrib,graph,rank \
		-max_graphs 50 \
		-ci ${CI} -l ${OL} -o ${POS}.tab
one_len:
	@echo
	@echo 
	@mkdir -p ${POS_DIR}
	@echo "${POS_CMD}"
	@${POS_CMD}
	@echo ${POS}.tab
	@text-to-html -i ${POS}.tab -o ${POS}.html
	@echo ${POS}.html

## Iterate over lengths
OLIGO_LENGTHS=1 2 3 4 5 6 7 8
all_len:
	@for l in ${OLIGO_LENGTHS}; do \
		${MAKE} one_len OL=$${l}; \
	done

## Iterate over sides
SIDE_TASK=all_len
all_sides:
	${MAKE} ${SIDE_TASK} ORI=start FROM=1 TO=200
	${MAKE} ${SIDE_TASK} ORI=end FROM=-1 TO=-200
