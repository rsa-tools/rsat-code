################################################################
## Some tests with the program matrix-distrib
## Given a PSSM, compare the theoretical distributions with various
## background models


include ${RSAT}/makefiles/util.mk
MAKEFILE=${RSAT}/makefiles/matrix_distrib.mk

all: list_params iid_model iid_distrib all_bgol compare_distrib

################################################################
## List the parameters
FACTOR=PHO4
MATRIX=${FACTOR}_matrix.tab
MATRIX_FORMAT=tab
ORG=Saccharomyces_cerevisiae
BG_MODEL=upstream-noorf
PSEUDO=1
DECIMALS=1
list_params:
	@echo "	FACTOR  	${FACTOR}"
	@echo "	ORG      	${ORG}"
	@echo "	BG_MODEL	${BG_MODEL}"
	@echo "	MATRIX    	${MATRIX}"
	@echo "	MATRIX_FORMAT	${MATRIX_FORMAT}"
	@echo "	PSEUDO    	${PSEUDO}"
	@echo "	DECIMALS    	${DECIMALS}"

################################################################
## Compute the theoretical distribution of scores for one background model
BGOL=1
BG_FILE=${RSAT}/data/genomes/${ORG}/oligo-frequencies/${BGOL}nt_${BG_MODEL}_${ORG}-noov-1str.freq.gz
DISTRIB_FILE=distrib_${BGOL}nt_${BG_MODEL}.tab
one_distrib:
	matrix-distrib -v 1 -m ${MATRIX} -matrix_format ${MATRIX_FORMAT} -pseudo ${PSEUDO} -bgfile ${BG_FILE} -o ${DISTRIB_FILE}
	@echo ${DISTRIB_FILE}

################################################################
## Compute theoretical distribution with markov chains of increasing orders
all_bgol:
	@for ol in 1 2 3 4 5 6 ; do \
		${MAKE} one_distrib BGOL=$${ol}; \
	done

################################################################
## Create a file with IID (independently and identically distributed)
## model
iid_model:
	@echo "A	0.25" > iid_model.tab
	@echo "C	0.25" >> iid_model.tab
	@echo "G	0.25" >> iid_model.tab
	@echo "T	0.25" >> iid_model.tab
	@echo iid_model.tab

################################################################
## Compute theoretical distribution with an IID model
iid_distrib:
	${MAKE} one_distrib BG_MODEL=iid BGOL=1 BG_FILE=iid_model.tab 

################################################################
## Compare the theoretical distributions
IMG_FORMAT=png
compare_distrib:
	compare-scores -sc 4 -numeric -suppress distrib_ -suppress .tab -o distrib_comparison.tab -files distrib_*nt*
	@echo  distrib_comparison.tab
	XYgraph -i distrib_comparison.tab -title '${ORG} ${FACTOR} matrix' \
		-xcol 1 -ycol 2-10 -pointsize 0 \
		-xleg1 'score' -yleg1 'P-value' \
		-legend -lines -xsize 600 -ysize 400 \
		-format ${IMG_FORMAT} -o distrib_comparison.${IMG_FORMAT}
	@echo  distrib_comparison.${IMG_FORMAT}

