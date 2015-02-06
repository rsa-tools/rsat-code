################################################################
## Quick tester for the RSAT tool convert-matrix
##
## Authors: Jacques van Helden
## Date: Feb 2015


include ${RSAT}/makefiles/util.mk
MAKEFILE=${RSAT}/makefiles/convert-matrix_test.mk

## Choose logo program (weblogo|seqlogo)
LOGO_PROGRAM=weblogo

#################################################################
## Generate logos from the matrices
PREFIX=Dmelanogaster_segmentation_12matrices
IN_MATRICES=${RSAT}/public_html/demo_files/${PREFIX}
RES_DIR=results/convert-matrix_test/${PREFIX}
LOGO_DIR=${RES_DIR}/${LOGO_PROGRAM}
OUT_MATRICES=${RES_DIR}/${PREFIX}
FROM=tf
TO=tab
LOGO_FORMAT=png
logos:
	@echo
	@echo "Converting ${IN_MATRICES}"
	@echo "	from ${FROM} to ${TO}"
	@echo "	LOGO_PROGRAM	${LOGO_PROGRAM}"
	@echo "	LOGO_FORMAT	${LOGO_FORMAT}"
	@echo "	RES_DIR		${RES_DIR}"
	@echo "	LOGO_DIR	${LOGO_DIR}"
	@mkdir -p ${RES_DIR} ${LOGO_DIR}
	convert-matrix -v ${V} -i ${IN_MATRICES}.${FROM} -from ${FROM} \
		-return counts,parameters,logo -logo_dir ${LOGO_DIR} \
		-logo_program ${LOGO_PROGRAM} \
		-logo_format ${LOGO_FORMAT} \
		-to ${TO} -o ${OUT_MATRICES}.${TO}
	@echo "	${OUT_MATRICES}.${TO}"

## Convert 12 matrices involved in Drosophila segmentation
droso_segm:
	@${MAKE} logos PREFIX=Dmelanogaster_segmentation_12matrices

## Convert LexA matrix (E.coli)
lexa:
	@${MAKE} logos PREFIX=LexA.2nt_upstream-noorf-ovlp-2str.20
