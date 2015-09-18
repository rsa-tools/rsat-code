################################################################
## Run infer-operons to detect operons and directons

include ${RSAT}/makefiles/util.mk
MAKEFILE=${RSAT}/makefiles/infer_operons_tester.mk

ORG=Escherichia_coli_K_12_substr__MG1655_uid57779
OPERON_DIR=operons/${ORG}
TYPE=operons
OPERONS=${OPERON_DIR}/${ORG}_${TYPE}_dist${DIST}_mingenes${MIN_GN}
FIELDS=name,operon,query,leader,upstr_dist,gene_nb
DIST=55
MIN_GN=1
infer_operons:
	@echo
	@echo "Inferring operons	${ORG}"
	@mkdir -p ${OPERON_DIR}
	infer-operons -v 1 -dist ${DIST} -min_gene_nb ${MIN_GN} -return ${FIELDS} -org ${ORG} -all -o ${OPERONS}.tab
	@echo "	${OPERONS}.tab"
	@text-to-html -i ${OPERONS}.tab -o ${OPERONS}.html
	@echo "	${OPERONS}.html"

infer_directons:
	${MAKE} infer_operons DIST=1000000 TYPE=directons
