################################################################
## Compare all operons with all gene ontology classes


include ${RSAT}/makefiles/util.mk
MAKEFILE=makefiles/operons_vs_go.mk

## Define parameters
ORG=Escherichia_coli_K12
MIN_GENES=2
OPERONS=results/GO/operons_Escherichia_coli_K12.tab
GENE_OPERON=results/GO/gene_operon_${ORG}.tab
GENE_GO=results/GO/gene_GO_name_ids_${ORG}.tab
GO_DESC=results/GO/GO_description.tab
OPERONS_VS_GO=results/GO/operons_vs_GO_${ORG}_QR${MIN_GENES}

operons:
	@echo
	@echo "Inferring operons"
	infer-operon -v 1 -dist 55 -min_gene_nb 2 \
		-return query,name,leader,operon,upstr_dist,q_info,gene_nb \
		-all  -org ${ORG} -o ${OPERONS}
	@echo "	${OPERONS}"

## Operons versus GO annotations
compa:
	@echo
	@echo "Comparing operons with GO annotations"
	compare-classes  -v 2 \
		-q ${GENE_OPERON} -r ${GENE_GO} \
		-lth Q ${MIN_GENES} -lth R ${MIN_GENES} -lth QR ${MIN_GENES} \
		-return occ,freq,jac_sim,proba,rank,common,Q_only \
		-dot ${OPERONS_VS_GO}.dot \
		-gml ${OPERONS_VS_GO}.gml \
		-rnames ${GO_DESC} \
		-o ${OPERONS_VS_GO}.tab
	@echo "	${OPERONS_VS_GO}.tab"
	@echo "	${OPERONS_VS_GO}.gml"
	@echo "	${OPERONS_VS_GO}.dot"


it:
	@text-to-html -i ${OPERONS_VS_GO}.tab -o ${OPERONS_VS_GO}.html
	@echo "	${OPERONS_VS_GO}.html"
