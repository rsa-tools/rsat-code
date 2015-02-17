################################################################
## Demonstration for the too footprint-discovery

include ${RSAT}/makefiles/util.mk
MAKEFILE=${RSAT}/makefiles/footprint-discovery_demo.mk

ORG=Escherichia_coli_K_12_substr__MG1655_uid57779
TAXON=Enterobacteriales
GENE=lexA

TASK=query_seq,filter_dyads,orthologs,ortho_seq,purge,dyads,maps,gene_index,index

## Run footprint-discovery for a single gene of interest.
FP_DISCO_DIR=results/footprint-discovery_demo/${ORG}/${TAXON}/${GENE}
disco:
	@echo
	@echo "Running footprint-discovery 	${GEBE}	${ORG}	${TAXON}"
	footprint-discovery  -v 1 -org ${ORG} -taxon ${TAXON} \
		-q ${GENE} \
		-lth occ 1 \
		-lth occ_sig 0 \
		-uth rank 50 \
		-return occ,proba,rank \
		-filter \
		-bg_model taxfreq \
		-task ${TASK} ${OPT} \
		-o ${FP_DISCO_DIR}
