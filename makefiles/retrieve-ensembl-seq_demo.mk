################################################################
## Demos for the tool retrieve-ensembl-seq.pl

include ${RSAT}/makefiles/util.mk
MAKEFILE=makefiles/retrieve-ensembl-seq_demo.mk

RES_DIR=results/retrieve-ensembl-seq_demo
ORG=Homo_sapiens
GENE=ENSG00000139618
TYPE=upstream
FEATTYPE=mRNA
TYPE=upstream
FROM=-2000
TO=-1
SEQ=${RES_DIR}/${GENE}_${TYPE}_${FEATTYPE}${FROM}_${TO}.fasta
one_gene:
	@echo "Testing retrieve-ensembl-seq.pl"
	@mkdir -p ${RES_DIR}
	@echo
	@echo "Retrieving ${TYPE} sequences for ${ORG} gene ${GENE}"
	retrieve-ensembl-seq.pl -v ${V} \
		-org ${ORG} \
		-q ${GENE} \
		-from ${FROM} -to ${TO} \
		-feattype ${FEATTYPE} -type ${TYPE} \
		-lw '60' -alltranscripts -utr 'all' -header_org 'scientific' \
		-o ${SEQ}
	@echo "Sequence file"
	@echo "	${SEQ}"
	@echo "Sequence lengths"
	@sequence-lengths -i ${SEQ} -unit bp
