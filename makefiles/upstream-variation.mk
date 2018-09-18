#############
#
# This makefile takes as input a retrieve-seq file with upstream sequences and a vcf file of the same species.
# Plots the density of polymorphisms (SNPs and indels) overlapping the retrieved sequences.
# Optionally it supports different chromosome names in both input files.
#
# Author: Chesco Montardit Tarda 2018 (edited by BContreras)
#
#############

include ${RSAT}/makefiles/util.mk
MAKEFILE=${RSAT}/makefiles/upstream-variation.mk

## Define parameters


DIR = ${DIR}
VCF_INPUT = ${VCF_INPUT}
RET_INPUT = ${RET_INPUT}
RSAT_ORG = ${RSAT_ORG}
CHR_NAME = ${CHR_NAME}
VCF_POSITIONS = ${RSAT_ORG}_SNPs_position.tab
RET_RANGES = ${RSAT_ORG}_ranges.tab

#### vcf_pos parses a VCF-file to output a list of polymorphisms positions
vcf_pos:
	@echo
	@echo	"Collecting polymorphisms positions from ${VCF_INPUT}" 
	@perl -F'\t' -lane 'print "@F[0..1]"' ${DIR}/${VCF_INPUT} > ${DIR}/${VCF_POSITIONS}
	@perl -i -ne 's/${CHR_NAME}//; print' ${DIR}/${VCF_POSITIONS}
	@sed -i '/scaffold/d' ${DIR}/${VCF_POSITIONS}
	@echo	"List of positions can be found at  ${DIR}/${VCF_POSITIONS}"

### sequence_pos gathers genomic ranges and strand of retrieved sequences
sequence_pos:
	@echo
	@echo "Parsing ${RET_INPUT} to acquire genomic ranges of sequences"	
	@awk -F';' '/location:/{sub(/.*location: /,""); print $1}' ${RET_INPUT} > ${DIR}/${RET_RANGES}
	@perl -i -ne 's/${RSAT_ORG}://, print' ${DIR}/${RET_RANGES}
	@sed -i '/\;/ s/;.*//' ${DIR}/${RET_RANGES}
	@echo "Genomic ranges can be found at ${DIR}/${RET_RANGES}"

### This command needs the next order of variables when called: VCF_POSITIONS=X RET_RANGES=Y DIR=Z etc
snp_vs_sequence:
	@echo
	@echo "Overlapping polymorphisms in upstream sequences"
	@Rscript ${RSAT}/R-scripts/upstream_SNP_distribution.R ${DIR}/${VCF_POSITIONS} ${DIR}/${RET_RANGES}
	@rm tmp_points.tab
	@echo "Done"
