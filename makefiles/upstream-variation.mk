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

VCF_INPUT = ${VCF_INPUT}    # might be GZIP-compressed
RET_INPUT = ${VCF_INPUT}    # might be GZIP-compressed

# default chr prefixes
ifeq ($(CHR_ORI),)
	CHR_ORI=""
	CHR_NEW=""
endif

# default location and names of temp files
TMPDIR        = /tmp
VCF_POSITIONS = ${TMPDIR}/$(basename VCF_INPUT).variants.tab
RET_RANGES    = ${TMPDIR}/$(basename RET_INPUT).ranges.tab

#this could be used to call this makefile from others such as ensemblgenomes_FTP_client.mk
#SPECIES_DIR=${RSAT}/data/genomes/${SPECIES_RSAT_ID} instead of RSAT_ORG
#GENOME_DIR=${SPECIES_DIR}/genome
#VARIATIONS_DIR=${SPECIES_DIR}/variations
#data/genomes/Oryza_sativa.IRGSP-1.0.38/genome/Oryza_sativa.IRGSP-1.0.38_upstream-noorf.fasta.gz

#### Parses a VCF-file and produces a temporary list of tab-separated polymorphism positions.
#### Optionanly uses $CHR_NEW to replace $CHR_ORI chr names in VCF so that they match those in FASTA file
#### NOTE: only first base of variant is taken
_vcf_pos:
	@echo
	@echo "Collecting polymorphisms positions from ${VCF_INPUT}" 
	@gzip -cdfq ${VCF_INPUT} | perl -lne 's/^${CHR_ORI}/${CHR_NEW}/; print' | \
		perl -F'\t' -lane 'next if(/^#/); print "@F[0..1]"' > ${VCF_POSITIONS}
	@echo "List of positions stored in ${VCF_POSITIONS}"

### Sequence_pos gathers genomic ranges and strand of retrieved sequences
_sequence_pos:
	@echo
	@echo "Parsing sequence range coordinates from ${RET_INPUT}"
	@gzip -cdfq ${RET_INPUT} | perl -lne 'print $$1 if(/location: \S+?:(\S+?);/)' > ${RET_RANGES} 
	@echo "Genomic ranges can be found at ${RET_RANGES}"

_clean:
	@echo
	@echo "Clean temporary files ${VCF_POSITIONS}, ${RET_RANGES}"
	@rm ${VCF_POSITIONS} ${RET_RANGES}

dist: _vcf_pos _sequence_pos
	@echo
	@echo "Computing distribution of polymorphisms in upstream sequences"
	@Rscript ${RSAT}/R-scripts/upstream_SNP_distribution.R ${VCF_POSITIONS} ${RET_RANGES}
	$(MAKE) _clean
	@echo "Done"

VCF_INPUT = ${RSAT}/public_html/demo_files/Ppersica_Varieties.chr1.vcf.gz
RET_INPUT = ${RSAT}/public_html/demo_files/Prunus_persica.Prunus_persica_NCBIv2.38_upstream.chr1.fasta.gz
CHR_ORI=Pp0
CHR_NEW=G
demo:
	$(MAKE) dist
