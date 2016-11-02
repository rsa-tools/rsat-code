################################################################
## Test the level of compression of different genomes.

include ${RSAT}/makefiles/util.mk
MAKEFILE=${RSAT}/makefiles/genome_compression_tests.mk

ORG=Escherichia_coli_K12
GENOME=${RSAT}/data/genomes/${ORG}/genome/contigs.txt

CHROM_SIZES=chrom_sizes_${ORG}.tab
chrom_lengths:
	@echo
	@echo "${DATE}	Measuring chromosome sizes	${ORG}"
	sequence-lengths -i ${GENOME} -in_format filelist -o ${CHROM_SIZES}
	@echo "	${CHROM_SIZES}"

GENOME_SIZE=`convert-seq -i ${GENOME} -from filelist -to fasta | wc -c`
uncompressed_genome_size:
	@echo
	@echo "${DATE}	Measuring uncompressed genome size	${ORG}"
	@echo "	${GENOME}"
	@echo "	GENOME_SIZE	${GENOME_SIZE}"

GZIPPED_GENOME_SIZE=`convert-seq -i ${GENOME} -from filelist -to fasta | gzip --to-stdout | wc -c`
compressed_genome_size:
	@echo
	@echo "${DATE}	Compressing genome	${ORG}"
	@echo "	${GENOME}"
	@echo "	GZIPPED_GENOME_SIZE	${GZIPPED_GENOME_SIZE}"

RAND_EQUI_SIZE=`random-seq -i ${CHROM_SIZES} -template_format len | gzip --to-stdout | wc -c`
rand_equi:

SIZE_TABLE=genome_sizes.tab
init_size_table:
	@echo
	@echo "${DATE}	Initiating genome size table"
	@echo "Organism	genome	gzipped_genome	rand_bernoulli" > ${SIZE_TABLE}
	@echo "	${SIZE_TABLE}"

TAXONOMY=`supported-organisms -return ID,taxonomy | awk -F '\t' '$$1=="${ORG}" {print $$2}'`
one_genome: chrom_lengths
	@echo "${DATE}	Analyzing genome size	${ORG}"
	@echo "${ORG}	${GENOME_SIZE}	${GZIPPED_GENOME_SIZE}	${RAND_EQUI_SIZE}	${TAXONOMY}" >> ${SIZE_TABLE}
	@tail -n 1 ${SIZE_TABLE}
	@echo "${DATE}	Done	${ORG}"
	@echo

ORGANISMS=Escherichia_coli_K12 \
	Saccharomyces_cerevisiae \
	Plasmodium_falciparum \
	Trypanosoma_brucei_TryBru_Apr2005_chr11 \
	Arabidopsis_thaliana_TAIR10.29 \
	Drosophila_melanogaster \
	Homo_sapiens_GRCh37 \
	Mus_musculus_GRCm38

#	${PLANTS} ${FUNGI} ${BACTERIA} ${ARCHAEA} ${METAZOA}
PLANTS=`supported-organisms  -unique_genus -taxon Viridiplantae`
FUNGI=`supported-organisms -unique_genus -taxon Fungi`
BACTERIA=`supported-organisms -unique_genus -taxon Bacteria`
ARCHAEA=`supported-organisms -unique_genus -taxon Archaea`
METAZOA=`supported-organisms -unique_genus -taxon Archaea`
some_genomes:
	@echo "ORGANISMS	${ORGANISMS}"
	@for org in ${ORGANISMS}; do \
		${MAKE} one_genome ORG=$${org}; \
	done

TAXON=Viridiplantae
one_taxon:
	@${MAKE} some_genomes ORGANISMS="`supported-organisms -unique_genus -taxon ${TAXON} | xargs`" SIZE_TABLE=genome_sizes_${TAXON}.tab


compression_plot:
	cat genome_sizes*.tab | grep -v 'genome'| \
		awk -F '\t' '$$4 >0 {print $$2"\t"$$3"\t"$$4"\t"$$4/$$3"\t"$$1"\t"$$5}'|\
		sort -rn -k 3 > genome_compression_stats.tab
	@echo "	genome_compression_stats.tab"
	XYgraph -v ${V} -i genome_compression_stats.tab \
		-xcol 1 -ycol 4 -title "Genome compression levels (gzip)" \
		-xleg1 "Uncompressed genome size" -yleg1 "Compressed genome vs random seq" \
		-xlog -hline red 1 -ygstep 0.1 ${OPT} \
		-o genome_compression_stats.png
	@echo "	genome_compression_stats.png"

