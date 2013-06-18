################################################################
## Install a genome from a gff file (specification of features) + a
## fasta file (contig sequences).

include ${RSAT}/makefiles/util.mk
MAKEFILE=makefiles/install_genome_from_gff.mk

## Input parameters
ORG=Acanthamoeba_castellanii
TAXID=5755
TAXONOMY=Eukaryota; Amoebozoa; Centramoebida; Acanthamoebidae; Acanthamoeba
DIR_SRC=${RSAT}/temp_data/new
CONTIGS_FASTA=${DIR_SRC}/acast_assembly_v5.fasta
ANNOT_GFF=${DIR_SRC}/acast_assembly_v5.genes.annotated.gff

## Output parameters
GENOME_DIR=${RSAT}/public_html/data/genomes/${ORG}/genome

################################################################
## Convert th fasta sequences into a series of raw sequence files + a
## file containing the list of raw files.
CONTIGS_LIST=${GENOME_DIR}/contigs.txt
format_contig_seq:
	@echo
	@echo "Formatting and installing contig sequences"
	@echo "Genome directory	${GENOME_DIR}"
	@mkdir -p ${GENOME_DIR}
	convert-seq -i ${CONTIGS_FASTA} -from fasta -to filelist -o ${CONTIGS_LIST}

################################################################
## Convert features from gff to RSAT-compatible genomic features
GENOME_FEATURES=${GENOME_DIR}/feature.tab
NB_FEATURES_GFF=`wc -l ${ANNOT_GFF} | awk '{print $$1}'`
NB_FEATURES_PARSED=`grep -v '^;' ${GENOME_FEATURES} | grep -v '^\#' | wc -l`
convert_gff:
	@echo
	@echo "Converting gff to RSAT genome features	${ANNOT_GFF}"
	@echo "GFF file contains ${NB_FEATURES_GFF} lines"
#	convert-features -from gff -to gft  -i ${ANNOT_GFF} -o ${GENOME_FEATURES}
	@echo "Genome feature file	${GENOME_FEATURES}"
	@echo "Parsed ${NB_FEATURES_PARSED} features"


################################################################
## Extract one fetaure type from the complete feature file
FEATTYPE=gene
GENOME_FEATTYPE=`echo ${GENOME_DIR}/${FEATTYPE}.tab | tr 'A-Z' 'a-z'`
extract_one_feattype:
	@echo
	@echo "Extracting ${FEATTYPE} from feature file	${GENOME_FEATURES}"
	@grep '^;' ${GENOME_FEATURES} > ${GENOME_FEATTYPE}
	@awk -F'\t' '$$2 == "${FEATTYPE}"' ${GENOME_FEATURES} >> ${GENOME_FEATTYPE}
	@echo "${FEATTYPE} file	${GENOME_FEATTYPE}"
	@wc -l ${GENOME_FEATTYPE}

extract_cds:
	${MAKE} extract_one_feattype FEATTYPE=CDS

extract_genes:
	${MAKE} extract_one_feattype FEATTYPE=gene

################################################################
## Add the genome to the list of genomes supported in RSAT
config:
	@echo "${TAXID}	${TAXONOMY}" > ${GENOME_DIR}/organism.tab
	@echo "Created organism file	${GENOME_DIR}/organism.tab"
	install-organism -org ${ORG} -task config


