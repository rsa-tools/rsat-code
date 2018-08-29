################################################################
## This makefile contains some targets to download genome seqs and
## annotations from ensemblgenome FTP site, parse and install them on
## an RSAT server.
##
## Authors:
##   Bruno Contreras Moreira <bcontreras@eead.csic.es>
##   Jacques van Helden <Jacques.van-Helden@univ-amu.fr>

## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
## Currently does not work for Bacteria besides Ecoli K12 and for most
## Fungi, as these are further grouped by collections, and hence
## stored in bacteria_NN_collection subfolders. However, the
## species_EnsemblBacteria.txt lists all available genomes and the
## collection they belong to. (BCM)
## 
## Does not work for the main ensembl ftp site, because there is no
## organism table as in ensemblgenomes. I (JvH) need to contact Dan
## Staines (ensemblgenomes) to see how we can manage this.
## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

include ${RSAT}/makefiles/util.mk
MAKEFILE=${RSAT}/makefiles/ensemblgenomes_FTP_client.mk

## Define parameters
V=2

################################################################
## The variable GROUP is defined in RSAT_config.mk, since it depends
## on the server; It can however be overwritten on the command
## line when calling the makefile , for example 
##   make -f makefiles/ensemblgenomes_FTP_client.mk GROUP=Plants
GROUP_LC=$(shell echo $(GROUP) | tr A-Z a-z)
ifeq ($(GROUP),Fungi)
  ## For Fungal genomes, we have to extract the collection (sub-folder
  ## of the ftp site) from the organism table.
  COLLECTION_FROM_TABLE=`awk -F'\t' '$$2 =="${SPECIES}" {print $$13}' ${ORGANISM_TABLE} | perl -pe 's|_core_${RELEASE}_.*||'`

  ## Trick: some Fungi are in the root folder rather than in a
  ## collection sub-folder.  These species however have a collection
  ## field in the species table, but its value equals the species
  ## name.
  ifeq (${COLLECTION_FROM_TABLE}, ${SPECIES})
    COLLECTION=TRICK
  else 
    COLLECTION=${COLLECTION_FROM_TABLE}
  endif

else ifeq ($(GROUP),Bacteria)
  ## For Fungal genomes, we have to extract the collection (sub-folder
  ## of the ftp site) from the organism table.
  COLLECTION_FROM_TABLE=`awk -F'\t' '$$2 =="${SPECIES}" {print $$13}' ${ORGANISM_TABLE} | perl -pe 's|_core_${RELEASE}_.*||'`
  COLLECTION=${COLLECTION_FROM_TABLE}
else
  COLLECTION=
endif


RELEASE=${ENSEMBLGENOMES_RELEASE}
# should be set in RSAT_config.props
SERVER_URL=ftp://ftp.ensemblgenomes.org/pub/${GROUP_LC}
DATABASE=${SERVER_URL}/release-${RELEASE}
# as of june2018 this works, but default seems to be:
#SERVER_URL=ftp://ftp.ensemblgenomes.org/pub/release-${RELEASE}
#DATABASE=${SERVER_URL}/${GROUP_LC}


#name preffix hard-coded, might change in future
SERVERLIST=${DATABASE}/species_Ensembl${GROUP}.txt

ORGANISM_DIR=${RSAT}/data/ensemblgenomes/${GROUP_LC}/release-${RELEASE}
ORGANISM_TABLE=${ORGANISM_DIR}/species_Ensembl${GROUP}.txt
## SPECIES=arabidopsis_thaliana ## The default species is coupled to the group specificity of the server -> I (JvH) move it to RSAT_config.mk

## Note (JvH 2015-11-06) I change SPECIES DIR to directly download
## fasta and gtf in the genome dir, since we will use it for vairous
## purposes.
SPECIES_UCFIRST=$(shell perl -e 'print ucfirst ${SPECIES}')
#SPECIES_ID=${SPECIES_UCFIRST}.${ASSEMBLY_ID}
SPECIES_ID=${SPECIES_UCFIRST}.${GCA_ID}
SPECIES_RSAT_ID=${SPECIES_UCFIRST}.${ASSEMBLY_ID}.${RELEASE}
# SPECIES_DIR=${ORGANISM_DIR}/${SPECIES}
SPECIES_DIR=${RSAT}/data/genomes/${SPECIES_RSAT_ID}
GENOME_DIR=${SPECIES_DIR}/genome
VARIATIONS_DIR=${SPECIES_DIR}/variations

###############################################################
## Get all supported organisms in an ensemblgenome release and store them 
organisms:
	@echo
	@mkdir -p ${ORGANISM_DIR}
	@echo "Getting list of organisms from ${DATABASE}"
	@echo "	${ORGANISM_DIR}"
	@wget -Ncnv ${SERVERLIST} -P ${ORGANISM_DIR}
	@echo
	@echo "	${ORGANISM_TABLE}"

list_param:
	@echo
	@echo "Parameters"
	@echo "	GROUP   		${GROUP} (${GROUP_LC})"
	@echo "	SPECIES			${SPECIES}"
	@echo "	SPECIES_UCFIRST		${SPECIES_UCFIRST}"
	@echo "	TAXON_ID 		${TAXON_ID}"
	@echo "	ASSEMBLY_ID 		${ASSEMBLY_ID}"
	@echo "	RELEASE 		${RELEASE}"
	@echo " SPECIES_ID		${SPECIES_ID}"
	@echo "	SPECIES_RSAT_ID		${SPECIES_RSAT_ID}"
	@echo "Files to download"
	@echo "	DOWNLOAD_TASKS		${DOWNLOAD_TASKS}"
	@echo "	COLLECTION_FROM_TABLE	${COLLECTION_FROM_TABLE}"
	@echo "	COLLECTION		${COLLECTION}"
	@echo "	GTF_FTP_URL		${GTF_FTP_URL}"
	@echo "	FASTA_RAW_FTP_URL	${FASTA_RAW_FTP_URL}"
	@echo "	FASTA_MSK_FTP_URL	${FASTA_MSK_FTP_URL}"
	@echo "	FASTA_PEP_FTP_URL	${FASTA_PEP_FTP_URL}"
	@echo "	SERVER_COMPARA_FILE	${SERVER_COMPARA_FILE}"
	@echo "LOCAL_FILES"
	@echo "	ORGANISM_DIR		${ORGANISM_DIR}"
	@echo "	ORGANISM_TABLE		${ORGANISM_TABLE}"
	@echo "	ALL_SPECIES_NB		${ALL_SPECIES_NB}"
	@echo "	SPECIES_DIR		${SPECIES_DIR}"
	@echo "	GENOME_DIR		${GENOME_DIR}"
	@echo "	PARSE_DIR		${PARSE_DIR}"
	@echo "	GO_DIR			${GO_DIR}"
	@echo "	GTF_LOCAL		${GTF_LOCAL}"
	@echo "	FASTA_RAW_LOCAL		${FASTA_RAW_LOCAL}"
	@echo "	FASTA_MSK_LOCAL		${FASTA_MSK_LOCAL}"
	@echo "	FASTA_PEP_LOCAL		${FASTA_PEP_LOCAL}"
	@echo "	GTF_LOCAL_GZ		${GTF_LOCAL_GZ}"
	@echo "	FASTA_RAW_LOCAL_GZ	${FASTA_RAW_LOCAL_GZ}"
	@echo "	FASTA_MSK_LOCAL_GZ	${FASTA_MSK_LOCAL_GZ}"
	@echo "	FASTA_PEP_LOCAL_GZ	${FASTA_PEP_LOCAL_GZ}"
	@echo "	CMP_GZ			${CMP_GZ}"

################################################################
## Download required files for all organisms
ALL_SPECIES=$(shell cut -f 2 ${ORGANISM_TABLE} | grep -v '^species')
ALL_SPECIES_NB=$(shell cut -f 2 ${ORGANISM_TABLE} | grep -v '^species' | wc -l)
list_all_species:
	@echo 
	@echo "All species from ${ORGANISM_TABLE}"
	@echo ${ALL_SPECIES} | perl -pe 's/\s+/\n/g' |add-linenb -before


ORG_TASKS=organisms
DOWNLOAD_TASKS=download_gtf download_fasta gunzip_downloads 
#INSTALL_TASKS=install_from_gtf index_fasta_downloads install_go_annotations
INSTALL_TASKS=install_from_gtf index_fasta_downloads
COMPARA_TASKS=organisms download_compara parse_compara install_compara
VARIANT_TASKS=download_vcf
# not used
#ALL_TASKS=${ORG_TASKS} ${DOWNLOAD_TASKS} ${INSTALL_TASKS} ${COMPARA_TASKS}

download_one_species: ${DOWNLOAD_TASKS}

install_one_species: ${INSTALL_TASKS}

variations_one_species: ${VARIANT_TASKS}

download_all_species: organisms
	@echo
	@echo Downloading all species in GROUP=${GROUP} RELEASE=${RELEASE}
	@for org in $(ALL_SPECIES); do \
		$(MAKE) download_one_species SPECIES=$$org; \
	done
	@${MAKE} download_go

################################################################
## Install files required for all organisms
install_all_species:
	@echo WARNING: Make sure you run organisms before install_all_species
	@echo
	@echo Installing all species in GROUP=${GROUP} RELEASE=${RELEASE}
	for org in $(ALL_SPECIES); do \
		$(MAKE) install_one_species SPECIES=$$org; \
	done

################################################################
## Download and install compara. This should be done once (concerns
## all species). In addition, it is not supported for all the groups
## -> should remain a separate task.
download_and_install_compara:
	@${MAKE} ${COMPARA_TASKS}


################################################################
## Calculate descriptive stats of installed genomes
STATSDIR=${RSAT}/data/stats/
calc_stats:
	@echo Calculating stats of installed organisms
	@mkdir -p ${STATSDIR}
	@supported-organisms-plots -o data/stats/ -ref thaliana

################################################################
## Check upstrean sequences of all installed species
check_all_species:
	@echo WARNING: Make sure you install_all_species before check_all_species
	@echo
	@echo Checking upstream sequences of all species in GROUP=${GROUP} RELEASE=${RELEASE}
	for org in $(ALL_SPECIES); do \
		$(MAKE) check_sequences SPECIES=$$org; \
	done

################################################################
## Download GTF files from ensemblgenomes
GTF_FTP_URL=${DATABASE}/gtf/${COLLECTION}/${SPECIES}/*${RELEASE}.gtf.gz
download_gtf:
	@echo
	@mkdir -p ${GENOME_DIR}	
	@echo "Downloading GTF file of ${SPECIES}"
	@echo "	GENOME_DIR	${GENOME_DIR}"
	@echo "	GTF_FTP_URL	${GTF_FTP_URL}"
#	@wget -Ncnv ${GTF_FTP_URL} -P ${GENOME_DIR}
#	@wget -cnv ${GTF_FTP_URL} -O ${GTF_LOCAL_GZ}
	@if test -s ${GTF_LOCAL}; then \
		echo "	Uncompressed file exists; skipping	${GTF_LOCAL}"; \
	else \
		wget -cnv ${GTF_FTP_URL} -O ${GTF_LOCAL_GZ} ; \
		echo; \
		ls -1 ${GENOME_DIR}/*.gtf.gz; \
	fi

################################################################
## Download VCF files from ensemblgenomes
VCF_FTP_URL=${DATABASE}/vcf/${COLLECTION}/${SPECIES}/
VCF_SERVER_GZ=${VCF_FTP_URL}/${SPECIES}.vcf.gz
VCF_SERVER_TBI=${VCF_FTP_URL}/${SPECIES}.vcf.gz.tbi
VCF_SERVER_README=${VCF_FTP_URL}/README
VCF_LOCAL_GZ=${VARIATIONS_DIR}/${SPECIES_RSAT_ID}.vcf.gz
VCF_LOCAL_TBI=${VARIATIONS_DIR}/${SPECIES_RSAT_ID}.vcf.gz.tbi
VCF_LOCAL_README=${VARIATIONS_DIR}/README
download_vcf:
	@echo
	@mkdir -p ${VARIATIONS_DIR}
	@echo " VARIATIONS_DIR  ${VARIATIONS_DIR}"	
	@echo
	@echo "Downloading VCF file for species ${SPECIES}"
	wget -cnv ${VCF_SERVER_GZ} -O ${VCF_LOCAL_GZ}; \
	echo "  VCF_LOCAL_GZ  ${VCF_LOCAL_GZ}"; \
	wget -cnv ${VCF_SERVER_TBI} -O ${VCF_LOCAL_TBI}; \
    echo "  VCF_LOCAL_TBI  ${VCF_LOCAL_TBI}"; \
	wget -cnv ${VCF_SERVER_README} -O ${VCF_LOCAL_README};


################################################################
## Download FASTA files with genomic sequences (raw and masked)
## and peptidic sequences
# release <32
#FASTA_RAW_SUFFIX=*${RELEASE}.dna.genome.fa
#FASTA_RAW_FTP_URL=${DATABASE}/fasta/${COLLECTION}/${SPECIES}/dna/${FASTA_RAW_SUFFIX}.gz
#FASTA_MSK_SUFFIX=*${RELEASE}.dna_rm.genome.fa
#FASTA_MSK_FTP_URL=${DATABASE}/fasta/${COLLECTION}/${SPECIES}/dna/${FASTA_MSK_SUFFIX}.gz
#FASTA_PEP_SUFFIX=*${RELEASE}.pep.all.fa
#FASTA_PEP_FTP_URL=${DATABASE}/fasta/${COLLECTION}/${SPECIES}/pep/${FASTA_PEP_SUFFIX}.gz
# releases >32
FTP_SPECIES_PREFIX=${SPECIES_UCFIRST}.*
#FTP_SPECIES_PREFIX=${SPECIES_ID}
FASTA_RAW_SUFFIX=${FTP_SPECIES_PREFIX}.dna.toplevel.fa
FASTA_RAW_FTP_URL=${DATABASE}/fasta/${COLLECTION}/${SPECIES}/dna/${FASTA_RAW_SUFFIX}.gz
FASTA_MSK_SUFFIX=${FTP_SPECIES_PREFIX}.dna_rm.toplevel.fa
FASTA_MSK_FTP_URL=${DATABASE}/fasta/${COLLECTION}/${SPECIES}/dna/${FASTA_MSK_SUFFIX}.gz
FASTA_PEP_SUFFIX=${FTP_SPECIES_PREFIX}.pep.all.fa
FASTA_PEP_FTP_URL=${DATABASE}/fasta/${COLLECTION}/${SPECIES}/pep/${FASTA_PEP_SUFFIX}.gz


## Define local files corresponding to the FTP-downloaded files.
## Note that only the first gtf file is considered
#FASTA_RAW_LOCAL=`ls -1 ${GENOME_DIR}/${FASTA_RAW_SUFFIX} | grep -v '.gz$$'| head -1`
#FASTA_RAW_LOCAL_GZ=`ls -1 ${GENOME_DIR}/${FASTA_RAW_SUFFIX}.gz | head -1`
# _OLD is needed for backwards compatibility
FASTA_RAW_LOCAL_OLD=${GENOME_DIR}/${SPECIES_RSAT_ID}.dna.genome.fa
FASTA_RAW_LOCAL=${GENOME_DIR}/${SPECIES_RSAT_ID}.dna.toplevel.fa
FASTA_RAW_LOCAL_GZ=${FASTA_RAW_LOCAL}.gz
#FASTA_MSK_LOCAL=${GENOME_DIR}/${SPECIES_RSAT_ID}.dna_rm.genome.fa
FASTA_MSK_LOCAL=${GENOME_DIR}/${SPECIES_RSAT_ID}.dna_rm.genome.fa
FASTA_MSK_LOCAL_GZ=${FASTA_MSK_LOCAL}.gz
FASTA_PEP_LOCAL=${GENOME_DIR}/${SPECIES_RSAT_ID}.pep.all.fa
FASTA_PEP_LOCAL_GZ=${FASTA_PEP_LOCAL}.gz
GTF_LOCAL=${GENOME_DIR}/${SPECIES_RSAT_ID}.gtf
GTF_LOCAL_GZ=${GTF_LOCAL}.gz
#FASTA_MSK_LOCAL=`ls -1 ${GENOME_DIR}/${FASTA_MSK_SUFFIX} | grep -v '.gz$$' | head -1`
#FASTA_PEP_LOCAL_GZ=`ls -1 ${GENOME_DIR}/${FASTA_PEP_SUFFIX}.gz | head -1`
#FASTA_PEP_LOCAL=`ls -1 ${GENOME_DIR}/${FASTA_PEP_SUFFIX} | grep -v '.gz$$' | head -1`
#GTF_LOCAL_GZ=$(shell ls -1 ${GENOME_DIR}/*${RELEASE}.gtf.gz | head -1)
#GTF_LOCAL=$(shell ls -1 ${GENOME_DIR}/*${RELEASE}.gtf | head -1)

download_fasta:
	@echo
	@mkdir -p ${GENOME_DIR}
	@echo "	GENOME_DIR	${GENOME_DIR}"
	@echo
	@echo "Downloading raw FASTA genome for species ${SPECIES}"
	@if test -s ${FASTA_RAW_LOCAL}; then \
		echo "	Uncompressed file exists; skipping	${FASTA_RAW_LOCAL}"; \
	else \
		wget -cnv ${FASTA_RAW_FTP_URL} -O ${FASTA_RAW_LOCAL_GZ}; \
		echo "	FASTA_RAW_LOCAL_GZ	${FASTA_RAW_LOCAL_GZ}"; \
		echo "Removing previous fasta index file (.fai)"; \
		echo "	${FASTA_RAW_LOCAL}.fai"; \
		rm -f ${FASTA_RAW_LOCAL}.fai; \
	fi
	@echo
	@echo "Downloading repeat-masked FASTA genome for species ${SPECIES}"
	@if test -s ${FASTA_MSK_LOCAL}; then \
		echo "	Uncompressed file exists; skipping	${FASTA_MSK_LOCAL}"; \
	else \
		wget -cnv ${FASTA_MSK_FTP_URL} -O ${FASTA_MSK_LOCAL_GZ}; \
		echo "	FASTA_MSK_LOCAL_GZ	${FASTA_MSK_LOCAL_GZ}"; \
		echo "Removing previous fasta index file (.fai)"; \
		echo "	${FASTA_MSK_LOCAL}.fai"; \
		rm -f ${FASTA_MSK_LOCAL}.fai; \
	fi
	@echo
	@echo "Downloading FASTA peptidic sequences for species ${SPECIES}"
	@if test -s ${FASTA_PEP_LOCAL}; then \
		echo "	Uncompressed file exists; skipping	${FASTA_PEP_LOCAL}"; \
	else \
		wget -cnv ${FASTA_PEP_FTP_URL} -O ${FASTA_PEP_LOCAL_GZ}; \
		echo "	FASTA_PEP_LOCAL_GZ	${FASTA_PEP_LOCAL_GZ}"; \
		echo "Removing previous fasta index file (.fai)"; \
		echo "	${FASTA_PEP_LOCAL}.fai"; \
		rm -f ${FASTA_PEP_LOCAL}.fai; \
	fi


################################################################
## Download a sample of eg upstream sequences to check installed sequences
check_sequences:
	@echo
	@check-retrieve-seq-rest -v ${V} \
        -org ${SPECIES} \

#################################################################
## Download group COMPARA files from eg
#SERVER_COMPARA_FILE=${DATABASE}/tsv/ensembl-compara/Compara.homologies.${RELEASE}.tsv.gz
SERVER_COMPARA_FILE=${DATABASE}/tsv/ensembl-compara/homologies/Compara.${ENSEMBL_RELEASE}.protein_default.homologies.tsv.gz
download_compara:
	@echo
	@mkdir -p ${ORGANISM_DIR}
	@echo "Downloading COMPARA file of ${GROUP}"
	@echo "	${SERVER_COMPARA_FILE}"
	@wget -Ncnv ${SERVER_COMPARA_FILE} -P ${ORGANISM_DIR}
	@echo
	@ls -1 ${ORGANISM_DIR}/Compara.${ENSEMBL_RELEASE}.protein_default.homologies.tsv.gz

##################################################################
## Download GO ontology file and parse it for server use
GO_DIR=${RSAT}/data/genomes/GO
download_go:
	@echo
	@mkdir -p ${GO_DIR}
	@echo "Downloading and parsing Gene Ontology"
	@make -f ${RSAT}/makefiles/go_analysis.mk GO_DIR=${GO_DIR} download_go parse_go

#################################################################
## Get & install GO annotations for a given species
GO_ANNOT_DIR=${RSAT}/data/genomes/${SPECIES_RSAT_ID}/
GO_ANNOT_FILE=${GO_ANNOT_DIR}/go_annotations.tsv
GO_ANNOT_LINK=go_annotations.tsv
GO_DESC=${GO_DIR}/GO_description.tab
GO_REL=${GO_DIR}/GO_relations.tab
GO_EXPANDED_FILE=expanded_${GO_ANNOT_LINK}
install_go_annotations:
	@echo
	@mkdir -p ${GO_ANNOT_DIR}
	@echo "Downloading GO annotations of ${SPECIES}"
	@echo download-ensembl-go-annotations-biomart -o ${GO_ANNOT_FILE} -org ${SPECIES} \
        -release ${RELEASE} -list ${ORGANISM_TABLE} 
	@download-ensembl-go-annotations-biomart -o ${GO_ANNOT_FILE} -org ${SPECIES} \
		-release ${RELEASE}	-list ${ORGANISM_TABLE}
	@echo "Expanding GO annotations of ${SPECIES}" 
	@rm -f ${GO_ANNOT_LINK}
	@ln -s ${GO_ANNOT_FILE} ${GO_ANNOT_LINK}
	@python2.7 ${RSAT}/python-scripts/go_analysis.py expand -a ${GO_ANNOT_LINK} -d ${GO_DESC} -r ${GO_REL} 
	@mv ${GO_EXPANDED_FILE} ${GO_ANNOT_DIR}
	@rm -f ${GO_ANNOT_LINK}



##################################################################
## Parse GTF file to extract gene, transcripts and cds coordinates.
##
## Beware: this task can be sent to the job queue (qsub) by adding the
## options WHEN=queue to the make command. This is convenient for the
## installation of all genomes.
##
## Examples
##   make -f ${RSAT}/makefiles/ensemblgenomes_FTP_client.mk \
##      install_all_species WHEN=queue
##
## Each species installation will be executed as a job for the
## cluster

## Attempt to get TAXON_ID from $ORGANISM_TABLE
#ifeq (${TAXON_ID},)
TAXON_ID=$(shell grep -w ${SPECIES} ${ORGANISM_TABLE} | cut -f 4)
#endif

## The Assembly ID is important for some model organisms (the
## community relies on some particular assemblies) but is sometimes
## not defined in the table.
ASSEMBLY_ID=$(shell grep -w ${SPECIES} ${ORGANISM_TABLE} | cut -f 5 | perl -lne 's/\s+//g; print')

## The GCA ID is now (2017) recognized by NCBI as well as
## EnsemblGenomes. We should decide if we include it in
## SPECIES_RSAT_ID. 
GCA_ID=$(shell grep -w ${SPECIES} ${ORGANISM_TABLE} | cut -f 6)
PARSE_DIR=${GENOME_DIR}
PARSE_TASK=all
GTF_SOURCE=ensemblgenomes
# Set the following option to -batch in order to dispatch the computation of oligo and dyad frequencies to the job scheduler
PARSE_GTF_OPT=
PARSE_GTF_CMD=parse-gtf -v ${V} -i ${GTF_LOCAL} \
		-fasta ${FASTA_RAW_LOCAL} \
		-fasta_rm ${FASTA_MSK_LOCAL} \
		-fasta_pep ${FASTA_PEP_LOCAL} \
		-org_name ${SPECIES_RSAT_ID} \
		-task ${PARSE_TASK} ${OPT} \
		-taxid ${TAXON_ID} \
		-gtf_source ${GTF_SOURCE} \
		${PARSE_GTF_OPT} -o ${PARSE_DIR} 
parse_gtf:
	@echo
	@echo "Parsing GTF file	${GTF_LOCAL}"
	@echo "TaxonID = ${TAXON_ID}"
	@${MAKE} my_command MY_COMMAND="${PARSE_GTF_CMD}"
	@echo "	${PARSE_DIR}"

###############################################################
## parse gtf and then install organism
# ${FASTA_RAW_LOCAL} is symb linked to ${FASTA_RAW_LOCAL_OLD} for compatibility
install_from_gtf:
	@echo
	@echo "Parsing and installing in RSAT	${SPECIES}"
	@ln -sf ${FASTA_RAW_LOCAL} ${FASTA_RAW_LOCAL_OLD}
	@${MAKE} parse_gtf PARSE_DIR=${RSAT}/public_html/data/genomes/${SPECIES_RSAT_ID}/genome

## Run some test for the GTF parsing result
parse_gtf_test:
	retrieve-seq -org ${SPECIES_RSAT_ID} -from 0 -to 3 -feattype gene | oligo-analysis -v 1 -l 3 -return occ,freq -sort 

###############################################################
## Uncompress GTF and genomic fasta files for beedtools
gunzip_downloads:
	@echo
	@echo "Uncompressing the downloaded GTF and fasta files."
	@if test -s ${GTF_LOCAL_GZ}; then echo "	${GTF_LOCAL_GZ}"; gunzip -qf ${GTF_LOCAL_GZ}; else echo "	skipping GTF_LOCAL_GZ ${GTF_LOCAL_GZ}"; fi;
	@if test -s ${FASTA_RAW_LOCAL_GZ}; then echo "	${FASTA_RAW_LOCAL_GZ}"; gunzip -qf ${FASTA_RAW_LOCAL_GZ}; else echo "	skipping FASTA_RAW_LOCAL_GZ ${FASTA_RAW_LOCAL_GZ}"; fi;
	@if test -s ${FASTA_MSK_LOCAL_GZ}; then echo "	${FASTA_MSK_LOCAL_GZ}"; gunzip -qf ${FASTA_MSK_LOCAL_GZ}; else echo "	skipping FASTA_MSK_LOCAL_GZ ${FASTA_MSK_LOCAL_GZ}"; fi;
	@if test -s ${FASTA_PEP_LOCAL_GZ}; then echo "	${FASTA_PEP_LOCAL_GZ}"; gunzip -qf ${FASTA_PEP_LOCAL_GZ}; else echo "	skipping FASTA_PEP_LOCAL_GZ ${FASTA_PEP_LOCAL_GZ}"; fi;

###############################################################
## Index fasta files for beedtools
index_fasta_downloads:
	@echo
	@echo "Indexing the downloaded fasta files."
	@if test -s ${FASTA_RAW_LOCAL}; then echo "	${FASTA_RAW_LOCAL}.fai"; samtools faidx ${FASTA_RAW_LOCAL}; else echo "	missing FASTA_RAW_LOCAL ${FASTA_RAW_LOCAL}"; fi;
	@if test -s ${FASTA_MSK_LOCAL}; then echo "	${FASTA_MSK_LOCAL}.fai"; samtools faidx ${FASTA_MSK_LOCAL}; else echo "	missing FASTA_MSK_LOCAL ${FASTA_MSK_LOCAL}"; fi;
	@if test -s ${FASTA_PEP_LOCAL}; then echo "	${FASTA_PEP_LOCAL}.fai"; samtools faidx ${FASTA_PEP_LOCAL}; else echo "	missing FASTA_PEP_LOCAL ${FASTA_PEP_LOCAL}"; fi;

###############################################################
## (Re)compress GTF and genomic fasta files for beedtools
gzip_downloads:
	@echo
	@echo "(Re)compressing the downloaded GTF and fasta files."
	@if test -s ${GTF_LOCAL}; then echo "	${GTF_LOCAL}"; gzip -f ${GTF_LOCAL}; else echo "	skipping GTF_LOCAL ${GTF_LOCAL}"; fi;
	@if test -s ${FASTA_RAW_LOCAL}; then echo "	${FASTA_RAW_LOCAL}"; gzip -f ${FASTA_RAW_LOCAL}; else echo "	skipping FASTA_RAW_LOCAL ${FASTA_RAW_LOCAL}"; fi;
	@if test -s ${FASTA_MSK_LOCAL}; then echo "	${FASTA_MSK_LOCAL}"; gzip -f ${FASTA_MSK_LOCAL}; else echo "	skipping FASTA_MSK_LOCAL ${FASTA_MSK_LOCAL}"; fi;
	@if test -s ${FASTA_PEP_LOCAL}; then echo "	${FASTA_PEP_LOCAL}"; gzip -f ${FASTA_PEP_LOCAL}; else echo "	skipping FASTA_PEP_LOCAL ${FASTA_PEP_LOCAL}"; fi;

################################################################
## Initialize the fasta indexes for bedtools getfasta.
##
## THIS TARGET IS NOT REQUIRED ANYMORE: replaced by samtools faidx in the script parse-gtf.
## It can be used as a test
RSAT_GTF=${RSAT}/public_html/data/genomes/${SPECIES_RSAT_ID}/genome/${SPECIES_RSAT_ID}.gtf
START_CODONS=${RSAT}/public_html/data/genomes/${SPECIES_RSAT_ID}/genome/${SPECIES_RSAT_ID}_start_codons
init_getfasta:
	@echo
	@echo "Initializing genomic fasta index	${SPECIES_RSAT_ID}"
	@grep start_codon ${RSAT_GTF} > ${START_CODONS}.gtf || true
	@echo "	Start codon coordinates	${START_CODONS}.gtf"
	@retrieve-seq-bed -i ${START_CODONS}.gtf -o ${START_CODONS}.fasta -org ${SPECIES_RSAT_ID}
	@echo "	Start codon sequences	${START_CODONS}.fasta"
	@oligo-analysis -v 1 -i ${START_CODONS}.fasta -l 3 -1str -return occ,freq -o ${START_CODONS}_3nt_freq.tab
	@echo "	Start codon frequencies	${START_CODONS}_3nt_freq.tab"
	@retrieve-seq-bed -i ${START_CODONS}.gtf -o ${START_CODONS}-rm.fasta -org ${SPECIES_RSAT_ID} -rm
	@echo "	Repeat-masked start codon sequences	${START_CODONS}-rm.fasta"
	@oligo-analysis -v 1 -i ${START_CODONS}-rm.fasta -l 3 -1str -return occ,freq -o ${START_CODONS}-rm_3nt_freq.tab
	@echo "	Repeat-masked start codon frequencies	${START_CODONS}-rm_3nt_freq.tab"

################################################################
## Install some pet genomes

## Arabidopsis thaliana (Plant)
install_thaliana:
	${MAKE} GROUP=Plants SPECIES=arabidopsis_thaliana organisms ${DOWNLOAD_TASKS} ${INSTALL_TASKS}

## Saccharomyces cerevisiae (Fungus)
install_yeast:
	${MAKE} GROUP=Fungi SPECIES=saccharomyces_cerevisiae COLLECTION= organisms ${DOWNLOAD_TASKS} ${INSTALL_TASKS}


## Mus musculus (Metazoa)
## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
## NOT YET WORKING, because of inconsistencies between ensembl and ensemblgenomes FTP servers:
## - ensembl FTP site has not the species table which we use to get assembly name
## - ensembl FTP site does not contain the single file per genome for DNA sequences. We should write a specific target to concatenate all the chromosome files.
## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
install_mouse:
	${MAKE} GROUP=Metazoa SPECIES=mus_musculus SERVER_URL=ftp://ftp.ensembl.org/pub \
		RELEASE=${ENSEMBL_RELEASE} ${DOWNLOAD_TASKS} ${INSTALL_TASKS}

## Drosophila melanogaster (Metazoa)
install_droso:
	${MAKE} GROUP=Metazoa SPECIES=drosophila_melanogaster organisms ${DOWNLOAD_TASKS} ${INSTALL_TASKS}

## Note: for bacteria we need to define a collection

## Escherichia coli (Bacteria)
install_ecoli:
	${MAKE} GROUP=Bacteria SPECIES=escherichia_coli_str_k_12_substr_mg1655 \
		COLLECTION=bacteria_0_collection organisms ${DOWNLOAD_TASKS} ${INSTALL_TASKS}


## Pseudomonas aeruginosa (Bacteria)
#install_pao1:
#	${MAKE} GROUP=Bacteria SPECIES=pseudomonas_aeruginosa_pao1_ve13 \
#		COLLECTION=bacteria_44_collection ${INSTALL_TASKS}

install_bsub:
	${MAKE} GROUP=Bacteria SPECIES=bacillus_subtilis_subsp_subtilis_str_168 \
		COLLECTION=bacteria_0_collection organisms ${DOWNLOAD_TASKS} ${INSTALL_TASKS}


##################################################################
## Parse Compara.homologies 
#CMP_GZ=$(shell ls -1 ${ORGANISM_DIR}/Compara.homologies*.gz)
CMP_GZ=${ORGANISM_DIR}/Compara.${ENSEMBL_RELEASE}.protein_default.homologies.tsv.gz
BDB_FILE=${ORGANISM_DIR}/compara.bdb
BDB_LOG=${ORGANISM_DIR}/compara.log
parse_compara:
	@echo
	@echo "Parsing Compara file ${CMP_GZ}"
	@echo
	@parse-compara -i ${CMP_GZ} -list ${ORGANISM_TABLE} -release ${RELEASE} \
		-o ${BDB_FILE} -log ${BDB_LOG} -v ${V}

##################################################################
## Parse Compara.homologies and match genomes names to support-organisms 
CMP_GZ=${ORGANISM_DIR}/Compara.${ENSEMBL_RELEASE}.protein_default.homologies.tsv.gz
BDB_FILE=${ORGANISM_DIR}/compara.bdb
BDB_LOG=${ORGANISM_DIR}/compara.log
parse_compara_match:
	@echo
	@echo "Parsing Compara file ${CMP_GZ}"
	@echo
	@parse-compara -i ${CMP_GZ} -list ${ORGANISM_TABLE} -match_genomes \
        -o ${BDB_FILE} -log ${BDB_LOG} -v ${V}

#################################################################
## Install Compara db
COMP_INSTALL_DIR=${RSAT}/public_html/data/genomes/
install_compara:
	@echo
	@echo "Installing Compara db ${BDB_FILE}"
	@echo
	@mv ${BDB_FILE} ${COMP_INSTALL_DIR}
	@echo
	@ls -1 ${COMP_INSTALL_DIR}/compara.bdb

##################################################################

# not used like that, requires previous download_go
#all: organisms \
	${DOWNLOAD_TASKS} \
	${INSTALL_TASKS} \
	${COMPARA_TASKS}

clean_compara:
	@echo
	@echo "Deleting ensemblgenomes Compara release ${RELEASE}"
	@[[ -e ${CMP_GZ} ]] && rm -f ${CMP_GZ}
	@echo

clean_all:
	@echo
	@echo "Deleting ensemblgenomes release ${RELEASE}" 
	@[[ -d ${ORGANISM_DIR} ]] && rm -rf ${ORGANISM_DIR}
	@echo	

clean:
	@echo
	@echo "Deleting ensemblgenomes species ${SPECIES} (release ${RELEASE})"
	@[[ -d ${SPECIES_DIR} ]] && rm -rf ${SPECIESS_DIR}
	@echo	

