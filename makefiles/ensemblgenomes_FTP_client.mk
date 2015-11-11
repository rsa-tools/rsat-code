################################################################
## This makefile contains some targets to download genome seqs and
## annotations from ensemblgenome FTP site, parseand install them
## Jacques Van Helden, Bruno Contreras Moreira

include ${RSAT}/makefiles/util.mk
MAKEFILE=${RSAT}/makefiles/ensemblgenomes_FTP_client.mk

## Define parameters
V=2


#should be set in env var {ORG_GROUP} ?
#currently does not work for Bacteria, as these are further grouped in bacteria_NN_collection subfolders
## GROUP=Plants I moved this variable to RSAT_config.mk, since it depends on the server;
GROUP_LC=$(shell echo $(GROUP) | tr A-Z a-z)
RELEASE=${ENSEMBLGENOMES_BRANCH}
# should be set in RSAT_config.props
SERVER_URL=ftp://ftp.ensemblgenomes.org/pub/${GROUP_LC}
DATABASE=${SERVER_URL}/release-${RELEASE}

#name preffix hard-coded, might change in future
SERVERLIST=${DATABASE}/species_Ensembl${GROUP}.txt

ORGANISMS_DIR=${RSAT}/data/ensemblgenomes/${GROUP_LC}/release-${RELEASE}
ORGANISMS_LIST=${ORGANISMS_DIR}/species_Ensembl${GROUP}.txt
## SPECIES=arabidopsis_thaliana ## The default species is coupled to the group specificity of the server -> I (JvH) move it to RSAT_config.mk

## Note (JvH 2015-11-06) I change SPECIES DIR to directly download
## fasta and gtf in the genome dir, since we will use it for vairous
## purposes.
SPECIES_UCFIRST=$(shell perl -e 'print ucfirst ${SPECIES}')
SPECIES_RSAT_ID=${SPECIES_UCFIRST}.${ASSEMBLY_ID}.${RELEASE}
# SPECIES_DIR=${ORGANISMS_DIR}/${SPECIES}
SPECIES_DIR=${RSAT}/data/genomes/${SPECIES_RSAT_ID}
GENOME_DIR=${SPECIES_DIR}/genome

###############################################################
## Get all supported organisms in an ensemblgenome release and store them in a file
organisms:
	@echo
	@mkdir -p ${ORGANISMS_DIR}
	@echo "Getting list of organisms from ${DATABASE}"
	@echo "	${ORGANISMS_DIR}"
	@wget -Ncnv ${SERVERLIST} -P ${ORGANISMS_DIR}
	@echo
	@echo "	${ORGANISMS_LIST}"

list_param:
	@echo
	@echo "Parameters"
	@echo "	GROUP   		${GROUP} (${GROUP_LC})"
	@echo "	SPECIES			${SPECIES} (${SPECIES_UCFIRST})"
	@echo "	TAXON_ID 		${TAXON_ID}"
	@echo "	ASSEMBLY_ID 		${ASSEMBLY_ID}"
	@echo "	RELEASE 		${RELEASE}"
	@echo "	SPECIES_RSAT_ID		${SPECIES_RSAT_ID}"
	@echo "Files to download"
	@echo "	GTF_FTP_URL		${GTF_FTP_URL}"
	@echo "	FASTA_RAW_FTP_URL	${FASTA_RAW_FTP_URL}"
	@echo "	FASTA_MSK_FTP_URL	${FASTA_MSK_FTP_URL}"
	@echo "	FASTA_PEP_FTP_URL	${FASTA_PEP_FTP_URL}"
	@echo "	SERVER_COMPARA_FILE	${SERVER_COMPARA_FILE}"
	@echo "LOCAL_FILES"
	@echo "	ORGANISMS_DIR		${ORGANISMS_DIR}"
	@echo "	SPECIES_DIR		${SPECIES_DIR}"
	@echo "	GENOME_DIR		${GENOME_DIR}"
	@echo "	PARSE_DIR		${PARSE_DIR}"
	@echo "	GO_DIR			${GO_DIR}"
	@echo "	GTF_LOCAL		${GTF_LOCAL}"
	@echo "	FASTA_RAW_LOCAL		${FASTA_RAW_LOCAL}"
	@echo "	FASTA_MSK_LOCAL		${FASTA_MSK_LOCAL}"
	@echo "	FASTA_PEP_LOCAL		${FASTA_PEP_LOCAL}"
	@echo "	CMP_GZ			${CMP_GZ}"

################################################################
## Download required files for all organisms
ALL_SPECIES=$(shell cut -f 2 ${ORGANISMS_LIST} | grep -v '^species')
list_all_species:
	@echo 
	@echo "All species from ${ORGANISMS_LIST}"
	@echo ${ALL_SPECIES} | perl -pe 's/\s+/\n/g' |add-linenb -before


COLLECTION=
ORG_TASKS=organisms
DOWNLOAD_TASKS=download_gtf download_fasta  gunzip_fasta 
INSTALL_TASKS=install_from_gtf init_getfasta install_go_annotations
COMPARA_TASKS=parse_compara install_compara
ALL_TASKS=${ORG_TASKS} ${DOWNLOAD_TASKS} ${INSTALL_TASKS} ${COMPARA_TASKS}

download_one_species: ${DOWNLOAD_TASKS}

install_one_species: ${INSTALL_TASKS}

download_all_species: organisms
	@echo WARNING: Make sure you run organisms before download_all_species
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
	@${MAKE} parse_compara
	@${MAKE} install_compara

################################################################
## Check upstrean sequences of all installed species
check_all_species:
	@echo WARNING: Make sure you install_all_species before check_all_species
	@echo
	@echo Checking upstraem sequences of all species in GROUP=${GROUP} RELEASE=${RELEASE}
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
	@wget -Ncnv ${GTF_FTP_URL} -P ${GENOME_DIR}
	@echo
	@ls -1 ${GENOME_DIR}/*.gtf.gz

################################################################
## Download FASTA files with genomic sequences (raw and masked)
## and peptidic sequences
FASTA_RAW_SUFFIX=*${RELEASE}.dna.genome.fa*
FASTA_RAW_FTP_URL=${DATABASE}/fasta/${COLLECTION}/${SPECIES}/dna/${FASTA_RAW_SUFFIX}
FASTA_MSK_SUFFIX=*${RELEASE}.dna_rm.genome.fa*
FASTA_MSK_FTP_URL=${DATABASE}/fasta/${COLLECTION}/${SPECIES}/dna/${FASTA_MSK_SUFFIX}
FASTA_PEP_SUFFIX=*${RELEASE}.pep.all.fa*
FASTA_PEP_FTP_URL=${DATABASE}/fasta/${COLLECTION}/${SPECIES}/pep/${FASTA_PEP_SUFFIX}
download_fasta:
	@echo
	@mkdir -p ${GENOME_DIR}
	@echo "Downloading raw FASTA genome for species ${SPECIES}"
	@echo "	GENOME_DIR	${GENOME_DIR}"
	@wget -Ncnv ${FASTA_RAW_FTP_URL} -P ${GENOME_DIR}
	@echo
	@echo "Downloading repeat-masked FASTA genome for species ${SPECIES}"
	@wget -Ncnv ${FASTA_MSK_FTP_URL} -P ${GENOME_DIR}
	@echo
	@echo
	@echo "Downloading FASTA peptidic sequences for species ${SPECIES}"
	@wget -Ncnv ${FASTA_PEP_FTP_URL} -P ${GENOME_DIR}
	@echo
	@ls -1 ${GENOME_DIR}/*.fa.gz


################################################################
## Download a sample of eg upstream sequences to check installed sequences
check_sequences:
	@echo
	@check-retrieve-seq-rest -v ${V} \
        -org ${SPECIES} \

#################################################################
## Download group COMPARA files from eg
SERVER_COMPARA_FILE=${DATABASE}/tsv/ensembl-compara/Compara.homologies.${RELEASE}.tsv.gz
download_compara:
	@echo
	@mkdir -p ${ORGANISMS_DIR}
	@echo "Downloading COMPARA file of ${GROUP}"
	@echo "	${SERVER_COMPARA_FILE}"
	@wget -Ncnv ${SERVER_COMPARA_FILE} -P ${ORGANISMS_DIR}
	@echo
	@ls -1 ${ORGANISMS_DIR}/Compara.homologies*gz

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
	@download-ensembl-go-annotations-biomart -o ${GO_ANNOT_FILE} -org ${SPECIES} \
		-release ${RELEASE}	-list ${ORGANISMS_LIST-}
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
## cluster.

FASTA_RAW_LOCAL=`ls -1 ${GENOME_DIR}/${FASTA_RAW_SUFFIX} | head -1`
FASTA_MSK_LOCAL=`ls -1 ${GENOME_DIR}/${FASTA_MSK_SUFFIX} | head -1`
FASTA_PEP_LOCAL=`ls -1 ${GENOME_DIR}/${FASTA_PEP_SUFFIX} | head -1`
# Note that only the first gtf file is considered
GTF_LOCAL=$(shell ls -1 ${GENOME_DIR}/*.gtf.gz ${GENOME_DIR}/*.gtf | head -1)
TAXON_ID=$(shell grep -w ${SPECIES} ${ORGANISMS_LIST} | cut -f 4)
ASSEMBLY_ID=$(shell grep -w ${SPECIES} ${ORGANISMS_LIST} | cut -f 5)
PARSE_DIR=${GENOME_DIR}
PARSE_TASK="parse_gtf,parse_fasta"
PARSE_GTF_CMD=parse-gtf -v ${V} -i ${GTF_LOCAL} \
		-fasta ${FASTA_RAW_LOCAL} \
		-fasta_rm ${FASTA_MSK_LOCAL} \
		-fasta_pep ${FASTA_PEP_LOCAL} \
		-org_name ${SPECIES_RSAT_ID} \
		-task ${PARSE_TASK} ${OPT} \
		-taxid ${TAXON_ID} \
		-gtf_source ensemblgenomes \
		-o ${PARSE_DIR} 
parse_gtf:
	@echo
	@echo "Parsing GTF file	${GTF_LOCAL}"
	@echo "TaxonID = ${TAXON_ID}"
	@${MAKE} my_command MY_COMMAND="${PARSE_GTF_CMD}"
	@echo "	${PARSE_DIR}"
#	@ls -1 ${PARSE_DIR}/*.tab

###############################################################
## parse gtf and then install organism
install_from_gtf:
	@echo
	@echo "Parsing and installing in RSAT	${SPECIES}"
	@${MAKE} parse_gtf PARSE_DIR=${RSAT}/public_html/data/genomes/${SPECIES_RSAT_ID}/genome PARSE_TASK="all"

## Run some test for the GTF parsing result
parse_gtf_test:
	retrieve-seq -org ${SPECIES} -from 0 -to 3 -feattype gene | oligo-analysis -v 1 -l 3 -return occ,freq -sort 
###############################################################
## Uncompress genomic fasta files for beedtools
gunzip_fasta:
	@echo
	@echo "Uncompressing the downloaded fasta files."
	@echo ${FASTA_RAW_LOCAL}
	@gunzip	${FASTA_RAW_LOCAL}
	@echo ${FASTA_MSK_LOCAL}
	@gunzip	${FASTA_MSK_LOCAL}
	@echo ${FASTA_PEP_LOCAL}
	@gunzip	${FASTA_PEP_LOCAL}

################################################################
## Initialize the fasta indexes for bedtools getfasta.
RSAT_GTF=${RSAT}/public_html/data/genomes/${SPECIES_RSAT_ID}/genome/${SPECIES_RSAT_ID}.gtf.gz
START_CODONS=${RSAT}/public_html/data/genomes/${SPECIES_RSAT_ID}/genome/${SPECIES_RSAT_ID}_start_codons
init_getfasta:
	@echo
	@echo "Initializing genomic fasta index	${SPECIES_RSAT_ID}"
	@zgrep start_codon ${RSAT_GTF} > ${START_CODONS}.gtf
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
	${MAKE} GROUP=Plants SPECIES=arabidopsis_thaliana ${INSTALL_TASKS}

## Zea mais (Plant)
install_zea:
	${MAKE} GROUP=Plants SPECIES=zea_mays ${INSTALL_TASKS}

## Saccharomyces cerevisiae (Fungus)
install_yeast:
	${MAKE} GROUP=Fungi SPECIES=saccharomyces_cerevisiae ${INSTALL_TASKS}

## Mus musculus (Metazoa)
## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
## NOT YET WORKING, because of inconsistencies between ensembl and ensemblgenomes FTP servers:
## - ensembl FTP site has not the species table which we use to get assembly name
## - ensembl FTP site does not contain the single file per genome for DNA sequences. We should write a specific target to concatenate all the chromosome files.
## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
install_mouse:
	${MAKE} GROUP=Metazoa SPECIES=mus_musculus SERVER_URL=ftp://ftp.ensembl.org/pub RELEASE=82 ${INSTALL_TASKS}

## Drosophila melanogaster (Metazoa)
install_droso:
	${MAKE} GROUP=Metazoa SPECIES=drosophila_melanogaster ${INSTALL_TASKS}

## Note: for bacteria we need to define a collection

## Escherichia coli (Bacteria)
install_ecoli:
	${MAKE} GROUP=Bacteria SPECIES=escherichia_coli_str_k_12_substr_mg1655 \
		COLLECTION=bacteria_0_collection ${INSTALL_TASKS}

## Pseudomonas aeruginosa (Bacteria)
#install_pao1:
#	${MAKE} GROUP=Bacteria SPECIES=pseudomonas_aeruginosa_pao1_ve13 \
#		COLLECTION=bacteria_44_collection ${INSTALL_TASKS}

install_bsub:
	${MAKE} GROUP=Bacteria SPECIES=bacillus_subtilis_subsp_subtilis_str_168 COLLECTION=bacteria_0_collection ${INSTALL_TASKS}


##################################################################
## Parse Compara.homologies 
CMP_GZ=$(shell ls -1 ${ORGANISMS_DIR}/Compara.homologies*.gz)
BDB_FILE=${ORGANISMS_DIR}/compara.bdb
BDB_LOG=${ORGANISMS_DIR}/compara.log
parse_compara:
	@echo
	@echo "Parsing Compara file ${CMP_GZ}"
	@echo
	@parse-compara -i ${CMP_GZ} -list ${ORGANISMS_LIST} -release ${RELEASE} \
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

all: organisms \
	${DOWNLOAD_TASKS} \
	${INSTALL_TASKS} \
	parse_compara install_compara

clean_compara:
	@echo
	@echo "Deleting ensemblgenomes Compara release ${RELEASE}"
	@[[ -e ${CMP_GZ} ]] && rm -f ${CMP_GZ}
	@echo

clean_all:
	@echo
	@echo "Deleting ensemblgenomes release ${RELEASE}" 
	@[[ -d ${ORGANISMS_DIR} ]] && rm -rf ${ORGANISMS_DIR}
	@echo	

clean:
	@echo
	@echo "Deleting ensemblgenomes species ${SPECIES} (release ${RELEASE})"
	@[[ -d ${SPECIES_DIR} ]] && rm -rf ${SPECIESS_DIR}
	@echo	

