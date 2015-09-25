################################################################
## This makefile contains some targets to download genome seqs and
## annotations from ensemblgenome FTP site
## Jacques Van Helden, Bruno Contreras Moreira

include ${RSAT}/makefiles/util.mk
MAKEFILE=${RSAT}/makefiles/ensemblgenomes_FTP_client.mk

## Define parameters
V=2

# should be set in RSAT_config.props
SERVERURL=ftp://ftp.ensemblgenomes.org

#should be set in env var {ORG_GROUP} ?
#currently does not work for Bacteria, as these are further grouped in bacteria_NN_collection subfolders
GROUP=Plants
GROUP_LC=$(shell echo $(GROUP) | tr A-Z a-z)
RELEASE=${ENSEMBLGENOMES_BRANCH}
DATABASE=${SERVERURL}/pub/${GROUP_LC}/release-${RELEASE}

#name preffix hard-coded, might change in future
SERVERLIST=${DATABASE}/species_Ensembl${GROUP}.txt

ORGANISMS_DIR=${RSAT}/data/ensemblgenomes/${GROUP_LC}/release-${RELEASE}
ORGANISMS_LIST=${ORGANISMS_DIR}/organisms.tab
SPECIES=chlamydomonas_reinhardtii
SPECIES_DIR=${ORGANISMS_DIR}/${SPECIES}

###############################################################
## Get all supported organisms in an eg release and store them in a file
organisms:
	@echo
	@mkdir -p ${ORGANISMS_DIR}
	@echo "Getting list if organisms from ${DATABASE}"
	@wget -Ncnv ${SERVERLIST} -O ${ORGANISMS_LIST} 
	@echo
	@echo "	${ORGANISMS_LIST}"

list_param:
	@echo
	@echo "Parameters"
	@echo "	GROUP   ${GROUP} (${GROUP_LC})"
	@echo "	SPECIES	${SPECIES}"
	@echo "	RELEASE ${RELEASE}"

################################################################
## Download GTF files from ensemblgenomes
SERVER_GTF_FILE=${DATABASE}/gtf/${SPECIES}/*${RELEASE}.gtf.gz

download_gtf:
	@echo
	@mkdir -p ${SPECIES_DIR}	
	@echo "Downloading GTF file of ${SPECIES}"	
	@wget -Ncnv ${SERVER_GTF_FILE} -P ${SPECIES_DIR}
	@echo
	@ls -1 ${SPECIES_DIR}/*.gtf.gz

################################################################
## Download genome FASTA files (raw and masked) from eg 
FASTA_RAW_SUFFIX=*${RELEASE}.dna.genome.fa.gz
FASTA_RM_SUFFIX=*${RELEASE}.dna_rm.genome.fa.gz
FASTA_RAW_FTP_URL=${DATABASE}/fasta/${SPECIES}/dna/${FASTA_RAW_SUFFIX}
FASTA_MSK_FTP_URL=${DATABASE}/fasta/${SPECIES}/dna/${FASTA_RM_SUFFIX}
download_fasta:
	@echo
	@mkdir -p ${SPECIES_DIR}
	@echo "Downloading FASTA genome files of ${SPECIES}"
	@wget -Ncnv ${FASTA_RAW_FTP_URL} -P ${SPECIES_DIR}
	@echo
	@wget -Ncnv ${FASTA_MSK_FTP_URL} -P ${SPECIES_DIR}
	@echo
	@ls -1 ${SPECIES_DIR}/*.genome.fa.gz

################################################################
## Download sequences of some eg genomic features to be used as control
## of RSAT scripts that slice sequences based on coordinates
SERVER_CDS_FILE=${DATABASE}/fasta/${SPECIES}/cds/*${RELEASE}.cds.all.fa.gz
download_feature_sequences:
	@echo
	@mkdir -p ${SPECIES_DIR}
	@echo "Downloading FASTA feature files of ${SPECIES}"
	@wget -Ncnv ${SERVER_CDS_FILE} -P ${SPECIES_DIR}
	@echo
	@ls -1 ${SPECIES_DIR}/*.cds.all.fa.gz

#################################################################
## Download group COMPARA files from eg
SERVER_COMPARA_FILE=${DATABASE}/tsv/ensembl-compara/Compara.homologies.${RELEASE}.tsv.gz
download_compara:
	@echo
	@mkdir -p ${ORGANISMS_DIR}
	@echo "Downloading COMPARA file of ${GROUP}"
	@wget -Ncnv ${SERVER_COMPARA_FILE} -P ${ORGANISMS_DIR}
	@echo
	@ls -1 ${ORGANISMS_DIR}/Compara.homologies*gz

##################################################################
## Parse GTF file to extract gene, transcripts and cds coords
FASTA_RAW=`ls -1 ${SPECIES_DIR}/${FASTA_RAW_SUFFIX} | head -1`
GTF_GZ=$(shell ls -1 ${SPECIES_DIR}/*.gtf.gz)
# Note that only the first file is considered
parse_gtf:
	@echo
	@echo "Parsing GTF file	${GTF_GZ}"
	parse-gtf -v ${V} -i ${GTF_GZ} -fasta ${FASTA_RAW} -o ${SPECIES_DIR} 
	@echo
	@ls -1 ${SPECIES_DIR}/*.tab

##################################################################
## Parse Compara.homologies 
CMP_GZ=$(shell ls -1 ${ORGANISMS_DIR}/Compara.homologies*.gz)
parse_compara:
	@echo
	@echo "Parsing Compara file ${CMP_GZ}"
	@echo
#@ls -1 ${SPECIES_DIR}/*.tab


##################################################################

all: download_fasta download_feature_sequences download_gtf download_compara \
	parse_gtf parse_compara

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



