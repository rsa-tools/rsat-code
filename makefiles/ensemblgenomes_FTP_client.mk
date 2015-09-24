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
GROUP=plants

#GROUP_LC=`echo $(GROUP) | tr A-Z a-z`
RELEASE=${ENSEMBLGENOMES_BRANCH}
DATABASE=${SERVERURL}/pub/${GROUP}/release-${RELEASE}

#name preffix hard-coded, might change in future
SERVERLIST=${DATABASE}/species_Ensembl${GROUP}.txt

ORGANISMS_DIR=${RSAT}/data/ensemblgenomes/${GROUP}/release-${RELEASE}
ORGANISMS_LIST=${ORGANISMS_DIR}/organisms.tab
SPECIES=chlamydomonas_reinhardtii

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
	@echo "	GROUP		${GROUP}"
	@echo "	SPECIES		${SPECIES}"
	@echo "	RELEASE		${RELEASE}"

################################################################
## Download GTF files from ensemblgenomes
SERVER_GTF_FILE=${DATABASE}/gtf/${SPECIES}/*${RELEASE}.gtf.gz

download_gtf:
	@echo
	@mkdir -p ${ORGANISMS_DIR}/${SPECIES}	
	@echo "Downloading GTF file of ${SPECIES}"	
	@wget -Ncnv ${SERVER_GTF_FILE} -P ${ORGANISMS_DIR}/${SPECIES}
	@echo
	@ls -1 ${ORGANISMS_DIR}/${SPECIES}/*.gtf.gz

################################################################
## Download genome FASTA files (raw and masked) from eg 
SERVER_RAW_FILE=${DATABASE}/fasta/${SPECIES}/dna/*${RELEASE}.dna.genome.fa.gz
SERVER_MSK_FILE=${DATABASE}/fasta/${SPECIES}/dna/*${RELEASE}.dna_rm.genome.fa.gz
download_fasta:
	@echo
	@mkdir -p ${ORGANISMS_DIR}/${SPECIES}
	@echo "Downloading FASTA genome files of ${SPECIES}"
	@wget -Ncnv ${SERVER_RAW_FILE} -P ${ORGANISMS_DIR}/${SPECIES}
	@echo
	@wget -Ncnv ${SERVER_MSK_FILE} -P ${ORGANISMS_DIR}/${SPECIES}
	@echo
	@ls -1 ${ORGANISMS_DIR}/${SPECIES}/*.fa.gz

################################################################
## Download sequences of some eg genomic features to be used as control
## of RSAT scripts that slice sequences based on coordinates
SERVER_CDS_FILE=${DATABASE}/fasta/${SPECIES}/cds/*${RELEASE}.cds.all.fa.gz
download_feature_sequences:
	@echo
	@mkdir -p ${ORGANISMS_DIR}/${SPECIES}
	@echo "Downloading FASTA feature files of ${SPECIES}"
	@wget -Ncnv ${SERVER_CDS_FILE} -P ${ORGANISMS_DIR}/${SPECIES}
	@echo
	@ls -1 ${ORGANISMS_DIR}/${SPECIES}/*.cds.all.fa.gz

#################################################################
## Download group COMPARA files from eg
SERVER_COMPARA_FILE=${DATABASE}/tsv/ensembl-compara/Compara.homologies.${RELEASE}.tsv.gz
download_compara:
	@echo
	@mkdir -p ${ORGANISMS_DIR}
	@echo "Downloading COMPARA file of ${GROUP}"
	@wget -Ncnv ${SERVER_COMPARA_FILE} -P ${ORGANISMS_DIR}
	@echo

##################################################################
## Parse GTF file to extract gene, transcripts a
GTF_GZ=`ls -1 ${ORGANISMS_DIR}/${SPECIES}/*.gtf.gz`
PARSE_DIR=${RSAT}/data/genomes/${SPECIES}/genome/
parse_gtf:
	@echo
	@echo "Parsing GTF file	${GTF_GZ}"
	@echo "SPECIES	${SPECIES}"
	parse-gtf -v ${V} -i ${GTF_GZ} -o ${PARSE_DIR}
	@echo "	${PARSE_DIR}"


##################################################################

all: download_fasta download_feature_sequences download_gtf download_compara

clean_all:
	@echo
	@echo "Deleting ensemblgenomes release ${RELEASE}" 
	@[[ -d ${ORGANISMS_DIR} ]] && rm -rf ${ORGANISMS_DIR}
	@echo	

clean:
	@echo
	@echo "Deleting ensemblgenomes species ${SPECIES} (release ${RELEASE})"
	@[[ -d ${ORGANISMS_DIR}/${SPECIES} ]] && rm -rf ${ORGANISMS_DIR}/${SPECIES}
	@echo	

