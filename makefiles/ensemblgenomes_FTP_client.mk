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
	@wget -cnv ${SERVERLIST} -O ${ORGANISMS_LIST} 
	@echo
	@echo "	${ORGANISMS_LIST}"

list_param:
	@echo
	@echo "Parameters"
	@echo "	GROUP   ${GROUP} (${GROUP_LC})"
	@echo "	SPECIES	${SPECIES}"
	@echo "	RELEASE ${RELEASE}"
	@echo "Files to download"
	@echo "	GTF_FTP_URL		${GTF_FTP_URL}"
	@echo "	FASTA_RAW_FTP_URL	${FASTA_RAW_FTP_URL}"
	@echo "	FASTA_MSK_FTP_URL	${FASTA_MSK_FTP_URL}"
	@echo "	FASTA_PEP_FTP_URL	${FASTA_PEP_FTP_URL}"
	@echo "LOCAL_FILES"
	@echo "	GTF_LOCAL		${GTF_LOCAL}"
	@echo "	FASTA_RAW_LOCAL		${FASTA_RAW_LOCAL}"
	@echo "	FASTA_MSK_LOCAL		${FASTA_MSK_LOCAL}"
	@echo "	FASTA_PEP_LOCAL		${FASTA_PEP_LOCAL}"

################################################################
## Download all required files
download_all: download_gtf download_fasta

################################################################
## Download GTF files from ensemblgenomes
GTF_FTP_URL=${DATABASE}/gtf/${SPECIES}/*${RELEASE}.gtf.gz
download_gtf:
	@echo
	@mkdir -p ${SPECIES_DIR}	
	@echo "Downloading GTF file of ${SPECIES}"	
	@wget -Ncnv ${GTF_FTP_URL} -P ${SPECIES_DIR}
	@echo
	@ls -1 ${SPECIES_DIR}/*.gtf.gz

################################################################
## Download FASTA files with genomic sequences (raw and masked)
## and peptidic sequences
FASTA_RAW_SUFFIX=*${RELEASE}.dna.genome.fa.gz
FASTA_RAW_FTP_URL=${DATABASE}/fasta/${SPECIES}/dna/${FASTA_RAW_SUFFIX}
FASTA_MSK_SUFFIX=*${RELEASE}.dna_rm.genome.fa.gz
FASTA_MSK_FTP_URL=${DATABASE}/fasta/${SPECIES}/dna/${FASTA_MSK_SUFFIX}
FASTA_PEP_SUFFIX=*${RELEASE}.pep.all.fa.gz
FASTA_PEP_FTP_URL=${DATABASE}/fasta/${SPECIES}/pep/${FASTA_PEP_SUFFIX}
download_fasta:
	@echo
	@mkdir -p ${SPECIES_DIR}
	@echo "Downloading raw FASTA genome for species ${SPECIES}"
	@wget -Ncnv ${FASTA_RAW_FTP_URL} -P ${SPECIES_DIR}
	@echo
	@echo "Downloading repeat-masked FASTA genome for species ${SPECIES}"
	@wget -Ncnv ${FASTA_MSK_FTP_URL} -P ${SPECIES_DIR}
	@echo
	@echo
	@echo "Downloading FASTA peptidic sequences for species ${SPECIES}"
	@wget -Ncnv ${FASTA_PEP_FTP_URL} -P ${SPECIES_DIR}
	@echo
	@ls -1 ${SPECIES_DIR}/*.fa.gz

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
FASTA_RAW_LOCAL=`ls -1 ${SPECIES_DIR}/${FASTA_RAW_SUFFIX} | head -1`
FASTA_MSK_LOCAL=`ls -1 ${SPECIES_DIR}/${FASTA_MSK_SUFFIX} | head -1`
FASTA_PEP_LOCAL=`ls -1 ${SPECIES_DIR}/${FASTA_PEP_SUFFIX} | head -1`
GTF_LOCAL=$(shell ls -1 ${SPECIES_DIR}/*.gtf.gz)
PARSE_DIR=${SPECIES_DIR}
PARSE_TASK=all
# Note that only the first file is considered
parse_gtf:
	@echo
	@echo "Parsing GTF file	${GTF_LOCAL}"
	parse-gtf -v ${V} -i ${GTF_LOCAL} \
		-fasta ${FASTA_RAW_LOCAL} \
		-fasta_rm ${FASTA_MSK_LOCAL} \
		-fasta_pep ${FASTA_PEP_LOCAL} \
		-org_name ${SPECIES} \
		-task ${PARSE_TASK} ${OPT} \
		-o ${PARSE_DIR} 
	@echo "	${PARSE_DIR}"
#	@ls -1 ${PARSE_DIR}/*.tab


install_from_gtf:
	@echo
	@echo "Parsing and installing in RSAT	${SPECIES}"
	@${MAKE} parse_gtf PARSE_DIR=${RSAT}/public_html/data/genomes/${SPECIES}/genome

## Run some test for the GTF parsing result
parse_gtf_test:
	retrieve-seq -org ${SPECIES} -from 0 -to 3 -feattype gene | oligo-analysis -v 1 -l 3 -return occ,freq -sort 

install_yeast:
	${MAKE} GROUP=Fungi SPECIES=saccharomyces_cerevisiae download_gtf download_fasta install_from_gtf

##################################################################
## Parse Compara.homologies 
CMP_GZ=$(shell ls -1 ${ORGANISMS_DIR}/Compara.homologies*.gz)
parse_compara:
	@echo
	@echo "Parsing Compara file ${CMP_GZ}"
	@echo
#@ls -1 ${PARSE_DIR}/*.tab


##################################################################

all: organisms download_fasta download_gtf download_compara \
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



