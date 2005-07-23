############################################################
#
# $Id: install_genomes.mk,v 1.20 2005/07/23 05:57:17 rsat Exp $
#
# Time-stamp: <2003-10-10 22:49:55 jvanheld>
#
############################################################

include ${RSAT}/makefiles/util.mk

DATE = `date +%Y%m%d_%H%M%S`

################################################################
#### Directories
GENBANK_DIR=${RSAT}/downloads/ftp.ncbi.nih.gov/genbank/genomes
NCBI_DIR=${RSAT}/downloads/ftp.ncbi.nih.gov/genomes

#################################################################
#### Programs

WGET = wget -np -rNL 
MAKEFILE=${RSAT}/makefiles/install_genomes.mk
#MAKE=nice -n 19 make -f ${MAKEFILE}
RSYNC_OPT = -ruptvl ${OPT}
SSH=-e 'ssh -x'
RSYNC = rsync ${RSYNC_OPT} ${SSH}

V=1 

################################################################
### Targets

### Install one organism
ORG=Arabidopsis_thaliana
ORG_DIR=${NCBI_DIR}/${ORG}
INSTALL_TASK=allup,config,dyads,oligos,parse,start_stop,upstream_freq,phylogeny
INSTALL_CMD=install-organism -v ${V}		\
		-genbank ${NCBI_DIR}		\
		-org ${ORG}			\
		-task ${INSTALL_TASK}		\
		${OPT}

install_one_organism:
	@echo "install log	${INSTALL_LOG}"
	@echo "Installing organism ${ORG}" 
	${MAKE} my_command MY_COMMAND="${INSTALL_CMD}" JOB_PREFIX=install_${ORG}

### Bacteria with a small genome for quick testing
BACT=Mycoplasma_genitalium

### All the bacteria in NCBI genome directory
BACTERIA = `ls -1 ${NCBI_DIR}/Bacteria | grep _ | sort -u | xargs `
list_bacteria:
	@echo "Bacteria to install	${BACTERIA}"


### Install all bacterial genomes on RSAT
install_all_bacteria:
	@for bact in ${BACTERIA}; do				\
		${MAKE} install_one_bacteria BACT=$${bact} ;	\
	done

### Install a single bacteria genome
install_one_bacteria:
	@echo
	@echo "${DATE}	Installing bacteria ${BACT}"
	@${MAKE} install_one_organism ORG=${BACT}		\
		ORG_DIR=${NCBI_DIR}/Bacteria/${BACT}

### Parse one organism
parse_organism:
	@echo "Parsing organism ${ORG}"
	${MAKE} install_one_organism INSTALL_TASK=parse


ENSEMBL_DIR=${RSAT}/downloads/ftp.ensembl.org/pub/current_worm/data/flatfiles/genbank
parse_organism_ensembl:
	parse-genbank.pl -v 1 -source ensembl -ext dat -i ${ENSEMBL_DIR} -org ${ORG}_ENSEMBL

install_organism_ensembl:
	${MAKE} install_one_organism  ORG=${ORG}_ENSEMBL OPT='-source ensembl -ensembl ${ENSEMBL_DIR}'

################################################################
#### Install all eukaryote genomes
#### Genomes are selected manually because NCBI directories are
#### a bit messy for eukaryotes.
EUKARYOTES=					\
	Saccharomyces_cerevisiae		\
	Schizosaccharomyces_pombe		\
	Encephalitozoon_cuniculi		\
	Plasmodium_falciparum			\
	Apis_mellifera				\
	Drosophila_melanogaster			\
	Arabidopsis_thaliana			\
	Caenorhabditis_elegans			\
	Gallus_gallus				\
	Homo_sapiens				\
	Mus_musculus				\
	Rattus_norvegicus			\
	Canis_familiaris			\
	Pan_troglodytes
install_all_eukaryotes:
	for org in ${EUKARYOTES} ; do \
		${MAKE} install_one_organism ORG=$${org} INSTALL_TASK=${INSTALL_TASK},clean; \
	done

FULL_ORG=${ORG}
LINK_DIR=${RSAT}/genome_installations
LINK_DIR_ORG=${LINK_DIR}/${FULL_ORG}
install_one_eukaryote:
	@mkdir -p ${LINK_DIR_ORG}
	@rm -rf ${LINK_DIR_ORG}/*
	@(cd ${LINK_DIR_ORG}; ln -s ${ORG_DIR}/CHR_*/*.gbk.gz .)
	@echo ${LINK_DIR_ORG}
	${MAKE} install_one_organism NCBI_DIR=${LINK_DIR} ORG=${FULL_ORG}
