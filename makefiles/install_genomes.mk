############################################################
#
# $Id: install_genomes.mk,v 1.4 2004/06/12 17:33:36 jvanheld Exp $
#
# Time-stamp: <2003-10-10 22:49:55 jvanheld>
#
############################################################

DATE = `date +%Y%m%d_%H%M%S`

################################################################
#### Directories
RSAT=${HOME}/rsa-tools/
GENBANK_DIR=${RSAT}/downloads/ftp.ncbi.nih.gov/genbank/genomes
NCBI_DIR=${RSAT}/downloads/ftp.ncbi.nih.gov/genomes

#################################################################
#### Programs

WGET = wget -np -rNL 
MAKEFILE=${RSAT}/makefiles/install_genomes.mk
MAKE=nice -n 19 make -s -f ${MAKEFILE}
RSYNC_OPT = -ruptvl ${OPT}
SSH=-e 'ssh -x'
RSYNC = rsync ${RSYNC_OPT} ${SSH}

V=1 

################################################################
### Targets

#### list the supported tasks
usage:
	@echo "usage: make [-OPT='options'] target"
	@echo "implemented targets"
	@perl -ne 'if (/^(\S+):/){ print "\t$$1\n"}' ${MAKEFILE}

#BACTERIA = `cat TO_INSTALL.txt| sort -ru | xargs `
#BACTERIA =					\
#	Clostridium_perfringens			\
#	Pyrobaculum_aerophilum			\
#	Pyrococcus_furiosus

### Pet bacteria for quick testing
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
	@${MAKE} install_organism ORGANISM=${BACT}		\
		ORGANISM_DIR=${NCBI_DIR}/Bacteria/${BACT}

### Parse one organism
parse_organism:
	@${MAKE} install_organism INSTALL_TASK=parse

################################################################
#### Install all eukaryote genomes
#### Genomes are selected manually because NCBI directories are
#### a bit messy for eukaryotes.
EUKARYOTES=					\
	Apis_mellifera				\
	Arabidopsis_thaliana			\
	Caenorhabditis_elegans			\
	Drosophila_melanogaster			\
	Encephalitozoon_cuniculi		\
	Homo_sapiens				\
	Mus_musculus				\
	Plasmodium_falciparum			\
	Rattus_norvegicus			\
	Saccharomyces_cerevisiae		\
	Schizosaccharomyces_pombe 
install_all_eukaryotes:
	for org in ${EUKARYOTES} ; do				\
		${MAKE} install_one_organism ORGANISM=$${org} ;	\
	done

### Install one organism
ORGANISM=Arabidopsis_thaliana
ORGANISM_DIR=${NCBI_DIR}/${ORGANISM}
INSTALL_TASK=allup,clean,config,dyads,ncf,intergenic_freq,oligos,parse,start_stop,upstream_freq
install_organism:
	@echo "install log	${INSTALL_LOG}"
	@echo "Parsing organism ${ORGANISM}" 
	install-organism -v ${V}								\
		-org ${ORGANISM}							\
		-task  ${INSTALL_TASK}

