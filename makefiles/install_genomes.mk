############################################################
#
# $Id: install_genomes.mk,v 1.9 2005/01/18 15:25:08 jvanheld Exp $
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
INSTALL_TASK=allup,clean,config,dyads,ncf,intergenic_freq,oligos,parse,start_stop,upstream_freq
INSTALL_CMD=install-organism -v ${V}		\
		-org ${ORG}			\
		-task  ${INSTALL_TASK}
install_one_organism:
	@echo "install log	${INSTALL_LOG}"
	@echo "Parsing organism ${ORG}" 
	${MAKE} my_command MY_COMMAND="${INSTALL_CMD}"

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
	@${MAKE} install_one_organism ORG=${BACT}		\
		ORG_DIR=${NCBI_DIR}/Bacteria/${BACT}

### Parse one organism
parse_organism:
	@echo "Parsing organism ${ORG}"
	${MAKE} install_one_organism INSTALL_TASK=parse

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
	for org in ${EUKARYOTES} ; do ${MAKE} install_one_organism	\
		ORG=$${org} ; done


################################################################
#
# Clean obsolete genome files (those with uppercases)
ORGANISMS=`ls -d1 data/genomes/*_* | xargs basename | xargs`
clean_genomes:
	${MAKE} iterate_organisms ORG_TASK=clean_one_genome

ORG=Mycoplasma_genitalium
clean_one_genome:
	@echo "Cleaning genome ${ORG}"
	@echo "Before cleaning"
	ls -ltr ${RSAT}/data/genomes/${ORG}/genome/
	rm -f ${RSAT}/data/genomes/${ORG}/genome/Feature*
	rm -f ${RSAT}/data/genomes/${ORG}/genome/Gene*
	rm -f ${RSAT}/data/genomes/${ORG}/genome/CDS*.ctl
	rm -f ${RSAT}/data/genomes/${ORG}/genome/CDS*
	rm -f ${RSAT}/data/genomes/${ORG}/genome/mRNA*
	rm -f ${RSAT}/data/genomes/${ORG}/genome/tRNA*.ctl
	rm -f ${RSAT}/data/genomes/${ORG}/genome/tRNA*
	rm -f ${RSAT}/data/genomes/${ORG}/genome/rRNA*
	rm -f ${RSAT}/data/genomes/${ORG}/genome/*RNA*
	rm -f ${RSAT}/data/genomes/${ORG}/genome/Organism*
	rm -f ${RSAT}/data/genomes/${ORG}/genome/Contig*
	rm -f ${RSAT}/data/genomes/${ORG}/genome/misc_RNA*
	rm -f ${RSAT}/data/genomes/${ORG}/genome/Source*
	rm -f ${RSAT}/data/genomes/${ORG}/genome/*.wc
	rm -f ${RSAT}/data/genomes/${ORG}/genome/*.wc.gz
	@echo "After cleaning"
	ls -ltr ${RSAT}/data/genomes/${ORG}/genome/
	${MAKE} clean_sql

################################################################
# Clean obsolete SQL files (before splitting into subfolders for
# mysql, postgresql and oracle)
clean_sql:
	@echo "Cleaning SQL files ${ORG}"
	@echo "Before cleaning"
	ls -ltr ${RSAT}/data/genomes/${ORG}/genome/sql_scripts/
	rm -f ${RSAT}/data/genomes/${ORG}/genome/sql_scripts/makefile
	rm -f ${RSAT}/data/genomes/${ORG}/genome/sql_scripts/*.sql
	rm -f ${RSAT}/data/genomes/${ORG}/genome/sql_scripts/*.ctl
	rm -f ${RSAT}/data/genomes/${ORG}/genome/sql_scripts/*.mk
	@echo "After cleaning"
	ls -ltr ${RSAT}/data/genomes/${ORG}/genome/sql_scripts/


