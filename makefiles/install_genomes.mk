############################################################
#
# $Id: install_genomes.mk,v 1.1 2003/12/12 11:29:12 jvanheld Exp $
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
	@make install_organism INSTALL_TASK=parse

################################################################
#### Install all eukaryote genomes
#### Genomes are selected manually because NCBI directories are
#### a bit messy for eukaryotes.
EUKARYOTES=					\
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
ORGANISM_DIR=${GENBANK_DIR}/${ORGANISM}
INSTALL_TASK=allup,clean,config,dyads,ncf,intergenic_freq,oligos,parse,start_stop,upstream_freq
install_organism:
	@echo "install log	${INSTALL_LOG}"
	@echo "Parsing organism ${ORGANISM}" 
	install-organism -v 1								\
		-org ${ORGANISM}							\
		-task  ${INSTALL_TASK}

# ### Install S.pombe genome (obsolete)
# POMBE_DIR=/win/databases/downloads/ftp.sanger.ac.uk/pub/yeast/Pombe/CONTIGS/
# install_pombe:
# 	echo "Parsing organism Schizosaccharomyces pombe" ;
# #	parse-embl.pl -i ${POMBE_DIR} -org 'Schizosaccharomyces pombe' -v 1
# 	install-organism -v 1											\
# 		-org Schizosaccharomyces_pombe									\
# 		-features ${RSA}/data/Schizosaccharomyces_pombe/genome/Gene_Schizosaccharomyces_pombe.tab	\
# 		-genome ${RSA}/data/genome/Contigs_Schizosaccharomyces_pombe.txt				\
# 		-format filelist										\
# 		-source genbank											\
# 		-step config -step start_stop -step ncf -step oligos -step dyads;
