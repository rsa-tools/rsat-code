############################################################
#
# $Id: install_rsat.mk,v 1.1 2003/10/29 00:35:14 jvanheld Exp $
#
# Time-stamp: <2003-05-23 09:36:00 jvanheld>
#
############################################################

RSA=${HOME}/rsa-tools/
#RSA=${HOME}/rsa/rsa-tools/
GENBANK_DIR=/home/rsa/downloads/ftp.ncbi.nih.gov/genbank/genomes
NCBI_DIR=/home/rsa/downloads/ftp.ncbi.nih.gov/genomes

DATE = `date +%Y%m%d_%H%M%S`


#################################################################
# programs

WGET = wget -np -rNL 
MAKEFILE=${RSA}/makefiles/install_rsat.mk
MAKE=nice -n 19 make -s -f ${MAKEFILE}
RSYNC_OPT = -ruptvl ${OPT}
SSH=-e 'ssh -x'
RSYNC = rsync ${RSYNC_OPT} ${SSH}

################################################################
#
# Servers

CIFN = jvanheld@embnet.cifn.unam.mx:rsa-tools/
PAULUS = jvanheld@paulus.ulb.ac.be:rsa-tools/
UPPSALA = jvanheld@bioinformatics.bmc.uu.se:rsa-tools
LIV = jvanheld@liv.bmc.uu.se:rsa-tools
RSAT = jvanheld@164.15.61.35:/home/rsa/rsa-tools
MEDICEL = root@grimsel.co.helsinki.fi:/work/programs/rsa-tools
SERVERS = ${CIFN} ${UPPSALA} ${LIV}

#SERVER=${RSAT}
SERVER=${UPPSALA}

### tags
usage:
	@echo "usage: make [-OPT='options'] target"
	@echo "implemented targets"
	@perl -ne 'if (/^(\S+):/){ print "\t$$1\n"}' ${MAKEFILE}

################################################################
#### from brol to servers
################################################################
DIR=perl-scripts
DIRS=perl-scripts public_html doc
rsync_servers:
	for server in ${SERVERS} ; do					\
		${MAKE} rsync_server SERVER=$${server} DIR=$${dir} ;	\
	done

rsync_server:
	@for dir in ${DIRS}; do							\
		${MAKE} rsync_dir  SERVER=${SERVER} DIR=$${dir} ;	\
	done

rsync_dir:
	echo "Synchronizing dir ${DIR} to server ${SERVER}"
	${RSYNC}  ${OPT} --exclude data --exclude tmp --exclude logs --exclude perl-scripts/lib/arch --exclude qd.pl ${DIR} ${SERVER}/

RSYNC_DATA_CMD=${RSYNC} --exclude 'Mus_musculus*' --exclude 'Homo_sapiens*' public_html/data ${SERVER}/public_html/ 
rsync_data:
	@for server in ${SERVERS} ; do					\
		${MAKE} rsync_data_one_server SERVER=$${server} ;	\
	done

rsync_data_one_server:
	@echo "Synchronizing data to server ${SERVER}" 
	@echo ${RSYNC_DATA_CMD} ;			
	${RSYNC_DATA_CMD};				


ORGS=Saccharomyces_cerevisiae Escherichia_coli_K12 Bacillus_subtilis
medicel:
	${RSYNC} config/medicel.config ${MEDICEL}/config/
	${RSYNC} doc/*.pdf ${MEDICEL}/doc/
	rsync ${SSH} -ruptvL distrib/* ${MEDICEL}/perl-scripts
	for org in ${ORGS}; do				\
		${RSYNC} data/$${org} ${MEDICEL}/data/;	\
	done

rsync_archives:
	@for server in ${SERVERS} ; do				\
		${RSYNC} archives/* $${server}/archives/ ;	\
	done

################################################################
#### from servers to brol
################################################################
rsync_logs:
	@for server in ${SERVERS} ; do					\
		echo "${RSYNC} $${server}/logs/log-file_* logs/" ;	\
		${RSYNC} $${server}/logs/log-file_* logs/ ;		\
	done

rsync_config:
	@for server in ${SERVERS} ; do \
		${RSYNC} $${server}/config/*.config config/ ;\
	done

FOLDERS=data pdf_files
from_ucmb:
	@for folder in ${FOLDERS}; do \
		${RSYNC} jvanheld@${PAULUS}:rsa-tools/$${folder} . ; \
	done

FOLDERS=data 
from_cifn:
	@for folder in ${FOLDERS}; do \
		${RSYNC} jvanheld@${CIFN}:rsa-tools/$${folder}/* ./$${folder} ; \
	done

SRC=perl-scripts
COMPIL=compil/
PROGRAMS=	\
	oligo-analysis	\
	retrieve-seq	\
	dyad-analysis	\
	dna-pattern	\
	orf-info
LIBRARIES=\
	RSA.stat.lib 	\
#	RSA.classes	\
#	RSA.seq.lib	\
#	RSA.cgi.lib	\
#	RSA.lib 	
compile:
	@mkdir -p ${COMPIL}/lib
	@mkdir -p ${COMPIL}/bin
	@(cd  ${COMPIL}/bin; ln -fs ../lib)
	@cp -f config/default.config ${COMPIL}/RSA.config

	@for lb in ${LIBRARIES}; do \
		echo "compiling library $${lb}"; \
		cp -f ${SRC}/lib/$${lb} ${COMPIL}/lib/$${lb}.pl ; \
		(cd ${COMPIL}/lib; pwd; perlcc $${lb}.pl && rm -f $${lb}.pl); \
	done

	@for pgm in ${PROGRAMS}; do \
		echo "compiling program $${pgm}"; \
		cp -f ${SRC}/$${pgm} ${COMPIL}/bin/$${pgm}.pl ; \
		(cd ${COMPIL}/bin; pwd; perlcc $${pgm}.pl && rm -f $${pgm}.pl); \
	dgone

#BACTERIA = `ls -1 genbank/genomes/Bacteria/ | grep _ | xargs`
#BACTERIAS = `ls -1 ${GENBANK_DIR}/Bacteria | grep -v Mesorhizobium_loti | grep _ | tail -37 | xargs `
#BACTERIA = `ls -1 ${NCBI_DIR}/Bacteria | grep _ | sort -ru | xargs `
BACTERIA = `cat TO_INSTALL.txt| sort -ru | xargs `
#BACTERIAS = Ralstonia_solanacearum
#BACTERIA =					\
#	Clostridium_perfringens			\
#	Pyrobaculum_aerophilum			\
#	Pyrococcus_furiosus
BACT=Mycoplasma_genitalium

list_bacteria:
	@echo "Bacteria to install	${BACTERIA}"

install_all_bacteria:
	@for bact in ${BACTERIA}; do				\
		${MAKE} install_one_bacteria BACT=$${bact} ;	\
	done

install_one_bacteria:
	@echo
	@echo "${DATE}	Installing bacteria ${BACT}"
	@${MAKE} install_organism ORGANISM=${BACT}		\
		ORGANISM_DIR=${NCBI_DIR}/Bacteria/${BACT}

#ORGANISM=Plasmodium_falciparum
ORGANISM=Homo_sapiens
ORGANISM_DIR=${GENBANK_DIR}/${ORGANISM}
INSTALL_TASK=allup,clean,config,dyads,ncf,intergenic_freq,oligos,parse,start_stop,upstream_freq
install_organism:
	@echo "install log	${INSTALL_LOG}"
	echo "Parsing organism ${ORGANISM}" 
	install-organism -v 1								\
		-org ${ORGANISM}							\
		-task  ${INSTALL_TASK}

POMBE_DIR=/win/databases/downloads/ftp.sanger.ac.uk/pub/yeast/Pombe/CONTIGS/
install_pombe:
	echo "Parsing organism Schizosaccharomyces pombe" ;
#	parse-embl.pl -i ${POMBE_DIR} -org 'Schizosaccharomyces pombe' -v 1
	install-organism -v 1											\
		-org Schizosaccharomyces_pombe									\
		-features ${RSA}/data/Schizosaccharomyces_pombe/genome/Gene_Schizosaccharomyces_pombe.tab	\
		-genome ${RSA}/data/genome/Contigs_Schizosaccharomyces_pombe.txt				\
		-format filelist										\
		-source genbank											\
		-step config -step start_stop -step ncf -step oligos -step dyads;

GD_DISTRIB_DIR=stein.cshl.org/WWW/software/GD/
GD_DISTRIB= http://${GD_DISTRIB_DIR}/GD.pm.tar.gz
get_gd:
	mkdir -p lib-sources
	(cd lib-sources; \
	${WGET} ${GD_DISTRIB})


GD_VERSION=2.06
uncompress_gd:
	(cd lib-sources/${GD_DISTRIB_DIR};		\
	gunzip -c GD.pm.tar.gz | tar -xpf - ;		\

install_gd:
	(cd lib-sources/${GD_DISTRIB_DIR}/GD-${GD_VERSION} ;	\
	perl Makefile.PL INSTALLDIRS=site			\
		INSTALLSITELIB=${RSA}/extlib			\
		INSTALLSITEARCH=${RSA}/extlib/arch ;		\
	make ;							\
	make install ;						\
	)


APP_DIR=${RSA}/applications
PROGRAM=consensus
PROGRAM_DIR=${APP_DIR}/${PROGRAM}
PROGRAM_ARCHIVE=`ls -1t ${RSA}/app_sources/${PROGRAM}* | head -1`
uncompress_program:
	@echo installing ${PROGRAM_ARCHIVE} in dir ${PROGRAM_DIR}
	@mkdir -p ${PROGRAM_DIR}
	@(cd ${PROGRAM_DIR} ;				\
	gunzip -c ${PROGRAM_ARCHIVE} | tar -xf - )

INSTALLED_PROGRAM=`ls -1t ${APP_DIR}/${PROGRAM}/${PROGRAM}*`
install_program:
	(cd ${PROGRAM_DIR}; make ${INSTALL_OPT})
	(cd bin; ln -fs ${INSTALLED_PROGRAM} ./${PROGRAM})

install_patser:
	${MAKE} uncompress_program PROGRAM=patser
	${MAKE} install_program PROGRAM=patser

install_consensus:
	${MAKE} uncompress_program PROGRAM=consensus
	${MAKE} install_program PROGRAM=consensus INSTALL_OPT='CPPFLAGS=""'

GIBBS_DIR=${APP_DIR}/gibbs/gibbs9_95
install_gibbs:
	${MAKE} uncompress_program PROGRAM=gibbs
	(cd ${GIBBS_DIR}; ./compile; cd ${GIBBS_DIR}/code; make clean)
	(cd bin; ln -fs ${GIBBS_DIR}/gibbs ./gibbs)

install_ext_apps:
	${MAKE} install_consensus
	${MAKE} install_patser
	${MAKE} install_gibbs
