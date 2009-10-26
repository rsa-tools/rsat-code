############################################################
#
# $Id: server.mk,v 1.22 2009/10/26 15:49:40 jvanheld Exp $
#
# Time-stamp: <2003-10-10 22:49:55 jvanheld>
#
############################################################

#RSAT=${HOME}/rsa-tools/
include ${RSAT}/makefiles/util.mk

GENBANK_DIR=/home/rsa/downloads/ftp.ncbi.nih.gov/genbank/genomes
NCBI_DIR=/home/rsa/downloads/ftp.ncbi.nih.gov/genomes

DATE = `date +%Y%m%d_%H%M%S`


#################################################################
# programs

WGET = wget -np -rNL 
MAKEFILE=${RSAT}/makefiles/server.mk
MAKE=nice -n 19 make -s -f ${MAKEFILE}
RSYNC_OPT = -ruptvl ${OPT}
SSH=-e 'ssh -x'
RSYNC = rsync ${RSYNC_OPT} ${SSH}

################################################################
# Mirrors
RSAT_SERVER = rsat@rsat.ulb.ac.be:rsa-tools
MAMAZE = jvanheld@164.15.109.52:rsa-tools
MERLIN = jvanheld@merlin.bigre.ulb.ac.be:rsa-tools
#FLYCHIP = jvanheld@flychip.org.uk:rsa-tools
CCG = jvanheld@kayab.ccg.unam.mx:rsa-tools
PORTPIN = rsat@portpin.univ-mrs.fr:rsa-tools
LIV = jvanheld@liv.bmc.uu.se:rsa-tools
TORONTO=jvanheld@ws03.ccb.sickkids.ca:rsa-tools
MILLIWAYS=jvanheld@milliways.bi.up.ac.za:rsa-tools
MIRROR_SERVERS = ${MAMAZE} ${MERLIN} ${LIV} ${TORONTO} ${CCG}  ${PORTPIN}
LOG_SERVERS= ${LIV} ${TORONTO} ${CCG} ${MILLIWAYS} ${RSAT_SERVER}  ${PORTPIN}

################################################################
## distribution
MEDICEL = root@grimsel.co.helsinki.fi:/work/programs/rsa-tools

MIRROR=${LIV}

################################################################
#### from brol to mirrors
################################################################
DIR=perl-scripts
DIRS=perl-scripts public_html doc
rsync_mirrors:
	@for mirror in ${MIRRORS} ; do					\
		${MAKE} rsync_one_mirror MIRROR=$${mirror} DIR=$${dir} ;	\
	done

rsync_one_mirror:
	@for dir in ${DIRS}; do							\
		${MAKE} rsync_dir  MIRROR=${MIRROR} DIR=$${dir} ;	\
	done

EXCLUDED=						\
	--exclude '*~'					\
	--exclude data					\
	--exclude tmp					\
	--exclude logs					\
	--exclude perl-scripts/lib/arch			\
	--exclude qd.pl				
rsync_dir:
	@echo "Synchronizing dir ${DIR} to mirror ${MIRROR}"
	${RSYNC} ${OPT} ${EXCLUDED}		\
		${DIR} ${MIRROR}/

DATA_EXCLUDED= --exclude 'Mus_*'		\
	--exclude 'Homo_*'			\
	--exclude 'Rattus_*'			\
	--exclude comparative_genomics		\
	--exclude upstream_calibrations

RSYNC_DATA_CMD=${RSYNC} ${DATA_EXCLUDED} \
	public_html/data  ${MIRROR}/public_html/ 
rsync_data:
	@for mirror in ${MIRRORS} ; do					\
		${MAKE} rsync_data_one_mirror MIRROR=$${mirror} ;	\
	done

rsync_data_one_mirror:
	@echo "Synchronizing data to mirror ${MIRROR}" 
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
	@for mirror in ${MIRRORS} ; do				\
		${RSYNC} archives/* $${mirror}/archives/ ;	\
	done

################################################################
#### from mirrors to brol
################################################################
rsync_logs:
	@for mirror in ${LOG_SERVERS} ; do					\
		echo "${RSYNC} $${mirror}/logs/log-file_* logs/" ;	\
		${RSYNC} $${mirror}/logs/log-file_* logs/ ;		\
	done

rsync_config:
	@for mirror in ${MIRROR_SERVERS} ; do \
		${RSYNC} $${mirror}/config/*.config config/ ;\
	done

FOLDERS=data pdf_files
from_ucmb:
	@for folder in ${FOLDERS}; do \
		${RSYNC} jvanheld@${PAULUS}:rsa-tools/$${folder} . ; \
	done

FOLDERS=data 
from_cifn:
	@for folder in ${FOLDERS}; do \
		${RSYNC} jvanheld@${CCG}:rsa-tools/$${folder}/* ./$${folder} ; \
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

# #### Obsolete
# #### Instal the GD library
# GD_DISTRIB_DIR=stein.cshl.org/WWW/software/GD/
# GD_DISTRIB= http://${GD_DISTRIB_DIR}/GD.pm.tar.gz
# get_gd:
# 	mkdir -p lib-sources
# 	(cd lib-sources; \
# 	${WGET} ${GD_DISTRIB})


# GD_VERSION=2.06
# uncompress_gd:
# 	(cd lib-sources/${GD_DISTRIB_DIR};		\
# 	gunzip -c GD.pm.tar.gz | tar -xpf - ;		\

# install_gd:
# 	(cd lib-sources/${GD_DISTRIB_DIR}/GD-${GD_VERSION} ;	\
# 	perl Makefile.PL INSTALLDIRS=site			\
# 		INSTALLSITELIB=${RSA}/extlib			\
# 		INSTALLSITEARCH=${RSA}/extlib/arch ;		\
# 	make ;							\
# 	make install ;						\
# 	)

################################################################
#### Install programs from third parties
APP_DIR=${RSA}/applications
PROGRAM=consensus
PROGRAM_DIR=${APP_DIR}/${PROGRAM}
PROGRAM_ARCHIVE=`ls -1t ${RSA}/app-sources/${PROGRAM}* | head -1`
uncompress_program:
	@echo installing ${PROGRAM_ARCHIVE} in dir ${PROGRAM_DIR}
	@mkdir -p ${PROGRAM_DIR}
	(cd  ${PROGRAM_DIR} ;		\
	gunzip -c ${PROGRAM_ARCHIVE} | tar -xvf - )

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
	(cd ${GIBBS_DIR}; ./compile)
	(cd bin; ln -fs ${GIBBS_DIR}/gibbs ./gibbs)

install_ext_apps:
	${MAKE} install_consensus
	${MAKE} install_patser
	${MAKE} install_gibbs

################################################################
#### clean temporary directory
CLEAN_DATE=3
clean_tmp:
	@echo `date` Cleaning temporary directory on `hostname` 
	@echo "Before cleaning	" `du -sk ${RSAT}/public_html/tmp`
	touch ${RSAT}/public_html/tmp
	find ${RSAT}/public_html/tmp/ -mtime +${CLEAN_DATE} -type f -exec rm -rf {} \;	
	@echo "After cleaning	" `du -sk ${RSAT}/public_html/tmp`
	@echo "Disk free" 
	@df `pwd`   
