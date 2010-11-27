############################################################
#
# $Id: server.mk,v 1.26 2010/11/27 16:26:10 rsat Exp $
#
# Time-stamp: <2003-10-10 22:49:55 jvanheld>
#
############################################################

#RSAT=${HOME}/rsa-tools/
include ${RSAT}/makefiles/util.mk
MAKEFILE=${RSAT}/makefiles/server.mk

GENBANK_DIR=/home/rsa/downloads/ftp.ncbi.nih.gov/genbank/genomes
NCBI_DIR=/home/rsa/downloads/ftp.ncbi.nih.gov/genomes

DATE = `date +%Y%m%d_%H%M%S`


#################################################################
# programs

WGET = wget -np -rNL 
MAKE=nice -n 19 make -s -f ${MAKEFILE}
RSYNC_OPT = -ruptvl ${OPT}
SSH=-e 'ssh -x'
RSYNC = rsync ${RSYNC_OPT} ${SSH}

################################################################
# Mirrors
BIGRE = rsat@rsat.ulb.ac.be:rsa-tools
BIGRE2 = rsat@wwwsup.scmbb.ulb.ac.be:rsa-tools
MAMAZE = jvanheld@164.15.109.52:rsa-tools
MERLIN = jvanheld@merlin.bigre.ulb.ac.be:rsa-tools
#FLYCHIP = jvanheld@flychip.org.uk:rsa-tools
CCG = jvanheld@kayab.ccg.unam.mx:rsa-tools
TAGC = jvanheld@139.124.66.43:rsa-tools
LIV = jvanheld@liv.bmc.uu.se:rsa-tools
TORONTO=jvanheld@ws03.ccb.sickkids.ca:rsa-tools
MILLIWAYS=jvanheld@milliways.bi.up.ac.za:rsa-tools
MIRROR_SERVERS = ${MAMAZE} ${MERLIN} ${LIV} ${TORONTO} ${CCG} 
LOG_SERVERS= ${BIGRE2} ${LIV} ${CCG} ${MILLIWAYS} ${BIGRE} ${TORONTO} 

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
	rsync -ruptvl -e 'ssh -p 24222'  jvanheld@139.124.66.43:rsa-tools/logs/log-file_* logs/
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
