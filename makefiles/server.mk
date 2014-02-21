############################################################
#
# $Id: server.mk,v 1.41 2012/08/05 13:13:10 rsat Exp $
#
# Time-stamp: <2003-10-10 22:49:55 jvanheld>
#
############################################################

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
#BIGRE=rsat@rsat.ulb.ac.be:rsa-tools
#WWWSUP=rsat@wwwsup.scmbb.ulb.ac.be:rsa-tools
#MAMAZE=rsat@mamaze.ulb.ac.be:rsa-tools


BIGRE=rsat@${MERLIN}:/rsat_servers/rsat
WWWSUP=rsat@${MERLIN}:/rsat_servers/wwwsup
MAMAZE=rsat@${MERLIN}:/rsat_servers/mamaze

#CCG=jvanheld@itzamna.ccg.unam.mx:rsa-tools
CCG=rsat@itzamna.ccg.unam.mx:rsa-tools
TAGC=jvanheld@pedagogix-tagc.univ-mrs.fr:rsa-tools
UPPSALA=jvanheld@bongcam1.hgen.slu.se:rsa-tools
PRETORIA=jvanheld@anjie.bi.up.ac.za:.
LOG_SERVERS= ${MAMAZE} ${WWWSUP}  ${CCG} ${UPPSALA} ${PRETORIA} ${BIGRE} 

################################################################
## OLD SERVERS, NOT MAINTAINED ANYMORE
#FLYCHIP=jvanheld@flychip.org.uk:rsa-tools
#TORONTO=jvanheld@ws03.ccb.sickkids.ca:rsa-tools
#PRETORIA=jvanheld@milliways.bi.up.ac.za:rsa-tools

################################################################
## distribution
MEDICEL=root@grimsel.co.helsinki.fi:/work/programs/rsa-tools

MIRROR=${UPPSALA}

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
	@echo
	@echo "Synchronizing mirror ${MIRROR}"
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
	@for mirror in ${LOG_SERVERS} ; do				\
		echo ;							\
		echo "Synchronizing mirror	$${mirror}";		\
		echo "${RSYNC} $${mirror}/logs/log-file_* logs/" ;	\
		${RSYNC} $${mirror}/logs/log-file_* logs/ ;		\
	done
	rsync -ruptvl -e 'ssh -p 24222'  jvanheld@pedagogix-tagc.univ-mrs.fr:rsa-tools/logs/log-file_* logs/


################################################################
## Clean temporary directory
CLEAN_LIMIT=3
clean_tmp:
	@echo "Cleaning temporary directory	`hostname` ${RSAT}/public_html/tmp/"
	@echo
	@date "+%Y/%m/%d %H:%M:%S"
	@echo "Free disk before cleaning" 
	@df -h ${RSAT}/public_html/tmp/
	@echo
	@date "+%Y/%m/%d %H:%M:%S"
	@echo "Measuring disk usage before cleaning"
	@echo "Before cleaning	" `du -sh public_html/tmp`
	@touch ${RSAT}/public_html/tmp/
	@echo
	@date "+%Y/%m/%d %H:%M:%S"
	@echo "Removing all files older than ${CLEAN_LIMIT} days"
	find ${RSAT}/public_html/tmp/ -mtime +${CLEAN_LIMIT} -type f -exec $(SUDO) rm -f {} \;	
	@echo
	@date "+%Y/%m/%d %H:%M:%S"
	@echo "Measuring disk usage after cleaning"
	@echo "After cleaning	" `du -sh public_html/tmp`
	@echo "Cleaned temporary directory" | mail -s 'cleaning tmp' ${RSAT_ADMIN_EMAIL}
	@echo
	@date "+%Y/%m/%d %H:%M:%S"
	@echo "Free disk after cleaning" 
	@df -h ${RSAT}/public_html/tmp/


################################################################
## Detect web spammers (2012). Several IP addresses are repeatedly sending Web
## spam to the Web interfaces of gene-info.cgi and convert-matrix.cgi.
YEAR=`date +%Y`
DENIED_IP_FILE=denied_IP_addresses_${RSAT_SITE}_${YEAR}.tab
FORM_DENIAL_THRESHOLD=500
TAG_DENIAL_THRESHOLD=5
ATTACKED_FORMS=gene-info.cgi convert-matrix.cgi RSAT_home.cgi
denied_ips:
	@echo 
	@echo "Detecting suspicious IP addresses (Web spammers)"
	@echo "	frequent HTML tags in queries (> ${TAG_DENIAL_THRESHOLD})"
	@cut -f 3 ${RSAT}/logs/web_attacks_log_${RSAT_SITE}_${YEAR}_*.txt \
			| perl -pe 's|\@||' \
			| perl -pe 's| \(\)||' \
			| contingency-table  -col1 1 -col2 1 -margin \
			| grep -v '^;' \
			| grep -v '^#' \
			| cut -f 1,2 \
			| awk '$$2 > ${TAG_DENIAL_THRESHOLD} {print $$1"\t"$$2"\tHTML_tags"}' \
			| sort > ${DENIED_IP_FILE}
	@for form in ${ATTACKED_FORMS}; do \
		echo "	abusive use of form $${form}  (> ${FORM_DENIAL_THRESHOLD})"; \
		${MAKE}  _denied_ips_one_script ATTACKED_FORM=$${form} ; \
	done
	@echo "	${DENIED_IP_FILE}"
	@wc -l ${DENIED_IP_FILE}

ATTACKED_FORM=gene-info.cgi
_denied_ips_one_script:
	grep ${ATTACKED_FORM} logs/log-file_${RSAT_SITE}_${YEAR}_*  \
		| cut -f 3 \
		| perl -pe 's|\@||' \
		| perl -pe 's| \(\)||' \
		| contingency-table  -col1 1 -col2 1 -margin \
		| grep -v '^;' \
		| grep -v '^#' \
		| cut -f 1,2 \
		| awk '$$2 > ${FORM_DENIAL_THRESHOLD} {print $$1"\t"$$2"\t${ATTACKED_FORM}"}' \
		| sort >> ${DENIED_IP_FILE} ; \

