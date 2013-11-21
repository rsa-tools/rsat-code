############################################################
#
# $Id: mirror.mk,v 1.54 2012/08/04 16:17:23 jvanheld Exp $
#
# Time-stamp: <2003-10-01 12:05:45 jvanheld>
#
############################################################

include ${RSAT}/makefiles/util.mk
MAKEFILE=${RSAT}/makefiles/mirror.mk

################################################################
#
# Server
RSAT_SERVER=merlin.bigre.ulb.ac.be
RSAT_SERVER_DIR=/rsat_servers/rsat
RSAT_SERVER_LOGIN=rsat
SERVER=${RSAT_SERVER_LOGIN}@${RSAT_SERVER}:${RSAT_SERVER_DIR}

DATE = `date +%Y%m%d_%H%M%S`


#################################################################
# programs
OPT=
MAKE=make -sk -f ${MAKEFILE}
RSYNC_OPT = -ruptvl ${OPT} --exclude '*~' --exclude jobs
SSH=-e ssh
RSYNC = rsync ${RSYNC_OPT} ${SSH}



################################################################
#### from local machine to servers
################################################################
DIR=perl-scripts
DIRS=perl-scripts public_html doc

rsync_to_server: script_to_server pub_to_server

scripts_to_server:
	echo "Synchronizing perl-scripts to server ${SERVER}"
	${RSYNC} --exclude perllib --exclude perl-scripts/lib/arch --exclude qd.pl perl-scripts ${SERVER}/

pub_to_server:
	echo "Synchronizing public_html to server ${SERVER}"
	${RSYNC} --exclude logs --exclude tmp --exclude data public_html ${SERVER}/  

data_to_server:
	echo "Synchronizing data to server ${SERVER}"
	${RSYNC} ${EXCLUDED_FILES}  public_html/data/* ${SERVER}/public_html/data/

genomes_to_server:
	echo "Synchronizing genomes to server ${SERVER}"
	${RSYNC} ${EXCLUDED_FILES} public_html/data/taxon_frequencies ${SERVER}/public_html/data/
	${RSYNC} ${EXCLUDED_FILES} public_html/data/genomes ${SERVER}/public_html/data/
	${RSYNC}  public_html/data/supported*.tab ${SERVER}/public_html/data/
	${RSYNC}  public_html/data/supported*.pl ${SERVER}/public_html/data/

doc_to_server:
	${MAKE} dir_to_server DIR=doc


dir_to_server:
	echo "Synchronizing dir ${DIR} to server ${SERVER}" 
	${RSYNC} ${DIR} ${SERVER}/ 


################################################################
#### From server to local machine
################################################################
all_from_server: scripts_from_server pub_from_server data_from_server supported_from_server logs_from_server

DIR=doc
TARGET_DIR=${RSAT}/
RSYNC_FROM_SERVER_CMD=${RSYNC} ${SERVER}/${DIR} ${TARGET_DIR}
dir_from_server:
	@echo ${RSYNC_FROM_SERVER_CMD}
	${RSYNC_FROM_SERVER_CMD}

doc_from_server:
	${MAKE} dir_from_server DIR=doc

logs_from_server:
	${MAKE} dir_from_server DIR='public_html/logs' RSYNC_OPT='-ruptvl ${OPT}' TARGET_DIR=${RSAT}/public_html/

scripts_from_server:
	${MAKE} dir_from_server DIR=perl-scripts

pub_from_server:
	${RSYNC}								\
		--exclude data							\
		--exclude logs							\
		--exclude tmp							\
		${SERVER}/public_html ${RSAT}/

EXCLUDED_GENOMES=					\
		--exclude phages 			\
		--exclude Bos_*				\
		--exclude Danio_rerio*			\
		--exclude Mus_musculus*			\
		--exclude Gallus_gallus*		\
		--exclude Canis_familiaris*		\
		--exclude Pan_troglodytes*		\
		--exclude Rattus_norvegicus*		\
		--exclude Danio_rerio*			\
		--exclude Drosophila_melanogaster*	\
		--exclude Anopheles_gambiae*		\
		--exclude Caenorhabditis_elegans*	\
		--exclude Tetraodon_nigroviridis*	\
		--exclude Oryzias_latipes*		\
		--exclude Arabidopsis_thaliana		\
		--exclude Apis_mellifera		\
		--exclude Homo_sapiens*			\
		--exclude *_EnsEMBL*

#EXCLUDED_BLAST= \
#		--exclude 'blastdb'		\
#		--exclude 'blast_hits'		

EXCLUDED_FILES=					\
		--exclude '*.wc'		\
		--exclude '*.wc.gz'		
#		--exclude '*.fasta'		\
#		--exclude '*.fasta.gz'

EXCLUDED_DIRS=					\
		--exclude 'published_data'	\
		--exclude embl_genomes		\
		--exclude previous_version	\
		--exclude tmp			\
		--exclude upstream_calibrations	\
		--exclude comparative_genomics

EXCLUDED=${EXCLUDED_GENOMES} ${EXCLUDED_DIRS} ${EXCLUDED_FILES} ${EXCLUDED_BLAST}
data_from_server:
	${RSYNC} ${EXCLUDED} ${OPT}							\
		${SERVER}/public_html/data/*	\
		${RSAT}/public_html/data/

## Synchronize the list of supported organisms from the server to the mirror
supported_from_server:
	rsync -ruptvl -e ssh rsat@rsat.ulb.ac.be:rsat/data/supported'*' data/

################################################################
## Use wget to synchronize this mirror from the main server
RSAT_HTTP=http://rsat.ulb.ac.be/rsat/
ORG=Mycoplasma_genitalium
wget_one_org:
	wget -rNL -P data/genomes/ ${RSAT_HTTP}/data/genomes/${ORG} --cut-dirs=3 -nH


# ################################################################
# #### Server in finland
# ################################################################
# ORGS=Saccharomyces_cerevisiae Escherichia_coli_K12 Bacillus_subtilis
# medicel:
# 	${RSYNC} config/medicel.config ${MEDICEL}/config/
# 	${RSYNC} doc/*.pdf ${MEDICEL}/doc/
# 	rsync ${SSH} -ruptvL distrib/* ${MEDICEL}/perl-scripts
# 	for org in ${ORGS}; do				\
# 		${RSYNC} data/$${org} ${MEDICEL}/data/;	\
# 	done


################################################################
## Clean temporary directory

## Days
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
	find ${RSAT}/public_html/tmp/ -mtime +${CLEAN_LIMIT} -type f -exec rm -f {} \;	
	@echo
	@date "+%Y/%m/%d %H:%M:%S"
	@echo "Measuring disk usage after cleaning"
	@echo "After cleaning	" `du -sh public_html/tmp`
	@echo "Cleaned temporary directory" | mail -s 'cleaning tmp' ${RSAT_ADMIN_EMAIL}
	@echo
	@date "+%Y/%m/%d %H:%M:%S"
	@echo "Free disk after cleaning" 
	@df -h ${RSAT}/public_html/tmp/
