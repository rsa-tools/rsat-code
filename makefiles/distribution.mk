################################################################
## Prepare a release version of rsa-tools
## distribution.mk makefile
## Usage: make -f distribution.mk

include ${RSAT}/makefiles/util.mk

MAKEFILE=${RSAT}/makefiles/distribution.mk
MAKE = make -sk -f ${MAKEFILE}

## Archive file
DATE=`date +%Y%m%d`
ARCHIVE_PREFIX=rsa-tools_${DATE}
ARCHIVE=rsa-tools/${ARCHIVE_PREFIX}

## Archive with tar
#TAR_EXCLUDE=-X CVS '*~' 
TAR_CREATE =tar ${TAR_EXCLUDE} -cpvf ${ARCHIVE}.tar rsa-tools/RSA.config.default 
TAR =tar ${TAR_EXCLUDE} -rpvf ${ARCHIVE}.tar 


## Manuals

manuals:
	(cd doc/manuals; make fullclean; make ig; make ug; make tex_clean)

## Archive with zip
ZIP_EXCLUDE=-x CVS '*~' 
ZIP =zip -ry ${ARCHIVE}.zip 

POST_CMD=


DISTRIB_FILES= rsa-tools/perl-scripts					\
	rsa-tools/RSA.config.default					\
	rsa-tools/public_html/data/supported_organisms_template.txt	\
	rsa-tools/public_html/data/supported_organisms.pl		\
	rsa-tools/makefiles						\
	rsa-tools/doc/manuals/*.pdf
fill_archive:
	(cd ..;						\
	for f in ${DISTRIB_FILES}; do			\
		${MAKE} add_one_file FILE=$${f};	\
	done)
	@echo "Archive created	${ARCHIVE}"

create_tar_archive:
	@echo ${TAR_CREATE} 
	(cd ..; ${TAR_CREATE})

FILE=rsa-tools/perl-scripts
add_one_file:
	@echo ${ARCHIVE_CMD} ${FILE} ${POST_CMD}
	${ARCHIVE_CMD} ${FILE}  ${POST_CMD}

tar_archive:
	${MAKE} create_tar_archive
	${MAKE} fill_archive ARCHIVE_CMD='${TAR}' POST_CMD=''
	(cd ..; gzip -f ${ARCHIVE}.tar)

zip_archive:
	${MAKE} fill_archive ARCHIVE_CMD='${ZIP}' POST_CMD='${ZIP_EXCLUDE}'


PUB_LOGIN=jvanheld
PUB_SERVER=rsat.scmbb.ulb.ac.be
PUB_DIR=/home/jvanheld/public_html/rsat_distrib/
PUB_FORMAT=zip
publish:
	rsync -ruptvl -e ssh ${ARCHIVE_PREFIX}.${PUB_FORMAT} ${PUB_LOGIN}@${PUB_SERVER}:${PUB_DIR}