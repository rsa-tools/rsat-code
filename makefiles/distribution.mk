################################################################
## Prepare a release version of rsa-tools
## distribution.mk makefile
## Usage: make -f distribution.mk

MAKEFILE=${RSAT}/makefiles/distribution.mk
MAKE = make -sk -f ${MAKEFILE}

## Archive file
DATE=`date +%Y%m%d`
ARCHIVE=rsa-tools/rsa-tools_${DATE}

## Archive with tar
TAR =tar ${TAR_EXCLUDE} -cpv ${ARCHIVE}.tar
TGZ =tar ${TAR_EXCLUDE} -cpvz ${ARCHIVE}.tgz

## Archive with zip
ZIP_EXCLUDE=-x CVS '*~' 
ZIP =zip  -ry ${ARCHIVE}.zip 

EXCLUDE=${ZIP_EXCLUDE}
ARCHIVE_CMD=${ZIP}

### tags
usage:
	@echo "usage: make [-OPT='options'] target"
	@echo "implemented targets"
	@perl -ne 'if (/^([a-z]\S+):/){ print "\t$$1\n";  }' ${MAKEFILE}

zip_distrib:
	(cd ..;								\
	${ARCHIVE_CMD} rsa-tools/perl-scripts ${EXCLUDE};		\
	${ARCHIVE_CMD} rsa-tools/doc/manuals/*.pdf ${EXCLUDE};		\
	${ARCHIVE_CMD} rsa-tools/makefiles ${EXCLUDE};			\
	${ARCHIVE_CMD} rsa-tools/config/RSA.config.default ${EXCLUDE};	\
	)
