############################################################
#
# $Id: install_EnsEMBL_genomes.mk,v 1.1 2006/10/16 14:17:06 rsat Exp $
#
# Time-stamp: <>
#
############################################################

include ${RSAT}/makefiles/util.mk

DATE = `date +%Y%m%d_%H%M%S`

MAKEFILE=${RSAT}/makefiles/install_EnsEMBL_genomes.mk

V=1

OPT=

INSTALL_CMD=install-organism -v ${V}		\
		-org ${ORG}			\
		-task ${INSTALL_TASK}		\
		${OPT}

VERSION=41

## One genome
ORGANISM=Danio_rerio
FROM=-2000

install_one_ensembl_genome:
	rsync -ruptlv rsat@rsat.scmbb.ulb.ac.be:rsa-tools/data/genomes/${ORG}_EnsEMBL_${VERSION} ./data/genomes/
	ln -s ${ORG}_EnsEMBL_${VERSION} ${ORG}_EnsEMBL
	${MAKE} my_command MY_COMMAND="${INSTALL_CMD}" ORG=${ORGANISM}_EnsEMBL INSTALL_TASK="config -up_from ${FROM}"
	${MAKE} my_command MY_COMMAND="${INSTALL_CMD}" ORG=${ORGANISM}_EnsEMBL	\
		INSTALL_TASK=allup,dyads,oligos,start_stop,upstream_freq OPT="-batch -rm"
		>& install_${ORGANISM}_${VERSION}.txt \&
