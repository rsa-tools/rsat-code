############################################################
#
# $Id: install_EnsEMBL_genomes.mk,v 1.2 2007/02/22 11:27:43 rsat Exp $
#
# Time-stamp: <>
#
############################################################

include ${RSAT}/makefiles/util.mk

DATE = `date +%Y%m%d_%H%M%S`

MAKEFILE=${RSAT}/makefiles/install_EnsEMBL_genomes.mk

V=1

VERSION=41

## One genome
ORG=Pan_troglodytes_EnsEMBL
FROM=-2000

install_one_ensembl_genome:
	install-organism -v ${V} -org ${ORG} -task config -up_from ${FROM}
	install-organism -v ${V} -org ${ORG} -task allup,dyads,oligos,start_stop,upstream_freq	\
		-batch -rm
