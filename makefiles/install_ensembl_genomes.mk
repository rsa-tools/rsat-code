############################################################
#
# $Id: install_EnsemblGenomes.mk,v 1.3 2011/04/10 13:29:22 rsat Exp $
#
# Time-stamp: <>
#
############################################################

## THIS MAKEFILE IS INCOMPLETE. BESIDES, IT SHOULD BE UPDATED TO TAKE
## INTO ACCOUNT THE NEW SCRIPTS DEVELOPED BY J. DELERCE ET
## AL. (2013-2014).

include ${RSAT}/makefiles/util.mk

DATE = `date +%Y%m%d_%H%M%S`

MAKEFILE=${RSAT}/makefiles/install_ensembl_genomes.mk

V=1

VERSION=41



################################################################
## Return the list of organisms supported at ensembl.org
list_supported:
	supported-organisms-ensembl


## One genome
ORG=Pan_troglodytes_EnsEMBL
FROM=-2000

install_one_ensembl_genome:
	install-organism -v ${V} -org ${ORG} -task config -up_from ${FROM}
	install-organism -v ${V} -org ${ORG} -task allup,dyads,oligos,start_stop,upstream_freq	\
		-batch -rm
