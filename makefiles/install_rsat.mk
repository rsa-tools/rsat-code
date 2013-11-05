############################################################
#
# $Id: install_rsat.mk,v 1.76 2013/07/19 06:29:14 jvanheld Exp $
#
# Time-stamp: <2003-05-23 09:36:00 jvanheld>
#
############################################################


################################################################
## This makefile serves to download and compile some C programs developed by
## third parties, and which are required for the Web site (e.g. RNSC, MCL) or
## can optionnally be used in some work flows (e.g. peak-motifs,
## cluster-motifs).

include ${RSAT}/makefiles/util.mk
MAKEFILE=${RSAT}/makefiles/install_rsat.mk

#################################################################
# Programs used for downloading and sycnrhonizing
SSH=-e 'ssh -x'

################################################################
## Install the RSAT package
install_rsat:
	make -f ${RSAT}/makefiles/init_RSAT.mk init
	make -f ${RSAT}/makefiles/init_RSAT.mk compile_all
	${MAKE} config_rsat
	make -f ${RSAT}/makefiles/install_software.mk


install_required_apps:
	make -f makefiles/install_software.mk install_seqlogo
	make -f makefiles/install_software.mk install_gnuplot


not_working:
	make -f makefiles/install_software.mk install_ghostscript

################################################################
## To do: write a Perl script that interactively prompts for
## configuration parameters and adapts the RSAT_config.mk and
## RSAT_config.props files.
config_rsat:
	@echo "INTERACTIVE CONFIGURATION NOT IMPLEMENTED YET" 
	@echo "Please edit manually the following files"
	@echo "	${RSAT}/RSAT_config.props"
	@echo "	${RSAT}/RSAT_config.mk"

