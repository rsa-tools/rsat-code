############################################################
#
# $Id: install_rsat.mk,v 1.73 2013/07/19 02:43:28 jvanheld Exp $
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
## Packages specifically required for Ubuntu systems
install_ubuntu_prereq:
	sudo apt-get install libgd2-xpm-dev
	sudo apt-get install libxml2-dev
	sudo apt-get install libmysqlclient15-dev
	sudo apt-get install libdb-dev
	sudo apt-get install libberkeleydb-perl
