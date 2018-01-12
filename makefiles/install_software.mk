############################################################
#
# $Id: install_software.mk,v 1.54 2013/11/03 19:35:38 jvanheld Exp $
#
# Time-stamp: <2017-03-17 10:41:52 eda>
#
############################################################


################################################################
## This makefile manages the installation of bioinformatics software
## and the configuration of the paths and environment variables for
## the users.
include ${RSAT}/makefiles/util.mk
MAKEFILE=${RSAT}/makefiles/install_software.mk
MAKE=make -f ${MAKEFILE}
V=1

#################################################################
# Programs used for downloading and sycnrhonizing
WGET = wget -np -rNL 
#MAKE=nice -n 19 make -s -f ${MAKEFILE}
RSYNC_OPT = -ruptvl ${OPT}
SSH=-e 'ssh -x'
RSYNC = rsync ${RSYNC_OPT} ${SSH}
PYTHON=python2.7

################################################################
## Install the software tools.
INSTALL_TASKS=`${MAKE} | grep 'install_'`
list_installation_targets:
	@echo "Supported installation tasks"
	@echo ${INSTALL_TASKS} | perl -pe 's|\s|\n|g'

all_installations:
	@for task in ${INSTALL_TASKS}; do \
		echo "Installation task	$${task}" ; \
		${MAKE} $${task} ; \
	done

VERSIONS=`grep '_VERSION=' ${MAKEFILE} | grep -v '^\#' | perl -pe 's|_VERSION=|\t|'`
list_versions:
	@echo
	@echo "Software versions"
	@echo "${VERSIONS}"

################################################################
## Install the applications developed by third-parties and which are required
## or useful for RSAT.
EXT_APP_TARGETS=\
	install_vmatch \
	install_seqlogo \
	install_d3 \
	install_bedtools \
	install_ensembl_api \
	install_ghostscript \
	install_gnuplot \
	install_weblogo3 \
	install_mcl \
	install_rnsc \
	install_ensembl_bioperl \
	install_blast
list_ext_apps:
	@echo
	@echo "External applications to install"
	@echo "	${EXT_APP_TARGETS}" | perl -pe 's| |\n\t|g'

install_ext_apps:
	@${MAKE} ${EXT_APP_TARGETS}


not_working:
	make -f makefiles/install_software.mk install_ghostscript_macosx

EXT_APP_TARGETS_OPTIONAL=install_gibbs \
	install_consensus \
	install_patser \
	install_meme \
	install_bedtools

install_ext_apps_optional:
	@${MAKE} ${EXT_APP_TARGETS_OPTIONAL}

################################################################
## Download the vmatch program
##
## IMPORTANT: this program is an indispensable companion for the motif
## discovery tools in RSAT. If is used to purge sequences from
## redundant fragments. The program requires a license, which can be
## obtained (free of charge for academics) at http://www.vmatch.de/
VMATCH_VERSION=2.3.0
VMATCH_VERSION_MACOSX=vmatch-${VMATCH_VERSION}-Darwin_i386-64bit
VMATCH_VERSION_LINUX=vmatch-${VMATCH_VERSION}-Linux_x86_64-64bit
VMATCH_ARCHIVE=${VMATCH_VERSION}.tar.gz
VMATCH_BASE_DIR=${SRC_DIR}/vmatch
#VMATCH_URL=ftp://lscsa.de/pub/lscsa/
VMATCH_URL=http://www.vmatch.de/distributions/
install_vmatch_param:
	@echo "Parameters for vmatch instlalation"
	@echo "	VMATCH_VERSION		${VMATCH_VERSION}"
	@echo "	OS			${OS}"
	@echo "	VMATCH_VERSION_MACOSX	${VMATCH_VERSION_MACOSX}"
	@echo "	VMATCH_VERSION_LINUX	${VMATCH_VERSION_LINUX}"
	@echo "	VMATCH_ARCHIVE		${VMATCH_ARCHIVE}"
	@echo "	VMATCH_BASE_DIR		${VMATCH_BASE_DIR}"
	@echo "	VMATCH_URL		${VMATCH_URL}"

install_vmatch:
	@echo
	@echo "Installing vmatch for operating system ${OS}"
	${MAKE} _install_vmatch_${OS}
#	${MAKE} _vmatch_warning

_install_vmatch_macosx:
	${MAKE} VMATCH_VERSION=${VMATCH_VERSION_MACOSX} _download_vmatch _install_vmatch 

_install_vmatch_linux:
	${MAKE} VMATCH_VERSION=${VMATCH_VERSION_LINUX} _download_vmatch _install_vmatch


#VMATCH_SOURCE_DIR=vmatch_latest
_download_vmatch: 
	@echo ""
	@echo "Downloading vmatch in folder"
	@echo "	${VMATCH_BASE_DIR}"
	@mkdir -p ${VMATCH_BASE_DIR}
	wget --no-clobber --no-directories --no-verbose  --directory-prefix ${VMATCH_BASE_DIR} ${VMATCH_URL}/${VMATCH_ARCHIVE}
	@ls ${VMATCH_BASE_DIR}/${VMATCH_ARCHIVE}

VMATCH_SOURCE_DIR=${VMATCH_BASE_DIR}/${VMATCH_VERSION}
_install_vmatch:
	@echo ""
	@echo "VMATCH_SOURCE_DIR	${VMATCH_SOURCE_DIR}"
	@echo "Uncompressing vmatch tar archive"
	@echo "	${VMATCH_BASE_DIR}/${VMATCH_ARCHIVE}"
	@tar -xzf  ${VMATCH_BASE_DIR}/${VMATCH_ARCHIVE} -C ${VMATCH_BASE_DIR}
	@echo "Synchronizing vmatch and mkvtree in RSAT_BIN	${RSAT_BIN}"
	@echo "	${RSAT_BIN}"
	@${SUDO} rsync -ruptl ${VMATCH_SOURCE_DIR}/vmatch ${RSAT_BIN}/
	@${SUDO} rsync -ruptl ${VMATCH_SOURCE_DIR}/mkvtree ${RSAT_BIN}/

# _vmatch_warning:
# 	@echo ""
# 	@echo ""
# 	@echo "vmatch has been installed in bin folder ${RSAT_BIN}"
# 	@echo "IN ORDER TO GET A FUNCTIONAL COPY, YOU NEED TO REQUEST A LICENSE"
# 	@echo "	http://www.vmatch.de/"
# 	@echo "AND PLACE THE FILE vmatch.lic IN FOLDER"
# 	@echo "	 ${RSAT_BIN}"
# 	@echo ""

################################################################
## LOGOL
LOGOL_ARCHIVE=logol_1.7.1.orig.tar.bz2
LOGOL_URL=https://gforge.inria.fr/frs/download.php/file/33588
LOGOL_DIR=${SRC_DIR}/logol
install_logol: _download_logol _compile_logol

_download_logol:
	@mkdir -p ${LOGOL_DIR}
	@echo
	@echo "Downloading logol	${LOGOL_URL}"
	(cd ${LOGOL_DIR}; \
		wget --timestamping --no-directories  ${LOGOL_URL}/${LOGOL_ARCHIVE}; \
		tar -xpzf ${LOGOL_ARCHIVE})
	@echo "logol dir	${LOGOL_DIR}"

_compile_logol:
	@echo "Installing logol in RSAT_BIN	${RSAT_BIN}"
	@${SUDO} rsync -ruptl ${LOGOL_DIR}/weblogo/logol ${RSAT_BIN}/
	@${SUDO} rsync -ruptl ${LOGOL_DIR}/weblogo/template.* ${RSAT_BIN}/
	@${SUDO} rsync -ruptl ${LOGOL_DIR}/weblogo/logo.pm ${RSAT_BIN}/

################################################################
## Get and install the program seqlogo
SEQLOGO_URL=http://weblogo.berkeley.edu/release
SEQLOGO_TAR=weblogo.2.8.2.tar.gz
SEQLOGO_DIR=${SRC_DIR}/seqlogo
install_seqlogo: _download_seqlogo _compile_seqlogo

_download_seqlogo:
	@mkdir -p ${SEQLOGO_DIR}
	@echo
	@echo "Downloading seqlogo	${SEQLOGO_URL}"
	(cd ${SEQLOGO_DIR}; wget --timestamping --no-verbose --no-directories  ${SEQLOGO_URL}/${SEQLOGO_TAR}; tar -xpzf ${SEQLOGO_TAR})
	@echo "seqlogo dir	${SEQLOGO_DIR}"

_compile_seqlogo:
	@echo "Installing seqlogo in RSAT_BIN	${RSAT_BIN}"
	@${SUDO} rsync -ruptl ${SEQLOGO_DIR}/weblogo/seqlogo ${RSAT_BIN}/
	@${SUDO} rsync -ruptl ${SEQLOGO_DIR}/weblogo/template.* ${RSAT_BIN}/
	@${SUDO} rsync -ruptl ${SEQLOGO_DIR}/weblogo/logo.pm ${RSAT_BIN}/

################################################################
## Get and install the program weblogo.  Weblogo is an upgrade from
## seqlogo. Seqlogo required sequences in input, whereas weblogo takes
## either sequences or matrices.
##
## Source: http://weblogo.threeplusone.com/

# ## Manual installation of Weblogo in ${RSAT_BIN}
# WEBLOGO3_VERSION=3.3
# WEBLOGO3_TAR=weblogo-${WEBLOGO3_VERSION}.tar.gz
# WEBLOGO3_URL=https://weblogo.googlecode.com/files/${WEBLOGO3_TAR}
# WEBLOGO3_DIR=${SRC_DIR}/weblogo3
# install_weblogo3: _download_weblogo3 _compile_weblogo3

# _download_weblogo3:
# 	@mkdir -p ${WEBLOGO3_DIR}
# 	@echo "Getting weblogo3 using wget"
# 	(cd ${WEBLOGO3_DIR}; wget --no-clobber -nv -nd ${WEBLOGO3_URL}; tar -xpzf ${WEBLOGO3_TAR})
# 	@echo "weblogo3 dir	${WEBLOGO3_DIR}"

# _compile_weblogo3:
# 	@echo "Building weblogo vesion ${WEBLOGO3_VERSION} and installing in ${RSAT_BIN}"
# 	(cd ${WEBLOGO3_DIR}/weblogo-${WEBLOGO3_VERSION}; \
# 	${SUDO} ${PYTHON} setup.py install --prefix ${RSAT})
# 	@echo "weblogo installed in ${RSAT_BIN}"



################################################################
## Install weblogo3
install_weblogo3: weblogo3_from_githhub


################################################################
## Clone weblogo3
WEBLOGO_SRC_DIR=${SRC_DIR}/weblogo
weblogo3_from_githhub:
	@echo "Installing weblogo3 from github"
	@if [ -d ${WEBLOGO_SRC_DIR} ] ; then \
		echo "	Updating weblogo with git pull"; \
		cd ${WEBLOGO_SRC_DIR} ]; git pull;  \
	else \
		echo "	Cloning weblogo" ; \
		(cd ${SRC_DIR}; git clone https://github.com/WebLogo/weblogo.git); \
	fi
	@echo "	WEBLOGO_SRC_DIR	${WEBLOGO_SRC_DIR})"
	ln -fs ${WEBLOGO_SRC_DIR}/weblogo ${RSAT_BIN}/

## Installation via pip is simpler, but cannot be done on all RSAT
## servers because it requires admin rights.
PIP=pip
install_weblogo3_pip:
	@echo "Installing weblogo3 using pip"
	${SUDO} ${PIP} install --install-option "--install-scripts=${RSAT_BIN}" weblogo
#	${SUDO} ${PIP} install --target ${RSAT_BIN} --install-option "--install-scripts=${RSAT_BIN}" weblogo
#	${PIP} install --target ${RSAT_BIN} weblogo



################################################################
## Get and install the program gnuplot
#GNUPLOT_VERSION=4.6.4
GNUPLOT_VERSION=5.0.0
GNUPLOT_TAR=gnuplot-${GNUPLOT_VERSION}.tar.gz
GNUPLOT_URL=http://sourceforge.net/projects/gnuplot/files/gnuplot/${GNUPLOT_VERSION}/${GNUPLOT_TAR}
GNUPLOT_DIR=${SRC_DIR}/gnuplot
install_gnuplot:
	@echo
	@echo "Installing gnuplot for operating system ${OS}"
	${MAKE}  install_gnuplot_${OS}

install_gnuplot_macosx:
	brew install gnuplot

install_gnuplot_linux: _download_gnuplot_linux _compile_gnuplot_linux

install_gnuplot_linux: _download_gnuplot_linux _compile_gnuplot_linux

_download_gnuplot_linux:
	@mkdir -p ${GNUPLOT_DIR}
	@echo "Getting gnuplot using wget"
	(cd ${GNUPLOT_DIR}; wget -nv -nd ${GNUPLOT_URL}; tar -xpzf ${GNUPLOT_TAR})
	@echo "gnuplot dir	${GNUPLOT_DIR}"

_compile_gnuplot_linux:
	@echo "Compiling and installing gnuplot"
	(cd ${GNUPLOT_DIR}/gnuplot-${GNUPLOT_VERSION}; \
	./configure --prefix ${GNUPLOT_DIR}/gnuplot-${GNUPLOT_VERSION} --bindir ${RSAT_BIN}  && make; ${SUDO} make install)

################################################################
## Get and install the program ghostscript
## Note: for Mac users, please go to the ghostscript Web site
##
GS_URL=http://downloads.ghostscript.com/public/binaries/
GS_VERSION=9
## Beware: I use an older version 9.07 because linux version 9.10 issues
## warnings Unrecoverable error: stackunderflow in .setdistillerparams"
GS_SUBVER=07
GS_BIN=gs-${GS_VERSION}${GS_SUBVER}-linux_x86_64
GS_DISTRIB=ghostscript-${GS_VERSION}.${GS_SUBVER}-linux-x86_64
GS_TAR=${GS_DISTRIB}.tgz
GS_DIR=${SRC_DIR}/ghostscript
install_ghostscript:
	@echo
	@echo "Installing ghostscript (gs) for operating system ${OS}"
	${MAKE} install_ghostscript_${OS}

install_ghostscript_macosx:
	brew install ghostscript

install_ghostscript_linux: _download_gs _install_gs

_download_gs:
	@mkdir -p ${GS_DIR}
	@echo "Getting gs using wget"
	(cd ${GS_DIR}; wget -nv -nd ${GS_URL}/${GS_TAR}; tar -xpzf ${GS_TAR})
	@echo "gs dir	${GS_DIR}"

_install_gs:
	@echo "Installing gs in RSAT_BIN	${RSAT_BIN}"
	(cd ${GS_DIR}/${GS_DISTRIB}; ${SUDO} rsync -ruptvl ${GS_BIN} ${RSAT_BIN}/; cd ${RSAT_BIN}; ${SUDO} rm -f gs; ${SUDO} ln -s ${GS_BIN} gs)

#_compile_gs:
#	@echo "Compiling gs"
#	(cd ${GS_DIR}/${GS_DISTRIB}; ./configure && make)


################################################################
## Install the EnsEMBL Perl API
##
## Important note
## (http://bacteria.ensembl.org/info/docs/api/api_cvs.html): you must
## install version 1.2.3, not a more recent version. Starting with
## 1.2.4, major changes were made to the BioPerl API which have made
## it incompatible with Ensembl
BIOPERL_VERSION=1-2-3
BIOPERL_DIR=${RSAT}/ext_lib/bioperl-release-${BIOPERL_VERSION}

install_ensembl_api: install_ensembl_api_git

## Note: In some cases, there are delays between Ensembl and
## EnsemblGenome releases. To ensure compatibility, please check the
## versions of both distributions on http://ensemblgenomes.org/.
ENSEMBL_API_DIR=${RSAT}/ext_lib/ensemblgenomes-${ENSEMBLGENOMES_RELEASE}-${ENSEMBL_RELEASE}
install_ensembl_api_param:
	@echo "	BIOPERL_VERSION		${BIOPERL_VERSION}"
	@echo "	BIOPERL_DIR		${BIOPERL_DIR}"
	@echo "	ENSEMBL_RELEASE		${ENSEMBL_RELEASE}"
	@echo "	ENSEMBLGENOMES_RELEASE	${ENSEMBLGENOMES_RELEASE}"
	@echo "	ENSEMBL_API_DIR		${ENSEMBL_API_DIR}"

## Install the required modules for Ensembl API
install_ensembl_api_cvs:
	@(cd ${BIOPERL_DIR}; \
		echo "" ; \
		echo "Installing ensembl branch ${ENSEMBLGENOMES_RELEASE} version ${ENSEMBL_RELEASE}"; \
		echo "	ENSEMBL_API_DIR		${ENSEMBL_API_DIR}" ; \
		mkdir -p "${ENSEMBL_API_DIR}"; \
		cd ${ENSEMBL_API_DIR}; \
		echo  "Password is 'CVSUSER'" ; \
		cvs -d :pserver:cvsuser@cvs.sanger.ac.uk:/cvsroot/ensembl login ; \
		cvs -d :pserver:cvsuser@cvs.sanger.ac.uk:/cvsroot/ensembl checkout -r branch-ensemblgenomes-${ENSEMBLGENOMES_RELEASE}-${ENSEMBL_RELEASE} ensembl-api)
	@echo
	@${MAKE} install_ensembl_api_env

install_ensembl_api_git:
	@echo ""
	@echo "	ENSEMBL_API_DIR		${ENSEMBL_API_DIR}"
	@mkdir -p "${ENSEMBL_API_DIR}"
	@echo ""
	@echo "Cloning git for ensemblgenomes API branch ${ENSEMBLGENOMES_RELEASE}"
	@(cd ${ENSEMBL_API_DIR}; git clone https://github.com/EnsemblGenomes/ensemblgenomes-api.git ; \
		cd ensemblgenomes-api/ ; \
		git checkout release/eg/${ENSEMBLGENOMES_RELEASE} )
	@echo ""
	@echo "Getting git clone for ensembl API release ${ENSEMBL_RELEASE}"
	@(cd ${ENSEMBL_API_DIR}; \
		git clone https://github.com/Ensembl/ensembl-git-tools.git; \
		export PATH=${ENSEMBL_API_DIR}/ensembl-git-tools/bin:${PATH}; \
		git ensembl --clone api; \
		git ensembl --checkout --branch release/${ENSEMBL_RELEASE} api)

## Install git without using https (not supported on all servers)
install_ensembl_api_git_git:
	@(cd ${ENSEMBL_API_DIR}; \
		git clone git://github.com/Ensembl/ensembl.git; \
		git clone git://github.com/Ensembl/ensembl-compara.git; \
		git clone git://github.com/Ensembl/ensembl-funcgen.git; \
		git clone git://github.com/Ensembl/ensembl-tools.git; \
		git clone git://github.com/Ensembl/ensembl-variation.git; \
		git clone git://github.com/Ensembl/ensembl-git-tools.git; \
	)


################################################################
## Ensembl API requires Bioperl version 1-2-3, as quoted in their
## installation page:
## 	http://www.ensembl.org/info/docs/api/api_git.html
##
## "Important note: you must install version 1.2.3, not a more recent
## version. Starting with 1.2.4, major changes were made to the
## BioPerl API which have made it incompatible with Ensembl."
##
## ISSUES ERROR MESSAGE (2016-10-18)
## IS IT STILL REQUIRED ?
install_ensembl_bioperl:
	@echo ""
	@echo "Installing bioperl release ${BIOPERL_VERSION} (required for ensembl)"
	@echo "	BIOPERL_DIR		${BIOPERL_DIR}"
	@mkdir -p "${BIOPERL_DIR}"
	@if [ -d ${BIOPERL_DIR}/bioperl-live ] ; then \
		echo "Bioperl already installed"; \
	else \
		echo "Cloning bioperl" ; \
		(cd ${BIOPERL_DIR}; git clone https://github.com/bioperl/bioperl-live.git); \
	fi
	@(cd ${BIOPERL_DIR}/bioperl-live; git checkout bioperl-release-${BIOPERL_VERSION})
	@echo "bioperl-release-${BIOPERL_VERSION} installed in ${BIOPERL_DIR}"


## Print the settings for including ensembl in the PERL5LIB environment variable
install_ensembl_api_env:
	@echo
	@echo "ENSEMBL Perl modules are installed in directory ${ENSEMBL_API_DIR}"
	@echo
	@echo "BEWARE !"
	@echo "You need to paste the following lines in the bash profile ${RSAT}/RSAT_config.bashrc"
	@echo
	@echo '################################################################'
	@echo '## Default path for the Ensembl Perl modules and sofwtare tools'
	@echo 'export ENSEMBL_RELEASE=${ENSEMBL_RELEASE}'
	@echo 'export ENSEMBLGENOMES_RELEASE=${ENSEMBLGENOMES_RELEASE}'
	@echo 'export PATH=$${RSAT}/ext_lib/ensemblgenomes-$${ENSEMBLGENOMES_RELEASE}-$${ENSEMBL_RELEASE}/ensembl-git-tools/bin:$${PATH}'
	@echo 'export PERL5LIB=$${RSAT}/ext_lib/bioperl-release-$${BIOPERL_VERSION}/bioperl-live::$${PERL5LIB}'
	@echo 'export PERL5LIB=$${RSAT}/ext_lib/ensemblgenomes-$${ENSEMBLGENOMES_RELEASE}-$${ENSEMBL_RELEASE}/ensembl/modules::$${PERL5LIB}'
	@echo 'export PERL5LIB=$${RSAT}/ext_lib/ensemblgenomes-$${ENSEMBLGENOMES_RELEASE}-$${ENSEMBL_RELEASE}/ensembl-compara/modules::$${PERL5LIB}'
	@echo 'export PERL5LIB=$${RSAT}/ext_lib/ensemblgenomes-$${ENSEMBLGENOMES_RELEASE}-$${ENSEMBL_RELEASE}/ensembl-external/modules::$${PERL5LIB}'
	@echo 'export PERL5LIB=$${RSAT}/ext_lib/ensemblgenomes-$${ENSEMBLGENOMES_RELEASE}-$${ENSEMBL_RELEASE}/ensembl-functgenomics/modules::$${PERL5LIB}'
	@echo 'export PERL5LIB=$${RSAT}/ext_lib/ensemblgenomes-$${ENSEMBLGENOMES_RELEASE}-$${ENSEMBL_RELEASE}/ensembl-tools/modules::$${PERL5LIB}'
	@echo 'export PERL5LIB=$${RSAT}/ext_lib/ensemblgenomes-$${ENSEMBLGENOMES_RELEASE}-$${ENSEMBL_RELEASE}/ensembl-variation/modules::$${PERL5LIB}'

################################################################
## Install biomart Perl libraries
# from: http://www.ensembl.org/info/data/biomart/biomart_perl_api.html#downloadbiomartperlapi
# cvs -d :pserver:cvsuser@cvs.sanger.ac.uk:/cvsroot/biomart login                                                                                
# # passwd: CVSUSER                                                                                                                             
# cvs -d :pserver:cvsuser@cvs.sanger.ac.uk:/cvsroot/biomart co -r release-0_7 biomart-perl                                                    

# # Paste the text obtained on
# # [http://plants.ensembl.org/biomart/martservice?type=registry] into
# # the biomart-perl/conf/martURLLocation.xml file
# # PERL5LIB=${PERL5LIB}:${RSAT}/ext_lib/biomart-perl/lib                                                                                         # export PERL5LIB                                                                                                                                
# # check missing dependencies
# perl bin/configure.pl -r conf/registryURLPointer.xml                                                                                        
# #install the required modules


################################################################
## Install the graph-based clustering algorithm MCL
MCL_BASE_DIR=${SRC_DIR}/mcl
#MCL_VERSION=12-135
MCL_VERSION=14-137
MCL_ARCHIVE=mcl-${MCL_VERSION}.tar.gz
MCL_URL=http://www.micans.org/mcl/src/${MCL_ARCHIVE}
MCL_DISTRIB_DIR=${MCL_BASE_DIR}/mcl-${MCL_VERSION}
install_mcl: _download_mcl _compile_mcl

_download_mcl:
	@echo
	@echo "Downloading MCL"
	@mkdir -p ${MCL_BASE_DIR}
	wget -nd  --directory-prefix ${MCL_BASE_DIR} -rNL ${MCL_URL}
	(cd ${MCL_BASE_DIR}; tar -xpzf ${MCL_ARCHIVE})
	@echo ${MCL_DISTRIB_DIR}

MCL_COMPILE_DIR=`dirname ${RSAT_BIN}`
MCL_BIN_DIR=${RSAT_BIN}
_compile_mcl:
	@echo
	@echo "Installing MCL in dir ${MCL_BIN_DIR}"
	@mkdir -p ${MCL_COMPILE_DIR}
	(cd ${MCL_DISTRIB_DIR}; ./configure --prefix=${MCL_COMPILE_DIR} ; \
	make clean; make ; ${SUDO} make install)
	@echo "Please check that MCL binary directory is in your PATH"
	@echo "	${MCL_BIN_DIR}"

################################################################
## Install the graph-based clustering algorithm RNSC
RNSC_BASE_DIR=${SRC_DIR}/rnsc
RNSC_ARCHIVE=rnsc.tar.gz
RNSC_URL=http://www.cs.utoronto.ca/~juris/data/rnsc/${RNSC_ARCHIVE}
install_rnsc: _download_rnsc _compile_rnsc

_download_rnsc:
	@echo
	@echo "Downloading RNSC"
	@mkdir -p ${RNSC_BASE_DIR}
	wget --no-directories  --directory-prefix ${RNSC_BASE_DIR} -rNL ${RNSC_URL}
	(cd ${RNSC_BASE_DIR}; tar -xpf ${RNSC_ARCHIVE})
	@echo ${RNSC_BASE_DIR}

_compile_rnsc:
	@echo
	@echo "Installing RNSC in RSAT_BIN	${RSAT_BIN}"
	@${SUDO} mkdir -p ${RSAT_BIN}
	(cd ${RNSC_BASE_DIR}; make ;  \
	${SUDO} rsync -ruptvl rnsc ${RSAT_BIN}; \
	${SUDO} rsync -ruptvl rnscfilter ${RSAT_BIN}; \
	)
#	${SUDO} rsync -ruptvl rnscconvert ${RSAT_BIN}/; \

################################################################
## Download and install NCBI BLAST "legacy" version, which was written in C
##
install_blast: list_blast_param _download_blast _install_blast

_download_blast: _download_blast_${OS}

list_blast_param:
	@echo "Downloading blast"
	@echo "	Operating system	${OS}"
	@echo "	BLAST_BASE_DIR		${BLAST_BASE_DIR}"
	@echo "	BLAST_LINUX_ARCHIVE	${BLAST_LINUX_ARCHIVE}"
	@echo "	BLAST_URL		${BLAST_URL}"
	@echo "	BLAST_SOURCE_DIR	${BLAST_SOURCE_DIR}"
	@echo "	BLAST_BASE_DIR		${BLAST_BASE_DIR}"
	@echo "	RSAT_BIN		${RSAT_BIN}"

## Download BLAST for linux
BLAST_BASE_DIR=${SRC_DIR}/blast
BLAST_VERSION=2.2.26
BLAST_LINUX_ARCHIVE=blast-${BLAST_VERSION}-${ARCHITECTURE}-linux.tar.gz
#BLAST_URL=ftp://ftp.ncbi.nih.gov/blast/executables/release/LATEST/
#BLAST_URL=ftp://ftp.ncbi.nih.gov/blast/executables/release/legacy/${BLAST_VERSION}
BLAST_URL=ftp://ftp.ncbi.nlm.nih.gov/blast/executables/legacy/${BLAST_VERSION}
BLAST_SOURCE_DIR=blast_latest
_download_blast_linux:
	@mkdir -p ${BLAST_BASE_DIR}
	wget --no-directories  --directory-prefix ${BLAST_BASE_DIR} -rNL ${BLAST_URL}/${BLAST_LINUX_ARCHIVE}
	(cd ${BLAST_BASE_DIR}; tar -xzf ${BLAST_LINUX_ARCHIVE}; )
	@echo ${BLAST_BASE_DIR}

## Download BLAST for Mac OS X
BLAST_BASE_DIR=${SRC_DIR}/blast
BLAST_MAC_ARCHIVE=blast-${BLAST_VERSION}-universal-macosx.tar.gz
BLAST_SOURCE_DIR=blast-${BLAST_VERSION}
_download_blast_macosx:
	@mkdir -p ${BLAST_BASE_DIR}
	wget --no-directories  --directory-prefix ${BLAST_BASE_DIR} -rNL ${BLAST_URL}/${BLAST_MAC_ARCHIVE}
	(cd ${BLAST_BASE_DIR}; tar -xzf ${BLAST_MAC_ARCHIVE})
	@echo ${BLAST_BASE_DIR}

## Install BLAST executablesin RSAT_BIN directory
_install_blast:
	@${SUDO} mkdir -p ${RSAT_BIN}
	${SUDO} rsync -ruptvl ${BLAST_BASE_DIR}/${BLAST_SOURCE_DIR}/bin/blastall ${RSAT_BIN}
	${SUDO} rsync -ruptvl ${BLAST_BASE_DIR}/${BLAST_SOURCE_DIR}/bin/formatdb ${RSAT_BIN}
	@echo "Please edit the RSAT configuration file"
	@echo "	${RSAT}/RSAT_config.props"
	@echo "and copy-paste the following line to specify the BLAST bin pathway"
	@echo "	blast_dir=${RSAT_BIN}"
	@echo "This will allow RSAT programs to idenfity BLAST path on this server."
	@echo
	@echo "You can also add the BLAST bin directory in your path."
	@echo "If your shell is bash"
	@echo "	export PATH=$${RSAT_BIN}:$${PATH}"
	@echo "If your shell is csh or tcsh"
	@echo "	setenv PATH $${RSAT_BIN}:$${PATH}"

################################################################
## Install blast+, the C++ version of BLAST
##
## Note that commands changed, formatdb and blastall do not exist
## anymore. The older versions (written in C) are found in the legacy
## directory of the FTP site.
install_blast+: list_blast+_param _download_blast+ _install_blast+

list_blast+_param:
	@echo "Downloading blast"
	@echo "	Operating system	${OS}"
	@echo "	BLASTPLUS_URL		${BLASTPLUS_URL}"
	@echo "	BLASTPLUS_ARCHIVE	${BLASTPLUS_ARCHIVE}"
	@echo "	BLASTPLUS_SOURCE_DIR	${BLASTPLUS_SOURCE_DIR}"
	@echo "	BLASTPLUS_BASE_DIR	${BLASTPLUS_BASE_DIR}"
	@echo "	BLASTPLUS_ARCHIVE	${BLASTPLUS_ARCHIVE}"
	@echo "	RSAT_BIN		${RSAT_BIN}"

## Download BLASTPLUS archive from NCBI ftp site
BLASTPLUS_BASE_DIR=${SRC_DIR}/blast
BLASTPLUS_VERSION=2.5.0
BLASTPLUS_SOURCE_DIR=ncbi-blast-${BLASTPLUS_VERSION}+
#BLASTPLUS_LINUX_ARCHIVE=blast-*-${ARCHITECTURE}-linux.tar.gz
BLASTPLUS_ARCHIVE=ncbi-blast-${BLASTPLUS_VERSION}+-${ARCHITECTURE}-${OS}.tar.gz
BLASTPLUS_URL=ftp://ftp.ncbi.nih.gov/blast/executables/blast+/${BLASTPLUS_VERSION}/
BLASTPLUS_BASE_DIR=${SRC_DIR}/blast
_download_blast+:
	@echo "Downloading blast+ from ${BLASTPLUS_URL}/${BLASTPLUS_ARCHIVE}"
	@mkdir -p ${BLASTPLUS_BASE_DIR}
	wget --no-clobber --no-directories  --directory-prefix ${BLASTPLUS_BASE_DIR} ${BLASTPLUS_URL}/${BLASTPLUS_ARCHIVE}
	(cd ${BLASTPLUS_BASE_DIR}; tar -xzf ${BLASTPLUS_ARCHIVE})
	@echo ${BLASTPLUS_BASE_DIR}

## Install BLASTPLUS from downloaded BLASTPLUS archive
_install_blast+:
	@echo "Installing blast+"
	@${SUDO} mkdir -p ${RSAT_BIN}
	${SUDO} rsync -ruptvl ${BLASTPLUS_BASE_DIR}/${BLASTPLUS_SOURCE_DIR}/bin/* ${RSAT_BIN}
	@echo "Please edit the RSAT configuration file"
	@echo "	${RSAT}/RSAT_config.props"
	@echo "and copy-paste the following line to specify the BLASTPLUS bin pathway"
	@echo "	blast_dir=${RSAT_BIN}"
	@echo "This will allow RSAT programs to idenfity BLASTPLUS path on this server."
	@echo
	@echo "You can also add the BLASTPLUS bin directory in your path."
	@echo "If your shell is bash"
	@echo "	export PATH=$${RSAT_BIN}:$${PATH}"
	@echo "If your shell is csh or tcsh"
	@echo "	setenv PATH $${RSAT_BIN}:$${PATH}"


################################################################
## Get the d3 javascript library, for motif clustering and dynamical display
install_d3: _download_d3

#D3_URL=http://d3js.org/d3.v3.zip
D3_URL=https://github.com/mbostock/d3/releases/download/v3.4.1/d3.v3.zip
D3_DIR=${SRC_DIR}/d3.${D3_VERSION}
D3_VERSION=v3
D3_ARCHIVE=d3.${D3_VERSION}.zip
_download_d3:
	@mkdir -p ${D3_DIR}
	(cd ${D3_DIR};  \
	wget http://mbostock.github.com/d3/d3.js ;  \
	wget http://mbostock.github.com/d3/d3.min.js )

#	wget http://mbostock.github.com/d3/d3.layout.js)


################################################################
## Install fastqc, a software tool to control the quality of read
## files (next generation sequencing).
## http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip
FASTQC_VER=0.11.5
FASTQC_ZIP=fastqc_v${FASTQC_VER}.zip
FASTQC_URL=http://www.bioinformatics.babraham.ac.uk/projects/fastqc/${FASTQC_ZIP}
FASTQC_DOWNLOAD_DIR=${SRC_DIR}/fastqc
FASTQC_INSTALL_DIR=${FASTQC_DOWNLOAD_DIR}
FASTQC_EXEC_DIR=${FASTQC_INSTALL_DIR}/FastQC
FASTQC=${FASTQC_EXEC_DIR}/fastqc
install_fastqc: _download_fastqc _install_fastqc

_download_fastqc:
	@mkdir -p ${FASTQC_DOWNLOAD_DIR}
	@echo "Getting fastqc using wget"
	(cd ${FASTQC_DOWNLOAD_DIR}; wget  --no-clobber -nv -nd ${FASTQC_URL}; unzip ${FASTQC_ZIP})
	@echo "	fastqc download dir	${FASTQC_DOWNLOAD_DIR}"
	@echo "	fastqc install dir	${FASTQC_INSTALL_DIR}"
	@chmod 755 ${FASTQC}
	@echo "	fastqc executable dir	${FASTQC_EXEC_DIR}"
	@echo "	fastqc script    	${FASTQC}"

_install_fastqc:
	@echo
	@ln -s -f ${FASTQC_INSTALL_DIR}/FastQC/fastqc ${RSAT_BIN}/fastqc
	@echo "	fastqc link    	 ${RSAT_BIN}/fastqc"
#	@echo "YOU NEED TO ADD THE FASTQC EXEC DIRECTORY TO YOUR PATH"
#	@echo "export PATH=$${PATH}:${FASTQC_EXEC_DIR}"

################################################################
## Install BEDTools
##
## BEDTools is a collection of utilities for comparing, summarizing,
## and intersecting genomic features in BED, GTF/GFF, VCF and BAM
## formats.
#bedtools: git_bedtools compile_bedtools _compile_bedtools
install_bedtools: _download_bedtools _compile_bedtools _install_bedtools

#http://bedtools.googlecode.com/files/BEDTools.v2.17.0.tar.gz
#BED_VERSION=2.13.3
#BED_VERSION=2.17.0
BED_VERSION=2.25.0
BED_RELEASE=2
BED_ARCHIVE=bedtools-${BED_VERSION}.tar.gz
BED_URL=https://github.com/arq5x/bedtools2/releases/download/v${BED_VERSION}/${BED_ARCHIVE}
BED_BASE_DIR=${SRC_DIR}/bedtools
#BED_URL=http://bedtools.googlecode.com/files/${BED_ARCHIVE} #deprecated
#BEDTOOL_MANUAL=http://bedtools.googlecode.com/files/BEDTools-User-Manual.v4.pdf
_download_bedtools:
	@echo
	@echo "Downloading bedtools ${BED_VERSION}"
	@echo
	@mkdir -p ${BED_BASE_DIR}
	(cd ${BED_BASE_DIR}; wget -nv -nd ${BED_URL} ; tar xpfz ${BED_ARCHIVE})

BED_GIT_DIR=${SRC_DIR}/bedtools
_git_bedtools:
	@mkdir -p ${BED_GIT_DIR}
	(cd ${SRC_DIR}; git clone git://github.com/arq5x/bedtools.git)

#BED_SRC_DIR=${BED_GIT_DIR}
BED_SRC_DIR=${BED_BASE_DIR}/bedtools${BED_RELEASE}
BED_BIN_DIR=${BED_SRC_DIR}/bin
_compile_bedtools:
	@echo
	@echo "Compiling bedtools from ${BED_SRC_DIR}"
	@echo
	@mkdir -p ${BED_SRC_DIR}
	(cd ${BED_SRC_DIR}; make all)

_install_bedtools:
	@echo
	@echo "Installing bedtools binaries from ${BEN_BIN_DIR} in RSAT_BIN	${RSAT_BIN}"
	@echo
	@mkdir -p ${RSAT_BIN}
	@${SUDO} rsync -ruptvl ${BED_BIN_DIR}/* ${RSAT_BIN}/

################################################################
## Install Biotoolbox
install_biotoolbox: _download_biotoolbox

BTB_VERSION=1.9.4
BTB_ARCHIVE=biotoolbox_v${BTB_VERSION}.tgz
BTB_URL=http://biotoolbox.googlecode.com/files/${BTB_ARCHIVE}
BTB_BASE_DIR=${SRC_DIR}/biotoolbox
BTB_DOWNLOAD_DIR=${BTB_BASE_DIR}/biotoolbox
_download_biotoolbox:
	@echo
	@echo "Downloading biotoolbox ${BTB_VERSION}"
	@echo
	@mkdir -p ${BTB_BASE_DIR}
	(cd ${BTB_BASE_DIR}; wget -nv -nd ${BTB_URL} ; tar -xpzf ${BTB_ARCHIVE})
	@echo ${BTB_DOWNLOAD_DIR}

################################################################
## Install MEME (Tim Bailey)
MEME_BASE_DIR=${SRC_DIR}/MEME
MEME_VERSION=4.10.2
#MEME_PATCH=_1 #Oct2015
#http://ebi.edu.au/ftp/software/MEME/4.9.1/meme_4.9.1_2.tar.gz
#http://meme-suite.org/meme-software/4.10.2/meme_4.10.2.tar.gz en Oct2015
MEME_ARCHIVE=meme_${MEME_VERSION}.tar.gz
MEME_URL=meme-suite.org/meme-software/${MEME_VERSION}/${MEME_ARCHIVE}
MEME_INSTALL_SUBDIR=${SRC_DIR}/MEME
MEME_INSTALL_DIR=${MEME_INSTALL_SUBDIR}/meme_${MEME_VERSION}
MEME_LOCAL_URL=http://localhost/meme
install_meme: _download_meme _compile_meme _after_meme

_download_meme:
	@echo
	@echo "Downloading MEME ${MEME_VERSION}"
	@echo
	@echo "MEME base directory	${MEME_BASE_DIR}"
	@mkdir -p ${MEME_BASE_DIR}
	(cd ${MEME_BASE_DIR}; wget -nv -nd ${MEME_URL})
	@echo "MEME base directory	${MEME_BASE_DIR}"
	@echo "MEME install directory	${MEME_INSTALL_DIR}"

## BEWARE, MEME creates a lot of folders and files, it should NOT be
## installed in /usr/local nor ${RSAT} directories
MEME_COMPILE_DIR=${MEME_INSTALL_DIR}
MEME_BIN_DIR=${MEME_COMPILE_DIR}/bin
_compile_meme:
	@echo
	@echo "Compiling MEME ${MEME_VERSION} in dir ${MEME_INSTALL_DIR}"
	@echo "	MEME_INSTALL_DIR	${MEME_INSTALL_DIR}"
	@echo "	MEME_INSTALL_SUBDIR	${MEME_INSTALL_SUBDIR}"
	@echo "	MEME_COMPILE_DIR	${MEME_COMPILE_DIR}"
	@mkdir -p ${MEME_INSTALL_DIR}
	(cd ${MEME_INSTALL_SUBDIR}; tar -xpzf ${MEME_BASE_DIR}/${MEME_ARCHIVE})
#	@echo "MEME configuration prefix	${MEME_CONFIG_PREFIX}"
	(cd ${MEME_INSTALL_DIR}; ./configure --prefix=${MEME_COMPILE_DIR} --with-url="${MEME_LOCAL_URL}")
	(cd ${MEME_INSTALL_DIR}; make clean; make ; make test; ${SUDO} make install)
	@echo "MEME installed in ${MEME_COMPILE_DIR}"

_after_meme:
	@echo "Creating links to meme"
	@echo "	MEME_BIN_DIR	${MEME_BIN_DIR}"
	@echo "	MEME_VERSION	${MEME_VERSION}"
	@echo "	RSAT_BIN	${RSAT_BIN}"
	cd ${RSAT_BIN}; rm -f meme; ln -s  ${MEME_BIN_DIR}/meme .
	@echo "Please edit the bashrc file"
	@echo "and copy-paste the following lines to specify the MEME bin pathway"
	@echo "	export PATH=${MEME_BIN_DIR}:$${PATH}"

################################################################
## Install a clustering algorithm "cluster"
CLUSTER_BASE_DIR=${SRC_DIR}/cluster
CLUSTER_VERSION=1.50
CLUSTER_ARCHIVE=cluster-${CLUSTER_VERSION}.tar.gz
CLUSTER_URL=http://bonsai.hgc.jp/~mdehoon/software/cluster/${CLUSTER_ARCHIVE}
CLUSTER_DISTRIB_DIR=${CLUSTER_BASE_DIR}/cluster-${CLUSTER_VERSION}
install_cluster: _download_cluster _compile_clustera

_download_cluster:
	@echo
	@echo "Downloading CLUSTER"
	@mkdir -p ${CLUSTER_BASE_DIR}
	wget -nd  --directory-prefix ${CLUSTER_BASE_DIR} -rNL ${CLUSTER_URL}
	(cd ${CLUSTER_BASE_DIR}; tar -xpzf ${CLUSTER_ARCHIVE})
	@echo ${CLUSTER_DISTRIB_DIR}

CLUSTER_COMPILE_DIR=`dirname ${RSAT_BIN}`
CLUSTER_BIN_DIR=${CLUSTER_COMPILE_DIR}/bin
_compile_cluster:
	@echo
	@echo "Installing CLUSTER in dir	${CLUSTER_COMPILE_DIR}"
	@mkdir -p ${CLUSTER_COMPILE_DIR}
	(cd ${CLUSTER_DISTRIB_DIR}; ./configure --without-x --prefix=${CLUSTER_COMPILE_DIR} ; \
	make clean; make ; ${SUDO} make install)

################################################################
## Generic call for installing a program. This tag is called with
## specific parameters for each program (consensus, patser, ...)
PROGRAM=consensus
PROGRAM_DIR=${SRC_DIR}/${PROGRAM}
PROGRAM_ARCHIVE=`ls -1t ${SRC_DIR}/${PROGRAM}* | head -1`
uncompress_program:
	@echo installing ${PROGRAM_ARCHIVE} in dir ${PROGRAM_DIR}
	@mkdir -p ${PROGRAM_DIR}
	@(cd ${PROGRAM_DIR} ;				\
	gunzip -c ${PROGRAM_ARCHIVE} | tar -xf - )

################################################################
## Common frame for installing programs
INSTALLED_PROGRAM=`ls -1t ${SRC_DIR}/${PROGRAM}/${PROGRAM}*`
_compile_program:
	(cd ${PROGRAM_DIR}; make ${COMPILE_OPT})
	(cd bin; ln -fs ${INSTALLED_PROGRAM} ./${PROGRAM})


################################################################
## Install Andrew Neuwald's gibbs sampler (1995 version)
GIBBS_DIR=${SRC_DIR}/gibbs/gibbs9_95
_compile_gibbs:
	${MAKE} uncompress_program PROGRAM=gibbs
	(cd ${GIBBS_DIR}; ./compile; cd ${GIBBS_DIR}/code; make clean)
	(cd bin; ln -fs ${GIBBS_DIR}/gibbs ./gibbs)

################################################################
## Get and install patser (matrix-based pattern matching)
PATSER_VERSION=patser-v3b.5
#PATSER_VERSION=patser-v3e.1
PATSER_TAR=${PATSER_VERSION}.tar.gz
PATSER_URL=http://stormo.wustl.edu/src/
#PATSER_URL=ftp://www.genetics.wustl.edu/pub/stormo/Consensus
PATSER_DIR=${SRC_DIR}/patser/${PATSER_VERSION}
PATSER_APP=`cd ${PATSER_DIR} ; ls -1tr patser-v* | grep -v .tar | tail -1 | xargs`
install_patser: _download_patser _compile_patser

_download_patser:
	@mkdir -p ${PATSER_DIR}
	@echo "Getting patser using wget"
	wget --no-directories  --directory-prefix ${PATSER_DIR} -rNL ${PATSER_URL}/${PATSER_TAR}
	(cd ${PATSER_DIR}; tar -xpzf ${PATSER_TAR})
#	(cd ${PATSER_DIR}; wget -nv  ${PATSER_URL}/${PATSER_TAR}; tar -xpzf ${PATSER_TAR})
	@echo "patser dir	${PATSER_DIR}"

_compile_patser:
	@echo "Compiling patser in RSAT_BIN	${PATSER_DIR}"
	(cd ${PATSER_DIR}; rm *.o; make)
	@echo "Installing patser in RSAT_BIN	${RSAT_BIN}"
	${SUDO} rsync -ruptvl ${PATSER_DIR}/${PATSER_APP} ${RSAT_BIN}
	(cd ${RSAT_BIN}; ${SUDO} ln -fs ${PATSER_APP} patser)
	@echo "ls -ltr ${RSAT_BIN}/patser*"


################################################################
## Install consensus (J.Hertz)
#CONSENSUS_VERSION=consensus-v6c.1 ## Not distributed anymore ?
CONSENSUS_VERSION=consensus-v6c
CONSENSUS_TAR=${CONSENSUS_VERSION}.tar.gz
#CONSENSUS_URL=ftp://www.genetics.wustl.edu/pub/stormo/Consensus
CONSENSUS_URL=http://stormo.wustl.edu/src/
#CONSENSUS_URL=http://gzhertz.home.comcast.net/~gzhertz/
CONSENSUS_DIR=${SRC_DIR}/consensus/${CONSENSUS_VERSION}
install_consensus: _download_consensus _compile_consensus

_download_consensus:
	@echo
	@echo "Downloading ${CONSENSUS_VERSION}"
	@echo
	@mkdir -p ${CONSENSUS_DIR}
	(cd ${SRC_DIR}/consensus; wget -nv -nd --no-clobber ${CONSENSUS_URL}/${CONSENSUS_TAR}; cd ${CONSENSUS_DIR}; tar -xpzf ../${CONSENSUS_TAR})
	@echo "consensus dir	${CONSENSUS_DIR}"

_compile_consensus:
	@echo "Compiling consensus in RSAT_BIN	${CONSENSUS_DIR}"
	(cd ${CONSENSUS_DIR}; rm *.o; make CPPFLAGS="")
	@echo "Installing consensus in RSAT_BIN	${RSAT_BIN}"
	${SUDO} rsync -ruptvl ${CONSENSUS_DIR}/${CONSENSUS_APP} ${RSAT_BIN}
	(cd ${RSAT_BIN}; ${SUDO} ln -fs ${CONSENSUS_APP} consensus)
	@echo "ls -ltr ${RSAT_BIN}/consensus*"
#	${MAKE} _compile_program PROGRAM=consensus COMPILE_OPT='CPPFLAGS=""'

################################################################
## UCSC tools (developed by Jim Kent)
install_ucsc_userapps: _download_ucsc_userapps _install_ucsc_userapps
UCSC_URL=http://hgdownload.cse.ucsc.edu/admin/exe/${UCSC_OS}/
UCSC_USERAPP_DIR=${SRC_DIR}/UCSC_userApps
_download_ucsc_userapps:
	@echo
	@echo "Downloading UCSC userApps from ${UCSC_URL}"
	(mkdir -p ${UCSC_USERAPP_DIR}; cd ${UCSC_USERAPP_DIR};  wget --no-verbose -rNL --no-directories --page-requisites --execute robots=off ${UCSC_URL})
	@echo "UCSC userApps downloaded to	${UCSC_USERAPP_DIR}"

_install_ucsc_userapps:
	(cd ${UCSC_USERAPP_DIR})


################################################################
################################################################
###########  SOFTWARE FOR HIGH-THROUGHPUT SEQUENCING ###########
################################################################
################################################################


################################################################
## TopHat - discovery splice junctions with RNA-seq
install_tophat: _download_tophat _compile_tophat

TOPHAT_BASE_DIR=${SRC_DIR}/TopHat
#TOPHAT_VERSION=1.2.0
#TOPHAT_ARCHIVE=tophat-${TOPHAT_VERSION}.tar.gz
#TOPHAT_URL=http://tophat.cbcb.umd.edu/downloads/${TOPHAT_ARCHIVE}
TOPHAT_VERSION=2.0.14
TOPHAT_ARCHIVE=tophat-${TOPHAT_VERSION}.tar.gz
TOPHAT_URL=https://ccb.jhu.edu/software/tophat/downloads/${TOPHAT_ARCHIVE}
TOPHAT_DISTRIB_DIR=${TOPHAT_BASE_DIR}/tophat-${TOPHAT_VERSION}
_download_tophat:
	@echo
	@echo "Downloading TopHat"
	@mkdir -p ${TOPHAT_BASE_DIR}
	@echo "TOPHAT_BASE_DIR	${TOPHAT_BASE_DIR}"
	wget -nd -N  --directory-prefix ${TOPHAT_BASE_DIR} -rNL ${TOPHAT_URL}
	(cd ${TOPHAT_BASE_DIR}; tar -xpzf ${TOPHAT_ARCHIVE})
	@echo ${TOPHAT_DISTRIB_DIR}

_clone_tophatgit:
	@echo
	@echo "cloning TopHat from git"
	(cd ${TOPHAT_BASE_DIR}; git clone https://github.com/infphilo/tophat.git)


## COMPILATION DOES NOT WORK - TO CHECK
TOPHAT_COMPILE_DIR=`dirname ${RSAT_BIN}`
TOPHAT_BIN_DIR=${TOPHAT_COMPILE_DIR}/bin
_compile_tophat:
	@echo
	@echo "Compiling and installing TopHat"
	@mkdir -p ${TOPHAT_COMPILE_DIR}
	(cd ${TOPHAT_DISTRIB_DIR}; ./configure --prefix=${TOPHAT_COMPILE_DIR}; \
	make ; ${SUDO} make install)


################################################################
## MACS, peak-calling program
MACS_BASE_DIR=${SRC_DIR}/MACS
MACS_VERSION=1.4.2
MACS_PATCH=-1
MACS_ARCHIVE=MACS-${MACS_VERSION}${MACS_PATCH}.tar.gz
MACS_URL=https://github.com/downloads/taoliu/MACS/${MACS_ARCHIVE}
MACS_DISTRIB_DIR=${MACS_BASE_DIR}/MACS-${MACS_VERSION}
install_macs: _download_macs _compile_macs

_download_macs:
	@echo
	@echo "Downloading MACS"
	@mkdir -p ${MACS_BASE_DIR}
	wget -nd  --directory-prefix ${MACS_BASE_DIR} -rNL ${MACS_URL}
	(cd ${MACS_BASE_DIR}; tar -xpzf ${MACS_ARCHIVE})
	@echo ${MACS_DISTRIB_DIR}

_compile_macs:
	(cd ${MACS_DISTRIB_DIR}; ${SUDO} python setup.py install )

## MACS version 2
install_macs2:
	${MAKE} _download_macs  MACS_VERSION=2.0.9 MACS_PATCH=-1
	${MAKE} _compile_macs  MACS_VERSION=2.0.9

## The python Cython library is required for installing macs2
install_macs2_dependencies:
	sudo easy_install cython

################################################################
## PeakSplitter, program for splitting the sometimes too large regions
## returned by MACS into "topological" peaks.
PEAKSPLITTER_BASE_DIR=${SRC_DIR}/PeakSplitter
PEAKSPLITTER_ARCHIVE=PeakSplitter_Cpp_1.0.tar.gz
PEAKSPLITTER_URL=http://www.ebi.ac.uk/sites/ebi.ac.uk/files/groups/bertone/software/${PEAKSPLITTER_ARCHIVE}
PEAKSPLITTER_DISTRIB_DIR=${PEAKSPLITTER_BASE_DIR}/PeakSplitter_Cpp
install_peaksplitter: _download_peaksplitter _compile_peaksplitter 

_download_peaksplitter:
	@echo
	@echo "Downloading PeakSplitter"
	@mkdir -p ${PEAKSPLITTER_BASE_DIR}
	wget -nd  --directory-prefix ${PEAKSPLITTER_BASE_DIR} -rNL ${PEAKSPLITTER_URL}
	(cd ${PEAKSPLITTER_BASE_DIR}; tar -xpzf ${PEAKSPLITTER_ARCHIVE})
	@echo ${PEAKSPLITTER_DISTRIB_DIR}

_compile_peaksplitter: _compile_peaksplitter_${OS}

_compile_peaksplitter_macosx:
	${MAKE} __compile_peaksplitter OS=MacOS

_compile_peaksplitter_linux:
	${MAKE} __compile_peaksplitter OS=Linux64

__compile_peaksplitter:
	(cd ${PEAKSPLITTER_DISTRIB_DIR}; ${SUDO} rsync -ruptvl -e ssh PeakSplitter_${OS}/PeakSplitter ${RSAT_BIN}/)

################################################################
## sickle
SICKLE_BASE_DIR=${SRC_DIR}/sickle
SICKLE_URL=https://github.com/ucdavis-bioinformatics/sickle.git
SICKLE_DISTRIB_DIR=${SICKLE_BASE_DIR}/SICKLE_V${SICKLE_VERSION}
install_sickle: _clone_sickle _compile_sickle

_clone_sickle:
	@if [ -d ${SICKLE_BASE_DIR} ] ; \
		then echo "Updating sickle"; \
		(cd ${SICKLE_BASE_DIR}; git pull) \
	else \
		echo "Cloning sickle" ; \
		(cd ${SRC_DIR}; git clone ${SICKLE_URL}) \
	fi

_compile_sickle:
		(cd ${SICKLE_BASE_DIR}; make clean; make; ${SUDO} rsync -ruptvl sickle ${RSAT_BIN})

################################################################
## FindPeaks
FINDPEAKS_VERSION=3-1-9-2
FINDPEAKS_VERSION_DIR=3.1.9.2
FINDPEAKS_URL=http://www.bcgsc.ca/platform/bioinfo/software/findpeaks/releases/${FINDPEAKS_VERSION_DIR}/findpeaks${FINDPEAKS_VERSION}-tar.gz
install_findpeaks:
	@echo "TO BE IMPLEMENTED"

################################################################
## SICER (peak calling for large regions e.g. methylation)
SICER_BASE_DIR=${SRC_DIR}/sicer
SICER_VERSION=1.1
SICER_ARCHIVE=SICER_V${SICER_VERSION}.tgz
SICER_URL=http://home.gwu.edu/~wpeng/${SICER_ARCHIVE}
SICER_DISTRIB_DIR=${SICER_BASE_DIR}/SICER_V${SICER_VERSION}
install_sicer: _download_sicer _compile_sicer numpy_and_scipy

_download_sicer:
	@echo
	@echo "Downloading SICER"
	@mkdir -p ${SICER_BASE_DIR}
	wget -nd  --directory-prefix ${SICER_BASE_DIR} -rNL ${SICER_URL}

## This installation is VERY tricky. The user has to replace the
## hard-coded path in 3 shell files
_compile_sicer:
	@echo
	@echo "Installing SICER in dir	${SICER_DISTRIB_DIR}"
	(cd ${SICER_BASE_DIR}; tar -xpzf ${SICER_ARCHIVE})
	@echo ${SICER_DISTRIB_DIR}
	@for f in SICER.sh SICER-rb.sh SICER-df.sh SICER-df-rb.sh; do \
		${MAKE} _sicer_path_one_file SICER_FILE=$${f} ; \
	done
	@echo "SICER requires python library Numpy	http://sourceforge.net/projects/numpy/files/NumPy/"
	@echo "SICER requires python library scipy	http://www.scipy.org/"
	@${MAKE} numpy_and_scipy

SICER_FILE=SICER.sh
_sicer_path_one_file:
	@if [ -f "${SICER_DISTRIB_DIR}/SICER/${SICER_FILE}.ori" ] ; then \
		echo "	backup copy already exists ${SICER_DISTRIB_DIR}/SICER//${SICER_FILE}.ori" ; \
	else \
		echo "	creating backup copy ${SICER_DISTRIB_DIR}/SICER/${SICER_FILE}.ori" ; \
		rsync -ruptl ${SICER_DISTRIB_DIR}/SICER/${SICER_FILE} ${SICER_DISTRIB_DIR}/SICER/${SICER_FILE}.ori ; \
	fi
	@perl -pe 's|/home/data/SICER${SICER_VERSION}|${SICER_DISTRIB_DIR}|g' -i ${SICER_DISTRIB_DIR}/SICER/${SICER_FILE}
	@echo '	specified SICER path in	${SICER_DISTRIB_DIR}/SICER/${SICER_FILE}'


## NUMPY is required for SICER installation.
##
## NUMPY requires two python libraries Numpy and Scipy
## Info for nose (test library): http://nose.readthedocs.org/en/latest/
## Info for NumPy and scipy: http://www.scipy.org/Installing_SciPy/Mac_OS_X
numpy_and_scipy:
	@echo
	@echo "Installing Python libraries nose, NumPy and scipy with easy_install"
	${SUDO} easy_install nose
	${SUDO} easy_install numpy
	${SUDO} easy_install scipy

################################################################
## SISSRS, peak-calling program. 
## References: 
##  Genome-wide identification of in vivo protein-DNA binding sites from ChIP-Seq data
##    Raja Jothi, Suresh Cuddapah, Artem Barski, Kairong Cui, Keji Zhao.
##    Nucleic Acids Research, 36(16):5221-31, 2008. [Pubmed] [PDF] [Text] [Download Sissrs]
##
##  ChIP-Seq data analysis: identification of protein-DNA binding sites with Sissrs peak-finder
##    Leelavati Narlikar, Raja Jothi.
##    Methods in Molecular Biology, 802:305-22, 2012. [Pubmed] [PDF]
SISSRS_BASE_DIR=${SRC_DIR}/SISSRS
SISSRS_VERSION=1.4
SISSRS_ARCHIVE=sissrs_v${SISSRS_VERSION}.tar.gz
SISSRS_URL=http://dir.nhlbi.nih.gov/papers/lmi/epigenomes/sissrs/${SISSRS_ARCHIVE}
install_sissrs: _download_sissrs _link_sissrs

_download_sissrs:
	@echo
	@echo "Downloading SISSRS"
	@mkdir -p ${SISSRS_BASE_DIR}
	wget -nd  --directory-prefix ${SISSRS_BASE_DIR} -rNL ${SISSRS_URL}

## This installation is VERY tricky. The user has to replace the
## hard-coded path in 3 shell files
_link_sissrs:
	@echo
	@echo "Installing SISSRS in dir	${SISSRS_BASE_DIR}"
	(cd ${SISSRS_BASE_DIR}; tar -xpzf ${SISSRS_ARCHIVE})
	@echo ${SISSRS_BASE_DIR}
	@echo "Linking sissrs in binary dir ${RSAT_BIN}"
	(cd ${RSAT_BIN}; ln -fs ${SISSRS_BASE_DIR}/sissrs.pl sissrs)


################################################################
## Install bfast
## 
## Prerequisites
## 
## BFAST requires the following packages to be installed: 
## - automake (part of GNU autotools) 
## - zlib-dev (ZLIB developer’s libraray)
## - libbz2-dev (BZIP2 developer’s library)
##
## The BZIP2 developer’s library may optionally not be installed, but
## the --disable-bz2 option must be used when running the configure
## script (see the next section).
BFAST_BASE_DIR=${SRC_DIR}/bfast
BFAST_VERSION=0.7.0
BFAST_ARCHIVE=bfast-${BFAST_VERSION}a.tar.gz
BFAST_URL=http://sourceforge.net/projects/bfast/files/bfast/${BFAST_VERSION}/${BFAST_ARCHIVE}
BFAST_DISTRIB_DIR=${BFAST_BASE_DIR}/bfast-${BFAST_VERSION}a
install_bfast: _download_bfast _compile_bfast 

_download_bfast:
	@echo
	@echo "Downloading BFAST"
	@mkdir -p ${BFAST_BASE_DIR}
	wget -nd  --directory-prefix ${BFAST_BASE_DIR} -rNL ${BFAST_URL}

_compile_bfast:
	@echo
	@echo "Installing BFAST in dir	${BFAST_DISTRIB_DIR}"
	(cd ${BFAST_BASE_DIR}; tar -xpzf ${BFAST_ARCHIVE})
	@echo ${BFAST_DISTRIB_DIR}
	(cd  ${BFAST_DISTRIB_DIR}; ./configure ; make; make check; ${SUDO} make install)

################################################################
## Install bowtie, read-mapping program
BOWTIE_BASE_DIR=${SRC_DIR}/bowtie
BOWTIE_VERSION=2.2.5
BOWTIE_ARCHIVE=bowtie2-${BOWTIE_VERSION}-${OS}.zip
BOWTIE_URL=http://sourceforge.net/projects/bowtie-bio/files/bowtie2/${BOWTIE_VERSION}/${BOWTIE_ARCHIVE}
#BOWTIE_URL=http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.5/bowtie2-2.2.5-linux-x86_64.zip/download
BOWTIE_DISTRIB_DIR=${BOWTIE_BASE_DIR}/bowtie2-${BOWTIE_VERSION}
install_bowtie: _download_bowtie _compile_bowtie 

_download_bowtie: _download_bowtie_${OS}

_download_bowtie_macosx:
	${MAKE} _download_bowtie_os OS=macos-x86_64

_download_bowtie_linux:
	${MAKE} _download_bowtie_os OS=linux-x86_64

_download_bowtie_os:
	@echo
	@echo "Downloading BOWTIE"
	@mkdir -p ${BOWTIE_BASE_DIR}
	wget -nd  --directory-prefix ${BOWTIE_BASE_DIR} -rNL ${BOWTIE_URL}

_install_bowtie: _install_bowtie_${OS}

_install_bowtie_macosx:
	${MAKE} _install_bowtie_os OS=macos-x86_64

_install_bowtie_linux:
	${MAKE} _install_bowtie_os OS=linux-x86_64


_install_bowtie_os:
	@echo
	@echo "Installing BOWTIE in dir	${BOWTIE_DISTRIB_DIR}"
	(cd ${BOWTIE_BASE_DIR}; unzip ${BOWTIE_ARCHIVE})
	@echo ${BOWTIE_DISTRIB_DIR}
	@chmod 755  ${BOWTIE_DISTRIB_DIR}/bowtie2*
	${SUDO} rsync -ruptvl ${BOWTIE_DISTRIB_DIR}/bowtie2* ${RSAT_BIN}/

################################################################
## Install  Cis-regulatory Element Annotation System  (CEAS)
CEAS_BASE_DIR=${SRC_DIR}/CEAS
CEAS_VERSION=1.0.2
CEAS_ARCHIVE=CEAS-Package-${CEAS_VERSION}.tar.gz
CEAS_URL=http://liulab.dfci.harvard.edu/CEAS/src/${CEAS_ARCHIVE}
CEAS_DISTRIB_DIR=${CEAS_BASE_DIR}/CEAS-Package-${CEAS_VERSION}
CEAS_COMPILE_DIR=`dirname ${RSAT_BIN}`
#make -f makefiles/install_software.mk RSAT_BIN=/usr/local/bin install_ceas
install_ceas: _download_ceas _compile_ceas 

_download_ceas:
	@echo
	@echo "Downloading CEAS"
	@mkdir -p ${CEAS_BASE_DIR}
	wget -nd  --directory-prefix ${CEAS_BASE_DIR} -rNL ${CEAS_URL}

#CEAS_COMPILE_DIR=/usr/local
_compile_ceas:
	@echo
	@echo "Installing CEAS in dir	${CEAS_DISTRIB_DIR}"
	@echo "CEAS_BASE_DIR	${CEAS_BASE_DIR}"
	@echo "CEAS_DISTRIB_DIR	${CEAS_DISTRIB_DIR}"
	@echo "CEAS_COMPILE_DIR	${CEAS_COMPILE_DIR}"
	(cd ${CEAS_BASE_DIR}; tar -xpzf ${CEAS_ARCHIVE})
	@echo ${CEAS_DISTRIB_DIR}
	(cd  ${CEAS_DISTRIB_DIR}; ${SUDO} python setup.py install --prefix=${CEAS_COMPILE_DIR})
	@echo "CEAS has been installed in dir ${CEAS_COMPILE_DIR}/bin"
	@echo "Before using CEAS, you need to add a line to the log-in shell script"
	@echo "(i.e. .bashrc in case of bash shell)"
	@echo "Adapt the python version in the path below"
	@echo 'export PYTHONPATH=$$PYTHONPATH:${RSAT_BIN}/ext_lib/${PYTHON}/site-packages'

################################################################
## Download data required to run CEAS with Human genome hg19
CEAS_DATA_DIR=${RSAT}/data/ceas
CEAS_GENOME=hg19
download_ceas_data:
	@echo
	@echo "Downloading data required to run CEAS with Human genome ${CEAS_GENOME}"
	mkdir -p ${CEAS_DATA_DIR}/${CEAS_GENOME}
	(cd ${CEAS_DATA_DIR}/${CEAS_GENOME}; \
		wget --no-clobber http://liulab.dfci.harvard.edu/CEAS/src/${CEAS_GENOME}.refGene.gz; \
		gunzip ${CEAS_GENOME}.refGene.gz)


################################################################
## Install  SAMTOOLS
SAMTOOLS_BASE_DIR=${SRC_DIR}/samtools
SAMTOOLS_VERSION=1.3
SAMTOOLS_ARCHIVE=samtools-${SAMTOOLS_VERSION}.tar.bz2
SAMTOOLS_URL=https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2
#SAMTOOLS_URL=https://github.com/samtools/samtools/releases/download/1.3/samtools-1.3.tar.bz2
#SAMTOOLS_URL=https://sourceforge.net/projects/samtools/files/samtools/${SAMTOOLS_VERSION}/${SAMTOOLS_ARCHIVE}/download
#SAMTOOLS_URL=https://sourceforge.net/projects/samtools/files/samtools/1.3/samtools-1.3.tar.bz2/download
SAMTOOLS_DISTRIB_DIR=${SAMTOOLS_BASE_DIR}/samtools-${SAMTOOLS_VERSION}
install_samtools: _download_samtools _compile_samtools _install_pysam

_download_samtools:
	@echo
	@echo "Downloading SAMTOOLS"
	@mkdir -p ${SAMTOOLS_BASE_DIR}
	wget -nd --no-check-certificate --directory-prefix ${SAMTOOLS_BASE_DIR} -rNL ${SAMTOOLS_URL}

_compile_samtools:
	@echo
	@echo "Installing SAMTOOLS in dir	${SAMTOOLS_DISTRIB_DIR}"
	(cd ${SAMTOOLS_BASE_DIR}; tar --bzip2 -xpf ${SAMTOOLS_ARCHIVE})
	@echo ${SAMTOOLS_DISTRIB_DIR}
	(cd ${SAMTOOLS_DISTRIB_DIR}; make)
	${SUDO} find  ${SAMTOOLS_DISTRIB_DIR} -maxdepth 1 -perm 775 -type f  -exec rsync -uptvL {} ${RSAT_BIN}/ \;

## Install a python library required for some samtool functionalities
_install_pysam:
	@echo 
	@echo "Installing pysam library for python"
	sudo pip install pysam


################################################################
## Install  SRA toolkit
SRA_OS=ubuntu64
SRA_BASE_DIR=${SRC_DIR}/sra
SRA_VERSION=2.4.5-2
SRA_ARCHIVE=sratoolkit.${SRA_VERSION}-${SRA_OS}.tar.gz
SRA_URL=ftp://ftp-private.ncbi.nlm.nih.gov/sra/sdk/${SRA_VERSION}/${SRA_ARCHIVE}
SRA_DISTRIB_DIR=${SRA_BASE_DIR}/sratoolkit.${SRA_VERSION}-${SRA_OS}
install_sra: _download_sra_${OS} _install_sra_${OS}

_download_sra_macosx:
	${MAKE} _download_sra SRA_OS=mac64

_download_sra_linux:
	${MAKE} _download_sra SRA_OS=ubuntu64

_download_sra:
	@echo
	@echo "Downloading SRA"
	@mkdir -p ${SRA_BASE_DIR}
	wget -nd  --directory-prefix ${SRA_BASE_DIR} -rNL ${SRA_URL}

_install_sra_macosx:
	${MAKE} _install_sra SRA_OS=mac64

_install_sra_linux:
	${MAKE} _install_sra SRA_OS=ubuntu64

_install_sra:
	@echo
	@echo "Installing SRA in dir	${SRA_DISTRIB_DIR}"
	(cd ${SRA_BASE_DIR}; tar -xpzf ${SRA_ARCHIVE})
	@echo ${SRA_DISTRIB_DIR}
	(cd ${RSAT}/bin; ln -s ${SRA_DISTRIB_DIR}/bin/* . )
#	${SUDO} find  ${SRA_DISTRIB_DIR} -maxdepth 1 -perm 755 -type f  -exec rsync -uptvL {} ${RSAT_BIN}/ \;

################################################################
## Install  SWEMBL
SWEMBL_BASE_DIR=${SRC_DIR}/SWEMBL
SWEMBL_VERSION=3.3.1
SWEMBL_ARCHIVE=SWEMBL.${SWEMBL_VERSION}.tar.bz2
SWEMBL_URL=http://www.ebi.ac.uk/~swilder/SWEMBL/${SWEMBL_ARCHIVE}
SWEMBL_DISTRIB_DIR=${SWEMBL_BASE_DIR}/SWEMBL.${SWEMBL_VERSION}
install_swembl: _download_swembl _compile_swembl 

_download_swembl:
	@echo
	@echo "Downloading SWEMBL"
	@mkdir -p ${SWEMBL_BASE_DIR}
	wget -nd  --directory-prefix ${SWEMBL_BASE_DIR} -rNL ${SWEMBL_URL}

_compile_swembl:
	@echo
	@echo "Installing SWEMBL in dir	${SWEMBL_DISTRIB_DIR}"
	(cd ${SWEMBL_BASE_DIR}; tar --bzip2 -xpf ${SWEMBL_ARCHIVE})
	@echo ${SWEMBL_DISTRIB_DIR}
	(cd ${SWEMBL_DISTRIB_DIR}; make; ${SUDO} rsync -ruptvl ${SWEMBL_DISTRIB_DIR}/SWEMBL ${RSAT_BIN}/)

################################################################
## Internet Genome Browser
IGB_VERSION=5GB
IGB_URL=http://bioviz.org/igb/releases/current/igb-${IGB_VERSION}.jnlp
IGB_BASE_DIR=${SRC_DIR}/IGB
install_igb: _download_igb

_download_igb:
	@echo
	@echo "Downloading IGB"
	@mkdir -p ${IGB_BASE_DIR}
	wget -nd  --directory-prefix ${IGB_BASE_DIR} -rNL ${IGB_URL}

################################################################
## Integrative Genomics Viewer (IGV) tools
IGV_BASE_DIR=${SRC_DIR}/IGV
IGV_VERSION=2.1.7
#IGV_VERSION=nogenomes_2.1.7
IGV_ARCHIVE=igvtools_${IGV_VERSION}.zip
IGV_URL=http://www.broadinstitute.org/igv/projects/downloads/${IGV_ARCHIVE}
IGV_DISTRIB_DIR=${IGV_BASE_DIR}/IGVTools
install_igv: _download_igv _compile_igv 

_download_igv:
	@echo
	@echo "Downloading IGV"
	@mkdir -p ${IGV_BASE_DIR}
	wget -nd  --directory-prefix ${IGV_BASE_DIR} -rNL ${IGV_URL}

_compile_igv:
	@echo
	@echo "Installing IGV in dir	${IGV_DISTRIB_DIR}"
	(cd ${IGV_BASE_DIR}; unzip ${IGV_ARCHIVE})
	@echo ${IGV_DISTRIB_DIR}
	@echo "Please add IGVTools folder to your path"
	@echo 'export PATH=$${PATH}:${IGV_DISTRIB_DIR}'

################################################################
## BLAT
## 
## Fast sequence similarity search program developed by Jim Kent to
## identify homogous segments of genomes. BLAT is required for Homer,
## which uses it to purge sequences from redundant fragments before
## motif discovery.  
##
## Info on BLAT:https://genome.ucsc.edu/FAQ/FAQblat.html
BLAT_VERSION=35
BLAT_ARCHIVE=blatSrc${BLAT_VERSION}.zip
BLAT_URL=https://users.soe.ucsc.edu/~kent/src/${BLAT_ARCHIVE}
BLAT_DOWNLOAD_DIR=${SRC_DIR}/BLAT
install_blat: _download_blat

_download_blat:
	@echo
	@echo "Downloading BLAT sources from	${BLAT_URL}"
	mkdir -p ${BLAT_DOWNLOAD_DIR}
	wget --no-clobber --no-directories --directory-prefix ${BLAT_DOWNLOAD_DIR} ${BLAT_URL}
	(cd ${BLAT_DOWNLOAD_DIR}; unzip -o  ${BLAT_ARCHIVE})
	@echo "BLAT sources downloaded in ${BLAT_DOWNLOAD_DIR}"

BLAT_SRC_DIR=${BLAT_DOWNLOAD_DIR}/blatSrc
_compile_blat:
	@echo
	@echo "Compiling BLAT sources in ${BLAT_SRC_DIR}"
	(cd ${BLAT_SRC_DIR}; make)
	@echo "BLAT executables are in ${BLAT_DOWNLOAD_DIR}"



################################################################
## HOMER
HOMER_CONFIG_URL=http://homer.salk.edu/homer/configureHomer.pl
#HOMER_CONFIG_URL=http://biowhat.ucsd.edu/homer/configureHomer.pl
HOMER_BASE_DIR=${SRC_DIR}/HOMER
install_homer: _download_homer _install_homer

_download_homer:
	@echo 
	@echo "Downloading HOMER"
	@mkdir -p ${HOMER_BASE_DIR}
	wget -nd  --directory-prefix ${HOMER_BASE_DIR} -rNL ${HOMER_CONFIG_URL}
	@echo "	${HOMER_BASE_DIR}"

_install_homer:
	@echo
	@echo "Installing HOMER in dir	${HOMER_BASE_DIR}"
	(cd ${HOMER_BASE_DIR}; perl ./configureHomer.pl -install)
	@echo "HOMER installed in dir	${HOMER_BASE_DIR}"
	@echo "Please add the three following lines to your .bashrc file in order to include HOMER programs in your path"
#	@echo 'export PATH=$${PATH}:${HOMER_BASE_DIR}/bin'
	@echo 'export HOMER=${HOMER_BASE_DIR}'
	@echo 'export PATH=$${PATH}:\$$HOMER/bin'
	@echo 'export PERL5LIB=$${PERL5LIB}:$$HOMER/bin'


HOMER_GENOME=mm9
install_homer_one_genome:
	@echo
	@echo "Installing HOMER genome	HOMER_GENOME=${HOMER_GENOME}"
	(cd ${HOMER_BASE_DIR}; perl ./configureHomer.pl -install ${HOMER_GENOME})

HOMER_GENOMES=dm3 mm8 mm9 hg18 hg19 rn4 rn5
install_homer_some_genomes:
	@echo
	@echo "Installing genomes for HOMER"
	@echo "	${HOMER_GENOMES}"
	@for g in ${HOMER_GENOMES}; do \
		${MAKE} install_homer_one_genome HOMER_GENOME=$${g} ; \
	done

################################################################
## clustalW (multiple alignment)
install_clustalw: _download_clustalw _compile_clustalw

CLUSTALW_BASE_DIR=${SRC_DIR}/clustalw
CLUSTALW_ARCHIVE=clustalw-${CLUSTALW_VERSION}.tar.gz
CLUSTALW_VERSION=2.1
CLUSTALW_URL=http://www.clustal.org/download/current/${CLUSTALW_ARCHIVE}
CLUSTALW_SOURCE_DIR=clustalw_latest
_download_clustalw:
	@mkdir -p ${CLUSTALW_BASE_DIR}
	wget --no-directories  --directory-prefix ${CLUSTALW_BASE_DIR} -rNL ${CLUSTALW_URL} -A "${CLUSTALW_ARCHIVE}"
	(cd ${CLUSTALW_BASE_DIR}; tar -xzf ${CLUSTALW_ARCHIVE})
	@echo ${CLUSTALW_BASE_DIR}

_compile_clustalw:
	(cd ${CLUSTALW_BASE_DIR}/${CLUSTALW_SOURCE_DIR}; ./configure; make clean ; make ; \
	${SUDO} rsync -ruptvl src/clustalw2 ${RSAT_BIN}/)
	@echo "	clustalw2 should now be executable from ${RSAT_BIN}";


################################################################
## Obsolete: compile some perl scripts to binaries.  This was a test and the
## results were not very good, the compiled programs were unstable.
SRC=perl-scripts
COMPIL=compil/
PROGRAMS=	\
	oligo-analysis	\
	retrieve-seq	\
	dyad-analysis	\
	dna-pattern	\
	orf-info
LIBRARIES=\
	RSA.stat.lib 	\
#	RSA.classes	\
#	RSA.seq.lib	\
#	RSA.cgi.lib	\
#	RSA.lib 	
_compile_perl_scripts:
	@mkdir -p ${COMPIL}/lib
	@mkdir -p ${COMPIL}/bin
	@(cd  ${COMPIL}/bin; ln -fs ../lib)
	@cp -f config/default.config ${COMPIL}/RSA.config

	@for lb in ${LIBRARIES}; do \
		echo "compiling library $${lb}"; \
		cp -f ${SRC}/lib/$${lb} ${COMPIL}/lib/$${lb}.pl ; \
		(cd ${COMPIL}/lib; pwd; perlcc $${lb}.pl && rm -f $${lb}.pl); \
	done

	@for pgm in ${PROGRAMS}; do \
		echo "compiling program $${pgm}"; \
		cp -f ${SRC}/$${pgm} ${COMPIL}/bin/$${pgm}.pl ; \
		(cd ${COMPIL}/bin; pwd; perlcc $${pgm}.pl && rm -f $${pgm}.pl); \
	dgone



################################################################
## This library allows you to install the Perl libraries locally, if you are not system administrator
LOCAL_LIB_URL=http://search.cpan.org/CPAN/authors/id/A/AP/APEIRON/local-lib-1.008004.tar.gz
LOCAL_LIB_DIR=lib/perl_lib/locallib

local_lib: download_local_lib install_local_lib config_local_lib

download_local_lib:
	@echo "Downloading Perl module local::lib"
	(mkdir -p ${LOCAL_LIB_DIR}; cd ${LOCAL_LIB_DIR}; wget ${LOCAL_LIB_URL}; tar -xzf local-lib-1.008004.tar.gz)

install_local_lib:
	@echo "Installing Perl module local::lib"
	(mkdir -p ${RSAT}/ext_lib/perl5; ln -s ${RSAT}/ext_lib/perl5 ${HOME}/perl5)
	(cd ${LOCAL_LIB_DIR}/local-lib-1.008004;  perl Makefile.PL --bootstrap; make; make test; make install)

config_local_lib:
	@echo "Adding path to Perl module local::lib in ${HOME}/.bashrc"
	@echo ''  >>~/.bashrc
	@echo '################################################################'  >>~/.bashrc
	@echo '## Perl local::lib module'  >>~/.bashrc
	@echo 'eval $$(perl -I$$HOME/perl5/ext_lib/perl5 -Mlocal::lib)' >>~/.bashrc

#install_one_perl_module_locally:
#	${MAKE} SUDO='' install_one_perl_module
#
#install_perl_modules_locally:
#	${MAKE} SUDO='' install_perl_modules

################################################################
## Install Python 2.7.  We deliberately chose version 2.7 (and not
## version 3.x) because some modules are not working with version 3.x.
PYTHON_COMPILE_DIR=${SRC_DIR}/python
install_python: _download_python _compile_python

_download_python:
	@echo "Downloading Python-2.7 to dir ${PYTHON_COMPILE_DIR}"
	(mkdir -p ${PYTHON_COMPILE_DIR}; cd ${PYTHON_COMPILE_DIR}; wget -NL http://www.python.org/ftp/python/2.7/Python-2.7.tgz; tar -xpzf Python-2.7.tgz)

_compile_python:
	@echo "Compiling python2.7"
	(cd ${PYTHON_COMPILE_DIR}/Python-2.7; ./configure; make; ${SUDO} make install)

install_python_suds: _download_python_suds _compile_python_suds

################################################################
## Install suds library for python2.7, required for the MICROME Web
## clients to connect Genoscope/Microscope Web services.
SUDS_VERSION=0.4
SUDS_TAR=python-suds-${SUDS_VERSION}.tar.gz
SUDS_URL=https://fedorahosted.org/releases/s/u/suds/${SUDS_ARCHIVE}
SUDS_DIR=${SRC_DIR}/suds
_download_python_suds:
	@mkdir -p ${SUDS_DIR}
	@echo "Getting suds (Python library) using wget"
	(cd ${SUDS_DIR}; wget -nv -nd ${SUDS_URL}/${SUDS_TAR}; tar -xpzf ${SUDS_TAR})
	@echo "suds dir	${SUDS_DIR}"

SUDS_INSTALL_DIR=${SUDS_DIR}/python-suds-${SUDS_VERSION}
_compile_python_suds:
	@echo "Installing suds"
	(cd ${SUDS_INSTALL_DIR}; python2.7 setup.py build; ${SUDO} python2.7 setup.py install)


################################################################
## Install STAMP (zip archive kindly sent by email by Shaun Mahony)
install_stamp: _clone_stamp _compile_stamp _install_stamp

_clone_stamp:
	@echo "Cloning stamp"
	(cd ${SRC_DIR}; git clone https://github.com/shaunmahony/stamp.git)

_compile_stamp:
	@echo "Compiling stamp"
	(cd ${SRC_DIR}/stamp/src; \
		g++ -O3 -o stamp Motif.cpp Alignment.cpp ColumnComp.cpp \
                PlatformSupport.cpp PlatformTesting.cpp Tree.cpp \
                NeuralTree.cpp MultipleAlignment.cpp RandPSSMGen.cpp \
                ProteinDomains.cpp main.cpp -lm -lgsl -lgslcblas)

_install_stamp:
	${SUDO} rsync -ruptvl ${SRC_DIR}/stamp/src/stamp ${RSAT_BIN}


################################################################
## Install MATLIGN

################################################################
## Get and install the program matlign
MATLIGN_ARCHIVE=matlign.tgz
MATLIGN_URL=http://ekhidna.biocenter.helsinki.fi/poxo/download_folder
MATLIGN_DOWNLOAD_DIR=${SRC_DIR}/matlign/
MATLIGN_COMPILE_DIR=${MATLIGN_DOWNLOAD_DIR}data/backup/zope_data/poxo/MATLIGN
install_matlign: _download_matlign _compile_matlign _install_matlign

_download_matlign:
	@mkdir -p ${MATLIGN_DOWNLOAD_DIR}
	@echo
	@echo "Downloading matlign	${MATLIGN_URL}"
	(cd ${MATLIGN_DOWNLOAD_DIR}; wget --timestamping ${MATLIGN_URL}/${MATLIGN_ARCHIVE}; tar -xpzf ${MATLIGN_ARCHIVE})
	@echo "MATLING_DOWNLOAD_DIR	${MATLIGN_DOWNLOAD_DIR}"
	@echo "MATLING_COMPILE_DIR	${MATLIGN_COMPILE_DIR}"

_compile_matlign:
	@echo
	@echo "Compiling matlign in MATLIGN_DIR	${MATLIGN_DIR}"
	(cd ${MATLIGN_COMPILE_DIR}; ./compile1)


_install_matlign:
	@echo
	@echo "Installing matlign in RSAT_BIN	${RSAT_BIN}"
	@${SUDO} rsync -ruptl ${MATLIGN_COMPILE_DIR}/matlign ${RSAT_BIN}/


################################################################
## plink - polymorphism linkage analysis
PLINK_RELEASE=160705
PLINK_URL=http://www.cog-genomics.org/static/bin/plink${PLINK_RELEASE}
PLINK_ARCHIVE_linux=plink_linux_x86_64.zip
PLINK_ARCHIVE_macosx=plink_mac.zip
PLINK_ARCHIVE=${PLINK_ARCHIVE_${OS}}
PLINK_MAC_URL=${PLINK_URL}/${PLINK_ARCHIVE}
PLINK_DIR=${SRC_DIR}/plink
install_plink:
	@echo
	@echo "Installing plink for operating system ${OS}"
	@echo "	PLINK_URL	${PLINK_URL}"
	@echo "	PLINK_DIR	${PLINK_DIR}"
	@echo "	PLINK_ARCHIVE	${PLINK_ARCHIVE}"
	${MAKE} _download_plink _install_plink

_download_plink:
	@mkdir -p ${PLINK_DIR}
	@echo "Downloading plink using wget"
	(cd ${PLINK_DIR}; wget -nv -nd ${PLINK_URL}/${PLINK_ARCHIVE}; unzip ${PLINK_ARCHIVE})
	@echo "plink dir	${PLINK_DIR}"
	@echo "plink zip	${PLINK_DIR}/${PLINK_ARCHIVE}"

_install_plink:
	@echo "Uncompressing PLINK_ACHIVE	${PLINK_DIR}/${PLINK_ARCHIVE}"
	@echo "Installing in RSA_BIN	${RSAT_BIN}"
	(cd ${PLINK_DIR}; rsync -ruptvl plink ${RSAT_BIN}/; rsync -ruptvl prettify ${RSAT_BIN}/)

#}; ${SUDO} rsync -ruptvl ${PLINK_BIN} ${RSAT_BIN}/; cd ${RSAT_BIN}; ${SUDO} rm -f plink; ${SUDO} ln -s ${PLINK_BIN} plink)

