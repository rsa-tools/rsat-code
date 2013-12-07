############################################################
#
# $Id: install_software.mk,v 1.54 2013/11/03 19:35:38 jvanheld Exp $
#
# Time-stamp: <2003-05-23 09:36:00 jvanheld>
#
############################################################


################################################################
## This makefile manages the installation of bioinformatics software
## and the configuration of the paths and environment variables for
## the users.

include ${RSAT}/makefiles/util.mk
MAKEFILE=makefiles/install_software.mk
MAKE=make -f ${MAKEFILE}
V=1

#################################################################
# Programs used for downloading and sycnrhonizing
WGET = wget -np -rNL 
#MAKE=nice -n 19 make -s -f ${MAKEFILE}
RSYNC_OPT = -ruptvl ${OPT}
SSH=-e 'ssh -x'
RSYNC = rsync ${RSYNC_OPT} ${SSH}

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
	(mkdir -p ${RSAT}/lib/perl5; ln -s ${RSAT}/lib/perl5 ${HOME}/perl5)
	(cd ${LOCAL_LIB_DIR}/local-lib-1.008004;  perl Makefile.PL --bootstrap; make; make test; make install)

config_local_lib:
	@echo "Adding path to Perl module local::lib in ${HOME}/.bashrc"
	@echo ''  >>~/.bashrc
	@echo '################################################################'  >>~/.bashrc
	@echo '## Perl local::lib module'  >>~/.bashrc
	@echo 'eval $$(perl -I$$HOME/perl5/lib/perl5 -Mlocal::lib)' >>~/.bashrc

install_one_perl_module_locally:
	${MAKE} SUDO='' install_one_perl_module

install_perl_modules_locally:
	${MAKE} SUDO='' install_perl_modules

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
## Install the applications developed by third-parties and which are required
## or useful for RSAT.
install_ext_apps:
	${MAKE} download_seqlogo install_seqlogo
	${MAKE} bedtools
	${MAKE} download_meme install_meme
	${MAKE} download_mcl install_mcl
	${MAKE} download_rnsc install_rnsc
#	${MAKE} download_blast install_blast
#	${MAKE} download_gs install_gs
#	${MAKE} download_gnuplot install_gnuplot
#	${MAKE} install_gibbs
#	${MAKE} download_consensus install_consensus
#	${MAKE} download_patser install_patser

# ################################################################
# ## Install the EnsEMBL Perl API
# ENSEMBL_BRANCH=74
# ensembl_api:	
# 	@echo  "Password is 'CVSUSER'"
# 	@cvs -d :pserver:cvsuser@cvs.sanger.ac.uk:/cvsroot/ensembl login
# 	@(cd ${RSAT}/lib; \
# 		cvs -d :pserver:cvsuser@cvs.sanger.ac.uk:/cvsroot/ensembl \
# 		checkout -r branch-ensembl-${ENSEMBL_BRANCH} ensembl ; \
# 		cvs -d :pserver:cvsuser@cvs.sanger.ac.uk:/cvsroot/ensembl \
# 		checkout -r branch-ensembl-${ENSEMBL_BRANCH} ensembl-compara)
# 	@echo "Don't forget to adapt the following lines in the file ${RSAT}/RSAT_config.props"
# 	@echo "ensembl=${RSAT}/lib/ensembl/modules"
# 	@echo "compara=${RSAT}/lib/ensembl-compara/modules"


################################################################
## Install the EnsEMBL Perl API
##
## Important note
## (http://bacteria.ensembl.org/info/docs/api/api_cvs.html): you must
## install version 1.2.3, not a more recent version. Starting with
## 1.2.4, major changes were made to the BioPerl API which have made
## it incompatible with Ensembl
BIOPERL_VERSION=1-2-3
BIOPERL_DIR=${PERLLIB_DIR}/bioperl-release-${BIOPERL_VERSION}

## Note: In some cases, there are delays between Ensembl and
## EnsemblGenome releases. To ensure compatibility, please check the
## versions of both distributions on http://ensemblgenomes.org/.
ENSEMBL_VERSION=73
ENSEMBL_BRANCH=ensemblgenomes-20
ENSEMBL_API_DIR=${PERLLIB_DIR}/${ENSEMBL_BRANCH}-${ENSEMBL_VERSION}
install_ensembl_api_param:
	@echo "	BIOPERL_VERSION		${BIOPERL_VERSION}"
	@echo "	BIOPERL_DIR		${BIOPERL_DIR}"
	@echo "	ENSEMBL_VERSION		${ENSEMBL_VERSION}"
	@echo "	ENSEMBL_BRANCH		${ENSEMBL_BRANCH}"
	@echo "	ENSEMBL_API_DIR		${ENSEMBL_API_DIR}"

## Install the required modules for Ensembl API
install_ensembl_api:
	@(cd ${BIOPERL_DIR}; \
		echo "" ; \
		echo "Installing ensembl branch ${ENSEMBL_BRANCH} version ${ENSEMBL_VERSION}"; \
		echo "	ENSEMBL_API_DIR		${ENSEMBL_API_DIR}" ; \
		mkdir -p "${ENSEMBL_API_DIR}"; \
		cd ${ENSEMBL_API_DIR}; \
		echo  "Password is 'CVSUSER'" ; \
		cvs -d :pserver:cvsuser@cvs.sanger.ac.uk:/cvsroot/ensembl login ; \
		cvs -d :pserver:cvsuser@cvs.sanger.ac.uk:/cvsroot/ensembl checkout -r branch-${ENSEMBL_BRANCH}-${ENSEMBL_VERSION} ensembl-api)
		echo "" ; \
		echo "Installing bioperl release ${BIOPERL_VERSION} (required for ensembl)"; \
		echo "	BIOPERL_DIR		${BIOPERL_DIR}" ; \
		mkdir -p "${BIOPERL_DIR}" ; \
		echo  "Password is 'cvs'" ; \
		cvs -d :pserver:cvs@code.open-bio.org:/home/repository/bioperl login ; \
		cvs -d :pserver:cvs@code.open-bio.org:/home/repository/bioperl checkout -r bioperl-release-1-2-3 bioperl-live ; \
	@echo
	@${MAKE} install_ensembl_api_env

## Print the settings for including ensembl in the PERL5LIB environment variable
install_ensembl_api_env:
	@echo
	@echo "ENSEMBL Perl modules are installed in directory ${ENSEMBL_API_DIR}"
	@echo "Don't forget to adapt the following lines in your bash profile"
	@echo 'export PERL5LIB=${ENSEMBL_API_DIR}/ensembl/modules::$${PERL5LIB}'
	@echo 'export PERL5LIB=${ENSEMBL_API_DIR}/ensembl-compara/modules::$${PERL5LIB}'
	@echo 'export PERL5LIB=${ENSEMBL_API_DIR}/ensembl-external/modules::$${PERL5LIB}'
	@echo 'export PERL5LIB=${ENSEMBL_API_DIR}/ensembl-functgenomics/modules::$${PERL5LIB}'
	@echo 'export PERL5LIB=${ENSEMBL_API_DIR}/ensembl-tools/modules::$${PERL5LIB}'
	@echo 'export PERL5LIB=${ENSEMBL_API_DIR}/ensembl-variation/modules::$${PERL5LIB}'
	@echo 'export PERL5LIB=${BIOPERL_DIR}/bioperl-live::$${PERL5LIB}'

################################################################
## Get and install the program seqlogo
SEQLOGO_URL=http://weblogo.berkeley.edu/release
SEQLOGO_TAR=weblogo.2.8.2.tar.gz
SEQLOGO_DIR=${SRC_DIR}/seqlogo
install_seqlogo: _download_seqlogo _compile_seqlogo

_download_seqlogo:
	@mkdir -p ${SEQLOGO_DIR}
	@echo "Getting seqlogo using wget"
	(cd ${SEQLOGO_DIR}; wget -nv -nd ${SEQLOGO_URL}/${SEQLOGO_TAR}; tar -xpzf ${SEQLOGO_TAR})
	@echo "seqlogo dir	${SEQLOGO_DIR}"

_compile_seqlogo:
	@echo "Installing seqlogo"
	@${SUDO} rsync -ruptl ${SEQLOGO_DIR}/weblogo/seqlogo ${BIN_DIR}
	@${SUDO} rsync -ruptl ${SEQLOGO_DIR}/weblogo/template.* ${BIN_DIR}
	@${SUDO} rsync -ruptl ${SEQLOGO_DIR}/weblogo/logo.pm ${BIN_DIR}

################################################################
## Get and install the program ghostscript
## Note: for Mac users, please go to the ghostscript Web site
##
## GS_URL=http://ghostscript.com/releases/
## GS_VER=ghostscript-8.64
## GS_TAR=${GS_VER}.tar.gz
GS_URL=http://downloads.ghostscript.com/public/binaries/
GS_VER=ghostscript-9.06-linux-x86_64
GS_TAR=${GS_VER}.tgz
GS_DIR=${SRC_DIR}/ghostscript
install_ghostscript: _download_gs _compile_gs

_download_gs:
	@mkdir -p ${GS_DIR}
	@echo "Getting gs using wget"
	(cd ${GS_DIR}; wget -nv -nd ${GS_URL}/${GS_TAR}; tar -xpzf ${GS_TAR})
	@echo "gs dir	${GS_DIR}"

_compile_gs:
	@echo "Installing gs"
	(cd ${GS_DIR}/${GS_VER}; ./configure && make)

################################################################
## Get and install the program gnuplot
GNUPLOT_VER=4.6.0
GNUPLOT_TAR=gnuplot-${GNUPLOT_VER}.tar.gz
GNUPLOT_URL=http://sourceforge.net/projects/gnuplot/files/gnuplot/${GNUPLOT_VER}/${GNUPLOT_TAR}
GNUPLOT_DIR=${SRC_DIR}/gnuplot
install_gnuplot: _download_gnuplot _compile_gnuplot

_download_gnuplot:
	@mkdir -p ${GNUPLOT_DIR}
	@echo "Getting gnuplot using wget"
	(cd ${GNUPLOT_DIR}; wget -nv -nd ${GNUPLOT_URL}; tar -xpzf ${GNUPLOT_TAR})
	@echo "gnuplot dir	${GNUPLOT_DIR}"

_compile_gnuplot:
	@echo "Installing gnuplot"
	(cd ${GNUPLOT_DIR}/gnuplot-${GNUPLOT_VER}; ./configure && make; ${SUDO} make install)


################################################################
## Install fastqc, a software tool to control the quality of read
## files (next generation sequencing).
FASTQC_VER=0.10.1
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
	(cd ${FASTQC_DOWNLOAD_DIR}; wget -nv -nd ${FASTQC_URL}; unzip ${FASTQC_ZIP})
	@echo "	fastqc download dir	${FASTQC_DOWNLOAD_DIR}"
	@echo "	fastqc install dir	${FASTQC_INSTALL_DIR}"
	@chmod 755 ${FASTQC}
	@echo "	fastqc executable dir	${FASTQC_EXEC_DIR}"
	@echo "	fastqc script    	${FASTQC}"
	@echo
	@echo "YOU NEED TO ADD THE FASTQC EXEC DIRECTORY TO YOUR PATH"
	@echo "export PATH=\$$PATH:${FASTQC_EXEC_DIR}"

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
BED_VERSION=2.17.0
BED_ARCHIVE=BEDTools.v${BED_VERSION}.tar.gz
BED_URL=http://bedtools.googlecode.com/files/${BED_ARCHIVE}
BEDTOOL_MANUAL=http://bedtools.googlecode.com/files/BEDTools-User-Manual.v4.pdf
BED_BASE_DIR=${SRC_DIR}/BEDTools
BED_DOWNLOAD_DIR=${BED_BASE_DIR}/bedtools-${BED_VERSION}
_download_bedtools:
	@echo
	@echo "Downloading BEDTools ${BED_VERSION}"
	@echo
	@mkdir -p ${BED_BASE_DIR}
	(cd ${BED_BASE_DIR}; wget -nv -nd ${BED_URL} ; tar -xpzf ${BED_ARCHIVE}; \
		wget ${BEDTOOL_MANUAL})
	@echo ${BED_DOWNLOAD_DIR}

BED_GIT_DIR=${SRC_DIR}/bedtools
_git_bedtools:
	@mkdir -p ${BED_GIT_DIR}
	(cd ${SRC_DIR}; git clone git://github.com/arq5x/bedtools.git)

#BED_SRC_DIR=${BED_GIT_DIR}
BED_SRC_DIR=${BED_DOWNLOAD_DIR}
BED_BIN_DIR=${BED_SRC_DIR}/bin
_compile_bedtools:
	@echo
	@echo "Installing bedtools from ${BED_SRC_DIR}"
	@echo
	@mkdir -p ${BED_SRC_DIR}
	(cd ${BED_SRC_DIR}; make clean; make)

_install_bedtools:
	@echo
	@echo "Installing bedtools binaries from ${BEN_BIN_DIR} to ${BIN_DIR}"
	@echo
	@mkdir -p ${BIN_DIR}
	@${SUDO} rsync -ruptvl ${BED_BIN_DIR}/* ${BIN_DIR}

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
MEME_VERSION=4.9.0
MEME_PATCH=_4
#MEME_VERSION=4.8.0
#MEME_PATCH=
#MEME_VERSION=current
MEME_ARCHIVE=meme_${MEME_VERSION}${MEME_PATCH}.tar.gz
##MEME_URL=http://meme.nbcr.net/downloads/${MEME_ARCHIVE}
#MEME_URL=ftp://ftp.ebi.edu.au/pub/software/MEME/r${MEME_VERSION}/rc5/${MEME_ARCHIVE}
#MEME_URL=http://ebi.edu.au/ftp/software/MEME/${MEME_VERSION}/${MEME_ARCHIVE}
MEME_URL=http://ebi.edu.au/ftp/software/MEME/${MEME_VERSION}/${MEME_ARCHIVE}
MEME_INSTALL_SUBDIR=${SOFT_DIR}/MEME
MEME_INSTALL_DIR=${MEME_INSTALL_SUBDIR}/meme_${MEME_VERSION}
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
	@echo "Installing MEME ${MEME_VERSION} in dir ${MEME_INSTALL_DIR}"
	@mkdir -p ${MEME_INSTALL_DIR}
	(cd ${MEME_INSTALL_SUBDIR}; tar -xpzf ${MEME_BASE_DIR}/${MEME_ARCHIVE})
#	@echo "MEME configuration prefix	${MEME_CONFIG_PREFIX}"
	(cd ${MEME_INSTALL_DIR}; ./configure --prefix=${MEME_COMPILE_DIR} --with-url="http://localhost/meme")
	(cd ${MEME_INSTALL_DIR}; make clean; make ; make test; ${SUDO} make install)
	@echo "MEME installed in ${MEME_COMPILE_DIR}"

_after_meme:
#	@cd ${MEME_BIN_DIR}; rm -f meme; ln -s meme_${MEME_VERSION} meme
	@echo "Please edit the bashrc file"
	@echo "and copy-paste the following lines to specify the MEME bin pathway"
	@echo "	export PATH=${MEME_BIN_DIR}:\$$PATH"

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

CLUSTER_COMPILE_DIR=`dirname ${BIN_DIR}`
CLUSTER_BIN_DIR=${CLUSTER_COMPILE_DIR}/bin
_compile_cluster:
	@echo
	@echo "Installing CLUSTER in dir	${CLUSTER_COMPILE_DIR}"
	@mkdir -p ${CLUSTER_COMPILE_DIR}
	(cd ${CLUSTER_DISTRIB_DIR}; ./configure --without-x --prefix=${CLUSTER_COMPILE_DIR} ; \
	make clean; make ; ${SUDO} make install)

################################################################
## Install the graph-based clustering algorithm MCL
MCL_BASE_DIR=${SRC_DIR}/mcl
MCL_VERSION=12-135
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

MCL_COMPILE_DIR=`dirname ${BIN_DIR}`
MCL_BIN_DIR=${BIN_DIR}
_compile_mcl:
	@echo
	@echo "Installing MCL"
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
	@echo "Installing RNSC"
	@${SUDO} mkdir -p ${BIN_DIR}
	(cd ${RNSC_BASE_DIR}; make ;  \
	${SUDO} rsync -ruptvl rnsc ${BIN_DIR}; \
	${SUDO} rsync -ruptvl rnscfilter ${BIN_DIR}; \
	)
	@echo "Please check that RNSC bin directory in your path."
	@echo "	${BIN_DIR}"
#	${SUDO} rsync -ruptvl rnscconvert ${BIN_DIR}; \

################################################################
## Install BLAST
install_blast: _list_blast_param _download_blast _install_blast

_download_blast: _download_blast_${OS}

#_install_blast: _install_blast_${OS}
_install_blast:
	@${SUDO} mkdir -p ${BIN_DIR}
	@echo "Please edit the RSAT configuration file"
	@echo "	${RSAT}/RSAT_config.props"
	@echo "and copy-paste the following line to specify the BLAST bin pathway"
	@echo "	blast_dir=${BIN_DIR}"
	@echo "This will allow RSAT programs to idenfity BLAST path on this server."
	@echo
	@echo "You can also add the BLAST bin directory in your path."
	@echo "If your shell is bash"
	@echo "	export PATH=${BIN_DIR}:\$$PATH"
	@echo "If your shell is csh or tcsh"
	@echo "	setenv PATH ${BIN_DIR}:\$$PATH"

_list_blast_param:
	@echo "Downloading blast"
	@echo "	BLAST_BASE_DIR		${BLAST_BASE_DIR}"
	@echo "	BLAST_LINUX_ARCHIVE	${BLAST_LINUX_ARCHIVE}"
	@echo "	BLAST_URL		${BLAST_URL}"
	@echo "	BLAST_SOURCE_DIR	${BLAST_SOURCE_DIR}"
	@echo "	BLAST_BASE_DIR		${BLAST_BASE_DIR}"
	@echo "	BIN_DIR			${BIN_DIR}"

################################################################
## Install the BLAST on linux
BLAST_BASE_DIR=${SRC_DIR}/blast
BLAST_LINUX_ARCHIVE=blast-*-${ARCHITECTURE}-linux.tar.gz
BLAST_URL=ftp://ftp.ncbi.nih.gov/blast/executables/release/LATEST/
BLAST_SOURCE_DIR=blast_latest
_download_blast_linux:
	@mkdir -p ${BLAST_BASE_DIR}
	wget --no-directories  --directory-prefix ${BLAST_BASE_DIR} -rNL ${BLAST_URL} -A "${BLAST_LINUX_ARCHIVE}"
	(cd ${BLAST_BASE_DIR}; tar -xvzf ${BLAST_LINUX_ARCHIVE}; \
		${SUDO} rm -rf ${BLAST_SOURCE_DIR}; \
		${SUDO} mv -f ${BLAST_LINUX_ARCHIVE} ..; \
		${SUDO} mv -f blast-*  ${BLAST_SOURCE_DIR} \
	)
	@echo ${BLAST_BASE_DIR}

#_install_blast_linux:
#	@${SUDO} mkdir -p ${BIN_DIR}
#	${SUDO} rsync -ruptvl ${BLAST_BASE_DIR}/${BLAST_SOURCE_DIR}/bin/blastall ${BIN_DIR}
#	${SUDO} rsync -ruptvl ${BLAST_BASE_DIR}/${BLAST_SOURCE_DIR}/bin/formatdb ${BIN_DIR}
#	@echo "Please check that the BLAST install directory is in your path"
#	@echo "	${BIN_DIR}"

################################################################
## Install the BLAST on MAC
BLAST_BASE_DIR=${SRC_DIR}/blast
BLAST_MAC_ARCHIVE=blast-*-universal-macosx.tar.gz
BLAST_URL=ftp://ftp.ncbi.nih.gov/blast/executables/release/LATEST/
BLAST_SOURCE_DIR=blast_latest
_download_blast_macosx:
	@mkdir -p ${BLAST_BASE_DIR}
	wget --no-directories  --directory-prefix ${BLAST_BASE_DIR} -rNL ${BLAST_URL} -A "${BLAST_MAC_ARCHIVE}"
	(cd ${BLAST_BASE_DIR}; tar -xvzf ${BLAST_MAC_ARCHIVE}; \
		${SUDO} rm -rf ${BLAST_SOURCE_DIR}; \
		${SUDO} mv -f ${BLAST_MAC_ARCHIVE} ..; \
		${SUDO} mv -f blast-*  ${BLAST_SOURCE_DIR})
	@echo ${BLAST_BASE_DIR}

# _install_blast_macosx:
# 	@${SUDO} mkdir -p ${BIN_DIR}
# 	${SUDO} rsync -ruptvl ${BLAST_BASE_DIR}/${BLAST_SOURCE_DIR}/bin/blastall ${BIN_DIR}
# 	${SUDO} rsync -ruptvl ${BLAST_BASE_DIR}/${BLAST_SOURCE_DIR}/bin/formatdb ${BIN_DIR}

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
	(cd ${PROGRAM_DIR}; make ${_COMPILE_OPT})
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
PATSER_URL=ftp://www.genetics.wustl.edu/pub/stormo/Consensus
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
	@echo "Installing patser"
	(cd ${PATSER_DIR}; rm *.o; make)
	${SUDO} rsync -ruptvl ${PATSER_DIR}/${PATSER_APP} ${BIN_DIR}
	(cd ${BIN_DIR}; ${SUDO} ln -fs ${PATSER_APP} patser)
	@echo "ls -ltr ${BIN_DIR}/patser*"


################################################################
## Install consensus (J.Hertz)
CONSENSUS_VERSION=consensus-v6c.1
CONSENSUS_TAR=${CONSENSUS_VERSION}.tar.gz
CONSENSUS_URL=ftp://www.genetics.wustl.edu/pub/stormo/Consensus
CONSENSUS_DIR=ext/consensus/${CONSENSUS_VERSION}
_download_consensus:
	@echo
	@echo "Downloading ${CONSENSUS_VERSION}"
	@echo
	@mkdir -p ${CONSENSUS_DIR}
	(cd ${CONSENSUS_DIR}; wget -nv -nd ${CONSENSUS_URL}/${CONSENSUS_TAR}; tar -xpzf ${CONSENSUS_TAR})
	@echo "consensus dir	${CONSENSUS_DIR}"

_compile_consensus:
#	${MAKE} uncompress_program PROGRAM=consensus
	${MAKE} _compile_program PROGRAM=consensus _COMPILE_OPT='CPPFLAGS=""'



################################################################
## UCSC tools (developed by Jim Kent)
## APPARENTLY THIS DOES NOT WORK YET
install_ucsc_tools: _download_ucsc_tools
UCSC_URL=http://hgdownload.cse.ucsc.edu/admin/exe/${UCSC_OS}/
UCSC_SRC_DIR=${SRC_DIR}/UCSC
_download_ucsc_tools:
	(mkdir -p ${UCSC_SRC_DIR}; cd ${SRC_DIR}/ucsc/;  wget  -nv -nd -rNL ${UCSC_URL})

################################################################
################################################################
###########  SOFTWARE FOR HIGH-THROUGHPUT SEQUENCING ###########
################################################################
################################################################

################################################################
## TopHat - discovery splice junctions with RNA-seq
install_tophat: _download_tophat _compile_tophat

TOPHAT_BASE_DIR=${SRC_DIR}/TopHat
TOPHAT_VERSION=1.2.0
TOPHAT_ARCHIVE=tophat-${TOPHAT_VERSION}.tar.gz
TOPHAT_URL=http://tophat.cbcb.umd.edu/downloads/${TOPHAT_ARCHIVE}
TOPHAT_DISTRIB_DIR=${TOPHAT_BASE_DIR}/tophat-${TOPHAT_VERSION}
_download_tophat:
	@echo
	@echo "Downloading TopHat"
	@mkdir -p ${TOPHAT_BASE_DIR}
	wget -nd  --directory-prefix ${TOPHAT_BASE_DIR} -rNL ${TOPHAT_URL}
	(cd ${TOPHAT_BASE_DIR}; tar -xpzf ${TOPHAT_ARCHIVE})
	@echo ${TOPHAT_DISTRIB_DIR}


## COMPILATION DOES NOT WORK - TO CHECK
TOPHAT_COMPILE_DIR=`dirname ${BIN_DIR}`
TOPHAT_BIN_DIR=${TOPHAT_COMPILE_DIR}/bin
_compile_tophat:
	@echo
	@echo "Installing TopHat"
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
PEAKSPLITTER_ARCHIVE=PeakSplitter_Cpp.tar.gz
PEAKSPLITTER_URL=http://www.ebi.ac.uk/bertone/software/${PEAKSPLITTER_ARCHIVE}
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
	(cd ${PEAKSPLITTER_DISTRIB_DIR}; ${SUDO} rsync -ruptvl -e ssh PeakSplitter_${OS}/PeakSplitter ${BIN_DIR})

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
	@echo "Linking sissrs in binary dir ${BIN_DIR}"
	(cd ${BIN_DIR}; ln -fs ${SISSRS_BASE_DIR}/sissrs.pl sissrs)


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
BOWTIE_VERSION=2.0.0-beta6
BOWTIE_ARCHIVE=bowtie2-${BOWTIE_VERSION}-${OS}.zip
BOWTIE_URL=http://sourceforge.net/projects/bowtie-bio/files/bowtie2/${BOWTIE_VERSION}/${BOWTIE_ARCHIVE}
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

_compile_bowtie: _compile_bowtie_${OS}

_compile_bowtie_macosx:
	${MAKE} _compile_bowtie_os OS=macos-x86_64

_compile_bowtie_linux:
	${MAKE} _compile_bowtie_os OS=linux-x86_64


_compile_bowtie_os:
	@echo
	@echo "Installing BOWTIE in dir	${BOWTIE_DISTRIB_DIR}"
	(cd ${BOWTIE_BASE_DIR}; unzip ${BOWTIE_ARCHIVE})
	@echo ${BOWTIE_DISTRIB_DIR}
	${SUDO} find  ${BOWTIE_DISTRIB_DIR} -maxdepth 1 -perm 755 -type f  -exec rsync -uptvL {} ${BIN_DIR} \;

################################################################
## Install  Cis-regulatory Element Annotation System  (CEAS)
CEAS_BASE_DIR=${SRC_DIR}/CEAS
CEAS_VERSION=1.0.2
CEAS_ARCHIVE=CEAS-Package-${CEAS_VERSION}.tar.gz
CEAS_URL=http://liulab.dfci.harvard.edu/CEAS/src/${CEAS_ARCHIVE}
CEAS_DISTRIB_DIR=${CEAS_BASE_DIR}/CEAS-Package-${CEAS_VERSION}
install_ceas: _download_ceas _compile_ceas 

_download_ceas:
	@echo
	@echo "Downloading CEAS"
	@mkdir -p ${CEAS_BASE_DIR}
	wget -nd  --directory-prefix ${CEAS_BASE_DIR} -rNL ${CEAS_URL}

CEAS_COMPILE_DIR=`dirname ${BIN_DIR}`
_compile_ceas:
	@echo
	@echo "Installing CEAS in dir	${CEAS_DISTRIB_DIR}"
	(cd ${CEAS_BASE_DIR}; tar -xpzf ${CEAS_ARCHIVE})
	@echo ${CEAS_DISTRIB_DIR}
	(cd  ${CEAS_DISTRIB_DIR}; ${SUDO} python setup.py install --prefix=${CEAS_COMPILE_DIR})
	@echo "CEAS was installed in dir ${BIN_DIR}"
	@echo "Before using CEAS, you need to add a line to the log-in shell script (i.e. .bashrc in case of bash shell)"
	@echo "Adapt the python version in the path below"
	@echo 'export PYTHONPATH=$$PYTHONPATH:${BIN_DIR}/lib/python2.7/site-packages'


################################################################
## Install  SAMTOOLS
SAMTOOLS_BASE_DIR=${SRC_DIR}/samtools
SAMTOOLS_VERSION=0.1.18
SAMTOOLS_ARCHIVE=samtools-${SAMTOOLS_VERSION}.tar.bz2
SAMTOOLS_URL=http://sourceforge.net/projects/samtools/files/samtools/${SAMTOOLS_VERSION}/${SAMTOOLS_ARCHIVE}
SAMTOOLS_DISTRIB_DIR=${SAMTOOLS_BASE_DIR}/samtools-${SAMTOOLS_VERSION}
install_samtools: _download_samtools _compile_samtools 

_download_samtools:
	@echo
	@echo "Downloading SAMTOOLS"
	@mkdir -p ${SAMTOOLS_BASE_DIR}
	wget -nd  --directory-prefix ${SAMTOOLS_BASE_DIR} -rNL ${SAMTOOLS_URL}

_compile_samtools:
	@echo
	@echo "Installing SAMTOOLS in dir	${SAMTOOLS_DISTRIB_DIR}"
	(cd ${SAMTOOLS_BASE_DIR}; tar --bzip2 -xpf ${SAMTOOLS_ARCHIVE})
	@echo ${SAMTOOLS_DISTRIB_DIR}
	(cd ${SAMTOOLS_DISTRIB_DIR}; make)
	${SUDO} find  ${SAMTOOLS_DISTRIB_DIR} -maxdepth 1 -perm 755 -type f  -exec rsync -uptvL {} ${BIN_DIR} \;

################################################################
## Install  SRA toolkit
SRA_OS=mac64
SRA_BASE_DIR=${SRC_DIR}/sra
SRA_VERSION=2.1.10
SRA_ARCHIVE=sratoolkit.${RSA_VERSION}-${SRA_OS}.tar.gz
SRA_URL=ftp-private.ncbi.nlm.nih.gov/sra/sdk/${RSA_VERSION}/${SRA_ARCHIVE}
SRA_DISTRIB_DIR=${SRA_BASE_DIR}/sra-${SRA_VERSION}
install_sra: _download_sra _compile_sra 

_download_sra:
	@echo
	@echo "Downloading SRA"
	@mkdir -p ${SRA_BASE_DIR}
	wget -nd  --directory-prefix ${SRA_BASE_DIR} -rNL ${SRA_URL}

_compile_sra:
	@echo
	@echo "Installing SRA in dir	${SRA_DISTRIB_DIR}"
	(cd ${SRA_BASE_DIR}; tar --bzip2 -xpf ${SRA_ARCHIVE})
	@echo ${SRA_DISTRIB_DIR}
	(cd ${SRA_DISTRIB_DIR}; make)
	${SUDO} find  ${SRA_DISTRIB_DIR} -maxdepth 1 -perm 755 -type f  -exec rsync -uptvL {} ${BIN_DIR} \;

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
	(cd ${SWEMBL_DISTRIB_DIR}; make; ${SUDO} rsync -ruptvl ${SWEMBL_DISTRIB_DIR}/SWEMBL ${BIN_DIR})

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
	@echo 'export PATH=$$PATH:${IGV_DISTRIB_DIR}'



################################################################
## HOMER
HOMER_CONFIG_URL=http://biowhat.ucsd.edu/homer/configureHomer.pl
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
#	@echo 'export PATH=$$PATH:${HOMER_BASE_DIR}/bin'
	@echo 'export HOMER=${HOMER_BASE_DIR}'
	@echo 'export PATH=$$PATH:$$HOMER/bin'
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
	(cd ${CLUSTALW_BASE_DIR}; tar -xvzf ${CLUSTALW_ARCHIVE})
	@echo ${CLUSTALW_BASE_DIR}

_compile_clustalw:
	(cd ${CLUSTALW_BASE_DIR}/${CLUSTALW_SOURCE_DIR}; ./configure; make clean ; make ; \
	${SUDO} rsync -ruptvl src/clustalw2 ${BIN_DIR}/)
	@echo "	clustalw2 should now be executable from ${BIN_DIR}";
