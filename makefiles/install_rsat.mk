############################################################
#
# $Id: install_rsat.mk,v 1.58 2012/03/03 15:15:26 jvanheld Exp $
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
V=1

## Operating system
## Supported: linux | mac
OS=linux

#################################################################
# Programs used for downloading and sycnrhonizing
WGET = wget -np -rNL 
#MAKE=nice -n 19 make -s -f ${MAKEFILE}
RSYNC_OPT = -ruptvl ${OPT}
SSH=-e 'ssh -x'
RSYNC = rsync ${RSYNC_OPT} ${SSH}
APP_SRC_DIR=${RSAT}/app_sources

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

################################################################
## Install perl modules
## 
## Modules are installed using cpan. Beware, this requires admin
## rights.
PERL_MODULES= \
	PostScript::Simple \
	Statistics::Distributions \
	File::Spec \
	POSIX \
	Data::Dumper \
	Util::Properties \
	Class::Std::Fast \
	XML::LibXML \
	DBD::mysql \
	DBI \
	XML::Compile::Cache \
	XML::Compile::SOAP11 \
	XML::Compile::WSDL11 \
	XML::Compile::Transport::SOAPHTTP \
	SOAP::WSDL \
	SOAP::Lite \
	Module::Build::Compat \
	GD \
	DB_File \
	LWP::Simple \
	Bio::Perl \
	Bio::Das
list_perl_modules:
	@echo
	@echo "Perl modules required for RSAT"
	@echo "------------------------------"
	@echo ${PERL_MODULES} | perl -pe 's|\s+|\n|g'
	@echo

install_perl_modules:
	@for module in ${PERL_MODULES} ; do \
		${MAKE} install_one_perl_module PERL_MODULE=$${module}; \
	done

## Install a single Perl module
PERL_MODULE=PostScript::Simple
#PERL=`which perl`
PERL='/usr/bin/perl'
SUDO=sudo
install_one_perl_module:
	@echo "Installing Perl module ${PERL_MODULE}"
	@${SUDO} ${PERL} -MCPAN -e 'install ${PERL_MODULE}'


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
## Install the BioPerl library
## For this example, we install Bioperl and EnsEMBL libraries 
## in $RSAT/lib, but you can install it in some other place
### (password is 'cvs')
_old_bioperl:
	@mkdir -p ${RSAT}/lib
	@echo "Password is 'cvs'"
	@cvs -d :pserver:cvs@code.open-bio.org:/home/repository/bioperl login
	(cd ${RSAT}/lib;  cvs -d :pserver:cvs@code.open-bio.org:/home/repository/bioperl checkout bioperl-live)

bioperl_git:
	@echo "This method is obsolete, BioPerl module can now be installed with cpan"
	@mkdir -p $RSAT/lib
	@cd $RSAT/lib
	git clone git://github.com/bioperl/bioperl-live.git

bioperl_test:
	perl -MBio::Perl -le 'print Bio::Perl->VERSION;'


################################################################
## Install the EnsEMBL Perl API
ENSEMBL_BRANCH=56
ensembl_api:	
	@echo  "Password is 'CVSUSER'"
	@cvs -d :pserver:cvsuser@cvs.sanger.ac.uk:/cvsroot/ensembl login
	@(cd ${RSAT}/lib; \
		cvs -d :pserver:cvsuser@cvs.sanger.ac.uk:/cvsroot/ensembl \
		checkout -r branch-ensembl-${ENSEMBL_BRANCH} ensembl ; \
		cvs -d :pserver:cvsuser@cvs.sanger.ac.uk:/cvsroot/ensembl \
		checkout -r branch-ensembl-${ENSEMBL_BRANCH} ensembl-compara)
	@echo "Don't forget to adapt the following lines in the file ${RSAT}/RSAT_config.props"
	@echo "ensembl=${RSAT}/lib/ensembl/modules"
	@echo "compara=${RSAT}/lib/ensembl-compara/modules"
	@echo "bioperl=${RSAT}/lib/bioperl-live"


################################################################
## Get and install the program seqlogo
SEQLOGO_URL=http://weblogo.berkeley.edu/release
SEQLOGO_TAR=weblogo.2.8.2.tar.gz
SEQLOGO_DIR=${RSAT}/ext/seqlogo
download_seqlogo:
	@mkdir -p ${SEQLOGO_DIR}
	@echo "Getting seqlogo using wget"
	(cd ${SEQLOGO_DIR}; wget -nv -nd ${SEQLOGO_URL}/${SEQLOGO_TAR}; tar -xpzf ${SEQLOGO_TAR})
	@echo "seqlogo dir	${SEQLOGO_DIR}"

install_seqlogo:
	@echo "Installing seqlogo"
	@rsync -ruptl ${SEQLOGO_DIR}/weblogo/seqlogo ${RSAT}/bin/
	@rsync -ruptl ${SEQLOGO_DIR}/weblogo/template.* ${RSAT}/bin/
	@rsync -ruptl ${SEQLOGO_DIR}/weblogo/logo.pm ${RSAT}/bin/

################################################################
## Get and install the program ghostscript
## Note: for Mac users, please go to the ghostscript Web site
GS_URL=http://ghostscript.com/releases/
GS_VER=ghostscript-8.64
GS_TAR=${GS_TAR}.tar.gz
GS_DIR=${RSAT}/ext/ghostscript
download_gs:
	@mkdir -p ${GS_DIR}
	@echo "Getting gs using wget"
	(cd ${GS_DIR}; wget -nv -nd ${GS_URL}/${GS_TAR}; tar -xpzf ${GS_TAR})
	@echo "gs dir	${GS_DIR}"

install_gs:
	@echo "Installing gs"
	(cd ${GS_DIR}; wget -nv -nd ${GS_URL}/${GS_TAR}; tar -xpzf ${GS_TAR})
	(cd ${GS_DIR}/${GS_VER}; ./configure && make)

################################################################
## Get and install the program gnuplot
GNUPLOT_URL=http://sourceforge.net/projects/gnuplot/files/gnuplot/4.2.6/
GNUPLOT_VER=gnuplot-4.2.5
GNUPLOT_TAR=${GNUPLOT_VER}.tar.gz
GNUPLOT_DIR=${RSAT}/ext/gnuplot
download_gnuplot:
	@mkdir -p ${GNUPLOT_DIR}
	@echo "Getting gnuplot using wget"
	(cd ${GNUPLOT_DIR}; wget -nv -nd ${GNUPLOT_URL}/${GNUPLOT_TAR}; tar -xpzf ${GNUPLOT_TAR})
	@echo "gnuplot dir	${GNUPLOT_DIR}"

install_gnuplot:
	@echo "Installing gnuplot"
	(cd ${GNUPLOT_DIR}/${GNUPLOT_VER}; ./configure && make)

################################################################
## Install BEDTools
##
## BEDTools is a collection of utilities for comparing, summarizing,
## and intersecting genomic features in BED, GTF/GFF, VCF and BAM
## formats.
#bedtools: git_bedtools compile_bedtools install_bedtools
bedtools: download_bedtools compile_bedtools install_bedtools

BED_VERSION=2.13.3
BED_ARCHIVE=BEDTools.v${BED_VERSION}.tar.gz
BED_URL=http://bedtools.googlecode.com/files/${BED_ARCHIVE}
BED_BASE_DIR=${APP_SRC_DIR}/BEDTools
BED_DOWNLOAD_DIR=${BED_BASE_DIR}/BEDTools-Version-${BED_VERSION}
download_bedtools:
	@echo
	@echo "Downloading BEDTools ${BED_VERSION}"
	@echo
	@mkdir -p ${BED_BASE_DIR}
	(cd ${BED_BASE_DIR}; wget -nv -nd ${BED_URL} ; tar -xpzf ${BED_ARCHIVE})
	@echo ${BED_DOWNLOAD_DIR}

BED_GIT_DIR=${APP_SRC_DIR}/bedtools
git_bedtools:
	@mkdir -p ${BED_GIT_DIR}
	(cd ${APP_SRC_DIR}; git clone git://github.com/arq5x/bedtools.git)

#BED_SRC_DIR=${BED_GIT_DIR}
BED_SRC_DIR=${BED_DOWNLOAD_DIR}
BED_BIN_DIR=${BED_SRC_DIR}/bin
compile_bedtools:
	@echo
	@echo "Installing bedtools from ${BED_SRC_DIR}"
	@echo
	@mkdir -p ${BED_SRC_DIR}
	(cd ${BED_SRC_DIR}; make clean; make all) #; make test; make install)

BIN=${RSAT}/bin
install_bedtools:
	@echo
	@echo "Synchronizing bedtools from ${BEN_BIN_DIR} to ${BIN}"
	@echo
	@mkdir -p ${BIN}
	rsync -ruptvl ${BED_BIN_DIR}/* ${BIN}

################################################################
## Install MEME (Tim Bailey)
MEME_BASE_DIR=${APP_SRC_DIR}/MEME
MEME_VERSION=4.7.0
#MEME_VERSION=current
MEME_ARCHIVE=meme_${MEME_VERSION}.tar.gz
MEME_URL=http://meme.nbcr.net/downloads/${MEME_ARCHIVE}
MEME_DISTRIB_DIR=${MEME_BASE_DIR}/meme_${MEME_VERSION}
download_meme:
	@echo
	@echo "Downloading MEME ${MEME_VERSION}"
	@echo
	@mkdir -p ${MEME_BASE_DIR}
	(cd ${MEME_BASE_DIR}; wget -nv -nd ${MEME_URL} ; tar -xpzf ${MEME_ARCHIVE})
	@echo ${MEME_DISTRIB_DIR}

MEME_INSTALL_DIR=${MEME_DISTRIB_DIR}_installed
MEME_BIN_DIR=${MEME_INSTALL_DIR}/bin
install_meme:
	@echo
	@echo "Installing MEME ${MEME_VERSION}"
	@echo
	@mkdir -p ${MEME_INSTALL_DIR}
	(cd ${MEME_DISTRIB_DIR}; ./configure --prefix=${MEME_INSTALL_DIR} --with-url="http://localhost/meme")
	(cd ${MEME_DISTRIB_DIR}; make clean; make ; make test; make install)
	@echo "Please edit the RSAT configuration file"
	@echo "	${RSAT}/RSAT_config.props"
	@echo "and copy-paste the following lines to specify the MEME bin pathway"
	@echo "	meme_dir=${MEME_BIN_DIR}"
	@echo "	MEME_DIRECTORY=${MEME_BIN_DIR}"
	@echo "This will allow RSAT programs to idenfity meme path on this server."
	@echo
	@echo "You can also add the MEME bin directory in your path."
	@echo "If your shell is bash"
	@echo "	export PATH=${MEME_BIN_DIR}:\$$PATH"
	@echo "If your shell is csh or tcsh"
	@echo "	setenv PATH ${MEME_BIN_DIR}:\$$PATH"

################################################################
## Install a clustering algorithm "cluster"

################################################################
## Install the graph-based clustering algorithm MCL
CLUSTER_BASE_DIR=${APP_SRC_DIR}/cluster
CLUSTER_VERSION=1.50
CLUSTER_ARCHIVE=cluster-${CLUSTER_VERSION}.tar.gz
CLUSTER_URL= http://bonsai.hgc.jp/~mdehoon/software/cluster/${CLUSTER_ARCHIVE}
CLUSTER_DISTRIB_DIR=${CLUSTER_BASE_DIR}/cluster-${CLUSTER_VERSION}
download_cluster:
	@echo
	@echo "Downloading CLUSTER"
	@mkdir -p ${CLUSTER_BASE_DIR}
	wget -nd  --directory-prefix ${CLUSTER_BASE_DIR} -rNL ${CLUSTER_URL}
	(cd ${CLUSTER_BASE_DIR}; tar -xpzf ${CLUSTER_ARCHIVE})
	@echo ${CLUSTER_DISTRIB_DIR}
#	tar xvfz cluster-1.50.tar.gz
#	(cd cluster-1.50; ./configure --without-x; make clean; make )

CLUSTER_INSTALL_DIR=${RSAT}
CLUSTER_BIN_DIR=${CLUSTER_INSTALL_DIR}/bin
install_cluster:
	@echo
	@echo "Installing CLUSTER"
	@mkdir -p ${CLUSTER_INSTALL_DIR}
	(cd ${CLUSTER_DISTRIB_DIR}; ./configure --without-x --prefix=${CLUSTER_INSTALL_DIR} ; \
	make clean; make ; make install)



################################################################
## Install the graph-based clustering algorithm MCL
MCL_BASE_DIR=${APP_SRC_DIR}/mcl
MCL_VERSION=09-308
MCL_ARCHIVE=mcl-${MCL_VERSION}.tar.gz
MCL_URL=http://www.micans.org/mcl/src/${MCL_ARCHIVE}
MCL_DISTRIB_DIR=${MCL_BASE_DIR}/mcl-${MCL_VERSION}
download_mcl:
	@echo
	@echo "Downloading MCL"
	@mkdir -p ${MCL_BASE_DIR}
	wget -nd  --directory-prefix ${MCL_BASE_DIR} -rNL ${MCL_URL}
	(cd ${MCL_BASE_DIR}; tar -xpzf ${MCL_ARCHIVE})
	@echo ${MCL_DISTRIB_DIR}

MCL_INSTALL_DIR=${RSAT}
MCL_BIN_DIR=${MCL_INSTALL_DIR}/bin
install_mcl:
	@echo
	@echo "Installing MCL"
	@mkdir -p ${MCL_INSTALL_DIR}
	(cd ${MCL_DISTRIB_DIR}; ./configure --prefix=${MCL_INSTALL_DIR} ; \
	make clean; make ; make install)
	@echo "Please edit the RSAT configuration file"
	@echo "	${RSAT}/RSAT_config.props"
	@echo "and copy-paste the following line to specify the MCL bin pathway"
	@echo "	mcl_dir=${MCL_BIN_DIR}"
	@echo "This will allow RSAT programs to idenfity mcl path on this server."
	@echo
	@echo "You can also add the MCL bin directory in your path."
	@echo "If your shell is bash"
	@echo "	export PATH=${MCL_BIN_DIR}:\$$PATH"
	@echo "If your shell is csh or tcsh"
	@echo "	setenv PATH ${MCL_BIN_DIR}:\$$PATH"

################################################################
## Install the graph-based clustering algorithm RNSC
RNSC_BASE_DIR=${APP_SRC_DIR}/rnsc
RNSC_VERSION=09-308
RNSC_ARCHIVE=rnsc.zip
RNSC_URL=http://www.cs.utoronto.ca/~juris/data/rnsc/rnsc.zip
download_rnsc:
	@echo
	@echo "Downloading RNSC"
	@mkdir -p ${RNSC_BASE_DIR}
	wget --no-directories  --directory-prefix ${RNSC_BASE_DIR} -rNL ${RNSC_URL}
	(cd ${RNSC_BASE_DIR}; unzip ${RNSC_ARCHIVE})
	@echo ${RNSC_BASE_DIR}

RNSC_BIN_DIR=${RSAT}/bin
install_rnsc:
	@echo
	@echo "Installing RNSC"
	@mkdir -p ${RNSC_BIN_DIR}
	(cd ${RNSC_BASE_DIR}; make ; \
	mv -f rnsc ${RNSC_BIN_DIR}; \
	mv -f rnscfilter ${RNSC_BIN_DIR}; \
	mv -f rnscconvert ${RNSC_BIN_DIR}; \
	)
	@echo "Please edit the RSAT configuration file"
	@echo "	${RSAT}/RSAT_config.props"
	@echo "and copy-paste the following line to specify the RNSC bin pathway"
	@echo "	rnsc_dir=${RNSC_BIN_DIR}"
	@echo "This will allow RSAT programs to idenfity rnsc path on this server."
	@echo
	@echo "You can also add the RNSC bin directory in your path."
	@echo "If your shell is bash"
	@echo "	export PATH=${RNSC_BIN_DIR}:\$$PATH"
	@echo "If your shell is csh or tcsh"
	@echo "	setenv PATH ${RNSC_BIN_DIR}:\$$PATH"


################################################################
## Install BLAST
download_blast: download_blast_${OS}

install_blast: install_blast_${OS}


################################################################
## Install the BLAST on linux
ARCHITECTURE=32
BLAST_BASE_DIR=${APP_SRC_DIR}/blast
BLAST_LINUX_ARCHIVE=blast-*-ia${ARCHITECTURE}-linux.tar.gz
BLAST_URL=ftp://ftp.ncbi.nih.gov/blast/executables/release/LATEST/
BLAST_BIN_DIR=${RSAT}/bin
BLAST_SOURCE_DIR=blast_latest
download_blast_linux:
	@mkdir -p ${BLAST_BASE_DIR}
	wget --no-directories  --directory-prefix ${BLAST_BASE_DIR} -rNL ${BLAST_URL} -A "${BLAST_LINUX_ARCHIVE}"
	(cd ${BLAST_BASE_DIR}; tar -xvzf ${BLAST_LINUX_ARCHIVE}; rm -r ${BLAST_SOURCE_DIR};mv ${BLAST_LINUX_ARCHIVE} ..;mv blast-*  ${BLAST_SOURCE_DIR})
	@echo ${BLAST_BASE_DIR}

install_blast_linux:
	@mkdir -p ${BLAST_BIN_DIR}
	(cd ${BLAST_BIN_DIR}; \
	ln -fs ${BLAST_BASE_DIR}/${BLAST_SOURCE_DIR}/bin/blastall blastall; \
	ln -fs ${BLAST_BASE_DIR}/${BLAST_SOURCE_DIR}/bin/formatdb formatdb; \
	)
	@echo "Please edit the RSAT configuration file"
	@echo "	${RSAT}/RSAT_config.props"
	@echo "and copy-paste the following line to specify the BLAST bin pathway"
	@echo "	blast_dir=${BLAST_BIN_DIR}"
	@echo "This will allow RSAT programs to idenfity BLAST path on this server."
	@echo
	@echo "You can also add the BLAST bin directory in your path."
	@echo "If your shell is bash"
	@echo "	export PATH=${BLAST_BIN_DIR}:\$$PATH"
	@echo "If your shell is csh or tcsh"
	@echo "	setenv PATH ${BLAST_BIN_DIR}:\$$PATH"

################################################################
## Install the BLAST on MAC
BLAST_BASE_DIR=${APP_SRC_DIR}/blast
BLAST_MAC_ARCHIVE=blast-*-universal-macosx.tar.gz
BLAST_URL=ftp://ftp.ncbi.nih.gov/blast/executables/release/LATEST/
BLAST_BIN_DIR=${RSAT}/bin
BLAST_SOURCE_DIR=blast_latest
download_blast_mac:
	@mkdir -p ${BLAST_BASE_DIR}
	wget --no-directories  --directory-prefix ${BLAST_BASE_DIR} -rNL ${BLAST_URL} -A "${BLAST_MAC_ARCHIVE}"
	(cd ${BLAST_BASE_DIR}; tar -xvzf ${BLAST_MAC_ARCHIVE}; rm -r ${BLAST_SOURCE_DIR};mv ${BLAST_MAC_ARCHIVE} ..;mv blast-*  ${BLAST_SOURCE_DIR})
	@echo ${BLAST_BASE_DIR}

install_blast_mac:
	@mkdir -p ${BLAST_BIN_DIR}
	(cd ${BLAST_BIN_DIR}; \
	ln -fs ${BLAST_BASE_DIR}/${BLAST_SOURCE_DIR}/bin/blastall blastall; \
	ln -fs ${BLAST_BASE_DIR}/${BLAST_SOURCE_DIR}/bin/formatdb formatdb; \
	)
	@echo "Please edit the RSAT configuration file"
	@echo "	${RSAT}/RSAT_config.props"
	@echo "and copy-paste the following line to specify the BLAST bin pathway"
	@echo "	blast_dir=${BLAST_BIN_DIR}"
	@echo "This will allow RSAT programs to idenfity BLAST path on this server."
	@echo
	@echo "You can also add the BLAST bin directory in your path."
	@echo "If your shell is bash"
	@echo "	export PATH=${BLAST_BIN_DIR}:\$$PATH"
	@echo "If your shell is csh or tcsh"
	@echo "	setenv PATH ${BLAST_BIN_DIR}:\$$PATH"



################################################################
## Generic call for installing a program. This tag is called with
## specific parameters for each program (consensus, patser, ...)
APP_DIR=${RSAT}/applications
PROGRAM=consensus
PROGRAM_DIR=${APP_DIR}/${PROGRAM}
PROGRAM_ARCHIVE=`ls -1t ${APP_SRC_DIR}/${PROGRAM}* | head -1`
uncompress_program:
	@echo installing ${PROGRAM_ARCHIVE} in dir ${PROGRAM_DIR}
	@mkdir -p ${PROGRAM_DIR}
	@(cd ${PROGRAM_DIR} ;				\
	gunzip -c ${PROGRAM_ARCHIVE} | tar -xf - )

################################################################
## Common frame for installing programs
INSTALLED_PROGRAM=`ls -1t ${APP_DIR}/${PROGRAM}/${PROGRAM}*`
install_program:
	(cd ${PROGRAM_DIR}; make ${INSTALL_OPT})
	(cd bin; ln -fs ${INSTALLED_PROGRAM} ./${PROGRAM})


################################################################
## Install Andrew Neuwald's gibbs sampler (1995 version)
GIBBS_DIR=${APP_DIR}/gibbs/gibbs9_95
install_gibbs:
	${MAKE} uncompress_program PROGRAM=gibbs
	(cd ${GIBBS_DIR}; ./compile; cd ${GIBBS_DIR}/code; make clean)
	(cd bin; ln -fs ${GIBBS_DIR}/gibbs ./gibbs)

################################################################
## Get and install patser (matrix-based pattern matching)
#PATSER_TAR=patser-v3e.1.tar.gz
PATSER_VERSION=patser-v3b.5
#PATSER_VERSION=patser-v3b.5
PATSER_TAR=${PATSER_VERSION}.tar.gz
PATSER_URL=ftp://www.genetics.wustl.edu/pub/stormo/Consensus
PATSER_DIR=${RSAT}/ext/patser/${PATSER_VERSION}
PATSER_APP=`cd ${PATSER_DIR} ; ls -1tr patser-v* | grep -v .tar | tail -1 | xargs`
download_patser:
	@mkdir -p ${PATSER_DIR}
	@echo "Getting patser using wget"
	wget --no-directories  --directory-prefix ${PATSER_DIR} -rNL ${PATSER_URL}/${PATSER_TAR}
	(cd ${PATSER_DIR}; tar -xpzf ${PATSER_TAR})
#	(cd ${PATSER_DIR}; wget -nv  ${PATSER_URL}/${PATSER_TAR}; tar -xpzf ${PATSER_TAR})
	@echo "patser dir	${PATSER_DIR}"

install_patser:
	@echo "Installing patser"
	(cd ${PATSER_DIR}; rm *.o; make)
	rsync -ruptvl ${PATSER_DIR}/${PATSER_APP} ${RSAT}/bin/
#	(cd ${RSAT}/bin; ln -fs ${PATSER_APP} patser)
	@echo "ls -ltr ${RSAT}/bin/patser*"
#	${MAKE} uncompress_program PROGRAM_DIR=${PATSER_DIR} PROGRAM=patser
#	${MAKE} install_program PROGRAM=patser


################################################################
## Install consensus (J.Hertz)
CONSENSUS_VERSION=consensus-v6c.1
CONSENSUS_TAR=${CONSENSUS_VERSION}.tar.gz
CONSENSUS_URL=ftp://www.genetics.wustl.edu/pub/stormo/Consensus
CONSENSUS_DIR=ext/consensus/${CONSENSUS_VERSION}
download_consensus:
	@echo
	@echo "Downloading ${CONSENSUS_VERSION}"
	@echo
	@mkdir -p ${CONSENSUS_DIR}
	(cd ${CONSENSUS_DIR}; wget -nv -nd ${CONSENSUS_URL}/${CONSENSUS_TAR}; tar -xpzf ${CONSENSUS_TAR})
	@echo "consensus dir	${CONSENSUS_DIR}"

install_consensus:
#	${MAKE} uncompress_program PROGRAM=consensus
	${MAKE} install_program PROGRAM=consensus INSTALL_OPT='CPPFLAGS=""'

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
################################################################
###########  SOFTWARE FOR HIGH-THROUGHPUT SEQUENCING ###########
################################################################
################################################################

## TopHat
tophat: download_tophat install_tophat


TOPHAT_BASE_DIR=${APP_SRC_DIR}/TopHat
TOPHAT_VERSION=1.2.0
TOPHAT_ARCHIVE=tophat-${TOPHAT_VERSION}.tar.gz
TOPHAT_URL=http://tophat.cbcb.umd.edu/downloads/${TOPHAT_ARCHIVE}
TOPHAT_DISTRIB_DIR=${TOPHAT_BASE_DIR}/tophat-${TOPHAT_VERSION}
download_tophat:
	@echo
	@echo "Downloading TopHat"
	@mkdir -p ${TOPHAT_BASE_DIR}
	wget -nd  --directory-prefix ${TOPHAT_BASE_DIR} -rNL ${TOPHAT_URL}
	(cd ${TOPHAT_BASE_DIR}; tar -xpzf ${TOPHAT_ARCHIVE})
	@echo ${TOPHAT_DISTRIB_DIR}

TOPHAT_INSTALL_DIR=${RSAT}
TOPHAT_BIN_DIR=${TOPHAT_INSTALL_DIR}/bin
install_tophat:
	@echo
	@echo "Installing TopHat"
	@mkdir -p ${TOPHAT_INSTALL_DIR}
	(cd ${TOPHAT_DISTRIB_DIR}; ./configure --prefix=${TOPHAT_INSTALL_DIR} ; \
	make clean; make ; make install)
