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
	make -f ${RSAT}/makefiles/init_rsat.mk init
	make -f ${RSAT}/makefiles/init_rsat.mk compile_all
	make -f ${RSAT}/makefiles/install_software.mk install_ext_apps

################################################################
## Install Unix packages required for RSAT
UNIX_PACKAGES_COMMON= \
	emacs \
	git \
	cvs \
	wget \
	ghostscript \
	gnuplot \
	graphviz \
	links \
	finger \
	zip \
	unzip \
	python3.2 \
	python3-setuptools \
	python2.7



UNIX_PACKAGES_CENTOS= \
	httpd \
	cpan \
	php \
	glibc.i686 \
	zlib.i686 \
	gd gd gd-devel php-gd perl-GD.x86_64 \
	tetex-latex tetex-doc tetex-fonts

## gfortran required for python scipy
UNIX_PACKAGES_MACOSX= \
	gd \
	gfortran

UNIX_PACKAGES_UBUNTU= \
	make \
	yum \
	apache2 \
	apache2-utils \
	php5 \
	libapache2-mod-php5 \
	php-elisp \
	texlive-latex-base \
	libgd2-xpm-dev \
	libgd-gd2-perl \
	python3 \
	python3-dev

unix_packages_list:
	@echo
	@echo "Required Unix packages"
	@echo "======================"
	@echo ${UNIX_PACKAGES_COMMON} | perl -pe 's|\s+|\n\t|g'
	@echo
	@echo "Additional packages for Ubuntu"
	@echo "=============================="
	@echo ${UNIX_PACKAGES_UBUNTU} | perl -pe 's|\s+|\n\t|g'
	@echo
	@echo "Additional packages for Centos"
	@echo "=============================="
	@echo ${UNIX_PACKAGES_CENTOS} | perl -pe 's|\s+|\n\t|g'
	@echo

PACKAGE_MANAGER_MAC=brew install
PACKAGE_MANAGER_CENTOS=yum install
PACKAGE_MANAGER_UBUNTU=get-apt
UNIX_PACKAGES_CMD=${PACKAGE_MANAGER} ${UNIX_PACKAGES}
unix_packages_cmd:
	@echo "${UNIX_PACKAGES_CMD}"

## Install required Unix packages
unix_packages_install:
	@${UNIX_PACKAGES_CMD}

## Install required Unix packages
unix_packages_install_centos:
	yes | yum upgrade
	@${MAKE} unix_packages_install PACKAGE_MANAGER=${PACKAGE_MANAGER_CENTOS} UNIX_PACKAGES="${UNIX_PACKAGES_COMMON} ${UNIX_PACKAGES_CENTOS}"

## Install required Unix packages
unix_packages_install_ubuntu:
	yes | apt-get upgrade
	@${MAKE} unix_packages_install PACKAGE_MANAGER=${PACKAGE_MANAGER_UBUNTU} UNIX_PACKAGES="${UNIX_PACKAGES_COMMON} ${UNIX_PACKAGES_UBUNTU}"


################################################################
## Install perl modules
## 
## Modules are installed using cpan. Beware, this requires admin
## rights.
PERL_MODULES= \
	YAML \
	CGI \
	MIME::Lite \
	PostScript::Simple \
	Statistics::Distributions \
	Algorithm::Cluster \
	File::Spec \
	POSIX \
	Data::Dumper \
	Digest::MD5::File \
	IO::All \
	LockFile::Simple \
	Object::InsideOut \
	Util::Properties \
	Class::Std::Fast  \
	GD \
	REST::Client \
	JSON \
	MIME::Base64 \
	XML::LibXML \
	XML::LibXML::Simple \
	XML::Compile \
	XML::Compile::Cache \
	XML::Compile::SOAP11 \
	XML::Compile::WSDL11 \
	XML::Parser::Expat \
	XML::Compile::Transport::SOAPHTTP \
	SOAP::WSDL \
	SOAP::Lite \
	SOAP::Transport::HTTP \
	Module::Build::Compat \
	DBI \
	DBD::mysql \
	DB_File \
	LWP::Simple \
	Bio::Perl \
	Bio::Das

## This module is problematic (not maintained anymore), and I am not
## sure it is required anymore. To be checked
PERL_MODULES_EXTRA = SOAP

PERL_MODULES_PROBLEMS= \

PERLMOD_TO_UPGRADE=Archive::Tar

perl_modules_list:
	@echo ${PERL_MODULES} | perl -pe 's|\s+|\n|g'

perl_modules_cmd:
	@echo "${CPAN_CMD} -i ${PERL_MODULES}"

CPAN_OPT=-T 
CPAN_CMD=cpan ${CPAN_OPT}
## Install all Perl modules in one short. Beware: depending on the
## configuration, cpan may ask you to answer y/n for each module and
## dependency.
perl_modules_install:
	@sudo ${CPAN_CMD} -i ${PERL_MODULES}

## This is a somewhat risky but less cumbersome way to install Perl
## modules: automatically send a carriage return to accept the default
## options for all the modules
perl_modules_install_noprompt:
	@yes '' | sudo ${CPAN_CMD} -i ${PERL_MODULES}

perl_modules_install_one_by_one:
	@for module in ${PERL_MODULES} ; do \
		${MAKE} _install_one_perl_module PERL_MODULE=$${module}; \
	done
	${MAKE} perl_modules_install_by_force

## Some Perl modules cannot be installed without force
perl_modules_install_by_force:
	@sudo ${PERL} -MCPAN -e 'force install SOAP::SWDL'
	@sudo ${PERL} -MCPAN -e 'force install SOAP::Transport::HTTP'

## Install a single Perl module
PERL_MODULE=PostScript::Simple
#PERL=`which perl`
PERL='/usr/bin/perl'
_install_one_perl_module:
	@echo "Installing Perl module ${PERL_MODULE}"
	@sudo ${PERL} -MCPAN -e 'install ${PERL_MODULE}'

## Check which modules are installed
perl_modules_check:
	@echo
	@echo "Checking perl modules"
	@echo `hostname` > perl_modules_check.txt
	@for module in ${PERL_MODULES} ; do \
		 perldoc -l $${module} >> check_perl_modules.txt; \
	done
	@echo "	check_perl_modules.txt"

################################################################
## Install modules required for python
PYTHON_MODULES=SUDS Rpy2 lxml SOAPpy
python_modules_list:
	@echo ${PYTHON_MODULES} | perl -pe 's|\s+|\n|g'

python_modules_install:
	@echo
	@echo "Installing modules for python"
	@for module in ${PYTHON_MODULES} ; do \
		${SUDO} easy_install $${module}; \
	PYTHON3

MODULES_done=numpy scipy soappy
python3_modules_install:
	@echo
	@echo "Installing modules for python3"
	@echo "${PYTHON3_MODULES}"
	@echo
	@for module in ${PY	THON3_MODULES} ; do \
		${SUDO} pip install -U $${module}; \
	done

################################################################
## Install R modules required for some RSAT scripts
R_MODULES=RJSONIO reshape plyr dendroextras dendextend
## Note: package  ctc does not exist. To check with Jaime Castro
r_modules_list:
	@echo ${R_MODULES} | perl -pe 's|\s+|\n|g'

r_modules_install_all:
	@echo
	@echo "Installing R modules"
	@for m in ${R_MODULES}; do \
		${MAKE} r_modules_install_one R_MODULE=$${m}; \
	done

R_MODULE=RJSONIO
r_modules_install_one:
	${SUDO} R CMD INSTALL ${R_MODULE}

BIOCONDUCTOR_MODULES=ctc
r_bioconductor_modules:
	for module in ${BIOCONDUCTOR_MODULES}; do \
		echo "Insalling bioconductor module	$${module}"; \
		${SUDO} echo "source('http://bioconductor.org/biocLite.R'); biocLite('"$${module}"')" \
		| R --slave --no-save --no-restore --no-environ ; \
	done

################################################################
## Install tex-live for generating the doc
##
## I found the instructions here:
##  http://tex.stackexchange.com/questions/1092/how-to-install-vanilla-texlive-on-debian-or-ubuntu
TL_VERSION=20131213
install_latex:
	wget http://mirror.ctan.org/systems/texlive/tlnet/install-tl-unx.tar.gz
	tar -xzf install-tl-unx.tar.gz
	cd install-tl-${TL_VERSION}
	sudo ./install-tl


LATEX_PACKAGES=pst-pdf ifplatform 
install_latex_packages:
	sudo tlmgr install ${LATEX_PACKAGES}

# ## Some modules must be upgraded befinre installing required ones
# upgrade_perl_modules:
# 	@for module in ${PERLMOD_TO_UPGRADE}; do \
# 		${MAKE} _upgrade_one_perl_module PERL_MODULE=$${module}; \
# 	done

# ## Upgrade a single Perl module
# _upgrade_one_perl_module:
# 	@echo "Upgrading Perl module ${PERL_MODULE}"
# 	@sudo ${PERL} -MCPAN -e 'upgrade ${PERL_MODULE}'


################################################################
## Install the BioPerl library
##
## For this example, we install Bioperl and EnsEMBL libraries 
## in $RSAT/lib, but you can install it in some other place
### (password is 'cvs')
_old_bioperl:
	@mkdir -p ${RSAT}/lib
	@echo "Password is 'cvs'"
	@cvs -d :pserver:cvs@code.open-bio.org:/home/repository/bioperl login
	(cd ${RSAT}/lib;  cvs -d :pserver:cvs@code.open-bio.org:/home/repository/bioperl checkout bioperl-live)

_old_bioperl_git:
	@echo "This method is obsolete, BioPerl module can now be installed with cpan"
	@mkdir -p $RSAT/lib
	@cd $RSAT/lib
	git clone git://github.com/bioperl/bioperl-live.git

bioperl_install:
	@${MAKE} _install_one_perl_module PERL_MODULE=Bio::Perl

bioperl_test:
	perl -MBio::Perl -le 'print Bio::Perl->VERSION;'


