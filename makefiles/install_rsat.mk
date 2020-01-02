############################################################
#
# $Id: install_rsat.mk,v 1.76 2013/07/19 06:29:14 jvanheld Exp $
#
# Time-stamp: <2017-03-17 10:39:28 eda>
#
############################################################


################################################################
## This makefile serves to download and compile some C programs developed by
## third parties, and which are required for the Web site (e.g. RNSC, MCL) or
## can optionnally be used in some work flows (e.g. peak-motifs,
## cluster-motifs).
#RSAT=$(CURDIR)

include ${RSAT}/makefiles/util.mk
MAKEFILE=${RSAT}/makefiles/install_rsat.mk
TEST_DATE=`date +%Y-%m-%d`
TEST_DIR=${RSAT}/install_tests_${TEST_DATE}
TEST_CHECK_DIR=${RSAT}/install_tests_correct


## BEWARE: the location of perl and cpan may vary bewteen computers.
## Even worse on Mac OSX there are 2 versions of perl and cpan.
## These two programs are found in both /usr/bin and /usr/local/bin.
## Make sure to use the same path for perl and cpan. 
PERL=`which perl`
CPAN=`which cpan`
#PERL='/usr/bin/perl'

#@echo "perl path	${PERL}"
#@echo "cpan path	${CPAN}"

#################################################################
# Programs used for downloading and sycnrhonizing
SSH=-e 'ssh -x'

################################################################
## Install the RSAT package
install:
	make -f ${RSAT}/makefiles/init_rsat.mk init
	make -f ${RSAT}/makefiles/init_rsat.mk compile_all ## Compile C programs
	make -f ${RSAT}/makefiles/init_rsat.mk ws_init ws_stub  ## Generate the stub for the Web services
	make -f ${RSAT}/makefiles/install_rsat.mk install_r_packages ## Install external and custom R packages
	make -f ${RSAT}/makefiles/install_software.mk install_ext_apps ## Instal external applications required for some RSAT functionalities

## Generate the pdf documentation
## This requires a LaTeX compiler
doc:
	(cd ${RSAT}/do/manuals/; make all_pdf)

################################################################
## Update RSAT from the git repository and check the dependent tasks
## (compilation of C programs, installation of R packages).
update_from_git:
	git pull
	${MAKE} _update_tasks

update_from_wget:
	@echo "DOWNLOAD THE LATEST TAR ARCHIVE: NOT AUTOMATED YET"
	${MAKE} _update_tasks
#	make -f ${RSAT}/makefiles/init_rsat.mk compile_all ## Compile C programs
#	make -f ${RSAT}/makefiles/init_rsat.mk ws_init ws_stub  ## Generate the stub for the Web services
#	make -f ${RSAT}/makefiles/install_rsat.mk install_r_packages ## Install external and custom R packages

_update_tasks:
	make -f ${RSAT}/makefiles/init_rsat.mk compile_all ## Compile C programs
	make -f ${RSAT}/makefiles/init_rsat.mk ws_init ws_stub  ## Generate the stub for the Web services
	make -f ${RSAT}/makefiles/install_rsat.mk install_r_packages ## Install external and custom R packages
	make -f ${RSAT}/makefiles/server.mk denied_ips

################################################################
## Install Unix packages required for RSAT
#
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
	unzip 

# june2018 BCM: these I had to install while installing some R packages
# libcurl-devel mariadb-devel cairo-devel postgresql-devel libXt-devel
UNIX_PACKAGES_CENTOS= \
	httpd \
	cpan \
	php \
	glibc.i686 \
	zlib.i686 \
	libxml2-devel \
	gd gd gd-devel php-gd \
	perl-GD.x86_64 perl-SOAP-WSDL \
	tetex-latex tetex-doc tetex-fonts \
        netcdf-devel \
	python2 \
	python34 \
    python34-setuptools

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
	python2.7 \
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

BREW_PACKAGES=gd python python3 gnuplot ghostscript
macosx_package_list:
	@echo "Packages to install for Mac OSX (brew install)"
	@echo ${BREW_PACKAGES} | perl -pe 's|\s+|\n|g'

macosx_package_install:
	@echo
	@echo "Installing packages for Mac OSX (brew install)"
	@for package in ${BREW_PACKAGES} ; do \
		echo "Installing brew package	$${package}" ; \
		brew install $${package}; \
	done

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
	@${MAKE} unix_packages_install PACKAGE_MANAGER="${PACKAGE_MANAGER_CENTOS}" UNIX_PACKAGES="${UNIX_PACKAGES_COMMON} ${UNIX_PACKAGES_CENTOS}"

## Install required Unix packages
unix_packages_install_ubuntu:
	yes | apt-get upgrade
	@${MAKE} unix_packages_install PACKAGE_MANAGER="${PACKAGE_MANAGER_UBUNTU}" UNIX_PACKAGES="${UNIX_PACKAGES_COMMON} ${UNIX_PACKAGES_UBUNTU}"


################################################################
## Install perl modules
## 
## Modules are installed using cpan. Beware, this requires admin
## rights.
PERL_MODULES= \
	MCE::Shared \
	YAML \
	Module::Build::Compat \
	CGI \
	Email::Sender \
	Email::Sender::Transport::SMTPS \
	Email::Simple \
	Email::Simple::Creator \
	PostScript::Simple	 \
	Statistics::Distributions \
	Math::CDF \
	Algorithm::Cluster \
	File::Spec \
	POSIX \
	Data::Dumper \
	Digest::MD5::File \
	IO::All \
	LockFile::Simple \
	Object::InsideOut Util::Properties \
	Class::Std::Fast  \
	GD \
	DBI \
	DBD::mysql \
	DB_File \
	LWP::Simple \
	REST::Client \
	JSON \
	HTTP::Tiny \
	LWP::UserAgent \
	XML::LibXML \
	XML::LibXML::Simple \
	XML::Parser::Expat \
	XML::Compile \
	XML::Compile::Cache \
	XML::Compile::SOAP11 \
	XML::Compile::WSDL11 \
	XML::Compile::Transport::SOAPHTTP \
	SOAP::Lite \
	SOAP::Packager \
	SOAP::Transport::HTTP \
	SOAP::WSDL \
	Bio::Perl \
	Bio::Das \
	XML::DOM \
	Spreadsheet::WriteExcel::Big \
	Spreadsheet::WriteExcel \
	Log::Log4perl \
	Number::Format \
	OLE::Storage_Lite \
	Template::Plugin::Number::Format \
	Readonly \
	Email::Sender::Transport::SMTPS \
	Parallel::ForkManager \
	MCE::Shared

#t/013_complexType.t ................................... 1/? Can't locate object method "new" via package "MyElement" (perhaps you forgot to load "MyElement"?) at lib/SOAP/WSDL/XSD/Typelib/ComplexType.pm line 213.



## To fix problem with SOAP::WSDL.
## Found at http://www.perlmonks.org/?node_id=823801
## But is apparently not sufficient
#	Module::Build \
#	Devel::Loaded \
#	File::Basename \

## Why was this library required ???

## This module is problematic (not maintained anymore), and I am not
## sure it is required anymore. To be checked
PERL_MODULES_EXTRA = SOAP

PERL_MODULES_PROBLEMS= \

PERLMOD_TO_UPGRADE=Archive::Tar

## Display the parameters for perl module installation
perl_modules_param:
	@echo
	@echo "Parameters for perl module installation"
	@echo "	PERL		${PERL}"
	@echo "	CPAN		${CPAN}"
	@echo "	CPAN_CMD	${CPAN_CMD}"
	@echo "	PERL_MODULES	${PERL_MODULES}"

perl_modules_list:
	@echo ${PERL_MODULES} | perl -pe 's|\s+|\n|g'

perl_modules_cmd:
	@echo "${CPAN_CMD} -i ${PERL_MODULES}"

## Do not test the modules, simply install them
CPAN_OPT=-T 
CPAN_CMD=${CPAN} ${CPAN_OPT}
## Install all Perl modules in one short. Beware: depending on the
## configuration, cpan may ask you to answer y/n for each module and
## dependency.
perl_modules_install:
	@${SUDO} ${CPAN_CMD} -i ${PERL_MODULES}

## This is a somewhat risky but less cumbersome way to install Perl
## modules: automatically send a carriage return to accept the default
## options for all the modules
perl_modules_install_noprompt:
	@yes '' | ${SUDO} ${CPAN_CMD} -i ${PERL_MODULES}

perl_modules_install_one_by_one:
	@for module in ${PERL_MODULES} ; do \
		${MAKE} perl_install_one_module PERL_MODULE=$${module}; \
	done
	${MAKE} perl_modules_install_by_force

## Some Perl modules cannot be installed without force
## About SOAP::Transport::HTTP, I think that there is no doc but the
## module is installed correctly.
PERL_MODULES_TO_FORCE=Object::InsideOut SOAP SOAP::Transport SOAP::WSDL
perl_modules_install_by_force:
	@for module in ${PERL_MODULES_TO_FORCE} ; do \
		${SUDO} ${CPAN_CMD} -f -i $${module}; \
	done

## Install a single Perl module
PERL_MODULE=PostScript::Simple
perl_install_one_module:
	@echo "Installing Perl module ${PERL_MODULE}"
	@${SUDO} ${PERL} -MCPAN -e 'install ${PERL_MODULE}'

## Check which modules are installed
PERL_MODULE_TEST=eval
PERL_MODULES_CHECK_FILE=check_perl_modules_${PERL_MODULE_TEST}.txt
perl_modules_check:
	@echo
	@echo "Checking perl modules ${PERL_MODULE_TEST}"
	@echo "; Checking perl modules ${PERL_MODULE_TEST}" > ${PERL_MODULES_CHECK_FILE}
	@echo "; Host: `hostname`" >> ${PERL_MODULES_CHECK_FILE}
	@echo "; CPAN		${CPAN}"
	@echo "; CPAN_CMD	${CPAN_CMD}"
	@for module in ${PERL_MODULES} ; do \
		 ${MAKE} perl_module_test_${PERL_MODULE_TEST} PERL_MODULE=$${module}; \
	done
	@echo "Report file for Perl modules test"
	@echo "	${PERL_MODULES_CHECK_FILE}"

perl_modules_check_version:
	@${MAKE} perl_modules_check PERL_MODULE_TEST=version

perl_modules_check_doc:
	@${MAKE} perl_modules_check PERL_MODULE_TEST=doc

perl_module_test_eval:
	@echo "	Checking perl module	${PERL_MODULE}	${PERL}"
	@echo "${PERL_MODULE}" | xargs -I MODULE ${PERL} -e  'print eval "use MODULE;1"?"OK\t${PERL_MODULE}\n":"Fail\t${PERL_MODULE}\n"' >> ${PERL_MODULES_CHECK_FILE}

perl_module_test_version:
	@echo "	Checking perl module version	${PERL_MODULE}"
	${PERL} -M${PERL_MODULE} -le 'print ${PERL_MODULE}->VERSION."\t".${PERL_MODULE};' >> ${PERL_MODULES_CHECK_FILE}

perl_module_test_doc:
	@echo "	Checking perl module doc	${PERL_MODULE}"
	perldoc -l ${PERL_MODULE} >> ${PERL_MODULES_CHECK_FILE}

################################################################
## Install modules required for python
PYTHON2_MODULES=numpy \
	scipy \
	SUDS \
	Rpy2 \
	lxml \
	SOAPpy \
	httplib2 \
	requests
python2_modules_list:
	@echo "Modules to install for Python 2"
	@echo ${PYTHON2_MODULES} | perl -pe 's|\s+|\n|g'

python2_modules_install:
	@echo
	@echo "Installing modules for python"
	@for module in ${PYTHON2_MODULES} ; do \
		echo "Installing python2 module	$${module}" ; \
		${SUDO} pip install $${module}; \
	done

PYTHON3_MODULES=numpy scipy snakemake docutils pyyaml
PYTHON3_NONSUPPORTED= soappy
python3_modules_list:
	@echo "Modules to install for Python 3"
	@echo ${PYTHON3_MODULES} | perl -pe 's|\s+|\n|g'

python3_modules_install:
	@echo
	@echo "Installing modules for python3"
	@echo "${PYTHON3_MODULES}"
	@echo
	@for module in ${PYTHON3_MODULES} ; do \
		echo "Installing python3 module	$${module}" ; \
		${SUDO} pip3 install $${module}; \
	done

################################################################
## Install R modules required for some RSAT scripts

## Install R from the command line
## Source: https://www.digitalocean.com/community/tutorials/how-to-set-up-r-on-ubuntu-14-04
install_r:
	@echo
	@echo "Installing R"
	${SUDO} sh -c 'echo "deb http://cran.rstudio.com/bin/linux/ubuntu trusty/" >> /etc/apt/sources.list'
	gpg --keyserver keyserver.ubuntu.com --recv-key E084DAB9
	gpg -a --export E084DAB9 | ${SUDO} apt-key add -
	${SUDO} apt-get update
	${SUDO} apt-get -y install r-base
	${SUDO} apt-get -y install libcurl4-gnutls-dev libxml2-dev libssl-dev
	${SUDO} su - -c "R -e \"install.packages('devtools', repos='http://cran.rstudio.com/')\""
	@echo "R installed"
	@R --version

install_r_packages:
	Rscript ${RSAT}/R-scripts/install_packages_for_rsat.R

R_PACKAGES=RColorBrewer devtools RJSONIO dendextend flux gplots RColorBrewer jpeg
install_r_packages_cmd:
	@for rpack in ${R_PACKAGES}; do \
		echo "Installing R package $${rpack}"; \
		${${SUDO}} R --slave -e "\"if (require('$${rpack}')) {message('$${rpack} already installed')} else {install.packages('$${rpack}', repos='http://cran.rstudio.com/')}\""; \
	done

update_r_packages:
	Rscript ${RSAT}/R-scripts/update_packages_for_rsat.R

##  --slave --no-save --no-restore --no-environ

# R_MODULES=RJSONIO dendextend Rcpp Rclusterpp gplots devtools
# ## reshape plyr: are these still required ?
# ## Note: package  ctc does not exist. To check with Jaime Castro
# r_modules_list:
# 	@echo ${R_MODULES} | perl -pe 's|\s+|\n|g'

# r_modules_install_all:
# 	@echo
# 	@echo "Installing R modules"
# 	@for m in ${R_MODULES}; do \
# 		${MAKE} r_modules_install_one R_MODULE=$${m}; \
# 	done

# R_MODULE=RJSONIO
# r_modules_install_one:
# 	${SUDO} echo "install.packages('${R_MODULE}', repos='http://cran.rstudio.com/', dependencies=TRUE)" | ${SUDO} R --slave --no-save --no-restore --no-environ
# #	${SUDO} R CMD INSTALL ${R_MODULE}

# BIOCONDUCTOR_MODULES=ctc
# r_bioconductor_modules:
# 	for module in ${BIOCONDUCTOR_MODULES}; do \
# 		echo "Installing bioconductor module	$${module}"; \
# 		${SUDO} (echo "source('http://bioconductor.org/biocLite.R'); biocLite('"$${module}"')" \
# 		| R --slave --no-save --no-restore --no-environ ;) \
# 	done

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
	${SUDO} ./install-tl


LATEX_PACKAGES=pst-pdf ifplatform 
install_latex_packages:
	${SUDO} tlmgr install ${LATEX_PACKAGES}


################################################################
## Install the BioPerl library
##
## For this example, we install Bioperl and EnsEMBL libraries 
## in $RSAT/ext_lib, but you can install it in some other place
### (password is 'cvs')
_old_bioperl:
	@mkdir -p ${RSAT}/ext_lib
	@echo "Password is 'cvs'"
	@cvs -d :pserver:cvs@code.open-bio.org:/home/repository/bioperl login
	(cd ${RSAT}/ext_lib;  cvs -d :pserver:cvs@code.open-bio.org:/home/repository/bioperl checkout bioperl-live)

_old_bioperl_git:
	@echo "This method is obsolete, BioPerl module can now be installed with cpan"
	@mkdir -p $RSAT/ext_lib
	@cd $RSAT/ext_lib
	git clone git://github.com/bioperl/bioperl-live.git

bioperl_install:
	@${MAKE} perl_install_one_module PERL_MODULE=Bio::Perl

bioperl_test:
	perl -MBio::Perl -le 'print Bio::Perl->VERSION;'


################################################################
## Installation tests

## Run a series of installation test for the different elements
## required to get a fully working RSAT instance.
all_tests:  test_dir \
	test_random-seq \
	test_purge-sequence

## Create day-specific test directory
test_dir:
	@mkdir -p ${TEST_DIR}
	@echo "Test directory	${TEST_DIR}"

## Check one test by comparing the result to the correct answer
TEST_RESULT=${TEST_DIR}/random-seq_header.txt
TEST_CORRECT=${TEST_CHECK_DIR}/random-seq_header.txt
TEST_NAME=test_random-seq
TEST_DIFF=`diff ${TEST_RESULT} ${TEST_CORRECT} | wc -l | xargs`
test_check:
#	@echo "Test check"
#	@echo "	TEST_NAME	${TEST_NAME}"
#	@echo "	TEST_RESULT	${TEST_RESULT}"
#	@echo "	TEST_CORRECT	${TEST_CORRECT}"
#	@echo "	TEST_DIFF	${TEST_DIFF}"
	if [ "${TEST_DIFF}" -eq "0" ]; \
	then echo "PASSED	${TEST_NAME}	${TEST_RESULT}"; \
	else echo "FAILED	${TEST_NAME}	${TEST_RESULT}	${TEST_CORRECT}"; fi

## A simple Perl script that does not require external revices
test_random-seq:
	@random-seq -l 100 | grep '^>' > ${TEST_DIR}/random-seq_header.txt
	@${MAKE} test_check  TEST_NAME=test_random-seq \
		TEST_RESULT=${TEST_DIR}/random-seq_header.txt \
		TEST_CORRECT=${TEST_CHECK_DIR}/random-seq_header.txt

## Check that the vmatch library is correctly installed
test_purge-sequence:
	@random-seq -l 100 | purge-sequence > ${TEST_DIR}/purge-seq_res.txt 2> ${TEST_DIR}/purge-seq_stderr.txt
	@${MAKE} test_check  TEST_NAME=test_purge-seq \
		TEST_RESULT=${TEST_DIR}/purge-seq_stderr.txt \
		TEST_CORRECT=${TEST_CHECK_DIR}/purge-seq_stderr.txt
## Test retrieve-ensembl-seq, which requires Ensembl API Perl libraries
test_retrieve-ensembl-seq:
	make -f ${RSAT}/makefiles/retrieve-ensembl-seq_demo.mk one_gene \
		> ${TEST_DIR}/retrieve-ensembl-test_res.tab \
		2> ${TEST_DIR}/retrieve-ensembl-test_stderr.txt
	@${MAKE} test_check  TEST_NAME=retrieve-ensembl-seq \
		TEST_RESULT=${TEST_DIR}/retrieve-ensembl-test_res.tab \
		TEST_CORRECT=${TEST_CHECK_DIR}/retrieve-ensembl-test_res.tab

## Test compilation of C programs by running the help message
test_c_compilations:
	info-gibbs -h > ${TEST_DIR}/info-gibbs_help.txt

test_ensembl_available_species:
	@make -f makefiles/install_genomes_from_ensemblgenomes.mk available_species DB=ensembl
	@grep -v '^;' results/ensemblgenomes/available_species_ensembl_release${ENSEMBL_RELEASE}_${TEST_DATE}.txt | wc -l > results/ensemblgenomes/available_species_ensembl_release${ENSEMBL_RELEASE}_${TEST_DATE}_nb_genomes.txt

