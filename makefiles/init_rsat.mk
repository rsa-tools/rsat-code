################################################################
## Initialize Regulatory Sequence Analysis Tools (in principle, this
## script should be used only once at installation).

MAKEFILE=${RSAT}/makefiles/init_rsat.mk
MAKE = make -sk -f ${MAKEFILE}
include ${RSAT}/RSAT_config.mk

### tags
usage:
	@echo "usage: make [-OPT='options'] target"
	@echo "implemented targets"
	@perl -ne 'if (/^([a-z]\S+):/){ print "\t$$1\n";  }' ${MAKEFILE}

################################################################
## Initialize directories and config files
SUPPORTED_ORGANISMS=public_html/data/supported_organisms.tab
LOGS_DIR=${RSAT}/public_html/logs
COUNT_FILE=${LOGS_DIR}/count-file

EXEC_FILES='bin/*' \
	'python-scripts/*' \
	'perl-scripts/*' \
	'perl-scripts/parsers/*' \
	'public_html/*.cgi' \
	'public_html/web_services/*.cgi' \
	'ws_clients/perl_clients/*.pl'
init:
	@echo ""
	@echo "Creating directories"
	@echo "	data	${RSAT}/public_html/data"
	@mkdir -p public_html/data
	@echo "Options +Indexes" > public_html/data/.htaccess
	@echo "	genomes	${RSAT}/public_html/data/genomes"
	@mkdir -p public_html/data/genomes
#	mkdir -p public_html/data/KEGG
#	mkdir -p public_html/data/metabolic_networks
#	mkdir -p public_html/data/metabolic_networks/GER_files
	@echo "	bin	${RSAT}/bin"
	@mkdir -p bin
	@echo "	lib	${RSAT}/lib"
	@mkdir -p lib
	@${MAKE} _create_download_dir

	@echo "	tmp	${RSAT}/public_html/tmp"
	@mkdir -p public_html/tmp
	@mkdir -p public_html/tmp/peak-footprints_output; chmod 777 public_html/tmp/peak-footprints_output
	@chmod 777 public_html/tmp
#	echo "Options -Indexes" > public_html/tmp/.htaccess
	@rm -f public_html/tmp/index.html
	@echo "<html><body><b>Forbidden</b></body></html>" > public_html/tmp/index.html
	@chmod 444 public_html/tmp/index.html

	@echo "	logs	${LOGS_DIR}"
	@mkdir -p ${LOGS_DIR}
	@chmod 777 ${LOGS_DIR}
	@mkdir -p ${LOGS_DIR}/peak-footprints_logs; chmod 777 ${LOGS_DIR}/peak-footprints_logs
#	echo "Options -Indexes" > ${LOGS_DIR}/.htaccess
	@rm -f ${LOGS_DIR}/index.html
	@echo "<html><body></b<Forbidden</b></body></html>" > ${LOGS_DIR}/index.html
	@chmod 444 ${LOGS_DIR}/index.html

	@echo
	@echo "Setting exec rights to executable files"
	@for f in ${EXEC_FILES} ; do \
		echo "	$${f}"; \
		chmod 755 $${f}; \
	done ; \

#	@chmod -R 755 bin
#	@chmod 755 python-scripts/*
#	@chmod 755 perl-scripts/*
#	@chmod 755 perl-scripts/parsers/*
#	@chmod 755 public_html/*.cgi
#	@chmod 755 public_html/web_services/*.cgi
#	@chmod 755 ws_clients/perl_clients/*.pl

	@echo ""
	@echo "Creating links"
	@echo "	data"
	@ln -fs public_html/data .
	@echo "	tmp"
	@ln -fs public_html/tmp .
	@echo "	logs"
	@ln -fs ${LOGS_DIR} logs

	@echo ""
	@echo "Checking config files"
	@if [ -f "${SUPPORTED_ORGANISMS}" ] ; then \
		echo "	File already exists	${SUPPORTED_ORGANISMS}" ; \
	else \
		echo "	Creating empty file with supported organisms ${SUPPORTED_ORGANISMS}" ; \
	fi
	@if [ -f "${COUNT_FILE}" ] ; then \
		echo "	File already exists	${COUNT_FILE}" ; \
	else \
		echo "	Creating count file ${COUNT_FILE}" ; \
		echo "0" > ${COUNT_FILE}; \
	fi
	@if [ -f "${RSAT}/RSAT_config.props" ] ; then \
		echo "	RSAT property file already exists	${RSAT}/RSAT_config.props" ; \
	else \
		echo "	Creating RSAT property file ${RSAT}/RSAT_config.props" ; \
		cp ${RSAT}/RSAT_config_default.props ${RSAT}/RSAT_config.props; \
	fi
	@if [ -f "${RSAT}/RSAT_config.mk" ] ; then \
		echo "	RSAT makefiles config already exists	${RSAT}/RSAT_config.mk" ; \
	else \
		echo "	Creating RSAT config for makefiles ${RSAT}/RSAT_config.mk" ; \
		cp ${RSAT}/RSAT_config_default.mk ${RSAT}/RSAT_config.mk; \
	fi
	@chmod a+w ${COUNT_FILE}



################################################################
## Create a directory for downloading genomes
_create_download_dir:
	cd ${RSAT}
	mkdir -p downloads
	(cd downloads; ln -fs $RSAT/makefiles/downloads.mk ./makefile)


## Init Web services
ws_stubb:
	@echo
	@echo "Initiating Web services at ${RSAT_WS}"
	(cd ${RSAT}/ws_clients/perl_clients/; chmod 755 *.pl; make stubb SERVER=${RSAT_WS})

ws_stubb_test:
	@echo
	@echo "Testing Web services at ${RSAT_WS}"
	(cd ${RSAT}/ws_clients/perl_clients/; make test SERVER=${RSAT_WS})


################################################################
## Compile and install C/C++ programs that are part of the RSAT
## distribution (since 2009).

## Compile all programs
compile_all: compile_info_gibbs compile_count_words compile_matrix_scan_quick  compile_compare_matrices_quick compile_pathway_tools

PROGRAM=info-gibbs
SRC_DIR=${RSAT}/contrib/${PROGRAM}
BIN=${RSAT}/bin
## It may be necessary to run the synchronization as super-user (su) with the command sudo.
## For this, type:
##   make -f makefiles/init_rsat.mk compile_all SUDO=sudo BIN=[target_dir]
## For instance, on the lab cluster Brussels I run
##   make -f makefiles/init_rsat.mk compile_all SUDO=sudo BIN=/usr/local/bin
SUDO=
compile_one_program:
	@echo "Compiling ${PROGRAM}"
	(cd ${SRC_DIR}; make all; ${SUDO} rsync -ruptL ${SRC_DIR}/${PROGRAM} ${BIN}/)
#	(cd ${SRC_DIR}; make all; ln -fs ${SRC_DIR}/${PROGRAM} ${RSAT}/bin/${PROGRAM})
	@echo ${BIN}/${PROGRAM}
	@echo ""


## Compile and install info-gibbs (developed by Matthieu Defrance)
compile_info_gibbs:
	@${MAKE} compile_one_program PROGRAM=info-gibbs

## Compile and install count-words (developed by Matthieu Defrance)
compile_count_words:
	@${MAKE} compile_one_program PROGRAM=count-words

compile_matrix_scan_quick:
	@${MAKE} compile_one_program PROGRAM=matrix-scan-quick

compile_compare_matrices_quick:
	@${MAKE} compile_one_program PROGRAM=compare-matrices-quick

## Install external tools useful for RSAT
install_rsat_extra:
	${MAKE} -f ${RSAT}/makefiles/install_software.mk install_seqlogo

################################################################
## Compile the NeAT tools (network analysis + pathway analysis)

install_neat_extra:
	${MAKE} -f ${RSAT}/makefiles/install_software.mk install_mcl
	${MAKE} -f ${RSAT}/makefiles/install_software.mk install_rnsc

compile_pathway_tools: compile_floydwarshall compile_kwalks compile_rea

################################################################
## Compile the program floydwarshall, located in the folder
## contrib/floydwarshall
compile_floydwarshall:
	@echo "Compiling tool floydwarshall required for pathway analysis tools"
	@gcc ${RSAT}/contrib/floydwarshall/floydwarshall.c -o ${BIN}/floydwarshall
	@echo "	${BIN}/floydwarshall"

################################################################
## Compile kwalks (subgraph extraction algorithm developed by Jerome
## Callut and Pierre Dupont).
##
## We also need to set read/write access to all users on the
## kwalks/bin directory because it is used to store temporary files
## during execution. Quick and dirty solution, will need to be revised
compile_kwalks:
	@echo "Compiling kwalks (software developed by Jerome Callut and Pierre Dupont, UCL, Belgium)"
	@(cd ${RSAT}/contrib/kwalks/src; make ; \
		echo "Installing lkwalk executable in bin directory ${BIN}"; \
		cd ../bin; rsync -ruptvl lkwalk ${BIN})
	@echo "Setting read/write access to ${RSAT}/contrib/kwalks for temporary files"
	@chmod 777 ${RSAT}/contrib/kwalks
	@chmod 777 ${RSAT}/contrib/kwalks/bin
	@echo "Executable	${RSAT}/contrib/kwalks/bin/lkwalk"
	@${MAKE} check_kwalks_config

check_kwalks_config:
	@echo
	@echo "Checking KWALKS_ROOT config in RSAT config file"
	@echo "Current value of KWALKS_ROOT in ${RSAT}/RSAT_config.props"
	@grep "KWALKS_ROOT" ${RSAT}/RSAT_config.props
	@echo
	@echo "Please check that the above line equals the following"
	@echo "KWALKS_ROOT=${RSAT}/contrib/kwalks/bin"

check_lkwalk_help:
	@echo "Checkin lkwalk help"
	${RSAT}/contrib/kwalks/bin/lkwalk

################################################################
## Compile REA (shortest path finding algorithm)

compile_rea:
	@echo "Compiling REA"
	@(cd ${RSAT}/contrib/REA/; \
		make; rsync -ruptvl REA ${BIN})
	@echo "Setting read/write access to  ${RSAT}/contrib/REA for temporary files"
	@chmod 777 ${RSAT}/contrib/REA
	@echo "Executable	 ${BIN}/REA"
	@${MAKE} check_rea_config

## Specific configuration to compile REA on Mac OSX
compile_rea_macosx:
	${MAKE} compile_rea CFLAGS='-O3 -Wall -I.'

check_rea_config:
	@echo
	@echo "Checking REA_ROOT config in RSAT config file"
	@echo "Current value of REA_ROOT in ${RSAT}/RSAT_config.props"
	@grep "REA_ROOT" ${RSAT}/RSAT_config.props
	@echo
	@echo "Please check that the above line equals the following"
	@echo "REA_ROOT=${RSAT}/contrib/REA"

check_rea_help:
	@echo "Checking rea help"
	${RSAT}/contrib/REA/REA -help

################################################################
##########   Ubuntu-specific installation
################################################################

################################################################
## install some modules required for proper functioning of RSAT o
## Ubuntu 12
ubuntu_addons:
	sudo apt-get install zlib1g-dev
	sudo apt-get install libgd2-xpm-dev
	sudo apt-get install libxml2-dev
	sudo apt-get install libmysqlclient15-dev
	sudo apt-get install libdb-dev
	sudo apt-get install libberkeleydb-perl
	sudo apt-get install ia32-libs

## Note Berkeley db package should be isntalled with synaptic
ubuntu_install_logwatch:
	sudo apt-get install logwatch
	sudo logwatch.conf

################################################################
## Install tomcat7 for Ubuntu server
ubuntu_tomcat:
	sudo apt-get install libtomcat7-java
	sudo apt-get install tomcat7-common
	sudo apt-get install tomcat7
	sudo apt-get install tomcat7-docs

################################################################
## Restart the pathway tools services
restart_pathway_tools:
	(su - ; su - postgres (pg_ctl -D data -l logfile restart))
	(su -; cd /usr/share/tomcat6/bin; ./catalina.sh start)

################################################################
## Packages specifically required for Ubuntu 14
install_ubuntu14_packages:
	sudo apt-get install python
	sudo apt-get install python-pip
	sudo apt-get install python-dev
	sudo apt-get install ipython
	sudo apt-get install ipython-notebook
	sudo apt-get install git
	sudo apt-get install emacs
	sudo apt-get install make
	sudo apt-get install g++
	sudo apt-get install yum
	sudo apt-get install apache2
	sudo apt-get install php5
	sudo apt-get install libapache2-mod-php5
	sudo apt-get install php-elisp
	sudo apt-get install perldoc
	sudo apt-get install texlive-latex-base
	sudo apt-get install cvs
	sudo apt-get install wget
	sudo apt-get install finger
	sudo apt-get install zip
	sudo apt-get install unzip
	sudo apt-get install libgd2-xpm-dev
	sudo apt-get install libgd-gd2-perl
	sudo apt-get install libxml2-dev
	sudo apt-get install libnet-ssleay-perl
	sudo apt-get install libcrypt-ssleay-perl
	sudo apt-get install libssl-dev
	sudo apt-get install ghostscript
	sudo apt-get install gnuplot
	sudo apt-get install graphviz
	sudo apt-get install links
	sudo apt-get install lib32z1 lib32ncurses5 lib32bz2-1.0

## Some linux packages required for R BioConductor
	sudo apt-get install -y make
	sudo apt-get install libc6-dev
	sudo apt-get install gfortran
	sudo apt-get install build-essential
	sudo apt-get install libreadline-gplv2-dev:i386 lib64readline-gplv2-dev:i386 libreadline-gplv2-dev
	sudo apt-get install libx11-dev
	sudo apt-get install libxt-dev
	sudo apt-get install libcurl4-openssl-dev
	sudo apt-get install libxml2-dev
	sudo apt-get install texlive-full
	sudo apt-get install tcl8.5-dev
	sudo apt-get install tk8.5-dev
	sudo apt-get install libxss-dev
	sudo apt-get install libpng12-dev
	sudo apt-get install libjpeg62-dev
	sudo apt-get install libcairo2-dev

## Unable to locate package. WHy do I have it in the list ? Was this
## package reomved after Ubuntu 12 ?
# 	sudo apt-get install apache2-util

## The following packages "have no candidates". Maybe they were required for Ubuntun 12 ?
##	gfortran-4.3 \
##	gcj
##
##  .. to be completed

################################################################
## Python modules
	sudo apt-get install python python-pip python-dev
	sudo pip install numpy
	sudo pip install scipy
	sudo pip install matplotlib
	sudo pip install soappy

## We need both python2.7 and python3 (for different scripts)
	sudo apt-get install python3 python3-pip 
	sudo pip3 install numpy
	sudo pip3 install scipy
	sudo pip3 install matplotlib
	sudo pip3 install soappy

## This line does not work, although I found it in some docs
## sudo apt-get install pyton3-dev


