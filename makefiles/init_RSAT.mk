################################################################
## Initialize Regulatory Sequence Analysis Tools (in principle, this
## script should be used only once at installation).

MAKEFILE=${RSAT}/makefiles/init_RSAT.mk
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
COUNT_FILE=public_html/logs/count-file
CONFIG_FILE=RSA.config
init:
	@echo ""
	@echo "Creating directories"

	@echo ""
	@echo "Initializing directory	public_html/data"
	mkdir -p public_html/data
	echo "Options +Indexes" > public_html/data/.htaccess
	mkdir -p public_html/data/genomes
#	mkdir -p public_html/data/KEGG
#	mkdir -p public_html/data/metabolic_networks
#	mkdir -p public_html/data/metabolic_networks/GER_files
	mkdir -p bin
	mkdir -p lib
	@${MAKE} _create_download_dir

	@echo ""
	@echo "Initializing directory	public_html/tmp"
	mkdir -p public_html/tmp
	mkdir -p public_html/tmp/peak-footprints_output; chmod 777 public_html/tmp/peak-footprints_output
	chmod 777 public_html/tmp
#	echo "Options -Indexes" > public_html/tmp/.htaccess
	echo "<html><body><b>Forbidden</b></body></html>" > public_html/tmp/index.html
	chmod 644 public_html/tmp/index.html

	@echo ""
	@echo "Initializing directory	public_html/logs"
	chmod 777 public_html/logs
	mkdir -p public_html/logs
	mkdir -p public_html/logs/peak-footprints_logs; chmod 777 public_html/logs/peak-footprints_logs
#	echo "Options -Indexes" > public_html/logs/.htaccess
	echo "<html><body></b<Forbidden</b></body></html>" > public_html/logs/index.html
	chmod 644 public_html/logs/index.html


	@echo
	@echo "Setting exec rights to script directories"
	chmod -R 755 bin
	chmod 755 python-scripts/*
	chmod 755 perl-scripts/*
	chmod 755 perl-scripts/parsers/*
	chmod 755 public_html/*.cgi
	chmod 755 public_html/web_services/*.cgi
	chmod 755 ws_clients/perl_clients/*.pl

	@echo ""
	@echo "Creating links"
	ln -fs public_html/data .
	ln -fs public_html/tmp .
	ln -fs public_html/logs .

	@echo ""
	@echo "Checking config files"
	@if [ -f "${SUPPORTED_ORGANISMS}" ] ; then \
		echo "file already exists	${SUPPORTED_ORGANISMS}" ; \
	else \
		echo "creating empty file with supported organisms ${SUPPORTED_ORGANISMS}" ; \
	fi
	@if [ -f "${COUNT_FILE}" ] ; then \
		echo "file already exists	${COUNT_FILE}" ; \
	else \
		echo "creating count file ${COUNT_FILE}" ; \
		echo "0" > ${COUNT_FILE}; \
	fi
	@if [ -f "${CONFIG_FILE}" ] ; then \
		echo "file already exists	${CONFIG_FILE}" ; \
	else \
		echo "creating config file ${CONFIG_FILE}" ; \
		cp RSA.config.default ${CONFIG_FILE}; \
	fi
	@if [ -f "${RSAT}/RSAT_config.props" ] ; then \
		echo "RSAT property file already exists	${RSAT}/RSAT_config.props" ; \
	else \
		echo "Creating RSAT property file ${RSAT}/RSAT_config.props" ; \
		cp ${RSAT}/RSAT_config_default.props ${RSAT}/RSAT_config.props; \
	fi
	@if [ -f "${RSAT}/RSAT_config.mk" ] ; then \
		echo "RSAT makefiles config already exists	${RSAT}/RSAT_config.mk" ; \
	else \
		echo "Creating RSAT config for makefiles ${RSAT}/RSAT_config.mk" ; \
		cp ${RSAT}/RSAT_config_default.mk ${RSAT}/RSAT_config.mk; \
	fi
	chmod a+w ${COUNT_FILE}



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
##   make -f makefiles/init_RSAT.mk compile_all SUDO=sudo BIN=[target_dir]
## For instance, on the lab cluster Brussels I run
##   make -f makefiles/init_RSAT.mk compile_all SUDO=sudo BIN=/usr/local/bin
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
	@(cd ${RSAT}/contrib/; \
		tar -xpzf kwalks/kwalks.tgz; \
		cd kwalks/src; make ; \
		echo "Installing lkwalk executable in bin directory ${BIN}"; \
		cd ../bin; rsync -ruptvl lkwalk ${BIN})
	@echo "Setting read/write access to ${RSAT}/contrib/kwalks for temporary files"
	@chmod 777 ${RSAT}/contrib/kwalks
#	@chmod 777 ${RSAT}/contrib/kwalks/bin
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
		tar -xpzf REA.tgz; \
		make; rsync -ruptvl REA ${BIN})
	@echo "Setting read/write access to  ${RSAT}/contrib/REA for temporary files"
	@chmod 777 ${RSAT}/contrib/REA
	@echo "Executable	 ${BIN}/rea"
	@${MAKE} check_rea_config

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
## install  some modules required for proper functioning of RSAT o Ubuntu
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
## Packages specifically required for Ubuntu systems
install_ubuntu_prereq:
	sudo apt-get install libgd2-xpm-dev
	sudo apt-get install libxml2-dev
	sudo apt-get install libmysqlclient15-dev
	sudo apt-get install libdb-dev
	sudo apt-get install libberkeleydb-perl

