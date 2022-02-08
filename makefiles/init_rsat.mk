################################################################
## Initialize Regulatory Sequence Analysis Tools (in principle, this
## script should be used only once at installation).

RSAT=$(CURDIR)
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
LOG_DIR=${RSAT}/logs
COUNT_FILE=${LOG_DIR}/count-file

#WS_PERL_CLIENT_DIR=public_html/web_services/clients/perl
#WS_PERL_CLIENT_DIR=ws_clients/perl_clients
WS_PERL_CLIENT_DIR=public_html/web_services/clients/perl_clients
WS_STUB_CLIENT_DIR=${WS_PERL_CLIENT_DIR}/soap-wsdl
EXEC_FILES='bin/*' \
	'python-scripts/*' \
	'perl-scripts/*' \
	'perl-scripts/parsers/*' \
	'public_html/*.cgi' \
	'public_html/web_services/*.cgi' \
	'${WS_PERL_CLIENT_DIR}/*.pl' \
	'${WS_STUB_CLIENT_DIR}/*.pl'
init:
	@echo ""
	@echo "Creating directories"
	@echo "	data	${RSAT}/public_html/data"
	@mkdir -p ${RSAT}/bin
	@mkdir -p ${RSAT}/public_html/data
	@echo "Options +Indexes" > public_html/data/.htaccess
	@echo "	genomes	${RSAT}/public_html/data/genomes"
	@mkdir -p ${RSAT}/public_html/data/genomes
	@mkdir -p ${RSAT}/public_html/data/KEGG
	@mkdir -p ${RSAT}/public_html/data/metabolic_networks
	@mkdir -p ${RSAT}/public_html/data/metabolic_networks/GER_files
	@echo "	bin	${RSAT}/bin"
	@mkdir -p ${RSAT}/bin
	@echo "	ext_lib	${RSAT}/ext_lib"
	@mkdir -p ${RSAT}/ext_lib
	@${MAKE} _create_download_dir
	@echo "	tmp	${RSAT}/public_html/tmp"
	@mkdir -p ${RSAT}/public_html/tmp
	@chmod 777 ${RSAT}/public_html/tmp
	@echo "	serialized_genomes	${RSAT}/public_html/tmp/serialized_genomes"
	@mkdir -p ${RSAT}/public_html/tmp/serialized_genomes
	@chmod 777 ${RSAT}/public_html/tmp/serialized_genomes
#	echo "Options -Indexes" > public_html/tmp/.htaccess
	@echo "	logs	${LOG_DIR}"
	@mkdir -p ${LOG_DIR}
	@chmod 777 ${LOG_DIR}
	@${MAKE} init_robots ROBOTS=${RSAT}/public_html/data/robots.txt
	@${MAKE} init_robots ROBOTS=${RSAT}/public_html/tmp/robots.txt
	@${MAKE} init_robots ROBOTS=${LOG_DIR}/robots.txt
	@rm -f public_html/tmp/index.html
	@echo "<html><body><b>Forbidden</b></body></html>" > public_html/tmp/index.html
	@chmod 444 public_html/tmp/index.html
#	@echo "	peak-footprints_logs	${RSAT}/${LOG_DIR}/peak-footprints_logs"
#	@mkdir -p ${LOG_DIR}/peak-footprints_logs; chmod 777 ${LOG_DIR}/peak-footprints_logs
#	echo "Options -Indexes" > ${LOG_DIR}/.htaccess
	@rm -f ${LOG_DIR}/index.html
	@echo "<html><body></b<Forbidden</b></body></html>" > ${LOG_DIR}/index.html
	@chmod 444 ${LOG_DIR}/index.html

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
#	@chmod 755 ${WS_PERL_CLIENT_DIR}/*.pl

	@echo ""
	@echo "Creating links"
	@echo "	data"
	@ln -fs public_html/data .
	@echo "	tmp"
	@ln -fs public_html/tmp .
#	@echo "	logs"
#	@ln -fs ${LOG_DIR} logs

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
	mkdir -p ${RSAT}/downloads
	(cd downloads; ln -fs ${RSAT}/makefiles/downloads.mk ./makefile)

################################################################
## Configuration for robots.
## Details can be found here: http://www.robotstxt.org/orig.html
ROBOTS=${RSAT}/public_html/robots.txt
init_robots:
	@echo
	@echo "Disabling robots from indexing data directory"
	@echo "	Robot file	${ROBOTS}"
	@echo "User-agent: *" > ${ROBOTS}
#	@echo "Disallow: ." >> ${ROBOTS}
#	@echo "Disallow: /rsat/data/" >> ${ROBOTS}
#	@echo "Disallow: /data/" >> ${ROBOTS}
#	@echo "Disallow: data/" >> ${ROBOTS}
#	@echo "Disallow: /rsat/tmp/" >> ${ROBOTS}
#	@echo "Disallow: /tmp/" >> ${ROBOTS}
#	@echo "Disallow: tmp/" >> ${ROBOTS}
#	@echo "	ROBOTS	${ROBOTS}"

## Adapt the IP address in the RSATWS.wsdl file
ws_init:
	@echo
	@echo "Initializing RSATWS.wsdl file with SERVER=RSAT_WS=${RSAT_WS}"
	@perl -pe 's|\[RSAT_WS\]|${RSAT_WS}|g' \
		${RSAT}/public_html/web_services/RSATWS_default.wsdl \
		> ${RSAT}/public_html/web_services/RSATWS.wsdl
	@cp ${RSAT}/public_html/web_services/RSATWS.wsdl ${RSAT}/public_html/web_services/RSATWS_documentation.xml
	@${MAKE} ws_param

## List the parameters for the WS stub
ws_param:
	@echo ""
	@echo "Web services parameters"
	@echo "	RSAT_WS		${RSAT_WS}"
	@echo "Local WSDL file"
	@echo "	${RSAT}/public_html/web_services/RSATWS.wsdl"
	@echo "WSDL access"
	@echo "	${RSAT_WS}/web_services/RSATWS.wsdl"
	@echo "WS location"
	@grep location ${RSAT}/public_html/web_services/RSATWS.wsdl

## Init Web services
ws_stub: ws_init
	@echo
	@echo "Initiating Web services at SERVER=RSAT_WS=${RSAT_WS}"
	(cd ${RSAT}/${WS_STUB_CLIENT_DIR}; chmod 755 *.pl; make stub SERVER=${RSAT_WS})

ws_stub_test:
	@echo
	@echo "Testing Web services at SERVER=RSAT_WS=${RSAT_WS} using the stub"
	(cd ${RSAT}/${WS_STUB_CLIENT_DIR}/; make stub_test)

ws_nostub_test:
	@echo
	@echo "Testing Web services at SERVER=RSAT_WS=${RSAT_WS} without stub"
	perl ${RSAT}/${WS_PERL_CLIENT_DIR}/supported-organisms.pl ${RSAT_WS}

## Test web services to all RSAT servers
ws_test_all_servers:
	@echo
	@echo "Testing Web services at various RSAT servers"
	(cd ${RSAT}/${WS_PERL_CLIENT_DIR}; make test_all_servers)

################################################################
## Compile and install C/C++ programs that are part of the RSAT
## distribution (since 2009).

## Compile all programs
compile_all:  compile_info_gibbs compile_count_words compile_matrix_scan_quick compile_compare_matrices_quick compile_retrieve_variation_seq compile_variation_scan

compile_all_bugged: compile_retrieve_variation_seq compile_variation_scan compile_info_gibbs compile_count_words compile_matrix_scan_quick compile_compare_matrices_quick

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
	(cd ${SRC_DIR}; make clean; make all; ${SUDO} rsync -uptl ${SRC_DIR}/${PROGRAM} ${BIN}/)
#	(cd ${SRC_DIR}; make clean; make all; ${SUDO} rsync -ruptL ${SRC_DIR}/${PROGRAM} ${BIN}/)
#	(cd ${SRC_DIR}; make all; ln -fs ${SRC_DIR}/${PROGRAM} ${RSAT}/bin/${PROGRAM})
	@echo ${BIN}/${PROGRAM}
	@echo ""



## Compile and install retrieve-variation-seq
compile_retrieve_variation_seq:
	@${MAKE} compile_one_program PROGRAM=retrieve-variation-seq

## Compile and install variation-scan
compile_variation_scan:
	@${MAKE} compile_one_program PROGRAM=variation-scan

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
	${MAKE} -f ${RSAT}/makefiles/install_software.mk install_weblogo

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
	@(cd ${RSAT}/contrib/kwalks/src; make clean; make; \
		echo "	Installing lkwalk executable in bin directory ${BIN}"; \
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
## Restart the pathway tools services
## 
## Not really working as such, commands have to be typed after having
## entered in su -. I still should improve this.
restart_pathway_tools:
	(su - ; su - postgres (pg_ctl -D data -l logfile restart))
	(su -; cd /usr/share/tomcat6/bin; ./catalina.sh start)


