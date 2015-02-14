################################################################ The
## command necessary to (re-)generate the Stub for ensuring
## compatibility with the Web services

include ${RSAT}/makefiles/util.mk
MAKEFILE=makefile

SERVER=http://pedagogix-tagc.univ-mrs.fr/rsat

## For these servers, the SOAP/WSDL should be working fine
SERVERS=http://pedagogix-tagc.univ-mrs.fr/rsat \
	http://tagc.univ-mrs.fr/rsa-tools \
	http://embnet.ccg.unam.mx/rsa-tools \
	http://rsat.ulb.ac.be/rsat \
	http://nexus.hgen.slu.se/rsa-tools

## This server should be working only if you have a local installation of RSAT
SERVERS_TO_CHECK = \
	http://localhost/rsat

## I have to understand why web sevices do not work whe I use
## redirections. It seems to be related to a .. in the path: since
## these servers have a direct link without the rsat or rsa-tools
## folder.
SERVERS_NOT_WORKING= \
	http://www.rsat.eu \
	http://prokaroyte.rsat.eu \
	http://fungi.rsat.eu \
	http://plants.rsat.eu \
	http://metazoa.rsat.eu \
	http://protists.rsat.eu \
	http://teaching.rsat.eu \
	http://training.rsat.eu

list_servers:
	@echo "SERVERS	${SERVERS}"
	@echo "SERVER	${SERVER}"
	@echo "WSDL_URL	${WSDL_URL}"
	@echo "RSAT_WS	${RSAT_WS}"

stub_param:
	@echo "Parameters for generating the stub"
	@echo "SERVERS	${SERVERS}"
	@echo "SERVER	${SERVER}"
	@echo "WSDL_URL	${WSDL_URL}"

WSDL_URL=${SERVER}/web_services/RSATWS.wsdl
stub:
	@echo
	@echo "Generating stub for server"
	@echo "	SERVER=${SERVER}"
	@echo "	WSDL_URL=${WSDL_URL}"
	perl ./wsdl2perl.pl -b RSATWS ${WSDL_URL}
	@echo "Generated stub for server"
	@echo "	SERVER=${SERVER}"
	@echo "	WSDL_URL=${WSDL_URL}"

stub_test:
	perl ./convert_matrix_soap_wsdl_2.pl 


nostub_test:
	@echo "Testing web services without stub	${RSAT_WS}"
	perl ./supported-organisms_client_nostub_wsdl.pl ${RSAT_WS}

nostub_test_all_servers:
	@for server in ${SERVERS}; do \
		${MAKE} RSAT_WS=$${server} nostub_test; \
	done

