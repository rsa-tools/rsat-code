################################################################
## Check that web services are working at different RSAT servers.

include ${RSAT}/makefiles/util.mk
MAKEFILE=makefile

SERVER=http://pedagogix-tagc.univ-mrs.fr/rsat

## For these servers, the SOAP/WSDL should be working fine
SERVERS=http://embnet.ccg.unam.mx/rsat \
	http://rsat.france-bioinformatique.fr/rsat \
	http://rsat.france-bioinformatique.fr/fungi \
	http://rsat.france-bioinformatique.fr/metazoa \
	http://rsat.france-bioinformatique.fr/protists \
	http://rsat.france-bioinformatique.fr/teaching \
	http://rsat01.biologie.ens.fr/rsat \
	http://floresta.eead.csic.es/rsat


# SERVERS=http://pedagogix-tagc.univ-mrs.fr/rsat \
# 	http://rsat-tagc.univ-mrs.fr/rsat \
# 	http://floresta.eead.csic.es/rsat \
# 	http://rsat01.biologie.ens.fr/rsat \
# 	http://embnet.ccg.unam.mx/rsa-tools \
# 	http://rsat.ulb.ac.be/rsat \
# 	http://wwwsup.scmbb.ulb.ac.be/rsat

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
	http://training.rsat.eu \
	http://rsat.sb-roscoff.fr \
	http://tagc.univ-mrs.fr/rsa-tools \
	http://nexus.hgen.slu.se/rsa-tools

WSDL_URL=${SERVER}/web_services/RSATWS.wsdl
list_param:
	@echo "SERVERS	${SERVERS}"
	@echo "SERVER	${SERVER}"
	@echo "RSAT_WS	${RSAT_WS}"
	@echo "WSDL_URL	${WSDL_URL}"

# stub_param:
# 	@echo "Parameters for generating the stub"
# 	@echo "SERVERS	${SERVERS}"
# 	@echo "SERVER	${SERVER}"
# 	@echo "WSDL_URL	${WSDL_URL}"

# stub:
# 	@echo
# 	@echo "Generating stub for server"
# 	@echo "	SERVER=${SERVER}"
# 	@echo "	WSDL_URL=${WSDL_URL}"
# 	perl ./wsdl2perl.pl -b RSATWS ${WSDL_URL}
# 	@echo "Generated stub for server"
# 	@echo "	SERVER=${SERVER}"
# 	@echo "	WSDL_URL=${WSDL_URL}"

# stub_test:
# 	perl ./convert_matrix_soap_wsdl_2.pl 


test_one_server:
	@echo
	@echo "Testing web services (without stub) at URL	${RSAT_WS}"
	perl ./supported-organisms.pl ${RSAT_WS}

test_all_servers:
	@for server in ${SERVERS}; do \
		${MAKE} RSAT_WS=$${server} test_one_server; \
	done

