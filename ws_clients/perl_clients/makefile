################################################################ The
## command necessary to (re-)generate the Stubb for ensuring
## compatibility with the Web services

include ${RSAT}/makefiles/util.mk
MAKEFILE=makefile

SERVER=http://www.rsat.eu
#SERVER=http://rsat.ulb.ac.be/rsat
#SERVER=http://139.124.66.4/rsat
#SERVER=http://liv.bmc.uu.se/rsa-tools
#SERVER=http://bongcam1.hgen.slu.se/rsat
#SERVER=http://tagc.univ-mrs.fr/rsa-tools
#SERVER=http://embnet.ccg.unam.mx/rsa-tools
#SERVER=http://localhost/rsat
SERVERS= \
	http://www.rsat.eu \
	http://rsat.ulb.ac.be/rsat \
	http://139.124.66.4/rsat \
	http://liv.bmc.uu.se/rsa-tools \
	http://tagc.univ-mrs.fr/rsa-tools \
	http://embnet.ccg.unam.mx/rsa-tools \
	http://localhost/rsat
list_servers:
	@echo "SERVERS	${SERVERS}"
	@echo "SERVER	${SERVER}"
	@echo "WSDL_URL	${WSDL_URL}"

stubb_param:
	@echo "Parameters for generating the stubb"
	@echo "SERVERS	${SERVERS}"
	@echo "SERVER	${SERVER}"
	@echo "WSDL_URL	${WSDL_URL}"

WSDL_URL=${SERVER}/web_services/RSATWS.wsdl
stubb:
	@echo
	@echo "Generating stubb for server	${SERVER}"
	@echo "WSDL URL				${WSDL_URL}"
	perl ./wsdl2perl.pl -b RSATWS ${WSDL_URL}
	@echo "Generated stubb for server	${SERVER}"

test:
	perl ./convert_matrix_soap_wsdl_2.pl 
