################################################################ The
## command necessary to (re-)generate the Stubb for ensuring
## compatibility with the Web services

include ${RSAT}/makefiles/util.mk
MAKEFILE=makefile

SERVER=http://139.124.66.4/rsat
#http://rsat.ulb.ac.be/rsat
#SERVER=http://mamaze.ulb.ac.be/rsat
#SERVER=http://liv.bmc.uu.se/rsa-tools
#SERVER=http://bongcam1.hgen.slu.se/rsat/
#SERVER=http://tagc.univ-mrs.fr/rsa-tools
#SERVER=http://embnet.ccg.unam.mx/rsa-tools
#SERVER=http://localhost/rsat
SERVERS= \
	http://rsat.ulb.ac.be/rsat \
	http://liv.bmc.uu.se/rsa-tools \
	http://tagc.univ-mrs.fr/rsa-tools \
	http://embnet.ccg.unam.mx/rsa-tools \
	http://localhost/rsat
list_servers:
	@echo "SERVER	${SERVER}"
	@echo "SERVERS	${SERVERS}"

stubb:
	@echo "Generating stubb for server	${SERVER}"
	perl ./wsdl2perl.pl -b RSATWS ${SERVER}/web_services/RSATWS.wsdl
	@echo "Generated stubb for server	${SERVER}"

test:
	./convert_matrix_soap_wsdl_2.pl 
