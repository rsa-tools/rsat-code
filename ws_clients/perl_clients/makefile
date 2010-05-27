################################################################ The
## command necessary to (re-)generate the Stubb for ensuring
## compatibility with the Web services

include ${RSAT}/makefiles/util.mk
MAKEFILE=makefile

SERVER=rsat.ulb.ac.be
#SERVER=localhost
stubb:
	wsdl2perl.pl -b RSATWS http://${SERVER}/rsat/web_services/RSATWS.wsdl


test:
	./convert_matrix_soap_wsdl_2.pl 
