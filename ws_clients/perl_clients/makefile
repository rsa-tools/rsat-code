################################################################ The
## command necessary to (re-)generate the Stubb for ensuring
## compatibility with the Web services

include ${RSAT}/makefiles/util.mk
MAKEFILE=makefile

stubb:
	wsdl2perl.pl -b RSATWS http://rsat.ulb.ac.be/rsat/web_services/RSATWS.wsdl
