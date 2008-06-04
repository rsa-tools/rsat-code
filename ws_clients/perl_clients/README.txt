This folder contains the Perl libraries and sample clients for using
the Regulatory Sequence Analysis Tools Web Services (RSATWS).


REQUIREMENTS
============
This version of the RSATWS clients require a recent version of the
SOAP:WSDL library (later than v2.00.01, released in 2008). To install
this library, log in as a super-user and run the following CPAN commands
under cpan.

cpan> install Module::Build::Compat
cpan> install SOAP::WSDL

If you are using an earlier version of SOAP::WSDL and, for some reason, do not want to update it, Perl clients scripts working with the older versions of SOAP::WSDL can be found in the included subfolder clients_for_SOAP_WSDL_v1.

Note for Mac users
------------------
The cpan installation tests for SOAP::WSDL fail on the Macinstosh
OS10.5. We needed to force the cpan installation as follows:

cpan> force install SOAP::WSDL


TESTING
=======

You can either run the scripts directly from the directory where
they are installed, or from anywhere else, but in such case you need
to include the RSATWS stub in your Perl library (the environent
variable is $PERL5LIB).

For a simple test, run the scripts directly from their directory. 

perl ./retrieve-seq_client_soap-wsdl.pl


FILES
=====
README.txt	this file
RSATWS	RSAT Stub for the SOAP::WSDL library
*.pl		Perl clients

				   
