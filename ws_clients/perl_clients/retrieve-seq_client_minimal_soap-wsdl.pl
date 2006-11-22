#!/usr/bin/perl
# retrieve-seq_client_soap-wsdl.pl - Client retrieve-seq using the SOAP::WSDL module

################################################################
##
## This script runs a simple demo of the web service inerface to the
## RSAT tool retrieve-seq. It sends a request to the server for
## obtaining the start codons of 3 E.coli genes.
##
################################################################

use strict;
use SOAP::WSDL;
use SOAP::LITE;
# import SOAP::Lite +trace;

## WSDL location
#my $WSDL = 'http://rsat.scmbb.ulb.ac.be/rsat/web_services/RSATWS.wsdl';
#my $proxy = 'http://rsat.scmbb.ulb.ac.be/rsat/web_services/RSATWS.cgi';
#my $WSDL = 'http://localhost/rsa-tools/web_services/RSATWS.wsdl';
#my $proxy = 'http://localhost/rsa-tools/web_services/RSATWS.cgi';
#my $WSDL = 'http://embnet.ccg.unam.mx/rsa-tools/web_services/RSATWS.wsdl';
#my $proxy = 'http://embnet.ccg.unam.mx/rsa-tools/web_services/RSATWS.cgi';
my $WSDL = 'http://kinich.ccg.unam.mx/rsa-tools/web_services/RSATWS.wsdl';
my $proxy = 'http://kinich.ccg.unam.mx/rsa-tools/web_services/RSATWS.cgi';

my $soap = SOAP::WSDL->new(wsdl => $WSDL)->proxy($proxy);
$soap->wsdlinit();
# $soap->wsdl_checkoccurs(0);

my %args = ('return' => 'both', ## Store the result of the server + return it directly
	    'organism' => '-org Saccharomyces_cerevisiae',
	    'query' => ['PHO5', 'PHO8', 'PHO11', 'PHO84', 'PHO80'], 
	    'noorf' => '-noorf');

## Send the request to the server
print "Sending request to the server\t",$proxy,"\n";
my $som = $soap->call('retrieve_seq' => 'request' => \%args);

## Get the result
if ($som->fault){ ## Report error if any
    printf "A fault (%s) occured: %s\n", $som->faultcode, $som->faultstring;
} else {
    my $results_ref = $som->result;  ## A reference to the result hash table
    my %results = %$results_ref;  ## Dereference the result hash table

    ## Report the remote command
    my $command = $results{'command'};
    print "Command used on the server:\n\t", $command, "\n";
    print "Result file on the server:\n\t",  $results{'file'};
    print "Retrieved sequence(s):\n\n", $results{'result'};
}
