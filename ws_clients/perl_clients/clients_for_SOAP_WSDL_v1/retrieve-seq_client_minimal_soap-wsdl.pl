#!/usr/bin/perl

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
#my $server = 'http://rsat.scmbb.ulb.ac.be/rsat/web_services/';
my $server = 'http://localhost/rsat/web_services/';

my $WSDL = $server.'RSATWS.wsdl';
my $proxy = $server.'RSATWS.cgi';
my $soap = SOAP::WSDL->new(wsdl => $WSDL)->proxy($proxy);
$soap->wsdlinit();
# $soap->wsdl_checkoccurs(0);

my %args = ('return' => 'both', ## Store the result of the server + return it directly
	    'organism' => '-org Escherichia_coli_K12',
	    'query' => ['CRP', 'FruR'], 
	    'noorf' => '-noorf');

## Send the request to the server
print "Sending request to the server\t",$server,"\n";
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
