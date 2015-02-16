#!/usr/bin/perl -w
# gene-info_client_soap-wsdl.pl - Client gene-info using the SOAP::WSDL module

################################################################
##
## This script runs a simple demo of the web service interface to the
## RSAT tool gene-info. It searches all the genes havine the words
## "mehtionine" or "purine" in their description.
##
################################################################

use strict;
use SOAP::WSDL;
#import SOAP::Lite+trace;

## Service location
my $server = 'http://rsat.scmbb.ulb.ac.be/rsat/web_services';
#my $server = 'http://localhost/rsat/web_services';
my $WSDL = $server.'/RSATWS.wsdl';
my $proxy = $server.'/RSATWS.cgi';

## Call the service
my $soap=SOAP::WSDL->new(wsdl => $WSDL)->proxy($proxy);
$soap->wsdlinit;

## Gene-info parameters
my $organism = 'Escherichia_coli_K12';  ## Name of the query organism
my @queries = ("methionine", "purine");  ## List of queries
my $full = 0;  ## Also looking for substring matches, not only full string matches.
my $descr = 1;  ## Also looking in description field, not just gene name

my %args = ('organism' => $organism,
	    'query' => \@queries,
	    'full' => $full,
	    'descr' => $descr);

## Send the request to the server
print "Sending request to the server $server\n";
my $som = $soap->call('gene_info' => 'request' => \%args);

## Get the result
if ($som->fault){ ## Report error if any
    printf "A fault (%s) occured: %s\n", $som->faultcode, $som->faultstring;
} else {
    my $results_ref = $som->result;  ## A reference to the result hash table
    my %results = %$results_ref;  ## Dereference the result hash table

    ## Report the remote command
    my $command = $results{'command'};
    print "Command used on the server: ".$command, "\n";

    ## Report the result
    my $result = $results{'client'};
    print "Gene(s) info(s): \n".$result;
}
