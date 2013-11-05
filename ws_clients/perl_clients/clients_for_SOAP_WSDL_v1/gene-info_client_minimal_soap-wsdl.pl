#!/usr/bin/perl -w
# gene-info_client_minimal_soap-wsdl.pl - Client gene-info using the SOAP::WSDL module.

################################################################
##
## This script runs a simple demo of the web service interface to the
## RSAT tool gene-info. It sends a query to the web serice to ask all
## the yeast genes whose name starts with MET, followed by one or
## several numbers.
##
################################################################

use strict;
use SOAP::WSDL;

## Service location
my $server = 'http://rsat.scmbb.ulb.ac.be/rsat/web_services';
my $WSDL = $server.'/RSATWS.wsdl';
my $proxy = $server.'/RSATWS.cgi';

## Call the service
my $soap=SOAP::WSDL->new(wsdl => $WSDL)->proxy($proxy);
$soap->wsdlinit;

## Gene-info parameters
my $organism = 'Saccharomyces_cerevisiae';  ## Name of the query organism
my @queries = ('MET\d+');  ## This query is a regular expression
my $full = 1;  ## Looking for full match, not substring match.

my %args = ('organism' => $organism,
	    'query' => \@queries,
	    'full' => $full);

## Send the request to the server
warn "\nSending request to the server $server\n";
my $call = $soap->call('gene_info' => 'request' => \%args);

## Get the result
if ($call->fault){ ## Report error if any
    die (sprintf "\nA fault (%s) occured: %s\n", $call->faultcode, $call->faultstring);
} else {
    my $results_ref = $call->result;  ## A reference to the result hash table
    my %results = %$results_ref;  ## Dereference the result hash table

    ## Report the remote command
    my $command = $results{'command'};
    print "\nCommand used on the server: ".$command, "\n";

    ## Report the result
    my $result = $results{'client'};
    print "\nGene(s) info(s): \n".$result;
}
