#!/usr/bin/perl -w
# gene-info_client_soap-wsdl.pl - Client gene-info using the SOAP::WSDL module

################################################################
##
## This script runs a simple demo of the web service interface to the
## RSAT tool gene-info. It sends a request to the server for
## obtaining information on 3 E. coli genes.
##
################################################################

use strict;
use SOAP::WSDL;

## Service location
my $WSDL = 'http://rsat.scmbb.ulb.ac.be/rsat/web_services/RSATWS.wsdl';
my $proxy = 'http://rsat.scmbb.ulb.ac.be/rsat/web_services/RSATWS.cgi';

## Call the service
my $soap=SOAP::WSDL->new(wsdl => $WSDL)->proxy($proxy);
$soap->wsdlinit;

## Gene-info parameters
my $organism = 'Escherichia_coli_K12';  ## Name of the query organism
my @gene = ("metA", "metB", "metC");  ## List of query genes
my $full = 'full';  ## Looking for full match, not substring match.
my $noquery = '';  ## Not used here.
my $descr = '';  ## Accepted value: 'descr'. Not used here.
my $feattype = '';  ## If the -feattype option value is not specified, the default is used (CDS for E. coli)

my %args = ('organism' => $organism,
	    'query' => \@gene,
	    'full' => $full,
	    'noquery' => $noquery,
	    'descr' => $descr,
	    'feattype' => $feattype);

## Send the request to the server
print "Sending request to the server\n";
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
