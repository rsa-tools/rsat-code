#!/usr/bin/perl -w
# gene-info_client_soap-wsdl.pl - Client gene-info using the SOAP::WSDL module

################################################################
##
## This script runs a simple demo of the web service inerface to the
## RSAT tool gene-info. It sends a request to the server for
## obtaining information on 3 E.coli genes.
##
################################################################

use strict;
use SOAP::WSDL;
# import SOAP::Lite +trace;

## WSDL location
my $WSDL = 'http://rsat.scmbb.ulb.ac.be/rsat/web_services/RSATWS.wsdl';
my $proxy = 'http://rsat.scmbb.ulb.ac.be/rsat/web_services/RSATWS.cgi';

my $soap=SOAP::WSDL->new(wsdl => $WSDL)->proxy($proxy);

$soap->wsdlinit;

# $soap->wsdl_checkoccurs(0);

## Return option
my $return_choice = 'both';  ## Accepted values: 'file', 'result', 'both'

## Gene-info parameters
my $organism = 'Escherichia_coli_K12';  ## Name of the query organism
#my $organism = 'Escherichia_colix_K12';
#my @gene = ("metA", "metB", "metC");  ## List of query genes
# my @gene = ("zorglB");
my $query = 'metA';
my $full = 'full';  ## The -full option.
my $noquery = '';  ## 
my $descr = '';  ## Accepted value: 'descr'
my $feattype = '';  ## If the -feattype option value is not specified, the default is used

my %args = ('organism' => $organism,
#	    'query' => \@gene,  ## An array in a hash has to be referenced (correct?)
	    'query' => $query,
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
#    if ($return_choice eq 'file') {
#	my $server_file = $results{'file'};
#	print "Result file on the server: ".$server_file;
#    } elsif ($return_choice eq 'result') {
	my $result = $results{'result'};
	print "Gene(s) info(s): \n".$result;
#    } elsif ($return_choice eq 'both') {
#	my $server_file = $results{'file'};
#	my $result = $results{'result'};
#	print "Result file on the server: ".$server_file;
#	print "Retrieved sequence(s): \n".$result;
#    }
}
