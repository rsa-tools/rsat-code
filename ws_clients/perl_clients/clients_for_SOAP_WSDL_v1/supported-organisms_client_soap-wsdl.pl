#!/usr/bin/perl -w
# supported-organisms_client_soap-wsdl.pl - Client supported-organisms using the SOAP::WSDL module

################################################################
##
## This script runs a simple demo of the web service interface to the
## RSAT tool supported-organisms.
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

## Supported-organisms parameters
my $format;  ## Output format (supported: html_list,html_table,array,text,keys,names,sizes,full,tree,html_tree)
my $taxon;  ## Root taxon

my %args = ('format' => $format,
	    'taxon' => $taxon,
	    );

## Send the request to the server
print "Sending request to the server $server\n";
my $som = $soap->call('supported_organisms' => 'request' => \%args);

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
