#!/usr/bin/perl -w
# gene-info_client.pl - Client gene-info using the SOAP::WSDL module and a property file

################################################################
##
## This script runs a simple demo of the web service interface to the
## RSAT tool gene-info. It searches all the genes havine the words
## "mehtionine" or "purine" in their description.
##
################################################################

use strict;
use SOAP::WSDL;
use Util::Properties;
#import SOAP::Lite+trace;

## Service location
my $server = 'http://rsat.scmbb.ulb.ac.be/rsat/web_services';
#my $server = 'http://localhost/rsat/web_services';
my $WSDL = $server.'/RSATWS.wsdl';
my $proxy = $server.'/RSATWS.cgi';

## Property file is firste argument
my $property_file = shift @ARGV;
die "\tYou must specify the property file as first argument\n"
  unless $property_file;

## Call the service
my $soap=SOAP::WSDL->new(wsdl => $WSDL)->proxy($proxy);
$soap->wsdlinit;

my $prop =  Util::Properties->new();
$prop->file_name($property_file);
$prop->load();
my %args = $prop->prop_list();

## Convert the query string into a list
my @queries = split(",", $args{query});
$args{query} = \@queries;

my $output_choice = $args{output_choice} || 'both';

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
