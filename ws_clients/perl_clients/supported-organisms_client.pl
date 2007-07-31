#!/usr/bin/perl -w
# supported-organisms_client.pl - Client supported-organisms using the SOAP::WSDL module and a property file

################################################################
##
## This script runs a simple demo of the web service interface to the
## RSAT tool supported-organisms.
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

## Property file is first argument
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
    print "Supported organisms: \n".$result;
}
