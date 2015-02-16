#!/usr/bin/perl -w
################################################################
##
## This script runs a simple demo of the web service inerface to the
## RSAT tool retrieve-seq. It sends a request to the server for
## obtaining the start codons of 3 E.coli genes.
##
################################################################

#use strict;
use SOAP::WSDL;
use Util::Properties;

## WSDL location
my $server = 'http://rsat.scmbb.ulb.ac.be/rsat/web_services';
my $property_file = shift @ARGV;
die "\tYou must specify the property file as first argument\n"
  unless $property_file;

## Service call
my $soap=SOAP::WSDL->new(wsdl => $server.'/RSATWS.wsdl')->proxy($server.'/RSATWS.cgi');
$soap->wsdlinit;

## Read properties
my $prop =  Util::Properties->new();
$prop->file_name($property_file);
$prop->load();
my %args = $prop->prop_list();
$args{output_choice} = 'both'; ## Force both outputs for this minimal script

## Convert the query string into a list
my @queries = split(",", $args{query});
$args{query} = \@queries;

## Send the request to the server
warn "Sending request to the server $server\n";
my $som = $soap->call('retrieve_seq' => 'request' => \%args);

## Get the result
if ($som->fault){ ## Report error if any
    die sprintf( "A fault (%s) occured: %s\n", $som->faultcode, $som->faultstring);
} else {
    my $results_ref = $som->result;  ## A reference to the result hash table
    my %results = %$results_ref;  ## Dereference the result hash table

    ## Report the remote command and results
    print "Command used on the server: ".$results{'command'}, "\n";
    print "Result file on the server: ", $results{'server'};
    print "Retrieved sequence(s): \n", $results{'client'};
}
