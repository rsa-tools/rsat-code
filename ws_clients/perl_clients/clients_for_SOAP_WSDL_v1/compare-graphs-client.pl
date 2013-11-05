#!/usr/bin/perl -w
# compare-graphs-client.pl is a perl example of a RSAT convert-graph using the SOAP::WSDL module

################################################################
##
## This script runs a simple demo of the web service interface to the
## RSAT tool compare-graphs. It sends a request to the server for
## comparing two graphs
##
################################################################

#use strict;
use SOAP::WSDL;
use Util::Properties;
# import SOAP::Lite + trace;

## Service location
my $server = 'http://rsat.scmbb.ulb.ac.be/rsat/web_services';
#my $server = 'http://localhost/rsat/web_services';
my $WSDL = $server.'/RSATWS.wsdl';
my $proxy = $server.'/RSATWS.cgi';

## Property file is first argument
my $property_file = shift @ARGV;
die "\tYou must specify the property file as first argument\n"
  unless $property_file;

# Service call
my $soap=SOAP::WSDL->new(wsdl => $WSDL)->proxy($proxy);
$soap->wsdlinit;

my $prop =  Util::Properties->new();
$prop->file_name($property_file);
$prop->load();
my %args = $prop->prop_list();

#ARGUMENTS
$args{Rinputgraph} = `cat $args{Rinputgraph}`;
$args{Qinputgraph} = `cat $args{Qinputgraph}`;

my $output_choice = $args{output_choice} || 'both';

warn "\nThis demo script returns the intersection of two graphs \n\n";

## Send the request to the server
print "Sending request to the server $server\n";
my $som = $soap->call('compare_graphs' => 'request' => \%args);

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
    if ($output_choice eq 'server') {
	my $server_file = $results{'server'};
	print "Result file on the server: ".$server_file;
    } elsif ($output_choice eq 'client') {
	my $result = $results{'client'};
	print "Clusters\n".$result;
    } elsif ($output_choice eq 'both') {
	my $server_file = $results{'server'};
	my $result = $results{'client'};
	print "Result file on the server: ".$server_file;
	print "Conversion \n".$result;
    }
}
