#!/usr/bin/perl -w
# display_graph_client.pl is a perl example of a RSAT display-graph using the SOAP::WSDL module.
# it creates a demo.jpg file representing the graph.

################################################################
##
## This script runs a simple demo of the web service interface to the
## RSAT tool display-graph. It sends a request to the server to display
## a graph in a user chosen format
## 
################################################################

#use strict;
use SOAP::WSDL;
use Util::Properties;
import SOAP::Lite + trace;

# # ## Service location
my $server = 'http://rsat.scmbb.ulb.ac.be/rsat/web_services';
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

$args{inputgraph} = `cat $args{inputgraph}`;


my $output_choice = $args{output_choice} || 'both';

warn "\nThis demo script displays a graph into a picture format (png/jpeg/ps) \n\n";

## Send the request to the server
print "Sending request to the server $server\n";
my $som = $soap->call('display_graph' => 'request' => \%args);

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
	open (DEMO, ">demo.jpg");
	print DEMO $result;
	close DEMO;
    } elsif ($output_choice eq 'both') {
	my $server_file = $results{'server'};
	my $result = $results{'client'};
	print "Result file on the server: ".$server_file;
	open (DEMO, ">demo.jpg");
	print DEMO $result;
	close DEMO;
    }
}
