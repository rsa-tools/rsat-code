#!/usr/bin/perl -w
# compare-classes_client.pl - Client compare-classes using the SOAP::WSDL module
# and a property file

################################################################
##
## This script runs a simple demo of the web service interface to the
## RSAT tool compare-classes. It sends a request to the server for
## obtaining comparaison between two files.
##
################################################################

use strict;
use SOAP::WSDL;
use Util::Properties;
#import SOAP::Lite + trace;

## WSDL location
my $server = 'http://rsat.scmbb.ulb.ac.be/rsat/web_services';
#my $server = 'http://localhost/rsat/web_services';
my $WSDL = $server.'/RSATWS.wsdl';
my $proxy = $server.'/RSATWS.cgi';
my $property_file = shift @ARGV;
die "\tYou must specify the property file as first argument\n"
  unless $property_file;

## Service call
my $soap=SOAP::WSDL->new(wsdl => $WSDL)->proxy($proxy);
$soap->wsdlinit;

# $soap->wsdl_checkoccurs(0);

my $prop =  Util::Properties->new();
$prop->file_name($property_file);
$prop->load();
my %args = $prop->prop_list();

if ($args{ref_classes}) {
$args{ref_classes} = `cat $args{ref_classes}`;
chomp($args{ref_classes});
}

if ($args{query_classes}) {
$args{query_classes} = `cat $args{query_classes}`;
chomp($args{query_classes});
}

if ($args{input_classes}) {
$args{input_classes} = `cat $args{input_classes}`;
chomp($args{input_classes});
}

my $output_choice = $args{output_choice} || 'both';

warn "\nThis demo script compares two class files(the query and the reference files) Each class of the query file is compared to each class of the reference file\n\n";

## Send the request to the server
print "Sending request to the server $server\n";
my $som = $soap->call('compare_classes' => 'request' => \%args);

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
	print "Compared classes(s): \n".$result;
    } elsif ($output_choice eq 'both') {
	my $server_file = $results{'server'};
	my $result = $results{'client'};
	print "Result file on the server: ".$server_file;
	print "Compared classes(s): \n".$result;
    }
}
