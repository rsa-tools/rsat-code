#!/usr/bin/perl -w
# convert-seq_client.pl - Client convert-sequence using the SOAP::WSDL module
# and a property file

################################################################
##
## This script runs a simple demo of the web service interface to the
## RSAT tool convert-sequence. It sends a request to the server for
## converting a sequence between two formats (e.g. fasta -> raw).
##
################################################################

use strict;
use SOAP::WSDL;
use Util::Properties;
#import SOAP::Lite +trace;

## WSDL location
my $server = 'http://rsat.scmbb.ulb.ac.be/rsat/web_services';
#my $server = 'http://localhost/rsat/web_services';
my $WSDL = $server.'/RSATWS.wsdl';
my $proxy = $server.'/RSATWS.cgi';

## Property file is firste argument
my $property_file = shift @ARGV;
die "\tYou must specify the property file as first argument\n"
  unless $property_file;

## Service call
my $soap=SOAP::WSDL->new(wsdl => $WSDL)->proxy($proxy);
$soap->wsdlinit;

my $prop =  Util::Properties->new();
$prop->file_name($property_file);
$prop->load();
my %args = $prop->prop_list();

## TEMPORARY: replace sequence file by sequence
$args{sequence} = `cat $args{sequence_file}`;
chomp($args{sequence});
delete($args{sequence_file});

my $output_choice = $args{output_choice} || 'both';

## Send the request to the server
print "Sending request to the server $server\n";
my $som = $soap -> call('convert_seq' => 'request' => \%args);

## Get the result
if ($som->fault){  ## Report error if any
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
	print "Converted sequence(s): \n".$result;
    } elsif ($output_choice eq 'both') {
	my $server_file = $results{'server'};
	my $result = $results{'client'};
	print "Result file on the server: ".$server_file;
	print "Converted sequence(s): \n".$result;
    }
}
