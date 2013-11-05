#!/usr/bin/perl -w
# oligos-analysis_client.pl - Client oligo-analysis using the SOAP::WSDL module and a property file

################################################################
##
## This script runs a simple demo of the web service interface to the
## RSAT tool oligo-analysis. It sends a request to the server for
## discovering 6 letter words in the upstream sequences of 5 yeast genes.
##
################################################################

use strict;
use SOAP::WSDL;
use Util::Properties;

warn "\nINFO: This demo script sends a set of sequences to the RSAT web service, and runs oligo-analysis to detect over-represented oligonuclotides\n\n";

## WSDL location
my $server = 'http://rsat.scmbb.ulb.ac.be/rsat/web_services';
#my $server = 'http://localhost/rsat/web_services';
my $WSDL = $server.'/RSATWS.wsdl';
my $proxy = $server.'/RSATWS.cgi';

## Property file is firste argument
my $property_file = shift @ARGV;
die "\tYou must specify the property file as first argument\n"
  unless $property_file;

my $soap=SOAP::WSDL->new(wsdl => $WSDL)->proxy($proxy);
$soap->wsdlinit;

my $prop =  Util::Properties->new();
$prop->file_name($property_file);
$prop->load();
my %args = $prop->prop_list();

my $sequence = `cat $args{sequence}`;
chomp($sequence);

$args{sequence} = $sequence;

## Convert the lth string into a list
my @lths = split(",", $args{lth});
$args{lth} = \@lths;

## Output option
my $output_choice = $args{output_choice} || 'both';

## Send request to the server
print "Sending request to the server $server\n";
my $som = $soap->call('oligo_analysis' => 'request' => \%args);

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
	print "Discovered oligo(s): \n".$result;
    } elsif ($output_choice eq 'both') {
	my $server_file = $results{'server'};
	my $result = $results{'client'};
	print "Result file on the server: ".$server_file;
	print "Discovered oligo(s): \n".$result;
    }
}
