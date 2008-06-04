#!/usr/bin/perl -w
# matrix-scan_client.pl - Client matrix-scan using the SOAP::WSDL module
# and a property file

################################################################
##
## This script runs a simple demo of the web service interface to the
## RSAT tool matrix-scan. It sends a request to the server for
## obtaining a tab-delimited file, with one row per match.
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

if ($args{sequence_file}) {
$args{sequence_file} = `cat $args{sequence_file}`;
chomp($args{sequence_file});
}

if ($args{matrix_file}) {
$args{matrix_file} = `cat $args{matrix_file}`;
chomp($args{matrix_file});
}

if ($args{matrix_list}) {
$args{matrix_list} = `cat $args{matrix_list}`;
chomp($args{matrix_list});
}

if ($args{background}) {
$args{background} = `cat $args{background}`;
chomp($args{background});
}

## Convert the lth string into a list
my @lths = split(",", $args{lth});
$args{lth} = \@lths;

## Convert the uth string into a list
my @uths = split(",", $args{uth});
$args{uth} = \@uths;

my $output_choice = $args{output_choice} || 'both';

warn "\nThis demo script scans sequences with one or several position-specific scoring matrices (PSSM) to identify instances of the corresponding motifs (putative sites).\n\n";

## Send the request to the server
print "Sending request to the server $server\n";
my $som = $soap->call('matrix_scan' => 'request' => \%args);

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
	print "Matrix scan(s): \n".$result;
    } elsif ($output_choice eq 'both') {
	my $server_file = $results{'server'};
	my $result = $results{'client'};
	print "Result file on the server: ".$server_file;
	print "Matrix scan(s): \n".$result;
    }
}
