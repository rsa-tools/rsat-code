#!/usr/bin/perl -w
# retrieve_purge_oligos_client.pl - Client retrieve-seq + oligo-analysis

################################################################
##
## This script runs a simple demo of the web service interface to the
## RSAT tools retrieve-seq, purge-sequence and oligo-analysis linked in a workflow.
##  It sends a request to the server for discovering 6 letter words
## in upstream sequences of 5 yeast genes. The sequences are first
## retrieved and purged for repeated segments
##
################################################################

use strict;
use SOAP::WSDL;
use Util::Properties;

warn "\nThis demo script illustrates a work flow combining three requests to the RSAT web services:\n\tretrieve-seq | purge-sequence | oligo-analysis\n\n";


## Service location
my $server = 'http://rsat.scmbb.ulb.ac.be/rsat/web_services';
my $WSDL = $server.'/RSATWS.wsdl';
my $proxy = $server.'/RSATWS.cgi';
my $property_file = shift @ARGV;
die "\tYou must specify the property file as first argument\n"
  unless $property_file;

## Service call
my $soap=SOAP::WSDL->new(wsdl => $WSDL)->proxy($proxy);
$soap->wsdlinit;

#################################################
## Retrieve-seq part

my $prop =  Util::Properties->new();
$prop->file_name($property_file);
$prop->load();
my %args = $prop->prop_list();

## Convert the query string into a list
my @queries = split(",", $args{query});
$args{query} = \@queries;

## Send request to the server
print "\nRetrieve-seq: sending request to the server\t", $server, "\n";
my $som = $soap->call('retrieve_seq' => 'request' => \%args);

## Get the result
my $server_file;  ## That variable needs to be declared outside the if..else block to be useable in the next part
if ($som->fault){  ## Report error if any
printf "A fault (%s) occured: %s\n", $som->faultcode, $som->faultstring;
} else {
my $results_ref = $som->result;  ## A reference to the result hash table
my %results = %$results_ref;  ## Dereference the result hash table

## Report the remote command
my $command = $results{'command'};
print "Command used on the server:\n\t".$command, "\n";

## Report the result file name on the server
$server_file = $results{'server'};
print "Result file on the server:\n\t".$server_file;
}

#################################################
## Purge-sequence part

$args{tmp_infile} = $server_file;

## Send the request to the server
print "\nPurge-sequence: sending request to the server\t", $server, "\n";
$som = $soap -> call('purge_seq' => 'request' => \%args);

## Get the result
if ($som->fault){  ## Report error if any
printf "A fault (%s) occured: %s\n", $som->faultcode, $som->faultstring;
} else {
my $results_ref = $som->result;  ## A reference to the result hash table
my %results = %$results_ref;  ## Dereference the result hash table

## Report the remote command
my $command = $results{'command'};
print "Command used on the server: \n\t".$command, "\n";

## Report the result file name on the server
$server_file = $results{'server'};
print "Result file on the server: \n\t".$server_file;
}
#################################################
## Oligo-analysis part

$args{output} = "both";
$args{tmp_infile} = $server_file;

## Convert the lth string into a list
my @lths = split(",", $args{lth});
$args{lth} = \@lths;

## Send request to the server
print "\nOligo-analysis: sending request to the server\t", $server, "\n";
$som = $soap->call('oligo_analysis' => 'request' => \%args);

## Get the result
if ($som->fault){  ## Report error if any
    printf "A fault (%s) occured: %s\n", $som->faultcode, $som->faultstring;
} else {
    my $results_ref = $som->result;
    my %results = %$results_ref;
    
    ## Report remote commande
    my $command = $results{'command'};
    print "Command used on the server: ".$command, "\n";
    
    ## Report the result
    $server_file = $results{'server'};
    my $result = $results{'client'};
    print "Result file on the server: \n\t".$server_file;
    print "Discovered oligo(s): \n".$result;
}
