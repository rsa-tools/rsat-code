#!/usr/bin/perl -w
# purge_client_soap-wsdl.pl - Client purge-sequence using the SOAP::WSDL module

################################################################
##
## This script runs a simple demo of the web service interface to the
## RSAT tool purge-sequence. It sends a request to the server for
## obtaining the purged sequences of 2 yeast upstream sequences.
## Segments repeated in the two sequences will be returned masked.
##
################################################################

use strict;
use SOAP::WSDL;
#import SOAP::Lite+trace;

## WSDL location
my $server = 'http://rsat.scmbb.ulb.ac.be/rsat/web_services';
#my $server = 'http://localhost/rsat/web_services';
my $WSDL = $server.'/RSATWS.wsdl';
my $proxy = $server.'/RSATWS.cgi';

## Service call
my $soap=SOAP::WSDL->new(wsdl => $WSDL)->proxy($proxy);
$soap->wsdlinit;

## Output option
my $output_choice = 'both';  ## Accepted values: 'server', 'client', 'both'

## Purge-sequence parameters
my $sequence = '>NP_009576.1    YBR020w; upstream from -800 to -1; size: 800; location: NC_001134.7 278221 279020 D
CAGGTTATCAGCAACAACACAGTCATATCCATTCTCAATTAGCTCTACCACAGTGTGTGAACCAATGTATCCAGCACCACCTGTAACCAAAACAATTTTAGAAGTACTTTCACTTTGTAACTGAGCTGTCATTTATATTGAATTTTCAAAAATTCTTACTTTTTTTTTGGATGGACGCAAAGAAGTTTAATAATCATATTACATGGCATTACCACCATATACATATCCATATCTAATCTTACTTATATGTTGTGGAAATGTAAAGAGCCCCATTATCTTAGCCTAAAAAAACCTTCTCTTTGGAACTTTCAGTAATACGCTTAACTGCTCATTGCTATATTGAAGTACGGATTAGAAGCCGCCGAGCGGGCGACAGCCCTCCGACGGAAGACTCTCCTCCGTGCGTCCTCGTCTTCACCGGTCGCGTTCCTGAAACGCAGATGTGCCTCGCGCCGCACTGCTCCGAACAATAAAGATTCTACAATACTAGCTTTTATGGTTATGAAGAGGAAAAATTGGCAGTAACCTGGCCCCACAAACCTTCAAATTAACGAATCAAATTAACAACCATAGGATGATAATGCGATTAGTTTTTTAGCCTTATTTCTGGGGTAATTAATCAGCGAAGCGATGATTTTTGATCTATTAACAGATATATAAATGGAAAAGCTGCATAACCACTTTAACTAATACTTTCAACATTTTCAGTTTGTATTACTTCTTATTCAAATGTCATAAAAGTATCAACAAAAAATTGTTAATATACCTCTATACTTTAACGTCAAGGAGAAAAAACTATA
>NP_009575.1    YBR019c; upstream from -800 to -1; size: 800; location: NC_001134.7 278353 279152 R
CGGTTTAGCATCATAAGCGCTTATAAATTTCTTAATTATGCTCGGGCACTTTTCGGCCAATGGTCTTGGTAATTCCTTTGCGCTAGAATTGAACTCAGGTACAATCACTTCTTCTGAATGAGATTTAGTCATTATAGTTTTTTCTCCTTGACGTTAAAGTATAGAGGTATATTAACAATTTTTTGTTGATACTTTTATGACATTTGAATAAGAAGTAATACAAACTGAAAATGTTGAAAGTATTAGTTAAAGTGGTTATGCAGCTTTTCCATTTATATATCTGTTAATAGATCAAAAATCATCGCTTCGCTGATTAATTACCCCAGAAATAAGGCTAAAAAACTAATCGCATTATCATCCTATGGTTGTTAATTTGATTCGTTAATTTGAAGGTTTGTGGGGCCAGGTTACTGCCAATTTTTCCTCTTCATAACCATAAAAGCTAGTATTGTAGAATCTTTATTGTTCGGAGCAGTGCGGCGCGAGGCACATCTGCGTTTCAGGAACGCGACCGGTGAAGACGAGGACGCACGGAGGAGAGTCTTCCGTCGGAGGGCTGTCGCCCGCTCGGCGGCTTCTAATCCGTACTTCAATATAGCAATGAGCAGTTAAGCGTATTACTGAAAGTTCCAAAGAGAAGGTTTTTTTAGGCTAAGATAATGGGGCTCTTTACATTTCCACAACATATAAGTAAGATTAGATATGGATATGTATATGGTGGTAATGCCATGTAATATGATTATTAAACTTCTTTGCGTCCATCCAAAAAAAAAGTAAGAATTTTTGAAAATTCAATATAA';
my $format = '';  ## Default sequence format is used
my $match_length;  ## Default match length (40) is used
my $mismatch;  ## Default number of mismatch (3) is used
my $str = '';  ## Discard duplications on direct or both strands; default is used (both strands)
my $delete = 0;  ## Other example: 1, to delete repeats instead of masking them
my $mask_short;  ## mask_short option to mask sequences shorter than specified length

my %args = ('output'=> $output_choice,
	    'sequence'=> $sequence,
	    'format'=> $format,
	    'match_length'=> $match_length,
	    'mismatch'=> $mismatch,
	    'str'=> $str,
	    'delete'=> $delete,
	    'mask_short'=> $mask_short);

## Send the request to the server
print "Sending request to the server $server\n";
my $som = $soap -> call('purge_seq' => 'request' => \%args);

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
	print "Purged sequence(s): \n".$result;
    } elsif ($output_choice eq 'both') {
	my $server_file = $results{'server'};
	my $result = $results{'client'};
	print "Result file on the server: ".$server_file;
	print "Purged sequence(s): \n".$result;
    }
}
