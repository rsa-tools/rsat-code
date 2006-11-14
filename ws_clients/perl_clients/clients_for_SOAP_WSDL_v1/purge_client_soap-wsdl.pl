#!/usr/bin/perl -w
# purge-seq_client_wsdl.pl - Client purge-sequence using the SOAP::WSDL module

################################################################
##
## This script runs a simple demo of the web service inerface to the
## RSAT tool purge-sequence. It sends a request to the server for
## obtaining the purged sequences of 2 yeast upstream sequences.
## Segments repeated in the two sequences will be returned masked.
##
################################################################

use strict;
use SOAP::WSDL;

## WSDL location
my $WSDL = 'http://rsat.scmbb.ulb.ac.be/rsat/web_services/RSATWS.wsdl';
my $proxy = 'http://rsat.scmbb.ulb.ac.be/rsat/web_services/RSATWS.cgi';

my $soap=SOAP::WSDL->new(wsdl => $WSDL)->proxy($proxy);

$soap->wsdlinit;

## Return option
my $return_choice = 'both';  ## Accepted values: 'file', 'result', 'both'

## Purge-sequence parameters
my $sequence = '>NP_009576.1    YBR020w; upstream from -800 to -1; size: 800; location: NC_001134.7 278221 279020 D
CAGGTTATCAGCAACAACACAGTCATATCCATTCTCAATTAGCTCTACCACAGTGTGTGAACCAATGTATCCAGCACCACCTGTAACCAAAACAATTTTAGAAGTACTTTCACTTTGTAACTGAGCTGTCATTTATATTGAATTTTCAAAAATTCTTACTTTTTTTTTGGATGGACGCAAAGAAGTTTAATAATCATATTACATGGCATTACCACCATATACATATCCATATCTAATCTTACTTATATGTTGTGGAAATGTAAAGAGCCCCATTATCTTAGCCTAAAAAAACCTTCTCTTTGGAACTTTCAGTAATACGCTTAACTGCTCATTGCTATATTGAAGTACGGATTAGAAGCCGCCGAGCGGGCGACAGCCCTCCGACGGAAGACTCTCCTCCGTGCGTCCTCGTCTTCACCGGTCGCGTTCCTGAAACGCAGATGTGCCTCGCGCCGCACTGCTCCGAACAATAAAGATTCTACAATACTAGCTTTTATGGTTATGAAGAGGAAAAATTGGCAGTAACCTGGCCCCACAAACCTTCAAATTAACGAATCAAATTAACAACCATAGGATGATAATGCGATTAGTTTTTTAGCCTTATTTCTGGGGTAATTAATCAGCGAAGCGATGATTTTTGATCTATTAACAGATATATAAATGGAAAAGCTGCATAACCACTTTAACTAATACTTTCAACATTTTCAGTTTGTATTACTTCTTATTCAAATGTCATAAAAGTATCAACAAAAAATTGTTAATATACCTCTATACTTTAACGTCAAGGAGAAAAAACTATA
>NP_009575.1    YBR019c; upstream from -800 to -1; size: 800; location: NC_001134.7 278353 279152 R
CGGTTTAGCATCATAAGCGCTTATAAATTTCTTAATTATGCTCGGGCACTTTTCGGCCAATGGTCTTGGTAATTCCTTTGCGCTAGAATTGAACTCAGGTACAATCACTTCTTCTGAATGAGATTTAGTCATTATAGTTTTTTCTCCTTGACGTTAAAGTATAGAGGTATATTAACAATTTTTTGTTGATACTTTTATGACATTTGAATAAGAAGTAATACAAACTGAAAATGTTGAAAGTATTAGTTAAAGTGGTTATGCAGCTTTTCCATTTATATATCTGTTAATAGATCAAAAATCATCGCTTCGCTGATTAATTACCCCAGAAATAAGGCTAAAAAACTAATCGCATTATCATCCTATGGTTGTTAATTTGATTCGTTAATTTGAAGGTTTGTGGGGCCAGGTTACTGCCAATTTTTCCTCTTCATAACCATAAAAGCTAGTATTGTAGAATCTTTATTGTTCGGAGCAGTGCGGCGCGAGGCACATCTGCGTTTCAGGAACGCGACCGGTGAAGACGAGGACGCACGGAGGAGAGTCTTCCGTCGGAGGGCTGTCGCCCGCTCGGCGGCTTCTAATCCGTACTTCAATATAGCAATGAGCAGTTAAGCGTATTACTGAAAGTTCCAAAGAGAAGGTTTTTTTAGGCTAAGATAATGGGGCTCTTTACATTTCCACAACATATAAGTAAGATTAGATATGGATATGTATATGGTGGTAATGCCATGTAATATGATTATTAAACTTCTTTGCGTCCATCCAAAAAAAAAGTAAGAATTTTTGAAAATTCAATATAA';
my $format = '';  ## Default sequence format is sused
my $match_length = '';  ## Default match length (40) is used
my $mismatch = '';  ## Default number of mismatch (3) is used
my $str = '';  ## Discard duplications on direct or both strands; default is used
my $delete = '';  ## Other example: '-del' to delete repeats instead of masking them
my $mask_short = '';  ## -mask_short option to mas sequences shorter than specified length

my %args = ('return'=> $return_choice,
	    'sequence'=> $sequence,
	    'format'=> $format,
	    'match_length'=> $match_length,
	    'mismatch'=> $mismatch,
	    'str'=> $str,
	    'delete'=> $delete,
	    'mask_short'=> $mask_short);

## Send the request to the server
print "Sending request to the server\n";
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
    if ($return_choice eq 'file') {
	my $server_file = $results{'file'};
	print "Result file on the server: ".$server_file;
    } elsif ($return_choice eq 'result') {
	my $result = $results{'result'};
	print "Purged sequence(s): \n".$result;
    } elsif ($return_choice eq 'both') {
	my $server_file = $results{'file'};
	my $result = $results{'result'};
	print "Result file on the server: ".$server_file;
	print "Purged sequence(s): \n".$result;
    }
}
