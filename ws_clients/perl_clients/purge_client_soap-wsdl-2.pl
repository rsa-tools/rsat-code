#!/usr/bin/perl -w
# purge_client_soap-wsdl-2.pl - Client purge-sequence using the SOAP::WSDL module

################################################################
##
################################################################

use strict;
use SOAP::WSDL; ## Requires version 2.0 or later of SOAP::WSDL
use lib 'RSATWS';
use MyInterfaces::RSATWebServices::RSATWSPortType;

warn "\nThis demo script retrieves the start codons for a set of query genes\n\n";

## WSDL location
my $server = 'http://rsat.scmbb.ulb.ac.be/rsat/web_services';

## Service call
my $soap=MyInterfaces::RSATWebServices::RSATWSPortType->new();

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
my $som = $soap->purge_seq({'request' => \%args});

## Get the result
unless ($som) {
	printf "A fault (%s) occured: %s\n", $som->get_faultcode(), $som->get_faultstring();
} else {
	my $results = $som->get_response();

    ## Report the remote command
    my $command = $results -> get_command();
    print "Command used on the server: ".$command, "\n";

    ## Report the result
    if ($output_choice eq 'server') {
 		my $server_file = $results -> get_server();
		print "Result file on the server: ".$server_file;
   } elsif ($output_choice eq 'client') {
		my $result = $results -> get_client();
		print "Retrieved sequence(s): \n".$result;
   } elsif ($output_choice eq 'both') {
		my $server_file = $results -> get_server();
		my $result = $results -> get_client();
		print "Result file on the server: ".$server_file."\n";
		print "Retrieved sequence(s): \n".$result;
    }
}
