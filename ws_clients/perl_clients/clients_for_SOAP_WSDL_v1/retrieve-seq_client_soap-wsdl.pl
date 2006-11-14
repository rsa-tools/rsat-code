#!/usr/bin/perl -w
# retrieve-seq_client_soap-wsdl.pl - Client retrieve-seq using the SOAP::WSDL module

################################################################
##
## This script runs a simple demo of the web service inerface to the
## RSAT tool retrieve-seq. It sends a request to the server for
## obtaining the start codons of 3 E.coli genes.
##
################################################################

use strict;
use SOAP::WSDL;
# import SOAP::Lite +trace;

## WSDL location
my $WSDL = 'http://rsat.scmbb.ulb.ac.be/rsat/web_services/RSATWS.wsdl';
my $proxy = 'http://rsat.scmbb.ulb.ac.be/rsat/web_services/RSATWS.cgi';

my $soap=SOAP::WSDL->new(wsdl => $WSDL)->proxy($proxy);

$soap->wsdlinit;

# $soap->wsdl_checkoccurs(0);

## Return option
my $return_choice = 'both';  ## Accepted values: 'file', 'result', 'both'

## Retrieve-seq parameters
my $organism = '-org Escherichia_coli_K12';  ## Name of the query organism
# my $organism = '-org Escherichia_colix_K12';
my @gene = ("metA", "metB", "metC");  ## List of query genes
# my @gene = ("zorglB");
my $all = '';  ## the -all option. This option is incompatible with the query list @gene (above)
my $noorf = '-noorf';  ## Clip sequences to avoid upstream ORFs
my $from = '-from 0';  ## Start position of the sequence
my $to = '-to 2';  ## End position of the sequence
my $feattype = '';  ## The -feattype option value is  not specified, the default is used
my $type = '';  ## The -type option value; other example:'-type downstream'
my $format = '';  ## The -format option value. We use the default (fasta), but other formats could be specified, for example 'multi'
my $label = '-label id,name';  ## Choice of label for the retrieved sequence(s)
my $label_sep = '';  ## Choice of separator for the label(s) of the retrieved sequence(s)
my $nocom = '';  ## Other possible value = '-nocom', to get sequence(s) whithout comments
my $repeat =  '';  ## Other possible value = '-rm', to have annotated repeat regions masked

my %args = ('return' => $return_choice,
	    'organism' => $organism,
	    'query' => \@gene,  ## An array in a hash has to be referenced (correct?)
	    'noorf' => $noorf,
	    'from' => $from,
	    'to' => $to,
	    'feattype' => $feattype,
	    'type' => $type,
	    'format' => $format,
	    'label' => $label,
	    'label_sep' => $label_sep,
	    'nocom' => $nocom,
	    'repeat' => $repeat);

## Send the request to the server
print "Sending request to the server\n";
my $som = $soap->call('retrieve_seq' => 'request' => \%args);

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
    if ($return_choice eq 'file') {
	my $server_file = $results{'file'};
	print "Result file on the server: ".$server_file;
    } elsif ($return_choice eq 'result') {
	my $result = $results{'result'};
	print "Retrieved sequence(s): \n".$result;
    } elsif ($return_choice eq 'both') {
	my $server_file = $results{'file'};
	my $result = $results{'result'};
	print "Result file on the server: ".$server_file;
	print "Retrieved sequence(s): \n".$result;
    }
}
