#!/usr/bin/perl -w
# retrieve-seq_client_soap-wsdl-2.pl - Client retrieve-seq using the SOAP::WSDL module

################################################################
##
## This script runs a simple demo of the web service inerface to the
## RSAT tool retrieve-seq. It sends a request to the server for
## obtaining the start codons of 3 E.coli genes.
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

## Retrieve-seq parameters
my $organism = 'Escherichia_coli_K12';  ## Name of the query organism
my @gene = ("metA", "metB", "metC");  ## List of query genes
my $all = 0;  ## the -all option (other accepted value = 1). This option is incompatible with the query list @gene (above)
my $noorf = 1;  ## Clip sequences to avoid upstream ORFs
my $from = 0;  ## Start position of the sequence
my $to = 2;  ## End position of the sequence
my $feattype = '';  ## The -feattype option value is  not specified, the default is used
my $type = '';  ## The -type option value; other example:'-type downstream'
my $format = '';  ## The -format option value. We use the default (fasta), but other formats could be specified, for example 'multi'
my $lw = 0;  ## Line width. 0 means all on one line
my $label = 'id,name';  ## Choice of label for the retrieved sequence(s)
my $label_sep = '';  ## Choice of separator for the label(s) of the retrieved sequence(s)
my $nocom = 0;  ## Other possible value = 1, to get sequence(s) whithout comments
my $repeat =  0;  ## Other possible value = 1, to have annotated repeat regions masked
my $imp_pos = 0;  ## Admit imprecise position (value = 1 to do so)

my %args = (
	    'output' => $output_choice,
	    'organism' => $organism,
	    'query' => \@gene,  ## An array in a hash has to be referenced (always?)
	    'noorf' => $noorf,
	    'from' => $from,
	    'to' => $to,
	    'feattype' => $feattype,
	    'type' => $type,
	    'format' => $format,
	    'lw' => $lw,
	    'label' => $label,
	    'label_sep' => $label_sep,
	    'nocom' => $nocom,
	    'repeat' => $repeat,
	    'imp_pos' => $imp_pos
	    );

## Send the request to the server
print "Sending request to the server $server\n";
my $som = $soap->retrieve_seq({'request' => \%args});

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
