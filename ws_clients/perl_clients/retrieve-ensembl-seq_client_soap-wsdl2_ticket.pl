#!/usr/bin/perl -w
# retrieve-ensembl-seq_client_soap-wsdl-2.pl - Client retrieve-seq using the SOAP::WSDL module

################################################################
##
##
################################################################

use strict;
use SOAP::WSDL; ## Requires version 2.0 or later of SOAP::WSDL
use lib 'RSATWS';
use MyInterfaces::RSATWebServices::RSATWSPortType;

## WSDL location
my $server = 'http://rsat.scmbb.ulb.ac.be/rsat/web_services';

## Service call
my $soap=MyInterfaces::RSATWebServices::RSATWSPortType->new();

## Output option
#my $output_choice = 'client';  ## Accepted values: 'server', 'client', 'both'
my $output_choice = 'ticket';  ## Accepted values: 'server', 'client', 'both'

## Retrieve-seq parameters
my $organism = 'Homo_sapiens';  ## Name of the query organism
#my $organism = 'Ciona_intestinalis';  ## Name of the query organism
my @gene = ("ENSG00000139618", "ENSG00000138411");  ## List of query genes
#my @gene = ("ENSCING00000006669");  ## List of query genes
my $alltranscripts = 1;  ## the -all option (other accepted value = 1). This option is incompatible with the query list @gene (above)
my $noorf = 0;  ## Clip sequences to avoid upstream ORFs
my $from = -2000;  ## Start position of the sequence
my $to = -1;  ## End position of the sequence
my $feattype = 'mRNA';  ## The -feattype option value is  not specified, the default is used
my $type = 'upstream';  ## The -type option value; other example:'-type downstream'
my $lw = 60;  ## Line width. 0 means all on one line
my $ortho = 1;
my $taxon = 'Mammalia';
my $homology_type = 'ortholog';
my $header_organism = 'scientific';

my %args = (
            'output' => $output_choice,
            'organism' => $organism,
            'query' => \@gene,  ## An array in a hash has to be referenced (always?)
            'noorf' => $noorf,
            'from' => $from,
            'to' => $to,
            'feattype' => $feattype,
           'type' => $type,
            'line_width' => $lw,
            'ortho' => $ortho,
            'homology_type' => $homology_type,
#            'taxon' => $taxon,
            'all_transcripts' => $alltranscripts,
            'header_organism' => $header_organism
    );

## Send the request to the server
warn "Sending request to the server via SOAP::WSDL\n";

my $som = $soap->retrieve_ensembl_seq({'request' => \%args});

## Get the result
unless ($som) {
    printf "A fault (%s) occured: %s\n", $som->get_faultcode(), $som->get_faultstring();
} else {
    my $results = $som->get_response();

    ## Report the remote command
    my $command = $results -> get_command();
    print "Command used on the server: ".$command, "\n";

    ## Report the result
    my $ticket = $results -> get_server();
    print "Ticket: ".$ticket."\n";
}
