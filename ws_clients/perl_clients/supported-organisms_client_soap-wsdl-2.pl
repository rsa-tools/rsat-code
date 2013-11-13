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

warn "\nThis demo script gets the list of supported organisms from the remote RSAT server\n\n";

## WSDL location
##
## Actually, the server location is not affected by the variable
## $server, it depends on the stubb, which has to be configured
## separately (comment by JvH, 2012-01-30) -> I suppress this misleading message
##
#my $server = $ARGV[0] || 'http://rsat.ulb.ac.be/rsat/web_services';

## Service call
my $soap=MyInterfaces::RSATWebServices::RSATWSPortType->new();

## Output option
my $output_choice = 'both';  ## Accepted values: 'server', 'client', 'both'

## Retrieve-seq parameters
#my $taxon='Saccharomycetales';
my $taxon='Enterobacteriales';
my $return='ID,taxonomy';
my %args = (
	    'taxon' => $taxon,
	    'return'=>$return
	    );

## Send the request to the server
warn "Sending request to the server\n";
warn "Taxon\t", $taxon, "\n";
my $som = $soap->supported_organisms({'request' => \%args});

## Get the result
unless ($som) {
  die (printf "A fault (%s) occured: %s\n", $som->get_faultcode(), $som->get_faultstring());
} else {
  my $results = $som->get_response();

  ## Report the remote command
  my $command = $results -> get_command();
  warn "Command used on the server: ".$command, "\n";

  ## Report the result
  if ($output_choice eq 'server') {
    my $server_file = $results -> get_server();
    print "Result file on the server: ".$server_file;
  } elsif ($output_choice eq 'client') {
    my $result = $results -> get_client();
    warn "Result: \n";
    print $result;
  } elsif ($output_choice eq 'both') {
    my $server_file = $results -> get_server();
    my $result = $results -> get_client();
    warn "Result file on the server: ".$server_file."\n";
    warn "Result: \n";
    print $result;
  }
}
