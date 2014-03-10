#!/usr/bin/perl

################################################################
## This Perl script shows how to use the Web services
## "fetch-sequences" from the Regulatory Sequence Analysis Tools
## (RSAT, http://www.rsat.eu/). 
##
## The interface relies on the SOAP/SWDL standard, which requires to
## pre-install some libraries.
##
## This script sends the URL of a bed file to the RSAT Web services,
## which run fetch-sequences in order to collect the sequences from
## UCSC. The sequences are then transferred from the RSAT Web services
## to the current client.
##
## This tester program only collects the two first sequences of the
## bed file (option -top 2), checks if the first sequence of the
## result matches the expected sequence and issue a message "Passed".

use strict;

# import the modules we need for this test; XML::Compile is included on the server
# by default.
use XML::Compile::SOAP11;
use XML::Compile::WSDL11;
use XML::Compile::Transport::SOAPHTTP;
use Data::Dumper;

eval
{
  ## Choosing the RSAT server
  my $wsdl  = XML::LibXML->new->parse_file('http://www.rsat.eu/web_services/RSATWS.wsdl');
#  my $wsdl  = XML::LibXML->new->parse_file('http://www.rsat.fr/web_services/RSATWS.wsdl');
#  my $wsdl  = XML::LibXML->new->parse_file('http://localhost/rsat/web_services/RSATWS.wsdl');
  
  ## Retriving and processing the WSDL
  my $proxy = XML::Compile::WSDL11->new($wsdl);
  
  ## Generating a request message based on the WSDL
  my $method = 'fetch_sequences';
  my $client = $proxy->compileClient($method);
  
  ## Defining the parameters for fetch-sequences
  my %args = ( 
      output => "both", ## Return both the URL of the result file on the server, and the result itself.  Alternatives: "server", "client"
      url => "http://www.rsat.eu/demo_files/fetch-sequences_Schmidt_2011_mm9_CEBPA_SWEMBL_R0.12_702peaks.bed",
      genome => "mm9",
      header_format => "", ## Header format
      upstr_ext => 0, ## Upstream extension
      downstr_ext => 0, ## Downstream extension
      extend => 0, ## Extension (alternative to the two previous arguments, since it applies downstream + upstream)
      reference => "segment",
      top => 2, ## For this test, we only fetch the two first sequences
      );
  
  # Calling the service and getting the response
  my $answer = $client->( request => {%args});

  # If the response arrived, look for a specific match
  # If the match is correct, return 0 because the test passed.
  # If the result is something else, return 2 to indicate a warning.
  # If no answer has arrived, return 1 to indicate the test failed.
  if ( defined $answer ) {
      warn ("Server command : ".$answer->{output}->{response}->{command}."\n");
      
      
      ## Note: this currently does not work. In debugging
      #    warn ("Server command : ".$answer->{output}->{response}->{server}."\n");
      #    my $response = $answer->get_response();
      #    y $server_file = $response->get_server();
      #    warn ("Result path on the server : ".$server_file."\n");
      
      if ($answer->{output}->{response}->{client} =~ 'CTGTCTATATGCCAC'){
	  print "Passed\n";
	  print "\nResult :\n\n", $answer->{output}->{response}->{client}, "\n";
	  print "\nServer :\n\n", $answer->{output}->{response}->{server}, "\n";
	  exit 0;
      } else {
	  print "Unexpected data\n";
	  print "\nResult : ", Dumper($answer), "\n";
	  exit 2;
      }
  } else {    
      print "Failed\n";
      exit 1;
  }
};

if ($@)
{
  print "Caught an exception\n";
  print $@."\n";
  exit 1;
}
