#!/usr/bin/env perl

use strict;
use XML::Compile::SOAP11;
use XML::Compile::WSDL11;
use XML::Compile::Transport::SOAPHTTP;
use Data::Dumper;

eval
{
    # Retriving and processing the WSDL
    my $wsdl  = XML::LibXML->new->parse_file('http://pedagogix-tagc.univ-mrs.fr/rsat/web_services/RSATWS.wsdl');
    my $proxy = XML::Compile::WSDL11->new($wsdl);
    
    # Generating a request message based on the WSDL
    my $client = $proxy->compileClient('get_result');

    my $ticket = $ARGV[0];

    my %args = (
	'ticket' => $ticket
	);

    # Calling the service and getting the response
    print "Sending request to the server\n";
    my $answer = $client->( request => {%args});

    if ( defined $answer ) {
      if ($answer->{output}->{response}->{client}){
#      if ($answer->{output}->{response}->{server}){
	print "\nResult :\n", $answer->{output}->{response}->{client}, "\n";
#	print "\nResult :\n", $answer->{output}->{response}->{server}, "\n";
	exit 0;
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
