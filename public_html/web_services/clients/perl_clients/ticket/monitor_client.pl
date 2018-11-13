#!/usr/bin/env perl -w

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
    my $client = $proxy->compileClient('monitor');

    my $ticket = $ARGV[0];

    my %args = (
	'ticket' => $ticket
	);

    # Calling the service and getting the response
    print "Sending request to the server\n";
    my $answer = $client->( request => {%args});


    if ( defined $answer ) {
      if ($answer->{output}->{response}->{status}){
	print "\nStatus : ", $answer->{output}->{response}->{status}, "\n";
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
