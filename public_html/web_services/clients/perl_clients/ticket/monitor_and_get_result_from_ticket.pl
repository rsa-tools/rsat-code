#!/usr/bin/env perl
# the first line of the script must tell us which language interpreter to use,
# in this case its perl

# how to test this program : 
# 1. run a client that returns a ticket (eg: matrix_scan_ticket.pl)
# 2. execute this code : monitor_and_get_result_from_ticket.pl ticket

use strict;

# import the modules we need for this test; XML::Compile is included on the server
# by default.
use XML::Compile::SOAP11;
use XML::Compile::WSDL11;
use XML::Compile::Transport::SOAPHTTP;
use Data::Dumper;

eval
{
    # Retrieving and processing the WSDL
    my $wsdl  = XML::LibXML->new->parse_file('http://pedagogix-tagc.univ-mrs.fr/rsat/web_services/RSATWS.wsdl');
    my $proxy = XML::Compile::WSDL11->new($wsdl);
    
    ##### 1.  Monitor to check if the job is finished
    
    # Generating a request message based on the WSDL
    my $client = $proxy->compileClient('monitor');
    
    #Defining a few parameters
    my $ticket = $ARGV[0];

    my %args = (
	'ticket' => $ticket
	);

    # Calling the service and getting the response
    my $answer = $client->( request => {%args});

    if ( defined $answer ) {
      print "\Status of job :", $answer->{output}->{response}->{status}, "\n";
    }
    
    ##### 2. Obtain result if job is finished
	if ($answer->{output}->{response}->{status} eq "Done") {
	
	my $client = $proxy->compileClient('get_result');
	
	# Calling the service and getting the response
    my $answer = $client->( request => {%args});
	
	if ( defined $answer ) {
     	print "\nResult :\n\n", $answer->{output}->{response}->{client}, "\n";
    }
    
    }
    
};

if ($@)
{
    print "Caught an exception\n";
    print $@."\n";
    exit 1;
}