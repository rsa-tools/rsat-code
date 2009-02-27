#!/usr/bin/perl -w

use strict;
use SOAP::WSDL; ## Requires version 2.0 or later of SOAP::WSDL
use lib 'RSATWS';
use MyInterfaces::RSATWebServices::RSATWSPortType;

## WSDL location
my $server = 'http://rsat.scmbb.ulb.ac.be/rsat/web_services';

## Service call
my $soap=MyInterfaces::RSATWebServices::RSATWSPortType->new();
my $ticket = $ARGV[0];

my %args = (
            'ticket' => $ticket
    );

## Send the request to the server
print "Sending request to the server $server\n";
my $som = $soap->monitor({'request' => \%args});

## Get the result
unless ($som) {
    printf "A fault (%s) occured: %s\n", $som->get_faultcode(), $som->get_faultstring();
} else {
    my $response = $som->get_response();

    ## Report the status
    my $status = $response -> get_status();
    print "Status: ".$status."\n";
}
