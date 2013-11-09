#!/usr/bin/perl
################################################################
#### use the soap interface to retrieve a list of oranisms supported in KEGG
###INSTAL SOAP!

use strict;
use SOAP::Lite;
my $def='http://soap.genome.ad.jp/KEGG.wsdl';
my $response = SOAP::Lite->service("$def")->list_organisms;

warn "I got a response\n", $response, "\n";

foreach(@{$response}) {
    my $orgs = SOAP::Lite->service("$def")->dbinfo($_);
    print $orgs, "\n";
}
