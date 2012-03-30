#!/usr/bin/perl -w
use SOAP::Lite;
use strict;
use Getopt::Std;

my %options=();
getopts("hu:o:i:",\%options);
# like the shell getopt, "d:" means d takes an argument
die "Unprocessed by Getopt::Std:\n" if $ARGV[0];

if ($options{h}) {
  print "pathway-extractor_soapclient OPTIONS\n";
  print "===================\n";
  print "-h displays this help message and exit\n";                                                                                                                             
  print "-u server url \n";
#   print "-i seeds file \n";
  exit(0);
}
my $url = $options{u} || "http://localhost";


# defining service
my $soap = SOAP::Lite
    -> service("$url/rsat/web_services/RSATWM.wsdl");
    
# executing request    
my $results =    $soap -> supported_networks_metab();
my %resulthash = %{$results};
    
print "Supported Networks: \n" . $resulthash{"networkslist"}."\n";