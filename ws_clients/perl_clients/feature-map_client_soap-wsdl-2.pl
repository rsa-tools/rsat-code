#!/usr/bin/perl -w
# feature-map client

use strict;
use SOAP::WSDL;
use MIME::Base64;
use lib 'RSATWS';
use MyInterfaces::RSATWebServices::RSATWSPortType;

## Service location
my $server = 'http://rsat.scmbb.ulb.ac.be/rsat/web_services';

## Service call
my $soap=MyInterfaces::RSATWebServices::RSATWSPortType->new();

## Fixed parameters
my $sequence = "";
my $features = "PHO5	dnapat	acgtgc|gcacgt	R	438	443	ACGTGC	4.37
PHO8	dnapat	acgtgc|gcacgt	D	268	273	ACGTGC	4.37";
my $legend = 1; ## do put a legend on the feature map
my $scalebar = 1; ## do put a scalebar on the feature map
my $scalestep = 50; ## set scalesteps at 50 for scalebar
my $map_from = -800; ## start scalebar at -800
my $map_to = 0; ## end scalebar at 0
my $scorethick = 1; ## set thickness of marks on the map proportionnal to score
# my $htmap = 1;
my $map_format = 'jpg'; ## Format of feature map image

my $output_choice = 'client';

my %args = (
           'output' => $output_choice,
           'features' => $features,
#           'sequence' => $sequence,
           'legend' => $legend,
           'scalebar' => $scalebar,
           'scalestep' => $scalestep,
           'scorethick' => $scorethick,
#           'from' => $map_from,
#           'to' => $map_to,
#          'htmap' => $htmap
           'format' => $map_format
          );

  ## Send request to the server
  print "\nFeature-map: sending request to the server\t", $server, "\n";
  my $som = $soap->feature_map({'request' => \%args});


  ## Get the result
unless ($som) {
        printf "A fault (%s) occured: %s\n", $som->get_faultcode(), $som->get_faultstring();
} else {
        my $results = $som->get_response();

    ## Report the remote command
    my $command = $results -> get_command();
    print "Command used on the server: ".$command, "\n";

    ## Report the result
    my $result = $results -> get_client();
    my $output_map_file = "/no_backup/rsa-tools/ws_clients/perl_clients/feature_map_test.".$map_format;
    open OUTPUT_MAP, ">$output_map_file" or die "Can't open output file: $!\n";
    binmode OUTPUT_MAP;
    print OUTPUT_MAP decode_base64($result);
    close OUTPUT_MAP;
  }
