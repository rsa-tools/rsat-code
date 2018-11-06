#!/usr/bin/perl

################################################################
## Connect a RSAT server and get a  list of supported organisms
##
## Usage:
##   perl supported-organisms_client_nostub.wsdl [server_URL]

use strict;

## Import the modules we need for this test.  XML::Compile enables to
## parse the WSDL (description of the web services) on the flight.
use REST::Client;
use JSON qw(decode_json encode_json);

my $client = REST::Client->new();

sub date {
    my ($sec, $min, $hour,$day,$month,$year) = localtime(time());
    my $my_date = sprintf("%02d-%02d-%02d.%02d%02d%02d", 1900+$year,$month+1,$day,$hour, $min, $sec);
    return $my_date;
}

## Specification of the server(s)
#my $server = $ARGV[0] || "http://rsat.bigre.ulb.ac.be/rsat";
my @servers =  $ARGV[0] ||
qw(http://fungi.rsat.eu
http://plants.rsat.eu
http://metazoa.rsat.eu
http://protists.rsat.eu
http://teaching.rsat.eu
http://dev.rsat.eu/
);


my @not_working = qw(
http://rsat.ulb.ac.be/rsat/
);

my @more_servers = qw(
http://rsat-tagc.univ-mrs.fr/rsat/
http://embnet.ccg.unam.mx/rsa-tools/
http://rsat.sb-roscoff.fr/
http://rsat01.biologie.ens.fr/rsa-tools/
http://floresta.eead.csic.es/rsat/
http://pedagogix-tagc.univ-mrs.fr/rsat/
http://rsat-tagc.univ-mrs.fr/rsat-dev/
http://tagc.univ-mrs.fr/rsa-tools
http://localhost/rsat
);

#        http://nexus.hgen.slu.se/rsat
#        http://anjie.bi.up.ac.za/rsa-tools
#        http://wwwsup.scmbb.ulb.ac.be/rsat

## Query parameters
my $taxon = 'Fungi';
my $return = 'ID,taxonomy';
my $format = "tab";
my $depth = '5';
my %args = (
#'return'=>$return,
'format'=>$format,
'group' => $taxon,
'depth'=>$depth,
'output'=>'email'
);


warn "Getting lists of supported organisms from server(s)\n\t";

eval
    {
        # Retrieving and processing the WSDL
        my $arg = encode_json(\%args);
        
        $client->POST('http://rsatlocal/rest/supported-organisms', $arg, {'Content-type' => 'application/json'});
        
        my $ret = $client->responseContent();
        my $decode = decode_json($ret);
        
        if($decode){
            
            my $result = $decode->{'server'};
            #my @files = split("\n", $result);
        
            warn (&date(), "\t", "Result stored in file\t", $result, "\n");
        
        ## Open output file
        # If the response arrived, look for a specific pattern
        # If the pattern is present, return 0 because the test passed.
        # If the result is something else, return 2 to indicate a warning.
        # If no answer has arrived, return 1 to indicate the test failed.
            my $answer = $decode->{'output'};
        if ( defined $answer ) {
            warn (&date(), "\t", "Server command : ".$decode->{command}."\n");
            print "; Server : ", 'http://rsatlocal/rest/', "\n";
            print "; Server command : ".$decode->{command}."\n";
            print "; Server file : ".$decode->{server}."\n";
            print $answer."\n";
        } else {
            print "No answer\n";
        }
        
        my $nb_organisms = `grep -v '^;' ${result} | grep -v '^#' | wc -l`;
        chomp($nb_organisms);
        print join("\t", $nb_organisms, "organisms at", $result), "\n";
    };
    
    ## Report exceptions
    if ($@) {
        warn (&date(), "\t", "Caught an exception\n");
        warn ($@."\n");
        print (&date(), "\t", "Caught an exception\n");
        print ($@."\n");
    }  
}
