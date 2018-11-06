#!/usr/bin/perl

################################################################
## Connect a RSAT server and get a  list of supported organisms

use strict;

## Import the modules we need for this test.
use REST::Client;
use JSON qw(decode_json encode_json);

my $client = REST::Client->new();

sub date {
    my ($sec, $min, $hour,$day,$month,$year) = localtime(time());
    my $my_date = sprintf("%02d-%02d-%02d.%02d%02d%02d", 1900+$year,$month+1,$day,$hour, $min, $sec);
    return $my_date;
}

## Query parameters
my $taxon = 'Fungi';
my $format = "tab";
my $depth = '5';
my %args = (
#'return'=>$return,
'format'=>$format,
'group' => $taxon,
'depth'=>$depth,
'output'=>'email'
);


#warn "Getting lists of supported organisms from server(s)\n\t";

eval
    {
        # send REST request
        my $arg = encode_json(\%args);
        $client->POST('http://rsatlocal/rest/supported-organisms', $arg, {'Content-type' => 'application/json'});
        
        # get response
        my $ret = $client->responseContent();
        my $decode = decode_json($ret);
        
        # process response
        if($decode){
            my $file = $decode->{'server'};
            warn (&date(), "\t", "Result stored in file\t", $file, "\n");
        
        ## Open output file
        
            my $answer = $decode->{'output'};
        if ( defined $answer ) {
            warn (&date(), "\t", "Server command : ".$decode->{command}."\n");
            print "; Server : ", 'http://rsatlocal/rest/', "\n";
            print "; Server command : ".$decode->{command}."\n";
            print "; Server file : ".$decode->{server}."\n";
            print "; RESULT :\n". $answer."\n";
        } else {
            print "No answer\n";
        }
        
        my $nb_organisms = `grep -v '^;' ${file} | grep -v '^#' | wc -l`;
        chomp($nb_organisms);
        print join("\t", $nb_organisms, "organisms at", $file), "\n";
    };
    
    ## Report exceptions
    if ($@) {
        warn (&date(), "\t", "Caught an exception\n");
        warn ($@."\n");
        print (&date(), "\t", "Caught an exception\n");
        print ($@."\n");
    }  
}
