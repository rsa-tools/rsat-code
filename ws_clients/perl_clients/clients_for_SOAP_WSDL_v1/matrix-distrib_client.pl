#!/usr/bin/perl -w
# matrix-distrib_client.pl - Client matrix-distrib using the SOAP::WSDL module

################################################################
##
## This script runs a simple demo of the web service interface to the
## RSAT tool matrix-distrib. It sends a request to the server for
## obtaining a tab-delimited file.
##
################################################################

#use strict;
use SOAP::WSDL;
import SOAP::Lite + trace;


my %args = ();

$args{matrix_file} = "
; MET4 matrix, from Gonze et al. (2005). Bioinformatics 21, 3490-500.
A |   7   9   0   0  16   0   1   0   0  11   6   9   6   1   8
C |   5   1   4  16   0  15   0   0   0   3   5   5   0   2   0
G |   4   4   1   0   0   0  15   0  16   0   3   0   0   2   0
T |   0   2  11   0   0   1   0  16   0   2   2   2  10  11   8";

$args{matrix_pseudo} = "1";

$args{background} = "
a       a       0.3211836168662 922942  
c       c       0.1827597426890 525172  
g       g       0.1764908745757 507158  
t       t       0.3195657658692 918293  
";

$args{background_pseudo} = "0.001";
$args{matrix_format} = "tab";
$args{background_format} ="oligo-analysis"; 
$args{decimals} = "1";

my $output_choice = $args{output_choice} || 'both';

## WSDL location
my $server = 'http://rsat.scmbb.ulb.ac.be/rsat/web_services';
#my $local_server = 'http://localhost/rsa-tools/web_services';
my $WSDL = $server.'/RSATWS.wsdl';
my $proxy = $server.'/RSATWS.cgi';


## Service call
my $soap=SOAP::WSDL->new(wsdl => $WSDL)->proxy($proxy);
$soap->wsdlinit;

# $soap->wsdl_checkoccurs(0);


warn "\nThis demo script returns the theoretical distribution of matrix weights within the given background model.\n\n";

## Send the request to the server
print "Sending request to the server $server\n";
my $som = $soap->call('matrix_distrib' => 'request' => \%args);

## Get the result
if ($som->fault){ ## Report error if any
    printf "A fault (%s) occured: %s\n", $som->faultcode, $som->faultstring;
} else {
    my $results_ref = $som->result;  ## A reference to the result hash table
    my %results = %$results_ref;  ## Dereference the result hash table

    ## Report the remote command
    my $command = $results{'command'};
    print "Command used on the server: ".$command, "\n";

    ## Report the result
    if ($output_choice eq 'server') {
	my $server_file = $results{'server'};
	print "Result file on the server: ".$server_file;
    } elsif ($output_choice eq 'client') {
	my $result = $results{'client'};
	print "Matrix distrib(s): \n".$result;
    } elsif ($output_choice eq 'both') {
	my $server_file = $results{'server'};
	my $result = $results{'client'};
	print "Result file on the server: ".$server_file;
	print "Matrix distrib(s): \n".$result;
    }
}
