#!/usr/bin/perl -w
# retrieve_purge_oligos_client_soap-wsdl.pl - Client retrieve-seq + oligo-analysis

################################################################
##
## This script runs a simple demo of the web service interface to the
## RSAT tools retrieve-seq, purge-sequence and oligo-analysis linked in a workflow.
##  It sends a request to the server for discovering 6 letter words
## in upstream sequences of 5 yeast genes. The sequences are first
## retrieved and purged for repeated segments
##
################################################################

use strict;
use SOAP::WSDL;

warn "\nThis demo script illustrates a work flow combining three requests to the RSAT web services:\n\tretrieve-seq | purge-sequence | oligo-analysis\n\n";


## Service location
#my $server = 'http://localhost/rsat/web_services';
my $server = 'http://rsat.scmbb.ulb.ac.be/rsat/web_services';
my $WSDL = $server.'/RSATWS.wsdl';
my $proxy = $server.'/RSATWS.cgi';

## Service call
my $soap=SOAP::WSDL->new(wsdl => $WSDL)->proxy($proxy);
$soap->wsdlinit;

#################################################
## Retrieve-seq part

## Output option
my $output_choice = 'server'; ## The result will stay in a file on the server

## Parameters
my $organism = 'Saccharomyces_cerevisiae';  ## Name of the query organism
my @gene = ("PHO5", "PHO8", "PHO11", "PHO81", "PHO84");  ## List of query genes
my $all = 0;  ## -all option. This option is incompatible with the query list @gene (above)
my $noorf = 1;  ## Clip sequences to avoid upstream ORFs
my $from;  ## Start position of the sequence. Default is used (-800).
my $to;  ## End position of te sequence. Default is used (-1).
my $feattype = '';  ## -feattype option value is not defined, default is used (CDS).
my $type = '';  ## -type option value; other example:'-type downstream'
my $format = 'fasta';  ## the format of the retrieved sequence(s)
my $label = '';  ## Choice of label for the retrieved sequence(s). Default is used.
my $label_sep = '';  ## Choice of separator for the label(s) of the retrieved sequence(s). Default is used.
my $nocom = 0;  ## Other possible value = 1.

my %args = ('output' => $output_choice,
    'organism' => $organism,
    'query' => \@gene,  ## An array in a hash has to be referenced
    'noorf' => $noorf,
    'from' => $from,
    'to' => $to,
    'feattype' => $feattype,
    'type' => $type,
    'format' => $format,
    'all' => $all,
    'label' => $label,
    'label_sep' => $label_sep,
    'nocom' => $nocom);

## Send request to the server
print "\nRetrieve-seq: sending request to the server\t", $server, "\n";
my $som = $soap->call('retrieve_seq' => 'request' => \%args);

## Get the result
my $server_file;  ## That variable needs to be declared outside the if..else block to be useable in the next part
if ($som->fault){  ## Report error if any
printf "A fault (%s) occured: %s\n", $som->faultcode, $som->faultstring;
} else {
my $results_ref = $som->result;  ## A reference to the result hash table
my %results = %$results_ref;  ## Dereference the result hash table

## Report the remote command
my $command = $results{'command'};
print "Command used on the server:\n\t".$command, "\n";

## Report the result file name on the server
$server_file = $results{'server'};
print "Result file on the server:\n\t".$server_file;
}

#################################################
## Purge-sequence part

## Define hash of parameters
%args = ('output' => $output_choice,  ## Same 'server' output option
 'tmp_infile' => $server_file);  ## Output from retrieve-seq part is used as input here

## Send the request to the server
print "\nPurge-sequence: sending request to the server\t", $server, "\n";
$som = $soap -> call('purge_seq' => 'request' => \%args);

## Get the result
if ($som->fault){  ## Report error if any
printf "A fault (%s) occured: %s\n", $som->faultcode, $som->faultstring;
} else {
my $results_ref = $som->result;  ## A reference to the result hash table
my %results = %$results_ref;  ## Dereference the result hash table

## Report the remote command
my $command = $results{'command'};
print "Command used on the server: \n\t".$command, "\n";

## Report the result file name on the server
$server_file = $results{'server'};
print "Result file on the server: \n\t".$server_file;
}
#################################################
## Oligo-analysis part

## Output option
$output_choice = 'both'; ## We want to get the result on the client side, as well as the server file name

## Parameters
$format = 'fasta';  ## The format of input sequences
my $length = 6;  ## Length of patterns to be discovered
my $background = 'upstream-noorf';  ## Type of background used
my $stats = 'occ,proba,rank';  ## Returned statistics
my $noov = 1;  ## Do not allow overlapping patterns
my $str = 2;  ## Search on both strands
my $sort = 1;  ## Sort the result according to score
my @lth = ('occ_sig 0');  ## Lower limit to score is 0, less significant patterns are not displayed

%args = ('output' => $output_choice, 
	 'tmp_infile' => $server_file, 
	 'format' => $format,
	 'length' => $length,
	 'organism' => $organism, 
	 'background' => $background,
	 'stats' => $stats,
	 'noov' => $noov,
	 'str' => $str,
	 'sort' => $sort,
	 'lth' => \@lth);

## Send request to the server
print "\nOligo-analysis: sending request to the server\t", $server, "\n";
$som = $soap->call('oligo_analysis' => 'request' => \%args);

## Get the result
if ($som->fault){  ## Report error if any
    printf "A fault (%s) occured: %s\n", $som->faultcode, $som->faultstring;
} else {
    my $results_ref = $som->result;
    my %results = %$results_ref;
    
    ## Report remote commande
    my $command = $results{'command'};
    print "Command used on the server: ".$command, "\n";
    
    ## Report the result
    if ($output_choice eq 'server') {
	$server_file = $results{'server'};
	print "Result file on the server: \n\t".$server_file;
    } elsif ($output_choice eq 'client') {
	my $result = $results{'client'};
	print "Discovered oligo(s): \n".$result;
    } elsif ($output_choice eq 'both') {
	$server_file = $results{'server'};
	my $result = $results{'client'};
	print "Result file on the server: \n\t".$server_file;
	print "Discovered oligo(s): \n".$result;
    }
}
