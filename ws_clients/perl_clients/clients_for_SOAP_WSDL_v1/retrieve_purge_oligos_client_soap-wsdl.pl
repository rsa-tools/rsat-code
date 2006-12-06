#!/usr/bin/perl -w
# retrieve_purge_oligos_client_wsdl.pl - Client retrieve-seq + oligo-analysis

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
#use SOAP::Lite +trace;

## Service location
my $server = 'http://localhost/rsat/web_services';
#my $server = 'http://rsat.scmbb.ulb.ac.be/rsat/web_services';
my $WSDL = $server.'/RSATWS.wsdl';
my $proxy = $server.'/RSATWS.cgi';

my $soap=SOAP::WSDL->new(wsdl => $WSDL)->proxy($proxy);

$soap->wsdlinit;

#################################################
## Retrieve-seq part

## Return option
my $return_choice = 'file';

## Parameters
my $organism = '-org Saccharomyces_cerevisiae';  ## Name of the query organism
my @gene = ("PHO5", "PHO8", "PHO11", "PHO81", "PHO84");  ## List of query genes
my $all = '';  ## -all option. This option is incompatible with the query list @gene (above)
my $noorf = '-noorf';  ## Clip sequences to avoid pstream ORFs
my $from;  ## Start position of the sequence
my $to;  ## End position of te sequence
my $feattype = '';  ## -feattype option value is not defined, default is used
my $type = '';  ## -type option value; other example:'-type downstream'
my $format = '-format fasta';  ## the format of the retrieved sequence(s)
my $label = '';  ## Choice of label for the retrieved sequence(s)
my $label_sep = '';  ## Choice of separator for the label(s) of the retrieved sequence(s)
my $nocom = '';  ## Other possible value = '-nocom'

my %args = ('return' => $return_choice,
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
my $server_file;  ## That variable needs to be declared outside the if..else block
if ($som->fault){  ## Report error if any
    printf "A fault (%s) occured: %s\n", $som->faultcode, $som->faultstring;
} else {
    my $results_ref = $som->result;  ## A reference to the result hash table
    my %results = %$results_ref;  ## Dereference the result hash table
    
    ## Report the remote command
    my $command = $results{'command'};
    print "Command used on the server:\n\t".$command, "\n";
    
    ## Report the result file name on the server
    $server_file = $results{'file'};
    print "Result file on the server:\n\t".$server_file;
}

#################################################
## Purge-sequence part

## Define hash of parameters
%args = ('return' => $return_choice,  ## Same 'file' return option
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
    $server_file = $results{'file'};
    print "Result file on the server: \n\t".$server_file;
}
#################################################
## Oligo-analysis part

## Return option
$return_choice = 'both';
## Parameters
$format = '-format fasta';  ## The format of input sequences
my $length = 6;  ## Length of patterns to be discovered
my $background = '-bg upstream-noorf';  ## Type of background used
my $stats = '-return occ,proba,rank';  ## Returned statistics
my $noov = '-noov';  ## Do not allow overlapping patterns
my $str = '-2str';  ## Search on both strands
my $sort = '-sort';  ## Sort the result according to score
my $lth = '-lth occ_sig 0';  ## Lower limit to score is 0, less significant patterns are not displayed

%args = ('return' => $return_choice, 
	 'tmp_infile' => $server_file, 
	 'format' => $format,
	 'length' => $length,
	 'organism' => $organism, 
	 'background' => $background,
	 'stats' => $stats,
	 'noov' => $noov,
	 'str' => $str,
	 'sort' => $sort,
	 'lth' => $lth);

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
    if ($return_choice eq 'file') {
	$server_file = $results{'file'};
	print "Result file on the server: \n\t".$server_file;
    } elsif ($return_choice eq 'result') {
	my $result = $results{'result'};
	print "Discovered oligo(s): \n".$result;
    } elsif ($return_choice eq 'both') {
	$server_file = $results{'file'};
	my $result = $results{'result'};
	print "Result file on the server: \n\t".$server_file;
	print "Discovered oligo(s): \n".$result;
    }
}
