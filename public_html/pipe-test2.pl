#!/usr/bin/perl -w
# pipe-test2.pl - Client retrieve-seq -> oligo-analysis -> feature-âp

################################################################
##
## This script runs a simple demo of the web service interface to the
## RSAT tools retrieve-ensembl-seq, purge-sequence and oligo-analysis linked in a workflow.
##  It sends a request to the server for discovering 6 letter words
## in upstream sequences of some genes. The sequences are first
## retrieved and purged for repeated segments
##
################################################################

## http://rsat.ulb.ac.be/rsat/pipe-test2.pl/q?taxonomy_id=83333&ensembl_gene_id=YKR034W,YJR152W,YKR039W,YGR121C,YNL142W,YPR138C,YOR348C
## http://rsat.ulb.ac.be/rsat/pipe-test2.pl/q?taxonomy_id=9606&ensembl_gene_id=ENSG00000067646

#use strict;
use SOAP::WSDL;

use CGI::Carp qw(fatalsToBrowser);
use Data::Dumper;
use File::Basename;

if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}

require "RSA.lib";
require "RSA2.cgi.lib";

## Service location
#my $server = 'http://localhost/rsat/web_services';
my $server = 'http://rsat.ulb.ac.be/rsat/web_services';
my $WSDL = $server.'/RSATWS.wsdl';
my $proxy = $server.'/RSATWS.cgi';

## Service call
my $soap=SOAP::WSDL->new(wsdl => $WSDL)->proxy($proxy);
$soap->wsdlinit;

## HTML page start
print "Content-type: text/html\n\n";
print '<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">', "\n";
print "<html><head><title>PIPE-TEST</title></head><body>\n";
#print '<H3>Contacting RSAT...</H3>';

#################################################
## Retrieve-seq part

# Parse the URL
my $query_string = $ENV{QUERY_STRING};
my @request = split (/\&/, $query_string);
my %params = ();
foreach my $param(@request){
    my @tuple = split (/=/, $param);
    $params{$tuple[0]} = $tuple[1];
}

# Construct the command
my $org = "";
my $organism = "";

# this is a trick just for tests
#if ($params{taxonomy_id} eq '9606'){
#    $organism = 'Homo_sapiens_EnsEMBL';
#    $org = 'Homo_sapiens';
#} elsif ($params{taxonomy_id} eq '83333'){
#    $organism = 'Saccharomyces_cerevisiae';
#    $org = 'Saccharomyces_cerevisiae';
#}

my @ensembl_gene_ids = split (/,/, $params{ensembl_gene_id});

if ($params{ensembl_gene_id} =~ /ENSG/) {
    $organism = 'Homo_sapiens_EnsEMBL'; ## Organism name in RSAT
    $org = 'Homo_sapiens';
} elsif ($params{ensembl_gene_id} =~ /Y/) {
    $organism = 'Saccharomyces_cerevisiae';
    $org = 'Saccharomyces_cerevisiae';
} elsif ($params{ensembl_gene_id} =~ /ENSMUSG/) {
    $organism = 'Mus_musculus_EnsEMBL';
    $org = 'Mus_musculus';
} elsif ($params{ensembl_gene_id} =~ /CG/) {
    $organism = 'Drosophila_melanogaster_EnsEMBL';
    $org = 'Drosophila_melanogaster';
} elsif ($params{ensembl_gene_id} =~ /ENSRNOG/) {
    $organism = 'Rattus_norvegicus_EnsEMBL';
    $org = 'Rattus_norvegicus';
}

## Output option
my $output_choice = 'server'; ## The result will stay in a file on the server

## Parameters
my @gene = @ensembl_gene_ids;  ## List of query genes
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
my $ensembl_host = 'xserve2.bigre.ulb.ac.be';

my %args = ('output' => $output_choice,
    'organism' => $org,
    'query' => \@gene,  ## An array in a hash has to be referenced
    'noorf' => $noorf,
    'from' => $from,
    'to' => $to,
    'feattype' => $feattype,
    'type' => $type,
    'format' => $format,
    'all' => $all,
#    'ensembl_host' => $ensembl_host
#    'label' => $label,
#    'label_sep' => $label_sep,
#    'nocom' => $nocom
);

## Send request to the server
my $som = $soap->call('retrieve_ensembl_seq' => 'request' => \%args);

## Get the result
my $server_file;  ## That variable needs to be declared outside the if..else block to be useable in the next part
if ($som->fault){  ## Report error if any
printf "A fault (%s) occured: %s<br/>", $som->faultcode, $som->faultstring;
} else {
my $results_ref = $som->result;  ## A reference to the result hash table
my %results = %$results_ref;  ## Dereference the result hash table

## Report the result file name on the server
$server_file = $results{'server'};
}

my $file_link = "http://rsat.ulb.ac.be/rsat/tmp/".basename($server_file);
print "<br/><b>Retrieved sequence(s):</b> <a href='".$file_link."'>".$file_link."</a><br/>";

my $sequence_file = $server_file;

#################################################
## Purge-sequence part

## Define hash of parameters
%args = ('output' => $output_choice,  ## Same 'server' output option
 'tmp_infile' => $server_file);  ## Output from retrieve-seq part is used as input here

## Send the request to the server
$som = $soap -> call('purge_seq' => 'request' => \%args);

## Get the result
if ($som->fault){  ## Report error if any
printf "A fault (%s) occured: %s<br/>", $som->faultcode, $som->faultstring;
} else {
my $results_ref = $som->result;  ## A reference to the result hash table
my %results = %$results_ref;  ## Dereference the result hash table

## Report the result file name on the server
$server_file = $results{'server'};
}

$file_link = "http://rsat.ulb.ac.be/rsat/tmp/".basename($server_file);
print "<br/><b>Purged sequence(s):</b> <a href='".$file_link."'>".$file_link."</a><br/>";

#################################################
## Oligo-analysis part

## Output option
#$output_choice = 'both'; ## We want to get the result on the client side, as well as the server file name

## Parameters
$verbosity = '1';
$format = 'fasta';  ## The format of input sequences
my $length = 6;  ## Length of patterns to be discovered
my $background = 'upstream-noorf';  ## Type of background used
my $stats = 'occ,proba,rank';  ## Returned statistics
my $noov = 1;  ## Do not allow overlapping patterns
my $str = 2;  ## Search on both strands
my $sort = 1;  ## Sort the result according to score
my $lth = 'occ_sig 0';  ## Lower limit to score is 0, less significant patterns are not displayed
my $pseudo = '0.05'; ## Pseudo-weight

%args = ('output' => $output_choice, 
	 'tmp_infile' => $server_file, 
	 'verbosity' => $verbosity,
	 'format' => $format,
	 'length' => $length,
	 'organism' => $organism, 
	 'background' => $background,
	 'stats' => $stats,
	 'noov' => $noov,
	 'str' => $str,
	 'sort' => $sort,
	 'lth' => $lth,
	 'pseudo' => $pseudo);

## Send request to the server
$som = $soap->call('oligo_analysis' => 'request' => \%args);

## Get the result
if ($som->fault){  ## Report error if any
    printf "A fault (%s) occured: %s<br/>", $som->faultcode, $som->faultstring;
} else {
    my $results_ref = $som->result;
    my %results = %$results_ref;
    
    ## Report the result
    $server_file = $results{'server'};
#    print $results{'client'},"<br/>";
}

$file_link = "http://rsat.ulb.ac.be/rsat/tmp/".basename($server_file);
print "<br/><b>Oligos:</b> <a href='".$file_link."'>".$file_link."</a><br/>";

###############################
# Dna-pattern

## Output option
#my $output_choice = 'both';  ## Accepted values: 'server', 'client', 'both'

my $format = 'fasta';  ## The format of input sequences
my $pattern = $server_file;  ## Pattern to be matched
my $subst;  ## Number of allowed substitutions
my $id;  ## id for the pattern
my $origin = '-0';  ## match position from end of sequence 
#my $noov = 1;  ## Do not allow overlapping patterns
my $str = 2;  ## Search on both strands
my $sort = 0;  ## Sort the result according to score when value = 1 (0 otherwise)
#my $th = 1;  ## Lower limit to score (matches count) is 1.
my $score = 8; ## Score column
my $return_field = 'sites,limits';

my %args = ('output' => $output_choice, 
	    'tmp_infile' => $sequence_file, 
	    'format' => $format,
	    'tmp_pattern_file' => $pattern,
	    'subst' => $subst, 
	    'id' => $id,
	    'origin' => $origin,
	    'score' => $score,
	    'str' => $str,
	    'sort' => $sort,
	    'return' => $return_field
#	    'th' => $th
	);

## Send request to the server
my $som = $soap->call('dna_pattern' => 'request' => \%args);

## Get the result
if ($som->fault){  ## Report error if any
    printf "A fault (%s) occured: %s<br/>", $som->faultcode, $som->faultstring;
} else {
    my $results_ref = $som->result;  ## A reference to the result hash table
    my %results = %$results_ref;  ## Dereference the result hash table

    ## Report the result
    $server_file = $results{'server'};
#    my $result = $results{'client'};
#    print "<PRE>";
#    print "<b>Matches:</b> <br/><br/>";
#    print "Pattern_ID    Strand    Pattern_sequence    Gene    Start    End    Match_sequence    Score<br/>";
#    print $result;
#    print "</PRE>";
#    print "<b>Command:</b> ".$results{'command'}."<br/>";
    $file_link = "http://rsat.ulb.ac.be/rsat/tmp/".basename($server_file);
    print "<br/><b>Feature(s):</b> <a href='".$file_link."'>".$file_link."</a><br/>";    
}

###############################                                                                         ## Convert-features

my $from = 'dnapat';
my $to = 'ft';

my %args = ('tmp_infile' => $server_file,
	    'from' => $from,
	    'to' => $to
	    );

## Send request to the server
my $som = $soap->call('convert_features' => 'request' => \%args);

## Get the result
if ($som->fault){  ## Report error if any
    printf "A fault (%s) occured: %s<br/>", $som->faultcode, $som->faultstring;
} else {
    my $results_ref = $som->result;  ## A reference to the result hash table
    my %results = %$results_ref;  ## Dereference the result hash table

    ## Report the result
    $server_file = $results{'server'};
    my $result = $results{'client'};
#    print "<b>Command:</b> ".$results{'command'}."<br/>";
    $file_link = "http://rsat.ulb.ac.be/rsat/tmp/".basename($server_file);
    print "<br/><b>Converted features:</b> <a href='".$file_link."'>".$file_link."</a><br/>";    
}

###############################
## Feature map

## Output option
$output_choice = 'both';  ## Accepted values: 'server', 'client', 'both'

my $format = 'jpg';  ## The format of output
my $tmp_infile = $server_file;  ## 
my $scorethick = 1;  ## 
my $legend = 1;  ## 
my $scalebar = 1;  ## 
my $from = 0;  ## 
my $to = 801;  ## 
my $origin = 801;  ## 
my $htmap = 1;  ## 
my $sequence_format = 'fasta';


my %args = ('output' => $output_choice, 
	    'tmp_infile' => $tmp_infile, 
	    'format' => $format, 
	    'legend' => $legend, 
	    'scorethick' => $scorethick, 
#	    'tmp_sequence_file' => $sequence_file, 
#	    'sequence_format' => $sequence_format,
	    'scalebar' => $scalebar
#	    'origin' => $origin,
#	    'from' => $from,
#	    'to' => $to, 
#	    'htmap' => $htmap
	);

## Send request to the server
my $som = $soap->call('feature_map' => 'request' => \%args);

## Get the result
if ($som->fault){  ## Report error if any
    printf "A fault (%s) occured: %s<br/>", $som->faultcode, $som->faultstring;
} else {
    my $results_ref = $som->result;  ## A reference to the result hash table
    my %results = %$results_ref;  ## Dereference the result hash table

    ## Report the result
    $server_file = $results{'server'};
    my $result = $results{'client'};
#    print "<b>Command:</b> ".$results{'command'}."<br/>";
    print "<br/><b>Feature map:</b><br/>";
    $file_link = "http://rsat.ulb.ac.be/rsat/tmp/".basename($server_file);
    print "<img src='".$file_link."'/>";
#    print "<br/><b>Feature map:</b> <a href='".$file_link."'>".$file_link."</a><br/>";    
}

print "</body></html>\n";
