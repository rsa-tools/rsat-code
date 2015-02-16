#!/usr/bin/perl -w
# workflow_nature_protocols.pl

################################################################
##
## This script runs a workflow of the web service interface to the
## RSAT tools retrieve-seq, purge-sequence, oligo-analysis, dna-pattern,
## convert-features and feature-map.
##  It sends requests to the server for discovering 6 letter words
## in upstream sequences of clusters of yeast genes. The sequences are first
## retrieved and purged for repeated segments. The discovered motifs as well
## as a map showing their positions on the sequences are printed to files
##
################################################################

use strict;
use SOAP::WSDL;

warn "\nThis demo script illustrates a work flow combining several requests to the RSAT web services:\n\tretrieve-seq | purge-sequence | oligo-analysis | dna_pattern | convert_features | feature_map\n\n";


## Service location
#my $server = 'http://localhost/rsat/web_services';
my $server = 'http://rsat.scmbb.ulb.ac.be/rsat/web_services';
my $WSDL = $server.'/RSATWS.wsdl';
my $proxy = $server.'/RSATWS.cgi';

## Service call
my $soap=SOAP::WSDL->new(wsdl => $WSDL)->proxy($proxy);
$soap->wsdlinit;

## Fixed parameters
my $organism = 'Saccharomyces_cerevisiae';  ## Name of the query organism
my $noorf = 1;  ## Clip sequences to avoid upstream ORFs
my $background = 'upstream-noorf';  ## Type of background used
my $length = 6;  ## Length of patterns to be discovered
my $stats = 'occ,proba,rank';  ## Returned statistics
my $noov = 1;  ## Do not allow overlapping patterns
my $str = 2;  ## Search on both strands
my $sort = 1;  ## Sort the result according to score
my $lth = 'occ_sig 0';  ## Lower limit to score is 0, less significant patterns are not displayed
my $origin = '-0';  ## match position from end of sequence
my $score_column = 8;
my $th = 1;  ## Lower limit to score (matches count) is 1.
my $feature_in = 'dnapat';
my $feature_out = 'ft';
my $legend = 1;
my $scalebar = 1;
my $scalestep = 50;
my $map_from = -800;
my $map_to = 0;
my $scorethick = 1;
# my $htmap = 1;
my $map_format = 'jpg';

## Read cluster file
open INPUT, "/Users/oly/Desktop/Nature_protocols/data/Harbison_2004_sig0.fam" or die "Can't open input file: $!\n";
# open INPUT, "/Users/oly/Desktop/Nature_protocols/data/Harbison_2004_sig0_select_sig4.fam" or die "Can't open input file: $!\n";

## Loop over gene clusters
my @genes;
my $cluster = "couille";
while (my $line = <INPUT>) {
  chomp $line;
  next if $line =~ /ID/;
  my @entry = split /\t/, $line;

  ## Put gene list in array
  if ($entry[1] eq $cluster) {
    push @genes, $entry[0];
#    print ("\tCluster= ", $cluster,"\n\tGenes= ", "@genes", "\n");
    next;
  } else {
    if (@genes) {
      &Analysis();
    }
    $cluster = $entry[1];
    @genes='';
    push @genes, $entry[0];
    next;
  }
}
&Analysis();  ## run with last cluster

close INPUT;

###################################################
sub Analysis {

#  die ("TEST FINISHED\n");
  #################################################
  ## Retrieve-seq part

  ## Output option
  my $output_choice = 'server'; ## The result will stay in a file on the server

  my %args = (
	      'output' => $output_choice,
	      'organism' => $organism,
	      'query' => \@genes,  ## An array in a hash has to be referenced
	      'noorf' => $noorf
	     );

  ## Send request to the server
  print "\nRetrieve-seq: sending request to the server\t", $server, "\n";
  my $som = $soap->call('retrieve_seq' => 'request' => \%args);

  ## Get the result
  my $sequence_file;  ## That variable needs to be declared outside the if..else block to be useable in the next part
  if ($som->fault){  ## Report error if any
    printf "A fault (%s) occured: %s\n", $som->faultcode, $som->faultstring;
  } else {
    my $results_ref = $som->result;  ## A reference to the result hash table
    my %results = %$results_ref;  ## Dereference the result hash table

    ## Report the remote command
    my $command = $results{'command'};
    print "Command used on the server:\n\t".$command, "\n";

    ## Report the result file name on the server
    $sequence_file = $results{'server'};
    print "Result file on the server:\n\t".$sequence_file;
  }
  #################################################
  ## Purge-sequence part

  %args = (
	   'output' => $output_choice,  ## Same 'server' output option
	   'tmp_infile' => $sequence_file  ## Output from retrieve-seq part is used as input here
	  );

  ## Send the request to the server
  print "\nPurge-sequence: sending request to the server\t", $server, "\n";
  $som = $soap -> call('purge_seq' => 'request' => \%args);

  ## Get the result
  my $purged_sequence_file;  ## That variable needs to be declared outside the if..else block to be useable in the next part
  if ($som->fault){  ## Report error if any
    printf "A fault (%s) occured: %s\n", $som->faultcode, $som->faultstring;
  } else {
    my $results_ref = $som->result;  ## A reference to the result hash table
	my %results = %$results_ref;  ## Dereference the result hash table

    ## Report the remote command
    my $command = $results{'command'};
    print "Command used on the server: \n\t".$command, "\n";

    ## Report the result file name on the server
    $purged_sequence_file = $results{'server'};
    print "Result file on the server: \n\t".$purged_sequence_file;
  }
  #################################################
  ## Oligo-analysis part

  $output_choice = 'both';

  %args = (
	   'output' => $output_choice,
	   'tmp_infile' => $purged_sequence_file,
	   'length' => $length,
	   'organism' => $organism,
	   'background' => $background,
	   'stats' => $stats,
	   'noov' => $noov,
	   'str' => $str,
	   'sort' => $sort,
	   'lth' => $lth,
	   'verbosity' => 1
	  );

  ## Send request to the server
  print "\nOligo-analysis: sending request to the server\t", $server, "\n";
  $som = $soap->call('oligo_analysis' => 'request' => \%args);

  ## Get the result
  my $oligos_file;  ## That variable needs to be declared outside the if..else block to be useable in the next part
  if ($som->fault){  ## Report error if any
    printf "A fault (%s) occured: %s\n", $som->faultcode, $som->faultstring;
  } else {
    my $results_ref = $som->result;
    my %results = %$results_ref;

    ## Report remote command
    my $command = $results{'command'};
    print "Command used on the server: \n\t".$command, "\n";

    ## Report the result file name on the server
    $oligos_file = $results{'server'};
    print "Result file on the server: \n\t".$oligos_file;

    ## Store the result in file
    my $result = $results{'client'};
    my $output_oligo_file = "/Users/oly/Desktop/Nature_protocols/results/sig0/".$cluster."_oligos.tab";
#    my $output_oligo_file = "/Users/oly/Desktop/Nature_protocols/results/sig4/".$cluster."_oligos.tab";
    open OUTPUT_OLIGO, ">$output_oligo_file" or die "Can't open input file: $!\n";
    print OUTPUT_OLIGO $result;
    close OUTPUT_OLIGO;
  }

  #################################################
  ## Dna-pattern part

  $output_choice = 'server';

  %args = ('output' => $output_choice,
	   'tmp_infile' => $sequence_file,
	   'tmp_pattern_file' => $oligos_file,
	   'origin' => $origin,
	   #	 'noov' => $noov,
	   'str' => $str,
	   'sort' => $sort,
	   'th' => $th,
	   'return' => 'sites,limits',
	   'score' => $score_column
	  );

  ## Send request to the server
  print "\nDna-pattern: sending request to the server\t", $server, "\n";
  $som = $soap->call('dna_pattern' => 'request' => \%args);

  ## Get the result
  my $feature_file;
  if ($som->fault){  ## Report error if any
    printf "A fault (%s) occured: %s\n", $som->faultcode, $som->faultstring;
  } else {
    my $results_ref = $som->result;
    my %results = %$results_ref;

    ## Report remote command
    my $command = $results{'command'};
    print "Command used on the server: \n\t".$command, "\n";

    ## Report the result file name on the server
    $feature_file = $results{'server'};
    print "Result file on the server: \n\t".$feature_file, "\n";

    ## Report the result
#    my $result = $results{'client'};
#    print "Features(s): \n\t".$result;
  }

  #################################################
  ## Convert-features part

  %args = ('output' => $output_choice,
	   'tmp_infile' => $feature_file,
	   'from' => $feature_in,
	   'to' => $feature_out
	  );

  ## Send request to the server
  print "\nConvert-features: sending request to the server\t", $server, "\n";
  $som = $soap->call('convert_features' => 'request' => \%args);

  ## Get the result
  my $converted_feature_file;
  if ($som->fault){  ## Report error if any
    printf "A fault (%s) occured: %s\n", $som->faultcode, $som->faultstring;
  } else {
    my $results_ref = $som->result;
    my %results = %$results_ref;

    ## Report remote command
    my $command = $results{'command'};
    print "Command used on the server: \n\t".$command, "\n";

    ## Report the result file name on the server
    $converted_feature_file = $results{'server'};
    print "Result file on the server: \n\t".$converted_feature_file, "\n";

    ## Report the result
    #  my $result = $results{'client'};
    #  print "Converted features(s): \n\t".$result;
  }

  #################################################
  ## Feature-map part
  chomp $sequence_file;

  $output_choice = 'both';

  %args = (
	   'output' => $output_choice,
	   'tmp_infile' => $converted_feature_file,
	   'tmp_sequence_file' => $sequence_file,
	   'legend' => $legend,
	   'scalebar' => $scalebar,
	   'scalestep' => $scalestep,
	   'scorethick' => $scorethick,
	   'from' => $map_from,
	   'to' => $map_to,
	   #	 'htmap' => $htmap
	   'format' => $map_format
	  );

  ## Send request to the server
  print "\nFeature-map: sending request to the server\t", $server, "\n";
  $som = $soap->call('feature_map' => 'request' => \%args);

  ## Get the result
  my $feature_map_file;
  if ($som->fault){  ## Report error if any
    printf "A fault (%s) occured: %s\n", $som->faultcode, $som->faultstring;
  } else {
    my $results_ref = $som->result;
    my %results = %$results_ref;

    ## Report remote command
    my $command = $results{'command'};
    print "Command used on the server: \n\t".$command, "\n";

    ## Report the result file name on the server
    $feature_map_file = $results{'server'};
    print "Result file on the server: \n\t".$feature_map_file, "\n";

    ## Store the result in file
    my $result = $results{'client'};
    my $output_map_file = "/Users/oly/Desktop/Nature_protocols/results/sig0/".$cluster."_map.".$map_format;
#    my $output_map_file = "/Users/oly/Desktop/Nature_protocols/results/sig4/".$cluster."_map.".$map_format;
    open OUTPUT_MAP, ">$output_map_file" or die "Can't open input file: $!\n";
    print OUTPUT_MAP $result;
    close OUTPUT_MAP;
  }
}
