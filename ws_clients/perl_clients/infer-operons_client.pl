#!/usr/bin/perl

################################################################
## Run the RSATtool infer-operon remotely, using the Web servicee
## (WSDL inferface).
##
## Usage:

use strict;
use Getopt::Long qw(:config bundling); ## Required for parsing command-line arguments

# import the modules we need for this test; XML::Compile is included
# on the server by default.
use XML::Compile::SOAP11;
use XML::Compile::WSDL11;
use XML::Compile::Transport::SOAPHTTP;

package main;
{
  ## Initialize variables
  my $server = "http://rsat.ulb.ac.be/rsat";
  my $organism = "";
  my $help = "";
  my $verbose = 1;
  my $output = "client";
  my $return = "query,name,leader,operon,q_info";
  my $infile = "",
  my $distance = 55;
  my $min_gene_nb = 2;
  my $query = "";
  my $all = 0;
  my $no_args = 0;

  ## Parse arguments from the command line
  my %opt = ('v|verbose=i'=> \$verbose,
	     's|server=s'=>\$server,
	     'o|organism=s'=>\$organism,
	     'd|distance=i'=>\$distance,
	     'q|query=s'=>\$query, ## Note: alternativey, the Web services also support an array of query genes.
	     'g|min_gene_nb=i'=>\$min_gene_nb,
	     'r|return=s'=>\$return,
	     'h|help'=>\$help,
	    );

  ## Print command-line arguments for the sake of reproducibility
  print "; Client command: ", join(" ", $0, @ARGV), "\n";

  ## Check arguments before running &GetOptions() because
  ## &GetOptions() apparently creates values in @ARGV
  if (scalar(@ARGV) < 1) {
    $no_args = 1;
  }

  ## Automatically parse the options with the Getopt module
  &GetOptions(%opt);

  ## Print help message
  if (($no_args) || ($help)) {
    &PrintHelp();
    exit(0) if ($help);
  }


  ## Organism is a mandatory argument
  unless ($organism) {
    die "\n\nError: Organism should be specified (option -o)\n\n";
  }

  ## If no query is specified, compute all operons
  unless ($query) {
    $all = 1;
    warn "No query specified, inferring operons for all genes.\n";
  }

  ## Specification of the server.
  unless ($server) {
    $server = "http://rsat.ulb.ac.be/rsat";
  }

  ## Query parameters
  my %args = (
	      'output' => "both",
	      'organism'=>$organism,
	      'query'=>$query,
	      'tmp_infile'=>$infile,
	      'all'=>$all,
	      'distance'=>$distance,
	      'min_gene_nb'=>$min_gene_nb,
	      'return'=>$return,
	     );

  warn "DEBUG\targs\t", join "; ", %args, "\n" if ($verbose >= 3);

  eval {
    # Retrieving and processing the WSDL
    my $wsdl_url = $server.'/web_services/RSATWS.wsdl';
    warn ("Parsing Web service description from WSDL", "\t", $wsdl_url, "\n") if ($verbose >= 2);
    my $wsdl  = XML::LibXML->new->parse_file($wsdl_url);
    my $proxy = XML::Compile::WSDL11->new($wsdl);

    ## Compiling the client for infer-operon
    warn ("Compiling client\n") if ($verbose >= 2);
    my $client = $proxy->compileClient('infer_operon');

    # Calling the service and getting the response
    warn ("Sending query to server", "\t", $server, "\n") if ($verbose >= 2);
    my $answer = $client->( request => {%args});


    ## Analyze the answer and print out the result
    if ( defined $answer ) {
      print "; Server : ", $server, "\n";
      print "; WSDL URL : ", $wsdl_url, "\n";
      print "; Server command : ".$answer->{output}->{response}->{command}."\n";
      print "; Server file : ".$answer->{output}->{response}->{server}."\n";
      print $answer->{output}->{response}->{client};

    } else {
      print "No answer\n";
    }

    ## Catch exceptions issued by the WS server
    if ($@) {
      warn "Caught an exception\n";
      warn $@."\n";
      print OUT "Caught an exception\n";
      print OUT $@."\n";
    }
  }
}

sub PrintHelp {
  print <<EndHelp;
NAME

infer-operons_client.pl

DESCRIPTION

This script runs the RSAT command infer-operons on a remote server,
using the Web Services interface.

AUTHOR

Jacques.van-Helden\@univ-amu.fr

USAGE

    perl infer-operons_client.pl -o organism [-q query] \
        [-v verbose_level] [-s server] [-d distance] [-g min_gene_number] \
        [-return return_fields]

EXAMPLE
    perl infer-operons_client.pl -o Escherichia_coli_K12 -q lysA

ARGUMENTS

  Mandatory arguments

    -o, --organism  organism_name

       The list of organisms supported on the server can be obtained with
       the script supported-organisms_client.pl.


  Optional argument
    -v, --verbose #
       Verbosity level.

    -s, --server server_url
        URL of the server (default: $main::server)

    -d, --distance
        Splitting distance between two genes: if intergenic regions is
        larger than this value, they two genes are considered to
        belong to distinct operons.

    -g, --min_gene_nb #
         Mininal number of genes (operons with less than g genes are
         not reported).

    -q, --query
        Query gene. 

        If no query is specified, infer-operon returns the predicted
        operons for all the genes in the genome.

    -r, --return return_fields
        List of fields to return, separated by commas.


    -i input_file
        Name of an input file contianing the queries (one query per
        row).

EndHelp
}
