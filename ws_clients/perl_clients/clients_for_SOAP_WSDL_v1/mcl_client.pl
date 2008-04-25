#!/usr/bin/perl -w

###############################################################################
##
## This script runs a simple demo of the web service interface to the
## NeAT tool MCL. It sends a request to the RSAT server 
## 1) Send a graph to MCL and clusters it
## 2) Converts the output of MCL into a tab delimited file with convert-classes
## 3) To get the help use mcl_client.pl -h
##
################################################################################

#use strict;
use SOAP::WSDL;
use Getopt::Std;
import SOAP::Lite;
# import SOAP::Lite + trace;

my %options=();
getopts("i:c:o:hv",\%options);


if (defined($options{h})) {
  print "mcl_client.pl options\n";
  print "=====================\n";
  print "-h : print this help message\n";
  print "-i : input file (tab-delimited file - required) \n";
  print "-o : output file (optional)\n"; 
  print "-c : inflation value (required)\n"; 
  print "-v : verbosity on\n"; 
  exit(0);
}
print "Unprocessed by Getopt::Std:\n" if $ARGV[0];

################################################
## Check parameters

my $input;
my $output = "mcl_clusters.tab";
my $inflation = "";
my $verbose = 0;

# check verbosity
if (exists($options{v})) {
  $verbose = 1;
}
# check inputfile
if (exists($options{i})) {
  $input = $options{i};
  if (! -e $input) {
    die("File $input does not exist\n");
  }
} else {
  die("You must specify an input file\n");
}
# check inflation
if (exists($options{c})) {
  $inflation = $options{c};
  if ($inflation !~ /^-?\d+\.?\d*$/) {
  die ("$inflation is not a valid inlation value\n");
  }
} else {
  die("You must specify an inflation value\n");
}

# check outputfile
if (!exists($options{o})) {
  $output = $input;
  $output =~ s/\.tab/_clusters.tab/;
  $output =~ s/\.csv/_clusters.tab/;
  $output =~ s/\.txt/_clusters.tab/;
  if ($output eq $input) {
    print "Rename the input file or specify an output file name\n";
    exit(0);
  }
} else {
  $output = $options{o};
}
###############################################
## Service location
my $server = 'http://rsat.scmbb.ulb.ac.be/rsat/web_services';
my $WSDL = $server.'/RSATWS.wsdl';
my $proxy = $server.'/RSATWS.cgi';
# Service call
my $soap=SOAP::WSDL->new(wsdl => $WSDL)->proxy($proxy);
$soap->wsdlinit;

if ($verbose) {
  print "Webservice properties :\n";
  print "=======================\n";
  print "Server : $server\n";
  print "WSDL File : $WSDL\n";
  print "Proxy : $proxy\n";
  print "Documentation for all RSAT / NeAT webservices can be found at http://rsat.scmbb.ulb.ac.be/rsat/web_services/RSATWS_documentation.xml\n";
}


################################################
## LAUNCH MCL
# MCL parameters
my $mcl_args = ();
# read input graph
my $input_graph = "";
open INPUTGRAPH,$input;
while (my $ligne = <INPUTGRAPH>) {
  $input_graph .= $ligne;
}

$mcl_args{inputgraph} = $input_graph;
$mcl_args{inflation} = $inflation;
$mcl_args{output_choice} = 'client';
my $output_choice = 'client';
## Send the MCL request to the server
print "Running MCL......\n" if ($verbose);
my $mcl_som = $soap->call('mcl' => 'request' => \%mcl_args);
my $mcl_result;

## Get the MCL result
if ($mcl_som->fault){ ## Report error if any
    printf "A fault (%s) occured: %s\n", $mcl_som->faultcode, $mcl_som->faultstring;
} else {
    my $mcl_results_ref = $mcl_som->result;  ## A reference to the result hash table
    my %mcl_results = %$mcl_results_ref;  ## Dereference the result hash table
    ## Report the remote command
    my $mcl_command = $mcl_results{'command'};
    print "Command used on the server: ".$mcl_command, "\n" if ($verbose);
    $mcl_result = $mcl_results{'client'};
}

################################################
## Send MCL results to convert-class
# convert-classes parameters
my $cc_args = ();
$cc_args{inputclasses} = $mcl_result;
$cc_args{informat} = 'mcl';
$cc_args{outformat} = 'tab';
$cc_args{output_choice} = 'client';


## Send the convert-classes request to the server
print "Running convert-classes.....\n" if ($verbose);
my $cc_som = $soap->call('convert_classes' => 'request' => \%cc_args);
my $cc_result;

## Get the convert-classes results
if ($cc_som->fault){ ## Report error if any
    printf "A fault (%s) occured: %s\n", $cc_som->faultcode, $cc_som->faultstring;
} else {
    my $cc_results_ref = $cc_som->result;  ## A reference to the result hash table
    my %cc_results = %$cc_results_ref;  ## Dereference the result hash table
    ## Report the remote command
    my $cc_command = $cc_results{'command'}."\n";
    print "Command used on the server: ".$cc_command, "\n" if ($verbose);
    $cc_result = $cc_results{'client'};
}

open OUTPUT, ">$output";

print OUTPUT "$cc_result";