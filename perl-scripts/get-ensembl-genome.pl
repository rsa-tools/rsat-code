#!/usr/bin/perl -w
############################################################
#
# $Id: get-ensembl-genome.pl,v 1.1 2005/02/21 16:48:56 jvanheld Exp $
#
# Time-stamp: <2003-07-04 12:48:55 jvanheld>
#
############################################################
#use strict;;
BEGIN {
    if ($0 =~ /([^(\/)]+)$/) {
	push (@INC, "$`lib/");
    }
}
require "RSA.lib";

#use DBI;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
#use Bio::EnsEMBL::DBSQL::DBConnection;
use Bio::EnsEMBL::DBSQL::SliceAdaptor;

################################################################
#### initialise parameters
my $start_time = &AlphaDate();


local %dir = ();
local %outfile = ();
$outfile{feature} = "feature.tab";

local $verbose = 0;
local $out = STDOUT;

my $host = 'ensembldb.ensembl.org';
my $user = "anonymous";
my $dbname = 'homo_sapiens_core_28_35a';

&ReadArguments();

################################################################
#### check argument values


################################################################
### open output streams

## log file
$out = &OpenOutputFile($outfile{output});

## feature file
$FT_TABLE =  &OpenOutputFile($outfile{feature});

################################################################
#### print verbose
&Verbose() if ($verbose);

################################################################
###### execute the command


## Connect to database:
my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host => $host, -user => $user, -dbname => $dbname);

my $slice_adaptor = $db->get_SliceAdaptor();

## Get one chromosome object
my $slice = $slice_adaptor->fetch_by_region('chromosome', 'X'); 
foreach my $gene (@{$slice->get_all_Genes()}) {
    my @feature = &get_feature($gene);
    print join("\t", @feature), "\n";
}


################################################################
###### print output


################################################################
###### finish verbose
if ($verbose >= 1) {
    my $done_time = &AlphaDate();
    print $out "; Job started $start_time\n";
    print $out "; Job done    $done_time\n";
}


################################################################
###### close output stream
close $out if ($outfile{output});
close $FT_TABLE if ($outfile{feature});


exit(0);


################################################################
################### subroutine definition ######################
################################################################


################################################################
#### display full help message 
sub PrintHelp {
  open HELP, "| more";
  print HELP <<End_of_help;
NAME
	get-ensembl-genome.pl

        2004 by Jacques van Helden (jvanheld\@scmbb.ulb.ac.be)
	
DESCRIPTION
	Get-Ensembl-Genome.Pl for writing new perl scripts

CATEGORY
	util

USAGE
        get-ensembl-genome.pl [-i inputfile] [-o outputfile] [-v]

OPTIONS
	-h	display full help message
	-help	display options
	-v	verbose
	-outdir output directory

   Connection to the ENSEMBL MYSQL server
   	-user	username (default: $user)
	-host	host (default: $host)
	-dbname	database name (default: $dbname)


End_of_help
  close HELP;
  exit;
}

################################################################
#### display short help message
sub PrintOptions {
  open HELP, "| more";
  print HELP <<End_short_help;
get-ensembl-genome.pl options
----------------
-h		(must be first argument) display full help message
-help		(must be first argument) display options
-i		input file
-outdir		output directory
-v		verbose
-user		username (default: $user)
-host		host (default: $host)
-dbname		database name (default: $dbname)
End_short_help
  close HELP;
  exit;
}


################################################################
#### read arguments 
sub ReadArguments {
    foreach my $a (0..$#ARGV) {
	### verbose  
	if ($ARGV[$a] eq "-v") {
	    if (&IsNatural($ARGV[$a+1])) {
		$verbose = $ARGV[$a+1];
	    } else {
		$verbose = 1;
	    }
	    
	    ### detailed help
	} elsif ($ARGV[$a] eq "-h") {
	    &PrintHelp();
	    
	    ### list of options
	} elsif ($ARGV[$a] eq "-help") {
	    &PrintOptions();
	    
	    ### output dir  
	} elsif ($ARGV[$a] eq "-outdir") {
	    $dir{output} = $ARGV[$a+1];
	    
	    ### DB user
	} elsif ($ARGV[$a] eq "-user") {
	    $user = $ARGV[$a+1];
	    
	    ### DB host
	} elsif ($ARGV[$a] eq "-host") {
	    $host = $ARGV[$a+1];
	    
	    ### DB user
	} elsif ($ARGV[$a] eq "-dbname") {
	    $dbname = $ARGV[$a+1];
	    
	}
    }
}

################################################################
#### verbose message
sub Verbose {
    print $out "; get-ensembl-genome.pl ";
    &PrintArguments($out);
    if (defined(%dir)) {
	print $out "; Directories\n";
	while (($key,$value) = each %dir) {
	    print $out ";\t$key\t$value\n";
	}
    }
    if (defined(%outfile)) {
	print $out "; Output files\n";
	while (($key,$value) = each %outfile) {
	    print $out ";\t$key\t$value\n";
	}
    }
}


################################################################
## Convert a feature to a string for export to the feature table
sub get_feature {
    my ($gene) = @_;
    my @feature = ();
   
    push @feature, $f->stable_id(); ## ID
    push @feature, "CDS"; ## We need here the feature type 
    push @feature, $f->slice->seq_region_name(); ## Chromosome name. Actually we need the contig ID
    push @feature, $f->start(); ## Start position
    push @feature, $f->end(); ## End position
    my $strand = "D";
    unless ($f->strand()) {
	$strand =  "R";
    }
    push @feature, $strand; ## Strand
    push @feature, $f->description(); ## Description
    push @feature, $f->external_name();
    
    return @feature;
}
