#!/usr/bin/perl -w
############################################################
#
# $Id: get-ensembl-genome.pl,v 1.5 2005/02/23 21:07:39 oly Exp $
#
# Time-stamp: <2003-07-04 12:48:55 jvanheld>
#
############################################################
#use strict;
BEGIN {
    if ($0 =~ /([^(\/)]+)$/) {
	push (@INC, "$`lib/");
    }
}
require "RSA.lib";

# use DBI;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
# use Bio::EnsEMBL::DBSQL::DBConnection;
use Bio::EnsEMBL::DBSQL::SliceAdaptor;

## TO DO
## Loop over all chromosomes: get the chromosome names from ENSEMBL, for a given organism
## Get the chromosome sequences 
## Get as much names as possible
## Get cross-references to various databases (Refseq, Uniprot, ...)
## Add an argument -org to specify the organism, and see how to obtain the current version
## Get intron and exon positions


################################################################
#### initialise parameters
my $start_time = &AlphaDate();


local %dir = ();

local %outfile = ();
$outfile{log}="get-ensembl-genome_log.txt";
$outfile{err}="get-ensembl-genome_err.txt";
$outfile{feature} = "feature.tab";
$outfile{xref} = "xref.tab";
$outfile{seq} = "sequence.raw";

local $verbose = 0;

## Connection to the ENSEMBL MYSQL database
my $host = 'ensembldb.ensembl.org';
my $user = "anonymous";
my $dbname = 'homo_sapiens_core_28_35a';

&ReadArguments();

################################################################
#### check argument values


################################################################
### open output streams
unless ($dir{output}) {
    $dir{output} = "ENSEMBL_".$dbname;
}
&CheckOutDir($dir{output});
chdir($dir{output});

## log file
open $log, ">".$outfile{log} || die "cannot open error log file".$outfile{log}."\n";

## error file
open ERR, ">".$outfile{err} || die "cannot open error log file".$outfile{err}."\n";

## feature file
$FT_TABLE =  &OpenOutputFile($outfile{feature});
&PrintFtHeader();

## xref file
# $XREF_TABLE = &OpenOutputFile($outfile{xref});

## sequence files???
# $SEQ = &OpenOutputFile($outfile{seq});

################################################################
#### print verbose
&Verbose() if ($verbose);

################################################################
###### execute the command


## Connect to database:
my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host => $host, -user => $user, -dbname => $dbname);
warn join ("\t", "; Db", $db), "\n" if ($main::verbose >= 3);

my $slice_adaptor = $db->get_SliceAdaptor();
warn join ("\t", "; Adaptor", $slice_adaptor), "\n" if ($main::verbose >= 3);


## TO DO: Get the list of chromosomes from the slice adaptor. 
## For the time begin I use a dirty trick
## my @chromosomes = 1..22;
## push @chromosomes, "X";

## foreach my $chromosome (@chromosomes) {
##    warn join ("\t", "; Getting chromosome", $chromosome), "\n" if ($main::verbose >= 0);
    
    ## Get one chromosome object
##    my $slice = $slice_adaptor->fetch_by_region('chromosome', $chromosome);
##    warn join ("\t", "; Slice", $slice), "\n" if ($main::verbose >= 3);
    
## Get all chromosome objects
my @slices = @{$slice_adaptor->fetch_all('chromosome')};
  
foreach $slice (@slices) {    
    
    foreach my $gene (@{$slice->get_all_Genes()}) {
	    warn join("\t", "gene", $gene), "\n" if ($main::verbose >= 5);
	    my @feature = &get_feature($gene);
	    print $FT_TABLE join("\t", @feature), "\n";
	    # print_DBEntries($gene->get_all_DBLinks());
    }
#	my $sequence = $slice->seq();
#	print $SEQ "$sequence";  ## TO DO: will need to print each sequence to a different file
}


################################################################
###### finish verbose
if ($verbose >= 1) {
    my $done_time = &AlphaDate();
    print $log "; Job started $start_time\n";
    print $log "; Job done    $done_time\n";
}


################################################################
###### close output stream
close $log if ($outfile{log});
close ERR if ($outfile{err});
close $FT_TABLE if ($outfile{feature});
# close $XREF_TABLE if ($outfile{xref});
# close $SEQ if ($outfile{seq});


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
    print $log "; get-ensembl-genome.pl ";
    &PrintArguments($log);
    if (defined(%dir)) {
	print $log "; Directories\n";
	while (($key,$value) = each %dir) {
	    print $log ";\t$key\t$value\n";
	}
    }
    if (defined(%outfile)) {
	print $log "; Output files\n";
	while (($key,$value) = each %outfile) {
	    print $log ";\t$key\t$value\n";
	}
    }
}


################################################################
## Convert a feature to a string for export to the feature table
sub get_feature {
    my ($gene) = @_;
    my @feature = ();
   
    ## ID
    my $id = $gene->stable_id(); 
    push @feature, $id;

    ## Type
    push @feature, $gene->type();


    ## Gene name
    my $name = $gene->external_name();
    unless ($name) {
	$name = $id;
	print ERR join ("\t", "No name for gene", $id, "Using ID instead"), "\n";
    }
    push @feature, $name;

    ## Chromosome name.
    push @feature, $gene->slice->seq_region_name(); 

    ## Start position
    push @feature, $gene->start();

    ## End position
    push @feature, $gene->end();

    ## Strand
    my $strand = "D";
    unless ($gene->strand()) {
	$strand =  "R";
    }
    push @feature, $strand; 

    ## Description
    my $description = $gene->description();
    if ($description) {
	push @feature, $description; 
    } else {
	push @feature, "";
	print ERR join ("\t", "No description for gene", $id), "\n";
    }
    
    return @feature;
}


################################################################
# Print cross-references !!!Problem = need to add feature id!!!
# sub print_DBEntries {
#        my $db_entries = shift;
#        foreach my $dbe (@$db_entries) {
#			print $XREF_TABLE $dbe->dbname(),"\t",$dbe->display_id(),"\n";
#        }
#}


################################################################
## Print header for the feature table
sub PrintFtHeader {
    print $FT_TABLE "-- dump date   	", &AlphaDate(), "\n";
    print $FT_TABLE "-- class       	ENSEMBL feature", "\n";
    print $FT_TABLE "-- table       	feature", "\n";
    print $FT_TABLE "-- table       	main", "\n";
    print $FT_TABLE "-- field 1	id", "\n";
    print $FT_TABLE "-- field 2	type", "\n";
    print $FT_TABLE "-- field 3	name", "\n";
    print $FT_TABLE "-- field 4	contig", "\n";
    print $FT_TABLE "-- field 5	start_pos", "\n";
    print $FT_TABLE "-- field 6	end_pos", "\n";
    print $FT_TABLE "-- field 7	strand", "\n";
    print $FT_TABLE "-- field 8	description", "\n";
    print $FT_TABLE "-- header", "\n";
    print $FT_TABLE "-- id	type	name	contig	start_pos	end_pos	strand	description", "\n";
}
