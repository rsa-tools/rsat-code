#!/usr/bin/perl -w
############################################################
#
# $Id: get-ensembl-genome.pl,v 1.8 2005/03/10 14:30:37 jvanheld Exp $
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
## Add an argument -org to specify the organism, and see how to obtain the current version
## (not sure it is possible; dbname seems mandatory)
## Get intron positions ? (we have transcripts and exons)

################################################################
#### initialise parameters
my $start_time = &AlphaDate();
local $slice_type = "chromosome";
local $noseq=0;

local %dir = ();

local %outfile = ();
$outfile{log}="get-ensembl-genome_log.txt";
$outfile{err}="get-ensembl-genome_err.txt";
$outfile{feature} = "feature.tab";
$outfile{xreference} = "xreference.tab";
$outfile{ft_name} = "feature_name.tab";

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
open $XREF_TABLE, ">".$outfile{xreference} || die "cannot open error log file".$outfile{xreference}."\n";

## feature_name table
open $FT_NAME_TABLE, ">".$outfile{ft_name} || die "cannot open error log file".$outfile{ft_name}."\n";


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

## Get the list of  slices (chromosomes by default)
warn join ("\t", "; Getting the list of slices of type", $slice_type), "\n" if ($main::verbose >= 1);
my @slices = @{$slice_adaptor->fetch_all($slice_type)};

warn join ("\n;\t", "; Number of slices", scalar(@slices)), "\n" if ($main::verbose >= 1);

foreach my $slice (@slices) {
    my $slicename = $slice->seq_region_name();
    warn join ("\t", "; getting features for slice", $slicename), "\n" if ($main::verbose >= 2);
    
    ## Get all Gene objects
    foreach my $gene (@{$slice->get_all_Genes()}) {
        warn join("\t", "gene", $gene), "\n" if ($main::verbose >= 5);
        @feature = &get_feature($gene);
        print $FT_TABLE join("\t", @feature), "\n";
        print_DBEntries($gene->get_all_DBLinks());
	
        ## Get all Transcript objects
        foreach my $trans (@{$gene->get_all_Transcripts()}) {
            warn join("\t", "transcript", $trans), "\n" if ($main::verbose >= 5);
            my @feature = &get_feature($trans);
            $feature[1] = "transcript";
            print $FT_TABLE join("\t", @feature), "\n";
            $transcriptID = $feature[0];
	    
	    ## Get CDS ID and coordinates (relative to chromosome) - there is a strand trick (see API doc)
	    ## Problem: corrdinates are strange
	    my $coding_region_start = $trans->coding_region_start;
	    my $coding_region_end = $trans->coding_region_end;
	    if($trans->translation()) {
		$feature[0] =  $trans->translation()->stable_id();
		$feature[1] = "CDS";
		if ($feature[6] eq 'D') {
		    $feature[4] = $coding_region_start;
		    $feature[5] = $coding_region_end;
		} else {
		    $feature[4] = $coding_region_end;
		    $feature[5] = $coding_region_start;
		}
		print $FT_TABLE join ("\t", @feature), "\n"; 
	    }
	    
## Tests to make sure a transcripts includes UTRs (are also in first and last exons!)
#        print $feature[0], " : ", $trans->spliced_seq(), "\n";
#        print $feature[0], " : ", $trans->translateable_seq(), "\n";
#        my $fiv_utr = $trans->five_prime_utr();
#        my $thr_utr = $trans->three_prime_utr();
#        print $feature[0], " : ", ($fiv_utr) ? $fiv_utr->seq() : 'No 5 prime UTR', "\n";
#        print $feature[0], " : ", ($thr_utr) ? $thr_utr->seq() : 'No 3 prime UTR', "\n";

            ## Get all Exon objects
            foreach my $exon (@{$trans->get_all_Exons()}) {
                my @exonfeature = &get_exonfeature($exon);
                print $FT_TABLE join("\t", @exonfeature), "\n";
            }
	    
            ## Get all Intron objects
            foreach my $intron (@{$trans->get_all_Introns()}) {
                my @intronfeature = &get_intronfeature($intron);
                $intronfeature[0] = "Trnscrpt - ".$transcriptID;
                print $FT_TABLE join("\t", @intronfeature), "\n";
            }
        }

    }
    
    ################################################################
    ## Export the raw sequence in a separate file.
    unless ($noseq) {
	my $sequence = $slice->seq();
	$outfile{sequence} = "$feature[3].raw";
	open $SEQ, ">".$outfile{sequence} || die "cannot open error log file".$outfile{sequence}."\n";
	print $SEQ $sequence;
	close $SEQ if ($outfile{sequence});
    }
}

## foreach my $chromosome (@chromosomes) {
##    warn join ("\t", "; Getting chromosome", $chromosome), "\n" if ($main::verbose >= 0);
    
    ## Get one chromosome object
##    my $slice = $slice_adaptor->fetch_by_region('chromosome', $chromosome);
##    warn join ("\t", "; Slice", $slice), "\n" if ($main::verbose >= 3);
    
# ## Get all chromosome objects
# my @slices = @{$slice_adaptor->fetch_all('chromosome')};
  
# foreach my $slice (@slices) {    
#     foreach my $gene (@{$slice->get_all_Genes()}) {
# 	    warn join("\t", "gene", $gene), "\n" if ($main::verbose >= 5);
# 	    @feature = &get_feature($gene);
# 	    print $FT_TABLE join("\t", @feature), "\n";
# 	    print_DBEntries($gene->get_all_DBLinks());
#     }
    
#  	my $sequence = $slice->seq();
#  	$outfile{sequence} = "$feature[3].raw";
#     open $SEQ, ">".$outfile{sequence} || die "cannot open error log file".$outfile{sequence}."\n";
#  	print $SEQ $sequence;
#     close $SEQ if ($outfile{sequence});

# }


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
close $XREF_TABLE if ($outfile{xreference});
close $FT_NAME_TABLE if ($outfile{ft_name});


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
	-noseq  do not export the sequence (only the features)

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
-v		verbose
-i		input file
-outdir		output directory
-noseq          do not export the sequence (only the features)
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
	    
	    ### do not import the sequence
	} elsif ($ARGV[$a] eq "-noseq") {
	    $noseq=1;
	    
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
## Convert an exon feature to a string for export to the feature table
sub get_exonfeature {
    my ($exon) = @_;
    my @exonfeature = ();
 
    ## ID
    my $id = $exon->stable_id();
    push @exonfeature, $id;

    ## Type
    push @exonfeature, "exon";
   
    ## Exon name
    push @exonfeature, "";

    ## Chromosome name.
    push @exonfeature, $exon->slice->seq_region_name();

    ## Start position
    push @exonfeature, $exon->start();

    ## End position
    push @exonfeature, $exon->end();
   
    ## Strand
    my $strand = "D";
    unless ($exon->strand()) {
        $strand =  "R";
    }
    push @exonfeature, $strand;

    ## Description
    push @exonfeature, "";
  
    return @exonfeature;
}


################################################################
## Convert an intron feature to a string for export to the feature table
sub get_intronfeature {
    my ($intron) = @_;
    my @intronfeature = ();

    ## ID (introns have no ID in EnsEMBL)
    push @intronfeature, "";

    ## Type
    push @intronfeature, "intron";

    ## Intron name
    push @intronfeature, "";

    ## Chromosome name.
    push @intronfeature, $intron->slice->seq_region_name();

    ## Start position
    push @intronfeature, $intron->start();

    ## End position
    push @intronfeature, $intron->end();

    ## Strand
    my $strand = "D";
    unless ($intron->strand() == 1) {
        $strand =  "R";
    }
    push @intronfeature, $strand;

    ## Description
    push @intronfeature, "";

    return @intronfeature;
}



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



# ################################################################
# ## Print header for the feature table
# sub PrintFtHeader {
#     print $FT_TABLE "-- dump date   	", &AlphaDate(), "\n";
#     print $FT_TABLE "-- class       	ENSEMBL feature", "\n";
#     print $FT_TABLE "-- table       	feature", "\n";
#     print $FT_TABLE "-- table       	main", "\n";
#     print $FT_TABLE "-- field 1	id", "\n";
#     print $FT_TABLE "-- field 2	type", "\n";
#     print $FT_TABLE "-- field 3	name", "\n";
#     print $FT_TABLE "-- field 4	contig", "\n";
#     print $FT_TABLE "-- field 5	start_pos", "\n";
#     print $FT_TABLE "-- field 6	end_pos", "\n";
#     print $FT_TABLE "-- field 7	strand", "\n";
#     print $FT_TABLE "-- field 8	description", "\n";
#     print $FT_TABLE "-- header", "\n";
#     print $FT_TABLE "-- id	type	name	contig	start_pos	end_pos	strand	description", "\n";
# }

################################################################
# Print cross-references and names to files
sub print_DBEntries {
    my $db_entries = shift;
    foreach my $dbe (@$db_entries) {
        print $XREF_TABLE $feature[0],"\t",$dbe->dbname(),"\t",$dbe->display_id(),"\n";
        if ($dbe->dbname() eq 'HUGO') {
            print $FT_NAME_TABLE $feature[0],"\t",$dbe->display_id(),"\n";
        }
    }
}
