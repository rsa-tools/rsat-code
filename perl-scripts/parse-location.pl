#!/usr/bin/perl -w
############################################################
#
# $Id: parse-location.pl,v 1.8 2009/11/05 00:32:07 jvanheld Exp $
#
# Time-stamp: <2002-06-06 14:06:08 jvanheld>
#
############################################################
#use strict;;
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}
require "RSA.lib";


#### initialise parameters ####
my $start_time = &AlphaDate;

local %infile = ();
local %outfile = ();

local $verbose = 0;

&ReadArguments;


#### check argument values ####


### open output files ###
$out = &OpenOutputFile($outfile{output});
print $out join ("\t", ";ORF", "EXON"), "\n";

##### read input #####
($in) = &OpenInputFile($infile{input});
while (<$in>) {
    chomp;
    next if (/^;/);
    next unless (/\S/);
    my @fields = split "\t";
    my $orf_id = $fields[0];
    my $location = $fields[1];
    my @exons = split ",", $location;
    my @lefts = ();
    my @rights = ();

#    print join ("\t", "HELLO", $orf_id, $#exons + 1, $location);
    foreach my $e (0..$#exons) {
	my $exon = $exons[$e];
  	my @limits = split "-", $exon;
  	if ($limits[0] < $limits[1] ) {
  	    $strand = "R";
  	} else {
  	    $strand = "D";
  	}
  	my $left = &min(@limits);
  	my $right = &max(@limits);
	$lefts[$e] = $left;
	$rights[$e] = $right;
	$exon =~ s/\-/\.\./;
	print $out join ("\t", $orf_id, $left."..".$right, 
#			$exon, $e+1, $left, $right, $strand
			), "\n";
    }

}

close $in if ($infile{input});

#### verbose ####
&Verbose if ($verbose);

###### execute the command #########


###### print output ######


###### verbose ######
if ($verbose) {
    my $done_time = &AlphaDate;
    print $out "; Job started $start_time\n";
    print $out "; Job done    $done_time\n";
}


###### close output file ######
close $out if ($outfile{output});


exit(0);


########################## subroutine definition ############################

sub PrintHelp {
#### display full help message #####
  open HELP, "| more";
  print HELP <<End_of_help;
NAME
	parse-location

        2001 by Jacques van Helden (jvanheld\@bigre.ulb.ac.be)
	
USAGE
        parse-location [-i inputfile] [-o outputfile] [-v]

DESCRIPTION
	Parse ORF locations and generates one file with introns 
	and one with exons.

CATEGORY
	parser

INPUT AND OUPUT FORMATS

	Tab-delimited text files with 2 columns. the first column
	contains the ORF identifier, and the second the ORF location.

	ORF location converts from the format

	    YDR424c	1319830-1319806,1319709-1319687,1319606-1319379

	to 
		YDR424c	1319806..1319830
		YDR424c	1319687..1319709
		YDR424c	1319379..1319606


OPTIONS
	-h	(must be first argument) display full help message
	-help	(must be first argument) display options
	-v	verbose
	-i inputfile
		if not specified, the standard input is used.
		This allows to place the command within a pipe.
	-out	output file
		if not specified, the standard output is used.
		This allows to place the command within a pipe.

End_of_help
  close HELP;
  exit;
}

sub PrintOptions {
#### display short help message #####
  open HELP, "| more";
  print HELP <<End_short_help;
parse-location options
----------------
-h	(must be first argument) display full help message
-help	(must be first argument) display options
-i	input file
-out	output file
-v	verbose
End_short_help
  close HELP;
  exit;
}


sub ReadArguments {
#### read arguments ####
    foreach my $a (0..$#ARGV) {
	### verbose ###
	if ($ARGV[$a] eq "-v") {
	    if (&IsNatural($ARGV[$a+1])) {
		$verbose = $ARGV[$a+1];
	    } else {
		$verbose = 1;
	    }
	    
	    ### detailed help
	} elsif ($ARGV[$a] eq "-h") {
	    &PrintHelp;
	    
	    ### list of options
	} elsif ($ARGV[$a] eq "-help") {
	    &PrintOptions;
	    
	    ### input file ###
	} elsif ($ARGV[$a] eq "-i") {
	    $infile{input} = $ARGV[$a+1];
	    
	    ### output file
	} elsif ($ARGV[$a] eq "-o") {
	    $outfile{output} = $ARGV[$a+1];
	    
	}
    }
}

sub Verbose {
    print $out "; parse-location ";
    &PrintArguments($out);
    if (defined(%infile)) {
	print $out "; Input files\n";
	while (($key,$value) = each %infile) {
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
