#!/usr/bin/perl -w
############################################################
#
# $Id: calibrate-oligos.pl,v 1.1 2003/12/23 02:55:40 jvanheld Exp $
#
# Time-stamp: <2003-07-04 12:48:55 jvanheld>
#
############################################################
#use strict;;
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}
require "RSA.lib";

################################################################
#### initialise parameters
my $start_time = &AlphaDate;
my $oligo_length = 6;
my $repet = 1000;
my $oligo_len = 6;
my $str = "-1str";
my $seq_len=800;
my $seq_nb=10;
my $ov="-noov";
my $organism_name = "Saccharomyces_cerevisiae";
my %dir = ();
$dir{output} = "oligo-calibrations";

local %infile = ();
local %outfile = ();

local $verbose = 0;
#local $in = STDIN;
local $out = STDOUT;

&ReadArguments();

################################################################
#### check argument values

## organism name
&CheckOrganismName($organism_name);

## output directory
unless (defined($dir{output})) {
    $dir{output} = "oligo_check_results";
}
&CheckOutputDir($dir{output});
chdir($dir{output});

#### output file
$outfile{output} = $organism_name."_".$oligo_len."nt_".$str.$ov."_n".$seq_nb."_l".$seq_len."_r".$repet."distrib.tab";

################################################################
### open output stream
$out = &OpenOutputFile($outfile{output});

################################################################
#### print verbose
&Verbose() if ($verbose);


################################################################
#### Retrieve all upstream sequences in wc format. This will allow to
#### load them rapidly for further analysis
warn "; ", &AlphaDate(), "\tRetrieving upstream sequences\n" if ($verbose >= 1);
$seq_file = $organism_name."_allup".$seq_len.".wc";
my $command = "retrieve-seq -org $organism_name ";
$command .= " -all -type upstream -nocomment -lw 0";
$command .= " -from -".$seq_len;
$command .= "  -to -1";
$command .= "  -format wc";
$command .= " -label orf";
$command .= " -o $seq_file";
print $out "; $command\n";
&doit($command);

################################################################
##### Random gene families
warn "; ", &AlphaDate(), "\tSelecting random gene families\n" if ($verbose >= 1);
$family_file = "random_genes.tab";
$command = "random-genes -org ".$organism_name;
$command .= " -r ".$repet;
$command .= " -n ".$seq_nb;
$command .= " -o ".$family_file;
print $out "; $command\n";
&doit($command);

################################################################
#### oligo-analysis
for my $r (1..$repet) {
    warn "; ", &AlphaDate(), "\toligo-analysis\trepetition\t",$r,"\n" if ($verbose >= 1);

    #### select one gene family
    my $one_family_file = "RAND_n".$seq_nb."_r".$r.".fam";
    $command = "grep 'RAND".$r."\$' ".$family_file;
    $command .= " | cut -f 1 > ".$one_family_file;
    print $out "; $command\n" if ($r == 1);
    &doit($command);

    #### select the corresponding upstream sequences
    my $oligo_file = "oligos_".$oligo_len."nt".$str.$ov."_n".$seq_nb."_r".$r;
    $command = "grep -f ".$one_family_file." ".$seq_file;
    $command .= " | oligo-analysis -format wc -return occ ".$str;
    $command .= " -v 1" if ($r == 1);
    $command .= " ".$ov;
    $command .= " -l ".$oligo_len;
    $command .= " -o ".$oligo_file;
    print $out "; $command\n" if ($r == 1);
    &doit($command);
}

################################################################
#### Regroup oligo-counts in a single table
#  my $oligo_table = $organism_name."_oligo_counts_".$oligo_len."nt_r".$repet.".tab";
#  warn "; ", &AlphaDate(), "\toligo table\t",$oligo_table,"\n" if ($verbose >= 1);
#  $command = "compare-scores -o ".$oligo_table;
#  $command .= " -null 0 -sc 3 -files oligos_".$oligo_len."nt_r*";
# print $out "; $command\n";
#&doit($command);

################################################################
#### Calculate distribution for each oligonucleotide
warn "; ", &AlphaDate(), "\toligo count distribution\t",$outfile{output},"\n" if ($verbose >= 1);
my $max_occ = 0;
my %counts = ();
my %count_sum = ();
for my $r (1..$repet) {
    my $oligo_file = "oligos_".$oligo_len."nt".$str.$ov."_n".$seq_nb."_r".$r;
    open OLIGOS, $oligo_file;
    while (<OLIGOS>) {
	next if (/^;/);
	chomp;
	my ($pattern, $id, $occ) = split "\t", $_;
#	warn $r, "\t", $pattern, "\t", $occ, "\n";
	$counts{$pattern}{$occ}++;
	$max_occ = $occ if ($max_occ < $occ);
	$count_sum{$pattern}++;
    }
    close OLIGOS;
}

#### print the distribution in the output file
print $out "; oligo count distribution\n";
print $out "; ", join ("\t", 0..$max_occ), "\n";
foreach my $pattern (sort keys %count_sum) {
    ## calculate number of families without any occurrence of this pattern
    $counts{$pattern}{0} = $repet - $count_sum{$pattern}; 
    print $out $pattern;
    for my $occ (0..$max_occ) {
	unless ($counts{$pattern}{$occ}) {
	    $counts{$pattern}{$occ} = 0;
	}
	printf $out "\t%d", $counts{$pattern}{$occ};
    }
    print $out "\n";
}


################################################################
###### execute the command

################################################################
###### print output


################################################################
###### finish verbose
if ($verbose) {
    my $done_time = &AlphaDate;
    print $out "; Job started $start_time\n";
    print $out "; Job done    $done_time\n";
}


################################################################
###### close output stream
close $out if ($outfile{output});


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
	calibrate-oligos.pl

        2002 by Jacques van Helden (jvanheld\@ucmb.ulb.ac.be)
	
DESCRIPTION

	Calibrate oligonucleotide frequencies in upstream sequences of
	a selected organism. The calibration can be gene-wise (count
	oligo frequencies in the upstream sequence of each gene) or
	cluster-wise (count oligo frequencies in upstream sequences of
	random gene selections).

CATEGORY
	util

USAGE
        calibrate-oligos.pl -org organism [-r #] [-start #]
	    [-sl #] [-sn #] [-ol #]

OPTIONS
	-h	(must be first argument) display full help message
	-help	(must be first argument) display options
	-v	verbose
	-outdir outputdir
	-org organism
	-r #	repetitions
	-start #	starting iteration (to pursue an interrupted test)
	-sl	sequence length
	-sn	sequence number
	-ol	oligo length
End_of_help
  close HELP;
  exit;
}

################################################################
#### display short help message
sub PrintOptions {
  open HELP, "| more";
  print HELP <<End_short_help;
calibrate-oligos.pl options
----------------
-h		(must be first argument) display full help message
-help		(must be first argument) display options
-i		input file
-o		output file
-outdir		output dir
-v		verbose
-org		organism
-r #		repetitions
-start #	starting iteration (to pursue an interrupted test)E
-sl		sequence length
-sn		sequence number
-ol		oligo length
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
	    
	    ### input file  
#	} elsif ($ARGV[$a] eq "-i") {
#	    $infile{input} = $ARGV[$a+1];
	    
#	    ### output file  
#	} elsif ($ARGV[$a] eq "-o") {
#	    $outfile{output} = $ARGV[$a+1];
#	    
	    #### organism
	} elsif ($ARGV[$a] eq "-org") {
	    $organism_name =$ARGV[$a+1];

	    ### repetitions
	} elsif ($ARGV[$a] eq "-r") {
	    $repet = $ARGV[$a+1];
	    
	    ### sequence length
	} elsif ($ARGV[$a] eq "-sl") {
	    $seq_len = $ARGV[$a+1];
	    
	    ### sequence number
	} elsif ($ARGV[$a] eq "-sn") {
	    $seq_nb = $ARGV[$a+1];
	    
	    ### oligo length
	} elsif ($ARGV[$a] eq "-ol") {
	    $oligo_len = $ARGV[$a+1];
	    
	    ### output file  
	} elsif ($ARGV[$a] eq "-outdir") {
	    $dir{output} = $ARGV[$a+1];
	    	
	    ### starting iteration
	} elsif ($ARGV[$a] eq "-start") {
	    $start = $ARGV[$a+1];
	    
	}
    }
}

################################################################
#### verbose message
sub Verbose {
    print $out "; calibrate-oligos.pl ";
    &PrintArguments($out);
    printf "; %-29s\t%s\n", "Output directory", $dir{output};
    printf "; %-29s\t%s\n", "Organism", $organism_name;
    printf "; %-29s\t%s\n", "Repetitions", $repet;
    printf "; %-29s\t%s\n", "Sequence length", $seq_len;
    printf "; %-29s\t%s\n", "Sequence number", $seq_nb;
    printf "; %-29s\t%s\n", "Sequence file", $seq_file;
    printf "; %-29s\t%s\n", "Random gene selections", $family_file;
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
