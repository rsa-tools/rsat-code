#!/usr/bin/perl -w
############################################################
#
# $Id: calibrate-oligos.pl,v 1.2 2003/12/25 19:52:18 jvanheld Exp $
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
my $start_time = &AlphaDate();
my $oligo_length = 6;
$repet = 1000;
my $oligo_len = 6;
my $str = "-2str";
my $seq_len=800;
my $seq_nb=10;
my $ov="-noov";
my $organism_name = "Saccharomyces_cerevisiae";
my %dir = ();
my $start = 1;

$seq_file = "";
$family_file = "";

@supported_tasks = qw (all upstream random oligos distrib clean_seq clean_oligos);
foreach my $task (@supported_tasks) {
    $supported_task{$task} = 1;
}
$supported_tasks = join ",", @supported_tasks;

local %infile = ();
local %outfile = ();

local $verbose = 0;

&ReadArguments();

################################################################
#### check argument values


#### check selected tasks
unless (defined(%task)) {
    &FatalError("You should select at least one task.");
}
if ($task{all}) {
    foreach my $t (@supported_tasks) {
	next if ($t =~ /^clean/); ## clean must be actively requested
	$task{$t} = 1;
    }
}

## organism name
&CheckOrganismName($organism_name);

## output directory
unless (defined($dir{output})) {
    $dir{output} = "results/".$organism_name."/rand_gene_selections/".$oligo_len."nt".$str.$ov."_N".$seq_nb."_L".$seq_len."_R".$repet;
}
&CheckOutputDir($dir{output});
$dir{oligos} = $dir{output}."/oligos";
&CheckOutputDir($dir{oligos});

chdir($dir{output});

################################################################
#### Open output stream
#### output files
$outfile{distrib} = $organism_name."_".$oligo_len."nt_".$str.$ov."_n".$seq_nb."_l".$seq_len."_r".$repet."_distrib.tab";
$outfile{stats} = $organism_name."_".$oligo_len."nt_".$str.$ov."_n".$seq_nb."_l".$seq_len."_r".$repet."_stats.tab";
if ($task{distrib}) {
    $out = &OpenOutputFile($outfile{distrib});
}

################################################################
#### print verbose
&Verbose() if ($verbose);

################################################################
#### Retrieve all upstream sequences in wc format. This will allow to
#### load them rapidly for further analysis
$seq_file = $organism_name."_allup".$seq_len.".wc";
if ($task{upstream}) {
    warn "; ", &AlphaDate(), "\tRetrieving upstream sequences\n" if ($verbose >= 1);
    my $command = "retrieve-seq -org $organism_name ";
    $command .= " -all -type upstream -nocomment -lw 0";
    $command .= " -from -".$seq_len;
    $command .= "  -to -1";
    $command .= "  -format wc";
    $command .= " -label orf";
    $command .= " -o $seq_file";
    print $out "; $command\n";
    &doit($command);
}

################################################################
##### Random gene families
$family_file = "random_genes.tab";
if ($task{random}) {
    warn "; ", &AlphaDate(), "\tSelecting random gene families\n" if ($verbose >= 1);
    $command = "random-genes -org ".$organism_name;
    $command .= " -r ".$repet;
    $command .= " -n ".$seq_nb;
    $command .= " -o ".$family_file;
    print $out "; $command\n";
    &doit($command);
}

################################################################
#### oligo-analysis
if ($task{oligos}) {
    for my $r ($start..$repet) {
	warn "; ", &AlphaDate(), "\toligo-analysis\trepetition\t",$r,"\n" if ($verbose >= 1);
	
	#### select one gene family
	my $one_family_file = "oligos/RAND_n".$seq_nb."_r".$r.".fam";
	$command = "grep 'RAND".$r."\$' ".$family_file;
	$command .= " | cut -f 1 > ".$one_family_file;
	print $out "; $command\n" if ($r == 1);
	&doit($command);
	
	#### select the corresponding upstream sequences
	my $oligo_file = "oligos/oligos_".$oligo_len."nt".$str.$ov."_L".$seq_len."_n".$seq_nb."_r".$r;
	$command = "grep -f ".$one_family_file." ".$seq_file;
	$command .= " | oligo-analysis -format wc -return occ ".$str;
	$command .= " -v 1" if ($r == 1);
	$command .= " ".$ov;
	$command .= " -l ".$oligo_len;
	$command .= " -o ".$oligo_file;
	print $out "; $command\n" if ($r == 1);
	&doit($command);
    }
}

################################################################
#### Regroup oligo-counts in a single table
#  my $oligo_table = $organism_name."_oligo_counts_".$oligo_len."nt_r".$repet.".tab";
#  warn "; ", &AlphaDate(), "\toligo table\t",$oligo_table,"\n" if ($verbose >= 1);
#  $command = "compare-scores -o ".$oligo_table;
#  $command .= " -null 0 -sc 3 -files oligos/oligos_".$oligo_len."nt_r*";
# print $out "; $command\n";
#&doit($command);

################################################################
#### Calculate distribution for each oligonucleotide
local $max_occ = 0;
local %counts = ();
local %count_sum = ();
local %overlaps = ();
local %sum = ();
local %ssq = ();
if ($task{distrib}) {
    warn "; ", &AlphaDate(), "\toligo count distribution\t",$outfile{distrib},"\n" if ($verbose >= 1);
    for my $r (1..$repet) {
	my $oligo_file = "oligos/oligos_".$oligo_len."nt".$str.$ov."_L".$seq_len."_n".$seq_nb."_r".$r;
#	my $oligo_file = "oligos/oligos_".$oligo_len."nt".$str.$ov."_L".$seq_len."_n".$seq_nb."_r".$r;
	&FatalError("File ".$oligo_file." does not exist in directory ".$dir{output}."\n") unless (-e $oligo_file);
	my ($oligos) = &OpenInputFile($oligo_file);
	while (<$oligos>) {
	    next if (/^;/);
	    chomp;
	    my ($pattern, $id, $occ, $ovl) = split "\t", $_;
#	warn $r, "\t", $pattern, "\t", $occ, "\n";
	    $counts{$pattern}{$occ}++;
	    $overlaps{$pattern}{$ovl}++;
	    $max_occ = $occ if ($max_occ < $occ);
	    $count_sum{$pattern}++;
	}
	close $oligos;
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
	    $sum{$pattern} += $occ*$counts{$pattern}{$occ};
	    $ssq{$pattern} += $occ*$occ*$counts{$pattern}{$occ};
	}
	print $out "\n";
    }
}

################################################################
#### Clean data files (just hold the distrib files)
if ($task{clean_seq}) {
    warn "; ", &AlphaDate(), "\tDeletingsequence file\n" if ($verbose >= 1);
    my $command = "rm -f ".$seq_file;
    &doit($command);
}
if ($task{clean_oligos}) {
    warn "; ", &AlphaDate(), "\tDeleting oligo files\n" if ($verbose >= 1);
    my $command = "rm -f ";
    $command .= " oligos/oligos_".$oligo_len."nt".$str.$ov."_L".$seq_len."_n".$seq_nb."_r*";
    $command .= " oligos/RAND*.fam";
    &doit($command);
}

################################################################
###### finish verbose
if ($verbose) {
    my $done_time = &AlphaDate();
    print $out "; Job started $start_time\n";
    print $out "; Job done    $done_time\n";
}

################################################################
#### Calculate statistics on each pattern count distribution
if ($task{distrib}) {
    warn "; ", &AlphaDate(), "\toligo count statistics\t",$outfile{stats},"\n" if ($verbose >= 1);
    $stats = &OpenOutputFile($outfile{stats});

    #### header
    print $stats join ("\t", 
		       "; pattern",
		       "sum",
		       "ssq",
		       "avg",
		       "var",
		       "std",
		       ), "\n";

    ## Pattern count statistics
    foreach my $pattern (sort keys %count_sum) {
	my $avg = $sum{$pattern}/$repet;
	my $var = $ssq{$pattern}/$repet - $avg*$avg;
	my $std = sqrt($var);
	print $stats join ("\t", 
			   $pattern,
			   $sum{$pattern},
			   $ssq{$pattern},
			   $avg,
			   $var,
			   $std,
			   ), "\n";
    }
}

warn "; Result stored in directory\t", $dir{output}, "\n";

################################################################
#### Close output streams
close $out;
close $stats;



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
	-task selected_task
		Select the tasks to be performed.
		Supported tasks: $supported_tasks

		Can be used iteratively on the same command line to 
		select multiple tasks.  

		Example:
		    -task upstream -task oligos -task distrib
		For a full analysis, simply type 
		    -task all
    Repetitions
	-r #	repetitions
	-start #	starting iteration (to pursue an interrupted test)

    Upstream sequences
	-org organism
	-sl	sequence length
	-sn	sequence number

    oligo-analysis
	-ol	oligo length
	-1str   strand-sensitive analysis
	-2str   strand-insensitive analysis
	-noov	prevent overlapping matches for self-overlapping patterms
		(default)
	-ovlp	allow overlapping matches for self-overlapping patterms
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
-outdir		output dir
-v		verbose
-task		selected task (supported: $supported_tasks)
-r #		repetitions
-start #	starting iteration (to pursue an interrupted test)E
-org		organism
-sl		sequence length
-sn		sequence number
-ol		oligo length
-1str   	strand-sensitive analysis
-2str   	strand-insensitive analysis
-noov		prevent overlapping matches for self-overlapping patterms
-ovlp		allow overlapping matches for self-overlapping patterms
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
	    
	    ### strands 
	} elsif ($ARGV[$a] eq "-1str") {
	    $force{strands} = "1str";
	} elsif ($ARGV[$a] eq "-2str") {
	    $force{strands} = "2str";

	    #### prevent self-overlap
	} elsif ($ARGV[$a] eq "-noov") {
	    $noov = 1;

	    #### allow self-overlap
	} elsif ($ARGV[$a] eq "-ovlp") {
	    $noov = 0;

	    #### prevent feature-map drawing

	    ### output directory  
	} elsif ($ARGV[$a] eq "-outdir") {
	    $dir{output} = $ARGV[$a+1];
	    	
	    ### starting iteration
	} elsif ($ARGV[$a] eq "-start") {
	    $start = $ARGV[$a+1];
	    
	    #### task selection
	} elsif ($ARGV[$a] eq "-task") {
	    my @requested_tasks = split ",", $ARGV[$a+1];
	    foreach my $task (@requested_tasks) {
		next unless $task;
		if ($supported_task{$task}) {
		    $task{$task} = 1;
		} else {
		    &FatalError("Unsupported task '$task'. \n\tSupported: $supported_tasks");
		}
	    }

	}
    }
}

################################################################
#### verbose message
sub Verbose {
    print $out "; calibrate-oligos.pl ";
    &PrintArguments($out);
    printf $out "; %-29s\t%s\n", "Output directory", $dir{output};
    printf $out "; %-29s\t%s\n", "Distribution", $outfile{distrib};
    printf $out "; %-29s\t%s\n", "Stats", $outfile{stats};
    printf $out "; %-29s\t%s\n", "Organism", $organism_name;
    printf $out "; %-29s\t%s\n", "Sequence length", $seq_len;
    printf $out "; %-29s\t%s\n", "Sequence number", $seq_nb;
    printf $out "; %-29s\t%s\n", "Repetitions", $repet;
    printf $out "; %-29s\t%s\n", "Sequence file", $seq_file;
    printf $out "; %-29s\t%s\n", "Random gene selections", $family_file;
    printf $out "; %-29s\t%s\n", "Oligonucleotide length", $oligo_len;
    printf $out "; %-29s\t%s\n", "Strands", $str;
    printf $out "; %-29s\t%s\n", "Overlap mode", $ov;
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
