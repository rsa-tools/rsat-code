#!/usr/bin/perl 
############################################################
#
# $Id: parse-yeast-from-mips.pl,v 1.3 2003/10/29 09:04:12 jvanheld Exp $
#
# Time-stamp: <2002-06-06 13:29:16 jvanheld>
#
############################################################
#use strict;;
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}
require "RSA.lib";


#### initialise parameters ####
my $start_time = &AlphaDate;

local %outfile = ();

local $verbose = 0;
local $out = STDOUT;

$dir{input} = "/win/databases/downloads/ftpmips.gsf.de/yeast";
$dir{output} = $RSA."/data/Saccharomyces_cerevisiae";
$dir{genome} = $dir{output}."/genome/";
$infile{orfs} = $dir{input}."/CYGD/all_orf_coor_descr";
$infile{synonyms} = $dir{input}."/CYGD/orf_gene_alias";
$chr_file = $dir{genome}."/Contigs_Saccharomyces_cerevisiae.txt";
$feature_file = $dir{genome}."Feature_Saccharomyces_cerevisiae.tab";
$synonym_file = $dir{genome}."Feature_names_Saccharomyces_cerevisiae.tab";

@chromosome_dirs = qw ( chri 
			chrii
			chriii
			chriv
			chrv
			chrvi
			chrvii
			chrviii
			chrix
			chrx
			chrxi
			chrxii
			chrxiii
			chrxiv
			chrxv
			chrxvi
			mito);

&ReadArguments;


#### check argument values ####

#### check for existence of input directory
unless (-d $dir{input}) {
    &FatalError("Input directory $dir{input} does not exist\n");
}
foreach my $dir (@chromosome_dirs) {
    unless (-d "$dir{input}/$dir") {
	&FatalError("Chromosome directory $dir does not exist\n");
    }
}

unless (-e $infile{orfs}) {
    &FatalError("ORF location file $infile{orfs} does not exist\n");
}

#### create output directory if necessary
unless (-d $dir{genome}) {
    system "mkdir -p $dir{genome}";
    unless (-d $dir{genome}) {
	&FatalError("Could not create output directory $dir{genome}\n");
    }
}



#### verbose ####
&Verbose if ($verbose);

&ConvertSequences();
&ParseFeatures();
&ParseSynonyms();

$command = "install-organism -v $verbose";
$command .= " -org Saccharomyces_cerevisiae";
$command .= " -organism 'Saccharomyces cerevisiae'";
$command .= " -source MIPS";
$command .= " -dir $dir{output}";
$command .= " -step config";
$command .= " -step start_stop";
$command .= " -step ncf";
$command .= " -step oligos";
$command .= " -step dyads";

warn $command, "\n" if ($verbose >= 1);
system $command;

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
	parse_yeast_genome_from_mips

        2001 by Jacques van Helden (jvanheld\@ucmb.ulb.ac.be)
	
USAGE
        parse_yeast_genome_from_mips [-i inputfile] [-o outputfile] [-v]

DESCRIPTION
	Parse the yeast genome imported from the MIPS. 

	Data source:
		ftp://ftpmips.gsf.de/yeast/

CATEGORY
	parser

OPTIONS
	-h	(must be first argument) display full help message
	-help	(must be first argument) display options
	-v	verbose
	-i inputfile
		if not specified, the standard input is used.
		This allows to place the command within a pipe.
	-o outputfile
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
parse_yeast_genome_from_mips options
----------------
-h	(must be first argument) display full help message
-help	(must be first argument) display options
-i	input file
-o	output file
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
	} elsif ($ARGV[0] eq "-help") {
	    &PrintOptions;
	    
	}
    }
}

sub Verbose {
    print $out "; parse_yeast_genome_from_mips ";
    &PrintArguments($out);
    printf $out "; %-29s\t%s\n", "Input directory", $dir{input};
    printf $out "; %-29s\t%s\n", "Output directory", $dir{genome};
}


sub ConvertSequences {
    #### chromosome list file
    open CHR, ">$chr_file"
	|| &FatalError("Cannot write chromosome sfile $chr_file\n");

    #### convert sequences into raw format (without any spacing) #########
    foreach my $chr (@chromosome_dirs) {
	my @ascii_files = glob($dir{input}."/".$chr."/${chr}*.ascii");
	if ($#ascii_files == -1) {
	    &FatalError("Cannot find the ascii sequence for chromosome $chr\n");
	} elsif ($#ascii_files > 0) {
	    &FatalError("Multiple ascii files in directory $chr\n");
	}
	my $raw_file= $dir{genome}."/$chr.raw";
	my $command = "convert-seq -i $ascii_files[0] -from raw -to raw -lw 0 -o $raw_file";
	warn $command, "\n" if ($verbose >= 1);;
	system $command;
	print CHR join ("\t",  
			$chr.raw,
			$chr,
			"linear"),"\n";
    }
    
    close CHR;
}

sub ParseFeatures {
    #### open output file
    open FEATURES, ">$feature_file"
	|| &FatalError("Cannot create feature file $feature_file\n");
    print FEATURES join ("\t",  
			 ";ID",	
			 "TYPE",
			 "NAME",
			 "CHROM",
			 "LEFT",
			 "RIGHT",
			 "STRAND",
			 "DESCR"), "\n" ; # header

    #### parse feature table
    $chr_conversion{A_CONTIG} = "chri";
    $chr_conversion{B_CONTIG} = "chrii";
    $chr_conversion{C_CONTIG} = "chriii";
    $chr_conversion{D_CONTIG} = "chriv";
    $chr_conversion{E_CONTIG} = "chrv";
    $chr_conversion{F_CONTIG} = "chrvi";
    $chr_conversion{G_CONTIG} = "chrvii";
    $chr_conversion{H_CONTIG} = "chrviii";
    $chr_conversion{I_CONTIG} = "chrix";
    $chr_conversion{J_CONTIG} = "chrx";
    $chr_conversion{K_CONTIG} = "chrxi";
    $chr_conversion{L_CONTIG} = "chrxii";
    $chr_conversion{M_CONTIG} = "chrxiii";
    $chr_conversion{N_CONTIG} = "chrxiv";
    $chr_conversion{O_CONTIG} = "chrxv";
    $chr_conversion{P_CONTIG} = "chrxvi";
    $chr_conversion{Q_CONTIG} = "mito";

    open ORFS, $infile{orfs} 
	|| &FatalError ("Cannot read ORF file $infile{orfs}") ;
    my $header = <ORFS>; # skip header row
    my $l = 0;
    
    while (<ORFS>) {
	chomp;
	$l++;
	my $type = "CDS";
	my @fields = split "\t";
	for my $f (0..$#fields) {
	    $fields[$f] =~ s/^\s+//;
	    $fields[$f] =~ s/\s+$//;
	}
	my $orf = $fields[0];
	my $name = $fields[1];
	$name = $orf unless ($name =~ /\S/);
	my $chr = $chr_conversion{$fields[2]};
	my $coordinates = $fields[3];
	my $classification = $fields[4];
	my $description = $fields[5];
	
	#### gene coordinates
	$coordinates =~ s/^\(//;
	$coordinates =~ s/\)$//;
	my $strand = "D";
	if ($coordinates =~ /\(C\)/) {
	    $strand = "R";
	    $coordinates = $`
	}
	my @exons = split ",", $coordinates;
	my @exon_starts = ();
	my @exon_ends = ();
	foreach my $e (@exons) {
	    my ($exon_start, $exon_end) = split "\-", $e;
	    push @exon_starts, $exon_start;
	    push @exon_ends, $exon_end;
	}
	my $start = $exon_starts[0];
	my $end = $exon_ends[$#exon_ends];
	my $left = &min(@exon_starts, @exon_ends);
	my $right = &max(@exon_starts, @exon_ends);

	#### correction for the missing stop codons
	if ($strand eq "D") {
	    $right += 3;
	} else {
	    $left -= 3;
	}


	print FEATURES join ("\t", 
			     $orf, 
			     $type,
			     $name,
			     $chr,
			     $left,
			     $right,
			     $strand,
			     $description,
			     $coordinates
			     ), "\n";
    }
    close ORFS;
    close FEATURES;
    
}

sub ParseSynonyms {
    open SYNONYMS, ">$synonym_file"
	|| &FatalError("Cannot create synonym file $synonym_file\n");
    open NAMES, "$infile{synonyms}" 
	|| &FatalError ("Cannot read synonym file $infile{synonyms}") ;
    my $header = <NAMES>; # skip header row
    my $l = 0;
    while (<NAMES>) {
	chomp;
	my @fields = split '\|';
	for my $f (0..$#fields) {
	    $fields[$f] =~ s/^\s+//;
	    $fields[$f] =~ s/\s+$//;
	}
	my $orf = $fields[0];
	my $gene = $fields[1];

	my @synonyms = ();
	push @synonyms, $gene if ($gene =~ /\S/);
	push @synonyms, (split ";", $fields[2]);
	#print "$orf\t",$#synonyms+1,"\n";
	foreach my $synonym (@synonyms) {
	    print SYNONYMS "$orf\t$synonym\n";
	}
    }
    close NAMES;
    close SYNONYMS;
    
}
