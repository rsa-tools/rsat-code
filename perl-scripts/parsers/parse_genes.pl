#!/usr/bin/perl
############################################################
#
# $Id: parse_genes.pl,v 1.11 2000/11/29 14:00:48 jvanheld Exp $
#
# Time-stamp: <2000-11-29 13:11:53 jvanheld>
#
############################################################

### parse_kegg.pl
### type parse_kegg.pl -h for info

#use strict;
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`"); ### add the program's directory to the lib path
}
require "PFBP_classes.pl";
require "PFBP_config.pl";
require "PFBP_util.pl";
require "PFBP_loading_util.pl"; ### for converting polypeptide IDs into ACs
require "PFBP_parsing_util.pl";


package main;

### files to parse
@selected_organisms= ();

$species_name{yeast} = "S.cerevisiae";
$species_name{human} = "H.sapiens";
$species_name{ecoli} = "E.coli";

foreach $org (keys %species_name) {
    my $data_file = $dir{KEGG}."/genomes/genes/".$species_name{$org}.".ent";
    if (-e $data_file) {
	$in_file{$org} = "cat ${data_file} | ";
    } elsif (-e "${data_file}.gz") {
	$in_file{$org} = "gunzip -c ${data_file}.gz | ";
    } else {
	die ("Error: cannot find data file for organism ", $org, "\n",
	     "\t", $data_file{$org}, "\n");
    }
}

$organism{yeast}->{name} = "Saccharomyces cerevisiae";
$organism{ecoli}->{name} = "Escherichia coli";
$organism{human}->{name} = "Homo sapiens";


$warn_level = 0;
$out_format = "obj";

#### classes and classholders
#@classes = qw( PFBP::Gene PFBP::Expression );
@classes = qw( PFBP::Gene );
$genes = PFBP::ClassFactory->new_class(object_type=>"PFBP::Gene",
				       prefix=>"gene_");
$expressions = PFBP::ClassFactory->new_class(object_type=>"PFBP::Expression",
					     prefix=>"expr_");

### default output fields for each class
@{$out_fields{'PFBP::Gene'}} = qw( id organism names position chromosome strand start end exons definition source);

@{$out_fields{'PFBP::Expression'}} = qw( id gene_id polypeptide_id );

&ReadArguments;

### outfile names
$suffix = "";
foreach $organism (@selected_organisms) {
  $suffix .= "_$organism";
}
if ($export{all}) {
  $suffix .= "_all";
} elsif ($export{enzymes}) {
  $suffix .= "_enz";
}
$suffix .= "_test" if ($test);

$dir{output} = $parsed_data."/kegg_parsed/".$delivery_date;
unless (-d $dir{output}) {
    warn "Creating output dir $dir{output}";
    mkdir $dir{output}, 0775 || die "Error: cannot create directory $dir\n";
}
chdir $dir{output};
$out_file{error} = "$dir{output}/Gene".$suffix.".errors.txt";
$out_file{stats} = "$dir{output}/Gene".$suffix.".stats.txt";
$out_file{genes} = "$dir{output}/Gene".$suffix.".obj";
$out_file{synonyms} = "$dir{output}/Gene".$suffix.".synonyms.tab";

### open error report file
open ERR, ">$out_file{error}" || die "Error: cannot write error file $out_file{error}\n";


### select all organisms if none was selected (-org)
unless ($#selected_organisms >= 0) {
  push @selected_organisms, "ecoli";
  push @selected_organisms, "human";
  push @selected_organisms, "yeast";
}

### test conditions
if ($test) {
  warn ";TEST\n" if ($warn_level >= 1);
  ### fast partial parsing for debugging
  foreach $key (keys %in_file) {
    $in_file{$key} .= " head -1000 |";
  }
}

&DefaultVerbose if ($warn_level >= 1);
warn "; Selected organisms\n;\t", join("\n;\t", @selected_organisms), "\n"
    if ($warn_level >= 1);


### parse data from original files
foreach $org (@selected_organisms) {
  &ParseKeggFile($in_file{$org}, 
		 $genes, 
		 source=>"KEGG", 
		 organism=>$organism{$org}->{name});
}
$genes->index_names();

&ParsePositions();

foreach $gene ($genes->get_objects()) {
    ### define a primary name (take the first value in the name list)
    if ($name = $gene->get_name()) {
	$gene->set_attribute("primary_name",$name);
    }
    #### check for genes without definition
    if ($gene->get_attribute("definition") eq "<UNDEF>") {
	$gene->set_attribute("definition",$gene->get_name());
    }
}


### print result
&ExportClasses($out_file{genes}, $out_format, PFBP::Gene);

$genes->dump_tables($suffix);

&PrintStats($out_file{stats}, @classes);

### report execution time
if ($warn_level >= 1) {
  $done_time = &AlphaDate;
  warn ";\n";
  warn "; job started $start_time\n";
  warn "; job done    $done_time\n";
}

close ERR;

system "gzip -f $dir{output}/*.tab $dir{output}/*.obj $dir{output}/*.txt";


exit(0);

### subroutines for the main package

### print the help message 
### when the program is called with the -h or -help option
sub PrintHelp {
  open HELP, "| more";
  print <<EndHelp;
NAME
	parse_kegg.pl

DESCRIPTION
	Parse genes from a KEGG file (http://www.genome.ad.jp/kegg/). 
	
AUTHOR
	Jacques van Helden (jvanheld\@ucmb.ulb.ac.be)  

VERSION
	0.01
	Created		1999/12/16
	Last modified	2000/01/08
	
SYNOPSIS	
	parse_kegg_genes.pl [-v] [-vv] [-i infile] [-format output_format]
		[-o outfile] 

OPTIONS	
	-h	detailed help
	-help	short list of options
	-test	fast parsing of partial data, for debugging
	-w #	warn level
		Warn level 1 corresponds to a restricted verbose
		Warn level 2 reports all polypeptide instantiations
		Warn level 3 reports failing get_attribute()
	-i	input file
		If ommited, STDIN is used
		This allows to insert the program within a unix pipe
	-o	output file
		If ommited, STDOUT is used. 
		This allows to insert the program within a unix pipe
	-enz	export enzymes only
	-all	export all polypeptides for the selected organisms
		(default)
	-org	select an organism for exportation
		can be used reiteratively in the command line 
		to select several organisms
		   Supported organisms :
			ecoli
			human
			yeast
		by default, all organisms are selected
EndHelp
  close HELP;
}

  



### read arguments from the command line
sub ReadArguments {
  for my $a (0..$#ARGV) {
    
    ### warn level
    if (($ARGV[$a] eq "-w" ) && 
	($ARGV[$a+1] =~ /^\d+$/)){
      $main::warn_level = $ARGV[$a+1];
      
      ### test run
    } elsif ($ARGV[$a] eq "-test") {
      $main::test = 1;
      
      ### input file
    } elsif ($ARGV[$a] eq "-i") {
      $a++;
      $main::in_file{pathway_index} = $ARGV[$a];
      
      ### output file
    } elsif ($ARGV[$a] eq "-o") {
      $a++;
      $main::out_file{pathways} = $ARGV[$a];
      
      ### help
    } elsif (($ARGV[$a] eq "-h") ||
	     ($ARGV[$a] eq "-help")) {
      &PrintHelp;
      exit(0);

      ### select enzymes for exportation
    } elsif ($ARGV[$a] =~ /^-org/) {
      push @selected_organisms, lc($ARGV[$a+1]);
      
      ### select enzymes for exportation
    } elsif ($ARGV[$a] =~ /^-enz/) {
      $main::export{enzymes} = 1;
    } elsif ($ARGV[$a] =~ /^-all/) {
      $main::export{all} = 1;
    }
    
  }
}





sub ParsePositions {
### clean up positions
    warn ("; parsing gene positions\n")
	if ($warn_level >= 1);

    foreach my $gene ($genes->get_objects()) {
	my $organism = $gene->get_attribute("organism");
	unless (defined($organism)) {
	    &ErrorMessage("Error: gene ", $gene->get_attribute("id"), " has no organism attribute\n");
	    next;
	}
	my $position = $gene->get_attribute("position");

	if ($position eq "<UNDEF>") {
	    $gene->set_attribute("position","<NULL>");
	    $gene->set_attribute("chromosome", "<NULL>");
	    $gene->set_attribute("strand","<NULL>");
	    $gene->set_attribute("start","<NULL>");
	    $gene->set_attribute("end","<NULL>");
	    &ErrorMessage("Error: gene ", $gene->get_attribute("id"), " has no position attribute\n");
	    next;
	}
	my $coord;
	my $chomosome;
	my $chrom_pos;
	my $strand;
	my $start;
	my $end;
	if ($organism eq "Escherichia coli") {
	    $chromosome = "genome";
	    $chrom_pos = $position;
	} elsif ($organism eq "Saccharomyces cerevisiae") {
	    if ($position =~ /^(\d+)\:(.*)/) {
		$chromosome = $1;
		$chrom_pos = $2;
	    } elsif ($position =~ /^mit\:(.*)/) {
		$chromosome = "mitochondrial";
		$chrom_pos = $1;
	    } else {
		&ErrorMessage("Error in gene ",$gene->get_attribute("id"),"\tinvalid position\t$position\n");
		next;
	    }
	} elsif ($organism eq "Homo sapiens") {
	    $gene->set_attribute("chromosome", "<NULL>");
	    $gene->set_attribute("strand","<NULL>");
	    $gene->set_attribute("start","<NULL>");
	    $gene->set_attribute("end","<NULL>");
	    next;
	}
	
	if ($chrom_pos =~ /complement\((.*)\)/) {
	    $strand = "R";
	    $coord = $1;
	} else {
	    $strand = "D";
	    $coord = $chrom_pos;
	}
	
	if ($coord =~ /^(\d+)\.\.(\d+)$/) {
	    $start = $1;
	    $end = $2;
	} elsif ($coord =~ /^join\((.*)\)$/) { ### exons
	    @exons = split ",", $1;
	    foreach $exon (@exons) {
		$gene->push_attribute("exons", $exon);
	    }
	    if ($exons[0] =~ /^(\d+)\.\.(\d+)$/) {
		$start = $1;
	    }
	    if ($exons[$#exons] =~ /^(\d+)\.\.(\d+)$/) {
		$end = $2;
	    }
	} else {
	    &ErrorMessage("Error : gene ",$gene->get_attribute("id"),"\tinvalid gene coordinate $coord\n");
	    next;
	}
	$gene->set_attribute("chromosome", "genome");
	$gene->set_attribute("strand",$strand);
	$gene->set_attribute("start",$start);
	$gene->set_attribute("end",$end);
    }
}

