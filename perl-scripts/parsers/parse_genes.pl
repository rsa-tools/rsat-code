#!/usr/bin/perl
############################################################
#
# $Id: parse_genes.pl,v 1.18 2000/12/29 21:12:58 jvanheld Exp $
#
# Time-stamp: <2000-12-29 22:12:07 jvanheld>
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

$dir{output} = $parsed_data."/kegg_parsed/".$delivery_date;
#$dir{genes} = $dir{KEGG}."/genomes/previous_genes/";
$dir{genes} = $dir{KEGG}."/genomes/genes/";

$kegg_file{yeast} = "S.cerevisiae.ent";
$kegg_file{human} = "H.sapiens.ent";
$kegg_file{ecoli} = "E.coli.ent";

@all_files = ();
push @all_files, glob($dir{genes}."/*\.ent");
push @all_files, glob($dir{genes}."/*\.ent\.gz");
foreach $file (@all_files) {
    if ($file =~ /\/([^\/]*)$/) {
	$file = $1;
    }
    my $organism = $file;
    $organism =~ s/\.ent$//;
    $organism =~ s/\.ent\.gz$//;
    $organism{$organism}->{name} = $organism;
    $kegg_file{lc($organism)} = $file;
    push @all_organisms, $organism;
}




$organism{yeast}->{name} = "Saccharomyces cerevisiae";
$organism{ecoli}->{name} = "Escherichia coli";
$organism{human}->{name} = "Homo sapiens";

$warn_level = 0;
$out_format = "obj";

#### classes and classholders
@classes = qw( PFBP::Gene );
$genes = PFBP::ClassFactory->new_class(object_type=>"PFBP::Gene",
				       prefix=>"gene_");
### default output fields for each class


&ReadArguments;

#### specific export format for RSA-tools
if ($rsa) {
    $single_name = 1;
    $genes->set_out_fields(qw( id type name chromosome start end strand description position names ));
} else {
    $genes->set_out_fields(qw( id source organism type position chromosome strand start end description names exons introns ));
    #@{$out_fields{'PFBP::Gene'}} = qw( id source organism raw_position chromosome strand start_base end_base description names exons );
}

### outfile names
unless (defined($suffix)) {
    foreach my $organism (@selected_organisms) {
	$suffix .= "_$organism";
    }
    if ($export{enzymes}) {
	$suffix .= "_enz";
    }
    $suffix .= "_test" if ($test);
}

unless (-d $dir{output}) {
    warn "Creating output dir $dir{output}\n";
    mkdir $dir{output}, 0775 || die "Error: cannot create directory $dir\n";
}
chdir $dir{output};
$out_file{error} = "$dir{output}/Gene".$suffix.".errors.txt";
$out_file{stats} = "$dir{output}/Gene".$suffix.".stats.txt";
$out_file{genes} = "$dir{output}/Gene".$suffix.".obj" if ($export{obj});

### open error report file
open ERR, ">$out_file{error}" || die "Error: cannot write error file $out_file{error}\n";


### select all organisms if none was selected (-org)
unless ($#selected_organisms >= 0) {
    @selected_organisms = @all_organisms;
}

foreach $org (@selected_organisms) {
    die "; Fatal Error: no data file for organism $org\n"
	unless (defined($kegg_file{lc($org)}));
    my $data_file = $dir{genes}.$kegg_file{lc($org)};
    if (-e $data_file) {
	if ($data_file =~ /\.gz$/) {
	    $in_file{$org} = "gunzip -c ${data_file} | ";
	} else {
	    $in_file{$org} = "cat ${data_file} | ";
	}
    } elsif (-e "${data_file}.gz") {
	$in_file{$org} = "gunzip -c ${data_file}.gz | ";
    } else {
	die ("Error: cannot find data file for organism ", $org, "\n",
	     "\t", $data_file{$org}, "\n");
    }
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
		   source=>"KEGG:".$kegg_file{lc($org)});

    ### check organism attribute
    foreach $gene ($genes->get_objects()) {
	unless ($gene->get_attribute("organism")) {
	    &ErrorMessage("Warning: gene ", $gene->get_attribute("id"), " has no organism attribute\n");
	    $gene->set_attribute("organism",$org);
	    next;
	}
    }
}

&ParsePositions();

foreach $gene ($genes->get_objects()) {
    ### define a single name  (take the first value in the name list)
    if ($single_name) {
	if ($name = $gene->get_name()) {
	    $gene->set_attribute("name",$name);
	} else {
	    $gene->set_attribute("name",$gene->get_id());
	}
    }

    #### check for genes without description
    if ($gene->get_attribute("description") eq "<UNDEF>") {
	$gene->set_attribute("description",$gene->get_name());
    }
}


### print result
&PrintStats($out_file{stats}, @classes);
$genes->dump_tables($suffix);
&ExportClasses($out_file{genes}, $out_format, PFBP::Gene) if $export{obj};


### report execution time
if ($warn_level >= 1) {
  $done_time = &AlphaDate;
  warn ";\n";
  warn "; job started $start_time";
  warn "; job done    $done_time\n";
}

close ERR;


warn "; compressing the files\n" if ($warn_level >= 1);
system "gzip -f $dir{output}/*.tab $dir{output}/*.txt";
system "gzip -f $dir{output}/*.obj" if ($export{obj});


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
	
OPTIONS	
	-h	detailed help
	-help	short list of options
	-test	fast parsing of partial data, for debugging
	-outdir output directory
	-w #	warn level
		Warn level 1 corresponds to a restricted verbose
		Warn level 2 reports all polypeptide instantiations
		Warn level 3 reports failing get_attribute()
	-obj    Export results in "obj" format (human readable)
	-enz	export enzymes only
	-all	export all polypeptides for the selected organisms
		(default)
	-org	select an organism for exportation
		can be used reiteratively in the command line 
		to select several organisms
		by default, all organisms found in the input directory
		are selected
	-suffix suffix
		add a suffix to output file names
		    -suffix '' prevents from adding any suffix
	-name
		exports a name as single value attribute in
		the main table (this is redundant but can be useful)

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
      
      ### export single name in main table
    } elsif ($ARGV[$a] eq "-name") {
      $main::single_name = 1;
      
      ### specific export format  for RSA-tools
    } elsif ($ARGV[$a] eq "-rsa") {
      $main::rsa = 1;
      
      ### suffix
    } elsif ($ARGV[$a] eq "-suffix") {
      $a++;
      $main::suffix = $ARGV[$a];
      
      ### output dir
    } elsif ($ARGV[$a] eq "-outdir") {
      $a++;
      $main::dir{output} = $ARGV[$a];
      
      ### help
    } elsif (($ARGV[$a] eq "-h") ||
	     ($ARGV[$a] eq "-help")) {
      &PrintHelp;
      exit(0);

      ### select enzymes for exportation
    } elsif ($ARGV[$a] =~ /^-org/) {
      push @selected_organisms, $ARGV[$a+1];
      
      ### select enzymes for exportation
    } elsif ($ARGV[$a] =~ /^-enz/) {
      $main::export{enzymes} = 1;

    } elsif ($ARGV[$a] =~ /^-all/) {
      $main::export{all} = 1;

      #### export object file
    } elsif ($ARGV[$a] =~ /^-obj/) {
      $main::export{obj} = 1;
    }
    
  }
}





sub ParsePositions {
### clean up positions
    warn ("; ",
	  &AlphaDate(),
	  "\tparsing gene positions\n")
	if ($warn_level >= 1);

    foreach my $gene ($genes->get_objects()) {
	my $position = $gene->get_attribute("position");

	if ($position eq "<UNDEF>") {
	    $gene->set_attribute("position","<NULL>");
	    $gene->set_attribute("chromosome", "<NULL>");
	    $gene->set_attribute("strand","<NULL>");
	    $gene->set_attribute("start","<NULL>");
	    $gene->set_attribute("end","<NULL>");
	    &ErrorMessage("Warning: gene ", $gene->get_attribute("id"), " has no position attribute\n");
	    next;
	}
	my $coord = "<NULL>";
	my $chomosome = "<NULL>";
	my $chrom_pos = "<NULL>";
	my $strand = "<NULL>";
	my $start = "<NULL>";
	my $end = "<NULL>";

	if ($position =~ /^(\S+)\:(.*)/) {
	    $chromosome = $1;
	    $chrom_pos = $2;
	} elsif ($position =~ /^([^\:]*)$/) {
	    $chromosome = "genome";
	    $chrom_pos = $1;
	} else {
	    &ErrorMessage("Warning: invalid position",
			  "\t", $gene->get_attribute("id"),
			  "\t", $gene->get_attribute("organism"),
			  "\t", $position, "\n");
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
	    my @exons = split ",", $1;
	    my @exon_starts = ();
	    my @exon_ends = ();

	    #### exon limits
	    foreach my $exon (@exons) {
		$exon =~ s/\s+//g;
		$gene->push_attribute("exons", $exon);
		if ($exon =~ /(\d+)\.\.(\d+)/) {
		    push @exon_starts, $1;
		    push @exon_ends, $2;
		} else {
		    &ErrorMessage("Error gene\t",$gene->get_id(),"\tinvalid exon\t$exon\n");
		}
	    }
	    
	    #### gene start and end 
	    $start = $exon_starts[0] || "<NULL>";
	    $end = $exon_ends[$#exons] || "<NULL>";

	    #### introns
	    my @introns = ();
	    for my $e (0..$#exon_starts - 1) {
		my $intron = $exon_ends[$e] + 1;
		$intron .= "..";
		$intron .= $exon_starts[$e+1] -1;
		$gene->push_attribute("introns", $intron);
	    }
	} else {
	    &ErrorMessage("Warning : gene ",$gene->get_attribute("id"),"\tinvalid gene position $position\n");
	    $strand = "<NULL>";
	    $start = "<NULL>";
	    $end = "<NULL>";
	}
	$gene->set_attribute("chromosome", $chromosome);
	$gene->set_attribute("strand",$strand);
	$gene->set_attribute("start",$start);
	$gene->set_attribute("end",$end);
    }
}

