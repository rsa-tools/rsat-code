#!/usr/bin/perl
############################################################
#
# $Id: parse_genes.pl,v 1.3 2000/03/10 14:22:59 jvanheld Exp $
#
# Time-stamp: <2000-03-10 02:17:58 jvanheld>
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

### initialization
$start_time = &AlphaDate;

### files to parse
$KEGG_dir = "~jvanheld/pfmp/Databases/KEGG/downloaded";


@selected_organisms= ();

$in_file{yeast} = "gunzip -c $KEGG_dir/genomes/genes/S.cerevisiae.ent.gz | ";
$in_file{human} = "gunzip -c $KEGG_dir/genomes/genes/H.sapiens.ent.gz | ";
$in_file{ecoli} = "gunzip -c $KEGG_dir/genomes/genes/E.coli.ent.gz | ";

$organism{yeast}->{name} = "Saccharomyces cerevisiae";
$organism{ecoli}->{name} = "Escherichia coli";
$organism{human}->{name} = "Homo sapiens";


$warn_level = 0;
$out_format = "obj";

#### classes and classholders
@classes = qw( PFBP::Gene PFBP::Expression );
$genes = PFBP::ClassFactory->new_class(object_type=>"PFBP::Gene",
				       prefix=>"gene_");
$expressions = PFBP::ClassFactory->new_class(object_type=>"PFBP::Expression",
					     prefix=>"expr_");


### open error report file
open ERR, ">$out_file{error}" || die "Error: cannot write error file $out_file{error}\n";

### default output fields for each class
@{$out_fields{'PFBP::Gene'}} = qw( id organism names position definition source);

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

$out_file{error} = "Gene".$suffix.".errors.txt";
$out_file{stats} = "Gene".$suffix.".stats.txt";
$out_file{genes} = "Gene".$suffix.".obj";
$out_file{synonyms} = "Gene".$suffix.".synonyms.tab";
$out_file{mldbm_genes} = "Gene".$suffix.".mldbm";
$out_file{expressions} = "Expression".$suffix.".obj";
$out_file{mldbm_expressions} = "Expression".$suffix.".mldbm";

### select all organisms if none was selected (-org)
unless ($#selected_organisms >= 0) {
  push @selected_organisms, "ecoli";
  push @selected_organisms, "human";
  push @selected_organisms, "yeast";
}

&DefaultVerbose if ($warn_level >= 1);
warn "; Selected organisms\n;\t", join("\n;\t", @selected_organisms), "\n"
    if ($warn_level >= 1);

### actualize parameters
if ($test) {
  print STDERR ";TEST\n" if ($warn_level >= 1);
  ### fast partial parsing for debugging
  foreach $key (keys %in_file) {
    $in_file{$key} .= " head -3000 |";
  }
}


### parse data from original files
foreach $org (@selected_organisms) {
  &ParseKeggFile($in_file{$org}, 
		 $genes, 
		 source=>"KEGG", 
		 organism=>$organism{$org}->{name});
}
$genes->index_names();


&CreateExpression;
$expressions->index_names();

### print result
&ExportClasses($out_file{genes}, $out_format, PFBP::Gene);
&ExportClasses($out_file{expressions}, $out_format, PFBP::Expression);


$genes->dump_tables($suffix);
$expressions->dump_tables($suffix);

$genes->export('MLDBM',$out_file{mldbm_genes});
$expressions->export('MLDBM',$out_file{mldbm_expressions});

&PrintStats($out_file{stats}, @classes);

### report execution time
if ($warn_level >= 1) {
  $done_time = &AlphaDate;
  print STDERR ";\n";
  print STDERR "; job started $start_time\n";
  print STDERR "; job done    $done_time\n";
}

close ERR;

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



sub CreateExpression {
### the field "definition" contains a reference to the Swissprot ID
### of the polypeptide that the gene codes for
### this sbroutine etxracts this information, and instantiates 
### objects in the PFBP::Expression class

  ### load polypeptide class 
  ### for converting Swissprot IDs in ACs
#  $class_holder{polypeptides} = &ImportClass("$data_dir/polypeptides/Polypeptide.mldbm", "PFBP::Polypeptide");
 
  
  $infile{polypeptide_index} = "$mldbm_dir/Polypeptide".$suffix."__name_index.mldbm";
  warn "; opening polpeptide index $infile{polypeptide_index}\n"
      if ($warn_level >= 1);

  ### open the persistent polypeptide index
  tie %index, 'MLDBM', $infile{polypeptide_index} ||
      die "Error: cannot tie file $infile{polypeptide_index}";
  
  foreach my $gene (PFBP::Gene->get_objects()) {
    my $gene_id = $gene->get_attribute("id");
    if (my $definition = join(" ", $gene->get_attribute("definition"))) {
      if ($definition =~ /\[SP:(\S+)\]/) {
	while ($definition =~ s/\[SP:(\S+)\]//) {
	  
	  ### identify the polypeptide
	  my $swissprot_id = $1;
	  my $polypeptide = undef;
	  my $polypeptide_id = undef;
	  
	  if (exists $index{uc($swissprot_id)}) {
	    $id_table = $index{uc($swissprot_id)};
	    $polypeptide_id = $$id_table[0];
	  } else {
	    $error_message = "ERROR : could not identify polypeptide $swissprot_id\n";
	    print ERR $error_message;
	    warn  $error_message if ($warn_level >= 2);
	    next;
	  }
	  
	  $expr_id = "expr_".$gene_id."_".$polypeptide_id;
	  if (my $expression = $expressions->new_object(id=>$expr_id,
							source=>"KEGG")) {
	    print STDERR "created expression from gene $gene_id to polypeptide $polypeptide_id\n" 
		if ($warn_level >= 2);
	    $expression->set_attribute("gene_id",$gene_id);
	    $expression->set_attribute("polypeptide_id",$polypeptide_id);
	  } else {
	    die "Error: could not instantiate expression\n";
	  }
	}
      } else {
	$error_message =  "WARNING: gene $gene_id could not be linked to swissprot\n";
	$error_message .= "\tdefinition\t$definition\n";
	print ERR $error_message;
	warn $error_message if ($warn_level >= 3);
      }
    } else {
      print ERR "WARNING: gene $gene_id has no definition field\n";
    }
    
    untie %polypeptides_name_index ||
	die "Error: cannot untie file $infile{polypeptide_index}";
    
  }
}

