#!/usr/bin/perl
############################################################
#
# $Id: parse_genes.pl,v 1.2 2000/03/08 00:22:15 jvanheld Exp $
#
# Time-stamp: <2000-03-08 01:21:40 jvanheld>
#
############################################################

### parse_kegg.pl
### type parse_kegg.pl -h for info

#use strict;
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`"); ### add the program's directory to the lib path
}
require "PFBP_classes.pl";
require "PFBP_parsing_util.pl";


package main;

### initialization
$start_time = `date +%Y%m%d_%H%M%S`;

### files to parse
$KEGG_dir = "~jvanheld/pfmp/Databases/KEGG/downloaded";

$in_file{yeast} = "gunzip -c $KEGG_dir/genomes/genes/S.cerevisiae.ent.gz | ";
$in_file{human} = "gunzip -c $KEGG_dir/genomes/genes/H.sapiens.ent.gz | ";
$in_file{coli} = "gunzip -c $KEGG_dir/genomes/genes/E.coli.ent.gz | ";

$organism{yeast}->{name} = "Saccharomyces cerevisiae";
$organism{coli}->{name} = "Escherichia coli";
$organism{human}->{name} = "Homo sapiens";

$out_file{error} = "Gene.errors.txt";
$out_file{stats} = "Gene.stats.txt";
$out_file{genes} = "Gene.obj";
$out_file{synonyms} = "Gene.synonyms.tab";
$out_file{mldbm_genes} = "Gene.mldbm";
$out_file{expressions} = "Expression.obj";
$out_file{mldbm_expressions} = "Expression.mldbm";
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

@{$out_fields{'PFBP::Expression'}} = qw( id input output );

&ReadArguments;
&DefaultVerbose if ($warn_level >= 1);

### actualize parameters
if ($test) {
  print STDERR ";TEST\n" if ($warn_level >= 1);
  ### fast partial parsing for debugging
  foreach $key (keys %in_file) {
    $in_file{$key} .= " head -3000 |";
  }
}

### parse data from original files
while (my ($org, $file) = each %in_file) {
  &ParseKeggFile($file, 
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

$genes->dump_tables();
$expressions->dump_tables();

$genes->export('MLDBM',$out_file{mldbm_genes});
$expressions->export('MLDBM',$out_file{mldbm_expressions});

&PrintStats($out_file{stats}, @classes);

### report execution time
if ($warn_level >= 1) {
  $done_time = `date +%Y%m%d_%H%M%S`;
  print STDERR ";\n";
  print STDERR "; job started $start_time";
  print STDERR "; job done    $done_time";
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
    }
  }
}



sub CreateExpression {
### the field "definition" contains a reference to the Swissprot ID
### of the polypeptide that the gene codes for
### this sbroutine etxracts this information, and instantiates 
### objects in the PFBP::Expression class
  foreach my $gene (PFBP::Gene->get_objects()) {
    my $gene_id = $gene->get_attribute("id");
    if (my $definition = join(" ", $gene->get_attribute("definition"))) {
      if ($definition =~ /\[SP:(\S+)\]/) {
	while ($definition =~ /\[SP:(\S+)\]/) {
	  my $prot_id = $1;
	  if (my $expression = $expressions->new_object()) {
	    print STDERR "created expression for gene $gene_id\n" 
	      if ($warn_level >= 2);
	    $expression->set_attribute("input",$gene_id);
	    $expression->set_attribute("output",$prot_id);
	    $definition = $';
	  } else {
	    die "Error: could not instantiated expression\n";
	  }
	}
      } else {
	print ERR "WARNING: gene $gene_id could not be linked to swissprot\n";
	print ERR "\tdefinition\t$definition\n";
      }
    } else {
      print ERR "WARNING: gene $gene_id has no definition field\n";
    }
  }
}

