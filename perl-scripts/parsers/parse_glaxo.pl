#!/usr/bin/perl
############################################################
#
# $Id: parse_glaxo.pl,v 1.1 2000/10/25 07:10:34 jvanheld Exp $
#
# Time-stamp: <2000-10-25 09:07:19 jvanheld>
#
############################################################
### parse_ligand.plt
### type parse_ligand.pl -h for info

#use strict;
if ($0 =~ /([^(\/)]+)$/) {
  push (@INC, "$`"); ### add the program's directory to the lib path
}
require "PFBP_classes.pl";
require "PFBP_config.pl";
require "PFBP_util.pl";
require "PFBP_parsing_util.pl";

package main;

### initialization
$start_time = `date +%Y-%m-%d.%H%M%S`;
$warn_level = 0;

### files to parse
$in_file{glaxo} = "gunzip -c  $data_glaxo/didier.fdt.gz | ";

$output_dir = "$parsed_data/glaxo";
unless (-d $output_dir) {
    warn "Creating output dir $output_dir";
    mkdir $output_dir, 0755 || die "Error: cannot create directory $dir\n";
}
chdir $output_dir;
$out_file{compounds} = "$output_dir/Compound.obj";
$out_file{stats} = "$output_dir/Compound.stats.txt";
$out_file{errors} = "$output_dir/Compound.errors.txt";

open ERR, ">$out_file{errors}" || die "Error: cannot write error report fle $$out_file{errors}\n";

$out_format = "obj";

@classes = qw( PFBP::Compound );

### default output fields for each class
#@{$out_fields{'PFBP::Compound'}} = qw( id names source );

&ReadArguments;

### testing mode
if ($test) {
  warn ";TEST\n" if ($warn_level >= 1);
  ### fast partial parsing for debugging
  $in_file{glaxo} .= " head -5000 |";
}

### default verbose message
&DefaultVerbose if ($warn_level >= 1);

### parse data from original files
$compounds = PFBP::ClassFactory->new_class(object_type=>"PFBP::Compound",
					   prefix=>"g_");
&ParseGlaxoFile($in_file{glaxo}, $compounds, source=>'LIGAND:compound');
$compounds->index_names();
### define a primary name (take the first value in the name list)
foreach $compound ($compounds->get_objects()) {
    @BRENDA_names = $compound->get_attribute("BRENDA_name");
    if ($#BRENDA_names >= 0) {
	$compound->set_attribute("primary_name",$BRENDA_names[0]);
    } elsif ($name = $compound->get_name()) {
	$compound->set_attribute("primary_name",$name);
    }
}


### print result
$compounds->dump_tables();
&ExportClasses($out_file{compounds}, $out_format, PFBP::Compound);

### print some stats after parsing
&PrintStats($out_file{stats}, @classes);

### report execution time
if ($warn_level >= 1) {
  $done_time = `date +%Y-%m-%d.%H%M%S`;
  warn ";\n";
  warn "; job started $start_time";
  warn "; job done    $done_time";
}

close ERR;

system "gzip -f $output_dir/*.tab $output_dir/*.obj $output_dir/*.txt";

exit(0);

### subroutines for the main package

### print the help message 
### when the program is called with the -h or -help option
sub PrintHelp {
  open HELP, "| more";
  print <<EndHelp;
NAME
        parse_ligand.pl

DESCRIPTION
	Parse	compounds from Glaxo
	
AUTHOR
	Jacques van Helden (jvanheld\@ucmb.ulb.ac.be)  

VERSION
	0.01
	Created		2000/10/25
	Last modified	2000/10/25
	
SYNOPSIS	
	parse_ligand.pl [-v] [-vv] [-i infile] [-format output_format]
		[-o outfile] [-smiles]

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
  my $a = "";
  for $a (0..$#ARGV) {

    ### warn level
    if (($ARGV[$a] eq "-w" ) && 
	     ($ARGV[$a+1] =~ /^\d+$/)){
      $main::warn_level = $ARGV[$a+1];
      
      ### fast test
    } elsif ($ARGV[$a] eq "-test") {
      $test = 1;

      ### input file
    } elsif ($ARGV[$a] eq "-i") {
      $a++;
      $in_file{compounds} = $ARGV[$a];
      
      ### output file
    } elsif ($ARGV[$a] eq "-o") {
      $a++;
      $out_file{compounds} = $ARGV[$a];

      ### help
    } elsif (($ARGV[$a] eq "-h") ||
	     ($ARGV[$a] eq "-help")) {
      &PrintHelp;
      exit(0);
    }
  }
}


sub ParseGlaxoFile{
  my ($input_file,$class_holder, %args) = @_;
  my $class = $class_holder->get_object_type();
  
  ### open the input stream
  if ($input_file) {
    warn (";\n; ", &AlphaDate,
	  " parsing class $class from $input_file\n")
      if ($warn_level >= 1);
    open COMP, $input_file 
      || die "Error: cannot $input_file\n";
    $in = COMP;
  } else {
    warn (";\n; ", &AlphaDate,
	  " parsing class $class from STDIN\n")
      if ($warn_level >= 1);
    $in = STDIN;
  } 
  
  my @multi_value_attributes = ();
  push @multi_value_attributes, "Remark";
  push @multi_value_attributes, "Parent Formula";
  push @multi_value_attributes, "BRENDA name";
  push @multi_value_attributes, "Parent Avg molecular weight";

  my @single_value_attributes = ();
  push @single_value_attributes, "Biosequence";
  push @single_value_attributes, "LIGAND AccNo";

  foreach $attr (@multi_value_attributes, @single_value_attributes) {
     $field_name{$attr} = $attr;
     $field_name{$attr} =~ s/\s/_/g;
  }

  ### read data
  while (<$in>) {
      chomp;
      if (/^SMILES\s*/) { ### end of record
	  my $smiles = $';

	  ### instantiate new object
	  warn ";creating object\t$current_id in class\t$class\n" 
	      if ($warn_level >= 2);
	  warn ("\t", $smiles, "\n")
 	      if ($warn_level >= 3);
	  $object = $class_holder->new_object();
	  $object->set_attribute("SMILES",$smiles);
	  $object->set_attribute("source","Glaxo");
	  
	  ### reeinitialize the arguments
	  my %attributes = ();
	  
	  ### read next line
	  next;
      }

      ### attributes
      my $found = 0;
      foreach $attr (@multi_value_attributes) {
	  if (/^\s+$attr\s*/) {
	      $object->new_attribute_value($field_name{$attr},$');
	      $found = 1;
	  }
      }
      foreach $attr (@single_value_attributes) {
	  if (/^\s+$attr\s*/) {
	      $object->set_attribute($field_name{$attr},$');
	      $found = 1;
	  }
      }
      next if ($found);
      if (/^\t\S+/) {
	  print ERR ("uknown attribute, not parsed\t",
		     $_, "\n");
      }

  }
  
  close $in if ($input_file);

  return $class_holder;
}
