#!/usr/bin/perl
############################################################
#
# $Id: parse_swissprot.pl,v 1.2 2000/11/29 14:00:48 jvanheld Exp $
#
# Time-stamp: <2000-11-29 10:39:55 jvanheld>
#
############################################################

### add the program's directory to the lib path
BEGIN {
    if ($0 =~ /([^(\/)]+)$/) {
	$script_dir = $` || ".";
	push @INC, "$script_dir";
	push @INC, "$script_dir/perllib/";
    }
}
require "PFBP_config.pl";
require "PFBP_classes.pl";
require "PFBP_util.pl";
require "PFBP_parsing_util.pl";
require "PFBP_loading_util.pl";
use SWISS::Entry;

package main;
{


  ### organism selection
  @selected_organisms = ();
  $full_name{yeast} = "Saccharomyces cerevisiae (Baker's yeast)";
  $full_name{ecoli} = "Escherichia coli";
  $full_name{human} = "Homo sapiens (Human)";

  ### input directories and files
  $dir{swissprot} = "/win/databases/downloads/ftp.ebi.ac.uk/pub/databases/swissprot/release_compressed/";
  $source{swissprot} = "sprot39";
  
  $dir{trembl} = "/win/databases/downloads/ftp.ebi.ac.uk/pub/databases/trembl/";
  $source{yeast} = "sptrembl/fun";
  $source{human} = "sptrembl/hum";
  $source{ecoli} = "sptrembl/pro";
#  &OpenIndex("ECSet");
  
  ### export classes
  @classes = qw( PFBP::Polypeptide PFBP::Catalysis );
  @{$out_fields{'PFBP::Polypeptide'}} = qw( id name description gene organisms swissprot_acs swissprot_ids names );
  @{$out_fields{'PFBP::Catalysis'}} = qw( id catalyst catalyzed source );
  
  &ReadArguments;
  
  #### output directory
  $out_format = "obj";
  $dir = $parsed_data."/swissprot_parsed";
  unless (-d $dir) {
      warn "Creating output dir $dir", "\n";
      mkdir $dir, 0775 || die "Error: cannot create directory $dir\n";
  }
  $dir{output} = $dir."/".$delivery_date;
  unless (-d $dir{output}) {
      warn "Creating output dir $dir{output}\n";
      mkdir $dir{output}, 0775 || die "Error: cannot create directory $dir\n";
  }
  die unless chdir $dir{output};
  
  #### output file names depend on the organisms to parse and source 
  #### -> defined after reading arguments
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
  
  $out_file{polypeptides} = $dir{output}."/Polypeptide".$suffix.".obj";
##  $out_file{mldbm_polypeptides} = $dir{output}."/Polypeptide".$suffix.".mldbm";
#  $out_file{mldbm_polypeptide_index} = $dir{output}."/Polypeptide".$suffix."__name_index.mldbm";
  $out_file{catalyses} = $dir{output}."/Catalysis".$suffix.".obj";
#  $out_file{mldbm_catalyses} = $dir{output}."/Catalysis".$suffix.".mldbm";
  $out_file{errors} = $dir{output}."/swissprot".$suffix.".errors.txt";
  $out_file{stats} = $dir{output}."/swissprot".$suffix.".stats.txt";
  
  ### open error report file
  open ERR, ">$out_file{errors}" || 
      die "Error: cannot write error file $out_file{errors}\n";
  
  ### select all organisms if none was selected (-org)
  unless ($#selected_organisms >= 0) {
    push @selected_organisms, "ecoli";
    push @selected_organisms, "human";
    push @selected_organisms, "yeast";
  }

  ### define regular expressions to match the selected organisms
  foreach $organism (@selected_organisms) {
      $regexp = $full_name{$organism};
      $regexp = $` if ($regexp =~ /\(/); 
      $regexp = $` if ($regexp =~ /\)/);
      $regexp = $` if ($regexp =~ /\'/);
      $selected_organism{uc($full_name{$organism})} = $regexp;
      push @regexps, qr /$regexp/; ### list of precompiled regular expressions
  }
  

  ### by default, parse swissprot only
  unless (defined(%data_sources)) {
      $data_sources{swissprot} = 1;
  }
  if ($data_sources{swissprot}) {
      $in_file{$source{swissprot}} = "uncompress -c ".$dir{swissprot}."/".$source{swissprot}.".dat.Z | ";
      $files_to_parse{$source{swissprot}} = 1;
  }
  if ($data_sources{trembl}) {
      foreach $organism (@selected_organisms) {
	  $in_file{${source{$organism}}} = "uncompress -c ".$dir{trembl}."/".$source{$organism}.".dat.Z | ";
	  $files_to_parse{$source{$organism}} = 1;
      }
  }

  ### select all polypeptide if no specific class was selected (-enz)
  unless (defined(%export)) {
      $export{all} = 1;
  }
  
  ### create class holders
  $polypeptides = PFBP::ClassFactory->new_class(object_type=>"PFBP::Polypeptide",
						prefix=>"spp_");
  $catalyses = PFBP::ClassFactory->new_class(object_type=>"PFBP::Catalysis",
					     prefix=>"sct_");
  
  ### testing mode
  if ($test) {
    warn ";TEST\n" if ($warn_level >= 1);
    ### fast partial parsing for debugging
    foreach $key (keys %in_file) {
	$in_file{$key} .= " head -100000 |";
    }
}
  
  
  ### default verbose message
  if ($warn_level >= 1) {
      warn "; Selected organisms\n;\t", join("\n;\t", @selected_organisms), "\n";
      warn "; Polypeptide classes\n;\t", join("\n;\t", keys %export), "\n";
      warn "; Data sources\n;\t", join("\n;\t",  keys %data_sources), "\n";
      &DefaultVerbose;
  }
  
  
  ### parse data from original files
  foreach $source (keys %files_to_parse) {
      &ParseSwissprot($in_file{$source}, $source);
  }
#  &CreateCatalyses;

  #### use first name as primary name
  foreach $polypeptide ($polypeptides->get_objects()) {
      if ($name = $polypeptide->get_name()) {
	  $polypeptide->set_attribute("primary_name",$name);
      }
  }
  

  #### print the result
  &ExportClasses($out_file{polypeptides}, $out_format,PFBP::Polypeptide);
  &ExportClasses($out_file{catalyses}, $out_format, PFBP::Catalysis);
  $polypeptides->dump_tables($suffix);
  $catalyses->dump_tables($suffix);
#  $catalyses->export('MLDBM',$out_file{mldbm_catalyses});
#  $polypeptides->export("MLDBM",$out_file{mldbm_polypeptides});
#  $polypeptides->export_name_index("MLDBM",$out_file{mldbm_polypeptide_index});
  
  &PrintStats($out_file{stats}, @classes);
  
  #if ($out_file{stats}) {
  #  open STDOUT, ">>$out_file{stats}"
  #    || die "Error: cannot write stat file $out_file{stats}\n";
  #  &SwissStats;
  #  close $out if ($out_file{stats});
  #}
  
  
  ### report execution time
  if ($warn_level >= 1) {
      $done_time = &AlphaDate;
      warn (";\n",
	    "; job started $start_time",
	    "; job done    $done_time\n")
      }
  
  &CloseIndex("ECSet");
  close ERR;
  
  system "gzip -f $dir{output}/*.tab $dir{output}/*.obj $dir{output}/*.txt";
  
  exit(0);
}

### subroutines for the main package

### print the help message 
### when the program is called with the -h or -help option
sub PrintHelp {
  open HELP, "| more";
  print <<EndHelp;
NAME
	parse_swissprot.pl

DESCRIPTION
	Parse polypeptides from a Swissprot flat file

AUTHOR
	Jacques van Helden (jvanheld\@ucmb.ulb.ac.be)  

VERSION
	0.01
	Created		2000/01/11
	Last modified	2000/01/11
	
SYNOPSIS	
	parse_swissprot.pl [-v] [-vv] [-i infile] [-o outfile] 
		[-w warn_level]


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
	-org	select an organism for exportation
		can be used reiteratively in the command line 
		to select several organisms
		   Supported organisms :
			ecoli
			human
			yeast
		by default, all organisms are selected
	-data	select and organism for exportation
		   Valid data sources:
			swissprot
			trembl

		by default, all sources are selected
EXAMPLE
	parse_polypeptides.pl -w 2 -org ecoli -data swissprot -enz
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
    } elsif ($ARGV[$a] =~ /^-data/) {
      $data_sources{lc($ARGV[$a+1])} = 1;


      ### select enzymes for exportation
    } elsif ($ARGV[$a] =~ /^-enz/) {
      $main::export{enzymes} = 1;
    } elsif ($ARGV[$a] =~ /^-all/) {
      $main::export{all} = 1;
    }


  }

}

#### this is the routine for actually parsing the swissprot or trembl file
sub ParseSwissprot {
    my ($input_file, $source) = @_;
    $source = $input_file unless ($source);
    my %attrib_keys = (
		       id=>"ID",
		       ac=>"AC",
		       gn=>"Gene"
		       );
    my $class = PFBP::Polypeptide;
    
    warn (";\n; ", &AlphaDate,  " parsing polypeptides from $input_file\n")
	if ($warn_level >= 1);
    
    open DATA, $input_file || 
	die "Error: cannot open data file $input_file\n";
    
    # Read an entire record at a time
    $/ = "\/\/\n";
    
    my $entries = 0;
    while ($text_entry = <DATA>){
	$entries++;
	
	my $export = 0;
	
	### check that the entry matches organism name 
	### before converting it to an object
	my $parse = 0;
	foreach $regexp (@regexps) {
	    if ($text_entry =~ /$regexp/i) {
		$parse = 1;
		last;
	    }
	}
	next unless $parse;
	
	# Read the entry
	my $object_entry = SWISS::Entry->fromText($text_entry);
	
	### check whether the polypeptide organism(s)
	### match the organism selection
	my @organisms = $object_entry->OSs->elements;
	foreach $organism (@organisms) {
	    if ($selected_organism{uc($organism)}) {
		$export = 1;
		last;
	    }
	}
	next unless $export;
	
	### resets the export flag to 0
	$export = 0 unless ($export{all});
	
	# get the polypeptide attributes
	my $swissprot_ac = $object_entry->AC;
	my @swissprot_ids = $object_entry->IDs->elements;
	my @swissprot_acs = $object_entry->ACs->elements;
	my $descr = $object_entry->DEs->text;
	my $geneNames = $object_entry->GNs->text;
	my @geneNames = split " OR ", $geneNames;
	### extract name from description
	my @names = ();
	if ($descr =~ / \(/) {
	    push @names, lc($`);
	    #if ($descr =~ /^([^( \()]+)/) {
	    #push @names, lc($1);
	} else {
	    push @names, lc($descr);
	}
	push @names, @geneNames;
	
	### extract ECs from description
	my @ECs = ();
	my $tmp = $descr;
	while ($tmp =~ /\(EC ([^\(]+)\)/) {
	    push @ECs, $1;
	    $tmp = $';
	}
	if (($export{enzymes})  &&
	    ($#ECs >= 0)){
	    $export = 1; ### select enzyme for export
	}
	
	next unless ($export);
	
	### create a new polypeptide
	warn "$source\tentry $entries\t$swissprot_ids[0]\t$names[0]\n"
	    if ($warn_level >= 2); 
	my $polypeptide = $polypeptides->new_object(id=>$swissprot_ac,
						    source=>$source);
	if (defined($geneNames[0])) {
	    $polypeptide->set_attribute("gene",$geneNames[0]);
	} else {
	    $polypeptide->set_attribute("gene","<NULL>");
	}

	foreach my $name (@names) {
	    $polypeptide->push_attribute("names",$name);
	}
	foreach my $id (@swissprot_ids) {
	    $polypeptide->push_attribute("swissprot_ids",$id);
	    $polypeptide->push_attribute("names",$id);
	}
	foreach my $ac (@swissprot_acs) {
	    $polypeptide->push_attribute("swissprot_acs",$ac);
	    $polypeptide->push_attribute("names",$ac);
	}
	$polypeptide->set_attribute("description", $descr);

	my $pp_id = $polypeptide->get_id();
	foreach my $ec (@ECs) {
	    my $catalysis = $catalyses->new_object(source=>      $source,
						   catalyst=>    $pp_id,
						   catalyzed=>   $ec);
#	    $polypeptide->push_attribute("ECs",$EC);
	}
	foreach my $organism (@organisms) {
	    $polypeptide->push_attribute("organisms",$organism);
	}
	
	
    }
    close DATA;
    
    $polypeptides->index_names();
    return 1;
}


#### create catalysis objects from the EC numbers found in the description field
#### Note that these EC numbers are not cross-validated with those from KEGG (entered in teh DB)
#### this cross-validation has tom be done in the oracle database itself
#  sub CreateCatalyses {
#      foreach my $polypeptide (PFBP::Polypeptide->get_objects()) {
#  	my $pp_id = $polypeptide->get_attribute("id");
#  	if (my @ECs = $polypeptide->get_attribute("ECs")) {
#  	    foreach my $ec (@ECs) {
#  #		unless (defined($ECSet_index{$ec})) {
#  #		    $error_message = "Error: in $pp_id\tEC number $ec is not defined\n";
#  #		    print ERR $error_message;
#  #		    warn $error_message if ($warn_level >= 2);
#  #		}
#  		my $source = $polypeptide->get_attribute("source");
#  		$catalysis = $catalyses->new_object(source=>      $source,
#  						    catalyst=>    $pp_id,
#  						    catalyzed=>   $ec);
#  	    }
#  	}
#      }
#  }

