#!/usr/bin/perl
############################################################
#
# $Id: parse_transfac.pl,v 1.3 2003/10/29 09:04:13 jvanheld Exp $
#
# Time-stamp: <2003-07-10 11:52:52 jvanheld>
#
############################################################

#use strict;
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`"); ### add the program's directory to the lib path
}
require "config.pl";
require "lib/load_classes.pl";
require "lib/util.pl";
require "lib/loading_util.pl"; ### for converting polypeptide IDs into ACs
require "lib/parsing_util.pl";



################################################################
## regulatory site
package classes::TransfacSite;
{
  @ISA = qw ( classes::DatabaseObject );
  ### class attributes
  $_count = 0;
  $_prefix = "site_";
  @_objects = ();
  %_name_index = ();
  %_id_index = ();
  %_attribute_count = ();
  %_attribute_cardinality = (id=>"SCALAR",           ### internal ID for the pathway diagram
			     transfac_id=>"SCALAR",      ### site ID from the Transfac entry
			     description=>"SCALAR",
			     gene_ac=>"SCALAR",      ### regulated gene AC
			     gene_name=>"SCALAR",    ### regulated gene name
#			     factor_ac=>"SCALAR",    ### transcription factor AC
#			     factor_name=>"SCALAR",  ### transcription factor name
#			     quality=>"SCALAR",      ### quality of the site-factor assignment
			     bound_factors=>"EXPANDED", 
			     site_sequence=>"SCALAR"      ### site sequence
			     );     
}


################################################################
## Transcription factor
package classes::TransfacFactor;
{
  @ISA = qw ( classes::DatabaseObject );
  ### class attributes
  $_count = 0;
  $_prefix = "fact_";
  @_objects = ();
  %_name_index = ();
  %_id_index = ();
  %_attribute_count = ();
  %_attribute_cardinality = (id=>"SCALAR",           ### internal ID for the pathway diagram
			     transfac_id=>"SCALAR",  ### transfac factor ID
			     organism=>"SCALAR",     ### OS field
			     name=>"SCALAR",         ### FA field
			     );     
}

################################################################
## Gene
package classes::TransfacGene;
{
  @ISA = qw ( classes::DatabaseObject );
  ### class attributes
  $_count = 0;
  $_prefix = "fact_";
  @_objects = ();
  %_name_index = ();
  %_id_index = ();
  %_attribute_count = ();
  %_attribute_cardinality = (id=>"SCALAR",           ### internal ID for the pathway diagram
			     transfac_id=>"SCALAR",  ### transfac gene ID
			     organism=>"SCALAR",     ### OS field
			     description=>"SCALAR",  ### DE field
			     );     
}


package main ;
{
    
    #### input directory
    $dir{transfac} = "${Databases}/transfac";


    #### default export directory
    $export_subdir = "transfac";
    $dir{output} = "$parsed_data/${export_subdir}/$delivery_date";

    #### check existence of input directories
    foreach $d (keys %dir) {
	unless (-d ($dir{$d}) ) {
	    die "Error: $d directory $dir{$d} does not exist\n";
	}
    }

    @files = qw (site factor gene);
#    @files = qw (site factor matrix class cell );
    foreach my $file (@files) {
	$file_name{$file} = "tf_$file.dat";
	$in_file{$file} = "$dir{transfac}/tf_$file.dat";
	unless (-e $in_file{$file}) {
	    die "Error: $file file does not exist\t$in_file{$file}\n";
	}
    }
    
    $verbose = 0;
    $out_format = "obj";

    #### classes and class_holders
    @classes = qw( classes::TransfacSite classes::TransfacFactor classes::TransfacGene );
    $sites = classes::ClassFactory->new_class(object_type=>"classes::TransfacSite");
    $factors = classes::ClassFactory->new_class(object_type=>"classes::TransfacFactor");
    $genes = classes::ClassFactory->new_class(object_type=>"classes::TransfacGene");


    &ReadArguments;

    #### output directory
    &CheckOutputDir();
    $out_file{error} = "$dir{output}/transfac.errors.txt";
    $out_file{stats} = "$dir{output}/transfac.stats.txt";
    $out_file{transfac} = "$dir{output}/transfac.obj" if ($export{obj});

    ### open error report file
    open ERR, ">$out_file{error}" || die "Error: cannot write error file $out_file{error}\n";


    ### test conditions
    if ($test) {
	warn ";TEST\n" if ($verbose >= 1);
	### fast partial parsing for debugging
	foreach $key (keys %in_file) {
	    $in_file{$key} = " head -5000 $in_file{$key} |";
	}
    }

    &DefaultVerbose if ($verbose >= 1);
 
    ### parse data from original files
    &ParseTransfacFile($in_file{site}, $sites);
    &ParseTransfacFile($in_file{factor}, $factors);
    &ParseTransfacFile($in_file{gene}, $genes);
    &TreatSites();
    &TreatFactors();
    &TreatGenes();


    ### print result
    &PrintStats($out_file{stats}, @classes);
    foreach $class_holder ($sites, $factors, $genes) {
	$class_holder->dump_tables();
	$class_holder->generate_sql(schema=>"jvanheld", 
				    dir=>"$dir{output}/sql_scripts", 
				    prefix=>"transfac_");
    }
    &ExportClasses($out_file{transfac}, $out_format, @classes) if $export{obj};


    ### report execution time
    if ($verbose >= 1) {
	$done_time = &AlphaDate;
	warn ";\n";
	warn "; job started $start_time";
	warn "; job done    $done_time\n";
    }

    close ERR;


    &deliver("transfac_parsed") if ($deliver);

    &CompressParsedData();

    exit(0);
}


################################################################
############## subroutines for the main package ################
################################################################

### print the help message 
### when the program is called with the -h or -help option
sub PrintHelp {
  open HELP, "| more";
  print <<EndHelp;
NAME
	parse_kegg.pl

DESCRIPTION
	Parse TRANSFAC regulatory sites.
	
AUTHOR
	Jacques van Helden (jvanheld\@ucmb.ulb.ac.be)  

OPTIONS	
	-h	detailed help
	-help	short list of options
	-test	fast parsing of partial data, for debugging
	-outdir output directory
	-v #	warn level
		Warn level 1 corresponds to a restricted verbose
		Warn level 2 reports all polypeptide instantiations
		Warn level 3 reports failing get_attribute()
	-obj    Export results in "obj" format (human readable)
	-clean	remove all files from the output directory before
		parsing
	-deliver 
		copy the relevant files in the delivery directory,
		for loading in the aMAZE database. 
EndHelp
  close HELP;
}

  



### read arguments from the command line
sub ReadArguments {
    for my $a (0..$#ARGV) {
	
	### warn level
	if (($ARGV[$a] eq "-v" ) && 
	    ($ARGV[$a+1] =~ /^\d+$/)){
	    $main::verbose = $ARGV[$a+1];
	    
	    ### test run
	} elsif ($ARGV[$a] eq "-test") {
	    $main::test = 1;
	    
	    ### clean
	} elsif ($ARGV[$a] eq "-clean") {
	    $main::clean = 1;
	    
	    ### deliver
	} elsif ($ARGV[$a] eq "-deliver") {
	    $main::deliver = 1;
	    

	    ### output dir
	} elsif ($ARGV[$a] eq "-outdir") {
	    $a++;
	    $main::dir{output} = $ARGV[$a];
	    
	    ### help
	} elsif (($ARGV[$a] eq "-h") ||
		 ($ARGV[$a] eq "-help")) {
	    &PrintHelp;
	    exit(0);

	    #### export object file
	} elsif ($ARGV[$a] =~ /^-obj/) {
	    $main::export{obj} = 1;
	}
	
    }
}




################################################################
# Parse TRANSFAC regulatory sites from sites.dat
#
sub ParseTransfacFile {
    my ($input_file, $class_holder, $source) = @_;
    
    warn (";\n; ", &AlphaDate,  "  parsing file $input_file\n")
	if ($verbose >= 1);

    open DATA, $input_file || 
	die "Error: cannot open data file $input_file\n";
    
    #### Read an entire record at a time
    $/ = "\/\/\n";

    my $entries = 0;
    while ($text_entry = <DATA>) {
	$entries++;
	warn "; $input_file\t$entries entries\n" if ($verbose >=2);

	if ($entries==1) {
	    #### check if there is a header
    	    if ($header =~ /TRANSFAC \S+ TABLE, V.(\S+)\s+(\S+)/i) {
		$transfac_version = $1;
		$release_date = $2;
		if ($verbose >= 1) {
		    warn "; Transfac version $transfac_version\n";
		    warn "; Release date     $release_date\n";
		}
		next;
	    } else {
		warn "Warning: this file does not contain the header of a TRANSFAC file\n";
	    }
	}

	unless ($source) {
	    $source = "Transfac";
	    $source .= " $transfac_version" if ($transfac_version);
	}

	#### parse the entry
	my @lines = split "\n", $text_entry;
	foreach my $line (@lines) {
	    chomp $line;
	    next unless ($line); # skip empty lines
	    next if ($line =~ /^XX$/); # skip separator lines
	    next if ($line =~ /^\/\/$/); # skip record separator
	    if ($line =~ /^(\S{2})  /) {
		my $key = $1;
		my $value = $';
		if ($key eq "AC") {
		    $current_obj = $class_holder->new_object(source=>$source, 
						      id=>$value, 
						      transfac_ac=>$value);
		} elsif ($key eq "ID") {
		    $current_obj->set_attribute("transfac_id", "\"$value\""); ### to prevent Oracle loading problems with the '$' in site IDs 
		} else {
		    $current_obj->new_attribute_value($key, $value);
		}
		
		$last_key = $key;
	    } else {
		&ErrorMessage("Unrecognized syntax, line skipped\t$line\n");
	    }
	}
    }
}


################################################################
# Treatment for sites
#
sub TreatSites {
    $sites->set_attribute_header("bound_factors", join("\t", "factor_ac", "factor_name", "quality"));
    $sites->set_attribute_header("modifications", join("\t", "date", "author"));

    my @out_fields = qw (id 
			 transfac_id
			 transfac_ac
			 gene_ac
			 gene_region
			 site_sequence
			 organism
			 denomination
			 description

			 bound_factors
			 modifications
			 );

    my @other_fields = qw (
			 matrices
			 );

    foreach my $site ($sites->get_objects()) {
	################################################################
	# specific treatment for some attributes
	
	#### description
	my $description = join (" ", $site->get_attribute("DE"));
	if ($description =~ /\; Gene\: (\S+)\./) {
	    $gene_ac = $1;
	} else {
	    &ErrorMessage (join ("\t", "Missing gene in site ",$site->get_attribute("ac") , $description), "\n");
	}
	$site->set_attribute("description", $description);
	$site->set_attribute("gene_ac", $gene_ac);
	
	#### date and author
	foreach my $dt ($site->get_attribute("DT")) {
	    my ($date, $author) = split "; ", $dt;
	    $site->push_expanded_attribute("modifications", $date, $author);
	}
	#### sequence
	my $sequence = join ("", $site->get_attribute("SQ")); ## join for possible multi-line sequence
	if ($sequence =~ /\.$/) {
	    $sequence = $`;
	} elsif ($sequence) {
	    &ErrorMessage (join ("\t", "unproperly terminated sequence (no dot)", $site->get_attribute("ac"), $sequence), "\n");
	}
	if ($sequence) {
	    $site->set_attribute("site_sequence", $sequence);
	}

	#### organism classification
	$site->set_attribute("organism_classification", join (" ", $site->get_attribute("OC"))); ## join for possible multi-line sequence

	#### bound factor
	my @BF = $site->get_attribute("BF");
	foreach my $BF (@BF) {
	    my ($factor_ac, $factor_name, $quality) = split ';', $BF;
	    $quality =~ s/Quality: //;
	    $site->push_expanded_attribute("bound_factors", $factor_ac, $factor_name, $quality);
	}

	#### published references
	my @medline_ids = $site->get_attribute("RX");
	my @medline_ids = $site->get_attribute("RX");

	################################################################
        # Give an explicit name to other attributes
	
	#### single-value attributes
	$site->set_attribute("organism", $site->get_attribute("OS"));
	$site->set_attribute("gene_region", $site->get_attribute("RE"));
	$site->set_attribute("denomination", $site->get_attribute("EL"));
	$site->set_attribute("first_position", $site->get_attribute("SF"));
	$site->set_attribute("last_position", $site->get_attribute("SL"));
	$site->set_attribute("position_ref", $site->get_attribute("S1"));

	#### multi-value attributes
	$site->set_array_attribute("deduced_matrices", $site->get_attribute("MX"));
	$site->set_array_attribute("factor_sources", $site->get_attribute("SO"));
	$site->set_array_attribute("comments", $site->get_attribute("CC"));
	$site->set_array_attribute("external_refs", $site->get_attribute("DR"));

    }

    $sites->set_out_fields(@out_fields);
}


################################################################
## post-parsing for genes (assign explicit field names)
sub TreatGenes {
    foreach my $gene ($genes->get_objects()) {
	
        ################################################################
        # Give an explicit name to other attributes
	
	#### single-value attributes
	$gene->set_attribute("organism", $gene->get_attribute("OS"));
	$gene->set_attribute("description", $gene->get_attribute("DE"));
    }
}

################################################################
## post-parsing for factors (assign explicit field names)
sub TreatFactors {
    foreach my $factor ($factors->get_objects()) {
	
        ################################################################
        # Give an explicit name to other attributes
	
	#### single-value attributes
	$factor->set_attribute("organism", $factor->get_attribute("OS"));
	$factor->set_attribute("name", $factor->get_attribute("FA"));
    }
}
