#!/usr/bin/perl
############################################################
#
# $Id: parse_transfac.pl,v 1.4 2004/07/06 17:45:47 jvanheld Exp $
#
# Time-stamp: <2003-07-10 11:52:52 jvanheld>
#
############################################################

#use strict;
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`"); ### add the program directory to the lib path
}
require "config.pl";
require "lib/load_classes.pl";
require "lib/util.pl";
require "lib/loading_util.pl"; ### for converting polypeptide IDs into ACs
require "lib/parsing_util.pl";



################################################################
## regulatory site
package classes::Site;
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
			     organism => "SCALAR",   ## OS field
#			     gene_name=>"SCALAR",    ### regulated gene name
#			     factor_ac=>"SCALAR",    ### transcription factor AC
#			     factor_name=>"SCALAR",  ### transcription factor name
#			     quality=>"SCALAR",      ### quality of the site-factor assignment
			     binding_factors_expanded=>"EXPANDED", 
			     site_sequence=>"SCALAR",      ### site sequence

			     type => "SCALAR", ## 
			     position_first => "SCALAR",
			     position_last => "SCALAR",
			     position_reference => "SCALAR",
			     organism_classification => "ARRAY",
			     gene_region => "SCALAR",
			     element => "SCALAR",
			     );     
}


################################################################
## Transcription factor
package classes::Factor;
{
  @ISA = qw ( classes::DatabaseObject );
  ### class attributes
  $_count = 0;
  $_prefix = "fact_";
  @_objects = ();
  %_name_index = ();
  %_id_index = ();
  %_attribute_count = ();
  %_attribute_cardinality = (id=>"SCALAR",              ## internal ID for the pathway diagram
			     transfac_id=>"SCALAR",     ## transfac factor ID
			     organism=>"SCALAR",        ## OS field
			     name=>"SCALAR",            ## FA field
			     gene_name => "SCALAR",     ## GE field
			     factor_class => "SCALAR",  ## CL field

			     names=>"ARRAY", 
			     );     
}

################################################################
## Gene
package classes::Gene;
{
  @ISA = qw ( classes::DatabaseObject );
  ### class attributes
  $_count = 0;
  $_prefix = "fact_";
  @_objects = ();
  %_name_index = ();
  %_id_index = ();
  %_attribute_count = ();
  %_attribute_cardinality = (id=>"SCALAR",           ## internal ID for the pathway diagram

			     transfac_id=>"SCALAR",  ## transfac gene ID
			     organism=>"SCALAR",     ## OS field
			     description=>"SCALAR",  ## DE field
			     short_description => "SCALAR", ## SD field
			     EPD_classification => "SCALAR", ## BC field
			     
			     taxonomy => "SCALAR",    ## concatanation of the OC field
			     names=>"ARRAY", 
			    );     
}

################################################################
#### Main package
package main ;
{
    $host= $default{'host'};
    $schema="transfac";
    $user="transfac";
    $password="transfac";
#    $special_field_size{description} = 10000;

    %key_alias = ();
    %transfac_index = (); ## Index transfac objects with transfac AC as key

    
    #### input directory
    $dir{transfac} = "${Databases}/transfac";

    #### default export directory
    $export_subdir = "transfac";
    $dir{output} = "$parsed_data/${export_subdir}/$delivery_date";

    @files = qw (site factor gene);
#    @files = qw (site factor matrix class cell );

    $verbose = 0;
    $out_format = "obj";

    #### classes and class_holders
    @classes = qw( classes::Site classes::Factor classes::Gene );
    $sites = classes::ClassFactory->new_class(object_type=>"classes::Site");
    $factors = classes::ClassFactory->new_class(object_type=>"classes::Factor");
    $genes = classes::ClassFactory->new_class(object_type=>"classes::Gene");

    $sites->set_attribute_header("modifications", join("\t", "date", "author"));
    $factors->set_attribute_header("modifications", join("\t", "date", "author"));
    $genes->set_attribute_header("modifications", join("\t", "date", "author"));
    $genes->set_attribute_header("xrefs", join("\t", "xdb", "xid"));

    &ReadArguments();

    #### output directory
    &CheckOutputDir();

    #### check existence of input and output directories
    foreach $d (keys %dir) {
	warn "Checking existence of $d directory $dir{$d}\n" if ($main::verbose >= 1);
	unless (-d ($dir{$d}) ) {
	    die "Error: $d directory $dir{$d} does not exist\n";
	}
    }

    #### Check existence of input files
    foreach my $file (@files) {
	$in_file{$file} = $dir{transfac}."/".$file.".dat";
	unless (-e $in_file{$file}) {
	    die "Error: $file file does not exist\t$in_file{$file}\n";
	}
    }
     
    ### Report files
    $out_file{error} = "$dir{output}/transfac.errors.txt";
    open ERR, ">$out_file{error}" || die "Error: cannot write error file $out_file{error}\n";
    $out_file{stats} = "$dir{output}/transfac.stats.txt";
    $out_file{transfac} = "$dir{output}/transfac.obj" if ($export{obj});

    ### test conditions
    if ($test) {
	warn ";TEST\n" if ($verbose >= 1);
	### fast partial parsing for debugging
	foreach $key (keys %in_file) {
	    $in_file{$key} = " head -5000 $in_file{$key} |";
	}
    }

    &DefaultVerbose() if ($verbose >= 1);
 
    ### parse data from original files
    &ParseTransfacFile($in_file{site}, $sites);
    &ParseTransfacFile($in_file{factor}, $factors);
    &ParseTransfacFile($in_file{gene}, $genes);

    &TreatSites();
    &TreatFactors();
    &TreatGenes();

    ### print result
    &PrintStats($out_file{stats}, @classes);
    foreach my $class_factory ($sites, $factors, $genes) {
	warn "; Dumping class $factory_name $class_factory\n" if ($verbose >= 1);
	$class_factory->dump_tables();
	$class_factory->generate_sql(dir=>"$dir{output}/sql_scripts",
				     prefix=>"$class_factory_",
				     schema=>$schema,
				     host=>$host,
				     dbms=>$dbms,
				     user=>$user,
				     password=>$password
				    );
    }
    &ExportMakefile(@classes);

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

#    &CompressParsedData();
    
    warn "; Result saved in directory $dir{output}\n" if ($main::verbose >= 1);

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
	-indir  input directory. 
		This directory contains the TRANSFAC flat files. 
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
   Options for the automaticaly generated SQL scripts
	-schema database schema (default: $schema)
	-host	database host (efault: $host)
	-user	database user (efault: $user)
	-password	
		database password (default: $password)
EndHelp
  close HELP;
  exit(0);
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
	    
	    ### input dir
	} elsif ($ARGV[$a] eq "-indir") {
	    $a++;
	    $main::dir{transfac} = $ARGV[$a];
	    
	    ### help
	} elsif (($ARGV[$a] eq "-h") ||
		 ($ARGV[$a] eq "-help")) {
	    &PrintHelp;
	    exit(0);

	    #### export object file
	} elsif ($ARGV[$a] =~ /^-obj/) {
	    $main::export{obj} = 1;

	    ### database schema
	} elsif ($ARGV[$a] eq "-schema") {
	    $schema = $ARGV[$a+1];
	    
	    ### database host
	} elsif ($ARGV[$a] eq "-host") {
	    $host = $ARGV[$a+1];
	    
	    ### database user
	} elsif ($ARGV[$a] eq "-user") {
	    $user = $ARGV[$a+1];
	    
	    ### password 
	} elsif ($ARGV[$a] eq "-password") {
	    $password = $ARGV[$a+1];
	    
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
    	    if ($text_entry =~ /TRANSFAC \S+ TABLE/i) {
		if ($' =~ /Release (\S+)/) { 
		    $transfac_version = $1;
		   }
		if ($' =~ /(\d+\-\d+\-\d+)/) { 
		    $release_date = $1;
		}
		if ($verbose >= 1) {
		    warn "; Transfac version $transfac_version\n";
		    warn "; Release date     $release_date\n";
		}
		next;
	    } else {
		warn "Warning: file $input_file does not contain the header of a TRANSFAC file\n";
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
		my $value = $'; #';

		################################################################
		#### Create new object 
		if ($key eq "AC") {
		    $current_obj = $class_holder->new_object(source=>$source, 
							     id=>$value, 
							     transfac_ac=>$value);
		    $transfac_index{$value} = $current_obj;
		    next;
		}

		## Replace some field names by aliases to avoid
		## conflicts with reserved SQL words (e.g. ID, IN)

		## Aliases for attributes common to several types of objects
		$key_alias{CO} = "copyrigth";
		$key_alias{ID} = "transfac_id";
		$key_alias{DT} = "date";
		$key_alias{SD} = "short_description";
		$key_alias{DE} = "descr";
		$key_alias{OS} = "organism";
		$key_alias{OC} = "organism_classification";
		$key_alias{DR} = "database_references";
		$key_alias{RN} = "reference_number";
		$key_alias{RX} = "reference_xref";
		$key_alias{RA} = "reference_authors";
		$key_alias{RT} = "reference_title";
		$key_alias{RL} = "reference_source";
		$key_alias{FA} = "factor_name";
		$key_alias{BS} = "binding_sites";
		$key_alias{MX} = "matrices";

		## Gene-specific aliases
		$key_alias{SY} = "synonyms";
		$key_alias{CH} = "chromosomal_location";
		$key_alias{BC} = "EPD_classification";
		$key_alias{RG} = "regulation";
		$key_alias{CE} = "composite_elements";

		## Factor-specific aliases
		$key_alias{GE} = "gene_name";
		$key_alias{HO} = "homologs";
		$key_alias{CL} = "factor_class";
		$key_alias{SZ} = "size";
		$key_alias{SQ} = "seq";
		$key_alias{SC} = "sequence_source";
		$key_alias{FT} = "feature_table";
		$key_alias{SF} = "structural_features";
		$key_alias{CP} = "cell_specificity_positive";
		$key_alias{CN} = "cell_specificity_negative";
		$key_alias{EX} = "expression_patttern";
		$key_alias{FF} = "functional_properties";
		$key_alias{IN} = "interacting_factor";

		## Site-specific aliases
		$key_alias{TY} = "sequence_type";
		$key_alias{RE} = "gene_region";
		$key_alias{EL} = "element";
		$key_alias{SF} = "position_first";
		$key_alias{ST} = "position_last";
		$key_alias{S1} = "position_reference";
		$key_alias{BF} = "binding_factors";
		$key_alias{SO} = "cellular_factor_source";
		$key_alias{MM} = "methods";
		$key_alias{CC} = "comments";
#		$key_alias{} = "";

		if (defined($key_alias{$key})) {
		    $key = $key_alias{$key};
		}
		$current_obj->new_attribute_value($key, $value);
	    } else {
		&ErrorMessage("Unrecognized syntax, line skipped\t$line\n");
	    }
	}
	&treat_generic_attributes($current_obj);
    }
}


################################################################
#### Treat some generic attributes
sub treat_generic_attributes {
    my ($current_obj) = @_;

    
    #### date and author
    foreach my $dt ($current_obj->get_attribute($key_alias{"DT"})) {
	my ($date, $author) = split "; ", $dt;
	$current_obj->push_expanded_attribute("modifications", $date, $author);
    }
    
    #### organism classification
    my $taxonomy = &join_attribute($current_obj, "organism_classification", " ", "taxonomy");
#    &join_$current_obj->set_attribute("taxonomy", join (" ", $current_obj->get_attribute("OC"))); ## join for possible multi-line sequence
    
    #### Add AC as name, in order to be able to retrieve it
#    $current_obj->push_attribute("names", $current_obj->get_attribute("transfac_ac"));
#    $current_obj->index_object_names();
#    $current_obj->set_attribute("name", $current_obj->get_name());

    #### published references
#    my @medline_ids = $current_obj->get_attribute("RX");
#    my @medline_ids = $current_obj->get_attribute("RX");
}   


################################################################
## Parse database references, and create a tab-delimited table with
## one column for the external database and one column ffor the
## external ID
sub parse_database_references {
    my ($class) = @_;
    
    foreach my $object ($class->get_objects()) {
	foreach my $db_ref ($object->get_attribute("database_references")) {
	    my @fields = split ":", $db_ref;
	    my $xdb = &trim(shift @fields);
	    my $xids = &trim(shift @fields);
	    $xids =~ s/\.$//; 
	    my @xids = split ";", $xids;
	    foreach my $xid (@xids) {
		$current_obj->push_expanded_attribute("xrefs", $xdb, &trim($xid));
	    }
	}
    }
}
    
################################################################
# Specific treatment for sites
#
sub TreatSites {
    $sites->set_attribute_header("binding_factors_expanded", join("\t", "factor_ac", "factor_name", "quality"));

    foreach my $site ($sites->get_objects()) {
	#### Description
	my $description = &join_attribute($site, "descr", "", "description");

	################################################################
	#### extract gene AC from description
	if ($description =~ /^consensus./i) {
	    $site->set_attribute("type", "consensus");
	} elsif ($description =~ /^artificial/i) {
	    $site->set_attribute("type", "artificial");
	} elsif ($description =~ /Gene\:\s*(\S+)\./) {
	    $site->set_attribute("type", "cis-acting element");
	    my $gene_ac = $1;
	    $site->set_attribute("gene_ac", $gene_ac);
	    if (defined($transfac_index{$gene_ac})) {
		my $gene_obj = $transfac_index{$gene_ac};
		my $gene_name = join ("", $gene_obj->get_attribute($key_alias{SD}));
		$site->set_attribute("gene_name", $gene_name);
		warn join ("\t", "; Site",  
			   $site->get_attribute("id"),
			   "Gene AC", $gene_ac,
			   "Gene name", $gene_name,
			  ), "\n" if($main::verbose >= 10);
	    } else {
		&ErrorMessage(join("\t", "; Site", 
				   $site->get_attribute('id'),
				   "Cannot identify gene", $gene_ac));
			     }
	} else {
	    &ErrorMessage (join ("\t", "Missing gene in site ",$site->get_attribute("ac") , $description), "\n");
	}
	
	#### sequence
	my $sequence = join ("", $site->get_attribute($key_alias{"SQ"})); ## join for possible multi-line sequence
	if ($sequence =~ /\.$/) {
	    $sequence = $`;
	} elsif ($sequence) {
	    &ErrorMessage (join ("\t", "unproperly terminated sequence (no dot)", $site->get_attribute("ac"), $sequence), "\n");
	}
	if ($sequence) {
	    $site->set_attribute("site_sequence", $sequence);
	}

	#### bound factor
	my @BF = $site->get_attribute($key_alias{"BF"});
	foreach my $BF (@BF) {
	    my ($factor_ac, $factor_name, $quality) = split ';', $BF;
	    $quality =~ s/Quality: //;
	    $site->push_expanded_attribute("binding_factors_expanded", $factor_ac, $factor_name, $quality);
	}


    }
}


################################################################
## post-parsing for genes (assign explicit field names)
sub TreatGenes {
    foreach my $gene ($genes->get_objects()) {

	#### Use short description as primary name
	$gene->set_attribute("name", $gene->get_attribute("short_description"));
	$gene->push_attribute("names", $gene->get_attribute("short_description"));

	#### Use description as name
	my $description = &join_attribute($gene, $key_alias{DE}, "", "description");
	$gene->push_attribute("names", $description);

	warn join ("\t", $gene, $gene->get_attribute('id'), 
		  $description,
		  join (";", $gene->get_attribute("names")),
		 ),"\n" if ($main::verbose >= 10);

    }

    &parse_database_references($genes);
}

################################################################
## post-parsing for factors (assign explicit field names)
sub TreatFactors {
    foreach my $factor ($factors->get_objects()) {
	my $primary_name = $factor->get_attribute($key_alias{FA});
	$factor->set_attribute("name", $primary_name);
	$factor->push_attribute("names", $primary_name);
    }

}


################################################################
## Join attribute fragmented on several lines
## Create a new single-value attribute with the result
sub join_attribute {
    my ($object, $attribute_name, $separator, $join_attribute_name) = @_;
    my @values = $object->get_attribute($attribute_name);
    my $value = join($separator, @values);
    warn join ("\t", 
	      "; Joining attribute", 
	      $object->get_attribute('id'), 
	      $attribute_name, 
	      $join_attribute_name, 
	      join (";", @values),
	      "'".$value."'"), "\n" if ($main::verbose >= 10);

    if ($join_attribute_name) {
	warn join ("\t", "; Creating new attribute", $join_attribute_name, $value), "\n" if ($main::verbose >= 4);
	$object->force_attribute($join_attribute_name, $value);
    }
    return $value;
}
