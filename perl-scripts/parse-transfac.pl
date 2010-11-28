#!/usr/bin/perl
############################################################
#
# $Id: parse-transfac.pl,v 1.14 2010/11/28 02:14:58 rsat Exp $
#
# Time-stamp: <2003-07-10 11:52:52 jvanheld>
#
############################################################

#use strict;
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}
require "RSA.lib";
require "RSA.seq.lib";
push @INC, "$ENV{RSAT}/perl-scripts/parsers/";
#require "config.pl";
require "lib/load_classes.pl";
require "lib/util.pl";
require "lib/loading_util.pl"; ### for converting polypeptide IDs into ACs
require "lib/parsing_util.pl";
require RSAT::util;
require RSAT::feature;

################################################################
## Matrix
package TRANSFAC::Matrix;
{
  @ISA = qw ( classes::DatabaseObject );
  ### class attributes
  $_count = 0;
  $_prefix = "matrix_";
  @_objects = ();
  %_name_index = ();
  %_id_index = ();
  %_attribute_count = ();
  %_attribute_cardinality = (id=>"SCALAR",           ### internal ID 
			     transfac_id=>"SCALAR",      ### site ID from the Transfac entry
			     description=>"SCALAR",
			     );     
}


################################################################
## regulatory site
package TRANSFAC::Site;
{
  @ISA = qw ( classes::DatabaseObject );
  ### class attributes
  $_count = 0;
  $_prefix = "site_";
  @_objects = ();
  %_name_index = ();
  %_id_index = ();
  %_attribute_count = ();
  %_attribute_cardinality = (id=>"SCALAR",           ### internal ID
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
package TRANSFAC::Factor;
{
  @ISA = qw ( classes::DatabaseObject );
  ### class attributes
  $_count = 0;
  $_prefix = "fact_";
  @_objects = ();
  %_name_index = ();
  %_id_index = ();
  %_attribute_count = ();
  %_attribute_cardinality = (id=>"SCALAR",              ## internal ID
			     transfac_id=>"SCALAR",     ## transfac factor ID
			     organism=>"ARRAY",        ## OS field
			     name=>"SCALAR",            ## FA field
			     gene_name => "ARRAY",     ## GE field became an array from version 11.2
			     factor_class => "SCALAR",  ## CL field

			     names=>"ARRAY",
			     );
}

################################################################
## Gene
package TRANSFAC::Gene;
{
  @ISA = qw ( classes::DatabaseObject );
  ### class attributes
  $_count = 0;
  $_prefix = "fact_";
  @_objects = ();
  %_name_index = ();
  %_id_index = ();
  %_attribute_count = ();
  %_attribute_cardinality = (id=>"SCALAR",           ## internal ID

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
## TransPro
package TRANSFAC::TransPro;
{
  @ISA = qw ( classes::DatabaseObject );
  ### class attributes
  $_count = 0;
  $_prefix = "prom_";
  @_objects = ();
  %_name_index = ();
  %_id_index = ();
  %_attribute_count = ();
  %_attribute_cardinality = (id=>"SCALAR",           ### internal ID 
			     transfac_id=>"SCALAR",      ### site ID from the Transfac entry
			     descr=>"SCALAR",
			     organism=>"SCALAR",
			     chromosomal_location=>"SCALAR",
			     );
}


################################################################
#### Main package
package main;
{

    ## Initialize parameters

#    $special_field_size{description} = 10000;
    $parse_transpro = 0;
    $host= $default{'host'};
    $schema="transfac";
    $user="transfac";
    $password="transfac";
    $pssm_sep="\t";
    %key_alias = ();
    %transfac_index = (); ## Index transfac objects with transfac AC as key
    $verbose = 0;
    $out_format = "obj";

    ## Read arguments
    &ReadArguments();

    #### default input and output directories
    if ($parse_transpro) {
	unless ($dir{transfac}) {$dir{transfac} = "${Databases}/transpro"};
	$export_subdir = "transpro";
    } else {
	unless ($dir{transfac}) {$dir{transfac} = "${Databases}/transfac"};
	$export_subdir = "transfac";
    }
    unless ($dir{output}) { $dir{output} = join('/', $parsed_data, ${export_subdir},$delivery_date) };

    #### classes and class_holders
    @data_types = ();
    @class_factories = ();
    @classes = ();
    if ($parse_transpro) {
	&RSAT::message::Info("Parsing TRANSPRO") if ($main::verbose >= 1);
	$transpro_holder = classes::ClassFactory->new_class(object_type=>"TRANSFAC::TransPro");
	$transpro_holder->set_attribute_header("modifications", join("\t", "date", "author"));

	## Create a class holder for parsed features
	$feature_holder = classes::ClassFactory->new_class(object_type=>"RSAT::feature");

	push @data_types, qw(transpro);
	push @classes, qw(TRANSFAC::TransPro);
	push @class_factories, ($transpro_holder, $feature_holder);

    } else {	
	$site_holder = classes::ClassFactory->new_class(object_type=>"TRANSFAC::Site");
	$site_holder->set_attribute_header("modifications", join("\t", "date", "author"));

	$factor_holder = classes::ClassFactory->new_class(object_type=>"TRANSFAC::Factor");
	$factor_holder->set_attribute_header("modifications", join("\t", "date", "author"));

	$gene_holder = classes::ClassFactory->new_class(object_type=>"TRANSFAC::Gene");
	$gene_holder->set_attribute_header("modifications", join("\t", "date", "author"));
	$gene_holder->set_attribute_header("xrefs", join("\t", "xdb", "xid"));

	$matrix_holder = classes::ClassFactory->new_class(object_type=>"TRANSFAC::Matrix");
	$matrix_holder->set_attribute_header("pssm", join("\t", "pos", "A", "C", "G", "T", "consensus"));

	push @data_types, qw (matrix site factor gene);
	push @classes, qw( TRANSFAC::Site TRANSFAC::Factor TRANSFAC::Gene TRANSFAC::Matrix );
	push @class_factories, ($site_holder, $factor_holder, $gene_holder, $matrix_holder);
    }

    ## Index classes to be parsed
    foreach my $type (@data_types) {
	$parse{$type} = 1;
    }

    #### output directory
    &RSAT::util::CheckOutDir($dir{output});
#    &CheckOutputDir();

    #### check existence of input and output directories
    foreach $d (keys %dir) {
	warn "; Checking existence of $d directory $dir{$d}\n" if ($main::verbose >= 1);
	unless (-d ($dir{$d}) ) {
	    die "Error: $d directory $dir{$d} does not exist\n";
	}
    }

    #### Check existence of input files
#    foreach my $file ("factor") {
    foreach my $type (@data_types) {
	$in_file{$type} = $dir{transfac}."/".$type.".dat";
	unless (-e $in_file{$type}) {
	    &FatalError(join ("\t", "Input file for data type", $type, "does not exist", $in_file{$type}));
	}
    }
     
    ### Report files
    $out_file{error} = "$dir{output}/transfac.errors.txt";
    open ERR, ">$out_file{error}" || &RSAT::error::FatalError("cannot write error file $out_file{error}");
    $out_file{stats} = "$dir{output}/transfac.stats.txt";
    $out_file{transfac} = "$dir{output}/transfac.obj" if ($export{obj});

    ### test conditions
    if ($test) {
	warn ";TEST\n" if ($main::verbose >= 1);
	### fast partial parsing for debugging
	foreach $key (keys %in_file) {
	    $in_file{$key} = " head -5000 $in_file{$key} |";
	}
    }

    if ($main::verbose >= 1) {
	&DefaultVerbose();
	warn sprintf "%-20s\t%d\n", "parse_transpro", $parse_transpro;
	warn sprintf "%-20s\t%s\n", "host", $host;
	warn sprintf "%-20s\t%s\n", "schema", $schema;
	warn sprintf "%-20s\t%s\n", "user", $user;
    }


    ### parse data from original files
    foreach my $type (@data_types) {
	my $holder = $type."_holder";
	&ParseTransfacFile($in_file{$type}, $$holder);
    } 

    ## Data type-specific treatment
    &TreatPromoters() if ($parse{transpro});
    &TreatSites() if ($parse{site});
    &TreatFactors() if ($parse{factor});
    &TreatGenes() if ($parse{gene});
    &TreatMatrices() if ($parse{matrix});

    ### print result
    &PrintStats($out_file{stats}, @classes);
    foreach my $class_factory (@class_factories) {
	warn (join "\t", "; Dumping class", $class_factory->get_object_type()), "\n" if ($main::verbose >= 1);
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


    ################################################################
    ### Report execution time
    my $exec_time = &RSAT::util::ReportExecutionTime($start_time); ## This has to be exectuted by all scripts
    warn $exec_time if ($main::verbose >= 1); ## only report exec time if verbosity is specified

    close ERR;
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
	parse_transfac.pl

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

	-outdir output 
		If not specified, a directory is automatically
		created in your home directory.
		${HOME}/parsed_data/transfac/[YYYYMMDD]
		where [YYYYMMDD] indicates the current date

	-v #	warn level
		Warn level 1 corresponds to a restricted verbose
		Warn level 2 reports all polypeptide instantiations
		Warn level 3 reports failing get_attribute()

	-obj    Export results in "obj" format (human readable)

	-clean	remove all files from the output directory before
		parsing
	-transpro
		parse TRANSPRO (promoter database) 
		This requires the commercial version of TRANSPRO(R).
		The input directrry must then contain a file called
		    transpro.tat

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
	    
# 	    ### deliver
# 	} elsif ($ARGV[$a] eq "-deliver") {
# 	    $main::deliver = 1;
	    

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
	    $main::schema = $ARGV[$a+1];
	    
	    ### database host
	} elsif ($ARGV[$a] eq "-host") {
	    $main::host = $ARGV[$a+1];
	    
	    ### database user
	} elsif ($ARGV[$a] eq "-user") {
	    $main::user = $ARGV[$a+1];
	    
	    ### password 
	} elsif ($ARGV[$a] eq "-password") {
	    $main::password = $ARGV[$a+1];
	    
	    ### transpro
	} elsif ($ARGV[$a] eq "-transpro") {
	    $main::parse_transpro = 1;

	}
	
    }
}




################################################################
# Parse TRANSFAC regulatory sites from sites.dat
#
sub ParseTransfacFile {
    my ($input_file, $class_holder, $source) = @_;
    
    &RSAT::message::TimeWarn(join ("\t","Parsing file", $input_file, "class holder", $class_holder->get_object_type(), "source", $source))
	if ($main::verbose >= 1);

    open DATA, $input_file || 
	die "Error: cannot open data file ", $input_file, "\n";

    #### Read an entire record at a time
    $/ = "\/\/\n";


    ## Parse the entry
    my $entries = 0;
    while ($text_entry = <DATA>) {
      $entries++;
      warn "; $input_file\t$entries entries\n" if ($main::verbose >=2);

      if ($entries==1) {
	#### check if there is a header
	if (($text_entry =~ /TRANSFAC \S+ TABLE/i) ||
	    ($text_entry =~ /TRANSPro/i))
	  {
	    if ($' =~ /Release (\S+)/) { 
	      $transfac_version = $1;
	    }
	    if ($' =~ /(\d+\-\d+\-\d+)/) { 
	      $release_date = $1;
	    }
	    if ($main::verbose >= 1) {
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


	## If entry is a matrix, save it in a separate file in TRANSFAC format
	if ($class_holder->get_object_type() eq "TRANSFAC::Matrix") {
	  chdir ($dir{output});
	  unless (-d "matrices_tf_format") {
	    system "mkdir -p matrices_tf_format";
	  }
	  if ($text_entry =~ /AC +(M\d+)\s*\n/) {
	    my $ac = $1;
	    my $matrix_file = "matrices_tf_format/".$ac.".tf";
	    my $out = &OpenOutputFile($matrix_file);
	    print $out $text_entry;
	    close  $out;
	    &RSAT::message::Info("Exported matrix in TRANSFAC format", $dir{output}."/".$matrix_file) if ($main::verbose >= 3);
	  } else {
	    &RSAT::error::FatalError("Matrix entry without AC field", "\n".$text_entry);
	  }
	}

	#### parse the entry
	my @lines = split "\n", $text_entry;
	foreach my $line (@lines) {
	    chomp $line;
	    next unless ($line); # skip empty lines
	    next if ($line =~ /^XX$/); # skip separator lines
	    next if ($line =~ /^\/\/$/); # skip record separator

	    if ($line =~ /^(\S{2})  /) {
		my $key = $1; ## Attribute name
		my $value = $'; #'; ## Attribute value

		################################################################
		#### Create new object 
		if ($key eq "AC") {
		    $current_obj = $class_holder->new_object(source=>$source, 
							     id=>$value, 
							     transfac_ac=>$value);
		    $transfac_index{$value} = $current_obj;
		    next;

		    ## Specific treatment for matrices
		} elsif ($class_holder->get_object_type() eq "TRANSFAC::Matrix") {
		    if ($key eq "P0") { ## Matrix header
			$value = &RSAT::util::trim($value);
			unless ($value =~ /^\s*A\s+C\s+G\s+T/i) {
			    ## Check if all matrices have the same header
			    &ErrorMessage("Non-conform header", $value, "\n");
			    my @header = split /\s+/, $value;
			    $current_obj->push_expanded_attribute("pssm_header", join "\t", $key, @header); 
			}
			next;
		    } elsif ($key =~ /^\d+$/) {
			$value = &RSAT::util::trim($value);
			my @values = split /\s+/, $value;
			warn join "\t", "PSSM", $key, @values, "\n" if ($main::verbose >= 10); 
			$current_obj->push_expanded_attribute("pssm", join "\t", $key, @values);
			next;
		    }
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
      my $xdb = &RSAT::util::trim(shift @fields);
      my $xids = &RSAT::util::trim(shift @fields);
      $xids =~ s/\.$//; 
      my @xids = split ";", $xids;
      foreach my $xid (@xids) {
	$current_obj->push_expanded_attribute("xrefs", $xdb, &RSAT::util::trim($xid));
      }
    }
  }
}

################################################################
# Specific treatment for promoters (TRANSPRO)
#
sub TreatPromoters {
    ## Open a separate file for storing the sequence in fasta format
    &RSAT::message::Info(join ("\t", "; Exporting promoter sequences to file", 
			       $out_file{promoter_sequences})) if ($main::verbose >= 1);
    
    $out_file{promoter_sequences} = $dir{output}."/promoter_sequences.fasta";
    open SEQ, ">$out_file{promoter_sequences}";
    $out_file{features} = $dir{output}."/promoter_features.ft";
    open FT, ">$out_file{features}";
    foreach my $promoter ($transpro_holder->get_objects()) {
	my $prom_ac = $promoter->get_attribute("id");
	my $prom_id = $promoter->get_attribute("transfac_id");
	my $f = 0; ## Count features per promoter
	my $synonyms = join (", ", $promoter->get_attribute('synonyms'));
	my @synonyms = split(", ", $synonyms);
	my $prom_label = join ("|", $prom_id, $prom_ac, @synonyms);

	################################################################
	## Export features in a .ft format (in addition to the field
	## transpro_features)
	foreach my $feature_line ($promoter->get_attribute("feature_table")) {
	    $f++;
	    my $feature_id;
	    if ($feature_line =~ /\:/) {
		my $feature_type = $`;
		my @feature_fields = split("\s*;\s*", "$'");
		unless ($feature_type =~ /\S+/) {
		    &ErrorMessage("Invalid feature format: missing feature type", "promoter ID", $prom_id, $prom_ac, $feature_line);
		    next;
		}
		
		## Create a new feature object
		my $feature = $feature_holder->new_object();
		$feature->set_attribute("score", 20); ## Aritrary score, but convenient since the feature will be visible above patser results (max score ~10)
		$feature->set_attribute("promoter", $prom_ac);
		$feature->set_attribute("seq_name", $prom_label);
		
		## Treat position
		my $position = &trim($feature_fields[$#feature_fields]);
		if ($position =~ /^(\d+)\.\.(\d+)\.$/) {
		    $feature->set_attribute("start", $1);
		    $feature->set_attribute("end", $2);
		} elsif ($position =~ /^(\d+)\.$/) {
		    $feature->set_attribute("start", $1);
		    $feature->set_attribute("end", $1);
		} else {
		    &ErrorMessage("Invalid feature position", "promoter ID", $prom_id, $prom_ac, $feature_line);
		}
		$feature->set_attribute("strand", "DR");
		
		## Treat separately TRANSFAC SITES and other features, because they have different formats !!!
		my $feature_name;
		my $feature_source;
		if ($feature_type =~ /transfac site/i) {
		    $feature_source = "TRANSFAC";
		    $feature_id = &trim($feature_fields[0]);
		    $feature_name =  $feature_id;
		} else {
		    $feature_source = $feature_fields[0];
		    $feature_id =  $feature_fields[1];
		    $feature_name = $feature_type;
		}
		$feature->set_attribute("id", $feature_id);
		$feature->set_attribute("ft_type", $feature_type);
		$feature->set_attribute("feature_name",$feature_name);
		$feature->set_attribute("description", join ("; ", $prom_id, $prom_ac, $feature_type, $feature_source, $feature_id, $synonyms));

		print FT $feature->to_text("ft", $null);
	    } else {
		&ErrorMessage("Invalid feature format: missing column", "promoter ID", $prom_id, $prom_ac, $feature_line);
	    }
	}

	################################################################
	## Export promoter sequences
	
	## Read the chromosomal location from the comments
	my @comments = $promoter->get_attribute("comments");
	foreach my $comment (@comments) {
	    if ($comment =~  /Sequence Fragment/) {
		if ($comment =~ /Sequence Fragment\s+\[(.*)\]/) {
		    $promoter->set_attribute("build", $1);
		}
		if ($comment =~ /\:\s+chr.(\S+)\s+(\d+)\.\.(\d+)/i) {
		    $promoter->set_attribute("chromosome", $1);
		    $promoter->set_attribute("left", $2);
		    $promoter->set_attribute("right", $3);
		    if ($comment =~ /FORWARD/i) {
			$promoter->set_attribute("strand", "D");
		    } elsif ($comment =~ /REVERSE/i) {
			$promoter->set_attribute("strand", "R");
		    } else {
			&ErrorMessage(join("\t", "Invalid strand specification for promoter", $prom_ac, $comment));
		    }
		}
	    }
	}

	## Concatenate the sequence fragments and export them
	my $seq_comment = $promoter->get_attribute("organism");
	$seq_comment .= "; ";
	$seq_comment .= $promoter->get_attribute("descr");
	my $sequence = join ("", $promoter->get_attribute($key_alias{"SQ"})); ## join for possible multi-line sequence
	
	if ($sequence =~ /\.$/) {
	    $sequence = $`;
#	} elsif ($sequence) {
#	    &ErrorMessage (join ("\t", "unproperly terminated sequence (no dot)", $promoter->get_attribute("ac"), $sequence), "\n");
	}
	if ($sequence) {
	    $promoter->set_array_attribute("seq"); ### empty the seq vector
#	    $promoter->push_attribute("sequence", $sequence);
	    &PrintNextSequence(SEQ,"fasta",60,$sequence,$prom_label, $seq_comment);
	}
    }
    close SEQ;
    close FT;

}

################################################################
# Specific treatment for sites
#
sub TreatSites {
    $site_holder->set_attribute_header("binding_factors_expanded", join("\t", "factor_ac", "factor_name", "quality"));

    foreach my $site ($site_holder->get_objects()) {
	#### Description
	my $description = &join_attribute($site, "descr", "", "description");

	#### organism classification
	my $taxonomy = &join_attribute($current_obj, "organism_classification", " ", "taxonomy");

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
		$site->force_attribute("gene_name", $gene_name);
		warn join ("\t", "; Site",  
			   $site->get_attribute("id"),
			   "Gene AC", $gene_ac,
			   "Gene name", $gene_name,
			  ), "\n" if($main::verbose >= 10);
	    } else {
		&ErrorMessage(join("\t", "; Site", 
				   $site->get_attribute('id'),
				   "Cannot identify gene", $gene_ac), "\n");
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
  foreach my $gene ($gene_holder->get_objects()) {
    #### organism classification
    my $taxonomy = &join_attribute($current_obj, "organism_classification", " ", "taxonomy");

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

  &parse_database_references($gene_holder);
}

################################################################
## post-parsing for factors (assign explicit field names)
sub TreatFactors {
  foreach my $factor ($factor_holder->get_objects()) {
    #### organism classification
    my $taxonomy = &join_attribute($current_obj, "organism_classification", " ", "taxonomy");

    ## Factor name(s)
    my @names = $factor->get_attribute($key_alias{FA});
    my $primary_name = $names[0];
    unless ($primary_name) {
      $primary_name = $factor->get_attribute("transfac_id");
    };
    $factor->set_attribute("name", $primary_name);
    $factor->push_attribute("names", $primary_name);
  }
}


################################################################
## post-parsing for matrices (assign explicit field names)
sub TreatMatrices {
    warn "; Treating matrices\n" if ($main::verbose >= 1);

    ## Description
    foreach my $matrix ($matrix_holder->get_objects()) {
	$matrix->set_attribute("description", $matrix->get_attribute("descr"));
    }

    ## Export the matrix in tab-delimited format and extract consensus
    chdir ($dir{output});
    unless (-d "matrices") {
	system "mkdir -p matrices";
    }
    foreach my $matrix ($matrix_holder->get_objects()) {
	warn "Exporting PSSM for matrix", $matrix->get_attribute("id"), "\n" if ($main::verbose >= 20);
	my $consensus  = "";
	my @pssm = ();
	my @columns = $matrix->get_attribute("pssm");
	my $max_pos = -1;
	foreach my $column_ref (@columns) {
	    my $fields = shift @{$column_ref};
	    my @fields = split "\t", $fields;
	    my $pos = shift @fields; 
	    for my $i (0..3) {
		$pssm[$pos][$i] = shift @fields;
		$max_pos = $pos if ($pos > $max_pos);
	    }
	    $consensus .= shift @fields;
#	    warn join ("\t", $matrix->get_attribute("id"), $column_ref, $pos, $consensus), "\n" if ($main::verbose >= 10)
	}

	## Export the matrix in a separate file
	my $file_name = $matrix->get_attribute("id");
	my $pssm = "";
	$file_name .= ".matrix";
	$file_name =~ s/\$/_/g;
	open MATRIX, ">matrices/".$file_name || die "Cannot write matrix file", $file_name, "\n";
	@residues = qw( A C G T ); 
	for my $r (0..$#residues) {
	    $pssm .= $residues[$r]." |";
	    for my $pos (1..$max_pos) {
		$pssm .= $pssm_sep.$pssm[$pos][$r];
	    }
	    $pssm .= "\n";
	}
	print MATRIX  $pssm;
	close MATRIX;

	## Set the consensus attribute
	$matrix->set_attribute("consensus", $consensus); 
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
