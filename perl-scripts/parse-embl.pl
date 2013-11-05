#!/usr/bin/perl
############################################################
#
# $Id: parse-embl.pl,v 1.30 2013/06/18 20:18:28 jvanheld Exp $
#
# Time-stamp: <2003-10-21 01:17:49 jvanheld>
#
############################################################
#use strict;;
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}

require "RSA.lib";
push @INC, "$ENV{RSAT}/perl-scripts/parsers/";
require "lib/load_classes.pl";
require "lib/parsing_util.pl";


################################################################
#### Class for EMBL feature
package EMBL::Feature;
{
  @ISA = qw ( classes::DatabaseObject );
  ### class attributes
  $_count = 0;
  $_prefix = "ft_";
  @_objects = ();
  %_name_index = ();
  %_id_index = ();
  %_attribute_count = ();
  %_attribute_cardinality = (id=>"SCALAR",
			     names=>"ARRAY",
			     organism=>"SCALAR",
			     type=>"SCALAR",
			     description=>"SCALAR",
			     position=>"SCALAR",
			     contig=>"SCALAR",
			     strand=>"SCALAR",
			     start_pos=>"SCALAR",
			     end_pos=>"SCALAR",
			     GeneID=>"SCALAR",
			     xrefs=>"EXPANDED"
			     );
}

################################################################
#### Class for CDS
package EMBL::CDS;
{
  @ISA = qw ( classes::DatabaseObject );
  ### class attributes
  $_count = 0;
  $_prefix = "cds_";
  @_objects = ();
  %_name_index = ();
  %_id_index = ();
  %_attribute_count = ();
  %_attribute_cardinality = (id=>"SCALAR",
			     names=>"ARRAY",
			     organism=>"SCALAR",
			     type=>"SCALAR",
			     description=>"SCALAR",
			     position=>"SCALAR",
			     contig=>"SCALAR",
			     strand=>"SCALAR",
			     start_pos=>"SCALAR",
			     end_pos=>"SCALAR",
			     GeneID=>"SCALAR",
			     xrefs=>"EXPANDED"
			    );
}

################################################################
#### Class for tRNA
package EMBL::tRNA;
{
  @ISA = qw ( classes::DatabaseObject );
  ### class attributes
  $_count = 0;
  $_prefix = "trna_";
  @_objects = ();
  %_name_index = ();
  %_id_index = ();
  %_attribute_count = ();
  %_attribute_cardinality = (id=>"SCALAR",
			     names=>"ARRAY",
			     organism=>"SCALAR",
			     type=>"SCALAR",
			     description=>"SCALAR",
			     position=>"SCALAR",
			     contig=>"SCALAR",
			     strand=>"SCALAR",
			     start_pos=>"SCALAR",
			     end_pos=>"SCALAR",
			     GeneID=>"SCALAR",
			     xrefs=>"EXPANDED"
			    );
}

################################################################
#### Class for rRNA
package EMBL::rRNA;
{
  @ISA = qw ( classes::DatabaseObject );
  ### class attributes
  $_count = 0;
  $_prefix = "rrna_";
  @_objects = ();
  %_name_index = ();
  %_id_index = ();
  %_attribute_count = ();
  %_attribute_cardinality = (id=>"SCALAR",
			     names=>"ARRAY",
			     organism=>"SCALAR",
			     type=>"SCALAR",
			     description=>"SCALAR",
			     position=>"SCALAR",
			     contig=>"SCALAR",
			     strand=>"SCALAR",
			     start_pos=>"SCALAR",
			     end_pos=>"SCALAR",
			     GeneID=>"SCALAR",
			     xrefs=>"EXPANDED"
			    );
}

################################################################
#### Class for mRNA
package EMBL::mRNA;
{
  @ISA = qw ( classes::DatabaseObject );
  ### class attributes
  $_count = 0;
  $_prefix = "mrna_";
  @_objects = ();
  %_name_index = ();
  %_id_index = ();
  %_attribute_count = ();
  %_attribute_cardinality = (id=>"SCALAR",
			     names=>"ARRAY",
			     organism=>"SCALAR",
			     type=>"SCALAR",
			     description=>"SCALAR",
			     position=>"SCALAR",
			     contig=>"SCALAR",
			     strand=>"SCALAR",
			     start_pos=>"SCALAR",
			     end_pos=>"SCALAR",
			     GeneID=>"SCALAR",
			     xrefs=>"EXPANDED"
			    );
}


################################################################
#### Class for EMBL organism
package EMBL::Organism;
{
  @ISA = qw ( classes::DatabaseObject );
  ### class attributes
  $_count = 0;
  $_prefix = "ft_";
  @_objects = ();
  %_name_index = ();
  %_id_index = ();
  %_attribute_count = ();
  %_attribute_cardinality = (id=>"SCALAR",
			     names=>"ARRAY",
			     taxonomy=>"SCALAR",
			     );
}

################################################################
#### Class for EMBL contig
package EMBL::Contig;
{
  @ISA = qw ( classes::DatabaseObject );
  ### class attributes
  $_count = 0;
  $_prefix = "ft_";
  @_objects = ();
  %_name_index = ();
  %_id_index = ();
  %_attribute_count = ();
  %_attribute_cardinality = (id=>"SCALAR",
			     accession=>"SCALAR",
			     version=>"SCALAR",
			     type=>"SCALAR",
			     length=>"SCALAR",
			     description=>"SCALAR",
			     );
}




################################################################
#### main package
package main;
{

    #### initialise parameters ####
    local $start_time = &RSAT::util::StartScript();

    ### initial directory
    $dir{main} = `pwd`; #### remember working directory
    chomp($dir{main});

    local $data_source = "embl";
    local %infile = ();
    local %outfile = ();
    local $verbose = 0;
    local $in = STDIN;
    local $out = STDOUT;
    local $test_lines = 1000;
    local $feature_count = 0;
    %feature_out_fields = (); #### index of feature out fields, for checking extra fields
    %feature_extra_fields = (); #### index of feature out fields, for checking extra fields


    #### SQL options
    $host= $default{'host'};
    $schema="embl";
    $user="embl";
    $password="embl";

    ## Extension for embl-formatted files
    $extension = ".embl";


#    my $genes = classes::ClassFactory->new_class(object_type=>"EMBL::Gene",
#						 prefix=>"gene_");

    ################################################################
    #### Open class factories

    #### Organism
    local $organisms = classes::ClassFactory->new_class(object_type=>"EMBL::Organism",prefix=>"org_");
    $organisms->set_out_fields(qw (id
				   taxonomy
				   name
				   names
				   ));

    #### Contig
    local $contigs = classes::ClassFactory->new_class(object_type=>"EMBL::Contig",prefix=>"ctg_");
    $contigs->set_out_fields(qw (id
				 accession
				 version
				 form
				 type
				 file
				 length
				 description
				 ));

    ## Specific feature types
    local $CDSs = classes::ClassFactory->new_class(object_type=>"EMBL::CDS",prefix=>"cds_", source=>$data_source);
    local $mRNAs = classes::ClassFactory->new_class(object_type=>"EMBL::mRNA",prefix=>"mrna_", source=>$data_source);
    local $tRNAs = classes::ClassFactory->new_class(object_type=>"EMBL::tRNA",prefix=>"trna_", source=>$data_source);
    local $rRNAs = classes::ClassFactory->new_class(object_type=>"EMBL::rRNA",prefix=>"rrna_", source=>$data_source);
    
    #### Features
    local $features = classes::ClassFactory->new_class(object_type=>"EMBL::Feature",prefix=>"ft_", source=>$data_source);
    @feature_out_fields = qw( id
			      type
			      name
			      contig
			      start_pos
			      end_pos
			      strand
			      description
			      GeneID
			      organism
			      position
			      names
			      db_xref
			      introns
			      exons

			      EC_number

			      codon_start
			      gene
			      protein_id
			      transl_table
			      function
			      isolate
			      mol_type
			      note
			      product
			      translation
			      strain
			      locus_tag
			      sub_strain
			      bound_moiety

			      db_xref_exp
			      ); 
    $features->set_out_fields(@feature_out_fields);


    #### index feature out fields for later checking
    foreach my $field (@feature_out_fields) {
	$feature_out_fields{$field}++;
    }
#die join "'\n'", keys %feature_out_fields, "\n";

    #### classes
    @feature_holders = ($features, $CDSs, $mRNAs, $tRNAs, $rRNAs);

    @classes = qw ( EMBL::Organism
		    EMBL::Contig
		    EMBL::Feature
		    EMBL::CDS		
		    EMBL::mRNA
		    EMBL::tRNA
		    EMBL::rRNA
  );

    #### working directory
    $wd = `pwd`;
    chomp($wd);

    &ReadArguments();

    ################################################################
    #### check argument values ####

    #### input directory
    unless (defined($dir{input})) {
	&RSAT::error::FatalError("You must specify the input directory.\n");
    }
    unless (-d $dir{input}) {
	&RSAT::error::FatalError("Input directory '$dir{input}' does not exist.\n");
    }

    #### First guess for organism
    unless (defined($org)) {
	$org = `basename $dir{input}`;
	chomp($org);
	&RSAT::message::Info("Auto selection of organism name", $org) if ($main::verbose >= 1);
    }

    #### find embl files in the input directory
    chdir ($dir{input});
    @embl_files =();
    push @embl_files, glob("*".$extension);
    push @embl_files, glob("*".$extension.".Z");
    push @embl_files, glob("*".$extension.".gz");
    if ($#embl_files < 0) {
	&RSAT::error::FatalError("There is no file with extension ${extension} in the input directory $dir{input}\n");
    } else {
	&RSAT::message::Info("EMBL files\n;\t", join("\n;\t", @embl_files)) if ($main::verbose >= 1);
    }
    #### come back to the starting directory
    chdir($wd);

    #### output directory
    unless (defined($dir{output})) {
	#### default output directory
	$export_subdir = "embl";
	$dir{output} = "$parsed_data/$export_subdir/$delivery_date/$org";
	$dir{output} = "parsed_embl/".$org."/".$delivery_date;
#	$dir{output} = "$ENV{RSAT}/public_html/data/embl_genomes/$org/genome";
	&RSAT::message::Warning("Auto selection of output dir", $dir{output}) if ($main::verbose >= 1);
    }
    &CheckOutputDir($dir{output});
    chdir $dir{output};
    $outfile{error} = "$dir{output}/genbank.errors.txt";
    $outfile{stats} = "$dir{output}/genbank.stats.txt";
    open ERR, ">$outfile{error}";

    #### come back to the starting directory
    chdir($wd);

    #### verbose ####
    &Verbose() if ($main::verbose >= 1);

    #### parse the embl files
    $contig_handle = &OpenOutputFile("$dir{output}/contigs.txt"); # file with contig IDs
    foreach my $file (@embl_files) {
      &ParseEMBLFile("$dir{input}/$file", 
		     source=>$file);
    }
    close $contig_handle;


    ################################################################
    #### Postprocessing : we interpret the features to extract more information
    #### check some feature attributes (name, description, ...)
    &RSAT::message::TimeWarn("Processing features") if ($main::verbose >= 1);
    my $organism_name = $org; #### first guess
    my $contig = "";
    foreach my $class_holder (@feature_holders) {
      &RSAT::message::TimeWarn("Processing features of class", $class_holder->get_object_type()) if ($main::verbose >= 1);
      foreach $feature ($class_holder->get_objects()) {
	################################################################
	## The organism name is annotated in the feature of type
	## "source"
	if ($feature->get_attribute("type") eq "source") {
	  $organism_name = $feature->get_attribute("organism");
	  if (($organism_name) && ($organism_name ne $null)) {
	    &CreateOrganism($organisms, $organism_name);
	  }

	  $contig = $feature->get_attribute("contig");
	  my @xrefs = $feature->get_attribute("db_xref");
	  foreach my $xref (@xrefs) {
	    if ($xref =~ /taxon:(\S+)/) {
	      $taxon = $1;
	    }
	  }
	  if (!($contig) || ($contig eq $null)) {
	    $contig = $feature->get_attribute("contig");
	    $feature->force_attribute("contig", $contig);
	  }
	} else {
	  $feature->set_attribute("organism", $organism_name);
	  #	   $feature->force_attribute("contig", $contig);
	}

        ### Add gene attribute as name
	my @gene = $feature->get_attribute("gene");
	push @gene, $feature->get_attribute("Gene"); ## In some files, the attribute "gene" starts with an uppercase
#	&RSAT::message::Debug( "Gene name(s)", $feature->get_attribute("type"), "id=".$feature->get_attribute('id'), "gene", join(";", @gene)) if ($main::verbose >= 5);
	if (scalar(@gene) >= 1) {
	  foreach my $gene (@gene) {
#	    &RSAT::message::Debug( "Adding gene name", $gene, $feature->get_attribute("type"), $feature->get_attribute('id')) if ($main::verbose >= 5);
	    if ($feature->get_attribute("type") eq "CDS") {
	      $feature->push_attribute("names",$gene);
	    } else {
	      $feature->push_attribute("names",$gene." ".$feature->get_attribute("type"));
	    }
	  }

	  ## Use gene name as main ID for the CDS 
	  ## (this is just a first guess, it will be replaced by protein_id or locus_tag if these are defined)
	  if ($gene[0] =~ /(\S+)/) {
	    $feature->force_attribute('id', $1);
	  }
	}

	#	## Use attribute "accession" as main ID if there is one
	#	my $accession = $feature->get_attribute("accession");
	#	if ($accession =~ /(\S+)/) {
	#	    $feature->force_attribute('id', $1);
	#	}

	### add protein_id as valid name
	@protein_ids = $feature->get_attribute("protein_id");
	if ($#protein_ids >= 0) {
	  foreach my $id (@protein_ids) {
	    $feature->push_attribute("names",$id);
	  }

	  ## Use protein_id as main ID for the CDS
	  ## (this is just a second guess, it will be replaced by locus_tag if it is defined)
	  if ($protein_ids[0] =~ /(\S+)/) {
	    $feature->force_attribute('id', $1);
	  }
	}

	### add locus_tag as valid name
	@locus_tags = $feature->get_attribute("locus_tag");
	if ($#locus_tags >= 0) {
	  foreach my $id (@locus_tags) {
	    $feature->push_attribute("names",$id);
	  }

	  ## Use locus_tag as main ID if there is one
	  if ($locus_tags[0] =~ /(\S+)/) {
	    $feature->force_attribute('id', $1);
	  }
	}

	#### add SWISS-PROT/Uniprot/TrEMBL ID as valid name
	my @xrefs = $feature->get_attribute("db_xref");
	my $gi = "";
	foreach my $xref (@xrefs) {
	  if (($xref =~ /SWISS-PROT:/i) ||
	      ($xref =~ /Uniprot\/TrEMBL:/i) ||
	      ($xref =~ /SPTREMBL:/i)) {
	    my $name = "$'";
	    $feature->push_attribute("names",$name);		
	  } 
	}
	
	### define a single name  (take the first value in the name list)
	if ($name = $feature->get_name()) {
	  $feature->set_attribute("name",$name);
	} else {
	  $feature->set_attribute("name",$feature->get_id());
	}
	
	#### check for features without description
	if (($feature->get_attribute("description") eq $null) 
	    || ($feature->get_attribute("description") eq "")) {
	  my $description =  join (";", $feature->get_attribute("product"));
	  unless ($description) {
	    $description =  join (";", $feature->get_attribute("note"));
	  }
	  $feature->set_attribute("description",$description);
        }
	
	
	#### parse cross-references
	my @xrefs = $feature->get_attribute("db_xref");
	$class_holder->set_attribute_header("db_xref_exp", "xref_db", "xref_id");
	my $gi = "";
	foreach my $xref (@xrefs) {
	  if ($xref =~ /:/) {
	    my $xref_db = $`;
	    my $xref_id = "$'";
	    $feature->push_expanded_attribute("db_xref_exp", $xref_db, $xref_id);
	  } else {
	    &ErrorMessage(join ("\t", 
				"Invalid cross reference", 
				$feature->get_attribute("id"), 
				$xref), "\n");
	  }
	}
      }

      ## Specify the GeneID (required for retrieve-seq)
      &RSAT::message::TimeWarn("Specifying GeneIDs for class", $class_holder->get_object_type()) if ($main::verbose >= 2);
      foreach my $feature ($class_holder->get_objects()) {
	$feature->set_attribute("GeneID", $feature->get_attribute("id"));
	#	&RSAT::message::Debug("Using ID as GeneID for feature", $feature->get_attribute("id"),$feature->get_attribute("GeneID")) if ($main::verbose >= 10);
      }
    }
    ################################################################
    ### Save result in tab files
    chdir $dir{main};
    #    chdir $dir{output};
    foreach $class_factory ($organisms, $contigs, @feature_holders) {
      $class_factory->dump_tables();
      $class_factory->generate_sql(dir=>"$dir{output}/sql_scripts",
				   schema=>$schema,
				   host=>$host,
				   user=>$user,
				   password=>$password
				  );
    }
    
    &ExportProteinSequences($features,$org);
    &ExportMakefile(@classes);

    chdir($dir{output});
    &PrintStats($outfile{stats}, @classes);

    #### document the ffields which were parsed but not exported
    open STATS, ">>$outfile{stats}";
    print STATS "; \n; Fields parsed but not exported\n";
    foreach my $key (sort keys %feature_extra_fields) {
	my $value = $feature_extra_fields{$key};
	printf STATS ";\t%-20s\t%d\n", $key, $value;
    }
    close STATS;

    ###### close output file ######
    my $exec_time = &RSAT::util::ReportExecutionTime($start_time);
    print $main::out $exec_time if ($main::verbose >= 1);
    close $out if ($outfile{output});
    close ERR;

    &RSAT::message::TimeWarn("Results exported in directory", $dir{output}) if ($main::verbose >= 1);

    exit(0);
}

########################## subroutine definition ############################

################################################################
#### display full help message #####
sub PrintHelp {
  open HELP, "| more";
  print HELP <<End_of_help;
NAME
	parse-embl

        2001 by Jacques van Helden (jvanheld\@bigre.ulb.ac.be)
	
USAGE
        parse-embl [-dir input_dir][-i inputfile] [-o outputfile] [-v]

DESCRIPTION
	Parse one or sveral EMBL files for extracting genome
	information.

CATEGORY
	parser

OPTIONS
	-h	(must be first argument) display full help message
	-help	(must be first argument) display options

	-v	verbose

	-test #	quick test (for debugging): only parse the # first
		lines of each Genabnk file (default $test_lines).

	-i	input directory
		input directory. This directory must contain one or
		several embl files (extension .contig). 

	-ext	extension for EMBL flat files (default: $extension)

	-o	output directory
		The parsing result will be saved in this directory. If
		the directory does not exist, it will be created.

	-org	organism name

	-noseq  do not export sequences in .raw files

	-source	data source (default: $data_source)

   Options for the automaticaly generated SQL scripts
	-schema database schema (default: $schema)
	-host	database host (default: $host)
	-user	database user (default: $user)
	-password
		database password (default: $password)

INPUT FILE

        This parser takes as input all the .embl and .embl.Z files
        found in the input directory. Each file is parsed and the
        results are exported in a single directory (the output
        directory). 

	One input directory is supposed corresponds to one genome. For
        eucaryotes (and some procaryotes), a directory can contain
        several .embl files, one per chromosome.Features found in
        different files are merged in the same output file
        (features.tab).

OBJECT TYPES

        The parser is based on a syntactic decomposition of the .embl
        files. Internally, it creates the following objects

	    EMBL::Organism
	    EMBL::Contig
	    EMBL::Feature

	Each feature has a type (CDS, mRNA, tRNA, rRNA, promoter,
	...), as found in the original EMBL file. The feature type si
	indicated by a line starting with FT, followed by 3 white
	spaces, the type, some additional spaces, and the position.

	There is thus one generic object called 'feature', which has a
	type attribute. An alternative possibility would have been to
	create one distinct object type for each feature
	type. However, since the types of features encountered vary
	from file to file, we preferred to leave the feature
	completely generic, so we avoid to skip features which would
	not have been encountered in the files used for the
	development of the parser.

	Each feature has several attributes. The attributes names are
	assigned dynamically on the basis of a syntactic
	analysis. This means that the parser might create attributes
	which had bever been assigned before. This makes a problem,
	since one organism would export some tables which were not
	found in another organism. Since tables are created only once
	(before loading the first organism), the new data could not be
	loaded. To avoid this, the attributes to be exported are
	specified in the parser. Unrecognized attributes are parsed
	and appear in the file embl_tats.txt, but the tables and SQL
	scripts for loading these attributes are not exported.

OUTPUT

	The parsed data is exported in several tab-delimited files:
	one main table with all the single-value attributes, and one
	separate table for each multi-value attribute.

	However, in this parser, I do not strictly follow this
	normalization rule, for the sae of space economy. Indeed, some
	single-value attributes are specific for some feature types
	(e.g. the attribute 'product' is found for the feature type
	'CDS' but not for feature type 'promoter'; the attribute
	'bound_moiety' is assigned to feature type 'protein_bind',
	which is annotated in some genomes like ecoli_K12, but absent
	from most other genome annotations). Thus, having these fields
	in the main table feature.tab would result in columns almost
	full of <NULL>. To avoid this, I export these fields in side
	tables, with a foreign key to the feature ID.


	The program also exports a file embl_stats.txt, with
	statistics about the number of objects parsed for each calss,
	and the number of attribute allocations. The stats file also
	contains a list of attributes which were found in the original
	EMBL file, but not exported because they were not part of the
	predefined list of attributes to export
	(\@feature_out_fields). This allows to complete this list when
	some features are annotated in some organisms, which were not
	identified in the organisms used for the development of this
	program.

End_of_help
  close HELP;
  exit;
}

################################################################
#### display short help message #####
sub PrintOptions {
  open HELP, "| more";
  print HELP <<End_short_help;
parse-embl options
----------------
-h	(must be first argument) display full help message
-help	(must be first argument) display options
-test #	quick test (for debugging)
-i	input dir
-o	output dir
-source	data source (default: $data_source)
-v	verbose
-ext	extension for EMBL flat files (default: $extension)
-org	organism name
-prefid preferred ID for a given feature type
-schema database schema (default: $schema)
-host	database host (default: $host)
-user	database user (default: $user)
-pass	database password (default: $password)
-noseq  do not export sequences in .raw files
End_short_help
  close HELP;
  exit();
}

################################################################
#### read arguments ####
sub ReadArguments {
    foreach my $a (0..$#ARGV) {
	### verbose ###
	if ($ARGV[$a] eq "-v") {
	    if (&IsNatural($ARGV[$a+1])) {
		$main::verbose = $ARGV[$a+1];
	    } else {
		$main::verbose = 1;
	    }

	    ### detailed help
	} elsif ($ARGV[$a] eq "-h") {
	    &PrintHelp();

	    ### list of options
	} elsif ($ARGV[$a] eq "-help") {
	    &PrintOptions();

	    ### quick test
	} elsif ($ARGV[$a] eq "-test") {
	    $test = 1;
	    if (&IsNatural($ARGV[$a+1])) {
		$test_lines = $ARGV[$a+1];
	    }

	    ### input file ###
	} elsif ($ARGV[$a] eq "-i") {
	    $dir{input} = $ARGV[$a+1];

	    ### extension
	} elsif ($ARGV[$a] eq "-ext") {
	    $extension = $ARGV[$a+1];

	    ### output file ###
	} elsif ($ARGV[$a] eq "-o") {
	    $dir{output} = $ARGV[$a+1];

	    ### organism ###
	} elsif ($ARGV[$a] eq "-org") {
	    $org = $ARGV[$a+1];
	    $org =~ s/\s+/_/g;

	    ### data source
	} elsif ($ARGV[$a] eq "-source") {
	    $data_source = $ARGV[$a+1];

	    ### do not export sequences
	} elsif ($ARGV[$a] eq "-noseq") {
	    $noseq = 1;

	    ################################################################
	    #### SQL database parameters

	    ### schema
	} elsif ($ARGV[$a] eq "-schema") {
	    $schema = $ARGV[$a+1];

	    ### host
	} elsif ($ARGV[$a] eq "-host") {
	    $host = $ARGV[$a+1];

	    ### user
	} elsif ($ARGV[$a] eq "-user") {
	    $user = $ARGV[$a+1];

	    ### password 
	} elsif ($ARGV[$a] =~ /^-pass/) {
	    $password = $ARGV[$a+1];

	}
    }
}


################################################################
#### print verbose message
sub Verbose {
    print $out "; parse-embl ";
    &PrintArguments($out);
    if (%dir) {
	print $out "; Directories\n";
	while (($key,$value) = each %dir) {
	    print $out ";\t$key\t$value\n";
	}
    }
    if (%main::infile) {
	print $out "; Input files\n";
	while (($key,$value) = each %infile) {
	    print $out ";\t$key\t$value\n";
	}
    }
    if (%outfile) {
	print $out "; Output files\n";
	while (($key,$value) = each %outfile) {
	    print $out ";\t$key\t$value\n";
	}
    }
    printf $out "; %-29s\t%s\n", "organism", $org;
    printf $out "; %-29s\t%s\n", "embl files", join (" ", @embl_files);
}



################################################################
#### parse one EMBL file
sub ParseEMBLFile {
    my ($input_file, %args) = @_;
    &RSAT::message::TimeWarn (join ("\t", "Parsing file", $input_file)) if ($main::verbose >= 1);
    
    my ($file,$dir) = &OpenInputFile($input_file);
#    open EMBL, $input_file 
#	|| die "Error: cannot open input file $input_file.\n";

    #### initialize parsing variables
    my $in_FT = 0; #### flag to indicate whether we are reading a feature
#    my $in_feature = 0;
#    my $in_cds = 0;

    my $current_feature = null;
    my $sequence = "";
    my $organism_name = "";
    my $organism;

    ################################################################
    #### parse the file
    my $l = 0; #### line counter
    while (my $line = <$file>) {
	$l++;
	if (($test) && ($l > $test_lines)) {
	    warn "Test: stopping at line $l\n";
	    last;
	}
	if (($main::verbose >= 1) && ($l%10000 == 0)) {
	  &RSAT::message::TimeWarn("Parsed", $l, "lines");
	}

#	warn $line if ($main::verbose >= 10);
	chomp $line;
	next unless ($line =~ /\S/);

	#### read the full sequence
	if  ($line =~ /^SQ/) {
	    $in_FT = 0;
	    $in_sequence = 1;
	    while (($in_sequence) && 
		   (my $line = <$file>)) {
		if ($line =~ /\d+\s*$/) {
		    $sequence .= $`;
		} elsif ($line =~ /^\/\/$/) {
		  ## The sequence is finished
		  $in_sequence = 0;
		  ## Store the sequence if required
		  unless ($noseq) {
		    my $seq_file = $contig->get_attribute("id");
		    $seq_file .= ".raw";
		    open RAW, ">$dir{output}/$seq_file";
		    &PrintNextSequence(RAW, "raw", 0, $sequence, $contig);
		    close RAW;
		    print $contig_handle $seq_file, "\t", $contig->get_attribute("id"), "\n";
		  }
		  undef($sequence);
		}
	    }
        }

	#### contig ID line
	if  ($line =~ /^ID\s+(\S+)\s*/) {
	    my $contig_id = $1;
	    $contig_id =~ s/\;$//;
	    if ($contig_id =~ /unknown/i) {
	      ## This is a bit tricky, but in Leishmania genome, most contigs havve as ID "unknown id"
	      $contig = $contigs->new_object();
	      $contig_id = $contig->get_attribute("id");
	    } else {
	      $contig = $contigs->new_object(id=>$contig_id);
	    }
	    &RSAT::message::Info("line", $l, "New contig", $contig->get_attribute("id")) if ($main::verbose >= 2);
	    my $contig_description = "$'"; 
	    $contig->set_attribute("description",$contig_description);
	    $contig->set_attribute("file",$input_file);

	    #### contig length
	    if ($contig_description =~ /(\d+)\s+BP\./i) {
		$contig->set_attribute("length", $1);
	    }

	    #### contig form (circular or linear)
	    if ($contig_description =~ /circular/i) {
		$contig->set_attribute("form", "circular");
	    } else {
		$contig->set_attribute("form", "linear");
	    }

	} elsif ($line =~ /^AC\s+/) {
	    #### accession number for the currrent contig
	    my $AC = $'; ## '
	    $AC =~ s/;$//;
	    $contig->set_attribute("accession", $AC);

	} elsif ($line =~ /^SV\s+/) {
	    #### sequence version for the currrent contig
	    my $version = $';
	    $contig->set_attribute("version", $version);

	} elsif ($line =~ /^OS\s+/) {
	    #### organism name
	    $organism_name = $';
	    &CreateOrganism($organisms, $organism_name);
	} elsif ($line =~ /^OC\s+/) {
	    my $taxonomy = "$'";

	    #### read the rest of the taxxonomy (it can be larger than
	    #### one line)
	    while ($line = <$file>) {
		chomp $line;
		$l++;
		if ($line =~ /^OC\s+/) {
		    $taxonomy .= "$'";
		} else {
		    last;
		}
	    }

	    #### organism taxonomy
	    if ($organism) {
		$taxonomy = &trim($taxonomy);
		$taxonomy =~ s/\.$//;
		$organism->set_attribute("taxonomy", $taxonomy);
	    }
	}


	#### start reading features
	unless ($in_FT) {
	    if ($line =~ /^FT\s+/) {
		$in_FT = 1 ;
		warn "; Reading features\n" if ($main::verbose >= 10);
	    }
	}


	#### new feature
	if ($line =~ /^FT   (\S+)\s+(.*)/) {

	    #### feature type
	    $feature_type = $1;

	    #### read the feature position
	    $position = &trim($2);
	    if ($position =~ /join\(/){
		### check that the position is complete
		my $start_line = $l;
		while ($line = <$file>) {
		    $l++;
		    chomp($line);
		    if (($line =~ /^FT   (\S+)\s+(.*)/) ||
			($line =~ /^FT    \s+\/(\S+)\=/)) {
			last;
		    } elsif ($line =~ /^FT    \s+([^\/]\S+)/) {
			$position .= $1;
		    } else {
			die "Error: feature position started at line $start_line not properly terminated at line $l\n";
		    }
		}
		
#  		unless ($' =~ /\)/) {
#  		    do {
#  			die "Error: position starting at line $l is not terminated properly.\n"
#  			    unless $position_suite = <$file>;
#  			$position_suite =~ s/^FT//;
#  			$position .= &trim($position_suite);
#  		    } until ($position =~ /\)/);
#  		}
	    }

	    ################################################################
	    ## Create a new object for the feature
	    $class_holder = $features;
	    if ($feature_type eq "CDS") {
	      $class_holder = $CDSs;
	    } elsif ($feature_type eq "mRNA") {
	      $class_holder = $mRNAs;
	    } elsif ($feature_type eq "tRNA") {
	      $class_holder = $tRNAs;
	    } elsif ($feature_type eq "rRNA") {
	      $class_holder = $rRNAs;
	    }
	      
	    $current_feature = $class_holder->new_object(%args);
	    $feature_count++;
	    &RSAT::message::Debug( "file", $input_file, "line", $l, "feature", $feature_count, $feature_type) 
	      if ($main::verbose >= 3);
	    $current_feature->set_attribute("type",$feature_type);
	    $current_feature->set_attribute("contig",$contig->get_attribute("id"));
#	    $current_feature->set_attribute("organism",$organism_name);
	    $current_feature->set_attribute("position",$position);
	    &ParsePositionEMBL($current_feature);

	    ### contig
#	    my $contig = $current_feature->get_attribute("source");
#	    $contig =~ s/\.embl.*//;
#	    $current_feature->set_attribute("contig",$contig);
	}



	#### new feature attribute
	if ($line =~ /^FT    \s+\/(\S+)\=/) {
	    my $key = $1;
	    #### check that this attribute is part of the exported ones
	    unless ($feature_out_fields{$key}) {
		$feature_extra_fields{$key}++;
	    }

	    my $value = "$'";
	    my $start_l = $l;

	    #### If the attribute starts with a quote, make sure to
	    #### have the closing quote, and trim the quotes
	    if  ($value =~ /^\"/) {
		$value = "$'"; #### trim the leading quote

		if ($value =~ /\"\s*$/) {
		    #### value terminates on the first line
		    $value = $`; #### retain what precedes the closing quote
		} else {
		    #### collect following lines until value is terminated
		    while ($line = <$file>) {
			chomp $line;
			$l++;
			if  ($line =~ /^FT\s+/) {
			    #### value is included in the first line
			    $value .= " $'";
			    if  ($value =~ /\"\s*$/) {
				#### value is terminated
				$value = $`;
				last;
			    }
			} else  {
			    die "Error: feature attribute started at line $start_l is not terminated at line $l\n";
			}
		    }
		}
	    }

	    #### remove spaces from sequences
	    if ($key eq "translation") {
		$value =~ s/\s//g;
	    }


	    #### add the new attribtue to the feature
	    if (($key eq "id") || ($key eq "type")) {
		$current_feature->force_attribute($key, $value);
	    } else {
		$current_feature->new_attribute_value($key, $value);
	    }
	    warn join ("\t", "; attribute", $feature_count, $start_l, $l, $key, $value), "\n" if ($main::verbose >= 4);

	}
	
    }
    close $file;
    return $sequence;
}


################################################################
#### parse the position of the current feature
sub ParsePositionEMBL {
    my ($current_feature) = @_;
    my $position = $current_feature->get_attribute("position");
    my $complete_position = $position;
    my $start_pos = $null;
    my $end_pos = $null;
    my $strand = $null;

    #### strand
    if ($position =~ /^complement\((.*)\)$/) {
	$position = $1;
	$strand = "R";
    } else {
	$strand = "D";
    }

    #### introns/exons
    if ($position =~ /^join\((.*)\)$/) {
	$position = $1;
    }
    my @exons = split ",", $position;
    

    #### The start of the first exon is the feature start, the end of
    #### the last exon is the feature end
    my @exon_starts = ();
    my @exon_ends = ();
    foreach my $exon (@exons) {

	#### strand !!! in some cases, the complement is WITHIN the join, and repeated for each exon 
	##  FT   CDS             join(complement(173191..173432),complement(172942..173107)
	##  FT                   )
	##  FT                   /db_xref="SPTREMBL:O11851"
	##  FT                   /note="ORF YCR028"
	##  FT                   /gene="RIM1"
	if ($exon =~ /^complement\((.*)\)$/) {
	    $exon = $1;
	    $strand = "R";
	}
	
	if ($exon =~ /^(<*\d+)..(>*\d+)$/) {
	    push @exon_starts, $1;
	    push @exon_ends, $2;
	} else {
	    &ErrorMessage("invalid exon specification\t", $exon, "\t", $complete_position, "\n");
	}
    }

    #### if the feature contains several exons infer introns
    if ($#exon_starts >= 1) {
	for my $e (0..$#exon_starts-1) {
	    my $intron = $exon_ends[$e]."..".$exon_starts[$e+1];
	    $current_feature->push_attribute("introns", $intron);
	}
    }

    if ($#exon_starts >= 0) {
	$current_feature->set_attribute("start_pos",&min(@exon_starts,@exon_ends));
	$current_feature->set_attribute("end_pos",&max(@exon_starts,@exon_ends));
    } else {
	$current_feature->set_attribute("start_pos",$null);
	$current_feature->set_attribute("end_pos",$null);
	&ErrorMessage(join ("\t", "Invalid position for feature", $feature_counts, $feature->get_attribute("id")), "\n");
    }
    $current_feature->set_attribute("strand",$strand);


    
    return();
} 


################################################################
## Export protein sequences in fasta format
sub ExportProteinSequences {
    my ($features, $org) = @_;
    $outfile{pp} = $dir{output}."/".$org."_aa.fasta";

    warn join ("\t", "; Exporting translated sequences to file", $outfile{pp}), "\n" if ($main::verbose >= 1);

    open PP, ">$outfile{pp}";
    foreach my $feature ($features->get_objects()) {
	next unless ($feature);
	if ($feature->get_attribute("type") eq "CDS") {
	    my ($translation) = $feature->get_attribute("translation");
	    next unless ($translation =~ /\S+/);

	    my $id = $feature->get_attribute("id");
	    my $pp_id = $id;
#	    my $pp_id = "id|".$id;
#	    $pp_id .= "|ref|".$feature->get_attribute("id");
#	    $pp_id .= "|name|".$feature->get_attribute("name");
#	    $pp_id .= "| ";

	    my $description;
	    $description .= $feature->get_attribute("description");
	    $description .= "; ".join ("|", $id, $feature->get_attribute("names"));


	    print PP $header, "\n";
#	    &PrintNextSequence(PP,"fasta",60,$translation,$pp_id);
	    &PrintNextSequence(PP,"fasta",60,$translation,$pp_id, $description);
	}
    }
    close PP;
}


################################################################
## Create the organism object if it has not ben created
sub CreateOrganism {
  my ($organisms, $organism_name) = @_;
  unless ($organism_created{$organism_name}) {
    &RSAT::message::Info("Creating organism Organism name", $organism_name) if ($main::verbose >= 1);
    $organism = $organisms->new_object(%args);
    $organism->push_attribute("names", $organism_name);
    $organism->set_attribute("name",$organism_name);
    $organism_created{$organism_name} = 1;
  }
}
