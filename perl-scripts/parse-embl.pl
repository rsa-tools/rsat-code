#!/usr/bin/perl 
############################################################
#
# $Id: parse-embl.pl,v 1.9 2004/03/29 12:41:47 jvanheld Exp $
#
# Time-stamp: <2003-10-21 01:17:49 jvanheld>
#
############################################################
#use strict;;
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}
#require "RSA.lib";
#push @INC, $ENV{PARSER};
#require "PFBP_classes.pl";
#require "PFBP_parsing_util.pl";

require "RSA.lib";
push @INC, "$RSA/perl-scripts/parsers/";
require "lib/load_classes.pl";
#require "lib/util.pl";
require "lib/parsing_util.pl";
#require "classes/Genbank_classes.pl";


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
			     position=>"SCALAR",
			     chromosome=>"SCALAR",
			     strand=>"SCALAR",
			     start_pos=>"SCALAR",
			     end_pos=>"SCALAR",
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
    my $start_time = &AlphaDate();

    ### initial directory
    $dir{main} = `pwd`; #### remember working directory
    chomp($dir{main});

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
    

#    my $genes = classes::ClassFactory->new_class(object_type=>"EMBL::Gene",
#						 prefix=>"gene_");

    ################################################################
    #### Open class factories

    #### Organism
    my $organisms = classes::ClassFactory->new_class(object_type=>"EMBL::Organism",prefix=>"org_");
    $organisms->set_out_fields(qw (id
				   name
				   names
				   taxonomy
				   ));
    

    #### Contig
    my $contigs = classes::ClassFactory->new_class(object_type=>"EMBL::Contig",prefix=>"ctg_");
    $contigs->set_out_fields(qw (id
				 accession
				 version
				 form
				 type
				 file
				 length
				 description
				 ));
    
    
    #### Features
    my $features = classes::ClassFactory->new_class(object_type=>"EMBL::Feature",prefix=>"ft_");
    @feature_out_fields = qw( id
			      type
			      name
			      contig
			      start_pos
			      end_pos
			      strand
			      description
			      organism
			      chromosome
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
    @classes = qw ( EMBL::Feature
		  EMBL::Contig
		  EMBL::Organism
		    );
    
    #### working directory
    $wd = `pwd`;
    chomp($wd);
    
    &ReadArguments();
    
    ################################################################
    #### check argument values ####
    
    #### input directory
    unless (defined($dir{input})) {
	&FatalError("You must specify the input directory.\n");
    }
    unless (-d $dir{input}) {
	&FatalError("Input directory '$dir{input}' does not exist.\n");
    }

    #### First guess for organism
    unless (defined($org)) {
	$org = `basename $dir{input}`;
	chomp($org);
	warn "; Auto selection of organism name\t$org\n" if ($verbose >= 1);
    }

    #### find embl files in the input directory
    chdir ($dir{input});
    @embl_files =();
    push @embl_files, glob("*.embl");
    push @embl_files, glob("*.embl.Z");
    if ($#embl_files < 0) {
	&FatalError("There is no embl file in the input directory $dir{input}\n");
    } else {
	warn "; EMBL files\n;\t", join("\n;\t", @embl_files), "\n" if ($verbose >= 1);
    }
    #### come back to the starting directory
    chdir($wd);

    #### output directory
    unless (defined($dir{output})) {
	#### default output directory
	$export_subdir = "embl";
	$dir{output} = "$parsed_data/$export_subdir/$delivery_date/$org";
#	$dir{output} = "$RSA/data/embl_genomes/$org/genome";
	warn "; Auto selection of output dir\t$dir{output}\n" if ($verbose >= 1);
    }
    &CheckOutputDir($dir{output});
    chdir $dir{output};
    $out_file{error} = "$dir{output}/genbank.errors.txt";
    $out_file{stats} = "$dir{output}/genbank.stats.txt";
    open ERR, ">$outfile{error}";

    #### come back to the starting directory
    chdir($wd);

    #### verbose ####
    &Verbose() if ($verbose);

    #### parse the embl files
    $contig_handle = &OpenOutputFile("$dir{output}/contigs.txt"); # file with chromosome IDs
    foreach my $file (@embl_files) {
	my $sequence = &ParseEMBLFile("$dir{input}/$file", 
				      $organisms, 
				      $contigs,
				      $features, 
				      source=>$file);

	unless ($noseq) {
	    my $seq_file = $contig->get_attribute("id");
	    $seq_file .= ".raw";
	    open RAW, ">$dir{output}/$seq_file";
	    &PrintNextSequence(RAW, "raw", 0, $sequence, $contig);
	    close RAW;
	    print $contig_handle $seq_file, "\t", $contig->get_attribute("id"), "\n";
	}
    }
    close $contig_handle;


    ################################################################
    #### Postprocessing : we interpret the features to extract more information
    #### check some feature attributes (name, description, ...)

    warn "; Processing features\n" if ($verbose >= 1);
    my $organism_name = $org; #### first guess
    my $chromosome = "";
    foreach $feature ($features->get_objects()) {

	#### the real organism name is annotated in the feature of type "source"
	if ($feature->get_attribute("type") eq "source") {
            $organism_name = $feature->get_attribute("organism");
            $chromosome = $feature->get_attribute("chromosome");
            if (($chromosome eq "") || ($chromosome eq $null)) {
		$chromosome = $feature->get_attribute("contig");
		$feature->set_attribute("chromosome", $chromosome);
	    }
        } else {
           $feature->set_attribute("organism", $organism_name);
	   $feature->set_attribute("chromosome", $chromosome);
        }
	
        ### use gene attribute as name
	foreach my $name ($feature->get_attribute("gene")) {
	    if ($feature->get_attribute("type") eq "CDS") {
		$feature->push_attribute("names",$name);
	    } else {
		$feature->push_attribute("names",$name." ".$feature->get_attribute("type"));
	    }
	}

	### add protein_id as valid name
	@protein_ids = $feature->get_attribute("protein_id");
	if ($#protein_ids >= 0) {
#	    $feature->force_attribute("id", $protein_ids[0]); #### Use protein ID as ID for the CDS
	    foreach my $id (@protein_ids) {
		$feature->push_attribute("names",$id);
	    }
	}
	

	#### add SWISS-PROT ID as valid name
	my @xrefs = $feature->get_attribute("db_xref");
	my $gi = "";
	foreach my $xref (@xrefs) {
	    if (($xref =~ /SWISS-PROT:/) ||
		($xref =~ /SPTREMBL:/)) {
		my $name = $';
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
	$features->set_attribute_header("db_xref_exp", "xref_db", "xref_id");
	my $gi = "";
	foreach my $xref (@xrefs) {
	    if ($xref =~ /:/) {
		my $xref_db = $`;
		my $xref_id = $';
		$feature->push_expanded_attribute("db_xref_exp", $xref_db, $xref_id);
	    } else {
		&ErrorMessage(join ("\t", 
				    "Invalid cross reference", 
				    $feature->get_attribute("id"), 
				    $xref), "\n");
	    }
	}
	

#  	#### use GI as feature identifier
#  	my @xrefs = $feature->get_attribute("db_xref");
#  	my $gi = "";
#  	foreach my $xref (@xrefs) {
#  	    if ($xref =~ /GI:/) {
#  		$gi = $';
#  		last;
#  	    } 
#  	}
#  	if ($gi) {
#  	    $feature->force_attribute("id",$gi);
#  	} else {
#  	    &Warning("feature\t".$feature->get_attribute("id")."\thas no GI.\n") if ($verbose >= 2); 
#  	}
	
#  	#### use embl name as chromosome name
#  	my $source = $feature->get_attribute("source");
#  	if ($source =~ /embl:/) {
#  	    my $chromosome = $';
#  	    $chromosome =~ s/\.gz$//;
#  	    $chromosome =~ s/\.contig$//;
#  	    $feature->force_attribute("chromosome",$chromosome);
#  	}
    }
    
    ################################################################
    ### Save result in tab files
    chdir $dir{main};
#    chdir $dir{output};
    foreach $class_factory ($organisms, $contigs, $features) {
	$class_factory->dump_tables();
	$class_factory->generate_sql(dir=>"$dir{output}/sql_scripts",
				schema=>$schema,
				host=>$host,
				user=>$user,
				password=>$password
				);
    }
    &ExportMakefile(@classes);
    &PrintStats($outfile{stats}, @classes);

    #### document the ffields which were parsed but not exported
    chdir($dir{output});
    open STATS, ">>$outfile{stats}";
    print STATS "; \n; Fields parsed but not exported\n";
    foreach my $key (sort keys %feature_extra_fields) {
	my $value = $feature_extra_fields{$key};
	printf STATS ";\t%-20s\t%d\n", $key, $value;
    }
    close STATS;

    ###### verbose ######
    if ($verbose) {
	my $done_time = &AlphaDate();
	print $out "; Job started $start_time\n";
	print $out "; Job done    $done_time\n";
    }

    ###### close output file ######
    close $out if ($outfile{output});
    close ERR;

    warn "; Results exported in directory\t", $dir{output}, "\n" if ($verbose >= 1);

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

        2001 by Jacques van Helden (jvanheld\@ucmb.ulb.ac.be)
	
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
	-o	output directory
		The parsing result will be saved in this directory. If
		the directory does not exist, it will be created.
	-org	organism name
	-noseq  do not export sequences in .raw files

   Options for the automaticaly generated SQL scripts
	-schema database schema (default: $schema)
	-host	database host (efault: $host)
	-user	database user (efault: $user)
	-password	
		database password (default: $password)

INPUT FILE

        This parser takes as input all the .embl and .embl.Z files
        found in the input directory. Each file is parsed and the
        results are exported in a single directory (the outpu
        directory). Features found in different files are merged in
        the same output file (features.tab). The idea is tha one
        parsing corresponds to one genome.

OBJECT TYPES

        The parser is based on a syntactic decomposition of the .embl
        files. Internally, it creates the following objects

	    EMBL::Organism
	    EMBL::Contig
	    EMBL::Feature

	Each feature has a type (CDS, tRNA, rRNA, promoter, ...), as
	found in the original EMBL file. The feature type si indicated
	by a line starting with FT, followed by 3 white spaces, the
	type, some additional spaces, and the position.

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
-v	verbose
-org	organism name
-schema database schema (default: $schema)
-host	database host (default: $host)
-user	database user (default: $user)
-password	database password (default: $password)
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
		$verbose = $ARGV[$a+1];
	    } else {
		$verbose = 1;
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

	    ### output file ###
	} elsif ($ARGV[$a] eq "-o") {
	    $dir{output} = $ARGV[$a+1];
	    
	    ### organism ###
	} elsif ($ARGV[$a] eq "-org") {
	    $org = $ARGV[$a+1];
	    $org =~ s/\s+/_/g;

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
	} elsif ($ARGV[$a] eq "-password") {
	    $password = $ARGV[$a+1];
	    
	}
    }
}


################################################################
#### print verbose message
sub Verbose {
    print $out "; parse-embl ";
    &PrintArguments($out);
    if (defined(%dir)) {
	print $out "; Directories\n";
	while (($key,$value) = each %dir) {
	    print $out ";\t$key\t$value\n";
	}
    }
    if (defined(%infile)) {
	print $out "; Input files\n";
	while (($key,$value) = each %infile) {
	    print $out ";\t$key\t$value\n";
	}
    }
    if (defined(%outfile)) {
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
    my ($input_file, $organisms, $contigs, $features, %args) = @_;
    warn ";\n; Parsing file $input_file\n" if ($verbose >= 1);
    
    my ($file,$dir) = &OpenInputFile($input_file);
#    open EMBL, $input_file 
#	|| die "Error: cannot open input file $input_file.\n";

    #### initialize parsing variables
    my $l = 0; #### line counter
    my $in_FT = 0; #### flag to indicate whether we are reading a feature
#    my $in_feature = 0;
#    my $in_cds = 0;

    my $current_feature = null;
    my $sequence = "";
    my $organism_name = "";
    my $organism;

    #### parse the file
    while (my $line = <$file>) {
	$l++;
	if (($test) && ($l > $test_lines)) {
	    warn "Test: stopping at line $l\n";
	    last; 
	}

	warn $line if ($verbose >= 10);
	chomp $line;
	next unless ($line =~ /\S/);

	#### read the full sequence
	if  ($line =~ /^SQ/) {
	    $in_FT = 0;
	    $in_sequence = 1;
	    while (my $line = <$file>) {
		if ($line =~ /\d+\s*$/) {
		    $sequence .= $`;
		} elsif ($line =~ /^\/\/$/) {
		    $in_sequence = 0;
		}
	    }
        }
	
	if  ($line =~ /^ID\s+(\S+)\s*/) {
	    #### contig ID line
	    $contig = $contigs->new_object(id=>$1);
	    my $contig_description = $';
	    $contig->set_attribute("description",$contig_description);
	    $contig->set_attribute("file",$input_file);

	    #### contig length
	    if ($contig_description =~ /(\d+)\s+BP\./i) {
		$contig->set_attribute("length", $1);
	    }

	    #### contig form
	    if ($contig_description =~ /circular/i) {
		$contig->set_attribute("form", "circular");
	    } else {
		$contig->set_attribute("form", "linear");
	    }

	} elsif ($line =~ /^AC\s+/) {
	    #### accession number for the currrent contig
	    my $AC = $';
	    $AC =~ s/;$//;
	    $contig->set_attribute("accession", $AC);

	} elsif ($line =~ /^SV\s+/) {
	    #### sequence version for the currrent contig
	    my $version = $';
	    $contig->set_attribute("version", $version);

	} elsif ($line =~ /^OS\s+/) {
	    #### organism name
	    $organism_name = $';
	    warn "; Organism name\t", $organism_name, "\n" if ($verbose >= 2);
	    $organism = $organisms->new_object(%args);
	    $organism->push_attribute("names", $organism_name);
	    $organism->set_attribute("name",$organism_name);

	} elsif ($line =~ /^OC\s+/) {
	    my $taxonomy = $';

	    #### read the rest of the taxxonomy (it can be larger than
	    #### one line)
	    while ($line = <$file>) {
		chomp $line;
		$l++;
		if ($line =~ /^OC\s+/) {
		    $taxonomy .= $';
		} else {
		    last;
		}
	    }
	    
	    #### organism taxonomy
	    if ($organism) {
		$organism->set_attribute("taxonomy", $taxxonomy);
	    }
	} 


	#### start reading features
	unless ($in_FT) {
	    if ($line =~ /^FT\s+/) {
		$in_FT = 1 ;
		warn "; Reading features\n" if ($verbose >= 10);
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

	    #### create a new object for the feature
	    $current_feature = $features->new_object(%args);
	    $feature_count++;
	    warn join "\t",  "; file $input_file", "line $l", "feature $feature_count", $feature_type, "\n" if ($verbose >= 2);
	    $current_feature->set_attribute("type",$feature_type);
	    $current_feature->set_attribute("contig",$contig->get_attribute("id"));
#	    $current_feature->set_attribute("organism",$organism_name);
	    $current_feature->set_attribute("position",$position);
	    &ParsePositionEMBL($current_feature);

	    ### chromosome
#	    my $chromosome = $current_feature->get_attribute("source");
#	    $chromosome =~ s/\.embl.*//;
#	    $current_feature->set_attribute("chromosome",$chromosome);
	}



	#### new feature attribute
	if ($line =~ /^FT    \s+\/(\S+)\=/) {
	    my $key = $1;
	    #### check that this attribute is part of the exported ones
	    unless ($feature_out_fields{$key}) {
		$feature_extra_fields{$key}++;
	    }

	    my $value = $';
	    my $start_l = $l;

	    #### If the attribute starts with a quote, make sure to
	    #### have the closing quote, and trim the quotes
	    if  ($value =~ /^\"/) {
		$value = $'; #### trim the leading quote

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
			    $value .= " ".$';
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
	    $current_feature->new_attribute_value($key, $value);
	    warn join ("\t", "; attribute", $feature_count, $start_l, $l, $key, $value), "\n" if ($verbose >= 4);

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
	} else {
	    $strand = "D";
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
	$current_feature->set_attribute("start_pos",$exon_starts[0]);
	$current_feature->set_attribute("end_pos",$exon_ends[$#exon_ends]);
    } else {
	$current_feature->set_attribute("start_pos",$null);
	$current_feature->set_attribute("end_pos",$null);
	&ErrorMessage(join ("\t", "Invalid position for feature", $feature_counts, $feature->get_attribute("id")), "\n");
    }
    $current_feature->set_attribute("strand",$strand);


    
    return();
} 




