#!/usr/bin/perl
############################################################
#
# $Id: parsing_util.pl,v 1.10 2004/01/23 16:25:50 oly Exp $
#
# Time-stamp: <2003-10-01 17:00:56 jvanheld>
#
############################################################
### util.pl
### utilities for the AMAZE project

require "lib/load_classes.pl";
require "lib/util.pl";
require Data::Dumper;

#### check output directory
sub CheckOutputDir {
    #### special dir for quick tests
    if ($main::test) {
	$dir{output} .= "_test";
    }
    unless (-d $dir{output}) {
	warn "; Creating output dir $dir{output}\n" if ($main::verbose >= 1);
	`mkdir -p $dir{output}`;
#	mkdir $dir{output}, 0775 || die "Error: cannot create directory $dir\n";
	unless (-d $dir{output}) {
	    die "ERROR\tCannot find directory $dir{output}\n";
	}
    }
    if ($clean) {
	system "\\rm -rf $dir{output}/*";
	system "\\rm -rf $dir{output}/sql_scripts/*";
    }
#    chdir $dir{output};
}

#### copy the required files to the delivery directory
sub deliver {
    my ($delivery_dir) = $_[0];
    unless ($delivery_target) {
	$delivery_target = "jvanheld\@paulus.ulb.ac.be:/rubens/dsk3/genomics/amaze/parsed_data/${delivery_dir}";
	
	$delivery_target =~ s|//|/|g;
    }
    
    my $delivery_source = $dir{output};
    $delivery_source =~ s|/+$||;
    unless ($delivery_source) {
	warn "ERROR: empty delivery source\n";
	return;
    }

    warn "; Delivering parsed data to directory\n;\t$delivery_target.\n" if ($verbose >=1); 
    $command = "rsync -e ssh -ruptvl $delivery_source $delivery_target";
    warn "; $command\n" if ($verbose >=1);
    system $command;
}


################################################################
### Export all classes from a list provided as argument
### usage: 
###    &ExportClasses($output_file, $out_format, @classes);
sub ExportClasses {
  local ($output_file, $out_format, @classes) = @_;
  if ($output_file) {
    open STDOUT, ">$output_file"  || die "Error : cannot write file $output_file\n"; 
  }
  foreach $class (@classes) {
    warn (";\n; ", &AlphaDate, " exporting class ", $class, " to file '$output_file'\n") 
      if ($verbose >= 1);
    local @selected = @{$out_fields{$class}};
    foreach $object ($class->get_objects()) {
      $object->print_attributes($out_format, @selected);
    }
    warn ("; ", &AlphaDate, " class ", $class, " exported\n") 
      if ($verbose >= 1);
  }
  close STDOUT if ($output_file);
  return 1;
}




################################################################
## Export a makefile which will call other makefiles for the loading
###    &ExportMakefile(@classes);
sub ExportMakefile {
    my @classes = @_;

    ## data loading
    foreach my $dbms (keys %main::supported_dbms) {
	my $makefile=$dir{output}."/sql_scripts/".$dbms."/makefile";
	warn "; Exporting makefile for $dbms in file $makefile\n" if ($main::verbose >= 2);
	open MAKEFILE, ">$makefile" || die "Cannot write makefile in dir $dir\n";
	#### usage
	print MAKEFILE "\nusage:\n";
	print MAKEFILE "\t", '@perl -ne \'if (/^([a-z]\S+):/){ print "\t$$1\n";  }\' ', "makefile\n";
        foreach my $target ("create", "uncompress", "load", "recompress", "all", "drop") {
	    print MAKEFILE "${target}:\n";
	    my $usage_done = 0;
	    for my $class (@classes) {
		(my $short_class = $class) =~ s/.*:://g;
		print MAKEFILE "\tmake -i -f ${short_class}.mk ${target}\n";
	    }
	}
	close MAKEFILE;
    }
}


################################################################
### parse a file in KEGG format (flat file)
### see KEGG ftp site for format specifications
### usage : $class_holder = &ParseKeggFile($kegg_file,
###                                        $class_holder,
###                                        %default_args)
sub ParseKeggFile {
    my ($kegg_file,$class_holder, %args) = @_;
    my $key = "";
    my $value = "";
    my $id = "";
    my $in;
    my %parse_fields = ();
    my $class = $class_holder->get_object_type();
    my %attr_types = $class->get_attribute_cardinalities();

    my $parse_all = 0;
    my @out_fields = $class_holder->get_out_fields();
    warn ("; Output fields\t", join (";", @out_fields), "\n") if ($verbose >=1);

    if ($#out_fields <0) {
	$parse_all = 1;
    } else {
	foreach $field (@out_fields) {
#    foreach $field (@{$out_fields{$class}}) {
	    $parse_fields{$field} = 1;
	}
    }
    
    
    ### open the input stream
    if ($kegg_file) {
	warn (";\n; ", &AlphaDate,
	      " parsing class $class from $kegg_file\n")
	    if ($verbose >= 1);
	open INFILE, $kegg_file 
	    || die "Error: cannot $kegg_file\n";
	$in = INFILE;
    } else {
	warn (";\n; ", &AlphaDate,
	      " parsing class $class from STDIN\n")
	    if ($verbose >= 1);
	$in = STDIN;
    } 
    
    ### read data
    my @names = ();
    my $l = 0;
    while (<$in>) {
	$l++;
	chomp;
	if (/\/\/\//) { ### end of record
	    if (defined($object)) {

		#### special treatment for yeast genes, where the ID can be considered as a name
		if (($class eq "classes::Gene") && 
		    ($object->get_attribute("organism") =~ /S.*cerevisiae/)) {
		    $object->push_attribute("names", $object->get_id());
		}
	    } else {
		&ErrorMessage("Error: record separator at line $l without any defined object");
	    }
	    undef($object);
	    next;
	}
	
	### split the line into field key and value
	$key = lc(substr($_,0,12));
	$value = substr($_,12);
    
	if ($key =~ /^(\S+)/) { ### new field
	    $current_key = $1; ### trim the white spaces
	    $current_key =~ s/^name$/names/i; 
	    $current_key =~ s/^definition$/description/i;
	    $current_key =~ s/^position$/chrom_position/i;
	    $current_key =~ s/^source$/ext_source/i;
	    if ($current_key =~ /entry/) { ### new entry
		my %gene_args = ();

		### specific re-id for KEGG genes, when their ID is ambiguous
		if ($class eq "classes::Gene") {
		    my @fields = split /\s+/, $value;
		    my $entry_id = $fields[0];
		    my $entry_type = $fields[1];
		    my $entry_organism = $fields[2];
		    $gene_args{type} = $entry_type;
		    $gene_args{organism} = $entry_organism;
#		    if ($main::rsa) {
			$current_id = $entry_id;
#		    } else {
#			$current_id = $entry_organism."_".$entry_id;
#		    }
		    warn "gene\t$current_id\n" if ($verbose >= 4);
		    $current_id =~ s/\s+/ /g;
		    $current_id =~ s/\s/_/g;
		} elsif ($value =~ /^(\S+)/) {
		    $current_id = $1;
		}	  
		
		
		### instantiate new object
		warn ";creating object\t$current_id in class\t$class\n" 
		    if ($verbose >= 2);
		$object = $class_holder->new_object(%args, 
						    %gene_args,
						    id=>$current_id);
	    } elsif ($current_key eq "names") {
		#### for KEGG genes only, split the string, which may contain multiple names
		if ($class_holder->get_object_type() eq "KEGG::Gene") {
		    my @new_names = split /\, +/, $value;
		    foreach my $name (@new_names) {
			$object->new_attribute_value($current_key,$name);
		    }
		} else {
		    $object->new_attribute_value($current_key,$value);
		}

	    } elsif ($current_key eq "chrom_position") {
		$object->new_attribute_value($current_key,$value);

	    } elsif ($current_key) {
		if (($parse_all) || 
		    ($parse_fields{$current_key})) { ### field to parse
		    ### set new attribute value
		    $object->new_attribute_value($current_key,$value);
		}
	    }
	    
	    ### continuation of the same field as previous line
	} elsif (($parse_all) || 
		 ($parse_fields{$current_key})) {
	    if (($value =~ /^\$(.+)/) && ($current_key)) {    ### continuation from the previous line
		### concatenate the new value the attribute previously stored
		$object->append_attribute($current_key,$1);
	    } elsif ($attr_types{$current_key} eq "SCALAR") { ### single-value attribute
		### concatenate the new value the attribute previously stored
		$object->append_attribute($current_key," ".$value);
	    } else { ### new value for a multiple-value field
		### push a new attribute value to the list
		$object->new_attribute_value($current_key,$value);
	    }
	}
	
    }
    
    
    ### create indexes
    #$class->index_ids();
    #$class->index_names();
    
    ### create indexes
    #$class_holder->index_ids();
    #$class_holder->index_names();
    
    close $in if ($kegg_file);

    return $class_holder;
}


################################################################
#### print command-line arguments
sub ArgString {
    my $arg_string = "";
    foreach my $a (@main::ARGV) {
	if (($a =~ /\s+/)  ||
	    ($a !~ /\S+/) ||
	    ($a =~ /[\(\)\>\<\&]/)) {
	    $arg_string .=  " '$a'";
	} else {
	    $arg_string .= " $a";
	}
    }
    return $arg_string;
}

################################################################
#### prints the number of objects and the number of attribute allocations
#### for each selected class 
#### usage : &PrintStats($file,@classes);
sub PrintStats {
  my ($stat_file, @classes) = @_;
  
  if ($stat_file) {
    open STATS, ">$stat_file" || die "Error: cannot write stat file $stat_file\n";
    $stats = STATS;
  } else {
    $stats = STDOUT;
  }
  
  ### print some stats
  $alpha_date =  &AlphaDate();
  warn (";\n; ", $alpha_date,
	" printing parsing statistics in file ", 
	$stat_file, "\n")  
      if ($verbose >= 1);
  
  printf $stats "; Parsing date\t$alpha_date\n";
  $program = join (" ", $main::0, &ArgString());
  printf $stats "; Parsing command\t",$program, "\n";
  printf $stats "; stats\n";
  foreach $class (@classes) {
    printf $stats ";\n; class\t$class\n";
    if ($class->can("get_count")) {
      my $count =  $class->get_count();
      printf $stats ";\t%d entries\n", $count;
      if ($count) {
	### print attribute counts
	print $stats ";\tAttributes\n";
	#@fields = @{$out_fields{$class}} || $class->get_attribute_names;
	%attr_counts = $class->get_attribute_counts();
	%attr_types = $class->get_attribute_cardinalities();
	
	foreach $attribute ($class->get_attribute_names()) {
	  $count = $attr_counts{$attribute};
	  $type = $attr_types{$attribute};
	  printf $stats ";\t\t%-15s\t%d\t%s\n", $attribute, $count, $type;
	}
      }
    }
    
  }
  
  close STATS if ($stat_file);
}

################################################################
### split a side of an equation into its reactants, and for each,
### separates the stoichiometry from the compound
sub SplitReactants {
  my @reactants = split / \+ /, $_[0];
  my %reactants;
  my $stoichio = "";
  my $compound = "";
  my $compound = "";
  
  foreach $reactant (@reactants) {
    if ($reactant =~ /^(.*)\s*(C\d+)/) {
      $stoichio = $1 || 1;
      $stoichio = &trim($stoichio);
      $compound = $2;
      $reactants{$compound} = $stoichio;
    }
  }
  return %reactants;
}


################################################################
#
# After having parsed genes, CDS and RNA from Genbank, create one
# feature for each CDS and another one for each RNA. This feature
# summarizes information (only retins selected fields) and reformats
# it (multiple gene names) for ease of manipulation.
sub CreateGenbankFeatures {
    my ($features, $genes, $mRNAs, $tRNAs, $rRNAs, $misc_RNAs, $misc_features, $CDSs, $sources) = @_;

    #### warning message
    warn ("; ",
	  &AlphaDate(),
	  "\tChecking features found in Genbank file\n")
	if ($verbose >= 1);



    #### assign gene name to the different types of features
    foreach my $source ($sources->get_objects()) {
	$source->get_taxid();
    }
    
    foreach my $object ($features->get_objects(),
			$genes->get_objects(),
			$mRNAs->get_objects(),
			$tRNAs->get_objects(),
			$rRNAs->get_objects(),
			$misc_RNAs->get_objects(),
			$misc_features->get_objects(),
			$CDSs->get_objects()
			) {
	if ($object->get_attribute("gene")) {
	    #### add gene name to the list of synonyms
	    $object->push_attribute("names", $object->get_attribute("gene"));
	    #### define this name as primary name
	    $object->set_attribute("name",  $object->get_attribute("gene"));
	}
    }
    
    ################################################################
    #### this only applies to mRNA and CDS because the other features
    #### (gene, tRNA, rRNA, misc_RNA, misc_feature) have no IDs in Genbank flat files !!!!
    foreach my $object ($mRNAs->get_objects(),
			$CDSs->get_objects()
			) {
	$object->UseGenbankGIasID();
    }
    
    
    #### initialize parameters
    my $gi_as_id = 1;
    warn "; Creating features from CDS and RNAs\n" if ($verbose >= 1);
    foreach my $parsed_feature ($CDSs->get_objects(),
				$mRNAs->get_objects(),
				$tRNAs->get_objects(),
				$rRNAs->get_objects(),
				$misc_RNAs->get_objects(),
				$misc_features->get_objects()
				) {

	

	$created_feature = $features->new_object(%args);
	
	$created_feature->set_attribute("type",$parsed_feature->get_attribute("type"));
	$created_feature->set_attribute("organism",$parsed_feature->get_attribute("organism"));
	$created_feature->set_attribute("contig",$parsed_feature->get_attribute("contig"));
	$created_feature->set_attribute("chrom_position",$parsed_feature->get_attribute("chrom_position"));
	$created_feature->set_attribute("start_pos",$parsed_feature->get_attribute("start_pos"));
	$created_feature->set_attribute("end_pos",$parsed_feature->get_attribute("end_pos"));
	$created_feature->set_attribute("strand",$parsed_feature->get_attribute("strand"));

	warn join ("\t", "feature", 
		   $created_feature->get_attribute("id"),
		   "type", $created_feature->get_attribute("type"),
		   "organism", $created_feature->get_attribute("organism"),
		   ), "\n" if ($verbose >= 2);


	################################################################
	#### Define names for the new feature
	warn join ("\t",  "feature",		   
		   $created_feature->get_attribute("id"),
		   "Adding gene names", 
		   ), "\n" if ($verbose >= 5);
	my $gene_name = "";
	my $gene_id = "";
	my $gene = "";
	if ($parsed_feature->get_attribute("gene")) {
	    $gene_name = $parsed_feature->get_attribute("gene");
	    $gene_id = $parsed_feature->get_attribute("gene_id");
	    warn join ("\t", "feature", $created_feature->get_attribute("id"),
		       "gene", 
		       "id", $gene_id,
		       "name", $gene_name,
		       ), "\n" if ($verbose >= 5);
	    
	    $created_feature->push_attribute("names",$gene_name);

	    $gene = $genes->get_object($gene_id);
	    if ($gene) {
		#### use gene name as primary name for the feature
		$created_feature->set_attribute("name", $gene->get_name());
		$created_feature->push_attribute("names", $gene->get_name());

	    	#### accept all gene names for the new feature
		foreach my $name ($gene->get_attribute("names")) {
		    $created_feature->push_attribute("names",$name) unless ($names{$name});
		}
	    }
	}

	#### add names from the parsed feature
	warn join ("\t",  "feature",		   
		   $created_feature->get_attribute("id"),
		   "Adding names from original feature", $parsed_feature->get_attribute("id")
		   ), "\n" if ($verbose >= 5);
	foreach my $name ($parsed_feature->get_attribute("names")) {
	    $created_feature->push_attribute("names",$name) unless ($name eq $gene_name);
	}
	
	################################################################
	#### cross-references
	warn join ("\t",  "feature",		   
		   $created_feature->get_attribute("id"),
		   "Cross references", 
		   ), "\n" if ($verbose >= 5);
	my @xrefs = $parsed_feature->get_attribute("db_xref");
	my $gi = "";
	foreach my $xref (@xrefs) {
	    #### extract GI from cross-references
	    if ($xref =~ /GI:/) {
		$gi = "$'";
		last;
	    } 
	}
	
	warn join ("\t", ";", 
		   "parsed feature", $parsed_feature->get_attribute("id"), 
		   "GI=$gi", 
		   "product=$products[0]",
		   ), "\n" if ($verbose >= 5);
	
	if ($gi) {
	    if ($gi_as_id) {
		$created_feature->force_attribute("id",$gi); #### use GI as dientifier
	    } else {
		$created_feature->push_attribute("names",$gi); #### accept GI as synonym
	    }
	} else {
	    &ErrorMessage("; Error\tfeature ".$created_feature->get_attribute("id")." has no GI.\n"); 
	}
	
	################################################################
	#### define a single name  (take the first value in the name list)
	if ($single_name) {
	    if ($name = $created_feature->get_name()) {
		$created_feature->force_attribute("name",$name);
	    } else {
		$created_feature->force_attribute("name",$created_feature->get_id());
	    }
	    warn join ("\t", "feature",		   
		       $created_feature->get_attribute("id"),
		       "single name", $name,, 
		       ), "\n" if ($verbose >= 5);
	}
	
	################################################################
	#### create a description for the new feature
	my @products = $parsed_feature->get_attribute("product");
	if (($products[0]) && ($products[0] ne $null)){
	    $created_feature->set_attribute("description",$products[0]);
	} else {
	
	    #### if there is no product field, use the parsed feature note
	    my @feature_notes = $parsed_feature->get_attribute("note");
	    if ($feature_notes[0]) {
		$created_feature->set_attribute("description",$feature_notes[0]);
	    } else {
		#### if there is no feature not, use the gene note
		my @gene_notes = $gene->get_attribute("note");
		if ($gene_notes[0]) {
		    $created_feature->set_attribute("description",$gene_notes[0]);
		}
	    }
	}
	warn join ("\t", "feature",		   
		   $created_feature->get_attribute("id"),
		   "description", $created_feature->get_attribute("description"), 
		   ), "\n" if ($verbose >= 5);
	
	################################################################
	#### protein ID (for CDS)
	if ($parsed_feature->get_attribute("type") eq "CDS") {
	    my @protein_ids= $parsed_feature->get_attribute("protein_id");
	    foreach my $protein_id (@protein_ids) {
		if (($protein_id) && ($protein_id ne $null)){
		    ### remove the version number
		    if ($protein_id =~ /(.*)\.\d+$/) {
			$no_version = $1;
			$created_feature->push_attribute("names", $no_version);
			$parsed_feature->push_attribute("names", $no_version);
		    } else {
			$created_feature->push_attribute("names", $protein_id);
			$parsed_feature->push_attribute("names", $protein_id);
		    }
		}
	    }
	
	################################################################
	#### transcript ID (for mRNA)
	} elsif ($parsed_feature->get_attribute("type") eq "mRNA") {
	    my @transcript_ids= $parsed_feature->get_attribute("transcript_id");
	    foreach my $transcript_id (@transcript_ids) {
		if (($transcript_id) && ($transcript_id ne $null)) {
		    ### remove the version number
		    if ($transcript_id =~ /(.*)\.\d+$/) {
			$no_version = $1;
			$created_feature->push_attribute("names", $no_version);
			$parsed_feature->push_attribute("names", $no_version);
		    } else {
			$created_feature->push_attribute("names", $transcript_id);
			$parsed_feature->push_attribute("names", $transcript_id);
		    }
		}
	    }
	}

	################################################################
	#### delete duplicate names
	$created_feature->unique_names();
    }
	
}


################################################################
#### parse position string and add the different fields to each object
#### - chromosome
#### - strand
#### - start_pos
#### - end_pos
sub ParsePositions {
    my ($features) = @_;

    warn ("; ",
	  &AlphaDate(),
	  "\tparsing feature positions\n")
	if ($verbose >= 1);

    foreach my $feature ($features->get_objects()) {
	my $position = $feature->get_attribute("chrom_position");

	if ($position eq $null) {
	    $feature->set_attribute("chrom_position",$null);
#	    $feature->set_attribute("chromosome", $null);
	    $feature->set_attribute("strand",$null);
	    $feature->set_attribute("start_pos",$null);
	    $feature->set_attribute("end_pos",$null);
	    &ErrorMessage("Warning: feature ", $feature->get_attribute("id"), " has no attribute chrom_position\n");
	    next;
	}
	my $coord = $null;
	my $chomosome = $null;
	my $chrom_pos = $null;
	my $strand = $null;
	my $start_pos = $null;
	my $end_pos = $null;


	################################################################
	#### separate chromosome from chromosomal position
	if ($position =~ /^(\S+)\:(.*)/) {
	    $chromosome = $1;
	    $chrom_pos = $2;
	} elsif ($position =~ /^([^\:]*)$/) {
	    $chromosome = $feature->get_attribute("source");
	    $chrom_pos = $1;
	} else {
	    &ErrorMessage("Warning: invalid position",
			  "\t", $feature->get_attribute("id"),
			  "\t", $feature->get_attribute("organism"),
			  "\t", $position, "\n");
	    $feature->force_attribute("chromosome", $null);
	    $feature->force_attribute("strand",$null);
	    $feature->force_attribute("start_pos",$null);
	    $feature->force_attribute("end_pos",$null);
	    next;
	}

	################################################################
	#### direct or reverse strand
	if ($chrom_pos =~ /complement\((.*)\)/) {
	    $strand = "R";
	    $coord = $1;
	} else {
	    $strand = "D";
	    $coord = $chrom_pos;
	}

	################################################################
	#### split exons

	
	if ($coord =~ /^join\((.*)\)$/) { ### exons
	    #### multiple segments
	    my @exons = split ",", $1;
	    my @exon_starts = ();
	    my @exon_ends = ();

	    #### exon limits
	    foreach my $exon (@exons) {
		$exon =~ s/\s+//g;
		$feature->push_attribute("exons", $exon);
		my ($exon_start, $exon_end) = &segment_limits($exon);
		push @exon_starts, $exon_start;
		push @exon_ends, $exon_end;
	    }
	    
	    #### feature start and end 
	    $start_pos = $exon_starts[0] || $null;
	    $end_pos = $exon_ends[$#exons] || $null;

	    #### introns
	    my @introns = ();
	    for my $e (0..$#exon_starts - 1) {
		my $intron = $exon_ends[$e] + 1;
		$intron .= "..";
		$intron .= $exon_starts[$e+1] -1;
		$feature->push_attribute("introns", $intron);
	    }

	} else {
	    #### a single segment
	    ($start_pos, $end_pos) = &segment_limits($coord);
	}


	#### a single segment
#	if ($coord =~ /^([\>\<]{0,1}\d+)\.\.([\>\<]{0,1}\d+)$/) {
#	    $start_pos = $1;
#	    $end_pos = $2;
#
#	} elsif ($coord =~ /^join\((.*)\)$/) { ### exons
#	    #### multiple segments
#	    my @exons = split ",", $1;
#	    my @exon_starts = ();
#	    my @exon_ends = ();
#
#	    #### exon limits
#	    foreach my $exon (@exons) {
#		$exon =~ s/\s+//g;
#		$feature->push_attribute("exons", $exon);
#		if ($exon =~ /([\>\<]{0,1}\d+)\.\.([\>\<]{0,1}\d+)/) {
#		    push @exon_starts, $1;
#		    push @exon_ends, $2;
#		    
#		} elsif ($exon =~ /([\>\<]{0,1}\d+)/) {
#		    #### some segments in Genbank annotations are made of a single nucleotide !
#		    push @exon_starts, $1;
#		    push @exon_ends, $1;
#
#		} else {
#		    &ErrorMessage("Error feature\t",$feature->get_id(),"\tinvalid exon\t$exon\n");
#		}
#	    }
#	    
#	    #### feature start and end 
#	    $start_pos = $exon_starts[0] || $null;
#	    $end_pos = $exon_ends[$#exons] || $null;
#
#	    #### introns
#	    my @introns = ();
#	    for my $e (0..$#exon_starts - 1) {
#		my $intron = $exon_ends[$e] + 1;
#		$intron .= "..";
#		$intron .= $exon_starts[$e+1] -1;
#		$feature->push_attribute("introns", $intron);
#	    }
#	} else {
#	    &ErrorMessage("Warning : feature ",$feature->get_attribute("id"),"\tinvalid feature position $position\n");
#	    $strand = $null;
#	    $start_pos = $null;
#	    $end_pos = $null;
#	}
	$feature->force_attribute("chromosome", $chromosome);
	$feature->force_attribute("strand",$strand);
	$feature->force_attribute("start_pos",$start_pos);
	$feature->force_attribute("end_pos",$end_pos);
    }
}


################################################################
#### return start and end positions of a segment
sub segment_limits {
    my ($segment) = @_;
    my $segment_start = $null;
    my $segment_end = $null;

    #### start and end positions are different
    if ($segment =~ /^([\>\<]{0,1}\d+)\.\.([\>\<]{0,1}\d+)$/) {
	$segment_start = $1;
	$segment_end = $2;

    #### single residue position
    } elsif ($segment =~ /^([\>\<]{0,1}\d+)$/) {
	$segment_start = $1;
	$segment_end = $1;
    } else {
	&ErrorMessage ("Invalid segment format\t$segment");
    }
    return ($segment_start, $segment_end);
}

################################################################
#
# Parse a Genbank genome file. 
# Create one object per feature.
# The sequence is either returned (default), or directly stored 
# in a file, when this file is specified with the argument
#     seq_file=>myfile
# Direct storage is particularly useful for big genomes (e.g. 
# human), where the size of some chromosomes would crate an 
# "out of memory" error. 
# The option 
#    no_seq=>1
# prevents from parsing the sequence.
sub ParseGenbankFile {
    my ($input_file, 
	$feature_holder, 
	$gene_holder, 
	$mRNA_holder, 
	$tRNA_holder, 
	$rRNA_holder, 
	$misc_RNA_holder, 
	$misc_feature_holder, 
	$CDS_holder, 
	$contig_holder, 
	$organism_holder, 
	$source_holder, 
	%args) = @_;

    my %features_to_parse = (CDS=>1,
			     source=>1,
			     mRNA=>1,
			     tRNA=>1,
			     rRNA=>1,
			     misc_RNA=>1,
			     misc_feature=>1,
			     gene=>1
			     );

    if ($main::verbose >= 1) {
	warn ";\n; Parsing file $input_file.\n"; 
	warn ";\n; Features to parse\t", join (",", keys %features_to_parse), "\n";
    }

    open GBK, $input_file 
	|| die "Error: cannot open input file $input_file.\n";
    my $l = 0;
    my $in_sequence = 0;
    my $in_features = 0;
    my $in_feature = 0;
    my $in_gene = 0;
    my $in_cds = 0;
    my $current_feature = null;
    my $current_contig = null;    
    my $taxid = "";
    my $organism = "";
    my $organism_name = "";
    my $sequence = "";

    %contig_keys = ("ACCESSION"=>1,
#		    "LOCUS"=>1,
		    "DEFINITION"=>1,
		    "VERSION"=>1,
		    "DBSOURCE"=>1,
		    "COMMENT"=>1,
		    "PUBMED"=>1,
		    "MEDLINE"=>1,
		    );

    while (my $line = <GBK>) {
	$l++;
	print STDERR "$l\t$line" if ($main::verbose >= 10);
	chomp $line;
	next unless ($line =~ /\S/);
	unless (($in_features) || 
		($in_sequence)){

	    #### organism name
	    if ($line =~ /^\s+ORGANISM\s+/) {
		$organism_name = "$'";
		warn "; Organism name\t\t$organism_name\n" if ($main::verbose >= 1);
		$current_contig->set_attribute("organism", $organism_name);		
		$organism = $organism_holder->get_object($organism_name);
		if ($organism) {
		    warn "; Organism $organism already created\n" if ($main::verbose >= 3);
		} else {
		    $organism = $organism_holder->new_object();
		    $organism->push_attribute("names", $organism_name);
		    $organism_holder->index_names(); ### required to prevent creating several objects for the same organism
		}

		#### collect the taxonomy
		my $taxonomy = "";
		while ($line = <GBK>) {
		    if ($line =~ /^\S/) {
			$taxonomy = &trim($taxonomy);
			$taxonomy =~ s/\s+/ /g;
			$taxonomy =~ s/\.\s*$//;
#			$current_contig->set_attribute("taxonomy", $taxonomy);
			$organism->force_attribute("taxonomy", $taxonomy);
			last;
		    } else {
			$taxonomy .= $line;
		    }
		}
	    }
	    if ($line =~ /^([A-Z]+)\s+/) {
		$current_contig_key = $1;
		$current_contig_value = "$'";
		chomp $current_contig_value;
		$current_contig_value =~ s/\.\s*$//;
		if ($contig_keys{$current_contig_key}) {
		    warn "parsing\t$current_contig_key\t$current_contig_value\n" if ($verbose >= 4);
		    $current_contig->set_attribute(lc($current_contig_key), $current_contig_value);
		}
		
	    } elsif ($line =~ /^ {12}/) {
		### suite of the current contig key
		$current_contig_value .= " ".$';
		warn "parsing\t$current_contig_key\t$current_contig_value\n" if ($verbose >= 4);
		$current_contig->force_attribute(lc($current_contig_key), $current_contig_value);
	    }

	    if ($line =~ /^FEATURES/) {
		$in_features = 1 ;
		undef($current_contig_key);
		undef($current_contig_value);
		warn "; Reading features\n" if ($main::verbose >= 1);
		next;
	    }
	}

	if  ($line =~ /^LOCUS/) {
	    #### new contig
	    @fields = split /\s+/, $line;
	    warn join ("\t", "; New contig", $line), "\n" if ($main::verbose >= 1);
	    my $contig_id = $fields[1];
	    my $length = $fields[2];
	    my $type = $fields[4];      #### DNA or peptide
	    my $form = $fields[5];      #### circular or linear
	    my $taxo_group = $fields[6];      #### short taxonomic group
	    my $contig_date = $fields[7];      #### date of the contig 
	    $current_contig = $contig_holder->new_object(id=>$contig_id);
	    if ($input_file =~ /[^\/]+\.gbk/) {
		$current_contig->set_attribute("genbank_file", $&);
	    }
	    $current_contig->set_attribute("length", $length);
	    $current_contig->set_attribute("type", $type);
	    $current_contig->set_attribute("form", $form);
	    $current_contig->set_attribute("taxo_group", $taxo_group);
	    $current_contig->set_attribute("date", $contig_date);

	} elsif  ($line =~ /^BASE COUNT/) {
	    #### base count
	    #### currently ignored

	    #### read the full sequence
	} elsif  ($line =~ /^ORIGIN/) {
	    warn "; Reading sequence\n" if ($main::verbose >= 1);
	    $in_features = 0;
	    $in_sequence = 1;
	    if ($args{no_seq}) {
		#### skip the sequence
		while (<GBK>) {
		    if (/^\/\/$/) {
			$current_contig = null;
			last;
		    }
		}

		################################################################
		#### save the whole contig sequence in a file
	    } elsif ($args{seq_dir}) {
		my $seq_file = $current_contig->get_attribute("id").".raw";
		$current_contig->set_attribute("seq_dir", $args{seq_dir});
		$current_contig->set_attribute("file", $seq_file);
		
#		die join ("\t", $current_contig, $seq_file), "\n";
		
		warn ("; Storing sequence ",
		      $current_contig->get_attribute("id"), 
		      " in file $args{seq_dir}/$seq_file\n" )
		    if ($main::verbose >= 1);

		open SEQ, ">$args{seq_dir}/$seq_file" 
		    || die "Error: cannot write sequence file $args{seq_file}\n";
		while ($line = <GBK>) {
		    if ($line =~ /^\s+\d+\s+/) {
			$sequence = "$'";
			$sequence =~ s/\s//g;
			print SEQ $sequence;
		    } elsif ($line =~ /^\/\/$/) {
			print SEQ "\n";
			$in_sequence = 0;
			$current_contig = null;
			last;
		    } else {
			&ErrorMessage("Invalid sequence format, skipped\t$line\n");
		    }
		}
		close SEQ;
	    } else {
		#### load sequence in memory and return it
		while ($line = <GBK>) {
		    if ($line =~ /^\s+\d+\s+/) {
			$sequence .= "$'";
		    } elsif ($line =~ /^\/\/$/) {
			$in_sequence = 0;
			$current_contig = null;
			last;
		    } else {
			&ErrorMessage("Invalid sequence format, skipped\t$line\n");
		    }
		}
	    }
	
	#### new feature
	} elsif ($line =~ /^ {5}(\S+)\s+/) {
	    $feature_type = $1;
	    $value = "$'";

	    #### parse feature  position
	    $position = &trim($value);
	    if ($position =~ /join\(/){
		### check that the position is complete
		my $start_line = $l;
		my $next = "$'";;
		unless ($next =~ /\)/) {
		    do {
			$l++;
			die "Error: position starting at line $l is not terminated properly.\n"
			    unless $position_suite = <GBK>;
			$position_suite =~ s/^FT//;
			$position .= &trim($position_suite);
		    } until ($position =~ /\)/);
		}
	    }

	    #### create an object for the new feature
	    if ($features_to_parse{$feature_type}) {
		warn join "\t", $l, "new feature", $feature_type, $position, "\n" if ($main::verbose >= 2);

		if ($feature_type eq "gene") {
		    $holder = $gene_holder;
		} elsif ($feature_type eq "mRNA") {
		    $holder = $mRNA_holder;
 		} elsif ($feature_type eq "tRNA") {
 		    $holder = $tRNA_holder;
 		} elsif ($feature_type eq "rRNA") {
 		    $holder = $rRNA_holder;
 		} elsif ($feature_type eq "misc_RNA") {
 		    $holder = $misc_RNA_holder;
 		} elsif ($feature_type eq "misc_feature") {
 		    $holder = $misc_feature_holder;
		} elsif ($feature_type eq "CDS") {
		    $holder = $CDS_holder;
		} elsif ($feature_type eq "source") {
		    $holder = $source_holder;
		} else {
		    $holder = $feature_holder;
		}
		$current_feature = $holder->new_object();

		$current_feature->set_attribute("type",$feature_type);
		unless ($feature_type eq "source") {
		    $current_feature->set_attribute("organism",$organism_name);
		    $current_feature->set_attribute("taxid",$taxid);
		}
#		$current_feature->force_attribute("id",$organism_name."_".$id); ### add organism name as prefix
		$current_feature->set_attribute("contig",$current_contig->get_attribute("id"));
		$current_feature->set_attribute("chrom_position",$position);

		#### gene
		if ($feature_type eq 'gene') {
		    $last_gene = $current_feature;

		    #### WARNING : this is VERY tricky : some gene
		    #### names are annotated in the gene, and some in
		    #### the CDS. I automatically add all the names of
		    #### the gene to the CDS. This relies on the
		    #### assumption that each CDS is preceeded by the
		    #### corresponding gene.
		} elsif (($feature_type eq 'CDS') ||
			 ($feature_type eq 'mRNA') ||
			 ($feature_type eq 'tRNA') ||
			 ($feature_type eq 'rRNA') ||
			 ($feature_type eq 'misc_RNA')) {
		    if ($last_gene) {
			#### addd the names of the last gene to the current feature
			foreach my $gene_name ($last_gene->get_attribute("names")) {
			    $current_feature->push_attribute("names", $gene_name);
			}
			my $gene_id =  $last_gene->get_attribute("id");
			$current_feature->set_attribute("gene_id", $gene_id);
#                   print join ("\t", $feature_type, join (";", $last_gene->get_attribute("names")), $gene_id), "\n";
#		    $current_feature->set_attribute("gene_id", $last_gene->get_attribute("id"));

			#### set the type of the last gene to the current feature type, if it has not yet been done
			my $last_gene_type = $last_gene->get_attribute("type");
			if (($last_gene_type eq "") ||
			    ($last_gene_type eq $null) ||
			    ($last_gene_type eq "gene")) {
			    $last_gene->force_attribute("type", $feature_type);
			}
		    }
		}


	    } else {
		warn join "\t", $l, "ignoring ft", $feature_type, $position, "\n" if ($main::verbose >= 3);
		undef($current_feature);
	    }


	    #### new feature attribute
	} elsif ($line =~ /^ {21}\/(\S+)=/) {
	    $attribute_type=$1;
	    $attribute_value=&trim("$'");

	    #### string attributes : make sure that the entire string is parsed
	    if ($attribute_value =~ /^\"/) {
		unless ($attribute_value =~ /\"$/) {
		    do {
			$line = <GBK>;
			chomp $line;
			$attribute_value .= " ";
			$attribute_value .= &trim($line);
		    } until ($line =~ /\"$/);
		}
		$attribute_value =~ s/^\"//;
		$attribute_value =~ s/\"$//;
	    }
	    if ($current_feature) {
		$current_feature->new_attribute_value($attribute_type, $attribute_value);
		warn join "\t", $l, "\tadding attribute", $attribute_type, $attribute_value, "\n" if ($main::verbose >=3);


		#### detect taxid
		if (($current_feature->get_attribute("type") eq "source") ) {
		    if (($attribute_type eq "db_xref") && ($attribute_value =~ /taxon:(\d+)/)){
			$taxid = $1;
			$current_contig->force_attribute("taxid", $taxid);
			$organism->force_attribute("id", $taxid);
			foreach my $class_holder ($contigs,
						  $features, 
						  $genes, 
						  $mRNAs, 
						  $tRNAs, 
						  $rRNAs, 
						  $misc_RNAs, 
						  $misc_features, 
						  $CDSs, 
						  $sources, 
						  $organisms) {
			    $class_holder->set_prefix($taxid);
			}
		    }
		}
		
	    } else {
		warn join "\t", $l, "\tignoring attibute", $attribute_type, $attribute_value, "\n" if ($main::verbose >=5);
	    }
	    

	} else {
	    warn ("file ".$short_name{$org}."\tline $l\tnot parsed\t$line\n") if ($main::verbose >= 2);
	}
    }
    close GBK;

    ################################################################
    ##### parse some names from the notes
    #### guess which "notes" actually
    #### correspond to synonyms.

    #### The basic guessing rule is : any note made of a single word
    #### is considered as a synonym (this is of course not optimal)
    #### For some genomes, additional identifiers/synonyms have
    #### specific formats (locus_tag, Accession, ...).

    warn ("; ",
	  &AlphaDate(),
	  "\tguessing feature synonyms from genbank notes\n")
	if ($verbose >= 1);

    foreach my $feature ($genes->get_objects(),
			 $mRNAs->get_objects(),
			 $CDSs->get_objects()
			 ) {
	
	my @notes = $feature->get_attribute("note");
	foreach my $note (@notes) {
	    $note = &trim($note); ### remove leading and trailing spaces
	    warn "NOTE\t'$note'\n" if ($verbose >= 4);
	    unless ($note =~ /\S/) {
		$feature->push_attribute("names", $note);
		next;
	    }
	    
	    
	    if (
#		($note =~ /^locus_tag\: ([\w_\-\.]+)/i)  ||
		($note =~ /^locus_tag\: ([a-z0-9_\-\.]+)/i)  ||
		($note =~ /^Accession ([\w_\-\.]+)/) ||
		($note =~ /^(\S+)\,\s+len:/)
		){
		
		$new_name = $1;
		
#		#### check whether th name has already been declared
#		my @names = $feature->get_attribute("names");
#		foreach my $name (@names) {
#		    if ($new_names eq $name) {
#			next;
#		    }
#		}
		
		#### add the new name
		warn join( "\t", "Adding name to feature", 
			   $feature->get_attribute("id"), 
			   $new_name
			   ), "\n" if ($verbose >= 4);
		$feature->push_attribute("names", $new_name);
	    }
	}
    }


    return ($file_description, $sequence);
}




1;


