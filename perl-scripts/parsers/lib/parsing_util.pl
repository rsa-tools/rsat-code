#!/usr/bin/perl
############################################################
#
# $Id: parsing_util.pl,v 1.27 2010/11/12 22:09:23 rsat Exp $
#
# Time-stamp: <2003-10-01 17:00:56 jvanheld>
#
############################################################

require "lib/load_classes.pl";
require "lib/util.pl";
require Data::Dumper;

=pod

=head1 NAME

    parsing_util.pl

=head1 DESCRIPTION

    Utilities for parsing various types of flat files. 

=over

=cut


################################################################
=pod

=item 

Check for the existence of a directory, and create it if if does not exist.

=cut
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

################################################################
=pod

=item deliver()

Copy the required files to the delivery directory

=cut

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
=pod

=item ExportClasses()

Export all classes from a list provided as argument.

Usage: &ExportClasses($output_file, $out_format, @classes);

=cut

sub ExportClasses {
  local ($output_file, $out_format, @classes) = @_;
  if ($output_file) {
    open STDOUT, ">$output_file"  || die "Error : cannot write file $output_file\n"; 
  }
  foreach $class (@classes) {
    warn (";\n; ", &RSAT::util::AlphaDate, " exporting class ", $class, " to file '$output_file'\n") 
      if ($verbose >= 1);
    local @selected = @{$out_fields{$class}};
    foreach $object ($class->get_objects()) {
      $object->print_attributes($out_format, @selected);
    }
    warn ("; ", &RSAT::util::AlphaDate, " class ", $class, " exported\n") 
      if ($verbose >= 1);
  }
  close STDOUT if ($output_file);
  return 1;
}




################################################################
=pod

=item ExportMakeFile()

Export a makefile which will call other makefiles for the loading
&ExportMakefile(@classes);

=cut

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
        foreach my $target ("create", "alter", "uncompress", "load", "recompress", "all", "drop") {
	    print MAKEFILE "${target}:\n";
	    my $usage_done = 0;
	    for my $class (@classes) {
		(my $short_class = $class) =~ s/.*:://g;
		my $class_makefile = lc($short_class);
		$class_makefile .= ".mk";
		print MAKEFILE "\tmake -i -f ", $class_makefile, " ", ${target}, "\n";
	    }
	}
	close MAKEFILE;
    }
}


################################################################
=pod

Parse a file in KEGG format (flat file)

=item ParseKeggFile

See KEGG ftp site for format specifications

Usage : $class_holder = &ParseKeggFile($kegg_file,
				       $class_holder,
				       %default_args)

=cut

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
	warn (";\n; ", &RSAT::util::AlphaDate,
	      " parsing class $class from $kegg_file\n")
	    if ($verbose >= 1);
	open INFILE, $kegg_file 
	    || die "Error: cannot $kegg_file\n";
	$in = INFILE;
    } else {
	warn (";\n; ", &RSAT::util::AlphaDate,
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
		    $current_id = $entry_id;
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
=pod

=item ArgString()

Print command-line arguments

=cut

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
=pod

=item PrintStats()

Prints the number of objects and the number of attribute allocations
for each selected class 

Usage : &PrintStats($file,@classes);

=cut

sub PrintStats {
  my ($stat_file, @classes) = @_;
  
  if ($stat_file) {
    open STATS, ">$stat_file" || die "Error: cannot write stat file $stat_file\n";
    $stats = STATS;
  } else {
    $stats = STDOUT;
  }
  
  ### print some stats
  my $alpha_date =  &RSAT::util::AlphaDate();
  &RSAT::message::TimeWarn("Printing parsing statistics in file", $stat_file)
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
	  $count = $attr_counts{$attribute} || 0;
	  $type = $attr_types{$attribute} || "NA";
	  printf $stats ";\t\t%-15s\t%d\t%s\n", $attribute, $count, $type;
	}
      }
    }
  }
  close STATS if ($stat_file);
}

################################################################
=pod

=item SplitReactants()

Split a side of an equation into its reactants, and for each,
separates the stoichiometry from the compound.

=cut

sub SplitReactants {
  my @reactants = split / \+ /, $_[0];
  my %reactants;
  my $stoichio = "";
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
=pod

=item ParsePositions()

Parse position string and add the different fields to each object
(chromosome,  strand,  start_pos, end_pos).

=cut

sub ParsePositions {
    my ($features) = @_;

    &RSAT::message::TimeWarn("Parsing feature positions", $features->get_object_type())
	if ($verbose >= 2);

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
	my $chromosome = $null;
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
	$coord = $chrom_pos;
	if ($chrom_pos =~ /complement\((.*)\)/) {
	    $strand = "R";
	    $coord =~ s/complement\(//g;
	    $coord =~ s/\)//g;
	} else {
	    $strand = "D";
	}

	################################################################
	#### split exons or treat genes accross replication origin (circular genomes)
	if ($coord =~ /^join\(/) { ### exons or genes accross replication origin (circular genomes)
	    $coord =~ s/^join\(//;
	    $coord =~ s/\)$//;

	    #### genes accross replication origin (circular genomes)
	    if ($coord =~ /^([\>\<]{0,1}\d+)\.\.([\>\<]{0,1}\d+),(1)\.\.([\>\<]{0,1}\d+)$/) {
		$start_pos = $1;
		$end_pos = $4;

		warn (join ("\t", "DEBUG", "ParsePositions", $position, , "\n", "chrom_pos", $chrom_pos, "\n", "coord", $strand, $coord, $start_pos, $end_pos, "\n")) if ($main::verbose >= 1);

	    #### exons
	    } else {

	    #### multiple segments
	    my @exons = split ",", $coord;
	    my @exon_starts = ();
	    my @exon_ends = ();

	    #### exon limits
	    foreach my $exon (@exons) {
		$exon =~ s/\s+//g;
		$feature->push_attribute("exons", $exon);
		my ($exon_start, $exon_end) = &segment_limits($exon, $position);
		push @exon_starts, $exon_start;
		push @exon_ends, $exon_end;
	    }
	    
	    #### feature start and end 
	    $start_pos = &min(@exon_starts) || $null;
	    $end_pos = &max(@exon_ends) || $null;

	    #### introns
	    my @introns = ();
	    for my $e (0..$#exon_starts - 1) {
		my $intron = $exon_ends[$e] + 1;
		$intron .= "..";
		$intron .= $exon_starts[$e+1] -1;
		$feature->push_attribute("introns", $intron);
	    }


	}

	} else {
	    #### a single segment
	    ($start_pos, $end_pos) = &segment_limits($coord);
	}

	warn (join ("\t", "DEBUG", "ParsePositions", $position, , "\n", "chrom_pos", $chrom_pos, "\n", "coord", $strand, $coord, "\n")) if ($main::verbose >= 10);


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


# 	&RSAT::message::Debug("ParsePositions", 
# 			      $feature->get_attribute("id"),
# 			      $feature->get_attribute("type"),
# 			      $feature->get_attribute("position"),
# 			      $feature->get_attribute("chromosome"),
# 			      $feature->get_attribute("end_pos"),
# 			      $feature->get_attribute("start_pos"),
# 			      $feature->get_attribute("strand"),
# 			      ) if ($main::verbose >= 0);


    }

}


################################################################
=pod

=item segment_limits()

Return start and end positions of a segment

=cut

sub segment_limits {
    my ($segment, $position) = @_;
    my $segment_start = $null;
    my $segment_end = $null;

    warn "; Segment limits $segment\n" if ($main::verbose >= 10);

    #### start and end positions are different
    if ($segment =~ /^([\>\<]{0,1}\d+)\.\.([\>\<]{0,1}\d+)$/) {
	$segment_start = $1;
	$segment_end = $2;

    #### single residue position
    } elsif ($segment =~ /^([\>\<]{0,1}\d+)$/) {
	$segment_start = $1;
	$segment_end = $1;
    } else {
	&ErrorMessage (join "\t", "Invalid segment format", $segment, $position, "\n");
    }
    return ($segment_start, $segment_end);
}


1;


__END__

=pod

=back
