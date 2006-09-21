###############################################################
#
# Class to handle organisms
#

package RSAT::organism;

use RSAT::GenericObject;
use RSAT::error;
use RSAT::message;
use RSAT::SequenceOnDisk;
use RSAT::GenomeFeature;
use RSAT::Index;
use RSAT::stats;

@ISA = qw( RSAT::GenericObject );

=pod

=head1 NAME

    RSAT::organism

=head1 DESCRIPTION

Object for manipulating organisms.

=cut


################################################################
=pod

=item new()

Create a new organism.

=cut
sub new {
    my ($class, %args) = @_;

    my $self = bless {
	}, $class;
    my $index = new RSAT::Index();
    $self->set_attribute("name_index", $index);

#    die join ("\t", $self, $index, $self->get_attribute("name_index"));

    return $self;
}

################################################################
=pod

=item to_text($format)

Converts the organism in a string for exporting it in the specified
format.

=cut
sub to_text {
    my ($self, $format, $null) = @_;
    my $text = "";

    return $text;
}

return 1;


################################################################
=pod

=item get_contigs()

Usage: my %contigs = $organism->get_contigs();

Return the has table with indexed contigs (key=contig ID, value =
contig object).

=cut
sub get_contigs {
    my ($self) = @_;
    return $self->get_attribute("contigs");
}


################################################################
=pod

=item add_contigs

Add a contig to the list.

=cut
sub add_contig {
    my ($self, $contig_id, $contig_obj) = @_;
    &RSAT::message::Info(join("\t", "Adding contig", $contig_id, $contig_obj)) if ($main::verbose >= 10);
    $self->add_hash_attribute("contigs", $contig_id, $contig_obj);
}


################################################################
=pod

=item is_supported

Indicates whether a given organism name is supported on this RSAT
site.  Return value is boolean (1=true, 0=false).

=cut
sub is_supported {
    my ($self, $organism_name) = @_;
    if (defined($main::supported_organism{$organism_name})) {
	return 1;
    } else {
	return 0;
    }
}

################################################################
=pod

=item

Check if the specified organism has been installed on this RSAT
site. If not, die.

=cut
sub check_name {
    my ($self, $organism_name) = @_;
    unless ($organism_name) {
	$organism_name = $self->get_attribute("name");
    }
    unless ($organism_name) {
	&RSAT::error::FatalError("you should specify an organism name");
    }
    my $supported = $self->is_supported($organism_name);
    unless ($supported) {
	 &RSAT::error::FatalError("Organism $organism_name is not supported.",
		     "Use the command supported-organisms for a list of supported organisms");

     }
}



################################################################
=pod

=item OpenContigs

Check if an organism is supported on the current installation,
and open streams to read its genome sequence.

Usage:
   Automatic selection of genome and feature file :
       $organism->OpenContigs($organism_name);

   Manual specification of input files :
       $organism->OpenContigs($organism_name,
                              $annotation_table,
                              $input_sequence_file,
                              $input_sequence_format);

=cut
sub OpenContigs {
    my ($self,
	$organism_name,
	$annotation_table,
	$input_sequence_file,
	$input_sequence_format,
	%args) = @_;

    my $genome_file;
    my $genome_file_format;

    if ($main::verbose >= 2) {
	&RSAT::message::Info(join ("\t", "Opening contigs", @_));
	&RSAT::message::Info(join ("\t",
		       "Manually specified input sequence file",
		       $input_sequence_file,
		       $input_sequence_format))
	    if ($input_sequence_file);

    }

    ################################################################
    #### check input sequence file and format
    if ($input_sequence_file) {
	$genome_file = $input_sequence_file;
	$genome_seq_format = $input_sequence_format;
    } elsif ($organism_name) {
	$self->check_name($organism_name);
	$genome_file = $main::supported_organism{$organism_name}->{'genome'};
	$genome_seq_format = $main::supported_organism{$organism_name}->{'seq_format'};
    } else {
	&RSAT::error::FatalError("You should specify an organism.", "Use the command supported-organisms to get a list of suported organisms");
    }
    if ($main::verbose >= 4) {
	&RSAT::message::Info("genome file\t".$genome_file);
	&RSAT::message::Info("genome format\t".$genome_seq_format) ;
    }

    #### check feature file
    unless ($annotation_table) {
	$annotation_table = $main::supported_organism{$organism_name}->{'features'};
    }
    if ($main::verbose >= 4) {
	&RSAT::message::Info("feature table\t".$annotation_table);
    }

    ################################################################
    ### open a direct read access to the contig files
    if ($genome_seq_format eq "raw") {
	$contig_id = $organism_name;
	$genome_file =~ s/\.fna$/\.raw/;
	$contig_seq{$contig_id} = new RSAT::SequenceOnDisk (filename=>  $genome_file,
							 id=>        $contig_id,
							 circular=>  1, ### for bacterial genomes
							 organism=>  $organism_name);
	$contig_obj = new RSAT::contig();
	$contig_obj->push_attribute("names", $contig_id);
	$contig_obj->set_sequence_object($contig_seq{$contig_id});
	$contig_obj->set_organism($organism_name);
	$self->add_contig($contig_id, $contig_obj);

    } elsif ($genome_seq_format eq "filelist") {

	## Check the existence of genome file
	unless (-e $genome_file) {
	    &RSAT::error::FatalError("Genome file $genome_file does not exist.");
	}

	## Separate the genome directory and file name
	my ($seq_dir, $seq_list) = &RSAT::util::SplitFileName($genome_file);

	if ($main::verbose >= 4) {
	    &RSAT::message::Info(join ("\t", "Genome dir", $seq_dir));
	    &RSAT::message::Info(join ("\t", "Sequence list", $seq_list));
	}

	## Read the list of contigs for this genome, and open a stream for each contig
	open FILES, $genome_file ||
	    &RSAT::error::FatalError("Cannot open genome file $genome_file for reading.");
	while (<FILES>) {
	    chomp;
	    my ($seq_file, $contig_id, $type) = split "\t";

	    ## Contig ID
	    $contig_id = $seq_file unless ($contig_id);

	    ## Specify whether the contig is circular
	    my $circular = 0;
	    if ((defined($type)) && ($type eq "circular")) {
		$circular = 1;
	    }

	    ## Repeat masked
	    if ($args{rm}) {
		my $masked_seq_file = $seq_file;
		$masked_seq_file =~ s/\.raw/_repeat_masked.raw/;
		if (-e $seq_dir.$masked_seq_file) {
		    $seq_file = $masked_seq_file;

		    &RSAT::message::Warning(join("\t",
						 "Using masked repeat version of contig",
						 $contig_id,
						 $seq_dir.$seq_file,
						 )) if ($main::verbose >= 2);
		} else {
		    &RSAT::message::Warning(join("\t",
						 "There is no masked repeat version of contig",
						 $contig_id,
						 "missing file", $seq_dir.$masked_seq_file,
						 ));
		}
	    }

	    &RSAT::message::Info(join ("\t", "Contig sequence file", $seq_file, "Contig", $contig_id, "circular: $circular") )
		if ($main::verbose >= 2);

	    $contig_seq{$contig_id} = new RSAT::SequenceOnDisk(filename=>  $seq_dir.$seq_file,
							       id=>        $contig_id,
							       circular=>  $circular,
							       organism=>  $organism_name);
	    $contig_obj = new RSAT::contig();
	    $contig_obj->push_attribute("names",$contig_id);
	    $contig_obj->set_sequence_object($contig_seq{$contig_id});
	    $contig_obj->set_organism($organism_name);
	    $self->add_contig($contig_id, $contig_obj);

 	    warn join ("\t","CONTIG",
 		       $seq_file,
 		       $contig_id,
 		       $contig_obj,
 		       $contig_obj->get_attribute("sequence"),
 		       $contig_obj->get_attribute("sequence")->get_length()
		      ), "\n"  if ($main::verbose >= 10);

 	    warn join ("\t", "SEQUENCE",
 		       $contig_id,
 		       $contig_seq{$contig_id}->get_attribute("filename"),
 		       $contig_seq{$contig_id}->get_length(),
 		       $contig_seq{$contig_id}->get_attribute("circular"),
 		       $contig_seq{$contig_id}->get_sequence(1,12,"D")
		      ), "\n" if ($main::verbose >= 10);

	}
	close FILES;

    } elsif ($accepted_input_seq{$genome_seq_format}) {
	unless (-e $genome_file) {
	    &RSAT::error::FatalError("Genome file $genome_file does not exist.");
	}
	my ($in, $input_dir) = &OpenInputFile($genome_file, $genome_seq_format);

	while ((($current_seq, $current_id, @comments) = &ReadNextSequence($in, $genome_seq_format, $input_dir)) &&
	       (($current_seq) || ($current_id))) {
	    $contig_id = $current_id;
	    $contig_seq{$contig_id} = new RSAT::Sequence (id=>        $contig_id,
						       sequence=>  $current_seq,
						       organism=>  $organism_name);
	    $contig_obj = new RSAT::contig();
	    $contig_obj->push_attribute("names",$contig_id);
	    $contig_obj->set_sequence_object($contig_seq{$contig_id});
	    $contig_obj->set_organism($organism_name);
	    $self->add_contig($contig_id, $contig_obj);
	}
    } else {
	&RSAT::error::FatalError("Sequence format $genome_seq_format is not supported for genomes.");
    }
}


################################################################
=pod

=item index_attribute_by_feature

Return a hash table with an index of a specified feature attribute.

Usage: my %attr_index = $organism->index_attribute_by_feature($attribute);

Example: my %left  = $organism->index_attribute_by_feature("left");
Returns a hashtabe with left positions as keys, and RSAT::GenomeFeature as object

=cut
sub index_attribute_by_feature {
    my ($self, $attribute) = @_;
    my %index = ();
    my $null = "<NULL>";
    foreach my $feature ($self->get_attribute("features")) {
	my $value  = $feature->get_attribute($attribute);
	if (($value) && ($value ne $null)) {
	    $index{$feature} = $value;
	} else {
	    &RSAT::message::Warning(join "\t",
				    "Attribute", $attribute,
				    "is not defined for feature",
				    $feature->get_attribute("id")) if ($main::verbose >= 3);
	}
    }

#    &RSAT::message::Debug("Indexed values", $attribute, scalar(keys %index)) if ($main::verbose >= 10);

    return %index;
}

################################################################
=pod

=item DefineAcceptedFeatureTypes(@accepted_feature_types)

Define the accepted feature types.
If the parameter @accepted_feature_table is empty, the default feature type is set to cds.

=cut
sub DefineAcceptedFeatureTypes {
  my ($self, @accepted_feature_types) = @_;

  unless (scalar(@accepted_feature_types) > 0) {
    @accepted_feature_types = 'cds';
  }


  foreach my $feature_type (@accepted_feature_types) {
    &RSAT::message::Info(join("\t", "Adding feature type", $feature_type)) if ($main::verbose >= 3);
    $self->push_attribute("feature_types", $feature_type);
  }
}


################################################################
=pod

=item LoadFeatures()

Read annotations for the current organism.

Usage: &LoadFeatures($annotation_table, $feature_types)

@param $annotation_table  name of the file containing the annotations
@param $feature_types
@param $imp_pos accept to load imperfectly specified positions

The features are embedded in the organism object, and indexed by contig.

=cut
sub LoadFeatures {
  my ($self,
      $annotation_table,
#      $feature_types,
      $imp_pos)= @_;

  ## Organism name
  my $organism_name = $self->get_attribute("name");
  my %contig = $self->get_contigs();
  my $name_index = $self->get_attribute("name_index");

  ## Parse feature types
#  my @feature_types = ();
#  my %accepted_feature_types = ();
#  if ($feature_types) {
#    @feature_types = split ",", $feature_types;
#  } else {
#    @feature_types = qw (cds);
#  }
  my @feature_types = $self->get_attribute("feature_types");
  if (scalar(@feature_types) < 1) {
    @feature_types = qw (cds);
  }

  foreach my $feature_type (@feature_types) {
    warn join ("\t", ";", "accepting feature type", $feature_type), "\n" if ($main::verbose >= 10);
    $accepted_feature_types{$feature_type}++;
  }
  warn join ("\t", "; Accepted feature types",
	     join( ",", keys %accepted_feature_types)), "\n"
	       if ($main::verbose >= 2);

  ## Annotation table
  if ($annotation_table) {
    $self->push_attribute("annotation_tables", $annotation_table);
  } else {
    #	$annotation_table = $main::supported_organism{$organism_name}->{'features'};
    foreach my $type (@feature_types) {
      $annotation_table = join("", $main::supported_organism{$organism_name}->{'data'}, "/genome/", $type, ".tab");
      $self->push_attribute("annotation_tables", $annotation_table);
    }
  }

  foreach my $annotation_table ($self->get_attribute("annotation_tables")) {

    &RSAT::message::Info(join ("\t",
			       &RSAT::util::AlphaDate(),
			       "Loading annotation table",
			       $self->get_attribute("name"),
			       $annotation_table)
			) if ($main::verbose >= 2);

    ## Default column order for genomic features
    ## Note that this order can be redefined by the header of the
    ## annotation table (see below)
    my %col = ();
#     $col{'id'} = 0;
#     $col{'type'} = 1;
#     $col{'name'} = 2;
#     $col{'ctg'} = 3;
#     $col{'left'} = 4;
#     $col{'right'} = 5;
#     $col{'strand'} = 6;
#     $col{'descr'} = 7;
#     $col{'location'} = 8;

    ### Read feature positions (genome annotations)
    my ($annot, $annot_dir) = &main::OpenInputFile($annotation_table);
    my %type = ();
    my $linenb = 0;
    while (my $line = <$annot>) {
      $linenb++;

      if (($main::verbose >= 2) && ($linenb % 1000 == 1)) {
	&RSAT::message::psWarn("Loaded features", $linenb);
      }

      chomp($line);
      next unless ($line =~ /\S/);
      if ($line =~ /^\-\-/) {
	## Internal colunm specification in tables resulting from RSAT parsers
	if ($line =~ /^-- field (\d+)\t(\S+)/) {
	  my $field_column = $1;
	  my $field = lc($2);
	  ## Convert field names
	  $field =~ s/contig/ctg/;
	  $field =~ s/start_pos/left/;
	  $field =~ s/end_pos/right/;
	  $field =~ s/description/descr/;
	  $field =~ s/chrom_position/location/;
	  $col{$field} = $field_column - 1;
	  &RSAT::message::Info(join("\t", "Column specification", $field_column, $field)) if ($main::verbose >= 3);
	}
	next;
      }
      next if (($line =~ /^;/) || ($line =~ /^\#/));

      ## Split the columns of the input row
      my @fields = split "\t", $line;
      foreach $f (keys %col) {
	$$f = $fields[$col{$f}];
#	&RSAT::message::Debug("field",$f, $$f, "column", $col{$f}) if ($main::verbose >= 10);
      }

      ## Check if the type of this feature is accepted
      unless ($accepted_feature_types{lc($type)}) {
	&RSAT::message::Info(join( "\t","skipping feature", $id, "Non-selected feature type",$type))
	  if ($main::verbose >= 3);
	next;
      }

      ################################################################
      #### check mandatory attributes

      ## Check ID
      unless ($id) {
	&RSAT::message::Warning("invalid orf identifier specification in the feature table line $linenb\n;\t",join "\t", @fields) if ($main::verbose >= 3);
	next;
      }

      ## Check contig
      unless ($ctg) {
	&RSAT::message::Warning("invalid contig specification in the feature table line $linenb\n;\t",join "\t", @fields) if ($main::verbose >= 3);
	next;
      }

      ## Check left
      unless ($left) {
	&RSAT::message::Warning("left position not specified in the feature table line $linenb\n;\t",join "\t", @fields) if ($main::verbose >= 3);
	next;
      }
      unless (&RSAT::util::IsNatural($left) ) {
	if ($imp_pos) {
	  &RSAT::message::Warning("imprecise specification of the left position for gene $id\n;\t",join "\t", @fields) if ($main::verbose >= 2);
	  $left =~ s/\>//;
	  $left =~ s/\<//;
	} else {
	  &RSAT::message::Warning("invalid left position specification in the feature table line $linenb\n;\t",join "\t", @fields) if ($main::verbose >= 3);
	  next;
	}
      }

      ## Check right
      unless ($right) {
	&RSAT::message::Warning("Right position not specified in the feature table line $linenb\n;\t",join "\t", @fields) if ($main::verbose >= 3);
	next;
      }
      unless (&RSAT::util::IsNatural($right) ) {
	if ($imp_pos) {
	  &RSAT::message::Warning("imprecise specification of the right position for gene $id\n;\t",join "\t", @fields) if ($main::verbose >= 2);
	  $right =~ s/\>//;
	  $right =~ s/\<//;
	} else {
	  &RSAT::message::Warning("invalid right position specification in the feature table line $linenb\n;\t",join "\t", @fields) if ($main::verbose >= 3);
	  next;
	}
      }

      unless ($left < $right) {
	&RSAT::message::Warning("left should be smaller than right position specification in in  feature table line $linenb\n;\t",join "\t", @fields) if ($main::verbose >= 3);
	next;
      }

      ## Check strand
      unless ($strand) {
	&RSAT::message::Warning("invalid strand specification in the feature table line $linenb\n;\t",join "\t", @fields) if ($main::verbose >= 3);
	next;
      }

      #### make sure the strand format is correct
      $strand = &RSAT::util::ConvertStrand($strand);

      ### Create a new genomic feature
      my $feature = new RSAT::GenomeFeature();
      foreach $f (keys %col) {
	$feature->set_attribute($f, $$f);
      }
      $self->push_attribute("features", $feature);

      ## Accept ID as name
      $feature->push_attribute("names", $id);
      $name = $id unless defined($name);
      if (($name eq "<NULL>") || ($name eq "")) {
	$name = $id;
	$feature->force_attribute("name", $id);
      }
      if ($name ne $id) {
	$feature->push_attribute("names", $name);
      }
      $type{$id} = $type;	## For the loading statistics

      ################################################################
      ## Index genome features by names and ID
      ## The key is the name, the value is the RSAT::GenomeFeature object
      ## The key is converted to uppercases to allow case-insensitive searches
      my $null = "<NULL>";
      $name_index->add_value(uc($name), $feature) if (($name) && ($name ne $null)) ;
      $name_index->add_value(uc($id), $feature) if (($id) && ($id ne $null)) ;


      ## Add the new feature to the contig
      unless (defined($contig{$ctg})) {
	&RSAT::message::Info(join("\t", "Creating new contig", $ctg)) if ($main::verbose >= 3);
	$contig{$ctg} = new RSAT::contig(id=>$ctg);
	$contig{$ctg}->set_organism($organism_name);
      }
      $contig{$ctg}->add_gene($feature);

#      &RSAT::message::Debug("Loaded feature",
#			    $contig{$ctg}->count_genes(),
#			    $id,$ctg,$strand,$left,$right,
#			    $feature->get_attribute("id"),
#			    $feature->get_attribute("ctg"),
#			    $feature->get_attribute("strand"),
#			    $feature->get_attribute("left"),
#			    $feature->get_attribute("right"),
#			    $feature->get_attribute("descr"),
#			    $descr{$id},
#			    ) if ($main::verbose >= 10);
    }
    close $annot if ($annotation_table);
  }

  ## Check the number of features
  if (scalar($self->get_attribute("features")) < 1) {
    &RSAT::message::Warning("There is no annotated feature of type ".$feature_types." in the genome of ".$organism_name);
  } elsif ($main::verbose >= 2) {
    ## Print stats on the features
    &RSAT::message::Info(join ("\t",
			       "Loaded",
			       scalar($self->get_attribute("features")),
			       "features for organism",
			       $self->get_attribute("name"),
			      ));
    my %stats = ();
    foreach my $id (keys %type) {
      $stats{$type{$id}}++;
    }
    foreach my $type (sort keys %stats) {
      &RSAT::message::Info(join ("\t", "", $stats{$type}, $type));
    }
  }
  #    my %test = $self->get_attribute("feature_id");
  #    &RSAT::message::Debug("Feature keys", scalar(keys(%test))) if ($main::verbose >= 10);
}


################################################################
=pod

=item CalcNeighbourLimits()

Calculate the limits of the non-coding region on both sides of each
gene

=cut

sub CalcNeighbourLimits {
  my ($self) = @_;

  my %contig = $self->get_contigs();

  &RSAT::message::TimeWarn("Calculating neighbour limits") if ($main::verbose >= 1);

  my %left = $self->index_attribute_by_feature("left");
  my %right = $self->index_attribute_by_feature("right");
  foreach my $ctg (sort keys %contig) {
    my @genes = sort { $left{$a} <=> $left{$b} } $contig{$ctg}->get_genes();
    my $contig_length = $contig{$ctg}->get_length();

    &RSAT::message::Info(join ("\t", "Features per contig", $ctg, scalar(@genes))) if ($main::verbose >= 2);
    @ctg_lefts = sort {$a <=> $b} @left{@genes};
    @ctg_rights = sort {$a <=> $b} @right{@genes};


    ################################################################
    ## Identify the left neighbour of each feature
    ## TO CHECK: does it correctly treat the case of completely overlapping genes, embedded genes
    for my $g (0..$#genes) {
      my $gene = $genes[$g];
      my $ln = $g -1;		### first guess for left neighbour
      my $found = 0;

      do {
	&RSAT::message::Debug("contig", $ctg, "gene",
			      "g=".$g,
			      "geneId=".$gene->get_attribute("geneid"),
			      "name=".$gene->get_attribute("name"),
			      "candidate left neighbour",
			      "ln=".$ln,
			      "GeneID:".$genes[$ln]->get_attribute("geneid"),
			      "name=".$genes[$ln]->get_attribute("name"),
			     ) if ($main::verbose >= 4);

	if ($ln < 0) {
	  ## Leftmost gene of a contig

	  &RSAT::message::Debug("No left neighbour for genomic feature",
				"g=".$g,
				"ID=".$genes[$g]->get_attribute("id"),
				"GeneID=".$genes[$g]->get_attribute("geneid"),
				"name=".$genes[$g]->get_attribute("name"),
			       ) if ($main::verbose >= 10);

	  $gene->set_attribute("left_neighbour", "<NULL>");
	  $gene->set_attribute("left_limit", 1);
	  $gene->set_attribute("left_size", $ctg_lefts[$g]);
	  next;

	} elsif ($gene->get_attribute("geneid") eq $genes[$ln]->get_attribute("geneid")) {

	  ## Skip features having the same GeneID as the current
	  ## feature (e.g. multiple mRNA for a same gene, due to
	  ## the presence of alternative transcription start
	  ## sites (TSS) or to alternative splicing.

	  &RSAT::message::Debug("genomic feature",
				"ln=".$ln,
				"ID=".$genes[$ln]->get_attribute("id"),
				"GeneID=".$genes[$ln]->get_attribute("geneid"),
				"name=".$genes[$ln]->get_attribute("name"),
				"has same GeneID as genomic feature",
				"g=".$g,
				"ID=".$gene->get_attribute("id"),
				"GeneID=".$gene->get_attribute("geneid"),
				"name=".$gene->get_attribute("name"),
			       ) if ($main::verbose >= 3);

	  $ln -= 1;
	  next;

	} elsif (($ctg_rights[$ln] > $left{$gene}) && ($ctg_lefts[$ln] > $left{$gene})) {
	  ## candidate left neighour gene is completely embedded in current gene -> select the next left candidate
	  $ln--;
	  &RSAT::message::Debug("genomic feature",
				$ln,
				$genes[$ln]->get_attribute("geneid"),
				$genes[$ln]->get_attribute("name"),
				"is embedded in genomic feature", $g,
				$gene->get_attribute("geneid"),
				$gene->get_attribute("name"),
			       ) if ($main::verbose >= 4);

	} elsif ($ctg_rights[$ln+1] < $left{$gene}) {
	  $ln++;

	} else {
	  $found = 1;
	}
      } until (($found) || ($ln < 0) || ($ln > $#ctg_rights));

      if ($found) {
	$neighb_left_limit = $ctg_rights[$ln];
      } elsif ($ln < 0) {
	$neighb_left_limit = 0;
      } else {
	$neighb_left_limit = undef;
      }
      $neighb_left_size = &RSAT::stats::max(0, $left{$gene} - $neighb_left_limit -1);

#      &RSAT::message::Debug("Identified left neighbour", $ln,
# 			    $genes[$ln]->get_attribute("geneid"),
# 			    $genes[$ln]->get_attribute("name"),
# 			    "for gene", "g=".$g,
# 			    $gene->get_attribute("geneid"),
# 			    $gene->get_attribute("name"),
# 			   ) if ($main::verbose >= 10);

      $gene->set_attribute("left_neighbour", $genes[$ln]);
      $gene->set_attribute("left_limit", $neighb_left_limit);
      $gene->set_attribute("left_size", $neighb_left_size);
    }


    ################################################################
    ## Identify the right neighbour for each feature
    for my $g (0..$#genes) {
      my $gene = $genes[$g];
      my $rn = $g +1;		### first guess for right neighbour
      my $found = 0;

      #### Calculate intergenic limit on the right side of the gene
      $found = 0;
      do {
	if ($rn > $#genes) {
	  ## Rightmost gene of a contig
# 	  &RSAT::message::Debug("No right neighbour for genomic feature",
# 				"g=".$g,
# 				"ID=".$genes[$g]->get_attribute("id"),
# 				"GeneID=".$genes[$g]->get_attribute("geneid"),
# 				"name=".$genes[$g]->get_attribute("name"),
# 			       ) if ($main::verbose >= 10);
	  $gene->set_attribute("right_neighbour", "<NULL>");
	  $gene->set_attribute("right_limit", $contig_length);
	  $gene->set_attribute("right_size", $contig_lenth - $ctg_rights[$g]);
	  next;
	} elsif ($gene->get_attribute("geneid") eq $genes[$rn]->get_attribute("geneid")) {

	  ## Skip features having the same GeneID as the current
	  ## feature (e.g. multiple mRNA for a same gene, due to
	  ## the presence of alternative transcription start
	  ## sites (TSS) or to alternative splicing.

	  &RSAT::message::Debug("genomic feature",
				"rn=".$rn,
				"ID=".$genes[$rn]->get_attribute("id"),
				"GeneID=".$genes[$rn]->get_attribute("geneid"),
				"name=".$genes[$rn]->get_attribute("name"),
				"has same GeneID as genomic feature",
				"g=".$g,
				"ID=".$gene->get_attribute("id"),
				"GeneID=".$gene->get_attribute("geneid"),
				"name=".$gene->get_attribute("name"),
			       ) if ($main::verbose >= 10);

	  $rn += 1;
	  next;

	} elsif (($ctg_lefts[$rn] < $right{$gene}) && ($ctg_rights[$rn] < $right{$gene})) {
	  $rn++;
	} elsif ($ctg_lefts[$rn-1] > $right{$gene}) {
	  $rn--;
	} else {
	  $found = 1;
	}
      } until (($found) || ($rn < 0) || ($rn > $#ctg_lefts));
      if ($found) {
	$neighb_right_limit = $ctg_lefts[$rn];
      } elsif ($rn < 0) {
	$neighb_right_limit = undef;
      } else {
	$neighb_right_limit = $contig_seq{$ctg}->get_length() + 1;
      }
      $neighb_right_size = &RSAT::stats::max(0, $neighb_right_limit - $right{$gene} -1);

      $gene->set_attribute("right_neighbour", $genes[$rn]);
      $gene->set_attribute("right_limit", $neighb_right_limit);
      $gene->set_attribute("right_size", $neighb_right_size);
    }

    ################################################################
    ## Calculate upstream and downstream limits (according to the strand)
    for my $g (0..$#genes) {
      my $gene = $genes[$g];
      my $un;			## upstream neighbour index
      my $dn;			## upstream neighbour index

      if ($gene->get_attribute("strand") eq "R") {

	## Upstream neighbour is on the right side
	$gene->set_attribute("upstr_neighbour", $gene->get_attribute("right_neighbour"));
	$gene->set_attribute("upstr_limit", $gene->get_attribute("right_limit"));
	$gene->set_attribute("upstr_size", $gene->get_attribute("right_size"));

	## Downstream neighbour is on the left side
	$gene->set_attribute("downstr_neighbour", $gene->get_attribute("left_neighbour"));
	$gene->set_attribute("downstr_limit", $gene->get_attribute("left_limit"));
	$gene->set_attribute("downstr_size", $gene->get_attribute("left_size"));

      } else {
	## Upstream neighbour is on the left side
	$gene->set_attribute("upstr_neighbour", $gene->get_attribute("left_neighbour"));
	$gene->set_attribute("upstr_limit", $gene->get_attribute("left_limit"));
	$gene->set_attribute("upstr_size", $gene->get_attribute("left_size"));

	## Downstream neighbour is on the right side
	$gene->set_attribute("downstr_neighbour", $gene->get_attribute("right_neighbour"));
	$gene->set_attribute("downstr_limit", $gene->get_attribute("right_limit"));
	$gene->set_attribute("downstr_size", $gene->get_attribute("right_size"));
      }

#      &RSAT::message::Debug ("neighbours",
#			     $g, $ctg, $gene->get_attribute("id"), $left{$gene}, $right{$gene}, $gene->get_attribute("strand"),
#			     "left", "gene=".$gene->get_attribute("left_neighbour"), "limit=".$gene->get_attribute("left_limit"), "limit=".$gene->get_attribute("left_size"),
#			     "right", "gene=".$gene->get_attribute("right_neighbour"), "limit=".$gene->get_attribute("right_limit"), "limit=".$gene->get_attribute("right_size"),
#			     "upstr", "gene=".$gene->get_attribute("upstr_neighbour"), "limit=".$gene->get_attribute("upstr_limit"), "limit=".$gene->get_attribute("upstr_size"),
#			     "downstr", "gene=".$gene->get_attribute("downstr_neighbour"), "limit=".$gene->get_attribute("downstr_limit"), "limit=".$gene->get_attribute("downstr_size"),
#			    ) if ($main::verbose >= 10);

    }
  }
  &RSAT::message::TimeWarn("Calculated neighbour limits") if ($main::verbose >= 2);
  &RSAT::message::psWarn("Calculated neighbour limits") if ($main::verbose >= 2);
}


################################################################
=pod

=item LoadSynonyms

Read a list of synonyms for gene names, given an organism name.

Usage : $organism->LoadSynonyms();

=cut
sub LoadSynonyms {
  my ($self) = @_;
  my $organism_name = $self->get_attribute("name");
  my $name_index = $self->get_attribute("name_index");


  ## Annotation table
  my @feature_types = $self->get_attribute("feature_types");
  if (scalar(@feature_types) < 1) {
    @feature_types = qw (cds);
  }
  foreach my $type (@feature_types) {
    $synonym_file = join("", $main::supported_organism{$organism_name}->{'data'}, "/genome/", $type, "_names.tab");
    if (-e $synonym_file) {
	$self->push_attribute("synonym_files", $synonym_file);
    } else {
	&RSAT::message::Warning(join("\t", "synonym file does not exist", $synonym_file));
    }
  }

  foreach my $synonym_file ($self->get_attribute("synonym_files")) {
      &RSAT::message::Info(join("\t","Loading synonyms from file", $synonym_file)) if ($main::verbose >= 2);

    #    my $synonym_file = $main::supported_organism{$organism_name}->{'synonyms'};
#    unless ($synonym_file) {
      #	&RSAT::message::Warning(join "\t", "no synonym table for organism",$organism_name)
      #	  if ($main::verbose >= 1);
      #	return;
      #      }

    my ($syn) = &RSAT::util::OpenInputFile( $synonym_file);
    while (<$syn>) {
      next if ((/^;/) || (/^\#/) || (/^--/)); ## skip comment lines
      next unless (/\S/);			## skip empty lines
      chomp();
      my @fields = split "\t";
      my $id = uc($fields[0]);
      my $name = $fields[1];
      my $feature = $self->get_feature_for_name($id);
      if ($feature) {
	$feature->push_attribute("names", $name);
	$name_index->add_value(uc($name), $feature);
	#	    &RSAT::message::Debug(join ("\t", "Added synonym", $id, $name)) if ($main::verbose >= 10);
	} else {
	  &RSAT::message::Warning(join ("\t", "Cannot add synonym", "no feature with ID", $id)) if ($main::verbose >= 4);
	}
    }
    close $syn if ($synonym_file);
  }
}


################################################################
=pod

=item SelectRandomORFs

return a random ORF selection

Usage: @random_orfs = $organism->SelectRandomORFs($n, $replace, $init, $seed);

@param $replace  if true, the selections is made with replacement
@param $init     if true, the random generator is initialized (srand)
@param $seed     specify the seed for the random generator

=cut
sub SelectRandomGenes {
  my ($self, $rand_gene_nb, $replace, $init, $seed) = @_;
  my @random_genes = ();
  my @random_ids = ();

  #### initalize the random generator
  if ($init) {
      if ($seed) {
	  ## User-defined seed
	  srand($seed);
      } else {
	  ## Use current time as random seed
	  srand (time);
      }
  }

  my @genes = $self->get_attribute("features");
  for my $i (1..$rand_gene_nb) {
      my $remaining_genes = scalar(@genes);
      my $selected = int(rand($remaining_genes));
      warn ";", join ("\t", "Selected gene", $selected, $remaining_genes), "\n" if ($main::verbose >= 2);
      if ($replace) {
	  push @random_genes, $genes[$selected];
      } else {
	  push @random_genes, splice(@genes, $selected, 1);
      }
  }

  foreach my $g (@random_genes) {
      push @random_ids, $g->get_attribute("id"),
  }

  return @random_ids;
}


################################################################
=pod

=item get_features_for_name

Return genomic features (RSAT::GenomicFeature objects) given an ID or
a name. Due to homonyms, a single name can return several features.

Usage: my @ft = $organism->get_features_for_name($query);

=cut
sub get_features_for_name {
    my ($self, $query) = @_;
    my $name_index = $self->get_attribute("name_index");
    my @values = $name_index->get_values(uc($query));
#    &RSAT::message::Debug("getting features for name", $query, scalar(@values), join(";", @values));
    return  @values;
}

################################################################
=pod

=item get_feature_for_name

Return a genomic feature (RSAT::GenomicFeature objects) given an ID or
a name. In case of homonyms, only the first value is returned.

Usage: my $ft = $organism->get_feature_for_name($query);

=cut
sub get_feature_for_name {
    my ($self, $query) = @_;
    my @result = $self->get_features_for_name($query);
    return $result[0];
}

return 1;


__END__


