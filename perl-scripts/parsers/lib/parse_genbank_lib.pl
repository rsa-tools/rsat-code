#!/usr/bin/perl
############################################################
#
# $Id: parse_genbank_lib.pl,v 1.46 2012/10/03 05:21:24 jvanheld Exp $
#
# Time-stamp: <2003-10-01 17:00:56 jvanheld>
#
############################################################


=pod

=head1 NAME

    parse_genbank_lib.pl

=head1 DESCRIPTION

    Utilities for parsing NCBI/Genbank flat files.

=over

=cut


################################################################
=pod

=item ReadNextLine()

Read one line of the Genbank file, increment line counter, and chomp.

=cut

sub ReadNextLine {
    my $line = <GBK>;
    print STDERR "$l\t$line" if ($main::verbose >= 10);
    $l++;
    $line =~ s/\r//;
    chomp($line);
    return ($line);
}



################################################################
=pod

=item ParseAllGenbankFiles

Parse a list of genbank files (each file normally contains
information about one chromosome).

This procedure calls &ParseGenbankFile on each input file, and then
collects the names for each feature type.

=cut

sub ParseAllGenbankFiles {
    my (@genbank_files) = @_;
    foreach my $file (@genbank_files) {
	$file{input} = $file;

	#### check whether the file has to be uncompressed on the fly
	if ($file{input} =~ /.gz/) {
	    $file{input} = "gunzip -c $file{input} |";
	} else {
	    $file{input} = "cat $file{input} |";
	}

	#### for quick testing; only parse the first  lines
	if ($main::test) {
	    $file{input} .= " head -$test_lines | ";
	}

	&ParseGenbankFile($file{input},
			  $features,
			  $genes,
			  $mRNAs,
			  $scRNAs,
			  $tRNAs,
			  $rRNAs,
			  $repeat_regions,
			  $misc_RNAs,
			  $misc_features,
			  $CDSs,
			  $contigs,
			  $organisms,
			  $sources,
			  source=>"Genbank",
			  seq_dir=>$main::dir{sequences}
			 );
    }

#    ## For Genes, use the GeneID as unique identifier
#    foreach my $object ($genes->get_objects()) {
#	$object->UseGeneIDasID();
#    }


    ## DEBUG NOTE: ParseFeatureNames and CheckObjectNames are
    ## apparently partly redundant. I should check and remove
    ## redundancy.

    ## Parse feature names
    &ParseFeatureNames($genes, $mRNAs, $scRNAs, $tRNAs, $rRNAs, $CDSs);

    ## Extract cross-references
    &ExtractCrossReferences($genes, $mRNAs, $scRNAs, $tRNAs, $rRNAs, $CDSs);

    ## Update name index for genes (since the original IDs have been changed for GeneIDs)
    $genes->index_ids();
    $genes->index_names();

    ################################################################
    ## Replace IDs for the different feature types.
    ##
    ## This is still a bit tricky (2006/01/06), because in NCBI/Genbank there
    ## is no suitable ID which would be defined for all the feature types.
    ## For example
    ## - the locus_tag is not unique : when a gene contains several mRNA and
    ##   CDS (alternative splicing), they all bear the same locus_tag (see
    ##   Arabidopsis for an example).
    ## - transcript_id and GI are unique identifiers for mRNA
    ## - protein_id and GI are unique identifiers for CDS
    ## - tRNA has no GI !
    ##
    ## Temporarily, we use
    ## - GeneID for genes
    ## - protein_id for CDS
    ## - transcript_id for mRNA
    ## - locus_tag for all the other feature types
    ##
    ## In addition, each genome has its specific problems.
    ## For instance, in Debaryomyces_hansenii_CBS767, some CDSs have
    ## no protein_id, whereas some other CDS do.

    ## For CDS, use the GI or the protein_id the locus_tag as unique identifier
    my @preferred_ids_cds;
    if (defined($main::preferred_id{cds})) {
      @preferred_ids_cds = $main::preferred_id{cds};
    }
    push @preferred_ids_cds, qw (locus_tag protein_id GI GeneID gene); ## These IDs will be used if the preferred ID is not found

    &RSAT::message::TimeWarn("Finding preferred IDs for CDS", join(",", @preferred_ids_cds)) if ($main::verbose >= 3);
    foreach my $object ($CDSs->get_objects()) {
      $object->UseAttributesAsID(@preferred_ids_cds);
    }
    $CDSs->index_ids();

    ## For mRNA, use the GI or the transcript_id the locus_tag as unique identifier
    my @preferred_ids_mrna;
    if (defined($preferred_id{mrna})) {
      @preferred_ids_mrna = $preferred_id{mrna};
    } else {
      @preferred_ids_mrna = qw (locus_tag transcript_id GI);
    }
    &RSAT::message::TimeWarn("Finding preferred IDs for mRNA", join(",", @preferred_ids_mrna)) if ($main::verbose >= 3);
    foreach my $object ($mRNAs->get_objects()) {
      $object->UseAttributesAsID(@preferred_ids_mrna);
    }
    $mRNAs->index_ids();

    ## Use the locus_tag as unique identifier for the other feature types
    foreach my $holder ($rRNAs, $misc_RNAs, $scRNAs, $tRNAs, $genes) {
      my $object_type = lc($holder->get_object_type());
      $object_type =~ s/genbank:://;
      my $preferred_id = $preferred_id{$object_type} || "locus_tag";
      &RSAT::message::TimeWarn("Finding preferred IDs for class", $object_type, $preferred_id) if ($main::verbose >= 3);
      foreach my $object ($holder->get_objects()) {
	$object->UseAttributesAsID($preferred_id);
      }
      $holder->index_ids();
    }

    ## Check object names for all the parsed features, before building the RSAT features from it
    &CheckObjectNames($genes, $mRNAs, $scRNAs, $tRNAs, $rRNAs, $misc_RNAs, $misc_features, $CDSs);

    &InheritGeneNames($CDSs, $mRNAs);

    ## Find a description for the different classes

    ## VERY TRICKY: for the yeast Saccharomyces cerevisiae, the fields
    ## "note" and "product" are intermingled (28 March 2010 and
    ## before). In the same file (chromosome), some genes have the
    ## product name in the field "description" and the descriptin
    ## inthe field "product", and some other genes the opposite. I
    ## TEMPORARILY merge note and product to ensure the full
    ## description.
    &SetDescriptions($CDSs, "note+product","product", "locus_tag", "gene", "name");
    &SetDescriptions($mRNAs, "locus_tag", "gene", "name");
    &SetDescriptions($scRNAs, "product", "locus_tag", "note");
    &SetDescriptions($tRNAs, "product", "locus_tag", "note");
    &SetDescriptions($rRNAs, "product", "locus_tag", "note");

    ## Parse Gene Ontology terms
    &ParseGO($CDSs);

    ################################################################
    ## index names
    for my $holder ($genes,$mRNAs, $scRNAs,$tRNAs,$rRNAs,$misc_RNAs,$misc_features,$CDSs) {
	&RSAT::message::TimeWarn("Indexing IDs and names for class" , $holder->get_object_type()) if ($main::verbose >= 2);
	$holder->index_ids();
	$holder->index_names();
    }

    ################################################################
    ## parse chromosomal positions
    for my $holder ($genes,$mRNAs, $scRNAs,$tRNAs,$rRNAs, $repeat_regions, $misc_RNAs,$misc_features,$CDSs) {
	&ParsePositions($holder);
    }

    ################################################################
    ## Specific treatment for refseq files (probably obsolete)
    if ($data_type eq "refseq") {
	&RefseqPostProcessing();
    } else {
	#### Create features from CDSs and RNAs (obsolete)
	&CreateGenbankFeatures($features, $genes, $mRNAs, $scRNAs, $tRNAs, $rRNAs, $misc_RNAs, $misc_features, $CDSs, $sources, $contigs);
    }
}

################################################################
=pod

=item ParseGenbankFile()

Parse a single Genbank genome file. Create one object per feature.

The sequence is either returned (default), or directly stored in a
file, when this file is specified with the argument seq_file=>myfile.
Direct storage is particularly useful for big genomes (e.g.  human),
where the size of some chromosomes would create an "out of memory"
error.

The option no_seq=>1 prevents from parsing the sequence.

=cut

  sub ParseGenbankFile {
    my ($input_file,
	$features,
	$genes,
	$mRNAs,
	$scRNAs,
	$tRNAs,
	$rRNAs,
	$repeat_regions,
	$misc_RNAs,
	$misc_features,
	$CDSs,
	$contigs,
	$organisms,
	$sources,
	%args) = @_;

    my %features_to_parse = (CDS=>1,
			     source=>1,
			     mRNA=>1,
			     scRNA=>1,
			     tRNA=>1,
			     rRNA=>1,
			     repeat_region=>1,
			     misc_RNA=>1,
			     misc_feature=>1,
			     gene=>1
			    );

    &RSAT::message::TimeWarn("Parsing gbk file",$input_file)
      if ($main::verbose >= 1);
    &RSAT::message::TimeWarn("Parsing features", join (",", keys %features_to_parse))
      if ($main::verbose >= 2);
    &RSAT::message::TimeWarn(join("\t", "Sequence directory", $args{sequence_dir}))
      if ($main::verbose >= 2);

    open GBK, $input_file
      || die "Error: cannot open input file $input_file.\n";
    my $l = 0;
    my $in_sequence = 0;
    my $in_features = 0;
    my $in_feature = 0;
    my $in_gene = 0;
    my $in_cds = 0;
    my $current_feature;
    my $current_contig;
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

    my $l = 0;
    while (my $line = &ReadNextLine()) {
      $l++;
#      &RSAT::message::Debug("line", $l, $line) if ($main::verbose >= 10);

      next unless ($line =~ /\S/);
      unless (($in_features) ||
	      ($in_sequence)){
	#### organism name
	if ($line =~ /^\s+ORGANISM\s+/) {
	  $organism_name = "$'";
	  &RSAT::message::Info("Organism name", $organism_name) if ($main::verbose >= 2);
	  $current_contig->set_attribute("organism", $organism_name);
	  $organism = $organisms->get_object($organism_name);
	  if ($organism) {
	    warn "; Organism $organism already created\n" if ($main::verbose >= 3);
	  } else {
	    $organism = $organisms->new_object();
	    $organism->push_attribute("names", $organism_name);
	    $organism->force_attribute("source", $args{source} || "Genbank");
	    $organisms->index_names(); ### required to prevent creating several objects for the same organism
	  }

	  #### collect the taxonomy
	  my $taxonomy = "";
	  while ($line = &ReadNextLine()) {
	    if ($line =~ /^\S/) {
	      $taxonomy = &trim($taxonomy);
	      $taxonomy =~ s/\s+/ /g;
	      $taxonomy =~ s/\.\s*$//;
	      #			$current_contig->set_attribute("taxonomy", $taxonomy);
	      $organism->force_attribute("taxonomy", $taxonomy);
	      last;
	    } else {
	      $taxonomy .= $line;
	      ## A horrible fix for the fact that there is no clear
	      ## separation between organism name and taxonomy
	      ## For example
	      ##      ORGANISM  Salmonella enterica subsp. enterica serovar Choleraesuis str.
	      ##           SC-B67
	      ##           Bacteria; Proteobacteria; Gammaproteobacteria; Enterobacteriales;
	      ##			 Enterobacteriaceae;
	      ##
	      ## I delete the taxonomy until there is at least one
	      ## semilcolon in the string.
	      unless ($taxonomy =~ /;/) {
		$taxonomy = "";
	      }
	    }
	  }
	}


	if ($line =~ /^([A-Z]+)\s+/) {
	  $current_contig_key = $1;
	  $current_contig_value = "$'";
	  chomp($current_contig_value);
	  $current_contig_value =~ s/\.\s*$//;
	  if ($contig_keys{$current_contig_key}) {
	    warn "parsing\t$current_contig_key\t$current_contig_value\n" if ($verbose >= 4);
	    $current_contig->push_attribute(lc($current_contig_key), $current_contig_value);
	  }

	  ## Use VERSION number as ID
	  if (lc($current_contig_key) eq "version") {
	    @current_contig_values = split / +/, $current_contig_value;
	    $current_contig->force_attribute("id", $current_contig_values[0]);
	  }

	} elsif ($line =~ /^ {12}/) {
	  ### suite of the current contig value
	  $current_contig_value .= " ".$'; ##'
	  warn "parsing\t$current_contig_key\t$current_contig_value\n" if ($verbose >= 4);
	  $current_contig->append_attribute(lc($current_contig_key), $current_contig_value);
	}

	if ($line =~ /^FEATURES/) {
	  $in_features = 1 ;
	  undef($current_contig_key);
	  undef($current_contig_value);
	  warn "; Reading features\n" if ($main::verbose >= 2);
	  next;
	}
      }

      if ($line =~ /^LOCUS/) {
	#### new contig
	@fields = split /\s+/, $line;
	&RSAT::message::Info("New contig", $line) if ($main::verbose >= 2);
	my $contig_name = $fields[1];
	my $length = $fields[2];
	my $type = $fields[4];		#### DNA or peptide
	my $form = $fields[5];		#### circular or linear
	my $taxo_group = $fields[6];	#### short taxonomic group
	my $contig_date = $fields[7];	#### date of the contig
	$current_contig = $contigs->new_object();
	$current_contig->push_attribute("names",$contig_name);
	if ($input_file =~ /[^\/]+\.gbk/) {
	  $current_contig->set_attribute("genbank_file", $&);
	}
	$current_contig->set_attribute("length", $length);
	$current_contig->set_attribute("type", $type);
	$current_contig->set_attribute("form", $form);
	$current_contig->set_attribute("taxo_group", $taxo_group);
	$current_contig->set_attribute("date", $contig_date);

      } elsif ($line =~ /^BASE COUNT/) {
	#### base count
	#### currently ignored

	################################################################
	#### read the full sequence (at the end of the contig)
      } elsif ($line =~ /^ORIGIN/) {
	&RSAT::message::TimeWarn("Reading sequence") if ($main::verbose >= 2);
	$in_features = 0;
	$in_sequence = 1;
	if ($args{no_seq}) {
	  #### skip the sequence
	  while (($line = &readNextLine()) && ($current_contig)) {
	    if ($line =~ /^\/\/$/) {
	      $current_contig = "";
	      last;
	    }
	  }

	  ################################################################
	  ## Save the whole contig sequence in a file
	} elsif ($args{seq_dir}) {
	  my $seq_file = $current_contig->get_attribute("id").".raw";
	  $seq_file =~ s/:/_/g;
	  $current_contig->set_attribute("seq_dir", $args{seq_dir});
	  $current_contig->set_attribute("file", $seq_file);
	  my $seq_file_path = $args{seq_dir}."/".$seq_file;

	  &RSAT::message::Debug ("Contig sequence file", $current_contig, $contig_id, $seq_file) if ($main::verbose >= 2);

	  &RSAT::message::Info (join("\t", "Storing sequence ",
				     $current_contig->get_attribute("id"),
				     "in file", $seq_file_path))
	    if ($main::verbose >= 2);

	  open SEQ, ">".$seq_file_path
	    || die "Error: cannot write sequence file $seq_file_path\n";

	  while (($line = &ReadNextLine()) && ($current_contig)) {
	    if ($line =~ /^\s+\d+\s+/) {
	      $sequence = "$'";
	      $sequence =~ s/\s//g;
	      print SEQ $sequence;
	    } elsif ($line =~ /^\/\/$/) {
	      print SEQ "\n";
	      $in_sequence = 0;
	      $current_contig = "";
#	      last;
	    } else {
	      &ErrorMessage("Invalid sequence format, skipped\t$line\n");
	    }

	    #		    &RSAT::message::Debug("Sequence parsing", length($sequence), $sequence) if (main::verbose >= 10);
	  }
	  close(SEQ);
	  $pwd = `pwd`;
	  chomp($pwd);
	  &RSAT::message::Debug("Working dir", $pwd,  "Sequence saved in file", $seq_file_path) if ($main::verbose >= 2);

	} else {
	  #### load sequence in memory and return it
	  while ($line = &ReadNextLine()) {
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

	################################################################
	#### new feature
      } elsif ($line =~ /^ {5}(\S+)\s+/) {

	$feature_type = $1;
	$value = "$'";


	#### parse feature  position
	$position = &trim($value);
	if ($position =~ /join\(/) {
	  ### check that the position is complete
	  my $start_line = $l;
	  my $next = "$'";;
	  my $end_expression = '\)';
	  if ($position =~ /join\(complement\(/) {
	    $end_expression = '\)\)';
	  }

	  unless ($next =~ /${end_expression}/) {
	    do {
	      die "Error: position starting at line $l is not terminated properly.\n"
		unless $position_suite = &ReadNextLine();
	      $position_suite =~ s/^FT//;
	      $position .= &trim($position_suite);
	    } until ($position =~ /${end_expression}/);
	  }
	}

	#	    &RSAT::message::Debug("feature type", $feature_type, $value, $position) if ($main::verbose >= 10);

	#### create an object for the new feature
	if ($features_to_parse{$feature_type}) {
	  &RSAT::message::Debug($l, "new feature", $feature_type, $position) if ($main::verbose >= 3);

	  if ($feature_type eq "gene") {
	    $holder = $genes;
	  } elsif ($feature_type eq "mRNA") {
	    $holder = $mRNAs;
	  } elsif ($feature_type eq "scRNA") {
	    $holder = $scRNAs;
	  } elsif ($feature_type eq "tRNA") {
	    $holder = $tRNAs;
	  } elsif ($feature_type eq "rRNA") {
	    $holder = $rRNAs;
	  } elsif ($feature_type eq "repeat_region") {
	    $holder = $repeat_regions;
	    #		    &RSAT::message::Debug("repeat region", $feature_type) if ($main::verbose >= 10);
	    #		    die "HELLO";
	  } elsif ($feature_type eq "misc_RNA") {
	    $holder = $misc_RNAs;
	  } elsif ($feature_type eq "misc_feature") {
	    $holder = $misc_features;
	  } elsif ($feature_type eq "CDS") {
	    $holder = $CDSs;
	  } elsif ($feature_type eq "source") {
	    $holder = $sources;
	  } else {
	    $holder = $features;
	  }
	  $current_feature = $holder->new_object();

	  $current_feature->set_attribute("type",$feature_type) ;

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
		   ($feature_type eq 'scRNA') ||
		   ($feature_type eq 'tRNA') ||
		   ($feature_type eq 'rRNA') ||
		   ($feature_type eq 'misc_RNA')) {

	    ## Link the current feature to its parent gene
	    if ($last_gene) {
	      my $last_gene_id =  $last_gene->get_attribute("id");
	      $current_feature->set_attribute("gene_id", $last_gene_id);

	      #			&RSAT::message::Debug("\t", "feature", $current_feature->get_attribute("id"),
	      #				   "parent gene", $last_gene->get_attribute("id"),
	      #				   "gene_id", $last_gene_id,
	      #				  ), "\n"if ($main::verbose >= 5);


	      #### set the type of the last gene to the
	      #### current feature type, if it has not yet
	      #### been done
	      my $last_gene_type = $last_gene->get_attribute("type");
	      if (($last_gene_type eq "") ||
		  ($last_gene_type eq $main::null) ||
		  ($last_gene_type eq "gene")) {
		$last_gene->force_attribute("type", $feature_type);
	      }
	    }
	  }


	  ## Other feature types are currently ignored
	} else {
	  &RSAT::message::Info(join("\t", $l, "Ignoring feature", $feature_type, $position)) if ($main::verbose >= 2);
	  undef($current_feature);
	}

	#### new feature attribute

	### A new feature attribute is detected by
	### - 21 spaces
	### - followed by a slash (/),
	### - followed by a word (the name of the attribute),
	### - followed by the attribute value

      } elsif ($line =~ /^ {21}\/(\S+)=/) {
	$attribute_type=$1;
	$attribute_value=&trim("$'");

	#### string attributes : make sure that the entire string is parsed
	if ($attribute_value =~ /^\"/) {
	  unless ($attribute_value =~ /\"$/) {
	    do {
	      $line = &ReadNextLine();
	      $line = &trim($line);
	      $attribute_value .= " ";
	      $attribute_value .= &trim($line);
	      #			&RSAT::message::Debug("\t", $l, "Collecting string attribute", $attribute_type, $attribute_value), "\n" if ($main::verbose >= 10);
	    } until ($line =~ /\"$/);
	  }
	  $attribute_value =~ s/^\"//;
	  $attribute_value =~ s/\"$//;
	}
	if ($current_feature) {


#	  &RSAT::message::Debug("new attribute", $attribute_type, $attribute_value) if ($main::verbose >= 10);

	  ## the attribute "type" is reserved, since all the features have a type (source, gene, CDS, ...).
	  if ($attribute_type eq "type") {
	      $attribute_type = "type2";
	  }
	  $current_feature->new_attribute_value($attribute_type, $attribute_value);
	  #		&RSAT::message::Debug( "\t", $l, "\tadding attribute", $attribute_type, $attribute_value) if ($main::verbose >=3);


	  ## Detect GeneID
	  if (($attribute_type eq "db_xref") && ($attribute_value =~ /^GeneID:(\d+)/i)) {
	    my $GeneID = $1;
	    if ($GeneID) {
	      $current_feature->force_attribute("GeneID", $GeneID);
	      &RSAT::message::Warning(join("\t", "Assigning GeneID",$GeneID,
					   "to feature", $current_feature->get_attribute("id"),
					   "type",  $current_feature->get_attribute("type"),
					  )) if ($main::verbose >= 5);
	      if (($current_feature->get_attribute("type") eq "gene") ) {
		$current_feature->ReplaceID("GeneID");
	      }
	    } else {
	      &RSAT::message::Warning(join("\t", "Empty GeneID for feature", $current_feature->get_attribute("id")));
	    }
	  }

	  #### detect taxid
	  if (($current_feature->get_attribute("type") eq "source") ) {
	    if (($attribute_type eq "db_xref") && ($attribute_value =~ /taxon:(\d+)/)) {
	      $taxid = $1;
	      $current_contig->force_attribute("taxid", $taxid);
	      $current_contig->force_attribute("taxid", $taxid);
	      $organism->force_attribute("id", $taxid);
	      foreach my $class ($contigs,
				 $features,
				 $genes,
				 $mRNAs,
				 $scRNAs,
				 $tRNAs,
				 $rRNAs,
				 $repeat_regions,
				 $misc_RNAs,
				 $misc_features,
				 $CDSs,
				 $sources,
				 $organisms) {
		$class->set_prefix($taxid);
	      }
	    }
	  }

	} else {
	  &RSAT::message::Warning(join("\t", $l, "\tignoring attibute", $attribute_type, $attribute_value)) if ($main::verbose >=5);
	}


      } else {
	&RSAT::message::Warning("file ".$short_name{$org}."line",$l,"\tnot parsed",$line) if ($main::verbose >= 3);
      }
    }
    close GBK;


    return ($file_description, $sequence);
  }


################################################################
=pod

=item ParseGO()

Extract GO identifiers from the feature notes

=cut

sub ParseGO {
    my ($CDSs) = @_;
    &RSAT::message::TimeWarn("Parsing Gene Ontology for CDS") if ($main::verbose >= 1);
    foreach my $cds ($CDSs->get_objects()) {
	foreach my $note ($cds->get_attribute("note")) {
	    if ($note =~ /\[goid (\S+)]/) {
		my $goid = $1;

		## The GO: is present some genomes
		## (e.g. S.cerevisiae), absent from others
		## (e.g. E.coli K12). I am pretty sure that it will
		## sometimes be in lowercases, sometimes in
		## uppercases. I check them for the sake of robustness.
		if ($goid =~ /^\d+$/) {
		  $goid = "GO:".$goid;
		} elsif ($goid =~ /^GO:(\d+)/i) {
		  $goid = uc($goid);
		} else {
		  &RSAT::message::Warning(join("\t", "Invalid GO identifier", $goid));
		}
#		$goid =~ s/GO\://i;
		$cds->push_attribute("GO", $goid);
	    }
	}
    }

}

################################################################
=pod

=item ParseFeatureNames()

Parse some names from the notes. Guess which "notes" actually
correspond to synonyms.

The basic guessing rule is : any note made of a single word is
considered as a synonym (this is of course not optimal)

For some genomes, additional identifiers/synonyms have specific
formats (synonyms, locus_tag, Accession, ...).

=cut
sub ParseFeatureNames {
    my (@class_holders) = @_;

    foreach my $class_holder (@class_holders) {
	&RSAT::message::TimeWarn("Parsing feature names for class", $class_holder->get_object_type())
	    if ($verbose >= 2);

	foreach my $feature ($class_holder->get_objects()) {

	    ################################################################
	    ## Attribute "synonym" (e.g. C.elegans genome version July 2003)
	    foreach my $synonym_line ($feature->get_attribute("synonym")) {
		my @synonyms = split /\,\s+/, $synonym_line;
#		&RSAT::message::Debug("feature synonym field", $feature->get_attribute("id"), scalar(@synonyms), join ("; ", @synonyms)) if ($$main::verbose >= 10);
		foreach my $synonym (@synonyms) {
		    $feature->push_attribute("names", $synonym);
		}
	    }

	    ################################################################
	    ## Attribute "gene"
	    my @genes = $feature->get_attribute("gene");
	    foreach my $gene_name (@genes) {
		$new_name = &trim($gene_name);
		&RSAT::message::Debug( "\t", "Adding gene name as name to feature",
			   $feature->get_attribute("id"),
			   $new_name
			   ), "\n" if ($verbose >= 4);
		$feature->push_attribute("names", $new_name);
	    }


	    ################################################################
	    ## Attribute "locus_tag"
	    my ($locus_tag) = $feature->get_attribute("locus_tag");
	    $locus_tag = &trim($locus_tag); ## suppress leading and trainling spaces
	    &RSAT::message::Debug( "\t", "Adding locus tag as name to feature",
				   $feature->get_attribute("id"),
				   $locus_tag
				   ), "\n" if ($verbose >= 4);
	    $feature->push_attribute("names", $locus_tag);

	    ################################################################
	    ## Attribute "GeneID"
	    my ($GeneID) = $feature->get_attribute("GeneID");
	    $GeneID = &trim($GeneID); ## suppress leading and trainling spaces
	    &RSAT::message::Debug( "\t", "Adding locus tag as name to feature",
				   $feature->get_attribute("id"),
				   $GeneID
				   ), "\n" if ($verbose >= 4);
	    $feature->push_attribute("names", $GeneID);




	    ################################################################
	    ## Some types of notes contain identifiers or synonyms
	    my @notes = $feature->get_attribute("note");
	    foreach my $note (@notes) {
		$note = &trim($note); ### remove leading and trailing spaces
		warn "NOTE\t'$note'\n" if ($verbose >= 4);

		if ($note =~ /synonyms:/) {
		    ## For some genomes there is a note of type 'synonyms:' (e.g. Saccharomyces cerevisiae),
		    my $synonyms = $'; ##'
		    my @synonyms = split /, /, $synonyms;
		    foreach my $new_name (@synonyms) {
			$new_name = &trim($new_name);
			&RSAT::message::Debug( "\t", "Adding synonym to feature",
				   $feature->get_attribute("id"),
				   $new_name
				   ), "\n" if ($verbose >= 4);
			$feature->push_attribute("names", $new_name);
		    }

		} elsif ($note =~ /synonym:/) {
		    ## For other genomes there is a note of type 'synonym:' (e.g. Arabidopsis thaliana),
		    my $synonyms = $';  ##'
		    $synonyms =~ s/;.*//; ## In A.thaliana there are comments after the synonym

		    my @synonyms = split /, /, $synonyms;
		    foreach my $new_name (@synonyms) {
			$new_name = &trim($new_name);
			&RSAT::message::Debug( "\t", "Adding synonym to feature",
				   $feature->get_attribute("id"),
				   $new_name
				   ), "\n" if ($verbose >= 4);
			$feature->push_attribute("names", $new_name);
		    }

		} elsif ($note !~ /\S/) {
		    ## A single-word note is usually (but not always, I guess) a synonym
		    $feature->push_attribute("names", $note);

		} elsif ($org eq "Plasmodium_falciparum") {
		    my @notes = $feature->get_attribute("note");
		    foreach my $note (@notes) {
			next if ($note =~/score/);
			my @names = split /\;\s*/, $note;
			foreach my $name (@names) {
			    $name = &RSAT::util::trim($name);
			    if (($name =~ /^(PF\S+)/) || ($name =~ /^(MAL\S+)/)) {
				$feature->push_attribute("names", $1);
# 				&RSAT::message::Debug ("\t",
# 					   "Plasmodium_falciparum",
# 					   $feature->get_attribute("contig"),
# 					   $feature->get_attribute("id"),
# 					   $feature->get_attribute("type"),
# 					   , "adding name", $1, "from note", $note), "\n" if ($main::verbose >= 10);
			    }
			}
		    }
#		    die "Specific treatment for Plasmodium_falciparum";


		} elsif (
		    ## For some genomes, the locus_tag, Accession or other IDs have been annotated as "note" !
		    ($note =~ /^locus_tag\: (\S+)/i)  ||
		    ($note =~ /^Accession ([\w_\-\.]+)/) ||
		    ($note =~ /^(\S+)\,\s+len:/)
		){

		    $new_name = $1;
		    $new_name =~ s/\;$//;
		    $new_name =~ s/\:$//;

		    #### add the new name
		    &RSAT::message::Debug( "\t", "Adding name to feature",
			       $feature->get_attribute("id"),
			       $new_name
			       ), "\n" if ($verbose >= 4);
		    $feature->push_attribute("names", $new_name);

		}

	    }

	}
    }
}





################################################################

=pod

=item  CreateGenbankFeatures()

After having features of different types from Genbank, create
RSAT-formatted features for specific types (CDS, mRNA, tRNA, ...).

=cut
sub CreateGenbankFeatures {
  my ($features, $genes, $mRNAs, $scRNA, $tRNAs, $rRNAs, $misc_RNAs, $misc_features, $CDSs, $sources,
      #	$repeat_regions,
      $contigs) = @_;

  ## Index gene names
  #    &RSAT::message::TimeWarn("Indexing gene names") if ($main::verbose >= 2);
  #    $genes->index_ids();
  #    $genes->index_names();


  #### extract taxid for each source object
  foreach my $source ($sources->get_objects()) {
    $source->get_taxid();
  }

  #### Create features for contig limits
  &RSAT::message::TimeWarn("Creating features for contig limits")
    if ($main::verbose >= 3);
  foreach my $contig ($contigs->get_objects()) {
    my $contig_length = $contig->get_attribute("length");

    #  	## Contig start position
    #  	my $start_feature = $features->new_object(%args);
    #  	$start_feature->set_attribute("type","SEQ_START");
    #  	$start_feature->set_attribute("name","SEQ_START");
    #  	$start_feature->set_attribute("description",$contig->get_attribute("id")."; contig start");
    #  	$start_feature->set_attribute("contig",$contig->get_attribute("id"));
    #  	$start_feature->set_attribute("organism",$contig->get_attribute("organism"));
    #  	$start_feature->set_attribute("chrom_position","1..1");
    #  	$start_feature->set_attribute("start_pos",1);
    #  	$start_feature->set_attribute("end_pos",1);
    #  	$start_feature->set_attribute("strand","DR");

    ## Contig start position
    my $end_feature = $features->new_object(%args);
    $end_feature->set_attribute("type","SEQ_END");
    $end_feature->set_attribute("name","SEQ_END");
    $end_feature->set_attribute("description",$contig->get_attribute("id")."; contig end");
    $end_feature->set_attribute("contig",$contig->get_attribute("id"));
    $end_feature->set_attribute("organism",$contig->get_attribute("organism"));
    $end_feature->set_attribute("chrom_position",$contig_length."..".$contig_length);
    $end_feature->set_attribute("start_pos",$contig_length);
    $end_feature->set_attribute("end_pos",$contig_length);
    $end_feature->set_attribute("strand","DR");

  }

  ################################################################
  ## Create unified features
  &RSAT::message::TimeWarn("Creating unified features from parsed features")
    if ($main::verbose >= 2);

  foreach my $parent_feature ($CDSs->get_objects(),
			      $mRNAs->get_objects(),
			      $scRNAs->get_objects(),
			      $tRNAs->get_objects(),
			      $rRNAs->get_objects(),
			      #				$misc_RNAs->get_objects(),
			      #				$misc_features->get_objects()
			     ) {
    &RSAT::message::TimeWarn(join ("\t", "Creating feature for parsed feature",
				   $parent_feature->get_attribute("id"),
				   "type", $parent_feature->get_attribute("type"),
				  ))
      if ($main::verbose >= 3);

    ## Create a new feature from the parsed feature
    $created_feature = $features->new_object(%args);
    $created_feature->force_attribute("id",$parent_feature->get_attribute("id"));
    $created_feature->set_attribute("type",$parent_feature->get_attribute("type"));
    $created_feature->set_attribute("organism",$parent_feature->get_attribute("organism"));
    $created_feature->set_attribute("contig",$parent_feature->get_attribute("contig"));
    $created_feature->set_attribute("chrom_position",$parent_feature->get_attribute("chrom_position"));
    $created_feature->set_attribute("start_pos",$parent_feature->get_attribute("start_pos"));
    $created_feature->set_attribute("end_pos",$parent_feature->get_attribute("end_pos"));
    $created_feature->set_attribute("strand",$parent_feature->get_attribute("strand"));

    $created_feature->set_attribute("locus_tag",$parent_feature->get_attribute("locus_tag"));
    $created_feature->set_attribute("GeneID",$parent_feature->get_attribute("GeneID"));

    &RSAT::message::Debug ("\t", "feature",
			   $created_feature->get_attribute("id"),
			   "type", $created_feature->get_attribute("type"),
			   "organism", $created_feature->get_attribute("organism"),
			  ), "\n" if ($verbose >= 4);


    # 	################################################################
    # 	#### Define names for the new feature
    # 	&RSAT::message::Debug ("feature",
    # 			       $created_feature->get_attribute("id"),
    # 			       "Adding gene names",
    # 			      ), "\n" if ($verbose >= 4);
    # 	my $gene_name = "";
    # 	my $ParsedGeneID = $parent_feature->get_attribute("GeneID");
    # 	my $gene = "";
    # 	if ($ParsedGeneID) {
    # 	    ## Primary gene name is the one documented as "gene" attribute in the feature iself
    # 	    $gene_name = $parent_feature->get_attribute("gene");
    # 	    if ($gene_name) {
    # 		$created_feature->force_attribute("name",$gene_name);
    # 		$created_feature->push_attribute("names",$gene_name);
    # 	    }

    # 	    ## Identify the parent gene
    # 	    $ParsedGeneID = $parent_feature->get_attribute("GeneID");
    # 	    &RSAT::message::Debug("feature",
    # 				  $parent_feature->get_attribute("id"),
    # 				  "type", $parent_feature->get_attribute("type"),
    # 				  "GeneID",$ParsedGeneID,
    # 				 ) if ($main::verbose >= 4);
    # 	    unless ($ParsedGeneID) {
    # 		&RSAT::message::Warning(join("\t", "There is no GeneID for feature",
    # 					     $parent_feature->get_attribute("id"),
    # 					     "type", $parent_feature->get_attribute("type")));
    # 		next;
    # 	    }
    # 	    if ($ParsedGeneID eq $main::null) {
    # 		&RSAT::message::Warning(join("\t", "GeneID is null for feature",
    # 					     $parent_feature->get_attribute("id"),
    # 					     "type", $parent_feature->get_attribute("type")));
    # 		next;
    # 	    }

    # 	    ## Add parent gene names to the current feature
    # 	    $gene = $genes->get_object($ParsedGeneID);
    # 	    if ($gene) {
    # 		&RSAT::message::Debug("Feature parent gene", $created_feature->get_attribute("id"),
    # 				      "gene name", $gene_name,
    # 				      "GeneID", $ParsedGeneID,
    # 				      "gene ID", $gene->get_attribute("id"),
    # 				     ), "\n" if ($verbose >= 10);

    # 		#### use gene name as primary name for the feature
    # 		my $gene_name = $gene->get_name();
    # 		if (($gene_name) && ($gene_name ne $main::null)) {
    # 		    $created_feature->force_attribute("name", $gene_name);
    # 		    $created_feature->push_attribute("names", $gene_name);
    # 		}

    # 	    	#### accept all gene names for the new feature
    # 		foreach my $name ($gene->get_attribute("names")) {
    # 		    $created_feature->push_attribute("names",$name) unless ($names{$name});
    # 		}

    # 	    } else {
    # 		&RSAT::message::Warning(join("\t", "Cannot identify gene",$ParsedGeneID,
    # 					     "for feature", $parent_feature->get_attribute("id"),
    # 					     "type", $parent_feature->get_attribute("type")));
    # 		next;
    # 	    }
    # 	}


    ################################################################
    #### Inherit names from the parent feature
    foreach my $name ($parent_feature->get_attribute("names")) {



      $created_feature->push_attribute("names",$name);
    }
    &RSAT::message::Debug ("\t",  "feature",
			   $created_feature->get_attribute("id"),
			   "Added names from original feature",
			   $parent_feature->get_attribute("id"),
			   join(";", $parent_feature->get_attribute("names"))
			  ) if ($verbose >= 5);

    ################################################################
    #### Inherit cross-references from parent feature
    &RSAT::message::Debug ("\t",  "feature",
			   $created_feature->get_attribute("id"),
			   "Cross references",
			  ), "\n" if ($verbose >= 5);



    my @xrefs = $parent_feature->get_attribute("db_xref");
    $created_feature->set_array_attribute("db_xref", @xrefs);



    &ExtractCrossReferencesForFeature($created_feature);



    ################################################################
    ## Inherit locus tag from parent feature
    my @locus_tags = $parent_feature->get_attribute("locus_tag");




    &RSAT::message::Debug( "\t",
			   "parsed feature", $parent_feature->get_attribute("id"),
			   "locus tag", join ";", @locus_tags) if ($main::verbose >= 5);

    ## Add locus tags to the list of synonyms
    foreach my $locus_tag (@locus_tags) {
      &RSAT::message::Debug ("\t", "Feature", $created_feature->get_attribute("id"),
			     "Adding locus tag as synonym", $locus_tag),  "\n"
			       if ($main::verbose >= 5);
      $created_feature->push_attribute("locus_tags", $locus_tag);
      $created_feature->push_attribute("names", $locus_tag);

    }

    ################################################################
    #### Define a single name  (take the first value in the name list)
    $single_name = 1;
    if ($single_name) {
      my ($name) = $created_feature->get_attribute("names"); ## Use first name as primary
      if ($name) {
	$created_feature->force_attribute("name",$name);
      } else {
	$created_feature->force_attribute("name",$created_feature->get_attribute("id"));
      }
      &RSAT::message::Debug ("\t", "feature",
			     $created_feature->get_attribute("id"),
			     "single name", $name,,
			    ), "\n" if ($verbose >= 5);
    }

    ################################################################
    ## Create a description for the new feature

    ## Use parent feature description
    my $description = $parent_feature->get_attribute("description");
    if (($description) && ($description ne $main::null)) {
      $created_feature->set_attribute("description", $description);
    } else {
      ## Use the parent feature product field if any
      my @products = $parent_feature->get_attribute("product");
      if (($products[0]) && ($products[0] ne $main::null)) {
	$created_feature->set_attribute("description",$products[0]);
      } else {
	#### use the first parsed feature note
	my @feature_notes = $parent_feature->get_attribute("note");
	if ($feature_notes[0]) {
	  $created_feature->set_attribute("description",$feature_notes[0]);
	} else {
	  #### if there is no feature note, use the gene note
	  my $GeneID = $parent_feature->get_attribute("DeneID");
	  if (($GeneID) && ($GeneID ne $main::null)) {
	    my $gene = $genes->get_object($GeneID);
	    my @gene_notes = $gene->get_attribute("note");
	    if ($gene_notes[0]) {
	      $created_feature->set_attribute("description",$gene_notes[0]);
	    }
	  }
	}
      }
    }
    &RSAT::message::Debug ("\t", "feature",
			   $created_feature->get_attribute("id"),
			   "description", $created_feature->get_attribute("description"),
			  ), "\n" if ($verbose >= 5);

    ################################################################
    #### protein ID (for CDS)
    if ($parent_feature->get_attribute("type") eq "CDS") {
      my @protein_ids= $parent_feature->get_attribute("protein_id");
      foreach my $protein_id (@protein_ids) {
	if (($protein_id) && ($protein_id ne $main::null)) {
	  ### remove the version number
	  if ($protein_id =~ /(.*)\.\d+$/) {
	    $no_version = $1;
	    $created_feature->push_attribute("names", $no_version);
	    $parent_feature->push_attribute("names", $no_version);
	  } else {
	    $created_feature->push_attribute("names", $protein_id);
	    $parent_feature->push_attribute("names", $protein_id);
	  }
	}
      }

      ################################################################
      #### transcript ID (for mRNA)
    } elsif ($parent_feature->get_attribute("type") eq "mRNA") {
      my @transcript_ids= $parent_feature->get_attribute("transcript_id");
      foreach my $transcript_id (@transcript_ids) {
	if (($transcript_id) && ($transcript_id ne $main::null)) {
	  ### remove the version number
	  if ($transcript_id =~ /(.*)\.\d+$/) {
	    $no_version = $1;
	    $created_feature->push_attribute("names", $no_version);


	    $parent_feature->push_attribute("names", $no_version);
	  } else {
	    $created_feature->push_attribute("names", $transcript_id);
	    $parent_feature->push_attribute("names", $transcript_id);
	  }
	}
      }
    }
  }
}

################################################################

=pod

=item  CheckObjectNames()

Try to extract as many names as possible, make sure that each obect
has at least one name, and that there are no duplicate names for the
same object.

=cut
sub CheckObjectNames {
    my (@holders) = @_;


    ## Each feature inherits the names of its parent gene
    foreach my $holder (@holders) {
	&RSAT::message::TimeWarn("Checking object names for class",
				 $holder->get_object_type()) if ($main::verbose >= 2);

	foreach my $object ($holder->get_objects()) {
#	  &RSAT::message::Debug("Checking object names for object",
#				$holder->get_object_type(),
#				"id=".$object->get_attribute("id"),
#				"name=".$object->get_attribute("name"),
#				"gene=".$object->get_attribute("gene"),
#				"names=", join(";", $object->get_attribute("names")),
#				  ) if ($main::verbose >= 10);

	    ## Make sure that the object has no null name
	    my $primary_given = 0;
	    my $primary_name = $object->get_attribute("name");

	    if (($primary_name) && ($primary_name ne $main::null)) {
		$primary_given = 1;

	    }



	  ## Test different attributes which can be considered as primary name
	  for my $attribute ("gene", "synonym", "synonyms", "locus_tag", "GeneID", "id") {
	    my @names = $object->get_attribute($attribute);
	    foreach my $name (@names) {
	      if (($name) && ($name ne $main::null)) {
		#### add gene name to the list of synonyms
		$object->push_attribute("names", $name);

		unless ($primary_given) {
		  #### define this name as primary name
		  $object->set_attribute("name",  $name);
#		  &RSAT::message::Debug("Assigned field", $attribute, "value", $name,
#					"as primary name for object",
#					$object->get_attribute("id")) if ($main::verbose >= 10);
		  $primary_given = 1;
		}
	      }
	    }
	  }

# 	    ## The feature attribute "gene" is used as name
# 	    my $gene_attr = $object->get_attribute("gene");
# 	    if ($gene_attr) {
# 		#### add gene name to the list of synonyms
# 		$object->push_attribute("names", $gene_attr);

# 		unless ($primary_given) {
# 		    #### define this name as primary name
# 		    $object->set_attribute("name",  $gene_attr);
# 		    $primary_given = 1;
# 		}
# 	    }

# 	    ## The feature attribute "locus_tag" is used as name and
# 	    ## becomes primary name if there is no attribute "gene".
# 	    my @locus_tags = $object->get_attribute("locus_tag");
# 	    foreach my $locus_tag_attr (@locus_tags) {
# 		#### add locus_tag name to the list of synonyms
# 		$object->push_attribute("names", $locus_tag_attr);

# 		unless ($primary_given) {
# 		    #### define this name as primary name
# 		    $object->set_attribute("name",  $locus_tag_attr);
# 		    $primary_given = 1;
# 		}
# 	    }

	    ################################################################
# 	    #### Add all the names of the parent gene to the current feature
# 	    my $ParentGeneID = $object->get_attribute("GeneID");
# 	    if (($ParentGeneID) && ($ParentGeneID ne $main::null)) {
# 		my $gene = $genes->get_object($ParentGeneID);
# 		if ($gene) {
# 		    foreach my $gene_name ($gene->get_attribute("names")) {
# 			$object->push_attribute("names", $gene_name);
# 		    }
# 		} else {
# 		    &ErrorMessage("Cannot identify gene $ParentGeneID for object\t", $object->get_attribute("id"), "\n");
# 		}
# 	    }

# 	    #### use ID as primary name
# 	    unless ($primary_given) {
# 		$object->set_attribute("name",  $object->get_attribute("id"));
# 		$primary_given = 1;
# 	    }

	    ## Make sure that object has no duplicate names
	    $object->unique_names();
	}
    }

}


################################################################
## Set the description for a given class of objects, by taking the
## value of one or several specified field. If several fields are
## specified, the first field having a non-null value is used as
## description.
sub SetDescriptions {
  my ($holder, @fields) = @_;
  &RSAT::message::TimeWarn("Setting descriptions for class",
			   $holder->get_object_type()) if ($main::verbose >= 2);
  foreach my $object ($holder->get_objects()) {
    #    &RSAT::message::Debug("Searching description  for object",
    #			  $holder->get_object_type(),
    #			  "id=".$object->get_attribute("id"),
    #			  "name=".$object->get_attribute("name"),
    #			 ) if ($main::verbose >= 10);

    ## Make sure that the object has no null name
    my $description_given = 0;
    my $description = $object->get_attribute("description");
    if (($description) && ($description ne $main::null)) {
      $description_given = 1;
    }

    ## Test different attributes which can be used as description
    my $id = $object->get_attribute("id");
#    &RSAT::message::Debug("SetDescriptions", $id) if ($main::verbose >= 10);
    for my $field (@fields) {
      my @fields_to_merge = split /\+/, $field;
#      &RSAT::message::Debug("", "Field", $field, join (" +++ ", @fields_to_merge)) if ($main::verbose >= 10);
      my @values = ();
      foreach my $to_merge (@fields_to_merge) {
	push @values, $object->get_attribute($to_merge);
      }
      my $value = join("; ", @values); ## Merge multi-value fields into a single string
      #      &RSAT::message::Debug("", "attribute=".$field, $value) if ($main::verbose >= 10);
      if (($value) && ($value ne $main::null)) {
	unless ($description_given) {
	  #### define this field as description
	  $object->force_attribute("description",  $value);
#	  &RSAT::message::Debug("Assigned description", $object->get_attribute("id"), "field", $field, "value", $value) if ($main::verbose >= 10);
	  $description_given = 1;
	}
      }
    }
  }
}


################################################################
=pod

=itme B<ExtractCrossReferences>

Extract information from the cross-references for al the objects of a
class holder.

=cut
sub ExtractCrossReferences {
  my (@class_holders) = @_;

  foreach my $class_holder (@class_holders) {
    &RSAT::message::TimeWarn("Extracting cross-references for class", $class_holder->get_object_type())
      if ($verbose >= 2);

    foreach my $feature ($class_holder->get_objects()) {
      &ExtractCrossReferencesForFeature($feature);
    }
  }
}


################################################################
=pod

=itme B<ExtractCrossReferences>

Extract information from the cross-references for a single feature
object.

=cut
sub ExtractCrossReferencesForFeature {
  my ($feature) = @_;
  &RSAT::message::Info("Extracting cross-references for feature", $feature->get_attribute("id"))
    if ($verbose >= 3);
  my $gi;
  my $GeneID;
  my $LocusID;
  my $locus_tag;
  my @xrefs = $feature->get_attribute("db_xref");
  foreach my $xref (@xrefs) {
    #### extract GI from cross-references
    if ($xref =~ /GI:/) {
      $gi = $';  ##'

      #### accept GI as synonym
      $feature->push_attribute("names", $gi);
      $feature->force_attribute("GI", $gi);

    } elsif ($xref =~ /GeneID:/) {
      $GeneID = $';  ##'
      $feature->force_attribute("GeneID", $GeneID);
      #### accept GeneID as synonym
      $feature->push_attribute("names", $GeneID);

    } elsif ($xref =~ /LocusID:/) {
      $LocusID = $';  ##'
      #### accept LocusID as synonym
      $feature->push_attribute("names", $LocusID);

    } elsif ($xref =~ /locus_tag:/) {
      ## Locus tags were previously annotated as cross-references rather than direct attributes
      $locus_tag = $';  ##'
      #### accept locus_tag as synonym
      $feature->push_attribute("locus_tags", $locus_tag);
      $feature->push_attribute("names", $locus_tag);
    }
  }

  &RSAT::message::Debug ("\t", ";",
			 "parsed feature", $feature->get_attribute("id"),
			 "GI=$gi",
			 "locus_tag=$locus_tag",
			 "LocusID=$LocusID",
			 "product=$products[0]",
			), "\n" if ($verbose >= 3);

}

################################################################
## CDSs and mRNAs inherit names from their parent gene
sub InheritGeneNames {
  my (@class_holders) = @_;

  ## Index genes in order to be able to retrieve them
  $genes->index_ids();
  $genes->index_names();

  foreach my $class_holder (@class_holders) {
    &RSAT::message::TimeWarn("Extracting cross-references for class", $class_holder->get_object_type())
      if ($verbose >= 2);

    foreach my $feature ($class_holder->get_objects()) {
      my $GeneID = $feature->get_attribute("GeneID");

      if (($GeneID eq "") || ($GeneID eq $main::null)) {
	&RSAT::message::Warning("No GeneID for feature", $feature->get_attribute("id")) if ($main::verbose >= 2);
      } else {
	my $gene = $genes->get_object($GeneID);
	if ($gene) {
	  my @gene_names = $gene->get_attribute("names");
	  foreach my $name (@gene_names) {
	    $feature->push_attribute("names", $name);
	    &RSAT::message::Info("Added gene name", $name, "to feature", $feature->get_attribute("id")) if ($main::verbose >= 4);
	  }
	} else {
	  &RSAT::message::Warning("No gene corresponding to GeneID", $GeneID, "for feature", $feature->get_attribute("id")) if ($main::verbose >= 2);
	}
      }
    }
  }
}

1;


__END__

=pod

=back
