#!/usr/bin/perl
############################################################
#
# $Id: parse_genbank_lib.pl,v 1.2 2004/04/01 22:14:42 jvanheld Exp $
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

=item ParseGenbankFile()

Parse a Genbank genome file. Create one object per feature.

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
	$tRNAs, 
	$rRNAs, 
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
		$organism = $organisms->get_object($organism_name);
		if ($organism) {
		    warn "; Organism $organism already created\n" if ($main::verbose >= 3);
		} else {
		    $organism = $organisms->new_object();
		    $organism->push_attribute("names", $organism_name);
		    $organisms->index_names(); ### required to prevent creating several objects for the same organism
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
	    $current_contig = $contigs->new_object(id=>$contig_id);
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
		    $holder = $genes;
		} elsif ($feature_type eq "mRNA") {
		    $holder = $mRNAs;
 		} elsif ($feature_type eq "tRNA") {
 		    $holder = $tRNAs;
 		} elsif ($feature_type eq "rRNA") {
 		    $holder = $rRNAs;
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
		    
		    ## Link the current feature to its parent gene
		    if ($last_gene) {
			my $gene_id =  $last_gene->get_attribute("id");
			$current_feature->set_attribute("gene_id", $gene_id);
			
			warn join ("\t", "feature", $current_feature->get_attribute("id"),
				   "parent gene", $last_gene->get_attribute("id"),
				   "gene_id", $gene_id,
				   ), "\n"if ($main::verbose >= 5);

			
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
		warn join "\t", $l, "ignoring ft", $feature_type, $position, "\n" if ($main::verbose >= 2);
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
			foreach my $classs ($contigs,
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
			    $classs->set_prefix($taxid);
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
    #### this only applies to mRNA and CDS because the other features
    #### (gene, tRNA, rRNA, misc_RNA, misc_feature) have no IDs in Genbank flat files !!!!
    foreach my $object ($mRNAs->get_objects(),
			$CDSs->get_objects()
			) {
	$object->UseGenbankGIasID();
    }

    &ParseFeatureNames($genes, $mRNAs, $CDSs);
    &ParseFeatureNames($CDSs);
    
    return ($file_description, $sequence);
}


################################################################
=pod

=item ParseGO()

Extract GO identifiers from the feature notes

=cut

sub ParseGO {
    my ($CDSs) = @_;
    warn "ParseGO: TO BE IMPLEMENTED";
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
    warn ("; ",
	  &AlphaDate(),
	  "\tGuessing feature synonyms\n")
	if ($verbose >= 1);
    
    foreach my $class_holder (@class_holders) {

	warn "; Adding names to features of type ", $class_holder->get_object_type(), "\n"
	    if ($main::verbose >= 1);

	foreach my $feature ($class_holder->get_objects()) {
	    
	    ## Attribute "gene"
	    my @genes = $feature->get_attribute("gene");
	    foreach my $gene_name (@genes) {
		$new_name = &trim($gene_name);
		warn join( "\t", "Adding gene name as name to feature", 
			   $feature->get_attribute("id"), 
			   $new_name
			   ), "\n" if ($verbose >= 4);
		$feature->push_attribute("names", $new_name);
	    }
	    
	    
	    ## Attribute "locus_tag"
	    my @locus_tags = $feature->get_attribute("locus_tag");
	    foreach my $locus_tag (@locus_tags) {
		$new_name = &trim($locus_tag);
		warn join( "\t", "Adding locus tag as name to feature", 
			   $feature->get_attribute("id"), 
			   $new_name
			   ), "\n" if ($verbose >= 4);
		$feature->push_attribute("names", $new_name);
	    }
	    
	    ## Some types of notes
	    my @notes = $feature->get_attribute("note");
	    foreach my $note (@notes) {
		$note = &trim($note); ### remove leading and trailing spaces
		warn "NOTE\t'$note'\n" if ($verbose >= 4);
		
		if ($note =~ /synonyms:/) {
		    ## For some genomes (e.g. Saccharomyces cerevisiae),
		    ## there is a note of type 'synonyms'
		    my @synonyms = split ',', $';
		    foreach my $new_name (@synonyms) {
			$new_name = &trim($new_name);
			warn join( "\t", "Adding synonym to feature", 
				   $feature->get_attribute("id"), 
				   $new_name
				   ), "\n" if ($verbose >= 4);
			$feature->push_attribute("names", $new_name);
		    }
		    
		} elsif ($note !~ /\S/) {
		    ## A single-word note is usually (but not always, I guess) a synonym
		    $feature->push_attribute("names", $note);
		    
		} elsif (
#		($note =~ /^locus_tag\: ([\w_\-\.]+)/i)  ||
			 ($note =~ /^locus_tag\: (\S+)/i)  ||
#		($note =~ /^locus_tag\: ([a-z0-9_\-\.]+)/i)  ||
			 ($note =~ /^Accession ([\w_\-\.]+)/) ||
			 ($note =~ /^(\S+)\,\s+len:/)
			 ){
		    
		    $new_name = $1;
		    
		    #### add the new name
		    warn join( "\t", "Adding name to feature", 
			       $feature->get_attribute("id"), 
			       $new_name
			       ), "\n" if ($verbose >= 4);
		    $feature->push_attribute("names", $new_name);
		}
		
	    }	
	    
	    warn join ("\t", "names done",
		       $class_holder->get_object_type(),
		       $feature->get_attribute("id"),
		       $feature->get_attribute('gene_id'),
		       $feature->get_attribute("type"),
		       join (":", $feature->get_attribute("names")),
		       ), "\n" if ($main::verbose >= 10);
	}
	
    }
}

################################################################

=pod

=item  CreateGenbankFeatures()

After having parsed genes, CDS and RNA from Genbank, create one
feature for each CDS and another one for each RNA. This feature
summarizes information (only retins selected fields) and reformats it
(multiple gene names) for ease of manipulation.

=cut
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
	if ($parsed_feature->get_attribute("gene_id")) {
	    ## Primary gene name is the one documented as "gene" attribute in the feature iself
	    $gene_name = $parsed_feature->get_attribute("gene");
	    if ($gene_name) {
		$created_feature->push_attribute("names",$gene_name);
	    }
	    
	    ## Identify the parent gene
	    $gene_id = $parsed_feature->get_attribute("gene_id");
	    next unless ($gene_id);
	    next if ($gene_id eq $main::null);
	    $gene = $genes->get_object($gene_id);


	    ## Add parent gene names to the current feature
	    if ($gene) {
		warn join ("\t", "feature parent gene", $created_feature->get_attribute("id"),
			   "gene name", $gene_name,
			   "gene_id", $gene_id,
			   "gene ID", $gene->get_attribute("id"),
			   ), "\n" if ($verbose >= 5);
		#### use gene name as primary name for the feature
		$created_feature->set_attribute("name", $gene->get_name());
		$created_feature->push_attribute("names", $gene->get_name());
		
	    	#### accept all gene names for the new feature
		foreach my $name ($gene->get_attribute("names")) {
		    $created_feature->push_attribute("names",$name) unless ($names{$name});
		}
	    } else {
		&ErrorMessage( "Cannot identify gene $gene_id for feature ", $parsed_feature, "\n");
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
	if (($products[0]) && ($products[0] ne $main::null)){
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
		if (($protein_id) && ($protein_id ne $main::null)){
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
		if (($transcript_id) && ($transcript_id ne $main::null)) {
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

    }

    ## Check object for all the parsed features, before building the RSAT features from it
    &CheckObjectNames($features, $genes, $mRNAs, $tRNAs, $rRNAs, $misc_RNAs, $misc_features, $CDSs);

}

################################################################

=pod

=item  CheckObjectNames()

Try to extract as many names as possible, make sure that each obect
has at least one name, and that there are no duplicate names for the
same object.

=cut
sub CheckObjectNames {
    my ($feature, $genes, $mRNAs, $tRNAs, $rRNAs, $misc_RNAs, $misc_features, $CDSs) = @_;


    ## Each feature inherits the names of its parent gene
    warn "; Treating object names\n" if ($main::verbose >= 10);
    foreach my $object (
			$genes->get_objects(),
			$mRNAs->get_objects(),
			$tRNAs->get_objects(),
			$rRNAs->get_objects(),
			$misc_RNAs->get_objects(),
			$misc_features->get_objects(),
			$CDSs->get_objects(),
			$features->get_objects(),
			) {

	## Make sure that the object has no null name
	my $primary_given = 0;
	my $primary_name = $object->get_name();
	if (($primary_name) && ($primary_name ne $main::null)) {
	    $primary_given = 1;
	}

	## The feature attribute "gene" is used as name
	my $gene_attr = $object->get_attribute("gene");
	if ($gene_attr) {
	    #### add gene name to the list of synonyms
	    $object->push_attribute("names", $gene_attr);

	    unless ($primary_given) {
		#### define this name as primary name
		$object->set_attribute("name",  $gene_attr);
		$primary_given = 1;
	    }
	}

	## The feature attribute "locus_tag" is used as name and
	## becomes primary name if there is no attribute "gene".
	my @locus_tags = $object->get_attribute("locus_tag");
	foreach my $locus_tag_attr (@locus_tags) {
	    #### add locus_tag name to the list of synonyms
	    $object->push_attribute("names", $locus_tag_attr);
	    
	    unless ($primary_given) {
		#### define this name as primary name
		$object->set_attribute("name",  $locus_tag_attr);
		$primary_given = 1;
	    }
	}

	################################################################
	#### Add all the names of the parent gene to the current feature
	my $gene_id = $object->get_attribute("gene_id");
	if ($gene_id) {
	    my $gene = $genes->get_object($gene_id);
	    if ($gene) {
		foreach my $gene_name ($gene->get_attribute("names")) {
		    $object->push_attribute("names", $gene_name);
		}
	    } else {
		&ErrorMessage("Cannot identify gene $gene_id for object ", $object, "\n");
	    }
	}

	unless ($primary_given) {
	    #### use ID as primary name
	    $object->set_attribute("name",  $object->get_attribute("id"));
	    $primary_given = 1;
	}

	## Make sure that object has no duplicate names
	$object->unique_names();
    }
    
}


1;


__END__

=pod

=back
