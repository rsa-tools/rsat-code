###############################################################
#
# Manipulation of features
#
package RSAT::feature;

%strand_index = (D=>0,
		 R=>1,
		 DR=>2
		 );

%default = (strand =>"DR",
	    frame =>".",
	    score =>0
	    );

# RSAT feature-map format (ft)
@{$columns{ft}} = qw (
		       seq_name
		       ft_type
                       feature_name
		       strand
		       start
		       end
		       description
                       score
		       );
@{$strands{ft}} = ("D", "R", "DR");
$comment_char{ft} = "; ";

# RSAT dna-pattern format (dnapat)
@{$columns{dnapat}} = qw (
			  feature_name
			  strand
			  pattern_sequence
			  seq_name
			  start
			  end
			  description
			  score
		       );
@{$strands{dnapat}} = ("D", "R", "DR");
$comment_char{dnapat} = "; ";

# RSAT Genome feature format (gft)
@{$columns{gft}} = qw (ft_id
		       ft_type
		       feature_name
		       seq_name
		       start
		       end
		       strand
		       description
		       );
@{$strands{gft}} = ("D", "R", "DR");
$comment_char{gft} = "; ";

## Sanger generic feature format (gff)
## http://www.sanger.ac.uk/Software/formats/GFF/GFF_Spec.shtml
@{$columns{gff}} = qw (seq_name
		       source
		       ft_type
		       start
		       end
		       score
		       strand
		       frame
		       attribute
		       );
@{$strands{gff}} = ("+", "-", ".");
$comment_char{gff} = "## ";

## Generic Feature Format version 3 (gff3)
## http://flybase.net/annot/gff3.html
@{$columns{gff3}} = qw (seq_name
			source
			ft_type
			start
			end
			score
			strand
			frame
			attribute
		       );
@{$strands{gff3}} = ("+", "-", ".");
$comment_char{gff3} = "## ";
@gff3_attributes = ("ID", #    name of the feature
		    "Name", #	 display name for the feature
		    "Alias", #	 secondary name for the feature
		    "Parent", #parent of the feature
		    "Target", #target (of an alignment)
		    "Gap", #   alignment of the feature to the target
		    "Note", #  A free text note
		    "Dbxref", # database cross reference
		    "Ontology_term"	# cross-reference to an ontology term
		   );

## UCSC BED
## http://genome.ucsc.edu/goldenPath/help/customTrack.html#BED
@{$columns{bed}} = qw (seq_name
		       start
		       end
		       feature_name
		       score
		       strand
		       thickStart
		       thickEnd
		       itemRGB
		       blockCount
		       blockSizes
		       blckStarts
		      );
@{$strands{bed}} = ("+", "-", ".");
$comment_char{bed} = "## ";

require "RSA.seq.lib";

use RSAT::GenericObject;
@ISA = qw( RSAT::GenericObject );

=pod

=head1 NAME

    RSAT::feature

=head1 DESCRIPTION

Object for manipulating features.

=head1 OUTPUT FORMATS

Feature formats are generally tab-delimited files. For historical
reasons, different formats have been used at different
sites. Originally, these formats were conceived to represent different
types of information, and, for this reason, the contents of their
column slightly differs. However, the main information is similar, and
the differing columns can be easily extrapolated or skipped.

=head2 ft

RSAT feature format for the program feature-map.

=over

=item col 1: map label (eg gene name)

=item col 2: feature type

=item col 3: feature identifier (ex: GATAbox, Abf1_site)

=item col 4: strand (D for Direct, R for Reverse),

=item col 5: feature start position

=item col 6: feature end position

=item col 7: (optional) description 

=item col 8: (optional) score

=back

=head2 gft

RSAT Genome features (file Feature.tab in the directory data/genomes).q

=over 

=item column 1: id

=item column 2:	type

=item column 3: name

=item column 4: contig

=item column 5: start

=item column 6: end

=item column 7: strand

=item column 8: description

=item column 9: chrom_position

=item column 10: organism


=back


=head2  gff: Sanger general feature file (extension .gff)

Information: http://www.sanger.ac.uk/Software/formats/GFF/GFF_Spec.shtml

=over

=item column 1: seqname (contig or sequence ID)

=item column 2: source

=item column 3: feature (the deature type name)

=item column 4: start

=item column 5: end

=item column 6: score

=item column 7: strand (+, - or .)

=item column 8: frame (0, 1, 2 or .)

=item column 9: attribute

=back 


=head3 Example of gff1:

 SEQ1	EMBL	atg	103	105	.	+	0
 SEQ1	EMBL	exon	103	172	.	+	0
 SEQ1	EMBL	splice5	172	173	.	+	.
 SEQ1	netgene	splice5	172	173	0.94	+	.
 SEQ1	genie	sp5-20	163	182	2.3	+	.
 SEQ1	genie	sp5-10	168	177	2.1	+	.
 SEQ2	grail	ATG	17	19	2.1	-	0

=head3 Example of gff2:

 seq1     BLASTX  similarity   101  235 87.1 + 0   Target "HBA_HUMAN" 11 55 ; E_value 0.0003
 dJ102G20 GD_mRNA coding_exon 7105 7201   .  - 2   Sequence "dJ102G20.C1.1"

=head2  gff3: Generic Feature Format version 3

Information: http://flybase.net/annot/gff3.html

=over

=item column 1: seqname (contig or sequence ID)

=item column 2: source

=item column 3: feature (the deature type name)

=item column 4: start

=item column 5: end

=item column 6: score

=item column 7: strand (+, - or .)

=item column 8: frame (0, 1, 2 or .)

=item column 9: attribute

 ID:     name of the feature
 Name: 	 display name for the feature
 Alias:	 secondary name for the feature
 Parent: parent of the feature
 Target: target (of an alignment)
 Gap:    alignment of the feature to the target
 Note:   A free text note
 Dbxref: database cross reference
 Ontology_term:	cross-reference to an ontology term

=back

=head3 Example of gff3

 ##gff-version 3 
 ##sequence-region ctg123 1 1497228 
 ctg123	.	gene	1000	9000	.	+	.	ID=gene00001;Name=EDEN	
 ctg123	.	TF_binding_site	1000	1012	.	+	.	ID=tfbs00001;Parent=gene00001	
 ctg123	.	mRNA	1050	9000	.	+	.	ID=mRNA00001;Parent=gene00001;Name=EDEN.1	
 ctg123	.	mRNA	1050	9000	.	+	.	ID=mRNA00002;Parent=gene00001;Name=EDEN.2
 ctg123	.	mRNA	1300	9000	.	+	.	ID=mRNA00003;Parent=gene00001;Name=EDEN.3	
 ctg123	.	exon	1300	1500	.	+	.	ID=exon00001;Parent=mRNA00003	
 ctg123	.	exon	1050	1500	.	+	.	ID=exon00002;Parent=mRNA00001,mRNA00002	
 ctg123	.	exon	3000	3902	.	+	.	ID=exon00003;Parent=mRNA00001,mRNA00003	
 ctg123	.	exon	5000	5500	.	+	.	ID=exon00004;Parent=mRNA00001,mRNA00002,mRNA00003	
 ctg123	.	exon	7000	9000	.	+	.	ID=exon00005;Parent=mRNA00001,mRNA00002,mRNA00003	
 ctg123	.	CDS	1201	1500	.	+	0	ID=cds00001;Parent=mRNA00001;Name=edenprotein.1	
 ctg123	.	CDS	3000	3902	.	+	0	ID=cds00001;Parent=mRNA00001;Name=edenprotein.1	
 ctg123	.	CDS	5000	5500	.	+	0	ID=cds00001;Parent=mRNA00001;Name=edenprotein.1	
 ctg123	.	CDS	7000	7600	.	+	0	ID=cds00001;Parent=mRNA00001;Name=edenprotein.1	
 ctg123	.	CDS	1201	1500	.	+	0	ID=cds00002;Parent=mRNA00002;Name=edenprotein.2	
 ctg123	.	CDS	5000	5500	.	+	0	ID=cds00002;Parent=mRNA00002;Name=edenprotein.2	
 ctg123	.	CDS	7000	7600	.	+	0	ID=cds00002;Parent=mRNA00002;Name=edenprotein.2	
 ctg123	.	CDS	3301	3902	.	+	0	ID=cds00003;Parent=mRNA00003;Name=edenprotein.3	
 ctg123	.	CDS	5000	5500	.	+	2	ID=cds00003;Parent=mRNA00003;Name=edenprotein.3	
 ctg123	.	CDS	7000	7600	.	+	2	ID=cds00003;Parent=mRNA00003;Name=edenprotein.3	
 ctg123	.	CDS	3391	3902	.	+	0	ID=cds00004;Parent=mRNA00003;Name=edenprotein.4	
 ctg123	.	CDS	5000	5500	.	+	2	ID=cds00004;Parent=mRNA00003;Name=edenprotein.4	
 ctg123	.	CDS	7000	7600	.	+	2 	ID=cds00004;Parent=mRNA00003;Name=edenprotein.4

=head2 bed

Genomic features in the UCSC format.  Beware, this format assumes that
features are described with chromosomal positions.

Information: http://genome.ucsc.edu/goldenPath/help/hgTracksHelp.html#BED

=over

=item 1. chrom

The name of chromosome (e.g. chr3, chrY, chr2_random) or scaffold
(e.g. scaffold10671).

=item 2. chromStart

The starting position of the feature in the chromosome or
scaffold. The first base in a chromosome is numbered 0.

=item 3. chromEnd

The ending position of the feature in the chromosome or scaffold. The
chromEnd base is not included in the display of the feature. For
example, the first 100 bases of a chromosome are defined as
chromStart=0, chromEnd=100, and span the bases numbered 0-99.

=item 4. name

Defines the name of the BED line. This label is displayed to the left
of the BED line in the Genome Browser window when the track is open to
full display mode or directly to the left of the item in pack mode.

=item 5. score

A score between 0 and 1000.

=item 6. strand

Defines the strand - either '+' or '-'.

=item 7. thickStart

The starting position at which the feature is drawn thickly (for
example, the start codon in gene displays).

=item 8. thickEnd

The ending position at which the feature is drawn thickly (for
example, the stop codon in gene displays).

   
=item 9. itemRgb

An RGB value of the form R,G,B (e.g. 255,0,0). If the track line
itemRgb attribute is set to "On", this RBG value will determine the
display color of the data contained in this BED line. NOTE: It is
recommended that a simple color scheme (eight colors or less) be used
with this attribute to avoid overwhelming the color resources of the
Genome Browser and your Internet browser.

=item 10. blockCount

The number of blocks (exons) in the BED line.

=item 11. blockSizes

A comma-separated list of the block sizes. The number of items in this
list should correspond to blockCount.

=item 12. blockStarts

A comma-separated list of block starts. All of the blockStart
positions should be calculated relative to chromStart. The number of
items in this list should correspond to blockCount.

=back

=head3 Example of bed format

 browser position chr7:127471196-127495720
 browser hide all
 track name="ItemRGBDemo" description="Item RGB demonstration" visibility=2 
 itemRgb="On" 
 chr7	127471196  127472363  Pos1  0  +  127471196  127472363  255,0,0
 chr7	127472363  127473530  Pos2  0  +  127472363  127473530  255,0,0
 chr7	127473530  127474697  Pos3  0  +  127473530  127474697  255,0,0
 chr7	127474697  127475864  Pos4  0  +  127474697  127475864  255,0,0
 chr7	127475864  127477031  Neg1  0  -  127475864  127477031  0,0,255
 chr7	127477031  127478198  Neg2  0  -  127477031  127478198  0,0,255
 chr7	127478198  127479365  Neg3  0  -  127478198  127479365  0,0,255
 chr7	127479365  127480532  Pos5  0  +  127479365  127480532  255,0,0
 chr7	127480532  127481699  Neg4  0  -  127480532  127481699  0,0,255


=head1 METHODS

=cut


################################################################
=pod

=over

=item new()

Create a new feature.

=cut
sub new {
    my ($class, %args) = @_;
    my $feature = bless {
	}, $class;
    return $feature;
}

################################################################
=pod

=item parse_one_row($row, $format)

Parse the feature from a text row.

=cut
sub parse_from_row {
  my ($self, $row, $format) = @_;
  chomp($row);
  $row =~ s/\r//g;
  my @fields = split("\t", $row);
  warn join( "\t", "parsing from ", $format, @fields), "\n" if ($main::verbose >= 10);

  my @cols = @{$columns{$format}};
  foreach my $c (0..$#cols) {
    my $attr = $cols[$c];
    if (defined($fields[$c])) {
      my $value = $fields[$c];
      &RSAT::message::Debug("column ".($c+1), "attr:".$attr, "value:".$value) if ($main::verbose >= 10);
      $self->set_attribute($attr, $value);
    } else {
      &RSAT::message::Warning("Missing attribute ".$attr, "column:".($c+1)) if ($main::verbose >= 3);
    }
  }

  ## Convert strand format
  my $strand = $self->get_attribute("strand");
#  &RSAT::message::Debug("strand", $strand) if ($main::verbose >= 0);
  if ($strand) {
    $strand =~ s/\+/D/;
    $strand =~ s/\-/R/;
    $strand =~ s/\./D/;
  }
  $self->force_attribute("strand", $strand);

  ## Format-specific conversions
  if (($format eq "gff") || ($format eq "gff3")) {
    $self->set_attribute("feature_name", $self->get_attribute("source"));
    $self->set_attribute("description", $self->get_attribute("attribute"));
  }

  ## parse attributes from the attribute/description field
  my $description = "";
  if (($format eq "gff") || ($format eq "gff3")) {
    $description = $self->get_attribute("attribute");
    ## Convert single-note attribute into description (suppress "Note=" from the beginning)
    if ($description =~ /^Note=([^;]+)$/) {
      $self->force_attribute("description",  $1);
#      die "HELLO";
    }

  } else {
    $description = $self->get_attribute("description");
  }

  ## Parse attributes from the description field
  if ($description) {

    ## Parse attributes from gff/gff3 format
    my @attributes = split /; */, $description;
	&RSAT::message::Debug( $description, join(";", @attributes)) if ($main::verbose >= 3);
    if (scalar(@attributes) > 0) {
      my $gff_attributes_found = 1;
      foreach my $a (0..$#attributes) {
	my $attribute = $attributes[$a];
	if (($attribute =~ /(\S+) +(\S.+)/) ||
	    ($attribute =~ /(\S+)\s*=\s*(\S+)/)) {
	  $gff_attributes_found = 1;
	  my $attr = $1;
	  my $value = $2;
	  $value =~ s/^\"//; 
	  $value =~ s/\"$//;
	  $self->force_attribute($attr, $value);
	  if (lc($attr) eq "id") {
	    $self->force_attribute("ft_id", $value);
	    ## Use ID as name unless name has already been defined
	    #		      die "HELLO";
	    $self->force_attribute("feature_name", $value);
	    ##			$self->force_attribute("id", $value);
	  }
	  if (lc($attr) eq "name") {
	    $self->force_attribute("feature_name", $value);
	    ##			$self->force_attribute("name", $value);
	  }
	  if (lc($attr) eq "site") {
	    $self->force_attribute("pattern_sequence", $value);
	  }
	  &RSAT::message::Debug("Feature",$self->get_attribute("ID"), 
				"Parsed attribute", $a."/".$#attributes, $attr, $value) 
	    if ($main::verbose >= 5);
	}
      }

       ## convert description into gff note
#       unless ($gff_attribute_found) {
# 	$self->force_attribute("Note", $description);
#       }
      #    } else {
    }

  }

  ## dna-pattern
  if ($format eq "dnapat") {
    $self->force_attribute("type", "dnapat");
  }
  &RSAT::message::Info(join("\t",
			    "Parsed new feature",
			    $self->get_attribute("seq_name"),
			    $self->get_attribute("feature_type"),
			    $self->get_attribute("feature_name"),
			    $self->get_attribute("id"),
			    $self->get_attribute("start"),
			    $self->get_attribute("end"),
			    $self->get_attribute("strand"),
			    $self->get_attribute("description"),
			    $self->get_attribute("score"))
		      ) if ($main::verbose >= 3);
}

################################################################
=pod

=item to_text($format)

Converts the feature in a single-row string for exporting it in the
specified format.

=cut

sub to_fasta {
  my ($self) = @_;
  my $feature_id = join( ":", 
			 $self->get_attribute("seq_name"),
			 $self->get_attribute("start"),
			 $self->get_attribute("end"),
			 $self->get_attribute("strand"),
		       );
  my $fasta_string = ">".$feature_id."\n";
  $fasta_string .= $self->get_attribute("description");
  $fasta_string .= "\n";
#  &RSAT::message::Debug($fasta_string);
  return($fasta_string);
}

################################################################
=pod

=item to_text($format)

Converts the feature in a single-row string for exporting it in the
specified format.

=cut
sub to_text {
    my ($self, $format, $null) = @_;
    $null = "" unless (defined($null));


    if ($format eq "fasta") {
      return $self->to_fasta();
    }

    my @cols = @{$columns{$format}};


    ## Index column number by contents
    my ${col_index} = ();
    foreach my $c (0..$#cols) {
	$col_index{$cols[$c]} = $c;
    }

    ## Select the fields
    my @fields = ();
    foreach my $c (0..$#cols) {
      $attr = $cols[$c];
      $field_value = $self->get_attribute($attr);
      unless ($field_value) {
	if (defined($default{$attr})) {
	  $field_value = $default{$attr};
	} elsif ($attr eq "source") {
	  $field_value = $main::input_format;
	} else {   
	  $field_value = $null;
	}
      }
      $fields[$c] =  $field_value;
      &RSAT::message::Debug("field", $c, $fields[$c])if ($main::verbose >= 10);
    }

    ################################################################
    ## Format-specific attributes

    ## Collect attributes for gff and gff3 formats
    if ($format =~ /gff/) {
#       unless ($self->get_attribute("gene")) {
# 	if ($self->get_attribute("feature_name")) {
# 	  $self->set_attribute("gene", $self->get_attribute("feature_name"));
# 	}
#       }

      my @attributes = ();
      foreach $attr (@gff3_attributes) {
	my $value = "";
	if (defined($self->get_attribute($attr))) {
	  $value = $self->get_attribute($attr);
	  push @attributes, $attr."=".$value;
	  &RSAT::message::Debug("export attribute", $attr, $value) if ($main::verbose >= 3);
	} else {
	  &RSAT::message::Debug("feature", $self->get_attribute("ID"), "undefined gff3 attribute", $attr) if ($main::verbose >= 3);
	}
      }

      ## If no info has been found, use the description field as note
      if (scalar(@attributes) == 0) {
	my $description = $self->get_attribute("description");
	if ($description) {
	  if ($format =~ /gff/) {
	    push @attributes, "Note=".$description;
	  }
	}
      }

      ## Concatenate the attributes in a string
      my $attribute = join ";", @attributes;
      $fields[$col_index{attribute}] = $attribute;
    }

    ## Format-specific treatment for the strand
    my @strands = @{$strands{$format}};
    my $strand = $self->get_attribute("strand") || $default{strand};
    my $f = $col_index{"strand"};
    if ($strand) {
	$s = $strand_index{$strand};
    }
    $fields[$f] = $strands[$s];
#    &RSAT::message::Debug( "strand", $strand, "f=$f", 
#			   "index:".join(";", %strand_index),
#			   $strand_index{$strand},
#			   "format:".join(";", @strands), 
#			   "s=".$s, $strands[$s], "field=$fields[$f]") if ($main::verbose >= 10);

    ## Generate the row to be printed
    my $row = join ("\t", @fields);
    $row .= "\n";

#    &RSAT::message::Debig ("printing in format ", $format, $row) if ($main::verbose >= 10);

    return($row);
}

################################################################
=pod

=item header($format)

Print the header in the specified format.

=cut

sub header {
    my ($format) = @_;
    if ($format eq "fasta") {
      return();
    }

    ## Print format
    my $header = $comment_char{$format};
    if ($format eq "gff3") {
      $header .= "gff-version\t3";
    } else {
      $header .= $format;
    }
    $header .= "\n";

    ## Print column content
    my @cols = @{$columns{$format}};
    $header .= $comment_char{$format};
    $header .= join ("\t", @cols);
    $header .= "\n";
    return $header;
}


################################################################
=pod

=item full_id

Return a full identifier for the feature

=cut

sub full_id {
    my ($self)  = @_;
    return (join (":", 
		  $self->get_attribute("filename"),
		  $self->get_attribute("id"),
		  $self->get_attribute("feature_name")
		 ));

}


return 1;


__END__

=pod

=back

