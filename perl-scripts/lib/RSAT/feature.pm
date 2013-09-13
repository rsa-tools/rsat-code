###############################################################
#
# Manipulation of features
#
package RSAT::feature;

%supported_input_format =(ft=>1,
			  gft=>1,
			  gff=>1,
			  gff3=>1,
			  dnapat=>1,
			  bed=>1,
			  galaxy_seq=>1,
			  ucsc_seq=>1,
			  swembl=>1,
			 );
$supported_input_formats = join (",", keys %supported_input_format);

%supported_output_format =(ft=>1,
			   fasta=>1,
			   gft=>1,
			   gff=>1,
			   gff3=>1,
			   dnapat=>1,
			   bed=>1,
			  );
$supported_output_formats = join (",", keys %supported_output_format);

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
$header_char{ft} = "# ";

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
$header_char{dnapat} = "# ";

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
$header_char{gft} = "# ";

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
$header_char{gff} = "## ";

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
$header_char{gff3} = "## ";
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
## warning: UCSC BED files should be zero-based
@{$columns{bed}} = qw (seq_name
		       start
		       end
		       feature_name
		       score
		       strand
		       thickStart
		       thickEnd
		       itemRgb
		       blockCount
		       blockSizes
		       blckStarts
		      );
@{$strands{bed}} = ("+", "-", ".");
$comment_char{bed} = "## ";
$header_char{bed} = "## ";

## SWEMBL
## http://www.ebi.ac.uk/~swilder/SWEMBL/
## SWEMBL exports a bed-like format for the 3 first columns, and custom info in the other columns
@{$columns{swembl}} = qw (seq_name
		       start
		       end
		       count
		       length
		       unique_pos
		       score
		       ref_count
		       max_coverage
		       summit
		      );
@{$strands{swembl}} = ("+", "-", ".");
$comment_char{swembl} = "#";
$header_char{swembl} = "#";
#$header_char{swembl} = "Region";


## Galaxy sequences
## http://main.g2.bx.psu.edu/tool_runner?tool_id=Extract+genomic+DNA+1
## warning: the coordinates are zero-based
@{$columns{galaxy_seq}} = qw (assembly
			      seq_name
			      start
			      end
			      strand
			     );
@{$strands{galaxy_seq}} = ("+", "-");
$comment_char{galaxy_seq} = "#";
$header_char{galaxy_seq} = "#";


## UCSC sequences
## warning: the coordinates are zero-based
@{$columns{ucsc_seq}} = qw (seq_name
			    start
			    end
			    strand
			   );
@{$strands{ucsc_seq}} = ("1", "2");

## Define specific formats
our %format = ();
$format{"start"} = '%d';
$format{"end"} = '%d';


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


=item I<column 1: seqname>

The name of the sequence. Having an explicit sequence name allows a
feature file to be prepared for a data set of multiple
sequences. Normally the seqname will be the identifier of the sequence
in an accompanying fasta format file. An alternative is that <seqname>
is the identifier for a sequence in a public database, such as an
EMBL/Genbank/DDBJ accession number. Which is the case, and which file
or database to use, should be explained in accompanying information.

=item I<column 2: source>

The source of this feature. This field will normally be used to
indicate the program making the prediction, or if it comes from public
database annotation, or is experimentally verified, etc.

=item I<column 3: feature>

The feature type name. We hope to suggest a standard set of features,
to facilitate import/export, comparison etc.. Of course, people are
free to define new ones as needed. For example, Genie splice detectors
account for a region of DNA, and multiple detectors may be available
for the same site, as shown above.

We would like to enforce a standard nomenclature for common GFF
features. This does not forbid the use of other features, rather, just
that if the feature is obviously described in the standard list, that
the standard label should be used. For this standard table we propose
to fall back on the international public standards for genomic
database feature annotation, specifically, the DDBJ/EMBL/GenBank
feature table documentation).

=item I<column 4: start>

=item I<column 5: end>

Integers. <start> must be less than or equal to <end>. Sequence
numbering starts at 1, so these numbers should be between 1 and the
length of the relevant sequence, inclusive. (Version 2 change: version
2 condones values of <start> and <end> that extend outside the
reference sequence. This is often more natural when dumping from
acedb, rather than clipping. It means that some software using the
files may need to clip for itself.)

=item I<column 6: score>

A floating point value. When there is no score (i.e. for a sensor that
just records the possible presence of a signal, as for the EMBL
features above) you should use '.'. (Version 2 change: in version 1 of
GFF you had to write 0 in such circumstances.)

=item I<column 7: strand>

One of '+', '-' or '.'. 

'.' should be used when strand is not relevant, e.g. for dinucleotide
repeats. Version 2 change: This field is left empty '.' for RNA and
protein features.

=item I<column 8: frame>

One of '0', '1', '2' or '.'. 

'0' indicates that the specified region is in frame, i.e. that its
first base corresponds to the first base of a codon. '1' indicates
that there is one extra base, i.e. that the second base of the region
corresponds to the first base of a codon, and '2' means that the third
base of the region is the first base of a codon. If the strand is '-',
then the first base of the region is value of <end>, because the
corresponding coding region will run from <end> to <start> on the
reverse strand. As with <strand>, if the frame is not relevant then
set <frame> to '.'. It has been pointed out that "phase" might be a
better descriptor than "frame" for this field. Version 2 change: This
field is left empty '.' for RNA and protein features.

=item I<column 9: attribute>

From version 2 onwards, the attribute field must have an tag value
structure following the syntax used within objects in a .ace file,
flattened onto one line by semicolon separators. Tags must be standard
identifiers ([A-Za-z][A-Za-z0-9_]*). Free text values must be quoted
with double quotes. Note: all non-printing characters in such free
text value strings (e.g. newlines, tabs, control characters, etc) must
be explicitly represented by their C (UNIX) style backslash-escaped
representation (e.g. newlines as '\n', tabs as '\t'). As in ACEDB,
multiple values can follow a specific tag. The aim is to establish
consistent use of particular tags, corresponding to an underlying
implied ACEDB model if you want to think that way (but acedb is not
required). Examples of these would be:

 seq1     BLASTX  similarity   101  235 87.1 + 0	Target "HBA_HUMAN" 11 55 ; E_value 0.0003
 dJ102G20 GD_mRNA coding_exon 7105 7201   .  - 2 Sequence "dJ102G20.C1.1"

The semantics of tags in attribute field tag-values pairs has
intentionally not been formalized. Two useful guidelines are to use
DDBJ/EMBL/GenBank feature 'qualifiers' (see DDBJ/EMBL/GenBank feature
table documentation), or the features that ACEDB generates when it
dumps GFF.

Version 1 note In version 1 the attribute field was called the group
field, with the following specification: An optional string-valued
field that can be used as a name to group together a set of
records. Typical uses might be to group the introns and exons in one
gene prediction (or experimentally verified gene structure), or to
group multiple regions of match to another sequence, such as an EST or
a protein.

=back


B<Description of gfffrom UCSC genome browser help page>

http://genome.ucsc.edu/goldenPath/help/customTrack.html#GFF

=over

=item column 1: seqname (contig or sequence ID)

The name of the sequence. Must be a chromosome or scaffold.

=item column 2: source

The program that generated this feature.

=item column 3: feature

The name of this type of feature. Some examples of standard feature
types are "CDS", "start_codon", "stop_codon", and "exon".


=item column 4: start

The starting position of the feature in the sequence. The first base
is numbered 1.


=item column 5: end

The ending position of the feature (inclusive).

=item column 6: score

A score between 0 and 1000. If the track line useScore attribute is
set to 1 for this annotation data set, the score value will determine
the level of gray in which this feature is displayed (higher numbers =
darker gray). If there is no score value, enter ".".


=item column 7: strand (+, - or .)

Valid entries include '+', '-', or '.' (for don't know/don't care).


=item column 8: frame (0, 1, 2 or .)

If the feature is a coding exon, frame should be a number between 0-2
that represents the reading frame of the first base. If the feature is
not a coding exon, the value should be '.'.

=item column 9: attribute

All lines with the same group are linked together into a single item.

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

Genomic features in the UCSC format.  Warning: this format assumes that
features are described with chromosomal positions, and should be zero-based (meaning that
the first position is 0, not 1).

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

=head2 galaxy_seq

Fasta sequences retrieved from the website Galaxy
(http://main.g2.bx.psu.edu/tool_runner?tool_id=Extract+genomic+DNA+1)
Warning: this format assumes that features are described with
chromosomal positions, and should be zero-based (meaning that the
first position is 0, not 1).

>hg17_chr7_127475281_127475310_+

=over

=item 1. assembly

UCSC-style description of the assembly (e.g. hg19,mm9)

=item 2. chromosome

The name of chromosome (e.g. chr3, chrY, chr2_random) or scaffold
(e.g. scaffold10671).

=item 3. chromStart

The starting position of the feature in the chromosome or
scaffold. The first base in a chromosome is numbered 0.

=item 4. chromEnd

The ending position of the feature in the chromosome or scaffold. The
chromEnd base is not included in the display of the feature. For
example, the first 100 bases of a chromosome are defined as
chromStart=0, chromEnd=100, and span the bases numbered 0-99.

=item 5. strand

Defines the strand - either '+' or '-'.

=back

=head2 ucsc_seq

Fasta sequences retrieved from the UCSC server.  Warning: this format
assumes that features are described with chromosomal positions, and
should be zero-based (meaning that the first position is 0, not 1).

Coordinates are parsed from the fasta header.

>11:3052022..3052331:1


=head3 Example of galaxy_seq format

>hg17_chr7_127475281_127475310_+
GTAGGAATCGCAGCGCCAGCGGTTGCAAG
>hg17_chr7_127485994_127486166_+
GCCCAAGAAGCCCATCCTGGGAAGGAAAATGCATTGGGGAACCCTGTGCG
GATTCTTGTGGCTTTGGCCCTATCTTTTCTATGTCCAAGCTGTGCCCATC
CAAAAAGTCCAAGATGACACCAAAACCCTCATCAAGACAATTGTCACCAG
GATCAATGACATTTCACACACG
>hg17_chr7_127486011_127486166_+
TGGGAAGGAAAATGCATTGGGGAACCCTGTGCGGATTCTTGTGGCTTTGG
CCCTATCTTTTCTATGTCCAAGCTGTGCCCATCCAAAAAGTCCAAGATGA
CACCAAAACCCTCATCAAGACAATTGTCACCAGGATCAATGACATTTCAC
ACACG


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

=item parse_one_row($row, $in_format)

Parse the feature from a text row.

=cut
sub parse_from_row {
  my ($self, $row, $in_format, $out_format) = @_;
  chomp($row);
  $row =~ s/\r//g;

  ## Split the row into fields (tab-delimited columns)
  my @fields = ();
  if ($in_format eq "galaxy_seq"){
    $row =~ s/^\s*>//;

    ## PROBLEM HERE: DOES NOT WORK IF THE ID CONTAINS "_" characters
    #    @fields = split("_", $row);
    if ($row =~ /(\S+)*_(\S+)_(\d+)_(\d+)_([+-])$/) {
      @fields = ($1, $2, $3, $4, $5);
    } else {
      &RSAT::message::Warning("Invalid galaxy fasta header for feature extraction", $row) if ($main::verbose >= 0);
      return();
    }

  } elsif ($in_format eq "ucsc_seq") {
    if ($row =~ /^>(\S+):(\S+)\.\.(\S+)\:(\S+)/) {
#    if ($row =~ />11:3052022..3052331:1	/target='11:3052022..3052331' /seq_id='11' /strand='+' /type='region' /original_strand='+' /end='3052331' /start='3052022'/) {
      @fields = ($1, $2, $3, $4);
    } else {
      &RSAT::message::Warning("Invalid galaxy fasta header for feature extraction", $row) if ($main::verbose >= 0);
      return();
    }

  } else {
    @fields = split("\t", $row);
  }
  &RSAT::message::Debug("parsing from ", $in_format, @fields) if ($main::verbose >= 10);

  ## Identify attributes in columns
  my @cols = @{$columns{$in_format}};
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
#  &RSAT::message::Debug("strand", $strand) if ($main::verbose >= 10);
  if ($strand) {
    $strand =~ s/\+/D/;
    $strand =~ s/\-/R/;
    $strand =~ s/\./D/;
  } else {
    $strand = "DR";
  }
  $self->force_attribute("strand", $strand);

  ## Format-specific conversions
  if (($in_format eq "gff") || ($in_format eq "gff3")) {
    $self->set_attribute("feature_name", $self->get_attribute("source"));
    $self->set_attribute("description", $self->get_attribute("attribute"));

  } elsif ($in_format eq "swembl")  {
    ## SWEMBL occasionally returns peak with a negative start
    ## coordinate ! I guess this is due to the algorithm for defining
    ## the peak width, but it creates obvious problems with the
    ## programs used to further analyze the peaks. We circumvent this
    ## by replacing negative and null values by 1.
    if ($self->get_attribute("start") < 0) {
      $self->force_attribute("start", 1)
    }

    my $name = join("_",
		    $self->get_attribute("seq_name"),
		    $self->get_attribute("start"),
		    $self->get_attribute("end"),
		    "+"
		    );
    $self->set_attribute("feature_name", $name);
    $self->set_attribute("description", $self->get_attribute("feature_name"));

  } elsif ($in_format eq "galaxy_seq")  {
    $self->set_attribute("ft_type", "");
    $row =~ s/^\s*>//;
    $self->set_attribute("feature_name", $row);
  }

  ## Convert name
  if (defined($out_format) && ($out_format =~ /gff/)) {
    if (($self->get_attribute('feature_name')) && (!$self->get_attribute("Name"))) {
      $self->force_attribute("Name", $self->get_attribute("feature_name"));
    }
    if (($self->get_attribute("description")) &&
	($self->get_attribute("description") !~ /=/) &&
	(!$self->get_attribute("Note"))) {
      $self->force_attribute("Note", $self->get_attribute("description"));
    }
  }

  ## Parse attributes from the attribute/description field
  my $description = "";
  if (($in_format eq "gff") || ($in_format eq "gff3")) {
    $description = $self->get_attribute("attribute");
    ## Convert single-note attribute into description (suppress "Note=" from the beginning)
    if ($description =~ /^Note=([^;]+)$/) {
      $self->force_attribute("description",  $1);
    }

  } else {
    $description = $self->get_attribute("description");
  }

  ## Parse attributes from the description field
  if ($description) {

    ## Parse attributes from gff/gff3 format
    my @attributes = split /; */, $description;
	&RSAT::message::Debug( $description, join(";", @attributes)) if ($main::verbose >= 4);
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
  if ($in_format eq "dnapat") {
    $self->force_attribute("ft_type", "pattern");
  }

  ## bed format: convert start position from 0-based to 1-based coordinate
  ##
  ## See http://genome.ucsc.edu/FAQ/FAQformat.html#format1
  ##
  ##  chromStart - The starting position of the feature in the
  ##    chromosome or scaffold. The first base in a chromosome is
  ##    numbered 0.
  ##  chromEnd - The ending position of the feature in the
  ##     chromosome or scaffold. The chromEnd base is not included
  ##     in the display of the feature. For example, the first 100
  ##     bases of a chromosome are defined as chromStart=0,
  ##     chromEnd=100, and span the bases numbered 0-99.
  if (($in_format eq "bed") || ($in_format eq "galaxy_seq")) {
     my $start_0_based = $self->get_attribute("start");
     $self->force_attribute("start",$start_0_based+1); ## only the fist base is shifted
  }

  ## parsed row
  &RSAT::message::Info("Parsed new feature",
		       $self->get_attribute("seq_name"),
		       $self->get_attribute("ft_type"),
		       $self->get_attribute("feature_name"),
		       $self->get_attribute("id"),
		       $self->get_attribute("start"),
		       $self->get_attribute("end"),
		       $self->get_attribute("strand"),
		       $self->get_attribute("description"),
		       $self->get_attribute("score"),
		      ) if ($main::verbose >= 4);

  return();
}

################################################################

=pod

=item to_text($out_format)

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

=item to_text($out_format)

Converts the feature in a single-row string for exporting it in the
specified format.

=cut
sub to_text {
  my ($self, $out_format, $null) = @_;
  $null = "" unless (defined($null));

  ## For the BED format
  if ($out_format eq "bed") {

    ## Suppress sequence start and end features (temporary fix for
    ## UCSC genome browser)
    my $ft_type = $self->get_attribute("ft_type") || "feature";
    if ($ft_type eq "limit") {
      if ($self->get_attribute("feature_name") eq "START_END") {
	my $seq_name = $self->get_attribute("seq_name");
	my $seq_start = $self->get_attribute("start");
	my $seq_end = $self->get_attribute("end");
	my $string = "browser position ".${seq_name}.":".${seq_start}."-".${seq_end}."\n";
	return($string);
      } else {
	&RSAT::message::Warning("Skipping feature", $self->get_attribute("feature_name"), "for BED compatibility") if ($main::verbose >= 2);
	return();
      }
    }

    ## Transform the start position into zero-based coordinate.
    ##
    ## See http://genome.ucsc.edu/FAQ/FAQformat.html#format1
    ##
    ##  chromStart - The starting position of the feature in the
    ##    chromosome or scaffold. The first base in a chromosome is
    ##    numbered 0.
    ##  chromEnd - The ending position of the feature in the
    ##     chromosome or scaffold. The chromEnd base is not included
    ##     in the display of the feature. For example, the first 100
    ##     bases of a chromosome are defined as chromStart=0,
    ##     chromEnd=100, and span the bases numbered 0-99.
    my $start_1_based = $self->get_attribute("start");
    $self->force_attribute("start",$start_1_based-1); ## only the fist base is shifted
  }

  ## Treat the thickStart and thickEnd attributes
  unless (defined($self->get_attribute("thickStart"))) {
    $self->set_attribute("thickStart",$self->get_attribute("start"));
  }
  unless (defined($self->get_attribute("thickEnd"))) {
    $self->set_attribute("thickEnd",$self->get_attribute("end"));
  }

  ## Fasta format
  if ($out_format eq "fasta") {
    return $self->to_fasta();
  }


  ## Tab-delimited column files
  my @cols = @{$columns{$out_format}};


  ## Index column number by contents
  my ${col_index} = ();
  foreach my $c (0..$#cols) {
    $col_index{$cols[$c]} = $c;
  }

  ## Select the fields
  my @fields = ();
  foreach my $c (0..$#cols) {
    my $attr = $cols[$c];
    my $field_value = $self->get_attribute($attr);

    ## Check null attributes
    unless ($field_value) {
      if (defined($default{$attr})) {
	$field_value = $default{$attr};
      } elsif ($attr eq "source") {
	$field_value = $main::input_format;
      } elsif (!defined($format{$attr})) {
	$field_value = $null;
      }
    }

    ## Check attribute formats
    if (defined($format{$attr})) {
#      &RSAT::message::Warning("Formatting attriubte", $attr, $format{$attr}, $field_value);
      $field_value = sprintf $format{$attr}, $field_value;
    }
    $fields[$c] =  $field_value;
    #      &RSAT::message::Debug("field", $c, sprintf("%-15s", $attr), $fields[$c])if ($main::verbose >= 10);
  }

  ################################################################
  ## Format-specific attributes

  ## Collect attributes for gff and gff3 formats
  if ($out_format =~ /gff/) {
    #       unless ($self->get_attribute("gene")) {
    # 	if ($self->get_attribute("feature_name")) {
    # 	  $self->set_attribute("gene", $self->get_attribute("feature_name"));
    # 	}
    #       }

    my @attributes = ();
    my %attributes = ();	## Index for further tests
    foreach $attr (@gff3_attributes) {
      my $value = "";
      $attributes{$attr} =  $value; ## index for further tests
      if (defined($self->get_attribute($attr))) {
	$value = $self->get_attribute($attr);
	push @attributes, $attr."=".$value;
	&RSAT::message::Debug("export attribute", $attr, $value) if ($main::verbose >= 4);
      } else {
	&RSAT::message::Debug("feature", $self->get_attribute("ID"), "undefined gff3 attribute", $attr) if ($main::verbose >= 3);
      }
    }

    ## If no info has been found, use the description field as note
    if (scalar(@attributes) == 0) {
      my $description = $self->get_attribute("description");
      if ($description) {
	if ($out_format =~ /gff/) {
	  push @attributes, "Note=".$description;
	}
      }
    }

    ## Concatenate the attributes in a string
    my $attribute = join ";", @attributes;
    $fields[$col_index{attribute}] = $attribute;
  }

  ## Format-specific treatment for the strand
  my @strands = @{$strands{$out_format}};
  my $strand = $self->get_attribute("strand") || $default{strand};
  my $f = $col_index{"strand"};
  my $s;
  if ($strand) {
    if (defined($strand_index{$strand})) {
      $s = $strand_index{$strand};
    } else {
      $s = $strand_index{'DR'};
    }
  } else {
    $s = $strand_index{'DR'};
  }
#  &RSAT::message::Debug($f, $strand, $s, @strands, %strand_index) if ($main::verbose >= 10);
  $fields[$f] = $strands[$s];
  #    &RSAT::message::Debug( "strand", $strand, "f=$f", 
  #			   "index:".join(";", %strand_index),
  #			   $strand_index{$strand},
  #			   "format:".join(";", @strands), 
  #			   "s=".$s, $strands[$s], "field=$fields[$f]") if ($main::verbose >= 10);

  ## Generate the row to be printed
  my $row = join ("\t", @fields);
  $row .= "\n";

  #    &RSAT::message::Debig ("printing in format ", $out_format, $row) if ($main::verbose >= 10);

  return($row);
}

################################################################

=pod

=item header($out_format)

Print the header in the specified format.

=cut

sub header {
    my ($out_format) = @_;
    if ($out_format eq "fasta") {
      return();
    }

    ## Display options
    my $use_scores = 0; ## For the time being, we set this option to 0 because it supposes that scores are comprized between 160 and 1000, which is not the case for the features produced by RSAT (typically, weight scores are comprized between 0 and 20)

    ## Print format
    my $header = "";
    if ($out_format eq "gff3") {
      $header .= $comment_char{$out_format};
      $header .= "gff-version\t3";
      $header .= "\n";
    } elsif ($out_format eq "bed") {
      my $track_name =  "RSAT_features";
      if (defined($main::infile{input})) {
	$track_name = $main::infile{input};
      }
      $header .= "track name='".$track_name."' description='".$track_name."' visibility=3 itemRgb='On' use_score=".$use_scores."\n";
      $header .= "browser dense ".$track_name."\n";
#      $header .= "browser dense all\n";
    } else {
      $header .= $comment_char{$out_format};
      $header .= "Feature format:".$out_format."\n";
    }

    ## Print column content
    my @cols = @{$columns{$out_format}};
    $header .= $header_char{$out_format};
    $header .= join ("\t", @cols);
    $header .= "\n";
    return $header;
}


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

