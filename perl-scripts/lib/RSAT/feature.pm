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

# RSAT feature-map format
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

# RSAT dna-pattern format
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

# RSAT Genome feature format
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

# Sanger general feature format
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

=over

=item gft

RSAT Genome features (file Feature.tab in the directory data/genomes)

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


=item  gff: Sanger general feature file (extension .gff)

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

=head1 METHODS

=over

=cut


################################################################
=pod

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
	$self->set_attribute($cols[$c], $fields[$c]);
    }

    ## Convert strand format
    my $strand = $self->get_attribute("strand");
    $strand =~ s/\+/D/;
    $strand =~ s/\-/R/;
    $strand =~ s/\./DR/;
    $self->force_attribute("strand", $strand);

    ## Format-specific conversions

    ## dna-pattern
    if ($format eq "dnapat") {
	$self->force_attribute("type", "dnapat");
    }
    &RSAT::message::Info(join("\t",
			      "new feature",			      
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
sub to_text {
    my ($self, $format, $null) = @_;

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
	warn join ( "\t", "field", $c, $fields[$c]), "\n" if ($main::verbose >= 10);
    }

    ## Format-specific features

    ## GFF 
    if ($format eq "gff") {
	my $attribute = "gene \"".$self->get_attribute("feature_name")."\"";
	$attribute .= "; note \"".$self->get_attribute("description")."\"";
	$fields[$col_index{attribute}] = $attribute;
    }


    ## Specific treatment for the strand
    my @strands = @{$strands{$format}};
    my $strand = $self->get_attribute("strand");
    my $f = $col_index{"strand"};
    my $s = $strand_index{$strand};
    $fields[$f] = $strands[$s];
#    die join( "\t", "f=$f", "s=$s", $strands[$s], "field=$fields[$f]");

    ## Generate the row to be printed
    my $row = join ("\t", @fields);
    $row .= "\n";

    warn (join "\t", "printing in format ", $format, $row) if ($main::verbose >= 10);

    return $row;
}

################################################################
=pod

=item header($format)

Print the header in the specified format.

=cut

sub header {
    my ($self, $format) = @_;

    my @cols = @{$columns{$format}};

    my $header = $comment_char{$format};
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

