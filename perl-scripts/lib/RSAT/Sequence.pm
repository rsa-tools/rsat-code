###############################################################
#
# Class Sequence
#
package RSAT::Sequence;

use RSAT::GenericObject;
use RSAT::error;
@ISA = qw( RSAT::GenericObject );

=pod

=head1 NAME

    RSAT::Sequence

    =head1 DESCRIPTION

    Class for managing Sequence.

=cut


################################################################

=pod

=item new()

Create a new Sequence.

=cut
sub new {
    my ($class, %args) = @_;
    my $self = bless {
	id=>$args{id} || $class->auto_id(),
	sequence=>$args{sequence},
	type=>lc($args{type}) || "dna", ### dna (default) or aa
	source=>$args{source},
	description=>$args{description},
    }, $class;
    $self->init();
    return $self;
}


################################################################

=pod

=item Initiator.

initialize the sequence. Suppress all the spaces.

=cut
sub init {
    my ($self) = @_;
    $self->{sequence} =~ s/\s//g;
}


sub get_id {
    my ($self) = @_;
    return $self->{id};
}
sub get_type {
    my ($self) = @_;
    return $self->{type};
}


################################################################

=pod

=item  get_sequence($from, $to, $strand)

Return the whole sequence or a fragment (if arguements "from" and "to"
are specified).  The reverse complement can be obtained by specifying
"R" for the argument "strand".

=cut
sub get_sequence {
  my ($self, $from, $to, $strand) = @_;
  if (($from) && ($to) && ($strand)) {
    if ($to < $from) {
      my $tmp = $from;
      $from = $to;
      $to = $tmp;
    }

    my $fragment_length;
    my $fragment;
    if ($from < 0) {
      ## Negative coordinates indicate sequence fragment from the end
      ## of the sequence.  In this case, the fragment is either
      ## located entirely at the end (if $to <=0), or includes a piece
      ## from the end ($from < 0) and another piece from the beginning
      ## ($to > 0) of the sequence.

      ## Select the end fragment
      if ($to <= 0) {
	$fragment_length = $to - $from + 1;
      } else {
	$fragment_length = -$from;
      }
      my $seq_len = $self->get_length();
      $fragment =  substr($self->{sequence}, $seq_len + $from - 1, $fragment_length);

      ## Append start fragment if required
      if ($to > 0) {
	$fragment_length = $to;
	$fragment .= substr($self->{sequence}, $from-1, $fragment_length);
      }

    } else {
      ## Simple case: both from and to are positive
      $fragment_length = $to - $from + 1;
      $fragment = substr($self->{sequence}, $from-1, $fragment_length);
    }

    if ($strand eq "R") { ### reverse complement
      return (&main::ReverseComplement($fragment));
    } else {
      return $fragment;
    }
  } else {
    return $self->{sequence};
  }
}

################################################################

=pod

=item get_length

Return the sequence length

=cut
sub get_length {
    my ($self) = @_;
    return length($self->{sequence});
}

################################################################

=pod

=item get_source

Return the sequence source

=cut
sub get_source {
    my ($self) = @_;
    return $self->{source};
}


################################################################

=pod

=item get_description

Return the sequence description

=cut
sub get_description {
    my ($self) = @_;
    return $self->{description};
}

return 1;


__END__

