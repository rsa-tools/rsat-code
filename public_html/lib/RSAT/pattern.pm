###############################################################
#
# Class pattern
#
package RSAT::pattern;

use RSAT::GenericObject;
use RSAT::Sequence;
use RSAT::error;

@ISA = qw( RSAT::GenericObject RSAT::Sequence );

=pod

=head1 NAME

    RSAT::pattern

=head1 DESCRIPTION

Class for managing patterns.

=cut



=pod

=item new()

Create a new pattern.

=cut

### creator
sub new {
    my ($class, %args) = @_;
    my $self = bless {
	%args,
	id=>$args{id} || $args{seq},
	type=>"pattern",
	sequence=>$args{sequence},
	source=>$args{source},
	description=>$args{description},
	score=>$args{score},
    }, $class;
    foreach my $key (keys %args) {
	$self->force_attribute($key, $args{$key});
    }
    return $self;
}


=pod

=item contains()

Pattern comparison method. 

Returns 1 if the pattern contains a match with the other
pattern. Matching includes degenerate code, for instance W matches M
since they have A as a possible intersection

Usage: if ($pattern1->contains($pattern2)) { ... } if
    ($pattern1->contains($pattern2), rc=>1) { ... }

Options: 
    rc=>1: try to match reverse complement


=cut
sub contains {
    my ($pattern1, $pattern2,%args) = @_;
    my $seq1 = $pattern1->get_sequence();
    my $seq2 = $pattern2->get_sequence();
    my $len1 = length($seq1);
    my $len2 = length($seq2);
    my $max_offset = $len1 - $len2;
    my @seqs = ();
    my $min_matches = $len2;
    my $min_score = 0;
    if ($args{min_score}) {
	$min_score = $args{min_score};
    } elsif ($args{perfect}) {
	$min_score = $len2;
    }

    return 0 if ($len2 > $len1);
    push @seqs, $seq2;
    unless ($args{single_strand}){ ### reverse complement
	push @seqs, &main::ReverseComplement($seq2) 
				 }

    foreach my $seq (@seqs) {
	for my $i (0..$max_offset) {
	    my $subseq = substr($seq1,$i,$len2);
	    my ($matches, $score) = &main::CountMatches($subseq,$seq);
	    return $score if (($matches >= $min_matches) &&
			      ($score >= $min_score)); 
	}
    }
    return 0;
}


return 1;


__END__

