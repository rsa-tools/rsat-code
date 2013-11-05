###############################################################
#
# Class match
#
package RSAT::match;

use RSAT::GenericObject;
use RSAT::error;
@ISA = qw( RSAT::GenericObject  RSA::Sequence );

=pod

=head1 NAME

    RSAT::match

=head1 DESCRIPTION

Match of a pattern against a sequence.

=cut


################################################################
=pod

=item new()

Create a new match.

=cut
sub new {
    my ($class, %args) = @_;
    my $self = bless {
	id=>$args{id} || $args{seq},
	start_pos=>$args{start_pos},
	end_pos=>$args{end_pos},
	strand=>$args{strand},
	sequence=>$args{sequence},
	score=>$args{score},
    }, $class;
	return $self;
}

################################################################
=pod

=item get_length

Calculate match length from start and end positions

=cut
sub get_length {
    my ($self) = @_;
    my $length = abs($self->get_attribute("end_pos") - $self->get_attribute("start_pos") + 1);
    return $length;
}

return 1;


__END__

