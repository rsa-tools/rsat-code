###############################################################
#
# Class blast_hit
#
package RSAT::blast_hit;

use RSAT::GenericObject;
@ISA = qw( RSAT::GenericObject);

=pod

=head1 NAME

    RSAT::blast_hit

=head1 DESCRIPTION

BLAST hit

=cut


################################################################
=pod

=item new()

Create a new blast_hit.

=cut
sub new {
    my ($class, %args) = @_;
    my $self = bless {
	%args
    }, $class;
    return $self;
}


return 1;


__END__

