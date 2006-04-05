###############################################################
#
# Class Template
#
package RSAT::Template;

use RSAT::GenericObject;
use RSAT::error;
@ISA = qw( RSAT::GenericObject );

### class attributes

=pod

=head1 NAME

    RSAT::Template

=head1 DESCRIPTION

Template class. 

=cut


################################################################
=pod

=item B<new()>

Create a new Template.

=cut
sub new {
    my ($class, %args) = @_;
    my $object = bless {
	}, $class;
    return $object;
}

return 1;

__END__

