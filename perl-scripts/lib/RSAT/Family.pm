###############################################################
#
# Class Family
#
package RSAT::Family;

use RSAT::GenericObject;
use RSAT::error;
@ISA = qw( RSAT::GenericObject );

### class attributes
$_count = 0;
$_prefix = "ctg_";
@_objects = ();
%_name_index = ();
%_id_index = ();
%_attribute_count = ();
%_attribute_cardinality = (id=>"SCALAR",
			   name=>"SCALAR",
			   organism=>"SCALAR",
			   size=>"SCALAR"
			  );

=pod

=head1 NAME

    RSAT::Family

=head1 DESCRIPTION

Class used to store a family (cluster) of genes.

=cut



################################################################
=pod

get_members()

Return the list of members

=cut

sub get_members {
    my ($self) = @_;
    return @{$self->{members}};
}

################################################################
=pod

get_size()

return the number of members

=cut
sub get_size {
    my ($self) = @_;
    my $size = scalar(@{$self->{members}});
    return $size;
}

################################################################
=pod

new_member($member)

Add a member to the family

=cut
sub new_member {
    my ($self, $new_member) = @_;
    push @{$self->{members}}, $new_member;
}


return 1;

__END__

