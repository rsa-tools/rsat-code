###############################################################
#
# Class GraphNode
#
package RSAT::GraphNode;

use RSAT::GenericObject;
use RSAT::error;
@ISA = qw( RSAT::GenericObject );

### class attributes

=pod

=head1 NAME

    RSAT::GraphNode

=head1 DESCRIPTION

GraphNode class. 

=cut

################################################################
=pod 

=item B<node_neighbours>

Return the direct neighbours of the node.

Usage: my @neighbours = $node->get_neighbours();

=cut
sub get_neighbours {
  my ($self) = @_;
  my @neighbours = ();
  push @neighbours, $self->get_attribute("in_nodes");
  push @neighbours, $self->get_attribute("out_nodes");
  return @neighbours;
}



return 1;

__END__

