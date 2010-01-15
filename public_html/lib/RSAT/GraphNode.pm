###############################################################
#
# Class GraphNode
#
package RSAT::GraphNode;

use RSAT::GenericObject;
use RSAT::GraphNeighbour;
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

=item B<get_neighbours>

Return the direct neighbours of the node.

Usage: my @neighbours = $node->get_neighbours();

The output object is a list containing objects of the class
RSAT::GraphNeighbour. 

=cut
sub get_neighbours {
  my ($self, $steps) = @_;
  my @neighbours = ();

  foreach my $arc ($self->get_attribute("in_arcs")) {
    my $neighbour = new RSAT::GraphNeighbour();

    my $node = $arc->get_attribute("source");
    $neighbour->set_attribute("node", $node);
    $neighbour->set_attribute("direction", "in");

    $neighbour->set_attribute("steps", 1);
#    my $weight = $arc->get_attribute("weight");
#    $neighbour->set_attribute("weight", $weight);
    push @neighbours, $neighbour;
  }

  foreach my $arc ($self->get_attribute("out_arcs")) {
    my $neighbour = new RSAT::GraphNeighbour();
    my $node = $arc->get_attribute("target");
    $neighbour->set_attribute("node", $node);
    $neighbour->set_attribute("direction", "out");

    $neighbour->set_attribute("steps", 1);
#    my $weight = $arc->get_attribute("weight");
#    $neighbour->set_attribute("weight", $weight);
    push @neighbours, $neighbour;
  }


  return @neighbours;
}


################################################################
=pod

=item B<get_neighbours_step>


Return the direct and indirect neighbours of the node up to a certain
number of steps.

Usage: my @neighbours = $node->get_neighbours_step($steps, $self_included);

Parameters: 

$steps: maximal number of steps from the seed node

$self_included: if this attribute is set, the seed node is included in
the neighbourhood (with a distance 0), even if there is no self-loop.

The output object is a list containing objects of the class
RSAT::GraphNeighbour.

=cut
sub get_neighbours_step {

  my ($self, $steps, $self_included) = @_;

  #  my @neighbours = ();
  my %neighbour_ids = ();
  my @seed_nodes = ();
  my @neighbours = ();

  $neighbour_ids{$self->get_attribute("id")} = 1; ## Avoid including the seed node in the neighbour node

  ## initialize the neighborhood by including the node itself, with a distance of 0
  my $self_neighbour = new RSAT::GraphNeighbour();
  my $self_id = $self->get_attribute("id");
  $self_neighbour->set_attribute("node", $self);
  $self_neighbour->set_attribute("direction", "self");
  $self_neighbour->set_attribute("steps", 0);
  #   $self_neighbour->set_attribute("weight", $weight, 0);
  $neighbours_ids{$self_id} = $self_neighbour;
  if ($self_included) {
    push @neighbours, $self_neighbour;
  }

  for my $s (1..$steps) {
    if ($s == 1) {
      @seed_nodes = ($self);
    } else {
      @seed_nodes = @new_neighbours;
    }
    #    &RSAT::message::TimeWarn("Getting neighbours for node", $self->get_attribute("id"), scalar(@seed_nodes), "seed nodes at step", $s) if ($main::verbose >= 0);
    @new_neighbours = ();
    foreach my $n (0..$#seed_nodes) {
      $node = $seed_nodes[$n];
      &RSAT::message::Debug("step", $s, "Getting neighbours for node", ($n+1)."/".scalar(@seed_nodes), $node->get_attribute("id")) if ($main::verbose >= 5);
      my @next_neighbours = $node->get_neighbours();
      foreach my $neighbour (@next_neighbours) {
	$neighb_node = $neighbour->get_attribute("node");
	## If a node is already in the neighborhood, avoid to incorporate it twice
	my $neighb_id = $neighb_node->get_attribute("id");
	if (defined($neighbour_ids{$neighb_id})) {
	  &RSAT::message::Info("node already included in neighbourhood", $neighb_id) if ($main::verbose >= 5);
	} else {
	  $neighbour->force_attribute("steps", $s);
	  $neighbour->force_attribute("direction", "na") if ($s > 1);
	  #	  $neighbour->force_attribute("weight", $s);
	  push @neighbours, $neighbour;
	  push @new_neighbours, $neighb_node;
	  $neighbour_ids{$neighb_id}++;
	}
      }
    }
    &RSAT::message::Debug("GraphNode::get_neighbours_step", 
			  "node", $self->get_attribute("id"),
			  "step", $s,
			  "seed nodes", scalar(@seed_nodes),
			  "new neighbours", scalar(@new_neighbours),
			  "total neighbours", scalar(@neighbours),
			 ) if ($main::verbose >= 4);
  }
  return @neighbours;
}


return 1;

__END__

