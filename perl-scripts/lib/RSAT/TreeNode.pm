=head1 NAME

RSAT::TreeNode - A Tree object

=head1 SYNOPSIS


=head1 DESCRIPTION

This module make a Bio::Tree object from a hash asoociating organism name and their taxonomy.

=head1 AUTHOR

Email rekins@scmbb.ulb.ac.be

=cut

package RSAT::TreeNode;
use vars qw(@ISA);
#use RSAT::util;
use RSAT::GenericObject;
use RSAT::error;
use RSAT::message;
use RSAT::Tree;
#use RSAT::Object;
@ISA = qw( RSAT::GenericObject RSAT::Tree );

=head2 getid

 Title   : getid()

 Usage   : my $node = new RSAT::TreeNode("id"=>$nodeid,
					   "name"=>$nodename);
          $node->getid();

 Function: Get id of a tree node
 Returns : L<RSAT::TreeNode>

=cut


sub getid{
  my $self = shift;
  return $self->{id};
}

=head2 get_name

 Title   : get_name()
 Usage   : my $node = new RSAT::TreeNode("id"=>$nodeid,
					   "name"=>$nodename);
          $node->get_name();
 Function: Get name of a tree node
 Returns : L<RSAT::TreeNode>

=cut


sub get_name{
  my $self = shift;
  return $self->{name};
}

=head2 get_type

 Title   : get_type()
 Usage   : my $node = new RSAT::TreeNode("id"=>$nodeid,
                                            "name"=>$node_name,
					   "type"=>$node_type);
          $node->get_type();
 Function: Get type of a tree node
 Returns : string can be ["leaf","node","root"]

=cut

sub get_type{
  my $self = shift;
  return $self->{type};
}

=head2 get_level

 Title   : get_level()
 Usage   : my $node = new RSAT::TreeNode("id"=>$nodeid,
                                            "name"=>$node_name,
					   "level"=>$node_level);
          $node->get_level();
 Function: Get level of a tree node
 Returns : integer

=cut


sub get_level{
  my $self = shift;
  return $self->{level};
}

=head2 add_child

 Title   : add_child()
 Usage   : $parent_node->add_child($child_node);
 Function: add a child to a tree node
 Returns : L<RSAT::TreeNode>

=cut

sub add_child  {
  my ($self,$node) = @_;
  # add parent to the child node
#  $node->{"parent"=>$self};
  # add child 
  $self->{'child'}->{$node->getid()} = $node;
  return scalar keys %{$self->{'child'}};
}

=head2 get_children

 Title   : get_children()
 Usage   : @children =get_children($node);
 Function: get children of a tree node
 Returns : Array

=cut

sub get_children  {
  my ($self) = @_;
  my @children=();
  if (defined $self->{child}){
    @children = values %{$self->{child}};
    return @children;
  }else{
    return;
  }
}

=head2 is this node a leaf ?

 Title   : is_leaf()
 Usage   : $node->is_leaf
 Function: check if this node is a leaf
 Returns : true if node is a leaf

=cut

sub is_leaf  {
  my ($self) = shift;
  my $isleaf=0;
  if ((defined $self->{'child'} && (keys %{$self->{'child'}} > 0))&&
      ($self->get_type eq "leaf")){
    $isleaf = 1;
  }
  return $isleaf;
}

=head2 set level

 Title   : set_level()
 Usage   : $node->set_level($level)
 Function: attribute a level to this node
 Returns : level attributed to this node
 Argument: $level [int]

=cut

sub set_level  {
  my $self = shift;
  my $level = shift;
  $self->{"level"}= $level;
  RSAT::message::Warning(join("\t","Attribute level $level to node",$self->getid())) if ($main::verbose >=10);
  return $self->get_level();
}


=head2 set children levels

 Title   : set_children_levels()
 Usage   : $node->set_children_levels($pevious_level)
 Function: Attribute a level to each children node of a given node
 Returns : level
 Args    : $previous_level [int]

=cut

sub set_children_levels{
  my $self =shift;
  my $level=shift;
  $level++;
  foreach my $child ($self->get_children()){
    if ($self->is_leaf){
      next;
    }else{
      $child->set_level($level);
      $child->set_children_levels($level);
    }
  }
 return ();
}

=head2 get all nodes

 Title   : get_all_nodes()
 Usage   : my @descendants = $node->get_all_nodes()
 Function: get all descendant nodes from this node by a depht-first-search algorithm (DFS)
 Returns : Array of nodes

=cut

sub get_all_nodes{
  my $self=shift;
  my $nodes = shift;
  foreach my $child ($self->get_children()){
    if ($self->is_leaf){
      next;
    }else{
      push @{$nodes},$child;
      $child->get_all_nodes($nodes);
    }
  }
 return ($nodes);
}


1;
