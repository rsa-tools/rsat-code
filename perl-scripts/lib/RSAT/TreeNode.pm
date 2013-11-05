=pod

=head1 NAME

RSAT::TreeNode - A Tree object

=head1 AUTHOR

Email rekins@bigre.ulb.ac.be

=cut

package RSAT::TreeNode;
use vars qw(@ISA);
use RSAT::GenericObject;
use RSAT::error;
use RSAT::message;
use RSAT::Tree;

@ISA = qw( RSAT::GenericObject RSAT::Tree );

=pod

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

=pod

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

=pod

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


=pod

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

=pod

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

=pod

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

=pod

=head2 is this node a leaf ?

 Title   : is_leaf()
 Usage   : $node->is_leaf
 Function: check if this node is a leaf
 Returns : true if node is a leaf

=cut

sub is_leaf  {
  my ($self) = shift;
  my $isleaf=0;
#  if ((defined $self->{'child'} && (keys %{$self->{'child'}} > 0))&&
  if ($self->get_type eq "leaf"){
    RSAT::message::Warning(join("\t","Node",$self->getid(),"identified as leaf.")) if ($main::verbose >=10);
    $isleaf = 1;
  }
  return $isleaf;
}

=pod

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


=pod

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

=pod

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

=pod

=head2 get all descendents

 Title   : get_all_descendents()
 Usage   : my @descendants = $node->get_all_descendents()
 Function: get all descendant descendants from this node
 Returns : reference to an Array of descendent nodes

=cut

sub get_all_descendents{
  my $self=shift;
  my $order=shift||"DFS"; # BFS DFS
  my $type=shift||"all"; # all, leave, node
  my $max_depth=shift;
  my $max_leaves=shift;
  my @descendents =();

  if ($order eq 'DFS'){
    (@descendents) = $self->get_all_descendents_by_DFS($type,$max_depth,$max_leaves);
#  }elsif ($order eq 'BFS'){
#    (@descendents) = $self->get_all_descendents_by_BFS($type,$max_depth,$max_leaves);
  }else{
#    RSAT::message::Warning(0,join("\t","Order error : Please specify a order (DFS or BFS).")) if ($main::verbose >=0);
    RSAT::message::Warning(0,join("\t","Order error : Please specify a DFS order (BFS not yet implemented).")) if ($main::verbose >=0);
  }
  return (@descendents);
}


=head2 get all descendents DFS

  Title   : get_all_descendents_by_DFS()
  Usage   : my @descendants = $node->get_all_desceneants_by_DFS()
  Function: get all descendant descendants from this node by a breadth-first-search algorithm (DFS)
  Returns : reference to an Array of descendent nodes

=cut

sub get_all_descendents_by_DFS{
  
  my $self=shift;
  my $type=shift; # all, leave, node
  my $max_depth=shift;
  my $max_leaves=shift;
  my $depth=shift||0;
  my (@descendents) =();
  if ($max_depth){
    RSAT::message::Warning(1,join("\t",
				  "Node",$self->get_name(),
				  "Level",$self->get_level(),
				  "Depth",$depth,
				  "max_depth",$max_depth)) if ($main::verbose >= 3);
    if ($depth >= $max_depth){
      return (@descendents);
    }
  }
  $depth++;
  foreach my $child (sort {$a->get_name() cmp $b->get_name()} $self->get_children()) {
    if ($type eq "all"){
      push @descendents,$child,($child->get_all_descendents_by_DFS($type,$max_depth,$max_leaves,$depth));
      
    }elsif ($child->get_type() eq $type){
      if ($child->get_type() eq "leaf"){
	push @descendents,$child;
      }else{

	push @descendents,$child,($child->get_all_descendents_by_DFS($type,$max_depth,$max_leaves,$depth));
      }
 
    }else{
      if($child->get_type() eq "node"){
	push @descendents,($child->get_all_descendents_by_DFS($type,$max_depth,$max_leaves,$depth));
      }
      next;
    }
  }
  return (@descendents);
}

=pod

=head2 get all descendents by BFS (TO BE IMPLEMENTED)

 Title   : get_all_descendents_by_BFS()
 Usage   : my @descendants = $node->get_all_desceneants_by_BFS()
 Function: get all descendant descendants from this node by a breadth-first-search algorithm (BFS)
 Returns : reference to an Array of descendent nodes

=cut

### TO BE IMPLEMENTED


=pod

=head2 get leaves per node

 Title   : get_leaves
 Usage   : 
 Function: get the leaves under each node in a tree

=cut

sub get_leaves  {
   my ($self) = shift;
   my @leaves=();
   foreach my $node ( $self->get_all_descendents() ) {
     RSAT::message::Warning(join("\t","Descendant node",$node->get_name())) if ($main::verbose >=10);
     if ($node->is_leaf){
       push @leaves,$node;
     }
   }
   $self->{leaves}=@leaves;
   return (@leaves);
}

=pod

=head2 get_leaves_names()

 Title    : get_leaves_names()
 Usage    : my @leaves_labels = $node->get_leaves_names()
 Function : returns a list of node labels corresponding to the leaves
 Returns  : @leaves_labels

=cut

sub get_leaves_names {
  my $self = shift;
  my @leaves_labels=();
  foreach my $leaf ( $self->get_leaves()){
    push @leaves_labels, $leaf->get_name();
  }
  return (@leaves_labels);
}

1;
