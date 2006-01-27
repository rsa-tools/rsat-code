=head1 NAME

RSAT::Tree - A Tree object

=head1 SYNOPSIS

{
   my %organisms=();
   my $select = "SELECT names, taxonomy FROM organism";
   my $result = $dba->execute_SQL($select);
   while (my($name,$taxonomy) = $result->fetchrow_array() ) {
     $organisms{$name}=$taxonomy;
   }
   my $tree = &RSAT::Tree::makeTree("organism",%organisms);
}

=head1 DESCRIPTION

This module make a Bio::Tree object from a hash asoociating organism name and their taxonomy.

=head1 AUTHOR

Email rekins@scmbb.ulb.ac.be

=cut

package RSAT::Tree;
use vars qw(@ISA);
#use RSAT::util;
use RSAT::GenericObject;
use RSAT::error;
use RSAT::message;

#use RSAT::Object;
use RSAT::TreeNode;
use Data::Dumper; # useful for debugging (print contents of hashes, objects)
@ISA = qw( RSAT::GenericObject );
$default_indent_string = ":-";

################################################################
#### METHODS TO USE TREE
################################################################

=pod

=head2 set root node

 Title   : set_root_node($node)
 Usage   : $root = $tree->set_root_node($node)
 Function: Set the root node
 Returns : Root node object
 Args    : RSAT::TreeNode object

=cut

sub set_root_node  {
  my $self = shift;
  my $value = shift;
  $self->{'rootnode'} = $value;
  return $self->get_root_node;
}

=pod

=head2 get root node

 Title   : get_root_node($node)
 Usage   : $root = $tree->get_root_node($node)
 Function: Get the root node
 Returns : Root node object

=cut

sub get_root_node {
  my $self = shift;
  return $self->{'rootnode'};
}

=pod

=head2 set all levels

 Title   : set_all_levels()
 Usage   : my %node2level = $tree->set_all_levels($root_level)
 Function: Attribute a level to each node starting from the root ()
 Returns : Hash
 Args    : $root_level [int] (default:1)

=cut

sub set_all_levels{
  my $self =shift;
  my $level=1;
  my $root_node=$self->get_root_node();
  $root_node->RSAT::TreeNode::set_level($level);
  $root_node->set_children_levels($level);
  return();
}

=pod

=head2 get all nodes

 Title   : get_all_nodes()
 Usage   : my @descendants = $tree->get_all_nodes()
 Function: Get all nodes of the tree from the root by DFS algorithm
 Returns : Array of nodes

=cut

sub get_all_nodes{
  my $self =shift;
  my $root_node=$self->get_root_node();
  my ($descendants) = $root_node->get_all_nodes();
  return ($root_node,@{$descendants});
}

################################################################
#### IMPORT METHODS
################################################################

=pod

=head2 make a tree from ncbi taxonomy

 Title   : parse_NCBI_taxonomy()

 Usage   : my $tree = RSAT::Tree::parse_NCBI_taxonomy($rootname,\%supported_organisms)

 Function: Make a tree object from a hash 

 Returns : L<Bio::Tree::Tree>

 Args    :

      $rootname  [string] Name to be attributed to the root

      %supported_organisms [hash]   ( '$organism_name' => '$taxonomy')

=cut

sub parse_NCBI_taxonomy {
  my ($self,$root_name,$supported_organism)=@_;
  my %supported_organism=%{$supported_organism};
  my %nodes = (); # node index

  ## Initiate the root of the taxonomy
  my $root_node = new RSAT::TreeNode("id"=>$root_name,
					"name"=>$root_name,
					"type"=>"root"
				       );
  $nodes{$root_name} = $root_node;
  my $root=$self->set_root_node($root_node);
  RSAT::message::Warning("Root node :\t",$root->getid()) if ($main::verbose >= 3);
  
  ## get taxonomy
  my $c = 0;
  foreach my $org (sort {$supported_organism{$a}->{"taxonomy"} cmp $supported_organism{$b}->{"taxonomy"}} 
		   keys (%supported_organism)) {
    $c++;
    my @taxons = split /\s*;\s*/, $supported_organism{$org}->{"taxonomy"};
    &RSAT::message::Warning(join ("\t", $c, $org,scalar(@taxons),"taxons"), "\n") if ($main::verbose >=3);
    &RSAT::message::Warning(join ("\t","taxons",(@taxons)), "\n") if ($main::verbose >= 4);
    my $root_found=0;

    # initiate the leaf
    my $leaf = new RSAT::TreeNode(id=>$org,
				     name=>$org,
				     type=>"leaf"
				    );
    RSAT::message::Warning(join("\t","Initiate leaf",$leaf->get_name())) if ($main::verbose >= 4);

    for my $t (0..$#taxons) {
      RSAT::message::Warning(3,join("\t","Compare taxon",$taxons[$t],"with root name",$root->get_name())) if ($main::verbose >=3);
      ## identify root taxon
      if (($taxons[$t] eq $root->get_name())&&($root_found==0)){
	RSAT::message::Warning(3,"Taxon identified as root for\t",$org) if ($main::verbose >=3);
	$root_found=1;
	next;
      }
      # start top->down to increase the tree
      if ($root_found==1){
	if (defined $nodes{$taxons[$t-1]}){
	  my $node = new RSAT::TreeNode(id=>$taxons[$t],
					   name=>$taxons[$t],
					   type=>"node",
					   all_leaves=>[$org]
					  );
	  $nodes{$taxons[$t]}=$node;
	  RSAT::message::Warning(3,join("\t","Adding node",$node->get_name(),
					"to node",$nodes{$taxons[$t-1]}->get_name())) if ($main::verbose >=3);
	  $nodes{$taxons[$t-1]}->add_child($node);
	  # attach organism as leaf if it is the last taxon
	  if ($t == $#taxons){
	    RSAT::message::Warning(3,join("\t","Adding leaf",$leaf->get_name(),
					  "to node",$node->get_name())) if ($main::verbose >=3);
	    $node->add_child($leaf);
	  }
	}else{
	  next;
	}
      }
    }
  }
  return $self;
}


## ##############################################################
=pod

=item LoadSupportedTaxonomy

Fill a tree (RSAT::Tree) with the taxonomy of supported organisms on RSAT

Usage:  my $tree = SupportedOrganismTree($no_species);

Parameters:

=over

=item $no_species

do not create a node for the species, but only for

=back

=cut

sub LoadSupportedTaxonomy {
    my ($self, $no_species) = @_;

    my %nodes = (); # node index

    ## Initiate the root of the taxonomy
    my $root_node = new RSAT::TreeNode();
    $root_node->force_attribute("id", "Organism");
    $root_node->set_attribute("name", "Organism");
    $root_node->set_attribute("description", "Organism");
    $nodes{organism} = $root_node;
    $self->set_root_node($root_node);

    ## Iterate over all supported organisms
    my $org_counter = 0;
    foreach my $org (keys %main::supported_organism) {
	&RSAT::message::Info(join("\t", "adding organism", $org)) if ($main::verbose >= 3);
	$org_counter++;
	my $org_node = new RSAT::TreeNode();
	$org_node->force_attribute("id", $org);
	$org_node->set_attribute("name", $org);
	$org_node->set_attribute("description", $org);
	$nodes{$org} = $org_node; # index the new node

	my $taxonomy = $main::supported_organism{$org}->{taxonomy};

	## Replace prolematic characters by _
	$taxonomy = &RSAT::util::trim($taxonomy);
	$taxonomy =~ s|/|_|g; ## / are reserved in phylip format
	$taxonomy =~ s|; +|;|g; ## 
	$taxonomy =~ s| |_|g; ## 
	$taxonomy =~ s|\(|_|g; ## / are reserved in phylip format
	$taxonomy =~ s|\)|_|g; ## / are reserved in phylip format

	my @taxonomy = split /\s*;\s*/, $taxonomy;

	&RSAT::message::Info(join ("\t", $org_counter, $org, $taxonomy)), "\n" if ($main::verbose >= 4);
	
	## Initiate child to the level of the organism
	my $child = $org;
	$child_node = $org_node; 
	
	## ##############################################################
	## Traverse the current taxonomy bottom -> up (from species to
	## phyllum) and create nodes if they don't exist yet
	for my $tr (0..$#taxonomy) {
	    my $t = $#taxonomy -$tr;
	    my $parent = $taxonomy[$t];
	    if (defined $nodes{$parent}) {
		$nodes{$parent}->add_child($child_node) unless (($no_species) && ($child_node eq $org_node));
		warn join("\t", ";\t", "parent found", $t, $parent, $nodes{$parent}), "\n" if ($main::verbose >= 4);
		$child_node = $nodes{$parent};
		last;
	    } else {
		$parent_node = new RSAT::TreeNode();
		$parent_node->force_attribute("id", $parent);
		$parent_node->set_attribute("name", $parent);
		$parent_node->set_attribute("description", $parent);

		$nodes{$parent} = $parent_node;
		$nodes{$parent}->add_child($child_node ) unless (($no_species) && ($child_node eq $org_node));
		warn join("\t", ";\t", 
			  "new parent", $t, $parent, 
			  "child", $child_node->id(), 
			 ), "\n" if ($main::verbose >= 4);

		$child_node = $nodes{$parent};


		## Attach the top node to the root
		if ($t == 0) {
		    $root_node->add_child($child_node);
		}
	    }
	}
    }
}

################################################################
#### EXPORT METHODS
################################################################


################################################################
=pod

=head2 node_names()

 Title    : node_names()
 Usage    : my @node_labels = $tree->as_list()
 Function : returnsa a list of node labels
 Returns  : @node_labels

=cut

sub node_names {
    my ($self) = @_;
    my @node_names = ();
    my @nodes = $self->get_all_nodes();
    foreach my $node (@nodes) {
	push @node_labels, $node->get_attribute("name");
    }
    return @node_labels;
}


################################################################
=pod

=head2 as_indented_text()

 Title   : as_indented_text()
 Usage   : $tree->as_indented_text($indent_string)
 Function: Export tree as indented text.
           You can specify the string character to use for the indentation.
 Returns : $text_to_print
 Argument: $indent [string]

=cut
sub as_indented_text{
    my ($self,$indent_string) = @_;
    unless (defined($indent_string)) {
	$indent_string = $default_indent_string;;
    }
    my $output ="";
    $self->set_all_levels();
    foreach my $n ($self->get_all_nodes()){
	$output .= join(" ",$indent_string x $n->get_level(),$n->getid())."\n";
    }
    return ($output);
}

################################################################
#### CGI METHOD

=pod

=head2 export tree as a hash

 Title   : as_indented_hash()
 Usage   : $tree->as_indented_hash($indent_string)
 Function: Export tree as indented hash.
           You can specify the string character to use for the indentation.
           To be used in CGI form
 Returns : $hash (key=taxon, value=indented_taxon)
 Argument: $indent [string]

=cut

sub as_indented_hash{
  my ($self,$indent_string) = @_;
  unless (defined($indent_string)) {
      $indent_string = $default_indent_string;;
  }
  my %taxons =();
  $self->set_all_levels();
  foreach my $n ($self->get_all_nodes()){
    $taxons{$n->getid()} = join(" ",$indent_string x $n->get_level(),$n->getid())."\n";
  }
  return (%taxons);
}

1;
