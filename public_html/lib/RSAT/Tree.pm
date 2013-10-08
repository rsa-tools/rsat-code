package RSAT::Tree;

=pod

=head1 NAME

RSAT::Tree - Class for handling taxonomic trees in RSAT. 

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

Class for the maipulation of taxonomic trees in RSAT.

=head1 AUTHOR

Email rekins@bigre.ulb.ac.be

=cut

use vars qw(@ISA);
use RSAT::GenericObject;
use RSAT::error;
use RSAT::message;
use RSAT::TreeNode;

# useful for debugging (print contents of hashes, objects)
use Data::Dumper; 

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
    my ($self) =@_;
    my $level=1;
    my $root_node=$self->get_root_node();
    $root_node->set_level($level);
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
    my ($self) = @_;
    my $root_node=$self->get_root_node();
    my ($descendants) = $root_node->get_all_nodes();
    return ($root_node,@{$descendants});
}

=pod

=head2 get all descendents

 Title   : get_all_descendents()
 Usage   : my @descendants = $tree->get_all_descendents()
 Function: Get all descendents of the tree from the root by DFS algorithm
 Returns : Array of nodes

=cut

sub get_all_descendents{
    my ($self, $order, $type, $max_depth, $max_leaves) = @_;
#     my $self = =shift;
#     my $order=shift;
#     my $type=shift; # all, leave, node
#     my $max_depth=shift;
#     my $max_leaves=shift;
    my $root_node=$self->get_root_node();
    my (@descendents) = $root_node->get_all_descendents($order,$type,$max_depth,$max_leaves);
    return ($root_node,@descendents);
}

=pod

=head2 get node descendents

 Title   : get_node_descendents()
 Usage   : my @descendants = $tree->get_node_descendents($taxon, $order, $type, $max_depth, $max_leaves);
 Example : my @descendants = $tree->get_node_descendents("Gammaproteobacteria", "DFS", "all");
 Function: Get node descendents of the tree from the root by DFS algorithm
 Returns : Array of nodes

=cut

sub get_node_descendents{
  my ($self, $node_id, $order, $type, $max_depth, $max_leaves) = @_;
#  &RSAT::message::Info("RSAT::Tree", $self, "Getting node descendents", $node_id, $order, $type) if ($main::verbose >= 0);
  if ($node_id){
    my $node = $self->get_node_by_id($node_id);
#    &RSAT::message::Debug("RSAT::Tree", $node_id, "node", $node) if ($main::verbose >= 10);
    my (@descendents) = $node->get_all_descendents($order,$type,$max_depth,$max_leaves);
    return ($node,@descendents);
  }else{
    &RSAT::error::FatalError("&RSAT::Tree::get_node_descendents() was called without specifying \$node_id");
  }
}

=pod

=head2 get node descendents names

 Title   : get_node_descendents()
 Usage   : my @descendants = $tree->get_node_descendents(@args)
 Arguments: all arguments are passed to RSAT::Tree::get_node_descendents
 Function: Get node descendents of the tree from the root by DFS algorithm
 Returns : Array of nodes

=cut

sub get_node_descendents_names {
  my ($self, @args) = @_;
  &RSAT::message::Info("RSAT::Tree", $self, "Getting names of node descendents") if ($main::verbose >= 4);
#  my $node_id=shift;
#  my $order=shift;
#  my $type=shift; # all, leave, node
#  my $max_depth=shift;
#  my $max_leaves=shift;
#  my (@descendents) = $self->get_node_descendents($node_id,$order,$type,$max_depth,$max_leaves);
  my (@descendents) = $self->get_node_descendents(@args);
  my @node_names=();
  foreach my $n (@descendents){
    push @node_names, $n->getid();
  }
  return (@node_names);
}


=pod

=head2 get node by id

 Title   : get_node_by_id()
 Usage   : my $node = $tree->get_node_by_id($id)
 Function: Return a node object if exists in the tree.
 Returns : node object

=cut

sub get_node_by_id {
  my ($self, $id) = @_;
#  my $self = shift;
#  my $id = shift;
  my $rootnode = $self->get_root_node();

#  &RSAT::message::Debug("RSAT::Tree", $self, "Getting node by id", $id, $rootnode) if ($main::verbose >= 5);

  if ( ($rootnode->getid) && ($rootnode->getid eq $id) ) {
    return $rootnode;
  }
  foreach my $node ($rootnode->get_all_descendents(undef,"all") ) { ## SB: The type was set to "node" but I changed it to all which seemes more logical to me!
    
    my $current_id = $node->get_id();
    if ($current_id eq $id) {
      
      return $node;
    }
#    &RSAT::message::Debug("&RSAT::Tree::get_node_by_id()",$current_id, "differs from", $node_id) if ($main::verbose >= 10);
  }
  &RSAT::message::Warning("&RSAT::Tree::get_node_by_id()", "no node with id", $id);
  return(undef);
}

# =pod 

# =head2 get all descendents

#  Title   : get_all_descendents()
#  Usage   : my @descendants = $tree->get_all_descendents()
#  Function: Get all descendents of the tree from the root by DFS algorithm
#  Returns : Array of nodes

# =cut

# sub get_all_descendents{
#   my $self =shift;
#   my $order=shift;
#   my $type=shift; # all, leave, node
#   my $max_depth=shift;
#   my $max_leaves=shift;
#   my $root_node=$self->get_root_node();
#   my (@descendents) = $root_node->get_all_descendents($order,$type,$max_depth,$max_leaves);
#   return ($root_node,@descendents);
# }

# =pod

# =head2 get node descendents

#  Title   : get_node_descendents()
#  Usage   : my @descendants = $tree->get_node_descendents("Gammaproteobacteria")
#  Function: Get node descendents of the tree from the root by DFS algorithm
#  Returns : Array of nodes

# =cut

# sub get_node_descendents{
#   my $self =shift;
#   my $node_id=shift;
#   my $order=shift;
#   my $type=shift; # all, leave, node
#   my $max_depth=shift;
#   my $max_leaves=shift;
#   if ($node_id){
#     my $node=$self->get_node_by_id($node_id);
#     my (@descendents) = $node->get_all_descendents($order,$type,$max_depth,$max_leaves);
#     return ($node,@descendents);
#   }else{
#     die("No valid node id !");
#   }
# }


# =pod

# =head2 get node by id

#  Title   : get_node_by_id()
#  Usage   : my $node = $tree->get_node_by_id($id)
#  Function: Return a node object if exists in the tree.
#  Returns : node object

# =cut

# sub get_node_by_id  {
#   my $self = shift;
#   my $id = shift;
#   my $rootnode = $self->get_root_node();
#   if ( ($rootnode->getid) && ($rootnode->getid eq $id) ) {
#     return $rootnode;
#   }
#   foreach my $node ( $rootnode->get_all_descendents(undef,"node") ) {
#     if ( ($node->getid) and ($node->getid eq $id ) ) {
#       return $node;
#     }
#   }
#   return(0);
# }

################################################################
=pod

=head2 get_leaves_names()

 Title    : get_leaves_names()
 Usage    : my @leaves_labels = $tree->get_leaves_names()
 Function : returns a list of node labels corresponding to the leaves
 Returns  : @leaves_labels

=cut

sub get_leaves_names {
    my $self = shift;
    my $root_node=$self->get_root_node();
    my (@leaves_labels) = $root_node->get_leaves_names();
    return (@leaves_labels);
}

################################################################
#### IMPORT METHODS
################################################################

## ##############################################################
## TEMP: allows to choose between two versions of LoadSupportedTaxonomy
## In order to compare the results.
sub LoadSupportedTaxonomy {
  my ($self) = @_;
  &RSAT::message::Info("RSAT::Tree", $self, "Loading supported taxonomy") if ($main::verbose >= 3);
  #  &LoadSupportedTaxonomy_jvh(@_);
  &LoadSupportedTaxonomy_rj(@_);
  $self->force_attribute("loaded", 1)
;
}


=pod

=head2 make a tree from taxonomy

 Title   : loadSupportedTaxonomy()
 Usage   : my $tree = RSAT::Tree::loadSupportedTaxonomy($rootname,\%supported_organisms)
 Function: Make a tree object from a hash
 Returns : L<RSAT::Tree>
 Args    :
      $rootname  [string] Name to be attributed to the root
      %supported_organisms [hash]   ( '$organism_name' => '$taxonomy')

=cut

sub LoadSupportedTaxonomy_rj {
    my ($self,$root_name,$supported_organism) = @_;
    unless ($root_name) {
	$root_name = "Organism";
    }

    &RSAT::message::Info("RSAT::Tree", "Loading supported taxonomy with root", $root_name)
	if ($main::verbose >= 3);

    my %supported_organism=%{$supported_organism};
    my %nodes = ();		# node index

    ## Instantiate the root of the taxonomy
    my $root_node = new RSAT::TreeNode("id"=>$root_name,
				       "name"=>$root_name,
				       "type"=>"root"
	);
    $nodes{$root_name} = $root_node;
    my $root=$self->set_root_node($root_node);

    ## Get  the taxonomy
    my $c = 0;
    foreach my $org (sort {$supported_organism{$a}->{"taxonomy"} cmp $supported_organism{$b}->{"taxonomy"}}
		     keys (%supported_organism)) {

	## Define organism identifier
	my $org_id = $org;
	$org_id =~ s/\s+/_/g;

	$c++;
	my $taxonomy = $supported_organism{$org}->{"taxonomy"};
	$taxonomy = &RSAT::util::trim($taxonomy);
	$taxonomy =~ s/\s+/ /g;
	$taxonomy =~ s/\.$//g;
	my @taxa = split /\s*;\s*/, $taxonomy;
	@taxa = map {$_ =~ s/\s+/\_/g;$_} @taxa; # removing spaces in the top and bottom taxa names
	@taxa = map {$_ =~ s/(\(|\))/\_/g;$_} @taxa; # removing spaces in the top and bottom taxa names
	&RSAT::message::Warning($c, $org,scalar(@taxa),"taxa"),  if ($main::verbose >= 5);
	&RSAT::message::Warning($c, $org_id, "taxa", join ("::",(@taxa))) if ($main::verbose >= 0);

	# Instantiate the leaf
	my $leaf = new RSAT::TreeNode(id=>$org_id,
				      name=>$org,
				      type=>"leaf"
	    );
	$nodes{$org_id}=$leaf; ## Index the leaves by organism ID
	&RSAT::message::Warning(join("\t","Instantiated leaf",$leaf->get_name())) if ($main::verbose >= 5);

#    &RSAT::message::Debug("Created TreeNode for organism", $org, $org_id, $species_node) if ($main::verbose >= 10);



	for my $t (0..$#taxa) {

#	    # TEMPORARY PATCH
#	    # correct the taxon name for weird taxon name due to parsing error (cases of Salmonella enterica)
#	    if (($taxa[$t] =~ "^SC-B67")||($taxa[$t] =~ "^9150")) {
#		$taxa[$t]="Bacteria";
#	    }

	    # start top->down to build the tree
	    if (defined $nodes{$taxa[$t]}) {
		if ($t == $#taxa) {
		    ## Link the lowest-level taxon to the organism node (the leaf)
		    $nodes{$taxa[$t]}->add_child($leaf);
		} else {
		    next;
		}
	    } else {
		my $node = new RSAT::TreeNode(id=>$taxa[$t],
					      name=>$taxa[$t],
					      type=>"node",
					      ## all_leaves=>[$org]
		    );
		$nodes{$taxa[$t]}=$node; ## Index taxa by name
		&RSAT::message::Debug("Created TreeNode for taxon", $t, $taxa[$t], $node) 
		    if ($main::verbose >= 5);
		
		if ((defined $nodes{$taxa[$t-1]})&&($t-1>=0)) {
		    $nodes{$taxa[$t-1]}->add_child($node);
		} else {
		    # attach first taxon to the root
		    $nodes{$root_name}->add_child($node);
		}

		# attach organism as leaf if it is the last taxon
		if ($t == $#taxa) {
		    $node->add_child($leaf);
		}
	    }
	}
    }
    return $self;
}


=pod

=item LoadSupportedTaxonomy

Fill a tree (RSAT::Tree) with the taxonomy of supported organisms on RSAT

#Usage:  my $tree = &SupportedOrganismTree($no_species);
 Usage:  $tree->LoadSupportedTaxonomy_jvh($no_species);

Parameters:

=over

=item $no_species

Do not create a node for the species, but only for higher taxonomical levels.

=back

=cut

sub LoadSupportedTaxonomy_jvh {
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
    &RSAT::message::Info(join("\t", "adding organism", $org)) if ($main::verbose >= 4);
    $org_counter++;
    my $org_node = new RSAT::TreeNode();
    $org_node->force_attribute("id", $org);
    $org_node->set_attribute("name", $org);
    $org_node->set_attribute("description", $org);
    $nodes{$org} = $org_node; # index the new node

    my $taxonomy = $main::supported_organism{$org}->{taxonomy};

    ## Replace problematic characters by _
    $taxonomy = &RSAT::util::trim($taxonomy);
    $taxonomy =~ s|/|_|g; ## / are reserved in phylip format
    $taxonomy =~ s|; +|;|g; ## 
    $taxonomy =~ s| |_|g; ## 
    $taxonomy =~ s|\(|_|g; ## / are reserved in phylip format
    $taxonomy =~ s|\)|_|g; ## / are reserved in phylip format

    my @taxonomy = split /\s*;\s*/, $taxonomy;

    &RSAT::message::Info(join ("\t", $org_counter, $org, $taxonomy), "\n") if ($main::verbose >= 4);

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
	&RSAT::message::Debug("\t\tparent found", $t, $parent, $nodes{$parent}) if ($main::verbose >= 4);
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
#### EXPORT METHODS
################################################################

=pod

=head2 export tree as indented text

 Title   : as_indented_text()
 Usage   : $tree->as_indented_text($indent_string,$start_node_id)
 Function: Export tree as indented text.
           You can specify the string character to use for the indentation.
 Returns : $text_to_print
 Argument: $indent [string]
           $start_node_id [string]

=cut

sub as_indented_text{
  my $self = shift;
  my $indent_string=shift||"-";
  my $start_node_id = shift||$self->get_root_node()->getid();
  my $format = shift||"";
  my $node_type=shift||"all";
  my $max_depth=shift;
  my $output ="";
  $output = "<HTML><HEAD><TITLE>Taxonomic Tree - $start_node_id</TITLE></HEAD><BODY><PRE>\n"  if ($format =~ /^html/i);
  $self->set_all_levels();
  my $start_node=$self->get_node_by_id($start_node_id);

  if (! $start_node){
    die("No node with this id in the tree : \"$start_node_id\" !");
  }
  my $initlevel = $start_node->get_level();

  foreach my $n ($start_node,$start_node->get_all_descendents("DFS",$node_type,$max_depth,undef)){
    if (($n->is_leaf())&&($format =~ /^HTML/i)){
      $output .= join(" ",$indent_string x ($n->get_level() - $initlevel),"<i>",$n->getid())."</i>\n";
    }elsif($format =~ /^HTML/i){
      $output .= "<b>".join(" ",$indent_string x ($n->get_level() - $initlevel),$n->getid())."</b>\n";
    }else{
      $output .= join(" ",$indent_string x ($n->get_level() - $initlevel),$n->getid())."\n";

    }
  }
  $output.= "</PRE></BODY></HTML>\n"  if ($format =~ /^HTML/i);
  return ($output);
}

################################################################
#### EXPORT METHODS
################################################################

=pod

=head2 export tree as newick format (format used to represent taxonomical trees)

 Title   : as_newick()
 Usage   : $tree->as_newick($start_node_id)
 Function: Export tree as newick 
 Returns : $text_to_print
 Argument: $indent [string]
           $start_node_id [string]

=cut

sub as_newick  {
  my $self = shift;
  my $taxon = shift || "Organisms";
  our %newick_results = ();
  our %parents = ();
  $root = $self->get_node_by_id($taxon);
  
  
  my $root_id = $root->getid;
  &create_newick($self, $root_id, 1);
  my $output = "((".(join ",", @{$newick_results{"$taxon"}}).")$taxon)";
  return ($output);
}



sub create_newick {
  my ($self, $taxonid, $depth) = @_;
  my ($what, @descendants) = $self->get_node_descendents($taxonid, "DFS", "all", 1, undef);
  my $descendants_nb = scalar @descendants;
  if ($descendants_nb > 0) {
    foreach my $child (@descendants) {
      my $childid = $child->getid;
      $parents{$childid} = $taxonid;
      &create_newick($self,$childid, 1);
    }
  }
  my @species = @{$newick_results{$taxonid}};


  my $species_group = "";

  if (scalar @species > 0) {
    $species_group = "(".(join ",", @species).")$taxonid";
  } else {
    $species_group = $taxonid;
  }
  my $parent = $parents{$taxonid};
  push @{$newick_results{$parent}}, $species_group if $parent;

}


################################################################
#### CGI METHOD

=pod

=head2 export tree as a hash

 Title   : as_indented_hash()
 Usage   : $tree->as_indented_hash($indent_string,$start_node_id)
 Function: Export tree as indented hash.
           You can specify the string character to use for the indentation.
           To be used in CGI form
 Returns : $hash (key=taxon, value=indented_taxon)
 Argument: $indent [string]
           $start_node_id [string]

=cut

sub as_indented_hash{
  my ($self, $indent_string, $start_node_id) = @_;
  unless (defined($indent_string)) {
    $indent_string = $default_indent_string;;
  }
  my %taxa =();
  $self->set_all_levels();

  my $start_node=$self->get_node_by_id($start_node_id);
  if (! $start_node){
      die("No node with this id in the tree : \"$start_node_id\" !");
  }
  my $initlevel = $start_node->get_level();
  foreach my $n ($start_node,$start_node->get_all_descendents("DFS","node",undef,undef)){
      if ($n->is_leaf()){
	  die("This node must not be a leaf ! ".$n->getid());
      }else{
	  $taxa{$n->getid()} = join(" ",$indent_string x ($n->get_level() - $initlevel),$n->getid())."\n";
      }
  }
  return (%taxa);
}

1;
