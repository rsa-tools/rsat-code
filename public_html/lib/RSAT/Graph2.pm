###############################################################
#
# Class Graph2
#
package RSAT::Graph2;

use RSAT::GenericObject;
use RSAT::error;
use RSAT::util;



### class attributes
@ISA = qw( RSAT::GenericObject );
#$node_color="#000088";
#$arc_color="#000044";

=pod

=head1 NAME

    RSAT::Graph

=head1 DESCRIPTION

Implementation of basic graph functions. This class allows to create
directed graphs and export them in different formats.
Graph class. 


=head1 METHODS

=cut

################################################################
=pod

=item B<new()>

Creates a new graph.

=cut

sub new {
     my ($class, %args) = @_;
     my $self = bless {
	}, ref($class) || $class;
    $self->init(%args);
    my @out_neighbours =();
    my @in_neighbours = ();
    my @arc_in_label = ();
    my @arc_out_label = ();
    my @arc_in_color = ();
    my @arc_out_color = ();
    my @arcs = ();
    
    my $max_arc = 1;
    
    my %arcs_name_id = ();
    my %nodes_color = ();
    my %nodes_label = ();
    my %node_name_id = ();
    my %node_id_name = ();
    # @out_neighbours is a multidimensional array : the index corresponds to the id of the node 
    # @{$out_neighbours[x]} contains a list of the id of the out neighbours of the nodes having id x
    $self->set_array_attribute("out_neighbours", @out_neighbours);
    # @in_neighbours is a multidimensional array : the index corresponds to the id of the node 
    # @{$in_neighbours[x]} contains a list of the id of the in neighbours of the nodes having id x
    $self->set_array_attribute("in_neighbours", @in_neighbours);
    # @arc_in_label, @arc_out_label, @arc_in_color, @arc_out_color are multidimensional arrays : the index corresponds to the id of the node 
    # @{$arc_in_label[x]} contains a list of the label (weight), colors of the in_arcs / out_arcs corresponding to the 
    # in_neighbours / out_neigbours
    $self->set_array_attribute("in_label", @arc_in_label);
    # @arc_out_label is a multidimensional array : the index corresponds to the id of the node 
    # @{$arc_out_label[x]} contains a list of the label (weight) of the out_arcs corresponding to the 
    # out_neighbours
    $self->set_array_attribute("out_label", @arc_out_label);
    $self->set_array_attribute("in_color", @arc_in_color);
    $self->set_array_attribute("out_color", @arc_out_color);
    # @arcs is a multidimensional array[nX3], with n the number of edges of the graph.
    # this array is redundant with the above information but help in sparing a lot of computation time
    # $arcs[x][0] : sourcenode
    # $arcs[x][1] : targetnode
    # $arcs[x][2] : label
    # $arcs[x][3] : color
    $self->set_array_attribute("arcs", @arcs);
    # $maxarc represents the number of time that an arc may be duplicated.
    $self->force_attribute("nb_arc_bw_node", $max_arc);
    # %nodes_name_id correspondance between the name of an arc and its internal id in the @arcs array.
    # the key of the hash is source_target_x, with x an integer from 1 to $max_arc;
    $self->set_hash_attribute("arcs_name_id", %arcs_name_id);
   
    # %nodes_name_id correspondance between name and id
    $self->set_hash_attribute("nodes_name_id", %nodes_name_id);
    # %nodes_name_id correspondance between id and name
    $self->set_hash_attribute("nodes_id_name", %nodes_id_name);
    # %nodes_name_id correspondance between node_name and color
    $self->set_hash_attribute("nodes_color", %nodes_color);
    # %nodes_name_id correspondance between node_name and label
    $self->set_hash_attribute("nodes_label", %nodes_label);
    return $self;
}


################################################################
=pod

=item B<get_node_color()>

get the color of the node

=cut
sub get_node_color {
    my ($self, $node_name) = @_;
    my $numId = $self->node_by_name($node_name);
    my %nodes_color = $self->get_attribute("nodes_color");
    my $color = $nodes_color{$numId};
    return $color;
}

################################################################
=pod

=item B<get_out_labels()>

returns the labels of the out neighbours of a node

=cut
sub get_out_labels {
    my ($self, $node_name) = @_;
    my $numId = $self->node_by_name($node_name);
    my @out_label = $self->get_attribute("out_label");
    @node_out_label = ();
    if (defined(@{$out_label[$numId]})) {
      @node_out_label = @{$out_label[$numId]};
    }
    return @node_out_label;
}

################################################################
=pod

=item B<get_out_labels_id()>

returns the labels of the out neighbours of a node given its internal id

=cut
sub get_out_labels_id {
    my ($self, $numId) = @_;
    my @out_label = $self->get_attribute("out_label");
    @node_out_label = ();
    if (defined(@{$out_label[$numId]})) {
      @node_out_label = @{$out_label[$numId]};
    }
    return @node_out_label;
}

################################################################
=pod

=item B<get_out_color()>

returns the color of the out arcs of a node

=cut
sub get_out_colors {
    my ($self, $node_name) = @_;
    my $numId = $self->node_by_name($node_name);
    my @out_color = $self->get_attribute("out_color");
    my @node_out_color = ();
    if (defined(@{$out_color[$numId]})) {
      @node_out_color = @{$out_color[$numId]};
    }
    return @node_out_color;
}


################################################################
=pod

=item B<get_out_color()>

returns the color of the out arcs of a node given its id

=cut
sub get_out_colors_id {
    my ($self, $numId) = @_;
    my @out_color = $self->get_attribute("out_color");
    my @node_out_color = ();
    if (defined(@{$out_color[$numId]})) {
      @node_out_color = @{$out_color[$numId]};
    }
    return @node_out_color;
}


################################################################
=pod

=item B<contains_node>

returns 1 if the graph contains a node having name $node_name, returns 0 otherwise

=cut
sub contains_node {
    my ($self, $node_name) = @_;
    my $numId = $self->node_by_name($node_name);
    my $contains = 0;
    if (defined($numId)) {
      $contains = 1;
    }
    return $contains;

}





################################################################
=pod

=item B<get_out_neighbours()>

returns the out labels of a node

=cut
sub get_out_neighbours {
    my ($self, $node_name) = @_;
    my $numId = $self->node_by_name($node_name);
    my @out_neighbours = $self->get_attribute("out_neighbours");
    my @out_neighbours_names = ();
    if (defined($out_neighbours[$numId])) {
      my @out_neighbours_indices = @{$out_neighbours[$numId]};
      foreach my $out_neighbour_id(@out_neighbours_indices) {
        my $out_neighbour_name = $self->node_by_id($out_neighbour_id);
        push @out_neighbours_names, $out_neighbour_name;
      }
    } else {
      &RSAT::message::Warning("Node $node_name has no out neighbours") if $main::verbose >= 3;
    }
    return @out_neighbours_names;
}

################################################################
=pod

=item B<get_neighbours_id()>

returns an array contening the out neighbours of a node to a certain step self included or not
1st col : neighbour
2nd col : 


=cut
# sub get_neighbours {
#     my ($self, $node_name) = @_;
#     my $numId = $self->node_by_name($node_name);
#     my @out_neighbours = $self->get_attribute("out_neighbours");
#     my @out_neighbours_names = ();
#     if (defined($out_neighbours[$numId])) {
#       my @out_neighbours_indices = @{$out_neighbours[$numId]};
#       foreach my $out_neighbour_id(@out_neighbours_indices) {
#         my $out_neighbour_name = $self->node_by_id($out_neighbour_id);
#         push @out_neighbours_names, $out_neighbour_name;
#       }
#     } else {
#       &RSAT::message::Warning("Node $node_name has no out neighbours") if $main::verbose >= 3;
#     }
#     return @out_neighbours_names;
# }

################################################################
=pod

=item B<properties()>

Return properties of the graph (edges, nodes).

Supported formats: ???

=cut
sub properties {
    my ($self) = @_;
    my @arcs = $self->get_attribute("arcs");
    my $nbarcs = scalar @arcs;
    my %nodes_name_id = $self->get_attribute("nodes_name_id");
    my $nbnodes = scalar (keys(%nodes_name_id));    
    
return ($nbnodes, $nbarcs);
}
################################################################
=pod

=item B<get_out_neighbours()>

returns the out labels of a node given its internal id

=cut
sub get_out_neighbours_id {
    my ($self, $numId) = @_;
    my @out_neighbours = $self->get_attribute("out_neighbours");
    my @out_neighbours_names = ();
    if (defined($out_neighbours[$numId])) {
      my @out_neighbours_indices = @{$out_neighbours[$numId]};
      foreach my $out_neighbour_id(@out_neighbours_indices) {
        my $out_neighbour_name = $self->node_by_id($out_neighbour_id);
        push @out_neighbours_names, $out_neighbour_name;
      }
    } else {
      my $node_name = $self->node_by_id($numId);
      &RSAT::message::Warning("Node $node_name has no out neighbours") if $main::verbose >= 3;;
    }
    return @out_neighbours_names;
}



################################################################
=pod

=item B<get_node_label()>

get the color of the node

=cut
sub get_node_label {
    my ($self, $node_name) = @_;
    my $numId = $self->node_by_name($node_name);
    my %nodes_label = $self->get_attribute("nodes_label");
    my $label = $nodes_label{$numId};
    return $label;
}




################################################################
=pod

=item B<create_node()>

Create and add a node to the graph. 

=cut
sub create_node {
    my ($self, $node_name, $label, $color) = @_;
    my %nodes_name_id = $self->get_attribute("nodes_name_id");
    my %nodes_id_name = $self->get_attribute("nodes_id_name");
    my %nodes_color = $self->get_attribute("nodes_color");
    my %nodes_label = $self->get_attribute("nodes_label");
    my $numId = scalar (keys (%nodes_name_id)); 
    if (!defined($label) || $label eq "") {
      $label = $node_name;
    }
    if (!defined($color) || $color eq "") {
      $color = "#000088";
    }
    $nodes_name_id{$node_name} = $numId;
    $nodes_id_name{$numId} = $node_name;
    $nodes_color{$numId} = $color;
    $nodes_label{$numId} = $label;
    $self->set_hash_attribute("nodes_name_id", %nodes_name_id);
    $self->set_hash_attribute("nodes_id_name", %nodes_id_name);
    $self->set_hash_attribute("nodes_color", %nodes_color);
    $self->set_hash_attribute("nodes_label", %nodes_label);
}

################################################################
=pod

=item B<create_arc()>

Create an arc between to nodes

=cut
sub create_arc {
    my ($self, $source_node_name, $target_node_name, $label, $color) = @_;
    my @out_neighbours = $self->get_attribute("out_neighbours");
    my @in_neighbours = $self->get_attribute("in_neighbours");
    my @arc_in_label = $self->get_attribute("in_label");
    my @arc_out_label = $self->get_attribute("out_label");
    my @arc_in_color = $self->get_attribute("in_color");
    my @arc_out_color = $self->get_attribute("out_color");
    my @arcs = $self->get_attribute("arcs");
    my %arcs_name_id = $self->get_attribute("arcs_name_id");
    my $max_arc_nb = $self->get_attribute("nb_arc_bw_node");
    
    my $numId = scalar (@arcs); 
    if (!defined($label) || $label eq "") {
      $label = $source_node_name."_".$target_node_name;
    }
    if (!defined($color) || $color eq "") {
      $color = "#000044";
    }
    my $source_node_index = $self->node_by_name($source_node_name);
    my $target_node_index = $self->node_by_name($target_node_name);
    my $error =  "Could not create an arc between $source_node_name and $target_node_name\n" ;
    if (defined($source_node_index) && defined($target_node_index)) {
      my $exist = 0;
      my $arc_id = "";
      my $i;
      for ($i = 1; $i <= $max_arc_nb; $i++) {
        $arc_id = $source_node_name."_".$target_node_name."_".$i;
	$exist = exists($arcs_name_id{$arc_id});
      }
      if ($exist) {
	$arc_id = $source_node_name."_".$target_node_name."_".($i);
	$max_arc_nb++;
      } else {
	$arc_id = $source_node_name."_".$target_node_name."_".($i-1);
      }
      $arcs_name_id{$arc_id} = $numId;
      push @{$out_neighbours[$source_node_index]}, $target_node_index;
      push @{$in_neighbours[$target_node_index]}, $source_node_index;
      push @{$arc_out_label[$source_node_index]}, $label;
      push @{$arc_in_label[$target_node_index]}, $label;
      push @{$arc_out_color[$source_node_index]}, $color;
      push @{$arc_in_color[$target_node_index]}, $color;    
      $arcs[$numId][0] = $source_node_name;
      $arcs[$numId][1] = $target_node_name;
      $arcs[$numId][2] = $label;
      $arcs[$numId][3] = $color;
      $self->set_array_attribute("out_neighbours", @out_neighbours);
      $self->set_array_attribute("in_neighbours", @in_neighbours);
      $self->set_array_attribute("in_label", @arc_in_label);
      $self->set_array_attribute("out_label", @arc_out_label);
      $self->set_array_attribute("in_color", @arc_in_color);
      $self->set_array_attribute("out_color", @arc_out_color);
      $self->set_array_attribute("arcs", @arcs);
      $self->set_hash_attribute("arcs_name_id", %arcs_name_id);
      $self->force_attribute("nb_arc_bw_node", $max_arc_nb);  
    } 
    if (!defined($source_node_index)) {
      $error .= "\t$target_node_name not found in the graph\n";
    }
    if (!defined($target_node_index)) {
      $error .= "\t$target_node_name not found in the graph\n";
    }
    if (!defined($target_node_index) || !defined($source_node_index)) {
      &RSAT::message::Warning($error);
    }
}





################################################################
=pod

=item B<get_nodes()>

Return the list of nodes namse of the graph. 

=cut
sub get_nodes {
  my ($self) = @_;
  my %nodes_name_id = $self->get_attribute("nodes_name_id");
  return keys %nodes_name_id;
}

################################################################
=pod

=item B<get_nodes()>

Return a reference to the list of arcs

=cut
sub get_arcs_ref {
  my ($self) = @_;
  my @arcs = $self->get_attribute("arcs");
  my $arcsRef = \@arcs;
  return $arcsRef
}

################################################################
=pod

=item B<get_size()>

Return the number of nodes and arcs of the graph.

Usage: my ($nodes, $arcs) = $graph->get_size();

=cut
sub get_size {
  my ($self) = @_;
  my @nodes = $self->get_nodes();
  my @arcs = $self->get_attribute("arcs");
  my $nodes_nb = scalar(@nodes);
  my $arcs_nb = scalar(@arcs);
  return ($nodes_nb, $arcs_nb);
}


################################################################
=pod

=item B<get_nodes_clusters()>

Return the clusters to which the node specified by its name belongs

=cut
sub get_nodes_clusters {
    my ($self, $node_name) = @_;
    my $numId = $self->node_by_name($node_name);
    my @nodes_clusters = $self->get_attribute("nodes_clusters");
    my @node_clusters = ();
    if (defined($numId) && defined(@{$nodes_clusters[$numId]})) {
      @node_clusters = @{$nodes_clusters[$numId]};
    } else {
      &RSAT::message::Warning("\t","Node",$node_name,"does not belong to any cluster") if ($main::verbose >= 4);;
    }
    return @node_clusters;
}

################################################################
=pod

=item B<node_by_name()>

Returns the id of the node having name $name

=cut
sub node_by_name {
    my ($self, $name) = @_;
    #print "$self";
    my %node_names_id = $self->get_attribute("nodes_name_id");
    if (defined($node_names_id{$name})) {
	return($node_names_id{$name});
    } else {
#	&RSAT::message::Debug("The graph does not contain a node with name ".$name) if ($main::verbose >= 5);
	return();
    }
}

################################################################
=pod

=item B<arcs_by_name()>

Returns the internal ids (in the @arcs array) of the arc having $source as source node and $target as target node

=cut
sub arcs_by_name {
    my ($self, $source, $target) = @_;
    my $max_arc = $self->get_attribute("nb_arc_bw_node");
    my %arcs_name_id = $self->get_attribute("arcs_name_id");
    my @arc_ids = ();
    for (my $i = 1; $i <= $max_arc; $i++) {
      my $arc_id = join ("_", $source, $target, $i);
      my $id = $arcs_name_id{$arc_id};
      if (defined($arcs_name_id{$arc_id})) {
        push @arc_ids, $id;
      }
    }
    return @arc_ids;
}


################################################################
=pod

=item B<node_by_id()>

Returns the name of the node having name $id

=cut
sub node_by_id {
    my ($self, $id) = @_;
    #print "$self";
    my %node_id_names = $self->get_attribute("nodes_id_name");
    if (defined($node_id_names{$id})) {
	return($node_id_names{$id});
    } else {
#	&RSAT::message::Debug("The graph does not contain a node with ID ".$id) if ($main::verbose >= 5);
	return();
    }
}




################################################################
=pod

=item B<read_from_table>

Read the graph from a tab-delimited text file.


 Title    : from_table
 Usage    : $graph->->read_from_table($input_file)
 Function : Read the graph from a tab-delimited text file.
 Returns  : void

=cut

sub read_from_table2 {
    ################################"
    # Define variables
    my ($self, $inputfile, $source_col, $target_col, $weight_col) = @_;
    my @out_neighbours = $self->get_attribute("out_neighbours");
    my @in_neighbours = $self->get_attribute("in_neighbours");
    my @arc_out_label = $self->get_attribute("out_label");
    my @arc_in_label = $self->get_attribute("in_label");
    my @arc_in_color = $self->get_attribute("in_color");
    my @arc_out_color = $self->get_attribute("out_color");
    my @arcs = $self->get_attribute("arcs");
    my %arcs_name_id = $self->get_attribute("arcs_name_id");
 
    
    my %nodes_name_id = $self->get_attribute("nodes_names_id");
    my %nodes_id_name = $self->get_attribute("nodes_id_names");
    my %nodes_color = $self->get_attribute("nodes_color");
    my %nodes_label = $self->get_attribute("nodes_label");
    
    my $max_arc_nb = $self->get_attribute("nb_arc_bw_node");
   
    my $nodecpt = 0;
    my $arccpt = 0;
    &RSAT::message::TimeWarn("Loading graph from tab file", $inputfile) if ($main::verbose >= 2);
    ($main::in) = &RSAT::util::OpenInputFile($inputfile); 
    my $no_weight = 0;
    my $default_weight = 1;
    my $node_default_color = $self->get_attribute("node_color") || "#000088";
    my $arc_default_color = $self->get_attribute("arc_color") || "#000044";

    ## Check input parameters
    unless (&RSAT::util::IsNatural($source_col) && ($source_col > 0)) {
	&RSAT::error::FatalError(join("\t", $source_col, "Invalid source column secification for graph loading. Should be a strictly positive natural number."));
    }
    unless (&RSAT::util::IsNatural($target_col) && ($target_col > 0)) {
	&FatalError(join("\t", $target_col, "Invalid target column specification for graph loading. Should be a strictly positive natural number."));
    }
    unless (&RSAT::util::IsNatural($weight_col) && ($weight_col > 0)) {
	$no_weigth = 1;
    }

    ## Load the graph
    my $l = 0;
    while (<$main::in>) {
	$l++;
	if(($main::verbose >= 2) && ($l % 1000 == 0)) {
	    &RSAT::message::TimeWarn("\tLoaded", $l, "lines from file", $inputfile);
	}
	next if (/^--/); # Skip mysql-like comments
	next if (/^;/); # Skip RSAT comments
	next if (/^#/); # Skip comments and header
	next unless (/\S/); # Skip empty rows
	chomp;
	my @fields = split("\t");
	my $source_name = $fields[$source_col-1];
	my $target_name = $fields[$target_col-1];
	my $weight = $default_weight;
	unless ($no_weight) {
	  if ($weight_col > 0) {
	    $weight = $fields[$weight_col-1];
	  }
	}

	## Source node
	my $source_node_index = $nodes_name_id{$source_name};
	if (!defined($source_node_index)) {
	    my $node_label = $source_name;
	    $source_node_index = $nodecpt;
	    $nodes_color{$source_node_index} = $node_default_color;
	    $nodes_label{$source_node_index} = $node_label;
	    $nodes_name_id{$source_name} = $nodecpt;
	    $nodes_id_name{$nodecpt} = $source_name;
	    $nodecpt++;
	    &RSAT::message::Info(join("\t", "Created source node", 
				      $source_name,
				      $source_node_index, 
				      $node_label)
				     ) if ($main::verbose >= 3);
	}

	## Target node
	my $target_node_index = $nodes_name_id{$target_name};
	if (!defined($target_node_index)) {
	    my $node_label = $target_name;
	    $target_node_index = $nodecpt;
	    $nodes_color{$target_node_index} = $node_default_color;
	    $nodes_label{$target_node_index} = $node_label;
	    $nodes_name_id{$target_name} = $nodecpt;
	    $nodes_id_name{$nodecpt} = $target_name;
	    $nodecpt++;
	    &RSAT::message::Info(join("\t", "Created target node", 
				      $target_name,
				      $target_node_index, 
				      $node_label)
				     ) if ($main::verbose >= 3);
	}
	## Create the arc
	my $arc_label = "";
	if ($no_weight) {
	    $arc_label = join ("_", $source_name, $target_name);
	} else {
	    $arc_label = $weight;
	}
	push @{$out_neighbours[$source_node_index]}, $target_node_index;
	push @{$in_neighbours[$target_node_index]}, $source_node_index;
	push @{$arc_out_label[$source_node_index]}, $arc_label;
	push @{$arc_in_label[$target_node_index]}, $arc_label;
	push @{$arc_out_color[$source_node_index]}, $arc_default_color;
	push @{$arc_in_color[$target_node_index]}, $arc_default_color;
        my $exist = 0;
        my $arc_id = "";
        my $i;
	for ($i = 1; $i <= $max_arc_nb; $i++) {
	  $arc_id = $source_name."_".$target_name."_".$i;
	  $exist = exists($arcs_name_id{$arc_id});
	}
	if ($exist) {
	  $arc_id = $source_name."_".$target_name."_".($i);
	  $max_arc_nb++;
	} else {
	  $arc_id = $source_name."_".$target_name."_".($i-1);
	}
	$arcs_name_id{$arc_id} = $arccpt;
	$arcs[$arccpt][0] = $source_name;
	$arcs[$arccpt][1] = $target_name;
	$arcs[$arccpt][2] = $arc_label;
	$arcs[$arccpt][3] = $arc_default_color;
	$arccpt++;
	&RSAT::message::Info(join("\t", "Created arc", 
				  $source_name, $target_name
				 )) if ($main::verbose >= 4);
    }
    close $main::in if ($inputfile);
    ################################"
    # Save the tables
    $self->set_array_attribute("out_neighbours", @out_neighbours);
    $self->set_array_attribute("in_neighbours", @in_neighbours);
    $self->set_array_attribute("in_label", @arc_in_label);
    $self->set_array_attribute("out_label", @arc_out_label);
    $self->set_array_attribute("in_color", @arc_in_color);
    $self->set_array_attribute("out_color", @arc_out_color);
    $self->set_array_attribute("arcs", @arcs);
    
    $self->set_hash_attribute("nodes_name_id", %nodes_name_id);
    $self->set_hash_attribute("nodes_id_name", %nodes_id_name);
    $self->set_hash_attribute("nodes_color", %nodes_color);
    $self->set_hash_attribute("nodes_label", %nodes_label);
    $self->set_hash_attribute("arcs_name_id", %arcs_name_id);
    
    $self->force_attribute("nb_arc_bw_node", $max_arc_nb);   
}

################################################################

sub read_from_table {
    ################################"
    # Define variables
    my ($self, $inputfile, $source_col, $target_col, $weight_col, $source_color_col, $target_color_col, $edge_color_col) = @_;
    &RSAT::message::TimeWarn("Loading graph from tab file", $inputfile) if ($main::verbose >= 2);
    ($main::in) = &RSAT::util::OpenInputFile($inputfile); 
    my $weight = 0;
    my $default_weight = 1;
    my @array = ();
    my $default_node_color = "#000088";
    my $default_edge_color = "#000044";
    ## Check input parameters
    unless (&RSAT::util::IsNatural($source_col) && ($source_col > 0)) {
	&RSAT::error::FatalError(join("\t", $source_col, "Invalid source column secification for graph loading. Should be a strictly positive natural number."));
    }
    unless (&RSAT::util::IsNatural($target_col) && ($target_col > 0)) {
	&FatalError(join("\t", $target_col, "Invalid target column specification for graph loading. Should be a strictly positive natural number."));
    }
    if (&RSAT::util::IsNatural($weight_col) && ($weight_col > 0)) {
	$weight = 1;
    }

    ## Load the graph
    $cpt = 0;
    while (my $ligne = <$main::in>) {
      next if ($ligne =~ /^--/); # Skip mysql-like comments
      next if ($ligne =~ /^;/); # Skip RSAT comments
      next if ($ligne =~ /^#/); # Skip comments and header
      next unless ($ligne =~ /\S/); # Skip empty rows
      chomp ($ligne);
      my @lignecp = split "\t", $ligne;
      $array[$cpt][0] = $lignecp[$source_col-1];
      $array[$cpt][1] = $lignecp[$target_col-1];
      if ($weight) {
        $array[$cpt][2] = $lignecp[$weight_col-1];
      } else {
        $array[$cpt][2] = join("_",$lignecp[$source_col-1],$lignecp[$target_col-1]);
        
      }
      if (defined($source_color_col)) {
        $array[$cpt][3] = $lignecp[$source_color_col-1] || $default_node_color;
      } else {
        $array[$cpt][3] = $default_node_color;
      }
      if (defined($target_color_col)) {
        $array[$cpt][4] = $lignecp[$target_color_col-1] || $default_node_color;
      } else {
        $array[$cpt][4] = $default_node_color;
      }
      if (defined($edge_color_col)) {
        $array[$cpt][5] = $lignecp[$edge_color_col-1] || $default_edge_color;
      } else {
        $array[$cpt][5] = $default_edge_color;
      }
      $cpt++;
    }
    $self->load_from_array(@array);
    return $self;
}

################################################################

=pod

=item B<remove_duplicated_arcs>

From a graph where arcs may be duplicated, create a graph with only unique arcs.
Usage	:	$graph->remove_duplicated_arcs($directed)

=cut

sub remove_duplicated_arcs {
  my ($self, $directed) = @_;
  if ($main::verbose >= 2) {
    &RSAT::message::TimeWarn("Remove duplicated edges");
  }
  my %seen;
  my %nodes_name_id = $self->get_attribute("nodes_name_id");
  my %nodes_color = $self->get_attribute("nodes_color");
  my @unique_array = ();
  my @arcs = $self->get_attribute("arcs");
  my $arccpt = 0;
  for (my $i = 0; $i < scalar(@arcs); $i++) {
    my @nodes = ();
    push @nodes, $arcs[$i][0], $arcs[$i][1];
    if (!$directed) {
      @nodes = sort(@nodes);
    }
    my $arc_id = $nodes[0]."_".$nodes[1];
    if (!exists($seen{$arc_id})) {
      my $source_id = $nodes_name_id{$arcs[$i][0]};
      my $target_id = $nodes_name_id{$arcs[$i][1]};
      $unique_array[$arccpt][0] = $arcs[$i][0];
      $unique_array[$arccpt][1] = $arcs[$i][1];
      $unique_array[$arccpt][2] = $arcs[$i][2];
      $unique_array[$arccpt][3] = $nodes_color{$source_id};
      $unique_array[$arccpt][4] = $nodes_color{$target_id};
      $unique_array[$arccpt][5] = $arcs[$i][3];
      $arccpt++;
      $seen{$arc_id}++;      
    }
  }
  $self->reload_graph;
  $self->load_from_array(@unique_array);
  return $self;
}


################################################################

=pod

=item B<reload_graph>

empty all tables composing the graph

=cut

sub reload_graph {
  if ($main::verbose >= 2) {
    &RSAT::message::TimeWarn("Reload graph");
  }
  my ($self) = @_;
  my @out_neighbours =();
  my @in_neighbours = ();
  my @arc_in_label = ();
  my @arc_out_label = ();
  my @arc_in_color = ();
  my @arc_out_color = ();
  my @arcs = ();
  my $max_arc = 1;
  my %arcs_name_id = ();
  my %nodes_color = ();
  my %nodes_label = ();
  my %node_name_id = ();
  my %node_id_name = ();
  $self->set_array_attribute("out_neighbours", @out_neighbours);
  $self->set_array_attribute("in_neighbours", @in_neighbours);
  $self->set_array_attribute("in_label", @arc_in_label);
  $self->set_array_attribute("out_label", @arc_out_label);
  $self->set_array_attribute("in_color", @arc_in_color);
  $self->set_array_attribute("out_color", @arc_out_color);
  $self->set_array_attribute("arcs", @arcs);
  $self->force_attribute("nb_arc_bw_node", $max_arc);
  $self->set_hash_attribute("arcs_name_id", %arcs_name_id);
  $self->set_hash_attribute("nodes_name_id", %nodes_name_id);
  $self->set_hash_attribute("nodes_id_name", %nodes_id_name);
  $self->set_hash_attribute("nodes_color", %nodes_color);
  $self->set_hash_attribute("nodes_label", %nodes_label);
  return $self;
}


################################################################



sub load_from_gml {
    ################################
    # Define variables
    my ($self, $inputfile) = @_;
    my @out_neighbours = $self->get_attribute("out_neighbours");
    my @in_neighbours = $self->get_attribute("in_neighbours");
    my @arc_out_label = $self->get_attribute("out_label");
    my @arc_in_label = $self->get_attribute("in_label");
    my @arc_in_color = $self->get_attribute("in_color");
    my @arc_out_color = $self->get_attribute("out_color");
    my @arcs = $self->get_attribute("arcs");
    my %arcs_name_id = $self->get_attribute("arcs_name_id");
 
    
    my %nodes_name_id = $self->get_attribute("nodes_name_id");
    my %nodes_id_name = $self->get_attribute("nodes_id_name");
    my %nodes_color = $self->get_attribute("nodes_color");
    my %nodes_label = $self->get_attribute("nodes_label");
    my %gml_id = ();
    
    my $max_arc_nb = $self->get_attribute("nb_arc_bw_node");
   
    my $nodecpt = 0;
    my $arccpt = 0;
    &RSAT::message::TimeWarn("Loading graph from gml file", $inputfile) if ($main::verbose >= 2);
    ($main::in) = &RSAT::util::OpenInputFile($inputfile); 
    
    my $fichier;

    while (my $ligne = <$main::in>) {
      chomp($ligne);
      $ligne =~ s/\t/ /g;
      $fichier .= $ligne;
      
    }
    
    
    my @fichier_node = split /node( |\t|)/, $fichier;
    my @fichier_edge = split /edge( |\t|)/, $fichier;
  
    foreach my $node (@fichier_node) {
      if ($node ne " ") {
        $node =~ s/edge.*\[.*//;
        my $node_id =  "NA";
        my $node_label = "NA";
        my $node_color  = "#000088";
# NODE ID
        if ($node =~ /id/) {
          my @node_cp = split /.*id /, $node;
          $node_id = $node_cp[1];
          $node_id = substr($node_id,0, index($node_id, " "));
        }
# NODE LABEL
        if ($node =~ /label/) {
          my @label_cp = split /.*label /, $node;
          $node_label = $label_cp[1];
          $node_label = substr($node_label,1, index($node_label, "\" "));
          $node_label =~ s/\"//;
        }
# NODE COLOR
        if ($node =~ /outline/) {
          my @color_cp = split /.*outline /, $node;
          $node_color = $color_cp[1];
          $node_color = substr($node_color,1, index($node_color, "\" "));
          $node_color =~ s/\"//;
        }
	
	
        if ($node_id ne "NA") {
          if ($node_label eq "NA") {
            $node_label = $node_id;
          }
          my $node_index = $nodes_name_id{$node_label};
          if (!defined($node_index)) {
	    $node_index = $nodecpt;
	    $gml_id{$node_id} = $node_index;
	    $nodes_color{$node_index} = $node_color;
	    $nodes_label{$node_index} = $node_label;
	    $nodes_name_id{$node_label} = $nodecpt;
	    $nodes_id_name{$nodecpt} = $node_label;
	    $nodecpt++;
	    &RSAT::message::Info(join("\t", "Created node", 
	     $node_label,
	     $node_index, 
	     $node_label)
	     ) if ($main::verbose >= 3);
	  }
        }
      }
    }

    
    

    foreach my $edge (@fichier_edge) {
      if ($edge ne " ") {
        my $source_edge =  "NA";
        my $target_edge = "NA";
        my $label_edge  = "NA";
        my $color_edge = "#000044";
    
# SOURCE EDGE
        if ($edge =~ /source/) {
          my @source_cp = split /.*source /, $edge;
          $source_edge = $source_cp[1];
          $source_edge = substr($source_edge,0, index($source_edge, " "));
        }
# TARGET EDGE
        if ($edge =~ /source/) {
          my @target_cp = split /.*target /, $edge;
          $target_edge = $target_cp[1];
          $target_edge = substr($target_edge,0, index($target_edge, " "));
        }
# LABEL EDGE
        if ($edge =~ /label/) {
          my @label_cp = split /.*label /, $edge;
          $label_edge = $label_cp[1];
          $label_edge = substr($label_edge,1, index($label_edge, "\" "));
          $label_edge =~ s/\"//;
        }
# COLOR EDGE
        if ($edge =~ /label/) {
          my @color_cp = split /.*fill /, $edge;
          $color_edge = $color_cp[1];
          $color_edge = substr($color_edge,1, index($color_edge, "\" "));
          $color_edge =~ s/\"//;
        }
	
        if ($source_edge ne "NA" || $target_edge ne "NA") {
          $source_node_id = $gml_id{$source_edge};
	  $source_node = $source_edge;
	  $target_node_id = $gml_id{$target_edge};
	  $target_node = $target_edge;
	  $source_name = $nodes_id_name{$source_node_id};
	  $target_name = $nodes_id_name{$target_node_id}; 
	  if ($label_edge eq "NA") {
	    $label_edge = $node_id_name{$source_node_id}."_".$node_id_name{$target_node_id};
	  }
	    $nodes_label = $label_edge;
	    $nodes_color = $color_edge;
	    
	  push @{$out_neighbours[$source_node_id]}, $target_node_id;
	  push @{$in_neighbours[$target_node_id]}, $source_node_id;
	  push @{$arc_out_label[$source_node_id]}, $label_edge;
	  push @{$arc_in_label[$target_node_id]}, $label_edge;
	  push @{$arc_out_color[$source_node_id]}, $color_edge;
	  push @{$arc_in_color[$target_node_id]}, $color_edge;
          my $exist = 0;
          my $arc_id = "";
          my $i;
	  
	  for ($i = 1; $i <= $max_arc_nb; $i++) {
	    $arc_id = $source_node."_".$target_node."_".$i;
	    $exist = exists($arcs_name_id{$arc_id});
	  }
	  if ($exist) {
	    $arc_id = $source_node."_".$target_node."_".($i);
	    $max_arc_nb++;
	  } else {
	    $arc_id = $source_node."_".$target_node."_".($i-1);
	  }
	  $arcs_name_id{$arc_id} = $arccpt;
	  $arcs[$arccpt][0] = $source_name;
	  $arcs[$arccpt][1] = $target_name;
	  $arcs[$arccpt][2] = $label_edge;
	  $arcs[$arccpt][3] = $color_edge;
	  $arccpt++;
	  
	  &RSAT::message::Info(join("\t", "Created arc", 
	  $label_edge,
	  $arccpt)
	  ) if ($main::verbose >= 3);
        }
      }
    }
    $self->set_array_attribute("out_neighbours", @out_neighbours);
    $self->set_array_attribute("in_neighbours", @in_neighbours);
    $self->set_array_attribute("in_label", @arc_in_label);
    $self->set_array_attribute("out_label", @arc_out_label);
    $self->set_array_attribute("in_color", @arc_in_color);
    $self->set_array_attribute("out_color", @arc_out_color);
    $self->set_array_attribute("arcs", @arcs);
    $self->force_attribute("nb_arc_bw_node", $max_arc);
    $self->set_hash_attribute("arcs_name_id", %arcs_name_id);
    $self->set_hash_attribute("nodes_name_id", %nodes_name_id);
    $self->set_hash_attribute("nodes_id_name", %nodes_id_name);
    $self->set_hash_attribute("nodes_color", %nodes_color);
    $self->set_hash_attribute("nodes_label", %nodes_label);
    return $self;
}
    
    


################################################################

=pod

=item B<load_from_array>

Read the graph from a array
  where col1 = source node
  	col2 = target node
	col3 = weight or edge label
	col4 = source color
	col5 = target color
	col6 = arc color


 Title    : load_from_array
 Usage    : $graph->load_from_array(@input_array)
 Function : Read the graph from a array
 Returns  : a graph

=cut

sub load_from_array {
    ################################"
    # Define variables
    my ($self, @array) = @_;
    my @out_neighbours = $self->get_attribute("out_neighbours");
    my @in_neighbours = $self->get_attribute("in_neighbours");
    my @arc_out_label = $self->get_attribute("out_label");
    my @arc_in_label = $self->get_attribute("in_label");
    my @arc_in_color = $self->get_attribute("in_color");
    my @arc_out_color = $self->get_attribute("out_color");
    my @arcs = $self->get_attribute("arcs");
    my %arcs_name_id = $self->get_attribute("arcs_name_id");
 
    
    my %nodes_name_id = $self->get_attribute("nodes_name_id");
    my %nodes_id_name = $self->get_attribute("nodes_id_name");
    my %nodes_color = $self->get_attribute("nodes_color");
    my %nodes_label = $self->get_attribute("nodes_label");
    
    my $max_arc_nb = $self->get_attribute("nb_arc_bw_node");
   
    my $nodecpt = 0;
    my $arccpt = 0;
    ($main::in) = &RSAT::util::OpenInputFile($inputfile); 
    

    ## Load the graph
    for ($l = 0; $l < scalar(@array); $l++){
	if(($main::verbose >= 3) && ($l % 1000 == 0)) {
	    &RSAT::message::TimeWarn("\tLoaded", $l, "edges");
	}
	my $source_name = $array[$l][0];
	my $target_name = $array[$l][1];
	my $weight = $array[$l][2] || $default_weight;
	if (!defined($array[$l][2])) {
	  $no_weight = 0;
	}
	my $source_color = $array[$l][3];
	my $target_color = $array[$l][4];
	my $edge_color = $array[$l][5];
	
	## Source node
	my $source_node_index = $nodes_name_id{$source_name};
	if (!defined($source_node_index)) {
	    my $node_label = $source_name;
	    $source_node_index = $nodecpt;
	    $nodes_color{$source_node_index} = $source_color;
	    $nodes_label{$source_node_index} = $node_label;
	    $nodes_name_id{$source_name} = $nodecpt;
	    $nodes_id_name{$nodecpt} = $source_name;
	    $nodecpt++;
	    &RSAT::message::Info(join("\t", "Created source node", 
				      $source_name,
				      $source_node_index, 
				      $node_label)
				     ) if ($main::verbose >= 3);
	}

	## Target node
	my $target_node_index = $nodes_name_id{$target_name};
	if (!defined($target_node_index)) {
	    my $node_label = $target_name;
	    $target_node_index = $nodecpt;
	    $nodes_color{$target_node_index} = $target_color;
	    $nodes_label{$target_node_index} = $node_label;
	    $nodes_name_id{$target_name} = $nodecpt;
	    $nodes_id_name{$nodecpt} = $target_name;
	    $nodecpt++;
	    &RSAT::message::Info(join("\t", "Created target node", 
				      $target_name,
				      $target_node_index, 
				      $node_label)
				     ) if ($main::verbose >= 3);
	}
	## Create the arc
	my $arc_label = "";
	if ($no_weight) {
	    $arc_label = join ("_", $source_name, $target_name);
	} else {
	    $arc_label = $weight;
	}
	push @{$out_neighbours[$source_node_index]}, $target_node_index;
	push @{$in_neighbours[$target_node_index]}, $source_node_index;
	push @{$arc_out_label[$source_node_index]}, $arc_label;
	push @{$arc_in_label[$target_node_index]}, $arc_label;
	push @{$arc_out_color[$source_node_index]}, $edge_color;
	push @{$arc_in_color[$target_node_index]}, $edge_color;
        my $exist = 0;
        my $arc_id = "";
        my $i;
	for ($i = 1; $i <= $max_arc_nb; $i++) {
	  $arc_id = $source_name."_".$target_name."_".$i;
	  $exist = exists($arcs_name_id{$arc_id});
	}
	if ($exist) {
	  $arc_id = $source_name."_".$target_name."_".($i);
	  $max_arc_nb++;
	} else {
	  $arc_id = $source_name."_".$target_name."_".($i-1);
	}
	$arcs_name_id{$arc_id} = $arccpt;
	$arcs[$arccpt][0] = $source_name;
	$arcs[$arccpt][1] = $target_name;
	$arcs[$arccpt][2] = $arc_label;
	$arcs[$arccpt][3] = $edge_color;
	$arccpt++;
	&RSAT::message::Info(join("\t", "Created arc", 
				  $source_name, $target_name, $arc_label
				 )) if ($main::verbose >= 4);
    }
    close $main::in if ($inputfile);
    ################################"
    # Save the tables
    $self->set_array_attribute("out_neighbours", @out_neighbours);
    $self->set_array_attribute("in_neighbours", @in_neighbours);
    $self->set_array_attribute("in_label", @arc_in_label);
    $self->set_array_attribute("out_label", @arc_out_label);
    $self->set_array_attribute("in_color", @arc_in_color);
    $self->set_array_attribute("out_color", @arc_out_color);
    $self->set_array_attribute("arcs", @arcs);
    
    $self->set_hash_attribute("nodes_name_id", %nodes_name_id);
    $self->set_hash_attribute("nodes_id_name", %nodes_id_name);
    $self->set_hash_attribute("nodes_color", %nodes_color);
    $self->set_hash_attribute("nodes_label", %nodes_label);
    $self->set_hash_attribute("arcs_name_id", %arcs_name_id);
    
    $self->force_attribute("nb_arc_bw_node", $max_arc_nb);   
}

################################################################
=pod

=item B<graph_from_text()>

Supported formats: gml, tab

=cut
sub graph_from_text {
    my ($self, $in_format, @args) = @_;
    if ($in_format eq "gml") {
	return $self->load_from_gml(@args);
    } elsif ($in_format eq "tab") {
	return $self->read_from_table(@args);
    } else {
	&RSAT::error::FatalError(join ("\t", $in_format, "Invalid format"));
    }
}

################################################################



=pod

=item B<to_text()>

Return the graph in various format.

Supported formats: dot, gml, gdl, tab

=cut
sub to_text {
    my ($self, $out_format, @args) = @_;
    if ($out_format eq "dot") {
	return $self->to_dot(@args);
    } elsif ($out_format eq "gml") {
	return $self->to_gml(@args);
    } elsif ($out_format eq "tab") {
	return $self->to_tab(@args);
    } elsif ($out_format eq "node_table") {
	return $self->to_node_table(@args);
    } else {
	&RSAT::error::FatalError(join ("\t", $out_format, "Invalid format for a graph."));
    }
}

################################################################
=pod

=item B<to_dot()>

Return the graph in dot format. This format can be read by the
GraphViz suite.

=cut
sub to_dot {
    my ($self) = @_;    
    my @nodes = $self->get_nodes();
    my @arcs = $self->get_attribute("arcs");
    
    #print "NODES @nodes\n\n";
    my $dot = "graph G {\n";
    $dot .= "overlap=scale;\n";
    $dot .= "size=\"7,10\";\n";

    foreach my $node (@nodes) {
      my $nodeid = $self->node_by_name($node);
      my $color = $self->get_node_color($node);
      my $label = $self->get_node_label($node);
      #print "NODEID"." ".$node." ".$nodeid." ".$color." ".$label."\n";
      $dot .= join("", 
		   "\"", $node, "\"",  
		   " [color=\"$color\"",
		   ",label=\"$label",
		   "\"]",
		   "\n");
    }
    
    for (my $i = 0; $i < scalar(@arcs); $i++) {
      my $source = $arcs[$i][0];
      my $target = $arcs[$i][1];
      my $label = $arcs[$i][2];
      $dot .=  join ("", "\"", $source, "\" -- \"", $target, "\" [label=\"", $label,"\"]", "\n");
    }
    
    $dot .= "}\n";

    return $dot;
}




################################################################
=pod

=item B<to_gml()>

Return the graph in gml format. 

=cut
sub to_gml {
    my ($self) = @_;    
    my $gml = "";
    my %nodes_id_name = $self->get_attribute("nodes_id_name");
    my @out_neighbours = $self->get_attribute("out_neighbours");
    my @out_label = $self->get_attribute("out_label");
    my @out_color = $self->get_attribute("out_color");
    ## Graph description
    my $graph_label = $self->get_attribute("label") || "graph";
    $gml .= "Creator \"RSAT\"\n";
    $gml .= "Version 1.0\n";
    $gml .= "graph\n";
    $gml .= "[\n";
    $gml .= "	label	\"".$graph_label."\"\n";
    $gml .= "	directed	1\n";

    ## Export nodes
    while (($id, $node_name) = each %nodes_id_name) {
	my $label = $self->get_node_label($node_name); ## node label
	my $w = length($label)*10; ## label width
	my $h = 16; ## label height
	my $x = $id*10;
	my $y = $id*20;
	my $box_color = $self->get_node_color($node_name) || "#0000EE"; ## color for the box around the node
	$gml .= "\t"."node\n";
	$gml .= "\t"."[\n";
	$gml .= "\t\t"."id	".$id."\n";
	$gml .= "\t\t"."label	\"".$label."\"\n";
	$gml .= "\t\t"."graphics\n";
	$gml .= "\t\t"."[\n";
	$gml .= "\t\t\t"."x	".$x."\n";
	$gml .= "\t\t\t"."y	".$y."\n";
	$gml .= "\t\t\t"."w	".$w."\n";
	$gml .= "\t\t\t"."h	".$h."\n";
#	$gml .= "\t\t\t"."width	1.00000\n";
	$gml .= "\t\t\t"."fill	\"\#ffffff\"\n";
	$gml .= "\t\t\t"."outline	\"".$box_color."\"\n";
	$gml .= "\t\t\t"."type	\"rectangle\"\n";
	$gml .= "\t\t"."]\n";
	$gml .= "\t"."]\n";
    }
    ## Export arcs
    for (my $i = 0; $i < scalar(@out_neighbours); $i++) {
      if (defined (@{$out_neighbours[$i]})) {
        my $source_id = $i;
        my $nodename = $nodes_id_name{$source_id};
        my @source_out_neighbours = @{$out_neighbours[$i]};
        my @source_out_label = @{$out_label[$i]};
        my @source_out_color = @{$out_color[$i]};

        for (my $j = 0; $j < scalar(@source_out_neighbours); $j++) {
          my $target_id = $source_out_neighbours[$j];
          my $arc_color = $source_out_color[$j];
          my $arc_label = $source_out_label[$j];
          $gml .= "\tedge\n";
	  $gml .= "\t"."[\n";
  	  $gml .= "\t\t"."source\t".$source_id."\n";
	  $gml .= "\t\t"."target\t".$target_id."\n";
	  $gml .= "\t\t"."label\t\"".$arc_label."\"\n" if (defined($arc_label));
	  $gml .= "\t\t"."graphics\n";
	  $gml .= "\t\t"."[\n";
	  $gml .= "\t\t\t"."width\t2\n";
	  $gml .= "\t\t\t"."type\t\"line\"\n";
	  $gml .= "\t\t\t"."fill\t\"".$arc_color."\"\n";
	  $gml .= "\t\t]\n";
	  $gml .= "\t]\n";
        }
      }
    }
    ## Close the graph
    $gml .= "]\n";

    return $gml;
}

################################################################
=pod

=item B<to_tab()>

Return the graph in a tab-delimited format. 

=cut
sub to_tab {
    my ($self) = @_;    
    my @arcs = $self->get_attribute("arcs");
    my @arcs_attributes = $self->get_attribute("arcs_attribute");
    my $tab = join("\t","#source","target","label","color");
    $tab .= "\n";
    if (@arcs_attributes && scalar(@arcs_attributes) > 0) {
      $tab = $self->to_tab_arcs_attribute();
    } else {
      for (my $i = 0; $i < scalar(@arcs); $i++) {
        $tab .= $arcs[$i][0]."\t";
        $tab .= $arcs[$i][1]."\t";
        $tab .= $arcs[$i][2]."\t";
        $tab .= $arcs[$i][3]."\n";
      }
    }
    return $tab;
}

################################################################
=pod

=item B<to_tab_index($arcs_multiple_attribute)>

Return the graph in a tab-delimited format. 
$arcs_multiple_attribute is an Index object, having
source_target_label_color as key

=cut
sub to_tab_arcs_attribute {
    my ($self) = @_;    
    my @arcs = $self->get_attribute("arcs");
    my @arcs_attributes = $self->get_attribute("arcs_attribute");
    my @arcs_attribute_header = $self->get_attribute("arcs_attribute_header");
    my $tab = join("\t","#source", "target", "label", "color");
    if (@arcs_attribute_header) { 
      $tab .= "\t".join("\t",@arcs_attribute_header);
    } else {
      $tab .= "\tattribute";
    }
    $tab .= "\n";
    # if @arcs_attribute_header is not defined or has scalar = 1, then one attribute for each row 
    # else one attribute by tab
    for (my $i = 0; $i < scalar(@arcs); $i++) {
      my $source = $arcs[$i][0];
      my $target = $arcs[$i][1];
      my $label = $arcs[$i][2];
      my $color =  $arcs[$i][3];
      my $attribute = $arcs_attributes[$i];
      if (defined($attribute)) {
        my @clusters = @{$attribute};
	if (@clusters && (!@arcs_attribute_header || scalar(@arcs_attribute_header) == 1)) { 
          foreach my $cluster (@clusters) {
            $tab .= $source."\t";
            $tab .= $target."\t";
            $tab .= $label."\t";
	    $tab .= $color."\t";
            $tab .= $cluster."\n";
	  }
	} elsif (@clusters && scalar(@arcs_attribute_header) >= 1) {
	  $tab .= $source."\t";
          $tab .= $target."\t";
          $tab .= $label."\t";
	  $tab .= $color;
	  foreach my $cluster (@clusters) {
            $tab .= "\t".$cluster;
          }
          $tab .= "\n";
        } else {
            $tab .= $source."\t";
            $tab .= $target."\t";
            $tab .= $label."\t";
	    $tab .= $color."\t";
            $tab .= $attribute."\n";	  
	}
      } 
    }
    return $tab;
}




################################################################
=pod

=item B<to_node_table(@out_fields)>

Returns a tab-delimited string with one row per node, and one column
per attribute. 

Parameters

=over

=item I<@out_fields>

List of attributes to export for each node. These must be
single-valued fields.

=back

=cut
sub to_node_table {
  my ($self, @out_fields) = @_;
  &RSAT::message::TimeWarn("Exporting node tables with fields", join(",", @out_fields)) 
    if ($main::verbose >= 2);

  my $node_table  = "";
  my $null = $main::null || "<NULL>";

  ## Print the header
  $node_table .= '#';
  $node_table .= join ("\t", @out_fields);
  $node_table .= "\n";

  ## Check if there is a selection of nodes to export
  my @to_export = ();
  if (defined($self->{nodes_to_export})) {
    @to_export = $self->get_attribute("nodes_to_export");
  } else {
    @to_export = $self->get_nodes();
  }

  ## Export one row per node
  foreach my $node (@to_export) {
    my @out_values = ();
    foreach my $field (@out_fields) {
      my $new_value = $null;
      my @new_values =  $node->get_attribute($field);
      if (scalar(@new_values) > 0) {
	$new_value = $new_values[0];
      }
      push @out_values, $new_value;
    }
    $node_table .= join ("\t", @out_values);
    $node_table .= "\n";
  }

  return $node_table;
}
################################################################
=pod

=item B<load_classes($class_file)>

Load the $class_file by adding an array containing having as coordinate the internal index of the nodes and as component the class_names.
Class names are stored within the attribute cluster_list of the graph object.

Parameters

=over

=item I<@out_fields>

The class file

=back

=cut

sub load_classes {
  my ($self, $inputfile) = @_;
  &RSAT::message::TimeWarn("Loading class information", $inputfile) if ($main::verbose >= 2);
  ($main::in) = &RSAT::util::OpenInputFile($inputfile);
  my %cluster_list;
  my @nodes_clusters;
  while (<$main::in>) {
	next if (/^--/); # Skip mysql-like comments
	next if (/^;/); # Skip RSAT comments
	next if (/^#/); # Skip comments and header
	next unless (/\S/); # Skip empty rows
	chomp;
	my @fields = split("\t");
	my $node_name = $fields[0];
        my $family_name = $fields[1];
        my $node_index = $self->node_by_name($node_name);
        if (defined($node_index)) {
          push @{$nodes_clusters[$node_index]}, $family_name;
          $cluster_list{$family_name} = 1;
        } else {
          #&RSAT::message::TimeWarn("Node $node_id does not exist in the graph") if ($main::verbose >= 2);
        }
  }
  @cluster_list_array = sort(keys(%cluster_list));
  $self->set_array_attribute("cluster_list", @cluster_list_array);
  $self->set_array_attribute("nodes_clusters", @nodes_clusters);

}
return 1;

__END__

