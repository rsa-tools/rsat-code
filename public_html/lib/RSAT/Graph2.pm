###############################################################
#
# Class Graph2
#
package RSAT::Graph2;

use RSAT::GenericObject;
use RSAT::error;
use RSAT::util;
use RSAT::stats;
use List::Util 'shuffle';
use POSIX qw(log10);

require "RSA.lib";


### class attributes
@ISA = qw( RSAT::GenericObject );
#$node_color="#000088";
#$arc_color="#000044";

=pod

=head1 NAME

    RSAT::Graph2

=head1 DESCRIPTION

Implementation of basic graph functions. This class allows (among other things) to create
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
    
    
    my $mean_weight = "null";
    my $sd_weight = "null";
    my $max_weight = "null";
    my $min_weight = "null";
    my $real = "null";
    
    my %arcs_name_id = ();
    my %nodes_color = ();
    my %nodes_label = ();
    my %nodes_name_id = ();
    my %nodes_id_name = ();
    my %nodes_id_xpos = ();
    my %nodes_id_ypos = ();
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
    # Properties of the weight (may be calculated only if all weights are real)
    $self->force_attribute("mean_weight", $mean_weight);
    $self->force_attribute("mean_sd", $mean_sd);
    $self->force_attribute("max_weight", $max_weight);
    $self->force_attribute("min_weight", $min_weight);
    $self->force_attribute("real", $real);
    # $maxarc represents the number of time that an arc may be duplicated.
    $self->force_attribute("nb_arc_bw_node", $max_arc);
    # %arcs_name_id correspondance between the name of an arc and its internal id in the @arcs array.
    # the key of the hash is source_target_x, with x an integer from 1 to $max_arc;
    $self->set_hash_attribute("arcs_name_id", %arcs_name_id);
   
    # %nodes_name_id correspondance between name and id
    $self->set_hash_attribute("nodes_name_id", %nodes_name_id);
    # %nodes_id_name correspondance between id and name
    $self->set_hash_attribute("nodes_id_name", %nodes_id_name);
    # %nodes_color correspondance between node_name and color
    $self->set_hash_attribute("nodes_color", %nodes_color);
    # %nodes_label correspondance between node_name and label
    $self->set_hash_attribute("nodes_label", %nodes_label);
    # node x and y positions
    $self->set_hash_attribute("nodes_id_xpos", %nodes_id_xpos);
    $self->set_hash_attribute("nodes_id_ypos", %nodes_id_ypos);
    return $self;
}

################################################################
=pod

=item B<empty_graph()>

=cut
sub empty_graph {
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
    @out_neighbours = ();
    @in_neighbours = ();
    @arc_out_label = ();
    @arc_in_label = ();
    @arc_in_color = ();
    @arc_out_color = ();
    @arcs = ();
    %arcs_name_id = ();
    %nodes_name_id = ();
    %nodes_id_name =  ();
    %nodes_color = ();
    %nodes_label =  ();
    %gml_id = ();
    $self->set_array_attribute("out_neighbours", @out_neighbours);
    $self->set_array_attribute("in_neighbours", @in_neighbours);
    $self->set_array_attribute("in_label", @arc_in_label);
    $self->set_array_attribute("out_label", @arc_out_label);
    $self->set_array_attribute("in_color", @arc_in_color);
    $self->set_array_attribute("out_color", @arc_out_color);
    $self->set_array_attribute("arcs", @arcs);
    $self->force_attribute("nb_arc_bw_node", $max_arc_nb);
    $self->set_hash_attribute("arcs_name_id", %arcs_name_id);
    $self->set_hash_attribute("nodes_name_id", %nodes_name_id);
    $self->set_hash_attribute("nodes_id_name", %nodes_id_name);
    $self->set_hash_attribute("nodes_color", %nodes_color);
    $self->set_hash_attribute("nodes_label", %nodes_label);
    return $self;
}

################################################################
=pod

=item B<set_seed_nodes>

Define node filter for the loading. When seed nodes are defined, the
method load_from_table() only accepts loads the induced subgraph,
i.e. arcs for which both source and target nodes belong to the set of
seeds.

Usage:
$graph->set_seed_nodes(@seed_ids);

=cut
sub set_seed_nodes {
  my ($self, @seed_ids) = @_;

  ## Create an empty seed index
  %{$self->{seed_index}} = ();
  foreach my $seed (@seed_ids) {
    $self->{seed_index}->{$seed}++;
  }
  $self->{seed_nodes} = 1;
}

################################################################
=pod

=item B<randomize()> 

returns a randomized graph having the same number of edges, each node having the same number of neighbours

=cut
sub randomize {

    my ($self, $self_loop, $duplicated, $directed, $iter) = @_;
    # get arcs and nodes of $graph
    my @arcs = $self->get_attribute("arcs");
    my %nodes_name_id = $self->get_attribute("nodes_name_id");
    my $nodes_number = scalar keys %nodes_name_id;
    my $req_edges = scalar @arcs;
    my %nodes_label = $self->get_attribute("nodes_label");
    my %nodes_color = $self->get_attribute("nodes_color");
    my $max_arc_number;
    
    ## Computation of the maximum number of edges
    if (!$duplicated) { 
      if (!$directed && !$self_loop) {
        $max_arc_number = ($nodes_number*($nodes_number-1))/2;
      } elsif ($directed && $self_loop) {
        $max_arc_number = ($nodes_number*$nodes_number);
      } elsif ($directed && !$self_loop) {
        $max_arc_number = ($nodes_number*($nodes_number-1));
      } elsif (!$directed && $self_loop) {
        $max_arc_number = ($nodes_number*($nodes_number+1))/2;
      }
      if ($max_arc_number < $req_edges) {
        &RSAT::error::FatalError("\t","More requested edges than possible edges", "requested", $req_edges, "available", scalar $max_arc_number);
      }
    }
    
    # create a new graph object;
    my $rdm_graph = new RSAT::Graph2;
    my @rdm_graph_array = ();
    # graph containing the non relevant edges (self-loop, ...)
    my @rdm_graph_array_err = (); # graph containing the non relevant edges (self-loop, ...)
    my %seen = ();
    my @arcs_to_shuffle = ();
    for (my $i = 0; $i < scalar(@arcs); $i++) {
      $arcs_to_shuffle[$i][0] = $arcs[$i][0];
      $arcs_to_shuffle[$i][1] = $arcs[$i][1];
      $arcs_to_shuffle[$i][2] = $arcs[$i][2];
      $arcs_to_shuffle[$i][3] = $arcs[$i][3];
    }
    my $invalid = 0;
    my $invalid_old = 0;
    my $non_decreasing_count = 0;
    my $loopcount = 0;
    while (scalar(@arcs) != scalar(@rdm_graph_array)) {
      $loopcount++;
      my @arcs_size = 0 .. (scalar(@arcs_to_shuffle)-1);
      my @shuffled_arcs_size = shuffle(@arcs_size);
      for (my $i = 0; $i < scalar(@arcs_to_shuffle); $i++) {
        my $source = $arcs_to_shuffle[$i][0];
        my $target = $arcs_to_shuffle[$shuffled_arcs_size[$i]][1];
        my $label = join("_", $source, $target);
        if (&RSAT::util::IsReal($arcs_to_shuffle[$i][2]) || $main::weight_col) {
          $label = $arcs[$i][2];
        }
        my $arc_color = $arcs[$i][3];
        my $source_id = $nodes_name_id{$source};
        my $target_id = $nodes_name_id{$target};
        my $source_color = $nodes_color{$source};
        my $target_color = $nodes_color{$target};
        my $interaction1 = join("_",$source,$target);
        my $interaction2 = $interaction1;
        if (!$directed) {
          $interaction2 = join("_",$target,$source);
        }

        if ((!$self_loop && ($source_id eq $target_id)) || (!$duplicated && ($seen{$interaction1} || $seen{$interaction2}))) {
          if ($loopcount <= $iter) {
            $rdm_graph_array_err[$invalid][0] = $source;
            $rdm_graph_array_err[$invalid][1] = $target;
            $rdm_graph_array_err[$invalid][2] = $label;
            $rdm_graph_array_err[$invalid][3] = $source_color;
            $rdm_graph_array_err[$invalid][4] = $target_color;
            $rdm_graph_array_err[$invalid][5] = $arc_color;
            $invalid++;
            next;
          } else {
            $invalid++;
          }
        }
        $cpt = scalar(@rdm_graph_array);
        $rdm_graph_array[$cpt][0] = $source;
        $rdm_graph_array[$cpt][1] = $target;
        $rdm_graph_array[$cpt][2] = $label;
        $rdm_graph_array[$cpt][3] = $source_color;
        $rdm_graph_array[$cpt][4] = $target_color;
        $rdm_graph_array[$cpt][5] = $arc_color;
        $seen{$interaction1}++;
        $seen{$interaction2}++;
      }
      if ($loopcount > $iter) {
        last;
      }
      @arcs_to_shuffle = ();
      for (my $j = 0; $j < $invalid; $j++) {
        $arcs_to_shuffle[$j][0] = $rdm_graph_array_err[$j][0];
        $arcs_to_shuffle[$j][1] = $rdm_graph_array_err[$j][1];
        $arcs_to_shuffle[$j][2] = $rdm_graph_array_err[$j][2];
        $arcs_to_shuffle[$j][3] = $rdm_graph_array_err[$j][3];
      }

      my $val = (int (scalar(@arcs_to_shuffle))/20)+1;
      if ($invalid == $invalid_old) {
        ## if the number of invalid edges did not decrease in the 10 
        ## last shuffling, a  60d of the number of the remaining arcs to shuffle ((scalar(@arcs_to_shuffle))/60))
        ## are removed from the already shuffled arcs (@rdm_graph_array) and added to the 
        ## arcs that need to be shuffled again 
        ## these arcs are also removed from the %seen hash
        $non_decreasing_count++;
        if ($non_decreasing_count >= 20 && $val < scalar(@rdm_graph_array)) {
          my @to_remove = @rdm_graph_array[0..$val-1];
          splice(@rdm_graph_array, 0, $val);
          for (my $z = 0; $z < scalar(@to_remove); $z++) {
            my $label1 = join("_", $to_remove[$z][0], $to_remove[$z][1]);

            delete $seen{$label1};
            if (!$directed) {
              my $label2 = join("_", $to_remove[$z][1], $to_remove[$z][0]);
              delete $seen{$label2};
              if (!exists($seen{$label2})) {
              }
            }
            $arcs_to_shuffle[$z+$invalid][0] = $to_remove[$z][0];
            $arcs_to_shuffle[$z+$invalid][1] = $to_remove[$z][1];
            $arcs_to_shuffle[$z+$invalid][2] = $to_remove[$z][2];
            $arcs_to_shuffle[$z+$invalid][3] = $to_remove[$z][3];
          }
          $invalid += $val;
        }
      } else {
        $non_decreasing_count = 0;
      }
      if ($main::verbose >= 3) {
        &RSAT::message::TimeWarn("\t",$invalid, "arcs invalid... Shuffling procedure starts again on those arcs", $loopcount, "iteration") if ($main::verbose >= 3);
      }

      
      $invalid_old = $invalid;
      $invalid = 0;
      @rdm_graph_array_err = ();
      
      
    }
    $rdm_graph->load_from_array(@rdm_graph_array);
    my %rdm_nodes_name_id = $rdm_graph->get_attribute("nodes_name_id");
    
    if ((scalar(keys %nodes_name_id)) != scalar(keys(%rdm_nodes_name_id))) {
      my %rdm_nodes_id_name = $rdm_graph->get_attribute("nodes_id_name");
      my %rdm_nodes_label = $rdm_graph->get_attribute("nodes_label");
      my %rdm_nodes_color = $rdm_graph->get_attribute("nodes_color");
      my $node_id = scalar(keys(%rdm_nodes_name_id));
      foreach my $node_name (keys %nodes_name_id) {
        if (!exists($rdm_nodes_name_id{$node_name})) {
          my $node_color = $nodes_color{$nodes_name_id{$node_name}};
          my $node_label = $nodes_label{$nodes_name_id{$node_name}};
          $rdm_nodes_name_id{$node_name} = $node_id;
          $rdm_nodes_id_name{$node_id} = $node_name;
          $rdm_nodes_label{$node_id} = $node_label;
          $rdm_nodes_color{$node_id} = $node_color;
        }
      }
      $rdm_graph->set_hash_attribute("nodes_id_name", %rdm_nodes_id_name);
      $rdm_graph->set_hash_attribute("nodes_name_id", %rdm_nodes_name_id);
      $rdm_graph->set_hash_attribute("nodes_label", %rdm_nodes_label);
      $rdm_graph->set_hash_attribute("nodes_color", %rdm_nodes_color);
    }
    return ($rdm_graph);
}


################################################################
=pod

=item B<create_random_graph()> 
Usage : $graph->create_random_graph(
      $nodes_ref,
      $req_nodes,
      $req_edges,
      $self_loops,
      $duplicated,
      $directed,
      $max_degree,
      $single,
      $mean,
      $sd,
      $normal,
      $column,
      $weights_ref,
      $source_nodes_ref,
      $target_nodes_ref,
      $node_prefix,
      $edge_prefix);


Create a random graph from the nodes in @nodes having $req_nodes nodes
of maximum degree $max_degree and $req_edges edges.

This supports multi-edges ($duplicated = 1) or not ($duplicated = 0)
and self loops ($self_loops = 1) or not ($self_loops = 0). A weight is
calculated according to the normal distribution if the $mean and $sd
value are given as argument.

@source_nodes and @target_nodes contain the source nodes and the
target nodes of the original graph respectively (if exists), if the
$column boolean is set to 1, then source node will remain source node
and target node will remain target nodes (useful for bipartite
graphs).

=cut
sub create_random_graph {
  my ($self,
      $nodes_ref,
      $req_nodes,
      $req_edges,
      $self_loops,
      $duplicated,
      $directed,
      $max_degree,
      $single,
      $mean,
      $sd,
      $normal,
      $column,
      $weights_ref,
      $source_nodes_ref,
      $target_nodes_ref,
      $node_prefix,
      $edge_prefix) = @_;
  my $rdm_graph = new RSAT::Graph2();
  my @rdm_graph_array = ();
  my $max_arc_number = 10000000;
  my @nodes = @{$nodes_ref};
  my @source_nodes = @{$source_nodes_ref};
  my @target_nodes = @{$target_nodes_ref};
  my $req_source_nodes = $req_nodes;
  my $req_target_nodes = $req_nodes;
  my @labels = @{$weights_ref}; 
  $node_prefix = "" if (!defined($node_prefix));
  $edge_prefix = "" if (!defined($egde_prefix));


  ## creation of the list of nodes
  if (scalar(@nodes) > 0) {
    if (scalar(@nodes) < $req_nodes) {
      &RSAT::error::FatalError("\t","More requested nodes than available nodes", "requested", $req_nodes, "available", scalar @nodes);
    } elsif (!$column) {
      my @indices = 0 .. (scalar(@nodes)-1);
      my @shuffle_indices = &shuffle(@indices);
      my @random_nodes = ();
      for (my $i = 0; $i < $req_nodes; $i++) {
	push @random_nodes, $nodes[$shuffle_indices[$i]];
      }
      @source_nodes = @random_nodes;
      @target_nodes = @random_nodes;
    }
  } else {
    for (my $i = 1; $i <= $req_nodes; $i++) {
      my $node_id = $node_prefix."n".$i;
      push @source_nodes, $node_id;
      push @target_nodes, $node_id;
      push @nodes, $node_id;
    }
  }
  if ($column) {
    $req_source_nodes = scalar(@source_nodes);
    $req_target_nodes = scalar(@target_nodes);
  }
  ## Computation of the maximum number of edges
  if (!$duplicated) { 
    if (!$directed && !$self_loops) {
      $max_arc_number = ($req_nodes*($req_nodes-1))/2;
#       print "REQ NODES".$req_nodes."\n";
    } elsif ($directed && $self_loops) {
      $max_arc_number = ($req_nodes*$req_nodes);
    } elsif ($directed && !$self_loops) {
      $max_arc_number = ($req_nodes*($req_nodes-1));
    } elsif (!$directed && $self_loops) {
      $max_arc_number = ($req_nodes*($req_nodes+1))/2;
    }
  }
  if ($max_arc_number < $req_edges) {
    &RSAT::error::FatalError("\t","More requested edges than possible edges", "requested", $req_edges, "available", $max_arc_number);
  }
  my @possible_source = ();
  my @possible_target = ();
  my %degree;
  my %graph_node;
  my %seen;
  my $k = 0;
  my @random_source = 0 .. ($req_source_nodes-1);
  my @random_target = 0 .. ($req_target_nodes-1);
  @random_source = &shuffle(@random_source);
  for (my $i = 0; $i < $req_source_nodes; $i++) {
    @random_target = &shuffle(@random_target);
    
    for (my $j = 0; $j < $req_target_nodes; $j++) {
      my $source_index = $random_source[$i%($req_source_nodes)];
      my $target_index = $random_target[$j%($req_target_nodes)];
      my $source = $source_nodes[$source_index];
      my $target = $target_nodes[$target_index];
      my $label = join("_", $source, $target);
      my $inv_label = join("_", $target, $source);

      if (($source eq $target) && !$self_loops) {
	next;
      }
      if (exists($seen{$label}) && !$duplicated) {
	next;
      }
      if ((exists($seen{$label}) || exists($seen{$inv_label})) && !$directed && !$duplicated) {
	next;
      } 
      if ((exists($seen{$label}) && exists($seen{$inv_label})) && !$duplicated) {
	next;
      }
      if ($max_degree > 0 && exists($degree{$source}) && exists($degree{$target})) {
	if ($degree{$source} > ($max_degree-1)) {
	  next;
	}
	if ($degree{$target} > ($max_degree-1)) {
	  next;
	}
      }
      $degree{$source}++;
      $degree{$target}++;
      $seen{$label}++;
      push @possible_source, $source;
      push @possible_target, $target;
      $k++;
      if (($k % 100000 == 0) && ($main::verbose >= 3)) {
	&RSAT::message::psWarn("\t","$k" ,"potential edges created.");
      }
    }
    if ($k > (50*$req_edges) || $k > $max_arc_number) {
      last;
    }
  }
  if ($duplicated) {
    @possible_target = &shuffle(@possible_target)
  }
  &RSAT::message::Info("\t",scalar(@possible_target) ,"potential edges created.") if ($main::verbose >= 3);
  ## In case the maximum degree specification prevented to create enough edges : error.
  if ((scalar(@possible_target)) < $req_edges) {
    &RSAT::error::FatalError("\t","Maximal node degree specification is not compatible with the number of requested edges");
  }
  my @random_edges = 0 .. (scalar(@possible_source)-1);
  my @random_weights = 0 .. $req_edges-1;
  @random_edges = &shuffle(@random_edges);
  @random_weights = &shuffle(@random_weights);
  my $weightcpt = 0;
  my $count = 0;
  while (scalar(@rdm_graph_array) < $req_edges) {
    for (my $i = 0; $i < scalar(@random_edges); $i++) {
      if (scalar(@rdm_graph_array) == $req_edges) {
        last;
      }
      my $source = $possible_source[$random_edges[$i]];
      my $target = $possible_target[$random_edges[$i]];
      my $label = $edge_prefix.$source."_".$target;
      my $cpt = scalar(@rdm_graph_array); ## Number of arcs in the table
      if ($mean ne 'null' && $sd ne 'null' && $normal) {
        $label = &gaussian_rand();
        $label = ($label * $sd) + $mean;
      } elsif (!$normal) {
        my $weight = $labels[$random_weights[$cpt]];
        if (&RSAT::util::IsReal($arcs_to_shuffle[$i][2]) || $main::weight_col) {
	  $label = $weight;
        } else {
# 	  &RSAT::message::Warning("Not a numeric weight on edge between", $source, "and", $target) if ($main::verbose >= 3);
        }
      }

      if (($single && $count == 0) && (exists($graph_node{$source}) && exists($graph_node{$target}))) {
        next;
      }
      $graph_node{$source}++;
      $graph_node{$target}++;
      $rdm_graph_array[$cpt][0] = $source;
      $rdm_graph_array[$cpt][1] = $target;
      $rdm_graph_array[$cpt][2] = $label;
      $rdm_graph_array[$cpt][3] = "#000088";
      $rdm_graph_array[$cpt][4] = "#000088";
      $rdm_graph_array[$cpt][5] = "#000044";
      &RSAT::message::Info("\t", "Random edge created between", $source, "and", $target, "with label", $label) if ($main::verbose >= 3);
    }
    $count++;
  }
  $rdm_graph->load_from_array(@rdm_graph_array);
    
  ## Add nodes that have degree 0 (if %graph_node < req_nodes)
  
  
  ## MAY BE ERRONEOUS
  if (scalar(keys(%graph_node)) < $req_nodes) {
    my %rdm_nodes_id_name = $rdm_graph->get_attribute("nodes_id_name");
    my %rdm_nodes_name_id = $rdm_graph->get_attribute("nodes_name_id");
    my %rdm_nodes_label = $rdm_graph->get_attribute("nodes_label");
    my %rdm_nodes_color = $rdm_graph->get_attribute("nodes_color");
    my $node_id = scalar(keys(%rdm_nodes_name_id));
    foreach my $node (@nodes) {
      if (!exists($graph_node{$node})) {
	my $node_color = "#000088";
	my $node_label = $node;
	$rdm_nodes_name_id{$node} = $node_id;
	$rdm_nodes_id_name{$node_id} = $node;
	$rdm_nodes_label{$node_id} = $node_label;
	$rdm_nodes_color{$node_id} = $node_color;
	$node_id++;
      }
    }
    $rdm_graph->set_hash_attribute("nodes_id_name", %rdm_nodes_id_name);
    $rdm_graph->set_hash_attribute("nodes_name_id", %rdm_nodes_name_id);
    $rdm_graph->set_hash_attribute("nodes_label", %rdm_nodes_label);
    $rdm_graph->set_hash_attribute("nodes_color", %rdm_nodes_color);
  }
  return ($rdm_graph);
}



################################################################
=pod

=item B<random_graph_degree_distrib()> 
Usage : $graph->random_graph_degree_distrib();

create a random graph from another graph. The global connectivity distribution being conserved.
This function works by permuting node names.


=cut
sub random_graph_degree_distrib {
    my ($self) = @_;
    my $rdm_graph = new RSAT::Graph2();
    my @rdm_graph_array = ();
    my @arcs = $self->get_attribute("arcs");
    my @new_arcs; 
    my @nodes = $self->get_nodes;
    my @new_nodes = &shuffle(@nodes);
    my %nodes_new_nodes = ();
    my %nodes_name_id = $self->get_attribute("nodes_name_id");
    my %nodes_label = $self->get_attribute("nodes_label");
    my %nodes_color = $self->get_attribute("nodes_color");
    for (my $i = 0; $i < scalar(@nodes); $i++) {
      $node_new_nodes{$nodes[$i]} = $new_nodes[$i]; 
    }
    my %graph_node = ();
    for (my $i = 0; $i < scalar(@arcs); $i++) {
      $rdm_graph_array[$i][0] = $node_new_nodes{$arcs[$i][0]};
      $rdm_graph_array[$i][1] = $node_new_nodes{$arcs[$i][1]};
      my $label = join("_", $rdm_graph_array[$i][0], $rdm_graph_array[$i][1]);
      if (&RSAT::util::IsReal($arcs[$i][2]) || $main::weight_col) {
        $label = $arcs[$i][2];
      }
      my $source_id = $nodes_name_id{$rdm_graph_array[$i][0]};
      my $target_id = $nodes_name_id{$rdm_graph_array[$i][1]};
      my $source_color = $nodes_color{$source_id};
      my $target_color = $nodes_color{$target_id};
      $graph_node{$arcs[$i][0]}++;
      $graph_node{$arcs[$i][1]}++;
      $rdm_graph_array[$i][2] = $label;
      $rdm_graph_array[$i][3] = $source_color;
      $rdm_graph_array[$i][4] = $target_color;
      $rdm_graph_array[$i][5] = "#000044";
    }
    $rdm_graph->load_from_array(@rdm_graph_array);
    ## Add nodes that have degree 0
    if (scalar(keys(%graph_node)) != scalar(keys(%nodes_name_id))) {
      my %rdm_nodes_id_name = $rdm_graph->get_attribute("nodes_id_name");
      my %rdm_nodes_name_id = $rdm_graph->get_attribute("nodes_name_id");
      my %rdm_nodes_label = $rdm_graph->get_attribute("nodes_label");
      my %rdm_nodes_color = $rdm_graph->get_attribute("nodes_color");
      my $node_id = scalar(keys(%rdm_nodes_name_id));
      foreach my $node (keys(%nodes_name_id)) {
        if (!exists($graph_node{$node})) {
          my $node_id = $nodes_name_id{$rdm_graph_array[$i][0]};
          my $node_color = $nodes_color{$node_id};
          my $node_label = $node;
          $rdm_nodes_name_id{$node} = $node_id;
          $rdm_nodes_id_name{$node_id} = $node;
          $rdm_nodes_label{$node_id} = $node_label;
          $rdm_nodes_color{$node_id} = $node_color;
          $node_id++;
        }
      }
      $rdm_graph->set_hash_attribute("nodes_id_name", %rdm_nodes_id_name);
      $rdm_graph->set_hash_attribute("nodes_name_id", %rdm_nodes_name_id);
      $rdm_graph->set_hash_attribute("nodes_label", %rdm_nodes_label);
      $rdm_graph->set_hash_attribute("nodes_color", %rdm_nodes_color);
    }
    return ($rdm_graph);    
    
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

returns an array contening the out neighbours of a node to a certain step self, the node being included or not

1st col : neighbour_node
2nd col : seed_node
3d col  : direction
4d col  : steps
5th col : weight (only if step == 1)


=cut
sub get_neighbours_id {
  my ($self, $node_id, $step, $included, $weight) = @_;
  &RSAT::message::Info("\t","Looking for neighbours of node", $node_id) if $main::verbose > 2;
  my %node_id_name = $self->get_attribute("nodes_id_name");
  my @out_neighbours = $self->get_attribute("out_neighbours");
  my @in_neighbours = $self->get_attribute("in_neighbours");
  my @arc_out_label = $self->get_attribute("out_label");
  my @arc_in_label = $self->get_attribute("in_label");
  my %direction = ();
  my %weight = ();
  my %seen_nodes = ();
  my @result = ();
  if ($included) {
    $result[0][0] = $node_id_name{$node_id};
    $result[0][1] = $node_id_name{$node_id};
    $result[0][2] = 0;
    $result[0][3] = "self";
    if ($weight) {
      $result[0][4] = "NA";
    }

  }
  $seen_nodes{$node_id} = 0;
  my $neighbours_cpt = scalar(@result);
  
  for (my $i = 1; $i <= $step; $i++) {
    my @to_find_nodes = keys (%seen_nodes);
    my @to_add = ();
    ## look for all out neighbours of the nodes that are in %seen_nodes and collect their weight if step = 1
    foreach my $to_find_id (@to_find_nodes) {
      if (defined($out_neighbours[$to_find_id])) {
        push @to_add, @{$out_neighbours[$to_find_id]};
        if ($i == 1) {
          for (my $j = 0; $j < scalar(@{$out_neighbours[$to_find_id]}); $j ++) {
            $direction{$out_neighbours[$to_find_id][$j]} = "out";
            $weight{$out_neighbours[$to_find_id][$j]}  = $arc_out_label[$to_find_id][$j];
          }
        }
      }
      ## look for all in neighbours of the nodes that are in %seen_nodes and collect their weight if step = 1
      if (defined($in_neighbours[$to_find_id])) {
        push @to_add, @{$in_neighbours[$to_find_id]};
        if ($i == 1) {
          for (my $j = 0; $j < scalar(@{$in_neighbours[$to_find_id]}); $j ++) {
            $direction{$in_neighbours[$to_find_id][$j]} = "in";
            $weight{$in_neighbours[$to_find_id][$j]}  = $arc_in_label[$to_find_id][$j];
          }
        }
      }
    }
    foreach my $node_to_add (@to_add) {
      if (!exists($seen_nodes{$node_to_add})) {
        $seen_nodes{$node_to_add} = $step;
        $result[$neighbours_cpt][0] = $node_id_name{$node_to_add};
        $result[$neighbours_cpt][1] = $node_id_name{$node_id};
        $result[$neighbours_cpt][2] = $i;
        $result[$neighbours_cpt][3] = $direction{$node_to_add} || "NA";
        if ($weight) {
          $result[$neighbours_cpt][4] = $weight{$node_to_add};
        }
        $seen_nodes{$node_to_add}++;
        $neighbours_cpt++;
      }
    }
  }

  return @result;
}
################################################################
=pod

=item B<get_clust_coef(@nodes)>

Return the clustering coefficient of a group of nodes.
The clustering coefficient consist in the number of edges among the nodes of the
group divided by the maximum number of edges.
If the graph is directed of may contain self-loop the maximum number of edges is different.
Duplicated edges will be counted only once.

=cut

sub get_clust_coef {
  my ($self, @nodes) = @_;
  my %arcs_name_id = $self->get_attribute("arcs_name_id");
  my $max_arc = $self->get_attribute("nb_arc_bw_node");
  my $group_size = scalar(@nodes);
  my $max_arc_number = (scalar(@nodes)*(scalar(@nodes)-1))/2;
  my $arc_cpt = 0;
  for (my $i = 0; $i < scalar(@nodes); $i++) {
    for (my $j = $i+1; $j < scalar(@nodes); $j++) {
      next if (!$self_loops && ($i == $j));
      my $label = $nodes[$i]."_".$nodes[$j]."_1";
      my $invlabel = $nodes[$j]."_".$nodes[$i]."_1";
      if (exists($arcs_name_id{$label}) || exists($arcs_name_id{$invlabel})) {
        $arc_cpt++;
      }
    }
  }
  my $result = $arc_cpt / $max_arc_number;
  return $result;
}






################################################################
=pod

=item B<properties()>

Return properties of the graph (nodes, edges).

=cut
sub properties {
    my ($self) = @_;
    my @arcs = $self->get_attribute("arcs");
    my $nbarcs = scalar @arcs;
    my %nodes_id_name = $self->get_attribute("nodes_id_name");
    my $nbnodes = scalar (keys(%nodes_id_name));
    return ($nbnodes, $nbarcs);
}
################################################################
=pod

=item B<get_weights()>

Check if the weights are real numbers (or labels)
=cut
sub get_weights() {
  my ($self) = @_;
  my @arcs = $self->get_attribute("arcs");
  my $real = 1;
  for (my $i = 0; $i < scalar(@arcs); $i++) {
    if (!&RSAT::util::IsReal($arcs[$i][2])) {
      $real = 0;
      last;
    } 
  }
  $self->force_attribute("real", $real);
  return $real;
}


################################################################
=pod

=item B<weight_properties()>

Return the properties of the weight distribution

=cut
sub weight_properties {
    my ($self, @weights) = @_;
    my $mean = "null";
    my $sd = "null";
    my $min = "null";
    my $max = "null";
    my $real = $self->get_attribute("real");
    my @arcs = $self->get_attribute("arcs");
    if ($real eq "null") {
      $real = $self->get_weights();
    } 
    if (!@weights) {
      @weights = map $_->[ 2 ], @arcs;
    }
    if ($real) {
      my %summary = &RSAT::stats::summary(@weights);
      $mean = $summary{mean};
      $sd = $summary{sd};
      $min = $summary{min};
      $max = $summary{max};
    } else {
      &RSAT::message::Warning("Cannot compute the mean and standard deviation of the edges : edge weights contain\n\tat least one non real value"."\n") if ($main::verbose >= 5);
    }
    $self->force_attribute("mean_weight", $mean);
    $self->force_attribute("mean_sd", $sd);
    $self->force_attribute("max_weight", $max);
    $self->force_attribute("min_weight", $min);
    return ($mean, $sd, $min, $max);
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
    my @nodes_clusters = $self->get_attribute("nodes_clusters");
    my %node_names_id = $self->get_attribute("nodes_name_id");
    my $numId = $node_names_id{$node_name};
    if (defined($numId) && defined(@{$nodes_clusters[$numId]})) {
      @node_clusters = @{$nodes_clusters[$numId]};
    } else {
      &RSAT::message::Warning("\t","Node",$node_name,"does not belong to any cluster") if ($main::verbose >= 4);;
      @node_clusters = ();
    }
    return @node_clusters;
}


################################################################
=pod

=item B<get_nodes_clusters()>

Return the clusters to which the node specified by its name belongs

=cut
sub get_node_id_clusters {
  my ($self, $numId, @nodes_clusters) = @_;
  @node_clusters = ();
  if (defined(@{$nodes_clusters[$numId]})) {
    @node_clusters = @{$nodes_clusters[$numId]};
  } else {
    &RSAT::message::Warning("\t","Node",$numId,"does not belong to any cluster") if ($main::verbose >= 4);
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
## Load a graph from a tab-delimited text file
sub read_from_table {

    my ($self, $inputfile, $source_col, $target_col, $weight_col,$source_color_col, $target_color_col, $edge_color_col, $source_xpos_col, $source_ypos_col, $target_xpos_col, $target_ypos_col) = @_;

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
    unless (&RSAT::util::IsNatural($source_color_col) && ($source_color_col > 0)) {
      undef $source_color_col; 
    }
    unless (&RSAT::util::IsNatural($target_color_col) && ($target_color_col > 0)) {
      undef $target_color_col;
    }
    unless (&RSAT::util::IsNatural($edge_color_col) && ($edge_color_col > 0)) {
      undef $edge_color_col;
    }
    unless (&RSAT::util::IsNatural($source_xpos_col) && ($source_xpos_col > 0)) {
      undef $source_xpos_col;
    }
    unless (&RSAT::util::IsNatural($source_ypos_col) && ($source_ypos_col > 0)) {
      undef $source_ypos_col;
    }
    unless (&RSAT::util::IsNatural($target_xpos_col) && ($target_xpos_col > 0)) {
      undef $target_xpos_col;
    }
    unless (&RSAT::util::IsNatural($target_ypos_col) && ($target_ypos_col > 0)) {
      undef $target_ypos_col;
    }
    
    ## Load the graph
    $cpt = 0;
    $line_cpt = 0;
    while (my $line = <$main::in>) {
      $linecpt++;
      next if ($line =~ /^--/); # Skip mysql-like comments
      next if ($line =~ /^;/); # Skip RSAT comments
      next if ($line =~ /^#/); # Skip comments and header
      next unless ($line =~ /\S/); # Skip empty rows
      chomp ($line);
      my @linecp = split "\t", $line;
      ## Filter on node names (for induced graphs and graph-get-clusters)
      my $source_id = $linecp[$source_col-1];
      chomp $source_id;
      if ($linecp[$target_col-1]) {
        my $target_id = $linecp[$target_col-1];
        chomp $target_id;
      }
      if ($self->{seed_nodes}) {
	next unless $self->{seed_index}->{$source_id};
	next unless $self->{seed_index}->{$target_id};
      }
      $array[$cpt][0] = $linecp[$source_col-1];
      $array[$cpt][1] = $linecp[$target_col-1];
      ## If there is only a source node, the target node is called ###NANODE###
      ## This term will be recognized by the function load_from_array in order no to create an edge
      ## for this node
      if ((!defined ($array[$cpt][1]) || $array[$cpt][1] eq "")  && $array[$cpt][0] ne "") {
        $array[$cpt][1] = "###NANODE###";
      }
      # Edge weight
      $array[$cpt][2] = join("_",$array[$cpt][0],$array[$cpt][1]);
      if ($weight && $array[$cpt][1] ne "###NANODE###") {
        if (defined ($linecp[$weight_col-1])) {
          $array[$cpt][2] = $linecp[$weight_col-1];
          $array[$cpt][2] =~ s/^\s*//; # remove ending space
         } else {
           &RSAT::message::Warning("No label or weight in column $weight_col on line $linecpt");# if ($main::verbose >= 1);
         }
      }
      # Source Node color
      $array[$cpt][3] = $default_node_color;
      if (defined($source_color_col)) {
        if (defined ($linecp[$source_color_col-1]) && $array[$cpt][1] ne "###NANODE###") { 
        $array[$cpt][3] = $linecp[$source_color_col-1] || $default_node_color;
        } else {
          &RSAT::message::Warning("No source node color in column $source_color_col on line $linecpt");# if ($main::verbose >= 1);
        }
      }
      # Target Node color
      $array[$cpt][4] = $default_node_color;
      if (defined($target_color_col) && $array[$cpt][1] ne "###NANODE###") {
        if (defined ($linecp[$target_color_col-1])) { 
        $array[$cpt][4] = $linecp[$target_color_col-1] || $default_node_color;
        } else {
          &RSAT::message::Warning("No target node color in column $target_color_col on line $linecpt");# if ($main::verbose >= 1);
        }
      }
      # Edge color
      $array[$cpt][5] = $default_edge_color;
      if (defined($edge_color_col) && $array[$cpt][1] ne "###NANODE###") {
        if (defined($linecp[$edge_color_col-1])) {
          $array[$cpt][5] = $linecp[$edge_color_col-1] || $default_edge_color;
        } else {
          &RSAT::message::Warning("No edge color in column $default_edge_color on line $linecpt");# if ($main::verbose >= 1);
        }
      }
      # Source node X position 
      $array[$cpt][6] = undef;
      
      if (defined($source_xpos_col)) {
        if (defined($linecp[$source_xpos_col-1])) {
          $array[$cpt][6] = $linecp[$source_xpos_col-1];
        } else {
          &RSAT::message::Warning("No valid X position for source node in column $source_xpos_col on line $linecpt");# if ($main::verbose >= 1);
        }
      }
      # Source node Y position 
      $array[$cpt][7] = undef;
      if (defined($source_ypos_col)) {
        if (defined($linecp[$source_ypos_col-1])) {
          $array[$cpt][7] = $linecp[$source_ypos_col-1];
        } else {
          &RSAT::message::Warning("No valid Y position for source node in column $source_ypos_col on line $linecpt");# if ($main::verbose >= 1);
        }
      }
      # target node X position 
      $array[$cpt][8] = undef;
      if (defined($target_xpos_col)) {
        if (defined($linecp[$target_xpos_col-1])) {
          $array[$cpt][8] = $linecp[$target_xpos_col-1];
        } else {
          &RSAT::message::Warning("No valid X position for target node in column $target_xpos_col on line $linecpt");# if ($main::verbose >= 1);
        }
      }
      # target node Y position 
      $array[$cpt][9] = undef;
      if (defined($target_ypos_col)) {
        if (defined($linecp[$target_ypos_col-1])) {
          $array[$cpt][9] = $linecp[$target_ypos_col-1];
        } else {
          &RSAT::message::Warning("No valid Y position for target node in column $target_ypos_col on line $linecpt");# if ($main::verbose >= 1);
        }
      }
      $cpt++;
    }
    
#     for (my $i = 0; $i < scalar @array; $i++) {
#       for (my $j = 0; $j <= 9; $j++) {
#         print $array[$i][$j]." ";
#       }
#       print "\n";
#     }
    
    
    
    $self->load_from_array(@array);
    return $self;
}
######################################################################################
## Load a graph from a tab-delimited path file (returned by NeAT PathFinder algorithm)

sub read_from_paths {

    my ($self,$inputfile, $path_col, $distinct_path) = @_;
    $path_col = 7;
    &RSAT::message::TimeWarn("Loading path(s) from tab file", $inputfile) if ($main::verbose >= 2);
    ($main::in) = &RSAT::util::OpenInputFile($inputfile); 
    my $weight = 0;
    my $default_weight = 1;
    my @array = ();
    my $default_node_color = "#000088";
    my $default_edge_color = "#000044";

    ## Check input parameters
    unless (&RSAT::util::IsNatural($path_col) && ($path_col > 0)) {
	&RSAT::error::FatalError(join("\t", $path_col, "Invalid path column secification for paths loading. Should be a strictly positive natural number."));
    }
    
    ## Load the graph
    $cpt = 0;
    $path_cpt = 1;
    my %seen_nodes = ();
    while (my $line = <$main::in>) {
      next if ($line =~ /^--/); # Skip mysql-like comments
      next if ($line =~ /^;/); # Skip RSAT comments
      next if ($line =~ /^#/); # Skip comments and header
      next unless ($line =~ /\S/); # Skip empty rows
      chomp ($line);
      my @linecp = split "\t", $line;
      my $path_arrows = $linecp[$path_col-1];
      &RSAT::error::FatalError("\t","Column $path_col does not contain any path in a valid format") if ($path_arrows !~ /\-\>/);
      my @path = split /\-\>/, $path_arrows;
      my $path_cpt++;
      for (my $i = 0; $i < ((scalar @path) -1); $i++) {
        my $source_node_name = $path[$i];
        my $target_node_name = $path[$i+1];
        $source_node_name = join ("_", $path[$i], $path_cpt) if (defined $seen_nodes{$path[$i]} && $distinct_path);
        $target_node_name = join ("_", $path[$i+1], $path_cpt) if (defined $seen_nodes{$path[$i+1]} && $distinct_path);
        $array[$cpt][0] = $source_node_name;
        $array[$cpt][1] = $target_node_name;
        $array[$cpt][2] = $default_weight;
        $array[$cpt][3] = $default_node_color;
        $array[$cpt][4] = $default_node_color;
        $array[$cpt][5] = $default_edge_color;
        $seen_nodes{$path[$i]}++;
        $seen_nodes{$path[$i+1]}++ if (($i+1) == (scalar (@path) -1));
        $cpt++;
      } 
    }
    $self->load_from_array(@array);
    return $self;
}
################################################################
## Load a graph from a tab delimited adjacency matrix

#   where col1 = source node
#   	col2 = target node
# 	col3 = weight or edge label
# 	col4 = source color
# 	col5 = target color
# 	col6 = arc color


sub read_from_adj_matrix {
    my ($self, $inputfile, $directed) = @_;
    &RSAT::message::TimeWarn("Loading graph from tab delimited adjacency matrix file", $inputfile) if ($main::verbose >= 2);
    ($main::in) = &RSAT::util::OpenInputFile($inputfile); 
    my $weight = 0;
    my $default_weight = 1;
    my @array = ();
    my $default_node_color = "#000088";
    my $default_edge_color = "#000044";
    my %edge_list = ();
    my %nodes_id_name = ();
    my %nodes_name_id = ();


    ## Load the graph
    # find the first line containing the nodes names (starting with a <TAB>)
    my $line = "";
    while ($line = <$main::in> ) {
#       next if ($line !~ /^\t/);
      last if $line =~ /^\t/;
      
    }
    ## if no header -> ERROR
    chomp $line;
    if ($line !~ /^\t/) {
      &RSAT::error::FatalError("\t","Missing header line in input adjacency table",$inputfile);
    }
    my @linecp = split "\t", $line;
    for (my $i = 1; $i < scalar @linecp; $i++) {
      $nodes_id_name{$i-1} = $linecp[$i];
      $nodes_name_id{$linecp[$i]} = $i-1;
    }
    my $cptarray = 0;
    my $cpt = 0;
    while ($line = <$main::in>) {
      next if ($line =~ /^--/); # Skip mysql-like comments
      next if ($line =~ /^;/); # Skip RSAT comments
      next if ($line =~ /^#/); # Skip comments and header
      next unless ($line =~ /\S/); # Skip empty rows
      chomp ($line);
      my @linecp = split "\t", $line;
      if ($nodes_id_name{$cpt} ne $linecp[0]) {
        
        &RSAT::error::FatalError("\t", "Rows and columns are not in the same order");
        
      }
      $cpt++;
      for (my $i = 1; $i < scalar (@linecp); $i++) {
        if ($linecp[$i] ne '0') {
          my $direct_edge_name = join("_", $linecp[0], $nodes_id_name{$i-1});
          my $reverse_edge_name = join("_", $nodes_id_name{$i-1}, $linecp[0]);
          if (!exists($edge_list{$direct_edge_name}) && !exists($edge_list{$reverse_edge_name})) {
            $array[$cptarray][0] = $linecp[0];
            $array[$cptarray][1] = $nodes_id_name{$i-1}; 
            $array[$cptarray][2] = $linecp[$i];
            $array[$cptarray][3] = $default_node_color;
            $array[$cptarray][4] = $default_node_color;
            $array[$cptarray][5] = $default_edge_color;
            if (!$directed) {
              $edge_list{$direct_edge_name}++;
              $edge_list{$reverse_edge_name}++;
            }
            $cptarray++;
          }
        }
      }
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
  my %nodes_id_name = $self->get_attribute("nodes_id_name");
  my %nodes_color = $self->get_attribute("nodes_color");
  my %old_nodes_label = $self->get_attribute("nodes_label");
  my %degree_0_nodes_name_id = $self->get_attribute("nodes_name_id");
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
      delete $degree_0_nodes_name_id{$arcs[$i][0]};
      delete $degree_0_nodes_name_id{$arcs[$i][1]};
    }
  }
  $self->reload_graph;
  $self->load_from_array(@unique_array);
  ## Add nodes that have degree 0
  if (scalar(keys (%degree_0_nodes_name_id)) > 0) {
    my %old_nodes_id_name = %nodes_id_name;
    my %old_nodes_color = %nodes_color;
    my %nodes_id_name = $self->get_attribute("nodes_id_name");
    my %nodes_name_id = $self->get_attribute("nodes_name_id");
    my %nodes_label = $self->get_attribute("nodes_label");
    my %nodes_color = $self->get_attribute("nodes_color");
    my $node_cpt = scalar keys %nodes_id_name;
    while (my ($name, $id) = each (%degree_0_nodes_name_id)) {
      $nodes_name_id{$name} = $node_cpt;
      $nodes_id_name{$node_cpt} = $name;
      $nodes_color{$node_cpt} = $old_nodes_color{$id};
      $nodes_label{$node_cpt} = $old_nodes_label{$id};
      $node_cpt++;
    }
  }

  $self->set_hash_attribute("nodes_id_name", %nodes_id_name);
  $self->set_hash_attribute("nodes_name_id", %nodes_name_id);
  $self->set_hash_attribute("nodes_label", %nodes_label);
  $self->set_hash_attribute("nodes_color", %nodes_color);
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
    my %nodes_id_xpos = $self->get_attribute("nodes_id_xpos");
    my %nodes_id_ypos = $self->get_attribute("nodes_id_ypos");
    my %gml_id = ();
    my $max_arc_nb = $self->get_attribute("nb_arc_bw_node");
    my $nodecpt = 0;
    my $arccpt = 0;
    &RSAT::message::TimeWarn("Loading graph from gml file", $inputfile) if ($main::verbose >= 2);
    ($main::in) = &RSAT::util::OpenInputFile($inputfile); 
    
    my $fichier;

    while (my $line = <$main::in>) {
      chomp($line);
      next if ($line =~ /^#/);
      $line .= " ";
      $line =~ s/\t/ /g;
      $fichier .= $line;
      
    }
    
    
    my @fichier_node = split /node( |\t|)/, $fichier;
    my @fichier_edge = split /edge( |\t|)/, $fichier;

    my %discarded_nodes = (); # if $self->{seed_nodes} is defined, this hash will be filled with the gml_id of the nodes that must no be taken into account when parsing the edges
  
    foreach my $node (@fichier_node) {
      if ($node ne " ") {
        $node =~ s/edge.*\[.*//;
        my $node_id =  "NA#";
        my $node_label = "NA#";
        my $node_color  = "#000088";
        my $node_xpos = "NA#";
        my $node_ypos = "NA#";
# NODE ID
        if ($node =~ /id/) {
          my @node_cp = split /.*id /, $node;
          $node_id = $node_cp[1];
          $node_id = substr($node_id, 0, index($node_id, " "));
        }

# NODE LABEL
        if ($node =~ /label/) {
          my @label_cp = split /.*label /, $node;
          $node_label = $label_cp[1];
          $node_label = substr($node_label, 1, index($node_label, "\" "));
          $node_label =~ s/\"$//;
        }
        if ($self->{seed_nodes} && !$self->{seed_index}->{$node_label}) {
          $discarded_nodes{$node_id}++;
	  next;
        }
# NODE COLOR
        if ($node =~ /outline/) {
          my @color_cp = split /.*outline /, $node;
          $node_color = $color_cp[1];
          $node_color = substr($node_color,1, index($node_color, "\" "));
          $node_color =~ s/\"//;
        }
# NODE X POSITION
        if ($node =~ /x /) {
          my @xpos_cp = split /.*x /, $node;
          $node_xpos = $xpos_cp[1];
          $node_xpos = substr($node_xpos,0, index($node_xpos, " "));
          $node_xpos =~ s/\"//;
         }
# NODE Y POSITION
        if ($node =~ /y /) {
          my @ypos_cp = split /.*y /, $node;
          $node_ypos = $ypos_cp[1];
          $node_ypos = substr($node_ypos,0, index($node_ypos, " "));
          $node_ypos =~ s/\"//;
         }
	
        if ($node_id ne "NA#") {
          if ($node_label eq "NA#") {
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
	    if ($node_xpos ne "NA#") {
	      $nodes_id_xpos{$node_id} = $node_xpos;
	    }
	    if ($node_ypos ne "NA#") {
	      $nodes_id_ypos{$node_id} = $node_ypos;
	    }	    
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
        my $source_edge =  "NA#";
        my $target_edge = "NA#";
        my $label_edge  = "NA#";
        my $color_edge = "#000044";
    
# SOURCE EDGE
        if ($edge =~ /source/) {
          my @source_cp = split /.*source /, $edge;
          $source_edge = $source_cp[1];
          $source_edge = substr($source_edge,0, index($source_edge, " "));
        }
# TARGET EDGE
        if ($edge =~ /target/) {
          my @target_cp = split /.*target /, $edge;
          $target_edge = $target_cp[1];
          $target_edge = substr($target_edge,0, index($target_edge, " "));
        }
# LABEL EDGE
        if ($edge =~ /label/) {
          my @label_cp = split /.*label /, $edge;
          $label_edge = $label_cp[1];
          $label_edge = substr($label_edge,1, index($label_edge, "\" "));
          $label_edge =~ s/\"$//;
        }
# COLOR EDGE
        if ($edge =~ /fill/) {
          my @color_cp = split /.*fill /, $edge;
          $color_edge = $color_cp[1];
          $color_edge = substr($color_edge,1, index($color_edge, "\" "));
          $color_edge =~ s/\"//;
        }
	
	
        if ($source_edge ne "NA#" || $target_edge ne "NA#") {
          next if $discarded_nodes{$source_edge};
          next if $discarded_nodes{$target_edge};
          $source_node_id = $gml_id{$source_edge};
	  $source_node = $source_edge;
	  $target_node_id = $gml_id{$target_edge};
	  $target_node = $target_edge;
	  $source_name = $nodes_id_name{$source_node_id};
	  $target_name = $nodes_id_name{$target_node_id}; 
	  if ($label_edge eq "NA#") {
	    $label_edge = $nodes_id_name{$source_node_id}."_".$nodes_id_name{$target_node_id};
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
    $self->force_attribute("nb_arc_bw_node", $max_arc_nb);
    $self->set_hash_attribute("arcs_name_id", %arcs_name_id);
    $self->set_hash_attribute("nodes_name_id", %nodes_name_id);
    $self->set_hash_attribute("nodes_id_name", %nodes_id_name);
    $self->set_hash_attribute("nodes_color", %nodes_color);
    $self->set_hash_attribute("nodes_label", %nodes_label);
    $self->set_hash_attribute("nodes_id_xpos", %nodes_id_xpos);    
    $self->set_hash_attribute("nodes_id_ypos", %nodes_id_ypos);
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
	col7 = source xpos
	col8 = source ypos
	col9 = target xpos
	col10 = target ypos

The four columns specifying the nodes position must not all be filled. In case the position is not a real value, the value is replaced by a random value.


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
    my ($mean, $sd, $min, $max);
    my $nodecpt = 0;
    my $arccpt = 0;
    
    # if $main::edge_colors (color gradient) or $main::edge_width is defined in convert-graph
    # then, extract the third column of the array of edges 
    # that contains the weight and 
    # that it does not contain not real values.
    # However check first that there is more than one weight.
    # If there is less than 2 edges, $main::edge_colors and $main::edge_width are set to 0
    if (scalar @array < 2 && ($main::edge_colors || $main::edge_width)) {
      $main::edge_colors = 0;
      $main::edge_width = 0;
      &RSAT::message::Warning("The graph has less than 2 edges. -ewidth and -ecolors options will be ignored");
#       print "merde";
    }
    
    
    if ($main::edge_colors || $main::edge_width) {
      $edge_weight_colors = $main::edge_colors;
      my $real = $self->get_weights();
      if ($real) {
        my @weights = map $_->[ 2 ], @array;
        ($mean, $sd, $min, $max) = $self->weight_properties(@weights);
      }
      # if $main:$min_value a,d $main:max_value are defined in convert-graph
      # then use these value as minimum or maximum to compute the edge color gradient
      $min = $main::min_value if (defined $main::min_value && $main::min_value <= $min);
      $max = $main::max_value if (defined $main::max_value && $main::max_value >= $max);      
    }
    

    
    
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
	my $source_xpos = $array[$l][6];
	my $source_ypos = $array[$l][7];
	my $target_xpos = $array[$l][8];
	my $target_ypos = $array[$l][9];
	
	
	## Source node
	my $source_node_index = $nodes_name_id{$source_name};
	if (!defined($source_node_index)) {
	    my $node_label = $source_name;
	    $source_node_index = $nodecpt;
	    $nodes_color{$source_node_index} = $source_color;
	    $nodes_label{$source_node_index} = $node_label;
	    $nodes_name_id{$source_name} = $nodecpt;
	    $nodes_id_name{$nodecpt} = $source_name;
	    
	    
	    
	    $source_xpos = $nodecpt*10 if (!defined($source_xpos) || !&RSAT::util::IsReal($source_xpos));
	    $source_ypos = $nodecpt*10 if (!defined($source_ypos) || !&RSAT::util::IsReal($source_ypos));
	    $nodes_id_xpos{$nodecpt} = $source_xpos;
            $nodes_id_ypos{$nodecpt} = $source_ypos;
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
	    next if ($target_name eq "###NANODE###");
	    my $node_label = $target_name;
	    $target_node_index = $nodecpt;
	    $nodes_color{$target_node_index} = $target_color;
	    $nodes_label{$target_node_index} = $node_label;
	    $nodes_name_id{$target_name} = $nodecpt;
	    $nodes_id_name{$nodecpt} = $target_name;
	    $target_xpos = $nodecpt*10 if (!defined($target_xpos) || !&RSAT::util::IsReal($target_xpos));
	    $target_ypos = $nodecpt*10 if (!defined($target_ypos) || !&RSAT::util::IsReal($target_ypos));
	    $nodes_id_xpos{$nodecpt} = $target_xpos;
            $nodes_id_ypos{$nodecpt} = $target_ypos;
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
	if ($edge_weight_colors) {
	  $edge_color = &RSAT::util::getBgColorFromOneScore($arc_label, $min, $max, 0, $edge_weight_colors);
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
    $self->set_hash_attribute("nodes_id_xpos", %nodes_id_xpos);
    $self->set_hash_attribute("nodes_id_ypos", %nodes_id_ypos);    
    
    
    $self->force_attribute("nb_arc_bw_node", $max_arc_nb);
    

}

################################################################
=pod

=item B<graph_from_text()>

Supported formats: gml, tab, adj_matrix

=cut
sub graph_from_text {
    my ($self, $in_format, @args) = @_;
    if ($in_format eq "gml") {
	return $self->load_from_gml(@args);
    } elsif ($in_format eq "tab") {
        my $inputfile =  $args[0];
        my $scol = $args[1];
        my $tcol = $args[2];
        my $wcol = $args[3] || 0;
        my $sccol = $args[4] || 0;
        my $tccol = $args[5] || 0;
        my $ecol = $args[6] || 0;
        my $sxcol = $args[7] || 0;
        my $sycol = $args[8] || 0;
        my $txcol = $args[9] || 0;
        my $tycol = $args[10] || 0;
        
        
	return $self->read_from_table($inputfile, $scol, $tcol, $wcol, $sccol, $tccol, $ecol, $sxcol, $sycol, $txcol, $tycol);
    } elsif ($in_format eq "adj_matrix") {
	return $self->read_from_adj_matrix($args[0], $args[11]);
    } elsif ($in_format eq "path") {
      my $inputfile = $args[0];
      my $path_col = $args[12];
      my $distinct_path = $args[13] || 0;
      return $self->read_from_paths($inputfile, $path_col, $distinct_path);
    } else {
	&RSAT::error::FatalError(join ("\t", $in_format, "Invalid format"));
    }
}

################################################################



=pod

=item B<to_text()>

Return the graph in various format.

Supported formats: dot, gml, tab, rnsc, tab_numeric

=cut
sub to_text {
    my ($self, $out_format, @args) = @_;
    if ($out_format eq "dot") {
	return $self->to_dot(@args);
    } elsif ($out_format eq "gml") {
        shift @args;
        shift @args;
	return $self->to_gml(@args);
    } elsif ($out_format eq "tab") {
	return $self->to_tab(@args);
    } elsif ($out_format eq "tab_java") {
	return $self->to_tab_java(@args);
    } elsif ($out_format eq "adj_matrix") {
        shift @args;
	return $self->to_adj_matrix(@args);
    } elsif ($out_format eq "rnsc") {
	return $self->to_rnsc(@args);
    } elsif ($out_format eq "tab_numeric") {
	return $self->to_tab_numeric(@args);
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
    my $min = $self->get_attribute("min_weight");
    my $max = $self->get_attribute("max_weight");    
    
    my $edge_width_calc = 0;
    if ($min ne "null" && $max ne "null" && $main::edge_width) {
      $edge_width_calc = 1;
      ## attribution of the minimal and maximal value if specified as arguments
      $min = $min_value if (defined $min_value && $min_value <= $min);
      $max = $max_value if (defined $max_value && $max_value >= $max);
#       print "MIN $min_value $min\n";
    }     
    
    
    #print "NODES @nodes\n\n";
    my $dot = "graph G {\n";
    $dot .= "overlap=scale;\n";
    $dot .= "size=\"36,40\";\n";

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
      my $color = $arcs[$i][3];
      if ($label eq $source."_".$target) {
        $label = "";
      }
      my $edge_width= 2;
      if ($edge_width_calc) {
        $edge_width = ((($label-$min)/($max+1-$min))*6.5)+0.5;
      }
      $dot .=  join ("", "\"", $source, "\" -- \"", $target, "\" [label=\"\", weight=\"$label\", color=\"$color\"]", "\n");
    }
    
    $dot .= "}\n";

    return $dot;
}

################################################################
=pod

=item B<to_adj_matrix()>

Return the graph as an adjacency matrix.

=cut
sub to_adj_matrix {
    my ($self) = shift;
    my $directed = shift;
    my $max_arc = $self->get_attribute("nb_arc_bw_node");
    ## Graph having more than one edge between cannot be
    ## exported in the adjacency matrix format
    if ($max_arc > 1) {
      &RSAT::error::FatalError("\t","This graph is not a simple graph");
    }
    my %nodes_id_name = $self->get_attribute("nodes_id_name");
    my @out_neighbours = $self->get_attribute("out_neighbours");
    my @in_neighbours = $self->get_attribute("in_neighbours");
    my @arc_in_label = $self->get_attribute("in_label");
    my @arc_out_label = $self->get_attribute("out_label");
    my @nodes = $self->get_nodes();
    my $adj_matrix = "";
    my @empty_row = ();
    ## print the header 
    for (my $i = 0; $i < scalar @nodes; $i++) {
      $adj_matrix .= "\t";
      $adj_matrix .= $nodes_id_name{$i};
      push @empty_row, 0;
    }
    $adj_matrix .= "\n";
    ## print the rows
    for (my $i = 0; $i < scalar @nodes; $i++) {
      my @row = @empty_row;
      ## out_neighbours
      if (defined @{$out_neighbours[$i]}) {
        for (my $j = 0; $j < scalar @{$out_neighbours[$i]}; $j++) {
          my $out_id = $out_neighbours[$i][$j];
          my $out_label = $arc_out_label[$i][$j];
          if (&RSAT::util::IsReal($out_label)) {
            $row[$out_id] = $out_label;
          } else {
            $row[$out_id] = 1;
          }
        }
      }
      ## in_neighbours
      if (defined @{$in_neighbours[$i]} && !$directed) {
        for (my $j = 0; $j < scalar @{$in_neighbours[$i]}; $j++) {
          my $in_id = $in_neighbours[$i][$j];
          my $in_label = $arc_in_label[$i][$j];
          if (&RSAT::util::IsReal($in_label)) {
            $row[$in_id] = $in_label;
          } else {
            $row[$in_id] = 1;
          }
        }
      }      
      $adj_matrix .= $nodes_id_name{$i}."\t".join("\t", @row)."\n";
    }
    return $adj_matrix;
}

################################################################
=pod

=item B<get_position()>

Fill the hashes %node_id_xpos and %node_id_ypos using the program
$RSAT/bin/fr_layout. This program based on the boost graph library, computes
the node position using the Fruchterman and Reingold algorithm.
If the program is not available, no changes are brought to the hashes



=cut

sub get_position {
  my ($self) = @_;
  my %nodes_id_xpos = $self->get_attribute("nodes_id_xpos");
  my %nodes_id_ypos = $self->get_attribute("nodes_id_ypos");
  my %nodes_name_id = $self->get_attribute("nodes_name_id");
  my @arcs = $self->get_attribute("arcs");
  my $fr_layout = $ENV{RSAT}."/bin/fr_layout";
  my $nodes_nb = scalar keys %nodes_name_id;
  my $layout_size = 1000;
  my $outfile = $main::outfile{output};

  my $dir = ".";
  if ($outfile) {
    my ($outdir, $outfile) = &RSAT::util::SplitFileName($outfile);
    $dir = $outdir;
    $dir = "." if ($dir eq "");
  }
  $dir .= "/";
  if ($nodes_nb > 6000) {
    $layout_size = 3000;
  } elsif ($nodes_nb < 30) {
    $layout_size = 1000;
  } else {
    $layout_size = int(1000+($nodes_nb/(2)))
  }
  if (!-e $fr_layout) {
    &RSAT::message::Warning("Layout calculator program $fr_layout missing\nLayout will not be computed");
  } else {
    ## The string $arcs_list will be submitted to fr_layout
    ## This C++ script must be located in $RSAT/bin
    ## It takes as argument a list of edges under the form 
    ## source_node1 target_node1
    ## source_node2 target_node2
    ## ... ...
    ## As we may work with large graphs, this function creates a temporary file
    ## where the edges are stored.
    ## This file is then removed.
    
    my $arcs_list = "";
    my $tempfilecmd = "mktemp $dir"."graph_fr.XXXXXXXXXX";
    my $tempfile = `$tempfilecmd`; 
    chomp $tempfile;
    open (TMP, ">$tempfile");
    for (my $i = 0; $i < scalar(@arcs); $i++) {
      my $source = $arcs[$i][0];
      my $target = $arcs[$i][1];
      $source =~ s/ /####space####/g;
      $target =~ s/ /####space####/g;
      $arcs_list .= join(" ",$source, $target, "\n");
    }

    print TMP $arcs_list;
    my $command = "cat $tempfile | $fr_layout --iterations 400 $layout_size $layout_size";
    &RSAT::message::TimeWarn("Calculating the layout with $fr_layout") if ($main::verbose >= 2);
    my $coordinates = `$command`;
    my @lignes = split /\n/, $coordinates;
    system ("rm $tempfile");
    foreach my $ligne (@lignes) {
      ## The four first line of the output of fr_layout displays the iterations
      ## of the program and so must no be parsed, so we skip them.
      next unless ($ligne =~  /\S/);
      next if ($ligne =~ /^0\%\ /); # Skip mysql-like comments
      next if ($ligne =~ /^\|\-/); # Skip RSAT comments
      next if ($ligne =~ /^\*\*/); # Skip comments and header
      my @lignecp = split "\t", "$ligne";
      my $node = $lignecp[0];
      $node =~ s/####space####/ /g; 
      my $id = $nodes_name_id{$node};
      $nodes_id_xpos{$id} = int($lignecp[1]+($layout_size/2));
      $nodes_id_ypos{$id} = int($lignecp[2]+($layout_size/2));
    }
  }
  $self->set_hash_attribute("nodes_id_xpos", %nodes_id_xpos);
  $self->set_hash_attribute("nodes_id_ypos", %nodes_id_ypos);

}

# ################################################################
# =pod
# 
# =item B<get_position()>
# 
# Fill the hashes %node_id_xpos and %node_id_ypos using a spring embedding method.
# 
# 
# 
# =cut
# 
# sub get_position {
#   my ($self) = @_;
#   my %nodes_id_xpos = $self->get_attribute("nodes_id_xpos");
#   my %nodes_id_ypos = $self->get_attribute("nodes_id_ypos");
#   my %nodes_name_id = $self->get_attribute("nodes_name_id");
#   my %arc_name_id = $self->get_attribute("arcs_name_id");
#   my %nodes_id_force = ();
#   my @arcs = $self->get_attribute("arcs");
#   my $nodes_nb = scalar keys %nodes_name_id;
#   my $layout_size = 1000;
#   my $attraction_1 = 2;
#   my $attraction_2 = 1;
#   my $repulsion = 2;
#   my $node_mobility = 0.1;
#   my $max_it = 100;
#   if ($nodes_nb > 6000) {
#     $layout_size = 3000;
#   } elsif ($nodes_nb < 30) {
#     $layout_size = 1000;
#   } else {
#     $layout_size = int(1000+($nodes_nb/(2)))
#   }
#   my @random_position = 0..$layout_size;
#   @random_position = shuffle(@random_position);
#   my $i = 0;
#   # randomize node positions
#   while (($node_name, $id) = each %nodes_name_id) {
#     $nodes_id_xpos{$id} = $random_position[$i];
#     $nodes_id_ypos{$id} = $random_position[($layout_size-$i)];
#     $i++;
#   }
#   my @nodes_names = keys %nodes_name_id;
#   my $mean_mov = 30;
#   for (my $iteration = 0; $iteration < $max_it; $iteration++) {
# #     if ($mean_mov < 17) {
# #       last;
# #     }
#     $mean_mov = 0;
#     my %nodes_force_x = ();
#     my %nodes_force_y = ();
#     for (my $i = 0; $i < scalar(@nodes_names); $i++) {
#       @nodeA_pos = ($nodes_id_xpos{$nodes_name_id{$nodes_names[$i]}}, $nodes_id_ypos{$nodes_name_id{$nodes_names[$i]}});
#       my %already_calculated;
#       my $force_x = 0;
#       my $force_y = 0;
#       for (my $j = 0; $j < scalar(@nodes_names); $j++) {
#         my $arc_id = join("_", $nodes_names[$i], $nodes_names[$j], "1");
#         my $arc_id_rev = join("_", $nodes_names[$j], $nodes_names[$i], "1");
#         if ($i != $j) {
#           @nodeB_pos = ($nodes_id_xpos{$nodes_name_id{$nodes_names[$j]}}, $nodes_id_ypos{$nodes_name_id{$nodes_names[$j]}});
#           # distance calculation between nodeA and nodeB
#           $x_dir = $nodeA_pos[0] - $nodeB_pos[0];
#           print "XDIR $x_dir ";
#           $y_dir = $nodeA_pos[1] - $nodeB_pos[1];
#           $dist_x = abs($x_dir);
#           $dist_y = abs($y_dir);
#           if ($x_dir < 0) {
#             $x_dir = +1;
#           } else {
#             $x_dir = -1;
#           }
#           if ($y_dir < 0) {
#             $y_dir = +1;
#           } else {
#             $y_dir = -1;
#           }          
# 
#           print $iteration." ".$nodes_names[$i]." ".$nodes_names[$j]." P".join(" ", @nodeA_pos)." ".join(" ", @nodeB_pos)." D $dist_x $dist_y";
#           if (defined($arc_name_id{$arc_id}) || defined($arc_name_id{$arc_id_rev})) {
#             # Both nodes are neighbours
#             print " CONNECTED ";
#             $force_x += $x_dir*($attraction_1*log10($dist_x/$attraction_2));
#             $force_y += $y_dir*($attraction_1*log10($dist_y/$attraction_2));
#             print " F ".$x_dir*($attraction_1*log10($dist_x/$attraction_2))." ".$y_dir*($attraction_1*log10($dist_y/$attraction_2))."\n";
#           } else {
#             print " NOT CONNECTED ";
#             # Both nodes are not connected
# #             print "DISTX ".$dist_x."\n";
#             $force_x += -1*($x_dir*($repulsion/sqrt($dist_x)));
# #             print "DISTY ".$dist_y."\n";
#             $force_y += -1*($y_dir*($repulsion/sqrt($dist_y)));
#             print " F".-1*($x_dir*($repulsion/sqrt($dist_x)))." ".-1*($y_dir*($repulsion/sqrt($dist_y)))."\n";
#           }
# #           $already_calculated{$arc_id}++;
# #           $already_calculated{$arc_id_rev}++;
#             
#         
#         }
#       }
#       $nodes_force_x{$nodes_names[$i]} = $force_x;
#       $nodes_force_y{$nodes_names[$i]} = $force_y;
# 
#     }
#     for (my $k = 0; $k < scalar(@nodes_names); $k++ ) {
#       $nodes_id_xpos{$nodes_name_id{$nodes_names[$k]}} += $node_mobility*$nodes_force_x{$nodes_names[$k]};
#       $nodes_id_ypos{$nodes_name_id{$nodes_names[$k]}} += $node_mobility*$nodes_force_y{$nodes_names[$k]};
#       $mean_mov += ($node_mobility*$nodes_force_y{$nodes_names[$k]})/scalar(@nodes_names)^2;
# #       print $c4*$nodes_force_y{$nodes_names[$k]}."\n";
#     }
# #     print "MEAN MOV $mean_mov\n";
# #     print "ITER $iteration\n";
#   }
#   $self->set_hash_attribute("nodes_id_xpos", %nodes_id_xpos);
#   $self->set_hash_attribute("nodes_id_ypos", %nodes_id_ypos);
# }
# 



################################################################
=pod

=item B<to_gml()>

Return the graph in gml format. 

=cut
sub to_gml {
    my ($self, $edge_width,$min_value, $max_value) = @_;
    my $gml = "";
    my %nodes_id_xpos = $self->get_attribute("nodes_id_xpos");
    my %nodes_id_ypos = $self->get_attribute("nodes_id_ypos");
    my %nodes_id_name = $self->get_attribute("nodes_id_name");
    my %nodes_name_id = $self->get_attribute("nodes_name_id");
    my %nodes_color = $self->get_attribute("nodes_color");
    my %nodes_label = $self->get_attribute("nodes_label");
    my @out_neighbours = $self->get_attribute("out_neighbours");
    my @out_label = $self->get_attribute("out_label");
    my @out_color = $self->get_attribute("out_color");
    my @arcs = $self->get_attribute("arcs");
    my $min = $self->get_attribute("min_weight");
    my $max = $self->get_attribute("max_weight");
    my $edge_width_calc = 0;
    if ($min ne "null" && $max ne "null" && $main::edge_width) {
      $edge_width_calc = 1;
      ## attribution of the minimal and maximal value if specified as arguments
      $min = $min_value if (defined $min_value && $min_value <= $min);
      $max = $max_value if (defined $max_value && $max_value >= $max);
#       print "MIN $min_value $min\n";
    } 
    
    

    ## Graph description
    my $graph_label = $self->get_attribute("label") || "graph";
    $gml .= "Creator \"RSAT\"\n";
    $gml .= "Version 1.0\n";
    $gml .= "graph\n";
    $gml .= "[\n";
    $gml .= "	label	\"".$graph_label."\"\n";
    $gml .= "	directed	1\n";
    
    ## Export nodes
    &RSAT::message::Info("Exporting nodes") if $main::verbose >= 3;
    while (($id, $node_name) = each %nodes_id_name) {
        my $label = $nodes_label{$id};
	my $w = length($label)*10; ## label width
	my $h = 16; ## label height
	my $x = $nodes_id_xpos{$id} || $id*10;
	my $y = $nodes_id_ypos{$id} || $id*10;
	my $box_color = $nodes_color{$id} || "#0000EE"; ## color for the box around the node
	
	
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
# 	print "noeud $label\n";
    }
    ## Export edges
    &RSAT::message::Info("Exporting edges") if ($main::verbose >= 3);
    for (my $j = 0; $j < scalar(@arcs); $j++) {
      my $source_name = $arcs[$j][0];
      my $target_name = $arcs[$j][1];
      my $edge_label = $arcs[$j][2];
      my $edge_color = $arcs[$j][3];
      my $source_id = $nodes_name_id{$source_name};
      my $target_id = $nodes_name_id{$target_name};
      $edge_label =~ s/^\s*// if (defined ($edge_label));
      my $edge_width= 2;
      if ($edge_width_calc) {
        $edge_width = ((($edge_label-$min)/($max+1-$min))*6.5)+0.5;
      }
      $gml .= "\tedge\n";
      $gml .= "\t"."[\n";
      $gml .= "\t\t"."source\t".$source_id."\n";
      $gml .= "\t\t"."target\t".$target_id."\n";
      $gml .= "\t\t"."label\t\"".$edge_label."\"\n" if (defined($edge_label));
      $gml .= "\t\t"."graphics\n";
      $gml .= "\t\t"."[\n";
      $gml .= "\t\t\t"."width\t$edge_width\n";
      $gml .= "\t\t\t"."type\t\"line\"\n";
      $gml .= "\t\t\t"."fill\t\"".$edge_color."\"\n";
      $gml .= "\t\t]\n";
      $gml .= "\t]\n";      
    }
    
    
    
    
    ## Close the graph
    $gml .= "]\n";

    return $gml;
}

################################################################
=pod

=item B<to_tab($arc_id)>

Return the graph in a tab-delimited format. 

=cut
sub to_tab {
    
    my ($self, $arc_id) = @_;
    my @arcs = $self->get_attribute("arcs");
    my @arcs_attributes = $self->get_attribute("arcs_attribute");
    my %nodes_name_id = $self->get_attribute("nodes_name_id");
    my $tab = "";
    if (@arcs_attributes && scalar(@arcs_attributes) > 0) {
     $tab = $self->to_tab_arcs_attribute($arc_id);
    } else {
      if (!$arc_id) {
        $tab = join("\t","#source","target","label","color");
        $tab .= "\n";
        for (my $i = 0; $i < scalar(@arcs); $i++) {
          $tab .= $arcs[$i][0]."\t";
          $tab .= $arcs[$i][1]."\t";
          $tab .= $arcs[$i][2]."\t";
          $tab .= $arcs[$i][3]."\n";
          delete $nodes_name_id{$arcs[$i][0]};
          delete $nodes_name_id{$arcs[$i][1]};
        }
        $tab .= join("\n", keys %nodes_name_id)."\n";
      } else {
        my %arcs_name_id = $self->get_attribute("arcs_name_id");
        my $max_arc = $self->get_attribute("nb_arc_bw_node");
        foreach my $arc_name (keys %arcs_name_id) {
          my $id = $arcs_name_id{$arc_name};
          if ($max_arc == 1) {
            $arc_name =~ s/(_1)$//;
          }
          $tab .= $arcs[$id][0]."\t";
          $tab .= $arcs[$id][1]."\t";
          $tab .= $arcs[$id][2]."\t";
          $tab .= $arcs[$id][3]."\t";
          $tab .= $arc_name."\n";
          delete $nodes_name_id{$arcs[$i][0]};
          delete $nodes_name_id{$arcs[$i][1]};
        }
        $tab .= join("\n", keys %nodes_name_id)."\n";
      }
    }
    return $tab;
}
################################################################
################################################################
=pod

=item B<to_rnsc()>

Return the graph in the format used by RNSC. This will consist of one string having separator
'#####NODES_NAME#####'. This string has then to be cut into two files 
- adjacency table of id
- list of id and label association

=cut
sub to_rnsc {
    
    my ($self) = shift;
    my @out_neighbours = $self->get_attribute("out_neighbours");
    my @in_neighbours = $self->get_attribute("in_neighbours");
    my %nodes_name_id = $self->get_attribute("nodes_name_id");
    my $rnsc = "";
    # print the adjacency table of label
    for (my $i = 0; $i < scalar (keys %nodes_name_id); $i++) {
      my @i_neighbours = ();
      @i_neighbours = (@i_neighbours, @{$out_neighbours[$i]}) if (defined $out_neighbours[$i]);
      @i_neighbours = (@i_neighbours, @{$in_neighbours[$i]}) if (defined $in_neighbours[$i]);
      @i_neighbours_sorted = sort {$a <=> $b} @i_neighbours;
      @i_neighbours = @i_neighbours_sorted;
      my $j = 0;
      for ($j = 0; $j < scalar (@i_neighbours); $j++) {
        last if ($i_neighbours[$j] > $i);
      }
      splice (@i_neighbours, 0, $j);
      $rnsc .= $i." ".join (" ",@i_neighbours). " -1\n";
    }
    $rnsc .= "#####NODES_NAME#####";
    while (($name, $id) = each %nodes_name_id) {
     $rnsc .= join("\t", $id, $name)."\n";
    } 
    return $rnsc;
}

################################################################
################################################################
=pod

=item B<to_tab_numeric()>

Return a two columns string. This will consist of one string having separator
'#####NODES_NAME#####'. Before this separator, the nodes names are replaced by their numerical id (two-column tab delimited graph file). In the second part, the string consists in a list of id and their label association.

=cut
sub to_tab_numeric {
    
    my ($self) = shift;
    my @arcs = $self->get_attribute("arcs");
    my %nodes_name_id = $self->get_attribute("nodes_name_id");
    my $output = "";
    for (my $i = 0; $i < scalar @arcs; $i++) {
      $output .= join "\t", (($nodes_name_id{$arcs[$i][0]})+1), (($nodes_name_id{$arcs[$i][1]})+1)."\n";
    }
    $output .= "#####NODES_NAME#####";
    while (($name, $id) = each %nodes_name_id) {
     $output .= join("\t", ($id+1), $name)."\n";
    } 
    return $output;
}

################################################################
=pod

=item B<to_tab_index($arcs_multiple_attribute)>

Return the graph in a tab-delimited format. 
$arcs_multiple_attribute is an Index object, having
source_target_label_color as key

=cut
sub to_tab_arcs_attribute {
    my ($self, $arc_id) = @_;    
    my @arcs = $self->get_attribute("arcs");
    my @arcs_attributes = $self->get_attribute("arcs_attribute");
    my @arcs_attribute_header = $self->get_attribute("arcs_attribute_header");
    my %nodes_name_id = $self->get_attribute("nodes_name_id");
    my $tab = "";
    if (!$arc_id) { 
      $tab = join("\t","#source", "target", "label", "color");
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
	  } elsif (@clusters && (scalar(@arcs_attribute_header) > 1)) {
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
        delete $nodes_name_id{$source};
        delete $nodes_name_id{$target};
      }
    } else {
      $tab = join("\t","#source", "target", "label", "color", "arc_id");
      if (@arcs_attribute_header) { 
        $tab .= "\t".join("\t",@arcs_attribute_header);
     } else {
        $tab .= "\tattribute";
      }
      $tab .= "\n";
      # if @arcs_attribute_header is not defined or has scalar = 1, then one attribute for each row 
      # else one attribute by tab
      my %arcs_name_id = $self->get_attribute("arcs_name_id");
      my $max_arc = $self->get_attribute("nb_arc_bw_node");
        foreach my $arc_name (keys %arcs_name_id) {
        my $id = $arcs_name_id{$arc_name};
        if ($max_arc == 1) {
          $arc_name =~ s/(_1)$//;
        }
        my $source = $arcs[$id][0];
        my $target = $arcs[$id][1];
        my $label = $arcs[$id][2];
        my $color =  $arcs[$id][3];
        my $attribute = $arcs_attributes[$id];
        if (defined($attribute)) {
          my @clusters = @{$attribute};
	  if (@clusters && (!@arcs_attribute_header || scalar(@arcs_attribute_header) == 1)) { 
            foreach my $cluster (@clusters) {
              $tab .= $source."\t";
              $tab .= $target."\t";
              $tab .= $label."\t";
	      $tab .= $color."\t";
	      $tab .= $arc_name."\t";
              $tab .= $cluster."\n";
	    }
	  } elsif (@clusters && scalar(@arcs_attribute_header) >= 1) {
	    $tab .= $source."\t";
            $tab .= $target."\t";
            $tab .= $label."\t";
	    $tab .= $color."\t";
	    $tab .= $arc_name;
	    foreach my $cluster (@clusters) {
              $tab .= "\t".$cluster;
            }
            $tab .= "\n";
          } else {
              $tab .= $source."\t";
              $tab .= $target."\t";
              $tab .= $label."\t";
	      $tab .= $color."\t";
	      $tab .= $arc_name."\t";
              $tab .= $attribute."\n";	  
	  }
        }
        delete $nodes_name_id{$source};
        delete $nodes_name_id{$target};
      }
    }
    $tab .= join("\n", keys %nodes_name_id)."\n";
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

Load the $class_file by adding an array having as
coordinate the internal index of the nodes and as component the
class_names. Class names are stored within the attribute cluster_list
of the graph object.
Nodes belonging to the clusters are placed within the attribute cluster_node_list of
the graph object.

Parameters

=over

=item I<@out_fields>

The class file

=back

=cut
sub load_classes {
  my ($self, $inputfile, $inducted) = @_;
  &RSAT::message::TimeWarn("Loading class information", $inputfile) if ($main::verbose >= 2);
  ($main::in) = &RSAT::util::OpenInputFile($inputfile);
  my %cluster_list;
#   my %cluster_node_list;
  my @nodes_clusters;
  while (<$main::in>) {
	next if (/^--/); # Skip mysql-like comments
	next if (/^;/); # Skip RSAT comments
	next if (/^#/); # Skip comments and header
	next unless (/\S/); # Skip empty rows
	chomp;
	my @fields = split("\t");
	my $node_name = $fields[0]; 
	$family_name = "graph";
	if (!$inducted) {
          $family_name = $fields[1];
        }
        my $node_index = $self->node_by_name($node_name);
        my $line = $node_name."\t".$family_name;
        if (defined($node_index)) {
          if (!exists($lines{$line})) {
            push @{$nodes_clusters[$node_index]}, $family_name;
            $cluster_list{$family_name} = 1;
            $lines{$line}++;
          } else {
            &RSAT::message::Warning("Node\t$node_name\talready in class\t$family_name\t.Skipped") if ($main::verbose >= 3);
          }
        } else {
          &RSAT::message::Warning("Node\t$node_name\tdoes not exist in the graph") if ($main::verbose >= 3);
        }
#         $cluster_node_list{$node_name}++;
        
  }
  @cluster_list_array = sort(keys(%cluster_list));
  $self->set_array_attribute("cluster_list", @cluster_list_array);
  $self->set_array_attribute("nodes_clusters", @nodes_clusters);
#   $self->set_array_attribute("cluster_node_list", keys %cluster_node_list);
}

################################################################
# Returns an array of  normally distributed random numbers
#
# The gaussian_rand function implements the polar Box Muller method
# for turning two independent uniformly distributed random numbers
# between 0 and 1 (such as rand returns) into two numbers with a mean
# of 0 and a standard deviation of 1 (i.e., a Gaussian
# distribution). To generate numbers with a different mean and
# standard deviation, multiply the output of gaussian_rand by the new
# standard deviation, and then add the new mean
sub gaussian_rand {
    my ($u1, $u2);  # uniformly distributed random numbers
    my $w;          # variance, then a weight
    my ($g1, $g2);  # gaussian-distributed numbers

    do {
        $u1 = 2 * rand() - 1;
        $u2 = 2 * rand() - 1;
        $w = $u1*$u1 + $u2*$u2;
    } while ( $w >= 1 );

    $w = sqrt( (-2 * log($w))  / $w );
    $g2 = $u1 * $w;
    $g1 = $u2 * $w;
    # return both if wanted, else just one
    return wantarray ? ($g1, $g2) : $g1;
}

return 1;

#############################################################
################################################################
=pod

=item B<c_topology($arc_id)>

# Calculates the closeness and the betweenness of all nodes.
# It first runs a C program ($RSAT/bin/floydwarshall) to calculate 
# the shortest path between all nodes in a network


=cut 

sub c_topology {
  my ($self, $directed) = (shift, shift);
  my $floydwarshall = $ENV{RSAT}."/bin/floydwarshall";
#   my $floydwarshall = $ENV{RSAT}."/bin/fw-test";
  my %nodes_id_name = $self->get_attribute("nodes_id_name");
  
  my @out_neighbours = $self->get_attribute("out_neighbours");
  my @in_neighbours = $self->get_attribute("in_neighbours");
  my @out_labels = $self->get_attribute("out_label");
  my @in_labels = $self->get_attribute("in_label");
  my %betweenness = ();
  my %distances = ();
  my %closeness = ();
  
  my $real = $self->get_weights();
  
#   print join("\n",%nodes_id_name);
  
  my $outfile = $main::outfile{output};
  my $dir = ".";
  if ($outfile) {
    my ($outdir, $outfile) = &RSAT::util::SplitFileName($outfile);
    $dir = $outdir;
    $dir = "." if ($dir eq "");
  }
  $dir .= "/";
  my $fw_input_file_cmd = "mktemp $dir"."/floydwarshall.input.XXXXX";
  my $fw_output_file_cmd = "mktemp $dir"."/floydwarshall.output.XXXXX";
  my $fw_input_file = `$fw_input_file_cmd`;
  my $fw_output_file = `$fw_output_file_cmd`;
  chomp $fw_input_file;
  chomp $fw_output_file;
  my @empty_array = ();
  for ($i = 0; $i < scalar (keys (%nodes_id_name)); $i++) {
    push @empty_array, 0;
  }
  # Write the adjacency matrix in the format required by $ENV{RSAT}."/bin/floydwarshall
  open FWINPUTFILE, "> $fw_input_file";
  # add the number of nodes
  print FWINPUTFILE scalar (keys %nodes_id_name)."\n";
  # specifies whether the graph is directed or not
  print FWINPUTFILE $directed."\n" if ($floydwarshall !~ /ndir/);
  my $node_nb = scalar (keys %nodes_id_name);
  for (my $i = 0; $i < $node_nb; $i++) {
    &RSAT::message::Info("Computing topology for node", $i."/".$node_nb) if ($main::verbose >= 3);
    @weight_list = @empty_array;
    @i_neighbours = ();
    @i_neighbours = (@i_neighbours, @{$out_neighbours[$i]}) if (defined $out_neighbours[$i]);
#     print "SOURCE $i OUT :";
    for (my $j = 0; $j < scalar @i_neighbours; $j++) {
#       print "$i_neighbours[$j] ";
      my $neighbour = $i_neighbours[$j];
      $label = 1;
      $label = $out_labels[$i][$j] if ($real);
      $weight_list[$neighbour] = $label;
      
    }
#     print "\n###lala####\n";
    if (!$directed) {
      @i_neighbours = ();
      @i_neighbours = (@i_neighbours, @{$in_neighbours[$i]}) if (defined $in_neighbours[$i]);
      for (my $j = 0; $j < scalar @i_neighbours; $j++) {
        my $neighbour = $i_neighbours[$j];
        $label = 1;
        $label = $in_labels[$i][$j] if ($real);
        $weight_list[$neighbour] = $label;
      }
    }
    $weight_list[$i] = 0;
    print FWINPUTFILE join (" ",@weight_list)."\n";
  }
  # Create and launch the command
  my $fw_command = "$floydwarshall $fw_input_file > $fw_output_file";
  &RSAT::message::Info("\t", "Computing all shortest paths") if ($main::verbose >= 2);
  system ("$fw_command");
  # Read the file computed via the C algorithm
  open (FWINPUTFILE, $fw_output_file);
  my $pathnb = 0;
  my %nodes = ();
  &RSAT::message::Info("\t", "Computing all closeness and betweenness") if ($main::verbose >= 2);

  while (my $ligne = <FWINPUTFILE>) {
    next if ($ligne =~ /No route/);
    my @lignecp = split / /, $ligne;
    my $source = $nodes_id_name{$lignecp[0]};
    my $target = $nodes_id_name{$lignecp[1]};
    my $path = $lignecp[2];
    my $dist = $lignecp[3]; 
    my @pathcp = split /,/, $path;
    for (my $i = 0; $i < scalar @pathcp; $i++) {
      my $inter_node = $nodes_id_name{$pathcp[$i]};
      $betweenness{$inter_node}++;
      $nodes{$inter_node}++;
    }
    ${$distances{$source}}[0] += $dist;
    ${$distances{$source}}[1] ++;
    if (!$directed) {
      ${$distances{$target}}[0] += $dist;
      ${$distances{$target}}[1] ++;
    }
    $pathnb++;
    $nodes{$source}++;
    $nodes{$target}++;
  }
  foreach my $node (keys %nodes) {
    if (defined($betweenness{$node})) {
      $betweenness{$node} = ($betweenness{$node}/$pathnb);
    } 
    $closeness{$node} = (${$distances{$node}}[0] / ${$distances{$node}}[1])  if defined($distances{$node});
  }
#   system ("rm $fw_input_file $fw_output_file");
  return (\%betweenness, \%closeness);
}




################################################################
=pod

=item B<to_tab_java($arc_id)>

# Return the graph in a tab-delimited format fully compatible with 
# Neat Java tools (pathfinder ...)

=cut
sub to_tab_java {
    
    my ($self) = @_;
    my @arcs = $self->get_attribute("arcs");
    my %nodes_name_id = $self->get_attribute("nodes_name_id");
    my %nodes_id_name = $self->get_attribute("nodes_id_name");
    my %nodes_color = $self->get_attribute("nodes_color");
    my $tab = "";
    my $real = $self->get_attribute("real");
    if ($real eq "null") {
      $real = $self->get_weights();
    }
    
    #### Nodes definition
    # Nodes list ($tab)+ color ($tab_color)
    $tab_color = ";NODES\trgb_color\n";
    while (my ($id, $color) = each(%nodes_color)) {
      my $name = $nodes_id_name{$id};
      $tab_color .= $name."\t".$color."\n";
    }
    $tab .= $tab_color;
    $tab .= ";ARCS\trgb_color\tlabel\n";
    # Edge list ($tab color + weight)
    for (my $i = 0; $i < scalar(@arcs); $i++) {
      my $source = $arcs[$i][0];
      my $target = $arcs[$i][1];
      my $label = $arcs[$i][2];
      my $color = $arcs[$i][3];
      if (!$real) {
        $label = 1;
      }
      $tab .= join("\t", $source, $target, $color, $label)."\n";
    }
    return $tab;
}

################################################################


__END__
