###############################################################
#
# Class Graph
#
package RSAT::Graph;

use RSAT::GenericObject;
use RSAT::error;
use RSAT::GraphNode;
use RSAT::GraphArc;
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

=item B<create_node()>

Create a node and add it to the graph.

=cut
sub create_node {
    my ($self, %args) = @_;
    my $node = new RSAT::GraphNode(%args);
    &RSAT::message::Info(join ("\t", "Created node", 
			       $node->get_attribute("id"), 
			       $node->get_attribute("label"), 
			       $node)) if ($main::verbose >= 4);
    $self->add_node($node);
    return $node;
}

################################################################
=pod

=item B<add_node()>

Add a (previously created) node to the graph.

=cut
sub add_node {
    my ($self, $node) = @_;    
    &RSAT::message::Info(join ("\t", "Added node", 
			       $node->get_attribute("id"), 
			       $node->get_attribute("label"), 
			       $node)) if ($main::verbose >= 4);
    $self->push_attribute("nodes", $node);
    $self->add_hash_attribute("nodes_by_id", $node->get_attribute("id"), $node);
    return $node;
}

################################################################
=pod

=item B<get_nodes()>

Return the list of nodes of the graph. 

=cut
sub get_nodes {
  my ($self) = @_;
  return $self->get_attribute("nodes");
}

################################################################
=pod

=item B<get_size()>

Return the number of nodes and arcs of the graph.

Usage: my ($nodes, $arcs) = $graph->get_size();

=cut
sub get_size {
  my ($self) = @_;
  my $nodes = scalar(@{$self->{nodes}});
  my $arcs = scalar(@{$self->{arcs}});
  return ($nodes, $arcs);
}



################################################################
=pod

=item B<node_by_id()>

Returns a node specified by its ID. 

=cut
sub node_by_id {
    my ($self, $node_id) = @_;    
    my %node_index = $self->get_attribute("nodes_by_id");
    if (defined($node_index{$node_id})) {
	return($node_index{$node_id});
    } else {
#	&RSAT::message::Debug("The graph does not contain a node with ID ".$node_id) if ($main::verbose >= 5);
	return();
    }
    return $node;
}

################################################################
=pod

=item B<add_arc()>

Create an arc between two existing nodes, and add the new 
graph.

=cut
sub create_arc {
    my ($self, $source_node, $target_node, %args) = @_;    
    my $arc = new RSAT::GraphArc(source=>$source_node, target=>$target_node, %args);
    &RSAT::message::Info(join ("\t", "Created arc from", $source_node->get_attribute("id"),
			       "to", $source_node->get_attribute("id"),
			       $arc)) if ($main::verbose >= 4);
    $self->push_attribute("arcs", $arc);

    ## Add the arc to the source and target nodes
    $source_node->push_attribute("out_arcs", $arc);
    $target_node->push_attribute("in_arcs", $arc);

    ## Cross-linking between source and target nodes
    $source_node->push_attribute("out_nodes", $target_node);
    $target_node->push_attribute("in_nodes", $source_node);

    return $arc;
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

sub read_from_table {
    my ($self, $inputfile, $source_col, $target_col, $weight_col) = @_;
    &RSAT::message::TimeWarn("Loading graph from tab file", $inputfile) if ($main::verbose >= 2);
    ($main::in) = &RSAT::util::OpenInputFile($inputfile); 
    my $no_weight = 0;
    my $default_weight = 1;
    my $node_color = $self->get_attribute("node_color") || "#000088";
    my $arc_color = $self->get_attribute("arc_color") || "#000044";

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
	my $source_node = $self->node_by_id($source_name);
	unless ($source_node) {
	    $source_node = $self->create_node(id=>$source_name, 
					      label=>$source_name,
					      color=>$node_color,
					     );
	    &RSAT::message::Info(join("\t", "Created source node", 
				      $source_name,
				      "id=".$source_node->get_attribute("id"), 
				      "label=".$source_node->get_attribute("label"),
				     )) if ($main::verbose >= 5);
	}

	## Target node
	my $target_node = $self->node_by_id($target_name);
	unless ($target_node) {
	    $target_node = $self->create_node(id=>$target_name, 
					      label=>$target_name,
					      color=>$node_color,
					     );
	    &RSAT::message::Info(join("\t", "Created target node", 
				      $target_name,
				      "id=".$target_node->get_attribute("id"), 
				      "label=".$target_node->get_attribute("label"),
				     )) if ($main::verbose >= 5);
	}

	## Create the arc
	my $arc_label = "";
	if ($no_weight) {
	    $arc_label = join ("_", $source_name, $target_name);
	} else {
	    $arc_label = $weight;
	}
	my $arc = $self->create_arc($source_node, 
				    $target_node, 
				    color=>$arc_color,
				    label=>$arc_label,
				    weight=>$weight);
	&RSAT::message::Info(join("\t", "Created arc", 
				  $source_node->get_attribute("id"), 
				  $target_node->get_attribute("id"), 
				  $arc,
				  $arc->get_attribute("label"),
				  $arc->get_attribute("weight"),
				 )) if ($main::verbose >= 4);
    }
    close $main::in if ($inputfile);    
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
    } elsif ($out_format eq "gdl") {
	return $self->to_gdl(@args);
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
    my $dot = "graph G {;\n";
    $dot .= "overlap=scale;\n";
    $dot .= "size=\"7,10\";\n";

    foreach my $node ($self->get_attribute("nodes")) {
      $dot .= join("", 
		   "\"", $node->get_attribute("id"), "\"",  
		   " [color=",$node->get_attribute("color"),
		   ",fontsize=",$node->get_attribute("fontsize"),
		   ",fontcolor=", $node->get_attribute("fontcolor"),
		   ",shape=",$node->get_attribute("shape"),
		   ",label=\"", $node->get_attribute("label"),
		   "\"]",
		   "\n");
    }

    foreach my $arc ($self->get_attribute("arcs")) {
	my $source = $arc->get_attribute("source");
	my $source_id = $source->get_attribute("id");
	my $target = $arc->get_attribute("target");
	my $target_id = $target->get_attribute("id");
	my $label = $arc->get_attribute("label") || "";
	$dot .=  join ("", "\"", $source_id, "\" -- \"", $target_id, "\" [label=\"", $label,"\"]", "\n");
    }
    
    $dot .= "}\n";

    return $dot;
}

################################################################
=pod

=item B<to_gdl()>

Return the graph in gdl (Graph Data Linker) format.  This was an XML
description of the graph object used by the Snow system. It is now
obsolete, since the Snow team has stopped activities in 2008.

=cut
sub to_gdl {
  require Scmbb::Snow::Graph;
  require Scmbb::Snow::GraphDataLinker;
  my ($self) = @_;
  my $gdl = "";
  &RSAT::message::TimeWarn("Exporting graph in GDL format") if ($main::verbose >= 0); 
  return $gdl;
}


################################################################
=pod

=item B<to_gml()>

Return the graph in gml format. 

=cut
sub to_gml {
    my ($self) = @_;
    my $gml = "";

    ## Graph description
    my $graph_label = $self->get_attribute("label") || "graph";
    $gml .= "Creator \"RSAT\"\n";
    $gml .= "Version 1.0\n";
    $gml .= "graph\n";
    $gml .= "[\n";
    $gml .= "	label	\"".$graph_label."\"\n";
    $gml .= "	directed	1\n";

    ## Export nodes
    my $n = 0;
    my %node = ();
    foreach my $node ($self->get_attribute("nodes")) {
	$n++; 
	$node_internal_id = $n;
#	$node_internal_id = '"'.$node->get_attribute("id").'"';
	$node_internal_id{$node->get_attribute("id")} = $node_internal_id;
	my $label = $node->get_attribute("label"); ## node label
	my $w = length($label)*10; ## label width
	my $h = 16; ## label height
	my $x = $node->get_attribute("x") || $n*10;
	my $y = $node->get_attribute("y") || $n*20;
	my $box_color = $node->get_attribute("color") || "#0000EE"; ## color for the box around the node
	$gml .= "\t"."node\n";
	$gml .= "\t"."[\n";
	$gml .= "\t\t"."id	".$node_internal_id."\n";
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
    foreach my $arc ($self->get_attribute("arcs")) {
	my $source_node = $arc->get_attribute("source");
	my $target_node = $arc->get_attribute("target");
	my $source_id = $node_internal_id{$source_node->get_attribute("id")};
	my $target_id = $node_internal_id{$target_node->get_attribute("id")};

	my $arc_label = $arc->get_attribute("label");
	unless ($arc_label) {
	  if ($weight = $arc->get_attribute("weight")) {
	    $arc_label = $weight;
	  }
	}

	my $arc_color = $arc->get_attribute("color") || "#000000";
	$gml .= "\tedge\n";
	$gml .= "\t"."[\n";
	$gml .= "\t\t"."source\t".$source_id."\n";
	$gml .= "\t\t"."target\t".$target_id."\n";
	$gml .= "\t\t"."label\t\"".$arc_label."\"\n" if ($arc_label);
	$gml .= "\t\t"."graphics\n";
	$gml .= "\t\t"."[\n";
	my $width = 2;
	$gml .= "\t\t\t"."width\t".$width."\n";
	$gml .= "\t\t\t"."type\t\"line\"\n";
	$gml .= "\t\t\t"."fill\t\"".$arc_color."\"\n";
	$gml .= "\t\t]\n";
	$gml .= "\t]\n";
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
    my $tab = "";
    

    ## Graph description
    my $graph_label = $self->get_attribute("label") || "graph";
    $tab .= ";source\ttarget\n";

    ## Export arcs
    foreach my $arc ($self->get_attribute("arcs")) {
	my $source_node = $arc->get_attribute("source");
	my $target_node = $arc->get_attribute("target");
	my $arc_label = $arc->get_attribute("label");
	unless ($arc_label) {
	  if ($weight = $arc->get_attribute("weight")) {
	    $arc_label = $weight;
	  }
	}
	my $arc_color = $arc->get_attribute("color") || "#000000";
        my $sourceid = $source_node->get_attribute("label");
        my $targetid = $target_node->get_attribute("label");
	$tab .= "$sourceid\t$targetid";
        if ($arc_label) {
          $tab .= "\t$arc_label\n";
        } else {
          $tab .= "\n";
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

Load the $class_file by adding to each node the cluster(s) to which it belongs 
Clusters are stored within the attribute cluster_list of the graph object.

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
  while (<$main::in>) {
	next if (/^--/); # Skip mysql-like comments
	next if (/^;/); # Skip RSAT comments
	next if (/^#/); # Skip comments and header
	next unless (/\S/); # Skip empty rows
	chomp;
	my @fields = split("\t");
	my $node_id = $fields[0];
        my $family_name = $fields[1];
        my $node = $self->node_by_id($node_id);
        if ($node) {
          $node->push_attribute("clusters", $family_name);
          $cluster_list{$family_name} = 1;
        } else {
          #&RSAT::message::TimeWarn("Node $node_id does not exist in the graph") if ($main::verbose >= 2);
        }
  }
  @cluster_list = sort(keys(%cluster_list));
  $self->set_array_attribute("cluster_list", @cluster_list);

}
return 1;

__END__

