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

Impementation of basic graph functions. This class allows to create
directed graphs and export them in different formats.

Graph class. 


=head1 METHODS

=cut



################################################################
=pod

=item B<add_node()>

Create a node and add it to the graph. 

=cut
sub create_node {
    my ($self, %args) = @_;    
    my $node = new RSAT::GraphNode(%args);
    &RSAT::message::Info(join ("\t", "Created node", 
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
    while (<$main::in>) {
	next if (/^;/);
	next if (/^#/);
	next unless (/\S/);
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
				     )) if ($main::verbose >= 3);
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
				     )) if ($main::verbose >= 3);
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

Supported formats: dot, gml,gdl

=cut
sub to_text {
    my ($self, $out_format) = @_;
    if ($out_format eq "dot") {
	return $self->to_dot();
    } elsif ($out_format eq "gdl") {
	return $self->to_gdl();
    } elsif ($out_format eq "gml") {
	return $self->to_gml();
    } else {
	&RSAT::error::FatalError(join ("\t", $out_format, "Invlid format for a graph."));
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

Return the graph in gdl (Graph Data Linker) format.  This is an XML
description of the graph object used by the Snow system
(http://www.scmbb.ulb.ac.be/biomaze/).

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
	my $box_color = $node->get_attribute("fontcolor") || "#0000EE"; ## color for the box around the node
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
	my $arc_label = $arc->get_attribute("label");
	my $source_node = $arc->get_attribute("source");
	my $target_node = $arc->get_attribute("target");
	my $source_id = $node_internal_id{$source_node->get_attribute("id")};
	my $target_id = $node_internal_id{$target_node->get_attribute("id")};
	my $arc_color = $arc->get_attribute("color") || "#000000";
	$gml .= "\tedge\n";
	$gml .= "\t"."[\n";
	$gml .= "\t\t"."source\t".$source_id."\n";
	$gml .= "\t\t"."target\t".$target_id."\n";
	$gml .= "\t\t"."label\t\"".$arc_label."\"\n";
	$gml .= "\t\t"."graphics\n";
	$gml .= "\t\t"."[\n";
	$gml .= "\t\t\t"."width\t2\n";
	$gml .= "\t\t\t"."type\t\"line\"\n";
	$gml .= "\t\t\t"."fill\t\"".$arc_color."\"\n";
	$gml .= "\t\t]\n";
	$gml .= "\t]\n";
    }

    ## Close the graph
    $gml .= "]\n";

    
    return $gml;
}


return 1;

__END__

