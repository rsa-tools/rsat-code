###############################################################
#
# Class Graph
#
package RSAT::Graph;

use RSAT::GenericObject;
use RSAT::error;
use RSAT::GraphNode;
use RSAT::GraphArc;

### class attributes
@ISA = qw( RSAT::GenericObject );


=pod

=head1 NAME

    RSAT::Graph

=head1 DESCRIPTION

Impementation of basic graph functions. This class allows to create
directed graphs and export them in different formats.

Graph class. 

=cut



################################################################
=pod

=item add_node()

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
    return $node;
}

################################################################
=pod

=item add_arc()

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
    return $arc;
}

################################################################
=pod

=item to_dot()

Return the graph indot format.

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
	my $label = $arc->get_attribute("label");
	$dot .=  join ("", "\"", $source_id, "\" -- \"", $target_id, "\" [label=\"", $label,"\"]", "\n");
    }
    
    $dot .= "}\n";

    return $dot;
}

################################################################
## Export the graph in GML format
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

