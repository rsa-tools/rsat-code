###############################################################
#
# Class GraphArc
#
package RSAT::GraphArc;

use RSAT::GenericObject;
use RSAT::message;
use RSAT::error;
@ISA = qw( RSAT::GenericObject );

### class attributes

=pod

=head1 NAME

    RSAT::GraphArc

=head1 DESCRIPTION

GraphArc class. 

=cut


################################################################
=pod

=item new()

Create a new GraphArc.

=cut
sub new {
    my ($class, %args) = @_;
    my $self = bless {
	}, $class;

    $self->init(%args);

    ## check source node
    my $source;
    my $source_id;
    if($args{source}) {
	$source = $self->get_attribute("source");
	unless ($source->isa("RSAT::GraphNode")) {
	    &RSAT::error::FatalError("Source node does not belong to the class RSAT::GraphNode");
	}
	$source_id = $source->get_attribute("id");
    } else {
	&RSAT::error::FatalError( "Cannot create an arc without source node");
    }

    ## check target node
    my $target;
    my $target_id;
    if($args{target}) {
	$target = $self->get_attribute("target");
	unless ($target->isa("RSAT::GraphNode")) {
	    &RSAT::error::FatalError("Target node does not belong to the class RSAT::GraphNode");
	}
	$target_id = $target->get_attribute("id");
    } else {
	&RSAT::error::FatalError( "Cannot create an arc without target node");
    }


    ## Report the creation of the arc
    &RSAT::message::Info(join("\t", "Created arc",
			      "source",
			      $source,
			      $source_id,
			      "target",
			      $target,
			      $target_id,
			     )) if ($main::verbose >= 5);
    return $self;
}

return 1;

__END__

