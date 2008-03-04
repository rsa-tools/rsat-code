###############################################################
#
# Class GenomeFeature
#

## CVS: implemented the method matches

package RSAT::GenomeFeature;

use RSAT::GenericObject;
use RSAT::error;
@ISA = qw( RSAT::GenericObject );

### class attributes

=pod

=head1 NAME

    RSAT::GenomeFeature

=head1 DESCRIPTION

GenomeFeature class, for manipulating genome features (CDS, mRNA,
tRNA, ...).

=cut

################################################################
=pod

=item new()

Create a new GenomeFeature.

=cut
sub new {
    my ($class, %args) = @_;
    my $self = bless {
	}, $class;
    return $self;
}



################################################################
=pod

=item matches

Match a string pattern against the ID, names and (optionnally)
description of the GenomicFeature. Return 1 if a match is found, 0
otherwise.

Usage: $feature->matches($query, $match_description]);

The argument $match_description (boolean) specifies if the query has
to be matched against the feature descrioption, in addition to the
ID and names. The comparison is case-insensitive. 

The argument $full specifies if the query must match the whole ID,
name or description (full=1), or if partial matcehs are allowed
(full=0).

=cut
sub matches {
    my ($self, $query, $match_description, $full) = @_;
    my $match = 0;
    my @to_match = ($self->get_attribute("id"),
		    $self->get_attribute("names"));
    if ($match_description) {
	push @to_match, $self->get_attribute("descr");
    }
    
    ## Regular expression
    my $expression;
    if ($full) {
	## Full match
	$expression = "^".$query."\$";
    } else {
	## Partial match
	$expression = $query;
    }

    foreach my $to_match (@to_match) {
	if ($to_match=~ /$expression/i) {
	    $match = 1;
	}
    }

#    &RSAT::message::Debug($self->get_attribute("id"), $query, $match, join (";", @to_match)) if ($main::verbose >= 0);
    return $match;
}


return 1;

__END__

