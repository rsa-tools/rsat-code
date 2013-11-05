###############################################################
#
# Class Analysis
#
package RSAT::Analysis;

use RSAT::GenericObject;
use RSAT::error;
@ISA = qw( RSAT::GenericObject );

### class attributes
$_count = 0;
$_prefix = "ctg_";
@_objects = ();
%_name_index = ();
%_id_index = ();
%_attribute_count = ();
%_attribute_cardinality = (id=>"SCALAR");

=pod

=head1 NAME

    RSAT::Analysis

=head1 DESCRIPTION

This class is used by multiple-family-analysis to store the parameters
of one analysis.

=cut

return 1;


__END__

