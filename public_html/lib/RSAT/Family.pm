###############################################################
#
# Class Family
#
package RSAT::Family;

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
%_attribute_cardinality = (id=>"SCALAR",
			   name=>"SCALAR",
			   organism=>"SCALAR",
			   size=>"SCALAR"
			  );

=pod

=head1 NAME

    RSAT::Family

=head1 DESCRIPTION

Class used to store a family (cluster) of genes.

=cut


return 1;

__END__

