
######################################################################################
##################################### REFERENCES #####################################
######################################################################################



package PFBP::Reference;
{
  @ISA = qw ( PFBP::DatabaseObject );
  ### class attributes
  $_count = 0;
  $_prefix = "refr_";
  @_objects = ();
  %_name_index = ();
  %_id_index = ();
  %_attribute_count = ();
  %_attribute_cardinality = (
		      id=>"SCALAR"
		     );
}

package PFBP::Paper;
{
  @ISA = qw ( PFBP::DatabaseObject );
  ### class attributes
  $_count = 0;
  $_prefix = "papr_";
  @_objects = ();
  %_name_index = ();
  %_id_index = ();
  %_attribute_count = ();
  %_attribute_cardinality = (
		      id=>"SCALAR"
		     );
}

package PFBP::Book;
{
  @ISA = qw ( PFBP::DatabaseObject );
  ### class attributes
  $_count = 0;
  $_prefix = "book_";
  @_objects = ();
  %_name_index = ();
  %_id_index = ();
  %_attribute_count = ();
  %_attribute_cardinality = (
		      id=>"SCALAR"
		     );
}


return 1;
