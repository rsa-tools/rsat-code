

######################################################################################
######################### PATHWAYS AND GENERIC OBJECTS ###############################
######################################################################################

package PFBP::ProcessNode;
### A class to treat process nodes
{
    @ISA = qw ( PFBP::ObjectSet );
    ### class attributes
    $_count = 0;
    $_prefix = "pthw_";
    @_objects = ();
    %_name_index = ();
    %_id_index = ();
    %_attribute_count = ();
    %_attribute_cardinality = (id=>"SCALAR");
}

package PFBP::Process;
### A class to treat processes
{
    @ISA = qw ( PFBP::ProcessNode );
    ### class attributes
    $_count = 0;
    $_prefix = "proc_";
    @_objects = ();
    %_name_index = ();
    %_id_index = ();
    %_attribute_count = ();
    %_attribute_cardinality = (id=>"SCALAR",
			       source=>"SCALAR",
			       description=>"SCALAR",
			       names=>"ARRAY",
			       nodes=>"ARRAY",
			       arcs=>"EXPANDED"
			       );
}


package PFBP::ProcessLeaf;
### A class to treat process leaves
{
    @ISA = qw ( PFBP::ProcessNode );
    ### class attributes
    $_count = 0;
    $_prefix = "leaf_";
    @_objects = ();
    %_name_index = ();
    %_id_index = ();
    %_attribute_count = ();
    %_attribute_cardinality = (id=>"SCALAR",
			       interaction=>"SCALAR",
			       reverse=>"SCALAR");
}


package PFBP::Pathway;
#### OBSOLETE: maintained only for backward compatibility
### A class to treat Pathways
{
  @ISA = qw ( PFBP::ObjectSet );
  ### class attributes
  $_count = 0;
  $_prefix = "pthw_";
  @_objects = ();
  %_name_index = ();
  %_id_index = ();
  %_attribute_count = ();
  %_attribute_cardinality = (id=>"SCALAR",
		      description=>"SCALAR",
		      source=>"SCALAR",
		      names=>"ARRAY",
		      #element_class=>"SCALAR",
		      elements=>"ARRAY",
		      subsets=>"ARRAY",
		      supersets=>"ARRAY"
		     );
}

package PFBP::ECSet;
### A class to treat EC numenclature
{
  @ISA = qw ( PFBP::ObjectSet );
  ### class attributes
  $_count = 0;
  $_prefix = "ecst_";
  @_objects = ();
  %_name_index = ();
  %_id_index = ();
  %_attribute_count = ();
  %_attribute_cardinality = (id=>"SCALAR",
		      names=>"ARRAY",
		      #element_class=>"SCALAR",
		      elements=>"ARRAY",
		      ec=>"SCALAR",
		      description=>"SCALAR",
		      parent=>"SCALAR",
#		      subsets=>"ARRAY",
#		      supersets=>"ARRAY"
		     );
}


return 1;
