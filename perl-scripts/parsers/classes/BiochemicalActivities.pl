###############################################################################################
################################## BIOCHEMICAL ACTIVITIES  ####################################
###############################################################################################


package classes::BiochemicalActivity;
{
  @ISA = qw ( classes::DatabaseObject );
  ### this is a super-class for all the biochemical activities
  ### it should be a method class only, instantiation necessitates 
  ### to specify the sub-class

  ### return the list of inputs
  ### the generic method takes the attribute called "input"
  sub get_inputs {
    my ($object) =  @_;
    my @pseudo_pionters = ();
    foreach my $input ($object->get_attribute("inputs")) {
      push @pseudo_pointers, $input->get_pseudo_pointer();
    }
    return @pseudo_pointers;
  }

  ### return the list of outputs
  ### the generic method takes the attribute called "output"
  sub get_outputs {
    my ($object) =  @_;
    my @pseudo_pionters = ();
    foreach my $output ($object->get_attribute("outputs")) {
      push @pseudo_pointers, $output->get_pseudo_pointer();
    }
    return @pseudo_pointers;
  }

}

package classes::Expression;
### activity connecting a gene to the polypeptide it codes for
{
  @ISA = qw ( classes::BiochemicalActivity );
  ### class attributes
  $_count = 0;
  $_prefix = "expr_";
  @_objects = ();
  %_name_index = ();
  %_id_index = ();
  %_attribute_count = ();
  %_attribute_cardinality = (id=>"SCALAR",
		      gene_id=>"SCALAR",
		      polypeptide_id=>"SCALAR",
		      source=>"SCALAR");

  sub get_inputs {
    my ($object) =  @_;
    my $gene_id = $object->get_attribute("gene_id");
    my $pseudo_pointer = "classes::Gene::".$gene_id;
    return $pseudo_pointer;
  }

  sub get_outputs {
    my ($object) =  @_;
    my $polypeptide_id = $object->get_attribute("polypeptide_id");
    my $pseudo_pointer = "classes::Polypeptide::".$polypeptide_id;
    return $pseudo_pointer;
  }

}

### reactant is just a class to store a hash attribute for Reaction
### it describes 
### - the compound ID, 
### - the stoichiometry 
### - the validity of the reactant as intermediate betweeen 2 successive reactions in a pathway
package classes::Reactant;
{
  @ISA = qw ( classes::DatabaseObject );
  ### class attributes
  $_count = 0;
  $_prefix = "rctt_";
  @_objects = ();
  %_name_index = ();
  %_id_index = ();
  %_attribute_count = ();
  %_attribute_cardinality = (id=>"SCALAR",
		      compound_id=>"SCALAR",
		      stoichio=>"SCALAR",
		      valid_interm=>"SCALAR"
		     );
}

package classes::HasReactants;
{
  @ISA = qw ( classes::BiochemicalActivity );
  ### return the list of reactant objects
  sub get_reactants {
    my ($activity, $reactant_type) =  @_;
    my @reactants = ();
    if ($reactant_type =~ /^substrate/i) {
      @reactant_types = ("substrates");
    } elsif ($reactant_type =~ /^product/i) {
      @reactant_types = ("products");
    } else {
      @reactant_types = ("substrates", "products");
    }
    
    foreach $reactant_type (@reactant_types) {
      foreach my $reactant_id ($activity->get_attribute($reactant_type)) {
	my $reactant = undef;
	if (
	    ($reactant = classes::Reactant->get_object($reactant_id)) ||
	    ((defined($main::reaction_reactants)) && 
	     ($reactant = $main::reaction_reactants->get_object($reactant_id))) ||
	    ((defined($main::assembly_reactants)) && 
	     ($reactant = $main::assembly_reactants->get_object($reactant_id))) ||
	    ((defined($main::reactants)) && 
	     ($reactant = $main::reactants->get_object($reactant_id))) ||
	    ((defined($main::database)) && 
	     ($reactant = $main::database->get_object("classes::Reactant::$reactant_id"))) 
	    ) {
#print STDERR "get_reactants\t$activity\t$reactant_type\t$reactant_id\t$reactant\n";
	  push @reactants, $reactant;
	} else {
	  warn ("Error: cannot retrieve reactant $reactant_id ",
		"for activity $activity ", $activity->get_attribute("id"),
		"\n")
	      if ($main::verbose >= 1);
	}
      }
    }
    return @reactants;
  }

  ### overloading of the method get_input, in order to treat the special case of reactants
  ### the pseudo_pointers returned refer to compounds rather than reactant objects
  sub get_inputs {
    my ($activity, %args) =  @_;
    my @pseudo_pionters = ();
    foreach my $reactant ($activity->get_reactants("substrates")) {
      $comp_id = $reactant->get_attribute("comp_id");
#print "hello get_inputs\t$activity\t$reactant\t$comp_id\n";
      
      if ($args{return} =~ "id") {
	### return input ids
	push @pseudo_pointers, $comp_id;
	
      } else  {
	### return input objects
	if (($compound = $main::database->get_object("classes::Compound::$comp_id")) ||
	    ($compound = classes::Compound->get_object($comp_id))) {
	  push @pseudo_pointers, $compound->get_pseudo_pointer();
	} else {
	  warn ("Error: cannot retrieve compound $comp_id ",
		"for reactant $reactant ", $reactant->get_attribute("id"),
		"\n")
	      if ($main::verbose >= 1);
	}
      }
    }
    return @pseudo_pointers;
  }

  sub calc_equation {
    ### calculates a activity equation, from the list of substrates and products
    ### each side of the equation is sorted alphabetically, which 
    ### allows fast equation matching by a simple string matching
    ### (idea from Christian Lemer)
    ### supported arguments:
    ### - reverse: calculate the reverse equation
    ### - names:   use compound primary names instead of compound IDs in the equation
    my ($activity, %args) =  @_;
    local @substrates = ();
    local @products = ();
    my $equation = "";
    foreach $side ("substrates", "products") {
      foreach my $reactant ($activity->get_reactants($side)) {
	my $term = "";
	my $stoichio = "";
	if (($stoichio = $reactant->get_attribute("stoichio")) &&
	    ($stoichio ne "1")) {
	  $term = "$stoichio ";
	}
	if (my $comp_id = $reactant->get_attribute("comp_id")) {
	  if ($args{names}) {
	    ### try to identify compound name
	    if (($compound = classes::Compound->get_object($comp_id)) ||
		((defined($main::compounds)) && 
		 ($compound = $main::compounds->get_object("classes::Compound::$comp_id"))) ||
		((defined($main::database)) && 
		 ($compound = $main::database->get_object("classes::Compound::$comp_id")))) {
	      
	      $comp_name = $compound->get_name();
	      $term .= $comp_name;
	    } else {
	      warn ("Error: cannot identify compound $comp_id\n")
		  if ($main::verbose > 2);
	      $term .= $comp_id
		}
	    ### use compound ID for calculating equation (default behaviour) 
	  } else {
	    $term .= $comp_id;
	  }
	  push @{$side}, $term;
	}
      }
    }
    my $substrates = join " + ", sort @substrates;
    my $products= join " + ", sort @products;
    if ($args{reverse}) { ### reverse equation only
      $equation = join (" -> ", $products, $substrates);
    } elsif ($args{direct}) { ### direct equation only
	$equation = join (" -> ", $substrates, $products);
    } else { ### direction inensitive equation  -> sort sides alphabetically as well
      $equation = join (" <=> ", 
			sort ($substrates,$products));
    }
    return $equation;
  }
}

package classes::Reaction;
{
  @ISA = qw ( classes::HasReactants );
  ### class attributes
  $_count = 0;
  $_prefix = "rctn_";
  @_objects = ();
  %_name_index = ();
  %_id_index = ();
  %_attribute_count = ();
  %_attribute_cardinality = (id=>"SCALAR",
		      equation=>"SCALAR",
		      description=>"SCALAR",
#		      substrates=>"ARRAY",
#		      products=>"ARRAY",
		      source=>"SCALAR");
}

package classes::Assembly;
{
  @ISA = qw ( classes::HasReactants );
  ### class attributes
  $_count = 0;
  $_prefix = "asmb_";
  @_objects = ();
  %_name_index = ();
  %_id_index = ();
  %_attribute_count = ();
  %_attribute_cardinality = (id=>"SCALAR",
		      equation=>"SCALAR",
		      substrates=>"ARRAY",
		      products=>"ARRAY",
		      source=>"SCALAR");
}

package classes::Translocation;
{
  @ISA = qw ( classes::HasReactants );
  ### class attributes
  $_count = 0;
  $_prefix = "trlc_";
  @_objects = ();
  %_name_index = ();
  %_id_index = ();
  %_attribute_count = ();
  %_attribute_cardinality = (id=>"SCALAR",
		      equation=>"SCALAR",
		      substrates=>"ARRAY",
		      products=>"ARRAY",
		      source=>"SCALAR");
}


package classes::IndirectInteraction;
{
  @ISA = qw ( classes::BiochemicalActivity );
  ### class attributes
  $_count = 0;
  $_prefix = "treg_";
  @_objects = ();
  %_name_index = ();
  %_id_index = ();
  %_attribute_count = ();
  %_attribute_cardinality = (id=>"SCALAR",
			     equation=>"SCALAR",
			     sign=>"SCALAR",
			     inputs=>"ARRAY",
			     output=>"SCALAR",
			     source=>"SCALAR");
}

package classes::Induction;
{
  @ISA = qw ( classes::BiochemicalActivity );
  ### class attributes
  $_count = 0;
  $_prefix = "treg_";
  @_objects = ();
  %_name_index = ();
  %_id_index = ();
  %_attribute_count = ();
  %_attribute_cardinality = (id=>"SCALAR",
		      equation=>"SCALAR",
		      sign=>"SCALAR",
		      inputs=>"ARRAY",
		      outputs=>"ARRAY",
		      source=>"SCALAR");
}


package classes::TranscriptRegul;
{
  @ISA = qw ( classes::BiochemicalActivity );
  ### class attributes
  $_count = 0;
  $_prefix = "treg_";
  @_objects = ();
  %_name_index = ();
  %_id_index = ();
  %_attribute_count = ();
  %_attribute_cardinality = (id=>"SCALAR",
		      equation=>"SCALAR",
		      sign=>"SCALAR",
		      inputs=>"ARRAY",
		      outputs=>"ARRAY",
		      source=>"SCALAR");
}


package classes::TransportFacilitation;
{
  @ISA = qw ( classes::BiochemicalActivity );
  ### class attributes
  $_count = 0;
  $_prefix = "trsp_";
  @_objects = ();
  %_name_index = ();
  %_id_index = ();
  %_attribute_count = ();
  %_attribute_cardinality = (id=>"SCALAR",
		      equation=>"SCALAR",
		      inputs=>"ARRAY",


		      outputs=>"ARRAY",
		      source=>"SCALAR");
}

package classes::Catalysis;
### activity connecting a proteinaceousentity (polypeptide or protein)
### to an EC number
{
  @ISA = qw ( classes::BiochemicalActivity );
  ### class attributes
  $_count = 0;
  $_prefix = "cata_";
  @_objects = ();
  %_name_index = ();
  %_id_index = ();
  %_attribute_count = ();
  %_attribute_cardinality = (id=>"SCALAR",
		      catalyst=>"SCALAR",
		      catalyzed=>"SCALAR",
		      source=>"SCALAR");
}



package classes::ControlOfControl;
{
  @ISA = qw ( classes::BiochemicalActivity );
  ### class attributes
  $_count = 0;
  $_prefix = "ctrl_";
  @_objects = ();
  %_name_index = ();
  %_id_index = ();
  %_attribute_count = ();
  %_attribute_cardinality = (id=>"SCALAR",
		      equation=>"SCALAR",
		      sign=>"SCALAR",
		      inputs=>"ARRAY",
		      outputs=>"ARRAY",
		      source=>"SCALAR");

}

return 1;
