### A generic class to treat all taxonomic classifications
### Examples of derived classes : 
### - SystematicGroup (organism classification)
### - ECnumber (enzyme classification)
### - FunctionalClass (gene function classification; to import from the MIPS for yeast)
package classes::ObjectSet;
{
  @ISA = qw ( classes::DatabaseObject );
  ### class attributes
  $_count = 0;
  $_prefix = "oset_";
  @_objects = ();
  %_name_index = ();
  %_id_index = ();
  %_attribute_count = ();
  %_attribute_cardinality = (
		      element_class=>"SCALAR",
		      id=>"SCALAR",
		      names=>"ARRAY",
		      supersets=>"ARRAY",
		      subsets=>"ARRAY",
		      elements=>"ARRAY",
		     );

  sub get_elements {
    ### return all elements comprized in a set by 
    ### iteratively collecting subset elements 
    ### (transitive closure)
    ### usage: @elements = $set->get_elements;
    my ($self) = @_;
    my @elements = ();
    push @elements, @{$self->{elements}};
    foreach my $sub_id (@{$self->{subsets}}) {
      if (my $subset = $class->get_object($sub_id)) {
	if (my @sub_elements = @{$subset->get_elements()}) {
	  push @elements, @sub_elements;
	}
      }
    }
    return @elements;
  }

  sub count_elements {
    ### count all elements comprized in a set and its subsets
    ### (transitive closure)
    ### usage : $count = $set->count_elements();
    my ($self) = @_;
    my $class = ref($self) || $self;
    my $element_count = 0;

    $element_count = $#{$self->{elements}} +1;
    foreach my $sub_id (@{$self->{subsets}}) {
      if (my $subset = $class->get_object($sub_id)) {
	if (my $subcount = $subset->count_elements()) {
	  $element_count += $subcount;
	}
      }
    }
    return $element_count;
  }

  sub expand {
    ### prints a tree in an expanded form, one line per set
    ### each set is followed by the list of its subsets,
    ### right-indented one step further
    ### parameters : 
    ### indent_level:   starting indent level (default = 1)
    ### indent_string:  string used for indentation (default = "  ")
    ### count_elements: return the number of elements in each set
    ### usage: $set->expand(indent_level=>$number,indent_string=>$str,count_elements=>1);
    my ($self, %args) = @_;
    my $indent_level = $args{indent_level} || 1;
    my $indent_string = $args{indent_string} || "  ";
    my $class = ref($self) || $self;
    
    print $indent_string x $indent_level; 
    
    print $self->get_id();
    if ($args{count_elements}) {
      printf "\t%6d", $self->count_elements;
    }
    print "\t", $self->get_name();
    print "\n";
    
    foreach my $sub_id (@{$self->{subsets}}) {
      if (my $subset = $class->get_object($sub_id)) {
	$subset->expand(%args,indent_level=>$indent_level+1);
      }
    }
    return 1;
  }

}


return 1;
