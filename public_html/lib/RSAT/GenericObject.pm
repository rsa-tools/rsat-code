################################################################
### The parent class for all objects
### include generic accessor methods
package RSAT::GenericObject;
$_count = 0;
$_prefix = "";
# use Data::Dumper

################################################################
#### Instance methods
################################################################

################################################################
## Instantiate a new object
sub new {
    my ($class,%args) = @_;
#    warn join (" ", keys(%args));
    ### bless the new object
    my $self = bless {
    }, ref($class) || $class;
    $self->init(%args) if $self->can("init"); ### initialization
    return $self;
}


################################################################
### create object without initialization
### WARNING: the class attributes are thus not updated
### and some class methods (e.g. $class->get_object($id)) do not function
sub fast_new {
    my ($class,%args) = @_;
    ### bless the new object
    my $self = bless {}, $class;

    ### add the new object to the class list
    push @{$class."::_objects"}, $self;

    return $self;
}


################################################################
=pod

=item auto_id

Automatically assign a unique identifier to the object

=cut

sub auto_id {
  my ($self, $new_prefix) = @_;
  my $class = ref($self) || $self;
  my $id;
  if ($new_prefix) {
    $id = $new_prefix;
  } else {
    $id = $class->{_prefix} || $class;
  }

  ## Add a counter
  my $id_nb = sprintf "%6s", $class->{_count}++;
  $id_nb =~ s/ /0/g;
  $id .= $id_nb;

  return $id;
}

################################################################
## Initialize the object
sub init {
    my ($self, %args) = @_;
    my $class = ref($self);

    ### increment class counters
    $self->_incr_count;
    foreach $super (@{$class."::ISA"}) {
	$super->_incr_count() if $super->can("_incr_count");
    }
    
    #### initialize SCALAR variables to $null
    my %attribute_cardinalities = $class->get_attribute_cardinalities();
    my @scalar_fields = ();
    while (my ($attribute,$cardinality) = each %attribute_cardinalities) {
	unless ($attribute eq "id") { ### ID  treated separately because comes as first column
	    if ($cardinality eq "SCALAR") {
		$self->init_attribute($attribute, $main::null);
	    }
	}
    }
    
    ### assign attribute values
    while (($key, $value) = each %args) {
	$self->set_attribute($key,$value);
    }
    
    ### make sure the new object has an ID
    unless ($args{id}) {
# 	my $auto_id = sprintf "%s%6d", $self->get_prefix, $class->get_count;
# 	$auto_id =~ s/ /0/g;
	my $auto_id = $class->auto_id();
 	$self->set_attribute("id",$auto_id);
    }
    
    ### add the new object to the class list
    push @{$class."::_objects"}, $self;
    
}


################################################################
## return the value of the specified attribute the result is either a
## SCALAR or and ARRAY, depending on the attribute cardinality.
sub get_attribute {
    my ($self,$attr) = @_;
    my $class = ref($self);
    unless (defined($self->{$attr})) {
	warn("WARNING: object $self of class $class has no attribute named '$attr'\n") if ($main::verbose >= 6);
	return;
    }
    if (ref($self->{$attr}) eq "ARRAY") {
	return @{$self->{$attr}};
    } elsif (ref($self->{$attr}) eq "EXPANDED") { #### an array of arrays
	return @{$self->{$attr}};
    } elsif (ref($self->{$attr}) eq "HASH") {
	return %{$self->{$attr}};
    } else {
	return $self->{$attr};
    }
}


################################################################
## set the initial value for a SCALAR attribute do not increment the
## attribute counter
sub init_attribute {
    my ($self,$attr,$value) = @_;
    my $class = ref($self);
    $self->_set_attribute_cardinality($attr, "SCALAR");
    $self->{$attr} = $value;
}


################################################################
## set the value for a SCALAR attribute
sub set_attribute {
    my ($self,$attr,$value) = @_;
    my $class = ref($self);
    if (($self->get_attribute($attr)) &&
	($self->get_attribute($attr) ne $main::null)) {
	die ("\n\tError: object ", $self->get_attribute("id"),
	     "\n\talready has an attribute '$attr' with value ", $self->get_attribute($attr),
	     "\n\tand can not be reset to value $value", 
	     "\n\tto reset an attribute value, you should use the method force_attribute()", 
	     "\n");
    } else {
	$self->_set_attribute_cardinality($attr, "SCALAR");
	$self->_incr_attribute_count($attr);
	$self->{$attr} = $value;
    }
}

################################################################
## set the value for a SCALAR attribute without checking if the
## attribute already had a value (contrarily to set_attribute_value
sub force_attribute {
    my ($self,$attr,$value) = @_;
    my $class = ref($self);    
    $self->_set_attribute_cardinality($attr, "SCALAR");
    #$self->_incr_attribute_count($attr);
    $self->{$attr} = $value;
}

################################################################
## same as set_attribute but does not increment attribute counter
sub reset_attribute {
    my ($self,$attr,$value) = @_;
    $self->{$attr} = $value;
}

################################################################
## Set a whole set of key-value pairs for a HASH attribute
sub set_hash_attribute {
    my ($self, $attr, %hash) = @_;
    my $class = ref($self);
    $self->_set_attribute_cardinality($attr, "HASH");
    %{$self->{$attr}} = %hash;
#    die join "\t", "POUET", keys %hash;
}

################################################################
## add a key-value pair to a HASH attribute
sub add_hash_attribute {
    my ($self, $attr, $key, $value) = @_;
    my $class = ref($self);
    $self->_set_attribute_cardinality($attr, "HASH");
    $self->{$attr}->{$key} = $value;
}


################################################################
## concatenate the input string to the previous content of the variable
sub append_attribute {
    my ($self,$attr,$value) = @_;
    if ($self->get_attribute_cardinality($attr) eq "ARRAY") {
	### concatenate with the last array element
	${$self->{$attr}}[$#{$self->{$attr}}] .= $value;
    } else {  
	### concatenate with the scalar variable
	$self->_set_attribute_cardinality($attr, "SCALAR");
	$self->{$attr} .= $value;
    }
}

################################################################
## add a value to an ARRAY attribute
sub push_attribute {
    my ($self,$attr,@values) = @_;
    $self->_set_attribute_cardinality($attr, "ARRAY");
    $self->_incr_attribute_count($attr);
    push @{$self->{$attr}}, @values;
}

################################################################
## Set the value of an array attribute to the specified array
sub set_array_attribute {
    my ($self,$attr,@values) = @_;
    my $class = ref($self) || $self;
    $self->_set_attribute_cardinality($attr, "ARRAY");
    $self->_incr_attribute_count($attr);
    @{$self->{$attr}} = @values;
}

################################################################
## Set the value of an array attribute to the specified array
sub delete_array_attribute {
    my ($self,$attr) = @_;
    my $class = ref($self) || $self;
    $self->_set_attribute_cardinality($attr, "ARRAY");
    $self->_incr_attribute_count($attr);
    @{$self->{$attr}} = ();
}

################################################################
## add an entry to an EXPANDED attribute an EXPANDED is an array of
## entries, where #### each entry is itself an array of values
sub push_expanded_attribute {
    my ($self,$attr,@value_array) = @_;
    warn join "\t", "pushing expanded attribute", $self->get_attribute("id"), @value_array, "\n" if ($main::verbose >=3);
    $self->_set_attribute_cardinality($attr, "EXPANDED");
    $self->_incr_attribute_count($attr);
    push @{$self->{$attr}}, \@value_array;
}

################################################################
## select the appropriate method depending on the attribute cardinality
sub new_attribute_value {
    ### chooses between set_attribute or push_attribute 
    ### depending on the existence of a predefined type
    ### for that attribute
    my ($self,$attr,@value) = @_;
    my $class = ref($self) || $self;
    my $attr_cardinality = $class->get_attribute_cardinality($attr);
    if ($attr_cardinality eq "SCALAR") {
	$self->set_attribute($attr,$value[0]);
    } elsif ($attr_cardinality eq "EXPANDED") {
	$self->push_expanded_attribute($attr,@value);
    } else { ### ARRAY by default
	$self->push_attribute($attr,@value);
    }
}


################################################################
## set the value of an array attribute to an empty list
sub empty_array_attribute {
    my ($self,$attr) = @_;
    my $class = ref($self) || $self;
    $self->_set_attribute_cardinality($attr, "ARRAY");
    @{$self->{$attr}} = ();
}


################################################################
## print all attribute values 
## usage 
##    $self->print_attributes($format,@attrs);
## example
##    $self->print_attributes(tab,names,formula);
## supported formats:
##      tab    tab-dlimited text, one column per attribute
##      obj    objects embedded within curly braces
##      dump   use perl Data::Dumper library
##
## if @attrs is left empty, all attributes are printed
sub print_attributes {
    my ($self, $format, @attrs) = @_;
    my $class = ref($self) || $self;
    
    if ($#attrs < 0) {
	@attrs = $class->get_attribute_names();
    }
    
    ### tab-delimited text, one column per attribute
    if ($format =~ /^tab/i) {
	printf "%s", ref($self); 
	foreach $attr (@attrs) {
	    print "\t";
	    if (ref($self->{$attr}) eq "ARRAY") {
		print "[";
		foreach $i (0..$#{$self->{$attr}}) {
		    $value = ${$self->{$attr}}[$i];
		print "|" if ($i > 0);
		printf "%s", $value;
	    }
	    print "]";
	} else {
	    print "\t", $self->{$attr};
	}
    }
    print "\n";
    
    ### objects embedded within curly braces
} elsif ($format =~ /^obj/i) {
    printf "%s {\n", ref($self); 
    foreach $attr (@attrs) {
	my $attribute_cardinality = $class->get_attribute_cardinality($attr);
	
	### attribute cardinality ARRAY
	if ($attribute_cardinality eq "ARRAY") {
	    foreach $value ($self->get_attribute($attr)) {
		printf "    %-15s %s\n", $attr, $value;
	    }
	    
	    ### attribute cardinality HASH
	} elsif ($attribute_cardinality eq "HASH") {
	    my @pairs = ();
	    my %attr = $self->get_attribute($attr);
	    while (($key,$value) = each (%attr)) {
		push @pairs, "$key=>$value";
	    }
	    printf "    %-15s { %s }\n", $attr, join(" , ", @pairs);
	    
	    
	    ### attribute cardinality EXPANDED
	} elsif  ($attribute_cardinality eq "EPXANDED") {
	    foreach $array_pointer ($self->get_attribute($attr)) {
		print join "\t", @{$array_pointer};
		print "\n";
	    }

	    ### attribute cardinality SCALAR
	} else {
	    $value = $self->get_attribute($attr);
	    unless ($value eq "") {
		printf "    %-15s %s\n", $attr, $value;
	    }
	    
	}
    }
    print "}\n\n";
} else {
    die "Error: '$format' invalid format for print_attributes\n";
}
}

################################################################
## Get an object name, i.e. the first entry of the object attribute
## called "names" (it is a Vector) if there is not a single name for
## this object, return its ID
sub get_name {
    my ($self) = @_;
    if (@names = $self->get_attribute("names")) {
	return $names[0];
    } else {
	return $self->get_attribute("id");
    }
}


################################################################
## Check that names are unique preserve the original order of the
## names
sub unique_names {
    my ($self) = @_;
    my @new_names = ();
    my %names = ();
    if (@names = $self->get_attribute("names")) {
	foreach my $name (@names) {
	    unless ($names{lc($name)}) {
		push @new_names, $name;
		$names{lc($name)}++; 
	    }
	}
	@{$self->{names}} = @new_names;
    }
}

################################################################
## Return a pseudo_pointer to the object
sub get_pseudo_pointer {
    my ($object) = @_;
    my $class = ref($object);
    my $id = $object->get_attribute("id");
    my $pseudo_pointer = $class."::".$id;
    return $pseudo_pointer;
}

################################################################
## return an object's ID
sub get_id {
    my ($self) = @_;
    if ($id = $self->get_attribute("id")) {
	return $id;
    }
}


################################################################
################################################################
#### Class methods
################################################################
################################################################

################################################################
## accessors for class attributes
sub get_count { 
    my $class = ref($_[0]) || $_[0];
    return ${$class."::_count"}; 
}

################################################################
## Set the class counter to a specified value
sub set_count { 
    my ($class, $new_count) = @_;
    ${$class."::_count"} = $new_count; 
}

################################################################
## Increment the class counter
sub _incr_count {
    my $class = ref($_[0]) || $_[0];
    return ++${$class."::_count"}; 
}

################################################################
## Get rteh value of the class prefix
sub get_prefix { 
    my $class = ref($_[0]) || $_[0];
    return ${$class."::_prefix"}; 
}

################################################################
## Gets the list of objects of this class
sub get_objects { 
    my $class = ref($_[0]) || $_[0];
    return @{$class."::_objects"};
}

################################################################
## Set the list of objects for this class
sub set_objects { 
    my ($class,@new_objects) = @_;
    @{$class."::_objects"} = @new_objects;
}


################################################################
## Return an object given its ID or one of its names. The query is
## case insensitive: all keys are converted to uppercases
sub get_object { 
    my $class = ref($_[0]) || $_[0];
    my $key = uc($_[1]);
    my $id = "";
    my $obj = "";

    ### find key in id index
    if ($obj = ${$class."::_id_index"}{$key}) {
    #warn "get object\tclass : $class\tkey:$key\t$obj\n";
        return $obj;
    }

    ### find key in name index
    if (($id = ${$class."::_name_index"}{$key}) &&
        ($obj = $class->get_object($id))) {
    #warn "get object\tclass : $class\tkey:$key\t$id\t$obj\n";
        return $obj;
    }

    ### failed to identify the object
    return 0; 
}

################################################################
## Get the number of assignemtns per attribute
sub get_attribute_counts {
    my $class = ref($_[0]) || $_[0];
    return %{$class."::_attribute_count"};
}

################################################################
## get the index of object names
sub get_name_index {
    my $class = ref($_[0]) || $_[0];
    return %{$class."::_name_index"};
}

################################################################
## Get the index of object identifiers
sub get_id_index {
    my $class = ref($_[0]) || $_[0];
    return %{$class."::_id_index"};
}

################################################################
## Return the headers of all expanded attributes 
sub get_attribute_headers {
    my $class = ref($_[0]) || $_[0];
    return %{$class."::_attribute_header"};
}

################################################################
## Return the header for a specified expanded attribute
##
## usage: $header = $class->get_attribute_header($attr_name)
sub get_attribute_header { 
    my ($self, $attr) = @_;
    my $class = ref($self) || $self;
    return ${$class."::_attribute_header"}{$attr};
}


################################################################
## Get attribute cardinalities  
sub get_attribute_cardinalities {
    my $class = ref($_[0]) || $_[0];
    return %{$class."::_attribute_cardinality"};
}

################################################################
## get the cardinality of a specified attribute
sub get_attribute_cardinality { 
    ### usage: $cardinality = $class->get_attribute_cardinality($attr_name)
    my ($self, $attr) = @_;
    my $class = ref($self) || $self;
    return ${$class."::_attribute_cardinality"}{$attr};
}

#  sub get_attribute_type { 
#      ### usage: $type = $class->get_attribute_type($attr_name)
#      my ($self, $key) = @_;
#      my $class = ref($self) || $self;
#      return ${$class."::_attribute_type"}{$key};
#  }


################################################################
## Get the names of all attributes
sub get_attribute_names {
    ### usage : @names = $class->get_attribute_names();
    ### usage : @names = $object->get_attribute_names();
    my $class = ref($_[0]) || $_[0];
    return sort keys %{$class."::_attribute_count"};
}





################################################################
## Index all objects by name. This method uses the object ID as value
## index the "names" attribute in the %name_index hash the "id"
## attribute is indexed as well, so each object can be called by
## either its id or one of its names name_index allows to retrieve
## the object id on basis of any of the indexed names id_index
## allows then to retrieve the object reference. 
## Indexing is case-insensitive, i.e. all keys are converted to
## uppercases in the hash table
##
## usage: $class->index_names(); 
sub index_object_names {
    my ($object) = @_;
    my $class = ref($object);
    my @keys = ();
    $index_value = $object->get_attribute("id");
    push @keys, $index_value;
    push @keys, $object->get_attribute("names");
    foreach $key (@keys) {
	${$class."::_name_index"}{uc($key)} = $index_value 
	    if (($key) && ($index_value));
}
}

################################################################
## Index all objects on basis of their identifier
sub index_ids {
    my ($class) = @_;
    foreach $object ($class->get_objects()) {
	$object->index_id();
    }
}


################################################################
## Retrieve an object on basis of its identifier.
##   usage: $object->index_id($name); ### uses the object reference as value
##   usage: $object->index_id($name, $value); ### use specified value
## WARNING: indexing is case-insensitive, 
## i.e. all keys are converted to uppercases in the hash table
sub index_id {
    my ($self, $key, $value) = @_;
    my $class = ref($self) || $self;
    $key = $self->get_attribute("id") unless ($key);
    $value = $self unless ($value);
    ${$class."::_id_index"}{uc($key)} = $value if (($key) && ($value));
}


################################################################
## Index the "names" attribute in the %name_index hash the "id"
## attribute is indexed as well, so each object can be called by
## either its id or one of its names name_index allows to retrieve the
## object id on basis of any of the indexed names.
##
## Uses the object ID as value.
##
## $class->id_index() allows then to retrieve the object reference 
##
## Indexing is case-insensitive, i.e. all keys are converted
## to uppercases in the hash table
##
## usage: $class->index_names();
## 
sub index_names {
    my $class = ref($_[0]) || $_[0];
    foreach $object ($class->get_objects()) {
	$object->unique_names(); #### make sure a name only appears once
	my @keys = ();
	$index_value = $object->get_attribute("id");
	push @keys, $index_value;
	push @keys, $object->get_attribute("names");
	foreach $key (@keys) {
	    ${$class."::_name_index"}{uc($key)} = $index_value 
		if (($key) && ($index_value));
    }
}
}

################################################################
## Increment the counter of assignments for a given attribute
sub _incr_attribute_count { 
    my $class = ref($_[0]) || $_[0];
    my $key = $_[1];
    ${$class."::_attribute_count"}{$key}++ if ($key);
}

################################################################
## Set the cardinality of an attribute
sub _set_attribute_cardinality { 
    my ($self, $key, $cardinality) = @_;
    my $class = ref($self) || $self;
    ${$class."::_attribute_cardinality"}{$key} = $cardinality;
}

################################################################
## Set the header for an expanded attribute
sub _set_attribute_header { 
    my ($self, $attr, $header) = @_;
    my $class = ref($self) || $self;
    ${$class."::_attribute_header"}{$attr} = $header;
}

################################################################
#  sub _set_attribute_type { 
#    my ($self, $key, $type) = @_;
#    my $class = ref($self) || $self;
#    ${$class."::_attribute_type"}{$key} = $type;
#  }


return 1;
