
### a persistent class consists basically in an index hash and a class
### hash the index hash provides, given a name, th list of IDs that
### correspond to that name. There might be several IDs associated to
### a single name (homonymy), or, at reverse, several names associated
### to a same ID (synonymy). 
package classes::PersistentClass; 
{ 
  ### open a new persistenc class and associate it to
  ### - a name index file
  ### - a class file (containig the objects)
  sub new {
    my ($class,%args) = @_;
    my %_class_hash = ();
    my %_name_index = ();
    die "Error: the class file should be specified to open a persistent class\n"
	unless (defined($args{class_file}));
    die "Error: the name index file should be specified to open a persistent class\n"
	unless (defined($args{name_index_file}));

    warn (";\n;\topening class \t", $args{object_type},"\n",
	  ";\tclass file     \t",$args{class_file},"\n",
	  ";\tname index file\t",$args{name_index_file},"\n")
	if ($main::verbose >= 1);


    ### open class object file
    die ("Error: cannot read class file ", $args{class_file}, "\n") 
	unless (-r $args{class_file});
    tie %_class_hash, 'MLDBM', $args{class_file} ||
	die "Error: cannot tie to ", $args{class_file},"\n";
    
    ### open name index file
    die ("Error: cannot read class file ", $args{name_index_file}, "\n") 
	unless (-r $args{name_index_file});
    tie %_name_index, 'MLDBM', $args{name_index_file} ||
	die "Error: cannot tie to ", $args{name_index_file},"\n";
      
    ### create the persistent class holder
    $class_holder = bless {
      object_type=> $args{object_type},
      class_file=>      $args{class_file},
      name_index_file=> $args{name_index_file},
      _class_hash=>     \%_class_hash,
      _name_index=>     \%_name_index
    }, $class;

    return $class_holder;
  }
  
  sub open_equation_index {
    my ($class_holder,$index_file) = @_;
    ### open equation index file
    die ("Error: cannot read class file ", $index_file, "\n") 
	unless (-r $index_file);
    tie %_equation_index, 'MLDBM', $index_file ||
	die "Error: cannot tie to ", $index_file,"\n";
    $class_holder->{equation_index_file} = $index_file;
    $class_holder->{_equation_index} = \%_equation_index;
  }

  sub open_input_index {
    my ($class_holder,$index_file) = @_;
    ### open input index file
    die ("Error: cannot read class file ", $index_file, "\n") 
	unless (-r $index_file);
    tie %_input_index, 'MLDBM', $index_file ||
	die "Error: cannot tie to ", $index_file,"\n";
    $class_holder->{input_index_file} = $index_file;
    $class_holder->{_input_index} = \%_input_index;
  }

  sub open_output_index {
    my ($class_holder,$index_file) = @_;
    ### open output index file
    die ("Error: cannot read class file ", $index_file, "\n") 
	unless (-r $index_file);
    tie %_output_index, 'MLDBM', $index_file ||
	die "Error: cannot tie to ", $index_file,"\n";
    $class_holder->{output_index_file} = $index_file;
    $class_holder->{_output_index} = \%_output_index;
  }

  ### return the type (class) of objects contained in the class holder
  sub get_object_type {
    my ($class_holder) = @_;
    return $class_holder->{object_type};
  }

  ### given an object name, return its ID
  sub get_id {
    my ($class_holder, $query) = @_;
    my $key = uc($query);
    
    if (${$class_holder->{_name_index}}{$key}) {
      my @IDs = @{${$class_holder->{_name_index}}{$key}};
      return $IDs[0];
    } else {
      return undef;
    }
  }

  ### return the list of IDs associated to a given name
  ### remind that indexes are hash of arrays, i.e.
  ### any index key can be associated to a list of objects
  ### this allows to deal with homonymy
  sub get_ids {
    my ($class_holder, $query) = @_;
    my $key = uc($query);
    
    if (${$class_holder->{_name_index}}{$key}) {
      my @IDs = @{${$class_holder->{_name_index}}{$key}};
      return @IDs;
    } else {
      return undef;
    }
  }

  ### given an object ID or name, return the object
  sub get_object {
    my ($class_holder, $query) = @_;
    if ($object = ${$class_holder->{_class_hash}}{$query}) {
      return $object;
    } elsif (($id = $class_holder->get_id($query)) &&
	     ($object = ${$class_holder->{_class_hash}}{$id})) {
      return $object;
    } else {
      return undef;
    }
  }

  sub get_objects_by_equation {
    ### works only for BiochemicalActivities !
    ### return all the activities in the class holder
    ### which match the equation (perfect matching is required)
    ### usage: 
    ###    @IDs = $class_holder->get_objects_by_equation($query_pseudo_pointer);
    my ($class_holder, $query) = @_;
    my $index = $class_holder->{_equation_index};
    my @pseudo_pointers = @{$index->{$query}};
    return @pseudo_pointers;
  }

  sub get_objects_by_input {
    ### works only for BiochemicalActivities !
    ### given a pseudo_pointer, find the referred object and
    ### return all the activities in the class holder
    ### for which the object is an input
    ### usage: 
    ###    @IDs = $class_holder->get_objects_by_input($query_pseudo_pointer);
    my ($class_holder, $query) = @_;
    my $index = $class_holder->{_input_index};
    my @pseudo_pointers = @{$index->{$query}};
    return @pseudo_pointers;
  }

  sub get_objects_by_output {
    ### works only for BiochemicalActivities !
    ### given a pseudo_pointer, find the referred object and
    ### return all the activities in the class holder
    ### for which the object is an output
    ### usage: 
    ###    @IDs = $class_holder->get_objects_by_output($query_pseudo_pointer);
    my ($class_holder, $query) = @_;
    my $index = $class_holder->{_output_index};
    my @pseudo_pointers = @{$index->{$query}};
    return @pseudo_pointers;
  }

  ### return the list of IDs for all the objects contained in the class holder
  sub get_all_ids {
    my ($class_holder) = @_;
    return keys (%{$class_holder->{_class_hash}});
  }
}



return 1;
