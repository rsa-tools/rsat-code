### the package Database allows to open the persistent classes in read mode
### and to retrieve objects
#### This requires a prior installation of the MLDBM  module
#### This routine is probably obsolete, I cannot test it anymore
#### because I don't have this module on my machine anymore.
package classes::Database;
{
  my $data_dir = $parsed_data || $ENV{HOME}."parsed_data";
  my $mldbm_dir = "$data_dir/mldbm/";
  
  my %name_index_file = ();
  $name_index_file{"classes::ECSet"} = "$data_dir/ec/ECSet__name_index.mldbm";
  $name_index_file{"classes::Compound"} = "$data_dir/compounds/Compound__name_index.mldbm";
  $name_index_file{"classes::Reaction"} = "$data_dir/reactions/Reaction__name_index.mldbm";
  $name_index_file{"classes::Reactant"} = "$data_dir/reactions/Reactant__name_index.mldbm";

  $name_index_file{"classes::Gene"} = "$data_dir/genes/Gene__name_index.mldbm";
  $name_index_file{"classes::Expression"} = "$data_dir/genes/Expression__name_index.mldbm";

  $name_index_file{"classes::Polypeptide"} = "$data_dir/polypeptides/Polypeptide__name_index.mldbm";
#  $name_index_file{"classes::Catalysis"} = "$data_dir/polypeptides/Catalysis__name_index.mldbm";

  my %input_index_file = ();
  $input_index_file{"classes::Expression"} = "$data_dir/genes/Expression__input_index.mldbm";
  
  my %output_index_file = ();
  $output_index_file{"classes::Expression"} = "$data_dir/genes/Expression__output_index.mldbm";
  
  my %equaton_index_file = ();
  $equation_index_file{"classes::Reaction"} = "$data_dir/reactions/Reaction__equation_index.mldbm";
  
  my %class_file = ();
  $class_file{"classes::ECSet"} = "$data_dir/ec/ECSet.mldbm";
  $class_file{"classes::Compound"} = "$data_dir/compounds/Compound.mldbm";
  $class_file{"classes::Reaction"} = "$data_dir/reactions/Reaction.mldbm";
  $class_file{"classes::Reactant"} = "$data_dir/reactions/Reactant.mldbm";

  $class_file{"classes::Gene"} = "$data_dir/genes/Gene.mldbm";
  $class_file{"classes::Expression"} = "$data_dir/genes/Expression.mldbm";

  $class_file{"classes::Polypeptide"} = "$data_dir/polypeptides/Polypeptide.mldbm";
#  $class_file{"classes::Catalysis"} = "$data_dir/polypeptides/Catalysis.mldbm";
  
  
  sub new {
    ### open an empty database
    my ($class) = @_;
    
    $database = bless { 
    }, $class;
    return $database;
  }

  sub open_all_classes {
    my ($database) = @_;
    foreach my $class (keys %class_file) {
      $database->open_class($class);
    }
  }

  sub open_class {
    my ($database, $class) = @_;
    $class_holder = classes::PersistentClass->new(class_file=>       $class_file{$class},
					       name_index_file=>  $name_index_file{$class},
					       object_type=>  $class); 
    ### open various indexes (by input, output, equation)
    $class_holder->open_input_index($input_index_file{$class})
	if (defined($input_index_file{$class}));
    $class_holder->open_output_index($output_index_file{$class})
	if (defined($output_index_file{$class}));
    $class_holder->open_equation_index($equation_index_file{$class})
	if (defined($equation_index_file{$class}));
    $database->{$class} = $class_holder;
  }

  sub get_class_holder {
    my ($database, $class) = @_;
    return $database->{$class};
  }

  sub get_id {
    ### return the first ID associated to the ley in the name index hash
    ### remind that indexes are hash of arrays, i.e.
    ### any index key can be associated to a list of objects
    ### this allows to deal with homonymy
    my ($database, $query) = @_;
    my ($class, $key) = classes::DatabaseManagement->split_pseudo_pointer($query);
    my $class_holder = $database->get_class_holder($class);
    return $class_holder->get_id($key);
  }

  sub get_ids {
    ### return the complete list of IDs associated to the ley in the name index hash
    ### remind that indexes are hash of arrays, i.e.
    ### any index key can be associated to a list of objects
    ### this allows to deal with homonymy
    my ($database, $query) = @_;
    my ($class, $key) = classes::DatabaseManagement->split_pseudo_pointer($query);
    my $class_holder = $database->get_class_holder($class);
    return $class_holder->get_ids($key);
  }

  sub get_object {
    my ($database, $query) = @_;
    my ($class, $key) = classes::DatabaseManagement->split_pseudo_pointer($query);
    if (defined($database->{$class})) {
      return $database->{$class}->get_object($key);
    } else {
      warn "Error: get_object called with query '$query'\tdatabase does not contain a class called $class\n";
      return undef;
    }
  }

  sub get_objects_by_equation {
    my ($database, $class, $query) = @_;
    return $database->{$class}->get_objects_by_equation($query);
  }

  sub get_objects_by_input {
    my ($database, $class, $query) = @_;
    return $database->{$class}->get_objects_by_input($query);
  }

  sub get_objects_by_output {
    my ($database, $class, $query) = @_;
    return $database->{$class}->get_objects_by_output($query);
  }


}
return 1;
