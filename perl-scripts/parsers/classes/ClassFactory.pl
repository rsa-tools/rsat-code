################################################################
### the class factory manages class holders 
### a class holder contains a set of objects belonging to the same class
### (_object_type)
### it allows to automatically : 
### - access all objects of a class (_objects)
### - counters (_count) and ids, 
### - indexes (_id_index, _name_index),

package classes::ClassFactory;
use Data::Dumper;
{
   
  sub new_class {
    ### create an empty class factory
    ### usage:
    ###    $myclass = classes::ClassFactory->new_class(object_type=>$object_type,prefix=>$prefix);
    ### example:
    ###    $compounds = classes::ClassFactory->new_class(object_type=>"classes::Compound",prefix=>"comp_");

    my ($class,%args) = @_;
    my $object_type = "";

    ### check the object type 
    unless (defined($args{object_type})) {
	die "Error : cannot create a classholder without specified object type\n";
    }
    $object_type = $args{object_type};

#    die "HELLO\t'$object_type'";
    unless ($object_type->can("new")) {
      die "Error: $class->new_class $object_type is not a valid object type\n";
    }

    my $input_index = new classes::Index;
    my $output_index = new classes::Index;

    ### instantiate the new class holder
    my $class_holder = bless {
			      _object_type =>  $object_type,
			      _prefix =>  $args{prefix},
			      _count => 0,
			      _objects => [],
			      _id_index => {},
			      _name_index => {},
			      _input_index => $input_index,
			      _output_index => $output_index,
			      _attribute_count => {},
			      _attribute_cardinality => {},
			      _attribute_header => {},
			      _out_fields => []
			     }, $class;
    $class_holder->init();

    return $class_holder;
  }


  sub init() {
      my ($class_holder) = @_;
      $class_holder->set_attribute_header("xrefs", join("\t", "external_db", "external_id"));
  }

  ################################################################
  ### Create a new object and store its references in the class holder
  ### usage : $object = $class_holder->new_object(id=>$id,%other_args);
  sub new_object {
    my ($class_holder, %args) = @_;
    my $object_type = $class_holder->get_object_type();

    ### make sure the new object has an ID
    unless ($args{id}) {
	my $auto_id = sprintf "%s%6d", $class_holder->get_prefix(), $class_holder->get_count + 1;
	$auto_id =~ s/ /0/g;
	$args{id} = $auto_id;
    }
    if ($object = $object_type->new(%args)) {
	$class_holder->{_count}++;
	push @{$class_holder->{_objects}}, $object;
	${$class_holder->{_id_index}}{$args{id}} = $object;
    } else {
      die "Error: could not create object of type $object_type\n";
    }
    $class_holder->index_object_id($object);
    $class_holder->index_object_names($object);

    return $object;
  }

  ################################################################
  ### Add a previously created  object and store its references in the class holder
  ### usage : $object = $class_holder->add_object(id=>$id,%other_args);
  sub add_object {
    my ($class_holder, $object, %args) = @_;
    my $object_type = $class_holder->get_object_type();

    ### make sure the new object has an ID
    unless ($object->get_attribute("id")) {
      my $auto_id = sprintf "%s%6d", $class_holder->get_prefix(), $class_holder->get_count + 1;
      $auto_id =~ s/ /0/g;
      $args{id} = $auto_id;
    }
    $class_holder->{_count}++;
    push @{$class_holder->{_objects}}, $object;
    ${$class_holder->{_id_index}}{$args{id}} = $object;
    $class_holder->index_object_id($object);
    $class_holder->index_object_names($object);
    return $object;
  }

  ### accessors for the class holder attributes
  sub get_object_type {
    my ($class_holder) = @_;
    return $class_holder->{_object_type};
  }

  sub get_prefix {
    my ($class_holder) = @_;
    return $class_holder->{_prefix};
  }

  sub set_prefix {
    my ($class_holder, $new_prefix) = @_;
    $class_holder->{_prefix} = $new_prefix;
  }

  sub get_count {
    my ($class_holder) = @_;
    return $class_holder->{_count};
  }
  
  sub get_id_index {
    my ($class_holder) = @_;
    return %{$class_holder->{_id_index}};
  }
  
  sub get_name_index {
    my ($class_holder) = @_;
    return %{$class_holder->{_name_index}};
  }
  
  ################################################################
  ### returns all objects of a class holder
  sub get_objects {
    my ($class_holder) = @_;
    return @{$class_holder->{_objects}};
  }

  sub get_out_fields {
    ### returns all objects of a class holder
    my ($class_holder) = @_;
    return @{$class_holder->{_out_fields}};
  }

  sub set_out_fields {
    ### returns all objects of a class holder
    my ($class_holder, @out_fields) = @_;
    @{$class_holder->{_out_fields}} = @out_fields;
  }

  ## #############################################################
  ## Return an object given its ID or one of its names.
  ## Case insensitive : all keys are converted to uppsercase.
  sub get_object { 
    my ($class_holder, $key) = @_;
    my $id = undef;
    my $obj = undef;
    
    ### find key in id index
    if ($obj =  ${$class_holder->{_id_index}}{uc($key)}) {
      return $obj;
    }

    ### find key in name index
    if ($r_ids = ${$class_holder->{_name_index}}{uc($key)}) {
      my $id = $$r_ids[0];
       if ($obj = $class_holder->get_object($id)) {
	 return $obj;
       } else {
	 warn  ("Error: ",
		"class ", $class_holder->get_object_type(),
		"object indexed by name ($key) but not by ID ($id)\n"
		)
	     if ($main::verbose >= 1);
	 
	 return undef;
       }
    }

    ### fail to identify object
    return undef; 
  }

  #### set the multi-column header for expanded attributes. The $heade string should contain \t as field separator 
  sub set_attribute_header {
      my ($class_holder, $attr, $header) = @_;
      my $object_type = $class_holder->get_object_type();
      $object_type->_set_attribute_header($attr,$header);
  }

#  sub get_objects_by_name {
#    ### returns all objects whose name matches a given string
#    my ($class_holder, $query) = @_;
#    my @matching_objects = ();
#    $query = uc($query); ### case-insensitive
#    if (defined($$class_holder->{_name_index}->{$query})) {
#      return @$$class_holder->{_name_index}->{$query};
#    } else {
#      warn ("; WARGNING: could not identify object $query in class ",
#	    $class_holder->get_object_type(),
#	    "\n") 
#	  if ($main::verbose >= 1);
#      return;
#    }
#  }
  
  ### reset the list of objects inluded in the class
  ### this is mainly for having the appropriate pseudo_pointers after reloading a DB
  sub set_objects {
    my ($class_holder, @new_objects) = @_;
    @{$class_holder->{_objects}} = @new_objects;
    return @{$class_holder->{_objects}};
  }
  
  ################################################################
  ## index the names for a single object in the class_holder
  sub index_object_id {
    my ($class_holder, $object) = @_;
    unless ($id = $object->get_attribute("id")) {
	warn "Object $object has no ID\n";
	return;
    }
    ${$class_holder->{_id_index}}{uc($id)} = $object;
  }


  ################################################################
  ## index the names for a single object in the class_holder
  sub index_object_names {
      my ($class_holder, $object) = @_;
      unless (my $id = $object->get_attribute("id")) {
	  warn "Object $object has no ID\n";	
	  return;
      }
      foreach $name  ($id, $object->get_attribute("names")) {
	  my $index_key = uc($name); ### indexes ar case-insensitive
	  $index_key =~ s/^\s+//; ### suppress leading spaces
	  $index_key =~ s/\s+$//; ### suppress trailing spaces
	  push @{$class_holder->{_name_index}{$index_key}}, $id; 
      }
  }

  ################################################################
  ## delete the previous ID index and recalculates it from 
  ## from all objects contained in the class holder
  sub index_ids {
    my ($class_holder) = @_;
    %{$class_holder->{_id_index}} = ();
    foreach $object ($class_holder->get_objects()) {
      $class_holder->index_object_id($object);
#      next unless ($id = $object->get_attribute("id"));
#      ${$class_holder->{_id_index}}{$args{id}} = $object;
    }
    return %{$class_holder->{_id_index}};
  }

  
  ################################################################
  ## index object names for all the objects of the embedded class 
  ## Delete the previous ID index and recalculates it 
  ## from all objects contained in the class holder
  ## the index is now a has of lists
  ## this allows to stand for potential homonymy
  sub index_names {
      my ($class_holder) = @_;
      
      warn ("; indexing names for class ", 
	    $class_holder->get_object_type(), "\n") 
      if ($main::verbose >= 4);
      %{$class_holder->{_name_index}} = ();
      
      foreach $object ($class_holder->get_objects()) {
	  $object->unique_names(); #### make sure a name only appears once
	  $class_holder->index_object_names($object);
      }
      return %{$class_holder->{_name_index}};
  }
  
  ### generate an index with all input objects
  ### the input is a classes::Index, where 
  ### - keys are pseudo_pointers to the input object (ex classes::Gene::GeneID)
  ### - values are pseudo_pointers to the Activity
  sub index_inputs {
    my ($class_holder) = @_;
    my $class = $class_holder->get_object_type();
    die "Error: get_inputs method is not implemented for class $class\n"
	unless ($class->can("get_inputs"));
    my $index = $class_holder->{_input_index};
    foreach my $object ($class_holder->get_objects()) {
      my $pseudo_pointer = $object->get_pseudo_pointer();
      foreach my $input ($object->get_inputs()) {
        $index->add_value($input,$pseudo_pointer);
      }
    }
  }


  ### generate an index with all output objects
  ### the output is a classes::Index, where 
  ### - keys are pseudo_pointers to the output object (ex classes::Gene::GeneID)
  ### - values are pseudo_pointers to the Activity
  sub index_outputs {
    my ($class_holder) = @_;
    my $class = $class_holder->get_object_type();
    die "Error: get_outputs method is not implemented for class $class\n"
	unless ($class->can("get_outputs"));
    my $index = $class_holder->{_output_index};
    foreach my $object ($class_holder->get_objects()) {
      my $pseudo_pointer = $object->get_pseudo_pointer();
      foreach my $output ($object->get_outputs()) {
        $index->add_value($output,$pseudo_pointer);
      }
    }
  }



  ### export a class in different formats
  ### usage : 
  ### %tables = $class_holder->export("tab");
  ### returns a hash of strings, each hash entry representing a table
  ### - all unique (SCALAR) arguments are returned in a $tables{description}
  ### - multiple (ARRAY) arguments are returned in separate tables
  sub export {
      my ($class_holder,$format,$outfile) = @_;
      my $object_type = $class_holder->get_object_type();


      (my $short_class = $object_type) =~ s/.*:://g;
      $short_class = lc($short_class);

#      if ($main::dbms eq "mysql") {
#	  $comment_symbol = "#";
#      } else {
	  $comment_symbol = "--";
#      }

      ################################################################
      #### export the data in tab-delimited format
      if ($format eq "tab") {
	  my %tables = ();
	  my %attribute_cardinalities = $object_type->get_attribute_cardinalities();
#	  my %attribute_headers = $object_type->get_attribute_headers();
	  
	  #### out fields
	  my @out_fields = $class_holder->get_out_fields();
	  if ($#out_fields < 0) {
	      @out_fields = sort keys %attribute_cardinalities;
	  }

	  ### header of the main table, indicating the column content
	  my @scalar_fields = ();
	  foreach my $attribute (@out_fields) { 
	      my $cardinality = $attribute_cardinalities{$attribute} || "ARRAY";
	      my $header = $attribute_header{$attribute} || $attribute;
	      unless ($attribute eq "id") { ### ID  treated separately because comes as first column
		  if ($cardinality eq "SCALAR") {
		      push @scalar_fields, $header;
		  }
	      }
	  }
	  $tables{$short_class} = sprintf "$comment_symbol %-12s\t%s\n", "table", "main";
	  $tables{$short_class} .= sprintf "$comment_symbol %s %d\t%s\n", "field", 1, "id";
	  for $f (0..$#scalar_fields) {
	      $tables{$short_class} .= sprintf "$comment_symbol %s %d\t%s\n", "field", $f+2, $scalar_fields[$f];
	  }
	  
	  ### header line (all field names on the same line, separated by tabs; this is convenient to have head associate with column data with cut)
	  $tables{$short_class} .= sprintf "$comment_symbol header\n";
	  $tables{$short_class} .= join("\t",
					"$comment_symbol id", 
					@scalar_fields);
	  $tables{$short_class} .= "\n";
	  
	  ### dump object content 
	  foreach my $object ($class_holder->get_objects()) {
	      my $id = $object->get_attribute("id");
	      my @scalar_fields = ();
	      foreach my $attribute (@out_fields) { 
		  my $cardinality = $attribute_cardinalities{$attribute} || "ARRAY";
		  my $header = $object_type->get_attribute_header($attribute) || $attribute;
		  unless ($attribute eq "id") {
		      
		      #### SCALAR attibute
		      if ($cardinality eq "SCALAR") {
			  $value = $object->get_attribute($attribute);
			  push @scalar_fields, $value;
			  
		      #### ARRAY attribute
		      } elsif ($cardinality eq "ARRAY") {
			  my $table_name = lc($short_class."_".$attribute);
			  unless (defined($tables{$table_name})) {
			      $tables{$table_name} = "$comment_symbol id\t".$header."\n";
			  }
			  if (my @values = $object->get_attribute($attribute)) {


			      ### special treatment for names : add a
			      ### column with name qualifier. The
			      ### first name is labeled 'primary', the
			      ### next names 'alternate'
			      if ($attribute eq "names") {
				  $tables{$table_name} .= join("\t",
							       $id,
							       $values[0],
							       "primary");
				  $tables{$table_name} .= "\n";
				  foreach my $v (1..$#values) {
				      my $cell = $values[$v]; 
				      $tables{$table_name} .= join("\t", 
								   $id,
								   $cell,
								   "alternate");
				      $tables{$table_name} .= "\n"; 
				  }
			      } else {
				  foreach my $value (@values) {
				      my $cell = $value; 
				      $tables{$table_name} .= join("\t", 
								   $id,
								   $cell);
				      $tables{$table_name} .= "\n";
				  }
			      }
			  }
			  
			  #### EXPANDED attribute
		      } elsif ($cardinality eq "EXPANDED") {
			  my @expanded_fields = split "\t", $header;
			  my $table_name = lc($short_class."_".$attribute);
			  unless (defined($tables{$table_name})) {
			      $tables{$table_name} = "$comment_symbol id\t".$header."\n";
			  }
			  if (my @value_array = $object->get_attribute($attribute)) {
			      foreach my $array_pointer (@value_array) {
				  my @array_values = @{$array_pointer};
				  #### check the size of the array of values
				  if ($#array_values < $#expanded_fields) {
				      for $f ($#array_values+1..$#expanded_fields) {
					  $array_values[$f] = $main::null;
				      }
				  } elsif ($#array_values > $#expanded_fields) {
				      &main::ErrorMessage(join "\t", "Too many values for attribute", 
							  $attribute,
							  "fields", scalar(@expanded_fields),
							  "values", scalar(@array_values),
							  join ";", @expanded_fields, 
							  join ";", @array_values, 
							  "supplementary values are ignored",
							  "\n"
							  );
				      @array_values = @array_values[0..$#expanded_fields];
				  }
				  for my $f (0..$#array_values) {
				      if ($array_values[$f] eq "") {
					  $array_values[$f] = $main::null;
				      }
				  }
				  $tables{$table_name} .= join("\t", 
							       $id,
							       @array_values);
				  $tables{$table_name} .= "\n";
			      }
			  }
			  
		      }
		  }
	      }
	      
	      $tables{$short_class} .= $id;
	      if (scalar(@scalar_fields) > 0) {
		  $tables{$short_class} .= join("\t","", @scalar_fields);
	      }
	      $tables{$short_class} .= "\n";
	  }
	  return %tables;
	  
	  ### MLDBM export
#      } elsif (lc($format) eq "mldbm") {
#	  use MLDBM qw(DB_File Storable);
#	  die "Error: export called with empty file name\n"
#	      unless ($outfile);
#	  warn (";\n; ", &main::AlphaDate, " exporting class $object_type in format MLDBM to file $outfile\n")
#	      if ($main::verbose >= 1);
#	  tie %db, 'MLDBM', $outfile ||
#	      die "Cannnot tie to $outfile\n";
#	  ### previous version: serialize the whole class holder
#	  ### $db{$object_type} = $class_holder;
#	  
#	  ### serialize each object belonging to the class holder
#	  foreach my $object ($class_holder->get_objects()) {
#	      if (my $id = $object->get_attribute("id")) {
#		  $db{$id} = $object;
#	      }
#	  }
	  
	  untie %db;
	  warn ("; ", &main::AlphaDate, " export done\n")
	      if ($main::verbose >= 2);
	  return();
	  
	  
	  ### export on "obj" format
      } elsif ($format eq "obj") {
	  if ($outfile) {
	      open STDOUT, ">$outfile"  || die "Error : cannot write file $outfile\n"; 
	  }
	  warn (";\n; ", &main::AlphaDate, " exporting class ", $object_type, " to file '$outfile'\n") 
	      if ($main::verbose >= 4);
	  my @selected = $class_holder->get_out_fields();
#	  my @selected = @{$out_fields{$object_type}};
	  foreach my $object ($class_holder->get_objects()) {
	      $object->print_attributes($format, @selected);
	  }
	  warn ("; ", &main::AlphaDate, " class ", $object_type, " exported\n") 
	      if ($main::verbose >= 2);
	  close STDOUT if ($outfile);
	  return(1);
	  
	  
	  ### unsupported format
      } else {
	  die "Error: format $format is not supported for exporting a class\n";
      }
  }


  sub export_name_index {
      ### exports the name index hash as a separate file, 
      ### which can be used for consistency checkings
      ### usage: 
      ###   $class_holder->export_name_index($format,$outfile);
      ### supported formats: MLDBM
      
      my ($class_holder, $format, $outfile) = @_;
      my $object_type = $class_holder->get_object_type();
      
      if (lc($format) eq 'mldbm') {
	  ### MLDBM export
#	  use MLDBM qw(DB_File Storable);
	  
	  die "Error: export called with empty file name\n"
	      unless ($outfile);
	  warn (";\n; ", &main::AlphaDate, 
		" exporting name index for class $object_type in format MLDBM to file $outfile\n")
	      if ($main::verbose >= 4);
	  
	  tie %db, 'MLDBM', $outfile ||
	      die "Cannnot tie to $outfile\n";
	  
	  
	  
	  my %name_index = $class_holder->get_name_index();
	  
	  while (($key, $value) = each %name_index) {
	      $db{$key} = $value;
	  }
      
	  untie %db;
	  warn ("; ", &main::AlphaDate, " export done\n")
	      if ($main::verbose >= 2);
	  return;
	  
	  ### unsupported format
      } else {
	  die "Error: format '$format' is not supported for export_name_index\n";
      }
  }

  sub export_input_index {
    ### exports the input index hash as a separate file, 
    ### which can be used for cross-reference searching
    ### usage: 
    ###   $class_holder->export_input_index($format,$outfile);
    ### supported formats: MLDBM
    my ($class_holder, $format, $outfile) = @_;
    my $object_type = $class_holder->get_object_type();

    ### MLDBM export
    if (lc($format) eq 'mldbm') {
#      use MLDBM qw(DB_File Storable);
      die "Error: export called with empty file name\n"
	  unless ($outfile);
      warn (";\n; ", &main::AlphaDate, 
	    " exporting input index for class $object_type in format MLDBM to file $outfile\n")
	  if ($main::verbose >= 4);
      tie %db, 'MLDBM', $outfile ||
	  die "Cannnot tie to $outfile\n";
      
      my $input_index = $class_holder->{_input_index};
      foreach my $key ($input_index->get_keys()) {
	my @values = $input_index->get_values($key);
	$db{$key} = \@values;
      }

      untie %db;
      warn ("; ", &main::AlphaDate, " export done\n")
	  if ($main::verbose >= 2);
      return;
      
      ### unsupported format
    } else {
      die "Error: format '$format' is not supported for export_input_index\n";
    }
  }

  sub export_output_index {
    ### exports the output index hash as a separate file, 
    ### which can be used for cross-reference searching
    ### usage: 
    ###   $class_holder->export_output_index($format,$outfile);
    ### supported formats: MLDBM
    my ($class_holder, $format, $outfile) = @_;
    my $object_type = $class_holder->get_object_type();

    ### MLDBM export
    if (lc($format) eq 'mldbm') {
#      use MLDBM qw(DB_File Storable);
      die "Error: export called with empty file name\n"
	  unless ($outfile);
      warn (";\n; ", &main::AlphaDate, 
	    " exporting output index for class $object_type in format MLDBM to file $outfile\n")
	  if ($main::verbose >= 4);
      tie %db, 'MLDBM', $outfile ||
	  die "Cannnot tie to $outfile\n";
      
      my $output_index = $class_holder->{_output_index};
      foreach my $key ($output_index->get_keys()) {
	my @values = $output_index->get_values($key);
	$db{$key} = \@values;
      }

      untie %db;
      warn ("; ", &main::AlphaDate, " export done\n")
	  if ($main::verbose >= 2);
      return;
      
      ### unsupported format
    } else {
      die "Error: format '$format' is not supported for export_output_index\n";
    }
  }

  ### dumps the class in a tabular format
  ### uses the export("tab") function
  sub dump_tables {
      my ($class_holder,$file_suffix, $export_indexes) = @_;
      $file_suffix = "" unless (defined($file_suffix));
      my $object_type = $class_holder->get_object_type();
      (my $short_class = $object_type) =~ s/.*:://g;
      $short_class = lc($short_class);

      my $pwd = `pwd`;

      warn (";\n; ", &main::AlphaDate, " dumping class ", $class_holder->get_object_type(),
	    " to tables\n") if ($main::verbose >= 2);
      
      ### dump the class content in tables 
      ### according to relational normalization standards
      my %tables = $class_holder->export("tab");
      foreach my $table_name (keys %tables) {
	  my $file_name = $main::dir{output}."/".$table_name.$file_suffix.".tab";
	  open TABLE, ">$file_name" || die "Error: cannot write file $file_name\n";
	  my $dump_date = `date +%Y%m%d_%H%M%S`;
	  chomp($dump_date);
	  ### print version
	  printf TABLE "$comment_symbol %-12s\t%s\n", "dump date", $dump_date;
	  printf TABLE "$comment_symbol %-12s\t%s\n", "class", $class_holder->get_object_type(); 
	  printf TABLE "$comment_symbol %-12s\t%s\n", "table", $table_name; 
	  ### print content of the table
	  print TABLE $tables{$table_name};
	  close TABLE;
	  warn "; dumping table: $table_name in file $file_name\n" 
	      if ($main::verbose >= 4);
      }
      
      ### dump separate tables with indexes 
      if ($export_indexes) {

	  ### name index table
	  ### this table is a N to N correspondance
	  ### between names (all in uppercase) and object IDs
	  ### a priori, a class holder can contain hmonyms, i.e. several objects 
	  ### referred to by the same name
	  ### e.g.: a holder for BiochemicalEntities can contain 
	  ### a polypeptide and the gene coding for it that have the same name
	  my $file_name = lc($short_class.$file_suffix."__name_index.tab");
	  
	  warn "; dumping $object_type name index in file $file_name\n"
	      if ($main::verbose >= 2);
	  
	  open INDEX, ">$file_name" || 
	      die "Error: cannot write name index file '$file_name'\n";
	  
	  while (my ($name, $r_ids) = each %{$class_holder->{_name_index}}) {
	      foreach my $id (@$r_ids) {
		  print INDEX "$name\t$id\n";
	      }
	  }
	  close INDEX;
      
	  ### dump a separate table with the input_index index
	  if (($object_type->isa("classes::BiochemicalActivity")) &&
	      (defined( $class_holder->{_input_index}))){
	      my $file_name = lc($short_class.$file_suffix."__input_index.tab");
	      warn "; dumping $object_type input index in file $file_name\n"
		  if ($main::verbose >= 2);
	      open INDEX, ">$file_name" || 
		  die "Error: cannot write input index file '$file_name'\n";
	      my $input_index = $class_holder->{_input_index};
	      foreach my $key ($input_index->get_keys()) {
		  foreach my $value ($input_index->get_values($key)) {
		      print INDEX "$key\t$value\n";
		  }
	      }
	      close INDEX;
	  }
	  
	  
	  ### dump a separate table with the output index
	  if (($object_type->isa("classes::BiochemicalActivity")) &&
	      (defined( $class_holder->{_output_index}))){
	      my $file_name = lc($short_class.$file_suffix."__output_index.tab");
	      warn "; dumping $object_type output index in file $file_name\n"
		  if ($main::verbose >= 2);
	      open INDEX, ">$file_name" || 
		  die "Error: cannot write output index file '$file_name'\n";
	      my $output_index = $class_holder->{_output_index};
	      foreach my $key ($output_index->get_keys()) {
		  foreach my $value ($output_index->get_values($key)) {
		      print INDEX "$key\t$value\n";
		  }
	      }
	      close INDEX;
	  }
      }

      warn ("; ", &main::AlphaDate, " class ", $object_type, " dumped\n") 
	  if ($main::verbose >= 2);

      
      return;
  }

  ################################################################
  #### format a string into SQL comment
  sub print_sql_header {
      my ($sql_header) = @_;
      my $rep=64/length($comment_symbol);
      print SQL "\n", ${comment_symbol}x$rep, "\n";
      print SQL $comment_symbol, "\n";
      print SQL $comment_symbol, " ", $sql_header, "\n";
      print SQL $comment_symbol, "\n";
  }

  ################################################################
  ### Automatic generation of the SQL scripts for creating tables and loading the data
  ### usage
  ### $holder->generate_sql(schema=>"amaze",
  ###		            grant=>"granted_reader", 
  ###		            grant_write=>"granted_writer",
  ###		            dir=>"$dir{output}/sql_scripts", 
  ###		            prefix=>"prefix_");
  sub generate_sql {
      my ($class_holder,%args) = @_;
      foreach my $dbms (keys %main::supported_dbms) {
	  warn "; Exporting SQL scripts for $dbms\n" if ($main::verbose >= 4);
	  $class_holder->generate_sql_one_dbms(%args, dbms=>$dbms,
					       dir=>$main::dir{output}."/sql_scripts/".$dbms);
      }
  }


  ################################################################
  ## Calculate table path with or without schema as prefix
  sub get_table_path {
      my ($table_name, $schema, $full_path) = @_;
      my $table_path;
      if ($full_path) {
	  $table_path = $schema.".".$table_name;
      } else {
	  $table_path = $table_name;
      }
      return($table_path);
  }

  ################################################################
  #### Generate SQL scripts for one specific DBMS. This method is
  #### called iteratively by generate_sql with the different DBMS. 
  sub generate_sql_one_dbms {
      my ($class_holder,%args) = @_;
      my $object_type = $class_holder->get_object_type();
      (my $short_class = $object_type) =~ s/.*:://g;
      $short_class = lc($short_class);
      my $table_prefix = $short_class;
      if ($args{prefix}) {
	  $table_prefix = $args{prefix}.$table_prefix;
      }
      $table_prefix = lc($table_prefix); #### all in lowercase

      $max_col_length = 255; #beyond this length fields are of type long/text/longtext

     #### Database management system
      my $dbms = $args{dbms} || $main::default{dbms};
      my $host = $args{host} || $main::default{host};
      my $user = $args{user} || $main::default{user};
      my $password = $args{password} || $main::default{password};
      my $full_path = $args{full_path} || $main::default{full_path};
      my $schema = $args{schema} || $main::default{schema};
      my $sql_dir = $args{dir} || $main::dir{output}."/sql_scripts/$dbms";
      
      warn join ("\t", ";\n;", &main::AlphaDate(), 
		 " Generating SQL scripts for class ", $class_holder->get_object_type(),
		 "sql_dir: $sql_dir"),
		 "\n" if ($main::verbose >= 4);
      
      #### create SQL export directory if required
      unless (-d $sql_dir) {
	  warn "; Creating SQL dir $sql_dir" if ($main::verbose > 2);
	  system "mkdir -p $sql_dir";
	  unless (-e $sql_dir) {
	      die "Error: cannot create directory $sql_dir for exporting SQL scripts\n";
	  }
      }

      my $comment_symbol = "--";

      #### DBMS-specific options
      if ($dbms eq "mysql") {
	    $long_format = "LONGTEXT";
      } elsif ($dbms eq "postgresql") {
	    $long_format = "TEXT";
#	  $comment_symbol = "#";
      } else {
	    $long_format = "LONG";
      }

      #### grants
      my @granted_readers = ();
      my @granted_writers = ();
      if ($dbms eq "oracle") {
	  push @granted_readers, "amaze";
#	  push @granted_readers, $user."_reader";
#	  push @granted_writers, $user."_writer";
	  if ($args{grant}) {
	      push @granted_readers, $args{grant};
	  }
	  if ($args{grant_write}) {
	      push @granted_writers, $args{grant_write};
	  }
      }
      #### report DBMS options
      warn (";\tdbms\t", $dbms, "\n",
	    ";\thost\t", $host, "\n",
	    ";\tschema\t", $schema, "\n",
	    ";\tuser\t", $user, "\n",
	    ";\tpassword\t", $password, "\n",
	    ";\tgranted readers\t", join (",",@granted_readers), "\n",
	    ";\tgranted writers\t", join (",",@granted_writers), "\n",
	    ";\tpassword\t", $password, "\n",
	    ) if ($main::verbose >= 4);



      #### field sizes
      my $default_id_size = 32;
      my $default_field_size = 255;
      my $default_field_format = "VARCHAR($default_field_size)"; ### temporary : all fields are strings
      

      ################################################################
      #### get attribute cardinalities
      my %attribute_cardinalities = $object_type->get_attribute_cardinalities();
      
      #### out fields
      my @out_fields = $class_holder->get_out_fields();
      if ($#out_fields < 0) {
	  @out_fields = sort keys %attribute_cardinalities;
      }
      
      ### header of the main table, indicating the column content
      my @scalar_fields = ();
      my @array_fields = ();
      my @expanded_fields = ();
      foreach my $attribute (@out_fields) { 
	  my $cardinality = $attribute_cardinalities{$attribute} || "ARRAY";
	  my $header = $attribute_header{$attribute} || $attribute;
#	  while (my ($attribute,$cardinality) = each %attribute_cardinalities) {
	  unless ($attribute eq "id") { ### ID  treated separately because comes as first column
	      if ($cardinality eq "SCALAR") {
		  push @scalar_fields, $header;
	      } elsif ($cardinality eq "ARRAY") {
		  push @array_fields, $header;
	      } elsif ($cardinality eq "EXPANDED") {
		  push @expanded_fields, $header;
	      } 
	  }
      }


      ################################################################
      #### Table creation
      my $create_file = "${table_prefix}_table_create.sql";
      warn ";\ttable creation scripts to file $create_file\n" 
	  if ($main::verbose >= 4);
      open SQL, "> $sql_dir/$create_file" || die "Error: cannot write file $create_file\n";

      ################################################################
      #### Alter table
      my $alter_file = "${table_prefix}_table_alter.sql";
      warn ";\talter table scripts to file $alter_file\n" 
	  if ($main::verbose >= 4);
      open ALTER, "> $sql_dir/$alter_file" || die "Error: cannot write file $alter_file\n";

      print_sql_header ("Table creation scripts for class $table_prefix");

      #### schema
      if ($schema) {
	  print_sql_header ("Schema");
	  if ($dbms eq "oracle") {
	      print SQL "alter session set current_schema=", $schema, ";\n\n";
#	  }  elsif ($dbms eq "mysql") {
#	      print SQL "use ", $schema, ";\n\n";
	  }
      }

      #### name and path for the main table
      my $main_table_name = $table_name = $table_prefix;
      my $main_table_path = $table_path = &get_table_path($table_name, $schema, $full_path);
      print_sql_header ("Main table - $table_name");
      print SQL "CREATE TABLE $table_path", "\n\t(", "\n";
      my @field_defs = ();
      push @field_defs, sprintf "\t\t%-33s\t%s\t%s", "id", "VARCHAR($default_id_size)", "NOT NULL PRIMARY KEY";
      foreach my $field (@scalar_fields) {
	  $field = lc($field);
	  my $field_size = $main::special_field_size{$field} || $default_field_size;
	  my $field_format;
	  if ($field_size <= $max_col_length) {
	      $field_format  = "VARCHAR($field_size)";
	  } else {
	      $field_format  = $long_format;	      
	  }
#	  print STDERR join "\t", "HELLO SCALAR FIELD", $field, $default_field_size,  $main::special_field_size{$field}, $field_size, "\n";
	  push @field_defs, sprintf "\t\t%-33s\t%s", $field, $field_format;
      }
      print SQL join (",\n", @field_defs), "\n";
      if ($dbms eq "mysql") {
	  print SQL "\t) TYPE=INNODB", "\n", ";", "\n";
      } else {
	  print SQL "\t)", "\n", ";", "\n";
      }

      #### grant 
      foreach $g (@granted_readers) {
	  print SQL "GRANT select ON $table_path TO $g;\n";
      }
      foreach $g (@granted_writers) {
	  print SQL "GRANT select, insert, update ON $table_path TO $g;\n";
      }

      #### multivalue attributes
      foreach my $field (@array_fields) {
	  $field = lc($field);
	  my $table_name = $table_prefix."_".$field;
	  my $table_path = &get_table_path($table_name, $schema, $full_path);
	  print_sql_header ("Multivalue field - $table_name");
	  print SQL "CREATE TABLE $table_path", "\n\t(", "\n";
	  my @field_defs = ();
	  push @field_defs, sprintf "\t\t%-33s\t%s\t%s", "id", "VARCHAR($default_id_size)", "NOT NULL";
	  my $field_size = $main::special_field_size{$field} || $default_field_size;
	  my $field_format;
	  if ($field_size <= $max_col_length) {
	      $field_format  = "VARCHAR($field_size)";
	  } else {
	      $field_format  = $long_format;	      
	  }
#	  my $field_size = $main::special_field_size{$field} || $default_field_size;
#	  my $field_format = "VARCHAR($field_size)";
#	  print STDERR join "\t", "HELLO ARRAY FIELD", $field, $default_field_size,  $main::special_field_size{$field}, $field_size, $field_format, "\n";
	  push @field_defs, sprintf "\t\t%-33s\t%s", $field, $field_format;
	  #### quick and dirty fix for the fact that names are still not real expanded fields
	  if ($field eq "names") {
	      push @field_defs, sprintf "\t\t%-33s\t%s", "qualifier", "VARCHAR(20)";
	  }
	  print SQL join (",\n", @field_defs, 
#			  "INDEX(id)",
#			  "FOREIGN KEY (id) REFERENCES ${main_table_name}(id) ON DELETE CASCADE",
			  ), "\n";
	  if ($dbms eq "mysql") {
	      print SQL "\t) TYPE=INNODB", "\n", ";", "\n";
	  } else {
	      print SQL "\t)", "\n", ";", "\n";
	  }
	  print SQL "CREATE INDEX ".$table_name."_id_index ON ".$table_path." (id);\n";
	  print ALTER "ALTER TABLE ${table_path} ADD FOREIGN KEY (id) REFERENCES ${main_table_path}(id) ON DELETE CASCADE;\n";
	  foreach $g (@granted_readers) {
	      print SQL "GRANT select ON $table_path TO $g;\n";
	  }
	  foreach $g (@granted_writers) {
	      print SQL "GRANT select, insert, update ON $table_path TO $g;\n";
	  }
      }

      #### expanded attributes
      foreach my $field (@expanded_fields) {
	  $field = lc($field);
	  my $header = $object_type->get_attribute_header($field) || $field;
	  my @fields = split "\t", $header;
#	  shift @fields;
#	  warn "HELLO\tEXPANDED\t$field\t$cardinality\t$header\t", join (";", @fields), "\n";
	  my $table_name = $table_prefix."_".$field;
	  my $table_path = &get_table_path($table_name, $schema, $full_path);
	  print_sql_header ("Expanded attribute - $table_name");
	  print SQL "CREATE TABLE $table_path", "\n\t(", "\n";
	  my @field_defs = ();
	  push @field_defs, sprintf "\t\t%-33s\t%s\t%s", "id", "VARCHAR($default_id_size)", "NOT NULL";
	  my $field_size = $main::special_field_size{$field} || $default_field_size;
	  my $field_format;
	  if ($field_size <= $max_col_length) {
	      $field_format  = "VARCHAR($field_size)";
	  } else {
	      $field_format  = $long_format;	      
	  }
#	  my $field_size = $main::special_field_size{$field} || $default_field_size;
#	  my $field_format = "VARCHAR($field_size)";
#	  print STDERR join( "\t", "HELLO EXPANDED FIELD", 
#			     $field, 
#			     $default_field_size,  
#			     "HI", $main::special_field_size{features}, 
#			     $main::special_field_size{$field}, 
#			     $field_size,
#			     $field_format), "\n";
	  foreach my $field (@fields) {
	      $field = lc($field);
	      push @field_defs, sprintf "\t\t%-33s\t%s", $field, $field_format;
	  }
	  print SQL join (",\n", @field_defs, 
#			  "INDEX(id)",
#			  "FOREIGN KEY (id) REFERENCES ${main_table_name}(id) ON DELETE CASCADE",
			  ), "\n";
	  if ($dbms eq "mysql") {
	      print SQL "\t) TYPE=INNODB", "\n", ";", "\n";
	  } else {
	      print SQL "\t)", "\n", ";", "\n";
	  }
	  print SQL "CREATE INDEX ".$table_name."_id_index ON ".$table_path." (id);\n";
	  print ALTER "ALTER TABLE ${table_path} ADD FOREIGN KEY (id) REFERENCES ${main_table_path}(id) ON DELETE CASCADE;\n";
	  foreach $g (@granted_readers) {
	      print SQL "GRANT select ON $table_path TO $g;\n";
	  }
	  foreach $g (@granted_writers) {
	      print SQL "GRANT select, insert, update ON $table_path TO $g;\n";
	  }
      }

      #### close the table creation and alter script files
      if ($dbms eq "oracle") {
	  print_sql_header ("Quit");
	  print SQL "quit;";
      } 
      
      close SQL;
      close ALTER;

      ################################################################
      #### Table loading
      if ($dbms eq "oracle") {
	  #### load main table
	  my $table_name = $table_prefix;
	  my $table_path = &get_table_path($table_name, $schema, $full_path);
	  my $load_file = "${table_prefix}_table_load.ctl";
	  warn ";\ttable loading scripts to file $load_file\n" 
	      if ($main::verbose >= 4);
	  open SQL, "> $sql_dir/$load_file" || die "Error: cannot write file $load_file\n";
	  print_sql_header ("Table loading scripts for class $table_prefix");
	  print_sql_header ("Main table - $table_prefix");
	  print SQL "LOAD DATA", "\n";
	  print SQL "INFILE '", "../../${short_class}.tab", "'\n";
	  print SQL "APPEND INTO TABLE $table_path", "\n";
	  print SQL "REENABLE DISABLED_CONSTRAINTS", "\n";
	  print SQL "FIELDS TERMINATED BY X'09'", "\n";
#      print SQL "TRAILING NULLCOLS", "\n";
	  print SQL "(", "\n";
	  my @load_defs = ();
	  push @load_defs, sprintf "\t%-33s\t%s", "id", "CHAR($default_id_size)";
	  foreach my $field (@scalar_fields) {
	      $field = lc($field);
	      my $field_size = $main::special_field_size{$field} || $default_field_size;
	      my $load_field_format;
	      if ($field_size <= $max_col_length) {
		  $load_field_format  = "CHAR($field_size)";
	      } else {
		  $load_field_format  = "VARRAW";	      
	      }
#	      my  $field_size = $main::special_field_size{$field} || $default_field_size;
#	      my $load_field_format = "CHAR(${field_size})";
	      if ($main::null) {
		  push @load_defs, sprintf "\t%-33s\t%s\tNULLIF %s=\"%s\"", $field, $load_field_format, $field, $main::null;
	      } else {
		  push @load_defs, sprintf "\t%-33s\t%s", $field, $load_field_format;
	      }
	  }
	  print SQL join(",\n", @load_defs);
	  print SQL "\n";
#      print SQL "\tTERMINATED BY X'10'", "\n";
	  print SQL ")", "\n";
	  close SQL;
	  
	  #### load multivalue attributes
	  foreach my $field (@array_fields) {
	      $field = lc($field);
	      my $table_name = $table_prefix."_".$field;
  	      my $table_path = &get_table_path($table_name, $schema, $full_path);
	      my $load_file = "${table_prefix}_${field}_table_load.ctl";
	      warn ";\ttable loading scripts to file $load_file\n" 
		  if ($main::verbose >= 4);
	      open SQL, "> $sql_dir/$load_file" || die "Error: cannot write file $load_file\n";
	      print_sql_header ("Table loading scripts for class $table_prefix");
	      print_sql_header ("Multivalue attribute - $table_name");
	      print SQL "LOAD DATA", "\n";
	      print SQL "INFILE '", "../../${short_class}_${field}.tab", "'\n";
	      print SQL "APPEND INTO TABLE $table_path", "\n";
	      print SQL "REENABLE DISABLED_CONSTRAINTS", "\n";
	      print SQL "FIELDS TERMINATED BY X'09'", "\n";
#         print SQL "TRAILING NULLCOLS", "\n";
	      print SQL "(", "\n";
	      my @load_defs = ();
	      push @load_defs, sprintf "\t%-33s\t%s", "id", "CHAR($default_id_size)";
	      my $field_size = $main::special_field_size{$field} || $default_field_size;
	      my $load_field_format;
	      if ($field_size <= $max_col_length) {
		  $load_field_format  = "CHAR($field_size)";
	      } else {
		  $load_field_format  = "VARRAW";	      
	      }
#	      my $load_field_format = "CHAR(${field_size})";
#	      my $default_field_format = "CHAR(${default_field_size})";
	      my @fields = ();
	      if ($field eq "names") {
		  @fields = qw (names qualifier);
	      } else {
		  @fields = ($field);
	      }
	      foreach my $field (@fields) {
		  $field = lc($field);
		  if ($main::null) {
		      push @load_defs, sprintf "\t%-33s\t%s\tNULLIF %s=\"%s\"", $field, $load_field_format, $field, $main::null;
		  } else {
		      push @load_defs, sprintf "\t%-33s\t%s", $field, $load_field_format;
		  }
	      }
	      print SQL join(",\n", @load_defs);
	      print SQL "\n";
#         print SQL "\tTERMINATED BY X'10'", "\n";
	      print SQL ")", "\n";
	      close SQL;
	  }
	  
	  #### load expanded attributes
	  foreach $field (@expanded_fields) {
	      $field = lc($field);
	      my $table_name = $table_prefix."_".$field;
  	      my $table_path = &get_table_path($table_name, $schema, $full_path);
	      my $load_file = "${table_prefix}_${field}_table_load.ctl";
	      my $field_size = $main::special_field_size{$field} || $default_field_size;
	      my $load_field_format;
	      if ($field_size <= $max_col_length) {
		  $load_field_format  = "CHAR($field_size)";
	      } else {
		  $load_field_format  = "VARRAW";	      
	      }
#	      my $load_field_format = "CHAR(${field_size})";
#	      my $default_field_format = "CHAR(${default_field_size})";
	      my $header = $object_type->get_attribute_header($field) || $field;
	      my @fields = split "\t", $header;
#	      shift @fields;
	      warn ";\ttable loading scripts to file $load_file\n" 
		  if ($main::verbose >= 4);
	      open SQL, "> $sql_dir/$load_file" || die "Error: cannot write file $load_file\n";
	      print_sql_header ("Table loading scripts for class $table_prefix");
	      print_sql_header ("Multivalue attribute - $table_name");
	      print SQL "LOAD DATA", "\n";
	      print SQL "INFILE '", "../../${short_class}_${field}.tab", "'\n";
	      print SQL "APPEND INTO TABLE $table_path", "\n";
	      print SQL "REENABLE DISABLED_CONSTRAINTS", "\n";
	      print SQL "FIELDS TERMINATED BY X'09'", "\n";
#         print SQL "TRAILING NULLCOLS", "\n";
	      print SQL "(", "\n";
	      my @load_defs = ();
	      push @load_defs, sprintf "\t%-33s\t%s", "id", "CHAR($default_id_size)";
	      foreach my $field (@fields) {
		  $field = lc($field);
		  if ($main::null) {
		      push @load_defs, sprintf "\t%-33s\t%s\tNULLIF %s=\"%s\"", $field, $load_field_format, $field, $main::null;
		  } else {
		      push @load_defs, sprintf "\t%-33s\t%s", $field, $load_field_format;
		  }
	      }
	      print SQL join(",\n", @load_defs);
	      print SQL "\n";
#         print SQL "\tTERMINATED BY X'10'", "\n";
	      print SQL ")", "\n";
	      close SQL;
	  }

	  #### postgresql loader
      } elsif ($dbms eq "postgresql") {
	  my $table_name = $table_prefix;
	  my $table_path = &get_table_path($table_name, $schema, $full_path);
	  my $load_file = "${table_prefix}_table_load.ctl";
	  warn ";\ttable loading scripts to file $load_file\n" 
	      if ($main::verbose >= 4);
	  open SQL, "> $sql_dir/$load_file" || die "Error: cannot write file $load_file\n";
	  print_sql_header ("Table loading scripts for class $table_prefix");
	  print_sql_header ("Main table - $table_prefix");
#	  print SQL "copy $table_path from ../../${short_class}.tab with null as '", $main::null, "'\n";
	  print SQL "copy $table_path from stdin with null as '", $main::null, "'\n";
	  close SQL;
	  
	  foreach $field (@array_fields, @expanded_fields) {
	      $field = lc($field);
	      my $table_name = $table_prefix."_".$field;
	      my $table_path = &get_table_path($table_name, $schema, $full_path);
	      my $load_file = "${table_prefix}_${field}_table_load.ctl";
	      warn ";\ttable loading scripts to file $load_file\n" 
		  if ($main::verbose >= 4);
	      open SQL, "> $sql_dir/$load_file" || die "Error: cannot write file $load_file\n";
	      print_sql_header ("Table loading scripts for class $table_prefix");
	      print_sql_header ("Multivalue attribute - $table_name");
#	      print SQL "copy $table_path from ../../${short_class}_${field}.tab with null as '", $main::null, "'\n";
	      print SQL "copy $table_path from stdin with null as '", $main::null, "'\n";
	      close SQL;
	  }
	  
	  #### mysql loader
      } elsif ($dbms eq "mysql") {
	  my $table_name = $table_prefix;
	  my $table_path = &get_table_path($table_name, $schema, $full_path);
	  my $load_file = "${table_prefix}_table_load.ctl";
	  warn ";\ttable loading scripts to file $load_file\n" 
	      if ($main::verbose >= 4);
	  open SQL, "> $sql_dir/$load_file" || die "Error: cannot write file $load_file\n";
	  print_sql_header ("Table loading scripts for class $table_prefix");
#	  print SQL "use ", $schema, ";\n\n";
	  print_sql_header ("Main table - $table_prefix");
	  print SQL "LOAD DATA LOCAL INFILE \"../../${short_class}_filtered.tab\" INTO TABLE $table_path;\n";
	  close SQL;

	  
	  foreach $field (@array_fields, @expanded_fields) {
	      $field = lc($field);
	      my $table_name = $table_prefix."_".$field;
	      my $table_path = &get_table_path($table_name, $schema, $full_path);
	      my $load_file = "${table_prefix}_${field}_table_load.ctl";
	      warn ";\ttable loading scripts to file $load_file\n" 
		  if ($main::verbose >= 4);
	      open SQL, "> $sql_dir/$load_file" || die "Error: cannot write file $load_file\n";
	      print_sql_header ("Table loading scripts for class $table_prefix");
#	      print SQL "use ", $schema, ";\n\n";
	      print_sql_header ("Multivalue attribute - $table_name");
	      print SQL "LOAD DATA LOCAL INFILE \"../../${short_class}_${field}_filtered.tab\" INTO TABLE $table_path;\n";
	      close SQL;
	  }
	  
	  

      }
      
      ################################################################
      #### Table dropping
      my $drop_file = "${table_prefix}_table_drop.sql";
      warn ";\ttable dropping scripts to file $drop_file\n" 
	  if ($main::verbose >= 4);
      open SQL, "> $sql_dir/$drop_file" || die "Error: cannot write file $drop_file\n";
      
      print_sql_header ("Table droping scripts for class $table_prefix");

      #### schema
      if  ($schema) {
	  print_sql_header ("Schema");
	  if ($dbms eq "oracle") {
	      print SQL "alter session set current_schema=", $schema, ";\n\n";
#	  } elsif ($dbms eq "mysql") {
#	      print SQL "use ", $schema, ";\n\n";	      
	  }
      }

      #### main table
      my $table_name = $table_prefix;
      my $table_path = &get_table_path($table_name, $schema, $full_path);
      print_sql_header ("Main table - $table_prefix");
      print SQL "DROP TABLE ", $table_path, ";\n";

      #### multivalue attributes
      print_sql_header ("Multivalue attributes - $table_prefix");
      foreach $field (@array_fields) {
	  $field = lc($field);
	  my $table_name = $table_prefix."_".$field;
	  my $table_path = &get_table_path($table_name, $schema, $full_path);
	  print SQL "DROP TABLE ", $table_path, ";\n";
      }

      #### expanded attributes
      print_sql_header ("Expanded attributes - $table_prefix");
      foreach $field (@expanded_fields) {
	  $field = lc($field);
	  my $table_name = $table_prefix."_".$field;
	  my $table_path = &get_table_path($table_name, $schema, $full_path);
	  print SQL "DROP TABLE ", $table_path, ";\n";
      }

      #### close the table creation script file
      if ($dbms eq "oracle") {
	  print_sql_header ("Quit");
	  print SQL "quit;";
	  close SQL;
      }

      ################################################################
      #### generate the makefile
      my $makefile;
      if ($args{makefile}) {
	  $makefile = $args{makefile};
      } else {
	  $makefile = $table_prefix.".mk";
      }
      my $skip_lines= 7 + scalar(@scalar_fields);
      open MK, "> $sql_dir/${makefile}" || die "Error: cannot create the makefile\n";
      #### variables
      if ($dbms eq "oracle") {
	  print MK 'SQLLOGIN=',$user,'/',$password,'@',$host,"\n";
	  print MK 'SQLLOADER=sqlldr userid=${SQLLOGIN}', "\n";
	  print MK 'SQLPLUS=sqlplus ${SQLLOGIN}', "\n";
      } elsif ($dbms eq "mysql") {
#	  print MK 'MYSQL=mysql', "\n";
	  my $mysql = "mysql";
	  $mysql .= " --local-infile";
	  $mysql .= " -u $user" if $user;
	  $mysql .= " -D $schema" if $schema;
	  $mysql .= " -p$password" if $password;
	  $mysql .= " -h $host" if $host;
	  print MK "MYSQL=$mysql\n";
      }

      #### usage
      print MK "\nusage:\n";
      print MK "\t", '@perl -ne \'if (/^([a-z]\S+):/){ print "\t$$1\n";  }\' ', ${makefile}, "\n";

      #### all
      print MK "\nall: create alter load recompress\n";

      #### uncompression
      print MK "\nuncompress:\n";
      print MK "\tgunzip -f ../../${short_class}*.tab.gz\n";

      #### recompression
      print MK "\nrecompress:\n";
      print MK "\tgzip -f ../../${short_class}*.tab\n";

      #### table creation
      print MK "\ncreate:\n";
      print MK "\t", '@echo "--- Creating ',$table_prefix,'"', "\n";
      if ($dbms eq "oracle") {
	  print MK "\t", '${SQLPLUS} < ', $create_file, "\n";
      } elsif ($dbms eq "postgresql") {
	  print MK "\t", 'psql -f ', $create_file;
	  print MK " -d ", $schema if ($schema);
	  print MK "\n";
      } elsif ($dbms eq "mysql") {
	  print MK "\tcat $create_file | \${MYSQL}\n";
      }

      #### alter table 
      print MK "\nalter:\n";
      print MK "\t", '@echo "--- Alter ',$table_prefix,'"', "\n";
      if ($dbms eq "oracle") {
	  print MK "\t", '${SQLPLUS} < ', $alter_file, "\n";
      } elsif ($dbms eq "postgresql") {
	  print MK "\t", 'psql -f ', $alter_file;
	  print MK " -d ", $schema if ($schema);
	  print MK "\n";
      } elsif ($dbms eq "mysql") {
	  print MK "\tcat $alter_file | \${MYSQL}\n";
      }

      #### table loading
      my $load_file = "${table_prefix}_table_load.ctl";
      print MK "\nload:\n";
      print MK "\t", '@echo "--- Loading ',$short_class,'"', "\n";

      if ($dbms eq "oracle") {
	 print MK "\t", '${SQLLOADER} control=', $load_file, " skip=",$skip_lines,"\n";
	 foreach $field (@array_fields,@expanded_fields) {
	     $field = lc($field);
	     my $load_file = "${table_prefix}_${field}_table_load.ctl";
	     my $skip_lines=4;
	     print MK "\t", '${SQLLOADER} control=', ${load_file}, " skip=",$skip_lines,"\n";
         }
      } elsif ($dbms eq "postgresql") {
         my $table_file = "../../".$short_class.".tab";
	 print MK "\t", "grep -v '^$comment_symbol' ${table_file} | psql -f ", $load_file;
         print MK " -d ", $schema if ($schema);
         print MK "\n";
	 foreach $field (@array_fields,@expanded_fields) {
	     $field = lc($field);
             my $table_file = lc("../../".$short_class."_".$field.".tab");
	     my $load_file = "${table_prefix}_${field}_table_load.ctl";
	     my $skip_lines=4;

	     print MK "\t", "grep -v '^$comment_symbol' ${table_file} | psql -f ", $load_file;
             print MK " -d ", $schema if ($schema);
             print MK "\n";
         }
     } elsif ($dbms eq "mysql") {
	 print MK "\t", "\@if [ -f \"../../${short_class}.tab\" ] ; then ";
	 print MK "grep -v '^$comment_symbol' ../../${short_class}.tab > ../../${short_class}_filtered.tab"," ; ";
	 print MK "cat $load_file | \${MYSQL}"," ; ";
	 print MK "rm ../../${short_class}_filtered.tab"," ; ";
	 print MK "fi\n";
	 foreach $field (@array_fields,@expanded_fields) {
	     $field = lc($field);
	     my $load_file = "${table_prefix}_${field}_table_load.ctl";
	     print MK "\t", "\@if [ -f \"../../${short_class}_${field}.tab\" ] ; then ";
	     print MK "grep -v '^$comment_symbol' ../../${short_class}_${field}.tab > ../../${short_class}_${field}_filtered.tab"," ; ";
	     print MK "cat ${load_file} | \${MYSQL}"," ; ";
	     print MK "rm ../../${short_class}_${field}_filtered.tab"," ; ";
	     print MK "fi\n";
         }
      }

      #### table dropping
      print MK "\ndrop:\n";
      print MK "\t", '@echo "--- Dropping ',$table_prefix,'"', "\n";
      if ($dbms eq "oracle") {
          print MK "\t", '${SQLPLUS} < ', $drop_file, "\n";
      } elsif  ($dbms eq "postgresql") {
          print MK "\t", 'psql -f ', $drop_file;
          print MK " -d ", $schema if ($schema);
          print MK "\n";
      } elsif ($dbms eq "mysql") {
          print MK "\tcat $drop_file | \${MYSQL}\n";
      }
      close MK;


      return();
  }
  
}




return 1;
