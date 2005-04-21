#!/usr/bin/perl
############################################################
#
# $Id: loading_util.pl,v 1.2 2005/04/21 22:55:08 jvanheld Exp $
#
# Time-stamp: <2003-07-10 11:52:50 jvanheld>
#
############################################################
### util.pl
### utilities for the AMAZE project


#require "config.pl";
require "lib/load_classes.pl";
require "lib/util.pl";
require Data::Dumper;

### initialization
$index_file{"ECSet"} = "$parsed_data/ec/ECSet__name_index.mldbm";
$index_file{"Compound"} = "$parsed_data/compounds/Compound__name_index.mldbm";
$index_file{"Reaction"} = "$parsed_data/reactions/Reaction__name_index.mldbm";
$index_file{"Reactant"} = "$parsed_data/reactions/Reactant__name_index.mldbm";

$class_file{"ECSet"} = "$parsed_data/ec/ECSet.mldbm";
$class_file{"Compound"} = "$parsed_data/compounds/Compound.mldbm";
$class_file{"Reaction"} = "$parsed_data/reactions/Reaction.mldbm";
$class_file{"Reactant"} = "$parsed_data/reactions/Reactant.mldbm";

sub OpenClass {
  ### usage:
  ###     &OpenClass($class);
  ### example:
  ###     &OpenClass("Compound");
  ### Any object of the selected class can then be accessed via its ID  
  ### ex: my $compound = $Compound_class{$comp_id};


  my @classes = @_;
  my $filename = undef;

  foreach $class (@classes) {
    if (defined ($class_file{$class})) {
      $filename = $class_file{$class};
    } else {
      die "Error: there is no class file for class $class\n";
    }
    die "Error: file $filename does not exist\n"
	unless (-e $filename);
    
    tie %{"main::${class}_class"}, 'MLDBM', $filename ||
	die "Error: cannot tie to $filename\n";
  }
  return;
}

sub CloseClass {
  my @classes = @_;
  foreach $class (@classes) {
    untie %{"main::${class}_class"}    ||
	die "Error: cannot untie $class class\n";
  }
  return;
}

sub OpenIndex {
  ### usage:
  ###     &OpenIndex($class);
  ### example:
  ###     &OpenIndex("Compound");
  ### Any object of the selected cass can then be identified with the hash
  ### %main::${class}_index
  ### ex: $name = "h2o"; print ${$main::Compound_index{uc($name)}}[0];


  my @classes = @_;
  my $filename = undef;

  foreach $class (@classes) {
    if (defined ($index_file{$class})) {
      $filename = $index_file{$class};
    } else {
      die "Error: there is no index file for class $class\n";
    }
    die "Error: file $filename does not exist\n"
	unless (-e $filename);
    
    tie %{"main::${class}_index"}, 'MLDBM', $filename ||
	die "Error: cannot tie to $filename\n";
  }
  return;
}

sub CloseIndex {
  my @classes = @_;
  foreach $class (@classes) {
    untie %{"main::${class}_index"}    ||
	die "Error: cannot untie $class index\n";
  }
  return;
}

sub ImportHash {
  my ($filename) = @_;
  my %hash = ();
  my $to_compress = 0;
  
  unless (-e  $filename) {
    if (-e "$filename.gz")  { 
      ### temporarily uncompress file
      warn ("; uncompressing file $filename.gz\n")
	  if ($verbose >= 1);
      system "gunzip $filename.gz";
      $to_compress = 1;
    } else {
      die "Error: file $filename does not exist\n";
    }
  }
  
  warn (";\n; ", &AlphaDate, " importing hash from $filename\n")
      if ($verbose >= 1);
  
  my %db = ();
  tie %db, 'MLDBM', $filename
      || die "Error: cannot tie to $filename\n";
  
  %hash = %db;
  

  untie %db;
  
  warn ("; ", &AlphaDate, " imported\n")
      if ($verbose >= 1);

  if ($to_compress) {
    warn ("; compressing file $filename\n")
	if ($verbose >= 1);
    system "gzip $filename";
  }

  return %hash;

}


sub ImportClass {
  ### import a class holder and its content from a MLDBM file
  ### the class holder must have been previously exported with 
  ### the method $class_holder->export('MDBLM',$file)
  ###
  ### usage: 
  ###     $class_holder = &ImportClass($file, $object_type);
  ###     $class_holder = &ImportClass($file);
  ###
  ### example: 
  ###     $commpounds = &ImportClass("$mldbm_dir/Compound.mldbm","classes::Compound");
  ###
  ### If object type is specified, checks that there is a classholder
  ### containing this type and return it
  ### If object type is not provided, return the first (supposedly unique)
  ### element of the MLDBM file
  my ($filename, $object_type) = @_;
  my $class_holder = undef;
  my $to_compress = 0;

  warn (";\n; ", &AlphaDate, " importing $object_type from $filename\n")
      if ($verbose >= 1);
  
  unless (-e  $filename) {
    if (-e "$filename.gz")  { 
      ### temporarily uncompress file
      warn ("; uncompressing file $filename.gz\n")
	  if ($verbose >= 1);
      system "gunzip $filename.gz";
      $to_compress = 1;
    } else {
      die "Error: file $filename does not exist\n";
    }
  }
  
  my %db = ();
  tie %db, 'MLDBM', $filename
      || die "Error: cannot tie to $filename\n";
  
  if ($object_type) {
    ### get the class holder containing the requested object type
    $class_holder = exists $db{$object_type} ? $db{$object_type} : undef;
  } else {
    ### take the first entry of the tied hash
    my ($object_type, $class_holder) = each %db;
  }


  untie %db;

  if ($to_compress) {
    warn ("; compressing file $filename\n")
	if ($verbose >= 1);
    system "gzip $filename";
  }

  ### check that the object is well a class holder
  die "failed to load a class holder from file $filename\n" 
      unless ref($class_holder) eq "classes::ClassFactory";

  warn ("; ", &AlphaDate, " $object_type imported\n")
      if ($verbose >= 1);
  
  return $class_holder;
}


sub LoadClass {
  ### usage : &LoadClass($file,$class);
  ### loads a class previously stored using the routine DumpClass;
  ### This routine :
  ### - loads the data in the file
  ### - updates the class variable @_objects
  ### - updates ID and name indexes
  ### This is necessary to use methods like get_objects, get_object, ...
  my ($file,$class) = @_;
  my $short_class = $class;
  my $to_compress = 0;
  $short_class =~ s/.*:://g;
  unless (-r $file) {
    if (-e $file.".gz") {
      warn "; uncompressing file $file.gz\n";
      system "gunzip $file.gz";
      $to_compress = 1;
    } else {
      die "Error: cannot load class $class from file '$file'\n";
    }
  }
  warn (";\n; ", &AlphaDate, " loading class $class ", "from file $file\n")
      if ($verbose >= 1);
  do $file;
  $class->set_objects(@{$short_class});
  $class->set_count($#{$short_class} + 1);
  $class->index_ids();
  $class->index_names();
  warn (";\n; ", &AlphaDate, 
	" class $short_class loaded\n") 
    if ($verbose >= 1);
  if ($to_compress) {
    warn "compressing file $file\n";
    system "gzip $file";
  }
}

sub LoadAMAZEfile {
  ### load a data file in the .obj format (created with lib/load_classes.pl)
  ### usage : &LoadAMAZEfile($data_file);
  ### usage : @classes = &LoadAMAZEfile($data_file); 
  ###         returns the list of instantiated classes
  my ($data_file) = @_;
  my $self = "";
  my $class = "";
  my %classes = ();

  warn (";\n; ", &AlphaDate,
	" Loading file $data_file\n")
    if ($verbose >= 1);
  unless (open DATA, $data_file) {
    warn ";WARNING: could not load data file '$data_file'\n";
    return 0;
  }
  
  while (<DATA>) {
    chomp;
    if (/^(classes::\S*)\s+\{/) { ### new object
      $class = $1;
      $classes{$class}++;
      $self = $class->fast_new();
    } elsif (/\}/) { ### object is terminated
      printf STDERR "; loaded object %s in class %s\n", $self->get_attribute("id"),$class if ($verbose >= 3);
      $class = $self = "";
    } elsif (($self) && (/^\s+(\S*)\s*(.*)/)) {
      $key = $1;
      $value = $2;
      $self->new_attribute_value($key,$value);
    }
  }
  close DATA;
  
  foreach $class (keys %classes) {
    $class->index_ids();
    $class->index_names();
  }
  
  return keys %classes;
}


1;

