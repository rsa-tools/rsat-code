### the package classes::Index allows to manipulate a special kind of index, where 
### each key is associated to a list of values rather than a single one
### This allows for example to deal with homonymy (several IDs associated to the same name)
package classes::Index;
{
    ### creator
    sub new {
	my ($class) = @_;
	$index = bless  {
	}, $class;
	return $index;
    }

   ### accessors

    ### given a key, get the list of associated values 
  sub get_values {
      my ($index, $key) = @_;
      return @{$index->{$key}};
  }

  #### returns the first value associated to the key
  sub get_first_value {
      my ($index, $key) = @_;
      my @values = @{$index->{$key}};
      return $values[0];
  }

  #### returns true if the index contains the specified key.
  sub contains {
      my ($index, $key) = @_;
#      if (defined(@{$index->{$key}})) {
      if (defined($index->{$key})) {
	  return 1;
      } else {
	  return 0;
      }
  }

  #### returns true if the index contain the specified key-value pair. 
  sub contains_pair {
      my ($index, $key, $value) = @_;
#      if (defined(@{$index->{$key}})) {
      if (defined($index->{$key})) {
	  foreach my $indexed_value ($index->get_value($key)) {
	      return 1 if ($value eq $indexed_value);
	  }
      } 
      return 0;
  }

  sub set_values {
      my ($index, $key, @values) = @_;
      @{$index->{$key}} = @values;
  }

  sub add_value {
    my ($index, $key, $new_value) = @_;
    my @values = $index->get_values($key);
    $index->set_values($key, sort (@values, $new_value));
  }

  ### return all index keys
    sub get_keys {
	my ($index) = @_;
	my %hash = %{$index};
	return keys(%{$index});
    }

  #### return the number of keys in the index
  sub get_size {
      my ($index) = @_;
      my @keys = $index->get_keys;
      return $#keys+1;
  }

  ### exports the index hash as a separate file, 
  ### which can be used for cross-reference searching
  ### usage: 
  ###   $index->export($format,$outfile);
  ### supported formats: MLDBM, tab
  sub export {
    my ($index, $format, $outfile) = @_;

    ### MLDBM export
    if (lc($format) eq 'mldbm') {
#      use MLDBM qw(DB_File Storable);
      die "Error: export called with empty file name\n"
	  unless ($outfile);
      warn (";\n; ", &main::AlphaDate, 
	    " exporting index in format MLDBM to file $outfile\n")
	  if ($main::verbose >= 1);
      tie %db, 'MLDBM', $outfile ||
	  die "Cannnot tie to $outfile\n";
      
      foreach my $key ($index->get_keys()) {
	my @values = $index->get_values($key);
	$db{$key} = \@values;
      }
      
      untie %db;
      warn ("; ", &main::AlphaDate, " export done\n")
	  if ($main::verbose >= 1);
      return;
      
      ### tab export
    } elsif (lc($format) eq 'tab') {
	die "Error: export called with empty file name\n"
	    unless ($outfile);
	warn (";\n; ", &main::AlphaDate, 
	      " exporting index in tab format to file $outfile\n")
	    if ($main::verbose >= 1);
	open INDEX, ">$outfile" ||
	    die "Cannnot write file $outfile\n";
	
	foreach my $key ($index->get_keys()) {
	    my @values = $index->get_values($key);
	    foreach my $value (@values) {
		print INDEX "$key\t$value\n";
	    }
	}
	
	close INDEX;
	warn ("; ", &main::AlphaDate, " export done\n")
	    if ($main::verbose >= 1);
	return;
	
	### unsupported format
    } else {
	die "Error: format '$format' is not supported for index export\n";
    }
  }

  sub load {
      #### Loads an index from a tab-delimited text file.
      #### The first column contains the key, the second the value.
      #### Lines starting with '--' are ignored. 

      #### The second and third arguments, specify whether the keys
      #### and values have to be converted to lower- or
      #### upper-cases. They can take the values "lc" (lowercases), "uc"
      #### (uppercases), or be ommited.

      my ($index, $infile, $standard_key, $standard_value, %args) = @_;
      
      die "Error: a file name should be specified for loading the index.\n"
	  unless $infile;

      unless ($infile =~ /|\s+$/) {
	  #### test the existenc of the file
	  die "Error: index file $infile does not exist.\n"
	      unless (-e $infile);
	  die "Error: cannot read index file $infile.\n"
	      unless (-r $infile);
      }
      open INDEX, $infile ||
	  die "Error: cannot open index file $infile.\n";
      while (<INDEX>) {
	  chomp;
	  s/\r//g;
	  next if (/^--/);
	  next unless (/\S/);
	  #### read key-value pair
	  if ($args{reverse}==1) {
	      ($value, $key) = split "\t";
	  } else {
	      ($key, $value) = split "\t";
	  }
	  #### standardize the key if required
	  if ($standard_key) {
	      $key = &main::standardize($key);
	  }
	  #### standardize the value if required
	  if ($standard_value) {
	      $value = &main::standardize($value);
	  }
	  if (($key) && ($value)) {
	      $index->add_value($key,$value);
	  }
#	  print STDERR join "\t", $standard_key, "key=$key", $standard_value, "value=$value", "\n";
      }
      close INDEX;
      return;
  }

  #### returns an index with the current data in the inverse mapping :
  #### keys become valuaes and values become keys
  #### usage:
  ####   $reversed_index = $index->reverse();
  sub reverse {
      my ($index) = @_;
      my $reversed_index = new classes::Index;
      foreach my $key ($index->get_keys) {
	  my @values = $index->get_values($key);
	  foreach my $value (@values) {
	      $reversed_index->add_value($value, $key);
	  }
      }
      return $reversed_index;
  }

}

return 1;
