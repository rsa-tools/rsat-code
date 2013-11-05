################################################################
### DATABASE MANAGEMENT

package classes::DatabaseManagement;
{
  ### convert a query into a key for indexes
  ### - convert all letters to uppercases
  ### - convert white spaces and dashes into underscores
  sub get_index_key {
    my ($self, $query) = @_;
    my $index_key = uc($query);
    my @blanks = ();
    push @blanks, " ";
    push @blanks, "-";
    foreach $blank (@blanks) {
      $index_key =~ s/$blank/_/g;
    }
    return $index_key;
  }


  sub split_pseudo_pointer {
    ### a pseudo_pointer is a string that refers to a database object (e.g. classes::Gene::b0001)
    ### it consists in two parts :
    ### - the class (e.g. classes::Gene)
    ### - the object ID (e.g. b00001)
    ### this allows a system-independent cross-referencing between classes
    ###
    ### usage: ($class, $key) = classes::DatabaseManagement->split_pseudo_pointer($pseudo_pointer);
    my ($database, $pseudo_pointer) = @_;
    my $class;
    my $key;
    if ($pseudo_pointer =~ /(.*)::/) {
      $class = $1;
#      unless (defined($database->{$class})) {
#	die "Error: pseudo_pointer '$pseudo_pointer': class '$class' not found in the database\n";
#      }
      $key = "$'";
    } else {
      die "Error: invalid database pseudo_pointer $pseudo_pointer\n";
    }
    return  ($class, $key);
  }


}

return 1;
