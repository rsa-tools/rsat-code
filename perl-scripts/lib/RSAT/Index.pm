###############################################################
#
# Class Index
#
package RSAT::Index;

use RSAT::GenericObject;
use RSAT::error;
@ISA = qw( RSAT::GenericObject );

### class attributes

=pod

=head1 NAME

    RSAT::Index

=head1 DESCRIPTION

Manipulation an "index", i.e. a hashtable where each key is associated
to a list of values.

This allows for example to deal with homonymy (several IDs associated
to the same name).

=cut



=pod

=item get_values

Given a key, return the list of associated values.

=cut
sub get_values {
    my ($self, $key) = @_;
    my @values = ();
    if (defined($self->{$key})) {
	@values = @{$self->{$key}};
#    } else {
#	&RSAT::message::Warning("Index contains no value associated to key $key") if ($main::verbose >= 0);
    }

#    &RSAT::message::Debug("RSAT::index::get_values()", scalar(@values), "values for key", 
#			  "'".$key."'", join(";", @values)) if ($main::verbose >= 10);

    return @values;
}



=pod

=item get_first_value

returns the first value associated to the key

=cut
sub get_first_value {
    my ($self, $key) = @_;
    my @values = @{$self->{$key}};
    return $values[0];
}


=pod

=item contains

Returns true if the index contains the specified key.

=cut
sub contains {
    my ($self, $key) = @_;
    if (defined($self->{$key})) {
	return 1;
    } else {
	return 0;
    }
}

=pod

=item contains_pair

Returns true if the index contain the specified key-value pair.

Usage: if ($index->contains_key($key, $value)) { ... }

=cut
sub contains_pair {
    my ($self, $key, $value) = @_;
    if (defined($self->{$key})) {
	foreach my $indexed_value ($self->get_value($key)) {
	    return 1 if ($value eq $indexed_value);
	}
    }
    return 0;
}

=pod

=item set_values

Assigns a list of value to one key.

Usage: $index->set_values($key, @values);

=cut
sub set_values {
    my ($self, $key, @values) = @_;
    @{$self->{$key}} = @values;
}

=pod

=item add_value

Adds a value to the list of values associated to a given key. 

Usage: $index->add_value($key, $value);

=cut
sub add_value {
    my ($self, $key, $new_value) = @_;
    my @values = $self->get_values($key);

#    &RSAT::message::Debug("indexing", $key, $new_value) if ($main::verbose >= 10);

    $self->set_values($key, sort (@values, $new_value));
}

################################################################
=pod

=item get_keys

Return all the keys in the index.

Usage: my @keys = $index->get_keys();

=cut
sub get_keys {
    my ($self) = @_;
    my %hash = %{$self};
    return keys(%{$self});
}

################################################################
=pod

=item get_size

Returns the number of keys in the index.

Usage: my $index_size = $index->get_size();

return the number of keys in the index
=cut
sub get_size {
    my ($self) = @_;
    my @keys = $self->get_keys;
    return $#keys+1;
}

################################################################
=pod

=item export

Export the index in different formats.

Supported formats: mldbm (requires MLDBM probably not working anymore), tab.

Usage: $index->export("tab", $outfile);

=cut
sub export {
    my ($self, $format, $outfile) = @_;

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
	
	foreach my $key ($self->get_keys()) {
	    my @values = $self->get_values($key);
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
	
	foreach my $key ($self->get_keys()) {
	    my @values = $self->get_values($key);
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

################################################################
=pod

=item load

Loads an index from a tab-delimited text file.  The first column
contains the key, the second the value.  Lines starting with '--' are
ignored.

The second and third arguments, specify whether the keys and values
have to be converted to lower- or upper-cases. They can take the
values "lc" (lowercases), "uc" (uppercases), or be ommited.

Usage: $index->load($infile, $standard_key, $standard_value, %args);

Arguments:
    $infile: input file

    $standard_key: specifies whether the keys have to be standardized
                   (converted to lowercases)

    $standard_value: specifies whether the values have to be standardized
                   (converted to lowercases)

    $reverse: if true, values are taken from the first column of
              the file, and keys from the second column.

=cut
sub load {
    my ($self, $infile, $standard_key, $standard_value, %args) = @_;
    
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
	if ($reverse) {
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
	    $self->add_value($key,$value);
	}
#	  print STDERR join "\t", $standard_key, "key=$key", $standard_value, "value=$value", "\n";
    }
    close INDEX;
    return;
}

################################################################
=pod

=item reverse

Returns an index with the current data in the inverse mapping : keys
become values and values become keys.

Usage: $reversed_index = $index->reverse();

=cut
sub reverse {
    my ($self) = @_;
    my $reversed_index = new classes::Index;
    foreach my $key ($self->get_keys) {
	my @values = $self->get_values($key);
	foreach my $value (@values) {
	    $reversed_index->add_value($value, $key);
	}
    }
    return $reversed_index;
}




return 1;


__END__


