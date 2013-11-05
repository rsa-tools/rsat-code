###############################################################
#
# Manipulation of Position-Specific Scoring Matrices (PSSM)
#
package RSAT::table;

use RSAT::GenericObject;
use RSAT::util;
use RSAT::stats;
@ISA = qw( RSAT::GenericObject );

=pod

=head1 NAME

    RSAT::table

=head1 DESCRIPTION

Generic class for manipuating tables. 

=head1 METHODS

=cut


################################################################

=pod

=item B<new()>

Create an empty table.

=cut
sub new {
    my ($class, %args) = @_;
    my $table = bless {
	nrow=>0,
	ncol=>0
	}, $class;
    return $table;
}


################################################################

=pod

=item B<transpose($nrow, $ncol, @table)>

transpose a table (rows become columns and reciprocally)

=cut
sub transpose {
    my ($nrow, $ncol, @table) = @_;
    my @transposed = ();
    for my $c (1..$ncol) {
	for my $r (1..$nrow) {
	    $transposed[$r][$c] = $table[$c][$r];
	}
    }    
    return @transposed;
}


################################################################

=pod

=item B<col_freq>

Return the frequency per column. Computation is cached (computed only
once for the life of the object).

=cut
sub col_freq {
  my ($self) = @_;
  unless ($self->get_attribute("col_freq_calculated")) {
    $self->calc_col_freq();
  }
  return $self->get_attribute("col_freq");
}

################################################################

=pod

=item B<calc_col_freq>

Calculate the frequency per column, i.e. the cell counts divided by
sum per column.

=cut
sub calc_col_freq {
  my ($self) = @_;
  my @col_sum = $self->col_sum();;
  my @col_freq = ();
  my $ncol = $self->ncol();
  my $nrow = $self->nrow();
  my @table = $self->getTable();
  &RSAT::message::Info(join("\t", "Calculating frequencies per column","rows:".$nrow, "columns:".$ncol))
    if ($main::verbose >= 3);
  foreach my $c (0..($ncol-1)) {
    for my $r (0..($nrow-1)) {
      if ($col_sum[$c] > 0) {
	$col_freq[$c][$r] = $table[$c][$r]/$col_sum[$c];
      } else {
	$col_freq[$c][$r] = 0;
      }
    }
  }
  @{$self->{col_freq}} = @col_freq;
  $self->force_attribute("col_freq_calculated", 1);
}


################################################################

=pod

=item B<row_freq>

Return the frequency per row. Computation is cached (computed only
once for the life of the object).

=cut
sub row_freq {
  my ($self) = @_;
  unless ($self->get_attribute("row_freq_calculated")) {
    $self->calc_row_freq();
  }
  return $self->get_attribute("row_freq");
}

################################################################

=pod

=item B<calc_row_freq>

Calculate the frequency per row, i.e. the cell counts divided by
sum per row.

=cut
sub calc_row_freq {
  my ($self) = @_;
  my @row_sum = $self->row_sum();;
  my @row_freq = ();
  my $ncol = $self->ncol();
  my $nrow = $self->nrow();
  my @table = $self->getTable();
  &RSAT::message::Info(join("\t", "Calculating frequencies per row","rows:".$nrow, "columns:".$ncol))
    if ($main::verbose >= 3);
  foreach my $r (0..($nrow-1)) {
    for my $c (0..($ncol-1)) {
      if ($row_sum[$r] > 0) {
	$row_freq[$c][$r] = $table[$c][$r]/$row_sum[$r];
      } else {
	$row_freq[$c][$r] = 0;
      }
    }
  }
  @{$self->{row_freq}} = @row_freq;
  $self->force_attribute("row_freq_calculated", 1);
}


################################################################

=pod

=item B<col_sum>

Return the sum of each column. Computation is cached (computed only
once for the life of the object).

=cut
sub col_sum {
  my ($self) = @_;
  unless ($self->get_attribute("col_sum_calculated")) {
    $self->calc_col_sum();
  }
  return $self->get_attribute("col_sum");
}

################################################################

=pod

=item B<calc_col_sum>

Calculate the sum of each column.

=cut
sub calc_col_sum {
    my ($self) = @_;
    my @col_sum = ();
    my $ncol = $self->ncol();
    my $nrow = $self->nrow();
    my @table = $self->getTable();
    &RSAT::message::Info(join("\t", "Calculating sum per column","rows:".$nrow, "columns:".$ncol))
	if ($main::verbose >= 5);
    foreach my $c (0..($ncol-1)) {
	my $col_sum = 0;
	for my $r (0..($nrow-1)) {
	    $col_sum += $table[$c][$r];
	}
	push @col_sum, $col_sum;
#	&RSAT::message::Debug("column", $c, "sum", $col_sum)  if ($main::verbose >= 10);
    }
    @{$self->{col_sum}} = @col_sum;
    $self->force_attribute("col_sum_calculated", 1);
}

################################################################

=pod

=item B<col_max>

Return the max of each column. Computation is cached (computed only
once for the life of the object).

=cut
sub col_max {
  my ($self) = @_;
  unless ($self->get_attribute("col_max_calculated")) {
    $self->calc_col_max();
  }
  return $self->get_attribute("col_max");
}

################################################################

=pod

=item B<calc_col_max>

Calculate the max of each column.

=cut
sub calc_col_max {
    my ($self) = @_;
    my @col_max = ();
    my $ncol = $self->ncol();
    my $nrow = $self->nrow();
    my @table = $self->getTable();
    &RSAT::message::Info(join("\t", "Calculating max per column","rows:".$nrow, "columns:".$ncol))
	if ($main::verbose >= 4);
    foreach my $c (0..($ncol-1)) {
	my $col_max;
	for my $r (0..($nrow-1)) {
	    $col_max = &RSAT::stats::checked_max($col_max, $table[$c][$r]);
	}
	push @col_max, $col_max;
#	warn join ("\t", ";", "column", $c, "max", $col_max), "\n"  if ($main::verbose >= 5);
    }
    @{$self->{col_max}} = @col_max;
    $self->force_attribute("col_max_calculated", 1);
}

################################################################

=pod

=item B<col_min>

Return the min of each column. Computation is cached (computed only
once for the life of the object).

=cut
sub col_min {
  my ($self) = @_;
  unless ($self->get_attribute("col_min_calculated")) {
    $self->calc_col_min();
  }
  return $self->get_attribute("col_min");
}

################################################################

=pod

=item B<calc_col_min>

Calculate the min of each column.

=cut
sub calc_col_min {
    my ($self) = @_;
    my @col_min = ();
    my $ncol = $self->ncol();
    my $nrow = $self->nrow();
    my @table = $self->getTable();
    &RSAT::message::Info(join("\t", "Calculating min per column","rows:".$nrow, "columns:".$ncol))
	if ($main::verbose >= 4);
    foreach my $c (0..($ncol-1)) {
	my $col_min;
	for my $r (0..($nrow-1)) {
	    $col_min = &RSAT::stats::checked_min($col_min, $table[$c][$r]);
	}
	push @col_min, $col_min;
#	warn join ("\t", ";", "column", $c, "min", $col_min), "\n"  if ($main::verbose >= 5);
    }
    @{$self->{col_min}} = @col_min;
    $self->force_attribute("col_min_calculated", 1);
}


################################################################

=pod

=item B<row_sum>

Return the sum of each row. Computation is cached (computed only
once for the life of the object).

=cut
sub row_sum {
  my ($self) = @_;
  unless ($self->get_attribute("row_sum_calculated")) {
    $self->calc_row_sum();
  }
  return $self->get_attribute("row_sum");
}

################################################################

=pod

=item B<calc_row_sum>

Calculate the sum of each row.

=cut
sub calc_row_sum {
    my ($self) = @_;
    my @row_sum = ();
    my $ncol = $self->ncol();
    my $nrow = $self->nrow();
    my @table = $self->getTable();
    &RSAT::message::Info(join("\t", "Calculating sum per row","rows:".$nrow, "columns:".$ncol))
	if ($main::verbose >= 5);
    foreach my $r (0..($nrow-1)) {
	my $row_sum = 0;
	for my $c (0..($ncol-1)) {
	    $row_sum += $table[$c][$r];
	}
	push @row_sum, $row_sum;
#	&RSAT::message::Debug( "row", $r, "sum", $row_sum)  if ($main::verbose >= 10);
    }
    @{$self->{row_sum}} = @row_sum;
    $self->force_attribute("row_sum_calculated", 1);
}

################################################################

=pod

=item B<force_row_sum>

Force the row sums to take specified values rather than calculating
them from the table itself.

=cut
sub force_row_sum {
  my ($self, @row_sum) = @_;
  my $nrow = $self->nrow();
  unless (scalar(@row_sum) == $nrow) {
    &RSAT::error::FatalError(join("\t", "RSAT::table::force_row_sum", 
				  "number of specified values", scalar(@row_sum),
				  "differs from the number of rows", $nrow));
  }
  @{$self->{row_sum}} = @row_sum;
  $self->force_attribute("row_sum_calculated", 1);
}

################################################################

=pod

=item B<force_col_sum>

Force the column sums to take specified values rather than calculating
them from the table itself.

=cut
sub force_col_sum {
  my ($self, @col_sum) = @_;
  my $ncol = $self->ncol();
  unless (scalar(@col_sum) == $ncol) {
    &RSAT::error::FatalError(join("\t", "RSAT::table::force_col_sum", 
				  "number of specified values", scalar(@col_sum),
				  "differs from the number of cols", $ncol));
  }
  @{$self->{col_sum}} = @col_sum;
  $self->force_attribute("col_sum_calculated", 1);
}


################################################################

=pod

=item B<row_mean>

Return the mean of each row. Computation is cached (computed only
once for the life of the object).

=cut
sub row_mean {
  my ($self) = @_;
  unless ($self->get_attribute("row_mean_calculated")) {
    $self->calc_row_mean();
  }
  return $self->get_attribute("row_mean");
}

################################################################

=pod

=item B<calc_row_mean>

Calculate the mean of each row.

=cut
sub calc_row_mean {
    my ($self) = @_;
    my $ncol = $self->ncol();
    my $nrow = $self->nrow();
    &RSAT::message::Info(join("\t", "Calculating mean per row","rows:".$nrow, "columns:".$ncol))
	if ($main::verbose >= 4);
    my @row_mean = ();
    my @row_sum = $self->row_sum();
    foreach my $r (0..($nrow-1)) {
	my $row_mean = $row_sum[$r]/$ncol;
	push @row_mean, $row_mean;
    }
    @{$self->{row_mean}} = @row_mean;
    $self->force_attribute("row_mean_calculated", 1);
}

################################################################

=pod

=item B<col_mean>

Return the mean of each column. Computation is cached (computed only
once for the life of the object).

=cut
sub col_mean {
  my ($self) = @_;
  unless ($self->get_attribute("col_mean_calculated")) {
    $self->calc_col_mean();
  }
  return $self->get_attribute("col_mean");
}

################################################################

=pod

=item B<calc_col_mean>

Calculate the mean of each column.

=cut
sub calc_col_mean {
    my ($self) = @_;
    my @col_mean = ();
    my $ncol = $self->ncol();
    my $nrow = $self->nrow();
    &RSAT::message::Info(join("\t", "Calculating mean per column","rows:".$nrow, "columns:".$ncol))
	if ($main::verbose >= 4);
    my @col_sum = $self->col_sum();
    foreach my $c (0..($ncol-1)) {
	my $col_mean = $col_sum[$c]/$nrow;
	push @col_mean, $col_mean;
    }
    @{$self->{col_mean}} = @col_mean;
    $self->force_attribute("col_mean_calculated", 1);
}

################################################################

=pod

=item B<row_max>

Return the max of each row. Computation is cached (computed only
once for the life of the object).

=cut
sub row_max {
  my ($self) = @_;
  unless ($self->get_attribute("row_max_calculated")) {
    $self->calc_row_max();
  }
  return $self->get_attribute("row_max");
}

################################################################

=pod

=item B<calc_row_max>

Calculate the max of each row.

=cut
sub calc_row_max {
    my ($self) = @_;
    my @row_max = ();
    my $ncol = $self->ncol();
    my $nrow = $self->nrow();
    my @table = $self->getTable();
    &RSAT::message::Info(join("\t", "Calculating max per row","rows:".$nrow, "columns:".$ncol))
	if ($main::verbose >= 4);
    foreach my $r (0..($nrow-1)) {
	my $row_max;
	for my $c (0..($ncol-1)) {
	    $row_max = &RSAT::stats::checked_max($row_max, $table[$c][$r]);
	}
	push @row_max, $row_max;
#	warn join ("\t", ";", "column", $c, "max", $col_max), "\n"  if ($main::verbose >= 5);
    }
    @{$self->{row_max}} = @row_max;
    $self->force_attribute("row_max_calculated", 1);
}

################################################################

=pod

=item B<row_min>

Return the min of each row. Computation is cached (computed only
once for the life of the object).

=cut
sub row_min {
  my ($self) = @_;
  unless ($self->get_attribute("row_min_calculated")) {
    $self->calc_row_min();
  }
  return $self->get_attribute("row_min");
}

################################################################

=pod

=item B<calc_row_min>

Calculate the min of each row.

=cut
sub calc_row_min {
    my ($self) = @_;
    my @row_min = ();
    my $ncol = $self->ncol();
    my $nrow = $self->nrow();
    my @table = $self->getTable();
    &RSAT::message::Info(join("\t", "Calculating min per row","rows:".$nrow, "columns:".$ncol))
	if ($main::verbose >= 4);
    foreach my $r (0..($nrow-1)) {
	my $row_min;
	for my $c (0..($ncol-1)) {
	    $row_min = &RSAT::stats::checked_min($row_min, $table[$c][$r]);
	}
	push @row_min, $row_min;
#	warn join ("\t", ";", "column", $c, "min", $col_min), "\n"  if ($main::verbose >= 5);
    }
    @{$self->{row_min}} = @row_min;
    $self->force_attribute("row_min_calculated", 1);
}



################################################################

=pod

=item B<addColumn()>

Add a new column to the table

=cut
sub addColumn {
    my ($self,@new_col) = @_;
#	push @{$self->{table}}, [@new_col];
    
    
    &RSAT::message::Debug("Table: adding column", join (" ", @new_col)) if ($main::verbose >= 5);
    
    ## Update number of columns
    my $ncol = $self->ncol()+1;
    $self->force_attribute("ncol", $ncol);
    &RSAT::message::Debug("Table: updating number of columns", $self->ncol()) if ($main::verbose >= 5);
    
    ## update number of rows
    my $column_size = scalar(@new_col);
    if ($column_size >= $self->nrow()) {
	&RSAT::message::Debug ("Table: updating number of rows", $column_size) if ($main::verbose >= 5);
	$self->force_attribute("nrow", scalar(@new_col));
    }
    
    ## update table content
    for my $r (0..$#new_col) {
	${$self->{table}}[$ncol-1][$r] = $new_col[$r];
    }
}


################################################################

=pod

=item B<addRow(@new_row)>

Add a new row to the table

=cut
sub addRow {
    my ($self,@new_row) = @_;
    
    ## Update number of rows
    my $nrow = $self->nrow()+1;
	$self->force_attribute("nrow", $nrow);
    &RSAT::message::Debug("Table: updating number of rows", $self->nrow()) if ($main::verbose >= 5);
    
    ## update number of colmuns
    my $row_size = scalar(@new_row);
    if ($row_size >= $self->ncol()) {
	&RSAT::message::Debug("Table: updating number of columns", $row_size) if ($main::verbose >= 5);
	$self->force_attribute("ncol", scalar(@new_row));
    }
    
    ## update table content
    for my $c (0..$#new_row) {
	${$self->{table}}[$c][$nrow-1] = $new_row[$c];
    }
}


################################################################

=pod

=item B<setAlphabet(@alphabet)>

Specify the alphabet (i.e. the list of valid letters) for the table.

=cut
sub setAlphabet {
    my ($self, @new_alphabet) = @_;
    @{$self->{alphabet}} = @new_alphabet;

    ## update the number of columns
    $self->force_attribute("nrow", scalar(@new_alphabet));
#    &RSAT::message::Debug("&RSAT::table::setAlphabet()", "new alphabet", $self->getAlphabet())) if ($main::verbose >= 10);
}


################################################################

=pod

=item B<setAlphabet_uc(@alphabet)>

Same as setAlphabet(), but first converts the alphabet to uppercases,
to ensure case-insensitivvity.

=cut
sub setAlphabet_uc {
    my ($self, @new_alphabet) = @_;

    ## Convert alphabet to uppercases
    for my $i (0..$#new_alphabet) {
	$new_alphabet[$i] = uc($new_alphabet[$i]);
    }

    $self->setAlphabet(@new_alphabet);
}


################################################################

=pod

=item B<setAlphabet_lc(@alphabet)>

Same as setAlphabet(), but first converts the alphabet to lowercases,
to ensure case-insensitivvity.

=cut
sub setAlphabet_lc {
    my ($self, @new_alphabet) = @_;

    ## Convert alphabet to uppercases
    for my $i (0..$#new_alphabet) {
	$new_alphabet[$i] = lc($new_alphabet[$i]);
    }
    $self->setAlphabet(@new_alphabet);
}


################################################################

=pod

=item B<getAlphabet()>

Return the list of valid letters for the table

=cut
sub getAlphabet {
    my ($self) = @_;
    return @{$self->{alphabet}};
}

################################################################

=pod

=item B<addIndexedRow($index, @new_row)>

Add a new row and append its symbol to the alphabet.

=cut
sub addIndexedRow {
  my ($self, $index, @new_row) = @_;
  $self->addRow(@new_row);
  $self->push_attribute("alphabet", $index);

  ## Check that the row contains numerical values
  for my $value (@new_row) {
    unless (&RSAT::util::IsReal($value)) {
      &RSAT::error::FatalError("Invalid matrix row: all values after the index ($index) should be numbers");
    }
  }
}


################################################################

=pod

=item B<ncol()>

Return number of columns

=cut
sub ncol {
    my ($self) = @_;
    return $self->get_attribute("ncol");
}

################################################################

=pod

=item B<nrow()>

Return number of rows

=cut
sub nrow {
    my ($self) = @_;
    return $self->get_attribute("nrow");
}

################################################################

=pod

=item B<size()>

Return table size

=cut
sub size {
    my ($self) = @_;
    return ($self->nrow(), $self->ncol());
}


################################################################

=pod

=item B<getTable()>

Return the whole table (a vector of vectors)

=cut
sub getTable {
    my ($self) = @_;
    return @{$self->{table}};
}

################################################################

=pod

=item B<setTable($nrow, $ncol, @table)>

Specify the whole table

=cut
sub setTable {
    my ($self,$nrow, $ncol, @table) = @_;
    $self->force_attribute("nrow", $nrow);
    $self->force_attribute("ncol", $ncol);
    @{$self->{table}} = @table;
}


################################################################

=pod

=item B<setCell($row, $col, $value)>

Specify the content of a single cell. 

=cut
sub setCell {
    my ($self,$row, $col, $value) = @_;
#    warn join("\t", "Setting cell", 
#	      "row", $row, 
#	      "column", $col, 
#	      "value", $value), "\n" 
#		  if (main::verbose >= 10); 
    ${$self->{table}}[$col-1][$row-1] = $value;
}


################################################################

=pod

=item B<getCell($row, $col)>

Return the content of a single cell. 

=cut
sub getCell {
    my ($self,$row, $col) = @_;
    return ${$self->{table}}[$col-1][$row-1];
}



################################################################

=pod

=item B<readFromFile($file, $format)>

Read a table from a file

=cut
sub readFromFile {
    my ($self, $file, $format, %args) = @_;
    if ($format =~ /tab/i) {
	$self->_readFromTabFile($file, %args);
    } else {
	&main::FatalError("Invalid format for reading table\t$format");
    }
    
    ## Check that the table contains at least one row and one col
    if (($self->nrow() > 0) && ($self->ncol() > 0)) {
      &RSAT::message::Info(join("\t", "Table read", 
				"nrow = ".$self->nrow(),
				"ncol = ".$self->ncol(),
			       )) if ($main::verbose >= 2);
    } else {
	&main::FatalError("The file $file does not seem to contain a table in format $format. Please check the file format and contents.");
    }
}


################################################################

=pod
    
=item B<_readFromTabFile($file)>

Read a table from a tab-delimited file. This method is called by the
method C<readFromFile($file, "tab")>.

=cut
sub _readFromTabFile {
    my ($self, $file, %args) = @_;
    &RSAT::message::Info(join("\t", "Reading table from tab file",$file)) if ($main::verbose >= 2);

    ## open input stream
    my ($in, $dir) = &RSAT::util::OpenInputFile($file);
    my $current_table_nb = 0;
    my $l=0;
    ## read header
    my $header = "";
    if ($args{header}) {
      do {
	$header = <$in>;
      } until ($header !~ /^;/);
      chomp ($header);
      $l++;
    }
    while (<$in>) {
	$l++;
	next unless (/\S/);
	next if (/^;/);
	chomp();
	if (/^#/) {
	  $header = $1;
	  next;
	}
	s/\r//;
	if (/^\s*(\S+)\s+/) {
	    my @fields = split /\t/, $_;

#	    warn join("\t", @fields), "\n" if ($main::verbose >= 10);

	    ## residue associated to the row
	    my $residue = shift @fields;

#	    warn join("\t", "line", $l, "row name", $residue), "\n" if ($main::verbose >= 10);

	    ## skip the | between residue and numbers
	    shift @fields unless &main::IsReal($fields[0]);	
            
	    $self->addIndexedRow($residue, @fields);
	}
    }
    close $in if ($file);

    ## Assign the header
    if ($header) {
      @header = split "\t", $header;
      shift @header;
      $self->push_attribute("header", @header);
    }
}


################################################################

=pod

=item B<init>

Initialize the table.

=cut

sub init {
    my ($self) = @_;
    
    ## initialize the table
    my $nrow = $self->nrow();
    my $ncol = $self->ncol();
    warn "Initializing the table $nrow rows, $ncol columns\n" if ($main::verbose >= 5);
    foreach my $r (1..$nrow) {
	foreach my $c (1..$ncol) {
	    $self->setCell($r,$c,0);
	}
    }
}

################################################################

=pod

=item B<_printSeparator($ncol)>

Print a separator between header/footer and table

=cut

sub _printSeparator {
    my ($self, $ncol) = @_;
    my $sep = $self->get_attribute("sep") || "\t";
    my $col_width = $self->get_attribute("col_width");
    my $separator = "";
    
    if (($col_width) && ($col_width < 6)){
	$separator .= ";-";
    } else {
	$separator .= "; -----";
    }
    $separator .= $sep."|-";
    for $c (0..($ncol-1)) {
	if ($col_width) {
	    $separator .= "-"x($col_width-1);
	} else {
	    $separator .= "-"x7;
	}
	$separator .= "|";
    }
    $separator .= "\n";
    return $separator;
}



################################################################

=pod

=item B<get_row()>

Return a row of the table as a list.

=cut
sub get_row {
    my ($self, $row_nb) = @_;
    my @row = ();
    my $ncol = $self->ncol();
    my $nrow = $self->nrow();
    &RSAT::error::FatalError(join ("\t", $row_nb, "invalid row number (max=".$nrow.")")) if ($row_nb > $nrow);
    &RSAT::error::FatalError(join ("\t", $row_nb, "row numbers cannot be negative")) if ($row_nb < 1);
    for my $c (0..($ncol-1)) {
	push @row, ${$self->{table}}[$c][$row_nb-1];
    }
    return @row;
}

################################################################

=pod

=item B<get_column($col_nb, $nrow, @table)>

Return a column of the table as a list.

=cut
sub get_column {
    my ($self, $col_nb) = @_;
    my @col = ();
    my $ncol = $self->ncol();
    my $nrow = $self->nrow();
    &RSAT::error::FatalError(join ("\t", $col_nb, "invalid column number (max=".$ncol.")")) if ($col_nb > $ncol);
    &RSAT::error::FatalError(join ("\t", $col_nb, "column numbers cannot be negative")) if ($col_nb < 1);
    for my $r (0..($nrow-1)) {
	push @col, ${$self->{table}}[$col_nb-1][$r];
    }
    return @col;
}


################################################################

=pod

=item B<toString(sep=>$sep, col_width=>$col_width, $decimals=>$decimals)>

Return a string representing the table. 

Arguments:

=over 4

=item I<title>

Text printed above the table.

=item I<corner>

Text printed in the left top corner of the tabe.

=item I<decimals>

Number of decimals.

=back

=cut
sub toString {
    my ($self, %args) = @_;
    my $to_print = "";


    ## Set formatting parameters provided in arguments as table attribute
    foreach my $key ("sep", "col_width", "decimals") {
	if (defined($args{$key})) {
	    $self->force_attribute($key, $args{$key});
	}
    }

    ## Format for the table entries
    my $sep = $self->get_attribute("sep") || "\t";
    my $col_width = $self->get_attribute("col_width");
    my $decimals = $self->get_attribute("decimals");

    ## Calculate number width
    my $number_width = 0;
    if ($col_width) {
	$number_width = $col_width - 1;
    }

    ## Number of decimal digits 
    unless ($decimals) {
	$decimals = $number_width - 2;
    }
    
    ################################################################
    ## Print a table
    my @table = @{$self->{table}};
    my @alphabet = $self->getAlphabet();

    my $ncol = $self->ncol();
    my $nrow = $self->nrow();
    
    ## Header for the table
    my @header = $self->get_attribute("header");
    my $header_given = 1;
    if (scalar(@header) <= 0) {
	for my $c (0..($ncol-1)) {
	    my $pos = $c+1;
	    push @header, $pos;
	  }
    }

    if ($self->get_attribute("margins")) {
      push @header, "|", "r.mean", "r.sum", "r.max", "r.min";
    }

    if ($main::verbose >= 1) {
	$to_print .= ";\n";
	if ($args{"title"}) {
	  $to_print .= "; $args{title}\n";
	}
	if (($col_width) && ($col_width < 6)) {
	    $to_print .= ";";
	} else {
	    $to_print .= "; $args{corner}";
	}
	foreach my $header_item (@header) {
	    if ($col_width) {
	      $to_print .= sprintf "%${col_width}s", $header_item;
	    } else {
	      $to_print .= $sep;
	      $to_print .= $header_item;
	    }
	}
	$to_print .= "\n";
	
#	$to_print .= $self->_printSeparator($ncol, $to_print);
    }
    
    ## Print the table
    my @row_mean;
    my @row_sum;
    my @row_max;
    my @row_min;
    if ($self->get_attribute("margins")) {
      @row_mean = $self->row_mean();
      @row_sum = $self->row_sum();
      @row_max = $self->row_max();
      @row_min = $self->row_min();
    }
    for $a (0..$#alphabet) {
	my @row = $self->get_row($a+1, $ncol, @table);
	if ($self->get_attribute("margins")) {
	  push @row, "|", sprintf("%.3g", $row_mean[$a]), $row_sum[$a], $row_max[$a], $row_min[$a];
	}
	$to_print .= $self->_printTableRow($alphabet[$a], @row);
    }

    ################################################################
    ##Print column statistics
    if ($self->get_attribute("margins")) {
	$prefix_letter = substr($type, 0, 1) || "c";
	$to_print .= $self->_printSeparator($ncol+4, $to_print);

	## Mean per column
	my @col_mean = $self->col_mean($nrow, $ncol, @table);
	push @col_mean, "|", &RSAT::stats::mean(@col_mean);
	$to_print .= $self->_printTableRow("; ".$prefix_letter.".mea", @col_mean);

	## Sum per column
	my @col_sum = $self->col_sum($nrow, $ncol, @table);
	push @col_sum, "|", "", &RSAT::stats::sum(@col_sum);
	$to_print .= $self->_printTableRow("; ".$prefix_letter.".sum", @col_sum);

	## Maximum per column
	my @col_max = $self->col_max($nrow, $ncol, @table);
	push @col_max, "|", "", "", &RSAT::stats::checked_max(@col_max);
	$to_print .= $self->_printTableRow("; ".$prefix_letter.".max", @col_max);

	## Minimum per column
	my @col_min = $self->col_min($nrow, $ncol, @table);
	push @col_min, "|", "", "", "", &RSAT::stats::checked_min(@col_min);
	$to_print .= $self->_printTableRow("; ".$prefix_letter.".min", @col_min);
    }

    return $to_print;
}

################################################################

=pod

=item B<_printTableRow($row_name, @values)>

Print a row for the table output.

=cut

sub _printTableRow {
  my ($self, $row_name, @values) = @_;
  my $row_string = $row_name;
  my $ncol = scalar(@values);

  #	my ($decimals, $sep, $col_width, $number_width) = $self->_get_format;
	
  my $decimals = $self->get_attribute("decimals");

  my $sep = $self->get_attribute("sep") || "\t";
  my $col_width = $self->get_attribute("col_width");
    
  ## Format for the table entries
  my $number_width = 0;
  if ($col_width) {
    $number_width = $col_width - 1;
  }
  unless ($decimals) {
    if ($type eq "counts") {
      $decimals = 0;
    } else {
      $decimals = $number_width - 2;
    }
  }
    
  ## Print the table row
  #    $row_string .= $sep."|";
  for $c (0..($ncol-1)) {
    my $value = $values[$c];
    if ($col_width) {
      my $value_format = "%${number_width}s";
      if (&main::IsReal($value)) {
	if ($type eq "counts") {
	  $value_format = "%${number_width}d";
	} else {
	  $value_format= "%${number_width}.${decimals}f";
	}
      }
      $row_string .= sprintf " ${value_format}", $value;
    } elsif (($decimals > 0) && (&main::IsReal($value))) {
      $row_string .= $sep;
      $row_string .= sprintf "%.${decimals}g", $value;
    } else {
      $row_string .= $sep.$value;
    }
  }
  $row_string .= "\n";
  return $row_string;
}


################################################################

=pod

=item B<_get_format()>

Return the parameters for formatting output columns.

Usage:  my ($decimals, $sep, $col_width, $number_width) = $self->_get_format;

=cut
sub _get_format {
    my ($self) = @_;

    ## Column separator
    my $sep = $self->get_attribute("sep") || "\t";

    ## Column width
    my $col_width = $self->get_attribute("col_width");

    ## Format for the table entries
    my $number_width = 0;
    if ($col_width) {
	$number_width = $col_width - 1;
    }

    ## Number of decimals
    my $decimals = $self->get_attribute("decimals");
    unless ($decimals) {
	if ($type eq "counts") {
	    $decimals = 0;
	} else {
	    $decimals = $number_width - 2;
	}
    }

    return ($decimals, $sep, $col_width, $number_width);
}


return 1;


__END__

=pod

=back

