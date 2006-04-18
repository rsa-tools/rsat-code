###############################################################
#
# Manipulation of Position-Specific Scoring Matrices (PSSM)
#
package RSAT::table;

use RSAT::GenericObject;
use RSAT::util;
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

=item B<calc_col_sums>

Calculate the sum of each column

=cut
sub calc_col_sums {
    my ($self) = @_;
    my @col_sums = ();
    my $ncol = $self->ncol();
    my $nrow = $self->nrow();
    my @table = $self->getTable();
    warn join("\t", "; Calculating sum per column","rows:".$nrow, "columns:".$ncol),"\n"
	if ($main::verbose > 2);
    foreach my $c (0..($ncol-1)) {
	my $col_sum = 0;
	for my $r (0..($nrow-1)) {
	    $col_sum += $table[$c][$r];
	}
	push @col_sums, $col_sum;
#	warn join ("\t", ";", "column", $c, "sum", $col_sum), "\n"  if ($main::verbose >= 5);
    }
    @{$self->{col_sums}} = @col_sums;
}


################################################################
=pod

=item B<calc_row_sums>

Calculate the sum of each rowumn

=cut
sub calc_row_sums {
    my ($self) = @_;
    my @row_sums = ();
    my $ncol = $self->ncol();
    my $nrow = $self->nrow();
    my @table = $self->getTable();
    warn join("\t", "; Calculating sum per row","rows:".$nrow, "columns:".$ncol),"\n"
	if ($main::verbose > 2);
    for my $r (0..($nrow-1)) {
	my $row_sum = 0;
	foreach my $c (0..($ncol-1)) {
	    $row_sum += $table[$c][$r];
	}
	push @row_sums, $row_sum;
#	warn join ("\t", ";", "row", $r, "sum", $row_sum), "\n"  if ($main::verbose >= 5);
    }
    @{$self->{row_sums}} = @row_sums;
}


################################################################
=pod

=item B<addColumn()>

Add a new column to the table

=cut
sub addColumn {
    my ($self,@new_col) = @_;
#	push @{$self->{table}}, [@new_col];
    
    
    warn ("Table: adding column\t", join (" ", @new_col), "\n") if ($main::verbose >= 5);
    
    ## Update number of columns
    my $ncol = $self->ncol()+1;
    $self->force_attribute("ncol", $ncol);
    warn ("Table: updating number of columns\t", $self->ncol(), "\n") if ($main::verbose >= 5);
    
    ## update number of rows
    my $column_size = scalar(@new_col);
    if ($column_size >= $self->nrow()) {
	warn ("Table: updating number of rows\t", $column_size, "\n") if ($main::verbose >= 5);
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
    warn ("Table: updating number of rows\t", $self->nrow(), "\n") if ($main::verbose >= 5);
    
    ## update number of colmuns
    my $row_size = scalar(@new_row);
    if ($row_size >= $self->ncol()) {
	warn ("Table: updating number of columns\t", $row_size, "\n") if ($main::verbose >= 5);
	$self->force_attribute("ncol", scalar(@new_row));
    }
    
    ## update table content
    for my $c (0..$#new_row) {
	${$self->{table}}[$c][$nrow-1] = $new_row[$c];
    }
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

=item B<setAlphabet(@alphabet)>

Specify the alphabet (i.e. the list of valid letters) for the table.

=cut
sub setAlphabet {
    my ($self, @new_alphabet) = @_;

    @{$self->{alphabet}} = @new_alphabet;	

#    warn join("\t", "; Alphabet", $self->getAlphabet()), "\n" if ($main::verbose >= 10);
    
    ## update the number of columns
    $self->force_attribute("nrow", scalar(@alphabet));
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
	warn join("\t", "; Table read", 
		  "nrow = ".$self->nrow(),
		  "ncol = ".$self->ncol(),
		 ), "\n" if ($main::verbose >= 2);
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
    warn ("; Reading table from tab file\t",$file, "\n") if ($main::verbose >= 2);
    
#	($in, $dir) = &main::OpenInputFile($file);
#
    
    ## open input stream
    my ($in, $dir) = &RSAT::util::OpenInputFile($file);
#     my $in = STDIN;
#     if ($file) {
#	 open $in, $file;
# 	open INPUT, $file;
# 	$in = INPUT;
#     }
    my $current_table_nb = 0;
    my $l=0;
    ## read header
    if ($args{header}) {
	$header = <$in>;
	chomp ($header);
	@header = split "\t", $header;
	$self->push_attribute("header", @header);
	$l++
    }
    while (<$in>) {
	$l++;
	next unless (/\S/);
	s/\r//;
	chomp();
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
#    die;
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
    if ($main::verbose >= 1) {
	$to_print .= ";\n";
	if (($col_width) && ($col_width < 6)) {
	    $to_print .= ";";
	} else {
	    $to_print .= "; Table";
	}
#	$to_print .= $sep."|";
	for my $c (0..($ncol-1)) {
	    my $pos = $c+1;
	    if ($col_width) {
		$to_print .= sprintf "%${col_width}s", $pos;
	    } else {
		$to_print .= $sep;
		$to_print .= $pos;
	    }
	}
	$to_print .= "\n";
	
#	$to_print .= $self->_printSeparator($ncol, $to_print);
    }
    
    ## Print the table
    for $a (0..$#alphabet) {
	my @row = $self->get_row($a+1, $ncol, @table);
	$to_print .= $self->_printTableRow($alphabet[$a], @row);
    }

#    die "HELLO";
    ################################################################
    ##Print column statistics
    if ($self->get_attribute("margins")) {
	$prefix_letter = substr($type, 0, 1);
	$to_print .= $self->_printSeparator($ncol, $to_print);
	
	## Sum per column
#	my @col_sum = &col_sum($nrow, $ncol, @table);
	my @col_sum = $self->col_sum($nrow, $ncol, @table);
	push @col_sum, &main::sum(@col_sum);
	$to_print .= $self->_printTableRow("; ".$prefix_letter.".sum", @col_sum);
	
	## Maximum per column
#	my @col_max = &col_max($nrow, $ncol, @table);
	my @col_max = $self->col_max($nrow, $ncol, @table);
	push @col_max, &main::max(@col_max);
	$to_print .= $self->_printTableRow("; ".$prefix_letter.".max", @col_max);
	
	## Minimum per column
#	my @col_min = &col_min($nrow, $ncol, @table);
	my @col_min = $self->col_min($nrow, $ncol, @table);
	push @col_min, &main::min(@col_min);
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
	unless ($decimals) {
	    if ($type eq "counts") {
		$decimals = 0;
  	} else {
  	    $decimals = $number_width - 2;
  	}
    }
    my $sep = $self->get_attribute("sep") || "\t";
    my $col_width = $self->get_attribute("col_width");
    
    ## Format for the table entries
    my $number_width = 0;
    if ($col_width) {
  	$number_width = $col_width - 1;
    }
    
    ## Print the table row
#    $row_string .= $sep."|";
    for $c (0..($ncol-1)) {
	my $value = $values[$c];
	if ($col_width) {
	    my $value_format = "%${number_width}s";
	    if (&main::IsReal($value)){
		if ($type eq "counts") {
		    $value_format = "%${number_width}d";
		} else {
		    $value_format= "%${number_width}.${decimals}f";
		}
	    }
	    $row_string .= sprintf " ${value_format}", $value;
	} else {
	    $row_string .= $sep.$value;
	}
    }
    $row_string .= "\n";
    return $row_string
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

