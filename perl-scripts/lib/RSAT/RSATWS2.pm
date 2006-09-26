# RSATWS2.pm - rsa-tools web services module

package RSATWS2;

my $RSAT = $0; $RSAT =~ s|/public_html/+web_services/.*||;
my $SCRIPTS = $RSAT.'/perl-scripts';
my $TMP = $RSAT.'/public_html/tmp';

=pod

=head2 retrieve-seq

WS interface to the script B<I<retrieve-seq>>.

=head3 Usage

my ($server_file, $command, $result) = $service->retrieve_seq(SOAP::Data->name('param' => \SOAP::Data->value([parameters]))); 

=head3 Parameters

=over

=item I<organism>

organism name

=item I<...>

=back

=cut

sub retrieve_seq {
    my ($self, @args) = @_;
    my $return_choice = $args[0];
    my $command = $self->retrieve_seq_cmd(@args);
    my $result = `$command`;
    my $tmp_outfile = `mktemp $TMP/retrieve-seq.XXXXXXXXXX`;
    open TMP, ">".$tmp_outfile or die "cannot open temp file ".$tmp_outfile."\n";
    print TMP $result;
    close TMP;
    if ($return_choice eq 'file') {
	return {'command' => $command, 
		'file' => $tmp_outfile};
    } elsif ($return_choice eq 'result') {
	return {'command' => $command,
		'result' => $result};
    } elsif ($return_choice eq 'both') {
	return {'file' => $tmp_outfile,
		'command' => $command, 
		'result' => $result};
    }
}

=pod

=head2 retrieve_seq_cmd

Generates the appropriate retrieve-seq command given a hash of parameters.

=cut

sub retrieve_seq_cmd {
    my ($self,@args) = @_;
    my $organism = $args[1];
    my $noorf = $args[3];
    my $from = $args[4];
    my $to = $args[5];

    ## List of query genes
    my $query_ref = $args[2];
    my @query = @$query_ref;
    my $query = "-q ";
    $query .= join " -q ", @query;

    my $feattype = $args[6];
    my $type = $args[7];
    my $format = $args[8];
    my $all = $args[9];
    my $n = $args[10];
    my $label = $args[11];
    my $label_sep = $args[12];
    my $nocom = $args[13];

    my $command = "$SCRIPTS/retrieve-seq $organism $query -lw 0 $noorf $from $to $feattype $type $format $all $n $label $label_sep $nocom";
    return $command;
}

sub purge_seq {
    my ($self, $args_ref) = @_;
    my %args = %$args_ref;
    my $return_choice = $args{"return"};
    my $command = $self->purge_seq_cmd(%args);
    my $result = `$command`;
    my $tmp_outfile = `mktemp $TMP/purge-seq.XXXXXXXXXX`;
    open TMP_OUT, ">".$tmp_outfile or die "cannot open temp file ".$tmp_outfile."\n";
    print TMP_OUT $result;
    close TMP_OUT;
    if ($return_choice eq 'file') {
	return {'command' => $command, 
		'file' => $tmp_outfile};
    } elsif ($return_choice eq 'result') {
	return {'command' => $command,
		'result' => $result};
    } elsif ($return_choice eq 'both') {
	return {'file' => $tmp_outfile,
		'command' => $command, 
		'result' => $result};
    }
}

sub purge_seq_cmd {
    my ($self, %args) = @_;
    if ($args{"sequence"}) {
	my $sequence = $args{"sequence"};
	chomp $sequence;
	$tmp_infile = `mktemp $TMP/purge-seq.XXXXXXXXXX`;
	open TMP_IN, ">".$tmp_infile or die "cannot open temp file ".$tmp_infile."\n";
	print TMP_IN $sequence;
	close TMP_IN;
    } elsif ($args{"tmp_infile"}) {
	$tmp_infile = $args{"tmp_infile"};
	chomp $tmp_infile;
    }
    my $format = $args{"format"};
    my $match_length = $args{"match_length"};
    my $mismatch = $args{"mismatch"};
    my $str = $args{"str"};
    my $delete = $args{"delete"};
    my $mask_short = $args{"mask_short"};

    my $command = "$SCRIPTS/purge-sequence $format $match_length $mismatch $str $delete $mask_short -i $tmp_infile";
    return $command;
}

sub oligo_analysis {
    my ($self, $args_ref) = @_;
    my %args = %$args_ref;
    my $return_choice = $args{"return"};
    my $command = $self->oligo_analysis_cmd(%args);
    my $result = `$command`;
    my $tmp_outfile = `mktemp $TMP/oligo.XXXXXXXXXX`;
    open TMP_OUT, ">".$tmp_outfile or die "cannot open temp file ".$tmp_outfile."\n";
    print TMP_OUT $result;
    close TMP_OUT;
    if ($return_choice eq 'file') {
	return {'command' => $command, 
		'file' => $tmp_outfile};
    } elsif ($return_choice eq 'result') {
	return {'command' => $command,
		'result' => $result};
    } elsif ($return_choice eq 'both') {
	return {'file' => $tmp_outfile,
		'command' => $command, 
		'result' => $result};
    }
}

sub oligo_analysis_cmd {
    my ($self, %args) =@_;
    if ($args{"sequence"}) {
	my $sequence = $args{"sequence"};
	chomp $sequence;
	$tmp_infile = `mktemp $TMP/oligo.XXXXXXXXXX`;
	open TMP_IN, ">".$tmp_infile or die "cannot open temp file ".$tmp_infile."\n";
	print TMP_IN $sequence;
	close TMP_IN;
    } elsif ($args{"tmp_infile"}) {
	$tmp_infile = $args{"tmp_infile"};
	chomp $tmp_infile;
    }
    my $format = $args{"format"};
    my $length = $args{"length"};
    my $organism = $args{"organism"};
    my $background = $args{"background"};
    my $stats = $args{"stats"};
    my $noov = $args{"noov"};
    my $str = $args{"str"};
    my $sort = $args{"sort"};
    my $lth = $args{"lth"};

    my $command = "$SCRIPTS/oligo-analysis $format $length $organism $background $stats $noov $str $sort $lth -i $tmp_infile";
    return $command;
}

=pod

=back

=cut

1;
