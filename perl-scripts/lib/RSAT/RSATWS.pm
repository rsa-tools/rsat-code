# RSATWS.pm - rsa-tools web services module

package RSATWS;

use SOAP::Lite;

my $RSAT = $0; $RSAT =~ s|/public_html/+web_services/.*||;
my $SCRIPTS = $RSAT.'/perl-scripts';
my $TMP = $RSAT.'/public_html/tmp';

=pod

=head2 retrieve_seq

WS interface to the script B<I<retrieve-seq>>.

=head3 Usage
my %args = ('return' => $return_choice,
            'organism' => $organism,
            'query' => \@gene,  # an array in a hash has to be referenced (correct?)
            'noorf' => $noorf,
            'from' => $from,
            'to' => $to,
            'feattype' => $feattype,
            'type' => $type,
            'format' => $format,
            'all' => $all,
            'number' => $n,
            'label' => $label,
            'label_sep' => $label_sep,
            'nocom' => $nocom,
            'repeat' => $repeat);
my $results_ref = $soap->call('retrieve_seq' => 'request' => \%args)->result;
my %results = %$results_ref;
my $command = $results{'command'};
my $server_file = $results{'file'};
my $result = $results{'result'};

=head3 Parameters

=over

=item I<return>

return choice (file,result,both)

=item I<organism>

organism name

=item I<query>

array of gene names

=item I<noorf>

cut overlapping orf

=item I<from>

upstream coordinate

=item I<to>

downstream coordinate

=item I<feattype>

feature type (CDS,mRNA,tRNA,rRNA,scRNA)

=item I<type>

sequence type (upstream,downstream,orf,random)

=item I<format>

output format (IG,WC,raw,FastA))

=item I<all>

all genomic upstream regions

=item I<number>

number of sequences (only with type random)

=item I<label>

field used to label sequence (id,name,organism_name,sequence_type,current_from,current_to,ctg,orf_strand,reg_left,reg_right)

=item I<label_sep>

separator between label fields

=item I<nocom>

no comments, only identifier

=item I<repeat>

repeat masked

=back

=cut

sub retrieve_seq {
    my ($self, $args_ref) = @_;
    my %args = %$args_ref;
    my $return_choice = $args{"return"};
    my $command = $self->retrieve_seq_cmd(%args);
    my $stderr = `$command 2>&1 1>/dev/null`;
    if ($stderr) {
	die SOAP::Fault -> faultcode('Server.ExecError') -> faultstring("Execution error: $stderr");
    }
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
    my ($self,%args) = @_;
    my $organism = $args{"organism"};
    my $noorf = $args{"noorf"};
    my $from = $args{"from"};
    my $to = $args{"to"};

    ## List of query genes
    my $query_ref = $args{"query"};
    my @query = @{$query_ref};
    my $query = "-q ";
    $query .= join " -q ", @query;

    my $feattype = $args{"feattype"};
    my $type = $args{"type"};
    my $format = $args{"format"};
    my $all = $args{"all"};
    my $n = $args{"number"};
    my $label = $args{"label"};
    my $label_sep = $args{"label_sep"};
    my $nocom = $args{"nocom"};
    my $repeat = $args{'repeat'};

    my $command = "$SCRIPTS/retrieve-seq $organism $query -lw 0 $noorf $from $to $feattype $type $format $all $n $label $label_sep $nocom $repeat";
    return $command;
}

=pod

=head2 purge_seq

WS interface to the script B<I<purge-sequence>>.

=head3 Usage
my %args = ('return'=> $return_choice,
            'sequence'=> $sequence,
            'format'=> $format,
            'match_length'=> $match_length,
            'mismatch'=> $mismatch,
            'str'=> $str,
            'delete'=> $delete,
            'mask_short'=> $mask_short);
my $results_ref = $soap->call('purge_seq' => 'request' => \%args)->result;
my %results = %$results_ref;
my $command = $results{'command'};
my $server_file = $results{'file'};
my $result = $results{'result'};

=head3 Parameters

=over

=item I<return>

return choice (file,result,both)

=item I<sequence>

sequence(s)

=item I<format>

input sequence format

=item I<match_length>

match length

=item I<mismatch>

number of mismatches

=item I<str>

discard duplications on direct only or both strands (-1str,-2str)

=item I<delete>

delete repeats

=item I<mask_short>

mask sequences shorter than specified length

=back

=cut

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

=pod

=head2 purge_seq_cmd

Generates the appropriate purge-seq command given a hash of parameters.

=cut

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

=pod

=head2 oligo_analysis

WS interface to the script B<I<oligo-analysis>>.

=head3 Usage
my %args = ('return' => $return_choice, 
            'sequence' => $sequence, 
            'format' => $format,
            'length' => $length,
            'organism' => $organism, 
            'background' => $background,
            'stats' => $stats,
            'noov' => $noov,
            'str' => $str,
            'sort' => $sort,
            'lth' => $lth);
my $results_ref = $soap->call('purge_seq' => 'request' => \%args)->result;
my %results = %$results_ref;
my $command = $results{'command'};
my $server_file = $results{'file'};
my $result = $results{'result'};

=head3 Parameters

=over

=item I<return>

return choice (file,result,both)

=item I<sequence>

sequence(s)

=item I<format>

input sequence format

=item I<length>

oligo length

=item I<organism>

organism name

=item I<background>

background model (upstream,upstream-noorf,intergenic,upstreamL,input)

=item I<stats>

list of statistics to return (occ,mseq,freq,proba,ratio,zscore,like,pos,rank)

=item I<noov>

no overlapping

=item I<str>

oligo occurences on one or both strands (-1str,-2str)

=item I<sort>

sort oligos according to overrepresentation

=item I<lth>

lower threshold parameter (occ,occ_P,occ_E,occ_sig,observed_freq,exp_freq,zscore,mseq,ms_P,ms_E,ms_sig,ratio,rank) value

=back

=cut

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

=pod

=head2 oligo_analysis_cmd

Generates the appropriate oligo-analysis command given a hash of parameters.

=cut

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
