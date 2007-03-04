# RSATWS.pm - rsa-tools web services module

package RSATWS;

use SOAP::Lite;

my $RSAT = $0; $RSAT =~ s|/public_html/+web_services/.*||;
my $SCRIPTS = $RSAT.'/perl-scripts';
my $TMP = $RSAT.'/public_html/tmp';

sub retrieve_seq {
    my ($self, $args_ref) = @_;
    my %args = %$args_ref;
    my $output_choice = $args{"output"};
    unless ($output_choice) {
	$output_choice = 'both';
    }
    my $command = $self->retrieve_seq_cmd(%args);
    my $stderr = `$command 2>&1 1>/dev/null`;
    if ($stderr) {
	die SOAP::Fault -> faultcode('Server.ExecError') -> faultstring("Execution error: $stderr\ncommand: $command");
    }
    my $result = `$command`;
    my $tmp_outfile = `mktemp $TMP/retrieve-seq.XXXXXXXXXX`;
    open TMP, ">".$tmp_outfile or die "cannot open temp file ".$tmp_outfile."\n";
    print TMP $result;
    close TMP;
    if ($output_choice eq 'server') {
	return {'command' => $command, 
		'server' => $tmp_outfile};
    } elsif ($output_choice eq 'client') {
	return {'command' => $command,
		'client' => $result};
    } elsif ($output_choice eq 'both') {
	return {'server' => $tmp_outfile,
		'command' => $command, 
		'client' => $result};
    }
}

sub retrieve_seq_cmd {
    my ($self,%args) = @_;
    my $organism = $args{"organism"};
    my $noorf = $args{"noorf"};
    my $from = $args{"from"};
    my $to = $args{"to"};

    ## List of query genes
    my $query_ref = $args{"query"};
    my $query = "";
    if ($query_ref) {
	my @query = @{$query_ref};
	foreach $q (@query) {
	    $q =~s/\'//g;
	    $q =~s/\"//g;
	}
	$query = " -q '";
	$query .= join "' -q '", @query;
	$query .= "'";
    }

    my $feattype = $args{"feattype"};
    my $type = $args{"type"};
    my $format = $args{"format"};
    my $all = $args{"all"};
    my $lw = $args{"lw"};
    my $label = $args{"label"};
    my $label_sep = $args{"label_sep"};
    my $nocom = $args{"nocom"};
    my $repeat = $args{'repeat'};
    my $imp_pos = $args{'imp_pos'};

    my $command = "$SCRIPTS/retrieve-seq";

    if ($organism) {
      $organism =~ s/\'//g;
      $organism =~ s/\"//g;
      $command .= " -org '".$organism."'";
    }
    if ($query) {
	$command .= $query;
    }

    if ($noorf == 1) {
	$command .= " -noorf";
    }
    if ($from =~ /\d/) { ## This is to make the difference between unspecified parameter and value 0
	$from =~ s/\'//g;
	$from =~ s/\"//g;
	$command .= " -from '".$from."'";
    }
    if ($to =~ /\d/) { ## This is to make the difference between unspecified parameter and value 0
	$to =~ s/\'//g;
	$to =~ s/\"//g;
	$command .= " -to '".$to."'";
    }
    if ($feattype) {
	$feattype =~ s/\'//g;
	$feattype =~ s/\"//g;
	$command .= " -feattype '".$feattype."'";
    }
    if ($type) {
	$type =~ s/\'//g;
	$type =~ s/\"//g;
	$command .= " -type '".$type."'";
    }
    if ($format) {
	$format =~ s/\'//g;
	$format =~ s/\"//g;
	$command .= " -format '".$format."'";
    }
    if ($all == 1) {
	$command .= " -all";
    }
    if (defined($lw)) {
	$lw =~ s/\'//g;
	$lw =~ s/\"//g;
	$command .= " -lw '".$lw."'";
    }
    if ($label) {
	$label =~ s/\'//g;
	$label =~ s/\"//g;
	$command .= " -label '".$label."'";
    }
    if ($label_sep) {
	$label_sep =~ s/\'//g;
	$label_sep =~ s/\"//g;
	$command .= " -labelsep '".$label_sep."'";
    }
    if ($nocom == 1) {
	$command .= " -nocom";
    }
    if ($repeat == 1) {
	$command .= " -rm";
    }
    if ($imp_pos == 1) {
	$command .= " -imp_pos";
    }

    return $command;
}

sub purge_seq {
    my ($self, $args_ref) = @_;
    my %args = %$args_ref;
    my $output_choice = $args{"output"};
    unless ($output_choice) {
	$output_choice = 'both';
    }
    my $command = $self->purge_seq_cmd(%args);
    my $stderr = `$command 2>&1 1>/dev/null`;
#    if ($stderr) {
#	die SOAP::Fault -> faultcode('Server.ExecError') -> faultstring("Execution error: $stderr\ncommand: $command");
#    }
    my $result = `$command`;
    my $tmp_outfile = `mktemp $TMP/purge-seq.XXXXXXXXXX`;
    open TMP_OUT, ">".$tmp_outfile or die "cannot open temp file ".$tmp_outfile."\n";
    print TMP_OUT $result;
    close TMP_OUT;
    if ($output_choice eq 'server') {
	return {'command' => $command, 
		'server' => $tmp_outfile};
    } elsif ($output_choice eq 'client') {
	return {'command' => $command,
		'client' => $result};
    } elsif ($output_choice eq 'both') {
	return {'server' => $tmp_outfile,
		'command' => $command,
		'client' => $result};
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
	$tmp_infile =~ s/\'//g;
	$tmp_infile =~ s/\"//g;
	chomp $tmp_infile;
    }
    my $format = $args{"format"};
    my $match_length = $args{"match_length"};
    my $mismatch = $args{"mismatch"};
    my $str = $args{"str"};
    my $delete = $args{"delete"};
    my $mask_short = $args{"mask_short"};

    my $command = "$SCRIPTS/purge-sequence";

    if ($format) {
      $format =~ s/\'//g;
      $format =~ s/\"//g;
      $command .= " -format '".$format."'";
    }

    if ($str) {
	if ($str == 1 || $str == 2) {
	    $command .= " -".$str."str";
	} else {
	    die "str value must 1 or 2";
	}
    }

    if ($delete == 1) {
      $command .= " -del";
    }

    if (defined($mask_short)) {
      $mask_short =~ s/\'//g;
      $mask_short =~ s/\"//g;
      $command .= " -mask_short '".$mask_short."'";
    }

    if (defined($match_length)) {
	$match_length =~ s/\'//g;
	$match_length =~ s/\"//g;
	$command .= " -ml '".$match_length."'";
    }
    if (defined($mismatch)) {
	$mismatch =~ s/\'//g;
	$mismatch =~ s/\"//g;
	$command .= " -mis '".$mismatch."'";
    }

    $command .= " -i '".$tmp_infile."'";

    return $command;
}

sub oligo_analysis {
    my ($self, $args_ref) = @_;
    my %args = %$args_ref;
    my $output_choice = $args{"output"};
    unless ($output_choice) {
	$output_choice = 'both';
    }
    my $command = $self->oligo_analysis_cmd(%args);
#    my $stderr = `$command 2>&1 1>/dev/null`;
#    if ($stderr) {
#	die SOAP::Fault -> faultcode('Server.ExecError') -> faultstring("Execution error: $stderr\ncommand: $command");
#    }
    my $result = `$command`;
    my $tmp_outfile = `mktemp $TMP/oligo.XXXXXXXXXX`;
    open TMP_OUT, ">".$tmp_outfile or die "cannot open temp file ".$tmp_outfile."\n";
    print TMP_OUT $result;
    close TMP_OUT;
    if ($output_choice eq 'server') {
	return {'command' => $command, 
		'server' => $tmp_outfile};
    } elsif ($output_choice eq 'client') {
	return {'command' => $command,
		'client' => $result};
    } elsif ($output_choice eq 'both') {
	return {'server' => $tmp_outfile,
		'command' => $command, 
		'client' => $result};
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
    }
    chomp $tmp_infile;
    my $format = $args{"format"};
    my $length = $args{"length"};
    my $organism = $args{"organism"};
    my $background = $args{"background"};
    my $stats = $args{"stats"};
    my $noov = $args{"noov"};
    my $str = $args{"str"};
    my $sort = $args{"sort"};
    my $lth = $args{"lth"};

    my $command = "$SCRIPTS/oligo-analysis";

    if ($format) {
      $format =~ s/\'//g;
      $format =~ s/\"//g;
      $command .= " -format '".$format."'";
    }

    if ($organism) {
      $organism =~ s/\'//g;
      $organism =~ s/\"//g;
      $command .= " -org '".$organism."'";
    }

    if ($background) {
      $background =~ s/\'//g;
      $background =~ s/\"//g;
      $command .= " -bg '".$background."'";
    }

    if ($stats) {
      $stats =~ s/\'//g;
      $stats =~ s/\"//g;
      $command .= " -return '".$stats."'";
    }

    if ($noov == 1) {
      $command .= " -noov";
    }

    if ($str) {
	if ($str == 1 || $str == 2) {
	    $command .= " -".$str."str";
	} else {
	    die "str value must 1 or 2";
	}
    }

    if ($sort == 1) {
      $command .= " -sort";
    }

    if ($lth) {
      $lth =~ s/\'//g;
      $lth =~ s/\"//g;
      @lth = split / /, $lth;
      $command .= " -lth '".$lth[0]."' '".$lth[1]."'";
    }

    if (defined($length)) {
	$length =~ s/\'//g;
	$length =~ s/\"//g;
	$command .= " -l '".$length."'";
    }

    $command .= " -i '".$tmp_infile."'";

    return $command;
}

sub dna_pattern {
    my ($self, $args_ref) = @_;
    my %args = %$args_ref;
    my $output_choice = $args{"output"};
    unless ($output_choice) {
	$output_choice = 'both';
    }
    my $command = $self->dna_pattern_cmd(%args);
    my $stderr = `$command 2>&1 1>/dev/null`;
    if ($stderr) {
	die SOAP::Fault -> faultcode('Server.ExecError') -> faultstring("Execution error: $stderr\ncommand: $command");
    }
    my $result = `$command`;
    my $tmp_outfile = `mktemp $TMP/dna_pattern.XXXXXXXXXX`;
    open TMP_OUT, ">".$tmp_outfile or die "cannot open temp file ".$tmp_outfile."\n";
    print TMP_OUT $result;
    close TMP_OUT;
    if ($output_choice eq 'server') {
	return {'command' => $command, 
		'server' => $tmp_outfile};
    } elsif ($output_choice eq 'client') {
	return {'command' => $command,
		'client' => $result};
    } elsif ($output_choice eq 'both') {
	return {'server' => $tmp_outfile,
		'command' => $command, 
		'client' => $result};
    }
}

sub dna_pattern_cmd {
    my ($self, %args) =@_;
    if ($args{"sequence"}) {
	my $sequence = $args{"sequence"};
	chomp $sequence;
	$tmp_infile = `mktemp $TMP/dna_pattern.XXXXXXXXXX`;
	open TMP_IN, ">".$tmp_infile or die "cannot open temp file ".$tmp_infile."\n";
	print TMP_IN $sequence;
	close TMP_IN;
    } elsif ($args{"tmp_infile"}) {
	$tmp_infile = $args{"tmp_infile"};
    }
    chomp $tmp_infile;
    my $format = $args{"format"};
    my $pattern = $args{"pattern"};
    my $subst = $args{"subst"};
    my $id = $args{"id"};
    my $origin = $args{"origin"};
    my $noov = $args{"noov"};
    my $str = $args{"str"};
    my $sort = $args{"sort"};
    my $th = $args{"th"};

    my $command = "$SCRIPTS/dna-pattern";

    if ($format) {
      $format =~ s/\'//g;
      $format =~ s/\"//g;
      $command .= " -format '".$format."'";
    }

    if ($pattern) {
      $pattern =~ s/\'//g;
      $pattern =~ s/\"//g;
      $command .= " -p '".$pattern."'";
    }

    if (defined($subst)) {
      $subst =~ s/\'//g;
      $subst =~ s/\"//g;
      $command .= " -subst '".$subst."'";
    }

    if ($id) {
      $id =~ s/\'//g;
      $id =~ s/\"//g;
      $command .= " -id '".$id."'";
    }

    if ($noov == 1) {
      $command .= " -noov";
    }

    if ($str) {
	if ($str == 1 || $str == 2) {
	    $command .= " -".$str."str";
	} else {
	    die "str value must 1 or 2";
	}
    }

    if ($sort == 1) {
      $command .= " -sort";
    }

    if (defined($th)) {
      $th =~ s/\'//g;
      $th =~ s/\"//g;
      $command .= " -th '".$th."'";
    }

    if ($origin) {
	$origin =~ s/\'//g;
	$origin =~ s/\"//g;
	$command .= " -origin '".$origin."'";
    }

    $command .= " -i '".$tmp_infile."'";

    return $command;
}

sub gene_info {
    my ($self, $args_ref) = @_;
    my %args = %$args_ref;
    my $output_choice = $args{"output"};
    unless ($output_choice) {
	$output_choice = 'both';
    }
    my $command = $self->gene_info_cmd(%args);
    my $result = `$command`;
    my $stderr = `$command 2>&1 1>/dev/null`;
    if ($stderr) {
	die SOAP::Fault -> faultcode('Server.ExecError') -> faultstring("Execution error: $stderr\ncommand: $command");
    }
    my $tmp_outfile = `mktemp $TMP/oligo.XXXXXXXXXX`;
    open TMP_OUT, ">".$tmp_outfile or die "cannot open temp file ".$tmp_outfile."\n";
    print TMP_OUT $result;
    close TMP_OUT;
    if ($output_choice eq 'server') {
	return {'command' => $command, 
		'server' => $tmp_outfile};
    } elsif ($output_choice eq 'client') {
	return {'command' => $command,
		'client' => $result};
    } elsif ($output_choice eq 'both') {
	return {'server' => $tmp_outfile,
		'command' => $command, 
		'client' => $result};
    }
}

sub gene_info_cmd {
  my ($self, %args) =@_;

  ## List of queries
  my $query_ref = $args{"query"};
  my $query = "";
  if ($query_ref) {
    my @query = @{$query_ref};
    foreach $q (@query) {
	$q =~s/\'//g;
	$q =~s/\"//g;
    }
    $query = " -q '";
    $query .= join "' -q '", @query;
    $query .= "'";
  }

  my $command = "$SCRIPTS/gene-info";

  if ($args{organism}) {
      $args{organism} =~ s/\'//g;
      $args{organism} =~ s/\"//g;
      $command .= " -org '".$args{organism}."'";
  }
  if ($query) {
    $command .= $query;
  }
  if ($args{full} == 1) {
    $command .= " -full";
  }
  if ($args{noquery} == 1){
    $command .= " -noquery";
  }
  if ($args{descr} == 1) {
    $command .= " -descr";
  }
  if ($args{feattype}) {
      $args{feattype} =~ s/\'//g;
      $args{feattype} =~ s/\"//g;
      $command .= " -feattype '".$args{feattype}."'";
  }

  return $command;
}

sub supported_organisms {
    my ($self, $args_ref) = @_;
    my %args = %$args_ref;
    my $output_choice = $args{"output"};
    unless ($output_choice) {
	$output_choice = 'both';
    }
    my $command = $self->supported_organisms_cmd(%args);
    my $result = `$command`;
    my $stderr = `$command 2>&1 1>/dev/null`;
    if ($stderr) {
	die SOAP::Fault -> faultcode('Server.ExecError') -> faultstring("Execution error: $stderr\ncommand: $command");
    }
    my $tmp_outfile = `mktemp $TMP/oligo.XXXXXXXXXX`;
    open TMP_OUT, ">".$tmp_outfile or die "cannot open temp file ".$tmp_outfile."\n";
    print TMP_OUT $result;
    close TMP_OUT;
    if ($output_choice eq 'server') {
	return {'command' => $command, 
		'server' => $tmp_outfile};
    } elsif ($output_choice eq 'client') {
	return {'command' => $command,
		'client' => $result};
    } elsif ($output_choice eq 'both') {
	return {'server' => $tmp_outfile,
		'command' => $command, 
		'client' => $result};
    }
}

sub supported_organisms_cmd {
  my ($self, %args) =@_;

  my $command = "$SCRIPTS/supported-organisms";

  if ($args{format}) {
      $args{format} =~ s/\'//g;
      $args{format} =~ s/\"//g;
      $command .= " -format '".$args{format}."'";
  }
  if ($args{taxon}) {
      $args{taxon} =~ s/\'//g;
      $args{taxon} =~ s/\"//g;
      $command .= " -taxon '".$args{taxon}."'";
  }

  return $command;
}

1;
