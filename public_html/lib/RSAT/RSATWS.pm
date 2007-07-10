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
    if ($lw) {
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
    if ($stderr) {
	die SOAP::Fault -> faultcode('Server.ExecError') -> faultstring("Execution error: $stderr\ncommand: $command");
    }
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
    }
    chomp $tmp_infile;
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
	    die "str value must be 1 or 2";
	}
    }

    if ($delete == 1) {
      $command .= " -del";
    }

    if ($mask_short) {
      $mask_short =~ s/\'//g;
      $mask_short =~ s/\"//g;
      $command .= " -mask_short '".$mask_short."'";
    }

    if ($match_length) {
	$match_length =~ s/\'//g;
	$match_length =~ s/\"//g;
	$command .= " -ml '".$match_length."'";
    }

    if ($mismatch) {
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
    my $stderr = `$command 2>&1 1>/dev/null`;
    if ($stderr) {
	die SOAP::Fault -> faultcode('Server.ExecError') -> faultstring("Execution error: $stderr\ncommand: $command");
    }
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
	    die "str value must be 1 or 2";
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

    if ($length) {
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
      chomp $pattern;
      $pattern =~ s/\'//g;
      $pattern =~ s/\"//g;
      $command .= " -pl '".$pattern."'";
    }

    if ($subst) {
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

    if ($th) {
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

sub convert_seq {
    my ($self, $args_ref) = @_;
    my %args = %$args_ref;
    my $output_choice = $args{"output"};
    unless ($output_choice) {
	$output_choice = 'both';
    }
    my $command = $self->convert_seq_cmd(%args);
    my $result = `$command`;
    my $stderr = `$command 2>&1 1>/dev/null`;
    if ($stderr) {
	die SOAP::Fault -> faultcode('Server.ExecError') -> faultstring("Execution error: $stderr\ncommand: $command");
    }
    my $tmp_outfile = `mktemp $TMP/convert-seq.XXXXXXXXXX`;
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

sub convert_seq_cmd {
  my ($self, %args) =@_;
  if ($args{"sequence"}) {
    my $sequence = $args{"sequence"};
    chomp $sequence;
    $tmp_infile = `mktemp $TMP/convert-seq.XXXXXXXXXX`;
    open TMP_IN, ">".$tmp_infile or die "cannot open temp file ".$tmp_infile."\n";
    print TMP_IN $sequence;
    close TMP_IN;
  } elsif ($args{"tmp_infile"}) {
    $tmp_infile = $args{"tmp_infile"};
    $tmp_infile =~ s/\'//g;
    $tmp_infile =~ s/\"//g;
  }
  chomp $tmp_infile;
  my $command = "$SCRIPTS/convert-seq";

  if ($args{from}) {
      $args{from} =~ s/\'//g;
      $args{from} =~ s/\"//g;
      $command .= " -from '".$args{from}."'";
  }
  if ($args{to}) {
      $args{to} =~ s/\'//g;
      $args{to} =~ s/\"//g;
      $command .= " -to '".$args{to}."'";
  }

  $command .= " -i '".$tmp_infile."'";

  return $command;
}

sub compare_classes {
  my ($self, $args_ref) = @_;
  my %args = %$args_ref;
  my $output_choice = $args{"output"};
  unless ($output_choice) {
    $output_choice = 'both';
  }
  my $command = $self->compare_classes_cmd(%args);
  my $stderr = `$command 2>&1 1>/dev/null`;
  if ($stderr) {
    die SOAP::Fault -> faultcode('Server.ExecError') -> faultstring("Execution error: $stderr\ncommand: $command");
  }
  my $result = `$command`;
  my $tmp_outfile = `mktemp $TMP/compare-classes.XXXXXXXXXX`;
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

sub compare_classes_cmd {
  my ($self,%args) = @_;
  #creation d'un fichier temporaire qui sera intégré dans la commande
  if ($args{"ref_classes"}) {
    my $reference = $args{"ref_classes"};
    chomp $reference;
    $tmp_ref = `mktemp $TMP/compare-ref-classes.XXXXXXXXXX`;
    open TMP_IN, ">".$tmp_ref or die "cannot open temp file ".$tmp_ref."\n";
    print TMP_IN $reference;
    close TMP_IN;
  }
  #idem
  if ($args{"query_classes"}) {
    my $query = $args{"query_classes"};
    chomp $query;
    $tmp_query = `mktemp $TMP/compare-query-classes.XXXXXXXXXX`;
    open TMP_IN, ">".$tmp_query or die "cannot open temp file ".$tmp_query."\n";
    print TMP_IN $query;
    close TMP_IN;
  }
  if ($args{"input_classes"}) {
    my $input = $args{"input_classes"};
    chomp $input;
    $tmp_input = `mktemp $TMP/compare-input-classes.XXXXXXXXXX`;
    open TMP_IN, ">".$tmp_input or die "cannot open temp file ".$tmp_input."\n";
    print TMP_IN $input;
    close TMP_IN;
  }

  # my $ref_classes = $args{"ref_classes"};
  # my $query_classes = $args{"query_classes"};
  my $return_fields = $args{"return_fields"};
  my $score_column = $args{"score_column"};
  my $upper_threshold_field = $args{"upper_threshold_field"};
  my $upper_threshold_value = $args{"upper_threshold_value"};
  my $lower_threshold_field = $args{"lower_threshold_field"};
  my $lower_threshold_value = $args{"lower_threshold_value"};
  my $sort = $args{"sort"};
  my $distinct = $args{"distinct"};
  my $triangle = $args{"triangle"};
  my $matrix = $args{"matrix"};

  my $command = "$SCRIPTS/compare-classes";

  #pas d'utilité directe de "nettoyage" de la commande sauf si l'on rajoute un elsif...  
  if ($tmp_ref) {
    $tmp_ref =~ s/\'//g;
    $tmp_ref =~ s/\"//g;
    chomp $tmp_ref;
    $command .= " -r '".$tmp_ref."'";
  }

  if ($tmp_query) {
    $tmp_query =~ s/\'//g;
    $tmp_query =~ s/\"//g;
    chomp $tmp_query;
    $command .= " -q '".$tmp_query."'";
  }

  if ($tmp_input) {
    $tmp_input =~ s/\'//g;
    $tmp_input =~ s/\"//g;
    chomp $tmp_input;
    $command .= " -i '".$tmp_input."'";
  }

  if ($return_fields) {
    $return_fields =~ s/\'//g;
    $return_fields =~ s/\"//g;
    $command .= " -return '".$return_fields."'";
  }

  if ($score_column) {
    $score_column =~ s/\'//g;
    $score_column =~ s/\"//g;
    $command .= " -sc '".$score_column."'";
  }

  if ($upper_threshold_field && $upper_threshold_value =~ /\d/)  {
    $upper_threshold_field =~ s/\'//g;
    $upper_threshold_field =~ s/\"//g;
    $upper_threshold_value =~ s/\'//g;
    $upper_threshold_value =~ s/\"//g;
    $command .= " -uth '".$upper_threshold_field."' '".$upper_threshold_value."'";
  }

  if ($lower_threshold_field && $lower_threshold_value =~ /\d/) {
    $lower_threshold_field =~ s/\'//g;
    $lower_threshold_field =~ s/\"//g;
    $lower_threshold_value =~ s/\'//g;
    $lower_threshold_value =~ s/\"//g;
    $command .= " -lth '".$lower_threshold_field."' '".$lower_threshold_value."'";
  }

  if ($sort) {
    $sort =~ s/\'//g;
    $sort =~ s/\"//g;
    $command .= " -sort '".$sort."'";
  }

  if ($distinct == 1) {
    $command .= " -distinct";
  }

  if ($triangle == 1) {
    $command .= " -triangle";
  }

  if ($matrix) {
    $matrix =~ s/\'//g;
    $matrix =~ s/\"//g;
    $command .= " -matrix '".$matrix."'";
  }

  return $command;
}

sub matrix_scan {
  my ($self, $args_ref) = @_;
  my %args = %$args_ref;
  my $output_choice = $args{"output"};
  unless ($output_choice) {
    $output_choice = 'both';
  }
  my $command = $self->matrix_scan_cmd(%args);
#  my $stderr = `$command 2>&1 1>/dev/null`;  ####cette gestion des erreurs est incompatible avec le fonctionnement de matrix-scan dans RSAT #######
  if ($stderr) {
    die SOAP::Fault -> faultcode('Server.ExecError') -> faultstring("Execution error: $stderr\ncommand: $command");
  }
  my $result = `$command`;
  my $tmp_outfile = `mktemp $TMP/matrix-scan.XXXXXXXXXX`;
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

sub matrix_scan_cmd {
  my ($self,%args) = @_;
  #creation d'un fichier temporaire qui sera intégré dans la commande
  if ($args{"sequence_file"}) {
    my $sequence = $args{"sequence_file"};
    chomp $sequence;
    $tmp_sequence_file = `mktemp $TMP/matscan-sequence_file.XXXXXXXXXX`;
    open TMP_IN, ">".$tmp_sequence_file or die "cannot open temp file ".$tmp_seqence_file."\n";
    print TMP_IN $sequence;
    close TMP_IN;
  }

    #idem
    if ($args{"matrix_file"}) {
      my $input_matrix = $args{"matrix_file"};
      chomp $input_matrix;
      $tmp_input_matrix = `mktemp $TMP/matscan-matrix_file.XXXXXXXXXX`;
      open TMP_IN, ">".$tmp_input_matrix or die "cannot open temp file ".$tmp_input_matrix."\n";
      print TMP_IN $input_matrix;
      close TMP_IN;
    }

    if ($args{"matrix_list"}) {
      my $input_list = $args{"matrix_list"};
      chomp $input_list;
      $tmp_input_list = `mktemp $TMP/matscan-matrix_list.XXXXXXXXXX`;
      open TMP_IN, ">".$tmp_input_list or die "cannot open temp file ".$tmp_input_list."\n";
      print TMP_IN $input_list;
      close TMP_IN;
    }

    if ($args{"background"}) {
      my $background = $args{"background"};
      chomp $background;
      $tmp_background = `mktemp $TMP/matscan-background.XXXXXXXXXX`;
      open TMP_IN, ">".$tmp_background or die "cannot open temp file ".$tmp_background ."\n";
      print TMP_IN $background;
      close TMP_IN; 
    }

    my $matrix_format = $args{"matrix_format"}; 
    my $top_matrices = $args{"top_matrices"};
    my $background_input = $args{"background_input"};
    my $background_window = $args{"background_window"};
    my $markov = $args{"markov"};
    my $background_pseudo = $args{"background_pseudo"};
    my $return_fields = $args{"return_fields"};
    my $upper_threshold_field = $args{"upper_threshold_field"};
    my $upper_threshold_value = $args{"upper_threshold_value"};
    my $lower_threshold_field = $args{"lower_threshold_field"};
    my $lower_threshold_value = $args{"lower_threshold_value"};
    my $both_strand = $args{"both_strand"};
    my $single_strand = $args{"single_strand"};

    my $command = "$SCRIPTS/matrix-scan";

 #pas d'utilite directe de "nettoyage" de la commande sauf si l'on rajoute un elsif...
    if ($tmp_sequence_file) {
      $tmp_sequence_file =~ s/\'//g;
      $tmp_sequence_file =~ s/\"//g;
      chomp $tmp_sequence_file;
      $command .= " -i '".$tmp_sequence_file."'";
    }

    if ($tmp_input_matrix) {
      $tmp_input_matrix =~ s/\'//g;
      $tmp_input_matrix =~ s/\"//g;
      chomp $tmp_input_matrix;
      $command .= " -m '".$tmp_input_matrix."'";
    }

    if ($matrix_format) {
      $matrix_format =~ s/\'//g;
      $matrix_format=~ s/\"//g;
      chomp $matrix_format;
      $command .= " -matrix_format '".$matrix_format."'";
    }

    if ($tmp_input_list) {
      $tmp_input_list =~ s/\'//g;
      $tmp_input_list =~ s/\"//g;
      $command .= " -mlist '".$tmp_input_list."'";
    }

    if ($top_matrices ) {
      $top_matrices  =~ s/\'//g;
      $top_matrices  =~ s/\"//g;
      $command .= " -top_matrices '".$top_matrices."'";
    }

    if ($tmp_background) {
      $tmp_background  =~ s/\'//g;
      $tmp_background  =~ s/\"//g;
      chomp $tmp_background;
      $command .= " -bgfile '".$tmp_background."'";
    }

     if ($background_input == 1 ) {
      $command .= " -bginput";
    }

     if ($background_window) {
      $background_window  =~ s/\'//g;
      $background_window=~ s/\"//g;
      $command .= " -window '".$background_window."'";
    }

     if ($markov =~ /\d/) {
      $markov  =~ s/\'//g;
      $markov =~ s/\"//g;
      $command .= " -markov '".$markov."'";
    }

     if ($background_pseudo) {
      $background_pseudo =~ s/\'//g;
      $background_pseudo =~ s/\"//g;
      $command .= " -bg_pseudo '".$background_pseudo."'";
    }

     if ($return_fields) {
      $return_fields =~ s/\'//g;
      $return_fields =~ s/\"//g;
      $command .= " -return '".$return_fields."'";
    }

    if ($upper_threshold_field && $upper_threshold_value =~ /\d/)  {
      $upper_threshold_field =~ s/\'//g;
      $upper_threshold_field =~ s/\"//g;
      $upper_threshold_value =~ s/\'//g;
      $upper_threshold_value =~ s/\"//g;
      $command .= " -uth '".$upper_threshold_field."' '".$upper_threshold_value."'";
    }

    if ($lower_threshold_field && $lower_threshold_value =~ /\d/) {
      $lower_threshold_field =~ s/\'//g;
      $lower_threshold_field =~ s/\"//g;
      $lower_threshold_value =~ s/\'//g;
      $lower_threshold_value =~ s/\"//g;

      $command .= " -lth '".$lower_threshold_field."' '".$lower_threshold_value."'";
    }

    if ($both_strand  == 1) {
      $command .= " -2str";
    }

    if ($single_strand == 1) {
       $command .= " -1str";
    }

    return $command;
}

1;
