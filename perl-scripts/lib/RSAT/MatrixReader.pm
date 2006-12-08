###############################################################
#
# Class MatrixReader
#
package RSAT::MatrixReader;

use RSAT::GenericObject;
use RSAT::matrix;
use RSAT::error;
@ISA = qw( RSAT::GenericObject );

### class attributes

=pod

=head1 NAME

RSAT::MatrixReader

=head1 DESCRIPTION

A class for reading position-specific scoring matrices from different
formats.

=cut


################################################################
=pod

=item readFromFile($file, $format)

Read a matrix from a file

=cut
sub readFromFile {
    my ($file, $format, %args) = @_;
    my @matrices = ();

    if ((lc($format) eq "consensus") || ($format =~ /^wc/i)) {
	@matrices = _readFromConsensusFile($file);
    } elsif (lc($format) eq "assembly") {
	@matrices = _readFromAssemblyFile($file);
    } elsif (lc($format) eq "gibbs") {
	@matrices = _readFromGibbsFile($file);
    } elsif (lc($format) eq "tab") {
	@matrices = _readFromTabFile($file, %args);
    } elsif (lc($format) eq "MotifSampler") {
	@matrices = _readFromMotifSamplerFile($file);
    } elsif (lc($format) eq "meme") {
	@matrices = _readFromMEMEFile($file);
    } elsif (lc($format) eq "clustal") {
	@matrices = _readFromClustalFile($file);
    } else {
	&main::FatalError("&RSAT::matrix::readFromFile", "Invalid format for reading matrix\t$format");
    }

    ## Check that there was at least one matrix in the file
    if (scalar(@matrices) == 0) {
      &RSAT::message::Warning("File",  $file, "does not contain any matrix in format", $format);
    }  else {
      &RSAT::message::Info("Read ".scalar(@matrices)." matrices from file ", $file) if ($main::verbose >= 3);
    }

    foreach my $matrix (@matrices) {
      ## Check that each matrix contains at least one row and one col
      if (($matrix->nrow() > 0) && ($matrix->ncol() > 0)) {
	&RSAT::message::Info("Matrix read", 
			     "nrow = ".$matrix->nrow(),
			     "ncol = ".$matrix->ncol(),
			     "prior : ".join (" ", $matrix->getPrior()),
			    ) if ($main::verbose >= 3);
      } else {
	&RSAT::message::Warning("The file $file does not seem to contain a matrix in format $format. Please check the file format and contents.");
      }

      ## Replace undefined values by 0
      $matrix->treat_null_values();

    }

    return @matrices;
}


################################################################
=pod

=item _readFromGibbsFile($file)

Read a matrix from a gibbs file. This method is called by the method 
C<readFromFile($file, "gibbs")>.

=cut
sub _readFromGibbsFile {
    my ($file) = @_;
    
    ## open input stream
    my $in = STDIN;
    if ($file) {
	open INPUT, $file;
	$in = INPUT;
    }
    $in_matrix = 0;

    my @matrices = ();
    my $matrix = new RSAT::matrix();
    push @matrices, $matrix;

    my @matrix = ();
    my @alphabet = ();
    my $ncol = 0;
    my $nrow = 0;
    my $last_ncol = 0;
    my $last_nrow = 0;
    while (<$in>) {
	next unless (/\S/);
	s/\r//;
	chomp();
	if (/Information \(relative entropy\) contribution in tenth bits\:/) {
	    $in_matrix = 1;
	    # default nucletodide alphabet
	    $matrix->setAlphabet_lc("A","C","G","T");   
	    next;

	} elsif (/site/) {
	    ### Empty the previous matrix because it was not the definitive result
	    $in_matrix = 0;
	    @last_matrix = @matrix;
	    $last_nrow = $nrow;
	    $last_ncol = $ncol;
	    @matrix = ();
	    $nrow = 0;
	    $ncol = 0;

	    next;
	} elsif ((/^\s*POS/) && ($in_matrix)) {
	    s/\r//;
	    chomp;
	    @header = split " +";
	    @alphabet = @header[1..$#header-1];
#		$matrix->setAlphabet_lc(@alphabet);
	} elsif (/model map = (\S+); betaprior map = (\S+)/) {
	    $matrix->set_parameter("model.map", $1);
	    $matrix->set_parameter("betaprior.map", $2);
	} elsif (/MAP = (\S+)/) {
	    $matrix->set_parameter("MAP", $1);
	} elsif (/seed: (\S+)/) {
	    $matrix->set_parameter("seed", $1);
	} elsif (/^gibbs /) {
	    $matrix->set_parameter("command", $_);
	} elsif ($in_matrix) {
	    ## Add a column to the matrix (gibbs rows correspond to our columns)
	    s/\r//;
	    chomp;
	    s/^\s+//;
	    @fields = split " +";
	    @values = @fields[1..$#header-1];
	    $nrow = scalar(@values);
	    foreach my $v (0..$#values) {
		$values[$v] =~ s/^\.$/0/;
		$matrix[$ncol][$v] = $values[$v];
	    }
	    $ncol++;
	}
    }
    close $in if ($file);

    $matrix->setAlphabet_lc (@alphabet);
    $matrix->force_attribute("nrow", $last_nrow);
    $matrix->force_attribute("ncol", $last_ncol);
    $matrix->setMatrix ($last_nrow, $last_ncol, @last_matrix);
    return @matrices;
}


################################################################
=pod

=item _readFromConsensusFile($file)

Read a matrix from a consensus file. This method is called by the
method C<readFromFile($file, "consensus")>.

=cut
sub _readFromConsensusFile {
  my ($file) = @_;
  &RSAT::message::Info ("Reading matrix from consensus file", $file) if ($main::verbose >= 3);
  
  #    ($in, $dir) = &main::OpenInputFile($file);
  
  ## open input stream
  my $in = STDIN;
  if ($file) {
    open INPUT, $file;
    $in = INPUT;
  }
  my $current_matrix_nb = 0;
  my @matrices = ();
  my $matrix;
  my $command = "";

  my %prior = ();
  my $l = 0;
  while (<$in>) {
    $l++;
    next unless (/\S/);
    s/\r//;
    chomp();
    
    ## The following information (final cycle) is only exported when
    ## the number of cycles is automatic. I don't understand the
    ## reason for this. I need to ask Jerry. In the mean time, I
    ## always use the same information (THE LIST OF TOP MATRICES
    ## FROM EACH CYCLE).
    last if (/THE LIST OF MATRICES FROM FINAL CYCLE/);
    
    ## Read the command line
    if (/COMMAND LINE: /) {
      $command = $';		# '
      
      ## Start a new matrix (one consensus file contains several matrices)
    } elsif (/MATRIX\s(\d+)/) {
      $current_matrix_nb = $1;
      $matrix = new RSAT::matrix();
      push @matrices, $matrix;
      $matrix->set_attribute("number", $current_matrix_nb);
      $matrix->set_parameter("command", $command);
      $matrix->setPrior(%prior);
      next;

      ## Read prior frequency for one residue in the consensus header
    } elsif (/letter\s+\d:\s+(\S+).+prior frequency =\s+(\S+)/) {
      my $letter = lc($1);
      my $prior = $2;
      &RSAT::message::Info ("Prior from consensus file", $letter, $prior) if ($main::verbose >= 3);
      $prior{$letter} = $prior;
      
    } elsif ($current_matrix_nb >= 1) {
      
      ## Matrix content (counts) for one residue
      if (/^\s*(\S+)\s+\|/) {
	my @fields = split / +/, $_;
	## residue associated to the row
	my $residue = lc(shift @fields);
	
	## skip the | between residue and numbers
	shift @fields unless &main::IsReal($fields[0]);	
	
	$matrix->addIndexedRow($residue, @fields);
	
	## Sites used to build the matrix
      } elsif (/(\d+)\|(\d+)\s*\:\s*(-){0,1}(\d+)\/(\d+)\s+(\S+)/) {
	my $site_nb = $1;
	my $site_cycle = $2;
	my $site_strand = $3;
	my $site_seq_nb = $4;
	my $site_pos = $5;
	my $site_sequence = $6;
	$matrix->push_attribute("sequences", $site_sequence);
	&RSAT::message::Debug("line", $l, "site", $site_sequence) if ($main::verbose >= 4);
	
	## Other matrix parameters
      } elsif (/number of sequences = (\d+)/) {
	$matrix->set_parameter("nb.sequences", $1); 
      } elsif (/unadjusted information = (\S+)/) {
	$matrix->set_parameter("unadjusted.information", $1); 
      } elsif (/sample size adjusted information = (\S+)/) {
	$matrix->set_parameter("adjusted.information", $1); 
      } elsif (/ln\(p\-value\) = (\S+)   p\-value = (\S+)/) {
	$matrix->set_parameter("ln.Pval", $1); 
	$matrix->set_parameter("Pval", $2); 
      } elsif (/ln\(expected frequency\) = (\S+)   expected frequency = (\S+)/) {
	$matrix->set_parameter("ln.exp", $1); 
	$matrix->set_parameter("exp", $2); 
      }
    }
  }
  close $in if ($file);
  return @matrices;
}

################################################################
=pod

=item _readFromAssemblyFile($file)

Read a matrix from the output file of pattern-assembly. This method is
called by the method C<readFromFile($file, "assembly")>.

=cut
sub _readFromAssemblyFile {
  my ($file) = @_;
  &RSAT::message::Info ("Reading matrix from pattern-assembly file", $file) if ($main::verbose >= 3);
  
  #    ($in, $dir) = &main::OpenInputFile($file);
  
  ## open input stream
  my $in = STDIN;
  if ($file) {
    open INPUT, $file;
    $in = INPUT;
  }
  my $current_matrix_nb = 0;
  my @matrices = ();
  my $matrix;
  my $command = "";

#  my %prior = ();
  my $l = 0;
  while (my $line = <$in>) {
    $l++;
    next unless ($line =~  /\S/);
    chomp($line);

    ## Read the command line
    if ($line =~ /;assembly # (\d+)\s+seed:\s+(\S+)/) {
      $current_matrix_nb = $1;
      my $seed = $2;
      $matrix = new RSAT::matrix();
      push @matrices, $matrix;
      $matrix->setAlphabet_lc("A","C","G","T");
      $matrix->set_attribute("nrow", 4);
      $matrix->set_attribute("number", $current_matrix_nb);
      $matrix->set_parameter("seed", $seed);
      &RSAT::message::Debug("New matrix from assembly", $current_matrix_nb."/".scalar(@matrices), "seed", $seed) if ($main::verbose >= 4);

    } elsif ($line =~ /^(\S+)\t(\S+)\s+(\S+)\s+isol/) {
      $current_matrix_nb++;
      my $pattern = $1; 
      my $pattern_rc = $2; 
      my $score = $3;
      $matrix = _from_isolated($pattern, $pattern_rc, $score, @matrices);
      push @matrices, $matrix;

    } elsif ($line =~ /^(\S+)\s+(\S+)\s+isol/) {
      $current_matrix_nb++;
      my $pattern = $1; 
      my $score = $2;
      $matrix = _from_isolated($pattern, "", $score, @matrices);
      push @matrices, $matrix;


      ## Consensus from a 2-strand assembly
    } elsif ($line =~ /^(\S+)\s+(\S+)\s+(\S+)\s+best consensus/) {
      $matrix->set_attribute("consensus.assembly", $1);
      $matrix->set_attribute("consensus.assembly.rc", $2);
      $matrix->set_attribute("assembly.top.score", $3);
#      &RSAT::message::Debug("Consensus for matrix", $current_matrix_nb, $1) if ($main::verbose >= 5);

      ## Consensus from a 1-strand assembly
    } elsif ($line =~ /^(\S+)\s+(\S+)\s+best consensus/) {
      $matrix->set_attribute("consensus.assembly", $1);
      $matrix->set_attribute("assembly.top.score", $3);
#      &RSAT::message::Debug("Consensus for matrix", $current_matrix_nb, $1) if ($main::verbose >= 5);

    } elsif ($line =~ /^;/) {
      next;

      ## New site from a two-strand assembly
    } elsif ($line =~ /^(\S+)\t(\S+)\s+(\S+)/) {
      my $pattern = $1; 
      my $pattern_rc = $2;
      my $score = $3;
      my $pattern_id = $pattern."|";
      $pattern =~ s/\./n/g;
      $pattern_rc =~ s/\./n/g;
#      &RSAT::message::Debug("ASSEMBLY LINE", $l, $pattern, $pattern_rc, $score) if ($main::verbose >= 5);
      $matrix->add_site($pattern, , score=>$score, id=>$pattern_id, max_score=>1);

      ## New site from a single-strand assembly
    } elsif ($line =~ /^(\S+)\s+(\S+)/) {
      my $pattern = $1; 
      my $score = $2;
      $pattern =~ s/\./n/g;
#      &RSAT::message::Debug("ASSEMBLY LINE", $l, $pattern, $pattern_rc, $score) if ($main::verbose >= 5);
      $matrix->add_site($pattern, score=>$score, id=>$pattern, max_score=>1);

    } else {
      &RSAT::message::Warning("&RSAT::Matrixreader::_readFromAssemblyFile", "line", $l, "not parsed", $_) if ($main::verbose >= 4);
    }
  }
  close $in if ($file);

  return @matrices;
}

sub _from_isolated {
  my ($pattern, $pattern_rc, $score, @matrices) = @_;
  my $pattern_id = $pattern;
  $pattern =~ s/\./n/g;
  if ($pattern_rc) {
    $pattern_rc =~ s/\./n/g;
    $pattern_id = $pattern."|".$pattern_rc;
  }
  $matrix = new RSAT::matrix();
  $matrix->setAlphabet_lc("A","C","G","T");
  $matrix->set_attribute("nrow", 4);
  $matrix->set_attribute("number", $current_matrix_nb);
  $matrix->set_parameter("seed", $pattern);
  $matrix->set_attribute("consensus.assembly", $pattern);
  $matrix->set_attribute("consensus.assembly.rc", $pattern_rc);
  $matrix->set_attribute("assembly.top.score", $score);
  $matrix->add_site($pattern, score=>$score, id=>$pattern_id,, max_score=>1);
  &RSAT::message::Debug("New matrix from isolated pattern", $current_matrix_nb."/".scalar(@matrices), "seed", $seed) if ($main::verbose >= 4);
  return $matrix;
}

################################################################
=pod

=item _readFromTabFile($file)

Read a matrix from a tab-delimited file. This method is called by the
method C<readFromFile($file, "tab")>.

=cut
sub _readFromTabFile {
    my ($file, %args) = @_;
    &RSAT::message::Info(join("\t", "Reading matrix from tab file\t",$file)) if ($main::verbose >= 3);


    ## open input stream
    my ($in, $dir) = &main::OpenInputFile($file);
    if ($file) {
	open INPUT, $file;
	$in = INPUT;
    }


    ## read header
    if ($args{header}) {
	$header = <$in>;
	$header =~ s/\r//;
	chomp ($header);
	$header =~ s/\s+/\t/g;
	@header = split "\t", $header;
	$matrix->push_attribute("header", @header);
    }


    ## Initialize the matrix list
    my @matrices = ();
    my $matrix = new RSAT::matrix();
    push @matrices, $matrix;
    my $current_matrix_nb = 1;
    my $l = 0;
    while ($line = <$in>) {
      $l++;
      next unless ($line =~ /\S/); ## Skip empty lines
      chomp($line); ## Suppress newline
      $line =~ s/\r//; ## Suppress carriage return
      $line =~ s/\s+/\t/g; ## Replace spaces by tabulation

      next if ($line =~ /^;/) ; # skip comment lines

#	&RSAT::message::Debug("line", $l, $line) if ($main::verbose >= 10);
	## Create a new matrix if required
	if  ($line =~ /\/\//) {
	  $matrix = new RSAT::matrix();
	  push @matrices, $matrix;
	  $current_matrix_nb++;
	  &RSAT::message::Info("line", $l, "new matrix", $current_matrix_number) if ($main::verbose >= 5);
	  next;
	}

      if ($line =~ /^\s*(\S+)\s+/) {

	my @fields = split /\t/, $line;


	## residue associated to the row
	my $residue = lc(shift @fields);

	## skip the | between residue and numbers
	shift @fields unless &main::IsReal($fields[0]);

	$matrix->addIndexedRow($residue, @fields);
      }
    }
    close $in if ($file);


    ## Initialize prior as equiprobable alphabet
    foreach my $matrix (@matrices) {
      my @alphabet = $matrix->getAlphabet();
      my %tmp_prior = ();
      my $prior = 1/scalar(@alphabet);
      foreach my $residue (@alphabet) {
	$tmp_prior{$residue} = $prior;
	#	&RSAT::message::Debug("initial prior", $residue, $prior) if ($main::verbose >= 10);
      }
      $matrix->setPrior(%tmp_prior);
      
      if ($main::verbose >= 3) {
	&RSAT::message::Debug("Read matrix with alphabet", join(":", $matrix->getAlphabet()));
	&RSAT::message::Debug("Initialized prior as equiprobable", join(":", $matrix->getPrior()));
	&RSAT::message::Debug("Matrix size", $matrix->nrow()." rows",  $matrix->ncol()." columns");
      }
    }

    return (@matrices);
}


################################################################
=pod

=item _readFromMEMEFile($file)

Read a matrix from a MEME file. This method is called by the
method C<readFromFile($file, "MEME")>.

=cut
sub _readFromMEMEFile {
  my ($file) = @_;
  &RSAT::message::Info("Reading matrix from consensus file\t", $file) if ($main::verbose >= 3);
    
  ## open input stream
  #    ($in, $dir) = &main::OpenInputFile($file);
  my $in = STDIN;
  if ($file) {
    open INPUT, $file;
    $in = INPUT;
  }
  my @matrices = ();
  my $current_matrix_nb = 0;
  my $matrix;
  my $current_col = 0;
  my $in_proba_matrix = 0;
  my $in_blocks = 0;
  my $width_to_parse = 0;
  my %alphabet = ();
  my %residue_frequencies = ();
  my @alphabet = ();
  my @frequencies = ();
  my $parsed_width = 0;
  my $l = 0;
  while (<$in>) {
    $l++;
    next unless (/\S/);
    s/\r//;
    chomp();
    $_ = &main::trim($_);
    if (/MOTIF\s+(\d+)\s+width =\s+(\d+)\s+sites =\s+(\d+)\s+llr =\s+(\d+)\s+E-value =\s+(\S+)/) {
      &RSAT::message::Debug("line", $l, "Parsing matrix parameters") if ($main::verbose >= 5);

      $current_matrix_nb = $1;
      $width_to_parse = $2;
      $matrix = new RSAT::matrix();
      $matrix->init();
      $matrix->set_attribute("number", $current_matrix_nb);
      $matrix->set_attribute("ncol", $2);
      $matrix->set_parameter("sites", $3);
      $matrix->set_parameter("llr", $4);
      $matrix->set_parameter("E-value", $5);
      $matrix->setPrior(%residue_frequencies);
#      &RSAT::message::Debug("line", $l, "Read letter frequencies", %residue_frequencies) if ($main::verbose >= 10);
      $matrix->setAlphabet_lc(@alphabet);
      $matrix->force_attribute("nrow", scalar(@alphabet)); ## Specify the number of rows of the matrix
      push @matrices, $matrix;

      ## Parse alphabet
    } elsif (/Letter frequencies in dataset/) {
      my $alphabet = <$in>;
      $alphabet = lc($alphabet);
      $alphabet = &main::trim($alphabet);
      %residue_frequencies = split /\s+/, $alphabet;
      @alphabet = sort (keys %residue_frequencies);
#      &RSAT::message::Debug("line", $l, "Read letter frequencies", %residue_frequencies) if ($main::verbose >= 10);

      ## Index the alphabet
      foreach my $l (0..$#alphabet) {
	$alphabet{$alphabet[$l]} = $l;
      }

      ## Parse BLOCKS format
    } elsif (/Motif (\d+) in BLOCKS format/) {
      $current_matrix_nb = $1;
      $in_blocks = 1;
      &RSAT::message::Debug("line", $l, "Starting to parse BLOCKS format") if ($main::verbose >= 5);

    } elsif ($in_blocks) {
      if (/(\S+)\s+\(\s*\d+\)\s+(\S+)/) {
	my $seq_id = $1;
	my $seq = lc($2);
	my $seq_len =  length($seq);
	if ($seq_len > 0) {
	  $parsed_width = &main::max($parsed_width, $seq_len);
	  $matrix->add_site($seq, score=>1, id=>$seq_id, max_score=>0);
	}

      } elsif (/\/\//) {
	&RSAT::message::Debug("line", $l, "BLOCKS format parsed") if ($main::verbose >= 5);
	$in_blocks = 0;

      }
    }
  }
  close $in if ($file);
  return @matrices;
#  return $matrices[0];
}

################################################################
=pod

=item _readFromMotifSamplerFile($file)

Read a matrix from a MotifSampler file. This method is called by the
method C<readFromFile($file, "MotifSampler")>.

TO BE IMPLEMENTED

=cut

sub _readFromMotifSamplerFile {
    &RSAT::error::FatalError("The MotifSampler format is not yet supported in this version of the program.");

}




################################################################
=pod

=item _readFromClustalFile($file)

Read a matrix from a multiple alignment in clustal format (extension
 .aln).  This method is called by the method C<readFromFile($file,
 "clustal")>.

=cut
sub _readFromClustalFile {
    my ($file) = @_;

    my @matrices = ();
    my $matrix = new RSAT::matrix();
    push @matrices, $matrix;
    
    ## open input stream
    my $in = STDIN;
    if ($file) {
	open INPUT, $file;
	$in = INPUT;
    }

    ## Check the header
    my $header = <$in>;
    unless ($header =~ /clustal/i) {
	&main::Warning("This file does not contain the clustal header");
    }

    ## Read the sequences
    my %sequences = ();
    warn "; Reading sequences\n" if ($main::verbose >= 3);
    while (<$in>) {
	next unless (/\S/);
	s/\r//;
	chomp();
	if (/^\s*(\S+)\s+(.+)$/) {
	    my $seq_id = $1;
	    next if ($seq_id eq "*"); ## asterisks are used to mark conservation
	    my $new_seq = $2;
	    
	    ## index the new sequence
	    $sequences{$seq_id} .= $new_seq;
	    warn join ("\t", ";", "Sequence", $seq_id, 
		       length($new_seq), length($sequences{$seq_id}),
		       ),"\n" if ($main::verbose >= 5);
	}
    }
    
    ## Calculate count matrix
    my %matrix = ();
    my %prior = ();
    my $ncol = 0;
    my $nrow = 0;
    &RSAT::message::Info("Calculating profile matrix from sequences") if ($main::verbose >= 3);
    foreach my $seq_id (sort keys %sequences) {
	my $sequence = $sequences{$seq_id};
	$sequence =~ s/\s+//g;

	################################################################
	## Distinguish between insertions and leading/trailing gaps
	$terminal_gap_char = ".";

	## Substitute leading gaps
	if ($sequence =~ /^(\-+)/) {
	    $leading_gap_len = length($1);
	    my $leading_gap = ${terminal_gap_char}x$leading_gap_len;
	    $sequence =~ s|^(\-+)|${leading_gap}|;
	}
	## Substitute trailing gaps
	if ($sequence =~ /(\-+)$/) {
	    $trailing_gap_len = length($1);
	    my $trailing_gap = ${terminal_gap_char}x$trailing_gap_len;
	    $sequence =~ s|(\-+)$|${trailing_gap}|;
	}
	warn join ("\t",";", $seq_id,$sequence), "\n" if ($main::verbose >= 5);
	    
	$ncol = &main::max($ncol, length($sequence));
	warn join ("\t", ";", "Sequence", $seq_id, length($sequence)),"\n" if ($main::verbose >= 5);
	my @sequence = split '|', $sequence;
	foreach my $i (0..$#sequence) {
	    my $res = lc($sequence[$i]);
	    next if ($res eq "N"); ## BEWARE: THIS IS FOR DNA ONLY
#	    next if ($res eq "-");
	    next if ($res eq "."); ## leading and trailing gaps
	    next if ($res eq "*");
	    $prior{$res}++;
	    $matrix{$res}->[$i] += 1;
	}
    }
    $matrix->set_attribute("ncol", $ncol);

    ## Define prior probabilities, alphabet, and matrix size
    my @alphabet = sort keys %prior;
    my $alpha_sum = 0;
    foreach my $res (@alphabet) {
	$alpha_sum += $prior{$res};
    }
    foreach my $res (@alphabet) {
	if ($alpha_sum > 0) {
	    $prior{$res} /= $alpha_sum;
	} else {
	    $prior{$res} = 0;
	}
#	warn join "\t", $res, $alpha_sum, $prior{$res};
    }
    $matrix->setPrior(%prior);

    ## Store the matrix
    my @matrix = ();
    foreach my $r (0..$#alphabet) {
	my $res = $alphabet[$r];
	my @row = @{$matrix{$res}};
	$nrow++;
	foreach $i (0..($ncol-1)) {
	    $row[$i] = 0 unless (defined($row[$i]));
	}
	$matrix->addRow(@row);
	warn join ("\t", "Adding row", $r, $res, join ":", @row, "\n"), "\n" if ($main::verbose >= 4); 
    }
    $matrix->setAlphabet_lc(@alphabet);
    $matrix->force_attribute("ncol", $ncol);
    $matrix->force_attribute("nrow", $nrow);

    warn join ("\t", "; Matrix size",  
	       $nrow,
	       $ncol,
	       $matrix->nrow(), 
	       $matrix->ncol()), "\n" 
		  if ($main::verbose >= 3);
    close $in if ($file);

    return (@matrices);
}



return 1;

__END__

