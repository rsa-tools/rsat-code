###############################################################
#
# Class MarkovModel
#
package RSAT::MarkovModel;

## TO DO : 
## - fix bugs with the method average_strands()
## - for the model update, take into account the strand-insensitive models

use RSAT::GenericObject;
use RSAT::SeqUtil;
use RSAT::error;
use Data::Dumper;

@ISA = qw( RSAT::GenericObject );

### class attributes
%supported_input_formats = ("oligo-analysis"=>1, 
			    "oligos"=>1, ## Abbreviation for oligo-analysis
			    "dyads"=>1, ## Abbreviation for dyad-analysis
			    "motifsampler"=>1, 
			    "ms"=>1, ## Abbreviation for MotifSampler
			    "meme"=>1,
			   );
%supported_output_formats = ("tab"=>1, ## transition_freq
			     "transitions"=>1, ## transition table
			     "tables"=>1, ## frequency + transition tables
			     "patser"=>1,
			     "oligo-analysis"=>1,
			     "oligos"=>1, ## Abbreviation for oligo-analysis
			     "meme"=>1,
			     "motifsampler"=>1,
			     "ms"=>1, ## Abbreviation for MotifSampler
			    );
@supported_input_formats = sort keys %supported_input_formats;
@supported_output_formats = sort keys %supported_output_formats;
$supported_input_formats = join(",", @supported_input_formats);
$supported_output_formats = join(",", @supported_output_formats);


=pod

=head1 NAME

    RSAT::MarkovModel

=head1 DESCRIPTION

Class for handling a Markov model. 

=cut


################################################################
=pod

=item B<new()>

Create a new MarkovModel.

=cut
sub new {
    my ($class, %args) = @_;
    my $object = bless {
	}, $class;
    return $object;
}



################################################################
=pod

=item B<get_supported_input_formats()>

Return input formats supported by the function load_from_file().

=cut
sub get_supported_input_formats { 
    my $class = ref($_[0]) || $_[0];
    return %{$class."::supported_input_formats"}; 
}

################################################################
=pod

=item B<get_supported_output_formats()>

Return output formats supported by the function to_string().

=cut
sub get_supported_output_formats { 
    my $class = ref($_[0]) || $_[0];
    return %{$class."::supported_output_formats"}; 
}




################################################################
=pod

=item B<load_from_file($bg_file, $format)>

Load the markov model from a file. 

In order to ensure case-insensivity, Markov models are automatically
converted to lowercases.

=cut
sub load_from_file {
  my ($self, $bg_file, $format) = @_;

  $format =~ s/oligo-analysis/oligos/;
  $format =~ s/dyad-analysis/dyads/;
  $format =~ s/^ms$/motifsampler/;

  $self->force_attribute("strand", "undef");

  if ($format eq "oligos") {
    $self->load_from_file_oligos($bg_file);
  } elsif ($format eq "dyads") {
    $self->load_from_file_oligos($bg_file, 1); ## Dyads with spacing 0 are converted into oligos
  } elsif ($format eq "meme") {
    $self->load_from_file_meme($bg_file);
  } elsif ($format eq "motifsampler") {
    $self->load_from_file_MotifSampler($bg_file);
  } elsif ($format eq "") {
    &RSAT::error::FatalError(join("\t", "MarkovModel::load_from_file", 
				  "the format of the backgroun model must be specified"));
  } else {

    &RSAT::error::FatalError(join("\t", "MarkovModel::load_from_file", 
				  $format, "invalid format. Supported:",
				  join(",", @supported_input_formats)));
  }
  $self->normalize_transition_frequencies();
}



################################################################
=pod

=item B<load_from_file_oligos($bg_file)>

Load the markov model from a result file from oligo-analysis. 

The model file is a tab-delimited text file indicating the expected
frequency of each word of length m+1, for a model of order m. These
files can be produced with the ccmmand oligo-analysis. When new
genomes are installed, background files are automatically calculated
for upstream sequences, and stored in
$RSAT/data/genomes/[Organism_name]/oligo-frequencies.

In order to ensure case-insensivity, Markov models are automatically
converted to lowercases.

=cut
sub load_from_file_oligos {
  my ($self, $bg_file, $dyad_conversion) = @_;

  &RSAT::message::TimeWarn("Loading Markov model from oligo file $bg_file") if ($main::verbose >= 2);

  my ($file_type, %patterns) = &main::ReadPatternFrequencies($bg_file) ;

  ################################################################
  ## For ydad input format, convert dyads into oligos: only retain
  ## dyads with spacing 0
  if ($dyad_conversion) {
    my %oligo_patterns;
    foreach my $dyad_seq (keys %patterns) {
      ## Dyads with spacing 0 are converted into oligos
      if ($dyad_seq =~ /n\{0\}/) {
	my $oligo_seq = $dyad_seq;
	$oligo_seq =~ s/n\{0\}//;
	$patterns{$oligo_seq} = $patterns{$dyad_seq};
      }
      delete $patterns{$dyad_seq};
    }
  }
  ## This is a bit tricky: ReadExpectedFrequencies sets a global variable
  ## $file_type to "2str" if the model is strand-insensitive/
  if ($file_type eq "2str") {
    $self->force_attribute("strand", "insensitive");
  } elsif ($file_type eq "1str") {
    $self->force_attribute("strand", "sensitive");
  } else {
    $self->force_attribute("strand", "undef");
  }

#   if ($main::verbose >= 10) {
#     foreach my $oligo_seq (sort(keys(%patterns))) {
#       &RSAT::message::Debug("Loaded oligo frequency", $oligo_seq,  $patterns{$oligo_seq}->{exp_freq});
#     }
#   }

  $self->oligos_to_frequency_table(%patterns);
  #    &RSAT::message::Debug("MARKOV MODEL", $order, join (' ', @patterns)) if ($main::verbose >= 5);
}

################################################################

=pod

=item B<oligos_to_frequency_table>

Convert oligomer frequencies into a prefix/suffix table (one row per
prefix, one column per suffix, each cell indicates the absolute
frequency of the prefix->suffix transition).

=cut
sub oligos_to_frequency_table {
  my ($self, %patterns) = @_;
  my @patterns = sort(keys(%patterns));
  my $pattern_len = length($patterns[0]);
  my $pattern_nb = scalar(keys(%patterns));
  my $order =  $pattern_len -1 ;	## Use the first pattern to calculate model order
  $self->force_attribute("order", $order);
  &RSAT::message::Info("Converting oligos into frequency tables", $pattern_nb." patterns", "length", $pattern_len, "order", $order) 
    if ($main::verbose >= 3);

  foreach my $pattern_seq (@patterns) {
    my $pattern_freq =  $patterns{$pattern_seq}->{exp_freq};
    ## Check pattern length
    $pattern_seq = lc($pattern_seq); ## Ensure case-insensitivity
    $pattern_len = length($pattern_seq);
    #    &RSAT::message::Debug($pattern_seq, $pattern_len, $pattern_freq, 
    #    			      join(";", keys %{$patterns{$pattern_seq}})) if ($main::verbose >= 0);
    &RSAT::error::FatalError("All patterns should have the same length in a Markov model file.") 
      unless $pattern_len = $order+1;
    my $prefix = substr($pattern_seq,0,$order);
    my $suffix = substr($pattern_seq,$order, 1);
    $self->{oligo_freq}->{$prefix}->{$suffix} = $pattern_freq;
    #	&RSAT::message::Debug("transition count", $prefix.".".$suffix, $pattern_freq, $self->{oligo_freq}->{$prefix}->{$suffix}) if ($main::verbose >= 10);
  }
}

################################################################
=pod 

=item B<transitions_to_oligo_frequencies>


Convert a prefix/suffix table into oligo frequencies.

The frequency table has one row per prefix, one column per suffix,
each cell indicates the absolute frequency of the prefix->suffix
transition. Oligos are obtained by concatenating each prefix with each suffix

=cut
sub transitions_to_oligo_frequencies {
  my ($self) = @_;
  &RSAT::message::Info("Converting transitions into oligo frequencies") if ($main::verbose >= 3);
  my %patterns = ();
  my @prefix = sort($self->get_attribute("prefixes"));
  my @suffix = sort($self->get_attribute("suffixes"));
  foreach my $prefix (sort (@prefix)) {
    foreach my $suffix (sort(@suffix)) {
      my $oligo_seq = $prefix.$suffix;
      my $oligo_freq = $self->{prefix_proba}->{$prefix}*$self->{transitions}->{$prefix}->{$suffix};
      $patterns{$oligo_seq}->{exp_freq} = $oligo_freq;
#      &RSAT::message::Debug($prefix, $suffix, $oligo_seq, $oligo_freq) if ($main::verbose >= 10);
    }
  }
  return %patterns;
}


################################################################
=pod

=item B<load_from_file_meme($bg_file)>

Load the Markov model from a MEME background file.

MEME (http://meme.sdsc.edu/meme/) is a matrix-based pattern discovery
program developed by Tim Bailey, and supporting Marckov chain
background models. 

The MEME suite includes

 - I<fasta-get-markov>, which esimtates a background model
   from a sequence file.

 - I<meme>, for discovering new motifs from a sequence set.

 - I<mast>, for localizing known motifs in a sequence set.

 - a set of other utilities described in MEME documentation.

In order to ensure case-insensivity, Markov models are automatically
converted to lowercases.

=cut
sub load_from_file_meme {
  my ($self, $bg_file) = @_;

  &RSAT::message::TimeWarn("Loading Markov model from MEME background file $bg_file") if ($main::verbose >= 2);

  #    my ($file_type, %patterns) = &main::ReadPatternFrequencies($bg_file) ;

  my ($in, $dir) = &main::OpenInputFile($bg_file);
  my $max_order = 0;
  my $order  =0;
  my %freq = ();
  while (<$in>) {
    if (/^# order (\d+)/) {
      $order = $1;
      if ($order > $max_order) {
	## Meme bg files contain frequencies for all oligo lengths
	## from 0 to k=m+1 (where m the order of the Markov
	## chain). We only need the largest patterns, so we
	## reinitialize the frequency table after smaller words.
	%freq = ();
      }
      $max_order = &RSAT::stats::max($max_order, $order);
    } elsif (/^(\S+)\s+(\S+)/) {
      my $pattern_seq = $1;
      my $freq = $2;
      $freq{$pattern_seq} = $freq;
    }
  }
  close $in if ($bg_file);

  $order = $max_order;
  $self->force_attribute("order", $order);
  &RSAT::message::Info("Loading Markov model of order", $order) if ($main::verbose >= 3);

  ## Calculate alphabet from expected frequency keys
  foreach my $pattern_seq (keys %freq) {
    ## Check pattern length
    my $pattern_freq =  $freq{$pattern_seq};
    $pattern_seq = lc($pattern_seq);
    my $pattern_len = length($pattern_seq);
    #	&RSAT::message::Debug($pattern_seq, $pattern_len, $pattern_freq, 
    #			      join(";", keys %{$patterns{$pattern_seq}})) if ($main::verbose >= 0);
    &RSAT::error::FatalError("All patterns should have the same length in a Markov model file.") 
      unless $pattern_len = $order+1;
    my $prefix = substr($pattern_seq,0,$order);
    my $suffix = substr($pattern_seq,$order, 1);
    $self->{oligo_freq}->{$prefix}->{$suffix} = $pattern_freq;
    #	&RSAT::message::Debug("transition count", $prefix.".".$suffix, $pattern_freq, $self->{oligo_freq}->{$prefix}->{$suffix}) if ($main::verbose >= 10);

  }
  #    &RSAT::message::Debug("MARKOV MODEL", $order, join (' ', @patterns)) if ($main::verbose >= 5);
}


################################################################
=pod

=item B<load_from_file_MotifSampler($bg_file)>

Load the markov model from a result file from MotifSampler.

These files are used by the various programs of the software suite
INCLUSive (http://homes.esat.kuleuven.be/~thijs/download.html),
developed by Gert Thijs.

This suite includes the program CreateBackgroundModel, which creates a
Markov model from a background sequence.

=cut

sub load_from_file_MotifSampler {
  my ($self, $bg_file) = @_;

  &RSAT::message::TimeWarn("Loading Markov model from MotifSampler file $bg_file") if ($main::verbose >= 2);

  my @prefixes = ();
  my @suffixes = qw(a c g t);
  my %prefix_proba = ();
  my %suffix_proba = ();

  my ($in, $dir) = &main::OpenInputFile($bg_file);
  while (<$in>) {
    next unless /\S/;		## Skip empty lines

    if (/^#Order\s*=\s*(\d+)/i) {
      ## Markov order
      my $order = $1;
      $self->force_attribute("order", $order);
      if ($order == 0) {
	@prefixes = ("");
      } else {
	@prefixes = &RSAT::SeqUtil::all_possible_oligos($order);
      }
      &RSAT::message::Info("Markov order", $order, scalar(@prefixes)." prefixes")
	if ($main::verbose >= 2);

    } elsif (/^#snf/i) {
      ## Single nucleotide frequencies
      &RSAT::message::Info("Reading single nucleotide probabilities") if ($main::verbose >= 3);
      my $line = (<$in>);
      chomp($line);
      my @probas = split (/\s+/, $line);
      foreach my $i (0..$#suffixes) {
	my $suffix= $suffixes[$i];
	my $proba= $probas[$i];
	$suffix_proba{$suffix} = $proba;
	&RSAT::message::Debug("letter proba", $suffix, $proba) if ($main::verbose >= 5);
      }
      $self->set_hash_attribute("suffix_count", %suffix_proba);
      $self->set_hash_attribute("suffix_proba", %suffix_proba);

    } elsif (/^#oligo frequency/i) {
      ## Prefix probabilities
      &RSAT::message::Info("Reading prefix probabilities") if ($main::verbose >= 3);
      foreach my $prefix (@prefixes) {
	my $proba = (<$in>);
	chomp($proba);
	$proba = &RSAT::util::trim($proba); ## Remove leading and trailing spaces
	$prefix_proba{$prefix} = $proba;
	&RSAT::message::Debug("prefix proba", $prefix, $proba) if ($main::verbose >= 5);
      }
      $self->set_hash_attribute("prefix_count", %prefix_proba);
      $self->set_hash_attribute("prefix_proba", %prefix_proba);


    } elsif (/^#transition matrix/i) {
      ## transition probabilities
#      &RSAT::message::Debug("Reading transition matrix") if ($main::verbose >= 5);
      foreach my $prefix (@prefixes) {
	my $line = (<$in>);
	chomp($line);
	@transition_proba = split(/\s+/, $line);
	foreach my $i (0..$#suffixes) {
	  my $suffix= $suffixes[$i];
	  my $proba= $transition_proba[$i];
	  $self->{oligo_freq}->{$prefix}->{$suffix} = $proba*$prefix_proba{$prefix};
	  $self->{transitions}->{$prefix}->{$suffix} = $proba*$prefix_proba{$prefix};
	  &RSAT::message::Debug("transition", $prefix, $prefix_proba{$prefix}, $suffix, $suffix_proba{$suffix}, $proba, ) if ($main::verbose >= 5);
	}
      }
    }
  }
  close $in if ($bg_file);
}



############################################################
=pod

=item B<add_pseudo_freq>

Add pseudo-frequencies to transition frequencies. This ensures that
each possible transition will have a non-null probability, even though
this transition might be absent from the sequence on which the Markov
model was trained.

=cut
sub add_pseudo_freq {
  my ($self) = @_;

  my $pseudo_freq = $self->get_attribute("bg_pseudo") || 0;

  ## get alphabet
  my @dna_alphabet =  qw (a c g t);

  ################################################################
  ## Ensure that all possible suffixes are taken into consideration
  my @possible_suffixes = @dna_alphabet;
  $self->set_array_attribute("suffixes", @possible_suffixes);
  my $alpha_size = scalar(@possible_suffixes);

  ################################################################
  ## Calculate all possible prefix oligomers, even those
  ## which are not observed in the input background model
  my @possible_prefixes;
  if ($self->{order} > 0) {
    @possible_prefixes  = &RSAT::SeqUtil::all_possible_oligos($self->{order}, @dna_alphabet);
  } elsif ($self->{order} == 0) {
    @possible_prefixes  = ("");
  }
  $self->set_array_attribute("prefixes", @possible_prefixes);

  ################################################################
  ## Update prefix sum
  my $total_prefix_sum = 0;
  foreach my $prefix (@possible_prefixes) {
    $total_prefix_sum += $self->{prefix_sum}->{$prefix};
  }
#  &RSAT::message::Debug("Total prefix sum", $total_prefix_sum) if ($main::verbose >= 0);

  my $prefix_pseudo_count = $total_prefix_sum*$pseudo_freq/scalar(@possible_prefixes);
  foreach my $prefix (@possible_prefixes) {
    $self->{prefix_sum}->{$prefix} = $self->{prefix_sum}->{$prefix}*(1-$pseudo_freq) + $prefix_pseudo_count;
    #    &RSAT::message::Debug("Corrected prefix sum", $prefix, $self->{prefix_sum}->{$prefix}) if ($main::verbose >= 10);
  }

  ################################################################
  ## Update suffix sum
  my $total_suffix_sum = 0;
  foreach my $suffix (@possible_suffixes) {
    $total_suffix_sum += $self->{suffix_sum}->{$suffix};
  }
#  &RSAT::message::Debug("Total suffix sum", $total_suffix_sum) if ($main::verbose >= 0);
  my $suffix_pseudo_count = $total_suffix_sum*$pseudo_freq/scalar(@possible_suffixes);
  foreach my $suffix (@possible_suffixes) {
    $self->{suffix_sum}->{$suffix} = $self->{suffix_sum}->{$suffix}*(1-$pseudo_freq) + $suffix_pseudo_count;
    #    &RSAT::message::Debug("Corrected suffix sum", $suffix, $self->{suffix_sum}->{$suffix}) if ($main::verbose >= 10);
  }


  ################################################################
  ## Update transition frequencies
  foreach my $prefix (@possible_prefixes) {

    ## Compute the sum of transition counts for the current prefix
#    my $current_prefix_sum = 0;
#    foreach my $suffix (@possible_suffixes) {
#      $current_prefix_sum += $self->{oligo_freq}->{$prefix}->{$suffix}
#	if (defined($self->{oligo_freq}->{$prefix}->{$suffix}));
#    }
    my $current_prefix_sum = $self->{prefix_sum}->{$prefix};

    ## pseudo count varies depending on the prefix sum
    my $pseudo_count = $pseudo_freq*$current_prefix_sum;
    if ($current_prefix_sum == 0) {
      $pseudo_count = $pseudo_freq;
    }
    &RSAT::message::Debug("Prefix", $prefix, $current_prefix_sum, "Pseudo count=".$pseudo_count) if ($main::verbose >= 0);

    foreach my $suffix (@possible_suffixes) {
      if (defined($self->{oligo_freq}->{$prefix}->{$suffix})) {
	## Adding the pseudo-freq on the background model.
	my $oligo_freq_corrected = 
	  (1-$pseudo_freq)*$self->{oligo_freq}->{$prefix}->{$suffix} + $pseudo_freq/$alpha_size;
	$self->{oligo_freq}->{$prefix}->{$suffix} = $oligo_freq_corrected;
      } else {			## missing transitions
	$self->{oligo_freq}->{$prefix}->{$suffix} = $pseudo_freq/$alpha_size;
      }
      &RSAT::message::Debug("Corrected transition count", $prefix.".".$suffix, 
			    "pseudo_freq=".$pseudo_freq, 
			    "prefix_sum=".$current_prefix_sum,
			    "pseudo_count=".$pseudo_count,
			    $self->{oligo_freq}->{$prefix}->{$suffix}) if ($main::verbose >= 0);

    }
  }

}

################################################################
=pod

=item B<normalize_transition_frequencies>

Normalize transition frequencies, i.e. make sure that, for each prefix, the sum
of transition frequencies (frequency of all possible suffixes given the
prefix) is 1.

=cut
sub normalize_transition_frequencies {
  my ($self, %args) = @_;

  &RSAT::message::TimeWarn(join("\t", "MarkovModel", "Normalizing transition frequencies")) if ($main::verbose >= 3);

  ################################################################
  ## Calculate sum of counts (frequencies) per prefix and suffix by
  ## summing over all transition frequencies
  &RSAT::message::TimeWarn("Computing prefix and suffix sums") if ($main::verbose >= 4);
  my %prefix_sum = ();
  my %suffix_sum = ();
  my $freq_sum = 0;
  my $p=0;			## Number of prefixes
  my $s=0;			## Number of suffixes
  my $strand = $self->get_attribute("strand") || "undef";

  foreach my $prefix (sort keys(%{$self->{oligo_freq}})) {
    $p++;
    #	foreach my $suffix (keys (%suffix_sum)) {
    foreach my $suffix (sort keys (%{$self->{oligo_freq}->{$prefix}})) {
      $s++;
      my $pattern_count = $self->{oligo_freq}->{$prefix}->{$suffix};

      ## Particular treatment for case-insensitive models: the
      ## oligo_freq of each non-rc-palindromic pattern is tken into
      ## account in two entries of the hash table (direct + rc
      ## pattern) -> to compute the sum of frequencies, on has to count it once only
      my $to_sum = 1;
      if ($strand eq "insensitive") {
	my $oligo_seq = $prefix.$suffix;
	my $oligo_rc = lc(&main::SmartRC($oligo_seq));
	if ($oligo_rc lt $oligo_seq) {
	  $to_sum = 0;
	}
      }

      ## Increment frequencies
      if ($to_sum) {
	$freq_sum += $pattern_count;
	$suffix_sum{$suffix} += $pattern_count;
	$prefix_sum{$prefix} += $pattern_count;
      }

#      &RSAT::message::Debug("Count sum", $p, $prefix.".".$suffix, $pattern_count, $to_sum, $freq_sum) if ($main::verbose >= 0);
    }
  }

  $self->set_hash_attribute("prefix_sum", %prefix_sum);
  $self->set_hash_attribute("suffix_sum", %suffix_sum);

  ## Store prefixes and suffixes in arrays for quick access
  $self->set_array_attribute("prefixes", sort(keys(%prefix_sum)));
  $self->set_array_attribute("suffixes", sort(keys(%suffix_sum)));

  ## Add pseudo-counts to transition counts. This has to be done after the computation of prefix counts
  unless ($args{no_pseudo}) {
    $self->add_pseudo_freq();	### adding the pseudo-freq
  }

  ## Retrieved the updated prefix and suffix sums
  %prefix_sum = $self->get_attribute("prefix_sum");
  %suffix_sum = $self->get_attribute("suffix_sum");

  ## Calculate transition frequencies == Transition table with
  ## conditional probabilities
  my $missing_transitions = 0;
  foreach my $prefix ($self->get_attribute("prefixes")) {
    foreach my $suffix ($self->get_attribute("suffixes")) {
      #	foreach my $suffix (keys (%suffix_sum)) {
      if (defined($self->{oligo_freq}->{$prefix}->{$suffix})) {
	if ($prefix_sum{$prefix} > 0) {
	  $self->{transitions}->{$prefix}->{$suffix} = 
	    $self->{oligo_freq}->{$prefix}->{$suffix}/$prefix_sum{$prefix};
	} else {
	  $self->{transitions}->{$prefix}->{$suffix} =  0;
	}
      } else {
	$missing_transitions++;
	$self->{transitions}->{$prefix}->{$suffix} = 0;
	#		&RSAT::message::Warning(join(" ",
	#					     "No transition between prefix",$prefix, 
	#					     "and suffix", $suffix)) if ($main::verbose >= 3);
      }
      ;
    }
  }
  $self->force_attribute("missing_transitions", $missing_transitions);


  ## Calculate prefix probabilities
  my %prefix_proba = ();
  foreach my $prefix ($self->get_attribute("prefixes")) {
    #    foreach my $prefix (keys(%prefix_sum)) {
    $prefix_proba{$prefix} = $prefix_sum{$prefix}/$freq_sum;
  }

  ## Calculate suffix probabilities
  my %suffix_proba = ();
  foreach my $suffix ($self->get_attribute("suffixes")) {
    $suffix_proba{$suffix} = $suffix_sum{$suffix}/$freq_sum;
  }

  ## store prefix sums and suffix sums
  $self->set_hash_attribute("prefix_proba", %prefix_proba);
  $self->set_hash_attribute("suffix_proba", %suffix_proba);

  ## Calculate relative frequencies
  foreach my $prefix ($self->get_attribute("prefixes")) {
    foreach my $suffix ($self->get_attribute("suffixes")) {
      $self->{oligo_freq_rel}->{$prefix}->{$suffix} = $self->{oligo_freq}->{$prefix}->{$suffix}/$freq_sum;
    }
  }

  #    &RSAT::message::Debug("SUFFIX PROBA", $self->set_hash_attribute("suffix_proba", %suffix_proba)) if ($main::verbose >= 5);
  ## Average counts and frequencies on both strands if required
  #       my $strand = $self->get_attribute("strand") || "sensitive";
  #     if ($strand eq "insensitive") {
  # 	$self->average_strands();
  #     }

  &RSAT::message::TimeWarn(join("\t", 
				"Normalized background model", 
				"prefixes: ".$p,
				"transitions: ".$s)) if ($main::verbose >= 3);
}

################################################################
=pod

=item B<check_missing_transitions>

Check that the transition matrix contains values >0 in all cells. If
not, report the number of missing transitions.

=cut 
sub check_missing_transitions {
    my ($self) = @_;

    my $missing_transitions = 0;
    foreach my $prefix ($self->get_attribute("prefixes")) {
	foreach my $suffix ($self->get_attribute("suffixes")) {
	    unless ((defined($self->{transitions}->{$prefix}->{$suffix})) 
		    && ($self->{transitions}->{$prefix}->{$suffix} > 0)) {
		$missing_transitions++;
		&RSAT::message::Warning(join(" ",
					     "No transition between prefix",$prefix, 
					     "and suffix", $suffix)) if ($main::verbose >= 3);
	    };
	}
    }
    $self->force_attribute("missing_transitions", $missing_transitions);
    if ($missing_transitions > 0) {
	&RSAT::message::Warning(join(" ", $missing_transitions,
				     "missing transitions in the transition matrix.",
				     "Over-fitting risk.You should better sequences or a lower order Markov model. ")) if
					 ($main::verbose >= 2);
    }
}

################################################################
=pod

=item B<check_transition_alphabet>

Check the transition matrix to suppress residues which are not part of
the accepted alphabet. Such residues can be present in the sequence
for various reasons. For example, the sequence can contain stretches
of N corresponding to masked repeats.

=cut
sub check_transition_alphabet {
    my ($self) = @_;

    my $seq_type = $self->get_attribute("seq_type") || "dna";
    my %accepted_residues = &RSAT::SeqUtil::get_accepted_residues($seq_type);
    &RSAT::message::Info("Accepted residues", $seq_type, %accepted_residues) if ($main::verbose >= 3);

    my $accepted_expression = join ('', '^[', sort(keys(%accepted_residues)), ']*$');

    my %suppressed_prefixes = ();
    my %suppressed_suffixes = ();
    my %checked_prefixes = ();
    my %checked_suffixes = ();
    foreach my $prefix (keys(%{$self->{oligo_freq}})) {
#    foreach my $prefix ($self->get_attribute("prefixes")) {
      if ($prefix =~ /${accepted_expression}/) {
	$checked_prefixes{$prefix} = 1;
      } else {
	$suppressed_prefixes{$prefix}++;
#	delete $self->{transitions}->{$prefix};
	delete $self->{oligo_freq}->{$prefix};
	&RSAT::message::Info(join(" ",
				     "Supressing prefix",$prefix, 
				     "from transition matrix")) if ($main::verbose >= 5);
	next;
      }
      foreach my $suffix (keys(%{$self->{oligo_freq}->{$prefix}})) {
#      foreach my $suffix ($self->get_attribute("suffixes")) {
	$suppressed_suffixes{$suffix}++;
	if ($suffix =~ /${accepted_expression}/) {
	  $checked_suffixes{$suffix} = 1;
	} else {
	  delete $self->{transitions}->{$prefix}->{$suffix};
	  delete $self->{oligo_freq}->{$prefix}->{$suffix};
	  &RSAT::message::Info(join(" ",
				    "Supressing transition from",$prefix, "to", $suffix, 
				    "from transition matrix")) if ($main::verbose >= 5);
	}
      }
    }

    if (scalar(keys(%suppressed_prefixes)) > 0) {
      $self->set_array_attribute("prefixes", sort(keys(%checked_prefixes)));
      &RSAT::message::Info("Suppressed", scalar(keys(%suppressed_prefixes)),"prefixes, incompatible with sequence type", $seq_type)
	if ($main::verbose >= 3);
      &RSAT::message::Debug("Checked prefixes", join (",", $self->get_attribute("prefixes")),
			   ) if ($main::verbose >= 5);
    }
    if (scalar(keys(%suppressed_suffixes)) > 0) {
      $self->set_array_attribute("suffixes", sort(keys(%checked_suffixes)));
      &RSAT::message::Info("Suppressed", scalar(keys(%suppressed_suffixes)),"suffixes, incompatible with sequence type", $seq_type)
	if ($main::verbose >= 3);
    }
}


################################################################
=pod

=item B<calc_from_seq($sequence, [add=>0|1])>

Calculate background model from a sequence.

If the argument add=>1 is specified, the new sequence is added to the
background model.

=cut
sub calc_from_seq {
  my ($self, $sequence, %args) = @_;
  my $seq_len = length($sequence);
  my $order = $self->get_attribute("order");
  $sequence = lc($sequence);
  my $previous_words = $self->get_attribute("training_words") || 0;

  if ($args{add}) {
    &RSAT::message::TimeWarn(join(" ", 
				  "Updating markov model (order ".$order.")",
				  "by adding sequence of length", $seq_len)) if ($main::verbose >= 2);

  } else {
    $self->set_hash_attribute("oligo_counts",()); ## contains the words occurences, not frequencies
    $self->force_attributes("oligo_counts_defined", 1);
    if ($main::verbose >= 2) {
      &RSAT::message::TimeWarn(join(" ", 
				    "Calculating markov model (order ".$order.")",
				    "from sequence of length", $seq_len));
      &RSAT::message::Warning(join("", "RSAT::MarkovModel->calc_from_seq: sequence is too short (len=",$seq_len,
				   "). It must be larger than Markov order + 1 = ", $order+1))
	unless ($seq_len > $order + 1);
    }
  }

  ## Transition counts : count the occurences of the words on 1 strand
  my $last_pos = $seq_len - $order;
  for my $offset (0..($last_pos-1)) {
    my $prefix = substr($sequence, $offset, $order);
    my $suffix = substr($sequence, $offset + $order,1);
    $self->{oligo_counts}->{$prefix}->{$suffix}++;
  }
  $self->force_attribute("training_words", $last_pos + $previous_words);

  ## Delete transitions between letters which do not belong to the accepted alphabet
  $self->check_transition_alphabet();

  ## Convert occurences (word count) into relative frequencies in order to add the pseudo-frequencies
  ## as done when the bg model is read from a bg file, where relative frequencies (not word count) are
  ## read from the file and used to calculate transition matrix.
  foreach my $prefix (sort keys(%{$self->{oligo_counts}})) {
    foreach my $suffix (sort keys(%{$self->{oligo_counts}->{$prefix}})) {
      my $pattern_count = $self->{oligo_counts}->{$prefix}->{$suffix};
      my $pattern_rel_freq = $pattern_count / $self->{training_words};
      $self->{oligo_freq}->{$prefix}->{$suffix} = $pattern_rel_freq;
    }
  }

  ## Initialize transition frequencies
  $self->set_hash_attribute("transition_freq", %transition_count);

  ## Convert counts to transition frequencies
  $self->normalize_transition_frequencies();
}

################################################################
=pod

=item B<two_words_update($added_word, $deleted_word)>

Update the model by adding a word and deleting another word. 

This invovles to update only a few transition frequencies : those including
the prefix of the added and deleted words. 

=cut
sub two_words_update {
    my ($self, $added_word, $deleted_word, $window_offset) = @_;

	my $pseudo_freq = $self->get_attribute("bg_pseudo");
	my @dna_alphabet =  qw (a c g t);
	
    ## No need to update if added word equald deleted word
    return(0) if ($added_word eq $deleted_word);
    
    ## Update transition count for the added word
    my $added_prefix = substr($added_word, 0, $self->{order});
    my $added_suffix = substr($added_word, $self->{order}, 1);
    $self->{oligo_counts}->{$added_prefix}->{$added_suffix}++;
     if (($self->{oligo_counts}->{$added_prefix}->{$added_suffix} == 1) 
	&&($main::verbose >= 4)){
	&RSAT::message::Warning(join (" ", "Model update:", $added_word, 
				      "appeared in updated window starting at", $window_offset));
    }
   

    ## Update transition count for the deleted word
    my $deleted_prefix = substr($deleted_word, 0, $self->{order});
    my $deleted_suffix = substr($deleted_word, $self->{order}, 1);
    $self->{oligo_counts}->{$deleted_prefix}->{$deleted_suffix}--;
     if (($self->{oligo_counts}->{$deleted_prefix}->{$deleted_suffix} == 0) 
	&&($main::verbose >= 4)){
	&RSAT::message::Warning(join (" ", "Model update:", $deleted_word, 
				      "disappeared from updated window starting at", $window_offset));
    }


	## Update relative frequencies for the added and deleted prefix
	## added prefix
	my $added_pattern_count = $self->{oligo_counts}->{$added_prefix}->{$added_suffix};	
	my $added_pattern_rel_freq = $added_pattern_count / $self->{training_words};
	$self->{oligo_freq}->{$added_prefix}->{$added_suffix} = $added_pattern_rel_freq;
	##pseudo-count
	my $added_pattern_pseudo_freq = 
	    ((1 - $pseudo_freq)*$self->{oligo_freq}->{$added_prefix}->{$added_suffix}) + $pseudo_freq/scalar(@dna_alphabet);
	$self->{oligo_freq}->{$added_prefix}->{$added_suffix} = $added_pattern_pseudo_freq;
	## prefix sum
	$self->{prefix_sum}->{$added_prefix} = 0;
	foreach my $suffix (sort keys (%{$self->{oligo_freq}->{$added_prefix}})) {
	    my $pattern_count = $self->{oligo_freq}->{$added_prefix}->{$suffix};	  
	    $self->{prefix_sum}->{$added_prefix} += $pattern_count;
	}
	
	## deleted prefix
	my $deleted_pattern_count = $self->{oligo_counts}->{$deleted_prefix}->{$deleted_suffix};	    
	my $deleted_pattern_rel_freq = $deleted_pattern_count / $self->{training_words};
	$self->{oligo_freq}->{$deleted_prefix}->{$deleted_suffix} = $deleted_pattern_rel_freq;
	##pseudo-count
	my $deleted_pattern_pseudo_freq = 
	    ((1 - $pseudo_freq)*$self->{oligo_freq}->{$deleted_prefix}->{$deleted_suffix}) + $pseudo_freq/scalar(@dna_alphabet);
	$self->{oligo_freq}->{$deleted_prefix}->{$deleted_suffix} = $deleted_pattern_pseudo_freq;
	## prefix sum
	$self->{prefix_sum}->{$deleted_prefix} = 0;
	foreach my $suffix (sort keys (%{$self->{oligo_freq}->{$deleted_prefix}})) {
	    my $pattern_count = $self->{oligo_freq}->{$deleted_prefix}->{$suffix};	    
	    $self->{prefix_sum}->{$deleted_prefix} += $pattern_count;
	}

    ## Update transition frequencies for the added and deleted prefix 
    foreach my $suffix (sort keys (%{$self->{oligo_freq}->{$added_prefix}})) {
    	$self->{transitions}->{$added_prefix}->{$suffix} = 
			$self->{oligo_freq}->{$added_prefix}->{$suffix}/$self->{prefix_sum}->{$added_prefix};	  
	}
	foreach my $suffix (sort keys (%{$self->{oligo_freq}->{$deleted_prefix}})) {
    	$self->{transitions}->{$deleted_prefix}->{$suffix} = 
			$self->{oligo_freq}->{$deleted_prefix}->{$suffix}/$self->{prefix_sum}->{$deleted_prefix};	  
	}
    
#     &RSAT::message::Debug("Updated model", 
# 			  "added",$added_word,
# 			  $added_prefix, 
# 			  $self->{prefix_sum}->{$added_prefix},
# 			  $added_suffix,
# 			  $self->{transitions}->{$added_prefix}->{$added_suffix},
# 			  "deleted", $deleted_word,
# 			  $deleted_prefix, 
# 			  $self->{prefix_sum}->{$deleted_prefix},
# 			  $deleted_suffix,
# 			  $self->{transitions}->{$deleted_prefix}->{$deleted_suffix},
# 			  ) if ($main::verbose >= 0);
}

################################################################
=pod

=item B<one_word_update($word, $mode,$window_offset)>

Update the model by adding or deleting a word depending on the choen mode. 

This invovles to update only a few transition frequencies : those including
the prefix of the word involved. 


=cut
sub one_word_update {
    my ($self, $word, $mode,$window_offset) = @_;

	my $pseudo_freq = $self->get_attribute("bg_pseudo");
	my @dna_alphabet =  qw (a c g t);
	
    ## Update transition count for the added word
    my $curr_prefix = substr($word, 0, $self->{order});
    my $curr_suffix = substr($word, $self->{order}, 1);
  
    if ($mode eq "add"){
    	$self->{oligo_counts}->{$curr_prefix}->{$curr_suffix}++;
     if (($self->{oligo_counts}->{$curr_prefix}->{$curr_suffix} == 1) 
		&&($main::verbose >= 4)){
		&RSAT::message::Warning(join (" ", "Model update: add", $word, 
				      "appeared in updated window starting at", $window_offset));
    }
    } elsif ($mode eq "delete"){
    	$self->{oligo_counts}->{$curr_prefix}->{$curr_suffix}--;
     if (($self->{oligo_counts}->{$curr_prefix}->{$curr_suffix} == 0) 
	&&($main::verbose >= 4)){
	&RSAT::message::Warning(join (" ", "Model update: delete", $word, 
				      "disappeared from updated window starting at", $window_offset));
    }
    }
    
    &RSAT::message::Debug("updated absolute count",Dumper($self->{oligo_counts})) if ($main::verbose >= 5);
	
	## Update relative frequencies for the added and deleted prefix
	## added prefix
	my $curr_pattern_count = $self->{oligo_counts}->{$curr_prefix}->{$curr_suffix};	
	my $curr_pattern_rel_freq = $curr_pattern_count / $self->{training_words};
	$self->{oligo_freq}->{$curr_prefix}->{$curr_suffix} = $curr_pattern_rel_freq;
	##pseudo-count
	my $curr_pattern_pseudo_freq = 
	    ((1 - $pseudo_freq)*$self->{oligo_freq}->{$curr_prefix}->{$curr_suffix}) + $pseudo_freq/scalar(@dna_alphabet);
	$self->{oligo_freq}->{$curr_prefix}->{$curr_suffix} = $curr_pattern_pseudo_freq;
	## prefix sum
	$self->{prefix_sum}->{$curr_prefix} = 0;
	foreach my $suffix (sort keys (%{$self->{oligo_freq}->{$curr_prefix}})) {
	    my $pattern_count = $self->{oligo_freq}->{$curr_prefix}->{$suffix};	  
	    $self->{prefix_sum}->{$curr_prefix} += $pattern_count;
	}


    ## Update transition frequencies for the added and deleted prefix 
    foreach my $suffix (sort keys (%{$self->{oligo_freq}->{$curr_prefix}})) {
    	$self->{transitions}->{$curr_prefix}->{$suffix} = 
			$self->{oligo_freq}->{$curr_prefix}->{$suffix}/$self->{prefix_sum}->{$curr_prefix};	  
	}
    
#     &RSAT::message::Debug("Updated model", 
# 			  "added",$added_word,
# 			  $added_prefix, 
# 			  $self->{prefix_sum}->{$added_prefix},
# 			  $added_suffix,
# 			  $self->{transitions}->{$added_prefix}->{$added_suffix},
# 			  "deleted", $deleted_word,
# 			  $deleted_prefix, 
# 			  $self->{prefix_sum}->{$deleted_prefix},
# 			  $deleted_suffix,
# 			  $self->{transitions}->{$deleted_prefix}->{$deleted_suffix},
# 			  ) if ($main::verbose >= 0);
}


################################################################
=pod

=item B<to_string($format, %args)>

Return a string describing the transition matrix.

Supported formats:

=over

=item tab

tab-delimited format

=item patser

Format supported by the programs patser and consensus, developed by
Gerald Z. Hertz.

=item MotifSampler

See input formats for description. 

=item oligos

See input formats for description. 

=back

=cut

sub to_string {
    my ($self, $format, %args) = @_;

    $format =~ s/oligo-analysis/oligos/;
    $format =~ s/^ms$/motifsampler/;

    if ($format eq ("tab")) {
      &RSAT::message::Warning("Output format tab is deprecated, please use format transitions instead.");
      $self->to_prefix_suffix_table(%args, type=>"transitions");
    } elsif ($format eq ("transitions")) {
      $self->to_prefix_suffix_table(%args, type=>"transitions");
    } elsif ($format eq ("tables")) {
      my $string = "";

      ## Original frequencies
      $string .= ";\n; Original oligomer frequencies\n";
      $string .= $self->to_prefix_suffix_table(%args, type=>"oligo_freq");

      ## Relative frequencies
      $string .= ";\n; Relative oligmer frequencies\n";
      $string .= $self->to_prefix_suffix_table(%args, type=>"oligo_freq_rel");

      ## Oligo counts
      if ($self->get_attribute("oligo_counts_defined")) {
	$string .= ";\n; Oligo counts\n";
	$string .= $self->to_prefix_suffix_table(%args, type=>"oligo_counts");
      }

      ## Transition table
      $string .= ";\n; Transition frequencies\n";
      $string .= $self->to_prefix_suffix_table(%args, type=>"transitions");

      ## Oligo frequencies derived from transitions
#      $string .= ";\n; Oligomer frequencies derived from transitions\n";
#      $string .= $self->to_string_oligos(%args);

    } elsif ($format eq ("oligos")) {
	$self->to_string_oligos(%args); 
    } elsif ($format eq ("meme")) {
	$self->to_string_meme(%args); 
    } elsif ($format eq ("motifsampler")) {
	$self->to_string_MotifSampler(%args); 
    } elsif ($format eq ("patser")) {
	$self->to_string_patser(%args); 
    } elsif ($format eq "") {
	&RSAT::error::FatalError(join("\t", "MarkovModel::to_string()", 
				      "the format of the backgroun model must be specified"));
    } else {
	&RSAT::error::FatalError(join("\t", "MarkovModel::to_string()", 
				      $format, "invalid output format for a background model. Supported:",
				      join(",", @supported_output_formats)));
    }
}


################################################################
=pod

=item B<to_prefix_suffix_table(%args)>

Export the Markov model in a tab-delimited format.

Supported arguments: comment_string, decimals

Usage:
  my $transition_table = $bg_model->to_prefix_suffix_table();
  my $transition_table = $bg_model->to_prefix_suffix_table(decimals=>3,comment_string="#");
  my $transition_table = $bg_model->to_prefix_suffix_table(type=>"transitions");
  my $transition_table = $bg_model->to_prefix_suffix_table(type=>"frequencies");

=cut
sub to_prefix_suffix_table{
    my ($self, %args) = @_;
    my $decimals = $args{decimals} || "5";
    my $string = "";
    my $type = $args{type} || "transitions";
    my %prefix_proba = $self->get_attribute("prefix_proba");
    my @prefix = sort($self->get_attribute("prefixes"));
    my %suffix_sum = $self->get_attribute("suffix_sum");
    my @suffix = sort($self->get_attribute("suffixes"));

    ## Print header
    $string .= join ("\t", ";pr\\suf",
		     @suffix,
		     "P(pre)",
		     "N(pre)",
		    );
    $string .= "\n";

    ## Print transition frequencies and sum and proba per prefix
    foreach my $prefix (@prefix) {
	if (defined($args{comment_string})) {
	    $string .= $args{comment_string};
	}
	$string .= $prefix;
	foreach my $suffix (@suffix) {
	  my $value =  $self->{$type}->{$prefix}->{$suffix};
	  if (&RSAT::util::IsInteger($value)) {
	    $string .= sprintf "\t%d", $value;
	  } else {
	    $string .= sprintf "\t%.${decimals}f", $value;
	  }
	}
	$string .= sprintf "\t%.${decimals}f", $self->{prefix_proba}->{$prefix};
	my $count = $self->{prefix_sum}->{$prefix};;
	if (&RSAT::util::IsNatural($count)) {
	  $string .= sprintf "\t%d", $count;
	} else {
	  $string .= sprintf "\t%.${decimals}f", $count;
	}
	$string .= "\n";
    }

    ## Print suffix probabilities
    my $row_name_len = &RSAT::stats::max(5,$self->get_attribute("order"));
    $string .= sprintf("; %${row_name_len}s", "P(su)");
#    $string .= "; P(su)";
    foreach my $suffix (@suffix) {
	$string .= sprintf "\t%.${decimals}f",  $self->{suffix_proba}->{$suffix};
    }
    $string .= "\n";

    ## Print suffix counts
    $string .= sprintf("; %${row_name_len}s", "N(su)");
#    $string .= "; N(su)";
    foreach my $suffix (@suffix) {
      my $count =   $self->{suffix_sum}->{$suffix};
      if (&RSAT::util::IsNatural($count)) {
	$string .= sprintf "\t%d",  $count;
      } else {
	$string .= sprintf "\t%.${decimals}f",  $count;
      }
    }
    $string .= "\n";


    return $string;
}


################################################################
=pod

=item B<to_string_oligos(%args)>

Convert the Markov model into oligomer frequencies.

A Markov model of order m will generate a frequency table for
oligomers of length k=m+1.

=cut
sub to_string_oligos {
  my ($self, %args) = @_;
  my $decimals = $args{decimals} || "8";
  my $string = "";
#  my %prefix_proba = $self->get_attribute("prefix_proba");
#  my %suffix_sum = $self->get_attribute("suffix_sum");

  my $order = $self->get_attribute("order");
  my $strand = $self->get_attribute("strand") || "sensitive";
  if ($strand eq "insensitive") {
    $string .= join("\t", "#seq", "seq|revcpl", "freq");
  } else {
    $string .= join("\t", "#seq", "id", "freq");
  }
  $string .= "\n";


  &RSAT::message::Info("Exporting Markov model", "order", $order, "strand", $strand) 
    if ($main::verbose >= 2);

  ## Compute oligomer frequencies from the prefix prior proba + transition proba
#  my @prefix = sort($self->get_attribute("prefixes"));
#  my @suffix = sort($self->get_attribute("suffixes"));
#  foreach my $prefix (sort (@prefix)) {
#    foreach my $suffix (sort(@suffix)) {
#      my $oligo_freq = $self->{prefix_proba}->{$prefix}*$self->{transitions}->{$prefix}->{$suffix};
  my %patterns = $self->transitions_to_oligo_frequencies();
  foreach my $oligo_seq (sort (keys(%patterns))) {
    if ($strand eq "insensitive") {
      my $oligo_rc = lc(&RSAT::SeqUtil::ReverseComplement($oligo_seq));;
#      next if ($oligo_rc lt $oligo_seq);
      $string .= $oligo_seq;
      $string .= "\t".$oligo_seq;
      $string .= "|";
      $string .= $oligo_rc;
    } else {
      $string .= $oligo_seq;
      $string .= "\t".$oligo_seq;
    }
    $string .= sprintf "\t%.${decimals}f",  $patterns{$oligo_seq}->{exp_freq};
    $string .= "\n";
  }
  return $string;
}


################################################################
=pod

=item B<to_string_meme(%args)>

Convert the Markov model into oligomer frequencies in the MEME
format. 

A Markov model of order m corresponds a frequency table for oligomers
of length k=m+1. The MEME format includes oligomer frequency tables
for all k-mers from k=1 to k=m+1.

=cut
sub to_string_meme {
  my ($self, %args) = @_;
  my $decimals = $args{decimals} || "3";
  my $string = "";
  my %prefix_proba = $self->get_attribute("prefix_proba");
  my @prefix = sort($self->get_attribute("prefixes"));
  my %suffix_sum = $self->get_attribute("suffix_sum");
  my @suffix = sort($self->get_attribute("suffixes"));

  my $order = $self->get_attribute("order");

  ## Print residue frequencies
  $string .= "# order 0\n";
  foreach my $suffix (sort @suffix) {
    $string .= sprintf "%s %.${decimals}e\n", uc($suffix),  $self->{suffix_proba}->{$suffix};
  }

  if ($order > 0) {
    ## Estimate subword frequencies from prefix frequencies
    foreach my $i (2..$order) {
      $string .= "# order ".($i-1)."\n";
      my %freq = ();
      foreach my $prefix (@prefix) {
	my $prefix_proba = $self->{prefix_proba}->{$prefix};
	my $subword = substr($prefix,0, $i);
	$freq{$subword} += $prefix_proba;
      }
      foreach my $oligo (sort keys %freq) {
	$string .= uc($oligo);
	$string .= sprintf " %.${decimals}e",  $freq{$oligo};
	$string .= "\n";
      }
    }

    ## Calculate oligomer frequencies from prefix+suffix frequencies
    $string .= "# order ".$order."\n";
    foreach my $prefix (sort (@prefix)) {
      foreach my $suffix (sort(@suffix)) {
	my $oligo = $prefix.$suffix;
	my $oligo_freq = $self->{prefix_proba}->{$prefix}*$self->{transitions}->{$prefix}->{$suffix};
	$string .= uc($oligo);
	$string .= sprintf " %.${decimals}e",  $oligo_freq;
	$string .= "\n";
      }
    }
  }
  return $string;
}


################################################################
=pod

=item B<to_string_MotifSampler(%args)>

Export the Markov model in  MotifSampler format.

Supported arguments: comment_string, decimals

=cut
sub to_string_MotifSampler {
    my ($self, %args) = @_;
    my $decimals = $args{decimals} || "5";
    my $string = "";
    my %prefix_proba = $self->get_attribute("prefix_proba");
    my @prefix = sort($self->get_attribute("prefixes"));

    my %suffix_sum = $self->get_attribute("suffix_sum");
    my @suffix = sort($self->get_attribute("suffixes"));

    ## Print header
    $string .= join ("\n", 
		     "#INCLUSive Background Model v1.0",
		     "#",
		     "#Order = ".$self->get_attribute("order"),
		     "#Organism = unknown",
		     "#Sequences = ",
		     "#Path = ",
		     "#");
    $string .= "\n";

    ## Single nucleotide frequencies
    $string .= "\n#snf\n";
    foreach my $suffix (@suffix) {
      push @suffix_proba, sprintf "%.${decimals}f",  $self->{suffix_proba}->{$suffix};
    }
    $string .= join ("\t", @suffix_proba);
    $string .= "\n";

    ## Prefix probabilities

    ## Particularity of MotifSampler format for oder 0: all oligos are exported as prefix
    $string .= "\n#oligo frequency\n";
    if ($self->{order}==0) {
      foreach my $suffix (@suffix) {
	push @prefix_proba, sprintf "%.${decimals}f",  $self->{suffix_proba}->{$suffix};
      }
    } else {
      foreach my $prefix (@prefix) {
	push @prefix_proba, sprintf "%.${decimals}f",  $self->{prefix_proba}->{$prefix};
      }
    }
    $string .= join ("\n", @prefix_proba);
    $string .= "\n";

    ## Print transition frequencies and sum and proba per prefix

    ## Particularity of MotifSampler format for oder 0: the oligo frequencies are repeated 4 times
    if ($self->{order}==0) {
      @prefix = ("","","","");
    }

    $string .= "\n#transition matrix\n";
    foreach my $prefix (@prefix) {
      my @transitions = ();
      foreach my $suffix (@suffix) {
	push @transitions, sprintf "%.${decimals}f",  $self->{transitions}->{$prefix}->{$suffix};
      }
      $string .= join ("\t", @transitions);
      $string .= "\n";
    }

    return $string;
}


################################################################
=pod

=item B<to_string_patser>

Export the Markov model in patser format.  

WARNING: patser only supports Bernoulli models ! If the Markov order is
superior to 0, this function issues a fatal error.

=cut
sub to_string_patser {
    my ($self, %args) = @_;
    my $string = "";
    my $decimals = $args{decimals} || "5";
    my @prefix = sort($self->get_attribute("prefixes"));
    my @suffix = sort($self->get_attribute("suffixes"));
    
    ## Check that the model is Bernoulli (Markov order = 0)
    &RSAT::error::FatalError(join("\t",
				  "MarkovModel::to_string_patser()",
				  "The patser format only supports Bernoulli models (Markov order must be 0)"
				  ))
	unless ($self->get_attribute("order") == 0);
    
    ## The output ormat differs between strand-sensitive and strand-insensitive models
    my $strand = $self->get_attribute("strand");
    if ($strand eq "insensitive") {
	foreach my $prefix (@prefix) {
 	    foreach my $suffix (@suffix) {
		my $suffix_rc = lc(&main::SmartRC($suffix));

		next if ($suffix gt $suffix_rc);
		if (defined($args{comment_string})) {
		    $string .= $args{comment_string};
		}
		$string .= $suffix.":".$suffix_rc;
		$string .= sprintf ("\t%.${decimals}f\n",$self->{transitions}->{$prefix}->{$suffix});
	    }
	}
	
    } else {
	foreach my $prefix (@prefix) {
	    foreach my $suffix (@suffix) {
		if (defined($args{comment_string})) {
		    $string .= $args{comment_string};
		}
		$string .= $suffix;
		$string .= sprintf ("\t%.${decimals}f\n",$self->{transitions}->{$prefix}->{$suffix});
	    }
	}
	
    }
    return $string;

}


################################################################
=pod

=item B<average_strands>

Average transition frequencies between pairs of reverse complements,
in order to obtain a strand-insensitive model.

The averaging is performed in 3 steps

 - convert the transition matrix (order m) into oligo frequencies
   (k=m+1).
 - average oligo frequencies
 - convert averaged oligo frequencies into transition matrix

=cut
sub average_strands {
  my ($self) = @_;
    &RSAT::message::Info(join("\t", "MarkovModel", "Averaging transition frequencies on both strands."))
	if ($main::verbose >= 0);

  ## Get (supposedly strand-sensitive) oligo frequencies from the
  ## transition matrix
  my %patterns_1str = $self->transitions_to_oligo_frequencies();

  ## Compute strand-insensitive oligo frequencies
  my %patterns_2str = ();
  foreach my $pattern (keys(%patterns_1str)) {
    my $pattern_rc = lc(&main::SmartRC($pattern));
    if ($pattern lt $pattern_rc) {
      $patterns_2str{$pattern}->{exp_freq} =
	$patterns_2str{$pattern_rc}->{exp_freq} =
	($patterns_1str{$pattern}->{exp_freq} +
	  $patterns_1str{$pattern_rc}->{exp_freq})/2;
    } elsif ($pattern eq $pattern_rc) {
      $patterns_2str{$pattern}->{exp_freq} = $patterns_1str{$pattern}->{exp_freq};
    }
  }
  $self->oligos_to_frequency_table(%patterns_2str);
  $self->force_attribute("strand", "insensitive");
  $self->normalize_transition_frequencies(no_pseudo=>1);
}

################################################################
=pod

=item B<average_strands>

Average transition frequencies between pairs of reverse complements,
in order to obtain a strand-insensitive model.

=cut
# sub average_strands_bof {
#     my ($self) = @_;

#     &RSAT::message::Info(join("\t", "MarkovModel", "Averaging transition frequencies on both strands."))
# 	if ($main::verbose >= 0);

#     ## Sum per prefix
#     my %prefixes_2str = ();
#     foreach my $prefix ($self->get_attribute("prefixes")) {
# 	next if (defined($prefixes_2str{$prefix}));
# 	my $prefix_rc = lc(&main::SmartRC($prefix));

# 	## Make sure that both the prefix and its reverse complement are indexed
# 	$prefixes_2str{$prefix} = 1;
# 	$prefixes_2str{$prefix_rc} = 1;
# #	&RSAT::message::Debug("average_strands sum for prefix", $prefix, $prefix_rc) if ($main::verbose >= 0);

# 	unless (defined($self->{prefix_sum}->{$prefix} )) {
# 	    $self->{prefix_sum}->{$prefix}  = 0;
# 	}
# 	unless (defined($self->{prefix_sum}->{$prefix_rc} )) {
# 	    $self->{prefix_sum}->{$prefix_rc}  = 0;
# 	}
# 	$self->{prefix_sum}->{$prefix} 
# 	= $self->{prefix_sum}->{$prefix_rc}
# 	= ($self->{prefix_sum}->{$prefix} +$self->{prefix_sum}->{$prefix_rc})/2 ;

# 	unless (defined($self->{prefix_proba}->{$prefix} )) {
# 	    $self->{prefix_proba}->{$prefix}  = 0;
# 	}
# 	unless (defined($self->{prefix_proba}->{$prefix_rc} )) {
# 	    $self->{prefix_proba}->{$prefix_rc}  = 0;
# 	}
# 	$self->{prefix_proba}->{$prefix} 
# 	= $self->{prefix_proba}->{$prefix_rc}
# 	= ($self->{prefix_proba}->{$prefix} + $self->{prefix_proba}->{$prefix_rc})/2 ;
#     }
#     $self->set_array_attribute("prefixes", sort keys %prefixes_2str);

#     ## Sum per suffix
#     my %suffixes_2str = ();
#     foreach my $suffix ($self->get_attribute("suffixes")) {
# 	next if (defined($suffixes_2str{$suffix}));
# 	my $suffix_rc = lc(&main::SmartRC($suffix));

# 	## Make sure that both the suffix and its reverse complement are indexed
# 	$suffixes_2str{$suffix} = 1;
# 	$suffixes_2str{$suffix_rc} = 1;
# #	&RSAT::message::Debug("average_strands sum for suffix", $suffix, $suffix_rc) if ($main::verbose >= 0);

# 	$self->{suffix_sum}->{$suffix} 
# 	= $self->{suffix_sum}->{$suffix_rc}
# 	= ($self->{suffix_sum}->{$suffix} +$self->{suffix_sum}->{$suffix_rc})/2 ;

# 	$self->{suffix_proba}->{$suffix} 
# 	= $self->{suffix_proba}->{$suffix_rc}
# 	= ($self->{suffix_proba}->{$suffix} +$self->{suffix_proba}->{$suffix_rc})/2 ;
#     }
#     $self->set_array_attribute("suffixes", sort keys %suffixes_2str);

#     ## Transition matrix
#     foreach my $prefix ($self->get_attribute("prefixes")) {
# #	my $prefix_rc = lc(&main::SmartRC($prefix));
# #	next if ($prefix_rc ge $prefix);
# 	foreach my $suffix ($self->get_attribute("suffixes")) {
# 	    my $word = $prefix.$suffix;
# 	    my $rc_word = lc(&main::SmartRC($word));
# 	    next if ($rc_word ge $word);
# 	    my $rc_prefix = substr($rc_word, 0, $self->{order});
# 	    my $rc_suffix = substr($rc_word, $self->{order}, 1);
# 	    unless (defined($self->{transitions}->{$prefix}->{$suffix})) {
# 		$self->{transitions}->{$prefix}->{$suffix} = 0;
# 	    }
# 	    unless (defined($self->{transitions}->{$rc_prefix}->{$rc_suffix})) {
# 		$self->{transitions}->{$rc_prefix}->{$rc_suffix} = 0;
# 	    }
# 	    $self->{transitions}->{$prefix}->{$suffix} 
# 	    = $self->{transitions}->{$rc_prefix}->{$rc_suffix} 
# 	    = ($self->{transitions}->{$prefix}->{$suffix} +
# 	       $self->{transitions}->{$rc_prefix}->{$rc_suffix})/2;
	    
# 	    unless (defined($self->{oligo_freq}->{$prefix}->{$suffix})) {
# 		$self->{oligo_freq}->{$prefix}->{$suffix} = 0;
# 	    }
# 	    unless (defined($self->{oligo_freq}->{$rc_prefix}->{$rc_suffix})) {
# 		$self->{oligo_freq}->{$rc_prefix}->{$rc_suffix} = 0;
# 	    }
# 	    $self->{oligo_freq}->{$prefix}->{$suffix} 
# 	    = $self->{oligo_freq}->{$rc_prefix}->{$rc_suffix} 
# 	    = ($self->{oligo_freq}->{$prefix}->{$suffix} +
# 	       $self->{oligo_freq}->{$rc_prefix}->{$rc_suffix})/2;

# #	    &RSAT::message::Debug("average_strands transitions", $prefix, $suffix, $rc_word, $rc_prefix, $rc_suffix) if ($main::verbose >= 0);
# 	}
#     }
#     $self->force_attribute("strand", "insensitive");
# #    die "HELLO";
# }


################################################################
=pod 

=item B<segment_proba($segment)>

Calculate the probability of a segment of sequence. The sequence segment must
by larger than the Markov order + 1.

=cut

sub segment_proba {
    my ($self, $segment) = @_;

#    return(1);

    my $seq_len = length($segment);
    my $order = $self->get_attribute("order");
    $segment =  lc($segment);

    &RSAT::error::FatalError("&RSAT::MarkovModel->segment_proba. The segment ($segment) length ($seq_len) must be larger than the markov order ($order) + 1.") 
	if ($seq_len < $order + 1);

    my $prefix = substr($segment,0,$order);
    my $segment_proba = 0;
    my $c = 0;
    if (defined($self->{prefix_proba}->{$prefix})) {
      $segment_proba = $self->{prefix_proba}->{$prefix};
      &RSAT::message::Info(#"MarkovModel::segment_proba()",
			   "offset:".$c,
			   "i=".($c+1),
			   "P(".$prefix.")", $self->{prefix_proba}->{$prefix},
			   "P(S)", $segment_proba,
			   $prefix.uc($suffix),
			   $prefix.uc($suffix),
			  ) if ($main::verbose >= 4);
    } elsif ($self->get_attribute("n_treatment") eq "score") {
      $segment_proba = 1;
      foreach my $i (1..length($prefix)) {
	my $residue = substr($prefix,$i-1,1);
	if (defined($self->{suffix_proba}->{$residue})) {
	  $segment_proba *= $self->{suffix_proba}->{$residue};
	}
      }

      #	&RSAT::message::Debug("Ignoring undefined prefix", $prefix,  "proba set to 1") if ($main::verbose >= #0);
    } else {
      &RSAT::error::FatalError("\t", "MarkovModel::segment_proba",
			       "Invalid prefix for the selected sequence type", $prefix);
    }
    
    for $c ($order..($seq_len-1)) {
	my $residue_proba = 0;
	
	my $suffix = substr($segment, $c, 1);
	my $prefix = substr($segment,($c-$order),$order);
	if (defined($self->{transitions}->{$prefix}->{$suffix})) {
	  $residue_proba = $self->{transitions}->{$prefix}->{$suffix};
	  if ($residue_proba <= 0) {
	    &RSAT::error::FatalError(join("\t", "MarkovModel::segment_proba",
					  "null transition between prefix ", $prefix, " and suffix", $suffix)) if ($main::verbose >= 0);
	  }
	} elsif ($self->get_attribute("n_treatment") eq "score") {
	  if (defined($self->{suffix_proba}->{$suffix})) {
	    $residue_proba = $self->{suffix_proba}->{$suffix};
	  } else {
	    $residue_proba = 1;
	  }
	  #	    &RSAT::message::Debug("Ignoring undefined transition", $prefix, $suffix, "proba set to 1") if ($main::verbose >= 10);
	} else { 
	  &RSAT::error::FatalError(join("\t", "MarkovModel::segment_proba",
					"undefined transition between prefix ", 
					$prefix, "and suffix", 
					$suffix));
	}
	$segment_proba *= $residue_proba;
	
	&RSAT::message::Info(#"MarkovModel::segment_proba()",
			     "offset:".$c,
			     "i=".($c+1),
			     "P(".$suffix."|".$prefix.")", $residue_proba,
			     "P(S)", $segment_proba,
			     $prefix.uc($suffix),
			     substr($segment,0,$c+1),
			    ) if ($main::verbose >= 4);
#	&RSAT::message::Info("segment_proba", 
#			      "prefix=".$word, 
#			      "prefix=".$prefix, 
#			      "suffix:".$suffix, 
#			      "offset:".$c, 
#			      "P(letter)=".$residue_proba, 
#			      "P(S)=".$segment_proba) if ($main::verbose >= 4);
    }
    
#    &RSAT::message::Debug("segment_proba", $segment, "P(segm)=".$segment_proba) if ($main::verbose >= 0);
#    die;
    return $segment_proba;
}


return 1;

__END__

