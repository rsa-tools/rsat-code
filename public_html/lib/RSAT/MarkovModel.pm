###############################################################
#
# Class MarkovModel
#
package RSAT::MarkovModel;

## CVS: fixed bugs with the method average_strands()

## TO DO : 
## - for the model update, take into account the strand-insensitive models
## - Would it be appropriate to have a pseudo-count for the estimation
##  of markov chain models from input sequenes ?

use RSAT::GenericObject;
use RSAT::SeqUtil;
use RSAT::error;
@ISA = qw( RSAT::GenericObject );

### class attributes
%supported_input_formats = ("oligo-analysis"=>1);
%supported_output_formats = ("tab"=>1, 
			     "patser"=>1);
@supported_input_formats = sort keys %supported_input_formats;
@supported_output_formats = sort keys %supported_output_formats;

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
    if ($format eq "oligo-analysis") {
	$self->load_from_file_oligos($bg_file);
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

The model file is a tab-delimited text file indicating the expected frequency
of each word of length m+1, for a model of order m. These files can be
produced with the ccmmand oligo-analysis. When new genomes are installed,
background files are automatically calculated for upstream sequences, and
stored in $RSAT/data/genomes/[Organism_name]/oligo-frequencies.

In order to ensure case-insensivity, Markov models are automatically
converted to lowercases.

=cut
sub load_from_file_oligos {
    my ($self, $bg_file) = @_;

    &RSAT::message::TimeWarn("Loading Markov model from file $bg_file") if ($main::verbose >= 2);

    my ($file_type, %patterns) = &main::ReadPatternFrequencies($bg_file) ;
    my @patterns = keys(%patterns);


    ## This is a bit tricky: ReadExpectedFrequencies sets a global variable
    ## $file_type to "2str" if the model is strand-insensitive/
    if ($file_type eq "2str") {
	$self->force_attribute("strand", "insensitive");
    }

    my $order = length($patterns[0]) -1 ; ## Use the first pattern to calculate model order
    $self->force_attribute("order", $order);


    ## Calculate alphabet from expected frequency keys
    foreach my $pattern_seq (keys %patterns) {
	## Check pattern length
	$pattern_seq = lc($pattern_seq);
	my $pattern_len = length($pattern_seq);
	my $pattern_freq =  $patterns{$pattern_seq}->{exp_freq};

#	&RSAT::message::Debug($pattern_seq, $pattern_len, $pattern_freq, 
#			      join(";", keys %{$patterns{$pattern_seq}})) if ($main::verbose >= 0);
	&RSAT::error::FatalError("All patterns should have the same length in a Markov model file.") 
	    unless $pattern_len = $order+1;
	my $prefix = substr($pattern_seq,0,$order);
	my $suffix = substr($pattern_seq,$order, 1);
	$self->{transition_count}->{$prefix}->{$suffix} = $pattern_freq;
#	&RSAT::message::Debug("transition count", $prefix.".".$suffix, $pattern_freq, $self->{transition_count}->{$prefix}->{$suffix}) if ($main::verbose >= 10);
    }
    
#    &RSAT::message::Debug("MARKOV MODEL", $order, join (' ', @patterns)) if ($main::verbose >= 5);
}


################################################################
=pod

=item B<normalize_transition_frequencies>

Normalize transition frequencies, i.e. make sure that, for each prfix, the sum
of transition frequencies (frequency of all possible siffuxed given the
prefix) is 1.

=cut
    
sub normalize_transition_frequencies {
    my ($self) = @_;
    
    &RSAT::message::TimeWarn(join("\t", "MarkovModel", "Normalizing transition frequencies")) if ($main::verbose >= 3);

    ## Calculate sum of counts (frequencies) per prefix and suffix
    my %prefix_sum = ();
    my %suffix_sum = ();
    my $freq_sum = 0;
    my $p=0;
    my $s=0;


#    foreach my $prefix (keys(%prefix_sum)) {
    foreach my $prefix (sort keys(%{$self->{transition_count}})) {
	$p++;
#	foreach my $suffix (keys (%suffix_sum)) {
	foreach my $suffix (sort keys (%{$self->{transition_count}->{$prefix}})) {
	    $s++;
	    my $pattern_count = $self->{transition_count}->{$prefix}->{$suffix};	    
	    $prefix_sum{$prefix} += $pattern_count;
	    $freq_sum += $pattern_count;
	    $suffix_sum{$suffix} += $pattern_count;
#	    &RSAT::message::Debug("Count sum", $p, $prefix.".".$suffix, $pattern_count, $freq_sum) if ($main::verbose >= 0);
	}
    }
    $self->set_hash_attribute("prefix_sum", %prefix_sum);
    $self->set_hash_attribute("suffix_sum", %suffix_sum);

    ## Store prefixes and suffixes in arrays for quick access
    $self->set_array_attribute("prefixes", sort(keys(%prefix_sum)));
    $self->set_array_attribute("suffixes", sort(keys(%suffix_sum)));
    
    ## Calculate transition frequencies
    my $missing_transitions = 0;
    foreach my $prefix ($self->get_attribute("prefixes")) {
	foreach my $suffix ($self->get_attribute("suffixes")) {
#	foreach my $suffix (keys (%suffix_sum)) {
	    if (defined($self->{transition_count}->{$prefix}->{$suffix})) {	
		if ($prefix_sum{$prefix} > 0) {
		    $self->{transition_freq}->{$prefix}->{$suffix} =  
			$self->{transition_count}->{$prefix}->{$suffix}/$prefix_sum{$prefix};
		} else {
		    $self->{transition_freq}->{$prefix}->{$suffix} =  0;
		}
	    } else {
		$missing_transitions++;
		$self->{transition_freq}->{$prefix}->{$suffix} = 0;
#		&RSAT::message::Warning(join(" ",
#					     "No transition between prefix",$prefix, 
#					     "and suffix", $suffix)) if ($main::verbose >= 3);
	    };
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

#    &RSAT::message::Debug("SUFFIX PROBA", $self->set_hash_attribute("suffix_proba", %suffix_proba)) if ($main::verbose >= 5);

    ## Average counts and frequencies on both strands if required
    my $strand = $self->get_attribute("strand") || "sensitive";
    if ($strand eq "insensitive") {
	$self->average_strands();
    }

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
	    unless ((defined($self->{transition_freq}->{$prefix}->{$suffix})) 
		    && ($self->{transition_freq}->{$prefix}->{$suffix} > 0)) {
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
    foreach my $prefix (keys(%{$self->{transition_count}})) {
#    foreach my $prefix ($self->get_attribute("prefixes")) {
      if ($prefix =~ /${accepted_expression}/) {
	$checked_prefixes{$prefix} = 1;
      } else {
	$suppressed_prefixes{$prefix}++;
#	delete $self->{transition_freq}->{$prefix};
	delete $self->{transition_count}->{$prefix};
	&RSAT::message::Info(join(" ",
				     "Supressing prefix",$prefix, 
				     "from transition matrix")) if ($main::verbose >= 5);
	next;
      }
      foreach my $suffix (keys(%{$self->{transition_count}->{$prefix}})) {
#      foreach my $suffix ($self->get_attribute("suffixes")) {
	$suppressed_suffixes{$suffix}++;
	if ($suffix =~ /${accepted_expression}/) {
	  $checked_suffixes{$suffix} = 1;
	} else {
	  delete $self->{transition_freq}->{$prefix}->{$suffix};
	  delete $self->{transition_count}->{$prefix}->{$suffix};
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


#    &RSAT::message::Debug("Current prefixes", join (",", $self->get_attribute("prefixes")),
#			  "freq", join(",", sort keys(%{$self->{transition_freq}})),
#			  "counts", join(",", sort keys(%{$self->{transition_count}})),
#			  "Current suffixes", join (",", $self->get_attribute("suffixes")),
#			   ) if ($main::verbose >= 10);
}


################################################################
=pod

=item B<calc_from_seq($sequence, [add=>0|1])>

Calculate background model from a sequence.

If the argument add=>1 is specified, the new sequence is added to the background model.

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
	$self->set_hash_attribute("transition_count",());
	if ($main::verbose >= 2) {
	  &RSAT::message::TimeWarn(join(" ", 
					"Calculating markov model (order ".$order.")",
					"from sequence of length", $seq_len));
	  &RSAT::message::Warning(join("", "RSAT::MarkovModel->calc_from_seq: sequence is too short (len=",$seq_len,
				       "). It must be larger than Markov order + 1 = ", $order+1))
	    unless ($seq_len > $order + 1);
	}
      }

    
    ## Transition counts
    my $last_pos = $seq_len - $order;
    for my $offset (0..($last_pos-1)) {
	my $prefix = substr($sequence, $offset, $order);
	my $suffix = substr($sequence, $offset + $order,1);
	$self->{transition_count}->{$prefix}->{$suffix}++;
    }
    $self->force_attribute("training_words", $last_pos + $previous_words);

    ## Delete transitions between letters which do not belong to the accepted alphabet
    $self->check_transition_alphabet();

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

    ## No need to update if added word equald deleted word
    return(0) if ($added_word eq $deleted_word) ;
    
    ## Update transition count for the added word
    my $added_prefix = substr($added_word, 0, $self->{order});
    my $added_suffix = substr($added_word, $self->{order}, 1);
    $self->{transition_count}->{$added_prefix}->{$added_suffix}++;
    $self->{prefix_sum}->{$added_prefix}++;
    if (($self->{transition_count}->{$added_prefix}->{$added_suffix} == 1) 
	&&($main::verbose >= 4)){
	&RSAT::message::Warning(join (" ", "Model update:", $added_word, 
				      "appeared in updated window starting at", $window_offset));
    }

    ## Update transition count for the deleted word
    my $deleted_prefix = substr($deleted_word, 0, $self->{order});
    my $deleted_suffix = substr($deleted_word, $self->{order}, 1);
    $self->{transition_count}->{$deleted_prefix}->{$deleted_suffix}--;
    if (($self->{transition_count}->{$deleted_prefix}->{$deleted_suffix} == 0) 
	&&($main::verbose >= 4)){
	&RSAT::message::Warning(join (" ", "Model update:", $deleted_word, 
				      "disappeared from updated window starting at", $window_offset));
    }
    $self->{prefix_sum}->{$deleted_prefix}--;

    ## Update transition frequencies for the added and deleted prefix
    foreach my $suffix ($self->{suffixes}) {
	$self->{transition_freq}->{$added_prefix}->{$suffix} = 
	    $self->{transition_count}->{$added_prefix}->{$added_suffix}/$self->{prefix_sum}->{$added_prefix}++;
	$self->{transition_freq}->{$deleted_prefix}->{$suffix} = 
	    $self->{transition_count}->{$deleted_prefix}->{$deleted_suffix}/$self->{prefix_sum}->{$deleted_prefix}++;
    }

    
    
#     &RSAT::message::Debug("Updated model", 
# 			  "added",$added_word,
# 			  $added_prefix, 
# 			  $self->{prefix_sum}->{$added_prefix},
# 			  $added_suffix,
# 			  $self->{transition_freq}->{$added_prefix}->{$added_suffix},
# 			  "deleted", $deleted_word,
# 			  $deleted_prefix, 
# 			  $self->{prefix_sum}->{$deleted_prefix},
# 			  $deleted_suffix,
# 			  $self->{transition_freq}->{$deleted_prefix}->{$deleted_suffix},
# 			  ) if ($main::verbose >= 10);
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

=back

=cut

sub to_string {
    my ($self, $format, %args) = @_;
    if ($format eq ("tab")) {
	$self->to_string_tab(%args);
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

=item B<to_string_tab(%args)>

Export the Markov model in a tab-delimited format.

Supported arguments: comment_string, decimals

=cut
sub to_string_tab {
    my ($self, %args) = @_;
    my $decimals = $args{decimals} || "5";
    my $string = "";
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
	    $string .= sprintf "\t%.${decimals}f",  $self->{transition_freq}->{$prefix}->{$suffix};
	}
	$string .= sprintf "\t%.${decimals}f", $self->{prefix_proba}->{$prefix};
	$string .= sprintf "\t%.${decimals}f", $self->{prefix_sum}->{$prefix};
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
	$string .= sprintf "\t%.${decimals}f",  $self->{suffix_sum}->{$suffix};
    }
    $string .= "\n";


    return $string;
}


################################################################
## Export the Markov model in patser format
## Beware:  patser only supports Bernoulli models !
## If the Markov order is superior to 0, this function issues a fatal error.
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
		$string .= sprintf ("\t%.${decimals}f\n",$self->{transition_freq}->{$prefix}->{$suffix});
	    }
	}
	
    } else {
	foreach my $prefix (@prefix) {
	    foreach my $suffix (@suffix) {
		if (defined($args{comment_string})) {
		    $string .= $args{comment_string};
		}
		$string .= $suffix;
		$string .= sprintf ("\t%.${decimals}f\n",$self->{transition_freq}->{$prefix}->{$suffix});
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

=cut

sub average_strands {
    my ($self) = @_;

    &RSAT::message::Info(join("\t", "MarkovModel", "Averaging transition frequencies on both strands."))
	if ($main::verbose >= 3);

    ## Sum per prefix
    my %prefixes_2str = ();
    foreach my $prefix ($self->get_attribute("prefixes")) {
	next if (defined($prefixes_2str{$prefix}));
	my $prefix_rc = lc(&main::SmartRC($prefix));

	## Make sure that both the prefix and its reverse complement are indexed
	$prefixes_2str{$prefix} = 1;
	$prefixes_2str{$prefix_rc} = 1;
#	&RSAT::message::Debug("average_strands sum for prefix", $prefix, $prefix_rc) if ($main::verbose >= 0);

	unless (defined($self->{prefix_sum}->{$prefix} )) {
	    $self->{prefix_sum}->{$prefix}  = 0;
	}
	unless (defined($self->{prefix_sum}->{$prefix_rc} )) {
	    $self->{prefix_sum}->{$prefix_rc}  = 0;
	}
	$self->{prefix_sum}->{$prefix} 
	= $self->{prefix_sum}->{$prefix_rc}
	= ($self->{prefix_sum}->{$prefix} +$self->{prefix_sum}->{$prefix_rc})/2 ;

	unless (defined($self->{prefix_proba}->{$prefix} )) {
	    $self->{prefix_proba}->{$prefix}  = 0;
	}
	unless (defined($self->{prefix_proba}->{$prefix_rc} )) {
	    $self->{prefix_proba}->{$prefix_rc}  = 0;
	}
	$self->{prefix_proba}->{$prefix} 
	= $self->{prefix_proba}->{$prefix_rc}
	= ($self->{prefix_proba}->{$prefix} + $self->{prefix_proba}->{$prefix_rc})/2 ;
    }
    $self->set_array_attribute("prefixes", sort keys %prefixes_2str);

    ## Sum per suffix
    my %suffixes_2str = ();
    foreach my $suffix ($self->get_attribute("suffixes")) {
	next if (defined($suffixes_2str{$suffix}));
	my $suffix_rc = lc(&main::SmartRC($suffix));

	## Make sure that both the suffix and its reverse complement are indexed
	$suffixes_2str{$suffix} = 1;
	$suffixes_2str{$suffix_rc} = 1;
#	&RSAT::message::Debug("average_strands sum for suffix", $suffix, $suffix_rc) if ($main::verbose >= 0);

	$self->{suffix_sum}->{$suffix} 
	= $self->{suffix_sum}->{$suffix_rc}
	= ($self->{suffix_sum}->{$suffix} +$self->{suffix_sum}->{$suffix_rc})/2 ;

	$self->{suffix_proba}->{$suffix} 
	= $self->{suffix_proba}->{$suffix_rc}
	= ($self->{suffix_proba}->{$suffix} +$self->{suffix_proba}->{$suffix_rc})/2 ;
    }
    $self->set_array_attribute("suffixes", sort keys %suffixes_2str);

    ## Transition matrix
    foreach my $prefix ($self->get_attribute("prefixes")) {
#	my $prefix_rc = lc(&main::SmartRC($prefix));
#	next if ($prefix_rc ge $prefix);
	foreach my $suffix ($self->get_attribute("suffixes")) {
	    my $word = $prefix.$suffix;
	    my $rc_word = lc(&main::SmartRC($word));
	    next if ($rc_word ge $word);
	    my $rc_prefix = substr($rc_word, 0, $self->{order});
	    my $rc_suffix = substr($rc_word, $self->{order}, 1);
	    unless (defined($self->{transition_freq}->{$prefix}->{$suffix})) {
		$self->{transition_freq}->{$prefix}->{$suffix} = 0;
	    }
	    unless (defined($self->{transition_freq}->{$rc_prefix}->{$rc_suffix})) {
		$self->{transition_freq}->{$rc_prefix}->{$rc_suffix} = 0;
	    }
	    $self->{transition_freq}->{$prefix}->{$suffix} 
	    = $self->{transition_freq}->{$rc_prefix}->{$rc_suffix} 
	    = ($self->{transition_freq}->{$prefix}->{$suffix} +
	       $self->{transition_freq}->{$rc_prefix}->{$rc_suffix})/2;
	    
	    unless (defined($self->{transition_count}->{$prefix}->{$suffix})) {
		$self->{transition_count}->{$prefix}->{$suffix} = 0;
	    }
	    unless (defined($self->{transition_count}->{$rc_prefix}->{$rc_suffix})) {
		$self->{transition_count}->{$rc_prefix}->{$rc_suffix} = 0;
	    }
	    $self->{transition_count}->{$prefix}->{$suffix} 
	    = $self->{transition_count}->{$rc_prefix}->{$rc_suffix} 
	    = ($self->{transition_count}->{$prefix}->{$suffix} +
	       $self->{transition_count}->{$rc_prefix}->{$rc_suffix})/2;

#	    &RSAT::message::Debug("average_strands transitions", $prefix, $suffix, $rc_word, $rc_prefix, $rc_suffix) if ($main::verbose >= 0);
	}
    }
#    die "HELLO";
}


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
	if (defined($self->{transition_freq}->{$prefix}->{$suffix})) {
	  $residue_proba = $self->{transition_freq}->{$prefix}->{$suffix};
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

