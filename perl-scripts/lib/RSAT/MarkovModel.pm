###############################################################
#
# Class MarkovModel
#
package RSAT::MarkovModel;

## TO DO : for the model update, take into account the strand-insensitive models

use RSAT::GenericObject;
use RSAT::error;
@ISA = qw( RSAT::GenericObject );

### class attributes

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

=item B<load_from_file($bg_file)>

Load the markov model from a file. 

The model file is a tab-delimited text file indicating the expected frequency
of each word of length m+1, for a model of order m. These files can be
produced with the ccmmand oligo-analysis. When new genomes are installed,
background files are automatically calculated for upstream sequences, and
stored in $RSAT/data/genomes/[Organism_name]/oligo-frequencies.


=cut
sub load_from_file {
    my ($self, $bg_file) = @_;
    

    &RSAT::message::TimeWarn("Loading Markov model from file $bg_file") if ($main::verbose >= 2);

    my %patterns = &main::ReadExpectedFrequencies($bg_file) ;
    my @patterns = keys(%patterns);

    my $order = length($patterns[0]) -1 ; ## Use the first pattern to calculate model order
    $self->force_attribute("order", $order);

    ## Calculate alphabet from expected frequency keys
    foreach my $pattern_seq (keys %patterns) {
	## Check pattern length
	$pattern_seq = lc($pattern_seq);
	my $pattern_len = length($pattern_seq);
	my $pattern_freq =  $patterns{$pattern_seq}->{exp_freq};
	&RSAT::error::FatalError("All patterns should have the same length in a Markov model file.") 
	    unless $pattern_len = $order+1;
	my $prefix = substr($pattern_seq,0,$order);
	my $suffix = substr($pattern_seq,$order, 1);
	$self->{transition_count}->{$prefix}->{$suffix} = $pattern_freq;
    }
    
    $self->normalize_transition_frequencies();
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

    ## Calculate sum of counts (frequencies) per prefix and suffix
    my %prefix_sum = ();
    my %suffix_sum = ();
    my $freq_sum = 0;
    my $p=0;
    my $s=0;
#    foreach my $prefix (keys(%prefix_sum)) {
    foreach my $prefix (keys(%{$self->{transition_count}})) {
	$p++;
#	foreach my $suffix (keys (%suffix_sum)) {
	foreach my $suffix (keys (%{$self->{transition_count}->{$prefix}})) {
	    $s++;
	    my $pattern_count = $self->{transition_count}->{$prefix}->{$suffix};	    
	    $prefix_sum{$prefix} += $pattern_count;
	    $freq_sum += $pattern_count;
	    $suffix_sum{$suffix} += $pattern_count;
	}
    }
    $self->set_hash_attribute("prefix_sum", %prefix_sum);
    $self->set_hash_attribute("suffix_sum", %suffix_sum);

    ## Store prefixes and suffixes in arrays for quick access
    $self->set_array_attribute("prefixes", sort(keys(%prefix_sum)));
    $self->set_array_attribute("suffixes", sort(keys(%suffix_sum)));
    

    ## Calculate transition frequencies
    foreach my $prefix ($self->get_attribute("prefixes")) {
	foreach my $suffix ($self->get_attribute("suffixes")) {
#	foreach my $suffix (keys (%suffix_sum)) {
	    if (defined($self->{transition_count}->{$prefix}->{$suffix})) {	
		$self->{transition_freq}->{$prefix}->{$suffix} =  
		    $self->{transition_count}->{$prefix}->{$suffix}/$prefix_sum{$prefix};
	    } else {
		$self->{transition_freq}->{$prefix}->{$suffix} = 0;
		&RSAT::message::Warning(join(" ",
					     "No transition between prefix",$prefix, 
					     "and suffix", $suffix)) if ($main::verbose >= 1);
	    };
	}
    }


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

#    &RSAT::message::Debug("SUFFIX PROBA", $self->set_hash_attribute("suffix_proba", %suffix_proba));

    ## Average counts and frequencies on both strands if required
    my $strand = $self->get_attribute("strand");
    if ($strand eq "insensitive") {
	$self->average_strands();
    }

    &RSAT::message::TimeWarn(join("\t", 
				  "Normalized background model", 
				  "prefixes: ".$p,
				  "transitions: ".$s)) if ($main::verbose >= 2);
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
    
    if ($args{add}) {
	&RSAT::message::TimeWarn(join(" ", 
				      "Updating markov model (order ".$order.")",
				      "by adding sequence of length", $seq_len)) if ($main::verbose >= 2);
    } else {
	$self->set_hash_attribute("transition_count",());
	&RSAT::message::TimeWarn(join(" ", 
				      "Calculating markov model (order ".$order.")",
				      "from sequence of length", $seq_len)) if ($main::verbose >= 2);
    }


    &RSAT::error::FatalErrorr("RSAT::MarkovModel->calc_from_seq: sequence must be MUCH larger than order +1")
	unless ($seq_len > $order + 1);
    
    ## Transition counts
    my $last_pos = $seq_len - $order;
    for my $offset (0..($last_pos-1)) {
	my $prefix = substr($sequence, $offset, $order);
	my $suffix = substr($sequence, $offset + $order,1);
	$self->{transition_count}->{$prefix}->{$suffix}++;
    }

    ## Initialize transition frequencies
    $self->set_hash_attribute("transition_freq", %transition_count);

    ## Convert counts to transition frequencies
    $self->normalize_transition_frequencies();
   
    $self->force_attribute("training_words", $last_pos);
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
	&&($main::verbose >= 0)){
	&RSAT::message::Warning(join (" ", "Model update:", $added_word, 
				      "appeared in updated window starting at", $window_offset));
    }
    
    ## Update transition count for the deleted word
    my $deleted_prefix = substr($deleted_word, 0, $self->{order});
    my $deleted_suffix = substr($deleted_word, $self->{order}, 1);
    $self->{transition_count}->{$deleted_prefix}->{$deleted_suffix}--;
    if (($self->{transition_count}->{$deleted_prefix}->{$deleted_suffix} == 0) 
	&&($main::verbose >= 0)){
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

=item B<to_string()>

Return a string describing the transition matrix

=cut

sub to_string {
    my ($self, %args) = @_;
    my $string = "";
    my %prefix_proba = $self->get_attribute("prefix_proba");
    my @prefix = sort($self->get_attribute("prefixes"));
    my %suffix_sum = $self->get_attribute("suffix_sum");
    my @suffix = sort($self->get_attribute("suffixes"));


    ## Print header
    $string .= join ("\t", ";prefix",
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
	    $string .= sprintf "\t%.5f",  $self->{transition_freq}->{$prefix}->{$suffix};
	}
	$string .= sprintf "\t%.5f", $self->{prefix_proba}->{$prefix};
	$string .= sprintf "\t%.5f", $self->{prefix_sum}->{$prefix};
	$string .= "\n";
    }

    ## Print suffix probabilities
    my $row_name_len = &RSAT::stats::max(5,$self->get_attribute("order"));
    $string .= sprintf("; %${row_name_len}s", "P(su)");
#    $string .= "; P(su)";
    foreach my $suffix (@suffix) {
	$string .= sprintf "\t%.5f",  $self->{suffix_proba}->{$suffix};
    }
    $string .= "\n";

    ## Print suffix counts
    $string .= sprintf("; %${row_name_len}s", "N(su)");
#    $string .= "; N(su)";
    foreach my $suffix (@suffix) {
	$string .= sprintf "\t%.5f",  $self->{suffix_sum}->{$suffix};
    }
    $string .= "\n";


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
    
    &RSAT::message::Info(join("\t", "Matrix", $self->get_attribute("id"), "Averaging transition frequencies on both strands."))
	if ($main::verbose >= 2);

    ## Sum per prefix
    foreach my $prefix ($self->get_attribute("prefixes")) {
	my $prefix_rc = &main::SmartRC($prefix);	
	next if ($prefix_rc ge $prefix);	

	$self->{prefix_sum}->{$prefix} 
	= $self->{prefix_sum}->{$prefix_rc}
	= ($self->{prefix_sum}->{$prefix} +$self->{prefix_sum}->{$prefix_rc})/2 ;

	$self->{prefix_proba}->{$prefix} 
	= $self->{prefix_proba}->{$prefix_rc}
	= ($self->{prefix_proba}->{$prefix} +$self->{prefix_proba}->{$prefix_rc})/2 ;
    }

    ## Sum per suffix
    foreach my $suffix ($self->get_attribute("suffixes")) {
	my $suffix_rc = &main::SmartRC($suffix);	
	next if ($suffix_rc ge $suffix);

	$self->{suffix_sum}->{$suffix} 
	= $self->{suffix_sum}->{$suffix_rc}
	= ($self->{suffix_sum}->{$suffix} +$self->{suffix_sum}->{$suffix_rc})/2 ;

	$self->{suffix_proba}->{$suffix} 
	= $self->{suffix_proba}->{$suffix_rc}
	= ($self->{suffix_proba}->{$suffix} +$self->{suffix_proba}->{$suffix_rc})/2 ;
    }

    ## Transition matrix
    foreach my $prefix ($self->get_attribute("prefixes")) {
	my $prefix_rc = &main::SmartRC($prefix);
#	&RSAT::message::Debug("average_strands for prefix", $prefix, $prefix_rc) if ($main::verbose >= 10);
	next if ($prefix_rc ge $prefix);
	foreach my $suffix ($self->get_attribute("suffixes")) {
	    my $suffix_rc = &main::SmartRC($suffix);
	    next if ($suffix_rc ge $suffix);
	    $self->{transition_freq}->{$prefix}->{$suffix} 
	    = $self->{transition_freq}->{$prefix_rc}->{$suffix_rc} 
	    = ($self->{transition_freq}->{$prefix}->{$suffix} +
	       $self->{transition_freq}->{$prefix_rc}->{$suffix_rc})/2;
	    
	    $self->{transition_count}->{$prefix}->{$suffix} 
	    = $self->{transition_count}->{$prefix_rc}->{$suffix_rc} 
	    = ($self->{transition_count}->{$prefix}->{$suffix} +
	       $self->{transition_count}->{$prefix_rc}->{$suffix_rc})/2;
#	    &RSAT::message::Debug("average_strands", $prefix, $suffix, $prefix_rc, $suffix_rc) if ($main::verbose >= 10);
	}
    }
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
    if (defined($self->{prefix_proba}->{$prefix})) {
	$segment_proba = $self->{prefix_proba}->{$prefix};
    }
    
    for my $c ($order..($seq_len-1)) {
	my $letter_proba = 0;

	my $suffix = substr($segment, $c, 1);
	my $prefix = substr($segment,($c-$order),$order);
	if (defined($self->{transition_freq}->{$prefix}->{$suffix})) {
	    $letter_proba = $self->{transition_freq}->{$prefix}->{$suffix};
	}
	
#	my $word = substr($segment,($c-$order),$order+1);
#	if (defined($self->{transition_quick}->{$word})) {
#	    $letter_proba = $self->{transition_quick}->{$word};
#	}
	
#	&RSAT::message::Debug("letter proba", $word,$letter_proba) if ($main::verbose >= 0);

	$segment_proba *= $letter_proba;

#	&RSAT::message::Debug("segment_proba", 
#			      "prefix=".$word, 
#			      "prefix=".$prefix, 
#			      "suffix:".$suffix, 
#			      "offset:".$c, 
#			      "P(letter)=".$letter_proba, 
#			      "P(segm)=".$segment_proba) if ($main::verbose >= 0);
    }
    
#    &RSAT::message::Debug("segment_proba", $segment, "P(segm)=".$segment_proba) if ($main::verbose >= 10);
#    die;
    return $segment_proba;
}


return 1;

__END__

