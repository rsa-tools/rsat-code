###############################################################
#
# Class MarkovModel
#
package RSAT::MarkovModel;

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
    my %prefix_sum = ();
    my %prefix_proba = ();
    my %suffix_sum = ();
    my @patterns = keys(%patterns);
    my $freq_sum = 0;

    my $order = length($patterns[0]) -1 ; ## Use the first pattern to calculate model order
    $self->force_attribute("order", $order);

    #### Calculate alphabet from expected frequency keys
    foreach my $pattern_seq (keys %patterns) {
	## Check pattern length
	$pattern_seq = lc($pattern_seq);
	my $pattern_len = length($pattern_seq);
	my $pattern_freq =  $patterns{$pattern_seq}->{exp_freq};
	&RSAT::error::FatalError("All patterns should have the same length in a Markov model file.") 
	    unless $pattern_len = $order+1;
	my $prefix = substr($pattern_seq,0,$order);
	my $suffix = substr($pattern_seq,$order, 1);
	$self->{transition}->{$prefix}->{$suffix} = $pattern_freq;
	$prefix_sum{$prefix} += $pattern_freq;
	$freq_sum += $pattern_freq;
	$suffix_sum{$suffix} += $pattern_freq;
    }
    
    ## Calculate transition probabilities
    foreach my $prefix (keys(%prefix_sum)) {
	foreach my $suffix (keys (%suffix_sum)) {
	    if (defined($self->{transition}->{$prefix}->{$suffix})) {	
		$self->{transition}->{$prefix}->{$suffix} /= $prefix_sum{$prefix};
	    } else {
		$self->{transition}->{$prefix}->{$suffix} = 0;
	    };
	}
    }

    ## Calculate prefix probabilities
    foreach my $prefix (keys(%prefix_sum)) {
	$prefix_proba{$prefix} = $prefix_sum{$prefix}/$freq_sum;
    }    

    ## store prefix sums and suffix sums
    $self->set_hash_attribute("prefix_proba", %prefix_proba);
    $self->set_hash_attribute("prefix_sum", %prefix_sum);
    $self->set_hash_attribute("suffix_sum", %suffix_sum);
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
    my @prefix = sort(keys(%prefix_proba));
    my %suffix_sum = $self->get_attribute("suffix_sum");
    my @suffix = sort(keys(%suffix_sum));
    
    $string .= join ("\t", ";prefix",
		     @suffix,
		     "P(pref)");
    $string .= "\n";

    foreach my $prefix (@prefix) {
	if (defined($args{comment_string})) {
	    $string .= $args{comment_string};
	}
	$string .= $prefix;
	foreach my $suffix (@suffix) {
	    $string .= sprintf "\t%.5f",  $self->{transition}->{$prefix}->{$suffix};
	}
	$string .= sprintf "\t%.5f", $prefix_proba{$prefix};
	$string .= "\n";
    }

    return $string;
}



################################################################
=pod 

=item B<segment_proba($segment)>

Calculate the probability of a segment of sequence. The sequence segment must
by larger than the Markov order + 1.

=cut

sub segment_proba {
    my ($self, $segment) = @_;

    my $seq_len = length($segment);
    my $order = $self->get_attribute("order");

    &RSAT::error::FatalError("&RSAT::MarkovModel->segment_proba. The segment ($segment) length ($seq_len) must be larger than the markov order ($order) + 1.") 
	if ($seq_len < $order + 1);

    my $prefix = substr($segment,0,$order);
    my $segment_proba = 0;
    if (defined($self->{prefix_proba}->{$prefix})) {
	$segment_proba = $self->{prefix_proba}->{$prefix};
    }
    
    for my $c ($order..($seq_len-1)) {
	my $suffix = lc(substr($segment, $c, 1));
	my $prefix = substr($segment,($c-$order),$order);
	my $letter_proba = 0;
	if (defined($self->{transition}->{$prefix}->{$suffix})) {
	    $letter_proba = $self->{transition}->{$prefix}->{$suffix};
	}
	$segment_proba *= $letter_proba;
#	&RSAT::message::Debug("segment_proba", 
#			      "prefix=".$prefix, 
#			      "suffix:".$suffix, 
#			      "offset:".$c, 
#			      "P(letter)=".$letter_proba, 
#			      "P(segm)=".$segment_proba) if ($main::verbose >= 0);
    }
    
    &RSAT::message::Debug("segment_proba", $segment, "P(segm)=".$segment_proba) if ($main::verbose >= 10);
#    die;
    return $segment_proba;
}


return 1;

__END__

