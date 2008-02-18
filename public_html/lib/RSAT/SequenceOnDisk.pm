###############################################################
#
# Class SequenceOnDisk
#
package RSAT::SequenceOnDisk;

use RSAT::GenericObject;
use RSAT::error;
#@ISA = qw( RSAT::GenericObject );
@ISA = qw( RSAT::GenericObject RSAT::Sequence );

=pod

=head1 NAME

    RSAT::SequenceOnDisk

=head1 DESCRIPTION

SequenceOnDisk provides direct read access to a sequence on disk this
is particularly important for retrieving sub-sequences of large
contigs whose size could exceed the RAM (example whole human
chromosomes).

=cut


################################################################
=pod

=item new()

Create a new SequenceOnDisk.

=cut
sub new {
    my ($class, %args) = @_;
    
    my %new_seq = ();
    unless (defined($args{filename})) {
	die "Error: cannot open a sequence on disk without a specified filename\n";
    }
    
    %new_seq = %args;
    unless (defined($new_seq{id})) {
	$new_seq{id} = $class->auto_id();
    }
    my $self = bless \%new_seq, $class;
    $self->init();
    
    return $self;
}
  
  sub init {
    my ($self) = @_;
    my $filename = $self->get_attribute("filename");

    if ($filename) {
      die "Error: file $filename does not exist\n"
	unless (-e $filename);
      die "Error: cannot read file $filename\n"
	unless (-r $filename);

      ### open a stream to the file
## TO DEBUG: this does not work !

#      open FILE, $filename ||
#	die "Error: cannot open file ", $filename, "\n";
#      #    no strict;
#      my $filehandle = FILE;
#      $self->set_attribute("filehandle", $filehandle);      
#      return $filehandle;

    } else { 
      &RSAT::error::FatalError("A file name has to be specified");
    }
  }

  sub DESTROY {
    my ($self) = @_;
### TO DEBUG: THIS DOES NOT WORK !
#    my $filehandle =  $self->get_attribute("filehandle");
#    my $filename =  $self->get_attribute("filename");
#    close $filehandle ||
#	die "Error: cannot close file $filename\n";;
  }

  sub get_length() {
    my ($self) = @_;
    my $filename = $self->get_attribute("filename");
    unless (defined($self->{length})) {
    ## here changed by Morgane : with -1, the last nucleotide is missed
    #  my $length = (-s $filename) -1;
    ##
    my $length = (-s $filename);
      $self->set_attribute("length", $length);
    }
    return $self->get_attribute("length");
  }


  ### retrieve a sub-sequence
  sub get_sequence {
    my ($self, $from, $to, $strand) = @_;

#    no strict;
#    my $filehandle =  $self->get_attribute("filehandle");

    my $sequence = "";

    ### check from and to
    $from = 1 unless (defined($from));
    $to = $self->get_length() unless (defined($to));
    if ($from > $to) { 
      my $tmp = $to;
      $to = $from;
      $from = $tmp;
    }

    warn (join "\t",
	  "; get_sequence",
	  $self->get_attribute("filename"),
	  $from,
	  $to,
	  $strand,
	  $sequence,
	  "\n")
	if ($main::verbose >= 3);


    ### treat out-of-bonds cases

    ### overlap sequence origin
    if (($from < 0) && 
	($to > 0)) {
      if ($self->get_attribute("circular") == 1) { 
	### concatenate end and beginning of circular sequence
	$sequence = $self->get_sequence($self->get_length() + $from, $self->get_length());
	$sequence .= $self->get_sequence(1,$to);
      } else {
	### truncate and warn
	warn join "\t", "; WARNING: cannot retrieve sub-sequence with negative limits", $self->get_id(), $from, $to, "\n";
	$sequence = $self->get_sequence(1,$to);
      }
      
    ### treat out-of-bonds cases      ### overlap sequence end
    } elsif (($from < $self->get_length()) &&
	     ($to > $self->get_length())) {
    
      if ($self->get_attribute("circular") == 1) { 
	### concatenate end and beginning of circular sequence
	$sequence = $self->get_sequence($from,$self->get_length());
	$sequence .= $self->get_sequence(1,$to - $self->get_length());
      } else {
	### truncate and warn
	warn join "\t", "; WARNING: cannot retrieve sub-sequence with limits larger than sequence length", $self->get_id(), $from, $to, "\n";
	$sequence = $self->get_sequence($from,$self->get_length());
      }


      ### negative coordinates
    } elsif (($from < 0) &&
	     ($to < 0)) {
      if ($self->get_attribute("circular") == 1) { 
	$from += $self->get_length();
	$to += $self->get_length();
	$sequence =  $self->get_sequence($from,$to);
      } else {
	warn "Error: cannot retrieve sub-sequence with negative limits\n";
	return undef;
      }

      ### coordinates larger than sequence length
    } elsif (($from > $self->get_length()) && 
	     ($to > $self->get_length())) {
      if ($self->get_attribute("circular") == 1) { 
	$from -= $self->get_length();
	$to -= $self->get_length();
	$sequence = $self->get_sequence($from,$to);
      } else {
	warn "Error: cannot retrieve sub-sequence with limits larger than whole sequence size\n";
	return undef;
      }

      ### normal case
    } else { 
      my $len = $to - $from + 1;
      my $offset = $from -1;
      
      ### open a stream to the file
### TO DEBUG: there is a problem with file hndles: I cannot attach a file handle to the object, I don't understand why

#      $filehandle =  $self->get_attribute("filehandle");
#      unless ($filehandle) {
	  warn "Opening file handle\n" if ($main::verbose >= 5);
	  my $filename = $self->get_attribute("filename");
	  open FH, $filename ||
	      die "Error: cannot open file ", $filename, "\n";
	  $filehandle = FH;
#	  $self->set_attribute("filehandle", FH);
#	  $self->set_attribute("filehandle", $filehandle);
#      }
      
#      sysseek($filehandle, $offset, 0) &&
#	  sysread($filehandle, $sequence, $len);
#      close $filehandle;

      sysseek(FH, $offset, 0) &&
	  sysread(FH, $sequence, $len);
      close FH;

      warn (join "\t",
	    "; get_sequence",
	    FH,
	    $from,
	    $to,
	    $len,
	    $offset,
	    "\n")
	  if ($main::verbose >= 3);

    }
    

    if ($strand eq "R") { ### reverse complement
      return (&main::ReverseComplement($sequence));
    } else {
      return $sequence;
    }
  }

################################################################
=pod

=item get_file_name

Return the name of the file containing the sequence.

=cut
sub get_file_name {
    my ($self) = @_;
    return $self->{name};
}

################################################################
=pod

=item set_file_name

Specify the name of the file containing the sequence.

=cut
sub set_file_name {
    my ($self, $new_name) = @_;
    $self->{name} = $new_name;
}

return 1;


__END__

