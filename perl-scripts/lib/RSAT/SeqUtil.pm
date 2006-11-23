###############################################################
#
# Class SeqUtil
#
package RSAT::SeqUtil;

use RSAT::GenericObject;
use RSAT::error;
@ISA = qw( RSAT::GenericObject );

### class attributes

=pod

=head1 NAME

    RSAT::SeqUtil

=head1 DESCRIPTION

Abstract class containing a set of methods for trating sequences.

=cut


################################################################
=pod

=item get_accepted_residues()

Returns a hash table with the accepted residues, given the sequence
type (DNA or protein). 

Usage:  my %accepted_residues = &RSAT::SeqUtil::get_accepted_residues($seq_type);

=cut
sub get_accepted_residues {
  my ($seq_type) = @_;
  &RSAT::error::FatalError("The method &RSAT::SeqUtil::get_accepted_residues() requires to specify a sequence type.") 
    unless ($seq_type);

  my @dna_alphabet = qw (a c g t);
  my @protein_alphabet = qw(a c d e f g h i k l m n p q r s t v w y);
  my @alphabet = ();
  my %accepted_residues = ();

  if (lc($seq_type) eq "dna") {
    @alphabet = @dna_alphabet;
  } elsif (lc($seq_type) eq "protein") {
    @alphabet = @protein_alphabet;
  } else {
    &RSAT::error::FatalError("The method &RSAT::SeqUtil::get_accepted_residues()", "Invalid sequence type.", "Supported: dna,protein");
  }

  foreach my $residue (@alphabet) {
    $accepted_residues{$residue} = 0.25;
  }
  
  return(%accepted_residues);
}


return 1;

__END__

